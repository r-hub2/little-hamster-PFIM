// PGBO kernel (Le Nagard et al., Population Genetics Based Optimizer).
//
// Contract with R:
//   - `pos` is a flat sampling-time vector
//   - `eval_d(pos)` returns the D-criterion (higher is better); R owns FIM assembly
//   - cost = 1/D; fitness = 10^(-theta * cost)  (same cost geometry as PSO)
//   - `check_valid_group(trial, group_1based)` enforces arm / window constraints
//   - mutate one sampling group per attempt; Metropolis-like accept on fitness
//   - purge every purge_iteration: rescale theta from current fit_a, reset fit_a
//     so acceptance does not freeze as the search improves
//
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;

// --- Helpers ---------------------------------------------------------------

// Sort the coordinates of one sampling group (1-based indices from R).
static void sort_group(NumericVector& pos, const IntegerVector& idx_one_based) {
  const int n = idx_one_based.size();
  std::vector<double> buf(static_cast<size_t>(n));
  for (int k = 0; k < n; ++k) buf[static_cast<size_t>(k)] = pos[idx_one_based[k] - 1];
  std::sort(buf.begin(), buf.end());
  for (int k = 0; k < n; ++k) pos[idx_one_based[k] - 1] = buf[static_cast<size_t>(k)];
}

// Load and validate one sorting group (1-based indices into pos).
static IntegerVector load_sorting_group(SEXP group_sexp, int n_pos) {
  IntegerVector grp = group_sexp;
  if (grp.size() == 0)
    stop("PGBOAlgorithm: sorting_groups must not contain empty groups.");
  for (int k = 0; k < grp.size(); ++k) {
    const int idx = grp[k];
    if (IntegerVector::is_na(idx) || idx < 1 || idx > n_pos)
      stop("PGBOAlgorithm: sorting_groups indices must be in 1..length(pos).");
  }
  return grp;
}

// Uniform draw of a group index in [0, n_groups).
static int random_group_index(int n_groups) {
  const int g = static_cast<int>(std::floor(R::runif(0.0, static_cast<double>(n_groups))));
  return (g >= n_groups) ? (n_groups - 1) : g;
}

// --- Exported PGBO loop ----------------------------------------------------
// [[Rcpp::export]]
List pgbo_optimize_Rcpp(NumericVector initial_pos,
                        List sorting_groups,
                        int max_iteration,
                        int N,
                        double mute_effect,
                        int purge_iteration,
                        double fit_base,
                        int max_attempts,
                        double cauchy_prob,
                        bool show_process,
                        Function eval_d,
                        Function check_valid_group,
                        double ftol = 0.0,
                        int stall_iterations = 5) {
  if (!R_finite(cauchy_prob) || cauchy_prob < 0.0 || cauchy_prob > 1.0)
    stop("PGBOAlgorithm: cauchy_prob must be in [0, 1].");
  NumericVector pos = clone(initial_pos);
  const int n_groups = sorting_groups.size();
  if (n_groups <= 0)
    stop("PGBOAlgorithm: sorting_groups must be non-empty.");

  std::vector<IntegerVector> group_indices(static_cast<size_t>(n_groups));
  for (int g = 0; g < n_groups; ++g)
    group_indices[static_cast<size_t>(g)] =
      load_sorting_group(sorting_groups[g], pos.size());

  const double initial_d = as<double>(eval_d(pos));
  if (!R_finite(initial_d) || initial_d <= 0.0)
    stop("PGBOAlgorithm: initial design has invalid D-criterion.");

  // Fitness uses cost = 1/D (same convention as PSO).
  // fit = 10^(-theta * cost); theta chosen so initial fit equals fit_base.
  const double initial_cost = 1.0 / initial_d;
  double theta = -std::log10(fit_base) / initial_cost;
  double fit_a = std::pow(10.0, -theta * initial_cost);
  double best_d = initial_d;

  NumericVector best_pos = clone(pos);

  std::vector<char> has_group(static_cast<size_t>(n_groups), 0);

  bool converged = false;
  int iterations_done = 0;
  int stall_count = 0;
  const int stall_need = std::max(1, stall_iterations);

  for (int iteration = 1; iteration <= max_iteration; ++iteration) {
    Rcpp::checkUserInterrupt();
    iterations_done = iteration;
    const double iter_start_best = best_d;

    // --- Collect one valid mutation per arm group --------------------------
    std::fill(has_group.begin(), has_group.end(), 0);
    std::vector<NumericVector> group_values(static_cast<size_t>(n_groups));

    bool found = false;
    int attempts = 0;

    while (!found) {
      ++attempts;
      if (attempts > max_attempts) {
        stop("PGBOAlgorithm: could not find a valid sampling configuration within %d attempts.\n"
             "Check that the sampling constraints allow enough freedom, or reduce muteEffect (current: %g).",
             max_attempts, mute_effect);
      }

      const int group_idx = random_group_index(n_groups);
      const IntegerVector& grp_idx = group_indices[static_cast<size_t>(group_idx)];

      NumericVector trial = clone(pos);

      // Mostly Cauchy jump on one time point; otherwise Gaussian on the whole group.
      if (R::runif(0.0, 1.0) < cauchy_prob) {
        const int local = static_cast<int>(std::floor(R::runif(0.0, static_cast<double>(grp_idx.size()))));
        const int pick = (local >= static_cast<int>(grp_idx.size())) ? (static_cast<int>(grp_idx.size()) - 1) : local;
        trial[grp_idx[pick] - 1] += R::rcauchy(0.0, 1.0) * mute_effect;
      } else {
        for (int k = 0; k < grp_idx.size(); ++k) {
          trial[grp_idx[k] - 1] += R::rnorm(0.0, 1.0) * mute_effect;
        }
      }
      sort_group(trial, grp_idx);

      // Accept this group's mutation if R says the constraints are satisfied.
      if (as<bool>(check_valid_group(trial, group_idx + 1))) {
        NumericVector stored(grp_idx.size());
        for (int k = 0; k < grp_idx.size(); ++k) stored[k] = trial[grp_idx[k] - 1];
        group_values[static_cast<size_t>(group_idx)] = stored;
        has_group[static_cast<size_t>(group_idx)] = 1;
        found = std::all_of(has_group.begin(), has_group.end(),
                            [](char v) { return v != 0; });
      }
    }

    // --- Assemble full candidate from accepted per-group mutations ---------
    NumericVector candidate = clone(pos);
    for (int g = 0; g < n_groups; ++g) {
      const IntegerVector& grp_idx = group_indices[static_cast<size_t>(g)];
      const NumericVector& stored = group_values[static_cast<size_t>(g)];
      for (int k = 0; k < grp_idx.size(); ++k) {
        candidate[grp_idx[k] - 1] = stored[k];
      }
    }

    const double d_b = as<double>(eval_d(candidate));

    // Track the best D seen so far (for reporting / return value).
    if (R_finite(d_b) && d_b > best_d) {
      best_d = d_b;
      best_pos = clone(candidate);
      if (show_process)
        Rcout << "Iteration = " << iteration << " | Criterion = " << best_d << "\n";
    }

    // --- Metropolis-like acceptance on fitness (not raw D) -----------------
    if (R_finite(d_b) && d_b > 0.0) {
      const double cost_b = 1.0 / d_b;
      const double fit_b = std::pow(10.0, -theta * cost_b);

      if (fit_b > 0.0) {
        double proba = 0.0;
        if (fit_a == fit_b) {
          proba = 1.0 / static_cast<double>(N);
        } else {
          // Population-genetics acceptance probability (Le Nagard).
          const double f = fit_a / fit_b;
          proba = (1.0 - f * f) / (1.0 - std::pow(f, 2.0 * static_cast<double>(N)));
          if (R_IsNaN(proba)) proba = 0.0;
        }

        if (R::runif(0.0, 1.0) < proba) {
          fit_a = fit_b;
          pos = clone(candidate);
        }
      }
    }

    // --- Periodic "purge": rescale theta so acceptance does not freeze -----
    // Recover cost from fit_a = 10^(-theta * cost), then rescale; reset fit_a.
    if ( purge_iteration > 0 && ( iteration % purge_iteration ) == 0 ) {
      if ( fit_a > 0.0 && fit_a < 1.0 ) {
        const double cost_tmp = -std::log10( fit_a ) / theta;
        if ( std::isfinite( cost_tmp ) && cost_tmp > 0.0 )
          theta = -std::log10( fit_base ) / cost_tmp;
      }
      fit_a = fit_base;
    }

    // Relative improvement stall stop on best_d within this iteration.
    if (ftol > 0.0 && R_finite(best_d) && R_finite(iter_start_best) && iter_start_best > 0.0) {
      const double rel = std::abs(best_d - iter_start_best) / iter_start_best;
      if (rel < ftol)
        ++stall_count;
      else
        stall_count = 0;
      if (stall_count >= stall_need) {
        converged = true;
        break;
      }
    }
  }

  const bool improved = R_finite(best_d) && best_d > initial_d * (1.0 + 1e-12);

  // converged is only set by the stall stop (ftol > 0); never equate to improved.
  return List::create(
    Named("bestDesign") = best_pos,
    Named("bestD")      = best_d,
    Named("initialD")   = initial_d,
    Named("converged")  = converged,
    Named("iterations") = iterations_done,
    Named("improved")   = improved
  );
}
