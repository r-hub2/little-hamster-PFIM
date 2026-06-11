// PGBO kernel (Le Nagard et al.): mutate one sampling group per attempt,
// accept via Metropolis-like rule on fitness = 10^(-theta * D).
// eval_d returns D; check_valid_group enforces arm constraints from R.
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;

static void sort_group(NumericVector& pos, const IntegerVector& idx_one_based) {
  const int n = idx_one_based.size();
  std::vector<double> buf(static_cast<size_t>(n));
  for (int k = 0; k < n; ++k) buf[static_cast<size_t>(k)] = pos[idx_one_based[k] - 1];
  std::sort(buf.begin(), buf.end());
  for (int k = 0; k < n; ++k) pos[idx_one_based[k] - 1] = buf[static_cast<size_t>(k)];
}

static int random_group_index(int n_groups) {
  const int g = static_cast<int>(std::floor(R::runif(0.0, static_cast<double>(n_groups))));
  return (g >= n_groups) ? (n_groups - 1) : g;
}

// [[Rcpp::export]]
List pgbo_optimize_Rcpp(NumericVector initial_pos,
                        List sorting_groups,
                        int max_iteration,
                        int N,
                        double mute_effect,
                        int purge_iteration,
                        double fit_base,
                        int max_attempts,
                        bool show_process,
                        Function eval_d,
                        Function check_valid_group) {
  NumericVector pos = clone(initial_pos);
  const int n_groups = sorting_groups.size();

  std::vector<IntegerVector> group_indices(static_cast<size_t>(n_groups));
  for (int g = 0; g < n_groups; ++g)
    group_indices[static_cast<size_t>(g)] = sorting_groups[g];

  const double initial_d = as<double>(eval_d(pos));
  double theta = -std::log10(fit_base) / initial_d;
  double fit_a = std::pow(10.0, -theta * initial_d);
  double fit_best = fit_a;

  NumericVector best_pos = clone(pos);

  std::vector<char> has_group(static_cast<size_t>(n_groups), 0);

  for (int iteration = 1; iteration <= max_iteration; ++iteration) {
    Rcpp::checkUserInterrupt();

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
      if (R::runif(0.0, 1.0) < 0.8) {
        const int local = static_cast<int>(std::floor(R::runif(0.0, static_cast<double>(grp_idx.size()))));
        const int pick = (local >= static_cast<int>(grp_idx.size())) ? (static_cast<int>(grp_idx.size()) - 1) : local;
        trial[grp_idx[pick] - 1] += R::rcauchy(0.0, 1.0) * mute_effect;
      } else {
        for (int k = 0; k < grp_idx.size(); ++k) {
          trial[grp_idx[k] - 1] += R::rnorm(0.0, 1.0) * mute_effect;
        }
      }
      sort_group(trial, grp_idx);

      if (as<bool>(check_valid_group(trial, group_idx + 1))) {
        NumericVector stored(grp_idx.size());
        for (int k = 0; k < grp_idx.size(); ++k) stored[k] = trial[grp_idx[k] - 1];
        group_values[static_cast<size_t>(group_idx)] = stored;
        has_group[static_cast<size_t>(group_idx)] = 1;
        found = std::all_of(has_group.begin(), has_group.end(),
                            [](char v) { return v != 0; });
      }
    }

    NumericVector candidate = clone(pos);
    for (int g = 0; g < n_groups; ++g) {
      const IntegerVector& grp_idx = group_indices[static_cast<size_t>(g)];
      const NumericVector& stored = group_values[static_cast<size_t>(g)];
      for (int k = 0; k < grp_idx.size(); ++k) {
        candidate[grp_idx[k] - 1] = stored[k];
      }
    }

    const double d_b = as<double>(eval_d(candidate));
    const double fit_b = std::pow(10.0, -theta * d_b);

    if (!R_IsNaN(fit_b) && fit_b > 0.0) {
      double proba = 0.0;
      if (fit_a == fit_b) {
        proba = 1.0 / static_cast<double>(N);
      } else {
        const double f = fit_a / fit_b;
        proba = (1.0 - f * f) / (1.0 - std::pow(f, 2.0 * static_cast<double>(N)));
        if (R_IsNaN(proba)) proba = 0.0;
      }

      if (R::runif(0.0, 1.0) < proba) {
        fit_a = fit_b;
        pos = clone(candidate);

        const double inv_d = 1.0 / d_b;
        if (fit_best < inv_d) {
          fit_best = inv_d;
          best_pos = clone(candidate);
          if (show_process) {
            Rcout << "Iteration = " << iteration
                  << " | Criterion = " << inv_d << "\n";
          }
        }
      }
    }

    // Rescale theta so the acceptance temperature does not freeze (fit_a reset).
    if (purge_iteration > 0 && (iteration % purge_iteration) == 0) {
      const double d_tmp = -std::log10(fit_a) / theta;
      theta = -std::log10(fit_base) / d_tmp;
      fit_a = fit_base;
    }
  }

  return List::create(
    Named("bestDesign") = best_pos,
    Named("fitBest")    = fit_best
  );
}
