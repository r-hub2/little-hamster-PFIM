// Particle Swarm Optimization on flat sampling-time vectors.
// Clerc & Kennedy constriction PSO.
//
// Contract with R:
//   - each particle is a row of the position matrix (flat sampling times)
//   - cost = eval_fitness(row) = 1/D  (lower is better); optional batch API
//   - after every velocity update: project onto sampling windows, then
//     re-sort times within each arm group so chronological order is kept
//   - empty window matrix for a coordinate = fixed (pinned to initial_pos)
//   - FIM / Bayesian prior live in the R fitness callback, not here
//
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <vector>

using namespace Rcpp;

// --- Constraint helpers ----------------------------------------------------

// Project a scalar onto the nearest admissible sampling window.
// Each row of `bounds` is [w_min, w_max]; if already inside any window, return as-is.
// Empty bounds here mean "no clamp" — callers treat empty windows_list entries as
// pinned coordinates and rewrite those entries to initial_pos after clamp/sort.
static double clamp_to_windows(double val, const arma::mat& bounds) {
  if ( bounds.n_rows == 0 )
    return val;

  double min_dist = arma::datum::inf;
  double clamped_val = val;

  for (arma::uword k = 0; k < bounds.n_rows; ++k) {
    const double w_min = bounds(k, 0);
    const double w_max = bounds(k, 1);

    if (val >= w_min && val <= w_max) return val;

    const double d1 = std::abs(val - w_min);
    const double d2 = std::abs(val - w_max);

    if (d1 < min_dist) { min_dist = d1; clamped_val = w_min; }
    if (d2 < min_dist) { min_dist = d2; clamped_val = w_max; }
  }
  return clamped_val;
}

// Within-arm sampling times must stay ordered; indices are 1-based in R input,
// converted to 0-based Armadillo indices before calling this helper.
static void sort_groups(arma::rowvec& row, const std::vector<arma::uvec>& sorting_indices) {
  for (size_t g = 0; g < sorting_indices.size(); ++g) {
    const arma::uvec& idx = sorting_indices[g];
    arma::rowvec sub_pos(idx.n_elem);
    for (arma::uword k = 0; k < idx.n_elem; ++k) sub_pos(k) = row(idx(k));
    sub_pos = arma::sort(sub_pos);
    for (arma::uword k = 0; k < idx.n_elem; ++k) row(idx(k)) = sub_pos(k);
  }
}

// Print iteration progress; show D = 1/cost when cost is a valid positive finite value.
static void report_pso_progress(int iteration, double gbest_cost) {
  Rcout << "Iteration = " << iteration;
  if (R_finite(gbest_cost) && gbest_cost > 0.0 && gbest_cost < 1e6) {
    Rcout << " | Criterion = " << (1.0 / gbest_cost) << "\n";
  } else {
    Rcout << " | Criterion = N/A (no valid design yet)\n";
  }
}

// Prefer eval_fitness_batch (one R call for the whole swarm) when provided;
// otherwise fall back to per-row scalar callbacks.
static arma::vec eval_fitness_rows(const Nullable<Function>& eval_batch,
                                   const Function& eval_scalar,
                                   const arma::mat& positions) {
  if (eval_batch.isNotNull()) {
    Function batch_fn = eval_batch.get();
    arma::vec costs = as<arma::vec>(batch_fn(wrap(positions)));
    if (static_cast<arma::uword>(costs.n_elem) != positions.n_rows)
      stop("pso_optimize_Rcpp: batch fitness length must equal population size.");
    return costs;
  }
  arma::vec costs(positions.n_rows);
  for (arma::uword i = 0; i < positions.n_rows; ++i) {
    costs(i) = as<double>(eval_scalar(wrap(positions.row(i))));
  }
  return costs;
}

// Update personal bests and global best from a freshly evaluated cost vector.
static void update_best_from_costs(const arma::vec& costs,
                                   const arma::mat& positions,
                                   arma::mat& pbest_pos,
                                   arma::vec& pbest_cost,
                                   arma::rowvec& gbest_pos,
                                   double& gbest_cost) {
  for (arma::uword i = 0; i < costs.n_elem; ++i) {
    const double cost = costs(i);
    if (cost < pbest_cost(i)) {
      pbest_cost(i) = cost;
      pbest_pos.row(i) = positions.row(i);
    }
    if (cost < gbest_cost) {
      gbest_cost = cost;
      gbest_pos  = positions.row(i);
    }
  }
}

// Convert R 1-based sorting group to validated 0-based Armadillo indices.
static arma::uvec sorting_group_indices(const IntegerVector& raw, arma::uword n_dims) {
  if (raw.size() == 0)
    stop("pso_optimize_Rcpp: sorting_groups must not contain empty groups.");
  arma::uvec idx(static_cast<arma::uword>(raw.size()));
  for (int k = 0; k < raw.size(); ++k) {
    const int one_based = raw[k];
    if (IntegerVector::is_na(one_based) || one_based < 1 ||
        static_cast<arma::uword>(one_based) > n_dims)
      stop("pso_optimize_Rcpp: sorting_groups indices must be in 1..n_dims.");
    idx(static_cast<arma::uword>(k)) = static_cast<arma::uword>(one_based - 1);
  }
  return idx;
}

// --- Exported constriction PSO ---------------------------------------------
// [[Rcpp::export]]
List pso_optimize_Rcpp(int n_pop_in,
                       int max_iter,
                       arma::rowvec initial_pos,
                       List windows_list,
                       List sorting_groups,
                       double phi1,
                       double phi2,
                       double constriction,
                       bool show_process,
                       Function eval_fitness,
                       Nullable<Function> sample_valid_pos = R_NilValue,
                       Nullable<Function> eval_fitness_batch = R_NilValue,
                       double ftol = 0.0,
                       int stall_iterations = 5) {
  if (n_pop_in <= 0)
    stop("pso_optimize_Rcpp: n_pop_in must be positive.");
  const arma::uword n_pop  = static_cast<arma::uword>(n_pop_in);
  const arma::uword n_dims = initial_pos.n_elem;
  const bool use_valid_sampler = sample_valid_pos.isNotNull();

  if (windows_list.size() != static_cast<int>(n_dims))
    stop("pso_optimize_Rcpp: windows_list length must equal length(initial_pos).");

  // Per-dimension window matrices (empty = fixed coordinate pinned to initial).
  std::vector<arma::mat> windows(static_cast<size_t>(n_dims));
  for (arma::uword j = 0; j < n_dims; ++j) {
    const arma::mat bounds = as<arma::mat>(windows_list[static_cast<int>(j)]);
    if (bounds.n_rows > 0 && bounds.n_cols < 2)
      stop("pso_optimize_Rcpp: each non-empty window matrix must have >= 2 columns.");
    windows[static_cast<size_t>(j)] = bounds;
  }

  // Convert R 1-based group indices to 0-based Armadillo indices.
  std::vector<arma::uvec> sorting_indices(static_cast<size_t>(sorting_groups.size()));
  for (int g = 0; g < sorting_groups.size(); ++g)
    sorting_indices[static_cast<size_t>(g)] =
      sorting_group_indices(sorting_groups[g], n_dims);

  // Swarm state: positions, velocities, personal / global bests.
  arma::mat pos = arma::repmat(arma::mat(initial_pos), n_pop, 1);
  arma::mat vel = arma::zeros<arma::mat>(n_pop, n_dims);
  arma::mat pbest_pos = pos;
  arma::vec pbest_cost(n_pop);
  pbest_cost.fill(arma::datum::inf);

  arma::rowvec gbest_pos = initial_pos;
  double gbest_cost = arma::datum::inf;
  double initial_gbest_cost = arma::datum::inf;

  // --- Seed global best from the ordered initial design --------------------
  {
    arma::rowvec init_row = initial_pos;
    sort_groups(init_row, sorting_indices);
    arma::mat init_mat(1, n_dims, arma::fill::zeros);
    init_mat.row(0) = init_row;
    const arma::vec init_costs = eval_fitness_rows(
      eval_fitness_batch, eval_fitness, init_mat
    );
    if (init_costs(0) < gbest_cost) {
      gbest_cost = init_costs(0);
      gbest_pos  = init_row;
    }
  }

  // --- Initialise the swarm ------------------------------------------------
  // Optional R sampler supplies constraint-feasible designs; otherwise draw
  // uniformly inside each coordinate's windows, then clamp + sort.
  for (arma::uword i = 0; i < n_pop; ++i) {
    if (use_valid_sampler) {
      Function sampler = sample_valid_pos.get();
      pos.row(i) = as<arma::rowvec>(sampler());
    } else {
      for (arma::uword j = 0; j < n_dims; ++j) {
        const arma::mat& bounds = windows[static_cast<size_t>(j)];
        if ( bounds.n_rows == 0 ) {
          pos(i, j) = initial_pos(j);
        } else {
          const double w_min = bounds.col(0).min();
          const double w_max = bounds.col(1).max();
          pos(i, j) = w_min + R::runif(0.0, 1.0) * (w_max - w_min);
          pos(i, j) = clamp_to_windows(pos(i, j), bounds);
        }
      }
    }

    arma::rowvec current_row = pos.row(i);
    sort_groups(current_row, sorting_indices);
    for (arma::uword j = 0; j < n_dims; ++j) {
      if ( windows[static_cast<size_t>(j)].n_rows == 0 )
        current_row(j) = initial_pos(j);
    }
    pos.row(i) = current_row;
  }

  // Evaluate initial swarm and set personal / global bests.
  {
    const arma::vec init_costs = eval_fitness_rows(eval_fitness_batch, eval_fitness, pos);
    pbest_cost = init_costs;
    pbest_pos  = pos;
    update_best_from_costs(init_costs, pos, pbest_pos, pbest_cost, gbest_pos, gbest_cost);
    initial_gbest_cost = gbest_cost;
  }

  bool converged = false;
  int iterations_done = 0;
  int stall_count = 0;
  const int stall_need = std::max( 1, stall_iterations );
  double prev_gbest = gbest_cost;

  // Reused each iteration: R1/R2 from R's RNG (seedable); gbest broadcast without repmat alloc churn.
  arma::mat R1(n_pop, n_dims);
  arma::mat R2(n_pop, n_dims);
  arma::mat attract_g(n_pop, n_dims);

  // --- Main constriction PSO loop ------------------------------------------
  for (int iter = 0; iter < max_iter; ++iter) {
    // Checked once per iteration; fitness callbacks may also be interrupted from R.
    Rcpp::checkUserInterrupt();
    iterations_done = iter + 1;

    if (show_process) {
      report_pso_progress(iter + 1, gbest_cost);
    }

    // Velocity: v = chi * (v + phi1*U1*(pbest-x) + phi2*U2*(gbest-x))
    // Draw via R::runif so set.seed / optimizerParameters$seed control the swarm.
    for (arma::uword i = 0; i < R1.n_elem; ++i) {
      R1(i) = R::runif(0.0, 1.0);
      R2(i) = R::runif(0.0, 1.0);
    }
    attract_g.each_row() = gbest_pos;
    attract_g -= pos;

    vel = constriction * (vel +
      phi1 * R1 % (pbest_pos - pos) +
      phi2 * R2 % attract_g);
    pos += vel;

    // Feasibility: clamp to windows, then re-sort within each arm group.
    for (arma::uword i = 0; i < n_pop; ++i) {
      for (arma::uword j = 0; j < n_dims; ++j) {
        const arma::mat& bounds = windows[static_cast<size_t>(j)];
        if ( bounds.n_rows == 0 )
          pos(i, j) = initial_pos(j);
        else
          pos(i, j) = clamp_to_windows(pos(i, j), bounds);
      }

      arma::rowvec current_row = pos.row(i);
      sort_groups(current_row, sorting_indices);
      for (arma::uword j = 0; j < n_dims; ++j) {
        if ( windows[static_cast<size_t>(j)].n_rows == 0 )
          current_row(j) = initial_pos(j);
      }
      pos.row(i) = current_row;
    }

    // Evaluate and refresh personal / global bests.
    const arma::vec costs = eval_fitness_rows(eval_fitness_batch, eval_fitness, pos);
    for (arma::uword i = 0; i < n_pop; ++i) {
      if (costs(i) < pbest_cost(i)) {
        pbest_cost(i) = costs(i);
        pbest_pos.row(i) = pos.row(i);
      }
      if (costs(i) < gbest_cost) {
        gbest_cost = costs(i);
        gbest_pos  = pos.row(i);
      }
    }

    // Relative stall stop on D-criterion (same geometry as PGBO).
    // cost = 1/D → D = 1/cost when cost is finite and positive.
    if ( ftol > 0.0 && R_finite( gbest_cost ) && gbest_cost > 0.0 &&
         R_finite( prev_gbest ) && prev_gbest > 0.0 ) {
      const double gbest_d = 1.0 / gbest_cost;
      const double prev_d  = 1.0 / prev_gbest;
      const double rel = std::abs( gbest_d - prev_d ) / std::max( prev_d, 1e-12 );
      if ( rel < ftol )
        ++stall_count;
      else
        stall_count = 0;
      prev_gbest = gbest_cost;
      if ( stall_count >= stall_need ) {
        converged = true;
        break;
      }
    }
  }

  if ( !converged && iterations_done == 0 )
    iterations_done = max_iter;

  const bool improved = R_finite( gbest_cost ) && R_finite( initial_gbest_cost ) &&
    gbest_cost < initial_gbest_cost;

  // converged is only set by the stall stop (ftol > 0); never equate to improved.
  return List::create(
    Named("globalBestDesign") = gbest_pos,
    Named("globalBestCost")   = gbest_cost,
    Named("converged")        = converged,
    Named("iterations")       = iterations_done,
    Named("improved")         = improved
  );
}
