// PSO on flat sampling-time vectors (Clerc & Kennedy constriction PSO).
// R side supplies eval_fitness / eval_fitness_batch: cost = 1/D, lower is better.
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <vector>

using namespace Rcpp;

// Project particle coordinates onto the nearest admissible window (per dimension).
static double clamp_to_windows(double val, const arma::mat& bounds) {
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

// Within-arm sampling times must stay ordered; indices are 1-based in R input.
static void sort_groups(arma::rowvec& row, const std::vector<arma::uvec>& sorting_indices) {
  for (size_t g = 0; g < sorting_indices.size(); ++g) {
    const arma::uvec& idx = sorting_indices[g];
    arma::rowvec sub_pos(idx.n_elem);
    for (arma::uword k = 0; k < idx.n_elem; ++k) sub_pos(k) = row(idx(k));
    sub_pos = arma::sort(sub_pos);
    for (arma::uword k = 0; k < idx.n_elem; ++k) row(idx(k)) = sub_pos(k);
  }
}

static void report_pso_progress(int iteration, double gbest_cost) {
  Rcout << "Iteration = " << iteration;
  if (R_finite(gbest_cost) && gbest_cost > 0.0 && gbest_cost < 1e6) {
    Rcout << " | Criterion = " << (1.0 / gbest_cost) << "\n";
  } else {
    Rcout << " | Criterion = N/A (no valid design yet)\n";
  }
}

// Prefer eval_fitness_batch (one R call per swarm state) when provided;
// otherwise fall back to per-row scalar callbacks.
static arma::vec eval_fitness_rows(const Nullable<Function>& eval_batch,
                                   const Function& eval_scalar,
                                   const arma::mat& positions) {
  if (eval_batch.isNotNull()) {
    Function batch_fn = eval_batch.get();
    return as<arma::vec>(batch_fn(wrap(positions)));
  }
  arma::vec costs(positions.n_rows);
  for (arma::uword i = 0; i < positions.n_rows; ++i) {
    costs(i) = as<double>(eval_scalar(wrap(positions.row(i))));
  }
  return costs;
}

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
                       Nullable<Function> eval_fitness_batch = R_NilValue) {
  const arma::uword n_pop  = static_cast<arma::uword>(n_pop_in);
  const arma::uword n_dims = initial_pos.n_elem;
  const bool use_valid_sampler = sample_valid_pos.isNotNull();

  std::vector<arma::mat> windows(static_cast<size_t>(n_dims));
  for (arma::uword j = 0; j < n_dims; ++j)
    windows[static_cast<size_t>(j)] = as<arma::mat>(windows_list[static_cast<int>(j)]);

  std::vector<arma::uvec> sorting_indices(static_cast<size_t>(sorting_groups.size()));
  for (int g = 0; g < sorting_groups.size(); ++g)
    sorting_indices[static_cast<size_t>(g)] = as<arma::uvec>(sorting_groups[g]) - 1;

  arma::mat pos = arma::repmat(arma::mat(initial_pos), n_pop, 1);
  arma::mat vel = arma::zeros<arma::mat>(n_pop, n_dims);
  arma::mat pbest_pos = pos;
  arma::vec pbest_cost(n_pop);
  pbest_cost.fill(arma::datum::inf);

  arma::rowvec gbest_pos = initial_pos;
  double gbest_cost = arma::datum::inf;

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

  for (arma::uword i = 0; i < n_pop; ++i) {
    if (use_valid_sampler) {
      Function sampler = sample_valid_pos.get();
      pos.row(i) = as<arma::rowvec>(sampler());
    } else {
      for (arma::uword j = 0; j < n_dims; ++j) {
        const arma::mat& bounds = windows[static_cast<size_t>(j)];
        const double w_min = bounds.col(0).min();
        const double w_max = bounds.col(1).max();
        pos(i, j) = w_min + R::runif(0.0, 1.0) * (w_max - w_min);
        pos(i, j) = clamp_to_windows(pos(i, j), bounds);
      }
    }

    arma::rowvec current_row = pos.row(i);
    sort_groups(current_row, sorting_indices);
    pos.row(i) = current_row;
  }

  {
    const arma::vec init_costs = eval_fitness_rows(eval_fitness_batch, eval_fitness, pos);
    pbest_cost = init_costs;
    pbest_pos  = pos;
    update_best_from_costs(init_costs, pos, pbest_pos, pbest_cost, gbest_pos, gbest_cost);
  }

  for (int iter = 0; iter < max_iter; ++iter) {
    // Checked once per iteration; fitness callbacks may also be interrupted from R.
    Rcpp::checkUserInterrupt();

    if (show_process) {
      report_pso_progress(iter + 1, gbest_cost);
    }

    const arma::mat R1 = arma::randu<arma::mat>(n_pop, n_dims);
    const arma::mat R2 = arma::randu<arma::mat>(n_pop, n_dims);
    const arma::mat gbest_mat = arma::repmat(arma::mat(gbest_pos), n_pop, 1);

    vel = constriction * (vel +
      phi1 * R1 % (pbest_pos - pos) +
      phi2 * R2 % (gbest_mat - pos));
    pos += vel;

    for (arma::uword i = 0; i < n_pop; ++i) {
      for (arma::uword j = 0; j < n_dims; ++j) {
        pos(i, j) = clamp_to_windows(pos(i, j), windows[static_cast<size_t>(j)]);
      }

      arma::rowvec current_row = pos.row(i);
      sort_groups(current_row, sorting_indices);
      pos.row(i) = current_row;
    }

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
  }

  return List::create(
    Named("globalBestDesign") = gbest_pos,
    Named("globalBestCost")   = gbest_cost
  );
}
