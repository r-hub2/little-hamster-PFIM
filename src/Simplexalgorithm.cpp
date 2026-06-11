// Nelder-Mead simplex on sampling-time parameters (Press et al. amoeba).
// funk(data, named_vector, outcomes) returns 1/D from R; we minimise y.
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <limits>
#include <numeric>
#include <vector>

using namespace Rcpp;

static void fill_named_vector(NumericVector& out, const NumericMatrix& p, int row,
                              const CharacterVector& param_names) {
  for (int j = 0; j < p.ncol(); ++j) out[j] = p(row, j);
  if (param_names.size() == p.ncol()) out.attr("names") = param_names;
}

// ilo = lowest y (best design), ihi = highest (worst vertex).
static void sort_indices_by_y(std::vector<int>& ord, const NumericVector& y) {
  ord.resize(y.size());
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(),
            [&](int a, int b) { return y[a] < y[b]; });
}

// [[Rcpp::export]]
List fun_amoeba_Rcpp(NumericMatrix p,
                     NumericVector y,
                     double ftol,
                     int itmax,
                     Function funk,
                     SEXP data,
                     SEXP outcomes,
                     bool show_process) {
  const double alpha = 1.0;
  const double beta  = 0.5;
  const double gamma = 2.0;
  const double eps   = 1e-10;

  const int mpts  = p.nrow();
  const int ncols = p.ncol();
  CharacterVector param_names = colnames(p);
  const bool has_names = param_names.size() == ncols;

  int iter = 0;
  bool contin = true;
  bool converge = false;

  NumericMatrix results_mat(itmax + 1, 2);
  std::vector<int> ord(mpts);
  NumericVector pbar(ncols);
  NumericVector pr(ncols);
  NumericVector prr(ncols);
  NumericVector row_vec(ncols);
  if (has_names) {
    pr.attr("names") = param_names;
    prr.attr("names") = param_names;
    row_vec.attr("names") = param_names;
  }

  while (contin) {
    const double current_min_y = min(y);
    const double current_criterion = 1.0 / current_min_y;

    if (show_process) {
      Rcout << "iter = " << iter << "\n";
      Rcout << "Criterion = " << current_criterion << "\n";
    }

    results_mat(iter, 0) = iter;
    results_mat(iter, 1) = current_criterion;

    sort_indices_by_y(ord, y);
    const int ilo  = ord.front();
    const int ihi  = ord.back();
    const int inhi = ord[mpts - 2];

    // Relative spread of vertex costs — same stopping rule as legacy PFIM simplex.
    const double rtol = 2.0 * std::abs(y[ihi] - y[ilo]) /
                      (std::abs(y[ihi]) + std::abs(y[ilo]) + eps);

    if (rtol < ftol || iter == itmax) {
      contin = false;
      converge = (iter != itmax);
    } else {
      ++iter;
      Rcpp::checkUserInterrupt();

      for (int j = 0; j < ncols; ++j) {
        double sum = 0.0;
        for (int i = 0; i < mpts; ++i) {
          if (i != ihi) sum += p(i, j);
        }
        pbar[j] = sum / (mpts - 1);
      }

      for (int j = 0; j < ncols; ++j) {
        pr[j] = (1.0 + alpha) * pbar[j] - alpha * p(ihi, j);
      }
      const double ypr = as<double>(funk(data, pr, outcomes));
      const double ypr_safe = (R_IsNaN(ypr) || !R_finite(ypr))
        ? std::numeric_limits<double>::infinity() : ypr;

      if (ypr_safe <= y[ilo]) {
        for (int j = 0; j < ncols; ++j) {
          prr[j] = gamma * pr[j] + (1.0 - gamma) * pbar[j];
        }
        const double yprr = as<double>(funk(data, prr, outcomes));
        const double yprr_safe = (R_IsNaN(yprr) || !R_finite(yprr))
          ? std::numeric_limits<double>::infinity() : yprr;

        if (yprr_safe < y[ilo]) {
          for (int j = 0; j < ncols; ++j) p(ihi, j) = prr[j];
          y[ihi] = yprr_safe;
        } else {
          for (int j = 0; j < ncols; ++j) p(ihi, j) = pr[j];
          y[ihi] = ypr_safe;
        }
      } else {
        if (ypr_safe >= y[inhi]) {
          if (ypr_safe < y[ihi]) {
            for (int j = 0; j < ncols; ++j) p(ihi, j) = pr[j];
            y[ihi] = ypr_safe;
          }

          for (int j = 0; j < ncols; ++j) {
            prr[j] = beta * p(ihi, j) + (1.0 - beta) * pbar[j];
          }
          const double yprr2 = as<double>(funk(data, prr, outcomes));
          const double yprr2_safe = (R_IsNaN(yprr2) || !R_finite(yprr2))
            ? std::numeric_limits<double>::infinity() : yprr2;

          if (yprr2_safe < y[ihi]) {
            for (int j = 0; j < ncols; ++j) p(ihi, j) = prr[j];
            y[ihi] = yprr2_safe;
          } else {
            for (int i = 0; i < mpts; ++i) {
              if (i != ilo) {
                for (int j = 0; j < ncols; ++j) {
                  p(i, j) = 0.5 * (p(i, j) + p(ilo, j));
                }
                fill_named_vector(row_vec, p, i, param_names);
                y[i] = as<double>(funk(data, row_vec, outcomes));
              }
            }
          }
        } else {
          for (int j = 0; j < ncols; ++j) p(ihi, j) = pr[j];
          y[ihi] = ypr_safe;
        }
      }
    }
  }

  const int n_rows = iter + 1;
  NumericVector iter_col(n_rows);
  NumericVector crit_col(n_rows);
  for (int r = 0; r < n_rows; ++r) {
    iter_col[r] = results_mat(n_rows - 1 - r, 0);
    crit_col[r] = results_mat(n_rows - 1 - r, 1);
  }
  DataFrame results = DataFrame::create(
    Named("iter")      = iter_col,
    Named("criterion") = crit_col
  );

  if (!converge && iter == itmax && show_process) {
    Rcout << "Amoeba exceeding maximum iterations.\n";
  }

  return List::create(
    Named("p")        = p,
    Named("y")        = y,
    Named("iter")     = iter,
    Named("converge")  = converge,
    Named("results")  = results
  );
}
