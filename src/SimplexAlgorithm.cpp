// Nelder–Mead simplex ("amoeba") on flat sampling-time parameters.
// Reference: Press et al., Numerical Recipes — classic reflection / expansion /
// contraction / shrink geometry.
//
// Contract with R:
//   - vertices live in rows of `p`; `y[i]` = cost of vertex i
//   - cost = funk(data, named_vector, outcomes) = 1/D  (lower is better)
//   - we minimise y; reported criterion / history is D = 1/min(y)
//   - column names on `p` (when present) are forwarded to `funk`
//   - no window projection here — R's funk / constraint layer owns feasibility
//   - non-finite funk returns map to +Inf so ranking and rtol stay defined
//
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <limits>
#include <numeric>
#include <vector>

using namespace Rcpp;

// --- Helpers ---------------------------------------------------------------

// Copy one simplex vertex into `out`, preserving column names when present.
static void fill_named_vector(NumericVector& out, const NumericMatrix& p, int row,
                              const CharacterVector& param_names) {
  for (int j = 0; j < p.ncol(); ++j) out[j] = p(row, j);
  if (param_names.size() == p.ncol()) out.attr("names") = param_names;
}

// Sort vertex indices by ascending cost: ilo = best (lowest y), ihi = worst.
static void sort_indices_by_y(std::vector<int>& ord, const NumericVector& y) {
  ord.resize(y.size());
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(),
            [&](int a, int b) { return y[a] < y[b]; });
}

// Map non-finite costs to +Inf so ranking / rtol stay well-defined.
static double finite_or_inf(double y) {
  return (R_IsNaN(y) || !R_finite(y))
    ? std::numeric_limits<double>::infinity() : y;
}

// --- Exported Nelder–Mead loop ---------------------------------------------
// [[Rcpp::export]]
List fun_amoeba_Rcpp(NumericMatrix p,
                     NumericVector y,
                     double ftol,
                     int itmax,
                     Function funk,
                     SEXP data,
                     SEXP outcomes,
                     bool show_process) {
  // Standard Nelder–Mead move coefficients.
  const double alpha = 1.0;   // reflection
  const double beta  = 0.5;   // contraction
  const double gamma = 2.0;   // expansion
  const double eps   = 1e-10; // floor in relative tolerance

  const int mpts  = p.nrow(); // number of simplex vertices (dim + 1 typically)
  const int ncols = p.ncol(); // flat sampling-time dimension
  if (mpts < 2)
    stop("Simplex: at least 2 simplex vertices are required.");
  if (y.size() != mpts)
    stop("Simplex: length(y) must equal nrow(p).");
  if (itmax < 0)
    stop("Simplex: itmax must be non-negative.");
  CharacterVector param_names = colnames(p);
  const bool has_names = param_names.size() == ncols;

  int iter = 0;
  bool contin = true;
  bool converge = false;

  // Iteration history: (iter, D-criterion) for the R report.
  NumericMatrix results_mat(itmax + 1, 2);
  std::vector<int> ord(mpts);
  NumericVector pbar(ncols);   // centroid of all but the worst vertex
  NumericVector pr(ncols);     // reflected trial point
  NumericVector prr(ncols);    // expanded / contracted trial point
  NumericVector row_vec(ncols);
  if (has_names) {
    pr.attr("names") = param_names;
    prr.attr("names") = param_names;
    row_vec.attr("names") = param_names;
  }

  while (contin) {
    // --- Progress / history ------------------------------------------------
    const double current_min_y = min(y);
    const double current_criterion = 1.0 / current_min_y;

    if (show_process) {
      Rcout << "iter = " << iter << "\n";
      Rcout << "Criterion = " << current_criterion << "\n";
    }

    results_mat(iter, 0) = iter;
    results_mat(iter, 1) = current_criterion;

    // --- Rank vertices -----------------------------------------------------
    sort_indices_by_y(ord, y);
    const int ilo  = ord.front();       // best
    const int ihi  = ord.back();        // worst
    const int inhi = ord[mpts - 2];     // second-worst

    // Relative spread of vertex costs — same stopping rule as legacy PFIM simplex.
    const double rtol = 2.0 * std::abs(y[ihi] - y[ilo]) /
                      (std::abs(y[ihi]) + std::abs(y[ilo]) + eps);

    if (rtol < ftol || iter == itmax) {
      contin = false;
      converge = (rtol < ftol);
    } else {
      ++iter;
      Rcpp::checkUserInterrupt();

      // --- Centroid of all vertices except the worst -----------------------
      for (int j = 0; j < ncols; ++j) {
        double sum = 0.0;
        for (int i = 0; i < mpts; ++i) {
          if (i != ihi) sum += p(i, j);
        }
        pbar[j] = sum / (mpts - 1);
      }

      // --- Reflect the worst point through the centroid --------------------
      // pr = (1+alpha)*pbar - alpha*p_worst
      for (int j = 0; j < ncols; ++j) {
        pr[j] = (1.0 + alpha) * pbar[j] - alpha * p(ihi, j);
      }
      const double ypr = as<double>(funk(data, pr, outcomes));
      const double ypr_safe = finite_or_inf(ypr);

      // --- Accept / expand / contract / shrink -----------------------------
      if (ypr_safe <= y[ilo]) {
        // Reflection is a new best: try expansion beyond pr.
        for (int j = 0; j < ncols; ++j) {
          prr[j] = gamma * pr[j] + (1.0 - gamma) * pbar[j];
        }
        const double yprr = as<double>(funk(data, prr, outcomes));
        const double yprr_safe = finite_or_inf(yprr);

        if (yprr_safe < y[ilo]) {
          for (int j = 0; j < ncols; ++j) p(ihi, j) = prr[j];
          y[ihi] = yprr_safe;
        } else {
          for (int j = 0; j < ncols; ++j) p(ihi, j) = pr[j];
          y[ihi] = ypr_safe;
        }
      } else {
        if (ypr_safe >= y[inhi]) {
          // Reflection is not better than second-worst: try contraction.
          if (ypr_safe < y[ihi]) {
            // Keep the reflected point if it beats the current worst.
            for (int j = 0; j < ncols; ++j) p(ihi, j) = pr[j];
            y[ihi] = ypr_safe;
          }

          // Contract toward the centroid from the (possibly updated) worst.
          for (int j = 0; j < ncols; ++j) {
            prr[j] = beta * p(ihi, j) + (1.0 - beta) * pbar[j];
          }
          const double yprr2 = as<double>(funk(data, prr, outcomes));
          const double yprr2_safe = finite_or_inf(yprr2);

          if (yprr2_safe < y[ihi]) {
            for (int j = 0; j < ncols; ++j) p(ihi, j) = prr[j];
            y[ihi] = yprr2_safe;
          } else {
            // Shrink all non-best vertices halfway toward the best.
            for (int i = 0; i < mpts; ++i) {
              if (i != ilo) {
                for (int j = 0; j < ncols; ++j) {
                  p(i, j) = 0.5 * (p(i, j) + p(ilo, j));
                }
                fill_named_vector(row_vec, p, i, param_names);
                y[i] = finite_or_inf(as<double>(funk(data, row_vec, outcomes)));
              }
            }
          }
        } else {
          // Intermediate reflection: replace the worst vertex with pr.
          for (int j = 0; j < ncols; ++j) p(ihi, j) = pr[j];
          y[ihi] = ypr_safe;
        }
      }
    }
  }

  // --- Pack reverse-chronology history for the R report --------------------
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
    Named("p")          = p,
    Named("y")          = y,
    Named("iterations") = iter,
    Named("converged")  = converge,
    Named("results")    = results
  );
}
