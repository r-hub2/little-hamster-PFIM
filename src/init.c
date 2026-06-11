/* .Call registration for PFIM C++ optimizers (paired with R/RcppExports.R). */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP _PFIM_FedorovWynnAlgorithm_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PFIM_MultiplicativeAlgorithm_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PFIM_fun_amoeba_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PFIM_pso_optimize_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PFIM_pgbo_optimize_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_PFIM_FedorovWynnAlgorithm_Rcpp", (DL_FUNC) &_PFIM_FedorovWynnAlgorithm_Rcpp, 12},
    {"_PFIM_MultiplicativeAlgorithm_Rcpp", (DL_FUNC) &_PFIM_MultiplicativeAlgorithm_Rcpp, 8},
    {"_PFIM_fun_amoeba_Rcpp", (DL_FUNC) &_PFIM_fun_amoeba_Rcpp, 8},
    {"_PFIM_pso_optimize_Rcpp", (DL_FUNC) &_PFIM_pso_optimize_Rcpp, 12},
    {"_PFIM_pgbo_optimize_Rcpp", (DL_FUNC) &_PFIM_pgbo_optimize_Rcpp, 11},
    {NULL, NULL, 0}
};

void R_init_PFIM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
