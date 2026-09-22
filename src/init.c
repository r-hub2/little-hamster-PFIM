/* Native routine registration. Keep in sync with R/RcppExports.R.
   Rcpp::compileAttributes() does not write this file. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP _PFIM_FedorovWynnAlgorithm_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PFIM_MultiplicativeAlgorithm_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PFIM_pgbo_optimize_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PFIM_pso_optimize_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PFIM_fun_amoeba_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PFIM_chol_inv_Rcpp(SEXP);
extern SEXP _PFIM_safe_solve_Rcpp(SEXP);
extern SEXP _PFIM_computeMFVar_Rcpp(SEXP, SEXP);
extern SEXP _PFIM_computeMFVar_mixed_Rcpp(SEXP, SEXP, SEXP);
extern SEXP _PFIM_computePopFimCombo_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PFIM_residualErrorDerivatives_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PFIM_pfimCacheHash_Rcpp(SEXP);
extern SEXP _PFIM_pfimEnvCacheKey_Rcpp(SEXP);
extern SEXP _PFIM_pfimOdeSimTimesCached_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PFIM_pfimOdeSimTimesCacheClear_Rcpp(void);
extern SEXP _PFIM_pfimOdeSimTimesCacheSize_Rcpp(void);

static const R_CallMethodDef CallEntries[] = {
    {"_PFIM_FedorovWynnAlgorithm_Rcpp", (DL_FUNC) &_PFIM_FedorovWynnAlgorithm_Rcpp, 13},
    {"_PFIM_MultiplicativeAlgorithm_Rcpp", (DL_FUNC) &_PFIM_MultiplicativeAlgorithm_Rcpp, 8},
    {"_PFIM_pgbo_optimize_Rcpp", (DL_FUNC) &_PFIM_pgbo_optimize_Rcpp, 14},
    {"_PFIM_pso_optimize_Rcpp", (DL_FUNC) &_PFIM_pso_optimize_Rcpp, 14},
    {"_PFIM_fun_amoeba_Rcpp", (DL_FUNC) &_PFIM_fun_amoeba_Rcpp, 8},
    {"_PFIM_chol_inv_Rcpp", (DL_FUNC) &_PFIM_chol_inv_Rcpp, 1},
    {"_PFIM_safe_solve_Rcpp", (DL_FUNC) &_PFIM_safe_solve_Rcpp, 1},
    {"_PFIM_computeMFVar_Rcpp", (DL_FUNC) &_PFIM_computeMFVar_Rcpp, 2},
    {"_PFIM_computeMFVar_mixed_Rcpp", (DL_FUNC) &_PFIM_computeMFVar_mixed_Rcpp, 3},
    {"_PFIM_computePopFimCombo_Rcpp", (DL_FUNC) &_PFIM_computePopFimCombo_Rcpp, 8},
    {"_PFIM_residualErrorDerivatives_Rcpp", (DL_FUNC) &_PFIM_residualErrorDerivatives_Rcpp, 7},
    {"_PFIM_pfimCacheHash_Rcpp", (DL_FUNC) &_PFIM_pfimCacheHash_Rcpp, 1},
    {"_PFIM_pfimEnvCacheKey_Rcpp", (DL_FUNC) &_PFIM_pfimEnvCacheKey_Rcpp, 1},
    {"_PFIM_pfimOdeSimTimesCached_Rcpp", (DL_FUNC) &_PFIM_pfimOdeSimTimesCached_Rcpp, 5},
    {"_PFIM_pfimOdeSimTimesCacheClear_Rcpp", (DL_FUNC) &_PFIM_pfimOdeSimTimesCacheClear_Rcpp, 0},
    {"_PFIM_pfimOdeSimTimesCacheSize_Rcpp", (DL_FUNC) &_PFIM_pfimOdeSimTimesCacheSize_Rcpp, 0},
    {NULL, NULL, 0}
};

void R_init_PFIM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

