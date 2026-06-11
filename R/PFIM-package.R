#' @name PFIM-package
#' @aliases PFIM
#' @aliases package-PFIM
#' @docType package
#' @encoding UTF-8
#' @title Population Fisher information for design evaluation and optimization in NLME models.
#' @section Description:
#' Nonlinear mixed-effects models (NLMEM) are widely used in model-based drug development.
#' The population Fisher information matrix (FIM) is an efficient alternative to clinical
#' trial simulation for optimizing study designs. **PFIM 7.1** is an R package using the
#' **S7** object system to evaluate and optimize population designs from the FIM.
#'
#' PFIM includes libraries of PK and PD models (S7 classes): bolus, infusion,
#' first-order absorption, one- to three-compartment structures, linear or
#' Michaelis-Menten elimination, direct and turnover PD models, and combined
#' PK/PD models. Users may also supply custom analytical or ODE models.
#' The FIM is computed by first-order linearization with a block-diagonal structure
#' (Mentré et al., 1997). A Bayesian FIM provides shrinkage-related predictions
#' (Combes et al., 2013). Design optimization under the D-criterion uses:
#' the simplex algorithm (Nelder-Mead) (Nelder & Mead, 1965), the
#' multiplicative algorithm (Seurat et al., 2021), the Fedorov-Wynn algorithm (Fedorov, 1972), PSO (*Particle Swarm Optimization*) and PGBO (*Population Genetics Based
#' Optimizer*) (Le Nagard et al., 2011).
#'
#' @section Documentation:
#' Package source, issues, and vignettes are available at
#' \url{https://github.com/packagePFIM/PFIM}
#'
#' @section Session options:
#' Runtime toggles live in an internal session environment. Use
#' \code{pfim_set_option(fim.cache = FALSE)} and \code{pfim_get_option("fim.cache")}.
#' \describe{
#'   \item{\code{fim.cache}}{Reuse evaluated designs during optimization (default \code{TRUE}).}
#'   \item{\code{fim.cache.hits}}{Cache hit counter (reset at each constraint enumeration).}
#'   \item{\code{fim.cache.scope}}{Active cache namespace (set when optimization starts).}
#'   \item{\code{constraints.maxTasks}}{Cap dose x sampling cells on the constraint grid.}
#'   \item{\code{verbose}}{Progress messages during long FIM enumerations.}
#'   \item{\code{eval.batch}}{Batch evaluation mode during constraint-grid enumeration.}
#'   \item{\code{perf.adminCache}, \code{perf.fdCache}, \code{perf.odeTimesCache}}{Gradient
#'     performance caches for ODE models (default \code{TRUE}).}
#'   \item{\code{model.cache.signature}}{Content-based key for \code{rebuildEvalModel}.}
#'   \item{\code{devPath}}{Package root when loading from source (\code{.onLoad}).}
#' }
#'
#' @section Validation:
#' Results are shown via \code{show()}, plots (\code{ggplot2}), and HTML reports.
#' Covariates, inter-occasion variability (IOV), and occasion covariates are supported
#' for population FIM; the number of occasions is inferred automatically on the model.
#'
#' @references
#' Dumont C, Lestini G, Le Nagard H, Mentré F, Comets E, Nguyen TT, et al. PFIM 4.0, an extended R program for design evaluation and optimization in
#' nonlinear mixed-effect models. Comput Methods Programs Biomed. 2018;156:217-29.
#'
#' Chambers JM. Object-Oriented Programming, Functional Programming and R. Stat Sci. 2014;29:167-80.
#'
#' Mentré F, Mallet A, Baccar D. Optimal Design in Random-Effects Regression Models. Biometrika. 1997;84:429-42.
#'
#' Combes FP, Retout S, Frey N, Mentré F. Prediction of shrinkage of individual parameters using the Bayesian information matrix in nonlinear mixed effect models with evaluation in pharmacokinetics. Pharm Res. 2013;30:2355-67.
#'
#' Nelder JA, Mead R. A simplex method for function minimization. Comput J. 1965;7:308-13.
#'
#' Seurat J, Tang Y, Mentré F, Nguyen, TT. Finding optimal design in nonlinear mixed effect models using multiplicative algorithms. Computer Methods and Programs in Biomedicine, 2021.
#'
#' Fedorov VV. Theory of Optimal Experiments. Academic Press, New York, 1972.
#'
#' Eberhart RC, Kennedy J. A new optimizer using particle swarm theory. Proc. of the Sixth International Symposium on Micro Machine and Human Science, Nagoya, 4-6 October 1995, 39-43.
#'
#' Le Nagard H, Chao L, Tenaillon O. The emergence of complexity and restricted pleiotropy in adapting networks. BMC Evol Biol. 2011;11:326.
#'
#' Wickham H. ggplot2: Elegant Graphics for Data Analysis, Springer-Verlag New York, 2016.
#'
#' @importFrom Matrix bdiag
#' @importFrom rlang obj_address
#'
# @section Content of the source code and files in the \code{/R} folder:
#'
"_PACKAGE"
