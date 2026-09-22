# Package-level Rd documentation (`_PACKAGE`). Describes PFIM purpose, session
# options, validation, authors, and Imports / useDynLib declarations.

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
#' first-order absorption, one- and two-compartment structures, linear or
#' Michaelis-Menten elimination, direct and turnover PD models, and combined
#' PK/PD models. Users may also supply custom analytical or ODE models; set
#' \code{modelClass} explicitly or use \code{\link{pfim_resolve_model_class}}.
#' The FIM is computed by first-order linearization with a block-diagonal
#' structure. A Bayesian FIM provides shrinkage-related predictions.
#' Design optimization under the D-criterion uses the simplex algorithm
#' (Nelder-Mead), the multiplicative algorithm, Fedorov-Wynn, PSO
#' (Particle Swarm Optimization) and PGBO (Population Genetics Based Optimizer).
#'
#' @section User API:
#' Constructors: \code{\link{Evaluation}}, \code{\link{Optimization}},
#' \code{\link{Arm}}, \code{\link{Design}}, \code{\link{ModelParameter}},
#' residual-error and distribution classes, FIM types, and algorithm classes.
#' Verbs: \code{\link{run}}, \code{\link{Report}}, \code{\link{defineFim}},
#' \code{\link{definePKModel}}, \code{\link{definePKPDModel}},
#' \code{getSE}/\code{getRSE}/\code{getFisherMatrix}/\code{getDcriterion}/\code{getShrinkage},
#' \code{plotEvaluation}/\code{plotSE}/\code{plotRSE}/\code{plotSensitivityIndices}/\code{plotWeights}/\code{plotFrequencies},
#' \code{\link{showFIM}}, \code{\link{covariateTest}}.
#' Session: \code{\link{pfim_set_option}}, \code{\link{pfim_reset_session}}.
#' Pipeline helpers (\code{evaluateDesign}, \code{optimizeDesign}, \ldots) are
#' internal; they remain callable from tests via the package namespace.
#'
#' @section Documentation:
#' Package source, issues, and vignettes are available at
#' \url{https://github.com/packagePFIM/PFIM}
#'
#' @section Session options:
#' \code{\link{pfim_set_option}} / \code{\link{pfim_get_option}}; \code{\link{pfim_reset_session}} clears caches.
#' See \code{?pfim_set_option} for \code{fim.cache}, \code{constraints.maxTasks}, and performance options.
#'
#' @section Validation:
#' Results are shown via \code{show()}, plots (\code{ggplot2}), and HTML reports.
#' Covariates, inter-occasion variability (IOV), and occasion covariates are supported
#' for population FIM; the number of occasions is inferred automatically on the model.
#'
#' @references
#' Dumont C, Lestini G, Le Nagard H, Comets E, Nguyen TT, et al. PFIM 4.0, an extended R program for design evaluation and optimization in
#' nonlinear mixed-effect models. Comput Methods Programs Biomed. 2018;156:217-29.
#'
#' Chambers JM. Object-Oriented Programming, Functional Programming and R. Stat Sci. 2014;29:167-80.
#'
#' Nelder JA, Mead R. A simplex method for function minimization. Comput J. 1965;7:308-13.
#'
#' Seurat J, Tang Y, Nguyen TT. Finding optimal design in nonlinear mixed effect models using multiplicative algorithms. Computer Methods and Programs in Biomedicine, 2021.
#'
#' Fedorov VV. Theory of Optimal Experiments. Academic Press, New York, 1972.
#'
#' Eberhart RC, Kennedy J. A new optimizer using particle swarm theory. Proc. of the Sixth International Symposium on Micro Machine and Human Science, Nagoya, 4-6 October 1995, 39-43.
#'
#' Le Nagard H, Chao L, Tenaillon O. The emergence of complexity and restricted pleiotropy in adapting networks. BMC Evol Biol. 2011;11:326.
#'
#' Wickham H. ggplot2: Elegant Graphics for Data Analysis, Springer-Verlag New York, 2016.
#'
#' @import S7
#' @import Rcpp
#' @importFrom methods show
#' @importFrom deSolve ode
#' @importFrom purrr accumulate compact detect_index flatten imap iwalk keep list_c list_flatten list_rbind map map2 map2_chr map_chr map_dbl map_int map_lgl pluck pmap reduce reduce2 set_names walk walk2
#' @importFrom stringr fixed regex str_c str_detect str_remove str_replace str_replace_all str_split str_starts
#' @importFrom ggplot2 .data aes coord_flip element_blank element_line element_rect element_text expansion facet_wrap geom_bar geom_col geom_line geom_point ggplot labs margin scale_x_continuous scale_x_discrete scale_y_continuous sec_axis theme theme_bw theme_grey theme_minimal unit
#' @importFrom rmarkdown render
#' @importFrom tibble as_tibble
#' @importFrom utils combn tail head flush.console capture.output
#' @importFrom Matrix bdiag
#' @importFrom Deriv Simplify
#' @importFrom scales pretty_breaks
#' @importFrom kableExtra kable_styling kbl add_header_above footnote
#' @importFrom knitr asis_output
#' @importFrom rlang hash `%||%`
#' @importFrom stats D cov2cor pnorm qnorm reorder uniroot
#' @useDynLib PFIM, .registration = TRUE
#'
#' @author
#' \strong{Maintainer}: France Mentré \email{PFIM@inserm.fr}
#' (\href{https://orcid.org/0000-0002-7045-1275}{ORCID}).
#'
#' \strong{Author, creator of PFIM 7.1}: Romain Leroux \email{romain.leroux@inserm.fr}
#' (\href{https://orcid.org/0009-0009-5779-5303}{ORCID}) (PFIM, 2020--2026).
#'
#' \strong{Other contributors}:
#' \strong{Contributor} Jérémy Seurat \email{jeremy.seurat@prive.fr},
#' \strong{Contributor} Antoine Croxo \email{antoine.croxo@inserm.fr},
#' \strong{Contributor} Lucie Fayette \email{lucie.fayette@inserm.fr}.
#'
"_PACKAGE"
