#' CovariateTest Class
#'
#' @description
#' \code{CovariateTest} stores covariate-test results.
#'
#' @details
#' \describe{
#'   \item{covariate}{Significance on covariate effects \eqn{\beta}}
#'   \item{nonRelevance}{TOST non-relevance on \eqn{\beta}}
#'   \item{relevance}{Clinical relevance on \eqn{\beta}}
#'   \item{settings}{Run options (\code{tests}, \code{thetaL/U}, \code{target_power}, \code{alpha}, \code{N0})}
#' }
#'
#' @param covariate    \code{data.frame} for \eqn{\beta} significance.
#' @param parameter    Deprecated; kept empty (typical \eqn{\mu} are not reported).
#' @param nonRelevance \code{data.frame} for TOST non-relevance on \eqn{\beta}.
#' @param relevance    \code{data.frame} for clinical relevance on \eqn{\beta}.
#' @param settings     List of run/display options.
#'
#' @return An S7 object of class \code{CovariateTest}.
#' @examples
#' \donttest{
#' source(system.file("examples", "covariate-test-minimal.R", package = "PFIM"))
#' class(ct)
#' }
#'
#' @name CovariateTest-class
#' @rdname CovariateTest-class
#' @aliases CovariateTest
#' @usage NULL
#' @export
CovariateTest = new_class( "CovariateTest", package = "PFIM",
                           properties = list(
                             covariate    = new_property( class_any, default = data.frame() ),
                             parameter    = new_property( class_any, default = data.frame() ),
                             nonRelevance = new_property( class_any, default = data.frame() ),
                             relevance    = new_property( class_any, default = data.frame() ),
                             settings     = new_property( class_any, default = list() )
                           )
)
S4_register( CovariateTest )

#' covariateTest
#'
#' @description
#' Assess significance and clinical relevance from the evaluation FIM.
#'
#' @details
#' Covariate effects \eqn{\beta} only. Use \code{tests} to select which tables
#' are computed and displayed (\code{"significance"}, \code{"nonRelevance"},
#' \code{"relevance"}). Requires a \strong{population} FIM: individual and
#' Bayesian FIMs omit \eqn{\beta} (subject-level precision conditions on
#' covariates).
#'
#' Relative-effect display (\code{Ratio}, TOST IC) uses \eqn{e^{\beta}} for the
#' exponential covariate link and \eqn{1+\beta} for the additive link. The
#' clinical window \code{[thetaL, thetaU]} remains on the log-ratio scale
#' (default bioequivalence \code{[log(0.80), log(1.25)]}).
#'
#' Sample-size scaling: population FIM recovers unit variance as
#' \code{SE^2 * N0}; individual / Bayesian FIM (not multiplied by arm size)
#' use \code{N_scale = 1}, so \code{N_Required} is the number of independent
#' replications of the subject-level design.
#'
#' Output column order (PFIM4-compatible):
#' \code{Value}, \code{SE}, \code{RSE}, then (TOST) \code{Ratio} and IC bounds,
#' then \code{Power}, \code{N_Required}. Each table is preceded by the target
#' power for \code{N_Required} and, for TOST tables, the equivalence IC bounds.
#'
#' @param pfimproject   An \code{Evaluation} object (after \code{run()}).
#' @param thetaL        Lower bound on the log-ratio scale (default: \code{log(0.80)}).
#' @param thetaU        Upper bound on the log-ratio scale (default: \code{log(1.25)}).
#' @param target_power  Target power (default: 0.90).
#' @param alpha         Nominal type-I error rate (default: 0.05).
#' @param tests         Character vector: \code{"significance"},
#'   \code{"nonRelevance"}, \code{"relevance"} (subset allowed).
#'
#' @return A \code{CovariateTest} object.
#'
#' @examples
#' \donttest{
#' source(system.file("examples", "covariate-test-minimal.R", package = "PFIM"))
#' nrow(prop(ct, "covariate"))
#' ct_nr = covariateTest(ev, tests = "nonRelevance")
#' }
#'
#' @name covariateTest
#' @rdname covariateTest
#' @aliases covariateTest tost
#' @include Evaluation.R
#' @export
covariateTest = new_generic( "covariateTest", "pfimproject",
                             fun = function( pfimproject,
                                             thetaL       = log( 0.80 ),
                                             thetaU       = log( 1.25 ),
                                             target_power = 0.90,
                                             alpha        = 0.05,
                                             tests        = c(
                                               "significance", "nonRelevance", "relevance"
                                             ) ) {
                               S7_dispatch()
                             }
)

# Allowed `tests=` tokens (print / HTML / file writers share this set).
.VALID_COVARIATE_TESTS = c( "significance", "nonRelevance", "relevance" )

#' Reject unknown \code{tests} tokens early (shared by run + display helpers).
#' @noRd
#' @keywords internal
.normalizeCovariateTests = function( tests ) {
  tests = unique( as.character( tests ) )
  bad   = setdiff( tests, .VALID_COVARIATE_TESTS )
  if ( length( bad ) )
    stop(
      "Invalid tests: ", paste( bad, collapse = ", " ),
      ". Choose from: ", paste( .VALID_COVARIATE_TESTS, collapse = ", " ),
      call. = FALSE
    )
  tests
}

#' Ceiling for finite N; keep NA / Inf as-is for “unreachable power” cases.
#' @noRd
#' @keywords internal
.formatNRequired = function( N_req ) {
  if ( is.na( N_req ) ) return( NA_real_ )
  if ( is.infinite( N_req ) ) return( Inf )
  ceiling( N_req )
}

#' Banner line for console / file tables (power, N0, optional TOST window).
#' @noRd
#' @keywords internal
.covariateTestHeader = function( kind, settings ) {
  # Empty CovariateTest() has settings = list(); avoid NULL arithmetic / exp().
  tp = settings$target_power %||% NA_real_
  N0 = settings$N0 %||% NA_real_
  hdr = sprintf(
    "Target power: %.0f%% | Current sample size N = %.0f",
    as.numeric( tp ) * 100, as.numeric( N0 )
  )
  if ( kind %in% c( "nonRelevance", "relevance" ) ) {
    # TOST uses one-sided alpha per bound => displayed CI level is 1 - 2*alpha
    # (e.g. alpha = 0.05 -> 90% IC). Bounds already use z_one = qnorm(1 - alpha).
    alpha  = settings$alpha  %||% NA_real_
    thetaL = settings$thetaL %||% NA_real_
    thetaU = settings$thetaU %||% NA_real_
    ic_pct = ( 1 - 2 * as.numeric( alpha ) ) * 100
    hdr = paste0(
      hdr,
      sprintf(
        " | Equivalence IC (%.0f%%) on ratio: [%.2f, %.2f]",
        ic_pct, exp( as.numeric( thetaL ) ), exp( as.numeric( thetaU ) )
      )
    )
  }
  hdr
}

#' PFIM4 column order for display; drop missing cols when a table is empty/partial.
#' @noRd
#' @keywords internal
.covariateTestDisplayDf = function( df, kind ) {
  if ( is.null( df ) || !nrow( df ) ) return( df )
  if ( kind == "significance" ) {
    cols = c( "Parameter", "Value", "SE", "RSE", "Power", "N_Required" )
  } else {
    cols = c(
      "Parameter", "Value", "SE", "RSE", "Ratio",
      "IC_Inf", "IC_Sup", "Power", "N_Required"
    )
  }
  df[ , intersect( cols, names( df ) ), drop = FALSE ]
}

#' Human titles for the three beta-test blocks (print / HTML / file).
#' @noRd
#' @keywords internal
.covariateTestSlotTitle = function( kind ) {
  switch( kind,
    significance = "Statistical significance (\u03b2, bilateral Wald test)",
    nonRelevance = "Clinical non-relevance (TOST on \u03b2)",
    relevance    = "Clinical relevance (\u03b2 outside equivalence bounds)"
  )
}

#' Prefer explicit \code{tests=} override; else settings from the last run.
#' @noRd
#' @keywords internal
.covariateTestResolvedTests = function( covariateTestResult, tests = NULL ) {
  if ( !is.null( tests ) ) return( .normalizeCovariateTests( tests ) )
  st = prop( covariateTestResult, "settings" )$tests
  if ( length( st ) ) return( .normalizeCovariateTests( st ) )
  .VALID_COVARIATE_TESTS
}

#' Map display kind -> S7 slot (\code{significance} is stored as \code{covariate}).
#' @noRd
#' @keywords internal
.covariateTestSlotDf = function( covariateTestResult, kind ) {
  slot_name = switch( kind,
    significance = "covariate",
    nonRelevance = "nonRelevance",
    relevance    = "relevance"
  )
  .covariateTestDisplayDf( prop( covariateTestResult, slot_name ), kind )
}

#' @noRd
#' @keywords internal
.printCovariateTestSlot = function( object, kind ) {
  settings = prop( object, "settings" )
  df       = .covariateTestSlotDf( object, kind )
  cat( "\n", strrep( "=", 60 ), "\n", sep = "" )
  cat( " ", .covariateTestSlotTitle( kind ), "\n", sep = "" )
  cat( " ", .covariateTestHeader( kind, settings ), "\n", sep = "" )
  cat( strrep( "=", 60 ), "\n\n", sep = "" )
  if ( !nrow( df ) ) {
    cat( "  (none)\n" )
  } else {
    print( df, row.names = FALSE )
  }
}

#' @noRd
#' @keywords internal
.formatCovariateTestSlotLines = function( object, kind ) {
  settings = prop( object, "settings" )
  df       = .covariateTestSlotDf( object, kind )
  c(
    strrep( "=", 60 ),
    paste0( "  ", .covariateTestSlotTitle( kind ) ),
    paste0( "  ", .covariateTestHeader( kind, settings ) ),
    strrep( "=", 60 ),
    "",
    if ( !nrow( df ) ) {
      "  (none)"
    } else {
      utils::capture.output( print( df, row.names = FALSE ) )
    },
    ""
  )
}

#' Two-sided Wald power for testing \eqn{H_0: \beta = 0}.
#' @noRd
#' @keywords internal
.powerSignificance = function( beta, SE, z_half ) {
  1 - pnorm( z_half - beta / SE ) + pnorm( -z_half - beta / SE )
}

#' Sample size so that Wald power for \eqn{\beta} reaches target \code{PS}.
#'
#' Inverts the critical-value equation for the SE that yields power \code{PS},
#' then converts with \eqn{N = \sigma^2_{\mathrm{unit}} / \mathrm{SE}^2}.
#' @noRd
#' @keywords internal
.nRequiredSignificance = function( beta, sigma2_unit, z_half, PS ) {
  if ( beta == 0 ) return( NA_real_ )
  seS = if ( beta > 0 ) {
    beta  / ( z_half - qnorm( 1 - PS ) )
  } else {
    -beta / ( z_half + qnorm( PS ) )
  }
  if ( seS <= 0 ) return( NA_real_ )
  sigma2_unit / seS^2
}

#' TOST power for clinical non-relevance (\eqn{\beta} inside [\code{Binf}, \code{Bsup}]).
#' @noRd
#' @keywords internal
.powerNonRelevance = function( beta, SE, Binf, Bsup, z_one ) {
  # When the CI width needed for TOST exceeds the equivalence window, power is 0.
  if ( 2 * z_one >= ( Bsup - Binf ) / SE ) return( 0 )
  pnorm( -z_one + ( Bsup - beta ) / SE ) -
    pnorm(  z_one + ( Binf - beta ) / SE )
}

#' Expand an upper search bound until \code{f(N) >= 0} (sample-size root finding).
#' @noRd
#' @keywords internal
.nBracketSampleSizeRoot = function( f, N_min = 2, N_max = 1e7, growth = 2 ) {
  if ( f( N_min ) >= 0 ) return( N_min )
  N_hi = N_min
  repeat {
    N_next = min( N_hi * growth, N_max )
    if ( N_next <= N_hi ) break
    if ( f( N_next ) >= 0 ) return( N_next )
    if ( N_next >= N_max ) break
    N_hi = N_next
  }
  N_max
}

#' Sample size for TOST non-relevance power (uniroot on power − target).
#' @noRd
#' @keywords internal
.nRequiredNonRelevance = function( beta, sigma2_unit, Binf, Bsup, z_one, PS ) {
  if ( beta <= Binf || beta >= Bsup ) return( NA_real_ )

  f = function( N ) {
    .powerNonRelevance( beta, sqrt( sigma2_unit / N ), Binf, Bsup, z_one ) - PS
  }
  p_max = .powerNonRelevance( beta, 1e-12, Binf, Bsup, z_one )

  if ( f( 2 ) >= 0 ) return( 2 )
  if ( p_max < PS ) return( NA_real_ )
  # Near-asymptotic power: report Inf rather than a huge finite N.
  if ( p_max - PS < 0.005 ) return( Inf )

  N_hi = .nBracketSampleSizeRoot( f )
  if ( f( N_hi ) < 0 ) return( NA_real_ )

  tryCatch(
    uniroot( f, interval = c( 2, N_hi ) )$root,
    error = function( e ) NA_real_
  )
}

#' Power that \eqn{\beta} lies outside the equivalence window (clinical relevance).
#' @noRd
#' @keywords internal
.powerRelevance = function( beta, SE, Binf, Bsup, z_one ) {
  pnorm( -z_one + ( Binf - beta ) / SE ) +
    1 - pnorm(  z_one + ( Bsup - beta ) / SE )
}

#' Sample size for clinical-relevance power (beta outside [\code{Binf}, \code{Bsup}]).
#' @noRd
#' @keywords internal
.nRequiredRelevance = function( beta, sigma2_unit, Binf, Bsup, z_one, PS ) {
  if ( beta >= Binf && beta <= Bsup ) return( NA_real_ )
  seR = if ( beta > Bsup ) {
    ( Bsup - beta ) / ( qnorm( 1 - PS ) - z_one )
  } else {
    ( Binf - beta ) / ( qnorm( PS ) + z_one )
  }
  if ( is.na( seR ) || seR <= 0 ) return( NA_real_ )
  sigma2_unit / seR^2
}

#' Wald significance rows: power at current SE, N to reach target power.
#'
#' \code{Power} uses the design SE; \code{N_Required} uses unit variance
#' \code{sigma2} (already scaled by FIM type in \code{.covariateTestSeData}).
#' @noRd
#' @keywords internal
.significanceTable = function( idx, names, values, SE, RSE, sigma2, z_half,
                               target_power, N0 ) {
  if ( !length( idx ) ) return( data.frame() )
  map( idx, function( i ) {
    b     = values[i]
    se_i  = SE[i]
    N_req = .nRequiredSignificance( b, sigma2[i], z_half, target_power )
    p_N0  = .powerSignificance( b, se_i, z_half )
    data.frame(
      Parameter  = names[i],
      Value      = round( b,     4 ),
      SE         = round( se_i,  4 ),
      RSE        = round( RSE[i], 2 ),
      Power      = round( p_N0 * 100, 1 ),
      N_Required = .formatNRequired( N_req ),
      stringsAsFactors = FALSE,
      row.names        = NULL
    )
  }) |> list_rbind()
}

#' Whether covariate effects use the additive link \eqn{\theta=\mu(1+\beta\cdot cov)}.
#' @noRd
#' @keywords internal
.covariateTestIsAdditive = function( pfimproject ) {
  eq = prop( pfimproject, "modelCovariatesEquation" )
  if ( length( eq ) == 0L || is.null( eq ) ) return( FALSE )
  if ( is.character( eq ) )
    return( identical( tolower( eq[[ 1L ]] ), "additive" ) )
  inherits( eq, "S7_object" ) && S7::S7_inherits( eq, Additive )
}

#' Extract SE/RSE rows and current sample size for covariate tests.
#'
#' For \strong{population} FIM, \code{sigma2 = SE^2 * N0} recovers the unit
#' variance because the FIM is scaled by arm size. For \strong{individual} /
#' \strong{Bayesian} FIM (not scaled by \eqn{N}), \code{N_scale = 1} so that
#' \code{N_Required} counts independent replications of the subject-level design.
#' @noRd
#' @keywords internal
.covariateTestSeData = function( pfimproject ) {
  fim  = setEvaluationFim( prop( pfimproject, "fim" ), pfimproject )
  seDF = prop( fim, "SEAndRSE" )$SEAndRSE
  names = rownames( seDF )
  N0 = sum( map_dbl( prop( prop( pfimproject, "designs" )[[ 1L ]], "arms" ),
                     ~ prop( .x, "size" ) ) )
  fimType = tolower( as.character( prop( pfimproject, "fimType" ) ) )
  is_pop  = identical( fimType, "population" ) ||
    S7::S7_inherits( fim, PopulationFim )
  list(
    names     = names,
    values    = set_names( seDF$parametersValues, names ),
    SE        = set_names( seDF$SE, names ),
    RSE       = set_names( seDF$RSE, names ),
    N0        = N0,
    N_scale   = if ( is_pop ) N0 else 1,
    fimType   = fimType,
    additive  = .covariateTestIsAdditive( pfimproject ),
    idx_beta  = which( str_starts( names, "\u03b2_" ) | str_starts( names, "beta_" ) )
  )
}

#' Map clinical log-ratio window to the beta / display scale for TOST tables.
#' @noRd
#' @keywords internal
.covariateTestEffectScale = function( values, SE, z_one, thetaL, thetaU, additive ) {
  if ( additive ) {
    # Additive link: relative factor is (1 + beta); keep power formulas on beta with
    # bounds mapped from the same clinical ratio window [e^{theta_L}, e^{theta_U}].
    list(
      effect  = 1 + values,
      lowerCI = 1 + values - z_one * SE,
      upperCI = 1 + values + z_one * SE,
      Binf    = exp( thetaL ) - 1,
      Bsup    = exp( thetaU ) - 1
    )
  } else {
    # Exponential link: display e^beta; TOST power stays on the log-ratio (beta) scale.
    list(
      effect  = exp( values ),
      lowerCI = exp( values - z_one * SE ),
      upperCI = exp( values + z_one * SE ),
      Binf    = thetaL,
      Bsup    = thetaU
    )
  }
}

#' Run significance / TOST non-relevance / relevance tables on \eqn{\beta} effects.
#'
#' Pipeline: pull SE/RSE -> classify betas inside/outside the equivalence window
#' on the relative-effect scale -> fill requested tables -> wrap in \code{CovariateTest}.
#' @name covariateTest
#' @keywords internal
method( covariateTest, Evaluation ) = function( pfimproject,
                                                thetaL       = log( 0.80 ),
                                                thetaU       = log( 1.25 ),
                                                target_power = 0.90,
                                                alpha        = 0.05,
                                                tests        = c(
                                                  "significance",
                                                  "nonRelevance", "relevance"
                                                ) ) {

  if ( thetaL >= thetaU ) {
    stop( "thetaL must be strictly less than thetaU on the log-ratio scale.", call. = FALSE )
  }

  tests_norm = .normalizeCovariateTests( tests )
  z_half     = qnorm( 1 - alpha / 2 )
  z_one      = qnorm( 1 - alpha       )
  d          = .covariateTestSeData( pfimproject )
  # Unit variance for N-scaling: population uses SE^2*N0; Ind/Bayes use SE^2.
  sigma2     = d$SE^2 * d$N_scale

  sc = .covariateTestEffectScale(
    d$values, d$SE, z_one, thetaL, thetaU, d$additive
  )
  effectBeta  = sc$effect[ d$idx_beta ]
  ratioLo     = exp( thetaL )
  ratioHi     = exp( thetaU )
  # Non-relevance: point relative effect inside the clinical ratio window.
  idx_beta_NR = d$idx_beta[ effectBeta >= ratioLo & effectBeta <= ratioHi ]
  idx_beta_R  = d$idx_beta[ effectBeta <  ratioLo | effectBeta >  ratioHi ]
  lowerCI     = sc$lowerCI
  upperCI     = sc$upperCI
  Binf        = sc$Binf
  Bsup        = sc$Bsup

  settings = list(
    thetaL       = thetaL,
    thetaU       = thetaU,
    target_power = target_power,
    alpha        = alpha,
    N0           = d$N0,
    N_scale      = d$N_scale,
    fimType      = d$fimType,
    additive     = d$additive,
    tests        = tests_norm
  )

  df_cov = if ( "significance" %in% tests_norm ) {
    .significanceTable(
      d$idx_beta, d$names, d$values, d$SE, d$RSE, sigma2,
      z_half, target_power, d$N0
    )
  } else {
    data.frame()
  }

  df_NR = if ( "nonRelevance" %in% tests_norm ) {
    map( idx_beta_NR, function( i ) {
      b     = d$values[i]
      se_i  = d$SE[i]
      ratio = sc$effect[i]
      N_req = .nRequiredNonRelevance(
        b, sigma2[i], Binf, Bsup, z_one, target_power
      )
      p_N0  = .powerNonRelevance( b, se_i, Binf, Bsup, z_one )
      data.frame(
        Parameter  = d$names[i],
        Value      = round( b,          4 ),
        SE         = round( se_i,       4 ),
        RSE        = round( d$RSE[i],   2 ),
        Ratio      = round( ratio,      4 ),
        IC_Inf     = round( lowerCI[i], 4 ),
        IC_Sup     = round( upperCI[i], 4 ),
        Power      = round( p_N0 * 100, 1 ),
        N_Required = .formatNRequired( N_req ),
        stringsAsFactors = FALSE,
        row.names        = NULL
      )
    }) |> list_rbind()
  } else {
    data.frame()
  }

  df_R = if ( "relevance" %in% tests_norm ) {
    map( idx_beta_R, function( i ) {
      b     = d$values[i]
      se_i  = d$SE[i]
      ratio = sc$effect[i]
      N_req = .nRequiredRelevance( b, sigma2[i], Binf, Bsup, z_one, target_power )
      p_N0  = .powerRelevance( b, se_i, Binf, Bsup, z_one )
      data.frame(
        Parameter  = d$names[i],
        Value      = round( b,          4 ),
        SE         = round( se_i,       4 ),
        RSE        = round( d$RSE[i],   2 ),
        Ratio      = round( ratio,      4 ),
        IC_Inf     = round( lowerCI[i], 4 ),
        IC_Sup     = round( upperCI[i], 4 ),
        Power      = round( p_N0 * 100, 1 ),
        N_Required = .formatNRequired( N_req ),
        stringsAsFactors = FALSE,
        row.names        = NULL
      )
    }) |> list_rbind()
  } else {
    data.frame()
  }

  CovariateTest(
    covariate    = df_cov,
    parameter    = data.frame(),
    nonRelevance = df_NR,
    relevance    = df_R,
    settings     = settings
  )
}

#' Build HTML tables from a covariate test result
#'
#' @param covariateTestResult A \code{CovariateTest} object.
#' @param tests Optional character vector overriding \code{settings$tests}.
#' @return Named list of \code{knitr_kable} objects (or \code{NULL} if empty).
#'   Element \code{significance} is the \code{covariate} slot (report alias).
#' @export
getCovariateTestTables = function( covariateTestResult, tests = NULL ) {

  tests = .covariateTestResolvedTests( covariateTestResult, tests )
  settings = prop( covariateTestResult, "settings" )

  .makeKable = function( df, kind ) {
    if ( is.null( df ) || !nrow( df ) ) return( NULL )
    caption = paste(
      .covariateTestSlotTitle( kind ),
      .covariateTestHeader( kind, settings ),
      sep = " - "
    )
    kbl( df, align = "c", caption = caption ) |>
      kable_styling(
        bootstrap_options = c( "hover", "striped", "bordered" ),
        full_width        = FALSE,
        position          = "center",
        font_size         = 12
      )
  }

  out = list(
    significance = NULL,
    parameter    = NULL,
    nonRelevance = NULL,
    relevance    = NULL
  )

  slots = intersect( c( "significance", "nonRelevance", "relevance" ), tests )
  kables = compact( set_names(
    map( slots, ~ .makeKable( .covariateTestSlotDf( covariateTestResult, .x ), .x ) ),
    slots
  ) )
  out[ names( kables ) ] = kables
  out
}

#' @keywords internal
method( show, CovariateTest ) = function( object ) {
  walk( .covariateTestResolvedTests( object ), ~ .printCovariateTestSlot( object, .x ) )
  invisible( object )
}

#' Save covariate test results to a text file
#'
#' @param relevanceResults A \code{CovariateTest} object.
#' @param file_name Output file name (or full path when \code{folder} is \code{NULL}).
#' @param folder Optional directory; created if missing.
#' @param tests Optional character vector overriding \code{settings$tests}.
#' @return Invisibly, the output file path.
#' @export
saveCovariateTest = function( relevanceResults, file_name, folder = NULL,
                              tests = NULL ) {

  tests = .covariateTestResolvedTests( relevanceResults, tests )
  body  = unlist( map( tests, ~ .formatCovariateTestSlotLines( relevanceResults, .x ) ) )

  output_lines = c(
    strrep( "=", 60 ),
    "   PFIM -- Covariate test",
    strrep( "=", 60 ),
    "",
    body,
    strrep( "=", 60 ),
    "   END OF FILE",
    strrep( "=", 60 )
  )

  if ( !is.null( folder ) ) {
    if ( !dir.exists( folder ) ) {
      warning( sprintf( "saveCovariateTest: directory '%s' did not exist and was created.", folder ) )
      dir.create( folder, recursive = TRUE )
    }
    file_name = file.path( folder, file_name )
  }

  writeLines( output_lines, file_name )
  message( "Saved: ", file_name )
  invisible( file_name )
}

#' @rdname covariateTest
#' @export
tost = covariateTest
