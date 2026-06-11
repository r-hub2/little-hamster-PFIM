#' CovariateTest Class
#'
#' @description
#' \code{CovariateTest} is the S7 class storing the covariate test results.
#'
#' @details
#' Three slots are available:
#' \describe{
#'   \item{significance}{Significance test (two-sided Wald, H0: beta = 0)}
#'   \item{nonRelevance}{Non-relevance / TOST test (H1: ratio in [0.80, 1.25])}
#'   \item{relevance}{Relevance test (H1: ratio outside [0.80, 1.25])}
#' }
#'
#' @param significance A \code{data.frame} for the significance test results.
#' @param nonRelevance A \code{data.frame} for the non-relevance test results.
#' @param relevance    A \code{data.frame} for the relevance test results.
#'
#' @name CovariateTest-class
#' @rdname CovariateTest-class
#' @aliases CovariateTest
#' @usage NULL
#' @export
CovariateTest = new_class( "CovariateTest", package = "PFIM",
                           properties = list(
                             significance = new_property( class_any, default = data.frame() ),
                             nonRelevance = new_property( class_any, default = data.frame() ),
                             relevance    = new_property( class_any, default = data.frame() )
                           )
)
S4_register( CovariateTest )

#' covariateTest
#'
#' @description
#' Assess statistical significance and clinical relevance of covariate effects from the FIM.
#'
#' @details
#' Three tests are computed:
#' \enumerate{
#'   \item \strong{Significance} (slot \code{significance}) -- covariate effects \eqn{\beta}
#'     (and population typical values \eqn{\mu} when \code{fimType = "population"}).
#'     Two-sided Wald test (H0: effect = 0 on the log scale for \eqn{\beta}).
#'     \eqn{PS = 1 - \Phi(z_{1-\alpha/2} - \beta/SE) + \Phi(-z_{1-\alpha/2} - \beta/SE)}.
#'   \item \strong{Non-relevance / TOST} (slot \code{nonRelevance}) -- \eqn{\beta} with ratio in [0.80, 1.25].
#'     \eqn{PNR = \Phi(-z_{1-\alpha} + (B_{sup}-\beta)/SE) - \Phi(z_{1-\alpha} + (B_{inf}-\beta)/SE)}.
#'   \item \strong{Relevance} (slot \code{relevance}) -- \eqn{\beta} with ratio outside [0.80, 1.25].
#'     \eqn{PR = \Phi(-z_{1-\alpha} + (B_{inf}-\beta)/SE) + 1 - \Phi(z_{1-\alpha} + (B_{sup}-\beta)/SE)}.
#' }
#'
#' @section Standard errors and interpretation:
#' Standard errors (SE) and power calculations depend on \code{fimType}:
#' \describe{
#'   \item{\code{population}}{
#'     FIM includes \eqn{\mu}, \eqn{\beta}, IIV (\eqn{\omega^2}), IOV (\eqn{\gamma^2}), and residual
#'     \eqn{\sigma}. SEs are appropriate for population design and sample-size formulas on all blocks.
#'   }
#'   \item{\code{individual}}{
#'     Subject-level FIM: block-diagonal \eqn{(\mu,\beta)} and \eqn{\sigma} without \eqn{\omega^2}/\eqn{\gamma^2}
#'     in the matrix. SEs for \eqn{\mu} are often large (no shrinkage); \strong{clinical tests (slots 2--3)
#'     should focus on \eqn{\beta}}. Slot 1 includes \eqn{\mu} only for \code{population}.
#'   }
#'   \item{\code{Bayesian}}{
#'     Prior-shrunk FIM on typical PK parameters (\eqn{\mu} block) plus a separate \eqn{\beta} block when
#'     covariates are present; residual \eqn{\sigma} is not in this matrix. SEs reflect design information
#'     with population prior on \eqn{\mu}; use slots 2--3 for covariate ratios.
#'   }
#' }
#'
#' @param pfimproject    An \code{Evaluation} object (after \code{run()}).
#' @param thetaL         Lower bound on the log-ratio scale (default: log(0.80)).
#' @param thetaU         Upper bound on the log-ratio scale (default: log(1.25)).
#' @param target_power   Target power (default: 0.90).
#' @param alpha          Nominal type-I error rate (default: 0.05).
#'
#' @return A \code{CovariateTest} object containing the three test dataframes.
#'
#' @name covariateTest
#' @rdname covariateTest
#' @aliases covariateTest
#' @include Evaluation.R
#' @export
covariateTest = new_generic( "covariateTest", "pfimproject",
                             fun = function( pfimproject,
                                             thetaL       = log( 0.80 ),
                                             thetaU       = log( 1.25 ),
                                             target_power = 0.90,
                                             alpha        = 0.05 ) {
                               S7_dispatch()
                             }
)

# Private pure functions -- package level

# -- Status string helpers (%.0f: R %d rejects non-integer doubles / Inf) ------
.statusSignificanceN = function( N_req, N0, p_N0 ) {
  if ( is.infinite( N_req ) ) {
    return( sprintf(
      "Asymptotic (N -> Inf) -- power at N = %.0f: %.1f%%", N0, p_N0 * 100
    ) )
  }
  if ( N_req <= 2 ) {
    return( sprintf( "OK (N = 2) -- power at N = %.0f: %.1f%%", N0, p_N0 * 100 ) )
  }
  sprintf(
    "N = %.0f required -- power at N = %.0f: %.1f%%",
    ceiling( N_req ), N0, p_N0 * 100
  )
}

.statusSampleSizeN = function( N_req, N0, p_N0 ) {
  if ( is.infinite( N_req ) ) {
    return( sprintf(
      "Asymptotic (N -> Inf) -- power at N = %.0f: %.1f%%", N0, p_N0 * 100
    ) )
  }
  if ( N_req <= 2 ) {
    return( sprintf( "OK (N = 2) -- power at N = %.0f: %.1f%%", N0, p_N0 * 100 ) )
  }
  sprintf(
    "N = %.0f required -- power at N = %.0f: %.1f%%",
    ceiling( N_req ), N0, p_N0 * 100
  )
}

# -- Significance --------------------------------------------------------------
.powerSignificance = function( beta, SE, z_half ) {
  1 - pnorm( z_half - beta / SE ) + pnorm( -z_half - beta / SE )
}

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

# -- Non-relevance / TOST ------------------------------------------------------
.powerNonRelevance = function( beta, SE, Binf, Bsup, z_one ) {
  if ( 2 * z_one >= ( Bsup - Binf ) / SE ) return( 0 )
  pnorm( -z_one + ( Bsup - beta ) / SE ) -
    pnorm(  z_one + ( Binf - beta ) / SE )
}

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

.nRequiredNonRelevance = function( beta, sigma2_unit, Binf, Bsup, z_one, PS ) {
  if ( beta <= Binf || beta >= Bsup ) return( NA_real_ )

  f = function( N ) {
    .powerNonRelevance( beta, sqrt( sigma2_unit / N ), Binf, Bsup, z_one ) - PS
  }
  p_max = .powerNonRelevance( beta, 1e-12, Binf, Bsup, z_one )

  if ( f( 2 ) >= 0 ) return( 2 )
  if ( p_max < PS ) return( NA_real_ )
  if ( p_max - PS < 0.005 ) return( Inf )

  N_hi = .nBracketSampleSizeRoot( f )
  if ( f( N_hi ) < 0 ) return( NA_real_ )

  tryCatch(
    uniroot( f, interval = c( 2, N_hi ) )$root,
    error = function( e ) NA_real_
  )
}

# -- Relevance -----------------------------------------------------------------
.powerRelevance = function( beta, SE, Binf, Bsup, z_one ) {
  pnorm( -z_one + ( Binf - beta ) / SE ) +
    1 - pnorm(  z_one + ( Bsup - beta ) / SE )
}

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

# Method dispatch
method( covariateTest, Evaluation ) = function( pfimproject,
                                                thetaL       = log( 0.80 ),
                                                thetaU       = log( 1.25 ),
                                                target_power = 0.90,
                                                alpha        = 0.05 ) {

  z_half = qnorm( 1 - alpha / 2 )
  z_one  = qnorm( 1 - alpha       )

  fim  = prop( pfimproject, "fim" )
  fim  = setEvaluationFim( fim, pfimproject )
  seDF = prop( fim, "SEAndRSE" )$SEAndRSE

  allParametersNames  = rownames( seDF )
  betaHat             = set_names( seDF$parametersValues, allParametersNames )
  standardErrorsValue = set_names( seDF$SE,               allParametersNames )
  RSE                 = set_names( seDF$RSE,               allParametersNames )

  N0 = sum( map_dbl( prop( prop( pfimproject, "designs" )[[1L]], "arms" ),
                     ~ prop( .x, "size" ) ) )

  idx_mu   = which( str_starts( allParametersNames, "\u03bc_" ) )
  idx_beta = which( str_starts( allParametersNames, "\u03b2_" ) |
                      str_starts( allParametersNames, "beta_"  ) )

  fim_type   = tolower( prop( pfimproject, "fimType" ) )
  idx_signif = if ( fim_type %in% c( "individual", "bayesian" ) ) {
    idx_beta
  } else {
    c( idx_mu, idx_beta )
  }

  ratioBeta   = exp( betaHat[ idx_beta ] )
  idx_beta_NR = idx_beta[ ratioBeta >= 0.80 & ratioBeta <= 1.25 ]
  idx_beta_R  = idx_beta[ ratioBeta <  0.80 | ratioBeta >  1.25 ]

  lowerCI      = exp( betaHat - z_one * standardErrorsValue )
  upperCI      = exp( betaHat + z_one * standardErrorsValue )
  sigma2_units = standardErrorsValue^2 * N0

  # -- Slot 1: Significance (beta; mu only for population FIM) -------------------
  df_signif = map( idx_signif, function( i ) {
    b     = betaHat[i]
    SE    = standardErrorsValue[i]
    N_req = .nRequiredSignificance( b, sigma2_units[i], z_half, target_power )
    p_N0  = .powerSignificance( b, SE, z_half )

    status = if ( is.na( N_req ) ) {
      "Impossible (beta = 0)"
    } else {
      .statusSignificanceN( N_req, N0, p_N0 )
    }

    data.frame(
      Parameter  = allParametersNames[i],
      Value      = round( b,      4 ),
      SE         = round( SE,     4 ),
      RSE        = round( RSE[i], 2 ),
      Power_N0   = round( p_N0 * 100, 1 ),
      N_Required = ceiling( N_req ),
      Status     = status,
      stringsAsFactors = FALSE,
      row.names        = NULL
    )
  }) |> list_rbind()

  # -- Slot 2: Non-relevance (beta with ratio in [0.80, 1.25]) ------------------
  df_NR = map( idx_beta_NR, function( i ) {
    b     = betaHat[i]
    SE    = standardErrorsValue[i]
    ratio = exp( b )
    p_max = .powerNonRelevance( b, 1e-12, thetaL, thetaU, z_one )
    N_req = .nRequiredNonRelevance( b, sigma2_units[i], thetaL, thetaU, z_one, target_power )
    p_N0  = .powerNonRelevance( b, SE, thetaL, thetaU, z_one )

    status = if ( is.na( N_req ) ) {
      sprintf( "Impossible -- max. power (N->Inf): %.1f%%", p_max * 100 )
    } else if ( is.infinite( N_req ) ) {
      sprintf( "Asymptotic -- max. power (N->Inf): %.1f%% (target: %.0f%%)",
               p_max * 100, target_power * 100 )
    } else {
      .statusSampleSizeN( N_req, N0, p_N0 )
    }

    data.frame(
      Parameter  = allParametersNames[i],
      Value      = round( b,          4 ),
      SE         = round( SE,         4 ),
      RSE        = round( RSE[i],     2 ),
      Ratio      = round( ratio,      4 ),
      IC90_Inf   = round( lowerCI[i], 4 ),
      IC90_Sup   = round( upperCI[i], 4 ),
      Power_N0   = round( p_N0  * 100, 1 ),
      Power_max  = round( p_max * 100, 1 ),
      N_Required = ceiling( N_req ),
      Status     = status,
      stringsAsFactors = FALSE,
      row.names        = NULL
    )
  }) |> list_rbind()

  # -- Slot 3: Relevance (beta with ratio outside [0.80, 1.25]) -----------------
  df_R = map( idx_beta_R, function( i ) {
    b     = betaHat[i]
    SE    = standardErrorsValue[i]
    ratio = exp( b )
    N_req = .nRequiredRelevance( b, sigma2_units[i], thetaL, thetaU, z_one, target_power )
    p_N0  = .powerRelevance( b, SE, thetaL, thetaU, z_one )

    status = if ( is.na( N_req ) ) {
      "Impossible -- max. power insufficient"
    } else {
      .statusSampleSizeN( N_req, N0, p_N0 )
    }

    data.frame(
      Parameter  = allParametersNames[i],
      Value      = round( b,          4 ),
      SE         = round( SE,         4 ),
      RSE        = round( RSE[i],     2 ),
      Ratio      = round( ratio,      4 ),
      IC90_Inf   = round( lowerCI[i], 4 ),
      IC90_Sup   = round( upperCI[i], 4 ),
      Power_N0   = round( p_N0 * 100, 1 ),
      N_Required = ceiling( N_req ),
      Status     = status,
      stringsAsFactors = FALSE,
      row.names        = NULL
    )
  }) |> list_rbind()

  CovariateTest(
    significance = df_signif,
    nonRelevance = df_NR,
    relevance    = df_R
  )
}

#' getCovariateTestTables: build kableExtra tables from a CovariateTest object.
#'
#' Returns a named list of three kable objects (NULL when a slot is empty).
#' Intended for use inside \code{Report} methods.
#'
#' @param covariateTestResult A \code{CovariateTest} object.
#' @return Named list: \code{significance}, \code{nonRelevance}, \code{relevance}.
#' @export
getCovariateTestTables = function( covariateTestResult ) {

  .makeKable = function( df, caption ) {
    if ( is.null( df ) || nrow( df ) == 0L ) return( NULL )
    kbl( df, align = "c", caption = caption ) |>
      kable_styling(
        bootstrap_options = c( "hover", "striped", "bordered" ),
        full_width        = FALSE,
        position          = "center",
        font_size         = 12
      )
  }

  list(
    significance = .makeKable(
      prop( covariateTestResult, "significance" ),
      paste0( "Statistical Significance \u2014 Wald test" )
    ),
    nonRelevance = .makeKable(
      prop( covariateTestResult, "nonRelevance" ),
      paste0( "Clinical Non-Relevance \u2014 TOST" )
    ),
    relevance = .makeKable(
      prop( covariateTestResult, "relevance" ),
      paste0( "Clinical Relevance" )
    )
  )
}

method( show, CovariateTest ) = function( object ) {

  .print_slot = function( df, title ) {
    cat( "\n", strrep( "=", 60 ), "\n", sep = "" )
    cat( " ", title, "\n" )
    cat( strrep( "=", 60 ), "\n\n", sep = "" )
    if ( nrow( df ) == 0L ) {
      cat( "  (no parameter in this slot)\n" )
    } else {
      print( df, row.names = FALSE )
    }
  }

  .print_slot( prop( object, "significance" ), "SLOT 1 -- Statistical Significance"      )
  .print_slot( prop( object, "nonRelevance" ), "SLOT 2 -- Clinical Non-Relevance (TOST)" )
  .print_slot( prop( object, "relevance"    ), "SLOT 3 -- Clinical Relevance"             )

  invisible( object )
}

#' saveCovariateTest: save the three slots to a text file.
#' @name saveCovariateTest
#' @param relevanceResults A \code{CovariateTest} object.
#' @param file_name        Output file name.
#' @param folder           Destination folder (optional).
#' @export
saveCovariateTest = function( relevanceResults, file_name, folder = NULL ) {

  .format_slot = function( df, title ) {
    c(
      strrep( "=", 60 ),
      paste0( "  ", title ),
      strrep( "=", 60 ),
      "",
      if ( nrow( df ) == 0L ) {
        "  (no parameter in this slot)"
      } else {
        utils::capture.output( print( df, row.names = FALSE ) )
      },
      ""
    )
  }

  output_lines = c(
    strrep( "=", 60 ),
    "   PFIM -- Clinical Relevance",
    strrep( "=", 60 ),
    "",
    .format_slot( prop( relevanceResults, "significance" ), "SLOT 1 -- Statistical Significance"      ),
    .format_slot( prop( relevanceResults, "nonRelevance" ), "SLOT 2 -- Clinical Non-Relevance (TOST)" ),
    .format_slot( prop( relevanceResults, "relevance"    ), "SLOT 3 -- Clinical Relevance"             ),
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
}

#' tost
#'
#' Backward-compatible wrapper for \code{covariateTest}.
#'
#' @param pfimproject    An \code{Evaluation} object.
#' @param thetaL         Lower bound on log-ratio scale (default: log(0.80)).
#' @param thetaU         Upper bound on log-ratio scale (default: log(1.25)).
#' @param target_power   Target power (default: 0.90).
#' @param alpha          Nominal type-I error rate (default: 0.05).
#' @return A \code{CovariateTest} object.
#' @export
tost = covariateTest
