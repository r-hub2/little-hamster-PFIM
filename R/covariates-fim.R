# Covariate x occasion FIM assembly.
#
# Population mixes Fisher blocks: M = Sum_c pi_c M_c.
# Individual / Bayesian mix covariances: C̄ = Sum_c pi_c M_c^{-1}, M_eff = C̄^{-1}.
# Occasions inside a combination belong to one subject (stacked / independent IOV).
# Plot helpers average gradients and variances; they must not feed Fisher blocks.

#' Detect nested covariatexoccasion arm evaluations (gradients or variances).
#' @param x \code{evaluationGradients} or \code{evaluationVariance} from an arm.
#' @return Logical scalar.
#' @noRd
#' @keywords internal
.isNestedArmEvaluation = function( x ) {
  length( x ) > 0L && is.list( x[[1L]] ) &&
    ( !is.null( x[[1L]]$gradients ) || !is.null( x[[1L]]$variances ) )
}

#' Detect nested \code{evaluateModel()} output (covariate x occasion combinations).
#' @param x Arm \code{evaluationModel} structure.
#' @return Logical scalar.
#' @noRd
#' @keywords internal
.isNestedEvaluationModel = function( x ) {
  length( x ) > 0L && is.list( x[[1L]] ) && !is.null( x[[1L]]$evaluations )
}

#' Proportion-weighted mean model curves for report plots.
#'
#' Averages over occasions within each combination, then weights by combination
#' proportion so plotted curves match the mixture used in the population FIM.
#' Occasions inside a combo get equal weight (\code{pi_c / n_occ}); this is for
#' display only - pop FIM still sums \eqn{\pi_c M_c} at the Fisher-block level.
#' @param nested_eval Nested evaluation structure from covariatexoccasion path.
#' @return Named list of \code{data.frame}s (\code{time}, output) per outcome.
#' @noRd
#' @keywords internal
.aggregateEvaluationModelForPlot = function( nested_eval ) {
  outputNames = names( pluck( nested_eval, 1L, "evaluations", 1L, "evaluation" ) )
  map( outputNames, function( out_name ) {
    # Equal weight within combo, then pi_c across combinations.
    weighted = nested_eval |>
      map( function( combo ) {
        w_occ = combo$proportion / length( combo$evaluations )
        map( combo$evaluations, function( occ ) {
          df = occ$evaluation[[ out_name ]]
          list( time = df$time, y = w_occ * df[[ out_name ]] )
        } )
      } ) |>
      list_flatten()
    acc = reduce(
      weighted[-1L],
      function( acc, piece ) {
        idx = match( acc$time, piece$time )
        acc$y = acc$y + piece$y[ idx ]
        acc
      },
      .init = weighted[[ 1L ]]
    )
    out = data.frame( time = acc$time )
    out[[ out_name ]] = acc$y
    out
  } ) |> set_names( outputNames )
}

#' Names of declared covariate effects, including zeros (null-hypothesis \eqn{\beta}).
#' Overlay vectors from \code{getCovariateEffects()} pad every model parameter
#' with 0, so callers must pass the user \code{effects} list, not the overlay.
#' @noRd
#' @keywords internal
.namedEffectParams = function( effectVec ) {
  if ( is.null( effectVec ) || !length( effectVec ) )
    return( character( 0L ) )
  nms = names( effectVec )
  if ( is.null( nms ) )
    return( character( 0L ) )
  nms[ !is.na( nms ) & nzchar( nms ) ]
}

#' Build a dictionary of declared covariate beta effect sizes.
#'
#' Keys follow \code{beta_<param>_<covariate>_<category>} (reference category skipped).
#' Declared zeros are kept so the SE of \eqn{\beta} under the null remains defined.
#' @param modelCovariates List of covariate objects.
#' @return Named numeric vector of beta values.
#' @noRd
#' @keywords internal
.betaDictFromCovariates = function( modelCovariates ) {
  if ( length( modelCovariates ) == 0L )
    return( set_names( numeric( 0L ), character( 0L ) ) )
  lst = modelCovariates |>
    map( function( cov ) {
      covName    = prop( cov, "name" )
      categories = prop( cov, "categories" )
      effects    = prop( cov, "effects" )
      ( seq_len( length( categories ) - 1L ) + 1L ) |>
        map( function( icat ) {
          cat = categories[[ icat ]]
          if ( !cat %in% names( effects ) ) return( list() )
          eff = effects[[ cat ]]
          .namedEffectParams( eff ) |>
            map( ~ list(
              key = paste0( "beta_", .x, "_", covName, "_", cat ),
              val = as.numeric( eff[[ .x ]] )
            ) )
        }) |> list_flatten()
    }) |> list_flatten()
  if ( length( lst ) == 0L )
    return( set_names( numeric( 0L ), character( 0L ) ) )
  set_names( map_dbl( lst, "val" ), map_chr( lst, "key" ) )
}

#' Map internal beta names to numeric values from covariate definitions.
#'
#' @param modelCovariates List of covariate objects.
#' @param betaInternal Character vector of \code{beta_*} names
#' @return Numeric vector aligned with \code{betaInternal}.
#' @noRd
#' @keywords internal
.betaValuesFromCovariates = function( modelCovariates, betaInternal ) {
  if ( length( betaInternal ) == 0L ) return( numeric( 0L ) )
  dict = .betaDictFromCovariates( modelCovariates )
  vals = map_dbl( betaInternal, ~ dict[[ .x ]] )
  vals
}

#' Collapse nested combination-by-occasion gradients to one matrix (all outputs stacked).
#'
#' Plot / label helper only: \eqn{E[G] = \sum_c \pi_c \mathrm{mean}_k G_{c,k}}.
#' Individual / Bayesian FIM must \emph{not} use this for Fisher blocks - use
#' \code{.evaluateIndBayesMixtureFim()} (\eqn{\sum_c \pi_c \sum_k G^\top R^{-1} G}).
#'
#' @param model A \code{Model} object.
#' @param arm An \code{Arm} object with nested \code{evaluationGradients}.
#' @return Numeric matrix (stacked outputs x parameters).
#' @keywords internal
aggregateGradientsWithCovariates = function( model, arm ) {

  allGradientsData = prop( arm, "evaluationGradients" )
  outputNames      = prop( model, "outputNames" )

  if ( length( outputNames ) == 0L && length( allGradientsData ) > 0L ) {
    outputNames = names( pluck( allGradientsData, 1L, "gradients", 1L, "gradient" ) )
  }

  gradDfs = map( outputNames, function( outName ) {
    reduce(
      allGradientsData,
      function( acc, combinationData ) {
        grads = map( combinationData$gradients, ~ .x$gradient[[ outName ]] )
        meanGrad = reduce( grads, `+` ) / length( grads )
        piece = combinationData$proportion * meanGrad
        if ( is.null( acc ) ) piece else acc + piece
      },
      .init = NULL
    )
  })

  mat = list_rbind( gradDfs ) |> as.matrix()
  if ( !is.null( gradDfs[[1L]] ) && ncol( gradDfs[[1L]] ) > 0L )
    colnames( mat ) = colnames( gradDfs[[1L]] )
  mat
}

#' Expectation of residual variance over combinations / occasions (plots / flat V).
#'
#' Same averaging contract as \code{aggregateGradientsWithCovariates}. FIM assembly
#' for Individual / Bayesian uses \code{.evaluateIndBayesMixtureFim()} instead.
#'
#' @param arm An \code{Arm} object with nested \code{evaluationVariance}.
#' @return List with \code{errorVariance} and \code{sigmaDerivatives}.
#' @keywords internal
aggregateVarianceWithCovariates = function( arm ) {

  varianceResults = prop( arm, "evaluationVariance" )

  meanOccasionVariance = function( occVars ) {
    n = length( occVars )
    err = reduce(
      map( occVars, ~ as.matrix( .x$variance$errorVariance ) ),
      `+`
    ) / n
    sig = pluck( occVars, 1L, "variance", "sigmaDerivatives" )
    if ( n > 1L && length( sig ) ) {
      sig = map( seq_along( sig ), function( j ) {
        reduce(
          map( occVars, ~ pluck( .x, "variance", "sigmaDerivatives", j ) ),
          `+`
        ) / n
      } )
    }
    list( errorVariance = err, sigmaDerivatives = sig )
  }

  reduce(
    varianceResults,
    function( out, combo ) {
      p = combo$proportion
      v = meanOccasionVariance( combo$variances )
      if ( is.null( out ) ) {
        list(
          errorVariance    = p * v$errorVariance,
          sigmaDerivatives = map( v$sigmaDerivatives, ~ p * .x )
        )
      } else {
        list(
          errorVariance    = out$errorVariance + p * v$errorVariance,
          sigmaDerivatives = map2( out$sigmaDerivatives, v$sigmaDerivatives, `+` )
        )
      }
    },
    .init = NULL
  )
}

#' Stack per-output gradient frames for one (combo, occasion) leaf.
#' @noRd
#' @keywords internal
.indBayesOccasionGradientMat = function( comboGrad, occ, outputNames, feCols ) {
  dfs = map( outputNames, ~ comboGrad$gradients[[ occ ]]$gradient[[ .x ]] )
  mat = as.matrix( list_rbind( dfs ) )
  mat[ , feCols, drop = FALSE ]
}

#' Gradient columns needed to inflate residual V by IOV (every gamma>0).
#' @noRd
#' @keywords internal
.indBayesIovGradientCols = function( iovIdx, sampleGrad, bareNames ) {
  if ( !length( iovIdx ) || !length( sampleGrad ) )
    return( character( 0L ) )
  hits = map_chr( iovIdx, function( ii ) {
    muCol = paste0( "mu_", bareNames[[ ii ]] )
    if ( muCol %in% sampleGrad )
      return( muCol )
    if ( bareNames[[ ii ]] %in% sampleGrad )
      return( bareNames[[ ii ]] )
    NA_character_
  } )
  unique( hits[ !is.na( hits ) ] )
}

#' FO IOV inflation: \eqn{V = R + \sum_j \gamma_j^2 (\mu_j g_j)(\mu_j g_j)^\top}.
#'
#' One BLAS \code{tcrossprod} after column scaling; same algebra as a loop of
#' rank-1 updates.
#' @noRd
#' @keywords internal
.indBayesIovInflatedV = function( R, G_all, gradCols, bareNames, gamma_all, muChain ) {
  bare = sub( "^mu_", "", gradCols )
  ii   = match( bare, bareNames )
  keep = !is.na( ii ) & gamma_all[ ii ] > 0
  if ( !any( keep ) )
    return( R )
  scale = muChain[ ii[ keep ] ] * gamma_all[ ii[ keep ] ]
  R + tcrossprod( sweep( G_all[ , keep, drop = FALSE ], 2L, scale, `*` ) )
}

#' Reorder occasion gradient columns to the requested FE layout.
#' @noRd
#' @keywords internal
.indBayesAlignFeGradient = function( G_all, gradCols, feCols ) {
  fe_in = which( gradCols %in% feCols )
  G = G_all[ , fe_in, drop = FALSE ]
  G[ , match( feCols, gradCols[ fe_in ] ), drop = FALSE ]
}

#' Individual / Bayesian Fisher pieces for every occasion of one combination.
#' @noRd
#' @keywords internal
.indBayesOccasionBlocks = function( comboGrad, comboVar, nOcc, outputNames,
                                    gradCols, feCols, has_iov, bareNames,
                                    gamma_all, muChain, bayesian ) {
  map( seq_len( nOcc ), function( k ) {
    G_all = .indBayesOccasionGradientMat( comboGrad, k, outputNames, gradCols )
    R = as.matrix( pluck( comboVar, "variances", k, "variance", "errorVariance" ) )
    V = if ( has_iov )
      .indBayesIovInflatedV( R, G_all, gradCols, bareNames, gamma_all, muChain )
    else
      R
    V_inv = .safeCholInv( V )
    G = .indBayesAlignFeGradient( G_all, gradCols, feCols )
    Mb = crossprod( G, V_inv ) %*% G
    Mv = if ( !bayesian ) {
      dV = pluck( comboVar, "variances", k, "variance", "sigmaDerivatives" )
      if ( length( dV ) ) .computeMFVar( V_inv, dV ) else NULL
    } else {
      NULL
    }
    list( Mb = Mb, Mv = Mv )
  } )
}

#' Stack occasion Fisher blocks into one subject FIM for a covariate combination.
#' @noRd
#' @keywords internal
.indBayesComboSubjectFim = function( blocks, model, bayesian ) {
  Mc_b = reduce( map( blocks, "Mb" ), `+` )
  if ( bayesian )
    return( .pfimBayesianSubjectFim( Mc_b, model ) )
  mvs = compact( map( blocks, "Mv" ) )
  Mc_v = if ( length( mvs ) )
    reduce( mvs, `+` )
  else
    matrix( numeric( 0L ), 0L, 0L )
  as.matrix( Matrix::bdiag( Mc_b, Mc_v ) )
}

#' Individual / Bayesian subject FIM: harmonic mean over covariate strata.
#'
#' Per covariate combination \eqn{c}, occasions are stacked (sum of Fisher
#' blocks for one subject). Across combinations, precision is averaged:
#' \eqn{\bar C = \sum_c \pi_c M_c^{-1}}, \eqn{M_{\mathrm{eff}} = \bar C^{-1}}.
#' Covariate \code{beta} columns are dropped (MAP conditions on \eqn{\beta}).
#'
#' When \code{bayesian = TRUE}, each stratum uses
#' \code{.pfimBayesianSubjectFim()} (data FIM + \eqn{\Omega^{-1}}).
#' With IOV (\eqn{\gamma > 0}), both individual and Bayesian paths use
#' \eqn{V_k = R_k + F_k \mathrm{diag}(\gamma^2) F_k^\top}. The sum runs over
#' every parameter with \eqn{\gamma > 0} (even without IIV); \eqn{M_{\mathrm{data}}}
#' keeps only the requested \code{feCols} (eta / mu columns).
#'
#' @param model \code{Model} with nested arm evaluations.
#' @param arm Evaluated \code{Arm}.
#' @param feCols Character vector of estimable FE column names (beta stripped).
#' @param bayesian If \code{TRUE}, build Bayesian subject FIMs.
#' @return List with \code{fisherMatrix} (and \code{covariance}); legacy
#'   \code{MFbeta}/\code{MFVar} when \code{bayesian = FALSE} for callers that
#'   still split FE / residual blocks after a flat path.
#' @noRd
#' @keywords internal
.evaluateIndBayesMixtureFim = function( model, arm, feCols, bayesian = FALSE ) {

  feCols = .pfimSubjectFeCols( feCols )
  allGradientsData = prop( arm, "evaluationGradients" )
  varianceResults  = prop( arm, "evaluationVariance" )
  outputNames      = prop( model, "outputNames" )
  if ( length( outputNames ) == 0L && length( allGradientsData ) > 0L )
    outputNames = names( pluck( allGradientsData, 1L, "gradients", 1L, "gradient" ) )

  parameters = prop( model, "modelParameters" )
  muChain    = .pfimPopMuChainFactors( parameters )
  bareNames  = map_chr( parameters, ~ prop( .x, "name" ) )
  gamma_all  = vapply(
    parameters,
    function( x ) pluck( x, "gamma", .default = 0 ),
    numeric( 1L ),
    USE.NAMES = FALSE
  )
  iovIdx  = which( gamma_all > 0 )
  has_iov = length( iovIdx ) > 0L

  # Column names available in the nested gradient frames.
  sampleGrad = character( 0L )
  if ( length( allGradientsData ) && length( outputNames ) ) {
    sampleGrad = tryCatch( {
      colnames( as.matrix( list_rbind( map(
        outputNames,
        ~ allGradientsData[[ 1L ]]$gradients[[ 1L ]]$gradient[[ .x ]]
      ) ) ) )
    }, error = function( e ) character( 0L ) )
  }

  # IOV noise needs gradients for every gamma>0 parameter (even without eta / IIV).
  iovCols = .indBayesIovGradientCols( iovIdx, sampleGrad, bareNames )
  gradCols = unique( c( feCols, iovCols ) )
  if ( length( sampleGrad ) )
    gradCols = gradCols[ gradCols %in% sampleGrad ]
  if ( !length( gradCols ) )
    gradCols = feCols

  pieces = compact( map( seq_along( allGradientsData ), function( iter ) {
    comboGrad = allGradientsData[[ iter ]]
    nOcc = length( comboGrad$gradients )
    if ( !nOcc )
      return( NULL )
    blocks = .indBayesOccasionBlocks(
      comboGrad, varianceResults[[ iter ]], nOcc, outputNames,
      gradCols, feCols, has_iov, bareNames, gamma_all, muChain, bayesian
    )
    list(
      mat    = .indBayesComboSubjectFim( blocks, model, bayesian ),
      weight = comboGrad$proportion
    )
  } ) )
  mats    = map( pieces, "mat" )
  weights = map_dbl( pieces, "weight" )

  if ( !length( mats ) ) {
    empty = matrix( numeric( 0L ), 0L, 0L )
    return( list(
      fisherMatrix = empty, covariance = empty,
      MFbeta = empty, MFVar = empty
    ) )
  }

  harm = .pfimHarmonicMeanFim( mats, weights )
  # Split FE / residual for Individual callers that still expect two blocks.
  if ( !bayesian && length( mats ) ) {
    n_fe = length( feCols )
    M = harm$fisherMatrix
    if ( nrow( M ) > n_fe ) {
      harm$MFbeta = M[ seq_len( n_fe ), seq_len( n_fe ), drop = FALSE ]
      harm$MFVar  = M[ ( n_fe + 1L ):nrow( M ), ( n_fe + 1L ):ncol( M ), drop = FALSE ]
    } else {
      harm$MFbeta = M
      harm$MFVar  = matrix( numeric( 0L ), 0L, 0L )
    }
  }
  harm
}

#' Flat gradient matrix for Individual/Bayesian FIM (aggregates nested evaluations when needed).
#'
#' @param model A \code{Model} object.
#' @param arm An \code{Arm} object.
#' @return Numeric matrix (observation rows x parameter columns).
#' @keywords internal
getArmEvaluationGradientsMatrix = function( model, arm, pfimproject = NULL, evalModel = NULL ) {
  if ( !is.null( evalModel ) ) {
    model = evalModel
  } else if ( !is.null( pfimproject ) ) {
    model = rebuildEvalModel( pfimproject, finiteDifference = FALSE )
  }
  raw = prop( arm, "evaluationGradients" )
  if ( .isNestedArmEvaluation( raw ) ) {
    aggregateGradientsWithCovariates( model, arm )
  } else if ( is.data.frame( raw ) ) {
    as.matrix( raw )
  } else {
    do.call( rbind, raw ) |> as.matrix()
  }
}

#' Flat residual variance for Individual/Bayesian FIM.
#'
#' @param arm An \code{Arm} object.
#' @return List with \code{errorVariance} and \code{sigmaDerivatives}.
#' @keywords internal
getArmEvaluationVarianceFlat = function( arm ) {
  raw = prop( arm, "evaluationVariance" )
  if ( .isNestedArmEvaluation( raw ) ) aggregateVarianceWithCovariates( arm ) else raw
}

#' Internal beta column names (beta_param_covariate_category) from covariate definitions.
#'
#' @param modelCovariates List of covariate objects.
#' @return Character vector.
#' @noRd
#' @keywords internal
.betaInternalNamesFromCovariates = function( modelCovariates ) {
  if ( length( modelCovariates ) == 0L ) return( character( 0L ) )
  extractNames = function( covList ) {
    if ( length( covList ) == 0L ) return( character( 0L ) )
    covList |>
      map( function( cov ) {
        covName    = prop( cov, "name" )
        categories = prop( cov, "categories" )
        effects    = prop( cov, "effects" )
        ( seq_len( length( categories ) - 1L ) + 1L ) |>
          map( function( icat ) {
            cat = categories[[ icat ]]
            if ( !cat %in% names( effects ) ) return( character( 0L ) )
            paste0(
              "beta_", .namedEffectParams( effects[[ cat ]] ),
              "_", covName, "_", cat
            )
          }) |> unlist( use.names = FALSE )
      }) |> unlist( use.names = FALSE )
  }
  byClass = .splitCovariatesByClass( modelCovariates )
  unique( c(
    extractNames( pluck( byClass, "CategoricalCovariate",        .default = list() ) ),
    extractNames( pluck( byClass, "CategoricalCovariateWithIOV", .default = list() ) )
  ) )
}

#' Drop fixed-\eqn{\mu} columns from a fixed-effects name vector.
#'
#' Gradient matrices may still carry zero columns for \code{fixedMu} parameters;
#' Population already removes them when building \code{MFbeta}. Individual /
#' Bayesian must use the same estimable-only set so the FE block stays
#' non-singular.
#'
#' @param cn Character vector (\code{mu_*}, \code{beta_*}, or bare mu names).
#' @param parameters List of \code{ModelParameter} objects.
#' @return \code{cn} without non-estimable \eqn{\mu} entries (\eqn{\beta} kept).
#' @noRd
#' @keywords internal
.fimDropFixedMuColumns = function( cn, parameters ) {
  if ( !length( cn ) ) return( cn )
  estimable_mu   = .estimableMuNames( parameters, "mu_" )
  estimable_bare = .estimableMuNames( parameters, "" )
  keep = vapply( cn, function( nm ) {
    if ( startsWith( nm, "beta_" ) ) return( TRUE )
    if ( startsWith( nm, "mu_" ) ) return( nm %in% estimable_mu )
    # Bare parameter name (no-covariate Individual / Bayesian path).
    if ( nm %in% estimable_bare ) return( TRUE )
    bare_all = map_chr( parameters, ~ prop( .x, "name" ) )
    if ( nm %in% bare_all ) return( FALSE )
    TRUE
  }, logical( 1L ) )
  cn[ keep ]
}

#' Column names for the fixed-effects block (mu_* and beta_* when covariates/IOV apply).
#'
#' @param model A \code{Model} object.
#' @param arm An \code{Arm} object (required when covariates or IOV are present).
#' @return Character vector of estimable FE columns (subset of
#'   \code{getArmEvaluationGradientsMatrix()} when covariates/IOV apply).
#' @noRd
#' @keywords internal
.fimFixedEffectColumnNames = function( model, arm = NULL, pfimproject = NULL, evalModel = NULL ) {
  if ( !is.null( evalModel ) ) {
    model = evalModel
  } else if ( !is.null( pfimproject ) ) {
    model = rebuildEvalModel( pfimproject, finiteDifference = FALSE )
  }
  parameters = prop( model, "modelParameters" )
  if ( !usesCovariateOccasionStructure( model ) ) {
    .estimableMuNames( parameters )
  } else {
    if ( is.null( arm ) )
      stop( "arm is required to resolve FIM columns with covariates or IOV.", call. = FALSE )
    cn = colnames( getArmEvaluationGradientsMatrix( model, arm, evalModel = model ) )
    if ( is.null( cn ) || !length( cn ) ) {
      cn = c(
        .estimableMuNames( parameters, "mu_" ),
        .betaInternalNamesFromCovariates( prop( model, "modelCovariates" ) )
      )
    }
    .fimDropFixedMuColumns( cn, parameters )
  }
}

#' Fixed-effects labels after \code{run()} (Individual + Bayesian FIM).
#'
#' Prefer column order from the estimable-filtered gradient layout when
#' covariates/IOV are present so SE/RSE rows align with \code{fisherMatrix}.
#' Population \code{setEvaluationFim} builds mu/omega names from parameters
#' directly; beta still comes through this helper.
#' @noRd
#' @keywords internal
.fimFixedEffectLabels = function( evaluation, greek = .greekConsole ) {

  parameters      = prop( evaluation, "modelParameters" )
  modelCovariates = prop( evaluation, "modelCovariates" )
  has_cov         = length( modelCovariates ) > 0L

  evalDesign = pluck( prop( evaluation, "evaluationDesign" ), 1L )
  evalArm    = pluck( prop( evalDesign, "evaluationArms" ), 1L )
  evalModel  = rebuildEvalModel( evaluation, finiteDifference = FALSE )
  hasComplex = usesCovariateOccasionStructure( evalModel )

  feInternal = if ( hasComplex ) {
    .fimFixedEffectColumnNames( evalModel, evalArm, evalModel = evalModel )
  } else {
    character( 0L )
  }

  # Prefer FE column order from the (already estimable-filtered) gradient layout.
  muBare = if ( hasComplex && length( feInternal ) > 0L ) {
    sub( "^mu_", "", feInternal[ startsWith( feInternal, "mu_" ) ] )
  } else {
    .estimableMuNames( parameters, "" )
  }
  columnNamesMu = paste0( greek[ "mu" ], muBare )

  betaInternal = if ( !has_cov ) {
    character( 0L )
  } else if ( length( feInternal ) > 0L ) {
    feInternal[ startsWith( feInternal, "beta_" ) ]
  } else {
    .betaInternalNamesFromCovariates( modelCovariates )
  }

  columnNamesBeta = if ( length( betaInternal ) > 0L ) {
    paste0( greek[ "beta" ], sub( "^beta_", "", betaInternal ) )
  } else {
    character( 0L )
  }

  paramByName = set_names( parameters, map_chr( parameters, ~ prop( .x, "name" ) ) )
  muValues = if ( length( muBare ) ) {
    map_dbl( muBare, function( nm ) .paramDist( paramByName[[ nm ]], "mu" ) )
  } else {
    numeric( 0L )
  }

  betaValues = if ( length( betaInternal ) > 0L ) {
    .betaValuesFromCovariates( modelCovariates, betaInternal )
  } else {
    numeric( 0L )
  }

  list(
    feInternal      = feInternal,
    columnNamesMu   = columnNamesMu,
    columnNamesBeta = columnNamesBeta,
    betaInternal    = betaInternal,
    muValues        = muValues,
    betaValues      = betaValues,
    has_cov         = has_cov,
    hasComplex      = hasComplex,
    evalModel       = evalModel,
    evalArm         = evalArm
  )
}

#' Display names and values for the residual-variance block (Individual FIM).
#'
#' @param evaluation A \code{Evaluation} object.
#' @param greek Named character vector (default \code{.greekConsole}).
#' @return List with \code{columnNamesSigma} and \code{sigmaValues}.
#' @noRd
#' @keywords internal
.fimSigmaBlockLabels = function( evaluation, greek = .greekConsole ) {
  modelError = prop( evaluation, "modelError" )
  list(
    columnNamesSigma = .sigmaNames( modelError, greek[ "sigma" ] ),
    sigmaValues      = .sigmaValues( modelError )
  )
}

#' Diagonal of a Moore-Penrose covariance from a possibly singular FIM.
#'
#' Used when Cholesky inversion fails (e.g. inactive covariate effect when the
#' global covariate value is 0). Directions with no information get \code{Inf}.
#' Does not alter D-criterion / \code{singularFim}: callers set the flag when
#' this path is taken from \code{.fimBuildSeAndRse}.
#' @noRd
#' @keywords internal
.fimPinvCovarianceDiagonal = function( M ) {
  pmax( diag( .pfimPsdPseudoInverse( M ) ), 0 )
}

#' Build SE / RSE data frames from a labelled FIM matrix.
#'
#' @param M Full Fisher matrix with \code{dimnames} set.
#' @param allNames Row/column names (same order as \code{pVals}).
#' @param pVals Parameter values for RSE denominators.
#' @param abs_denominator If \code{TRUE} (default), RSE uses \code{abs(pVals)}.
#' @return List with \code{SE}, \code{RSE}, \code{table}/\code{SEAndRSE}, and
#'   \code{singular} (logical: Cholesky failed, pseudo-inverse used).
#' @noRd
#' @keywords internal
.fimBuildSeAndRse = function( M, allNames, pVals, abs_denominator = TRUE ) {
  # Singular FIM: keep run()/show() alive via pseudo-inverse SE;
  # non-informative directions (e.g. beta_SEX when SEX = 0) get Inf.
  # Contract: length(allNames) == ncol(M) == length(pVals); callers label first.
  singular = FALSE
  SE = tryCatch(
    sqrt( pmax( diag( .safeCholInv( M ) ), 0 ) ),
    error = function( e ) NULL
  )
  if ( is.null( SE ) ) {
    singular = TRUE
    SE = sqrt( .fimPinvCovarianceDiagonal( M ) )
  }
  if ( singular ) {
    n_inf = sum( !is.finite( SE ) )
    warning(
      "Fisher information matrix is singular or not positive definite; ",
      "SE/RSE used a Moore-Penrose pseudo-inverse",
      if ( n_inf > 0L )
        paste0( " (", n_inf, " non-informative direction(s) have Inf SE)" )
      else
        "",
      ". Interpret SE/RSE and design criteria with caution.",
      call. = FALSE
    )
  }
  # abs() keeps RSE finite for negative betas without changing |value| scale.
  den = if ( abs_denominator ) abs( pVals ) else pVals
  RSE = SE / den * 100
  seDF = data.frame( parametersValues = pVals, SE = SE, RSE = RSE )
  rownames( seDF ) = allNames
  list(
    SE       = seDF[ , c( "parametersValues", "SE"  ), drop = FALSE ],
    RSE      = seDF[ , c( "parametersValues", "RSE" ), drop = FALSE ],
    table    = seDF,
    SEAndRSE = seDF,
    singular = singular
  )
}
