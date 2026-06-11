#' isNestedArmEvaluation
#' @param x \code{evaluationGradients} or \code{evaluationVariance} from an arm.
#' @return Logical scalar.
#' @keywords internal
.isNestedArmEvaluation = function( x ) {
  length( x ) > 0L && is.list( x[[1L]] ) &&
    ( !is.null( x[[1L]]$gradients ) || !is.null( x[[1L]]$variances ) )
}

#' numeric vector of covariate effect sizes
#'
#' @param modelCovariates List of covariate objects.
#' @return Named numeric vector.
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
          names( eff )[ eff != 0 ] |>
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
#' @keywords internal
.betaValuesFromCovariates = function( modelCovariates, betaInternal ) {
  if ( length( betaInternal ) == 0L ) return( numeric( 0L ) )
  dict = .betaDictFromCovariates( modelCovariates )
  vals = map_dbl( betaInternal, ~ dict[[ .x ]] )
  vals
}

#' Collapse nested combination-by-occasion gradients to one matrix (all outputs stacked).
#'
#' @param model A \code{Model} object.
#' @param arm An \code{Arm} object with nested \code{evaluationGradients}.
#' @return Numeric matrix (stacked outputs x parameters).
#' @keywords internal
aggregateGradientsWithCovariates = function( model, arm ) {

  allGradientsData = prop( arm, "evaluationGradients" )
  outputNames      = prop( model, "outputNames" )

  if ( length( outputNames ) == 0L && length( allGradientsData ) > 0L ) {
    outputNames = names( allGradientsData[[1L]]$gradients[[1L]]$gradient )
  }

  gradDfs = map( outputNames, function( outName ) {
    zeroGrad = allGradientsData[[1L]]$gradients[[1L]]$gradient[[ outName ]] * 0

    matList = map( allGradientsData, function( combinationData ) {
      nOcc     = length( combinationData$gradients )
      meanGrad = reduce( map( combinationData$gradients, ~ .x$gradient[[ outName ]] ), `+` ) / nOcc
      combinationData$proportion * meanGrad
    })

    reduce( matList, `+`, .init = zeroGrad )
  })

  mat = list_rbind( gradDfs ) |> as.matrix()
  if ( !is.null( gradDfs[[1L]] ) && ncol( gradDfs[[1L]] ) > 0L )
    colnames( mat ) = colnames( gradDfs[[1L]] )
  mat
}

#' Expectation of the residual variance structure over covariate combinations and occasions.
#'
#' @param arm An \code{Arm} object with nested \code{evaluationVariance}.
#' @return List with \code{errorVariance} and \code{sigmaDerivatives}.
#' @keywords internal
aggregateVarianceWithCovariates = function( arm ) {

  varianceResults = prop( arm, "evaluationVariance" )

  meanOccasionVariance = function( occVars ) {
    n = length( occVars )
    err = Reduce( `+`, map( occVars, ~ as.matrix( .x$variance$errorVariance ) ) ) / n
    sig = occVars[[1L]]$variance$sigmaDerivatives
    if ( n > 1L ) {
      sig = map( seq_along( sig ), function( k ) {
        Reduce( `+`, map( occVars, ~ .x$variance$sigmaDerivatives[[ k ]] ) ) / n
      })
    }
    list( errorVariance = err, sigmaDerivatives = sig )
  }

  out = NULL
  for ( combo in varianceResults ) {
    p = combo$proportion
    v = meanOccasionVariance( combo$variances )
    if ( is.null( out ) ) {
      out = list(
        errorVariance      = p * v$errorVariance,
        sigmaDerivatives   = map( v$sigmaDerivatives, ~ p * .x )
      )
    } else {
      out$errorVariance    = out$errorVariance + p * v$errorVariance
      out$sigmaDerivatives = map2( out$sigmaDerivatives, v$sigmaDerivatives, `+` )
    }
  }
  out
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
              "beta_", names( effects[[ cat ]] )[ effects[[ cat ]] != 0 ],
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

#' Column names for the fixed-effects block (mu_* and beta_* when covariates/IOV apply).
#'
#' @param model A \code{Model} object.
#' @param arm An \code{Arm} object (required when covariates or IOV are present).
#' @return Character vector matching columns of \code{getArmEvaluationGradientsMatrix()}.
#' @keywords internal
.fimFixedEffectColumnNames = function( model, arm = NULL, pfimproject = NULL, evalModel = NULL ) {
  if ( !is.null( evalModel ) ) {
    model = evalModel
  } else if ( !is.null( pfimproject ) ) {
    model = rebuildEvalModel( pfimproject, finiteDifference = FALSE )
  }
  if ( !usesCovariateOccasionStructure( model ) ) {
    prop( model, "modelParameters" ) |>
      keep( ~ !isTRUE( prop( .x, "fixedMu" ) ) ) |>
      keep( ~ prop( prop( .x, "distribution" ), "mu" ) != 0 ) |>
      map_chr( ~ prop( .x, "name" ) )
  } else {
    if ( is.null( arm ) )
      stop( "arm is required to resolve FIM columns with covariates or IOV.", call. = FALSE )
    cn = colnames( getArmEvaluationGradientsMatrix( model, arm, evalModel = model ) )
    if ( !is.null( cn ) && length( cn ) > 0L ) return( cn )
    c(
      paste0( "mu_", prop( model, "modelParameters" ) |>
                keep( ~ !isTRUE( prop( .x, "fixedMu" ) ) ) |>
                keep( ~ prop( prop( .x, "distribution" ), "mu" ) != 0 ) |>
                map_chr( ~ prop( .x, "name" ) ) ),
      .betaInternalNamesFromCovariates( prop( model, "modelCovariates" ) )
    )
  }
}

#' Display names and values for the fixed-effects block after \code{run()}.
#'
#' Shared by \code{IndividualFim} and \code{BayesianFim} \code{setEvaluationFim}.
#'
#' @param evaluation A \code{Evaluation} object after \code{run()}.
#' @param greek Named character vector (default \code{.greekConsole}).
#' @return List with \code{feInternal}, \code{columnNamesMu}, \code{columnNamesBeta},
#'   \code{betaInternal}, \code{muValues}, \code{betaValues}, \code{has_cov}, \code{hasComplex}.
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

  columnNamesMu = if ( hasComplex && length( feInternal ) > 0L ) {
    paste0( greek[ "mu" ], sub( "^mu_", "", feInternal[ startsWith( feInternal, "mu_" ) ] ) )
  } else {
    parameters |>
      keep( ~ !isTRUE( prop( .x, "fixedMu" ) ) ) |>
      keep( ~ prop( prop( .x, "distribution" ), "mu" ) != 0 ) |>
      map_chr( ~ paste0( greek[ "mu" ], prop( .x, "name" ) ) )
  }

  betaInternal = if ( has_cov ) {
    if ( length( feInternal ) > 0L ) {
      feInternal[ startsWith( feInternal, "beta_" ) ]
    } else {
      .betaInternalNamesFromCovariates( modelCovariates )
    }
  } else {
    character( 0L )
  }

  columnNamesBeta = if ( length( betaInternal ) > 0L ) {
    paste0( greek[ "beta" ], sub( "^beta_", "", betaInternal ) )
  } else {
    character( 0L )
  }

  muValues = parameters |>
    keep( ~ !isTRUE( prop( .x, "fixedMu" ) ) ) |>
    keep( ~ prop( prop( .x, "distribution" ), "mu" ) != 0 ) |>
    map_dbl( ~ prop( prop( .x, "distribution" ), "mu" ) )

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
#' @keywords internal
.fimSigmaBlockLabels = function( evaluation, greek = .greekConsole ) {

  modelError = prop( evaluation, "modelError" )

  columnNamesSigma = modelError |>
    map( function( err ) {
      out = prop( err, "output" )
      c(
        if ( prop( err, "sigmaInter" ) != 0 && !prop( err, "sigmaInterFixed" ) )
          paste0( greek[ "sigma" ], "_inter_", out ),
        if ( prop( err, "sigmaSlope"  ) != 0 && !prop( err, "sigmaSlopeFixed"  ) )
          paste0( greek[ "sigma" ], "_slope_", out )
      )
    }) |> unlist( use.names = FALSE )

  sigmaValues = modelError |>
    map( function( err ) {
      v = c()
      if ( prop( err, "sigmaInter" ) != 0 && !prop( err, "sigmaInterFixed" ) )
        v = c( v, prop( err, "sigmaInter" ) )
      if ( prop( err, "sigmaSlope"  ) != 0 && !prop( err, "sigmaSlopeFixed"  ) )
        v = c( v, prop( err, "sigmaSlope" ) )
      v
    }) |> unlist( use.names = FALSE )

  list(
    columnNamesSigma = columnNamesSigma,
    sigmaValues      = sigmaValues
  )
}

#' Build SE / RSE data frames from a labelled FIM matrix.
#'
#' @param M Full Fisher matrix with \code{dimnames} set.
#' @param allNames Row/column names (same order as \code{pVals}).
#' @param pVals Parameter values for RSE denominators.
#' @param abs_denominator If \code{TRUE}, use \code{abs(pVals)} (Bayesian).
#' @return List with \code{SE}, \code{RSE}, and \code{SEAndRSE} data frames.
#' @keywords internal
.fimBuildSeAndRse = function( M, allNames, pVals, abs_denominator = FALSE ) {
  SE  = sqrt( diag( .safeCholInv( M ) ) )
  den = if ( abs_denominator ) abs( pVals ) else pVals
  RSE = SE / den * 100
  seDF = data.frame( parametersValues = pVals, SE = SE, RSE = RSE )
  rownames( seDF ) = allNames
  list(
    SE       = seDF[ , c( "parametersValues", "SE"  ), drop = FALSE ],
    RSE      = seDF[ , c( "parametersValues", "RSE" ), drop = FALSE ],
    table    = seDF,
    SEAndRSE = seDF
  )
}
