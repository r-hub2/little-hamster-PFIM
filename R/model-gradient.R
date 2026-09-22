#' Evaluate model outputs over finite-difference parameter shifts.
#' @param model \code{Model} with precomputed FD shift matrix.
#' @param arm \code{Arm} object used for model evaluation.
#' @return List of model evaluations, one per FD shift column.
#' @noRd
#' @keywords internal
.fdModelEvaluations = function( model, arm ) {
  parameters = prop( model, "modelParameters" )
  gradParams = prop( model, "parametersForComputingGradient" )
  shifted    = gradParams$shifted
  if ( is.null( shifted ) )
    stop( "parametersForComputingGradient$shifted is missing.", call. = FALSE )
  if ( !is.matrix( shifted ) )
    shifted = matrix( shifted, ncol = 1L )

  # Strip covariate/occasion structure: FD shifts act on theta at one occasion.
  workModel = .pfimCloneS7( model )
  prop( workModel, "modelCovariates" ) = list()
  prop( workModel, "numberOfOccasions" ) = 1L

  # Column 1 = nominal; remaining columns = ±ε (or polynomial) parameter shifts.
  # Fresh parameter clones per column (S7 prop<- copies lists; avoid shared mu).
  map( seq_len( ncol( shifted ) ), function( iter ) {
    m = workModel
    prop( m, "modelParameters" ) = .pfimShiftModelParameters( parameters, shifted[ , iter ] )
    evaluateModel( .pfimPrepareModelForEvaluation( m, arm ), arm )
  } )
}

#' Reconstruct parameter gradients from FD evaluations.
#' @param model \code{Model} providing parameter names and FD coefficients.
#' @param evaluations List of model outputs from \code{.fdModelEvaluations()}.
#' @return Named list of gradient data frames by output.
#' @noRd
#' @keywords internal
.gradientsFromFdEvaluations = function( model, evaluations ) {
  parameters  = prop( model, "modelParameters" )
  npar        = length( parameters )
  paramNames  = map_chr( parameters, ~ prop( .x, "name" ) )
  outputNames = prop( model, "outputNames" )
  gradParams  = prop( model, "parametersForComputingGradient" )
  XcolsInv    = gradParams$XcolsInv
  frac        = gradParams$frac
  freeIdx     = gradParams$freeIdx %||% seq_len( npar )
  nFree       = length( freeIdx )

  # Stencil column contract (matches .pfimFdStencilScheme):
  #   linearOnly: evals = [nominal | +e_i | -e_i]  ->  central (f+ - f-) / (2 incr_i)
  #   quadratic:  XcolsInv %*% f / frac; columns 2..(1+nFree) are the first derivatives
  # Parameters not in freeIdx stay as exact zeros. freeIdx includes estimable mu,
  # any omega>0, and any gamma>0 (IOV still needs df/dmu).
  map( outputNames, function( outName ) {
    if ( isTRUE( gradParams$linearOnly ) ) {
      incr  = gradParams$incr
      nRow  = length( evaluations[[ 1L ]][[ outName ]][ , 2L ] )
      gradsFree = vapply(
        seq_len( nFree ),
        function( i ) {
          fp = evaluations[[ 1L + i ]][[ outName ]][ , 2L ]
          fm = evaluations[[ 1L + nFree + i ]][[ outName ]][ , 2L ]
          ( fp - fm ) / ( 2 * incr[ i ] )
        },
        numeric( nRow )
      )
      if ( nFree == 1L )
        gradsFree = matrix( gradsFree, ncol = 1L )
      # One sampling time: vapply drops to a vector; force (nRow x nFree).
      gradsFree = matrix( gradsFree, nrow = nRow, ncol = nFree )
    } else {
      # XcolsInv is NULL when linearOnly; that branch is handled above.
      cols = map( evaluations, ~ .x[[ outName ]][ , 2L ] )
      outputMatrix = t( do.call( cbind, cols ) )

      # frac undoes the incr / incr^2 / cross-product scaling baked into Xcols.
      raw       = XcolsInv %*% outputMatrix / frac
      gradsFree = matrix(
        t( raw )[ , 2L:( 1L + nFree ), drop = FALSE ],
        ncol = nFree
      )
    }

    grads = matrix( 0, nrow = nrow( gradsFree ), ncol = npar )
    grads[ , freeIdx ] = gradsFree
    as.data.frame( grads ) |> stats::setNames( paramNames )
  } ) |> set_names( outputNames )
}

#' Compute FD gradients and nominal evaluation together.
#' @param model \code{Model} configured for gradient computation.
#' @param arm \code{Arm} object used for model evaluation.
#' @return List with \code{gradTheta} and \code{nominalEvaluation}.
#' @noRd
#' @keywords internal
.evaluateModelGradientCoreInner = function( model, arm ) {
  # column 0 of the FD grid is the nominal evaluation point
  outputNames = prop( model, "outputNames" )
  evals       = .fdModelEvaluations( model, arm )
  list(
    gradTheta         = .gradientsFromFdEvaluations( model, evals ),
    nominalEvaluation = map( outputNames, ~ evals[[ 1L ]][[ .x ]] ) |> set_names( outputNames )
  )
}

#' Evaluate gradient core without covariate expansion.
#' @param model \code{Model} configured for gradient computation.
#' @param arm \code{Arm} object used for model evaluation.
#' @return Named list of gradient data frames by output.
#' @keywords internal
evaluateModelGradientCore = function( model, arm ) {
  .evaluateModelGradientCoreInner( model, arm )$gradTheta
}

#' Dispatch gradient evaluation: simple path or covariatexoccasion expansion.
#' @return Gradient structure for the model, with covariate expansion when needed.
#' @name evaluateModelGradient
#' @keywords internal
method( evaluateModelGradient, Model ) = function( model, arm ) {
  if ( !usesCovariateOccasionStructure( model ) ) {
    evaluateModelGradientCore( model, arm )
  } else {
    evaluateModelGradientWithCovariates( model, arm, evaluateModelGradientCore )
  }
}

#' Map over data-frame columns with names preserved.
#' @param .x Data frame or matrix-like object.
#' @param .f Function applied to each column and its name.
#' @return Data frame assembled from mapped columns.
#' @noRd
#' @keywords internal
.imap_dfc = function( .x, .f ) {
  .x = as.data.frame( .x, check.names = FALSE, stringsAsFactors = FALSE )
  nm = names( .x )
  if ( is.null( nm ) ) nm = paste0( "V", seq_len( ncol( .x ) ) )
  cols = Map( .f, .x, nm )
  as.data.frame( cols, check.names = FALSE, stringsAsFactors = FALSE )
}

#' Evaluate gradients with covariate and occasion structure.
#' @param model \code{Model} containing covariate combinations.
#' @param arm \code{Arm} object used for evaluation.
#' @param evaluateModelGradientCore Function evaluating gradients for one occasion model.
#' @return Nested list of gradients by combination and occasion.
#' @keywords internal
evaluateModelGradientWithCovariates = function( model, arm, evaluateModelGradientCore ) {
  .evaluateCovariateOccasions(
    model, arm, evaluateModelGradientCore = evaluateModelGradientCore
  )$evaluationGradients
}

#' Precompute shared inputs for covariate gradient evaluation.
#' @param model \code{Model} object with covariate definitions.
#' @return List of cached values used by occasion gradient routines.
#' @noRd
#' @keywords internal
.covariateGradientSetup = function( model ) {
  covariatesEffect  = prop( model, "covariatesEffect" )
  modelCovariatesEq = prop( model, "modelCovariatesEquation" )
  modelParameters   = prop( model, "modelParameters" )
  modelCovariates   = prop( model, "modelCovariates" )

  paramNames = map_chr( modelParameters, ~ prop( .x, "name" ) )
  muValues   = map_dbl( modelParameters, ~ prop( prop( .x, "distribution" ), "mu" ) ) |>
    set_names( paramNames )

  covEqType = if ( S7::S7_inherits( modelCovariatesEq, Additive ) ) modelCovEqAdditive
  else                                              modelCovEqExponential

  covariateWithoutIov = .filterCovariatesByClass( modelCovariates, "CategoricalCovariate"        )
  covariateWithIov    = .filterCovariatesByClass( modelCovariates, "CategoricalCovariateWithIOV" )

  list(
    covariatesEffect      = covariatesEffect,
    paramNames            = paramNames,
    muValues              = muValues,
    covEqType             = covEqType,
    covariateWithoutIov   = covariateWithoutIov,
    covariateWithIov      = covariateWithIov,
    betaList              = .buildBetaList( covariateWithoutIov, covariateWithIov )
  )
}

#' Check whether a beta effect applies to an occasion.
#' @param betaInfo List describing beta-covariate mapping.
#' @param combinationInfo Named list from parsed combination name.
#' @param occasionIndex Integer occasion index.
#' @param setup Setup list returned by \code{.covariateGradientSetup()}.
#' @return Logical scalar indicating activity for the occasion.
#' @noRd
#' @keywords internal
.isBetaActiveForOccasion = function( betaInfo, combinationInfo, occasionIndex, setup ) {
  if ( !betaInfo$isIOV ) {
    identical( combinationInfo[[ betaInfo$covName ]], betaInfo$category )
  } else {
    seqLabel = combinationInfo[[ betaInfo$covName ]]
    if ( is.null( seqLabel ) ) return( FALSE )
    cov      = setup$covariateWithIov[[ betaInfo$covIndex ]]
    seqs     = prop( cov, "sequences" )
    seqNames = names( seqs ) %||% paste0( "sequence_", seq_along( seqs ) )
    seqIdx   = match( seqLabel, seqNames )
    if ( is.na( seqIdx ) ) return( FALSE )
    seqValues = seqs[[ seqIdx ]]
    if ( occasionIndex > length( seqValues ) ) return( FALSE )
    identical( seqValues[[ occasionIndex ]], betaInfo$category )
  }
}

#' Build gradient matrix for one combination and occasion.
#' @param model \code{Model} object defining outputs and covariates.
#' @param arm \code{Arm} object used in model evaluation.
#' @param tempModel Temporary \code{Model} configured for one occasion.
#' @param setup Setup list returned by \code{.covariateGradientSetup()}.
#' @param combinationInfo Named list of covariate levels for the combination.
#' @param occIndex Integer occasion index.
#' @param occasionParams Named numeric vector of occasion parameters.
#' @param evaluateModelGradientCore Gradient core evaluator function.
#' @param reuseFdNominal Logical scalar; include nominal output from FD run.
#' @return Gradient list, or list with gradient and nominal evaluation.
#' @noRd
#' @keywords internal
.gradientWithCovariatesForOccasion = function(
    model, arm, tempModel, setup, combinationInfo, occIndex, occasionParams,
    evaluateModelGradientCore, reuseFdNominal = FALSE ) {

  # reuseFdNominal: FD column 1 is the nominal evaluation - avoid a second
  # evaluateModel() when both gradient and model output are requested.
  if ( reuseFdNominal ) {
    fdOut     = .evaluateModelGradientCoreInner( tempModel, arm )
    gradTheta = fdOut$gradTheta
    nominal   = fdOut$nominalEvaluation
  } else {
    gradTheta = evaluateModelGradientCore( tempModel, arm )
    nominal   = NULL
  }

  # FD grads are df/dtheta at the occasion theta. Convert to df/dmu via the
  # covariate link (not LogNormal/Normal adjustGradient - that eta-scale lives
  # in the FIM assembly: .pfimPopMuChainFactors).
  # Both Additive (theta = mu(1+Sumbeta)) and Exponential (theta = mu exp(beta'z)) share:
  #   df/dmu = (theta/mu) df/dtheta
  # Using theta/mu from occasionParams avoids the old Pi(1+beta) additive bug and the
  # IOV effects[[category]] lookup that returned NULL on sequence-indexed lists.
  gradMu = map( gradTheta, function( gradDf ) {
    .imap_dfc( gradDf, function( gradVec, pName ) {
      mu_i    = setup$muValues[[ pName ]]
      theta_i = occasionParams[[ pName ]]
      if ( is.null( mu_i ) || mu_i == 0 || is.na( mu_i ) ) gradVec
      else gradVec * ( theta_i / mu_i )
    })
  })

  outputNames = prop( model, "outputNames" )
  gradient = map( outputNames, function( outName ) {
    gradMu_out    = gradMu[[ outName ]]
    gradTheta_out = gradTheta[[ outName ]]
    nRows         = nrow( gradMu_out )

    muCols = map( setup$paramNames, ~ gradMu_out[[ .x ]] ) |>
      set_names( paste0( "mu_", setup$paramNames ) )

    # Inactive betas are exact zeros (reference category / wrong occasion).
    # Chain rule uses dtheta/dbeta: exp link -> theta; additive -> mu.
    betaCols = imap( setup$betaList, function( betaInfo, betaName ) {
      if ( !.isBetaActiveForOccasion( betaInfo, combinationInfo, occIndex, setup ) )
        return( rep( 0, nRows ) )
      pName = betaInfo$param
      if ( setup$covEqType == modelCovEqExponential )
        gradTheta_out[[ pName ]] * occasionParams[[ pName ]]
      else
        gradTheta_out[[ pName ]] * setup$muValues[[ pName ]]
    })

    as.data.frame( c( muCols, betaCols ) )
  }) |> set_names( outputNames )

  if ( reuseFdNominal && !is.null( nominal ) )
    return( list( gradient = gradient, nominalEvaluation = nominal ) )
  gradient
}

#' Evaluate one covariate-occasion block with caching.
#' @param model \code{Model} object to evaluate.
#' @param arm \code{Arm} object providing administration and sampling data.
#' @param combinationName Character combination identifier.
#' @param occasion Character occasion label.
#' @param occasionParams Named numeric vector of occasion-specific parameters.
#' @param occIndex Integer occasion index.
#' @param combinationInfo Parsed combination metadata.
#' @param baseModelParameters List of base \code{ModelParameter} objects.
#' @param setup Gradient setup list or \code{NULL}.
#' @param evaluateModelCore Optional model core evaluator.
#' @param evaluateModelGradientCore Optional gradient core evaluator.
#' @param wantModel Logical scalar requesting model outputs.
#' @param wantGrad Logical scalar requesting gradients.
#' @return List containing requested evaluation and/or gradient fields.
#' @noRd
#' @keywords internal
.evalCovOccasion = function(
    model, arm, combinationName, occasion, occasionParams, occIndex,
    combinationInfo, baseModelParameters, setup,
    evaluateModelCore, evaluateModelGradientCore, wantModel, wantGrad ) {

  # Cache key: arm + combination + occasion + occasionParams (+ want flags).
  # Hit skips both FD and model evaluation for this (combo, occ) node.
  cached = .pfimCovOccasionCacheGet(
    arm, combinationName, occasion, occasionParams, wantModel, wantGrad, model
  )
  if ( !is.null( cached ) )
    return( cached )

  occ  = list( occasion = occasion )
  both = wantModel && wantGrad

  if ( both ) {
    # One FD pass: gradient + nominal (reuseFdNominal); variance from that f.
    tempModel = .modelForCovariateOccasion(
      model, arm, baseModelParameters, occasionParams, needFd = TRUE
    )
    gradOut = .gradientWithCovariatesForOccasion(
      model, arm, tempModel, setup, combinationInfo, occIndex, occasionParams,
      evaluateModelGradientCore, reuseFdNominal = TRUE
    )
    occ$gradient   = gradOut$gradient
    occ$evaluation = gradOut$nominalEvaluation
    occ$variance   = evaluateModelVarianceCore( model, gradOut$nominalEvaluation )
  } else {
    if ( wantGrad ) {
      tempModel = .modelForCovariateOccasion(
        model, arm, baseModelParameters, occasionParams, needFd = TRUE
      )
      occ$gradient = .gradientWithCovariatesForOccasion(
        model, arm, tempModel, setup, combinationInfo, occIndex, occasionParams,
        evaluateModelGradientCore
      )
    }
    if ( wantModel ) {
      tempModel = .modelForCovariateOccasion(
        model, arm, baseModelParameters, occasionParams
      )
      occ$evaluation = evaluateModelCore( tempModel, arm )
      occ$variance   = evaluateModelVarianceCore( model, occ$evaluation )
    }
  }

  .pfimCovOccasionCacheSet(
    arm, combinationName, occasion, occasionParams, occ, wantModel, wantGrad, model
  )
  occ
}

#' Extract a nested covariate-occasion field structure.
#' @param perCombination List of per-combination evaluation records.
#' @param field Character field name to extract from each occasion node.
#' @param subfield Character name for the nested output list.
#' @return List grouped by combination with selected occasion fields.
#' @noRd
#' @keywords internal
.covOccasionField = function( perCombination, field, subfield ) {
  map( perCombination, function( x ) {
    nested = map( x$byOccasion, function( occ ) {
      out = list( occasion = occ$occasion )
      out[[ field ]] = occ[[ field ]]
      out
    })
    out = list( combination = x$combination, proportion = x$proportion )
    out[[ subfield ]] = nested
    out
  })
}

#' Evaluate model and gradients across covariate occasions.
#' @param model \code{Model} object with covariate combination data.
#' @param arm \code{Arm} object used for all evaluations.
#' @param evaluateModelCore Optional function evaluating model outputs.
#' @param evaluateModelGradientCore Optional function evaluating gradients.
#' @return List containing model, variance, and/or gradient results.
#' @noRd
#' @keywords internal
.evaluateCovariateOccasions = function(
    model, arm, evaluateModelCore = NULL, evaluateModelGradientCore = NULL ) {

  wantModel = !is.null( evaluateModelCore )
  wantGrad  = !is.null( evaluateModelGradientCore )
  if ( !wantModel && !wantGrad )
    stop( "evaluateModelCore or evaluateModelGradientCore required.", call. = FALSE )

  combinations        = prop( model, "covariatesCombination" )$combinations
  modelParamsWithCov  = prop( model, "modelParametersWithCovariates" )
  baseModelParameters = prop( model, "modelParameters" )
  modelCovariates     = prop( model, "modelCovariates" )
  setup               = if ( wantGrad ) .covariateGradientSetup( model ) else NULL

  perCombination = imap( combinations$name, function( combinationName, i ) {
    paramsForCombo  = modelParamsWithCov[[ combinationName ]]
    combinationInfo = parseCombinationName( combinationName, modelCovariates )

    evaluationByOccasion = imap( paramsForCombo, function( occasionParams, occasion ) {
      .evalCovOccasion(
        model, arm, combinationName, occasion, occasionParams,
        as.integer( str_remove( occasion, "^occ" ) ), combinationInfo,
        baseModelParameters, setup, evaluateModelCore, evaluateModelGradientCore,
        wantModel, wantGrad
      )
    })

    list(
      combination = combinationName,
      proportion  = combinations$proportion[[ i ]],
      byOccasion  = evaluationByOccasion
    )
  })

  list(
    evaluationModel     = if ( wantModel )     .covOccasionField( perCombination, "evaluation", "evaluations" )     else NULL,
    evaluationVariance  = if ( wantModel )     .covOccasionField( perCombination, "variance",   "variances" )       else NULL,
    evaluationGradients = if ( wantGrad )      .covOccasionField( perCombination, "gradient",   "gradients" )       else NULL
  )
}

#' Select base evaluator for covariate occasion computations.
#' @param model \code{Model} object whose class determines evaluator.
#' @return Function for core evaluation, or \code{NULL} when unsupported.
#' @noRd
#' @keywords internal
.covariateEvaluationCore = function( model ) {
  if ( S7::S7_inherits( model, ModelAnalyticInfusionSteadyState ) )
    return( evaluateAnalyticInfusionSteadyStateCore )
  if ( S7::S7_inherits( model, ModelAnalyticInfusion ) )
    return( evaluateAnalyticInfusionCore )
  if ( S7::S7_inherits( model, ModelAnalyticSteadyState ) )
    return( evaluateAnalyticSteadyStateCore )
  if ( S7::S7_inherits( model, ModelAnalytic ) )
    return( evaluateAnalyticCore )
  if ( S7::S7_inherits( model, ModelODEInfusionDoseInEquation ) )
    return( .odeInfusionEvaluateModelCore )
  if ( S7::S7_inherits( model, ModelODEBolus ) ||
       S7::S7_inherits( model, ModelODEDoseInEquations ) ||
       S7::S7_inherits( model, ModelODEDoseNotInEquations ) )
    return( .odeEvaluateModelCore )
  NULL
}

#' Test availability of a covariate evaluation core.
#' @param model \code{Model} object to test.
#' @return Logical scalar indicating evaluator availability.
#' @noRd
#' @keywords internal
.hasCovariateEvaluationCore = function( model ) !is.null( .covariateEvaluationCore( model ) )

# Beta FIM columns: one per (covariate, non-reference category, affected parameter).
# categories[[1]] is the reference (no beta); IOV betas activate per occasion sequence.

#' Name beta metadata entries for gradient column assembly.
#' @param entries List of \code{list(name = ..., value = ...)} beta descriptors.
#' @return Named list of beta metadata values.
#' @noRd
#' @keywords internal
.nameBetaList = function( entries ) {
  if ( !length( entries ) ) return( list() )
  set_names( map( entries, "value" ), map_chr( entries, "name" ) )
}

#' Build beta-column metadata from declared covariate effects.
#'
#' Uses the user \code{effects} list (including zeros). Overlay vectors pad
#' every parameter with 0 and must not be used to decide which \eqn{\beta}
#' columns exist.
#' @param covariateWithoutIov List of \code{CategoricalCovariate} objects.
#' @param covariateWithIov List of \code{CategoricalCovariateWithIOV} objects.
#' @return Named list describing beta columns and activation rules.
#' @noRd
#' @keywords internal
.buildBetaList = function( covariateWithoutIov, covariateWithIov ) {

  betasWithoutIov = if ( length( covariateWithoutIov ) > 0L ) {
    covariateWithoutIov |>
      map( function( cov ) {
        covName    = prop( cov, "name" )
        categories = prop( cov, "categories" )
        declared   = prop( cov, "effects" )

        if ( length( categories ) <= 1L ) return( list() )

        ( seq_len( length( categories ) - 1L ) + 1L ) |>
          map( function( icat ) {
            cat = categories[[ icat ]]
            map( .namedEffectParams( declared[[ cat ]] ), ~ list(
              name  = paste0( "beta_", .x, "_", covName, "_", cat ),
              value = list( covName  = covName,
                            category = cat,
                            param    = .x,
                            isIOV    = FALSE )
            ))
          }) |> list_flatten()
      }) |> list_flatten() |>
      .nameBetaList()
  } else list()

  betasWithIov = if ( length( covariateWithIov ) > 0L ) {
    imap( covariateWithIov, function( cov, covIndex ) {
      covName    = prop( cov, "name" )
      categories = prop( cov, "categories" )
      declared   = prop( cov, "effects" )

      if ( length( categories ) <= 1L ) return( list() )

      ( seq_len( length( categories ) - 1L ) + 1L ) |>
        map( function( icat ) {
          category = categories[[ icat ]]
          map( .namedEffectParams( declared[[ category ]] ), ~ list(
            name  = paste0( "beta_", .x, "_", covName, "_", category ),
            value = list(
              covName    = covName,
              covIndex   = covIndex,
              category   = category,
              categories = categories,
              param      = .x,
              isIOV      = TRUE
            )
          ))
        }) |> list_flatten()
      }) |> list_flatten() |>
      .nameBetaList()
  } else list()

  c( betasWithoutIov, betasWithIov )
}
