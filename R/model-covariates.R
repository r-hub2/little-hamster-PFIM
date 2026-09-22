# Covariate combinations and FD Hessian for population FIM.
NULL

#' Evaluate model outputs over every covariate combination and occasion.
#' @return Nested list of evaluations keyed by combination and occasion.
#' @name evaluateModelWithCovariates
#' @keywords internal
method( evaluateModelWithCovariates, Model ) = function( model, arm, evaluateModelCore ) {
  .evaluateCovariateOccasions( model, arm, evaluateModelCore = evaluateModelCore )$evaluationModel
}

#' Report whether the model defines any covariates.
#' @param model First argument of generic.
#' @return Logical scalar indicating whether covariates are configured.
#' @name hasCovariates
#' @keywords internal
method( hasCovariates, Model ) = function( model ) {
  length( prop( model, "modelCovariates" ) ) > 0L
}

#' Prepare all covariate-derived structures used by evaluation and the FIM.
#'
#' Pipeline: effect vectors -> combination grid -> occasion-specific mus -> omega.
#' @param model First argument of generic.
#' @return \code{Model} with covariate effects, combinations, and omega matrix prepared.
#' @name defineCovariatesData
#' @keywords internal
method( defineCovariatesData, Model ) = function( model ) {
  model |>
    evaluateCovariatesEffects()         |>
    generateCovariatesCombination()     |>
    modelParametersWithCovariates()     |>
    evaluateOmegaMatrixFromCovariates()
}

#' Build per-covariate effect vectors aligned with model parameter names.
#' @param model First argument of generic.
#' @return \code{Model} with grouped covariate effect vectors.
#' @name evaluateCovariatesEffects
#' @keywords internal
method( evaluateCovariatesEffects, Model ) = function( model ) {
  covariates = prop( model, "modelCovariates" )

  if ( length( covariates ) == 0L ) {
    prop( model, "covariatesEffect" ) = list()
    return( model )
  }

  zeroEffectVector = prop( model, "modelParameters" ) |>
    map_chr( ~ prop( .x, "name" ) ) |>
    set_names() |>
    map_dbl( ~ 0 )

  covariatesEffects = imap( .splitCovariatesByClass( covariates ), function( covList, className ) {
    effects = map( covList, ~ getCovariateEffects( .x, zeroEffectVector ) )
    set_names( effects, map_chr( covList, ~ prop( .x, "name" ) ) )
  })

  prop( model, "covariatesEffect" ) = covariatesEffects
  model
}

#' Build the Cartesian product of categorical and IOV covariate levels.
#'
#' Proportions multiply across covariates; names encode \code{cov=level} pairs.
#' @param model First argument of generic.
#' @return \code{Model} with covariate combination grid and proportions.
#' @name generateCovariatesCombination
#' @keywords internal
method( generateCovariatesCombination, Model ) = function( model ) {
  covariates = prop( model, "modelCovariates" )

  if ( length( covariates ) == 0L ) {
    prop( model, "covariatesCombination" ) = list(
      combinations = data.frame(
        combinationIndex               = 1L,
        covariateWithoutIovCombination = 1L,
        covariateWithIovCombination    = 1L,
        proportion                     = 1,
        name                           = "Reference",
        stringsAsFactors               = FALSE
      ),
      covariateWithoutIovGridValues = NULL,
      covariateWithIovGridValues    = NULL
    )
    return( model )
  }

  covariateWithoutIov = .filterCovariatesByClass( covariates, "CategoricalCovariate"        )
  covariateWithIov    = .filterCovariatesByClass( covariates, "CategoricalCovariateWithIOV" )

  # Build an index grid (data.frame) for one list of covariates.
  makeIndexGrid = function( covList, propName ) {
    if ( length( covList ) == 0L ) return( NULL )
    grids = covList |>
      map( ~ seq_along( prop( .x, propName ) ) ) |>
      set_names( map_chr( covList, ~ prop( .x, "name" ) ) )
    do.call( expand.grid, c( grids, stringsAsFactors = FALSE ) )
  }

  covWithoutIovGrid = makeIndexGrid( covariateWithoutIov, "categories" )
  covWithIovGrid    = makeIndexGrid( covariateWithIov,    "sequences"  )

  hasWithoutIov = !is.null( covWithoutIovGrid )
  hasWithIov    = !is.null( covWithIovGrid    )

  # Unified cross-product: absent type contributes a single pseudo-index 1L.
  withoutIovIdx = if ( hasWithoutIov ) seq_len( nrow( covWithoutIovGrid ) ) else 1L
  withIovIdx    = if ( hasWithIov    ) seq_len( nrow( covWithIovGrid    ) ) else 1L

  fullCombinations = expand.grid(
    covariateWithoutIovIndex = withoutIovIdx,
    covariateWithIovIndex    = withIovIdx,
    stringsAsFactors         = FALSE
  )

  combinations = pmap( fullCombinations, function( covariateWithoutIovIndex,
                                                    covariateWithIovIndex ) {
    dataWithoutIov = if ( hasWithoutIov ) {
      imap( covariateWithoutIov, function( cov, i ) {
        catIdx = covWithoutIovGrid[ covariateWithoutIovIndex, i ]
        list(
          prop = prop( cov, "categoriesProportions" )[[ catIdx ]],
          key  = prop( cov, "name" ),
          val  = prop( cov, "categories"            )[[ catIdx ]]
        )
      })
    } else list()

    dataWithIov = if ( hasWithIov ) {
      imap( covariateWithIov, function( cov, i ) {
        seqIdx   = covWithIovGrid[ covariateWithIovIndex, i ]
        seqs     = prop( cov, "sequences"            )
        seqProps = prop( cov, "sequencesProportions" )
        seqNames = if ( is.null( names( seqs ) ) ) paste0( "sequence_", seq_along( seqs ) ) else names( seqs )
        list(
          prop = seqProps[[ seqIdx ]],
          key  = prop( cov, "name"    ),
          val  = seqNames[[ seqIdx ]]
        )
      })
    } else list()

    allData         = c( dataWithoutIov, dataWithIov )
    proportion      = reduce( map_dbl( allData, "prop" ), `*`, .init = 1 )
    nameParts       = set_names( map_chr( allData, "val" ), map_chr( allData, "key" ) )
    combinationName = .encodeCombinationName( as.list( nameParts ) )

    data.frame(
      combinationIndex               = NA_integer_,
      covariateWithoutIovCombination = covariateWithoutIovIndex,
      covariateWithIovCombination    = covariateWithIovIndex,
      proportion                     = proportion,
      name                           = combinationName,
      stringsAsFactors               = FALSE
    )
  }) |> list_rbind()

  combinations$combinationIndex = seq_len( nrow( combinations ) )

  prop( model, "covariatesCombination" ) = list(
    combinations                  = combinations,
    covariateWithoutIovGridValues = covWithoutIovGrid,
    covariateWithIovGridValues    = covWithIovGrid
  )
  model
}

#' Compute occasion-specific typical values for every covariate combination.
#'
#' Applies additive or exponential covariate links to mu for each occasion.
#' @param model First argument of generic.
#' @return \code{Model} with occasion-specific parameter sets by combination.
#' @name modelParametersWithCovariates
#' @keywords internal
method( modelParametersWithCovariates, Model ) = function( model ) {

  modelCovariatesEquation = prop( model, "modelCovariatesEquation" )
  covariatesEffect        = prop( model, "covariatesEffect"        )
  modelParameters         = prop( model, "modelParameters"         )

  muValues = map_dbl( modelParameters, ~ prop( prop( .x, "distribution" ), "mu" ) ) |>
    set_names( map_chr( modelParameters, ~ prop( .x, "name" ) ) )

  combinations      = prop( model, "covariatesCombination" )$combinations
  covWithoutIovGrid = prop( model, "covariatesCombination" )$covariateWithoutIovGridValues
  covWithIovGrid    = prop( model, "covariatesCombination" )$covariateWithIovGridValues

  covariateWithoutIov = .filterCovariatesByClass( prop( model, "modelCovariates" ), "CategoricalCovariate"        )
  covariateWithIov    = .filterCovariatesByClass( prop( model, "modelCovariates" ), "CategoricalCovariateWithIOV" )

  covWithoutIovNames = map_chr( covariateWithoutIov, ~ prop( .x, "name" ) )
  covWithIovNames    = map_chr( covariateWithIov,    ~ prop( .x, "name" ) )

  maxOccasions = getNumberOfOccasionsForModel( model )
  zeroEffect   = set_names( rep( 0, length( muValues ) ), names( muValues ) )

  # cumulated covariate effect (without IOV) for one combination row.
  effectWithoutIov = function( rowIndex ) {
    if ( length( covWithoutIovNames ) == 0L ) return( zeroEffect )
    imap( covWithoutIovNames, function( covName, j ) {
      catIdx  = as.integer( covWithoutIovGrid[ rowIndex, j ] )
      catName = prop( covariateWithoutIov[[ j ]], "categories" )[[ catIdx ]]
      pluck( covariatesEffect, "CategoricalCovariate", covName, catName )
    }) |> reduce( `+`, .init = zeroEffect )
  }

  # cumulated covariate effect (with IOV) for one combination row and occasion.
  effectWithIov = function( rowIndex, occIndex ) {
    if ( length( covWithIovNames ) == 0L ) return( zeroEffect )
    occName = paste0( "occasion_", occIndex )
    imap( covWithIovNames, function( covName, j ) {
      seqIdx    = as.integer( covWithIovGrid[ rowIndex, j ] )
      sequences = prop( covariateWithIov[[ j ]], "sequences" )
      seqNames  = names( sequences ) %||% paste0( "sequence_", seq_along( sequences ) )
      seqName   = seqNames[[ seqIdx ]]
      effList   = pluck( covariatesEffect, "CategoricalCovariateWithIOV", covName, seqName )
      if ( is.null( effList ) )
        stop(
          sprintf(
            "No IOV effect list for covariate '%s', sequence '%s'.",
            covName, seqName
          ),
          call. = FALSE
        )
      if ( !occName %in% names( effList ) )
        stop(
          sprintf(
            "No IOV effect for covariate '%s', sequence '%s', %s.",
            covName, seqName, occName
          ),
          call. = FALSE
        )
      effList[[ occName ]]
    }) |> reduce( `+`, .init = zeroEffect )
  }

  modelParamsForCov = pmap( combinations, function( name,
                                                    covariateWithoutIovCombination,
                                                    covariateWithIovCombination, ... ) {
    baseEffect = effectWithoutIov( covariateWithoutIovCombination )

    map( seq_len( maxOccasions ), function( occ ) {
      combinedEffect = baseEffect + effectWithIov( covariateWithIovCombination, occ )

      if ( is.null( modelCovariatesEquation ) ) {
        if ( !all( combinedEffect == 0 ) )
          stop( "`modelCovariatesEquation` must be defined when covariate effects are present." )
        return( muValues )
      }

      prop(
        computeCovariateValue( modelCovariatesEquation, beta = muValues, combinedEffect = combinedEffect ),
        "value"
      )
    }) |> set_names( paste0( "occ", seq_len( maxOccasions ) ) )

  }) |> set_names( combinations$name )

  prop( model, "modelParametersWithCovariates" ) = modelParamsForCov
  model
}

#' Cache the IIV variance diagonal (\eqn{\omega^2}) on the model.
#' @param model First argument of generic.
#' @return \code{Model} with \code{omegaWithIOV} variance matrix populated.
#' @name evaluateOmegaMatrixFromCovariates
#' @keywords internal
method( evaluateOmegaMatrixFromCovariates, Model ) = function( model ) {

  modelParameters = prop( model, "modelParameters" )
  paramNames      = map_chr( modelParameters, ~ prop( .x, "name" ) )

  # omega stores standard deviations; cache squared IIV diagonal for population FIM.
  # IOV gamma remains on parameters (occasion structure is in V, not here).
  omegaDiag = map_dbl( modelParameters, ~ prop( prop( .x, "distribution" ), "omega" ) )^2

  omegaWithIOV = diag( omegaDiag, nrow = length( paramNames ) )
  dimnames( omegaWithIOV ) = list( paramNames, paramNames )

  prop( model, "omegaWithIOV" ) = omegaWithIOV
  model
}

#' Cached IIV variance diagonal (\eqn{\omega^2}) from model parameters.
#' @param model \code{Model} object, optionally with \code{omegaWithIOV} populated.
#' @return Numeric vector of IIV variances aligned with \code{modelParameters}.
#' @noRd
#' @keywords internal
.modelOmegaIIVVariance = function( model ) {
  omegaMat = prop( model, "omegaWithIOV" )
  npar = length( prop( model, "modelParameters" ) )
  if ( is.matrix( omegaMat ) && nrow( omegaMat ) == npar && ncol( omegaMat ) == npar )
    return( diag( omegaMat ) )
  map_dbl(
    prop( model, "modelParameters" ),
    ~ prop( prop( .x, "distribution" ), "omega" )^2
  )
}

#' Precompute the finite-difference shift grid for gradient evaluation.
#'
#' Relative steps are \eqn{\varepsilon^{1/3}\max(|\mu|, 10^{-4})}. Scale mus to
#' O(1) so the absolute floor does not dominate the truncation error.
#' @param model First argument of generic.
#' @return \code{Model} with finite-difference gradient precomputation fields.
#' @name finiteDifferenceHessian
#' @keywords internal
method( finiteDifferenceHessian, Model ) = function( model ) {

  parameters = prop( model, "modelParameters" )
  pars       = map_dbl( parameters, ~ prop( prop( .x, "distribution" ), "mu" ) )
  freeIdx    = .pfimFdFreeIndices( parameters )

  prop( model, "parametersForComputingGradient" ) =
    .pfimFiniteDifferenceGrid(
      pars, freeIdx,
      odeSolverParameters = prop( model, "odeSolverParameters" ),
      scaleToOdeTol = .pfimIsOdeModel( model )
    )
  model
}

#' Update parameter mus for a specific occasion.
#' @param baseModelParameters List of \code{ModelParameter} objects to copy and update.
#' @param occasionParams Named numeric vector of occasion-specific \code{mu} values.
#' @param respectFixedMu Logical scalar; keep fixed parameters unchanged when \code{TRUE}.
#' @return Updated list of \code{ModelParameter} objects.
#' @noRd
#' @keywords internal
.updateParamsForOccasion = function( baseModelParameters, occasionParams, respectFixedMu = FALSE ) {
  map( baseModelParameters, function( param ) {
    param = .pfimCloneS7( param )
    if ( respectFixedMu && prop( param, "fixedMu" ) ) return( param )
    pName = prop( param, "name" )
    newMu = unname( occasionParams[[ pName ]] )
    if ( is.null( newMu ) || length( newMu ) != 1L ) return( param )
    distribution = prop( param, "distribution" )
    prop( distribution, "mu" ) = newMu
    prop( param, "distribution" ) = distribution
    param
  })
}

#' Build a temporary model for one covariate occasion.
#' @param model \code{Model} object used as template.
#' @param arm \code{Arm} object used to define administration.
#' @param baseModelParameters List of base \code{ModelParameter} objects.
#' @param occasionParams Named numeric vector of occasion-specific parameters.
#' @param needFd Logical scalar indicating whether FD setup is required.
#' @return \code{Model} configured for one covariate-occasion evaluation.
#' @noRd
#' @keywords internal
.modelForCovariateOccasion = function(
    model, arm, baseModelParameters, occasionParams, needFd = FALSE ) {
  tempModel = .pfimCloneS7( model )
  # needFd: overwrite every mu with occasion theta so the FD grid includes
  # fixedMu params at their covariate-adjusted value. Plain evaluation keeps
  # fixedMu at the base (reference) mu via respectFixedMu = TRUE.
  prop( tempModel, "modelParameters" ) = .updateParamsForOccasion(
    baseModelParameters, occasionParams, respectFixedMu = !needFd
  )
  if ( needFd ) {
    # FD is on theta for this occasion only; strip covxocc so .fdModelEvaluations
    # does not recurse into evaluateModelGradientWithCovariates.
    prop( tempModel, "modelCovariates" ) = list()
    prop( tempModel, "numberOfOccasions" ) = 1L
    tempModel = .pfimFiniteDifferenceHessianCached( tempModel )
  }
  .pfimPrepareModelForEvaluation( tempModel, arm )
}
