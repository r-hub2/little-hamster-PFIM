# Covariate/occasion gradient assembly and finite-difference Hessian for population FIM.

# Covariate combination and parameter methods for Model.
NULL

# S7 Method Implementations

method( evaluateModelWithCovariates, Model ) = function( model, arm, evaluateModelCore ) {

  covariatesCombinations = prop( model, "covariatesCombination" )$combinations
  modelParamsWithCov     = prop( model, "modelParametersWithCovariates" )
  baseModelParameters    = prop( model, "modelParameters" )

  map( seq_len( nrow( covariatesCombinations ) ), function( i ) {
    combinationName    = covariatesCombinations$name[[i]]
    proportion         = covariatesCombinations$proportion[[i]]
    parametersForCombo = modelParamsWithCov[[ combinationName ]]

    evaluationByOccasion = map( names( parametersForCombo ), function( occasion ) {
      occasionParams = parametersForCombo[[ occasion ]]

      updatedParameters = map2( baseModelParameters, occasionParams, function( param, newMu ) {
        if ( prop( param, "fixedMu" ) ) return( param )
        distribution = prop( param, "distribution" )
        prop( distribution, "mu" ) = newMu
        prop( param, "distribution" ) = distribution
        param
      })

      tempModel = model
      prop( tempModel, "modelParameters" ) = updatedParameters
      tempModel = .pfimPrepareModelForEvaluation( tempModel, arm )
      list( occasion = occasion, evaluation = evaluateModelCore( tempModel, arm ) )
    })

    list( combination = combinationName, proportion = proportion, evaluations = evaluationByOccasion )
  })
}

method( hasCovariates, Model ) = function( model ) {
  length( prop( model, "modelCovariates" ) ) > 0L
}

method( defineCovariatesData, Model ) = function( model ) {
  model |>
    evaluateCovariatesEffects()         |>
    generateCovariatesCombination()     |>
    modelParametersWithCovariates()     |>
    evaluateOmegaMatrixFromCovariates()
}

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
      catIdx = as.integer( covWithoutIovGrid[ rowIndex, j ] )
      covariatesEffect$CategoricalCovariate[[ covName ]][[ catIdx ]]
    }) |> reduce( `+`, .init = zeroEffect )
  }

  # cumulated covariate effect (with IOV) for one combination row and occasion.
  effectWithIov = function( rowIndex, occIndex ) {
    if ( length( covWithIovNames ) == 0L ) return( zeroEffect )
    imap( covWithIovNames, function( covName, j ) {
      seqIdx    = as.integer( covWithIovGrid[ rowIndex, j ] )
      sequences = prop( covariateWithIov[[ j ]], "sequences" )
      seqName   = if ( is.null( names( sequences ) ) ) paste0( "sequence_", seqIdx ) else names( sequences )[[ seqIdx ]]
      effList   = covariatesEffect$CategoricalCovariateWithIOV[[ covName ]][[ seqName ]]
      if ( is.null( effList ) || occIndex > length( effList ) ) zeroEffect else effList[[ occIndex ]]
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

method( evaluateOmegaMatrixFromCovariates, Model ) = function( model ) {

  modelParameters = prop( model, "modelParameters" )
  paramNames      = map_chr( modelParameters, ~ prop( .x, "name" ) )

  # omega / gamma store standard deviations; square them for the variance matrix.
  omegaDiag = map_dbl( modelParameters, ~ prop( prop( .x, "distribution" ), "omega" ) )
  gammaDiag = map_dbl( modelParameters, ~ prop( .x, "gamma" ) )

  omegaMat = diag( omegaDiag^2 )
  gammaMat = diag( gammaDiag^2 )

  maxOccasions = getNumberOfOccasionsForModel( model )
  nIov         = max( maxOccasions - 1L, 0L )
  nIiv         = length( paramNames )
  nTotal       = nIiv + nIov * nIiv

  # IOV block indices: occasion k occupies rows/cols [nIiv + (k-1)*nIiv + 1 : nIiv + k*nIiv].
  # reduce() threads gammaMat into each block; with an empty list it returns .init unchanged.
  iovIndices   = map( seq_len( nIov ), ~ nIiv + ( .x - 1L ) * nIiv + seq_len( nIiv ) )
  base         = matrix( 0, nrow = nTotal, ncol = nTotal )
  base[ seq_len( nIiv ), seq_len( nIiv ) ] = omegaMat

  omegaWithIOV = reduce( iovIndices,
                         function( mat, idx ) { mat[ idx, idx ] = gammaMat; mat },
                         .init = base
  )

  iovNames = if ( nIov > 0L ) {
    paste0( rep( paramNames, nIov ), "_occ", rep( seq_len( nIov ), each = nIiv ) )
  } else character( 0L )

  dimnames( omegaWithIOV ) = list( c( paramNames, iovNames ), c( paramNames, iovNames ) )

  prop( model, "omegaWithIOV" ) = omegaWithIOV
  model
}

method( finiteDifferenceHessian, Model ) = function( model ) {

  pars    = map_dbl( prop( model, "modelParameters" ), ~ prop( prop( .x, "distribution" ), "mu" ) )
  npar    = length( pars )
  relStep = .Machine$double.eps^( 1 / 3 )
  incr    = pmax( abs( pars ), 1e-4 ) * relStep
  baseInd = diag( npar )

  # Build shift columns and quadratic-approximation fraction vector.
  extraCols = map( seq_len( npar - 1L ), ~ baseInd[ , .x ] + baseInd[ , -seq_len( .x ) ] )
  extraFrac = map( seq_len( npar - 1L ), ~ incr[ .x ] * incr[ -seq_len( .x ) ]            )

  cols    = c( list( 0, baseInd, -baseInd ), extraCols )
  frac    = c( 1, incr, incr^2, unlist( extraFrac ) )
  indMat  = do.call( cbind, cols )
  shifted = pars + incr * indMat

  indMatT = t( indMat )
  Xcols   = c(
    list( 1, indMatT, indMatT^2 ),
    map( seq_len( npar - 1L ), ~ indMatT[ , .x ] * indMatT[ , -seq_len( .x ) ] )
  )

  prop( model, "parametersForComputingGradient" ) = list(
    XcolsInv = .safeSolve( do.call( cbind, Xcols ) ),
    shifted  = shifted,
    frac     = frac
  )
  model
}
