NULL

evaluateModelGradientCore = function( model, arm ) {

  parameters  = prop( model, "modelParameters" )
  paramNames  = map_chr( parameters, ~ prop( .x, "name" ) )
  outputNames = prop( model, "outputNames" )
  gradParams  = prop( model, "parametersForComputingGradient" )

  XcolsInv = gradParams$XcolsInv
  shifted  = gradParams$shifted
  frac     = gradParams$frac

  # Strip covariates / force single occasion so the core stays atomic.
  tempModel = model
  prop( tempModel, "modelCovariates" ) = list()
  prop( tempModel, "numberOfOccasions" ) = 1L

  # Each column of `shifted` is one finite-difference perturbation → one ODE solve.
  evaluations = map(
    seq_len( ncol( shifted ) ),
    function( iter ) {
      shiftedParams = .pfimShiftModelParameters( parameters, shifted[ , iter ] )
      iterModel = tempModel
      prop( iterModel, "modelParameters" ) = shiftedParams
      iterModel = .pfimPrepareModelForEvaluation( iterModel, arm )
      evaluateModel( iterModel, arm )
    }
  )

  map( outputNames, function( outName ) {
    outputMatrix = evaluations |>
      map( ~ .x[[ outName ]][ , 2L ] ) |>
      as.data.frame() |>
      t()

    raw   = XcolsInv %*% outputMatrix / frac
    grads = as.data.frame(
      matrix( t( raw )[ , 2L:( 1L + length( parameters ) ) ], ncol = length( parameters ) )
    )
    colnames( grads ) = paramNames
    grads
  } ) |> set_names( outputNames )
}

method( evaluateModelGradient, Model ) = function( model, arm ) {
  if ( !usesCovariateOccasionStructure( model ) ) {
    evaluateModelGradientCore( model, arm )
  } else {
    evaluateModelGradientWithCovariates( model, arm, evaluateModelGradientCore )
  }
}

.imap_dfc = function( .x, .f ) {
  .x = as.data.frame( .x, check.names = FALSE, stringsAsFactors = FALSE )
  nm = names( .x )
  if ( is.null( nm ) ) nm = paste0( "V", seq_len( ncol( .x ) ) )
  cols = Map( .f, .x, nm )
  as.data.frame( cols, check.names = FALSE, stringsAsFactors = FALSE )
}

evaluateModelGradientWithCovariates = function( model, arm, evaluateModelGradientCore ) {

  covariatesCombinations = prop( model, "covariatesCombination" )$combinations
  modelParamsWithCov     = prop( model, "modelParametersWithCovariates" )
  covariatesEffect       = prop( model, "covariatesEffect" )
  modelCovariatesEq      = prop( model, "modelCovariatesEquation" )
  modelParameters        = prop( model, "modelParameters" )
  modelCovariates        = prop( model, "modelCovariates" )

  paramNames = map_chr( modelParameters, ~ prop( .x, "name" ) )
  muValues   = map_dbl( modelParameters, ~ prop( prop( .x, "distribution" ), "mu" ) ) |>
    set_names( paramNames )

  covEqType = if ( S7::S7_inherits( modelCovariatesEq, Additive ) ) modelCovEqAdditive
  else                                              modelCovEqExponential

  covariateWithoutIov = .filterCovariatesByClass( modelCovariates, "CategoricalCovariate"        )
  covariateWithIov    = .filterCovariatesByClass( modelCovariates, "CategoricalCovariateWithIOV" )

  # Built once, reused across all combination x occasion iterations.
  betaList = .buildBetaList( covariateWithoutIov, covariateWithIov, covariatesEffect )

  # predicates and correction helpers (closed over local env)

  .isBetaActive = function( betaInfo, combinationInfo, occasionIndex ) {
    if ( !betaInfo$isIOV ) {
      identical( combinationInfo[[ betaInfo$covName ]], betaInfo$category )
    } else {
      seqLabel = combinationInfo[[ betaInfo$covName ]]
      if ( is.null( seqLabel ) ) return( FALSE )
      cov      = covariateWithIov[[ betaInfo$covIndex ]]
      seqs     = prop( cov, "sequences" )
      seqNames = if ( is.null( names( seqs ) ) ) paste0( "sequence_", seq_along( seqs ) ) else names( seqs )
      seqIdx   = match( seqLabel, seqNames )
      if ( is.na( seqIdx ) ) return( FALSE )
      seqValues = seqs[[ seqIdx ]]
      if ( occasionIndex > length( seqValues ) ) return( FALSE )
      identical( seqValues[[ occasionIndex ]], betaInfo$category )
    }
  }

  # Multiplicative correction per parameter for the additive model:
  # d(theta)/d(mu) = (1 + beta * cov).
  .additiveMuCorrection = function( combinationInfo, occIndex ) {
    correction = set_names( rep( 1, length( muValues ) ), names( muValues ) )

    applyCorrection = function( covList, effectClass ) {
      reduce( covList, function( corr, cov ) {
        covName = prop( cov, "name" )
        val     = combinationInfo[[ covName ]]
        if ( is.null( val ) ) return( corr )

        activeCat = if ( S7::S7_inherits( cov, CategoricalCovariateWithIOV ) ) {
          seqs     = prop( cov, "sequences" )
          seqNames = names( seqs ) %||% paste0( "sequence_", seq_along( seqs ) )
          seqIdx   = match( val, seqNames )
          if ( is.na( seqIdx ) || occIndex > length( seqs[[ seqIdx ]] ) ) return( corr )
          seqs[[ seqIdx ]][[ occIndex ]]
        } else val

        if ( identical( activeCat, prop( cov, "categories" )[[1L]] ) ) return( corr )

        effects = covariatesEffect[[ effectClass ]][[ covName ]]
        catIdx  = match( activeCat, prop( cov, "categories" ) )
        if ( is.na( catIdx ) || is.null( effects[[ catIdx ]] ) ) return( corr )

        eff    = effects[[ catIdx ]]
        active = names( eff )[ eff != 0 ]
        reduce( active,
                function( cr, pName ) { cr[[ pName ]] = cr[[ pName ]] * ( 1 + eff[[ pName ]] ); cr },
                .init = corr )
      }, .init = correction )
    }

    correction = applyCorrection( covariateWithoutIov, "CategoricalCovariate"        )
    correction = applyCorrection( covariateWithIov,    "CategoricalCovariateWithIOV" )
    correction
  }

  # Main loop

  map( covariatesCombinations$name, function( combinationName ) {
    paramsForCombo  = modelParamsWithCov[[ combinationName ]]
    combinationInfo = parseCombinationName( combinationName, modelCovariates )
    proportion      = covariatesCombinations$proportion[[ match( combinationName,
                                                                 covariatesCombinations$name ) ]]

    evaluationByOccasion = imap( paramsForCombo, function( occasionParams, occasion ) {
      occIndex = as.integer( str_remove( occasion, "^occ" ) )

      updatedParams = map2( modelParameters, occasionParams, function( param, newMu ) {
        distr = prop( param, "distribution" )
        prop( distr, "mu" ) = newMu
        prop( param, "distribution" ) = distr
        param
      })

      tempModel = model
      prop( tempModel, "modelParameters" ) = updatedParams
      prop( tempModel, "modelCovariates" ) = list()
      prop( tempModel, "numberOfOccasions" ) = 1L
      tempModel = .pfimFiniteDifferenceHessianCached( tempModel )
      tempModel = .pfimPrepareModelForEvaluation( tempModel, arm )

      gradTheta = evaluateModelGradientCore( tempModel, arm )

      # Gradient w.r.t. mu
      gradMu = if ( covEqType == modelCovEqExponential ) {
        # d(f)/d(mu) = d(f)/d(theta) x (theta / mu)
        map( gradTheta, function( gradDf ) {
          .imap_dfc( gradDf, function( gradVec, pName ) {
            mu_i    = muValues[[ pName ]]
            theta_i = occasionParams[[ pName ]]
            if ( is.null( mu_i ) || mu_i == 0 || is.na( mu_i ) ) gradVec
            else gradVec * ( theta_i / mu_i )
          })
        })
      } else {
        # d(f)/d(mu) = d(f)/d(theta) x (1 + beta x cov)
        correction = .additiveMuCorrection( combinationInfo, occIndex )
        map( gradTheta, function( gradDf ) {
          .imap_dfc( gradDf, function( gradVec, pName ) gradVec * correction[[ pName ]] )
        })
      }

      # Gradient w.r.t. beta
      outputNames = prop( model, "outputNames" )

      gradWithCovariates = map( outputNames, function( outName ) {
        gradMu_out    = gradMu[[ outName ]]
        gradTheta_out = gradTheta[[ outName ]]
        nRows         = nrow( gradMu_out )

        muCols = map( paramNames, ~ gradMu_out[[ .x ]] ) |>
          set_names( paste0( "mu_", paramNames ) )

        betaCols = imap( betaList, function( betaInfo, betaName ) {
          if ( !.isBetaActive( betaInfo, combinationInfo, occIndex ) )
            return( rep( 0, nRows ) )
          pName = betaInfo$param
          if ( covEqType == modelCovEqExponential )
            gradTheta_out[[ pName ]] * occasionParams[[ pName ]]   # d(f)/d(beta) = d(f)/d(theta) x theta
          else
            gradTheta_out[[ pName ]] * muValues[[ pName ]]         # d(f)/d(beta) = d(f)/d(theta) x mu
        })

        as.data.frame( c( muCols, betaCols ) )
      }) |> set_names( outputNames )

      list( occasion = occasion, gradient = gradWithCovariates )
    })

    list( combination = combinationName, proportion = proportion, gradients = evaluationByOccasion )
  })
}

# One entry per (covariate, non-reference category, affected parameter) triple.
.buildBetaList = function( covariateWithoutIov, covariateWithIov, covariatesEffect ) {

  betasWithoutIov = if ( length( covariateWithoutIov ) > 0L ) {
    covariateWithoutIov |>
      map( function( cov ) {
        covName    = prop( cov, "name"       )
        categories = prop( cov, "categories" )
        effects    = covariatesEffect$CategoricalCovariate[[ covName ]]

        if ( length( categories ) <= 1L ) return( list() )

        ( seq_len( length( categories ) - 1L ) + 1L ) |>
          map( function( icat ) {
            affectedParams = names( effects[[ icat ]] )[ effects[[ icat ]] != 0 ]
            map( affectedParams, ~ list(
              name  = paste0( "beta_", .x, "_", covName, "_", categories[[ icat ]] ),
              value = list( covName  = covName,
                            category = categories[[ icat ]],
                            param    = .x,
                            isIOV    = FALSE )
            ))
          }) |> list_flatten()
      }) |> list_flatten() |>
      (\( lst ) set_names( lst, map_chr( lst, function( x ) x[["name"]] ) ))() |>
      map( "value" )
  } else list()

  betasWithIov = if ( length( covariateWithIov ) > 0L ) {
    imap( covariateWithIov, function( cov, covIndex ) {
      covName    = prop( cov, "name"       )
      sequences  = prop( cov, "sequences"  )
      categories = prop( cov, "categories" )
      effects    = covariatesEffect$CategoricalCovariateWithIOV[[ covName ]]

      if ( length( categories ) <= 1L ) return( list() )

      ( seq_len( length( categories ) - 1L ) + 1L ) |>
        map( function( icat ) {
          category = categories[[ icat ]]

          affectedParams = seq_along( sequences ) |>
            map( function( iseq ) {
              effSeq = effects[[ iseq ]]
              if ( is.null( effSeq ) ) return( character( 0L ) )
              seqVals     = sequences[[ iseq ]]
              matchingOcc = detect_index( seq_along( seqVals ), ~ seqVals[[ .x ]] == category )
              if ( matchingOcc == 0L ) return( character( 0L ) )
              occEff = effSeq[[ matchingOcc ]]
              names( occEff )[ occEff != 0 ]
            }) |> list_c() |> unique()

          map( affectedParams, ~ list(
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
      (\( lst ) set_names( lst, map_chr( lst, function( x ) x[["name"]] ) ))() |>
      map( "value" )
  } else list()

  c( betasWithoutIov, betasWithIov )
}
