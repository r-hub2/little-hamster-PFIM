#' @title Optimization
#' @description
#' Design optimization: search over sampling times (or related design variables)
#' using a metaheuristic or exchange algorithm. Holds a \code{PFIMProject} in the
#' \code{project} slot (composition) plus optimization results.
#' @inheritParams PFIMProject
#' @param project Nested \code{PFIMProject} (filled by the constructor).
#' @param optimisationDesign List with initial and optimal design evaluations.
#' @param optimisationAlgorithmOutputs Raw outputs from the optimization algorithm.
#' @include PFIMProject.R
#' @include pfim-project-access.R
#' @export

Optimization = new_class( "Optimization",
                          package    = "PFIM",
                          properties = list(
                            project                      = new_property( PFIMProject ),
                            optimisationDesign           = new_property( class_list, default = list() ),
                            optimisationAlgorithmOutputs = new_property( class_list, default = list() )
                          ),
                          constructor = function(
                              name = character(0),
                              modelEquations = list(),
                              modelCovariatesEquation = character(0),
                              modelFromLibrary = list(),
                              modelParameters = list(),
                              modelCovariates = list(),
                              modelError = list(),
                              optimizer = character(0),
                              optimizerParameters = list(),
                              outputs = list(),
                              designs = list(),
                              fimType = character(0),
                              fim = NULL,
                              odeSolverParameters = list(),
                              project = NULL,
                              optimisationDesign = list(),
                              optimisationAlgorithmOutputs = list() ) {
                            if ( is.null( project ) ) {
                              if ( is.null( fim ) )
                                fim = .placeholderFim( fimType )
                              project = PFIMProject(
                                name = name,
                                modelEquations = modelEquations,
                                modelCovariatesEquation = modelCovariatesEquation,
                                modelFromLibrary = modelFromLibrary,
                                modelParameters = modelParameters,
                                modelCovariates = modelCovariates,
                                modelError = modelError,
                                optimizer = optimizer,
                                optimizerParameters = optimizerParameters,
                                outputs = outputs,
                                designs = designs,
                                fimType = fimType,
                                fim = fim,
                                odeSolverParameters = odeSolverParameters
                              )
                            }
                            new_object(
                              .parent                      = S7_object(),
                              project                      = project,
                              optimisationDesign           = optimisationDesign,
                              optimisationAlgorithmOutputs = optimisationAlgorithmOutputs
                            )
                          }
)
S4_register( Optimization )

method( defineFim, Optimization ) = function( pfimproject ) {
  defineFim( projectOf( pfimproject ) )
}

method( defineModelEquationsFromLibraryOfModel, Optimization ) = function( pfimproject ) {
  res = defineModelEquationsFromLibraryOfModel( projectOf( pfimproject ) )
  projectProp( pfimproject, "modelEquations" ) <- res
  res
}

method( defineModelType, Optimization ) = function( pfimproject ) {
  defineModelType( projectOf( pfimproject ) )
}

defineOptimizationAlgorithm  = new_generic( "defineOptimizationAlgorithm", c( "optimization" ) )
generateFimsFromConstraints  = new_generic( "generateFimsFromConstraints",  c( "optimization" ) )
plotWeights                  = new_generic( "plotWeights",                  c( "optimization" ) )
plotFrequencies              = new_generic( "plotFrequencies",              c( "optimization" ) )
optimizeDesign               = new_generic( "optimizeDesign",  c( "optimizationObject", "optimizationAlgorithm" ) )
constraintsTableForReport    = new_generic( "constraintsTableForReport",    c( "optimizationAlgorithm" ) )

.getOptimalEval = function( optimization ) {
  prop( optimization, "optimisationDesign" )$evaluationOptimalDesign
}

.getOptimalFim = function( optimization ) {
  eval  = .getOptimalEval( optimization )
  fim   = prop( eval, "fim" )
  setEvaluationFim( fim, eval )
}

#' Instantiate the optimizer from project settings
#' @name defineOptimizationAlgorithm
#' @export
method( defineOptimizationAlgorithm, Optimization ) = function( optimization ) {
  optimizerParameters = projectProp( optimization, "optimizer" )

  switch( optimizerParameters,
          MultiplicativeAlgorithm = MultiplicativeAlgorithm(),
          FedorovWynnAlgorithm    = FedorovWynnAlgorithm(),
          PSOAlgorithm            = PSOAlgorithm(),
          PGBOAlgorithm           = PGBOAlgorithm(),
          SimplexAlgorithm        = SimplexAlgorithm(),
          stop( sprintf( "Unknown optimizer: '%s'", optimizerParameters ) )
  )
}

#' Run design optimization
#' @name run
#' @export
method( run, Optimization ) = function( pfimproject ) {
  .invalidateEvalModelCache( pfimproject )
  optimizationAlgorithm = defineOptimizationAlgorithm( pfimproject )
  projectProp( pfimproject, "fim" ) <- defineFim( pfimproject )

  if ( length( projectProp( pfimproject, "modelFromLibrary" ) ) != 0L )
    projectProp( pfimproject, "modelEquations" ) <-
      defineModelEquationsFromLibraryOfModel( pfimproject )

  optimizeDesign( pfimproject, optimizationAlgorithm )
}

method( show, Optimization ) = function( object ) {
  optimisationDesign      = prop( object, "optimisationDesign" )
  evaluationInitialDesign = optimisationDesign$evaluationInitialDesign
  evaluationOptimalDesign = optimisationDesign$evaluationOptimalDesign
  .armTable = function( evaluation ) {
    designs  = prop( evaluation, "designs" )
    armsData = list_flatten( map( pluck( map( designs, ~ prop( .x, "arms" ) ), 1L ),
                                  getArmData ) )
    df = map( armsData, ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |> list_rbind()
    colnames( df ) = c( "Arms name", "Number of subjects", "Outcome", "Dose", "Sampling times" )
    df
  }

  fimInitialDesign = setEvaluationFim( prop( evaluationInitialDesign, "fim" ), evaluationInitialDesign )
  fimOptimalDesign = setEvaluationFim( prop( evaluationOptimalDesign, "fim" ), evaluationOptimalDesign )

  cat( "\n--- Initial design ---\n\n" )
  print( .armTable( evaluationInitialDesign ) )
  showFIM( fimInitialDesign )
  cat( "\n--- Optimal design ---\n\n" )
  print( .armTable( evaluationOptimalDesign ) )
  showFIM( fimOptimalDesign )

  invisible( object )
}

#' Fisher matrix blocks for the optimal design
#' @name getFisherMatrix
#' @export
method( getFisherMatrix, Optimization ) = function( pfimproject ) {
  fim = .getOptimalFim( pfimproject )
  list(
    fisherMatrix    = prop( fim, "fisherMatrix"    ),
    fixedEffects    = prop( fim, "fixedEffects"    ),
    varianceEffects = prop( fim, "varianceEffects" )
  )
}

#' Standard errors from the optimal design FIM
#' @name getSE
#' @export
method( getSE, Optimization ) = function( pfimproject ) {
  prop( .getOptimalFim( pfimproject ), "SEAndRSE" )$SE
}

#' Relative standard errors from the optimal design FIM
#' @name getRSE
#' @export
method( getRSE, Optimization ) = function( pfimproject ) {
  prop( .getOptimalFim( pfimproject ), "SEAndRSE" )$RSE
}

#' Parameter shrinkage from the Bayesian FIM
#' @name getShrinkage
#' @export
method( getShrinkage, Optimization ) = function( pfimproject ) {
  prop( .getOptimalFim( pfimproject ), "shrinkage" )
}

#' Determinant of the optimal design FIM
#' @name getDeterminant
#' @export
method( getDeterminant, Optimization ) = function( pfimproject ) {
  det( prop( .getOptimalFim( pfimproject ), "fisherMatrix" ) )
}

#' Parameter correlations from the optimal design FIM
#' @name getCorrelationMatrix
#' @export
method( getCorrelationMatrix, Optimization ) = function( pfimproject ) {
  cor( prop( .getOptimalFim( pfimproject ), "fisherMatrix" ) )
}

#' D-criterion of the optimal design
#' @name getDcriterion
#' @export
method( getDcriterion, Optimization ) = function( pfimproject ) {
  Dcriterion( .getOptimalFim( pfimproject ) )
}

#' Sensitivity indices for the optimal design
#' @name plotSensitivityIndices
#' @export
method( plotSensitivityIndices, Optimization ) = function( pfimproject ) {
  plotSensitivityIndices( .getOptimalEval( pfimproject ) )
}

#' SE barplot for the optimal design
#' @name plotSE
#' @export
method( plotSE, Optimization ) = function( pfimproject ) {
  plotSEFIM( .getOptimalFim( pfimproject ), .getOptimalEval( pfimproject ) )
}

#' RSE barplot for the optimal design
#' @name plotRSE
#' @export
method( plotRSE, Optimization ) = function( pfimproject ) {
  plotRSEFIM( .getOptimalFim( pfimproject ), .getOptimalEval( pfimproject ) )
}

#' Multiplicative algorithm weights by iteration
#' @name plotWeights
#' @export
method( plotWeights, Optimization ) = function( optimization ) {
  plotWeightsMultiplicativeAlgorithm(
    optimization,
    prop( optimization, "optimisationAlgorithmOutputs" )$optimizationAlgorithm
  )
}

#' Fedorov-Wynn optimal frequencies
#' @name plotFrequencies
#' @export
method( plotFrequencies, Optimization ) = function( optimization ) {
  plotFrequenciesFedorovWynnAlgorithm(
    optimization,
    prop( optimization, "optimisationAlgorithmOutputs" )$optimizationAlgorithm
  )
}

#' Stratified subsample of constraint-grid cells
#' @name .pfimConstraintTaskIndices
#' @keywords internal
.pfimConstraintTaskIndices = function( totalIterations, numberOfDoses, nCombinations ) {
  max_tasks = pfim_get_option( "constraints.maxTasks", NULL )
  if ( is.null( max_tasks ) || !is.finite( max_tasks ) || max_tasks >= totalIterations )
    return( seq_len( totalIterations ) )
  max_tasks = max( 1L, as.integer( max_tasks ) )
  base      = max_tasks %/% numberOfDoses
  remainder = max_tasks %% numberOfDoses
  indices   = integer( 0 )
  for ( dose_idx in seq_len( numberOfDoses ) ) {
    pool   = ( ( dose_idx - 1L ) * nCombinations + 1L ):( dose_idx * nCombinations )
    n_pick = base + as.integer( dose_idx <= remainder )
    n_pick = min( length( pool ), n_pick )
    if ( n_pick > 0L ) indices = c( indices, sample( pool, n_pick ) )
  }
  sort( unique( indices ) )
}

#' Enumerate FIMs on the dose and sampling grid
#' @name generateFimsFromConstraints
#' @export
method( generateFimsFromConstraints, Optimization ) = function( optimization ) {
  show_progress = isTRUE( projectProp( optimization, "optimizerParameters" )$showProcess ) ||
    isTRUE( pfim_get_option( "verbose", FALSE ) )
  .pfimFimCacheBegin( optimization )
  pfim_set_option( fim.cache.hits = 0L )
  old_batch = pfim_get_option( "eval.batch", FALSE )
  on.exit( pfim_set_option( eval.batch = old_batch ), add = TRUE )
  pfim_set_option( eval.batch = TRUE )
  evaluation  = .evaluationFromProject( optimization )
  baseModel   = rebuildEvalModel( evaluation, finiteDifference = TRUE )
  baseFim     = defineFim( evaluation )
  designs     = projectProp( optimization, "designs" )
  designNames = map_chr( designs, ~ prop( .x, "name" ) )
  dosesForFIMs = map( designs, ~ generateDosesCombination( .x ) ) |> set_names( designNames )
  samplingsForFIMs = map( designs, ~ generateSamplingTimesCombination( .x ) ) |>
    set_names( designNames )
  allResults = imap( set_names( designs, designNames ), function( design, designName ) {
    arms = prop( design, "arms" )
    if ( length( arms ) > 1L )
      stop(
        "generateFimsFromConstraints: single-arm designs only; design '",
        designName, "' has ", length( arms ), " arm(s).", call. = FALSE
      )
    dosesForDesign = dosesForFIMs[[ designName ]]
    numberOfDoses  = dosesForDesign$numberOfDoses
    combinationGrid = expand.grid( map( samplingsForFIMs[[ designName ]], seq_along ) )
    nCombinations   = nrow( combinationGrid )
    totalIterations = numberOfDoses * nCombinations
    taskIndices     = .pfimConstraintTaskIndices( totalIterations, numberOfDoses, nCombinations )
    nTasks          = length( taskIndices )
    if ( nTasks < totalIterations && show_progress )
      message( sprintf(
        "FIM evaluation: sampling %d / %d constraint-grid cells (PFIM.constraints.maxTasks)",
        nTasks, totalIterations
      ) )
    eval_cell = get(
      ".evaluateFimConstraintsCell",
      envir = asNamespace( "PFIM" ),
      inherits = FALSE
    )
    perDesignResults = map( taskIndices, function( fimIndex ) {
      iterDose = ( fimIndex - 1L ) %/% nCombinations + 1L
      iterComb = ( fimIndex - 1L ) %% nCombinations + 1L
      eval_cell(
        fimIndex         = match( fimIndex, taskIndices ),
        iterDose         = iterDose,
        iterComb         = iterComb,
        totalIterations  = nTasks,
        show_progress    = show_progress,
        evaluation       = evaluation,
        design           = design,
        arms             = arms,
        dosesForDesign   = dosesForDesign,
        samplingsForFIMs = samplingsForFIMs,
        designName       = designName,
        combinationGrid  = combinationGrid,
        baseModel        = baseModel,
        baseFim          = baseFim
      )
    } )
    for ( cell in perDesignResults )
      if ( !is.null( cell$cachedEvaluation ) )
        .pfimFimCacheRegister( cell$cachedEvaluation )
    perDesignResults
  } )
  list(
    listArms                    = map( allResults, ~ map( .x, "armResult" ) ),
    dimFim                      = allResults[[ 1L ]][[ 1L ]]$dimFim,
    listFimsAlgoFW              = map( allResults, ~ map( .x, "fisherMatrixForAlgoFW" ) ),
    listFimsAlgoMult            = map( allResults, ~ map( .x, "fisherMatrix" ) ),
    samplingsForFedorovWynnAlgo = map( allResults, ~ map( .x, "samplingsForFW" ) )
  )
}

#' HTML report for an optimization run
#' @name Report
#' @export
method( Report, Optimization ) = function( pfimproject, outputPath, outputFile, plotOptions ) {
  projectName = projectProp( pfimproject, "name" )
  optimisationDesign      = prop( pfimproject, "optimisationDesign" )
  evaluationInitialDesign = optimisationDesign$evaluationInitialDesign
  evaluationOptimalDesign = optimisationDesign$evaluationOptimalDesign
  evaluationOutputs       = projectProp( pfimproject, "outputs" )
  model          = rebuildEvalModel( evaluationInitialDesign, finiteDifference = FALSE )
  modelEquations = prop( model, "modelEquations" )
  modelErrorData = prop( evaluationInitialDesign, "modelError" ) |>
    map( getModelErrorData ) |>
    map( ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |>
    list_rbind()
  colnames( modelErrorData ) = c( "Output", "Type", "$\\sigma_{slope}$", "$\\sigma_{inter}$" )
  modelErrorTable = kbl( modelErrorData, align = "c" ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )
  modelParametersTable = .buildModelParametersKable( prop( evaluationInitialDesign, "modelParameters" ) )
  modelCovariates = prop( evaluationInitialDesign, "modelCovariates" )
  hasCov          = length( modelCovariates ) > 0L
  covariatesTable     = .buildCovariatesKable( modelCovariates )
  covariateTestTables = if ( hasCov ) .buildCovariateTestSection( evaluationOptimalDesign ) else NULL
  .buildArmTable = function( evaluation, colnamesVec ) {
    designs  = prop( evaluation, "designs" )
    armsData = list_flatten( map( pluck( map( designs, ~ prop( .x, "arms" ) ), 1L ), getArmData ) )
    df = map( armsData, ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |> list_rbind()
    colnames( df ) = colnamesVec
    df
  }
  administrationData = prop( evaluationInitialDesign, "designs" ) |>
    map( ~ prop( .x, "arms" ) ) |> pluck( 1L ) |>
    map( armAdministration ) |> list_flatten() |>
    map( ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |> list_rbind()
  colnames( administrationData ) = c(
    "Design name", "Arms name", "Number of subject",
    "Outcome", "Dose", "Time dose", "$\\tau$", "$T_{inf}$"
  )
  administrationTable = kbl( administrationData, align = c( "l", "l", "l", "c", "c", "c", "c" ) ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )
  armCols            = c( "Arms name", "Number of subjects", "Outcome", "Dose", "Sampling times" )
  initialDesignData  = .buildArmTable( evaluationInitialDesign, armCols )
  initialDesignTable = kbl( initialDesignData, align = c( "l", "c", "c", "c" ) ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )
  fimInitial            = setEvaluationFim( prop( evaluationInitialDesign, "fim" ), evaluationInitialDesign )
  fimInitialDesignTable = tablesForReport( fimInitial, evaluationInitialDesign )
  optimAlgoOutputs      = prop( pfimproject, "optimisationAlgorithmOutputs" )
  optimizationAlgorithm = optimAlgoOutputs$optimizationAlgorithm
  constraintsData       = constraintsTableForReport(
    optimizationAlgorithm,
    map( prop( evaluationInitialDesign, "designs" ), ~ prop( .x, "arms" ) )
  )
  optimalDesignData  = .buildArmTable( evaluationOptimalDesign, armCols )
  optimalDesignTable = kbl( optimalDesignData, align = c( "l", "c", "c", "c", "c" ) ) |>
    kable_styling( bootstrap_options = "hover", full_width = FALSE,
                   position = "center", font_size = 13 )
  fimOptimal      = setEvaluationFim( prop( evaluationOptimalDesign, "fim" ), evaluationOptimalDesign )
  fimOptimalTable = tablesForReport( fimOptimal, evaluationOptimalDesign )
  plotsEvaluationData        = plotEvaluation( evaluationOptimalDesign, plotOptions )
  plotSensitivityIndicesData = plotSensitivityIndices( evaluationOptimalDesign, plotOptions )
  plotSEData                 = plotSEFIM( fimOptimal, evaluationOptimalDesign )
  plotRSEData                = plotRSEFIM( fimOptimal, evaluationOptimalDesign )
  reportTables = list(
    evaluationOutputs         = evaluationOutputs,
    modelEquations            = modelEquations,
    modelErrorTable           = modelErrorTable,
    modelParametersTable      = modelParametersTable,
    covariatesTable           = covariatesTable,
    covariateTestTables       = covariateTestTables,
    administrationTable       = administrationTable,
    initialDesignTable        = initialDesignTable,
    constraintsTableForReport = constraintsData,
    armsConstraintsTable      = constraintsData,
    fimInitialDesignTable     = fimInitialDesignTable,
    optimalDesignTable        = optimalDesignTable,
    fimOptimalTable           = fimOptimalTable,
    plotsEvaluation           = plotsEvaluationData,
    plotSensitivityIndices    = plotSensitivityIndicesData,
    plotSE                    = plotSEData,
    plotRSE                   = plotRSEData,
    fim                       = fimOptimal,
    evaluationForPlot         = evaluationOptimalDesign,
    pfimproject               = pfimproject,
    projectName               = projectName
  )
  generateReportOptimization( fimOptimal, optimizationAlgorithm,
                              reportTables, outputFile, outputPath )
}

#' @keywords internal
.buildFlatSamplingLayout = function( design ) {
  arms           = prop( design, "arms" )
  initial_flat   = numeric()
  windows_list   = list()
  sorting_groups = list()
  group_specs    = list()
  flat_dim       = 0L
  group_idx      = 0L
  for ( arm in arms ) {
    armName = prop( arm, "name" )
    for ( st in prop( arm, "samplingTimes" ) ) {
      outcome = prop( st, "outcome" )
      samps   = prop( st, "samplings" )
      n_samps = length( samps )
      group_idx = group_idx + 1L
      initial_flat = c( initial_flat, samps )
      constraint = pluck(
        keep( prop( arm, "samplingTimesConstraints" ), ~ prop( .x, "outcome" ) == outcome ),
        1L
      )
      win        = prop( constraint, "samplingsWindows" )
      win_matrix = do.call( rbind, lapply( win, unlist ) )
      for ( i in seq_len( n_samps ) ) windows_list = append( windows_list, list( win_matrix ) )
      sorting_groups = append( sorting_groups, list( seq( flat_dim + 1L, flat_dim + n_samps ) ) )
      group_specs[[ group_idx ]] = list(
        arm = arm, armName = armName, outcome = outcome, constraint = constraint
      )
      flat_dim = flat_dim + n_samps
    }
  }
  list(
    initial_flat   = initial_flat,
    windows_list   = windows_list,
    sorting_groups = sorting_groups,
    group_specs    = group_specs
  )
}

#' @keywords internal
.sampleFlatFromConstraints = function( layout ) {
  flat = layout$initial_flat
  for ( gi in seq_along( layout$group_specs ) ) {
    spec = layout$group_specs[[ gi ]]
    idx  = layout$sorting_groups[[ gi ]]
    flat[ idx ] = sort( generateSamplingsFromSamplingConstraints( spec$constraint ) )
  }
  flat
}

#' @keywords internal
.checkFlatValidGroup = function( layout, flat_pos, group_id ) {
  spec      = layout$group_specs[[ group_id ]]
  idx       = layout$sorting_groups[[ group_id ]]
  samplings = flat_pos[ idx ]
  all( unlist( checkSamplingTimeConstraintsForMetaheuristic(
    spec$constraint, spec$arm, samplings, spec$outcome
  ) ) )
}

#' @keywords internal
.isFlatValid = function( layout, flat_pos ) {
  all( vapply(
    seq_along( layout$group_specs ),
    function( g ) .checkFlatValidGroup( layout, flat_pos, g ),
    logical( 1L )
  ) )
}

.metaheuristicFitnessPenalty = 1e7
