#' @title Design
#' @description
#' Experimental design: one or more arms, each with dosing and sampling schedules.
#' Filled by \code{evaluateDesign()} with per-arm FIM results.
#' @param name Character string: design name.
#' @param size Total number of subjects (sum of arm sizes).
#' @param arms List of \code{Arm} objects.
#' @param numberOfArms Optional total study size for discrete designs (historical
#'   CRAN field; Mult/FW prefer \code{size} / arm sizes / \code{numberOfSubjects}).
#'   Not the count of \code{Arm} objects.
#' @param evaluationArms Evaluated arms after \code{evaluateDesign()}.
#' @param fim Aggregated \code{Fim} for the design.
#' @include Fim.R
#' @include Arm.R
#' @return An S7 object of class \code{Design}.
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' design = prop(ev, "designs")[[1L]]
#' prop(prop(design, "arms")[[1L]], "size")
#' }
#' @export

Design = new_class("Design", package = "PFIM",
                   properties = list(
                     name = new_property(class_character, default = character(0)),
                     size = new_property(class_double, default = 0.0),
                     arms = new_property(class_list, default = list()),
                     evaluationArms = new_property(class_list, default = list()),
                     numberOfArms = new_property(class_double, default = 0.0),
                     fim = new_property(NULL | Fim, default = NULL)
                   ),
                   validator = function( self ) {
                     n = prop( self, "name" )
                     if ( length( n ) > 1L )
                       return( "Design: name must be a single string." )
                     if ( .pfimIsBlankScalar( n ) )
                       return( "Design: name must be a non-empty string." )
                     size = prop( self, "size" )
                     if ( length( size ) == 1L && ( !is.finite( size ) || size < 0 ) )
                       return( "Design: size must be non-negative." )
                     nArms = prop( self, "numberOfArms" )
                     if ( length( nArms ) == 1L && ( !is.finite( nArms ) || nArms < 0 ) )
                       return( "Design: numberOfArms must be non-negative." )
                     msg = .validateS7List(
                       prop( self, "arms" ), Arm, "Design:arms", "Arm"
                     )
                     if ( !is.null( msg ) )
                       return( msg )
                     msg = .validateS7List(
                       prop( self, "evaluationArms" ), Arm, "Design:evaluationArms", "Arm"
                     )
                     if ( !is.null( msg ) )
                       return( msg )
                     NULL
                   })


#' Sampling times for one outcome on an arm.
#' @param samplingTimes List of \code{SamplingTimes} objects.
#' @param outcome Character outcome name.
#' @return Numeric vector of sampling times.
#' @noRd
#' @keywords internal
.samplingsForOutcome = function( samplingTimes, outcome ) {
  keep( samplingTimes, ~ prop( .x, "outcome" ) == outcome ) |>
    map( ~ prop( .x, "samplings" ) ) |>
    unlist()
}


#' Arm sizes for design FIM mixture weights (\code{<= 0} -> drop).
#' @noRd
#' @keywords internal
.designArmSizes = function( evaluationArms ) {
  map_dbl( evaluationArms, function( arm ) {
    s = prop( arm, "size" )
    if ( length( s ) != 1L || !is.finite( s ) ) 0 else as.numeric( s )
  } )
}

#' Assemble the design-level FIM from evaluated arms.
#'
#' \strong{Population:} arm FIMs are already \eqn{\times n_a}; the design matrix is
#' their sum.
#'
#' \strong{Individual / Bayesian:} each arm stores a \emph{per-subject} FIM
#' (Bayesian arms keep their prior). Across arms, average covariances then invert:
#' \deqn{\bar C = \sum_a (n_a / N)\, M_a^{-1},\quad M_{\mathrm{eff}} = \bar C^{-1}.}
#' Arms with \code{size <= 0} are ignored. Shrinkage uses \eqn{M_{\mathrm{eff}}}.
#' @noRd
#' @keywords internal
.assembleDesignFim = function( fimPrototype, evaluationArms, model = NULL ) {

  designFim = .duplicateFim( fimPrototype )
  if ( !length( evaluationArms ) ) {
    prop( designFim, "fisherMatrix" ) = matrix( numeric( 0L ), 0L, 0L )
    return( designFim )
  }

  sizes = .designArmSizes( evaluationArms )
  keep  = sizes > 0
  # If every size is missing/zero (malformed design), fall back to equal weights.
  if ( !any( keep ) ) {
    keep  = rep( TRUE, length( evaluationArms ) )
    sizes = rep( 1, length( evaluationArms ) )
  }
  evaluationArms = evaluationArms[ keep ]
  sizes          = sizes[ keep ]
  N              = sum( sizes )
  weights        = sizes / N

  fisherBlocks = map(
    evaluationArms,
    ~ prop( prop( .x, "evaluationFim" ), "fisherMatrix" )
  )

  isBayes = S7::S7_inherits( fimPrototype, BayesianFim )
  isPop   = S7::S7_inherits( fimPrototype, PopulationFim )

  if ( isPop ) {
    # Population: each arm FIM is already x n_a.
    M = Reduce( `+`, fisherBlocks )
  } else {
    # Individual / Bayesian: arms keep prior; aggregate by covariance mixture.
    M = .pfimHarmonicMeanFim( fisherBlocks, weights )$fisherMatrix
  }

  if ( is.matrix( M ) && nrow( M ) == ncol( M ) && length( M ) )
    M = 0.5 * ( M + t( M ) )
  prop( designFim, "fisherMatrix" ) = M

  # Shrinkage from the design FIM (same matrix as SE/RSE).
  if ( isBayes && !is.null( model ) ) {
    prop( designFim, "shrinkage" ) = .bayesianShrinkage(
      M, model, evaluationArms[[ 1L ]]
    )
  }

  designFim
}

#' Evaluate a design for a model and FIM type
#' @param design A \code{Design} object.
#' @param ... Method-specific arguments. For \code{Design}, \code{model} and \code{fim} objects.
#' @usage evaluateDesign(design, ...)
#' @return The updated \code{Design} object after evaluation.
#' @examples
#' \donttest{
#' source(system.file("examples", "evaluation-minimal.R", package = "PFIM"))
#' length(prop(ev, "evaluationDesign"))
#' }
#' @name evaluateDesign
#' @keywords internal
evaluateDesign = new_generic( "evaluateDesign", c( "design" ) )

#' Enumerate dose combinations under design constraints
#' @param design A \code{Design} object.
#' @param ... Optional method arguments.
#' @name generateDosesCombination
#' @return A list or matrix of generated dose combinations.
#' @keywords internal
generateDosesCombination = new_generic( "generateDosesCombination", c( "design" ) )

#' Enumerate sampling-time combinations under constraints
#' @param design A \code{Design} object.
#' @param ... Optional method arguments.
#' @name generateSamplingTimesCombination
#' @return A list or matrix of generated sampling-time combinations.
#' @keywords internal
generateSamplingTimesCombination = new_generic( "generateSamplingTimesCombination", c( "design" ) )

#' Hard cap on discrete dose/sampling Cartesian products.
#'
#' Enumerating \code{combn()} / \code{expand.grid()} without a bound can exhaust
#' memory on large constraint grids. This limit applies to the *enumeration*
#' step only; it is independent of \code{constraints.maxTasks}, which only
#' subsamples FIM evaluations *after* the grid exists.
#' @noRd
#' @keywords internal
.pfimMaxCombinationCount = 1000000L

#' Stop when a discrete enumeration would exceed \code{.pfimMaxCombinationCount}.
#'
#' @param count Number of combinations about to be materialised (must be
#'   \code{nrow} for dose grids, or product of per-outcome list lengths for
#'   sampling grids - never \code{prod(dim(...))}).
#' @param context Short label included in the error (e.g. \code{"Dose enumeration"}).
#' @noRd
#' @keywords internal
.pfimStopIfCombinationOverflow = function( count, context ) {
  if ( !is.finite( count ) || count <= .pfimMaxCombinationCount )
    return( invisible( NULL ) )
  stop(
    context, " would enumerate ", format( count, scientific = FALSE ), " combinations ",
    "(limit ", format( .pfimMaxCombinationCount, scientific = FALSE ), "). ",
    "Reduce the number of candidate doses or sampling times.",
    call. = FALSE
  )
}

#' Validate generated sampling constraints for feasibility
#' @param design A \code{Design} object.
#' @param ... Optional method arguments.
#' @name checkValiditySamplingConstraint
#' @return Logical value indicating whether sampling constraints are valid.
#' @keywords internal
checkValiditySamplingConstraint = new_generic( "checkValiditySamplingConstraint", c( "design" ) )

#' Apply sampling constraints for optimization algorithms
#' @param design A \code{Design} object.
#' @param ... Optional method arguments.
#' @name setSamplingConstraintForOptimization
#' @return The modified \code{Design} object.
#' @keywords internal
setSamplingConstraintForOptimization = new_generic( "setSamplingConstraintForOptimization", c( "design" ) )

#' Evaluate every arm and assemble the design-level FIM.
#'
#' Runs \code{evaluateArm()} on each arm (the arm method clones \code{model}
#' before binding administration), then \code{.assembleDesignFim()}.
#' @return The same \code{Design} with \code{evaluationArms} and \code{fim} set.
#' @name evaluateDesign
#' @usage NULL
method( evaluateDesign, Design ) = function( design, model, fim ) {

  arms = prop( design, "arms" )
  prop( design, "evaluationArms" ) = map(
    arms,
    function( arm ) evaluateArm( arm, model, fim )
  )
  prop( design, "fim" ) = .assembleDesignFim( fim, prop( design, "evaluationArms" ), model )

  design
}

#' Enumerate dose levels under administration constraints.
#'
#' Builds the Cartesian product of dose lists across arms/outcomes, guarded by
#' \code{.pfimMaxCombinationCount}. Returns a nested list keyed by arm and
#' outcome, plus \code{numberOfDoses}.
#' @name generateDosesCombination
#' @keywords internal

method( generateDosesCombination, Design ) = function( design ) {

  arms = prop( design, "arms" )

  # Cartesian product of dose levels across administrations (one row = one combo).
  # Missing AdministrationConstraints: singleton grid from the arm's current doses
  # (optimize sampling times only).
  dosesForFIMsTmp = map( arms, function( arm ) {
    cons = prop( arm, "administrationsConstraints" )
    adms = prop( arm, "administrations" )
    if ( length( cons ) ) {
      if ( length( cons ) != length( adms ) )
        .pfimStop(
          "Arm '", prop( arm, "name" ),
          "': administrationsConstraints must have one entry per Administration (",
          length( adms ), " administration(s), ", length( cons ), " constraint(s))."
        )
      return( map( cons, ~ prop( .x, "doses" ) ) )
    }
    map( adms, function( adm ) {
      d = as.numeric( prop( adm, "dose" ) )
      if ( !length( d ) || any( !is.finite( d ) ) )
        .pfimStop(
          "Arm '", prop( arm, "name" ),
          "': Administration(outcome = \"", prop( adm, "outcome" ),
          "\") has no finite dose. Set dose= or add AdministrationConstraints."
        )
      as.list( d[[ 1L ]] )
    } )
  } ) |> flatten()

  if ( !length( dosesForFIMsTmp ) )
    .pfimStop(
      "Fedorov-Wynn / Multiplicative: no doses to enumerate. ",
      "Add Administration() on the arm, or AdministrationConstraints(doses = ...)."
    )

  dosesForFIMsTmp = expand.grid( dosesForFIMsTmp )

  # Guard uses nrow (= number of dose combinations), not prod(dim) = nrow*ncol.
  .pfimStopIfCombinationOverflow(
    nrow( dosesForFIMsTmp ),
    "Dose enumeration"
  )

  admin_specs = map( arms, function( arm ) {
    armName = prop( arm, "name" )
    map( prop( arm, "administrations" ), function( administration ) {
      list( armName = armName, outcome = prop( administration, "outcome" ) )
    } )
  } ) |> flatten()

  dosesForFIMs = reduce(
    imap( admin_specs, ~ list( spec = .x, col = .y ) ),
    function( acc, item ) {
      acc[[ item$spec$armName ]][[ item$spec$outcome ]] =
        unlist( dosesForFIMsTmp[, item$col, drop = TRUE ] )
      acc
    },
    .init = list()
  )
  c( dosesForFIMs, numberOfDoses = dim( dosesForFIMsTmp )[1L] )
}

#' Enumerate sampling-time combinations under sampling constraints.
#'
#' For each arm and outcome: \code{combn} of free times plus fixed times, then
#' \code{expand.grid} across outcomes (guarded before materialising). Returns a
#' named list of arms -> list of sampling schedules for FIM evaluation.
#' @name generateSamplingTimesCombination
#' @keywords internal

method( generateSamplingTimesCombination, Design ) = function( design ) {

  arms = prop( design, "arms" )

  armResults = map( arms, function( arm ) {

    armName = prop( arm, "name" )
    samplingTimesConstraints = prop( arm, "samplingTimesConstraints" )
    samplingTimes = prop( arm, "samplingTimes")
    outcomeNames = map_chr( samplingTimesConstraints, ~ prop( .x, "outcome" ) )

    # Per outcome: choose k free times from the candidate grid (+ fixed times).
    samplingTimesCombinations = map( samplingTimesConstraints, function( samplingTimesConstraint ) {
      initialSamplings = prop( samplingTimesConstraint, "initialSamplings" )
      fixedTimes = prop( samplingTimesConstraint, "fixedTimes" )
      numberOfsamplingsOptimisable = prop( samplingTimesConstraint, "numberOfsamplingsOptimisable" )
      availableSamplings = setdiff( initialSamplings, fixedTimes )
      outcome = prop( samplingTimesConstraint, "outcome" )
      k = as.integer( numberOfsamplingsOptimisable - length( fixedTimes ) )
      if ( k < 0L )
        stop(
          "numberOfsamplingsOptimisable cannot be less than the number of fixed times for outcome ",
          outcome, ".", call. = FALSE
        )
      if ( k > length( availableSamplings ) )
        stop(
          "Not enough available sampling times for outcome ", outcome,
          " (need ", k, ", have ", length( availableSamplings ), ").", call. = FALSE
        )
      combinations = if ( k == 0L ) list( fixedTimes ) else combn( availableSamplings, k, simplify = FALSE )
      map( combinations, ~ c( fixedTimes, .x ) )
    })

    # Guard before expand.grid (Cartesian product across outcomes).
    n_combinations = prod( vapply( samplingTimesCombinations, length, integer( 1L ) ) )
    .pfimStopIfCombinationOverflow(
      n_combinations,
      paste0( "Sampling enumeration for arm '", armName, "'" )
    )

    samplingTimesCombinations = expand.grid( samplingTimesCombinations )
    colnames( samplingTimesCombinations ) = outcomeNames
    samplingTimesCombinations = pmap( samplingTimesCombinations, ~ list( ... ) )

    # Materialise SamplingTimes clones with the chosen times (sorted).
    samplingsForFIM = map( samplingTimesCombinations, function( samplingTimeCombination ) {
      map2( samplingTimes, outcomeNames, function( samplingTime, outcomeName ) {
        st = .pfimCloneS7( samplingTime )
        prop( st, "samplings" ) = sort( samplingTimeCombination[[outcomeName]] )
        st
      })
    })
  })
  set_names( armResults, map_chr( arms, ~ prop( .x, "name" ) ) )
}

#' Fail fast if window counts / spacing cannot fit the arm's sampling times.
#'
#' Same slack test as \code{generateSamplingsFromSamplingConstraints}:
#' \code{max - min - (n-1)*delta >= 0}. Also requires
#' \code{sum(numberOfTimesByWindows) == length(samplings)} per outcome, and
#' that the current sampling times actually occupy the declared windows
#' (counts and \code{minSampling}).
#' @name checkValiditySamplingConstraint
#' @keywords internal

method( checkValiditySamplingConstraint, Design ) = function( design ) {

  arms = prop( design, "arms" )

  walk ( arms, function( arm ) {
    armName = prop( arm, "name" )

    samplingTimesConstraints = prop( arm, "samplingTimesConstraints" )
    outcomes = map( samplingTimesConstraints, ~ prop( .x, "outcome" ) ) |> unlist()

    samplingsWindow = map( samplingTimesConstraints, ~ prop( .x, "samplingsWindows" ) ) |> stats::setNames( outcomes )
    numberOfTimesByWindows = map( samplingTimesConstraints, ~ prop( .x, "numberOfTimesByWindows" ) ) |> stats::setNames( outcomes )
    minSampling = map( samplingTimesConstraints, ~ prop( .x, "minSampling" ) ) |> stats::setNames( outcomes )

    inputRandomSpaced = list()
    samplingTimesArms = list()

    walk ( outcomes, function( outcome ) {
      intervalsConstraints = list()

      samplingTimes = prop( arm, "samplingTimes")
      samplings = .samplingsForOutcome( samplingTimes, outcome )

      minSamplingAndNumberOfTimesByWindows = as.data.frame( list( minSampling[[outcome]], numberOfTimesByWindows[[outcome]] ) )
      tmp = t( as.data.frame(samplingsWindow[[outcome]] ) )
      inputRandomSpaced[[outcome]] = as.data.frame( do.call( "cbind", list( tmp, minSamplingAndNumberOfTimesByWindows ) ) )

      colnames( inputRandomSpaced[[outcome]] ) = c("min","max","delta","n")
      rownames( inputRandomSpaced[[outcome]] ) = NULL

      # Window counts must partition the current sampling vector exactly.
      if ( sum( numberOfTimesByWindows[[outcome]] ) != length( samplings ) ) {
        stop(
          "The sampling times constraint is not possible for arm ",
          armName, " and outcome ", outcome, ".", call. = FALSE
        )
      }

      walk( seq_len( length( inputRandomSpaced[[outcome]]$n)), function( iter ) {
        min = inputRandomSpaced[[outcome]]$min[iter]
        max = inputRandomSpaced[[outcome]]$max[iter]
        delta = inputRandomSpaced[[outcome]]$delta[iter]
        n = inputRandomSpaced[[outcome]]$n[iter]

        # Negative slack => n points cannot fit with min gap delta in [min, max].
        distance = max-min-(n-1)*delta

        if ( distance < 0 ) {
          stop(
            "The sampling times constraint is not possible for arm ",
            armName, " and outcome ", outcome, ".", call. = FALSE
          )
        }
      })

      sc = pluck(
        keep( samplingTimesConstraints, ~ prop( .x, "outcome" ) == outcome ),
        1L
      )
      occ = .pfimTimesWindowError( sc, samplings, label = "sampling times" )
      if ( !is.null( occ ) )
        .pfimStop( "arm '", armName, "', outcome '", outcome, "': ", occ )
    })
  })
}

#' Fill missing sampling constraints from the arm's current sampling times.
#'
#' Outcomes that have \code{SamplingTimes} but no \code{SamplingTimeConstraints}
#' get a single window \code{[min, max]}, all points free (\code{minSampling = 0}).
#' Optimizers then see one constraint object per sampled outcome.
#' @name setSamplingConstraintForOptimization
#' @keywords internal

method( setSamplingConstraintForOptimization, Design ) = function( design ) {

  arms = prop( design, "arms" )

  arms = map ( arms, function( arm )
  {
    samplingTimes = prop( arm, "samplingTimes" )
    outcomes = map_chr( samplingTimes, ~ prop( .x, "outcome" ) )

    samplingTimesConstraints = prop( arm, "samplingTimesConstraints" )
    outcomesSamplingTimesConstraints = map_chr( samplingTimesConstraints, ~ prop( .x, "outcome" ) )
    outcomesSamplingNotInTimesConstraints = outcomes[!outcomes %in% outcomesSamplingTimesConstraints]

    if ( length( outcomesSamplingNotInTimesConstraints ) !=0 )
    {
      samplingTimesConstraints = reduce(
        outcomesSamplingNotInTimesConstraints,
        function( stc, outcomeSamplingNotInTimesConstraints ) {
          samplings = .samplingsForOutcome( samplingTimes, outcomeSamplingNotInTimesConstraints )
          # Default: one window over the existing times, no min-gap (fully free).
          newSamplingTimeConstraints = SamplingTimeConstraints( outcome = outcomeSamplingNotInTimesConstraints,
                                                                initialSamplings = samplings,
                                                                samplingsWindows = list( c( min( samplings ), max( samplings ) ) ),
                                                                numberOfTimesByWindows = length( samplings ),
                                                                minSampling = 0 )
          append( stc, newSamplingTimeConstraints )
        },
        .init = samplingTimesConstraints
      )
    }
    prop( arm, "samplingTimesConstraints" ) = samplingTimesConstraints
    arm
  })

  prop( design, "arms" ) = arms
  design
}






