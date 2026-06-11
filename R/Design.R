#' @title Design
#' @description
#' Experimental design: one or more arms, each with dosing and sampling schedules.
#' Filled by \code{evaluateDesign()} with per-arm FIM results.
#' @param name Character string: design name.
#' @param size Total number of subjects (sum of arm sizes).
#' @param arms List of \code{Arm} objects.
#' @param numberOfArms Number of arms (computed).
#' @param evaluationArms Evaluated arms after \code{evaluateDesign()}.
#' @param fim Aggregated \code{Fim} for the design.
#' @include Fim.R
#' @export

Design = new_class("Design", package = "PFIM",
                   properties = list(
                     name = new_property(class_character, default = character(0)),
                     size = new_property(class_double, default = 0.0),
                     arms = new_property(class_list, default = list()),
                     evaluationArms = new_property(class_list, default = list()),
                     numberOfArms = new_property(class_double, default = 0.0),
                     fim = new_property(Fim, default = NULL)
                   ))


#' Sum arm-level FIM matrices into a new design-level \code{Fim} object.
#'
#' @param fimPrototype Template from \code{defineFim()} (class only; not mutated).
#' @param evaluationArms List of evaluated \code{Arm} objects.
#' @return Aggregated \code{Fim} for the design.
#' @keywords internal
.assembleDesignFim = function( fimPrototype, evaluationArms ) {

  designFim = .duplicateFim( fimPrototype )
  armFims   = map( evaluationArms, ~ prop( .x, "evaluationFim" ) )

  fisherBlocks = map( armFims, ~ prop( .x, "fisherMatrix" ) )
  prop( designFim, "fisherMatrix" ) = Reduce( `+`, fisherBlocks )

  # Shrinkage
  shrinkage = prop( armFims[[ 1L ]], "shrinkage" )
  if ( length( shrinkage ) > 0L )
    prop( designFim, "shrinkage" ) = shrinkage

  designFim
}

evaluateDesign = new_generic( "evaluateDesign", c( "design" ) )
generateDosesCombination = new_generic( "generateDosesCombination", c( "design" ) )
generateSamplingTimesCombination = new_generic( "generateSamplingTimesCombination", c( "design" ) )
checkValiditySamplingConstraint = new_generic( "checkValiditySamplingConstraint", c( "design" ) )
setSamplingConstraintForOptimization = new_generic( "setSamplingConstraintForOptimization", c( "design" ) )

#' Evaluate Fisher information for a design
#' @name evaluateDesign
#' @export

method( evaluateDesign, Design ) = function( design, model, fim ) {

  arms = prop( design, "arms" )
  prop( design, "evaluationArms" ) = map(
    arms,
    function( arm ) evaluateArm( arm, model, fim )
  )
  prop( design, "fim" ) = .assembleDesignFim( fim, prop( design, "evaluationArms" ) )

  design
}

#' Enumerate dose levels under constraints
#' @name generateDosesCombination
#' @export

method( generateDosesCombination, Design ) = function( design ) {

  arms = prop( design, "arms" )
  armNames = map_chr( arms, ~ prop( .x, "name" ) )
  outcomes = map( arms, ~ map_chr(prop( .x, "administrationsConstraints" ), ~ prop( .x,"outcome" ) ) )

  # combination of the doses
  dosesForFIMsTmp = map( arms, function( arm ) {
    map( prop( arm, "administrationsConstraints" ), ~ prop( .x, "doses" ) )
  } ) |> flatten() |> expand.grid()

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

#' Enumerate sampling combinations under constraints
#' @name generateSamplingTimesCombination
#' @export

method( generateSamplingTimesCombination, Design ) = function( design ) {

  arms = prop( design, "arms" )
  armNames = map_chr( arms, ~ prop( .x, "name" ) )

  armResults = map( arms, function( arm ) {

    # Extract the relevant properties for each arm
    armName = prop( arm, "name" )
    samplingTimesConstraints = prop( arm, "samplingTimesConstraints" )
    samplingTimes = prop( arm, "samplingTimes")
    outcomeNames = map_chr( samplingTimesConstraints, ~ prop( .x, "outcome" ) )

    # Generate all combinations of sampling times for the arm
    samplingTimesCombinations = map( samplingTimesConstraints, function( samplingTimesConstraint ) {
      initialSamplings = prop( samplingTimesConstraint, "initialSamplings" )
      fixedTimes = prop( samplingTimesConstraint, "fixedTimes" )
      numberOfsamplingsOptimisable = prop( samplingTimesConstraint, "numberOfsamplingsOptimisable" )
      availableSamplings = setdiff( initialSamplings, fixedTimes )
      combinations = combn( availableSamplings, numberOfsamplingsOptimisable - length( fixedTimes ), simplify = FALSE )

      # Generate all combinations of sampling times
      map( combinations, ~ c( fixedTimes, .x ) )
    })

    # Flatten the list and convert to a data frame
    samplingTimesCombinations = expand.grid( samplingTimesCombinations )
    colnames( samplingTimesCombinations ) = outcomeNames

    # Convert to a list of named lists
    samplingTimesCombinations = pmap( samplingTimesCombinations, ~ list( ... ) )

    # Map the sampling times combinations to the sampling times for each outcome
    samplingsForFIM = map( samplingTimesCombinations, function( samplingTimeCombination ) {
      map2( samplingTimes, outcomeNames, function( samplingTime, outcomeName ) {
        prop( samplingTime, "samplings" ) = sort( samplingTimeCombination[[outcomeName]] )
        samplingTime
      })
    })
  })
  set_names( armResults, armNames )
}

#' Check feasibility of sampling constraints
#' @name checkValiditySamplingConstraint
#' @export

method( checkValiditySamplingConstraint, Design ) = function( design ) {

  arms = prop( design, "arms" )

  walk ( arms, function( arm ) {
    armName = prop( arm, "name" )

    # get the outcomes
    samplingTimesConstraints = prop( arm, "samplingTimesConstraints" )
    outcomes = map( samplingTimesConstraints, ~ prop( .x, "outcome" ) ) |> unlist()

    # get samplings window constraints
    samplingsWindow = map( samplingTimesConstraints, ~ prop( .x, "samplingsWindows" ) ) |> setNames( outcomes )

    # get numberOfTimesByWindows constraints
    numberOfTimesByWindows = map( samplingTimesConstraints, ~ prop( .x, "numberOfTimesByWindows" ) ) |> setNames( outcomes )

    # get minimal time step for each windows
    minSampling = map( samplingTimesConstraints, ~ prop( .x, "minSampling" ) ) |> setNames( outcomes )

    inputRandomSpaced = list()
    samplingTimesArms = list()

    walk ( outcomes, function( outcome ) {
      intervalsConstraints = list()

      # get samplingTimes and samplings
      samplingTimes = prop( arm, "samplingTimes")
      samplings = map( samplingTimes, ~ if ( prop( .x, "outcome") == outcome ) prop(.x ,"samplings" ) ) |> compact() |> unlist()

      minSamplingAndNumberOfTimesByWindows = as.data.frame( list( minSampling[[outcome]], numberOfTimesByWindows[[outcome]] ) )
      tmp = t( as.data.frame(samplingsWindow[[outcome]] ) )
      inputRandomSpaced[[outcome]] = as.data.frame( do.call( "cbind", list( tmp, minSamplingAndNumberOfTimesByWindows ) ) )

      colnames( inputRandomSpaced[[outcome]] ) = c("min","max","delta","n")
      rownames( inputRandomSpaced[[outcome]] ) = NULL

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

        distance = max-min-(n-1)*delta

        if ( distance < 0 ) {
          stop(
            "The sampling times constraint is not possible for arm ",
            armName, " and outcome ", outcome, ".", call. = FALSE
          )
        }
      })
    })
  })
}

#' Prepare sampling constraints for optimization
#' @name setSamplingConstraintForOptimization
#' @export

method( setSamplingConstraintForOptimization, Design ) = function( design ) {

  arms = prop( design, "arms" )

  arms = map ( arms, function( arm )
  {
    # get the outcomes in the sampling times
    samplingTimes = prop( arm, "samplingTimes" )
    outcomes = map_chr( samplingTimes, ~ prop( .x, "outcome" ) )

    # set the sampling time constraints for missing outcomes ie from its sampling times
    samplingTimesConstraints = prop( arm, "samplingTimesConstraints" )
    outcomesSamplingTimesConstraints = map_chr( samplingTimesConstraints, ~ prop( .x, "outcome" ) )
    outcomesSamplingNotInTimesConstraints = outcomes[!outcomes %in% outcomesSamplingTimesConstraints]

    if ( length( outcomesSamplingNotInTimesConstraints ) !=0 )
    {
      samplingTimesConstraints = reduce(
        outcomesSamplingNotInTimesConstraints,
        function( stc, outcomeSamplingNotInTimesConstraints ) {
          samplings = map( samplingTimes, ~ if ( prop( .x, "outcome") == outcomeSamplingNotInTimesConstraints ) prop(.x ,"samplings" ) ) |> compact() |> unlist()
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
    return(arm)
  })

  # set the design with arms
  prop( design, "arms" ) = arms
  return( design )
}






