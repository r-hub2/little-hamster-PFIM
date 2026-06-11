# Report helpers shared by Evaluation and Optimization.

.buildModelParametersKable = function( modelParameters ) {
  df = modelParameters |>
    map( getModelParametersData ) |>
    map( ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |>
    list_rbind()
  colnames( df ) = c(
    "Parameter",
    "$\\mu$",
    "$\\omega^2$",
    "$\\gamma^2$",
    "Distribution",
    paste0( "$\\mu$",      " fixed" ),
    paste0( "$\\omega^2$", " fixed" )
  )
  kbl( df, align = c( "l","l","l","l","c","c","c" ) ) |>
    kable_styling( bootstrap_options = c( "hover" ),
                   full_width = FALSE, position = "center", font_size = 13 )
}

# Build the covariates structure kable.
.buildCovariatesKable = function( modelCovariates ) {
  if ( length( modelCovariates ) == 0L ) return( NULL )

  df = map( modelCovariates, function( cov ) {
    nm    = prop( cov, "name" )
    isIOV = S7::S7_inherits( cov, CategoricalCovariateWithIOV )

    if ( isIOV ) {
      seqs      = prop( cov, "sequences"            )
      seqProps  = prop( cov, "sequencesProportions" )
      seqNames  = if ( is.null( names( seqs ) ) ) paste0( "seq", seq_along( seqs ) ) else names( seqs )

      modalStr  = paste(
        map2_chr( seqNames, seqs,
                  ~ sprintf( "%s: %s", .x, paste( unlist( .y ), collapse = "\u2192" ) ) ),
        collapse = "  |  "
      )
      propStr   = paste( round( unlist( seqProps ), 3 ), collapse = " | " )
      data.frame( Covariate = nm, Type = "IOV",
                  Modalities = modalStr, Proportions = propStr,
                  stringsAsFactors = FALSE )
    } else {
      cats  = prop( cov, "categories"             )
      props = prop( cov, "categoriesProportions"  )
      data.frame( Covariate   = nm,
                  Type        = "Categorical",
                  Modalities  = paste( unlist( cats ),                    collapse = " | " ),
                  Proportions = paste( round( unlist( props ), 3 ),       collapse = " | " ),
                  stringsAsFactors = FALSE )
    }
  } ) |> list_rbind()

  kbl( df, align = c( "l","c","l","l" ) ) |>
    kable_styling( bootstrap_options = c( "hover", "bordered" ),
                   full_width = FALSE, position = "center", font_size = 13 )
}

# Run covariateTest() on pfimproject and return the kable tables
.buildCovariateTestSection = function( pfimproject ) {
  tryCatch(
    getCovariateTestTables( covariateTest( pfimproject ) ),
    error = function( e ) NULL
  )
}
