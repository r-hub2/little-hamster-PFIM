# kable helpers for Evaluation / Optimization HTML reports.
# Escape is off so LaTeX math in column headers (\mu, \sigma, ...) renders in HTML.

#' Styled \code{kableExtra} table for HTML evaluation/optimization reports.
#' @param x Data frame or matrix passed to \code{kable()}.
#' @param bootstrap_options Bootstrap table classes for \code{kable_styling()}.
#' @param ... Additional arguments forwarded to \code{kable()}.
#' @return A \code{kableExtra} object with striped hover styling.
#' @noRd
#' @keywords internal
.kblReportStyled = function(
    x,
    ...,
    bootstrap_options = c( "striped", "hover" ) ) {
  kbl( x, escape = FALSE, ... ) |>
    kable_styling(
      bootstrap_options = bootstrap_options,
      full_width        = FALSE,
      position          = "center",
      font_size         = 13
    )
}

#' Column align vector for report kables: \code{n_left} left, rest centered.
#' @noRd
#' @keywords internal
.kblAlign = function( n_col, n_left = 1L ) {
  n_col  = as.integer( n_col )
  n_left = max( 0L, min( as.integer( n_left ), n_col ) )
  c( rep( "l", n_left ), rep( "c", max( n_col - n_left, 0L ) ) )
}

#' Drop row names before report kable (SE/RSE tables).
#' @noRd
#' @keywords internal
.pfimKableDataFrame = function( df ) {
  row.names( df ) = NULL
  df
}

#' Log-determinant / D-criterion / condition-number table for evaluation reports.
#'
#' Bayesian FIMs omit the variance condition-number column (\code{condVariance}
#' left \code{NULL}); population/individual reports pass both condition numbers.
#' @noRd
#' @keywords internal
.fimCriteriaKable = function( logDeterminant, dcriterion, condFixed,
                              condVariance = NULL, singularFim = FALSE ) {
  if ( is.null( condVariance ) ) {
    df = data.frame( ld = logDeterminant, dcrit = dcriterion, fe = condFixed )
    out = kbl(
      df,
      col.names = c( "", "", "Fixed effects" ),
      align     = "c",
      format    = "html"
    ) |>
      add_header_above( c( "log-Determinant" = 1, "D-criterion" = 1, "Condition number" = 1 ) )
  } else {
    df = data.frame(
      ld = logDeterminant, dcrit = dcriterion,
      fe = condFixed, ve = condVariance
    )
    out = kbl(
      df,
      col.names = c( "", "", "Fixed effects", "Variance effects" ),
      align     = "c",
      format    = "html"
    ) |>
      add_header_above( c( "log-Determinant" = 1, "D-criterion" = 1, "Condition number" = 2 ) )
  }
  out = out |>
    kable_styling(
      bootstrap_options = c( "striped", "hover" ),
      full_width        = FALSE,
      position          = "center",
      font_size         = 13
    )
  if ( isTRUE( singularFim ) ) {
    out = footnote(
      out,
      general = paste(
        "singularFim = TRUE: Fisher information matrix was singular or not",
        "positive definite; SE/RSE used a Moore-Penrose pseudo-inverse.",
        "Interpret SE/RSE and design criteria with caution."
      ),
      general_title = "Note: "
    )
    # Visible HTML banner above the criteria table (all Evaluation / Optimization Rmds).
    banner = paste0(
      '<div class="pfim-singular-fim" style="margin:1em 0;padding:0.75em;',
      'border-left:4px solid #b00020;background:#fff5f5;">',
      '<strong>singularFim = TRUE</strong> \u2014 Fisher information matrix was ',
      'singular or not positive definite; SE/RSE used a Moore-Penrose ',
      'pseudo-inverse. Interpret SE/RSE and design criteria with caution.',
      '</div>\n'
    )
    out = structure(
      paste0( banner, paste( as.character( out ), collapse = "\n" ) ),
      class = "knit_asis",
      knit_cacheable = NA
    )
  }
  out
}

#' SE/RSE (and optional shrinkage) table for evaluation reports.
#' @noRd
#' @keywords internal
.fimSeRseKable = function( paramLabels, se_rse_df, shrinkage = NULL ) {
  df = data.frame(
    Parameters = paramLabels,
    round( se_rse_df, 3 ),
    stringsAsFactors = FALSE
  )
  cols = c( "Parameters", "Parameter values", "SE", "RSE (%)" )
  if ( !is.null( shrinkage ) ) {
    df$Shrinkage = round( as.numeric( shrinkage ), 3 )
    cols = c( cols, "Shrinkage (%)" )
  }
  .kblReportStyled(
    .pfimKableDataFrame( df ),
    col.names = cols,
    align     = "c",
    row.names = FALSE
  )
}

#' kable table of residual error model settings.
#' @param modelError List of \code{ModelError} objects.
#' @return A styled \code{kableExtra} table.
#' @noRd
#' @keywords internal
.buildModelErrorKable = function( modelError ) {
  modelErrorData = modelError |>
    map( getModelErrorData ) |>
    map( ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |>
    list_rbind()
  modelErrorData$cError = NULL
  modelErrorData$varianceForm = NULL
  colnames( modelErrorData ) = c(
    "Output", "Type",
    "$\\sigma_{slope}$", "$\\sigma_{inter}$"
  )
  .kblReportStyled( modelErrorData, align = c( "c", "c", "c", "c" ) )
}

#' kable table of population model parameters (mu, omega, optional gamma, fixed flags).
#' @param modelParameters List of \code{ModelParameter} objects.
#' @return A styled \code{kableExtra} table.
#' @noRd
#' @keywords internal
.buildModelParametersKable = function( modelParameters ) {
  has_iov = .hasIovParameters( modelParameters )
  df = modelParameters |>
    map( getModelParametersData ) |>
    map( ~ as.data.frame( .x, stringsAsFactors = FALSE ) ) |>
    list_rbind()
  if ( !has_iov ) df$gamma2 = NULL
  colnames( df ) = if ( has_iov ) {
    c(
      "Parameter", "$\\mu$", "$\\omega^2$", "$\\gamma^2$",
      "Distribution",
      paste0( "$\\mu$", " fixed" ), paste0( "$\\omega^2$", " fixed" )
    )
  } else {
    c(
      "Parameter", "$\\mu$", "$\\omega^2$", "Distribution",
      paste0( "$\\mu$", " fixed" ), paste0( "$\\omega^2$", " fixed" )
    )
  }
  align = if ( has_iov ) c( "l", "l", "l", "l", "c", "c", "c" ) else c( "l", "l", "l", "c", "c", "c" )
  .kblReportStyled( df, align = align )
}

#' Report chunk flags for covariate sections (omit empty headings in HTML).
#' @param hasCov Whether the evaluation defines covariates.
#' @param covariateTestTables Optional list from \code{getCovariateTestTables()}.
#' @param covariatesTable Optional covariate-structure kable from \code{.buildCovariatesKable()}.
#' @return Named logical vector for R Markdown \code{eval} conditions.
#' @noRd
#' @keywords internal
.covReportFlags = function( hasCov, covariateTestTables = NULL, covariatesTable = NULL ) {
  if ( !isTRUE( hasCov ) ) {
    return( list(
      showCovariates          = FALSE,
      showCovTestSignificance = FALSE,
      showCovTestNonRelevance = FALSE,
      showCovTestRelevance    = FALSE
    ) )
  }
  ct = covariateTestTables %||% list()
  list(
    showCovariates          = !is.null( covariatesTable ),
    showCovTestSignificance = !is.null( ct$significance ),
    showCovTestNonRelevance = !is.null( ct$nonRelevance ),
    showCovTestRelevance    = !is.null( ct$relevance )
  )
}

#' kable table of categorical and IOV covariate structure.
#' @param modelCovariates List of covariate objects (may be empty).
#' @return A \code{kableExtra} table, or \code{NULL} when no covariates are defined.
#' @noRd
#' @keywords internal
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

  .kblReportStyled(
    df,
    align = c( "l", "c", "l", "l" ),
    bootstrap_options = c( "striped", "hover", "bordered" )
  )
}

#' Run \code{covariateTest()} and return kable tables for the evaluation report.
#' @param pfimproject A finished \code{Evaluation} \code{PFIMProject} object.
#' @return List of kable tables from \code{getCovariateTestTables()}, or \code{NULL} on error.
#' @noRd
#' @keywords internal
.buildCovariateTestSection = function( pfimproject ) {
  tryCatch(
    getCovariateTestTables( covariateTest( pfimproject ) ),
    error = function( e ) {
      warning(
        "Covariate-test section omitted from report: ",
        conditionMessage( e ),
        call. = FALSE
      )
      NULL
    }
  )
}

#' Convert console FIM dimnames to LaTeX labels for HTML reports.
#'
#' Replaces Greek console prefixes (\code{.greekConsole}) with LaTeX fragments
#' (\code{.greekLatex}) and closes the math mode started by those fragments.
#' @param consoleLabels Character vector from \code{rownames(fixedEffects)} etc.
#' @return A character vector of LaTeX parameter labels for report tables.
#' @noRd
#' @keywords internal
.fimConsoleToLatexLabels = function( consoleLabels ) {
  if ( !length( consoleLabels ) )
    return( character( 0L ) )
  out = reduce(
    names( .greekConsole ),
    function( acc, k )
      gsub( .greekConsole[ k ], .greekLatex[ k ], acc, fixed = TRUE ),
    .init = consoleLabels
  )
  paste0( out, "}$" )
}

#' LaTeX row labels for the fixed-effects FIM block (mu and beta columns).
#' @param evaluation A \code{PFIMProject} object (typically \code{Evaluation}).
#' @param fixedEffects Optional fixed-effects block; when set, labels match its \code{rownames}.
#' @return Character vector of LaTeX labels in column order.
#' @noRd
#' @keywords internal
.fimFixedEffectLatexLabels = function( evaluation, fixedEffects = NULL ) {
  if ( !is.null( fixedEffects ) ) {
    rn = rownames( fixedEffects )
    if ( length( rn ) == nrow( fixedEffects ) )
      return( .fimConsoleToLatexLabels( rn ) )
  }
  fe = .fimFixedEffectLabels( evaluation, .greekLatex )
  c(
    paste0( fe$columnNamesMu,   "}$" ),
    paste0( fe$columnNamesBeta, "}$" )
  )
}

#' LaTeX labels for the full SE/RSE report table (fixed + variance blocks).
#' @param evaluation A \code{PFIMProject} object (typically \code{Evaluation}).
#' @param SEAndRSE SE/RSE data frame from \code{setEvaluationFim()}.
#' @param fixedEffects Optional fixed-effects FIM block.
#' @param varianceEffects Optional variance-effects FIM block.
#' @return Character vector of LaTeX parameter labels.
#' @noRd
#' @keywords internal
.fimSeRseReportLabels = function( evaluation,
                                  SEAndRSE,
                                  fixedEffects    = NULL,
                                  varianceEffects = NULL ) {
  se_rn = rownames( SEAndRSE )
  if ( length( se_rn ) == nrow( SEAndRSE ) )
    return( .fimConsoleToLatexLabels( se_rn ) )

  feL = .fimFixedEffectLatexLabels( evaluation, fixedEffects = fixedEffects )
  veL = .fimConsoleToLatexLabels( rownames( as.matrix( varianceEffects ) ) )
  labels = c( feL, veL )
  if ( length( labels ) != nrow( SEAndRSE ) ) {
    stop( sprintf(
      "Report SE/RSE labels (%d) do not match table rows (%d).",
      length( labels ), nrow( SEAndRSE )
    ), call. = FALSE )
  }
  labels
}

#' Wrap a \code{kableExtra} object for \code{results = "asis"} R Markdown chunks.
#' @param x A \code{kableExtra} object or \code{NULL}.
#' @return \code{knitr::asis_output()} HTML or \code{NULL}.
#' @noRd
#' @keywords internal
.pfimKableAsis = function( x ) {
  if ( is.null( x ) )
    return( NULL )
  html = if ( inherits( x, "knitr_kable" ) ) {
    paste( as.character( x ), collapse = "\n" )
  } else if ( inherits( x, "knit_asis" ) ) {
    as.character( x )
  } else {
    paste( utils::capture.output( print( x ) ), collapse = "\n" )
  }
  if ( !.pfimIsNonEmptyScalar( html ) )
    return( NULL )
  asis_output( html )
}

#' Pre-render kable tables that appear in \code{results = "asis"} report chunks.
#'
#' Converts kable objects to \code{knitr::asis_output()} and clears section flags
#' when the corresponding table is missing, so empty headings are omitted.
#' @param reportTables List passed as an R Markdown \code{params} object.
#' @return Updated \code{reportTables}.
#' @noRd
#' @keywords internal
.pfimPrepareReportAsisTables = function( reportTables ) {
  ct = reportTables$covariateTestTables
  if ( !is.null( reportTables$covariatesTable ) )
    reportTables$covariatesTable = .pfimKableAsis( reportTables$covariatesTable )
  if ( !is.null( ct ) ) {
    ct = imap( ct, function( x, .nm ) {
      if ( is.null( x ) ) x else .pfimKableAsis( x )
    } )
    reportTables$covariateTestTables = ct
  }
  if ( !is.null( reportTables$showCovariates ) )
    reportTables$showCovariates = isTRUE( reportTables$showCovariates ) &&
      !is.null( reportTables$covariatesTable )
  flag_slots = c(
    showCovTestSignificance = "significance",
    showCovTestNonRelevance = "nonRelevance",
    showCovTestRelevance    = "relevance"
  )
  reportTables = reduce( names( flag_slots ), function( rt, flag ) {
    if ( !is.null( rt[[ flag ]] ) )
      rt[[ flag ]] = isTRUE( rt[[ flag ]] ) &&
        !is.null( ct[[ flag_slots[[ flag ]] ]] )
    rt
  }, .init = reportTables )
  reportTables$printPlots = .pfimPrintReportPlots
  reportTables
}
