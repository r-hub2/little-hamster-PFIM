# Shared rmarkdown render path for evaluation / optimization HTML reports.

#' @include Fim.R
#' @importFrom tools file_path_sans_ext
#' @noRd
#' @keywords internal
NULL

#' @noRd
#' @keywords internal
.reportTemplatePath = function( filename ) {
  dir = file.path( system.file( package = "PFIM" ),
                    "rmarkdown", "templates", "skeleton" )
  path = file.path( dir, filename )
  if ( file.exists( path ) ) return( path )
  # Case-insensitive fallback (.Rmd vs .rmd) for installed package layouts.
  stem = sub( "\\.[rR]md$", "", filename )
  hits = list.files( dir,
                      pattern = paste0( "^", stem, "\\.[rR]md$" ),
                      full.names = TRUE, ignore.case = TRUE )
  if ( length( hits ) >= 1L ) return( hits[[ 1L ]] )
  path
}

#' Close graphics devices opened after \code{open_at_start}.
#'
#' Does not close the null device or devices the caller already had open.
#' @noRd
#' @keywords internal
.pfimCloseOwnedGraphicsDevices = function( open_at_start ) {
  now = grDevices::dev.list()
  if ( is.null( now ) ) return( invisible( NULL ) )
  extra = setdiff( now, open_at_start )
  for ( d in extra )
    try( grDevices::dev.off( d ), silent = TRUE )
  invisible( NULL )
}

#' Close every open device except the null device (opt-in via \code{pfim_reset_session}).
#' @noRd
#' @keywords internal
.pfimCloseGraphicsDevices = function() {
  while ( grDevices::dev.cur() > 1L )
    try( grDevices::dev.off(), silent = TRUE )
  invisible( NULL )
}

#' @noRd
#' @keywords internal
.pfimIsKnitting = function() {
  if ( !requireNamespace( "knitr", quietly = TRUE ) ) return( FALSE )
  # knitr::is_knitting exists only in recent knitr; look it up without attaching.
  fn = get0( "is_knitting", envir = asNamespace( "knitr" ), inherits = FALSE )
  is.function( fn ) && isTRUE( fn() )
}

#' Environment where user symbols (global covariates) are resolved.
#'
#' Analytic/ODE wrappers look up free variables such as \code{logtWT} or
#' \code{SEX} here. Prefers \code{knitr::knit_global()} when knitr is available
#' (chunk objects during R Markdown), otherwise \code{globalenv()}.
#' Does not rely on \code{knitr::is_knitting()} (absent in older knitr).
#' @noRd
#' @keywords internal
.pfimUserSymbolEnv = function() {
  if ( requireNamespace( "knitr", quietly = TRUE ) ) {
    kg = get0( "knit_global", envir = asNamespace( "knitr" ), inherits = FALSE )
    if ( is.function( kg ) ) {
      env = tryCatch( kg(), error = function( e ) NULL )
      if ( is.environment( env ) ) return( env )
    }
  }
  globalenv()
}

#' @noRd
#' @keywords internal
.pfimPrepareReportInputForRender = function( template_path, render_dir ) {
  lines = readLines( template_path, warn = FALSE, encoding = "UTF-8" )
  # Nested knit: rename the template's setup chunk to avoid knitr name clash.
  if ( .pfimIsKnitting() && any( grepl( "`\\{r setup([,}])", lines, perl = TRUE ) ) )
    lines = gsub( "`\\{r setup([,}])", "`{r pfim_report_setup\\1", lines, perl = TRUE )
  dest = file.path( render_dir, basename( template_path ) )
  writeLines( lines, dest, useBytes = TRUE )
  dest
}

#' Remove knitr/rmarkdown sidecars for one report only.
#'
#' Deletes \code{\{stem\}_files/} (media folder for this HTML) and optional
#' \code{.knit.md} leftovers. Intentionally does \emph{not} delete every
#' \code{.png} under \code{outputPath}: that would wipe unrelated plots when
#' several reports share a directory.
#' @noRd
#' @keywords internal
.pfimCleanupReportSidecars = function( outputPath, outputFile, template = NULL ) {
  stem  = tools::file_path_sans_ext( basename( outputFile ) )
  # Only this report's media folder --- not every PNG under outputPath.
  media = file.path( outputPath, paste0( stem, "_files" ) )
  if ( dir.exists( media ) )
    unlink( media, recursive = TRUE )
  if ( !is.null( template ) ) {
    knit_md = paste0( tools::file_path_sans_ext( basename( template ) ), ".knit.md" )
    walk(
      unique( c(
        file.path( dirname( .reportTemplatePath( template ) ), knit_md ),
        file.path( outputPath, knit_md )
      ) ),
      ~ if ( file.exists( .x ) ) unlink( .x )
    )
  }
  invisible( NULL )
}

#' @noRd
#' @keywords internal
.pfimRenderReport = function( template, outputFile, outputPath, reportTables ) {
  if ( !requireNamespace( "rmarkdown", quietly = TRUE ) )
    .pfimStop( "Report() needs the 'rmarkdown' package." )
  if ( !rmarkdown::pandoc_available() )
    .pfimStop( "Report() needs pandoc (https://pandoc.org)." )
  reportTables = .pfimPrepareReportAsisTables( reportTables )
  pn = reportTables$projectName
  if ( !.pfimIsNonEmptyScalar( pn ) )
    reportTables$projectName = "PFIM Report"
  open_at_start = grDevices::dev.list()
  chunk_old = NULL
  if ( requireNamespace( "knitr", quietly = TRUE ) ) {
    chunk_old = knitr::opts_chunk$get( c( "dpi", "fig.retina" ) )
    knitr::opts_chunk$set( dpi = 96, fig.retina = 1 )
  }
  render_dir = tempfile( pattern = "pfim_report" )
  dir.create( render_dir, showWarnings = FALSE )
  on.exit( {
    unlink( render_dir, recursive = TRUE )
    if ( !is.null( chunk_old ) )
      do.call( knitr::opts_chunk$set, chunk_old )
    .pfimCloseOwnedGraphicsDevices( open_at_start )
  }, add = TRUE )
  template_path = .reportTemplatePath( template )
  input = .pfimPrepareReportInputForRender( template_path, render_dir )
  # When already knitting, pin the knit root so relative paths stay inside render_dir.
  knit_root_dir = if ( .pfimIsKnitting() ) render_dir else NULL
  out = rmarkdown::render(
    input          = input,
    output_file    = outputFile,
    output_dir     = render_dir,
    knit_root_dir  = knit_root_dir,
    params         = list( reportTables = reportTables ),
    # Self-contained HTML embeds figures; no external _files/ dependency for users.
    output_options = list( html_document = list( self_contained = TRUE ) ),
    envir          = new.env( parent = globalenv() ),
    quiet          = TRUE
  )
  dir.create( outputPath, recursive = TRUE, showWarnings = FALSE )
  file.copy(
    file.path( render_dir, outputFile ),
    file.path( outputPath, outputFile ),
    overwrite = TRUE
  )
  .pfimCleanupReportSidecars( render_dir, outputFile, template )
  .pfimCleanupReportSidecars( outputPath, outputFile )
  invisible( out )
}

#' @noRd
#' @keywords internal
.renderReport = function( template ) {
  # Capture template name; return an S7 method body for generateReportOptimization.
  force( template )
  function( fim, optimizationAlgorithm, tablesForReport, outputFile, outputPath )
    .pfimRenderReport( template, outputFile, outputPath, tablesForReport )
}

#' @noRd
#' @keywords internal
.renderEvalReport = function( template ) {
  # Same factory as .renderReport, but for generateReportEvaluation (no optimizer arg).
  force( template )
  function( fim, tablesForReport, outputFile, outputPath )
    .pfimRenderReport( template, outputFile, outputPath, tablesForReport )
}

#' @noRd
#' @keywords internal
.pfimRegisterOptimizationReports = function( fimClass, fimSuffix ) {
  # Wire each optimizer class to its matching Rmd template for this FIM type.
  algos = list(
    MultiplicativeAlgorithm, FedorovWynnAlgorithm, SimplexAlgorithm,
    PSOAlgorithm, PGBOAlgorithm
  )
  stems = c(
    "OptimizationMultiplicativeAlgorithm",
    "OptimizationFedorovWynnAlgorithm",
    "OptimizationSimplexAlgorithm",
    "OptimizationPSOAlgorithm",
    "OptimizationPGBOAlgorithm"
  )
  walk2( algos, stems, function( algo, stem ) {
    tpl = paste0( stem, fimSuffix, "FIM.Rmd" )
    method( generateReportOptimization, list( fimClass, algo ) ) =
      .renderReport( tpl )
  } )
}
