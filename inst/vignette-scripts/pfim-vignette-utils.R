# Shared helpers for Example01--04 vignette execute scripts.

pfimEnsurePandoc = function() {
  if (requireNamespace("rmarkdown", quietly = TRUE) &&
      rmarkdown::pandoc_available()) return(invisible(TRUE))
  candidates = unique(normalizePath(c(
    Sys.getenv("RSTUDIO_PANDOC", unset = ""),
    file.path(Sys.getenv("ProgramFiles", unset = "C:/Program Files"),
              "RStudio", "resources", "app", "bin", "quarto", "bin", "tools"),
    file.path(Sys.getenv("ProgramFiles", unset = "C:/Program Files"), "Pandoc")
  ), winslash = "/", mustWork = FALSE))
  for (dir in candidates) {
    exe = file.path(dir, if (.Platform$OS.type == "windows") "pandoc.exe" else "pandoc")
    if (file.exists(exe)) {
      Sys.setenv(RSTUDIO_PANDOC = dir)
      if (requireNamespace("rmarkdown", quietly = TRUE) &&
          rmarkdown::pandoc_available()) return(invisible(TRUE))
    }
  }
  invisible(FALSE)
}

pfimSafeReport = function(...) {
  if (!pfimEnsurePandoc()) return(invisible(NULL))
  tryCatch(invisible(Report(...)), error = function(e) {
    message("[PFIM vignette] Report() skipped: ", conditionMessage(e))
    invisible(NULL)
  })
}

pfimCapture = function(expr) paste(capture.output(expr), collapse = "\n")

pfimSafePlot = function(expr) {
  tryCatch(expr, error = function(e) {
    message("[PFIM vignette] plot skipped: ", conditionMessage(e))
    NULL
  })
}

pfimFinalizeOpt = function(opt, paths, reportName, showName, reportOpts) {
  show_file = file.path(paths$outputs, showName)
  out = pfimOptShow(opt, show_file)
  writeLines(out, show_file)
  invisible(pfimSafeReport(opt, paths$reports, reportName, reportOpts))
  attr(opt, "showOutput") = out
  opt
}

pfimOptShow = function(opt, show_file) {
  tryCatch(
    pfimCapture({
      show(opt)
      getFisherMatrix(opt)
      getCorrelationMatrix(opt)
      getSE(opt)
      getRSE(opt)
      getShrinkage(opt)
      getDeterminant(opt)
      getDcriterion(opt)
    }),
    error = function(e) {
      if (file.exists(show_file)) {
        paste(readLines(show_file, warn = FALSE), collapse = "\n")
      } else {
        paste0(
          "[Loaded optimization — show() output unavailable for this RDS object.\n",
          "Run run() and saveRDS() with the current PFIM version to regenerate.\n",
          conditionMessage(e)
        )
      }
    }
  )
}

# Load a shipped Optimization RDS. During R CMD check, missing caches stop;
# rebuild locally only when NOT_CRAN=true and the file is absent.
pfimRunOrLoad = function(rds_path, runner) {
  if (file.exists(rds_path)) {
    message("[PFIM vignette] Loading ", basename(rds_path))
    return(readRDS(rds_path))
  }
  under_check = nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")) ||
    !identical(Sys.getenv("NOT_CRAN"), "true")
  if (under_check) {
    stop(
      "Missing vignette cache '", basename(rds_path),
      "' under vignettes/data/. Ship the RDS with the package; ",
      "do not rebuild optimizations during R CMD check.",
      call. = FALSE
    )
  }
  message("[PFIM vignette] Running optimization -> ", basename(rds_path))
  dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)
  opt = runner()
  saveRDS(opt, rds_path)
  message("[PFIM vignette] Saved ", basename(rds_path))
  opt
}

# Vignette HTML / figures / text dumps. Default is tempdir() so R CMD check
# never writes into the package or check tree (CRAN: only the session temp dir).
# Set PFIM_VIGNETTE_RESULTS to an absolute path to keep a local results/ copy.
pfimResultsRoot = function() {
  override = Sys.getenv("PFIM_VIGNETTE_RESULTS", unset = "")
  if (nzchar(override)) return(override)
  file.path(tempdir(), "PFIM_vignette_results")
}

pfimVignettePaths = function() {
  root = pfimResultsRoot()
  data_candidates = unique(c(
    file.path(getwd(), "data"),
    file.path(getwd(), "vignettes", "data"),
    system.file("vignette-data", package = "PFIM"),
    file.path(getwd(), "..", "vignettes", "data"),
    file.path(getwd(), "..", "..", "vignettes", "data")
  ))
  hit = Find(
    function(d) {
      nzchar(d) && dir.exists(d) &&
        length(list.files(d, pattern = "[.]RDS$")) > 0L
    },
    data_candidates
  )
  list(
    data    = if (is.null(hit)) data_candidates[[1L]] else
      normalizePath(hit, winslash = "/", mustWork = FALSE),
    figures = file.path(root, "figures"),
    reports = root,
    outputs = file.path(root, "outputs")
  )
}

# Create figure / report / output dirs under pfimResultsRoot(). Do not mkdir
# the data cache (shipped, read-only).
pfimVignetteSetupPaths = function() {
  paths = pfimVignettePaths()
  lapply(
    paths[c("figures", "reports", "outputs")],
    dir.create, recursive = TRUE, showWarnings = FALSE
  )
  paths
}
