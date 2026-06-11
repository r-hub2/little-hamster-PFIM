args <- commandArgs(trailingOnly = TRUE)
pkg <- if (length(args)) normalizePath(args[1], winslash = "/") else {
  script_dir <- getwd()
  if (file.exists(file.path(script_dir, "DESCRIPTION"))) {
    script_dir
  } else {
    normalizePath(file.path(script_dir, "..", ".."), winslash = "/")
  }
}
stopifnot(file.exists(file.path(pkg, "DESCRIPTION")))

rcmd_bin <- file.path(R.home("bin"), "Rcmd.exe")
if (!file.exists(rcmd_bin)) rcmd_bin <- file.path(R.home("bin"), "x64", "Rcmd.exe")
objdump_cfg <- system2(rcmd_bin, "config OBJDUMP", stdout = TRUE, stderr = TRUE)
message("Rcmd OBJDUMP: ", paste(objdump_cfg, collapse = " "))
message("which objdump: ", Sys.which("objdump"))
message(
  "src: init.c=", file.exists(file.path(pkg, "src/init.c")),
  ", RcppExports has R_init_PFIM=",
  any(grepl(
    "R_init_PFIM",
    readLines(file.path(pkg, "src/RcppExports.cpp"), warn = FALSE)
  ))
)

tmpdir <- tempfile("pfimchk")
dir.create(tmpdir)
rcmd <- rcmd_bin
inst <- system2(rcmd, c("INSTALL", pkg, "-l", tmpdir), stdout = TRUE, stderr = TRUE)
message("R CMD INSTALL exit status: ", attr(inst, "status"))

dll <- file.path(tmpdir, "PFIM", "libs", "x64", "PFIM.dll")
message("DLL: ", dll, " (exists=", file.exists(dll), ")")

if (file.exists(dll)) {
  syms <- tryCatch(
    suppressWarnings(tools:::read_symbols_from_dll(dll, "x86_64-w64-mingw32")),
    error = function(e) paste0("ERROR:", conditionMessage(e)),
    warning = function(w) {
      invokeRestart("muffleWarning")
    }
  )
  if (is.character(syms) && length(syms) == 1L && grepl("^ERROR:", syms)) {
    message("read_symbols_from_dll failed: ", syms)
  } else {
    message(
      "read_symbols_from_dll: n=", length(syms),
      ", R_init_PFIM=", any(grepl("R_init_PFIM", syms)),
      ", R_registerRoutines=", any(grepl("R_registerRoutines", syms)),
      ", head=", paste(head(syms, 20), collapse = ",")
    )
  }
}

nm_paths <- c(
  "C:/RBuildTools/4.4/x86_64-w64-mingw32.static.posix/bin/nm.exe",
  "C:/rtools44/mingw64/bin/nm.exe",
  "C:/Program Files/Git/usr/bin/nm.exe"
)
for (p in nm_paths) {
  if (file.exists(p) && file.exists(dll)) {
    out <- system2(p, c("-g", dll), stdout = TRUE, stderr = TRUE)
    message(
      "nm (", p, "): R_init_PFIM=", any(grepl("R_init_PFIM", out)),
      ", R_registerRoutines=", any(grepl("R_registerRoutines", out)),
      ", R_useDynamicSymbols=", any(grepl("R_useDynamicSymbols", out))
    )
  }
}
