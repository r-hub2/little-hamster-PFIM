# G-CSF / filgrastim PK-PD population FIM (PopED K/PD example).
# Canonical script: inst/vignette-scripts/example04_execute.R (vignette Example04).
#
# Usage:
#   source(system.file("examples", "gcsf-poped-evaluation.R", package = "PFIM"))

script = system.file("vignette-scripts", "example04_execute.R", package = "PFIM")
if (!nzchar(script) || !file.exists(script))
  script = file.path("inst", "vignette-scripts", "example04_execute.R")
if (!file.exists(script))
  stop("example04_execute.R not found.", call. = FALSE)

env = new.env(parent = globalenv())
source(script, local = env)
invisible(list(
  evaluation = env$evaluationPop,
  rse = env$rse,
  cmp_mu = env$cmp_mu,
  cmp_d = env$cmp_d
))
