## ----global_options, echo = FALSE, include = FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(tangle = FALSE)
backup_options = options()
library(PFIM)
set.seed(42)
options(width = 200)
utils = system.file("vignette-scripts", "pfim-vignette-utils.R", package = "PFIM")
if (!nzchar(utils)) stop("pfim-vignette-utils.R not found.", call. = FALSE)
source(utils, local = knitr::knit_global())
paths = pfimVignetteSetupPaths()
plotOptions = list(unitTime = c("hour"), unitOutcomes = c("mcg/mL", "DI%"))
.pfimVignetteHas = function( name ) {
  exists( name, inherits = TRUE ) && {
    val = get( name, inherits = TRUE )
    !is.null( val ) && ( !is.character( val ) || any( nzchar( val ) ) )
  }
}
knitr::opts_chunk$set(purl = FALSE, collapse = TRUE,
                      comment = "#>", echo = FALSE, warning = FALSE, message = FALSE,
                      cache = FALSE, tidy = FALSE,
                      fig.align = "center", out.width = "62%", dpi = 110,
                      fig.width = 5, fig.height = 4, dev = "png")

