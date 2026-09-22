# G-CSF PopED Design 1 RSE agreement.
# Always runs under R CMD check (~20s linear-FD FIM). Plots stay skip_on_cran.

.gcsfExample04Script = function() {
  script = file.path( "inst", "vignette-scripts", "example04_execute.R" )
  if ( !file.exists( script ) )
    script = system.file( "vignette-scripts", "example04_execute.R", package = "PFIM" )
  script
}

.gcsfRestoreEnv = function( name, old ) {
  if ( is.na( old ) )
    Sys.unsetenv( name )
  else
    do.call( Sys.setenv, stats::setNames( list( old ), name ) )
}

.gcsfSourceExample04 = function( cmp_only ) {
  script = .gcsfExample04Script()
  skip_if_not( nzchar( script ) && file.exists( script ), "example04_execute.R not found" )

  old_force = Sys.getenv( "PFIM_GCSF_FORCE_RUN", unset = NA_character_ )
  old_cmp = Sys.getenv( "PFIM_GCSF_CMP_ONLY", unset = NA_character_ )
  Sys.setenv( PFIM_GCSF_FORCE_RUN = "true" )
  if ( isTRUE( cmp_only ) )
    Sys.setenv( PFIM_GCSF_CMP_ONLY = "true" )
  else
    Sys.unsetenv( "PFIM_GCSF_CMP_ONLY" )
  on.exit( {
    .gcsfRestoreEnv( "PFIM_GCSF_FORCE_RUN", old_force )
    .gcsfRestoreEnv( "PFIM_GCSF_CMP_ONLY", old_cmp )
  }, add = TRUE )

  env = new.env( parent = globalenv() )
  source( script, local = env )
  env
}

test_that( "G-CSF PopED Design 1 RSE matches reference (Example04)", {
  env = .gcsfSourceExample04( cmp_only = TRUE )

  expect_true( is.data.frame( env$cmp_mu ) )
  expect_true( is.data.frame( env$cmp_d ) )
  expect_true( all( is.finite( env$cmp_mu$pfim ) ) )
  expect_equal( env$cmp_mu$pfim, env$cmp_mu$poped, tolerance = 1e-3 )
  expect_equal( env$cmp_d$pfim, env$cmp_d$poped, tolerance = 1e-3 )
  expect_true( all( env$cmp_mu$rel_pct >= 0 ) )
  expect_true( all( env$cmp_d$rel_pct >= 0 ) )
  expect_true( all( env$cmp_mu$rel_pct < 0.05 ) )
  expect_true( all( env$cmp_d$rel_pct < 0.05 ) )
  expect_true( is.data.frame( env$cmp_sigma ) )
  expect_true( all( is.finite( env$cmp_sigma$pfim_rse_sd ) ) )
  # 2 * RSE(SD) ≈ PopED RSE(variance)
  expect_equal(
    env$cmp_sigma$pfim_rse_var_equiv,
    env$cmp_sigma$poped_rse_var,
    tolerance = 1e-3
  )
  expect_true( all( env$cmp_sigma$rel_pct >= 0 ) )
  expect_true( all( env$cmp_sigma$rel_pct < 0.5 ) )
  expect_equal( env$cmp_sigma$pfim_sd, sqrt( env$cmp_sigma$poped_var ), tolerance = 1e-6 )
} )

test_that( "G-CSF PopED overlay is a true ODE (not FIM-cache interpolation)", {
  skip_on_cran()
  env = .gcsfSourceExample04( cmp_only = FALSE )

  expect_true( inherits( env$plotOutcomesEvaluationPopedStyle, "ggplot" ) )
  pk = subset(
    env$plotOutcomesEvaluationPopedStyle$data,
    Model == "Model: PK" & Group == "Group 3"
  )
  expect_gt( nrow( pk ), 100L )
  dpred = diff( pk$pred[ order( pk$time ) ] )
  n_peaks = sum( dpred[ -length( dpred ) ] > 0 & dpred[ -1L ] < 0 )
  expect_gte( n_peaks, 5L )
} )
