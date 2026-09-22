# Population mu chain rule (adjustGradient in variance block).

test_that("adjustGradient drives population mu chain factors", {
  mp_ln = ModelParameter( name = "k", distribution = LogNormal( mu = 2, omega = 0.3 ) )
  mp_n  = ModelParameter( name = "k", distribution = Normal( mu = 2, omega = 0.3 ) )
  expect_equal( PFIM:::.pfimPopMuChainFactors( list( mp_ln ) ), 2 )
  expect_equal( PFIM:::.pfimPopMuChainFactors( list( mp_n ) ), 1 )
})

test_that("Normal chain factor does not scale population combo block like LogNormal", {
  G = matrix( c( 1, 0.5 ), nrow = 2L, ncol = 1L )
  base = list(
    gradients         = G,
    omega_iiv         = c( 0.01, 0.01 ),
    gamma_values      = c( 0, 0 ),
    error_variance    = matrix( 0.1, 1L, 1L ),
    occ_col_widths    = 1L,
    sigma_derivatives = list(),
    has_iov           = FALSE
  )
  normal = do.call( .computePopFimCombo_R, c( base, list( mu_values = c( 1, 1 ) ) ) )
  wrong  = do.call( .computePopFimCombo_R, c( base, list( mu_values = c( 10, 10 ) ) ) )
  expect_false( isTRUE( all.equal( normal, wrong, tolerance = 1e-12 ) ) )
})
