# Frozen gold table: do not "fix" numbers to match a new run.

test_that( "PFIM gold table is complete and frozen", {
  expect_true( is.list( .pfimGold ) )
  expect_equal(
    names( .pfimGold$gcsf_design1$rse_mu ),
    c( "KA", "KEL", "VD", "KD", "KINT", "KSI", "KMT", "KTT", "NB0", "SC1", "SM1", "SM2" )
  )
  expect_equal( unname( .pfimGold$gcsf_design1$rse_mu[[ "KA" ]] ), 3.289140 )
  expect_equal( .pfimGold$cas1_noCov_noIOV$D, 731.067615428391 )
  expect_equal( .pfimGold$cas6_occasion_iov$D, 1837.7481985995 )
} )
