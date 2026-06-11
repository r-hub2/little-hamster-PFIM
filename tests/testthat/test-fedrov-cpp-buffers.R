# Fedorov-Wynn C++ integration: buffer sizing and clean Rcpp calls.

test_that("Fedorov-Wynn buffer sizes follow FIM dimension, not protocol count", {
  buf = getFromNamespace( ".fedorovWynnBufferSizes", "PFIM" )( 4L, 3L )
  expect_equal( buf$nFisher, 10L )
  expect_equal( buf$nMaxPop, 11L )
  expect_equal( buf$nBuf, 11L )
  expect_gt( buf$nFisher, 3L * 4L / 2L ) # old nProtocols-based fisher length was 6
})

test_that("FedorovWynnAlgorithm run() has no C++ subscript warnings", {
  expect_no_warning(
    run( .minimal_discrete_opt( "FedorovWynnAlgorithm", name = "fw_cpp_buffers" ) )
  )
})

test_that("FedorovWynnAlgorithm converges with ndimFim > nProtocols grid", {
  opt = expect_no_warning(
    run( .minimal_discrete_opt( "FedorovWynnAlgorithm", name = "fw_cpp_conv" ) )
  )
  od = prop( opt, "optimisationDesign" )
  expect_s7_class( od$evaluationOptimalDesign, Evaluation )
  expect_gt( getDeterminant( od$evaluationOptimalDesign ), 0 )
  expect_gt( getDcriterion( opt ), 0 )
})
