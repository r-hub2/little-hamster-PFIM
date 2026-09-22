# Evaluation_tests / Optimisation_tests: D-criterion and FIM vs baked references.

purrr::walk( eval_opt_reference_cases(), function( case ) {
  local( {
    ref_case = case
    test_that( paste0( "eval-opt reference: ", ref_case$id ), {
      if ( isTRUE( ref_case$slow ) ) skip_on_cran()

      if ( isTRUE( ref_case$expect_converge_warning ) ) {
        # testthat 3e: expect_warning() returns the condition, not run()'s value.
        prom = testthat::evaluate_promise( ref_case$run() )
        expect_true(
          any( grepl( "converge", prom$warnings ) ),
          info = paste( prom$warnings, collapse = "\n" )
        )
        obj = prom$result
      } else {
        obj = ref_case$run()
      }
      stats = ref_case$summarize( obj )
      expect_fim_summary_matches_reference( stats, ref_case$reference, ref_case$tol, label = ref_case$id )

      M = stats$fisherMatrix
      expect_true( is.matrix( M ) && isSymmetric( M, tol = 1e-8 ) )
      expect_gt( det( M ), 0 )
      expect_equal( nrow( M ), ncol( M ) )

      if ( !is.null( ref_case$ref_html ) && !is.na( ref_case$ref_html ) && file.exists( ref_case$ref_html ) ) {
        html_parts = Filter( Negate( is.null ), parse_report_fim_summaries( ref_case$ref_html ) )
        idx = ref_case$html_index %||% 1L
        if ( length( html_parts ) >= idx ) {
          expect_equal(
            stats$d, html_parts[[ idx ]]$d,
            tolerance = ( ref_case$tol$d %||% 1e-3 ) + 1e-6,
            info = paste0( ref_case$id, ": D vs ", basename( ref_case$ref_html ) )
          )
        }
      }
    } )
  } )
})

test_that( "parse_report_fim_summaries reads mult opt popFIM.html", {
  ref_root = eval_opt_reference_root()
  skip_if( is.na( ref_root ), "tests_PFIM/resultats_de_references not found" )
  html = file.path(
    ref_root,
    "Optimisation_tests/discrete/MultiplicativeAlgorithm",
    "multiplicative_Algorithm_PK_ode_dose_not_in_eqs_results/popFIM.html"
  )
  skip_if_not( file.exists( html ), "reference mult opt HTML missing" )

  parts = Filter( Negate( is.null ), parse_report_fim_summaries( html ) )
  expect_gte( length( parts ), 2L )
  expect_equal( parts[[ 1L ]]$d, 6810.422, tolerance = 1e-3 )
} )
