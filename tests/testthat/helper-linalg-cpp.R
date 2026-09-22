# R references for C++ linalg kernels (chol_inv, safe_solve).

.safeCholInv_ref = function( V ) {
  M = 0.5 * ( V + t( V ) )
  chol2inv( chol( M ) )
}

.safeSolve_ref = function( M ) solve( M )

.random_spd = function( p, seed ) {
  set.seed( seed )
  crossprod( matrix( rnorm( p * p ), p, p ) ) + diag( p ) * 1e-6
}

.pop_fim_simple_V = function( evaluation ) {
  model = rebuildEvalModel( evaluation, finiteDifference = FALSE )
  arm   = prop( prop( evaluation, "evaluationDesign" )[[ 1L ]], "evaluationArms" )[[ 1L ]]
  parameters = prop( model, "modelParameters" )
  outputNames = prop( model, "outputNames" )
  allGradientsData = prop( arm, "evaluationGradients" )
  varianceResults  = prop( arm, "evaluationVariance" )
  distributions = map( parameters, ~ prop( .x, "distribution" ) )
  omega_IIV = vapply( distributions, function( x ) prop( x, "omega" )^2, numeric( 1L ) )
  muChain   = PFIM:::.pfimPopMuChainFactors( parameters )
  gradients = do.call( cbind, map( outputNames, ~ t( allGradientsData[[ .x ]] ) ) )
  G = gradients * muChain
  crossprod( diag( omega_IIV ) %*% G, G ) + as.matrix( varianceResults$errorVariance )
}
