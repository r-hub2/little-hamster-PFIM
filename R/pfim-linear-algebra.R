# Cholesky-based inverses with jitter when V is near-singular (population FIM blocks).

.safeCholInv = function( V ) {
  tryCatch(
    chol2inv( chol( V ) ),
    error = function( e ) {
      eps = .Machine$double.eps^0.5 * max( abs( diag( V ) ) )
      chol2inv( chol( V + diag( eps, nrow( V ) ) ) )
    }
  )
}

.safeSolve = function( M ) {
  tryCatch(
    solve( M ),
    error = function( e ) {
      warning(
        "Matrix is singular or ill-conditioned; using Tikhonov-regularized inverse.",
        call. = FALSE
      )
      eps = .Machine$double.eps^0.5 * max( abs( diag( M ) ) )
      solve( M + diag( eps, nrow( M ) ) )
    }
  )
}
