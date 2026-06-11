// Fedorov-Wynn exchange on discrete protocol weights (Fedorov 1972).
// fisher buffers hold ndim*(ndim+1)/2 packed elements — sized from R via ndimFim.
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

using namespace Rcpp;

#define CHECK_ALLOC(ptr, label) do { \
  if ( (ptr) == NULL ) Rcpp::stop("Fedorov-Wynn: memory allocation failed (%s)", label); \
} while (0)

// ── Tolerances / constants ──────────────────────────────────────────────────
#define FTOL    1e-6
#define FREQMIN 1e-8
#define EPSILON DBL_EPSILON
#define CRI_MAX 1e30
#define TINY    1.0e-20
#define SMALL   1.0e-7

// ── Data structures ─────────────────────────────────────────────────────────

typedef struct PROTOC {
  int     ntps;
  double *tps;
  int     ndim;
  double *fisher;
} PROTOC;

typedef struct POPROT {
  int     maxnp;
  int     np;
  int     ndim;
  int    *num;
  PROTOC *pind;
  double *freq;
  double *fisher;
  double *finv;
  double  det;
} POPROT;

typedef struct matrix_t {
  int     nrow;
  int     ncol;
  double **m;
} matrix_t;

// ── Allocation helpers ───────────────────────────────────────────────────────

static void PROTOC_alloc( PROTOC *p, int ntps, int ndim ) {
  p->ntps   = ntps;
  p->tps    = (double *)calloc( p->ntps, sizeof(double) );
  CHECK_ALLOC( p->tps, "PROTOC.tps" );
  p->ndim   = ndim;
  p->fisher = (double *)calloc( p->ndim*(p->ndim+1)/2, sizeof(double) );
  CHECK_ALLOC( p->fisher, "PROTOC.fisher" );
}

static POPROT *POPROT_alloc( int np, int ndim, int maxnp ) {
  POPROT *p = (POPROT *)malloc( sizeof(POPROT) );
  CHECK_ALLOC( p, "POPROT" );
  p->maxnp  = maxnp;
  p->np     = np;
  p->ndim   = ndim;
  p->num    = (int    *)calloc( p->maxnp, sizeof(int)    );  CHECK_ALLOC( p->num,    "POPROT.num"    );
  p->pind   = (PROTOC *)calloc( p->maxnp, sizeof(PROTOC) );  CHECK_ALLOC( p->pind,   "POPROT.pind"   );
  p->freq   = (double *)calloc( p->maxnp, sizeof(double) );  CHECK_ALLOC( p->freq,   "POPROT.freq"   );
  p->fisher = (double *)calloc( p->ndim*(p->ndim+1)/2, sizeof(double) );  CHECK_ALLOC( p->fisher, "POPROT.fisher" );
  p->finv   = (double *)calloc( p->ndim*(p->ndim+1)/2, sizeof(double) );  CHECK_ALLOC( p->finv,   "POPROT.finv"   );
  p->det    = 1.0;
  return p;
}

static void PROTOC_copy( PROTOC *p, PROTOC *p1 ) {
  int ndat = p->ndim*(p->ndim+1)/2;
  for ( int i = 0; i < p->ntps; ++i ) p1->tps[i]    = p->tps[i];
  for ( int i = 0; i < ndat;    ++i ) p1->fisher[i]  = p->fisher[i];
}

static matrix_t *matrix_create( int nrow, int ncol ) {
  matrix_t *mat = (matrix_t *)malloc( sizeof(matrix_t) );
  CHECK_ALLOC( mat, "matrix_t" );
  mat->nrow = nrow;
  mat->ncol = ncol;
  mat->m    = (double **)calloc( nrow, sizeof(double *) );
  CHECK_ALLOC( mat->m, "matrix_t.m" );
  for ( int i = 0; i < nrow; ++i ) {
    mat->m[i] = (double *)calloc( ncol, sizeof(double) );
    CHECK_ALLOC( mat->m[i], "matrix_t.row" );
  }
  return mat;
}

static void PROTOC_destroy( PROTOC *p ) { free(p->fisher); free(p->tps); }
static void POPROT_destroy( POPROT *p ) {
  free(p->finv); free(p->fisher); free(p->freq);
  free(p->pind); free(p->num); free(p);
}
static void matrix_destroy( matrix_t *mat ) {
  for ( int i = 0; i < mat->nrow; ++i ) free( mat->m[i] );
  free( mat->m ); free( mat );
}

// ── Numerical helpers ────────────────────────────────────────────────────────

static void lubksb( double **a, int n, int *indx, double b[] ) {
  int ii = -1;
  for ( int i = 0; i < n; ++i ) {
    int    ip  = indx[i];
    double sum = b[ip];
    b[ip] = b[i];
    if ( ii >= 0 )
      for ( int j = ii; j <= i-1; ++j ) sum -= a[i][j]*b[j];
    else if ( sum ) ii = i;
    b[i] = sum;
  }
  for ( int i = n-1; i >= 0; --i ) {
    double sum = b[i];
    for ( int j = i+1; j < n; ++j ) sum -= a[i][j]*b[j];
    b[i] = sum / a[i][i];
  }
}

static int ludcmp( double **a, int n, int *indx, double *d ) {
  double *vv = (double *)calloc( n, sizeof(double) );
  *d = 1.0;
  for ( int i = 0; i < n; ++i ) {
    double big = 0.0;
    for ( int j = 0; j < n; ++j ) {
      double temp = fabs( a[i][j] );
      if ( temp > big ) big = temp;
    }
    if ( big == 0.0 ) { free(vv); return 1; }
    vv[i] = 1.0 / big;
  }
  int imax = 0;
  for ( int j = 0; j < n; ++j ) {
    for ( int i = 0; i < j; ++i ) {
      double sum = a[i][j];
      for ( int k = 0; k < i; ++k ) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }
    double big = 0.0;
    for ( int i = j; i < n; ++i ) {
      double sum = a[i][j];
      for ( int k = 0; k < j; ++k ) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      double dum = vv[i]*fabs(sum);
      if ( dum >= big ) { big = dum; imax = i; }
    }
    if ( j != imax ) {
      for ( int k = 0; k < n; ++k ) {
        double dum   = a[imax][k];
        a[imax][k]   = a[j][k];
        a[j][k]      = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if ( a[j][j] == 0.0 ) a[j][j] = TINY;
    if ( j != n-1 ) {
      double dum = 1.0 / a[j][j];
      for ( int i = j+1; i < n; ++i ) a[i][j] *= dum;
    }
  }
  free( vv );
  return 0;
}

// ── FIM helpers ──────────────────────────────────────────────────────────────

static double lik( POPROT *pop ) {
  int    ndim  = pop->ndim;
  int    ncase = ndim*(ndim+1)/2;
  int   *indx  = (int    *)calloc( ndim, sizeof(int)    );
  double *col  = (double *)calloc( ndim, sizeof(double) );
  matrix_t *xa = matrix_create( ndim, ndim );
  double cri;

  for ( int i = 0; i < ncase; ++i ) {
    pop->fisher[i] = 0.0;
    for ( int j = 0; j < pop->np; ++j )
      pop->fisher[i] += pop->freq[j] * pop->pind[j].fisher[i];
  }

  int jj = 0;
  for ( int i = 0; i < ndim; ++i )
    for ( int j = 0; j <= i; ++j ) {
      xa->m[i][j] = pop->fisher[jj];
      xa->m[j][i] = pop->fisher[jj];
      ++jj;
    }

  int ifail = ludcmp( xa->m, ndim, indx, &cri );
  if ( ifail == 1 ) { matrix_destroy(xa); free(indx); free(col); return CRI_MAX; }

  for ( int i = 0; i < ndim; ++i ) cri *= xa->m[i][i];
  pop->det = cri;

  for ( int i = 0; i < ndim; ++i ) {
    for ( int j = 0; j < ndim; ++j ) col[j] = 0.0;
    col[i] = 1.0;
    lubksb( xa->m, ndim, indx, col );
    for ( int j = 0; j < ndim; ++j )
      pop->finv[ i*(i+1)/2 + j ] = col[j];   // upper triangle only (i>=j)
  }

  matrix_destroy(xa); free(indx); free(col);
  return cri;
}

static double matrace( int iprot, PROTOC *prot, POPROT *pop ) {
  double xcal = 0.0;
  for ( int i = 0; i < prot->ndim; ++i ) {
    for ( int j = 0; j < i; ++j ) {
      int ij = i*(i+1)/2 + j;
      xcal  += 2.0 * prot[iprot].fisher[ij] * pop->finv[ij];
    }
    int ij = i*(i+3)/2;
    xcal += prot[iprot].fisher[ij] * pop->finv[ij];
  }
  return xcal;
}

// ── Optimisation ─────────────────────────────────────────────────────────────

static void tassement( POPROT *pop ) {
  int idec = 0, np = pop->np;
  for ( int i = 0; i < pop->np; ++i ) {
    if ( pop->freq[i] < FREQMIN ) {
      ++idec; --np;
    } else if ( idec > 0 ) {
      pop->freq[ i-idec ] = pop->freq[i];
      pop->num[  i-idec ] = pop->num[i];
      PROTOC_destroy( &pop->pind[i-idec] );
      PROTOC_alloc(   &pop->pind[i-idec], pop->pind[i].ntps, pop->pind[i].ndim );
      PROTOC_copy(    &pop->pind[i],      &pop->pind[i-idec] );
    }
  }
  for ( int i = np; i < pop->np; ++i ) PROTOC_destroy( &pop->pind[i] );
  pop->np = np;
}

static int ajout( PROTOC *allprot, POPROT *pop, int nprot ) {
  double xtest = (double)pop->ndim, xdim = xtest;
  int    qajout = -1;

  for ( int iprot = 0; iprot < nprot; ++iprot ) {
    int deja = 0;
    for ( int j = 0; j < pop->np; ++j )
      if ( iprot == pop->num[j] ) deja = 1;
    if ( fabs(deja) < EPSILON ) {
      double xmul = matrace( iprot, allprot, pop );
      if ( xmul >= xtest ) { xtest = xmul; qajout = iprot; }
    }
  }

  if ( fabs(xtest - xdim) < EPSILON ) return 1;   // converged

  double rcalc = (xtest - xdim) / (xdim * (xtest - 1.0));
  PROTOC_alloc( &pop->pind[ pop->np ], allprot[qajout].ntps, allprot[qajout].ndim );
  PROTOC_copy(  &allprot[qajout], &pop->pind[ pop->np ] );
  pop->num[  pop->np ] = qajout;
  pop->freq[ pop->np ] = rcalc;
  for ( int j = 0; j < pop->np; ++j ) pop->freq[j] *= (1.0 - rcalc);
  ++pop->np;

  lik( pop );
  return 0;
}

static int takeout_k( POPROT *pop, PROTOC *allprot, int *nnul ) {
  lik( pop );
  for ( int j = 0; j < pop->np; ++j ) {
    if ( pop->freq[j] <= FREQMIN ) {
      double xmul = matrace( pop->num[j], allprot, pop );
      if ( xmul > (double)pop->ndim ) { *nnul -= 1; return j; }
    }
  }
  return -1;
}

static int project_grad( POPROT *pop, PROTOC *allprot, double *gal, int ir,
                         double *vmgal, int nnul ) {
  int in = -1;
  double sg = 0.0;
  lik( pop );
  sg = 0.0;
  for ( int j = 0; j < pop->np; ++j ) {
    gal[j] = 0.0;
    if ( pop->freq[j] >= FREQMIN || j == ir )
      gal[j] = matrace( pop->num[j], allprot, pop );
    sg += gal[j];
  }
  *vmgal = 0.0;
  for ( int j = 0; j < pop->np; ++j ) {
    if ( pop->freq[j] >= FREQMIN || j == ir ) {
      gal[j] -= sg / (double)(pop->np - nnul);
      if ( fabs(gal[j]) > *vmgal ) *vmgal = fabs(gal[j]);
      if ( gal[j] < 0.0 ) in = j;
    }
  }
  return in;
}

static int optim_lambda( POPROT *pop, double *gal, double *alkeep, double *crit,
                         double dmax, double crit2, double vmgal ) {
  double ro = 1.0, crit1;
  do {
    ro *= 0.5;
    for ( int j = 0; j < pop->np; ++j )
      pop->freq[j] = alkeep[j] + ro*dmax*gal[j];
    crit1 = lik( pop );
    if ( crit1 > *crit ) {
      for ( int j = 0; j < pop->np; ++j ) alkeep[j] = pop->freq[j];
      *crit = crit1;
      return ( fabs((crit2 - *crit)/crit2) > FTOL ) ? 0 : 1;
    }
  } while ( dmax*ro*vmgal > 1e-4 );
  for ( int j = 0; j < pop->np; ++j ) pop->freq[j] = alkeep[j];
  return 1;
}

static int doptimal( POPROT *pop, PROTOC *allprot ) {
  int     nnul = 0, ir = -1, in = -1;
  double  crit, crit1, crit2, vmgal, dmax;
  double *gal    = (double *)calloc( pop->np, sizeof(double) );
  double *alkeep = (double *)calloc( pop->np, sizeof(double) );
  crit = lik( pop );

  for ( int nit = 0; nit < 500; ++nit ) {
    for ( int j = 0; j < pop->np; ++j ) alkeep[j] = pop->freq[j];
    crit2 = crit;

    if ( ir == -1 ) {
      nnul = 0;
      for ( int j = 0; j < pop->np; ++j )
        if ( pop->freq[j] < FREQMIN ) ++nnul;
    }

    in = project_grad( pop, allprot, gal, ir, &vmgal, nnul );
    if ( in < 0 || gal[in] == 0.0 ) { free(gal); free(alkeep); return nnul; }

    dmax = -pop->freq[in] / gal[in];
    int imax = in;
    if ( ir != -1 && gal[ir] <= 0.0 ) { free(gal); free(alkeep); return nnul; }
    ir = -1;
    for ( int j = 0; j < pop->np; ++j ) {
      if ( gal[j] < 0.0 && fabs(pop->freq[j]/gal[j]) < dmax ) {
        dmax = -pop->freq[j] / gal[j]; imax = j;
      }
    }

    double sal = 0.0;
    for ( int j = 0; j < pop->np; ++j ) {
      pop->freq[j] += dmax * gal[j]; sal += pop->freq[j];
    }
    pop->freq[imax] = 0.0;
    crit1 = lik( pop );

    if ( crit1 > crit ) {
      for ( int j = 0; j < pop->np; ++j ) alkeep[j] = pop->freq[j];
      crit = crit1; ++nnul;
    } else {
      for ( int j = 0; j < pop->np; ++j ) pop->freq[j] = alkeep[j];
      int iopt = optim_lambda( pop, gal, alkeep, &crit, dmax, crit2, vmgal );
      if ( iopt == 1 ) {
        ir = takeout_k( pop, allprot, &nnul );
        if ( ir == -1 ) { free(gal); free(alkeep); return nnul; }
      }
    }
  }
  free(gal); free(alkeep);
  return -10;
}

static POPROT *initprot( PROTOC *allprot, IntegerVector protdep,
                         NumericVector freqdep ) {
  int np   = protdep[0];
  int nmax = allprot[0].ndim*(allprot[0].ndim+1)/2 + 1;
  POPROT *mypop = POPROT_alloc( np, allprot[0].ndim, nmax );
  for ( int i = 0; i < np; ++i ) {
    mypop->freq[i] = freqdep[i];
    mypop->num[i]  = protdep[i+1] - 1;
    int ij = mypop->num[i];
    PROTOC_alloc( &mypop->pind[i], allprot[ij].ntps, allprot[ij].ndim );
    PROTOC_copy(  &allprot[ij],    &mypop->pind[i] );
  }
  return mypop;
}

static int main_opt( PROTOC *allprot, POPROT *pop, IntegerVector nprot ) {
  double cri = lik( pop );
  for ( int nit = 0; nit < 100; ++nit ) {
    int nstop = 0, ifin = 0;
    do {
      ifin = ajout( allprot, pop, nprot[0] ); ++nstop;
    } while ( cri < -1e10 && nstop < 100 );
    if ( nstop > 99 && cri < -1e10 ) return -10;
    if ( ifin == 1 ) break;
    int ntest = doptimal( pop, allprot );
    if ( ntest < 0 ) return -100;
    tassement( pop );
    if ( pop->np > pop->ndim*(pop->ndim+1)/2 ) return -1;
    cri = lik( pop );
  }
  return 0;
}

// ── Exported Rcpp entry point ────────────────────────────────────────────────
// P1-02/03 fixes applied:
//   - Removed global `seed` variable (was long seed=-10 at file scope).
//     The RNG (ran1/gauss1) is not used in the D-optimality loop; the functions
//     have been removed.  If stochastic initialisation is ever needed, use
//     R's RNG via GetRNGstate()/PutRNGstate().
//   - Removed all static local variables from helper functions (gauss1, ran1
//     were dead code — not called from main_opt or doptimal).

// [[Rcpp::export]]
Rcpp::List FedorovWynnAlgorithm_Rcpp(
    Rcpp::List    protocols,
    IntegerVector ndimen,
    IntegerVector nbprot,
    IntegerVector numprot,
    NumericVector freq,
    IntegerVector nbdata,
    NumericVector vectps,
    NumericVector fisher,
    IntegerVector error,
    IntegerVector protdep,
    NumericVector freqdep,
    bool          show_process = false )
{
  IntegerVector nb_protocols   = Rcpp::as<IntegerVector>( protocols[0] );
  IntegerVector nb_times       = Rcpp::as<IntegerVector>( protocols[1] );
  IntegerVector nb_dimensions  = Rcpp::as<IntegerVector>( protocols[2] );
  IntegerVector total_cost     = Rcpp::as<IntegerVector>( ndimen );
  NumericMatrix samplingTimes  = Rcpp::as<NumericMatrix>( protocols[4] );
  NumericMatrix fisher_matrices = Rcpp::as<NumericMatrix>( protocols[5] );

  IntegerVector nprot = nb_protocols;
  IntegerVector ndim  = nb_dimensions;
  int Ntot            = total_cost[2];

  NumericMatrix optimal_sampling_times( nprot[0], nb_times[0] );

  PROTOC *allprot = (PROTOC *)calloc( nprot[0], sizeof(PROTOC) );

  for ( int iprot = 0; iprot < nprot[0]; ++iprot ) {
    PROTOC_alloc( &allprot[iprot], nb_times[0], ndim[0] );
    NumericVector row_st = samplingTimes.row( iprot );
    NumericVector row_fm = fisher_matrices.row( iprot );
    for ( int j = 0; j < allprot[iprot].ntps; ++j )
      allprot[iprot].tps[j] = row_st[j];
    int icas = 0;
    for ( int ij = 0; ij < ndim[0]; ++ij )
      for ( int ik = 0; ik <= ij; ++ik ) {
        allprot[iprot].fisher[icas] = row_fm[icas] * Ntot / nb_times[0];
        ++icas;
      }
  }

  POPROT *pop = initprot( allprot, protdep, freqdep );
  int nok     = main_opt( allprot, pop, nprot );
  (void) show_process;
  double sumf = 0.0;

  if ( nok >= 0 ) {
    for ( int i = 0; i < pop->np; ++i ) {
      numprot[i] = pop->num[i] + 1;
      nbdata[i]  = pop->pind[i].ntps;
      freq[i]    = pop->freq[i] * Ntot / nb_times[0];
      sumf       += freq[i];
      freqdep[i]  = pop->freq[i];
    }
    for ( int i = 0; i < pop->np; ++i ) {
      freq[i] /= sumf;
      for ( int j = 0; j < pop->pind[i].ntps; ++j ) {
        vectps[j] = pop->pind[i].tps[j];
        optimal_sampling_times( i, j ) = pop->pind[i].tps[j];
      }
    }
    for ( int i = 0; i < pop->ndim*(pop->ndim+1)/2; ++i )
      fisher[i] = pop->fisher[i];
  }

  POPROT_destroy( pop );
  for ( int i = 0; i < nprot[0]; ++i ) PROTOC_destroy( &allprot[i] );
  free( allprot );

  return Rcpp::List::create(
    Rcpp::Named("freq")                    = freq,
    Rcpp::Named("optimal_sampling_times")  = optimal_sampling_times,
    Rcpp::Named("fisher")                  = fisher,
    Rcpp::Named("numprot")                 = numprot
  );
}
