#include "RngStream.h"
#include <stdio.h>
#include <assert.h>
#include <R.h>
#include <Rinternals.h>
#include <omp.h>

SEXP noise(SEXP a) 
{
  SEXP c;
  double *xa, *xb, *xc;
  int rank, np, na, i;
  
  na = length(a);

  a = PROTECT(coerceVector(a, REALSXP));
  c = PROTECT(allocVector(REALSXP, na));
    
  xa = REAL(a); 
  xc = REAL(c);
  
  np = omp_get_num_procs();
  omp_set_num_threads(np);  
  printf("number of processors: %d\n", np);

  assert(na % np == 0);

  unsigned long seed[6] = {1806547166, 3311292359,  
			   643431772, 1162448557, 
			   3335719306, 4161054083};
  RngStream_SetPackageSeed(seed);
  
  RngStream RngArray[na];
  for (i = 0; i < na; ++i) {
    RngArray[i] = RngStream_CreateStream(NULL);
  }

  #pragma omp parallel default(shared) private(rank, i, xb)
  {
    rank = omp_get_thread_num();
  
    xb = (double*) malloc(na*sizeof(double));

    for (i = 0; i < na/np; ++i) {
      xb[i] = RngStream_RandU01(RngArray[rank*na/np+i]);
    }

    for (i = 0; i < na; ++i) {
      xc[i] = xa[i] + xb[i];
    }

    free(xb);
  }

  UNPROTECT(2);
  return c;
}
