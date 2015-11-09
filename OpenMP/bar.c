#include <stdio.h>
#include <assert.h>
#include <R.h>
#include <Rinternals.h>
#include <omp.h>

SEXP dotProd(SEXP a, SEXP b) 
{
  SEXP c;
  double *xa, *xb, *xc;
  int np, na, nb, i;
  double dot = 0.0;
  
  a = PROTECT(coerceVector(a, REALSXP));
  b = PROTECT(coerceVector(b, REALSXP));
  c = PROTECT(allocVector(REALSXP, 1));
    
  xa = REAL(a); 
  xb = REAL(b); 
  xc = REAL(c);
  
  np = omp_get_num_procs();
  omp_set_num_threads(np);  
  printf("number of processors: %d\n", np);

  na = length(a);
  nb = length(b);
  assert(na == nb);

  #pragma omp parallel for \
          private(i) \
          shared(xa, xb, na) \
          reduction(+:dot)
  for (i = 0; i < na; ++i)
    dot += xa[i] * xb[i];
    
  xc[0] = dot;
  UNPROTECT(3);
  return c;
}
