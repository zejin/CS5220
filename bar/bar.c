#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <omp.h>

SEXP convolve2(SEXP a, SEXP b) 
{
  int na, nb, nab;
  double *xa, *xb, *xab;
  SEXP ab;

  a = PROTECT(coerceVector(a, REALSXP));
  b = PROTECT(coerceVector(b, REALSXP));
  na = length(a); 
  nb = length(b); 
  nab = na + nb - 1;
  ab = PROTECT(allocVector(REALSXP, nab));
    
  xa = REAL(a); 
  xb = REAL(b); 
  xab = REAL(ab);

  for (int i = 0; i < nab; i++) 
    xab[i] = 0.0;

  int np = omp_get_num_procs();
  omp_set_num_threads(np);//set number of threads   
  printf("number of threads: %d\n", np);

  #pragma omp parallel for
  for (int i = 0; i < na; i++)
    for (int j = 0; j < nb; j++) 
      xab[i + j] += xa[i] * xb[j];
    
  UNPROTECT(3);
  return ab;
}
