#include "RngStream.h"
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <math.h>
#include <omp.h>

static const double PI = 3.14159265358979323846;

//
double RngStream_RandNormal(RngStream g) 
{
  double u1 = RngStream_RandU01(g);
  double u2 = RngStream_RandU01(g);

  return sqrt(-2*log(u1)) * cos(2*PI*u2);
}

//
void RngStream_RandShuffle(RngStream g, int* index, int n) 
{
  int i, j, temp;
  for (int i = n - 1; i > 0; --i) {
    j = RngStream_RandInt(g, 0, i);

    temp = index[i];
    index[i] = index[j];
    index[j] = temp;
  }
}

//
int main(int argc, char** argv)
{
  //
  double t0 = omp_get_wtime();

  int rank, i, j, n, m;
  double a, b;
  int nobs = 25;
  int ndim = 5;
  int nrep = 16;
  int nproc = 8;

  extern char* optarg;
  const char* optstring = "o:d:r:p:";
  int c;
  while ((c = getopt(argc, argv, optstring)) != -1) {
    switch (c) {
    case 'o': nobs = atoi(optarg); break; // number of observations
    case 'd': ndim = atoi(optarg); break; // number of dimensions
    case 'r': nrep = atoi(optarg); break; // number of repetitions
    case 'p': nproc = atoi(optarg); break; // number of processors
    }
  }

  assert(nrep % nproc == 0);

  printf("nobs: %d\n", nobs);
  printf("ndim: %d\n", ndim);
  printf("nrep: %d\n", nrep);
  printf("nproc: %d\n", nproc);

  //  
  unsigned long seed[6] = {1806547166, 3311292359,  
			   643431772, 1162448557, 
			   3335719306, 4161054083};
  RngStream_SetPackageSeed(seed);
  
  RngStream RngArray[nrep];
  for (i = 0; i < nrep; ++i) {
    RngArray[i] = RngStream_CreateStream(NULL);
  }

  //
  double* x;
  double* y;
  double* local;
  int* index;
  double* global = (double*) malloc(nrep*sizeof(double));

  omp_set_num_threads(nproc);

  #pragma omp parallel default(shared) \
  private(rank, i, j, n, m, a, b, x, y, local, index) 
  {
    rank = omp_get_thread_num();

    x = (double*) malloc(nobs*ndim*nrep/nproc*sizeof(double));
    y = (double*) malloc(nobs*ndim*nrep/nproc*sizeof(double));
    local = (double*) malloc(nrep/nproc*sizeof(double));
    index = (int*) malloc(nobs*sizeof(int));

    for (n = 0; n < nrep/nproc; ++n) {

      for (i = 0; i < nobs*ndim; ++i) {
        x[i] = RngStream_RandNormal(RngArray[rank*nrep/nproc+n]);
      }

      for (i = 0; i < nobs*ndim; ++i) {
        y[i] = RngStream_RandNormal(RngArray[rank*nrep/nproc+n]);
      }

      for (i = 0; i < nobs; ++i) {
        b = 0.0;
        for (j = 0; j < ndim; ++j) {
          a = x[i*ndim+j] - y[i*ndim+j];
          b += a * a;
        }
        b = sqrt(b);
        local[n] += b;
      }

      for (i = 0; i < nobs; ++i) {
        index[i] = i;
      }
      RngStream_RandShuffle(RngArray[rank*nrep/nproc+n], index, nobs);

      //







    }

    for (i = 0; i < nrep/nproc; ++i) {
      global[rank*nrep/nproc+i] = local[i];
    }

    if (rank == 0) {
      for (i = 0; i < nobs; ++i) {
        printf("%d\n", index[i]);
      }
    }

    free(x);
    free(y);
    free(local);
    free(index);
  }

  for (i = 0; i < nrep; ++i) {
    printf("%g\n", global[i]);
  }

  free(global);

  double t1 = omp_get_wtime();
  
  printf("time: %g\n", t1-t0);

  return 0;
}
