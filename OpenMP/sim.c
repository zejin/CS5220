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

int main(int argc, char** argv)
{
  //
  double t0 = omp_get_wtime();

  int rank, i;
  int nobs = 10;
  int nrep = 1;
  int nproc = 1;

  extern char* optarg;
  const char* optstring = "o:r:p:";
  int c;
  while ((c = getopt(argc, argv, optstring)) != -1) {
    switch (c) {
    case 'o': nobs = atoi(optarg); break; // number of observations
    case 'r': nrep = atoi(optarg); break; // number of repetitions
    case 'p': nproc = atoi(optarg); break; // number of processors
    }
  }

  assert(nobs % nproc == 0);

  printf("nobs: %d\n", nobs);
  printf("nrep: %d\n", nrep);
  printf("nproc: %d\n", nproc);

  //  
  unsigned long seed[6] = {1806547166, 3311292359,  
			   643431772, 1162448557, 
			   3335719306, 4161054083};
  RngStream_SetPackageSeed(seed);
  
  RngStream RngArray[nproc];
  for (i = 0; i < nproc; ++i) {
    RngArray[i] = RngStream_CreateStream(NULL);
  }

  //
  double* global = (double*) malloc(nobs*sizeof(double));
  double* local;
  int* index;

  omp_set_num_threads(nproc);

#pragma omp parallel default(shared) private(rank, i, index, local) 
  {
    rank = omp_get_thread_num();
    local = (double*) malloc(nobs/nproc*sizeof(double));
    index = (int*) malloc(nobs/nproc*sizeof(int));

    for (i = 0; i < nobs/nproc; ++i) {
      global[rank*nobs/nproc+i] = RngStream_RandNormal(RngArray[rank]);
    }

    for (i = 0; i < nobs/nproc; ++i) {
      index[i] = i;
    }
    RngStream_RandShuffle(RngArray[rank], index, nobs/nproc); 

    if (rank == 0) {
      for (i = 0; i < nobs/nproc; ++i) {
        printf("%d\n", index[i]);
      }
    }

    free(local);
    free(index);
  }

  for (i = 0; i < nobs; ++i) {
    printf("%g\n", global[i]);
  }

  free(global);

  double t1 = omp_get_wtime();
  
  printf("time: %g\n", t1-t0);

  return 0;
}
