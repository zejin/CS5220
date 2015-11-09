#include "RngStream.h"
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <omp.h>

int main(int argc, char** argv)
{
  //
  int c, i, rank;
  int nobs = 10;
  int nrep = 1;
  int nproc = 1;

  extern char* optarg;
  const char* optstring = "o:r:p:";
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

  double t0 = omp_get_wtime();

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

  omp_set_num_threads(nproc);

  #pragma omp parallel private(rank, i) 
  {
    rank = omp_get_thread_num();
    for (i = 0; i < nobs/nproc; ++i) {
      global[rank*nobs/nproc+i] = RngStream_RandU01(RngArray[rank]);
    }  
  }

  for (i = 0; i < nobs; ++i) {
    printf("%g\n", global[i]);
  }

  double t1 = omp_get_wtime();
  
  printf("time: %g\n", t1-t0);

  return 0;
}
