#include "RngStream.h"
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <mpi.h>

int main(int argc, char** argv)
{
  //
  int c, i, rank, nproc;
  int nobs = 10;
  int nrep = 1;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  extern char* optarg;
  const char* optstring = "o:r:";
  while ((c = getopt(argc, argv, optstring)) != -1) {
    switch (c) {
    case 'o': nobs = atoi(optarg); break; // number of observations
    case 'r': nrep = atoi(optarg); break; // number of repetitions
    }
  }

  if (rank == 0) {
    assert(nobs % nproc == 0);
    
    printf("nobs: %d\n", nobs);
    printf("nrep: %d\n", nrep);
    printf("nproc: %d\n", nproc);
  }

  double t0 = MPI_Wtime();

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
  double* local = (double*) malloc(nobs/nproc*sizeof(double));
  double* global = (double*) malloc(nobs*sizeof(double));

  for (i = 0; i < nobs/nproc; ++i) {
    local[i] = RngStream_RandU01(RngArray[rank]);
  }

  MPI_Gather(local, nobs/nproc, MPI_DOUBLE, 
	     global, nobs/nproc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  double t1 = MPI_Wtime();

  if (rank == 0) {
    for (i = 0; i < nobs; ++i) {
      printf("%g\n", global[i]);
    }

    printf("time: %g\n", t1-t0);
  }
  
  free(local);
  free(global);

  MPI_Finalize();
  
  return 0;
}
