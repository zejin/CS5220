#include "RngStream.h"
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

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
  MPI_Init(&argc, &argv);

  double t0 = MPI_Wtime();

  int nproc, rank, i;
  int nobs = 8;
  int nrep = 1;

  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  extern char* optarg;
  const char* optstring = "o:r:";
  int c;
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
  double* global;
  double* local = (double*) malloc(nobs/nproc*sizeof(double));
  int* index = (int*) malloc(nobs/nproc*sizeof(int));
  
  if (rank == 0) {
    global = (double*) malloc(nobs*sizeof(double));
  } 

  for (i = 0; i < nobs/nproc; ++i) {
    local[i] = RngStream_RandNormal(RngArray[rank]);
  }

  for (i = 0; i < nobs/nproc; ++i) {
    index[i] = i;
  }
  RngStream_RandShuffle(RngArray[rank], index, nobs/nproc);

  MPI_Gather(local, nobs/nproc, MPI_DOUBLE, 
	     global, nobs/nproc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    for (i = 0; i < nobs/nproc; ++i) {
      printf("%d\n", index[i]);
    }
    
    for (i = 0; i < nobs; ++i) {
      printf("%g\n", global[i]);
    }

    free(global);
  }
  
  free(local);
  free(index);

  double t1 = MPI_Wtime();

  if (rank == 0) {
    printf("time: %g\n", t1-t0);
  }

  MPI_Finalize();
  
  return 0;
}
