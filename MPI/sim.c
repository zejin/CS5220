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

  int nproc, rank, i, j, k;
  double a, b;
  int nobs = 25;
  int ndim = 5;
  int nrep = 16;

  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  extern char* optarg;
  const char* optstring = "o:d:r:";
  int c;
  while ((c = getopt(argc, argv, optstring)) != -1) {
    switch (c) {
    case 'o': nobs = atoi(optarg); break; // number of observations
    case 'd': ndim = atoi(optarg); break; // number of dimensions
    case 'r': nrep = atoi(optarg); break; // number of repetitions
    }
  }

  if (rank == 0) {
    assert(nrep % nproc == 0);
    
    printf("nobs: %d\n", nobs);
    printf("ndim: %d\n", ndim);
    printf("nrep: %d\n", nrep);
    printf("nproc: %d\n", nproc);
  }

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
  double* x = (double*) malloc(nobs*ndim*nrep/nproc*sizeof(double));
  double* y = (double*) malloc(nobs*ndim*nrep/nproc*sizeof(double));
  double* local = (double*) malloc(nrep/nproc*sizeof(double));
  int* index = (int*) malloc(nobs*sizeof(int));
  double* global;

  for (i = 0; i < nobs*ndim*nrep/nproc; ++i) {
    x[i] = RngStream_RandNormal(RngArray[rank]);
  }

  for (i = 0; i < nobs*ndim*nrep/nproc; ++i) {
    y[i] = RngStream_RandNormal(RngArray[rank]);
  }

  for (i = 0; i < nrep/nproc; ++i) {
    for (j = 0; j < nobs; ++j) {
      b = 0.0;
      for (k = 0; k < ndim; ++k) {
        a = x[i*nobs*ndim+j*ndim+k] - y[i*nobs*ndim+j*ndim+k];
        b += a * a;
      }
      b = sqrt(b);
      local[i] += b;
    }
  }

  for (i = 0; i < nobs; ++i) {
    index[i] = i;
  }
  RngStream_RandShuffle(RngArray[rank], index, nobs);

  if (rank == 0) {
    global = (double*) malloc(nrep*sizeof(double));
  } 

  MPI_Gather(local, nrep/nproc, MPI_DOUBLE, 
	     global, nrep/nproc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    for (i = 0; i < nobs; ++i) {
      printf("%d\n", index[i]);
    }
    
    for (i = 0; i < nrep; ++i) {
      printf("%g\n", global[i]);
    }

    free(global);
  }
  
  free(x);
  free(y);
  free(local);
  free(index);

  double t1 = MPI_Wtime();

  if (rank == 0) {
    printf("time: %g\n", t1-t0);
  }

  MPI_Finalize();
  
  return 0;
}
