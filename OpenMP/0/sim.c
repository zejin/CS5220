#include "RngStream.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

// X is a pxn matrix, XX is a nxn matrix
void Double_Center(int n, int p, double *X, double *XX) {
  double* row_sum = (double*) calloc(n, sizeof(double));
  double* col_sum = (double*) calloc(n, sizeof(double));
  
  double total_sum = 0.0;
  double elem, part_sum;
  int i, j, k;

  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      if (i != j) {
        part_sum = 0.0;

        //XX[i, j] = |X[i, ] - X[j, ]|
        for (k = 0; k < p; ++k) {
          elem = X[i*p+k] - X[j*p+k];
          part_sum += elem * elem;
        }

        part_sum = sqrt(part_sum);

        XX[i+j*n] = part_sum;
        row_sum[i] += part_sum;
        col_sum[j] += part_sum;
        total_sum += part_sum;
      } else {
        XX[i+j*n] = 0.0;
      }
    }
  }

  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      XX[i+j*n] -= row_sum[i] / n + col_sum[j] / n - total_sum / n / n;
    }
  }

  free(row_sum);
  free(col_sum);
}

// XX is a nxn matrix, YY is a nxn matrix
double Inner_Prod(int n, double *XX, double *YY) {
  double sum = 0.0; 
  int i, j;

  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      // XX[i, j] * YY[i, j]
      sum += XX[i+j*n] * YY[i+j*n];
    }
  }

  return sum / n / n;
}

// XX is a nxn matrix, YY is a nxn matrix
double Inner_Prod_Perm(int n, int *P, double *XX, double *YY) {
  double sum = 0.0; 
  int i, j;

  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      // XX[i, j] * YY[P[i], P[j]]
      sum += XX[i+j*n] * YY[P[i]+P[j]*n];
    }
  }

  return sum / n / n;
}

//
int main(int argc, char** argv)
{
  //
  double t0 = omp_get_wtime();

  int i, j;
  int nobs = 25;
  int ndim = 5;
  int nrep = 24;
  int nperm = 100;
  double alpha = 0.1;

  extern char* optarg;
  const char* optstring = "o:d:r:p:a:";
  int c;
  while ((c = getopt(argc, argv, optstring)) != -1) {
    switch (c) {
    case 'o': nobs = atoi(optarg); break; // number of observations
    case 'd': ndim = atoi(optarg); break; // number of dimensions
    case 'r': nrep = atoi(optarg); break; // number of repetitions
    case 'p': nperm = atoi(optarg); break; // number of permutations
    case 'a': alpha = atof(optarg); break; // significance level
    }
  }
  
  printf("====================\n");
  printf("nobs: %d\n", nobs);
  printf("ndim: %d\n", ndim);
  printf("nrep: %d\n", nrep);
  printf("nperm: %d\n", nperm);
  printf("alpha: %g\n", alpha);

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
  int* index;
  double *x, *y, *xx, *yy; 
  double stat, stat_perm;
  int count, local;
  int global = 0;

  index = (int*) malloc(nobs*sizeof(int));
  x = (double*) malloc(nobs*ndim*sizeof(double));
  y = (double*) malloc(nobs*ndim*sizeof(double));
  xx = (double*) malloc(nobs*nobs*sizeof(double));
  yy = (double*) malloc(nobs*nobs*sizeof(double));
  local = 0;

  for (i = 0; i < nobs; ++i) {
    index[i] = i;
  }

  for (j = 0; j < nrep; ++j) {
    for (i = 0; i < nobs*ndim; ++i) {
      x[i] = RngStream_RandNormal(RngArray[j]);
    }

    for (i = 0; i < nobs*ndim; ++i) {
      y[i] = RngStream_RandNormal(RngArray[j]);
    }

    Double_Center(nobs, ndim, x, xx);
    Double_Center(nobs, ndim, y, yy);

    stat = Inner_Prod(nobs, xx, yy);
    count = 0;

    for (i = 0; i < nperm; ++i) {
      RngStream_RandShuffle(RngArray[j], index, nobs);
      stat_perm = Inner_Prod_Perm(nobs, index, xx, yy);
      if (stat_perm > stat) {
        count += 1;
      }
    }

    if ((double) count / nperm < alpha) {
      local += 1;  
    }  
  }

  global += local;

  free(index);
  free(x);
  free(y);
  free(xx);
  free(yy);

  double t1 = omp_get_wtime();
  
  printf("====================\n");
  printf("size: %d / %d = %g\n", global, nrep, (double) global / nrep);
  printf("time: %g\n", t1-t0);
  printf("====================\n");

  return 0;
}
