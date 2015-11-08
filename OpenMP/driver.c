#include "RngStream.h"
#include <omp.h>
#include <stdio.h>

int main(){
 
  int np = omp_get_num_procs();
  omp_set_num_threads(np);//set number of threads                              

  unsigned long seed[6] = {1806547166, 3311292359,  
			  643431772, 1162448557, 
			  3335719306, 4161054083};
  RngStream_SetPackageSeed(seed);
  RngStream RngArray[np];//array of RngStream objects

  int myRank;
  #pragma omp parallel private(myRank)
  {
    myRank = omp_get_thread_num();
    #pragma omp critical
    {
      printf ( "RandU01 (g1) = %16.12f\n", RngStream_RandU01(RngArray[myRank]));
    }
  }

  return 0;
}