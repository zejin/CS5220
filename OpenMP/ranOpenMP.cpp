//ranOpenMP.cpp
//Using RngStreams for Parallel Random Number Generation in C++ and R
//Computational Statistics
//Andrew Karl, Randy Eubank, Jelena Milovanovic, Mark Reiser, Dennis Young
//akarl@asu.edu
#include <omp.h>
#include "RngStream.h"
#include <iostream>

int main(){
 
  int nP = omp_get_num_procs();
  omp_set_num_threads(nP);//set number of threads                              

  unsigned long seed[6] ={1806547166, 3311292359,  
			  643431772, 1162448557, 
			  3335719306, 4161054083};
  RngStream::SetPackageSeed (seed);
  RngStream RngArray[nP];//array of RngStream objects

  int myRank;
#pragma omp parallel private(myRank)
  {
    myRank = omp_get_thread_num();
#pragma omp critical
    {
      std::cout << RngArray[myRank].RandU01() << " "; 
    }
  }
  std::cout << std::endl;
  return 0;
}
