// driverLTM_openmp.cpp
// Using RngStreams for Parallel Random Number Generation in C++ and R
// Computational Statistics
// Andrew Karl, Randy Eubank, Jelena Milovanovic, Mark
// Reiser, Dennis Young akarl@asu.edu
#include <iostream>
#include <iomanip>
#include <omp.h>
#include "tnt_array1d.h"
#include "tnt_array1d_utils.h"
#include "tnt_array2d.h"
#include "tnt_array2d_utils.h"
#include "tnt_math_utils.h"
#include "funclass.h"
#include "RngStream.h"

int main(int argc, char** argv){
  //wallclock time
  double wtime = omp_get_wtime();

  //openMP related parameters//
  int nP = atoi(argv[1]);
  omp_set_num_threads(nP);//set number of threads                              
  int myRank;

  //data features//
  //generate all possible response patterns
  TNT::Array2D<int> patn(32, 5, 0);
  funClass::genPats(0, 31, 0, patn);
  //pointer for intercepts
  double* pAlpha = new double[5];
  //intercept values
  pAlpha[0] = -1.868265; pAlpha[1] = -0.7910062; pAlpha[2] = -1.460980;
  pAlpha[3] = -0.521506; pAlpha[4] = -1.992975;
  //slope
  double beta = 1.011268;

  //simulation parameters//
  int nRep = 1e5;
  int* nSim = new int[nP];
  //divide up the work to be done for each thread
  int rem = nRep%nP;
  for(int i = 0; i < nP; i++){
    nSim[i] = nRep/nP;
    if(rem != 0) {
      nSim[i]++;
      rem--;
    }
  }

  //initialize RngStream array
  unsigned long seed[6] = {1806547166, 3311292359,  
			  643431772, 1162448557, 
			  3335719306, 4161054083};
  RngStream::SetPackageSeed(seed);
  RngStream RngArray[nP];

  //pointer to Array1D objects that will hold output from each process
  TNT::Array1D<double>* tempVec = new TNT::Array1D<double>[nP];

  //pointer for random numbers and funClass object
  double* pZ;
  funClass* pf;
  //initialize array for cell probabilities
  TNT::Array1D<double> prob(patn.dim1(), 0.);

 //approximate the cell probabilities  
#pragma omp parallel private(myRank, pZ, pf)
  {
    myRank = omp_get_thread_num();
    pZ = new double[nSim[myRank]];
    pf = new funClass(pAlpha, beta, patn.dim2());
    for(int i = 0;i < nSim[myRank]; i++)
      pZ[i] =  2.*pf->qnorm(RngArray[myRank].RandU01());
    tempVec[myRank] 
      = pf->probFunc(TNT::Array1D<double>(nSim[myRank], 
					  pZ), patn); 
    delete[] pZ;
    delete pf;
  }

  for(int i = 0; i < nP; i++){
    prob += tempVec[i];
  }
  std::cout << "Integration results:"<<std::endl; 
  std::cout << std::fixed;
  for(int i = 0; i < patn.dim1(); i++){
    prob[i] /= (double)nP;
    std::cout<<std::setprecision(9)<<prob[i]<<std::endl;
  }

  std::cout <<"time with " << nP 
    << " processes was " << omp_get_wtime() - wtime 
    << std::endl;
  delete[] pAlpha;
  delete[] nSim;  
  delete[] tempVec;
  return 0;
}

