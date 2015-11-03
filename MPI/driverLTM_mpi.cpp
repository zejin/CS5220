// driverLTM_mpi.cpp
// Using RngStreams for Parallel Random Number Generation in C++ and R
// Computational Statistics
// Andrew Karl, Randy Eubank, Jelena Milovanovic, Mark
// Reiser, Dennis Young akarl@asu.edu
#include <iostream>
#include <mpi.h>
#include "tnt_array1d.h"
#include "tnt_array1d_utils.h"
#include "tnt_array2d.h"
#include "tnt_array2d_utils.h"
#include "tnt_math_utils.h"
#include "funclass.h"
#include "RngStream.h"

int main(){

  MPI::Init();
  //wallclock time
  double wtime = MPI::Wtime();

  int myRank = MPI::COMM_WORLD.Get_rank();
  int nP = MPI::COMM_WORLD.Get_size();

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
  //divide up the work to be done for each process
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

  //array for random numbers
  TNT::Array1D<double> Z(nSim[myRank], 0.);
  //funClass object
  funClass f(pAlpha, beta, patn.dim2());

  //array that will hold approximate probabilities
  TNT::Array1D<double> prob(patn.dim1(), 0.);

 //approximate the cell probabilities on each process  
  for(int i = 0;i < nSim[myRank]; i++)
    Z[i] =  2.*f.qnorm(RngArray[myRank].RandU01());
  prob = f.probFunc(Z, patn); 

  //now gather up the results on the master process
  double* temp;
  temp = new double[nP*patn.dim1()];

  double* pData_ = new double[patn.dim1()];
  for(int i = 0; i < patn.dim1(); i++)
    pData_[i] = prob[i];

  MPI::COMM_WORLD.Gather(pData_, patn.dim1(), MPI::DOUBLE, 
			 temp, patn.dim1(), MPI::DOUBLE, 0);
  if(myRank == 0){  
    std::cout << "Integration results:"<<std::endl; 
    for(int i = 0; i < prob.dim(); i++){
      for(int j = 1; j < nP; j++)
         prob[i] += temp[j*prob.dim() + i];
      prob[i] /= (double)nP;
      std::cout<<prob[i]<<std::endl;
    }
    std::cout << "time with " << nP << " processes was " << 
      MPI::Wtime() - wtime << std::endl;
  }

  delete[] pAlpha;
  delete[] nSim;
  delete[] temp;
  delete[] pData_; 
  MPI::Finalize();
  return 0;
}

