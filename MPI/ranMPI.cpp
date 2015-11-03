//ranMPI.cpp
//Using RngStreams for Parallel Random Number Generation in C++ and R
//Computational Statistics
//Andrew Karl, Randy Eubank, Jelena Milovanovic, Mark Reiser, Dennis Young
//akarl@asu.edu
#include <iostream>
#include "mpi.h"
#include "RngStream.h"
#include "RngStreamSupp.h"

int main(){
  unsigned long seed[6] ={1806547166, 3311292359,  
			  643431772, 1162448557, 
			  3335719306, 4161054083};
  MPI::Init();//start MPI
  int myRank = MPI::COMM_WORLD.Get_rank();//get process rank
  int nP = MPI::COMM_WORLD.Get_size();//get number of processes
  if(myRank == 0){
    //generate deviate for master process
    RngStream::SetPackageSeed(seed);
    RngStream Rng;
    double U = Rng.RandU01();
    //now send seeds to the other processes
    unsigned long tempSeed[6];
    for(int i = 1; i < nP; i++){
      RngStreamSupp::AdvanceSeed(seed, tempSeed);
      MPI::COMM_WORLD.Send(tempSeed, 6, MPI::UNSIGNED_LONG, 
			   i, 0);
      for(int j = 0; j < 6; j++)
        seed[j] = tempSeed[j];
    }
    std::cout << U << " ";
    //collect the other random deviates
    for(int i = 1; i < nP; i++){
      MPI::COMM_WORLD.Recv(&U, 1, MPI::DOUBLE, i, 0);
      std::cout << U << " ";
    }      
    std::cout << std::endl;
  }     
  else{
    unsigned long mySeed[6];
    //receive your seed if you are not the master process
    MPI::COMM_WORLD.Recv(mySeed, 6, MPI::UNSIGNED_LONG, 0, 0);
    RngStream::SetPackageSeed(mySeed);
    RngStream Rng;
    double U = Rng.RandU01();
    //send deviate to master process
    MPI::COMM_WORLD.Send(&U, 1, MPI::DOUBLE, 0, 0);
  }
  MPI::Finalize();//end MPI
  return 0;
}
