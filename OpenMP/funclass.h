// funclass.h
// Using RngStreams for Parallel Random Number Generation in C++ and R
// Computational Statistics
// Andrew Karl, Randy Eubank, Jelena Milovanovic, Mark
// Reiser, Dennis Young akarl@asu.edu
#ifndef FUNCLASS_H
#define FUNCLASS_H

namespace TNT{
  template <class T> class Array1D;
  template <class T> class Array2D;
}

#include <iostream>

class funClass{
  double* pAlpha; 
  double beta;

  TNT::Array2D<double> probMatFunc(double z, const TNT::Array2D<int>& patn);
 
 public:
  
  funClass(const double* pa, double b, int dim);
  ~funClass(){
    delete[] pAlpha;
  }

  static double qnorm(double u);
  static double pnorm(double x);
  static void genPats(int nrowLow, int nrowUp, int ncol, TNT::Array2D<int>& patn);

  TNT::Array1D<double> probFunc(const TNT::Array1D<double>& Z, const TNT::Array2D<int>& patn);
};
#endif
