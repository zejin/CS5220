// funclass.cpp
// Using RngStreams for Parallel Random Number Generation in C++ and R
// Computational Statistics
// Andrew Karl, Randy Eubank, Jelena Milovanovic, Mark
// Reiser, Dennis Young akarl@asu.edu
#include <cmath>
#include "tnt_array1d.h"
#include "tnt_array1d_utils.h"
#include "tnt_array2d.h"
#include "tnt_array2d_utils.h"
#include "tnt_math_utils.h"
#include "funclass.h"

//anonymous namespace for constants
namespace{
  // Coefficients in rational approximations
  static const double a[6] = {-3.969683028665376e+01,  2.209460984245205e+02,
		       -2.759285104469687e+02,  1.383577518672690e+02,
		       -3.066479806614716e+01,  2.506628277459239e+00};

  static const double b[5] = {-5.447609879822406e+01,  1.615858368580409e+02,
		       -1.556989798598866e+02,  6.680131188771972e+01,
		       -1.328068155288572e+01};
  
  static const double c[6] = {-7.784894002430293e-03, -3.223964580411365e-01,
		       -2.400758277161838e+00, -2.549732539343734e+00,
		       4.374664141464968e+00,  2.938163982698783e+00};

  static const double d[5] = {7.784695709041462e-03, 3.224671290700398e-01,
		       2.445134137142996e+00,  3.754408661907416e+00};

  // Define break-points.
  static const double ulow  = 0.02425;
  static const double uhigh = 1 - ulow;

  static const double pi = 2.*acos(0.);
}

funClass::funClass(const double* pa, double b, int dim){
  beta = b;
  pAlpha = new double[dim];
  for(int i = 0; i < dim; i++){
    pAlpha[i] = pa[i];
  }
}

TNT::Array2D<double> 
funClass::probMatFunc(double z, const TNT::Array2D<int>& patn){
  //temporary storage array
  TNT::Array1D<double> temp(patn.dim2(), 0.);
  //calculate the chance of a success for each variable at this deviate value
  for(int j = 0; j < patn.dim2(); j++){
    temp[j] = 1./(1. + exp(-(-pAlpha[j] + beta*z)));
  }
  
  //evaluate the probability for each variable on every pattern at this deviate value
  TNT::Array2D<double> probMat(patn.dim1(), patn.dim2(), 0.);
  for(int j = 0; j < patn.dim2(); j++){
    for(int i = 0; i < patn.dim1(); i++){
      if(patn[i][j] == 0) probMat[i][j] = 1. - temp[j];
      if(patn[i][j] == 1) probMat[i][j] = temp[j];
    }
  }
  return probMat;  
}

TNT::Array1D<double> 
funClass::probFunc(const TNT::Array1D<double>& Z, 
		   const TNT::Array2D<int>& patn){
  //pattern probability array used in intermediate calculations
  TNT::Array2D<double>  probMat(patn.dim1(), patn.dim2(), 0.);
  //work array
  TNT::Array1D<double> temp(patn.dim1(), 0.);
  //array that will hold the cell probabilities
  TNT::Array1D<double> prob(patn.dim1(), 0.);

  //sum across normal deviates
  for(int k = 0; k < Z.dim(); k++){
    probMat = probMatFunc(Z[k], patn); 
      
    //multiply to get the probability for each pattern
    for(int i = 0; i < patn.dim1(); i++){
      temp[i] = 1.;
      for(int j = 0; j < patn.dim2(); j++)
	temp[i] *= probMat[i][j];
      temp[i] *= 2.*pnorm(Z[k])
	/(pnorm(Z[k]/2.)*(double) Z.dim());
    }
    prob += temp;
  }
  return prob;
}

//utility functions
//C++ translation of Java code by Peter John Acklam 
//(http://home.online.no/~pjacklam) for inverse normal cdf
double funClass::qnorm(double u){
  // Rational approximation for lower region:
  if (u < ulow) {
    double q = sqrt(-2*log(u));
    return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
      ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1);
  }
  
  // Rational approximation for upper region:
  if (uhigh < u) {
    double q = sqrt(-2*log(1 - u));
    return -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
      ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1);
  }
  
  // Rational approximation for central region:
  double q = u - 0.5;
  double r = q*q;
  return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5])*q /
    (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1);
}

//standard normal density
double funClass::pnorm(double x){
  return exp(-x*x/2.)/sqrt(2*pi);
}

void funClass::genPats(int nrowLow, int nrowUp, int ncol, TNT::Array2D<int>& patn){
  int nrows = nrowUp - nrowLow + 1;
  int nrow2 = nrows/2;
  if(ncol < patn.dim2()){
    for(int i = nrowLow; i < nrowLow + nrow2; i++)
      patn[i][ncol] = 0;
    for(int i = nrowLow + nrow2; i <= nrowUp; i++)
      patn[i][ncol] = 1;
    genPats(nrowLow, nrowLow + nrow2 - 1, ncol + 1, patn);
    genPats(nrowLow + nrow2, nrowUp, ncol + 1, patn);
  }
}
