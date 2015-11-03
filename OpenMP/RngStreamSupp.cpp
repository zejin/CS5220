// RngStreamSupp.cpp
// supplement to RngStream that allows for manual advancement of 
// the generator's seeds
// Using RngStreams for Parallel Random Number Generation in C++ and R
// Computational Statistics
// Andrew Karl, Randy Eubank, Jelena Milovanovic, Mark
// Reiser, Dennis Young akarl@asu.edu

#include "RngStreamSupp.h"
using namespace std;

//anonymous namespace that contains the necessary code from RngStream.cpp
namespace
{
const double m1   =       4294967087.0;
const double m2   =       4294944443.0;
const double two17 =      131072.0;
const double two53 =      9007199254740992.0;

const double A1p127[3][3] = {
       {    2427906178.0, 3580155704.0,  949770784.0 },
       {     226153695.0, 1230515664.0, 3580155704.0 },
       {    1988835001.0,  986791581.0, 1230515664.0 }
       };

const double A2p127[3][3] = {
       {    1464411153.0,  277697599.0, 1610723613.0 },
       {      32183930.0, 1464411153.0, 1022607788.0 },
       {    2824425944.0,   32183930.0, 2093834863.0 }
       };



//-------------------------------------------------------------------------
double MultModM (double a, double s, double c, double m)
{
    double v;
    long a1;

    v = a * s + c;

    if (v >= two53 || v <= -two53) {
        a1 = static_cast<long> (a / two17);    a -= a1 * two17;
        v  = a1 * s;
        a1 = static_cast<long> (v / m);     v -= a1 * m;
        v = v * two17 + a * s + c;
    }

    a1 = static_cast<long> (v / m);
    /* in case v < 0)*/
    if ((v -= a1 * m) < 0.0) return v += m;   else return v;
}


//-------------------------------------------------------------------------
void MatVecModM (const double A[3][3], const double s[3], double v[3],
                 double m)
{
    int i;
    double x[3];               // Necessary if v = s

    for (i = 0; i < 3; ++i) {
        x[i] = MultModM (A[i][0], s[0], 0.0, m);
        x[i] = MultModM (A[i][1], s[1], x[i], m);
        x[i] = MultModM (A[i][2], s[2], x[i], m);
    }
    for (i = 0; i < 3; ++i)
        v[i] = x[i];
}

}

//-------------------------------------------------------------------------
// Generate the seed for the next stream
//
void RngStreamSupp::AdvanceSeed (unsigned long seedIn[6], unsigned long seedOut[6])
{
  double tempIn[6]; double tempOut[6];
  for(int i = 0; i < 6; i++)
    tempIn[i] = seedIn[i];
  MatVecModM (A1p127, tempIn, tempOut, m1);
  MatVecModM (A2p127, &tempIn[3], &tempOut[3], m2);
  for(int i = 0; i < 6; i++)
    seedOut[i] = tempOut[i];

};
