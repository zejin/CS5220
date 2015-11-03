// RngStreamSupp.h
// Using RngStreams for Parallel Random Number Generation in C++ and R
// Computational Statistics
// Andrew Karl, Randy Eubank, Jelena Milovanovic, Mark
// Reiser, Dennis Young akarl@asu.edu
#ifndef RNGSTREAMSUPP_H
#define RNGSTREAMSUPP_H

class RngStreamSupp{

 public:

  static void AdvanceSeed (unsigned long seedIn[6], unsigned long seedOut[6]);

};

#endif
