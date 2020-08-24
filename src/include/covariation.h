#ifndef COVARIATION
#define COVARIATION
#include <sys/types.h>
#include <string>
#include <vector>


double covariation(uint *dataVector, uint *targetVector, int vectorLength);

std::string covVectorStack(std::vector<std::string> seqs);

std::string covVector(std::vector<std::string> seqs);

unsigned nChoosek( unsigned n, unsigned k );




#endif