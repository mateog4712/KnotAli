#ifndef MUTUAL_INFORMATION
#define MUTUAL_INFORMATION
#include <vector>
#include <string>
#include "CalculateProbability.h"
#include <sys/types.h>

std::string MIVector(std::vector<std::string> seqs);
std::string MIVectorStack(std::vector<std::string> seqs);

/*******************************************************************************
** calculateMutualInformation returns the log base LOG_BASE mutual information between
** dataVector and targetVector, I(X;Y)

*******************************************************************************/
double calcMutualInformation(uint *dataVector, uint *targetVector, int vectorLength);
double discAndCalcMutualInformation(double *dataVector, double *targetVector, int vectorLength);

/*******************************************************************************
** Inner functions which operate on state structs.
*******************************************************************************/
double mi(JointProbabilityState state);




#endif