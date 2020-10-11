#ifndef MUTUAL_INFORMATION
#define MUTUAL_INFORMATION
#include "CalculateProbability.h"
#include <vector>
#include <string>
#include <tuple>

// structure holding the score and location of pair
struct Hotspot { 
    std::tuple<int,int> pair;
    double score;
};

std::string MIVector(std::vector<std::string> seqs, bool stack = false);
bool check_Pseudoknot(std::vector<std::tuple<int,int> > used, Hotspot hotspot);

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