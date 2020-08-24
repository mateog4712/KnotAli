#ifndef __Renyi_MutualInformation
#define __Renyi_MutualInformation
#include "CalculateProbability.h"
/*******************************************************************************
** calculateRenyiMIDivergence returns the log base LOG_BASE Renyi mutual information
** between dataVector and targetVector, I_{\alpha}(X;Y), for \alpha != 1
** This uses Renyi's generalised alpha divergence as the difference measure
** instead of the KL-divergence as in Shannon's Mutual Information
*******************************************************************************/
double calcRenyiMIDivergence(double alpha, uint *dataVector, uint *targetVector, int vectorLength);
double discAndCalcRenyiMIDivergence(double alpha, double *dataVector, double *targetVector, int vectorLength);






/*******************************************************************************
** Inner functions which operate on state structs.
*******************************************************************************/
double renyiMI(JointProbabilityState state, double alpha);

#endif