#include <string>
#include <vector>
#include <math.h>       
#include <iostream>
#include <numeric>
#include "/home/mgray7/KnotAli/src/include/RenyiMutualInformation.h"
#include "../include/CalculateProbability.h"

#include <map>
#include <tuple>
#include <algorithm>



double renyiMI(JointProbabilityState state, double alpha) {
  int firstIndex,secondIndex;
  int i;
  double jointTemp, marginalTemp;
  double invAlpha = 1.0 - alpha;
  double mutualInformation = 0.0;

  
  for (i = 0; i < state.numJointStates; i++) {
    firstIndex = i % state.numFirstStates;
    secondIndex = i / state.numFirstStates;
    
    if ((state.jointProbabilityVector[i] > 0) && (state.firstProbabilityVector[firstIndex] > 0) && (state.secondProbabilityVector[secondIndex] > 0)) {
      jointTemp = pow(state.jointProbabilityVector[i],alpha);
      marginalTemp = state.firstProbabilityVector[firstIndex] * state.secondProbabilityVector[secondIndex];
      marginalTemp = pow(marginalTemp,invAlpha);
      mutualInformation += (jointTemp * marginalTemp);
    }
  }

  mutualInformation = log(mutualInformation);
  mutualInformation /= log(2);  
  mutualInformation /= (alpha-1.0); 

  return mutualInformation;
}


double calcRenyiMIDivergence(double alpha, uint *dataVector, uint *targetVector, int vectorLength) {
  JointProbabilityState state = calculateJointProbability(dataVector,targetVector,vectorLength);
  double mutualInformation = renyiMI(state,alpha);
  
  freeJointProbabilityState(state);
  
  return mutualInformation;
}

double discAndCalcRenyiMIDivergence(double alpha, double *dataVector, double *targetVector, int vectorLength) {
  JointProbabilityState state = discAndCalcJointProbability(dataVector,targetVector,vectorLength);
  double mutualInformation = renyiMI(state,alpha);
  
  freeJointProbabilityState(state);
  
  return mutualInformation;
}