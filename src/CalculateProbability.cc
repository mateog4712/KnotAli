#include "CalculateProbability.hh"
#include "ArrayOperations.hh"
#include <iostream>
#include <stdlib.h>     /* calloc, exit, free */
#include <stdio.h>      /* printf, scanf, NULL */


JointProbabilityState calculateJointProbability(uint *firstVector, uint *secondVector, int vectorLength) {

  int *firstStateCounts;
  int *secondStateCounts;
  int *jointStateCounts;
  double *firstStateProbs;
  double *secondStateProbs;
  double *jointStateProbs;
  int firstNumStates;
  int secondNumStates;
  int jointNumStates;
  int i;
  double length = vectorLength;
  JointProbabilityState state;

  firstNumStates = maxState(firstVector,vectorLength);
  secondNumStates = maxState(secondVector,vectorLength);
  jointNumStates = firstNumStates * secondNumStates;
  
  firstStateCounts = (int *) calloc(firstNumStates,sizeof(int));
  secondStateCounts = (int *) calloc(secondNumStates,sizeof(int));
  jointStateCounts = (int *) calloc(jointNumStates,sizeof(int));
  
  firstStateProbs = (double *) calloc(firstNumStates,sizeof(double));
  secondStateProbs = (double *) calloc(secondNumStates,sizeof(double));
  jointStateProbs = (double *) calloc(jointNumStates,sizeof(double));
    
  /* Optimised for number of FP operations now O(states) instead of O(vectorLength) */
  for (i = 0; i < vectorLength; i++) {
    firstStateCounts[firstVector[i]] += 1;
    secondStateCounts[secondVector[i]] += 1;
    jointStateCounts[secondVector[i] * firstNumStates + firstVector[i]] += 1;
  }
  
  for (i = 0; i < firstNumStates; i++) {
    firstStateProbs[i] = firstStateCounts[i] / length;
  }
  
  for (i = 0; i < secondNumStates; i++) {
    secondStateProbs[i] = secondStateCounts[i] / length;
  }
  
  for (i = 0; i < jointNumStates; i++) {
    jointStateProbs[i] = jointStateCounts[i] / length;
  }
  
  state.jointProbabilityVector = jointStateProbs;
  state.numJointStates = jointNumStates;
  state.firstProbabilityVector = firstStateProbs;
  state.numFirstStates = firstNumStates;
  state.secondProbabilityVector = secondStateProbs;
  state.numSecondStates = secondNumStates;
  
  free(firstStateCounts);
  free(secondStateCounts);
  free(jointStateCounts);

  return state;
}
JointProbabilityState discAndCalcJointProbability(double *firstVector, double *secondVector, int vectorLength) {
  uint *firstNormalisedVector = (uint *) calloc(vectorLength,sizeof(uint));
  uint *secondNormalisedVector = (uint *) calloc(vectorLength,sizeof(uint));
  JointProbabilityState state;
  
  normaliseArray(firstVector,firstNormalisedVector,vectorLength);
  normaliseArray(secondVector,secondNormalisedVector,vectorLength);

  state = calculateJointProbability(firstNormalisedVector,secondNormalisedVector,vectorLength);

  free(firstNormalisedVector);
  free(secondNormalisedVector);

  return state;
}






ProbabilityState calculateProbability(uint *dataVector, int vectorLength) {
  int numStates;
  int *stateCounts;
  double *stateProbs;
  ProbabilityState state;
  int i;
  double length = vectorLength;

  numStates = maxState(dataVector,vectorLength);
  
  stateCounts = (int *) calloc(numStates,sizeof(int));
  stateProbs = (double *) calloc(numStates,sizeof(double));
  
  /* Optimised for number of FP operations now O(states) instead of O(vectorLength) */
  for (i = 0; i < vectorLength; i++) {
    stateCounts[dataVector[i]] += 1;
  }
  
  for (i = 0; i < numStates; i++) {
    stateProbs[i] = stateCounts[i] / length;
  }
  
  free(stateCounts);
  stateCounts = NULL;
  
  state.probabilityVector = stateProbs;
  state.numStates = numStates;

  return state;
}









/////////////////////////// FREE FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

void freeProbabilityState(ProbabilityState state) {
    free(state.probabilityVector);
    state.probabilityVector = NULL;
}

void freeJointProbabilityState(JointProbabilityState state) {
    free(state.firstProbabilityVector);
    state.firstProbabilityVector = NULL;
    free(state.secondProbabilityVector);
    state.secondProbabilityVector = NULL;
    free(state.jointProbabilityVector);
    state.jointProbabilityVector = NULL;
}