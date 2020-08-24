#ifndef STRUCTS
#define STRUCTS
#include <sys/types.h>

typedef struct jpstate{

  double *jointProbabilityVector;
  int numJointStates;
  double *firstProbabilityVector;
  int numFirstStates;
  double *secondProbabilityVector;
  int numSecondStates;
} JointProbabilityState;

typedef struct pState
{
  double *probabilityVector;
  int numStates;
} ProbabilityState;

/*******************************************************************************
** calculateJointProbability returns the joint probability vector of two vectors
** and the marginal probability vectors in a struct.
*******************************************************************************/
JointProbabilityState calculateJointProbability(uint *firstVector, uint *secondVector, int vectorLength);
/*******************************************************************************
** discAndCalcJointProbability discretises the double vectors into int vectors,
** then generates a JointProbabilityState.
*******************************************************************************/
JointProbabilityState discAndCalcJointProbability(double *firstVector, double *secondVector, int vectorLength);










/*******************************************************************************
** calculateProbability returns the probability vector from one vector.
** It is the base operation for all information theory calculations involving 
** one variable
*******************************************************************************/
ProbabilityState calculateProbability(uint *dataVector, int vectorLength);





/*******************************************************************************
** Frees the struct members and sets all pointers to NULL.
*******************************************************************************/
void freeProbabilityState(ProbabilityState state);
void freeJointProbabilityState(JointProbabilityState state);


#endif