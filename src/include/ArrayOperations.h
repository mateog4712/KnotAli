#ifndef ARRAY_OPERATIONS
#define ARRAY_OPERATIONS
#include <sys/types.h>

// /*******************************************************************************
// ** Finds the maximum state of an int array.
// *******************************************************************************/
int maxState(uint *vector, int length);

/******************************************************************************* 
** normaliseArray takes an input vector and writes an output vector
** which is a normalised version of the input, and returns the number of states
** A normalised array has min value = 0, max value = number of states
** and all values are integers
*******************************************************************************/
int normaliseArray(double *inputVector, uint *outputVector, int vectorLength);

#endif