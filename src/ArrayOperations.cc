#include "ArrayOperations.hh"
#include <errno.h>
#include <cmath>
#include <stdlib.h>     /* calloc, exit, free */


int maxState(uint *vector, int vectorLength) {
    uint max = 0;
    for (int i = 0; i < vectorLength; i++) {
        if (vector[i] > max) {
            max = vector[i];
        }
    }
    return max + 1;
}
/******************************************************************************* 
 ** normaliseArray takes an input vector and writes an output vector
 ** which is a normalised version of the input, and returns the number of states
 ** A normalised array has min value = 0, max value = number of states
 ** and all values are integers
 **
 ** length(inputVector) == length(outputVector) == vectorLength otherwise there
 ** is a memory leak
 *******************************************************************************/
int normaliseArray(double *inputVector, uint *outputVector, int vectorLength) {
    int maxVal = 0;

    if (vectorLength > 0) {
        int* tempVector = (int*) calloc(vectorLength,sizeof(int));
        int minVal = (int) floor(inputVector[0]);
        maxVal = (int) floor(inputVector[0]);

        for (int i = 0; i < vectorLength; i++) {
            int currentValue = (int) floor(inputVector[i]);
            tempVector[i] = currentValue;

            if (currentValue < minVal) {
                minVal = currentValue;
            } else if (currentValue > maxVal) {
                maxVal = currentValue;
            }
        }/*for loop over vector*/

        for (int i = 0; i < vectorLength; i++) {
            outputVector[i] = tempVector[i] - minVal;
        }

        maxVal = (maxVal - minVal) + 1;

        free(tempVector);
        
    }

    return maxVal;
}