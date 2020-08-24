#ifndef ALLFOLD
#define ALLFOLD
#include <vector>
#include <string>

std::string returnUngapped(std::string input_sequence, std::string consensus_structure);

bool canMatch(char x, char y);

bool call_simfold2 (char *programPath, char *input_sequence, char *output_structure, double *output_energy);

bool call_simfold3 (char *programPath, char *input_sequence, char *output_structure, double *output_energy, double *scores, int n);

std::string iterativeFold(std::string seq, std::string str);


#endif