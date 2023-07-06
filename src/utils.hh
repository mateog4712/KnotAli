#ifndef UTILS
#define UTILS
#include <vector>
#include <string>

char* replaceChar(char *str, char ch1, char ch2);

bool exists (const std::string& name);

std::string returnUngapped(std::string input_sequence, std::string consensus_structure);

bool canMatch(char x, char y);

bool call_simfold2 (char *programPath, char *input_sequence, char *output_structure, double *output_energy);

bool call_simfold3 (char *programPath, char *input_sequence, char *output_structure, double *output_energy, double *scores, int n);

std::string iterativeFold(std::string seq, std::string str, double &en);

#endif