#ifndef ALLFOLD
#define ALLFOLD
#include <vector>
#include <string>

std::string returnUngapped(std::string input_sequence, std::string consensus_structure);

bool canMatch(char x, char y);

double MI(std::vector<std::string> seqs, int i, int j);

std::string get_consensus(std::vector<std::string> vect, bool mis);

bool call_simfold2 (char *programPath, char *input_sequence, char *output_structure, double *output_energy);

bool call_simfold3 (char *programPath, char *input_sequence, char *output_structure, double *output_energy, double *scores, int n);

void rnaalifold(std::string sequences, std::string outputFile, std::string &structure);

std::vector <std::string> returnSeq(std::string address);

void runPerl(std::string family);

void tCoffee(std::string family, std::string dir);

std::string iterativeFold(std::string seq, std::string str);

void convertPS(std::string family, std::string name);


#endif