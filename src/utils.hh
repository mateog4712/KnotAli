#ifndef UTILS
#define UTILS
#include <vector>
#include <string>

std::string removeIUPAC(std::string sequence);

char* replaceChar(char *str, char ch1, char ch2);

bool exists (const std::string& name);

std::string returnUngapped(std::string input_sequence, std::string consensus_structure);

bool canMatch(char x, char y);

#endif