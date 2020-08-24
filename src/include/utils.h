#ifndef UTILS
#define UTILS
#include <vector>
#include <string>

int* toBinary(int n, int ns);

int hamming(int i, int j, int x, int y);

void ct2dotSin(std::string &family, std::vector <std::string> &list);

void ct2dot2(std::string &family, std::vector <std::string> &lists);

char* replaceChar(char *str, char ch1, char ch2);

bool getFileContent(std::string fileName, std::vector <std::string> &vecOfStrs);

bool fexists(const char *filename);

bool canPair(int n);




#endif