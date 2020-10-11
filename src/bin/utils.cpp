#include "../include/utils.h"
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <vector>

using namespace std;


bool canPair(int n){
    if(n == 9 || n == 13 || n == 17 || n == 19 || n == 21 || n == 23) return true;
    return false;
}

/** This functions takes a string and two chars. If char1 appears, it is replaced with char2**/
char* replaceChar(char *stri, char ch1, char ch2) {
    string seq(stri);
    // cout << seq.length() << endl;
    for (int i = 0; i < seq.length(); ++i) {
        if (seq[i] == ch1)
            seq[i] = ch2;
    }
    // char output[seq.length()+1];
    strcpy(stri, seq.c_str());
    return stri;
}