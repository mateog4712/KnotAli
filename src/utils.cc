#include "utils.hh"
#include <sys/stat.h>
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <vector>
#include <tuple>

std::string removeIUPAC(std::string sequence){
    int n = sequence.length();
    for(int i = 0; i<n; ++i){
        if(sequence[i] == 'R') sequence[i] = 'G';
        else if(sequence[i] == 'Y') sequence[i] = 'C';
        else if (sequence[i] == 'S') sequence[i] = 'G';
        else if (sequence[i] == 'W') sequence[i] = 'A';
        else if (sequence[i] == 'K') sequence[i] = 'T';
        else if (sequence[i] == 'M') sequence[i] = 'A';
        else if (sequence[i] == 'B') sequence[i] = 'C';
        else if (sequence[i] == 'D') sequence[i] = 'T';
        else if (sequence[i] == 'H') sequence[i] = 'A';
        else if (sequence[i] == 'V') sequence[i] = 'G';
        else if (sequence[i] == 'N') sequence[i] = 'A';

    }
    return sequence;
}


bool exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

/** This functions takes a string and two chars. If char1 appears, it is replaced with char2**/
char* replaceChar(char *stri, char ch1, char ch2) {
    std::string seq(stri);
    // cout << seq.length() << endl;
    for (int i = 0; i < seq.length(); ++i) {
        if (seq[i] == ch1)
            seq[i] = ch2;
    }
    // char output[seq.length()+1];
    strcpy(stri, seq.c_str());
    return stri;
}

/**
 *  Takes in the consensus sequence and structure and the original sequence and returns the ungapped version 
 * of the structure to be used in simfold
 *  Works with pseudoknotted structures
 **/
std::string returnUngapped(std::string input_sequence, std::string consensus_structure){
    int length = consensus_structure.length();

    std::vector<std::tuple<char,int> > paren;
    std::vector<std::tuple<char,int> > sb;
    std::vector<std::tuple<char,int> > cb;
    std::vector<std::tuple<char,int> > lts;
    for(int i=0;i<length;i++){
        if(consensus_structure[i] == '(') {
            paren.push_back(std::make_tuple(consensus_structure[i],i));
            continue;
        }
        else if(consensus_structure[i] == '<') {
            lts.push_back(std::make_tuple(consensus_structure[i],i));
            continue;
        }
        else if(consensus_structure[i] == '[') {
            sb.push_back(std::make_tuple(consensus_structure[i],i));
            continue;
        }
        else if(consensus_structure[i] == '{') {
            cb.push_back(std::make_tuple(consensus_structure[i],i));
            continue;
        }
        std::tuple<char,int> x;
        bool close = false;
        if (consensus_structure[i] == ')' && !paren.empty()){
            x = paren[paren.size()-1];
            paren.erase(paren.begin()+(paren.size()-1));
            close = true;

        }
        else if (consensus_structure[i] == '>' && !lts.empty()){
            x = lts[lts.size()-1];
            lts.erase(lts.begin()+(lts.size()-1));
            close = true;
        }
        else if (consensus_structure[i] == ']' && !sb.empty()){
            x = sb[sb.size()-1];
            sb.erase(sb.begin()+(sb.size()-1));
            close = true;
        }
        else if (consensus_structure[i] == '}' && !cb.empty()){
            x = cb[cb.size()-1];
            cb.erase(cb.begin()+(cb.size()-1));
            close = true;
        }
            
        // Checks to see if the pair can match
        if(close){
            if(canMatch(input_sequence[i],input_sequence[std::get<1>(x)])){
                
                continue;
            }
            else{
                if(input_sequence[i] == '-'){
                    consensus_structure[std::get<1>(x)] = '.';
                    continue;
                }else if(input_sequence[std::get<1>(x)] == '-'){
                    consensus_structure[i] = '.';
                    continue;   
                }
                else{
                    consensus_structure[i] = '.';
                    consensus_structure[std::get<1>(x)] = '.';
                    continue;
                }
            }
        }    
    }
    // Reports if not all pairs were closed
    if(!paren.empty() || !lts.empty() || !sb.empty() || !cb.empty()){
        std::cout << "Error in reported structure. More left pairings than right pairings" << std::endl;
        exit(0); 
    }
    // Erase Gaps
    for(int i= input_sequence.length()-1; i>=0;--i){
        if(input_sequence[i] == '-') consensus_structure.erase(i,1);
    }

    // Checks for pairs less than 3 distance away
    for(int i=0;i<length;i++){
        if(consensus_structure[i] == '(') {
            paren.push_back(std::make_tuple(consensus_structure[i],i));
            continue;
        }
        else if(consensus_structure[i] == '<') {
            lts.push_back(std::make_tuple(consensus_structure[i],i));
            continue;
        }
        else if(consensus_structure[i] == '[') {
            sb.push_back(std::make_tuple(consensus_structure[i],i));
            continue;
        }
        else if(consensus_structure[i] == '{') {
            cb.push_back(std::make_tuple(consensus_structure[i],i));
            continue;
        }
        std::tuple<char,int> x;
        bool close = false;
        if (consensus_structure[i] == ')' && !paren.empty()){
            x = paren[paren.size()-1];
            paren.erase(paren.begin()+(paren.size()-1));
            close = true;

        }
        else if (consensus_structure[i] == '>' && !lts.empty()){
            x = lts[lts.size()-1];
            lts.erase(lts.begin()+(lts.size()-1));
            close = true;
        }
        else if (consensus_structure[i] == ']' && !sb.empty()){
            x = sb[sb.size()-1];
            sb.erase(sb.begin()+(sb.size()-1));
            close = true;
        }
        else if (consensus_structure[i] == '}' && !cb.empty()){
            x = cb[cb.size()-1];
            cb.erase(cb.begin()+(cb.size()-1));
            close = true;
        }
            
        if(close){
            if(i-std::get<1>(x) < 4){
                consensus_structure[i] = '.';
                consensus_structure[std::get<1>(x)] = '.';
                continue;  
            }
        }
    }

return consensus_structure;

}
bool canMatch(char x, char y){

if((x == 'A' && y == 'T') || (x == 'T' && y == 'A')) {return true;}
else if((x == 'C' && y == 'G') || (x == 'G' && y == 'C')) {return true;}
else if((x == 'A' && y == 'U') || (x == 'G' && y == 'U') || (x == 'U' && y == 'G') || (x == 'U' && y == 'A')) {return true;}
else{return false;}

}