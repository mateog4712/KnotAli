#include "../include/allFold.h"
#include "../include/mutual_information.h"
#include "../include/covariation.h"
#include "../include/utils.h"
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include "constants.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>


using namespace std;

char Simfold[10] = "./simfold";
int main(int argc, char *argv[]) {
    
    //int c;
    vector <string> list;
    string family;
    string sequence;
    string fastaFile;
    string clustalFile;
    string file;
    bool make = false;
    bool gap = false;
    bool perl = false;
    int c;

    opterr = 0;

    while ((c = getopt (argc, argv, "f:c:")) != -1)
        switch (c)
        {
            case 'f':
                fastaFile = optarg;
                break;
            case 'c':
                clustalFile = optarg;
                break; 

            default:
                abort ();
        }
    // Get the arguments
    if((fastaFile.length()!=0)){
        vector <string> seqs;
        vector <string> seqs2;
        vector <string> names;
        
        ifstream in(fastaFile.c_str());
        string str;
        int i = -1;
        bool newSeq = false;
        while (getline(in, str)) {

            if(str[0] == '>'){
                newSeq = true;
                names.push_back(str.substr(1,str.length()));
                continue;
            } 
            
            if(newSeq == true) {
                string str2 = str;
                for(int j=str2.length()-1;j>=0;--j){
                    if(str2[j] == '-') str2.erase(j,1);
                }
                seqs2.push_back(str2);
                seqs.push_back(str);
                i++;
                newSeq = false;
            }
            else{
               string str2 = str;
                for(int j=str2.length()-1;j>=0;--j){
                    if(str2[j] == '-') str2.erase(j,1);
                }
                string temp = seqs2[i] + str2;
                seqs2[i] = temp;
                string temp2 = seqs[i] + str;
                seqs[i] = temp2;
            }
        }
    in.close();
    
    string structure = MIVector(seqs);
    
    ofstream out("../results.fa", ofstream::trunc);
    for(int i=0; i<seqs.size(); ++i){
        
        string consensusCh = returnUngapped(seqs[i],structure);
        
        // run it
        string final = iterativeFold(seqs2[i],consensusCh);
        
        // makes sure name is in correct format
        if(names[i].substr(names[i].length()-4,4) == ".seq") names[i] = names[i].substr(0,names[i].length()-4);
        if(names[i].substr(names[i].length()-3,3) == ".ct") names[i] = names[i].substr(0,names[i].length()-3);
        istringstream ss(names[i]);
        ss >> names[i];
        

        /**   Outputting to file     **/
        out << ">" + names[i] << endl;
        out << seqs2[i] << endl;
        out << final << endl;
        out << endl;
    }
    out.close();

    
    }else if(clustalFile.length()!= 0){

    }
    return 0;
}