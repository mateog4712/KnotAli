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
    string type;
    string fastaFile;
    string clustalFile;
    int c;

    opterr = 0;

    while ((c = getopt (argc, argv, "f:c:o:")) != -1)
        switch (c)
        {
            case 'f':
                fastaFile = optarg;
                break;
            case 'c':
                clustalFile = optarg;
                break;
            case 'o':
                type = optarg;
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
    
    ofstream out("results.fa", ofstream::trunc);
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
        vector <string> seqs;
        vector <string> seqs2;
        vector <string> names;
        
        ifstream in(clustalFile.c_str());
        string str;
        int i = 0;

        bool first = true;
        getline(in, str);
        if(str.substr(0,7) != "CLUSTAL"){
            cout << "Not a valid CLUSTAL file!" << endl;
            exit (EXIT_FAILURE);
        }
        getline(in, str);
        while (getline(in, str)) {
            if(names.size() == 0 && str == "") continue;

            if(str == "" || str[0] == ' '){
                first = false;
                i = 0;
                continue;
            }

            if(first){
                istringstream ss(str);
                ss >> str;
                names.push_back(str);
                ss >> str;
                seqs.push_back(str);
                for(int j=str.length()-1;j>=0;--j){
                    if(str[j] == '-') str.erase(j,1);
                }
                seqs2.push_back(str);

            }
            else{
                istringstream ss(str);
                ss >> str;
                ss >> str;
                seqs[i] = seqs[i] + str;
                for(int j=str.length()-1;j>=0;--j){
                    if(str[j] == '-') str.erase(j,1);
                }
                seqs2[i] = seqs2[i] + str;
                ++i;
            }
        }
    in.close();
    
    string structure = MIVector(seqs);
    
    ofstream out("results.fa", ofstream::trunc);
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
    
    
    }
    return 0;
}