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

    while ((c = getopt (argc, argv, "s:f:c:mg:p:")) != -1)
        switch (c)
        {
            case 'g':
                gap = true;
                file = optarg;
                break;
            case 'p':
                perl = true;
                family = optarg;
                break;
            case 'm':
               make = true;
               break;
            case 's':
                sequence = optarg;
                break;
            case 'f':
                fastaFile = optarg;
                break;
            case 'c':
                clustalFile = optarg;
                break; 

            default:
                abort ();
        }
   

    //convertPS("../hxmatch-1.2.1/test.ps");
   
    if(make == true){
        ct2dot2(family, list);
    }
    else if(perl == true){
        runPerl(family);
    }
    else if(gap == true){
        /** Create the gapped version of the sequences **/
        ostringstream oss;
        oss << "../muscle3.8.31_i86linux64 -in " << file << " -out out1.afa"; 
        string command = oss.str();
        system(command.c_str());
    }
    // Get the arguments
    else if((fastaFile.length()!=0)){
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
    // Get the family for output info
    cout << "Put Family: ";
    cin >> family;

    

    // Define base variables
    int n = seqs[0].length();
    int n_seq = seqs.size();
    
    // string structure = MIVector(seqs);
    // cout << structure << endl;
    
    for(int i=0; i<seqs.size(); ++i){
        ostringstream oss;
        if(names[i].substr(names[i].length()-4,4) == ".seq") names[i] = names[i].substr(0,names[i].length()-4);
        if(names[i].substr(names[i].length()-3,3) == ".ct") names[i] = names[i].substr(0,names[i].length()-3);
        istringstream ss(names[i]);
        
        ss >> names[i];

        if((names[i].substr(0,family.length())) == family){
            
            }
        else{
            names[i] = (family + "_" + names[i]);
            }
        
        convertPS(family,names[i]);

        // oss << "../EMBOSS-6.6.0/emboss/seqret -sequence " << "../output/sequences/" << family << "/" << names[i] << ".txt " <<  "-outseq "  << "../output/clustal/" << family << "/" << names[i] << ".txt " <<  "-osformat2 clustal";
        // string command = oss.str();
        // system(command.c_str());

        // string file = "../output/clustal/" + family + "/" + names[i] + ".txt";
        // ifstream in(file);
        // //string line;
        // vector<string> lines;
        // int x = 0;
        // string extra = "                   *";

        // while (getline(in, str)) {
        //     // cout << "line: " << str << endl;
        //     if(x%3 == 1 && x!=1){
        //         lines.push_back(extra);
        //     }else{
        //         lines.push_back(str);
        //     }
        //     ++x;
        // }
        // in.close();
        // // cout << lines[0]<< endl;

        // ofstream out(file.c_str(),ofstream::trunc);
        // for(int j = 0; j<lines.size();++j){
        //     out << lines[j] << endl;
        // }
        // out.close();
        // string fileO = "../output/hxmatch/" + family + "/" + names[i] + ".txt";
        // string fileI = "../output/clustal/" + family + "/" + names[i] + ".txt";
        // oss << "../hxmatch-1.2.1/hxmatch < " << fileI << " > " << fileO;
        // string command = oss.str();
        // system(command.c_str());


    }

    
    
    // for(int i=0; i<seqs.size(); ++i){
        
    //     string consensusCh = returnUngapped(seqs[i],structure);
        
    //     // run it
    //     cout << i << endl;
    //     string final = consensusCh;//iterativeFold(seqs2[i],consensusCh);
    //     cout << final << endl;
        
    //     /**   Outputting to file     **/
    //     string File;
        
    //     if(names[i].substr(names[i].length()-4,4) == ".seq") names[i] = names[i].substr(0,names[i].length()-4);
    //     if(names[i].substr(names[i].length()-3,3) == ".ct") names[i] = names[i].substr(0,names[i].length()-3);
    //     istringstream ss(names[i]);
        
    //     ss >> names[i];
        
    //     if((names[i].substr(0,family.length())) == family){
    //         File = "../output/hfold/" + family + "/" + names[i] + ".txt";
            
    //         }
    //     else{
    //         File = "../output/hfold/" + family + "/" + family + "_" + names[i] + ".txt";
    //         }
    //     ofstream out(File.c_str(), ofstream::trunc);
    //     out << ">" + names[i] << endl;
    //     out << seqs2[i] << endl;
    //     out << final << endl;
    //     out.close();
    //     cout << endl;
    // }

    
    }else if(clustalFile.length()!= 0){

    }
    else if(sequence != ""){
        double energy = 0;
        int x = sequence.length();
        char output[x+1];
        char consens[x+1];
        strcpy(consens, sequence.c_str());
        // call simfold
    
        call_simfold2(Simfold,consens,output,&energy);
        // replace the . with _
        char* output2 = replaceChar(output,'.','_');
        string seqStr(output2);
        string final = iterativeFold(sequence,seqStr);
        cout << final << endl;
    }
    
    //string command = "export DATAPATH=~/RNAstructure/data_tables/";
    return 0;
}