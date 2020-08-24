#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <stack>
#include <ctime>
#include <signal.h>
#include <pthread.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <sstream>
#include <fstream>
#include <vector>
#include <iterator>
#include <cstdlib>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
//#include "simfold.h"
#include "externs.h"
#include "../h_globals.h"
#include "constants.h"
#include "params.h"
//#include "hfold.h"
#include "../hfold_iterative.h"
#include "../HFold_iterative.cpp"
#include <assert.h>
#include "../include/allFold.h"
#include "../include/utils.h"
//#include "include/mutual_information.h"
#include <tuple>



using namespace std;

/** Takes in the consensus sequence and structure and the original sequence and returns the ungapped version 
of the structure to be used in simfold **/
string returnUngapped(string input_sequence, string consensus_structure){
int length = consensus_structure.length();

int sublength=0;
int j = length;
vector<tuple<char,int> > s;
for(int i=0;i<length;i++){
    if(consensus_structure[i] == '(') {
        s.push_back(make_tuple(consensus_structure[i],i));
        continue;
    }
    if (consensus_structure[i] == ')'){
        tuple<char,int> x = s[s.size()-1];
        s.erase(s.end());
        if(get<0>(x) == '('){
            if(canMatch(input_sequence[i],input_sequence[get<1>(x)])){
                
                continue;
            }
            else{
                if(input_sequence[i] == '-'){
                    consensus_structure[get<1>(x)] = '_';
                    continue;
                }else if(input_sequence[get<1>(x)] == '-'){
                    consensus_structure[i] = '_';
                    continue;   
                }
                else{
                    consensus_structure[i] = '_';
                    consensus_structure[get<1>(x)] = '_';
                    continue;
                }
            }
        }
    }
}
for(int i= input_sequence.length()-1; i>=0;--i){
    if(input_sequence[i] == '-') consensus_structure.erase(i,1);
}
return consensus_structure;
}
bool canMatch(char x, char y){

if((x == 'A' && y == 'T') || (x == 'T' && y == 'A')) {return true;}
else if((x == 'C' && y == 'G') || (x == 'G' && y == 'C')) {return true;}
else if((x == 'A' && y == 'U') || (x == 'G' && y == 'U') || (x == 'U' && y == 'G') || (x == 'U' && y == 'A')) {return true;}
else{return false;}

}



void convertPS(string family, string name){
string fileI = "../output/hxmatch/" + family + "/" + name + ".txt";
string fileO = "../output/hxbp/" + family + "/" + name + ".bpseq";
string fileSeq = "../output/sequences/" + family + "/" + name + ".txt";

string seq;
string str;
ifstream in(fileSeq);
getline(in, str);
getline(in, str);
seq = str;
in.close();

vector< tuple<int,int> > pairs;
ifstream in1(fileI);
getline(in1,str);
getline(in1,str);
while(getline(in1,str)){
istringstream ss(str);
int i;
ss >> i;
int j;
ss >> j;
pairs.push_back(make_tuple(i,j));
}
in1.close();
int j = 0;
ofstream out(fileO, ofstream::trunc);
for(int i=0;i<seq.length();++i){
if(get<0>(pairs[j]) == i) {
    out << (i+1) << " " << seq[i] << " " << (get<1>(pairs[j])+1) << endl;
    ++j;
}
else{
    out << (i+1) << " " << seq[i] << " " << "0" << endl;  
}

}
out.close();
}
bool call_simfold2 (char *programPath, char *input_sequence, char *output_structure, double *output_energy) {
        

	char config_file[200] = SIMFOLD_HOME "/params/multirnafold.conf";

	double temperature;
	temperature = 37;
	init_data ("./simfold", config_file, RNA, temperature);

    fill_data_structures_with_new_parameters (SIMFOLD_HOME "/params/turner_parameters_fm363_constrdangles.txt");
	// when I fill the structures with DP09 parameters, I get a segmentation fault for 108 base sequence!!!!
	// So I chopped the parameter set to only hold the exact number as the turner_parameters_fm363_constrdangles.txt,
	// but still getting seg fault!
	fill_data_structures_with_new_parameters (SIMFOLD_HOME "/params/parameters_DP09_chopped.txt");

	*output_energy = simfold (input_sequence, output_structure);
	//*output_energy = simfold_restricted (input_sequence, output_structure);
//	printf ("Call_Simfold_RES( can be called by different methods): %s  %.2lf\n", output_structure, output_energy);
	return true;
}
bool call_simfold3 (char *programPath, char *input_sequence, char *output_structure, double *output_energy, double *scores ,int n) {
        

	char config_file[200] = SIMFOLD_HOME "/params/multirnafold.conf";

	double temperature;
	temperature = 37;
	init_data ("./simfold", config_file, RNA, temperature);

    fill_data_structures_with_new_parameters (SIMFOLD_HOME "/params/turner_parameters_fm363_constrdangles.txt");
	// when I fill the structures with DP09 parameters, I get a segmentation fault for 108 base sequence!!!!
	// So I chopped the parameter set to only hold the exact number as the turner_parameters_fm363_constrdangles.txt,
	// but still getting seg fault!
	fill_data_structures_with_new_parameters (SIMFOLD_HOME "/params/parameters_DP09_chopped.txt");

	*output_energy = simfold_new (input_sequence, output_structure,scores,n);
	//*output_energy = simfold_restricted (input_sequence, output_structure);
//	printf ("Call_Simfold_RES( can be called by different methods): %s  %.2lf\n", output_structure, output_energy);
	return true;
}

string get_consensus(vector<string> vects, bool mis){
    int n = vects[0].length();
    int n_seq = vects.size();
    string consensus = "";
    if(mis){
        string IUP2 = "-ACCGGCGUUCCGGCG";
        string IUP = "-ACMGRSVUWYHKDBN";
        vector<int> bgfreq = {0,0,0,0,0};
        for(int i=0;i<n;i++){
            for(int s =0; s<n_seq;s++){
                if(vects[s][i] == '-') bgfreq[0]++;
                if(vects[s][i] == 'A') bgfreq[1]++;
                if(vects[s][i] == 'C') bgfreq[2]++;
                if(vects[s][i] == 'G') bgfreq[3]++;
                if(vects[s][i] == 'U') bgfreq[4]++;
            }
        }

        for(int i=0;i<n;i++){
            int code  = 0; 
            vector<int> freq = {0,0,0,0,0};
            for(int s =0; s<n_seq;s++){
                if(vects[s][i] == '-') freq[0]++;
                if(vects[s][i] == 'A') freq[1]++;
                if(vects[s][i] == 'C') freq[2]++;
                if(vects[s][i] == 'G') freq[3]++;
                if(vects[s][i] == 'U') freq[4]++;
            }
            for (int c = 4; c > 0; c--) {
                code <<= 1;
                if (freq[c] * n >= bgfreq[c]) code++;
            }
            consensus += IUP2[code];
            //if (freq[0] * n > bgfreq[0]) consensus[i] = tolower(IUP2[code]);
        }
    } else{
        string basic = "-ACGU";
        unsigned int s;
        for(int i=0;i<n;i++){
            vector<int> freq = {0,0,0,0,0};
            for(int s =0; s<n_seq;s++){
                if(vects[s][i] == '-') freq[0]++;
                if(vects[s][i] == 'A') freq[1]++;
                if(vects[s][i] == 'C') freq[2]++;
                if(vects[s][i] == 'G') freq[3]++;
                if(vects[s][i] == 'U') freq[4]++;
            }
            int c, fm; 
            for (s = c = fm = 0; s < 5; s++) {/* find the most frequent char */
                if (freq[s] > fm) {
                    c   = s;
                    fm  = freq[c];
                }
            }
            consensus += basic[c];
        }
    }
return consensus;
}






void runPerl(string family) {
/** Updates the perl file for BPseq and runs it on the files from the family inputted **/
    string line;
    string folderI;
    string folderO;
    string dirN;

    // Runs the following parts twice: once for predictions and once for the correct info
    for (int i = 0; i < 2; i++) {
        // Takes the converted perl file and opens it for reading and a new perl file for writing
        // When done as a perl file, the getline function messes up and skips lines
        string bpseq = "../bpseq.txt";
        string bpseqU = "../bpseqU.pl";
        ifstream bp(bpseq.c_str());
        ofstream bpU(bpseqU.c_str(), ofstream::trunc);

        // defines the folders that are going to be outputted to
        if (i == 0) {
            folderI = "hfold";
            folderO = "predbp";
        } else {
            folderI = "rnastructure";
            folderO = "corrbp";
        }
        // Creates the folders if there are not currently there
        dirN = "../output/" + folderO + "/" + family + "/";
        struct stat info;
        if (stat(dirN.c_str(), &info) != 0) {
            system(("mkdir " + dirN).c_str());
        }

        // Starts a counter for reading the file
        int counter = 1;
        while (getline(bp, line)) {
            // If we've reached the 7th line of the txt file, place the line which look at the correct folder
            if (counter == 7) {
                bpU << "$inputpath = \"/home/mgray7/output/" << folderI << "/" << family << "/\";" << "#input path" << endl;
            } // The same is done with the 8th line. It is replaced with the correct folder location
            else if (counter == 8) {
                bpU << "$outputpath = \"/home/mgray7/output/" << folderO << "/" << family << "/\";" << "#output path"
                    << endl;
            } else {
                // Send the exact line from the txt file into the perl file
                bpU << line << endl;
            }
            // Increment the counter
            counter++;
        }
        // Close the files
        bp.close();
        bpU.close();
        // Run the updated perl file
        string command = "perl ../bpseqU.pl";
        // system(command.c_str());
    }

    string calcf = "../calculate_f_measure.txt";
    string calcfU = "../calculate_f_measureU.pl";
    ifstream cf(calcf.c_str());
    ofstream cfU(calcfU.c_str(), ofstream::trunc);

    // Starts a counter for reading the file
    int counter = 1;
    while (getline(cf, line)) {
        // If we've reached the 7th line of the txt file, place the line which look at the correct folder
        if (counter == 98) {
            cfU << "$crr_path = \"/home/mgray7/output/corrbp/" << family << "/\";" << endl;
        } // The same is done with the 8th line. It is replaced with the correct folder location
        else if (counter == 99) {
            //change this line to fix problem
            cfU << "$chk_path = \"/home/mgray7/output/hxbp/" << family << "/\";" << endl;
        } else if (counter == 152) {
            cfU << "$experiment = \"/home/mgray7/output/results/" << family << ".txt\";" << endl;
        } else {
            // Send the exact line from the txt file into the perl file
            cfU << line << endl;
        }
        // Increment the counter
        counter++;

    }
    // Closes the files
    cf.close();
    cfU.close();

    // Runs the updated perl file
    string command = "perl ../calculate_f_measureU.pl";
    system(command.c_str());

}



/** This function takes a set of sequences and find the consensus structure and sequence for it **/
void rnaalifold(string sequences, string outputFile, string &structure) {
    ostringstream oss;
    // Takes the sequence and runs it through RNAalifold
    oss << "./" << "ViennaRNA/bin/" << "RNAalifold   " << sequences << " > " << outputFile;
    string command = oss.str();
    system(command.c_str());

    // Gets the content of the file that RNAalifold outputted to and takes the structure from the file and saves it to a variable
    // The variable is returned to used later through a pointer
    vector <string> lines;
    getFileContent(outputFile, lines);
    structure = lines[1];

    // Does something but im unsure what currently
    istringstream iss(structure);
    vector <string> results((istream_iterator<string>(iss)), istream_iterator<string>());
    structure = results[0];

}

// Runs iterative HFold
string iterativeFold(string seq, string str){

    void *res;

    char sequence[seq.length()+1];
    char structure[str.length()+1];
    strcpy(sequence, seq.c_str());
    strcpy(structure, str.c_str());

    char *output_path;
    char *method1_structure = (char*) malloc(sizeof(char) * MAXSLEN);
    char *method2_structure = (char*) malloc(sizeof(char) * MAXSLEN);
    char *method3_structure = (char*) malloc(sizeof(char) * MAXSLEN);
    char *method4_structure = (char*) malloc(sizeof(char) * MAXSLEN);
    char final_structure[MAXSLEN];

    double *method1_energy = (double*) malloc(sizeof(double) * INF);
    double *method2_energy = (double*) malloc(sizeof(double) * INF);
    double *method3_energy = (double*) malloc(sizeof(double) * INF);
    double *method4_energy = (double*) malloc(sizeof(double) * INF);
    double final_energy = INF;
    int method_chosen = -1;

    *method1_energy = INF;
    *method2_energy = INF;
    *method3_energy = INF;
    *method4_energy = INF;

    method1_structure[0] = '\0';
    method2_structure[0] = '\0';
    method3_structure[0] = '\0';
    method4_structure[0] = '\0';
    final_structure[0] = '\0';

    //printf("method1\n");
    *method1_energy = method1(sequence, structure, method1_structure);
    //printf("method2\n");
    *method2_energy = method2(sequence, structure, method2_structure);
    //printf("method3\n");
    

    *method3_energy = method3(sequence, structure, method3_structure);
    //printf("method4\n");
    *method4_energy = method4(sequence, structure, method4_structure);
  
  
  //double energy = 0;
	//int length = strlen(sequence);
	//char simfold_structure[length];
  //char restricted[length];
  //strcpy(restricted,structure);
  //call_simfold(SIMFOLD, sequence, restricted, simfold_structure, &energy);
  


    //We ignore non-negetive energy, only if the energy of the input sequnces are non-positive!
    if (*method1_energy < final_energy) {
        //if (*method1_energy < final_energy && *method1_energy != 0) {
        final_energy = *method1_energy;
        strcpy(final_structure, method1_structure);
        method_chosen = 1;
    }

    if (*method2_energy < final_energy) {
        //if (*method2_energy < final_energy && *method2_energy != 0) {
        final_energy = *method2_energy;
        strcpy(final_structure, method2_structure);
        method_chosen = 2;
    }

    if (*method3_energy < final_energy) {
        //if (*method3_energy < final_energy && *method3_energy != 0) {
        final_energy = *method3_energy;
        strcpy(final_structure, method3_structure);
        method_chosen = 3;
    }

    if (*method4_energy < final_energy) {
        //if (*method4_energy < final_energy && *method4_energy != 0) {
        final_energy = *method4_energy;
        strcpy(final_structure, method4_structure);
        method_chosen = 4;
    } 


    if (final_energy == INF || method_chosen == -1) {
        fprintf(stderr, "ERROR: could not find energy\n");
        fprintf(stderr, "SEQ: %s\n",sequence);
        fprintf(stderr, "Structure: %s\n",structure);
    }
    //cout << sequence << endl << final_structure << endl << final_energy;
    free(method1_energy); free(method1_structure);
    free(method2_energy); free(method2_structure);
    free(method3_energy); free(method3_structure);
    free(method4_energy); free(method4_structure);
    return final_structure;
    

}


