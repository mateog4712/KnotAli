#include "utils.hh"
#include "Iterative/HFold_iterative.cpp"
#include <sys/stat.h>
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <vector>
#include <tuple>


bool exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}


bool canPair(int n){
    if(n == 9 || n == 13 || n == 17 || n == 19 || n == 21 || n == 23) return true;
    return false;
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

    int sublength=0;
    int j = length;
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
            paren.erase(paren.end());
            close = true;

        }
        else if (consensus_structure[i] == '>' && !lts.empty()){
            x = lts[lts.size()-1];
            lts.erase(lts.end());
            close = true;
        }
        else if (consensus_structure[i] == ']' && !sb.empty()){
            x = sb[sb.size()-1];
            sb.erase(sb.end());
            close = true;
        }
        else if (consensus_structure[i] == '}' && !cb.empty()){
            x = cb[cb.size()-1];
            cb.erase(cb.end());
            close = true;
        }
            
        // Checks to see if the pair can match
        if(close){
            if(canMatch(input_sequence[i],input_sequence[std::get<1>(x)])){
                
                continue;
            }
            else{
                if(input_sequence[i] == '-'){
                    consensus_structure[std::get<1>(x)] = '_';
                    continue;
                }else if(input_sequence[std::get<1>(x)] == '-'){
                    consensus_structure[i] = '_';
                    continue;   
                }
                else{
                    consensus_structure[i] = '_';
                    consensus_structure[std::get<1>(x)] = '_';
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
            paren.erase(paren.end());
            close = true;

        }
        else if (consensus_structure[i] == '>' && !lts.empty()){
            x = lts[lts.size()-1];
            lts.erase(lts.end());
            close = true;
        }
        else if (consensus_structure[i] == ']' && !sb.empty()){
            x = sb[sb.size()-1];
            sb.erase(sb.end());
            close = true;
        }
        else if (consensus_structure[i] == '}' && !cb.empty()){
            x = cb[cb.size()-1];
            cb.erase(cb.end());
            close = true;
        }
            
        if(close){
            if(i-std::get<1>(x) < 4){
                consensus_structure[i] = '_';
                consensus_structure[std::get<1>(x)] = '_';
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

std::string iterativeFold(std::string seq, std::string str, double &en){

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

    en = final_energy;

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