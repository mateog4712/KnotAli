// a simple driver for the HFold
#include <iostream>
#include <algorithm>
#include <iterator>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <stack>
#include <ctime>
#include <dirent.h>
#include <signal.h>
#include <pthread.h>

// include the simfold header files
#include "simfold.h"
#include "externs.h"
#include "h_globals.h"
//#include "h_externs.h"
#include "constants.h"
#include "params.h"
#include "common.h"

#include <sys/types.h>
#include <sys/wait.h>
#include <string>

// Hosna June 20th, 2007
//#include "W_final.h"
//#include "hfold.h"
#include "hfold_iterative.h"

#include <ctime>
//#include <regex>
//#include <regex.h>

//kevin June 22 2017
#include "hfold_validation.h"
#include <getopt.h>

//#define HFOLD 						"./HFold"
//#define HFOLD_PKONLY 				"./HFold_pkonly"
//#define HFOLD_INTERACTING 			"./HFold_interacting"
//#define HFOLD_INTERACTING_PKONLY 	"./HFold_interacting_pkonly"
//#define SIMFOLD 					"./simfold"

// this causes less compilation warnings than the #defines
char HFOLD[8] =						"./HFold";
char HFOLD_PKONLY[15] =				"./HFold_pkonly";
char HFOLD_INTERACTING[20] =	    "./HFold_interacting";
char HFOLD_INTERACTING_PKONLY[27] =	"./HFold_interacting_pkonly";
char SIMFOLD[10] =                  "./simfold";

#define LOGFILEPATH 		"./logfile.txt"

FILE *logFile;
char *file;
int timeout = 1200;
//static std::smatch match;
//static std::regex sequenceRegex;
//static std::regex structureRegex;
//static std::regex replacedRegex;
static std::string sequenceString;
static std::string structureString;
static std::string replacementText;
static std::string energyString;
static std::stack<int> customStack;


/* I am changing iterative HFold so that we run 4 methods and choose the structure with the lowest energy as the winner
1) run HFold_PKonly on input structure, if pked, keep the pk bases as input structure and run HFold on them and get the result
2) run HFold on the original input and get the result
3) get the beginning and end of the given structure, give the subsequence and the substructure to simfold restricted as input, get the simfold structure, then give that as input to HFold_PKonly, get the pked bases as new input to give HFold and get the result
4) run simfold restricted with the given input sequence and structure, and get the simfold structure. Only choose part of the structure that contains the original given input structure, then give that to HFold_PKonly as input, get pked bases and run HFold on them and get the result. */




void printUsage(){
	/*
	printf ("\nUsage: HFold_iterative <sequence> <structure>\n");
	printf ("Example: ./HFold_iterative \"GCAACGAUGACAUACAUCGCUAGUCGACGC\" \"(____________________________)\" \n");
	printf ("Or \nUsage: HFold_iterative <path to input file> <path to output file>\n");
	printf ("Example: ./HFold_iterative \"/home/username/Desktop/inputFile.txt\" -o \"/home/username/Desktop/outFile.txt\" \n");
	printf ("\tRestricted structure symbols:\n");
	printf ("\t\t() restricted base pair\n");
	printf ("\t\t _ no restriction\n");
*/
	printf("Usage ./HFold_iterative --s <sequence> --r <structure> [--o </path/to/file>]\n");
	printf("or\n");
	printf("Usage ./HFold_iterative --i </path/to/file> [--o </path/to/file>]\n");
	printf ("  Restricted structure symbols:\n");
	printf ("    () restricted base pair\n");
	printf ("    _ no restriction\n");
	printf("Example:\n");
	printf("./HFold_iterative --s \"GCAACGAUGACAUACAUCGCUAGUCGACGC\" --r \"(____________________________)\"\n");
	printf("./HFold_iterative --i \"/home/username/Desktop/myinputfile.txt\" --o \"/home/username/Desktop/some_folder/outputfile.txt\"\n");
	printf("Please read README for more details\n");
}


void segfault_sigaction(int signal, siginfo_t *si, void *arg) {
	//Need to expand on this.
    fprintf(stderr, "Caught segfault at address %p\n", si->si_addr);
    exit(si->si_errno);
}



//Calling HFold: programPath = HFOLD
//Calling HFold_PKonly: programPath = HFOLD_PKONLY
bool call_HFold (char *programPath, char *input_sequence, char *input_structure, char *output_structure, double *output_energy) {
	char config_file[200];
	strcpy (config_file, SIMFOLD_HOME "/params/multirnafold.conf");

	//what to fold: RNA or DNA
	int dna_or_rna;
	dna_or_rna = RNA;
	//temperature: any integer or real number between 0 and 100
	// represents degrees Celsius
	double temperature = 37.0;

	// initialize the thermodynamic parameters
	// call init_data only once for the same dna_or_rna and same temperature
	// if one of them changes, call init_data again

	//init_data ("./HFold", config_file, dna_or_rna, temperature);
	init_data (programPath, config_file, dna_or_rna, temperature);

	fill_data_structures_with_new_parameters ( SIMFOLD_HOME "/params/turner_parameters_fm363_constrdangles.txt");

	// in HotKnots and ComputeEnergy package the most up-to-date parameters set is DP09.txt
	// so we add it here
	fill_data_structures_with_new_parameters ( SIMFOLD_HOME "/params/parameters_DP09.txt");

        if (strcmp(programPath, HFOLD) == 0) {
            *output_energy = hfold(input_sequence, input_structure, output_structure);
        }
        else if (strcmp(programPath, HFOLD_PKONLY) == 0){
                *output_energy = hfold_pkonly(input_sequence, input_structure, output_structure);
        }else{
                printf("Error: invalid arguments are given: %s \nValid aurgumnets are: HFOLD and HFOLD_PKONLY\n",programPath);
                return false;
         }
	if(is_invalid_restriction(input_structure,output_structure)){
		fprintf(stderr,"ERROR!!! There is something wrong with the structure, doesn't match restricted\n");
        fprintf(stderr,"  %s\n  %s\n  %s\t%.2lf\n", input_sequence, input_structure, output_structure, *output_energy);
        fprintf(stderr,"ERROR!!! There is something wrong with the structure, doesn't match restricted\n");
		exit(11);
	}
	return true;

}

bool call_simfold (char *programPath, char *input_sequence, char *input_structure, char *output_structure, double *output_energy) {
        

	char config_file[200] = SIMFOLD_HOME "/params/multirnafold.conf";

	double temperature;
	temperature = 37;
	init_data ("./simfold", config_file, RNA, temperature);

    fill_data_structures_with_new_parameters (SIMFOLD_HOME "/params/turner_parameters_fm363_constrdangles.txt");
	// when I fill the structures with DP09 parameters, I get a segmentation fault for 108 base sequence!!!!
	// So I chopped the parameter set to only hold the exact number as the turner_parameters_fm363_constrdangles.txt,
	// but still getting seg fault!
	fill_data_structures_with_new_parameters (SIMFOLD_HOME "/params/parameters_DP09_chopped.txt");

	*output_energy = simfold_restricted (input_sequence, input_structure, output_structure);
	//*output_energy = simfold_restricted (input_sequence, output_structure);
//	printf ("Call_Simfold_RES( can be called by different methods): %s  %.2lf\n", output_structure, output_energy);
	if(is_invalid_restriction(input_structure,output_structure)){
			fprintf(stderr,"ERROR!!! There is something wrong with the structure, doesn't match restricted\n");
			fprintf(stderr,"  %s\n  %s\n  %s\t%.2lf\n", input_sequence, input_structure, output_structure, *output_energy);
			fprintf(stderr,"ERROR!!! There is something wrong with the structure, doesn't match restricted\n");
			exit(11);
	}
	return true;
}


//kevin 28 Aug
// Function to remove all spaces from a given string
void removeSpaces(char *str)
{
    // To keep track of non-space character count
    int count = 0;

    // Traverse the given string. If current character
    // is not space, then place it at index 'count++'
    for (int i = 0; str[i]; i++)
        if (str[i] != ' ' || str[i] != '\n' || str[i] != '\t' || str[i] != '\r')
            str[count++] = str[i]; // here count is incremented
    str[count-1] = '\0';
}


bool get_sequence_structure (char *fileName, char *sequence, char *structure) {
        FILE *ioFile;
        size_t len = 0;
        ssize_t read = 0;
        char *sequenceBuffer = NULL;
        char *structureBuffer = NULL;

        //replacementText = "";
        //sequenceRegex = "(\\s|\n|\r)";
        //structureRegex = "(\\s|\n|\r|\Z)";

        // Get structure and sequence from file
        ioFile = fopen(fileName, "r");
        if (ioFile != NULL) {
                if((read = getline(&sequenceBuffer, &len, ioFile)) == -1) {
                        write_log_file("Could not read the first sequence", fileName, 'E');
                        return false;
                }
                //sequenceString.assign(sequenceBuffer, read);

                if((read = getline(&structureBuffer, &len, ioFile)) == -1) {
                        write_log_file("Could not read the first structure", fileName, 'E');
                        return false;
                }
                //structureString.assign(structureBuffer, read);

        } else {
                write_log_file("Could not open file", fileName, 'E');
                return false;
        }

		//kevin 28 Aug 2017
        //added removeSpaces to replace the following regex things
        removeSpaces(sequenceBuffer);
        //printf("new: %s \n",sequenceBuffer);
        removeSpaces(structureBuffer);
        //printf("new: %s \n",structureBuffer);
        strcpy(sequence, sequenceBuffer);
        strcpy(structure, structureBuffer);

/*
        sequenceString = std::regex_replace(sequenceString, sequenceRegex, replacementText);
        std::cout << "Sequence:  " << sequenceString << " Length: " << sequenceString.length() << '\n' << std::flush;

        // Format structure
        structureString = std::regex_replace(structureString, structureRegex, replacementText);
        std::cout << "Structure: " << structureString << " Length: " << structureString.length() << '\n' << std::flush;

        // Copy strings
        strcpy(sequence, sequenceString.c_str());
        strcpy(structure, structureString.c_str());
*/
        // Clean up
        fclose(ioFile);
        if (sequenceBuffer)
        free(sequenceBuffer);
        if (structureBuffer)
        free(structureBuffer);
/*
        if(sequenceString.length() != structureString.length()) {
                write_log_file("Sequence and structure length do not match", fileName, 'E');
                return false;
        }
*/
        //kevin 28 Aug 2017
        //changed from sequenceString.length to strlen(sequence)
        if(strlen(sequence) != strlen(structure)) {
                write_log_file("Sequence and structure length do not match", fileName, 'E');
                return false;
        }

        return true;
}


bool save_file (const char *fileName, char *outputPath, const char *sequence, char *restricted, char *structure, double energy, int chosen_method) {
	FILE *ioFile;
	char filePath[MAXSLEN];
	strcpy(filePath, outputPath);
	strcat(filePath, fileName);

	std::cout << "Output File: " << filePath << '\n' << std::flush;

	ioFile = fopen(filePath, "w");
	if (ioFile != NULL) {
		fprintf(ioFile, "Sequence: %s\nInput_structure: %s\nOutput_structure: %s\nEnergy: %.2f\nMethod: %i", sequence, restricted, structure, energy, chosen_method);
	} else {
		write_log_file("Could not open file", fileName, 'E');
		perror("Write to file error");
		return false;
	}

	// Clean up
	fclose(ioFile);

	return true;
}

void write_log_file(const char *message, const char *fileName, const char option) {
	time_t now = time(0);
	struct tm tstruct = *localtime(&now);
	char buf[80];

	strftime(buf, sizeof(buf), "%c", &tstruct);

	switch (option) {
		case 'i':
		case 'I':
			fprintf(logFile, "[%s] Info: %s\n", buf, message);
			break;

		case 'a':
		case 'A':
			fprintf(logFile, "[%s] Action: %s in file %s\n", buf, message, fileName);
			break;

		case 'e':
		case 'E':
			fprintf(logFile, "[%s] Error: %s in file %s\n", buf, message, fileName);
			break;
	}

	fflush(logFile);
}



//---------------------------------------this function is suppose to be the same as the one in Hfold_interacting, if any changes are made, please change that one too--------------------
//kevin 19 july (verified by Mahyar 19 july 2017)
//count+1 when open bracket, count-1 when close bracket
//whgen count is 0, we have a substructure
//assume valid structure
void find_disjoint_substructure(char* structure, std::vector< std::pair<int,int> > &pair_vector){
	int length = strlen(structure);
	int count = 0;
	int first_time = 1; //flag for first time getting a open bracket for substructure
	int i = 0;
	int j = 0;
	for(int k=0; k<length;k++){
		if(structure[k] == '(' || structure[k] == '['){
			if(first_time && count == 0){
				first_time = 0;
				i = k;
			}
			count += 1;

		}else if(structure[k] == ')' || structure[k] == ']'){
			count -= 1;
			j = k;
			if(count == 0){
				std::pair <int,int> ij_pair (i,j);
				pair_vector.push_back(ij_pair);
				first_time = 1;
			}
		}
	}
}

//---------------------------------------this function is suppose to be the same as the one in Hfold_interacting, if any changes are made, please change that one too--------------------
//31 Aug 2017 kevin and Mahyar
//input[i] is _ or . and it did not turn into a . in the output structure, then it is not empty
int is_empty_structure(char* input_structure, char* output_structure){
	for(int i=0; i<strlen(input_structure);i++){
		if(input_structure[i] != output_structure[i]){
			if((input_structure[i] == '_' || input_structure[i] == '.' ) && output_structure[i] != '.'){
				return 0;
			}
		}
	}
	return 1;
}


//---------------------------------------this function is suppose to be the same as the one in Hfold_interacting, if any changes are made, please change that one too--------------------
//kevin 18 July
int paired_structure(int i, int j, int *pair_index, int length){
	if(i >= 0 && j < length && (pair_index[i] == j) && (pair_index[j] == i) ){
		return 1;
	}
	return 0;
}

//---------------------------------------this function is suppose to be the same as the one in Hfold_interacting, if any changes are made, please change that one too--------------------
//kevin 18 July
void obtainRelaxedStems(char* G1, char* G2, char* Gresult){
	int length = strlen(G1);
	int G1_pair[length];
	int G2_pair[length];

	//Gresult <- G1
	strcpy(Gresult,G1);

	detect_original_pairs(G1, G1_pair);
	detect_original_pairs(G2, G2_pair);

	//for(int d=0;d<length;d++){
	//	printf("%c %d\n",G2[d],G2_pair[d]);
	//}

	int i = 0;
	int j = 0;

	for(int k=0;k<length;k++){
		if(G2_pair[k] > -1){
			i = k;
			j = G2_pair[k];
			if(i < j){ //for each ij in G2
				if( (G1[i] != G2[i]) && (G1[j] != G2[j]) ){//if ij not in G1
					//include bulges of size 1
					if(paired_structure(i-1,j+1,G1_pair,length) || paired_structure(i+1,j-1,G1_pair,length) ){
						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
					//include loops of size 1x1
					}else if( paired_structure(i-2,j+1,G1_pair,length) || paired_structure(i-1,j+2,G1_pair,length) || \
							paired_structure(i+1,j-2,G1_pair,length) || paired_structure(i+2,j-1,G1_pair,length) ){
						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
					//include loops of size 1x2 or 2x1
					}else if( paired_structure(i-2,j+2,G1_pair,length) || paired_structure(i+2,j-2,G1_pair,length) ){
						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
					}else if( paired_structure(i-3,j+2,G1_pair,length) || paired_structure(i-2,j+3,G1_pair,length) || \
							paired_structure(i+2,j-3,G1_pair,length) || paired_structure(i+3,j-2,G1_pair,length) ){

						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
					}
				}
			}
		}
	}
}




//---------------------------------------this function is suppose to be the same as the one in Hfold_interacting, if any changes are made, please change that one too--------------------
//kevin 30 Aug 2017
//check if the computed structure matches the restricted structure
int is_invalid_restriction(char* restricted_structure, char* current_structure){
	std::string openBracketArray ("({[");
	std::string closeBracketArray (")}]");

	for (int i=0; i < strlen(restricted_structure); i++){
        if(restricted_structure[i] != '_' && restricted_structure[i] != current_structure[i]){
			if( (openBracketArray.find_first_of(restricted_structure[i]) != -1) && ((openBracketArray.find_first_of(current_structure[i]) != -1)) ){
				continue;
			}else if ( (closeBracketArray.find_first_of(restricted_structure[i]) != -1) && ((closeBracketArray.find_first_of(current_structure[i]) != -1)) ){
				continue;
			}else{
				return 1;
			}
		}

    }
	return 0;
}


//30 Aug 2017 kevin and Mahyar
double method1(char *sequence, char *restricted, char *structure){
	//printf("inside method1 restricrted:\n%s\n",restricted);
	double energy = 0;
	call_HFold(HFOLD, sequence, restricted, structure, &energy);
	//printf("method1: %s\n",structure);
	return energy;
}


//30 Aug 2017 kevin and Mahyar
double method2(char *sequence, char *restricted, char *structure){
	//printf("inside method2 restricrted:\n%s\n",restricted);
	double energy = 0;

	call_HFold(HFOLD_PKONLY, sequence, restricted, structure, &energy);
	//printf("inside method2 structure:\n%s\n",structure);
	//printf("restricted: %s\nstructure: %s\n",restricted,structure);
	if(is_empty_structure(restricted,structure)){
		//printf("is empty\n");
		return energy;
	}else{
		//printf("is not empty\n");
		char G_prime[strlen(structure)];
		remove_structure_intersection(structure,restricted, G_prime);
		//printf("G_prime: %s\n",G_prime);
		energy = method1(sequence, G_prime, structure);
		return energy;
	}

}

//30 Aug 2017 kevin and Mahyar
double method3(char *sequence, char *restricted, char *structure){
	double energy = 0;
	int length = strlen(sequence);
	char simfold_structure[length];

	call_simfold(SIMFOLD, sequence, restricted, simfold_structure, &energy);

	//^ G' simfold_structure <- SimFold(S sequence, G restricted)
	char G_updated[length+1];
	G_updated[length] = '\0';
	// std::cout << "restricted is: " << restricted << std::endl;
	obtainRelaxedStems(restricted ,simfold_structure, G_updated);
	// printf("the structure is: %s\n",G_updated);
	//^Gupdated G_updated<- ObtainRelaxedStems(G restricted,G' simfold_structure)
	energy = method2(sequence, G_updated, structure);

	return energy;
}

//30 Aug 2017 kevin and Mahyar
double method4(char *sequence, char *restricted, char *structure){
	int KEVIN_DEBUG = 0;
	double energy = 0;
	int length = strlen(sequence);
	char G_updated[length+1];
	int k = 1;
	//^k <- 1
	strcpy(G_updated, restricted);
	//^Gupdated <- G
	std::vector< std::pair<int,int> > disjoint_substructure_index; //contain pair(i,j) in a vector, ij are index of begin and end of substructure
	find_disjoint_substructure(restricted, disjoint_substructure_index);
	//^get disjoint substructure
	int i = 0;
	int j = 0;

	//31 Aug 2017 kevin and Mahyar
	if(disjoint_substructure_index.size() == 0){
		return INF;
	}

	for(auto current_substructure_index : disjoint_substructure_index){
		i = current_substructure_index.first;
		j = current_substructure_index.second;
		char subsequence[length+1];
		char substructure[length+1];
		char simfold_structure[length+1];

		strncpy(subsequence, sequence+i,j-i+1);
		subsequence[j-i+1] = '\0';
		//^Sk
		strncpy(substructure, restricted+i,j-i+1);
		substructure[j-i+1] = '\0';
		//^Gk

    
		call_simfold(SIMFOLD, subsequence, substructure, simfold_structure, &energy);

		//^ SimFold(Sk,Gk,Gk')
		char Gp_k_updated[length];
		obtainRelaxedStems(substructure, simfold_structure, Gp_k_updated);
		//^obtainRelaxedStems(Gk,Gk',G'kupdated)

		int m = 0; //index for going through Gp_k_updated
        for(int k =i;k<j;k++){
			if(G_updated[k] != Gp_k_updated[m]){
				G_updated[k] = Gp_k_updated[m];
			}
			m++;
		}
		//^Gupdated <- Gupdated U G'kupdated
	}
	energy = method2(sequence, G_updated, structure);

	return energy;
}

//---------------------------------------this function is suppose to be the same as the one in Hfold_interacting, if any changes are made, please change that one too--------------------
// Aug 31, 2017 kevin and Mahyar
//does G_p = G1-G
void remove_structure_intersection(char* G1, char* G0, char* G_p){
	strcpy(G_p,G1);
	for(int i=0; i< strlen(G1); i++){
		if(G_p[i] == '.'){
			G_p[i] = '_';
		}
		if (G1[i] == G0[i]){
			G_p[i] = '_';
		}
		if (G_p[i] == '[' || G_p[i] == '{'){
			G_p[i] = '(';
		}
		if (G_p[i] == ']' || G_p[i] == '}'){
			G_p[i] = ')';
		}
	}
}
