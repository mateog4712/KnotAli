#include "base_types.hh"
#include "mutual_information.hh"
#include "utils.hh"
#include "Iterative-HFold.hh"
#include "cmdline.hh"
#include "h_globals.hh"
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <unistd.h>
#include <sstream>

void seqtoRNA(std::string &sequence){
    for (char &c : sequence) {
      	if (c == 'T') c = 'U';
    }
}

void updateVectors(std::vector<std::string> & seqs, std::vector<std::string> &seqs2, std::vector<std::string> & names, const std::string & input_file, const std::string& input_type){
    std::string type = input_type;
    std::transform(type.begin(), type.end(), type.begin(), ::toupper);
    // Get the arguments
    if(type == "FASTA"){

        std::ifstream in(input_file.c_str());
        std::string str;
        cand_pos_t i = -1;
        bool newSeq = false;
        while (std::getline(in, str)) {

            if(str[0] == '>'){
                newSeq = true;
                names.push_back(str.substr(1,str.length()));
                continue;
            } 
            
            if(newSeq == true) {
                std::transform(str.begin(), str.end(), str.begin(), ::toupper);

                std::string str2 = str;
                cand_pos_t n = str2.length();
                for(cand_pos_t j=n-1;j>=0;--j){
                    if(str2[j] == '-') str2.erase(j,1);
                }
                seqs2.push_back(str2);
                seqs.push_back(str);
                i++;
                newSeq = false;
            }
            else{
                std::transform(str.begin(), str.end(), str.begin(), ::toupper);
                std::string str2 = str;
                cand_pos_t n = str2.length();
                for(cand_pos_t j=n-1;j>=0;--j){
                    if(str2[j] == '-') str2.erase(j,1);
                }
                std::string temp = seqs2[i] + str2;
                seqs2[i] = temp;
                std::string temp2 = seqs[i] + str;
                seqs[i] = temp2;
            }
        }
        in.close();

    }else if(type == "CLUSTAL"){
        
        std::ifstream in(input_file.c_str());
        std::string str;
        cand_pos_t i = 0;

        bool first = true;
        std::getline(in, str);
        if(str.substr(0,7) != "CLUSTAL"){
            std::cout << "Not a valid CLUSTAL file!" << std::endl;
            exit (EXIT_FAILURE);
        }
        std::getline(in, str);
        while (std::getline(in, str)) {
            if(names.size() == 0 && str == "") continue;

            if(str == "" || str[0] == ' '){
                first = false;
                i = 0;
                continue;
            }

            if(first){
                std::istringstream ss(str);
                ss >> str;
                names.push_back(str);
                ss >> str;
                std::transform(str.begin(), str.end(), str.begin(), ::toupper);
                seqs.push_back(str);
                cand_pos_t n = str.length();
                for(cand_pos_t j=n-1;j>=0;--j){
                    if(str[j] == '-') str.erase(j,1);
                }
                seqs2.push_back(str);

            }
            else{
                std::istringstream ss(str);
                ss >> str;
                ss >> str;
                std::transform(str.begin(), str.end(), str.begin(), ::toupper);
                seqs[i] = seqs[i] + str;
                cand_pos_t n = str.length();
                for(cand_pos_t j=n-1;j>=0;--j){
                    if(str[j] == '-') str.erase(j,1);
                }
                seqs2[i] = seqs2[i] + str;
                ++i;
            }
        }
        in.close();
    }
    else if(type == "STOCKHOLM"){
        std::ifstream in(input_file.c_str());
        std::string str;
        cand_pos_t i = 0;
        bool first = true;
        std::getline(in,str);
        getline(in,str);
        while(str.substr(0,4) == "#=GF"){
            getline(in,str);
        }
        while(getline(in,str)){

            if(str.substr(0,4) == "#=GC") break;
            if(str == "" || str[0] == ' '){
                first = false;
                i = 0;
                continue;
            }

            if(first){
                std::istringstream ss(str);
                ss >> str;
                names.push_back(str);
                ss >> str;
                std::transform(str.begin(), str.end(), str.begin(), ::toupper);
                seqs.push_back(str);
                cand_pos_t n = str.length();
                for(cand_pos_t j=n-1;j>=0;--j){
                    if(str[j] == '-') str.erase(j,1);
                }
                seqs2.push_back(str);

            }
            else{
                std::istringstream ss(str);
                ss >> str;
                ss >> str;
                std::transform(str.begin(), str.end(), str.begin(), ::toupper);
                seqs[i] = seqs[i] + str;
                cand_pos_t n = str.length();
                for(cand_pos_t j=n-1;j>=0;--j){
                    if(str[j] == '-') str.erase(j,1);
                }
                seqs2[i] = seqs2[i] + str;
                ++i;
            }
        }

    }
    else{
        std::cout << "Please give valid input type" << std::endl;
        exit (EXIT_FAILURE);
    }
    cand_pos_t n = seqs[0].length();
    cand_pos_t n_seq = seqs.size();
    for(cand_pos_t i=0;i< n_seq;++i){
        cand_pos_t n1 = seqs[i].length();
        if (n1 != n){
            std::cout << "All sequences must be the same length in the alignment" << std::endl;
            exit(0);
        }
    }
}


int main(int argc, char *argv[]) {


    args_info args_info;

	// get options (call gengetopt command line parser)
	if (cmdline_parser (argc, argv, &args_info) != 0) {
	exit(1);
	}

	std::string input_file;
	if (args_info.inputs_num>0) {
	input_file=args_info.inputs[0];
	} else {
	std::getline(std::cin,input_file);
	}
    bool verbose;
	verbose = args_info.verbose_given;

    // bool stacking;
    // stacking = args_info.stacking_given;

    std::string input_type;
    args_info.input_type_given ? input_type = type : input_type = "FASTA";

    std::string output_file_k;
    args_info.output_file_given ? output_file_k = output_file : output_file_k = "results.afa";

    std::string file= args_info.paramFile_given ? parameter_file : "params/rna_DirksPierce09.par";
	if(exists(file)){
		vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
	}

    int dangles = args_info.dangles_given ? dangle_model : 2;

    cmdline_parser_free(&args_info);

    if(!exists(input_file)){
        std::cout << "Input File does not exist!" << std::endl;
        exit (EXIT_FAILURE);
    }



    // Define the data structures for the sequences
    std::vector <std::string> seqs;
    std::vector <std::string> seqs2;
    std::vector <std::string> names;

    updateVectors(seqs,seqs2,names,input_file, input_type);

    // number of sequences
    int n_seq = seqs.size();

    // Converts any possible T's to U's in the sequences
    for(std::string &seq : seqs) seqtoRNA(seq);
    for(std::string &seq : seqs2) seqtoRNA(seq);

    // Uses covariation/ Mutual Information to find probable structurally important base pairs
    std::string structure = MIVector(seqs);

    if(verbose){
      printf ("\t The number of sequences are %d\n", n_seq);
      printf ("\t The structure found through covariation of the alignment is: \n\n%s\n", structure.c_str());  
    }
    
    // The output file
    std::ofstream out(output_file_k, std::ofstream::trunc);
    for(int i=0; i<n_seq; ++i){
    
        std::string consensus = returnUngapped(seqs[i],structure);
        std::string seq = removeIUPAC(seqs2[i]);
        double energy;        

        std::string final = Iterative_Fold(seq,consensus, energy, dangles);
    
        // makes sure name is in correct format
        if(names[i].substr(names[i].length()-4,4) == ".seq") names[i] = names[i].substr(0,names[i].length()-4);
        if(names[i].substr(names[i].length()-3,3) == ".ct") names[i] = names[i].substr(0,names[i].length()-3);
        std::istringstream ss(names[i]);
        ss >> names[i];

        /**   Outputting to file     **/
        out << ">" + names[i] << std::endl;
        out << seqs2[i] << std::endl;
        out << final << std::endl;
        out << energy << std::endl;
        out << std::endl;
    }
    out.close();
    
    return 0;
}