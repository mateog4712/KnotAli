#include "mutual_information.hh"
#include "utils.hh"
#include "cmdline.hh"
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <unistd.h>

using namespace std;

int main(int argc, char *argv[]) {


    args_info args_info;

	// get options (call gengetopt command line parser)
	if (cmdline_parser (argc, argv, &args_info) != 0) {
	exit(1);
	}

	string input_file;
	if (args_info.inputs_num>0) {
	input_file=args_info.inputs[0];
	} else {
	getline(cin,input_file);
	}
    

    bool verbose;
	verbose = args_info.verbose_given;

    bool stacking;
    stacking = args_info.stacking_given;

    string input_type;
    args_info.input_type_given ? input_type = type : input_type = "FASTA";

    int threads;
    args_info.threads_given ? threads = numThreads : threads = 1;

    string output_file;
    args_info.output_file_given ? output_file = file : output_file = "results.afa";

    cmdline_parser_free(&args_info);

    if(!exists(input_file)){
        cout << "Input File does not exist!" << endl;
        exit (EXIT_FAILURE);
    }

    

    // Define the data structures for the sequences
    vector <string> seqs;
    vector <string> seqs2;
    vector <string> names;
    // Get the arguments
    if(input_type == "FASTA"){

        ifstream in(input_file.c_str());
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

    }else if(input_type == "CLUSTAL"){
        
        ifstream in(input_file.c_str());
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
    }
    else{
        cout << "Please give valid input type" << endl;
        exit (EXIT_FAILURE);
    }
    // number of sequences
    int n_seq = seqs.size();

    // Uses covariation/ Mutual Information to find probable structurally important base pairs
    string structure = MIVector(seqs,stacking);
    if(verbose){
      printf ("\t The number of sequences is %d\n", n_seq);
      printf ("\t The structure found through covariation of the alignment is: \n\n%s\n", structure.c_str());  
    }
    
    // The output file
    ofstream out(output_file, ofstream::trunc);
    for(int i=0; i<n_seq; ++i){
    
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
    
    return 0;
}