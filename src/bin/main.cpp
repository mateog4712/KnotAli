#include "../include/allFold.h"
#include "../include/mutual_information.h"
#include "../include/utils.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <unistd.h>

using namespace std;

#define MAXSLEN 5000
void usage (const char *prog)
// PRE:  None
// POST: Prints the usage message and exits
{
    printf ("\nUsage: %s -f <location of alignment> [options]\n\n", prog);
    printf ("  -f <location of alignment>\n");
    printf ("\t The alignment to fold, max length is %d\n\n", MAXSLEN);
    printf ("Options:\n");
    printf ("  -i <Alignment Type> \n");
    printf ("\t The type of alignment: FASTA or CLUSTAL (default is FASTA) \n\n");
    printf ("  -p <number of threads>\n");
    printf ("\t Runs the program in parallel using the specified number of threads\n\n");
    printf ("  -s\n\tTurn stacking on which takes into account surrounding base pairs while running, by default it is false\n\n");
    printf ("  -v\n\tGives a verbose output\n\n");
    printf ("  -h\n\tPrint this help message\n\n");
    printf ("Examples:\n");
    printf ("\t%s -f sample.afa -i FASTA \n", prog);
    printf ("\t%s -f sample.aln -i CLUSTAL \n", prog);
    printf ("\t%s -f sample.aln -i CLUSTAL -s\n", prog);
    exit (0);
}
int main(int argc, char *argv[]) {
    
    //int c;
    vector <string> list;
    string typeI = "FASTA";
    string inputFile;
    int threads = 1;
    bool stack = false;
    bool verbose = false;
    int c;

    opterr = 0;

    while ((c = getopt (argc, argv, "f:i:p:svh")) != -1)
        switch (c)
        {
            case 'f':
                inputFile = optarg;
                break;
            case 'i':
                typeI = optarg;
                break;
            case 'p':
                threads = atoi(optarg);
                break;
            case 's':
                stack = true;
                break;
            case 'v':
                verbose = true;
                break;    
            case 'h':
                usage(argv[0]);
                break;

            default:
                abort ();
        }

    // Define the data structures for the sequences
    vector <string> seqs;
    vector <string> seqs2;
    vector <string> names;
    // Get the arguments
    if(typeI == "FASTA"){

        ifstream in(inputFile.c_str());
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

    }else if(typeI == "CLUSTAL"){
        
        ifstream in(inputFile.c_str());
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
    string structure = MIVector(seqs,stack);
    if(verbose){
      printf ("\t The number of sequences is %d\n", n_seq);
      printf ("\t The structure found through covariation of the alignment is: \n\n%s\n", structure.c_str());  
    }
    
    // The output file
    ofstream out("results.afa", ofstream::trunc);
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