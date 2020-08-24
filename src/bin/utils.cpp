#include "../include/utils.h"
#include <algorithm>
#include <stdlib.h>
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
#include <cmath>

using namespace std;


bool canPair(int n){
    if(n == 9 || n == 13 || n == 17 || n == 19 || n == 21 || n == 23) return true;
    return false;
}



// Converts an integer to binary with number of bits = ns+1
int* toBinary(int n, int ns){
    ns +=1;
    int* a =  (int *) calloc(ns,sizeof(int));
    int i;
    for(i=0; n>0; i++)    
    {
        a[ns-1-i]= n%2;    
        n= n/2;  
    } 
return a;    
}
// Returns the hamming distance between i and j with i' and j'
int hamming(int i, int j, int x, int y){

int m = max(j,y);
int ns = log2(m);
int* ii = toBinary(i,ns);
int* jj = toBinary(j,ns);
int* xx = toBinary(x,ns);
int* yy = toBinary(y,ns);
int distance = 0;
for(int k =0;k<ns+1;++k){
    if(*(ii+k) != *(xx+k)) distance+=1;
    if(*(jj+k) != *(yy+k)) distance+=1;
}
    
free(ii);
free(jj);
free(xx);
free(yy);
    
return distance;
}



/** This function is used to get the lines within a txt file **/
bool getFileContent(string fileName, vector <string> &vecOfStrs) {
    // Open the File
    ifstream in(fileName.c_str());

    // Check if object is valid
    if (!in) {
        cerr << "Cannot open the File : " << fileName << endl;
        return false;
    }
    string str;
    // Read the next line from File untill it reaches the end.
    while (getline(in, str)) {
        // Line contains string of length > 0 then save it in vector
        if (str.size() > 0)
            vecOfStrs.push_back(str);
    }
    //Close The File
    in.close();
    return true;
}

bool fexists(const char *filename) {
    ifstream ifile(filename);
    return bool(ifile);
}

/** This functions takes a string and two chars. If char1 appears, it is replaced with char2**/
char* replaceChar(char *stri, char ch1, char ch2) {
    string seq(stri);
    // cout << seq.length() << endl;
    for (int i = 0; i < seq.length(); ++i) {
        if (seq[i] == ch1)
            seq[i] = ch2;
    }
    // char output[seq.length()+1];
    strcpy(stri, seq.c_str());
    return stri;
}



void ct2dotSin(string &family, vector <string> &lists) {

    // Will only add sequences that are the same length as the first one
    int seqLength;
    // Sets up the variables for converting to dot format
    vector <string> direct;
    string dir, filepath;
    DIR *dp;
    struct dirent *dirp;
    struct stat filestat;

    // User gives the folder to look in
    cout << "dir to get files of: " << flush;
    getline(cin, dir);
    

    // User gives the family to look at
    cout << "Family to look at: " << flush;
    cin >> family;

    // The while loop opens the directory and goes through each file
    dp = opendir(dir.c_str());
    while ((dirp = readdir(dp))) {
        filepath = dir + "/" + dirp->d_name;
        

        // If the file being looked at is a folder, it is skipped
        if (stat(filepath.c_str(), &filestat)) continue;
        if (S_ISDIR(filestat.st_mode)) continue;



        // if the file does not contain the family we are looking for, skip it
        int x = filepath.find(family);
        if (x == -1) continue;
        // If the file is not a .seq file, skip it
        int check = filepath.find(".seq");
        if (check == -1) continue;
        // Define the temporary variables
        string fileTemp;
        vector <string> lines;


        // Create the string for the .ct file
        fileTemp = filepath.substr(0, filepath.length() - 4) + ".ct";
        // Add the size of the .seq file's sequence to the variable seqLength for comparing for the first sequence
        if (direct.empty()) {
            getFileContent(filepath, lines);
            seqLength = lines[2].length();
        }
        // if a .ct file exists for the .seq file, check the length and add it if it is the same length
        if (fexists(fileTemp.c_str())) {
            getFileContent(filepath, lines);
            int length = lines[2].length();
            if (76 == length) {
                direct.push_back(fileTemp);
        
            }
        }


    }
    closedir(dp);
    // Deletes the . at the end of the name so as to not mess with the creation of text files
    if (family[family.length() - 1] == '.') family.erase(family.length() - 1, 1);


    // creates the directory where ct2dot will place its files
    string dir1 = "../output/rnastructure/" + family + "/";
    struct stat info;
    if (stat(dir1.c_str(), &info) != 0) {
        system(("mkdir " + dir1).c_str());
    }

    // Creates the directory where HFold will later put its files
    string dir2 = "../output/hfold/" + family + "/";
    struct stat info2;
    if (stat(dir2.c_str(), &info2) != 0) {
        system(("mkdir " + dir2).c_str());
    }
    // Creates the directory where the fasta file will go
    string dir3 = "../output/fasta/" + family + "/";
    if (stat(dir3.c_str(), &info2) != 0) {
        system(("mkdir " + dir3).c_str());
    }

    // Creates new file in case it already exists

    ostringstream oss;
    ofstream out(("../output/fasta/" + family + "/" + family + ".txt").c_str(), ofstream::trunc);
    out.close();
    ofstream out2(("../output/fasta/" + family + "/" + family + "Seq.txt").c_str(), ofstream::trunc);
    out2.close();
    for (int i = 0; i < 100;++i){//direct.size(); i++) {
        // Resets the oss stream
        oss.str("");
        // Gets the file name from direct
        string file = direct[i];
        // Takes the filenames from direct and erases the .ct and archiveII/ from the string
        string output = direct[i];
        output.erase(output.length() - 3, 3);
        int loc = file.find("/");
        output = output.substr(loc + 1, output.length());
        loc = output.find("/");
        output = output.substr(loc + 1, output.length());
        
        // Pushes back the locations of the files into a list for Hfold later
        lists.push_back(dir2 + output + ".txt");
        string outer = dir1 + output + ".txt";
        // Runs RNAstructure's CT2dot conversion
        
        
        oss << "./" << "../RNAstructure/exe/ct2dot " << file << " ALL " << outer;
        string command = oss.str();
        system(command.c_str());

        // Opens the current family file
        ifstream in(outer.c_str());
        //create an ofstream for appending
        ofstream out(("../output/fasta/" + family + "/" + family + ".txt").c_str(), ofstream::app);
        ofstream out1(("../output/fasta/" + family + "/" + family + "Seq.txt").c_str(), ofstream::app);
        string str;
        // Read the next line
        int j = 1;
        while (getline(in, str)) {
            // Reads line into fasta File
            out << str << endl;
            if(j%3 != 0) out1 << str << endl;
            ++j;
        }
        //Close The File
        in.close();
        out.close();
        out1.close();
    }
}
void ct2dot2(string &family, vector <string> &lists) {

    // Sets up the variables for converting to dot format
    vector <string> direct;
    string dir, filepath;
    DIR *dp;
    struct dirent *dirp;
    struct stat filestat;

    // User gives the folder to look in
    cout << "dir to get files of: " << flush;
    getline(cin, dir);
    

    // User gives the family to look at
    cout << "Family to look at: " << flush;
    cin >> family;

    // The while loop opens the directory and goes through each file
    dp = opendir(dir.c_str());
    while ((dirp = readdir(dp))) {
        filepath = dir + "/" + dirp->d_name;
        

        // If the file being looked at is a folder, it is skipped
        if (stat(filepath.c_str(), &filestat)) continue;
        if (S_ISDIR(filestat.st_mode)) continue;



        // if the file does not contain the family we are looking for, skip it
        int x = filepath.find(family);
        if (x == -1) continue;
        // // If the file is not a .seq file, skip it
        // int check = filepath.find(".seq");
        // if (check == -1) continue;
        // Define the temporary variables
        string fileTemp;
        vector <string> lines;


        // Create the string for the .ct file
        fileTemp = filepath.substr(0, filepath.length() - 4) + ".ct";
        // Add the size of the .seq file's sequence to the variable seqLength for comparing for the first sequence
        if (direct.empty()) {
            getFileContent(filepath, lines);
        }
        // if a .ct file exists for the .seq file, check the length and add it if it is the same length
        if (fexists(fileTemp.c_str())) {
            getFileContent(filepath, lines);
            
            direct.push_back(fileTemp);
        
        }


    }
    closedir(dp);
    // Deletes the . at the end of the name so as to not mess with the creation of text files
    if (family[family.length() - 1] == '.') family.erase(family.length() - 1, 1);


    // creates the directory where ct2dot will place its files
    string dir1 = "../output/rnastructure/" + family + "/";
    struct stat info;
    if (stat(dir1.c_str(), &info) != 0) {
        system(("mkdir " + dir1).c_str());
    }

    // Creates the directory where HFold will later put its files
    string dir2 = "../output/hfold/" + family + "/";
    struct stat info2;
    if (stat(dir2.c_str(), &info2) != 0) {
        system(("mkdir " + dir2).c_str());
    }
    // Creates the directory where the fasta file will go
    string dir3 = "../output/fasta/" + family + "/";
    if (stat(dir3.c_str(), &info2) != 0) {
        system(("mkdir " + dir3).c_str());
    }
    // Creates new file in case it already exists
    
    ostringstream oss;
    ofstream out(("../output/fasta/" + family + "/" + family + ".txt").c_str(), ofstream::trunc);
    out.close();
    ofstream out2(("../output/fasta/" + family + "/" + family + "Seq.txt").c_str(), ofstream::trunc);
    out2.close();
    for (int i = 0; i < direct.size(); i++) {
        // Resets the oss stream
        oss.str("");
        // Gets the file name from direct
        string file = direct[i];
        // Takes the filenames from direct and erases the .ct and archiveII/ from the string
        string output = direct[i];
        output.erase(output.length() - 3, 3);
        
        int loc = file.find("/");
        output = output.substr(loc + 1, output.length());
        loc = output.find("/");
        output = output.substr(loc + 1, output.length());
        
        // Pushes back the locations of the files into a list for Hfold later
        lists.push_back(dir2 + output + ".txt");
        string outer = dir1 + output + ".txt";
        // Runs RNAstructure's CT2dot conversion
        
        
        oss << "./" << "../RNAstructure/exe/ct2dot " << file << " ALL " << outer;
        string command = oss.str();
        system(command.c_str());
        

        // Opens the current family file
        ifstream in(outer.c_str());
        //create an ofstream for appending
        ofstream out(("../output/fasta/" + family + "/" + family + ".txt").c_str(), ofstream::app);
        ofstream out1(("../output/fasta/" + family + "/" + family + "Seq.txt").c_str(), ofstream::app);
        string str;
        // Read the next line
        int j = 1;
        
        while (getline(in, str)) {
            // Reads line into fasta File
            //if(j%3 == 1)cout << str << endl;
            out << str << endl;
            if(j%3 != 0) {
            out1 << str << endl;
            
            }
            ++j;
        }
        //Close The File
        in.close();
        out.close();
        out1.close();
    }
}



