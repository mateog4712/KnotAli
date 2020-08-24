#include "../include/CalculateProbability.h"
#include "../include/ArrayOperations.h"
#include "../include/covariation.h"
#include "../include/utils.h"
#include <iostream>
#include <stdlib.h>     /* calloc, exit, free */
#include <stdio.h>      /* printf, scanf, NULL */
#include <algorithm> // for std::find
#include <iterator> // for std::begin, std::end
#include <map>
#include <string>
#include <vector>
#include <iomanip>

using namespace std;
// structure holding the score and location of pair
struct Hotspot { 
    tuple<int,int> pair;
    double score;
};

string covVector(vector<string> seqs){
    vector<Hotspot> hotspots;


    // Make a map to easily index bases
    std::map<char,uint> base;
    base['-']=0;
    base['A']=1;
    base['C']=2;
    base['G']=3;
    base['U']=4;



    // For easier naming
    int n_seq = seqs.size();
    int n = seqs[0].length();

    uint cols[n][n_seq];
    // Find the columns for the list of seqs
    for(int i = 0; i<n;++i){
        for(int j = 0;j<n_seq;++j){

            cols[i][j] = base[seqs[j][i]]; 

        }
    }

    // Calculate the entropy values
    for(int i = 0; i<n;++i){
        for(int j = 4;(j<30 && (i+j)<n);++j){
            double score = covariation(&cols[i][0],&cols[i+j][0], n_seq); 
            if(score > .54) {
                Hotspot hotspot;
                hotspot.pair = make_tuple(i,i+j);
                hotspot.score = score;
                hotspots.push_back(hotspot);
            }
        }

    }
    // Sorts by decreasing score
    sort(hotspots.begin(), hotspots.end(), [](auto const &a, auto const &b) { return a.score > b.score; });

    // a string of unpaired bases
    string structure(n, '_');

    // vector of pairs that have been used
    vector<tuple<int,int> > used;

    // Go through each hotspot pair and change the structure accordingly
    for(int i = 0; i<hotspots.size();++i){
        bool pseudoknot = false;
        // If pair wont work due to previous pair
        if(structure[get<0>(hotspots[i].pair)] != '_' || structure[get<1>(hotspots[i].pair)] != '_') continue;

        // Checks to see if pair is pseudoknotted
        for(int j = 0; j<used.size();++j){
            if((get<0>(hotspots[i].pair) < get<0>(used[j])  && get<1>(hotspots[i].pair) >  get<0>(used[j])) || (get<0>(hotspots[i].pair) < get<1>(used[j])  && get<1>(hotspots[i].pair) >  get<1>(used[j]))) pseudoknot = true;
        }
        if(pseudoknot) continue;

        // Replaces blank structure based on pairs
        structure[get<0>(hotspots[i].pair)] = '(';
        structure[get<1>(hotspots[i].pair)] = ')';

        // Pushes back used pairs
        used.push_back(hotspots[i].pair);
    }

return structure;
}

double covariation(uint *dataVector, uint *targetVector, int vectorLength){
    int sum = 0;
    double penalty = 0;
    for(int i =0; i<vectorLength;++i){
        bool pairI = canPair( *(dataVector+i)*5+*(targetVector+i));
        penalty += !pairI - (*(dataVector+i)*5+*(targetVector+i) == 0)*.75;

        for(int j=i+1; j<vectorLength;++j){
            int hamming = 0;
            bool pairJ = canPair( *(dataVector+j)*5+*(targetVector+j));
            if(pairI && pairJ ) hamming = (*(dataVector+i) != *(dataVector+j)) + (*(targetVector+i) != *(targetVector+j));
            sum += hamming;
        }
    }
    return 1.0/nChoosek(vectorLength,2)*sum*.5 - 1.0/vectorLength*penalty;
    
}

unsigned nChoosek( unsigned n, unsigned k )
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

string covVectorStack(vector<string> seqs){
    vector<Hotspot> hotspots;


    // Make a map to easily index bases
    std::map<char,uint> base;
    base['-']=0;
    base['A']=1;
    base['C']=2;
    base['G']=3;
    base['U']=4;



    // For easier naming
    int n_seq = seqs.size();
    int n = seqs[0].length();

    uint cols[n][n_seq];
    // Find the columns for the list of seqs
    for(int i = 0; i<n;++i){
        for(int j = 0;j<n_seq;++j){

            cols[i][j] = base[seqs[j][i]]; 

        }
    }

    // Calculate the entropy values
    for(int i = 0; i<n;++i){
        for(int j = 4;(j<30 && (i+j)<n);++j){
            double score = covariation(&cols[i][0],&cols[i+j][0], n_seq);
            double prev = 0;
            double foll = 0;
            if(i-1 <=0 || i+j+1 >= n){
                prev = score;
            }
            else{
                prev = covariation(&cols[i-1][0],&cols[i+j+1][0], n_seq);
            }
            if(i+1 >=n && j-i-1>4){
                foll = score;
            }
            else{
                foll = covariation(&cols[i+1][0],&cols[i+j-1][0], n_seq);
            }
         
        
            score = (prev + 2*score+ foll)/4; 
            if(score > .28) {
                Hotspot hotspot;
                hotspot.pair = make_tuple(i,i+j);
                hotspot.score = score;
                hotspots.push_back(hotspot);
            }
        }

    }
    // Sorts by decreasing score
    sort(hotspots.begin(), hotspots.end(), [](auto const &a, auto const &b) { return a.score > b.score; });

    // a string of unpaired bases
    string structure(n, '_');

    // vector of pairs that have been used
    vector<tuple<int,int> > used;

    // Go through each hotspot pair and change the structure accordingly
    for(int i = 0; i<hotspots.size();++i){
        bool pseudoknot = false;
        // If pair wont work due to previous pair
        if(structure[get<0>(hotspots[i].pair)] != '_' || structure[get<1>(hotspots[i].pair)] != '_') continue;

        // Checks to see if pair is pseudoknotted
        for(int j = 0; j<used.size();++j){
            if((get<0>(hotspots[i].pair) < get<0>(used[j])  && get<1>(hotspots[i].pair) >  get<0>(used[j])) || (get<0>(hotspots[i].pair) < get<1>(used[j])  && get<1>(hotspots[i].pair) >  get<1>(used[j]))) pseudoknot = true;
        }
        if(pseudoknot) continue;

        // Replaces blank structure based on pairs
        structure[get<0>(hotspots[i].pair)] = '(';
        structure[get<1>(hotspots[i].pair)] = ')';

        // Pushes back used pairs
        used.push_back(hotspots[i].pair);
    }

return structure;
}