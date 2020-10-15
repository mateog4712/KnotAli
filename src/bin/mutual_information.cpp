#include "../include/mutual_information.h"
#include "../include/CalculateProbability.h"
#include <vector>
#include <math.h>
#include <iostream>
#include <map>
#include <tuple>
#include <algorithm>

using namespace std;

string MIVector(vector<string> seqs, bool stack){
  

  // Initialize our vector of pairs with score > mean
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


  double column_max[n] = {0};
  double column_sum[n] = {0};
  double scores[n][n] = {0};
  double sum = 0;
  int count = 0;
  double maxim = 0;
  


  // Calculate the entropy values
  for(int i = 0; i<n;++i){
    for(int j = 4;(i+j)<n;++j){
        if(scores[i][i+j] == 0) scores[i][i+j] = calcMutualInformation(&cols[i][0],&cols[i+j][0], n_seq);
        if(i-1 > 0 && i+j+1 < n && scores[i-1][i+j+1] == 0) scores[i-1][i+j+1] = calcMutualInformation(&cols[i-1][0],&cols[i+j+1][0], n_seq); 
        if(scores[i+1][i+j-1] == 0 && (i+1 < n)) scores[i+1][i+j-1] = calcMutualInformation(&cols[i+1][0],&cols[i+j-1][0], n_seq); 
        
        double score;
        if(!stack) score = scores[i][i+j];
        else {
          double prev = (i-1 <=0 || i+j+1 >= n) ? scores[i][i+j] : scores[i-1][i+j+1];
          double foll = (i+1 >=n && j-i-1>4) ? scores[i][i+j] : scores[i+1][i+j-1];
          score = (prev + 2*scores[i][i+j]+ foll)/4;;
        }
        sum+=score;
        ++count;
        maxim =max(maxim,score); 
        column_max[i] =max(column_max[i],score);
        column_max[i+j] =max(column_max[i+j],score);  
        column_sum[i] += score;
        column_sum[i+j] += score;
    }
  }
  double mean = sum/count;
  for(int i = 0; i<n;++i){
    double column_mean_i = column_sum[i]/n;
    for(int j = 4;(i+j)<n;++j){
      double column_mean_ij = column_sum[i+j]/n;
      if(scores[i][i+j] == 0) continue;
      double score = scores[i][i+j]-(column_mean_i*column_mean_ij)/mean;
      if(score > .80){
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
    // If pair wont work due to previous pair
    if(structure[get<0>(hotspots[i].pair)] != '_' || structure[get<1>(hotspots[i].pair)] != '_') continue;

    // checks for pseudknots
    bool pseudoknot = check_Pseudoknot(used,hotspots[i]);
    if(!pseudoknot){
      // Replaces blank structure based on pairs
      structure[get<0>(hotspots[i].pair)] = '(';
      structure[get<1>(hotspots[i].pair)] = ')';
      // Pushes back used pairs
      used.push_back(hotspots[i].pair);
    }
  }
  for(int i = 0; i<n; ++i){
    if(structure[i] == '_' && column_max[i] < mean){
      structure[i] = '.';
    }
  }


  return structure;
}
double mi(JointProbabilityState state) {
  double mutualInformation = 0.0;
  double penalty = 1;
  int firstIndex,secondIndex;
  int i;
    
  /*
  ** I(X;Y) = \sum_x \sum_y p(x,y) * \log (p(x,y)/p(x)p(y))
  */
  for (i = 0; i < state.numJointStates; i++) {
    if(i==0) penalty-=state.jointProbabilityVector[i];
    if(i == 9 || i== 13 || i== 17 || i==19 || i==21 || i== 23){
    penalty-=state.jointProbabilityVector[i];
    firstIndex = i % state.numFirstStates;
    secondIndex = i / state.numFirstStates;
    
    if ((state.jointProbabilityVector[i] > 0) && (state.firstProbabilityVector[firstIndex] > 0) && (state.secondProbabilityVector[secondIndex] > 0)) {
      
      mutualInformation += state.jointProbabilityVector[i] * log(state.jointProbabilityVector[i] / state.firstProbabilityVector[firstIndex] / state.secondProbabilityVector[secondIndex]);
    }
    }
  }
  

  mutualInformation /= log(2);
  // Gap penalty added
  mutualInformation-=penalty;
 
  return mutualInformation;
}

double calcMutualInformation(uint *dataVector, uint *targetVector, int vectorLength) {
    JointProbabilityState state = calculateJointProbability(dataVector,targetVector,vectorLength);
    
    double mutualInformation = mi(state);
    freeJointProbabilityState(state);
    return mutualInformation;
}

double discAndCalcMutualInformation(double *dataVector, double *targetVector, int vectorLength) {
  JointProbabilityState state = discAndCalcJointProbability(dataVector,targetVector,vectorLength);
    
  double mutualInformation = mi(state);
  freeJointProbabilityState(state);
  return mutualInformation;
}

// Checks to see if pair is pseudoknotted
bool check_Pseudoknot(vector<tuple<int,int> > used, Hotspot hotspot){
  for(int j = 0; j<used.size();++j){
      if((get<0>(hotspot.pair) < get<0>(used[j])  && get<1>(hotspot.pair) >  get<0>(used[j]) && get<1>(hotspot.pair) <  get<1>(used[j])) || (get<0>(hotspot.pair) < get<1>(used[j])  && get<1>(hotspot.pair) >  get<1>(used[j]) && get<0>(hotspot.pair) >  get<0>(used[j]))) return true;
  }
  return false;
}




