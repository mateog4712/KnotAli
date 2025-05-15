#include <string>
#include <vector>
#include <math.h>
#include <iostream>
#include <numeric>
#include "mutual_information.hh"
#include "utils.hh"
#include "CalculateProbability.hh"
#include <map>
#include <tuple>
#include <algorithm>
#include <cstring>
#include <fstream>

// Checks to see if pair is pseudoknotted
bool check_Pseudoknot(std::vector<std::tuple<int,int> > const& used, const Hotspot& hotspot){
  cand_pos_t n = used.size();
  for(cand_pos_t j = 0; j<n;++j){
      if((std::get<0>(hotspot.pair) < std::get<0>(used[j])  && std::get<1>(hotspot.pair) >  std::get<0>(used[j]) && std::get<1>(hotspot.pair) <  std::get<1>(used[j])) || (std::get<0>(hotspot.pair) < std::get<1>(used[j])  && std::get<1>(hotspot.pair) >  std::get<1>(used[j]) && std::get<0>(hotspot.pair) >  std::get<0>(used[j]))) return true;
  }
  return false;
  
}

// Calculates the APC value for finding MIp
double APC(double const& col_i, double const& col_j, double const& mean){

return (col_i*col_j)/mean;

}

std::string MIVector(std::vector<std::string> &seqs){
  

  // Initialize our vector of pairs with score > mean
  std::vector<Hotspot> hotspots;

  // Make a map to easily index bases
  std::map<char,uint> base;
  base['-']=0;
  base['A'] = base['W'] = base['M'] = base['H'] = base['N'] = 1;
  base['C'] = base['Y'] = base['B'] = 2;
  base['G'] = base['R'] = base['S'] = base['V'] = 3;
  base['U']=4;
  base['T'] = base['K'] = base['D'] = 5;


  // For easier naming
  int n_seq = seqs.size();
  int n = seqs[0].length();
  

  std::vector<uint> cols;
  std::vector<double> column_max;
  std::vector<double> column_sum;
  std::vector<double> cnt;



  cols.resize(n_seq*n,0);
  column_max.resize(n,0);
  column_sum.resize(n,0);
  cnt.resize(n,0);
  // Find the columns for the list of seqs
  for(int i = 0; i<n;++i){
    
    for(int j = 0;j<n_seq;++j){

      cols[i*n_seq+j] = base[seqs[j][i]];

    }
  }

  double sum = 0;
  int count = 0;

  // Calculate the entropy values
  for(int i = 0; i<n;++i){
    for(int j = 0;(i+j)<n && j<2000;++j){
        
      double score = calcMutualInformation(&cols[i*n_seq],&cols[(i+j)*n_seq], n_seq);
      sum+=score;
      ++count;
      column_max[i] =std::max(column_max[i],score);
      column_max[i+j] =std::max(column_max[i+j],score);  
      column_sum[i] += score;
      column_sum[i+j] += score;
      cnt[i]+=1;
      cnt[i+j]+=1;
    }
  }
  double mean = sum/count;


  for(int i = 0; i<n;++i){
    double column_mean_i = column_sum[i]/(cnt[i]); 
    for(int j = 4;(i+j)<n && j<2000;++j){
      double column_mean_ij = column_sum[i+j]/(cnt[i+j]);
      double colij = calcMutualInformation(&cols[i*n_seq],&cols[(i+j)*n_seq], n_seq);
      double score = colij - APC(column_mean_i,column_mean_ij,mean);

      if(score > .4){
        Hotspot hotspot;
        hotspot.pair = std::make_tuple(i,i+j);
        hotspot.score = score;
        hotspots.push_back(hotspot);
      }
    }
  }
  // Sorts by decreasing score
  std::sort(hotspots.begin(), hotspots.end(), [](auto const &a, auto const &b) { return a.score > b.score; });

  // a string of unpaired bases
  std::string structure(n, '.');

  // vector of pairs that have been used
  std::vector<std::tuple<int,int> > used;

  // Go through each hotspot pair and change the structure accordingly
  cand_pos_t n_hot = hotspots.size();
  for(cand_pos_t i = 0; i<n_hot;++i){
    // If pair wont work due to previous pair
    if(structure[std::get<0>(hotspots[i].pair)] != '.' || structure[std::get<1>(hotspots[i].pair)] != '.') continue;

    // checks for pseudoknots
    bool pseudoknot = check_Pseudoknot(used,hotspots[i]);
    if(!pseudoknot){
      // Replaces blank structure based on pairs
      structure[std::get<0>(hotspots[i].pair)] = '(';
      structure[std::get<1>(hotspots[i].pair)] = ')';
      // Pushes back used pairs
      used.push_back(hotspots[i].pair);
    }
  }
  int infoLoss=0;
  for(cand_pos_t i = 0; i<n; ++i){
    if(structure[i] == '.' && column_max[i] < mean){
      structure[i] = 'x';
      infoLoss++;
    }
  }

  if(infoLoss>2*n/3 ){
    std::cout << "Not enough info from sequence" << std::endl;
    std::cout << structure << std::endl;
    exit (EXIT_FAILURE);
  }
  return structure;
}

double mi(JointProbabilityState state) {
  double mutualInformation = 0.0;
  /*
  ** I(X;Y) = \sum_x \sum_y p(x,y) * \log (p(x,y)/p(x)p(y))
  */
  
  for (int i = 0; i < state.numJointStates; i++) {
    if(i == 10 || i==11 || i== 15 || i== 10 || i==22 || i==25 || i== 27 || i==31){
    
      int firstIndex = i % state.numFirstStates;
      int secondIndex = i / state.numFirstStates;

    
      if ((state.jointProbabilityVector[i] > 0) && (state.firstProbabilityVector[firstIndex] > 0) && (state.secondProbabilityVector[secondIndex] > 0)) {
        
        mutualInformation += (state.jointProbabilityVector[i]) * log((state.jointProbabilityVector[i]) / state.firstProbabilityVector[firstIndex] / state.secondProbabilityVector[secondIndex]);
      }
    }
  }
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