#include "mutual_information.hh"
#include "CalculateProbability.hh"
#include <vector>
#include <math.h>
#include <iostream>
#include <map>
#include <tuple>
#include <algorithm>
#include <cstring>

// Checks to see if pair is pseudoknotted
auto const check_Pseudoknot(auto const& used, auto const& hotspot){
  for(int j = 0; j<used.size();++j){
      if((std::get<0>(hotspot.pair) < std::get<0>(used[j])  && std::get<1>(hotspot.pair) >  std::get<0>(used[j]) && std::get<1>(hotspot.pair) <  std::get<1>(used[j])) || (std::get<0>(hotspot.pair) < std::get<1>(used[j])  && std::get<1>(hotspot.pair) >  std::get<1>(used[j]) && std::get<0>(hotspot.pair) >  std::get<0>(used[j]))) return true;
  }

  // std::for_each( std::cbegin( used ),std::cend( used ),[&]( auto const tuple ) {
  //             auto const [first,second] =tuple;
  //             if((std::get<0>(hotspot.pair) < first  && std::get<1>(hotspot.pair) >  first && std::get<1>(hotspot.pair) <  second) || (std::get<0>(hotspot.pair) < second  && std::get<1>(hotspot.pair) >  second && std::get<0>(hotspot.pair) >  first)){
  //               return true;
  //             } 
  //           } );
  
  
  return false;
  
}




std::string MIVector(std::vector<std::string> seqs, bool stack){
  

  // Initialize our vector of pairs with score > mean
  std::vector<Hotspot> hotspots;

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


  uint cols[n][n_seq] = {0};
  // Find the columns for the list of seqs
  for(int i = 0; i<n;++i){
    for(int j = 0;j<n_seq;++j){

       cols[i][j] = base[seqs[j][i]]; 

    }
  }


  double column_max[n] = {0};
  double column_sum[n] = {0};
  double scores[n][n] = {0};
  memset(scores,0,n*n);
  double sum = 0;
  int count = 0;
  double maxim = 0;
  double APC[n][n] = {0};


  
  

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
          score = (prev + 2*scores[i][i+j]+ foll)/4;
          scores[i][i+j] = score;
        }
        sum+=score;
        ++count;
        maxim =std::max(maxim,score); 
        column_max[i] =std::max(column_max[i],score);
        column_max[i+j] =std::max(column_max[i+j],score);  
        column_sum[i] += score;
        column_sum[i+j] += score;
    }
  }
  std::vector <std::tuple<int,int> > stacks;
  double mean = sum/count;
  for(int i = 0; i<n;++i){
    double column_mean_i = column_sum[i]/n;
    for(int j = 4;(i+j)<n;++j){
      double column_mean_ij = column_sum[i+j]/n;
      APC[i][i+j] = scores[i][i+j]-(column_mean_i*column_mean_ij)/mean;
      
      if(APC[i][i+j] > .4){
        Hotspot hotspot;
        
        hotspot.pair = std::make_tuple(i,i+j);
        hotspot.score = APC[i][i+j];
        hotspots.push_back(hotspot);
      }
    }
  }
  // Sorts by decreasing score
  std::sort(hotspots.begin(), hotspots.end(), [](auto const &a, auto const &b) { return a.score > b.score; });

  // a string of unpaired bases
  std::string structure(n, '_');

  // vector of pairs that have been used
  std::vector<std::tuple<int,int> > used;

  // Go through each hotspot pair and change the structure accordingly
  for(int i = 0; i<hotspots.size();++i){
    // If pair wont work due to previous pair
    if(structure[std::get<0>(hotspots[i].pair)] != '_' || structure[std::get<1>(hotspots[i].pair)] != '_') continue;

    // checks for pseudknots
    auto const pseudoknot = check_Pseudoknot(used,hotspots[i]);
    if(!pseudoknot){
      // Replaces blank structure based on pairs
      structure[std::get<0>(hotspots[i].pair)] = '(';
      structure[std::get<1>(hotspots[i].pair)] = ')';
      // Pushes back used pairs
      used.push_back(hotspots[i].pair);
    }
  }

  int infoLoss=0;
  for(int i = 0; i<n; ++i){
    if(structure[i] == '_' && column_max[i] < mean){
      structure[i] = '.';
      infoLoss++;
    }
  }
  if(infoLoss>n/2){
    std::cout << "Not enough info from sequence" << std::endl;
    exit (EXIT_FAILURE);
  }


  return structure;
}
double mi(JointProbabilityState state) {
  double mutualInformation = 0.0;
  // double penalty = 0;
  int firstIndex,secondIndex;
  int i;
  /*
  ** I(X;Y) = \sum_x \sum_y p(x,y) * \log (p(x,y)/p(x)p(y))
  */
  for (i = 0; i < state.numJointStates; i++) {
    // if(i==0 || i==5 || i==10 || i==15 || i==20) penalty+=state.jointProbabilityVector[i];
    if(i == 9 || i== 13 || i== 17 || i==19 || i==21 || i== 23){
    
    firstIndex = i % state.numFirstStates;
    secondIndex = i / state.numFirstStates;

    
    if ((state.jointProbabilityVector[i] > 0) && (state.firstProbabilityVector[firstIndex] > 0) && (state.secondProbabilityVector[secondIndex] > 0)) {
      
      mutualInformation += state.jointProbabilityVector[i] * log(state.jointProbabilityVector[i] / state.firstProbabilityVector[firstIndex] / state.secondProbabilityVector[secondIndex]);
    }
    }
  }
  

  mutualInformation /= log(2);
  // Gap penalty added
  // mutualInformation-=(penalty/2);
 
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






