//
//  MC_random.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 9/20/16.
//
//

#ifndef MC_random_hpp
#define MC_random_hpp

#include <stdio.h>
#include <random>
#include <cmath>
#include <iostream>
#include <vector>
#include <sys/time.h>
#include "MC_clusterRunQ.hpp"

using namespace std;

// ---------------------------
// INITIALIZE WITH RANDOM SEED
// ---------------------------
void InitializeRandomNumberGenerator(void);

// ---------------
// RANDOM INTEGERS
// ---------------
int RandomInteger(int iMin, int iMax);                     // creates a single uniformly distr. integer in the range [min, max]
vector<int> RandomInteger(int iMin, int iMax, int n);      // creates n independent uniformly distr. integers in the range [min, max]

// ---------------
// RANDOM REALS
// ---------------
// creates a single uniformly distr. number in the range [min, max]
inline double RandomReal(double dMin, double dMax) {
    double r = ((double)rand() / (double) RAND_MAX);
    return dMin + r * (dMax - dMin);
} // Tested and works

// creates n independent uniformly distr numbers in the range [min, max]
inline vector<double> RandomReal(double dMin, double dMax, int n) {
    vector<double> numbers(n);
    double d = dMax-dMin;
    for(int i=0; i<n; i++) {
        numbers[i] = dMin + d * (((double)rand()) / ((double) RAND_MAX));
    }
    return numbers;
} // Tested and works


// ---------------------------
// RANDOM POISSON DISTRIBUTION
// ---------------------------
// ------------------------------------------------------------------------------
// THE FOLLOWING TWO FUNCTIONS (two version of RandomPoisson) NEED C++11 FEATURES
// SINCE THESE ARE NOT SUPPORTED BY THE CLUSTER, I WILL COMMENT THESE FUNCTIONS OUT.
// HOWEVER, THEY WORK PERFECTLY, AND CAN BE USED LATER ON OTHER PLATFORMS
// TO SUBSTITUTE THEM, I WROTE RandomPoissonDummy, THAT WORKS IDENTICALLY TO THEM
// ------------------------------------------------------------------------------
#ifndef RUN_ON_CLUSTER
int RandomPoisson(double mean);                     // creates a single integer with Poisson distribution P_i = (mean^i/i!) Exp[-mean]
vector<int> RandomPoisson(double mean, int n);      // creates n independent integers with Poisson distribution P_i = (mean^i/i!) Exp[-mean]
#endif
int RandomPoissonDummy(double mean);                // creates a random Poisson number, without using C++11 features


// ----------------
// CREATE HISTOGRAM
// ----------------
vector< vector<double> > Histogram(vector<double> const& randomNumbers, int nBins);
vector< vector<double> > Histogram(vector<int> const& randomNumbers, int nBins);
// creates a histogram from the elements of randomNumbers, with the number of bins being nBins
// Histogram returns triplets {binMin, binMax, freq}, where (binMin, binMax) denotes the range of the bin, and freq is its frequency among the elements in randomNumbers


// ---------------------------------------
// RANDOMLY SHUFFLE AND SELECT FROM ARRAYS
// ---------------------------------------
template <typename arrayType>
void RandomShuffle(arrayType * a, int n); // Randomly shuffles elements of the array of length n

// Randomly shuffles elements of the array of length n
template <typename myType>
void RandomShuffle(myType * a, int n){
    int i,j;
    myType tmp;
    for(i=0; i<(n-1); i++) {
        j = RandomInteger(i, n-1);
        tmp  = a[i];
        a[i] = a[j];
        a[j] = tmp;
    }
}
// Randomly shuffles elements of the vector of length n
template <typename myType>
void RandomShuffle(vector<myType> & a){
    int i,j;
    int n = a.size();
    myType tmp;
    for(i=0; i<(n-1); i++) {
        j = RandomInteger(i, n-1);
        tmp  = a[i];
        a[i] = a[j];
        a[j] = tmp;
    }
}

// Randomly selects 'n' elements from 'choices', with repetitions
template <typename T>
vector<T> RandomChoice(vector<T> choices, int n) {
    int len = choices.size();
    int i,j;
    vector<T> sample(n);
    for (i=0; i<n; i++) {
        j = RandomInteger(0,len-1);
        sample[i] = choices[j];
    }
    return sample;
}

#endif /* MC_random_hpp */


// ALL OF THESE ROUTINES HAVE BEEN TESTED AND THEY WORK APPROPRIATELY
