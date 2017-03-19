//
//  MC_random.cpp
//  
//
//  Created by Marton Kanasz-Nagy on 9/20/16.
//
//

#include "MC_random.hpp"


// ---------------------------
// INITIALIZE WITH RANDOM SEED
// ---------------------------
void InitializeRandomNumberGenerator(void){
    struct timeval tp;
    gettimeofday(&tp, NULL);
    long int time_in_microsec = tp.tv_sec * 1000000 + tp.tv_usec;
    srand (abs(time_in_microsec) % RAND_MAX);
}


// ---------------
// RANDOM INTEGERS
// ---------------
// creates a single uniformly distr. integer in the range [min, max]
int RandomInteger(int iMin, int iMax) {
    int d = abs(iMax-iMin)+1; // range of random numbers
    int b = min(iMin,iMax); // bottom of the range
    
    // If the range is zero
    if(iMin==iMax) {
        return iMin;
    }
    
    return ((rand() % d) + b);
} // Tested and works

// creates n independent uniformly distr. integer in the range [min, max]
vector<int> RandomInteger(int iMin, int iMax, int n)  {
    vector<int> numbers(n);
    int i;
    int d = abs(iMax-iMin)+1; // range of random numbers
    int b = min(iMin,iMax); // bottom of the range
    
    // If the range is zero
    if(iMin==iMax) {
        for(i=0; i<n; i++) {
            numbers[i] = iMin;
        }
        return numbers;
    }
    
    for(i=0; i<n; i++) {
        numbers[i] = (rand() % d) + b;
    }
    return numbers;
} // Tested and works





// ---------------------------
// RANDOM POISSON DISTRIBUTION
// ---------------------------
// Creates a random Poisson number, without using C++11 features
// Shuold work identically to RandomPoisson
int RandomPoissonDummy(double mean) {
    double x = RandomReal(0., 1.);      // We draw a random number, and check in which bin of the
                                        // cumulative distribution function it falls
                                        // The maximal value of n that can be drawn is maxn.
    
    // We only allow n-s below 10 standard deviations (10*sqrt(mean)) from the mean
    int maxn =(int) max(mean + 10.0 * sqrt(mean), 15.);
    
    int n=0;                    // n will be drawn with Poisson distribution
    double pn = exp(-mean);     // pn = (mean^n/n!) Exp[-mean] probability of drawing n
    double CDFn = pn;           // cumulative distribution function \sum_{i=0}^n (mean^i/i!) Exp[-i]
    while(x>CDFn && n<=maxn) {
        n++;
        pn *= mean/n;
        CDFn += pn;
    }
    return n;
} // Tested and works

#ifndef RUN_ON_CLUSTER
// ------------------------------------------------------------------------------
// THE FOLLOWING TWO FUNCTIONS (two version of RandomPoisson) NEED C++11 FEATURES
// SINCE THESE ARE NOT SUPPORTED BY THE CLUSTER, I WILL COMMENT THESE FUNCTIONS OUT
// HOWEVER, THEY WORK PERFECTLY, AND CAN BE USED LATER ON OTHER PLATFORMS
// TO SUBSTITUTE THEM, I WROTE RandomPoissonDummy, THAT WORKS IDENTICALLY TO THEM
// ------------------------------------------------------------------------------
// creates a single integer with Poisson distribution P_i = (mean^i/i!) Exp[-mean]
int RandomPoisson(double mean) {
    // initialize random number generator with a random seed
    default_random_engine generator(rand());
    
    // create class that generates random poisson distributed numbers
    poisson_distribution<int> distribution(mean);
    return distribution(generator);
} // Tested and works

// creates n independent integers with Poisson distribution P_i = (mean^i/i!) Exp[-mean]
vector<int> RandomPoisson(double mean, int n) {
    vector<int> numbers(n);
    
    // initialize random number generator with a random seed
    default_random_engine generator(rand());
    
    // create class that generates random poisson distributed numbers
    poisson_distribution<int> distribution(mean);
    for(int i=0; i<n; i++) {
        numbers[i] = distribution(generator);
    }
    return numbers;
} // Tested and works

#endif





// ----------------
// CREATE HISTOGRAM
// ----------------
// creates a histogram from the elements of randomNumbers, with the number of bins being nBins
// returns triplets {binMin, binMax, freq}, where (binMin, binMax) denotes the range of the bin, and freq is its frequency among the elements in randomNumbers
vector< vector<double> > Histogram(vector<double> const& randomNumbers, int nBins) {
    double minValue, maxValue;      // minimal and maximal value in randomNumbers
    int len = randomNumbers.size(); // number of elements in the randomNumbers
    int i;                          // iterator
    int iBin;                       // bin index
    
    // determine minimum and maximum in randomNumbers
    minValue = (double) randomNumbers[0];
    maxValue = (double) randomNumbers[0];
    for(i=1; i<len; i++) {
        if(minValue > ((double) (randomNumbers[i])) ) {
            minValue = (double) randomNumbers[i];
        }
        if(maxValue < ((double) (randomNumbers[i])) ) {
            maxValue = (double) randomNumbers[i];
        }
    }
    
    // if there is only a single element in randomNumbers, we will return a single element bin
    if(len==1 || minValue==maxValue) {
        vector< vector<double> > hist(1,vector<double>(3));
        hist[0][0]=(double) randomNumbers[0];
        hist[0][1]=(double) randomNumbers[0];
        hist[0][2]=1.;
        return hist;
    }
    
    // determine the ranges of bins
    vector< vector<double> > hist(nBins, vector<double>(3));
    double d = (maxValue-minValue) / ((double) nBins);  // range of each bin
    double invd = 1.0/d;                                  // inverse of the range
    double invlen = 1.0/((double) len);                 // 1.0/(number of elements)
    hist[0][0] = minValue;
    hist[0][1] = minValue + d;
    for(i=1; i<nBins; i++) {
        hist[i][0] = hist[i-1][1];
        hist[i][1] = hist[i][0] + d;
        hist[i][2] = 0.;
    }
    
    // count the frequency of elements in randomNumbers
    for(i=0; i<len; i++) {
        iBin = floor( (((double) randomNumbers[i]) - minValue)*invd );
        if(iBin==nBins) {
            iBin--;
        }
        hist[iBin][2] += invlen;
    }
    
    return hist;
} // Tested and works

vector< vector<double> > Histogram(vector<int> const& randomNumbers, int nBins) {
    vector<double> randomNumbersDouble(randomNumbers.begin(), randomNumbers.end());
    return Histogram(randomNumbersDouble, nBins);
} // Tested and works
