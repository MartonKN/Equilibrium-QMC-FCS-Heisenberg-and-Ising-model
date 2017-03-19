//
//  MC_histogram.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/28/16.
//
//


#ifndef MC_histogram_hpp
#define MC_histogram_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// Simple histogram implementation assuming uniform bin width
// All routines have been tested and they work appropriately
class histogram {
protected:
    double minValue;                // lowest value in the histogram
    double maxValue;                // largest value in the histogram
    double binWidth;                // width of each bins (we choose uniform widths)
    int nBins;                      // numebr of bins
    vector<long int> bins;          // bins with element counts
    long int nLowerThanMinValue;    // number of counts lower than this->minValue
    long int nGreaterThanMaxValue;  // number of counts graeter than this->maxValue
    
public:
    // constructors and destructor
    histogram(void);
    histogram(double minimumValue, double maximumValue, int numberOfBins=100);
    ~histogram(void);
    
    // clear all bins
    void clear(void);
    
    // reset the number of bins, the limits and clear everything
    void reset(double minimumValue, double maximumValue, int numberOfBins=100);
    
    // add an element to the histogram
    histogram & operator+=(const int & val);
    histogram & operator+=(const double & val);
    histogram & operator+=(const vector<int> & vals);
    histogram & operator+=(const vector<double> & vals);
    histogram & operator+=(const histogram & hist);
    
    // assignment operator
    histogram & operator=(const histogram & hist);
    
    // get basic data
    inline int getNumberOfBins(void) const {return this->nBins;};
    inline double getMinValue(void) const {return this->minValue;};
    inline double getMaxValue(void) const {return this->maxValue;};
    
    // priting
    void printCounts(void) const;
    void printProbabilities(void) const;
    void printPDF(void) const;
    friend std::ostream& operator<<(std::ostream & out, const histogram & hist);
};


#endif /* MC_histogram_hpp */
