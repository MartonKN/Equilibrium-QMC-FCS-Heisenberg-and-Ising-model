//
//  MC_histogram.cpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/28/16.
//
//

#include "MC_histogram.hpp"


// constructor
histogram::histogram(void) {
    this->minValue=0.;
    this->maxValue=0.;
    this->binWidth=0.;
    this->nBins=0;
    this->bins.resize(0);
    this->nLowerThanMinValue=0;
    this->nGreaterThanMaxValue=0;
}; // Tested and works

// constructor
histogram::histogram(double minimumValue, double maximumValue, int numberOfBins) {
    this->reset(minimumValue, maximumValue, numberOfBins);
}; // Tested and works

// destructor
histogram::~histogram() {
    
}; // Tested and works

// clear all bins
void histogram::clear(void) {
    for (int i=0; i<this->bins.size(); i++) {
        this->bins[i] = 0;
    }
    this->nLowerThanMinValue=0;
    this->nGreaterThanMaxValue=0;
}; // Tested and works

// reset the number of bins, the limits and clear everything
void histogram::reset(double minimumValue, double maximumValue, int numberOfBins) {
    if(minimumValue < maximumValue) {
        this->minValue=minimumValue;
        this->maxValue=maximumValue;
    } else {
        this->minValue=maximumValue;
        this->maxValue=minimumValue;
    }
    
    if(numberOfBins < 1) {
        cerr << "Error in histogram::reset. Cannot create " << numberOfBins << " bins. ";
        cerr << "Bin number has been changed to 1." << endl;
        this->nBins = 1;
    } else {
        this->nBins=numberOfBins;
    }
    
    this->binWidth = (this->maxValue - this->minValue) / ((double) this->nBins);
    
    this->bins.resize(nBins);
    for(int i=0; i<nBins; i++) {
        this->bins[i]=0;
    }
    this->nLowerThanMinValue=0;
    this->nGreaterThanMaxValue=0;
}; // Tested and works

// add an element to the histogram
histogram & histogram::operator+=(const int & val) {
    (*this) += (double) val;
    return *this;
}; // Tested and works

histogram & histogram::operator+=(const double & val) {
    if(val < this->minValue) {
        this->nLowerThanMinValue++;
        return *this;
    }
    
    if(val >= this->maxValue) {
        this->nGreaterThanMaxValue++;
        return *this;
    }
    
    int i = floor((val - this->minValue) / this->binWidth);
    this->bins[i]++;
    return *this;
}; // Tested and works

histogram & histogram::operator+=(const vector<int> & vals) {
    for(int i=0; i<vals.size(); i++) {
        (*this) += (double) vals[i];
    }
    return *this;
}; // Tested and works

histogram & histogram::operator+=(const vector<double> & vals) {
    for(int i=0; i<vals.size(); i++) {
        (*this) += vals[i];
    }
    return *this;
}; // Tested and works

histogram & histogram::operator+=(const histogram & hist) {
    if(this->minValue != hist.minValue || this->maxValue != hist.maxValue || this->nBins != hist.nBins) {
        cerr << "Error in histogram::operator+=(const histogram &). The two histograms have different parameters." << endl;
        cerr << "We did not update the values of *this." << endl;
        return *this;
    }
    
    this->nLowerThanMinValue += hist.nLowerThanMinValue;
    this->nGreaterThanMaxValue += hist.nGreaterThanMaxValue;
    for(int i=0; i<bins.size(); i++) {
        this->bins[i] += hist.bins[i];
    }
    
    return *this;
}; // Tested and works

// assignment operator
histogram & histogram::operator=(const histogram & hist) {
    this->minValue = hist.minValue;
    this->maxValue = hist.maxValue;
    this->binWidth = hist.binWidth;
    this->nBins = hist.nBins;
    this->bins.resize(hist.bins.size());
    for(int i=0; i<this->bins.size(); i++) {
        this->bins[i] = hist.bins[i];
    }
    this->nLowerThanMinValue = hist.nLowerThanMinValue;
    this->nGreaterThanMaxValue = hist.nGreaterThanMaxValue;
    
    return *this;
}; // Tested and works


// printing
void histogram::printCounts(void) const {
    cout << *this;
}; // Tested and works

void histogram::printProbabilities(void) const {
    int i;
    long int totalCount;
    double binmin;
    double binmax;
    
    totalCount = nLowerThanMinValue + nGreaterThanMaxValue;
    for(i=0; i<this->bins.size(); i++) {
        totalCount += bins[i];
    }
    
    cout << "bin_min  bin_max  probs" << endl;
    binmin = minValue;
    binmax = minValue + binWidth;
    if(totalCount != 0) {
        for(i=0; i<this->bins.size(); i++) {
            cout << setprecision(6) << fixed << binmin << "\t" << binmax << "\t";
            cout << (double) this->bins[i] / ((double) totalCount) << endl;
            binmin = binmax;
            binmax += binWidth;
        }
    } else {
        for(i=0; i<this->bins.size(); i++) {
            cout << setprecision(6) << fixed << binmin << "\t" << binmax << "\t";
            cout << 0. << endl;
            binmin = binmax;
            binmax += binWidth;
        }
    }
}; // Tested and works

void histogram::printPDF(void) const {
    int i;
    long int totalCount;
    double binmin;
    double binmax;
    
    totalCount = nLowerThanMinValue + nGreaterThanMaxValue;
    for(i=0; i<this->bins.size(); i++) {
        totalCount += bins[i];
    }
    
    cout << "bin_min  bin_max  probs" << endl;
    binmin = minValue;
    binmax = minValue + binWidth;
    if(totalCount != 0) {
        for(i=0; i<this->bins.size(); i++) {
            cout << setprecision(6) << fixed << binmin << "\t" << binmax << "\t";
            cout << ((double) this->bins[i]) / this->binWidth / ((double) totalCount) << endl;
            binmin = binmax;
            binmax += binWidth;
        };
    } else {
        for(i=0; i<this->bins.size(); i++) {
            cout << setprecision(6) << fixed << binmin << "\t" << binmax << "\t";
            cout << 0. << endl;
            binmin = binmax;
            binmax += binWidth;
        };
    }
}; // Tested and works

std::ostream& operator<<(std::ostream & out, const histogram & hist) {
    out << "bin_min  bin_max  counts" << endl;
    double binmin = hist.minValue;
    double binmax = hist.minValue + hist.binWidth;
    
    for(int i=0; i<hist.bins.size(); i++) {
        out << setprecision(6) << fixed << binmin << "\t" << binmax << "\t";
        out << hist.bins[i] << endl;
        binmin = binmax;
        binmax += hist.binWidth;
    }
    
    return out;
}; // Tested and works


