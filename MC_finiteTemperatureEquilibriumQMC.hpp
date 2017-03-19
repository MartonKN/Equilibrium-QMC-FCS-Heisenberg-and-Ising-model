//
//  MC_finiteTemperatureEquilibriumQMC.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/28/16.
//
//

#ifndef MC_finiteTemperatureEquilibriumQMC_hpp
#define MC_finiteTemperatureEquilibriumQMC_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <utility>
#include "MC_coord.hpp"
#include "MC_histogram.hpp"
#include "MC_random.hpp"

using namespace std;


// --------------------------------------------------------------------------------------
// Parent class of equilibrium quantum Monte Carlo simulations that we are going to write
// -- probably the Ising model and the Heisenberg model in zero magnetic field.
// We assume that these are going to be spin QMCs on a bipartite lattice.
// --------------------------------------------------------------------------------------

class finiteTemperatureEquilibriumQMC {
public:
    // Basic parameters of the lattice
    vector<int> spins;              // spin states of the lattice, at each lattice site, +1 / -1
    vector< vector<int> > sites;    // array storing the coordinates of the sites
    vector< pair<int,int> > bonds;  // array containing the endpoints of each bond 'b', accessed as 'bonds[b].first' and 'bonds[b].second'
    vector<int> NeelPattern;        // (+1) on even sites, (-1) on odd sites. This array is needed to measure the staggered magnetization
    
    // which quantities' thermal average to measure
    bool measureSpinZCorrelatorsQ;
    bool measureMagnetizationQ;
    bool measureMagnetizationSquaredQ;
    bool measureMagnetizationFullCountingStatisticsQ;
    bool measureStaggeredMagnetizationQ;
    bool measureStaggeredMagnetizationSquaredQ;
    bool measureStaggeredMagnetizationFullCountingStatisticsQ;
    bool measureAvgEnergyQ;
    bool measureSpecificHeatQ;
    bool measureMagneticSusceptibilityQ;
    
    // variables to store measured quantities
    long double partitionFunction;
    long double magnetization;
    long double magnetizationSquared;
    long double avgEnergy;
    long double avgEnergySquared;
    long double specificHeat;
    long double magneticSusceptibility;
    long double staggeredMagnetization;
    long double staggeredMagnetizationSquared;
    histogram magnetizationFullCountingStatistics;
    histogram staggeredMagnetizationFullCountingStatistics;
    vector< pair<int, int> > spinZCorrelatorSites;  // which pairs of sites to measure the spin correlators at
    vector<double> spinZCorrelator;                 // values of the spin correlator at these sites
    
    
public:
    // Constructor and destructor
    finiteTemperatureEquilibriumQMC(void);          // constructor, setting measure...Q values to false
    ~finiteTemperatureEquilibriumQMC(void);         // destructor
    
    
    // ----------------------------------------------
    // ROUTINES TO BE OVERWRITTEN BY THE CHILD CLASS
    // In addtion, 'sites', 'bonds' and 'NeelPattern'
    // need to be specified.
    // ----------------------------------------------
public:
    // QMC routines (these are the ones that will be overwritten at each child of this class)
    virtual void thermalize(int nBurnInSteps) {};
    virtual void run(int nMonteCarloSteps) {};
    
protected:
    // update average energy of the system and its magnetic susceptibility
    virtual void updateAvgEnergy(void) {};
    virtual void updateSpecificHeat(void) {};
    virtual void updateMagneticSusceptibility(void) {};
    
    
    // -----------------------------------------
    // ROUTINES DEFINIED WITHIN THE PARENT CLASS
    // -----------------------------------------
public:
    // clear all physical quantities, but preserve the spin configuration and the lattice structure
    void clearMeasurements(void);
    
    // check if the arrays we use are of appropriate size
    bool check(void) const;
    
    // update the bins of measured quantities with their expectation values in the current quantum state
    void updateMeasuredQuantities(void);
    
    // turn on switches to measure physical quantities. Should be called before the simulation is started.
    // when a switch is flipped, all previous measurement data is deleted
    void measureMagnetization(bool Q = true);
    void measureMagnetizationSquared(bool Q = true);
    void measureAvgEnergy(bool Q = true);
    void measureSpecificHeat(bool Q=true);
    void measureMagneticSusceptibility(bool Q = true);
    void measureStaggeredMagnetization(bool Q = true);
    void measureStaggeredMagnetizationSquared(bool Q = true);
    void measureSpinZCorrelators(bool Q, vector< pair<int,int> > sitePairs);
    void measureMagnetizationFullCountingStatistics(bool Q = true, int nBins=100);
    void measureStaggeredMagnetizationFullCountingStatistics(bool Q = true, int nBins=100);
    virtual void measureAllQuantities(void);
    
    // get the current spin state of the system as a triple data: (x, y, spin), where x and y denote the coordinates of sites
    vector< vector<int> > getCurrentSpinState(void) const;
    
    // get the current value of the staggered magnetization
    double getCurrentStaggeredMagnetization(void) const;
    
    // return values of the measured quantities
    inline double getMagnetization(void) const {return (this->magnetization/this->partitionFunction);};
    inline double getMagnetizationSquared(void) const {return (this->magnetizationSquared/this->partitionFunction);};
    inline double getAvgEnergy(void) const {return (this->avgEnergy/this->partitionFunction);};
    inline double getSpecificHeat(void) const {return (this->specificHeat/this->partitionFunction);};
    inline double getMagneticSusceptibility(void) const {return (this->magneticSusceptibility/this->partitionFunction);};
    inline double getStaggeredMagnetization(void) const {return (this->staggeredMagnetization/this->partitionFunction);};
    inline double getStaggeredMagnetizationSquared(void) const {return (this->staggeredMagnetizationSquared/this->partitionFunction);};
    inline histogram getMagnetizationFullCountingStatistics(void) const {return this->magnetizationFullCountingStatistics;};
    inline histogram getStaggeredMagnetizationFullCountingStatistics(void) const{return this->staggeredMagnetizationFullCountingStatistics;};
    vector< vector<double> > getSpinZCorrelators(void) const;
    
    // print out all data
    friend ostream& operator<<(ostream &, const finiteTemperatureEquilibriumQMC &);
    
protected:
    // update values of the physical quantities we measure
    inline void updatePartitionFunction(void) {
        partitionFunction++;
    };
    void updateMagnetization(void);
    void updateStaggeredMagnetization(void);
    void updateSpinZCorrelators(void);
};






// ------------------------------------------------------------------------------------------
// Child class of finiteTemperatureEquilibriumQMC implementing Monte Carlo class on a square
// lattice of size (2*LATTICE_SIZE)x(2*LATTICE_SIZE), where LATTICE_SIZE is a global constant
// defined in MC_coord.hpp.
// ------------------------------------------------------------------------------------------
class squareLatticeFiniteTemperatureEquilibriumQMC : public finiteTemperatureEquilibriumQMC {
public:
    squareLatticeFiniteTemperatureEquilibriumQMC();
    ~squareLatticeFiniteTemperatureEquilibriumQMC();
    void measureAllQuantities(void) {
        this->finiteTemperatureEquilibriumQMC::measureAllQuantities();
        this->measureStaggeredMagnetizationFullCountingStatistics(true, this->spins.size());
        this->measureMagnetizationFullCountingStatistics(true, this->spins.size());
    };
};




#endif /* MC_finiteTemperatureEquilibriumQMC_hpp */
