//
//  MC_finiteTemperatureEquilibriumQMC.cpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/28/16.
//
//

#include "MC_finiteTemperatureEquilibriumQMC.hpp"




// -------------------------------------------------------------------------------------
// Parent class of equilibrium quantum Monte Carlo simulations that we are going
// to write -- probably the Ising model and the Heisenberg model in zero magnetic field.
// We assume that these are going to be spin QMCs on a bipartite lattice.
// -------------------------------------------------------------------------------------

// constructor, setting measure...Q values to false
finiteTemperatureEquilibriumQMC::finiteTemperatureEquilibriumQMC(void) {
    InitializeRandomNumberGenerator();
    
    this->spins.resize(0);
    this->sites.resize(0);
    this->bonds.resize(0);
    this->NeelPattern.resize(0);
    
    this->clearMeasurements();
    
    this->measureSpinZCorrelatorsQ = false;
    this->measureMagnetizationQ = false;
    this->measureMagnetizationSquaredQ = false;
    this->measureMagnetizationFullCountingStatisticsQ = false;
    this->measureStaggeredMagnetizationQ = false;
    this->measureStaggeredMagnetizationSquaredQ = false;
    this->measureStaggeredMagnetizationFullCountingStatisticsQ = false;
    this->measureAvgEnergyQ = false;
    this->measureSpecificHeatQ = false;
    this->measureMagneticSusceptibilityQ = false;
};  // Tested and works

// destructor
finiteTemperatureEquilibriumQMC::~finiteTemperatureEquilibriumQMC(void) {
};  // Tested and works

// clear all physical quantities, but preserve the spin configuration and the lattice structure
void finiteTemperatureEquilibriumQMC::clearMeasurements(void) {
    int i;
    this->partitionFunction = 0.;
    this->magnetization = 0.;
    this->magnetizationSquared = 0.;
    this->magnetizationFullCountingStatistics.clear();
    this->staggeredMagnetization = 0.;
    this->staggeredMagnetizationSquared = 0.;
    this->staggeredMagnetizationFullCountingStatistics.clear();
    this->avgEnergy = 0.;
    this->avgEnergySquared = 0.;
    this->specificHeat = 0.;
    this->magneticSusceptibility = 0.;
    for(i=0; i<this->spinZCorrelator.size(); i++) {
        this->spinZCorrelator[i] = 0.;
    }
};  // Tested and works

// Check if the arrays we use are of appropriate size
bool finiteTemperatureEquilibriumQMC::check(void) const {
    if(this->spins.size() != this->sites.size()) {
        cerr << "Error in finiteTemperatureEquilibriumQMC::check(). 'spins' and 'sites' arrays are of different size." << endl;
        return false;
    }
    if(this->spins.size() != this->NeelPattern.size()) {
        cerr << "Error in finiteTemperatureEquilibriumQMC::check(). 'spins' and 'NeelPattern' arrays are of different size." << endl;
        return false;
    }
    return true;
};  // Looked the code through, looks OK

// update the bins of measured quantities with their expectation values in the current quantum state
void finiteTemperatureEquilibriumQMC::updateMeasuredQuantities(void) {
    this->updatePartitionFunction();
    
    if(this->measureMagnetizationQ || this->measureMagnetizationSquaredQ || this->measureMagnetizationFullCountingStatisticsQ) {
        this->updateMagnetization();
    };
    
    if(this->measureStaggeredMagnetizationQ || this->measureStaggeredMagnetizationSquaredQ || this->measureStaggeredMagnetizationFullCountingStatisticsQ) {
        this->updateStaggeredMagnetization();
    };
    
    if(measureSpinZCorrelatorsQ) {
        this->updateSpinZCorrelators();
    };
    
    if(this->measureAvgEnergyQ || measureSpecificHeatQ) {
        this->updateAvgEnergy();
        this->updateSpecificHeat();
    };
    
    if(this->measureMagneticSusceptibilityQ) {
        this->updateMagneticSusceptibility();
    };
};  // Tested and works

// turn on switches to measure physical quantities. Should be called before the simulation is started.
void finiteTemperatureEquilibriumQMC::measureMagnetization(bool Q) {
    if(this->measureMagnetizationQ != Q) {
        this->clearMeasurements();
        this->measureMagnetizationQ = Q;
    }
};  // Tested and works

void finiteTemperatureEquilibriumQMC::measureMagnetizationSquared(bool Q) {
    if(this->measureMagnetizationSquaredQ != Q) {
        this->clearMeasurements();
        this->measureMagnetizationSquaredQ = Q;
    }
};  // Tested and works

void finiteTemperatureEquilibriumQMC::measureMagnetizationFullCountingStatistics(bool Q, int nBins) {
    if((this->measureMagnetizationFullCountingStatisticsQ != Q) || (this->magnetizationFullCountingStatistics.getNumberOfBins() != nBins)) {
        this->clearMeasurements();
        this->measureMagnetizationFullCountingStatisticsQ = Q;
        this->magnetizationFullCountingStatistics.reset(-0.5 - 0.5/nBins, 0.5 + 0.5/nBins, nBins + 1);
    }
};  // Tested and works

void finiteTemperatureEquilibriumQMC::measureStaggeredMagnetization(bool Q) {
    if(this->measureStaggeredMagnetizationQ != Q) {
        this->clearMeasurements();
        this->measureStaggeredMagnetizationQ = Q;
    }
};  // Tested and works

void finiteTemperatureEquilibriumQMC::measureStaggeredMagnetizationSquared(bool Q) {
    if(this->measureStaggeredMagnetizationSquaredQ != Q) {
        this->clearMeasurements();
        this->measureStaggeredMagnetizationSquaredQ = Q;
    }
};  // Tested and works

void finiteTemperatureEquilibriumQMC::measureStaggeredMagnetizationFullCountingStatistics(bool Q, int nBins) {
    if(this->measureStaggeredMagnetizationFullCountingStatisticsQ != Q  || (this->staggeredMagnetizationFullCountingStatistics.getNumberOfBins() != nBins)) {
        this->clearMeasurements();
        this->measureStaggeredMagnetizationFullCountingStatisticsQ = Q;
        this->staggeredMagnetizationFullCountingStatistics.reset(-0.5 - 0.5/nBins, 0.5 + 0.5/nBins, nBins + 1);
    }
};  // Tested and works

void finiteTemperatureEquilibriumQMC::measureAvgEnergy(bool Q) {
    if(this->measureAvgEnergyQ != Q) {
        this->clearMeasurements();
        this->measureAvgEnergyQ = Q;
    }
};  // Tested and works

void finiteTemperatureEquilibriumQMC::measureSpecificHeat(bool Q) {
    if(this->measureSpecificHeatQ != Q) {
        this->clearMeasurements();
        this->measureSpecificHeatQ = Q;
    }
};  // Tested and works

void finiteTemperatureEquilibriumQMC::measureMagneticSusceptibility(bool Q) {
    if(this->measureMagneticSusceptibilityQ != Q) {
        this->clearMeasurements();
        this->measureMagneticSusceptibilityQ = Q;
    }
};  // Tested and works

void finiteTemperatureEquilibriumQMC::measureSpinZCorrelators(bool Q, vector< pair<int,int> > sitePairs) {
    bool changeQ = false;
    int i;
    
    if((this->measureSpinZCorrelatorsQ != Q) || (this->spinZCorrelatorSites.size() != sitePairs.size()) || (this->spinZCorrelator.size() != sitePairs.size())) {
        changeQ = true;
    } else {
        for(i=0; i<sitePairs.size(); i++) {
            if(this->spinZCorrelatorSites[i] != sitePairs[i]) {
                changeQ = true;
                break;
            }
        }
    }
    
    if(changeQ == true) {
        this->clearMeasurements();
        this->measureSpinZCorrelatorsQ = Q;
        this->spinZCorrelatorSites.resize(sitePairs.size());
        this->spinZCorrelator.resize(sitePairs.size());
        for(i=0; i<sitePairs.size(); i++) {
            this->spinZCorrelatorSites[i] = sitePairs[i];
            this->spinZCorrelator[i] = 0.0;
        };
    };
};  // Tested and works

void finiteTemperatureEquilibriumQMC::measureAllQuantities(void) {
    this->measureMagnetization(true);
    this->measureMagnetizationSquared(true);
    this->measureAvgEnergy(true);
    this->measureSpecificHeat(true);
    this->measureMagneticSusceptibility(true);
    this->measureStaggeredMagnetization(true);
    this->measureStaggeredMagnetizationSquared(true);
    this->measureMagnetizationFullCountingStatistics(true);
    this->measureStaggeredMagnetizationFullCountingStatistics(true);
    
    vector< pair<int,int> > sitePairs(this->sites.size());
    for(int i=0; i<this->sites.size(); i++) {
        sitePairs[i].first  = 0;
        sitePairs[i].second = i;
    }
    this->measureSpinZCorrelators(true, sitePairs);
};  // Tested and works

// get the current spin state of the system as a triple data: (x, y, spin), where x and y denote the coordinates of sites
vector< vector<int> > finiteTemperatureEquilibriumQMC::getCurrentSpinState(void) const {
    vector< vector<int> > currentSpinState(this->sites.size(), vector<int> (3));
    for(int i=0; i<this->sites.size(); i++) {
        currentSpinState[i][0] = this->sites[i][0];
        currentSpinState[i][1] = this->sites[i][1];
        currentSpinState[i][2] = this->spins[i];
    };
    return currentSpinState;
};  // Looked the code through, looks OK

// get the current value of the staggered magnetization
double finiteTemperatureEquilibriumQMC::getCurrentStaggeredMagnetization(void) const {
    if(spins.size() == 0) {
        return 0.;
    }
    
    double staggeredMagnetization;
    int i;
    for(staggeredMagnetization=0., i=0; i<this->spins.size(); i++) {
        staggeredMagnetization += (double) (spins[i] * NeelPattern[i]);
    };
    
    return (staggeredMagnetization / ((double) spins.size()));
};  // Looked the code through, looks OK

// returns an array of spin correlators and the site coordinates
// (x1, y1, x2, y2, <S_z(x1,y1) S_z(x2,y2)>)
vector< vector<double> > finiteTemperatureEquilibriumQMC::getSpinZCorrelators(void) const {
    vector< vector<double> > corr(spinZCorrelatorSites.size(), vector<double> (5));
    double x1, y1, x2, y2;
    double i1,i2;
    
    if(spinZCorrelatorSites.size() == 0) {
        return corr;
    }
    for(int i=0; i<spinZCorrelatorSites.size(); i++) {
        i1 = this->spinZCorrelatorSites[i].first;
        i2 = this->spinZCorrelatorSites[i].second;
        x1 = this->sites[i1][0];
        y1 = this->sites[i1][1];
        x2 = this->sites[i2][0];
        y2 = this->sites[i2][1];
        
        corr[i][0] = (double) x1;
        corr[i][1] = (double) y1;
        corr[i][2] = (double) x2;
        corr[i][3] = (double) y2;
        corr[i][4] = this->spinZCorrelator[i] / this->partitionFunction;
    };
    return corr;
};  // Looked the code through, looks OK

void finiteTemperatureEquilibriumQMC::updateMagnetization(void) {
    if(spins.size() == 0) {
        return;
    }
    
    double currentMagn;
    int i;
    
    for(currentMagn=0., i=0; i<this->spins.size(); i++) {
        currentMagn += (double) spins[i];
    };
    currentMagn /= (double) spins.size();
    
    // spins are stored as +1 / -1, whereas they should be +1/2 and -1/2
    // we correct for this when we calculate the magnetization
    currentMagn *= 0.5;
    
    this->magnetization        += currentMagn;
    this->magnetizationSquared += currentMagn * currentMagn;
    this->magnetizationFullCountingStatistics += currentMagn;
};

void finiteTemperatureEquilibriumQMC::updateStaggeredMagnetization(void) {
    if(spins.size() == 0) {
        return;
    }
    
    double currentStaggeredMagn;
    int i;
    
    for(currentStaggeredMagn=0., i=0; i<this->spins.size(); i++) {
        currentStaggeredMagn += (double) (spins[i] * NeelPattern[i]);
    };
    currentStaggeredMagn /= (double) spins.size();
    
    // spins are stored as +1 / -1, whereas they should be +1/2 and -1/2
    // we correct for this when we calculate the staggered magnetization
    currentStaggeredMagn *= 0.5;

    
    this->staggeredMagnetization        += currentStaggeredMagn;
    this->staggeredMagnetizationSquared += currentStaggeredMagn * currentStaggeredMagn;
    this->staggeredMagnetizationFullCountingStatistics += currentStaggeredMagn;
};  // Looked the code through, looks OK

void finiteTemperatureEquilibriumQMC::updateSpinZCorrelators(void) {
    int i1,i2,j;
    for(j=0; j<this->spinZCorrelatorSites.size(); j++) {
        i1 = this->spinZCorrelatorSites[j].first;
        i2 = this->spinZCorrelatorSites[j].second;
        this->spinZCorrelator[j] += 0.25 * spins[i1] * spins[i2];
        // spins are stored as +1 / -1, whereas they should be +1/2 and -1/2.
        // that's why we need the 0.25 factor
    };
};  // Looked the code through, looks OK

ostream & operator<<(ostream & out, const finiteTemperatureEquilibriumQMC & sim) {
    int i;
    
    out << endl;
    out << "Which quantities to measure: ";
    out << endl;
    out << "measureSpinZCorrelatorsQ: " << (sim.measureSpinZCorrelatorsQ ? "true" : "false") << endl;
    out << "measureMagnetizationQ: " << (sim.measureMagnetizationQ ? "true" : "false") << endl;
    out << "measureMagnetizationSquaredQ: " << (sim.measureMagnetizationSquaredQ ? "true" : "false") << endl;
    out << "measureMagnetizationFullCountingStatisticsQ: " << (sim.measureMagnetizationFullCountingStatisticsQ ? "true" : "false") << endl;
    out << "measureStaggeredMagnetizationQ: " << (sim.measureStaggeredMagnetizationQ ? "true" : "false") << endl;
    out << "measureStaggeredMagnetizationSquaredQ: " << (sim.measureStaggeredMagnetizationSquaredQ ? "true" : "false") << endl;
    out << "measureStaggeredMagnetizationFullCountingStatisticsQ: " << (sim.measureStaggeredMagnetizationFullCountingStatisticsQ ? "true" : "false") << endl;
    out << "measureAvgEnergyQ: " << (sim.measureAvgEnergyQ ? "true" : "false") << endl;
    out << "measureSpecificHeatQ: " << (sim.measureSpecificHeatQ ? "true" : "false") << endl;
    out << "measureMagneticSusceptibilityQ: " << (sim.measureMagneticSusceptibilityQ ? "true" : "false") << endl;
    out << endl;
    
    
    out << "Lattice parameters" << endl;
    out << "sites: " << endl;
    for(i=0; i<sim.sites.size(); i++) {
        out << "(" << sim.sites[i][0] << ", " << sim.sites[i][1] << "), ";
    };
    out << endl << endl;
    out << "bonds: " << endl;
    for(i=0; i<sim.bonds.size(); i++) {
        out << "(";
        out << sim.sites[sim.bonds[i].first][0] << ", " <<  sim.sites[sim.bonds[i].first][1] << "; ";
        out << sim.sites[sim.bonds[i].second][0] << ", " <<  sim.sites[sim.bonds[i].second][1] << "), ";
    };
    out << endl << endl;
    out << "spins: " << endl;
    for(i=0; i<sim.spins.size(); i++) {
        out << "(" << sim.sites[i][0] << ", " << sim.sites[i][1] << "; " << sim.spins[i] << "), ";
    };
    out << endl << endl;
    out << "Neel pattern: " << endl;
    for(i=0; i<sim.NeelPattern.size(); i++) {
        out << "(" << sim.sites[i][0] << ", " << sim.sites[i][1] << "; " << sim.NeelPattern[i] << "), ";
    };
    out << endl << endl;
    
    
    out << "Bins of measurable quantities: " << endl;
    out << "partitionFunction = " << sim.partitionFunction << endl;
    out << "magnetization = " << sim.magnetization/sim.partitionFunction << endl;
    out << "magnetizationSquared = " << sim.magnetizationSquared/sim.partitionFunction << endl;
    out << "magnetizationFullCountingStatistics = " << endl << sim.magnetizationFullCountingStatistics << endl << endl;
    out << "staggeredMagnetization = " << sim.staggeredMagnetization/sim.partitionFunction << endl;
    out << "staggeredMagnetizationSquared = " << sim.staggeredMagnetizationSquared/sim.partitionFunction << endl;
    out << "staggeredMagnetizationFullCountingStatistics = " << endl << sim.staggeredMagnetizationFullCountingStatistics << endl << endl;
    out << "avgEnergy = " << sim.avgEnergy/sim.partitionFunction << endl;
    out << "avgEnergySquared = " << sim.avgEnergySquared/sim.partitionFunction << endl;
    out << "specificHeat = " << sim.specificHeat/sim.partitionFunction << endl;
    out << "magneticSusceptibility = " << sim.magneticSusceptibility/sim.partitionFunction << endl;
    out << "spinZCorrelator = "<< endl;
    for(i=0; i<sim.spinZCorrelator.size(); i++) {
        out << "(" << sim.sites[sim.spinZCorrelatorSites[i].first][0] << ",";
        out << sim.sites[sim.spinZCorrelatorSites[i].first][1] << ", ";
        out << sim.sites[sim.spinZCorrelatorSites[i].second][0] << ",";
        out << sim.sites[sim.spinZCorrelatorSites[i].second][1] << ", ";
        out << sim.spinZCorrelator[i]/sim.partitionFunction << ") ";
        out << endl;
    }
    out << endl << endl;
    
    return out;
};  // Tested and works








// ------------------------------------------------------------------------------------------
// Child class of finiteTemperatureEquilibriumQMC implementing Monte Carlo class on a square
// lattice of size (2*LATTICE_SIZE)x(2*LATTICE_SIZE), where LATTICE_SIZE is a global constant
// defined in MC_coord.hpp.
// ------------------------------------------------------------------------------------------

squareLatticeFiniteTemperatureEquilibriumQMC::squareLatticeFiniteTemperatureEquilibriumQMC()
: finiteTemperatureEquilibriumQMC::finiteTemperatureEquilibriumQMC() {
    int i1;
    int x, y;
    int nSites = 4 * LATTICE_SIZE * LATTICE_SIZE;
    
    // set the indices of the sites and the bonds connecting them
    this->sites.resize(nSites);
    this->spins.resize(nSites);
    this->bonds.resize(2*nSites);
    this->NeelPattern.resize(nSites);
    
    for(i1=0; i1<nSites; i1++) {
        x = indexToX(i1);
        y = indexToY(i1);
        
        this->sites[i1].resize(2);
        this->sites[i1][0] = x;
        this->sites[i1][1] = y;
        
        this->bonds[i1].first = i1;
        this->bonds[i1].second = coordToIndex(x+1, y);
        this->bonds[i1 + nSites].first = i1;
        this->bonds[i1 + nSites].second = coordToIndex(x, y+1);
        
        this->NeelPattern[i1] = (mod(x+y,2) == 0 ? 1 : -1);
    }
    
    // randomize initial spin configurations
    vector<int> spinChoices = {-1,1};
    this->spins = RandomChoice(spinChoices, nSites);
} // Tested and works

squareLatticeFiniteTemperatureEquilibriumQMC::~squareLatticeFiniteTemperatureEquilibriumQMC(){
} // Tested and works
