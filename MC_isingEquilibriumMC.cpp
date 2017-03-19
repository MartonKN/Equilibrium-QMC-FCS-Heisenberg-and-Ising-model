//
//  MC_isingEquilibriumMC.cpp
//  
//
//  Created by Marton Kanasz-Nagy on 11/2/16.
//
//

#include "MC_isingEquilibriumMC.hpp"


// constructor and destructor
isingEquilibriumMC::isingEquilibriumMC(double Jval, double Tval)
: squareLatticeFiniteTemperatureEquilibriumQMC::squareLatticeFiniteTemperatureEquilibriumQMC() {
    this->J = abs(Jval);
    this->T = abs(Tval);
    this->initialize();
};

isingEquilibriumMC::~isingEquilibriumMC() {
};

// reset data and reinitialize simulation
void isingEquilibriumMC::reset(double Jval, double Tval) {
    this->J = abs(Jval);
    this->T = abs(Tval);
    this->initialize();
};

// QMC routines
void isingEquilibriumMC::thermalize(int nBurnInSteps) {
    for(int i=0; i<nBurnInSteps; i++) {
        this->wangUpdate();
        // this->singleSpinFlipUpdate();
    }
};
void isingEquilibriumMC::run(int nMonteCarloSteps) {
    for(int i=0; i<nMonteCarloSteps; i++) {
        this->wangUpdate();
        // this->singleSpinFlipUpdate();
        this->updateMeasuredQuantities();
    }
};

// print out all data, using the operator << of the parent class as well
ostream & operator<<(ostream & out, const isingEquilibriumMC & sim) {
    out << endl << "**************************" << endl;
    out << "Ising model simulation" << endl;
    out << "J = " << sim.J << endl;
    out << "T = " << sim.T << endl;
    out << endl;
    out << static_cast<const finiteTemperatureEquilibriumQMC&>(sim);
    out << endl << "**************************" << endl << endl;
    return out;
};


// initialize the private variables
void isingEquilibriumMC::initialize(void) {
    int i, j;
    
    // clear measurement bins
    this->clearMeasurements();
    
    // set the probability of adding two AFM aligned spins to the cluster
    this->pAdd = 1. - exp(-0.5 * this->J / this->T);
    
    // determine which sites are nearest neighbors, using the arrays 'bonds' and 'sites' in the parent class
    this->nearestNeighbors.resize(this->sites.size());
    this->nearestNeighbors.shrink_to_fit();
    for(i=0; i<this->sites.size(); i++) {
        this->nearestNeighbors[i].resize(0);
        for(j=0; j<this->bonds.size(); j++) {
            if(this->bonds[j].first == i) {
                this->nearestNeighbors[i].push_back(this->bonds[j].second);
            } else if (this->bonds[j].second == i) {
                this->nearestNeighbors[i].push_back(this->bonds[j].first);
            }
        }
        this->nearestNeighbors[i].shrink_to_fit();
    }
    
    // allocate memory for 'isMember'
    this->isMember.reserve(sites.size());
    this->isMember.resize(sites.size());
    
    // randomize spins in the initial state, and calculate the energy corresponding to this state
    vector<int> spinChoices = {1,-1};
    this->spins = RandomChoice(spinChoices, this->sites.size());
    this->energy = this->getCurrentEnergy();
};

// check consistency of inner variables
bool isingEquilibriumMC::check(void) const {
    int i,j;
    int i1,i2;
    
    // Call the 'check()' routine of the parent class to see if everything is OK there
    if(!finiteTemperatureEquilibriumQMC::check()) {
        return false;
    }
    
    // Check if the coupling and temperature are positive
    if(this->J < 0) {
        cerr << "Error in isingEquilibriumMC::check(). The AFM coupling is negative." << endl;
        return false;
    }
    if(this->T < 0) {
        cerr << "Error in isingEquilibriumMC::check(). The temperature is negative." << endl;
        return false;
    }
    
    // Check if the probability of adding a bond to the cluster is appropriate
    if(this->pAdd != (1. - exp(-0.5 * this->J / this->T))) {
        cerr << "Error in isingEquilibriumMC::check(). pAdd != 1 - exp(-J/2T)." << endl;
        return false;
    }
    
    // Check if the energy of the current state is appropriate
    if(this->energy != this->getCurrentEnergy()) {
        cerr << "Error in isingEquilibriumMC::check(). Energy of the current spin state is not correct." << endl;
        return false;
    }
    
    // Check if the private arrays of this class have an appropriate sign
    if(this->isMember.size() != this->sites.size()) {
        cerr << "Error in isingEquilibriumMC::check(). The size of the array isMember is not right." << endl;
        return false;
    }
    if(this->nearestNeighbors.size() != this->sites.size()) {
        cerr << "Error in isingEquilibriumMC::check(). The size of the array nearestNeighbors is not right." << endl;
        return false;
    }    
    
    // Check if the array 'nearestNeighbors' indeed contains the nearest neighbors of the lattice
    bool foundIt;
    bool foundPair;
    for(j=0; j<this->bonds.size(); j++) {
        i1 = this->bonds[j].first;
        i2 = this->bonds[j].second;
        for(foundIt=false, i=0; i<this->nearestNeighbors[i1].size(); i++) {
            if(this->nearestNeighbors[i1][i] == i2) {
                foundIt = true;
                break;
            }
        }
        if(!foundIt) {
            cerr << "Error in isingEquilibriumMC::check(). nearestNeighbors was set up inappropriately." << endl;
            return false;
        }
        for(foundIt=false, i=0; i<this->nearestNeighbors[i2].size(); i++) {
            if(nearestNeighbors[i2][i] == i1) {
                foundIt = true;
                break;
            }
        }
        if(!foundIt) {
            cerr << "Error in isingEquilibriumMC::check(). nearestNeighbors was set up inappropriately." << endl;
            return false;
        }
    }
    for(i1=0; i1<this->sites.size(); i1++) {
        foundPair = false;
        for(j=0; j<this->bonds.size(); j++) {
            if(this->bonds[j].first == i1) {
                i2 = this->bonds[j].second;
                foundPair = true;
            } else if(this->bonds[j].second == i1) {
                i2 = this->bonds[j].first;
                foundPair = true;
            }
        }
        if(foundPair) {
            for(foundIt=false, i=0; i<this->nearestNeighbors[i1].size(); i++) {
                if(this->nearestNeighbors[i1][i] == i2) {
                    foundIt = true;
                    break;
                }
            }
            if(!foundIt) {
                cerr << "Error in isingEquilibriumMC::check(). nearestNeighbors was set up inappropriately." << endl;
                return false;
            }
            for(foundIt=false, i=0; i<nearestNeighbors[i2].size(); i++) {
                if(nearestNeighbors[i2][i] == i1) {
                    foundIt = true;
                    break;
                }
            }
            if(!foundIt) {
                cerr << "Error in isingEquilibriumMC::check(). nearestNeighbors was set up inappropriately." << endl;
                return false;
            }
        }
    }
    
    return true;
};

// one Monte Carlo step of the Wang cluster algorithm
void isingEquilibriumMC::wangUpdate(void) {
    int i1;
    int seedSite;
    
    // initialize the array storing whether a given site is a member of the flipped cluster
    for(i1=0; i1<this->isMember.size(); i1++) {
        this->isMember[i1] = false;
    }
    
    // Choose an initial site where the cluster is grown from
    seedSite = RandomInteger(0, this->sites.size()-1);
    
    // Grow the cluster iteratively, and flip the spins in the meantime
    this->growCluster(seedSite);
    
    // After the update has been made, let us check if the inner variables are consistent
    // This is only for testing the simulation, and will be commented out later on
    this->check();
};

// add AFM ordered neighboring sites of site 'site' to the cluster, with probability pAdd (Wang algorithm)
void isingEquilibriumMC::growCluster(const int siteInd) {
    // Flip the spin of 'siteInd', and indicate that is in the cluster
    this->isMember[siteInd] = true;
    this->spins[siteInd] *= -1;
    
    // Go over all the neighbors of 'siteInd', and try to add them to the cluster if they are not part of it already
    // Also update the change of the energy of the state due to the flip of the spin at 'siteInd'
    int neighborSiteInd;
    for(int i=0; i<this->nearestNeighbors[siteInd].size(); i++) {
        neighborSiteInd = this->nearestNeighbors[siteInd][i];
        this->energy += 0.5 * this->J * this->spins[siteInd] * this->spins[neighborSiteInd];
    };
    
    // We add 'neighborSiteInd' to the cluster with probability 'this->pAdd' if and only if
    // its spin in the initial spin configuration was the opposite as that of 'this->spins[siteInd]'.
    // However, since we already flipped the spin 'this->spins[siteInd]', we need 'spins[siteInd]==spins[neighborSiteInd]'
    // in the current spin state. Note that 'spins[neighborSiteInd]' was not flipped, since it is not part of the cluster.
    for(int i=0; i<this->nearestNeighbors[siteInd].size(); i++) {
        neighborSiteInd = this->nearestNeighbors[siteInd][i];
        if((!this->isMember[neighborSiteInd]) && (this->spins[siteInd]==this->spins[neighborSiteInd])) {
            if(RandomReal(0.,1.) < this->pAdd) {
                this->growCluster(neighborSiteInd);
            }
        }
    }
};

// one Monte Carlo step using the Metropolis update of a single spin
// much slower dynamics than that generated by 'wangUpdate', especially close to the ordering temperature.
// in the current implementation, it is not even optimized.
// it will be useful to test the results of the algorithm using the routine 'wangUpdate'.
void isingEquilibriumMC::singleSpinFlipUpdate(void) {
    int siteInd = RandomInteger(0, this->sites.size()-1);
    int neighborSiteInd;
    
    // Determine the energy change the spin flip would cause
    double dE = 0.;
    for(int i=0; i<this->nearestNeighbors[siteInd].size(); i++) {
        neighborSiteInd = this->nearestNeighbors[siteInd][i];
        dE -= (double) (this->spins[siteInd] * this->spins[neighborSiteInd]);
    }
    dE *= this->J * 0.5;
    
    if(dE <= 0.) {
        this->spins[siteInd] *= -1;
        this->energy += dE;
    } else if(RandomReal(0., 1.) < exp(-dE/this->T)) {
        this->spins[siteInd] *= -1;
        this->energy += dE;
    }
    
};

// calculate the energy of the current spin state
// energy = sum_<ij> J*spins[i]*spins[j]/4, where the 1/4 factor is there since the array 'spins' stores the spin
// state as (+1)/(-1)-s, not as (+1/2)/(-1/2)-s
double isingEquilibriumMC::getCurrentEnergy(void) const {
    double currentEnergy;
    int i1,i2, j;
    for(currentEnergy=0., i1=0; i1<sites.size(); i1++) {
        for(j=0; j<nearestNeighbors[i1].size(); j++) {
            i2 = this->nearestNeighbors[i1][j];
            currentEnergy += (double) (this->spins[i1] * this->spins[i2]);
        }
    }
    currentEnergy *= J;
    currentEnergy *= 0.25; // because we store the spins as (+1)/(-1)
    currentEnergy *= 0.5;  // because we counted each bond twice
    
    return currentEnergy;
};

// update average energy and magnetic susceptibility
void isingEquilibriumMC::updateAvgEnergy(void) {
    this->avgEnergy += this->energy / ((double) this->spins.size());
    this->avgEnergySquared += (this->energy * this->energy) / ((double) this->spins.size() * this->spins.size());
};

void isingEquilibriumMC::updateMagneticSusceptibility(void) {
    double totalCurrentMagnetization;
    int i;
    
    for(totalCurrentMagnetization=0., i=0; i<this->spins.size(); i++) {
        totalCurrentMagnetization += spins[i];
    }
    totalCurrentMagnetization *= 0.5;
    // spins are stored as +1 / -1, whereas they should be +1/2 and -1/2
    // we correct for this when we calculate the magnetization
    
    this->magneticSusceptibility += (totalCurrentMagnetization * totalCurrentMagnetization) / this->T / this->spins.size();
};

// update the specific heat C = T * (<E^2> - <E>^2) / N^2
void isingEquilibriumMC::updateSpecificHeat(void) {
    this->specificHeat = this->T * (this->avgEnergySquared - (this->avgEnergy * this->avgEnergy / this->partitionFunction));
}

