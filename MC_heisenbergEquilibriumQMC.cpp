//
//  MC_heisenbergEquilibriumQMC.cpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/28/16.
//
//

#include "MC_heisenbergEquilibriumQMC.hpp"

// constructor
heisenbergEquilibriumQMC::heisenbergEquilibriumQMC(double Jval, double Tval) : squareLatticeFiniteTemperatureEquilibriumQMC::squareLatticeFiniteTemperatureEquilibriumQMC() {
    this->J = abs(Jval);
    this->T = abs(Tval);
    this->initialize();
};  // Looked the code through, looks good

// destructor
heisenbergEquilibriumQMC::~heisenbergEquilibriumQMC() {
};  // Looked the code through, looks good


// reset data and reinitialize simulation
void heisenbergEquilibriumQMC::reset(double Jval, double Tval) {
    this->J = abs(Jval);
    this->T = abs(Tval);
    this->initialize();
};  // Looked the code through, looks good

// QMC routines
void heisenbergEquilibriumQMC::thermalize(int nBurnInSteps) {
    for(int i=0; i<nBurnInSteps; i++) {
        this->diagonalUpdate();
        this->operatorLoopUpdate();
        this->updatePowerCutoff();
    }
};  // Looked the code through, looks good

void heisenbergEquilibriumQMC::run(int nMonteCarloSteps) {
    for(int i=0; i<nMonteCarloSteps; i++) {
        this->diagonalUpdate();
        this->operatorLoopUpdate();
        this->updateMeasuredQuantities();
        this->updatePowerCutoff();
    }
};  // Looked the code through, looks good

// print out all data, using the operator << of the parent class as well
ostream & operator<<(ostream & out, const heisenbergEquilibriumQMC & sim) {
    out << endl << "**************************" << endl;
    out << "AFM Heisenberg model simulation" << endl;
    out << "J = " << sim.J << endl;
    out << "T = " << sim.T << endl;
    out << endl;
    out << static_cast<const finiteTemperatureEquilibriumQMC &>(sim);
    out << endl << "**************************" << endl << endl;
    return out;
};

// initialize private variables
void heisenbergEquilibriumQMC::initialize(void) {
    int i;
    
    // clear measurement bins
    this->clearMeasurements();
    
    // randomize spins
    vector<int> spinChoices = {1,-1};
    this->spins = RandomChoice(spinChoices, this->spins.size());
    
    // the simulation is started from a state with no Hamiltonian insertions
    this->hamiltonianPower = 0;
    this->hamiltonianPowerMax = 0;
    this->hamiltonianPowerCutoff = floor(this->sites.size() * 0.5 * abs(this->J / this->T));
    if(this->hamiltonianPowerCutoff < 100) {
        this->hamiltonianPowerCutoff = 100;
    }
    
    // we have no Hamiltonian insertions at the beginning, therefore oparray[] is zero
    this->oparray.resize(this->hamiltonianPowerCutoff);
    for(i=0; i<this->hamiltonianPowerCutoff; i++) {
        oparray[i] = 0;
    }
    
    // resize our auxiliary arrays
    this->linkedVertexList.resize(4 * hamiltonianPowerCutoff);
    this->vfirst.resize(sites.size());
    this->vlast.resize(sites.size());
    for(i=0; i<this->linkedVertexList.size(); i++) {
        this->linkedVertexList[i] = 0;
    }
    for(i=0; i<this->sites.size(); i++) {
        this->vfirst[i] = 0;
        this->vlast[i] = 0;
    }
    
    // determine the diagonal update probabilities
    this->pInsertDiagonal.resize(this->hamiltonianPowerCutoff);
    this->pDeleteDiagonal.resize(this->hamiltonianPowerCutoff);
    this->calculateDiagonalUpdateProbabilities();
};  // Looked the code through, looks good

// insertion or deletion of an S_z(i)S_z(j) diagonal operator, see Sandwik's notes
void heisenbergEquilibriumQMC::diagonalUpdate(void) {
    int p;      // index of the "time slice"
    int b;      // bond index
    int i,j;    // indices of the ends of the bonds
    
    for(p=0; p<this->hamiltonianPowerCutoff; p++) {
        // if there is no operator at time slice 'p', add a diagonal one with the
        // probability specified by pInsertDiagonal
        if(oparray[p] == 0) {
            b = RandomInteger(0, this->bonds.size() - 1);
            i = this->bonds[b].first;
            j = this->bonds[b].second;
            if(this->spins[i] == this->spins[j]) {
                // if the spins are the same at the ends of this bond, the matrix element
                // will be zero. Therefore, in this case, do nothing, go back to the for loop
            } else if(RandomReal(0., 1.) < pInsertDiagonal[hamiltonianPower]){
                oparray[p] = 2*b + 2;
                hamiltonianPower++;
                if(hamiltonianPowerMax < hamiltonianPower) {
                    hamiltonianPowerMax = hamiltonianPower;
                }
            }
        } else if (mod(oparray[p],2) == 0) { // if there is a diagonal operator, remove it
                                             // with probability pDeleteDiagonal
            if(RandomReal(0., 1.) < pDeleteDiagonal[hamiltonianPower]) {
                oparray[p] = 0;
                hamiltonianPower--;
            }
        } else { // if there is a spin-flip operator, we do nothing, just update
                 // the spins so that they encode the actual spin pattern at each
                 // imaginary time slice
            b = (oparray[p]/2) - 1;     // bond of this operator
            i = this->bonds[b].first;   // i and j are the two ends of the bond
            j = this->bonds[b].second;
            spins[i] = -spins[i];       // flip spins on these sites
            spins[j] = -spins[j];
        }
    }
};  // Looked the code through, looks good

// flip all spins and operators in each operator loop woth probability 1/2, see Sandwik's notes
void heisenbergEquilibriumQMC::operatorLoopUpdate(void) {
    // make sure that auxiliary arrays are of appropriate size
    this->vfirst.resize(sites.size());
    this->vlast.resize(sites.size());
    this->linkedVertexList.resize(4 * this->hamiltonianPowerCutoff);
    
    // define variables
    int i, s1, s2;  // site indices
    int b;          // bond index
    int p;          // imaginary time index
    int v0, v1, v2; // vertex indices
    
    // 1. we construct the linked vertex list
    // vertex legs are encoded as v = 4*p + l, where 'p' is the imaginary time index of the
    // operator, and 'l' is the leg index, l=0,1,2,3.
    // These legs are connected in the imaginary time direction at each site, see Sandwik's notes.
    // linkedVertexList contains the connectedness of these legs.
    
    // first, let us fill arrays vfirst and vlast with (-1)-s.
    // an element of these being -1 at the end of the construction of the linked vertex list
    // will indicate that the given site has no operator vertex attached to it
    for(i=0; i<this->vfirst.size(); i++) {
        this->vfirst[i] = -1;
        this->vlast[i] = -1;
    }
    
    // although Sandwik's notes do not talk about this step, I tink this is necessary.
    // we set all values of the linkedVertexList to (-1).
    // at the end of the construction of the linked vertex list, one element being -1
    // will indicate that there is no operator at the corresponding time slice v0/4.
    // later on, these will be the operator vertices, that will be left out of the loop update.
    for(v0=0; v0<4*this->hamiltonianPowerCutoff; v0++) {
        this->linkedVertexList[v0] = -1;
    }
    
    // go over the imaginary time direction, and construct the linked vertex list
    for(p=0; p<this->hamiltonianPowerCutoff; p++) {
        if(oparray[p]==0) { // if there is no operator at imaginary time slice 'p',
                            // then do nothing, go back to the for loop
        } else { // if there is an operator at imaginary time slice 'p' ...
            v0 = 4*p;               // "lower left leg" of the operator vertex
            b = (oparray[p]/2) - 1; // bond on which the operator acts
            s1 = bonds[b].first;    // s1 and s2 are endpoints of that bond
            s2 = bonds[b].second;
            v1 = vlast[s1];         // v1 and v2 are the vertex legs to which s1 and s2
            v2 = vlast[s2];         // are connected at an earlier time slice
            
            if(v1 != -1) {          // if there indeed is an earlier vertex on our records, we note
                                    // in linkedVertexList that these vertex legs are connected
                linkedVertexList[v1] = v0;
                linkedVertexList[v0] = v1;
            } else {                // if there was no operator legs at site s1 at earlier imaginary
                                    // times, we record that we have now visited vfirst
                vfirst[s1] = v0;
            }
            
            if(v2 != -1) {          // do the same thing for site s2
                                    // (this step is written incorrectly in Sandwik's notes)
                linkedVertexList[v2] = v0 + 1;
                linkedVertexList[v0 + 1] = v2;
            } else {
                vfirst[s2] = v0 + 1;
            }
            
            vlast[s1] = v0 + 2;
            vlast[s2] = v0 + 3;
        }
    }
    
    // update those connections that go through the imaginary time boundary
    for(i=0; i<sites.size(); i++) {
        if(vfirst[i] != -1) {
            linkedVertexList[vfirst[i]] = vlast [i];
            linkedVertexList[vlast[i]] = vfirst[i];
        }
    }
    // at this point, the connectedness of all vertex legs is stored in 'linkedVertexList'.
    // those sites 'i' that do not have an operator leg attached to them will have
    // vfirst[i] == vlast[i] == -1
    
    // 2. find all "operator loops" (see Sandwik's notes), and flip all of them with probability 1/2
    for(v0=0; v0<4*this->hamiltonianPowerCutoff; v0+=2) {
        if(linkedVertexList[v0]<0) {
            // Do nothing, go back to the for loop
        } else {
            v1 = v0;
            if(RandomReal(0., 1.) < 0.5) {      // with probability 1/2, traverse the loop, but do not
                                                // flip the spins or the operators (i.e. make
                                                // linkedVertexList[v1] = -1 for all vertex legs v1)
                do {
                    linkedVertexList[v1] = -1;              // the value -1 indicates that this vertex will not be flipped
                    p = v1/4;                               // time slice of the vertex
                    v2 = (mod(v1,2) == 0 ? v1+1 : v1-1);    // go to the neighboring leg of the vertex -> v2
                    v1 = linkedVertexList[v2];              // now go to the vertex that v2 is connected to
                    linkedVertexList[v2] = -1;
                } while (v1 != v0);
                
            } else {                            // with probability 1/2, traverse the loop, and flip the spins
                                                // and the operators (i.e. set linkedVertexList[v1] = -2)
                do {
                    linkedVertexList[v1] = -2;              // the value -2 indicates that this vertex will be flipped
                    p = v1/4;                               // time slice of the vertex
                    oparray[p] = (mod(oparray[p],2)==0 ? oparray[p]+1 : oparray[p]-1); // flip operator
                    v2 = (mod(v1,2) == 0 ? v1+1 : v1-1);    // go to the neighboring leg of the vertex -> v2
                    v1 = linkedVertexList[v2];              // now go to the vertex that v2 is connected to
                    linkedVertexList[v2] = -2;
                } while (v1 != v0);
            }
        }
    }
    
    // 3. flip spins in the selected operator loops
    // we can use the information stored in vfirst, vlast and linkedVertexList to figure out which spins to flip
    for(i=0; i<sites.size(); i++) {
        v0 = vfirst[i];
        if(v0 == -1) {  // if this spin has no operator leg attached to it, then flip it with probability 1/2
            if(RandomReal(0., 1.) < 0.5) {
                spins[i] = -spins[i];
            }
        } else {        // if there is an operator attached to it, then flip it only if we decided to do so
                        // at the construction of the operator loops
            if(linkedVertexList[v0]==-2) {
                spins[i] = -spins[i];
            }
        }
    }
};  // Looked the code through, looks good

// auxiliary function to update the values of pInsertDiagonal and pDeleteDiagonal
// needed when the maximal cutoff of powers of the Hamiltonian is changed
void heisenbergEquilibriumQMC::calculateDiagonalUpdateProbabilities(void) {
    int n;
    int Nb = this->bonds.size();
    int L = this->hamiltonianPowerCutoff;
    
    this->pInsertDiagonal.resize(this->hamiltonianPowerCutoff);
    this->pDeleteDiagonal.resize(this->hamiltonianPowerCutoff);
    
    for(n=0; n<this->hamiltonianPowerCutoff; n++) {
        pInsertDiagonal[n] = J * ((double) Nb) / (T * 2. * (L-n));
        pDeleteDiagonal[n] = 2. * T * ((double) L - n + 1) / (J * ((double) Nb));
        if(pInsertDiagonal[n] > 1.) {
            pInsertDiagonal[n] = 1.;
        }
        if(pDeleteDiagonal[n] > 1.) {
            pDeleteDiagonal[n] = 1.;
        }
    }
};  // Looked the code through, looks good

// auxiliary routine to update L during the thermalization process, if nmax becomes too large
// According to Sandwik, the ideal cutoff is L = nmax + nmax/3
void heisenbergEquilibriumQMC::updatePowerCutoff(void) {
    // if hamiltonianPowerMax becomes too large, we reset it
    int formerHamiltonianPowerCutoff;
    if(this->hamiltonianPowerMax > (int) (0.8 * this->hamiltonianPowerCutoff)) {
        formerHamiltonianPowerCutoff = this->hamiltonianPowerCutoff;
        this->hamiltonianPowerCutoff = this->hamiltonianPowerMax + this->hamiltonianPowerMax/3;
        if(this->hamiltonianPowerCutoff <= formerHamiltonianPowerCutoff) {
            this->hamiltonianPowerCutoff = formerHamiltonianPowerCutoff;
            return;
        }
    } else {
        return;
    }
    
    
    int i1;
    
    // when hamiltonianPowerCutoff is reset, we need to recalculate the diagonal update probabilities.
    // as the added imaginary time slices contain no operators, we need to set oparray in these slices to zero
    oparray.resize(this->hamiltonianPowerCutoff);
    for(i1=formerHamiltonianPowerCutoff; i1<this->hamiltonianPowerCutoff; i1++) {
        oparray[i1] = 0;
    }
    
    // we also need to resize some of our arrays
    this->linkedVertexList.resize(4 * this->hamiltonianPowerCutoff);
    this->pInsertDiagonal.resize(this->hamiltonianPowerCutoff);
    this->pDeleteDiagonal.resize(this->hamiltonianPowerCutoff);
    
    // after the probability arrays were made longer, we need to calculate them again
    this->calculateDiagonalUpdateProbabilities();
};  // Looked the code through, looks good

// update average energy (see Sandwik's notes about why avgEnergy = <-hamiltonianPower> * T).
// The minus sign is there because the Hamiltonian was multiplied by -1, so that the QMC measure be positive.
void heisenbergEquilibriumQMC::updateAvgEnergy(void) {
    this->avgEnergy += (-this->hamiltonianPower) * this->T / ((double) this->sites.size());
    this->avgEnergySquared += this->hamiltonianPower * (this->hamiltonianPower-1) * this->T * this->T / ((double) this->sites.size() * this->sites.size());
};  // Looked the code through, looks good

// update magnetic susceptibility X = <(S^z_tot)^2> / T / N
void heisenbergEquilibriumQMC::updateMagneticSusceptibility(void) {
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
void heisenbergEquilibriumQMC::updateSpecificHeat(void) {
    this->specificHeat = this->T * (this->avgEnergySquared - (this->avgEnergy * this->avgEnergy / this->partitionFunction));
}

