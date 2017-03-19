//
//  MC_heisenbergEquilibriumQMC.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/28/16.
//
//

#ifndef MC_heisenbergEquilibriumQMC_hpp
#define MC_heisenbergEquilibriumQMC_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include "MC_random.hpp"
#include "MC_finiteTemperatureEquilibriumQMC.hpp"

using namespace std;

// --------------------------------------------------------------------------------------------
// Child class of squareLatticeFiniteTemperatureEquilibriumQMC implementing an equilibrium
// quantum Monte Carlo calculation of the SU(2) symmetric AFM quantum Heisenberg model, without
// magnetic field. Lattice size is set in MC_coord.hpp by the global variable LATTICE_SIZE.
// --------------------------------------------------------------------------------------------
// This algorithm was written using the tutorials on Anders Sandwik's website and his papers.
// http://physics.bu.edu/~sandvik/vietri/sse.pdf
// Since most details are on the website, I will only write short comments to the code.
// ------------------------------------------------------------------------------------------
class heisenbergEquilibriumQMC : public squareLatticeFiniteTemperatureEquilibriumQMC {
protected:
    double T;                       // temperature
    double J;                       // spin coupling (positive since we assume antiferromagnetic ordering)
    
private:
    int hamiltonianPower;           // power of the Hamiltonian in its Taylor expansion
    int hamiltonianPowerMax;        // maximal Hamiltonian power used so far (to update hamiltonianPowerCutoff)
    int hamiltonianPowerCutoff;     // maximal allowed power of the Hamiltonian
    vector<int> oparray;            // array containing which operator we have at time p = 0...L-1
                                    // oparray[p] == 0: no operator
                                    // oparray[p] >  0: oparray[p] = a(p) + 2*bondindex(p) + 1
                                    // denoted s(p) in Sandwik's notes
                                    // Sandwik uses s(p) = a(p) + 2*bondindex(p) - 1, but in our case, the sign
                                    // shoud be +1, since we number the bonds from 0, whereas he numbers them
                                    // starting from 1
    vector<double> pInsertDiagonal; // pInsertDiagonal[n] = probability of accepting a diagonal
                                    // insertion update when the Hamiltonian is on the power n
    vector<double> pDeleteDiagonal; // pDeleteDiagonal[n] = probability of accepting a diagonal
                                    // deletion update when the Hamiltonian is on the power n
    vector<int> linkedVertexList;   // auxiliary array containing the connectivity of operator vertex legs.
                                    // Has as many elements as many vertex legs we have.
                                    // Note: v == linkedVertexList[linkedVertexList[v]]
    vector<int> vfirst;             // auxiliary array: vfirst[i] == the first operator vertex leg at site i
    vector<int> vlast;              // auxiliary array: vlast[i] == the last operator vertex leg at site i
    
public:
    // constructor and destructor
    heisenbergEquilibriumQMC(double Jval=1., double Tval=100.);
    ~heisenbergEquilibriumQMC();
    
    // reset data and reinitialize simulation
    void reset(double Jval, double Tval);
    
    // QMC routines
    void thermalize(int nBurnInSteps);
    void run(int nMonteCarloSteps);
    
    // print out all data, using the operator << of the parent class as well
    friend ostream& operator<<(ostream &, const heisenbergEquilibriumQMC &);
    
private:
    // initialize the private variables
    void initialize(void);
    
    // insertion or deletion of an S_z(i)*S_z(j) diagonal operator, see Sandwik's notes
    void diagonalUpdate(void);
    
    // flip all spins and operators in each operator loop with probability 1/2, see Sandwik's notes
    void operatorLoopUpdate(void);
    
    // auxiliary function to update the values of pInsertDiagonal and pDeleteDiagonal
    // needed when the maximal cutoff of powers of the Hamiltonian is changed
    void calculateDiagonalUpdateProbabilities(void);
    
    // auxiliary routine to update L during the thermalization process, if nmax becomes too large
    // According to Sandwik, the ideal cutoff is L = nmax + nmax/3
    void updatePowerCutoff(void);
    
    // update average energy and magnetic susceptibility
    void updateAvgEnergy(void);
    void updateMagneticSusceptibility(void);
    void updateSpecificHeat(void);
};


#endif /* MC_heisenbergEquilibriumQMC_hpp */
