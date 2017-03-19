//
//  MC_main.cpp
//  
//
//  Created by Marton Kanasz-Nagy on 9/19/16.
//
//

#include "MC_main.hpp"


// Arguments of the main function:
// 1. T/J = temperature in units of the coupling
// 2. nMCruns = number of independent Monte Carlo runs we start (usually 1-100)
// 3. nMCsteps = number of Monte Carlo steps per run (usually 100K)
// 4. nThermalizationSteps = number of thermalization steps (Monte Carlo steps during which we do not measure quantities at the beginning of the simulation)
int main(int argc,char *argv[])
{
    double J = 1.;
    double T;
    int nMCruns;
    int nMCsteps;
    int nThermalizationSteps;
    
    
    // Check if there is enough arguments
    if(argc!=5) {
        cerr << "Error in main(): Incorrect number of arguments." << endl;
        return 1;
    }
    
    // Get inputs
    T = atof(argv[1]) * J;
    nMCruns = atol(argv[2]);
    nMCsteps = atol(argv[3]);
    nThermalizationSteps = atol(argv[4]);
    
/*
    // Set up spin correlation measurement
    vector< vector<double> > spinCorr;
    int nSites = 4*LATTICE_SIZE*LATTICE_SIZE;
    vector< pair<int,int> > sitePairs(nSites*nSites);
    int i,j,k=0;
    for(i=0; i<nSites; i++) {
        for(j=0; j<nSites; j++) {
            sitePairs[k].first = i;
            sitePairs[k].second= j;
            k++;
        }
    }
*/
    // Set up simulation
//    isingEquilibriumMC simulation(J, T);
    heisenbergEquilibriumQMC simulation(J, T);
    simulation.measureAllQuantities();
    // simulation.measureStaggeredMagnetizationFullCountingStatistics(true, 2*LATTICE_SIZE*LATTICE_SIZE);
//    simulation.measureSpinZCorrelators(true, sitePairs);
    
    // Run quantum Monte Carlo
    simulation.thermalize(nThermalizationSteps);
    simulation.run(nMCsteps);
//    spinCorr = simulation.getSpinZCorrelators();
    
    cout << simulation << endl;
    
    histogram hist = simulation.getStaggeredMagnetizationFullCountingStatistics();
    hist.printProbabilities();
    

    // Print out full counting statistics of staggered magnetization:
//    cout << "Full counting statistics of the staggered magnetization: " << endl;
//    cout << simulation.getStaggeredMagnetizationFullCountingStatistics() << endl;

    
/*
 // Print out full counting statistics of magnetization:
 cout << "Full counting statistics of the staggered magnetization: " << endl;
 cout << simulation.getStaggeredMagnetizationFullCountingStatistics() << endl;
*/
    
/*
 // Print out spin correlations
 for(i=0; i<spinCorr.size(); i++) {
 cout << spinCorr[i][0] << "\t" << spinCorr[i][1] << "\t" << spinCorr[i][2] << "\t" << spinCorr[i][3] << "\t" << spinCorr[i][4] << endl;
 }
*/
/*
    // Print out magnetic susceptibility
    cout << simulation.getMagneticSusceptibility() << endl;
 */
    return 0;
}

/*
// Measure staggered magnetization full counting statistics
 double J = 1.;
 double T;
 int nMCruns;
 int nMCsteps;
 int nThermalizationSteps;
 
 
 // Check if there is enough arguments
 if(argc!=5) {
 cerr << "Error in main(): Incorrect number of arguments." << endl;
 return 1;
 }
 
 // Get inputs
 T = atof(argv[1]) * J;
 nMCruns = atol(argv[2]);
 nMCsteps = atol(argv[3]);
 nThermalizationSteps = atol(argv[4]);
 
 // Set up simulation
 heisenbergEquilibriumQMC simulation(J, T);
 histogram fullCountingStatistics(-0.5, 0.5, 2*LATTICE_SIZE*LATTICE_SIZE);
 simulation.measureAllQuantities();
 simulation.measureStaggeredMagnetizationFullCountingStatistics(true, 2*LATTICE_SIZE*LATTICE_SIZE);
 
 // Run quantum Monte Carlo
 simulation.thermalize(nThermalizationSteps);
 simulation.run(nMCsteps);
 fullCountingStatistics = simulation.getStaggeredMagnetizationFullCountingStatistics();
 simulation.clearMeasurements();
 
 for(int i=1; i<nMCruns; i++) {
 simulation.thermalize(nThermalizationSteps);
 simulation.run(nMCsteps);
 fullCountingStatistics += simulation.getStaggeredMagnetizationFullCountingStatistics();
 simulation.clearMeasurements();
 }
 
 cout << "full counting statistics: " << endl << fullCountingStatistics << endl;
*/
