//
//  MC_coord.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 9/19/16.
//
//
//  We perform our numerics on a two-dimensional lattice. There are two ways of indexing each lattice site:
//      - by its (x,y) coordinates:
//           (-N,N)   (-N+1,-N)     (-N+2,-N)  ... (N-1,-N)
//             ...       ...           ...     ...    ...
//             ...       ...           ...     ...    ...
//          (-N,-N+2) (-N+1,-N+2)  (-N+2,-N+2) ... (N-1,-N+2)
//          (-N,-N+1) (-N+1,-N+1)  (-N+2,-N+1) ... (N-1,-N+1)
//           (-N,-N)   (-N+1,-N)    (-N+2,-N)  ...  (N-1,-N)
//      - a single site index, identifying the sites:
//           4N   4N+1 4N+2 4N+3 4N+4  ...     6N-1
//           2N   2N+1 2N+2 2N+3 2N+4  ...     4N-1
//            0    1    2    3    4    ...     2N-1
//  where N=LATTICE_SIZE is the size of the sides of the 2D lattice
//
//  The functions below simply convert the site indices between these two representations


#ifndef MC_coord_hpp
#define MC_coord_hpp

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include "MC_clusterRunQ.hpp"


static const int LATTICE_SIZE = 1;     // global constant
                                       // size of the sides of the 2D lattice: x and y go from -LATTICE_SIZE to LATTICE_SIZE-1.
                                       // we use periodic boundary conditions

// modulo function (not identical to the symbol %, which does not work properly )
inline int mod(int x, int N) {
    return (x>=0) ? (x % N) : (N - 1 - ((-x-1) % N));
}

// converts x and y coordinates into site index
inline int coordToIndex(int x, int y) {
    return mod(x+LATTICE_SIZE, 2*LATTICE_SIZE) + mod(y + LATTICE_SIZE,2*LATTICE_SIZE)*(2*LATTICE_SIZE);
}

// converts site index into x coordinate
inline int indexToX(int ind) {
    return (mod(ind, 2*LATTICE_SIZE) - LATTICE_SIZE);
}

// converts site index into y coordinate
inline int indexToY(int ind) {
    return (mod(ind/(2*LATTICE_SIZE), 2*LATTICE_SIZE) - LATTICE_SIZE);
}

#endif /* MC_coord_hpp */
