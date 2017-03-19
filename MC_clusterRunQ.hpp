//
//  MC_clusterRunQ.hpp
//  
//
//  Created by Marton Kanasz-Nagy on 10/20/16.
//
//
//  This file contains a single pragma, that tells the compiler whether I want to compile the code on my laptop
//  or on the cluster. The reason for this solution is that the cluster has a different version of gcc giving all
//  sorts of errors, whereas my Mac compiles the function seemlessly. This is especially true for C++11 features.
//  I use this pragma to erease those parts of the code the cluster has problems with. This is a cheap solution,
//  but I don't have time to rewrite all my code in a cluster-friendly manner.

#ifndef MC_clusterRunQ_hpp
#define MC_clusterRunQ_hpp


#define RUN_ON_CLUSTER



#endif /* MC_clusterRunQ_hpp */
