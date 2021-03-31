//
//  mcmcmc.h
//  HaploCount
//
//  Created by Thomas Wong on 6/6/19.
//

#ifndef mcmcmc_h
#define mcmcmc_h

#include <cstdlib>
#include <pthread.h>
#include <vector>
#include <set>
#include <pthread.h>
#include "phyloMCMC.h"
#include "definitions.h"

using namespace std;

typedef struct ThreadArguments {
    PhyloMCMC* phyloMCMC;
    int nGenerations;
    int nPrintSteps;
} ThreadArguments;

// metropolis coupled
class MCMCMC {
public:
    
    vector<PhyloMCMC*> mcmc_set;
    int nChain;
    vector<pthread_t> thread_set;
    vector<ThreadArguments> args_threads;

    MCMCMC(int numChains); // constructor
    
    ~MCMCMC(); // destructor
    
    // load samfile
    void load_sam(char* samFile, int numHap, set<int>& show_pos);

    // run -- the main function to call for analysis (for sam file)
    void run(char* samFile, int numHap, int nGenerations, int nPrintSteps, int toConsiderAdj, char* startTreeFile, set<int>& show_pos, int toResume, int backMutateMethod, string& prefix, int nCheckChains);

    // initialization
    void init();
    
    // run multi-thread MCMCMC
    void runAllThreads(int nGenerations, int nPrintSteps);

    // swap j-chain with cold chain
    void swapColdChain(int chainJ);

    // randomly select one of the hot chains, compare its posterior probability with the one
    // in cold chain, and swap with the two chains if necessary
    void checkBetweenChains();

    // obtain the max startIteration
    int getmaxstartIteration();
    
};

#endif /* mcmcmc_h */
