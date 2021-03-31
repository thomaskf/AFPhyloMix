//
//  mcmc.h
//  mcmc
//
//  Created by Thomas Wong on 30/11/18.
//  Copyright © 2018 Thomas Wong. All rights reserved.
//

#ifndef mcmc_h
#define mcmc_h

// show detailed steps
// #define SHOW_DETAILS

#include <iostream>
#include <cstdlib>
#include "myRand.h"
#include "definitions.h"

using namespace std;

class MCMC {
    
public:

    // random
    MyRand myRand;
    
    // thread ID
    int threadID;
    
    // isColdChain
    int isColdChain;
    
    // temperature
    long double temperature;
    
    long double currLnP;
    int nSteps;
    int newStart;

    // constructor
    // also initialize the seed
    MCMC();
    
    //destructor
    virtual ~MCMC() {}
    
    // do/continue MCMC
    // parameters: number of generations
    void processMCMC(int nGenerations, int nPrintSteps);
    void processMCMC(int nGenerations, int nPrintSteps, int newStart);

    // initialize the state
    // and get the log-posterior value of the starting state
    // and the number of steps in each generation
    virtual void init(long double& currLnP, int& numSteps) {}
    
    // propose a new state
    // output the corresponding log-posterior and the hasting ratio
    virtual void proposeNewState(int whichStep, long double& logPost, long double& logHastRatio) {}

    // accept the new state
    virtual void acceptNewState(int whichStep) {}

    // reject the new state
    virtual void rejectNewState(int whichStep) {}
    
    // print the current state (without end of line character)
    virtual void printState(int round, int needComputeLogLike) {}
    
    // print the new propose state  (without end of line character)
    virtual void printNewProposeState(int round) {}

    // resume from the previous log file
    virtual void resumefrPreLog(char* preLogFile) {}

    // save the best result so far
    virtual void saveBestResult() {}

    // print the best result
    virtual void printMaxState() {}

// private:

    // how many times the proposal is accepted
    int numAccepts;
    
    // total number of generations
    // int totNumGenerations;
    
    // starting iteration
    int startIteration;
    
    // compute the temperature
    void computeTemperature(int numChains);
};

#endif /* mcmc_h */
