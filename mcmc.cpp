//
//  mcmc.cpp
//  mcmc
//
//  Created by Thomas Wong on 3/12/18.
//  Copyright © 2018 Thomas Wong. All rights reserved.
//

#include "mcmc.h"

// constructor
// also initialize the seed
MCMC::MCMC() {
    numAccepts = 0;
    startIteration = 1;
    temperature = 1.0;
    isColdChain = 1;
    threadID = 1;
    newStart = 1;
}

/*
// for debugging
// do/continue MCMC
// parameters: number of generations
void MCMC::processMCMC(int nGenerations, int nPrintSteps) {
    
    int i,k;
    int nSteps;
    long double currLnP;
    long double newLnP;
    long double lnHastRatio;
    long double lnAlpha;
    // double unitRand;
    // double lnUnitRand;
    bool stateAccepted;

    // initial the state
    // and get the log-posterior value of the starting state
    // and the number of steps in each generation
    init(currLnP, nSteps);
    cout << 0 << "\t"; printState(); cout << endl;
    
    for (i=0; i<=nGenerations; i++) {
        // cout << "generation " << i << endl << flush;
        for (k=0; k<nSteps; k++) {
            proposeNewState(k, newLnP, lnHastRatio);
            lnAlpha = newLnP - currLnP + lnHastRatio;
            stateAccepted = false;
            cout << i*nSteps + k + 1 << "\t"; printNewProposeState(); cout << endl;
            if ((lnAlpha >= 0) || (log(myRand.runif()) < lnAlpha)) {
                // accept the proposal
                acceptNewState(k);
                currLnP = newLnP;
                numAccepts++;
                stateAccepted = true;
            } else {
                // reject the proposal
                rejectNewState(k);
            }
        }
    }
    
    // print the last state
    if ((i-1) % nPrintSteps != 0) {
        cout << i-1 << "\t"; printState(); cout << endl;
    }
    
    // print out the ratio
    cout << "Acceptance proportion:" << (double) numAccepts / (nGenerations * nSteps) << endl;
}
*/

// do/continue MCMC
// parameters: number of generations
void MCMC::processMCMC(int nGenerations, int nPrintSteps) {
    
    int i,k;
    long double newLnP;
    long double lnHastRatio;
    long double lnAlpha;
    // double unitRand;
    // double lnUnitRand;
    bool stateAccepted;
    
    // initial the state
    // and get the log-posterior value of the starting state
    // and the number of steps in each generation
    if (newStart==1) {
        init(currLnP, nSteps);
        newStart = 0;
        printState(startIteration-1,1);
    }
    // cout << "startIteration = " << startIteration << endl;
    // cout << "nGenerations = " << nGenerations << endl;
    for (i=startIteration; i<=nGenerations; i++) {
        // cout << "generation " << i << endl << flush;
        for (k=0; k<nSteps; k++) {
            // cout << "k=" << k << endl << flush;
            proposeNewState(k, newLnP, lnHastRatio);
            lnAlpha = temperature * (newLnP - currLnP) + lnHastRatio;
            stateAccepted = false;
            if ((lnAlpha >= 0) || (log(myRand.runif()) < lnAlpha)) {
                // accept the proposal
                // cout << "accept the proposal" << endl << flush;
                acceptNewState(k);
                currLnP = newLnP;
                numAccepts++;
                stateAccepted = true;
            } else {
                // reject the proposal
                // cout << "reject the proposal" << endl << flush;
                rejectNewState(k);
            }
            
#ifdef SHOW_DETAILS
            printNewProposeState(i*nSteps + k);
            if (stateAccepted)
                cout << " ACCEPTED" << endl;
            else
                cout << " NOT" << endl;
#endif

            // save the best result so far
            // cout << "save the best result so far" << endl << flush;
            if (isColdChain) {
                saveBestResult();
            }
        }
        
        if (i > 0 && (i % nPrintSteps == 0)) {
            // print the current state
            printState(i,1);
            if (isColdChain) {
                // print the best result
                printMaxState();
            }
        }
    }
    
    /*
    // print the last state
    if ((i-1) % nPrintSteps != 0) {
        printState(i-1,0);
        if (isColdChain) {
            // print the best result
            printMaxState();
        }
    }
    */
    
    startIteration = nGenerations+1;
    
    // print out the ratio
    // cout << "Acceptance proportion:" << (double) numAccepts / (nGenerations * nSteps) << endl;
}

// compute the temperature
void MCMC::computeTemperature(int numChains) {
    
    // double delta_t = 1.0;
    // temperature = 1.0 / (1.0 + (threadID - 1.0) * delta_t);
    temperature = 1.0 - ((threadID - 1.0) / (double) numChains);
    
    cout << "Thread " << threadID << "'s temperature = " << temperature << endl;
}
