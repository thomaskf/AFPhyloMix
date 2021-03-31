//
//  mcmcmc.cpp
//  HaploCount
//
//  Created by Thomas Wong on 6/6/19.
//

#include "mcmcmc.h"

MCMCMC::MCMCMC(int numChains) {
    // constructor
    
    int i;
    
    if (numChains < 1) {
        cerr << "Error! Number of chains cannot be smaller than 1" << endl;
        exit(1);
    }
    
    for (i=0; i<numChains; i++) {
        mcmc_set.push_back(new PhyloMCMC());
    }
    nChain = numChains;
    
    // assign the thread IDs, temperature and isColdChain
    for (i=0; i<numChains; i++) {
        mcmc_set[i]->threadID = i+1;
        mcmc_set[i]->computeTemperature(numChains);
        mcmc_set[i]->isColdChain = 0;
    }
    mcmc_set[0]->isColdChain = 1;
    
}

MCMCMC::~MCMCMC() {
    // destructor
    int i;
    for (i=0; i<nChain; i++) {
        delete mcmc_set[i];
    }
    mcmc_set.clear();
}

// load samfile
void MCMCMC::load_sam(char* samFile, int numHap, set<int>& show_pos) {
    
    int i;
    mcmc_set[0]->load_sam(samFile, numHap, show_pos);
    for (i=1; i<nChain; i++) {
        mcmc_set[i]->load_sam(mcmc_set[0], numHap);
    }
}

// run -- the main function to call for analysis (for sam file)
// run -- the main function to call for analysis (for sam file)
void MCMCMC::run(char* samFile, int numHap, int nGenerations, int nPrintSteps, int toConsiderAdj, char* startTreeFile, set<int>& show_pos, int toResume, int backMutateMethod, string& prefix, int nCheckChains) {
    
    int i,k;
    string inputFileName;
    vector<int> resumes;
    
    cout << "Mathematical method for likelihood computation regarding the mutation location: ";
    if (METHOD_EQU == 1) {
        cout << "MAX";
    } else {
        cout << "SUM";
    }
    cout << endl;
    
    // setting parameters
    for (i=0; i<nChain; i++) {
        mcmc_set[i]->phyloRead.back_mutate_method = backMutateMethod;
        mcmc_set[i]->considerAdjSNPs = toConsiderAdj;
        mcmc_set[i]->fileLog = prefix + "." + int2str(i) + ".log";
        mcmc_set[i]->fileMaxTree = prefix + ".max.tree";
    }

    // check whether the files of previous run exist
    if (toResume) {
        for (i=0; i<nChain; i++) {
            resumes.push_back(mcmc_set[i]->checkWhetherPreRunFilesExist());
        }
    }

    // load the genome file
    cout << "loading the sam/bam file ... " << endl;
    load_sam(samFile, numHap, show_pos);

    // resume from the previous run
    if (toResume) {
        for (i=0; i<nChain; i++) {
            if (resumes[i]) {
                mcmc_set[i]->resumeFrPreRun();
            }
        }
    }
    
    // initialization
    init();

    // open the output file streams
    for (i=0; i<nChain; i++) {
        mcmc_set[i]->flog = new ofstream();
        mcmc_set[i]->flog->open(mcmc_set[i]->fileLog.c_str(), ofstream::out | ofstream::app);
    }

    // allow all chains to reach the same place
    k = getmaxstartIteration()-1;
    if (k > 0) {
        cout << "Number of generations being processed in the previous run: " << k << endl;
        runAllThreads(k, nPrintSteps);
    }
    
    // start MCMC from k-th generation
    cout << "Start MCMCMC" << endl << flush;
    while (k < nGenerations) {
        k += nCheckChains;
        if (k > nGenerations)
            k = nGenerations;
        runAllThreads(k, nPrintSteps);
        if (k < nGenerations)
            checkBetweenChains();
    }
    
    // close the output file streams
    for (i=0; i<nChain; i++)
        mcmc_set[i]->flog->close();
}

// initialization
void MCMCMC::init() {
    int i;
    thread_set.resize(nChain);
    args_threads.resize(nChain);
    for (i=0; i<nChain; i++) {
        args_threads[i].phyloMCMC = mcmc_set[i];
    }
}

// obtain the max startIteration
int MCMCMC::getmaxstartIteration() {
    int i;
    int k = mcmc_set[0]->startIteration;
    for (i=1; i<nChain; i++) {
        if (k < mcmc_set[i]->startIteration) {
            k = mcmc_set[i]->startIteration;
        }
    }
    return k;
}

void * runThread(void * arg) {
    ThreadArguments* threadArg = (ThreadArguments*) arg;
    threadArg->phyloMCMC->processMCMC(threadArg->nGenerations, threadArg->nPrintSteps);
    pthread_exit ( 0 );
    return 0;
}

// run multi-thread MCMCMC
void MCMCMC::runAllThreads(int nGenerations, int nPrintSteps) {
    int threadId, i;
    for (i=0; i<nChain; i++) {
        args_threads[i].nGenerations = nGenerations;
        args_threads[i].nPrintSteps = nPrintSteps;
    }
    for (threadId=0; threadId<nChain; threadId++) {
        if ( pthread_create ( & ( thread_set[threadId] ), NULL, runThread, ( void * ) & ( args_threads[threadId] ) ) ) {
            cerr << "[mcmcmc: ThreadId " << threadId << "] Cannot create ""runThread""\n" << endl << flush;
            exit(1);
        }
    }
    // args_threads[0].phyloMCMC->processMCMC(args_threads[0].nGenerations, args_threads[0].nPrintSteps);
    for (threadId=0; threadId<nChain; threadId++) {
        if ( thread_set[threadId] != 0 ) {
            if ( pthread_join ( thread_set[threadId], NULL ) ) {
                cerr << "[performRAL: ThreadId " << threadId << "] Crash!\n" << endl << flush;
                exit(1);
            }
            // pthread_detach(thread_set[threadId]);
            thread_set[threadId] = 0;
        }
    }
}

void MCMCMC::swapColdChain(int chainJ) {
    // swap j-chain with cold chain
    ofstream* tflog;
    long double ttemp;
    PhyloMCMC* t;
    
    // swap log files
    tflog = mcmc_set[0]->flog;
    mcmc_set[0]->flog = mcmc_set[chainJ]->flog;
    mcmc_set[chainJ]->flog = tflog;
    
    // swap temperature
    ttemp = mcmc_set[0]->temperature;
    mcmc_set[0]->temperature = mcmc_set[chainJ]->temperature;
    mcmc_set[chainJ]->temperature = ttemp;
    
    // swap the chain type
    mcmc_set[0]->isColdChain = 0;
    mcmc_set[chainJ]->isColdChain = 1;
    
    // swap the chain ID
    mcmc_set[0]->threadID = chainJ;
    mcmc_set[chainJ]->threadID = 0;
    
    // transfer the result with maximum posterior probability to the new chain
    mcmc_set[chainJ]->maxResultExist = mcmc_set[0]->maxResultExist;
    if (mcmc_set[chainJ]->maxResultExist) {
        mcmc_set[chainJ]->maxFreqs.clear();
        mcmc_set[chainJ]->maxFreqs.insert(mcmc_set[chainJ]->maxFreqs.begin(), mcmc_set[0]->maxFreqs.begin(), mcmc_set[0]->maxFreqs.end());
        mcmc_set[chainJ]->maxTree = mcmc_set[0]->maxTree;
        mcmc_set[chainJ]->maxErrRate = mcmc_set[0]->maxErrRate;
        mcmc_set[chainJ]->maxLogPost = mcmc_set[0]->maxLogPost;
    }

    // swap the object in the array
    t = mcmc_set[0];
    mcmc_set[0] = mcmc_set[chainJ];
    mcmc_set[chainJ] = t;
}

// randomly select one of the hot chains, compare its posterior probability with the one
// in cold chain, and swap with the two chains if necessary
void MCMCMC::checkBetweenChains() {
    long double logpost0 = mcmc_set[0]->logPost_curr;
    long double t0 = mcmc_set[0]->temperature;
    int j = mcmc_set[0]->myRand.iunif(1, nChain-1);
    long double logpost1 = mcmc_set[j]->logPost_curr;
    long double t1 = mcmc_set[j]->temperature;
    long double lnAlpha = (t0 - t1) * (logpost1 - logpost0);
    if ((lnAlpha >= 0) || (log(mcmc_set[0]->myRand.runif()) < lnAlpha)) {
        // swap between the cold chain and chain j
        cout << "chain " << j << " swap with the cold chain" << endl;
        swapColdChain(j);
    }
}

