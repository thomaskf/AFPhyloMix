//
//  phyloMCMC.h
//  HaploCount
//
//  Created by Thomas Wong on 11/2/19.
//

#ifndef phyloMCMC_h
#define phyloMCMC_h


#include <stdio.h>
#include <set>
#include "tree.h"
#include "phyloRead.h"
#include "mcmc.h"
#include "myRand.h"
#include "definitions.h"

class PhyloMCMC : public MCMC {
    
public:
    
    // to consider adjacent SNP positions
    int considerAdjSNPs;
    // 0 - No
    // 1 - Yes
    
    double freq_win_size;
    double edge_win_size;
    double err_win_size;
    double max_err;
    
    // the current frequencies of the haplotypes
    vector<double> haploFreqs_curr;

    // for the current topology
    vector<int> topologies_curr;
    vector<bool> isDescendant_curr;
    
    // the current log-likelihood
    long double logLike_curr;
    
    // the current log-prior
    long double logPrior_curr;
    
    // the current logPost
    long double logPost_curr;
    
    // the current state of the root for each position
    int* rootStates_curr;
    
    // the current mutation location for each position
    int* mutationLocs_curr;
    
    // the current expected Log frequencies and which edges have mutations
    double* expecLogFreqs_curr;
    
    // the current sequencing error
    double seqErr_curr;
    
    // the current function id and current k
    int func_id_curr, k_curr, pos_curr;

    // the proposed log-likelihood
    long double logLike_p;
    
    // the proposed log-prior
    long double logPrior_p;

    // the proposed logPost
    long double logPost_p;
    
    // prior for gamma
    long double freq_shape;
    long double freq_rate;
    long double edge_beta;
    
    // dimension
    int numFreqs;
    int numEdges;
    
    // random
    // MyRand myRand; (included in the base class)
    
    // log file
    string fileLog;
    // output file for logging sequences
    string fileSeq;
    // tree with max posterior probability
    string fileMaxTree;
    // sequences with max posterior probability
    string fileMaxSeq;
    
    // corresponding output file streams
    ofstream* flog;
    ofstream fseq;

    // the result with maximum posterior probability
    vector<string> maxSeqs;
    vector<double> maxFreqs;
    string maxTree;
    double maxErrRate;
    long double maxLogPost;
    bool maxResultExist; // false if there is no max result

    PhyloRead phyloRead;
    
    // constructor
    PhyloMCMC();
    
    // destructor
    ~PhyloMCMC();
    
    // get the information of the alignment
    void get_align_info(char* samFile, double& avgCover, int& numSNPs, int& seqLen);

    // load samfile
    void load_sam(char* samFile, int numHap, set<int>& show_pos);
    
    // load from other object
    void load_sam(PhyloMCMC* phyloMCMC, int numHap);
    
    // load genome file
    // void load_genomes(char* genomeFile, int numHap);
    
    // load genome file
    // void load_genomes(char* genomeFile, int numHap, vector<int>& hFreqs);
    
    // initialize the state
    // and get the log-posterior value of the starting state
    // and the number of steps in each generation
    void init(long double& currLnP, int& numSteps);
    
    // propose a new state
    // output the corresponding log-posterior and the hasting ratio
    void proposeNewState(int whichStep, long double& logPost, long double& logHastRatio);

    // check whether the files of previous run exist
    int checkWhetherPreRunFilesExist();
    
    // resume from the previous run
    void resumeFrPreRun();
    
    // resume from the previous log file
    // void resumefrPreLog(char* preLogFile);
    
    // load the starting tree file
    void loadStartTreeFile(char* startTreeFile);
    
    // accept the new state
    void acceptNewState(int whichStep);
    
    // reject the new state
    void rejectNewState(int whichStep);

    // print the result with max posterior probability
    void printMaxState();
    
    // print the current state
    void printState(int round, int needComputeLogLike);
    
    // print the new propose state  (without end of line character)
    // void printNewProposeState(int round);

    // run -- the main function to call for analysis (for sam file)
    void run(char* samFile, int numHap, int nGenerations, int nPrintSteps, int toConsiderAdj, char* startTreeFile, set<int>& show_pos, int toResume, int backMutateMethod, string& prefix);

    // run -- the main function to call for analysis (for genome file)
    // void runG(char* genomeFile, int numHap, int nGenerations, int nPrintSteps);

    // run -- the main function to call for analysis (for genome file)
    // void run(char* genomeFile, int numHap, vector<int>& hFreqs, int nGenerations, int nPrintSteps);

    // the log posterior for gamma distribution
    long double logPriorGamma(double x, double shape, double rate);

    // the log posterior for beta distribution
    long double logPriorBeta(double x, double beta);
    
    // compute the log-posterior of the current state / propose state
    void logPost(int isPropose);

    // do NNI
    // return 1 if there is change in the tree topology, but no change in the edge length array
    // return 2 if there are changes in both the tree topology and the edge length array
    int doNNI();

    // swap subtrees
    // return 0 if there is no change in the tree topology
    // return 1 if there is change in the tree topology, but no change in the edge length array
    // return 2 if there are changes in both the tree topology and the edge length array
    int swapSubTree();
    
    // distribute frequencies between two nodes
    void distributeFreqs(int& node1, int& node2);
    
    // swap frequencies between two nodes
    void swapFreqs(int& node1, int& node2);

    // merge a subtree of 2 tips into a single tip, and break one tip into a subtree of 2 tips
    int mergeAndBreak();

    // combine moves
    // randomly select two of actions: 0 - NNI; 1 - swapping subtrees; 2 - merge & break
    void combineMoves();
    
    // update the location of the mutations
    void updateMutateLocate(int vpos);

    // update the root state
    void updateRootState(int vpos);
    
    // update the haplotype frequency
    void updateHapFreq(int hap);

    // update the error rate
    void updateErrorRate();

    // add move
    void addMove(int func_id, int k, int weight);
    
    // save the best result so far
    void saveBestResult();

private:
    // for NNI computation
    vector<pair<int,int> > changeNodeIDs;
    vector<double> pr_tt;
    // for even frequency distribution
    int node1, node2;
    // for merge and break
    int node1a, node1b;
    int isChanged;
    // choice of moves
    vector<int> moveFuncPool;
    vector<int> moveKPool;
    // get the index of the node from the topology matrix
    int getIdxFrTop(int node);
    
};


#endif /* phyloMCMC_h */
