//
//  phyloRead.h
//  HaploCount
//
//  Created by Thomas Wong on 15/1/19.
//

#ifndef phyloRead_h
#define phyloRead_h

#include <cstdlib>
#include <cstring>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include "reads.h"
#include "tree.h"
#include "gtree.h"
#include "myRand.h"
#include "definitions.h"

using namespace std;

struct LogFreqs {
    double v[16]; // log-freqs
    // Symbol: 0 - root state; 1 - another state; 2 - error state; 3 - another error state
    // The sixteen states are:
    // 0 0
    // 0 1
    // 0 2
    // 0 3
    // 1 0
    // 1 1
    // 1 2
    // 1 3
    // 2 0
    // 2 1
    // 2 2
    // 2 3
    // 3 0
    // 3 1
    // 3 2
    // 3 3
};

class PhyloRead {
public:
    
    int numTips;
    int numNodes;
    int numInNodes; // number of internal nodes
    int numEdges; // number of edges in a rooted tree

    // the current tree topology
    vector<int> currTopology;

    // relationship between two nodes
    // isDescendant[i*numNodes+j] = true if node j is descendant of node i
    vector<bool> isDescendant;
    vector<int> refPos; // reference position for each position
    vector<int> refLst; // list of the reference positions

	// reads
	Reads rds;

    // data matrix (for each single position)
    int* sdMat;
    // data matrix (for adjacent snp positions)
    int* dMat;
    // dimension: numPos x 16 x 4
    int numPos; // number of positions
    int* charOrder; // the characters (i.e. 0,1,2,3) in decreasing order for each position
    // dimension: numPos * 4
    vector<long double> logStateFreqs; // freqs of A, C, G, and T

    bool* isInVar; // whether the site is invariable
    bool* isDiscard; // whether the site is discarded

    // vector<pair<int,int> > winLists; // start and end positions for each window

    long double logCombin1; // for no considering adjacent snp positions
    
    long double logCombin;
    
    vector<int> variablePos; // variable positions
    int numVarPos; // number of variable positions
    
    //--------------------------------------
    // Parameters
    //--------------------------------------

    // the frequencies of the haplotypes
    vector<double> haploFreqs_t;
    
    // the frequencies of the haplotypes
    vector<double> haploFreqs; // normalised
    
    // the weights of each node
    vector<double> nodeWeights;
    
    // the state of the root for each position
    int* rootStates;
    
    // the mutation location for each position
    int* mutationLocs;

    // the expected Log frequencies and which edges have mutations
    double* expecLogFreqs; // for two snp positions
    double* expecLogFreq1; // for one snp position
    
    // sequencing error
    double seqErr;
    
    // edge lengths
    vector<long double> edgeLen;
    vector<long double> logEdgeLen;

    //--------------------------------------
    // status
    //--------------------------------------
    
    bool isDescendantUpdated;
    bool isHapFreqNormalized;
    bool isWeightUpdated;
    bool isLogFreqUpdated;
    bool needUpdated;
    
    //--------------------------------------
    // additional parameters for likelihood calculations
    //--------------------------------------
    
    double* log_phi;
    // dimension: numPos * numEdges * numEdges * 2
    // phi[i*(numEdges*numEdges*2) + y*numEdges*2 + x*2 + a = phi(i,y,a,x)
    //
    // for i != ref-pos
    //     phi(i,y,a,x) = sum_b L(D | b,y,a,x) * Pr(b)
    // for i = ref-pos
    //     phi(i,a, x) = sum_b L(D | a,x) * Pr(a)
    // where 0<=x,y<=2N-3 : indicate the mutation location on the position
    //       a,b \in {A,C,G,T} : indicate the root state on the position
    
    double* lpr_ui;
    // LOG Probability of mutation location on each position
    // i.e. lpr_ui[i*numEdges + j] = log prob of mutation location = j on pos i
    
    //--------------------------------------
    // parameters for haplotype construction
    //--------------------------------------

    double* charMatHaplo;
    // the log-probability matrix for characters A,C,G,T for each position of each haplotype
    // dimension: numPos * numTips * 4
    // charMatHaplo[i * numTips * 4 + j * 4 + k] = the log-prob of char k on position i of haplotype j
    
    char* isValidCharMat;
    // referring to the above matrix
    // 1 if the probability is non-zero, 0 otherwise
    
    double* lprobRcMc;
    // a temporary array for computation of the log-probability Pr(r_c, m_c | W) for the ref pos c inside a window W

    //------------------------------------------------------
    // which method to estimate the back-mutation positions
    //------------------------------------------------------

    int back_mutate_method;
    
    //------------------------------------------------------
    // additional parameters for MCMCMC
    //------------------------------------------------------

    // Is using the matrices of other object?
    bool useOtherMat;
    
    //--------------------------------------
    // functions
    //--------------------------------------

    // constructor
    PhyloRead();
    // PhyloRead(int seed_t);
    
    // destructor
    ~PhyloRead();
    
    // build the array "isDescendant" that indicates the relationship between two nodes
    // isDescendant[i*numNodes+j] = true if node j is descendant of node i
    void buildIsDescendant();

    // randomly generate a random tree
    void computeRandTopology(int nTips, MyRand& myRand);

    // (for debugging) to get the actual back mutation snp positions
    void getActualBackMutate(char* samFile, vector<vector<int> >& c_sets, int numHaps);
    
    // from the genome alignments
    // compute the data matrix
    void loadData(char* samFile, set<int>& show_pos);
    
    // from another copy of phyloRead
    // the data matrix will be pointing to this object
    void loadData(PhyloRead& aPhyloRead);
    
    // get alignment information
    // Average coverage, number of SNP positions, and sequence length
    void getAlgnInfo(char* samFile, double& avgCover, int& numSNPs, int& seqLen);

    // (for debugging) check the snp positions
    void checkSNPPos(vector<int>& snpPosSet, vector<vector<int> >& c_sets, int numHaps);
    
    // identify the bad mutations
    void findBadMutations(vector<int>& snpPosSet, vector<pair<int,int> >& snpChars, set<int>& show_pos);
    
    // set all parameters to random values
    void randParameters(MyRand& myRand);

    // set up the variables:
    void setupVariables();
    
    // normalize haplotype frequencies
    void normalizeHapFreqs();

    // update the weight of the nodes
    // according to the updated haplotype frequencies
    // assume the haplotype frequencies to be normalized
    void updateWeights();

	// update the expected frequencies for mutations on different pairs of edges
	// update the arrays expectLogFreqs
	// expectLogFreqs[i*numEdges + j] = the expected log frequencies when the mutations happens on edges i+1 and j+1
	void updateExpectLogFreqs();

    // update the expected frequencies for mutations on different pairs of edges
    // update the arrays expectLogFreqs
    // expectLogFreqs[i*numEdges + j] = the expected log frequencies when the mutations happens on edges i+1 and j+1
    void updateExpectLogFreq1();

    // compute the logCombin
    void computelogCombin(int* dmatrix);

    // compute the value of logCombin (for data set without considering adjacent snp positions)
    void computelogCombin1(int* sdmatrix);

    // compute log-likelihood of the whole alignment without consider the adjacent snp positions
    long double logLikeS();
    // obtainMutatePos=1 if want to obtain the optimal mutation positions
    long double logLikeS(int obtainMutatePos);

    // compute log-likelihood of the whole alignment
    long double logLike();

    // compute log-likelihood of the whole alignment
    // obtainMutatePos=1 if want to obtain the optimal mutation positions
    long double logLike(int obtainMutatePos);
    
    // compute the log-likelihood of the whole alignment
    // by considering the edge lengths
    long double logLike2();
    
    // compute the array phi
    void computePhi();
    
    // compute the edge lengths
    void computeEdgeLen();

    // show the current tree
    string currTree(int num_digits);
    
    // show the current tree from the array edgeLen
    string currTree2(int num_digits);

    // show the current tree with all edge length = 1.0
    string currTree1(int num_digits);

    // get the current sequences
    void getCurrSeqs(vector<string>& seqs);
    
    // show the current sequences
    void showCurrSeqs(ofstream& fout, int round);
    
    // compute the loglikelihood given mi and ri for all positions
    void computeLogLikeRiMi();
    
    // compute the log-probability matrix for character 0 and 1 for each position of each tip (i.e. haplotype)
    void computeCharMatHaplo();

    // report the set of characters with maximum probability for each haplotype
    void reportMaxCharSet(vector<string>& haplos);
    
    // report the set of characters for each haplotype
    void reportHaploSeqs(vector<string>& haplos);

private:

    int* nucl2int;
    
    // for display purpose
    UnlabelTree tree;
    vector<double> edgelens;
    
    // for computation of loglikelihoods
    long double* tloglike;
    // int* tloglikeValid;
    // dimension: numEdges * 2
    
    // for storing the maximum of c2 (i.e. the root state on the current position)
    int* maxStateRoot;
    // dimension: numPos * numEdges * 2
    // maxStateRoot[pos * numEdges * 2 + m1*2 + c1] = the optimal root state on pos
    // when the mutation location on the ref position = m1 and the root state on the ref position = c1
    
    int* stateRoot;
    // dimension: numPos
    // stateRoot[pos] = the optimal root state on pos

    // for storing the maximum of m2 (i.e. the mutation location on the current position)
    int* maxMutatePos;
    // dimension: numPos * numEdges * 2
    // maxMutatePos[pos * numEdges * 2 + m1*2 + c1] = the optimal mutation location on pos
    // when the mutation location on the ref position = m1 and the root state on the ref position = c1
    
    int* mutatePos;
    // dimension: numPos
    // mutatePos[pos] = the optimal mutation location on pos
    
    double* logLikeRiMi;
    // dimension: numPos * numEdges * 2
    // logLikeRiMi[pos * numEdges * 2 + ri * numEdges + m1] = the log-likelihood given mi and ri for the position pos

    void startup();

    // compute whether the adjacent columns should be combined
    // should the adjacent positions be combined
    // vector<int> combinePos;
    // if position i can combine with position i+1, then combinePos[i] = 1
    // if position i can combine with position i+2, then combinePos[i] = 2
    // otherwise, combinePos[i] = 0
    void getCombinePos(Reads& rds, int numPos, vector<int>& combinePos);
    
    // get information of the data
    // construct the arrays chars0, chars1, isInVar, isDiscard
    void getDataInfo(int* posMatrix, vector<int>& coverage, vector<int>& coverage_with_gap);

    // characters of each tip
    // 0 - same as the root state; 1 - different from the root state
    // size of tipChars is initialized as the number of tips
    void getChars(int mutateLoc, vector<int>& tipChars);

    // obtain the state frequencies from posMatrix
    void getStateFreqs(int* posMatrix);

    // compute LOG Pr(u_c = x | theta)
    // i.e. LOG Pr(location of mutation on the reference position c = x | the parameters)
    void computeLPrUc();
    
    // compute LOG Pr(u_i = x | theta)
    // i.e. log Pr(location of mutation on the position i = x | the parameters)
    void computeLPrUi();

    // get the bad positions
    // pairstate: nxn dimension; 1 - good; 2 - bad; 0 - not enough coverage
    void selectBadPos(vector<char>& pairstate, int n, vector<int>& badPos, vector<double>& badPairRatio, vector<int>& numPair);
};


#endif /* phyloRead_h */
