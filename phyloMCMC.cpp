//
//  phyloMCMC.cpp
//  HaploCount
//
//  Created by Thomas Wong on 11/2/19.
//

#include "phyloMCMC.h"

// check balanced parentheses
bool isBalParentheses(string& s) {
	int i,k;
	k=0;
	for (i=0; i<s.length(); i++) {
		if (s[i] == '(') {
			k++;
		} else if (s[i] == ')') {
			if (k > 0)
				k--;
			else
				return false;
		}
	}
	if (k == 0)
		return true;
	else
		return false;
}

//==============================================================
// for MCMC
//==============================================================

// constructor
PhyloMCMC::PhyloMCMC() : MCMC() {
    // the current state of the root for each position
    rootStates_curr = NULL;
    // the current mutation location for each position
    mutationLocs_curr = NULL;
    // the current expected Log frequencies and which edges have mutations
    expecLogFreqs_curr = NULL;
    // max result not exist yet
    maxResultExist = false;
    flog = NULL;
}

// destructor
PhyloMCMC::~PhyloMCMC() {
    // the current state of the root for each position
    if (rootStates_curr != NULL) {
        delete[] rootStates_curr;
    }
    // the current mutation location for each position
    if (mutationLocs_curr != NULL) {
        delete[] mutationLocs_curr;
    }
    // the current expected Log frequencies and which edges have mutations
    if (expecLogFreqs_curr != NULL) {
        delete[] expecLogFreqs_curr;
    }
    if (flog != NULL) {
        delete flog;
    }
}

// load the starting tree file
void PhyloMCMC::loadStartTreeFile(char* startTreeFile) {
    
    ifstream fin;
    string strTree;
    string strErr;
    bool isTreeStr;
    UnlabelTree startTree;
    double maxErr;
    bool isErrInfoAvail;
    vector<int> topInt;
    vector<double> haploFreq;
    vector<double> edgeLen;
    
    if (!isColdChain)
        return;
    
    // the prefix of the line with error rate
    string prefix = "# error rate: ";
    int prefix_len = prefix.length();
    int i;
    
    isErrInfoAvail = false;
    cout << "Loading the starting tree file: " << startTreeFile << endl;
    // read the max-tree-file
    // cout << "read the max-tree-file" << endl << flush;
    fin.open(startTreeFile);
    isTreeStr = false;
    while ((!isTreeStr) && getline(fin,strTree)) {
        if (strTree.length() == 0)
            continue;
        // get the corresponding error rate
        if (strTree.length() > prefix_len && strTree.substr(0,prefix_len) == prefix) {
            strErr = strTree.substr(prefix_len);
            cout << "Error rate: " << strErr << endl;
            maxErr = atof(strErr.c_str());
            isErrInfoAvail = true;
        }
        // read the tree
        else if (strTree[0] != '#') {
            isTreeStr = true;
            startTree.topTxt2Int(strTree, topInt, haploFreq, edgeLen);
        }
    }
    fin.close();
    if (!isTreeStr) {
        cerr << "Error! The file: " << startTreeFile << " cannot be read" << endl;
        exit(1);
    }
    
    if (phyloRead.numTips != haploFreq.size()) {
        cerr << "Error! The number of tips in the starting tree file does not match with the input number of haplotypes" << endl;
        exit(1);
    }

    phyloRead.currTopology.clear();
    for (i=0; i<topInt.size(); i++)
        phyloRead.currTopology.push_back(topInt[i]);
    phyloRead.isDescendantUpdated = false;
    phyloRead.needUpdated = true;
    phyloRead.buildIsDescendant();
    for (i=0; i<phyloRead.numTips; i++)
        phyloRead.haploFreqs_t[i] = haploFreq[i];
    if (isErrInfoAvail)
        phyloRead.seqErr = maxErr;
    else
        phyloRead.seqErr = 0.002;
}

// run -- the main function to call for analysis (for SAM file)
void PhyloMCMC::run(char* samFile, int numHap, int nGenerations, int nPrintSteps, int toConsiderAdj, char* startTreeFile, set<int>& show_pos, int toResume, int backMutateMethod, string& prefix) {
    
    if (isColdChain) {
        cout << "Mathematical method for likelihood computation regarding the mutation location: ";
        if (METHOD_EQU == 1) {
            cout << "MAX";
        } else {
            cout << "SUM";
        }
        cout << endl;
    }
    
    phyloRead.back_mutate_method = backMutateMethod;
    considerAdjSNPs = toConsiderAdj;
    
    fileLog = prefix + "." + int2str(threadID-1) + ".log";
    // output file for logging sequences
    fileSeq = prefix + "." + int2str(threadID-1) + ".seq";
    // max tree file
    fileMaxTree = prefix + ".max.tree";
    // max seq file
    fileMaxSeq = prefix + ".max.seq";

    // check whether the files of previous run exist
    if (toResume)
        toResume = checkWhetherPreRunFilesExist();
    
    // load the genome file
    cout << "loading the sam/bam file ... " << endl;
    load_sam(samFile, numHap, show_pos);

    // load the starting tree file
    if (startTreeFile != NULL)
        loadStartTreeFile(startTreeFile);

    // resume from the previous run
    else if (toResume)
        resumeFrPreRun();

    // open the output file streams
    flog = new ofstream();
    flog->open(fileLog.c_str(), ofstream::out | ofstream::app);
    // fseq.open(fileSeq.c_str());
    
    cout << "start processMCMC()" << endl << flush;
    processMCMC(nGenerations, nPrintSteps);
    
    // open the output file streams
    flog->close();
    // fseq.close();
}

// get the information of the alignment
void PhyloMCMC::get_align_info(char* samFile, double& avgCover, int& numSNPs, int& seqLen) {
    phyloRead.getAlgnInfo(samFile, avgCover, numSNPs, seqLen);
}

// load samfile
void PhyloMCMC::load_sam(char* samFile, int numHap, set<int>& show_pos) {
    
    // randomly generate a random tree
    // cout << "randomly generate a random tree" << endl << flush;
    phyloRead.computeRandTopology(numHap, myRand);
    phyloRead.isDescendantUpdated = false;
    
    // compute the data matrix
    // cout << "compute the data matrix" << endl << flush;
    phyloRead.loadData(samFile, show_pos);
    
    // set up the variables:
    // cout << "set up the variables" << endl << flush;
    phyloRead.setupVariables();
    
    // initialize the parameters
    // cout << "initialize the parameters" << endl << flush;
    phyloRead.randParameters(myRand);
    
    numFreqs = (int) phyloRead.haploFreqs.size();
    // cout << "numFreqs = " << numFreqs << endl << flush;
}

// load from other object
void PhyloMCMC::load_sam(PhyloMCMC* phyloMCMC, int numHap) {

    // randomly generate a random tree
    phyloRead.computeRandTopology(numHap, myRand);
    phyloRead.isDescendantUpdated = false;
    
    // compute the data matrix
    phyloRead.loadData(phyloMCMC->phyloRead);
    
    // set up the variables:
    phyloRead.setupVariables();
    
    // initialize the parameters
    phyloRead.randParameters(myRand);
    
    numFreqs = (int) phyloRead.haploFreqs.size();
}

// initialize the state
// and get the log-posterior value of the starting state
// and the number of steps in each generation
void PhyloMCMC::init(long double& currLnP, int& numSteps) {
    int i;

    // allocate memory to the arrays
    haploFreqs_curr.clear();
    topologies_curr.clear();
    isDescendant_curr.clear();
    
    // allocate memories to the arrays
    expecLogFreqs_curr = new double[phyloRead.numEdges * phyloRead.numEdges * 16];

    // gamma prior for freqsuencies
    freq_shape = FREQ_SHAPE;
    freq_rate = FREQ_RATE;

    // compute the log-posterior of the current state
    int isPropose = 0;
    logPost(isPropose);
    currLnP = logPost_curr;
 
    // save all the current parameters
    for (i=0; i<numFreqs; i++)
        haploFreqs_curr.push_back(phyloRead.haploFreqs_t[i]);
    for (i=0; i<phyloRead.currTopology.size(); i++)
        topologies_curr.push_back(phyloRead.currTopology.at(i));
    for (i=0; i<phyloRead.isDescendant.size(); i++)
        isDescendant_curr.push_back(phyloRead.isDescendant[i]);

    // the current expected Log frequencies
    memcpy(expecLogFreqs_curr, phyloRead.expecLogFreqs, phyloRead.numEdges * phyloRead.numEdges * 16 * sizeof(double));
    
    // the current sequencing error
    seqErr_curr = phyloRead.seqErr;
    
    // window size for move
    freq_win_size = FREQ_WINSIZE;
    edge_win_size = EDGE_WINSIZE;
    err_win_size = ERR_WINSIZE;
    max_err = MAX_SEQ_ERR;
    // cout << "Window size for haplotype estimation = " << freq_win_size << endl;
    // cout << "edge_win_size = " << edge_win_size << endl;
    // cout << "Window size of sequencing error estimation = " << err_win_size << endl;
    
    if (isColdChain)
        cout << "Maximum sequencing error = " << max_err << endl;

    // add the set of moves
    moveFuncPool.clear();
    moveKPool.clear();
    
    for (i=0; i<phyloRead.numTips; i++)
        addMove(2, i, 1); // update the haplotype freqs
    addMove(3, 0, 1); // update the error rate
    addMove(4, 0, 1); // distribute the frequencies between two nodes
    addMove(5, 0, 1); // swap between frequencies of two nodes
    addMove(6, 0, 1); // NNI
    addMove(7, 0 ,1); // swap between two subtrees
    addMove(8, 0, 2); // merge and break
    addMove(9, 0, 1); // combine moves

    numSteps = moveFuncPool.size();
    if (isColdChain)
        cout << "Number of steps in each generation: " << numSteps << endl;
}

// compute the log-posterior of the current state / propose state
void PhyloMCMC::logPost(int isPropose) {
    int i;
    long double t;
    if (isPropose==1) {
        if (!considerAdjSNPs) {
            logLike_p = phyloRead.logLikeS();
        } else {
            if (METHOD_EQU == 1) {
                // max
                logLike_p = phyloRead.logLike();
            } else {
                // sum
                logLike_p = phyloRead.logLike2();
            }
        }
        logPrior_p = 0.0;
        for (i=0; i<numFreqs; i++) {
            t = logPriorGamma(phyloRead.haploFreqs[i], freq_shape, freq_rate);
            logPrior_p += t;
        }
        logPost_p = logPrior_p + logLike_p;
    } else {
        if (!considerAdjSNPs) {
            logLike_curr = phyloRead.logLikeS();
        } else {
            if (METHOD_EQU == 1) {
                // max
                logLike_curr = phyloRead.logLike();
            } else {
                // sum
                logLike_curr = phyloRead.logLike2();
            }
        }
        logPrior_curr = 0.0;
        for (i=0; i<numFreqs; i++) {
            t = logPriorGamma(phyloRead.haploFreqs[i], freq_shape, freq_rate);
            logPrior_curr += t;
        }
        logPost_curr = logPrior_curr + logLike_curr;
    }
}

// propose a new state
// output the corresponding log-posterior and the hasting ratio
void PhyloMCMC::proposeNewState(int whichStep, long double& logPostValue, long double& logHastRatio) {
    
    int i;
    int isPropose = 1;

    i = myRand.iunif(0, moveKPool.size()-1);
    func_id_curr = moveFuncPool[i];
    k_curr = moveKPool[i];
    switch(func_id_curr) {
        case 0:
            // update the root state
#ifdef SHOW_DETAILS
            cout << "update the root state" << k_curr << endl << flush;
#endif
            updateRootState(k_curr);
            break;
        case 1:
            // update the location of the mutations
#ifdef SHOW_DETAILS
            cout << "update the location of the mutations" << k_curr << endl << flush;
#endif
            updateMutateLocate(k_curr);
            break;
        case 2:
            // update the haplotype frequency
#ifdef SHOW_DETAILS
            cout << "update the haplotype frequency " << k_curr << " ";
            cout << phyloRead.haploFreqs[k_curr];
#endif
            updateHapFreq(k_curr);
#ifdef SHOW_DETAILS
            cout << " -> " << phyloRead.haploFreqs[k_curr] << endl << flush;
#endif
            break;
        case 3:
            // update the error rate
#ifdef SHOW_DETAILS
            cout << "update the error rate ";
            cout << phyloRead.seqErr;
#endif
            updateErrorRate();
#ifdef SHOW_DETAILS
            cout << " -> " << phyloRead.seqErr << endl << flush;
#endif
            break;
        case 4:
            // distribute the frequencies between two nodes
            distributeFreqs(node1, node2);
#ifdef SHOW_DETAILS
            cout << "distribute the frequencies between two nodes" << endl << flush;
            cout << node1 << ":" << haploFreqs_curr[node1] << "->" << phyloRead.haploFreqs[node1] << " ";
            cout << node2 << ":" << haploFreqs_curr[node2] << "->" << phyloRead.haploFreqs[node2] << endl << flush;
#endif
            break;
        case 5:
            // swap frequencies between two nodes
            swapFreqs(node1, node2);
#ifdef SHOW_DETAILS
            cout << "swap frequencies between two nodes" << endl << flush;
            cout << node1 << ":" << haploFreqs_curr[node1] << "->" << phyloRead.haploFreqs[node1] << " ";
            cout << node2 << ":" << haploFreqs_curr[node2] << "->" << phyloRead.haploFreqs[node2] << endl << flush;
#endif
            break;
        case 6:
            // NNI
            i = doNNI();
#ifdef SHOW_DETAILS
            cout << "NNI (" << i << ")" << endl << flush;
#endif
            break;
        case 7:
            // swap subtree
            i = swapSubTree();
#ifdef SHOW_DETAILS
            cout << "swap subtree (" << i << ")" << endl << flush;
#endif
            break;
        case 8:
            // merge and break
            isChanged = mergeAndBreak();
#ifdef SHOW_DETAILS
            cout << "mergeAndBreak (" << isChanged << ")" << endl << flush;
#endif
            break;
        case 9:
            // combine moves
            combineMoves();
#ifdef SHOW_DETAILS
            cout << "combine move" << endl << flush;
#endif
            break;

    }
    logHastRatio = 0.0; // always 0 for this
    logPost(isPropose);
    logPostValue = logPost_p;
    // cout << "logPostValue = " << logPostValue << endl << flush;
}

// check whether the files of previous run exist
int PhyloMCMC::checkWhetherPreRunFilesExist() {
    ifstream fin, fin2;
    // cout << "Checking whether the log files from the previous run exist..." << endl;
    /*
    if (isColdChain) {
        fin.open(fileMaxTree);
        if (!fin.good()) {
            cerr << "Error! The tree file from the previous run does not exist: " << fileMaxTree << endl;
            exit(1);
        }
        fin.close();
    }
    */
    fin2.open(fileLog);
    if (!fin2.good()) {
    	if (isColdChain) {
		cout << "The log file from the previous run does not exist: " << fileLog << endl;
		cout << "Thus the option '-r' is ignored" << endl;
        }
        fin2.close();
        return 0;
    }
    fin2.close();
    return 1;
}

// resume from the previous run
void PhyloMCMC::resumeFrPreRun() {
     // load the max tree and log file

    UnlabelTree utree;

    ifstream fin,fin2;
    string strTree;
    string strErr = "";
    string strLogP = "";
    bool isTreeStr;
    bool maxTreeFileValid;
    vector<int> topInt;
    vector<double> haploFreq;
    vector<double> edgeLen;

    string lastIter;
    string lastTree;
    string currTree;
    string lastErrRate;
    string lastLogP;
    string aline;
    double currLogP;
    vector<string> token;

    int numCol = 0;

    // the prefix of the line with error rate
    string err_prefix = "# error rate: ";
    int err_prefix_len = err_prefix.length();

    // the prefix of the line with posterior probability
    string posterior_prefix = "# max posterior probability: ";
    int posterior_prefix_len = posterior_prefix.length();

    int i,n;
    
    
    if (isColdChain) {
        cout << "loading from the log files of previous run...." << endl;
        // read the max-tree-file
        // cout << "read the max-tree-file" << endl << flush;
        fin.open(fileMaxTree);
        isTreeStr = false;
        maxTreeFileValid = false;

        while ((!isTreeStr) && getline(fin,strTree)) {
            if (strTree.length() == 0)
                continue;
            // get the corresponding error rate
            if (strTree.length() > err_prefix_len && strTree.substr(0,err_prefix_len) == err_prefix) {
                strErr = strTree.substr(err_prefix_len);
                cout << "Max error rate: " << strErr << endl;
                maxErrRate = atof(strErr.c_str());
            }
            // get the corresponding posterior probability
            else if (strTree.length() > posterior_prefix_len && strTree.substr(0,posterior_prefix_len) == posterior_prefix) {
                strLogP = strTree.substr(posterior_prefix_len);
                cout << "Max posterior probability: " << strLogP << endl;
                maxLogPost = atof(strLogP.c_str());
            }
            // read the tree
            else if (strTree[0] != '#') {
                isTreeStr = true;
                maxTree = strTree;
                cout << "Best tree: " << maxTree << endl;
                maxTreeFileValid = true;
            }
        }
        fin.close();
        
	/*
        if (strErr.length() == 0) {
            cerr << "Error! Error rate is not found in " << fileMaxTree << endl;
            exit(1);
        }

        if (strLogP.length() == 0) {
            cerr << "Error! Posterior probability is not found in " << fileMaxTree << endl;
            exit(1);
        }

        if (!isTreeStr) {
            cerr << "Error! The file: " << fileMaxTree << " cannot be read" << endl;
            exit(1);
        }
	*/
    }
    
    fin2.open(fileLog);
    while (getline(fin2,aline)) {
        tokenizer(aline, "\t", &token);
        n = token.size();
        if (n>4 && n>=numCol) {
        	currTree = token[n-1];
            if (currTree.length() > 0 && currTree[currTree.length()-1] == ')' && isBalParentheses(currTree)  && atof(token[1].c_str()) < 0.0) {
				numCol = n;
				lastIter = token[0];
				lastLogP = token[1];
				lastTree = currTree;
				lastErrRate = token[n-2];
				if ((!maxTreeFileValid) && isColdChain) {
					currLogP = atof(lastLogP.c_str());
					if ((!isTreeStr) || (currLogP>maxLogPost)) {
						isTreeStr = true;
						maxTree = lastTree;
						maxErrRate = atof(lastErrRate.c_str());
						maxLogPost = currLogP;
					}
				}
			}
        }
    }
    fin2.close();
    
    if (isColdChain) {
        cout << "last iteration: " << lastIter << endl;
        cout << "last log posterior prob: " << lastLogP << endl;
        cout << "last tree: " << lastTree << endl;
        cout << "last error rate: " << lastErrRate << endl;
        cout << "best log posterior prob: " << maxLogPost << endl;
        cout << "best tree: " << maxTree << endl;
        cout << "best error rate: " << maxErrRate << endl;
    }
    utree.topTxt2Int(lastTree, topInt, haploFreq, edgeLen);
    
    // save all the current parameters
    startIteration = atoi(lastIter.c_str())+1;
    for (i=0; i<numFreqs; i++)
        phyloRead.haploFreqs_t[i] = haploFreq[i];

    phyloRead.currTopology.clear();
    for (i=0; i<topInt.size(); i++)
        phyloRead.currTopology.push_back(topInt[i]);

    phyloRead.seqErr = atof(lastErrRate.c_str());
    
    phyloRead.isDescendantUpdated = false;
    phyloRead.needUpdated = true;
    phyloRead.buildIsDescendant();

}

// accept the new state
void PhyloMCMC::acceptNewState(int whichStep) {
    int i;
    switch(func_id_curr) {
        case 0:
            // update the root state
            rootStates_curr[pos_curr] = phyloRead.rootStates[pos_curr];
            break;
        case 1:
            // update the location of the mutations
            mutationLocs_curr[pos_curr] = phyloRead.mutationLocs[pos_curr];
            break;
        case 2:
            // update the haplotype frequency
            haploFreqs_curr[k_curr] = phyloRead.haploFreqs_t[k_curr];
            break;
        case 3:
            // update the error rate
            seqErr_curr = phyloRead.seqErr;
            break;
        case 4:
            // distribute the frequencies between two nodes
            haploFreqs_curr[node1] = phyloRead.haploFreqs_t[node1];
            haploFreqs_curr[node2] = phyloRead.haploFreqs_t[node2];
            break;
        case 5:
            // swap frequencies between two nodes
            haploFreqs_curr[node1] = phyloRead.haploFreqs_t[node1];
            haploFreqs_curr[node2] = phyloRead.haploFreqs_t[node2];
            break;
        case 6:
        case 7:
            /*
            cout << "size of topologies_curr: " << topologies_curr.size() << endl;
            cout << "size of isDescendant_curr: " << isDescendant_curr.size() << endl;
            cout << "size of phyloRead.currTopology: " << phyloRead.currTopology.size() << endl;
            cout << "size of phyloRead.isDescendant: " << phyloRead.isDescendant.size() << endl;*/
            // NNI / swap subtrees
            for (i=0; i<phyloRead.currTopology.size(); i++)
                topologies_curr[i] = phyloRead.currTopology[i];
            for (i=0; i<phyloRead.isDescendant.size(); i++)
                isDescendant_curr[i] = phyloRead.isDescendant[i];
            break;
        case 8:
            // merge and break
            if (isChanged) {
                for (i=0; i<phyloRead.currTopology.size(); i++)
                    topologies_curr[i] = phyloRead.currTopology[i];
                for (i=0; i<phyloRead.isDescendant.size(); i++)
                    isDescendant_curr[i] = phyloRead.isDescendant[i];
                haploFreqs_curr[node1a] = phyloRead.haploFreqs_t[node1a];
                haploFreqs_curr[node1b] = phyloRead.haploFreqs_t[node1b];
                haploFreqs_curr[node2] = phyloRead.haploFreqs_t[node2];
            }
            break;
        case 9:
            // combine moves
            for (i=0; i<phyloRead.currTopology.size(); i++)
                topologies_curr[i] = phyloRead.currTopology[i];
            for (i=0; i<phyloRead.isDescendant.size(); i++)
                isDescendant_curr[i] = phyloRead.isDescendant[i];
            for (i=0; i<phyloRead.numTips; i++)
                haploFreqs_curr[i] = phyloRead.haploFreqs_t[i];
            break;
            
    }

    // the current expected Log frequencies and which edges have mutations
    if (!considerAdjSNPs)
        memcpy(expecLogFreqs_curr, phyloRead.expecLogFreq1, phyloRead.numEdges * 4 * sizeof(double));
    else
        memcpy(expecLogFreqs_curr, phyloRead.expecLogFreqs, phyloRead.numEdges * phyloRead.numEdges * 16 * sizeof(double));

    logLike_curr = logLike_p;
    logPrior_curr = logPrior_p;
    logPost_curr = logPost_p;
}

// reject the new state
void PhyloMCMC::rejectNewState(int whichStep) {
    int i;
    switch(func_id_curr) {
        case 0:
            // update the root state
            phyloRead.rootStates[pos_curr] = rootStates_curr[pos_curr];
            break;
        case 1:
            // update the location of the mutations
            phyloRead.mutationLocs[pos_curr] = mutationLocs_curr[pos_curr];
            break;
        case 2:
            // update the haplotype frequency
            phyloRead.haploFreqs_t[k_curr] = haploFreqs_curr[k_curr];
            phyloRead.isHapFreqNormalized = false;
            phyloRead.needUpdated = true;
            break;
        case 3:
            // update the error rate
            phyloRead.seqErr = seqErr_curr;
            phyloRead.isLogFreqUpdated = false;
            phyloRead.needUpdated = true;
            break;
        case 4:
            // distribute the frequencies between two nodes
            phyloRead.haploFreqs_t[node1] = haploFreqs_curr[node1];
            phyloRead.haploFreqs_t[node2] = haploFreqs_curr[node2];
            phyloRead.isHapFreqNormalized = false;
            phyloRead.needUpdated = true;
            break;
        case 5:
            // swap frequencies between two nodes
            phyloRead.haploFreqs_t[node1] = haploFreqs_curr[node1];
            phyloRead.haploFreqs_t[node2] = haploFreqs_curr[node2];
            phyloRead.isHapFreqNormalized = false;
            phyloRead.needUpdated = true;
            break;
        case 6:
        case 7:
            // NNI / swap subtrees
            for (i=0; i<phyloRead.currTopology.size(); i++)
                phyloRead.currTopology[i] = topologies_curr[i];
            for (i=0; i<phyloRead.isDescendant.size(); i++)
                phyloRead.isDescendant[i] = isDescendant_curr[i];
            phyloRead.isWeightUpdated = false;
            phyloRead.needUpdated = true;
            break;
        case 8:
            // merge and break
            if (isChanged) {
                for (i=0; i<phyloRead.currTopology.size(); i++)
                    phyloRead.currTopology[i] = topologies_curr[i];
                for (i=0; i<phyloRead.isDescendant.size(); i++)
                    phyloRead.isDescendant[i] = isDescendant_curr[i];
                phyloRead.haploFreqs_t[node1a] = haploFreqs_curr[node1a];
                phyloRead.haploFreqs_t[node1b] = haploFreqs_curr[node1b];
                phyloRead.haploFreqs_t[node2] = haploFreqs_curr[node2];
                phyloRead.isWeightUpdated = false;
                phyloRead.isHapFreqNormalized = false;
                phyloRead.needUpdated = true;
            }
            break;
        case 9:
            // combine moves
            for (i=0; i<phyloRead.currTopology.size(); i++)
                phyloRead.currTopology[i] = topologies_curr[i];
            for (i=0; i<phyloRead.isDescendant.size(); i++)
                phyloRead.isDescendant[i] = isDescendant_curr[i];
            for (i=0; i<phyloRead.numTips; i++)
                phyloRead.haploFreqs_t[i] = haploFreqs_curr[i];
            phyloRead.isWeightUpdated = false;
            phyloRead.isHapFreqNormalized = false;
            phyloRead.needUpdated = true;
            break;
    }

    // the current expected Log frequencies and which edges have mutations
    if (!considerAdjSNPs)
        memcpy(phyloRead.expecLogFreq1, expecLogFreqs_curr, phyloRead.numEdges * 4 * sizeof(double));
    else
        memcpy(phyloRead.expecLogFreqs, expecLogFreqs_curr, phyloRead.numEdges * phyloRead.numEdges * 16 * sizeof(double));
}

// print the result with max posterior probability
void PhyloMCMC::printMaxState() {
    ofstream fout;
    ofstream fout2;

    fout.open(fileMaxTree.c_str());
    fout.precision(NUM_DIGITS);
    fout << "# max posterior probability: " << maxLogPost << endl;
    fout << "# error rate: " << maxErrRate << endl;
    fout << maxTree << endl;
    fout.close();
}

// print the current state
void PhyloMCMC::printState(int round, int needComputeLogLike) {
    int i;
    double sum;
    
    if (needComputeLogLike) {
        phyloRead.needUpdated = true;
        phyloRead.isHapFreqNormalized = false;
        if (!considerAdjSNPs) {
            phyloRead.logLikeS(1);
        } else {
            if (METHOD_EQU == 1) {
                // MAX
                phyloRead.logLike(1);
            } else {
                // SUM
                phyloRead.logLike2();
            }
        }
    }
    
    flog->precision(NUM_DIGITS);
    // cout.precision(NUM_DIGITS);
    
    (*flog) << round << "\t" << logPost_curr << "\t" << logLike_curr << "\t" << logPrior_curr;
    // cout << round << "\t" << logPost_curr << "\t" << logLike_curr << "\t" << logPrior_curr;
    sum = 0.0;
    for (i=0; i<numFreqs; i++) {
        sum += phyloRead.haploFreqs_t[i];
    }
    for (i=0; i<numFreqs; i++) {
        (*flog) << "\t" << phyloRead.haploFreqs_t[i] / sum;
        // cout << "\t" << phyloRead.haploFreqs_t[i] / sum;
    }
    (*flog) << "\t" << phyloRead.seqErr;
    // cout << "\t" << phyloRead.seqErr;
    
    if (METHOD_EQU == 1) {
        // MAX
        (*flog) << "\t" << phyloRead.currTree(NUM_DIGITS);
        // cout << "\t" << phyloRead.currTree(NUM_DIGITS);
    } else {
        // SUM
        (*flog) << "\t" << phyloRead.currTree2(NUM_DIGITS);
        // cout << "\t" << phyloRead.currTree2(NUM_DIGITS);
    }
    
    (*flog) << endl;
    // cout << endl;
    
    // show the sequence information
    // phyloRead.showCurrSeqs(fseq, round);
}

/*
// print the new propose state  (without end of line character)
void PhyloMCMC::printNewProposeState(int round) {
    int i;
    double sum;
    
    std::cout.precision(NUM_DIGITS);
    
    cout << round << "\t" << logPost_p << "\t" << logLike_p << "\t" << logPrior_p;
    sum = 0.0;
    for (i=0; i<numFreqs; i++)
        sum += phyloRead.haploFreqs_t[i];
    for (i=0; i<numFreqs; i++)
        cout << "\t" << phyloRead.haploFreqs_t[i] / sum;
    cout << "\t" << phyloRead.seqErr;
    
    if (METHOD_EQU == 1) {
        // MAX
        cout << "\t" << phyloRead.currTree1(NUM_DIGITS) << endl;
    } else {
        // SUM
        cout << "\t" << phyloRead.currTree2(NUM_DIGITS) << endl;
    }
}
*/

// the log posterior for gamma distribution
long double PhyloMCMC::logPriorGamma(double x, double shape, double rate) {
    // The normalizing constant in the prior densities can be ignored:
    return (shape - 1.0) * logl(x) - rate * x;
}

// the log posterior for beta distribution
long double PhyloMCMC::logPriorBeta(double x, double beta) {
    // The normalizing constant in the prior densities can be ignored:
    return - (beta * x);
}

// do NNI
// return 0 if there is no change in the tree topology
// return 1 if there is change in the tree topology
int PhyloMCMC::doNNI() {
    // NNI
    UnlabelTree uTree;
    bool treeChanged;
    vector<double> fs;

    // cout << "before NNI" << endl;
    // uTree.showTopology(phyloRead.topologies[0]);
    // cout << phyloRead.currTree() << endl;
    treeChanged = uTree.doNNI(phyloRead.currTopology, phyloRead.numTips, myRand, changeNodeIDs);
    if (!treeChanged)
        return 0;
    
    phyloRead.needUpdated = true;
    phyloRead.buildIsDescendant();
    return 1;
}

// swap subtrees
// return 0 if there is no change in the tree topology
// return 1 if there is change in the tree topology
int PhyloMCMC::swapSubTree() {
    UnlabelTree uTree;
    bool treeChanged;
    vector<double> fs;
    
    treeChanged = uTree.swapSubTree(phyloRead.currTopology, phyloRead.numTips, phyloRead.numNodes, myRand, phyloRead.isDescendant, changeNodeIDs);
    if (!treeChanged)
        return 0;
    
    phyloRead.needUpdated = true;
    phyloRead.buildIsDescendant();
    return 1;
}

// distribute the frequencies between two nodes
// hasting ratio is 1
void PhyloMCMC::distributeFreqs(int& node1, int& node2) {
    double newFreq1, newFreq2, sumFreq;
    node1 = myRand.iunif(0, numFreqs-1);
    node2 = myRand.iunif(0, numFreqs-2);
    if (node2 >= node1)
        node2++;
    sumFreq = phyloRead.haploFreqs_t[node1] + phyloRead.haploFreqs_t[node2];
    newFreq1 = myRand.runif(0, sumFreq);
    newFreq2 = sumFreq-newFreq1;
    phyloRead.haploFreqs_t[node1] = newFreq1;
    phyloRead.haploFreqs_t[node2] = newFreq2;
    phyloRead.isHapFreqNormalized = false;
    phyloRead.needUpdated = true;
}

// swap frequencies between two nodes
void PhyloMCMC::swapFreqs(int& node1, int& node2) {
    double newFreq1, newFreq2;
    node1 = myRand.iunif(0, numFreqs-1);
    node2 = myRand.iunif(0, numFreqs-2);
    if (node2 >= node1)
        node2++;
    newFreq1 = phyloRead.haploFreqs_t[node2];
    newFreq2 = phyloRead.haploFreqs_t[node1];
    phyloRead.haploFreqs_t[node1] = newFreq1;
    phyloRead.haploFreqs_t[node2] = newFreq2;
    phyloRead.isHapFreqNormalized = false;
    phyloRead.needUpdated = true;
}

int PhyloMCMC::getIdxFrTop(int node) {
    int i;
    for (i=0; i<phyloRead.currTopology.size(); i++) {
        if (phyloRead.currTopology[i] == node)
            return i;
    }
    return -1;
}

// merge a subtree of 2 tips into a single tip, and break one tip into a subtree of 2 tips
int PhyloMCMC::mergeAndBreak() {
    UnlabelTree utree;
    vector<pair<int,int> > changeNodeIDs;
    vector<int> rowsWithTwoLeaves;
    int i,n,t;
    int internode1, sister_internode1;
    int idx_node2, idx_inter1;
    int node_a, node_b, node_c;
    double new_freq2, new_freq1b;
    // collect all rows with two leaves
    n = phyloRead.currTopology.size();
    for (i=0; i<n; i+=2) {
        if (phyloRead.currTopology[i] >= 0 && phyloRead.currTopology[i+1] >=0
            && phyloRead.currTopology[i] < numFreqs && phyloRead.currTopology[i+1] < numFreqs) {
            rowsWithTwoLeaves.push_back(i/2);
        }
    }
    // randomly pick one subtree with two leaves
    i = myRand.iunif(0, rowsWithTwoLeaves.size()-1);
    internode1 = rowsWithTwoLeaves[i];
    // randomly assign node1a and node1b
    i = myRand.iunif(0,1);
    if (i==0) {
        node1a = phyloRead.currTopology[2*internode1];
        node1b = phyloRead.currTopology[2*internode1+1];
    } else {
        node1a = phyloRead.currTopology[2*internode1+1];
        node1b = phyloRead.currTopology[2*internode1];
    }
    if (node1a < node1b) {
        node_a = node1a;
        node_b = node1b;
    } else {
        node_a = node1b;
        node_b = node1a;
    }
    // get the sister node of internode1
    sister_internode1 = utree.getSisterNode(internode1, phyloRead.currTopology);
    if (sister_internode1 < node_a) {
        node_c = node_b;
        node_b = node_a;
        node_a = sister_internode1;
    } else if (sister_internode1 < node_b) {
        node_c = node_b;
        node_b = sister_internode1;
    } else {
        node_c = sister_internode1;
    }
    // randomly pick one tip not inside the selected subtree
    i = myRand.iunif(0, numFreqs-4);
    if (i>=node_a)
        i++;
    if (i>=node_b)
        i++;
    if (i>=node_c)
        i++;
    node2=i;
    // check whether the condition is satified
    if (phyloRead.haploFreqs_t[node2] - phyloRead.haploFreqs_t[node1a] >= MIN_FREQ_DIFF) {
        // show the topology
        // swapping between the node2 and internal node 1
        // cout << "swapping between node2=" << node2 << " and internal node 1=" << internode1 << endl << flush;
        idx_node2 = getIdxFrTop(node2);
        idx_inter1 = getIdxFrTop(internode1);
        t = phyloRead.currTopology[idx_inter1];
        phyloRead.currTopology[idx_inter1] = phyloRead.currTopology[idx_node2];
        phyloRead.currTopology[idx_node2] = t;
        // update the frequencies
        // node1b's frequency = node2's frequency - node1a's frequency
        // node2's frequency = node1a's frequency + node1b's frequency
        new_freq1b = phyloRead.haploFreqs_t[node2] - phyloRead.haploFreqs_t[node1a];
        new_freq2 = phyloRead.haploFreqs_t[node1a] + phyloRead.haploFreqs_t[node1b];
        phyloRead.haploFreqs_t[node1b] = new_freq1b;
        phyloRead.haploFreqs_t[node2] = new_freq2;
        // tidy up the messy topology
        utree.tidyUpMessyTopology(numFreqs, phyloRead.currTopology, changeNodeIDs);
        phyloRead.isHapFreqNormalized = false;
        phyloRead.needUpdated = true;
        phyloRead.buildIsDescendant();
        return 1;
    }
    return 0;
}

// combine moves
void PhyloMCMC::combineMoves() {
    // randomly select two of actions: 0 - NNI; 1 - swapping subtrees; 2 - merge & break
    // the two actions may be the same
    int action1, action2;
    action1 = myRand.iunif(0, 2);
    action2 = myRand.iunif(0, 2);

    switch(action1) {
        case 0:
#ifdef SHOW_DETAILS
            cout << "1: NNI" << endl << flush;
#endif
            doNNI();
            break;
        case 1:
#ifdef SHOW_DETAILS
            cout << "1: swapSubtree" << endl << flush;
#endif
            swapSubTree();
            break;
        case 2:
#ifdef SHOW_DETAILS
            cout << "1: merge & break" << endl << flush;
#endif
            mergeAndBreak();
            break;
    }
    
    switch(action2) {
        case 0:
#ifdef SHOW_DETAILS
            cout << "2: NNI" << endl << flush;
#endif
            doNNI();
            break;
        case 1:
#ifdef SHOW_DETAILS
            cout << "2: swapSubtree" << endl << flush;
#endif
            swapSubTree();
            break;
        case 2:
#ifdef SHOW_DETAILS
            cout << "2: merge & break" << endl << flush;
#endif
            mergeAndBreak();
            break;
    }
    // phyloRead.normalizeHapFreqs();
    phyloRead.isHapFreqNormalized = false;
    phyloRead.isDescendantUpdated = false;
    phyloRead.needUpdated = true;
}

// update the location of the mutations
void PhyloMCMC::updateMutateLocate(int vpos) {
    int r;
    pos_curr = phyloRead.variablePos[vpos];
    r = myRand.iunif(0,phyloRead.numEdges-2);
    if (r >= phyloRead.mutationLocs[pos_curr]) {
        r++;
    }
    phyloRead.mutationLocs[pos_curr] = r;
}

// update the root state
void PhyloMCMC::updateRootState(int vpos) {
    pos_curr = phyloRead.variablePos[vpos];
    phyloRead.rootStates[pos_curr] = (phyloRead.rootStates[pos_curr]+1) % 2;
    
}

// update the haplotype frequency
void PhyloMCMC::updateHapFreq(int hap) {
    double newValue;
    newValue = phyloRead.haploFreqs_t[hap] + myRand.runif(-freq_win_size/2.0, freq_win_size/2.0);
    if (newValue < 0)
        newValue = -newValue;
    phyloRead.haploFreqs_t[hap] = newValue;
    phyloRead.isHapFreqNormalized = false;
    phyloRead.needUpdated = true;
    phyloRead.normalizeHapFreqs();
}

// update the error rate
void PhyloMCMC::updateErrorRate() {
    double newValue;
    newValue = phyloRead.seqErr + myRand.runif(-err_win_size/2.0, err_win_size/2.0);
    if (newValue < 0)
        newValue = -newValue;
    if (newValue > MAX_SEQ_ERR)
        newValue = 2.0 * MAX_SEQ_ERR - newValue;
    phyloRead.seqErr = newValue;
    phyloRead.isLogFreqUpdated = false;
    phyloRead.needUpdated = true;
}

// add move
void PhyloMCMC::addMove(int func_id, int k, int weight) {
    int i;
    for (i=0; i<weight; i++) {
        moveFuncPool.push_back(func_id);
        moveKPool.push_back(k);
    }
}

// save the best result so far
void PhyloMCMC::saveBestResult() {
    int i;
    if (!maxResultExist || logPost_curr > maxLogPost) {
        // save the result

        phyloRead.needUpdated = true;
        phyloRead.isHapFreqNormalized = false;
        
        if (!considerAdjSNPs) {
            phyloRead.logLikeS(1);
        } else {
        if (METHOD_EQU == 1) {
            // max
            phyloRead.logLike(1);
        } else {
            // sum
            phyloRead.logLike2();
        }
        }
        
        maxResultExist = true;
        maxLogPost = logPost_curr;
        if (METHOD_EQU == 1) {
            // MAX
            maxTree = phyloRead.currTree(NUM_DIGITS);
        } else {
            // SUM
            maxTree = phyloRead.currTree2(NUM_DIGITS);
        }
        maxErrRate = phyloRead.seqErr;
        // phyloRead.getCurrSeqs(maxSeqs);
        maxFreqs.clear();
        for (i=0; i<phyloRead.haploFreqs.size(); i++) {
            maxFreqs.push_back(phyloRead.haploFreqs[i]);
        }
    }
}
