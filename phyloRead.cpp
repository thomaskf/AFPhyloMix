//
//  phyloRead.cpp
//  HaploCount
//
//  Created by Thomas Wong on 15/1/19.
//

#include "phyloRead.h"

// constructor
PhyloRead::PhyloRead() {
    startup();
}

// destructor
PhyloRead::~PhyloRead() {
    if (nucl2int != NULL)
        delete[] nucl2int;
    if (expecLogFreqs != NULL)
        delete[] expecLogFreqs;
    if (expecLogFreq1 != NULL)
        delete[] expecLogFreq1;
    if (rootStates != NULL)
        delete[] rootStates;
    if (mutationLocs != NULL)
        delete[] mutationLocs;
    if (tloglike != NULL)
        delete[] tloglike;
    // if (tloglikeValid != NULL)
    //    delete[] tloglikeValid;
    if (maxMutatePos != NULL)
        delete[] maxMutatePos;
    if (mutatePos != NULL)
        delete[] mutatePos;
    if (maxStateRoot != NULL)
        delete[] maxStateRoot;
    if (stateRoot != NULL)
        delete[] stateRoot;
    if (log_phi != NULL)
        delete[] log_phi;
    if (lpr_ui != NULL)
        delete[] lpr_ui;
    if (charMatHaplo != NULL)
        delete[] charMatHaplo;
    if (isValidCharMat != NULL)
        delete[] isValidCharMat;
    if (lprobRcMc != NULL)
        delete[] lprobRcMc;
    if (logLikeRiMi != NULL)
        delete[] logLikeRiMi;

    if (!useOtherMat) {
        if (dMat != NULL)
            delete[] dMat;
        if (sdMat != NULL)
            delete[] sdMat;
        if (charOrder != NULL)
            delete[] charOrder;
        if (isInVar != NULL)
            delete[] isInVar;
        if (isDiscard != NULL)
            delete[] isDiscard;
    }
}

// compute whether the adjacent columns should be combined
// should the adjacent positions be combined
// vector<int> combinePos;
// if position i can combine with position i+1, then combinePos[i] = 1
// if position i can combine with position i+2, then combinePos[i] = 2
// otherwise, combinePos[i] = 0
void PhyloRead::getCombinePos(Reads& rds, int numPos, vector<int>& combinePos) {
    Read* rd;
    int* distMat;
    int* distMat2;
    int* cov;
    int* cov2;
    int mpos;
    int i,j,c,c1,c2;
    vector<pair<int,int> > grps;
    int canCombine;
    distMat = new int[numPos*16];
    distMat2 = new int[numPos*16];
    cov = new int[numPos];
    cov2 = new int[numPos];
    
    memset(distMat, 0, numPos*16*sizeof(int));
    memset(distMat2, 0, numPos*16*sizeof(int));
    memset(cov, 0, numPos*sizeof(int));
    memset(cov2, 0, numPos*sizeof(int));
    for (i=0; i<(int)rds.readlist.size(); i++) {
        rd = rds.readlist[i];
        mpos = rd->mapPos;
        for (j=1; j<rd->seq.length(); j++) {
            c1 = nucl2int[(int)rd->seq[j-1]];
            c = nucl2int[(int)rd->seq[j]];
            if (c1!=0 && c!=0) {
                c1--; c--;
                distMat[(mpos+j)*16 + c1*4 + c]++;
                cov[mpos+j]++;
            }
            if (j>1) {
                c2 = nucl2int[(int)rd->seq[j-2]];
                c = nucl2int[(int)rd->seq[j]];
                if (c2!=0 && c!=0) {
                    c2--; c--;
                    distMat2[(mpos+j)*16 + c2*4 + c]++;
                    cov2[mpos+j]++;
                }
            }
        }
    }
    // compute the array "combinePos"
    combinePos.clear();
    for (i=0; i<numPos; i++) {
        grps.clear();
        for (j=0; j<16; j++) {
            if ((double)distMat[i*16+j] / (double)cov[i] > 2.0*MAX_SEQ_ERR) {
                grps.push_back(pair<int,int>(j/4,j%4));
            }
        }
        if (grps.size()==1) {
            canCombine = 1;
        } else if (grps.size()==2 && grps[0].first != grps[1].first && grps[0].second != grps[1].second) {
            canCombine = 1;
        } else {
            grps.clear();
            for (j=0; j<16; j++) {
                if ((double)distMat2[i*16+j] / (double)cov2[i] > 2.0*MAX_SEQ_ERR) {
                    grps.push_back(pair<int,int>(j/4,j%4));
                }
            }
            if (grps.size()==1) {
                canCombine = 2;
            } else if (grps.size()==2 && grps[0].first != grps[1].first && grps[0].second != grps[1].second) {
                canCombine = 2;
            } else {
                canCombine = 0;
            }
        }
        combinePos.push_back(canCombine);
    }
    delete[] distMat;
    delete[] distMat2;
    delete[] cov;
    delete[] cov2;
}

// obtain the state frequencies from posMatrix
void PhyloRead::getStateFreqs(int* posMatrix) {
    int i,j;
    double tot;
    logStateFreqs.clear();
    for (i=0; i<4; i++) {
        logStateFreqs.push_back(0.0);
    }
    for (i=0; i<numPos; i++) {
        for (j=1; j<=4; j++) {
            logStateFreqs[j-1] += posMatrix[i*5+j];
        }
    }
    // set the minimum value be 1
    for (i=0; i<4; i++) {
        if (logStateFreqs[i] == 0)
            logStateFreqs[i] = 1;
    }
    tot = logStateFreqs[0];
    for (i=1; i<4; i++) {
        tot += logStateFreqs[i];
    }
    for (i=0; i<4; i++) {
        logStateFreqs[i] = logStateFreqs[i] / tot;
    }
    // show the stateFreqs (before taking log)
    cout << "Freqs of A, C, G, T: ";
    for (i=0; i<4; i++) {
        if (i>0)
            cout << ", ";
        cout << logStateFreqs[i];
    }
    cout << endl;
    for (i=0; i<4; i++) {
        logStateFreqs[i] = log(logStateFreqs[i]);
    }
}

// get alignment information
// Average coverage, number of SNP positions, and sequence length
void PhyloRead::getAlgnInfo(char* samFile, double& avgCover, int& numSNPs, int& seqLen) {
    
    int* posMatrix;
    vector<int> combinePos;
    vector<int> coverage;
    vector<int> coverage_w_gap;
    vector<int> snpPosSet;
    vector<int> winRefPos;
    int tot;
    int i,j;
    double coverage_avg = 0.0;
    
    rds.readSamFile(samFile);
    posMatrix = rds.getPosMatrix(numPos);
    getStateFreqs(posMatrix);
    getCombinePos(rds, numPos, combinePos);
    
    // compute the coverages
    for (i=0; i<numPos; i++) {
        tot = 0;
        for (j=1; j<=4; j++) {
            tot += posMatrix[i*5+j];
        }
        coverage.push_back(tot);
        coverage_avg += tot;
        // cout << i << "\t" << tot;
        tot += posMatrix[i*5];
        coverage_w_gap.push_back(tot);
        // cout << "\t" << tot << endl;
    }
    avgCover = coverage_avg / (double) numPos;
    
    // set up the arrays
    charOrder = new int[numPos * 4];
    isInVar = new bool[numPos]; // whether the site is invariable
    isDiscard = new bool[numPos]; // whether the site is discarded
    
    // construct the arrays chars0, chars1, isInvar, isDiscard
    getDataInfo(posMatrix, coverage, coverage_w_gap);
    
    // collect a set of potential snp positions
    // i.e. the positions which have enough coverage and cannot be merged with the other positions
    for (i=0; i<numPos; i++) {
        if (coverage[i] >= COVER_THRES && combinePos[i]==0 && (!isDiscard[i]) && (!isInVar[i])) {
            snpPosSet.push_back(i);
        }
    }
    
    numSNPs = (int) snpPosSet.size();
    seqLen = numPos;

    // clear the memory
    delete[] posMatrix;
}

// (for debugging) to get the actual back mutation snp positions
void PhyloRead::getActualBackMutate(char* samFile, vector<vector<int> >& c_sets, int numHaps) {
    vector<int> combinePos;
    vector<int> coverage;
    vector<int> coverage_w_gap;
    vector<int> snpPosSet;
    vector<pair<int,int> > snpChars;
    vector<int> winRefPos;
    int tot;
    int i,j;
    int* posMatrix;
    double coverage_avg = 0.0;
    
    rds.readSamFile(samFile);
    posMatrix = rds.getPosMatrix(numPos);
    getStateFreqs(posMatrix);
    getCombinePos(rds, numPos, combinePos);
    
    // compute the coverages
    for (i=0; i<numPos; i++) {
        tot = 0;
        for (j=1; j<=4; j++) {
            tot += posMatrix[i*5+j];
        }
        coverage.push_back(tot);
        coverage_avg += tot;
        // cout << i << "\t" << tot;
        tot += posMatrix[i*5];
        coverage_w_gap.push_back(tot);
        // cout << "\t" << tot << endl;
    }
    coverage_avg = coverage_avg / (double) numPos;
    cout << "Average coverage: " << coverage_avg << endl;
    
    // set up the arrays
    charOrder = new int[numPos * 4];
    isInVar = new bool[numPos]; // whether the site is invariable
    isDiscard = new bool[numPos]; // whether the site is discarded
    
    // construct the arrays chars0, chars1, isInvar, isDiscard
    getDataInfo(posMatrix, coverage, coverage_w_gap);
    
    // collect a set of potential snp positions
    // i.e. the positions which have enough coverage and cannot be merged with the other positions
    for (i=0; i<numPos; i++) {
        // if (coverage[i] >= COVER_THRES && combinePos[i]==0 && (!isDiscard[i]) && (!isInVar[i])) {
        if (coverage[i] >= COVER_THRES && coverage[i] >= (COVER_RATIO_THRES * coverage_avg) && (!isDiscard[i]) && (!isInVar[i])) {
            snpPosSet.push_back(i);
            snpChars.push_back(pair<int,int>(charOrder[i*4],charOrder[i*4+1]));
        }
    }
    
    // cout << "Number of potential snp positions: " << snpPosSet.size() << endl;
    // cout << "Length of sequences: " << numPos << endl;
    
    if (snpPosSet.size() == 0) {
        cout << "Error! There is no snp position!" << endl;
        exit(1);
    }
    
    // (for debugging) check the snp positions
    checkSNPPos(snpPosSet, c_sets, numHaps);
}

void PhyloRead::loadData(char* samFile, set<int>& show_pos) {
    
    Read* rd;
    vector<int> combinePos;
    vector<int> coverage;
    vector<int> coverage_w_gap;
    vector<int> snpPosSet;
    vector<pair<int,int> > snpChars;
    vector<int> winRefPos;
    vector<pair<int,int> > dropRegions;
    int tot;
    int c,c2;
    int i,j,k,p,refpos,currpos,refj,idx,pre_j;
    int mpos,endpos;
    int* dMatrix;
    int* posMatrix;
    int t1[4], t2[4];
    int n[64];
    double coverage_avg = 0.0;
    
    rds.readSamFile(samFile);
    posMatrix = rds.getPosMatrix(numPos);
    getStateFreqs(posMatrix);
    getCombinePos(rds, numPos, combinePos);
    
    // compute the coverages
    for (i=0; i<numPos; i++) {
        tot = 0;
        for (j=1; j<=4; j++) {
            tot += posMatrix[i*5+j];
        }
        coverage.push_back(tot);
        coverage_avg += tot;
        // cout << i << "\t" << tot;
        tot += posMatrix[i*5];
        coverage_w_gap.push_back(tot);
        // cout << "\t" << tot << endl;
    }
    coverage_avg = coverage_avg / (double) numPos;
    cout << "Average coverage: " << coverage_avg << endl;
    
    // set up the arrays
    charOrder = new int[numPos * 4];
    isInVar = new bool[numPos]; // whether the site is invariable
    isDiscard = new bool[numPos]; // whether the site is discarded
    
    // construct the arrays chars0, chars1, isInvar, isDiscard
    getDataInfo(posMatrix, coverage, coverage_w_gap);
    
    // get the dropping regions
    rds.getCoverDropRegs(dropRegions, numPos, posMatrix);
    for (i=0; i<dropRegions.size(); i++) {
        for (j=dropRegions[i].first; j<=dropRegions[i].second; j++) {
            isDiscard[j] = true;
        }
    }
    
    // collect a set of potential snp positions
    // i.e. the positions which have enough coverage and cannot be merged with the other positions
    for (i=0; i<numPos; i++) {
        // if (coverage[i] >= COVER_THRES && combinePos[i]==0 && (!isDiscard[i]) && (!isInVar[i])) {
        if (coverage[i] >= COVER_THRES && coverage[i] >= (COVER_RATIO_THRES * coverage_avg) && (!isDiscard[i]) && (!isInVar[i])) {
            snpPosSet.push_back(i);
            snpChars.push_back(pair<int,int>(charOrder[i*4],charOrder[i*4+1]));
            // cout << i << "\t" << charOrder[i*4] << "\t" << charOrder[i*4+1] << endl;
        }
    }
    cout << "Number of potential snp positions: " << snpPosSet.size() << endl;
    cout << "Length of sequences: " << numPos << endl;
    
    if (snpPosSet.size() == 0) {
        cout << "Error! There is no snp position!" << endl;
        exit(1);
    }
    
    // identify the bad mutation
    cout << "Identifying the bad mutation..." << endl;
    findBadMutations(snpPosSet, snpChars, show_pos);
    
    /*
    // show the snp positions
    cout << (int) snpPosSet.size() << " SNP positions to consider" << endl;
    for (i=0; i<(int) snpPosSet.size(); i++) {
        cout << snpPosSet[i] << endl;
    }
    */
    // exit(1);
    
    // select the reference position for each window
    // consider the set of potential snp positions
    
    // get the first snp position
    int preRefPos = snpPosSet[0];
    winRefPos.push_back(preRefPos);
    i=1;
    
    while (i<snpPosSet.size()) {
        if (snpPosSet[i] - preRefPos > 2 * WINDOW_SIZE) {
            preRefPos = snpPosSet[i];
            winRefPos.push_back(preRefPos);
        }
        i++;
    }
    
    // get the last snp position
    if (winRefPos[winRefPos.size()-1] != snpPosSet[snpPosSet.size()-1]) {
        winRefPos.push_back(snpPosSet[snpPosSet.size()-1]);
    }
    
    // for each position, get the nearest reference position before and after
    vector<int> refBeforePos;
    vector<int> refAfterPos;
    pre_j = -1;
    for (i=0; i<winRefPos.size(); i++) {
        j = winRefPos[i];
        for (k=pre_j+1; k<j; k++) {
            refBeforePos.push_back(pre_j);
            refAfterPos.push_back(j);
        }
        refBeforePos.push_back(j);
        refAfterPos.push_back(j);
        pre_j = j;
    }
    for (k=pre_j+1; k<numPos; k++) {
        refBeforePos.push_back(pre_j);
        refAfterPos.push_back(-1);
    }
    
    // get the total coverage of reads covering the region from refBeforePos[i] to i
    // and the total coverage of reads covering the region from i to refAfterPos[i]
    vector<int> covBeforeRef (numPos, 0);
    vector<int> covAfterRef (numPos, 0);
    for (i=0; i<(int)rds.readlist.size(); i++) {
        rd = rds.readlist[i];
        mpos = rd->mapPos;
        endpos = mpos + rd->seq.length() - 1;
        for (j=0; j<rd->seq.length(); j++) {
            if (rd->seq.at(j) == 'N')
                continue;
            currpos = j+mpos;
            // check the reference position before
            refpos = refBeforePos[currpos];
            if (refpos >= mpos && refpos <= endpos) {
                covBeforeRef[currpos]++;
            }
            if (refpos != currpos) {
                // check the reference position after
                refpos = refAfterPos[currpos];
                if (refpos >= mpos && refpos <= endpos) {
                    covAfterRef[currpos]++;
                }
            }
        }
    }
    
    // set the refPos array
    refPos.clear();
    pre_j = -1;
    for (i=0; i<numPos; i++) {
        if (covBeforeRef[i] > covAfterRef[i]) {
            j = refBeforePos[i];
        } else {
            j = refAfterPos[i];
        }
        if (j < pre_j)
            j = pre_j;
        refPos.push_back(j);
        pre_j = j;
    }
    
    // set the refLst array
    refLst.clear();
    pre_j = -1;
    for (i=0; i<numPos; i++) {
        if (refPos[i] != pre_j) {
            pre_j = refPos[i];
            refLst.push_back(pre_j);
        }
    }
    
    // get the dMatrix
    dMatrix = new int[16 * numPos];
    memset(dMatrix, 0, 16 * numPos * sizeof(int));
    for (i=0; i<(int)rds.readlist.size(); i++) {
        rd = rds.readlist[i];
        mpos = rd->mapPos;
        endpos = mpos + rd->seq.length() - 1;
        for (j=0; j<rd->seq.length(); j++) {
            currpos = j+mpos;
            refpos = refPos[currpos];
            if (refpos >= mpos && refpos <= endpos) {
                refj = refpos - mpos;
                c = nucl2int[(int)rd->seq[refj]];
                c2 = nucl2int[(int)rd->seq[j]];
                if (c==0 || c2==0)
                    continue;
                c--; c2--;
                dMatrix[currpos*16 + c*4 + c2]++;
            }
        }
    }
    
    // compute the logCombins
    computelogCombin1(posMatrix);
    computelogCombin(dMatrix);

    // construct the sdMat
    sdMat = new int[8*numPos];
    
    for (i=0; i<4; i++)
        n[i] = i;
    n[4]=1;
    n[5]=0;
    n[6]=2;
    n[7]=3;
    p=0;
    for (i=0; i<numPos; i++) {
        for (j=0; j<8; j++) {
            c = charOrder[i*4 + n[j]];
            sdMat[p++] = posMatrix[i*5 + c + 1];
        }
    }

    // construct the dMat
    dMat = new int[64*numPos];
    
    //00
    t1[0] = 0; t1[1] = 1; t1[2] = 2; t1[3] = 3;
    t2[0] = 0; t2[1] = 1; t2[2] = 2; t2[3] = 3;
    k=0;
    for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
            n[k] = t1[i]*4 + t2[j];
            k++;
        }
    }
    //01
    t1[0] = 0; t1[1] = 1; t1[2] = 2; t1[3] = 3;
    t2[0] = 1; t2[1] = 0; t2[2] = 2; t2[3] = 3;
    k=0;
    for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
            n[16+k] = t1[i]*4 + t2[j];
            k++;
        }
    }
    //10
    t1[0] = 1; t1[1] = 0; t1[2] = 2; t1[3] = 3;
    t2[0] = 0; t2[1] = 1; t2[2] = 2; t2[3] = 3;
    k=0;
    for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
            n[32+k] = t1[i]*4 + t2[j];
            k++;
        }
    }
    //11
    t1[0] = 1; t1[1] = 0; t1[2] = 2; t1[3] = 3;
    t2[0] = 1; t2[1] = 0; t2[2] = 2; t2[3] = 3;
    k=0;
    for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
            n[48+k] = t1[i]*4 + t2[j];
            k++;
        }
    }
    
    // for each position:
    // int[16] for 00
    // int[16] for 01
    // int[16] for 10
    // int[16] for 11
    p=0;
    for (i=0; i<numPos; i++) {
        refj = refPos[i];
        if (i == refj) {
            for (j=0; j<4; j++) {
                // 0 1 2 3
                c = charOrder[i*4 + j];
                dMat[p++] = dMatrix[i*16+c*5];
            }
            // 1
            c = charOrder[i*4 + 1];
            dMat[p++] = dMatrix[i*16+c*5];
            // 0
            c = charOrder[i*4];
            dMat[p++] = dMatrix[i*16+c*5];
            for (j=2; j<4; j++) {
                // 2 3
                c = charOrder[i*4 + j];
                dMat[p++] = dMatrix[i*16+c*5];
            }
            p+=56;
        } else {
            for (j=0; j<4; j++) {
                for (k=0; k<16; k++) {
                    idx = n[16*j+k];
                    c = charOrder[refj*4 + idx/4];
                    c2 = charOrder[i*4 + idx%4];
                    dMat[p] = dMatrix[i*16 + c*4 + c2];
                    p++;
                }
            }
        }
    }
    
    // clear the memory
    delete[] dMatrix;
    delete[] posMatrix;
}

// from another copy of phyloRead
// the data matrix will be pointing to this object
void PhyloRead::loadData(PhyloRead& aPhyloRead) {

    numPos = aPhyloRead.numPos;
    useOtherMat = true;
    dMat = aPhyloRead.dMat;
    sdMat = aPhyloRead.sdMat;
    charOrder = aPhyloRead.charOrder;
    isInVar = aPhyloRead.isInVar;
    isDiscard = aPhyloRead.isDiscard;
    logCombin1 = aPhyloRead.logCombin1;
    logCombin = aPhyloRead.logCombin;
    // refPos and refLst
    refPos.clear();
    refLst.clear();
    logStateFreqs.clear();
    refPos.insert(refPos.begin(),aPhyloRead.refPos.begin(),aPhyloRead.refPos.end());
    refLst.insert(refLst.begin(),aPhyloRead.refLst.begin(),aPhyloRead.refLst.end());
    logStateFreqs.insert(logStateFreqs.begin(),aPhyloRead.logStateFreqs.begin(),aPhyloRead.logStateFreqs.end());
}

bool cmpPairInt (pair<int,int> p1,pair<int,int> p2) { return (p1.first > p2.first); }

// get information of the data
// construct the arrays chars0, chars1, isInVar, isDiscard
void PhyloRead::getDataInfo(int* posMatrix, vector<int>& coverage, vector<int>& coverage_with_gap) {
    vector<pair<int,int> > nuclpairs;
    int i,j,k;
    variablePos.clear();
    k=0;
    for (i=0; i<numPos; i++) {
        if (coverage[i] == 0 || ((double)posMatrix[i*5] / coverage_with_gap[i] > MAX_GAP_RATE)) {
            // too many reads with gap at this site
            isDiscard[i] = true;
            for (j=0; j<4; j++)
                charOrder[i*4+j] = 0;
            isInVar[i] = false;
        } else {
            nuclpairs.clear();
            for (j=1; j<=4; j++) {
                nuclpairs.push_back(pair<int,int>(posMatrix[i*5+j],j-1));
            }
            sort(nuclpairs.begin(),nuclpairs.end(), cmpPairInt);
            for (j=0; j<4; j++)
                charOrder[i*4+j] = nuclpairs[j].second;
            if ((double)nuclpairs[2].first / (double)coverage[i] > MAX_SEQ_ERR) {
                // more than one changes on the position
                isDiscard[i] = true;
                isInVar[i] = false;
                k++;
            } else if ((double)nuclpairs[1].first / (double)coverage[i] > MAX_SEQ_ERR) {
                // a snp position
                isDiscard[i] = false;
                isInVar[i] = false;
                variablePos.push_back(i);
            } else {
                // an invariable site
                isDiscard[i] = false;
                isInVar[i] = true;
            }
        }
    }
    numVarPos = variablePos.size();
    cout << k << " positions with more than one (but not back) changes" << endl;
}

void backmutate0(vector<char>& pairstate, int n, vector<int>& badPos, vector<double>& badPairRatio, vector<int>& numPair) {
    
    vector<char> isSelected (n, 0);
    int i,j;
    int numBad, numGood;
    int state;
    int max_i, max_k;
    double r;
    double max_ratio;
    badPos.clear();
    badPairRatio.clear();
    numPair.clear();
    do {
        max_ratio = 0.0;
        for (i=0; i<n; i++) {
            if (isSelected[i])
                continue;
            numBad = numGood = 0;
            for (j=0; j<n; j++) {
                if (i==j || isSelected[j])
                    continue;
                if (i<j)
                    state = pairstate[i*n + j];
                    else
                        state = pairstate[j*n + i];
                        if (state == 1) // good
                            numGood++;
                        else if (state == 2) // bad
                            numBad++;
            }
            if (numBad < BAD_MIN_PAIR)
                continue;
            r = (double) numBad / (numBad + numGood);
            if (r > max_ratio) {
                max_i = i;
                max_ratio = r;
                max_k = numBad + numGood;
            }
        }
        if (max_ratio > badPairRatioTHRES) {
            isSelected[max_i] = 1;
            badPos.push_back(max_i);
            badPairRatio.push_back(max_ratio);
            numPair.push_back(max_k);
        }
    } while (max_ratio > badPairRatioTHRES);
}

void backmutate1(vector<char>& pairstate, int n, vector<int>& badPos, vector<double>& badPairRatio, vector<int>& numPair, int numBadThres) {
    
    // report the positions with number of bad >= numBadThres
    
    int i,j;
    int numBad, numGood;
    int state;
    badPos.clear();
    badPairRatio.clear();
    numPair.clear();
    
    for (i=0; i<n; i++) {
        numBad = 0;
        numGood = 0;
        for (j=0; j<n; j++) {
            if (i==j)
                continue;
            if (i<j)
                state = pairstate[i*n + j];
            else
                state = pairstate[j*n + i];
            if (state == 1) // good
                numGood++;
            else if (state == 2) // bad
                numBad++;
        }
        if (numBad >= numBadThres) {
            badPos.push_back(i);
            badPairRatio.push_back((double)numBad / (numBad + numGood));
            numPair.push_back(numBad + numGood);
        }
    }
}


// get the bad positions
// pairstate: nxn dimension; 1 - good; 2 - bad; 0 - not enough coverage
void PhyloRead::selectBadPos(vector<char>& pairstate, int n, vector<int>& badPos, vector<double>& badPairRatio, vector<int>& numPair) {
    
    if (back_mutate_method==0) {
        backmutate0(pairstate, n, badPos, badPairRatio, numPair);
    } else {
        // report the positions with number of bad mutations >= back_mutate_method
        backmutate1(pairstate, n, badPos, badPairRatio, numPair, back_mutate_method);
    }
}

int whichHaplo(string& desc) {
    int p1;
    p1 = desc.find_first_of("-_");
    if (p1!=string::npos) {
        return desc[p1+1]-'A';
    }
    return -1;
}

bool hasMoreThanOneMutation(vector<int>& char_set, vector<vector<int> >& c_sets) {
    int cc_num; // total number of distinct characters
    int num_1;
    int c,c1,c2;
    int i,j;
    vector<int> cc_set;
    bool match;
    c1=c2=-1;
    cc_num=0;
    num_1 = 0;
    if (char_set.size() == 0)
        return false;
    c1 = char_set[0];
    cc_num++;
    cc_set.push_back(0);
    for (i=1; i<char_set.size(); i++) {
        c = char_set[i];
        if (c==c1) {
            cc_set.push_back(0);
        } else if (cc_num==1) {
            cc_num++;
            c2=c;
            cc_set.push_back(1);
            num_1++;
        } else if (c==c2) {
            cc_set.push_back(1);
            num_1++;
        } else {
            // more than one mutation
            return true;
        }
    }
    if (num_1 == 1 || num_1 == char_set.size()-1) {
        return false;
    }
    
    for (i=0; i<c_sets.size(); i++) {
        match = true;
        for (j=0; j<c_sets[i].size(); j++) {
            if (cc_set[j] != c_sets[i].at(j))
                match = false;
        }
        if (match)
            return false;
    }
    return true;
}

int diffCombin(vector<vector<int> >& char_sets, int i, int j) {
    set<pair<int,int> > p;
    int n,k;
    p.clear();
    n = char_sets[i].size();
    for (k=0; k<n; k++) {
        p.insert(pair<int,int>(char_sets[i].at(k),char_sets[j].at(k)));
    }
    return p.size();
}

// (for debugging) check the snp positions
void PhyloRead::checkSNPPos(vector<int>& snpPosSet, vector<vector<int> >& c_sets, int numhap) {

    int n = snpPosSet.size(); // number of snp positions
    vector<int> whichSnp (numPos, -1);
    
    // number of reads for each char for each haplotype for each snp position
    // haploCharSnpPos[i * numhap * 5 + j * 5 + k] = number of reads for haplotype j with character k at snp position i
    vector<int> haploCharSnpPos (numhap * 5 * n,0);
    
    Read *rd1, *rd2;
    int mpos1, mpos2, endpos1, endpos2;
    int currhap;
    int maxChar, maxCharCov;
    int i,j,k,c;
    
    vector<int> char_set;

    for (i=0; i<n; i++) {
        whichSnp[snpPosSet[i]] = i;
    }

    for (i=0; i<(int)rds.pairReads.size(); i++) {
        rd1 = rds.pairReads[i].first;
        currhap = whichHaplo(rd1->desc);
        if (currhap==-1) {
            cerr << rd1->desc << "\t" << currhap << endl;
            exit(1);
            continue;
        }
        rd2 = rds.pairReads[i].second;
        mpos1 = rd1->mapPos;
        mpos2 = rd2->mapPos;
        endpos1 = rd1->mapEndPos;
        endpos2 = rd2->mapEndPos;
        for (j=0; j<rd1->seq.length(); j++) {
            k = whichSnp[mpos1 + j];
            if (k > -1) {
                c = nucl2int[(int)rd1->seq[j]];
                haploCharSnpPos[k * numhap * 5 + currhap * 5 + c]++;
            }
        }
        for (j=0; j<rd2->seq.length(); j++) {
            k = whichSnp[mpos2 + j];
            if (k > -1) {
                c = nucl2int[(int)rd2->seq[j]];
                haploCharSnpPos[k * numhap * 5 + currhap * 5 + c]++;
            }
        }
    }
    
    /*
    // for checking
    vector<vector<int> > char_sets;
    for (i=0; i<snpPosSet.size(); i++) {
        char_set.clear();
        for (j=0; j<numhap; j++) {
            maxChar = 0;
            maxCharCov = haploCharSnpPos[i * numhap * 5 + j * 5];
            for (k=1; k<5; k++) {
                if (haploCharSnpPos[i * numhap * 5 + j * 5 + k] > maxCharCov) {
                    maxChar = k;
                    maxCharCov = haploCharSnpPos[i * numhap * 5 + j * 5 + k];
                }
            }
            char_set.push_back(maxChar);
        }
        char_sets.push_back(char_set);
    }
    int max_dist = 1000;
    bool found;
    for (i=0; i<snpPosSet.size(); i++) {
        if (snpPosSet[i] != 1373)
            continue;
        found = false;
        for (j=i-1; j>=0; j--) {
            if (snpPosSet[i] - snpPosSet[j] > max_dist)
                break;
            if (diffCombin(char_sets, j, i) > 3) {
                found = true;
                cout << snpPosSet[i] << "\t" << snpPosSet[j] << endl;
                // break;
            }
        }
        // if (found)
        //    continue;
        for (j=i+1; j<snpPosSet.size(); j++) {
            if (snpPosSet[j] - snpPosSet[i] > max_dist)
                break;
            if (diffCombin(char_sets, i, j) > 3) {
                found = true;
                cout << snpPosSet[i] << "\t" << snpPosSet[j] << endl;
                // break;
            }
        }
    }
     */

    // check whether the position has more than one mutation
    cout << "Actual positions with back mutation:" << endl;
    for (i=0; i<snpPosSet.size(); i++) {
        char_set.clear();
        for (j=0; j<numhap; j++) {
            maxChar = 0;
            maxCharCov = haploCharSnpPos[i * numhap * 5 + j * 5];
            for (k=1; k<5; k++) {
                if (haploCharSnpPos[i * numhap * 5 + j * 5 + k] > maxCharCov) {
                    maxChar = k;
                    maxCharCov = haploCharSnpPos[i * numhap * 5 + j * 5 + k];
                }
            }
            char_set.push_back(maxChar);
        }
        if (hasMoreThanOneMutation(char_set, c_sets)) {
            cout << snpPosSet[i];
            for (j=0; j<numhap; j++) {
                cout << "\t" << char_set[j];
            }
            cout << endl;
        }
    }
    
    /*
    // list out the statistics
    for (i=0; i<snpPosSet.size(); i++) {
        cout << snpPosSet[i];
        for (j=0; j<numhap; j++) {
            for (k=0; k<5; k++) {
                cout << "\t" << haploCharSnpPos[i * numhap * 5 + j * 5 + k];
            }
        }
        cout << endl;
    }
    */
}


// identify the bad mutations
void PhyloRead::findBadMutations(vector<int>& snpPosSet, vector<pair<int,int> >& snpChars, set<int>& show_pos) {

    long n = (long) snpPosSet.size();
    if (n < 2)
        return;

    vector<int> whichSnp (numPos, -1);
    vector<char> pairstate (n*n, 0); // 1 - good; 2 - bad; 0 - not enough coverage
    map<int,int> readSnp;
    map<int,int>::iterator itr, itr2;
    vector<int> badPos;
    vector<double> badPairRatio;
    vector<int> numPair;
    vector<int> tmp;
    int* nuclPair = new int[(long) n * n * 16];
    int* coverage = new int[(long) n * n];
    Read *rd1, *rd2;
    int mpos1, endpos1;
    int mpos2, endpos2;
    long i,j,k,ii,jj,jc,p,pc;
    int snp_id, snp_c;
    double r;
    double avg_pair_coverage;
    int pair_n;

    for (i=0; i<n; i++) {
        whichSnp[snpPosSet[i]] = i;
    }
    
    memset(nuclPair, 0, (long) n * n * 16 * sizeof(int));
    memset(coverage, 0, (long) n * n * sizeof(int));
    
    for (i=0; i<(int)rds.pairReads.size(); i++) {
        rd1 = rds.pairReads[i].first;
        rd2 = rds.pairReads[i].second;
        mpos1 = rd1->mapPos;
        mpos2 = rd2->mapPos;
        endpos1 = rd1->mapEndPos;
        endpos2 = rd2->mapEndPos;
        readSnp.clear();
        for (j=0; j<rd1->seq.length(); j++) {
            snp_id = whichSnp[mpos1 + j];
            if (snp_id > -1) {
                snp_c = nucl2int[(int)rd1->seq[j]]-1;
                if (snp_c == snpChars[snp_id].first || snp_c == snpChars[snp_id].second) {
                    readSnp.insert(pair<int,int>(snp_id, snp_c));
                }
            }
        }
        for (j=0; j<rd2->seq.length(); j++) {
            snp_id = whichSnp[mpos2 + j];
            if (snp_id > -1) {
                snp_c = nucl2int[(int)rd2->seq[j]]-1;
                if (snp_c == snpChars[snp_id].first || snp_c == snpChars[snp_id].second) {
                    itr = readSnp.find(snp_id);
                    if (itr == readSnp.end()) {
                        // not found
                        readSnp.insert(pair<int,int>(snp_id, snp_c));
                    } else {
                        if (itr->second != snp_c) {
                            itr->second = -1;
                        }
                    }
                }
            }
        }
        
        for (itr=readSnp.begin(); itr!=readSnp.end(); itr++) {
            if (itr->second == -1)
                continue;
            jj = itr->first * n;
            jc = itr->second * 4;
            for (itr2=itr; itr2!=readSnp.end(); itr2++) {
                if (itr2 == itr)
                    continue;
                if (itr2->second == -1)
                    continue;
                p = jj + itr2->first;
                pc = jc + itr2->second;
                nuclPair[(long) p * 16 + pc]++;
                coverage[p]++;
            }
        }
    }
    
    /*
    // collect the information for the pair of positions 1181 and 1275
    cout << "Information for the pair of the positions: 1181 and 1275" << endl;
    for (i=0; i<n-1; i++) {
        if (snpPosSet[i] == 1181) {
            for (j=i+1; j<n; j++) {
                p = i*n + j;
                if (snpPosSet[j] == 1275) {
                    for (k=0; k<16; k++) {
                        if (nuclPair[p * 16 + k] > 0) {
                            cout << nuclPair[p * 16 + k] << endl;
                        }
                    }
                }
            }
        }
    }
    exit(1);
    */
    
    // compute the average coverage for a pair of positions
    avg_pair_coverage = 0.0;
    pair_n = 0;
    for (i=0; i<n-1; i++) {
        ii = i*n;
        for (j=i+1; j<n; j++) {
            p = ii+j;
            if (coverage[p] >= COVER_THRES) {
                avg_pair_coverage += coverage[p];
                pair_n++;
            }
        }
    }
    avg_pair_coverage = avg_pair_coverage / (double) pair_n;
    cout << "Average coverage of a pair of positions: " << avg_pair_coverage << endl;
    double min_r;
    bool first;
    int i_pos, j_pos;
#ifdef SHOWPOS
    bool i_pos_show, j_pos_show;
#endif
    for (i=0; i<n-1; i++) {
        ii = i*n;
        i_pos = snpPosSet[i];
#ifdef SHOWPOS
        i_pos_show = (show_pos.find(i_pos)!=show_pos.end());
#endif
        for (j=i+1; j<n; j++) {
            p = ii+j;
            j_pos = snpPosSet[j];
#ifdef SHOWPOS
            j_pos_show = (show_pos.find(j_pos)!=show_pos.end());
#endif
            if (coverage[p] >= COVER_THRES && coverage[p] >= avg_pair_coverage * PAIR_COVER_RATIO_THRES) {
                jj=0;
                first = true;
                for (k=0; k<16; k++) {
                    r = (double) nuclPair[(long) p * 16 + k] / coverage[p];
                    if (nuclPair[(long) p * 16 + k] > 0 && first) {
                        min_r = r;
                        first = false;
                    } else if (nuclPair[(long) p * 16 + k] > 0 && min_r > r) {
                        min_r = r;
                    }
#ifdef SHOWPOS
                    if ((i_pos_show || j_pos_show) && (nuclPair[p * 16 + k] > 0)) {
                        cout << i_pos << "\t" << j_pos << "\t" << nuclPair[p * 16 + k] << "\t" << coverage[p] << "\t" << r << endl;
                    }
#endif
                    if (r >= BAD_MIN_RATIO) {
                        jj++;
                    }
                }
                if (jj>3) {
                    pairstate[ii+j] = 2; // bad
#ifdef SHOWPOS
                    if (i_pos_show || j_pos_show)
                        cout << i_pos << "\t" << j_pos << "\t" << "bad" << "\tmin_r" << min_r << endl;
#endif
                } else {
                    pairstate[ii+j] = 1; // good
                }
            }
        }
    }
#ifdef SHOWPOS
    exit(1);
#endif

#ifdef LISTBADPOS
    cout << "Number of SNP positions before bad position detection: " << snpPosSet.size() << endl;
    /*
    for (i=0; i<snpPosSet.size(); i++) {
        cout << snpPosSet[i] << endl;
    }*/
#endif

    // get the bad positions
    selectBadPos(pairstate, n, badPos, badPairRatio, numPair);
    
    cout << "Number of bad positions detected: " << badPos.size() << endl;
    
    // discard those bad positions
    for (i=0; i<badPos.size(); i++) {
        isDiscard[snpPosSet[badPos[i]]] = 1;
    }
    
#ifdef LISTBADPOS
    // list out all the bad positions
    cout << "Bad positions:" << endl;
    for (i=0; i<badPos.size(); i++) {
        cout << snpPosSet[badPos[i]] << endl; // "\t" << badPairRatio[i] << "\t" << numPair[i] << endl;
    }
    exit(1);
#endif
    
    // update the array snpPosSet
    k=0;
    for (i=0; i<snpPosSet.size(); i++) {
        if (!isDiscard[snpPosSet[i]]) {
            if (k < i) {
                snpPosSet[k] = snpPosSet[i];
            }
            k++;
        }
    }
    if (k > 0 && k < snpPosSet.size()) {
        cout << (int)snpPosSet.size() - k << " snp positions have been ignored due to bad mutations." << endl;
        snpPosSet.resize(k);
    } else if (k==0) {
        snpPosSet.clear();
    }
    
    delete[] nuclPair;
    delete[] coverage;
}

// randomly generate a random tree
void PhyloRead::computeRandTopology(int nTips, MyRand& myRand) {
    UnlabelTree utree;
    utree.randRootedTree(nTips, currTopology, myRand);
    numTips = nTips;
    numNodes = 2 * numTips - 1;
    numInNodes = numTips - 1;
    numEdges = numNodes - 1;
    /*
    // show current topology
    cout << "current topology: ";
    utree.showTopology(currTopology);
    cout << endl;*/
    isDescendantUpdated = false;
    needUpdated = true;
    buildIsDescendant();
}

// build the array "isDescendant" that indicates the relationship between two nodes
// isDescendant[i*numNodes+j] = true if node j is descendant of node i
void PhyloRead::buildIsDescendant() {
    int lchild, rchild;
    int i,j;
    isDescendant.clear();
    for (i=0; i<numNodes*numNodes; i++)
        isDescendant.push_back(false);
    for (i=numTips; i<numNodes; i++) {
        lchild = currTopology.at(2*i);
        rchild = currTopology.at(2*i+1);
        if (lchild == -1 || rchild == -1) {
            cerr << "[buildIsDescendant] Error! The value of lchild or rchild should not be -1" << endl;
            exit(1);
        }
        isDescendant[i*numNodes + lchild] = true;
        // consider all the descendants of lchild
        for (j=0; j<numNodes; j++) {
            isDescendant[i*numNodes + j] = isDescendant[i*numNodes + j] | isDescendant[lchild*numNodes + j];
        }
        isDescendant[i*numNodes + rchild] = true;
        // consider all the descendants of rchild
        for (j=0; j<numNodes; j++) {
            isDescendant[i*numNodes + j] = isDescendant[i*numNodes + j] | isDescendant[rchild*numNodes + j];
        }
    }
    isDescendantUpdated = true;
    isWeightUpdated = false;
}

// set all parameters to random values
void PhyloRead::randParameters(MyRand& myRand) {
    double avg_hap_freq;
    double sum_hap_freq, fr_hap_freq, to_hap_freq;

    avg_hap_freq = 1.0 / (double) numTips;

    // get a set of random haplotype frequencies
    sum_hap_freq = 1.0;
    fr_hap_freq = avg_hap_freq / 2.0;
    to_hap_freq = avg_hap_freq * 2.0;
    genRandom(fr_hap_freq, to_hap_freq, sum_hap_freq, numTips, haploFreqs_t);
    isHapFreqNormalized = false;
    needUpdated = true;
    
    /*
    haploFreqs_t[0] = 0.035;
    haploFreqs_t[1] = 0.24;
    haploFreqs_t[2] = 0.575;
    haploFreqs_t[3] = 0.06;
    haploFreqs_t[4] = 0.09;
    */
    

}

// set up the variables:
void PhyloRead::setupVariables() {
    nodeWeights.resize(numNodes);
    haploFreqs.resize(numTips);
    haploFreqs_t.resize(numTips);
    rootStates = new int[numPos];
    mutationLocs = new int[numPos];
    expecLogFreqs = new double[numEdges * numEdges * 16];
    expecLogFreq1 = new double[numEdges * 4];
    tloglike = new long double[numEdges * 2];
    // tloglikeValid = new int[numEdges * 2];
    maxMutatePos = new int[numPos*numEdges * 2];
    mutatePos = new int[numPos];
    maxStateRoot = new int[numPos*numEdges * 2];
    stateRoot = new int[numPos];
    log_phi = new double[numPos*numEdges*numEdges * 2];
    lpr_ui = new double[numPos*numEdges];
    charMatHaplo = new double[numPos*numTips*4];
    isValidCharMat = new char[numPos*numTips*4];
    lprobRcMc = new double[numEdges * 2];
    logLikeRiMi = new double[numPos*numEdges*2];
}

// normalize haplotype frequencies
void PhyloRead::normalizeHapFreqs() {
    int i;
    double haploFreqs_sum = 0.0;
    for (i=0; i<numTips; i++)
        haploFreqs_sum += haploFreqs_t[i];
    for (i=0; i<numTips; i++) {
        haploFreqs[i] = haploFreqs_t[i] / haploFreqs_sum;
    }
    isHapFreqNormalized = true;
    isWeightUpdated = false;
}

// update the weight of the nodes
// according to the updated haplotype frequencies
// assume the haplotype frequencies to be normalized
void PhyloRead::updateWeights() {
    int i;
    int leftChild, rightChild;
    // double t;
    
    // update the weights of the nodes
    for (i=0; i<numTips; i++)
        nodeWeights[i] = haploFreqs[i];
    
    for (i=numTips; i<numNodes; i++) {
        leftChild = currTopology.at(2*i);
        rightChild = currTopology.at(2*i+1);
        if (leftChild == -1 || rightChild == -1) {
            cerr << "[updateWeights] Error! The value of leftChild or rightChild should not be -1" << endl;
            exit(1);
        }
        nodeWeights[i] = nodeWeights[leftChild] + nodeWeights[rightChild];
    }
    isWeightUpdated = true;
    isLogFreqUpdated = false;
}

// update the expected frequencies for single mutation on an edge for ONE snp position
// update the array expecLogFreq1
// expecLogFreq1[i*4] = the expected log frequencies when a mutation happens on edge i+1
void PhyloRead::updateExpectLogFreq1() {
    int i;
    int c1,k;
    double* currExpecLogFreq1;
    double f[2];
    double p1;
    memset(expecLogFreq1, 0, numEdges * 4 * sizeof(double));
    
    // cout << "seq error rate: " << seqErr << endl;
    
    for (i=0; i<numEdges; i++) {
        currExpecLogFreq1 = &(expecLogFreq1[i*4]);
        f[0]=f[1]=0.0;
        if (i==0) {
            // no mutation
            f[0] = 1.0;
        } else {
            f[0] = 1.0 - nodeWeights[i-1];
            f[1] = nodeWeights[i-1];
        }
        // cout << f[0] << "," << f[1] << endl;
        // consider the sequencing error
        for (c1=0; c1<2; c1++) {
            for (k=0; k<4; k++) {
                if (k == c1) {
                    p1 = 1.0 - seqErr; // no seq error on snp pos 1
                } else {
                    p1 = seqErr / 3.0; // seq error on snp pos 1
                }
                currExpecLogFreq1[k] += f[c1] * p1;
                // cout << c1 << "," << k << "," << f[c1] * p1 << endl;
            }
        }
        
        // compute the log values
        for (k=0; k<4; k++)
            currExpecLogFreq1[k] = logl(currExpecLogFreq1[k]);
    }
    /*
    // show the array
    for (i=0; i<numEdges; i++) {
        currExpecLogFreq1 = &(expecLogFreq1[i*4]);
        for (k=0; k<4; k++) {
            if (k>0)
                cout << ",";
            cout << currExpecLogFreq1[k];
        }
        cout << endl;
    }
    exit(1);
     */
}

// update the expected frequencies for mutations on different pairs of edges
// update the array expectLogFreqs
// expectLogFreqs[i*numEdges*16 + j*16] = the expected log frequencies when the mutations happen on edges i+1 and j+1
void PhyloRead::updateExpectLogFreqs() {
    double f[4]; // frequencies when no sequencing error
    double p1,p2;
    int i,j,k,l,c1,c2;
    double* currExpecLogFreqs;
    
    // for mutations on single edge
    updateExpectLogFreq1();
    
    memset(expecLogFreqs, 0, numEdges * numEdges * 16 * sizeof(double));
    l=0;
    for (i=0; i<numEdges; i++) {
        for (j=0; j<numEdges; j++) {
            // reset the frequencies
            currExpecLogFreqs = &(expecLogFreqs[l*16]);
            l++;
            memset(f, 0, 4 * sizeof(double));
            if (i==0) {
                if (j==0) {
                    // no mutation on both snp positions
                    f[0] = 1.0;
                } else {
                    // mutation only on the second snp position
                    f[0] = 1 - nodeWeights[j-1];
                    f[1] = nodeWeights[j-1];
                }
            } else {
                if (j==0) {
                    // mutation only on the first snp position
                    f[0] = 1 - nodeWeights[i-1];
                    f[2] = nodeWeights[i-1];
                } else if (i==j) {
                    // mutation on the same edge for both snp positions
                    f[0] = 1 - nodeWeights[i-1];
                    f[3] = nodeWeights[i-1];
                } else {
                    // mutation on different edges for both snp positions
                    if (isDescendant[(i-1)*numNodes+(j-1)]) {
                        // node j-1 is a descendant of node i-1
                        f[0] = 1 - nodeWeights[i-1]; // 00
                        f[2] = nodeWeights[i-1] - nodeWeights[j-1]; // 10
                        f[3] = nodeWeights[j-1]; // 11
                    } else if (isDescendant[(j-1)*numNodes+(i-1)]) {
                        // node i-1 is a descendant of node j-1
                        f[0] = 1 - nodeWeights[j-1]; // 00
                        f[1] = nodeWeights[j-1] - nodeWeights[i-1]; // 01
                        f[3] = nodeWeights[i-1]; // 11
                    } else {
                        f[0] = 1 - nodeWeights[i-1] - nodeWeights[j-1]; // 00
                        f[1] = nodeWeights[j-1]; // 01
                        f[2] = nodeWeights[i-1]; // 10
                    }
                }
            }
            
            /*
            // show the content of f
            cout << "i=" << i << " j=" << j << " f:";
            for (k=0; k<4; k++) {
                if (k>0)
                    cout << " ";
                cout << f[k];
            }
            cout << endl;
            */

            for (c1=0; c1<2; c1++) {
                for (c2=0; c2<2; c2++) {
                    for (k=0; k<16; k++) {
                        if (k/4 == c1) {
                            p1 = 1.0 - seqErr; // no seq error on snp pos 1
                        } else {
                            p1 = seqErr / 3.0; // seq error on snp pos 1
                        }
                        if (k%4 == c2) {
                            p2 = 1.0 - seqErr; // no seq error on snp pos 2
                        } else {
                            p2 = seqErr / 3.0; // seq error on snp pos 2
                        }
                        currExpecLogFreqs[k] += f[c1*2+c2] * p1 * p2;
                    }
                }
            }
            
            // compute the log values
            for (k=0; k<16; k++)
                currExpecLogFreqs[k] = logl(currExpecLogFreqs[k]);
            
            // cout << "i=" << i << " j=" << j << " currExpecLogFreqs:" << endl;
            /*
            for (k=0; k<16; k++) {
                if (k>0)
                    cout << " ";
                cout << currExpecLogFreqs[k];
            }
            cout << endl;*/
        }
    }
    isLogFreqUpdated = true;
}

// compute the value of logCombin (for data set without considering adjacent snp positions)
void PhyloRead::computelogCombin1(int* posmatrix) {
    int p,j;
    vector<int> freqs;
    logCombin1 = 0.0;
    for (p=0; p<numPos; p++) {
        if (!isDiscard[p] && !isInVar[p]) {
        // if (!isDiscard[p]) {
            freqs.clear();
            for (j=0; j<4; j++)
                freqs.push_back(posmatrix[p*5 + j + 1]);
            logCombin1 += logMultiCoeff(freqs);
        }
    }
}

// compute the value of logCombin
void PhyloRead::computelogCombin(int* dmatrix) {
    int p,j;
    vector<int> freqs;
    logCombin = 0.0;
    for (p=0; p<numPos; p++) {
        if (!isDiscard[p] && !isInVar[p]) {
            freqs.clear();
            for (j=0; j<16; j++)
                freqs.push_back(dmatrix[p*16 + j]);
            logCombin += logMultiCoeff(freqs);
        }
    }
}

// compute the log-likelihood of the whole alignment
// by considering the edge lengths
long double PhyloRead::logLike2() {
    long double result = 0.0;
    // long double sub;
    long double value;
    long double c1_value;
    long double m1_value;
    long double sumValue;
    // bool isFirst;
    // bool sub_valid;
    int i,j,k,ii;
    int c1;
    int m1,m2;
    int c1_to,c2_to,m1_to,m2_to;
    int* currDMat;
    double* currExpecLogFrqs;
    int pre_ref = -1;
    
    if (needUpdated) {
        if (!isDescendantUpdated)
            buildIsDescendant();
        
        if (!isHapFreqNormalized)
            normalizeHapFreqs();
        
        if (!isWeightUpdated)
            updateWeights();
        
        if (!isLogFreqUpdated)
            updateExpectLogFreqs();
    }
    needUpdated = false;
    
    // compute the phi array
    computePhi();
    
    // compute the edge lengths
    computeEdgeLen();
    
    // memset(tloglikeValid, 0, numEdges*2*sizeof(int));
    c1_to = m1_to = 0;
    for (i=0; i<numPos; i++) {
        if (!isDiscard[i] && !isInVar[i]) {
            ii = i * numEdges * numEdges * 2;
            j = refPos[i];
            if (j != pre_ref) {
                if (pre_ref != -1) {
                    // max value out of tloglike[1...n] where n=2*numEdges
                    // c1*numEdges + m1
                    for (m1=0; m1<m1_to; m1++) {
                        for (c1=0; c1<c1_to; c1++) {
                            k = c1*numEdges + m1;
                            value = tloglike[k];
                            if (c1==0) {
                                c1_value = value;
                            } else {
                                c1_value = log_x_plus_y(c1_value, value);
                            }
                        }
                        if (m1==0) {
                            m1_value = c1_value + logEdgeLen[m1];
                        } else {
                            m1_value = log_x_plus_y(m1_value, c1_value + logEdgeLen[m1]);
                        }
                    }
                    result += m1_value;
                }
                // memset(tloglikeValid, 0, numEdges*2*sizeof(int));
                memset(tloglike, 0, numEdges*2*sizeof(long double));
                // set up m1_to and c1_to
                if (isInVar[j]) {
                    m1_to = c1_to = 1;
                } else {
                    m1_to = numEdges;
                    c1_to = 2;
                }
                pre_ref = j;
            }
            // let m1 = the mutation location on the ref position
            // let m2 = the mutation location on the current position
            // let c1 = the root state on the ref position
            // let c2 = the root state on the current position
            
            if (i==j) {
                for (c1=0; c1<c1_to; c1++) {
                    for (m1=0; m1<m1_to; m1++) {
                        currExpecLogFrqs = &(expecLogFreq1[m1*4]);
                        currDMat = &(dMat[i*64 + c1*4]);
                        // tloglikeValid[c1*numEdges + m1] = 1;
                        tloglike[c1*numEdges + m1] += log_phi[ii + m1*numEdges*2 + c1];
                    }
                }
            } else {
                
                // set up m2_to and c2_to
                if (isInVar[i]) {
                    m2_to = c2_to = 1;
                } else {
                    m2_to = numEdges;
                    c2_to = 2;
                }
                for (c1=0; c1<c1_to; c1++) {
                    for (m1=0; m1<m1_to; m1++) {
                        for (m2=0; m2<m2_to; m2++) {
                            value = log_phi[ii + m2*numEdges*2 + m1*2 + c1] + logEdgeLen[m2];
                            if (m2==0) {
                                sumValue = value;
                            } else {
                                sumValue = log_x_plus_y(sumValue, value);
                            }
                        }
                        tloglike[c1*numEdges + m1] += sumValue;
                        // tloglikeValid[c1*numEdges + m1] = 1;
                    }
                }
            }
        }
    }
    
    if (pre_ref != -1) {
        // max value out of tloglike[1...n] where n=2*numEdges
        // c1*numEdges + m1
        for (m1=0; m1<m1_to; m1++) {
            for (c1=0; c1<c1_to; c1++) {
                k = c1*numEdges + m1;
                value = tloglike[k];
                if (c1==0) {
                    c1_value = value;
                } else {
                    c1_value = log_x_plus_y(c1_value, value);
                }
            }
            if (m1==0) {
                m1_value = c1_value + logEdgeLen[m1];
            } else {
                m1_value = log_x_plus_y(m1_value, c1_value + logEdgeLen[m1]);
            }
        }
        result += m1_value;
    }
    
    return result + logCombin;
}

long double PhyloRead::logLikeS() {
    return logLikeS(0);
}

long double PhyloRead::logLikeS(int obtainMutatePos) {
    // long double result = 0.0;
    double* currExpecLogFrqs;
    int* currDMat;
    int i, k, c1, m1;
    long double maxLogLike;
    long double currLogLike, maxLogLike_m1, result;
    long double sumLogLike;

    /*
    if (needUpdated) {
        if (!isDescendantUpdated)
            buildIsDescendant();
        
        if (!isHapFreqNormalized)
            normalizeHapFreqs();
        
        if (!isWeightUpdated)
            updateWeights();
        
        if (!isLogFreqUpdated) {
            updateExpectLogFreq1();
            isLogFreqUpdated = true;
        }
    }
     */

    buildIsDescendant();
    normalizeHapFreqs();
    updateWeights();
    updateExpectLogFreq1();

    // memset(tloglikeValid, 0, numEdges*2*sizeof(int));
    if (obtainMutatePos) {
        memset(mutatePos, 0, numPos*sizeof(int));
        memset(stateRoot, 0, numPos*sizeof(int));
    }

    needUpdated = false;
    result = 0.0;
    for (i=0; i<numPos; i++) {
        if (!isDiscard[i] && !isInVar[i]) {
            for (c1=0; c1<2; c1++) {
                currDMat = &(sdMat[i*8 + c1*4]);
                maxLogLike_m1 = 0.0;
                for (m1=0; m1<numEdges; m1++) {
                    currExpecLogFrqs = &(expecLogFreq1[m1*4]);
                    currLogLike = 0.0;
                    for (k=0; k<4; k++) {
                        currLogLike += currDMat[k] * currExpecLogFrqs[k];
                    }
                    if ((c1==0 && m1==0) || (currLogLike > maxLogLike)) {
                        maxLogLike = currLogLike;
                        if (obtainMutatePos) {
                            stateRoot[i] = c1;
                            mutatePos[i] = m1;
                        }
                    }
                    if (m1==0 || currLogLike > maxLogLike_m1) {
                        maxLogLike_m1 = currLogLike;
                    }
                }
                if (c1==0) {
                    sumLogLike = maxLogLike_m1 + logStateFreqs[charOrder[i*4+c1]];
                } else {
                    sumLogLike = log_x_plus_y(sumLogLike, maxLogLike_m1 + logStateFreqs[charOrder[i*4+c1]]);
                }
            }
            result += sumLogLike;
            // result += maxLogLike;
        /*
        } if (isInVar[i]) {
            currDMat = &(sdMat[i*8]);
            currExpecLogFrqs = &(expecLogFreq1[0]);
            currLogLike = 0.0;
            for (k=0; k<4; k++) {
                currLogLike += currDMat[k] * currExpecLogFrqs[k];
            }
            result += currLogLike + logStateFreqs[charOrder[i*4]];*/
        } else if (isDiscard[i] && obtainMutatePos) {
            mutatePos[i] = -1;
            stateRoot[i] = -1;
        }
    }
    // for debugging
    /*
    if (obtainMutatePos) {
        for (i=0; i<numPos; i++) {
            if (i>0)
                cout << ",";
            cout << mutatePos[i];
        }
        cout << endl;
    }*/
    // cout << "final sum = " << result + logCombin1 << endl;
    return result + logCombin1;
}

// compute the log-likelihood for the whole alignment
// obtainMutatePos=1 if want to obtain the optimal mutation positions
long double PhyloRead::logLike() {
    return logLike(0);
}

// compute the log-likelihood for the whole alignment
// obtainMutatePos=1 if want to obtain the optimal mutation positions
long double PhyloRead::logLike(int obtainMutatePos) {
    long double result = 0.0;
    // long double sub;
    long double value;
    long double m1_value;
    long double m2_value;
    long double maxValue;
    long double sumValue;
    bool isFirst;
    // bool sub_valid;
    int i,j,k,pre_i,ii;
    int c1,c2,max_c1,max_c2;
    int m1,m2,max_m1,max_m2;
    int c1_to,c2_to,m1_to,m2_to;
    int* currDMat;
    double* currExpecLogFrqs;
    int pre_ref = -1;

    if (needUpdated) {
        if (!isDescendantUpdated)
            buildIsDescendant();
        
        if (!isHapFreqNormalized)
            normalizeHapFreqs();
        
        if (!isWeightUpdated)
            updateWeights();
        
        if (!isLogFreqUpdated)
            updateExpectLogFreqs();
    }

    needUpdated = false;
    
    // memset(tloglikeValid, 0, numEdges*2*sizeof(int));
    if (obtainMutatePos) {
        memset(maxMutatePos, 0, numPos*numEdges*2*sizeof(int));
        memset(mutatePos, 0, numPos*sizeof(int));
        memset(maxStateRoot, 0, numPos*numEdges*2*sizeof(int));
        memset(stateRoot, 0, numPos*sizeof(int));
    }
    c1_to = m1_to = 0;
    pre_i = 0;
    for (i=0; i<numPos; i++) {
        if (!isDiscard[i] && !isInVar[i]) {
            j = refPos[i];
            if (j != pre_ref) {
                if (pre_ref != -1) {
                    // max value out of tloglike[1...n] where n=2*numEdges
                    // c1*numEdges + m1
                    sumValue = 0.0;
                    isFirst = true;
                    for (c1=0; c1<c1_to; c1++) {
                        for (m1=0; m1<m1_to; m1++) {
                            k = c1*numEdges + m1;
                            value = tloglike[k];
                            if (m1 == 0) {
                                m1_value = value;
                            } else if (value > m1_value) {
                                m1_value = value;
                            }
                            if (obtainMutatePos && (isFirst || value > maxValue)) {
                                maxValue = value;
                                max_m1 = m1;
                                max_c1 = c1;
                                isFirst = false;
                            }
                        }
                        if (c1 == 0)
                            sumValue = m1_value + logStateFreqs[charOrder[pre_ref*4+c1]];
                        else
                            sumValue = log_x_plus_y(sumValue, m1_value + logStateFreqs[charOrder[pre_ref*4+c1]]);
                    }
                    result += sumValue;
                    if (obtainMutatePos) {
                        for (ii=pre_i; ii<i; ii++) {
                            if (isInVar[ii]) {
                                mutatePos[ii] = 0;
                                stateRoot[ii] = 0;
                            } else if (isDiscard[ii]) {
                                mutatePos[ii] = -1;
                                stateRoot[ii] = -1;
                            } else if (refPos[ii] != pre_ref) {
                                mutatePos[ii] = -1;
                                stateRoot[ii] = -1;
                            } else if (ii == pre_ref) {
                                mutatePos[ii] = max_m1;
                                stateRoot[ii] = max_c1;
                            } else {
                                mutatePos[ii] = maxMutatePos[ii * numEdges * 2 + max_m1*2 + max_c1];
                                stateRoot[ii] = maxStateRoot[ii * numEdges * 2 + max_m1*2 + max_c1];
                            }
                        }
                    }
                    pre_i = i;
                }
                // memset(tloglikeValid, 0, numEdges*2*sizeof(int));
                memset(tloglike, 0, numEdges*2*sizeof(long double));
                // set up m1_to and c1_to
                if (isInVar[j]) {
                    m1_to = c1_to = 1;
                } else {
                    m1_to = numEdges;
                    c1_to = 2;
                }
                pre_ref = j;
            }
            // let m1 = the mutation location on the ref position
            // let m2 = the mutation location on the current position
            // let c1 = the root state on the ref position
            // let c2 = the root state on the current position

            if (i==j) {
                for (c1=0; c1<c1_to; c1++) {
                    for (m1=0; m1<m1_to; m1++) {
                        currExpecLogFrqs = &(expecLogFreq1[m1*4]);
                        currDMat = &(dMat[i*64 + c1*4]);
                        for (k=0; k<4; k++) {
                            tloglike[c1*numEdges + m1] += currDMat[k] * currExpecLogFrqs[k];
                        }
                    }
                }
            } else {

                // set up m2_to and c2_to
                if (isInVar[i]) {
                    m2_to = c2_to = 1;
                } else {
                    m2_to = numEdges;
                    c2_to = 2;
                }
                for (c1=0; c1<c1_to; c1++) {
                    for (m1=0; m1<m1_to; m1++) {
                        isFirst = true;
                        sumValue = 0.0;
                        for (c2=0; c2<c2_to; c2++) {
                            currDMat = &(dMat[i*64 + c1*32 + c2*16]);
                            for (m2=0; m2<m2_to; m2++) {
                                currExpecLogFrqs = &(expecLogFreqs[m1*numEdges*16+m2*16]);
                                value = 0.0;
                                for (k=0; k<16; k++) {
                                    value += currExpecLogFrqs[k] * currDMat[k];
                                }
                                if (m2 == 0) {
                                    m2_value = value;
                                } else if (value > m2_value) {
                                    m2_value = value;
                                }
                                if (obtainMutatePos && (isFirst || value > maxValue)) {
                                    maxValue = value;
                                    max_c2 = c2;
                                    max_m2 = m2;
                                    isFirst = false;
                                }
                            }
                            if (c2 == 0)
                                sumValue = m2_value + logStateFreqs[charOrder[i*4+c2]];
                            else
                                sumValue = log_x_plus_y(sumValue, m2_value + logStateFreqs[charOrder[i*4+c2]]);
                        }
                        tloglike[c1*numEdges + m1] += sumValue;
                        if (obtainMutatePos) {
                            maxMutatePos[i * numEdges * 2 + m1*2 + c1] = max_m2;
                            maxStateRoot[i * numEdges * 2 + m1*2 + c1] = max_c2;
                        }
                    }
                }
            }
        }
    }

    if (pre_ref != -1) {
        sumValue = 0.0;
        isFirst = true;
        for (c1=0; c1<c1_to; c1++) {
            for (m1=0; m1<m1_to; m1++) {
                k = c1*numEdges + m1;
                value = tloglike[k];
                if (m1 == 0) {
                    m1_value = value;
                } else if (value > m1_value) {
                    m1_value = value;
                }
                if (obtainMutatePos && (isFirst || value > maxValue)) {
                    maxValue = value;
                    max_m1 = m1;
                    max_c1 = c1;
                    isFirst = false;
                }
            }
            if (c1 == 0)
                sumValue = m1_value + logStateFreqs[charOrder[pre_ref*4+c1]];
            else
                sumValue = log_x_plus_y(sumValue, m1_value + logStateFreqs[charOrder[pre_ref*4+c1]]);
        }
        result += sumValue;
    }
    
    if (obtainMutatePos) {
        for (ii=pre_i; ii<numPos; ii++) {
            if (isInVar[ii]) {
                mutatePos[ii] = 0;
                stateRoot[ii] = 0;
            } else if (isDiscard[ii]) {
                mutatePos[ii] = -1;
                stateRoot[ii] = -1;
            } else if (refPos[ii] != pre_ref) {
                mutatePos[ii] = -1;
                stateRoot[ii] = -1;
            } else if (ii == pre_ref) {
                mutatePos[ii] = max_m1;
                stateRoot[ii] = max_c1;
            } else {
                mutatePos[ii] = maxMutatePos[ii * numEdges * 2 + max_m1*2 + max_c1];
                stateRoot[ii] = maxStateRoot[ii * numEdges * 2 + max_m1*2 + max_c1];
            }
        }
    }
    return result + logCombin;
}

// compute the array phi
void PhyloRead::computePhi() {
    int m1, c1, m1_to, c1_to, m2_to, c2_to;
    // let m1 = the mutation location on the ref position
    // let c1 = the root state on the ref position
    // phi[i*numEdges*numEdges*2 + y*numEdges*2 + a*2 + x] = phi(i,y,a,x)
    int i,j,k,ii;
    int m2,c2;
    int pre_ref;
    int* currDMat;
    double* currExpecLogFrqs;
    long double currLogLike;
    long double logLike_c2;
    
    // for the ref positions (i.e. i = ref(H) )
    m1_to = numEdges;
    c1_to = 2;
    for (j=0; j<refLst.size(); j++) {
        i = refLst[j]; // ref pos
        ii = i * numEdges * numEdges * 2;
        if (isInVar[i]) {
            log_phi[ii] = 0.0;
            continue;
        }
        // variable site
        for (m1=0; m1<m1_to; m1++) {
            for (c1=0; c1<c1_to; c1++) {
                currExpecLogFrqs = &(expecLogFreq1[m1*4]);
                currDMat = &(dMat[i*64 + c1*4]);
                currLogLike = 0.0;
                for (k=0; k<4; k++) {
                    currLogLike += currDMat[k] * currExpecLogFrqs[k];
                }
                currLogLike += logStateFreqs[charOrder[i*4+c1]];
                log_phi[ii + m1*numEdges*2 + c1] = currLogLike;
            }
        }
    }
    
    // for the other positions (i.e. i != ref(H))
    pre_ref = -1;
    m2_to = numEdges;
    c2_to = 2;
    for (i=0; i<numPos; i++) {
        if (isInVar[i] || isDiscard[i]) {
            continue;
        }
        j = refPos[i];
        if (j==i) {
            // this is a reference position, skip it
            continue;
        }
        if (j != pre_ref) {
            // set up m1_to and c1_to
            if (isInVar[j]) {
                m1_to = c1_to = 1;
            } else {
                m1_to = numEdges;
                c1_to = 2;
            }
            pre_ref = j;
        }

        ii = i * numEdges * numEdges * 2;
        for (m2=0; m2<m2_to; m2++) {
            for (m1=0; m1<m1_to; m1++) {
                for (c1=0; c1<c1_to; c1++) {
                    for (c2=0; c2<c2_to; c2++) {
                        currExpecLogFrqs = &(expecLogFreqs[m1*numEdges*16+m2*16]);
                        currDMat = &(dMat[i*64 + c1*32 + c2*16]);
                        currLogLike = 0.0;
                        for (k=0; k<16; k++) {
                            currLogLike += currExpecLogFrqs[k] * currDMat[k];
                        }
                        currLogLike += logStateFreqs[charOrder[i*4+c2]];
                        if (c2==0) {
                            logLike_c2 = currLogLike;
                        } else {
                            logLike_c2 = log_x_plus_y(logLike_c2, currLogLike);
                        }
                    }
                    log_phi[ii + m2*numEdges*2 + m1*2 + c1] = logLike_c2;
                }
            }
        }
    }
}

// compute the edge lengths
void PhyloRead::computeEdgeLen() {
    int i, m;
    int firstVarSite;
    int numSites;
    long double sumEdgeLen;
    edgeLen.clear();
    logEdgeLen.clear();
    for (m=0; m<numEdges; m++) {
        logEdgeLen.push_back(0.0);
    }
    
    // for the first variable site
    firstVarSite=0;
    while (firstVarSite<numPos && (isInVar[firstVarSite] || isDiscard[firstVarSite]))
        firstVarSite++;
    if (firstVarSite>=numPos) {
        cerr << "Error! There is no variable site" << endl;
        exit(1);
    }
    
    // compute LOG Pr(u_c = x | theta)
    // i.e. log Pr(location of mutation on the reference position c = x | the parameters)
    computeLPrUc();

    // compute LOG Pr(u_i = x | theta)
    // i.e. log Pr(location of mutation on the position i = x | the parameters)
    computeLPrUi();

    for (m=0; m<numEdges; m++) {
        logEdgeLen[m] = lpr_ui[firstVarSite*numEdges+m];
    }
    numSites = 1;
    
    for (i=0; i<numPos; i++) {
        if (i == firstVarSite)
            continue;
        if (isDiscard[i])
            continue;
        numSites++;
        if (isInVar[i]) {
            logEdgeLen[0] = log_x_plus_y(logEdgeLen[0], 0);
        } else {
            for (m=0; m<numEdges; m++) {
                logEdgeLen[m] = log_x_plus_y(logEdgeLen[m], lpr_ui[i*numEdges+m]);
            }
        }
    }
    
    sumEdgeLen = logEdgeLen[0];
    for (m=1; m<numEdges; m++) {
        sumEdgeLen = log_x_plus_y(logEdgeLen[m], sumEdgeLen);
    }
    
    // normalize
    for (m=0; m<numEdges; m++) {
        logEdgeLen[m] = logEdgeLen[m] - sumEdgeLen;
    }
    
    // remove the first one and set the last one to zero
    for (m=0; m<numEdges-1; m++) {
        edgeLen.push_back(expl(logEdgeLen[m+1]));
    }
    edgeLen.push_back(0.0);
}

// characters of each tip
// 0 - same as the root state; 1 - different from the root state
// size of tipChars is initialized as the number of tips
void PhyloRead::getChars(int mutateLoc, vector<int>& tipChars) {
    vector<int> tList;
    int i,k;
    // reset all items in tipChars to 0
    for (i=0; i<tipChars.size(); i++) {
        tipChars[i] = 0;
    }
    if (mutateLoc == 0) {
        // no mutation
        return;
    }
    tList.push_back(mutateLoc-1);
    i=0;
    while (i < tList.size()) {
        k = tList[i];
        if (k < tipChars.size()) {
            if (k >= 0) {
                tipChars[k] = 1;
            }
        } else {
            tList.push_back(currTopology[2*k]);
            tList.push_back(currTopology[2*k+1]);
        }
        i++;
    }
}

// get the current sequences
void PhyloRead::getCurrSeqs(vector<string>& seqs) {
    int i,j,k1,k2;
    char c,c2;
    char int2nucl[4] = {'A','C','G','T'};
    vector<int> tipChars;
    seqs.clear();
    for (i=0; i<numTips; i++) {
        seqs.push_back("");
    }
    tipChars.resize(numTips);
    for (i=0; i<numPos; i++) {
        if (isInVar[i] || isDiscard[i] || mutatePos[i] < 0) {
            if (isInVar[i])
                c = int2nucl[charOrder[i*4]];
            else
                c = 'N';
            for (j=0; j<numTips; j++) {
                seqs[j].append(1,c);
            }
            // cout << i << "\t" << c << endl;
        } else {
            // root state
            if (stateRoot[i] == 0) {
                k1=charOrder[i*4];
                k2=charOrder[i*4 + 1];
            } else {
                k2=charOrder[i*4];
                k1=charOrder[i*4 + 1];
            }
            if (k1 < 4)
                c = int2nucl[k1];
            else
                c = 'N';
            if (k2 < 4)
                c2 = int2nucl[k2];
            else
                c2 = 'N';
            // get the character for each tip
            getChars(mutatePos[i], tipChars);
            // cout << i << "\t" << mutatePos[i] << endl;
            for (j=0; j<numTips; j++) {
                if (tipChars[j] == 0) {
                    // same as the root character
                    seqs[j].append(1,c);
                } else {
                    // different from the root character
                    seqs[j].append(1,c2);
                }
            }
        }
    }
}

// show the current sequences
void PhyloRead::showCurrSeqs(ofstream& fout, int round) {
    int i,j,l;
    vector<string> seqs;
    
    getCurrSeqs(seqs);
    // show the sequences
    for (i=0; i<numTips; i++) {
        fout << ">" << round << "_" << i+1 << "_f=" << haploFreqs[i] << endl;
        for (j=0; j<seqs[i].length(); j+=60) {
            l = 60;
            if (j+l > seqs[i].length()) {
                l = seqs[i].length() - j;
            }
            fout << seqs[i].substr(j, l) << endl;
        }
    }
}

// show the current tree
string PhyloRead::currTree(int num_digits) {
    int i;
    int tot = 0;
    vector<int> sum;
    for (i=0; i<numEdges; i++)
        sum.push_back(0);
    for (i=0; i<numPos; i++) {
        if (!isDiscard[i]) {
            sum[mutatePos[i]]++;
            tot++;
        }
    }
    edgelens.clear();
    for (i=1; i<numEdges; i++)
        edgelens.push_back((double)sum[i] / tot);
    edgelens.push_back(0.0);
    normalizeHapFreqs();
    return tree.topInt2Txt(currTopology, haploFreqs, edgelens, num_digits);
}

// show the current tree from the array edgeLen
string PhyloRead::currTree2(int num_digits) {
    /*
    cout << "currTopology" << endl << flush;
    for (int i=0; i<currTopology.size(); i+=2) {
        cout << currTopology[i] << "," << currTopology[i+1] << endl << flush;
    }
    cout << "haploFreqs" << endl << flush;
    for (int i=0; i<haploFreqs.size(); i++) {
        cout << haploFreqs[i] << endl;
    }*/
    return tree.topInt2Txt(currTopology, haploFreqs, edgeLen, num_digits);
}

// show the current tree with all edge length = 1.0
string PhyloRead::currTree1(int num_digits) {
    int i;
    string s;
    edgelens.clear();
    for (i=1; i<numEdges; i++)
        edgelens.push_back(1.0);
    edgelens.push_back(0.0);
    s = tree.topInt2Txt(currTopology, haploFreqs, edgelens, num_digits);
    return s;
}

//----------------------------
// supplementary functions
//----------------------------

void PhyloRead::startup() {
    
    nucl2int = new int[256];
    memset(nucl2int, 0, 256*sizeof(int));
    nucl2int[(int)'A'] = nucl2int[(int)'a'] = 1;
    nucl2int[(int)'C'] = nucl2int[(int)'c'] = 2;
    nucl2int[(int)'G'] = nucl2int[(int)'g'] = 3;
    nucl2int[(int)'T'] = nucl2int[(int)'t'] = 4;
    
    seqErr = INIT_SEQ_ERR;
    
    // initialize the arrays
    charOrder = NULL; // the characters (i.e. 0,1,2,3) in decreasing order according to freq for each position
    // dimension: numPos * 4
    isInVar = NULL; // whether the site is invariable
    isDiscard = NULL; // whether the site is discarded
    expecLogFreqs = NULL; // the expected log frequencies for different locations of the mutations on TWO snp positions
    expecLogFreq1 = NULL; // the expected log frequencies for different location of the mutation on ONE snp position
    // the state of the root for each position
    rootStates = NULL;
    // the mutation location for each position
    mutationLocs = NULL;
    
    // initialize the arrays for computation of likelihood
    tloglike = NULL;
    // tloglikeValid = NULL;
    maxMutatePos = NULL;
    mutatePos = NULL;
    maxStateRoot = NULL;
    stateRoot = NULL;
    log_phi = NULL;
    lpr_ui = NULL;
    charMatHaplo = NULL;
    isValidCharMat = NULL;
    lprobRcMc = NULL;
    logLikeRiMi = NULL;
    dMat = NULL;
    sdMat = NULL;

    //--------------------------------------
    // status
    //--------------------------------------
    isDescendantUpdated = false;
    isHapFreqNormalized = false;
    isWeightUpdated = false;
    isLogFreqUpdated = false;
    needUpdated = true;
    
    //------------------------------------------------------
    // which method to estimate the back-mutation positions
    //------------------------------------------------------
    back_mutate_method = DEFAULT_BACK_MUTATE_METHOD;
    
    // using the matrices of other object
    useOtherMat = false;

}

// compute LOG Pr(u_c = x | theta)
// i.e. log Pr(location of mutation on the reference position c = x | the parameters)
void PhyloRead::computeLPrUc() {
    int m1, c1, m1_to, c1_to;
    // let m1 = the mutation location on the ref position
    // let c1 = the root state on the ref position
    int i,j,ii;
    long double currLogLike;
    long double logLike;
    long double logSumLike;
    
    m1_to = numEdges;
    c1_to = 2;
    for (j=0; j<refLst.size(); j++) {
        i = refLst[j]; // ref pos
        ii = i * numEdges * numEdges * 2;
        if (isInVar[i]) {
            lpr_ui[i*numEdges] = 0.0;
            continue;
        }
        // compute likelihood
        for (m1=0; m1<m1_to; m1++) {
            for (c1=0; c1<c1_to; c1++) {
                currLogLike = log_phi[ii + m1*numEdges*2 + c1];
                if (c1==0) {
                    logLike = currLogLike;
                } else {
                    logLike = log_x_plus_y(logLike, currLogLike);
                }
            }
            lpr_ui[i*numEdges + m1] = logLike;
            if (m1==0) {
                logSumLike = logLike;
            } else {
                logSumLike = log_x_plus_y(logSumLike, logLike);
            }
        }
        // normalize
        for (m1=0; m1<m1_to; m1++) {
            lpr_ui[i*numEdges + m1] -= logSumLike;
        }
    }
}

// compute LOG Pr(u_i = x | theta)
// i.e. log Pr(location of mutation on the position i = x | the parameters)
void PhyloRead::computeLPrUi() {
    int c1,m1,m2;
    // let m1 = the mutation location on the ref position
    // let c1 = the root state on the ref position
    // let m2 = the mutation location on the curr position
    // let c2 = the root state on the curr position
    int c1_to,c2_to,m1_to,m2_to;
    int i,j,ii;
    int pre_ref;
    long double logLike_c2;
    long double logLike_c1;
    long double logLike_m1;
    long double logSumLike;
    
    pre_ref = -1;
    m2_to = numEdges;
    c2_to = 2;
    for (i=0; i<numPos; i++) {
        if (isInVar[i] || isDiscard[i]) {
            continue;
        }
        j = refPos[i];
        if (j==i) {
            // this is a reference position, skip it
            continue;
        }
        if (j != pre_ref) {
            // set up m1_to and c1_to
            if (isInVar[j]) {
                m1_to = c1_to = 1;
            } else {
                m1_to = numEdges;
                c1_to = 2;
            }
            pre_ref = j;
        }
        ii = i * numEdges * numEdges * 2;
        for (m2=0; m2<m2_to; m2++) {
            for (m1=0; m1<m1_to; m1++) {
                for (c1=0; c1<c1_to; c1++) {
                    /*
                     for (c2=0; c2<c2_to; c2++) {
                     currExpecLogFrqs = &(expecLogFreqs[m1*numEdges*16+m2*16]);
                     currDMat = &(dMat[i*64 + c1*32 + c2*16]);
                     currLogLike = 0.0;
                     for (k=0; k<16; k++) {
                     currLogLike += currExpecLogFrqs[k] * currDMat[k];
                     }
                     currLogLike += logStateFreqs[charOrder[i*4+c2]];
                     if (c2==0) {
                     logLike_c2 = currLogLike;
                     } else {
                     logLike_c2 = log_x_plus_y(logLike_c2, currLogLike);
                     }
                     }
                     */
                    logLike_c2 = log_phi[ii + m2*numEdges*2 + m1*2 + c1];
                    
                    logLike_c2 += logStateFreqs[charOrder[j*4+c1]];
                    
                    if (c1==0) {
                        logLike_c1 = logLike_c2;
                    } else {
                        logLike_c1 = log_x_plus_y(logLike_c1, logLike_c2);
                    }
                }
                logLike_c1 += lpr_ui[j*numEdges+m1];
                if (m1==0) {
                    logLike_m1 = logLike_c1;
                } else {
                    logLike_m1 = log_x_plus_y(logLike_m1, logLike_c1);
                }
            }
            lpr_ui[i*numEdges+m2] = logLike_m1;
            if (m2==0) {
                logSumLike = logLike_m1;
            } else {
                logSumLike = log_x_plus_y(logSumLike, logLike_m1);
            }
        }
        // normalize
        for (m2=0; m2<m2_to; m2++) {
            lpr_ui[i*numEdges+m2] -= logSumLike;
        }
    }
}

// compute the loglikelihood given mi and ri for all positions
void PhyloRead::computeLogLikeRiMi() {
    int* currDMat;
    double* currExpecLogFrqs;
    int i,k;
    int c,pre_c; // the current and previous ref pos
    int r1,m1,r1_to,m1_to; // the root state and the mutation location of the reference position
    int r2,m2,r2_to,m2_to; // the root state and the mutation location of the current position
    int max_r2,max_m2;
    double max_lprob2;
    double currLogLike;
    double sumLogLike;
    
    // initialize logLikeRiMi
    memset(logLikeRiMi, 0, numPos * numEdges * 2 * sizeof(double));

    r1_to = m1_to = 0;
    pre_c = -1;
    
    for (i=0; i<numPos; i++) {
        if (isDiscard[i]) {
            continue;
        }
        if (isInVar[i]) {
            continue;
        }
        c = refPos[i];
        if (c!=pre_c) {
            if (pre_c != -1) {
                // compute the log likelihood for the previous window
                for (r1=0; r1<r1_to; r1++) {
                    for (m1=0; m1<m1_to; m1++) {
                        currExpecLogFrqs = &(expecLogFreq1[m1*4]);
                        currDMat = &(dMat[pre_c*64 + r1*4]);
                        currLogLike = 0.0;
                        for (k=0; k<4; k++) {
                            currLogLike += currDMat[k] * currExpecLogFrqs[k];
                        }
                        logLikeRiMi[pre_c * numEdges * 2 + r1 * numEdges + m1] += currLogLike;
                        if (r1==0 && m1==0) {
                            sumLogLike = logLikeRiMi[pre_c * numEdges * 2 + r1 * numEdges + m1];
                        } else {
                            sumLogLike = log_x_plus_y(sumLogLike, logLikeRiMi[pre_c * numEdges * 2 + r1 * numEdges + m1]);
                        }
                    }
                }
                // normalize
                for (r1=0; r1<r1_to; r1++) {
                    for (m1=0; m1<m1_to; m1++) {
                        logLikeRiMi[pre_c * numEdges * 2 + r1 * numEdges + m1] -= sumLogLike;
                    }
                }
            }
            // a new window
            // set up m1_to and r1_to
            if (isInVar[c]) {
                r1_to = m1_to = 1;
            } else {
                r1_to = 2;
                m1_to = numEdges;
            }
            pre_c = c;
        }
        
        if (i==c) {
            continue;
        } else {
            if (isInVar[i]) {
                m2_to = r2_to = 1;
            } else {
                m2_to = numEdges;
                r2_to = 2;
            }
            for (r1=0; r1<r1_to; r1++) {
                for (m1=0; m1<m1_to; m1++) {
                    for (r2=0; r2<r2_to; r2++) {
                        for (m2=0; m2<m2_to; m2++) {
                            currDMat = &(dMat[i*64 + r1*32 + r2*16]);
                            currExpecLogFrqs = &(expecLogFreqs[m1*numEdges*16+m2*16]);
                            currLogLike = 0.0;
                            for (k=0; k<16; k++) {
                                currLogLike += currExpecLogFrqs[k] * currDMat[k];
                            }
                            if ((r2==0 && m2==0) || (max_lprob2 < currLogLike)) {
                                max_r2 = r2;
                                max_m2 = m2;
                                max_lprob2 = currLogLike;
                            }
                        }
                    }
                    logLikeRiMi[c * numEdges * 2 + r1 * numEdges + m1] += max_lprob2;
                }
            }
        }
    }
    
    // compute the log likelihood for the last previous window
    if (pre_c != -1) {
        for (r1=0; r1<r1_to; r1++) {
            for (m1=0; m1<m1_to; m1++) {
                currExpecLogFrqs = &(expecLogFreq1[m1*4]);
                currDMat = &(dMat[pre_c*64 + r1*4]);
                currLogLike = 0.0;
                for (k=0; k<4; k++) {
                    currLogLike += currDMat[k] * currExpecLogFrqs[k];
                }
                logLikeRiMi[pre_c * numEdges * 2 + r1 * numEdges + m1] += currLogLike;
                if (r1==0 && m1==0) {
                    sumLogLike = logLikeRiMi[pre_c * numEdges * 2 + r1 * numEdges + m1];
                } else {
                    sumLogLike = log_x_plus_y(sumLogLike, logLikeRiMi[pre_c * numEdges * 2 + r1 * numEdges + m1]);
                }
            }
        }
        // normalize
        for (r1=0; r1<r1_to; r1++) {
            for (m1=0; m1<m1_to; m1++) {
                logLikeRiMi[pre_c * numEdges * 2 + r1 * numEdges + m1] -= sumLogLike;
            }
        }
    }
    
    // compute the log likelihood for the other positions
    pre_c = -1;
    for (i=0; i<numPos; i++) {
        if (isDiscard[i]) {
            continue;
        }
        if (isInVar[i]) {
            continue;
        }
        c = refPos[i];
        if (i == c)
            continue;
        if (c != pre_c) {
            // a new window
            // set up m1_to and r1_to
            if (isInVar[c]) {
                r1_to = m1_to = 1;
            } else {
                r1_to = 2;
                m1_to = numEdges;
            }
            pre_c = c;
        }
        for (r2=0; r2<r2_to; r2++) {
            for (m2=0; m2<m2_to; m2++) {
                for (r1=0; r1<r1_to; r1++) {
                    for (m1=0; m1<m1_to; m1++) {
                        currDMat = &(dMat[i*64 + r1*32 + r2*16]);
                        currExpecLogFrqs = &(expecLogFreqs[m1*numEdges*16+m2*16]);
                        currLogLike = 0.0;
                        for (k=0; k<16; k++) {
                            currLogLike += currExpecLogFrqs[k] * currDMat[k];
                        }
                        currLogLike += logLikeRiMi[c * numEdges * 2 + r1 * numEdges + m1];
                        if (r1==0 && m1==0) {
                            logLikeRiMi[i * numEdges * 2 + r2 * numEdges + m2] = currLogLike;
                        } else {
                            logLikeRiMi[i * numEdges * 2 + r2 * numEdges + m2] = log_x_plus_y(logLikeRiMi[i * numEdges * 2 + r2 * numEdges + m2], currLogLike);
                        }
                    }
                }
                if (r2==0 && m2==0) {
                    sumLogLike = logLikeRiMi[i * numEdges * 2 + r2 * numEdges + m2];
                } else {
                    sumLogLike = log_x_plus_y(sumLogLike, logLikeRiMi[i * numEdges * 2 + r2 * numEdges + m2]);
                }
            }
        }
        // normalize
        for (r2=0; r2<r2_to; r2++) {
            for (m2=0; m2<m2_to; m2++) {
                logLikeRiMi[i * numEdges * 2 + r2 * numEdges + m2] -= sumLogLike;
            }
        }
    }
}

// compute the log-probability matrix for character A,C,G,T for each position of each tip (i.e. haplotype)
void PhyloRead::computeCharMatHaplo() {
    int i,j,k,m,r;
    // int max_m, max_r;
    // double max_logL;
    double curr_logL;
    double sum_logL;
    int ch, curr_char;
    vector<char> mutLocTipChars;
    // mutLocTipChars[i * numTips + j] = the char of j-th haplotype if the mutation location = i
    vector<int> tipChars;
    bool isFirst;
    
    // initialization
    for (i=0; i<numTips; i++)
        tipChars.push_back(0);
    mutLocTipChars.clear();
    for (m=0; m<numEdges; m++) {
        getChars(m, tipChars);
        for (i=0; i<numTips; i++) {
            mutLocTipChars.push_back(tipChars[i]);
        }
    }
    memset(isValidCharMat, 0, numPos * numTips * 4 * sizeof(char));
    
    for (i=0; i<numPos; i++) {
        // cout << "i=" << i << endl << flush;
        if (isInVar[i]) {
            curr_char = charOrder[i * 4];
            for (j=0; j<numTips; j++) {
                charMatHaplo[i * numTips * 4 + j * 4 + curr_char] = 0.0;
                isValidCharMat[i * numTips * 4 + j * 4 + curr_char] = 1;
            }
        } else {
            for (r=0; r<2; r++) {
                for (m=0; m<numEdges; m++) {
                    curr_logL = logLikeRiMi[i * numEdges * 2 + r * numEdges + m];
                    for (j=0; j<numTips; j++) {
                        ch = (r + mutLocTipChars[m * numTips + j]) % 2;
                        curr_char = charOrder[i * 4 + ch];
                        if (!isValidCharMat[i * numTips * 4 + j * 4 + curr_char]) {
                            charMatHaplo[i * numTips * 4 + j * 4 + curr_char] = curr_logL;
                            isValidCharMat[i * numTips * 4 + j * 4 + curr_char] = 1;
                        } else {
                            charMatHaplo[i * numTips * 4 + j * 4 + curr_char] = log_x_plus_y(charMatHaplo[i * numTips * 4 + j * 4 + curr_char], curr_logL);
                        }
                    }
                }
            }
            // normalize
            for (j=0; j<numTips; j++) {
                isFirst = true;
                for (k=0; k<4; k++) {
                    if (isValidCharMat[i * numTips * 4 + j * 4 + k]) {
                        if (isFirst) {
                            sum_logL = charMatHaplo[i * numTips * 4 + j * 4 + k];
                            isFirst = false;
                        } else {
                            sum_logL = log_x_plus_y(sum_logL, charMatHaplo[i * numTips * 4 + j * 4 + k]);
                        }
                    }
                }
                for (k=0; k<4; k++) {
                    if (isValidCharMat[i * numTips * 4 + j * 4 + k]) {
                        charMatHaplo[i * numTips * 4 + j * 4 + k] -= sum_logL;
                    }
                }
            }
        }
    }
}

// report the set of characters with maximum probability for each haplotype
void PhyloRead::reportMaxCharSet(vector<string>& haplos) {
    int i,j,k;
    double opt_logProb;
    int opt_char;
    bool isFirst;
    char int2nucl[5] = {'A','C','G','T','N'};
    
    // initialization
    haplos.clear();
    for (j=0; j<numTips; j++) {
        haplos.push_back("");
    }

    for (i=0; i<numPos; i++) {
        for (j=0; j<numTips; j++) {
            isFirst = true;
            opt_char = 4; // N
            for (k=0; k<4; k++) {
                if (isValidCharMat[i * numTips * 4 + j * 4 + k]) {
                    if (isFirst) {
                        opt_logProb = charMatHaplo[i * numTips * 4 + j * 4 + k];
                        opt_char = k;
                        isFirst = false;
                    } else if (charMatHaplo[i * numTips * 4 + j * 4 + k] > opt_logProb) {
                        opt_logProb = charMatHaplo[i * numTips * 4 + j * 4 + k];
                        opt_char = k;
                    }
                }
            }
            haplos[j].append(1, int2nucl[opt_char]);
        }
    }
}

#define LOG_RATIO_THRES 1.0
// report the set of characters for each haplotype
void PhyloRead::reportHaploSeqs(vector<string>& haplos) {

    double* probs;
    // array for storing the probabilities for each position of each haplotype based on the reads
    // dimension: numTips * numPos * 4
    // probs[i*numPos*4 + j*4 + k] = probability for position j of haplotype i with nucleotide k based on the reads
    
    double* readLogProbs;
    // dimension: numTips
    // readLogProbs[i] = log-prob(current read | haplotype i)
    
    // tip considered
    char* tipConsidered;
    
    double* allReadLogProbs;
    // dimension: numReadPairs * numTips
    // allReadLogProbs[i*numTips + j] = log-prob(read pair i | haplotype j)
    
    char* allTipConsidered;
    // dimension: numReadPairs * numTips
    // allTipConsidered[i*numTips + j] = 1 if the read pair i is considered for haplotype j, 0 otherwise
    
    // integer to nucletide
    char int2nucl[6] = {'-','A','C','G','T','N'};
    
    Read* rd;
    int i,j,k,nucl,p,p_g;
    double currLogProb, currProb;
    double maxLogProb;
    double sumLogProb;
    double noErrTranLogProb, errTranLogProb;
    double maxProb;
    int maxNucl;
    bool isFirst;
    
    // initialize the array
    readLogProbs = new double[numTips];
    tipConsidered = new char[numTips];
    
    // for each read R = r1 r2 ... rl, compute the loglike(read R | haplotype k)
    // loglike(read R | haplotype k) = sum_i of log-prob(nucl = ri | haplotype k)

    noErrTranLogProb = log(1.0 - seqErr);
    errTranLogProb = log((1.0 - seqErr) / 3.0);

    // ----------------------------
    // only consider pair-end reads
    // ----------------------------
    
    allReadLogProbs = new double[(int)rds.pairReads.size() * numTips];
    // dimension: numReadPairs * numTips
    // allReadLogProbs[i*numTips + j] = log-prob(read pair i | haplotype j)
    
    allTipConsidered = new char[(int)rds.pairReads.size() * numTips];
    memset(allTipConsidered, 0, (int)rds.pairReads.size() * numTips * sizeof(char));
    // dimension: numReadPairs * numTips
    // allTipConsidered[i*numTips + j] = 1 if the read pair i is considered for haplotype j, 0 otherwise

    
    for (i=0; i<(int)rds.pairReads.size(); i++) {
        // cout << "i==" << i << endl << flush;
        memset(readLogProbs, 0, numTips * sizeof(double));
        for (j=0; j<2; j++) {
            // cout << "j==" << j << endl << flush;
            if (j==0) {
                // consider the first read
                rd = rds.pairReads[i].first;
            } else {
                // consider the second read
                rd = rds.pairReads[i].second;
            }
            for (p=0; p<rd->seq.length(); p++) {
                p_g = rd->mapPos + p;
                if (isDiscard[p_g])
                    continue;
                for (k=0; k<numTips; k++) {
                    isFirst = true;
                    for (nucl=0; nucl<4; nucl++) {
                        if (isValidCharMat[p_g * numTips * 4 + k * 4 + nucl]) {
                            if (isFirst) {
                                if (nucl == nucl2int[(int)rd->seq[p]]-1) {
                                    currLogProb = charMatHaplo[p_g * numTips * 4 + k * 4 + nucl] + noErrTranLogProb;
                                } else {
                                    currLogProb = charMatHaplo[p_g * numTips * 4 + k * 4 + nucl] + errTranLogProb;
                                }
                                isFirst = false;
                            } else {
                                if (nucl == nucl2int[(int)rd->seq[p]]-1) {
                                    currLogProb = log_x_plus_y(currLogProb, charMatHaplo[p_g * numTips * 4 + k * 4 + nucl] + noErrTranLogProb);
                                } else {
                                    currLogProb = log_x_plus_y(currLogProb, charMatHaplo[p_g * numTips * 4 + k * 4 + nucl] + errTranLogProb);
                                }
                            }
                        }
                    }
                    readLogProbs[k] += currLogProb;
                }
            }
        }
        // get the max log prob
        maxLogProb = readLogProbs[0];
        for (k=1; k<numTips; k++) {
            if (readLogProbs[k] > maxLogProb) {
                maxLogProb = readLogProbs[k];
            }
        }
        // check which tip should be considered
        memset(tipConsidered, 0, numTips * sizeof(char));
        for (k=0; k<numTips; k++) {
            if (maxLogProb - readLogProbs[k] <= LOG_RATIO_THRES)
                tipConsidered[k] = 1;
        }
        // normalized the probabilities
        isFirst = true;
        for (k=0; k<numTips; k++) {
            if (tipConsidered[k]) {
                if (isFirst) {
                    sumLogProb = readLogProbs[k];
                    isFirst = false;
                } else {
                    sumLogProb = log_x_plus_y(sumLogProb, readLogProbs[k]);
                }
            }
        }
        for (k=0; k<numTips; k++) {
            if (tipConsidered[k]) {
                readLogProbs[k] -= sumLogProb;
                allTipConsidered[i*numTips + k] = 1;
                allReadLogProbs[i*numTips + k] = readLogProbs[k];
            }
        }
    }

    // compute the consensus
    // before do so, need to get back all the inserts
    rds.updateForInserts();
    numPos = rds.lastMapPos + 1;
    probs = new double[numTips * numPos * 5];
    memset(probs, 0, numTips * numPos * 5 * sizeof(double));

    for (i=0; i<(int)rds.pairReads.size(); i++) {
        for (k=0; k<numTips; k++) {
            if (allTipConsidered[i*numTips + k]) {
                currProb = exp(allReadLogProbs[i*numTips + k]);
                for (j=0; j<2; j++) {
                    // cout << "j==" << j << endl << flush;
                    if (j==0) {
                        // consider the first read
                        rd = rds.pairReads[i].first;
                    } else {
                        // consider the second read
                        rd = rds.pairReads[i].second;
                    }
                    for (p=0; p<rd->seq.length(); p++) {
                        p_g = rd->mapPos + p;
                        probs[k*numPos*5 + p_g*5 + nucl2int[(int)rd->seq[p]]] += currProb;
                    }
                }
            }
        }
    }
    
    // report the haplotypes
    haplos.clear();
    for (k=0; k<numTips; k++) {
        haplos.push_back("");
    }
    for (k=0; k<numTips; k++) {
        for (p=0; p<numPos; p++) {
            maxProb = 0.0;
            maxNucl = 5;
            for (j=0; j<5; j++) {
                if (probs[k*numPos*5 + p*5 + j] > maxProb) {
                    maxProb = probs[k*numPos*5 + p*5 + j];
                    maxNucl = j;
                }
            }
            haplos[k].append(1, int2nucl[maxNucl]);
        }
    }
    for (k=0; k<numTips; k++) {
        // remove the gaps
        i=0;
        for (p=0; p<numPos; p++) {
            if (haplos[k].at(p) != '-') {
                if (p > i) {
                    haplos[k].at(i) = haplos[k].at(p);
                }
                i++;
            }
        }
        if (i > 0 && i < numPos)
            haplos[k] = haplos[k].substr(0,i);
        else if (i==0)
            haplos[k] = "";
        // remove the N's in the front or at the end
        i=0;
        while (i<haplos[k].length() && haplos[k].at(i)=='N')
            i++;
        p=haplos[k].length()-1;
        while (p>=0 && haplos[k].at(p)=='N')
            p--;
        if (i <= p) {
            haplos[k] = haplos[k].substr(i,p-i+1);
        }
            
    }
    delete[] probs;
    delete[] readLogProbs;
    delete[] tipConsidered;
    delete[] allReadLogProbs;
    delete[] allTipConsidered;
}

