//
//  constructHaplo.cpp
//  HaploCount
//
//  Created by Thomas Wong on 10/6/19.
//

#include "constructHaplo.h"

void Haplotypes::loadMaxTreeFile(char* maxTreeFile) {
    ifstream fin;
    string strTree;
    string strErr;
    bool isTreeStr;
    UnlabelTree maxTree;
    double maxErr;
    vector<int> topInt;
    vector<double> haploFreq;
    vector<double> edgeLen;
    
    // the prefix of the line with error rate
    string prefix = "# error rate: ";
    int prefix_len = prefix.length();
    
    int i;
    
    // read the max-tree-file
    // cout << "read the max-tree-file" << endl << flush;
    fin.open(maxTreeFile);
    isTreeStr = false;
    while ((!isTreeStr) && getline(fin,strTree)) {
        if (strTree.length() == 0)
            continue;
        // get the corresponding error rate
        if (strTree.length() > prefix_len && strTree.substr(0,prefix_len) == prefix) {
            strErr = strTree.substr(prefix_len);
            cout << "Error rate: " << strErr << endl;
            maxErr = atof(strErr.c_str());
        }
        // read the tree
        else if (strTree[0] != '#') {
            isTreeStr = true;
            maxTree.topTxt2Int(strTree, topInt, haploFreq, edgeLen);
        }
    }
    fin.close();
    if (!isTreeStr) {
        cerr << "Error! The file: " << maxTreeFile << " cannot be read" << endl;
        exit(1);
    }
    
    // save the topology
    // cout << "save the topology" << endl << flush;

    phyloRead.numTips = haploFreq.size();
    phyloRead.numNodes = 2 * phyloRead.numTips - 1;
    phyloRead.numInNodes = phyloRead.numTips - 1;
    phyloRead.numEdges = phyloRead.numNodes - 1;
    phyloRead.currTopology.clear();
    for (i=0; i<topInt.size(); i++)
        phyloRead.currTopology.push_back(topInt[i]);
    phyloRead.isDescendantUpdated = false;
    phyloRead.needUpdated = true;
    
    // cout << "buildIsDescendant()" << endl << flush;
    phyloRead.buildIsDescendant();

    // set up the variables
    // cout << "set up the variables" << endl << flush;
    phyloRead.setupVariables();
    
    // save all the frequencies
    // cout << "save all the frequencies" << endl << flush;
    for (i=0; i<phyloRead.numTips; i++)
        phyloRead.haploFreqs_t[i] = haploFreq[i];

    // save the sequencing error
    // cout << "save the sequencing error" << endl << flush;
    phyloRead.seqErr = maxErr;
}

// load samfile
void Haplotypes::load_sam(char* samFile) {
    
    // compute the data matrix
    set<int> t;
    phyloRead.loadData(samFile, t);
}

// construct the haplotypes based on the connectivity of the columns inside reads
void Haplotypes::construct(char* samFile, char* maxTreeFile, char* outFile) {
    
    double logL;
    vector<long double> edgeLen;
    string treeStr;
    int num_digits = 7;
    vector<string> haplos;
    
    cout << "samFile = " << samFile << endl << flush;
    cout << "maxTreeFile = " << maxTreeFile << endl << flush;
    
    cout << "loading the samFile..." << endl << flush;
    load_sam(samFile);

    cout << "loading the max tree file..." << endl << flush;
    loadMaxTreeFile(maxTreeFile);
    
    cout << "set up variables..." << endl << flush;
    phyloRead.setupVariables();

    cout << "normalize haplotype frequecies..." << endl << flush;
    phyloRead.normalizeHapFreqs();

    cout << "update weights..." << endl << flush;
    phyloRead.updateWeights();

    cout << "update expected log freqs..." << endl << flush;
    phyloRead.updateExpectLogFreqs();
    
    cout << "compute edge lengths..." << endl << flush;
    phyloRead.computeEdgeLen();
    
    cout << "compute the loglikelihood given mi and ri for all positions..." << endl;
    phyloRead.computeLogLikeRiMi();

    cout << "compute the log-probability matrix for each character on each position of each tip (i.e. haplotype)..." << endl << flush;
    phyloRead.computeCharMatHaplo();

    // cout << "compute the haplotypes with max likelihoods..." << endl << flush;
    // phyloRead.reportMaxCharSet(haplos);

    cout << "compute the haplotypes ..." << endl << flush;
    phyloRead.reportHaploSeqs(haplos);
    
    cout << "output the haplotypes to a file" << endl << flush;
    outHaplo(outFile, haplos, phyloRead.haploFreqs);
}

// output the haplotypes
void Haplotypes::outHaplo(char* outFile, vector<string>& haplos, vector<double>& freqs) {
    
    ofstream fout;
    int i;
    
    fout.open(outFile);
    for (i=0; i<haplos.size(); i++) {
        fout << ">haplo" << i+1 << "_" << freqs[i] << endl;
        fout << haplos[i] << endl;
    }
    fout.close();
}

