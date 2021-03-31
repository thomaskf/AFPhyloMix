#ifndef _GTREE_
#define _GTREE_

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <set>
#include <algorithm>
#include <vector>
#include <string>
#include <math.h>
#include "mylib.h"
#include "definitions.h"

using namespace std;
class gNode {
public:
    vector<gNode*> children;
    vector<double> edgeLens;
    vector<set<string> > childValues;
    string name;
    string value;
    double freq;
    bool isLeaf;
    int id;
    int parent_id;
    bool consider;
    
    void addChild(gNode* child);
    
    void getChildValues(gNode* parent, set<string>& values);
    
    string getPattern(gNode* parent);

    string getNodeStr(gNode* parent);
    
    void replaceChild(gNode* newChild, gNode* childToReplace, double newEdgeLen);
    
	// compute the sum_of_square difference between two nodes
	bool SSDiff(gNode* aNode, gNode* thisParent, gNode* aParent, double& ssdiff, double diffThres, string& algn);
};

// rooted or unrooted tree
class gTree {
public:

    // root / an internal node
    gNode* root;
    
    // leave nodes
    vector<gNode*> leaves;
    
    // all nodes
    vector<gNode*> allNodes;
    
    // edges
    vector<pair<gNode*,gNode*> > edges;
    
    // destructor
    ~gTree();

	// load the topology file
	void loadTopFile(char* fileName);
    
    // load the tree sequence
    void loadTreeStr(string treeStr);
    
    // load the topology matrix
    void loadTopMatrix(vector<int>& topInt, vector<double>& haploFreq, vector<double>& edgeLen);
    
    // clear
    void clear();
    
    // print the topology
    void printTopology();
    
    // get tree in string format
    string getTreeStr();
    
    // reset all children values
    void resetAllChildValues();
    
    // check whether more than one mutation
    bool isMoreThanOneMutate();
    
    // get the pattern
    string getPattern();
    
    // change root to a terminal edge connecting to a leaf node
    void changeRoot(gNode* leafNode);
    
	// RMS difference between another rooted tree
	// report false if another rooted tree and this rooted tree are not in the same topology
	bool RMSDiff(gTree* aTree, double& diff, double diffThres, string& algn);

// supporting function

    int getCloseBracketPos(string& str, int openBracPos);

    // create a node for str[startPos ... endPos]
    // this node has to be enclosed by a pair of brackets
    // i.e. str[startPos] == '(' and str[endPos] == ')'
    gNode* createInterNode(string& str, int startPos, int endPos);
    
	// assign the IDs of the nodes
	void assignNodeID();
	
	// get all the nodes
	void getAllNodes();
    
    // remove an internal node from the tree
    // assume no children array containing that node
    void rmInterNode(gNode* nodeToRemove);
    
    // remove an edge
    void rmEdge(gNode* n1, gNode* n2);
    
    // get all the frequencies
    void getFreqs();

    // create an internal node from a single roow in the topology matrix
    gNode* createNodeFrRow(int rowID, vector<int>& topInt, vector<double>& haploFreq, vector<double>& edgeLen);
};

#endif
