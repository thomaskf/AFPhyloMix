#ifndef _TREE_
#define _TREE_

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <set>
#include <map>
#include <vector>
#include <stack>
#include <string>
#include <algorithm>
#include "mylib.h"
#include "rawMatrix.h"
#include "myRand.h"

using namespace std;

class UnlabelTree {
public:

    // randomly generate a rooted tree
    void randRootedTree(int numTips, vector<int>& topology, MyRand& myrand);
    
	// generate all possible rooted trees
	// output format: topologies[i] = the i-th topology (in integer format)
	void genRootedTrees(int numTips, vector<vector<int> >& topologies);

	// generate all possible unrooted trees
	// output format: topologies[i] = the i-th topology (in integer format)
	void genUnrootedTrees(int numTips, vector<vector<int> >& topologies);

	// only for the number of haplotypes between 2 and 20
	// quickly retrieve all possible unrooted trees
	// output format: topologies[i] = the i-th topology (in integer format)
	void getUnrootedTrees(int numHap, vector<vector<int> >& topologies);

	// tidy up the topology format such that all the leaf nodes appear in the beginning
	void tidyUpFormat(vector<int>& topInt, vector<int>& updateTop);
	
	// convert the topology from integer format into text format
	string topInt2Txt(vector<int>& topInt);

    // convert the topology from integer format into text format
    // for labelled tree
    string topInt2TxtLabelTree(vector<int>& topInt);

	// convert the topology from integer format into text format
	// the topology begins with leaves
	string topInt2Txt(vector<int>& topInt, vector<double>& haploFreq, vector<double>& edgeLen, int numDigits);
    string topInt2Txt(vector<int>& topInt, vector<double>& haploFreq, vector<long double>& edgeLen, int numDigits);


	// convert the topology from integer format into text format
	string topInt2Txt(vector<int>& topInt, vector<int>& nodeOrder);
	
    // convert the toplogy from integer format into text format (for labeled tree)
    string topInt2TxtLabelTree(vector<int>& topInt, vector<int>& nodeOrder);

    // convert the toplogy from text format into integer format
    void topTxt2Int(string& topStr, vector<int>& topInt, vector<double>& haploFreq, vector<double>& edgeLen);
	
	// show the topologies
	void showTopologies(vector<vector<int> >& topologies);
	
	// get the tree depth
	int treeDepth(vector<int>& topInt);

	// get the tree depth
	int treeDepth(vector<int>& topInt, vector<int>& nodeOrder);
	
	// change the root position to the edge connecting to the node x (i.e. newRoot)
	void changeRoot(vector<int>& topInt, int newRoot, vector<int>& newTop, vector<int>& newEdge, vector<int>& newNodeOrder);
	
	// collect all topologies (in text version) when rooted at leaf nodes
	// void collectLeafRootedTopologies(vector<int>& topInt, set<string>& allTopologies);

	// collect the rooted topologies with minimum tree height
	// void collectMinHeightTopologies(vector<int>& topInt, set<string>& minHeightTopologies, int& tree_depth);

	// obtain a reprsentative for the rooted topologies
	// (the rooted topologies with min height height and smallest in string)
	string getRepresent(vector<int>& topInt, vector<int>& representTop);

	// update the topology so that the node order becomes 0,1,2,3,....
	void updateTopology(vector<int>& topInt, vector<int>& nodeOrder, vector<int>& newTop);
	
	// show the topology
	void showTopology(vector<int>& topInt);

	// show the topology
	void showTopology(vector<int>& topInt, vector<int>& newNodeOrder);

	// show the topology
	void showTopology(vector<int>& topInt, vector<int>& newEdge, vector<int>& newNodeOrder);

	// show the topology in one line
	void showTopologyInOneLine(vector<int>& topInt);
    
    // generate K possible labelled topologies based on an input topology
    // void genKLabelTopoloy(int numTips, vector<int>& topInt, int K, vector<vector<int> >& newTopInts);
    
    // get the sister node
    int getSisterNode(int node1, vector<int>& topInt);
    
    // NNI
    // update: topology is always changed
    bool doNNI(vector<int>& topInt, int numTips, MyRand& myrand, vector<pair<int,int> >& changeNodeIDs);
    
    // subtree swapping
    // isDescendant[i*numNodes+j] = true if node j is descendant of node i (relationship between two nodes)
    // return true if the topology is changed
    bool swapSubTree(vector<int>& topInt, int numTips, int numNodes, MyRand& myrand, vector<bool>& isDescendant, vector<pair<int,int> >& changeNodeIDs);
    
    // tidy up the messy topology
    void tidyUpMessyTopology(int numTips, vector<int>& topIntMessy, vector<pair<int,int> >& changeNodeIDs);
    
private:
    void subrandRootedTree(int numTips, vector<int>& topology, MyRand& myrand, int& rootID);
};

// the leaf group with the same number of subgroups
class LeafGrpSameNum {
public:
    vector<vector<int> > leafGrps;
    void clear();
    int size();
    void show();
};

// all the leaf groups
class LeafGrps {
public:
    vector<LeafGrpSameNum> allLeafGrps;
    
    void addLeafGrp(vector<int>& leafGrp, int grp_num);
    void clear();
    void showSummary();
    void show();
};

class Node {
public:
    Node* parent;
    vector<Node*> children;
    int value;
    int label;
    
    // constructor
    Node(int v);
    Node(int v, Node* leftChild, Node* rightChild);

    // set the labels of the node and its children
    void setLabel(int parentLabel, int& lastLabel, vector<bool>& mutationOccur);
};

class RootTree {
public:

    // root
    Node* rootNode;
    
    // list of nodes (including leaves)
    vector<Node*> nodeList;
    
    // list of Leaves
    vector<Node*> leaveList;
    
    // destructor
    ~RootTree();
    
    // load the topology
    void loadTopology(vector<int>& topInt);
    
    // load the topology from a newick tree file
    // and then obtain the mutation sets
    void get_mutation_sets(char* treeFile, vector<vector<int> >& c_sets, int& numHaps);
    
    // clear
    void clear();
    
    // print the topology
    void printTopology();
    
    // get all possible groups of the leaves due to the mutations on different edges
    void getAllPossibleLeafGroups(LeafGrps& leafGrps);

private:
    // return the labels of the leaves given the positions of the mutations
    void getLabels(vector<bool>& mutationOccur, vector<int>& leaveLabels, string& s, int& numGrps);

    // next possible set of mutation positions
    bool nextSet(vector<bool>& mutationOccur);
};

#endif
