#include "tree.h"

// randomly generate a rooted tree
// numTips <= 1
void UnlabelTree::randRootedTree(int numTips, vector<int>& topology, MyRand& myrand) {
    int rootID;
    vector<int> top;
    subrandRootedTree(numTips, top, myrand, rootID);
    topology.clear();
    tidyUpFormat(top, topology);;
}

void UnlabelTree::subrandRootedTree(int numTips, vector<int>& topology, MyRand& myrand, int& rootID) {
    int leftRootID, rightRootID;
    int i;
    if (numTips==1) {
        topology.push_back(-1);
        topology.push_back(-1);
    } else {
        i = myrand.iunif(1, numTips-1);
        // left child has i tips and right child has numTips-i tips
        subrandRootedTree(i, topology, myrand, leftRootID);
        subrandRootedTree(numTips-i, topology, myrand, rightRootID);
        topology.push_back(leftRootID);
        topology.push_back(rightRootID);
    }
    rootID = (topology.size() - 2) / 2;
}

// generate all possible rooted trees
// output format: topologies[i] = the i-th topology (in integer format)
void UnlabelTree::genRootedTrees(int numTips, vector<vector<int> >& topologies) {
	int i,j,k,l,s,m,p,q;
	int numTipLeftTree, numTipRightTree;
	
	vector<vector<vector<int> > > treeList;
	vector<int> topology;
	vector<int> tidyTop;
	treeList.resize(numTips);
	// initalize the first two cases
	// tree with zero tip, which is empty
	treeList[0].clear();
	// tree with one tip
	treeList[1].resize(1);
	treeList[1].at(0).clear();
	treeList[1].at(0).push_back(-1);
	treeList[1].at(0).push_back(-1);
	for (i=2; i<numTips; i++) {
		m=0;
		for (j=1; j<=i/2; j++) {
			numTipLeftTree = 2 * j - 1;
			numTipRightTree = 2 * (i-j) - 1;
			for (k=0; k<treeList[j].size(); k++) {
				if (i-j==j)
					s = k;
				else
					s = 0;
				for (l=s; l<treeList[i-j].size(); l++) {
					topology.clear();
					// insert the left part of the tree
					topology.insert(topology.begin(), treeList[j].at(k).begin(), treeList[j].at(k).end());
					// insert the right part of the tree
					for (p=0; p<treeList[i-j].at(l).size(); p++) {
						q = treeList[i-j].at(l).at(p);
						if (q == -1)
							topology.push_back(q);
						else
							topology.push_back(q + numTipLeftTree);
					}
					// insert the root
					topology.push_back(numTipLeftTree-1);
					topology.push_back(numTipLeftTree+numTipRightTree-1);
					// put into treeList[i]
					treeList[i].resize(m+1);
					treeList[i].at(m).clear();
					treeList[i].at(m).insert(treeList[i].at(m).begin(),topology.begin(),topology.end());
					m++;
				}
			}
		}
		cerr << i << " tips: " << treeList[i].size() << " rooted unlabeled trees" << endl;
	}
	// for i = numTips
	i = numTips;
	m=0;
	topologies.clear();
	for (j=1; j<=i/2; j++) {
		numTipLeftTree = 2 * j - 1;
		numTipRightTree = 2 * (i-j) - 1;
		for (k=0; k<treeList[j].size(); k++) {
			if (i-j==j)
				s = k;
			else
				s = 0;
			for (l=s; l<treeList[i-j].size(); l++) {
				topology.clear();
				// insert the left part of the tree
				topology.insert(topology.begin(), treeList[j].at(k).begin(), treeList[j].at(k).end());
				// insert the right part of the tree
				for (p=0; p<treeList[i-j].at(l).size(); p++) {
					q = treeList[i-j].at(l).at(p);
					if (q == -1)
						topology.push_back(q);
					else
						topology.push_back(q + numTipLeftTree);
				}
				// insert the root
				topology.push_back(numTipLeftTree-1);
				topology.push_back(numTipLeftTree+numTipRightTree-1);
				// tidy up the topology
				tidyUpFormat(topology, tidyTop);
				// put into treeList[i]
				topologies.resize(m+1);
				topologies.at(m).clear();
				topologies.at(m).insert(topologies.at(m).begin(),tidyTop.begin(),tidyTop.end());
				m++;
			}
		}
	}
	cerr << i << " tips: " << topologies.size() << " unlabeled trees" << endl;
}

// generate all possible unrooted trees
// output format: topologies[i] = the i-th topology (in integer format)
void UnlabelTree::genUnrootedTrees(int numTips, vector<vector<int> >& topologies) {
	int i,j,k,l,s,m,p,q;
	int numTipLeftTree, numTipRightTree;
	set<string> unroot_trees;
	set<string>::iterator itr;
	vector<vector<vector<int> > > treeList;
	vector<int> topology;
	vector<int> tidyTop;
	string treeStr;
	treeList.resize(numTips);
	// initalize the first two cases
	// tree with zero tip, which is empty
	treeList[0].clear();
	// tree with one tip
	treeList[1].resize(1);
	treeList[1].at(0).clear();
	treeList[1].at(0).push_back(-1);
	treeList[1].at(0).push_back(-1);
	for (i=2; i<=numTips; i++) {
		m=0;
		for (j=1; j<=i/2; j++) {
			numTipLeftTree = 2 * j - 1;
			numTipRightTree = 2 * (i-j) - 1;
			for (k=0; k<treeList[j].size(); k++) {
				if (i-j==j)
					s = k;
				else
					s = 0;
				for (l=s; l<treeList[i-j].size(); l++) {
					topology.clear();
					// insert the left part of the tree
					topology.insert(topology.begin(), treeList[j].at(k).begin(), treeList[j].at(k).end());
					// insert the right part of the tree
					for (p=0; p<treeList[i-j].at(l).size(); p++) {
						q = treeList[i-j].at(l).at(p);
						if (q == -1)
							topology.push_back(q);
						else
							topology.push_back(q + numTipLeftTree);
					}
					// insert the root
					topology.push_back(numTipLeftTree-1);
					topology.push_back(numTipLeftTree+numTipRightTree-1);
					// put into treeList[i]
					if (i < numTips) {
						treeList[i].resize(m+1);
						treeList[i].at(m).clear();
						treeList[i].at(m).insert(treeList[i].at(m).begin(),topology.begin(),topology.end());
						m++;
					} else {
						tidyUpFormat(topology, tidyTop);
						treeStr = getRepresent(tidyTop,topology);
						itr = unroot_trees.find(treeStr);
						if (itr == unroot_trees.end()) {
							unroot_trees.insert(treeStr);
							
							topologies.resize(m+1);
							topologies.at(m).clear();
							topologies.at(m).insert(topologies.at(m).begin(),tidyTop.begin(),tidyTop.end());
							m++;
							
						}
					}
				}
			}
		}
	}
}

// only for the number of haplotypes between 2 and 20
// quickly retrieve all possible unrooted trees
// output format: topologies[i] = the i-th topology (in integer format)
void UnlabelTree::getUnrootedTrees(int numHap, vector<vector<int> >& topologies) {
	int numMatrix[] = {0, 1, 1, 1, 1, 1, 2, 2, 4, 6, 11, 18, 37, 66, 135, 265, 552, 1132, 2410, 5098, 11020};
	int i, j, k;
	if (numHap < 2 || numHap > 20) {
		cerr << "Number of haplotypes has to be between 2 and 20" << endl;
		exit(1);
	}
	topologies.resize(numMatrix[numHap]);
	k=0;
	for (i=2; i<numHap; i++) {
		k += numMatrix[i] * (i-1) * 2;
	}
	for (i=0; i<numMatrix[numHap]; i++) {
		topologies[i].clear();
		for (j=0; j<2*numHap; j++) {
			topologies[i].push_back(-1);
		}
		for (j=0; j<(numHap-1)*2; j++) {
			topologies[i].push_back(rawMatrix[k++]);
		}
	}
}

// tidy up the topology format such that all the leaf nodes appear in the beginning
void UnlabelTree::tidyUpFormat(vector<int>& topInt, vector<int>& updateTop) {
	vector<int> newID;
	int nodeNum;
	int i,k;
	
	k=0;
	nodeNum = ((int)topInt.size()) / 2;
	newID.resize(nodeNum);
	updateTop.clear();
	// assigning a new ID for all the leaves
	for (i=0; i<nodeNum; i++) {
		if (topInt[2*i]==-1 || topInt[2*i+1]==-1) {
			// node i is a leaf
			// its new ID is k
			newID[i] = k;
			k++;
		}
	}
	for (i=0; i<k; i++) {
		updateTop.push_back(-1);
		updateTop.push_back(-1);
	}
	// now consider the internal nodes
	for (i=0; i<nodeNum; i++) {
		if (!(topInt[2*i]==-1 || topInt[2*i+1]==-1)) {
			// an internal node
			newID[i] = k;
			k++;
			updateTop.push_back(newID[topInt[2*i]]);
			updateTop.push_back(newID[topInt[2*i+1]]);
		}
	}
}

// convert the toplogy from integer format into text format
string UnlabelTree::topInt2Txt(vector<int>& topInt) {
	vector<string> nodeStr;
	int leftID, rightID;
	string leftTxt, rightTxt;
	int k;
	
	nodeStr.clear();
	for (k=0; k<topInt.size()/2; k++) {
		leftID = topInt[2*k];
		rightID = topInt[2*k+1];
		if (leftID==-1 || rightID==-1) {
			// a leaf
			nodeStr.push_back("x");
		} else {
			// an internal node
			leftTxt = nodeStr[leftID];
			rightTxt = nodeStr[rightID];
			if (leftTxt <= rightTxt)
				nodeStr.push_back("(" + leftTxt + "," + rightTxt + ")");
			else
				nodeStr.push_back("(" + rightTxt + "," + leftTxt + ")");
		}
	}
	return nodeStr[nodeStr.size()-1];
}

// convert the toplogy from integer format into text format
// for labelled tree
string UnlabelTree::topInt2TxtLabelTree(vector<int>& topInt) {
    vector<string> nodeStr;
    int leftID, rightID;
    string leftTxt, rightTxt;
    int k;
    
    nodeStr.clear();
    for (k=0; k<topInt.size()/2; k++) {
        leftID = topInt[2*k];
        rightID = topInt[2*k+1];
        if (leftID==-1 || rightID==-1) {
            // a leaf
            nodeStr.push_back(int2str(k));
        } else {
            // an internal node
            leftTxt = nodeStr[leftID];
            rightTxt = nodeStr[rightID];
            if (leftTxt <= rightTxt)
                nodeStr.push_back("(" + leftTxt + "," + rightTxt + ")");
            else
                nodeStr.push_back("(" + rightTxt + "," + leftTxt + ")");
        }
    }
    return nodeStr[nodeStr.size()-1];
}

// convert the toplogy from integer format into text format
// the topology begins with leaves
string UnlabelTree::topInt2Txt(vector<int>& topInt, vector<double>& haploFreq, vector<double>& edgeLen, int numDigits) {
	vector<string> nodeStr;
	int leftID, rightID;
	string leftTxt, rightTxt;
	int k;
    double leftLen, rightLen;
    string s;
	
    // cout << "[enter UnlabelTree::topInt2Txt]" << endl << flush;
	nodeStr.clear();
	for (k=0; k<topInt.size()/2; k++) {
		leftID = topInt[2*k];
		rightID = topInt[2*k+1];
		if (leftID==-1 || rightID==-1) {
			// a leaf
			nodeStr.push_back("f=" + doublToStr(haploFreq[k], numDigits));
		} else {
			// an internal node
			leftTxt = nodeStr[leftID];
			rightTxt = nodeStr[rightID];
            
            if (leftID >= edgeLen.size())
                leftLen = 0.0;
            else
                leftLen = edgeLen[leftID];
            
            if (rightID >= edgeLen.size())
                rightLen = 0.0;
            else
                rightLen = edgeLen[rightID];

            if (leftTxt <= rightTxt) {
				nodeStr.push_back("(" + leftTxt + ":" + doublToStr(leftLen, numDigits) + "," + rightTxt + ":" + doublToStr(rightLen, numDigits) + ")");
            } else {
				nodeStr.push_back("(" + rightTxt + ":" + doublToStr(rightLen, numDigits) + "," + leftTxt + ":" + doublToStr(leftLen, numDigits) + ")");
            }
		}
	}
    s = nodeStr[nodeStr.size()-1];
    // cout << "[leave UnlabelTree::topInt2Txt]" << endl << flush;
    return s;
}

// convert the toplogy from integer format into text format
// the topology begins with leaves
string UnlabelTree::topInt2Txt(vector<int>& topInt, vector<double>& haploFreq, vector<long double>& edgeLen, int numDigits) {
    vector<string> nodeStr;
    int leftID, rightID;
    string leftTxt, rightTxt;
    int k;
    double leftLen, rightLen;
    string s;
    
    nodeStr.clear();
    for (k=0; k<topInt.size()/2; k++) {
        leftID = topInt[2*k];
        rightID = topInt[2*k+1];
        if (leftID==-1 || rightID==-1) {
            // a leaf
            nodeStr.push_back("f=" + doublToStr(haploFreq[k], numDigits));
        } else {
            // an internal node
            leftTxt = nodeStr[leftID];
            rightTxt = nodeStr[rightID];
            
            if (leftID >= edgeLen.size())
                leftLen = 0.0;
            else
                leftLen = edgeLen[leftID];
            
            if (rightID >= edgeLen.size())
                rightLen = 0.0;
            else
                rightLen = edgeLen[rightID];
            
            if (leftTxt <= rightTxt) {
                nodeStr.push_back("(" + leftTxt + ":" + doublToStr(leftLen, numDigits) + "," + rightTxt + ":" + doublToStr(rightLen, numDigits) + ")");
            } else {
                nodeStr.push_back("(" + rightTxt + ":" + doublToStr(rightLen, numDigits) + "," + leftTxt + ":" + doublToStr(leftLen, numDigits) + ")");
            }
        }
    }
    s = nodeStr[nodeStr.size()-1];
    return s;
}


// convert the toplogy from integer format into text format
string UnlabelTree::topInt2Txt(vector<int>& topInt, vector<int>& nodeOrder) {
	vector<string> nodeStr;
	int leftID, rightID;
	string leftTxt, rightTxt;
	int numNodes;
	int i,k;
	
	numNodes = ((int)topInt.size()) / 2;
	nodeStr.resize(numNodes);
	for (i=0; i<numNodes; i++) {
		k = nodeOrder[i];
		leftID = topInt[2*k];
		rightID = topInt[2*k+1];
		if (leftID==-1 || rightID==-1) {
			// a leaf
			nodeStr[k]="x";
		} else {
			// an internal node
			leftTxt = nodeStr[leftID];
			rightTxt = nodeStr[rightID];
			if (leftTxt <= rightTxt)
				nodeStr[k]="(" + leftTxt + "," + rightTxt + ")";
			else
				nodeStr[k]="(" + rightTxt + "," + leftTxt + ")";
		}
	}
	return nodeStr[nodeStr.size()-1];
}

// convert the toplogy from integer format into text format (for labeled tree)
string UnlabelTree::topInt2TxtLabelTree(vector<int>& topInt, vector<int>& tipOrder) {
    vector<string> nodeStr;
    int leftID, rightID;
    string leftTxt, rightTxt;
    int numNodes;
    int i;
    
    numNodes = ((int)topInt.size()) / 2;
    nodeStr.resize(numNodes);
    for (i=0; i<numNodes; i++) {
        leftID = topInt[2*i];
        rightID = topInt[2*i+1];
        if (leftID==-1 || rightID==-1) {
            // a leaf
            // the array begins with the leaves
            nodeStr[i]=int2str(tipOrder[i]);
        } else {
            // an internal node
            leftTxt = nodeStr[leftID];
            rightTxt = nodeStr[rightID];
            if (leftTxt <= rightTxt)
                nodeStr[i]="(" + leftTxt + "," + rightTxt + ")";
            else
                nodeStr[i]="(" + rightTxt + "," + leftTxt + ")";
        }
    }
    return nodeStr[nodeStr.size()-1];
}

// get the corresponding pairs of starting and ending positions from topology sequence
// assume no space in the topStr
void getPosPair(string& topStr, vector<pair<int,int> >& ipos, vector<int>& iparents, vector<pair<int,int> >& ichildren) {
    int p; // position on the sequence
    int k; // k-th bracket
    int t,n;
    string tStr;
    ipos.clear();
    iparents.clear();
    ichildren.clear();
    
    vector<pair<int,int> > tmpPos;
    vector<int> tmpParents;
    vector<bool> isLeaf;
    vector<int> changeTo;
    vector<int> changeFr;
    vector<pair<int,int> > tpos;
    vector<int> tparents;
    vector<pair<int,int> > tchildren;
    
    k=-1;
    for (p=0; p<(int)topStr.length(); p++) {
        if (topStr[p]=='(' || (p-1>=0 && topStr[p-1]=='(') || (p-1>=0 && topStr[p-1]==',')) {
            k = (int)tmpPos.size();
            tmpPos.push_back(pair<int,int>(p,-1));
            tmpParents.push_back(-1);
        } else if ((p+1<topStr.length() && topStr[p+1] == ')')|| (p+1<topStr.length() && topStr[p+1] == ',')) {
            if (k==-1) {
                cerr << "[getBracketPair] Error! The parentheses are not balanced" << endl;
                exit(1);
            }
            tmpPos[k].second = p;
            t = k;
            k--;
            while (k>=0 && tmpPos[k].second != -1)
                k--;
            tmpParents[t] = k;
        }
    }
    if (k==0) {
        tmpPos[k].second = topStr.length();
        k--;
    }
    if (k!=-1) {
        cerr << "[getBracketPair] Error! The parentheses are not balanced" << endl;
        exit(1);
    }
    
    // output the information reversely
    n = (int)tmpPos.size();
    t = n - 1;
    for (k=t; k>=0; k--) {
        tpos.push_back(pair<int,int>(tmpPos[k].first,tmpPos[k].second));
        if (tmpParents[k]==-1)
            tparents.push_back(-1);
        else
            tparents.push_back(t - tmpParents[k]);
    }
    // set the children
    for (k=0; k<n; k++) {
        tchildren.push_back(pair<int,int>(-1,-1));
    }
    for (k=0; k<n; k++) {
        t = tparents[k];
        if (t != -1) {
            if (tchildren[t].first == -1)
                tchildren[t].first = k;
            else if (tchildren[t].second == -1)
                tchildren[t].second = k;
            else {
                cerr << "[getBracketPair] Error! The input newick tree is not binary" << endl;
                exit(1);
            }
        }
    }

    // identify which are leaves
    for (k=0; k<n; k++) {
        if (tchildren[k].first == -1 && tchildren[k].second == -1) {
            isLeaf.push_back(true);
        } else {
            isLeaf.push_back(false);
        }
    }
    
    // reorder the nodes such that all leaves appear first
    changeTo.resize(n);
    changeFr.resize(n);
    t=0;
    // for all leaves
    for (k=0; k<n; k++) {
        if (isLeaf[k]) {
            changeTo[t]=k;
            changeFr[k]=t;
            t++;
        }
    }
    // for all internal nodes
    for (k=0; k<n; k++) {
        if (!isLeaf[k]) {
            changeTo[t]=k;
            changeFr[k]=t;
            t++;
        }
    }
    // update the arrays
    for (k=0; k<n; k++) {
        t = changeTo[k];
        ipos.push_back(tpos[t]);
        if (tparents[t] != -1)
            iparents.push_back(changeFr[tparents[t]]);
        else
            iparents.push_back(-1);
        if (tchildren[t].first != -1 && tchildren[t].second != -1)
            ichildren.push_back(pair<int,int>(changeFr[tchildren[t].first],changeFr[tchildren[t].second]));
        else
            ichildren.push_back(pair<int,int>(-1,-1));
    }
/*
    for (k=0; k<n; k++) {
        cout << k << "\t" << tpos[k].first << "\t" << tpos[k].second << "\t" << tparents[k] << "\t" << tchildren[k].first << "\t" << tchildren[k].second << "\t" << isLeaf[k];
        if (tchildren[k].first == -1 || tchildren[k].second==-1)
            cout << "\t" << topStr.substr(tpos[k].first, tpos[k].second-tpos[k].first+1);
        cout << endl;
    }
    
    cout << "changeTo:" << endl;
    for (k=0; k<changeTo.size(); k++) {
        cout << k << " -> " << changeTo[k] << endl;
    }
    
    for (k=0; k<n; k++) {
        cout << k << "\t" << ipos[k].first << "\t" << ipos[k].second << "\t" << iparents[k] << "\t" << ichildren[k].first << "\t" << ichildren[k].second;
        if (ichildren[k].first == -1 || ichildren[k].second==-1)
            cout << "\t" << topStr.substr(ipos[k].first, ipos[k].second-ipos[k].first+1);
        cout << endl;
    }
*/
    
}

// convert the binary toplogy from text format into integer format
void UnlabelTree::topTxt2Int(string& topStr, vector<int>& topInt, vector<double>& haploFreq, vector<double>& edgeLen) {
    
    // get corresponding pairs of brackets
    vector<pair<int,int> > pos;
    vector<int> parents;
    vector<pair<int,int> > children;
    vector<string> token;
    int k,s,t;
    string seq;
    double l;

    topInt.clear();
    haploFreq.clear();
    edgeLen.clear();
    getPosPair(topStr, pos, parents, children);
    for (k=0; k<(int)pos.size(); k++) {
        topInt.push_back(children[k].first);
        topInt.push_back(children[k].second);
        s = pos[k].first;
        t = pos[k].second;
        seq = topStr.substr(s,t-s+1);
        if (seq[seq.length()-1] == ')') {
            // root
            l = 0.0;
        } else {
            tokenizer(topStr.substr(s,t-s+1), "=:", &token);
            l = atof(token[token.size()-1].c_str());
            if (children[k].first == -1 && children[k].second==-1) {
                haploFreq.push_back(atof(token[1].c_str()));
            }
        }
        edgeLen.push_back(l);
    }
}


// show the topologies
void UnlabelTree::showTopologies(vector<vector<int> >& topologies) {
	int i;
	for (i=0; i<topologies.size(); i++) {
		cout << topInt2Txt(topologies.at(i)) << endl;
	}
}

// get the tree depth
int UnlabelTree::treeDepth(vector<int>& topInt) {
	int numNodes = (int)topInt.size() / 2;
	int* nodeDepths = new int[numNodes];
	int lnode, rnode;
	int treeDepth;
	int i;
	
	for (i=0; i<numNodes; i++) {
		lnode = topInt[i*2];
		rnode = topInt[i*2+1];
		if (lnode == -1 || rnode == -1) {
			nodeDepths[i] = 0;
		} else {
			nodeDepths[i] = maxInt(nodeDepths[lnode], nodeDepths[rnode])+1;
		}
	}
	treeDepth = nodeDepths[numNodes-1];
	delete[] nodeDepths;
	return treeDepth;
}

// get the tree depth
int UnlabelTree::treeDepth(vector<int>& topInt, vector<int>& nodeOrder) {
	int numNodes = (int)topInt.size() / 2;
	int* nodeDepths = new int[numNodes];
	int lnode, rnode;
	int treeDepth;
	int i,k;
	
	for (k=0; k<numNodes; k++) {
		i = nodeOrder[k];
		lnode = topInt[i*2];
		rnode = topInt[i*2+1];
		if (lnode == -1 || rnode == -1) {
			nodeDepths[i] = 0;
		} else {
			nodeDepths[i] = maxInt(nodeDepths[lnode], nodeDepths[rnode])+1;
		}
	}
	treeDepth = nodeDepths[numNodes-1];
	delete[] nodeDepths;
	return treeDepth;
}

// change the root position to the edge connecting to the node x (i.e. newRoot)
void UnlabelTree::changeRoot(vector<int>& topInt, int newRoot, vector<int>& newTop, vector<int>& newEdge, vector<int>& newNodeOrder) {
	int numNodes = (int)topInt.size() / 2;
	int numTips = (numNodes + 1) / 2;
	int* parents = new int[numNodes-1];
	int* newParents = new int[numNodes-1];
	int lnode, rnode, cnode;
	vector<int> nodesToProcess;
	int i,k;
	
	// check the value of newRoot
	if (newRoot >= numNodes-1) {
		cerr << "[UnlabelTree::changeRoot] Error! The value of newRoot >= the number of nodes - 1 (i.e. " << numNodes-1 << ")" << endl;
		exit(1);
	}
	
	// for all the nodes except root
	for (i=0; i<numNodes; i++) {
		lnode = topInt[2*i];
		rnode = topInt[2*i+1];
		if (lnode == -1 || rnode == -1)
			continue;
		if (i < numNodes-1) {
			// not a root
			parents[lnode] = i;
			parents[rnode] = i;
		} else {
			// a root
			parents[lnode] = rnode;
			parents[rnode] = lnode;
		}
	}
	
	// initialize the new topology
	newTop.clear();
	newTop.insert(newTop.begin(), topInt.begin(), topInt.end());
	
	// initialize the new node order
	newNodeOrder.resize(numNodes);
	for (i=0; i<numTips; i++) {
		newNodeOrder[i] = i;
	}
	
	// initialize the new edge
	newEdge.resize(numNodes-1);
	
	// build the new topology according to the new position of the root
	lnode = newRoot;
	rnode = parents[newRoot];
	cnode = numNodes-1;
	newTop[cnode*2] = lnode;
	newTop[cnode*2+1] = rnode;
	newParents[lnode] = rnode;
	newParents[rnode] = lnode;
	nodesToProcess.push_back(lnode);
	nodesToProcess.push_back(rnode);
	newNodeOrder[numNodes-1] = cnode;
	newEdge[lnode] = -1; // the edge length of the new root's left child is always ZERO
	newEdge[rnode] = minInt(lnode, rnode);
	i=0;k=1;
	while (i<nodesToProcess.size()) {
		cnode = nodesToProcess[i];
		if (cnode >= numTips) {
			// not a leaf
			if (newParents[cnode] != parents[cnode]) {
				// build the new subtree
				// get the new children
				if (topInt[cnode*2] == newParents[cnode]) {
					lnode = topInt[cnode*2+1];
					rnode = parents[cnode];
				} else {
					lnode = topInt[cnode*2];
					rnode = parents[cnode];
				}
				newTop[cnode*2] = lnode;
				newTop[cnode*2+1] = rnode;
			} else {
				lnode = topInt[cnode*2];
				rnode = topInt[cnode*2+1];
			}
			newParents[lnode] = cnode;
			newParents[rnode] = cnode;
			nodesToProcess.push_back(lnode);
			nodesToProcess.push_back(rnode);
			// set the new node order
			k++;
			newNodeOrder[numNodes-k] = cnode;
			// set the new edge
			newEdge[lnode] = minInt(lnode, cnode);
			newEdge[rnode] = minInt(rnode, cnode);
		}
		i++;
	}
	
	// clear the memory
	delete[] parents;
	delete[] newParents;
}

// obtain a reprsentative for the rooted topologies
// only consider the topologies when rooted on the leaf edge
// and pick the smallest one (w.r.t. the alphabetical order)
string UnlabelTree::getRepresent(vector<int>& topInt, vector<int>& representTop) {
	int i, besti;
	int numNodes = (int)topInt.size() / 2;
	int minTreeDepth = numNodes;
	int cTreeDepth;
	string represent, ctxt;
	vector<int> newTop;
	vector<int> newNodeOrder;
	vector<int> newEdge;
	represent = "";
	besti = 0;
	// cout << "input tree: " << topInt2Txt(topInt) << endl;
	for (i=0; i<numNodes-1; i++) {
		changeRoot(topInt, i, newTop, newEdge, newNodeOrder);
		// showTopology(newTop);
		cTreeDepth = treeDepth(newTop, newNodeOrder);
		if (cTreeDepth > minTreeDepth)
			continue;
		ctxt = topInt2Txt(newTop, newNodeOrder);
		// cout << ctxt << endl;
		if (cTreeDepth<minTreeDepth || ctxt<represent) {
			represent = ctxt;
			minTreeDepth = cTreeDepth;
			besti = i;
		}
	}
	changeRoot(topInt, besti, newTop, newEdge, newNodeOrder);
	updateTopology(newTop, newNodeOrder, representTop);
	return represent;
}

// update the topology so that the node order becomes 0,1,2,3,....
void UnlabelTree::updateTopology(vector<int>& topInt, vector<int>& nodeOrder, vector<int>& newTop) {
	vector<int> changeTo;
	int i,k;
	int left,right;
	newTop.clear();
	changeTo.resize(nodeOrder.size());
	for (i=0; i<nodeOrder.size(); i++) {
		changeTo[nodeOrder[i]] = i;
	}
	for (k=0; k<topInt.size()/2; k++) {
		i = nodeOrder[k];
		left = topInt[2*i];
		right = topInt[2*i+1];
		if (left != -1)
			left = changeTo[left];
		if (right != -1)
			right = changeTo[right];
		if (left <= right) {
			newTop.push_back(left);
			newTop.push_back(right);
		} else {
			newTop.push_back(right);
			newTop.push_back(left);
		}
	}
}

// show the topology
void UnlabelTree::showTopology(vector<int>& topInt) {
	int i;
	for (i=0; i<topInt.size(); i+=2) {
		cout << topInt[i] << "," << topInt[i+1] << endl;
	}
}

// show the topology
void UnlabelTree::showTopology(vector<int>& topInt, vector<int>& newNodeOrder) {
	int i,k;
	for (i=0; i<topInt.size()/2; i++) {
		k = newNodeOrder[i];
		cout << k << "," << topInt[2*k] << "," << topInt[2*k+1] << endl;
	}
}

// show the topology
void UnlabelTree::showTopology(vector<int>& topInt, vector<int>& newEdge, vector<int>& newNodeOrder) {
	int i,k;
	for (i=0; i<topInt.size()/2; i++) {
		k = newNodeOrder[i];
		cout << k << "," << topInt[2*k];
		if (topInt[2*k] >= 0)
			cout << "(edge:" << newEdge[topInt[2*k]] << ")";
		cout << "," << topInt[2*k+1];
		if (topInt[2*k+1] >= 0)
			cout << "(edge:" << newEdge[topInt[2*k+1]] << ")";
		cout << endl;
	}
}

// show the topology in one line
void UnlabelTree::showTopologyInOneLine(vector<int>& topInt) {
	int i;
	for (i=0; i<topInt.size(); i++) {
		if (i > 0)
			cout << ",";
		cout << topInt[i];
	}
	cout << endl;
}

// get the sister node
// return -1 if not found
int UnlabelTree::getSisterNode(int node1, vector<int>& topInt) {
    int i;
    for (i=0; i<topInt.size(); i++) {
        if (topInt[i] == node1) {
            if (i%2==0)
                return topInt[i+1];
            else
                return topInt[i-1];
        }
    }
    return -1;
}

// NNI
// update: topology is always changed
bool UnlabelTree::doNNI(vector<int>& topInt, int numTips, MyRand& myrand, vector<pair<int,int> >& changeNodeIDs) {
    
    int i;
    changeNodeIDs.clear();
    int which = myrand.iunif(0,1);
    if (which==2) {
        // no change, exit
        return false;
    }
    int line = myrand.iunif(numTips, topInt.size()/2-3);
    int selectNode1Idx = 2*line + which;
    int selectNode2Idx = -1;
    for (i=(line+1)*2; i<topInt.size(); i++) {
        if (topInt[i] == line) {
            if (i%2==0)
                selectNode2Idx = i+1;
            else
                selectNode2Idx = i-1;
            break;
        }
    }
    if (selectNode2Idx == -1) {
        cerr << "[UnlabelTree::doNNI] Error! selectNode2Idx == -1" << endl;
        exit(1);
    }
    // swap between "topInt[selectNode1Idx]" and "topInt[selectNode2Idx]"
    i = topInt[selectNode1Idx];
    topInt[selectNode1Idx] = topInt[selectNode2Idx];
    topInt[selectNode2Idx] = i;
    /*
    cout << "line = " << line << endl;
    cout << "which = " << which << endl;
    cout << "selectNode1Idx = " << selectNode1Idx << endl;
    cout << "selectNode2Idx = " << selectNode2Idx << endl;
    
    // show the topology
    cout << "topology before tidy up:" << endl;
    showTopology(topInt);
    */
    // tidy up the messy topology
    tidyUpMessyTopology(numTips, topInt, changeNodeIDs);
    
    return true;
}

// subtree swapping
// isDescendant[i*numNodes+j] = true if node j is descendant of node i (relationship between two nodes)
// return true if the topology is changed
bool UnlabelTree::swapSubTree(vector<int>& topInt, int numTips, int numNodes, MyRand& myrand, vector<bool>& isDescendant, vector<pair<int,int> >& changeNodeIDs) {
    
    vector<int> availableToChoose;
    int node1, node2, sister_node1;
    int i;
    node1 =myrand.iunif(0, 2*numTips-2);
    sister_node1 = getSisterNode(node1, topInt);
    changeNodeIDs.clear();
    if (sister_node1==-1)
        return false;
    
    for (i=0; i<=2*numTips-2; i++) {
        if ((i!=node1) && (i!=sister_node1) && (!isDescendant[node1*numNodes+i]) && (!isDescendant[i*numNodes+node1])) {
            availableToChoose.push_back(i);
        }
    }
    if (availableToChoose.size()==0)
        return false;
    node2 = availableToChoose[myrand.iunif(0,availableToChoose.size()-1)];
    // swap between node1 and node2
    for (i=0; i<topInt.size(); i++) {
        if (topInt[i] == node1)
            topInt[i] = node2;
        else if (topInt[i] == node2)
            topInt[i] = node1;
    }
    // tidy up the messy topology
    tidyUpMessyTopology(numTips, topInt, changeNodeIDs);
    return true;
}

// tidy up the messy topology
void UnlabelTree::tidyUpMessyTopology(int numTips, vector<int>& topIntMessy, vector<pair<int,int> >& changeNodeIDs) {
    
    int numLines = topIntMessy.size() / 2;
    if (numLines != (numTips*2-1)) {
        cerr << "[UnlabelTree::tidyUpMessyTopology] Error! The size of 'topIntMessy' does not match with the value of 'numTips'" << endl;
        exit(1);
    }
    vector<int> old2newNodes;
    vector<bool> lineExamined;
    vector<int> topInt;
    int leftChild, rightChild, interNode;
    bool leftChildExamined, rightChildExamined;
    int newLeftChildID, newRightChildID, newInterNodeID;
    int i;
    
    changeNodeIDs.clear();
    for (i=0; i<numLines; i++)
        lineExamined.push_back(false);
    for (i=0; i<numLines; i++)
        old2newNodes.push_back(-1);
    bool finish = false;
    while (!finish) {
        finish = true;
        for (i=0; i<numLines; i++) {
            if (!lineExamined[i]) {
                leftChild = topIntMessy[i*2];
                rightChild = topIntMessy[i*2+1];
                leftChildExamined = (leftChild==-1 || old2newNodes[leftChild]>=0);
                rightChildExamined = (rightChild==-1 || old2newNodes[rightChild]>=0);
                if (!leftChildExamined || !rightChildExamined) {
                    // skip the line
                    finish = false;
                } else {
                    if (leftChild==-1)
                        newLeftChildID=-1;
                    else
                        newLeftChildID=old2newNodes[leftChild];
                    if (rightChild==-1)
                        newRightChildID=-1;
                    else
                        newRightChildID=old2newNodes[rightChild];
                    interNode = i;
                    newInterNodeID = topInt.size()/2;
                    // insert the line into topInt
                    topInt.push_back(newLeftChildID);
                    topInt.push_back(newRightChildID);
                    old2newNodes[interNode] = newInterNodeID;
                    lineExamined[i] = true;
                }
            }
        }
    }
    // get all the updated nodes
    for (i=0; i<old2newNodes.size(); i++) {
        if (old2newNodes[i] != i) {
            changeNodeIDs.push_back(pair<int,int>(i,old2newNodes[i]));
        }
    }
    // update the topology
    for (i=0; i<topIntMessy.size(); i++) {
        topIntMessy[i] = topInt[i];
    }
}

/*
// generate K possible labelled topologies based on an input topology
void UnlabelTree::genKLabelTopoloy(int numTips, vector<int>& topInt, int K, vector<vector<int> >& newTopInts) {
    
    long totalPossibles = fact(numTips);
    cout << totalPossibles << endl;
}*/

// constructors
Node::Node(int v) {
    value = v;
}
Node::Node(int v, Node* leftChild, Node* rightChild) {
    value = v;
    children.push_back(leftChild);
    children.push_back(rightChild);
}

// set the labels of the node and its children
void Node::setLabel(int parentLabel, int& lastLabel, vector<bool>& mutationOccur) {
    int i;
    if (mutationOccur[value]) {
        // mutation on the edge between the parent and this node
        lastLabel++;
        label = lastLabel;
    } else {
        // no mutation on the edge between the parent and this node
        label = parentLabel;
    }
    for (i=0; i<children.size(); i++) {
        children[i]->setLabel(label, lastLabel, mutationOccur);
    }
}

// load the topology
void RootTree::loadTopology(vector<int>& topInt) {
    int i;
    int numLines = topInt.size() / 2;
    int leftID, rightID;
    Node* leftChild;
    Node* rightChild;
    Node* newNode;
    for (i=0; i<numLines; i++) {
        leftID = topInt[2*i];
        rightID = topInt[2*i+1];
        if (leftID==-1 || rightID==-1) {
            // it is a leaf
            newNode = new Node(i);
            nodeList.push_back(newNode);
            leaveList.push_back(newNode);
        } else {
            // it is an internal node
            leftChild = nodeList[leftID];
            rightChild = nodeList[rightID];
            newNode = new Node(i, leftChild, rightChild);
            nodeList.push_back(newNode);
            leftChild->parent = newNode;
            rightChild->parent = newNode;
            if (i == numLines-1) {
                // this is a root
                rootNode = newNode;
            }
        }
    }
}

void rmInvisibleChar(string& s) {
    int i,k;
    k=0;
    for (i=0;i<s.length();i++) {
        if (s[i]>='!' && s[i]<='~') {
            if (i > k)
                s[k] = s[i];
            k++;
        }
    }
    if (k > 0 && k < s.length())
        s.resize(k);
    if (k == 0)
        s = "";
}

void collectLeaves(vector<pair<int, int> >& children, vector<int>& leafid, int interid, vector<int>& c_set) {
    stack<int> q;
    int i,x,y;
    
    // initialized c_set
    for (i=0; i<c_set.size(); i++)
        c_set[i] = 0;
    
    q.push(interid);
    while (!q.empty()) {
        i = q.top();
        q.pop();
        x = children[i].first;
        y = children[i].second;
        if (x==-1 && y==-1) {
            // it is a leaf
            c_set[leafid[i]] = 1;
            // cout << interid << "\t" << leafid[i] << endl;
        } else {
            q.push(x);
            q.push(y);
        }
    }
    
    // update c_set if c_set[0] != 0
    if (c_set[0] == 1) {
        for (i=0; i<c_set.size(); i++) {
            c_set[i] = (c_set[i] + 1) % 2;
        }
    }
}

// load the topology from a newick tree file
// and then obtain the mutation sets
void RootTree::get_mutation_sets(char* treeFile, vector<vector<int> >& c_sets, int& numHaps) {
    vector<int> topInt;
    vector<double> haploFreq;
    vector<double> edgeLen;
    vector<int> parents;
    vector<pair<int,int> > children;
    vector<pair<int,int> > pos;
    vector<int> leafid;
    vector<int> c_set;
    string treeStr;
    ifstream fin;
    int numNodes;
    int k,l;
    char c;
    fin.open(treeFile);
    while (getline(fin, treeStr)) {
        rmInvisibleChar(treeStr);
        if (treeStr.length() > 0 && treeStr[0] != '#') {
            break;
        }
    }
    fin.close();
    
    // cout << "treeStr: " << treeStr << endl;
    pos.clear();
    parents.clear();
    children.clear();
    leafid.clear();
    c_sets.clear();
    
    getPosPair(treeStr, pos, parents, children);
    
    // collect all the leaf ID (A:0;B:1;C:2;...)
    numNodes = pos.size();
    numHaps = (numNodes + 1) / 2;
    for (k=0; k<pos.size(); k++) {
        if (children[k].first == -1 && children[k].second == -1) {
            c = treeStr[pos[k].first];
            leafid.push_back(c-'A');
        }
    }
    
    // collect the possible mutation sets
    c_set.resize(numHaps);
    for (k=numHaps; k<pos.size()-2; k++) {
        collectLeaves(children, leafid, k, c_set);
        c_sets.push_back(c_set);
    }
}


// destructor
RootTree::~RootTree() {
    clear();
}

// clear
void RootTree::clear() {
    int i;
    for (i=0; i<nodeList.size(); i++)
        delete nodeList[i];
    nodeList.clear();
    leaveList.clear();
}

// print the topology
void RootTree::printTopology() {
    vector<Node*> toProcess;
    Node* cNode;
    int i;
    toProcess.push_back(rootNode);
    i=0;
    while (i < toProcess.size()) {
        cNode = toProcess[i];
        if (cNode->children.size() > 0) {
            cout << cNode->value << ": " << cNode->children[0]->value << ", " << cNode->children[1]->value << endl;
            toProcess.push_back(cNode->children[0]);
            toProcess.push_back(cNode->children[1]);
        }
        i++;
    }
}

// return the labels of the leaves given the positions of the mutations
void RootTree::getLabels(vector<bool>& mutationOccur, vector<int>& leaveLabels, string& s, int& numGrps) {
    int i;
    int lastLabel;
    set<int> oldInt;
    set<int>::iterator itr;
    vector<int> old2New;
    lastLabel = 0;
    rootNode->setLabel(0, lastLabel, mutationOccur);
    leaveLabels.clear();
    for (i=0; i<leaveList.size(); i++) {
        leaveLabels.push_back(leaveList[i]->label);
    }
    for (i=0; i<leaveLabels.size(); i++) {
        itr=oldInt.find(leaveLabels[i]);
        if (itr == oldInt.end()) {
            // not found
            if (old2New.size() < leaveLabels[i]+1) {
                old2New.resize(leaveLabels[i]+1);
            }
            old2New[leaveLabels[i]] = oldInt.size();
            oldInt.insert(leaveLabels[i]);
        }
    }
    numGrps = oldInt.size();
    s = "";
    for (i=0; i<leaveLabels.size(); i++) {
        leaveLabels[i] = old2New[leaveLabels[i]];
        if (i>0) {
            s.append(",");
        }
        s.append(int2str(leaveLabels[i]));
    }
}

// get all possible groups of the leaves due to the mutations on different edges
void RootTree::getAllPossibleLeafGroups(LeafGrps& leafGrps) {
    // initialize the boolean array "mutationOccur"
    
    vector<bool> mutationOccur;
    vector<int> leafGrp;
    set<string> grpStrs;
    set<string>::iterator itr;
    int numGrps;
    string s;
    int i;
    for (i=0; i<nodeList.size(); i++) {
        mutationOccur.push_back(false);
    }
    while (nextSet(mutationOccur)) {
        getLabels(mutationOccur, leafGrp, s, numGrps);
        itr = grpStrs.find(s);
        if (itr==grpStrs.end()) {
            // not found
            grpStrs.insert(s);
            leafGrps.addLeafGrp(leafGrp, numGrps);
        }
    }
}

// next possible set of mutation positions
bool RootTree::nextSet(vector<bool>& mutationOccur) {
    int i,p;
    p = mutationOccur.size();
    for (i=0; i<mutationOccur.size(); i++) {
        if (mutationOccur[i]) {
            mutationOccur[i]=false;
        } else {
            mutationOccur[i]=true;
            p=i;
            break;
        }
    }
    if (p >= mutationOccur.size()-2)
        return false;
    else
        return true;
}

void LeafGrpSameNum::clear() {
    leafGrps.clear();
}

int LeafGrpSameNum::size() {
    return leafGrps.size();
}

void LeafGrpSameNum::show() {
    int i,j;
    for (i=0; i<leafGrps.size(); i++) {
        for (j=0; j<leafGrps[i].size(); j++) {
            if (j>0)
                cout << ",";
            cout << leafGrps[i].at(j);
        }
        cout << endl;
    }
}

void LeafGrps::addLeafGrp(vector<int>& leafGrp, int grp_num) {
    if (allLeafGrps.size() < grp_num) {
        allLeafGrps.resize(grp_num);
    }
    allLeafGrps[grp_num-1].leafGrps.push_back(leafGrp);
}

void LeafGrps::clear() {
    int i;
    for (i=0; i<allLeafGrps.size(); i++)
        allLeafGrps[i].clear();
}

void LeafGrps::show() {
    int i;
    cout << "--------------------------" << endl;
    for (i=1; i<allLeafGrps.size(); i++) {
        cout << "Number of groups with " << i+1 << " subgroups: " << allLeafGrps[i].size() << endl;
        allLeafGrps[i].show();
        cout << "--------------------------" << endl;
    }
}

void LeafGrps::showSummary() {
    int i;
    for (i=1; i<allLeafGrps.size(); i++) {
        if (i>1)
            cout << ",";
        cout << allLeafGrps[i].size();
    }
    cout << endl;
}
