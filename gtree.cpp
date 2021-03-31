#include "gtree.h"

void gNode::addChild(gNode* child) {
    children.push_back(child);
}

void gNode::getChildValues(gNode* parent, set<string>& values) {
	int i;
	for (i=0; i<children.size(); i++) {
    	if (children[i]!=parent) {
            if (childValues[i].size() == 0) {
				children[i]->getChildValues(this, childValues[i]);
			}
			values.insert(childValues[i].begin(), childValues[i].end());
    	}
    	if (isLeaf)
    		values.insert(value);
    }
}

string gNode::getPattern(gNode* parent) {
    vector<string> all_p;
    string s;
    int i;
    if (children.size() == 1) {
        // is a leaf
        return "x";
    } else {
        for (i=0; i<children.size(); i++) {
            if (children[i] != parent)
                all_p.push_back(children[i]->getPattern(this));
        }
        sort(all_p.begin(), all_p.end());
        s.append("(");
        for (i=0; i<all_p.size(); i++) {
            if (i>0)
                s.append(",");
            s.append(all_p[i]);
        }
        s.append(")");
    }
    return s;
}

string gNode::getNodeStr(gNode* parent) {
    vector<string> all_p;
    string s;
    int i;
    if (children.size() == 1) {
        // is a leaf
        return "f=" + doublToStr(freq, DECIMAL_PLACE);
    } else {
        for (i=0; i<children.size(); i++) {
            if (children[i] != parent)
                all_p.push_back(children[i]->getNodeStr(this) + ":" + doublToStr(edgeLens[i], DECIMAL_PLACE));
        }
        sort(all_p.begin(), all_p.end());
        s.append("(");
        for (i=0; i<all_p.size(); i++) {
            if (i>0)
                s.append(",");
            s.append(all_p[i]);
        }
        s.append(")");
    }
    return s;
}

void gNode::replaceChild(gNode* newChild, gNode* childToReplace, double newEdgeLen) {
    int i;
    for (i=0; i<children.size(); i++) {
        if (children[i] == childToReplace) {
            children[i] = newChild;
            edgeLens[i] = newEdgeLen;
            break;
        }
    }
}

// compute the sum_of_square difference between two nodes
bool gNode::SSDiff(gNode* aNode, gNode* thisParent, gNode* aParent, double& ssdiff, double diffThres, string& algn) {
    int i, j;
    double d1, d2, diff1, diff2;
    gNode *child1, *child2;
    gNode *aChild1, *aChild2;
    string algn1, algn2;
    bool hasValue;
    int state;
    if (children.size() != aNode->children.size())
        return false;
    if (children.size() == 1) {
        // this node and aNode are leaves
        d1 = fabs(freq - aNode->freq);
        if (d1 > diffThres)
            return false;
        ssdiff = d1 * d1;
        algn = name + "<->" + aNode->name;
        return true;
    } else {
        ssdiff = 0.0;
        child1=child2=aChild1=aChild2=NULL;
        for (i=0; i<children.size(); i++) {
            if (children[i] == thisParent)
                continue;
            if (child1==NULL)
                child1=children[i];
            else
                child2=children[i];
        }
        for (i=0; i<aNode->children.size(); i++) {
            if (aNode->children[i] == aParent)
                continue;
            if (aChild1==NULL)
                aChild1=aNode->children[i];
            else
                aChild2=aNode->children[i];
        }
        if (child2==NULL || aChild2==NULL) {
            // something wrong here
            cerr << "something wrong here" << endl;
            return false;
        }
        hasValue = false;
        if (child1->SSDiff(aChild1,this,aNode,d1,diffThres,algn1) && child2->SSDiff(aChild2,this,aNode,d2,diffThres,algn2)) {
            ssdiff = d1 + d2;
            state = 1;
            algn = algn1 + " " + algn2;
            hasValue = true;
        }
        if (child1->SSDiff(aChild2,this,aNode,d1,diffThres,algn1) && child2->SSDiff(aChild1,this,aNode,d2,diffThres,algn2)) {
            if (!hasValue || d1+d2 < ssdiff) {
                ssdiff = d1+d2;
                state = 2;
                algn = algn1 + " " + algn2;
                hasValue = true;
            }
        }
        return hasValue;
    }
}

// supporting function
int gTree::getCloseBracketPos(string& str, int openBracPos) {
    int i,k;
    k=0;
    for (i=openBracPos+1; i<str.length(); i++) {
    	if (str[i] == '(') {
    	    k++;
    	} else if (str[i] == ')') {
    	    if (k==0)
    	        return i;
    	    else
    	        k--;
    	}
    }
    return -1; // not a balanced brackets
}

// create a node for str[startPos ... endPos]
// this node has to be enclosed by a pair of brackets
// i.e. str[startPos] == '(' and str[endPos] == ')'
gNode* gTree::createInterNode(string& str, int startPos, int endPos) {
    int p, s, t;
    gNode* cNode;
    gNode* childNode;
    bool isLeaf;
    string v;
    double edgeLen;
    cNode = new gNode();
    childNode = NULL;
    isLeaf = (str[startPos] != '('); // if the first character is '(', then it is not a leaf
    
    p=startPos+1;
    while (p<=endPos) {
         if (str[p] == '(') {
             // an internal node
             s = p;
             t = getCloseBracketPos(str, s);
             if (t==-1) {
                 cerr << "Error! The brackets are not balanced" << endl;
                 exit(1);
             }
             childNode = createInterNode(str, s, t);
             p=t+1;
         } else {
             // a leaf node
             v = "";
             while (str[p]!=',' && str[p]!=':' && str[p]!=')') {
                 v.append(1, str[p]);
                 p++;
             }
             childNode = new gNode();
             childNode->name = v;
         }
         // get the edge length if available
         edgeLen = 0.0;
         if (str[p]==':') {
             v = "";
             p++;
             while (str[p]!=',' && str[p]!=')') {
                 v.append(1, str[p]);
                 p++;
             }
             edgeLen = atof(v.c_str());
         }
         // add the child
         if (str[p]==',' || str[p]==')') {
             cNode->addChild(childNode);
             childNode->addChild(cNode);
             cNode->edgeLens.push_back(edgeLen);
             childNode->edgeLens.push_back(edgeLen);
         }
         p++;
    }
    return cNode;
}

// create an internal node from a single roow in the topology matrix
gNode* gTree::createNodeFrRow(int rowID, vector<int>& topInt, vector<double>& haploFreq, vector<double>& edgeLen) {
    gNode* cNode;
    gNode* leftChildNode;
    gNode* rightChildNode;
    int leftChildID, rightChildID;
    cNode = new gNode();
    leftChildID = topInt[rowID*2];
    rightChildID = topInt[rowID*2+1];
    if (leftChildID == -1 || rightChildID == -1) {
        // this is a leaf
        cNode->freq = haploFreq[rowID];
    } else {
        // this is an internal node
        leftChildNode = createNodeFrRow(leftChildID, topInt, haploFreq, edgeLen);
        rightChildNode = createNodeFrRow(rightChildID, topInt, haploFreq, edgeLen);
        cNode->addChild(leftChildNode);
        cNode->addChild(rightChildNode);
        cNode->edgeLens.push_back(edgeLen[leftChildID]);
        cNode->edgeLens.push_back(edgeLen[rightChildID]);
        leftChildNode->addChild(cNode);
        rightChildNode->addChild(cNode);
        leftChildNode->edgeLens.push_back(edgeLen[leftChildID]);
        rightChildNode->edgeLens.push_back(edgeLen[rightChildID]);
    }
    return cNode;
}
// load the topology file
void gTree::loadTopFile(char* fileName) {
    int s, t;
    gNode* aNode;
    ifstream fin;
    string str;
    
    fin.open(fileName);
    getline(fin, str); // only one line in the file
    fin.close();
    
    loadTreeStr(str);
}

// load the tree sequence
void gTree::loadTreeStr(string str) {
    int s, t;
    gNode* aNode;

    s = 0;
    while (s < str.length()) {
        if (str[s] == '(')
            break;
        s++;
    }
    t = str.length() - 1;
    while (t >= 0) {
        if (str[t] == ')')
            break;
        t--;
    }
    aNode = createInterNode(str, s, t);
    root = aNode;
    
    // label all the internal nodes
    assignNodeID();
    
    // get all the nodes
    getAllNodes();
    
    // get all the frequencies
    getFreqs();
    
    // printTopology();
}



// load the topology matrix
void gTree::loadTopMatrix(vector<int>& topInt, vector<double>& haploFreq, vector<double>& edgeLen) {
    int numRow = topInt.size() / 2;
    root = createNodeFrRow(numRow-1, topInt, haploFreq, edgeLen);
    
    assignNodeID();
    
    getAllNodes();
}

// destructor
gTree::~gTree() {
    clear();
}

// clear
void gTree::clear() {
    vector<gNode*> toProcess;
    set<gNode*> processed;
    set<gNode*>::iterator itr;
    gNode* cNode;
    gNode* childNode;
    int i,j;
    
    i=0;
    toProcess.push_back(root);
    while (i < toProcess.size()) {
        cNode = toProcess[i];
        for (j=0; j<cNode->children.size(); j++) {
        	childNode = cNode->children[j];
        	itr = processed.find(childNode);
        	if (itr == processed.end()) {
                toProcess.push_back(childNode);
            }
        }
        processed.insert(cNode);
        i++;
    }
    for (i=0; i<toProcess.size(); i++)
        delete toProcess[i];
}

// print the topology
void gTree::printTopology() {
    vector<gNode*> toProcess;
    set<gNode*> processed;
    set<gNode*>::iterator itr;
    gNode* cNode;
    gNode* childNode;
    double edgeLen;
    int i,j;
    
    i=0;
    toProcess.push_back(root);
    while (i < toProcess.size()) {
        cNode = toProcess[i];
        if (cNode->children.size() == 1) {
            // is a leaf
            // cout << cNode->id << " is a leaf " << cNode->name << endl;
            cout << cNode->id << " is a leaf f=" << cNode->freq << endl;
        } else {
			cout << cNode->id << ":";
			for (j=0; j<cNode->children.size(); j++) {
				childNode = cNode->children[j];
				edgeLen = cNode->edgeLens[j];
				itr = processed.find(childNode);
				if (itr == processed.end()) {
					toProcess.push_back(childNode);
					cout << " " << childNode->id << "(" << edgeLen << ")";
				}
			}
			cout << endl;
			processed.insert(cNode);
		}
        i++;
    }
    // show the edges
    cout << "Edges:" << endl;
    for (i=0; i<edges.size(); i++) {
        cout << (edges[i].first)->id << "--" << (edges[i].second)->id << endl;
    }
    // show the leaves
    cout << "Leaves:" << endl;
    for (i=0; i<leaves.size(); i++) {
        cout << " " << leaves[i]->id;
    }
    cout << endl;
    // show allNodes
    cout << "All Nodes:" << endl;
    for (i=0; i<allNodes.size(); i++) {
        cout << " " << allNodes[i]->id;
    }
    cout << endl;
}

// get tree in string format
string gTree::getTreeStr() {
    return root->getNodeStr(NULL);
}


// reset all children values
void gTree::resetAllChildValues() {
    int i,j,n;
    gNode* cNode;
    for (i=0; i<allNodes.size(); i++) {
         cNode = allNodes[i];
         n = cNode->children.size();
         cNode->childValues.resize(n);
         for (j=0; j<n; j++) {
              cNode->childValues[j].clear();
         }
    }
}

// check whether more than one mutation
bool gTree::isMoreThanOneMutate() {
    gNode* node1;
    gNode* node2;
    set<string> values1, values2;
    set<string>::iterator itr;
    int i;
    resetAllChildValues();
    for (i=0; i<edges.size(); i++) {
    	node1 = edges[i].first;
    	node2 = edges[i].second;
    	values1.clear();
    	node1->getChildValues(node2, values1);
    	/*
    	cout << "Node 1's id: " << node1->id << " ";
    	cout << "values1:";
    	for (itr=values1.begin(); itr!=values1.end(); itr++)
    	    cout << " " << *itr;*/
    	values2.clear();
    	node2->getChildValues(node1, values2);
    	/*
    	cout << " Node 2's id: " << node2->id << " ";
    	cout << "values2:";
    	for (itr=values2.begin(); itr!=values2.end(); itr++)
    	    cout << " " << *itr;
    	cout << endl;
    	cout << "size of values1: " << values1.size();
    	cout << " size of values2: " << values2.size() << endl;*/
    	if (values1.size() == 1 && values2.size() == 1)
    	    return false;
    }
    return true;
}

// get the pattern
string gTree::getPattern() {
    return root->getPattern(NULL);
}

// change root to a terminal edge connecting to a leaf node
void gTree::changeRoot(gNode* leafNode) {
    gNode* n1;
    gNode* n2;
    gNode* newRoot;
    double edgeLen1, edgeLen2, newEdgeLen;
    
    // if this is a rooted tree, remove the root first
    if (root->children.size()==2) {
        n1 = root->children[0];
        n2 = root->children[1];
        edgeLen1 = root->edgeLens[0];
        edgeLen2 = root->edgeLens[1];
        newEdgeLen = edgeLen1 + edgeLen2;
        n1->replaceChild(n2, root, newEdgeLen);
        n2->replaceChild(n1, root, newEdgeLen);
        rmEdge(root, n1);
        rmEdge(root, n2);
        rmInterNode(root);
        edges.push_back(pair<gNode*,gNode*>(n1,n2));
    }
    
    // put a root on the terminal edge
    newRoot = new gNode();
    newRoot->isLeaf = false;
    n1 = leafNode->children[0]; // leaf node only has one child
    edgeLen1 = leafNode->edgeLens[0];
    newEdgeLen = edgeLen1/2.0;
    leafNode->replaceChild(newRoot, n1, newEdgeLen);
    n1->replaceChild(newRoot, leafNode, newEdgeLen);
    rmEdge(leafNode, n1);
    newRoot->addChild(leafNode);
    newRoot->addChild(n1);
    newRoot->edgeLens.push_back(newEdgeLen);
    newRoot->edgeLens.push_back(newEdgeLen);
    edges.push_back(pair<gNode*,gNode*>(newRoot,leafNode));
    edges.push_back(pair<gNode*,gNode*>(newRoot,n1));
    allNodes.push_back(newRoot);
    newRoot->id = root->id;
    delete root;
    root = newRoot;
}

// RMS difference between another rooted tree
// report false if another rooted tree and this rooted tree are not in the same topology
bool gTree::RMSDiff(gTree* aTree, double& diff, double diffThres, string& algn) {
    double d;
    if (root->SSDiff(aTree->root, NULL, NULL, d, diffThres, algn)) {
        diff = sqrt(d / (double)leaves.size());
        return true;
    } else {
        return false;
    }
}

// remove an internal node from the tree
// assume no children array containing that node
void gTree::rmInterNode(gNode* nodeToRemove) {
    int i,k;
    k=0;
    // remove that node from the array "allNodes"
    for (i=0; i<allNodes.size(); i++) {
        if (allNodes[i] != nodeToRemove) {
            if (k < i) {
                allNodes[k] = allNodes[i];
            }
            k++;
        }
    }
    allNodes.resize(k);
    // remove the edges with that node
    k=0;
    for (i=0; i<edges.size(); i++) {
        if (edges[i].first != nodeToRemove && edges[i].second != nodeToRemove) {
            if (k < i) {
                edges[k] = edges[i];
            }
            k++;
        }
    }
    edges.resize(k);
}

// remove an edge
void gTree::rmEdge(gNode* n1, gNode* n2) {
    int i,k;
    k=0;
    for (i=0; i<edges.size(); i++) {
        if ((edges[i].first == n1 && edges[i].second==n2) || (edges[i].first == n2 && edges[i].second==n1)) {
            // remove the edge
        } else {
            if (k < i) {
                edges[k] = edges[i];
            }
            k++;
        }
    }
    edges.resize(k);
}

// assign the IDs of the nodes
void gTree::assignNodeID() {
    vector<gNode*> toProcess;
    set<gNode*> processed;
    set<gNode*>::iterator itr;
    gNode* cNode;
    gNode* childNode;
    int i,j;
    
    i=0;
    toProcess.push_back(root);
    while (i < toProcess.size()) {
        cNode = toProcess[i];
        cNode->id = i;
        for (j=0; j<cNode->children.size(); j++) {
        	childNode = cNode->children[j];
        	itr = processed.find(childNode);
        	if (itr == processed.end()) {
            	toProcess.push_back(childNode);
            }
        }
        processed.insert(cNode);
        i++;
    }
}

// get all the nodes
void gTree::getAllNodes() {
    vector<gNode*> toProcess;
    set<gNode*> processed;
    set<gNode*>::iterator itr;
    gNode* cNode;
    gNode* childNode;
    double edgeLen;
    int i,j;
    
    i=0;
    leaves.clear();
    toProcess.push_back(root);
    while (i < toProcess.size()) {
        cNode = toProcess[i];
        allNodes.push_back(cNode);
        if (cNode->children.size() == 1) {
            // is a leaf
            leaves.push_back(cNode);
            cNode->isLeaf = true;
        } else {
            cNode->isLeaf = false;
			for (j=0; j<cNode->children.size(); j++) {
				childNode = cNode->children[j];
				edgeLen = cNode->edgeLens[j];
				itr = processed.find(childNode);
				if (itr == processed.end()) {
					toProcess.push_back(childNode);
					edges.push_back(pair<gNode*,gNode*>(cNode,childNode));
				}
			}
			processed.insert(cNode);
		}
        i++;
    }
}

 // get all the frequencies
void gTree::getFreqs() {
    int i;
    for (i=0; i<leaves.size(); i++) {
        // cout << leaves[i]->name << endl << flush;
        leaves[i]->freq = atof(leaves[i]->name.substr(2).c_str());
    }
}



