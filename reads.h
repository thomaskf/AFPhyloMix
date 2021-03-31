#ifndef _READS_
#define _READS_

#include<string>
#include<vector>
#include<map>
#include<fstream>
#include<iostream>
#include<cstdlib>
#include "fileHandler.h"
#include "mylib.h"
#include "definitions.h"

using namespace std;

class Nucl2Int {
public:
	char c[256];
	// constructor
	Nucl2Int();
};

class Insert {
public:
	int col;
	int readPos;
	string x;
	int readid;
	int len;
	Insert(int c, int p, string s, int id);
	void show();
};

class Read {
public:
	string seq;
    // string origSeq;
	int mapPos;
	int mapEndPos;
	string desc;
	int mapq;
	int frontTrimLen;
	int endTrimLen;
	bool deleted;
	string qualseq;

	// constructor
	Read(string s, int mPos, int mapQ, string des, int frontTrim, int endTrim, string qs);
    // Read(string s, string o_s, int mPos, int mapQ, string des, int frontTrim, int endTrim, string qs);
};

class Reads {
public:
	int pairReadNum;
	int singleReadNum;
	int readNum;
	int lastMapPos;
	vector<pair<Read*,Read*> > pairReads; // paired-end reads
	vector<Read*> unpairReads; // single-end reads;

	// load reads from the SAM file
	void readSamFile(char* samFile);

	// obtain the matrix for all positions
	// for each position, 0: '-'; 1: 'A'; 2: 'C'; 3: 'G'; 4: 'T'
	int* getPosMatrix(int& dim);

    ~Reads(); // destructor

	vector<Insert*> inserts;
	vector<Read*> readlist;
	// update the seq according to the cigar string
	// and store the insertions if exist
	// return false if the cigarStr is not valid
	bool updateSeq(string& seq, string& cigarStr, int readId, int mapPos, int& frontTrim, int& endTrim, string& qualseq);

	// update according to the insertions
	void updateForInserts();

    // check whether the trimming start position create issue
    // bool checkTrimmingStartPos(int pos, int len, double minRatio, int minFreq);
    
    // check whether the trimming end position create issue
    // bool checkTrimmingEndPos(int pos, int len, double minRatio, int minFreq);

    // get statistics of the trimming positions
    void getTrimmingStat(vector<int>& startTrimStat, vector<int>& endTrimStat, int numCol);
    
    // get the special region with significant coverage drop
    void getCoverDropRegs(vector<pair<int,int> >& dropRegions, int numCol, int* colStat);

    // change the bases to 'N' if its quality value is lower than the threshold
    void checkQual(int thres);
    
};

#endif
