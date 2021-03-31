#include "reads.h"

//=====================================================
// For the class Nucl2Int
//=====================================================

// constructor
Nucl2Int::Nucl2Int() {
	// '_'->0, A->1, C->2, G->3, T->4
	memset(c, 0, 256);
	c[(int)'A'] = c[(int)'a'] = 1;
	c[(int)'C'] = c[(int)'c'] = 2;
	c[(int)'G'] = c[(int)'g'] = 3;
	c[(int)'T'] = c[(int)'t'] = 4;
}

//=====================================================
// For the class Insert
//=====================================================

bool cmpInsert(Insert* i1, Insert* i2) {
	if (i1->col == i2->col) {
		if (i1->len == i2->len) {
			return i1->readid < i2->readid;
		} else {
			return i1->len > i2->len;
		}
	}
	return i1->col > i2->col;
}

Insert::Insert(int c, int p, string s, int id) {
	col = c;
	readPos = p;
	x = s;
	readid = id;
	len = (int) s.length();
}

void Insert::show() {
	cout << col << "\t" << readPos << "\t" << x << "\t" << readid << endl;
}

//=====================================================
// For the class Read
//=====================================================

// constructor
// Read::Read(string s, string o_s, int mPos, int mapQ, string des, int frontTrim, int endTrim, string qs) {
Read::Read(string s, int mPos, int mapQ, string des, int frontTrim, int endTrim, string qs) {
	seq = s;
    // origSeq = o_s;
	mapPos = mPos;
	mapEndPos = mPos + (int) s.length() - 1;
	mapq = mapQ;
	desc = des;
	frontTrimLen = frontTrim;
	endTrimLen = endTrim;
	deleted = false;
	qualseq = qs;
}


// load reads from the SAM/BAM file
void Reads::readSamFile(char* samFile) {
	string aline;
	int mpos;
	int mapq;
	int flag;
	int indx;
	int frontTrim;
	int endTrim;
	string cigar;
    // string o_seq;
	string seq;
	string qualseq;
	vector<string> token;
	SamBamFileHander fileHander;
	string preSeqName = "";
	string currSeqName;
	Read* preRead = NULL;
	Read* currRead;
	indx = 0;
	lastMapPos = 0;
	// clear the arrays
    readlist.clear();
	pairReads.clear();
	unpairReads.clear();
	// open the file
	fileHander.openFile(samFile);
	while (fileHander.getNextSeq(aline)) {
		if (aline.length() > 0 && aline[0]!='@') {
			tokenizer(aline, "\t", &token);
			if (token.size() > 9) {
				currSeqName = token[0];
				flag = atoi(token[1].c_str());
				if ((flag&2048)!=0) {
					// this is a supplementry alignment
					// discard it
					continue;
				}
				mpos = atoi(token[3].c_str());
				if (mpos==0) {
					// this read is not aligned
					continue;
				}
				mpos = mpos-1; // change to zero-based
				mapq = atoi(token[4].c_str());
				cigar = token[5];
				if (cigar == "*") {
					// this read is not aligned
					continue;
				}
				seq = token[9];
                // o_seq = seq;
				if (mapq < MAPQ_THRES) {
					// the mapping quality is too low
					// cout << aline << endl;
					continue;
				}
				qualseq = token[10];
				if (updateSeq(seq, cigar, indx, mpos, frontTrim, endTrim, qualseq)) {
					// alignment is valid
					/*
					// check whether the aligned region is too short
					if (seq.length() < origSeqLen * LEN_THRES) {
						// the read is too short
						// cout << seq << "\t" << seq.length() << "\t" << origSeqLen << endl;
						continue;
					}*/
					// currRead = new Read(seq, o_seq, mpos, mapq, currSeqName, frontTrim, endTrim, qualseq);
                    currRead = new Read(seq, mpos, mapq, currSeqName, frontTrim, endTrim, qualseq);
					readlist.push_back(currRead);
					indx++;
					if (lastMapPos < mpos + (int) seq.length() - 1) {
						lastMapPos = mpos + (int) seq.length() - 1;
					}
					if (preRead!=NULL) {
						if (currSeqName == preSeqName) {
							// pair-ended reads
							if (preRead->mapPos <= currRead->mapPos) {
								pairReads.push_back(pair<Read*,Read*>(preRead,currRead));
							} else {
								pairReads.push_back(pair<Read*,Read*>(currRead,preRead));
							}
							preRead = NULL;
							preSeqName = "";
						} else {
							// single-end read
							unpairReads.push_back(preRead);
							preRead = currRead;
							preSeqName = currSeqName;
						}
					} else {
						// first read (of the pair)
						preRead = currRead;
						preSeqName = currSeqName;
					}
				}
			}
		}
	}
	if (preRead!=NULL) {
		// single-end read
		unpairReads.push_back(preRead);
	}
    
	fileHander.closeFile();

	// set the number
	readNum = (int) readlist.size();
	pairReadNum = (int) pairReads.size();
	singleReadNum = (int) unpairReads.size();

	// update the inserts
	// updateForInserts();
	
	// cerr << "lastMapPos = " << lastMapPos << endl;
    
    // change the bases to 'N' if its quality value is lower than the threshold
    checkQual(QV_THRES);
	
}

// update the seq according to the cigar string
// and store the insertions if exist
// return false if the cigarStr is not valid
bool Reads::updateSeq(string& seq, string& cigarStr, int readId, int mapPos, int& frontTrim, int& endTrim, string& qualseq) {
	int i,j,k;
	int seqI;
	
	i=0;
	seqI = 0;
	frontTrim = 0;
	endTrim = 0;
	for (j=0; j<cigarStr.length(); j++) {
		// should begin with a number, and followed by a char
		if (isdigit(cigarStr[j]))
			continue;
		else {
			if (j==i) {
				// cout << "Warning! Invalid cigar string: " << cigarStr << endl;
				return false;
			}
			k = atoi(cigarStr.substr(i,j-i).c_str());
			switch (cigarStr[j]) {
				case 'I':
					inserts.push_back(new Insert(mapPos+seqI,seqI,seq.substr(seqI,k),readId));
				case 'S':
					seq.erase(seqI,k);
					qualseq.erase(seqI,k);
					if (seqI == 0)
						frontTrim = k;
					if (j==cigarStr.length()-1)
						endTrim = k;
					break;
				case 'D':
				case 'N':
					seq.insert(seqI,k,'-');
					qualseq.insert(seqI,k,0);
					seqI += k;
					break;
				case 'M':
					seqI += k;
					break;
				case 'P':
				case 'H':
					// do nothing
					break;
			}
			i=j+1;
		}
	}
	return true;
}

// update according to the insertions
void Reads::updateForInserts() {
	
	// sort the insertions
	sort(inserts.begin(), inserts.end(), cmpInsert);
	
	cerr << "Number of insertions: " << inserts.size() << endl;
	
	// update the positions and sequences
	int i,j;
	int cpos,rpos;
	string s;
	int rid;
	int l;
	int pre_cpos = -1;
	Read* curRd;
	for (i=0; i<inserts.size(); i++) {
		cpos = inserts[i]->col;
		rpos = inserts[i]->readPos;
		s = inserts[i]->x;
		rid = inserts[i]->readid;
		l = inserts[i]->len;
		if (cpos != pre_cpos) {
			// create l columns
			for (j=0; j<readlist.size(); j++) {
				curRd = readlist[j];
				if (curRd->mapPos >= cpos) {
					curRd->mapPos += l;
				} else if (curRd->mapPos < cpos && curRd->seq.length()+curRd->mapPos > cpos) {
					curRd->seq.insert(cpos - curRd->mapPos, l, '-');
				}
			}
			pre_cpos = cpos;
		}
		// put the characters into the corresponding positions of the read
		for (j=0; j<l; j++)
			readlist[rid]->seq.at(rpos+j) = s[j];
	}
    // update the lastMapPos
    for (i=0; i<readlist.size(); i++) {
        if (lastMapPos < readlist[i]->mapPos + (int) readlist[i]->seq.length() - 1) {
            lastMapPos = readlist[i]->mapPos + (int) readlist[i]->seq.length() - 1;
        }
    }
}

// obtain the matrix for all positions
// for each position, 0: '-'; 1: 'A'; 2: 'C'; 3: 'G'; 4: 'T'
int* Reads::getPosMatrix(int& dim) {
	int i,j,k;
	int mpos,pos;
	char ch;
	int* posMatrix;
	Read* rd;
	Nucl2Int nucl2int;
	
	dim = lastMapPos + 1;
	posMatrix = new int[dim * 5];
	memset(posMatrix, 0, dim * 5 * sizeof(int));
	
	for (i=0; i<(int)readlist.size(); i++) {
		rd = readlist[i];
		mpos = rd->mapPos;
		for (j=0; j<rd->seq.length(); j++) {
			pos = mpos + j;
			ch = rd->seq.at(j);
            if (ch != 'N') {
                k = nucl2int.c[(int)ch];
                posMatrix[pos * 5 + k]++;
            }
		}
	}
	return posMatrix;
}

// get statistics of the trimming positions
void Reads::getTrimmingStat(vector<int>& startTrimStat, vector<int>& endTrimStat, int numCol) {
    startTrimStat.clear();
    startTrimStat.insert(startTrimStat.begin(), numCol+1, 0);
    endTrimStat.clear();
    endTrimStat.insert(endTrimStat.begin(), numCol+1, 0);
    int i,p;
    Read* currRd;
    for (i=0; i<readlist.size(); i++) {
        currRd = readlist[i];
        if (currRd->frontTrimLen > TOO_LONG_TRIM) {
            startTrimStat[currRd->mapPos]++;
        }
        if (currRd->endTrimLen > TOO_LONG_TRIM) {
            p = currRd->mapPos + currRd->seq.length();
            endTrimStat[p]++;
        }
    }
}

// get the special region with significant coverage drop
void Reads::getCoverDropRegs(vector<pair<int,int> >& dropRegions, int numCol, int* colStat) {
    int i,j;
    vector<int> startTrimStat;
    vector<int> endTrimStat;
    vector<int> candiStart;
    vector<int> candiEnd;
    // vector<int> potStart;
    // vector<int> potEnd;
    vector<pair<int,int> > candidates;
    vector<int> totColStat;
    int subtot;
    double avgCoverage = 0.0;
    double thres;
    // bool hasIssue;
    int s,t;
    
    // initialization
    dropRegions.clear();
    
    // get the subtotal of coverage for each column
    for (i=0; i<numCol; i++) {
        subtot = 0;
        for (j=0; j<5; j++) {
            subtot += colStat[i*5+j];
        }
        totColStat.push_back(subtot);
        avgCoverage+=subtot;
    }
    avgCoverage = avgCoverage / numCol;
    // cerr << "avgCoverage = " << avgCoverage << endl << flush;
    thres = avgCoverage * MINDROPRATIO;
    if (thres < MINDROPFREQ)
        thres = MINDROPFREQ;
    // cerr << "thres = " << thres << endl << flush;
    
    getTrimmingStat(startTrimStat, endTrimStat, numCol);
    
    /*
    // list out the startTrimStat
    cerr << "startTrimStat" << endl << flush;
    for (i=0; i<numCol; i++)
    cout << i << "\t" << startTrimStat[i] << endl << flush;
    
    // list out the endTrimStat
    cerr << "endTrimStat" << endl << flush;
    for (i=0; i<numCol; i++)
    cout << i << "\t" << endTrimStat[i] << endl << flush;
    */
    
    for (i=1; i<numCol; i++) {
        // if (startTrimStat[i] > thres && totColStat[i-1] < totColStat[i] * (1.0 - DROP_LEVEL)) {
        if (startTrimStat[i] > thres) {
            candiStart.push_back(i);
        }
    }
    
    for (i=1; i<numCol; i++) {
        // if (endTrimStat[i] > thres && totColStat[i-1] * (1.0 - DROP_LEVEL) > totColStat[i]) {
        if (endTrimStat[i] > thres) {
            candiEnd.push_back(i);
        }
    }
    
    // combine if they are too close
    j=0;
    for (i=0; i<candiStart.size(); i++) {
        if (i+1 >= candiStart.size() || candiStart[i+1]-candiStart[i] > TOO_CLOSE_DROP) {
            // not too close
            if (i > j) {
                candiStart[j] = candiStart[i];
            }
            j++;
        }
    }
    if (j < candiStart.size())
        candiStart.resize(j);
    j=0;
    for (i=0; i<candiEnd.size(); i++) {
        if (i-1 < 0 || candiEnd[i]-candiEnd[i-1] > TOO_CLOSE_DROP) {
            // not too close
            if (i > j) {
                candiEnd[j] = candiEnd[i];
            }
            j++;
        }
    }
    if (j < candiEnd.size())
        candiEnd.resize(j);

    /*
    cerr << "candiStart" << endl << flush;
    for (i=0; i<candiStart.size(); i++)
        cout << candiStart[i] << endl << flush;
    cerr << "candiEnd" << endl << flush;
    for (i=0; i<candiEnd.size(); i++)
        cout << candiEnd[i] << endl << flush;

    exit(1);
     */
    
    /*
     // list out the potStart and potEnd
     cout << "potStart" << endl;
     for (i=0; i<potStart.size(); i++)
     cout << potStart[i] << endl;
     
     cout << "potEnd" << endl;
     for (i=0; i<potEnd.size(); i++)
     cout << potEnd[i] << endl;
     exit(1); */
    int preEnd = -1;
    for (i=0; i<candiEnd.size(); i++) {
        s = candiEnd[i];
        if (s <= preEnd)
            continue;
        for (j=0; j<candiStart.size(); j++) {
            t = candiStart[j];
            if (t > s && t-s+1 <= MAX_DROP_SIZE) {
                // if (i+1 >= potStart.size() || t < potStart[i+1] )
                dropRegions.push_back(pair<int,int>(s,t));
                preEnd = t;
                break;
            }
        }
    }
    
    // list out the drop regions
    cout << "Drop regions:" << endl;
    for (i=0; i<dropRegions.size(); i++) {
        cout << "[" << dropRegions[i].first << "," << dropRegions[i].second << "]" << endl;
    }
}

// change the bases to 'N' if its quality value is lower than the threshold
void Reads::checkQual(int thres) {
    // first guess the sanger format
    char highest_c = '!';
    char thres_c;
    int i,j;
    for (i=0; i<5000 && i<readlist.size(); i++) {
        // cout << readlist[i]->qualseq << endl;
        for (j=0; j<readlist[i]->qualseq.length(); j++) {
            if (readlist[i]->qualseq.at(j) > highest_c)
                highest_c = readlist[i]->qualseq.at(j);
        }
    }
    thres_c = highest_c - (40-thres);
    for (i=0; i<readlist.size(); i++) {
        for (j=0; j<readlist[i]->qualseq.length(); j++) {
            if (readlist[i]->seq[j] != '-' && readlist[i]->seq[j] != '_' && readlist[i]->qualseq.at(j) < thres_c){
                readlist[i]->seq[j] = 'N';
            }
        }
    }
}

Reads::~Reads() {
    // destructor
    int i;
    for (i=0; i<(int) readlist.size(); i++) {
        delete readlist[i];
    }
    readlist.clear();
}
