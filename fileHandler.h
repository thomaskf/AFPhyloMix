#ifndef _FILEHANDER_
#define _FILEHANDER_

#include<string>
#include<iostream>
#include<fstream>
#include<zlib.h>
#include "samtools-0.1.18/bam.h"

using namespace std;

class SamBamFileHander {
public:
	int type; // 1: BAM; 2: SAM
	void openFile(char* fileName);
	void closeFile();
	bool getNextSeq(string& line);
	
private:
	// for BAM file
	bamFile bamQueryFile;
	bam_header_t * bamHeader;
	bam1_t * bam;
	// char* rec_seq;
	// for SAM file
	ifstream fin;
};

#endif
