//
//  constructHaplo.h
//  HaploCount
//
//  Created by Thomas Wong on 10/6/19.
//

#ifndef constructHaplo_h
#define constructHaplo_h

#include <stdio.h>
#include <iomanip>
#include <set>
#include "reads.h"
#include "tree.h"
#include "gtree.h"
#include "tree.h"
#include "phyloRead.h"

class Haplotypes {
public:
    
    PhyloRead phyloRead;

    // load the max-tree file
    void loadMaxTreeFile(char* maxTreeFile);

    // load samfile
    void load_sam(char* samFile);
    
    // construct the haplotypes based on the connectivity of the columns inside reads
    void construct(char* samFile, char* maxTreeFile, char* outFile);

private:
    
    // output the haplotypes
    void outHaplo(char* outFile, vector<string>& haplos, vector<double>& freqs);
};

#endif /* constructHaplo_hpp */
