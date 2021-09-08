#include "reads.h"
#include "tree.h"
#include "phyloMCMC.h"
#include "constructHaplo.h"
#include "mcmcmc.h"
#include "mylib.h"
#include <math.h>
#include <time.h>
#include <fstream>

#define VERSION "c1.10"

void showHelpMenu(char** argv) {
    cerr << "AFPhyloMix (Assembly-Free Phylogenetics for Mixtures) is a method to recover the phylogeny of haplotypes from short-read sequences obtained using pooled amplicons from a mixture of individuals, without barcoding." << endl;
    cerr << endl;
    cerr << "Please note that AFPhyloMix only works on the reads sequenced under the same sequencing run." << endl;
    cerr << endl;
    cerr << "Syntax:" << endl;
    cerr << "To run MCMCMC from SAM/BAM file" << endl;
    cerr << "  " << argv[0] << " mcmcmc2 [SAM/BAM file] [# of haplos >= 2] [# of threads] <options>" << endl;
    cerr << endl;
    
    // cerr << "To show the details of the potential snp positions" << endl;
    // cerr << "  " << argv[0] << " potentialsnp [SAM/BAM file]" << endl << endl;

    cerr << "Options:" << endl;
    cerr << "-r : resume MCMCMC from the previous run" << endl;
    cerr << "-e [prefix of output file]" << endl;
    cerr << "-n [number of generations] (default: 50k)" << endl;
    cerr << "-p [print every number of cycles] (default: 100)" << endl;
    cerr << "-b [check hot chains every number of cycles] (default: 5)" << endl;
    cerr << endl;

    cerr << "Output fles:" << endl;
    cerr << "[SAM/BAM file w/o ext].max.tree : the tree with tip frequencies having the maximum posterior probability along the cold chain" << endl;
    cerr << "[SAM/BAM file w/o ext].x.log    : the log file listing the likelihod, posterior probability, and the parameter values every, by default, 100 cycles along the x-th chain" << endl;
    cerr << endl;

    cerr << "Example (# of tips = 5, # of CPU threads = 8):" << endl;
    cerr << "  " << argv[0] << " mcmcmc2 example.bam 5 8" << endl;
    cerr << endl;
    cerr << "The true tree is: example.real.tree" << endl;
    cerr << "You may compare the resulting estimated tree - example.max.tree with the true tree" << endl;
    cerr << endl;
}

void showDetailHelpMenu(char** argv) {
    cerr << "AFPhyloMix (Assembly-Free Phylogenetics for Mixtures) is developed to build a phylogenetic tree directly from the mixture of reads which cannot be separated to different individuals, without the assembly process of genomic sequences." << endl;
    cerr << endl;
    cerr << "Syntax:" << endl;
    cerr << "To run order-1 MCMC from SAM/BAM file" << endl;
    cerr << "  " << argv[0] << " mcmc1 [SAM/BAM file] [# of haplos >= 2] <options>" << endl;
    cerr << "To run order-2 MCMC from SAM/BAM file" << endl;
    cerr << "  " << argv[0] << " mcmc2 [SAM/BAM file] [# of haplos >= 2] <options>" << endl;
    cerr << "To run order-1 MCMCMC from SAM/BAM file" << endl;
    cerr << "  " << argv[0] << " mcmcmc1 [SAM/BAM file] [# of haplos >= 2] [# of threads] <options>" << endl;
    cerr << "To run order-2 MCMCMC from SAM/BAM file" << endl;
    cerr << "  " << argv[0] << " mcmcmc2 [SAM/BAM file] [# of haplos >= 2] [# of threads] <options>" << endl;
    cerr << "To constuct haplotypes from mcmc results" << endl;
    cerr << "  " << argv[0] << " haplo [SAM/BAM file] [mcmc max tree file] [haplo output file]" << endl;
    cerr << "To get the snp information from the alignment" << endl;
    cerr << "  " << argv[0] << " snp [SAM/BAM file]" << endl;
    cerr << "To obtain the read coverage information" << endl;
    cerr << "  " << argv[0] << " cover [file including all SAM file names]" << endl << endl;
    
    cerr << "To show the detailed information for each position" << endl;
    cerr << "  " << argv[0] << " detail [file including all SAM file names]" << endl;
    cerr << "To identify the dropping region" << endl;
    cerr << "  " << argv[0] << " dropreg [SAM/BAM file]" << endl;
    cerr << "To show the details of the potential snp positions" << endl;
    cerr << "  " << argv[0] << " potentialsnp [SAM/BAM file]" << endl << endl;

    cerr << "Options:" << endl;
    cerr << "-r : resume MCMC from the previous run (cannot use it with -t)" << endl;
    cerr << "-t [start tree file] : specify the starting tree for MCMC run (cannot use it with -r or with MCMCMC)" << endl;
    cerr << "-c [show position]" << endl;
    cerr << "-m [back mutation method]" << endl;
    cerr << "-e [prefix of output file]" << endl;
    cerr << "-n [number of generations] (default: 50k)" << endl;
    cerr << "-p [print every number of cycles] (default: 100)" << endl;
    cerr << "-b [check hot chains every number of cycles] (default: 5)" << endl;
    
    cerr << endl;

    // for examining the efficiency of the back-mutation algorithm
    cerr << "To report the actual position with back mutations:" << endl;
    cerr << "  " << argv[0] << " backmutate [SAM/BAM file] [tree file with tip names A, B, C, ...]" << endl;
    
}

void readOtherOptions(int argc, char** argv, int startIndex, int& toResume, char* startTreeFile, set<int>& show_pos, int& backMutateMethod, int& nGenerations, int& nPrintSteps, string& prefix, int& nCheckChains) {
    
    int i,j;
    vector<string> token;
    show_pos.clear();
    for (i=startIndex; i<argc; i++) {
        if (strcmp(argv[i],"-r")==0) {
            toResume = 1;
            cout << "Resume from the previous run" << endl;
        } else if (strcmp(argv[i],"-t")==0 && i+1<argc) {
            startTreeFile = argv[i+1];
            cout << "Starting tree file: " << startTreeFile << endl;
            i++;
        } else if (strcmp(argv[i],"-c")==0 && i+1<argc) {
            tokenizer(argv[i+1], ",", &token);
            for (j=0; j<token.size(); j++)
                show_pos.insert(atoi(token[j].c_str()));
            cout << "Showing position: " << argv[i+1] << endl;
            i++;
        } else if (strcmp(argv[i],"-e")==0 && i+1<argc) {
            prefix = argv[i+1];
            i++;
        } else if (strcmp(argv[i],"-m")==0 && i+1<argc) {
            backMutateMethod = atoi(argv[i+1]);
            i++;
        } else if (strcmp(argv[i],"-n")==0 && i+1<argc) {
            nGenerations = atoi(argv[i+1]);
            i++;
        } else if (strcmp(argv[i],"-p")==0 && i+1<argc) {
            nPrintSteps = atoi(argv[i+1]);
            i++;
        } else if (strcmp(argv[i],"-b")==0 && i+1<argc) {
            nCheckChains = atoi(argv[i+1]);
            i++;
        }
    }
    // print out information
    cout << "Prefix: " << prefix << endl;
    cout << "Back mutation detection method: " << backMutateMethod << endl;
    cout << "Number of generations: " << nGenerations << endl;
    cout << "Print every number of steps: " << nPrintSteps << endl;
    cout << "Check different chain every number of steps (for MCMCMC): " << nCheckChains << endl;
}

int main(int argc, char** argv) {
    
	// initialize the clock
	clock_t start_t = clock();
	 
	// display version number
	cout << "Version " << VERSION << endl;
    // cout << "This version considers substrings with two snp positions" << endl;
	cout << endl;
	
    // check the syntax
    if (argc < 3) {
        showHelpMenu(argv);
        exit(1);
    }

    PhyloMCMC phyloReadMCMC;
    Haplotypes haplo;
    
    char* option = argv[1];
    char* inputFile = argv[2];
    // char* preLogFile = NULL;
    char* maxTreeFile = NULL;
    char* outHaploFile = NULL;
    char* startTreeFile = NULL;
    int nGenerations = 50000;
    int nPrintSteps = 100;
    int nCheckChains = 5;
    char int2nucl[] = {'_','A','C','G','T'};
    int numHap;
    int toConsiderAdj;
    int i;
    int startOptionIndex;
    set<int> show_pos;
    int backMutateMethod = 3; // default is 3
    int toResume = 0;
    int numThreads;

    // get the prefix of input file
    string inputFileName = inputFile;
    string prefix = inputFileName;
    i = inputFileName.find_last_of('.');
    if (i > 1)
        prefix = inputFileName.substr(0,i);
    
    /*
    // -------------------------
    // for checking......
    
    Reads reads;
    vector<string> fileNames;
    vector<vector<int> > coverages;
    vector<bool> problematic; // working
    ifstream fin;
    string aline;
    int tot_cover;
    double minSnpRatio = 0.05;
    int max_dim;
    int dim;
    int* posMatrix;
    int j, k, n, l;
    coverages.clear();
    fileNames.clear();
    fin.open(inputFile);
    k=0;
    max_dim = 0;
    while (getline(fin, aline)) {
        if (aline.length() > 0) {
            fileNames.push_back(aline);
            reads.readSamFile((char *)aline.c_str());
            posMatrix = reads.getPosMatrix(dim);
            if (max_dim < dim)
                max_dim = dim;
            problematic.resize(max_dim, false);
            // check whether there exists a problem
            for (i=0; i<dim; i++) {
                tot_cover = 0;
                for (j=0; j<5; j++) {
                    tot_cover += posMatrix[i*5 + j];
                }
                l=0;
                for (j=0; j<5; j++) {
                    if ((double)posMatrix[i*5 + j] / tot_cover >= minSnpRatio) {
                        l++;
                    }
                }
                if (l>=2) {
                    // reads from the same bar codes should not have snp
                    problematic[i] = true;
                }
            }
            // get the number of reads covering each position
            coverages.resize(k+5);
            coverages[k].clear();
            for (j=0; j<5; j++) {
                for (i=0; i<dim; i++) {
                    coverages[k].push_back(posMatrix[i*5 + j]);
                }
                k++;
            }
            delete(posMatrix);
        }
    }
    fin.close();
    // print out the details
    n = fileNames.size();
    cout << "pos";
    for (i=0; i<n; i++) {
        for (j=0; j<5; j++) {
            if (j == 0)
                cout << "," << fileNames[i];
            else
                cout << ",";
        }
    }
    cout << ",isProblematic" << endl;
    cout << ",";
    for (i=0; i<n; i++) {
        for (j=0; j<5; j++) {
            cout << "," << int2nucl[j];
        }
    }
    cout << "," << endl;
    for (j=0; j<max_dim; j++) {
        cout << j+1;
        for (i=0; i<5*n; i++) {
            if (j < coverages[i].size()) {
                cout << "," << coverages[i].at(j);
            } else {
                cout << ",0";
            }
        }
        cout << "," << (int) problematic[j];
        cout << endl;
    }
    
    // -------------------------
    exit(1);
    */
    
    if (strcmp(option,"mcmc1") == 0) {
        numHap = atoi(argv[3]);
        startOptionIndex = 4;
        readOtherOptions(argc, argv, startOptionIndex, toResume, startTreeFile, show_pos, backMutateMethod, nGenerations, nPrintSteps, prefix, nCheckChains);
        toConsiderAdj = 0; // not considering adjacent SNP positions
        phyloReadMCMC.run(inputFile, numHap, nGenerations, nPrintSteps, toConsiderAdj, startTreeFile, show_pos, toResume, backMutateMethod, prefix);
    } else if (strcmp(option,"mcmc2") == 0) {
        numHap = atoi(argv[3]);
        startOptionIndex = 4;
        readOtherOptions(argc, argv, startOptionIndex, toResume, startTreeFile, show_pos, backMutateMethod, nGenerations, nPrintSteps, prefix, nCheckChains);
        toConsiderAdj = 1; // consider adjacent SNP positions
        phyloReadMCMC.run(inputFile, numHap, nGenerations, nPrintSteps, toConsiderAdj, startTreeFile, show_pos, toResume, backMutateMethod, prefix);
    } else if (strcmp(option, "mcmcmc1") == 0) {
        numHap = atoi(argv[3]);
        numThreads = atoi(argv[4]);
        startOptionIndex = 5;
        MCMCMC phyloReadMCMCMC(numThreads);
        readOtherOptions(argc, argv, startOptionIndex, toResume, startTreeFile, show_pos, backMutateMethod, nGenerations, nPrintSteps, prefix, nCheckChains);
        toConsiderAdj = 0; // consider adjacent SNP positions
        phyloReadMCMCMC.run(inputFile, numHap, nGenerations, nPrintSteps, toConsiderAdj, startTreeFile, show_pos, toResume, backMutateMethod, prefix, nCheckChains);
    } else if (strcmp(option, "mcmcmc2") == 0) {
        numHap = atoi(argv[3]);
        numThreads = atoi(argv[4]);
        startOptionIndex = 5;
        MCMCMC phyloReadMCMCMC(numThreads);
        readOtherOptions(argc, argv, startOptionIndex, toResume, startTreeFile, show_pos, backMutateMethod, nGenerations, nPrintSteps, prefix, nCheckChains);
        toConsiderAdj = 1; // consider adjacent SNP positions
        phyloReadMCMCMC.run(inputFile, numHap, nGenerations, nPrintSteps, toConsiderAdj, startTreeFile, show_pos, toResume, backMutateMethod, prefix, nCheckChains);
    } else if (strcmp(option, "haplo") == 0 && argc == 5) {
        maxTreeFile = argv[3];
        outHaploFile = argv[4];
        haplo.construct(inputFile, maxTreeFile, outHaploFile);
    } else if (strcmp(option, "snp") == 0) {
        double avgCover;
        int numSNPs, seqLen;
        // get the information of the alignment
        phyloReadMCMC.get_align_info(inputFile, avgCover, numSNPs, seqLen);
        cout << "Average coverage\tNumber of SNPs\tSequence length" << endl;
        cout << avgCover << "\t" << numSNPs << "\t" << seqLen << endl;
    } else if (strcmp(option, "backmutate") == 0) {
        RootTree rtree;
        vector<vector<int> > c_sets;
        PhyloRead phyloRead;
        int numHap;
        char* treeFile = argv[3];
        rtree.get_mutation_sets(treeFile, c_sets, numHap);
        phyloRead.getActualBackMutate(inputFile, c_sets, numHap);
    } else if (strcmp(option, "cover") == 0) {
        Reads reads;
        vector<string> fileNames;
        vector<vector<int> > coverages;
        ifstream fin;
        string aline;
        int max_dim;
        int dim;
        int* posMatrix;
        int i, j, k, n, cover;
        coverages.clear();
        fileNames.clear();
        fin.open(inputFile);
        k=0;
        max_dim = 0;
        while (getline(fin, aline)) {
            if (aline.length() > 0) {
                fileNames.push_back(aline);
                reads.readSamFile((char *)aline.c_str());
                posMatrix = reads.getPosMatrix(dim);
                if (max_dim < dim)
                    max_dim = dim;
                // get the number of reads covering each position
                coverages.resize(k+1);
                coverages[k].clear();
                for (i=0; i<dim; i++) {
                    cover = 0;
                    for (j=0; j<5; j++) {
                        cover += posMatrix[i*5 + j];
                    }
                    coverages[k].push_back(cover);
                }
                delete(posMatrix);
                k++;
            }
        }
        fin.close();
        // print out all the read coverages
        n = fileNames.size();
        cout << "pos";
        for (i=0; i<n; i++) {
            cout << " " << fileNames[i];
        }
        cout << endl;
        for (j=0; j<max_dim; j++) {
            cout << j+1;
            for (i=0; i<n; i++) {
                if (j < coverages[i].size()) {
                    cout << " " << coverages[i].at(j);
                } else {
                    cout << " -";
                }
            }
            cout << endl;
        }
    } else if (strcmp(option, "detail") == 0) {
        Reads reads;
        vector<string> fileNames;
        vector<vector<int> > coverages;
        ifstream fin;
        string aline;
        int max_dim;
        int dim;
        int* posMatrix;
        int i, j, k, n;
        coverages.clear();
        fileNames.clear();
        fin.open(inputFile);
        k=0;
        max_dim = 0;
        while (getline(fin, aline)) {
            if (aline.length() > 0) {
                fileNames.push_back(aline);
                reads.readSamFile((char *)aline.c_str());
                posMatrix = reads.getPosMatrix(dim);
                if (max_dim < dim)
                    max_dim = dim;
                // get the number of reads covering each position
                coverages.resize(k+5);
                coverages[k].clear();
                for (j=0; j<5; j++) {
                    for (i=0; i<dim; i++) {
                        coverages[k].push_back(posMatrix[i*5 + j]);
                    }
                    k++;
                }
                delete(posMatrix);
            }
        }
        fin.close();
        // print out the details
        n = fileNames.size();
        cout << "pos";
        for (i=0; i<n; i++) {
            for (j=0; j<5; j++) {
                if (j == 0)
                    cout << " " << fileNames[i];
                else
                    cout << " ";
            }
        }
        cout << endl;
        cout << " ";
        for (i=0; i<n; i++) {
            for (j=0; j<5; j++) {
                cout << " " << int2nucl[j];
            }
        }
        cout << endl;
        for (j=0; j<max_dim; j++) {
            cout << j+1;
            for (i=0; i<5*n; i++) {
                if (j < coverages[i].size()) {
                    cout << " " << coverages[i].at(j);
                } else {
                    cout << " 0";
                }
            }
            cout << endl;
        }
    } else if (strcmp(option, "potentialsnp") == 0) {
        Reads reads;
        int* posMatrix;
        int c, i, j, k, n;
        double r;
        double minSnpRatio = 0.05;
        int minCover = 1000;
        reads.readSamFile(inputFile);
        posMatrix = reads.getPosMatrix(n);
        cout << "pos - A C G T" << endl;
        for (i=0; i<n; i++) {
            c=0;
            for (j=0; j<5; j++) {
                c+=posMatrix[i*5 + j];
            }
            if (c < minCover)
                continue;
            if ((double) posMatrix[i*5] / c > minSnpRatio)
                continue; // too many gaps
            k=0;
            for (j=1; j<5; j++) {
                r = (double) posMatrix[i*5 + j] / c;
                if (r >= minSnpRatio)
                    k++;
            }   
            if (k>1) {
                // potential a snp position
                cout << i+1;
                for (j=0; j<5; j++) {
                    cout << " " << posMatrix[i*5 + j];
                }
                cout << endl;
            }
        }
    } else if (strcmp(option, "dropreg") == 0) {
        // read the sam/bam file
        Reads rds;
        int numcol;
        int* posmatrix;
        
        cerr << "Processing the sam/bam file..." << endl << flush;
        rds.readSamFile(inputFile);
        
        posmatrix = rds.getPosMatrix(numcol);
        
        // cerr << "Performing base error correction...." << endl << flush;
        // rds.baseCorrect(MINCOVER, MINRATIO);
        
        // cerr << "Performing kmer error correction...." << endl << flush;
        // rds.kmerCorrect(WINLEN, MINCOVER, MINFREQ, MINWINRATIO);
        
        vector<pair<int,int> > dropRegions;
        cerr << "Identifying the dropping regions...." << endl << flush;
        rds.getCoverDropRegs(dropRegions, numcol, posmatrix);
        
        delete[] posmatrix;

    } else {
        cerr << "Error! Unknown option: " << argv[1] << endl << endl;
        showHelpMenu(argv);
        exit(1);
    }
    
    cout << "finish!" << endl;
	
	// get the total time elapsed
	clock_t end_t = clock();
	timeElapsed(start_t, end_t, "total time elapsed");
	
}
