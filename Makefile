CC = g++
SAMLIB = samtools-0.1.18
SAMOBJLIBS = $(SAMLIB)/sam.o $(SAMLIB)/bam.o $(SAMLIB)/bgzf.o $(SAMLIB)/kstring.o $(SAMLIB)/bam_import.o $(SAMLIB)/faidx.o $(SAMLIB)/bam_pileup.o $(SAMLIB)/bam_aux.o $(SAMLIB)/sam_header.o $(SAMLIB)/razf.o
CFLAGS = -O3 -std=c++11
ZFLAGS = -lz -lpthread

all : AFPhyloMix

AFPhyloMix : myRand.h myRand.cpp mcmc.h mcmc.cpp phyloMCMC.h phyloMCMC.cpp gtree.cpp gtree.h fileHandler.cpp fileHandler.h main.cpp mylib.cpp mylib.h reads.cpp reads.h tree.cpp tree.h phyloRead.cpp phyloRead.h constructHaplo.h constructHaplo.cpp mcmcmc.cpp mcmcmc.h definitions.h $(SAMOBJLIBS)
	g++ -o AFPhyloMix myRand.cpp mcmc.cpp phyloMCMC.cpp gtree.cpp fileHandler.cpp main.cpp mylib.cpp reads.cpp tree.cpp phyloRead.cpp constructHaplo.cpp mcmcmc.cpp $(SAMOBJLIBS) $(ZFLAGS) $(CFLAGS)


clean :
	rm -rf AFPhyloMix $(SAMOBJLIBS)
