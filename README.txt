AFPhyloMix - Assembly-Free Phylogenetics for Mixtures

AFPhyloMix is a method to recover the phylogeny of haplotypes from short-read sequences obtained using pooled amplicons from a mixture of individuals, without barcoding. AFPhyloMix, accepts an alignment of the mixture of reads against a reference sequence, obtains the single-nucleotide-polymorphisms (SNP) patterns along the alignment, and constructs the phylogenetic tree according to the SNP patterns. AFPhyloMix adopts a Bayesian model of inference to estimates the phylogeny of the haplotypes and their relative frequencies, given that the number of haplotypes is known. 

INSTALLATION

The software was written in C++, and it has been tested under Linux and MacOS platform. You need to have a C++ compiler installed in the computer to compile the source codes. The compilation steps are shown as follows (xxx is the version of AFPhyloMix):

$ tar -zxvf AFPhyloMix-xxx.tar.gz
$ cd AFPhyloMix-xxx
$ make

Then an executable file: AFPhyloMix will appear.

SYNTAX

To run MCMCMC from SAM/BAM file

$ ./AFPhyloMix mcmcmc2 [SAM/BAM file] [# of haplos >= 2] [# of threads] <options>

To show the detailed information for each position

$ ./AFPhyloMix detail [file including all SAM file names]

To show the details of the potential snp positions

$ ./AFPhyloMix potentialsnp [SAM/BAM file]

Options

-r : resume MCMCMC from the previous run
-e [prefix of output file]
-n [number of generations] (default: 50k)
-p [print every number of cycles] (default: 100)
-b [check hot chains every number of cycles] (default: 5)
