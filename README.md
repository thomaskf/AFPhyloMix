# AFPhyloMix: Assembly-Free Phylogenetics for Mixtures

AFPhyloMix is a method to recover the phylogeny of haplotypes from short-read sequences obtained using pooled amplicons from a mixture of individuals, without barcoding. AFPhyloMix, accepts an alignment of the mixture of reads against a reference sequence, obtains the single-nucleotide-polymorphisms (SNP) patterns along the alignment, and constructs the phylogenetic tree according to the SNP patterns. AFPhyloMix adopts a Bayesian model of inference to estimates the phylogeny of the haplotypes and their relative frequencies, given that the number of haplotypes is known.

AFPhyloMix estimates the phylogeny of the haplotypes and their relative abundances, <b>without the need for reconstructing the haplotype sequences</b>.

Please note that AFPhyloMix only works on the reads sequenced under <b>the same sequencing run</b>. When all the haplotypes are sequenced under the same run (i.e. all haplotypes were pooled into the same library before sequencing), the read coverages of the haplotypes would have similar trends: all have relatively high (or low) read coverages at the same regions of the genome. The similar trends of the read coverages along the genome between the haplotypes led to a consistent distribution of the ratios of haplotypes along the genome. The consistent read coverage ratios along the genome worked in our method's favor.

Moreover, the assumption of an infinite sites model that is applied in AFPhyloMix is appropriate for <b>closely-related individuals</b>. To extend our algorithm to more divergent sequences will require a different model of mutation. This remains a work in progress.

*Installation*

The software was written in C++, and it has been tested under Linux and MacOS platform. You need
to have a C++ compiler installed in the computer to compile the source codes. The compilation steps
are shown as follows (xxx is the version of AFPhyloMix):

```
$ tar -zxvf AFPhyloMix-xxx.tar.gz
$ cd AFPhyloMix-xxx
$ make
```

Then an executable file: AFPhyloMix will appear. 

*Syntax*

To run MCMCMC from SAM/BAM file
```
$ ./AFPhyloMix mcmcmc2 [SAM/BAM file] [# of haplos >= 2] [# of threads] <options>
```

Options
```
-r : resume MCMCMC from the previous run
-e [prefix of output file]
-n [number of generations] (default: 50k)
-p [print every number of cycles] (default: 100)
-b [check hot chains every number of cycles] (default: 5)
```

Output files
```
[SAM/BAM file w/o ext].max.tree : the tree with tip frequencies having the maximum posterior probability along the cold chain
[SAM/BAM file w/o ext].x.log    : the log file listing the likelihod, posterior probability, and the parameter values every, by default, 100 cycles along the x-th chain
```

Example (# of tips = 5, # of CPU threads = 8)
```
$ ./AFPhyloMix mcmcmc2 example.bam 5 8
```

The true tree is: example.real.tree
You may compare the resulting estimated tree - example.max.tree with the true tree
