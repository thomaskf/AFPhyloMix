# AFPhyloMix

Background

Building a phylogenetic tree is an essential step in many researches in evolutionary biology. Methods of molecular phylogenetic reconstruction always requires the input of genomic sequences of individuals. Each genomic sequences can be produced through the assembly process on the reads belonging to the same individual. The phylogenetic reconstruction methods usually rely on the sequencing and the assembly processes of reads from the same individual. However, in some cases, the reads are mixed between multiple individuals and cannot be separated easily. For example, it is challenging to infer the phylogenetic tree for the rapidly evolving viruses like Human Immunodeficiency Virus (HIV) inside a host, because the viruses inside the host cannot be easily divided into isolated individual virus.

Method

AFPhyloMix (Assembly-Free Phylogenetics for Mixtures) is developed to build a phylogenetic tree directly from the mixture of reads which cannot be separated to different individuals, without the assembly process of genomic sequences. AFPhyloMix feeds in the alignments of the mixture of reads against a reference sequence, obtains the single-nucleotide-polymorphisms (SNP) patterns along the genomic positions, and constructs the phylogenetic tree according to the SNP patterns. AFPhyloMix adopts a Bayesian model of inference to estimates the phylogeny of the haplotypes and their relative frequencies, given that the number of haplotypes is known. 
