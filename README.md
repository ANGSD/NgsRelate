[![Build Status](https://travis-ci.org/ANGSD/NgsRelate.svg?branch=master)](https://travis-ci.org/ANGSD/NgsRelate)


30juni 2018. Added new version which does analysis from bcf/vcf and outputs all nine jacquards

This page refers to the new v2 of NgsRelate which coestimates relatedness and inbreeding. For the old version please use this link: http://www.popgen.dk/software/index.php?title=NgsRelate&oldid=694

# Brief description #

This page contains information about the program called NgsRelate, which can be used to infer relatedness, inbreeding coefficients and many other summary statistics for pairs of individuals from low coverage Next Generation Sequencing (NGS) data by using genotype likelihoods instead of called genotypes. To be able to infer the relatedness you will need to know the population allele frequencies and have genotype likelihoods. This can be obtained e.g. using the program ANGSD as shown in example 1 below. For more information about ANGSD see here: http://popgen.dk/angsd/index.php/Quick_Start. As of version 2, VCF/BCF files can also be parsed.

# Options #
```
$ ./ngsRelate 

Usage main analyses: ./ngsrelate  [options] 
Options:
   -f <filename>       Name of file with frequencies
   -O <filename>       Output filename
   -L <INT>            Number of genomic sites. Must be provided if -f (allele frequency file) is NOT provided 
   -m <INTEGER>        model 0=normalEM 1=acceleratedEM
   -i <UINTEGER>       Maximum number of EM iterations
   -t <FLOAT>          Tolerance for breaking EM
   -r <FLOAT>          Seed for rand
   -R <chr:from-to>    Region for analysis (only for bcf)
   -g gfile            Name of genotypellh file
   -p <INT>            threads (default 4)
   -c <INT>            Should call genotypes instead?
   -s <INT>            Should you swich the freq with 1-freq?
   -F <INT>            Estimate inbreeding instead of estimating the nine jacquard coefficients
   -o <INT>            estimating the 3 jacquard coefficient, assumming no inbreeding
   -v <INT>            Verbose. print like per iteration
   -e <INT>            Errorrates when calling genotypes?
   -a <INT>            First individual used for analysis? (zero offset)
   -b <INT>            Second individual used for analysis? (zero offset)
   -B <INT>            Number of bootstrap replicates for (only for single pairs)
   -N <INT>            How many times to start each pair with random seed?
   -n <INT>            Number of samples in glf.gz
   -l <INT>            minMaf or 1-Maf filter
   -z <INT>            Name of file with IDs (optional)
   -T <STRING>         For -h vcf use PL (default) or GT tag
   -A <STRING>         For -h vcf use allele frequency TAG e.g. AFngsrelate (default)
   -P <filename>       plink name of the binary plink file (excluding the .bed)

Or
 ./ngsrelate extract_freq_bim pos.glf.gz plink.bim plink.freq
Or
 ./ngsrelate extract_freq .mafs.gz .pos.glf.gz [-rmTrans]
Or
 ./ngsrelate -h my.bcf [DEVELOPMENT ONLY]
```


# How to download and install #
On a linux or mac system with curl and g++ installed NgsRelate can be downloaded and installed as follows:
``` bash
git clone https://github.com/SAMtools/htslib
git clone https://github.com/ANGSD/ngsRelate
cd htslib/;make -j2;cd ../ngsRelate;make HTSSRC=../htslib/
```

# Run examples #

## Run example 1: using only NGS data ##
Below is an example of how NgsRelate can be used to coestimate relatedness and inbreeding from NGS data.
Assume we have file (`filelist`) containing paths to 100 BAM/CRAM files; one line per BAM/CRAM file. Then we can use ANGSD to estimate allele frequencies and calculate genotype likelihoods while doing SNP calling and in the end produce the the two input files needed for the NgsRelate program as follows:
``` bash
### First we generate a file with allele frequencies (angsdput.mafs.gz) and a file with genotype likelihoods (angsdput.glf.gz).
./angsd -b filelist -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3

### Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
zcat angsdput.mafs.gz | cut -f5 |sed 1d >freq

### run NgsRelate
./ngsrelate  -g angsdput.glf.gz -n 100 -f freq  -O newres
```
The output should be a file called `newres` that contains the output for all pairs between six individuals.

## Run example 2: using only VCF/BCF files ##
As of version 2, NgsRelate can parse BCF/VCF files using [htslib](https://github.com/SAMtools/htslib) with the following command:
``` bash
./ngsrelate  -h my.VCF.gz -O vcf.res
```
By default, NgsRelate will estimate the allele frequencies using the individuals provided in the VCF files. Allele frequencies from the INFO field can used be used instead using `-A TAG`. The `TAG` usually take the form of `AF` or `AF1` but can be set to anything. By default the *PL* data (Phred-scaled likelihoods for genotypes) is parsed, however, the called genotypes can also be used instead with the `-T GT` option.
If called genotypes are being used, the software requires an additional argument (`-c 1`). If using `-c 2`, ngsRelate calls genotypes assuming hardy-weinberg.

# Input file format #

## Genotype likelihood input ##
NgsRelate takes two files as input: a file with genotype likelihoods `(-g)` and a file with population allele frequencies `(-f)` for the sites there are genotype likelihoods for.
The genotype likelihood file needs to contain a line for each site with 3 values for each individual (one log transformed genotype likelihood for each of the 3 possible genotypes encoded as 'double's) and it needs to be in binary format and gz compressed.
The frequency file needs to contain a line per site with the allele frequency of the site in it.

NgsRelate also calculates a few summary statistics based on the 2dsfs instead of known population allele frequencies. Even if the population allele frequencies are not available, these summary statistics can still be estimated with NgsRelate by providing it with the number of sites `(-L)` instead of the allele frequency file `(-f)`. 

## VCF input ##
NgsRelate takes a standard VCF file generated with e.g. [bcftools](http://samtools.github.io/bcftools/). By default, NgsRelate will estimate the allele frequencies using the individuals provided in the VCF files. As shown in example 2, external allele frequencies can be provided with `-f`. These frequencies will overwrite the ones estimated from the VCF file.


# Output format #
``` bash
a  b  nSites  J9        J8        J7        J6        J5        J4        J3        J2        J1        rab       Fa        Fb        theta     inbred_relatedness_1_2  inbred_relatedness_2_1  fraternity  identity  zygosity  2of3IDB   FDiff      loglh           nIter  coverage  2dsfs                                                                             R0        R1        KING       2dsfs_loglike   2dsfsf_niter
0  1  99927   0.384487  0.360978  0.001416  0.178610  0.071681  0.000617  0.002172  0.000034  0.000005  0.237300  0.002828  0.250330  0.127884  0.001091                0.035846                0.001451    0.000005  0.001456  0.345411  -0.088997  -341223.335664  103    0.999270  0.154920,0.087526,0.038724,0.143087,0.155155,0.139345,0.038473,0.087632,0.155138  0.497548  0.290124  0.000991   -356967.175857  7
```

The first two columns contain indices of the two individuals used for the analysis. The third column is the number of genomic sites considered. The following nine columns are the maximum likelihood (ML) estimates of the nine jacquard coefficients, where  K0==J9; K1==J8; K2==J7 in absence of inbreeding. 
Based on these Jacquard coefficients, NgsRelate calculates 11 summary statistics:

13. rab is the pairwise relatedness `(J1+J7+0.75*(J3+J5)+.5*J8)` [Hedrick et al](https://academic.oup.com/jhered/article/106/1/20/2961876)
14. Fa is the inbreeding coefficient of individual a `J1+J2+J3+J4` [Jacquard](https://www.springer.com/gp/book/9783642884177)
15. Fb is the inbreeding coefficient of individual b `J1+J2+J5+J6` [Jacquard](https://www.springer.com/gp/book/9783642884177)
16. theta is the coefficient of kinship `J1 + 0.5*(J3+J5+J7) + 0.25*J8)` [Jacquard](https://www.springer.com/gp/book/9783642884177)
17. inbred_relatedness_1_2 `J1+0.5*J3` [Ackerman et al](http://www.genetics.org/content/206/1/105)
18. inbred_relatedness_2_1 `J1+0.5*J5` [Ackerman et al](http://www.genetics.org/content/206/1/105)
19. fraternity `J2+J7` [Ackerman et al](http://www.genetics.org/content/206/1/105)
20. identity `J1` [Ackerman et al](http://www.genetics.org/content/206/1/105) 
21. zygosity `J1+J2+J7` [Ackerman et al](http://www.genetics.org/content/206/1/105)
22. Two-out-of-three IBD `J1+J2+J3+J5+J7+0.5*(J4+J6+J8)` [Miklos csuros](https://www.sciencedirect.com/science/article/pii/S0040580913001123)
23. Inbreeding difference `0.5*(J4-J6)` [Miklos csuros](https://www.sciencedirect.com/science/article/pii/S0040580913001123)
24. the log-likelihood of the ML estimate. 
25. number of EM iterations. If a `-1` is displayed. A boundary estimate had a higher likelihood. 
26. If differs from `-1`, a boundary estimate had a higher likelihood. Reported loglikelihood should be highly similar to the corresponding value reported in `loglh`
27. fraction of sites used for the ML estimate

The remaining columns relate to statistics based on a 2D-SFS. 

28. 2dsfs estimates using the same methodology as implemented in realSFS, see [ANGSD](https://github.com/ANGSD/angsd)
29. R0 [Waples et al](https://www.biorxiv.org/content/early/2018/08/31/260497)
30. R1 [Waples et al](https://www.biorxiv.org/content/early/2018/08/31/260497)
31. KING [Waples et al](https://www.biorxiv.org/content/early/2018/08/31/260497)
32. the log-likelihood of the 2dsfs estimate.
33. number of iterations for 2dsfs estimate

You can also input a file with the IDs of the individuals (on ID per line), using the `-z` option, in the same order as in the file `filelist` used to make the genotype likelihoods or the VCF file. If you do the output will also contain these IDs as column 3 and 4.

Note that in some cases nIter is -1. This indicates that values on the boundary of the parameter space had a higher likelihood than the values achieved using the EM-algorithm (ML methods sometimes have trouble finding the ML estimate when it is on the boundary of the parameter space, and we therefore test the boundary values explicitly and output these if these have the highest likelihood).

# Help and additional options #

To get help and a list of all options simply type

``` bash
./ngsrelate
```

Note that for the new version of ngsRelate it is not necessary to flip the allele frequencies if the input allele frequencies are for the minor allele.

## running version 1 of NgsRelate ##
To run the old version (V1) that assumes both individuals to be out-bred, just add `-o 1` to your command. This will force Jacqaurd coefficient 1-6 (the coefficients related to inbreeding) to zero.
``` bash
./ngsrelate -o 1
```

## Test for inbreeding only ##
To only estimate inbreeding, use the following command
``` bash
./ngsrelate -F 1
```


# Citing and references #

- Method (v1) is published here: http://bioinformatics.oxfordjournals.org/content/early/2015/08/29/bioinformatics.btv509.abstract
- Method (v2) is published here:
https://academic.oup.com/gigascience/article/8/5/giz034/5481763

# Changelog #

# Contacts #
thorfinn@binf.ku.dk, k.hanghoej@snm.ku.dk, and ida@binf.ku.dk

# Testdataset
```bash
bcftools call -v -m test.bcf -Ob >small.bcf
bcftools view small.bcf -Oz -o small.vcf.gz -m2 -M2 -v snps
plink --vcf small.vcf.gz --make-bed --out small 
bgzip small.bed
md5sum small.* >small.md5
```

# check glf file in `R`

```R
a <- file("test.glf", "rb")
nobs <- 1e8 # large number 3*N*Nsites
N <- 6 # number of individuals
mat <- matrix(readBin(a, "double", nobs ), byrow=T, ncol=3*N)
```
