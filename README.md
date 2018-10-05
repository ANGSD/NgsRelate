[![Build Status](https://travis-ci.org/ANGSD/NgsRelate.svg?branch=master)](https://travis-ci.org/ANGSD/NgsRelate)
30juni 2018. Added new version which does analysis from bcf/vcf and outputs all nine jacquards

This pages refers to the new v2 of ngsRelate which coestimates relatedness and inbreeding. For the old version please use this link: http://www.popgen.dk/software/index.php?title=NgsRelate&oldid=694

# Brief description #

This page contains information about the program called NgsRelate, which can be used to infer relatedness coefficients for pairs of individuals from low coverage Next Generation Sequencing (NGS) data by using genotype likelihoods instead of called genotypes. To be able to infer the relatedness you will need to know the population allele frequencies and have genotype likelihoods. This can be obtained e.g. using the program ANGSD as shown in the examples below. For more information about ANGSD see here: http://popgen.dk/angsd/index.php/Quick_Start.

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
Assume we have file containing paths to 100 BAM/CRAM files; one line per BAN/CRAM file. Then we can use ANGSD to estimate allele frequencies and calculate genotype likelihoods while doing SNP calling and in the end produce the the two input files needed for the NgsRelate program as follows:
``` bash
### First we generate a file with allele frequencies (angsdput.mafs.gz) and a file with genotype likelihoods (angsdput.glf.gz).
./angsd -b filelist -gl 1 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3

### Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
zcat angsdput.mafs.gz | cut -f5 |sed 1d >freq

### run NgsRelate
./ngsrelate  -g angsdput.glf.gz -n 6 -f freq  > newres
```
The output should be a file called `newres` that contains relatedness estimates for all pairs between six individuals.

## Run example 2: using only VCF/BCF files ##
As of version 2, NgsRelate can parse BCF/VCF files using `htslib` with the following command:
``` bash
./ngsrelate  -h my.VCF.gz > newres
```
By default, NgsRelate will estimate the allele frequencies using the individuals provided in the VCF files. External allele frequencies can be provided with the following command:
``` bash
./ngsrelate  -h my.VCF.gz -f freq > newres
```
The `freq` file should contain one allele frequency per line as shown in example 1. NOTE: The end-user must make sure that the allele frequencies overlap the sites provided in the VCF file.

# Input file format #

## Genotype likelihood input ##
NgsRelate takes two files as input: a file with genotype likelihoods and a file with frequencies for the sites there are genotype likelihoods for.
The genotype likelihood file needs to contain a line for each site with 3 values for each individual (one log transformed genotype likelihood for each of the 3 possible genotypes encoded as 'double's) and it needs to be in binary format and gz compressed.
The frequency file needs to contain a line per site with the allele frequency of the site in it.

## VCF input ##
NgsRelate takes a standard VCF file generated with e.g. (bcftools)[http://samtools.github.io/bcftools/]. By default, NgsRelate will estimate the allele frequencies using the individuals provided in the VCF files. As shown above, external allele frequencies can be provided.


# Output format #
COMING SOON

# Help and additional options #

To get help and a list of all options simply type

``` bash
./ngsrelate
```

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

Method (v1) is published here: http://bioinformatics.oxfordjournals.org/content/early/2015/08/29/bioinformatics.btv509.abstract
Method (v2) is currently under review.

# Changelog #
Important recent changes:
#The option -z has been added so one can get the sample IDs printed in the output (if one run the program with -z idfilename)
#We have fixed -m 1 so the estimates can no longer be negative

See github for the full change log.

=Bugs/Improvements=
-Make better output message if files doesn't exists when using the extract_freq option

# Contacts #
thorfinn@binf.ku.dk, k.hanghoej@snm.ku.dk, and ida@binf.ku.dk
