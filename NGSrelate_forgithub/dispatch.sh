#!/bin/bash

FILE=$1

#run NGSrelate
./NGSrelate -beagle $FILE -k2 0

#make angsd input
./beagleToTglfs $FILE  $FILE

#make plink input
./angsd0.551/angsd -sim1 ${FILE}.glf.gz -nInd 52 -out $FILE -doPlink 2 -doGeno -1 -doPost 1 -doMajorMinor 1 -doMaf 1

#make --genome analysis with plink
plink --noweb --tfile $FILE --genome --out $FILE
