#!/bin/bash


DNAME=$1

for i in `seq 20`;do sed -n '2p' ${DNAME}/fc$i.beagle.genome;done >$DNAME/fc.genomes
for i in `seq 20`;do sed -n '2p' ${DNAME}/sc$i.beagle.genome;done >${DNAME}/sc.genomes
for i in `seq 20`;do cat $DNAME/sc$i.beagle.res;done >${DNAME}/sc.ml
for i in `seq 20`;do cat ${DNAME}/fc$i.beagle.res;done >${DNAME}/fc.ml