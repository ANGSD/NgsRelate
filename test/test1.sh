#!/bin/bash
MD5=md5sum

if [ $# -eq 1 ] 
then
    PRG=$1
fi

if [[ "$OSTYPE" == "darwin"* ]]; then
    MD5=./md5osx.sh
fi


LOG=${0}.log
ODIR=odir
echo "Cleaning old output dir ${ODIR} PRG:${PRG}" >${LOG}
rm -rf ${ODIR} ${LOG}
mkdir -p ${ODIR}

gunzip -c small.bed.gz >small.bed

echo "1) part1" >>${LOG}
${PRG} -h small.vcf.gz -T GT -p 1 -c 1 -r 100 >${ODIR}/vcf.res 2>>${LOG} 
${PRG} -h small.bcf -T GT -p 1 -c 1 -r 100 >${ODIR}/vcf2.res 2>>${LOG} 


echo -e "\t test all done "  >>${LOG} 2>&1


##when generated the results:
##md5sum 
#${MD5}  -c tajima/md5/pestPG.md5sum >>${LOG} 2>&1 || exit 1
