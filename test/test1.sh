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
echo "rm -rf ${ODIR} ${LOG}"
rm -rf ${ODIR} ${LOG}
mkdir -p ${ODIR}

gunzip -c small.bed.gz >small.bed

echo "1) part1" >>${LOG}

for c in 1 2 ##with gc
do
    ${PRG} -h small.vcf.gz -T GT -p 1 -c ${c} -r 100 -O ${ODIR}/vcf1.c${c}.res 2>>${LOG}
    ${PRG} -h small.bcf    -T GT -p 1 -c ${c} -r 100 -O ${ODIR}/vcf2.c${c}.res 2>>${LOG} 
    ${PRG} -h small.vcf.gz       -p 1 -c ${c} -r 100 -O ${ODIR}/vcf3.c${c}.res 2>>${LOG} 
    ${PRG} -h small.bcf          -p 1 -c ${c} -r 100 -O ${ODIR}/vcf4.c${c}.res 2>>${LOG} 
    ${PRG} -P small.bed          -p 1 -c ${c} -r 100 -O ${ODIR}/vcf5.c${c}.res 2>>${LOG} 
    ${PRG} -h small.bcf          -p 1 -c ${c} -r 100 -o 1 -O ${ODIR}/vcf8.c${c}.res 2>>${LOG} 
    ${PRG} -h small.bcf          -p 1 -c ${c} -r 100 -o 1 -O ${ODIR}/vcf9.c${c}.res -R 18 2>>${LOG} 
done

for c in 0 ##without gc
do
    ${PRG} -h small.vcf.gz       -p 1 -c ${c} -r 100 -O ${ODIR}/vcf6.c${c}.res 2>>${LOG} 
    ${PRG} -h small.bcf          -p 1 -c ${c} -r 100 -O ${ODIR}/vcf7.c${c}.res 2>>${LOG} 
done

${PRG} -h small.bcf  -p 1 -r 100 -c 0 -F 1 -O ${ODIR}/vcf8.res  2>>${LOG} 


echo -e "\t test all done "  >>${LOG} 2>&1

#md5sum odir/* >test1.md5

#${MD5} -w -c test1.md5 >>${LOG} 2>&1
