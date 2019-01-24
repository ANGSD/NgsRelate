#!/bin/bash

PRG=""
if [ $# -eq 0 ] 
then
    exit 1;
fi

if [ $# -eq 1 ]
then
    PRG=$1
fi


echo "--------------------"
echo "Using PRG: '${PRG}'"
echo "--------------------"

WDIR=`dirname $PRG`

RVAL=0

echo "Testing bcf/vcf and plink"
./test1.sh ${PRG}
if [ ! $? -eq 0 ] ;then
    cat ./test1.sh.log
    RVAL=1
fi

exit ${RVAL}
