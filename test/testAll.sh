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
    echo "test1.sh: FAIL"
    cat ./test1.sh.log
    RVAL=1
else
    echo "test1.sh: PASS"
fi

exit ${RVAL}
