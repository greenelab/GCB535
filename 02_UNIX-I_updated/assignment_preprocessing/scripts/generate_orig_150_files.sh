#!/bin/bash

## alex amlie-wolf
## generate a tab separated file of original 150 pokemon

OUTFILE=$1
DATADIR=$2

printf "name\ttype1\ttype2\tnumber\n" > ${OUTFILE}

for NUM in {1..151}; do
    NAME=`grep "Name" ${DATADIR}/${NUM}.dat | cut -d'|' -f2`
    TYPE1=`grep "Type1" ${DATADIR}/${NUM}.dat | cut -d'|' -f2`
    TYPE2=`grep "Type2" ${DATADIR}/${NUM}.dat | cut -d'|' -f2`
    NUMBER=`grep "Number" ${DATADIR}/${NUM}.dat | cut -d'|' -f2`
    if [ -z "${TYPE2%$'\r'}" ]; then
	TYPE2="None"
    fi

    printf "%s\t%s\t%s\t%s\n" "${NAME%$'\r'}" "${TYPE1%$'\r'}" "${TYPE2%$'\r'}" "${NUMBER%$'\r'}" >> ${OUTFILE}
done
