#!/bin/bash

if [ $# -ne 3 ]
then
    echo ""
    echo "Usage: $0 <hg19.wgs|hg19.wes|hg19.haloplex> <output prefix> <sample1.bam>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Add all required binaries
export PATH=/g/funcgen/bin:${PATH}

# CMD parameters
ATYPE=${1}
OUTP=${2}
BAMFILE=${3}

# Calculate coverage in windows or exonic regions
if [[ ${ATYPE} = *"wgs"* ]]; then
    alfred count_dna -m 20 -s 10000 -t 10000 -o ${OUTP}.cov.gz ${BAMFILE}
else
    alfred count_dna -m 20 -i ${BASEDIR}/../R/${ATYPE}.bed.gz -o ${OUTP}.cov.gz ${BAMFILE}
fi
