#!/bin/bash

if [ $# -ne 5 ]
then
    echo ""
    echo "Usage: $0 <hg38.wgs|hg38.wes> <genome.fa> <map.fa.gz> <output prefix> <sample1.bam>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

# CMD parameters
ATYPE=${1}
GENOME=${2}
MAP=${3}
OUTP=${4}
BAMFILE=${5}
WIN=25000

# Calculate coverage
delly cnv -i ${WIN} -j ${WIN} -w ${WIN} -g ${GENOME} -m ${MAP} -o ${OUTP}.cov.bcf -c ${OUTP}.cov.gz ${BAMFILE}
