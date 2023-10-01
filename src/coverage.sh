#!/bin/bash

if [ $# -ne 5 ]
then
    echo ""
    echo "Usage: $0 <hg38.wgs|hg38.wes|hg38.haloplex> <genome.fa> <map.fa.gz> <output prefix> <sample1.bam>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

source activate align

# CMD parameters
ATYPE=${1}
GENOME=${2}
MAP=${3}
OUTP=${4}
BAMFILE=${5}

# Calculate coverage in windows or exonic regions
if [[ ${ATYPE} = *"wgs"* ]]; then
    WIN=10000
    delly cnv -i ${WIN} -j ${WIN} -w ${WIN} -g ${GENOME} -m ${MAP} -o ${OUTP}.cov.bcf -c ${OUTP}.cov.gz ${BAMFILE}
else
    delly cnv -r ${BASEDIR}/../genome/${ATYPE}.bed.gz -n -b ${BASEDIR}/../genome/${ATYPE}.bed.gz -g ${GENOME} -m ${MAP} -o ${OUTP}.cov.bcf -c ${OUTP}.cov.gz ${BAMFILE}
fi
