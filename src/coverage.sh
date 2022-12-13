#!/bin/bash

if [ $# -ne 5 ]
then
    echo ""
    echo "Usage: $0 <hg19.wgs|hg19.wes|hg19.haloplex> <genome.fa> <map.fa.gz> <output prefix> <sample1.bam>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate align

# CMD parameters
ATYPE=${1}
GENOME=${2}
MAP=${3}
OUTP=${4}
BAMFILE=${5}

# Calculate coverage in windows or exonic regions
if [[ ${ATYPE} = *"wgs"* ]]; then
    WIN=1000
    delly cnv -i ${WIN} -j ${WIN} -w ${WIN} -g ${GENOME} -m ${MAP} -o ${OUTP}.cov.bcf -c ${OUTP}.cov.gz ${BAMFILE}
else
    delly cnv -r <(zcat ${BASEDIR}/../genome/${ATYPE}.bed.gz) -n -b <(zcat ${BASEDIR}/../genome/${ATYPE}.bed.gz) -g ${GENOME} -m ${MAP} -o ${OUTP}.cov.bcf -c ${OUTP}.cov.gz ${BAMFILE}
fi
