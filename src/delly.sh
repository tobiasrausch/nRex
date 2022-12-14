#!/bin/bash

if [ $# -lt 4 ]
then
    echo ""
    echo "Usage: $0 <hg38.wgs|hg38.wes|hg38.haloplex> <genome.fa> <output prefix> <sample1.bam> <sample2.bam> ... <sampleN.bam>"
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
OUTP=${3}
EXCL=${BASEDIR}/../genome/human.hg38.excl.tsv
shift 3

# Delly
delly call -g ${GENOME} -x ${EXCL} -o ${OUTP}.delly.bcf $@
