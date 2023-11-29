#!/bin/bash

if [ $# -ne 3 ]
then
    echo ""
    echo "Usage: $0 <output prefix> <sample1.read1.fq.gz> <sample1.read2.fq.gz>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

# Input parameters
OUTP=${1}
FQ1=${2}
FQ2=${3}
THREADS=8

# FastQC
mkdir fastqc_${OUTP}
fastqc -t ${THREADS} -o fastqc_${OUTP}/ ${FQ1}
fastqc -t ${THREADS} -o fastqc_${OUTP}/ ${FQ2}

# Adapter trimming
trim_galore --paired --basename ${OUTP} ${FQ1} ${FQ2}
