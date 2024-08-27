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

## FASTQ pre-processing
fastp --thread ${THREADS} -i ${FQ1} -I ${FQ2} -o ${OUTP}.R1.fq.gz -O ${OUTP}.R2.fq.gz -j ${OUTP}.json -h ${OUTP}.html
