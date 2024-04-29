#!/bin/bash

if [ $# -ne 3 ]
then
    echo ""
    echo "Usage: $0 <hg38.wgs|hg38.wes> <genome.fa> <output prefix>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

# CMD parameters
ATYPE=${1}
GENOME=${2}
OUTP=${3}
THREADS=8

# Check BAM
if [ -f ${OUTP}.bam ]
then
    # Run stats using unfiltered BAM
    samtools idxstats -@ ${THREADS} ${OUTP}.bam > ${OUTP}.idxstats
    samtools flagstat -@ ${THREADS} ${OUTP}.bam > ${OUTP}.flagstat
    
    # Run alfred for BAM statistics
    alfred qc -i -b ${BASEDIR}/../genome/${ATYPE}.bed.gz -r ${GENOME} -j ${OUTP}.alfred.json.gz -o ${OUTP}.alfred.tsv.gz ${OUTP}.bam

    # NanoPlot
    NanoPlot -t ${THREADS} --bam ${OUTP}.bam -o ${OUTP}
fi
