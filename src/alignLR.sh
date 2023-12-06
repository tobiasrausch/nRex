#!/bin/bash

if [ $# -ne 4 ]
then
    echo ""
    echo "Usage: $0 <hg38.wgs|hg38.wes> <genome.fa> <output prefix> <sample1.fq.gz>"
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
FQ=${4}
THREADS=8

# BWA alignment
minimap2 -t ${THREADS} -a -x map-ont -L ${GENOME} ${FQ} | samtools sort -o ${OUTP}.bam
samtools index ${OUTP}.bam

# Run stats using unfiltered BAM
samtools idxstats -@ ${THREADS} ${OUTP}.bam > ${OUTP}.idxstats
samtools flagstat -@ ${THREADS} ${OUTP}.bam > ${OUTP}.flagstat

# Run alfred for BAM statistics
alfred qc -b ${BASEDIR}/../genome/${ATYPE}.bed.gz -r ${GENOME} -j ${OUTP}.alfred.json.gz -o ${OUTP}.alfred.tsv.gz ${OUTP}.bam

# NanoPlot
NanoPlot -t ${THREADS} --bam ${OUTP}.bam -o ${OUTP}
