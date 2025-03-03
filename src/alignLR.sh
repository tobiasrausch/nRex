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

source activate longread

# CMD parameters
ATYPE=${1}
GENOME=${2}
OUTP=${3}
FQ=${4}
THREADS=8

# BAM input?
samtools quickcheck ${FQ}
if [ $? -eq 0 ]
then
    # BAM input
    samtools fastq -@ ${THREADS} ${FQ} | gzip -c > ${OUTP}.fq.gz
    FQ=${OUTP}.fq.gz
    
    # Minimap2 alignment
    minimap2 -t ${THREADS} -a -x map-ont -L ${GENOME} ${FQ} | samtools sort -o ${OUTP}.bam
    samtools index ${OUTP}.bam

    # Clean-up
    rm ${OUTP}.fq.gz
else
    # Minimap2 alignment
    minimap2 -t ${THREADS} -a -x map-ont -L ${GENOME} ${FQ} | samtools sort -o ${OUTP}.bam
    samtools index ${OUTP}.bam
fi
