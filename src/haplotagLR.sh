#!/bin/bash

if [ $# -ne 4 ]
then
    echo ""
    echo "Usage: $0 <output prefix> <genome.fa> <phased.vcf> <sample.bam>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

# CMD params
THREADS=8
OUTP=${1}
GENOME=${2}
VCF=${3}
BAM=${4}

# Haplotag reads
whatshap haplotag -o ${OUTP}.haplotagged.bam --reference ${GENOME} --ignore-read-groups --tag-supplementary --skip-missing-contigs ${VCF} ${BAM}
samtools index ${OUTP}.haplotagged.bam
