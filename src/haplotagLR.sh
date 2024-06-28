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

# Generate phased blocks
whatshap phase --ignore-read-groups --reference ${GENOME} ${VCF} ${BAM} -o ${OUTP}.whatshap.vcf
bgzip ${OUTP}.whatshap.vcf
tabix ${OUTP}.whatshap.vcf.gz

# Haplotag reads
whatshap haplotag -o ${OUTP}.haplotagged.bam --reference ${GENOME} --ignore-read-groups --tag-supplementary --skip-missing-contigs ${OUTP}.whatshap.vcf.gz ${BAM}
samtools index ${OUTP}.haplotagged.bam
