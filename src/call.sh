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
export PATH=${BASEDIR}/../mamba/bin:${PATH}
source activate align

# CMD parameters
ATYPE=${1}
GENOME=${2}
OUTP=${3}
shift 3
THREADS=8

freebayes --targets <(zcat ${BASEDIR}/../genome/${ATYPE}.bed.gz) --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.15 --fasta-reference ${GENOME} --genotype-qualities $@ -v ${OUTP}.vcf
#freebayes --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.01 --fasta-reference ${GENOME} --genotype-qualities $@ -v ${OUTP}.vcf
bgzip ${OUTP}.vcf
tabix ${OUTP}.vcf.gz
