#!/bin/bash

if [ $# -ne 2 ]
then
    echo ""
    echo "Usage: $0 <output prefix> <input.hg19.vcf.gz>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

export PATH=/g/funcgen/bin/:${PATH}

# CMD params
THREADS=4
OP=${1}
INVCF=${2}

# Keep only ExAC regions
bcftools view -O b -o ${OP}.exac.bcf -T <(zcat ${BASEDIR}/../R/exac.hg19.bed.gz) ${INVCF}

# Remove low complexity regions
bcftools view -O b -o ${OP}.exac.sgdp.bcf -T ^<(zcat ${BASEDIR}/../R/sgdp.hg19.bed.gz) ${OP}.exac.bcf
bcftools index ${OP}.exac.sgdp.bcf
rm ${OP}.exac.bcf
