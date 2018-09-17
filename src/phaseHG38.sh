#!/bin/bash

if [ $# -ne 2 ]
then
    echo ""
    echo "Usage: $0 <output prefix> <input.hg38.vcf.gz>"
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

# Fetch input variants
bcftools annotate -O b -o ${OP}.input.bcf -x INFO,^FORMAT/GT ${INVCF}
bcftools index ${OP}.input.bcf

# Phase against 1kGP
FILES=""
for CHR in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
do
    echo "Eagle2 phasing ${CHR}"
    if [ `bcftools view ${OP}.input.bcf ${CHR} | grep -m 1 "^#CHROM" -A 1 | wc -l` -eq 2 ]
    then
	eagle --numThreads ${THREADS} --vcfRef ${BASEDIR}/../refpanelHG38/${CHR}.bcf --vcfTarget ${OP}.input.bcf --geneticMapFile ${BASEDIR}/../refpanelHG38/genetic_map_hg38_withX.txt.gz --outPrefix ${OP}.${CHR}.eagle2 --vcfOutFormat b --chrom ${CHR} 2>&1 | gzip -c > ${OP}.${CHR}.eagle2.log.gz
	bcftools index ${OP}.${CHR}.eagle2.bcf
	FILES=${FILES}" "${OP}.${CHR}.eagle2.bcf
    fi
done
rm ${OP}.input.bcf ${OP}.input.bcf.csi

# Concatenate chromosomes
bcftools concat -O b -o ${OP}.phased.bcf ${FILES}
bcftools index ${OP}.phased.bcf
rm ${OP}.chr*.eagle2.bcf ${OP}.chr*.eagle2.bcf.csi
