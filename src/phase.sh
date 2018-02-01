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

# Fetch input variants
bcftools annotate ${INVCF} -x INFO,^FORMAT/GT | bgzip > ${OP}.called.vcf.gz
tabix ${OP}.called.vcf.gz

# Rename chromosomes
rm -f ${OP}.rename.fwd.chrs ${OP}.rename.rev.chrs
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    echo chr${CHR} ${CHR} >> ${OP}.rename.fwd.chrs
    echo ${CHR} chr${CHR} >> ${OP}.rename.rev.chrs
done
bcftools annotate -O b -o ${OP}.input.bcf --rename-chrs ${OP}.rename.fwd.chrs ${OP}.called.vcf.gz
bcftools index ${OP}.input.bcf
rm ${OP}.called.vcf.gz ${OP}.called.vcf.gz.tbi

# Phase against 1kGP
FILES=""
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    echo "Eagle2 phasing chr${CHR}"
    if [ `bcftools view ${OP}.input.bcf ${CHR} | grep -m 1 "^#CHROM" -A 1 | wc -l` -eq 2 ]
    then
	eagle --numThreads ${THREADS} --vcfRef ${BASEDIR}/../refpanel/chr${CHR}.bcf --vcfTarget ${OP}.input.bcf --geneticMapFile ${BASEDIR}/../refpanel/genetic_map_hg19_withX.txt.gz --outPrefix ${OP}.chr${CHR}.eagle2 --vcfOutFormat b --chrom ${CHR} 2>&1 | gzip -c > ${OP}.chr${CHR}.eagle2.log.gz
	bcftools index ${OP}.chr${CHR}.eagle2.bcf
	FILES=${FILES}" "${OP}.chr${CHR}.eagle2.bcf
    fi
done
rm ${OP}.input.bcf ${OP}.input.bcf.csi

# Concatenate chromosomes
bcftools concat -O b -o ${OP}.eagle2join.bcf ${FILES}
bcftools index ${OP}.eagle2join.bcf

# Rename chromosomes
bcftools annotate -O b -o ${OP}.phased.bcf --rename-chrs ${OP}.rename.rev.chrs ${OP}.eagle2join.bcf
bcftools index ${OP}.phased.bcf

# Clean-up
rm ${OP}.eagle2join.bcf ${OP}.eagle2join.bcf.csi
rm ${OP}.rename.fwd.chrs ${OP}.rename.rev.chrs
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    if [ -f ${OP}.chr${CHR}.eagle2.bcf ]
    then
	rm ${OP}.chr${CHR}.eagle2.bcf ${OP}.chr${CHR}.eagle2.bcf.csi
    fi
done

