#!/bin/bash

if [ $# -ne 3 ]
then
    echo ""
    echo "Usage: $0 <output prefix> <genome.fa> <sample.bam>"
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
BAM=${3}

# Genotype known variants
longshot --no_haps --min_cov 3 --potential_variants ${BASEDIR}/../genomeLR/common_sites.vcf.gz --min_alt_count 2 --min_alt_frac 0.01 -s ${OUTP} --bam ${BAM} --ref ${GENOME} --out ${OUTP}.longshot.vcf
bgzip ${OUTP}.longshot.vcf
tabix ${OUTP}.longshot.vcf.gz

# Generate phased blocks
whatshap phase --ignore-read-groups --reference ${GENOME} ${OUTP}.longshot.vcf.gz ${BAM} -o ${OUTP}.whatshap.vcf
bgzip ${OUTP}.whatshap.vcf
tabix ${OUTP}.whatshap.vcf.gz

# Haplotag reads
whatshap haplotag -o ${OUTP}.haplotagged.bam --reference ${GENOME} --ignore-read-groups --tag-supplementary --skip-missing-contigs ${OUTP}.whatshap.vcf.gz ${BAM}
samtools index ${OUTP}.haplotagged.bam
