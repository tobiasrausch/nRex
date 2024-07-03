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

# Annotate common SNPs
bcftools mpileup --ignore-RG -a AD -B -Q 0 -q 1 -I -f ${GENOME} ${BAM} | bcftools call -P 0.01 -mv -V indels - | bcftools norm -f ${GENOME} -m -both -  | bcftools filter -O b -o ${OUTP}.tmp.bcf -i 'DP>=5 && AD[:1]>=2' -
bcftools index ${OUTP}.tmp.bcf
bcftools isec -O b -o ${OUTP}.common.bcf -n =2 -w 1 ${OUTP}.tmp.bcf ${BASEDIR}/../genomeLR/common_sites.vcf.gz
bcftools index ${OUTP}.common.bcf
rm ${OUTP}.tmp.bcf ${OUTP}.tmp.bcf.csi

# Fetch input variants
bcftools annotate ${OUTP}.common.bcf -x INFO,^FORMAT/GT | bgzip > ${OUTP}.called.vcf.gz
tabix ${OUTP}.called.vcf.gz
rm ${OUTP}.common.bcf ${OUTP}.common.bcf.csi

# Generate phased blocks
whatshap phase --ignore-read-groups --reference ${GENOME} ${OUTP}.called.vcf.gz ${BAM} -o ${OUTP}.whatshap.vcf
bgzip ${OUTP}.whatshap.vcf
tabix ${OUTP}.whatshap.vcf.gz
rm ${OUTP}.called.vcf.gz ${OUTP}.called.vcf.gz.tbi

# Haplotag reads
whatshap haplotag -o ${OUTP}.haplotagged.bam --reference ${GENOME} --ignore-read-groups --tag-supplementary --skip-missing-contigs ${OUTP}.whatshap.vcf.gz ${BAM}
samtools index ${OUTP}.haplotagged.bam
