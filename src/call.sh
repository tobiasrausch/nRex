#!/bin/bash

if [ $# -lt 4 ]
then
    echo ""
    echo "Usage: $0 <hg38.wgs|hg38.wes> <genome.fa> <output prefix> <sample1.bam> <sample2.bam> ... <sampleN.bam>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

source activate shortread

# CMD parameters
ATYPE=${1}
GENOME=${2}
OUTP=${3}
shift 3
THREADS=8

#freebayes --targets <(zcat ${BASEDIR}/../genome/${ATYPE}.bed.gz) --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.15 --fasta-reference ${GENOME} --genotype-qualities $@ -v ${OUTP}.vcf
freebayes --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.01 --fasta-reference ${GENOME} --genotype-qualities $@ -v ${OUTP}.vcf
bgzip ${OUTP}.vcf
tabix ${OUTP}.vcf.gz

# Normalize VCF
bcftools norm -O b -o ${OUTP}.norm.bcf -a -f ${GENOME} -m -both ${OUTP}.vcf.gz
bcftools index ${OUTP}.norm.bcf

# Filtering
bcftools filter -O b -o ${OUTP}.norm.filtered.bcf -e 'QUAL<=20 || QUAL/INFO/AO<=2 || SAF<=2 || SAR<=2 || RPR<=2 || RPL<=2' ${OUTP}.norm.bcf
bcftools index ${OUTP}.norm.filtered.bcf
rm ${OUTP}.norm.bcf ${OUTP}.norm.bcf.csi

# Subset to target regions
bcftools view -T <(zcat ${BASEDIR}/../genome/${ATYPE}.bed.gz) ${OUTP}.norm.filtered.bcf | bgzip > ${OUTP}.${ATYPE}.vcf.gz
tabix ${OUTP}.${ATYPE}.vcf.gz
rm ${OUTP}.norm.filtered.bcf ${OUTP}.norm.filtered.bcf.csi

source deactivate
