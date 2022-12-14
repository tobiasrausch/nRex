#!/bin/bash

if [ $# -lt 4 ]
then
    echo ""
    echo "Usage: $0 <hg19.wgs|hg19.wes|hg19.haloplex> <genome.fa> <output prefix> <sample1.bam> <sample2.bam> ... <sampleN.bam>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}
source activate align

# CMD parameters
ATYPE=${1}
GENOME=${2}
OUTP=${3}
shift 3
THREADS=8

freebayes --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.15 --fasta-reference ${GENOME} --genotype-qualities $@ -v ${OUTP}.vcf
#freebayes --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.01 --fasta-reference ${GENOME} --genotype-qualities $@ -v ${OUTP}.vcf
bgzip ${OUTP}.vcf
tabix ${OUTP}.vcf.gz

# Normalize VCF
bcftools norm -O b -o ${OUTP}.norm.bcf -f ${GENOME} -m -both ${OUTP}.vcf.gz
bcftools index ${OUTP}.norm.bcf
rm ${OUTP}.vcf.gz ${OUTP}.vcf.gz.tbi

# Filtering
if [[ ${ATYPE} = *"haloplex"* ]]; then
    bcftools filter -O b -o ${OUTP}.norm.filtered.bcf -e 'QUAL<=20 || QUAL/INFO/AO<=2 || SAF<=2 || SAR<=2' ${OUTP}.norm.bcf
else
    bcftools filter -O b -o ${OUTP}.norm.filtered.bcf -e 'QUAL<=20 || QUAL/INFO/AO<=2 || SAF<=2 || SAR<=2 || RPR<=2 || RPL<=2' ${OUTP}.norm.bcf
fi
bcftools index ${OUTP}.norm.filtered.bcf
rm ${OUTP}.norm.bcf ${OUTP}.norm.bcf.csi

# Subset to target regions
bcftools view -T <(zcat ${BASEDIR}/../genome/${ATYPE}.bed.gz) ${OUTP}.norm.filtered.bcf | bgzip > ${OUTP}.${ATYPE}.vcf.gz
tabix ${OUTP}.${ATYPE}.vcf.gz
rm ${OUTP}.norm.filtered.bcf ${OUTP}.norm.filtered.bcf.csi
