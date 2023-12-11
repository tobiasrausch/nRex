#!/bin/bash

if [ $# -ne 4 ]
then
    echo ""
    echo "Usage: $0 <hg38.wgs|hg38.wes> <genome.fa> <output prefix> <sample.bam"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

# CMD parameters
ATYPE=${1}
GENOME=${2}
OUTP=${3}
BAM=${4}
THREADS=8

# Pull DeepVariant image
BIN_VERSION=1.6.0
INPUT_DIR="."
OUTPUT_DIR="."
singularity pull docker://google/deepvariant:"${BIN_VERSION}"

# Run the container
cp ${GENOME} ${OUTP}.genome.fa
cp ${GENOME}.fai ${OUTP}.genome.fa.fai
singularity run --bind $(pwd) -B /usr/lib/locale/:/usr/lib/locale/ docker://google/deepvariant:"${BIN_VERSION}" /opt/deepvariant/bin/run_deepvariant --model_type=ONT_R104 --ref=${OUTP}.genome.fa --reads=${BAM} --output_vcf=${OUTP}.vcf.gz --output_gvcf=${OUTP}.g.vcf.gz --num_shards=${THREADS}
rm ${OUTP}.genome.fa ${OUTP}.genome.fa.fai

# Normalize VCF
bcftools norm -O b -o ${OUTP}.norm.bcf -a -f ${GENOME} -m -both ${OUTP}.vcf.gz
bcftools index ${OUTP}.norm.bcf

# Filtering
bcftools filter -O b -o ${OUTP}.norm.filtered.bcf -e 'QUAL<=10' ${OUTP}.norm.bcf
bcftools index ${OUTP}.norm.filtered.bcf
rm ${OUTP}.norm.bcf ${OUTP}.norm.bcf.csi

# Subset to target regions
bcftools view -T <(zcat ${BASEDIR}/../genomeLR/${ATYPE}.bed.gz) ${OUTP}.norm.filtered.bcf | bgzip > ${OUTP}.${ATYPE}.vcf.gz
tabix ${OUTP}.${ATYPE}.vcf.gz
rm ${OUTP}.norm.filtered.bcf ${OUTP}.norm.filtered.bcf.csi
