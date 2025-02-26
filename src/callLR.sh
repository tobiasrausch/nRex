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

source activate longread

# CMD parameters
ATYPE=${1}
GENOME=${2}
OUTP=${3}
shift 3
THREADS=8

source activate clair3
## Without phasing
run_clair3.sh --bam_fn=${@} --ref_fn=${GENOME} --threads=${THREADS} --platform="ont" --model_path=${BASEDIR}/../models/r1041_e82_260bps_sup_v400 --output=${OUTP}_clair3
## With phasing
#run_clair3.sh --bam_fn=${@} --ref_fn=${GENOME} --threads=${THREADS} --platform="ont" --model_path=${BASEDIR}/../models/r1041_e82_260bps_sup_v400 --output=${OUTP}_clair3 --use_whatshap_for_final_output_haplotagging --use_whatshap_for_final_output_phasing --enable_phasing
if [ -f ${OUTP}_clair3/phased_merge_output.vcf.gz ]
then
    cp ${OUTP}_clair3/phased_merge_output.vcf.gz ${OUTP}.vcf.gz
    cp ${OUTP}_clair3/phased_merge_output.vcf.gz.tbi ${OUTP}.vcf.gz.tbi
    cp ${OUTP}_clair3/phased_output.bam ${OUTP}.haplotagged.bam
    samtools index ${OUTP}.haplotagged.bam
else
    cp ${OUTP}_clair3/merge_output.vcf.gz ${OUTP}.vcf.gz
    cp ${OUTP}_clair3/merge_output.vcf.gz.tbi ${OUTP}.vcf.gz.tbi
fi
rm -rf ${OUTP}_clair3/
source deactivate

# Normalize VCF
bcftools norm -O b -o ${OUTP}.norm.bcf -a -f ${GENOME} -m -both ${OUTP}.vcf.gz
bcftools index ${OUTP}.norm.bcf

# Filtering
bcftools filter -O b -o ${OUTP}.norm.filtered.bcf -e 'QUAL<=5' ${OUTP}.norm.bcf
bcftools index ${OUTP}.norm.filtered.bcf
rm ${OUTP}.norm.bcf ${OUTP}.norm.bcf.csi

# Subset to target regions
bcftools view -T <(zcat ${BASEDIR}/../genome/${ATYPE}.bed.gz) ${OUTP}.norm.filtered.bcf | bgzip > ${OUTP}.${ATYPE}.vcf.gz
tabix ${OUTP}.${ATYPE}.vcf.gz
rm ${OUTP}.norm.filtered.bcf ${OUTP}.norm.filtered.bcf.csi
