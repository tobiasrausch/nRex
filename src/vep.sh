#! /bin/bash

if [ $# -ne 2 ]
then
    echo ""
    echo "Usage: $0 <outprefix> <sample.vcf.gz>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}
source activate align

OUTP=${1}
INVCF=${2}

# Create input and output directory for VEP
mkdir -p ${BASEDIR}/../vep_data/user
chmod a+rwx ${BASEDIR}/../vep_data/user

# Run VEP
bcftools view ${INVCF} | bgzip > ${BASEDIR}/../vep_data/user/${OUTP}.in.vcf.gz
docker run -v ${BASEDIR}/../vep_data:/opt/vep/.vep ensemblorg/ensembl-vep vep --assembly GRCh38 --everything --no_stats --cache --offline --compress_output bgzip --format vcf --vcf --force_overwrite --dir_cache /opt/vep/.vep/ --dir_plugins /opt/vep/.vep/Plugins/ --input_file /opt/vep/.vep/user/${OUTP}.in.vcf.gz --output_file /opt/vep/.vep/user/${OUTP}.out.vcf.gz
#--plugin dbNSFP,/opt/vep/.vep/anno/dbNSFP4.2a_grch37.gz,ALL

# Convert to BCF
bcftools view -O b -o ${OUTP}.vep.bcf ${BASEDIR}/../vep_data/user/${OUTP}.out.vcf.gz
bcftools index ${OUTP}.vep.bcf
rm /opt/vep/.vep/user/${OUTP}.in.vcf.gz /opt/vep/.vep/user/${OUTP}.out.vcf.gz
