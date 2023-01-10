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
export PATH=${BASEDIR}/../conda/bin:${PATH}
source activate align

# CMD params
THREADS=8
OUTP=${1}
INVCF=${2}

# Fetch input variants
bcftools annotate ${INVCF} -x INFO,^FORMAT/GT | bgzip > ${OUTP}.called.vcf.gz
tabix ${OUTP}.called.vcf.gz

# Phase against haplotype reference panel
FILES=""
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
    ls ${BASEDIR}/../refpanel/chr${CHR}.bcf
    if [ -f ${BASEDIR}/../refpanel/chr${CHR}.bcf ]
    then
	echo "Prepare chr${CHR}"
	bcftools view -O b -o ${OUTP}.chr${CHR}.in.bcf ${OUTP}.called.vcf.gz chr${CHR}
	bcftools index ${OUTP}.chr${CHR}.in.bcf

	echo "Shapeit4 phasing chr${CHR}"
	${BASEDIR}/../shapeit4-4.2.2/bin/shapeit4.2 --input ${OUTP}.chr${CHR}.in.bcf --thread ${THREADS} --map ${BASEDIR}/../shapeit4-4.2.2/maps/chr${CHR}.b38.gmap.gz --region chr${CHR} --reference ${BASEDIR}/../refpanel/chr${CHR}.bcf --output ${OUTP}.chr${CHR}.shapeit.bcf
	bcftools index ${OUTP}.chr${CHR}.shapeit.bcf
	rm ${OUTP}.chr${CHR}.in.bcf ${OUTP}.chr${CHR}.in.bcf.csi
	FILES=${FILES}" ${OUTP}.chr${CHR}.shapeit.bcf"
    fi
done

# Concatenate chromosomes
bcftools concat -O b -o ${OUTP}.shapeit.bcf ${FILES}
bcftools index ${OUTP}.shapeit.bcf
rm ${OUTP}.chr*.shapeit.bcf* ${OUTP}.called.vcf.gz ${OUTP}.called.vcf.gz.tbi
