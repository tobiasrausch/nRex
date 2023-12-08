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
bcftools mpileup -a AD -B -Q 0 -q 1 -I -f ${GENOME} ${BAM} | bcftools call -P 0.01 -mv -V indels - | bcftools norm -f ${GENOME} -m -both -  | bcftools filter -O b -o ${OUTP}.tmp.bcf -i 'DP>=5 && AD[:1]>=2' -
bcftools index ${OUTP}.tmp.bcf
bcftools isec -O b -o ${OUTP}.common.bcf -n =2 -w 1 ${OUTP}.tmp.bcf common_sites.vcf.gz
bcftools index ${OUTP}.common.bcf
rm ${OUTP}.tmp.bcf ${OUTP}.tmp.bcf.csi

# Fetch input variants
bcftools annotate ${OUTP}.common.bcf -x INFO,^FORMAT/GT | bgzip > ${OUTP}.called.vcf.gz
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
	shapeit4 --input ${OUTP}.chr${CHR}.in.bcf --thread ${THREADS} --map ${BASEDIR}/../maps/chr${CHR}.b38.gmap.gz --region chr${CHR} --reference ${BASEDIR}/../refpanel/chr${CHR}.bcf --output ${OUTP}.chr${CHR}.shapeit.bcf
	bcftools index ${OUTP}.chr${CHR}.shapeit.bcf
	rm ${OUTP}.chr${CHR}.in.bcf ${OUTP}.chr${CHR}.in.bcf.csi
	FILES=${FILES}" ${OUTP}.chr${CHR}.shapeit.bcf"
    fi
done

# Concatenate chromosomes
bcftools concat -O b -o ${OUTP}.shapeit.bcf ${FILES}
bcftools index ${OUTP}.shapeit.bcf
rm ${OUTP}.chr*.shapeit.bcf* ${OUTP}.called.vcf.gz ${OUTP}.called.vcf.gz.tbi

# Remove input variants
rm ${OUTP}.common.bcf ${OUTP}.common.bcf.csi
