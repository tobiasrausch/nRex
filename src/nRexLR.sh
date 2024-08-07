#!/bin/bash

if [ $# -ne 2 ]
then
    echo "**********************************************************************"
    echo "nRex: Single-nucleotide variant calling."
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "nRex (Version: 0.2.1)"
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <output prefix> [sample.fq.gz|sample.unaligned.bam]"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# CMD parameters
VERSION=38
ATYPE=hg${VERSION}.wgs    # [hg38.wgs|hg38.wes]
GENOME=${BASEDIR}/../genomeLR/hg${VERSION}.fa
MAP=${BASEDIR}/../genomeLR/Homo_sapiens.GRCh${VERSION}.dna.primary_assembly.fa.r101.s501.blacklist.gz
OUTP=${1}
FQ=${2}

# Alignment
if [ ! -f ${OUTP}.bam ]
then
    ${BASEDIR}/alignLR.sh ${ATYPE} ${GENOME} ${OUTP} ${FQ}
fi

# Alignment statistics
if [ ! -f ${OUTP}.alfred.tsv.gz ]
then
    ${BASEDIR}/statsLR.sh ${ATYPE} ${GENOME} ${OUTP}
fi

# SNP calling
if [ ! -f ${OUTP}.${ATYPE}.vcf.gz ]
then
    ${BASEDIR}/callLR.sh ${ATYPE} ${GENOME} ${OUTP} ${OUTP}.bam
fi

# Calculate coverage
if [ ! -f ${OUTP}.cov.gz ]
then
    ${BASEDIR}/coverageLR.sh ${ATYPE} ${GENOME} ${MAP} ${OUTP} ${OUTP}.bam
fi

# Haplotag BAM (if not done yet)
if [ ! -f ${OUTP}.haplotagged.bam ]
then
    ${BASEDIR}/haplotagLR.sh ${OUTP} ${GENOME} ${OUTP}.bam
fi

# Structural variants [can be jointly run on multiple BAM files]
if [ ! -f ${OUTP}.delly.bcf ]
then
    ${BASEDIR}/dellyLR.sh ${ATYPE} ${GENOME} ${OUTP} ${OUTP}.bam
fi
