#!/bin/bash

if [ $# -ne 2 ]
then
    echo "**********************************************************************"
    echo "nRex: Single-nucleotide variant calling."
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "nRex (Version: 0.1.8)"
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <output prefix> <sample.fq.gz>"
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

if [ ! -f ${OUTP}.bam ]
then
    ${BASEDIR}/alignLR.sh ${ATYPE} ${GENOME} ${OUTP} ${FQ}
fi

# Calculate coverage
if [ ! -f ${OUTP}.cov.gz ]
then
    ${BASEDIR}/coverage.sh ${ATYPE} ${GENOME} ${MAP} ${OUTP} ${OUTP}.bam
fi

# Structural variants [can be jointly run on multiple BAM files]
if [ ! -f ${OUTP}.delly.bcf ]
then
    ${BASEDIR}/dellyLR.sh ${ATYPE} ${GENOME} ${OUTP} ${OUTP}.bam
fi
