#!/bin/bash

if [ $# -ne 5 ]
then
    echo "**********************************************************************"
    echo "nRex: Single-nucleotide variant calling."
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "nRex (Version: 0.0.7)"
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <hg19.wgs|hg19.wes|hg19.haloplex> <genome.fa> <output prefix> <sample1.read1.fq.gz> <sample1.read2.fq.gz>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# CMD parameters
ATYPE=${1}
GENOME=${2}
OUTP=${3}

# Align the FASTQs
${BASEDIR}/align.sh $@

# Call variants
${BASEDIR}/call.sh ${ATYPE} ${GENOME} ${OUTP} ${OUTP}.bam

# Calculate coverage
${BASEDIR}/coverage.sh ${ATYPE} ${OUTP} ${OUTP}.bam

# Phase variants against 1000 Genomes reference panel
if [[ ${ATYPE} = *"hg19"* ]]; then
    ${BASEDIR}/phase.sh ${OUTP} ${OUTP}.${ATYPE}.vcf.gz

    # Further regional masks for hg19 exome data
    if [[ ${ATYPE} = *"wes"* ]]; then
	${BASEDIR}/exome_hg19.sh ${OUTP} ${OUTP}.${ATYPE}.vcf.gz
    fi
fi
