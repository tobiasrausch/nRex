#!/bin/bash

if [ $# -ne 5 ]
then
    echo "**********************************************************************"
    echo "nRex: Single-nucleotide variant calling."
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "nRex (Version: 0.0.6)"
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <wgs|wes|haloplex> <genome.fa> <output prefix> <sample1.read1.fq.gz> <sample1.read2.fq.gz>"
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

# Phase variants against 1000 Genomes reference panel
${BASEDIR}/phase.sh ${OUTP} ${OUTP}.${ATYPE}.vcf.gz
