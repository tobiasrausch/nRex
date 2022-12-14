#!/bin/bash

if [ $# -ne 3 ]
then
    echo "**********************************************************************"
    echo "nRex: Single-nucleotide variant calling."
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "nRex (Version: 0.1.8)"
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <output prefix> <sample1.read1.fq.gz> <sample1.read2.fq.gz>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# CMD parameters
VERSION=38
ATYPE=hg${VERSION}.wgs    # [hg38.wgs|hg38.wes|hg38.haloplex]
GENOME=${BASEDIR}/../genome/hg${VERSION}.fa
MAP=${BASEDIR}/../genome/Homo_sapiens.GRCh${VERSION}.dna.primary_assembly.fa.r101.s501.blacklist.gz
OUTP=${1}
FQ1=${2}
FQ2=${3}

if [ ! -f ${OUTP}.bam ]
then
    # Trim the FASTQs
    ${BASEDIR}/trim.sh ${OUTP} ${FQ1} ${FQ2}

    # Align the adapter-filtered FASTQs
    ${BASEDIR}/align.sh ${ATYPE} ${GENOME} ${OUTP} ${OUTP}_val_1.fq.gz ${OUTP}_val_2.fq.gz
    rm ${OUTP}_val_1.fq.gz ${OUTP}_val_2.fq.gz
fi

# Call variants
${BASEDIR}/call.sh ${ATYPE} ${GENOME} ${OUTP} ${OUTP}.bam

# Calculate coverage
${BASEDIR}/coverage.sh ${ATYPE} ${GENOME} ${MAP} ${OUTP} ${OUTP}.bam

# QC summary
python ${BASEDIR}/../scripts/qc.py -p ${OUTP} > ${OUTP}.qc.summary

# Phase variants against 1000 Genomes reference panel
if [[ ${ATYPE} = *"hg19"* ]]; then
    ${BASEDIR}/phase.sh ${OUTP} ${OUTP}.${ATYPE}.vcf.gz

    # Further regional masks for hg19 exome data
    if [[ ${ATYPE} = *"wes"* ]]; then
	${BASEDIR}/exome_hg19.sh ${OUTP} ${OUTP}.${ATYPE}.vcf.gz
    fi
fi
