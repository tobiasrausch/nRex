#!/bin/bash

if [ $# -ne 3 ]
then
    echo "**********************************************************************"
    echo "nRex: Single-nucleotide variant calling."
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "nRex (Version: 0.2.1)"
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
ATYPE=hg${VERSION}.wgs    # [hg38.wgs|hg38.wes]
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

# Read-depth
if [ ! -f ${OUTP}.cov.gz ]
then
    ${BASEDIR}/coverage.sh ${ATYPE} ${GENOME} ${MAP} ${OUTP} ${OUTP}.bam
fi

# Call variants [can be jointly run on multiple BAM files]
if [ ! -f ${OUTP}.vcf.gz ]
then
    ${BASEDIR}/call.sh ${ATYPE} ${GENOME} ${OUTP} ${OUTP}.bam
fi

# Phase variants
if [ ! -f ${OUTP}.shapeit.bcf ]
then
    if [[ ${ATYPE} = *"hg38"* ]]; then
	if [ -f ${BASEDIR}/../maps/chr21.b38.gmap.gz ]; then
	    if [ -f ${BASEDIR}/../refpanel/chr21.bcf ]; then
		${BASEDIR}/phase.sh ${OUTP} ${OUTP}.${ATYPE}.vcf.gz
	    fi
	fi
    fi
fi

# Structural variants [can be jointly run on multiple BAM files]
if [ ! -f ${OUTP}.delly.bcf ]
then
    ${BASEDIR}/delly.sh ${ATYPE} ${GENOME} ${OUTP} ${OUTP}.bam
fi

# QC summary
${BASEDIR}/qc.sh ${OUTP}

# Optional: Annotate variants
if [ ! -f ${OUTP}.vep.bcf ]
then
    if [ -d ${BASEDIR}/../vep_data ]
    then
	${BASEDIR}/vep.sh ${OUTP} ${OUTP}.${ATYPE}.vcf.gz
    fi
fi

