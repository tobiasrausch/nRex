#!/bin/bash

if [ $# -lt 3 ]
then
    echo "Usage: $0 <genome.fa> <tumor folder1> <normal folder2> ..."
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
GENOME=${1}
shift

# Programs
JAVABIN=java
PICARD=${BASEDIR}/picard/dist/picard.jar
SAM=${BASEDIR}/samtools/samtools
BCF=${BASEDIR}/bcftools/bcftools
BAMADDRG=/g/solexa/home/rausch/scripts/cpp/bamaddrg/bamaddrg
FREEBAYES=/g/solexa/home/rausch/scripts/cpp/freebayes/bin/freebayes

# Options
DSTR=$(date +'%a_%y%m%d_%H%M')
TMPDIR=/tmp
if [ -z "${LSF_TYPE}" ]
then
    TMPDIR=`pwd`/tmp_${DSTR}
else
    mkdir -p /scratch/rausch/tmp_${DSTR}
    TMPDIR=/scratch/rausch/tmp_${DSTR}
fi
export TMP=${TMPDIR}
JAVAOPT="-Xmx32g -Djava.io.tmpdir=${TMP}"
PICARDOPT="MAX_RECORDS_IN_RAM=5000000 TMP_DIR=${TMP} VALIDATION_STRINGENCY=SILENT"

# Create IDs and directories
mkdir -p ${TMP}

# Merge and remove duplicates
SAMPLEARR=( $@ )
BAMFILES=""
BAMLIST=""
BAYESLIST=""
for ((  i = 0 ;  i < ${#SAMPLEARR[@]};  i++  ))
do
    SAMPLEID=`echo ${SAMPLEARR[$i]} | sed 's/\/$//' | sed 's/^.*\///'`
    mkdir -p ${SAMPLEID}

    # Have we processed this sample already?
    if [ ! -f ${SAMPLEID}/${SAMPLEID}.rmdup.bam ]
    then
	if [ `ls -lh ${SAMPLEARR[$i]}/*.bam | wc -l | cut -f 1` -eq 1 ]
	then
	    ln -s ${SAMPLEARR[$i]}/*.bam ${SAMPLEID}/${SAMPLEID}.bam
	    ln -s ${SAMPLEARR[$i]}/*.bam.bai ${SAMPLEID}/${SAMPLEID}.bam.bai
	else
	    ${SAM} merge ${SAMPLEID}/${SAMPLEID}.bam ${SAMPLEARR[$i]}/*.bam
	    ${SAM} index ${SAMPLEID}/${SAMPLEID}.bam
	fi

        # Clean .bam file
	${JAVABIN} ${JAVAOPT} -jar ${PICARD}/CleanSam.jar I=${SAMPLEID}/${SAMPLEID}.bam O=${SAMPLEID}/${SAMPLEID}.clean.bam ${PICARDOPT}
	rm ${SAMPLEID}/${SAMPLEID}.bam ${SAMPLEID}/${SAMPLEID}.bam.bai
    
        # Adjusting .bam file read groups to match GATK requirements
	${JAVABIN} ${JAVAOPT} -jar ${PICARD}/AddOrReplaceReadGroups.jar I=${SAMPLEID}/${SAMPLEID}.clean.bam O=${SAMPLEID}/${SAMPLEID}.rg.bam SO=coordinate LB=EMBL RGPL=illumina PU=lane SM=${SAMPLEID} CN=GeneCore ${PICARDOPT}
	rm ${SAMPLEID}/${SAMPLEID}.clean.bam

        # Reorder resulting .bam file karyotypically
	${JAVABIN} ${JAVAOPT} -jar ${PICARD}/ReorderSam.jar I=${SAMPLEID}/${SAMPLEID}.rg.bam O=${SAMPLEID}/${SAMPLEID}.sort.bam REFERENCE=${GENOME} ${PICARDOPT}
	rm ${SAMPLEID}/${SAMPLEID}.rg.bam

        # Remove duplicates
	${JAVABIN} ${JAVAOPT} -jar ${PICARD}/MarkDuplicates.jar I=${SAMPLEID}/${SAMPLEID}.sort.bam O=${SAMPLEID}/${SAMPLEID}.rmdup.bam  M=${SAMPLEID}/${SAMPLEID}.metric.log PG=null MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 ${PICARDOPT}
	rm ${SAMPLEID}/${SAMPLEID}.sort.bam

        # Index resulting, sorted .bam file
	${SAM} index ${SAMPLEID}/${SAMPLEID}.rmdup.bam
    fi

    # Append BAM file
    BAMFILES=${BAMFILES}" -I "${SAMPLEID}/${SAMPLEID}.rmdup.bam" "
    BAMLIST=${BAMLIST}","${SAMPLEID}/${SAMPLEID}.rmdup.bam
    BAMLIST=`echo ${BAMLIST} | sed 's/^,//'`
    BAYESLIST=${BAYESLIST}" -b "${SAMPLEID}/${SAMPLEID}.rmdup.bam" -s "${SAMPLEID}" -r "${SAMPLEID}
done

# Realign indels
SAMPLEID=`echo ${SAMPLEARR[0]} | sed 's/\/$//' | sed 's/^.*\///'`
if [ ! -f ${SAMPLEID}/${SAMPLEID}.target_intervals.list ]
then
    # Create target intervals using GATK's RealignerTargetCreator
    ${JAVABIN} ${JAVAOPT} -jar ${GATK}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${GENOME} ${BAMFILES} -known ${GATK}/bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz -known ${GATK}/bundle/1000G_phase1.indels.hg19.vcf.gz -o ${SAMPLEID}/${SAMPLEID}.target_intervals.list
    cat ${SAMPLEID}/${SAMPLEID}.target_intervals.list | grep -v "^GL000" | grep -v "^MT" | grep -v "^NC_007605" | grep -v "hs37d5" > ${SAMPLEID}/${SAMPLEID}.target_intervals.list.tmp
    mv ${SAMPLEID}/${SAMPLEID}.target_intervals.list.tmp ${SAMPLEID}/${SAMPLEID}.target_intervals.list

    # Realign indels using GATK's IndelRealigner
    ${JAVABIN} ${JAVAOPT} -jar ${GATK}/GenomeAnalysisTK.jar -T IndelRealigner -R ${GENOME} ${BAMFILES} -targetIntervals ${SAMPLEID}/${SAMPLEID}.target_intervals.list -known ${GATK}/bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz -known ${GATK}/bundle/1000G_phase1.indels.hg19.vcf.gz -nWayOut ".bam"
fi

# Replace realigned files
for ((  i = 0 ;  i < ${#SAMPLEARR[@]};  i++  ))
do
    SAMPLEID=`echo ${SAMPLEARR[$i]} | sed 's/\/$//' | sed 's/^.*\///'`

    # Replace rmdup files with realigned files
    if [ -f ${SAMPLEID}.rmdup.bam ]
    then
	mv ${SAMPLEID}.rmdup.bam ${SAMPLEID}/${SAMPLEID}.rmdup.bam
	mv ${SAMPLEID}.rmdup.bai ${SAMPLEID}/${SAMPLEID}.rmdup.bam.bai
    fi
done

# Second realignment round: Initial mpileup indel list
SAMPLEID=`echo ${SAMPLEARR[0]} | sed 's/\/$//' | sed 's/^.*\///'`
if [ ! -f ${SAMPLEID}/${SAMPLEID}.mpileup.vcf ]
then
    ${SAM} mpileup -BDVS -ugf ${GENOME} `echo ${BAMLIST} | sed 's/,/ /g'` | ${BCF} view -cvg - > ${SAMPLEID}/${SAMPLEID}.mpileup.vcf
fi

# Realign indels a second time
if [ ! -f ${SAMPLEID}/${SAMPLEID}.indel.list ]
then
    # Create indel list
    cat ${SAMPLEID}/${SAMPLEID}.mpileup.vcf | grep "INDEL;" | cut -f 1,2 | awk '{print $1":"($2-25)"-"($2+25);}' > ${SAMPLEID}/${SAMPLEID}.indel.list

    # Realign indels using GATK's IndelRealigner
    ${JAVABIN} ${JAVAOPT} -jar ${GATK}/GenomeAnalysisTK.jar -T IndelRealigner -R ${GENOME} ${BAMFILES} -targetIntervals ${SAMPLEID}/${SAMPLEID}.indel.list -nWayOut ".bam"

    # Clean-up mpileup files
    rm ${SAMPLEID}/${SAMPLEID}.mpileup.vcf
fi

# Coverage
for ((  i = 0 ;  i < ${#SAMPLEARR[@]};  i++  ))
do
    SAMPLEID=`echo ${SAMPLEARR[$i]} | sed 's/\/$//' | sed 's/^.*\///'`

    # Replace rmdup files with realigned files
    if [ -f ${SAMPLEID}.rmdup.bam ]
    then
	mv ${SAMPLEID}.rmdup.bam ${SAMPLEID}/${SAMPLEID}.rmdup.bam
	mv ${SAMPLEID}.rmdup.bai ${SAMPLEID}/${SAMPLEID}.rmdup.bam.bai
    fi

    # Calculate Coverage
    if [ ! -f ${SAMPLEID}/${SAMPLEID}.${WIN}.cov ]
    then
	${BED} coverage -counts -abam ${SAMPLEID}/${SAMPLEID}.rmdup.bam -b ${GENOMEID}.${WIN}.win > ${SAMPLEID}/${SAMPLEID}.${WIN}.cov
    fi
done

# Mpileup
SAMPLEID=`echo ${SAMPLEARR[0]} | sed 's/\/$//' | sed 's/^.*\///'`
if [ ! -f ${SAMPLEID}/${SAMPLEID}.mpileup.vcf ]
then
    ${SAM} mpileup -BDVS -ugf ${GENOME} `echo ${BAMLIST} | sed 's/,/ /g'` | ${BCF} view -cvg - > ${SAMPLEID}/${SAMPLEID}.mpileup.vcf
fi

# Freebayes
if [ ! -f ${SAMPLEID}/${SAMPLEID}.freebayes.vcf ]
then
	${BAMADDRG} -c ${BAYESLIST} | ${FREEBAYES} --min-alternate-fraction 0.15 --fasta-reference ${GENOME} --genotype-qualities --stdin -v ${SAMPLEID}/${SAMPLEID}.freebayes.vcf
fi

# Variant calling using GATK's HaplotypeCaller
if [ ! -f ${SAMPLEID}/${SAMPLEID}.gatk.vcf ]
then
    for INTERVAL in `cut -f 1,2 ${GENOME}.fai | egrep -v "GL|NC|hs37d5|MT" | sed 's/\t/:1-/'`
    do
        # Wait for resources
	while [ `ps -a | grep "java" | wc -l | cut -f 1` -gt 12 ]
	do
	    sleep 20
	done
	INTID=`echo ${INTERVAL} | sed 's/:/_/' | sed 's/-/_/'`
	echo ${INTERVAL} ${INTID}
	${JAVABIN} ${JAVAOPT} -jar ${GATK}/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${GENOME} -L ${INTERVAL} ${BAMFILES} -o ${SAMPLEID}/${SAMPLEID}.interval.${INTID}.gatk.vcf > ${SAMPLEID}/${SAMPLEID}.interval.${INTID}.gatk.log 2> ${SAMPLEID}/${SAMPLEID}.interval.${INTID}.gatk.err &
	sleep 20
    done

    # Merge results
    while [ `tail -n 1 ${SAMPLEID}/*.gatk.log | grep "AWS" | wc -l | cut -f 1` -ne 24 ]
    do
	sleep 60
    done
    perl -I${VCFTOOLS} ${VCFTOOLS}/vcf-concat ${SAMPLEID}/${SAMPLEID}.interval.*.gatk.vcf > ${SAMPLEID}/${SAMPLEID}.gatk.vcf
    perl -I${VCFTOOLS} ${VCFTOOLS}/vcf-sort ${SAMPLEID}/${SAMPLEID}.gatk.vcf > ${SAMPLEID}/${SAMPLEID}.gatk.sort.vcf
    mv ${SAMPLEID}/${SAMPLEID}.gatk.sort.vcf ${SAMPLEID}/${SAMPLEID}.gatk.vcf
fi

