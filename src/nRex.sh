#!/bin/bash

if [ $# -lt 4 ]
then
    echo "Usage: $0 <genome.fa> <output prefix> <sample1.read1.fq.gz> <sample1.read2.fq.gz> ..."
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

export PERL5LIB=${BASEDIR}/perl/lib/perl5/:${BASEDIR}/perl/lib/5.24.0/
export PATH=${BASEDIR}/perl/bin:${PATH}

GENOME=${1}
shift
OUTP=${1}
shift
THREADS=4

# Programs
PICARD=${BASEDIR}/picard/picard.jar
SAM=${BASEDIR}/samtools/samtools
BCF=${BASEDIR}/bcftools/bcftools
BWA=${BASEDIR}/bwa/bwa
BAMSTATS=${BASEDIR}/bamStats/src/bamStats
BAMSTATR=${BASEDIR}/bamStats/R
FREEBAYES=${BASEDIR}/freebayes/bin/freebayes
TABIX=${BASEDIR}/htslib/tabix
BGZIP=${BASEDIR}/htslib/bgzip
VEP=${BASEDIR}/vep/vep.pl
VEP_DATA=${BASEDIR}/vepcache
JAVA=java

# Tmp directory
DSTR=$(date +'%a_%y%m%d_%H%M')
export TMP=/tmp/tmp_nrex_${DSTR}
mkdir -p ${TMP}
JAVAOPT="-Xms4g -Xmx16g -XX:ParallelGCThreads=4 -Djava.io.tmpdir=${TMP}"
PICARDOPT="MAX_RECORDS_IN_RAM=5000000 TMP_DIR=${TMP} VALIDATION_STRINGENCY=SILENT"

# Align
BAMLIST=""
mkdir -p ${OUTP}
SAMPLEARR=( $@ )
for ((  i = 0 ;  i < ${#SAMPLEARR[@]};  i = i + 2  ))
do
    j=`expr ${i} + 1`
    
    # Generate IDs
    BAMID=`echo ${OUTP} | sed "s/$/.align${i}/"`

    # BWA
    ${BWA} mem -t ${THREADS} -R "@RG\tID:${OUTP}\tSM:${OUTP}" ${GENOME} ${SAMPLEARR[$i]} ${SAMPLEARR[$j]} | ${SAM} view -bT ${GENOME} - > ${OUTP}/${BAMID}.bam

    # Sort & Index
    ${SAM} sort -o ${OUTP}/${BAMID}.srt.bam ${OUTP}/${BAMID}.bam && rm ${OUTP}/${BAMID}.bam && ${SAM} index ${OUTP}/${BAMID}.srt.bam

    # Clean .bam file
    ${JAVA} ${JAVAOPT} -jar ${PICARD} CleanSam I=${OUTP}/${BAMID}.srt.bam O=${OUTP}/${BAMID}.srt.clean.bam ${PICARDOPT} && rm ${OUTP}/${BAMID}.srt.bam*

    # Mark duplicates
    ${JAVA} ${JAVAOPT} -jar ${PICARD} MarkDuplicates I=${OUTP}/${BAMID}.srt.clean.bam O=${OUTP}/${BAMID}.srt.clean.rmdup.bam M=${OUTP}/${OUTP}.markdups.log PG=null MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 ${PICARDOPT} && rm ${OUTP}/${BAMID}.srt.clean.bam* && ${SAM} index ${OUTP}/${BAMID}.srt.clean.rmdup.bam
    # No mark duplicates (Haloplex)
    # mv ${OUTP}/${BAMID}.srt.clean.bam ${OUTP}/${BAMID}.srt.clean.rmdup.bam && ${SAM} index ${OUTP}/${BAMID}.srt.clean.rmdup.bam
    
    # Run stats using unfiltered BAM
    ${SAM} idxstats ${OUTP}/${BAMID}.srt.clean.rmdup.bam > ${OUTP}/${OUTP}.idxstats
    ${SAM} flagstat ${OUTP}/${BAMID}.srt.clean.rmdup.bam > ${OUTP}/${OUTP}.flagstat

    # Filter duplicates, unmapped reads, chrM and unplaced contigs
    #CHRS=`cat ${BASEDIR}/../bed/haloplex.bed | cut -f 1 | sort -k1,1V -k2,2n | uniq | tr '\n' ' '`
    CHRS=`cat ${BASEDIR}/../bed/exome.bed | cut -f 1 | sort -k1,1V -k2,2n | uniq | tr '\n' ' '`
    ${SAM} view -F 1024 -b ${OUTP}/${BAMID}.srt.clean.rmdup.bam ${CHRS} > ${OUTP}/${BAMID}.final.bam
    ${SAM} index ${OUTP}/${BAMID}.final.bam
    rm ${OUTP}/${BAMID}.srt.clean.rmdup.bam ${OUTP}/${BAMID}.srt.clean.rmdup.bam.bai

    # Run stats using filtered BAM, TSS enrichment, error rates, etc.
    ${BAMSTATS} -b ${BASEDIR}/../bed/exome.bed -r ${GENOME} -o ${OUTP}/${OUTP}.bamStats ${OUTP}/${BAMID}.final.bam
    #${BAMSTATS} -b ${BASEDIR}/../bed/haloplex.bed -r ${GENOME} -o ${OUTP}/${OUTP}.bamStats ${OUTP}/${BAMID}.final.bam

    # Collect BAMs
    BAMLIST=${BAMLIST}","${OUTP}/${BAMID}.final.bam
done

# Freebayes
${FREEBAYES} --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.15 --fasta-reference ${GENOME} --genotype-qualities `echo ${BAMLIST} | sed 's/,/ -b /g'` -v ${OUTP}/${BAMID}.vcf
${BGZIP} ${OUTP}/${BAMID}.vcf
${TABIX} ${OUTP}/${BAMID}.vcf.gz

# Fixed threshold filtering
${BCF} filter -O b -o ${OUTP}/${BAMID}.bcf -e '%QUAL<=20 || %QUAL/AO<=2 || SAF<=1 || SAR<=1 || RPR<=1 || RPL<=1' ${OUTP}/${BAMID}.vcf.gz
rm ${OUTP}/${BAMID}.vcf.gz ${OUTP}/${BAMID}.vcf.gz.tbi

# Normalize InDels
${BCF} norm -O z -o ${OUTP}/${BAMID}.norm.vcf.gz -c x -f ${GENOME} ${OUTP}/${BAMID}.bcf
${TABIX} ${OUTP}/${BAMID}.norm.vcf.gz
rm ${OUTP}/${BAMID}.bcf

# VEP
perl ${VEP} --species homo_sapiens --assembly GRCh37 --offline --no_progress --no_stats --sift b --polyphen b --ccds --hgvs --symbol --numbers --regulatory --canonical --protein --biotype --tsl --appris --gene_phenotype --af --af_1kg --af_esp --af_exac --pubmed --variant_class --no_escape --vcf --minimal --dir ${VEP_DATA} --fasta ${VEP_DATA}/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --input_file ${OUTP}/${BAMID}.norm.vcf.gz --output_file ${OUTP}/${BAMID}.vep.vcf --plugin ExAC,${VEP_DATA}/ExAC.r0.3.1.sites.vep.vcf.gz --plugin Blosum62 --plugin CSN --plugin Downstream --plugin GO --plugin LoFtool,${BASEDIR}/../bed/LoFtool_scores.txt --plugin TSSDistance
${BGZIP} ${OUTP}/${BAMID}.vep.vcf
${TABIX} ${OUTP}/${BAMID}.vep.vcf.gz

# Convert to BCF
${BCF} view -O b -o ${OUTP}/${BAMID}.vep.bcf ${OUTP}/${BAMID}.vep.vcf.gz
${BCF} index ${OUTP}/${BAMID}.vep.bcf
rm ${OUTP}/${BAMID}.vep.vcf.gz ${OUTP}/${BAMID}.vep.vcf.gz.tbi ${OUTP}/${BAMID}.norm.vcf.gz ${OUTP}/${BAMID}.norm.vcf.gz.tbi

# Clean-up
rm -rf ${TMP}
