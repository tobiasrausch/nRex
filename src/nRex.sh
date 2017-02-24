#!/bin/bash

if [ $# -lt 4 ]
then
    echo "**********************************************************************"
    echo "nRex: Single-nucleotide variant calling."
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "nRex (Version: 0.0.2)"
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <wgs|wex|haloplex> <genome.fa> <output prefix> <sample1.read1.fq.gz> <sample1.read2.fq.gz> ..."
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

export PERL5LIB=${BASEDIR}/perl/lib/perl5/:${BASEDIR}/perl/lib/5.24.0/
export PATH=${BASEDIR}/perl/bin:${PATH}

# CMD parameters
ATYPE=${1}
GENOME=${2}
OUTP=${3}
shift
shift
shift
THREADS=4

# Programs
PICARD=${BASEDIR}/picard/picard.jar
FASTQC=${BASEDIR}/FastQC/fastqc
FREEBAYES=${BASEDIR}/freebayes/bin/freebayes
VEP=${BASEDIR}/vep/vep.pl
VEP_DATA=${BASEDIR}/vepcache

# Tmp directory
DSTR=$(date +'%a_%y%m%d_%H%M')
if [ -n "${SCRATCHDIR}" ]
then
    export TMP=${SCRATCHDIR}
    echo "scratch directory" ${SCRATCHDIR}
else
    #export TMP=/tmp/tmp_atac_${DSTR}
    export TMP=/g/solexa/home/rausch/scripts/cpp/nRex/tmp/tmp_nrex_${DSTR}
    mkdir -p ${TMP}
fi
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
    FQ1ID=`echo ${OUTP} | sed 's/$/.fq${i}/'`
    FQ2ID=`echo ${OUTP} | sed 's/$/.fq${j}/'`
    BAMID=`echo ${OUTP} | sed "s/$/.align${i}/"`

    # Fastqc
    mkdir -p ${OUTP}/prefastqc/ && ${FASTQC} -t ${THREADS} -o ${OUTP}/prefastqc/ ${SAMPLEARR[$i]} && ${FASTQC} -t ${THREADS} -o ${OUTP}/prefastqc/ ${SAMPLEARR[$j]}
    
    # BWA
    bwa mem -t ${THREADS} -R "@RG\tID:${OUTP}\tSM:${OUTP}" ${GENOME} ${SAMPLEARR[$i]} ${SAMPLEARR[$j]} | samtools view -bT ${GENOME} - > ${OUTP}/${BAMID}.bam

    # Sort & Index
    samtools sort -o ${OUTP}/${BAMID}.srt.bam ${OUTP}/${BAMID}.bam && rm ${OUTP}/${BAMID}.bam && samtools index ${OUTP}/${BAMID}.srt.bam

    # Clean .bam file
    java ${JAVAOPT} -jar ${PICARD} CleanSam I=${OUTP}/${BAMID}.srt.bam O=${OUTP}/${BAMID}.srt.clean.bam ${PICARDOPT} && rm ${OUTP}/${BAMID}.srt.bam*

    # Mark duplicates
    if [ ${ATYPE} == "haloplex" ]
    then
	mv ${OUTP}/${BAMID}.srt.clean.bam ${OUTP}/${BAMID}.srt.clean.rmdup.bam && samtools index ${OUTP}/${BAMID}.srt.clean.rmdup.bam
    else
	java ${JAVAOPT} -jar ${PICARD} MarkDuplicates I=${OUTP}/${BAMID}.srt.clean.bam O=${OUTP}/${BAMID}.srt.clean.rmdup.bam M=${OUTP}/${OUTP}.markdups.log PG=null MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 ${PICARDOPT} && rm ${OUTP}/${BAMID}.srt.clean.bam* && ${SAM} index ${OUTP}/${BAMID}.srt.clean.rmdup.bam
    fi
    
    # Run stats using unfiltered BAM
    samtools idxstats ${OUTP}/${BAMID}.srt.clean.rmdup.bam > ${OUTP}/${OUTP}.idxstats
    samtools flagstat ${OUTP}/${BAMID}.srt.clean.rmdup.bam > ${OUTP}/${OUTP}.flagstat

    # Filter duplicates, unmapped reads, chrM and unplaced contigs
    if [ ${ATYPE} == "haloplex" ]
    then
	CHRS=`cat ${BASEDIR}/../bed/haloplex.bed | cut -f 1 | sort -k1,1V -k2,2n | uniq | tr '\n' ' '`
    else
	CHRS=`cat ${BASEDIR}/../bed/exome.bed | cut -f 1 | sort -k1,1V -k2,2n | uniq | tr '\n' ' '`
    fi
    samtools view -F 1024 -b ${OUTP}/${BAMID}.srt.clean.rmdup.bam ${CHRS} > ${OUTP}/${BAMID}.final.bam
    samtools index ${OUTP}/${BAMID}.final.bam
    rm ${OUTP}/${BAMID}.srt.clean.rmdup.bam ${OUTP}/${BAMID}.srt.clean.rmdup.bam.bai

    # Run stats using filtered BAM, TSS enrichment, error rates, etc.
    if [ ${ATYPE} == "haloplex" ]
    then
	alfred -b ${BASEDIR}/../bed/haloplex.bed -r ${GENOME} -o ${OUTP}/${OUTP}.bamStats ${OUTP}/${BAMID}.final.bam
    else
	alfred -b ${BASEDIR}/../bed/exome.bed -r ${GENOME} -o ${OUTP}/${OUTP}.bamStats ${OUTP}/${BAMID}.final.bam
    fi

    # Collect BAMs
    BAMLIST=${BAMLIST}","${OUTP}/${BAMID}.final.bam
done

# Freebayes
${FREEBAYES} --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.15 --fasta-reference ${GENOME} --genotype-qualities `echo ${BAMLIST} | sed 's/,/ -b /g'` -v ${OUTP}/${BAMID}.vcf
bgzip ${OUTP}/${BAMID}.vcf
tabix ${OUTP}/${BAMID}.vcf.gz

# Normalize VCF
vt normalize ${OUTP}/${BAMID}.vcf.gz -r ${GENOME} | vt decompose_blocksub - | vt decompose - | vt uniq - | bgzip > ${OUTP}/${BAMID}.norm.vcf.gz
tabix ${OUTP}/${BAMID}.norm.vcf.gz
rm ${OUTP}/${BAMID}.vcf.gz ${OUTP}/${BAMID}.vcf.gz.tbi

# Fixed threshold filtering
if [ ${ATYPE} == "haloplex" ]
then
    bcftools filter -O z -o ${OUTP}/${BAMID}.norm.filtered.vcf.gz -e '%QUAL<=20 || %QUAL/AO<=2 || SAF<=1 || SAR<=1' ${OUTP}/${BAMID}.norm.vcf.gz
else
    bcftools filter -O z -o ${OUTP}/${BAMID}.norm.filtered.vcf.gz -e '%QUAL<=20 || %QUAL/AO<=2 || SAF<=1 || SAR<=1 || RPR<=1 || RPL<=1' ${OUTP}/${BAMID}.norm.vcf.gz
fi
tabix ${OUTP}/${BAMID}.norm.filtered.vcf.gz
rm ${OUTP}/${BAMID}.norm.vcf.gz ${OUTP}/${BAMID}.norm.vcf.gz.tbi

# VEP
perl ${VEP} --species homo_sapiens --assembly GRCh37 --offline --no_progress --no_stats --sift b --polyphen b --ccds --hgvs --symbol --numbers --regulatory --canonical --protein --biotype --tsl --appris --gene_phenotype --af --af_1kg --af_esp --af_exac --pubmed --variant_class --no_escape --vcf --minimal --dir ${VEP_DATA} --fasta ${VEP_DATA}/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --input_file ${OUTP}/${BAMID}.norm.filtered.vcf.gz --output_file ${OUTP}/${BAMID}.vep.vcf --plugin ExAC,${VEP_DATA}/ExAC.r0.3.1.sites.vep.vcf.gz --plugin Blosum62 --plugin CSN --plugin Downstream --plugin GO --plugin LoFtool,${BASEDIR}/../bed/LoFtool_scores.txt --plugin TSSDistance --plugin MaxEntScan,${VEP_DATA}/Plugins/maxentscan/
bgzip ${OUTP}/${BAMID}.vep.vcf
tabix ${OUTP}/${BAMID}.vep.vcf.gz

# Convert to BCF
bcftools view -O b -o ${OUTP}/${BAMID}.vep.bcf ${OUTP}/${BAMID}.vep.vcf.gz
bcftools index ${OUTP}/${BAMID}.vep.bcf
rm ${OUTP}/${BAMID}.vep.vcf.gz ${OUTP}/${BAMID}.vep.vcf.gz.tbi ${OUTP}/${BAMID}.norm.filtered.vcf.gz ${OUTP}/${BAMID}.norm.filtered.vcf.gz.tbi

# Clean-up
if [ -n "${SCRATCHDIR}" ]
then
    ls ${SCRATCHDIR}
else
    rm -rf ${TMP}
fi
