#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

# GRCh38
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz
zcat GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz > hg38.fa
samtools faidx hg38.fa
rm GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz

# Index
bwa index hg38.fa

# Mappability and exclude maps
wget https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz
wget https://raw.githubusercontent.com/dellytools/delly/main/excludeTemplates/human.hg38.excl.tsv

# Genome bed
cat hg38.fa.fai  | awk '{print $1"\t1\t"$2;}' | grep -P "^chr[0-9XYM]*\t" > hg38.wgs.bed
bgzip hg38.wgs.bed 

# Exome bed
wget http://ftp.ensembl.org/pub/current_gff3/homo_sapiens/Homo_sapiens.GRCh38.110.gff3.gz
wget http://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
zcat Homo_sapiens.GRCh38.*.gtf.gz | grep "ensembl_havana" | grep "exon" | grep "protein_coding" | grep "gene_name" | cut -f 1,4,5,9 | sed 's/\t[^\t]*gene_name "/\t/' | sed 's/".*$//' | uniq | sed 's/^/chr/' | sort -k1,1V -k2,2n | uniq > hg38.wes.bed.tmp
bedtools merge -i hg38.wes.bed.tmp -c 4 -o mode > hg38.wes.bed
rm hg38.wes.bed.tmp
bgzip hg38.wes.bed
