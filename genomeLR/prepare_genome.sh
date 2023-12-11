#!/bin/bash

# GRCh38
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz
zcat GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz > hg38.fa
samtools faidx hg38.fa
rm GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz

# Genome bed
cat hg38.fa.fai  | awk '{print $1"\t1\t"$2;}' | grep -P "^chr[0-9XYM]*\t" > hg38.wgs.bed
bgzip hg38.wgs.bed 

# Mappability map
wget https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz
wget https://gear-genomics.embl.de/data/delly/common_sites.vcf.gz
tabix common_sites.vcf.gz

# T2T
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz
gunzip chm13v2.0_maskedY_rCRS.fa.gz
mv chm13v2.0_maskedY_rCRS.fa t2t.fa
samtools faidx t2t.fa
rm chm13v2.0_maskedY_rCRS.fa
