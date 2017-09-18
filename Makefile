SHELL := /bin/bash

# Sources
NREXSOURCES = $(wildcard src/*.cpp)
PBASE=$(shell pwd)

# Targets
TARGETS = .vep .biodbhts .gnomADg .maxentscan .egenome .loftool .exacpli .mpc src/nRex

all:   	$(TARGETS)

.biodbhts:
	module load Perl BioPerl DBD-mysql HTSlib git && cd src/BioDB && perl INSTALL.pl --prefix ${PBASE}/src/BioDbBuild && cd ../../ && touch .biodbhts

.vep:
	module load Perl BioPerl DBD-mysql HTSlib && mkdir -p src/vepcache/ && cd src/vep/ && perl INSTALL.pl -a acp -s homo_sapiens -y GRCh37 -v 90 -c ${PBASE}/src/vepcache -g Blosum62,CSN,Downstream,ExAC,GO,LoFtool,TSSDistance,MaxEntScan,ExACpLI,MPC -t --NO_HTSLIB --NO_BIOPERL && cd ../../ && touch .vep

.maxentscan:
	mkdir -p src/vepcache/Plugins/ && cd src/vepcache/Plugins/ && wget 'http://genes.mit.edu/burgelab/maxent/download/fordownload.tar.gz' && tar -xzf fordownload.tar.gz && rm fordownload.tar.gz && mv fordownload/ maxentscan/ && cd ../../../ && touch .maxentscan

.gnomADg:
	cd src/vepcache/ && wget --no-check-certificate 'ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz' && wget --no-check-certificate 'ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz.tbi' && cd ../../ && touch .gnomADg

.egenome:
	module load SAMtools && cd src/vepcache/ && wget --no-check-certificate 'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz' && gunzip Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz && samtools faidx Homo_sapiens.GRCh37.75.dna.primary_assembly.fa && cd ../../ && touch .egenome

.loftool:
	mkdir -p src/vepcache/Plugins/ && cd src/vepcache/Plugins/ && wget 'https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/90/LoFtool_scores.txt' && cd ../../../ && touch .loftool

.exacpli:
	mkdir -p src/vepcache/Plugins/ && cd src/vepcache/Plugins/ && wget 'https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/90/ExACpLI_values.txt' && cd ../../../ && touch .exacpli

.mpc:
	mkdir -p src/vepcache/Plugins/ && cd src/vepcache/Plugins/ && wget 'ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/regional_missense_constraint/fordist_constraint_official_mpc_values_v2.txt.gz' && cd ../../../ && touch .mpc

src/nRex: $(NREXSOURCES)
	${PBASE}/src/build_nRex

clean:
	rm -rf $(TARGETS) src/vepcache/
