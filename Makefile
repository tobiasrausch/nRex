SHELL := /bin/bash

# Sources
NREXSOURCES = $(wildcard src/*.cpp)
PBASE=$(shell pwd)

# Targets
TARGETS = .vep .gnomADg .maxentscan .egenome src/nRex

all:   	$(TARGETS)

.vep:
	module load Perl BioPerl DBD-mysql HTSlib && mkdir -p src/vepcache/ && cd src/vep/ && perl INSTALL.pl -a acp -s homo_sapiens -y GRCh37 -v 90 -c ${PBASE}/src/vepcache -g Blosum62,CSN,Downstream,ExAC,GO,LoFtool,TSSDistance,MaxEntScan,RankFilter -t --NO_HTSLIB --NO_BIOPERL && cd ../../ && touch .vep

.maxentscan: .vep
	mkdir -p src/vepcache/Plugins/ && cd src/vepcache/Plugins/ && wget 'http://genes.mit.edu/burgelab/maxent/download/fordownload.tar.gz' && tar -xzf fordownload.tar.gz && rm fordownload.tar.gz && mv fordownload/ maxentscan/ && cd ../../../ && touch .maxentscan

.gnomADg: .vep
	cd src/vepcache/ && wget --no-check-certificate 'ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz' && wget --no-check-certificate 'ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz.tbi' && cd ../../ && touch .gnomADg

.egenome: .vep
	module load SAMtools && cd src/vepcache/ && wget --no-check-certificate 'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz' && gunzip Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz && samtools faidx Homo_sapiens.GRCh37.75.dna.primary_assembly.fa && cd ../../ && touch .egenome

src/nRex: $(NREXSOURCES)
	${PBASE}/src/build_nRex

clean:
	rm -rf $(TARGETS) src/vepcache/
