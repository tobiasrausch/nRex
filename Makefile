SHELL := /bin/bash

CXX=g++
CXXFLAGS += -isystem ${EBROOTHTSLIB} -pedantic -W -Wall -Wno-unknown-pragmas -fno-strict-aliasing
LDFLAGS += -L${EBROOTHTSLIB}
LDFLAGS += -lhts -lz -Wl,-rpath,${EBROOTHTSLIB}
CXXFLAGS += -O3 -fno-tree-vectorize -DNDEBUG -D__STDC_LIMIT_MACROS

# External sources
NREXSOURCES = $(wildcard src/*.cpp)
PBASE=$(shell pwd)

# Targets
TARGETS = .vep .exac .maxentscan .egenome src/nRex

all:   	$(TARGETS)

.vep:
	module load Perl BioPerl DBD-mysql HTSlib && mkdir -p src/vepcache/ && cd src/vep/ && perl INSTALL.pl -a acp -s homo_sapiens -y GRCh37 -v 90 -c ${PBASE}/src/vepcache -g Blosum62,CSN,Downstream,ExAC,GO,LoFtool,TSSDistance,MaxEntScan,RankFilter -t --NO_HTSLIB --NO_BIOPERL && cd ../../ && touch .vep

.maxentscan: .vep
	cd src/vepcache/Plugins/ && wget 'http://genes.mit.edu/burgelab/maxent/download/fordownload.tar.gz' && tar -xzf fordownload.tar.gz && rm fordownload.tar.gz && mv fordownload/ maxentscan/ && cd ../../../ && touch .maxentscan

.exac: .vep
	cd src/vepcache/ && wget --no-check-certificate 'ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz' && wget --no-check-certificate 'ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz.tbi' && cd ../../ && touch .exac

.egenome: .vep
	module load SAMtools && cd src/vepcache/ && wget --no-check-certificate 'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz' && gunzip Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz && samtools faidx Homo_sapiens.GRCh37.75.dna.primary_assembly.fa && cd ../../ && touch .egenome

src/nRex: $(NREXSOURCES)
	module load foss Boost HTSlib && $(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

clean:
	rm -rf $(TARGETS) src/vepcache/
