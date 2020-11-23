SHELL := /bin/bash
DEBUG ?= 0
STATIC ?= 0

# Submodules
PWD = $(shell pwd)
EBROOTHTSLIB ?= ${PWD}/src/htslib/

# Sources
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
NREXSOURCES = $(wildcard src/*.h) $(wildcard src/*.cpp)
PBASE=$(shell pwd)

# Flags
CXX=g++
CXXFLAGS += -isystem ${EBROOTHTSLIB} -pedantic -W -Wall -Wno-unknown-pragmas -D__STDC_LIMIT_MACROS -fno-strict-aliasing -fpermissive
LDFLAGS += -L${EBROOTHTSLIB} -L${EBROOTHTSLIB}/lib -lboost_iostreams -lboost_filesystem -lboost_system -lboost_program_options -lboost_date_time 

# Flags for static compile
ifeq (${STATIC}, 1)
	LDFLAGS += -static -static-libgcc -pthread -lhts -lz -llzma -lbz2
else
	LDFLAGS += -lhts -lz -llzma -lbz2 -Wl,-rpath,${EBROOTHTSLIB}
endif

# Flags for debugging, profiling and releases
ifeq (${DEBUG}, 1)
	CXXFLAGS += -g -O0 -fno-inline -DDEBUG
else ifeq (${DEBUG}, 2)
	CXXFLAGS += -g -O0 -fno-inline -DPROFILE
	LDFLAGS += -lprofiler -ltcmalloc
else
	CXXFLAGS += -O3 -fno-tree-vectorize -DNDEBUG
endif
ifeq (${EBROOTHTSLIB}, ${PWD}/src/htslib/)
	SUBMODULES += .htslib
endif


# Targets
BUILT_PROGRAMS = src/nRex
TARGETS = ${SUBMODULES} ${BUILT_PROGRAMS} .vep .biodbhts .gnomADg .maxentscan .egenome .loftool .exacpli .mpc src/nRex

all:   	$(TARGETS)

.biodbhts:
	module load Perl BioPerl DBD-mysql HTSlib git && cd src/BioDB && perl INSTALL.pl --prefix ${PBASE}/src/BioDbBuild && cd ../../ && touch .biodbhts

.vep:
	module load Perl BioPerl DBD-mysql HTSlib && mkdir -p src/vepcache/ && cd src/vep/ && perl INSTALL.pl -a acp -s homo_sapiens -y GRCh37 -v 90 -c ${PBASE}/src/vepcache -g Blosum62,CSN,Downstream,ExAC,GO,LoFtool,TSSDistance,MaxEntScan,ExACpLI,MPC -t --NO_HTSLIB --NO_BIOPERL && cd ../../ && touch .vep

.maxentscan:
	mkdir -p src/vepcache/Plugins/ && cd src/vepcache/Plugins/ && wget 'http://genes.mit.edu/burgelab/maxent/download/fordownload.tar.gz' && tar -xzf fordownload.tar.gz && rm fordownload.tar.gz && mv fordownload/ maxentscan/ && cd ../../../ && touch .maxentscan

.gnomADg:
	module load HTSlib BCFtools && cd src/vepcache/ && wget --no-check-certificate 'ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz' && wget --no-check-certificate 'ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz.tbi' && bcftools norm -O z -o gnomad.tmp.vcf.gz -f Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -m -both gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz && tabix gnomad.tmp.vcf.gz && mv gnomad.tmp.vcf.gz gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz && mv gnomad.tmp.vcf.gz.tbi gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz.tbi && cd ../../ && touch .gnomADg

.egenome:
	module load SAMtools && cd src/vepcache/ && wget --no-check-certificate 'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz' && gunzip Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz && samtools faidx Homo_sapiens.GRCh37.75.dna.primary_assembly.fa && cd ../../ && touch .egenome

.loftool:
	mkdir -p src/vepcache/Plugins/ && cd src/vepcache/Plugins/ && wget 'https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/90/LoFtool_scores.txt' && cd ../../../ && touch .loftool

.exacpli:
	mkdir -p src/vepcache/Plugins/ && cd src/vepcache/Plugins/ && wget 'https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/90/ExACpLI_values.txt' && cd ../../../ && touch .exacpli

.mpc:
	mkdir -p src/vepcache/Plugins/ && cd src/vepcache/Plugins/ && wget 'ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/regional_missense_constraint/fordist_constraint_official_mpc_values_v2.txt.gz' && cd ../../../ && touch .mpc

.htslib: $(HTSLIBSOURCES)
	if [ -r src/htslib/Makefile ]; then cd src/htslib && autoheader && autoconf && ./configure --disable-s3 --disable-gcs --disable-libcurl --disable-plugins && $(MAKE) && $(MAKE) lib-static && cd ../../ && touch .htslib; fi

src/nRex: ${SUBMODULES} $(NREXSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

clean:
	if [ -r src/htslib/Makefile ]; then cd src/htslib && $(MAKE) clean; fi
	rm -rf $(TARGETS) src/vepcache/
