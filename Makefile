SHELL := /bin/bash

CXX=g++
CXXFLAGS += -isystem ${EBROOTHTSLIB} -pedantic -W -Wall -Wno-unknown-pragmas -fno-strict-aliasing
LDFLAGS += -L${EBROOTHTSLIB}
LDFLAGS += -lhts -lz -Wl,-rpath,${EBROOTHTSLIB}
CXXFLAGS += -O3 -fno-tree-vectorize -DNDEBUG -D__STDC_LIMIT_MACROS

# External sources
FREEBAYSOURCES = $(wildcard src/freebayes/src/*.cpp) $(wildcard src/freebayes/src/*.h)
NREXSOURCES = $(wildcard src/*.cpp)
PBASE=$(shell pwd)

# Targets
TARGETS = .fastqc .perl .vep .exac .maxentscan .egenome .picard .freebayes src/nRex

all:   	$(TARGETS)

.fastqc:
	module load Java && cd src && wget 'http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip' && unzip fastqc_v0.11.5.zip && chmod 755 FastQC/fastqc && rm fastqc_v0.11.5.zip && cd ../ && touch .fastqc

.perl:
	module load foss expat && mkdir src/perl/ && cd src/ && wget 'http://www.cpan.org/src/5.0/perl-5.24.0.tar.gz' && tar -xzf perl-5.24.0.tar.gz && rm perl-5.24.0.tar.gz && cd perl-5.24.0/ && ./Configure -des -Dprefix=${PBASE}/src/perl -Dotherlibdirs=${PBASE}/src/perl/lib/perl5/ && make && make install && cd ../../ && rm -rf src/perl-5.24.0/ && wget --no-check-certificate https://cpanmin.us/ -O src/perl/bin/cpanm && chmod +x src/perl/bin/cpanm && ${PBASE}/src/perl/bin/perl src/perl/bin/cpanm --notest -l ${PBASE}/src/perl App::cpanminus File::Find::Rule JSON JSON::Parse XML::Parser Set::IntervalTree XML::Simple LWP LWP::Simple LWP::Protocol::https Archive::Extract Archive::Tar Archive::Zip CGI DBI Time::HiRes DBD::mysql Encode File::Copy::Recursive PerlIO::gzip Perl::OSType Module::Metadata Statistics::Lite Tie::Autotie Tie::IxHash Log::Log4perl FindBin::Real Getopt::Long Catalyst::Runtime Catalyst::Devel List::Util Test::XML::Simple Test::XPath IO::String Bio::Perl version && touch .perl

.vep: .perl
	mkdir src/vepcache/ && cd src/vep/ && export PERL5LIB=${PBASE}/src/perl/lib/perl5/:${PBASE}/src/perl/lib/5.24.0/ && export PATH=${PBASE}/src/perl/bin:${PATH} && perl INSTALL.pl -a acp -s homo_sapiens -y GRCh37 -v 87 -c ${PBASE}/src/vepcache -g Blosum62,CSN,Downstream,ExAC,GO,LoFtool,TSSDistance,MaxEntScan,RankFilter -t && cd ../../ && touch .vep

.maxentscan: .perl .vep
	cd src/vepcache/Plugins/ && wget 'http://genes.mit.edu/burgelab/maxent/download/fordownload.tar.gz' && tar -xzf fordownload.tar.gz && rm fordownload.tar.gz && mv fordownload/ maxentscan/ && cd ../../../ && touch .maxentscan

.exac: .perl .vep
	cd src/vepcache/ && wget --no-check-certificate 'ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz' && wget --no-check-certificate 'ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz.tbi' && cd ../../ && touch .exac

.egenome: .perl .vep
	module load SAMtools && cd src/vepcache/ && wget --no-check-certificate 'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz' && gunzip Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz && samtools faidx Homo_sapiens.GRCh37.75.dna.primary_assembly.fa && cd ../../ && touch .egenome

.picard:
	module load Java && mkdir -p src/picard/ && cd src/picard && wget -O picard.jar 'https://github.com/broadinstitute/picard/releases/download/2.8.3/picard.jar' && cd ../../ && touch .picard

.freebayes: $(FREEBAYSOURCES)
	module load foss HTSlib CMake && cd src/freebayes && make && cd ../../ && touch .freebayes

src/nRex: $(NREXSOURCES)
	module load Boost HTSlib && $(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

clean:
	cd src/freebayes && make clean
	rm -rf $(TARGETS) src/perl/ src/vepcache/ src/picard/ src/FastQC
