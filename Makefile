SEQTK_ROOT ?= ${PWD}/src/htslib/

CXX=g++
CXXFLAGS += -isystem ${SEQTK_ROOT} -pedantic -W -Wall -Wno-unknown-pragmas -fno-strict-aliasing
LDFLAGS += -L${SEQTK_ROOT}
LDFLAGS += -lhts -lz -Wl,-rpath,${SEQTK_ROOT}
CXXFLAGS += -O3 -fno-tree-vectorize -DNDEBUG -D__STDC_LIMIT_MACROS

# External sources
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
BWASOURCES = ${wildcard src/bwa/*.c) $(wildcard src/bwa/*.h)
SAMSOURCES = $(wildcard src/samtools/*.c) $(wildcard src/samtools/*.h)
BCFSOURCES = $(wildcard src/bcftools/*.c) $(wildcard src/bcftools/*.h)
FREEBAYSOURCES = $(wildcard src/freebayes/src/*.cpp) $(wildcard src/freebayes/src/*.h)
BSTATSSOURCES = $(wildcard src/bamStats/src/*.h)
NREXSOURCES = $(wildcard src/*.cpp)
PBASE=$(shell pwd)

# Targets
TARGETS = .perl .vep .exac .egenome .picard .htslib .bwa .samtools .bcftools .freebayes .bamStats src/nRex

all:   	$(TARGETS)

.perl:
	mkdir src/perl/ && cd src/ && wget 'http://www.cpan.org/src/5.0/perl-5.24.0.tar.gz' && tar -xzf perl-5.24.0.tar.gz && rm perl-5.24.0.tar.gz && cd perl-5.24.0/ && ./Configure -des -Dprefix=${PBASE}/src/perl -Dotherlibdirs=${PBASE}/src/perl/lib/perl5/ && make && make install && cd ../../ && rm -rf src/perl-5.24.0/ && wget --no-check-certificate https://cpanmin.us/ -O src/perl/bin/cpanm && chmod +x src/perl/bin/cpanm && ${PBASE}/src/perl/bin/perl src/perl/bin/cpanm --notest -l ${PBASE}/src/perl App::cpanminus File::Find::Rule JSON JSON::Parse XML::Parser Set::IntervalTree XML::Simple LWP LWP::Simple LWP::Protocol::https Archive::Extract Archive::Tar Archive::Zip CGI DBI Time::HiRes DBD::mysql Encode File::Copy::Recursive PerlIO::gzip Perl::OSType Module::Metadata Statistics::Lite Tie::Autotie Tie::IxHash Log::Log4perl FindBin::Real Getopt::Long Catalyst::Runtime Catalyst::Devel List::Util Test::XML::Simple Test::XPath IO::String Bio::Perl version && touch .perl

.vep: .perl
	mkdir src/vepcache/ && cd src/vep/ && export PERL5LIB=${PBASE}/src/perl/lib/perl5/:${PBASE}/src/perl/lib/5.24.0/ && export PATH=${PBASE}/src/perl/bin:${PATH} && perl INSTALL.pl -a acp -s homo_sapiens -y GRCh37 -v 87 -c ${PBASE}/src/vepcache -g Blosum62,CSN,Downstream,ExAC,GO,LoFtool,TSSDistance -t && cd ../../ && touch .vep

.exac: .perl .vep
	cd src/vepcache/ && wget --no-check-certificate 'ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz' && wget --no-check-certificate 'ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz.tbi' && cd ../../ && touch .exac

.egenome: .perl .vep .samtools
	cd src/vepcache/ && wget --no-check-certificate 'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz' && gunzip Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz && ${PBASE}/src/samtools/samtools faidx Homo_sapiens.GRCh37.75.dna.primary_assembly.fa && cd ../../ && touch .egenome

.picard:
	mkdir -p src/picard/ && cd src/picard && wget -O picard.jar 'https://github.com/broadinstitute/picard/releases/download/2.7.1/picard.jar' && cd ../../ && touch .picard

.htslib: $(HTSLIBSOURCES)
	cd src/htslib && make && make lib-static && cd ../../ && touch .htslib

.bwa: $(HTSLIBSOURCES)
	cd src/bwa && make && cd ../../ && touch .bwa

.samtools: $(SAMSOURCES)
	cd src/samtools && make && cd ../../ && touch .samtools

.bcftools: $(BCFSOURCES)
	cd src/bcftools && make && cd ../../ && touch .bcftools

.freebayes: $(FREEBAYSOURCES)
	cd src/freebayes && make && cd ../../ && touch .freebayes

.bamStats: $(BSTATSSOURCES)
	cd src/bamStats && make all && cd ../../ && touch .bamStats

src/nRex: $(NREXSOURCES) .htslib
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

clean:
	cd src/htslib && make clean
	cd src/samtools && make clean
	cd src/bcftools && make clean
	cd src/freebayes && make clean
	cd src/bamStats && make clean
	rm -rf $(TARGETS) src/perl/ src/vepcache/ src/picard/
