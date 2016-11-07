# External sources
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
BWASOURCES = ${wildcard src/bwa/*.c) $(wildcard src/bwa/*.h)
SAMSOURCES = $(wildcard src/samtools/*.c) $(wildcard src/samtools/*.h)
BCFSOURCES = $(wildcard src/bcftools/*.c) $(wildcard src/bcftools/*.h)
FREEBAYSOURCES = $(wildcard src/freebayes/src/*.cpp) $(wildcard src/freebayes/src/*.h)
BSTATSSOURCES = $(wildcard src/bamStats/src/*.h)
PBASE=$(shell pwd)

# Targets
TARGETS = .perl .vep .picard .htslib .bwa .samtools .bcftools .freebayes .bamStats

all:   	$(TARGETS)

.perl:
	mkdir src/perl/ && cd src/ && wget 'http://www.cpan.org/src/5.0/perl-5.24.0.tar.gz' && tar -xzf perl-5.24.0.tar.gz && rm perl-5.24.0.tar.gz && cd perl-5.24.0/ && ./Configure -des -Dprefix=${PBASE}/src/perl -Dotherlibdirs=${PBASE}/src/perl/lib/perl5/ && make && make install && cd ../../ && rm -rf src/perl-5.24.0/ && wget --no-check-certificate https://cpanmin.us/ -O src/perl/bin/cpanm && chmod +x src/perl/bin/cpanm && ${PBASE}/src/perl/bin/perl src/perl/bin/cpanm --notest -l ${PBASE}/src/perl App::cpanminus File::Find::Rule JSON::Parse XML::Parser XML::Simple LWP LWP::Simple LWP::Protocol::https Archive::Extract Archive::Tar Archive::Zip CGI DBI Time::HiRes DBD::mysql Encode File::Copy::Recursive Perl::OSType Module::Metadata Statistics::Lite Tie::Autotie Tie::IxHash Log::Log4perl FindBin::Real Getopt::Long Catalyst::Runtime Catalyst::Devel List::Util Test::XML::Simple Test::XPath IO::String Bio::Perl version && touch .perl

.vep: .perl .htslib .samtools .bcftools
	mkdir -p src/vep/cache && cd src/vep/cache/ && wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz && wget ftp://ftp.ensembl.org/pub/release-86/variation/VEP/homo_sapiens_vep_86_GRCh37.tar.gz && tar -izxf homo_sapiens_vep_86_GRCh37.tar.gz && rm homo_sapiens_vep_86_GRCh37.tar.gz && cd ${PBASE}/src/vep/ && wget -O 86.tar.gz https://github.com/Ensembl/ensembl-tools/archive/release/86.tar.gz && tar -zxf 86.tar.gz --starting-file variant_effect_predictor --transform='s|.*/|./|g' && rm 86.tar.gz && export PERL5LIB=${PBASE}/src/perl/lib/perl5/:${PBASE}/src/perl/lib/5.24.0/ && export PATH=${PBASE}/src/htslib:${PBASE}/src/samtools:${PBASE}/src/bcftools:${PBASE}/src/perl/bin:${PATH} && sed -i '0,/ok = <>;/s//ok = "y";/' INSTALL.pl && perl INSTALL.pl -v 86 --NO_HTSLIB --CONVERT --AUTO afp --SPECIES homo_sapiens --ASSEMBLY GRCh37 --PLUGINS ExAC --DESTDIR ${PBASE}/src/vep/ --CACHEDIR ${PBASE}/src/vep/cache && tabix ${PBASE}/src/vep/cache/ExAC.r0.3.1.sites.vep.vcf.gz && cd ../../ && touch .vep

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

clean:
	cd src/htslib && make clean
	cd src/samtools && make clean
	cd src/bcftools && make clean
	cd src/freebayes && make clean
	cd src/bamStats && make clean
	rm -rf $(TARGETS) src/perl/ src/java/ src/vep/ src/picard/
