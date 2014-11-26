# External sources
PICARDSOURCES = $(wildcard src/picard/src/java/picard/*/*.java)
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
SAMSOURCES = $(wildcard src/samtools/*.c) $(wildcard src/samtools/*.h)
BCFSOURCES = $(wildcard src/bcftools/*.c) $(wildcard src/bcftools/*.h)
FREEBAYSOURCES = $(wildcard src/freebayes/src/*.cpp) $(wildcard src/freebayes/src/*.h)
BAMRGSOURCES = $(wildcard src/bamaddrg/*.cpp)

# Targets
TARGETS = .picard .htslib .samtools .bcftools .freebayes .bamaddrg

all:   	$(TARGETS)

.picard: $(PICARDSOURCES)
	cd src/picard && ant clone-htsjdk && ant -lib lib/ant/ all && cd ../../ && touch .picard

.htslib: $(HTSLIBSOURCES)
	cd src/htslib && make && cd ../../ && touch .htslib

.samtools: $(SAMSOURCES)
	cd src/samtools && make && cd ../../ && touch .samtools

.bcftools: $(BCFSOURCES)
	cd src/bcftools && make && cd ../../ && touch .bcftools

.freebayes: $(FREEBAYSOURCES)
	cd src/freebayes && make && cd ../../ && touch .freebayes

.bamaddrg: $(BAMRGSOURCES)
	cd src/bamaddrg && make && cd ../../ && touch .bamaddrg

clean:
	cd src/picard && ant clean
	cd src/htslib && make clean
	cd src/samtools && make clean
	cd src/bcftools && make clean
	cd src/freebayes && make clean
	cd src/bamaddrg && make clean
	rm -f .picard .htslib .samtools .bcftools .freebayes .bamaddrg
