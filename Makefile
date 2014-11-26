# External sources
PICARDSOURCES = $(wildcard src/picard/src/java/picard/*/*.java)

# Targets
TARGETS = .picard

all:   	$(TARGETS)

.picard: $(PICARDSOURCES)
	cd src/picard && ant clone-htsjdk && ant -lib lib/ant/ all && cd ../../ && touch .picard

clean:
	cd src/picard && ant clean
	rm -f .picard
