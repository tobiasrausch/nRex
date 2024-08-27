SHELL := /bin/bash

# Targets
TARGETS = .mamba .tools .clair3 .check .vep
PBASE=$(shell pwd)

all: ${TARGETS}

.mamba:
	curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(shell uname)-$(shell uname -m).sh" && bash Mambaforge-$(shell uname)-$(shell uname -m).sh -b -p mamba && rm "Mambaforge-$(shell uname)-$(shell uname -m).sh" && touch .mamba

.tools: .mamba
	export PATH=${PBASE}/mamba/bin:${PATH} && mamba install -y --override-channels -c conda-forge -c bioconda samtools bcftools bedtools htslib bwa fastp delly alfred freebayes wally minimap2 shapeit4 nanoplot sniffles whatshap longshot && touch .tools

.check: .mamba .tools
	export PATH=${PBASE}/mamba/bin:${PATH} && samtools --version && bcftools --version && bedtools --version && bgzip --version && tabix --version && delly --version && sniffles --version && NanoPlot --version && echo "Installation complete!" && touch .check

.clair3: .mamba .tools
	export PATH=${PBASE}/mamba/bin:${PATH} && mamba create -y -n clair3 --override-channels -c conda-forge -c bioconda clair3 python=3.9.0 && touch .clair3

.vep: .mamba .tools
	mkdir vep_data && chmod a+rwx vep_data && docker run -t -i -v ${PBASE}/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep INSTALL.pl -a cfp -s homo_sapiens -y GRCh38 -g all && touch .vep

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) mamba/

