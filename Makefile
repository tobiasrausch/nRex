SHELL := /bin/bash

# Targets
TARGETS = .mamba .tools .check .shapeit .vep
PBASE=$(shell pwd)

all: ${TARGETS}

.mamba:
	curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(shell uname)-$(shell uname -m).sh" && bash Mambaforge-$(shell uname)-$(shell uname -m).sh -b -p mamba && rm "Mambaforge-$(shell uname)-$(shell uname -m).sh" && touch .mamba

.tools: .mamba
	export PATH=${PBASE}/mamba/bin:${PATH} && mamba install -y -c conda-forge -c bioconda samtools bcftools bedtools htslib bwa trim-galore fastqc delly alfred freebayes wally minimap2 shapeit4 && pip install nanopack sniffles && touch .tools

.check: .mamba .tools
	export PATH=${PBASE}/mamba/bin:${PATH} && samtools --version && bcftools --version && bedtools --version && bgzip --version && tabix --version && trim_galore --version && delly --version && sniffles --version && echo "Installation complete!" && touch .check

.vep: .mamba .tools
	mkdir vep_data && chmod a+rwx vep_data && docker run -t -i -v ${PBASE}/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep INSTALL.pl -a cfp -s homo_sapiens -y GRCh38 -g all && touch .vep

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) mamba/

