SHELL := /bin/bash

# Targets
TARGETS = .check .conda .mamba .tools
PBASE=$(shell pwd)

all: ${TARGETS}

.conda:
	wget 'https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh' && bash Miniconda3-latest-Linux-x86_64.sh -b -p ${PBASE}/conda && rm -f Miniconda3-latest-Linux-x86_64.sh && touch .conda

.mamba: .conda
	export PATH=${PBASE}/conda/bin:${PATH} && conda install -y -n base -c conda-forge mamba && touch .mamba

.tools: .conda .mamba
	export PATH=${PBASE}/conda/bin:${PATH} && source activate base && mamba create -y -c conda-forge -c bioconda -n align samtools bcftools bedtools htslib bwa trim-galore fastqc delly alfred freebayes wally && touch .tools

.check: .conda .mamba .tools
	export PATH=${PBASE}/conda/bin:${PATH} && source activate align && samtools --version && bcftools --version && bedtools --version && bgzip --version && tabix --version && trim_galore --version && delly --version && echo "Installation complete!" && touch .check

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) conda/ pangolin/
