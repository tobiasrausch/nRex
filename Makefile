SHELL := /bin/bash

# Targets
TARGETS = .mamba .sr .lr .clair3 .check .vep
PBASE=$(shell pwd)

all: ${TARGETS}

.mamba:
	curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(shell uname)-$(shell uname -m).sh" && bash Miniforge3-$(shell uname)-$(shell uname -m).sh -b -p mamba && rm "Miniforge3-$(shell uname)-$(shell uname -m).sh" && touch .mamba

.sr: .mamba
	export PATH=${PBASE}/mamba/bin:${PATH} && mamba create -y -n shortread --override-channels -c conda-forge -c bioconda samtools bcftools bedtools htslib bwa fastp fastqc delly alfred freebayes wally && touch .sr

.lr: .mamba
	export PATH=${PBASE}/mamba/bin:${PATH} && mamba create -y -n longread --override-channels -c conda-forge -c bioconda samtools bcftools bedtools htslib delly alfred wally minimap2 shapeit4 nanoplot sniffles whatshap longshot && touch .lr

.check: .mamba .sr .lr
	export PATH=${PBASE}/mamba/bin:${PATH} && source activate shortread && samtools --version && bcftools --version && bedtools --version && bgzip --version && tabix --version && delly --version && source deactivate && source activate longread && sniffles --version && NanoPlot --version && echo "Installation complete!" && touch .check

.clair3: .mamba
	export PATH=${PBASE}/mamba/bin:${PATH} && mamba create -y -n clair3 --override-channels -c conda-forge -c bioconda clair3 python=3.9.0 && touch .clair3

.vep: .mamba
	mkdir vep_data && chmod a+rwx vep_data && docker run -t -i -v ${PBASE}/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep INSTALL.pl -a cfp -s homo_sapiens -y GRCh38 -g all && touch .vep

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) mamba/

