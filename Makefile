SHELL := /bin/bash

# Targets
TARGETS = .mamba .tools .check .shapeit .vep
PBASE=$(shell pwd)

all: ${TARGETS}

.mamba:
	curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(shell uname)-$(shell uname -m).sh" && bash Mambaforge-$(shell uname)-$(shell uname -m).sh -b -p mamba && rm "Mambaforge-$(shell uname)-$(shell uname -m).sh" && touch .mamba

.tools: .mamba
	export PATH=${PBASE}/mamba/bin:${PATH} && mamba create -y -c conda-forge -c bioconda -n align samtools bcftools bedtools htslib bwa trim-galore fastqc delly alfred freebayes wally && touch .tools

.check: .mamba .tools
	export PATH=${PBASE}/mamba/bin:${PATH} && source activate align && samtools --version && bcftools --version && bedtools --version && bgzip --version && tabix --version && trim_galore --version && delly --version && echo "Installation complete!" && touch .check

.shapeit: .mamba .tools
	wget 'https://github.com/odelaneau/shapeit4/archive/refs/tags/v4.2.2.tar.gz' && tar -xzf v4.2.2.tar.gz && rm v4.2.2.tar.gz && cd shapeit4-4.2.2/ && make all && cd .. && cd shapeit4-4.2.2/maps/ && tar -xzf genetic_maps.b38.tar.gz && cd ../../ && touch .shapeit

.vep: .mamba .tools
	mkdir vep_data && chmod a+rwx vep_data && docker run -t -i -v ${PBASE}/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep INSTALL.pl -a cfp -s homo_sapiens -y GRCh38 -g all && touch .vep

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) conda/ pangolin/

