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

.shapeit: .conda .mamba .tools
	wget 'https://github.com/odelaneau/shapeit4/archive/refs/tags/v4.2.2.tar.gz' && tar -xzf v4.2.2.tar.gz && rm v4.2.2.tar.gz && cd shapeit4-4.2.2/ && make all && cd .. && cd shapeit4-4.2.2/maps/ && tar -xzf genetic_maps.b38.tar.gz && cd ../../ && touch .shapeit

.vep: .conda .mamba .tools
	mkdir vep_data && chmod a+rwx vep_data && docker run -t -i -v ${PBASE}/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep INSTALL.pl -a cfp -s homo_sapiens -y GRCh38 -g all && touch .vep

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) conda/ pangolin/
