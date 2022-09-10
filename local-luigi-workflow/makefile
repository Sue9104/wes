SHELL := /bin/bash
PLATFORM := $(shell uname)
install: python

python:
	$(info "Install conda before installation!!!")
	conda create -y -n wes python=3
	conda activate wes && conda install luigi pandas numpy && conda install -c bioconda bwa samtools bcftools snpEFF multiqc qualimap trim-galore trimmomatic gatk gatktool freebayes picard sra-tools
	conda env export --file wes_env.yml --name wes

snpeff:
	source activate wes && snpEff download -v hg38
	source activate wes && snpEff download -v hg19

annovar:
	sh download_annovar_db.sh
