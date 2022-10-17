SHELL := /bin/bash
PLATFORM := $(shell uname)
install: python

python:
	$(info "Install conda before installation!!!")
	mamba create -y -n wes python==3.7
	mamba activate wes && mamba install luigi pandas numpy pyyaml && mamba install -c bioconda bwa samtools bcftools snpeff multiqc qualimap trim-galore trimmomatic gatk gatktool freebayes picard sra-tools
	mamba env export --file wes_env.yml --name wes

snpeff:
	source activate wes && snpEff download -v hg38
	source activate wes && snpEff download -v hg19

annovar:
	sh download_annovar_db.sh
