# TeXP
TeXP is a pipeline to evaluate the transcription level of transposable elements in short read RNA-seq data

#About
TeXP is a pipeline for quantifying abundances of Transposable Elements transcripts from RNA-Seq data. TeXP is based on the assumption that RNA-seq reads overlapping Transposable Elements is a composition of pervasive transcription signal and autonomous transcription of Transposable Elements.

[[REF]]

# Requirements
 - Bowtie2 (2.3+)
 - Bedtools (2.26+)
 - Fastx-toolkit (0.0.14+)
 - perl (5.24+)
 - python (2.7)
 - R (3.3+)
  - Penalized package (0.49+)
 - samtools (1.3+)
 - wgsim (a12da33 on Oct 17, 2011)
---
 - Bowtie2 hg38 reference index (http://homes.gersteinlab.org/people/fn64/TeXP/rep_annotation.hg38.tar.bz2)
 - Hg38 repetitive element annotation (http://homes.gersteinlab.org/people/fn64/TeXP/rep_annotation.hg38.tar.bz2)
 
# Download
 $> git clone https://github.com/fabiocpn/TeXP.git

 Edit TeXP.sh and Update INSTALL_DIR variable to the path where TeXP was cloned 

# Docker image
docker pull fnavarro/texp
https://hub.docker.com/r/fnavarro/texp/ for futher instructions

# TeXP config
 A few paramaters must be set in order to work outside the docker enviroment; Parameters are set on opts.mk

 LIBRARY_PATH: Should be pointing to TeXP library
 EXT_LIBRARY_PATH: Path containing the bowtie2 reference index and Transposable element annotation bed file
 EXE_DIR: Path to binaries


# Running TeXP
 $> ./TeXP.sh -f=[FILE_NAME] -t=[INT] -o=[OUTPUT_PATH] n=[SAMPLE_ID]

 -f: Input file (fastq,fastq.gz,sra)

 -t: Number of threads

 -o: Output path (i.e. ./ or ./processed)

 -n: Sample name (i.e. SAMPLE01)
