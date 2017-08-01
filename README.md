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
 
# Download
 $> git clone https://github.com/fabiocpn/TeXP.git

# Docker image
docker pull fnavarro/texp
https://hub.docker.com/r/fnavarro/texp/ for futher instructions


# Running TeXP
 $> make -f Makefile INPUT_FILE_PATH=[FILE_NAME] OUTPUT_DIR=[OUTPUT_PATH] N_THREADS=[INT] SAMPLE_NAME=[SAMPLE_ID]

 [FILE_NAME]: *.fastq/*.fq.gz/*.sra 

 [OUTPUT_PATH]: i.e. process/ (Path where results will be written)

 [N_THREADS]:  Number of threads TeXP will use

 [SAMPLE_ID]: Sample's name.


