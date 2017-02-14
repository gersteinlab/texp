# TeXP
TeXP is a pipeline to evaluate the transcription level of transposable elements in short read RNA-seq data

#About
TeXP is a pipeline for quantifying abundances of Transposable Elements transcripts from RNA-Seq data. TeXP is based on the assumption that RNA-seq reads overlapping Transposable Elements is a composition of pervasive transcription signal and autonomous transcription of Transposable Elements.

[[REF]]

# Requirements
 - Bowtie2 ()
 - Fastx-toolkit ()
 - perl ()
 - python ()
 - R ()
  - Penalized package ()
 - samtools ()
 - wgsim (a12da33  on Oct 17, 2011)
 - Bedtools
 
# Download
 $> git clone https://github.com/fabiocpn/TeXP.git

# Docker image


# Running TeXP
 $> make -f Makefile INPUT_FILE_PATH=[FILE_NAME] OUTPUT_DIR=[OUTPUT_PATH] N_THREADS=[INT] REFERENCE_GENOME=[REFERENCE_GENOME_ID] SAMPLE_NAME=[SAMPLE_ID]
 [FILE_NAME]: *.fastq/*.fq.gz/*.sra 
 [OUTPUT_PATH]: i.e. process/ (Path where results will be written)
 [N_THREADS]:  Number of threads TeXP will use
 [REFERENCE_GENOME]: hg19|hg38
 [SAMPLE_ID]: Sample's name.
