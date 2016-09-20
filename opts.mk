PIPELINE_NAME = TeXP

DATA_DIR          := NULL
OUTPUT_DIR        := NULL
INPUT_FILE_PATH   := NULL
SAMPLE_NAME       := NULL
REFERENCE_GENOME  := NULL
CONFIGURED        := NULL

##
## Use the input path to infer filetype and short name
##
INPUT_FILE_NAME := $(notdir $(INPUT_FILE_PATH))
INPUT_FILE_ID   := $(basename $(INPUT_FILE_NAME))

LIBRARY_PATH     := /scratch/fas/gerstein/fn64/TeExp/library/
EXE_DIR          := /home/fas/gerstein/fn64/tools/manual/

BOWTIE_BIN       := $(EXE_DIR)/bin/bowtie2
BOWTIE_PARAMS    := --sensitive-local -N1 --no-unal
BOWTIE_INDEX     := "/scratch/fas/gerstein/fn64/genome/Homo_sapiens/hg38/toBowtie2/hg38"

KALLISTO_BIN     := $(EXE_DIR)/bin/kallisto
KALLISTO_INDEX   := /scratch/fas/gerstein/fn64/genome/Homo_sapiens/gencode23/toKallisto/gencode.v23.transcripts.index

SAMTOOLS_BIN     := $(EXE_DIR)/bin/samtools
INTERSERC_BIN    := $(EXE_DIR)/bin/intersectBed
FASTX_FILTER_EXE := $(EXE_DIR)/bin/fastq_quality_filter
WGSIM_BIN        := $(EXE_DIR)/bin/wgsim
R_BIN            := $(EXE_DIR)/bin/R
PYTHON_BIN		 := /usr/bin/python

REPEAT_MASKER_BED          := /scratch/fas/gerstein/fn64/TeExp/hg38.rep.noexon.bed
REPEAT_MASKER_TOT_BED      := /scratch/fas/gerstein/fn64/TeExp/hg38.rep.bed

QFILTER_MIN_READ_FRAC      := 80
QFILTER_MIN_QUAL           := 20

COMMAND_HOMOPOL            := perl $(LIBRARY_PATH)/scripts/remove_homopol.pl
COMMAND_PARTIAL            := perl $(LIBRARY_PATH)/scripts/filter_qual.pl

##
## Simulation parameters
##
ERROR_RATE 		:= 0.1
NUMBER_OF_READS := 250000
NUMBER_OF_LOOPS := 100

