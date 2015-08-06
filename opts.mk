
PIPELINE_NAME     := TeXP

DATA_DIR          := NULL
OUTPUT_DIR        := NULL
INPUT_FILE_PATH   := NULL
SAMPLE_NAME       := NULL
REFERENCE_GENOME  := NULL

##
## Use the input path to infer filetype and short name
##
INPUT_FILE_NAME := $(notdir $(INPUT_FILE_PATH))
INPUT_FILE_ID   := $(basename $(INPUT_FILE_NAME))

LIBRARY_PATH     := /home2/fn64/projects/TeXP/library
EXE_DIR          := /home2/fn64/tools/manual

BOWTIE_BIN       := $(EXE_DIR)/bin/bowtie2
BOWTIE_PARAMS    := --sensitive-local -N1 --no-unal
BOWTIE_INDEX     := "/home2/fn64/genomes/Homo_sapiens/hg38/toBowtie2/hg38"

SAMTOOLS_BIN     := $(EXE_DIR)/bin/samtools
INTERSERC_BIN    := $(EXE_DIR)/bin/intersectBed
FASTX_FILTER_EXE := $(EXE_DIR)/bin/fastq_quality_filter
WGSIM_BIN        := $(EXE_DIR)/bin/wgsim
R_BIN            := $(EXE_DIR)/bin/R
KALLISTO_BIN     := $(EXE_DIR)/bin/kallisto

REPEAT_MASKER_OUT     := 
REPEAT_MASKER_BED     := ~/projects/hg38.rep.noexon.bed
REPEAT_MASKER_TOT_BED := ~/projects/hg38.rep.bed
KALLISTO_INDEX        := ~/projects/genomes/Homo_sapiens/gencode23/toKallisto/gencode.v23.transcripts.index

QFILTER_MIN_READ_FRAC           := 80
QFILTER_MIN_QUAL                := 20

COMMAND_HOMOPOL := perl $(LIBRARY_PATH)/scripts/remove_homopol.pl
COMMAND_PARTIAL := perl $(LIBRARY_PATH)/scripts/filter_qual.pl

##
## Simulation parameters
##
ERROR_RATE 		:= 0.1
NUMBER_OF_READS := 25000
NUMBER_OF_READS_SVA := 100000
NUMBER_OF_LOOPS := 100

LOG_FILE         := $(OUTPUT_DIR)/$(SAMPLE_ID).log

