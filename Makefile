include opts.mk

##
SAMPLE_ID := $(INPUT_FILE_ID)
ifneq ($(SAMPLE_NAME),NULL)
	SAMPLE_ID := $(SAMPLE_NAME)
endif

##
## Detect filetype and extract from SRA format if necessary
##
COMMAND_CONVERT_INPUT := cat $(INPUT_FILE_PATH)
ifeq ($(suffix $(INPUT_FILE_PATH)),.sra)
	COMMAND_CONVERT_INPUT := $(SRATOOLS_EXE) --split-spot --stdout $(INPUT_FILE_PATH) 
else ifeq ($(suffix $(INPUT_FILE_PATH)),.gz)
	COMMAND_CONVERT_INPUT := gunzip -c $(INPUT_FILE_PATH) 
else ifeq ($(suffix $(INPUT_FILE_PATH)),.tar.gz)
	COMMAND_CONVERT_INPUT := tar -xOzf $(INPUT_FILE_PATH) 
endif

USAGE := 
ifeq ($(INPUT_FILE_ID),NULL)
	USAGE := "make -f $PIPELINE_NAME 
		INPUT_FILE_PATH=[required: absolute/path/to/input/.fa|.fq|.sra|.fa.gz] 
		N_THREADS=[required: number of threads] 
		OUTPUT_DIR=<required: absolute/path/to/output> 
		INPUT_FILE_ID=[required: samplename] ADAPTER_SEQ=[optional: will guess sequence if not provided here; none, if already clipped input] 
		MAIN_ORGANISM=[optional: defaults to 'hsa'] 
		MAIN_ORGANISM_GENOME_ID=[optional: defaults to 'hg38']"
endif

ifneq (,$(wildcard $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).read_length))
OUT_LEN="$(timestamp) $(PIPELINE_NAME): Found $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).read_length \n" >> $(LOG_FILE)
$(eval MEAN_READ_LEN := $(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).read_length))
else
OUT_LEN="$(timestamp) $(PIPELINE_NAME): Did not found $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).read_length \n"
$(eval MEAN_READ_LEN := $(shell $(COMMAND_CONVERT_INPUT) | head -n 1000000 | awk '{if ((NR+2)%4==0) {i++; count[i] = length($1);}} END {asort(count); print count[int(i/2)];}'))
endif
#echo $(MEAN_READ_LEN)
$(info $(OUT_LEN))

LOG_FILE := $(OUTPUT_DIR)/$(SAMPLE_ID).log

##
## Bowtie2 command to align reads to a Reference genome
##
COMMAND_MAP := $(BOWTIE_BIN) -p $(N_THREADS) $(BOWTIE_PARAMS) -x $(BOWTIE_INDEX) -U $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).filtered.fastq 2>> $(LOG_FILE) | $(SAMTOOLS_BIN) view -Sb - 2>> $(LOG_FILE) > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).bam; $(SAMTOOLS_BIN) sort -@$(N_THREADS) $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).bam -o $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).sorted.bam; rm -R $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).bam; 



## Define current time
timestamp := `/bin/date "+%Y-%m-%d(%H:%M:%S)"`

##
## Main make target
##
.PHONY: all lock_L1 lock_SVA lock_LTR
.DEFAULT: all
.INTERMEDIATE: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).filtered.fastq $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).sorted.bam 

all: processSample 

##
## TO-DO:
##  - $(LIBRARY_PATH)/L1/ref/L1.bases.ref: Create this file dinamically based on a reference genome, instead of a static file?
##

##
## Make results directory & Convert INPUT if necessary
##
$(OUTPUT_DIR)/$(SAMPLE_ID): | $(INPUT_FILE_PATH)
	@echo -e "$(USAGE)"
	mkdir -p $(OUTPUT_DIR)/$(SAMPLE_ID)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Created results dir: $(OUTPUT_DIR)/$(SAMPLE_ID)\n" >> $(LOG_FILE)

##
## Guess meanread length
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).read_length: | $(INPUT_FILE_PATH) $(OUTPUT_DIR)/$(SAMPLE_ID)
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Guessing read legth based on fastq sequences:\n" >> $(LOG_FILE)
	$(eval MEAN_READ_LEN := $(shell $(COMMAND_CONVERT_INPUT) | head -n 1000000 | awk '{if ((NR+2)%4==0) {i++; count[i] = length($1);}} END {asort(count); print count[int(i/2)];}'))
	echo "$(MEAN_READ_LEN)" > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).read_length
	export MEAN_READ_LEN
	@echo -e "$(timestamp) $(PIPELINE_NAME): Finished guessing read legth based on fastq sequences:\n" >> $(LOG_FILE)

##
## Guess Fastq quality encoding
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).qualityEncoding: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).read_length
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Guessing encoding of fastq read-qualities:\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): $(COMMAND_CONVERT_INPUT) | head -n 400000 | awk '{if(NR%4==0) printf("%s",$$0);}' | od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($$i>max) max=$$i; if($$i<min) min=$$i;}}END{if(max<=74 && min<59) print "33"; else if(max>73 && min>=64) print "64"; else if(min>=59 && min<64 && max>73) print "64"; else print "64";}' > $@\n" >> $(LOG_FILE)
	$(COMMAND_CONVERT_INPUT) | head -n 400000 | awk '{if(NR%4==0) printf("%s",$$0);}' | od -A n -t u1 | grep -v "^*" | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($$i>max) max=$$i; if($$i<min) min=$$i;}}END{if(max<=74 && min<59) print "33"; else if(max>73 && min>=64) print "64"; else if(min>=59 && min<64 && max>73) print "64"; else print "Unknown";}' > $@
	@echo -e "$(timestamp) $(PIPELINE_NAME): Finished guessing encoding of fastq read-qualities:\n" >> $(LOG_FILE)


##
## FILTER clipped reads that have poor overall base quality  &  Remove homopolymer repeats
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).filtered.fastq: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).qualityEncoding 
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Filtering reads by base quality:\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): $(COMMAND_CONVERT_INPUT) | $(FASTX_FILTER_EXE) -v -Q$(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).qualityEncoding) -p $(QFILTER_MIN_READ_FRAC) -q $(QFILTER_MIN_QUAL) > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).filtered.fastq\n" >> $(LOG_FILE) 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	$(COMMAND_CONVERT_INPUT) | $(FASTX_FILTER_EXE) -v -Q$(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).qualityEncoding) -p $(QFILTER_MIN_READ_FRAC) -q $(QFILTER_MIN_QUAL) > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).filtered.fastq 2>>$(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Finished filtering reads by base quality\n" >> $(LOG_FILE)
	## Count reads that failed the quality filter
	grep "low-quality reads" $(LOG_FILE) | awk '{print "failed_quality_filter\t"$$2}' >> $(OUTPUT_DIR)/$(SAMPLE_ID).stats


##
## MAP filtered reads to a reference genome
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).sorted.bam: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).filtered.fastq
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Mapping reads to a reference genome:\n" >> $(LOG_FILE)
	echo -e "$(timestamp) $(PIPELINE_NAME): $(COMMAND_MAP)" >> $(LOG_FILE)
	$(COMMAND_MAP) 
	@echo -e "$(timestamp) $(PIPELINE_NAME): Indexing bam file:\n" >> $(LOG_FILE)
	$(SAMTOOLS_BIN) index $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).sorted.bam
	@echo -e "$(timestamp) $(PIPELINE_NAME): Counting total number of mapped reads:\n" >> $(LOG_FILE)
	samtools idxstats $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).sorted.bam | awk '{sum+=$$3} END{print sum}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).sorted.bam.tot


##
## FILTER repetitive element reads
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bed: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).sorted.bam
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Intersecting reads with repeat masked regions:\n" >> $(LOG_FILE)
	$(INTERSERC_BIN) -f 0.75 -a $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).sorted.bam -b $(REPEAT_MASKER_BED) > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.bam
	@echo -e "$(timestamp) $(PIPELINE_NAME): Removing Homopolymers and Partially mapped reads:\n" >> $(LOG_FILE)
	$(SAMTOOLS_BIN) view -h $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.bam | $(COMMAND_HOMOPOL) | $(COMMAND_PARTIAL) | $(SAMTOOLS_BIN) view -b - > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bam
	@echo -e "$(timestamp) $(PIPELINE_NAME): Creating BED with reads coordinates:\n" >> $(LOG_FILE)
	$(INTERSERC_BIN) -f 0.75 -a $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bam -b $(REPEAT_MASKER_TOT_BED) -bed -wo > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bed

##
## Quantifying transcripts annotated by GENCODE usign Kallisto 
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/transcripts_quantification/abundance.txt: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).filtered.fastq
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Fast quantify transcripts with Kallisto:\n" >> $(LOG_FILE)
	$(KALLISTO_BIN) quant -i $(KALLISTO_INDEX) -o $(OUTPUT_DIR)/$(SAMPLE_ID)/transcripts_quantification --single -l $(MEAN_READ_LEN) --plaintext $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).filtered.fastq

##
## Transcript count factor -- for TPM quantification
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).tpm.factor: $(OUTPUT_DIR)/$(SAMPLE_ID)/transcripts_quantification/abundance.txt
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Create fator file:\n" >> $(LOG_FILE)
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/transcripts_quantification/abundance.txt | grep -v "eff_length" | awk '{a+=$$4/$$3} END{print a}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).tpm.factor


include library/L1HS_hg38/ref/L1HS_hg38.makefile.sub
#include library/HERV_hg38/ref/HERV_hg38.makefile.sub
#include library/LINE2_hg38/ref/LINE2_hg38.makefile.sub
#include library/LTR_hg38/ref/LTR_hg38.makefile.sub
#include library/SVA_hg38/ref/SVA_hg38.makefile.sub

##
## Main sub-target
##
#!!!!!-----ALL REs:-----!!!!! processSample: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1HS_hg38.count.corrected $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).HERV_hg38.count.corrected $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LINE2_hg38.count.corrected $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR_hg38.count.corrected $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).SVA_hg38.count.corrected 
#processSample: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1HS_hg38.count.corrected $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).HERV_hg38.count.corrected $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LINE2_hg38.count.corrected $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).SVA_hg38.count.corrected $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR_hg38.count.corrected
processSample: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1HS_hg38.count.corrected
#	## Copy Output descriptions file
#	#cp $(SRNABENCH_LIBS)/sRNAbenchOutputDescription.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/sRNAbenchOutputDescription.txt 
#	## END PIPELINE
#	#@echo -e "$(ts) SMRNAPIPELINE: END smallRNA-seq Pipeline for sample $(SAMPLE_ID)\n======================\n" >> $(LOG_FILE)


