
PIPELINE_NAME = TeXP

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


EXE_DIR          := /home2/fn64/tools/manual

BOWTIE_BIN       := $(EXE_DIR)/bin/bowtie2
BOWTIE_PARAMS    := --sensitive-local -N1 --no-unal
SAMTOOLS_BIN     := $(EXE_DIR)/bin/samtools
INTERSERC_BIN    := $(EXE_DIR)/bin/intersectBed
BOWTIE_INDEX     := "/home2/fn64/genomes/Homo_sapiens/hg38/toBowtie2/hg38"
FASTX_FILTER_EXE := $(EXE_DIR)/bin/fastq_quality_filter
WGSIM_BIN        := /home2/fn64/tools/manual/bin/wgsim

LIBRARY_PATH     := /home2/fn64/projects/TeXP/library

REPEAT_MASKER_OUT     := 
REPEAT_MASKER_BED     := ~/projects/hg38.rep.noexon.bed
REPEAT_MASKER_TOT_BED := ~/projects/hg38.rep.bed

MEAN_READ_LEN                   := 100
QFILTER_MIN_READ_FRAC           := 80
QFILTER_MIN_QUAL                := 20

LOG_FILE         := $(OUTPUT_DIR)/$(SAMPLE_ID).log

##
SAMPLE_ID := $(INPUT_FILE_ID)
ifneq ($(SAMPLE_NAME),NULL)
  SAMPLE_ID := $(SAMPLE_ID)_$(SAMPLE_NAME)
endif

##
## Simulation parameters
##
ERROR_RATE 		:= 0.1
NUMBER_OF_READS := 25000
NUMBER_OF_LOOPS := 100

##
## Detect filetype and extract from SRA format if necessary
##
COMMAND_CONVERT_INPUT := cat $(INPUT_FILE_PATH) > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).fastq
ifeq ($(suffix $(INPUT_FILE_PATH)),.sra)
	COMMAND_CONVERT_INPUT := $(SRATOOLS_EXE) --split-spot --stdout $(INPUT_FILE_PATH) > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).fastq
else ifeq ($(suffix $(INPUT_FILE_PATH)),.gz)
	COMMAND_CONVERT_INPUT := gunzip -c $(INPUT_FILE_PATH) > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).fastq
endif


##
## Bowtie2 command to align reads to a Reference genome
##
COMMAND_MAP := $(BOWTIE_BIN) -p $(N_THREADS) $(BOWTIE_PARAMS) -x $(BOWTIE_INDEX) -U $(INPUT_FILE_PATH) 2>> $(LOG_FILE) | $(SAMTOOLS_BIN) view -Sb - 2>> $(LOG_FILE) > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).bam; $(SAMTOOLS_BIN) sort -@$(N_THREADS) -m 4G $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).bam $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).sorted; rm -R $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).bam; 
COMMAND_HOMOPOL := perl ~/projects/remove_homopol.pl
COMMAND_PARTIAL := perl ~/projects/filter_qual.pl

USAGE := 
ifeq ($(INPUT_FILE_ID),NULL)
  USAGE := "make -f $PIPELINE_NAME 
#  		INPUT_FILE_PATH=[required: absolute/path/to/input/.fa|.fq|.sra|.fa.gz] 
  		INPUT_FILE_PATH=[required: absolute/path/to/input/.fastq] 
  		N_THREADS=[required: number of threads] 
  		OUTPUT_DIR=<required: absolute/path/to/output> 
  		INPUT_FILE_ID=[required: samplename] ADAPTER_SEQ=[optional: will guess sequence if not provided here; none, if already clipped input] 
  		MAIN_ORGANISM=[optional: defaults to 'hsa'] 
  		MAIN_ORGANISM_GENOME_ID=[optional: defaults to 'hg38'] 
endif


## Define current time
timestamp := `/bin/date "+%Y-%m-%d(%H:%M:%S)"`


##
## Main make target
##
.PHONY: all
.DEFAULT: all
all: processSample

##
## TO-DO:
##  - ReadLength: Make it independent of the conversion process.
##  - Add a lock to the 
##

##
## Make results directory & Convert INPUT if necessary
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).fastq: $(INPUT_FILE_PATH)
	@echo -e "$(USAGE)"
	mkdir -p $(OUTPUT_DIR)/$(SAMPLE_ID)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Created results dir: $(OUTPUT_DIR)/$(SAMPLE_ID)\n" >> $(LOG_FILE)
	$(COMMAND_CONVERT_INPUT) 
	@echo -e "$(timestamp) $(PIPELINE_NAME): Converted INPUT_FILE: $(INPUT_FILE_NAME) to $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).fastq\n" >> $(LOG_FILE)

##
## Guess meanread length
##
update_read_lenth: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).fastq
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Guessing read legth based on fastq sequences:\n" >> $(LOG_FILE)
	$(eval MEAN_READ_LEN := $(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).fastq | head -n 40000 | awk '{if((NR+2)%4==0) {count++; sum+=length($$_)}} END{print sum/count}'))
	@echo -e "$(timestamp) $(PIPELINE_NAME): Finished guessing read legth based on fastq sequences:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log

##
## Guess Fastq quality encoding
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).qualityEncoding: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).fastq update_read_lenth $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Guessing encoding of fastq read-qualities:\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).fastq | head -n 40000 | awk '{if(NR%4==0) printf("%s",$$0);}' | od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($$i>max) max=$$i; if($$i<min) min=$$i;}}END{if(max<=74 && min<59) print "33"; else if(max>73 && min>=64) print "64"; else if(min>=59 && min<64 && max>73) print "64"; else print "64";}' > $@\n" >> $(LOG_FILE)
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).fastq | head -n 40000 | awk '{if(NR%4==0) printf("%s",$$0);}' | od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($$i>max) max=$$i; if($$i<min) min=$$i;}}END{if(max<=74 && min<59) print "33"; else if(max>73 && min>=64) print "64"; else if(min>=59 && min<64 && max>73) print "64"; else print "64";}' > $@
#	$(eval MEAN_READ_LEN := $(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).fastq | head -n 40000 | awk '{if((NR+2)%4==0) {count++; sum+=length($$_)}} END{print sum/count}'))
	@echo -e "$(timestamp) $(PIPELINE_NAME): Finished guessing encoding of fastq read-qualities:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log


##
## FILTER clipped reads that have poor overall base quality  &  Remove homopolymer repeats
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).filtered.fastq: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).qualityEncoding 
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Filtering reads by base quality:\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).fastq | $(FASTX_FILTER_EXE) -v -Q$(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).qualityEncoding) -p $(QFILTER_MIN_READ_FRAC) -q $(QFILTER_MIN_QUAL) > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).filtered.fastq\n" >> $(LOG_FILE) 2>> $(OUTPUT_DIR)/$(SAMPLE_ID).err
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).fastq | $(FASTX_FILTER_EXE) -v -Q$(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).qualityEncoding) -p $(QFILTER_MIN_READ_FRAC) -q $(QFILTER_MIN_QUAL) > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).filtered.fastq 2>>$(LOG_FILE)
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
	$(INTERSERC_BIN) -f 0.75 -a $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).sorted.bam -b $(REPEAT_MASKER_BED) -sorted > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.bam
	@echo -e "$(timestamp) $(PIPELINE_NAME): Removing Homopolymers and Partially mapped reads:\n" >> $(LOG_FILE)
	$(SAMTOOLS_BIN) view -h $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.bam | $(COMMAND_HOMOPOL) | $(COMMAND_PARTIAL) | $(SAMTOOLS_BIN) view -b - > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bam
	@echo -e "$(timestamp) $(PIPELINE_NAME): Creating BED with reads coordinates:\n" >> $(LOG_FILE)
	$(INTERSERC_BIN) -f 0.75 -a $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bam -b $(REPEAT_MASKER_TOT_BED) -sorted -bed -wo > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bed


##
## SIMULATE reads from L1 if they do not exists
##
$(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt: 
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): The profile for this study was not found at: $(LIBRARY_PATH)/L1/$(MEAN_READ_LEN).counts.txt\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Simulating reads with length equal to $(MEAN_READ_LEN)\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Creating reads from based on L1 reference sequence:\n" >> $(LOG_FILE)
	mkdir -p /home2/fn64/projects/TeXP/library/L1/simu/
	@for iter in $(shell seq 1 $(NUMBER_OF_LOOPS) ); do \
#		$(WGSIM_BIN) -1 $(MEAN_READ_LEN) -N $(NUMBER_OF_READS) -d0 -r $(ERROR_RATE) -e 0 -R 0 $(LIBRARY_PATH)/L1/ref/L1HS.ref.fa $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.simu /dev/null 2> /dev/null > /dev/null ; \
		$(WGSIM_BIN) -S $$(date "+%N") -1 $(MEAN_READ_LEN) -N $(NUMBER_OF_READS) -d0 -r$(ERROR_RATE) -e 0 -R 0 $(LIBRARY_PATH)/L1/ref/L1HS.ref.fa $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.simu /dev/null > /dev/null 2> /dev/null ; \
#		cat $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.simu | sed -n '1~4s/^@/>/p;2~4p' > $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.fasta; rm $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.simu; \
    done
#	@echo -e "$(timestamp) $(PIPELINE_NAME): Aligning simulated reads to the reference genome and counting L1 subfamilies:\n" >> $(LOG_FILE)
	@for iter in $(shell seq 1 $(NUMBER_OF_LOOPS) ); do \
		$(BOWTIE_BIN) -p $(N_THREADS) $(BOWTIE_PARAMS) -x $(BOWTIE_INDEX) -U $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.simu 2>> $(LOG_FILE) | $(SAMTOOLS_BIN) view -Sb - 2>> $(LOG_FILE) > $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.bam; \
		$(SAMTOOLS_BIN) sort -@$(N_THREADS) -m 4G $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.bam $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted; \
		rm -R $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.bam; \
		$(INTERSERC_BIN) -f 0.75 -a $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam -b $(LIBRARY_PATH)/L1/ref/L1.hg38.bed -sorted -bed -wo > $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.L1.bed; \
		cat $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.L1.bed | awk '{print $$(NF-1)}' | sort | uniq -c > $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.L1.bed.count; \
	done
	@echo -e "$(timestamp) $(PIPELINE_NAME): Calculating the expected number of reads on each subfamily:\n" >> $(LOG_FILE)
	cat $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_*.sorted.bam.L1.bed.count | sort -k2,2 | sed 's/^[ ]*//g' | awk 'BEGIN{first=1} {if ( first == 1 ) {id=$$2;first=0}; if ( id != $$2 ) {print sum/count,id;id=$$2;sum=0;count=0}; sum+=$$1;count++}; END{print sum/count,id;}' > $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt

##
## L1 Quantification repetitive element reads
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bed
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Counting the number of reads on each L1 subfamily:\n" >> $(LOG_FILE)
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bed | egrep -w "L1P1_LINE__L1|L1PA2_LINE__L1|L1PA3_LINE__L1|L1PA4_LINE__L1|L1HS_LINE__L1" | awk '{print $$(NF-1)}' | sort | uniq -c > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count

##
## Correcting the number of reads mapped to L1Hs
##	
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count.corrected: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Counting the number of reads on each L1 subfamily:\n" >> $(LOG_FILE)

##
## Main sub-target
##
processSample: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count.corrected
#	## Copy Output descriptions file
#	#cp $(SRNABENCH_LIBS)/sRNAbenchOutputDescription.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/sRNAbenchOutputDescription.txt 
#	## END PIPELINE
#	#@echo -e "$(ts) SMRNAPIPELINE: END smallRNA-seq Pipeline for sample $(SAMPLE_ID)\n======================\n" >> $(LOG_FILE)


