
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
R_BIN            := /home2/fn64/tools/manual/bin/R

LIBRARY_PATH     := /home2/fn64/projects/TeXP/library

REPEAT_MASKER_OUT     := 
REPEAT_MASKER_BED     := ~/projects/hg38.rep.noexon.bed
REPEAT_MASKER_TOT_BED := ~/projects/hg38.rep.bed

MEAN_READ_LEN                   := 100
QFILTER_MIN_READ_FRAC           := 80
QFILTER_MIN_QUAL                := 20

##
SAMPLE_ID := $(INPUT_FILE_ID)
ifneq ($(SAMPLE_NAME),NULL)
  SAMPLE_ID := $(SAMPLE_NAME)
endif

LOG_FILE         := $(OUTPUT_DIR)/$(SAMPLE_ID).log

##
## Simulation parameters
##
ERROR_RATE 		:= 0.1
NUMBER_OF_READS := 25000
NUMBER_OF_LOOPS := 100

##
## Detect filetype and extract from SRA format if necessary
##
COMMAND_CONVERT_INPUT := cat $(INPUT_FILE_PATH)
ifeq ($(suffix $(INPUT_FILE_PATH)),.sra)
	COMMAND_CONVERT_INPUT := $(SRATOOLS_EXE) --split-spot --stdout $(INPUT_FILE_PATH) 
else ifeq ($(suffix $(INPUT_FILE_PATH)),.gz)
	COMMAND_CONVERT_INPUT := gunzip -c $(INPUT_FILE_PATH) 
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
.PHONY: all update_read_lenth
.DEFAULT: all
all: processSample

##
## TO-DO:
##  - Add a lock in the simulation to avoid multiple simulations with the same parameters to run in parallel 
##  - $(LIBRARY_PATH)/L1/ref/L1.bases.ref: Create this file dinamically based on a reference genome, instead of a static file?

##
## Make results directory & Convert INPUT if necessary
##
$(OUTPUT_DIR)/$(SAMPLE_ID): $(INPUT_FILE_PATH)
	@echo -e "$(USAGE)"
	mkdir -p $(OUTPUT_DIR)/$(SAMPLE_ID)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Created results dir: $(OUTPUT_DIR)/$(SAMPLE_ID)\n" >> $(LOG_FILE)
#	$(COMMAND_CONVERT_INPUT) 
#	@echo -e "$(timestamp) $(PIPELINE_NAME): Converted INPUT_FILE: $(INPUT_FILE_NAME) to $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).fastq\n" >> $(LOG_FILE)

##
## Guess meanread length
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).read_lenth: $(INPUT_FILE_PATH) $(OUTPUT_DIR)/$(SAMPLE_ID)
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Guessing read legth based on fastq sequences:\n" >> $(LOG_FILE)
	$(eval MEAN_READ_LEN := $(shell $(COMMAND_CONVERT_INPUT) | head -n 40000 | awk '{if((NR+2)%4==0) {count++; sum+=length($$_)}} END{print sum/count}'))
	echo "$(MEAN_READ_LEN)" > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).read_lenth
	@echo -e "$(timestamp) $(PIPELINE_NAME): Finished guessing read legth based on fastq sequences:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log

##
## Guess Fastq quality encoding
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).qualityEncoding: $(INPUT_FILE_PATH) $(OUTPUT_DIR)/$(SAMPLE_ID) $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).read_lenth $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Guessing encoding of fastq read-qualities:\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): $(COMMAND_CONVERT_INPUT) | head -n 40000 | awk '{if(NR%4==0) printf("%s",$$0);}' | od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($$i>max) max=$$i; if($$i<min) min=$$i;}}END{if(max<=74 && min<59) print "33"; else if(max>73 && min>=64) print "64"; else if(min>=59 && min<64 && max>73) print "64"; else print "64";}' > $@\n" >> $(LOG_FILE)
#	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).fastq | head -n 40000 | awk '{if(NR%4==0) printf("%s",$$0);}' | od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($$i>max) max=$$i; if($$i<min) min=$$i;}}END{if(max<=74 && min<59) print "33"; else if(max>73 && min>=64) print "64"; else if(min>=59 && min<64 && max>73) print "64"; else print "64";}' > $@
	$(COMMAND_CONVERT_INPUT) | head -n 40000 | awk '{if(NR%4==0) printf("%s",$$0);}' | od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($$i>max) max=$$i; if($$i<min) min=$$i;}}END{if(max<=74 && min<59) print "33"; else if(max>73 && min>=64) print "64"; else if(min>=59 && min<64 && max>73) print "64"; else print "64";}' > $@
	@echo -e "$(timestamp) $(PIPELINE_NAME): Finished guessing encoding of fastq read-qualities:\n" >> $(OUTPUT_DIR)/$(SAMPLE_ID).log


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
		$(WGSIM_BIN) -S $$(date "+%N") -1 $(MEAN_READ_LEN) -N $(NUMBER_OF_READS) -d0 -r$(ERROR_RATE) -e 0 -R 0 $(LIBRARY_PATH)/L1/ref/L1HS.ref.fa $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.simu /dev/null > /dev/null 2> /dev/null ; \
    done
	@for iter in $(shell seq 1 $(NUMBER_OF_LOOPS) ); do \
		$(BOWTIE_BIN) -p $(N_THREADS) $(BOWTIE_PARAMS) -x $(BOWTIE_INDEX) -U $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.simu 2>> $(LOG_FILE) | $(SAMTOOLS_BIN) view -Sb - 2>> $(LOG_FILE) > $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.bam; \
		$(SAMTOOLS_BIN) sort -@$(N_THREADS) -m 4G $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.bam $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted; \
		rm -R $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.bam; \
		$(INTERSERC_BIN) -f 0.75 -a $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam -b $(LIBRARY_PATH)/L1/ref/L1.hg38.bed -sorted -bed -wo > $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.L1.bed; \
		cat $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.L1.bed | awk '{print $$(NF-1)}' | sort | uniq -c > $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.L1.bed.count; \
	done
	@echo -e "$(timestamp) $(PIPELINE_NAME): Calculating the expected number of reads on each subfamily:\n" >> $(LOG_FILE)
	echo "L1_Subfamily L1_Ref_bases" > $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt
	cat $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_*.sorted.bam.L1.bed.count | sort -k2,2 | sed 's/^[ ]*//g' | awk 'BEGIN{first=1} {if ( first == 1 ) {id=$$2;first=0}; if ( id != $$2 ) {print sum/count,id;id=$$2;sum=0;count=0}; sum+=$$1;count++}; END{print id,sum/count;}' | sed 's/_LINE__L1//g' >> $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt

##
## Create auxiliary file with proportion of simulated reads on each subfamily
##
$(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).prop.txt: $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt
	$(eval T_SUM := $(shell cat $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt | awk '{sum+=$$2} END{print sum}'))
	cat $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt | awk '{sum+=$$2} END{print sum}'
	cat $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt | awk -v sum=$(T_SUM) '{if ($$2 ~ /[0-9]*[.][0-9]*/ ) {print $$1,$$2/sum} else {print $$1,$$2} }' > $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).prop.txt

##
## Create signature file
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/L1.signatures.txt: $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).prop.txt $(LIBRARY_PATH)/L1/ref/L1.bases.ref
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Compiling L1 signature files:\n" >> $(LOG_FILE)
	paste $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).prop.txt $(LIBRARY_PATH)/L1/ref/L1.bases.ref | awk '{print $$1,$$2,$$4}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/L1.signatures.txt

##
## Quantification of L1 repetitive element reads
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bed
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Counting the number of reads on each L1 subfamily:\n" >> $(LOG_FILE)
	echo "L1_count L1_Subfamily" > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bed | egrep -w "L1P1_LINE__L1|L1PA2_LINE__L1|L1PA3_LINE__L1|L1PA4_LINE__L1|L1HS_LINE__L1" | awk '{print $$(NF-1)}' | sort | uniq -c | sed 's/_LINE__L1//g' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count | awk '{print $$2,$$1}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count.t
	mv $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count.t $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count 

##
## Correcting the number of reads mapped to L1Hs
##	
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count.corrected: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count $(OUTPUT_DIR)/$(SAMPLE_ID)/L1.signatures.txt
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Correcting the number of reads on L1Hs:\n" >> $(LOG_FILE)
	$(R_BIN) --no-restore --no-save --args $(OUTPUT_DIR)/$(SAMPLE_ID)/L1.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).sorted.bam.tot < $(LIBRARY_PATH)/L1/ref/lsei.template.r >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Writing L1 quantification files:" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME):  - $(OUTPUT_DIR)/$(SAMPLE_ID)/L1.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count.corrected" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME):  - $(OUTPUT_DIR)/$(SAMPLE_ID)/L1.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count.rpkm" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME):  - $(OUTPUT_DIR)/$(SAMPLE_ID)/L1.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count.rpkm.corrected" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME):  - $(OUTPUT_DIR)/$(SAMPLE_ID)/L1.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count.signal_proportions" >> $(LOG_FILE)


##
## Main sub-target
##
processSample: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count.corrected
#	## Copy Output descriptions file
#	#cp $(SRNABENCH_LIBS)/sRNAbenchOutputDescription.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/sRNAbenchOutputDescription.txt 
#	## END PIPELINE
#	#@echo -e "$(ts) SMRNAPIPELINE: END smallRNA-seq Pipeline for sample $(SAMPLE_ID)\n======================\n" >> $(LOG_FILE)


