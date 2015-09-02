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
endif

USAGE := 
ifeq ($(INPUT_FILE_ID),NULL)
  USAGE := "make -f $PIPELINE_NAME 
  		INPUT_FILE_PATH=[required: absolute/path/to/input/.fa|.fq|.sra|.fa.gz] 
  		N_THREADS=[required: number of threads] 
  		OUTPUT_DIR=<required: absolute/path/to/output> 
  		INPUT_FILE_ID=[required: samplename] ADAPTER_SEQ=[optional: will guess sequence if not provided here; none, if already clipped input] 
  		MAIN_ORGANISM=[optional: defaults to 'hsa'] 
  		MAIN_ORGANISM_GENOME_ID=[optional: defaults to 'hg38'] 
endif

ifneq ("$(wildcard $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).read_length)","")
	@echo -e "$(timestamp) $(PIPELINE_NAME): Did not found $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).read_length \n" >> $(LOG_FILE)
$(eval MEAN_READ_LEN := $(shell $(COMMAND_CONVERT_INPUT) | head -n 40000 | awk '{if((NR+2)%4==0) {count++; sum+=length($$_)}} END{print sum/count}'))
else
	@echo -e "$(timestamp) $(PIPELINE_NAME): Found $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).read_length \n" >> $(LOG_FILE)
$(eval MEAN_READ_LEN := $(shell cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).read_length)
endif

LOG_FILE := $(OUTPUT_DIR)/$(SAMPLE_ID).log

##
## Bowtie2 command to align reads to a Reference genome
##
COMMAND_MAP := $(BOWTIE_BIN) -p $(N_THREADS) $(BOWTIE_PARAMS) -x $(BOWTIE_INDEX) -U $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).filtered.fastq 2>> $(LOG_FILE) | $(SAMTOOLS_BIN) view -Sb - 2>> $(LOG_FILE) > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).bam; $(SAMTOOLS_BIN) sort -@$(N_THREADS) $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).bam $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).sorted; rm -R $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).bam; 



## Define current time
timestamp := `/bin/date "+%Y-%m-%d(%H:%M:%S)"`

##
## Main make target
##
.PHONY: all lock_L1 lock_SVA lock_LTR
.DEFAULT: all
.INTERMEDIATE: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).filtered.fastq $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).sorted.bam $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt.lock

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
	$(eval MEAN_READ_LEN := $(shell $(COMMAND_CONVERT_INPUT) | head -n 40000 | awk '{if((NR+2)%4==0) {count++; sum+=length($$_)}} END{print sum/count}'))
	echo "$(MEAN_READ_LEN)" > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).read_length
	export MEAN_READ_LEN
	@echo -e "$(timestamp) $(PIPELINE_NAME): Finished guessing read legth based on fastq sequences:\n" >> $(LOG_FILE)

##
## Guess Fastq quality encoding
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).qualityEncoding: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).read_length
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Guessing encoding of fastq read-qualities:\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): $(COMMAND_CONVERT_INPUT) | head -n 40000 | awk '{if(NR%4==0) printf("%s",$$0);}' | od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($$i>max) max=$$i; if($$i<min) min=$$i;}}END{if(max<=74 && min<59) print "33"; else if(max>73 && min>=64) print "64"; else if(min>=59 && min<64 && max>73) print "64"; else print "64";}' > $@\n" >> $(LOG_FILE)
	$(COMMAND_CONVERT_INPUT) | head -n 40000 | awk '{if(NR%4==0) printf("%s",$$0);}' | od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($$i>max) max=$$i; if($$i<min) min=$$i;}}END{if(max<=74 && min<59) print "33"; else if(max>73 && min>=64) print "64"; else if(min>=59 && min<64 && max>73) print "64"; else print "64";}' > $@
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
	$(INTERSERC_BIN) -f 0.75 -a $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).sorted.bam -b $(REPEAT_MASKER_BED) -sorted > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.bam
	@echo -e "$(timestamp) $(PIPELINE_NAME): Removing Homopolymers and Partially mapped reads:\n" >> $(LOG_FILE)
	$(SAMTOOLS_BIN) view -h $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.bam | $(COMMAND_HOMOPOL) | $(COMMAND_PARTIAL) | $(SAMTOOLS_BIN) view -b - > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bam
	@echo -e "$(timestamp) $(PIPELINE_NAME): Creating BED with reads coordinates:\n" >> $(LOG_FILE)
	$(INTERSERC_BIN) -f 0.75 -a $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bam -b $(REPEAT_MASKER_TOT_BED) -sorted -bed -wo > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bed

##
## Quantifying transcripts annotated by GENCODE usign Kallisto 
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/transcripts_quantification/abundance.txt: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).filtered.fastq
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Fast quantify transcripts with Kallisto:\n" >> $(LOG_FILE)
	$(KALLISTO_BIN) quant -i $(KALLISTO_INDEX) -o $(OUTPUT_DIR)/$(SAMPLE_ID)/transcripts_quantification --single -l 50 --plaintext $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).filtered.fastq

##
## Transcript count factor (for TPM quantification)
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).tpm.factor: $(OUTPUT_DIR)/$(SAMPLE_ID)/transcripts_quantification/abundance.txt
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Create fator file:\n" >> $(LOG_FILE)
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/transcripts_quantification/abundance.txt | grep -v "eff_length" | awk '{a+=$$4/$$3} END{print a}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).tpm.factor


##          ##   
##        ####   
##          ##   
##          ##   
##          ##   
##          ##   
########  ###### 


##
## SIMULATE reads from L1 if they do not exists
##
$(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt: $(LIBRARY_PATH)/L1/ref/L1HS.ref.fa
ifneq ("$(wildcard $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt.lock)","")
	@echo -e "$(timestamp) $(PIPELINE_NAME): There is another simulation running. Exiting without finishing."
	exit 1
endif

	touch $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt.lock

	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): The profile for this study was not found at: $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Simulating reads with length equal to $(MEAN_READ_LEN)\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Creating reads from based on L1 reference sequence:\n" >> $(LOG_FILE)
	@mkdir -p $(LIBRARY_PATH)/L1/simu/
	@for iter in $(shell seq 1 $(NUMBER_OF_LOOPS) ); do \
		$(WGSIM_BIN) -S $$(date "+%N") -1 $(MEAN_READ_LEN) -N $(NUMBER_OF_READS) -d0 -r$(ERROR_RATE) -e 0 -R 0 $(LIBRARY_PATH)/L1/ref/L1HS.ref.fa $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.simu /dev/null > /dev/null 2> /dev/null ; \
	done
	@echo -e "$(timestamp) $(PIPELINE_NAME): Aligning simulated reads to the reference genome:\n" >> $(LOG_FILE)
	@for iter in $(shell seq 1 $(NUMBER_OF_LOOPS) ); do \
		$(BOWTIE_BIN) -p $(N_THREADS) $(BOWTIE_PARAMS) -x $(BOWTIE_INDEX) -U $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.simu 2>> $(LOG_FILE) | $(SAMTOOLS_BIN) view -Sb - 2>> $(LOG_FILE) > $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.bam; \
		$(SAMTOOLS_BIN) sort -@$(N_THREADS) -m 4G $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.bam $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted; \
		rm -R $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.bam; \
		$(INTERSERC_BIN) -f 0.75 -a $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam -b $(LIBRARY_PATH)/L1/ref/L1.hg38.bed -sorted -bed -wo > $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.L1.bed; \
		cat $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.L1.bed | awk '{print $$(NF-1)}' | sort | uniq -c > $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.L1.bed.count; \
	done
	@echo -e "$(timestamp) $(PIPELINE_NAME): Calculating the expected number of reads on each subfamily:\n" >> $(LOG_FILE)
	echo "L1_Subfamily L1Hs_Transcript" > $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt
	cat $(LIBRARY_PATH)/L1/simu/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE)_*.sorted.bam.L1.bed.count | sort -k2,2 | sed 's/^[ ]*//g' | awk 'BEGIN{first=1} {if ( first == 1 ) {id=$$2;first=0}; if ( id != $$2 ) {print id,sum/count;id=$$2;sum=0;count=0}; sum+=$$1;count++}; END{print id,sum/count;}' | sed 's/_LINE__L1//g' >> $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt

	rm -Rf $(LIBRARY_PATH)/L1/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt.lock


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
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count.corrected: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).tpm.factor $(OUTPUT_DIR)/$(SAMPLE_ID)/L1.signatures.txt
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Correcting the number of reads on L1Hs:\n" >> $(LOG_FILE)
	$(R_BIN) --no-restore --no-save --args $(OUTPUT_DIR)/$(SAMPLE_ID)/L1.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).sorted.bam.tot $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).tpm.factor < $(LIBRARY_PATH)/L1/ref/lsei.template.r >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Writing L1 quantification files:" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): - $(OUTPUT_DIR)/$(SAMPLE_ID)/L1.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count.corrected" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): - $(OUTPUT_DIR)/$(SAMPLE_ID)/L1.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count.rpkm" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): - $(OUTPUT_DIR)/$(SAMPLE_ID)/L1.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count.rpkm.corrected" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): - $(OUTPUT_DIR)/$(SAMPLE_ID)/L1.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count.signal_proportions" >> $(LOG_FILE)


 ######  ##     ##    ###    
##    ## ##     ##   ## ##   
##       ##     ##  ##   ##  
 ######  ##     ## ##     ## 
      ##  ##   ##  ######### 
##    ##   ## ##   ##     ## 
 ######     ###    ##     ## 


##
## SIMULATE reads from SVA if they do not exists
##
$(LIBRARY_PATH)/SVA/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt: $(LIBRARY_PATH)/SVA/ref/SVA.ref.fa | lock_SVA
ifneq ("$(wildcard  $(LIBRARY_PATH)/SVA/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt.lock)","")
	@echo -e "$(timestamp) $(PIPELINE_NAME): There is another simulation running. Exiting without finishing."
	exit 1
endif

	touch $(LIBRARY_PATH)/SVA/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt.lock

	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): The profile for this study was not found at: $(LIBRARY_PATH)/SVA/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Simulating reads with length equal to $(MEAN_READ_LEN)\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Creating reads from based on SVA reference sequence:\n" >> $(LOG_FILE)

	mkdir -p $(LIBRARY_PATH)/SVA/simu/
	@for iter in $(shell seq 1 $(NUMBER_OF_LOOPS) ); do \
		$(WGSIM_BIN) -S $$(date "+%N") -1 $(MEAN_READ_LEN) -N $(NUMBER_OF_READS_SVA) -d0 -r$(ERROR_RATE) -e 0 -R 0 $(LIBRARY_PATH)/SVA/ref/SVA.ref.fa $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.simu /dev/null > /dev/null 2> /dev/null ; \
    done
	@echo -e "$(timestamp) $(PIPELINE_NAME): Aligning simulated reads to the reference genome:\n" >> $(LOG_FILE)
	@for iter in $(shell seq 1 $(NUMBER_OF_LOOPS) ); do \
		$(BOWTIE_BIN) -p $(N_THREADS) $(BOWTIE_PARAMS) -x $(BOWTIE_INDEX) -U $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.simu 2>> $(LOG_FILE) | $(SAMTOOLS_BIN) view -Sb - 2>> $(LOG_FILE) > $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.bam; \
		$(SAMTOOLS_BIN) sort -@$(N_THREADS) -m 4G $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.bam $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted; \
		rm -R $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.bam; \
		$(INTERSERC_BIN) -f 0.75 -a $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam -b $(LIBRARY_PATH)/SVA/ref/SVA.hg38.bed -sorted -bed -wo > $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.SVA.bed; \
		cat $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.SVA.bed | grep ")SVA_B_" | awk '{print $$(NF-1)}' | sort | uniq -c > $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.SVA_B.bed.count; \
		cat $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.SVA.bed | grep ")SVA_C_" | awk '{print $$(NF-1)}' | sort | uniq -c > $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.SVA_C.bed.count; \
		cat $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.SVA.bed | grep ")SVA_D_" | awk '{print $$(NF-1)}' | sort | uniq -c > $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.SVA_D.bed.count; \
		cat $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.SVA.bed | grep ")SVA_E_" | awk '{print $$(NF-1)}' | sort | uniq -c > $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.SVA_E.bed.count; \
		cat $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.SVA.bed | grep ")SVA_F_" | awk '{print $$(NF-1)}' | sort | uniq -c > $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.SVA_F.bed.count; \
	done
	@echo -e "$(timestamp) $(PIPELINE_NAME): Calculating the expected number of reads on each subfamily:\n" >> $(LOG_FILE)
	echo "SVA_Subfamily SVA_B_Transcript SVA_C_Transcript SVA_D_Transcript SVA_E_Transcript SVA_F_Transcript " > $(LIBRARY_PATH)/SVA/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt
	cat $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_*.sorted.bam.SVA_B.bed.count | sort -k2,2 | sed 's/^[ ]*//g' | awk 'BEGIN{first=1} {if ( first == 1 ) {id=$$2;first=0}; if ( id != $$2 ) {print id,sum/count;id=$$2;sum=0;count=0}; sum+=$$1;count++}; END{print id,sum/count;}' | sed 's/_Retroposon__SVA//g' > temp.SVA_B
	cat $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_*.sorted.bam.SVA_C.bed.count | sort -k2,2 | sed 's/^[ ]*//g' | awk 'BEGIN{first=1} {if ( first == 1 ) {id=$$2;first=0}; if ( id != $$2 ) {print id,sum/count;id=$$2;sum=0;count=0}; sum+=$$1;count++}; END{print id,sum/count;}' | sed 's/_Retroposon__SVA//g' > temp.SVA_C
	cat $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_*.sorted.bam.SVA_D.bed.count | sort -k2,2 | sed 's/^[ ]*//g' | awk 'BEGIN{first=1} {if ( first == 1 ) {id=$$2;first=0}; if ( id != $$2 ) {print id,sum/count;id=$$2;sum=0;count=0}; sum+=$$1;count++}; END{print id,sum/count;}' | sed 's/_Retroposon__SVA//g' > temp.SVA_D
	cat $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_*.sorted.bam.SVA_E.bed.count | sort -k2,2 | sed 's/^[ ]*//g' | awk 'BEGIN{first=1} {if ( first == 1 ) {id=$$2;first=0}; if ( id != $$2 ) {print id,sum/count;id=$$2;sum=0;count=0}; sum+=$$1;count++}; END{print id,sum/count;}' | sed 's/_Retroposon__SVA//g' > temp.SVA_E
	cat $(LIBRARY_PATH)/SVA/simu/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE)_*.sorted.bam.SVA_F.bed.count | sort -k2,2 | sed 's/^[ ]*//g' | awk 'BEGIN{first=1} {if ( first == 1 ) {id=$$2;first=0}; if ( id != $$2 ) {print id,sum/count;id=$$2;sum=0;count=0}; sum+=$$1;count++}; END{print id,sum/count;}' | sed 's/_Retroposon__SVA//g' > temp.SVA_F
	paste temp.SVA_B temp.SVA_C temp.SVA_D temp.SVA_E temp.SVA_F | awk '{print $$1,$$2,$$4,$$6,$$8,$$10}' >> $(LIBRARY_PATH)/SVA/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt
	rm temp.SVA_B temp.SVA_C temp.SVA_D temp.SVA_E temp.SVA_F

	rm -Rf $(LIBRARY_PATH)/SVA/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt.lock


##
## Create auxiliary file with proportion of simulated reads on each subfamily
##
$(LIBRARY_PATH)/SVA/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE).prop.txt: $(LIBRARY_PATH)/SVA/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Calculating simulation proportions:\n" >> $(LOG_FILE)
	echo -n "SVA_Subfamily " > $(LIBRARY_PATH)/SVA/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE).prop.txt
	$(R_BIN) --no-restore --no-save --args $(LIBRARY_PATH)/SVA/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt $(LIBRARY_PATH)/SVA/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE).prop.txt < $(LIBRARY_PATH)/SVA/ref/prop.template.r >> $(LOG_FILE)


##
## Create signature file
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/SVA.signatures.txt: $(LIBRARY_PATH)/SVA/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE).prop.txt $(LIBRARY_PATH)/SVA/ref/SVA.bases.ref
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Compiling SVA signature files:\n" >> $(LOG_FILE)
	paste $(LIBRARY_PATH)/SVA/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE).prop.txt $(LIBRARY_PATH)/SVA/ref/SVA.bases.ref | awk '{print $$1,$$5,$$6,$$8}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/SVA.signatures.txt


##
## Quantification of SVA repetitive element reads
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).SVA.count: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bed
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Counting the number of reads on each SVA subfamily:\n" >> $(LOG_FILE)
	echo "SVA_count SVA_Subfamily" > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).SVA.count
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bed | egrep -w "SVA_A_Retroposon__SVA|SVA_B_Retroposon__SVA|SVA_C_Retroposon__SVA|SVA_D_Retroposon__SVA|SVA_E_Retroposon__SVA|SVA_F_Retroposon__SVA" | awk '{print $$(NF-1)}' | sort | uniq -c | sed 's/_Retroposon__SVA//g' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).SVA.count
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).SVA.count | awk '{print $$2,$$1}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).SVA.count.t
	mv $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).SVA.count.t $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).SVA.count 


##
## Correcting the number of reads mapped to SVA
##	
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).SVA.count.corrected: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).SVA.count $(OUTPUT_DIR)/$(SAMPLE_ID)/SVA.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).tpm.factor
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Correcting the number of reads on SVA:\n" >> $(LOG_FILE)
	$(R_BIN) --no-restore --no-save --args $(OUTPUT_DIR)/$(SAMPLE_ID)/SVA.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).SVA.count $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).sorted.bam.tot $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).tpm.factor < $(LIBRARY_PATH)/SVA/ref/lsei.template.r >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Writing SVA quantification files:" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): - $(OUTPUT_DIR)/$(SAMPLE_ID)/SVA.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).SVA.count.corrected" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): - $(OUTPUT_DIR)/$(SAMPLE_ID)/SVA.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).SVA.count.rpkm" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): - $(OUTPUT_DIR)/$(SAMPLE_ID)/SVA.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).SVA.count.rpkm.corrected" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): - $(OUTPUT_DIR)/$(SAMPLE_ID)/SVA.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).SVA.count.signal_proportions" >> $(LOG_FILE)


##     ## ######## ########  ##     ## 
##     ## ##       ##     ## ##     ## 
##     ## ##       ##     ## ##     ## 
######### ######   ########  ##     ## 
##     ## ##       ##   ##    ##   ##  
##     ## ##       ##    ##    ## ##   
##     ## ######## ##     ##    ###    

$(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt: $(LIBRARY_PATH)/LTR/ref/LTR.ref.fa
ifneq ("$(wildcard $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt.lock)","")
	@echo -e "$(timestamp) $(PIPELINE_NAME): There is another simulation running. Exiting without finishing."
	exit 1
endif

	touch $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt.lock

	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): The profile for this study was not found at: $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Simulating reads with length equal to $(MEAN_READ_LEN)\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Creating reads from based on LTR reference sequence:\n" >> $(LOG_FILE)
	mkdir -p $(LIBRARY_PATH)/LTR/simu/
	@for iter in $(shell seq 1 $(NUMBER_OF_LOOPS) ); do \
		$(WGSIM_BIN) -S $$(date "+%N") -1 $(MEAN_READ_LEN) -N $(NUMBER_OF_READS_LTR) -d0 -r$(ERROR_RATE) -e 0 -R 0 $(LIBRARY_PATH)/LTR/ref/LTR.ref.fa $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.simu /dev/null > /dev/null 2> /dev/null ; \
    done
	@echo -e "$(timestamp) $(PIPELINE_NAME): Aligning simulated reads to the reference genome:\n" >> $(LOG_FILE)
	@for iter in $(shell seq 1 $(NUMBER_OF_LOOPS) ); do \
		$(BOWTIE_BIN) -p $(N_THREADS) $(BOWTIE_PARAMS) -x $(BOWTIE_INDEX) -U $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.simu 2>> $(LOG_FILE) | $(SAMTOOLS_BIN) view -Sb - 2>> $(LOG_FILE) > $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.bam; \
		$(SAMTOOLS_BIN) sort -@$(N_THREADS) -m 4G $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.bam $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted; \
		rm -R $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.bam; \
		$(INTERSERC_BIN) -f 0.75 -a $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam -b $(LIBRARY_PATH)/LTR/ref/LTR.hg38.bed -sorted -bed -wo > $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.LTR.bed; \
		cat $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.LTR.bed | awk -F "[$$\t ]" '{print $$4,$$20}' | sort -k1,1 -k2,2 | uniq -c > $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.LTR.bed.count; \
	done
	@echo -e "$(timestamp) $(PIPELINE_NAME): Calculating the expected number of reads on each subfamily:\n" >> $(LOG_FILE)
	cat $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_*.sorted.bam.LTR.bed.count | sort -k2,2 -k3,3 | sed 's/^[ ]*//g' | awk 'BEGIN{first=1} {if ( first == 1 ) {id=$$2"*"$$3;first=0}; if ( id != $$2"*"$$3 ) {print id,sum/$(NUMBER_OF_LOOPS);id=$$2"*"$$3;sum=0;count=0}; sum+=$$1;count++}; END{print id,sum/$(NUMBER_OF_LOOPS);}' > $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).means.txt
	$(PYTHON_BIN) $(LIBRARY_PATH)/scripts/complete_table.py -1 $(LIBRARY_PATH)/LTR/ref/LTR.bases.ref -2 $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).means.txt > $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt

	rm -Rf $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt.lock

##
## Create auxiliary file with proportion of simulated reads on each subfamily
##
$(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).prop.txt: $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Calculating simulation proportions:\n" >> $(LOG_FILE)
	echo -n "LTR_Subfamily " > $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).prop.txt
	$(R_BIN) --no-restore --no-save --args $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).prop.txt < $(LIBRARY_PATH)/SVA/ref/prop.template.r >> $(LOG_FILE)


##
## Create signature file
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/LTR.signatures.txt: $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).prop.txt $(LIBRARY_PATH)/LTR/ref/LTR.bases.ref
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Compiling SVA signature files:\n" >> $(LOG_FILE)
	cat $(LIBRARY_PATH)/LTR/ref/LTR.bases.ref | awk '{print $$2}' > $(LIBRARY_PATH)/LTR/ref/LTR.bases.ref.tmp
	paste $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).prop.txt $(LIBRARY_PATH)/LTR/ref/LTR.bases.ref.tmp | sed 's/[ \t][ \t]*/ /g' > $(OUTPUT_DIR)/$(SAMPLE_ID)/LTR.signatures.txt
	rm -Rf $(LIBRARY_PATH)/LTR/ref/LTR.bases.ref.tmp


##
## Quantification of LTR repetitive element reads
##
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bed
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Counting the number of reads on each LTR subfamily:\n" >> $(LOG_FILE)
	echo "LTR_count LTR_Subfamily" > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bed | grep "LTR" | awk '{print $$(NF-1)}' | sort | grep "^LTR" | uniq -c | sed 's/_LTR.*//g' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count
	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count | awk '{print $$2,$$1}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count.t
	mv $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count.t $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count 


##
## Correcting the number of reads mapped to LTR
##	
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count.corrected: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count $(OUTPUT_DIR)/$(SAMPLE_ID)/LTR.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).tpm.factor
	@echo -e "======================\n" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Correcting the number of reads on LTR:\n" >> $(LOG_FILE)
	$(R_BIN) --no-restore --no-save --args $(OUTPUT_DIR)/$(SAMPLE_ID)/LTR.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).sorted.bam.tot $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).tpm.factor < $(LIBRARY_PATH)/LTR/ref/lsei.template.r >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): Writing LTR quantification files:" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): - $(OUTPUT_DIR)/$(SAMPLE_ID)/LTR.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count.corrected" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): - $(OUTPUT_DIR)/$(SAMPLE_ID)/LTR.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count.rpkm" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): - $(OUTPUT_DIR)/$(SAMPLE_ID)/LTR.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count.rpkm.corrected" >> $(LOG_FILE)
	@echo -e "$(timestamp) $(PIPELINE_NAME): - $(OUTPUT_DIR)/$(SAMPLE_ID)/LTR.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count.signal_proportions" >> $(LOG_FILE)


##
## Main sub-target
##
#processSample: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count.corrected $(LIBRARY_PATH)/SVA/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt
processSample: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count.corrected $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).SVA.count.corrected $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count.corrected
#	## Copy Output descriptions file
#	#cp $(SRNABENCH_LIBS)/sRNAbenchOutputDescription.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/sRNAbenchOutputDescription.txt 
#	## END PIPELINE
#	#@echo -e "$(ts) SMRNAPIPELINE: END smallRNA-seq Pipeline for sample $(SAMPLE_ID)\n======================\n" >> $(LOG_FILE)


