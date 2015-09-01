PIPELINE_NAME     := TeXP

#DATA_DIR          := NULL
#OUTPUT_DIR        := NULL
#INPUT_FILE_PATH   := NULL
#SAMPLE_NAME       := NULL
REFERENCE_GENOME  := NULL

##
## Use the input path to infer filetype and short name
##
#INPUT_FILE_NAME := $(notdir $(INPUT_FILE_PATH))
#INPUT_FILE_ID   := $(basename $(INPUT_FILE_NAME))

LIBRARY_PATH     := /home2/fn64/projects/TeXP/library
EXE_DIR          := /home2/fn64/tools/manual

BOWTIE_BIN       := $(EXE_DIR)/bin/bowtie2
BOWTIE_PARAMS    := --sensitive-local -N1 --no-unal
BOWTIE_INDEX     := "/home2/fn64/genomes/Homo_sapiens/hg38/toBowtie2/hg38"

R_BIN            := $(EXE_DIR)/bin/R
BEDTOOLS_BIN     := $(EXE_DIR)/bin/bedtools
WGSIM_BIN        := $(EXE_DIR)/bin/wgsim
PYTHON_BIN       := /usr/bin/python

REPEAT_MASKER_OUT     := 
REPEAT_MASKER_BED     := ~/projects/hg38.rep.noexon.bed
REPEAT_MASKER_TOT_BED := ~/projects/hg38.rep.bed

#COMMAND_HOMOPOL := perl $(LIBRARY_PATH)/scripts/remove_homopol.pl
#COMMAND_PARTIAL := perl $(LIBRARY_PATH)/scripts/filter_qual.pl

##
## Simulation parameters
##
ERROR_RATE          := 0.1
NUMBER_OF_READS     := 25000
NUMBER_OF_LOOPS     := 100

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

LOG_FILE := $(ELE_NAME).log

## Define current time
timestamp := `/bin/date "+%Y-%m-%d(%H:%M:%S)"`

##
## Main make target
##
.PHONY: all
.DEFAULT: all
.INTERMEDIATE: $(ELE_NAME)/*.ref.bed.tmp

all: processElement

##
## Compile a sub-bed file containing only the desired elements.
##
$(ELE_NAME)/ref/$(ELE_NAME).bed: $(REPEAT_MASKER_TOT_BED)
	@echo -e "\n\n\n$(timestamp) **************************************************************" >> $(ELE_NAME).log
	@echo -e "$(timestamp) $(PIPELINE_NAME): Creating library folder: ./$(ELE_NAME)" >> $(ELE_NAME).log
	@mkdir -p $(ELE_NAME)/ref/
	@mkdir -p $(ELE_NAME)/info/
	@echo -e "$(timestamp) $(PIPELINE_NAME): Creating sub bed file containing only the element of interestet:" >> $(ELE_NAME).log
	egrep $(ELE_NAME) $(REPEAT_MASKER_TOT_BED) > $(ELE_NAME)/ref/$(ELE_NAME).bed

##
## Compile a sub-bed file containing only the desired elements.
##
$(ELE_NAME)/$(ELE_NAME).bed.elements: $(ELE_NAME)/ref/$(ELE_NAME).bed
	@echo -e "$(timestamp) $(PIPELINE_NAME): Elements to be processed:" >> $(ELE_NAME).log
	@cat $(ELE_NAME)/ref/$(ELE_NAME).bed | awk '{print $$NF}' | sort | uniq >> $(ELE_NAME).log 
	@cat $(ELE_NAME)/ref/$(ELE_NAME).bed | awk '{print $$NF}' | sort | uniq > $(ELE_NAME)/$(ELE_NAME).bed.elements
	@echo -e "$(timestamp) $(PIPELINE_NAME): *** WARNING ***:" >> $(ELE_NAME).log
	@echo -e "$(timestamp) $(PIPELINE_NAME): If there are undesired elements, remove it manually from ./$(ELE_NAME)/$(ELE_NAME).bed.elements and rerun librarian Makefile:" >> $(ELE_NAME).log
	@echo -e "$(timestamp) $(PIPELINE_NAME): ***         ***:" >> $(ELE_NAME).log

##
## Calculate the proportion of bases for each element.
##
$(ELE_NAME)/ref/$(ELE_NAME).ref.bases: $(ELE_NAME)/$(ELE_NAME).bed.elements
	@echo -e "$(timestamp) $(PIPELINE_NAME): Creating element summary:" >> $(ELE_NAME).log
	cat $(ELE_NAME)/ref/$(ELE_NAME).bed | sort -k4,4 | awk 'BEGIN{first=1} {if ( first == 1 ) {id=$$4;first=0}; if ( id != $$4 ) {print id,count,sum;id=$$4;sum=0;count=0}; sum+=($$3-$$2);count+=1}; END{print id,count,sum;}' | fgrep -w -f $(ELE_NAME)/$(ELE_NAME).bed.elements > $(ELE_NAME)/$(ELE_NAME).summary

	@echo -e "$(timestamp) $(PIPELINE_NAME): Calculating the proportion of bases of each subfamily:" >> $(ELE_NAME).log
	@echo -e "$(ELE_NAME)_Subfamily\t$(ELE_NAME)_Ref_bases" > $(ELE_NAME)/ref/$(ELE_NAME).ref.bases
	SUM=$$(cat $(ELE_NAME)/$(ELE_NAME).summary | awk '{sum+=$$3} END{print sum}'); \
	cat $(ELE_NAME)/$(ELE_NAME).summary | awk -v sum=$$SUM '{print $$1"\t"$$3/sum}' >> $(ELE_NAME)/ref/$(ELE_NAME).ref.bases

##
## Fetch bases the the reference genome and dump into the ref.fa file.
##
$(ELE_NAME)/ref/$(ELE_NAME).ref.fa: $(ELE_NAME)/$(ELE_NAME).bed.elements
	@echo -e "$(timestamp) $(PIPELINE_NAME): Fetching the reference sequence:" >> $(ELE_NAME).log
	for element in $$(cat $(ELE_NAME)/$(ELE_NAME).bed.elements); do \
	cat $(REPEAT_MASKER_TOT_BED) | awk -v element=$$element '{if ($$4 == element) {print}}' > $(ELE_NAME)/$$element.ref.bed.tmp; \
	$(BEDTOOLS_BIN) getfasta -fi $(REFERENCE_GENOME) -bed $(ELE_NAME)/$$element.ref.bed.tmp -fo $(ELE_NAME)/ref/$$element.ref.fasta.tmp -s -name; \
	done

	@cat $(ELE_NAME)/ref/*.ref.fasta.tmp > $(ELE_NAME)/ref/$(ELE_NAME).ref.fa
	rm -Rf $(ELE_NAME)/*.ref.bed.tmp

##
## Fetch bases the the reference genome and dump into the ref.fa file.
##
$(ELE_NAME)/info/copies_bases.plot.pdf: $(ELE_NAME)/$(ELE_NAME).summary
	sed 's/TEMPLATE/$(ELE_NAME)/g' templates/copies_bases.R.template > $(ELE_NAME)/info/copies_bases.R
	$(R_BIN) --no-restore --no-save --args $(ELE_NAME)/$(ELE_NAME).summary < $(ELE_NAME)/info/copies_bases.R >> $(LOG_FILE)


processElement: $(ELE_NAME)/ref/$(ELE_NAME).ref.bases $(ELE_NAME)/ref/$(ELE_NAME).ref.fa $(ELE_NAME)/info/copies_bases.plot.pdf
## $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt: $(LIBRARY_PATH)/LTR/ref/LTR.ref.fa
## ifneq ("$(wildcard $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt.lock)","")
## 	@echo -e "$(timestamp) $(PIPELINE_NAME): There is another simulation running. Exiting without finishing."
## 	exit 1
## endif
## 
## 	touch $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt.lock
## 
## 	@echo -e "======================\n" >> $(LOG_FILE)
## 	@echo -e "$(timestamp) $(PIPELINE_NAME): The profile for this study was not found at: $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt\n" >> $(LOG_FILE)
## 	@echo -e "$(timestamp) $(PIPELINE_NAME): Simulating reads with length equal to $(MEAN_READ_LEN)\n" >> $(LOG_FILE)
## 	@echo -e "$(timestamp) $(PIPELINE_NAME): Creating reads from based on LTR reference sequence:\n" >> $(LOG_FILE)
## 	mkdir -p $(LIBRARY_PATH)/LTR/simu/
## 	@for iter in $(shell seq 1 $(NUMBER_OF_LOOPS) ); do \
## 		$(WGSIM_BIN) -S $$(date "+%N") -1 $(MEAN_READ_LEN) -N $(NUMBER_OF_READS_LTR) -d0 -r$(ERROR_RATE) -e 0 -R 0 $(LIBRARY_PATH)/LTR/ref/LTR.ref.fa $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.simu /dev/null > /dev/null 2> /dev/null ; \
##     done
## 	@echo -e "$(timestamp) $(PIPELINE_NAME): Aligning simulated reads to the reference genome:\n" >> $(LOG_FILE)
## 	@for iter in $(shell seq 1 $(NUMBER_OF_LOOPS) ); do \
## 		$(BOWTIE_BIN) -p $(N_THREADS) $(BOWTIE_PARAMS) -x $(BOWTIE_INDEX) -U $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.simu 2>> $(LOG_FILE) | $(SAMTOOLS_BIN) view -Sb - 2>> $(LOG_FILE) > $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.bam; \
## 		$(SAMTOOLS_BIN) sort -@$(N_THREADS) -m 4G $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.bam $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted; \
## 		rm -R $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.bam; \
## 		$(INTERSERC_BIN) -f 0.75 -a $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam -b $(LIBRARY_PATH)/LTR/ref/LTR.hg38.bed -sorted -bed -wo > $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.LTR.bed; \
## 		cat $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.LTR.bed | awk -F "[$$\t ]" '{print $$4,$$20}' | sort -k1,1 -k2,2 | uniq -c > $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_$$iter.sorted.bam.LTR.bed.count; \
## 	done
## 	@echo -e "$(timestamp) $(PIPELINE_NAME): Calculating the expected number of reads on each subfamily:\n" >> $(LOG_FILE)
## 	cat $(LIBRARY_PATH)/LTR/simu/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE)_*.sorted.bam.LTR.bed.count | sort -k2,2 -k3,3 | sed 's/^[ ]*//g' | awk 'BEGIN{first=1} {if ( first == 1 ) {id=$$2"*"$$3;first=0}; if ( id != $$2"*"$$3 ) {print id,sum/$(NUMBER_OF_LOOPS);id=$$2"*"$$3;sum=0;count=0}; sum+=$$1;count++}; END{print id,sum/$(NUMBER_OF_LOOPS);}' > $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).means.txt
## 	$(PYTHON_BIN) $(LIBRARY_PATH)/scripts/complete_table.py -1 $(LIBRARY_PATH)/LTR/ref/LTR.bases.ref -2 $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS)_$(MEAN_READ_LEN)_$(ERROR_RATE).means.txt > $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt
## 
## 	rm -Rf $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt.lock
## 
## ##
## ## Create auxiliary file with proportion of simulated reads on each subfamily
## ##
## $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).prop.txt: $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt
## 	@echo -e "======================\n" >> $(LOG_FILE)
## 	@echo -e "$(timestamp) $(PIPELINE_NAME): Calculating simulation proportions:\n" >> $(LOG_FILE)
## 	echo -n "SVA_Subfamily " > $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).prop.txt
## 	$(R_BIN) --no-restore --no-save --args $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).prop.txt < $(LIBRARY_PATH)/SVA/ref/prop.template.r >> $(LOG_FILE)
## 
## 
## ##
## ## Create signature file
## ##
## $(OUTPUT_DIR)/$(SAMPLE_ID)/LTR.signatures.txt: $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).prop.txt $(LIBRARY_PATH)/LTR/ref/LTR.bases.ref
## 	@echo -e "======================\n" >> $(LOG_FILE)
## 	@echo -e "$(timestamp) $(PIPELINE_NAME): Compiling SVA signature files:\n" >> $(LOG_FILE)
## 	cat $(LIBRARY_PATH)/LTR/ref/LTR.bases.ref | awk '{print $$2}' > $(LIBRARY_PATH)/LTR/ref/LTR.bases.ref.tmp
## 	paste $(LIBRARY_PATH)/LTR/$(NUMBER_OF_READS_LTR)_$(MEAN_READ_LEN)_$(ERROR_RATE).prop.txt $(LIBRARY_PATH)/LTR/ref/LTR.bases.ref.tmp | sed 's/[ \t][ \t]*/ /g' > $(OUTPUT_DIR)/$(SAMPLE_ID)/LTR.signatures.txt
## 	rm -Rf $(LIBRARY_PATH)/LTR/ref/LTR.bases.ref.tmp
## 
## 
## ##
## ## Quantification of LTR repetitive element reads
## ##
## $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bed
## 	@echo -e "======================\n" >> $(LOG_FILE)
## 	@echo -e "$(timestamp) $(PIPELINE_NAME): Counting the number of reads on each LTR subfamily:\n" >> $(LOG_FILE)
## 	echo "LTR_count LTR_Subfamily" > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count
## 	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).re.filtered.bed | grep "LTR" | awk '{print $$(NF-1)}' | sort | grep "^LTR" | uniq -c | sed 's/_LTR.*//g' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count
## 	cat $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count | awk '{print $$2,$$1}' > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count.t
## 	mv $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count.t $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count 
## 
## 
## ##
## ## Correcting the number of reads mapped to LTR
## ##	
## $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count.corrected: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count $(OUTPUT_DIR)/$(SAMPLE_ID)/LTR.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).tpm.factor
## 	@echo -e "======================\n" >> $(LOG_FILE)
## 	@echo -e "$(timestamp) $(PIPELINE_NAME): Correcting the number of reads on LTR:\n" >> $(LOG_FILE)
## 	$(R_BIN) --no-restore --no-save --args $(OUTPUT_DIR)/$(SAMPLE_ID)/LTR.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).sorted.bam.tot $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).tpm.factor < $(LIBRARY_PATH)/LTR/ref/lsei.template.r >> $(LOG_FILE)
## 	@echo -e "$(timestamp) $(PIPELINE_NAME): Writing LTR quantification files:" >> $(LOG_FILE)
## 	@echo -e "$(timestamp) $(PIPELINE_NAME): - $(OUTPUT_DIR)/$(SAMPLE_ID)/LTR.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count.corrected" >> $(LOG_FILE)
## 	@echo -e "$(timestamp) $(PIPELINE_NAME): - $(OUTPUT_DIR)/$(SAMPLE_ID)/LTR.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count.rpkm" >> $(LOG_FILE)
## 	@echo -e "$(timestamp) $(PIPELINE_NAME): - $(OUTPUT_DIR)/$(SAMPLE_ID)/LTR.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count.rpkm.corrected" >> $(LOG_FILE)
## 	@echo -e "$(timestamp) $(PIPELINE_NAME): - $(OUTPUT_DIR)/$(SAMPLE_ID)/LTR.signatures.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count.signal_proportions" >> $(LOG_FILE)
## 
## 
## ##
## ## Main sub-target
## ##
## #processSample: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count.corrected $(LIBRARY_PATH)/SVA/$(NUMBER_OF_READS_SVA)_$(MEAN_READ_LEN)_$(ERROR_RATE).txt
## processSample: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).L1.count.corrected $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).SVA.count.corrected $(OUTPUT_DIR)/$(SAMPLE_ID)/$(SAMPLE_ID).LTR.count.corrected
## #	## Copy Output descriptions file
## #	#cp $(SRNABENCH_LIBS)/sRNAbenchOutputDescription.txt $(OUTPUT_DIR)/$(SAMPLE_ID)/sRNAbenchOutputDescription.txt 
## #	## END PIPELINE
## #	#@echo -e "$(ts) SMRNAPIPELINE: END smallRNA-seq Pipeline for sample $(SAMPLE_ID)\n======================\n" >> $(LOG_FILE)
## 
##
