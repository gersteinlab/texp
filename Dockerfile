FROM debian:stretch
RUN apt-get update

################
#Install binaries dependencies
################

RUN apt-get install -y \
	bedtools=2.26.0+dfsg-3 \
	bowtie2=2.3.0-2 \
	fastx-toolkit=0.0.14-3 \
	gawk=1:4.1.4+dfsg-1 \
	git \
	perl \
	python=2.7.13-2 \
	r-base=3.3.3-1 \
	r-base-dev=3.3.3-1 \
	samtools=1.3.1-3 \
	wget 


################
#Install Wgsim
################

RUN	mkdir -p /src; \ 
	cd /src ; \
	git clone https://github.com/lh3/wgsim.git; \
	cd wgsim; \
	gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm; \
	mv wgsim /usr/bin/; \
	cd /;

################
#Install TeXP
################

RUN	mkdir -p /src; \ 
	cd /src ; \
	git clone https://github.com/gersteinlab/texp.git
RUN	ln -s /src/TeXP/TeXP.sh /usr/bin/TeXP; ln -s /src/TeXP/TeXP.sh /usr/bin/TeXP.sh


################
#Download Libraries
################

RUN mkdir -p /data/library/rep_annotation; \
	cd /data/library/rep_annotation; \
	wget -c -t0 "http://homes.gersteinlab.org/people/fn64/TeXP/rep_annotation.hg38.tar.bz2" -O rep_annotation.hg38.tar.bz2; \
	tar xjvf rep_annotation.hg38.tar.bz2; \
	rm -Rf rep_annotation.hg38.tar.bz2
	
# DEPRECATED: 2017-02-14 - Only L1 subfamilies RPKM are estimated in the latest version		
#RUN mkdir -p /data/library/kallisto; \
#	cd /data/library/; \
#	wget -c -t0 "https://www.dropbox.com/s/bd0vqqyedx7p0jv/kallisto_gencode23.hg38.tar.bz2?dl=1" -O kallisto_gencode23.hg38.tar.bz2; \
#	tar xjvf kallisto_gencode23.hg38.tar.bz2; \
#	rm -Rf kallisto_gencode23.hg38.tar.bz2

RUN mkdir -p /data/library/bowtie2; \
	cd /data/library/bowtie2; \
	wget -c -t0 "http://homes.gersteinlab.org/people/fn64/TeXP/bowtie2.hg38.tar.bz2" -O bowtie2.hg38.tar.bz2; \
	tar xjvf bowtie2.hg38.tar.bz2; \
	rm -Rf bowtie2.hg38.tar.bz2


################
#Install R packages dependencies
################
#RUN cd /data/library/; \
#	wget -c -t0 "https://cloud.r-project.org/src/contrib/Archive/penalized/penalized_0.9-49.tar.gz" -O penalized_0.9-49.tar.gz; \
#	R CMD INSTALL penalized_0.9-49.tar.gz

RUN echo 'install.packages(c("penalized"), repos="http://cloud.r-project.org", dependencies=TRUE)' > /tmp/packages.R \
    && Rscript /tmp/packages.R


#make -f Makefile INPUT_FILE_PATH=/data/ENCFF000EBE.fastq.gz OUTPUT_DIR=/data/process/ N_THREADS=4 SAMPLE_NAME="t"

WORKDIR /src/TeXP/
CMD ["/src/TeXP/TeXP.sh"] 

