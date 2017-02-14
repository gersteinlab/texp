options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

#
sig_file = args[1]
L1HS_hg38sub_counts_file = args[2]
total_reads_file = args[3]
#
length = 6000

A=read.table(sig_file,header=T,row.names=1);  
B=read.table(L1HS_hg38sub_counts_file,header=T,row.names=1);   #Number of reads at each subfamily

#Select/order elements common in the experiment and simulations
names=c(rownames(A),rownames(B))
ok = names[duplicated(names)]
A=A[ok,,drop=FALSE]
B=B[ok,,drop=FALSE]

#Clean elements that were not contemplated in the simulations or wo reads
B2=merge(B,A,by=0,all=TRUE)[,1:3];
B2=B2[complete.cases(B2),,drop=FALSE]
B=B[B2$Row.names,,drop=FALSE]
A=A[B2$Row.names,,drop=FALSE]

if ( sum(B/length>=0.01) != 0 ) {
n_col_sig = ncol(A); 
library(penalized)
profile_lambda1=profL1(B[,1],penalized=A,unpenalized=~0,lambda2=0,positive=TRUE,standardize=TRUE,plot=TRUE,minsteps=10000,maxiter=1000); 
optimized_lambda1=optL1(B[,1],penalized=A,unpenalized=~0,lambda2=0,positive=TRUE,standardize=TRUE,fold=profile_lambda1$fold);

percentages = coefficients(optimized_lambda1$fullfit,"all")/sum(coefficients(optimized_lambda1$fullfit,"all"))
corrected_L1HS_hg38_reads = data.frame(coefficients(optimized_lambda1$fullfit,"all"))

L1HS_hg38_reads = B
#corrected_L1HS_hg38_reads = data.frame(percentages*sum(B))
rownames(corrected_L1HS_hg38_reads) = colnames(A)
colnames(corrected_L1HS_hg38_reads) = c("reads")

tot = scan(total_reads_file)

##
## Quantify element transcription
##
#L1HS_hg38_reads = L1HS_hg38_reads[1:dim(L1HS_hg38_reads)[1]-1,,drop=FALSE]
corrected_L1HS_hg38_reads.back = corrected_L1HS_hg38_reads
corrected_L1HS_hg38_reads = corrected_L1HS_hg38_reads[1:nrow(corrected_L1HS_hg38_reads)-1,,drop=FALSE]

rpkm = data.frame((L1HS_hg38_reads*10^9)/(length*tot))
colnames(rpkm)=c("RPKM")

rpkm.corrected = data.frame((corrected_L1HS_hg38_reads*10^9)/(length*tot))
colnames(rpkm.corrected)=c("RPKM")

##
## Dump results into file
##
write.table(percentages,file=paste(L1HS_hg38sub_counts_file,"signal_proportions",sep="."))

write.table(rpkm,file=paste(L1HS_hg38sub_counts_file,"rpkm",sep="."),quote=F)
write.table(rpkm.corrected,file=paste(L1HS_hg38sub_counts_file,"rpkm.corrected",sep="."),quote=F)

write.table(corrected_L1HS_hg38_reads.back,file=paste(L1HS_hg38sub_counts_file,"corrected",sep="."))
} else {
message("Not enough reads in L1HS_hg38 elements to perform the correction\n",appendLF=FALSE)
}
