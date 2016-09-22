options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

#
sig_file = args[1]
L1HS_hg38sub_counts_file = args[2]
total_reads_file = args[3]
factor_tpm_file = args[4]
#

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

n_col_sig = ncol(A); 
library(glmnet)
penalty.factor=rep(1,n_col_sig)
penalty.factor[n_col_sig]=0.01
fit = glmnet( x=as.matrix(A), y=c(B[,1]), lower.limits=rep(0,n_col_sig),intercept=FALSE,standardize=FALSE,alpha=1,penalty.factor=penalty.factor)
lambdas = seq(from=fit$lambda[1],to=100,by=-fit$lambda[1]*0.001)
fit = update(fit,nlambda=length(lambdas),lambda=lambdas)
coeflasso = as.vector(coef(fit,s=0))
percentages = coeflasso[-1]/sum(coeflasso[-1])

L1HS_hg38_reads = B
corrected_L1HS_hg38_reads = data.frame(percentages*sum(B))
rownames(corrected_L1HS_hg38_reads) = colnames(A)
colnames(corrected_L1HS_hg38_reads) = c("reads")

tot = scan(total_reads_file)
factor_tpm = scan(factor_tpm_file)
length = 6000

##
## Quantify element transcription
##
L1HS_hg38_reads = L1HS_hg38_reads[1:dim(L1HS_hg38_reads)[1]-1,,drop=FALSE]
corrected_L1HS_hg38_reads = corrected_L1HS_hg38_reads[1:nrow(corrected_L1HS_hg38_reads)-1,,drop=FALSE]

rpkm = data.frame((L1HS_hg38_reads*10^9)/(length*tot))
colnames(rpkm)=c("RPKM")

rpkm.corrected = data.frame((corrected_L1HS_hg38_reads*10^9)/(length*tot))
colnames(rpkm.corrected)=c("RPKM")

tpm = data.frame((L1HS_hg38_reads/length)*(1/factor_tpm)*10^6)
colnames(tpm)=c("TPM")

tpm.corrected = data.frame((corrected_L1HS_hg38_reads/length)*(1/factor_tpm)*10^6)
colnames(tpm.corrected)=c("TPM")

##
## Dump results into file
##
write.table(percentages,file=paste(L1HS_hg38sub_counts_file,"signal_proportions",sep="."))

write.table(rpkm,file=paste(L1HS_hg38sub_counts_file,"rpkm",sep="."),quote=F)
write.table(rpkm.corrected,file=paste(L1HS_hg38sub_counts_file,"rpkm.corrected",sep="."),quote=F)
write.table(tpm,file=paste(L1HS_hg38sub_counts_file,"tpm",sep="."),quote=F)
write.table(tpm.corrected,file=paste(L1HS_hg38sub_counts_file,"tpm.corrected",sep="."),quote=F)

write.table(corrected_L1HS_hg38_reads,file=paste(L1HS_hg38sub_counts_file,"corrected",sep="."))

