options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

#
signature_file = args[1]
svasub_counts_file = args[2]
total_reads_file = args[3]
factor_tpm_file = args[4]
#

A = read.table(signature_file,header=T,row.names=1);   #Expected number of reads based Assuming background transcription and simulation of L1Hs transcript;
A.scale=scale(A)
n_col_sig = ncol(A); 
B = read.table(svasub_counts_file,header=T,row.names=1);   #Number of reads on each L1 subfamily
B.scale = scale(B);
F = 1                 #Sum of percentages has to be equal to 1
G = diag(n_col_sig)   #
E = rep(1,n_col_sig)  #
H = rep(0,n_col_sig)  #
library(limSolve)
percentages = lsei(A.scale,B.scale,E,F,G,H)
B.corrected = B[5:6,]*percentages$X[1:2]

##
## Quantify L1Hs transcription
##
tot = scan(total_reads_file)
factor_tpm = scan(factor_tpm_file)
length = 1500
reads = B[5:6,]
svahs_percent_signal = percentages$X[1:2]
svahs_mismapping_percentage = cbind(A[5,1],A[6,2])    
corrected_sva_reads = reads*(1/svahs_mismapping_percentage)*svahs_percent_signal

rpkm = data.frame((reads*10^9)/(length*tot))
colnames(rpkm)=c("RPKM")
row.names(rpkm)=c("SVA_E","SVA_F")

rpkm.corrected = data.frame(t((corrected_sva_reads*10^9)/(length*tot)))
colnames(rpkm.corrected)=c("RPKM")
row.names(rpkm.corrected)=c("SVA_E","SVA_F")

tpm = data.frame((reads/length)*(1/factor_tpm)*10^6)
colnames(rpkm.corrected)=c("TPM")
row.names(rpkm.corrected)=c("SVA_E","SVA_F")

tpm.corrected = data.frame(t((corrected_sva_reads/length)*(1/factor_tpm)*10^6))
colnames(rpkm.corrected)=c("TPM")
row.names(rpkm.corrected)=c("SVA_E","SVA_F")

##
## Dump results into file
##
write.table(percentages$X,file=paste(svasub_counts_file,"signal_proportions",sep="."))

write.table(rpkm,file=paste(svasub_counts_file,"rpkm",sep="."),quote=F)
write.table(rpkm.corrected,file=paste(svasub_counts_file,"rpkm.corrected",sep="."),quote=F)
write.table(tpm,file=paste(svasub_counts_file,"tpm",sep="."),quote=F)
write.table(tpm.corrected,file=paste(svasub_counts_file,"tpm.corrected",sep="."),quote=F)

write.table(B.corrected,file=paste(svasub_counts_file,"corrected",sep="."))

