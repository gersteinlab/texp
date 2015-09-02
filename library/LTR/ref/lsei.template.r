options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

#
signature_file = args[1]
ltrsub_counts_file = args[2]
total_reads_file = args[3]
factor_tpm_file = args[4]
#

A = read.table(signature_file,header=T,row.names=1);   #Expected number of reads based Assuming background transcription and simulation of L1Hs transcript;
#A = A[,c("LTR_Ref_bases","LTR5_Hs_Transcript","LTR7Y_Transcript","LTR7C_Transcript","LTR7Y_Transcript","LTR7C_Transcript","LTR7B_Transcript","LTR12C_Transcript")]
incompatible_elements = match(NA, A$LTR5_Hs_Transcript)
incompatible_elements_id = row.names(A[incompatible_elements,])
if ( !is.na(incompatible_elements[1]) ){
A=A[-incompatible_elements,,drop=FALSE]
}
A=A[order(row.names(A)),,drop=FALSE]
A.scale=scale(A)

n_col_sig = ncol(A); 
B = read.table(ltrsub_counts_file,header=T,row.names=1);   #Number of reads on each L1 subfamily
B2=merge(B,A,by=0,all=TRUE)[,1:2];B2[is.na(B2)]=0;rownames(B2)=B2[,1];B2$Row.names=NULL;B=B2
B_incompatible_elements = match(incompatible_elements_id,rownames(B))
if ( !is.na(B_incompatible_elements[1]) ){
B=B[-B_incompatible_elements,,drop=FALSE]; 
}
B=B[order(row.names(B)),,drop=FALSE]
B.scale = scale(B);
F = 1                 #Sum of percentages has to be equal to 1
G = diag(n_col_sig)   #
E = rep(1,n_col_sig)  #
H = rep(0,n_col_sig)  #
library(limSolve)
percentages = lsei(A.scale,B.scale,E,F,G,H)
B.corrected = B*percentages$X

##
## Quantify L1Hs transcription
##
tot = scan(total_reads_file)
factor_tpm = scan(factor_tpm_file)
length = 1000
reads = B
ltr_percent_signal = percentages$X
ltr_mismapping_percentage=data.frame(x=numeric(dim(reads)[1])-1); count=0; for (i in rownames(A)){ j=paste(sub("-",".",i),"_Transcript",sep=""); ltr_mismapping_percentage[count,]=A[i,j]; count=count+1}    
corrected_ltr_reads = reads*(1/ltr_mismapping_percentage)*ltr_percent_signal

rpkm = data.frame((reads*10^9)/(length*tot))
colnames(rpkm)=c("RPKM")


rpkm.corrected = data.frame((corrected_ltr_reads*10^9)/(length*tot))
colnames(rpkm.corrected)=c("RPKM")

tpm = data.frame((reads/length)*(1/factor_tpm)*10^6)
colnames(tpm)=c("TPM")

tpm.corrected = data.frame((corrected_ltr_reads/length)*(1/factor_tpm)*10^6)
colnames(tpm.corrected)=c("TPM")

##
## Dump results into file
##
write.table(percentages$X,file=paste(ltrsub_counts_file,"signal_proportions",sep="."))

write.table(rpkm,file=paste(ltrsub_counts_file,"rpkm",sep="."),quote=F)
write.table(rpkm.corrected,file=paste(ltrsub_counts_file,"rpkm.corrected",sep="."),quote=F)
write.table(tpm,file=paste(ltrsub_counts_file,"tpm",sep="."),quote=F)
write.table(tpm.corrected,file=paste(ltrsub_counts_file,"tpm.corrected",sep="."),quote=F)

write.table(B.corrected,file=paste(ltrsub_counts_file,"corrected",sep="."))

