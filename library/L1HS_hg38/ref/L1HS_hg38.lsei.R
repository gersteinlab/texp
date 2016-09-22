options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

#
signature_file = args[1]
L1HS_hg38sub_counts_file = args[2]
total_reads_file = args[3]
factor_tpm_file = args[4]
#

A=read.table(signature_file,header=T,row.names=1);   
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

#Prepare/Run lsei
A.scale=scale(A,center=FALSE);
A.scale[is.nan(A.scale)] = 0
B.scale=scale(B,center=FALSE);

n_col_sig = ncol(A); 
E = rep(1,n_col_sig)  #
F = 1                 # (Ex=F) Sum of percentages has to be equal to 1

G = diag(n_col_sig)   #
H = rep(0,n_col_sig)  # Gx >= H

library(limSolve)
percentages = lsei(A.scale,B.scale,E,F,G,H)
#B.corrected = B*percentages$X

##
## Quantify element transcription
##
tot = scan(total_reads_file)
factor_tpm = scan(factor_tpm_file)
length = 6000
reads = B
L1HS_hg38_percent_signal = percentages$X[1:(length(percentages$X)-1)]

L1HS_hg38_reads = data.frame(x=numeric(dim(B)[1]-1));
corrected_L1HS_hg38_reads = data.frame(x=numeric(dim(B)[1]-1));


count=1; 
for (i in rownames(A)) { 
	j=paste(sub("-",".",i),"_Transcript",sep=""); 
	if ( !is.na(B[i,]) &&  !is.null(A[i,j]) && !is.na(L1HS_hg38_percent_signal[j]) ) {
		L1HS_hg38_reads[count,] =  B[i,]
		rownames(L1HS_hg38_reads)[count] = i
		signal = A.scale[,j]*L1HS_hg38_percent_signal[j]*attr(B.scale, 'scaled:scale')
		reads_from_signal = signal[i]*(1/A[i,j])
		corrected_L1HS_hg38_reads[count,] = reads_from_signal
#		corrected_L1HS_hg38_reads[count,] = B[i,]*(1/A[i,j])*
		rownames(corrected_L1HS_hg38_reads)[count] = i
		count=count+1;
	}
}
#L1HS_hg38_reads = L1HS_hg38_reads[1:count-1,,drop=FALSE]
#corrected_L1HS_hg38_reads = corrected_L1HS_hg38_reads[1:count-1,,drop=FALSE]

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
write.table(percentages$X,file=paste(L1HS_hg38sub_counts_file,"signal_proportions",sep="."))

write.table(rpkm,file=paste(L1HS_hg38sub_counts_file,"rpkm",sep="."),quote=F)
write.table(rpkm.corrected,file=paste(L1HS_hg38sub_counts_file,"rpkm.corrected",sep="."),quote=F)
write.table(tpm,file=paste(L1HS_hg38sub_counts_file,"tpm",sep="."),quote=F)
write.table(tpm.corrected,file=paste(L1HS_hg38sub_counts_file,"tpm.corrected",sep="."),quote=F)

write.table(corrected_L1HS_hg38_reads,file=paste(L1HS_hg38sub_counts_file,"corrected",sep="."))

