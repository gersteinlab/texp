options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

#
signature_file = args[1]
l1sub_counts_file = args[2]
total_reads_file = args[3]
#

A = read.table(signature_file,header=T,row.names=1);   #Expected number of reads based Assuming background transcription and simulation of L1Hs transcript;
A.scale=scale(A)
n_col_sig = ncol(A); 
B = read.table(l1sub_counts_file,header=T,row.names=1);   #Number of reads on each L1 subfamily
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
length = 6000
reads = B[1,1]
l1hs_percent_signal = percentages$X[1]
l1hs_mismapping_percentage = A[1,1]     
corrected_L1Hs_reads = reads*l1hs_percent_signal*(1/l1hs_mismapping_percentage)
rpkm = (reads*10^9)/(length*tot)
rpkm.corrected = (corrected_L1Hs_reads*10^9)/(length*tot)

##
## Dump results into file
##
write.table(percentages$X,file=paste(l1sub_counts_file,"signal_proportions",sep="."))

write(rpkm,file=paste(l1sub_counts_file,"rpkm",sep="."))
write(rpkm.corrected,file=paste(l1sub_counts_file,"rpkm.corrected",sep="."))

write.table(B.corrected,file=paste(l1sub_counts_file,"corrected",sep="."))

