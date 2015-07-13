options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

#
signature_file = args[1]
l1sub_counts_file = args[2]
total_reads_file = args[3]
#

A = read.table(signature_file);   #Expected number of reads based Assuming background transcription and simulation of L1Hs transcript;
n_col_sig = ncol(A); 
B = read.table(l1sub_counts_file);   #Number of reads on each L1 subfamily
F = 1                 #Sum of percentages has to be equal to 1
G = diag(n_col_sig)   #
E = rep(1,n_col_sig)  #
H = rep(0,n_col_sig)  #
library(limSolve)
percentages = lsei(A,B,E,F,G,H)

##
## Quantify L1Hs transcription
##
tot = scan(total_reads_file)
length = 6000
reads = B[2,1]
l1hs_percent_signal = percentages[1,]
l1hs_mismapping_percentage = A[1,2]     
corrected_L1Hs_reads = reads*l1hs_percent_signal*(1/l1hs_mismapping_percentage)
rpkm = (reads*10^9)/(length*tot)
rpkm.corrected = (corrected_L1Hs_reads*10^9)/(length*tot)
