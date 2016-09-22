options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

#
sig_file = args[1]
L1HS_hg38sub_counts_file = args[2]
total_reads_file = args[3]
factor_tpm_file = args[4]
#
library(mixedMem)

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
A=A[,colSums(A)!=0] 

A[A==0]=1e-10;  
sig=A
counts=B;   #Number of reads at each subfamily

#Number of samples
Total=1

# Number of Sub-populations
K = dim(sig)[2]

#Number of variables. (So if you think about a variable as a question in a survey, you can hypothetically think that you asked the machine for a gene, and the machine responded with a multinomial of size 1 by simply returning a single gene. Then you repeated/replicated this question many many times
J = 1

# Number of repetitions. This is now where we see exactly how many times the "machine responded"
Rj = c(sum(counts))

#Number of candidates (this is the number of categories in the multinomial data. so the number of options the machine could have returned for each replication)
Vj = dim(sig)[1]

#Distributions of each variable
dist = rep("multinomial",J)

dirichlet_parmeter = 0.1 #( 0 < x < 1 ), Where 0 is less certain about the priors and 1 is more certain about the priors.
#Inialization of populations proportions. 
alpha = rep(dirichlet_parmeter, K)
alpha[length(alpha)] = K*dirichlet_parmeter*8 #Mean of the initial dirichlet distribution of the last element (background) ~0.8 

#Start theta using the sig.
theta <- array(0, dim = c(J, K, max(Vj)))
for (j in 1:J) {
	theta[j, , ] = t(sig)/colSums(sig)
}

Nijr = array(1, dim = c(Total, J, max(Rj)))

#Reshape the observations
obs = array(-1,  dim = c(Total, J, max(Rj), max(Nijr)))
count_i=0;
count_j=0;
for (i in t(counts)) {
  for (j in seq(i)) { 
    if (i != 0 ) {
	    count_j = count_j + 1;
	    obs[1,1,count_j,1]=count_i;
	}
  } 
  count_i = count_i + 1; 
}

initial=mixedMemModel(Total = Total, J = J, Rj = Rj, Nijr = Nijr, K = K, Vj = Vj, alpha = alpha, theta = theta, dist = dist, obs = obs)
out = mmVarFit(initial, printStatus = 1, printMod = 25, stepType=0)

L1HS_hg38_reads = counts
corrected_L1HS_hg38_reads = t(out$phi)
rownames(corrected_L1HS_hg38_reads) = colnames(sig)
percentages = corrected_L1HS_hg38_reads/sum(corrected_L1HS_hg38_reads)

tot = scan(total_reads_file)
factor_tpm = scan(factor_tpm_file)
length = 6000

##
## Quantify element transcription
##
L1HS_hg38_reads = L1HS_hg38_reads[1:dim(L1HS_hg38_reads)[1]-1,,drop=FALSE]
corrected_L1HS_hg38_reads = corrected_L1HS_hg38_reads[1:length(corrected_L1HS_hg38_reads)-1,,drop=FALSE]

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

