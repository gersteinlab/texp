options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

#
input = args[1]
output = args[2]
#
temptable=read.table(file=input,row.names=1,header=T);
t=t(t(temptable)/colSums(temptable))
#t=temptable/rowSums(temptable)
write.table(file=output,x=t,quote=F,append=T);
