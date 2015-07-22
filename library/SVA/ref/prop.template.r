options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

#
input = args[1]
output = args[2]
#
temptable=read.table(file=input,row.names=1,header=T);
prop.table(temptable);
write.table(file=output,x=prop.table(temptable),quote=F,append=T);
