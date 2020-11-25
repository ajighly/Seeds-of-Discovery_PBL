library(data.table)
FullData=fread('Geno.txt',data.table=F)
map=FullData[,1:3]
row.names(FullData)=FullData[,1]
FullData=FullData[,-c(1:3)]
keep=read.table('keep',header=F)#subset of PBLs need to be analyzed
info=read.table('info',header=F)#pedigrees formated as the first column: the ID of the PBL, second column: the ID of the first parent, third column: the ID of the second parent
synthetics=FullData[,keep[,1]]
Contrib=matrix(NA,nrow=ncol(synthetics),ncol=2)
for (i in 1:ncol(synthetics)) {
  Geno=cbind(1:nrow(FullData),as.character(FullData[,info[i,1]]),as.character(FullData[,info[i,2]]),as.character(FullData[,info[i,3]]),as.character(synthetics[,i]))
  Geno=na.omit(Geno)
  Geno <- subset(Geno, Geno[,4] == Geno[,3] & Geno[,2] != Geno[,4])
  for (j in 2:5) Geno <- subset(Geno, Geno[,j] != "NN" )
  Contrib[i,1]=length(which(Geno[,2]==Geno[,5]))
  Contrib[i,2]=length(which(Geno[,4]==Geno[,5]))
}
row.names(Contrib)=colnames(synthetics)
write.table(Contrib,"contributions.txt",col.names=F,quote = F,sep = "\t")
