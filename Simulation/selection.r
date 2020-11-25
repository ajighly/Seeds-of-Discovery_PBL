library("data.table")
library("foreach")
library("doParallel")
library(Rcpp)
sourceCpp("test.cpp")
#rep
j=1

#PBLs with synthetic background
MYFULL=176
MYFULL1=22
MYFULL2=1665
Data=fread(paste("Tab_Short_Chr1_Rep-",j,"_GenerationNo-10000_out.txt",sep=''))#output of simulation software
Data=as.matrix(Data)
MYFULL3=matrix(0, nrow=(MYFULL+MYFULL1)*2, ncol=30000)
MYFULL5=sample(500,MYFULL)
MYFULL6=sample(500,MYFULL)
MYFULL7=sample(100,MYFULL1)
for (i in 1:MYFULL) {
	MYFULL3[i*2-1,1:10000]    = Data[1000+MYFULL6[i]*4-3,]
	MYFULL3[i*2,1:10000]      = Data[1000+MYFULL6[i]*4-2,]
	MYFULL3[i*2-1,10001:20000]= Data[1000+MYFULL6[i]*4-1,]
	MYFULL3[i*2,10001:20000]  = Data[1000+MYFULL6[i]*4  ,]
	MYFULL3[i*2-1,20001:30000]= Data[MYFULL5[i]*2-1,]
	MYFULL3[i*2,20001:30000]  = Data[MYFULL5[i]*2  ,]
}
for (i in 1:MYFULL1) {
	MYFULL3[MYFULL*2+i*2-1,1:10000]    =Data[3000+MYFULL7[i]*6-5,]
	MYFULL3[MYFULL*2+i*2,1:10000]      =Data[3000+MYFULL7[i]*6-4,]
	MYFULL3[MYFULL*2+i*2-1,10001:20000]=Data[3000+MYFULL7[i]*6-3,]
	MYFULL3[MYFULL*2+i*2,10001:20000]  =Data[3000+MYFULL7[i]*6-2,]
	MYFULL3[MYFULL*2+i*2-1,20001:30000]=Data[3000+MYFULL7[i]*6-1,]
	MYFULL3[MYFULL*2+i*2,20001:30000]  =Data[3000+MYFULL7[i]*6  ,]
}
ZMYFULL=cbind(sample(MYFULL,MYFULL2,replace=T),sample(MYFULL1,MYFULL2,replace=T)+MYFULL)
ZMYFULL1=MYFUNC(MYFULL3,3,10000,ZMYFULL[,1:2]-1)
ZMYFULL=cbind(1:MYFULL2,1:MYFULL2)
for (i in 1:6) ZMYFULL1=MYFUNC(ZMYFULL1,3,10000,ZMYFULL[,1:2]-1)
ZMYFULL1=rbind(ZMYFULL1,MYFULL3[(MYFULL*2+1):((MYFULL+MYFULL1)*2),])
ZMYFULL=cbind(1:MYFULL2,sample(MYFULL1,MYFULL2,replace=T)+MYFULL2)
ZMYFULL1=MYFUNC(ZMYFULL1,3,10000,ZMYFULL[,1:2]-1)
ZMYFULL=cbind(1:MYFULL2,1:MYFULL2)
for (i in 1:6) ZMYFULL1=MYFUNC(ZMYFULL1,3,10000,ZMYFULL[,1:2]-1)

MYFULL41XX <- matrix(nrow = nrow(ZMYFULL1)/2, ncol = ncol(ZMYFULL1), data = NA)
for (i in seq(2, nrow(ZMYFULL1), 2)) MYFULL41XX[(i/2),] <- colSums(ZMYFULL1[(i-2+1):i,])
MYFULL41XX=as.data.frame(MYFULL41XX)
MYFULL1XX <- sapply(MYFULL41XX ,function(x) sum(x)/(2*length(x)))
MYFULL4 <- as.matrix(MYFULL41XX[,MYFULL1XX>0 & MYFULL1XX<1])
MYFULL4 =t(MYFULL4)
MYFULL4[MYFULL4==0]="AA"
MYFULL4[MYFULL4==1]="AC"
MYFULL4[MYFULL4==2]="CC"
MYFULL1XX2=cbind(paste0('SNP',sprintf("%05d", 1:30000)),sort(rep(1:3,10000)),rep(seq(0.01,100,0.01)*1e6,3),rep(seq(0.01,100,0.01),3))
MYFULL1XX2=as.matrix(MYFULL1XX2[MYFULL1XX>0 & MYFULL1XX<1,])
MYFULL4=cbind(MYFULL1XX2,MYFULL4)
colnames(MYFULL4)=c('SNPname','chr','pos','genpos',paste0('gen',sprintf("%04d", 1:(ncol(MYFULL4)-4))))
fwrite(MYFULL4,paste("Syn_Rep_",j,".snp",sep=''),quote=F,row.names=F,col.names=T,sep="\t")



#PBLs with landrace background
MYFULL=190
MYFULL1=27
MYFULL2=1202
Data=fread(paste("Tab_Short_Chr1_Rep-",j,"_GenerationNo-10000_out.txt",sep=''))#output of PolySim
Data=as.matrix(Data)
MYFULL3=matrix(0, nrow=(MYFULL+MYFULL1)*2, ncol=30000)

Tri=matrix(0, nrow=200, ncol=30000)
for (i in 1:100) {
	ZMYFULL3[i*2-1,1:10000]    = Data[3000+i*6-5,]
	ZMYFULL3[i*2,1:10000]      = Data[3000+i*6-4,]
	ZMYFULL3[i*2-1,10001:20000]= Data[3000+i*6-3,]
	ZMYFULL3[i*2,10001:20000]  = Data[3000+i*6-2,]
	ZMYFULL3[i*2-1,20001:30000]= Data[3000+i*6-1,]
	ZMYFULL3[i*2,20001:30000]  = Data[3000+i*6  ,]
}
ZMYFULL=cbind(sample(100,MYFULL,replace=T),sample(100,MYFULL,replace=T))
ZMYFULL1=MYFUNC(MYFULL3,3,10000,ZMYFULL[,1:2]-1)
for (i in 1:50) {
	ZMYFULL=cbind(sample(MYFULL,MYFULL,replace=T),sample(MYFULL,MYFULL,replace=T))
	ZMYFULL1=MYFUNC(ZMYFULL1,3,10000,ZMYFULL[,1:2]-1)
}
ZMYFULL=cbind(1:MYFULL,1:MYFULL)
for (i in 1:6) ZMYFULL1=MYFUNC(ZMYFULL1,3,10000,ZMYFULL[,1:2]-1)
MYFULL3[1:(MYFULL*2),]=ZMYFULL1
MYFULL7=sample(100,MYFULL1)
for (i in 1:MYFULL1) {
	MYFULL3[MYFULL*2+i*2-1,1:10000]    =Data[3000+MYFULL7[i]*6-5,]
	MYFULL3[MYFULL*2+i*2,1:10000]      =Data[3000+MYFULL7[i]*6-4,]
	MYFULL3[MYFULL*2+i*2-1,10001:20000]=Data[3000+MYFULL7[i]*6-3,]
	MYFULL3[MYFULL*2+i*2,10001:20000]  =Data[3000+MYFULL7[i]*6-2,]
	MYFULL3[MYFULL*2+i*2-1,20001:30000]=Data[3000+MYFULL7[i]*6-1,]
	MYFULL3[MYFULL*2+i*2,20001:30000]  =Data[3000+MYFULL7[i]*6  ,]
}
ZMYFULL=cbind(sample(MYFULL,MYFULL2,replace=T),sample(MYFULL1,MYFULL2,replace=T)+MYFULL)
ZMYFULL1=MYFUNC(MYFULL3,3,10000,ZMYFULL[,1:2]-1)
ZMYFULL=cbind(1:MYFULL2,1:MYFULL2)
for (i in 1:6) ZMYFULL1=MYFUNC(ZMYFULL1,3,10000,ZMYFULL[,1:2]-1)
ZMYFULL1=rbind(ZMYFULL1,MYFULL3[(MYFULL*2+1):((MYFULL+MYFULL1)*2),])
ZMYFULL=cbind(1:MYFULL2,sample(MYFULL1,MYFULL2,replace=T)+MYFULL2)
ZMYFULL1=MYFUNC(ZMYFULL1,3,10000,ZMYFULL[,1:2]-1)
ZMYFULL=cbind(1:MYFULL2,1:MYFULL2)
for (i in 1:6) ZMYFULL1=MYFUNC(ZMYFULL1,3,10000,ZMYFULL[,1:2]-1)
MYFULL41XX <- matrix(nrow = nrow(ZMYFULL1)/2, ncol = ncol(ZMYFULL1), data = NA)
for (i in seq(2, nrow(ZMYFULL1), 2)) MYFULL41XX[(i/2),] <- colSums(ZMYFULL1[(i-2+1):i,])
MYFULL41XX=as.data.frame(MYFULL41XX)
MYFULL1XX <- sapply(MYFULL41XX ,function(x) sum(x)/(2*length(x)))
MYFULL4 <- as.matrix(MYFULL41XX[,MYFULL1XX>0 & MYFULL1XX<1])
MYFULL4 =t(MYFULL4)
MYFULL4[MYFULL4==0]="AA"
MYFULL4[MYFULL4==1]="AC"
MYFULL4[MYFULL4==2]="CC"

MYFULL1XX2=cbind(paste0('SNP',sprintf("%05d", 1:30000)),sort(rep(1:3,10000)),rep(seq(0.01,100,0.01)*1e6,3),rep(seq(0.01,100,0.01),3))
MYFULL1XX2=as.matrix(MYFULL1XX2[MYFULL1XX>0 & MYFULL1XX<1,])
MYFULL4=cbind(MYFULL1XX2,MYFULL4)
colnames(MYFULL4)=c('SNPname','chr','pos','genpos',paste0('gen',sprintf("%04d", 1:(ncol(MYFULL4)-4))))

fwrite(MYFULL4,paste("LR_Rep_",j,".snp",sep=''),quote=F,row.names=F,col.names=T,sep="\t")
