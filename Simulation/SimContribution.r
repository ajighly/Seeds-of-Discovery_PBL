library("data.table")
library("foreach")
library("doParallel")
library(Rcpp)
sourceCpp("test.cpp")


MYFULL9=matrix(0, nrow=700, ncol=1)
MYFULL=34
MYFULL1=17
MYFULL2=125

for (j in 1:700) {
	Data=fread(pMYFULL6e("../Selection/Tab_Short_Chr1_Rep-",j,"_GenerationNo-10000_out.txt",sep=''))
	Data=as.matrix(Data)
	MYFULL3=matrix(0, nrow=(MYFULL+MYFULL1)*2, ncol=30000)
	MYFULL4=sample(500,MYFULL)
	MYFULL5=sample(500,MYFULL)
	MYFULL6=sample(100,MYFULL1)
	for (i in 1:MYFULL) {
		MYFULL3[i*2-1,1:10000]    = Data[1000+MYFULL5[i]*4-3,]
		MYFULL3[i*2,1:10000]      = Data[1000+MYFULL5[i]*4-2,]
		MYFULL3[i*2-1,10001:20000]= Data[1000+MYFULL5[i]*4-1,]
		MYFULL3[i*2,10001:20000]  = Data[1000+MYFULL5[i]*4  ,]
		MYFULL3[i*2-1,20001:30000]= Data[MYFULL4[i]*2-1,]
		MYFULL3[i*2,20001:30000]  = Data[MYFULL4[i]*2  ,]
	}
	for (i in 1:MYFULL1) {
		MYFULL3[MYFULL*2+i*2-1,1:10000]    =Data[3000+MYFULL6[i]*6-5,]
		MYFULL3[MYFULL*2+i*2,1:10000]      =Data[3000+MYFULL6[i]*6-4,]
		MYFULL3[MYFULL*2+i*2-1,10001:20000]=Data[3000+MYFULL6[i]*6-3,]
		MYFULL3[MYFULL*2+i*2,10001:20000]  =Data[3000+MYFULL6[i]*6-2,]
		MYFULL3[MYFULL*2+i*2-1,20001:30000]=Data[3000+MYFULL6[i]*6-1,]
		MYFULL3[MYFULL*2+i*2,20001:30000]  =Data[3000+MYFULL6[i]*6  ,]
	}
	MYFULL8=cbind(sample(MYFULL,MYFULL2,replace=T),sample(MYFULL1,MYFULL2,replace=T)+MYFULL)
	MYFULL32=matrix(0, nrow=(MYFULL+MYFULL1)*2, ncol=30000)
	MYFULL32[(MYFULL*2+1):((MYFULL+MYFULL1)*2),]=2
	MYFULL7=MYFUNC(MYFULL32,3,10000,MYFULL8[,1:2]-1)
	MYFULL8=cbind(1:MYFULL2,1:MYFULL2)
	for (i in 1:6) MYFULL7=MYFUNC(MYFULL7,3,10000,MYFULL8[,1:2]-1)
	MYFULL7=rbind(MYFULL7,MYFULL32[(MYFULL*2+1):((MYFULL+MYFULL1)*2),])
	MYFULL8=cbind(1:MYFULL2,sample(MYFULL1,MYFULL2,replace=T)+MYFULL2)
	MYFULL7=MYFUNC(MYFULL7,3,10000,MYFULL8[,1:2]-1)
	MYFULL8=cbind(1:MYFULL2,1:MYFULL2)
	for (i in 1:6) MYFULL7=MYFUNC(MYFULL7,3,10000,MYFULL8[,1:2]-1)
	MYFULL9[j,1]=length(which(MYFULL7==0))/7500000
}
write.table(MYFULL9,'out.txt',quote=F)
