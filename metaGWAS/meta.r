#Provided by Reem Joukhadar
Ts=read.table('Ts.txt')
n=ncol(Ts)
pvalues=matrix(NA,nrow=nrow(Ts),ncol=1)
V=solve(cor(Ts))
for (i in 1:nrow(Ts)) pvalues[i,1]=pchisq(as.matrix(t(Ts[i,]))%*%V%*%as.matrix(Ts[i,]), df=n, lower.tail=FALSE)
write.table(pvalues,'metaAll',quote=F,sep='\t')
