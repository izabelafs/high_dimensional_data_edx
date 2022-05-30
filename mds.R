library(tissuesGeneExpression)
data(tissuesGeneExpression)

y = e - rowMeans(e)
s = svd(y)
z = s$d * t(s$v)

  library(rafalib)
ftissue = factor(tissue)
mypar(1,1)
plot(z[1,],z[2,],col=as.numeric(ftissue))
legend("topleft",levels(ftissue),col=seq_along(ftissue),pch=1)
d = dist(t(e))
mds = cmdscale(d)

cor(z[1,],mds[,1]) # Ans= 1 (-1)

cor(z[2,],mds[,2]) # Ans= 1 (-1)

library(rafalib)
ftissue = factor(tissue)
mypar(1,2)
plot(z[1,],z[2,],col=as.numeric(ftissue))
legend("topleft",levels(ftissue),col=seq_along(ftissue),pch=1)
plot(mds[,1],mds[,2],col=as.numeric(ftissue))

library(GSE5859Subset)
data(GSE5859Subset)

s = svd(geneExpression-rowMeans(geneExpression))
z = s$d * t(s$v)

ks = 1:nrow(z)
dimcors = sapply(ks,function(k){
  cor(z[k,],sampleInfo$group) 
})
plot(ks,dimcors) ##take a look
which.max( dimcors[] )

which.max(cor(sampleInfo$g,t(z))) # model answer

dimcors[which.max( dimcors[] )]

max(cor(sampleInfo$g,t(z))) # model answer

which.max(cor(sampleInfo$g,t(z))[-1]) + 1

sampleInfo$date

month = format( sampleInfo$date, "%m")
month = factor( month)

which.max(cor(as.numeric(month),t(z)))
max(cor(as.numeric(month),t(z)))

df<-data.frame(geneAnnotation$CHR,s$u[,6])
names(df)<-c("chr","ge")
df<-df[!df$chr == "chrUn", ]
boxplot(ge~chr,data=df)


result = split(s$u[,6],geneAnnotation$CHR)
result = result[ which(names(result)!="chrUn") ]
boxplot(result,range=0)
boxplot(result,range=0,ylim=c(-0.025,0.025))
medians = sapply(result,median)
names(result)[ which.max(abs(medians)) ]

