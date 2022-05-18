set.seed(1)
m = 10000
n = 24
x = matrix(rnorm(m*n),m,n)
colnames(x)=1:n

d <- dist( t(x) )
hc<-hclust(d)
plot(hc)

plot(hc, cex = 0.25)

cl <- cutree(hc, h = 10)

m = 10000
n = 24
x = matrix(rnorm(m*n),m,n)
cl <- cutree(hc, h = 10)

library(rafalib)
set.seed(1)
m = 10000
n = 24
nc = replicate(100,{
  x = matrix(rnorm(m*n),m,n)
  hc = hclust( dist( t(x)))
  length(unique(cutree(hc,h=143)))
})
plot(table(nc)) ## look at the distribution
popsd(nc) # unbiased estimate of sd

km <-  kmeans(t(e), centers = 7)

table(tissue, clusters = km$cluster)

d <- dist(t(e))
mds <- cmdscale(d)
plot(mds[,1],mds[,2],col = km$cluster)

set.seed(10)
library(GSE5859Subset)
data(GSE5859Subset)


km<-kmeans(t(geneExpression),centers=5)
mds=cmdscale(dist(t(geneExpression)))
mypar(1,1)
plot(mds,bg=km$cluster,pch=21)
table(sampleInfo$group,km$cluster)
table(sampleInfo$date,km$cluster)
##looks better if we re-order:
table(sampleInfo$date,km$cluster)[,c(4,1,5,3,2)]

