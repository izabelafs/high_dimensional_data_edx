library(devtools)
install_github("genomicsclass/tissuesGeneExpression")
install_github("genomicsclass/GSE5859Subset")

# The data represents RNA expression levels for seven tissues, each with several biological replicates. We call samples that we consider to be from the same population, such as liver tissue from different individuals, biological replicates:
library(tibble)
library(tissuesGeneExpression)

data(tissuesGeneExpression)
head(e)
head(tissue)
dist <- dist(t(e))

dist[3,45]

x<-e["210486_at",]
y<-e["200805_at",]
sqrt(crossprod(x-y))
d = dist(t(e))
length(d)

library(GSE5859Subset)
data(GSE5859Subset)
dim(geneExpression)
table(sampleInfo)
summary(sampleInfo)
sampleInfo
x <- geneExpression[3,]
y <- geneExpression[7,]
sqrt(crossprod(x-y))
x <- geneExpression[4,]
y <- geneExpression[14,]
sqrt(crossprod(x-y))
d = dist(t(geneExpression))
as.matrix(d)[3,7]
as.matrix(d)[4,14]


column = 1
y = nrow(geneExpression):1
x = 1:ncol(geneExpression)
dists = sapply(x, function(x){
  test = geneExpression[,x]
  sapply(y, function(y){
  target = geneExpression[y,]
  sqrt(crossprod(target-test))
  })
})
mean(dists)
z = 1:24
dist = sapply(z, function(z){
  d = mean(dists[,z])
})
which.max(dist)

d.ge <- dist(t(geneExpression))

x<-e["201371_s_at",]
y<-e["1007_s_at",]
sqrt(crossprod(x-y))

e <- as.matrix(geneExpression)

x1 <- e["202138_x_at",]
y1 <- e["202152_x_at",]
sqrt(crossprod(x1-y1))
ge <- dist(geneExpression)

