library(GSE5859Subset)
data(GSE5859Subset)

s <-  svd(geneExpression)

s$d[1]
s$d[1]^2/(sum(s$d^2))


x = 1:ncol(geneExpression)
m <- sapply(x, function(x){
  md = mean(geneExpression[x,])
  })

m = rowMeans(geneExpression)
cor(m, s$u[,1])

y <- geneExpression - rowMeans(geneExpression)

s <- svd(y)

sy$d[1]
sy$d[1]^2/(sum(sy$d^2))

x = 1:ncol(sy$u)
m <- sapply(x, function(x){
  md = sy$d[x]^2/(sum(sy$d^2))
})
length(m[m > 0.05])

sum(m[1:10])

y <- geneExpression

s <-  svd(geneExpression)

y2 = s$u %*% diag(s$d) %*% t(s$v)
resid = y - y2
max(abs(resid))

z = s$d * t(s$v)

x <- geneExpression[,1]
Y <- geneExpression[,2]
sqrt(crossprod(x-Y))

x <- y[,1]
Y <- y[,2]
sqrt(crossprod(x-Y))

x <- z[,1]
Y <- z[,2]
sqrt(crossprod(x-Y))

x <- z[1:10,1]
Y <- z[1:10,2]
sqrt(crossprod(x-Y))

d = dist(t(geneExpression))
mds = cmdscale(d)

fdate = factor(sampleInfo$date)
plot(mds,col = fdate)
