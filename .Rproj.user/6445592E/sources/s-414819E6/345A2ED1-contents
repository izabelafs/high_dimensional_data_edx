library(tissuesGeneExpression)
data(tissuesGeneExpression)
image(e[1:100,])
library(genefilter)
rv <-  rowVars(e)
idx <- order(-rv)[1:40]

heatmap((e[idx,]))

library(RColorBrewer)
library(gplots)
hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)

heatmap.2(e[idx,],col = hmcol)

cols <- palette(brewer.pal(7, "Dark2"))[as.fumeric(tissue)]
cbind(colnames(e),tissue,cols)

heatmap.2(e[idx,], labCol = tissue, trace = "none", ColSideColors = cols, col = hmcol)

library(GSE5859Subset)
data(GSE5859Subset)
# Pick the 25 genes with the highest across sample variance. This function might help

install.packages("matrixStats")
library(matrixStats)
?rowMads ##we use mads due to a outlier sample

heatmap.2(e[idx,], labCol = tissue, trace = "none", ColSideColors = cols, col = hmcol)

sampleInfo$group

hmcol<-colorRampPalette(brewer.pal(9,"GnBu"))(100)
cols<-colorRampPalette(brewer.pal(8,"Dark2"))(2)[as.factor(sampleInfo$group)] # assigns a colour to each group
#cbind(colnames(geneExpression),cols) 
heatmap.2(geneExpression[idx,],
          labRow=geneAnnotation$CHR[idx],
          labCol=gsub("2005-","",sampleInfo$date),
          trace="none", # Rafa never uses this
          ColSideColors=cols,
          col=hmcol,
          scale="row")
cols = colorRampPalette(rev(brewer.pal(11,"RdBu")))(25)
gcol=brewer.pal(3,"Dark2")
gcol=gcol[sampleInfo$g+1]
##make lables: remove 2005 since it's common to all
labcol= gsub("2005-","",sampleInfo$date)  
##pick highly variable genes:
sds =rowMads(geneExpression)
ind = order(sds,decreasing=TRUE)[1:25]
## make heatmap
heatmap.2(geneExpression[ind,],
          col=cols,
          trace="none",
          scale="row",
          labRow=geneAnnotation$CHR[ind],
          labCol=labcol,
          ColSideColors=gcol,
          key=FALSE)

set.seed(17)
m = nrow(geneExpression)
n = ncol(geneExpression)
x = matrix(rnorm(m*n),m,n)
g = factor(sampleInfo$g)

sds =rowttests(x)
ind = order(sds,decreasing=FALSE)[1:50]
# sds[ind,3]
##make colors
cols = colorRampPalette(rev(brewer.pal(11,"RdBu")))(25)
gcol=brewer.pal(3,"Dark2")
gcol=gcol[sampleInfo$g+1]
## make heatmap
heatmap.2(geneExpression[ind,],
          col=cols,
          trace="none",
          scale="row",
          labRow=geneAnnotation$CHR[ind],
          labCol=labcol,
          ColSideColors=gcol,
          key=FALSE)

##load libraries
library(rafalib)
library(gplots)
library(matrixStats)
library(RColorBrewer)
##pick genes with smallest p-values:
sds =rowVars(x)
ind = order(sds,decreasing=TRUE)[1:50]
# sds[ind,3]
##make colors
cols = colorRampPalette(rev(brewer.pal(11,"RdBu")))(25)
gcol=brewer.pal(3,"Dark2")
gcol=gcol[sampleInfo$g+1]
## make heatmap
heatmap.2(geneExpression[ind,],
          col=cols,
          trace="none",
          scale="row",
          labRow=geneAnnotation$CHR[ind],
          labCol=labcol,
          ColSideColors=gcol,
          key=FALSE)

library(gplots)
library(matrixStats)
library(genefilter)
library(RColorBrewer)
cols = colorRampPalette(rev(brewer.pal(11,"RdBu")))(25)
ttest = rowttests(x,g)
sds = rowSds(x)
Indexes = list(t=order(ttest$p.value)[1:50], s=order(-sds)[1:50])
for(ind in Indexes){
  heatmap.2(x[ind,],
            col=cols,
            trace="none",
            scale="row",
            labCol=g,
            key=FALSE)
}
