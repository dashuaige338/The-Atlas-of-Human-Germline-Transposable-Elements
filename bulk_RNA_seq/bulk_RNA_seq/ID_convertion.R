rm(list=ls())
data = read.csv("20210713Xinyu_geneexpression_RPKM_onlyTEwithReads.csv",header =TRUE)
library(clusterProfiler)
library(org.Hs.eg.db)
list <- bitr(data$Geneid,
             fromType = 'ENSEMBL',toType = c('SYMBOL'),
             OrgDb = 'org.Hs.eg.db')
write.csv(list,file = "list.csv")
test = data
changedata = read.csv('./list.csv',header = TRUE, stringsAsFactors = FALSE)

ddd = changedata[match(test$Geneid, changedata$ENSEMBL),]
test$symbol = ddd$SYMBOL

bb = na.omit(test)
symbol = as.data.frame(bb$symbol)
bb <- cbind(symbol,bb)

cc = bb[!duplicated(bb$symbol),]
write.csv(cc,file = "gene.csv")
