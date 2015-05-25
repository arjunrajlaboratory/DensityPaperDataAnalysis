setwd('~/Dropbox/Code/R/TxnSiteAnalysis/')
library(reshape)

alldata <- read.delim('alldata_crl.txt',header=T,stringsAsFactors=F)
alldata$id <- paste(alldata$cellType, alldata$gene, alldata$date, alldata$cellNumber, sep="_")

### cast way of doing things
tmp <- cast(alldata, id ~ variable)
for(i in which(!colnames(tmp) %in% c("id", "cellCycleStage"))) {
  tmp[,i] <- as.vector(tmp[,i])
  tmp[,i] <- as.numeric(tmp[,i])
}

### tapply way of doing things
#tmp <- tapply(alldata$value, list(alldata$id, alldata$variable), function(x) x)

idData <- do.call(rbind, strsplit(tmp$id, "_"))
idData <- as.data.frame(idData, stringsAsFactors=F)
colnames(idData) <- colnames(alldata)[1:4]
finalData <- cbind(idData, tmp)
finalData$normNumTxnSites <- finalData$numTxnSites
finalData$normNumTxnSites[finalData$cellCycleStage == 'G2'] <- 
  finalData$normNumTxnSites[finalData$cellCycleStage == 'G2']/2

geneData <- finalData[finalData$gene == "ACTA2",]
boxplot(geneData$numRnaPerTxnSiteAvgExon ~ geneData$cellCycleStage, na.rm=T)

pdf("trial.pdf")
for(gene in unique(finalData$gene)) {
  geneData <- finalData[finalData$gene == gene,]
  if(!all(is.na(geneData$volume))) {
    par(mar=c(4.5,5.1,0.5,0.5), oma=c(1,.5,2,2),mfrow=c(2,2))
    boxplot(geneData$numRnaPerTxnSiteAvgExon ~ geneData$cellCycleStage, na.rm=T,ylab='intensity')
    plot(geneData$volume,geneData$numRnaPerTxnSiteAvgExon,
     col = as.factor(geneData$cellCycleStage),ylab='intensity',xlab='vol')

    boxplot(geneData$numTxnSites ~ geneData$cellCycleStage, na.rm=T,ylab='num txn sites')
    boxplot(geneData$volume~geneData$numTxnSites,na.rm=T,xlab='num txn sites',ylab='vol')
    mtext(gene, outer = TRUE, cex = 1.5)
  }
}
dev.off()
