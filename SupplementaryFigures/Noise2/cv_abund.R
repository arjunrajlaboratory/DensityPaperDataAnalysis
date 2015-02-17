library(ggplot2)
library(plyr)

data <- read.delim('~/Dropbox/densitypaper/ExtractedData/CRL.txt',header=T,stringsAsFactors=F)

numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');


for (col in numericcols){
  data[[col]] <- as.numeric(data[[col]])
}

data$volume <- data$volume/1000

cvTable <- ddply(data,.(gene),summarize,
                 cv = sd(cytoRNA)/mean(cytoRNA),
                 mean = mean(cytoRNA),
                 seMean = sd(cytoRNA)/sqrt(length(cytoRNA)))

boot.cv <- function(mydata,i){
  d <- mydata[i,]
  cytoRNA <- d$cytoRNA
  volume <- d$volume
  
  out <- (sd(cytoRNA)/mean(cytoRNA))
  
  return(out)
}

cvR2Err <- ddply(data, .(gene), function(df) {
  results <- boot(data=df, statistic=boot.cv, R=5000)
  orderedResults <- results$t[order(results$t)]
  meas <- results$t0
  lower <- orderedResults[125]
  upper <- orderedResults[4875]
  
  out <- c(meas,lower, upper)
  names(out) <- c('value','lower','upper');
  out
})

outDat <- merge(cvTable,cvR2Err)

pdf('~/Dropbox/densitypaper/SupplementaryFigures/Noise2/cv_abund.pdf',width=3.3,height=3)
ggplot(outDat,aes(x=mean,y=cv)) +
  geom_point(size=1.3) +
  theme_classic() +
  geom_errorbar(aes(ymin=lower,ymax=upper),width=50) +
  ylab('mRNA CV') + xlab('Average mRNA count') +
  theme(axis.title=element_text(size='6'), axis.text=element_text(size='6'))
dev.off()
