
library(boot)
library(plyr)
library(ggplot2)

seedFromFile <- read.table('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/seed_marginalhist')
.Random.seed <- t(seedFromFile)

data <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/CRL.txt',header=T,stringsAsFactors=F)

numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');


for (col in numericcols){
  data[[col]] <- as.numeric(data[[col]])
}

dutyNoisefloorAbund <- ddply(data,.(gene),summarize,
                             dutycycle = mean(numTxnSites)/2,
                             noisefloor = 1/mean(cytoRNA)+.02,
                             meanRNA = mean(cytoRNA))

# This now agrees with Abhi's results!

boot.nm <- function(mydata,i){
  d <- mydata[i,]
  cytoRNA <- d$cytoRNA
  volume <- d$volume
  fit <- lm(cytoRNA ~ volume, data = d)
  slope <- fit$coeff[[2]]
  int <- fit$coeff[[1]]
  
  out <- (sd(cytoRNA)/mean(cytoRNA))^2 - 
    slope*mean(volume)/(int+slope*mean(volume)) * cov(cytoRNA,volume)/(mean(cytoRNA)*mean(volume))
  
  return(out)
}

boots <- ddply(data, .(gene), function(df) {
  results <- boot(data=df, statistic=boot.nm, R=5000)
  est <- mean(results$t)
  orderedResults <- results$t[order(results$t)]
  lower <- orderedResults[125]
  upper <- orderedResults[4875]
  out <- c(est,lower,upper)
  names(out) <- c('est','lower','upper');
  out
})

outData <- merge(dutyNoisefloorAbund,boots)
outData <- unique(outData)

cols = c('purple4','turquoise4','dodgerblue3')

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure6/Noise_Abund/noise_vs_abund.pdf',width=2.5,height=2.2)
ggplot(outData,aes(x=meanRNA,y=est,label=gene)) +
  geom_line(data=outData, aes(x=meanRNA, y=noisefloor), colour="gray78") +
  geom_line(data=outData, aes(x=meanRNA, y=1/meanRNA), col = 'gray78') +
  geom_point(size=1.2) +
  geom_errorbar(aes(ymin=lower,ymax=upper),width=.07) +
  scale_y_log10(limits = c(.01,4)) + scale_x_log10() +
  xlab('Mean RNA') + ylab('Noise Measure') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  theme(legend.position='')
dev.off()