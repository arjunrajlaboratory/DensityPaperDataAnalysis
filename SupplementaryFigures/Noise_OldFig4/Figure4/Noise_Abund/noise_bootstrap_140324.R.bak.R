library(boot)
library(plyr)
library(ggplot2)

seedFromFile <- read.table('~/Dropbox/densitypaper/densityfigures/seed_marginalhist')
.Random.seed <- t(seedFromFile)

data <- read.delim('~/Dropbox/ExtractedData_131216/alldata_CRL.txt',header=T,stringsAsFactors=F)

numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');


for (col in numericcols){
  data[[col]] <- as.numeric(data[[col]])
}

hl <- read.delim('~/Dropbox/Foreskin/half_life_with_gene_names.txt',header=F,stringsAsFactors=F)

colnames(hl) <- c('nm','gene','hl')
hl <- hl[,c('gene','hl')]
hl <- hl[hl$hl!='N.D.',]
hl$hl <- gsub(' ', '', hl$hl)
hl$hl <- gsub('>24', '24', hl$hl)
hl$hl <- as.numeric(hl$hl)

#data <- data[data$numCyclin<20 & data$numTxnSites<3,]
#data <- ddply(data,.(gene),subset,length(cytoRNA)>20)

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

outData <- merge(dutyNoisefloorAbund, hl)
outData <- merge(outData,boots)
outData <- unique(outData)

cols = c('purple4','turquoise4','dodgerblue3')

#pdf('~/Dropbox/densitypaper/densityfigures/Figure4/noise_vs_abund.pdf',width=2.5,height=2.2)
ggplot(outData,aes(x=meanRNA,y=est,label=gene)) +
  #geom_text() +
  geom_line(data=outData, aes(x=meanRNA, y=noisefloor), colour="gray78") +
  geom_line(data=outData, aes(x=meanRNA, y=1/meanRNA), col = 'gray78') +
  geom_point(size=1.2) +
  geom_errorbar(aes(ymin=lower,ymax=upper),width=.07) +
  scale_y_log10(limits = c(.01,4)) + scale_x_log10() +
  xlab('Mean RNA') + ylab('Noise Measure') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  theme(legend.position='')
#dev.off()