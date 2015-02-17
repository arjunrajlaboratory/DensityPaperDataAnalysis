library(boot)
library(plyr)
library(ggplot2)

seedFromFile <- read.table('~/Dropbox/densitypaper/densityfigures/seed_marginalhist')
.Random.seed <- t(seedFromFile)

crlDat <- read.delim('alldata_CRL.txt',header=T,stringsAsFactors=F)
crlDat$type <- 'cycling'

ssDat <- read.delim('alldata_SS.txt',header=T,stringsAsFactors=F)
ssDat$type <- 'quiescent'

senDat <- read.delim('alldata_Senescent.txt',header=T,stringsAsFactors=F)
senDat$type <- 'senescent'

a549Dat <- read.delim('alldata_A549.txt',header=T,stringsAsFactors=F)
a549Dat$type <- 'a549'

dat <- rbind(crlDat,ssDat,senDat,a549Dat)

numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');


for (col in numericcols){
  dat[[col]] <- as.numeric(dat[[col]])
}

dat$volume <- dat$volume/1000

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

boots <- ddply(dat, .(gene,cellType), function(df) {
  results <- boot(data=df, statistic=boot.nm, R=5000)
  est <- mean(results$t)
  orderedResults <- results$t[order(results$t)]
  lower <- orderedResults[125]
  upper <- orderedResults[4875]
  out <- c(est,lower,upper)
  names(out) <- c('est','lower','upper');
  out
})

crl <- subset(boots, cellType %in% 'CRL2097')
ss <- subset(boots, cellType %in% 'CRL2097_SerumStarved')
sen <- subset(boots, cellType %in% 'CRL2097_Senescent')
a549 <- subset(boots, cellType %in% 'A549')

cols <- c('deepskyblue4','darkorchid4','deeppink4','green4')

pdf('~/Dropbox/densitypaper/SupplementaryFigures/Noise/CRL.pdf',width=4,height=1.5)
ggplot(crl,aes(x=gene,y=est)) +
  geom_bar(stat='identity',fill=cols[1]) +
  #geom_errorbar(aes(ymin=lower,ymax=upper)) +
  theme_classic() +
  xlab('Gene') + ylab('Noise Measure') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf('~/Dropbox/densitypaper/SupplementaryFigures/Noise/Quiescent.pdf',width=4,height=1.5)
ggplot(ss,aes(x=gene,y=est)) +
  geom_bar(stat='identity',fill=cols[2]) +
  #geom_errorbar(aes(ymin=lower,ymax=upper)) +
  theme_classic() +
  xlab('Gene') + ylab('Noise Measure') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf('~/Dropbox/densitypaper/SupplementaryFigures/Noise/Senescent.pdf',width=4,height=1.5)
ggplot(sen,aes(x=gene,y=est)) +
  geom_bar(stat='identity',fill=cols[3]) +
  #geom_errorbar(aes(ymin=lower,ymax=upper)) +
  theme_classic() +
  xlab('Gene') + ylab('Noise Measure') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf('~/Dropbox/densitypaper/SupplementaryFigures/Noise/A549.pdf',width=4,height=1.5)
ggplot(a549,aes(x=gene,y=est)) +
  geom_bar(stat='identity',fill=cols[4]) +
  #geom_errorbar(aes(ymin=lower,ymax=upper)) +
  theme_classic() +
  xlab('Gene') + ylab('Noise Measure') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()