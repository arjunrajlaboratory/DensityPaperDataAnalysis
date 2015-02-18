library(boot)
library(plyr)
library(ggplot2)

seedFromFile <- read.table('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/seed_marginalhist')
.Random.seed <- t(seedFromFile)

crl <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/CRL.txt',header=T,stringsAsFactors=F)
a549 <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/A549.txt',header=T,stringsAsFactors=F)
quies <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/SerumStarved.txt',header=T,stringsAsFactors=F)

data <- rbind(crl,a549,quies)

numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');


for (col in numericcols){
  data[[col]] <- as.numeric(data[[col]])
}

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

boots <- ddply(data, .(gene,cellType), function(df) {
  results <- boot(data=df, statistic=boot.nm, R=5000)
  est <- mean(results$t)
  orderedResults <- results$t[order(results$t)]
  lower <- orderedResults[125]
  upper <- orderedResults[4875]
  mean <- mean(df$cytoRNA)
  sd <- sd(df$cytoRNA)
  out <- c(est,lower,upper,mean,sd)
  names(out) <- c('est','lower','upper','mean','sd');
  out
})

crl <- subset(boots,cellType %in% 'CRL2097')[,c('gene','est','lower','upper','mean','sd')]
quies <- subset(boots,cellType %in% 'CRL2097_SerumStarved')[,c('gene','est','lower','upper','mean','sd')]
a549 <- subset(boots,cellType %in% 'A549')[,c('gene','est','lower','upper','mean','sd')]

colnames(crl) <- c(colnames(crl)[1],paste(colnames(crl)[2:6],'crl',sep='.'))
colnames(quies) <- c(colnames(quies)[1],paste(colnames(quies)[2:6],'quies',sep='.'))
colnames(a549) <- c(colnames(a549)[1],paste(colnames(a549)[2:6],'a549',sep='.'))

dat2 <- merge(crl,quies)
dat3 <- merge(crl,a549)

### Now plot

cols = c('purple4','darkorchid3','magenta4')

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/SingleCellSeq_LowHi_SpecNon/Noise_CellTypes/cyc_vs_A549.pdf',width=2.2,height=2)
ggplot(dat3, aes(x=est.crl,y=est.a549,label=gene)) +
  #geom_text() +
  geom_abline(slope=1,intercept=0,col='gray78') +
  geom_point(size=1,col=cols[1]) +
  scale_y_log10() + scale_x_log10() +
  geom_errorbarh(aes(xmin=lower.crl,xmax=upper.crl),height=.1,col=cols[1]) +
  geom_errorbar(aes(ymin=lower.a549,ymax=upper.a549),width=.1,col=cols[1]) +
  expand_limits(x=c(.01,20),y=c(.01,20)) +
  xlab('Noise measure (CRL)') + ylab('Noise Measure (A549)') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()

dat3$meanlower.crl <- dat3$mean.crl-dat3$sd.crl
dat3$meanupper.crl <- dat3$mean.crl+dat3$sd.crl
dat3$meanlower.crl[dat3$meanlower.crl<0] <- dat3$mean.crl[dat3$meanlower.crl<0]

dat3$meanlower.a549 <- dat3$mean.a549-dat3$sd.a549
dat3$meanupper.a549 <- dat3$mean.a549+dat3$sd.a549
dat3$meanlower.a549[dat3$meanlower.a549<0] <- dat3$mean.a549[dat3$meanlower.a549<0]

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/SingleCellSeq_LowHi_SpecNon/Noise_CellTypes/cyc_vs_A549_abund.pdf',width=2.2,height=2)
ggplot(dat3, aes(x=mean.crl,y=mean.a549,label=gene)) +
  #geom_text() +
  geom_abline(slope=1,intercept=0,col='gray78') +
  geom_point(size=1,col=cols[1]) +
  scale_y_log10() + scale_x_log10() +
  geom_errorbarh(aes(xmin=meanlower.crl,xmax=meanupper.crl),height=.1,col=cols[1]) +
  geom_errorbar(aes(ymin=meanlower.a549,ymax=meanupper.a549),width=.1,col=cols[1]) +
  #expand_limits(x=c(1,4000),y=c(1,4000)) +
  xlab('Mean (CRL)') + ylab('Mean (A549)') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()

dat2$meanlower.crl <- dat2$mean.crl-dat2$sd.crl
dat2$meanupper.crl <- dat2$mean.crl+dat2$sd.crl
dat2$meanlower.crl[dat2$meanlower.crl<0] <- dat2$mean.crl[dat2$meanlower.crl<0]

dat2$meanlower.quies <- dat2$mean.quies-dat2$sd.quies
dat2$meanupper.quies <- dat2$mean.quies+dat2$sd.quies
dat2$meanlower.quies[dat2$meanlower.quies<0] <- dat2$mean.quies[dat2$meanlower.quies<0]
