library(data.table)
library(ggplot2)
library(plyr)

crlDat <- read.delim('~/Dropbox/densitypaper/ExtractedData/CRL.txt',header=T,stringsAsFactors=F)
crlDat$type <- 'cycling'

ssDat <- read.delim('~/Dropbox/densitypaper/ExtractedData/SerumStarved.txt',header=T,stringsAsFactors=F)
ssDat$type <- 'quiescent'

senDat <- read.delim('~/Dropbox/densitypaper/ExtractedData/Senescent.txt',header=T,stringsAsFactors=F)
senDat <- subset(senDat,!(date %in% 140128))
senDat$type <- 'senescent'

a549Dat <- read.delim('~/Dropbox/densitypaper/ExtractedData/A549.txt',header=T,stringsAsFactors=F)
a549Dat$type <- 'a549'

numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');

dat <- rbind(crlDat,ssDat,senDat,a549Dat)

for (col in numericcols){
  dat[[col]] <- as.numeric(dat[[col]])
}

dat$volume <- dat$volume/1000

densDat <- ddply(dat,c('gene','type'),summarize,
                 meanDens = mean(cytoRNA/volume))

dat <- merge(dat,densDat,by=c('gene','type')) 

crlDat <- subset(dat,type %in% 'cycling')
ssDat <- subset(dat,type %in% 'quiescent')
senDat <- subset(dat,type %in% 'senescent')
a549Dat <- subset(dat,type %in% 'a549')

pdf('~/Dropbox/densitypaper/SupplementaryFigures/DataDump/CRL_cycling_Conc.pdf',width=7.5,height=4)
ggplot(crlDat, 
       aes(x=volume,y=cytoRNA/volume)) + 
  geom_point(size=1) + facet_wrap(~gene, scales = 'free_y') + expand_limits(x=0,y=0) + 
  theme_classic() +
  geom_hline(aes(yintercept=meanDens),col='red') +
  theme(axis.title=element_text(size="10"), axis.text=element_text(size='10')) +
  xlab('Volume (picoliter)') + ylab('mRNA Concentration (count/pL)') + 
  theme(strip.background = element_rect(fill='white'),strip.text = element_text(face='italic',size=10))
dev.off()

pdf('~/Dropbox/densitypaper/SupplementaryFigures/DataDump/CRL_quiescent_Conc.pdf',width=7.5,height=4)
ggplot(ssDat, 
       aes(x=volume,y=cytoRNA/volume)) + 
  geom_point(size=1) + facet_wrap(~gene, scales = 'free_y') + expand_limits(x=0,y=0) + 
  theme_classic() +
  geom_hline(aes(yintercept=meanDens),col='red') +
  theme(axis.title=element_text(size="10"), axis.text=element_text(size='10')) +
  xlab('Volume (picoliter)') + ylab('mRNA Concentration (count/pL)') + 
  theme(strip.background = element_rect(fill='white'),strip.text = element_text(face='italic',size=10))
dev.off()

pdf('~/Dropbox/densitypaper/SupplementaryFigures/DataDump/CRL_senescent_Conc.pdf',width=7.5,height=4)
ggplot(senDat, 
       aes(x=volume,y=cytoRNA/volume)) + 
  geom_point(size=1.5) + facet_wrap(~gene, scales = 'free_y') + expand_limits(x=0,y=0) + 
  theme_classic() +
  geom_hline(aes(yintercept=meanDens),col='red') +
  theme(axis.title=element_text(size="10"), axis.text=element_text(size='10')) +
  xlab('Volume (picoliter)') + ylab('mRNA Concentration (count/pL)') + 
  theme(strip.background = element_rect(fill='white'),strip.text = element_text(face='italic',size=10))
dev.off()

pdf('~/Dropbox/densitypaper/SupplementaryFigures/DataDump/A549_cycling_Conc.pdf',width=7.5,height=4)
ggplot(a549Dat, aes(x=volume,y=cytoRNA/volume)) + 
  geom_point(size=1) + facet_wrap(~gene, scales = 'free_y') + expand_limits(x=0,y=0) + 
  theme_classic() +
  geom_hline(aes(yintercept=meanDens),col='red') +
  theme(axis.title=element_text(size="10"), axis.text=element_text(size='10')) +
  xlab('Volume (picoliter)') + ylab('mRNA Concentration (count/pL)') + 
  theme(strip.background = element_rect(fill='white'),strip.text = element_text(face='italic',size=10))
dev.off()