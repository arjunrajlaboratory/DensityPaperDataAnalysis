setwd('~/Dropbox/ExtractedData_131216/')
library(data.table)
library(ggplot2)
library(plyr)
library(boot)

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
  out <- c(est)
  names(out) <- c('est');
  out
})

abundDensTable <- ddply(dat,.(gene,cellType),summarize,
                        CVabund2 = (sd(cytoRNA)/mean(cytoRNA))^2,
                        CVdens2 = (sd(cytoRNA/volume)/mean(cytoRNA/volume))^2)

outTable <- merge(boots,abundDensTable)

meltTable <- melt(outTable,id.var=c('gene','cellType'))

pdf('~/Dropbox/densitypaper/SupplementaryFigures/Noise/CVabund_CVdens_noise_All.pdf',
    width=7.5,height=4)
ggplot(meltTable,aes(x=gene,y=value,fill=variable)) +
  geom_bar(stat='identity',position='dodge') +
  facet_wrap(~cellType, scales='free') +
  theme_classic() +
  xlab('Gene') + ylab('CV') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
dev.off()
