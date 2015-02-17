source('~/Dropbox/densitypaper/densitypaperdataanalysis/multiplot.R')
library(ggplot2)
library(plyr)

reg <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/CRL.txt',header=T,stringsAsFactors=F)
reg$volume <- reg$volume/1000
numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');


for (col in numericcols){
  reg[[col]] <- as.numeric(reg[[col]])
}
reg <- subset(reg,volume<8)

ss <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/SerumStarved.txt',header=T,stringsAsFactors=F)
ss <- ss[ss$numCyclin<20,]
ss$volume <- ss$volume/1000
numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');


for (col in numericcols){
  ss[[col]] <- as.numeric(ss[[col]])
}
ss <- subset(ss,volume<8)

seedFromFile <- read.table('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/seed_fig1_2')
.Random.seed <- t(seedFromFile)
regsub <- reg[runif(dim(reg)[1])<.08,]

cyclingCoeff <- lm(cytoGAPDH ~ volume, data=regsub)$coefficients

tmpCoeff <- lm(cytoGAPDH ~ volume, data=reg)$coefficients

seedFromFile <- read.table('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/seed_ss')
.Random.seed <- t(seedFromFile)
regsub <- reg[runif(dim(reg)[1])<.08,]

sssub <- ss[runif(dim(ss)[1])<.1,]

dat <- rbind(regsub,sssub)

maxVol <- max(dat$volume)
maxRna <- max(dat$cytoGAPDH)

colors <- c('lightseagreen','aquamarine3','seagreen3','gray87')

# CYCLING
p1 <- ggplot(regsub, aes(x=volume,y=cytoGAPDH)) +
  geom_abline(slope=cyclingCoeff[[2]],intercept=cyclingCoeff[[1]],linetype='dashed',col=colors[4]) +
  geom_point(size=1,col=colors[1]) +
  scale_x_continuous(limits = c(0, maxVol)) + scale_y_continuous(limits = c(0, maxRna)) +
  xlab('') + ylab('Cytoplasmic RNA') +
  theme_bw() + theme_classic() + theme(legend.position = "none") + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  theme(axis.text.y = element_text(angle = 90, hjust = 1))

# QUIESCENT
p2 <- ggplot(sssub, aes(x=volume,y=cytoGAPDH)) +
  geom_abline(slope=cyclingCoeff[[2]],intercept=cyclingCoeff[[1]],linetype='dashed',col=colors[4]) +
  geom_point(size=1,col=colors[2]) +
  scale_x_continuous(limits = c(0, maxVol)) + scale_y_continuous(limits = c(0, maxRna)) +
  xlab('Volume (picoliter)') + ylab('') +
  theme_bw() + theme_classic() + theme(legend.position = "none") + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  theme(axis.text.y = element_text(angle = 90, hjust = 1))



tmp <- ddply(dat,.(cellType),summarize,
             meanRNA = mean(cytoGAPDH),
             sterrRNA = sd(cytoGAPDH)/sqrt(length(cytoGAPDH)),
             sdRNA = sd(cytoGAPDH),
             meanDensity = mean(cytoGAPDH/volume),
             sterrDensity = sd(cytoGAPDH/volume)/sqrt(length(cytoGAPDH)),
             sdDensity = sd(cytoGAPDH/volume))

tmp$cellType <- c('Cycling','Quiescent')

m1 <- ggplot(tmp,aes(x=cellType,y=meanRNA)) +
  geom_bar(stat='identity',fill=colors[1:2]) +
  expand_limits(y=4850) +
  xlab('Cell Type') + ylab('Mean RNA') +
  geom_errorbar(aes(ymin=meanRNA-sterrRNA,ymax=meanRNA+sterrRNA),width=.25) +
  theme_bw() + theme_classic() + theme(legend.position = "none") + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

m2 <- ggplot(tmp,aes(x=cellType,y=meanDensity)) +
  geom_bar(stat='identity',fill=colors[1:2]) +
  expand_limits(y=2500) +
  geom_errorbar(aes(ymin=meanDensity-sterrDensity,ymax=meanDensity+sterrDensity),width=.25) +
  xlab('Cell Type') + ylab('Mean [RNA]') +
  theme_bw() + theme_classic() + theme(legend.position = "none") + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

pdf(file='~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure1/Cyc_Sen_SS/Cyc_Sen_SS.pdf',height=1.7,width=9)
multiplot(p1,p2,m1,m2,cols=4)
dev.off()
