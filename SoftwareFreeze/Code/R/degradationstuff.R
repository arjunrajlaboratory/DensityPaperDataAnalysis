setwd('~/Dropbox/ExtractedData_131216/')
library(data.table)
library(ggplot2)
library(plyr)

crlDat <- read.delim('alldata_CRL.txt',header=T,stringsAsFactors=F)
crlDat$type <- 'cycling'

ssDat <- read.delim('alldata_SS.txt',header=T,stringsAsFactors=F)
ssDat$type <- 'quiescent'

senDat <- read.delim('alldata_Senescent.txt',header=T,stringsAsFactors=F)
senDat$type <- 'senescent'

a549Dat <- read.delim('alldata_A549.txt',header=T,stringsAsFactors=F)
a549Dat$type <- 'a549'

halflife <- read.delim('~/Dropbox/Foreskin/half_life_rpkm_with_names.txt',header=T,stringsAsFactors=F)
hl <- halflife[,c('gene_name','halflife')]
colnames(hl) <- c('gene','halflife')
hl <- subset(hl,halflife!='N.D.')
hl$halflife[hl$halflife=='>24'] = 24
hl$halflife <- as.numeric(hl$halflife)

dat <- rbind(crlDat,ssDat,senDat,a549Dat)
dat <- merge(dat,hl)

numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');


for (col in numericcols){
  crlDat[[col]] <- as.numeric(crlDat[[col]])
  ssDat[[col]] <- as.numeric(ssDat[[col]])
  senDat[[col]] <- as.numeric(senDat[[col]])
  a549Dat[[col]] <- as.numeric(a549Dat[[col]])
}

crlDat$repNum <- NA

for(i in 1:length(unique(crlDat$gene))) {
  geneiter <- as.character(unique(crlDat$gene))[i]
  tmp <- subset(crlDat,gene %in% geneiter)
  for(j in 1:length(unique(tmp$date))){
    dateiter <- as.character(unique(tmp$date))[j]
    crlDat$repNum[(crlDat$gene==geneiter & crlDat$date ==dateiter)] <- j
  }
}

dat$mu <- dat$cytoRNA*log(2)/dat$halflife

lmDat <- ddply(dat,c('gene','type'),function(x) {
  tmp <- summary(lm(mu~volume,x))$coefficients; 
  out <- c(tmp[1,], tmp[2,]);
  names(out) <- paste(rep(c("inter", "slope"), each=ncol(tmp)), 
                      c('estimate','stderr','tval','pr'), sep=".");
  out
})

meanDat <- ddply(dat,~type,summarise,mean=mean(volume))

dat <- merge(dat,lmDat,by=c('gene','type')) 
dat <- merge(dat,meanDat)

dat$volindep.val <- dat$inter.estimate
dat$volindep.err <- dat$inter.stderr
dat$voldep.val <- dat$slope.estimate*dat$mean
dat$voldep.err <- dat$slope.stderr*dat$mean

indepDat <- dat[,c('type','gene','volindep.val','volindep.err')]
indepDat <- unique(indepDat)
colnames(indepDat) <- c('type','gene','val','err')
indepDat$voldep <- 'indep'
indepDat$graphTitle <- paste(indepDat$type,indepDat$voldep, sep=".")

depDat <- dat[,c('type','gene','voldep.val','voldep.err')]
depDat <- unique(depDat)
colnames(depDat) <- c('type','gene','val','err')
depDat$voldep <- 'dep'
depDat$graphTitle <- paste(depDat$type,depDat$voldep, sep=".")

meltDat <- rbind(depDat,indepDat)

# Trying this...
p <- ggplot(meltDat,aes(x=graphTitle,y=val,col=type))
p + geom_point() +
  geom_errorbar(aes(ymax=val+err,
                    ymin=val-err),width=.25) +
  facet_wrap(~gene,scales='free_y') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title=element_text(size="20"), axis.text=element_text(size='15')) +
  ggtitle("Volume-independent transcription")

# Volume-independent transcription
p <- ggplot(dat,aes(x=type,y=inter.estimate,col=type))
p + geom_point() +
  geom_errorbar(aes(ymax=inter.estimate+inter.stderr,
                    ymin=inter.estimate-inter.stderr),width=.25) +
  facet_wrap(~gene,scales='free_y') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title=element_text(size="20"), axis.text=element_text(size='15')) +
  ggtitle("Volume-independent transcription")

# Volume-dependent transcription
p <- ggplot(dat,aes(x=type,y=slope.estimate*mean,col=type))
p + geom_point() +
  geom_errorbar(aes(ymax=(slope.estimate+slope.stderr)*mean,
                    ymin=(slope.estimate-slope.stderr)*mean),width=.25) +
  facet_wrap(~gene,scales='free_y') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title=element_text(size="20"), axis.text=element_text(size='15')) +
  ggtitle("Volume-dependent transcription")

###

ggplot(crlDat, 
       aes(x=volume,y=totalRNA)) + 
  geom_point() + facet_wrap(~gene, scales = 'free_y') + expand_limits(y=0) + 
  stat_smooth(method="lm", se=FALSE, fullrange = TRUE,col='red') +
  xlab('Volume (femtoliter)') + ylab('mRNA Abundance') + ggtitle("Cycling")

#ggplot(crlDat, 
#       aes(x=log(volume),y=log(totalRNA))) + 
#  geom_point() + facet_wrap(~gene, scales = 'free_y') + expand_limits(y=0) + 
#  stat_smooth(method="lm", se=FALSE, fullrange = TRUE,col='red') +
#  xlab('log Volume (femtoliter)') + ylab('log mRNA Abundance') + ggtitle("Cycling - Log")

ggplot(ssDat, 
       aes(x=volume,y=totalRNA)) + 
  geom_point() + facet_wrap(~gene, scales = 'free_y') + expand_limits(y=0) + 
  stat_smooth(method="lm", se=FALSE, fullrange = TRUE,col='red') +
  xlab('Volume (femtoliter)') + ylab('mRNA Abundance') + ggtitle("Quiescent")

ggplot(senDat, 
       aes(x=volume,y=totalRNA)) + 
  geom_point() + facet_wrap(~gene, scales = 'free_y') + expand_limits(y=0) + 
  stat_smooth(method="lm", se=FALSE, fullrange = TRUE,col='red') +
  xlab('Volume (femtoliter)') + ylab('mRNA Abundance') + ggtitle("Senescent")

ggplot(a549Dat, 
       aes(x=volume,y=totalRNA)) + 
  geom_point() + facet_wrap(~gene, scales = 'free_y') + expand_limits(y=0) + 
  stat_smooth(method="lm", se=FALSE, fullrange = TRUE,col='red') +
  xlab('Volume (femtoliter)') + ylab('mRNA Abundance') + ggtitle("Senescent")

 crlTable <- data.table(crlDat)
crlCorr <- crlTable[ , list(crlDensity = mean(cytoRNA/volume),
                            crlCor = cor(volume,cytoRNA)), by = c('gene')]
crlCorr <- as.data.frame(crlCorr)
