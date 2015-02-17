# This code was modified from mu_analysis_2.R (OPM, 140330) #

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
  dat[[col]] <- as.numeric(dat[[col]])
}

dat$volume <- dat$volume/1000

lmDat <- ddply(dat,c('gene','type'),function(x) {
  tmp <- summary(lm(cytoRNA~volume,x))$coefficients; 
  out <- c(tmp[1,], tmp[2,]);
  names(out) <- paste(rep(c("inter", "slope"), each=ncol(tmp)), 
                      c('estimate','stderr','tval','pr'), sep=".");
  out
})

meanDat <- ddply(dat,~type,summarise,mean=mean(volume))

dat <- merge(lmDat,meanDat)

dat$volindep.val <- dat$inter.estimate
dat$volindep.err <- dat$inter.stderr
dat$voldep.val <- dat$slope.estimate*dat$mean
dat$voldep.err <- dat$slope.stderr*dat$mean

crlDat <- subset(dat,type %in% 'cycling')
ssDat <- subset(dat,type %in% 'quiescent')
senDat <- subset(dat,type %in% 'senescent')
a549Dat <- subset(dat,type %in% 'a549')

ggplot(crlDat,aes(x=voldep.val,y=volindep.val,label=gene)) +
  geom_point(size=1.5) +
  scale_x_log10() + scale_y_log10() +
  expand_limits(x=c(5,2000),y=c(5,2000)) +
  geom_abline(intercept=0,slope=1)

ggplot(ssDat,aes(x=voldep.val,y=volindep.val,label=gene)) +
  geom_point(size=1.5) +
  scale_x_log10() + scale_y_log10() +
  geom_abline(intercept=0,slope=1)

ggplot(senDat,aes(x=voldep.val,y=volindep.val,label=gene)) +
  geom_point() +
  geom_abline(intercept=0,slope=1)

ggplot(a549Dat,aes(x=voldep.val,y=volindep.val,label=gene)) +
  geom_point() +
  geom_abline(intercept=0,slope=1)

#indepDat <- dat[,c('type','gene','volindep.val','volindep.err')]
#indepDat <- unique(indepDat)
#colnames(indepDat) <- c('type','gene','val','err')
#indepDat$voldep <- 'indep'
#indepDat$graphTitle <- paste(indepDat$type,indepDat$voldep, sep=".")

#depDat <- dat[,c('type','gene','voldep.val','voldep.err')]
#depDat <- unique(depDat)
#colnames(depDat) <- c('type','gene','val','err')
#depDat$voldep <- 'dep'
#depDat$graphTitle <- paste(depDat$type,depDat$voldep, sep=".")

#meltDat <- rbind(depDat,indepDat)

# Trying this...
#p <- ggplot(meltDat,aes(x=graphTitle,y=val,col=type))
#p + geom_point() +
#  geom_errorbar(aes(ymax=val+err,
#                    ymin=val-err),width=.25) +
#  facet_wrap(~gene,scales='free_y') +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#        axis.title=element_text(size="20"), axis.text=element_text(size='15')) +
#  ggtitle("Volume-independent transcription")

# Volume-independent transcription
#p <- ggplot(dat,aes(x=type,y=inter.estimate,col=type))
#p + geom_point() +
#  geom_errorbar(aes(ymax=inter.estimate+inter.stderr,
#                    ymin=inter.estimate-inter.stderr),width=.25) +
#  facet_wrap(~gene,scales='free_y') +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#        axis.title=element_text(size="20"), axis.text=element_text(size='15')) +
#  ggtitle("Volume-independent transcription")

# Volume-dependent transcription
#p <- ggplot(dat,aes(x=type,y=slope.estimate*mean,col=type))
#p + geom_point() +
#  geom_errorbar(aes(ymax=(slope.estimate+slope.stderr)*mean,
#                    ymin=(slope.estimate-slope.stderr)*mean),width=.25) +
#  facet_wrap(~gene,scales='free_y') +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#        axis.title=element_text(size="20"), axis.text=element_text(size='15')) +
#  ggtitle("Volume-dependent transcription")