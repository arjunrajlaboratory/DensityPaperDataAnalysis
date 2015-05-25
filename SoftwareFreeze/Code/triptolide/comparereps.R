library(ggplot2)
library(plyr)

ctrl1 <- read.delim('~/Documents/Data/141113_CRL_MeOH_Triptolide/141113_Triptolide_Ctrl/141113_Triptolide_Ctrl.txt',
                   header=T,stringsAsFactors=F)
ctrl1$repNum <- 1
drug1 <- read.delim('~/Documents/Data/141113_CRL_MeOH_Triptolide/141113_Triptolide_100nM/141113_Triptolide_100nM.txt',
                   header=T,stringsAsFactors=F)
drug1$repNum <- 1

ctrl2 <- read.delim('~/Documents/Data/141114_CRL_MeOH_Triptolide/141114_Triptolide_Ctrl/141114_Triptolide_Ctrl.txt',
                   header=T,stringsAsFactors=F)
ctrl2$repNum <- 2
drug2 <- read.delim('~/Documents/Data/141114_CRL_MeOH_Triptolide/141114_Triptolide_100nM/141114_Triptolide_100nM.txt',
                   header=T,stringsAsFactors=F)
drug2$repNum <- 2

dat <- rbind(ctrl1,ctrl2,drug1,drug2)
dat$gene <- 'MYC'

ggplot(dat,aes(x=intensity,fill=drugConc)) + geom_density(alpha=.25)
ggplot(subset(dat,repNum==1),aes(x=intensity,fill=drugConc)) + geom_density(alpha=.25)
ggplot(subset(dat,repNum==2),aes(x=intensity,fill=drugConc)) + geom_density(alpha=.25)

t.test(subset(dat,drugConc %in% 'Ctrl')$intensity,subset(dat,drugConc %in% '100nM')$intensity)
t.test(subset(dat,drugConc %in% 'Ctrl' & repNum==1)$intensity,subset(dat,drugConc %in% '100nM' & repNum==1)$intensity)
t.test(subset(dat,drugConc %in% 'Ctrl' & repNum==2)$intensity,subset(dat,drugConc %in% '100nM' & repNum==2)$intensity)

ggplot(dat,aes(x=drugConc,y=intensity)) + geom_boxplot()

intSummary <- ddply(dat,.(drugConc,dataNum,objNum),summarise,maxInt=max(intensity))
ggplot(intSummary,aes(x=maxInt,fill=drugConc)) + geom_density(alpha=.25)



uniqueNumTxnSites <- unique(dat[,c('drugConc','dataNum','objNum','numTxnSites')])
tmp <- ddply(uniqueNumTxnSites,.(drugConc),summarise,freq=sum(numTxnSites)/length(numTxnSites))
