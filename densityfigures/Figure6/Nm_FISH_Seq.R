# Calculate Nm for FISH and seq 

library(ggplot2)
library(plyr)
library(reshape)
library(boot)

seedFromFile <- read.table('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/seed_marginalhist')
.Random.seed <- t(seedFromFile)

dataFish <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/CRL.txt',header=T,stringsAsFactors=F)
dataSeq <- read.table('~/Dropbox/densitypaper/densitypaperdataanalysis/SingleCellSeqData/seqData_volume_and_counts_fishGenes.txt',header=T,stringsAsFactors=F)

numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');

for (col in numericcols){
  dataFish[[col]] <- as.numeric(dataFish[[col]])
}

# Bootstrap for Seq

bootSeq.nm <- function(mydata,i){
  d <- mydata[i,]
  cytoRNA <- d$actualCount
  volume <- d$volume
  fit <- lm(cytoRNA ~ volume, data = d)
  slope <- fit$coeff[[2]]
  int <- fit$coeff[[1]]
  
  out <- (sd(cytoRNA)/mean(cytoRNA))^2 - 
    slope*mean(volume)/(int+slope*mean(volume)) * cov(cytoRNA,volume)/(mean(cytoRNA)*mean(volume))
  
  return(out)
}

bootsSeq <- ddply(dataSeq, .(gene_id), function(df) {
  results <- boot(data=df, statistic=bootSeq.nm, R=5000)
  est <- mean(results$t)
  orderedResults <- results$t[order(results$t)]
  lower <- orderedResults[125]
  upper <- orderedResults[4875]
  out <- c(est,lower,upper)
  names(out) <- c('nm','nm.lower','nm.upper');
  out
})

colnames(bootsSeq) <- c('gene','Nm_seq','lower_seq','upper_seq')

# Bootstrap for FISH

bootFish.nm <- function(mydata,i){
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

bootsFish <- ddply(dataFish, .(gene), function(df) {
  results <- boot(data=df, statistic=bootFish.nm, R=5000)
  est <- mean(results$t)
  orderedResults <- results$t[order(results$t)]
  lower <- orderedResults[125]
  upper <- orderedResults[4875]
  out <- c(est,lower,upper)
  names(out) <- c('nm','nm.lower','nm.upper');
  out
})

colnames(bootsFish) <- c('gene','Nm_fish','lower_fish','upper_fish')

#

data <- merge(bootsFish,bootsSeq)

tmp <- subset(data,!is.na(Nm_seq))

ggplot(tmp,aes(x=Nm_fish,y=Nm_seq)) + geom_point() + scale_x_log10() + scale_y_log10()

cor(log(tmp$Nm_fish),log(tmp$Nm_seq))

write.table(tmp,'~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure6/Nm_FISH_Seq_compare.txt',quote=F,sep='\t',row.names=F,col.names=T)


### Now for all the genes. (To do this, I changed R from 5000 to 1000)

dataSeqAll <- read.table('~/Dropbox/densitypaper/densitypaperdataanalysis/SingleCellSeqData/seqData_volume_and_counts_allGenes.txt',header=T,stringsAsFactors=F)

bootsSeqAll <- ddply(dataSeqAll, .(gene_id), function(df) {
  results <- boot(data=df, statistic=bootSeq.nm, R=1000)
  est <- mean(results$t)
  orderedResults <- results$t[order(results$t)]
  lower <- orderedResults[25]
  upper <- orderedResults[975]
  out <- c(est,lower,upper)
  names(out) <- c('nm','nm.lower','nm.upper');
  out
})

colnames(bootsSeqAll) <- c('gene','Nm_seq','lower_seq','upper_seq')

meanCountSeq <- ddply(dataSeqAll,.(gene_id),summarize,
                      meanCount = mean(actualCount),
                      seCount = sd(actualCount)/sqrt(length(actualCount)))
colnames(meanCountSeq) <- c('gene','meanCount','seCount')

NmAndCount <- merge(bootsSeqAll,meanCountSeq)
NmAndCount <- subset(NmAndCount, !(is.na(Nm_seq)))

ggplot(NmAndCount,aes(x=meanCount,y=Nm_seq)) + geom_point()

write.table(NmAndCount,'~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure6/Nm_all_genes.txt',quote=F,sep='\t',row.names=F,col.names=T)
