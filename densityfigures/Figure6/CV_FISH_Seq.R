# Calculate CV for FISH and seq 

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

bootSeq.cvR2.fpkm <- function(mydata,i){
  d <- mydata[i,]
  cytoRNA <- d$fpkm
  
  out <- (sd(cytoRNA)/mean(cytoRNA))^2
  
  return(out)
}

bootSeq.cvR2.count <- function(mydata,i){
  d <- mydata[i,]
  cytoRNA <- d$actualCount
  
  out <- (sd(cytoRNA)/mean(cytoRNA))^2
  
  return(out)
}

bootsSeqFpkm <- ddply(dataSeq, .(gene_id), function(df) {
  results <- boot(data=df, statistic=bootSeq.cvR2.fpkm, R=5000)
  est <- mean(results$t)
  orderedResults <- results$t[order(results$t)]
  lower <- orderedResults[125]
  upper <- orderedResults[4875]
  out <- c(est,lower,upper)
  names(out) <- c('nm.fpkm','nm.fpkm.lower','nm.fpkm.upper');
  out
})

colnames(bootsSeqFpkm) <- c('gene','CV2_seq_fpkm','lower_seq_fpkm','upper_seq_fpkm')

bootsSeqCount <- ddply(dataSeq, .(gene_id), function(df) {
  results <- boot(data=df, statistic=bootSeq.cvR2.count, R=5000)
  est <- mean(results$t)
  orderedResults <- results$t[order(results$t)]
  lower <- orderedResults[125]
  upper <- orderedResults[4875]
  out <- c(est,lower,upper)
  names(out) <- c('nm.count','nm.count.lower','nm.count.upper');
  out
})

colnames(bootsSeqCount) <- c('gene','CV2_seq_count','lower_seq_count','upper_seq_count')

# Bootstrap for FISH

bootFish.cvR2 <- function(mydata,i){
  d <- mydata[i,]
  cytoRNA <- d$totalRNA
  
  out <- (sd(cytoRNA)/mean(cytoRNA))^2
  
  return(out)
}

bootsFish <- ddply(dataFish, .(gene), function(df) {
  results <- boot(data=df, statistic=bootFish.cvR2, R=5000)
  est <- mean(results$t)
  orderedResults <- results$t[order(results$t)]
  lower <- orderedResults[125]
  upper <- orderedResults[4875]
  out <- c(est,lower,upper)
  names(out) <- c('nm','nm.lower','nm.upper');
  out
})

colnames(bootsFish) <- c('gene','CV2_fish','lower_fish','upper_fish')

#

data <- merge(bootsFish,bootsSeqFpkm)

data <- merge(data,bootsSeqCount)

tmp <- subset(data,!is.na(CV2_seq_fpkm))
tmp <- subset(tmp,!is.na(CV2_seq_count))

ggplot(tmp,aes(x=CV2_fish,y=CV2_seq_fpkm)) + geom_point() + scale_x_log10() + scale_y_log10()

cor(log(tmp$CV2_fish),log(tmp$CV2_seq_fpkm))

write.table(tmp,'~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure6/CV_FISH_Seq_compare.txt',quote=F,sep='\t',row.names=F,col.names=T)