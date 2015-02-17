library(ggplot2)
library(plyr)
library(reshape)
library(boot)

seedFromFile <- read.table('~/Dropbox/densitypaper/densityfigures/seed_marginalhist')
.Random.seed <- t(seedFromFile)

data <- read.delim('~/Dropbox/ExtractedData_131216/alldata_Senescent.txt',header=T,stringsAsFactors=F)

numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');


for (col in numericcols){
  data[[col]] <- as.numeric(data[[col]])
}

data$volume <- data$volume/1000

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

boot.cvR2 <- function(mydata,i){
  d <- mydata[i,]
  cytoRNA <- d$cytoRNA
  
  out <- (sd(cytoRNA)/mean(cytoRNA))^2
  
  return(out)
}

boot.cvC2 <- function(mydata,i){
  d <- mydata[i,]
  cytoRNA <- d$cytoRNA
  volume <- d$volume
  
  out <- (sd(cytoRNA/volume)/mean(cytoRNA/volume))^2
  
  return(out)
}

boots <- ddply(data, .(gene), function(df) {
  results <- boot(data=df, statistic=boot.nm, R=5000)
  est <- mean(results$t)
  orderedResults <- results$t[order(results$t)]
  lower <- orderedResults[125]
  upper <- orderedResults[4875]
  out <- c(est,lower,upper)
  names(out) <- c('nm','nm.lower','nm.upper');
  out
})

boots$measure <- 'Nm'
colnames(boots) <- c('gene','value','lower','upper','measure')

cvR2Err <- ddply(data, .(gene), function(df) {
  results <- boot(data=df, statistic=boot.cvR2, R=5000)
  orderedResults <- results$t[order(results$t)]
  meas <- results$t0
  lower <- orderedResults[125]
  upper <- orderedResults[4875]
  
  out <- c(meas,lower, upper)
  names(out) <- c('value','lower','upper');
  out
})

cvR2Err$measure <- 'cvR2'

cvC2Err <- ddply(data, .(gene), function(df) {
  
  results <- boot(data=df, statistic=boot.cvC2, R=5000)
  orderedResults <- results$t[order(results$t)]
  meas <- results$t0
  lower <- orderedResults[125]
  upper <- orderedResults[4875]
  
  out <- c(meas,lower,upper)
  names(out) <- c('value','lower','upper');
  out
})

cvC2Err$measure <- 'cvC2'

#outTable <- merge(cvC2Err,boots)
#outTable <- merge(outTable,cvR2Err)

outMelt <- rbind(cvR2Err,cvC2Err,boots)

#outMelt <- melt(outTable,idvar='gene')

outMelt$measure <- factor(outMelt$measure, levels = c('cvR2','cvC2','Nm'))

cols = c('slateblue4','turquoise4','seagreen3')
dodge <- position_dodge(.9)

pdf('~/Dropbox/densitypaper/SupplementaryFigures/Noise/Senescent.pdf',width=2.5,height=3)
ggplot(outMelt,aes(x=gene,y=value,fill=measure)) +
  geom_bar(position='dodge',stat='identity',col='black',width=.9) +
  geom_errorbar(aes(ymin=lower,ymax=upper),position=dodge,width=.3) +
  theme_classic() +
  scale_fill_manual(values=cols) +
  theme(legend.position = "none") + 
  xlab('Gene') + ylab('Noise Measure') +
  theme(axis.title=element_text(size="10"), axis.text=element_text(size='10'),
        axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

