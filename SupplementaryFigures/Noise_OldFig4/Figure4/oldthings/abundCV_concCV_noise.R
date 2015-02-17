library(ggplot2)
library(plyr)
library(reshape)

seedFromFile <- read.table('~/Dropbox/densitypaper/densityfigures/seed_marginalhist')
.Random.seed <- t(seedFromFile)

data <- read.delim('~/Dropbox/ExtractedData_131216/alldata_CRL.txt',header=T,stringsAsFactors=F)
data <- subset(data,gene %in% c('GAPDH','EEF2','MYC','LMNA','USF2','UBC'))

numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');


for (col in numericcols){
  data[[col]] <- as.numeric(data[[col]])
}

data$volume <- data$volume/1000

#noise <- ddply(data,.(gene),function(df) {
#  cytoRNA <- df$cytoRNA
#  volume <- df$volume
#  fit <- lm(cytoRNA ~ volume, data = df)
#  slope <- fit$coeff[[2]]
#  int <- fit$coeff[[1]]
#  
#  out <- (sd(cytoRNA)/mean(cytoRNA))^2 - 
#    slope*mean(volume)/(int+slope*mean(volume)) * cov(cytoRNA,volume)/(mean(cytoRNA)*mean(volume))
#  
#  return(out)
#})

#colnames(noise) <- c('gene','nm')

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

boots <- ddply(data, .(gene), function(df) {
  results <- boot(data=df, statistic=boot.nm, R=5000)
  est <- mean(results$t)
  orderedResults <- results$t[order(results$t)]
  lower <- orderedResults[125]
  upper <- orderedResults[4875]
  out <- c(est,lower,upper)
  names(out) <- c('nm','lower','upper');
  out
})

cvTable <- ddply(data,.(gene),summarize,
                 cvR2 = (sd(cytoRNA)/mean(cytoRNA))^2,
                 cvC2 = (sd(cytoRNA/volume)/mean(cytoRNA/volume))^2)

outTable <- merge(cvTable,boots)
outTable <- outTable[,c('gene','cvR2','cvC2','nm')]

outMelt <- melt(outTable,idvar='gene')

cols = c('slateblue4','turquoise4','seagreen3')

pdf('~/Dropbox/densitypaper/densityfigures/Figure4/gene_noise_measures.pdf',width=2.5,height=2)
ggplot(outMelt,aes(x=gene,y=value,fill=variable)) +
  geom_bar(position='dodge',stat='identity',col='black') +
  theme_classic() +
  scale_fill_manual(values=cols) +
  theme(legend.position = "none") + 
  xlab('Gene') + ylab('Noise Measure') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()

