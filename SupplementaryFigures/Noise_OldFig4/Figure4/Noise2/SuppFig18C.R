library(ggplot2)
library(plyr)
library(reshape)
library(boot)

seedFromFile <- read.table('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/seed_marginalhist')
.Random.seed <- t(seedFromFile)

data <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/CRL.txt',header=T,stringsAsFactors=F)

numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');


for (col in numericcols){
  data[[col]] <- as.numeric(data[[col]])
}

data$volume <- data$volume/1000

halflife <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/halflife.txt',header=T,stringsAsFactors=F)
hl <- halflife[,c('gene_name','halflife')]
colnames(hl) <- c('gene','halflife')
hl <- subset(hl,halflife!='N.D.')
hl$halflife[hl$halflife=='>24'] = 24
hl$halflife <- as.numeric(hl$halflife)

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

outDat <- merge(boots, hl)

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/Noise_OldFig4/Figure4/Noise2/halflife_noise.pdf',width=3.3,height=3)
ggplot(outDat,aes(x=halflife,y=value)) +
  geom_point(size=1.3) +
  theme_classic() +
  geom_errorbar(aes(ymin=lower,ymax=upper),width=0.5) +
  ylab('Volume-corrected noise measure') + xlab('mRNA Halflife') +
  theme(axis.title=element_text(size='6'), axis.text=element_text(size='6'))
dev.off()
