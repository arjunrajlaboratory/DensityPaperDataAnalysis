library(reshape2)
library(ggplot2)

data <- read.delim('~/Dropbox/densitypaper/ExtractedData/ActD.txt',header=T,stringsAsFactors=F)

numericcols <- c('timePt','volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');


for (col in numericcols){
  data[[col]] <- as.numeric(data[[col]])
}

data$volume <- data$volume/1000

# let's make a column called replicate number

data$repNum <- NA

for(i in 1:length(unique(data$gene))) {
  geneiter <- as.character(unique(data$gene))[i]
  tmp <- subset(data,gene %in% geneiter)
  for(j in 1:length(unique(tmp$date))){
    dateiter <- as.character(unique(tmp$date))[j]
    data$repNum[(data$gene==geneiter & data$date ==dateiter)] <- j
  }
}

###
# gene = 'UBC', rep = 3
# gene = 'IER2', rep = 1
selectGene <- 'IER2'
selectRep <- 1
selectData <- subset(data,gene %in% selectGene & repNum %in% selectRep)
timePts <- unique(selectData$timePt)

# Fit untreated line through zero
fit0h <- lm(cytoRNA ~ 0 + volume, data=subset(selectData, timePt %in% timePts[1]))
# Find slope and intercept
coeff0h <- coefficients(fit0h)

#Interpolated model
interpModel6h <- nls(cytoRNA ~ volume*coeff0h[1]*exp(-B*(1/volume)^alpha),
                     data=subset(selectData, timePt %in% timePts[2]),
                     start=list(B=1,alpha=0))
interpCoeff <- coefficients(interpModel6h)

#Linear model
linModel6h <- nls(cytoRNA ~ volume*coeff0h[1]*exp(-B),
                  data=subset(selectData, timePt %in% timePts[2]),
                  start=list(B=1))
linCoeff <- coefficients(linModel6h)

#1/V model
vModel6h <- nls(cytoRNA ~ volume*coeff0h[1]*exp(-B/volume),
                data=subset(selectData, timePt %in% timePts[2]),
                start=list(B=1))
vCoeff <- coefficients(vModel6h)

pdf('~/Dropbox/densitypaper/SupplementaryFigures/ActD/rawdata_IER2.pdf',width=2,height=1.8)
ggplot(selectData, aes(x=volume,y=cytoRNA,color=factor(timePt))) + 
  xlab('Volume (picoliter)') + ylab('Total RNA') +
  #stat_function(fun=function(x) x*coeff0h[1]*exp(-interpCoeff[1]*(1/x)^interpCoeff[2])) +
  stat_function(fun=function(x) x*coeff0h[1]*exp(-linCoeff[1]),col='deeppink4') +
  stat_function(fun=function(x) x*coeff0h[1]*exp(-vCoeff[1]/x),col='dodgerblue4') +
  geom_abline(intercept = 0, slope = coeff0h) +
  geom_point(size=.8) + 
  theme_classic() +
  theme(legend.position='none') +
  expand_limits(y=0) + expand_limits(x=0) + 
  #theme(axis.title=element_text(size="4"), axis.text=element_text(size='4')) +
  scale_colour_manual(values=c('0' = 'gray40',
                               '4' = 'gray40'))  
dev.off()

coeffTable <- summary(interpModel6h)$coefficients
minAlpha <- coeffTable['alpha','Estimate'] - 1.96*coeffTable['alpha','Std. Error']
maxAlpha <- coeffTable['alpha','Estimate'] + 1.96*coeffTable['alpha','Std. Error']
coeffTable['alpha','Estimate']
minAlpha
maxAlpha