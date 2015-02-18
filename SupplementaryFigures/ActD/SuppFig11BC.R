library(reshape2)
library(ggplot2)

data <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/ActD.txt',header=T,stringsAsFactors=F)

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
linModel6h <- nls(cytoRNA ~ volume*coeff0h[1]*exp(-B*timePts[2]),
                  data=subset(selectData, timePt %in% timePts[2]),
                  start=list(B=1))
linCoeff <- coefficients(linModel6h)

#1/V model
vModel6h <- nls(cytoRNA ~ volume*coeff0h[1]*exp(-B*timePts[2]/volume),
                data=subset(selectData, timePt %in% timePts[2]),
                start=list(B=1))
vCoeff <- coefficients(vModel6h)



selectData$tau <- timePts[2]/log(selectData$volume*coeff0h/selectData$cytoRNA)

tmp <- subset(selectData,timePt>0)

cols <- c('deeppink4','dodgerblue4')

test <- function(x) {vCoeff/x}

pdf(file='~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/ActD/Degradation_IER2.pdf',height=4,width=4.4)
ggplot(tmp,aes(x=volume,y=1/tau)) +
  geom_abline(slope=0,intercept=linCoeff,col=cols[1]) +
  stat_function(fun=test,col=cols[2]) +
  xlab('Volume (picoliter)') + ylab('Decay constant') +
  theme_bw() + theme_classic() + theme(legend.position = "none") + 
  theme(axis.title=element_text(size="12"), axis.text=element_text(size='12')) +
  geom_point(size=1.3)
dev.off()

pdf(file='~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/ActD/VolumeBeforeAfter.pdf',height=3,width=3.2)
ggplot(subset(data,timePt %in% c(0,4)),aes(x=volume,fill=factor(timePt))) +
  geom_density(alpha=.25) +
  xlab('Volume (picoliter)') + ylab('Density') +
  theme_bw() + theme_classic() + 
  theme(axis.title=element_text(size="12"), axis.text=element_text(size='12'),legend.position='none')
dev.off()