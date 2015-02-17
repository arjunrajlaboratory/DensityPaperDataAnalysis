source('~/Dropbox/densitypaper/densitypaperdataanalysis/multiplot.R')
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

### Read in click data

click <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/Click_SS.txt',header=T,stringsAsFactors=F)

numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'dapiIntensityTotal','dapiIntensityAvg','clickIntensityAvg',
                 'clickIntensityTotal','numCyclin');


for (col in numericcols){
  click[[col]] <- as.numeric(click[[col]])
}

click$volume <- click$volume/1000
click$clickIntensityTotal <- click$clickIntensityTotal/10^7

#dates = 130118, 130219, 130418
#dates (SS) = 121113, 130116, 130405

tmpclick <- subset(click,date %in% '121113')


###
# gene = 'UBC', rep = 3
# gene = 'IER2', rep = 1
selectGene <- 'UBC'
selectRep <- 3
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

deg <- ggplot(tmp,aes(x=volume,y=1/tau)) +
  geom_abline(slope=0,intercept=linCoeff,col=cols[1]) +
  stat_function(fun=test,col=cols[2]) +
  xlab('Volume (picoliter)') + ylab('Decay constant') +
  theme_bw() + theme_classic() + theme(legend.position = "none") + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  geom_point(size=1.3)

# now transcription

fittxn <- (lm(clickIntensityTotal ~ volume, data=tmpclick))
inttxn <- fittxn$coefficients[[1]]
slopetxn <- fittxn$coefficients[[2]]

txn <- ggplot(data=tmpclick,aes(x=volume,y=clickIntensityTotal)) +
  expand_limits(x=c(0,3.2),y=c(0,8.5)) +
  geom_abline(slope=slopetxn,intercept=inttxn,col=cols[1]) +
  geom_abline(slope=0,intercept=mean(tmpclick$clickIntensityTotal),col=cols[2]) +
  xlab('Volume (picoliter)') + ylab('Global transcription rate (arb. units)') +
  theme_bw() + theme_classic() + theme(legend.position = "none") + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  geom_point(size=1.3) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) 

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure2/TxnDeg/Txn_Deg.pdf',height=1.8,width=4)
multiplot(deg,txn,cols=2)
dev.off()

