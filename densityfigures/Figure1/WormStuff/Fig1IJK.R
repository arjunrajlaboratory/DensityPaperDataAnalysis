library("ggplot2")
library("plyr")
library(reshape2)
library(SDMTools)
library(grid)
source('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure1/WormStuff/loadMergeAndMakeDensities.R')
source('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure1/WormStuff/compareN2toCB502_Size.R')

areaRatio <- dataToPlot$mean[dataToPlot$Strain == 'N2']/dataToPlot$mean[dataToPlot$Strain == 'CB502']
volumeRatio <- areaRatio^(1.5) #Approximating the worm as a rectangular prism.
clrs <- c('royalblue','deepskyblue4')
clrs2 <- c('darkmagenta','darkorchid')

## Defining the general functions that will be used to produce the plots.

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

getThisDataReadyForPlotting <- function (data){
  meltedData <- melt(data, 
                     id.vars =c('objArrayNum', 'objNum', 'Strain', 'Volume'), 
                     measure.vars = c("VolDensityAlexa", "VolDensityTmr", 
                                      "PerNucDensityAlexa", "PerNucDensityTmr"))
  dataToPlot <- ddply(meltedData, c("Strain", "variable"), 
                      summarize, WeightedMean=weighted.mean(value, Volume), sd=sd(value), SEM=sd(value)/sqrt(length(value)),
                      weightedSD = wt.sd(value, Volume), weightedSEM = wt.sd(value, Volume)/sqrt(length(value)))
  return(dataToPlot)
}

makeVolComparisonPlot <- function(inputData){
  thisplot <- ggplot(inputData,aes(x=Strain, y=WeightedMean, fill=Strain)) + 
    geom_bar(stat="identity") + 
    geom_errorbar(aes(ymax=WeightedMean+weightedSEM, ymin=WeightedMean-weightedSEM), position=position_dodge(0.9), width=.35) +
    theme_classic() +
    scale_fill_manual(values=c('N2'=clrs[2], 'CB502'=clrs[1])) +
    theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
    ylab(yAxisName) +
    guides(fill=FALSE) 
  return(thisplot)
}

makePerNucComparisonPlot <- function(inputData){
  thisplot <- ggplot(inputData,aes(x=Strain, y=WeightedMean, fill=Strain)) + 
    geom_bar(stat="identity") + 
    geom_errorbar(aes(ymax=WeightedMean+weightedSEM, ymin=WeightedMean-weightedSEM), position=position_dodge(0.9), width=.35) +
    theme_classic() +
    theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
    scale_fill_manual(values=c('N2'=clrs[2], 'CB502'=clrs[1])) +
    ylab(yAxisName) +
    guides(fill=FALSE) 
  return(thisplot)
}

plotStrainComparisonsForAllVariables <- function(dataToPlot, plotFileName){
  pdf(plotFileName)
  d_ply(dataToPlot, 'variable', function(inputData){
    q <- makeStrainComparisonPlot(inputData)
    q <- q + ylab(inputData$variable[1])
    print(q)
  })
  dev.off()
}

## Specifying the input variables for the previously-defined functions.

pathToGraphs <- '~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure1/WormStuff/Gonad'

## Running the functions with the above specified inputs.

thisIsTheDataToBePlotted <- getThisDataReadyForPlotting(
  data =  subset(combinedData, anatomicalStructure %in% 'gonad'))

#To make a pdf with just the alexa volume graph
yAxisName <- "mRNA Volume Density, ama-1 (mRNA/fL)"

p1 <- makeVolComparisonPlot(inputData = subset(thisIsTheDataToBePlotted, variable %in% 'VolDensityAlexa'))

#To make a pdf with just the tmr volume graph
yAxisName <- "mRNA Volume Density, arf-3 (mRNA/fL)"

p2 <- makeVolComparisonPlot(inputData = subset(thisIsTheDataToBePlotted, variable %in% 'VolDensityTmr'))

#To make a pdf with just the alexa per nucleus graph
yAxisName <- "mRNA Molecules per Cell, ama-1"

p3 <- makePerNucComparisonPlot(inputData = subset(thisIsTheDataToBePlotted, variable %in% 'PerNucDensityAlexa'))

#To make a pdf with just the tmr per nucleus graph
yAxisName <- "mRNA Molecules per Cell, arf-3"

p4 <- makePerNucComparisonPlot(inputData = subset(thisIsTheDataToBePlotted, variable %in% 'PerNucDensityTmr'))

#To make all four plots in a 2x2 arrangement.
filename <- paste0('compareN2toCB502', '_', 'allNads', 'WeightedMeans', 'Alexa_vol.pdf')
plotFileName <- file.path(pathToGraphs, filename)
pdf(plotFileName,width=1.3,height=1.7)
p1
dev.off()

filename <- paste0('compareN2toCB502', '_', 'allNads', 'WeightedMeans', 'Tmr_vol.pdf')
plotFileName <- file.path(pathToGraphs, filename)
pdf(plotFileName,width=1.3,height=1.7)
p2
dev.off()

filename <- paste0('compareN2toCB502', '_', 'allNads', 'WeightedMeans', 'Alexa_nuc.pdf')
plotFileName <- file.path(pathToGraphs, filename)
pdf(plotFileName,width=1.3,height=1.7)
p3
dev.off()

filename <- paste0('compareN2toCB502', '_', 'allNads', 'WeightedMeans', 'Tmr_nuc.pdf')
plotFileName <- file.path(pathToGraphs, filename)
pdf(plotFileName,width=1.3,height=1.7)
p4
dev.off()

