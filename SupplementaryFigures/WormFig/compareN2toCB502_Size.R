##First, load the data and melt it into a useful form.
pathToData <- '~/Dropbox/densitypaper/WormData/140710_wormdensity_analysis_sizeComp'
relPathToN2 <- 'N2/totalextraction.csv'
relPathToCB502 <- 'CB502/totalextraction.csv'

big_size <- read.csv(file.path(pathToData, relPathToN2), header=T,stringsAsFactors=F)
small_size <- read.csv(file.path(pathToData, relPathToCB502),header=T,stringsAsFactors=F)

big_size <- cbind(big_size, Strain='N2')
small_size <- cbind(small_size, Strain='CB502')

sizeData <- rbind(big_size, small_size)

relPaths <- data.frame(Strain=c('N2', 
                                'CB502'), 
                       relPathTo=c('N2/totalextraction.csv', 
                                   'CB502/totalextraction.csv') )

sizeData <- ddply(relPaths, 'Strain', function(x){
  read.csv(file.path(pathToData, x$relPathTo), header=T,stringsAsFactors=F)
})

sizeData$area_sqMicrons <- sizeData$area*.125*.125

meltedSizeData <- melt(sizeData, 
                       id.vars =c('objArrayNum', 'Strain'), 
                       measure.vars = c("area_sqMicrons"))
dataToPlot <- ddply(meltedSizeData, c("Strain", "variable"), 
                    summarize, mean=mean(value), sd=sd(value), SEM=sd(value)/sqrt(length(value)))

## This is a function that is called inside of the next function, which generates and saves the graph in PDF form.
makeSizeComparisonPlot <- function(inputData){
  q <- ggplot(inputData,aes(x=Strain, y=mean, fill=Strain)) + 
    geom_bar(stat="identity") + 
    geom_errorbar(aes(ymax=mean+SEM, ymin=mean-SEM), position=position_dodge(0.9), width=.35)
return(q)
}

## Modified form of Gautham's function to save graphs as PDFs.
plotSizeComparison <- function(dataToPlot, plotFileName){
  pdf(plotFileName,height=2.2,width=1.5)
  d_ply(dataToPlot, 'variable', function(inputData){
    q <- makeSizeComparisonPlot(inputData)
    q <- q + ylab(inputData$variable[1]) + theme_classic() + theme(legend.position='none') +
      theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
    print(q)
  })
  dev.off()
}
 ## Define where the graph will go and what it will be called.
pathToGraphs <- '~/Dropbox/densitypaper/SupplementaryFigures/WormFig/'
fileName <- file.path(pathToGraphs, 'compareN2toCB502_size.pdf')

## Run the script.
plotSizeComparison(
  dataToPlot = dataToPlot, 
  plotFileName = fileName)