library(ggplot2)
library(reshape)

heights <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/FixationAndVolume/Heights_All.txt',
                    header=T, stringsAsFactors=F)

numericcols <- c('EtOH_Height','Fixed_Height','Live_Height');


for (col in numericcols){
  heights[[col]] <- as.numeric(heights[[col]])
}

heightmelt <- melt(heights,id.var='Classifier')

heightmelt$variable <- factor(heightmelt$variable, levels = c('Live_Height','Fixed_Height','EtOH_Height'))

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/FixationAndVolume/Height.pdf',width=7,height=3)
ggplot(heightmelt, aes(x=Classifier,y=value*.25,fill=variable)) +
  geom_bar(stat='identity',position="dodge") +
  theme_classic() +
  xlab('Cell') + ylab('Cell height (µm)') +
  theme(axis.title=element_text(size="12"), axis.text=element_text(size='12'),
        legend.position='none')
dev.off()
