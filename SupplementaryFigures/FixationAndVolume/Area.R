library(ggplot2)
library(reshape)

areas <- read.delim('~/Dropbox/densitypaper/SupplementaryFigures/FixationAndVolume/Areas_All.txt',
                    header=T, stringsAsFactors=F)

numericcols <- c('Ethanol','Fixed','Live');


for (col in numericcols){
  areas[[col]] <- as.numeric(areas[[col]])
}

areamelt <- melt(areas,id.var='Classifier')

areamelt$variable <- factor(areamelt$variable, levels = c('Live','Fixed','Ethanol'))

pdf('~/Dropbox/densitypaper/SupplementaryFigures/FixationAndVolume/Area.pdf',width=7,height=3)
ggplot(areamelt, aes(x=Classifier,y=value*.125*.125,fill=variable)) +
  geom_bar(stat='identity',position="dodge") +
  theme_classic() +
  xlab('Cell') + ylab('Cell area (µßm^2)') +
  theme(axis.title=element_text(size="12"), axis.text=element_text(size='12'),
        legend.position='none')
dev.off()
