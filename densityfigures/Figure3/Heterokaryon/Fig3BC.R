library(ggplot2)
library(plyr)

data <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/heterokaryon.txt',header=T,stringsAsFactors=F)

numericcols <- c('repNum','dataNum','objNum','volume','cytoGAPDH','cytoGFP','cytoGAS6');

for (col in numericcols){
  data[[col]] <- as.numeric(data[[col]])
}

data$volume <- data$volume/1000

# this gets rid of CRL controls (only keep WM controls)
data <- subset(data,!(class %in% 'ctrl' & cytoGFP < 100))

# this gets rid of WM-WM homokaryons
data <- subset(data,!(class %in% '12h' & cytoGAS6<50))

fm <- lm(cytoGFP ~ volume + 0, data=subset(data,class %in% 'ctrl' & volume < 5))
slope <- fm$coefficients[[1]]

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure3/Heterokaryon/heterokaryon.pdf',width=2.9,height=2.7)
ggplot(data, aes(x=volume,y=cytoGFP,color=factor(class))) + geom_point(size = 1.3) +
  expand_limits(x=c(0,10)) + expand_limits(y=c(0,4500)) + geom_abline(slope = slope, intercept = 0, linetype = 'dashed') +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  geom_abline(slope = slope/2, intercept = 0, linetype = 'dashed') +
  xlab('Volume') + ylab('GFP mRNA') +
  theme_bw() + theme_classic() + theme(legend.position = "none") + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  scale_colour_manual(values=c('darkcyan','darkslateblue'))
dev.off()

data$class <- factor(data$class, levels = c('ctrl','12h'))

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure3/Heterokaryon/rnaQuant.pdf',width=1,height=2.7)
ggplot(data, aes(x=class,y=cytoGFP,fill=class)) +
  geom_boxplot() +
  scale_fill_manual(values=c('darkcyan','darkslateblue')) +
  xlab('') + ylab('GFP mRNA') +
  theme_bw() + theme_classic() + theme(legend.position = "none") + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  expand_limits(y=c(0,4500)) +
  scale_y_continuous(expand = c(0, 0))
dev.off()

