library(ggplot2)

data <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/histonecyclin.txt',
                   header=T,stringsAsFactors=F)

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/CyclinHistone/CyclinHistone.pdf',width=5,height=3)
ggplot(data,aes(x=numCyclin,y=numHistone)) +
  geom_point(size=1.5) +
  theme_classic() +
  geom_vline(xintercept=20,col='gray60',linetype='dashed') +
  geom_vline(xintercept=230,col='gray60',linetype='dashed') +
  theme(axis.title=element_text(size="12"), axis.text=element_text(size='12')) +
  xlab('Cyclin mRNA') + ylab('Histone mRNA')
dev.off()
