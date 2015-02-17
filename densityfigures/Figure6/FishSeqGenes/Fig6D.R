library(ggplot2)

fishData <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/CRL.txt',header=T,stringsAsFactors=F)
fishData <- fishData[,c('gene','volume','totalRNA')]
colnames(fishData) <- c('gene','volume','count')
fishData$type <- 'fish'
seqData <- read.table('~/Dropbox/densitypaper/densitypaperdataanalysis/SingleCellSeqData/seqData_volume_and_counts_fishGenes.txt',header=T,stringsAsFactors=F)
seqData <- seqData[,c('gene_id','volume','actualCount')]
colnames(seqData) <- c('gene','volume','count')
seqData$type <- 'seq'

gene1 <- 'GAPDH'
gene2 <- 'MYC'

fishSubset <- subset(fishData,gene %in% c(gene1,gene2))
seqSubset <- subset(seqData,gene %in% c(gene1,gene2))

data <- rbind(fishSubset,seqSubset)

p1 <- ggplot(subset(data,gene %in% gene1 & type %in% 'fish'),aes(x=volume/1000,y=count)) +
  geom_point(size=1.1,col='mediumpurple3') + xlim(c(0,6)) + #ylim(c(0,250)) +
  theme_classic() + theme(legend.position='none') +
  xlab('Volume (picoliter)') + ylab('FISH count') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

p2 <- ggplot(subset(data,gene %in% gene1 & type %in% 'seq'),aes(x=volume/1000,y=count)) + 
  geom_point(size=1.1,col='mediumpurple3') + xlim(c(0,6)) + #ylim(c(0,250)) +
  theme_classic() + theme(legend.position='none') +
  xlab('Volume (picoliter)') + ylab('Estimated count (single-cell seq)') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

p3 <- ggplot(subset(data,gene %in% gene2 & type %in% 'fish'),aes(x=volume/1000,y=count)) + 
  geom_point(size=1.1,col='lightseagreen') + xlim(c(0,6)) + #ylim(c(0,320)) +
  theme_classic() + theme(legend.position='none') +
  xlab('Volume (picoliter)') + ylab('FISH count') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

p4 <- ggplot(subset(data,gene %in% gene2 & type %in% 'seq'),aes(x=volume/1000,y=count)) + 
  geom_point(size=1.1,col='lightseagreen') + xlim(c(0,6)) + #ylim(c(0,320)) +
  theme_classic() + theme(legend.position='none') +
  xlab('Volume (picoliter)') + ylab('Estimated count (single-cell seq)') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

pdf(file='~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure6/FishSeqGenes/GAPDH_FISH.pdf',height=1.5,width=1.7)
p1
dev.off()

pdf(file='~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure6/FishSeqGenes/GAPDH_Seq.pdf',height=1.5,width=1.7)
p2
dev.off()

pdf(file='~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure6/FishSeqGenes/MYC_FISH.pdf',height=1.5,width=1.7)
p3
dev.off()

pdf(file='~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure6/FishSeqGenes/MYC_Seq.pdf',height=1.5,width=1.7)
p4
dev.off()