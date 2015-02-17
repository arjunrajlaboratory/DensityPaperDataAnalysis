library(ggplot2)
library(plyr)

data <- read.delim('~/Dropbox/densitypaper/densityfigures/Figure5/Nm_all_genes.txt',header=T,stringsAsFactors=F)

bulkData <- read.delim('~/Dropbox/densitypaper/BulkSeqData/CRL_A549_Bulk_Fpkm.txt',header=T,stringsAsFactors=F)

seqDat <- read.table('~/Dropbox/densitypaper/SingleCellSeqData/seqData_volume_and_counts_allGenes.txt',header=T,stringsAsFactors=F)

bulkData$group <- 'NotConsidered'

bulkData$group[bulkData$fpkmCrl>5*bulkData$fpkmA549 & bulkData$fpkmCrl>5] <- 'crlSpecific'
bulkData$group[bulkData$fpkmA549>5*bulkData$fpkmCrl & bulkData$fpkmA549>5] <- 'a549Specific'
bulkData$group[bulkData$fpkmA549<2*bulkData$fpkmCrl & bulkData$fpkmCrl<2*bulkData$fpkmA549 & bulkData$fpkmCrl>5 & bulkData$fpkmA549>5] <- 'equalExpression'
colnames(bulkData) <- c('gene','fpkmCrl','fpkmA549','group')

data <- merge(data,bulkData[,c(1,4)])

crlSpecGenes <- subset(data,group %in% 'crlSpecific')$gene
nonspecGenes <- subset(data,group %in% 'equalExpression')$gene

seedFromFile <- read.table('~/Dropbox/densitypaper/densityfigures/seed_fig1_2')
.Random.seed <- t(seedFromFile)
tenRandomSpecificGenes <- crlSpecGenes[round(runif(5,min=1,max=length(crlSpecGenes)))]

#randomGeneData <- subset(seqDat,gene_id %in% tenRandomSpecificGenes)
tenRandomSpecificGenes <- c('ICAM1','GAS6','ACTA2')
randomGeneData <- subset(seqDat,gene_id %in% tenRandomSpecificGenes)
randomGeneData <- randomGeneData[,c('sampleID','gene_id','actualCount')]

randomGeneDataWide <- reshape(randomGeneData, idvar = 'sampleID',timevar = 'gene_id', direction = 'wide')

specSeqData <- subset(seqDat,gene_id %in% crlSpecGenes)
specSeqData <- specSeqData[,c('sampleID','gene_id','actualCount')]

nonspecSeqData <- subset(seqDat,gene_id %in% nonspecGenes)
nonspecSeqData <- nonspecSeqData[,c('sampleID','gene_id','actualCount')]

dataToPlayWith <- merge(specSeqData,randomGeneDataWide)

corrTable <- ddply(dataToPlayWith,.(gene_id),summarise,
                   cor1 = cor(actualCount,actualCount.ACTA2),
                   cor2 = cor(actualCount,actualCount.ICAM1),
                   cor3 = cor(actualCount,actualCount.GAS6))

mean(corrTable$cor1)
mean(corrTable$cor2)
mean(corrTable$cor3)
mean(corrTable$cor4)
mean(corrTable$cor5)

dataToPlayNonspec <- merge(nonspecSeqData,randomGeneDataWide)

corrTableNonspec <- ddply(dataToPlayNonspec,.(gene_id),summarise,
                          cor1 = cor(actualCount,parse(text = colnames(corrTableNonspec)[2])),
                          cor2 = cor(actualCount,parse(text = colnames(corrTableNonspec)[3])),
                          cor3 = cor(actualCount,parse(text = colnames(corrTableNonspec)[4])))

mean(corrTableNonspec$cor1)
mean(corrTableNonspec$cor2)
mean(corrTableNonspec$cor3)
mean(corrTableNonspec$cor4)
mean(corrTableNonspec$cor5)

twoTables <- rbind(corrTable,corrTableNonspec)

matrixThing <- as.matrix(corrTable[,c(2:4)])
heatmap(matrixThing)



eval(parse(text = "x"))