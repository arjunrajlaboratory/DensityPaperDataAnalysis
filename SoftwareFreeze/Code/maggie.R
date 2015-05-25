data <- read.delim('~/Documents/Data/131220_CRL_MeOH_UBC/TxnSites_IntensityNumCyclin.txt',header=F,stringsAsFactors=F)
colnames(data) <- c('datNum','objNum','intensity','numInt','cyclin')
data <- data[c('datNum','objNum','numInt')]

data <- unique(data)