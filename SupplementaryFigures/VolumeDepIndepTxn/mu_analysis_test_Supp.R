library(reshape2)
library(ggplot2)

data <- read.delim('~/Dropbox/ExtractedData_131216/alldata_CRL.txt',header=T,stringsAsFactors=F)
data$volume <- data$volume/1000

numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');


for (col in numericcols){
  data[[col]] <- as.numeric(data[[col]])
}

data <- data[,c('cellType','gene','date','dataNum','objNum','volume','cytoRNA')]

geneList <- c()
mu0List <- c()
muVList <- c()
se0List <- c()
seVList <- c()

for(i in 1:length(unique(data$gene))) {
  geneiter <- as.character(unique(data$gene))[i]
  tmp <- subset(data,gene %in% geneiter)
  
  fit <- lm(cytoRNA ~ volume, data = tmp)
  s <- summary(fit)
  mu0 <- coef(s)[1]
  muV <- coef(s)[2]
  se0 <- coef(s)[3]
  seV <- coef(s)[4]
  muVnorm <- muV*mean(tmp$volume)
  normFactor <- mu0+muVnorm
  
  geneList <- c(geneList,geneiter)
  mu0List <- c(mu0List,mu0/normFactor)
  muVList <- c(muVList,muVnorm/normFactor)
  se0List <- c(se0List,se0/normFactor)
  seVList <- c(seVList,seV*mean(tmp$volume)/normFactor)
  
}

dfmu <- data.frame(
  gene = geneList,
  mu0 = mu0List,
  muV = muVList
)

df <- data.frame(
  gene = geneList,
  mu0 = mu0List,
  muV = muVList,
  se0 = se0List,
  seV = seVList
  )

df <- subset(df,mu0>0)
dfmelt <- melt(df,id.var=c('gene','se0','seV'))

dfmelt$error <- NA
dfmelt$error[dfmelt$variable=='mu0'] <- dfmelt$se0[dfmelt$variable=='mu0']
dfmelt$error[dfmelt$variable=='muV'] <- dfmelt$seV[dfmelt$variable=='muV']
dfmelt <- dfmelt[,c('gene','variable','value','error')]

dfmelt$ystart <- NA
dfmelt$ystart[dfmelt$variable=='mu0'] <- dfmelt$value[dfmelt$variable=='mu0']
dfmelt$ystart[dfmelt$variable=='muV'] <- dfmelt$value[dfmelt$variable=='mu0']

dfmelt$yend <- dfmelt$ystart+dfmelt$error
dfmelt$yend[dfmelt$variable=='muV'] <- dfmelt$ystart[dfmelt$variable=='muV']-dfmelt$error[dfmelt$variable=='muV']
dfmelt$yend[dfmelt$yend<0] <- 0

dfmelt <- dfmelt[dfmelt$value>0 & dfmelt$value < 1,]

p <- ggplot(dfmelt,aes(x=gene,y=value,fill=variable))
#pdf('~/Dropbox/densitiyfigures/Figure2/Vol_Dep_Indep.pdf',width=3.2,height=2)
#pdf('~/Dropbox/LabMeeting140124/mu_crl.pdf',width=10,height=7)
p + geom_bar(stat='identity') + 
  #coord_flip() +
  geom_segment(aes(xend = gene,y = ystart,yend = yend)) + 
  geom_abline(slope = 0, intercept = 0.5,linetype = 'dashed',col='slategray4') +
  geom_point(aes(x = gene,y = yend),shape = 95,show_guide = FALSE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title=element_text(size="6"), axis.text=element_text(size='6'),
        legend.position = "none") +
  #theme_bw() + theme_classic() +
  xlab('Gene Name') + ylab('Fraction of total transcription') +
  scale_y_continuous(limits = c(0,1)) +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black')) +
  scale_fill_manual(values=c('mediumpurple3','lightskyblue3'))
#dev.off()