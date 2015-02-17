library(boot)
library(plyr)

data <- read.delim('~/Dropbox/ExtractedData_131216/alldata_CRL.txt',header=T,stringsAsFactors=F)

data <- data[data$numCyclin<20,]

data <- ddply(data,.(gene),subset,length(cytoRNA)>20)

noise <- ddply(data,~gene,summarize,
               cv = sd(cytoRNA)/mean(cytoRNA),
               cov = cov(cytoRNA,volume),
               meanM = mean(cytoRNA),
               meanV = mean(volume))

coefs <- ddply(data, .(gene), function(df) {
  m <- lm(cytoRNA ~ volume, data=df)
  data.frame(intercept = coef(m)[1], slope = coef(m)[2])
})

noiseDat <- merge(noise, coefs)

noiseDat$s <- noiseDat$slope*noiseDat$meanV/(noiseDat$intercept+noiseDat$slope*noiseDat$meanV)

noiseMeasure <- ddply(noiseDat,.(gene),summarize,
                      noise = cv^2 - s*cov/(meanM*meanV))