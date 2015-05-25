library("stringr")
library("plyr")
library("IRanges")

#pwd
setwd('~/Documents/Data/SNP stuff/')

#files
gtfAttr <- read.delim('refseq_genes_for_snp.txt', stringsAsFactors=F, header=T)
snpAttr<- read.delim('GM12878_Homozygous_NonReference_SNPs.txt', stringsAsFactors=F, header=T)

#find unique set of chromosomes included in each list
olistchrms <- unique(gtfAttr$chrom)
reflistchrms <- unique(snpAttr$chrom)

#placeholder variables for later use
ochrlist<-gtfAttr[gtfAttr$chrom=='chr1',]
refchrlist<-snpAttr[snpAttr$chrom=='chr1',]

#set mastertable for later use.
mastertable <- list()

#chromosome by chromosome, find all features in olivia's list on that chromosome
#and all SNPs on that chromosome in ref list. Then, set ranges as defined by starts
#and ends in the files. Then, cycle through snp list for given chromosome and find
#any overlaps on corresponding olivia feature list. Lastly, print to output file in new
#line all information from ref and olivia files on snp and feature, one line per snp.
for (ochr in reflistchrms){  

  ochrlist<-gtfAttr[gtfAttr$chrom==ochr,]
  refchrlist<-snpAttr[snpAttr$chrom==ochr,]
  
  ofeatsranges<-IRanges(start = ochrlist$start, end = ochrlist$end)
  refsnpranges<-IRanges(start = refchrlist$chromStart, end = refchrlist$chromEnd)
  
  overlaps<-as.matrix(findOverlaps(refsnpranges,ofeatsranges))
  
  chrhitstable<-cbind(ochrlist[,c(9,10,3,1,4,5,7)][overlaps[,2],], refchrlist[,c(4,3)][overlaps[,1],])
  
  mastertable[[ochr]] <- chrhitstable

}

finaltable<-do.call(rbind, mastertable)

colnames(finaltable)<-c('TranscriptID','GeneName','FeatureType','Chr','FeatureStart','FeatureEnd','FeatureStrand','SNP_ID','SNP_position')

outFile<-'genome_snp_aligns.txt'
write.table(finaltable,outFile,quote=F,sep='\t',row.names=F,col.names=T)
#end snp-feature alignment


#find genes of a certain abundance & halflife

seqAbund <- read.delim('~/Dropbox/NucSeqData/dataNewGtf/CRL_total1_FPKM.txt', 
                       stringsAsFactors=F, header=T)
colnames(seqAbund) <- c('GeneName', 'ExRPKM', 'IntRPKM')
seqAbund <- seqAbund[,1:2]

intAbund <- read.delim('~/Dropbox/NucSeqData/dataNewGtf/CRL_nuc1_FPKM.txt', 
                       stringsAsFactors=F, header=T)
colnames(intAbund) <- c('GeneName', 'ExRPKM', 'IntRPKM')
intAbund <- intAbund[,c(1,3)]

hlTable <- read.delim('~/Dropbox/Foreskin/half_life_rpkm_with_names.txt',
                      stringsAsFactors=F, header=T)
colnames(hlTable) <- c('GeneName','TxID','halflife','FPKM')
hlTable <- hlTable[,c(1,3)]
hlTable <- hlTable[hlTable$halflife!='N.D.',]
hlTable$halflife <- gsub(' ', '', hlTable$halflife)
hlTable$halflife <- as.numeric(hlTable$halflife)

abundTable <- merge(finaltable,seqAbund)
fpkmTable <- merge(abundTable,intAbund)
totalTable <- merge(fpkmTable,hlTable)
write.table(abundTable,'genome_snp_aligns_ex_int_hl.txt',quote=F,sep='\t',row.names=F,col.names=T)

tmp <- subset(totalTable,ExRPKM>100)
tmp <- subset(tmp,halflife>8)
length(unique(tmp$GeneName))

write.table(tmp,'candidates.txt',quote=F,sep='\t',row.names=F,col.names=T)

write.table(unique(tmp$SNP_ID),'candidate_snps.txt',quote=F,sep='\t',row.names=F,col.names=F)

#start probe design, GM first, hg19 second

#reading in maternal and paternal GM sequences, but only using snps homozygous over the two, so probes can
#be created from either.
targetchrs<-c('Chrx','Chry','Chrz') #set based on manually chosen target genes

#gm probe loop
for (chrm in targetchrs){
  
  GMinfileMat <- paste('GM12878/NA12878_diploid_genome_dec16_2013/', chrm , '_NA12878_maternal.fa', sep ="")
  GMinfilePat <- paste('GM12878/NA12878_diploid_genome_dec16_2013/', chrm , '_NA12878_paternal.fa', sep ="")
  GMinfileMap <- paste('GM12878/NA12878_diploid_genome_dec16_2013/', chrm , '_NA12878.map', sep ="")
  GMinfileSNP <- 'GM12878_HetSNPs_phased_pos_v2.txt' #need HOMOZYGOUS snps
  GMinFileGene <- 'Chrs_SKA3.txt' #rename
  GMoutFile <- 'Chrs_SKA3_snps.txt' #rename

  GMgeneList <- read.delim(GMinFileGene, stringsAsFactors=F, header=F)

  gmMAT <- read.delim(GMinfileMat, header = F, stringsAsFactors=F, colClasses = "character", col.names = "seq")
  gmPAT <- read.delim(GMinfilePat, header = F, stringsAsFactors=F, colClasses = "character", col.names = "seq")
  gmMAP <- read.table(GMinfileMap, header = F, stringsAsFactors=F, colClasses= "integer", 
                  col.names = c("REF","PAT","MAT") )
  gmSNPs <- read.table(GMinfileSNP, header = T, stringsAsFactors=F, 
                   colClasses=c("character", rep("integer", 2), "character") )

  gtfAttr <- read.delim('refSeq_refGene_12-13-12_hg19_with_group_columns.proteincoding.introns.txt', stringsAsFactors=F, header=T)
  chrgtfAttr <- gtfAttr[gtfAttr$seqname == chrm, ] #Gene model with introns for Chr

  gmMAT <- DNAString(str_c(gmMAT$seq[2:length(gmMAT$seq)], collapse =""))
  gmPAT <- DNAString(str_c(gmPAT$seq[2:length(gmPAT$seq)], collapse =""))

  GMchrSNPs <- gmSNPs[gmSNPs$chrom == chrm,]

}

gmfinaltable<-do.call(rbind, GMmastertable)

names(gmfinaltable)<-c('x','y','z')

gmoutfile<-'GMprobelist.txt'
write.table(gmfinaltable,gmoutFile,quote=F,sep='\t',row.names=F,col.names=T)