library(vcfR)
library(SNPRelate)
library(gdsfmt)
library(SeqArray)
#library(ggplot2)
#library(ggrepel)
#library(pheatmap)
#library(ggpubr)
#library(gridExtra)
#library(adegenet)
#library(RColorBrewer)
#library(reshape2)
#library(knitr)
#library(kableExtra)
#library(data.table)
#library(dplyr)
library(plyr)
library(Biostrings)
library(gtools)

homedir="/Users/Nenad"
projectdir=paste0(homedir,"/projects/dingo")
#projectdir=paste0(homedir,"/projects/Anthony")
scriptsdir=paste0(projectdir,"/scripts")
tabledir=paste0(projectdir,"/tables")
system(paste0("mkdir -p ",tabledir))
figdir=paste0(projectdir,"/figures")
system(paste0("mkdir -p ",figdir))
popfile<-paste0(projectdir,"/annotation/popmap2")
populationsdir<-paste0(projectdir,"/results/populations/")
annotationdir<-"../annotation"

#rm(list=ls())
missingness<-0.7
vcfFile<-paste0("../vcfs/final_dingo_dog_wolf_default_unique.3.dw.vcf.gz")
vcf<-read.vcfR(vcfFile)

head(vcf)

queryMETA(vcf)
gt <- extract.gt(vcf, IDtoRowNames = F)
fixed <- getFIX(vcf)
snps <- cbind(fixed[,1:5], gt)
head(snps)[,1:10]

snps.1 <- as.data.frame(as.matrix(snps))
snps.1$CHROM <- as.character(as.factor(snps.1$CHROM))
snps.1$POS <- as.character(as.factor(snps.1$POS))
snps.1$ID <- as.character(as.factor(snps.1$ID))
snps.1$identifier <- with(snps.1, paste(CHROM, POS, ID),sep=":")
snps.3 <- snps.1[,1:(ncol(snps.1)-1)]

######################3

read.depth <- extract.gt(vcf, element="AD")
length(unique(rownames(read.depth)))
nrow(read.depth)
read.depth.ref <- masplit(read.depth, record = 1, sort=0)
read.depth.snp <- masplit(read.depth, record = 2, sort=0)
#for reference allele coverage:
read.depth.ref.count<- rowSums(read.depth.ref, na.rm=T)
head(read.depth.ref.count)
read.depth.ref <- as.data.frame(read.depth.ref)
read.depth.ref$length <- rep(NA)
n <- ncol(read.depth.ref)-1
read.depth.ref[is.na(read.depth.ref)]<-0
read.depth.ref$length<-apply(read.depth.ref,1,function(x){sum(x[1:n]>0)})

read.depth.ref.avg <- read.depth.ref.count/read.depth.ref$length
read.depth.snp.count<- rowSums(read.depth.snp, na.rm=T)
head(read.depth.snp.count)
read.depth.snp <- as.data.frame(read.depth.snp)
read.depth.snp$length <- rep(NA)
n <- ncol(read.depth.snp)-1
read.depth.snp$length<-apply(read.depth.snp,1,function(x){sum(x[1:n]>0)})

read.depth.snp.avg <- read.depth.snp.count/read.depth.snp$length

coverage.rd <- cbind(read.depth.ref.avg, read.depth.snp.avg)
coverage.rd <- as.data.frame(coverage.rd)
coverage.rd$snp.index <- 1:nrow(coverage.rd)


index <- 1:nrow(snps.3)
snps.index <- cbind(index, snps.3)
snps.rd <- snps.index[which(snps.index[,1] %in% coverage.rd$snp.index),]
nrow(snps.rd)
callrate <- apply(snps.rd, 1, function(x) 100-(sum(is.na(x))/(ncol(snps.rd)-6))*100)



#We also need to filter by __coverage__. If the reference and SNP allele do not amplify at the same rate, 
#this may indicate potential bias (errors in calling). DArT PL includes AvgCountRef as the sum of the tag 
#read counts for all samples, divided by the number of samples with non-zero tag read counts for the reference 
#allele row, AvgCountSNP is the same but for the SNP allele row. 

#Coverage can be calculated as the absolute percentage difference between the AvgCountRef 
#and AvgCountSNP columns. vcf output gives Allele Depth, with coverage of reference, then 
#coverage of SNP allele e.g. AD = 4,3 means reference allele has coverage = 4, snp 
#allele coverage = 3. Count only non-zero reads (already filtered out anyway from above).

max_coverage <- pmax(read.depth.ref.avg, read.depth.snp.avg)
SNPcoverage <- ((abs(read.depth.ref.avg - read.depth.snp.avg))/(max_coverage))*100

#check if the format has to be switched to 49427264|F|0--18:A>G from 2975:12:+
final<-data.frame(AlleleID=snps.1$ID)
final$AvgCountRef<-read.depth.snp.avg

#DONE: CloneID: Unique identifier of the sequence tag   
#check if the format is always everything before the pipe sign
final$AlleleID<-gsub(":","\\|",final$AlleleID)
final$CloneID<-gsub("\\|.*","",final$AlleleID)

#not sure if this sequence is needed but check where it can be extracted from
fna<-readDNAStringSet("../seqs/catalog.default_unique.fa.gz")
seqnames<-names(fna)
annoSeqs<-strsplit(seqnames," ")
annoSeqsdf<-do.call("rbind",annoSeqs)
annoSeqsdf<-as.data.frame(annoSeqsdf)
annoSeqsdf$pos<-gsub("pos=","",annoSeqsdf$V2)
annoSeqsdf$seq<-as.character(fna)
#snps.1$identifier <- with(snps.1, paste(CHROM, POS, ID,sep=":"))

snps.1$SNPPos<-as.numeric(sapply(strsplit(snps.1$ID,":"),function(x){x[2]}))
snps.1$strand<-gsub(".*:","",snps.1$ID)

snps.1$startPos<-0
snps.1$startPos[grepl("+",snps.1$strand)]<-as.numeric(snps.1$POS[grepl("+",snps.1$strand)])-snps.1$SNPPos[grepl("+",snps.1$strand)]+1
snps.1$startPos[grepl("-",snps.1$strand)]<-as.numeric(snps.1$POS[grepl("-",snps.1$strand)])+snps.1$SNPPos[grepl("-",snps.1$strand)]-1

snps.1$identifier <- apply(snps.1,1,function(x){paste0(x["CHROM"],":",x["startPos"],":",x["strand"])})
snps.1$identifier<-gsub(" ","",snps.1$identifier)

#chrAnnotation
chrAnno<-read.table("../annotation/chrmap")
annoSeqsdf$chr<-gsub("\\:.*","",annoSeqsdf$pos)
annoSeqsdf$part2<-gsub(".*\\.3:","",annoSeqsdf$pos)
annoSeqsdf$chr2<-mapvalues(annoSeqsdf$chr,chrAnno$V1,chrAnno$V2)
annoSeqsdf$combined<-paste0(annoSeqsdf$chr2,":",annoSeqsdf$part2)
merged<-merge(snps.1,annoSeqsdf,by.x="identifier",by.y="combined",sort = F,all.x=T)

#
#merge(snps.1,annoSeqsdf,by.x="identifier",by.y="combined",sort = F)

final$AlleleSequence<-merged$seq
final$TrimmedSequence<-gsub("AGAT.*","",final$AlleleSequence)
final$TrimmedSequence[nchar(final$TrimmedSequence)<20]<-gsub("AGAT.*","",final$AlleleSequence[nchar(final$TrimmedSequence)<20])

#DONE: CallRate: Proportion of samples for which the genotype call is either "1" or "0", rather than "‐"    
final$CallRate<-callrate

#OneRatio: proportion of samples for which the genotype score is "1"
index <- 1:nrow(snps.3)
snps.index <- cbind(index, snps.3)
snps.rd <- snps.index[which(snps.index[,1] %in% coverage.rd$snp.index),]
snps.rdS<-snps.rd[,-(1:6)]
OneRatios<- apply(snps.rdS, 1, function(x){tabulate(factor(unlist(strsplit(x, "/", TRUE))))})
dfOneRatios<-apply(OneRatios,1,function(x){x/sum(x)})

#DONE: FreqHomRef	proportion of samples homozygous for the Reference allele 
homoHetcounts<- apply(snps.rdS, 1, function(x){tabulate(factor(x,levels=c("0/0","0/1","1/1")),nbins = 3)})
df<-t(apply(t(homoHetcounts),1,function(x){x/sum(x)}))
final$FreqHomRef<-df[,1]
#DONE: FreqHomSnp	
final$FreqHomSnp<-df[,2]





#DONE: Chrom
final$Chrom_DCan22_private_v1<-snps.1$CHROM
#DONE: Chrom Position
final$ChromPos_DCan22_private_v1<-snps.1$POS
#get read counts over SNP AlnCnt_DCan22_private_v1
final$AlnCnt_DCan22_private_v1<-""
#get alignment quality AlnEvalue_DCan22_private_v1
final$AlnEvalue_DCan22_private_v1<-""
#DONE: SNP
position<-sapply(strsplit(snps.1$ID,":"),function(x){x[2]})
final$SNP<-paste0(position,":",snps.1$REF,">",snps.1$ALT)
#DONE: SnpPosition
final$SnpPosition<-position
#DONE: CallRate: Proportion of samples for which the genotype call is either "1" or "0", rather than "‐"    
final$CallRate<-callrate

#DONE:OneRatioRef Proportion of (non -?) samples for which the genotype score is "0" in VCF (or 1 in final dart, first row for each SNP) 	
  #genotype encoded as alleles separated by /. Alleles are "0" for the reference allele, "1" for alternate, or "." if missing (GT). Heterozygote : 0/1

index <- 1:nrow(snps.3)
snps.index <- cbind(index, snps.3)
snps.rd <- snps.index[which(snps.index[,1] %in% coverage.rd$snp.index),]
snps.rdS<-snps.rd[,-(1:6)]
OneRatios<- apply(snps.rdS, 1, function(x){tabulate(factor(unlist(strsplit(x, "/", TRUE))))})

#DONE: OneRatioSnpProportion of (non -?) samples for which the genotype score is "1" in VCF (or 1 in final dart, second row for each SNP ) ) 	
df<-apply(t(OneRatios),1,function(x){x/sum(x)})
final$OneRatioRef<-df[1,]
final$OneRatioSnp<-df[2,]


#DONE: FreqHomRef	proportion of samples homozygous for the Reference allele 
homoHetcounts<- apply(snps.rdS, 1, function(x){tabulate(factor(x,levels=c("0/0","0/1","1/1")),nbins = 3)})
df<-t(apply(t(homoHetcounts),1,function(x){x/sum(x)}))
final$FreqHomRef<-df[,1]

#DONE: FreqHomSnp	
final$FreqHomSnp<-df[,2]

#DONE: FreqHets	proportion of samples which score as heterozygous, that is, scored as 1 
final$FreqHets<-df[,3]

#DONE: PICRef  Polymorphism Information Content extracted from the dartR script utils.recalc.avgpic


OneRatioRef <- final$OneRatioRef
OneRatioSnp <- final$OneRatioSnp
ZeroRatioRef <- 1 - OneRatioRef
ZeroRatioSnp <- 1 - OneRatioSnp
final$PICRef <-
  1 - ((OneRatioRef * OneRatioRef) + (ZeroRatioRef * ZeroRatioRef))

#DONE: PICSnp
final$PICSnp <-
  1 - ((OneRatioSnp * OneRatioSnp) + (ZeroRatioSnp * ZeroRatioSnp))

#DONE:AvgPIC	 average of the polymorphism information content (PIC) of the Reference and SNP alleles 
final$AvgPIC <-
  c(final$PICRef + final$PICSnp) / 2

#DONE:AvgCountRef	sum of the tag read counts for all samples, divided by the number of samples 
    #with non‐zero tag read counts, for the Reference allele row

final$AvgCountRef<-read.depth.ref.avg

#DONE: AvgCountSnp	sum of the tag read counts for all samples, divided by the number of samples 
    #with non‐zero tag read counts, for the Alternate (SNP) allele row
final$AvgCountSnp<-read.depth.snp.avg

################# rep avg

#now fix the rep avg

#import vcf file
vcf.fn<-"../vcfs/final_dingo_dog_wolf_default.3.dw.vcf.gz"
snpgdsVCF2GDS(vcf.fn, "ccm.gds")
genofile <- snpgdsOpen("ccm.gds")

#get genotypes and sample IDs
g <- read.gdsn(index.gdsn(genofile, "genotype"))
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

#load in the annotation to get replicate data from "target" files from dartseq
#files<-list.files(annotationdir,pattern="targets",full.names=T)
#results<-list()
#for(file in files){
#  results[[file]]<-read.csv(file)
#}
#targets<-do.call("rbind",results)
data<-read.table(paste0("../annotation/DCan22combined_annotation.txt"),sep="\t")

#clean targets
#targets$id<-paste0("a",targets$targetid)
#targets$cleanIDs<-targets$genotype
#targets$cleanIDs<-gsub("_rep.*","",targets$cleanIDs)
#merged<-merge(targets,data,by.x="cleanIDs",by.y="Sample.ID")

#clean genotypes
gS<-as.data.frame(g)
gS$id<-sample.id
gS<-gS[sample.id %in% as.character(data$targetid),]

data$genotypeNew<-data$genotype
data$genotypeNew<-gsub("_rep.*","",data$genotypeNew)
data$genotypeNew<-gsub("rep.*","",data$genotypeNew)
data$genotypeNew<-gsub("dog.*","",data$genotypeNew)

mergeG<-merge(gS,data[,c("genotypeNew","targetid")],by.x="id",by.y="targetid")


#find the proportion of technical pairs for which the marker score is consistent
#turn to matrix for faster processing
mergeGMat<-as.matrix(mergeG[,-c(1,dim(mergeG)[2])])
ids<-mergeG[,"genotypeNew"]

#split matrix by ID
mergeGL<-lapply(split(seq_along(ids), ids), #split indices by a
                function(m, ind) m[ind,], m = mergeGMat)[order(unique(ids))]

#extract only duplicates

#get only duplicate data
isDuplicate<-sapply(mergeGL,function(x){class(x)[1]!="integer"})
mergedDups<-mergeGL[isDuplicate]
consistentCount<-sapply(mergedDups,function(x){apply(x,2,function(x){c(abs(max(x)-min(x)))==0})})

#there were 85 duplicates, calculate sum of identical elements and divide by 85
RepAvg<-rowSums(consistentCount)/dim(consistentCount)[2]
snp.chromosome<-read.gdsn(index.gdsn(genofile, "snp.chromosome"))
snp.position<-read.gdsn(index.gdsn(genofile, "snp.position"))
snp.id<-read.gdsn(index.gdsn(genofile, "snp.rs.id"))

df<-data.frame(id=paste0("chr",snp.chromosome,snp.position),RepAvg=RepAvg)




final<-final[paste0(final$Chrom_DCan22_private_v1,final$ChromPos_DCan22_private_v1) %in% df$id,]
#RepAvg proportion of technical replicate assay pairs for which the marker score is consistent
final$RepAvg<-mapvalues(paste0(final$Chrom_DCan22_private_v1,final$ChromPos_DCan22_private_v1),df$id,df$RepAvg)

write.csv(final,"DCann_dogwolf.csv")
#final<-final[!is.na(final$RepAvg),]
#final<-final[final$RepAvg>0,]














