
#This reads the vcf files stored as names in samples files an makes a corrending data object with the name provided
## > sample.files
##                              snp                            indel 
##   "SKDP-FAM-26_All_snps.raw.vcf" "SKDP-FAM-26_All_DINDEL.raw.vcf"

library(foreach)  ## could not get this to work well
library(doMC)
num.cores<-7
registerDoMC(cores=num.cores)

for(i in 1:length(sample.files)){
################ LARGE FILES ########
  print(paste("Doing samples",sample.files[i],sep=" "))
setwd(eval(as.name(paste(names(sample.files)[i],"dir",sep=".")))) ## different directory for each indel type

######BELOW PROCESSING  this for snp for indel varient types in vcf4.0 or vcf 3.0 format

### get location of header:
chromo1<-try(scan(sample.files[i],what=character(),n=5000,sep="\n",skip=0,fill=TRUE)) ## find the start of the vcf file
skip.lines<-grep("^#CHROM",chromo1)
if(length(skip.lines)>1){print("ERROR multiple chrom lables found");skip.lines<-skip.lines[1]}
skip.lines<-skip.lines
 
options(show.error.messages = TRUE)
column.labels<-read.delim(sample.files[i],header=F,nrows=1,skip=(skip.lines-1),sep="\t",fill=TRUE,stringsAsFactors=FALSE)
num.vars<-dim(column.labels)[2]

vcf.head<-read.delim(sample.files[i],header=F,nrows=(skip.lines-1),skip=0,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
info.flag.string<-"^##INFO=<ID="
info.types<-grep(info.flag.string,vcf.head[,1])
info.class<-extract.value.from.format(vcf.head[info.types,1],"Type=")
info.description<-extract.value.from.format(vcf.head[info.types,1],"Description=")  
info.types<-extract.value.from.format(vcf.head[info.types,1],"ID=")


format.flag.string<-"^##FORMAT=<ID="
format.types<-grep(format.flag.string,vcf.head[,1])
format.class<-extract.value.from.format(vcf.head[format.types,1],"Type=")
format.description<-extract.value.from.format(vcf.head[format.types,1],"Description=")  
format.types<-extract.value.from.format(vcf.head[format.types,1],"ID=")  

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
##   library(VariantAnnotation)
##          fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
     
##        ## Subset on genome coordinates :
##        ## 'file' must have a Tabix index
##        rngs <- GRanges("20", IRanges(c(14370, 1110000), c(17330, 1234600)))
##        names(rngs) <- c("geneA", "geneB")
##        param <- ScanVcfParam(which=rngs,asGRanges=TRUE) 
##        compressVcf <- bgzip(fl, tempfile())
##        idx <- indexTabix(compressVcf, "vcf")
##        tab <- TabixFile(compressVcf, idx)
##        vcf <- readVcf(tab, "hg19", param)

## param <- ScanVcfParam(which = GRanges("chr1", IRanges(1, 1e8)),asGRanges=TRUE,fixed = c("ALT", "FILTER"),info = "DP",geno = "GT")
## param <- ScanVcfParam(which = GRanges("chr1", IRanges(1, 1e8)),asGRanges=TRUE,fixed = c("ALT", "FILTER"),info = "DP",geno = "GT")
## param <- ScanVcfParam(fixed="ALT", geno=c("GT", "HQ"), info=c("NS", "AF"))
## test<-readVcf("75_0080_All_snps.raw.vcf","hg19",param)

  ############# look for file corruption and find bad location
## test<- read.delim(sample.files[i],header=F,nrows=(72850),skip=skip.lines,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
## colnames(test)<-column.labels
## tapply(test[,1],test[,1],length)
## test[grep(FALSE,grepl("^chr",test[,1]))[1:5],1:10]
##   grep "91963371" 75_0080_All_VARIANTS.raw.vcf

 
indels<-try(scan(sample.files[i],what=character(num.vars),skip=skip.lines,sep="\t",fill=TRUE))
num.lines<-length(indels)/(num.vars)
dim(indels)<-c(num.vars,num.lines)
indels<-t(indels)

column.labels[match("#CHROM",column.labels)]<-"chr"
colnames(indels)<-column.labels 
## indels[1:5,]

####################################### FINISHED Read in data
indels<- process.GATK.indels(indels,vcf.type,format.types,info.types,info.labels,info.class,num.cores)
  
assign(sample.grs[i],value=indels)
  
## data.gr
print(paste("Done",sample.grs[i],sep=" "))
}




## require(data.table)

## computeAllPairSums <- function(filename, nbindiv,nrows.to.read)
## {
##    con <- file(filename, open="r")
##    on.exit(close(con))
##    ans <- matrix(numeric(nbindiv * nbindiv), nrow=nbindiv)
##    chunk <- 0L
##    while (TRUE) {
##        #read.table faster than scan
##        df0 <- read.table(con,col.names=c("ID1", "ID2", "ignored", "sharing"),
##                 colClasses=c("integer", "integer", "NULL",
## "numeric"),nrows=nrows.to.read,comment.char="")

##        DT <- data.table(df0)
##        setkey(DT,ID1,ID2)
##        ss <- DT[,sum(sharing),by="ID1,ID2"]

##        if (nrow(df0) == 0L)
##            break

##        chunk <- chunk + 1L
##        cat("Processing chunk", chunk, "... ")

##       idd <- as.matrix(subset(ss,select=1:2))
##       newvec <- as.vector(as.matrix(subset(ss,select=3)))
##       ans[idd] <- ans[idd] + newvec

##          cat("OK\n")
##      }
##    ans
##  }

