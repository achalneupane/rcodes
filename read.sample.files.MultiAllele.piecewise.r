
#This reads the vcf files stored as names in samples files an makes a corrending data object with the name provided
## > sample.files<-"2015-01-15_LeoPharma_Dec2014Freeze.chr21.output.recalibrated.filtered.vcf"

##     snps.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-01-15_LeoPharma_Dec2014Freeze/SNPs"                        snp                            indel 
##   "SKDP-FAM-26_All_snps.raw.vcf" "SKDP-FAM-26_All_DINDEL.raw.vcf"

library(foreach)  ## could not get this to work well
library(doMC)
if(!exists("num.cores")){num.cores<-4}
registerDoMC(cores=num.cores)

isamp<-1
for(isamp in 1:length(sample.files)){
################ LARGE FILES ########
print(paste("Doing samples",sample.files[isamp],sep=" "))
setwd(eval(as.name(paste(names(sample.files)[isamp],"dir",sep=".")))) ## different directory for each indel type

######BELOW PROCESSING  this for snp for indel varient types in vcf4.0 or vcf 3.0 format
  
################## get initial read information
start.data<-prepare.for.Vcf.file.read(sample.files[isamp])
for(i in 1:length(start.data)){assign( names(start.data)[i],value=start.data[[i]])}
###########################################
  
## ### get location of header:
## #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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

  
setwd(eval(as.name(paste(names(sample.files)[isamp],"dir",sep="."))))
con <- file(sample.files[isamp], open="r")  # close(con)
num.lines<-1 # so does while llop at least once
reads<-max.reads  #1.3M lines in snp file 50000 goes out to 24Gb without QC cgeck 
counter<- -1
indels.all<-{}
while (num.lines >0  ){
counter<-counter+1
print(counter)

if(counter==0){
indels<-try(scan(con,what=character(num.vars),skip=(reads*counter)+skip.lines,nlines=reads,sep="\t",fill=TRUE,na.strings="",quote="\""))
}else{
indels<-try(scan(con,what=character(num.vars),nlines=reads,sep="\t",fill=TRUE,na.strings="",quote="\""))
}
  

num.lines<-length(indels)/(num.vars)
dim(indels)<-c(num.vars,num.lines)
indels<-t(indels)
colnames(indels)<-column.labels 
print(dim(indels))
## indels[1:5,1:5]
if(num.lines==0){next}

####################################### FINISHED Read in data
indels<- process.GATK.indels(indels,sample.files[isamp],vcf.type,format.types,info.types,info.class,num.cores)

if(is.null(dim(indels.all))){
  if(length(remove.extensions)>1){
    for(ir in 1:length(remove.extensions)){
      remove.cols<-c(remove.cols,c(colnames(indels)[grepl(paste(remove.extensions[ir],"$",sep=""),colnames(indels))]))
                               }}
    cols.keep<-colnames(indels)[!(colnames(indels) %in% remove.cols)]
    indels.all<-indels[,cols.keep]
                     }else{
                       indels.all<-rbind(indels.all,indels[,cols.keep])
                     }
  

} ## loop over data chunks

close(con)
  
assign(sample.grs[isamp],value=indels.all)
print(dim(indels.all))  
## data.gr
print(paste("Done",sample.grs[isamp],sep=" "))
}


## remove.cols<-c(colnames(indels)[grepl(".PL$",colnames(indels))]

## require(data.table
