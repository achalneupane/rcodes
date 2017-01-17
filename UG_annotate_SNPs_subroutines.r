
############################START  FUNCTIONS USED TO READ AND PROCESS VCF
  #######--prototype to use instead of scan
## filena:me<-sample.files[i]
##   nrows.to.read<-10000
## fast.read <- function(filename,nrows.to.read, column.labels,skip.lines )
## {
##    con <- file(filename, open="r")
##    on.exit(close(con))
##    ans <- matrix(numeric(nbindiv * nbindiv), nrow=nbindiv)
##    chunk <- 0L
##    while (TRUE) {
##        #read.table faster than scan
## #     system.time(
##        indels1 <- read.table(con,col.names=column.labels,skip=skip.lines,fill=TRUE,stringsAsFactors=FALSE,
##                 colClasses="character",nrows=nrows.to.read,comment.char="",quote="")
## #)

##        if (nrow(indels) == 0L)
##            break


##      }
##    ans
##  }
 #######-----------------------
## the.indels<-g.indel
## root<-"test"
##################################### AOGC beta rescale
#res <- logscan(genotypes, a.fam,"pheno",num.samples)
# logistic.run(genotypes[41458:41460,], fam,"pheno",num.samples)

 ## logistic.run(genotypes, fam,"pheno",num.samples)

#genotypes, fam,"pheno",num.samples)
logistic.run<- function(gdata, pdata,pcolumn,num.samp) {
 ## gdata <- genotypes
 ## pdata <- fam
 ## pcolumn <- "pheno"
 ## num.samp <- num.samples

  
  ## print(dim(gdata))
  ## print(rownames(gdata))
  n <- dim(gdata)[1]
##   if(is.null(n)){ # only one element
##     gdata<-as.matrix(gdata)
##     if(dim(gdata)[1]!=1){gdata<-t(gdata)}
## }
#  n <- names(gdata)
  markers <- length(n) - 1

  all.missing<-apply(gdata,1,function(x) sum(is.na(as.numeric(x))))
  all.missing<-grep(dim(gdata)[2],all.missing,fixed=TRUE)
  do.rows<-1:dim(gdata)[1]
  if(length(all.missing)>0){do.rows<-do.rows[!(do.rows %in% all.missing)]}
  ## all.missing[1:5]
  ## sum(missing)

if(is.null(rownames(gdata))){rownames(gdata)<-1:dim(gdata)[1]}
  
  results <- matrix(data=NA,nrow=nrow(gdata),ncol=7)
  #results <- matrix(nrow=ncol(gdata),ncol=5)
  colnames(results) <- c("SNP","REFfreq","ngeno","Estimate", "Std Error", "Z", "pval" )
#  i<-2
  for( i in 1:length(do.rows)) {
    print(i)
    
 #   for( i in 2:length(n)) {  
   # pdata[,"missing"]<-as.numeric(gdata[i,])
    ## for aogc & hbm use both
    #  fit1 <- glm(form1,data=data,family=binomial(logit)) # logistic for case/control
  #  fit <- lm(gdata[i,] ~ pdata[,"PCA1"]+ pdata[,"PCA2"] + pdata[,"PCA3"] + pdata[,"PCA4"]+pdata[,pcolumn])
   fit <- lm(as.numeric(pdata[,pcolumn])~as.numeric(gdata[do.rows[i],]) ) # same results as plink
 #  fit <- lm(as.numeric(gdata[do.rows[i],]) ~ as.numeric(pdata[,pcolumn]) )
    a <- summary(fit)
    
    ### get coeff for pheno1 only
#     Estimate <- a$coefficients[2,1]
#     Std_error <- a$coefficients[2,2]
#     T <- a$coefficients[2,3]
#     pval <- a$coefficients[6,1:4]
      the.geno<-as.numeric(gdata[do.rows[i],])
      the.geno<-the.geno[!is.na(the.geno)]
      p= sum(the.geno)/(2*length(the.geno))
  # p= 1- sum(as.numeric(gdata[do.rows[i],]),na.rm=TRUE)/(2*num.samp)

    if( dim(a$coefficients)[1]==1){ # no slope (monomorphic typically)  so get P if intercept - bad! as dim(a$coefficients)[1]==1 NOT 2 (slope and intercept)
      results[do.rows[i],] <- c(rownames(gdata)[do.rows[i]],p,length(the.geno),rep(NA,times=4))
      }else{
    # results[do.rows[i],] <- c(rownames(gdata)[do.rows[i]],p,length(the.geno),a$coefficients[dim(a$coefficients)[1],1:4])    #a$coefficients[c(Estimate, Std_error, T, pval )
     results[do.rows[i],] <- c(rownames(gdata)[do.rows[i]],p,length(the.geno),a$coefficients["as.numeric(gdata[do.rows[i], ])",1:4]) ## just in case ordr changes
   
     }
  
  }

  if(length(all.missing)>0){results[all.missing,"SNP"]<-rownames(gdata)[all.missing]}
  #rownames(results) <- seq(from = 1, to = markers)
  #colnames(results) <- c("Estimate", "Std Error", "Z", "pval" )
 
  return(results)
  
}   






logscan<- function(gdata, pdata,pcolumn,num.samp) {
 ## gdata <- genotypes[c(41,54,76,77,81),]# [1,]
 ## pdata <- a.fam
 ## pcolumn <- "pheno"
 ## num.samp <- num.samples

  
  ## print(dim(gdata))
  ## print(rownames(gdata))
  n <- dim(gdata)[1]
##   if(is.null(n)){ # only one element
##     gdata<-as.matrix(gdata)
##     if(dim(gdata)[1]!=1){gdata<-t(gdata)}
## }
#  n <- names(gdata)
  markers <- length(n) - 1

  all.missing<-apply(gdata,1,function(x) sum(is.na(x)))
  all.missing<-grep(dim(gdata)[2],all.missing,fixed=TRUE)
  do.rows<-1:dim(gdata)[1]
  if(length(all.missing)>0){do.rows<-do.rows[!(do.rows %in% all.missing)]}
  ## all.missing[1:5]
  ## sum(missing)

if(is.null(rownames(gdata))){rownames(gdata)<-1:dim(gdata)[1]}
  
  results <- matrix(data=NA,nrow=nrow(gdata),ncol=7)
  #results <- matrix(nrow=ncol(gdata),ncol=5)
  colnames(results) <- c("SNP","REFfreq","ngeno","Estimate", "Std Error", "Z", "pval" )
#  i<- 1
  for( i in 1:length(do.rows)) {
 #   print(i)
    

  #  plot(as.numeric(gdata[do.rows[i],]), as.numeric(pdata[,pcolumn])) 
  #  fit <- lm(gdata[i,] ~ pdata[,"PCA1"]+ pdata[,"PCA2"] + pdata[,"PCA3"] + pdata[,"PCA4"]+pdata[,pcolumn])
 #  fit <- lm(as.numeric(pdata[,pcolumn])~as.numeric(gdata[do.rows[i],]) ) # same results as plink
   fit <- lm(as.numeric(gdata[do.rows[i],]) ~ as.numeric(pdata[,pcolumn]) )
    a <- summary(fit)
    
    ### get coeff for pheno1 only
#     Estimate <- a$coefficients[2,1]
#     Std_error <- a$coefficients[2,2]
#     T <- a$coefficients[2,3]
#     pval <- a$coefficients[6,1:4]
      the.geno<-as.numeric(gdata[do.rows[i],])
      the.geno<-the.geno[!is.na(the.geno)]
      p= sum(the.geno)/(2*length(the.geno))
    
  # p= 1- sum(as.numeric(gdata[do.rows[i],]),na.rm=TRUE)/(2*num.samp)
     if( dim(a$coefficients)[1]==1){ # no slope (monomorphic typically)  so get P if intercept - bad! as dim(a$coefficients)[1]==1 NOT 2 (slope and intercept)
      results[do.rows[i],] <- c(rownames(gdata)[do.rows[i]],p,length(the.geno),rep(NA,times=4))
      }else{
    # results[do.rows[i],] <- c(rownames(gdata)[do.rows[i]],p,length(the.geno),a$coefficients[dim(a$coefficients)[1],1:4])    #a$coefficients[c(Estimate, Std_error, T, pval )
     results[do.rows[i],] <- c(rownames(gdata)[do.rows[i]],p,length(the.geno),a$coefficients["as.numeric(pdata[, pcolumn])",1:4]) ## just in case ordr changes
   
     }
    results[do.rows[i],] <- c(rownames(gdata)[do.rows[i]],p,length(the.geno),a$coefficients[dim(a$coefficients)[1],1:4])    #a$coefficients[c(Estimate, Std_error, T, pval )
  
  }

  if(length(all.missing)>0){results[all.missing,"SNP"]<-rownames(gdata)[all.missing]}
  #rownames(results) <- seq(from = 1, to = markers)
  #colnames(results) <- c("Estimate", "Std Error", "Z", "pval" )
 
  return(results)
  
}   




get.forward.strand.allele<-function(indels,i,sample.files,snp.dir){
# c("chr","start","end","REF","ALT") columns required
  ## i = loop counter form read vcf projram
  ## sample.files = array with names of the file being read (used to name output list ) - ele set to data/ name (GENO)
  ## i,sample.files,snp.dir are optional and are set below if they do bot exist

  ################## get forward starnd allele function inputs i sindels output is indels and to.flip
  ################## done this way so can do 

library("Biostrings")
library("IRanges")

if(!exists("genome.build")){genome.build<-"hg19"}
if(genome.build=="hg19"){library("BSgenome.Hsapiens.UCSC.hg19"); have.chr<-seqnames(Hsapiens);chr.lengths<-seqlengths(Hsapiens)}
if(genome.build=="hg18"){library("BSgenome.Hsapiens.UCSC.hg18"); have.chr<-seqnames(Hsapiens);chr.lengths<-seqlengths(Hsapiens)}
if(genome.build=="mm9"){library("BSgenome.Mmusculus.UCSC.mm9"); have.chr<-seqnames(Mmusculus)}
if(genome.build=="mm10"){library("BSgenome.Mmusculus.UCSC.mm10"); have.chr<-seqnames(Mmusculus)}
######## CAREFUL this code used by read.plink.assoc.files.r
print(paste("Genome build is:",genome.build))
print("inside forward strand")
 print(dim(indels))
print("inside forward strand")

if(!grepl("^chr",indels[1,"chr"])){indels[,"chr"]<-paste("chr",indels[,"chr"],sep="")}

## print("test")
## print(!exists("sample.files"))
## print(!exists("sample.files",inherits=FALSE))
## print( !exists("sample.files",2,inherits=FALSE) )

if(!exists("sample.files",2,inherits=FALSE)){sample.files<-c("data");i<-1}
if(!exists("snp.dir",2,inherits=FALSE)){snp.dir<-getwd()}
#print(sample.files)
#print(sample.files[i])
## if(!grepl("^chr",indels[1,"chr"])){
## key.indels<-build.key(all,core.ann,add.chr.label=TRUE)
##   indels[,"chr"]<-paste("chr",indels[,"chr"],sep="")
## }else{key.indels<-build.key(indels,core.ann)      }
## indels[1:5,]
#print(class(indels))

## names(types)[grep("frameshift_variant",names(types))]
## types<-sort(tapply(poly[,"Consequence"],poly[,"Consequence"],length))
## to.get<-"Consequence"
## order.by<-vep.types
##  to.get.other<-c("Uploaded_variation","Gene","Feature","Protein_position","Amino_acids")
## to.get.other<-{}




#indels<-as.matrix(indels)  ## else get factor issues

#no.chr<-check.chr.and.positions(indels,chr.lengths)


  positions<-chr.lengths[indels[,"chr"]]
#as.integer(indels[,"start"])

indels[indels[,"chr"]=="chr23","chr"]<-"chrX"  ## no Y observed
indels[indels[,"chr"]=="chr24","chr"]<-"chrY"
indels[indels[,"chr"]=="chr25","chr"]<-"chrXY" ## chrXY is chr25 but coordinates/positions ARE on chrX (note no XY in UCSC build)
indels[indels[,"chr"]=="chr26","chr"]<-"chrM"
indels[indels[,"chr"]=="chrXY","chr"]<-"chrX"
## tapply(indels[,"chr"],indels[,"chr"],length)

no.chr<-check.chr.and.positions(indels,chr.lengths)
sum(no.chr)

if(sum(no.chr)>1){
print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
exclude.file<-paste("to.exclude",sample.files[i],"txt",sep=".")
print(paste("Missing chromosomes/off end of chromosomes  so won't be tested write to ",exclude.file,sep=""))
print(indels[no.chr,])
setwd(snp.dir)      
write.table(indels[no.chr,"SNP"],file=exclude.file,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
setwd(code.dir)
print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
}
## the.chrom[1:5]
## tapply(the.chrom,the.chrom,length)


    ##  X    X chromosome                    -> 23
    ##  Y    Y chromosome                    -> 24
    ##  XY   Pseudo-autosomal region of X    -> 25
    ##  MT   Mitochondrial                   -> 26
# Hsapiens


the.chrom<-as.character(indels[!no.chr,"chr"])
starts<-as.integer(as.character(indels[!no.chr,"start"]))
ends<-as.integer(as.character(indels[!no.chr,"end"]))
all.genomic.have<-getSeq(Hsapiens, the.chrom, starts, ends,as.character=TRUE) ## get forward stand alleles


all.genomic<-indels[,"REF"]
all.genomic[!no.chr]<-all.genomic.have #
all.genomic.ori<-all.genomic
## all.genomic[1:5]
## indels[,"REF.base"]<-
## alleles<-c("A","C","G","T")

ref.ori<-indels[,"REF"]
alt.ori<-indels[,"ALT"]

to.flip<-rep(0,times=length(ref.ori))
##TRY 
## (1) SWAP (+1 to to.flip)
## (2) Complement
## (3) SWAP  (-1 to to.flip)
## (4) Complent back so  no-fix's are back to original
## (5) assume what is left if forward strand alleles that are insertions
## to.flip !=0 need flipping
########### (1) ##########################################################################################
ref.diff<-indels[,"REF"]!=all.genomic & indels[,"REF"]!="-"  ## if a insertion can't tell if it match the genome
sum(ref.diff)

if(sum(ref.diff)!=0){
to.flip[ref.diff]<-to.flip[ref.diff]+1
## sum(ref.diff) #4601
print(paste(sum(ref.diff)," STEP1: strand differences or minor allele not REF"))
## cbind(indels,all.genomic)[ref.diff,][1:10,]

### swap REF AND ALT to see is just minor allele / major allele swap
the.ref<-indels[ref.diff,"REF"]
indels[ref.diff,"REF"]<-indels[ref.diff,"ALT"]
indels[ref.diff,"ALT"]<-the.ref
###### swap 0/1 --- keep track of times swapped

############ for indels must chenge start and end and the get new genomic sequence
an.indel<- (nchar(as.character(indels[ref.diff,"REF"])) >1) | (nchar(as.character(indels[ref.diff,"ALT"])) > 1) | indels[ref.diff,"REF"]=="-" | indels[ref.diff,"ALT"]=="-"
sum(an.indel)

if(sum(an.indel)>0){
  
  if(sum(ref.diff)==1){
temp<-  redo.start.end.annovar( subset(subset(indels,subset=ref.diff),subset=an.indel,select=c("chr","start","end","REF","ALT")) )
indels[ref.diff,"chr"][an.indel]<-temp[,"chr"]
indels[ref.diff,"start"][an.indel]<-temp[,"start"]
indels[ref.diff,"end"][an.indel]<-temp[,"end"]
indels[ref.diff,"REF"][an.indel]<-temp[,"REF"]
indels[ref.diff,"ALT"][an.indel]<-temp[,"ALT"]
}else{
  indels[ref.diff,][an.indel,c("chr","start","end","REF","ALT")]<-redo.start.end.annovar(subset(indels[ref.diff,],subset=an.indel,select=c("chr","start","end","REF","ALT")))
}

no.chr<-check.chr.and.positions(subset(subset(indels,subset=ref.diff),subset=an.indel,select=c("chr","start","end","REF","ALT")),chr.lengths)

the.chrom<-as.character(indels[ref.diff,"chr"][an.indel][!no.chr])
starts<-as.integer(as.character(indels[ref.diff,"start"][an.indel][!no.chr]))
ends<-as.integer(as.character(indels[ref.diff,"end"][an.indel][!no.chr]))
all.genomic.indels<-getSeq(Hsapiens, the.chrom, starts, ends,as.character=TRUE)
all.genomic[ref.diff][an.indel][!no.chr]<-all.genomic.indels
}
#cbind(all.genomic[!no.chr][ref.diff][an.indel],all.genomic.indels)
#### TOP /BOT diff
####################### issue here is that for INSERTIONS turned to DELETIONS may my chance match the reference mase cause insert/deletion is at a homopolymer run!
####################### so ca't tell IF it was more frequent or A1/ A2 or insetion or deletion


} # something to do in step 1

########### (2) ##########################################################################################

## ref.diff<-indels[,"REF"]!=all.genomic
## an.indel<- (nchar(as.character(indels[ref.diff,"REF"])) >1) | (nchar(as.character(indels[ref.diff,"ALT"])) > 1) | indels[ref.diff,"REF"]=="-" | indels[ref.diff,"ALT"]=="-"
## sum(ref.diff[an.indel]) # 6
## forward<-cbind(indels,all.genomic)[ref.diff,]
## save(list=c("forward"),file="forward.RData")
## cbind(forward[an.indel,],all.genomic.ori[ref.diff][an.indel],indels.ori[ref.diff,][an.indel,])

## cbind(forward,all.genomic.ori[ref.diff],indels[ref.diff,])[1:10,]
########################################################################################################
## IUPAC_CODE_MAP
   ##   A      C      G      T      M      R      W      S      Y      K      V      H      D      B      N 
   ## "A"    "C"    "G"    "T"   "AC"   "AG"   "AT"   "CG"   "CT"   "GT"  "ACG"  "ACT"  "AGT"  "CGT" "ACGT"

ref.diff<-indels[,"REF"]!=all.genomic
sum(ref.diff)


print(paste(sum(ref.diff),":STEP 2 Strand differences use complment DANGER with polymorphic alleles ASSUMING ONE IS THE REF ALLELE to flip written"))
print(paste("Flip file called:",paste("to.flip",sample.files[i],"txt",sep=".")," put in :",snp.dir))
 setwd(snp.dir)      
write.table(indels[ref.diff,"SNP"],file=paste("to.flip",sample.files[i],"txt",sep="."),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
setwd(code.dir)
## load("/media/scratch2/AOGC-NGS/ExomeChip/polymorphic.snps.RData") # save(list=c("poly.morphic.snps","all.data"),file="polymorphic.snps.RData")
## dim(all.data)
## cbind(forward,all.genomic.ori[ref.diff],indels[ref.diff,])[1:10,]
## posns<-match(indels[,"SNP"],poly.morphic.snps)
## poly.morphic.snps<-poly.morphic.snps[posns]
## cbind(forward,all.genomic.ori[ref.diff],indels[ref.diff,],poly.morphic.snps[ref.diff])[1:10,]
## test<- ref.diff & !is.na(poly.morphic.snps)
## sum(test)
## posns<-match(indels[,"SNP"],all.data[,"SNP"])
## all.data<-all.data[posns,]
## cbind(indels,all.genomic.ori,poly.morphic.snps,all.data)[test,]

if(sum(ref.diff)!=0){ ## something to do
  
######### this bad codes does not work for inserions or deletions 
#bad.codes<-(!(indels[ref.diff,"REF"] %in% names(IUPAC_CODE_MAP)) | !(indels[ref.diff,"ALT"] %in% names(IUPAC_CODE_MAP))) ## not a IUPAC id
## bad.codes2<-(nchar(indels[ref.diff,"REF"])==1 | nchar(indels[ref.diff,"REF"])==1 | indels[ref.diff,"REF"]!="-" | indels[ref.diff,"ALT"]!="-" ) ## An insert or del
## ## Can't handle indrts and del cause I have lost confortance with the real base, so only do Snps till this is fixed
## ##  [1,] "chr10" "61846"  "61846"  "+"    "snp" "chr10:61846"                        "C" "CTTT"          "-" "AAA"          "C"        
## ##  [2,] "chr10" "72170"  "72170"  "+"    "snp" "rs201275645"                        "C" "CT"            "-" "T"            "C"        
## ##  [3,] "chr10" "79178"  "79178"  "+"    "snp" "rs201215709"                        "A" "AT"            "-" "T"            "A"        
## ##  [4,] "chr10" "79179"  "79179"  "+"    "snp" "rs199708591"                        "T" "TG"            "-" "G"            "T"        
## ##  [5,] "chr10" "95429"  "95429"  "+"    "snp" "rs35442274"                         "C" "CA"            "-" "A"            "C"


## bad.codes<-bad.codes1 & !bad.codes2

## bad.codes.ref<-indels[ref.diff,"REF"]=="-"
## bad.codes.alt<-indels[ref.diff,"ALT"]=="-"
## bad.codes<-bad.codes.ref | bad.codes.alt
## sum(bad.codes)
## cbind(indels,all.genomic)[ref.diff,][bad.codes,]


#ref.diff[ref.diff] <- !bad.codes ## don't do ones that have bad codes


to.convert.ref<-as.character(indels[ref.diff,"REF"])
to.convert.alt<-as.character(indels[ref.diff,"ALT"])

one.base<-nchar(as.character(to.convert.ref))==1 &  nchar(as.character(to.convert.alt))==1
sum(one.base)
sum(ref.diff)

  
to.convert.ref[one.base]<-flip.one.base(to.convert.ref[one.base])       # handles "-"          
to.convert.alt[one.base]<-flip.one.base(to.convert.alt[one.base])      # handles "-"      

#cbind(forward[an.indel,],to.convert.ref[an.indel],to.convert.alt[an.indel],all.genomic.ori[ref.diff][an.indel],indels.ori[ref.diff,][an.indel,])

# icon<-1
if(length(to.convert.ref[!one.base])!=0){
for(icon in 1:length(to.convert.ref[!one.base])){
to.convert.ref[!one.base][icon]<-as.character(reverseComplement(DNAString(to.convert.ref[!one.base][icon])))
to.convert.alt[!one.base][icon]<-as.character(reverseComplement(DNAString(to.convert.alt[!one.base][icon])))
}
} ## more than one multibase



indels[ref.diff,"REF"]<-to.convert.ref
indels[ref.diff,"ALT"]<-to.convert.alt


} # somthing to do in step 2

########### (3) ##########################################################################################
#Migh have amajor / mino problem still 
ref.diff<-indels[,"REF"]!=all.genomic
sum(ref.diff)

if(sum(ref.diff)!=0){ ## something to do


to.flip[ref.diff]<-to.flip[ref.diff]-1

## sum(ref.diff) #4601
print(paste(sum(ref.diff),":STEP 3  minor allele not REF after complmenet swap these"))
## cbind(indels,all.genomic)[ref.diff,][1:10,]

### swap REF AND ALT to see is just minor allele / major allele swap
the.ref<-indels[ref.diff,"REF"]
indels[ref.diff,"REF"]<-indels[ref.diff,"ALT"]
indels[ref.diff,"ALT"]<-the.ref
###### swap 0/1 --- keep track of times swapped

############ for indels must chenge start and end and the get new genomic sequence
an.indel<- (nchar(as.character(indels[ref.diff,"REF"])) >1) | (nchar(as.character(indels[ref.diff,"ALT"])) > 1) | indels[ref.diff,"REF"]=="-" | indels[ref.diff,"ALT"]=="-"
sum(an.indel)





if(sum(an.indel)>0){
#

  ################# if ref.diff==1 then indels[ref.diff,] is a vector and indels[ref.diff,,drop=TRUE ] can't be used as input
if(sum(ref.diff)==1){
temp<-  redo.start.end.annovar( subset(subset(indels,subset=ref.diff),subset=an.indel,select=c("chr","start","end","REF","ALT")) )
indels[ref.diff,"chr"][an.indel]<-temp[,"chr"]
indels[ref.diff,"start"][an.indel]<-temp[,"start"]
indels[ref.diff,"end"][an.indel]<-temp[,"end"]
indels[ref.diff,"REF"][an.indel]<-temp[,"REF"]
indels[ref.diff,"ALT"][an.indel]<-temp[,"ALT"]
}else{
  indels[ref.diff,][an.indel,c("chr","start","end","REF","ALT")]<-redo.start.end.annovar(subset(indels[ref.diff,],subset=an.indel,select=c("chr","start","end","REF","ALT")))
}

no.chr<-check.chr.and.positions(subset(subset(indels,subset=ref.diff),subset=an.indel,select=c("chr","start","end","REF","ALT")),chr.lengths)
the.chrom<-as.character(indels[ref.diff,"chr"][an.indel][!no.chr])# this ok cause working with vectors from the start 
starts<-as.integer(as.character(indels[ref.diff,"start"][an.indel][!no.chr]))# this ok cause working with vectors from the start 
ends<-as.integer(as.character(indels[ref.diff,"end"][an.indel][!no.chr]))# this ok cause working with vectors from the start 
all.genomic.indels<-getSeq(Hsapiens, the.chrom, starts, ends,as.character=TRUE)
all.genomic[ref.diff][an.indel][!no.chr]<-all.genomic.indels # this ok cause working with vectors from the start 
}

#cbind(all.genomic[!no.chr][ref.diff][an.indel],all.genomic.indels)
#### TOP /BOT diff



} # something to do in step 3


## ref.diff<-indels[,"REF"]!=all.genomic
## sum(ref.diff)

## an.indel<- (nchar(as.character(indels[ref.diff,"REF"])) >1) | (nchar(as.character(indels[ref.diff,"ALT"])) > 1) | indels[ref.diff,"REF"]=="-" | indels[ref.diff,"ALT"]=="-"
## sum(ref.diff[an.indel]) # 21
## forward.remain<-cbind(indels,all.genomic)[ref.diff,]
## save(list=c("forward.remain"),file="forward.remain.RData")
## cbind(forward.remain[an.indel,],all.genomic.ori[ref.diff][an.indel],indels.ori[ref.diff,][an.indel,])
## cbind(forward.remain,all.genomic.ori[ref.diff],indels.ori[ref.diff,],indels.ori[ref.diff,])

ref.diff<-indels[,"REF"]!=all.genomic
sum(ref.diff)


########### (4) ##########################################################################################
if(sum(ref.diff)!=0){ ## something to do
  print(paste(sum(ref.diff),":STEP 4 cant be fixed after MINOR/COMPLEMENT/MINOR COMPLEMENT complement back to original state "))

print(paste("Un-fixable snp file called:",paste("un.fixable",sample.files[i],"txt",sep=".")," put in :",snp.dir))
setwd(snp.dir)      
write.table(indels[ref.diff,"SNP"],file=paste("un.fixable",sample.files[i],"txt",sep="."),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
setwd(code.dir)
  
######### this bad codes does not work for inserions or deletions 
#bad.codes<-(!(indels[ref.diff,"REF"] %in% names(IUPAC_CODE_MAP)) | !(indels[ref.diff,"ALT"] %in% names(IUPAC_CODE_MAP))) ## not a IUPAC id
## bad.codes2<-(nchar(indels[ref.diff,"REF"])==1 | nchar(indels[ref.diff,"REF"])==1 | indels[ref.diff,"REF"]!="-" | indels[ref.diff,"ALT"]!="-" ) ## An insert or del
## ## Can't handle indrts and del cause I have lost confortance with the real base, so only do Snps till this is fixed
## ##  [1,] "chr10" "61846"  "61846"  "+"    "snp" "chr10:61846"                        "C" "CTTT"          "-" "AAA"          "C"        
## ##  [2,] "chr10" "72170"  "72170"  "+"    "snp" "rs201275645"                        "C" "CT"            "-" "T"            "C"        
## bad.codes.ref<-indels[ref.diff,"REF"]=="-"
## bad.codes.alt<-indels[ref.diff,"ALT"]=="-"
## bad.codes<-bad.codes.ref | bad.codes.alt
## sum(bad.codes)
## cbind(indels,all.genomic)[ref.diff,][bad.codes,]


#ref.diff[ref.diff] <- !bad.codes ## don't do ones that have bad codes


to.convert.ref<-indels[ref.diff,"REF"]
to.convert.alt<-indels[ref.diff,"ALT"]

one.base<-nchar(as.character(to.convert.ref))==1 &  nchar(as.character(to.convert.alt))==1
sum(one.base)
sum(ref.diff)

  
to.convert.ref[one.base]<-flip.one.base(to.convert.ref[one.base])       # handles "-"          
to.convert.alt[one.base]<-flip.one.base(to.convert.alt[one.base])      # handles "-"      

## cbind(forward.remain[an.indel,],to.convert.ref[an.indel],to.convert.alt[an.indel],all.genomic.ori[ref.diff][an.indel],indels.ori[ref.diff,][an.indel,])

# icon<-1
if(length(to.convert.ref[!one.base])!=0){
for(icon in 1:length(to.convert.ref[!one.base])){
to.convert.ref[!one.base][icon]<-as.character(reverseComplement(DNAString(to.convert.ref[!one.base][icon])))
to.convert.alt[!one.base][icon]<-as.character(reverseComplement(DNAString(to.convert.alt[!one.base][icon])))
}
} ## more than one multibase



indels[ref.diff,"REF"]<-to.convert.ref
indels[ref.diff,"ALT"]<-to.convert.alt


} # somthing to do in step 4



ref.diff<-indels[,"REF"]!=all.genomic
sum(ref.diff)

########### (5) ##########################################################################################
if(sum(ref.diff)!=0){ ## something to do
#cbind(forward.remain,all.genomic.ori[ref.diff],indels[ref.diff,],all.genomic[ref.diff],indels.ori[ref.diff,])


## exm-IND15-41876706 out by one base +1 CT genomic wanted? 
## exm-IND16-87694958
## think both above are just denovo insertions

#### OF ones I have dome before like above they are either polymorphic OR the are native insertions and so am trying to
#### match the new ALT alleles with the reference which will never match. So in these cases make sure any indels are in fact INSERTIONS
### the insertion IS assumed now to be on the forward stand, it be be different which would change the protein....
### swap REF AND ALT to see is just minor allele / major allele swap

is.del<-indels[,"ALT"]=="-"
ref.diff<-ref.diff & is.del
sum(ref.diff)

################ make remaining fails insertions
if(sum(ref.diff)!=0){
to.flip[ref.diff]<-to.flip[ref.diff]+1
### swap REF AND ALT to see is just minor allele / major allele swap
the.ref<-indels[ref.diff,"REF"]
indels[ref.diff,"REF"]<-indels[ref.diff,"ALT"]
indels[ref.diff,"ALT"]<-the.ref
###### swap 0/1 --- keep track of times swapped

############ for indels must chenge start and end and the get new genomic sequence
an.indel<- (nchar(as.character(indels[ref.diff,"REF"])) >1) | (nchar(as.character(indels[ref.diff,"ALT"])) > 1) | indels[ref.diff,"REF"]=="-" | indels[ref.diff,"ALT"]=="-"
sum(an.indel)
if(sum(an.indel)>0){

if(sum(ref.diff)==1){
temp<-  redo.start.end.annovar( subset(subset(indels,subset=ref.diff),subset=an.indel,select=c("chr","start","end","REF","ALT")) )
indels[ref.diff,"chr"][an.indel]<-temp[,"chr"]
indels[ref.diff,"start"][an.indel]<-temp[,"start"]
indels[ref.diff,"end"][an.indel]<-temp[,"end"]
indels[ref.diff,"REF"][an.indel]<-temp[,"REF"]
indels[ref.diff,"ALT"][an.indel]<-temp[,"ALT"]
}else{
  indels[ref.diff,][an.indel,c("chr","start","end","REF","ALT")]<-redo.start.end.annovar(subset(indels[ref.diff,],subset=an.indel,select=c("chr","start","end","REF","ALT")))
}
  


}
############ 5 make reamining fails insertions 


}
} # finish step 5
ref.diff<-indels[,"REF"]!=all.genomic
sum(ref.diff)
to.flip<-to.flip!=0
#sum(to.flip)

rm(ref.ori)
rm(alt.ori)
rm(all.genomic)
rm(all.genomic.have)



the.samples<-colnames(indels)[grepl(".GT$",colnames(indels))]
if(length(the.samples)>0 & sum(to.flip)>0){
    geno.to.flip<-indels[to.flip,the.samples]
    geno.to.flip[geno.to.flip=="0/0"]<-"9/9"
    geno.to.flip[geno.to.flip=="1/1"]<-"0/0"
    geno.to.flip[geno.to.flip=="9/9"]<-"1/1"
    indels[to.flip,the.samples]<-geno.to.flip
  }

return(indels)

} ## function get.forward.start.allele


## out.list<-list(indels,to.flip)
## names(out.list)<-c("indels","to.flip")
## out.list
## ################## Unwind using code like
## start.data<-prepare.for.Vcf.file.read(sample.files[isamp])
## for(i in 1:length(start.data)){assign( names(start.data)[i],value=start.data[[i]])}
## ###########################################




##################################################################################################################################################



# the.indels<-a.indel[1:5,]

 write.plink<-function (the.indels,root) ## note this works with my read.plink function
{
  options("scipen"=300) ## if this does not work use "format" and digits
# c("chr","start",,"ALT","REF","TYPE") columns required : "end" "TYPE" added is needed

  
     ## readBin(con, what, n = 1L, size = NA_integer_, signed = TRUE,
     ##         endian = .Platform$endian)
     
     ## writeBin(object, con, size = NA_integer_,
     ##          endian = .Platform$endian, useBytes = FALSE)
## need a 2*(num of samples vector labelled (01 00 for missing) (01 01 for 1/1)  (00 01 for 0/1 ) (00 00  for 1/1) 
##  00 00 00 00 00 00#  at the end

          ## need to write 1 snp at a time
          ## need to run genotype. summary to get minor allele for write of bim ## is flip indels need to chnage the start/end
          ## wite the fam file
          ## write wite checck for same starts?
   
    #the.indels<-as.matrix(the.indels) ###### not SURE IF NEED THIS - yes makes replacement a lot faster
   core.ann<-c("chr","start","ALT","REF")
  if(sum(core.ann %in% colnames(the.indels))!=length(core.ann)){print(paste("Missing columns",core.ann[!(core.ann %in% colnames(the.indels))],"in input matrix",sep=" "))}

  
  if(is.null(dim(the.indels))){
    the.indels<-t(as.matrix(the.indels))
  }else{
    the.indels<-as.matrix(the.indels) ## Do this case way faster to to substitustion in a matrix
  }

  if(!("end" %in% colnames(the.indels))){
    end<-the.indels[,"start"]
    the.indels<-cbind(the.indels,end)
 
if(dim(the.indels)[i]==1){
temp<-  redo.start.end.annovar( subset(the.indels,select=c("chr","start","end","REF","ALT")) )
the.indels[,"chr"]  <-temp[,"chr"]
the.indels[,"start"]<-temp[,"start"]
the.indels[,"end"]  <-temp[,"end"]
the.indels[,"REF"]  <-temp[,"REF"]
the.indels[,"ALT"]  <-temp[,"ALT"]
}else{
  the.indels[,c("chr","start","end","REF","ALT")]<-redo.start.end.annovar( subset(the.indels,select=c("chr","start","end","REF","ALT")) )
}
 } # no end label


  if(!("TYPE" %in% colnames(the.indels))){
    TYPE="GENO"
    the.indels<-cbind(the.indels,TYPE)
  }

     the.indels[,"start"]<-gsub("^\\s+","",the.indels[,"start"])
     the.indels[,"end"]<-gsub("^\\s+","",the.indels[,"end"])

########### plink expect specifi encoding from chromosomes
    if(grepl("^chr",the.indels[1,"chr"])){the.indels[,"chr"]<-gsub("chr","",the.indels[,"chr"])}


    the.indels[the.indels[,"chr"]=="X","chr"]<-"23"  ## no Y observed
    the.indels[the.indels[,"chr"]=="Y","chr"]<-"24"
    the.indels[the.indels[,"chr"]=="XY","chr"]<-"25" ## chrXY is chr25 but coordinates/positions ARE on chrX (note no XY in UCSC build)
    the.indels[the.indels[,"chr"]=="M","chr"]<-"26"
 #   the.indels[the.indels[,"chr"]=="chrXY","chr"]<-"chrX"

    

      
    bed.file = paste(root, ".bed", sep = "")
    the.samples<- colnames(the.indels)[grepl(".GT$",colnames(the.indels))] #  gsub(".GT$","",colnames(indels)[grepl(".GT$",colnames(indels))])
    n.snps = dim(the.indels)[1]

    ## test.bytes = readBin(bin.connection, what = "raw", n = 3)
    ## if (!identical(as.character(test.bytes), c("6c", "1b", "01"))) {
    ##     stop("BED file not a v0.99 SNP-major BED file, please re-encode the data as v0.99 SNP-major file")
    ## }
    

    summary.geno<-genotype.summary(the.indels[,the.samples])
    colnames(summary.geno)<-c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO")



    ########## ALT most often th monor allele and set to 1. Plink expects MAJOR to be 1
    ######### 1) Ensure ALT is the mnor allele
    ######### 2) the recode the change

    
    ### plink A1 is the minor allele MAF is the minor allele  frequency
    ### Genotyping is w.r.t ALT. When is ALT NOT the minor allele
    ref.diff<-as.numeric(summary.geno[,"MAF"]) > 0.5 & !is.na(as.numeric(summary.geno[,"MAF"])) # say all those where ALT is minor.
    ###### plink BIM expects A1  (here use REF) as the A1 position
    ### swap REF AND ALT to see is just minor allele / major allele swap
    if(sum(ref.diff>0)){
    the.ref<-the.indels[ref.diff,"REF"]
    the.indels[ref.diff,"REF"]<-the.indels[ref.diff,"ALT"]
    the.indels[ref.diff,"ALT"]<-the.ref
    ###### swap 0/1 --- keep track of times swapped

    ############ for indels must chenge start and end and the get new genomic sequence
    an.indel<- (nchar(as.character(the.indels[ref.diff,"REF"])) >1) | (nchar(as.character(the.indels[ref.diff,"ALT"])) > 1) | the.indels[ref.diff,"REF"]=="-" | the.indels[ref.diff,"ALT"]=="-"
    sum(an.indel)

   if(sum(an.indel)>0){
#   the.indels[ref.diff,][an.indel,c("chr","start","end","REF","ALT")]<-redo.start.end.annovar(the.indels[ref.diff,][an.indel,c("chr","start","end","REF","ALT")])

if(sum(ref.diff)==1){
temp<-  redo.start.end.annovar( subset(subset(the.indels,subset=ref.diff),subset=an.indel,select=c("chr","start","end","REF","ALT")) )
the.indels[ref.diff,"chr"][an.indel]<-temp[,"chr"]
the.indels[ref.diff,"start"][an.indel]<-temp[,"start"]
the.indels[ref.diff,"end"][an.indel]<-temp[,"end"]
the.indels[ref.diff,"REF"][an.indel]<-temp[,"REF"]
the.indels[ref.diff,"ALT"][an.indel]<-temp[,"ALT"]
}else{
  the.indels[ref.diff,][an.indel,c("chr","start","end","REF","ALT")]<-redo.start.end.annovar(subset(the.indels[ref.diff,],subset=an.indel,select=c("chr","start","end","REF","ALT")))
}



   
  }
  } # sum(ref.diff)==0
   ######################################################

    ############## also flip the genotypes
    geno.to.flip<-the.indels[ref.diff,the.samples]
    geno.to.flip[geno.to.flip=="0/0"]<-"9/9"
    geno.to.flip[geno.to.flip=="1/1"]<-"0/0"
    geno.to.flip[geno.to.flip=="9/9"]<-"1/1"
    the.indels[ref.diff,the.samples]<-geno.to.flip

########## IF ALT was the minor allele I flipped it so REF is now the MINOR and 


    
  ## ALT is always the minor allele encode so major (REF) is 1 not 0 - like it is now.
    
    genotypes<-the.indels[,the.samples]
    genotypes<-as.matrix(genotypes)
    
    genotypes[genotypes=="NA" | is.na(genotypes)]<-"01 00"
    genotypes[genotypes=="0/0"]<-"01 01" ## encode 2
    genotypes[genotypes=="0/1"]<-"00 01" ## encode 1
    genotypes[genotypes=="1/1"]<-"00 00" ## encode 0

 


   the.key<-build.key(the.indels,c("chr","start","end","ALT","REF","TYPE"))

      if(!("SNP" %in% colnames(the.indels))){
        SNP<-the.key
      }else{
        SNP<-the.indels[,"SNP"]
      }

    
   write.table(cbind(the.indels[,c("chr")],SNP,0,the.indels[,c("start","ALT","REF")]),file=paste(root, ".bim", sep = ""),col.name=FALSE,row.names=FALSE,quote = FALSE) ## ALT is A1
    
   sample.labels<-gsub(".GT$","",the.samples)
   write.table(cbind(sample.labels,sample.labels,0,0,-9,-9),file=paste(root, ".fam", sep = ""),col.name=FALSE,row.names=FALSE,quote = FALSE) ## ALT is A1

    

    ## r.bin.snp.ori<-r.bin.snp
    ## test<-rawToBits(r.bin.snp)
    ## packBits(test,type="raw")
    ##  gen<-bin.snp[, 1] + bin.snp[, 2] - 10 * ((bin.snp[, 
    ##         1] == 1) & (bin.snp[, 2] == 0))
    ## packBits(intToBits(gen),type="raw")
    system(paste("rm -f ",bed.file))
    bin.connection = file(bed.file, "wb")
    plink.start<-packBits(as.raw(c("00","00","01","01","00","01","01","00","01","01","00","01","01","00","00","00","01","00","00","00","00","00","00","00")))
    writeBin(plink.start,bin.connection)
    
    geno.tail<-as.raw(c("00","00","00","00","00","00"))
    for (i in 1:n.snps) {
      geno.str<-paste(genotypes[i,],collapse=" ")
      geno.str<- unlist(strsplit(geno.str,split=" "))
      geno.str<-as.raw(geno.str)
      geno.str<-c(geno.str,geno.tail)
      geno.str<-packBits(geno.str)
      writeBin(geno.str,bin.connection)
    }

    close(bin.connection)
#    return(genotypes)
}



    ## genotypes[genotypes==0]<-NA
    ## genotypes[genotypes==0]<-"0/0"
    ## genotypes[genotypes==1]<-"0/1"
    ## genotypes[genotypes==2]<-"1/1"


## root<-paste(genotype.file.location,plink.file,sep="/")
## root<- "/media/UQCCG-Analysis/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL_chr21"

 read.plink<-function(root){ ## reda in a bed file add bim information do not mess with A1 or A2 just leave as is
   ##### WARNING  HERE A2 is labelled "1" and A1 is "0" A1 is the minor allele ## So NOT what I would expect most time arer allele ALT or A1 is usually "1"
   #### So Dosage is w.r.t A2 
    bed.file = paste(root, ".bed", sep = "")
    bed.file.size = file.info(bed.file)$size
    
    snp.names = read.table(paste(root, ".bim", sep = ""))
    colnames(snp.names)<-c("chr","SNP","cm","start","A1","A2")
    fam.names = read.table(paste(root, ".fam", sep = ""))
    colnames(fam.names)<-c("FID","IID","MOTHER","FATHER","SEX","AFF")
    
 
    #sample.size = dim(read.table(paste(root, ".fam", sep = "")))[1]
    sample.size = dim(fam.names)[1]
    snp.size = ceiling(sample.size/4)
    n.snps = round((bed.file.size - 3)/snp.size)
    bin.connection = file(bed.file, "rb")
    test.bytes = readBin(bin.connection, what = "raw", n = 3)
    if (!identical(as.character(test.bytes), c("6c", "1b", "01"))) {
        stop("BED file not a v0.99 SNP-major BED file, please re-encode the data as v0.99 SNP-major file")
    }
    genotypes = matrix(ncol = n.snps, nrow = sample.size)
    for (i in 1:n.snps) {
        r.bin.snp = readBin(bin.connection, what = "raw", n = snp.size)
        bin.snp = matrix(as.numeric(rawToBits(r.bin.snp)), ncol = 2,  byrow = TRUE)[1:sample.size, ]
        genotypes[, i] = bin.snp[, 1] + bin.snp[, 2] - 10 * ((bin.snp[,  1] == 1) & (bin.snp[, 2] == 0))
    }
    genotypes[genotypes == -9] = NA
     close(bin.connection)
    genotypes<-t(genotypes)
    ############# HERE A2 is labelled "1" and A1 is "0" A1 is the minor allele ## So NOT what I would expect most time ALT would be A1
    ############ the dosage is w.r.t A2 the most common allele. Change this labelling now so A1 is 1

    genotypes<-as.matrix(genotypes) ## Do this case way faster to to substitustion in a mtraix
    ## a.indel[,"start"]<-gsub("^\\s+","",a.indel[,"start"])
    

    genotypes[genotypes==0]<-"1/1"
    genotypes[genotypes==1]<-"0/1"
    genotypes[genotypes==2]<-"0/0"

     colnames(snp.names)[ colnames(snp.names)=="A1"]<-"ALT"
     colnames(snp.names)[ colnames(snp.names)=="A2"]<-"REF"
   SNP<-as.character(snp.names[,"SNP"])
    ## summary.geno<-genotype.summary(genotypes)
    ## colnames(summary.geno)<-c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO")
    
    colnames(genotypes) = paste(fam.names[,"IID"],".GT",sep="")
    the.samples<-colnames(genotypes)
    end<-snp.names[,"start"]
    TYPE<-rep("GENO",times=length(end))
    indels<-cbind(snp.names[,c("chr","start")],end,snp.names[,c("REF","ALT")],SNP,TYPE ,genotypes,stringsAsFactors = FALSE)

    indels<-as.matrix(indels) ## Do this case way faster to to substitustion in a matrix
    indels[,"start"]<-gsub("^\\s+","",indels[,"start"])
    
#
    indels[,c("chr","start","end","REF","ALT")]<-redo.start.end.annovar(indels[,c("chr","start","end","REF","ALT")]) ## put in annovar format


    ###### check have the forward startd allele
  
    print(dim(indels))
 #   setwd(code.dir)
    indels<-get.forward.strand.allele(indels) ## function conatins flip if samples are present (samples have a ".GT" extension in column name
   # source("get.forward.strand.allele.r") ## this should include genotype inversion

   ### flip genotypes of to.flip  to.flip HAve been flipped BUT genotypes have NOT BEEN

    ## geno.to.flip<-indels[to.flip,the.samples]
    ## geno.to.flip[geno.to.flip=="0/0"]<-"9/9"
    ## geno.to.flip[geno.to.flip=="1/1"]<-"0/0"
    ## geno.to.flip[geno.to.flip=="9/9"]<-"1/1"
    ## indels[to.flip,the.samples]<-geno.to.flip
    ################
    ## indels[to.flip,the.samples][1:5,1:5]
    ## rownames(genotypes) = snp.names
  
    return(indels)
}






 read.plink.raw<-function(root){ ## reda in a bed file add bim information do not mess with A1 or A2 just leave as is
   ##### WARNING  HERE A2 is labelled "1" and A1 is "0" A1 is the minor allele ## So NOT what I would expect most time arer allele ALT or A1 is usually "1"
   #### So Dosage is w.r.t A2 
    bed.file = paste(root, ".bed", sep = "")
    bed.file.size = file.info(bed.file)$size
    
    snp.names = read.table(paste(root, ".bim", sep = ""))
    colnames(snp.names)<-c("chr","SNP","cm","POS","A1","A2")
#    snp.names<-as.matrix(snp.names)
    fam.names = read.table(paste(root, ".fam", sep = ""))
    colnames(fam.names)<-c("FID","IID","MOTHER","FATHER","SEX","AFF")
    
 
    #sample.size = dim(read.table(paste(root, ".fam", sep = "")))[1]
    sample.size = dim(fam.names)[1]
    snp.size = ceiling(sample.size/4)
    n.snps = round((bed.file.size - 3)/snp.size)
    bin.connection = file(bed.file, "rb")
    test.bytes = readBin(bin.connection, what = "raw", n = 3)
    if (!identical(as.character(test.bytes), c("6c", "1b", "01"))) {
        stop("BED file not a v0.99 SNP-major BED file, please re-encode the data as v0.99 SNP-major file")
    }
    genotypes = matrix(ncol = n.snps, nrow = sample.size)


    
    
    for (i in 1:n.snps) {
 #       for (i in 1:16) {
        r.bin.snp = readBin(bin.connection, what = "raw", n = snp.size)
        bin.snp = matrix(as.numeric(rawToBits(r.bin.snp)), ncol = 2,  byrow = TRUE)[1:sample.size, ]
        genotypes[, i] = bin.snp[, 1] + bin.snp[, 2] - 10 * ((bin.snp[,  1] == 1) & (bin.snp[, 2] == 0))
    }
    genotypes[genotypes == -9] = NA

    genotypes<-t(genotypes)
    ############# HERE A2 is labelled "1" and A1 is "0" A1 is the minor allele ## So NOT what I would expect most time ALT would be A1

        ## a.indel<-as.matrix(a.indel) ## Do this case way faster to to substitustion in a matrix
    ## a.indel[,"start"]<-gsub("^\\s+","",a.indel[,"start"])
    
    genotypes<-as.matrix(genotypes) 
    genotypes[genotypes==0]<-"1/1"
    genotypes[genotypes==1]<-"0/1"
    genotypes[genotypes==2]<-"0/0"
    
    colnames(snp.names)[ colnames(snp.names)=="A1"]<-"ALT"
    colnames(snp.names)[ colnames(snp.names)=="A2"]<-"REF"
     colnames(snp.names)[ colnames(snp.names)=="POS"]<-"start"
    SNP<-as.character(snp.names[,"SNP"])
    
    colnames(genotypes) = paste(fam.names[,"IID"],".GT",sep="")
    genotypes<-cbind(snp.names[,!(colnames(snp.names) %in% "SNP")],SNP,genotypes,stringsAsFactors = FALSE)

    genotypes<-as.matrix(genotypes) ## Do this case way faster to to substitustion in a matrix
    genotypes[,"start"]<-gsub("^\\s+","",genotypes[,"start"])

    
    ## rownames(genotypes) = snp.names
    close(bin.connection)
    return(genotypes)
}



clean.unique.combine<-function(x,y){
  z<-unique(c(x,y))
  z[!is.na(z) & !grepl("^PLACE_",z)] #z[z!="-" & !is.na(z) & !grepl("^PLACE_",z)]
}


get.alleles.per.region<-function(indels,regions.gr,allele.counts,filter){
  ## regions.gr<-the.exons
  ## allele.counts<-genos.high[,"HET"]
  ## filter<-pass
  ## used in Analyse_group_project
  ## assumes indels has c("chr","start","end") and strand is +
  ## results counts of length regions
  core.ann<-c("chr","start","end")
  if(sum( (core.ann %in% colnames(indels)))!=length(core.ann)){print("ERROR missing core annoaation")}
  if(length(allele.counts)!=dim(indels)[1]){print("Error allele counts missing")}
  if(length(filter)!=dim(indels)[1]){filter<-rep(TRUE,times=dim(indels)[1]);print("WARNINQ no filtering")}
  
  indels.pass<-indels[filter,c("chr","start","end")]
  allele.counts<-as.numeric(allele.counts[filter])
 # allele.counts[is.na(allele.counts)]
  
  flat.index<-rep(1:dim(indels.pass)[1],times=allele.counts) ##works fine when allele counts has zeros
#  indels.pass.ori.size<-dim(indels.pass)[1]
  indels.pass<-indels.pass[flat.index,]
  indels.gr<-GRanges(seqnames =indels.pass[,"chr"],ranges = IRanges(start=as.numeric(indels.pass[,"start"]),end=as.numeric(indels.pass[,"end"])),strand="+")

  the.chromo<-as.character(unique(seqnames(indels.gr)))
#  human.chromlens<-the.chroms[the.chromo]
   regions.gr<-data.gr[seqnames(data.gr) %in% the.chromo,]
   countOverlaps(regions.gr,indels.gr,type="any")
}

get.genotype.counts<-function(genos){
  ## used in Analyse_greoup_project
  ## splits up indels[1,"GENO.HIGH"] == "1,4,496"
genos<-strsplit(genos,split=",")
genos.length<-length(genos)
genos<-as.numeric(unlist((genos)))

dim(genos)<-c(3,genos.length)
genos<-t(genos)
colnames(genos)<-c("ALT","HET","REF")
genos
}



combine.boolean<-function(a.matrix,columns.to.combine,combine.type,NA.set.to=TRUE){
  ### Combine columns in a matrix using OR or AND per specified column OR USE cbind(vector1,vector2) to make a.matrix  with columns.to.combine of length one
  ### Returns a bollean array
  ############ set any NA's to true
  if(length(columns.to.combine)==1){colnames(a.matrix)<-1:dim(a.matrix)[2];columns.to.combine<-colnames(a.matrix)}
  print(columns.to.combine)
  the.NAs<-is.na(a.matrix[,columns.to.combine])
  count.col.NAs<-apply(as.matrix(the.NAs),2,sum)# one filter .col gives a vector
  found.col.NAs<-count.col.NAs>0
  if(sum(found.col.NAs)>0){
    for(i in 1:length(found.col.NAs)){
      if(found.col.NAs[i]){
        a.matrix[,columns.to.combine[i]][is.na(a.matrix[,columns.to.combine[i]])]<-NA.set.to
      }}}
  ##########################
  
  
  num.trues<-apply(as.matrix(a.matrix[,columns.to.combine]),1,function(x) sum(x,na.rm=TRUE))
  if(combine.type=="OR"){the.filter<-num.trues >0} # ANY test is true the filter out 
  if(combine.type=="AND"){the.filter<-num.trues==length(columns.to.combine)} # ALL test must be true the filter out 
  the.filter
}


## data<-pholy.data
## score.col<-"PolyPhen.scores"
## desc.col<-"PolyPhen.desc"

functional.score.to.desc.calibrate<-function(data,score.col,desc.col){
  ## print(paste(score.col,desc.col))
  ## print(colnames(data))
  if(sum(!c(score.col,desc.col) %in% colnames(data)) >0){print(paste("ERROR columns not in table:",toString(colnames(data))));NA}else{
min.value<-tapply(as.numeric(data[,score.col]),data[,desc.col],function(x) min(x,na.rm=TRUE))
max.value<-tapply(as.numeric(data[,score.col]),data[,desc.col],function(x) max(x,na.rm=TRUE))
max.value<-max.value[names(min.value)]
cal.data<-cbind(min.value,max.value)
the.order<-order(cal.data[,"max.value"],decreasing=TRUE)
cal.data[the.order,]
}
}


apply.calibration<-function(scores,calibration){
  desc<-rep(NA,times=length(scores))
  for(i in 1:dim(calibration)[1]){
    found<-(as.numeric(scores)>=calibration[i,"min.value"]) & (as.numeric(scores)<=calibration[i,"max.value"])
    desc[found]<-rownames(calibration)[i]
  }
  desc
}

## pholy.format[a.deletion,][1:15,]
## pholy.format[a.insert,][1:15,]
## pholy.format[c(26,37),]
## the.indels<-indels
## the.indels<-subset(indels,missing)
##  check<-indels[,"REF"]=="-" | indels[,"ALT"]=="-"
## indels[check,][1:30,1:10]

#ok as of 9th may 2013  ### ADD TO CBIND ,stringsAsFactors = FALSE
to.pholy.format<-function(the.indels){
# see  http://www.ensembl.org/info/docs/variation/vep/vep_formats.html
  pholy.format<-cbind(subset(the.indels,select=c("chr","start","end")),paste(the.indels[,"REF"],"/",the.indels[,"ALT"],sep=""),"+")
  ## a.deletion<-the.indels[,"ALT"]=="-" deletions handled ok
  a.insert<-the.indels[,"REF"]=="-"
  pholy.format[a.insert,"start"]<-as.integer(as.integer(pholy.format[a.insert,"end"])+1) ## pholyphen has a strange format!
  pholy.format[,"start"]<-as.numeric(pholy.format[,"start"])
  pholy.format[,"end"]<-as.numeric(pholy.format[,"end"])
  pholy.format
}


process.format<-function(x) {
    posns<-match(format.string,unlist(strsplit(x[1],split=":")))
    unlist(strsplit(x[2],split=":"))[posns]
  }

process.format.list<-function(x) {
  x<-unlist(x)
    posns<-match(format.string,unlist(strsplit(x[1],split=":")))
    unlist(strsplit(x[2],split=":"))[posns]
  }

process.info.list<-function(x,info.types) {
         one<-strsplit(x,split="=")
         one<-sapply(one, function(x) c(x[1],x[2]))
         reorder<-match(info.types,one[1,])
         one[2,reorder]
       }

process.info<-function(x,info.types) {
      posns<-match(info.types,unlist(strsplit(x[1],split=":")))
       one<-unlist(strsplit(x,split=";"))
         one<-strsplit(one,split="=")
         one<-sapply(one, function(x) c(x[1],x[2]))
         reorder<-match(info.types,one[1,])
         one[2,reorder]
       }




## names(types)[grep("frameshift_variant",names(types))]/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/PCA_calc/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL_withRS.bed
## types<-sort(tapply(poly[,"Consequence"],poly[,"Consequence"],length))
## to.get<-"Consequence"
## order.by<-vep.types
##  to.get.other<-c("Uploaded_variation","Gene","Feature","Protein_position","Amino_acids")
## to.get.other<-{}
##indels<-subset(indels[ref.diff,,drop=F],subset=an.indel,select=c("chr","start","end","REF","ALT"))
#test<-indels[ref.diff,][an.indel,c("chr","start","end","REF","ALT")]/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.bed
## /media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL_withRS.bim
check.chr.and.positions<-function(indels,chr.lengths){
 ### used mainli for bim files to see if location has carp position, chromosome or length prior to getSeq
#  if(genome.build=="hg19"){library("BSgenome.Hsapiens.UCSC.hg19"); have.chr<-seqnames(Hsapiens);chr.lengths<-seqlengths(Hsapiens)}
  if( ( dim(indels)[1]>0 | !is.null(dim(indels))) ){
if(!grepl("^chr",indels[1,"chr"])){indels[,"chr"]<-paste("chr",indels[,"chr"],sep="")}
indels[indels[,"chr"]=="chr23","chr"]<-"chrX"  ## no Y observed
indels[indels[,"chr"]=="chr24","chr"]<-"chrY"
indels[indels[,"chr"]=="chr25","chr"]<-"chrXY"
indels[indels[,"chr"]=="chr26","chr"]<-"chrMT"
no.chr<-!(as.character(indels[,"chr"]) %in% names(chr.lengths))
if(sum(no.chr)>0){
  indels[indels[,"chr"]=="chrMT","chr"]<-"chrM"
}

  
positions<-chr.lengths[indels[,"chr"]]
positions[is.na(positions)]<-0
too.large<-as.integer(indels[,"start"]) > positions
no.start<-as.integer(indels[,"start"])==0

bad.locations<-too.large | no.chr | is.na(positions) | no.start
}else{
bad.locations<-{} }                                       # send in and empty array

bad.locations
}

flip.one.base<-function(to.flip){
  ### flip one base only to.flip is just an vector of A,T,G,C's
  As<-to.flip=="A"
  Ts<-to.flip=="T"
  Gs<-to.flip=="G"
  Cs<-to.flip=="C"
 # Dels<-to.flip=="-"
  
   to.flip[As]<-"T"
    to.flip[Ts]<-"A"
    to.flip[Gs]<-"C"
    to.flip[Cs]<-"G"
 #  to.flip[Dels]<-"-"
  to.flip
}


redo.start.end.annovar<-function(a.indel){ ## this funstion used in get.forward startnd alleles
  ## casue is you SWAP alleles you must change the start and end location for indels
########################### PUT ALLELES IN ANNOVAR FORMAT from convert2annovar.pl line 1083########
ref.length<-nchar(as.character(a.indel[,"REF"]))
alt.length<-nchar(as.character(a.indel[,"ALT"]))
is.snp<-(ref.length==1 & alt.length==1) ## indels that are "-" are included here start=end for those for both insertions and deletions

# sum(!is.snp)

POS.end<-as.numeric(a.indel[,"start"])
del<-ref.length > alt.length
ins<-(ref.length <= alt.length) & !is.snp
## a.indel[del,][1:5,]
## POS.end[del][1:5]
### deletion or block substitution
head<-substr(as.character(a.indel[del,"REF"]),1,alt.length[del])
head.is.mut<-(head==as.character(a.indel[del,"ALT"]))
a.indel[del,"REF"][head.is.mut]<-substr(as.character(a.indel[del,"REF"][head.is.mut]),(alt.length[del][head.is.mut]+1),ref.length[del][head.is.mut])
a.indel[del,"ALT"][head.is.mut]<-"-"
a.indel[del,"start"][head.is.mut]<-as.numeric(a.indel[del,"start"][head.is.mut]) + nchar(as.character(head[head.is.mut]))
#POS.end[del]<-POS.end[del]+ref.length[del]-1  # $$$ use this line in initial VCF read code when end does not exist 
a.indel[del,"end"]<-POS.end[del]+ref.length[del]-1 
## a.indel
## POS.end
### insertion or block substitution
head<-substr(as.character(a.indel[ins,"ALT"]),1,ref.length[ins])
head.is.ref<-(head==as.character(a.indel[ins,"REF"]))
a.indel[ins,"ALT"][head.is.ref]<-substr(as.character(a.indel[ins,"ALT"][head.is.ref]),(ref.length[ins][head.is.ref]+1),alt.length[ins][head.is.ref])
a.indel[ins,"REF"][head.is.ref]<-"-"
a.indel[ins,"start"][head.is.ref]<-as.numeric(a.indel[ins,"start"][head.is.ref]) + ref.length[ins][head.is.ref]-1
#POS.end[ins]<-POS.end[ins]+ref.length[ins]-1 # $$$ use this line in initial VCF read code when end does not exist 
a.indel[ins,"end"]<-POS.end[ins]+ref.length[ins]-1
a.indel
}



## table<-poly # key.pholy<-paste(poly[,"Location"],poly[,"Allele"],sep=":")
build.pholy.key<-function(table){  ### used to matcch indels across table use build.key(indels,c("chr","start","ALT")) for indels 
  if(is.null(dim(table))){table<-as.matrix(table)} # in casea vector sent
  if(sum(!(c("Location","Allele") %in% colnames(table)))>0){print("FAIL location and Allele not present");key.pholy<-1:dim(table)[1]}else{
    location<-strsplit(table[,"Location"],split="-")
    location<-unlist(lapply(location,function(x) x[1])) # mapply(function(x) x[1],location[1:20],SIMPLIFY=TRUE)    
    key.pholy<-paste(location,table[,"Allele"],sep=":")
}
key.pholy
}

#get.Consequence.from.VEP(indels,vep.types,poly,"Consequence",extra.vep.annotations)
get.Consequence.from.VEP<-function(indels,order.by,poly,to.get,to.get.other={}){
  ## pholy data is repeated many time so need to collapse
  ## this function gets the Polyphen and sift data which is encoded in a character string in the "Extra" column
  ## vep.types is rank ordered 
## poly.scores<-poly[,to.get]
 
  #### unwind the comma delimited Consequences


## get  "intron_variant,nc_transcript_variant"    splice_acceptor_variant,nc_transcript_variant  to getnode in from non-codiing gene
poly.scores<-poly[,to.get]
non.coding<-grep("nc_transcript_variant",poly.scores)
poly.scores[non.coding]<-paste("NC_",poly.scores[non.coding],sep="")      # poly.scores[non.coding][1:20]

poly.scores<-strsplit(poly.scores,split=",")
reps<-unlist(lapply(poly.scores,length))
key.pholy<-build.pholy.key(poly) # key.pholy<-paste(poly[,"Location"],poly[,"Allele"],sep=":")

index.use<-rep(c(1:length(reps)),times=reps) # used rep(key.pholy,times=reps) before by need index to get other cols
key.pholy<-key.pholy[index.use]  # key.pholy<-rep(key.pholy,times=reps)#  
poly.scores<-unlist(poly.scores)
if(length(to.get.other)>1){
poly.scores<-cbind(key.pholy ,poly.scores,poly[index.use,c(to.get,to.get.other)]) # ADD to.get here for keep original format
}else{
poly.scores<-cbind(key.pholy ,poly.scores)
}

############################################function(table,key.cols,add.chr.label=FALSE,delim=":"
key.indels.poly<-build.key(indels,c("chr","start","ALT"))
## key.indels.poly<-paste(build.key.delimit(indels,c("chr","start"),FALSE,"_"),paste(indels[,"REF"],"/",indels[,"ALT"],sep=""),sep="_")  # key.indels.poly<-build.key(indels,c("chr","start","ALT"))
## pholy.format<-paste(the.indels[,c("chr","start","end")],paste(the.indels[,"REF"],"/",the.indels[,"ALT"],sep=""),sep="_")
##poly scores comtains repeats for diff trnascripts choose the most damaging
if(dim(poly.scores)[1]>1){
  
possible<-unique(poly.scores[,"poly.scores"])

## the.types<-names(tapply(poly[,"Consequence"],poly[,"Consequence"],length))
#the.types<-unique(unlist(strsplit(the.types,split=";")))
if( sum( !(possible %in% vep.types))>0  ){print("WARNING VEP HAS NEW MUTATION TYPES DEFINED- REVIEW get.Cosequence.from.VEP in sunbroutines")
                                           print( possible [!(possible  %in% vep.types)])    }

new<-possible[!(possible %in% order.by)]
for(iu in 1:length(new)){
  poly.scores[poly.scores[,"poly.scores"]==new[iu],"poly.scores"]<-"not_assigned"
}
posns<-match(order.by,possible)
missing<-is.na(posns)
order.by<-order.by[!missing]
the.order.use<-match(poly.scores[,"poly.scores"],order.by) ## match position of verp muation that are ordered worst->ok
the.order<-order(the.order.use)
poly.scores<-poly.scores[the.order,]  ### with most damaging at beginning now 
dups<-duplicated(poly.scores[,"key.pholy"])
poly.scores<-poly.scores[!dups,]
}


## key.indels.poly[1:15]  
    ## list.element.lengths<-unlist(lapply(the.genes[[annotation.labels[i]]],length))        # list.lentghs
    ## indel.index<-rep(1:length(list.element.lengths),times=list.element.lengths)           # index
    ## flat.gene.list<-unlist(the.genes[[annotation.labels[i]]])   
if(is.null(dim(poly.scores))){poly.scores<-t(as.matrix(poly.scores))}
poly.scores<-match.cols.and.collapse(list(key.indels.poly),1,poly.scores,"key.pholy",colnames(poly.scores),"delimit")
colnames(poly.scores)<-gsub("poly.scores",paste(to.get,"SHORT",sep="_"),colnames(poly.scores))
poly.scores
}


## grep("156479569$",indels[,"start"])
## indels[612170:612179,1:10]

## key.indels.poly[indels[,"TYPE"]=="indel"][1:10]
## poly.scores[indels[,"TYPE"]=="indel",][1:10,]
## indels[indels[,"TYPE"]=="indel",1:7][1:10,]
## sort(tapply(poly.scores[,"poly.scores"],poly.scores[,"poly.scores"],length))

get.Column.from.VEP<-function(indels,order.by,poly,to.get){
  ## pholy data is repeated many time so need to collapse
  ## this function gets the Polyphen and sift data which is encoded in a character string in the "Extra" column
  ## vep.types is rank ordered 
poly.scores<-poly[,to.get]

key.pholy<-build.pholy.key(poly) # paste(poly[,"Location"],poly[,"Allele"],sep=":")
poly.scores<-cbind(key.pholy,poly.scores)
key.indels.poly<-build.key(indels,c("chr","start","ALT"))

##poly scores comtains repeats for diff trnascripts choose the most damaging
if(dim(poly.scores)[1]>1){
  
possible<-unique(poly.scores[,"poly.scores"])
new<-possible[!(possible %in% order.by)]
for(iu in 1:length(new)){
  poly.scores[poly.scores[,"poly.scores"]==new[iu],"poly.scores"]<-"not_assigned"
}
posns<-match(order.by,possible)
missing<-is.na(posns)
order.by<-order.by[!missing]
the.order.use<-match(poly.scores[,"poly.scores"],order.by)
the.order<-order(the.order.use)

poly.scores<-poly.scores[the.order,]  ### with most damaging at beginning now 
dups<-duplicated(poly.scores[,"key.pholy"])
poly.scores<-poly.scores[!dups,]
}
## key.indels.poly[1:15]
    ## list.element.lengths<-unlist(lapply(the.genes[[annotation.labels[i]]],length))        # list.lentghs
    ## indel.index<-rep(1:length(list.element.lengths),times=list.element.lengths)           # index
    ## flat.gene.list<-unlist(the.genes[[annotation.labels[i]]])   
if(is.null(dim(poly.scores))){poly.scores<-t(as.matrix(poly.scores))}
poly.scores<-match.cols.and.collapse(list(key.indels.poly),1,poly.scores,"key.pholy",colnames(poly.scores),"delimit")
colnames(poly.scores)<-gsub("poly.scores",gsub("=","",to.get),colnames(poly.scores))
poly.scores
}







## pholy.data<-get.function.from.VEP(indels,poly,"PolyPhen=")
## pholy.data<-get.function.from.VEP(indels,poly,"PolyPhen=")
## sift.data<-get.function.from.VEP(indels,poly,"SIFT=")  # to.get<-"PolyPhen="

## to.get<-"PolyPhen="
get.function.from.VEP<-function(indels,poly,to.get){
  ## this function gets the Polyphen and sift data which is encoded in a character string in the "Extra" column
has.pholy<-grepl(to.get,poly[,"Extra"]) # sum(has.pholy)
#poly.scores<-extract.value.from.format(poly[has.pholy,"Extra"],"PolyPhen=") #$ this just gives the predictions
poly.scores<-extract.value.from.poly(poly[has.pholy,"Extra"],to.get)
poly.scores<-gsub(")","",poly.scores) 
poly.scores<-strsplit(poly.scores,split="\\(")
#test<-unlist(lapply(poly.scores,length))
poly.desc<-unlist(lapply(poly.scores,function(x) x[1]))
poly.scores<-unlist(lapply(poly.scores,function(x) x[2]))
key.pholy<-build.pholy.key(poly[has.pholy,]) # paste(poly[has.pholy,"Location"],poly[has.pholy,"Allele"],sep=":")
poly.scores<-cbind(key.pholy,poly.desc,poly.scores)


##poly scores comtains repeats for diff trnascripts choose the most damaging
if(dim(poly.scores)[1]>1){
the.order<-order(as.numeric(poly.scores[,"poly.scores"]),decreasing=TRUE)
poly.scores<-poly.scores[the.order,]   
dups<-duplicated(poly.scores[,"key.pholy"])
poly.scores<-poly.scores[!dups,]
}

key.indels.poly<-build.key(indels,c("chr","start","ALT"))
## key.indels.poly[1:15]
    ## list.element.lengths<-unlist(lapply(the.genes[[annotation.labels[i]]],length))        # list.lentghs
    ## indel.index<-rep(1:length(list.element.lengths),times=list.element.lengths)           # index
    ## flat.gene.list<-unlist(the.genes[[annotation.labels[i]]])   
if(is.null(dim(poly.scores))){poly.scores<-t(as.matrix(poly.scores))}
poly.scores<-match.cols.and.collapse(list(key.indels.poly),1,poly.scores,"key.pholy",colnames(poly.scores),"delimit")
colnames(poly.scores)<-gsub("poly",gsub("=","",to.get),colnames(poly.scores))
poly.scores
}


#regulation.data<-get.regulation.from.VEP(indels,poly,c("MotifFeature","RegulatoryFeature")) # to.get<-c("MotifFeature","RegulatoryFeature")
get.regulation.from.VEP<-function(indels,poly,to.get){
has.pholy<-poly[,"Feature_type"] %in% to.get
#poly.scores<-extract.value.from.format(poly[has.pholy,"Extra"],"PolyPhen=") #$ this just gives the predictions
## poly.scores<-extract.value.from.poly(poly[has.pholy,"Extra"],to.get)
## poly.scores<-gsub(")","",poly.scores) 
## poly.scores<-strsplit(poly.scores,split="\\(")
#test<-unlist(lapply(poly.scores,length))
poly.desc<-poly[has.pholy,"Feature_type"]
poly.scores<-poly[has.pholy,"Feature"]
key.pholy<-build.pholy.key(poly[has.pholy,]) # paste(poly[has.pholy,"Location"],poly[has.pholy,"Allele"],sep=":")
poly.scores<-cbind(key.pholy,poly.desc,poly.scores)

# poly.scores[1:5,]
##poly scores comtains repeats for diff trnascripts choose the most damaging
## if(dim(poly.scores)[1]>1){
## the.order<-order(as.numeric(poly.scores[,"poly.scores"]),decreasing=TRUE)
## poly.scores<-poly.scores[the.order,]   
## dups<-duplicated(poly.scores[,"key.pholy"])
## poly.scores<-poly.scores[!dups,]
## }

key.indels.poly<-build.key(indels,c("chr","start","ALT"))
## key.indels.poly[1:15]
    ## list.element.lengths<-unlist(lapply(the.genes[[annotation.labels[i]]],length))        # list.lentghs
    ## indel.index<-rep(1:length(list.element.lengths),times=list.element.lengths)           # index
    ## flat.gene.list<-unlist(the.genes[[annotation.labels[i]]])   
if(is.null(dim(poly.scores))){poly.scores<-t(as.matrix(poly.scores))}
poly.scores<-match.cols.and.collapse(list(key.indels.poly),1,poly.scores,"key.pholy",colnames(poly.scores),"delimit")
colnames(poly.scores)<-gsub("poly","regulation",colnames(poly.scores))
colnames(poly.scores)<-gsub("scores","feature",colnames(poly.scores))
poly.scores
}


extract.value.from.info<-function(info,match.string){ #will be fooled if have GGMAF=  #ok if GMAF not present need [^;]
match.string.general<-paste(match.string,"[a-zA-Z0-9,\\.]*",sep="") #regexec(";GMAF=[a-zA-Z0-9,\\.]*;",indels[has.gmaf,"INFO"]) ";" not needed
match<- regexec(match.string.general,info)
start<-unlist(match)
length<-unlist(lapply(match,function(x)  attr(x,"match.length")))
end<-start+length-1 # minus one to get last position\
start<-start+nchar(as.character(match.string))
substr(info,start=start,stop=end)
}

# extract.value.from.format(vcf.head[info.types,1],"ID=")
## info<-vcf.head[info.types,1]
## match.string<-"ID="
extract.value.from.format<-function(info,match.string){ #will be fooled if have GGMAF=  #ok if GMAF not present need [^;]
match.string.general<-paste(match.string,"[a-zA-Z0-9_\\ \\.]*",sep="") #regexec(";GMAF=[a-zA-Z0-9,\\.]*;",indels[has.gmaf,"INFO"]) ";" not needed
match<- regexec(match.string.general,info)
start<-unlist(match)
length<-unlist(lapply(match,function(x)  attr(x,"match.length")))
end<-start+length-1 # minus one to get last position
start<-start+nchar(as.character(match.string))
substr(info,start=start,stop=end)
}

extract.value.from.DS<-function(info,match.string){ #will be fooled if have GGMAF=  #ok if GMAF not present need [^;]
match.string.general<-paste(match.string,"[a-zA-Z0-9,:\\+\\ \\.\\-]*",sep="")  #regexec(";GMAF=[a-zA-Z0-9,\\.]*;",indels[has.gmaf,"INFO"]) ";" not needed
match<- regexec(match.string.general,info)
start<-unlist(match)
length<-unlist(lapply(match,function(x)  attr(x,"match.length")))
end<-start+length-1 # minus one to get last position
start<-start+nchar(as.character(match.string))
substr(info,start=start,stop=end)
}

extract.value.from.poly<-function(info,match.string){ #will be fooled if have GGMAF=  #ok if GMAF not present need [^;]
match.string.general<-paste(match.string,"[a-zA-Z0-9_\\(\\)\\ \\.]*",sep="") #regexec(";GMAF=[a-zA-Z0-9,\\.]*;",indels[has.gmaf,"INFO"]) ";" not needed
match<- regexec(match.string.general,info)
start<-unlist(match)
length<-unlist(lapply(match,function(x)  attr(x,"match.length")))
end<-start+length-1 # minus one to get last position
start<-start+nchar(as.character(match.string))
substr(info,start=start,stop=end)
}


allele.summary.with.apply<-function(x){
  missing<-is.na(x) | x=="NA,NA" | x=="NA" | !grepl(",",x) # not an NA must have a comma
  x<-x[!missing]
  num.hetero<-length(x)
  if(num.hetero!=0){
  x<-as.numeric(unlist(strsplit(x,split=",")))
  dim(x)<-c(2,num.hetero) # x contains ref allele counts in first row
  the.counts<-apply(x,1,function(xx) sum(xx,na.rm=TRUE))
  c(the.counts[2],the.counts[1],signif((the.counts[2]/(the.counts[2]+the.counts[1]))*100,digits=3))
  }else{
    c(0,0,0) ## all values missing
  }
  
}
  
## x<-as.matrix(allele.depths.group)[,1]
allele.summary.individuals<-function(x){  ### work on COLUMNS
#  print(x)
  missing<-is.na(x) | x=="NA,NA" | x=="NA" | !grepl(",",x) # not an NA must have a comma
  het<-rep(NA,times=length(x))
  x<-x[!missing]
  x<-as.numeric(unlist(strsplit(x,split=",")))
  dim(x)<-c(2,sum(!missing)) 
  ref.hom<-x[2,]==0
  alt.hom<-x[1,]==0
  x<-x[2,]/(x[2,]+x[1,])

  het[!missing]<-x
  het[!missing][ref.hom]<-0
  het[!missing][alt.hom]<-1
  het
}

## x<-as.matrix(allele.depths.group)[,2]
allele.DP.from.AN<-function(x){  ### work on COLUMNS
#  print(x)
  missing<-is.na(x) | x=="NA,NA" | x=="NA" | !grepl(",",x) # not an NA must have a comma
  ## het<-rep(NA,times=length(x))
 #
  x[!missing]<-as.numeric(unlist(lapply(strsplit(x[!missing],split=","),function(x) sum(as.numeric(x)) ) ))
  x
}

ref.inversion.for.summary<-function(ref.inversion,summary.geno.group){
  num.inversions<-sum(ref.inversion)
  if(num.inversions>0){
  the.col.names<-colnames(summary.geno.group)
  alt.label<-the.col.names[grep("ALT.Alleles",the.col.names)]
  alt.homo.label<-the.col.names[grep("ALT_HOMO",the.col.names)]
  ref.label<-the.col.names[grep("REF.Alleles",the.col.names)]
  geno.label<-the.col.names[grep("GENO",the.col.names)]
  maf.label<-the.col.names[grep("MAF",the.col.names)]
  summary.geno.group[ref.inversion,maf.label]<- 1.0-as.numeric(summary.geno.group[ref.inversion,maf.label])
  tmp<- summary.geno.group[ref.inversion,ref.label]
  summary.geno.group[ref.inversion,ref.label]<- summary.geno.group[ref.inversion,alt.label]
  summary.geno.group[ref.inversion,alt.label]<-tmp
  tmp<-strsplit(summary.geno.group[ref.inversion,geno.label],split=",")
  
  tmp<-unlist(tmp)
  dim(tmp)<-c(3,num.inversions)
  
  tmp<-t(tmp)
  summary.geno.group[ref.inversion,alt.homo.label]<-tmp[,3] # make ALT.homo ref.homo
  tmp<-paste(tmp[,3],tmp[,2],tmp[,1],sep=",")
  summary.geno.group[ref.inversion,geno.label]<-tmp
}
   summary.geno.group
}

## sample.files<- sample.files[isamp]
prepare.for.Vcf.file.read<-function(sample.files){
print(paste("Doing samples",sample.files,sep=" "))
setwd(eval(as.name(paste(names(sample.files),"dir",sep=".")))) ## different directory for each indel type

######BELOW PROCESSING  this for snp for indel varient types in vcf4.0 or vcf 3.0 format

### get location of header:
chromo1<-try(scan(sample.files,what=character(),n=5000,sep="\n",skip=0,fill=TRUE,na.strings="",quote="\"")) ## find the start of the vcf file
skip.lines<-grep("^#CHROM",chromo1)
if(length(skip.lines)>1){print("ERROR multiple chrom lables found");skip.lines<-skip.lines[1]}
skip.lines<-skip.lines
 
options(show.error.messages = TRUE)
column.labels<-read.delim(sample.files,header=F,nrows=1,skip=(skip.lines-1),sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
num.vars<-dim(column.labels)[2]

vcf.head<-read.delim(sample.files,header=F,nrows=(skip.lines-1),skip=0,sep="\t",fill=TRUE,stringsAsFactors=FALSE,quote="\"")
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


column.labels[match("#CHROM",column.labels)]<-"chr"

out.list<-list(column.labels,skip.lines,num.vars,info.types,info.class,info.description,format.types,format.class,format.description)
names(out.list)<-c("column.labels","skip.lines","num.vars","info.types","info.class","info.description","format.types","format.class","format.description")
out.list
}
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



## indels[1:5,]

## process.GATK.indels(indels[1:50,],sample.files[isamp],vcf.type,format.types,info.types,info.class,num.cores)
##  ## indels<-  process.GATK.indels(indels,sample.files[isamp],vcf.type,format.types,info.types,info.class,num.cores)
##  indels<-indels.ori[1:50,]
##  the.sample.file <- sample.files[isamp]
##  vcf.type
##  format.types
##  info.types
##  info.class
##  num.cores

 process.GATK.indels<-function(indels,the.sample.file,vcf.type,format.types,info.types,info.class,num.cores){
#   print(dim(indels))
column.labels <- colnames(indels)
#######################Process the INFO section common to all indels
 ###
#print(num.cores)
 ###
 if(vcf.type=="v3"){
   hetero<-indels[,dim(indels)[2]]=="0/1"
  indels[!hetero,"INFO"]<- paste("AB=1",indels[!hetero,"INFO"],sep=";")
                }

  ########### to process using lists  use this method timing:  52.91s    0.05  ( {   } )


## system.time(
##  {

####multicore way WORKS! safe but slower  74.410  16.490  21.155 compare to  11.110   0.010  11.124 
## info<-strsplit(indels[,"INFO"],split=";")   
## key.indel<-build.key(indels,c("chr","POS","REF","ALT"))
## names(info)<-key.indel
## info<-mclapply(info,process.info.list,mc.cores=num.cores)
## posns<-match(key.indel,names(info))
## info<-info[posns]
## num.info<-length(info.types)
## len.info<-length(info)
## info<-unlist(info)
## dim(info)<-c(num.info,len.info)  
## info<-t(info) # info[1:5,]
## colnames(info)<-info.types

   ########vector.method WORKS in # 11.110   0.010  11.124 
info<-indels[,"INFO"]
#if(!skip.info.prosessing){  # skipping info processing causes error in fix multi alleles since that used AC anf AF
 
flag.items <- info.types[info.class=="Flag"]
 for(iflag in 1:length(flag.items)){
info<-gsub(paste(flag.items[iflag],";",sep=""),paste(flag.items[iflag],"=1;",sep=""),info) # is flag present just adds a one
}


info<-gsub("=",";",info) # get read to split
info<-strsplit(info,split=";") # unique(unlist(lapply(info,length)))
info.length<-unlist(lapply(info,length))

#abs(c(1,1,2.1,-2.6,0.001)-trunc(c(1,1,2.1,-2.6,0.001))
if(sum(abs(info.length-trunc(info.length))!=0)!=0){print("Error in processing INFO unpaired values")}

info<-unlist(info)
wanted<-seq(from=2, to=length(info),by=2)
wanted.labels<-seq(from=1, to=length(info),by=2)
  
info.labels<-info[wanted.labels]
info.num<-info[wanted]
info.length<-info.length/2
info.index<-rep(1:length(indels[,"INFO"]),times=info.length)
info<-matrix(data=NA,nrow=length(info.length),ncol=length(info.types))
colnames(info)<-info.types

# i.info.col<-1
for(i.info.col in 1:length(info.types)){
  the.col<-info.types[i.info.col]
  posns.in.flatten.info<-grep(paste("^",the.col,"$",sep=""),info.labels,fixed=FALSE) ## can have AC or AC1 etc correct make 01/21/2013
  posns<-info.index[posns.in.flatten.info]
  info[posns,the.col]<-info.num[posns.in.flatten.info]
}
                                        # numerical data wanted
#} # skip info 

## }) # system time test


##  info[1:5,]
## to.test<-15
## a.na<-is.na(info.ori[,to.test]) # | (is.na(info[,to.test]) | is.na(info.ori[,to.test]))
## sum(a.na)
## sum(  is.na(info[,to.test]))
## sum(  is.na(info.ori[,to.test]))
## sum(  info[!a.na,to.test]!=info.ori[!a.na,to.test]) 

####  Alternate mulitcore fails
##   library(foreach)  ## could not get this to work well
##    library(doMC)
## registerDoMC(cores=7)
## foreach(ifi=1:length(info),.combine=c,.multicombine=TRUE,.inorder=TRUE) %dopar% process.info.list(info[[ifi]],info.types)
###
  

 ########### to process using arrays use this method timing:  73.070s  ## may be abl to make multicore  
## system.time(
## {
## info<- apply(as.matrix(indels[,"INFO"]),1,process.info,info.types)
## info<-t(info) # info[1:5,]
## colnames(info)<-info.types
## }
## )
# info[1:5,] # info contains the common metrics for that location for all samples 
############################################

#################### Process the FORMAT sectionfor each sample
### assumes samples are listed after FORMAT COLUMN
info.sample<-{}
format.posn<-match("FORMAT",colnames(indels))

samples.processing<-length(column.labels)-format.posn # number of sample in the vcf file
samples.order.in.ALL<-column.labels[(format.posn+1):length(column.labels)] # samlple labels
  
#system.time({
format.string<-indels[1,format.posn]
complex.format<-sum(indels[,format.posn]!=format.string) !=0  # TRUE if all format are not the same order
format.string<-unlist(strsplit(format.string,split=":")) # used for column names below

if(complex.format){
  format.string<-indels[,format.posn]
  print("Comples Format found in vcf file")
}else{
  format.types<-format.string ## have same format labels but make sure in the correct order
} #use all possible attributes the function will then give the correct order




#  test2<-info.parse(indels[,c(format.posn+1)],format.types,format.string,complex.format) # protoy
# iter(indels[,c((format.posn+1):(format.posn+samples.processing))], by='col')
  
  
########Multicore run with foreach ## Done columnwise and Format the same for each column
print(paste("START processing",samples.processing,"samples",sep=" "))  
info.sample<-foreach(info2=iter(indels[,c((format.posn+1):(format.posn+samples.processing))], by='col',chunksize=1), .combine=cbind,.multicombine=TRUE,.inorder=TRUE) %dopar% info.parse(info2,format.types,format.string,complex.format)
###
print(paste("Finished processing",samples.processing,"samples",sep=" "))  

###construct the column names

if(vcf.type=="v3") {  
#  if(grepl(combined.extension,the.sample.file)){
  sample.class<-as.character( t(matrix(data=rep(samples.order.in.ALL,times=length(c("GT"))),nrow=length(samples.order.in.ALL),ncol=length(c("GT")),byrow=FALSE)) )
  colnames(info.sample)<-paste(sample.class,c("GT"),sep=".")
#}

#if(!grepl(combined.extension,the.sample.file)){ colnames(info.sample)<-c("GT")  }
  
}else{ # modern vcf 4.0 format
  
#if(grepl(combined.extension,the.sample.file)){
  sample.class<-as.character( t(matrix(data=rep(samples.order.in.ALL,times=length(format.types)),nrow=length(samples.order.in.ALL),ncol=length(format.types),byrow=FALSE)) )
  colnames(info.sample)<-paste(sample.class,format.types,sep=".")
#}

#if(!grepl(combined.extension,the.sample.file)){ colnames(info.sample)<-format.types  }

}

#################### FINISHED processing individual filters   indels<-indels.ori[c(1000,grep("1/2",indels.ori[,10])[1:3], grep("/3",indels.ori[,10])) ,]
TYPE<-rep(names(the.sample.file),times=dim(indels)[1])
filter.posn<-match("FILTER",colnames(indels))
  
indels<-cbind(subset(indels,select=1:filter.posn),TYPE,info,info.sample) ## use subset here incase indels has only 1 element

  
rm(TYPE)
rm(info)
rm(info.sample)
#colnames(data)<-c( "chr","start","end","strand",c("ID","REF","ALT","QUAL","FILTER"),"TYPE",colnames(info),colnames(info.sample))
############################### MAKE the final data file
  print("Have data now fix multi-alleles")

 # the.gen<-expand.labels.to.samples(c("GT"),samples.order.in.ALL)
############################## multi-allele genotypes need to fix
# poly position typicaly have no SB or HRun given

 #####################################################################uses just indels from this point ###################################

########Multicore run with foreach  #### WARNING PROBLEM IS LESS THAN num
if(dim(indels)[1]<1000){num.bits<-1}else{num.bits<-num.cores}  #### WARNING PROBLEM IS LESS indels than num.cores
print(paste("START processing muli-alleles",samples.processing,"samples",sep=" "))  
## indels<-foreach(indels.bit=iter(indels, by='row',chunksize=as.integer(dim(indels)[1]/num.bits)), .combine=rbind,.multicombine=TRUE,.inorder=TRUE) %dopar% correct.multi.alleles.flattern(indels.bit,samples.order.in.ALL)

indels<-foreach(indels.bit=iter(indels, by='row',chunksize=as.integer(dim(indels)[1]/num.bits)), .combine=rbind,.multicombine=TRUE,.inorder=TRUE) %dopar% correct.multi.alleles(indels.bit,samples.order.in.ALL)
###
print(paste("Finished processing muli-alleles",samples.processing,"samples",sep=" "))


########################### PUT ALLELES IN ANNOVAR FORMAT from convert2annovar.pl line 1083########

indels[is.na(indels[,"REF"]),"REF"]<-"NA" ## see a REF =N ALT= NA R confuses with "NA" with NA - kills substr below
indels[is.na(indels[,"ALT"]),"ALT"]<-"NA"


ref.length<-nchar(as.character(indels[,"REF"]))
alt.length<-nchar(as.character(indels[,"ALT"]))
is.snp<-(ref.length==1 & alt.length==1)

# sum(!is.snp)

POS.end<-as.numeric(indels[,"POS"])
del<-ref.length > alt.length
ins<-(ref.length <= alt.length) & !is.snp
## indels[del,][1:5,]
## POS.end[del][1:5]
### deletion or block substitution
head<-substr(as.character(indels[del,"REF"]),1,alt.length[del])
head.is.mut<-(head==as.character(indels[del,"ALT"]))
indels[del,"REF"][head.is.mut]<-substr(as.character(indels[del,"REF"][head.is.mut]),(alt.length[del][head.is.mut]+1),ref.length[del][head.is.mut])
indels[del,"ALT"][head.is.mut]<-"-"
indels[del,"POS"][head.is.mut]<-as.numeric(indels[del,"POS"][head.is.mut]) + nchar(as.character(head[head.is.mut]))
POS.end[del]<-POS.end[del]+ref.length[del]-1  # same for both head is mut and not head is mut

## indels
## POS.end
### insertion or block substitution
head<-substr(as.character(indels[ins,"ALT"]),1,ref.length[ins])
head.is.ref<-(head==as.character(indels[ins,"REF"]))
indels[ins,"ALT"][head.is.ref]<-substr(as.character(indels[ins,"ALT"][head.is.ref]),(ref.length[ins][head.is.ref]+1),alt.length[ins][head.is.ref])
indels[ins,"REF"][head.is.ref]<-"-"
indels[ins,"POS"][head.is.ref]<-as.numeric(indels[ins,"POS"][head.is.ref]) + ref.length[ins][head.is.ref]-1
POS.end[ins]<-POS.end[ins]+ref.length[ins]-1
########################################################################################################
#$$$  all.alleles
  

POS.end<-as.integer(POS.end) ### else can get 8e+7 type errors
begining.cols<-c("chr","POS")
strand<-rep("+",times=dim(indels)[1])
## indels<-cbind(indels[,begining.cols],POS.end,strand,indels[,colnames(indels)[!(colnames(indels) %in% begining.cols)]])
indels<-cbind(subset(indels,select=begining.cols),POS.end,strand,subset(indels,select=colnames(indels)[!(colnames(indels) %in% begining.cols)])) ##use incase indels 1 element long


colnames(indels)[colnames(indels)=="POS.end"]<-"end"
colnames(indels)[colnames(indels)=="POS"]<-"start"
  
rm(POS.end)
rm(strand)
indels[indels[,"chr"]=="chrMT","chr"]<-"chrM"  ### annotation algorithms expect MT
indels

   } # end function  process.GATK.indels




test.for.coding.type<-function(geneanno.table,transcripts,coding.type){
  #geneanno.table assumed has colums like: "refGene::location"
  #transcripts like:    refGene   knownGene     ensGene 
  #                    "refGene" "knownGene"   "ensGene"
  # coding.type like:  "nonsynonymous SNV"          "stopgain SNV" 
  
coding.type.found<-rep(FALSE,times=dim(geneanno.table)[1])
for(i in 1:length(transcripts)){
  the.test<-test.wanted.mutation(geneanno.table[,paste(transcripts[i],"location",sep="::")],coding.type,delimit.by=";")
  the.test[is.na(the.test)]<-FALSE # if no annotion set to FALSE for wanted   sum(the.test[is.na(the.test)])
  coding.type.found<- coding.type.found | the.test
}
 coding.type.found
} # end function test.for.coding.type



genotype.summary<-function(x){
  if(is.null(dim(x))){x<-matrix(data=x,nrow=1,ncol=length(x))}
  hetero.geno<-x=="0/1"
  ref.geno<- x=="0/0"
  alt.geno<-x=="1/1"
  alt.homo.counts<-apply(alt.geno,1,function(xx) sum(xx,na.rm=TRUE))
  hetero.counts<-apply(hetero.geno,1,function(xx) sum(xx,na.rm=TRUE))
  ref.homo.counts<-apply(ref.geno,1,function(xx) sum(xx,na.rm=TRUE))
  alt.counts<- alt.homo.counts*2+hetero.counts
  ref.counts<- ref.homo.counts*2+hetero.counts
  total.alleles<-ref.counts+alt.counts
  maf<-alt.counts/total.alleles
  genos.sum<-paste(alt.homo.counts,hetero.counts,ref.homo.counts,sep=",")
  cbind(signif(maf, digits = 6),alt.counts,ref.counts, total.alleles,(2*dim(x)[2]-total.alleles),alt.homo.counts,hetero.counts,genos.sum)
}
  
  
genotype.summary.with.apply<-function(x){
  counts<-tapply(x,x,length)
  missing.counts<-counts["NA"]
  genos<-counts[c("1/1","0/1","0/0")]
  genos[is.na(genos)]<-0 # case does not have that genotype
  alt.counts<-genos[1]*2+genos[2]
  ref.counts<-genos[2]+genos[3]*2
  total.alleles<-ref.counts+alt.counts
  maf<-alt.counts/total.alleles
  genos.sum<-paste(genos,collapse=",")
  c(signif(maf, digits = 6),alt.counts,ref.counts,total.alleles,missing.counts,genos[1],genos[2],genos.sum)
                              }

to.info.format<- function(x,labels=colnames(x)){ ##labels<- colnames(x) if not otherwise specified
  #print(labels)
  new.x<-""
  if(is.null(dim(x))){x<-as.matrix(x)}
 new.x<-paste(paste(labels[1],x[,1],sep="="),new.x,sep="")
 if(length(labels)>1){
 for ( i in 2:(length(labels))){
  new.x<-paste(paste(labels[i],x[,i],sep="="),new.x,sep=";")
}}
 new.x
}
#answer = which.isMatchingAt("N", seqs, at=6:1, follow.index=TRUE)




#test2<-info.parse(indels[,c(format.posn+1)],format.types,format.string,complex.format)
## info2<-indels[,c(format.posn+44)] # test<-indels[,c(format.posn+1)] # indels[385:390,c(format.posn+44)]
## format.types
## format.string<-format.string
## complex.format
## test<-strsplit(test,split=":")
## test<-unlist(lapply(test,function(x) x[1]))
## test[test=="./."]<-"NA"
## sum(test !=info2[,"GT"]) # & !(is.na(test) & is.na(info2[,"GT"]) ))

info.parse<-function(info2,format.types,format.string,complex.format=FALSE){
#print(complex.format)
if(vcf.type=="v3") {
no.reads<-grepl("./.",info2)  # info2=="./." # no reads so no information sum(no.reads) ## "./.:.:1" oct GATK recal versions 2.7-2-g6bda569
info2[no.reads]<-"NA/NA"
info2<-as.matrix(info2)
colnames(info2)<-c("GT")

}else{

#info2<-apply(info2,1,process.format) # original way

####multicore way
if(complex.format){  ## if complex format the the format column is not always the same

## info2<-cbind(format.string,info2)
## info2<-apply(info2,1,list)
## key.indel<-1:dim(info2)[1]
## names(info2)<-key.indel
## info2<-mclapply(info2,process.format.list,mc.cores=1)
## posns<-match(key.indel,names(info2))
## info2<-info2[posns]

############use same trick used to process info in process.GATK.indels

format.string<-strsplit(format.string,split=":") # unique(unlist(lapply(info,length)))
format.length<-unlist(lapply(format.string,length))


no.reads<-(grepl("./.",info2,fixed=TRUE) | (info2=="."))  #  (info2=="./.") | (info2==".") # no reads so no information sum(no.reads) changes Nov 4 2013 seen: "./.:.:1"
info2[no.reads]<-unlist(lapply(format.length[no.reads],function(x){ paste(rep(NA,times=x),collapse=":")}))


# WARNING info2 The format string can have MORE ITEMS than in sample string in merged VCF files
### fix this below make sure have the same length
format.num<-strsplit(info2,split=":")
format.num.length<-unlist(lapply(format.num,length))

## sum(format.num.length !=format.length) #
if(sum(format.num.length >format.length)){print("ERROR ERROR format and genotype mismatch")}
to.fix<-format.num.length !=format.length ## these have all got more formats than lengths
if(sum(to.fix)>0){
  to.strip<-format.length[to.fix] -format.num.length[to.fix]
  to.strip.types<-tapply(to.strip,to.strip,length)
  for(is in 1:length(to.strip.types)){
    
    the.strip.type<-as.integer(names(to.strip.types)[is])
    to.fix.subset<-to.strip==the.strip.type
    new.format.string<-lapply(format.string[to.fix][to.fix.subset],function(x,the.strip.type){x[1:(length(x)-1)]})
    format.string[to.fix][to.fix.subset]<-new.format.string
    format.length[to.fix][to.fix.subset]<-format.num.length[to.fix][to.fix.subset]
  } #strip.types many need to take off more than one element
} # to.fix required

########################


format.labels<-unlist(format.string)
format.num<-unlist(strsplit(info2,split=":"))

format.index<-rep(1:length(info2),times=format.length)
info2<-matrix(data=NA,nrow=length(info2),ncol=length(format.types))
colnames(info2)<-format.types

# i.info.col<-5
for(i.info.col in 1:length(format.types)){
  the.col<-format.types[i.info.col]
  posns.in.flatten.info<-grep(paste("^",the.col,"$",sep=""),format.labels,fixed=FALSE) ## can have AC or AC1 etc correct make 01/21/2013
  posns<-format.index[posns.in.flatten.info]
  info2[posns,the.col]<-format.num[posns.in.flatten.info]
}


}else{ # not complex format
no.reads<-(grepl("./.",info2,fixed=TRUE) | (info2=="."))  #   (info2=="./.") | (info2==".") # no reads so no information sum(no.reads) changes Nov 4 2013 seen: "./.:.:1"
info2[no.reads]<-gsub(", ",":",toString(rep("NA",times=length(format.types)))) #  this will be the longest required 
info2<-strsplit(info2,split=":")  # fastest  3.890   0.010   3.894
num.info2<-length(format.types)
len.info2<-length(info2)
info2<-unlist(info2)
dim(info2)<-c(num.info2,len.info2)  
info2<-t(info2)
colnames(info2)<-format.types
} ## not complex format

} ## not v3 VCF
info2
} ### end sub info.process
# sum(info2.ori[,1] != info2[,1])

genetic.sample.QC<-function(sample.index,genotypes,p,position.filter,all.sample.labels,on.X,parX){

  ## this routine required a sample index
##   sample.index<-matrix(data=NA,nrow=(length(all.sample.labels)*(length(all.sample.labels)+1)/2),ncol=2)
## k<-1
## for(i in 1:length(all.sample.labels)){  
##   for(j in i:length(all.sample.labels)){
##     sample.index[k,]<-c(i,j)
##     k<-k+1
## }}

  
print(dim(sample.index))

wanted.att<-c("sample_A","sample_B","IBS","sex_Predicted","sex.code_A","Het_chrX_A","ChrX-Het:Homo","Num_Good_SNPs_A")
related<-matrix(data=NA,nrow=dim(sample.index)[1],ncol=length(wanted.att))
colnames(related)<-wanted.att


for(k in 1:dim(sample.index)[1]){  

    
    ## all.sample.labels[i]
    ## all.sample.labels[j]
    i<-sample.index[k,1]
    j<-sample.index[k,2]

 #   print(paste(i,j,sep=":"))
  
    present<-!(is.na(genotypes[,i]) | is.na(genotypes[,j]))        
    ok<-present & position.filter
    
#sum(ok)
## y<-as.numeric(as.numeric(genotypes[ok,i])) 
## genotypes[,i]
## tapply(genotypes[,i],genotypes[,i],length)   
## genotypes[ok,c(i,j)]
  if(i==j){
   ## tapply(indels[ok,"chr"],indels[ok,"chr"],length)
   x.hetro<-sum(genotypes[ok & on.X & !parX,i]==1)
   x.homo<-sum(genotypes[ok & on.X & !parX,i]==2)
   x.het<-x.hetro/(x.hetro+x.homo)
   sex=2;sex.code="F";het.count<-paste(x.hetro,x.homo,sep=":")
   if(!is.finite(x.het)){sex=9;sex.code="Unknown"}else{
     if(x.het<0.3){sex=1;sex.code="M"}
     if(x.het<0.3){sex=1;sex.code="M"}
   }
  Aijk<- (as.numeric(genotypes[ok,i])^2 -(1.0 + 2*p[ok])*as.numeric(genotypes[ok,i]) + 2*p[ok]^2)/(2*p[ok]*(1-p[ok]))


 ## Aijk<- (as.numeric(genotypes[ok,i])^2 -(1.0 + 0)*as.numeric(genotypes[ok,i]) + 0)/(2*p[ok]*(1-p[ok]))


  Aij<-  1.0+sum(Aijk,na.rm=TRUE)/sum(!is.na(Aijk))
  Aij<-  1.0+(1-(1.0/sum(!is.na(Aijk))/var(Aijk,na.rm=TRUE)))*(Aij-1.0)

  }else{
    
  Aijk<- (as.numeric(genotypes[ok,i])-2*p[ok])*(as.numeric(genotypes[ok,j])-2*p[ok])/(2*p[ok]*(1-p[ok]))
  Aij<-  sum(Aijk,na.rm=TRUE)/sum(!is.na(Aijk))
  Aij<-  (1-(1.0/sum(!is.na(Aijk))/var(Aijk,na.rm=TRUE)))*Aij
}

    ## sum(ok)
    ## Aijc("sample_A","sample_B","IBS","sex_Predicted","sex.code_A","Het_chrX_A","ChrX-Het:Homo,""Num_Good_SNPs_A","Num_POLYmorph_SNPs_A")
    
  related[k,] <- c(all.sample.labels[i],all.sample.labels[j],signif(Aij,digits=3),sex,sex.code,signif(x.het,digits=3),het.count,sum(ok))
#    k<-k+1
  }
related
}
###################################################################


genetic.sample.QC.accumilate<-function(genotypes,p,position.filter,all.sample.labels,on.X,parX){

# position.filter<-position.filter.QC  
# genotypes<-all.genotypes$"genotypes"
## print(dim(genotypes))
## print(length(p))

if(is.null(dim(genotypes))){  # just in case not a matrix or data.frame
ncol<-length(genotypes)
nrow<-1
dim(genotypes)<-c(nrow,ncol)
}

Aij<-matrix(data=0,nrow=(length(all.sample.labels)+length(all.sample.labels)),ncol=(length(all.sample.labels)+2+length(all.sample.labels)))
colnames(Aij)<-c(all.sample.labels,"x.hetro","x.homo",paste(all.sample.labels,"C",sep="_"))
rownames(Aij)<-c(all.sample.labels,paste(all.sample.labels,"C",sep="_"))
#Aij

for(i in 1:length(all.sample.labels)){  # columns
  for(j in i:length(all.sample.labels)){  # rows
    
    ii<-i+2+length(all.sample.labels) # columns
    jj<-j+length(all.sample.labels)  # rows
    
     ## print(paste(i,j,sep=":"))
     ## print(paste(rownames(Aij)[j],colnames(Aij)[i],sep=" : "))
     ## print(paste(rownames(Aij)[jj],colnames(Aij)[ii],sep=" : "))
#   test=paste(rownames(Aij)[j],colnames(Aij)[i],sep=" : ") 

    present<-!(is.na(genotypes[,i]) | is.na(genotypes[,j]))        
    ok<-present & position.filter
    

  if(i==j){
   ## tapply(indels[ok,"chr"],indels[ok,"chr"],length)
   Aij[j,"x.hetro"]<-sum(genotypes[ok & on.X & !parX,i]==1)
   Aij[j,"x.homo"]<-sum(genotypes[ok & on.X & !parX,i]==2)
   Aijk<- (as.numeric(genotypes[ok,j])^2 -(1.0 + 2*p[ok])*as.numeric(genotypes[ok,j]) + 2*p[ok]^2)/(2*p[ok]*(1-p[ok]))

   
   Aij[j,i]<- sum( Aijk ,na.rm=TRUE)  #sum(Aijk,na.rm=TRUE)
   Aij[jj,ii]<- sum(!is.na(Aijk))
  }else{

  Aijk<- (as.numeric(genotypes[ok,i])-2*p[ok])*(as.numeric(genotypes[ok,j])-2*p[ok])/(2*p[ok]*(1-p[ok]))
  Aij[j,i]<- sum( Aijk ,na.rm=TRUE)  #sum(Aijk,na.rm=TRUE)
  Aij[jj,ii]<- sum(!is.na(Aijk))
  }

}}

Aij
}
###################################################################



sum.QC.matrix<-function(Aij,all.sample.labels){

wanted.att<-c("sample_A","sample_B","IBS","sex_Predicted","sex.code.Predicted","Het_chrX_A","ChrX-Het:Homo","Num_Good_SNPs_A")
related<-matrix(data=NA,nrow=(length(all.sample.labels)*(length(all.sample.labels)+1)/2),ncol=length(wanted.att))
colnames(related)<-wanted.att # set the columns of related; "sample_A","sample_B","IBS","sex_Predicted" can't be modified they are used in "Apply inheritance QC_MULTI_NEW.r "

k<-1
  for(i in 1:length(all.sample.labels)){  # columns
  for(j in i:length(all.sample.labels)){  # rows
    
    ii<-i+2+length(all.sample.labels) # columns
    jj<-j+length(all.sample.labels)  # rows
    
  if(i==j){
   ## tapply(indels[ok,"chr"],indels[ok,"chr"],length)
   x.het<- Aij[j,"x.hetro"]/( Aij[j,"x.hetro"]+ Aij[j,"x.homo"])
   sex=2;sex.code="F";het.count<-paste(Aij[j,"x.hetro"],Aij[j,"x.homo"],sep=":")
   if(!is.finite(x.het)){sex=9;sex.code="Unknown"}else{
     if(x.het<0.3){sex=1;sex.code="M"}
     if(x.het<0.3){sex=1;sex.code="M"}
   }


  A<-  1.0+ Aij[j,i]/Aij[jj,ii]
#  Aij<-  1.0+(1-(1.0/sum(!is.na(Aijk))/var(Aijk,na.rm=TRUE)))*(Aij-1.0)

  }else{
    
  A<-  Aij[j,i]/Aij[jj,ii]
#  Aij<-  (1-(1.0/sum(!is.na(Aijk))/var(Aijk,na.rm=TRUE)))*Aij
}

    ## sum(ok)
    ## Aijc("sample_A","sample_B","IBS","sex_Predicted","sex.code_A","Het_chrX_A","ChrX-Het:Homo,""Num_Good_SNPs_A","Num_POLYmorph_SNPs_A")
    
  related[k,] <- c(all.sample.labels[i],all.sample.labels[j],signif(A,digits=3),sex,sex.code,signif(x.het,digits=3),het.count,Aij[jj,ii])
    k<-k+1
  }}

related
}

## het.read.thresh<-1  # depth must exceed this to apply filer on genotypes
## het.low<-0
## het.high<-1.1
## het.low.ref<-0
## het.high.alt<-0.95
## depth.het.low<-10
## depth.homo.low<-5


filtered.genotype.summary<-function(indels,the.samples,prefix="",suffix="",het.read.thresh=20,het.low=0.2,het.high=0.80,het.low.ref=0.05,het.high.alt=0.95,depth.het.low=7,depth.homo.low=5){
########################### ALL ###############
## applied filters to genotypes and then makes a summary
## prefix<-""
## suffix<-".ALL"
## het.read.thresh<-20  # depth must exceed this to apply filer on genotypes
## het.low<-0.20
## het.high<-0.80
## het.low.ref<-0.05
## het.high.alt<-0.95
## depth.het.low<-10
## depth.homo.low<-5
### the.samples<-all.sample.labels
#the.samples
#if(prefix!=""){prefix<-paste(prefix,".",sep="")}
if(is.null(rownames(indels))){core.ann<-c("chr","start","end","REF","ALT","TYPE") ;rownames(indels)<-build.key(indels,core.ann,add.chr.label=FALSE)}

key.indels<-rownames(indels)


allele.depths.group<-indels[,paste(the.samples,"AD",sep=".")]
summary.het.indiv<-apply(as.matrix(allele.depths.group),2,allele.summary.individuals)
colnames(summary.het.indiv)<-gsub(".AD$",".HET",colnames(summary.het.indiv))


if( sum(paste(the.samples,"DP",sep=".") %in% colnames(indels))==length(the.samples)){
  allele.depths<-as.matrix(indels[,paste(the.samples,"DP",sep=".")])
   }else{
  allele.depths<-apply(as.matrix(allele.depths.group),2,allele.DP.from.AN)
  colnames(allele.depths)<-gsub(".AD",".DP",colnames(allele.depths))
}


#allele.depths<-as.matrix(indels[,paste(the.samples,"DP",sep=".")])

#allele.depths[5030,]
ncol<-length(the.samples)
nrow<-dim(indels)[1]
allele.depths<-as.numeric(allele.depths)
allele.depths[is.na(allele.depths)]<-0

if(is.null(dim(allele.depths))){
dim(allele.depths)<-c(nrow,ncol)
## allele.depths<-t(allele.depths)
colnames(allele.depths)<-paste(the.samples,"DP",sep=".")
}


genotypes<-as.matrix(indels[,paste(the.samples,"GT",sep=".")])
if(sum(colnames(summary.het.indiv)!=paste(the.samples,"HET",sep="."))!=0){summary.het.indiv<-summary.het.indiv[,paste(the.samples,"HET",sep=".")]}
if(sum(colnames(allele.depths)!=paste(the.samples,"DP",sep="."))!=0){summary.het.indiv<-summary.het.indiv[,paste(the.samples,"HET",sep=".")]}

## tapply(genotypes[,1],genotypes[,1],length)
passing.genos<- (genotypes=="1/1" & (( summary.het.indiv > het.high.alt) | allele.depths < het.read.thresh) &  allele.depths > depth.homo.low &  grepl("^snp",indels[,"TYPE"])) |
      (genotypes=="0/1" & ((summary.het.indiv >= het.low & summary.het.indiv <= het.high) | allele.depths < het.read.thresh)  &  allele.depths > depth.het.low &  grepl("^snp",indels[,"TYPE"])) |
      (genotypes=="0/0" & ((summary.het.indiv < het.low.ref) | allele.depths < het.read.thresh ) &  allele.depths > depth.homo.low  &  grepl("^snp",indels[,"TYPE"]))  |
      (genotypes=="1/1" &  allele.depths > depth.homo.low &  grepl("^indel",indels[,"TYPE"])) |
      (genotypes=="0/1" &  allele.depths > depth.het.low &  grepl("^indel",indels[,"TYPE"])) |
      (genotypes=="0/0" &  allele.depths > depth.homo.low &  grepl("^indel",indels[,"TYPE"]))


## test<- (genotypes=="0/0")& (summary.het.indiv < het.low.ref)

## &  allele.depths > depth.homo.low  &  grepl("^snp",indels[,"TYPE"]))
## apply(passing.genos,1,sum)

## allele.depths >= depth.homo.low

genotypes[!passing.genos]<-NA  ### exclude failing genotypes from the calculation

summary.geno.group<-genotype.summary(genotypes)
colnames(summary.geno.group)<-paste(prefix,c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),suffix,sep="")
rownames(summary.geno.group)<-key.indels

#genotypes[genotypes!="0/1"]<-NA

allele.depths.group[genotypes!="0/1" | genotypes=="NA" | is.na(genotypes) ] <- NA  # dim(allele.depths.group)
summary.depths.group<-apply(allele.depths.group,1,allele.summary.with.apply)
summary.depths.group<-t(summary.depths.group)
colnames(summary.depths.group)<-paste(prefix,c("Hetero.ALT.reads","Hetero.REF.reads","Hetero.Read.Balance"),suffix,sep="")
rownames(summary.depths.group)<-key.indels

colnames(summary.het.indiv)<-paste(prefix,colnames(summary.het.indiv),sep="")
rownames(summary.het.indiv)<-key.indels

 list(summary.geno.group=summary.geno.group,summary.depths.group=summary.depths.group,genotypes=genotypes)

}


#filtered.genotype.summary(a.indel[41458:41460,],fam[,1],prefix="",suffix="",20,0.2,0.80,0.05,0.95,10,5)

filtered.genotype<-function(indels,the.samples,prefix="",suffix="",het.read.thresh=20,het.low=0.2,het.high=0.80,het.low.ref=0.05,het.high.alt=0.95,depth.het.low=7,depth.homo.low=5){
########################### ALL ###############
## applied filters to genotypes and then makes a summary
## indels<-  a.indel[41458:41460,]
## prefix<-""
## suffix<-""
## het.read.thresh<-20  # depth must exceed this to apply filer on genotypes
## het.low<-0.20
## het.high<-0.80
## het.low.ref<-0.05
## het.high.alt<-0.95
## depth.het.low<-10
## depth.homo.low<-5
## the.samples<-fam[,1]
#the.samples
#if(prefix!=""){prefix<-paste(prefix,".",sep="")}
if(is.null(rownames(indels))){core.ann<-c("chr","start","end","REF","ALT","TYPE") ;rownames(indels)<-build.key(indels,core.ann,add.chr.label=FALSE)}

key.indels<-rownames(indels)


allele.depths.group<-indels[,paste(the.samples,"AD",sep=".")]
summary.het.indiv<-apply(as.matrix(allele.depths.group),2,allele.summary.individuals)
colnames(summary.het.indiv)<-gsub(".AD$",".HET",colnames(summary.het.indiv))


if( sum(paste(the.samples,"DP",sep=".") %in% colnames(indels))==length(the.samples)){
  allele.depths<-as.matrix(indels[,paste(the.samples,"DP",sep=".")])
   }else{
  allele.depths<-apply(as.matrix(allele.depths.group),2,allele.DP.from.AN)
  colnames(allele.depths)<-gsub(".AD",".DP",colnames(allele.depths))
}


#allele.depths<-as.matrix(indels[,paste(the.samples,"DP",sep=".")])

#allele.depths[5030,]
ncol<-length(the.samples)
nrow<-dim(indels)[1]
allele.depths<-as.numeric(allele.depths)
allele.depths[is.na(allele.depths)]<-0

if(is.null(dim(allele.depths))){
dim(allele.depths)<-c(nrow,ncol)
## allele.depths<-t(allele.depths)
colnames(allele.depths)<-paste(the.samples,"DP",sep=".")
}


genotypes<-as.matrix(indels[,paste(the.samples,"GT",sep=".")])
if(sum(colnames(summary.het.indiv)!=paste(the.samples,"HET",sep="."))!=0){summary.het.indiv<-summary.het.indiv[,paste(the.samples,"HET",sep=".")]}
if(sum(colnames(allele.depths)!=paste(the.samples,"DP",sep="."))!=0){summary.het.indiv<-summary.het.indiv[,paste(the.samples,"HET",sep=".")]}

## tapply(genotypes[,1],genotypes[,1],length)
passing.genos<- (genotypes=="1/1" & (( summary.het.indiv > het.high.alt) | allele.depths < het.read.thresh) &  allele.depths > depth.homo.low &  grepl("^snp",indels[,"TYPE"])) |
      (genotypes=="0/1" & ((summary.het.indiv >= het.low & summary.het.indiv <= het.high) | allele.depths < het.read.thresh)  &  allele.depths > depth.het.low &  grepl("^snp",indels[,"TYPE"])) |
      (genotypes=="0/0" & ((summary.het.indiv < het.low.ref) | allele.depths < het.read.thresh ) &  allele.depths > depth.homo.low  &  grepl("^snp",indels[,"TYPE"]))  |
      (genotypes=="1/1" &  allele.depths > depth.homo.low &  grepl("^indel",indels[,"TYPE"])) |
      (genotypes=="0/1" &  allele.depths > depth.het.low &  grepl("^indel",indels[,"TYPE"])) |
      (genotypes=="0/0" &  allele.depths > depth.homo.low &  grepl("^indel",indels[,"TYPE"]))


## test<- (genotypes=="0/0")& (summary.het.indiv < het.low.ref)

## &  allele.depths > depth.homo.low  &  grepl("^snp",indels[,"TYPE"]))
## apply(passing.genos,1,sum)

## allele.depths >= depth.homo.low

genotypes[!passing.genos]<-NA  ### exclude failing genotypes from the calculation


return(genotypes)

}



full.genotype.summary<-function(indels,the.samples,prefix="",suffix=""){
########################### ALL ###############
## just makes a summary , no filters so is faster


### the.samples<-all.sample.labels
#the.samples
#if(prefix!=""){prefix<-paste(prefix,".",sep="")}

if(is.null(rownames(indels))){core.ann<-c("chr","start","end","REF","ALT","TYPE") ;rownames(indels)<-build.key(indels,core.ann,add.chr.label=FALSE)}
key.indels<-rownames(indels)

allele.depths.group<-indels[,paste(the.samples,"AD",sep=".")]
summary.het.indiv<-apply(as.matrix(allele.depths.group),2,allele.summary.individuals)
colnames(summary.het.indiv)<-gsub(".AD$",".HET",colnames(summary.het.indiv))

if( sum(paste(the.samples,"DP",sep=".") %in% colnames(indels))==length(the.samples)){
  allele.depths<-as.matrix(indels[,paste(the.samples,"DP",sep=".")])
   }else{
  allele.depths<-apply(as.matrix(allele.depths.group),2,allele.DP.from.AN)
  colnames(allele.depths)<-gsub(".AD",".DP",colnames(allele.depths))
}


#allele.depths[1:10,1:10]
ncol<-length(the.samples)
nrow<-dim(indels)[1]
allele.depths<-as.numeric(allele.depths)
allele.depths[is.na(allele.depths)]<-0


if(is.null(dim(allele.depths))){
dim(allele.depths)<-c(nrow,ncol)
## allele.depths<-t(allele.depths)
colnames(allele.depths)<-paste(the.samples,"DP",sep=".")
}



genotypes<-as.matrix(indels[,paste(the.samples,"GT",sep=".")])
if(sum(colnames(summary.het.indiv)!=paste(the.samples,"HET",sep="."))!=0){summary.het.indiv<-summary.het.indiv[,paste(the.samples,"HET",sep=".")]}
if(sum(colnames(allele.depths)!=paste(the.samples,"DP",sep="."))!=0){summary.het.indiv<-summary.het.indiv[,paste(the.samples,"HET",sep=".")]}

## tapply(genotypes[,1],genotypes[,1],length)

summary.geno.group<-genotype.summary(genotypes)
colnames(summary.geno.group)<-paste(prefix,c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),suffix,sep="")
rownames(summary.geno.group)<-key.indels

#genotypes[genotypes!="0/1"]<-NA

allele.depths.group[genotypes!="0/1" | genotypes=="NA" | is.na(genotypes) ] <- NA
summary.depths.group<-apply(allele.depths.group,1,allele.summary.with.apply)
summary.depths.group<-t(summary.depths.group)
colnames(summary.depths.group)<-paste(prefix,c("Hetero.ALT.reads","Hetero.REF.reads","Hetero.Read.Balance"),suffix,sep="")
rownames(summary.depths.group)<-key.indels

colnames(summary.het.indiv)<-paste(prefix,colnames(summary.het.indiv),sep="")
rownames(summary.het.indiv)<-key.indels

 list(summary.geno.group=summary.geno.group,summary.depths.group=summary.depths.group,summary.het.indiv=summary.het.indiv)

}

#test<-cbind( indels[,colnames(indels) %in% global.labs], filter.table[,colnames(filter.table) %in% global.labs], filter.table.pholy[,colnames(filter.table.pholy) %in% global.labs] )


quality.filter<-function(indels,key.indels,quality.cut,quality.type,quality.dirn){
## sample.labels  ### this are in individual samples in ALL
## QUAL   QD HRun   SB 
##   50    5    5   -1 
## > global.quality.dirn
##      QUAL        QD      HRun        SB 
## "greater" "greater"    "less"    "less" 
## > global.quality.type
##      QUAL        QD      HRun        SB 
## "numeric" "numeric" "numeric" "numeric"

## quality.cut
## quality.type
## quality.dirn
###################################

quality.names<- names(quality.cut) ## qaulity.names is the column name in indel indels[1:10,quality.names[j]]
quality.measure<-names(quality.type) ## is the filter name :  X1000753W.Aff.DP_Low
quality.thresh<-matrix(data=TRUE,nrow=dim(indels)[1],ncol=length(quality.cut))
if(is.null(rownames(indels))){core.ann<-c("chr","start","end","REF","ALT","TYPE") ;rownames(indels)<-build.key(indels,core.ann,add.chr.label=FALSE)}
colnames(quality.thresh)<-quality.measure


 for(j in 1:length(quality.cut)){
  if(quality.type[j]=="factor"){

    if(quality.dirn[j]=="exact"){
                                quality.thresh[,quality.measure[j]] <- indels[,quality.names[j]]==quality.cut[j]
                                quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE
                                }

        if(quality.dirn[j]=="starts_with"){
                                quality.thresh[,quality.measure[j]] <- grepl(paste("^",quality.cut[j],sep=""),indels[,quality.names[j]]) 
                                quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE
                                }

        if(quality.dirn[j]=="ends_with"){
                                quality.thresh[,quality.measure[j]] <- grepl(paste(quality.cut[j],"$",sep=""),indels[,quality.names[j]]) 
                                quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE
                                }

                              }     ##character comparison NA gets FASLE automaticaly
  if(quality.type[j]=="numeric"){

    if(quality.dirn[j]=="greater"){
                         quality.thresh[,quality.measure[j]] <- as.integer( as.character(indels[,quality.names[j] ]) ) >= as.numeric(quality.cut[j]) ## get warning here of NAs intro by coersion
                         quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE  ## set NA to FALSE
                                  }
    if(quality.dirn[j]=="less"){
                         quality.thresh[,quality.measure[j]] <- as.integer( as.character(indels[,quality.names[j] ]) ) <  as.numeric(quality.cut[j]) ## get warning here of NAs intro by coersion
                         quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE  ## set NA to FALSE
                                  }
    if(quality.dirn[j]=="exact"){
                         quality.thresh[,quality.measure[j]] <- as.integer( as.character(indels[,quality.names[j] ]) ) ==  as.numeric(quality.cut[j]) ## get warning here of NAs intro by coersion
                         quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE  ## set NA to FALSE
                                  }

   
                              }
#  print(paste("passing filter",names(quality.type)[j],":",sum(quality.thresh[,quality.measure[j]] ),sep=" "))
}

## sum(test[[paste("All",names(quality.type)[j],sep="::")]] != quality.thresh[,names(quality.type)[j]])
quality.thresh

} # position.quality.filter subroutine end



## combine.position.indiv.filters<-function(sample.labels, global.quality,global.quality.names, global.quality.labs, global.quality.type, global.quality.dirn, indiv.quality, indiv.quality.names,   indiv.quality.labs,  indiv.quality.type, global.quality.dirn){

## names(indiv.quality.type)<-indiv.quality.labs                       
## names(indiv.quality)<-indiv.quality.names


## ###############################
## quality.labs<-expand.labels.to.samples(indiv.quality.labs,sample.labels)
## quality.names<-expand.labels.to.samples(indiv.quality.names,sample.labels)
## quality<-expand.quantity.to.samples(indiv.quality,sample.labels)
## quality.type<-expand.quantity.to.samples(indiv.quality.type,sample.labels)
## quality.dirn<-expand.quantity.to.samples(indiv.quality.dirn,sample.labels)

## ########################################################################
## ######################   Define global attributes ##################
## ######## Bulit ALL sample filter
## quality.names<-c(quality.names,global.quality.names)  # lables in the range list quality.names<-c("ID","DP","QUAL","QD")  #"QUAL < 30.0 || AB > 0.75 && DP > 25 || QD < 5.0 || HRun > 5 || SB > -0.10"
## quality.dirn<-c(quality.dirn,global.quality.dirn)

## quality.labs<-c(quality.labs,global.quality.labs)  # what to call them in the filter
## quality.cut<-c(quality,global.quality) ### remove the "." to keep the SNPs intact
## names(quality.cut)<-quality.labs
## names(quality.dirn)<-quality.labs
## quality.type<-c(quality.type,global.quality.type)
## names(quality.type)<-quality.names

## quality.names<- names(quality.cut) ## qaulity.names is the column name in indel indels[1:10,quality.names[j]]
## quality.measure<-names(quality.type)

## #### these two are in the same ORDER
## list(quality.cut=quality.cut,quality.type=quality.type,quality.dirn=quality.dirn)
## }  # unpack with # for(i in 1:length(start.data)){assign( names(start.data)[i],value=start.data[[i]])}



## quality.cut<-global.quality.cut
## quality.type<-global.quality.type
## quality.dirn<-global.quality.dirn   

position.quality.filter<-function(indels,quality.cut,quality.type,quality.dirn){

##### Missing or NA values are ALWAYS set to FALSE
## quality.cut
## quality.type
## quality.dirn
################################### only test attributes that are avaialable
avail<-names(quality.cut) %in% colnames(indels)
print(paste("Attributes not in indels file:",toString(names(quality.cut)[!avail]),sep=" "))
quality.cut<-quality.cut[avail]
quality.type<-quality.type[avail]
quality.dirn<-quality.dirn[avail]
  

quality.names<- names(quality.cut) ## qaulity.names is the column name in indel indels[1:10,quality.names[j]]
quality.measure<-names(quality.type) ## is the filter name :  X1000753W.Aff.DP_Low
## quality.thresh<-data.frame(key=key.indels,stringsAsFactors=FALSE)
quality.thresh<-matrix(data=TRUE,nrow=dim(indels)[1],ncol=length(quality.cut))
colnames(quality.thresh)<-quality.measure
#quality.thresh[1:5,]

 for(j in 1:length(quality.cut)){
  if(quality.type[j]=="factor"){

    if(quality.dirn[j]=="exact"){
                                quality.thresh[,quality.measure[j]] <- indels[,quality.names[j]]==quality.cut[j]
                                quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE  ## "J"=NA -> NA so NEEDED
                                }

        if(quality.dirn[j]=="starts_with"){
                                quality.thresh[,quality.measure[j]] <- grepl(paste("^",quality.cut[j],sep=""),indels[,quality.names[j]]) 
                                quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE  #  grepl("J",NA) -> FALSE so no really needed
                                }
    
        if(quality.dirn[j]=="not_starts_with"){
                                quality.thresh[,quality.measure[j]] <- grepl(paste("^",quality.cut[j],sep=""),indels[,quality.names[j]]) 
                                quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE  #  grepl("J",NA) -> FALSE so no really needed
                                }

        if(quality.dirn[j]=="ends_with"){
                                quality.thresh[,quality.measure[j]] <- grepl(paste(quality.cut[j],"$",sep=""),indels[,quality.names[j]]) 
                                quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE
                                }
    
        if(quality.dirn[j]=="not_ends_with"){
                                quality.thresh[,quality.measure[j]] <- !grepl(paste(quality.cut[j],"$",sep=""),indels[,quality.names[j]]) 
                                quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE
                                }

       if(quality.dirn[j]=="contains"){
                                quality.thresh[,quality.measure[j]] <- grepl(quality.cut[j],indels[,quality.names[j]]) 
                                quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE
                                }

      if(quality.dirn[j]=="not_contains"){
                                quality.thresh[,quality.measure[j]] <-  !grepl(quality.cut[j],indels[,quality.names[j]]) 
                                quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE # grepl("J",NA) -> FALSE so no really needed
                                quality.thresh[is.na(indels[,quality.names[j]]),quality.measure[j]] <- FALSE  #resquires additional check NA->FALSE ; FALSE->TRUE!!
                                }
    
            #  testing agaisnst more than one item c("CAT","DOG") ;-quality.cut[j] would need to be a list for this to work  properly
            ## if(quality.dirn[j]=="contains"){
            ##                     quality.thresh[,quality.measure[j]] <- indels[,quality.names[j]] %in% quality.cut[j] #  quality.cut[j] would need to be alist for this to work
            ##                     quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE #"J" %in% NA -> FALSE so not really needed
            ##                     }
    
            ##     if(quality.dirn[j]=="not_contains"){
            ##                     quality.thresh[,quality.measure[j]] <- !(indels[,quality.names[j]] %in% quality.cut[j])
            ##                     quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE #"J" %in% NA -> FALSE so not really needed
            ##                     quality.thresh[is.na(indels[,quality.names[j]]),quality.measure[j]] <- FALSE  #resquires additional check NA->FALSE ; FALSE->TRUE!!
            ##                     }

                              }     ## factor character comparison NA gets FASLE automaticaly
  if(quality.type[j]=="numeric"){

    if(quality.dirn[j]=="greater"){
                         quality.thresh[,quality.measure[j]] <- as.numeric( as.character(indels[,quality.names[j] ]) ) >= as.numeric(as.character(quality.cut[j])) ## get warning here of NAs intro by coersion
                         quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE  ## set NA to FALSE
                                  }
    if(quality.dirn[j]=="less"){
                         quality.thresh[,quality.measure[j]] <- as.numeric( as.character(indels[,quality.names[j] ]) ) <  as.numeric(as.character(quality.cut[j])) ## get warning here of NAs intro by coersion
                         quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE  ## set NA to FALSE
                                  }
    if(quality.dirn[j]=="exact"){
                         quality.thresh[,quality.measure[j]] <- as.numeric( as.character(indels[,quality.names[j] ]) ) ==  as.numeric(as.character(quality.cut[j])) ## get warning here of NAs intro by coersion
                         quality.thresh[is.na(quality.thresh[,quality.measure[j]]),quality.measure[j]] <- FALSE  ## set NA to FALSE
                                  }

   
                              }
#  print(paste("passing filter",names(quality.type)[j],":",sum(quality.thresh[,quality.measure[j]] ),sep=" "))
}

## sum(test[[paste("All",names(quality.type)[j],sep="::")]] != quality.thresh[,names(quality.type)[j]])
quality.thresh

} # position.quality.filter subroutine end

#quality.thresh[1:5,]




## indels.bit<-indels[1:4442,]

correct.multi.alleles<-function(indels.bit,samples.order.in.ALL){ ## add a flattern genotype as well to expanded indels

 #####################################################################uses just indels from this point ###################################
######## Check not multi ALT alleles if there are then I need to flatten that line- not done here see build.annotation.files.r
## cols.wanted<-colnames(indels.bit) %in% core.ann | grepl(".GT$",colnames(indels.bit))
## print(class(indels.bit))
## print(dim(indels.bit))
##    indels.bit[6607:6614,cols.wanted]
 add.flattern<-TRUE
           
 print(dim(indels.bit))
the.cols<-colnames(indels.bit)
has.comma.ref<-grepl(",",indels.bit[,"REF"])
if(sum(has.comma.ref)>0){print("ERROR comma in REF allele")}

  ###### expand indels.bit to make room for new data

  alt.list<-strsplit(indels.bit[,"ALT"],split=",")
  number.of.alleles<-unlist(lapply(alt.list,length))
#  number.of.alleles<- ((number.of.alleles+1)*(number.of.alleles+2)/2)- (number.of.alleles+1)  # extra +1 to count REf allele
  if(add.flattern){number.of.alleles[number.of.alleles>1]<-number.of.alleles[number.of.alleles>1]+1} ### so can have the collapsed one too

  flat.index<-rep(1:length(number.of.alleles),times=number.of.alleles)
  indels.bit.ori.size<-dim(indels.bit)[1]
  indels.bit<-indels.bit[flat.index,]

  if(is.null(dim(indels.bit))){ # rare case id indels.bit is only one element long so above line makes a col
  the.cols<-names(indels.bit)
  indels.bit<-matrix(data=indels.bit,nrow=1,ncol=length(indels.bit))
  colnames(indels.bit)<-the.cols
}
################### fix the genotypes

  
if(dim(indels.bit)[1]> indels.bit.ori.size){
null.genotype<-{} # posns in indels.bit where there are no genotypes  
## to.fix<-  c("AD","GT","AD")
fix.posns<-tapply(flat.index,flat.index,length)
fix.posns<-fix.posns[fix.posns>1]
location.in.data<-match(names(fix.posns),flat.index)
names(fix.posns)<-location.in.data
# fix.posns
##  ## Names(fix.posn) give the index location in data to strat and

## for(ifix in 1:length(to.fix)){ ##loop over items to fix then put into cols.wanted
##   a.string<-to.fix[ifix]
##   if(a.string %in% format.types){ cols.wanted<-paste(samples.order.in.ALL,a.string,sep=".")}
##   if(a.string %in% info.types ){cols.wanted<-a.string}
##   cols.wanted indels.bit[10847:10859,1:6]

# print(paste("fixing",length(fix.posns),"poly morphic sites",sep=" "))
## ipos<-1  grep("4799",names(fix.posns)) grep("rs141819080",indels.bit[,"ID"])
    for(ipos in 1:length(fix.posns)){
       print(ipos)
      a.posn<-as.integer(names(fix.posns[ipos]))
      num<-fix.posns[ipos]

      alleles <- unlist(strsplit(c(paste(indels.bit[a.posn,"REF"],indels.bit[a.posn,"ALT"], sep=",")),','))
      names(alleles)<-0:(length(alleles)-1)
      ## alleles
      ACs <- unlist(strsplit(indels.bit[a.posn,"AC"],','))
      names(ACs)<-1:(length(alleles)-1) # one AC for each ALT allele
      AFs <- unlist(strsplit(indels.bit[a.posn,"AF"],','))
      names(AFs)<-1:(length(alleles)-1) # one AF for each ALT allele

 

        ############### do AD for all samples are this is the sample for each posn, just replicate once
        ############### sort out the ref and alt alleles 1/1 could be from 0,1 OR 1,2
       ###############  use the PL to decide and test
     
         cols.wanted<-paste(samples.order.in.ALL,"GT",sep=".")
         the.genotypes<-indels.bit[a.posn,cols.wanted]
         the.alleles<-the.genotypes
         the.alleles[the.alleles=="NA" | is.na(the.alleles) ]<-"0/0" #just so next step does not overlop the.genotypes==NA are filled with 0,0,0 AD
         the.alleles<-strsplit(the.alleles,split="/")

         ref.allele<-as.integer(unlist(lapply(the.alleles,function(x) x[1]))) # +1 # COLUMN POSN IN AD.LIST
         ## ref.allele.ori<-ref.allele
         alt.allele<-as.integer(unlist(lapply(the.alleles,function(x) x[2]))) # +1 # COLUMN POSN IN AD.LIST

      ##################
      ## ref.allele
      ## ref.allele.ori
      ##  alt.allele
      ## a.posn
      ## indels.bit[a.posn,cols.wanted]
    

      ############################### define the PL array so can use to get most likely alleles for ties
      
       PL.combo<-matrix(data=NA,nrow=length(alleles),ncol=length(alleles))
      ############ set up the PL structure
      for(iref in 1:length(alleles)){
        for(ialt in iref:length(alleles)){
          i.ref<-as.character(iref) # cause sometimes want to refer to a position by name
          i.alt<-as.character(ialt) # cause sometimes want to refer to a position by name
          PL.combo[iref,ialt]<-paste(iref-1,ialt-1,sep=":") 
       }}
      PL.combo<-as.character(PL.combo)
      PL.combo<-PL.combo[!is.na(PL.combo)]
          ## PL.combo ### the v4 samtools 

         cols.wanted<-paste(samples.order.in.ALL,"PL",sep=".")
         the.PL<- indels.bit[a.posn,cols.wanted]
         the.PL[the.PL=="NA" | is.na(the.PL)]<-paste(rep("3000",times=length(PL.combo)),collapse=",")
         pl.list<-strsplit(the.PL,split=",")
         pl.list<-pl.list[cols.wanted]
         pl.list.length<-length(pl.list)
         pl.list<-as.integer(unlist(pl.list))
         dim(pl.list)<-c(length(PL.combo), pl.list.length)
         pl.list<-t(pl.list)
         colnames(pl.list)<-PL.combo ## p.list and array of sample by allele combos 



      
       ######### fix location where have   1/1 or 2/2 or 3/3 etc choose two most probable REF allele 0 if a tie..
       ### no loger need this the ref allele is always 0 now
       ref.allele<-rep(0,times=length(ref.allele))
      
         ## ties<-which(ref.allele==alt.allele)
         ## pl.ties<-pl.list[ties,] ## P
         ## if(length(ties)>0){  ## just in case is no ties which causes an error
         ##   # ir<-1
         ## for(ir in 1:length(ties)){
         ##   if(is.null(dim(pl.ties))){pl.ties<-as.matrix(pl.ties);pl.ties<-t(pl.ties)} ## in case inly one ties and get a vector!
         ##   best.pl<-sort(pl.ties[ir,]) ## this sort works as 0:1 reported before 1:2 is they are tied (see PL.combo order)
         ##   the.alleles<-{}
         ##   while(length(the.alleles)<2){
         ##     the.alleles<-sort(as.integer(unique(unlist(strsplit(paste(names(best.pl[1:2]),collapse=":"),split=":")))))
         ##   }
         ##   ref.allele[ties[ir]]<-the.alleles[1]
         ##   alt.allele[ties[ir]]<-the.alleles[2]
        
         ## }}
      ###################################################
       ## ref.allele ### these contain the best choese of alleles fo the genotypes listed
       ## alt.allele  ### these contain the best choese of alleles fo the genotypes listed 
      ############### do AD for all samples are this is the sample for each posn, just replicate once

      #####get The AD 
         cols.wanted<-paste(samples.order.in.ALL,"AD",sep=".")
         the.AN<- indels.bit[a.posn,cols.wanted]
         the.AN[the.AN=="NA" | is.na(the.AN) | the.AN=="."]<-paste(rep("0",times=length(alleles)),collapse=",")
         ad.list<-strsplit(the.AN,split=",")
         ad.list<-ad.list[cols.wanted]
         ad.list.length<-length(ad.list)
         ad.list<-as.integer(unlist(ad.list))
         dim(ad.list)<-c(length(alleles), ad.list.length)
         ad.list<-t(ad.list)

    
      ref.counts<-rep(NA,times=dim(ad.list)[1])
      alt.counts<-rep(NA,times=dim(ad.list)[1])
         for(ir in 1:dim(ad.list)[1]){
           ref.counts[ir]<-ad.list[ir,(ref.allele[ir]+1)] ## ad list is ordered in assending allele
           alt.counts[ir]<-ad.list[ir,(alt.allele[ir]+1)]
         }
         the.AN<-paste(ref.counts,alt.counts,sep=",")
         names(the.AN)<-cols.wanted
         
      #################################################################
        
      #####get The PL
      cols.wanted<-paste(samples.order.in.ALL,"PL",sep=".")
      the.PL<-rep(NA,times=dim(pl.list)[1])
         for(ir in 1:dim(pl.list)[1]){
            pl.wanted<-c(paste(ref.allele[ir],ref.allele[ir],sep=":"),
                         paste(ref.allele[ir],alt.allele[ir],sep=":"),
                         paste(alt.allele[ir],alt.allele[ir],sep=":"))
           
           the.PL[ir]<-paste(pl.list[ir,pl.wanted],collapse=":")
         }
         names(the.PL)<-cols.wanted
         
     #################################################################

      ## the.genotypes
      ## ACs
      ## AFs
      ## the.AN
      ## the.PL
      ## ref.allele
      ## ref.allele.ori
      ## alt.allele
      ## ad.list[1:5,]
      ## pl.list[1:5,]
      
      ## print(alleles)



###########################################################
iplace<-1 #################### now fix all positions that start with a.posn iplace goes from 1:num

      if(add.flattern){
        iref<-0
        ialt<-1
      ## for(iref in 0:(length(alleles)-1)){
      ##   for(ialt in iref:(length(alleles)-1)){
      ##     if(iref==ialt){next}
         i.ref<-as.character(iref) # cause sometimes want to refer to a position by name
         i.alt<-as.character(ialt) # cause sometimes want to refer to a position by name
      ##    print(paste(i.ref,i.alt,sep=":"))
      ##   }}
     
         a.place<-a.posn+iplace-1 # -1 because a.posn is the one I want to change

          ## print(a.place)
          ## indels.bit[a.place,cols.wanted]
          ## test[a.place,cols.wanted]
          ## indels.bit[a.place,1:10]

          
          indels.bit[a.place,"REF"]<-alleles[i.ref]
          indels.bit[a.place,"ALT"]<-alleles[i.alt]
          indels.bit[a.place,"AC"]<-sum(as.numeric(ACs),na.rm=TRUE) # 
          indels.bit[a.place,"AF"]<-sum(as.numeric(AFs),na.rm=TRUE) # AFs[i.alt]
          
          if("RPA" %in% the.cols){indels.bit[a.place,"RPA"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("HRun" %in% the.cols){indels.bit[a.place,"HRun"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("SB" %in% the.cols){indels.bit[a.place,"SB"]<--1} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("FS" %in% the.cols){indels.bit[a.place,"FS"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          
          indels.bit[a.place,"TYPE"]<-paste(indels.bit[a.place,"TYPE"],indels.bit[a.place,"POS"],"flat",sep=":")

         ########################################  FIX GT ################################

          new.gt<-the.genotypes
          target.GT<-!is.na(new.gt) | new.gt!="NA"

          
          new.gt[!target.GT]<-"NA"
          if(sum(target.GT)==0){null.genotype<-c(null.genotype,a.place)}
          

          ######### here convert all possible alternative alleles to one ( 2 stage so can handle 1/10 11/15 etc:
          to.one<-paste("^[",paste(names(alleles)[!(names(alleles) %in% c(i.ref))],collapse=" "),"]/",sep="") ## everying to aly except  0
          new.gt<-gsub(to.one,"1/",new.gt)  #convect all left alleles to zero
          to.one<-paste("/[",paste(names(alleles)[!(names(alleles) %in% c(i.ref))],collapse=" "),"]$",sep="") ## everying to aly except  0
          new.gt<-gsub(to.one,"/1",new.gt)  #convect all right alleles to zero
          
          new.gt<-gsub(i.alt,"1",new.gt)
          new.gt[new.gt=="1/0"]<-"0/1"
        
        
          cols.wanted<-paste(samples.order.in.ALL,"GT",sep=".")
          indels.bit[a.place,cols.wanted]<-new.gt

        ##################################################################################
          
        ########################################  FIX AD ################################

        cols.wanted<-paste(samples.order.in.ALL,"AD",sep=".")
        new.AD<-the.AN
        new.AD[!target.GT]<-"NA"
        indels.bit[a.place,cols.wanted]<-new.AD
 
        ##################################################################################

        ########################################  FIX PL ################################

        cols.wanted<-paste(samples.order.in.ALL,"PL",sep=".")
        new.PL<-the.PL
        new.PL[!target.GT]<-"NA"
        indels.bit[a.place,cols.wanted]<-new.PL
 
        ##################################################################################
        ################################# step one out to get the others

      
 
      iplace<-2 #################### now fix all positions that start with a.posn iplace goes from 1:num

      } # end add flattern if added flatten iplace is 2 otherwise is 1

      ## iref<-0
      ## ialt<-2
      for(iref in 0:0){  ## iref==0 always as now always have the reference alleles
        for(ialt in iref:(length(alleles)-1)){
          if(iref==ialt){next}
          i.ref<-as.character(iref) # cause sometimes want to refer to a position by name
          i.alt<-as.character(ialt) # cause sometimes want to refer to a position by name
         ## print(paste(i.ref,i.alt,sep=":"))
         ##  }}
     
         a.place<-a.posn+iplace-1 # -1 because a.posn is the one I want to change

           ## print(a.place)
          ## indels.bit[a.place,cols.wanted]
          ## test[a.place,cols.wanted]
          ## indels.bit[a.place,1:35]

          
          indels.bit[a.place,"REF"]<-alleles[i.ref]
          indels.bit[a.place,"ALT"]<-alleles[i.alt]
          indels.bit[a.place,"AC"]<-ACs[i.alt]
          indels.bit[a.place,"AF"]<-AFs[i.alt]
          
          if("RPA" %in% the.cols){indels.bit[a.place,"RPA"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("HRun" %in% the.cols){indels.bit[a.place,"HRun"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("SB" %in% the.cols){indels.bit[a.place,"SB"]<--1} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("FS" %in% the.cols){indels.bit[a.place,"FS"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          
          indels.bit[a.place,"TYPE"]<-paste(indels.bit[a.place,"TYPE"],indels.bit[a.place,"POS"],sep=":")

         ########################################  FIX GT ################################

          new.gt<-the.genotypes
          missing.genotype<-is.na(new.gt) | new.gt=="NA"
          target.GT<- !missing.genotype
          
          ## new.gt[!target.GT & !missing.genotype]<-"0/0"
          if(sum(target.GT)==0){null.genotype<-c(null.genotype,a.place)}
          

          ## new.gt<-gsub(i.alt,"1",new.gt) ## set the 0/2 to 0/1 --- fails in case where have 2/3 for i.alt=2
          
          ## new.gt<-gsub(i.ref,"0",new.gt) ## set the other alleles to reference now refernce alwasy zero
          ## Need a 2 step process for alleles like 0/10 11/11 etc 
          to.zero<-paste("^[",paste(names(alleles)[!(names(alleles) %in% c(i.alt))],collapse=" "),"]/",sep="") ## everying to aly except  0
          new.gt<-gsub(to.zero,"0/",new.gt)  #convect all left alleles to zero
          to.zero<-paste("/[",paste(names(alleles)[!(names(alleles) %in% c(i.alt))],collapse=" "),"]$",sep="") ## everying to aly except  0
          new.gt<-gsub(to.zero,"/0",new.gt)  #convect all right alleles to zero
          
          new.gt<-gsub(i.alt,"1",new.gt)
          new.gt[new.gt=="1/0"]<-"0/1"

          
          cols.wanted<-paste(samples.order.in.ALL,"GT",sep=".")
          indels.bit[a.place,cols.wanted]<-new.gt

        ##################################################################################
          
        ########################################  FIX AD ################################

        cols.wanted<-paste(samples.order.in.ALL,"AD",sep=".")
        new.AD<-the.AN
        new.AD[!target.GT]<-"NA"
        indels.bit[a.place,cols.wanted]<-new.AD
 
        ##################################################################################

        ########################################  FIX PL ################################

        cols.wanted<-paste(samples.order.in.ALL,"PL",sep=".")
        new.PL<-the.PL
        new.PL[!target.GT]<-"NA"
        indels.bit[a.place,cols.wanted]<-new.PL
 
        ##################################################################################
         
   
         iplace<-iplace+1
          }} #iref ialt iplace fill in same genotypes for location
      
    } #ipos new location
  
#} #ifix
    if(length(null.genotype)>1){  indels.bit<-indels.bit[-1*null.genotype,]}
    } # indels.bit.ori.size has to do work and fix alleles



 #####################################################################uses just indels.bit from this point ###################################

indels.bit



}  ## end function correct.muli.alleles





correct.multi.alleles.ref.change<-function(indels.bit,samples.order.in.ALL){ ## add a flattern genotype as well to expanded indels

 #####################################################################uses just indels from this point ###################################
######## Check not multi ALT alleles if there are then I need to flatten that line- not done here see build.annotation.files.r
## cols.wanted<-colnames(indels.bit) %in% core.ann | grepl(".GT$",colnames(indels.bit))
## indels.bit[
##    indels.bit[6607:6614,cols.wanted]

           
 ## print(dim(indels.bit))
the.cols<-colnames(indels.bit)
has.comma.ref<-grepl(",",indels.bit[,"REF"])
if(sum(has.comma.ref)>0){print("ERROR comma in REF allele")}

  ###### expand indels.bit to make room for new data

  alt.list<-strsplit(indels.bit[,"ALT"],split=",")
  number.of.alleles<-unlist(lapply(alt.list,length))
  number.of.alleles<- ((number.of.alleles+1)*(number.of.alleles+2)/2)- (number.of.alleles+1)  # extra +1 to count REf allele
  number.of.alleles[number.of.alleles>1]<-number.of.alleles[number.of.alleles>1]+1 ### so can have the collapsed one too

  flat.index<-rep(1:length(number.of.alleles),times=number.of.alleles)
  indels.bit.ori.size<-dim(indels.bit)[1]
  indels.bit<-indels.bit[flat.index,]

  if(is.null(dim(indels.bit))){ # rare case id indels.bit is only one element long so above line makes a col
  the.cols<-names(indels.bit)
  indels.bit<-matrix(data=indels.bit,nrow=1,ncol=length(indels.bit))
  colnames(indels.bit)<-the.cols
}
################### fix the genotypes

  
if(dim(indels.bit)[1]> indels.bit.ori.size){
null.genotype<-{} # posns in indels.bit where there are no genotypes  
## to.fix<-  c("AD","GT","AD")
fix.posns<-tapply(flat.index,flat.index,length)
fix.posns<-fix.posns[fix.posns>1]
location.in.data<-match(names(fix.posns),flat.index)
names(fix.posns)<-location.in.data
# fix.posns
##  ## Names(fix.posn) give the index location in data to strat and

## for(ifix in 1:length(to.fix)){ ##loop over items to fix then put into cols.wanted
##   a.string<-to.fix[ifix]
##   if(a.string %in% format.types){ cols.wanted<-paste(samples.order.in.ALL,a.string,sep=".")}
##   if(a.string %in% info.types ){cols.wanted<-a.string}
##   cols.wanted indels.bit[10847:10859,1:6]

# print(paste("fixing",length(fix.posns),"poly morphic sites",sep=" "))
## ipos<-892  grep("8126",names(fix.posns)) grep("rs201485735",indels.bit[,"ID"])
    for(ipos in 1:length(fix.posns)){
      ## print(ipos)
      a.posn<-as.integer(names(fix.posns[ipos]))
      num<-fix.posns[ipos]

      alleles <- unlist(strsplit(c(paste(indels.bit[a.posn,"REF"],indels.bit[a.posn,"ALT"], sep=",")),','))
      names(alleles)<-0:(length(alleles)-1)
      ## alleles
      ACs <- unlist(strsplit(indels.bit[a.posn,"AC"],','))
      names(ACs)<-1:(length(alleles)-1) # one AC for each ALT allele
      AFs <- unlist(strsplit(indels.bit[a.posn,"AF"],','))
      names(AFs)<-1:(length(alleles)-1) # one AF for each ALT allele

 

        ############### do AD for all samples are this is the sample for each posn, just replicate once
        ############### sort out the ref and alt alleles 1/1 could be from 0,1 OR 1,2
       ###############  use the PL to decide and test
     
         cols.wanted<-paste(samples.order.in.ALL,"GT",sep=".")
         the.genotypes<-indels.bit[a.posn,cols.wanted]
         the.alleles<-the.genotypes
         the.alleles[the.alleles=="NA" | is.na(alleles) ]<-"0/0" #just so next step does not overlop the.genotypes==NA are filled with 0,0,0 AD
         the.alleles<-strsplit(the.alleles,split="/")

         ref.allele<-as.integer(unlist(lapply(the.alleles,function(x) x[1]))) # +1 # COLUMN POSN IN AD.LIST
         alt.allele<-as.integer(unlist(lapply(the.alleles,function(x) x[2]))) # +1 # COLUMN POSN IN AD.LIST

      ##################
      ## ref.allele
      ##  alt.allele
      ## a.posn
      ## indels.bit[a.posn,cols.wanted]
      ## test[a.posn,1:10]

      ############################### define the PL array so can use to get most likely alleles for ties
      
       PL.combo<-matrix(data=NA,nrow=length(alleles),ncol=length(alleles))
      ############ set up the PL structure
      for(iref in 1:length(alleles)){
        for(ialt in iref:length(alleles)){
          i.ref<-as.character(iref) # cause sometimes want to refer to a position by name
          i.alt<-as.character(ialt) # cause sometimes want to refer to a position by name
          PL.combo[iref,ialt]<-paste(iref-1,ialt-1,sep=":") 
       }}
      PL.combo<-as.character(PL.combo)
      PL.combo<-PL.combo[!is.na(PL.combo)]
          ## PL.combo ### the v4 samtools 

         cols.wanted<-paste(samples.order.in.ALL,"PL",sep=".")
         the.PL<- indels.bit[a.posn,cols.wanted]
         the.PL[the.PL=="NA" | is.na(the.PL)]<-paste(rep("3000",times=length(PL.combo)),collapse=",")
         pl.list<-strsplit(the.PL,split=",")
         pl.list<-pl.list[cols.wanted]
         pl.list.length<-length(pl.list)
         pl.list<-as.integer(unlist(pl.list))
         dim(pl.list)<-c(length(PL.combo), pl.list.length)
         pl.list<-t(pl.list)
         colnames(pl.list)<-PL.combo ## p.list and array of sample by allele combos 



      
       ######### fix location where have   1/1 or 2/2 or 3/3 etc choose two most probable REF allele 0 if a tie..
         ties<-which(ref.allele==alt.allele)
         pl.ties<-pl.list[ties,] ## P
         if(length(ties)>0){  ## just in case is no ties which causes an error
           # ir<-1
         for(ir in 1:length(ties)){
           if(is.null(dim(pl.ties))){pl.ties<-as.matrix(pl.ties);pl.ties<-t(pl.ties)} ## in case inly one ties and get a vector!
           best.pl<-sort(pl.ties[ir,]) ## this sort works as 0:1 reported before 1:2 is they are tied (see PL.combo order)
           the.alleles<-{}
           while(length(the.alleles)<2){
             the.alleles<-sort(as.integer(unique(unlist(strsplit(paste(names(best.pl[1:2]),collapse=":"),split=":")))))
           }
           ref.allele[ties[ir]]<-the.alleles[1]
           alt.allele[ties[ir]]<-the.alleles[2]
        
         }}
      ###################################################
       ## ref.allele ### these contain the best choese of alleles fo the genotypes listed
       ## alt.allele  ### these contain the best choese of alleles fo the genotypes listed 
      ############### do AD for all samples are this is the sample for each posn, just replicate once

      #####get The AD 
         cols.wanted<-paste(samples.order.in.ALL,"AD",sep=".")
         the.AN<- indels.bit[a.posn,cols.wanted]
         the.AN[the.AN=="NA" | is.na(the.AN)]<-paste(rep("0",times=length(alleles)),collapse=",")
         ad.list<-strsplit(the.AN,split=",")
         ad.list<-ad.list[cols.wanted]
         ad.list.length<-length(ad.list)
         ad.list<-as.integer(unlist(ad.list))
         dim(ad.list)<-c(length(alleles), ad.list.length)
         ad.list<-t(ad.list)

    
      ref.counts<-rep(NA,times=dim(ad.list)[1])
      alt.counts<-rep(NA,times=dim(ad.list)[1])
         for(ir in 1:dim(ad.list)[1]){
           ref.counts[ir]<-ad.list[ir,(ref.allele[ir]+1)] ## ad list is ordered in assending allele
           alt.counts[ir]<-ad.list[ir,(alt.allele[ir]+1)]
         }
         the.AN<-paste(ref.counts,alt.counts,sep=",")
         names(the.AN)<-cols.wanted
         
      #################################################################
        
      #####get The PL
      cols.wanted<-paste(samples.order.in.ALL,"PL",sep=".")
      the.PL<-rep(NA,times=dim(pl.list)[1])
         for(ir in 1:dim(pl.list)[1]){
            pl.wanted<-c(paste(ref.allele[ir],ref.allele[ir],sep=":"),
                         paste(ref.allele[ir],alt.allele[ir],sep=":"),
                         paste(alt.allele[ir],alt.allele[ir],sep=":"))
           
           the.PL[ir]<-paste(pl.list[ir,pl.wanted],collapse=":")
         }
         names(the.PL)<-cols.wanted
         
     #################################################################

      the.genotypes
      ACs
      AFs
      the.AN
      the.PL
      ref.allele
      alt.allele
      ad.list[1:5,]
      pl.list[1:5,]
      
      ## print(alleles)



###########################################################
iplace<-1 #################### now fix all positions that start with a.posn iplace goes from 1:num
iref<-0
ialt<-1
      ## for(iref in 0:(length(alleles)-1)){
      ##   for(ialt in iref:(length(alleles)-1)){
      ##     if(iref==ialt){next}
         i.ref<-as.character(iref) # cause sometimes want to refer to a position by name
         i.alt<-as.character(ialt) # cause sometimes want to refer to a position by name
      ##    print(paste(i.ref,i.alt,sep=":"))
      ##   }}
     
         a.place<-a.posn+iplace-1 # -1 because a.posn is the one I want to change

          ## print(a.place)
          ## indels.bit[a.place,cols.wanted]
          ## test[a.place,cols.wanted]
          ## indels.bit[a.place,1:35]

          
          indels.bit[a.place,"REF"]<-alleles[i.ref]
          indels.bit[a.place,"ALT"]<-alleles[i.alt]
          indels.bit[a.place,"AC"]<-sum(as.numeric(ACs),na.rm=TRUE) # 
          indels.bit[a.place,"AF"]<-sum(as.numeric(AFs),na.rm=TRUE) # AFs[i.alt]
          
          if("RPA" %in% the.cols){indels.bit[a.place,"RPA"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("HRun" %in% the.cols){indels.bit[a.place,"HRun"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("SB" %in% the.cols){indels.bit[a.place,"SB"]<--1} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("FS" %in% the.cols){indels.bit[a.place,"FS"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          
          indels.bit[a.place,"TYPE"]<-paste(indels.bit[a.place,"TYPE"],indels.bit[a.place,"POS"],"flat",sep=":")

         ########################################  FIX GT ################################

          new.gt<-the.genotypes
          target.GT<-!is.na(new.gt) | new.gt!="NA"
          
          new.gt[!target.GT]<-"NA"
          if(sum(target.GT)==0){null.genotype<-c(null.genotype,a.place)}
          
          new.gt<-gsub(i.ref,"0",new.gt) ## set the other alleles to reference # ij<-2
          for(ij in 2:length(alleles)){
          new.gt<-gsub(names(alleles)[ij],"1",new.gt) ## set the 0/2 to 0/1
        }
          cols.wanted<-paste(samples.order.in.ALL,"GT",sep=".")
          indels.bit[a.place,cols.wanted]<-new.gt

        ##################################################################################
          
        ########################################  FIX AD ################################

        cols.wanted<-paste(samples.order.in.ALL,"AD",sep=".")
        new.AD<-the.AN
        new.AD[!target.GT]<-"NA"
        indels.bit[a.place,cols.wanted]<-new.AD
 
        ##################################################################################

        ########################################  FIX PL ################################

        cols.wanted<-paste(samples.order.in.ALL,"PL",sep=".")
        new.PL<-the.PL
        new.PL[!target.GT]<-"NA"
        indels.bit[a.place,cols.wanted]<-new.PL
 
        ##################################################################################
        ################################# step one out to get the others

      
      ## iref<-0
      ## ialt<-1
      iplace<-2 #################### now fix all positions that start with a.posn iplace goes from 1:num
      for(iref in 0:(length(alleles)-1)){
        for(ialt in iref:(length(alleles)-1)){
          if(iref==ialt){next}
          i.ref<-as.character(iref) # cause sometimes want to refer to a position by name
          i.alt<-as.character(ialt) # cause sometimes want to refer to a position by name
         ## print(paste(i.ref,i.alt,sep=":"))
         ## }}
     
         a.place<-a.posn+iplace-1 # -1 because a.posn is the one I want to change

          ## print(a.place)
          ## indels.bit[a.place,cols.wanted]
          ## test[a.place,cols.wanted]
          ## indels.bit[a.place,1:35]

          
          indels.bit[a.place,"REF"]<-alleles[i.ref]
          indels.bit[a.place,"ALT"]<-alleles[i.alt]
          indels.bit[a.place,"AC"]<-ACs[i.alt]
          indels.bit[a.place,"AF"]<-AFs[i.alt]
          
          if("RPA" %in% the.cols){indels.bit[a.place,"RPA"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("HRun" %in% the.cols){indels.bit[a.place,"HRun"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("SB" %in% the.cols){indels.bit[a.place,"SB"]<--1} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("FS" %in% the.cols){indels.bit[a.place,"FS"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          
          indels.bit[a.place,"TYPE"]<-paste(indels.bit[a.place,"TYPE"],indels.bit[a.place,"POS"],sep=":")

         ########################################  FIX GT ################################

          new.gt<-the.genotypes
          target.GT<- iref==ref.allele & ialt==alt.allele
          
          new.gt[!target.GT]<-"NA"
          if(sum(target.GT)==0){null.genotype<-c(null.genotype,a.place)}
          
          new.gt<-gsub(i.ref,"0",new.gt) ## set the other alleles to reference
          new.gt<-gsub(i.alt,"1",new.gt) ## set the 0/2 to 0/1
          cols.wanted<-paste(samples.order.in.ALL,"GT",sep=".")
          indels.bit[a.place,cols.wanted]<-new.gt

        ##################################################################################
          
        ########################################  FIX AD ################################

        cols.wanted<-paste(samples.order.in.ALL,"AD",sep=".")
        new.AD<-the.AN
        new.AD[!target.GT]<-"NA"
        indels.bit[a.place,cols.wanted]<-new.AD
 
        ##################################################################################

        ########################################  FIX PL ################################

        cols.wanted<-paste(samples.order.in.ALL,"PL",sep=".")
        new.PL<-the.PL
        new.PL[!target.GT]<-"NA"
        indels.bit[a.place,cols.wanted]<-new.PL
 
        ##################################################################################
         
   
         iplace<-iplace+1
          }} #iref ialt iplace fill in same genotypes for location
      
    } #ipos new location
  
#} #ifix
    if(length(null.genotype)>1){  indels.bit<-indels.bit[-1*null.genotype,]}
    } # indels.bit.ori.size has to do work and fix alleles



 #####################################################################uses just indels.bit from this point ###################################

indels.bit



}  ## end function correct.muli.alleles.ref.change



correct.multi.alleles.ref.change.no.flattern<-function(indels.bit,samples.order.in.ALL){ ##  one used to march 2013

 #####################################################################uses just indels from this point ###################################
######## Check not multi ALT alleles if there are then I need to flatten that line- not done here see build.annotation.files.r
## cols.wanted<-colnames(indels.bit) %in% core.ann | grepl(".GT$",colnames(indels.bit))
## indels.bit[
##    indels.bit[6607:6614,cols.wanted]

           
 ## print(dim(indels.bit))
the.cols<-colnames(indels.bit)
has.comma.ref<-grepl(",",indels.bit[,"REF"])
if(sum(has.comma.ref)>0){print("ERROR comma in REF allele")}

  ###### expand indels.bit to make room for new data

  alt.list<-strsplit(indels.bit[,"ALT"],split=",")
  number.of.alleles<-unlist(lapply(alt.list,length))
  number.of.alleles<- ((number.of.alleles+1)*(number.of.alleles+2)/2)- (number.of.alleles+1) # extra +1 to count REf allele

  flat.index<-rep(1:length(number.of.alleles),times=number.of.alleles)
  indels.bit.ori.size<-dim(indels.bit)[1]
  indels.bit<-indels.bit[flat.index,]

  if(is.null(dim(indels.bit))){ # rare case id indels.bit is only one element long so above line makes a col
  the.cols<-names(indels.bit)
  indels.bit<-matrix(data=indels.bit,nrow=1,ncol=length(indels.bit))
  colnames(indels.bit)<-the.cols
}
################### fix the genotypes

  
if(dim(indels.bit)[1]> indels.bit.ori.size){
null.genotype<-{} # posns in indels.bit where there are no genotypes  
## to.fix<-  c("AD","GT","AD")
fix.posns<-tapply(flat.index,flat.index,length)
fix.posns<-fix.posns[fix.posns>1]
location.in.data<-match(names(fix.posns),flat.index)
names(fix.posns)<-location.in.data
# fix.posns
##  ## Names(fix.posn) give the index location in data to strat and

## for(ifix in 1:length(to.fix)){ ##loop over items to fix then put into cols.wanted
##   a.string<-to.fix[ifix]
##   if(a.string %in% format.types){ cols.wanted<-paste(samples.order.in.ALL,a.string,sep=".")}
##   if(a.string %in% info.types ){cols.wanted<-a.string}
##   cols.wanted indels.bit[7235:7244,1:6]

# print(paste("fixing",length(fix.posns),"poly morphic sites",sep=" "))
## ipos<-892  grep("7235",names(fix.posns))
    for(ipos in 1:length(fix.posns)){
      ## print(ipos)
      a.posn<-as.integer(names(fix.posns[ipos]))
      num<-fix.posns[ipos]

      alleles <- unlist(strsplit(c(paste(indels.bit[a.posn,"REF"],indels.bit[a.posn,"ALT"], sep=",")),','))
      names(alleles)<-0:(length(alleles)-1)
      ## alleles
      ACs <- unlist(strsplit(indels.bit[a.posn,"AC"],','))
      names(ACs)<-1:(length(alleles)-1) # one AC for each ALT allele
      AFs <- unlist(strsplit(indels.bit[a.posn,"AF"],','))
      names(AFs)<-1:(length(alleles)-1) # one AF for each ALT allele

 

        ############### do AD for all samples are this is the sample for each posn, just replicate once
        ############### sort out the ref and alt alleles 1/1 could be from 0,1 OR 1,2
       ###############  use the PL to decide and test
     
         cols.wanted<-paste(samples.order.in.ALL,"GT",sep=".")
         the.genotypes<-indels.bit[a.posn,cols.wanted]
         the.alleles<-the.genotypes
         the.alleles[the.alleles=="NA" | is.na(alleles) ]<-"0/0" #just so next step does not overlop the.genotypes==NA are filled with 0,0,0 AD
         the.alleles<-strsplit(the.alleles,split="/")

         ref.allele<-as.integer(unlist(lapply(the.alleles,function(x) x[1]))) # +1 # COLUMN POSN IN AD.LIST
         alt.allele<-as.integer(unlist(lapply(the.alleles,function(x) x[2]))) # +1 # COLUMN POSN IN AD.LIST

      ##################
      ## ref.allele
      ##  alt.allele
      ## a.posn
      ## indels.bit[a.posn,cols.wanted]
      ## test[a.posn,1:10]

      ############################### define the PL array so can use to get most likely alleles for ties
      
       PL.combo<-matrix(data=NA,nrow=length(alleles),ncol=length(alleles))
      ############ set up the PL structure
      for(iref in 1:length(alleles)){
        for(ialt in iref:length(alleles)){
          i.ref<-as.character(iref) # cause sometimes want to refer to a position by name
          i.alt<-as.character(ialt) # cause sometimes want to refer to a position by name
          PL.combo[iref,ialt]<-paste(iref-1,ialt-1,sep=":") 
       }}
      PL.combo<-as.character(PL.combo)
      PL.combo<-PL.combo[!is.na(PL.combo)]
          ## PL.combo ### the v4 samtools 

         cols.wanted<-paste(samples.order.in.ALL,"PL",sep=".")
         the.PL<- indels.bit[a.posn,cols.wanted]
         the.PL[the.PL=="NA" | is.na(the.PL)]<-paste(rep("3000",times=length(PL.combo)),collapse=",")
         pl.list<-strsplit(the.PL,split=",")
         pl.list<-pl.list[cols.wanted]
         pl.list.length<-length(pl.list)
         pl.list<-as.integer(unlist(pl.list))
         dim(pl.list)<-c(length(PL.combo), pl.list.length)
         pl.list<-t(pl.list)
         colnames(pl.list)<-PL.combo ## p.list and array of sample by allele combos 



      
       ######### fix location where have   1/1 or 2/2 or 3/3 etc choose two most probable REF allele 0 if a tie..
         ties<-which(ref.allele==alt.allele)
         pl.ties<-pl.list[ties,] ## P
         if(length(ties)>0){  ## just in case is no ties which causes an error
           # ir<-1
         for(ir in 1:length(ties)){
           if(is.null(dim(pl.ties))){pl.ties<-as.matrix(pl.ties);pl.ties<-t(pl.ties)} ## in case inly one ties and get a vector!
           best.pl<-sort(pl.ties[ir,]) ## this sort works as 0:1 reported before 1:2 is they are tied (see PL.combo order)
           the.alleles<-{}
           while(length(the.alleles)<2){
             the.alleles<-sort(as.integer(unique(unlist(strsplit(paste(names(best.pl[1:2]),collapse=":"),split=":")))))
           }
           ref.allele[ties[ir]]<-the.alleles[1]
           alt.allele[ties[ir]]<-the.alleles[2]
        
         }}
      ###################################################
       ## ref.allele ### these contain the best choese of alleles fo the genotypes listed
       ## alt.allele  ### these contain the best choese of alleles fo the genotypes listed 
      ############### do AD for all samples are this is the sample for each posn, just replicate once

      #####get The AD 
         cols.wanted<-paste(samples.order.in.ALL,"AD",sep=".")
         the.AN<- indels.bit[a.posn,cols.wanted]
         the.AN[the.AN=="NA" | is.na(the.AN)]<-paste(rep("0",times=length(alleles)),collapse=",")
         ad.list<-strsplit(the.AN,split=",")
         ad.list<-ad.list[cols.wanted]
         ad.list.length<-length(ad.list)
         ad.list<-as.integer(unlist(ad.list))
         dim(ad.list)<-c(length(alleles), ad.list.length)
         ad.list<-t(ad.list)

    
      ref.counts<-rep(NA,times=dim(ad.list)[1])
      alt.counts<-rep(NA,times=dim(ad.list)[1])
         for(ir in 1:dim(ad.list)[1]){
           ref.counts[ir]<-ad.list[ir,(ref.allele[ir]+1)] ## ad list is ordered in assending allele
           alt.counts[ir]<-ad.list[ir,(alt.allele[ir]+1)]
         }
         the.AN<-paste(ref.counts,alt.counts,sep=",")
         names(the.AN)<-cols.wanted
         
      #################################################################
        
      #####get The PL
      cols.wanted<-paste(samples.order.in.ALL,"PL",sep=".")
      the.PL<-rep(NA,times=dim(pl.list)[1])
         for(ir in 1:dim(pl.list)[1]){
            pl.wanted<-c(paste(ref.allele[ir],ref.allele[ir],sep=":"),
                         paste(ref.allele[ir],alt.allele[ir],sep=":"),
                         paste(alt.allele[ir],alt.allele[ir],sep=":"))
           
           the.PL[ir]<-paste(pl.list[ir,pl.wanted],collapse=":")
         }
         names(the.PL)<-cols.wanted
         
     #################################################################

      the.genotypes
      the.AN
      the.PL
      ref.allele
      alt.allele
      ad.list[1:5,]
      pl.list[1:5,]
      
      ## print(alleles)
      iplace<-1 #################### now fix all positions that start with a.posn iplace goes from 1:num
      for(iref in 0:(length(alleles)-1)){
        for(ialt in iref:(length(alleles)-1)){
          if(iref==ialt){next}
          i.ref<-as.character(iref) # cause sometimes want to refer to a position by name
          i.alt<-as.character(ialt) # cause sometimes want to refer to a position by name
          ## print(paste(i.ref,i.alt,sep=":"))
         ## }}
     
         a.place<-a.posn+iplace-1 # -1 because a.posn is the one I want to change

          ## print(a.place)
          ## indels.bit[a.place,cols.wanted]
          ## test[a.place,cols.wanted]
          ## indels.bit[a.place,1:35]

          
          indels.bit[a.place,"REF"]<-alleles[i.ref]
          indels.bit[a.place,"ALT"]<-alleles[i.alt]
          indels.bit[a.place,"AC"]<-ACs[i.alt]
          indels.bit[a.place,"AF"]<-AFs[i.alt]
          
          if("RPA" %in% the.cols){indels.bit[a.place,"RPA"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("HRun" %in% the.cols){indels.bit[a.place,"HRun"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("SB" %in% the.cols){indels.bit[a.place,"SB"]<--1} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("FS" %in% the.cols){indels.bit[a.place,"FS"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          
          indels.bit[a.place,"TYPE"]<-paste(indels.bit[a.place,"TYPE"],indels.bit[a.place,"POS"],sep=":")

         ########################################  FIX GT ################################

          new.gt<-the.genotypes
          target.GT<- iref==ref.allele & ialt==alt.allele
          
          new.gt[!target.GT]<-"NA"
          if(sum(target.GT)==0){null.genotype<-c(null.genotype,a.place)}
          
          new.gt<-gsub(i.ref,"0",new.gt) ## set the other alleles to reference
          new.gt<-gsub(i.alt,"1",new.gt) ## set the 0/2 to 0/1
          cols.wanted<-paste(samples.order.in.ALL,"GT",sep=".")
          indels.bit[a.place,cols.wanted]<-new.gt

        ##################################################################################
          
        ########################################  FIX AD ################################

        cols.wanted<-paste(samples.order.in.ALL,"AD",sep=".")
        new.AD<-the.AN
        new.AD[!target.GT]<-"NA"
        indels.bit[a.place,cols.wanted]<-new.AD
 
        ##################################################################################

        ########################################  FIX PL ################################

        cols.wanted<-paste(samples.order.in.ALL,"PL",sep=".")
        new.PL<-the.PL
        new.PL[!target.GT]<-"NA"
        indels.bit[a.place,cols.wanted]<-new.PL
 
        ##################################################################################
         
   
         iplace<-iplace+1
          }} #iref ialt iplace fill in same genotypes for location
      
    } #ipos new location
  
#} #ifix
    if(length(null.genotype)>1){  indels.bit<-indels.bit[-1*null.genotype,]}
    } # indels.bit.ori.size has to do work and fix alleles



 #####################################################################uses just indels.bit from this point ###################################

indels.bit



}  ## end function correct.muli.alleles


## indels.bit<-indels
correct.multi.alleles.flattern<-function(indels.bit,samples.order.in.ALL){

 #####################################################################uses just indels from this point ###################################
######## Check not multi ALT alleles if there are then I need to flatten that line- not done here see build.annotation.files.r
skip<-TRUE
  the.cols<-colnames(indels.bit)
#  print(dim(indels.bit))
has.comma.ref<-grepl(",",indels.bit[,"REF"])
if(sum(has.comma.ref)>0){print("ERROR comma in REF allele")}

  ###### expand indels.bit to make room for new data

  alt.list<-strsplit(indels.bit[,"ALT"],split=",")
  number.of.alleles<-unlist(lapply(alt.list,length))
  

if(!skip){
  number.of.alleles<- ((number.of.alleles+1)*(number.of.alleles+2)/2)- (number.of.alleles+1) # extra +1 to count REf allele
  flat.index<-rep(1:length(number.of.alleles),times=number.of.alleles)
  indels.bit.ori.size<-dim(indels.bit)[1]
  indels.bit<-indels.bit[flat.index,]
}

  if(is.null(dim(indels.bit))){ # rare case id indels.bit is only one element long so above line makes a col
  the.cols<-names(indels.bit)
  indels.bit<-matrix(data=indels.bit,nrow=1,ncol=length(indels.bit))
  colnames(indels.bit)<-the.cols
}
################### fix the genotypes

  
if( sum(number.of.alleles>1)>1 ){ ## check have somthing to do
  
null.genotype<-{} # posns in indels.bit where there are no genotypes  
## to.fix<-  c("AD","GT","AD")
if(!skip){
fix.posns<-tapply(flat.index,flat.index,length)
fix.posns<-fix.posns[fix.posns>1]
location.in.data<-match(names(fix.posns),flat.index)
names(fix.posns)<-location.in.data
}else{
fix.posns<-(1:dim(indels.bit)[1])[number.of.alleles>1]
names(fix.posns)<-fix.posns
}
  
#fix.posns
##  ## Names(fix.posn) give the index location in data to strat and

## for(ifix in 1:length(to.fix)){ ##loop over items to fix then put into cols.wanted
##   a.string<-to.fix[ifix]
##   if(a.string %in% format.types){ cols.wanted<-paste(samples.order.in.ALL,a.string,sep=".")}
##   if(a.string %in% info.types ){cols.wanted<-a.string}
##   cols.wanted

# print(paste("fixing",length(fix.posns),"poly morphic sites",sep=" "))
# ipos<-1

    for(ipos in 1:length(fix.posns)){
     ## print(ipos)
      a.posn<-as.integer(names(fix.posns[ipos]))
      num<-fix.posns[ipos]

      alleles <- unlist(strsplit(c(paste(indels.bit[a.posn,"REF"],indels.bit[a.posn,"ALT"], sep=",")),','))
      names(alleles)<-0:(length(alleles)-1)
      ## alleles
      ACs <- unlist(strsplit(indels.bit[a.posn,"AC"],','))
      names(ACs)<-1:(length(alleles)-1) # one AC for each ALT allele
      AFs <- unlist(strsplit(indels.bit[a.posn,"AF"],','))
      names(AFs)<-1:(length(alleles)-1) # one AF for each ALT allele

 

        ############### do AD for all samples are this is the sample for each posn, just replicate once
        ############### sort out the ref and alt alleles 1/1 could be from 0,1 OR 1,2
       ###############  use the PL to decide and test
     
         cols.wanted<-paste(samples.order.in.ALL,"GT",sep=".")
         the.genotypes<-indels.bit[a.posn,cols.wanted]
         the.alleles<-the.genotypes
         the.alleles[the.alleles=="NA" | is.na(the.alleles)]<-"0/0" #just so next step does not overlop the.genotypes==NA are filled with 0,0,0 AD
         the.alleles<-strsplit(the.alleles,split="/")

         ref.allele<-as.integer(unlist(lapply(the.alleles,function(x) x[1]))) # +1 # COLUMN POSN IN AD.LIST
         alt.allele<-as.integer(unlist(lapply(the.alleles,function(x) x[2]))) # +1 # COLUMN POSN IN AD.LIST
      ## above are for the array of samples 
      ##################
      ## ref.allele
      ##  alt.allele
      ## a.posn
      ## indels.bit[a.posn,cols.wanted]
      ## test[a.posn,1:10]

      ############################### define the PL array so can use to get most likely alleles for ties
      
       PL.combo<-matrix(data=NA,nrow=length(alleles),ncol=length(alleles))
      ############ set up the PL structure
      for(iref in 1:length(alleles)){
        for(ialt in iref:length(alleles)){
          i.ref<-as.character(iref) # cause sometimes want to refer to a position by name
          i.alt<-as.character(ialt) # cause sometimes want to refer to a position by name
          PL.combo[iref,ialt]<-paste(iref-1,ialt-1,sep=":") 
       }}
      PL.combo<-as.character(PL.combo)
      PL.combo<-PL.combo[!is.na(PL.combo)]
          ## PL.combo ### the v4 samtools 

         cols.wanted<-paste(samples.order.in.ALL,"PL",sep=".")
         the.PL<- indels.bit[a.posn,cols.wanted]
         the.PL[the.PL=="NA" | is.na(the.PL)]<-paste(rep("3000",times=length(PL.combo)),collapse=",")
         pl.list<-strsplit(the.PL,split=",")
         pl.list<-pl.list[cols.wanted]
         pl.list.length<-length(pl.list)
         pl.list<-as.integer(unlist(pl.list))
         dim(pl.list)<-c(length(PL.combo), pl.list.length)
         pl.list<-t(pl.list)
         colnames(pl.list)<-PL.combo ## p.list and array of sample by allele combos 



       ###### This is done cause not sure what to do if no allel that matches reference
       ######### fix location where have 0/0 or  1/1 or 2/2 choose two most probable or 0 if a tie..
         ties<-which(ref.allele==alt.allele)
         pl.ties<-pl.list[ties,] ## P
         if(length(ties)>0){  ## just in case is no ties which causes an error
         for(ir in 1:length(ties)){
           if(is.null(dim(pl.ties))){pl.ties<-as.matrix(pl.ties);pl.ties<-t(pl.ties)} ## in case inly one ties and get a vector!
           best.pl<-sort(pl.ties[ir,]) ## this sort works as 0:1 reported before 1:2 is they are tied (see PL.combo order)
           the.alleles<-{}
           while(length(the.alleles)<2){
             the.alleles<-sort(as.integer(unique(unlist(strsplit(paste(names(best.pl[1:2]),collapse=":"),split=":")))))
           }
           ref.allele[ties[ir]]<-the.alleles[1]
           alt.allele[ties[ir]]<-the.alleles[2]
        
         }}
      ###################################################
       ## ref.allele ### these contain the best choese of alleles fo the genotypes listed
       ## alt.allele  ### these contain the best choese of alleles fo the genotypes listed 
      ############### do AD for all samples are this is the sample for each posn, just replicate once

      #####get The AD 
         cols.wanted<-paste(samples.order.in.ALL,"AD",sep=".")
         the.AN<- indels.bit[a.posn,cols.wanted]
         the.AN[the.AN=="NA" | is.na(the.AN)]<-paste(rep("0",times=length(alleles)),collapse=",")
         ad.list<-strsplit(the.AN,split=",")
         ad.list<-ad.list[cols.wanted]
         ad.list.length<-length(ad.list)
         ad.list<-as.integer(unlist(ad.list))
         dim(ad.list)<-c(length(alleles), ad.list.length)
         ad.list<-t(ad.list)

    
      ref.counts<-rep(NA,times=dim(ad.list)[1])
      alt.counts<-rep(NA,times=dim(ad.list)[1])
         for(ir in 1:dim(ad.list)[1]){
           ref.counts[ir]<-ad.list[ir,(ref.allele[ir]+1)] ## ad list is ordered in assending allele
           alt.counts[ir]<-ad.list[ir,(alt.allele[ir]+1)]
         }
         the.AN<-paste(ref.counts,alt.counts,sep=",")
         names(the.AN)<-cols.wanted
         
      #################################################################
        
      #####get The PL
      cols.wanted<-paste(samples.order.in.ALL,"PL",sep=".")
      the.PL<-rep(NA,times=dim(pl.list)[1])
         for(ir in 1:dim(pl.list)[1]){
            pl.wanted<-c(paste(ref.allele[ir],ref.allele[ir],sep=":"),
                         paste(ref.allele[ir],alt.allele[ir],sep=":"),
                         paste(alt.allele[ir],alt.allele[ir],sep=":"))
           
           the.PL[ir]<-paste(pl.list[ir,pl.wanted],collapse=":")
         }
         names(the.PL)<-cols.wanted
         
     #################################################################

      ## the.genotypes
      ## the.AN
      ## the.PL
      ## ref.allele
      ## alt.allele
      ## ad.list[1:5,]
      ## pl.list[1:5,]

      if(skip){

        the.cols<-colnames(indels.bit)
          ##### here just flattening
          a.place<-a.posn
          i.ref<-1
          i.alt<-which.max(nchar(as.character(alleles[2:length(alleles)])) )+1 # Choose ALT to be the longest allele
          indels.bit[a.place,"REF"]<-alleles[i.ref]
          indels.bit[a.place,"ALT"]<-alleles[i.alt]
          indels.bit[a.place,"AC"]<-sum(as.numeric(ACs),na.rm=TRUE)
          indels.bit[a.place,"AF"]<-sum(as.numeric(AFs),na.rm=TRUE)
        
          if("RPA" %in% the.cols){indels.bit[a.place,"RPA"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("HRun" %in% the.cols){indels.bit[a.place,"HRun"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("SB" %in% the.cols){indels.bit[a.place,"SB"]<--1} ## cause not defined for mulit alleles and NA set to FALSE in filtering
          if("FS" %in% the.cols){indels.bit[a.place,"FS"]<-0} ## cause not defined for mulit alleles and NA set to FALSE in filtering
         





         ########################################  FIX GT ################################
          
          cols.wanted<-paste(samples.order.in.ALL,"GT",sep=".")
          the.null.genotype<-is.na(the.genotypes)

         the.alleles.new<-the.genotypes
         the.alleles.new[the.alleles.new=="NA"]<-"0/0" #just so next step does not overlop the.genotypes==NA are filled with 0,0,0 AD
         the.alleles.new<-strsplit(the.alleles.new,split="/")

         #ref.allele and alt.allele are built to identify the best reference allele in case when have 2/2 or 3/3 here the ref alleles is as per the original
         ref.allele.new<-as.integer(unlist(lapply(the.alleles.new,function(x) x[1]))) # +1 # COLUMN POSN IN AD.LIST
         alt.allele.new<-as.integer(unlist(lapply(the.alleles.new,function(x) x[2]))) # +1 # COLUMN POSN IN AD.LIST

         ref.allele.new[ref.allele.new!=0]<-1  ## casue will have 2 or 1 for 1/1 or 2/2 2/3 etc
         alt.allele.new[alt.allele.new!=0]<-1

          new.gt<-paste(ref.allele.new,alt.allele.new,sep="/")
          names(new.gt)<-names(the.genotypes)

          
          new.gt[the.null.genotype]<-NA
          indels.bit[a.place,cols.wanted]<-new.gt  
        ##################################################################################
          
        ########################################  FIX AD ################################
        cols.wanted<-paste(samples.order.in.ALL,"AD",sep=".")
        new.AD<-the.AN
        indels.bit[a.place,cols.wanted]<-new.AD
        ##################################################################################

        ########################################  FIX PL ################################
        cols.wanted<-paste(samples.order.in.ALL,"PL",sep=".")
        new.PL<-the.PL
        indels.bit[a.place,cols.wanted]<-new.PL
        ##################################################################################
         




        

      }else{
      
      ## print(alleles)
      iplace<-1 #################### now fix all positions that start with a.posn iplace goes from 1:num
      for(iref in 0:(length(alleles)-1)){
        for(ialt in iref:(length(alleles)-1)){
          if(iref==ialt){next}
          i.ref<-as.character(iref) # cause sometimes want to refer to a position by name
          i.alt<-as.character(ialt) # cause sometimes want to refer to a position by name
          ## print(paste(i.ref,i.alt,sep=":"))
         ## }}
     
         a.place<-a.posn+iplace-1 # -1 because a.posn is the one I want to change

          ## print(a.place)
          ## indels.bit[a.place,cols.wanted]
          ## test[a.place,cols.wanted]
          ## indels.bit[a.place,1:35]

          
          indels.bit[a.place,"REF"]<-alleles[i.ref]
          indels.bit[a.place,"ALT"]<-alleles[i.alt]
          indels.bit[a.place,"AC"]<-ACs[i.alt]
          indels.bit[a.place,"AF"]<-AFs[i.alt]
          indels.bit[a.place,"HRun"]<-0 ## cause not defined for mulit alleles and NA set to FALSE in filtering
          indels.bit[a.place,"SB"]<--1 ## cause not defined for mulit alleles and NA set to FALSE in filtering
          indels.bit[a.place,"TYPE"]<-paste(indels.bit[a.place,"TYPE"],indels.bit[a.place,"POS"],sep=":")

         ########################################  FIX GT ################################

          new.gt<-the.genotypes
          target.GT<- iref==ref.allele & ialt==alt.allele
          
          new.gt[!target.GT]<-"NA"
          if(sum(target.GT)==0){null.genotype<-c(null.genotype,a.place)}
          
          new.gt<-gsub(i.ref,"0",new.gt) ## set the other alleles to reference
          new.gt<-gsub(i.alt,"1",new.gt) ## set the 0/2 to 0/1
          cols.wanted<-paste(samples.order.in.ALL,"GT",sep=".")
          indels.bit[a.place,cols.wanted]<-new.gt

        ##################################################################################
          
        ########################################  FIX AD ################################

        cols.wanted<-paste(samples.order.in.ALL,"AD",sep=".")
        new.AD<-the.AN
        new.AD[!target.GT]<-"NA"
        indels.bit[a.place,cols.wanted]<-new.AD
 
        ##################################################################################

        ########################################  FIX PL ################################

        cols.wanted<-paste(samples.order.in.ALL,"PL",sep=".")
        new.PL<-the.PL
        new.PL[!target.GT]<-"NA"
        indels.bit[a.place,cols.wanted]<-new.PL
 
        ##################################################################################
         
   
         iplace<-iplace+1
          }} #iref ialt iplace fill in same genotypes for location

    } #skip is true the just flattern

     print(a.posn)
     print(dim(indels.bit))
    } #ipos new location
  
#} #ifix
    #if(length(null.genotype)>1){  indels.bit<-indels.bit[-1*null.genotype,]}
    } # indels.bit.ori.size has to do work and fix alleles



 #####################################################################uses just indels.bit from this point ###################################

indels.bit



}  ## end function correct.muli.alleles






############################END FUNCTIONS USED TO READ AND PROCESS VCF


### Define functions ##########
clean.up.brackets<-function(x) {x<-gsub("\\(\\S*","",x,perl=TRUE)} # needed for function below
clean.up<-function(x) {x<-gsub("\\(\\S*","",x,perl=TRUE)}
split.genes<-function(x,delimit.by) {x<-unlist(strsplit(x,split=delimit.by))}


# get.gene.names(geneanno.table,ann.order)
## local.target.table<-geneanno.table
## col.names<-ann.order

get.gene.names<-function(local.target.table,col.names){  # give a table the the colnames to be used to get genes -> hget back a list names with col names vis a.list[["refGene::gene"]][1:50]
  ## remove the names on intergenic hits so dont greeop conserved distal regions of genes with the gene itself 
  a.list<-list()
  # i<-1
  for(i in 1:length(col.names)){
    local.gene.names<-local.target.table[,col.names[i]]
    to.fix1<-grepl("(",local.gene.names,fixed=TRUE)
    to.fix2<-grepl(",",local.gene.names,fixed=TRUE)
    to.fix3<-grepl(";",local.gene.names,fixed=TRUE)
    to.fix<-to.fix1 | to.fix2
    sum(to.fix)
    local.gene.names<-strsplit(local.gene.names,split="),") ### have Gene(transcript info, transcript info),Gene2(transcript info) OR  Gene(dist=XXX),Gene2(Dist=yyy)

    ####all lapply beow were sapply before but this was causing problems 
    local.gene.names[to.fix1]<-lapply(local.gene.names[to.fix1],clean.up.brackets)
    local.gene.names[to.fix2]<-lapply(local.gene.names[to.fix2],split.genes,delimit.by=",")
    local.gene.names[to.fix3]<-lapply(local.gene.names[to.fix3],split.genes,delimit.by=";")
    ################ remove below  to count transscript hitts and not gene hits
    poly.gene<-unlist(lapply(local.gene.names,length))>1
    local.gene.names[poly.gene]<-lapply(local.gene.names[poly.gene],unique)
    ################

    ########## cancel out intergenic hits
    location<-local.target.table[,gsub("gene$","location",col.names[i])]
    not.gene<-location %in% c("intergenic")
    labels<-paste("PLACE_",rownames(local.target.table)[not.gene],sep="")
    local.gene.names[not.gene]<-labels
    
    a.list[col.names[i]]<-list(local.gene.names)
  }
  a.list
}

##  strange<-"10a.872B9  10a)"
## posn<-grep("10a.872B9",geneanno.table[,"knownGene::gene"])

## posn<-grep(" \\(",geneanno.table[,"knownGene::gene"])

## geneanno.table[posn,]

## grep(strange,unlist(gene.names[["knownGene::gene"]]))
## unlist(the.genes[["knownGene::gene"]])[116511]

#### typical input data for above 
## [1] "CHURC1-FNTB,FNTB,MAX"

## [[2]]
## [1] "NONE(dist=NONE"     "MIR622(dist=83084)"

## [[3]]
## [1] "TMEM56(NM_001199679:exon4:c.246-4->T,NM_152487:exon4:c.246-4->T" "TMEM56-RWDD3(NM_001199691:exon4:c.246-4->T)"                    

## [[4]]
## [1] "BAGE,BAGE2,BAGE3,BAGE4,BAGE5"

## [[5]]
## [1] "BAGE(NM_001187:exon2:c.14+1A>G"   "BAGE4(NM_181704:exon2:c.14+1A>G"  "BAGE5(NM_182484:exon2:c.14+1A>G)"

######################### subroutines used in match.cols.and.collapse
delimit.collapse<-function(x,tab,tab.match.col){
    posns<-match(x,tab[,tab.match.col])
    missing<-is.na(posns)
    data<-as.matrix(tab[posns[!missing],])
    if(dim(data)[1]==0){data=matrix(data=NA,nrow=1,ncol=dim(tab)[2]);colnames(data)<-colnames(tab)} #no match found so return NA's data[1,]<-rep(NA,times=dim(tab)[2])
    apply(data,2,function(x) paste(x,collapse="::")) }

boolean.collapse<-function(x,tab,tab.match.col) {
    posns<-match(x,tab[,tab.match.col])
    missing<-is.na(posns)
    data<-tab[posns[!missing],]
    if(dim(data)[1]==0){data=matrix(data=NA,nrow=1,ncol=dim(tab)[2]);colnames(data)<-colnames(tab)} #no match found so return NA's
    the.dim<-dim(data)
    the.match.col<-paste(data[,1],collapse="::") # in begining forced first col to be the match column
    c(the.match.col, apply(as.matrix(data[,2:the.dim[2]]),2,function(x) sum(x,na.rm=TRUE)>0))
  }

  sum.collapse<-function(x,tab,tab.match.col) {
    posns<-match(x,tab[,tab.match.col])
    missing<-is.na(posns)
    data<-tab[posns[!missing],]
    if(dim(data)[1]==0){data=matrix(data=NA,nrow=1,ncol=dim(tab)[2]);colnames(data)<-colnames(tab)} #no match found so return NA's
    the.dim<-dim(data)
    the.match.col<-paste(data[,1],collapse="::") # in begining forced first col to be the match column
    c(the.match.col, apply(as.matrix(data[,2:the.dim[2]]),2,function(x) sum(as.numeric(x),na.rm=TRUE)))
  }

  max.collapse<-function(x,tab,tab.match.col) {
    posns<-match(x,tab[,tab.match.col])
    missing<-is.na(posns)
    data<-tab[posns[!missing],]
    if(dim(data)[1]==0){data=matrix(data=NA,nrow=1,ncol=dim(tab)[2]);colnames(data)<-colnames(tab)} #no match found so return NA's
    the.dim<-dim(data)
    the.match.col<-paste(data[,1],collapse="::") # in begining forced first col to be the match column
    c(the.match.col, apply(as.matrix(data[,2:the.dim[2]]),2,function(x) max(as.numeric(x),na.rm=TRUE)))
  }

  min.collapse<-function(x,tab,tab.match.col) {
    posns<-match(x,tab[,tab.match.col])
    missing<-is.na(posns)
    data<-tab[posns[!missing],]
    if(dim(data)[1]==0){data=matrix(data=NA,nrow=1,ncol=dim(tab)[2]);colnames(data)<-colnames(tab)} #no match found so return NA's
    the.dim<-dim(data)
    the.match.col<-paste(data[,1],collapse="::") # in begining forced first col to be the match column
    c(the.match.col, apply(as.matrix(data[,2:the.dim[2]]),2,function(x) min(as.numeric(x),na.rm=TRUE)))
  }


## list1<-gene.names
## list1.match.col<-"refGene::gene"
## delimit.by<-"::"
collapse.gene.names<-function(list1,list1.match.col,delimit.by){          ##  list1<-gene.names   list1.match.col<-"refGene::gene"
  if(is.null(names(list1))){local.list<-list1}else{local.list<-list1[[list1.match.col]]} # a single list has no name
  num.genes<-unlist(lapply(local.list,length))                            
  poly.genes<-num.genes>1 
  a.desc.table<-rep(NA,times=length(local.list))
  single.genes<-unlist(local.list[!poly.genes])
  multi.genes<-unlist(sapply(local.list[poly.genes],function(x) paste(x,collapse=delimit.by)))
  a.desc.table[!poly.genes]<-single.genes
  if(length(multi.genes)!=0){a.desc.table[poly.genes]<-multi.genes}
  a.desc.table
 }
                       

test.wanted.mutation<-function(type.vector,mutations,delimit.by){
  
  all.mutations<-rep(NA,times=length(type.vector))
  missing<-is.na(type.vector)
  type.vector<- type.vector[!missing]
  names(type.vector)<-1:length(type.vector)
  
  type.list<-strsplit(type.vector,split=delimit.by)
  list.element.lengths<-unlist(lapply(type.list,length)) ## check types
  sum(list.element.lengths>1)
  if( sum(list.element.lengths==0) ){
    bad.mutation<-type.vector %in% mutations
    all.mutations[!missing]<-bad.mutation
  }else{ # multipe types in one exon
    indel.index<-rep(1:length(list.element.lengths),times=list.element.lengths)
    flat.type.list<-unlist(type.list)
    bad.mutation<-flat.type.list %in% mutations
    bad.mutation<-tapply(bad.mutation,indel.index,function(x) sum(x,na.rm=TRUE)>0 )
    if(sum(names(bad.mutation)!=1:length( bad.mutation))!=0){ bad.mutation<- bad.mutation[as.character(1:length(bad.mutation))]} # reoder if necessary
    all.mutations[!missing]<-bad.mutation
  }
  all.mutations
}

#poly.scores<-match.cols.and.collapse(list(key.indels.poly),1,poly.scores,"key.pholy",colnames(poly.scores),"delimit")

## list1<-list(key.indels.poly)
## list1.match.col<-1
## tab2<-poly.scores
## tab2.match.col<-"key.pholy"
## tab2.collapse.cols<-colnames(poly.scores)
## gene.collapse<-"delimit"

#(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit")
#gene.names,"ensGene::gene",gene.desc,"ensembl_gene_id",gene.ann.wanted,"delimit"
## list1<-gene.names
## list1.match.col<-"refGene::gene"
## tab2<-a.bone.file
## tab2.match.col<-colnames(a.bone.file)[1]
## tab2.collapse.cols<-colnames(a.bone.file)
## gene.collapse<-"delimit"
####################### match genes to another column #########################################
### this method i faser when doing collapse with text and table where as match.cols.and.collapse.fast works better with interger collapse
match.cols.and.collapse<-function(list1,list1.match.col,tab2,tab2.match.col,tab2.collapse.cols,gene.collapse){ # match genes names (might be more than one) in list  to a table and collapse 
  num.genes<-unlist(lapply(list1[[list1.match.col]],length))                                     # gived back table of dimension of indel table gene.collapse="delimit","boolean","sum"
  poly.genes<-num.genes>1 ## multpe gene hits                                                    # IMPORTANT: tab2.collapse.cols also must contain the tab2.match.col !!
 #  print(gene.collapse)
 #  print(tab2.match.col)

                                        # is send it a vector then itwill assigne "genes" from names and values in vector
 ## if( (class(tab2)=="data.frame" | class(tab2)=="matrix")){if(dim(tab2)[2]<=1){tab2<-tab2<-data.frame(genes=tab2,value=tab2);colnames(tab2)<-c("genes","value");tab2.match.col<-"genes";tab2.collapse.cols<-c("value")}}
  
  if(class(tab2)=="logical" | class(tab2)=="integer" | class(tab2)=="character" ){tab2<-data.frame(genes=names(tab2),value=tab2);tab2.match.col<-"genes";tab2.collapse.cols<-c("genes","value")}

  if( !(gene.collapse %in% c("delimit","boolean","sum","max","min"))){gene.collapse<-"delimit"}
  # make sure the match column is the first columna as is in tab2 if not assume it is the rowanmes of the table
  if(tab2.match.col %in% colnames(tab2)){order<-c(tab2.match.col,colnames(tab2)[colnames(tab2)!=tab2.match.col]) 
                                        if(length(order)>1){tab2<-tab2[,order]}else{tab2<-data.frame(tab2[,order],stringsAsFactors=FALSE) ; colnames(tab2)<-order}
                                       }else{
                                      tab2.collapse.cols<-c(tab2.match.col,colnames(tab2))
                                      tab2<-cbind(rownames(tab2),tab2)
                                      colnames(tab2)<-tab2.collapse.cols
                                    }
                                         
# print(tab2[1:5,])
  if(is.null(dim(tab2))){tab2<-t(as.matrix(tab2))}
#map the single genes first
  a.desc.table<-matrix(data=NA,nrow=length(list1[[list1.match.col]]),ncol=length(tab2.collapse.cols))
  colnames(a.desc.table)<-tab2.collapse.cols
  a.desc.table<-as.data.frame(a.desc.table,stringsAsFactors=FALSE)

  posns<-match(unlist(list1[[list1.match.col]][!poly.genes]),tab2[,tab2.match.col]) ## match the single items first to get the out-of -the -way
  missing<-is.na(posns)
  #sum(missing)  sum(!missing)

 
  if(length(tab2.collapse.cols)==1){   ## have problem is gene.desc.tabel becomes a vector
     a.desc.table[!poly.genes,tab2.collapse.cols[1]][!missing]<-as.character(tab2[posns[!missing],tab2.collapse.cols[1]])
   }else{
     
  for(i in 1:length(tab2.collapse.cols)){  ## mave a matrix so no issues
    a.desc.table[!poly.genes,][!missing,tab2.collapse.cols[i]]<-as.character(tab2[posns[!missing],tab2.collapse.cols[i]])
        }
         }

  if(!sum(poly.genes)==0){
#print(gene.collapse)
  if(gene.collapse=="boolean"){
     other.genes<-sapply(list1[[list1.match.col]][poly.genes],boolean.collapse,tab=tab2,tab.match.col=tab2.match.col)
        }else if(gene.collapse=="max"){
     other.genes<-sapply(list1[[list1.match.col]][poly.genes],max.collapse,tab=tab2,tab.match.col=tab2.match.col)
             }else if(gene.collapse=="min"){
     other.genes<-sapply(list1[[list1.match.col]][poly.genes],min.collapse,tab=tab2,tab.match.col=tab2.match.col)
             }else if(gene.collapse=="sum"){
     other.genes<-sapply(list1[[list1.match.col]][poly.genes],sum.collapse,tab=tab2,tab.match.col=tab2.match.col)
   }else{
#     system.time(
     other.genes<-sapply(list1[[list1.match.col]][poly.genes],delimit.collapse,tab=tab2,tab.match.col=tab2.match.col)
#)   
   }  ## default is demilt collapse

                     
      other.genes<-t(other.genes)
      other.genes<-as.data.frame(other.genes,stringsAsFactors=FALSE)
      colnames(other.genes)<-colnames(tab2)

  if(sum(poly.genes)!=dim(other.genes)[1]){print("ERROR could not match polygenes")} ## must be true


  if(length(tab2.collapse.cols)==1){   ## have problem is gene.desc.tabel becomes a vector
     a.desc.table[poly.genes,tab2.collapse.cols[1]]<-other.genes[,tab2.collapse.cols[1]]
   }else{
     
  for(i in 1:length(tab2.collapse.cols)){
    a.desc.table[poly.genes,][,tab2.collapse.cols[i]]<-other.genes[,tab2.collapse.cols[i]]
  }
         } # length(tab2.collapse.cols)!=1

   } # no poly genes: sum(poly.genes)!=0

  a.desc.table
}

## a.desc.table[poly.genes,]
## bad1<-a.desc.table[,1]=="NA"
## bad2<-is.na(a.desc.table[,1])
## bad<-bad1 | bad2
## a.desc.table[!bad,]
####################### match genes to another column #########################################
## list.element.lengths,indel.index,flat.gene.list,counts.results[["num.wanted.muts"]],"genes",c("genes","value"),"max"
#match.cols.and.collapse.fast(list.element.lengths,indel.index,flat.gene.list,counts.results[["gene.counts"]],"genes",colnames(counts.results[["gene.counts"]]),"max")
#list.element.lengths,indel.index,flat.gene.list,counts.results[["num.wanted.muts"]],"genes",c("genes","value"),"max")


## list.length<-list.element.lengths
## index<-indel.index
## flat.list<-flat.gene.list
## tab2<-counts.results[["num.wanted.muts"]]
## tab2.match.col<-"genes"
## tab2.collapse.cols<-c("genes","value") # colnames(counts.results[["gene.counts"]])
## gene.collapse<-"max"

match.cols.and.collapse.fast<-function(list.lengths,index,flat.list,tab2,tab2.match.col,tab2.collapse.cols,gene.collapse){ # match genes names (might be more than one) in list  to a table and collapse 
### the idea here is used a flatted list and tapply
### IF A VECTOS IS USED THE THE NAMES OF THE VECTOR MUST BE THE GENE NAMES OR THE tab2.match.col
    ## list.element.lengths<-unlist(lapply(the.genes[[annotation.labels[i]]],length))        # list.lentghs
    ## indel.index<-rep(1:length(list.element.lengths),times=list.element.lengths)           # index
    ## flat.gene.list<-unlist(the.genes[[annotation.labels[i]]])                              # flat.list

#   if( (class(tab2)=="data.frame" | class(tab2)=="matrix")){if(dim(tab2)[2]<=1){tab2<-tab2<-data.frame(genes=tab2,value=tab2);colnames(tab2)<-c("genes","value");tab2.match.col<-"genes";tab2.collapse.cols<-c("value")}}
   
  if(class(tab2)=="logical" | class(tab2)=="integer" | class(tab2)=="character"){tab2<-data.frame(genes=names(tab2),value=tab2);tab2.match.col<-"genes";tab2.collapse.cols<-c("genes","value")}

  if( !(gene.collapse %in% c("delimit","boolean","sum","max","min","hit_count")) ){gene.collapse<-"delimit"}
  # make sure the match column is the first columna as is in tab2 if not assume it is the rowanmes of the table
  if(tab2.match.col %in% colnames(tab2)){order<-c(tab2.match.col,colnames(tab2)[colnames(tab2)!=tab2.match.col]) ;  tab2<-tab2[,order]
                                       }else{
                                      tab2.collapse.cols<-c(tab2.match.col,colnames(tab2))
                                      tab2<-cbind(rownames(tab2),tab2)
                                      colnames(tab2)<-tab2.collapse.cols
                                    }

   if(is.null(dim(tab2))){tab2<-t(as.matrix(tab2))} ## incase has become a vector
  
  posns<-match(flat.list,tab2[,tab2.match.col]) ## assumes tab2[,tab2.match.col] is a list of unique elements (no duplicates) 
  
  to.do<-tab2.collapse.cols[ tab2.collapse.cols !=tab2.match.col]
#  to.do.posns<-match(to.do,colnames(tab2))
#  a.desc.table<-matrix(data=,nrow=length(unique(index)),ncol=length(to.do))
   a.desc.table<-{}
# print(gene.collapse)
  
  if(gene.collapse=="boolean"){

a.desc.table<-foreach(info2=iter(subset(tab2[posns,],select=to.do), by='col',chunksize=1), .combine=cbind,.multicombine=TRUE,.inorder=TRUE) %dopar% tapply(info2,index,function(x) sum(x,na.rm=TRUE)>0 )

           }else if(gene.collapse=="max"){

a.desc.table<-foreach(info2=iter(subset(tab2[posns,],select=to.do), by='col',chunksize=1), .combine=cbind,.multicombine=TRUE,.inorder=TRUE) %dopar% tapply(info2,index,function(x) max(as.numeric(x),na.rm=TRUE))  

            }else if(gene.collapse=="hit_count"){
            
a.desc.table<-foreach(info2=iter(subset(tab2[posns,],select=to.do), by='col',chunksize=1), .combine=cbind,.multicombine=TRUE,.inorder=TRUE) %dopar% tapply(info2,index,function(x) sum((x==0 & !is.na(x)),na.rm=TRUE) ) 


             }else if(gene.collapse=="min"){

a.desc.table<-foreach(info2=iter(subset(tab2[posns,],select=to.do), by='col',chunksize=1), .combine=cbind,.multicombine=TRUE,.inorder=TRUE) %dopar% tapply(info2,index,function(x) min(as.numeric(x),na.rm=TRUE))
               
             }else if(gene.collapse=="sum"){

a.desc.table<-foreach(info2=iter(subset(tab2[posns,],select=to.do), by='col',chunksize=1), .combine=cbind,.multicombine=TRUE,.inorder=TRUE) %dopar% tapply(info2,index,function(x) sum(as.numeric(x),na.rm=TRUE))
                       
              }else{

a.desc.table<-foreach(info2=iter(subset(tab2[posns,],select=to.do), by='col',chunksize=1), .combine=cbind,.multicombine=TRUE,.inorder=TRUE) %dopar% tapply(info2,index,function(x) paste(as.character(x),collapse="::"))

                      }  ## default is demilt collapse

   if(length(to.do)>1){colnames(a.desc.table)<-to.do} # colnames(data)<-to.do

  a.desc.table
}
########################################################################################
########################################################################################
##
match.cols.and.collapse.fast.old<-function(list.lengths,index,flat.list,tab2,tab2.match.col,tab2.collapse.cols,gene.collapse){ # match genes names (might be more than one) in list  to a table and collapse 
### the idea here is used a flatted list and tapply
### IF A VECTOS IS USED THE THE NAMES OF THE VECTOR MUST BE THE GENE NAMES OR THE tab2.match.col
    ## list.element.lengths<-unlist(lapply(the.genes[[annotation.labels[i]]],length))        # list.lentghs
    ## indel.index<-rep(1:length(list.element.lengths),times=list.element.lengths)           # index
    ## flat.gene.list<-unlist(the.genes[[annotation.labels[i]]])                              # flat.list

#   if( (class(tab2)=="data.frame" | class(tab2)=="matrix")){if(dim(tab2)[2]<=1){tab2<-tab2<-data.frame(genes=tab2,value=tab2);colnames(tab2)<-c("genes","value");tab2.match.col<-"genes";tab2.collapse.cols<-c("value")}}
   
  if(class(tab2)=="logical" | class(tab2)=="integer" | class(tab2)=="character"){tab2<-data.frame(genes=names(tab2),value=tab2);tab2.match.col<-"genes";tab2.collapse.cols<-c("genes","value")}

  if( !(gene.collapse %in% c("delimit","boolean","sum","max","min","hit_count")) ){gene.collapse<-"delimit"}
  # make sure the match column is the first columna as is in tab2 if not assume it is the rowanmes of the table
  if(tab2.match.col %in% colnames(tab2)){order<-c(tab2.match.col,colnames(tab2)[colnames(tab2)!=tab2.match.col]) ;  tab2<-tab2[,order]
                                       }else{
                                      tab2.collapse.cols<-c(tab2.match.col,colnames(tab2))
                                      tab2<-cbind(rownames(tab2),tab2)
                                      colnames(tab2)<-tab2.collapse.cols
                                    }

   if(is.null(dim(tab2))){tab2<-t(as.matrix(tab2))} ## incase has become a vector
  
  posns<-match(flat.list,tab2[,tab2.match.col]) ## assumes tab2[,tab2.match.col] is a list of unique elements (no duplicates) 
  
  to.do<-tab2.collapse.cols[ tab2.collapse.cols !=tab2.match.col]
  a.desc.table<-{}
# print(gene.collapse)
  
  if(gene.collapse=="boolean"){
          for (i in 1:length(to.do)){    
         data<-tapply(tab2[posns,to.do[i]],index,function(x) sum(x,na.rm=TRUE)>0 )
         if(sum(names(data)!=1:length(data))!=0){data<-data[as.character(1:length(data))]}
         if(i==1){a.desc.table<-data}else{a.desc.table<-cbind(a.desc.table,data)}}

           }else if(gene.collapse=="max"){
          for (i in 1:length(to.do)){    
         data<-tapply(tab2[posns,to.do[i]],index,function(x) max(as.numeric(x),na.rm=TRUE))
         if(sum(names(data)!=1:length(data))!=0){data<-data[as.character(1:length(data))]}
         if(i==1){a.desc.table<-data}else{a.desc.table<-cbind(a.desc.table,data)} # a.desc.table<-cbind(a.desc.table,data)
       }

            }else if(gene.collapse=="hit_count"){
          for (i in 1:length(to.do)){    
         data<-tapply(tab2[posns,to.do[i]],index,function(x) sum((x==0 & !is.na(x)),na.rm=TRUE) )
         if(sum(names(data)!=1:length(data))!=0){data<-data[as.character(1:length(data))]}
         if(i==1){a.desc.table<-data}else{a.desc.table<-cbind(a.desc.table,data)}}

             }else if(gene.collapse=="min"){
          for (i in 1:length(to.do)){    
         data<-tapply(tab2[posns,to.do[i]],index,function(x) min(as.numeric(x),na.rm=TRUE))
         if(sum(names(data)!=1:length(data))!=0){data<-data[as.character(1:length(data))]}
         if(i==1){a.desc.table<-data}else{a.desc.table<-cbind(a.desc.table,data)}
       }
               
             }else if(gene.collapse=="sum"){
          for (i in 1:length(to.do)){    
         data<-tapply(tab2[posns,to.do[i]],index,function(x) sum(as.numeric(x),na.rm=TRUE))
         if(sum(names(data)!=1:length(data))!=0){data<-data[as.character(1:length(data))]}
         if(i==1){a.desc.table<-data}else{a.desc.table<-cbind(a.desc.table,data)}}
                       
              }else{
          for (i in 1:length(to.do)){    
         data<-tapply(tab2[posns,to.do[i]],index,function(x) paste(as.character(x),collapse="::"))
         if(sum(names(data)!=1:length(data))!=0){data<-data[as.character(1:length(data))]}
         if(i==1){a.desc.table<-data}else{a.desc.table<-cbind(a.desc.table,data)}}

                      }  ## default is demilt collapse

   if(length(to.do)>1){colnames(a.desc.table)<-to.do}

  a.desc.table
}
########################################################################################



expand.labels.to.samples.complex<-function(allele.labels,the.sample.order,paste.after=FALSE,seperator="."){
    sample.class<-as.character( t(matrix(data=rep(the.sample.order,times=length(allele.labels)),nrow=length(the.sample.order),ncol=length(allele.labels),byrow=FALSE)) )
    if(paste.after){sample.class<-paste(allele.labels,sample.class,sep=seperator)}else{
    sample.class<-paste(sample.class,allele.labels,sep=seperator)}
    sample.class}


expand.labels.to.samples<-function(allele.labels,the.sample.order,paste.after=FALSE){
    sample.class<-as.character( t(matrix(data=rep(the.sample.order,times=length(allele.labels)),nrow=length(the.sample.order),ncol=length(allele.labels),byrow=FALSE)) )
    if(paste.after){sample.class<-paste(allele.labels,sample.class,sep=".")}else{
    sample.class<-paste(sample.class,allele.labels,sep=".")}
    sample.class}

expand.quantity.to.samples<-function(allele.labels,the.sample.order){
    sample.class<-rep(allele.labels,times=length(the.sample.order))
    sample.class}

build.key<-function(table,key.cols,add.chr.label=FALSE){
  options(scipen=300)
  if(is.null(dim(table))){table<-as.matrix(table)} # in casea vector sent
  if(length(key.cols)<1){print("FAIL no keys columns specified");key<-1:dim(table)[1]}else{
    for (i in 1:length(key.cols)){
      if(i==1){key<-table[,key.cols[i]]}else{
      key<-paste(key,table[,key.cols[i]],sep=":")}
    }}
  if(add.chr.label){key<-paste("chr",key,sep="")}
         key}

build.key.delimit<-function(table,key.cols,add.chr.label=FALSE,delim=":"){
   options(scipen=300)
  if(is.null(dim(table))){table<-as.matrix(table)} # in casea vector sent
  if(length(key.cols)<1){print("FAIL no keys columns specified");key<-1:dim(table)[1]}else{
    for (i in 1:length(key.cols)){
      if(i==1){key<-table[,key.cols[i]]}else{
      key<-paste(key,table[,key.cols[i]],sep=delim)}
    }}
  if(add.chr.label){key<-paste("chr",key,sep="")}
         key}

#the.table,the.genes,annotation.labels[i],the.global,the.filter,the.samples,the.attributes

## a.table<- the.table
## the.local.genes<- the.genes
## gene.match.col<-annotation.labels[i]
## the.local.global<- the.global
##  the.local.filter<- the.filter
## the.local.samples<- the.samples
## the.local.attributes<-the.attributes


pheno.count.function<-function(x,the.local.global,the.local.filter,the.local.genes,gene.match.col){
    x<- x & the.local.global & the.local.filter  ## column wise set not interesting 
    poly.genes<-unlist(lapply(the.local.genes[[gene.match.col]],length))
    x<-rep(x,times=poly.genes)
    all.genes<-unlist(the.local.genes[[gene.match.col]])
    counts<-tapply(x,all.genes,sum)
    counts
  }


get.count.data<-function(a.table,the.local.genes,gene.match.col,the.local.global,the.local.filter,the.local.samples,the.local.attributes){
  
#  pheno.aff.gene.counts<-apply(as.matrix(a.table),2,pheno.count.function,the.local.global,the.local.filter,the.local.genes,gene.match.col)  ### this

   pheno.aff.gene.counts<-foreach(info2=iter(as.matrix(a.table), by='col',chunksize=1), .combine=cbind,.multicombine=TRUE,.inorder=TRUE) %dopar% pheno.count.function(info2,the.local.global,the.local.filter,the.local.genes,gene.match.col)

   colnames(pheno.aff.gene.counts)<- colnames(a.table)


  if(is.null(dim(pheno.aff.gene.counts))){
    the.col.names<-names(pheno.aff.gene.counts)
    the.row.names<-unique(unlist(the.local.genes[[gene.match.col]]))
    pheno.aff.gene.counts<-matrix(data=pheno.aff.gene.counts,nrow=1,ncol=length(pheno.aff.gene.counts))
    rownames(pheno.aff.gene.counts)<-the.row.names
    colnames(pheno.aff.gene.counts)<-the.col.names
  }
    
  gene.counts<-data.frame(key=rownames(pheno.aff.gene.counts),stringsAsFactors=FALSE)
  
  for(i in 1:length(the.local.samples)){
    one.sample<-expand.labels.to.samples(the.local.attributes,the.local.samples[i])
    one.count<-apply(as.matrix(pheno.aff.gene.counts[,one.sample]),1,sum)
    posns<-match(gene.counts[,"key"],names(one.count))
    missing<- is.na(posns)
    gene.counts[!missing,the.local.samples[i]]<-one.count[posns[!missing]]
  } #cum over attributes here so have number that pass per gene

## warning gene.counts can go from table to vector if only one column is left losing the rownames
  rownames.gene.counts<-gene.counts[,"key"]
  colnames.gene.counts<-colnames(gene.counts)[colnames(gene.counts) !="key"]
  gene.counts<-gene.counts[,colnames(gene.counts) != "key"] 
  gene.counts<-as.matrix(gene.counts)
  rownames(gene.counts)<-rownames.gene.counts 
  colnames(gene.counts)<-colnames.gene.counts

#  num.wanted.muts<-apply(gene.counts,1,sum) # total filtered hits per gene per sample


  num.wanted.muts<-apply(gene.counts,1,function(x){
    x<-x==0
    length(the.local.samples)-sum(x)
  })  # total filtered hits per gene per sample
  
  phenotype.hit<-apply(gene.counts,1,function(x){
    x<-x>0
    sum(x)==length(the.local.samples)
  }) 
 list(phenotype.hit=phenotype.hit,gene.counts=gene.counts,num.wanted.muts=num.wanted.muts)

}


 get.count.data.old<-function(a.table,the.local.genes,gene.match.col,the.local.global,the.local.filter,the.local.samples,the.local.attributes){
  
  pheno.aff.gene.counts<-apply(as.matrix(a.table),2,function(x){
    x<- x & the.local.global & the.local.filter  ## column wise set not interesting 
    poly.genes<-unlist(lapply(the.local.genes[[gene.match.col]],length))
    x<-rep(x,times=poly.genes)
    all.genes<-unlist(the.local.genes[[gene.match.col]])
    counts<-tapply(x,all.genes,sum)
    counts
  })  ### this

  

  if(is.null(dim(pheno.aff.gene.counts))){
    the.col.names<-names(pheno.aff.gene.counts)
    the.row.names<-unique(unlist(the.local.genes[[gene.match.col]]))
    pheno.aff.gene.counts<-matrix(data=pheno.aff.gene.counts,nrow=1,ncol=length(pheno.aff.gene.counts))
    rownames(pheno.aff.gene.counts)<-the.row.names
    colnames(pheno.aff.gene.counts)<-the.col.names
  }
    
  gene.counts<-data.frame(key=rownames(pheno.aff.gene.counts),stringsAsFactors=FALSE)
  
  for(i in 1:length(the.local.samples)){
    one.sample<-expand.labels.to.samples(the.local.attributes,the.local.samples[i])
    one.count<-apply(as.matrix(pheno.aff.gene.counts[,one.sample]),1,sum)
    posns<-match(gene.counts[,"key"],names(one.count))
    missing<- is.na(posns)
    gene.counts[!missing,the.local.samples[i]]<-one.count[posns[!missing]]
  } #cum over attributes here so have number that pass per gene

## warning gene.counts can go from table to vector if only one column is left losing the rownames
  rownames.gene.counts<-gene.counts[,"key"]
  colnames.gene.counts<-colnames(gene.counts)[colnames(gene.counts) !="key"]
  gene.counts<-gene.counts[,colnames(gene.counts) != "key"] 
  gene.counts<-as.matrix(gene.counts)
  rownames(gene.counts)<-rownames.gene.counts 
  colnames(gene.counts)<-colnames.gene.counts

#  num.wanted.muts<-apply(gene.counts,1,sum) # total filtered hits per gene per sample


  num.wanted.muts<-apply(gene.counts,1,function(x){
    x<-x==0
    length(the.local.samples)-sum(x)
  })  # total filtered hits per gene per sample
  
  phenotype.hit<-apply(gene.counts,1,function(x){
    x<-x>0
    sum(x)==length(the.local.samples)
  }) 
 list(phenotype.hit=phenotype.hit,gene.counts=gene.counts,num.wanted.muts=num.wanted.muts)

}
      
## pheno.qual.Unaff<-CountOverGene.withFilterAndGlobal(as.matrix(quality.thresh[,dom.Unaff.GT]),the.Unaffected,the.attributes,gene.names,filter.table,a.filter.cols.maf,"OR",filter.action="EXCLUDE",global.vars=c("wanted.muts","good.qual") )
## pheno.qual<-CountOverGene.withFilterAndGlobal(the.table<-as.matrix(quality.thresh[,dom.aff.GT]),the.affected,the.attributes,gene.names,filter.blank,c("all.true"),"OR",filter.action="INCLUDE",global.vars=c("good.qual") )filter.action="INCLUDE"
## the.genes<-gene.names
## the.samples<-the.affected
### the.filter.table<-filter.blank
## (as.matrix(quality.thresh[,dom.aff.GT]),the.affected,the.attributes,gene.names,filter.table,a.filter.cols.maf,"AND",filter.action="INCLUDE",global.vars=c("wanted.muts","good.qual") )

#input,the.affected,the.attributes,gene.names,filter.table,a.filter.cols.maf,"AND",filter.action="INCLUDE",global.vars=c("wanted.muts","good.qual")
## (input,the.affected,the.attributes,gene.names,filter.table,a.filter.cols.maf,"OR",filter.action="EXCLUDE",global.vars=c("wanted.muts","good.qual") )
## input,the.affected,the.attributes,gene.names,filter.table,a.filter.cols.maf,"OR",filter.action="EXCLUDE",global.vars=c("wanted.muts","good.qual.filter")
## input,the.affected,the.attributes,gene.names,filter.table,a.filter.cols.maf,"OR",filter.action="EXCLUDE",global.vars=c("wanted.muts","good.qual","ok.in.group")

## the.table<-input
## the.samples<-the.affected
## the.attributes<-the.attributes
## the.genes<-gene.names
## the.filter.table<-filter.table
## filter.cols<-a.filter.cols.maf
## filter.cols.combine<-"OR"
## filter.action="EXCLUDE"
## global.vars<-use.the.globals

#

CountOverGene.withFilterAndGlobal<-function(the.table,the.samples,the.attributes,the.genes,the.filter.table,filter.cols,filter.cols.combine="OR",filter.action="EXCLUDE",global.vars=c() ){


if(dim(the.table)[2]==0 | length(the.samples)==0){  #In case no sample
  nrows<-dim(the.table)[1]
  hit<-rep(FALSE,times=nrows)
  hit.count<-rep(0,times=nrows)
  gene.hit.count<-rep(0,times=nrows);gene.hit.count<-as.matrix(gene.hit.count)
}else{
  
  ## This function takes the quality threshold table made of samples and attributes filter with the filter table and applys & for the gloabl.vars
  ## Pheno hits are where the samples all have a filter hit in a gene
  ## the attributes are ones that would survive matrix mutliplication : ( quality.thresh[,dom.aff.GT] & quality.thresh[,dom.aff.GT.DP])
  ## filter.cols.combine="OR" or "AND"
  ## filter.action="EXCLUDE" or "INCLUDE
  ## global.vars<-c("wanted.muts.fil") ## names of columns vectors use & to combine
  ## include error if filter.cols not in the filter table
  ### eg exclude novel use OR and EXCLUDE and set filter.cols to tests is observed

  ############ set any NA's to true
  the.NAs<-is.na(the.filter.table[,filter.cols])
  count.col.NAs<-apply(as.matrix(the.NAs),2,sum)# one filter .col gives a vector
  found.col.NAs<-count.col.NAs>0
  if(sum(found.col.NAs)>0){
    for(i in 1:length(found.col.NAs)){
      if(found.col.NAs[i]){
        the.filter.table[,filter.cols[i]][is.na(the.filter.table[,filter.cols[i]])]<-TRUE
      }}}
  ##########################
  
  
  num.trues<-apply(as.matrix(the.filter.table[,filter.cols]),1,function(x) sum(x,na.rm=TRUE))
  if(filter.cols.combine=="OR"){the.filter<-num.trues >0} # ANY test is true the filter out 
  if(filter.cols.combine=="AND"){the.filter<-num.trues==length(filter.cols)} # ALL test must be true the filter out 

  if(filter.action=="EXCLUDE"){the.filter<-!the.filter}

  if(length(global.vars)>0){
    the.global<-eval(as.name(global.vars[1]))
    the.global[is.na(the.global)]<-TRUE
    if(length(global.vars)!=1){
      for(i in 2:(length(global.vars))){
        the.next<-eval(as.name(global.vars[i]))
        the.next[is.na(the.next)]<-TRUE
        the.global<-the.global & the.next
        }
  }} ## combine all the global variables
    

  annotation.labels<-names(the.genes)
  the.objects<-gsub("::gene","",annotation.labels)
#  the.objects<-paste( the.objects,"2",sep="") # for testing fast agaist other against 
# loop over   "refGene::gene"   "knownGene::gene" "ensGene::gene"   get hits for all annotion al the indel locations
# i<-1
  
  for(i in 1:length(annotation.labels)){
 # print(i)
    #annotation.labels[i]
#    system.time( #40.5     0.0    40.5 
    counts.results<-get.count.data(the.table,the.genes,annotation.labels[i],the.global,the.filter,the.samples,the.attributes)
#                )
    #names(counts.results)

    list.element.lengths<-unlist(lapply(the.genes[[annotation.labels[i]]],length))        
    indel.index<-rep(1:length(list.element.lengths),times=list.element.lengths)
    flat.gene.list<-unlist(the.genes[[annotation.labels[i]]])

#    system.time(
    hit<-match.cols.and.collapse.fast(list.element.lengths,indel.index,flat.gene.list,counts.results[["phenotype.hit"]],"genes",c("genes","value"),"boolean")
#                )
#    system.time( # 0.18    0.00    0.18 
    hit.count<-match.cols.and.collapse.fast(list.element.lengths,indel.index,flat.gene.list,counts.results[["num.wanted.muts"]],"genes",c("genes","value"),"max")
#                )
#    system.time( #161.760  15.580 177.342 .. This one HIt/GENE/INDIVIDUAL takes the longest too calculate and could skip
    gene.hit.count<-match.cols.and.collapse.fast(list.element.lengths,indel.index,flat.gene.list,counts.results[["gene.counts"]],"genes",colnames(counts.results[["gene.counts"]]),"max")
#    )
## hit<-match.cols.and.collapse(the.genes,annotation.labels[i],counts.results[["phenotype.hit"]],"genes",c("genes","value"),"boolean")
## hit.count<-match.cols.and.collapse(the.genes,annotation.labels[i],counts.results[["num.wanted.muts"]],"genes",c("genes","value"),"max")
## gene.hit.count<-match.cols.and.collapse(the.genes,annotation.labels[i],counts.results[["gene.counts"]],"genes",colnames(counts.results[["gene.counts"]]),"max")

  assign(paste("hit",the.objects[i],sep="."),value=hit)
  assign(paste("hit.count",the.objects[i],sep="."),value=hit.count)
  assign(paste("gene.hit.count",the.objects[i],sep="."),value=gene.hit.count)
}

### not using match.cols.and.collapse() get an extra dimension
## sum(hit.refGene[,2]!=hit.refGene2)
## sum(hit.count.refGene[,2]!=hit.count.refGene2)
## sum(gene.hit.count.refGene[,4]!=gene.hit.count.refGene2[,3])


for(i in 1:length(the.objects)){
  data<-eval(as.name(paste("hit",the.objects[i],sep=".") ))
  if(i==1){hit<-data}else{hit<- hit | data}
}

  for(i in 1:length(the.objects)){
    data<-eval(as.name(paste("hit.count",the.objects[i],sep=".") ))
    if(i==1){hit.count<-data}else{hit.count<- pmax.int(hit.count,data)}
  }


for(i in 1:length(the.objects)){
  data<-eval(as.name(paste("gene.hit.count",the.objects[i],sep=".") ))
  data<-as.matrix(data) # in case just one sample
  if(i==1){the.dim<-dim(data);the.colnames<-colnames(data);data<-as.vector(data);gene.hit.count<-data}else{data<-as.vector(data);gene.hit.count<- pmax.int(gene.hit.count,data)}
}
dim(gene.hit.count)<-the.dim
colnames(gene.hit.count)<-the.samples

## hit.count<-apply(as.matrix(gene.hit.count),1,function(x) sum(x==0,na.rm=TRUE))
## hit.count<-length(the.samples)-hit.count
                 

## ann<-collapse.gene.names(the.genes,"refGene::gene",delimit.by=";")
## ann[hit]
} # end of else statment of there are no samples
  
list(hit=hit,hit.count=hit.count,gene.counts=gene.hit.count)
}
############################## end function ###########
## a.test<-"MAFB"
## the.test<-geneanno.table[,"refGene::gene"]==a.test

## geneanno.table[the.test,]
## hit.count[the.test]

## pheno.maf.qual<-CountOverGene.withFilterAndGlobal(input,the.affected,the.attributes,gene.names,filter.table,a.filter.cols.maf,"AND",filter.action="INCLUDE",global.vars=use.the.globals )

CountOverGene.withFilterAndGlobal.OLD<-function(the.table,the.samples,the.attributes,the.genes,the.filter.table,filter.cols,filter.cols.combine="OR",filter.action="EXCLUDE",global.vars=c() ){

  ## This function takes the quality threshold table made of samples and attributes filter with the filter table and applys & for the gloabl.vars
  ## Pheno hits are where the samples all have a filter hit in a gene
  ## the attributes are ones that would survive matrix mutliplication : ( quality.thresh[,dom.aff.GT] & quality.thresh[,dom.aff.GT.DP])
  ## filter.cols.combine="OR" or "AND"
  ## filter.action="EXCLUDE" or "INCLUDE
  ## global.vars<-c("wanted.muts.fil") ## names of columns vectors use & to combine
  ## include error if filter.cols not in the filter table
  ############ set any NA's to true
  the.NAs<-is.na(the.filter.table[,filter.cols])
  count.col.NAs<-apply(as.matrix(the.NAs),2,sum)# one filter .col gives a vector
  found.col.NAs<-count.col.NAs>0
  if(sum(found.col.NAs)>0){
    for(i in 1:length(found.col.NAs)){
      if(found.col.NAs[i]){
        the.filter.table[,filter.cols[i]][is.na(the.filter.table[,filter.cols[i]])]<-TRUE
      }}}
  ##########################
  
  
  num.trues<-apply(as.matrix(the.filter.table[,filter.cols]),1,function(x) sum(x,na.rm=TRUE))
  if(filter.cols.combine=="OR"){the.filter<-num.trues >0}
  if(filter.cols.combine=="AND"){the.filter<-num.trues==length(filter.cols)}

  if(filter.action=="EXCLUDE"){the.filter<-!the.filter}

  if(length(global.vars)>0){
    the.global<-eval(as.name(global.vars[1]))
    the.global[is.na(the.global)]<-TRUE
    if(length(global.vars)!=1){
      for(i in 2:(length(global.vars))){
        print(i)
        the.next<-eval(as.name(global.vars[i]))
        the.next[is.na(the.next)]<-TRUE
        the.global<-the.global & the.next
        }
  }} ## combine all the global variables
    
  
## the.table<-( quality.thresh[,dom.aff.GT] & quality.thresh[,dom.aff.GT.DP])
  
pheno.aff.gene.counts<-apply(as.matrix(the.table),2,function(x){
    x<- x & the.global & the.filter  ## column wise set not interesting 
    tapply(x,the.genes,sum) ## count the hits
})  ### this is the number that pass for each sample-attribute

  
## cbind(the.table,the.global,the.filter,the.genes)
##   length(the.attributes)
  
gene.counts<-data.frame(key=rownames(pheno.aff.gene.counts),stringsAsFactors=FALSE)
for(i in 1:length(the.samples)){
one.sample<-expand.labels.to.samples(the.attributes,the.samples[i])
one.count<-apply(pheno.aff.gene.counts[,one.sample],1,sum)
posns<-match(gene.counts[,"key"],names(one.count))
missing<- is.na(posns)
gene.counts[!missing,the.samples[i]]<-one.count[posns[!missing]]
} #cum over attributes here so have number that pass per gene

## warning gene.counts can go from table to vector if only one column is left losing the rownames
rownames.gene.counts<-gene.counts[,"key"]
colnames.gene.counts<-colnames(gene.counts)[colnames(gene.counts) !="key"]
gene.counts<-gene.counts[,colnames(gene.counts) != "key"] 
gene.counts<-as.matrix(gene.counts)
rownames(gene.counts)<-rownames.gene.counts 
colnames(gene.counts)<-colnames.gene.counts

num.wanted.muts<-apply(gene.counts,1,sum) # total filtered hits per gene per sample
phenotype.hit<-apply(gene.counts,1,function(x){
  x<-x>0
  sum(x)==length(the.samples)
}) # pheno hit if there is a hit in that gene for all samples

#####################expand to full format
posns<-match(the.genes,names(phenotype.hit))
missing<-is.na(posns)
hit<-phenotype.hit[posns[!missing]]  ### marker for is a phenotype hit

posns<-match(the.genes,names(num.wanted.muts))
missing<-is.na(posns)
hit.count<-num.wanted.muts[posns[!missing]] ### number of hits in gene

posns<-match(the.genes,rownames(gene.counts))
missing<-is.na(posns)
sum(missing)
rownames.gene.counts<-rownames(gene.counts)
colnames.gene.counts<-colnames(gene.counts)
gene.counts<-gene.counts[posns[!missing],] ### number of hits in gene per affected individual again can be a problem if only one affected


if(is.null(dim(gene.counts))){
  rownames.gene.counts<-rownames.gene.counts[posns[!missing]]
  gene.counts<-as.matrix(gene.counts)
  rownames(gene.counts)<-rownames.gene.counts
  colnames(gene.counts)<-colnames.gene.counts
}
  
list(hit=hit,hit.count=hit.count,gene.counts=gene.counts)
}
############################## end functions ###########


################################testing junk below;

## list.lengths<-list.element.lengths
## index<-indel.index
## flat.list<-flat.gene.list
## tab2<-gene.desc
## tab2.match.col<-"ensembl_gene_id"
##  tab2.collapse.cols<-gene.ann.wanted
##  gene.collapse="delimit"
## (the.genes,"refGene::gene",gene.counts,"genes",colnames(gene.counts),gene.collapse="sum")
                       ##   list1<-the.genes
                       ##   list1.match.col<-annotation.labels[i]
                       ##  tab2<-counts.results[["gene.counts"]]
                       ##   tab2.match.col<-"genes"
                       ##   tab2.collapse.cols<-colnames(counts.results[["gene.counts"]])
                       ## gene.collapse<-"max"

                       ## ##   list1<-the.genes
                       ## ##   list1.match.col<-"refGene::gene"
                       ## ##  tab2<-num.wanted.muts
                       ## ##   tab2.match.col<-"genes"
                       ## ##   tab2.collapse.cols<-c("genes","value")
                       ## ## gene.collapse<-"max"

                        ##  list1<-gene.names
                        ##  list1.match.col<-"refGene::gene"
                        ## tab2<-counts.results[["phenotype.hit"]]
                        ##  tab2.match.col<-"genes"
                        ##  tab2.collapse.cols<-c("genes","value")
                        ## boolean.collapse<-TRUE


                        ## list1<-gene.names
                        ## list1.match.col<-"ensGene::gene"
                        ## tab2<-gene.desc
                        ## tab2.match.col<-"ensembl_gene_id"
                        ## tab2.collapse.cols<-gene.ann.wanted
                        ##  gene.collapse="delimit"

##                         list1<-omim.names
##                         list1.match.col<-1
##                         tab2<-omim.desc
##                         tab2.match.col<-"omim::name"
##                         tab2.collapse.cols<-colnames(omim.desc)
## length(tab2[posns[!missing],tab2.collapse.cols[i]])
##   length(gene.desc.table[!poly.genes,tab2.collapse.cols[i]][!missing]) # if only one dimension
## ## length( as.matrix(gene.desc.table[!poly.genes,])[!missing,tab2.collapse.cols[i]] )
## ## ##  length(a.desc.table[poly.genes,][,tab2.collapse.cols[i]])
## ## ## length(other.genes[,tab2.collapse.cols[i]])
## length(other.genes[,tab2.collapse.cols[1]])
## length(gene.desc.table[poly.genes,tab2.collapse.cols[1]])

##   class(gene.desc.table[!poly.genes,tab2.collapse.cols[i]])
##          [!missing,tab2.collapse.cols[i]])

  
##  dim(gene.desc.table[!poly.genes,])
## gene.desc.table@value
## [!poly.genes][1:5]  [!missing,tab2.collapse.cols[i]])


###http://www.biomedcentral.com/1471-2105/9/309#sec2
## perhaps MANTEL trend test better because it is not affected by departure fom hardy weinburg http://www.jstor.org/stable/2282717?origin=crossref&
## see library (ncf)  ade4 library cluster libraray __ I think I need all the genotypes to run 


 CATTs.row<-function ( controls,cases, scores = NULL, add.pval = TRUE)
    ## Cochran–Armitage test for trend
   ### modified rowCATTS from {scrime} it's statistics as stat^2 , now sign is correct
  ##3 swapped order of cases and controls so consistent with othe CATT defined below
  ## c(0,1,2)-additive c(0,1,1)-dominant  c(1,1,0)-recessive 
{
    if (!is.matrix(cases) | !is.matrix(controls)) 
        stop("cases and controls must be matrices.")
    if (any(dim(cases) != dim(controls))) 
        stop("cases and controls have not the same dimensions.")
    rn <- rownames(cases)
    if (any(rn != rownames(controls))) 
        stop("The row names differ between cases and controls.")
    cn <- colnames(cases)
    if (any(colnames(controls) != cn)) 
        stop("The column names differ between cases and controls.")
    if ("NA" %in% cn) {
        ids <- which(cn == "NA")
        cases <- cases[, -ids]
        controls <- controls[, -ids]
        warning("The column named NA is removed from both cases and controls.")
    }
    if (is.null(scores)) {
        if (is.null(cn)) 
            stop("Either scores must be specified or the matrices must have\n", 
                "(the same) column names specifying numeric scores.")
        scores <- as.numeric(cn)
        if (any(is.na(scores))) 
            stop("At least one of column names does not specify a numeric score.")
    }
    else {
        if (any(!is.numeric(scores))) 
            stop("scores must be numeric.")
        if (length(scores) != ncol(cases)) 
            stop("The number of scores must be equal to the number of columns.")
    }
    if (any(cases < 0)) 
        stop("All values in cases must be non-negative integers.")
    if (any(controls < 0)) 
        stop("All values in controls must be non-negative integers.")
    n.cases <- rowSums(cases)
    n.controls <- rowSums(controls)
    n <- n.cases + n.controls
    ybar <- n.cases/n
    xbar <- as.vector((cases + controls) %*% scores)
    xbar <- xbar/n
    mat.scores <- matrix(scores, length(n), length(scores), byrow = TRUE)
    mat.scores <- mat.scores - xbar
    cases <- cases * mat.scores
    controls <- controls * mat.scores
    num <- (-ybar) * rowSums(cases) + (1 - ybar) * rowSums(controls)
    cases <- cases + controls
    denom <- rowSums(cases * mat.scores)
    denom <- denom * ybar * n.controls
    stats <- num * num/denom
    stats <- stats * n
    if (!add.pval) 
        return(sqrt(stats)*sign(num))
    rawp <- pchisq(stats, 1, lower.tail = FALSE)
     structure(data.frame(stats =sqrt(stats)*sign(num), p.value = rawp),stringsAsFactors=FALSE)
#    structure(list(stats =sqrt(stats)*sign(num), rawp = rawp))
}

catt.rows<-function(data, score){
## Cochran–Armitage test for trend cbind(cases,controls)
## data c( controls,cases )   cols controls/cases are  aa,aA,aA A=effect allele
        data=matrix(data,nrow=2,byrow=TRUE)
	nn=apply(data,2,sum)
	n=sum(nn)

	Rbar=sum(nn*score)/n
	s2=sum(nn*(score-Rbar)^2)
	phi=sum(data[1,])/n

	catt=sum(data[1,]*(score-Rbar))/sqrt(phi*(1-phi)*s2)
	p.value = 2*(1-pnorm(abs(catt)))

	out=c(estimate = catt,p.value = p.value)
#	out=list(estimate = catt,p.value = p.value)
	return(out) 

} 

catt=function(data, score){
## Cochran–Armitage test for trend best for 2x3 table
## data 2x3 rows controls,cases  cols aa,aA,aA A=effect allele
	nn=apply(data,2,sum)
	n=sum(nn)

	Rbar=sum(nn*score)/n
	s2=sum(nn*(score-Rbar)^2)
	phi=sum(data[1,])/n

	catt=sum(data[1,]*(score-Rbar))/sqrt(phi*(1-phi)*s2)
	p.value = 2*(1-pnorm(abs(catt)))


	out=list(estimate = catt,p.value = p.value)
	return(out) 

} # end ftn_catt


## ########## TEST GERP
## library("DBI")
## library("RMySQL")
## db.details<-list(user="gerp2",pass="GeRpUQ8!6", dbname="gerp", host="ga-apps.di.uq.edu.au")
## drv <- dbDriver("MySQL") 
## con <- dbConnect(drv,user=db.details$user,pass=db.details$pass, dbname=db.details$dbname, host=db.details$host)
## dbname<-db.details$dbname
## test.mysql <- dbGetQuery(con,"select * from gerp where chr = \"20\" and (position > \"39316700\" and position < \"39316900\") ")
## test.mysql[1:50,]
## query<-paste("select position,score from ","gerp"," where chr = \"","20","\" and  ( ", " ( position >= \"", "39316730","\" and position <= \"",  "39316730" ,"\" )" , " )", sep="")
## query<-paste("select position,score from ","gerp"," where chr = \"","1","\" and  ( ", " ( position >= \"", "248785873","\" and position <= \"",  "248785873" ,"\" )" , " )", sep="")
## dbGetQuery(con,query)
## #############

## indels<-indels[,c("chr","start","end","REF")]
get.GERP.MULT<-function(db.details,indels,a.chunk=500){
  # gets gerp by connecting to a SQldb
  # a.chunk=500 seems fastest : trade off between qeury complexity and number of queries
  #db.details<-list(user="gerp2",pass="GeRpUQ8!6", dbname="gerp", host="ga-apps.di.uq.edu.au")
  # indels needs chr, start, end and REF (REF="-" the is deletion
  # 695017 vars  87.160    0.940 1981.781 for a.chunk=500 (33mins)
 print( dim(indels))
 print(a.chunk)
use.chrs<-sort(tapply(indels[,"chr"],indels[,"chr"],length))
if(grepl("^chr",names(use.chrs)[1])){names(use.chrs)<-gsub("^chr","",names(use.chrs)); indels[,"chr"]<-gsub("^chr","",indels[,"chr"])} ## gerp databse does not use chr



##### launch a parrallel process for each chromomes and collect data afterwards
result<-list()
 # ic<-1
for(ic in 1:length(use.chrs)){
  use.chr<-names(use.chrs)[ic]
  to.get<-indels[,"chr"]==use.chr 
  ## sum(to.get)
  print(paste("Starting GERP annotation for Chr",use.chr))
#chr start      end        strand ID            REF     ALT  
#"Y" "9307429"  "9307429"  "+"    "."           "-"     "GTGT
  the.starts<-as.numeric(indels[to.get,"start"])
  the.ends<-as.numeric(indels[to.get,"end"])
  a.insert<-indels[to.get,"REF"]=="-"

  result[[ic]]<-parallel(get.GERP.for.chromsome.MULTI(db.details,use.chr,the.starts,the.ends,a.insert,a.chunk))
}

print("Collecting Data...")
gerp.data<-rep(NA,times=dim(indels)[1])
#names(result)<-use.chrs

for(ic in 1:length(use.chrs)){
  use.chr<-names(use.chrs)[ic]
  to.get<-indels[,"chr"]==use.chr 
  chr.gerp.data<-collect(result[[ic]],wait=TRUE)
  gerp.data[to.get]<-chr.gerp.data[[1]]
}  # collect data loop
print("Got GERP!")
gerp.data<-as.numeric(gerp.data)
gerp.data[!is.finite(gerp.data)]<-0 # NA for snps -Inf for indels at end of chrosomes where GERP=0| is.na(gerp.data)
names(gerp.data)<-build.key(indels,c("chr","start","end","REF"))  
gerp.data

} # get function get gerp




#get.GERP.for.chromsome.MULTI(db.details,use.chr,the.starts,the.ends,a.insert,a.chunk

get.GERP.for.chromsome.MULTI<-function(db.details,use.chr,the.starts,the.ends,a.insert,a.chunk){

#### for insert also include the neighbouring bases
drv <- dbDriver("MySQL") 
con <- dbConnect(drv,user=db.details$user,pass=db.details$pass, dbname=db.details$dbname, host=db.details$host)
dbname<-db.details$dbname

to.get<-length(the.starts)
if(length(the.starts)==0){print("Error in get.GERP.for.chromsome.MULTI no data to get" )}

chr.gerp.data<-rep(NA,times=to.get)  
sizes<-the.ends-the.starts  
the.starts[a.insert]<-the.starts[a.insert]-1
the.ends[a.insert]<-the.ends[a.insert]+1

a.snp<-sizes==0 & !a.insert
## sum(a.snp)
## sum(!a.snp)

##a.chunk sweet spot:
##5000 steps 0.770   0.070  89.935
##1000  3.260   0.070  13.413 
##500  5.490   0.030  13.273
##256  11.860   0.200  19.144 
##100 30.410   0.470  40.035
the.seq<-seq(from=0,to=to.get,by=a.chunk)
the.seq<-unique(c(the.seq,to.get))

## ichunk<-100

 for(ichunk in 1:(length(the.seq)-1) ){
    the.start<-the.seq[ichunk]+1
    the.end<-the.seq[ichunk+1]

##     the.start
##     the.end
if(ichunk%%(a.chunk*10)==0){print(paste("iter:",ichunk," region:",the.start,the.end))}
####
    data.region<-the.start:the.end
    starts.chunk<-the.starts[data.region]
    ends.chunk<-the.ends[data.region]
    a.snp.chunk<-a.snp[data.region]



extras<-paste( " ( position >= \"", starts.chunk,"\" and position <= \"", ends.chunk,"\" ) ",sep="")
extras<-paste(extras,collapse=" or ")
 
## query<-paste("select position,score from ",dbname," where chr = \"",use.chr,"\" and  ( ", extras, " )", sep="")
query<-paste("select Distinct  position,score from ",dbname," where chr = \"",use.chr,"\" and  ( ", extras, " )", sep="")    
## print("start query")
# query<-paste("select position,score from ","gerp"," where chr = \"","1","\" and  ( ", " ( position >= \"", "249240334","\" and position <= \"",  "249240334" ,"\" )" , " )", sep="")
# query<-paste("select position,score from ","gerp"," where chr = \"","1","\" and  ( ", " ( position >= \"", "248785873","\" and position <= \"",  "248785873" ,"\" )" , " )", sep="")
  
gerp.chr <- dbGetQuery(con,query)
# print("end query")  
## dim(gerp.chr)
posns<-match(starts.chunk[a.snp.chunk],gerp.chr[,"position"])

## Below only needed for for diagnosic error correction
## missing<-is.na(posns)
## if(sum(missing)>0){print(paste("Error",ic,ichunk,sum(missing)));print(gerp.chr[1:100,]);save(list=c("gerp.chr","starts.chunk","a.snp.chunk","posns","missing","query","the.start","the.end","ichunk"),file="gerp.chr.RData")}  # load("gerp.chr.RData")
################################
    
chr.gerp.data[data.region[a.snp.chunk]]<-gerp.chr[posns,"score"]

    if(sum(!a.snp.chunk)>0){
    
ir.gerp<-IRanges(start=gerp.chr[,"position"],width=1)
cov<-coverage(ir.gerp,weight=as.numeric(gerp.chr[,"score"])) ## Problems exist here is the ir.gerp are not unique get coverage*weight 
myViews <- Views(cov,start=starts.chunk[!a.snp.chunk],end=ends.chunk[!a.snp.chunk] )
chr.gerp.data[data.region[!a.snp.chunk]]<-viewMaxs(myViews,na.rm=TRUE) # get "-Inf" in region where there is no data

## The below does the same as the above 4 lines does not requre SQL to be unique
## ir<-IRanges(start=starts.chunk[!a.snp.chunk],end=ends.chunk[!a.snp.chunk])
## ir.gerp<-IRanges(start=gerp.chr[,"position"],width=1)
## subject<-findOverlaps(ir,ir.gerp)
## query<-queryHits(subject)
## subject<-subjectHits(subject)
## ## must do below else can get cases at end of chr when there is no overlap.
## the.max<-tapply(as.numeric(gerp.chr[subject,"score"]),query,function(x) max(x,na.rm=TRUE))
## posns<-match(as.character(c(1:sum(!a.snp.chunk))),names(the.max))
## chr.gerp.data[data.region[!a.snp.chunk]]<-the.max[posns]
###############################################################
                       } # no indels to process
######## make an RleList for a score use "select Unique" in SQL : COVERAGE AND WEIGHTS
    

  } # loop over chunks
## print(paste("done chr:",ic))
dbDisconnect(con)
chr.gerp.data

} # end function get.GERP.for.chromsome.MULTI




get.data<-function(x,con){
unlist(dbGetQuery(con,x))
}



get.data.multi<-function(x){
con <- dbConnect(drv,user="gerp2",pass="GeRpUQ8!6", dbname="gerp", host="ga-apps.di.uq.edu.au")
unlist(dbGetQuery(con,x))
dbDisconnect(con)
}

