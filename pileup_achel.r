library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(XLConnect)
library("openxlsx")
library("data.table")
library("dplyr")
library("plyr")

code.dir <- "/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/Pvalue_enrichment_july9/"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")
core.ann <- c( "chr","start","end","ALT","REF")
core.ann.gr <- c( "chr","start","end","strand")
#mydf  <-  read.delim("mydf_enrichment_test_file.txt",header=TRUE,skip=0,sep="\t",check.names=FALSE,fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="")
mydf  <-  read.delim("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/Pvalue_enrichment_july9/2015-03-25_SKDPNexteraRun.THP.wanted.All-maf-filtered.txt",header=TRUE,skip=0,sep="\t",check.names=FALSE,fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="")
finaldf<-mydf  # to cbind the final coverage output

#####Check if chr is missing 
if (grepl("^chr", mydf[,"chr"])!=TRUE){
 mydf[,"chr"]<- paste("chr",mydf[,"chr"], sep="")
}

#cc<-mydf[!grepl("flat$",mydf[,"TYPE"]),] ##omitting all the flat ones
#mydf<-cc
bam.file <- "THP-1-AS703026resist.ReCal.sort.bam"
#bam.file <- "THP-1-Selumetinibresist.ReCal.sort.bam"
#bam.file <- "5_D1K4YACXX-5-01.ReCal.sort.bam"
#bam.dir <- "/media/UQCCG/Sequencing/Projects/AMAS/BAM/"

bam.dir <- "/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Projects/AMAS/BAM"

#bam.dir <- "/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Projects/TGCM-AML/BAM/"

#Selection regiions to check



##########################

#allkeys  <-  cat(paste(shQuote(all, type="cmd"), collapse=", "))
#dput(as.character(all))
#write.table(matrix(as.character(all),nrow=1), sep="'",row.names=FALSE, col.names=FALSE)
##########################
#coverage.bam  <-  read.delim("mycoverage.bam_enrichment_test_file.txt",header=TRUE,skip=0,sep="\t",check.names=FALSE,fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="")


#summary <- c("chr14","106066660","106066660","+","chr7","100680296","100680296","+")
#summary <- c("chr14","106066660","106066660","+","chr7","100680296","100680296","+")
#summary <- c("chr14","106066660","106066660","+")
#summary <- c("chr10","103310574","103310574","+","chr3","75790814","75790814","+")
#summary <- matrix(data=summary,nrow=2,ncol=4,byrow=TRUE)
#summary
all.keys <- cbind(mydf[,"chr"],mydf[,"start"],mydf[,"end"],mydf[,"strand"])
all.keys <- as.matrix(all.keys)
summary <- all.keys[c(4,6),]
summary <- all.keys
############################################ This concatenates all the elements in the column (don't need to use in here)
#all  <-  paste(apply(all.keys, 1, function(row) paste(dQuote(row), collapse=",")), collapse=',"+",')
#all.keys <- as.matrix(all)
#summary <- all.keys
#summary <- matrix(data=summary,nrow=813,ncol=4,byrow=TRUE)
#paste(t(all.keys), collapse=', ')
#paste(dQuote(t(df1)), collapse=', ')
############################################################

#summary <- summary[14:15,]

##      chr     start       end         strand
## [1,] "chr14" "106066660" "106066660" "+"   
## [2,] "chr7"  "100680296" "100680296" "+"


colnames(summary) <- core.ann.gr
summary #### a matrix that has position you want to look at
options(width=150) ## output length

data <- summary ### summary could be a subset of a.indel - . a.indel[wanted,core.ann.gr)
#'elementMetadata' cannot use "seqnames", "ranges", "strand", "seqlengths", "start", "end", "width", or "element" as column names
data.gr <- GRanges(seqnames =data[,"chr"],ranges = IRanges(start=as.numeric(data[,"start"]),end=as.numeric(data[,"end"])),strand=data[,"strand"])

## reserved.gr.names <- c("seqnames", "ranges", "strand", "seqlengths", "isCircular", "start", "end", "width", "element")
## names.have <- colnames(data)[!(colnames(data) %in% core.ann.gr)]
## names.have[names.have %in% reserved.gr.names] <- paste(names.have[names.have %in% reserved.gr.names],1,sep=".")
## colnames(data)[!(colnames(data) %in% core.ann.gr)] <- names.have
## values(data.gr) <- data[,(!(colnames(data) %in% core.ann.gr)) & !remove.cols.from.values ]

which <-   data.gr
which
#51623174

setwd(bam.dir)

files <- dir(getwd())
files[grepl("bam",files)]
files[match(bam.file,files)]  ### check  bam file must exist in specified locations
?ScanBamParam

params <- ScanBamParam(which=which,flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=FALSE,isNotPassingQualityControls=FALSE),simpleCigar = FALSE,reverseComplement = FALSE,what=c("qname","flag","rname","seq","strand","pos","qwidth","cigar","qual","mapq") )  ### NOTE isValidVendorRead=FALSE shoudl be TRUE


param.pile <- PileupParam(max_depth=2500, min_base_quality=15, min_mapq=0,min_nucleotide_depth=1, min_minor_allele_depth=0,distinguish_strands=TRUE, distinguish_nucleotides=TRUE,ignore_query_Ns=TRUE, include_deletions=TRUE,cycle_bins=numeric() )

test <- pileup(bam.file,scanBamParam=params,pileupParam=param.pile,ApplyPileupsParam=param)
setwd("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/Pvalue_enrichment_july9/")

#test[1:20,]
#dups<-duplicated(test[,"pos"])
#test[dups,]
#unique(test[,"which_label"])
#########################using tapply###############
#test.keys <- paste(test[,"seqnames"],test[,"pos"],sep=":")
newtest <- within(test, start  <-  paste(seqnames,pos, sep=':'))
newtest[1:5,]
#length(unique(newtest[,"start"]))
#crap<-strsplit(as.character(newtest[,"which_label"]),split="-")
#crap[1:5]
#crap<-unlist(lapply(crap,function(x) x[1]))
#sort(tapply(crap,crap,length),decreasing=TRUE)[1:10]

#newtest[newtest[,"start"] %in% names(sort(tapply(crap,crap,length),decreasing=TRUE)[1:2]),]

#sum(newtest[,"start"]!=crap)
#newtest[newtest[,"start"]!=crap,][1:10,]

yyy <- tapply(newtest[,"count"],list(newtest[,"nucleotide"],newtest[,"which_label"]),sum)
dim(yyy)
yyy[1:5,]

#yyy <- tapply(newtest[,"count"],list(newtest[,"nucleotide"],newtest[,"start"]),sum)
coverage.bam <- t(yyy)
coverage.bam <- cbind(rownames(coverage.bam),coverage.bam)
colnames(coverage.bam)[1] <- c("start")
rownames(coverage.bam) <- NULL
#coverage.bam<-(coverage.bam,c(gsub("*-.*", "\\1", coverage.bam[,"which_label"]),"newcol"))
coverage.bam[,"start"] <-c(gsub("*-.*", "\\1", coverage.bam[,"start"]))

#Since the bam file is missing some of the ncRNA_intronic, UTR, need to remove the difference between the file from the coverage.bam matrix  
#diff1 <- setdiff(mydf[,"start"],coverage.bam[,"start"])
mydf <- within(mydf, start  <-  paste(chr,start, sep=':'))
#No need to reduce the size
#mydf <- mydf[, c("chr","start", "end","REF","ALT","TYPE","refGene::location","refGene::type")]

###########

dim(mydf)  #6179,268
dim(coverage.bam) # 6017,8
missing <- !mydf[,"start"]%in%coverage.bam[,"start"]
sum(missing) # 148(excluding the duplicates)
#duplicated.poss <- duplicated(mydf[,"start"])
#sum(duplicated.poss) 
#missing.rows <- mydf[missing,]
#rows.present <- mydf[(mydf[,"start"]%in%coverage.bam[,"start"]),] #matching rows including duplicates
#mydf[,[match(mydf[,"start"],coverage.bam[,"start"])]]
#sum(grepl("^snp+:[0-9]", mydf[,"TYPE"]))
#wanted.snp <- mydf[(grepl("^snp+:[0-9]", mydf[,"TYPE"])),] # extracts all the rows that are 'snp:~'
#wanted.snp <- wanted.snp[!grepl("flat$",wanted.snp[,"TYPE"]),] ##omitting all the flat ones
#dupli.wanted <- duplicated(wanted.snp[,"start"])
#dupli.wanted <- wanted.snp[dupli.wanted,]
#one.set.dupli <- dupli.wanted
#sum(duplicated(wanted.snp[,"start"]))
#match(wanted.snp[,"start"],dupli.wanted[,"start"])
#dupli.wanted <- wanted.snp[duplicated.poss,]
#all.duplicates <- wanted.snp[wanted.snp[,"start"] %in% dupli.wanted[,"start"],]
#ttt <- wanted.snp[match(wanted.snp[,"start"],dupli.wanted[,"start"]),]
#dupli.start <- dupli.wanted[,"start"]
#dupli.rows <- mydf[(mydf[,"start"] %in% c(dupli.start)),] #the duplicate set of rows, 14 X 2 rows
#uniq.dupli.rows <- dupli.rows[!duplicated(dupli.rows$start),] 
#uniq.dupli.rows[, 'ALT']  <-  aggregate(ALT~start, data=dupli.rows, toString)[,2] ##################### #1   duplicates with aggregated ALTs
#sum(( mydf[,"start"] %in% c(dupli.start)))
rows.present <- mydf[(mydf[,"start"]%in%coverage.bam[,"start"]),] #matching rows including duplicates (i.e only -148, still 14X2 are duplicates)
dim(rows.present)

# rows.present.uniq <- !duplicated(rows.present[,"start"])
# rows.present.uniq <- rows.present[rows.present.uniq,]
# mydf.uniq<-rows.present.uniq
# dim(rows.present.uniq)

#matching.rows <- mydf[!missing&!duplicated.poss,] #these are the ones that are not duplicated by start(chr and positions) and also excludes the ones that are present in mydf but not in the bam.
######################
######################

#colnames(rows.present)
#foo  <-  aggregate(dupli.rows, by=list(Start=factor(dupli.rows$start)), FUN=function(x) x$ALT)
#concat.alt  <-  ddply(rows.present, .(ALT), summarise, alt=c(ALT))
#dim(concat.alt)
#head(concat.alt)
#foo  <-  ddply(rows.present, .(ALT), summarise, alt=ALT)
#head(concat.alt)
#concat.alt  <-  ddply(rows.present, .(start), summarise, alt=ALT)
#head(concat.alt)
#concat.alt  <-  ddply(rows.present, .(start), summarise, alt=c(ALT))
#head(concat.alt)
#concat.alt  <-  ddply(wanted.snp, .(start), summarise, alt=list(ALT))
#concat.alt  <-  ddply(rows.present, .(start), summarise, alt=paste(ALT, collapse=","))
#duplicated(rows.present$start)
#rows.present[!duplicated(rows.present$start), ]
#uniq.rows.present  <-  rows.present[!duplicated(rows.present$start), ]
#dim(uniq.rows.present)
#length(concat.alt)
#dim(concat.alt)
#concat.alt  <-  ddply(rows.present, .(start), summarise, ALT=paste(ALT, collapse=","))
#rows.concat.alt<-uniq.rows.present
#rows.concat.alt$ALT <- concat.alt$ALT[match(rows.concat.alt$start,concat.alt$start)]


#intersect.df <- intersect(mydf[,"start"],coverage.bam[,"start"])
dim(coverage.bam)
#dim(rows.concat.alt)
#rows.concat.alt <- rows.concat.alt[(!is.na(match(rows.concat.alt[,"start"] ,intersect.df))),]ix

#####################reordering according to the target vector#################
###Don't need to use this for now, because the coverage.bam is in the same order as the indel's input (GRranges)
#idx  <-  sapply(rows.present[,"start"], function(x) {
#  which(coverage.bam[,"start"] == x)
#})
#
#coverage.bam  <-  coverage.bam[unlist(idx),]
#coverage.bam

######################### difference between the dataframes############




#######################
#coverage.bam <- cbind(rows.concat.alt[,"chr"],rows.concat.alt[,"end"],coverage.bam)
#colnames(coverage.bam)[1:2] <- c("chr", "end")
coverage.bam <- coverage.bam[,c("start","A","C","G","T","N","=","-") ]
coverage.bam[1:5,]
#rows.concat.alt[1:5,]
###############split the contents of the ALT column at the delim ","   ##############
#split.col  <-  strsplit(rows.concat.alt[,"ALT"],",") ##not required!!
#split.col  <-  data.frame(do.call('rbind', strsplit(as.character(rows.concat.alt$ALT),',',fixed=TRUE)))
#split.col = transform(rows.concat.alt, ALT = colsplit(ALT, split = "\\,"))
#library(splitstackshape)
#split.col <- cSplit(rows.concat.alt, "ALT", ",")
#additional.cols <- rbind.fill.matrix(lapply(split.col, rbind))#not required
#colnames(additional.cols) <- c(paste("ALT",sequence(ncol(additional.cols)),sep=""))
#Not required!
#rows.concat.alt <- cbind(rows.concat.alt[,c("start","end","REF","ALT")],additional.cols, rows.concat.alt[, -which(names(rows.concat.alt) %in% c("chr","start","end","REF","ALT"))])

#result =merge(coverage.bam, rows.present, by="start") #merging bam data and mydf file
result<-cbind(coverage.bam,rows.present)
result<-result[!duplicated(result),]
cols <- names(result) == "start"
names(result)[cols] <- paste0("start", seq.int(sum(cols)))
#check if they are same
colnames(result)[which(names(result) == "start1")]<-"start"

sum(result[,"start"] %in% result[,"start2"])
wanted.pos<-!grepl("^indel",result[,"TYPE"])
sum(wanted.pos)
wanted.result<-result[wanted.pos,]
#newdf<-bamAD(wanted.result)
wanted.result[,"start2"]<- NULL
want.result<-wanted.result
new.file<-bamAD(want.result)
ggg<-new.file[,c(1:14,266)]
ggg[ggg[,"start"]=="chr10:17659149",]
which(colnames(result)=="THP-1-AS703026resist.AD")

########################################Function in annotate_SNPs_subroutine.r
normalCase <- function(x, ns) {
  ref.idx <- which(ns == "REF")
  ref.allele <- x[ref.idx]
  ref.count <- x[which(ns == ref.allele)]
  
  alt.idx <- which(ns == "ALT")
  alt.allele <- x[alt.idx]
  alt.count <- x[which(ns == alt.allele)]
  
  paste(ref.count, alt.count, sep=",")
}

flatCase <- function(x, mat, ns) {
  id.idx <- which(ns == 'start')
  type.idx <- which(ns == 'TYPE')
  ref.idx <- which(ns == 'REF')
  alt.idx <- which(ns == 'ALT')
  
  
  id <- x[id.idx]
  #m <- mat[mat[, id.idx] == id & mat[, type.idx] == "snp", ]
  #m <- mat[mat[, id.idx] == id & mat[, type.idx] == "snp", ]
  m<-mat[grepl(id,mat[, id.idx]) & grepl("snp:+[0-9]",mat[, type.idx]),]
  #flat<-mat[grepl("flat$",mat[, type.idx]),]
  ref.allele <- x[ref.idx]
  ref.count<-x[which(ns == ref.allele)]
  
  #indx <- grepl('flat', mat[,'TYPE'])
  #indx1 <- grepl('snp:\\d+', mat[,'TYPE'])
  #alt.count <- as.numeric(sub('.*,\\s+', '', v2[indx1]))
  #alt.allele <- x[alt.idx]
  alt.count <- sum(apply(m, 1, function(x) as.numeric(x[which(ns == x[alt.idx])])))
  paste(ref.count, alt.count, sep=",") 
}

calculateAD <- function(x, mat, ns) {
  if (grepl("flat$", x[which(ns == 'TYPE')])) {
    flatCase(x, mat, ns)
  } else {
    normalCase(x, ns)
  }
}

bamAD <- function(x) {
   new.x <- cbind(x, apply(x, 1, calculateAD, x, colnames(x)))
  colnames(new.x)[ncol(new.x)] <- "bam.AD"
  new.x      
}


x <- as.matrix(read.csv(text="start,A,T,G,C,REF,ALT,TYPE 
chr20:5363934,95,29,14,59,C,T,snp
chr5:8529759,24,1,28,41,G,C,snp
chr14:9620689,65,49,41,96,T,G,snp
chr18:547375,94,1,51,67,G,C,snp
chr8:5952145,27,80,25,96,T,T,snp
chr14:8694382,68,94,26,30,A,A,snp
chr16:2530921,49,15,79,72,A,T,snp:2530921
chr16:2530921,49,15,79,72,A,G,snp:2530921
chr16:2530921,49,15,79,72,A,T,snp:flat
chr14:4214117,73,49,18,77,G,A,snp
chr4:7799768,36,28,1,16,C,A,snp
chr3:9141263,27,41,93,90,A,A,snp", stringsAsFactors=FALSE))

