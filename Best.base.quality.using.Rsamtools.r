

## /media/UQCCG/Sequencing/Projects/AMAS/BAM/THP-1-AS703026resist.ReCal.sort.bam
## /media/UQCCG/Sequencing/Projects/AMAS/BAM/THP-1-normal-parent.ReCal.sort.bam
## /media/UQCCG/Sequencing/Projects/AMAS/BAM/THP-1-Selumetinibresist.ReCal.sort.bam
## THP-1-normal-parent.GT	THP-1-Selumetinibresist.GT	THP-1-AS703026resist.GT	THP-1-AS703026resist	THP-1-Selumetinibresist	ALT.reads	REF.reads	Read.Balance	ALT.reads.thresh	REF.reads.thresh	Read.Balance.thresh	Num.pos.thresh	Num.samples.thresh	THP-1-normal-parent.AD	THP-1-Selumetinibresist.AD	THP-1-AS703026resist.AD

## 0/0	0/1	0/0	2.60606541800293E-007	0.0003192576	1	29	3.33	0	0	0	0	0	29,1	53,7	40,8
## 0/0	0/1	0/0	5.47409567085124E-007	3.0251336455229E-033	0	25	0	0	0	0	0	0	25,0	25,2	36,1
## 1/1	1/1	1/1	0	0	0	0	0	0	0	0	0	0	0,2	0,11	0,4
## 1/1	NA	1/1	2.96230628546564E-219	1.0252397394	0	0	0	0	0	0	0	0	0,2	NA	0,1
## 1/1	1/1	1/1	0	0	0	0	0	0	0	0	0	0	0,53	0,97	0,84
## 0/0	0/1	0/0	3.88628235198148E-006	2.53671370192147E-082	0	25	0	0	0	0	0	0	25,0	21,3	42,1


## chr14	106066660	106066660	C	T
## chr7	100680296	100680296	A	C
## chr10	123845322	123845322	T	C
## chr10	124221227	124221227	C	T
## chr10	124271589	124271589	G	A
## chr8	101730326	101730326	G	C
## /media/UQCCG/Sequencing/Projects/AMAS/BAM/THP-1-AS703026resist.ReCal.sort.bam
## /media/UQCCG/Sequencing/Projects/AMAS/BAM/THP-1-normal-parent.ReCal.sort.bam
## /media/UQCCG/Sequencing/Projects/AMAS/BAM/THP-1-Selumetinibresist.ReCal.sort.bam


#########Define BAM files
library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
core.ann.gr<-c( "chr","start","end","strand")


bam.file<-"THP-1-AS703026resist.ReCal.sort.bam"
bam.dir<-"/media/UQCCG/Sequencing/Projects/AMAS/BAM/"


#Selection regiions to check
summary<-c("chr14","106066660","106066660","+","chr7","100680296","100680296","+")
summary<-matrix(data=summary,nrow=2,ncol=4,byrow=TRUE)
summary
##      chr     start       end         strand
## [1,] "chr14" "106066660" "106066660" "+"   
## [2,] "chr7"  "100680296" "100680296" "+"


colnames(summary)<-core.ann.gr
summary #### a matrix that has position you want to look at
options(width=150) ## output length

data<-summary ### summary could be a subset of a.indel - . a.indel[wanted,core.ann.gr)
 #'elementMetadata' cannot use "seqnames", "ranges", "strand", "seqlengths", "start", "end", "width", or "element" as column names
data.gr<-GRanges(seqnames =data[,"chr"],ranges = IRanges(start=as.numeric(data[,"start"]),end=as.numeric(data[,"end"])),strand=data[,"strand"])

reserved.gr.names<-c("seqnames", "ranges", "strand", "seqlengths", "isCircular", "start", "end", "width", "element")
names.have<-colnames(data)[!(colnames(data) %in% core.ann.gr)]
names.have[names.have %in% reserved.gr.names]<-paste(names.have[names.have %in% reserved.gr.names],1,sep=".")
colnames(data)[!(colnames(data) %in% core.ann.gr)]<-names.have
values(data.gr)<-data[,(!(colnames(data) %in% core.ann.gr)) & !remove.cols.from.values ]

which<-  data.gr
which


setwd(bam.dir)

files<-dir(getwd())
files[grepl("bam",files)]
files[match(bam.file,files)]  ### check  bam file must exist in specified locations


params<-ScanBamParam(which=which,flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=NA,isNotPassingQualityControls=FALSE),simpleCigar = FALSE,reverseComplement = FALSE,what=c("qname","flag","rname","seq","strand","pos","qwidth","cigar","qual","mapq") )  ### NOTE isValidVendorRead=FALSE shoudl be TRUE

## params<-ScanBamParam(which=which,flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=NA,isValidVendorRead=FALSE),simpleCigar = FALSE,reverseComplement = FALSE,what=c("seq","strand","pos","cigar") )
params
which
## params<-ScanBamParam(which=which,flag=scanBamFlag(isPaired=TRUE,isUnmappedQuery=FALSE,hasUnmappedMate=FALSE,isFirstMateRead=TRUE),what=c("qname","flag","strand","pos","qwidth","cigar","mrnm","mpos","isize") )

param.pile<-PileupParam(max_depth=2500, min_base_quality=10, min_mapq=13,
min_nucleotide_depth=1, min_minor_allele_depth=0,
distinguish_strands=TRUE, distinguish_nucleotides=TRUE,
ignore_query_Ns=TRUE, include_deletions=TRUE, include_insertions=FALSE,
cycle_bins=numeric() )

test<-pileup(bam.file,scanBamParam=params,pileupParam=param.pile)

######################## Phred scores DECODED are ascii 0-93 ### THIS CAN CAHNGES IN THE HISE 4000 or for a SOLID data 
str<-as.character(PhredQuality(0:93))
str<-unlist(strsplit(str,split = character(0)))
score<-0:93
names(score)<-str
 score
## the.qual<-unlist(strsplit(as.character(aln1[[1]]$qual[7]),split = character(0)))
## the.qual
## score[the.qual]  # provides the quality
## library("ChIPsim") ### needed for decodeQuality
## -10*log10(decodeQuality(aln1[[1]]$qual[7], type = c("Sanger")))

#####################################################

################## for eaxh BAM file DO:

aln1 <- scanBam(bam.file,param=params)
aln1
aln1[[1]]

regions<-names(aln1)
regions<-strsplit(gsub("-",":",regions),split=":")
regions
length(regions)
dim(summary)
bam.file  # bam.file <- "realign_cleaned_normal.sort.bam"

qual.filter<-15 # bases with a Phred less than 20 are ignored 20-99% 10-90% ok  15 appears to be about what GATK is using 
alleles<-{}
alleles.freq<-{}
coverage<-{}
i<-1

# based on the mpos, and cigar count along to extract the base at the reference location

for (i in 1:length(regions)){
   loc.wanted<-as.numeric(regions[[i]][3]) # loc.wanted<-101401645
   posns<-loc.wanted-aln1[[i]]$pos+1

   ## queryLocs2refLocs(aln1[[i]]$pos,aln1[[i]]$cigar, loc.wanted)  ##  not implemented
   ## queryLoc2refLoc(aln1[[i]]$pos[1],aln1[[i]]$cigar[1], loc.wanted)  ##  not implemented
   
  ######################### get poistion on READ that corresponds to the reference location wanted 
   posns
   test<-cigarOpTable(aln1[[i]]$cigar)
   print(paste("test no deletions: position:",i,":",sum(apply(test[,c("I","D","N","H","P")],2,sum))),sep="")

   shift<-rep(0,times=length(posns))
   length<-cigarWidthAlongQuerySpace(aln1[[i]]$cigar) #length<-cigarToQWidth(aln1[[i]]$cigar)
   true.length<-cigarWidthAlongQuerySpace(aln1[[i]]$cigar,after.soft.clipping=TRUE) #  true.length<-cigarToWidth(aln1[[i]]$cigar)
  

   true.length[true.length>length]<-length[true.length>length]
 
   soft.has.S<-grepl("S",aln1[[i]]$cigar) ## 2S  aln1[[i]]$cigar[soft.has.S]
   soft.at.end<-grepl("S$",aln1[[i]]$cigar,fixed=FALSE) # aln1[[i]]$cigar[soft.at.end]
   shift[soft.has.S & !soft.at.end]<-(length-true.length)[soft.has.S & !soft.at.end] # aln1[[i]]$cigar[soft.has.S & !soft.at.end]
   ## shift<-(length-true.length)
   ## shift[soft.at.end]<-0
   ## shift
  
   posns2<-cigarQNarrow(aln1[[i]]$cigar, start=posns,width=1)
   ## posns
   posns2<-attr(posns2,"rshift")
   ## posns2
   ## posns-posns2-1

   shift.other <- posns-posns2-1
   shift.other[posns2==0]<-shift[posns2==0] ### posn can be within the softclip which causes a problem

   posns<-posns+shift+(shift.other-shift)
   
############
   seq.at.posn<-subseq(aln1[[i]]$seq,start=posns,width=1)
   ## as.character(seq.at.posn)
   seq.at.posn<-as.character(seq.at.posn)
   
   qual.at.posn<-subseq(aln1[[i]]$qual,start=posns,width=1)
   ## as.character(qual.at.posn)
 
   qual.at.posn<-score[as.character(qual.at.posn)]
   bad.qual<-qual.at.posn <= qual.filter  # as.character(seq.at.posn)[bad.qual]
   ## sort(qual.at.posn,decreasing=TRUE)
   ## sum(qual.at.posn> qual.filter)

   to.print<-qual.at.posn
   names(to.print)<-seq.at.posn
   print(tapply(seq.at.posn,seq.at.posn,length))
   print(to.print)
   sort(aln1[[i]]$qname)

   #TEST a single read with IGV etc
   ## aln1[[i]]$qname[grep("2388",aln1[[i]]$qname)]
   ## to.print[grep("2388",aln1[[i]]$qname)]
   ## aln1[[i]]$pos[grep("2388",aln1[[i]]$qname)]

   ########### removes bad quality based < 15
   seq.at.posn<-seq.at.posn[!bad.qual]
   qual.at.posn<-qual.at.posn[!bad.qual]
   
   names(qual.at.posn)<-seq.at.posn
   mean.quals<-round(tapply(qual.at.posn, names(qual.at.posn),mean))

   counts<-tapply(seq.at.posn,seq.at.posn,length)
   counts<-sort(counts,decreasing=TRUE)
   mean.quals<-mean.quals[names(counts)]
   ## alleles[i]<-toString(paste(unlist(labels(counts)),counts,sep=":"))
   alleles[i]<-toString(paste(unlist(labels(counts)),":",counts," (",mean.quals,")",sep=""))
   if(length(counts)>1){
   alleles.freq[i]<-counts[2]/counts[1]
   coverage[i]<-counts[2]+counts[1]
 }else{
      alleles.freq[i]<-0
   coverage[i]<-counts[1]
    }
   ## alleles[i]
   ## toString(labels(counts))
   print(toString(paste(unlist(labels(counts)),":",counts," (",mean.quals,")",sep="")))
   

 }



names(alleles)<-names(aln1)
names(alleles.freq)<-names(aln1)
names(coverage)<-names(aln1)
summary.key<-paste(summary[,"chr"],":",summary[,"start"],"-",summary[,"end"],sep="")

alleles
coverage

posns<-match(summary.key,names(alleles))
missing<-is.na(posns)
sum(missing)
posns


alleles<-alleles[posns]
alleles.freq<-round(alleles.freq[posns],digits=3)
coverage<-coverage[posns]
summary.full<-cbind(summary,alleles,alleles.freq,coverage)

summary.full


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END
