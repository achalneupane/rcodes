

M

us musculus Ly-49Q mRNA for NK receptor Ly-49Q, complete cds.


Ly49a U34890 chr6:130,314,238-130,330,680
Ly49A (uc009eha.1) at chr6:129922100-130330735 - Mus musculus natural killer cell receptor Ly49A (Ly49A) mRNA, complete cds.



b AF253057  chr6:131,170,007-131,195,344
Klra2 (uc012eua.1) at chr6:131169253-131197380 - killer cell lectin-like receptor 2 isoform 1


c U49868  chr6:130,273,385-130,287,625
g AF307946 chr6:130,168,906-130,181,642
l AF307947  chr6:129,994,043-130,256,497
q AB193833 chr6:129,781,465-129,824,925

Ly-49Q mRNA for NK receptor Ly-49Q, complete cds. chr6:129781463-129824950


e AY620247 chr6:129,849,145-129,861,479
x
i AF237686  chr6:130,128,705-130,143,054
\alpha
\gamma


129,781,343-130,330,680

129,781,465-130,335,000


## Hi Paul and Mhairi,

## I have two sets of annotation files for mouse exome. They are saved under the HiSeq directory structure on the Brown share:

## B:\Brown Group\UQCCG-IlluminaHiSeqProcessing\Exome Design and Annotation Documents
##         \100803_MM9_exome_rebal_2_EZ_HX1
##                 \100803_MM9_exome_rebal_2_EZ_HX1.bed
##                 \100803_MM9_exome_rebal_2_EZ_HX1.gff
##         \101122_MM9_exome_L2R_D02_EZ_HX1
##                 \101122_MM9_exome_L2R_D02_EZ_HX1.gff
        

## The first set of probes we are using are from the 100803_MM9_exome_rebal_2_EZ_HX1 for which there are .bed and .gff files.

## Kind regards,



options(width=150,max.print=200)
library(GenomicFeatures)

#--------------- WARNING THIS on hg18
file.dir<-"/media/scratch/Genomes" # this is hg18
dir(getwd())
bed.input<-"2.1M_Human_Exome.gff"

### to use hg19 need to use the below- use auto_loop_over_UQCCG.r to process
## bed.input<-"/media/scratch/Genomes/Human_Exome_Targets_hg19Lift.list" ## above mapped to hg19 "Human_Exome_Targets_v1_hg19_targets.list.RData"
## output.bed<-"Human_Exome_Targets_v1_hg19_capture_targets.bed"
### these done in 
#--------
getwd()
load("Human_Exome_Targets_v1_hg19_targets.list.RData")

#----------
file.dir<-"/media/scratch/Genomes" # this is hg19
setwd(file.dir)
dir(getwd())
bed.input<-"Human_Exome_Targets_hg19Lift.list"
skip.lines<-0
dict.file<-"/media/scratch/Genomes/UCSC_ALL_FULL_CHROMS_HG19.dict"
output.file<-"Human_Exome_Targets_Nimble_v1_hg19_snp.calls.list"
output.RData<-"Human_Exome_Targets_Nimble_v1_hg19_targets"
#------------


#----------
file.dir<-"/media/scratch/Genomes/PN-FC-121-1008-to-FC-121-1960 Illumina TruSeqExome" # this is hg19
setwd(file.dir)
dir(getwd())
bed.input<-"TruSeq_exome_targeted_regions.hg19.bed.chr"
skip.lines<-0
dict.file<-"/media/scratch/Genomes/UCSC_ALL_FULL_CHROMS_HG19.dict"
output.file<-"Human_Exome_Targets_illumina_v2_hg19_snp.calls.list2"
output.RData<-"Human_Exome_Targets_illumina_v2_hg19_targets2"
#------------

#----------
file.dir<-"/media/scratch/Genomes/PN-05 836 611 001_NimbleGen SeqCap EZ Exome Library v2.0" # this is hg19
setwd(file.dir)
dir(getwd())
bed.input<-"SeqCap_EZ_Exome_v2.gff"
skip.lines<-1
dict.file<-"/media/scratch/Genomes/UCSC_ALL_FULL_CHROMS_HG19.dict"
output.file<-"Human_Exome_Targets_Nimble_v2_hg19_snp.calls.list"
output.RData<-"Human_Exome_Targets_Nimble_v2_hg19_targets"
#------------

#----------
file.dir<-"/media/scratch/Genomes/Exome Cature annotations/SeqCapEZ_Exome_v3.0_Design_Annotation_files" # this is hg19
setwd(file.dir)
dir(getwd())
bed.input<-"SeqCap_EZ_Exome_v3.gff"
skip.lines<-1
dict.file<-"/media/scratch/Genomes/UCSC_ALL_FULL_CHROMS_HG19.dict"
output.file<-"Human_Exome_Targets_Nimble_v3_hg19_snp.calls.list"
output.RData<-"Human_Exome_Targets_Nimble_v3_hg19_targets"
#------------


#----------
file.dir<-"/media/scratch/Genomes/100803_MM9_exome_rebal_2_EZ_HX1" # this is mm9
setwd(file.dir)
dir(getwd())
bed.input<-"100803_MM9_exome_rebal_2_EZ_HX1.gff"
skip.lines<-1
dict.file<-"/media/ga-apps/UQCCG/Data/Genomes/mm9/mm9.dict"
output.file<-"Mouse_Exome_Targets_Nimble_v100803.list"
output.RData<-"Mouse_Exome_Targets_Nimble_v100803_targets"
#------------


#----------
file.dir<-"/media/UQCCG/Data/Nextera" # this is hg19
setwd(file.dir)
dir(getwd())
bed.input<-"NexteraRapidCapture_Exome_TargetedRegions.bed"
skip.lines<-1
dict.file<-"/media/UQCCG/Sequencing/Data/Genomes/hg19/UCSC_ALL_FULL_CHROMS_HG19.dict"
output.file<-"NexteraRapidCapture_Exome_TargetedRegions_hg19_snp.calls.list2"
output.RData<-"NexteraRapidCapture_Exome_TargetedRegions_hg19_targets2"
#-----------

#----------
file.dir<-"/media/UQCCG/Data/Nextera" # this is hg19
setwd(file.dir)
dir(getwd())
bed.input<-"NexteraRapidCapture_ExpandedExome_TargetedRegions.bed"
skip.lines<-0
dict.file<-"/media/UQCCG/Sequencing/Data/Genomes/hg19/UCSC_ALL_FULL_CHROMS_HG19.dict"
output.file<-"NexteraRapidCapture_ExpandedExome_TargetedRegions_hg19_snp.calls.list"
output.RData<-"NexteraRapidCapture_ExpandedExome_TargetedRegions_hg19_targets"
#-----------

#----------
file.dir<-"/media/UQCCG/Sequencing/Data/Genomes/Exome Cature annotations/Netera/" # this is hg19
setwd(file.dir)
dir(getwd())
bed.input<-"/media/UQCCG/Sequencing/Data/Genomes/Exome Cature annotations/Netera/nexterarapidcapture_exome_targetedregions_v1.2.bed"
skip.lines<-0
dict.file<-"/media/UQCCG/Sequencing/Data/Genomes/hg19/UCSC_ALL_FULL_CHROMS_HG19.dict"
output.file<-"NexteraRapidCapture_Exome_TargetedRegions_v1.2_hg19_snp.calls.list2"
output.RData<-"NexteraRapidCapture_Exome_TargetedRegions_v1.2_hg19_targets"
#-----------

## load("/media/UQCCG/Sequencing/Data/Genomes/hg19/NexteraRapidCapture_Exome_TargetedRegions_hg19_targets.RData")
## data
## 37105383
## dim(data) 212158
## data.gr
data[1:5,]
sum(as.numeric(data[,"end"])-as.numeric(data[,"start"]) )

## /media/UQCCG/Sequencing/Data/Genomes/Exome Cature annotations/Netera/nexterarapidcapture_exome_targetedregions_v1.2.bed
## /media/UQCCG/Sequencing/Data/Nextera/NexteraRapidCapture_Exome_TargetedRegions.bed

dict<-read.delim(dict.file,header=F,skip=0,sep="\t",fill=TRUE)
dict
dict<-dict[grep("\\@SQ",dict[,1]),]
max.values<-as.numeric(gsub("LN:","",dict[,3]))
names(max.values)<-gsub("SN:","",dict[,2])

max.values

## dict.file<-"UCSC_ALL_FULL_CHROMS_HG18.dict"
## output.file<-"Human_Exome_Targets_hg18.list"
## skip.lines<-1

 skip.lines #check!
options(show.error.messages = TRUE)
chromo1<-read.delim(bed.input,header=F,nrows=1,skip=skip.lines,sep="\t",fill=TRUE)
num.vars<-dim(chromo1)[2]
chromo1

regions<-try(scan(bed.input,what=character(num.vars),skip=skip.lines,sep="\t",fill=TRUE))
num.lines<-length(regions)/(num.vars)
dim(regions)<-c(num.vars,num.lines)
regions<-t(regions)

regions[1:10,]
tail(regions)
######### FOR NIMBLEGEN ONLY 
wanted<-grep("capture_target",regions[,3])  ## also capture target are given"
regions<-regions[wanted,]
###############

if(bed.input=="Human_Exome_Targets_hg19Lift.list") # liftover file to parse
  {
    the.bed<-strsplit(regions,split=":")
    the.chr<-unlist(lapply(the.bed,function(x) x[1]))
    the.posn<-unlist(lapply(the.bed,function(x) x[2]))
    the.posn<-strsplit(the.posn,split="-")
    the.start<-unlist(lapply(the.posn,function(x) x[1]))
    the.end<-unlist(lapply(the.posn,function(x) x[2]))
    regions<-cbind(the.chr,the.chr,the.chr,the.start,the.end) # strange so in nimblegen format see LINE 120 Below
  }





regions[1:5,]

test.col<-1  # this should be the column that contains the chromosome name
tapply(regions[,test.col],regions[,test.col],length)
max.values


chromosomes<-tapply(regions[,test.col],regions[,test.col],length)
chromosomes[!(names(chromosomes) %in% names(max.values))] ## these are ones I don't want...


wanted<- regions[,test.col] %in% names(max.values)
sum(!wanted)
regions<-regions[wanted,]  ### just get the regular chromosomes
tapply(regions[,test.col],regions[,test.col],length)


regions[1:5,]

######### FOR NIMBLEGEN ONLY 
data<-cbind(regions[,1],regions[,4],regions[,5],"+") # chr start end "+"
###########################

######### FOR ILLUMINA ONLY 
data<-cbind(regions[,1],regions[,2],regions[,3],regions[,6])

data<-cbind(regions[,1],regions[,2],regions[,3],regions[,6],regions[,4])
###########################

######### FOR NEXTERA ONLY 
data<-cbind(regions[,1],regions[,2],regions[,3],"+",regions[,4])
###########################

######### FOR NEXTERA v1.2 ONLY 
data<-cbind(regions[,1],regions[,2],regions[,3],"+",NA)
###########################

colnames(data)<-c( "chr","start","end","strand","gene")
 data[1:5,]
##

sum( (as.numeric(data[,"end"])- as.numeric(data[,"start"]) ) < 0 )  # check should be zero

core.ann<-c( "chr","start","end","strand") #'elementMetadata' cannot use "seqnames", "ranges", "strand", "seqlengths", "start", "end", "width", or "element" as column names
data.gr<-GRanges(seqnames =data[,"chr"],ranges = IRanges(start=as.numeric(data[,"start"]),end=as.numeric(data[,"end"])),strand=data[,"strand"])

values(data.gr)<-data[,!(colnames(data) %in% core.ann)]

data.gr
paste(output.RData,"RData",sep=".")
save(list=c("data","data.gr"),file=paste(output.RData,"RData",sep="."))

data.gr<-data.gr+200 # expand regions all + stand so ok to to this way 100 for 2.1 200 for seqcap
length(data.gr) # 244619-seqcap
data.gr<-reduce(data.gr)
length(data.gr)  # 152953-secap

data<-as.data.frame(data.gr)

data[1:5,]
tapply(data[,"seqnames"],data[,"seqnames"],length)

## setwd(genome.dir)
## dir(getwd())

## dict<-read.delim(dict.file,header=F,skip=0,sep="\t",fill=TRUE)
## dict<-dict[grep("\\@SQ",dict[,1]),]
## max.values<-as.numeric(gsub("LN:","",dict[,3]))
## names(max.values)<-gsub("SN:","",dict[,2])


wanted.order<-names(max.values) # paste("chr",c(1:22,"X","Y","M"),sep="")

sum(! (unique(seqnames(data.gr)) %in% wanted.order)) # no mitochromdria
data<-{}

#### assumes start and end are on the "+" strand:
for(i in 1:length(wanted.order)){
temp<-as.data.frame( data.gr[seqnames(data.gr)==wanted.order[i],])
too.big<-temp[,"end"]>max.values[wanted.order[i]]
too.small<-temp[,"start"]<1

way.too.big<-temp[,"start"]>max.values[wanted.order[i]]
way.too.small<-temp[,"end"]<1

if(sum(too.big)>0){temp[too.big,"end"]<-max.values[wanted.order[i]]}
if(sum(too.small)>0){temp[too.small,"end"]<-1}


print(paste(wanted.order[i],":",sum(too.big),"end points exceeded chr length:",sum(way.too.big), "problems",sep=" "))
print(paste(wanted.order[i],":",sum(too.big),"  start  points   less  than 1:",sum(way.too.small), "problems",sep=" "))                   
## max(temp[,"end"])
data<-rbind(data,temp)
}


### chec
dim(data)
data[1:5,]
## tapply(data[,"seqnames"],data[,"seqnames"],length)
match(wanted.order,data[,"seqnames"])
unique(data[,"seqnames"])
sum(is.na(data[,"seqnames"]))
tapply(data[,"seqnames"],data[,"seqnames"],length)

data.write<-paste(data[,"seqnames"],paste(data[,"start"],data[,"end"],sep="-"),sep=":")
data.write[1:5]
write.table(data.write,file=output.file,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

########################## to combine nimblegen and ilumina


## data.nimble<-data.gr  ### expanded resgions for nimblegen 
## data.illumina<-data.gr

## data.gr<-c(data.illumina,data.nimble)
## data.gr<-reduce(data.gr)

## data<-as.data.frame(data.gr)
## data[1:5,]
## data.write<-paste(data[,"seqnames"],paste(data[,"start"],data[,"end"],sep="-"),sep=":")
## write.table(data.write,file="Human_Exome_Targets_Combines_Nimble_illumina_hg19_snp.calls.list",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

file.info(sample.bams)

##

sum( (as.numeric(data[,"end"])- as.numeric(data[,"start"]) ) < 0 )  # check

core.ann<-c( "chr","start","end","strand") #'elementMetadata' cannot use "seqnames", "ranges", "strand", "seqlengths", "start", "end", "width", or "element" as column names
data.gr<-GRanges(seqnames =data[,"chr"],ranges = IRanges(start=as.numeric(data[,"start"]),end=as.numeric(data[,"end"])),strand=data[,"strand"])

cov<-coverage(the.reads,width=human.chromlens) # readBamGappedAlignments use this is have a RangesList and RangedData objects where they are the length of each list

 hist(values(the.reads)[,"isize"],xlim=c(-800,800),breaks=seq(from=-1000,to=1000,by=10))

  median(abs(values(the.reads)[,"isize"]),na.rm=TRUE)
targets<-cov$chr1 > 10

 x<- slice(cov$chr1,lower=10,upper=Inf)
  cbind(start(x),end(x),width(x),viewMeans(x))



chr1:229995177-230070035 #hg18
chr1:231928554-232003412  #hg19
  chr1:231928064-232003412

  chr1:231928064 231934643  6580 2794.335410
  chr1:231935771 232004452 68682 1038.300064

output.file<-"QBI-DISC1"

output.file<-"Human_Exome_Targets_QBI_DISC1_hg19_snp.calls.list"
output.RData<-"Human_Exome_Targets_QBI_DISC1_hg19_targets"
data<-cbind("chr1",c(231928064,231935771),c(231934643,232004452),"+")

colnames(data)<-c( "chr","start","end","strand")
 data


  
 targets<- slice(cov,lower=50,upper=Inf)
 targets<-lapply(targets,function(x) x<-x+150)

  x<-
    targets[["chr1"]]+150

  lapply(targets, function(x) print(x))

  x<-x+150
x<-reduce(x)
 as.data.frame(x)




## my.loci<-IRanges(start=as.numeric(a.indel[,"start"]),end=as.numeric(a.indel[,"end"]))
## on.chr<-as.vector(seqnames(possible.loci.ori))==a.indel[1,"chr"] # paste("chr",a.indel[1,"chr"],sep="")
## possible.space<-IRanges(start=as.numeric(start(possible.loci.ori))[on.chr],end=as.numeric(end(possible.loci.ori))[on.chr])
## common.loci<-overlapsAny(my.loci,possible.space)



############################# regions in common

library("BSgenome.Hsapiens.UCSC.hg19")
the.chroms<-seqlengths(Hsapiens)
the.chromo<-names(cov.v2)
 human.chromlens<-the.chroms[the.chromo]


load("Human_Exome_Targets_Nimble_v1_hg19_targets.RData")

data.grv2<-data.gr
load("Human_Exome_Targets_Nimble_v3_hg19_targets.RData")



#      ir.gerp<-IRanges(start=gerp.chr[,"position"],width=1)
library("BSgenome.Hsapiens.UCSC.hg19")
the.chroms<-seqlengths(Hsapiens)
the.chromo<-names(cov.v2)
 human.chromlens<-the.chroms[the.chromo]


data.grv2<-reduce(data.grv2)
data.gr<-reduce(data.gr)
      cov.v2<-coverage(data.grv2,weight=5,width=human.chromlens) ## Problems exist here is the ir.gerp are not unique get coverage*weight
      cov<-coverage(data.gr,weight=5,width=human.chromlens) ## Problems exist here is the ir.gerp are not unique get coverage*weight

      cov<-cov[names(cov.v2)]
      cov.all<-cov+cov.v2

   x<- slice(cov.all,lower=9,upper=11) ### over is 10 otherwise is 5
 x[[5]][1:5]



########check
overlaps<-subsetByOverlaps(data.grv2,data.gr, type = "within")
data.grv2[1:50,]
data.gr[1:50,]
overlaps[[13]]

cbind(start(data.gr[seqnames(data.gr)=="chr13",])[50:55],end(data.gr[seqnames(data.gr)=="chr13",])[50:55])
data.grv2[seqnames(data.grv2)=="chr13",]
x[[5]][1:5]

 
regionViews<-x
regionViews[[chromo]][1]
the.counts<-{}
order.chromos<-names(regionViews)

# ik<-13
for (ik in 1:length(order.chromos)){
  chromo<-order.chromos[ik]
 ## print( chromo)

  a.chr<-chromo
  a.start<-start(regionViews[[chromo]])
  a.end<-end(regionViews[[chromo]])
  a.width<-width(regionViews[[chromo]])

 
  a.set<-cbind(a.chr,a.start,a.end,a.width)

  the.counts<-rbind(the.counts,a.set)
}
colnames(the.counts) <- c("chr","start","end","length")
## }) # system.time
the.counts[1:5,]

write.table(the.counts,file="Common_target_loci_between_Nimblegen_v2.and.v3.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
getwd()
