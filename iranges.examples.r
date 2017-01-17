
library("GenomicRanges")


load("/media/UQCCG/Sequencing/Data/Genomes/NexteraRapidCapture_Exome_TargetedRegions_v1.2_hg19_targets.RData")
## get cov data data.gr

data[1:5,]

### MAKE A GENOMIC RANGE
data.gr<-GRanges(seqnames =data[,"chr"],ranges = IRanges(start=as.numeric(data[,"start"]),end=as.numeric(data[,"end"])),strand=data[,"strand"])

## values(data.gr)<-data[,!(colnames(data) %in% core.ann)] ## add extra columns to a genomic ranges 
data.gr

data.gr<-data.gr+200 # expand regions all + stand so ok to to this way 100 for 2.1 200 for seqcap
length(data.gr) # 244619-seqcap
data.gr<-reduce(data.gr) # combines overlapping
length(data.gr)  # 152953-secap

data<-as.data.frame(data.gr)

data[1:5,]
# tapply(data[,"seqnames"],data[,"seqnames"],length)
--------------------------------------------------------------------

############################# regions in common between two regions:::

## colnames(data)<-c("chr","start","end")
## colnames(data2)<-c("chr","start","end")

## > data
##      chr    start end  
## [1,] "chr1" "10"  "30" 
## [2,] "chr1" "50"  "70" 
## [3,] "chr1" "50"  "100"
## > data2
##      chr    start end 
## [1,] "chr1" "10"  "30"
## [2,] "chr1" "60"  "70"
## [3,] "chr1" "10"  "40"

## library("GenomicRanges")

## data.gr<-GRanges(seqnames =data2[,"chr"],ranges = IRanges(start=as.numeric(data2[,"start"]),end=as.numeric(data2[,"end"])),strand="+")
## data.grv2<-GRanges(seqnames =data[,"chr"],ranges = IRanges(start=as.numeric(data[,"start"]),end=as.numeric(data[,"end"])),strand="+")

## data.grv2<-reduce(data.grv2)
## data.gr<-reduce(data.gr)




load("/media/UQCCG/Sequencing/Data/Genomes/hg19/Human_Exome_Targets_Nimble_v1_hg19_targets.RData")

data.grv2<-data.gr
load("/media/UQCCG/Sequencing/Data/Genomes/hg19/Human_Exome_Targets_Nimble_v3_hg19_targets.RData")

seqlevels(data.gr)
seqlevels(data.grv2)
#      ir.gerp<-IRanges(start=gerp.chr[,"position"],width=1)
library("BSgenome.Hsapiens.UCSC.hg19")
the.chroms<-seqlengths(Hsapiens)
the.chromo<-seqlevels(data.gr)
human.chromlens<-the.chroms[the.chromo]


data.grv2<-reduce(data.grv2)
data.gr<-reduce(data.gr)

data.gr  # a genomic range
data.grv2
########### make a coverge curve from a genomic range
      cov.v2<-coverage(data.grv2,weight=5,width=human.chromlens) ## Problems exist here is the ir.gerp are not unique get coverage*weight
      cov<-coverage(data.gr,weight=5,width=human.chromlens) ## Problems exist here is the ir.gerp are not unique get coverage*weight

cov

##### coverge curves don work if they contain different chromoeomse or the order is different (treat as a list)
      cov<-cov[names(cov.v2)]
      cov.all<-cov+cov.v2

   x<- slice(cov.all,lower=9,upper=11) ### over is 10 otherwise is 5

names(x)
 x[["chr10"]][1:5]
 x[[5]][1:5]


regionViews<-x
regionViews[[chromo]][1]
the.counts<-{}
order.chromos<-names(regionViews)

  ik<-1
#for (ik in 1:length(order.chromos)){
  chromo<-order.chromos[ik]
 ## print( chromo)

  a.chr<-chromo
  a.start<-start(regionViews[[chromo]])
  a.end<-end(regionViews[[chromo]])
  a.width<-width(regionViews[[chromo]])

 
  a.set<-cbind(a.chr,a.start,a.end,a.width)

  the.counts<-rbind(the.counts,a.set)
#
colnames(the.counts) <- c("chr","start","end","length")
## }) # system.time
the.counts[1:5,]


################################################ END


##############################################################################################
################# coverage like calcuation - specify a region

  rl.from.gr<-as(data.grv2, "RangesList") # COVERT  data.gr contains
rl.from.gr
  rl.from.gr.expanded<-rl.from.gr+30
  rl.from.gr.expanded<-reduce(rl.from.gr.expanded)



if(sum(names(rl.from.gr) %in% names(cov))==0){names(rl.from.gr)<-gsub("chr","",names(rl.from.gr))}


cov<-coverage(data.gr,weight=1,width=human.chromlens) 

cov<-cov[names(rl.from.gr)]  ## make sure in same order
total.coverage<-sum(as.numeric(sum(cov)))

#### this is important and useful
regionViews <- RleViewsList(rleList = cov, rangesList =rl.from.gr )
regionViews[1]
regionViews[[1]]

chromo<-order.chromos[ik]
#for (ik in 1:length(order.chromos)){
 a.start<-start(regionViews[[chromo]])
  a.end<-end(regionViews[[chromo]])
  a.width<-width(regionViews[[chromo]])

  ### multicore set up
basesAboveThreshold<-function(x,thresh=0){sum(x>=thresh)}
  ## m <- do.call(rbind, mclapply(d, function(x){ merge(x, y, by = c('a', 'b'))},mc.silent = TRUE))
  
  a.sum<-(viewSums(regionViews[[chromo]]))
  a.max<-(viewMaxs(regionViews[[chromo]]))
  a.mean<-(viewMeans(regionViews[[chromo]]))
  thr1<-(viewApply(regionViews[[chromo]],function(x) basesAboveThreshold(x,1)))
  thr5<-(viewApply(regionViews[[chromo]],function(x) basesAboveThreshold(x,5)))
  thr10<-(viewApply(regionViews[[chromo]],function(x) basesAboveThreshold(x,10)))
  thr15<-(viewApply(regionViews[[chromo]],function(x) basesAboveThreshold(x,15)))
  thr30<-(viewApply(regionViews[[chromo]],function(x) basesAboveThreshold(x,30)))
#}


################################### have a rangelist and see if has an overlapping region
############################ EXTRA THE flt# IDT FROM OUR STANDard ANNOATED FILES 

a.indel.gr<-GRanges(seqnames =a.indel[,"chr"],ranges = IRanges(start=as.numeric(a.indel[,"start"]),end=as.numeric(a.indel[,"end"])),strand="+")

the.chr<-"chr13"
begin<-28608214
end<-28608336
#########

indel.gr<-GRanges(seqnames =the.chr,ranges = IRanges(start=as.numeric(begin),end=as.numeric(end)),strand="+")

regionViews[[1]]
?findOverlaps
common<-overlapsAny(a.indel.gr,indel.gr) 
sum(common)
is.flat.geno<-grepl(":flat$",a.indel[,"TYPE"])

#### USE THE BOOLEAN ARRAY common in the orginal daat from...
a.indel[common & !is.flat.geno,1:4]  
a.indel[common & !is.flat.geno,the.samples.HC.GT][1:2,]
a.indel[common & !is.flat.geno,the.samples.HC.GT][,1]
a.indel[common & !is.flat.geno,the.samples.HC.GQ]
hits<-apply(a.indel[common & !is.flat.geno,the.samples.HC.GT],2,function(x) sum(x!="0/0"))

sum(common)


## hits<-apply(cg[common,the.samples],2,function(x) sum(x!="00"))


## hits<-apply(a.indel[common & !is.flat.geno,c(the.samples.HC.GT,the.samples.HC.GQ)],1,function(x){
##  # print(x)
##   a.geno<-x!="0/0" & grepl("GT$",names(x))
##   samples<-names(x)[a.geno]
##   sample.GQ<-gsub(".GT",".GQ",samples)
##   paste(x[c(samples,sample.GQ)],collapse="::")
## })



############################ using findOverlaps
 library("IRanges")

ir<-IRanges(start=starts.chunk,end=ends.chunk)
# ir<-ir+50
ir.gerp<-IRanges(start=gerp.chr[,"position"],width=1)

subject<-findOverlaps(ir,ir.gerp)
query<-queryHits(subject)
subject<-subjectHits(subject)


##############################



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
