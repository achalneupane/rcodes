




##################################################
#################################################
## interna code of auto_loop_over_UQCCG
## processes a cov object
## need rl.from.gr, and target file information
## rl.from.gr<-as(data.gr, "RangesList")
if(sum(names(rl.from.gr) %in% names(cov))==0){names(rl.from.gr)<-gsub("chr","",names(rl.from.gr))}
 
cov<-cov[names(rl.from.gr)]
total.coverage<-sum(as.numeric(sum(cov)))

if(is.na(total.coverage)){ # interger overlow occured
total.coverage<-sum(unlist(lapply(cov,function(x) sum(as.numeric(runValue(x)) * as.numeric(runLength(x)),na.rm = TRUE))),na.rm=TRUE)
}

#### excessive coverage can break this command
## > cov$chr1
## integer-Rle of length 249250621 with 13286012 runs
##   Lengths: 10000     3    22     2     5 ...    27    20    85    12 10025
##   Values :     0     2     8    10    14 ...     0     2     4     2     0
## > sum(cov$chr1)
## [1] NA
## Warning message:
## In sum(runValue(x) * runLength(x), ..., na.rm = na.rm) :
##   integer overflow - use runValue(.) <- as.numeric(runValue(.))
## > sum(cov$chr22)
## [1] 410699360
  
regionViews <- RleViewsList(rleList = cov, rangesList =rl.from.gr )# fast way to get info anout exons

regionViews.expanded <- RleViewsList(rleList = cov, rangesList =rl.from.gr.expanded )

chrom.cov<-viewSums(regionViews)
on.target.bases<-sum(as.numeric(unlist((lapply(chrom.cov,sum)))))
percent.on.target.bases<-100*on.target.bases/total.coverage

chrom.cov.expanded<-viewSums(regionViews.expanded)
on.target.bases.expanded<-sum(as.numeric(unlist((lapply(chrom.cov.expanded,function(x) sum(as.numeric(x),na.rm=TRUE))))))
percent.on.target.bases.expanded<-100*on.target.bases.expanded/total.coverage


print(on.target.bases.expanded)
print(percent.on.target.bases.expanded)
print(percent.on.target.bases)
print("-------------------")


########################report sum/man/mean counts in all regions (exons) without multicore
## system.time({
the.counts<-{}
order.chromos<-names(rl.from.gr)
for (ik in 1:length(order.chromos)){
  chromo<-order.chromos[ik]
 ## print( chromo)

  a.chr<-chromo
  a.start<-start(regionViews[[chromo]])
  a.end<-end(regionViews[[chromo]])
  a.width<-width(regionViews[[chromo]])

  ### multicore set up

  ## m <- do.call(rbind, mclapply(d, function(x){ merge(x, y, by = c('a', 'b'))},mc.silent = TRUE))
  
  a.sum<-mcparallel(viewSums(regionViews[[chromo]]))
  a.max<-mcparallel(viewMaxs(regionViews[[chromo]]))
  a.mean<-mcparallel(viewMeans(regionViews[[chromo]]))
  thr1<-mcparallel(viewApply(regionViews[[chromo]],function(x) basesAboveThreshold(x,1)))
  thr5<-mcparallel(viewApply(regionViews[[chromo]],function(x) basesAboveThreshold(x,5)))
  thr10<-mcparallel(viewApply(regionViews[[chromo]],function(x) basesAboveThreshold(x,10)))
  thr15<-mcparallel(viewApply(regionViews[[chromo]],function(x) basesAboveThreshold(x,15)))
  thr30<-mcparallel(viewApply(regionViews[[chromo]],function(x) basesAboveThreshold(x,30)))

 ## mulitcore collect data : data returns as a list
  a.sum<-unlist(mccollect(a.sum))
  a.max<-unlist(mccollect(a.max))
  a.mean<-unlist(mccollect(a.mean)) 
  thr1<-unlist(mccollect(thr1)) 
  thr5<-unlist(mccollect(thr5))
  thr10<-unlist(mccollect(thr10))
  thr15<-unlist(mccollect(thr15))
  thr30<-unlist(mccollect(thr30))

  a.set<-cbind(a.chr,a.start,a.end,a.width,a.sum,a.max,a.mean,thr1,thr5,thr10,thr15,thr30)

  the.counts<-rbind(the.counts,a.set)
}

colnames(the.counts) <- c("chr","start","end","length","Sum","Max","Mean","bp.gt.1","bp.gt.5","bp.gt.10","bp.gt.15","bp.gt.30")
## }) # system.time

 ccds.length<-sum(as.numeric(the.counts[,"length"]))
 ccds.gt.1<-sum(as.numeric(the.counts[,"bp.gt.1"]))
 ccds.gt.5<-sum(as.numeric(the.counts[,"bp.gt.5"]))
 ccds.gt.10<-sum(as.numeric(the.counts[,"bp.gt.10"]))
 ccds.gt.15<-sum(as.numeric(the.counts[,"bp.gt.15"]))
 ccds.gt.30<-sum(as.numeric(the.counts[,"bp.gt.30"]))
 on.target.bases<-sum(as.numeric(the.counts[,"Sum"]))
 percent.on.target.bases<-100*on.target.bases/total.coverage
print(percent.on.target.bases)

 percent.ccds.gt.1<-100*ccds.gt.1/ccds.length
 percent.ccds.gt.5<-100*ccds.gt.5/ccds.length
 percent.ccds.gt.10<-100*ccds.gt.10/ccds.length
 percent.ccds.gt.15<-100*ccds.gt.15/ccds.length
 percent.ccds.gt.30<-100*ccds.gt.30/ccds.length

 median.max.coverage<-median(as.numeric(the.counts[,"Max"]))
 median.mean.coverage<-median(as.numeric(the.counts[,"Mean"]))
 mean.mean.coverage<-mean(as.numeric(the.counts[,"Mean"]))
print("Mean")
 data<-data.frame(c(a.sample.ID,a.sample,a.recipe,a.capture,a.lane,a.run,targets.file,sum.total.reads,total.reads,total.dup.reads,total.unmapped.reads,total.coverage,percent.on.target.bases.expanded,on.target.bases, percent.on.target.bases, median.max.coverage, median.mean.coverage,mean.mean.coverage,percent.ccds.gt.1,percent.ccds.gt.5,percent.ccds.gt.10,percent.ccds.gt.15,percent.ccds.gt.30,100*total.dup.reads/(total.reads+total.unmapped.reads+total.dup.reads),100*total.unmapped.reads/(total.reads+total.unmapped.reads+total.dup.reads),a.DS)) # need c() abou all so get a column vector and not a row vector

dim(data)

### all.QC.column.labels<-c("ID","Sample","Recipe","Capture.Method","Lane","Run","target_file","total_reads","total_mapped_reads","total_dup_reads","unmapped_reads","total.bases",paste("on.target.bases.",extend.exons,"bp",sep=""),"on.target.bases", "percent.on.target.bases","median.max.coverage", "median.mean.coverage","mean.mean.coverage","percent.ccds.gt.1","percent.ccds.gt.5","percent.ccds.gt.10","percent.ccds.gt.15","percent.ccds.gt.30","percent.ccds.gt.40","percent.ccds.gt.50","percent.ccds.gt.60","percent.ccds.gt.70","percent.ccds.gt.80","percent.ccds.gt.90","percent.ccds.gt.100","percent_Duplicated","percent_Unmapped","Description")

length(all.QC.column.labels)
rownames(data)<-all.QC.column.labels
print("woser 3")
 print(all.QC.column.labels)



  
