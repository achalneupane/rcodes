




##################################################
#################################################
## interna code of auto_loop_over_UQCCG
## processes a cov object
## need rl.from.gr, and target file information rl.from.gr<-rl.from.gr[names(cov)] ## must have the same chromom names in both files
## rl.from.gr<-as(data.gr, "RangesList")
if(sum(names(rl.from.gr) %in% names(cov))==0){names(rl.from.gr)<-gsub("chr","",names(rl.from.gr))}
rl.from.gr<-rl.from.gr[names(rl.from.gr) %in% names(cov)]


cov<-cov[names(rl.from.gr)]
## total.coverage<-sum(as.numeric(sum(rl.from.gr)))

  
regionViews <- RleViewsList(rleList = cov, rangesList =rl.from.gr )# fast way to get info anout exons

## regionViews.expanded <- RleViewsList(rleList = cov, rangesList =rl.from.gr.expanded )

## chrom.cov<-viewSums(regionViews)
## on.target.bases<-sum(as.numeric(unlist((lapply(chrom.cov,sum)))))
## percent.on.target.bases<-100*on.target.bases/total.coverage

## chrom.cov.expanded<-viewSums(regionViews.expanded)
## on.target.bases.expanded<-sum(as.numeric(unlist((lapply(chrom.cov.expanded,function(x) sum(as.numeric(x),na.rm=TRUE))))))
## percent.on.target.bases.expanded<-100*on.target.bases.expanded/total.coverage
#rl.from.gr[[3]]



########################report sum/man/mean counts in all regions (exons) without multicore
## system.time({
the.counts<-{}
order.chromos<-names(rl.from.gr)
ik<-1

for (ik in 1:length(order.chromos)){
  chromo<-order.chromos[ik]
 print( chromo)

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

 ## mulitcore mccollect data : data returns as a list
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
## the.counts[1:5,]
## a.set[1:5,]

 ccds.length<-(as.numeric(the.counts[,"length"]))
 ccds.gt.1<-(as.numeric(the.counts[,"bp.gt.1"]))
 ccds.gt.5<-(as.numeric(the.counts[,"bp.gt.5"]))
 ccds.gt.10<-(as.numeric(the.counts[,"bp.gt.10"]))
 ccds.gt.15<-(as.numeric(the.counts[,"bp.gt.15"]))
 ccds.gt.30<-(as.numeric(the.counts[,"bp.gt.30"]))
 
 percent.ccds.gt.1<-100*ccds.gt.1/ccds.length
 percent.ccds.gt.5<-100*ccds.gt.5/ccds.length
 percent.ccds.gt.10<-100*ccds.gt.10/ccds.length
 percent.ccds.gt.15<-100*ccds.gt.15/ccds.length
 percent.ccds.gt.30<-100*ccds.gt.30/ccds.length

the.counts<-cbind(the.counts,percent.ccds.gt.1, percent.ccds.gt.5, percent.ccds.gt.10,percent.ccds.gt.15,percent.ccds.gt.30)
#ccds.gt.5[1:5]

#colnames(the.counts)<-output.labels

the.counts[1:5,]
## "chr,start,end,length,Sum,Max,Mean,bp.gt.1,bp.gt.5,bp.gt.10,bp.gt.15,bp.gt.30,percent.ccds.gt.1,percent.ccds.gt.5,percent.ccds.gt.10,percent.ccds.gt.15,percent.ccds.gt.30"
## paste(colnames(the.counts),collapse=",")
