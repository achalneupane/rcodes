####################################### build a sliging window
###################
################## take existing passing genotypes and make a window about then 



## snpinfo<-snpinfo.ori[pass,]
## sort(table(snpinfo[,"gene"]))

sliding.window.generation <- function(a.indel,radius=300){
## annotations.pass<-a.indel[pass,c("chr","start","end","Gene.Names")]

annotations.pass<-a.indel[,c("chr","start","end","Gene.Names")]
the.chrs.slide<-unique(annotations.pass[,"chr"])

window.all<-{}

iwindow<-1
for( iwindow in 1:length(the.chrs.slide)){
    print(paste0("Doing sliding window for : ",the.chrs.slide[iwindow] ))
sum(annotations.pass[,"chr"]!=the.chrs.slide[iwindow] )

annotations.pass.chr<-subset(annotations.pass,subset=(annotations.pass[,"chr"]==the.chrs.slide[iwindow])) #  annotations.pass.chr<-subset(annotations.pass,subset=(annotations.pass[,"chr"]=="test"))

if(dim(annotations.pass.chr)[1]>1){
order.by<-order(as.numeric(annotations.pass.chr[,"start"]))   
annotations.pass.chr<-annotations.pass.chr[order.by,]

print(dim(annotations.pass.chr))
loci = IRanges(start=as.numeric(annotations.pass.chr[,"start"]),end=as.numeric(annotations.pass.chr[,"end"]))
window<-loci+radius
#window<-reduce(window)
# a.hist<-hist(width(window))
# window[1:5,]

## counts<-countOverlaps(window,loci) ## could exclude windows with only one SNPs here is wanted
## length(counts) ## could exclude windows with only one SNPs here is wanted
# window<-window[counts>1,] ## could exclude windows with only one SNPs here is wanted

#overlap<-findOverlaps(loci,window)
overlap<-findOverlaps(window,loci)

## (overlap)[1:20]
## queryHits(overlap)[1:10]
## subjectHits(overlap)[1:10]
## length(queryHits(overlap))
## queryHits(overlap)

test.common<-tapply(subjectHits(overlap),queryHits(overlap),function(x) toString(sort(x)))
#test.common[1:5]

dups<-duplicated(test.common)
test.common<-test.common[!dups] ## remove windows with the same SNPs
#test.common[1:5]

use.windows<-queryHits(overlap) %in% names(test.common)
#use.windows[1:20]

gene<-paste(the.chrs.slide[iwindow],queryHits(overlap)[use.windows],sep="_")
cluster<-gene
Name<-rownames(annotations.pass.chr)[subjectHits(overlap)[use.windows]]

a.window<-cbind(Name,gene,cluster)
if(iwindow==1){window.all<-a.window}else{window.all<-rbind(window.all,a.window)}

}else{ ## only one SNP
    gene<- annotations.pass.chr[1,"Gene.Names"]
    cluster<-paste(the.chrs.slide[iwindow],"1",sep="_")
    Name<-rownames(annotations.pass.chr)
   a.window<-cbind(Name,gene,cluster) 
} 

}# iwindow loop

posns<-match(window.all[,"Name"],rownames(annotations.pass))
missing<-is.na(posns)
sum(missing)
#rownames(annotations.pass)[missing][1:10]
window.all[,"gene"]<-annotations.pass[posns,"Gene.Names"]
rownames(window.all)<-window.all[,"Name"]             


return(window.all)

}
## dim(window.all)
## length(window)
## table(counts)

## getwd()
## colnames(window.all)[1]<-"Name"
## save(list=c("window.all"),file=paste0("sliding_window",snp.file,".RData")
  
