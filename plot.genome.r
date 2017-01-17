###########################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"

the.dir<-"/media/scratch2/Scleroderma"
assoc.file<-"scl_related.genome" # "genome file 
#################################################

gen<-read.delim(paste(the.dir,assoc.file,sep="/"),header=T,sep="",fill=TRUE,stringsAsFactors=FALSE)

order.by<-order(gen[,"PI_HAT"],decreasing=TRUE)
gen<-gen[order.by,]


 par(mfrow=c(1,1),font=2,font.lab=2,font.axis=2,mgp=c(3.1,1,0),mar=c(5,5,4,2)+0.1)
######## all data : plot(as.numeric(vec[,eig1]),as.numeric(vec[,eig2]))
bp<-plot(1:dim(gen)[1],gen[,"PI_HAT"],cex=1.35,lwd=1.0,xlab="Samples", ylab="PI_HAT",cex.lab=2.0,cex.axis=2.0)

test<-identify(1:dim(gen)[1],gen[,"PI_HAT"],labels=as.character(gen[,"FID1"]))

setwd(the.dir)
savePlot(paste(assoc.file,"Related.tiff",sep="_"),,type="tiff")
savePlot(paste(assoc.file,"Related.png",sep="_"),type="png")


related<-gen[,"PI_HAT"]> 0.3
sum(related)

to.test<-gen[related,]
test1<- sort(tapply(to.test[,"FID1"],to.test[,"FID1"],length),decreasing=TRUE)
test2<- sort(tapply(to.test[,"FID2"],to.test[,"FID2"],length),decreasing=TRUE)
test1[1:50]
test2[1:50]
sum(test1>1)
sum(test2>1)
sum(test2>1)

the.related<-to.test

mutiples.to.remove<-unique(names(test)[test1>1])

the.related.new<-the.related[ !(the.related[,"FID1"] %in% mutiples.to.remove),]
dim(the.related.new)
dim(the.related)
length(mutiples.to.remove)

to.test<-the.related.new
test1<- sort(tapply(to.test[,"FID1"],to.test[,"FID1"],length),decreasing=TRUE)
test2<- sort(tapply(to.test[,"FID2"],to.test[,"FID2"],length),decreasing=TRUE)
test1[1:50]
test2[1:50]
sum(test1>1)
sum(test2>1)
######### if all zero the remove final samples

all.related<-unique(c(mutiples.to.remove,to.test[,"FID1"]))

length(all.related) # 295 for oats
getwd()
write.table(cbind(all.related,all.related),file=paste(assoc.file,"Related.txt",sep="_"),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
