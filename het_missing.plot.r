###########################################FILTER het vs miss #####################################


working.dir<-"/media/scratch2/Scleroderma"
missing.file<-"first.missing.imiss"
het.file<-"first.het.het"
fam.file<-"scl_controls.fam"

#################################################

setwd(working.dir)

miss.data<-read.table(missing.file,header=T,sep="",fill=TRUE)
het.data<-read.table(het.file,header=T,sep="",fill=TRUE)
fam<-read.table(fam.file,header=F,sep="",fill=TRUE)

colnames(fam)<-c("FID","IID","MOTHER","FATHER","SEX","STATUS")
table(fam[,6])

 het<- (as.numeric(het.data[,"N.NM."]) - as.numeric(het.data[,"O.HOM."]))/as.numeric(het.data[,"N.NM."])
 colnames(miss.data)[ colnames(miss.data)=="F_MISS"]<-"miss"

het.data<-cbind(het.data,het)

posns<-match(het.data[,"IID"],miss.data[,"IID"])
missing<-is.na(posns)
sum(missing)


het.miss<-cbind(miss.data[!missing,],het.data[posns[!missing],])

het.miss[1:5,]

rownames(het.miss)<-as.character(het.miss[,"IID"])
#colnames(het.miss)<-c("s1","s2","s3","s4","s5","s6","het","s8","miss")

posns<-match(het.data[,"IID"],fam[,"IID"])
missing<-is.na(posns)
sum(missing)

ann.het<-fam[posns[!missing],]
rownames(ann.het)<-ann.het[,"IID"]

######reorder sample.ann so same as vec
ann.het[1:5,]


 par(mfrow=c(1,1),font=2,font.lab=2,font.axis=2,mgp=c(3.1,1,0),mar=c(5,5,4,2)+0.1)

color.with<-"STATUS"
#color.with<-"Study"
#color.with<-"Product"
#color.with<-"Gender"
tapply(rownames(ann.het),ann.het[,color.with],length) # counts

class.array<-tapply(rownames(ann.het),ann.het[,color.with])
classes<-levels(as.factor(ann.het[,color.with]))  #order ok with the above
color.set<-rainbow(length(classes))          #get auto colors
color.array<-color.set[class.array] #array of colors of length of number of samples

range.eig1<-range(het.miss[,"het"]) #  range.eig1<-c(0.32,0.37)
range.eig2<-range(het.miss[,"miss"])
bp<-plot(het.miss[,"het"],het.miss[,"miss"],pch=20, cex=1,col = color.array, xlim=range.eig1,ylim=range.eig2,xlab="Heterozygosity",ylab="Missingness",main=paste("Heterozygosity vs. Missingness",color.with,sep=" : "))
legend(range.eig1[1],range.eig2[2],classes,col=color.set,pch=19,cex=0.8)    #,bg="white" # e1 vs e3
abline(v=0.343,lwd=2,col="red")
abline(v=0.36,lwd=2,col="red")
## abline(v=0.34,lwd=1,col="blue")
## abline(v=0.36,lwd=1,col="blue")
abline(h=0.025,lwd=2,col="red")

savePlot(paste(fam.file,"Het_Miss.tiff",sep="_"),,type="tiff")
savePlot(paste(fam.file,"Het_Miss.png",sep="_"),type="png")



miss.thresh.high<-0.025
miss.thresh.low<-0
het.thresh.low<-0.343
het.thresh.high<-0.36

samples.miss.cut.high<-rownames(het.miss)[het.miss[,"miss"]>miss.thresh.high]
samples.miss.cut.low<-rownames(het.miss)[het.miss[,"miss"]<miss.thresh.low]


sample.het.cut.low<-rownames(het.miss)[het.miss[,"het"]<het.thresh.low]
sample.het.cut.high<-rownames(het.miss)[het.miss[,"het"]>het.thresh.high]

length(samples.miss.cut.high)
length(samples.miss.cut.low)

length(sample.het.cut.high)
length(sample.het.cut.low)

sample.het.miss.cut<-unique(c(samples.miss.cut.high, samples.miss.cut.low,sample.het.cut.low, sample.het.cut.high))    # 123 removed      188 for David
length(sample.het.miss.cut)

getwd()
write.table(cbind(sample.het.miss.cut,sample.het.miss.cut),file=paste(fam.file,"het.miss.bad.txt",sep="_"),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

#################### NUMBERS FILTERED ##############################
