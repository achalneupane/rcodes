############################################################################################################################
############################################################################################################################
############################################################################################################################


options(width=200,mac.print=500)
# plink --bfile  Exome_Plus_AOGC_Gentrain_TOP --freq --remove   excluded.samples.txt  --out UQDI
# plink --bfile  AOGC_opticall_TOP --freq --remove   excluded.samples.txt  --out UQDI_optical

setwd("/media/scratch2/AOGC-NGS/ExomeChip")


gentrain<-read.table("UQDI.frq",header=T,fill=TRUE,stringsAsFactors=FALSE)
optical<-read.table("UQDI_optical.frq",header=T,fill=TRUE,stringsAsFactors=FALSE)

uqdi.snps<-read.table("UQDI.snps.txt",header=F,fill=TRUE,stringsAsFactors=FALSE)
uqdi.hwe<-read.table("UQDI.hwe",header=T,fill=TRUE,stringsAsFactors=FALSE)

uqdi.snps<-uqdi.snps[,1]
uqdi.snps[1:5]




dim(gentrain) # 23468     6
dim(optical) # 21898     6
tapply(optical[,"CHR"],optical[,"CHR"],length)
tapply(gentrain[,"CHR"],gentrain[,"CHR"],length)

posns<-match(gentrain[,"SNP"],optical[,"SNP"])
missing<-is.na(posns)
sum(missing) # 1570
gentrain[1:5,]
optical[1:5,]

colnames(optical)<-paste(colnames(optical),"_optical",sep="")

total<-cbind(gentrain,optical[posns,c("MAF_optical","NCHROBS_optical","SNP_optical","A1_optical","A2_optical")])

total[1:5,]

the.order<-order(total[,"MAF"],decreasing=TRUE)

total<-total[the.order,]


posns<-match(total[,"SNP"],uqdi.hwe[,"SNP"])
missing<-is.na(posns)
sum(missing) # 1570

uqdi.hwe[posns,][1:5,]


total<-cbind(total,uqdi.hwe[posns,c("GENO","O.HET.","E.HET.","P")])


write.table(total,file="AOGC_all_snp_frequencies.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

posns<-match(uqdi.snps,total[,"SNP"])
missing<-is.na(posns)
sum(missing)
total<-total[posns[!missing],]
dim(total) #  23468    11


write.table(total,file="AOGC_uqdi_snp_frequencies.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



in.common<-!is.na(total[,"SNP_optical"])
gentrain.only<-total[!in.common,]

common<-total[in.common,]


dim(common)
colnames(common)

plot(as.numeric(common[,"MAF"]),as.numeric(common[,"MAF_optical"]),type="p",xlab="Gentrain",ylab="Optical")

xrange<-c(0,0.005)
yrange<-c(0,0.005)

plot(as.numeric(common[,"MAF"]),as.numeric(common[,"MAF_optical"]),type="p",xlab="Gentrain",ylab="Optical",xlim=xrange,ylim=yrange)


test<-identify(as.numeric(common[,"MAF"]),as.numeric(common[,"MAF_optical"]),label=common[,"SNP"])
test<-identify(as.numeric(common[,"MAF"]),as.numeric(common[,"MAF_optical"]),label=paste("(",as.numeric(common[,"MAF"])," , ",as.numeric(common[,"MAF_optical"]),")",sep="") )

the.order<-order(common[test,"MAF"],decreasing=TRUE)
common[test[the.order],"SNP"]
common[test[the.order],1:7]

dim(common)
write.table(common[test[the.order],"SNP"],file="snp.list",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
getwd()
#y-y1=m(x-x1)
# y=m(x-x1)+y1
#intercept=0.5(*-x1) +y1


m<-1

## point.high<-common[test,]
## inter.high=as.numeric(point.high[,"MAF_optical"])-m*as.numeric(point.high[,"MAF"])
## inter.high #

plot(as.numeric(common[,"MAF"]),as.numeric(common[,"MAF_optical"]),type="p",xlab="Gentrain MAF",ylab="Optical MAF",cex.lab=1.5)

inter.high<-0.0003
abline(a=inter.high,b=m,col="magenta")

inter.high<-0.0003
abline(a=inter.high,b=m,col="magenta")


## point.high<-common[test,]
## inter.low=as.numeric(point.high[,"MAF_optical"])-m*as.numeric(point.high[,"MAF"])
## inter.low

inter.low<--0.0003
abline(a=inter.low,b=m,col="magenta")


inter.low<--0.001
abline(a=inter.low,b=m,col="magenta")

below.diag<-function(x,intercept,m){
  #
 y<-m*as.numeric(x[,"MAF"]) + intercept
# print(cbind(y,as.numeric(x[,"MAF_optical"]) ))
 y > as.numeric(x[,"MAF_optical"]) | is.na(as.numeric(x[,"MAF_optical"])) | is.na(as.numeric(x[,"MAF"]))
   }

above.diag<-function(x,intercept,m){
  #
 y<-m*as.numeric(x[,"MAF"]) + intercept
# print(cbind(y,as.numeric(x[,"MAF_optical"]) ))
 y < as.numeric(x[,"MAF_optical"]) | is.na(as.numeric(x[,"MAF_optical"])) | is.na(as.numeric(x[,"MAF"]))
   }


gentrain.thresh<-0.05
below<-below.diag(common,-0.0003,1)
x.thresh<-(as.numeric(common[,"MAF"]) < gentrain.thresh) & !is.na(as.numeric(common[,"MAF_optical"])) & !is.na(as.numeric(common[,"MAF"]))
sum(is.na(x.thresh))
sum(is.na(below))
sum(below & x.thresh)
below.low<-below & x.thresh
sum(below.low)
points(as.numeric(common[below.low,"MAF"]),as.numeric(common[below.low,"MAF_optical"]),col="blue")

below<-below.diag(common,-0.001,1)
x.thresh<-(as.numeric(common[,"MAF"]) > gentrain.thresh) & !is.na(as.numeric(common[,"MAF_optical"])) & !is.na(as.numeric(common[,"MAF"]))
sum(x.thresh)
below.high<-below & x.thresh
sum(below.high)
points(as.numeric(common[below.high,"MAF"]),as.numeric(common[below.high,"MAF_optical"]),col="green")




  

optical.thresh<-0.05

special<-(as.numeric(common[,"MAF_optical"]) > optical.thresh) & (as.numeric(common[,"MAF"]) < 0.002 ) & !is.na(as.numeric(common[,"MAF_optical"])) & !is.na(as.numeric(common[,"MAF"]))
sum(special)
points(as.numeric(common[special,"MAF"]),as.numeric(common[special,"MAF_optical"]),col="orange")

above<-above.diag(common,0.0003,1)
y.thresh<-(as.numeric(common[,"MAF_optical"]) < gentrain.thresh) & !is.na(as.numeric(common[,"MAF_optical"])) & !is.na(as.numeric(common[,"MAF"]))
sum(above & y.thresh)

above.low<-above & y.thresh & !special
sum(above.low)
points(as.numeric(common[above.low,"MAF"]),as.numeric(common[above.low,"MAF_optical"]),col="magenta")

above<-above.diag(common,0.001,1)
y.thresh<-(as.numeric(common[,"MAF_optical"]) > optical.thresh) & !is.na(as.numeric(common[,"MAF_optical"])) & !is.na(as.numeric(common[,"MAF"]))
sum(is.na(y.thresh))
sum(is.na(above))
sum(above & y.thresh)

above.high<-above & y.thresh & !special
sum()
points(as.numeric(common[above.high,"MAF"]),as.numeric(common[above.high,"MAF_optical"]),col="red")

low.thresh<-0.05
low.min<-0

y.thresh<-(as.numeric(common[,"MAF_optical"]) > low.min) & (as.numeric(common[,"MAF"]) > low.min) & (as.numeric(common[,"MAF_optical"]) < low.thresh) & (as.numeric(common[,"MAF"]) < low.thresh) & !above.high & !above.low & !special & !below.high & !below.low & !is.na(as.numeric(common[,"MAF_optical"])) & !is.na(as.numeric(common[,"MAF"]))


low<-y.thresh & (as.numeric(common[,"P"])< 1e-2)
sum(low)
points(as.numeric(common[low,"MAF"]),as.numeric(common[low,"MAF_optical"]),col="brown")
common[1:5,]



types<- c("above.high","above.low","special","below.high","below.low","low")


all.cases<-above.high | above.low | special | below.high | below.low | low
sum(all.cases)

root.dir<-"/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/cluster_viz/"
root.dir<-"/media/scratch2/AOGC-NGS/ExomeChip/cluster_viz"
# i<-1
 the.chrs<-unique(common[,"CHR"])

for(i in 1:length(types)){
   test<-eval(as.name(types[i]))
   print(types[i])
   print(sum(test))
  write.table(common[test,"SNP"],file=paste(root.dir,types[i],"_snp.list",sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

   the.chrs<-unique(common[,"CHR"])

   # ic<-1
   for(ic in 1:length(the.chrs)){
     on.chr<-common[,"CHR"]==the.chrs[ic]
     print(the.chrs[ic])
     print(sum(test & on.chr))
     write.table(common[test & on.chr,"SNP"],file=paste(root.dir,"/",types[i],".chr",the.chrs[ic],"_snp.list",sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
     ## write.table(common[test & on.chr,"SNP"],file=paste(root.dir,"chr",the.chrs[ic],"/",types[i],".chr",the.chrs[ic],"_snp.list",sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
   }

   
}

savePlot("AOGC_cluster_check.jpeg",type="jpeg")
savePlot("AOGC_cluster_check.png",type="png")
savePlot("AOGC_cluster_check.tiff",type="tiff")
## abline(v=0.004,col="red")
dim(common)

annotated

common<-cbind(common,above.high,above.low,special,below.high,below.low,low)

total[1:5,]

posns<-match(total[,"SNP"],common[,"SNP"])
dim(total)
length(posns)

common[1:5,]



total.ann<-cbind(total,common[posns,c("above.high","above.low","special","below.high","below.low","low","SNP")])
write.table(total.ann,file="uqdi.labeled.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)


## abline(h=0.2)
## abline(v=0.2)

options(width=2000,max.print=200)

sum(as.numeric(common[,"MAF"])!=0 & as.numeric(common[,"MAF_optical"])==0 & !is.na(as.numeric(common[,"MAF_optical"]))  )

sum(is.na(as.numeric(common[,"MAF_optical"])))


zero.test<-(as.numeric(common[,"MAF"])!=0 & as.numeric(common[,"MAF_optical"])==0 & !is.na(as.numeric(common[,"MAF_optical"]))  )
sum(a.test) ## 119
common[a.test,][1:5,1:7]

lower.test<-(as.numeric(common[,"MAF"])!=0 & as.numeric(common[,"MAF_optical"])==0 & !is.na(as.numeric(common[,"MAF_optical"]))  )


write.table(common[a.test,"SNP"],file="snp.list.zero.optical",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

7339

common[is.na(as.numeric(common[,"MAF_optical"])),]
