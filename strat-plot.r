
# load sample annotations
load("D:\\Research\\Matt Brown\\EIGENSOFT\\Merge ALL + HapMap Aug 2008\\sample.ann.RData")
load("sample.ann.RData")
#filtered data
pop.strat0<-read.delim("pop.strat.phase0.txt",header=F,sep="",fill=TRUE)
pop.strat0<-as.character(pop.strat0[,1])

related<-read.delim("related_samples.txt",header=F,sep="",fill=TRUE)
related<-as.character(related[,1])

related.miss<-read.delim("related_miss_samples.txt",header=F,sep="",fill=TRUE)
related.miss<-as.character(related.miss[,1])

filtered<-c(pop.strat0,related)


snps.p6.bad<-read.delim("p6_remove_bad_snps.txt",header=F,sep="",fill=TRUE)
snps.p6.bad<-as.character(snps.p6.bad[,1])

snps.bad<-read.delim("bad.snps.txt",header=F,sep="",fill=TRUE)
snps.bad<-as.character(snps.bad[,1])

snps.p6<-read.delim("p6_remove.txt",header=F,sep="",fill=TRUE)
snps.p6<-as.character(snps.p6[,1])

snps.p6<-setdiff(snps.p6.bad,snps.bad)

write.table(snps.p6,"p6_remove.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)



############################## READ IN DATA FILE
beta
file<-"AS_ALL_common.no6p.F.pca.evec"     # all including HapMap
file<-"AS_ALL_common.noR.no6p.F.pca.evec"
file<-"AS_ALL_common.noRHM.no6p.F.pca.evec"
file<-"AS_ALL_common.EIGEN.noLD.F.pca.evec"
file<-"AS_ALL_common.EIGEN.noLD.F.pca.evec"
file<-"cluster4_PCA.txt"
  #LOAD PLOT.EIGEN.RDATA in old merge all setwd("/media/scratch2/DORITHs_WORK/PCA STRANGE")


# > sample.ann[1:5,]
#                 Sample.ID Call.Rate Study         Product Ethnicity Gender     origin
# 1893411004_A 1893411004_A    0.9967   117 Human1Mv1_C.bpm       CEU      F US.control
# 1893411005_A 1893411005_A    0.9975   117 Human1Mv1_C.bpm       CEU      F US.control
# 1893411006_A 1893411006_A    0.9983   117 Human1Mv1_C.bpm       CEU      M US.control
# 1893429006_A 1893429006_A    0.9979   117 Human1Mv1_C.bpm       CEU      F US.control
# 1894819007_A 1894819007_A     0.998   117 Human1Mv1_C.bpm       CEU      F US.control
vec[1:5,]
vec<-read.delim(file,header=T,sep="",fill=TRUE)
rownames(vec)<-unlist(lapply(strsplit(rownames(vec),split=":"),function(x) x[1]))   # fix colnum names

# vec[,dim(vec)[2]+1]<-NA
# vec[,dim(vec)[2]+1]<-NA
# vec[,dim(vec)[2]+1]<-NA
#rownames(vec)<-vec[,1]
#vec<-vec[,-1]
###########################################################################################
###########################################################################################
######### need file for eigenstrat run:
write.table(cbind(rownames(vec),rownames(vec)),"eigenstrat_keep_samples.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
no.ibd<-read.delim("NO_IBD_CASES.txt",header=F,sep="\t",fill=TRUE)
no.ibd<-as.character(no.ibd[,1])
posns<-match(rownames(vec),no.ibd) # get no.ibd in rownames(vec)
no.ibd.keep<-no.ibd[posns[!is.na(posns)]]
no.ibd.case.fam<-cbind(no.ibd.keep,no.ibd.keep,0,0,0,1)


> length(no.ibd.keep)
[1] 1159
> length(no.ibd)
[1] 2964
## load in sample with controls
###########################################################################################
###########################################################################################
###########################################################################################


status<-read.delim("AS_ALL_common.EIGEN.with6p.F.ind",header=T,sep="",fill=TRUE)
status<-read.delim("AS_ALL_common.noIBD.with6p.F.ind",header=T,sep="",fill=TRUE)
controls<-as.character(status[status[,3]=="Control",1])
no.ibd.controls.fam<-cbind(controls,controls,0,0,0,2)

no.ibd.list<-c(no.ibd.keep,controls)
write.table(cbind(no.ibd.list,no.ibd.list),"no.ibd.keep_samples.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
                        
##############################################################

colnames(vec)<-c(paste("e.",1:10,sep=""),"status")   #,"origin","color","points")

dim(vec)
test<-grep("TRUE", rownames(sample.ann)==rownames(vec))

#########missing sample check
 posns<-match(rownames(sample.ann),rownames(vec))
 missing<-rownames(sample.ann[is.na(posns),])
length(missing)   # missing can be related, stratificatio,pop strat
setdiff(related.miss.1,missing)   # this should be sero
#############################

posns<-match(rownames(vec),rownames(sample.ann))
sum(is.na(posns)) # problem is sample id not found  as vec samples are a sunset of all known samples


######reorder sample.ann so same as vec
ann<-sample.ann[posns,]


############### vec and ann now in the same order
############################## READ IN DATA FILE
  ann[1:5,]
            Sample.ID Call.Rate Study       Product Ethnicity Gender     origin
WTCCC66069 WTCCC66069      <NA>   58C HumanHap550v1 Caucasian      M UK.control
 #par(mfrow=c(1,1),font=2,font.lab=2,font.axis=2,mgp=c(3.5,1,0),mar=c(5,5,4,2)+0.1)
 #par(mfrow=c(2,1),font=2,font.lab=2,font.axis=2,mar=c(1.0,4.5,1.25,2.1),mgp=c(1.5,0.5,0) )
 # par(mfrow=c(3,3),font=1,font.lab=2,font.axis=2,mgp=c(2,1,0),mar=c(5,3,4,2)+0.1)

eig1<-"e.1"
eig2<-"e.2"
color.with<-"origin"
color.with<-"Study"
color.with<-"Product"
color.with<-"Gender"
tapply(rownames(ann),ann[,color.with],length) # coulsnts

class.array<-tapply(rownames(ann),ann[,color.with])
classes<-levels(as.factor(ann[,color.with]))  #order ok with the above
color.set<-rainbow(length(classes))          #get auto colors
color.array<-color.set[class.array] #array of colors of length of number of samples


# pch=19: solid circle,
# pch=20: bullet (smaller circle),
# pch=21: circle,
# pch=22: square,
# pch=23: diamond,
# pch=24: triangle point-up,
# pch=25: triangle point down.


range.eig1<-range(vec[,eig1])
range.eig2<-range(vec[,eig2])
                                                                                           cut
# new.range<-locator()
# range.eig1<-new.range$x
# range.eig2<-new.range$y

eig1<-"e.2"
eig2<-"e.4"
#bp<-plot(range.eig1,range.eig2,xlab=eig1,ylab=eig2,main="Principle Components: Case+Control")

bp<-plot(vec[,eig1],vec[,eig2],pch=20, cex=1,col = color.array, xlim=range.eig1,ylim=range.eig2,xlab=eig1,ylab=eig2,main=paste("Principle Components",color.with,sep=" : "))

#bp<-plot(range(vec[,eig1]),range(vec[,eig2]),xlab="Eigenvector 1",ylab="Eigenvector 3",main="Principle Components: Case+Control")

#leg.txt<-c("US control","WTCCC control","UK affected","US affected","AUS affected")
#legend(-0.052,0.035,leg.txt,col=c("blue","light blue","orange","red","purple"),pch=19)    #,bg="white" # e1 vs e2
legend(range.eig1[1],-0.02,classes,col=color.set,pch=19,cex=0.9)
savePlot(filename="AS_noLD_e3.e4.origin.jpeg",type="jpeg")

legend(range.eig1[1],range.eig2[2],classes,col=color.set,pch=19,cex=0.9)    #,bg="white" # e1 vs e3


median(chi.data[!p6,"Chisq-EIG-2-corr"],na.rm=TRUE)/0.456
[1] 1.154841
> mean(chi.file[!p6,"Chisq_corr"],na.rm=TRUE)
[1] 1.153154
> mean(chi.data[!p6,"Chisq-EIG-4-corr"],na.rm=TRUE)
[1] 1.084971
> mean(chi.data[!p6,"Chisq-EIG-10-corr"],na.rm=TRUE)
[1] 1.085195

check<-rownames(vec)[vec[,"e.4"]< -0.02]

 check
  [1] "WTCCC66161"    "WTCCC66704"    "WTCCC66711"    "WTCCC88887"   
  [5] "WTCCC66238"    "WTCCC66277"    "WTCCC67004"    "WTCCC88319"   
  [9] "WTCCC88472"    "WTCCC88901"    "WTCCC88957"    "1485986170_A" 
 [13] "1485986215_A"  "1495201148_A"  "1518442198_A"  "1573151218_A" 
 [17] "1610071005_A"  "1508738273_A"  "1513810012_A"  "1516439174_A" 
 [21] "1573142584_A"  "1573151242_A"  "1573151421_A"  "1573169251_A" 
 [25] "1544641485_A"  "1629919076_A"  "1629919077_A"  "1629919166_A" 
 [29] "1629919265_A"  "1637564151_A"  "1656298517_A"  "1671121369_A" 
 [33] "1671121374_A"  "1671121503_A"  "1673012169_A"  "1673012178_A" 
 [37] "1673012323_A"  "1673012327_A"  "1673012477_A"  "1673047005_A" 
 [41] "1673047034_A"  "1673047092_A"  "1729987600_A"  "1739501496_A" 
 [45] "1740094352_A"  "1743164136_A"  "1544641164_A"  "1561862513_A" 
 [49] "1613642268_A"  "1629919143_A"  "1629919171_A"  "1629919212_A" 
 [53] "1629935295_A"  "1629935304_A"  "1637564280_A"  "1637564369_A" 
 [57] "1637564490_A"  "1663181366_A"  "1671105184_A"  "1673012054_A" 
 [61] "1673012133_A"  "1673012240_A"  "1673012242_A"  "1673047088_A" 
 [65] "1729987375_A"  "1739448376_A"  "1740086154_A"  "1740086412_A" 
 [69] "1739448472_A"  "1740086345_A"  "1767676095_A"  "1767676243_A" 
 [73] "1767676240_A"  "1767676324_A"  "1791043267_A"  "1792927150_A" 
 [77] "1792927216_A"  "1767676456_A"  "1791043133_A"  "1792927030_A" 
 [81] "1792927335_A"  "1793031180_A"  "1792927593_A"  "1793031430_A" 
 [85] "1803147084_A"  "1868063866_A"  "1868160458_A"  "1868160271_A" 
 [89] "1868160720_A"  "1868160778_A"  "US_01-0029-01" "US_01-0105-01"
 [93] "US_01-0233-01" "US_01-0397-01" "US_01-0424-01" "US_01-0441-01"
 [97] "US_01-0556-01" "US_02-0135-01" "US_13-N042-01" "US_14-S537-01"
[101] "US_355"        "US_815"        "US_PAX-002"    "US_SFX-062"   
[105] "UKB902"        "UKM0316"       "US_01-0198-01" "US_01-0567-01"
[109] "US_02-H143-01" "US_12-0009-01" "US_12-0017-02" "US_13-N103-01"
[113] "US_204"        "US_231"        "US_796"        "US_836"       
[117] "US_LAX-003"    "US_PAX-019"   


################################### DOUBLE strat plot ######################

#Read model filee
 CHR         SNP   A1   A2     TEST            AFF          UNAFF        CHISQ   DF            P
   1   RS3934834    A    G     GENO    32/347/1039    52/584/1747      0.02302    2       0.9886
   1   RS3934834    A    G    TREND       411/2425       688/4078     0.004571    1       0.9461

file<-"AS_ALL_common.8.controlUKvus.model"
file<-"AS_ALL_common.8.controlUKvus.assoc"
#################################################
################################ DO LAM vs STRAT plot
############################### get lambda data
#cutoffs<- seq(-0.04,0.005,0.005)
lam<-cutoffs
for(i in 1:length(cutoffs)) {
file<-paste("AS_ALL_common.",i,".model",sep="")
assoc.file<-read.delim(file,header=T,sep="",fill=TRUE)
#lam[i]<-mean(assoc.file[assoc.file[,"TEST"]=="ALLELIC","CHISQ"],na.rm=TRUE)
lam[i]<-median(assoc.file[assoc.file[,"TEST"]=="ALLELIC","CHISQ"],na.rm=TRUE)/0.456
}
#######################################
file<-paste("AS_ALL_common.",i,".assoc",sep="")


par(mfrow=c(2,1),font=2,font.lab=2,font.axis=2,mar=c(1.0,4.5,1.25,2.1),mgp=c(1.5,0.5,0) )
par(mfrow=c(2,1),font=2,font.lab=2,font.axis=2,mar=c(3.0,4.5,1.25,2.1),mgp=c(1.5,0.5,0) )
eig1<-"e.1"
eig2<-"e.2"
color.with<-"origin"
color.with<-"Study"
color.with<-"Product"
color.with<-"Gender"
tapply(rownames(ann),ann[,color.with],length) # coulsnts

class.array<-tapply(rownames(ann),ann[,color.with])
classes<-levels(as.factor(ann[,color.with]))  #order ok with the above
color.set<-rainbow(length(classes))          #get auto colors
color.array<-color.set[class.array] #array of colors of length of number of samples

range.eig1<-range(vec[,eig1])
range.eig2<-range(vec[,eig2])

eig1<-"e.1"
eig2<-"e.2"

bp<-plot(vec[,eig1],vec[,eig2],pch=20, cex=1,col = color.array, xlim=range.eig1,ylim=range.eig2,xlab="",ylab="Eigenvector 2",main="Principle Components")

legend(range.eig1[1],range.eig2[2],classes,col=color.set,pch=19,cex=0.9)    #,bg="white" # e1 vs e3

lam.cols<-rainbow(length(lam))
for(i in 1:length(lam)) {
abline(v=cutoffs[i],col=lam.cols[i],lwd=2) }
 # par(mar=c(3.0,4.5,1.25,2.1),mgp=c(1.5,0.5,0) )

bp2<-plot(cutoffs,lam,xlim=range.eig1,xlab="Eigenvector 1",ylab=expression(paste(bold(lambda),bold(" (mean(chi-squared values))"))),type="l")
points(cutoffs,lam,col=lam.cols,pch=19,cex=2.0 )
for(i in 1:length(lam)) {
abline(v=cutoffs[i],col=lam.cols[i],lwd=2) }

# count<-cbind(cutoffs,cutoffs) ; colnames(count)<-c("case","control")
# for(i in 1:length(cutoffs)) {
# temp<-vec[ vec[,"e.1"]> cutoffs[i],"status"]
# count[i,]<-tapply(temp,temp,length)  }

count[,1]<-paste("Cases=",count[,1],sep="")
count[,2]<-paste("Control=",count[,2],sep="")
leg2.txt<-paste(count[,1],count[,2],sep=" ")
leg2.txt<-count[,1]

legend(0.0052,2.10,leg2.txt,col=lam.cols,pch=19,cex=0.7)

count2<-cbind(cutoffs,cutoffs,cutoffs,cutoffs,cutoffs) ; colnames(count2)<-c("hapmap","UK.case","UK.control","US.case","US.control")
for(i in 1:length(cutoffs)) {
temp<- ann[vec[,"e.1"]> cutoffs[i],"origin"]
count2[i,]<-tapply(temp,temp,length)  }

cbind(cutoffs,cut.length,count,lam,count2)







###########################################FILTER het vs miss #####################################
het.miss<-read.delim("het_missing.txt",header=F,sep="",fill=TRUE)
het.miss<-read.delim("het_missing2.txt",header=F,sep="",fill=TRUE)
rownames(het.miss)<-as.character(het.miss[,1])
colnames(het.miss)<-c("s1","s2","s3","s4","s5","s6","het","s8","miss")

posns<-match(rownames(het.miss),rownames(sample.ann))
sum(is.na(posns)) # problem is sample id not found  as het.miss samples are a sunset of all known samples


######reorder sample.ann so same as vec
ann.het<-sample.ann[posns,]


par(mfrow=c(2,1),font=2,font.lab=2,font.axis=2,mar=c(1.0,4.5,1.25,2.1),mgp=c(1.5,0.5,0) )

color.with<-"origin"
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
abline(v=0.342,lwd=2,col="red")
abline(v=0.358,lwd=2,col="red")

abline(v=0.34,lwd=1,col="blue")
abline(v=0.36,lwd=1,col="blue")
abline(h=0.021,lwd=2,col="red")

samples.miss.cut<-rownames(het.miss)[het.miss[,"miss"]>0.021]
sample.het.cut.low<-rownames(het.miss)[het.miss[,"het"]<0.34]
sample.het.cut.high<-rownames(het.miss)[het.miss[,"het"]>0.36]

samples.miss.cut<-rownames(het.miss)[het.miss[,"miss"]>0.021]
sample.het.cut.low<-rownames(het.miss)[het.miss[,"het"]<0.342]
sample.het.cut.high<-rownames(het.miss)[het.miss[,"het"]>0.358]

sample.het.miss.cut<-unique(c( samples.miss.cut,sample.het.cut.low, sample.het.cut.high))    # 123 removed      188 for David
length(sample.het.miss.cut)

SWAP:
"UKM0139",
"UKM0164",  "UKS164"
"UKM0053",  "UK79.3"
"UKM0376", "UK156.3"
"UKS468" ,"UK674.3"
grep("UK674.3",sample.het.miss.cut)

#################### NUMBERS FILTERED ##############################

tapply(sample.ann[,"origin"],sample.ann[,"origin"],length)
tapply(sample.ann[,"Gender"],sample.ann[,"Gender"],length)
tapply(sample.ann[,"Product"],sample.ann[,"Product"],length)
tapply(sample.ann[sample.het.miss.cut,"origin"],sample.ann[sample.het.miss.cut,"origin"],length)
tapply(sample.ann[pop.strat0,"origin"],sample.ann[pop.strat0,"origin"],length)





rel.counts<-sample.ann[related,]
dim(rel.counts)
tapply(rel.counts[,"origin"],rel.counts[,"origin"],length)
> tapply(rel.counts[,"origin"],rel.counts[,"origin"],length)
    hapmap    UK.case UK.control    US.case US.control
       135         25          5         37        214 
       
miss.het.counts<-sample.ann[sample.het.miss.cut,]
dim(miss.het.counts)
tapply(miss.het.counts[,"origin"],miss.het.counts[,"origin"],length)

       
pop.only<-setdiff(pop.strat0,related)       ### David cuts
pop0.counts<-sample.ann[pop.only,]
tapply(pop0.counts[,"origin"],pop0.counts[,"origin"],length)
> tapply(pop0.counts[,"origin"],pop0.counts[,"origin"],length)
    hapmap    UK.case UK.control    US.case US.control
       127          7          1         33        111 
       

intersect(pop.strat1,pop.strat0)      6033
pop1.only<-setdiff(pop.strat1,related)   5908
intersect(pop1.only,sample.het.miss.cut)
pop1.only.miss<-setdiff(pop1.only,sample.het.miss.cut) 5774

fin.counts<-sample.ann[pop1.only.miss,]
dim(fin.counts)
tapply(fin.counts[,"origin"],fin.counts[,"origin"],length)

################################################

miss.counts<-sample.ann[sample.het.miss.cut,]




Davids cuts
tapply(sample.ann[,"origin"],sample.ann[,"origin"],length)
    hapmap    UK.case UK.control    US.case US.control
       262       1198       1436        983       4149
tapply(miss.counts[,"origin"],miss.counts[,"origin"],length)   # red/black lines used these
   UK.case UK.control    US.case US.control
        14         80         15         79

> tapply(miss.counts[,"origin"],miss.counts[,"origin"],length)  #blue/black line
   UK.case UK.control    US.case US.control 
         4         76          7         35

> intersect(sample.het.miss.cut,related)
0   # for both cuts
setdiff(sample.het.miss.cut,related)
setdiff(related,sample.het.miss.cut)


intersect(pop.strat0,related)

########################BASIC FILTER
related.miss<-unique(c(related,sample.het.miss.cut))
length(related.miss)
[1] 597



related.miss.frame<-data.frame( one=related.miss ,two=related.miss )
write.table(related.miss.frame,"related_miss_samples.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

length(related.miss)+length(pop.strat0)   #833
related.miss.1<-unique(c(related.miss,pop.strat0))
length(related.miss.1)    #800


########################Other filters

cutoffs<- seq(-0.04,0.005,0.005)
cut.length<-cutoffs
for(i in 1:length(cutoffs)) {
temp<-rownames(vec)[vec[,"e.1"]> cutoffs[i]]
cut.length[i]<-length(temp)
temp.frame<-data.frame( one=temp ,two=temp )
write.table(temp.frame,paste("keep_num_",i,".txt",sep=""),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
}



####################################### SNP filter ############################

file<-"AS_ALL_common.8.controlUKvus.model"
file<-"AS_ALL_common.8.controlUKvus.assoc"

file<-"AS_ALL_common.8.snpF.assoc"
############################### get lambda data
assoc.file<-read.delim(file,header=T,sep="",fill=TRUE)

 mean( assoc.file[,"CHISQ"])

assoc.file[,"P"]<-log10(assoc.file[,"P"])*-1
assoc.file[1:5,]

> sum(assoc.file[,"P"]>1.8)
[1] 11770
> sum(assoc.file[,"P"]>1.9)
[1] 9914
> sum(assoc.file[,"P"]>2)
[1] 8389
> sum(assoc.file[,"P"]>2.5)
[1] 3682
> sum(assoc.file[,"P"]>3)
[1] 1616
> sum(assoc.file[,"P"]>3.5)
[1] 784
> sum(assoc.file[,"P"]>4)
[1] 359
> sum(assoc.file[,"P"]>4.5)
[1] 171
> sum(assoc.file[,"P"]>5)
[1] 99


> mean( assoc.file[assoc.file[,"P"]<20000,"CHISQ"])
[1] 1.404701
> mean( assoc.file[assoc.file[,"P"]<5,"CHISQ"])
[1] 1.393783
> mean( assoc.file[assoc.file[,"P"]<4.5,"CHISQ"])
[1] 1.389465
> mean( assoc.file[assoc.file[,"P"]<4.0,"CHISQ"])
[1] 1.379779
> mean( assoc.file[assoc.file[,"P"]<3.5,"CHISQ"])
[1] 1.360914
> mean( assoc.file[assoc.file[,"P"]<3,"CHISQ"])
[1] 1.330221
> mean( assoc.file[assoc.file[,"P"]<2.5,"CHISQ"])
[1] 1.268842
> mean( assoc.file[assoc.file[,"P"]<2.0,"CHISQ"])
[1] 1.161507
> mean( assoc.file[assoc.file[,"P"]<1.5,"CHISQ"])
[1] 0.9767761
> mean( assoc.file[assoc.file[,"P"]<1.0,"CHISQ"])
[1] 0.6915331
> mean( assoc.file[assoc.file[,"P"]<1.9,"CHISQ"])
[1] 1.132132
> mean( assoc.file[assoc.file[,"P"]<1.8,"CHISQ"])
[1] 1.098634


mean( assoc.file[assoc.file[,"P"]<5,"CHISQ"])


markers<-read.delim("Confirmation set markers.txt" ,header=T,sep="",fill=TRUE)
markers.case<-as.character(markers[,"Discovery"])
markers.case<-gsub("rs","RS",markers.case)
markers.eth<-read.delim("all_the_used_strat_snps_ONLY.txt" ,header=F,sep="",fill=TRUE)
markers.eth<-as.character(markers.eth[,1])

intersect(markers.case,markers.eth)


snps.cut<- seq(1.3,4.5,0.025)
lam.val <-snps.cut
num.lost <-snps.cut
for(i in 1:length(snps.cut)) {
num.lost[i]<-sum(assoc.file[,"P"]>snps.cut[i])
lam.val[i]<-mean( assoc.file[assoc.file[,"P"]<snps.cut[i],"CHISQ"])
}


#par(mfrow=c(2,1),font=2,font.lab=2,font.axis=2 )
xlimits<-range(snps.cut)
xlimits<-c(1.5,2.0)
ylimits.lost<-range(num.lost[snps.cut>xlimits[1] &  snps.cut<xlimits[2]])
ylimits.lam<-range(lam.val[snps.cut>xlimits[1] &  snps.cut<xlimits[2]])

par(mfrow=c(2,1),font=2,font.lab=2,font.axis=2,mar=c(3.0,4.5,1.25,2.1),mgp=c(1.5,0.5,0) )
 plot(snps.cut,num.lost,xlim=xlimits,ylim=ylimits.lost,col="red",ylab="Number snps rejected",main="WTCCC vs US control",xlab="")
#for(i in 1:length(snps.cut:5)) {
#abline(v=snps.cut[i],lwd=1) }

plot(snps.cut,lam.val,xlim=xlimits,ylim=ylimits.lam,col="green",ylab="lam value",xlab="-log10(P) cutoff")

 for(i in 1:length(lam)) {
abline(v=cutoffs[i],col=lam.cols[i],lwd=2) }


 cbind(snps.cut,lam.val,num.lost)


[12,]    1.525 0.9880968    18844
 [13,]    1.550 0.9994821    18046
 [14,]    1.575 1.0105751    17284
 [15,]    1.600 1.0212070    16568
 [16,]    1.625 1.0314065    15895
 [17,]    1.650 1.0419727    15211
 [18,]    1.675 1.0519227    14579
 [19,]    1.700 1.0610911    14008
 [20,]    1.725 1.0715068    13371
 [21,]    1.750 1.0806643    12821
 [22,]    1.775 1.0900575    12267
 [23,]    1.800 1.0986338    11770
 [24,]    1.825 1.1065069    11322
 [25,]    1.850 1.1152013    10836
 [26,]    1.875 1.1233020    10391
 [27,]    1.900 1.1321323     9914
 [28,]    1.925 1.1401140     9490
 [29,]    1.950 1.1470440     9128


 bad.snps<-as.character(assoc.file[assoc.file[,"P"]>1.6,"SNP"])
length(bad.snps)
[1] 16568
intersect(bad.snps, markers.case)

assoc.file[match(intersect(bad.snps, markers.case),assoc.file[,"SNP"]) ,]


  CHR        SNP        BP A1     F_A     F_U A2  CHISQ         P     OR
6065     1 RS11590808  61571647  A 0.05677 0.03008  G 32.450  7.912574 1.9410   NFIA
40092    2  RS7580110 181765328  A 0.39610 0.42680  G  6.756  2.029514 0.8809   LOC729026
70873    4  RS1320213  39869933  G 0.12870 0.10670  A  8.266  2.393726 1.2370
88676    5  RS7720838  40522653  C 0.39120 0.41880  A  5.492  1.718739 0.8919
149284   9 RS12353357  11249755  A 0.04261 0.02433  G 19.350  4.963371 1.7850
149649   9  RS3931884  13792750  A 0.44020 0.47040  G  6.426  1.949234 0.8852
153736   9   RS348472  74710880  A 0.04211 0.05714  G  8.025  2.336017 0.7253    ALDH1A1
165237  10  RS1326986  19969519  G 0.05459 0.02123  A 58.970 13.794796 2.6620   LOC10012864
188932  11  RS2714068 122904751  A 0.38950 0.43540  G 15.110  3.993106 0.8274   GRAMD1B
236524  16  RS3093406  27360201  A 0.04586 0.02840  G 15.750  4.141824 1.6440   IL21R
272022  20  RS6027755  58702105  A 0.36800 0.34010  G  6.026  1.851089 1.1300

CHR	SNP	 POS	  A1	A2 CHISQ_plink	CHISQ_plink_perm  CHISQ_eigen	P_plink	   P_plink_perm	  P_eigen
1    RS11590808	61571647   A	G    15.42	15.42	          10.0711	8.59E-05    8.63E-05	0.00151
2    RS7580110	181765328  A	G    21.54	21.54	          18.2794	3.47E-06    2.00E-06	1.91E-05
4    RS1320213	39869933   G	A    16	        16	           9.586	6.35E-05    8.54E-05    0.00196
5    RS7720838	40522653   C	A    19.44	19.44	          13.9277	1.04E-05    1.30E-05	0.00019
9    RS12353357	11249755   A	G    18.12	18.12	          12.4387	2.07E-05    1.60E-05	0.000421
9    RS3931884	13792750   A	G    15.93	15.93	          17.868	6.56E-05    6.08E-05	2.37E-05
9    RS348472 74710880     A    G    19.13	19.13	           6.0125	1.22E-05    1.00E-05	0.0142
10   RS1326986	19969519   G	A    15.2	15.2	           12.0145	9.70E-05    8.26E-05	0.000528
11   RS2714068	122904751  A	G    15.21	15.21	           11.6422	9.60E-05    0.000118	0.000645
16   RS3093406	27360201   A	G    15.53	15.53	           14.7664	8.14E-05    0.0001275	0.000122
20   RS6027755	58702105   A	G    17.36	17.36	           10.7609	3.09E-05    4.21E-05	0.00104


4	RS1320213	39869933	G	A	16	16	9.586	6.35E-05	8.54E-05	0.00196

# snps.p6.bad<-read.delim("p6_remove_bad_snps.txt",header=F,sep="",fill=TRUE)
# snps.p6.bad<-as.character(snps.p6.bad[,1])
# 
# snps.bad<-read.delim("bad.snps.txt",header=F,sep="",fill=TRUE)
# snps.bad<-as.character(snps.bad[,1])
# 
# snps.p6<-read.delim("p6_remove.txt",header=F,sep="",fill=TRUE)
# snps.p6<-as.character(snps.p6[,1])
# 
# snps.p6<-setdiff(snps.p6.bad,snps.bad)

length(snps.p6.bad)
length( bad.snps)
snps.p6.bad.control<-unique(c(bad.snps,snps.p6.bad))
length(snps.p6.bad.control)

write.table(bad.snps,"bad_control_snps.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(snps.p6.bad.control,"p6_remove_bad_control_snps.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)








































for(i in 1:length(lambda)) {
abline(v=cut.e1e2[i],col=lam.cols[i],lwd=2) }
  par(mar=c(3.0,4.5,1.25,2.1),mgp=c(1.5,0.5,0) )

bp2<-plot(cut.e1e2,lambda,xlim=range.eig1,xlab="Eigenvector 1",ylab=expression(paste(bold(lambda),bold(" (mean(chi-squared values))"))),type="l")
points(cut.e1e2,lambda,col=lam.cols,pch=19,cex=2.0 )
leg2.txt<-c("Case=1870; Control=4306","Case=1833; Control=3960","Case=1821; Control=3748","Case=1788; Control=3488","Case=1708; Control=3088")
#legend(-0.052,0.035,leg.txt,col=c("blue","light blue","orange","red","purple"),pch=19)    #,bg="white" # e1 vs e2
points(vec["US_LAL-513",eig1],vec["US_LAL-513",eig2],pch=21,cex=2,col="magenta")

legend(0.002,1.45,leg2.txt,col=lam.cols,pch=19,cex=0.8)

for(i in 1:length(lambda)) {
abline(v=cut.e1e2[i],col=lam.cols[i],lwd=2) }

cut<-c(6sigma,keep2(e1e3),cut3, cut4
cut.e1e2<-c(-0.0505,-0.03, -0.01  ,-0.005,0.0)
lambda<-c(1.528005,1.284349,1.210552,1.174483,1.155952)
samples<-c( 6176,5793,5569, 5276, 4796 )
aff.samples<-c(1870,1833,1821,1788,1708)
#lam.cols<-c("magenta","red","orange","blue","green")
lam.cols<-c("green","blue","orange","magenta","red")

for(i in 1:length(lambda)) {
abline(v=cut.e1e2[i],col=lam.cols[i],lwd=2) }




 ###################### PLOT ALL
 
 
 ## Anova statistics for population differences along each eigenvector:
                                              p-value
             eigenvector_1_Control_Case_   2.11474e-16 +++
             eigenvector_2_Control_Case_      0.275883
             eigenvector_3_Control_Case_             0 +++
             eigenvector_4_Control_Case_     0.0103587 
             eigenvector_5_Control_Case_    1.3771e-12 +++
             eigenvector_6_Control_Case_    0.00278324 
             eigenvector_7_Control_Case_      0.206481 
             eigenvector_8_Control_Case_             0 +++
             eigenvector_9_Control_Case_      0.151565 
            eigenvector_10_Control_Case_   5.47718e-06 ***


# par(mfrow=c(2,1),font=2,font.lab=2,font.axis=2,mar=c(3.5,4.1,1.25,2.1),mgp=c(1.5,0.5,0) )
  par(mfrow=c(3,3),font=1,font.lab=2,font.axis=2,mgp=c(2,1,0),mar=c(5,3,4,2)+0.1)
 #  par(font=1,font.lab=1,font.axis=1,mgp=c(3.0,1,0),mar=c(5,4,4,2)+0.1) # default

vec[vec[,"origin"]=="us.con","color"]<-"blue"
vec[vec[,"origin"]=="wtccc","color"]<-"light blue"
vec[vec[,"origin"]=="hapmap","color"]<-"blue"
vec[vec[,"origin"]=="uk.aff","color"]<-"orange"
vec[vec[,"origin"]=="us.aff","color"]<-"red"
vec[vec[,"origin"]=="aus.aff","color"]<-"purple"

eig1<-"e.1"
eig2<-"e.2"
                    tapply(selected.ann[,"origin"],selected.ann[,"origin"],length)
    hapmap    UK.case    UK.conf UK.control    US.case    US.conf US.control 
       126          4         14          1         19        205         50 
# bp<-plot(c(-0.02,0.02),c(-0.025,0.0),xlab=eig1,ylab=eig2,main="Principle Components:Case+Control+HapMap")
  bp<-plot(range(vec[,eig1]),range(vec[,eig2]),xlab=eig1,ylab=eig2,main="EIGENSTRAT Principe Components")

#bp<-plot(c(-0.005,0.085),c(-0.035,0.095),xlab="Eigenvector 1",ylab=eig2,main="Principle Components:Case+Control+HapMap",las=2)
points(vec[,eig1],vec[,eig2],pch="+", col = vec[,"color"])  

 select<-vec[rownames(hapmap.eth[hapmap.eth[,"eth"]=="JPT",]),]
 select<-select[!is.na(select[,1]),]      # are NA's in the matrix since some CEUs already missing
 points(select[,eig1],select[,eig2],pch="+", col = "dark green")

 select<-vec[rownames(hapmap.eth[hapmap.eth[,"eth"]=="YRI",]),]
 select<-select[!is.na(select[,1]),]      # are NA's in the matrix since some CEUs already missing
 points(select[,eig1],select[,eig2],pch="+", col = "magenta")

 select<-vec[rownames(hapmap.eth[hapmap.eth[,"eth"]=="CEU",]),]      # just get out CEU hapmaps
 select<-select[!is.na(select[,1]),]      # are NA's in the matrix since some CEUs already missing
 points(select[,eig1],select[,eig2],pch=20, col = "black") 
 
  leg.txt<-c("icontrolDB","WTCCC control","UK affected","US affected","CHB","JPT","YRI","CEU")
  legend(0.05,0.07,leg.txt,col=c("blue","light blue","orange","red","green","dark green","magenta","black"),pch=20)

 ###################### PLOT ALL

 new.range1<-locator()   # right click and thes STOP
  new.range2<-locator()
 y=
m<-(new.range2$y-new.range1$y)/(new.range2$x-new.range1$x)
c<-new.range1$y - m*(new.range1$x)
abline(coef=c(c,m))
selected<-vec[vec[,eig2]> m*vec[,eig1] +c,]
#tapply(selected[,"origin"],selected[,"origin"],length)

 select<-vec[rownames(hapmap.eth[hapmap.eth[,"eth"]=="CEU",]),]      # just get out CEU hapmaps
 select<-select[!is.na(select[,1]),]      # are NA's in the matrix since some CEUs already missing
 points(select[,eig1],select[,eig2],pch=20, col = "black")









 #################### plot PLINK MDS DATA

file<-"merge_all_final_cSNP_no6pF.mds"
file<-"merge_all_final_cSNP_no6pF_mafF.mds"
file<-"merge_all_final_cSNP_no6pF_mafF3.mds"
file<-"merge_all_chr2F.mds"
file<-"merge_all_chr2F_noswaps.mds"
file<-"merge_subset_no6p_C.mds"
file<-"merge_subset_no6pALL_C.mds"


mds<-read.delim(file,header=T,sep="",fill=TRUE)
ori.size<-dim(mds)[2]-3
mds[,dim(mds)[2]+1]<-NA
mds[,dim(mds)[2]+1]<-NA
mds[,dim(mds)[2]+1]<-NA
mds[,dim(mds)[2]+1]<-NA
rownames(mds)<-mds[,1]
mds<-mds[,c(-1,-2,-3)]

colnames(mds)<-c(paste("e.",1:ori.size,sep=""),"status","origin","color","cluster")

mds[,"origin"]<-super[match(rownames(mds),rownames(super)),"origin"]
sum(is.na(mds[,"origin"]))     # check got all samples
mds[,"cluster"]<-cluster.type[match(rownames(mds),rownames(cluster.type)),"cluster"]
sum(!is.na(mds[,"cluster"])) # cluster type is for affected what cluster patrick used to genotype

affected<-grep(".aff",mds[,"origin"])
mds[affected,"status"]<-"Case"
mds[is.na(mds[,"status"]),"status"]<-"Control"

#mds[mds[,"status"]=="Case" & is.na(mds[,"cluster"]),]


sum(is.na(mds[,"status"]))

mds[mds[,"origin"]=="us.con","color"]<-"blue"
mds[mds[,"origin"]=="wtccc","color"]<-"light blue"
mds[mds[,"origin"]=="hapmap","color"]<-"blue"
mds[mds[,"origin"]=="uk.aff","color"]<-"orange"
mds[mds[,"origin"]=="us.aff","color"]<-"red"
mds[mds[,"origin"]=="aus.aff","color"]<-"purple"

sum(is.na(mds[,"color"]))

# pch=19: solid circle,
# pch=20: bullet (smaller circle),
# pch=21: circle,
# pch=22: square,
# pch=23: diamond,
# pch=24: triangle point-up,
# pch=25: triangle point down.

 par(mfrow=c(1,1),font=1,font.lab=2,font.axis=2,mgp=c(3.5,1,0),mar=c(5,5,4,2)+0.1)


#leg.txt<-c("UK Affected","US affected","AUS affected")

eig1<-"e.1"
eig2<-"e.2"

bp<-plot(range(mds[,eig1]),range(mds[,eig2]),xlab=eig1,ylab=eig2,main="PLINK Principle Components")

#bp<-plot(range(mds[,eig1]),range(mds[,eig2]),xlab="Eigenmdstor 1",ylab="Eigenmdstor 3",main="Principle Components: Case+Control")
#bp<-plot(range(mds[,eig1]),c(-0.0515,0.06),xlab=eig1,ylab=eig2,main="Principle Components:Affected only")

#points(sam_genes.pca$x[,1],sam_genes.pca$x[,2],col=colours[sam_genes.cl$cluster],pch=sam_genes.cl$cluster,cex=1.0,bg=colours[sam_genes.cl$cluster]) #colours and symbols
#text(mds[,eig1],mds[,eig2],label=mds[,"origin"],cex=0.5) #colours and symbols
points(mds[,eig1],mds[,eig2],pch="+", col = mds[,"color"])

mds.chr2<-mds
mds<-mds.chr2

select<-mds[,"origin"]=="us.aff"
points(mds[select,eig1],mds[select,eig2],pch="+", col = "red")

 select<-mds[rownames(hapmap.eth[hapmap.eth[,"eth"]=="JPT",]),]
  select<-select[!is.na(select[,1]),]      # are NA's in the matrix since some CEUs already missing
 points(select[,eig1],select[,eig2],pch="+", col = "dark green")

 select<-mds[rownames(hapmap.eth[hapmap.eth[,"eth"]=="YRI",]),]
  select<-select[!is.na(select[,1]),]      # are NA's in the matrix since some CEUs already missing
 points(select[,eig1],select[,eig2],pch="+", col = "magenta")
 
  select<-mds[rownames(hapmap.eth[hapmap.eth[,"eth"]=="CEU",]),]      # just get out CEU hapmaps
 select<-select[!is.na(select[,1]),]      # are NA's in the matrix since some CEUs already missing
 points(select[,eig1],select[,eig2],pch=20, col = "black") 
 
  leg.txt<-c("icontrolDB","WTCCC control","UK affected","US affected","CHB","JPT","YRI","CEU")
  legend(-1.5,0.15,leg.txt,col=c("blue","light blue","orange","red","green","dark green","magenta","black"),pch=20)


leg.txt<-c("US control","WTCCC control","UK affected","US affected","AUS affected")
 from<-"status"
 what<-"Control"

 from<-"origin"
 what<-"hapmap"
points(mds[mds[,from]==what,eig1],mds[mds[,from]==what,eig2],pch=19, col = "green")

select<-mds[,"e.1"]<0  & mds[,from]==what

select<-mds[,"e.1"]< -5
samples<-rownames(mds)[select]
mds<-mds[!select,]

type<-"uk.aff"
number<-100
temp<-mds[mds[,"origin"]==type,]
extra<-rownames(temp)[round(runif(number,1,dim(temp)[1]))]    # can get dupliactes 20.3 , 20.4 -> 20
samples<-c(samples,extra)
length(samples)
samples<-unique(samples)
set1<-data.frame( one=samples ,two=samples)
write.table(set1,"subset_mds_samples.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)


points(mds[select,eig1],mds[select,eig2],pch=21, col = "magenta")

> dim(mds)
[1] 6938   12
 > tapply(mds[mds[,"e.1"]>0,"status"],mds[mds[,"e.1"]>0,"status"],length)
   Case Control
    977    3434
> tapply(mds[mds[,"e.1"]<0,"status"],mds[mds[,"e.1"]<0,"status"],length)
   Case Control
   1091    1436
 tapply(mds[mds[,"e.1"]>0,"origin"],mds[mds[,"e.1"]>0,"origin"],length)
 aus.aff  hapmap  uk.aff  us.aff  us.con
     51     262     524     402    3172
 tapply(mds[mds[,"e.1"]<0,"origin"],mds[mds[,"e.1"]<0,"origin"],length)
 aus.aff  uk.aff  us.aff   wtccc
     18     492     581    1436
     
 # NOTE    clustering splits effected evenly
 # ALL US controls in one cluster , ALL WTCCC in another 
 # HAPMAP clusters with US contols
 select<-mds[,"e.1"]<0
 mds.set1<-mds[select,]
 dim(mds.set1)

 mds.set1[,"e.1"]
 mds.set1<-mds.set1[-1,] ## for lt 0 set

dim(mds.set1)
get.mds.set1<-data.frame( one=rownames(mds.set1) ,two=rownames(mds.set1))
write.table(get.mds.set1,"get_mds_set_lt0.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(get.mds.set1,"get_mds_set_gt0.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

US_LAL-513 -0.342521 2.4167000 -2.3811500 -0.02106980  0.00454281 #####strange one out of them all


 sample<-"UK1013"
 mds[sample,"e.1"]


























