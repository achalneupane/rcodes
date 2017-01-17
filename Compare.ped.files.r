
work.dir<-
  setwd("/media/scratch2/Uveitis/Problem Fix")
files<-dir(getwd())

files[grepl(".ped$",files)]


pedFile1<-"ANTXR2_sanger_on_tasc_F.ped"
pedFile2<-"ANTXR2.tasc.ped"

pedFile1<-"All_on_sanger.ped"
pedFile2<-"WTCC2_on_sanger.ped" 

pedFile1<-"All_on_WTCC.ped"
pedFile2<-"WTCC2_on_WTCC.ped" 

pedFile1<-"All_on_WTCC.ped"
pedFile2<-"Qced_on_WTCC_F.ped"

pedFile1<-"ALL_antrx.ped"
pedFile2<-"BC_antrx.ped"

pedFile1<-"All_test.ped"
pedFile2<-"Qced_test_F.ped"

pedFile1<-"Ichip_FFCS_QCids_test.ped"
pedFile2<-"total_common_QCids_test.ped" 

ped1<-read.delim(pedFile1,header=F,skip=0,sep="",fill=TRUE,stringsAsFactors=FALSE)
ped2<-read.delim(pedFile2,header=F,skip=0,sep="",fill=TRUE,stringsAsFactors=FALSE)

map1<-read.delim(gsub(".ped$",".map",pedFile1),header=F,skip=0,sep="",fill=TRUE,stringsAsFactors=FALSE)
map2<-read.delim(gsub(".ped$",".map",pedFile2),header=F,skip=0,sep="",fill=TRUE,stringsAsFactors=FALSE)


ped1
ped2

map1
map2
identical(map1[,4],map2[,4])
map1[,4]==map2[,4]

posns<-match(ped1[,1],ped2[,1])
missing<-is.na(posns)
sum(missing)

ped1<-ped1[!missing,]
ped2<-ped2[posns[!missing],]

posns<-match(map1[,2],map2[,2])
posns

if(sum(posns != c(1:length(posns)))>0){  ### need to change order of snps
reorder<-as.numeric(matrix(data=c(posns*2,posns*2+1),nrow=2,ncol=length(posns),byrow=TRUE))
reorder<-c(1:6,reorder+5)
missing<-is.na(reorder)
ped1<-ped1[,!missing]
ped2<-ped2[,reorder[!missing]]
}

colnames(ped1)<-1:dim(ped1)[2]
colnames(ped2)<-1:dim(ped2)[2]

rownames(ped1)<-1:dim(ped1)[1]
rownames(ped2)<-1:dim(ped2)[1]

ped1[ped1[,5]==0,]<-ped2[ped1[,5]==0,] # 
ped2[ped2[,5]==0,]<-ped1[ped2[,5]==0,]

identical(ped1,ped2,num.eq=TRUE)  ## This fails is rownames and colnames not identical
dim(ped1)
dim(ped2)
ped1[1:5,]
ped2[1:5,]
test<-TRUE
for(i in 1:dim(ped1)[1]){
 print(sum(ped1[i,] != ped2[i,])==0)
 test<-sum(ped1[i,] != ped2[i,])==0 & test
 }
test
 print(paste("Is identical:",test,sep=" "))


for(i in 1:dim(ped1)[2]){
 print(sum(ped1[,i] != ped2[,i])==0)
 test<-sum(ped1[,i] != ped2[,i])==0 & test
 }
test
 print(paste("Is identical:",test,sep=" "))
