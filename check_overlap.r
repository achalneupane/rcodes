

setwd("/media/UQCCG/GWAS/AS_Ichip_Study/HLA")
data<-read.table("Full_HLA_table_AS.txt",header=T,fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
samples<-read.delim("MBrown_Azores_HLA_Plate Layouts_20130924.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
ann<-read.delim("/media/UQCCG/GWAS/AS_Ichip_Study/Clinical_Data/UQCCG/clinical_data/allclinicialvariable_tasc_wttcc2.csv",header=T,sep=",",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
aogc<-read.table("/media/UQCCG/GWAS/AOGC/snp2hla/AOGC_IMPUTED.fam",header=F,fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")


dim(data)
dim(samples)

samples[1:5,]
ann[1:5,]
data[1:5,1:5]
aogc[1:5,]

samples$Sample.ID %in% ann$gentoype.ID
samples$Sample.ID %in% ann$BC.database
samples$Sample.ID %in% ann$phenotype.ID
the.missing<-samples$Sample.ID[!(samples$Sample.ID %in% data[,1])]

i<-1

for (i in 1:length(the.missing)){
found<-data[grep(the.missing[i],data[,1]),1]
if(length(found)>0){
a.set<-cbind(found,the.missing[i])
}
colnames(a.set)<-c("a","b")
if(i==1){set<-a.set}
set<-rbind(set,a.set)
}
set[1:5,]

posns<-match(data[,1],set[,1])
missing<-is.na(posns)
data[!missing,1]<-set[posns[!missing],2]

length(aogc$V1[!(aogc$V1 %in% data[,1])])
#934
#1002

data[grep("0004",data[,2]),2]

sum(aogc$V1=="AOGC-05-0004")

wanted<-c(samples$Sample.ID,aogc$V1)

posns<-match(data[,1],wanted)
missing<-is.na(posns)
data.sub<-data[!missing,]
dim(data.sub)

write.table(data.sub,file="HLA_Azores.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
