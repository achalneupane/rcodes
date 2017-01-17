As mentioned, ASAS have requested some additional ROC analyses.  I do not understand the reasons for their request but don’t think we can proceed without humouring them.  

 

T



 

I have gone back to them and confirmed that they won’t ask for additional analyses..

 

Thanks Paul.

 

Matt


setwd("/media/Bioinform-D/Research/ASAS/ASAS")
source("functions.R")
library(akima)
library(ROCR)
library(nlme)

## load clinical data
d <- read.table("asas_new_phenotypes_with_defnotAS.txt",header=T)
#d <- d[,c("ASAS_ID","HLAB27","Definite.AS","Def.Not.AS")]
table(d[,"ethnicity"])
table(d[,"Center"])


## 1-Asian
## 2-African
## 3-East Indian
## 4-Hispanic
## 5-Indigenous
## 6-Mixed
## 7-Turkish
## 8-Causian

#### old data from earily 2012;
	--bfile /home/acortes/Projects/Immunochip/ASAS/For_Dana/ASAS
	--extract confirmed_as_immunochip_loci.txt
	--recodeA
	--out ASAS


cd /media/Bioinform-D/Research/ASAS/ASAS/Sept12$
plink --bfile DATA_FOR_ASAS --recodeA --out ASAS.new
## load genotype data
g <- read.table("ASAS.raw",header=T)
snps.used<-colnames(g)[7:dim(g)[2]]
################################################################

### NEW DATA
snps.used.ori<-snps.used
g.ori<-g


g <- read.table("ASAS.new.raw",header=T)
#g <- g[,-c(1,3:6)]


snps.used<-colnames(g)[7:dim(g)[2]]

snps.used.ori %in% snps.used

tmp.id <- g$IID
tmp.id <- gsub("ASAS_","",tmp.id)
g$ASAS_ID <- tmp.id

new.status<-read.table("/media/Bioinform-D/Research/ASAS/ASAS/Sept12/training_data_MHC_snps.txt",header=T,stringsAsFactors=FALSE)
#ASAS.sample<-rep(FALSE,times=dim(new.status)[1])
ASAS.sample<-grepl("ASAS_",new.status$ID)
tmp.id <- new.status$ID
tmp.id <- gsub("ASAS_","",tmp.id)
new.status$ID <- tmp.id
new.status[1:5,]
new.status<-cbind(new.status,ASAS.sample)
colnames(new.status)<-paste("Ichip",colnames(new.status),sep=".")
new.status[1:5,]
g[1:5,]

### ALL ASAS 
sum(new.status$Ichip.PHEN[!ASAS.sample]==2 & new.status$Ichip.B27[!ASAS.sample]> 1.5 & !is.na(new.status$Ichip.B27[!ASAS.sample]))/sum(new.status$Ichip.PHEN[!ASAS.sample]==2)
### 14.6% cases not HLAB27
### 91 % controls not HLAB27
#### so this counts the reference allele
new.status$Ichip.B27<-abs(new.status$Ichip.B27-2) ### NOW it could the minor allele (that tags HLAB27)

sum(is.na(new.status$Ichip.B27)) #8
sum(new.status$Ichip.PHEN[!ASAS.sample]==2 & new.status$Ichip.B27[!ASAS.sample]> 0.5 & !is.na(new.status$Ichip.B27[!ASAS.sample] ))/
    sum(new.status$Ichip.PHEN[!ASAS.sample]==2 & !is.na(new.status$Ichip.B27[!ASAS.sample]))
14.9% cases 85% of cases are B27 positive
91% controls

> sum(new.status$Ichip.PHEN==1)
[1] 9644
> sum(new.status$Ichip.PHEN==2)
[1] 4775

g <- merge(g,new.status,by.x="ASAS_ID",by.y="Ichip.ID")
#############################################################





## load expression data
## Select genes that were found significan by Dana.
## e <- read.table("ExpressionData.txt")
## e <- e[,c("BZRAP1","COMT","FCGBP","FZD2")]

e <- read.table("ASAS_ALL_GENES_NORMALIZED.txt",stringsAsFactors=FALSE)
e<-t(e)
rownames(e)<-e[,1]
colnames(e)<-e[1,]
e<-e[,-1]
e<-e[-1,]
#
e[1:2,]

e<-as.data.frame(e,stringsAsFactors=FALSE)

gene.names<-colnames(e)


##get gene information
require("biomaRt")
genome.build<-"hg19"
## listMarts()  ; listMarts(host="july2009.archive.ensembl.org",archive=TRUE)
## ensembl=useMart("ensembl")  
## listDatasets(ensembl)  mart<-useMart("ensembl_mart_51",host="may2009.archive.ensembl.org",archive=TRUE); listDatasets(mart)
##  listFilters(mart)       mart<-useMart("ensembl_mart_51",host="may2009.archive.ensembl.org",dataset="hsapiens_gene_ensembl",archive=TRUE)

  
#if(!exists("mart")){
if(genome.build=="hg18"){mart<-useMart("ensembl_mart_51",host="may2009.archive.ensembl.org",dataset="hsapiens_gene_ensembl",archive=TRUE)
    }else if(genome.build=="hg19"){
      mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
    }else if(genome.build=="mm9"){
      mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl");print("here:")
 }
#} # is alread have mart then skip
  mart
## gene.names must be for each line, genes can be unordered 
genes<-unique(gene.names)
genes<-genes[!grepl("NONE",genes)]  #### NONE(dist=NONE exists)

#genes[1:10]
#USE the two beow for annotion with REFGENE
#gene.ann.wanted<-c("hgnc_symbol","gene_biotype","description")
#gene.desc<-getBM(attributes = gene.ann.wanted, filters = "hgnc_symbol", values=genes, mart = mart)
#USE the two beow for annotion with ENSGENE
if(genome.build=="hg19"){
gene.ann.wanted<-c("ensembl_gene_id","hgnc_symbol","gene_biotype","description")
}

if(genome.build=="hg18"){
gene.ann.wanted<-c("ensembl_gene_id","hgnc_symbol","biotype","description")
}
  
if(genome.build=="mm9"){
gene.ann.wanted<-c("ensembl_gene_id","mgi_symbol","gene_biotype","description")
}
  
gene.desc<-getBM(attributes = gene.ann.wanted, filters = "hgnc_symbol", values=genes, mart = mart)
#FIX THIS SO DOES THE CHROMOSOME REGION TO GET MiRNA and LOC names name use uscs ids if run with annovar
## gene.desc<-getBM( filters = "embl", values=genes, mart = mart)
## colnames(gene.desc)<- c("hgnc_symbol","gene_biotype","description")dbass3_name
gene.desc[1:2,]

######### get rid of dupicate entries with resplce to the fist column  dim(gene.desc)[1]==length(unique(gene.desc[,1]))
  dups<-duplicated(gene.desc[,2])
  gene.desc<-gene.desc[!dups,]
####################### match genes to another column #########################################

genes[!(genes %in% gene.desc[,2])]

write.table(gene.desc,file="table1.gene.list.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

#gene.desc.table<-match.cols.and.collapse(gene.names,"ensGene::gene",gene.desc,"ensembl_gene_id",gene.ann.wanted,"delimit")

## the above is much faster then below when doing delimited concatination compated to tapply
## list.element.lengths<-unlist(lapply(gene.names[["ensGene::gene"]],length))        
## indel.index<-rep(1:length(list.element.lengths),times=list.element.lengths)
## flat.gene.list<-unlist(the.genes[[annotation.labels[i]]])  
## gene.desc.table<-match.cols.and.collapse.fast(list.element.lengths,indel.index,flat.gene.list,gene.desc,"ensembl_gene_id",gene.ann.wanted,"delimit")


gene.desc[1:15,]

####################################################
##################################################3








tmp.id <- rownames(e)
tmp.id <- gsub("CtVals","",tmp.id)
tmp.id <- sapply(tmp.id,function(x) {
  out <- paste(substr(x,1,2),
               "-",
               substr(x,3,4),sep='')
  return(out)
})
e$ASAS_ID <- tmp.id
## ASAS_ID <- tmp.id
## e<-cbind(e,ASAS_ID)
## Merge datasets
data <- merge(d,g,by.x="ASAS_ID",by.y="ASAS_ID")
data <- merge(data,e,by.x="ASAS_ID",by.y="ASAS_ID")

data[1:5,]

## remove those without phen

num.cols<-c(gene.names,"HLAB27",snps.used) # num.cols<-c(gene.names,"HLAB27")
num.cols<-unique(num.cols)
num.cols
for( i in 1:length(num.cols)){
  ## a.line<-as.numeric(as.character( data[,num.cols[i]]))
  ## sum(is.na
  data[,num.cols[i]]<-as.numeric(as.character( data[,num.cols[i]]))
}

write.table(data,file="ASAS_all_data.used.in.new.analysis.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

data[,c("HLAB27","Ichip.B27","Ichip.RS2394250")]    #RS2394250 is the HLA_A tag snp


############### check for discordant HLAB27 typing
b27pos<-data[,"HLAB27"]==1
b27.ichip<-data[,"Ichip.B27"]>0.5

hla.conflict<-data[b27pos != b27.ichip,c("ASAS_ID","Center","sex","age", "ethnicity", "HLAB27","Ichip.B27","IID")]
sum(b27pos != b27.ichip)/dim(data)[1]
10% discordant ..... expect about 2
29 bad expected about 6
tapply(hla.conflict[,"Center"],hla.conflict[,"Center"],length)
 1  7 12 14 17 24 31 34 35 
 3  1  1  1  7  8  1  4  3

tapply(hla.conflict[,"ethnicity"],hla.conflict[,"ethnicity"],length)
hla.conflict

write.table(data[b27pos != b27.ichip,],file="HLA_conflict.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

dim(data)

################################################ wrong sex


bad.sex<-read.delim("wrong_sex.txt",header=F,skip=0,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
bad.sex<-unlist(bad.sex)

bad.samples<-c(bad.sex,as.character(hla.conflict[,"IID"]))
bad.samples<-unique(bad.samples) ## 36
bad.samples

###################################3

## data<-imm.data   all.data<-data
data<-all.data
colnames(data)[colnames(data)=="PHENOTYPE"]<-"PHEN"
#defin<-data[,"PHEN"]==1
#defin.NA<-is.na(data[,"PHEN"])
pos<-data[,"diagnosis"]==1
b27pos<-data[,"HLAB27"]==1
the.center<-data[,"Center"]==1

wanted<-data[,"ethnicity"] %in% c(8)
table(data[,"ethnicity"])

unique(data[,"Center"]) #  1  7 12 14 17 24 31 34 35
centers<-c(1,7,12,14,17,24,31,34,35)
ethnicitys<-c(1:9)
names(ethnicitys)<-c("Asian","African","East Indian","Hispanic","Indigenous","Mixed","Turkish","Causian","Unknown")

tapply(data[,"Center"],data[,"ethnicity"],length)



hla.aff<-matrix(data=NA,nrow=length(centers), ncol=length(ethnicitys))
hla.unaff<-matrix(data=NA,nrow=length(centers), ncol=length(ethnicitys))
pos.num<-matrix(data=NA,nrow=length(centers), ncol=length(ethnicitys))
notpos.num<-matrix(data=NA,nrow=length(centers), ncol=length(ethnicitys))

hla.pos.num<-matrix(data=NA,nrow=length(centers), ncol=length(ethnicitys))
hla.notpos.num<-matrix(data=NA,nrow=length(centers), ncol=length(ethnicitys))

summary<-matrix(data=NA,nrow=length(centers), ncol=length(ethnicitys))

pos<-data[,"diagnosis"]==1
b27pos<-data[,"HLAB27"]==1
n=3

i<-3
j<-1

sex<-data[,"sex"]
age<-data[,"age"]
defn.AS<-data[,"Definite.AS"]

male<-rep(NA,times=length(centers))
female<-rep(NA,times=length(centers))
mean.age<-rep(NA,times=length(centers))
sd.age<-rep(NA,times=length(centers))
range.age<-rep(NA,times=length(centers))
HLAB27.pos<-rep(NA,times=length(centers))
HLAB27.neg<-rep(NA,times=length(centers))
diagnosis.pos<-rep(NA,times=length(centers))
diagnosis.neg<-rep(NA,times=length(centers))
AS.NY.pos<-rep(NA,times=length(centers))
AS.NY.neg<-rep(NA,times=length(centers))

for (i in 1:length(centers)){
  the.center<-data[,"Center"]==centers[i]
  
for (j in 1:length(ethnicitys)){


the.ethnicitys<-data[,"ethnicity"]==ethnicitys[j]
hla.aff[i,j]<-100*sum(pos & b27pos & the.center & the.ethnicitys)/sum(pos  & the.center & the.ethnicitys)
hla.unaff[i,j]<-100*sum(!pos & b27pos & the.center & the.ethnicitys)/sum(!pos  & the.center & the.ethnicitys)

hla.pos.num[i,j]<-sum(pos & b27pos & the.center & the.ethnicitys)
hla.notpos.num[i,j]<-sum(!pos & b27pos & the.center & the.ethnicitys)

pos.num[i,j]<-sum(pos & the.center & the.ethnicitys)
notpos.num[i,j]<-sum(!pos & the.center & the.ethnicitys)

summary[i,j]<-paste(signif(hla.aff[i,j],digits=n),"% :",signif(hla.unaff[i,j],digits=n),"% [(",hla.pos.num[i,j],":",hla.notpos.num[i,j],")","(",pos.num[i,j],":",notpos.num[i,j],")]",paste="")
summary[i,j]<-gsub("NaN","-",summary[i,j])

}

#########Per centre stats
sex.c<-tapply(sex[the.center],sex[the.center],length)[c("1","2")]
  male[i]<-sex.c["1"]
  female[i]<-sex.c["2"]
  mean.age[i]<- mean(age[the.center],na.rm=TRUE)
  sd.age[i]<- sd(age[the.center],na.rm=TRUE)
  range.age[i]<- gsub(", ","-",toString(range(age[the.center],na.rm=TRUE)))
  HLAB27.pos[i]<-sum(b27pos[the.center],na.rm=TRUE)
  HLAB27.neg[i]<-sum(!b27pos[the.center],na.rm=TRUE)
  diagnosis.pos[i]<-sum(pos[the.center],na.rm=TRUE)
  diagnosis.neg[i]<-sum(!pos[the.center],na.rm=TRUE)     
  AS.NY.pos[i]<-sum(defn.AS[the.center],na.rm=TRUE)
  AS.NY.neg[i]<-sum(!defn.AS[the.center],na.rm=TRUE)                         
}


colnames(summary)<-paste(names(ethnicitys)," %HLAB27+ Aff : %HLAB27+ UnAff  Counts:[(HLAB27+ Aff:HLAB27+ UnAff) (Aff:Unaff)]",sep="")

summary<-cbind(centers,summary,male,female,mean.age,sd.age,range.age,HLAB27.pos,HLAB27.neg,diagnosis.pos, diagnosis.neg,AS.NY.pos,AS.NY.neg)
summary[1:5,]

tapply(data[,"Center"],data[,"Center"],length)


## write.table(summary,file="ASAS_HLA_status_Center_ethnicity_final.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

## write.table(summary,file="ASAS_HLA_status_Center_ethnicity_final_NEW2.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE) ### tatest OCT version 
summary <- read.delim("ASAS_HLA_status_Center_ethnicity_final_NEW2.txt",header=TRUE,sep="\t",skip=0,fill=TRUE,stringsAsFactors=FALSE)


colnames(hla.aff)<-names(ethnicitys)
colnames(hla.unaff)<-names(ethnicitys)

85.02687
8.570098

sum(pos & the.center & wanted)
sum(pos & b27pos & the.center & wanted)
sum(pos & b27pos & the.center & wanted)/sum(pos & the.center & wanted)

sum(!pos & b27pos & the.center & wanted)/sum(!pos & the.center & wanted)

sum(pos)
sum(b27pos)
sum(pos & b27pos)

sum(pos & b27pos & wanted)
sum(pos & wanted)

table(data[,"Center"])


100*sum(defin & !defin.NA & b27pos)/sum(defin & !defin.NA) 85%
100*sum(!defin & !defin.NA & b27pos)/sum(!defin & !defin.NA) 8.57

blank<-rep(FALSE,times=dim(data)[1])
for( i in 1:length(snps.used)){
  ## a.line<-as.numeric(as.character( data[,num.cols[i]]))
  is.blank<-is.na(data[,snps.used[i]])
  blank<-blank | is.blank
}
sum(blank)
sum(!defin & !defin.NA & blank)


#data.all<-data # data<-all.data
#write.table(data.all,file="Complete_normalised_dataset.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)




save.image("final.working.RData")
setwd("/media/Bioinform-D/Research/ASAS/ASAS/")


setwd("/media/scratch/ASAS/ASAS/")

source("functions.R")
library(akima)
library(ROCR)
library(nlme)


load("/media/Bioinform-D/Research/ASAS/ASAS/final.working.RData") # old computer
load("/media/scratch/ASAS/ASAS/final.working.RData")
load("/media/old-scratch/media/Bioinform-D/Research/ASAS/ASAS/final.working.RData")



hese are:

CASES	CONTROLS


IGAS AS CASES	ASAS IMAGING-NEGATIVE CONTROLS (those not meeting the ASAS-imaging positive criteria, irrespective if they meet other criteria)


ASAS IMAGING-POSITIVE CASES	ASAS IMAGING-NEGATIVE CONTROLS 


IGAS AS CASES	TOTAL ASAS CONTROLS (not meeting either imaging or non-imaging ASAS criteria)


ASAS IMAGING-POSITIVE CASES	TOTAL ASAS CONTROLS

############################################################## IMAGE POS? NEG

data<-data.all
#ori<- read.table("ASAS sample data.csv",header=T)
######################################################## SPA checks
ori<-read.delim("ASAS sample data.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
dim(ori)
dim(data)

colnames(ori)
colnames(data)
ori[,"ASAS_ID"]<-gsub(" ","",ori[,"ASAS_ID"])

posns<-match(as.character(data[,"ASAS_ID"]),ori[,"ASAS_ID"])
missing<-is.na(posns)
sum((missing))
posns

ori<-ori[posns[!missing],]

## sum(ori[,"ASAS_axial_spa_criteria"] != data[,"ASAS_axial_spa_criteria"]) #0
## sum(ori[,"diagnosis"] != data[,"diagnosis"]) #0
## sum(ori[,"HLAB27"] != data[,"HLAB27"]) #0
## sum(ori[,"MRI_SI_acute"] != data[,"MRI_SI_acute"]) #0

## pos<-data[,"ASAS_axial_spa_criteria"]==1 #137
## ori<-ori[pos,]
## data<-data[pos,]
## sum(pos)
## dim(ori)
## dim(data)

## risk<-data[,"HLAB27"]==1 | data[,"MRI_SI_acute"]==1
## sum(!risk) #140 # 16 not 

## strange<-pos & !pos
## sum(strange)

## risk<-data[,"HLAB27"]==1 | data[,"MRI_SI_acute"]==1 | data[,"Mod_NY_1"]==1
## sum(!risk) #140 # 16 not 
## ori[!risk,]

Mod_NY_1

################################ SET UP IMAGING ONLY :AXIAL SPA CRITERIA "Mod_NY_1" |    "MRI_SI_acute" AN ONE OTHER

colnames(ori)[c(39,44,46, 10, 13,  14,15, 22,23,  24,25,  26,27,28,29,41)]

test.1<-apply(ori[,c(39,10,13,14,15,22,23,24,25,26,27,28,29,41)],1,function(x) sum(as.numeric(x),na.rm=TRUE))
test.1<-test.1 >=1
sum(test.1)
length(test.1) #273 out of 286
test.1

colnames(ori)[c(39,10,13,14,15,22,23,24,25,26,27,28,29,41)]

colnames(ori)[c(44,46)]

spa.pos<-( ori[,44]==1 | ori[,46]==1 ) &  test.1
sum(spa.pos) # 92
table(ori[spa.pos,"HLAB27"]) #66%
##  0  1 
## 31 61
table(ori[ori[,"ASAS_axial_spa_criteria"]==1,"HLAB27"])
 0   1 
 31 106 

chk<-spa.pos & !ori[,"ASAS_axial_spa_criteria"]==1
sum(chk) #0

ori[chk,c(44,46, 10, 13, 14,15, 22,23,  24,25,  26,27,28,29,41)]
ori[chk,]
60/92

spa.pos.imaging<-spa.pos


spa.pos<-spa.pos.imaging
sum(spa.pos) # 92

######################################################################
###########image postive SPA
data[spa.pos,"PHEN"]<-1
data[!spa.pos,"PHEN"]<-0

###########image postive-controls SPA
data[spa.pos,"PHEN"]<-0
data[!spa.pos,"PHEN"]<-1
##################################

############## Image negative

colnames(ori)[c(13,26,27,28,29,41)]
colnames(ori)[c(14,15)]
colnames(ori)[c(22,23)]
colnames(ori)[c(24,25)]
colnames(ori)[c(39)]

test.A<-apply(ori[,c(13,26,27,28,29,41)],1,function(x) sum(as.numeric(x),na.rm=TRUE))

test.B<-apply(ori[,c(14,15)],1,function(x) sum(as.numeric(x),na.rm=TRUE))
test.B[test.B>=1]<-1
test.C<-apply(ori[,c(22,23)],1,function(x) sum(as.numeric(x),na.rm=TRUE))
test.C[test.C>=1]<-1
test.D<-apply(ori[,c(24,25)],1,function(x) sum(as.numeric(x),na.rm=TRUE))
test.D[test.D>=1]<-1




test.1<-test.A + test.B + test.C + test.D
test.1<-test.1 >=2
sum(test.1)
length(test.1) #176 out of 286
test.1
spa.pos<- ori[,39]==1 &  test.1

sum(spa.pos) # 87
table(ori[spa.pos,"HLAB27"]) #100%

chk<-spa.pos & !ori[,"ASAS_axial_spa_criteria"]==1
sum(chk) #1 07-04
ori[chk,]

spa.pos.no.imaging<-spa.pos
spa.pos<-spa.pos.no.imaging
#############################

spa.cases<-ori[,"ASAS_axial_spa_criteria"]==1
spa.controls<-ori[,"ASAS_axial_spa_criteria"]==0
sum(spa.pos.no.imaging & spa.pos.imaging) #48
sum(spa.pos.no.imaging & !spa.pos.imaging) #39
sum(!spa.pos.no.imaging & spa.pos.imaging)  #44
sum(!spa.pos.no.imaging & !spa.pos.imaging)

both<-spa.pos.no.imaging & spa.pos.imaging
either<-spa.pos.no.imaging | spa.pos.imaging
sum(either)
sum(both & ori[,39]==1) # chech 100
sum(spa.pos.imaging & ori[,39]==1 )
sum(!spa.pos.imaging)
sum(spa.cases)
sum(spa.controls)

sum(spa.pos.imaging & spa.controls  )

sum(ori[,"ASAS_axial_spa_criteria"] !=data[,"ASAS_axial_spa_criteria"]) # 0 so consistant

sum(data[,"ASAS_axial_spa_criteria"]==0)


sum(spa.pos.no.imaging) #87
sum(spa.pos.imaging  ) #92

# save.image("current.RData") load("current.RData")

###########image postive SPA
data[spa.pos,"PHEN"]<-1
data[!spa.pos,"PHEN"]<-0

chk<-(spa.pos.no.imaging | spa.pos.imaging) & !ori[,"ASAS_axial_spa_criteria"]==1
sum(chk) #1 07-04
ori[chk,]
sum(spa.pos.no.imaging | spa.pos.imaging)
sum(ori[,"ASAS_axial_spa_criteria"]==1)

chk<-ori[,"ASAS_axial_spa_criteria"]==1 & !(spa.pos.no.imaging | spa.pos.imaging) 
sum(chk) #1 07-04
ori[chk,]


data<-data.all
## snps.used <-colnames(g)[7:20]
snps.used<-colnames(data.all)[64:77]
snps.used
colnames(data)
missing.snps<-apply(data[,snps.used],2,function(x) sum(is.na(x)))
missing.snps[snps.used]
sort(missing.snps)
crap.snps<-names(missing.snps)[missing.snps>1]
good.snps<-snps.used[!(snps.used %in% crap.snps)]
## > (missing.snps)[missing.snps>2]
##     imm_1_25177701_T          rs4129267_T   X1kg_1_159746369_A 
##                    8                   18                    3 
##     imm_2_62421949_G    X1kg_3_27769911_0 ccc.5.96150086.T.C_T 
##                    3                  285                    9 
##    imm_5_158751323_A    imm_9_138373660_T   imm_10_101268715_T 
##                    6                   10                    7 
##   X1kg_14_87558574_C 
##                    5

gsub(", "," + ",toString(snps.used[!(snps.used %in% crap.snps)]))
gsub(", "," + ",toString(snps.used[!(snps.used %in% crap.snps)]))

############################
#just do genetic model
ASAS.sample<-grepl("ASAS_",g$FID)
data<-g[!ASAS.sample,]
colnames(data)[colnames(data)=="PHENOTYPE"]<-"PHEN"
colnames(data)[colnames(data)=="Ichip.B27"]<-"HLAB27"
data[data[,"PHEN"]==1,"PHEN"]<-0
data[data[,"PHEN"]==2,"PHEN"]<-1


sum(data[,"PHEN"]==0)
sum(data[,"PHEN"]==1)
sum(is.na(data[,"HLAB27"]))

data<-data[!is.na(data[,"HLAB27"]),]
tapply(data[,"PHEN"],data[,"PHEN"],length)
sum(is.na(data[,"PHEN"]))
colnames(data)
#########################################

## Phenotype based on Phil's Criteria  data<-all.data
b27.ichip<-data[,"Ichip.B27"]>0.5
sum(is.na(b27.ichip))
data[b27.ichip,"HLAB27"]<-1
data[!b27.ichip,"HLAB27"]<-0

phen <- rep(NA,nrow(data))
phen[data$Definite.AS] <- 1
phen[data$Def.Not.AS] <- 0
tapply(phen,phen,length)
data[,"PHEN"] <- phen
######################################


data[,"PHEN"]<-data[,"diagnosis"]


############### start here using ANY Axial Spa

data[,"PHEN"]<-data[,"ASAS_axial_spa_criteria"]

sum(data[,"PHEN"] != data[,"ASAS_axial_spa_criteria"])
sum(data[,"PHEN"] != data[,"diagnosis"])

sum(data[,"HLAB27"] != data[,"ASAS_axial_spa_criteria"])


 sum(data[,"diagnosis"] !=data[,"ASAS_axial_spa_criteria"])
68 discordatnt
sum(data[,"diagnosis"] ==data[,"ASAS_axial_spa_criteria"])
218 concordant

sum(data[,"diagnosis"] ==data[,"ASAS_axial_spa_criteria"])
####################################### 2013 test tp

data<-data.all
data[,"PHEN"]<-data[,"ASAS_axial_spa_criteria"]
data <- data[!is.na(data$PHEN),]
table(data[,"PHEN"])

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

############################################################################################################################################
##################################RESTART

data<-data.all
load("imm.data.RData")
colnames(imm.data)
table(imm.data[,"PHEN"])
table(data[,"PHEN"])

# data<-imm.data gentic only immunochip
#A ## befoe axial spa cases vs imm controls
colnames(imm.data) %in% colnames(data)

## data<-data[,colnames(imm.data)]
## data[1:5,]
## imm.data[1:5,]

## data<-data[data[,"PHEN"]==1,] # just the cases
## imm.data<-imm.data[imm.data[,"PHEN"]==1,] # just the cases
## imm.data[1:5,]

## table(imm.data[,"PHEN"])
## table(data[,"PHEN"])
## imm.data[,"PHEN"]<-0 # set imm AS cases to controls in this test

## data<-rbind(data,imm.data)
## table(data[,"PHEN"])

## data[data[,"PHEN"]==1,"PHEN"]<-2
## data[data[,"PHEN"]==0,"PHEN"]<-1
## data[data[,"PHEN"]==2,"PHEN"]<-0

## table(data[,"PHEN"])
## data<-data[data[,"PHEN"]==0,] # just the controls
## imm.data<-imm.data[imm.data[,"PHEN"]==0,] # just the controls
##################################################################
#1 ## befoe axial spa cases vs imm controls
colnames(imm.data) %in% colnames(data)

data<-data[,colnames(imm.data)]
data[1:5,]
imm.data[1:5,]

data<-data[data[,"PHEN"]==1,] # just the cases
imm.data<-imm.data[imm.data[,"PHEN"]==0,] # just the controls
data<-rbind(data,imm.data)
## data<-data[data[,"PHEN"]==0,] # just the controls
## imm.data<-imm.data[imm.data[,"PHEN"]==0,] # just the controls


sum(spa.pos.no.imaging & spa.pos.imaging) #48
sum(spa.pos.no.imaging & !spa.pos.imaging) #39
sum(!spa.pos.no.imaging & spa.pos.imaging)  #44
sum(spa.pos.imaging)
sum(spa.pos.no.imaging)
sum(spa.cases)
sum(spa.controls)
##################################################################
#A.1 ## ASAS cases vs  controls
colnames(imm.data) %in% colnames(data)

data<-data[,colnames(imm.data)]
data[1:5,]


data[spa.cases,"PHEN"]<-1 # just the cases
data[!spa.cases,"PHEN"]<-0
table(data[,"PHEN"])

##   0   1 
## 149 137

## data<-data[data[,"PHEN"]==0,] # just the controls
## imm.data<-imm.data[imm.data[,"PHEN"]==0,] # just the controls
##################################################################
#A.1 ## IGAS AS CASES	IGAS CONTROLS
load("imm.data.RData")
colnames(imm.data)



data<-data.all

table(imm.data[,"PHEN"])

colnames(imm.data) %in% colnames(data)

data<-imm.data 

table(data[,"PHEN"])


##################################################################


##################################################################
#A.21 ## IGAS AS CASES	ASAS IMAGING-POS CASES
load("imm.data.RData")
colnames(imm.data)



data<-data.all

table(imm.data[,"PHEN"])

colnames(imm.data) %in% colnames(data)

data<-data[,colnames(imm.data)]
data[1:5,]
imm.data[1:5,]
table(imm.data[,"PHEN"])

data<-data[spa.pos.imaging,] # ASAS IMAGING-POS cases
data[,"PHEN"]<-0
imm.data<-imm.data[imm.data[,"PHEN"]==1,] # just the cases
data<-rbind(data,imm.data)
table(data[,"PHEN"])
  0    1 
  92 4094 

##################################################################
#A.21 ## IGAS AS CASES	ASAS IMAGING-NEGATIVE CONTROLS
load("imm.data.RData")
colnames(imm.data)



data<-data.all

table(imm.data[,"PHEN"])

colnames(imm.data) %in% colnames(data)

data<-data[,colnames(imm.data)]
data[1:5,]
imm.data[1:5,]
table(imm.data[,"PHEN"])

data<-data[!spa.pos.no.imaging,] # ASAS IMAGING-NEGATIVE CONTROLS
data[,"PHEN"]<-0
imm.data<-imm.data[imm.data[,"PHEN"]==1,] # just the cases
data<-rbind(data,imm.data)
table(data[,"PHEN"])
  0    1 
 199 4094

##################################################################
#A.21 ## IGAS AS CASES	ASAS IMAGING-NEGATIVE cases
load("imm.data.RData")
colnames(imm.data)



data<-data.all

table(imm.data[,"PHEN"])

colnames(imm.data) %in% colnames(data)

data<-data[,colnames(imm.data)]
data[1:5,]
imm.data[1:5,]
table(imm.data[,"PHEN"])

data<-data[spa.pos.no.imaging,] # ASAS IMAGING-NEGATIVE cases
data[,"PHEN"]<-0
imm.data<-imm.data[imm.data[,"PHEN"]==1,] # just the cases
data<-rbind(data,imm.data)
table(data[,"PHEN"])
  ##  0    1 
  ## 87 4094


##################################################################
#A.21 ##ASAS IMAGING-POS CASES vs IGAS controls
load("imm.data.RData")
colnames(imm.data)



data<-data.all

table(imm.data[,"PHEN"])

colnames(imm.data) %in% colnames(data)

data<-data[,colnames(imm.data)]
data[1:5,]
imm.data[1:5,]
table(imm.data[,"PHEN"])

data<-data[spa.pos.imaging,] # ASAS IMAGING-POS cases
data[,"PHEN"]<-1
imm.data<-imm.data[imm.data[,"PHEN"]==0,] # just the controls
data<-rbind(data,imm.data)
table(data[,"PHEN"])
  1   0 
  92 4294 

##################################################################
#A.21 ##ASAS IMAGING-NEG CASES vs IGAS controls
load("imm.data.RData")
colnames(imm.data)



data<-data.all

table(imm.data[,"PHEN"])

colnames(imm.data) %in% colnames(data)

data<-data[,colnames(imm.data)]
data[1:5,]
imm.data[1:5,]
table(imm.data[,"PHEN"])

data<-data[spa.pos.no.imaging,] # ASAS IMAGING-POS cases
data[,"PHEN"]<-1
imm.data<-imm.data[imm.data[,"PHEN"]==0,] # just the controls
data<-rbind(data,imm.data)
table(data[,"PHEN"])
   0    1 
4294   87 
>


##################################################################
#A.21 ##ASAS IMAGING-NEG CONTROLS vs IGAS controls _2015 NEW
load("imm.data.RData")
colnames(imm.data)



data<-data.all

table(imm.data[,"PHEN"])

colnames(imm.data) %in% colnames(data)

data<-data[,colnames(imm.data)]
data[1:5,]
imm.data[1:5,]
table(imm.data[,"PHEN"])

data<-data[!spa.pos.no.imaging,] # ASAS IMAGING-POS cases
data[,"PHEN"]<-1
imm.data<-imm.data[imm.data[,"PHEN"]==0,] # just the controls
data<-rbind(data,imm.data)
table(data[,"PHEN"])
   0    1 
4294  199 
> 

##################################################################

##################################################################
#A.3 ##  ASAS IMAGING-POSITIVE CASES  ASAS IMAGING-NEGATIVE CONTROLS
colnames(imm.data) %in% colnames(data)

data<-data[,colnames(imm.data)]
data[1:5,]
sum(spa.pos.imaging)
sum(spa.pos.no.imaging)

data[spa.pos.imaging,"PHEN"]<-1 # just the cases
data[!spa.pos.imaging,"PHEN"]<-0
table(data[,"PHEN"])

##   0   1 
## 194  92

##   0   1 
## 192  91

## data<-data[data[,"PHEN"]==0,] # just the controls
## imm.data<-imm.data[imm.data[,"PHEN"]==0,] # just the controls






##################################################################
#A.4  ## IGAS AS CASES  TOTAL ASAS CONTROLS (not meeting
	either imaging or non-imaging ASAS
	criteria)
load("imm.data.RData")
colnames(imm.data)



data<-data.all

table(imm.data[,"PHEN"])

colnames(imm.data) %in% colnames(data)

data<-data[,colnames(imm.data)]
data[1:5,]
imm.data[1:5,]
table(imm.data[,"PHEN"])
table(data[,"PHEN"])

data<-data[spa.controls,] # ASAS IMAGING-NEGATIVE CONTROLS
data[,"PHEN"]<-0
imm.data<-imm.data[imm.data[,"PHEN"]==1,] # just the cases
data<-rbind(data,imm.data)
table(data[,"PHEN"])
 ##   0    1 
 ## 149 4094 


sum(spa.pos.imaging)
sum(spa.pos.no.imaging)
sum(spa.cases)
sum(spa.controls)


##################################################################
#A.5  ## ASAS IMAGING-POSITIVE CASES v TOTAL ASAS CONTROLS

load("imm.data.RData")
colnames(imm.data)



data<-data.all

table(imm.data[,"PHEN"])

colnames(imm.data) %in% colnames(data)

data<-data[,colnames(imm.data)]
data[1:5,]
sum(spa.pos.imaging)
sum(spa.controls)

sum(spa.pos.imaging & spa.controls)



data[spa.pos.imaging,"PHEN"]<-1 # just the cases
data[spa.controls,"PHEN"]<-0


data<-data[(spa.pos.imaging | spa.controls),]
table(data[,"PHEN"])
##   0   1 
## 149  92 


##################################################################
#A.6  ## ASAS IMAGING-NEGATIVE CASES v TOTAL ASAS CONTROLS

load("imm.data.RData")
colnames(imm.data)



data<-data.all

table(imm.data[,"PHEN"])

colnames(imm.data) %in% colnames(data)

data<-data[,colnames(imm.data)]
data[1:5,]
sum(spa.pos.no.imaging)
sum(spa.controls)

sum(spa.pos.no.imaging & spa.controls)



data[spa.pos.no.imaging,"PHEN"]<-1 # just the cases
data[spa.controls,"PHEN"]<-0


data<-data[(spa.pos.imaging | spa.controls),]
table(data[,"PHEN"])
 ##  0   1 
## 150  76

##################################################################


##################################################################
#1 ## befoe axial spa cases vs imm controls
colnames(imm.data) %in% colnames(data)

data<-data[,colnames(imm.data)]
data[1:5,]
imm.data[1:5,]

data<-data[data[,"PHEN"]==1,] # just the cases
imm.data<-imm.data[imm.data[,"PHEN"]==0,] # just the controls
data<-rbind(data,imm.data)
## data<-data[data[,"PHEN"]==0,] # just the controls
## imm.data<-imm.data[imm.data[,"PHEN"]==0,] # just the controls

#2 ######### non mhc-cases vs controls
table(data[,"PHEN"])
sum(data[,"PHEN"]==1 & data[,"HLAB27"]==0 ) # 613
cases<-data[ (data[,"PHEN"]==1 & data[,"HLAB27"]==0 ),]
controls<-data[ (data[,"PHEN"]==0 ),]
dim(cases)
dim(controls)
data<-rbind(cases,controls)

#3 ## befoe axial spa case HLA negative vs imm controls
colnames(imm.data) %in% colnames(data)

data<-data[,colnames(imm.data)]
data[1:5,]
imm.data[1:5,]

data<-data[data[,"PHEN"]==1 & data[,"HLAB27"]==0,] # just the cases
imm.data<-imm.data[imm.data[,"PHEN"]==0,] # just the controls
data<-rbind(data,imm.data)
## data<-data[data[,"PHEN"]==0,] # just the controls
## imm.data<-imm.data[imm.data[,"PHEN"]==0,] # just the controls

########## non mhc-cases vs controls
#4 ## befoe axial spa negative vs imm controls
colnames(imm.data) %in% colnames(data)

data<-data[,colnames(imm.data)]
data[1:5,]
imm.data[1:5,]

data<-data[data[,"PHEN"]==0,] # just the cases
data[,"PHEN"]<-1
imm.data<-imm.data[imm.data[,"PHEN"]==0,] # just the controls
data<-rbind(data,imm.data)
## data<-data[data[,"PHEN"]==0,] # just the controls
## imm.data<-imm.data[imm.data[,"PHEN"]==0,] # just the controls

########## non mhc-cases vs controls
#5 ## Spa positive and MRI positive  
"MRI_SI_acute"

colnames(imm.data) %in% colnames(data)
sum(data[,"PHEN"]==1 & data[,"MRI_SI_acute"]==1) #55  maybe chosen without B27
##                sum(data[,"PHEN"]==1 & data[,"MRI_SI_acute"]==1 & data[,"HLAB27"]==0) #15
##                sum(data[,"PHEN"]==1 & data[,"MRI_SI_acute"]==1 & data[,"HLAB27"]==1) #40 ## 72% B27+
## sum(data[,"PHEN"]==1 & data[,"MRI_SI_acute"]==0) #82 chosen using B27
##                sum(data[,"PHEN"]==1 & data[,"MRI_SI_acute"]==0 & data[,"HLAB27"]==0) #16
##                sum(data[,"PHEN"]==1 & data[,"MRI_SI_acute"]==0 & data[,"HLAB27"]==1) #66  ##80% B27+
## data<-data[data[,"PHEN"]==1 & data[,"HLAB27"]==0,] # just the cases


data<-data[(data[,"PHEN"]==1 & data[,"MRI_SI_acute"]==1),colnames(imm.data)] ## subset ROWS AND Columns
data[1:5,]
imm.data[1:5,]

imm.data<-imm.data[imm.data[,"PHEN"]==0,] # just the controls

table(data[,"PHEN"])
table(imm.data[,"PHEN"])
data<-rbind(data,imm.data)
table(data[,"PHEN"])
########## non mhc-cases vs controls
#6 ## Spa positive and MRI negative  
"MRI_SI_acute"

colnames(imm.data) %in% colnames(data)
sum(data[,"PHEN"]==1 & data[,"MRI_SI_acute"]==0) #82  maybe chosen without B27
##                sum(data[,"PHEN"]==1 & data[,"MRI_SI_acute"]==1 & data[,"HLAB27"]==0) #15
##                sum(data[,"PHEN"]==1 & data[,"MRI_SI_acute"]==1 & data[,"HLAB27"]==1) #40 ## 72% B27+
## sum(data[,"PHEN"]==1 & data[,"MRI_SI_acute"]==0) #82 chosen using B27
##                sum(data[,"PHEN"]==1 & data[,"MRI_SI_acute"]==0 & data[,"HLAB27"]==0) #16
##                sum(data[,"PHEN"]==1 & data[,"MRI_SI_acute"]==0 & data[,"HLAB27"]==1) #66  ##80% B27+
## data<-data[data[,"PHEN"]==1 & data[,"HLAB27"]==0,] # just the cases


data<-data[(data[,"PHEN"]==1 & data[,"MRI_SI_acute"]==0),colnames(imm.data)] ## subset ROWS AND Columns
data[1:5,]
imm.data[1:5,]



imm.data<-imm.data[imm.data[,"PHEN"]==0,] # just the controls

table(data[,"PHEN"])
table(imm.data[,"PHEN"])
data<-rbind(data,imm.data)
table(data[,"PHEN"])

################# Axal SPA and immunochop
################# Axal SPA and immunochop













table(data[,"PHEN"])

table(data[,"PHEN"])

   0    1 
4294  137


with no missing
> tapply(data[,"PHEN"],data[,"PHEN"],length)
   0    1 
4294  123
table(data[data[,"ASAS_axial_spa_criteria"]==1,"diagnosis"])
table(data[data[,"ASAS_axial_spa_criteria"]==0,"diagnosis"])

table(data[ (data[,"ASAS_axial_spa_criteria"]==1 & data[,"diagnosis"]==0 ),"HLAB27"])
table(data[ (data[,"ASAS_axial_spa_criteria"]==0 & data[,"diagnosis"]==1 ),"HLAB27"])

is.na(data[data[,"ASAS_axial_spa_criteria"]==1,"HLAB27"])
t <- chisq.test(make.mat,correct=TRUE)

table(data[data[,"Definite.AS"]==0,"HLAB27"])
  0   1 
142  85
table(data[data[,"Def.Not.AS"]==1,"HLAB27"])

table(data[data[,"PHEN"]==0,"HLAB27"])
is.na(data[data[,"PHEN"]==1,"HLAB27"])


#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
###############################################################################FIX SNPS

snps.used<-colnames(data)[2:14]
colnames(data)
snps.used[!(snps.used %in% colnames(data))] ## missing snps

snps.used<-snps.used[snps.used %in% colnames(data)]
snps.used

missing.snps<-apply(data[,snps.used],2,function(x) sum(is.na(x)))
missing.snps
sort(missing.snps)
crap.snps<-names(missing.snps)[missing.snps>1]
good.snps<-snps.used[!(snps.used %in% crap.snps)]
data<-data[,!(colnames(data) %in% crap.snps)]

crap.snps
good.snps
snps.used

#################################################
## ## use to subset by ethnicity or 
## wanted<-data[,"ethnicity"] %in% c(8) # CEU 
## wanted.c<-data[,"Center"] %in% c(1) # Asian
## sum(wanted & wanted.c)
## data<-data[wanted & wanted.c ,]
## table(data[,"PHEN"])
## data<-data[wanted,] # keep only core ethnicity

## ## use to sex swap samples for HLA status
## bad.sex
## bad.samples
## wanted<-!(data[,"IID"] %in% bad.samples)
## sum(wanted)
## data<-data[wanted,]
## ###############################################

## tapply(data[,"HLAB27"],data[,"PHEN"],length)
## sum(data[data[,"diagnosis"]==0,"HLAB27"])
## #####################################
## sum(is.na(data[,"HLAB27"]))


#################################################################################################################
#################################################################################################################
#################################################################################################################
################################################################################################################
#####################################################################
############ Genetic clean up 
missing.a.value<-apply(data[,c("HLAB27",good.snps)],1,function(x) sum(is.na(x)))
missing.a.value<-missing.a.value>0
sum(missing.a.value) # 9 # 14 for spa genetic 
data[missing.a.value,]
data<-data[!missing.a.value,]


tapply(data[,"PHEN"],data[,"PHEN"],length)
dim(data)
data <- data[!is.na(data$PHEN),]
dim(data)
table(data$PHEN)  ## must be 0 and 1 not (1 and 2 else routines below give problems)
data[1:5,]
##  0  1 
## 63 59

#### Diagnosis
##   0   1 
## 123 162

#### Diagnosis no bad samples
##   0   1 
## 105 144 

## Split data for CV nchooseK
Nsplits <- 300

## Nsplits=dim(data)[1]
## (1-1/Nsplits)*Nsplits
## data.list <- get_N_data_splits(data=data,trai.p=(1-1/Nsplits),N=Nsplits,phen.col='PHEN')


data.list <- get_N_data_splits(data=data,trai.p=.75,N=Nsplits,phen.col='PHEN') # using 0.7 as a slip made no difference
## highter numbers inclde more samples in training less in test

## Model to try
names(data.list[[1]])
data.list[[1]]$train[,"PHEN"]
sum(data.list[[1]]$train[,"PHEN"])
sum(data.list[[1]]$test[,"PHEN"])
#formula <- PHEN ~  HLAB27 + BZRAP1 +  EDG8 + GPR128 # originally
## USe SignifGroupsxls Yellow tagged  < 0.01
#formula <- PHEN ~  BZRAP1 +  EDG8 + GPR128  # originally and center # 2A "3 genes"

######################### Sept 4 2013 
formula <- PHEN ~ PSMG1 + ANTXR2 + SIRPA + CTNNB1 + TSC22D3 + SH3BP5L + INSIG1 + CREBZF + P2RY13 + ACTB + HLADRB5 + GIMAP4 + CD83 + SORT1 +SNX5 + TNFAIP8L2 + GZMK + TNFAIP3 + CXCR4 + DRAM + MOAP1 + BZRAP1 + RGS1 + NCF4 + TGFBR3 + TM9SF1 + HPRT1 + NAGK + RPL32 + GPR128 + HMGB1 + BUD13 +  FCGR1B ## Axial Spa plus HLAB27 genes

formula <- PHEN ~ PSMG1 + TNFAIP8L2 + INSIG1 ## Cau

formula <- PHEN ~  BZRAP1 + GPR128 + HPRT1 + BUD13 + CD83 + FCGR1B + CTNNB1 + ACTB + SNX5 + INSIG1 # in HLAB27

###


formula <- PHEN ~  HLAB27 + PSMG1 + ANTXR2 + SIRPA + CTNNB1 + TSC22D3 + SH3BP5L + INSIG1 + CREBZF + P2RY13 + ACTB + HLADRB5 + GIMAP4 + CD83 + SORT1 +SNX5 + TNFAIP8L2 + GZMK + TNFAIP3 + CXCR4 + DRAM + MOAP1 + BZRAP1 + RGS1 + NCF4 + TGFBR3 + TM9SF1 + HPRT1 + NAGK + RPL32 + GPR128 + HMGB1 + BUD13 +  FCGR1B ## Axial Spa plus HLAB27 genes

formula <- PHEN ~  HLAB27 + PSMG1 + TNFAIP8L2 + INSIG1 ## Cau

formula <- PHEN ~   HLAB27 + BZRAP1 + GPR128 + HPRT1 + BUD13 + CD83 + FCGR1B + CTNNB1 + ACTB + SNX5 + INSIG1 # in HLAB27


formula <- PHEN ~   HLAB27


formula <- PHEN ~  BZRAP1 +  EDG8 
 
#formula <- PHEN ~  BZRAP1 +  EDG8 + GPR128 + HLADRB5 + MOAP1 + PSMG1  #2B  "6 genes"
formula <- PHEN ~  BZRAP1 +  EDG8  + HLADRB5 + MOAP1 + PSMG1  #2B  "5 genes"

formula <- PHEN ~  BZRAP1 +  EDG8 + GPR128 + HLADRB5 + MOAP1 + PSMG1 + CD83 + NDUFB4 + RFWD2  #2C  "9 genes" 2D for CEU only


#formula <- PHEN ~  HLAB27 + BZRAP1 +  EDG8 + GPR128 # 3A
formula <- PHEN ~  HLAB27 + BZRAP1 +  EDG8 # 3A
#formula <- PHEN ~  HLAB27 + BZRAP1 +  EDG8 + GPR128 + HLADRB5 + MOAP1 + PSMG1 # 3B
formula <- PHEN ~  HLAB27 + BZRAP1 +  EDG8 + HLADRB5 + MOAP1 + PSMG1 # 3B
formula <- PHEN ~  HLAB27 + BZRAP1 +  EDG8 + GPR128 + HLADRB5 + MOAP1 + PSMG1 + CD83 + NDUFB4 + RFWD2 #3C # 3D for CEU only

####### Core Group use with NYAS group
formula <- PHEN ~  BZRAP1 + COMT + FCGBP + FZD2 + INSIG1 + MOAP1
formula <- PHEN ~ HLAB27 + BZRAP1 + COMT + FCGBP + FZD2 + INSIG1 + MOAP1 ## 0.70 gues to 0.72 with + TNFAIP3

formula <- PHEN ~ HLAB27


### genetic only tests current<-toString(good.snps)
toString(good.snps)
current==toString(good.snps)
length(good.snps)

formula <- PHEN ~ HLAB27 + imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A + imm_2_62404976_A + imm_2_102030060_A + rs11098964_A + imm_5_40526547_A + imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A ### for A.2 & A.3 & A.4 & A.5

formula <- PHEN ~ HLAB27 + imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A + imm_2_62404976_A + imm_2_102030060_A + rs11098964_A + imm_5_40526547_A + imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A ### for A.1 ASAS cases vs controls NOv 15th

formula <- PHEN ~ HLAB27 + imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A + imm_2_62404976_A + imm_2_102030060_A + rs11098964_A + imm_5_40526547_A + imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A ### for MY axial spa positive  image negative  vs controls

formula <- PHEN ~ HLAB27 + imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A + imm_2_62404976_A + imm_2_102030060_A + rs11098964_A + imm_5_40526547_A + imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A ### for MY axial spa not negative  vs controls

formula <- PHEN ~ HLAB27 + imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A + imm_2_62404976_A + imm_2_102030060_A + rs11098964_A + imm_5_40526547_A + imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A ### for MY axial spa positive  image pos  vs controls

formula <- PHEN ~ HLAB27 + imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A + imm_2_62404976_A + imm_2_102030060_A + rs11098964_A + imm_5_40526547_A + ccc.5.96150086.T.C_A + imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A ### for axial spa MRI+ vs controls

formula <- PHEN ~ HLAB27 + imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A + imm_2_62404976_A + imm_2_102030060_A + rs11098964_A + imm_5_40526547_A + imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A ### for axial spa MRI- vs controls

formula <- PHEN ~ HLAB27 + imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A + imm_2_62404976_A + imm_2_102030060_A + rs11098964_A + imm_5_40526547_A + imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A ### for axial spa pos no imaging - vs ASAS controls






formula <- PHEN ~ HLAB27 + imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A + imm_2_62404976_A + imm_2_102030060_A + rs11098964_A + imm_5_40526547_A + imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A ### for axial spa and immunochip

formula <- PHEN ~ imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A + imm_2_62404976_A + imm_2_102030060_A + rs11098964_A + imm_5_40526547_A + imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A ### for axial spa and immunochip

formula <- PHEN ~ HLAB27 +imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A + imm_2_62404976_A + imm_2_102030060_A + rs11098964_A + imm_5_40526547_A + ccc.5.96150086.T.C_A + imm_5_158751323_A + imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A + imm_21_39391390_G  ### immunochip


formula <- PHEN ~ imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A + imm_2_62404976_A + imm_2_102030060_A + rs11098964_A + imm_5_40526547_A + ccc.5.96150086.T.C_A + imm_5_158751323_A + imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A + imm_21_39391390_G
## 6 diagnosis genes + geneitic
## formula <- PHEN ~  HLAB27 + BZRAP1 +  EDG8 + GPR128 + HLADRB5 + MOAP1 + PSMG1 +  imm_1_67478546_A + imm_1_67520024_A + imm_1_159121069_T + imm_1_199144473_G + rs12758027_C + imm_2_102013732_A + ccc.2.102135805.G.T_G + imm_2_181660357_A + imm_2_181756697_C + X1kg_2_241212412_A + imm_5_35917200_A + imm_5_40560617_A + imm_5_96137839_G + imm_5_96278345_C + imm_5_158759370_A + rs17765610_G + imm_6_91047852_A + rs2402752_C + imm_10_80730323_A + imm_12_6317038_A + X1kg_12_6361386_A + imm_12_110346958_T + imm_16_28236248_A + imm_17_23120724_A + imm_17_23172294_G + rs9901869_G + imm_19_10386181_A + imm_19_10486067_A + imm_21_39388614_A + imm_21_44440169_G + imm_22_20286653_G

formula <- PHEN ~  HLAB27 + BZRAP1 +  EDG8 + GPR128 + HLADRB5 + MOAP1 + PSMG1 + CD83 + NDUFB4 + RFWD2 +  imm_1_67478546_A + imm_1_67520024_A + imm_1_159121069_T + imm_1_199144473_G + rs12758027_C + imm_2_102013732_A + ccc.2.102135805.G.T_G + imm_2_181660357_A + imm_2_181756697_C + X1kg_2_241212412_A + imm_5_35917200_A + imm_5_40560617_A + imm_5_96137839_G + imm_5_96278345_C + imm_5_158759370_A + rs17765610_G + imm_6_91047852_A + rs2402752_C + imm_10_80730323_A + imm_12_6317038_A + X1kg_12_6361386_A + imm_12_110346958_T + imm_16_28236248_A + imm_17_23120724_A + imm_17_23172294_G + rs9901869_G + imm_19_10386181_A + imm_19_10486067_A + imm_21_39388614_A + imm_21_44440169_G + imm_22_20286653_G

####
## formula <- PHEN ~  HLAB27 + BZRAP1 +  EDG8 + GPR128 + HLADRB5 + MOAP1 + PSMG1 + imm_2_102013732_A + imm_2_181756697_C  + imm_5_96137839_G  + imm_22_20286653_G + imm_5_158759370_A + imm_21_39388614_A + imm_1_67520024_A + imm_17_23172294_G + imm_1_159121069_T + imm_1_199144473_G

formula <- PHEN ~   HLAB27 + BZRAP1 +  EDG8 + GPR128 + HLADRB5 + MOAP1 + PSMG1 + CD83 + NDUFB4 + RFWD2 + imm_2_102013732_A + imm_2_181756697_C  + imm_5_96137839_G  + imm_22_20286653_G + imm_5_158759370_A + imm_21_39388614_A + imm_1_67520024_A + imm_17_23172294_G + imm_1_159121069_T + imm_1_199144473_G 
##########################

# NYAS 6 core genes + geneitic
formula <- PHEN ~  HLAB27 + BZRAP1 + COMT + FCGBP + FZD2 + INSIG1 + MOAP1 +  imm_1_67478546_A + imm_1_67520024_A + imm_1_159121069_T + imm_1_199144473_G + rs12758027_C + imm_2_102013732_A + ccc.2.102135805.G.T_G + imm_2_181660357_A + imm_2_181756697_C + X1kg_2_241212412_A + imm_5_35917200_A + imm_5_40560617_A + imm_5_96137839_G + imm_5_96278345_C + imm_5_158759370_A + rs17765610_G + imm_6_91047852_A + rs2402752_C + imm_10_80730323_A + imm_12_6317038_A + X1kg_12_6361386_A + imm_12_110346958_T + imm_16_28236248_A + imm_17_23120724_A + imm_17_23172294_G + rs9901869_G + imm_19_10386181_A + imm_19_10486067_A + imm_21_39388614_A + imm_21_44440169_G + imm_22_20286653_G

####
formula <- PHEN ~  HLAB27 + BZRAP1 + COMT + FCGBP + FZD2 + INSIG1 + MOAP1 + imm_19_10486067_A + imm_12_110346958_T  + rs17765610_G  + imm_2_102013732_A + ccc.2.102135805.G.T_G + X1kg_2_241212412_A + imm_6_91047852_A  + imm_5_40560617_A + rs12758027_C + rs9901869_G



## formula <- PHEN ~  BZRAP1 +  EDG8 + GPR128 + CD83 + NDUFB4 + RFWD2
##  formula <- PHEN ~ HLAB27 + BZRAP1 +  EDG8 + GPR128 + HLADRB5 + MOAP1 + PSMG1 + CD83 + NDUFB4 + RFWD2

## ## Use < 0.05  
## formula <- PHEN ~  BZRAP1 + CD83 + EDG8 + GPR128 +  MOAP1 + NDUFB4 + RFWD2 + RGS1
## formula <- PHEN ~  HLAB27 + BZRAP1 + CD83 + EDG8 + GPR128 +  MOAP1 + NDUFB4 + RFWD2 + RGS1  
## ## Use < 0.05  + whites
## formula <- PHEN ~  BZRAP1 + CD83 + EDG8 + GPR128 +  MOAP1 + NDUFB4 + RFWD2 + RGS1  + RNASE6 + FZD2 + FOSB
## ## whites only 
## formula <- PHEN ~  BZRAP1 + CD83 + EDG8 + GPR128 +  MOAP1 + NDUFB4  + RGS1  + RNASE6 + FZD2 + FOSB
## #data<- data[data[,"ethnicity"]==8,] # data<-all.data
## dim(data)
## #formula <- PHEN ~  HLAB27 + BZRAP1 + CD83 + EDG8 + GPR128 +  MOAP1 + NDUFB4  + RGS1  + RNASE6 + FZD2 + FOSB

## ### 5 gene model for diagnosis
## formula <- PHEN ~  HLAB27 + BZRAP1 +  EDG8 + GPR128 + CD83  + RNASE6

 
## ###########
## #diagnosis with HLAB27 ststus 
## formula <- PHEN ~  HLAB27+ ANTXR2 + CTNNB1 + DDIT4 + CXCR4 + GABARAPL1 + GIMAP4 + GZMK + HLADRB5 + MOAP1 + PSMG1 + PTRH2 + SH3BP5L + SIRPA + TGFBR3 + TSC22D3
## formula <- PHEN ~  ANTXR2 + CTNNB1 + DDIT4 + CXCR4 + GABARAPL1 + GIMAP4 + GZMK + HLADRB5 + MOAP1 + PSMG1 + PTRH2 + SH3BP5L + SIRPA + TGFBR3 + TSC22D3
## formula <- PHEN ~  ANTXR2 + HLADRB5 + MOAP1 + PSMG1  + SH3BP5L + SIRPA  + TSC22D3
## formula <- PHEN ~  HLAB27+ ANTXR2 + HLADRB5 + MOAP1 + PSMG1  + SH3BP5L + SIRPA  + TSC22D3


## ### from HALB27 controlled data
## formula <- PHEN ~  HLAB27 + TNFAIP3 + SH3BP5L + PTRH2 + PSMG1 + MOAP1 + INSIG1 + IL23R + HLADRB5 + ERAP2 + SNX5 + RFWD2 + DUSP2 + BATF + EDG8 + DDIT4 + GPR128
## formula <- PHEN ~  TNFAIP3 + SH3BP5L + PTRH2 + PSMG1 + MOAP1 + INSIG1 + IL23R + HLADRB5 + ERAP2 + SNX5 + RFWD2 + DUSP2 + BATF + EDG8 + DDIT4 + GPR128

## ## SNX5 + RFWD2 + DUSP2 + BATF  # Asian
## ## EDG8 + DDIT4 + GPR128  #HISPANIC

## ####### Core Group
## formula <- PHEN ~  BZRAP1 + COMT + FCGBP + FZD2 + INSIG1 + MOAP1
## formula <- PHEN ~ HLAB27 + BZRAP1 + COMT + FCGBP + FZD2 + INSIG1 + MOAP1 + TNFAIP3
## ######################################


## ############################# model choice: ### just run the linear fit select the top 5 preformers  and select the best ones (5-7) to make the final model)
## ###do le lot genes only
## formula <- PHEN ~  TNFAIP3 + SH3BP5L + PTRH2 + PSMG1 + MOAP1 + INSIG1 + IL23R + HLADRB5 + ERAP2 + SNX5 + RFWD2 + DUSP2 + BATF + EDG8 + DDIT4 + GPR128 + ERAP1 + COMT + FCGBP + FZD2 + CD83  + RNASE6 + SIRPA  + TSC22D3 + CXCR4 + GABARAPL1  + GZMK + TGFBR3
## ### select the best
## formula <- PHEN ~ EDG8 + MOAP1 + CD83 + SH3BP5L + HLADRB5 + COMT + FZD2 + GABARAPL1 + TNFAIP3 # AUC=0.62 so no better than other methods
## ###do le lot genes only HLA
## formula <- PHEN ~  HLAB27 + TNFAIP3 + SH3BP5L + PTRH2 + PSMG1 + MOAP1 + INSIG1 + IL23R + HLADRB5 + ERAP2 + SNX5 + RFWD2 + DUSP2 + BATF + EDG8 + DDIT4 + GPR128 + ERAP1 + COMT + FCGBP + FZD2 + CD83  + RNASE6 + SIRPA  + TSC22D3 + CXCR4 + GABARAPL1  + GZMK + TGFBR3
## ### select the best
## formula <- PHEN ~ HLAB27 + MOAP1 + EDG8 + TNFAIP3 + CD83 + GPR128 + FZD2 + COMT + TGFBR3 + PTRH2 # AUC=0.62 so no better than other methods
## ############################################

## ### Frome Gethin from NYAS
## formula <- PHEN ~ FCGBP+ NDUFB4 + ID3+ BZRAP1 + E2F2+ ANTXR2
## formula <- PHEN ~ HLAB27 + FCGBP+ NDUFB4 + ID3+ BZRAP1 + E2F2+ ANTXR2 # 0.73
## #E2F2+ ANTXR2 # from HLAB27 pos
## formula <- PHEN ~ HLAB27 + FCGBP+ NDUFB4 + ID3+ BZRAP1 + E2F2+ ANTXR2+ BZRAP1 + COMT  + FZD2 + INSIG1 + MOAP1 + TNFAIP3 ## No logical reason to try this and gives 0.70

## ## Use < 0.05  + whites + new genetic for diagnosis and HLAB27
## formula <- PHEN ~ HLAB27 + BZRAP1 + CD83 + EDG8 + GPR128 +  MOAP1 + NDUFB4 + RFWD2 + RGS1 +  imm_1_67478546_A + imm_1_67520024_A + imm_1_159121069_T + imm_1_199144473_G + rs12758027_C + imm_2_102013732_A + ccc.2.102135805.G.T_G + imm_2_181660357_A + imm_2_181756697_C + X1kg_2_241212412_A + imm_5_35917200_A + imm_5_40560617_A + imm_5_96137839_G + imm_5_96278345_C + imm_5_158759370_A + rs17765610_G + imm_6_91047852_A + rs2402752_C + imm_10_80730323_A + imm_12_6317038_A + X1kg_12_6361386_A + imm_12_110346958_T + imm_16_28236248_A + imm_17_23120724_A + imm_17_23172294_G + rs9901869_G + imm_19_10386181_A + imm_19_10486067_A + imm_21_39388614_A + imm_21_44440169_G + imm_22_20286653_G


## formula <- PHEN ~ HLAB27 + BZRAP1 + COMT + FCGBP + FZD2 + INSIG1 + MOAP1 + TNFAIP3 +  imm_1_67478546_A + imm_1_67520024_A + imm_1_159121069_T + imm_1_199144473_G + rs12758027_C + imm_2_102013732_A + ccc.2.102135805.G.T_G + imm_2_181660357_A + imm_2_181756697_C + X1kg_2_241212412_A + imm_5_35917200_A + imm_5_40560617_A + imm_5_96137839_G + imm_5_96278345_C + imm_5_158759370_A + rs17765610_G + imm_6_91047852_A + rs2402752_C + imm_10_80730323_A + imm_12_6317038_A + X1kg_12_6361386_A + imm_12_110346958_T + imm_16_28236248_A + imm_17_23120724_A + imm_17_23172294_G + rs9901869_G + imm_19_10386181_A + imm_19_10486067_A + imm_21_39388614_A + imm_21_44440169_G + imm_22_20286653_G


## ### Frome Gethin from NYAS + genetic 
## formula <- PHEN ~ HLAB27 + FCGBP + NDUFB4 + ID3 + BZRAP1 + E2F2 + ANTXR2

## + imm_1_67478546_A + imm_1_67520024_A + imm_1_159121069_T + imm_1_199144473_G + imm_2_102013732_A + ccc.2.102135805.G.T_G + imm_2_181660357_A + imm_2_181756697_C + X1kg_2_241212412_A + imm_5_40560617_A + imm_5_96137839_G + imm_5_96278345_C + imm_5_158759370_A + rs17765610_G + imm_6_91047852_A + rs2402752_C + imm_10_80730323_A + imm_12_6317038_A + X1kg_12_6361386_A + imm_12_110346958_T + imm_16_28236248_A + imm_17_23120724_A + rs9901869_G + imm_19_10386181_A + imm_21_39388614_A + imm_21_44440169_G + imm_22_20286653_G



## formula <- PHEN ~  HLAB27 + BZRAP1 +  EDG8 + GPR128 + CD83 + Hs18s # originally and center

## formula <- PHEN ~  HLAB27 + BZRAP1 +  COMT + FCGBP + FZD2  + ERAP1  + INSIG1 +  MOAP1 +  ERAP2 + TNFAIP3  # new genes


## formula <- PHEN ~  BZRAP1 +  EDG8 + GPR128 + CD83 + Hs18s # originally and center

## formula <- PHEN ~  BZRAP1 +  COMT + FCGBP + FZD2   + INSIG1 +  MOAP1 + TNFAIP3  # new genes

## formula <- PHEN ~ BZRAP1 +  COMT + FCGBP + FZD2  + ERAP1  + INSIG1 +  MOAP1 +  ERAP2 + TNFAIP3 + EDG8 + GPR128 + CD83 + Hs18s # All genes

## #formula <- PHEN ~  HLAB27 +  BZRAP1 + COMT + FCGBP + ERAP1 +  ERAP2+ FZD2 +  EDG8 + GPR128 + INSIG1 +  MOAP1

## #formula <- PHEN ~  BZRAP1 + COMT + FCGBP + ERAP1 +  ERAP2 + FZD2 +  EDG8 + GPR128 + INSIG1 +  MOAP1 # USE THIS IF HLAB27 a random variable

## formula <- PHEN ~ HLAB27 + BZRAP1 +  COMT + FCGBP + FZD2  + ERAP1  + INSIG1 +  MOAP1 +  ERAP2 + TNFAIP3 +
##   imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A +
##   imm_2_62404976_A + imm_2_102030060_A + rs11098964_A +
##   imm_5_40526547_A + ccc.5.96150086.T.C_A + imm_5_158751323_A +
##   imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A +
##   imm_21_39391390_G

## formula <- PHEN ~ HLAB27 + BZRAP1 +  COMT + FCGBP + FZD2  + ERAP1  + INSIG1 +  MOAP1 +  ERAP2 + TNFAIP3 + EDG8 +
##   GPR128 + CD83 + Hs18s +
##   imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A +
##   imm_2_62404976_A + imm_2_102030060_A + rs11098964_A +
##   imm_5_40526547_A + ccc.5.96150086.T.C_A + imm_5_158751323_A +
##   imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A +
##   imm_21_39391390_G

#formula <- PHEN ~ HLAB27 + imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A + imm_2_62404976_A + imm_2_102030060_A + rs11098964_A +
  ## imm_5_40526547_A + ccc.5.96150086.T.C_A + imm_5_158751323_A +
  ## imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A +
  ## imm_21_39391390_G

#######genetic only
formula <- PHEN ~ HLAB27 + 
imm_1_67478546_A + imm_1_67520024_A + imm_1_159121069_T + imm_1_199144473_G + rs12758027_C + imm_2_102013732_A + ccc.2.102135805.G.T_G + imm_2_181660357_A + imm_2_181756697_C + X1kg_2_241212412_A + imm_5_35917200_A + imm_5_40560617_A + imm_5_96137839_G + imm_5_96278345_C + imm_5_158759370_A + rs17765610_G + imm_6_91047852_A + rs2402752_C + imm_10_80730323_A + imm_12_6317038_A + X1kg_12_6361386_A + imm_12_110346958_T + imm_16_28236248_A + imm_17_23120724_A + imm_17_23172294_G + rs9901869_G + imm_19_10386181_A + imm_19_10486067_A + imm_21_39388614_A + imm_21_44440169_G + imm_22_20286653_G


### 
formula <- PHEN ~ imm_1_67478546_A + imm_1_67520024_A + imm_1_159121069_T + imm_1_199144473_G + imm_2_102013732_A + ccc.2.102135805.G.T_G + imm_2_181660357_A + imm_2_181756697_C + X1kg_2_241212412_A + imm_5_40560617_A + imm_5_96137839_G + imm_5_96278345_C + imm_5_158759370_A + rs17765610_G + imm_6_91047852_A + rs2402752_C + imm_10_80730323_A + imm_12_6317038_A + X1kg_12_6361386_A + imm_12_110346958_T  + imm_17_23120724_A + rs9901869_G + imm_19_10386181_A + imm_21_39388614_A + imm_21_44440169_G + imm_22_20286653_G + imm_16_28236248_A

formula <- PHEN ~ HLAB27 + BZRAP1 +  COMT + FCGBP + FZD2  + ERAP1  + INSIG1 +  MOAP1 +  ERAP2 + TNFAIP3
+
formula <- PHEN ~ HLAB27 + imm_1_67478546_A + imm_1_67520024_A + imm_1_159121069_T + imm_1_199144473_G + imm_2_102013732_A + ccc.2.102135805.G.T_G + imm_2_181660357_A + imm_2_181756697_C + X1kg_2_241212412_A + imm_5_40560617_A + imm_5_96137839_G + imm_5_96278345_C + imm_5_158759370_A + rs17765610_G + imm_6_91047852_A + rs2402752_C + imm_10_80730323_A + imm_12_6317038_A + X1kg_12_6361386_A + imm_12_110346958_T + imm_16_28236248_A + imm_17_23120724_A + rs9901869_G + imm_19_10386181_A + imm_21_39388614_A + imm_21_44440169_G + imm_22_20286653_G
#################################### RUN 

formula <- PHEN ~ HLAB27 +imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A + imm_2_62404976_A + imm_2_102030060_A + rs11098964_A + imm_5_40526547_A + ccc.5.96150086.T.C_A + imm_5_158751323_A + imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A + imm_21_39391390_G




############################################ RUN the formulae
#PHEN ~ something + (1|ethnicity) + (1|site)  
predictions0 <- list()   
labels <- list()
problem.boot<-{}
top.performers<-{}

for(i in 1:Nsplits) {
  train <- data.list[[i]]$train
  test <- data.list[[i]]$test
  ## table(train$PHEN)
  ## table(test$PHEN)
  ## table(train$imm_16_28236248_A)

  fit1 <- glm(formula,data=train, family=binomial(logit),control=glm.control(epsilon=1e-16,maxit=1000))
 #fit1 <- lm(formula,data=train) # gives the same results as above
 #  fit1 <- glmmPQL(formula , data=train, random= (~1 |ethnicity) , family=binomial(logit),control=glm.control(epsilon=1e-16,maxit=1000)) 
#  fit1 <- glmmPQL(formula , data=train, random= list((~1 | ethnicity),(~1 | HLAB27)) , family=binomial(logit),control=glm.control(epsilon=1e-16,maxit=1000))
 ##  print("---------------------------------------")
 ## print(i)
 ## print(summary(fit1))
 ##   print("---------------------------------------")
  
  pred <- predict(fit1,newdata=test,type='response')
  predictions0[[i]] <- pred

  
  top.performers<-c(top.performers,names(sort(summary(fit1)$coefficients[,"Pr(>|z|)"],decreasing=FALSE))[1:5])
 if(sum(abs(summary(fit1)$coefficients[,"z value"])< 1e-6)){if(length(problem.boot)==0){problem.boot<-i}else{problem.boot<-c(problem.boot,i)}}

  labels[[i]] <-  test$PHEN

  
}

sort(tapply(top.performers,top.performers,length),decreasing=TRUE)
names(sort(tapply(top.performers,top.performers,length),decreasing=TRUE))
gsub(", "," + ",toString(names(sort(tapply(top.performers,top.performers,length),decreasing=TRUE))))

##  sort(tapply(top.performers,top.performers,length),decreasing=TRUE)
##        EDG8       MOAP1        CD83     SH3BP5L     HLADRB5        COMT 
##         188         152         108          80          60          52 
##        FZD2   GABARAPL1     TNFAIP3        GZMK        SNX5      GPR128 
##          50          42          36          33          28          21 
##       SIRPA       PTRH2      TGFBR3       DDIT4        BATF       ERAP1 
##          17          14          14          11          10          10 
##      RNASE6       ERAP2       FCGBP     TSC22D3       IL23R (Intercept) 
##          10           9           9           9           7           7 
##       RFWD2       DUSP2      INSIG1       PSMG1       CXCR4 
##           7           6           5           4           1

## "EDG8 + MOAP1 + CD83 + SH3BP5L + HLADRB5 + COMT + FZD2 + GABARAPL1 + TNFAIP3 + GZMK + SNX5 + GPR128 + SIRPA + PTRH2 + TGFBR3 + DDIT4 + BATF + ERAP1 + RNASE6 + ERAP2 + FCGBP + TSC22D3 + IL23R + (Intercept) + RFWD2 + DUSP2 + INSIG1 + PSMG1 + CXCR4"
## model <- glm(formula,data=data,family=binomial(logit),control=glm.control(epsilon=1e-16,maxit=1000))
## error<-cv.glm(data,model)
## error$delta^2
## class(error)
length((predictions0[[1]]))
length((predictions0[[1]]))
## Get plots
length(predictions0)
length(problem.boot)
dim(data.list[[1]]$train)

if(length(problem.boot)>0){ ### these indicated instances where one covaraiate not polymorphic/ variable for a givem predict/train combination they do not affect the requst
predictions0<-predictions0[ -1*problem.boot ]
labels<-labels[ -1*problem.boot ]
}
predict.gene.definite<-predictions0
length(predictions0)

model<-"HLA_only"
model<-"Axial_SpA_33_genes"
model<-"Axial_SpA_33_CEU_genes"
model<-"Axial_SpA_33_B27_genes"

model<-"Axial_SpA_33_genes_HLA"
model<-"Axial_SpA_3_CEU_genes_HLA"
model<-"Axial_SpA_10_B27_genes_HLA"

 model<-"Diagnosis_2_genes_figure_2A"
 model<-"Diagnosis_2_genes_HLA_figure_3A"
# model<-"Diagnosis_3_genes_HLA_figure_3A"
# model<-"Diagnosis_6_genes_figure_2B"
 model<-"Diagnosis_5_genes_figure_2B"
# model<-"Diagnosis_9_genes_figure_2C"  #  model<-"Diagnosis_9_genes_figure_2C_NO_BAD_SAMPLES"
# model<-"Diagnosis_9_genes_CEU_only_figure_2D"
# model<-"Diagnosis_9_genes_HLA_CEU_only_figure_2D"
# model<-"Diagnosis_6_genes_HLA_figure_3B" # model<-"Diagnosis_6_genes_HLA_figure_3B_NO_BAD_SAMPLES"
model<-"Diagnosis_5_genes_HLA_figure_3B"
# model<-"Diagnosis_9_genes_HLA_figure_3C"

# model<-"Diagnosis_core_genes_figure_4D"
# model<-"Diagnosis_core_genes_HLA_figure_4B"
# model<-"NYAS_core_genes_figure_4A"
# model<-"NYAS_core_genes_HLA_figure_4C"

# model<-"Diagnosis_model_choice_figure_AA"
# model<-"Diagnosis_model_choice_HLA_figure_AB"
 model<-"Axial_SpA_9_genes_figure"
 model<-"Axial_SpA_9_genes_HLA_figure"
 model<-"Axial_SpA_HLA_only_figure"

 model<-"immuno_chip_genetc_only"
model<-"immuno_chip_genetc_only_noHLA"

model<-"immuno_chip_genetc_only_HLAneg_cases_vs_controls"

model<-"Axial_Spa_positive_immuno_chip_Controls_genetc_only"
model<-"Axial_Spa_positive_HLAneg_vs_immuno_chip_Controls_genetc_only"


model<-"Axial_Spa_positive_immuno_chip_Controls_genetc_only_noHLA"

model<-"Axial_Spa_positive_MRI_positive_immuno_chip_Controls_genetic"
model<-"Axial_Spa_positive_MRI_negative_immuno_chip_Controls_genetic"
model<-"IGNORE_B27_Axial_Spa_positive_immuno_chip_Controls_genetic"


model<-"Image_Positive_SPA_immuno_chip_Cases_genetic"
model<-"Image_Negative_SPA_immuno_chip_Cases_genetic"

model<-"immuno_chip_Cases_Image_Negative_SPA_genetic"

model<-"Image_Positive_SPA_immuno_chip_Controls_genetic"
model<-"Image_Negative_SPA_immuno_chip_Controls_genetic"

model<-"Not_Image_Positive_immuno_chip_Controls_genetic"

model<-"Axial_Spa_negative_immuno_chip_Controls_genetic_only"

 model<-"Diagnosis_9_genes_HLA_genetic_figure_3C"  #  model<-"Diagnosis_9_genes_HLA_genetic_figure_3C_NO_BAD_SAMPLES"
model<-"Diagnosis_9_genes_HLA_genetic_reduced_figure_5B2"
# model<-"Diagnosis_6_genes_HLA_genetic_reduced_figure_5B2"
# model<-"NYAS_core_genes_HLA_genetic_reduced_figure_5A2"
 model<-"test"

model<-"Genetics_axial_spa_wtcc_controls_genetic"

model<-"Axial_spa_cases_v_controls_genetic"

model<-"IGAS_AS_CASES_V_ ASAS_IMAGING_NEGATIVE_CONTROLS_genetic"

model<-"ASAS_IMAGING_POSITIVE_CASES_v_ASAS_IMAGING_NEGATIVE_CONTROLS_genetic"

model<-"IGAS_AS_CASES_v_TOTAL_ASAS_CONTROLS"

model<-"ASAS_IMAGING_POSITIVE_CASES_v_TOTAL_ASAS_CONTROLS"

model<-"IGAS_CASES_v_ASAS_IMAGING_POSITIVE_CASES"

model<-"IGAS_CASES_v_ASAS_IMAGING_POSITIVE_Controls_final"



model<-"IGAS_CASES_v_ASAS_IMAGING_POSITIVE_CASES_final"

model<-"ASAS_IMAGING_POSITIVE_CASES_v_ASAS_Controls_final"

model<-"ASAS_IMAGING_POSITIVE_CASES_v_IGAS_Controls_final"

model<-"ASAS_IMAGING_NEG_CASES_v_IGAS_Controls_final"




model<-"IGAS_CASES_v_ASAS_IMAGING_POSITIVE_CASES_genetic_final"

model<-"IGAS_CASES_v_ASAS_IMAGING_NEGATIVE_Controls_final_genetic"

model<-"ASAS_IMAGING_POSITIVE_CASES_v_IGAS_Controls_genetic_final"

model<-"ASAS_IMAGING_NEGATIVE_CONTROLS_v_IGAS_Controls_genetic_final"

formula
model


pred0 <- prediction(predictions0,labels)
perf0 <- performance(pred0, measure = "tpr", x.measure = "fpr")

length(pred0)
plot(perf0,avg="threshold",lwd=2,col=1,colorize=F)
plot(perf0,lwd=2,col=1,colorize=F,add=T)

## # it allows two different plots in the same frame
## par(mfrow = c(1,2))
## # plot a ROC curve for a single prediction run
## # and color the curve according to cutoff.
## library(ROCR)
## data(ROCR.simple)
## pred <- prediction(ROCR.simple$predictions, ROCR.simple$labels)
## perf <- performance(pred,"tpr", "fpr")
## plot(perf,colorize = TRUE)
## # plot a ROC curve for a single prediction run
## # with CI by bootstrapping and fitted curve
## library(verification)
## roc.plot(ROCR.simple$labels,ROCR.simple$predictions, xlab = "False positive rate",
## ylab = "True positive rate", main = NULL, CI = T, n.boot = 100, plot = "both", binormal = TRUE)
## (auc <- as.numeric(performance(pred, measure = "auc", x.measure = "cutoff")@y.values))


par(mar=c(4.5,4.5,2.5,2),mgp=c(3,1,0)) #c(bottom, left, top, right)
perf0.acc <- performance(pred0, "acc")
plot(perf0.acc, avg= "vertical", spread.estimate="boxplot", lwd=3,col='blue', show.spread.at= seq(0.1, 0.9, by=0.1),main= "",cex=1.5,cex.lab=1.75,font=2,xaxis.cex.axis=2,yaxis.cex.axis=2,xlab.cex.axis=2,colorkey.relwidth=0.5)
                                        # main= "Accuracy across the range of possible cutoffs"
savePlot(paste(model,"Accuracy_vs_Cutoff.png",sep=""),type="png")

plot(perf0, avg= "threshold", colorize=T, lwd= 5, main= "",colorkey.relwidth=0.5,colorkey.pos="right",cex.lab=1.75,font=2,xaxis.cex.axis=2,yaxis.cex.axis=2,xlab.cex.axis=2)
#print.cutoffs.at=c(0.4,0.5,0.6)
auc<-performance(pred0,'auc')
auc<- unlist(auc@y.values)
the.auc<-mean(auc,trim=0.1)
the.auc
sd.auc<-sd(auc)
sd.auc


x<-0.5 # x<-0.8

sen<-performance(pred0,'sens')
datax<-sen@x.values
datay<-sen@y.values

##
temp<-lapply(datax,length)
tapply(unlist(temp),unlist(temp),length)
test<-tapply(unlist(temp),unlist(temp),length)
to.keep<-max(test)
to.keep<-as.numeric(names(test)[test==max(test)])
to.keep
if(length(to.keep)>1){to.keep<-to.keep[1]} # in case multiple ones of the same length
datax<-datax[unlist(temp) ==to.keep]
datay<-datay[unlist(temp) ==to.keep]
length(datax)




av.x<-apply(do.call("rbind",datax),2,function(x) mean(x,trim=0.5))
av.y<-apply(do.call("rbind",datay),2,function(x) mean(x,trim=0.5)) # need trim to remove outlyers
sd.y<-apply(do.call("rbind",datay),2,function(x) sd(x))

trim<-1
the.length<-length(av.x)
av.sens<-aspline(as.numeric(av.x)[(1+trim):(the.length-trim)], as.numeric(av.y)[(1+trim):(the.length-trim)],x,method="improved",degree=3)$y
sd.sens<-aspline(as.numeric(av.x)[(1+trim):(the.length-trim)], as.numeric(sd.y)[(1+trim):(the.length-trim)],x,method="improved",degree=3)$y

## plot(av.x,av.y)
## plot(av.x,sd.y)
## length(as.numeric(av.x)[(1+trim):(the.length-trim)])
## length(as.numeric(av.y)[(1+trim):(the.length-trim)])


spec<-performance(pred0,'spec')
datax<-spec@x.values
datay<-spec@y.values

temp<-lapply(datax,length)
tapply(unlist(temp),unlist(temp),length)
test<-tapply(unlist(temp),unlist(temp),length)
to.keep<-max(test)
to.keep<-as.numeric(names(test)[test==max(test)])
if(length(to.keep)>1){to.keep<-to.keep[1]}
datax<-datax[unlist(temp) ==to.keep]
datay<-datay[unlist(temp) ==to.keep]
length(datax)


av.x<-apply(do.call("rbind",datax),2,function(x) mean(x,trim=0.5))
av.y<-apply(do.call("rbind",datay),2,function(x) mean(x,trim=0.5)) # need trim to remove outlyers
sd.y<-apply(do.call("rbind",datay),2,function(x) sd(x))
trim<-1
the.length<-length(av.x)
av.spec<-aspline(as.numeric(av.x)[(1+trim):(the.length-trim)], as.numeric(av.y)[(1+trim):(the.length-trim)],x,method="improved",degree=3)$y
sd.spec<-aspline(as.numeric(av.x)[(1+trim):(the.length-trim)], as.numeric(sd.y)[(1+trim):(the.length-trim)],x,method="improved",degree=3)$y




## plot(av.x,av.y,col="blue",cex=3)
## auc<- unlist(auc@y.values)
## the.auc<-mean(auc,trim=0.1)
## the.auc
## sd.auc<-sd(auc)
## sd.auc

leg.txt<-c(paste("Av. AUC        : ",signif(the.auc,digits=2), expression("\u00B1") ,signif(sd.auc,digits=1),sep=" "), # use character map to define
           paste("Av. Sensitivity: ",signif(av.sens,digits=2), expression("\u00B1") ,signif(abs(sd.sens),digits=1),sep=" "),
              paste("Av. Specificity: ",signif(av.spec,digits=2), expression("\u00B1") ,signif(abs(sd.spec),digits=1),sep=" ") )

leg.txt
# leg.txt<-gsub("0.7 ","0.70 ",leg.txt)
# leg.txt<-gsub("0.6 ","0.60 ",leg.txt)
#text(0.1,0.2,expression(paste("Av. AUC: ",eval(val),sep="") %+-% 1),cex=2)
#legend(-0.1,1.0,legend=leg.txt,cex=1.5,bty="n")

legend(0.3,0.6,legend=leg.txt,cex=1.5,bty="n")
#legend(0.2,0.8,legend=leg.txt,cex=1.5,bty="n")
#text(0,0.8,labels=leg.txt,cex=1.5,pos=4)
paste(model,"_ROC.png",sep="")
getwd()
savePlot(paste(model,"_ROC.png",sep=""),type="png")
savePlot(paste(model,"_ROC.tiff",sep=""),type="tiff")




###################################### SAVE a hi resolution epx figure###############
X11(width = 7, height = 7) 
plot(perf0, avg= "threshold", colorize=T, lwd= 5, main= "",colorkey.relwidth=0.5,colorkey.pos="right",cex.lab=1.75,font=2,xaxis.cex.axis=2,yaxis.cex.axis=2,xlab.cex.axis=2)

#legend(-0.1,1.0,legend=leg.txt,cex=1.5,bty="n")

legend(0.3,0.6,legend=leg.txt,cex=1.5,bty="n")

dev.copy2eps(file=paste(model,"_ROC.eps",sep=""))
dev.off()


##################### straight ROC no color
plot(perf0, avg= "threshold", colorize=F, lwd= 5, main= "",cex.lab=1.75,font=2,xaxis.cex.axis=2,yaxis.cex.axis=2,xlab.cex.axis=2)


#legend(-0.1,1.0,legend=leg.txt,cex=1.5,bty="n")
legend(0.3,0.6,legend=leg.txt,cex=1.5,bty="n")

savePlot(paste(model,"_ROC_noColour.png",sep=""),type="png")
savePlot(paste(model,"_ROC_noColour.tiff",sep=""),type="tiff")
dev.off()
########## EPS

X11(width = 7, height = 7) 
plot(perf0, avg= "threshold", colorize=F, lwd= 5, main= "",cex.lab=1.75,font=2,xaxis.cex.axis=2,yaxis.cex.axis=2,xlab.cex.axis=2)

#legend(-0.1,1.0,legend=leg.txt,cex=1.5,bty="n")
legend(0.3,0.6,legend=leg.txt,cex=1.5,bty="n")

dev.copy2eps(file=paste(model,"_ROC_noColour.eps",sep=""))
dev.off()

##########################################################################




auc<-performance(pred0,'auc')
auc<- unlist(auc@y.values)
the.auc2<-mean(auc,trim=0.1)
the.auc2
sd.auc2<-sd(auc)
sd.auc2




leg.txt<-c(paste("Av. AUC        : ",signif(the.auc,digits=2), expression("\u00B1") ,signif(sd.auc,digits=2),sep=" "), # use character map to define
           paste("Av. Sensitivity: ",signif(av.sens,digits=2), expression("\u00B1") ,signif(sd.sens,digits=2),sep=" "),
              paste("Av. Specificity: ",signif(av.spec,digits=2), expression("\u00B1") ,signif(sd.spec,digits=2),sep=" "),
              paste("HLAB27 only Av. AUC: ",signif(the.auc2,digits=2), expression("\u00B1") ,signif(sd.auc2,digits=2), "(-)",sep=" "))

legend(0.05,0.7,legend=leg.txt,cex=1.5,bty="n")



###########################################
datax<-perf0@x.values
datay<-perf0@y.values


#### excluding sex mimatched 
## sum<-{}
## for(i in 1:length(datax)){
##   if(i==1){
## sumx<-datax[[i]]
## sumy<-datay[[i]]
## }else{
##   sumx<-sumx+datax[[i]]
##   sumy<-sumy+datay[[i]]
## }
## }
## av.x<-sumx/length(datax)
## av.y<-sumy/length(datay)
## points(av.x,av.y,col="red",cex=3)
## av.x<-apply(do.call("rbind",datax),2,function(x) mean(x,trim=0.5))
## av.y<-apply(do.call("rbind",datay),2,function(x) mean(x,trim=0.5)) # need trim to remove outlyers
## points(av.x,av.y,col="blue",cex=3)
## ###########################
## library(verification)
## library(akima)
## table.predicts<-do.call("rbind",predictions0)
## table.labels<-do.call("rbind",labels)

## roc.plot(t(table.labels)[,1],t(table.predicts)[,5], xlab = "False positive rate",ylab = "True positive rate", main = NULL, CI = T, n.boot = 100, plot = "both", binormal = TRUE)

## roc.plot(t(table.labels)[,1],t(table.predicts), xlab = "False positive rate",ylab = "True positive rate", main = NULL, CI = T, n.boot = 100, plot = "both", binormal = TRUE)

#####################















x <- seq(-4, 4, len = 101)
y <- cbind(sin(x), cos(x))
matplot(x, y, type = "l", xaxt = "n",
        main = expression(paste(plain(sin) * phi, "  and  ",
                                plain(cos) * phi)),
        ylab = expression("sin" * phi, "cos" * phi), # only 1st is taken
        xlab = expression(paste("Phase Angle ", phi)),
        col.main = "blue")
axis(1, at = c(-pi, -pi/2, 0, pi/2, pi),
     labels = expression(-pi, -pi/2, 0, pi/2, pi))









savePlot(paste(model,"_cutoff.png",sep=""),type="png")

auc0 <- performance(pred0,'auc')
auc0 <- unlist(auc0@y.values)


mean(auc0)
sd(auc0)

auc0s <- performance(pred0,'sens')
auc0s <- unlist(auc0@y.values)


mean(auc0s)
sd(auc0s)

########################### immumo chip trained:
##########################################################################
#########################################################################
#########################################################################

############# ALLPY A KNOWN MODEL TO PREDICTE DATA:

#data.all<-data #
data<-data.all

############################

data[,"PHEN"]<-data[,"diagnosis"]
data[,"PHEN"]<-data[,"ASAS_axial_spa_criteria"]

sum(data[,"PHEN"] != data[,"ASAS_axial_spa_criteria"])
sum(data[,"PHEN"] != data[,"diagnosis"])


sum(data[,"Definite.AS"])
sum(data[,"Def.Not.AS"])
sum(data[,"Definite.AS"] & !data[,"Def.Not.AS"])

data[,"Definite.AS"]
data[,"Def.Not.AS"]

NYAS<-rep(NA,times=dim(data)[1])
NYAS[data[,"Definite.AS"]]<-1
NYAS[data[,"Def.Not.AS"]]<-0
sum(NYAS

sum(data[,"HLAB27"] != data[,"ASAS_axial_spa_criteria"])

wanted<-data[,"ethnicity"] %in% c(1,4,8)
data<-data[wanted,] # keep only core ethnicity



####################################### 2013 test tp 
data[,"PHEN"]<-data[,"ASAS_axial_spa_criteria"]
data <- data[!is.na(data$PHEN),]




library(verification)
load("imm.data.RData")
colnames(imm.data)
table(imm.data[,"PHEN"])
table(data[,"PHEN"])


    #######################
colnames(data)    
## colnames(imm.data)
## imm.aff<-imm.data[,"PHEN"]==1
## imm.b27<-imm.data[,"HLAB27"]==1

## sum(imm.aff)
## sum(imm.b27)
## sum(imm.aff & imm.b27)

## sum(imm.aff & imm.b27)/sum(imm.aff)

## sum(imm.data[imm.data[,"PHEN"]==1,"HLAB27"]==1


formula <- PHEN ~  BZRAP1 +  COMT + FCGBP + FZD2  + INSIG1 +  MOAP1  + TNFAIP3 + EDG8 +
  GPR128 + CD83 + Hs18s

formula <- PHEN ~  HLAB27+  BZRAP1 +  COMT + FCGBP + FZD2  + ERAP1  + INSIG1 +  MOAP1 +  ERAP2 + TNFAIP3 + EDG8 +
  GPR128 + CD83 + Hs18s


formula <- PHEN ~ HLAB27 + BZRAP1 +  COMT + FCGBP + FZD2  + ERAP1  + INSIG1 +  MOAP1 +  ERAP2 + TNFAIP3 + EDG8 +
  GPR128 + CD83 + Hs18s +
  imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A +
  imm_2_62404976_A + imm_2_102030060_A + rs11098964_A +
  imm_5_40526547_A + ccc.5.96150086.T.C_A + imm_5_158751323_A +
  imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A +
  imm_21_39391390_G

formula <- PHEN ~  BZRAP1 +  COMT + FCGBP + FZD2  + ERAP1  + INSIG1 +  MOAP1 +  ERAP2 + TNFAIP3 + EDG8 +
  GPR128 + CD83 + Hs18s +
  imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A +
  imm_2_62404976_A + imm_2_102030060_A + rs11098964_A +
  imm_5_40526547_A + ccc.5.96150086.T.C_A + imm_5_158751323_A +
  imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A +
  imm_21_39391390_G

formula <- PHEN ~ HLAB27 + imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A +
  imm_2_62404976_A + imm_2_102030060_A + rs11098964_A +
  imm_5_40526547_A + ccc.5.96150086.T.C_A + imm_5_158751323_A +
  imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A +
  imm_21_39391390_G

formula <- PHEN ~  imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A +
  imm_2_62404976_A + imm_2_102030060_A + rs11098964_A +
  imm_5_40526547_A + ccc.5.96150086.T.C_A + imm_5_158751323_A +
  imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A +
  imm_21_39391390_G

## formula.gen <- PHEN ~ imm_1_25169065_A + imm_1_67478546_A + imm_1_199226930_A +
##   imm_2_62404976_A + imm_2_102030060_A + rs11098964_A +
##   imm_5_40526547_A + ccc.5.96150086.T.C_A + imm_5_158751323_A +
##   imm_9_138389159_A + ccc.12.6373003.G.A_A + rs12943464_A +
##   imm_21_39391390_G

model <- glm(formula,data=data,family=binomial(logit),control=glm.control(epsilon=1e-16,maxit=1000))

model <- glmmPQL(formula , data=data, random= (~1 | ethnicity) , family=binomial(logit),control=glm.control(epsilon=1e-16,maxit=1000))

model <- glmmPQL(formula.gen , data=imm.data, random= (~1 | HLAB27) , family=binomial(logit),control=glm.control(epsilon=1e-16,maxit=1000))

model <- glm(formula,data=imm.data,family=binomial(logit),control=glm.control(epsilon=1e-16,maxit=1000))
model <- glmmPQL(formula , data=imm.data, random= (~1 | Centre) , family=binomial(logit),control=glm.control(epsilon=1e-16,maxit=1000))

###############$$$$$$### immunochip train on half 
## samples.labels<-permute(1:dim(imm.data)[1])
## totalN<-length(samples.labels)
## half<-totalN/2
## table(imm.data[samples.labels[1:half],"PHEN"])

## model <- glm(formula,data=imm.data[samples.labels[1:half],],family=binomial(logit),control=glm.control(epsilon=1e-16,maxit=1000))
## pred.probs <- predict(model,newdata=imm.data[samples.labels[half+1:totalN],],type='response')
## the.pred<-prediction(as.numeric(pred.probs),imm.data[samples.labels[half+1:totalN],"PHEN"])
##############$$$$$

colnames(data)
colnames(imm.data) %in% colnames(data)

## trained model with AS definate or spa_criate and apply to daignosis:
data<-data.all
data[,"PHEN"]<-data[,"diagnosis"]
data <- data[!is.na(data$PHEN),]

wanted<-data[,"ethnicity"] %in% c(8)
wanted.c<-data[,"Center"] %in% c(1)
sum(wanted & wanted.c)
data<-data[wanted & wanted.c ,]
table(data[,"PHEN"])

100*sum(data[,"PHEN"]==1 & data[,"HLAB27"]==1)/sum(data[,"PHEN"]==1)
###################################


pred.probs <- predict(model,newdata=data,type='response')
missing<-is.na(pred.probs)
the.pred<-prediction(as.numeric(pred.probs)[!missing],data[!missing,"PHEN"])


the.pred.perf <- performance(the.pred, measure = "tpr", x.measure = "fpr")
plot(the.pred.perf,avg="threshold",lwd=2,col="green",colorize=F)
plot(the.pred.perf,avg="threshold",lwd=2,col="green",colorize=F,add=T) # this same as below


perf.acc <- performance(the.pred, "acc")


plot(perf.acc, avg= "vertical", spread.estimate="boxplot", lwd=3,col='blue',
     show.spread.at= seq(0.1, 0.9, by=0.1),
     main= "Accuracy across the range of possible cutoffs")


par(mar=c(4.5,6.5,2.5,2),mgp=c(2.75,1,0)) #c(bottom, left, top, right)
roc.plot(data[,"PHEN"],pred.probs, xlab = "False positive rate",ylab = "True positive rate", main = NULL, CI = T, n.boot = 100, plot = "both", binormal = TRUE, plot.thres = seq(0.1,0.9, 0.1),cex.lab=1.75,font=1)


## a.test<-roc.plot(imm.data[samples.labels[half+1:totalN],"PHEN"],pred.probs, xlab = "False positive rate",ylab = "True positive rate", main = NULL, CI = T, n.boot = 100, plot = "both", binormal = TRUE, plot.thres = seq(0.1,0.9, 0.1),cex.lab=1.75,font=1,cex.axis=2,yaxis.cex.axis=2,xlab.cex.axis=2,col="green",add=T)


auc0 <- performance(the.pred, "auc")
the.auc<- unlist(auc0@y.values)  # area under curve

leg.txt<-c(paste("AUC : ",signif(the.auc,digits=3),sep=""))

legend(0.5,0.5,legend=leg.txt,cex=1.5,bty="n",text.col="red")
the.auc
formula

legend(0.5,0.4,legend=leg.txt,cex=1.5,bty="n",text.col="green")

savePlot("ImmunoChip_trained_to immunochip_NOHLA.png",type="png")
savePlot("ImmunoChip_trained_to immunochip.png",type="png")

savePlot("ImmunoChip_trained_to diagnosis_NOHLA.png",type="png")
savePlot("ImmunoChip_trained_to diagnosis.png",type="png")

savePlot("Spa_criteria_trained_gene_snps_NOHLA_to diagnosis.png",type="png")

savePlot("ASDefinite_trained_gene_to diagnosis.png",type="png")

savePlot("Spa_criteria_trained_gene_snps_to diagnosis.png",type="png")





performance(the.pred,'auc')@y.values # the area under the curve

performance(the.pred,measure='sens',x.measure="spec")

sum(is.na(pred.probs))
(pred.probs)

cutoff <- .2
pred.labs <- rep(0,length(pred.probs))
pred.labs[pred.probs > cutoff] <- 1
pred.labs



plot(perf0, avg= "threshold", colorize=T, lwd= 5, main= "",colorkey.relwidth=0.5,colorkey.pos="right",cex.lab=1.75,font=2,xaxis.cex.axis=2,yaxis.cex.axis=2,xlab.cex.axis=2)
