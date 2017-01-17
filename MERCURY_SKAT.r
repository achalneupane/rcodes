

###################### load the data file provided
library(skatMeta)  ### you need this R-library
load("MERCURY_BONE.RData") ## sent via cloudstor


######################  CARRIE YOU NEED TO ADD THE *RESIDUAL BMD* TO COLUMN **SampleProject** in pheno: (it's all NA at the moment)
pheno[1:5,]
##    SAMPLE PATIENT SampleProject FamilyCode ParticipantCode PaternalID MaternalID Sex AffectionStatus  all
## 1 MrOS001 MrOS001            NA        ALL         MrOS001          0          0   0               1 TRUE
## 2 MrOS002 MrOS002            NA        ALL         MrOS002          0          0   0               1 TRUE
## 3 MrOS003 MrOS003            NA        ALL         MrOS003          0          0   0               1 TRUE
## 4 MrOS004 MrOS004            NA        ALL         MrOS004          0          0   0               1 TRUE
## 5 MrOS005 MrOS005            NA        ALL         MrOS005          0          0   0               1 TRUE

covars<-"1"
target.pheno.col<-"SampleProject"
formula<-paste(target.pheno.col,"~",paste(covars,collapse="+"),sep="")
print(formula)
formula<-formula(formula)
#### ADD the residual BMD to the SampleProject column:
formula # this is the formule # skat will use
a.indel[1:5,1:50]
dim(a.indel) # 440130   1379 all the annotated genotypes

############################ this is for the coding mtation with MAF< 0.01
pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype & !are.repeats &
ok.missing & hw.controls.ok.filt & !no.genotypes &  rare.in.all #  & !on.x.y

sum(pass) # 56385

help<-cbind(full.qual,bad.coding,maf.filter,rare.in.group,no.genotypes,in.common.hit.gene ,hw.controls.ok,on.x.y,unannotated.hits,not.flat.genotype,are.repeats,ok.missing,ok.missing,is.unwound.geno,is.unwound.geno ,hw.p.control.filt, rare.in.all) ## summary of QC




################################# GEFOS FILTERING cause sending all

genotypes<-a.indel[pass,the.samples.use] ## ordered correctly for phenotypes
snp.names<-key[pass] ## GEFOS ony name with start

#### snpinfo now A different size than a.indel since added pathways!!!
snpinfo<-snpinfo.ori[snpinfo.ori[,"Name"] %in% snp.names,]
if( sum(!(snp.names %in% snpinfo.ori[,"Name"]))>0){print("WARINING snp.names not in snpinfo- unusual!")}
dim(snpinfo)
length(snp.names)
dim(genotypes)


dim(genotypes)

print("start QC")
genotypes[genotypes=="NA"]<-NA
genotypes[genotypes=="0/0"]<-0
genotypes[genotypes=="0/1"]<-1
genotypes[genotypes=="1/1"]<-2

########### prevent any averaging
dim(genotypes)
genotypes[is.na(genotypes)]<-0
dim(genotypes)
########### prevent any averaging


#tapply(as.vector(genotypes),as.vector(genotypes),length)
num.col<-dim(genotypes)[2]
num.row<-dim(genotypes)[1]
## genotypes[1:5,1:20]
genotypes<-as.numeric(as.matrix(genotypes))
dim(genotypes)<-c(num.row,num.col)
genotypes<-t(genotypes) # samples x SNPS

colnames(genotypes)<-snp.names
rownames(genotypes)<-gsub(".GT$","",the.samples.use)


#a.indel[pass,][1:5,1:10]
###  if(target.pheno.col %in% case.control){
###  cohort.seq <- skatCohort(Z=genotypes,formula, SNPInfo = snpinfo, data=pheno,aggregateBy="cluster",family=binomial(),verbose=FALSE)

### }else{
cohort.seq <- skatCohort(Z=genotypes,formula, SNPInfo = snpinfo, data=pheno,aggregateBy="cluster",family=gaussian(),verbose=FALSE) ## genes and clusters
#}

meta.results.burden<-burdenMeta(cohort.seq,wts=1,mafRange = c(0,1),SNPInfo = snpinfo,aggregateBy="cluster")
meta.results.skat<-skatMeta(cohort.seq,SNPInfo = snpinfo,aggregateBy="cluster")
meta.results.skatO<-skatOMeta(cohort.seq,burden.wts =1,SNPInfo = snpinfo,aggregateBy="cluster")


#meta.results.skat<-{}
#meta.results.skatO<-{}

the.order<-     order(meta.results.burden[,"p"])
sum(is.na(meta.results.burden[,"p"])) ## bad p-values shoudl not happen
meta.results.burden<-  meta.results.burden[the.order,]

meta.results.burden[1:50,]
meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]

the.order<-     order(meta.results.skat[,"p"])
meta.results.skat<-  meta.results.skat[the.order,]
meta.results.skat[1:50,]

the.order<-     order(meta.results.skatO[,"p"])
sum(is.na(meta.results.skatO[,"p"])) ## bad p-values shoudl not happen
meta.results.skatO<-  meta.results.skatO[the.order,]
meta.results.skatO[1:50,]

clusters[1:5,]

##### write results-----------------------

write.table(meta.results.burden[1:50,],file=paste("Burden",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skat,file=paste("SkaO",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skatO,file=paste("SkatO",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
##### write results-----------------------

