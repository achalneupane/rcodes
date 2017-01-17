
library(GenomicFeatures)
library(HardyWeinberg)
library(Biostrings)

code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir) # load in the UCSC tables these use the db file names and not their lable-names 
source("annotate_SNPs_subroutines.r")


options("width"=250,"max.print"=1000)



vcf.files<-"/media/UQCCG/Sequencing/CompleteGenomics/Diamantina_T_79/Diamantina_T_79-Small_Variant_Table.tsv"
output.name<-"cg_filtered"

names(vcf.files)<-"snp"
snp.dir<-"/media/UQCCG/Sequencing/CompleteGenomics/Diamantina_T_79"


 # .txt extension will be added


location<-snp.dir
skip.lines<-86

########################BEGIN
load("/media/UQCCG/Sequencing/Data/Genomes/hg19/NexteraRapidCapture_Exome_TargetedRegions_hg19_targets.RData")
length(data.gr)
dim(data)
data.gr<-data.gr+200
data.gr<-reduce(data.gr)

########################

setwd(location)
the.files<-dir(getwd())

if(paste(output.name,"_maf.txt",sep="") %in% the.files){print("WARNING output files exits they will be appended too!!")}

################ LARGE FILES ########
  print(vcf.files)
######BELOW PROCESSING  this for snp for indel varient types in vcf4.0 or vcf 3.0 format
close(con)

column.labels<-read.delim(vcf.files,header=F,nrows=1,skip=skip.lines,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
num.vars<-length(column.labels)
num.vars

con <- file(vcf.files, open="r")  # close(con)
num.lines<-1 # so does while llop at least once
reads<-1000000  #1.3M lines in snp file 50000 goes out to 24Gb without QC cgeck 
counter<- -1
while (num.lines >0  ){
counter<-counter+1
print(counter)

if(counter==0){
indels<-try(scan(con,what=character(num.vars),skip=(reads*counter)+skip.lines+1,nlines=reads,sep="\t",fill=TRUE,na.strings="",quote="\""))
}else{
indels<-try(scan(con,what=character(num.vars),nlines=reads,sep="\t",fill=TRUE,na.strings="",quote="\""))
}
 ## indels1 <- read.table(con,col.names=column.labels,sep="\t",skip=skip.lines,fill=TRUE,stringsAsFactors=FALSE,colClasses="character",nrows=reads,comment.char="",quote="")
             

num.lines<-length(indels)/(num.vars)
print(num.lines)
if(num.lines==0){next}
dim(indels)<-c(num.vars,num.lines)
indels<-t(indels)
colnames(indels)<-column.labels

indels[1:5,]

indel.gr<-GRanges(seqnames =indels[,"chromosome"],ranges = IRanges(start=as.numeric(indels[,"begin"]),end=as.numeric(indels[,"end"])),strand="+")

data.gr

common<-overlapsAny(indel.gr,data.gr)


if(sum(common)>0){
   indels<-indels[common,]
   
 if(counter==0){
 write.table(indels,file="cg_table1_Rapid_exomes.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)  
      }else{
write.table(indels,file="cg_table1_Rapid_exomes.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,,append=TRUE)     
             }
 }

} # loop over chunks






cg<-read.delim("/media/UQCCG/Sequencing/CompleteGenomics/Diamantina_T_79/cg_table1_Rapid_exomes.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)



##############

code<-read.delim("/media/UQCCG/Sequencing/CompleteGenomics/ID_conversion_sample_sheet_check.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
code2<-read.delim("/media/UQCCG/Sequencing/CompleteGenomics/917_Data_Delivery_101614.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
code<-code[,1:3]
code[1:5,]
code2[1:5,]
extra<-code2[!(code2[,"Assembly.ID"] %in%  code[,"assembly_id"]),]
dim(extra)

extra<-extra[,c("Assembly.ID","Cust..Subject.ID","Cust..Sample.ID")]
colnames(extra)<-colnames(code)

code<-rbind(code,extra)

recodes<-read.delim("/media/UQCCG/Sequencing/CompleteGenomics/RecodingSampleID.csv",header=T,sep=",",fill=TRUE,stringsAsFactors=FALSE) 
recodes[1:5,]

dim(cg)
cg[1:5,1:10]
code[!(code[,"customer_sample_id"] %in% all.possible.samples),]

code[,"exomes.ids"]<-code[,"customer_sample_id"]
code[1:5,]
recodes[1:5,]
posns<-match(code[,"exomes.ids"],recodes[,"ID"])
missing<-is.na(posns)
cbind(code[!missing,"exomes.ids"],recodes[posns[!missing],])

code[!missing,"exomes.ids"]<-recodes[posns[!missing],"MhairiRecoded"]

code[!(code[,"exomes.ids"] %in% all.possible.samples),]

recodes

write.table(indels,file="cg_table1_Rapid_exomes.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,,append=TRUE) 

## /media/UQCCG/Sequencing/Projects/TGCM-AML/Analysis/CREST_Analysis_FLT3_exome capture DNAs updated Jan 2013 _TJG and AB.csv

## /media/UQCCG/Sequencing/CompleteGenomics/Final 51 AMLM12 Pts - with seq+FLt3+ Blast and Cyto 140929.csv

## /media/UQCCG/Sequencing/CompleteGenomics/20140918_QCTB-UQDI_AML_Cytogenetics_etc.csv


code[1:5,]

cg[1:5,]
keep1<-code # code<-keep1


cyto<-read.delim("/media/UQCCG/Sequencing/CompleteGenomics/20140918_QCTB-UQDI_AML_Cytogenetics_etc.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
#cg.done<-read.delim("/media/UQCCG/Sequencing/CompleteGenomics/917_Data_Delivery_091214.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
cyto[1:5,]


code[1:5,]

#posns<-match(code[,"exomes.ids"],cyto[,"ID"])
posns<-match(code[,"customer_sample_id"],cyto[,"ID"])
posns
missing<-is.na(posns)
sum(missing)

code[!missing,"Blast"]<-cyto[posns[!missing],"Bone.marrow.blast....morphology."]
code[!missing,"FLT3-ITD"]<-cyto[posns[!missing],"FLT3.ITD"]
code[!missing,"FLT3.ITD.allelic.ratio"]<-cyto[posns[!missing],"FLT3.ITD.allelic.ratio"]

code[!missing,"FLT3.D835"]<-cyto[posns[!missing],"FLT3.D835.mutation"]
code[!missing,"FLT3.D835H.Y"]<-"as D835"
code[!missing,"FLT3.D835G.V"]<-"as D835"
code[!missing,"FLT3.D835E"]<-"as D835"
code[!missing,"FLT3.I836DEL"]<-"not tested"
code[!missing,"FLT3.I836INS"]<-"not tested"


code[!missing,"NPM1"]<-cyto[posns[!missing],"NPM1.Mutation"]
code[!missing,"PATH-main"]<-cyto[posns[!missing],"Main.Cytogenetic.Feature"]
code[!missing,"PATH-full"]<-cyto[posns[!missing],"Full.Pathology.Cytogenetics"]
code[1:15,]
sum(missing)
dim(cg.done)
sum(!missing)

#####################################

keep2<-code
cyto<-read.delim("/media/UQCCG/Sequencing/CompleteGenomics/Final 51 AMLM12 Pts - with seq+FLt3+ Blast and Cyto 140929.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
#cg.done<-read.delim("/media/UQCCG/Sequencing/CompleteGenomics/917_Data_Delivery_091214.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
cyto[1:5,]


code[1:5,]

#posns<-match(code[,"exomes.ids"],cyto[,"Registration.Number"])
posns<-match(code[,"customer_sample_id"],cyto[,"Registration.Number"])
posns
missing<-is.na(posns)
sum(missing)

code[!missing,"Blast"]<-cyto[posns[!missing],"Blast.."]
code[!missing,"FLT3-ITD"]<-cyto[posns[!missing],"FLT3.ITD"]
code[!missing,"FLT3.ITD.allelic.ratio"]<-cyto[posns[!missing],"Ratio"]

code[!missing,"FLT3.D835"]<-"as components"
code[!missing,"FLT3.D835H.Y"]<-cyto[posns[!missing],"D835H.Y"]
code[!missing,"FLT3.D835G.V"]<-cyto[posns[!missing],"D835G.V"]
code[!missing,"FLT3.D835E"]<-cyto[posns[!missing],"D835E"]
code[!missing,"FLT3.I836DEL"]<-cyto[posns[!missing],"I836DEL"]
code[!missing,"FLT3.I836INS"]<-cyto[posns[!missing],"I836INS"]


code[!missing,"NPM1"]<-cyto[posns[!missing],"W288C"]
code[!missing,"PATH-main"]<-cyto[posns[!missing],"Cytogenetics"]
code[!missing,"PATH-full"]<-cyto[posns[!missing],"Cytogenetics"]
code[1:15,]
sum(missing)
dim(cg.done)
sum(!missing)


#####################################

keep3<-code
cyto<-read.delim("/media/UQCCG/Sequencing/CompleteGenomics/CREST_Analysis_FLT3_exome.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
#cg.done<-read.delim("/media/UQCCG/Sequencing/CompleteGenomics/917_Data_Delivery_091214.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
cyto[1:5,]


code[1:5,]

posns<-match(code[,"exomes.ids"],cyto[,"c"])
#posns<-match(code[,"customer_sample_id"],cyto[,"c"])
posns
missing<-is.na(posns)
sum(missing)

code[!missing,"Blast"]<-cyto[posns[!missing],"BM.Blast.."]
code[!missing,"FLT3-ITD"]<-cyto[posns[!missing],"FLT3.ITD"]
code[!missing,"FLT3.ITD.allelic.ratio"]<-"NA"

code[!missing,"FLT3.D835"]<-"as components"
code[!missing,"FLT3.D835H.Y"]<-cyto[posns[!missing],"D835H.Y"]
code[!missing,"FLT3.D835G.V"]<-cyto[posns[!missing],"D835G.V"]
code[!missing,"FLT3.D835E"]<-cyto[posns[!missing],"D835E"]
code[!missing,"FLT3.I836DEL"]<-cyto[posns[!missing],"I836DEL"]
code[!missing,"FLT3.I836INS"]<-cyto[posns[!missing],"I836INS"]


code[!missing,"NPM1"]<-cyto[posns[!missing],"W288C"]
code[!missing,"PATH-main"]<-cyto[posns[!missing],"DIAGNOSIS"]
code[!missing,"PATH-full"]<-cyto[posns[!missing],"cytogenetics"]
code[1:15,]
sum(missing)
dim(cg.done)
sum(!missing)

getwd()

setwd("/media/UQCCG/Sequencing/CompleteGenomics")
write.table(code,file="core_Cytogenics.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE,,append=TRUE)


code<-read.delim("/media/UQCCG/Sequencing/CompleteGenomics/core.tytgenetics.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
code[1:5,]

coverage[1:5,]

sum(coverage[,"Sample"] %in% code[,"exomes.ids"])
aml.in.code<-coverage[,"Sample"] %in% code[,"exomes.ids"]

mean( coverage[aml.in.code,"mean.mean.coverage"])
range( coverage[aml.in.code,"mean.mean.coverage"])

mean( coverage[aml.in.code,"percent.ccds.gt.10"])
range( coverage[aml.in.code,"percent.ccds.gt.10"])


cyto<-read.delim("/media/UQCCG/Sequencing/CompleteGenomics/CREST_Analysis_FLT3_exome.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

code[1:3,]
cg[1:1,1:20]

code[,"assembly_id"]
the.colnames<-colnames(cg)
the.colnames<-gsub("\\.","-",the.colnames)

posns<-match(the.colnames,code[,"assembly_id"]) # posns<-match(the.colnames,code[,"exomes.ids"])
missing<-is.na(posns)
posns
sum(missing)
the.colnames[missing]
length(the.colnames)-20
dim(code)
the.colnames[!missing]<-code[posns[!missing],"exomes.ids"] # the.colnames[!missing]<-code[posns[!missing],"customer_sample_id"]

code[posns[!missing],1:3]

colnames(cg)<-the.colnames

order.wanted<-c(colnames(cg)[1:20],code[,"exomes.ids"])

wanted<-grepl("GATA1",cg[,"geneVar"]) #& grepl("COSMIC",cg[,"geneVar"]) &
wanted<- grepl("COSMIC",cg[,"xRef"]) & cg[,"1409201"]!="0" & cg[,"1409201"]!="00" & cg[,"1409201"]!="NN" & !grepl(":SYNONYMOUS$",cg[,"geneVar"]) & !grepl("dbsnp",cg[,"xRef"])
sum(wanted)
getwd()
cg[wanted,][1:10,c(colnames(cg)[1:20],"1409201")]
write.table(cg[wanted,],file="1409201_all.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE,,append=TRUE)

write.table(cg[wanted,c(colnames(cg)[1:20],"1409201")],file="1409201_all.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE,,append=TRUE)

length(order.wanted)
dim(cg)
cg<-cg[,order.wanted]
colnames(cg)
code[,"exomes.ids"]
the.samples<-code[,"exomes.ids"]

## code[grep("GS000038096",code[,"assembly_id"]),"assembly_id"]
## the.colnames[grep("GS000038096",the.colnames)]

cg.gr<-GRanges(seqnames =cg[,"chromosome"],ranges = IRanges(start=as.numeric(cg[,"begin"]),end=as.numeric(cg[,"end"])),strand="+")


### NPM1

target<-grepl("NPM1",cg[,"geneVar"])
sum(target)
cg[target,c(1:8,20)]

the.chr<-"chr5"
begin<-170837513
end<-170837670

the.chr<-"chr5"
begin<-170837542
end<-170837646

#########

### FLT3
colnames(cg)
target<-grepl("FLT3",cg[,"geneVar"])
sum(target)
cg[target,c(1:8,20)]
cg[target,"geneVar"]

the.chr<-"chr13"
begin<-28592639
end<-28592642

the.chr<-"chr13"
begin<-28608224
end<-28608267
end<-28608265
#########

indel.gr<-GRanges(seqnames =the.chr,ranges = IRanges(start=as.numeric(begin),end=as.numeric(end)),strand="+")

common<-overlapsAny(cg.gr,indel.gr)
sum(common)
cg[common,c(1:8,20)]

wanted<-c(colnames(cg)[1:6],"9")

cg[common,wanted]

hits<-apply(cg[common,the.samples],2,function(x) sum(x!="00"))
sort(hits)
hits<-hits[code[,"exomes.ids"]]
hits



########

code[,"NPM1-CG"]<-hits
code[1:5,]

code[,"FLT3-ITD-CG"]<-hits
code[1:5,]

code[,"FLT3-D835-CG"]<-hits


write.table(code,file="core.tytgenetics.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)









###########################################################

a.indel.gr<-GRanges(seqnames =a.indel[,"chr"],ranges = IRanges(start=as.numeric(a.indel[,"start"]),end=as.numeric(a.indel[,"end"])),strand="+")

is.flat.geno<-grepl(":flat$",a.indel[,"TYPE"])
code[1:5,]
getwd() #  "/media/UQCCG/Sequencing/CompleteGenomics"


colnames(a.indel)
target<-grepl("FLT3",a.indel[,"Gene.Names"])
sum(target)
a.indel[target,c(1:6)]


the.chr<-"chr13"
begin<-28608214
end<-28608336

the.chr<-"chr5"
begin<-170837522
end<-170837676

#########

AMLM12005R-G

grep("AMLM12005R-G",all.possible.samples)
the.samples<-aml.samples
all.possible.samples 
the.samples.HC<-all.possible.samples[all.possible.samples %in% the.samples]
the.samples.HC.GT<-paste(the.samples.HC,"GT",sep=".")
the.samples.HC.GQ<-paste(the.samples.HC,"GQ",sep=".")
indel.gr<-GRanges(seqnames =the.chr,ranges = IRanges(start=as.numeric(begin),end=as.numeric(end)),strand="+")

common<-overlapsAny(a.indel.gr,indel.gr)
sum(common)
a.indel[common & !is.flat.geno,1:9]
a.indel[common & !is.flat.geno,the.samples.HC.GT][1:2,]
a.indel[common & !is.flat.geno,the.samples.HC.GT][,1]
a.indel[common & !is.flat.geno,the.samples.HC.GQ]
hits<-apply(a.indel[common & !is.flat.geno,the.samples.HC.GT],2,function(x) sum(x!="0/0"))

x<-a.indel[common & !is.flat.geno,c(the.samples.HC.GT,the.samples.HC.GQ)][1,]
hits<-apply(a.indel[common & !is.flat.geno,c(the.samples.HC.GT,the.samples.HC.GQ)],1,function(x){
 # print(x)
  a.geno<-x!="0/0" & grepl("GT$",names(x))
  samples<-names(x)[a.geno]
  sample.GQ<-gsub(".GT",".GQ",samples)
  paste(x[c(samples,sample.GQ)],collapse="::")
})
  sum(hits!=0)

names(hits)<-gsub(".GT$","",names(hits))
hits<-hits[code[,"exomes.ids"]]
hits
sort(hits)
colnames(code)
code[,"FLT3-ITD-HC"]<-hits
code[,"NPM1-HC"]<-hits
code[1:5,]
write.table(code,file="core.tytgenetics.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
colnames(a.indel)[1:36]
wanted<-c(colnames(a.indel)[c(1:6,34)],"AMLM12005R-G.GT","AMLM12005R-G.GQ")
a.indel[common & !is.flat.geno,wanted]
a.indel[common & !is.flat.geno,][1,the.samples.HC.GT]

wanted<-c(colnames(a.indel)[1:69],the.samples.HC.GT,colnames(a.indel)[1279:1317])
write.table(a.indel[common & !is.flat.geno,wanted],file="FLT3.GENOTYPES.in.HC.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


load("/media/UQCCG/Sequencing/Data/Genomes/hg19/NexteraRapidCapture_Exome_TargetedRegions_hg19_targets.RData")

dim(cg)
dim(a.indel)

dim(cg)
common<-overlapsAny(cg.gr,data.gr)
sum(common)
cg.sub<-cg[common,]
dim(cg.sub)


common<-overlapsAny(cg.gr,data.gr)
sum(common)
cg.sub<-cg[common,]

all.possible.samples 
the.samples.HC<-all.possible.samples[all.possible.samples %in% the.samples]
the.samples.HC.GT<-paste(the.samples.HC,"GT",sep=".")


dim(a.indel)
common<-overlapsAny(a.indel.gr,data.gr)
is.flat.geno<-grepl(":flat$",a.indel[,"TYPE"])


#pass<- full.qual & !no.genotypes.AML  & !in.common.hit.gene  & !are.repeats & !are.in.repeats & !is.unwound.geno & ok.missing & ok.missing.filt  & hw.controls.ok.filt & !no.genotypes.filt
pass<- full.qual & !no.genotypes.AML   & !are.repeats & !are.in.repeats  & ok.missing & ok.missing.filt  & hw.controls.ok.filt & !no.genotypes.filt


sum(common & !is.flat.geno)

sum(common & !is.flat.geno & pass)

sum(common & !is.flat.geno & full.qual & !no.genotypes.AML)

sum(common & !is.flat.geno & full.qual & !no.genotypes.AML  & ok.missing & ok.missing.filt )

a.indel.sub<-a.indel[common & !is.flat.geno,]
dim(a.indel.sub)
a.indel.sub[1:5,1:10]

cg.sub[1:5,]
dim(cg.sub)
unique(cg.sub[,"varType"])
a.cg.snp<-cg.sub[,"varType"]=="snp"
sum(!a.cg.snp)
cg.sub.indel<-cg.sub[!a.cg.snp,]
cg.sub.indel[1:5,]
unique(cg.sub.indel[,55]) #[1] "00" "NN" "11" "0N" "01" "1N" "N" 
 place<-      grep("1",cg.sub.indel[,50])
 dim(cg.sub.indel)      
       sum(place)
cg.sub.indel["892376",]

x<-cg.sub.indel[3,the.samples]
 x      
cg.miss<-apply(cg.sub.indel[1:5,the.samples],1,function(x) sum(x=="NN")/sum(x=="01" | x=="1" |  x=="11" )   )
cg.miss<-apply(cg.sub.indel[,the.samples],1,function(x) sum(x=="NN" | x=="N")/sum(x!="00" | x!="0N")   )

cg.miss<-apply(cg.sub.indel[,the.samples],1,function(x) sum(x=="NN" | x=="N")/sum(x!="00")   )

hist(cg.miss)
cg.miss[1:50]
cg.miss.ok<-cg.miss<=0.1
sum(cg.miss.ok)
cg.sub.indel.pass<-cg.sub.indel[cg.miss.ok,]
dim(cg.sub.indel.pass)


cg.sub.indel.pass.gr<-GRanges(seqnames =cg.sub.indel.pass[,"chromosome"],ranges = IRanges(start=as.numeric(cg.sub.indel.pass[,"begin"]),end=as.numeric(cg.sub.indel.pass[,"end"])),strand="+")
a.indel.sub.gr<-GRanges(seqnames =a.indel.sub[,"chr"],ranges = IRanges(start=as.numeric(a.indel.sub[,"start"]),end=as.numeric(a.indel.sub[,"end"])),strand="+")

cg.sub.indel.pass.gr<-cg.sub.indel.pass.gr+10
a.indel.sub.gr<-a.indel.sub.gr+10
common<-overlapsAny(cg.sub.indel.pass.gr,a.indel.sub.gr)

sum(common)
length(a.indel.sub.gr)
length(cg.sub.indel.pass.gr)


> sum(common)
[1] 3281
> length(a.indel.sub.gr)
[1] 13097
> length(cg.sub.indel.pass.gr)
[1] 20638

exome.gt.ref<-a.indel.sub[,"AMAS-9.3-Diagnostic.GT"]=="0/0" | a.indel.sub[,"AMAS-9.3-Diagnostic.GT"]=="NA" 
sum(!exome.gt.ref)
length(!exome.gt.ref)
unique(cg.sub.indel.pass[,"AMAS-9.3-Diagnostic"])
cg.gt.ref<-cg.sub.indel.pass[,"AMAS-9.3-Diagnostic"]=="00" | cg.sub.indel.pass[,"AMAS-9.3-Diagnostic"]=="0N" 
sum(!cg.gt.ref)
length(cg.gt.ref)


cg.sub.indel.pass.gr<-GRanges(seqnames =cg.sub.indel.pass[!cg.gt.ref,"chromosome"],ranges = IRanges(start=as.numeric(cg.sub.indel.pass[!cg.gt.ref,"begin"]),end=as.numeric(cg.sub.indel.pass[!cg.gt.ref,"end"])),strand="+")
a.indel.sub.gr<-GRanges(seqnames =a.indel.sub[!exome.gt.ref,"chr"],ranges = IRanges(start=as.numeric(a.indel.sub[!exome.gt.ref,"start"]),end=as.numeric(a.indel.sub[!exome.gt.ref,"end"])),strand="+")

cg.sub.indel.pass.gr<-cg.sub.indel.pass.gr+10
a.indel.sub.gr<-a.indel.sub.gr+10
common<-overlapsAny(cg.sub.indel.pass.gr,a.indel.sub.gr)

sum(common)
length(a.indel.sub.gr)
length(cg.sub.indel.pass.gr)

length(cg.sub.indel.pass.gr)
> [1] 494
> [1] 1089
> [1] 1181

dim(cg.sub.indel.pass[!cg.gt.ref,])
dim(a.indel.sub)

length(common)
cg.sub.indel.pass
common<-overlapsAny(a.indel.sub.gr,cg.sub.indel.pass.gr)
sum(common)
length(a.indel.sub.gr)
length(cg.sub.indel.pass.gr)
median(as.numeric(a.indel.sub[!exome.gt.ref,][!common,"AMAS-9.3-Diagnostic.GQ"]))

hist(as.numeric(a.indel.sub[!exome.gt.ref,][!common,"AMAS-9.3-Diagnostic.GQ"]))
sum(as.numeric(a.indel.sub[!exome.gt.ref,][common,"AMAS-9.3-Diagnostic.GQ"])<20)



##                chr13:28608214:28608214:-:TTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCA:indel                               chr13:28608228:28608228:-:ATTTTCTCTTGGAAACTCCCATTTGAGATCATATTC:indel 
##                                                                                          "1/1::39"                                                                                          "0/1::99" 
##                                  chr13:28608230:28608230:-:TTTCTCTTGGAAACTCCCATTTGAGAGGGTATC:indel                               chr13:28608231:28608231:-:TTCTCTTGGAAACTCCCATTTGAGATCATATTCTTG:indel 
##                                                                                          "0/1::99"                                                                                                 "" 
##                                  chr13:28608239:28608239:-:GAAACTCCCATTTGAGATCATATTCAGGGGTCC:indel                chr13:28608242:28608242:-:ACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGCCATTGGTAG:indel 
##                                                                                                 ""                                                                                          "0/1::99" 
##                                           chr13:28608249:28608249:-:TTTGAGATCATATTCATATTCTCT:indel                         chr13:28608250:28608250:-:TTGAGATCATATTCATATTCTCTGAAATCAACG:indel:28608250 
##                                                                                                 ""                                                                                                 "" 
##                            chr13:28608250:28608250:-:TTGAGATCATATTCATATTCTCTGAAATCC:indel:28608250       chr13:28608252:28608252:-:GAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGG:indel 
##                                                                                                 ""                                                                                          "0/1::99" 
##                                           chr13:28608254:28608254:-:GATCATATTCATATTCTCTGAAAT:indel                                              chr13:28608255:28608255:-:ATCATATTCATATTCTCTGAA:indel 
##                                                                                          "0/1::99"                                                                                                 "" 
##                                                 chr13:28608257:28608257:-:CATATTCATATTCTCTGA:indel                   chr13:28608260:28608260:-:ATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGC:indel 
##                                                                                                 ""                                                                                          "0/1::99" 
##                   chr13:28608261:28608261:-:TTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCC:indel                                        chr13:28608262:28608262:-:TCATATTCTCGGGAATAA:indel:28608262 
##                                                                                                 ""                                                                                                 "" 
##                         chr13:28608262:28608262:-:TCATATTCTCTGAAATCAACGTAGAAGTACTCA:indel:28608262                                  chr13:28608262:28608262:-:TCATATTCTCTGAAATCAACGTAG:indel:28608262 
##                                                                                          "0/1::99"                                                                                          "0/1::99" 
## chr13:28608273:28608273:-:GAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCTGTACCATCTGTAGCTGGCTTTC:indel                               chr13:28608274:28608274:-:AAATCAACGTAGAAGTACTCATTATCTGAGGAGCCT:indel 
##                                                                                                 ""                                                                                          "0/1::99" 
##                                                  chr13:28608275:28608275:-:TTGCGTTC:indel:28608275                                                       chr13:28608275:28608275:-:CCT:indel:28608275 
##                                                                                                 ""                                                                                          "0/1::99" 
##                               chr13:28608280:28608280:-:ACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACC:indel                                                        chr13:28608310:28608310:-:ACATATTCTCT:indel 
##                                                                                                 ""                                                                                          "0/1::99" 
##                                         chr13:28608311:28608311:-:AGTCAAACTCTAAATTT:indel:28608311                                        chr13:28608311:28608311:-:AAATCAACGTAGAAGTAC:indel:28608311 
##                                                                                                 ""                                                                                          "0/1::99" 
##                                                    chr13:28608336:28608336:-:CTCATCACATCCTAA:indel 
##                                                                                                 ""


                                                                                                   FILTER                          AMLM12037SAT.GT AMLM12037SAT.GQ
chr13:28608214:28608214:-:TTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCA:indel                "PASS"                          "0/0"           "0"            
chr13:28608228:28608228:-:ATTTTCTCTTGGAAACTCCCATTTGAGATCATATTC:indel                               "VQSRTrancheINDEL99.00to99.90"  "0/1"           "99"           
chr13:28608230:28608230:-:TTTCTCTTGGAAACTCCCATTTGAGAGGGTATC:indel                                  "VQSRTrancheINDEL99.00to99.90"  "0/0"           "99"           
chr13:28608231:28608231:-:TTCTCTTGGAAACTCCCATTTGAGATCATATTCTTG:indel                               "VQSRTrancheINDEL99.00to99.90"  "0/0"           "99"           
chr13:28608239:28608239:-:GAAACTCCCATTTGAGATCATATTCAGGGGTCC:indel                                  "PASS"                          "0/0"           "99"           
chr13:28608242:28608242:-:ACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGCCATTGGTAG:indel                "VQSRTrancheINDEL99.00to99.90"  "0/0"           "99"           
chr13:28608249:28608249:-:TTTGAGATCATATTCATATTCTCT:indel                                           "VQSRTrancheINDEL99.00to99.90"  "0/0"           "99"           
chr13:28608250:28608250:-:TTGAGATCATATTCATATTCTCTGAAATCAACG:indel:28608250                         "VQSRTrancheINDEL99.00to99.90"  "0/0"           "99"           
chr13:28608250:28608250:-:TTGAGATCATATTCATATTCTCTGAAATCC:indel:28608250                            "VQSRTrancheINDEL99.00to99.90"  "0/0"           "99"           
chr13:28608252:28608252:-:GAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGG:indel       "VQSRTrancheINDEL99.00to99.90"  "0/0"           "99"           
chr13:28608254:28608254:-:GATCATATTCATATTCTCTGAAAT:indel                                           "VQSRTrancheINDEL99.00to99.90"  "0/0"           "99"           
chr13:28608255:28608255:-:ATCATATTCATATTCTCTGAA:indel                                              "VQSRTrancheINDEL99.00to99.90"  "0/0"           "99"           
chr13:28608257:28608257:-:CATATTCATATTCTCTGA:indel                                                 "VQSRTrancheINDEL99.00to99.90"  "0/0"           "99"           
chr13:28608260:28608260:-:ATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGC:indel                   "VQSRTrancheINDEL99.00to99.90"  "0/0"           "99"           
chr13:28608261:28608261:-:TTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCC:indel                   "VQSRTrancheINDEL99.00to99.90"  "0/0"           "99"           
chr13:28608262:28608262:-:TCATATTCTCGGGAATAA:indel:28608262                                        "VQSRTrancheINDEL99.00to99.90"  "0/0"           "99"           
chr13:28608262:28608262:-:TCATATTCTCTGAAATCAACGTAGAAGTACTCA:indel:28608262                         "VQSRTrancheINDEL99.00to99.90"  "0/0"           "99"           
chr13:28608262:28608262:-:TCATATTCTCTGAAATCAACGTAG:indel:28608262                                  "VQSRTrancheINDEL99.00to99.90"  "0/0"           "99"           
chr13:28608273:28608273:-:GAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCTGTACCATCTGTAGCTGGCTTTC:indel "PASS"                          "0/0"           "62"           
chr13:28608274:28608274:-:AAATCAACGTAGAAGTACTCATTATCTGAGGAGCCT:indel                               "VQSRTrancheINDEL99.00to99.90"  "0/0"           "62"           
chr13:28608275:28608275:-:TTGCGTTC:indel:28608275                                                  "VQSRTrancheINDEL99.00to99.90"  "0/0"           "62"           
chr13:28608275:28608275:-:CCT:indel:28608275                                                       "VQSRTrancheINDEL99.00to99.90"  "0/0"           "62"           
chr13:28608280:28608280:-:ACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACC:indel                               "VQSRTrancheINDEL99.00to99.90"  "0/0"           "52"           
chr13:28608310:28608310:-:ACATATTCTCT:indel                                                        "VQSRTrancheINDEL99.00to99.90"  "0/0"           "70"           
chr13:28608311:28608311:-:AGTCAAACTCTAAATTT:indel:28608311                                         "VQSRTrancheINDEL99.00to99.90"  "0/0"           "70"           
chr13:28608311:28608311:-:AAATCAACGTAGAAGTAC:indel:28608311                                        "VQSRTrancheINDEL99.00to99.90"  "0/0"           "70"           
chr13:28608336:28608336:-:CTCATCACATCCTAA:indel                                                    "VQSRTrancheINDEL99.90to100.00" "0/0"           "70"           
> hits<-apply(a.indel[common & !is.flat.geno,c(the.samples.HC.GT,the.samples.HC.GQ)],1,function(x){
+   a.geno<-x!="0/0" & grepl("GT$",names(x))
  samples<-names(x)[a.geno]
+ +   sample.GQ<-gsub(".GT",".GQ",samples)
+   paste(x[c(samples,sample.GQ)],collapse="::")
})
+ > hits
               chr13:28608214:28608214:-:TTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCA:indel                               chr13:28608228:28608228:-:ATTTTCTCTTGGAAACTCCCATTTGAGATCATATTC:indel 
                                                                                         "1/1::39"                                                                                          "0/1::99" 
                                 chr13:28608230:28608230:-:TTTCTCTTGGAAACTCCCATTTGAGAGGGTATC:indel                               chr13:28608231:28608231:-:TTCTCTTGGAAACTCCCATTTGAGATCATATTCTTG:indel 
                                                                                         "0/1::99"                                                                                                 "" 
                                 chr13:28608239:28608239:-:GAAACTCCCATTTGAGATCATATTCAGGGGTCC:indel                chr13:28608242:28608242:-:ACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGCCATTGGTAG:indel 
                                                                                                ""                                                                                          "0/1::99" 
                                          chr13:28608249:28608249:-:TTTGAGATCATATTCATATTCTCT:indel                         chr13:28608250:28608250:-:TTGAGATCATATTCATATTCTCTGAAATCAACG:indel:28608250 
                                                                                                ""                                                                                                 "" 
                           chr13:28608250:28608250:-:TTGAGATCATATTCATATTCTCTGAAATCC:indel:28608250       chr13:28608252:28608252:-:GAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGG:indel 
                                                                                                ""                                                                                          "0/1::99" 
                                          chr13:28608254:28608254:-:GATCATATTCATATTCTCTGAAAT:indel                                              chr13:28608255:28608255:-:ATCATATTCATATTCTCTGAA:indel 
                                                                                         "0/1::99"                                                                                                 "" 
                                                chr13:28608257:28608257:-:CATATTCATATTCTCTGA:indel                   chr13:28608260:28608260:-:ATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGC:indel 
                                                                                                ""                                                                                          "0/1::99" 
                  chr13:28608261:28608261:-:TTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCC:indel                                        chr13:28608262:28608262:-:TCATATTCTCGGGAATAA:indel:28608262 
                                                                                                ""                                                                                                 "" 
                        chr13:28608262:28608262:-:TCATATTCTCTGAAATCAACGTAGAAGTACTCA:indel:28608262                                  chr13:28608262:28608262:-:TCATATTCTCTGAAATCAACGTAG:indel:28608262 
                                                                                         "0/1::99"                                                                                          "0/1::99" 
chr13:28608273:28608273:-:GAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCTGTACCATCTGTAGCTGGCTTTC:indel                               chr13:28608274:28608274:-:AAATCAACGTAGAAGTACTCATTATCTGAGGAGCCT:indel 
                                                                                                ""                                                                                          "0/1::99" 
                                                 chr13:28608275:28608275:-:TTGCGTTC:indel:28608275                                                       chr13:28608275:28608275:-:CCT:indel:28608275 
                                                                                                ""                                                                                          "0/1::99" 
                              chr13:28608280:28608280:-:ACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACC:indel                                                        chr13:28608310:28608310:-:ACATATTCTCT:indel 
                                                                                                ""                                                                                          "0/1::99" 
                                        chr13:28608311:28608311:-:AGTCAAACTCTAAATTT:indel:28608311                                        chr13:28608311:28608311:-:AAATCAACGTAGAAGTAC:indel:28608311 
                                                                                                ""                                                                                          "0/1::99" 
                                                   chr13:28608336:28608336:-:CTCATCACATCCTAA:indel 
                                                                                                "" 
