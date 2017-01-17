






#answer = which.isMatchingAt("N", seqs, at=6:1, follow.index=TRUE)

code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir) # load in the UCSC tables these use the db file names and not their lable-names 
source("annotate_SNPs_subroutines.r")


options("width"=250,"max.print"=1000)

vcf.files<-"00-All.vcf"
output.name<-"hg19_snp137"

vcf.files<-"clinvar_20120616_2.vcf"
output.name<-"hg19_snp137_clinical"




names(vcf.files)<-"snp"
snp.dir<-"/media/scratch2/dbSNP"


 # .txt extension will be added


location<-snp.dir
unwanted.cols<-c("QUAL","FILTER")

########################BEGIN
setwd(location)
the.files<-dir(getwd())

if(paste(output.name,"_maf.txt",sep="") %in% the.files){print("WARNING output files exits they will be appended too!!")}

################ LARGE FILES ########
  print(vcf.files)
######BELOW PROCESSING  this for snp for indel varient types in vcf4.0 or vcf 3.0 format
rm(start.data)
start.data<-prepare.for.Vcf.file.read(vcf.files)
for(i in 1:length(start.data)){assign( names(start.data)[i],value=start.data[[i]])}

cbind(info.types,info.description)


con <- file(vcf.files, open="r")  # close(con)
num.lines<-1 # so does while llop at least once
reads<-1000000  #1.3M lines in snp file 50000 goes out to 24Gb without QC cgeck 
counter<- -1
while (num.lines >0  ){
counter<-counter+1
print(counter)

if(counter==0){
indels<-try(scan(con,what=character(num.vars),skip=(reads*counter)+skip.lines,nlines=reads,sep="\t",fill=TRUE,na.strings="",quote="\""))
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




#indels[1:5,]

##PM,Number=0,Type=Flag,Description="Variant is Precious(Clinical,Pubmed Cited)">
##INFO=<ID=TPA,Number=0,Type=Flag,Description="Provisional Third Party Annotation(TPA) (currently rs from PHARMGKB who will give phenotype data)">
 ##INFO=<ID=CDA,Number=0,Type=Flag,Description="Variation is interrogated in a clinical diagnostic assay">

####################################### FINISHED Read in data DO PROCESSIng below
###################################################################################################

alt.list<-strsplit(indels[,"ALT"],split=",")

has.g5.or.g5A<-grepl(";G5",indels[,"INFO"]) # catch G5 and G5A
## has.g5<-grepl(";G5;",indels[,"INFO"])
## has.g5a<-grepl(";G5A;",indels[,"INFO"])
has.gmaf<-grepl(";GMAF=",indels[,"INFO"])
#########extract A quatity from INFO GMAF DONE HERE

#indels[has.gmaf,"INFO"][1:20]

the.gmaf<-extract.value.from.info(indels[has.gmaf,"INFO"],"GMAF=")

###########################################


# make MAF 
the.af<-rep(0,times=dim(indels)[1]) 
the.af[!has.gmaf & has.g5.or.g5A]<-0.5  ### maf > 0.5% so set to 0.5 so it is excluded entirely
the.af[has.gmaf]<-the.gmaf

  
################### fltten the data and previos tests

number.of.alleles<-unlist(lapply(alt.list,length))
flat.index<-rep(1:length(number.of.alleles),times=number.of.alleles)
indels<-indels[flat.index,]
the.af<-the.af[flat.index]
has.g5.or.g5A<-has.g5.or.g5A[flat.index]
indels[,"ALT"]<-unlist(alt.list) # they are unlisted in the same order as they appear


#################### FLAG interesting stuff


has.OM<-grepl(";OM",indels[,"INFO"])
has.GNO<-grepl(";GNO",indels[,"INFO"])
has.CLN<-grepl(";CLN",indels[,"INFO"]) # variant is clinical (but of that is benighn - taste)
has.CDA<-grepl(";CDA",indels[,"INFO"]) #  "Variation is interrogated in a clinical diagnostic assay"
has.PM<-grepl(";PM",indels[,"INFO"]) #  ""Variant is Precious(Clinical,Pubmed Cited)""
has.MUT<-grepl(";MUT",indels[,"INFO"]) #  Is mutation (journal citation, explicit fact): a low frequency variation that is cited in journal and other reputable sources">
has.SCS<-grepl(";SCS",indels[,"INFO"]) # variant is clinical significance
#0 - unknown, 1 - untested, 2 - non-pathogenic, 3 - probable-non-pathogenic, 4 - probable-pathogenic, 5 - pathogenic, 6 - drug-response, 7 - histocompatibility, 255 - other">

## has.SCS.prob.patho<-grepl(";SCS=4;",indels[,"INFO"])
## has.SCS.patho<-grepl(";SCS=5;",indels[,"INFO"])
## has.SCS.drug<-grepl(";SCS=6;",indels[,"INFO"])
## has.SCS.histo<-grepl(";SCS=7;",indels[,"INFO"])
## has.bad.path<-(has.SCS.prob.patho | has.SCS.patho | has.SCS.drug | has.SCS.histo)
## to.test<- has.g5.or.g5A & has.bad.path
## if(sum(to.test)>0){print("common pathogenic allele");print(indels[to.test,])} # these appear to be GWAS hits


## pathological<-vector(mode="character",length=dim(indels)[1])
## clinical<-vector(mode="character",length=dim(indels)[1])

## pathological[has.SCS.patho]<-"pathogenic"
## pathological[has.SCS.prob.patho]<-"probable-pathogenic"
## pathological[has.SCS.drug]<-"drug-response"
## pathological[has.SCS.histo]<-"histocompatibility"

## clinical.event<-extract.value.from.info(indels[has.CLN,"INFO"],"CLN=")
## clinical[has.CLN]<-clinical.event


########################### PUT ALLELES IN ANNOVAR FORMAT from convert2annovar.pl line 1083########
ref.length<-nchar(as.character(indels[,"REF"]))
alt.length<-nchar(as.character(indels[,"ALT"]))
is.snp<-(ref.length==1 & alt.length==1)

POS.end<-as.numeric(indels[,"POS"])
del<-ref.length > alt.length
ins<-(ref.length <= alt.length) & !is.snp
#indels[del,][1:5,]
#POS.end[del][1:5]
### deletion or block substitution
head<-substr(as.character(indels[del,"REF"]),1,alt.length[del])
head.is.mut<-(head==as.character(indels[del,"ALT"]))
indels[del,"REF"][head.is.mut]<-substr(as.character(indels[del,"REF"][head.is.mut]),(alt.length[del][head.is.mut]+1),ref.length[del][head.is.mut])
indels[del,"ALT"][head.is.mut]<-"-"
indels[del,"POS"][head.is.mut]<-as.numeric(indels[del,"POS"][head.is.mut]) + nchar(as.character(head[head.is.mut]))
POS.end[del]<-POS.end[del]+ref.length[del]-1  # same for both head is mut and not head is mut

## indels
## POS.end
### insertion or block substitution
head<-substr(as.character(indels[ins,"ALT"]),1,ref.length[ins])
head.is.ref<-(head==as.character(indels[ins,"REF"]))
indels[ins,"ALT"][head.is.ref]<-substr(as.character(indels[ins,"ALT"][head.is.ref]),(ref.length[ins][head.is.ref]+1),alt.length[ins][head.is.ref])
indels[ins,"REF"][head.is.ref]<-"-"
indels[ins,"POS"][head.is.ref]<-as.numeric(indels[ins,"POS"][head.is.ref]) + ref.length[ins][head.is.ref]-1
POS.end[ins]<-POS.end[ins]+ref.length[ins]-1
########################################################################################################
#indels[ins,]
#POS.end[ins]
###wait here

indels<-cbind(indels[,c("chr","POS")],POS.end,indels[,c("REF","ALT")],the.af,indels[,c("ID","INFO")])
#indels[1:5,]

write.table(indels,file=paste(output.name,"_maf.txt",sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=TRUE)

## write.table(indels[has.bad.path,],file=paste(output.name,"_pathalogical_maf.txt",sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=TRUE)

write.table(indels[has.OM,],file=paste(output.name,"_omim_maf.txt",sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
write.table(indels[has.CDA,],file=paste(output.name,"_clinical_assay_maf.txt",sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
write.table(indels[has.MUT,],file=paste(output.name,"_mutation_maf.txt",sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
write.table(indels[has.PM,],file=paste(output.name,"_pubmed_maf.txt",sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
} ## loop over data chunks 


## Error in indels[del, "REF"][head.is.mut] <- substr(as.character(indels[del,  : 
##   NAs are not allowed in subscripted assignments
#This reads the vcf files stored as names in samples files an makes a corrending data object with the name provided
## > sample.files
##                              snp                            indel 
##   "SKDP-FAM-26_All_snps.raw.vcf" "SKDP-FAM-26_All_DINDEL.raw.vcf"


## http://www.ncbi.nlm.nih.gov/projects/SNP/docs/rs_attributes.html#gmaf
## Global minor allele frequency (MAF):  dbSNP is reporting the minor allele frequency for each rs included in  a default global population. Since this is being provided to distinguish common polymorphism from rare variants, the MAF is actually the second most frequent allele value. In other words, if there are 3 alleles, with frequencies of 0.50, 0.49, and 0.01, the MAF will be reported as 0.49. The current default global population is 1000Genome phase 1 genotype data from 1094 worldwide individuals, released in the May 2011 dataset.

## For example, refSNP page for rs222 reports: "MAF/MinorAlleleCount:G=0.262/330". This means that for rs222, minor allele is 'G' and has a frequency of 26.2% in the 1000Genome phase 1 population and that 'G' is observed 330 times in the sample population of 629 people (or 1258 chromosomes). 
## ###################################to a read in a lrage file

