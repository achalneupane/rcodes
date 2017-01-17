
################# recode GWAS to phenotype IDS:
  #working with
fam.template.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL_chr1.fam"
fam<-read.table(fam.template.file,header=F,fill=TRUE,stringsAsFactors=FALSE)


### exclude<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.recode.txt",header=F,sep=",",fill=TRUE,stringsAsFactors=FALSE)
related.to.remove<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.recode.txt",header=F,fill=TRUE,sep=",",stringsAsFactors=FALSE)
related.to.remove<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.txt",header=F,fill=TRUE,sep=",",stringsAsFactors=FALSE)
related.to.remove[1:5,]

/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/cluster_viz/excluded.samples.txt
evoker.remove<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/cluster_viz/excluded.samples.txt",header=F,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

dim(exclude)
dim(related.to.remove)
identical(exclude,related.to.remove) #TRUE

annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

fam.new<-read.table("/media/UQCCG/UQCCG-Projects/PAUL_LEO/AOGC exome chip core genotyping/all.lists_filtered-clean.f.fam",header=F,fill=TRUE,stringsAsFactors=FALSE)
# new checked genotypes in:"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP.list.re-checked.ED.chr1-19"
fam.new<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/BMD_EFF_STD_HIP.fam",header=F,fill=TRUE,stringsAsFactors=FALSE)
removes.new<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/remove_from.BMD_EFF_STD_HIP",header=F,fill=TRUE,stringsAsFactors=FALSE)


fam.new<-fam.new[!(fam.new[,1] %in% removes.new[,1]),]
dim(fam.new) ## this are samples actually used 

sum(related.to.remove[,1] %in% fam.new[,1])

fam.new[fam.new[,1] %in% related.to.remove[,1],]


missing[,1] %in% fam[,1]
sum(evoker.remove[,1] %in% fam.new[,1])
sum(!(evoker.remove[,1] %in% related.to.remove[,1])) ## all samples in evoker removed are in related.to.remove.txt
related.to.remove[1:5,]




recodes<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/pheno_gwas_id_missmatches.csv",,header=T,fill=TRUE,stringsAsFactors=FALSE)
recodes[1:5,]


dim(fam)
dim(fam.new)


sum(!(fam[,1] %in% fam.new[,1]))
fam[!(fam[,1] %in% fam.new[,1]),][1:10,]

recodes[1:5,]
recodes[,2] %in% fam[,1] # fam uses GWASID
missing<-recodes[!(recodes[,2] %in% fam[,1]),]
missing[,2] %in% fam[,1]
missing[,1] %in% fam[,1] ## missing ones are no present anywgere

recodes[,2] %in% exclude[,1] # exclude used GWASID
recodes[,2] %in% ann[,"PATIENT"]  # ann used GWASID


recodes[,1] %in% fam.new[,1] # fam.new (original data  uses PhenoID
recodes[1:5,]
# fam.new (original data  uses PhenoID

### convert from 
posns<-match(fam.new[,1],recodes[,"PhenoID"])
missing<-is.na(posns)
sum(missing)
fam.new[!missing,1]<-recodes[posns[!missing],"GWASID"]
fam.new[!missing,2]<-recodes[posns[!missing],"GWASID"]
sum(duplicated(fam.new[,1]))

## replace fam file with recoded
#setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP.list.re-checked.ED.chr1-19")
write.table(fam.new,file="all.lists_filtered.fam",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)


######################################################### single point for exome chip
######################################################### single point for exome chip
######################################################### single point for exome chip

fam.template.file<-"recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL"

fam.template.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL_chr1.fam"
annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"


setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/")
projects<-"recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL"

fam<-read.table(fam.template.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

exclude<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.recode.txt",header=F,sep=",",fill=TRUE,stringsAsFactors=FALSE)

related.to.remove<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.recode.txt",header=F,fill=TRUE,sep=",",stringsAsFactors=FALSE)
related.to.remove[1:5,]

annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)



########### Get Hardy for new tests

fam.template.file # set above
projects # set above

traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")

#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
names(traits)<-rep("logistic",times=length(traits))
names(traits)[1:3]<-rep("assoc",times=3)
traits
traits %in% colnames(ann)


fam[1:5,]


i<-1
for (i in 1:length(traits)){

  ### alredy got remove file from above
## write.table(fam[missing,c(1,1)],file=paste("remove_from",traits[i],sep="."),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

## fam[!missing,6]<-ann[posns[!missing],traits[i]]
the.fam<-paste(traits[i],"fam",sep=".")
### already wrote fam
## write.table(fam,file=the.fam,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)


ip<-1
projects
for(ip in 1:length(projects)){
print(projects[ip])
#setwd(projects[ip])

the.bed<-paste(projects[ip],"bed",sep=".")
the.bim<-paste(projects[ip],"bim",sep=".")

  system( paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam,"--remove",paste("remove_from",traits[i],sep="."),"--hardy","--allow-no-sex","--out",paste(traits[i],projects[ip],sep="."),"--noweb",sep=" ") )

} # ip loop over projects / chromosomes
} # i loop over traits

   







############old fro


















                                        #------------------------------------------------------------------
# loccation of rescored files  #1  
ROOT.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP.list.re-checked.ED.chr1-19" # location of recorces beds

##IMPORTANT  The bed/list/score file must contain the chromosome name prefix chromsome sufix
## names of the list/score/bed files must be the chromosome
### set up the list names
check.flag<-"ED.modified"
    
list.file.prefix<-"meta.snp.test_"
list.file.sufix<-".txt"
### set up the score file names    
score.file.prefix<-list.file.prefix
score.file.sufix<-paste(list.file.sufix,".scores",sep="")
# set up to get the bed.file.names
bed.file.prefix<-"meta.snp.test_"
bed.file.sufix<-".bed"    
#---------------------------------------------------------------------------------    

    

############ location of fam and bim files to use GENERIC
the.bim.root<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Cluster_vis_clean/Cluster_viz/" # .chr.bim added later
the.fam<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Cluster_vis_clean/Cluster_viz/chr1/AOGC_gentrain.fam"

 fam.template.file<-"all.lists_filtered-clean.f.fam" ## has been flipped to match
    

    






##  NEED TO GET BEFORE FILTERED
##
## data<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/final.meta.analysis.single.point.Filtered2.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)











  ############ insert into previos data chip P, meta P x2 hwe for chip


########################## get processed bim files that are all on the forward strand in annovar format"
code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")

######


annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
  
traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")

names(traits)<-rep("logistic",times=length(traits))
names(traits)[1:3]<-rep("assoc",times=3)
traits
traits %in% colnames(ann)


#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
names(traits)<-rep(".assoc.logistic",times=length(traits))
names(traits)[1:3]<-rep(".qassoc",times=3)
traits


## fam.template.file
## # fam.template.file set above in cases where flipping data is required
## fam<-read.table(fam.template.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
## projects<-gsub(".fam$","",fam.template.file)


    
the.seq.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/AOGC_vcf_to_plink" # dirname(fam.template.file)
files<-dir(the.seq.dir)  
projects<-gsub(".bed","",files[grepl(".bed$",files)])  # files[grepl(".bed$",files)]


project.prefix<-"AOGC_to_vcf_filtered_" ## outpt used same name as above
project.prefix.old<-"AOGC_vcf_final_chr"

projects.old<-projects[grepl(paste("^",project.prefix.old,sep=""),projects)]
projects<-projects[grepl(paste("^",project.prefix,sep=""),projects)] #

projects
projects.old

the.chr.old<-gsub(project.prefix.old,"",projects.old)
the.chr<-gsub(project.prefix,"",projects)

sum(the.chr.old !=the.chr)


# Corresdonding in same order for exome chip
 logistic.project.dir<-ROOT.dir
 logistic.project<-gsub(".fam$","",fam.template.file)
 if(!grepl("/$",logistic.project.dir)){logistic.project.dir<-paste(logistic.project.dir,"/",sep="")}

#logistic.files<-paste(logistic.project.dir,traits,".",logistic.project,names(traits),sep="")
hwe.files.chip<-paste(logistic.project.dir,traits,".",logistic.project,".hwe",sep="")
logistic.files<-paste(logistic.project.dir,"chip.common.extra.",traits,sep="")
logistic.seq.files<-paste(logistic.project.dir,"seq.common.extra.",traits,sep="")

#"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/filtered/"
meta.dir<-logistic.project.dir
meta.stderr.files<-paste(meta.dir,"META.STDERR.common.extra.",traits,"1.tbl",sep="")
meta.size.files<-paste(meta.dir,"META.SAMPLE.SIZE.common.extra.",traits,"1.tbl",sep="")

#data<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/final.meta.analysis.single.point.significant.reorg.evoker.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

#write.table(all.sig.final,file="final.meta.analysis.single.point.Filtered2.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
getwd()
setwd(logistic.project.dir)



ann.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/annoatation_short_exome_chip.txt"
ann.short.ori<-read.delim(ann.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)






############
########### want some extra info in data: RS.id, p.hardy.chip Stderr.chip  Stderr.seq 

##   colnames(data)[2]
## colnames(data)[colnames(data)=="traits.i."]<-"traits"
## data[1:20,1:10]
  
  i<-1

 all.meta<-{}
for (i in 1:length(traits)){


print(logistic.files[i])
chip<-read.table(logistic.files[i],header=T,fill=TRUE,stringsAsFactors=FALSE)
seq<-read.table(logistic.seq.files[i],header=T,fill=TRUE,stringsAsFactors=FALSE)
hwe.chip<-read.table(hwe.files.chip[i],header=T,fill=TRUE,stringsAsFactors=FALSE)
stderr<-read.table(meta.stderr.files[i],header=T,fill=TRUE,stringsAsFactors=FALSE)
sample.size<-read.table(meta.size.files[i],header=T,fill=TRUE,stringsAsFactors=FALSE)

ip<-1
for(ip in 1:length(projects)){
setwd(the.seq.dir)
############ read the hwe data
print(paste((traits)[i],projects[ip],"hwe",sep="."))

a.hwe<- read.table(paste((traits)[i],projects[ip],"hwe",sep="."),header=T,fill=TRUE,stringsAsFactors=FALSE)
a.hwe.old<- read.table(paste((traits)[i],projects.old[ip],"hwe",sep="."),header=T,fill=TRUE,stringsAsFactors=FALSE)

a.hwe<-as.matrix(a.hwe) ## colClasses="character"
a.hwe.old<-as.matrix(a.hwe.old)  
##############

if( sum( c("AFF","UNAFF") %in% a.hwe[,"TEST"] )==2 ){ # case/control unwind AFF and UNAFF
  aff<-a.hwe[,"TEST"]=="AFF"
  unaff<-a.hwe[,"TEST"]=="UNAFF"
  GENO.aff<-a.hwe[aff,"GENO"]
  GENO.unaff<-a.hwe[unaff,"GENO"]
 sum(aff)
 sum(unaff)
 if(sum(a.hwe[unaff,"SNP"]!=a.hwe[unaff,"SNP"])>0){print("Error in hwe alignments")}
  a.hwe<-cbind(a.hwe[unaff,c("CHR","SNP","A1","A2","P")],GENO.aff,GENO.unaff,stringsAsFactors=FALSE)
  aff<-a.hwe.old[,"TEST"]=="AFF"
  unaff<-a.hwe.old[,"TEST"]=="UNAFF"
  GENO.aff<-a.hwe.old[aff,"GENO"]
  GENO.unaff<-a.hwe.old[unaff,"GENO"]
 sum(aff)
 sum(unaff)
 if(sum(a.hwe.old[unaff,"SNP"]!=a.hwe.old[unaff,"SNP"])>0){print("Error in hwe alignments")}
  a.hwe.old<-cbind(a.hwe.old[unaff,c("CHR","SNP","A1","A2","P")],GENO.aff,GENO.unaff,stringsAsFactors=FALSE)
}

if(ip==1){hwe<-a.hwe}else{hwe<-rbind(hwe,a.hwe)}
if(ip==1){hwe.old<-a.hwe.old}else{hwe.old<-rbind(hwe.old,a.hwe.old)}

#dim(a.logistic.seq)
} ## loop over ip


if( sum( c("AFF","UNAFF") %in% hwe.chip[,"TEST"] )==2 ){ # case/control unwind AFF and UNAFF
  aff<-hwe.chip[,"TEST"]=="AFF"
  unaff<-hwe.chip[,"TEST"]=="UNAFF"
  GENO.aff<-hwe.chip[aff,"GENO"]
  GENO.unaff<-hwe.chip[unaff,"GENO"]
 sum(aff)
 sum(unaff)
 if(sum(hwe.chip[unaff,"SNP"]!=hwe.chip[unaff,"SNP"])>0){print("Error in hwe alignments")}
  hwe.chip<-cbind(hwe.chip[unaff,c("CHR","SNP","A1","A2","P")],GENO.aff,GENO.unaff,stringsAsFactors=FALSE)
}


setwd(logistic.project.dir)



chip[1:5,]
seq[1:5,]
#seq.unfilt[1:5,]
sample.size[1:5,]
stderr[1:5,]
hwe.chip[1:5,]
hwe[1:5,]
hwe.old[1:5,]

core.ann<-c("chr","start","REF","ALT")
seq.key<-build.key(seq,core.ann)


############## Sort out hwe for sequencing
posns<-match(hwe.old[,"SNP"],hwe[,"SNP"])
hwe<-hwe[posns,]

grep("^UQDI",hwe.old[,"SNP"])

hwe.old.key<-strsplit(hwe.old[,"SNP"],split=":")
start<-unlist(lapply(hwe.old.key,function(x) x[2]))
start[1:10]
hwe.old<-cbind(hwe.old,start)

hwe.key<-build.key(hwe.old,c("CHR","start","A2","A1"))

hwe.key[1:10]
posns<-match(seq.key,hwe.key)
missing<-is.na(posns)
hwe<-hwe[posns,]
hwe.old<-hwe.old[posns,]
hwe.old[,"SNP"]<-seq[,"SNP"]
hwe[,"SNP"]<-seq[,"SNP"]
###############################

dim(chip)
dim(seq)
dim(stderr)
## dim(logistic.seq)


## logistic.seq[1:5,]
## dim(logistic.seq)
## dim(seq)

######### alignen all to chip

posns<-match(chip[,"SNP"],seq[,"SNP"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
seq<-seq[posns,]

posns<-match(chip[,"SNP"],ann.short.ori[,"SNP"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
ann.short<-ann.short.ori[posns,]

posns<-match(chip[,"SNP"],stderr[,"MarkerName"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
stderr<-stderr[posns,]

posns<-match(chip[,"SNP"],sample.size[,"MarkerName"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
sample.size<-sample.size[posns,]

posns<-match(chip[,"SNP"],hwe[,"SNP"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
hwe<-hwe[posns,]

posns<-match(chip[,"SNP"],hwe.old[,"SNP"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
hwe.old<-hwe.old[posns,]

posns<-match(chip[,"SNP"],hwe.chip[,"SNP"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
hwe.chip<-hwe.chip[posns,]


ann.short[1:5,]


dim(seq)
dim(chip)
dim(stderr)
dim(sample.size)
dim(ann.short)
dim(hwe)
dim(hwe.old)
dim(hwe.chip)


P.value.StdErr<-stderr[,"P.value"]
P.value.Size<-sample.size[,"P.value"]
Weight<-sample.size[,"Weight"]

P.Seq.Filt<-seq[,"P"]
P.Seq<-rep(NA,times=length(P.Seq.Filt))
P.Chip<-chip[,"P"]

BETA.chip<-chip[,"BETA"]
BETA.seq<-seq[,"BETA"]

SE.chip<-chip[,"SE"]
SE.seq<-seq[,"SE"]

#### hwe by trait and pheno so can modify hwe.ori<-hwe ; hwe.old.ori<-hwe.old # hwe<-hwe.ori ; hwe.old<-hwe.old.ori
hwe<-hwe[,(grepl("^GENO",colnames(hwe)) |  colnames(hwe)=="P")]
hwe.old<-hwe.old[,(grepl("^GENO",colnames(hwe.old)) |  colnames(hwe.old)=="P")]
hwe.chip<-hwe.chip[,(grepl("^GENO",colnames(hwe.chip)) |  colnames(hwe.chip)=="P")]

colnames(hwe)[colnames(hwe)=="P"]<-"P.hardy"
colnames(hwe.old)[colnames(hwe.old)=="P"]<-"P.hardy"
colnames(hwe.chip)[colnames(hwe.chip)=="P"]<-"P.hardy"

colnames(hwe)<-paste(colnames(hwe),"Seq.Filt",sep=".")
colnames(hwe.old)<-paste(colnames(hwe.old),"Seq",sep=".")
colnames(hwe.chip)<-paste(colnames(hwe.chip),"Chip",sep=".")

all.hwe<-cbind(hwe,hwe.old,hwe.chip)

all.hwe[1:5,]
all.hwe<-as.matrix(all.hwe)

extra.cols<-c("GENO.Chip","GENO.Seq.Filt","GENO.Seq","GENO.aff.Chip","GENO.aff.Seq.Filt","GENO.unaff.Chip","GENO.unaff.Seq.Filt","GENO.aff.Seq","GENO.unaff.Seq","P.hardy.Chip","P.hardy.Seq.Filt","P.hardy.Seq")
extra<-matrix(data=NA,nrow=dim(all.hwe)[1],ncol=length(extra.cols))
colnames(extra)<-extra.cols

extra[,extra.cols[extra.cols %in% colnames(all.hwe)]]<-all.hwe[,extra.cols[extra.cols %in% colnames(all.hwe)]]
                                        # reorganise HWE
extra[1:5,]

dim(extra)
length(BETA.chip)


wanted.cols<-c("rs.id","design","refGene..type","Gene.Names","description","gerp.scores","PolyPhen.desc","PolyPhen.scores","SIFT.desc","SIFT.scores","skeletome","mouse.defect","sewell.cycling","Dequeant.cycling","ingenuity.bone.genes","Consequence.Embl","Amino_acids.Embl","MAF.lt.0.001","MAF.lt.0.5","Final_Score")


#wanted<-c("CHR","SNP","POS","A1","TEST","NMISS","OR","SE","L95","U95","STAT","P","design")

## cols.to.halpview.fix<-c("GENO","design","refGene..type","Gene.Names","description","gerp.scores","PolyPhen.desc","PolyPhen.scores","SIFT.desc","SIFT.scores","skeletome","mouse.defect","sewell.cycling","Dequeant.cycling","ingenuity.bone.genes","Consequence.Embl","Amino_acids.Embl","MAF","MAF.lt.0.001","MAF.lt.0.5")

wanted.cols %in% colnames(ann.short)

meta<-cbind(traits[i],stderr[,c("MarkerName","Allele1","Allele2","Effect","StdErr","Direction")],P.value.StdErr,P.value.Size,Weight,P.Chip,P.Seq.Filt,P.Seq,extra,ann.short[,wanted.cols],chip[,c("chr","start","REF","ALT","SNP")],BETA.chip,BETA.seq,SE.chip,SE.seq,stringsAsFactors=FALSE)

dim(meta)
colnames(meta)

if(is.null(dim(all.meta))){
  all.meta<-meta
}else{
  all.meta<-rbind(all.meta,meta)
}


} ## loop over traits











####################old from line 3380

#########################  ANNOTATE META ANALYSIS/media/UQCCG-Analysis
#########################  ANNOTATE META ANALYSIS/media/UQCCG-Analysis



#traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")

code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")

traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")
#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
names(traits)<-rep(".assoc.logistic",times=length(traits))
names(traits)[1:3]<-rep(".qassoc",times=3)
traits




hwe.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/zCall_AOGC.with.evoker_corrected_clean_FINAL.hwe"
hwe<-read.table(hwe.file,header=T,fill=TRUE,stringsAsFactors=FALSE)
wanted<-grepl("ALL",hwe[,"TEST"])
hwe.ori<-hwe[wanted,]
hwe.ori[1:5,]

frq.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/zCall_AOGC.with.evoker_corrected_clean_FINAL.frq"
frq.ori<-read.table(frq.file,header=T,fill=TRUE,stringsAsFactors=FALSE)
frq.ori[1:5,]

ann.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/annoatation_short_exome_chip.txt"
ann.short<-read.delim(ann.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

posns<-match(ann.short[,"SNP"],hwe.ori[,"SNP"])
missing<-is.na(posns)
sum(missing)
hwe.ori<-hwe.ori[posns,]

posns<-match(ann.short[,"SNP"],frq.ori[,"SNP"])
missing<-is.na(posns)
sum(missing)
frq.ori<-frq.ori[posns,]

ann.short.ori<-cbind(ann.short,hwe.ori,frq.ori,stringsAsFactors=FALSE)
ann.short[1:5,]


the.seq.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/AOGC_vcf_to_plink" # dirname(fam.template.file)
files<-dir(the.seq.dir)  
projects<-gsub(".bed","",files[grepl(".bed$",files)])  # files[grepl(".bed$",files)]
project.prefix.old<-"AOGC_vcf_final_chr"
project.prefix<-"AOGC_to_vcf_filtered_"
projects.old<-projects[grepl(paste("^",project.prefix.old,sep=""),projects)]
projects<-projects[grepl(paste("^",project.prefix,sep=""),projects)] #

projects
projects.old

the.chr.old<-gsub(project.prefix.old,"",projects.old)
the.chr<-gsub(project.prefix,"",projects)

sum(the.chr.old !=the.chr)
the.chr<-the.chr[!(the.chr %in% c(23,24))]


### to use this would need to add reat the the single point sequening results to get the key right from the SNP name and then align those.
## ref.cohort<-paste("/media/UQCCG/Sequencing/Projects/AOGC-NGS/Australian Reference Genome Cohort/AOGC_Controls/AOGC-Genotyping.output.chr",the.chr,".AOGC_ALL.controls.txt",sep="")

## ref.cohort
## ip<-1
## for(ip in 1:length(ref.cohort)){
##   print(ref.cohort[ip])
## a.ref<- read.table(ref.cohort[ip],header=T,fill=TRUE,stringsAsFactors=FALSE,sep="\t",quote="")
## if(ip==1){ref.ann<-a.ref}else{ref.ann<-rbind(ref.ann,a.ref)}
## }
  
## ref.ann.key<-build.key(ref.ann,c("chr","start","end","REF","ALT","TYPE"))
## ref.ann.key[1:20]

########### need to load in seq from one of the runs below
## seq.key<-seq[,c("chr","start","start","REF","ALT","SNP")]
## colnames(seq.key)<-c("chr","start","end","REF","ALT","SNP")

## seq.key[,c("chr","start","end","REF","ALT")]<-redo.start.end.annovar(seq.key[,c("chr","start","end","REF","ALT")]) ## put in annovar format
## seq.key<-get.forward.strand.allele(seq.key)
## dim(seq.key)
## seq.key[1:5,]
## seq.key[,"chr"]<-gsub("^chr","",seq.key[,"chr"])
## seq.key.use<-build.key(seq.key,c("chr","start","end","REF","ALT"))
## ref.ann.key.use<-build.key(ref.ann,c("chr","start","end","REF","ALT"))

## seq.key.use[1:5]
## ref.ann.key.use[1:5]
## posns<-match(seq.key.use,ref.ann.key.use)
## missing<-is.na(posns)
## sum(missing)
## save(list=c("ref.ann.key","ref.ann.key.use","ref.ann","seq.key","seq.key.use"),file="/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/ref.ann.RData")

setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point")

filtered.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/filtered"
unfiltered.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/notFiltered"

#load("/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward.RData")
#load("/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_fromVCF.RData")
## load("/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_A1_A2.RData") # contains "bim.chip" "bim.seq"

## bim.seq[1:5,]
## bim.chip[1:5,]
## chipkey<-build.key(bim.chip,c("chr","start","REF","ALT"))
## seqkey<-build.key(bim.seq,c("chr","start","REF","ALT"))
## posns<-match(chipkey,seqkey)
## missing<-is.na(posns)
## sum(!missing)
## chipSeqkey<-cbind(bim.chip[!missing,],bim.seq[posns[!missing],"SNP"])
## chipSeqkey[1:5,]
## colnames(chipSeqkey)[7]<-"marker"
## save(list=c("bim.chip","bim.seq","chipSeqkey"),file="/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_A1_A2_withKEY.RData")

load("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_A1_A2_withKEY.RData") # contains "bim.chip" "bim.seq" "chipSeqkey" key using posn and ref/alt

load("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/ref.ann.RData")
bim.chip.ori<-bim.chip
bim.seq.ori<-bim.seq
ref.ann.key.ori<-ref.ann.key
                                 
i<-1
for (i in 1:length(traits)){

  ## for (i in 12:17){
ann.short<-ann.short.ori
ref.ann.key<-ref.ann.key.ori
bim.seq<-bim.seq.ori
setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/AOGC_vcf_to_plink")
ip<-1

for(ip in 1:length(projects)){


if(i %in% 1:3){ # use rescalled data
#  NewBETA
logistic.file.seq<-paste("NewBETA",traits[i],projects[ip],sep=".")
print(logistic.file.seq)
setwd(the.seq.dir)
a.logistic.seq<-read.table(logistic.file.seq,header=T,fill=TRUE,stringsAsFactors=FALSE)
}else{
logistic.file.seq<-paste(traits[i],projects[ip],sep=".")
print(logistic.file.seq)
setwd(the.seq.dir)
a.logistic.seq<-read.table(paste(logistic.file.seq,names(traits)[i],sep=""),header=T,fill=TRUE,stringsAsFactors=FALSE)
}

colnames(a.logistic.seq)[colnames(a.logistic.seq)=="BP"]<-"POS"

if("TEST" %in% colnames(a.logistic.seq)){
wanted<-grepl("ADD",a.logistic.seq[,"TEST"])
a.logistic.seq<-a.logistic.seq[wanted,]
}

if(ip==1){logistic.seq<-a.logistic.seq}else{logistic.seq<-rbind(logistic.seq,a.logistic.seq)}

############ read the hwe data
print(paste((traits)[i],projects[ip],"hwe",sep="."))

a.hwe<- read.table(paste((traits)[i],projects[ip],"hwe",sep="."),header=T,fill=TRUE,stringsAsFactors=FALSE)
a.hwe.old<- read.table(paste((traits)[i],projects.old[ip],"hwe",sep="."),header=T,fill=TRUE,stringsAsFactors=FALSE)

a.hwe<-as.matrix(a.hwe) ## colClasses="character"
a.hwe.old<-as.matrix(a.hwe.old)  
##############

if( sum( c("AFF","UNAFF") %in% a.hwe[,"TEST"] )==2 ){ # case/control unwind AFF and UNAFF
  aff<-a.hwe[,"TEST"]=="AFF"
  unaff<-a.hwe[,"TEST"]=="UNAFF"
  GENO.aff<-a.hwe[aff,"GENO"]
  GENO.unaff<-a.hwe[unaff,"GENO"]
 sum(aff)
 sum(unaff)
 if(sum(a.hwe[unaff,"SNP"]!=a.hwe[unaff,"SNP"])>0){print("Error in hwe alignments")}
  a.hwe<-cbind(a.hwe[unaff,c("CHR","SNP","A1","A2","P")],GENO.aff,GENO.unaff,stringsAsFactors=FALSE)
  aff<-a.hwe.old[,"TEST"]=="AFF"
  unaff<-a.hwe.old[,"TEST"]=="UNAFF"
  GENO.aff<-a.hwe.old[aff,"GENO"]
  GENO.unaff<-a.hwe.old[unaff,"GENO"]
 sum(aff)
 sum(unaff)
 if(sum(a.hwe.old[unaff,"SNP"]!=a.hwe.old[unaff,"SNP"])>0){print("Error in hwe alignments")}
  a.hwe.old<-cbind(a.hwe.old[unaff,c("CHR","SNP","A1","A2","P")],GENO.aff,GENO.unaff,stringsAsFactors=FALSE)
}




if(ip==1){hwe<-a.hwe}else{hwe<-rbind(hwe,a.hwe)}
if(ip==1){hwe.old<-a.hwe.old}else{hwe.old<-rbind(hwe.old,a.hwe.old)}

#dim(a.logistic.seq)
} ## loop over ip






dim(logistic.seq)
dim(bim.seq)
#colnames(logistic.seq)[colnames(logistic.seq)=="BP"]<-"POS"
#gefos[,1] %in%  data[,"MarkerName"]
## logistic.seq[1:5,]
## logistic[1:5,]
## bim.seq[1:5,]

## posns<-match(logistic.seq[,"SNP"],bim.seq[,"SNP"]) # for annotation
## missing<-is.na(posns)
## sum(missing)
## logistic.seq<-logistic.seq[!missing,]
## bim.seq<-bim.seq[posns[!missing],]
###################################################



dim(hwe)
dim(hwe.old)
hwe[1:5,]
hwe.old[1:5,]


out.file<-paste("META.BOTH.common.RESCALED.",traits[i],"1.tbl",sep="")
out.file.sig<-paste("META.BOTH.common.RESCALED.SIGNIF2.",traits[i],"1.tbl",sep="")
out.file.ann<-paste("META.BOTH.common.RESCALED.",traits[i],"1.tbl.ann",sep="")


  
setwd(filtered.dir)
## sample.size.file<-paste("META.SAMPLE.SIZE.common.",traits[i],"1.tbl",sep="")
## stderr.file<-paste("META.STDERR.common.",traits[i],"1.tbl",sep="")

sample.size.file<-paste("META.SAMPLE.SIZE.",traits[i],"1.tbl",sep="")
stderr.file<-paste("META.STDERR.",traits[i],"1.tbl",sep="")

chip<-read.delim(paste("chip",traits[i],sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
seq<-read.delim(paste("seq",traits[i],sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

sample.size<-read.delim(sample.size.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
stderr<-read.delim(stderr.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)


setwd(unfiltered.dir)
### chip.unfilt<-read.delim(paste("chip",traits[i],sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) # chip  should be the same
seq.unfilt<-read.delim(paste("seq",traits[i],sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)




chip[1:5,]
seq[1:5,]
seq.unfilt[1:5,]
sample.size[1:5,]
stderr[1:5,]
hwe[1:5,]
hwe.old[1:5,]


dim(seq)
dim(stderr)
dim(logistic.seq)
length(ref.ann.key.ori)

logistic.seq[1:5,]
dim(logistic.seq)
dim(seq)

posns<-match(chip[,"SNP"],seq[,"SNP"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
common<-chip[!missing,"SNP"]

gefos[,1] %in%  chip[,"SNP"]
gefos[,1] %in%  common
#Collect stuff that is
## 1)significant ORIGINALLY or in meta analysis
## 2)in common no matter what

## posns<-match(chip[,"SNP"],seq[,"SNP"]) ## make chip the default id name before.
## missing<-is.na(posns)
## sum(!missing)
## common<-chip[!missing,"SNP"]

## hist(as.numeric(seq[,"NMISS"]))
## hist(as.numeric(chip[,"NMISS"]))

num.samples<-max(as.numeric(seq[,"NMISS"]),na.rm=TRUE)
min.miss<-num.samples-0.2*num.samples
min(as.numeric(seq[,"NMISS"],na.rm=TRUE))
sig.seq<-as.numeric(as.character(seq[,"P"]))< 0.00001 & !is.na(as.numeric(as.character(seq[,"P"]))) & as.numeric(seq[,"NMISS"]) > min.miss
sum(sig.seq)
sig.seq<-seq[sig.seq,"SNP"]

num.samples<-max(as.numeric(chip[,"NMISS"]),na.rm=TRUE)
min.miss<-num.samples-0.2*num.samples
min(as.numeric(chip[,"NMISS"],na.rm=TRUE))
sig.chip<-as.numeric(chip[,"P"])< 0.0001 & !is.na(as.numeric(chip[,"P"])) & as.numeric(chip[,"NMISS"]) > min.miss
sum(sig.chip)
sig.chip<-chip[sig.chip,"SNP"]


to.check<-unique(c(common,sig.seq,sig.chip))

############ put everything in the order of to.check


## seq<-seq[posns,]
## BETA.chip.seq<-paste(signif(chip[,"BETA"],4),signif(seq[,"BETA"],4),sep=":")
## P.chip.seq<-paste(signif(chip[,"P"],4),signif(seq[,"P"],4),sep=":")

##  chip.seq<-cbind(chip[,"SNP"],BETA.chip.seq,P.chip.seq)
## colnames(chip.seq)[1]<-"SNP"

## seq[1:5,]
## chip[1:5,]
## sample.size[1:5,]
## stderr[1:5,]
## ann.short[1:5,]
## ann.short[1:5,]
## hwe[1:5,]
#hwe.ori<-hwe

## grep("UQDIexomechipWithAOGC_chr7:141619286:141619286:T:G-1_B_F_2098213710",stderr[,"MarkerName"]) #yes
## grep("UQDIexomechipWithAOGC_chr7:141619286:141619286:T:G-1_B_F_2098213710",seq[,"SNP"]) #yes
## grep("UQDIexomechipWithAOGC_chr7:141619286:141619286:T:G-1_B_F_2098213710",seq.unfilt[,"SNP"])
## seq.unfilt[,"SNP"]
## grep("UQDIexomechipWithAOGC_chr7:141619286:141619286:T:G-1_B_F_2098213710",bim.seq[,"SNP"])
## grep("UQDIexomechipWithAOGC_chr7:141619286:141619286:T:G-1_B_F_2098213710",seq.key[,"SNP"]) #yes
## grep("UQDIexomechipWithAOGC_chr7:141619286:141619286:T:G-1_B_F_2098213710",seq.key.use) #no
## grep("UQDIexomechipWithAOGC_chr7:141619286:141619286:T:G-1_B_F_2098213710",ref.ann.key.use) #no
## grep("UQDIexomechipWithAOGC_chr7:141619286:141619286:T:G-1_B_F_2098213710",hwe[,"SNP"]) #no
## grep("UQDIexomechipWithAOGC_chr7:141619286:141619286:T:G-1_B_F_2098213710",chipSeqkey[,"SNP"]) # yes WITH KEY
## chipSeqkey[39197,]
# chipSeqkey matched the chip to the sequencing ID


seq.key[1:5,]
ref.ann.key.use[1:5]
seq.key.use[1:5]
names(seq.key.use)[1:5]
names(ref.ann.key.use)[1:5]

posns<-match(chipSeqkey[,"marker"],hwe[,"SNP"])
missing<-is.na(posns)
hwe[posns[!missing],"SNP"]<-chipSeqkey[!missing,"SNP"]

posns<-match(chipSeqkey[,"marker"],hwe.old[,"SNP"])
missing<-is.na(posns)
hwe.old[posns[!missing],"SNP"]<-chipSeqkey[!missing,"SNP"]

## chipSeqkey[!missing,"SNP"][1:10]
## hwe[posns[!missing],"SNP"][1:10]

## grep("7:141619286",hwe[,"SNP"]) # 1494718
## grep("7:141619286",ref.ann.key.use)
## seq.key.use[grep("7:141619286",ref.ann.key.use)]
## seq.key[grep("7:141619286",ref.ann.key.use),]



posns<-match(stderr[,"MarkerName"],seq.key[,"SNP"])
sum(is.na(posns))
targets<-seq.key.use[posns]
posns<-match(targets,ref.ann.key.use)

## ref.ann.key.ori[1:50]
## ## get meta in same order
## posns<-match(seq[,"SNP"],ref.ann.key.ori)
## missing<-is.na(posns)
## sum(!missing)

## stderr[missing,"MarkerName"][1:20]



##     ###### check have the forward startd allele
  
##     print(dim(indels))
##  #   setwd(code.dir)
##     indels<-get.forward.strand.allele(indels)



                                             

posns<-match(stderr[,"MarkerName"],sample.size[,"MarkerName"])
missing<-is.na(posns)
sum(missing) # sample.size[missing,]
sample.size<-sample.size[posns,]

posns<-match(stderr[,"MarkerName"],chip[,"SNP"])
missing<-is.na(posns)
sum(!missing)
chip<-chip[posns,]

posns<-match(stderr[,"MarkerName"],seq[,"SNP"])
missing<-is.na(posns)
sum(missing)
seq<-seq[posns,]

posns<-match(stderr[,"MarkerName"],seq.unfilt[,"SNP"])
missing<-is.na(posns)
sum(missing)
#stderr[missing,][1:20,]
seq.unfilt<-seq.unfilt[posns,]

## posns<-match(stderr[,"MarkerName"],chip.unfilt[,"SNP"])
## missing<-is.na(posns)
## sum(missing)
## #stderr[missing,][1:20,]
## chip.unfilt<-chip.unfilt[posns,]

posns<-match(stderr[,"MarkerName"],ann.short[,"SNP"])
missing<-is.na(posns)
sum(!missing)
ann.short<-ann.short[posns,]

posns<-match(stderr[,"MarkerName"],hwe[,"SNP"])
missing<-is.na(posns)
sum(!missing)
hwe<-hwe[posns,]

posns<-match(stderr[,"MarkerName"],hwe.old[,"SNP"])
missing<-is.na(posns)
sum(!missing)
hwe.old<-hwe.old[posns,]

dim(seq)
dim(chip)
dim(stderr)
dim(sample.size)
dim(ann.short)
dim(hwe)
dim(hwe.old)

P.value.StdErr<-stderr[,"P.value"]
P.value.Size<-sample.size[,"P.value"]
Weight<-sample.size[,"Weight"]

P.Seq.Filt<-seq[,"P"]
P.Seq<-seq.unfilt[,"P"]
P.Chip<-chip[,"P"]


#### hwe by trait and pheno so can modify hwe.ori<-hwe ; hwe.old.ori<-hwe.old # hwe<-hwe.ori ; hwe.old<-hwe.old.ori
hwe<-hwe[,(grepl("^GENO",colnames(hwe)) |  colnames(hwe)=="P")]
hwe.old<-hwe.old[,(grepl("^GENO",colnames(hwe.old)) |  colnames(hwe.old)=="P")]

colnames(hwe)[colnames(hwe)=="P"]<-"P.hardy"
colnames(hwe.old)[colnames(hwe.old)=="P"]<-"P.hardy"

colnames(hwe)<-paste(colnames(hwe),"Seq.Filt",sep=".")
colnames(hwe.old)<-paste(colnames(hwe.old),"Seq",sep=".")


all.hwe<-cbind(hwe,hwe.old)

all.hwe[1:5,]


extra.cols<-c("GENO.Seq.Filt","GENO.Seq","GENO.aff.Seq.Filt","GENO.unaff.Seq.Filt","GENO.aff.Seq","GENO.unaff.Seq","P.hardy.Seq.Filt","P.hardy.Seq")
extra<-matrix(data=NA,nrow=dim(all.hwe)[1],ncol=length(extra.cols))
colnames(extra)<-extra.cols

extra[,extra.cols[extra.cols %in% colnames(all.hwe)]]<-all.hwe[,extra.cols[extra.cols %in% colnames(all.hwe)]]
                                        # reorganise HWE





wanted.cols<-c("GENO","design","refGene..type","Gene.Names","description","gerp.scores","PolyPhen.desc","PolyPhen.scores","SIFT.desc","SIFT.scores","skeletome","mouse.defect","sewell.cycling","Dequeant.cycling","ingenuity.bone.genes","Consequence.Embl","Amino_acids.Embl","MAF","MAF.lt.0.001","MAF.lt.0.5")


#wanted<-c("CHR","SNP","POS","A1","TEST","NMISS","OR","SE","L95","U95","STAT","P","design")

cols.to.halpview.fix<-c("GENO","design","refGene..type","Gene.Names","description","gerp.scores","PolyPhen.desc","PolyPhen.scores","SIFT.desc","SIFT.scores","skeletome","mouse.defect","sewell.cycling","Dequeant.cycling","ingenuity.bone.genes","Consequence.Embl","Amino_acids.Embl","MAF","MAF.lt.0.001","MAF.lt.0.5")



meta<-cbind(stderr[,c("MarkerName","Allele1","Allele2","Effect","StdErr","Direction")],P.value.StdErr,P.value.Size,Weight,P.Chip,P.Seq.Filt,P.Seq,extra,ann.short[,wanted.cols],seq[,c("chr","start","REF","ALT","SNP")],stringsAsFactors=FALSE)

write.table(meta,file=paste(out.file.sig,"ALL",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


posns<-match(to.check,meta[,"MarkerName"])
missing<-is.na(posns)
sum(missing)


meta<-meta[posns[!missing],]

setwd("/media/UQCCG-Analysis/AOGC_exome_chip/meta_analysis_single_point")
write.table(meta,file=out.file.sig,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
}



## hwe.old<-hwe.old[posns,]
## to.check
## meta[1:2,]











########################### aligned refann
sum(is.na(a.sig[,"snps"]))

seq.key[1:5,]

ref.ann[1:5,]
## seq.key<-build.key(ref.ann,c("chr","start","end","REF","ALT","TYPE"))
## seq.key.s<-build.key(ref.ann,c("chr","start","REF","ALT"))
snps<-a.sig[,"snps"]


#### these two keys in the same order make seq.key.s the best key
seq.key.new<-build.key(ref.ann,c("chr","start","end","ALT","REF","TYPE"))
seq.key.s<-build.key(ref.ann,c("chr","start","ALT","REF"))
sum(!(seq.key[,"SNP"] %in% seq.key.s))

seq.key.new<-build.key(ref.ann,c("chr","start","end","ALT","REF","TYPE"))



posns<-match(snps,seq.key.new)
missing<-is.na(posns)
sum(!missing)
sum(missing)
### thes ones will match snps
## snps[!missing][1:10]
## seq.key.new[posns[!missing]][1:10]

seq.key.s[posns[!missing]]<-snps[missing]

                                         
a.sig[missing,][1:5,1:10]
types<-table(a.sig[missing,"design"])
types
sum(types) ### most are  chip ids
sum(missing)
sum(is.na(a.sig[missing,"chr"])) 


######## build new key for exome chip ids
snps.s<-build.key(a.sig[missing,],c("chr","start","ALT","REF"))
a.sig[1:5,]
snps.s[1:5]
seq.key.s[1:5]

posns.s<-match(snps.s,seq.key.s)
missing.s<-is.na(posns.s)
sum(!missing.s)
snps.s[!missing.s][1:10]
seq.key.s[posns.s[!missing.s]][1:10]
snps[missing][!missing.s][1:10]
seq.key.s[posns.s[!missing.s]]<-snps[missing][!missing.s]


posns<-match(snps,seq.key.s)
missing<-is.na(posns)
sum(missing)
############################ try flipping remaining

a.sig[missing,][1:5,]
a.sig[missing,][1:5,1:10]
types<-table(a.sig[missing,"design"])
types
sum(types) ### most are  chip ids
sum(missing)
sum(is.na(a.sig[missing,"chr"]))

snps.s<-build.key(a.sig[missing,],c("chr","start","REF","ALT"))
a.sig[1:5,]
snps.s[1:5]
seq.key.s[1:5]

posns.s<-match(snps.s,seq.key.s)
missing.s<-is.na(posns.s)
sum(!missing.s)
snps.s[!missing.s][1:10]
seq.key.s[posns.s[!missing.s]][1:10]
snps[missing][!missing.s][1:10]
seq.key.s[posns.s[!missing.s]]<-snps[missing][!missing.s]

posns<-match(snps,seq.key.s)
missing<-is.na(posns)
sum(missing) #23205



## extra<-ref.ann[posns,c( "chr","start","end","Hetero.ALT.reads","Hetero.REF.reads","Hetero.Read.Balance")]
## extra<-ref.ann[posns,c( "Hetero.ALT.reads","Hetero.REF.reads","Hetero.Read.Balance","Gene.Names","description","Amino_acids.Embl","FILTER","gerp.scores")]
extra<-ref.ann[posns,c("Hetero.ALT.reads","Hetero.REF.reads","Hetero.Read.Balance","FILTER","gerp.scores")]
 ## if(sum(c("GENO.aff.Seq.Filt","GENO.unaff.Seq.Filt","GENO.aff.Seq","GENO.unaff.Seq") %in% colnames(a.sig))==4){
 ##   a.sig<-


HI Matt,



1.       AOGC final SNP selection

Been going flat out on this, Ive reformatted the output to make the process a bit easier and Emma  I can finalize these in the next few days.

Also need to send the genewise data for the whole genome extremes

2.       Cervical cancer outstanding stuff

Be another week 

3.       HLA-B imputation from AOGC?

This has finished, and I've asked Felicity finish up the QC script only takes an hour or two. The output is dosage and plink (with best guess genotypes). What format do you think they can work with? 


HI Adrian,

Some of Matt's collaborators are worth with AOGC and mostly the Australia's AS samples.

Can you point me to the HLA imputation of AS immunochip? of AS samples in general

Dosage (and plink if you have it ) with a lits of the outliers.

I can pull what I need form there, well full dataset in any case. I'll put into BC.

Thanks
Paul

PS I moved your  Chinese TB folder to /UQCCG/GWAS  (it was a little too bing to copy) lefet a README in the original folder in case you forget
