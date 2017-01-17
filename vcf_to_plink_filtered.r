





options(width=250,max.print=5000)

code.dir<-"/media/Bioinform-D/Research/AML sequencing"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")


annotation.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)


#traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
## traits<-c("TOT_HIP_GCM","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")

## traits<-c("BMD_EFF_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")


traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")

names(traits)<-rep("logistic",times=length(traits))
names(traits)[1:3]<-rep("assoc",times=3)
traits
traits %in% colnames(ann)




#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")


input.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2012-06-18_All_AOGC/Analysis"
files<-dir(input.dir)  
project.files<-files[ grepl("^AOGC-Genotyping.output.",files) & grepl(".AOGC_ALL.analysis.txt$",files)]
project.files
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2012-06-18_All_AOGC/Analysis/AOGC-Genotyping.output.chr8.AOGC_ALL.analysis.txt

ichr<-1
i<-1

library(doMC)
num.cores<-7
registerDoMC(cores=num.cores)


for(ichr in 2:length(project.files)){
################## fast read ###########
  setwd(input.dir)
column.labels<-read.delim(project.files[ichr],header=F,nrows=1,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
num.vars<-dim(column.labels)[2]
a.indel<-scan(project.files[ichr],what=character(num.vars),skip=1,sep="\t",fill=TRUE,na.strings="",quote="\"")
num.lines<-length(a.indel)/(num.vars)
dim(a.indel)<-c(num.vars,num.lines)
a.indel<-t(a.indel)
colnames(a.indel)<-column.labels
########################################
the.chr<-a.indel[1,"chr"]
print(paste("Doing Chromosome ",the.chr))

if(!grepl("^chr",the.chr)){
a.indel[,"chr"]<-paste("chr",a.indel[,"chr"],sep="")
}
core.ann<-c("chr","start","end","REF","ALT","TYPE")
key<-build.key(a.indel,core.ann)

rownames(a.indel)<-key
plink.chr<-the.chr


if(plink.chr=="X"){plink.chr<-23}
if(plink.chr=="Y"){plink.chr<-24}
if(plink.chr=="XY"){plink.chr<-25}
if(plink.chr=="M"){plink.chr<-26}
plink.chr

the.samples<-colnames(a.indel)[grepl(".GT$",colnames(a.indel))]
length(the.samples)
the.samples.fam<-gsub(".GT$","",the.samples)


####
a.fam<-cbind(the.samples.fam,the.samples.fam,0,0,9,-9)
colnames(a.fam)<-c("PATIENT","SAMPLE","mother","father","sex","pheno")


the.fam<-paste(traits[i],"fam",sep=".")
a.fam<-read.table(the.fam,header=F,fill=TRUE,stringsAsFactors=FALSE)
colnames(a.fam)<-c("PATIENT","SAMPLE","mother","father","sex","pheno")
### colnames(fam)<-c("PATIENT","SAMPLE","mother","father","sex","pheno")


if(dim(a.indel)[1]<200){
num.bits<-1; print(paste("Only",dim(a.indel)[1],"genotypes remaining",sep=" ")) }else{num.bits<-num.cores}

#chunksize=as.integer(dim(a.indel)[1]/num.bits) )

fil.genotypes<-foreach(a.indel.bit=iter(a.indel,by='row',chunksize=5000 ), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype(a.indel.bit,the.samples.fam,prefix="",suffix="",20,0.2,0.80,0.05,0.95,10,5)

dim(fil.genotypes)

a.indel<-cbind(a.indel[,c("chr","start","end","REF","ALT","TYPE")],fil.genotypes)

root<-paste("AOGC_to_vcf_filtered",plink.chr,sep="_")
write.plink(a.indel,root)

} # short_cut loop to write filtered genotypes


############ below for testing:

g.indel<-read.plink.raw("/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink/AOGC_to_vcf_filtered_chr10") #A1->ALT POS->start ALT is 1 REF 0
dim(g.indel)
g.indel[1:15,1:10]
rownames(g.indel)<-build.key(g.indel,c("chr","start","REF","ALT"))


/media/UQCCG/Sequencing/Projects/AOGC-NGS/Australian Reference Genome Cohort/AOGC_Controls.zip
