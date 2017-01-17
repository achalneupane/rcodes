######################################################### single point for exome chip
######################################################### single point for exome chip
######################################################### single point for exome chip
#fam.template.file<-"recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL"

fam.template.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.fam"
#fam.template.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL_chr1.fam"
#fam.template.file<-"/media/UQCCG/UQCCG-Projects/PAUL_LEO/AOGC exome chip core genotyping/keep.txt"
annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"

#recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL
#the.chr<-basename(fam.template.file)
## the.fam.dir<-dirname(fam.template.file)
## files<-dir(the.fam.dir)  
## projects<-gsub(".bed","",files[grepl(".bed$",files)])  # files[grepl(".bed$",files)]
setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/")
projects<-"recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL"

fam<-read.table(fam.template.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

related.to.remove<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.recode.txt",header=F,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
related.to.remove<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.txt",header=F,fill=TRUE,sep=",",stringsAsFactors=FALSE)
related.to.remove[1:5,]

sum(related.to.remove[,1] %in% fam[,1])


traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")

#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
names(traits)<-rep("logistic",times=length(traits))
names(traits)[1:3]<-rep("assoc",times=3)
traits
traits %in% colnames(ann)


fam[1:5,]

i<-1
the.fam.dir<-dirname(fam.template.file)
setwd(the.fam.dir)


for (i in 4:length(traits)){

  
fam[,6]<--9
posns<-match(fam[,1],ann[,"PATIENT"])
missing<-is.na(posns)
print(sum(missing))

no.pheno<-ann[posns,traits[i]]==3 | ann[posns,traits[i]]==-9 | is.na(ann[posns,traits[i]])  # ann[posns,traits[i]]

print(sum(no.pheno))


missing<-missing | no.pheno
#cbind(fam[no.pheno,],ann[posns[no.pheno],c("PATIENT",traits[i])])
print(sum(missing))

if(names(traits)[i]=="logistic"){
 using.covars<-c("AGE_SCAN","PCA1","PCA2","PCA3","PCA4")
  covars<- ann[posns,c("PATIENT","PATIENT", using.covars)]
 colnames(covars)[1:2]<-c("FID","IID")
 covars[1:5,]
  missing.covars<-is.na(covars[,using.covars])
   missing.covars<-apply(missing.covars,1,sum,na.rm=TRUE)
 missing.covars<-  missing.covars>0
missing<-missing | missing.covars
write.table(covars[!missing,],file=paste("covars",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
}


has.missing<-sum(missing)>0
has.missing
write.table(fam[missing,c(1,1)],file=paste("remove_from",traits[i],sep="."),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

fam[!missing,6]<-ann[posns[!missing],traits[i]]
the.fam<-paste(traits[i],"fam",sep=".")

write.table(fam,file=the.fam,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)



ip<-1
projects
for(ip in 1:length(projects)){
print(projects[ip])
#setwd(projects[ip])

the.bed<-paste(projects[ip],"bed",sep=".")
the.bim<-paste(projects[ip],"bim",sep=".")

#  system( paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam,"--remove",paste("remove_from",traits[i],sep="."),"--hardy","--allow-no-sex","--out",paste(traits[i],projects[ip],sep="."),"--noweb",sep=" ") )


if(names(traits)[i]=="assoc"){
  system( paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam,"--remove",paste("remove_from",traits[i],sep="."),"--assoc","--allow-no-sex","--out",paste(traits[i],projects[ip],sep="."),"--noweb",sep=" ") )
}

if(names(traits)[i]=="logistic"){
  system( paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam,"--remove",paste("remove_from",traits[i],sep="."),"--covar",paste("covars",traits[i],sep="."),"--covar-name AGE_SCAN,PCA1,PCA2,PCA3,PCA4","--logistic --ci 0.95 ","--allow-no-sex","--out",paste(traits[i],projects[ip],sep="."),"--noweb",sep=" ") )
}
} # ip loop over chromosomes
} # i loop over traits

