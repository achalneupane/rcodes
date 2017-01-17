################################ prepare for META analysis of extra
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
#### all these in

######


##############################these bim files have the minor allele as A1 



code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")






output.dir<-"/media/UQCCG/GWAS/Scleroderma/scl_2014/metal"


setwd(output.dir)

bim.file<-"/media/UQCCG/GWAS/Scleroderma/scl_2014/metal/scl_case_control_common_final_gender.bim"
bim<-read.table(bim.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
colnames(bim)<-c("chr","SNP","cm","start","ALT","REF")
bim.chip<-bim


bim.file<-"/media/UQCCG/GWAS/Scleroderma/scl_2014/metal/dbgap_aogc_10kmerge_imputed_omni.bim"
bim<-read.table(bim.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
colnames(bim)<-c("chr","SNP","cm","start","ALT","REF") # A1 is labelled ALT here
bim.seq<-bim

bim.chip[1:5,]
bim.seq[1:5,]
dim(bim.seq)
dim(bim.chip)


#################################################

logistic.seq<-read.table("/media/UQCCG/GWAS/Scleroderma/scl_2014/metal/dbgap_aogc_10kmerge_imputed_omni_assoc.assoc.logistic",header=T,fill=TRUE,stringsAsFactors=FALSE)
logistic.seq[1:5,]
colnames(logistic.seq)[colnames(logistic.seq)=="BP"]<-"POS"

if("TEST" %in% colnames(logistic.seq)){
wanted<-grepl("ADD",logistic.seq[,"TEST"])
logistic.seq<-logistic.seq[wanted,]
}
dim(logistic.seq)
dim(bim.seq)


posns<-match(logistic.seq[,"SNP"],bim.seq[,"SNP"])
missing<-is.na(posns)
sum(missing)
bim.seq<-bim.seq[posns,]

junk<-bim.seq[,"start"]==0 | bim.seq[,"chr"]==0 | is.na(bim.seq[,"chr"])
sum(junk)
      
bim.seq<-bim.seq[!junk,]
logistic.seq<-logistic.seq[!junk,]


      
key.logistic.seq<-build.key(bim.seq,core.ann)
dups<-duplicated(key.logistic.seq)
      sum(dups)

rownames(logistic.seq)<-key.logistic.seq
TYPE<-rep("SEQ",times=length(key.logistic.seq))
#logistic.seq<-cbind(logistic.seq,TYPE,bim.seq[, core.ann],stringsAsFactors=FALSE)
logistic.seq<-cbind(logistic.seq,TYPE,bim.seq,stringsAsFactors=FALSE)
logistic.seq[1:5,]

#############################################
 logistic<-read.table("/media/UQCCG/GWAS/Scleroderma/scl_2014/metal/scl_case_control_common_final_gender.assoc.logistic",header=T,fill=TRUE,stringsAsFactors=FALSE)
logistic[1:5,]
colnames(logistic)[colnames(logistic)=="BP"]<-"POS"

if("TEST" %in% colnames(logistic)){
wanted<-grepl("ADD",logistic[,"TEST"])
logistic<-logistic[wanted,]
}
dim(logistic)
dim(bim.chip)


posns<-match(logistic[,"SNP"],bim.chip[,"SNP"])
missing<-is.na(posns)
sum(missing)
bim.chip<-bim.chip[posns,]

junk<-bim.chip[,"start"]==0 | bim.chip[,"chr"]==0 | is.na(bim.chip[,"chr"])
sum(junk)

bim.chip[junk,][1:5,]
logistic[junk,][1:5,]
      
bim.chip<-bim.chip[!junk,]
logistic<-logistic[!junk,]


      
key.logistic<-build.key(bim.chip,core.ann)
dups<-duplicated(key.logistic)
      sum(dups)

rownames(logistic)<-key.logistic
TYPE<-rep("GENO",times=length(key.logistic))
#logistic<-cbind(logistic,TYPE,bim.chip[, core.ann],stringsAsFactors=FALSE)
logistic<-cbind(logistic,TYPE,bim.chip,stringsAsFactors=FALSE)

      #############################


dim(logistic.seq)
dim(logistic)


logistic.seq[1:5,]
logistic[1:5,]

#################### get orderso can compare efefct directions

## posns<-match(logistic.seq[,"SNP"],logistic[,"SNP"])
##    missing<-is.na(posns)
## sum(missing)

## ## logistic.miss<-logistic[posns[missing],]
## ## logistic<-logistic[posns[!missing],]
## ## logistic.seq.miss<-logistic.seq[missing,]
## ## logistic.seq<-logistic.seq[!missing,]      

## dim(logistic.seq)
## dim(logistic)

##  logistic.seq[1:5,]
## logistic[1:5,]     
## ##################
#########
if("A1" %in% colnames(logistic)){ ### make sure ALT is A1 so effect direction is consistant - dont need to change sign of dirn then


  ################ set ALT and A! for the same
ref.diff<-logistic[,"A1"]!=logistic[,"ALT"]
sum(ref.diff)
logistic[ref.diff,][1:10,]
the.ref<-logistic[ref.diff,"REF"]
logistic[ref.diff,"REF"]<-logistic[ref.diff,"ALT"]
logistic[ref.diff,"ALT"]<-the.ref
ref.diff<-logistic[,"A1"]!=logistic[,"ALT"]
sum(ref.diff)

ref.diff<-logistic.seq[,"A1"]!=logistic.seq[,"ALT"]
sum(ref.diff)
logistic.seq[ref.diff,][1:10,]
the.ref<-logistic.seq[ref.diff,"REF"]
logistic.seq[ref.diff,"REF"]<-logistic.seq[ref.diff,"ALT"]
logistic.seq[ref.diff,"ALT"]<-the.ref
ref.diff<-logistic.seq[,"A1"]!=logistic.seq[,"ALT"]
sum(ref.diff)
#########################

############## don't chane the effect direction just let it go... meatl can sort out
## sum(is.na(logistic.seq[,"A1"]))
## remove<-is.na(logistic.seq[,"A1"]) | is.na(logistic.seq[,"ALT"]) |  is.na(logistic.seq[,"REF"]) # crap results
## sum(remove)
## logistic.seq<-logistic.seq[!remove,]

## ref.diff<-logistic.seq[,"A1"]!=logistic[,"A1"] 
## sum(ref.diff)
## logistic.seq[ref.diff,][1:10,]
## the.ref<-logistic.seq[ref.diff,"REF"]
## logistic.seq[ref.diff,"REF"]<-logistic.seq[ref.diff,"ALT"]
## logistic.seq[ref.diff,"ALT"]<-the.ref
## ref.diff<-logistic.seq[,"A1"]!=logistic.seq[,"ALT"]
## sum(ref.diff)

# redo.start.end.annovar 

}

### using minor allele now so should be fineetal will sort out differences




dim(logistic.seq)
dim(logistic)

logistic.seq[1:5,]
logistic[1:5,]
################### match up snp ids:
## posns<-match(logistic.seq[,"SNP"],logistic[,"SNP"])
## missing<-is.na(posns)
## sum(missing)
##  test<-logistic[,"P"]< 1e-6
##       sum(test)
## logistic[test,]
      
 


if(!("BETA" %in% colnames(logistic))){
BETA=log(as.numeric(logistic[,"OR"]))
logistic<-cbind(logistic,BETA,stringsAsFactors=FALSE)

## exp(BETA +- 1.96SE)=95% CI
## so
## se<--(log(as.numeric(logistic[,"L95"]))-BETA)/1.96
## logistic<-cbind(logistic,BETA,se,stringsAsFactors=FALSE)
## and se = SE in plink so !! SE Is the standard error of BETA!! no need to mess with
}

if(!("BETA" %in% colnames(logistic.seq))){
BETA=log(as.numeric(logistic.seq[,"OR"]))
logistic.seq<-cbind(logistic.seq,BETA,stringsAsFactors=FALSE)

#logistic.seq[,"SE"]<-log(as.numeric(logistic.seq[,"SE"]))
}
print("PAST BETA")
logistic[1:5,]
logistic.seq[1:5,]

setwd(output.dir)

### below is commented as just doing common
## write.table(logistic,file=paste("chip",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(logistic.seq,file=paste("seq",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

## write.table(logistic[posns[!missing],],file=paste("chip.common.extra",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(logistic.seq[!missing,],file=paste("seq.common.extra",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

write.table(logistic,file=paste("chip.common.extra","logistic",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(logistic.seq,file=paste("seq.common.extra","logistic",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



 logistic.seq[1:5,]
logistic[1:5,]   
##print(doing merge)
## common.cols<-colnames(logistic)[colnames(logistic) %in% colnames(logistic.seq)]
## the.merge<-rbind(logistic[,common.cols],logistic.seq[,common.cols])
## the.merge[1:5,]
## order.by<-order(the.merge[,"CHR"],the.merge[,"start"])
## the.merge<-the.merge[order.by,]
## the.merge[1:5,]
## write.table(the.merge,file=paste("merge.extra",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
##

## system(paste("sed s/TRAIT/",traits[i],"/ run_config_for_P_and_BETA_TEMPLATE.txt > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
## system(paste("./metal ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))

system(paste("sed s/TRAIT/common.extra.","logistic","/ run_config_for_P_and_BETA_TEMPLATE.txt > ",paste("sample.size.CONFIG.common.extra","logistic","txt",sep="."),sep=""))
system(paste("./metal ",paste("sample.size.CONFIG.common.extra","logistic","txt",sep="."),sep=""))


## system(paste("sed s/TRAIT/",traits[i],"/ run_config_for_inverse_varience_TEMPLATE.txt > ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))
## system(paste("./metal ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))

system(paste("sed s/TRAIT/common.extra.","logistic","/ run_config_for_inverse_varience_TEMPLATE.txt > ",paste("STDERR.CONFIG.common.extra","logistic","txt",sep="."),sep=""))
system(paste("./metal ",paste("STDERR.CONFIG.common.extra","logistic","txt",sep="."),sep=""))


} ## loop over traits


  ############# some testing

################# annotate ##########
combined.bim<-rbind(bim.chip,bim.seq)
dups<-duplicated(combined.bim[,"SNP"])
combined.bim<-combined.bim[!dups,]
combined.bim[1:5,]
write.table(combined.bim,file="combined.bim",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

###!!!!!!!  run /media/scratch/GWAS_Figures/Annotate_bim_file.r now seperately

################# annotate ##########





write.table(assoc.ann,file="final_combined_annotated.assoc",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
#write.table(assoc.ann[1:10,],file="final_annotated2.assoc",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

ann<-read.table("/media/UQCCG/GWAS/Scleroderma/scl_2014/metal/final_combined_annotated.assoc",header=T,fill=TRUE,stringsAsFactors=FALSE)
ann[1:5,]
  
meta.stderr<-read.table("/media/UQCCG/GWAS/Scleroderma/scl_2014/metal/META.STDERR.common.extra.logistic1.tbl",header=T,fill=TRUE,stringsAsFactors=FALSE)
meta.size<-read.table("/media/UQCCG/GWAS/Scleroderma/scl_2014/metal/META.SAMPLE.SIZE.common.extra.logistic1.tbl",header=T,fill=TRUE,stringsAsFactors=FALSE)

meta.stderr[1:5,]
meta.size[1:5,]
logistic[1:5,]
logistic.seq[1:5,]

data<-meta.stderr

extra<-logistic[,c("CHR","SNP","POS","A1","OR","SE","P")]
extra<-logistic.seq[,c("CHR","SNP","POS","A1","OR","SE","P")]

data[1:5,]


posns<-match(data[,"MarkerName"],logistic[,"SNP"])
missing<-is.na(posns)
sum(missing)
extra<-logistic[posns,c("CHR","SNP","POS","A1","OR","SE","P")]

posns<-match(data[,"MarkerName"],logistic.seq[,"SNP"])
missing<-is.na(posns)
sum(missing)
extra.seq<-logistic.seq[posns,c("CHR","SNP","POS","A1","OR","SE","P")]



data[1:5,]
extra[1:5,]
extra.seq[1:5,]

summary<-cbind(extra[,c("CHR","POS")],data,extra[,"P"],extra.seq[,"P"],extra[,"OR"],extra.seq[,"OR"])


dim(summary)
colnames(summary)[10:13]<-c("P.GWAS","P.Rep","OR.GWAS","OR.rep")

summary[1:5,]
extra.seq[missing.in.extra,][1:5,]
missing.in.extra<-is.na(summary[,"CHR"])
sum(missing.in.extra)
summary[missing.in.extra,c("CHR","POS")]<-extra.seq[missing.in.extra,c("CHR","POS")]
missing.in.extra<-is.na(summary[,"CHR"])
sum(missing.in.extra)

colnames(summary)[colnames(summary)=="MarkerName"]<-"SNP"


posns<-match(summary[,"SNP"],ann[,"ID"])
missing<-is.na(posns)
sum(missing)
ann[1:5,]
extra.ann<-ann[posns, c("Func.refGene","Gene.refGene","ExonicFunc.refGene","AAChange.refGene")]
summary<-cbind(summary,extra.ann)

summary[1:5,]

sig<-summary[,"P.value"]<1e-3 # | summary[,"P.GWAS"]<1e-4 | summary[,"P.Rep"]<1e-4 &   is.na(summary[,"P.GWAS"]) | is.na(summary[,"P.Rep"])
sum(sig)

summary[sig,][1:10,]



meta.stderr2<-read.table("/media/UQCCG/GWAS/Scleroderma/scl_2014/metal/",header=T,fill=TRUE,stringsAsFactors=FALSE)

posns<-match(summary[,"SNP"],meta.stderr2[,"SNP"])
missing<-is.na(posns)
sum(missing)
extra.seq<-logistic.seq[posns,c("CHR","SNP","POS","A1","OR","SE","P")]



  
write.table(summary[sig,],file="meta.signif.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(summary,file="meta.all.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)






















output.dir<-"/media/UQCCG/GWAS/Scleroderma/scl_2014/metal"


setwd(output.dir)

bim.file<-"/media/UQCCG/GWAS/Scleroderma/scl_2014/metal/scl_wtcc_merge_common.bim"
bim<-read.table(bim.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
colnames(bim)<-c("chr","SNP","cm","start","ALT","REF")
bim.chip<-bim


bim.file<-"/media/UQCCG/GWAS/Scleroderma/scl_2014/metal/dbgap_aogc_10kmerge_imputed_omni.bim"
bim<-read.table(bim.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
colnames(bim)<-c("chr","SNP","cm","start","ALT","REF") # A1 is labelled ALT here
bim.seq<-bim

bim.chip[1:5,]
bim.seq[1:5,]
dim(bim.seq)
dim(bim.chip)


#################################################

logistic.seq<-read.table("/media/UQCCG/GWAS/Scleroderma/scl_2014/metal/scl_wtcc_merge_common_assoc.assoc.logistic",header=T,fill=TRUE,stringsAsFactors=FALSE)
logistic.seq[1:5,]
colnames(logistic.seq)[colnames(logistic.seq)=="BP"]<-"POS"

if("TEST" %in% colnames(logistic.seq)){
wanted<-grepl("ADD",logistic.seq[,"TEST"])
logistic.seq<-logistic.seq[wanted,]
}
dim(logistic.seq)
dim(bim.seq)


posns<-match(logistic.seq[,"SNP"],bim.seq[,"SNP"])
missing<-is.na(posns)
sum(missing)
bim.seq<-bim.seq[posns,]

junk<-bim.seq[,"start"]==0 | bim.seq[,"chr"]==0 | is.na(bim.seq[,"chr"])
sum(junk)
      
bim.seq<-bim.seq[!junk,]
logistic.seq<-logistic.seq[!junk,]


      
key.logistic.seq<-build.key(bim.seq,core.ann)
dups<-duplicated(key.logistic.seq)
      sum(dups)

rownames(logistic.seq)<-key.logistic.seq
TYPE<-rep("SEQ",times=length(key.logistic.seq))
#logistic.seq<-cbind(logistic.seq,TYPE,bim.seq[, core.ann],stringsAsFactors=FALSE)
logistic.seq<-cbind(logistic.seq,TYPE,bim.seq,stringsAsFactors=FALSE)
logistic.seq[1:5,]

#############################################
 logistic<-read.table("/media/UQCCG/GWAS/Scleroderma/scl_2014/metal/scl_wtcc_merge_common_assoc.assoc.logistic",header=T,fill=TRUE,stringsAsFactors=FALSE)
logistic[1:5,]
colnames(logistic)[colnames(logistic)=="BP"]<-"POS"

if("TEST" %in% colnames(logistic)){
wanted<-grepl("ADD",logistic[,"TEST"])
logistic<-logistic[wanted,]
}
dim(logistic)
dim(bim.chip)


posns<-match(logistic[,"SNP"],bim.chip[,"SNP"])
missing<-is.na(posns)
sum(missing)
bim.chip<-bim.chip[posns,]

junk<-bim.chip[,"start"]==0 | bim.chip[,"chr"]==0 | is.na(bim.chip[,"chr"])
sum(junk)

bim.chip[junk,][1:5,]
logistic[junk,][1:5,]
      
bim.chip<-bim.chip[!junk,]
logistic<-logistic[!junk,]


      
key.logistic<-build.key(bim.chip,core.ann)
dups<-duplicated(key.logistic)
      sum(dups)

rownames(logistic)<-key.logistic
TYPE<-rep("GENO",times=length(key.logistic))
#logistic<-cbind(logistic,TYPE,bim.chip[, core.ann],stringsAsFactors=FALSE)
logistic<-cbind(logistic,TYPE,bim.chip,stringsAsFactors=FALSE)

      #############################


dim(logistic.seq)
dim(logistic)


logistic.seq[1:5,]
logistic[1:5,]

#################### get orderso can compare efefct directions

## posns<-match(logistic.seq[,"SNP"],logistic[,"SNP"])
##    missing<-is.na(posns)
## sum(missing)

## ## logistic.miss<-logistic[posns[missing],]
## ## logistic<-logistic[posns[!missing],]
## ## logistic.seq.miss<-logistic.seq[missing,]
## ## logistic.seq<-logistic.seq[!missing,]      

## dim(logistic.seq)
## dim(logistic)

##  logistic.seq[1:5,]
## logistic[1:5,]     
## ##################
#########
if("A1" %in% colnames(logistic)){ ### make sure ALT is A1 so effect direction is consistant - dont need to change sign of dirn then


  ################ set ALT and A! for the same
ref.diff<-logistic[,"A1"]!=logistic[,"ALT"]
sum(ref.diff)
logistic[ref.diff,][1:10,]
the.ref<-logistic[ref.diff,"REF"]
logistic[ref.diff,"REF"]<-logistic[ref.diff,"ALT"]
logistic[ref.diff,"ALT"]<-the.ref
ref.diff<-logistic[,"A1"]!=logistic[,"ALT"]
sum(ref.diff)

ref.diff<-logistic.seq[,"A1"]!=logistic.seq[,"ALT"]
sum(ref.diff)
logistic.seq[ref.diff,][1:10,]
the.ref<-logistic.seq[ref.diff,"REF"]
logistic.seq[ref.diff,"REF"]<-logistic.seq[ref.diff,"ALT"]
logistic.seq[ref.diff,"ALT"]<-the.ref
ref.diff<-logistic.seq[,"A1"]!=logistic.seq[,"ALT"]
sum(ref.diff)
#########################

############## don't chane the effect direction just let it go... meatl can sort out
## sum(is.na(logistic.seq[,"A1"]))
## remove<-is.na(logistic.seq[,"A1"]) | is.na(logistic.seq[,"ALT"]) |  is.na(logistic.seq[,"REF"]) # crap results
## sum(remove)
## logistic.seq<-logistic.seq[!remove,]

## ref.diff<-logistic.seq[,"A1"]!=logistic[,"A1"] 
## sum(ref.diff)
## logistic.seq[ref.diff,][1:10,]
## the.ref<-logistic.seq[ref.diff,"REF"]
## logistic.seq[ref.diff,"REF"]<-logistic.seq[ref.diff,"ALT"]
## logistic.seq[ref.diff,"ALT"]<-the.ref
## ref.diff<-logistic.seq[,"A1"]!=logistic.seq[,"ALT"]
## sum(ref.diff)

# redo.start.end.annovar 

}

### using minor allele now so should be fineetal will sort out differences




dim(logistic.seq)
dim(logistic)

logistic.seq[1:5,]
logistic[1:5,]
################### match up snp ids:
## posns<-match(logistic.seq[,"SNP"],logistic[,"SNP"])
## missing<-is.na(posns)
## sum(missing)
##  test<-logistic[,"P"]< 1e-6
##       sum(test)
## logistic[test,]
      
 


if(!("BETA" %in% colnames(logistic))){
BETA=log(as.numeric(logistic[,"OR"]))
logistic<-cbind(logistic,BETA,stringsAsFactors=FALSE)

## exp(BETA +- 1.96SE)=95% CI
## so
## se<--(log(as.numeric(logistic[,"L95"]))-BETA)/1.96
## logistic<-cbind(logistic,BETA,se,stringsAsFactors=FALSE)
## and se = SE in plink so !! SE Is the standard error of BETA!! no need to mess with
}

if(!("BETA" %in% colnames(logistic.seq))){
BETA=log(as.numeric(logistic.seq[,"OR"]))
logistic.seq<-cbind(logistic.seq,BETA,stringsAsFactors=FALSE)

#logistic.seq[,"SE"]<-log(as.numeric(logistic.seq[,"SE"]))
}
print("PAST BETA")
logistic[1:5,]
logistic.seq[1:5,]

setwd(output.dir)

### below is commented as just doing common
## write.table(logistic,file=paste("chip",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(logistic.seq,file=paste("seq",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

## write.table(logistic[posns[!missing],],file=paste("chip.common.extra",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(logistic.seq[!missing,],file=paste("seq.common.extra",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

write.table(logistic,file=paste("chip.common.extra","logistic2",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(logistic.seq,file=paste("seq.common.extra","logistic2",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



 logistic.seq[1:5,]
logistic[1:5,]   
##print(doing merge)
## common.cols<-colnames(logistic)[colnames(logistic) %in% colnames(logistic.seq)]
## the.merge<-rbind(logistic[,common.cols],logistic.seq[,common.cols])
## the.merge[1:5,]
## order.by<-order(the.merge[,"CHR"],the.merge[,"start"])
## the.merge<-the.merge[order.by,]
## the.merge[1:5,]
## write.table(the.merge,file=paste("merge.extra",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
##

## system(paste("sed s/TRAIT/",traits[i],"/ run_config_for_P_and_BETA_TEMPLATE.txt > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
## system(paste("./metal ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))

system(paste("sed s/TRAIT/common.extra.","logistic2","/ run_config_for_P_and_BETA_TEMPLATE.txt > ",paste("sample.size.CONFIG.common.extra","logistic2","txt",sep="."),sep=""))
system(paste("./metal ",paste("sample.size.CONFIG.common.extra","logistic2","txt",sep="."),sep=""))


## system(paste("sed s/TRAIT/",traits[i],"/ run_config_for_inverse_varience_TEMPLATE.txt > ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))
## system(paste("./metal ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))

system(paste("sed s/TRAIT/common.extra.","logistic2","/ run_config_for_inverse_varience_TEMPLATE.txt > ",paste("STDERR.CONFIG.common.extra","logistic2","txt",sep="."),sep=""))
system(paste("./metal ",paste("STDERR.CONFIG.common.extra","logistic2","txt",sep="."),sep=""))













write.table(summary,file="meta.all.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

summary<-read.table("/media/UQCCG/GWAS/Scleroderma/scl_2014/metal/meta.all.txt",header=T,fill=TRUE,stringsAsFactors=FALSE)


meta.stderr2<-read.table("/media/UQCCG/GWAS/Scleroderma/scl_2014/metal/META.STDERR.common.extra.logistic21.tbl",header=T,fill=TRUE,stringsAsFactors=FALSE)
meta.stderr2[1:5,]
posns<-match(summary[,"SNP"],meta.stderr2[,"MarkerName"])
missing<-is.na(posns)
sum(missing)
extra.seq2<-meta.stderr2[posns,]




posns<-match(summary[,"SNP"],logistic[,"SNP"])
missing<-is.na(posns)
sum(missing)
logistic<-logistic[posns,]

extra.seq2[1:5,]
summary[1:5,]
logistic[1:5,]

replace<-!is.na(extra.seq2[,"P.value"]) & (extra.seq2[,"P.value"] < summary[,"P.value"] ) & !is.na(logistic[,"P"])

reject<-extra.seq2[,"P.value"] < 1e-20
sum(replace & !reject)
replace<-replace & !reject
cbind(summary,extra.seq2)[replace,][1:10,]
P.old<-summary[,"P.value"]
P.GWAS.old<-summary[,"P.GWAS"]


summary[replace,"P.value"]<-extra.seq2[replace,"P.value"]
summary[replace,"P.GWAS"]<-logistic[replace,"P"]

summ




summary[replace,"P.value"]<-extra.seq2[replace,"P.value"]
summary2<-cbind(summary,P.old,P.GWAS.old)

sig<-summary2[,"P.value"]<1e-3 # | summary[,"P.GWAS"]<1e-4 | summary[,"P.Rep"]<1e-4 &   is.na(summary[,"P.GWAS"]) | is.na(summary[,"P.Rep"])
sum(sig)

summary2[sig,][1:10,]

write.table(summary2[sig,],file="meta.signif.final.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(summary2,file="meta.all.final.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
sig<-summary2[,"P.value"]<0.1 # | summary[,"P.GWAS"]<1e-4 | summary[,"P.Rep"]<1e-4 &   is.na(summary[,"P.GWAS"]) | is.na(summary[,"P.Rep"])
sum(sig)

write.table(summary2[sig,],file="meta.signif.final.less.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



info<-read.table("/media/UQCCG/GWAS/Scleroderma/scl_2014/metal/AOGC_uk10k_info.txt",header=T,fill=TRUE,stringsAsFactors=FALSE)

info[1:5,]
posns<-match(summary[,"SNP"],info[,"SNP"])
missing<-is.na(posns)
sum(missing)

info
