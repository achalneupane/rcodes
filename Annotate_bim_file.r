


###########################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-TRUE
plot.dir<-"/media/scratch2/GBS"
bim.file<-"/media/scratch2/GBS/GBS_hg19.bim"  ### bim file that is the superset of all genotypes
bim.file.cols<-c("chr","SNP","AF","POS","A1","A2")
bim.has.header<-FALSE

#################################################


###########################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-TRUE
plot.dir<-"/media/UQCCG/GWAS/Scleroderma/scl_2014/metal"
bim.file<-"/media/UQCCG/GWAS/Scleroderma/scl_2014/metal/combined.bim"  ### bim file that is the superset of all genotypes
bim.file.cols<-c("chr","SNP","AF","POS","A1","A2")
bim.has.header<-FALSE

#################################################

###########################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-FALSE
plot.dir<-"/media/scratch/immunoChip"
bim.file<-"/media/scratch/immunoChip/QC1_QC2_pass.bim"  ### bim file that is the superset of all genotypes
bim.file.cols<-c("chr","SNP","AF","POS","A1","A2")
bim.has.header<-FALSE

#################################################
###########################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-FALSE
plot.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes"
bim.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.bim"  ### bim file that is the superset of all genotypes
bim.file.cols<-c("chr","SNP","AF","POS","A1","A2")
bim.has.header<-FALSE

#################################################

bim<-read.table(bim.file,header=bim.has.header,fill=TRUE,stringsAsFactors=FALSE)

if(!bim.has.header){
colnames(bim)<-bim.file.cols
}


dim(bim)
bim[1:2,]



############### STOP give "final.eigen.assoc" to katie


##################################
core.ann<-c("chr","start","end","REF","ALT","TYPE")
code.dir<-"/mnt/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"## location of "annotate_SNPs_subroutines.r" & ("get.forward.strand.allele.r")
genome.build<-"hg19"
target<-"temp"

bim[1:5,]
indels<-cbind(bim[,c("chr","POS","POS","A1","A2")],"snp",bim[,"SNP"])
indels[1:5,]
colnames(indels)<-c(core.ann,"SNP")
indels[1:5,]
## > indels[1:5,]
##    chr  start    end REF ALT TYPE
## 1 chr1  72017  72017   0   A  snp
## 2 chr1 752566 752566   G   A  snp
## 3 chr1 752721 752721   A   G  snp
## 4 chr1 768448 768448   A   G  snp
## 5 chr1 776546 776546   G   A  snp

geneanno.DB<-c("refGene")  # returns 2 extra columns c("refGene","knownGene","ensGene")
names(geneanno.DB)<-c("refGene") # c("refGene","knownGene","ensGene")

filter.DB<-c("snp138")  # returns 2  extra columns (DB,score) c("snp132","1000g2011may_all")
names(filter.DB)<-c("ID")

## setwd(code.dir)
## source("ucsc.table.names.r")   # load in the UCSC tables these use the db file names and not their lable-names 
## source("ucsc.table.names.processor.r")
setwd(code.dir)
source("annotate_SNPs_subroutines.r")
source("get.forward.strand.allele.r") # this one takes a little while to run

annovar.loc<-"/media/UQCCG/Software/annovar"
anno.DB.location.core<-"/media/UQCCG/Software/annovar/humandb" 
anno.DB.location<-paste(anno.DB.location.core,genome.build,sep="/") #  "/media/Bioinform-D/Research/annovar/humandb/hg19"


indels[,"start"]<-gsub("^\\s+","",indels[,"start"])
indels[,"end"]<-gsub("^\\s+","",indels[,"end"])
indels[,"chr"]<-gsub("^\\s+","",indels[,"chr"])

setwd(plot.dir)
write.table(indels,file=target,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)


#table_annovar.pl ex1.human humandb/ -protocol refGene,phastConsElements44way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp135,avsift,ljb_all -operation g,r,r,f,f,f,f,f -nastring NA
system(paste(annovar.loc,"/","table_annovar.pl ",target,"  ",anno.DB.location," -buildver ",genome.build," -protocol refGene,snp1write.table(assoc.ann[1:10,],file="final_annotated2.assoc",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)










