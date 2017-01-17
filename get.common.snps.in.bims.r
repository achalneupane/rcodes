

## file<-"/media/old-scratch/media/scratch2/X_chromo/ImmunochipXYcs/brisbane_data_gencall_23112012.fam"

## file<-"/media/old-scratch/media/scratch2/X_chromo/ImmunochipXYcs/PostQC_Immunochip_fam_file.fam"

## fam2<-read.table(file,header=F,fill=TRUE,stringsAsFactors=FALSE)

## gsub("$ ","",fam[1

## posns<-match(fam1[,1],fam2[,1])
## missing<-is.na(posns)
## sum(missing)

## fam1[missing,][1:500,]

## strange<-(fam1[,1]!=fam1[,2])
## fam1[strange,]
## table(fam1[,3])





#####################################################
work.dir<-"/media/scratch2/GBS"
setwd(work.dir)
bim.file.list<-c("wtcc2_forward_f.bim","GBS_hg19_clean_f.bim") ## list all bim files to get common snp basis
bim.file.list<-c("GBS_wtccc_commom_ld_f_thin_f_clean.bim","ALL_eth_650Y_forward_hg19_clean.bim") ## list all bim files to get common snp basis
##################################################

#####################################################
work.dir<-"/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/PCA_calc"
setwd(work.dir)
bim.file.list<-c("AOGC_merge_common._Best.bim","recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL_withRS.bim","ALL_eth_650Y_forward_hg19_Best.bim") ## list all bim files to get common snp basis
##################################################

#####################################################
work.dir<-"/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/PCA_calc"
setwd(work.dir)
bim.file.list<-c("AOGC_merge_common._Best.bim","recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL_withRS.bim","HBM_Best_Unique.noExomeC.bim") ## list all bim files to get common snp basis

##################################################
work.dir<-"/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/PCA_calc"
setwd(work.dir)
bim.file.list<-c("ALL_eth_650Y_forward_hg19_Best.bim","AOGC_merge_common._Best.bim","recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL_withRS.bim","HBM_Best_Unique.noExomeC.bim") ## list all bim files to get common snp basis
##################################################ALL_eth_650Y_forward_hg19_Best

##################################################
work.dir<-"/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/PCA_calc"
setwd(work.dir)
bim.file.list<-c("ALL_eth_650Y_forward_hg19_Best.bim","exomesChip.rs.f.ld.c.bim") ## list all bim files to get common snp basis
##################################################ALL_eth_650Y_forward_hg19_Best

##################################################
work.dir<-"/media/UQCCG/GWAS/Cervical_Cancer/Imputation_1000G_Impute2_common_snp/Strand_fixed_QC"
setwd(work.dir)
bim.file.list<-c("ALL_eth_650Y_forward_hg19_Best.bim","wtcc.omni.660.QC.final.bim") ## list all bim files to get common snp basis
##################################################ALL_eth_650Y_forward_hg19_Best


num.vars<-6; skip.lines<-0
all.snps<-{}
# i<-2

for(i in 1:length(bim.file.list)){

a.bim<-try(scan(bim.file.list[i],what=character(num.vars),skip=skip.lines,fill=TRUE,))
num.lines<-length(a.bim)/(num.vars)
dim(a.bim)<-c(num.vars,num.lines)
a.bim<-t(a.bim)
if(sum(grepl("^RS",a.bim[,2]))>0){
  print(paste("warning",bim.file.list[i]," contains RS NOT rs ids use -sed -i '/RS/rs/' file.bim . Converting...",sep=" "))
  a.bim[,2]<-gsub("^RS","rs",a.bim[,2])
      }
all.snps<-c(all.snps,a.bim[,2])
a.bim[1:5,]
}

counts<-tapply(all.snps,all.snps,length)
in.common<-counts==length(bim.file.list)
print(paste(sum(in.common)," snps in common for all bims"))

write.table(names(counts[in.common]),"common_snps_Chip.650.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t",append=FALSE)
