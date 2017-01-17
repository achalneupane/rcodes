
#####################################################
work.dir<-"/media/scratch2/GBS"
setwd(work.dir)
bim.file.list<-c("wtcc2_forward_f.bim","GBS_hg19_clean_f.bim") ## list all bim files to get common snp basis
bim.file.list<-c("GBS_wtccc_commom_ld_f_thin_f_clean.bim","ALL_eth_650Y_forward_hg19_clean.bim") ## list all bim files to get common snp basis
##################################################

#####################################################
work.dir<-"/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/PCA_calc2"
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


num.vars<-6; skip.lines<-0
all.snps<-{}
all.snps.kry<-{}
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


####
colnames(a.bim)<-c("chr","start","cm","A1","A2")
a.key<-build.key(a.bim,c(c("chr","start","A1","A2"))
all.snps.key<-c(all.snps.key,a.key)                 
###                 

all.snps<-c(all.snps,a.bim[,2])

                 
a.bim[1:5,]
}

counts<-tapply(all.snps,all.snps,length)


in.common<-counts==length(bim.file.list)
print(paste(sum(in.common)," snps in common for all bims"))

write.table(names(counts[in.common]),"common_snps_Chip.650.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t",append=FALSE)



build.key<-function(table,key.cols,add.chr.label=FALSE){
  options(scipen=300)
  if(is.null(dim(table))){table<-as.matrix(table)} # in casea vector sent
  if(length(key.cols)<1){print("FAIL no keys columns specified");key<-1:dim(table)[1]}else{
    for (i in 1:length(key.cols)){
      if(i==1){key<-table[,key.cols[i]]}else{
      key<-paste(key,table[,key.cols[i]],sep=":")}
    }}
  if(add.chr.label){key<-paste("chr",key,sep="")}
         key}
