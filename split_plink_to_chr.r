
##################################################
work.dir<-"/media/UQCCG-Analysis/AOGC_exome_chip/working_genotypes"
setwd(work.dir)
bim.file<-c("recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL_withRS.bim") ## list all bim files to get common snp basis
##################################################ALL_eth_650Y_forward_hg19_Best

#/media/UQCCG-Analysis/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.bim
num.vars<-6; skip.lines<-0
all.snps<-{}
# i<-2

a.bim<-try(scan(bim.file,what=character(num.vars),skip=skip.lines,fill=TRUE,))
num.lines<-length(a.bim)/(num.vars)
dim(a.bim)<-c(num.vars,num.lines)
a.bim<-t(a.bim)

the.chr<-unique(a.bim[,1])
file.prefix<-gsub(".bim$","",bim.file)
i<-1
print(the.chr)
print(dim(a.bim))
for(i in 1:length(the.chr)){

system( paste("plink --bfile ",file.prefix , " --chr ",the.chr[i]," --make-bed --out ", paste(file.prefix,"_chr",the.chr[i],sep=""),sep="") )

}



