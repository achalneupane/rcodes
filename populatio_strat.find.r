

1)
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr20.output.recalibrated.filtered.vcf --remove-filtered-geno-all --plink  --out chr20
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr21.output.recalibrated.filtered.vcf  --remove-filtered-geno-all  --plink  --out chr21
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr22.output.recalibrated.filtered.vcf  --remove-filtered-geno-all  --plink  --out chr22
## vcftools --vcf 2015-08-14_AML_mixedAligners.chrX.output.recalibrated.filtered.vcf --remove-filtered-geno-all   --plink  --out chrX
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr1.output.recalibrated.filtered.vcf  --remove-filtered-geno-all  --plink  --out chr1
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr2.output.recalibrated.filtered.vcf  --remove-filtered-geno-all  --plink  --out chr2
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr3.output.recalibrated.filtered.vcf  --remove-filtered-geno-all  --plink  --out chr3
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr4.output.recalibrated.filtered.vcf  --remove-filtered-geno-all  --plink  --out chr4
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr5.output.recalibrated.filtered.vcf  --remove-filtered-geno-all  --plink  --out chr5
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr6.output.recalibrated.filtered.vcf  --remove-filtered-geno-all  --plink  --out chr6
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr7.output.recalibrated.filtered.vcf  --remove-filtered-geno-all  --plink  --out chr7
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr8.output.recalibrated.filtered.vcf  --remove-filtered-geno-all  --plink  --out chr8
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr9.output.recalibrated.filtered.vcf  --remove-filtered-geno-all  --plink  --out chr9
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr10.output.recalibrated.filtered.vcf  --remove-filtered-geno-all   --plink  --out chr10
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr11.output.recalibrated.filtered.vcf --remove-filtered-geno-all   --plink  --out chr11
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr12.output.recalibrated.filtered.vcf  --remove-filtered-geno-all  --plink  --out chr12
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr13.output.recalibrated.filtered.vcf  --remove-filtered-geno-all  --plink  --out chr13
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr14.output.recalibrated.filtered.vcf  --remove-filteredR:4-geno-all   --plink  --out chr14
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr15.output.recalibrated.filtered.vcf --remove-filtered-geno-all   --plink  --out chr15
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr16.output.recalibrated.filtered.vcf  --remove-filtered-geno-all  --plink  --out chr16
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr17.output.recalibrated.filtered.vcf  --remove-filtered-geno-all  --plink  --out chr17
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr18.output.recalibrated.filtered.vcf  --remove-filtered-geno-all  --plink  --out chr18
## vcftools --vcf 2015-08-14_AML_mixedAligners.chr19.output.recalibrated.filtered.vcf  --remove-filtered-geno-all   --plink  --out chr19
## vcftools --vcf 2015-08-14_AML_mixedAligners.chrY.output.recalibrated.filtered.vcf  --remove-filtered-geno-all  --plink  --out chrY
## /media/UQCCG/GWAS/Popularion stratification/ALL_eth_650Y_forward_hg19_Best.bim


plink --file chr10 --merge-list file.csv --make-bed --out aml  ## merge multiple map/ped files into a bed/bim/fam

 plink --bfile aml  --maf 0.02 --geno 0.1 --make-bed --ou t aml.f  ## merge multiple file 
plink --bfile aml.f  --freq  # 219560

plink --bfile ALL_eth_650Y_forward_hg19_Best  --freq


.bim




#########################################################
code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")
source("hwe.r")

#######Do one at a time

bim.file.list<-c("ALL_eth_650Y_forward_hg19_Best.bim")
work.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/SNPs"
setwd(work.dir)

bim.file.list<-c("aml.f.bim")
work.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/SNPs"
setwd(work.dir)

#####################



getwd()

i<-1
num.vars<-6; skip.lines<-0
all.snps<-{}
a.bim<-try(scan(bim.filels.list[i],what=character(num.vars),skip=skip.lines,fill=TRUE,))
num.lines<-length(a.bim)/(num.vars)
dim(a.bim)<-c(num.vars,num.lines)
a.bim<-t(a.bim)

colnames(a.bim)<-c("chr","SNP","CM","start","A1","A2")

a.bim[1:5,]
key<-build.key(a.bim,c("chr","start"),add.chr=TRUE)

key[1:5]
a.bim[,"SNP"]<-key
a.bim[1:5,]
write.table(a.bim,bim.file.list[i],col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t",append=FALSE)

#################################################


work.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/SNPs"
setwd(work.dir)
bim.file.list<-c("ALL_eth_650Y_forward_hg19_Best.bim","aml.f.bim") ## list all bim files to get common snp basis
##################################################ALL_eth_650Y_forward_hg19_Best

##################################################
work.dir<-"/media/UQCCG/GWAS/Cervical_Cancer/Imputation_1000G_Impute2_common_snp/Strand_fixed_QC"
setwd(work.dir)
bim.file.list<-c("ALL_eth_650Y_forward_hg19_Best.bim","wtcc.omni.660.QC.final.bim") ## list all bim files to get common snp basis
##################################################ALL_eth_650Y_forward_hg19_Best


num.vars<-6; skip.lines<-0
all.snps<-{}
# i<-2
combined.bim<-{}
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
if(is.null(dim(combined.bim))){ combined.bim<-a.bim}else{combined.bim<-rbind(combined.bim,a.bim)}
}

colnames(combined.bim)<-c("chr","SNP","CM","start","A1","A2")
combined.bim[1:5,]
counts<-tapply(all.snps,all.snps,length)
in.common<-counts==length(bim.file.list)
print(paste(sum(in.common)," snps in common for all bims"))

combined.bim.common<-combined.bim[ combined.bim[,"SNP"] %in% names(counts[in.common]) ,]
dim(combined.bim.common)


################################# extra code I needed to get rid of 
counts[1:5]
all.snps[1:5,]
combined.bim[1:5,]
colnames(combined.bim)<-c("chr","SNP","CM","start","A1","A2")

dim(combined.bim.common)


order.by<-order(combined.bim.common[,"SNP"])
combined.bim.common<-combined.bim.common[order.by,]
combined.bim.common[1:6,]

unique.snps<-unique(combined.bim.common[,"SNP"])
dups<-duplicated(combined.bim.common[,"SNP"])
dups[1:5]

table.common<-cbind(combined.bim.common[!dups,],combined.bim.common[dups,c("A1","A2")])
table.common[1:5,]
test<-apply(table.common[,5:8],1,function(x) length(unique(x)) )
test[1:10]
names(test)<-table.common[,"SNP"]
bad.alleles<-names(test)[test>2]

keep<-names(counts[in.common])[!(names(counts[in.common]) %in% bad.alleles)]
length(keep)
combined.bim[,"SNP"]
####################################
keep[1:50]

file.list<-gsub(".bim$","",bim.file.list)
file.list
fam<-read.table(paste(file.list,"fam",sep=".")[1],header=F,fill=TRUE,stringsAsFactors=FALSE)

fam[1:5,]
table(fam[,6])
## posns<-match(fam[,1],score[,"IID")
## missing<-is.na(posns)
## sum(missing)

write.table(keep,"common_snps_Chip.650.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t",append=FALSE)

#### extract coommon snps
bim.file.list<-c("ALL_eth_650Y_forward_hg19_Best.bim","wtcc.omni.660.QC.final.bim")
# c("ALL_eth_650Y_forward_hg19_Best.bim","wtcc.omni.660.QC.final.bim")

i<-1
system(paste("plink --bfile", file.list[i]," --extract common_snps_Chip.650.txt --allow-no-sex --make-bed --out  ",paste(file.list[i],"C",sep=".")))
i<-2
system(paste("plink --bfile", file.list[i]," --extract common_snps_Chip.650.txt --allow-no-sex --make-bed --out  ",paste(file.list[i],"C",sep=".")))



### merge files
system(paste("plink --bfile",  paste(file.list[1],"C",sep=".")," --bmerge", paste(file.list[2],"C","bed",sep="."), paste(file.list[2],"C","bim",sep="."), paste(file.list[2],"C","fam",sep="."),"--allow-no-sex --flip-scan",sep=" "))

system("plink --bfile",paste(file.list[1],"C",sep="."),"--flip --allow-no-sex plink.missnp --make-bed --out ",paste(file.list[1],"C","F",sep=".")

system(paste("plink --bfile",  paste(file.list[1],"C",sep=".")," --bmerge", paste(file.list[2],"C","bed",sep="."), paste(file.list[2],"C","bim",sep="."), paste(file.list[2],"C","fam",sep="."),"--allow-no-sex --make-bed  --out merge ",sep=" "))



merge.name<-"merge"

bim<-read.table(paste(merge.name,"bim",sep="."),header=F,fill=TRUE,stringsAsFactors=FALSE)
table(bim[,1])
bim[1:5,]
write.table(bim[bim[,1] %in% 23,2],"chrX.snps.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t",append=FALSE)
#### extract long range LD positions 
system(paste("plink --bfile", merge.name," --allow-no-sex --exclude chrX.snps.txt --make-bed --out ",paste(merge.name,"noX",sep="."),sep=" "))

system(paste("plink --bfile",paste(merge.name,"noX",sep="."),"--noweb --allow-no-sex --exclude Price2008_hg19_paul.txt --range --make-bed --out ",paste(merge.name,"noX","F",sep="."),sep=" "))

system(paste("plink --bfile",paste(merge.name,"noX","F",sep=".")," --indep 50 5 1.5 --out ldtrimset",sep=" "))

### LD prune
system(paste("plink --bfile",paste(merge.name,"noX","F",sep=".")," --noweb --allow-no-sex --exclude ldtrimset.prune.out --make-bed --out ",paste(merge.name,"noX","F","P",sep="."),sep=" "))
## THin
system(paste("plink --bfile",paste(merge.name,"noX","F","P",sep=".")," --noweb --allow-no-sex --thin .9  --make-bed --out ",paste(merge.name,"noX","F","P","T",sep="."),sep=" "))



system(paste("plink --bfile",paste(merge.name,"noX","F","P","T",sep=".")," --noweb --allow-no-sex --maf 0.05 --geno 0.03 --hwe 0.0000001 --mind 0.1 --make-bed --out ",paste(merge.name,"noX","F","P","T","QC",sep="."),sep=" "))


system(paste("/media/scratch/Shellfish/shellfish/shellfish.py --pca --numpcs 10 --maxprocs 8 --file", paste(merge.name,"noX","F","P","T","QC",sep=".")," --out ",paste(merge.name,"QCed",sep="."),sep=" ") )

#####################run plot_650_strat.r

"merge.QCed.evecs_keep.txt"

system(paste("plink --bfile",paste(merge.name,"noX","F","P","T","QC",sep="."),"--noweb --allow-no-sex --keep merge.QCed.evecs_keep.txt --make-bed --out ",paste(merge.name,"QCed","FINAL",sep="."),sep=" "))

system(paste("plink --bfile",paste(merge.name,"QCed","FINAL",sep=".")," --noweb --allow-no-sex --maf 0.05 --geno 0.03 --hwe 0.0000001 --mind 0.1 --make-bed --out ",paste(merge.name,"QCed","FINAL","QC",sep="."),sep=" ")) #2683 are cases and 6435 are controls


system(paste("/media/scratch/Shellfish/shellfish/shellfish.py --pca --numpcs 10 --maxprocs 8 --file", paste(merge.name,"QCed","FINAL","QC",sep=".")," --out ",paste(merge.name,"QCed","FINAL","QC",sep="."),sep=" ") )




2683 are cases and 6435 are controls.
"merge.QCed.evecs_keep.txt"

$$$
/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/SNPs/aml.evecs_keep_6SD.txt

plink --bfile data_WTCCC_f_650   --keep aml.evecs_keep_6SD.txt --make-bed --out data_WTCCC_f_650.ETH --noweb

/media/scratch/Shellfish/shellfish/shellfish.py --pca --numpcs 10 --maxprocs 8 --file  data_WTCCC_f_650.ETH --out aml.ETH
