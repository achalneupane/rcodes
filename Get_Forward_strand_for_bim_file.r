


###########################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-TRUE
plot.dir<-"/media/scratch2/GBS"
bim.file<-"/media/scratch2/GBS/GBS_hg19.bim"  ### bim file that has the path as well
bim.file.cols<-c("chr","SNP","AF","POS","A1","A2")
bim.has.header<-FALSE
genome.build<-"hg19"
#################################################
code.dir<-"/media/Bioinform-D/Research/AML sequencing"
source("annotate_SNPs_subroutines.r")

sample.files<-basename(bim.file)
names(sample.files)<-"GENO"
snp.dir<-dirname(bim.file)

sample.grs<-paste(names(sample.files),"grs",sep=".")
setwd(code.dir)

source("read.plink.bim.files.r")
indels<-GENO.grs
indels[1:5,]

