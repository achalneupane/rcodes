#!/usr/bin/Rscript --vanilla --slave

## Achal-2021-06-15 - UPDATED and MODIFIED for FASe Replication paper from Vicky's -EOAD analysis


### Script to run SKAT over a set of genes WITHIN A CHROMOSOME DIRECTORY
### USAGE: ./SKAT-unrelated-loop-OPTIMAL.R  ${WORKDIR} ${CHR} ${RAW} ${SET} ${COVARS}



##### 1 - Use PLINK to generate raw files,a VCF file and the covars file for the given dataset; this VCF file will go into annotation
#CHR_TO_PROCESS="$(seq 1 22) X Y"
#THREADS=12
#BFILE="FASe-1374-NOLATINS-SNPs-INDELs"
# parallel -j ${THREADS} plink1.9 --bfile ${BFILE} --max-maf 0.05 --keep-allele-order --allow-no-sex --recode A --prune --output-chr MT --chr {} --out chr{}/${BFILE}-MAF5P-chr{} ::: ${CHR_TO_PROCESS}
####  2 - annotate VCF file
#$ SNPEff + SNPSift annotate 
####  3 - generate list of variants of interest
####  4 - generate gene sets Within each chr directory



library(SKAT)
library(tools)
library(data.table)
# Define arguments
args <- commandArgs(TRUE)

# workdir <- "/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/02-gene-based"
workdir<-args[1]

# chr<-"X"
chr<-args[2]

# raw <- paste0("AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean3-chr") 
raw<-args[3]

# set<-"geneset-MAF1PEXAC"
set<-args[4]

# covarsfile <- "/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/02-gene-based/Replication_1590.covars"
covarsfile<-args[5]

# set up directories:
chrdir<-paste0(workdir,"/chr", chr)
# rawfile<-paste0(chrdir,"/",raw,chr,".raw")
rawfile<-paste0(chrdir,"/",raw,chr,".raw")
setwd<-paste0(workdir, "/chr",chr)

#### LOAD RAW file
cat(paste0("\nUploading raw file ",rawfile," \n"))
system.time(raw <- read.table(rawfile, header=T, check.names=F) )
cat(dim(raw))

# arrange header of raw file
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
colnames(raw) = HEADER

## LOAD COVARS file
cat(paste0("\nUploading covariates file ",covarsfile,"\n"))
covars <-read.table(covarsfile, head=T, check.names=F)

#ARRANGE covars so missing data is NA, PCs dont need to be corrected. For this analyisis we are dropping AGE
#covars$AGE[covars$AGE==-9]=NA
covars$APOE[covars$APOE==-9]=NA
covars$APOE[covars$APOE==-0]=NA

#MERGE1=merge(covars,raw, by="IID" ) #use this way and follwoing STEPS IF adding further covariants
com.data=merge(covars,raw, by="IID" )
cat("\n\ndimensions of com.data matrix are:\n")
dim(com.data)

#ARRANGE PHENO and GENDER so values are within  [0,1]
#GENDER
com.data$SEX.y[com.data$SEX.y==1]=0
com.data$SEX.y[com.data$SEX.y==2]=1
com.data$SEX.y[com.data$SEX.y==-9]=NA

#Pheno
com.data$STATUS[com.data$STATUS==-9]=NA
com.data$STATUS[com.data$STATUS==1]=0
com.data$STATUS[com.data$STATUS==2]=1

#################################
#########START LOOP PART 
#################################
## CHANGE DIR
#set<-"geneset-MAF-0.005"
dir<-paste0(chrdir,"/",set)
extension<-".gene"

## Initialize values:
cat("\nInitialize values for table generation:\nGENE\tSNPID\tN0\tN1\tSKAT\tSKATO\tSKAT_125\n")
genelist<-"GENE"
SNPID<-"SNPID"
N0<-"N0"
N1<-"N1"
SKAT<-"SKAT"
SKATO<-"SKATO"
SKAT_125<-"SKAT_125"
SKAT_C <- "SKAT_C"
### START loop over list of genes
cat(paste0("\nStart loop over list of genes in chr ",chr,"in geneset driectory",dir,"\n"))
datafiles<-Sys.glob(paste(dir,"/*",extension,sep=""))
for (i in 1:length(datafiles)) {
  print(i)
  #print(datafiles[i])
  genename=basename(file_path_sans_ext(datafiles[i]))
  print(genename)
  
  ### GET GENO DATA
  SET<-read.table(datafiles[i],head=F, as.is=T, check.names=F,sep="\n")
  SETlist<-c(SET$V1)
  SETlist <- SETlist[SETlist %in% colnames(com.data)]
  if (length(SETlist) == 0) {
   next
  }
  snpID<-SET[1,1]
  SETgeno<-as.matrix(com.data[,SETlist])
  
  ### Define y
  y=as.matrix(com.data$STATUS)
  table(y)
  
  # MODEL adjusted by SEX and fisrt three PCs
  covariates = c("SEX.y","PC1","PC2","PC3")
  index= which(colnames(com.data)%in%covariates)
  Xm2=as.matrix(com.data[,index])
  obj<-SKAT_Null_Model(y ~ Xm2,out_type="D")
  n0<-SKAT(SETgeno, obj, kernel = "linear.weighted")$param$n.marker	# n input SNPs
  n1<-SKAT(SETgeno, obj, kernel = "linear.weighted")$param$n.marker.test # n ouptut SNPs
  skat<-SKAT(SETgeno, obj, kernel = "linear.weighted")$p.value
  skato<-SKAT(SETgeno, obj, kernel = "linear.weighted", method="optimal.adj")$p.value
  skat_125<-SKAT(SETgeno,obj, kernel="linear.weighted", weights.beta=c(1,25))$p.value
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3675243/
  # SKAT-C	SC	combined sum test with SKAT tests for rare and common variants
  skat_c<-SKAT_CommonRare(SETgeno,obj, method="C")$p.value
  
  ## add values to lists
  genelist<-rbind(genelist,genename)
  SNPID<-rbind(SNPID,snpID)
  N0<-rbind(N0,n0)
  N1<-rbind(N1,n1)
  
  SKAT<-rbind(SKATO,skat)
  SKATO<-rbind(SKATO,skato)
  SKAT_125<-rbind(SKAT_125,skat_125)
  SKAT_C <- rbind(SKAT_C, skat_c)
}

write.table(cbind(genelist,SNPID,N0,N1,SKAT,SKATO,SKAT_125, SKAT_C), paste(dir,"_chr",chr,"_SKAT.pval",sep=""), sep="\t", col.names=F, row.names=F, quote=F)
