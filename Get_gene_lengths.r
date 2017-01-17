



####################################################
########### END OF FUNCTIONS #############
############# ############# ############# #############
##### ---- START HERE #####################
############# ############# ############# #############
## setwd(BAM.directory)
## files<-dir(getwd())
## sample.bams<-files[grepl(".bam$",files)] # sample.bams<-sample.bams.ori[48:132] # sample.bams<-sample.bams.ori[1:47] # sample.bams.ori<-sample.bams
## sample.bams


# PARAMETERS TO MODIFY 
### AIDEEN TO CHANGE THE PATH
#genes.file <- "/media/scratch/Aideens_program/common.genes.txt.csv"
genes.file <- "/media/scratch/Aideens_program/aml.txt.csv"
genes.file <- "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/clusters_as.using.txt"

#Sample.name <- c("AMAS-7.3")
conversion.file <- "/media/scratch/Aideens_program/GeneNameConversion.csv"
out.file<- "/media/UQCCG/Sequencing/Projects/Aglient_clinical_exomes/Analysis/aml_all.txt"

path.type<-"/"  # "\\" use 

## Preset parameters
the.genome<-"hg19" 
genomes.path<-"/mnt/UQCCG/Sequencing/Data/Genomes/hg19/" ## path to genomes files
code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
source(paste(code.dir,"annotate_SNPs_subroutines.r",sep=path.type))


force.redo.with.new.ranges<-FALSE # force.redo.with.new.ranges<-TRUE # set true is want to recalculate : IF FALSE it will over-write out of date files: see all.QC.column.labels
extend.exons<-0 # all to exons start/end to get additional coverage -200 
##### need to choose one of the below:

######################################## no need to change below ###################
all.QC.column.labels<-c("ID","Sample","Recipe","Capture.Method","Lane","Run","target_file","total_reads","total_mapped_reads","total_dup_reads","unmapped_reads","total.bases",paste("on.target.bases.",extend.exons,"bp",sep=""),"on.target.bases", "percent.on.target.bases","median.max.coverage", "median.mean.coverage","mean.mean.coverage","percent.ccds.gt.1","percent.ccds.gt.5","percent.ccds.gt.10","percent.ccds.gt.15","percent.ccds.gt.30","percent_Duplicated","percent_Unmapped","Description")






library(SKAT)
library(GenomicFeatures)
##library(multicore)
library(Rsamtools)
library(BSgenome)
library(multicore)

package<-"BSgenome.Hsapiens.UCSC.hg19"
library(package,character.only=TRUE)
the.chroms<-seqlengths(Hsapiens)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
                                        #
## setRepositories()
## install.packages("TxDb.Hsapiens.UCSC.hg19.knownGene")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

########### just select and run

 




basesAboveThreshold<-function(x,thresh=0){sum(x>=thresh)}

##################################################
# STEP 1 - READ IN THE GENES LIST
##################################################

#-----------------------
## Load in the genes
#-----------------------

## Read in the genes list
genes<- {}
genes.lists<-read.delim(genes.file, header=T,sep="\t",skip=0,fill=TRUE,stringsAsFactors=FALSE)

genes.lists[1:5,]

genes<-genes.lists[,"FANC_complex.all"]
genes<-genes.lists[,"Checkpoint.Proteins"]

### Make sure no spaces - bad spaces!
genes <- gsub("\\s+","",genes)
genes <- t(genes)
genes <- toupper(genes)
genes<- genes[genes!=""]
genes

## or you can override genes and use a list
##genes <- c("COL1A1","COL1A2","FBN1")



#-----------------------
## Convert the gene names to unique gene_id
#-----------------------
options(max.print=300)
conversion <- read.delim(conversion.file, header=F,sep=",",skip=0,fill=TRUE,stringsAsFactors=FALSE)
conversion <- t(conversion)
conversion <- t(conversion)
dim(conversion)
conversion[1:5,1:2]
colnames(conversion)<-conversion[1,]
conversion<-conversion[-1,]
conversion[1:5,]




print(paste("Missing genes ",paste(genes[!(genes %in% conversion[,"Approved Symbol"])],collpase=", "), "check HUGO gene symbol"))
### error => need to check that the gene exists
genes.hugo <- conversion[conversion[,"Approved Symbol"] %in% genes,"Entrez Gene ID (mapped data supplied by NCBI)"]



  cds <- {}

  columns(txdb) # what you can get
  cds <- cds(txdb, vals=list(gene_id=genes.hugo), columns=c( "gene_id",  "tx_name", "GENEID","CDSID","EXONRANK","EXONID"))
  cds

seqlevels(cds)
length(cds)
length(width(cds))

total.length<-sum(width(cds))
total.length
  
rl.from.gr<-as(cds, "RangesList") # data.gr contains
cds<-reduce(cds)
 # rl.from.gr<-reduce(rl.from.gr) # use to collapse regions
 # names(rl.from.gr)
empty.chromo<-lapply(rl.from.gr,length)==0
rl.from.gr<-rl.from.gr[!empty.chromo]



write.table( cbind(the.loci,Gene.Name,all.counts[,order.wanted]),file=out.file,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE,append=FALSE)


out.file <- "/media/scratch/Aideens_program/Russel.txt"

library("SKAT")

data(SKAT.haplotypes)
names(SKAT.haplotypes)
attach(SKAT.haplotypes)


clinical<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/Clinical.GENOTYPE.conponents:.SkatO.clusters.coding.0.001.all.geno.all.filters_no.imput_paper.txt", header=T,sep="\t",skip=0,fill=TRUE,stringsAsFactors=FALSE)
clinical[1:5,1:10]

dim(clinical)
sum(clinical[,"p"]<=0.05)
17/153

set.seed(500)
# out.c<-Power_Continuous(Haplotype,SNPInfo$CHROM_POS, SubRegion.Length=total.length,Causal.Percent= 20, N.Sim=10, MaxBeta=2,Negative.Percent=20)
[1] "10/10"
 out.c

total.length

out.b<-Power_Logistic(Haplotype,SNPInfo$CHROM_POS, SubRegion.Length=total.length,Prevalence=0.005,Causal.Percent= 10,Causal.MAF.Cutoff=0.01, N.Sim=1 ,MaxOR=7,N.Sample.ALL=c(200), Negative.Percent=0)
out.b
out.b<-Power_Logistic(Haplotype,SNPInfo$CHROM_POS, SubRegion.Length=total.length,Prevalence=0.005,Causal.Percent= 10,Causal.MAF.Cutoff=0.01, N.Sim=1 ,MaxOR=7,N.Sample.ALL=c(400), Negative.Percent=0)
out.b

out.b<-Power_Logistic(Haplotype,SNPInfo$CHROM_POS, SubRegion.Length=157000,Prevalence=0.005,Causal.Percent= 15,Causal.MAF.Cutoff=0.001, N.Sim=10 ,MaxOR=7,N.Sample.ALL=c(200), Negative.Percent=0)
out.b
> $Power
    0.01     0.001       1e-06
200    1 0.6254334 0.002003634


out.b<-Power_Logistic(Haplotype,SNPInfo$CHROM_POS, SubRegion.Length=157000,Prevalence=0.005,Causal.Percent= 15,Causal.MAF.Cutoff=0.001, N.Sim=50 ,MaxOR=7,N.Sample.ALL=c(200), Negative.Percent=0)
out.b



out.b<-Power_Logistic(Haplotype,SNPInfo$CHROM_POS, SubRegion.Length=157000,Prevalence=0.005,Causal.Percent= 10,Causal.MAF.Cutoff=0.001, N.Sim=50 ,MaxOR=7,N.Sample.ALL=c(200), Negative.Percent=0)
out.b


$Power
         0.01   0.001      1e-06
200 0.9027881 0.66911 0.07783071


out.b<-Power_Logistic(Haplotype,SNPInfo$CHROM_POS, SubRegion.Length=50000,Prevalence=0.005,Causal.Percent= 10,Causal.MAF.Cutoff=0.001, N.Sim=1 ,MaxOR=7,N.Sample.ALL=c(400), Negative.Percent=0)
out.b



out.b<-Power_Logistic(Haplotype,SNPInfo$CHROM_POS, SubRegion.Length=50000,Prevalence=0.005,Causal.Percent= 15,Causal.MAF.Cutoff=0.001, N.Sim=1 ,MaxOR=7,N.Sample.ALL=c(200), Negative.Percent=0)
out.b


$Power
         0.01   0.001      1e-06
200 0.9027881 0.66911 0.07783071




         0.01       0.001        1e-06
200 0.02228399 0.001993131 1.284074e-06

        0.01      0.001        1e-06
200 0.1371193 0.03030285 0.0003588669

## > out.b<-Power_Logistic(Haplotype,SNPInfo$CHROM_POS, SubRegion.Length=157000,Prevalence=0.005,Causal.Percent= 10,Causal.MAF.Cutoff=0.001, N.Sim=10 ,MaxOR=7,N.Sample.ALL=c(200), Negative.Percent=0)
##          0.01    0.001     1e-06
## 200 0.9419695 0.800915 0.3618237



##  > out.b<-Power_Logistic(Haplotype,SNPInfo$CHROM_POS, SubRegion.Length=157000,Prevalence=0.005,Causal.Percent= 10,Causal.MAF.Cutoff=0.001, N.Sim=50 ,MaxOR=7,N.Sample.ALL=c(200), Negative.Percent=0)
## out.b
## [1] "10/50"
## [1] "20/50"
## [1] "30/50"
## [1] "40/50"
## [1] "50/50"
## > $Power
##          0.01     0.001     1e-06
## 200 0.9102846 0.8261971 0.3698216


## > out.b<-Power_Logistic(Haplotype,SNPInfo$CHROM_POS, SubRegion.Length=157000,Prevalence=0.005,Causal.Percent= 10,Causal.MAF.Cutoff=0.001, N.Sim=50 ,MaxOR=7,N.Sample.ALL=c(200), Negative.Percent=0)
## out.b
## [1] "10/50"
## [1] "20/50"
## [1] "30/50"
## [1] "40/50"
## [1] "50/50"
## > $Power
##          0.01     0.001     1e-06
## 200 0.8976724 0.7950395 0.3556721


##  out.b<-Power_Logistic(Haplotype,SNPInfo$CHROM_POS, SubRegion.Length=157000,Prevalence=0.005,Causal.Percent= 15,Causal.MAF.Cutoff=0.001, N.Sim=10 ,MaxOR=7,N.Sample.ALL=c(200), Negative.Percent=0)
## out.b
## [1] "10/10"
## > $Power
##     0.01 0.001     1e-06
## 200    1     1 0.9113089




## > out.b<-Power_Logistic(Haplotype,SNPInfo$CHROM_POS, SubRegion.Length=157000,Prevalence=0.005,Causal.Percent= 15,Causal.MAF.Cutoff=0.001, N.Sim=50 ,MaxOR=7,N.Sample.ALL=c(200), Negative.Percent=0)
## out.b
## [1] "10/50"
## [1] "20/50"
## [1] "30/50"
## [1] "40/50"
## [1] "50/50"
## > $Power
##          0.01     0.001     1e-06
## 200 0.9916834 0.9474374 0.7033643
