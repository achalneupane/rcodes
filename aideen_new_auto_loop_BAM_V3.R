


library(GenomicFeatures)
##library(multicore)
library(Rsamtools)
library(BSgenome)
library(parallel)

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
code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")

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
## genes.file <- "/media/scratch/Aideens_program/genescdp.txt"
## Project.data<-"/media/UQCCG/Sequencing/Projects/SKDP"
## Sample.name <- c("NSGC-9.3","SKDP-55.3","SKDP-6.1") # "ALL" use "ALL" is want to get all in the /media/UQCCG/Sequencing/Projects/SKDP directory!
## #Sample.name <- c("NSGC-9.3_C4GD4ACXX-8-15","SKDP-55.3_C0J43ACXX-2-16","SKDP-6.1")
## conversion.file <- "/media/scratch/Aideens_program/GeneNameConversion.csv"
## out.file<- "/media/scratch/Aideens_program/out.txt"

## ### AIDEEN TO CHANGE THE PATH
## genes.file <- "/media/UQCCG/Sequencing/Projects/Aglient_clinical_exomes/List of genes for testing agilent.csv"
## Project.data<-"/media/UQCCG/Sequencing/Projects/Aglient_clinical_exomes"
## Sample.name <- "ALL" # "ALL" use "ALL" is want to get all in the /media/UQCCG/Sequencing/Projects/SKDP directory!
## #Sample.name <- c("NSGC-9.3_C4GD4ACXX-8-15","SKDP-55.3_C0J43ACXX-2-16","SKDP-6.1")
## conversion.file <- "/media/UQCCG/Sequencing/Data/Genomes/GeneNameConversion.csv"
## out.file<- "/media/UQCCG/Sequencing/Projects/Aglient_clinical_exomes/out.txt"







## ### AIDEEN TO CHANGE THE PATH
## ## [1] "Missing genes  CYP26 ,  check HUGO gene symbol"  CYP26A1
## ## [2] "Missing genes  EPH4A ,  check HUGO gene symbol"  EPHA4
## ## [3] "Missing genes  PARAXIS ,  check HUGO gene symbol" TCF15 
## ## [4] "Missing genes  RALDH2 ,  check HUGO gene symbol" ALDH1A2
## ## [5] "Missing genes  SNAIL1 ,  check HUGO gene symbol" SNAI1
## ## [6] "Missing genes  SNAIL2 ,  check HUGO gene symbol" SNAI2

#ann<-read.delim("/media/scratch/Aideens_program/Copy of Duncan and Sally's samples for coverage analysis v2.csv", header=T,sep="\t",skip=0,fill=TRUE,stringsAsFactors=FALSE)
ann<-read.delim("/media/scratch/Aideens_program/Copy of Duncan and Sally's samples for coverage analysis v2.csv", header=T,sep="\t",skip=0,fill=TRUE,stringsAsFactors=FALSE)
/media/scratch/Aideens_program/duncans.extended.list.csv
table(ann$Project)
## SDDS SKDP 
##   90   74
the.project<-"SDDS"
the.project<-"SKDP"

## genes.file <- "/media/scratch/Aideens_program/Duncan's candidate gene list.csv"
genes.file <- "/media/scratch/Aideens_program/duncans.extended.list.csv"

Project.data<-paste("/media/UQCCG/Sequencing/Projects",the.project,sep="/") 
Sample.name <- ann[ann$Project==the.project,"Sample"] # "ALL" # "ALL" use "ALL" is want to get all in the /media/UQCCG/Sequencing/Projects/SKDP directory!
#Sample.name <- c("NSGC-9.3_C4GD4ACXX-8-15","SKDP-55.3_C0J43ACXX-2-16","SKDP-6.1")
conversion.file <- "/media/UQCCG/Sequencing/Data/Genomes/GeneNameConversion.csv"
out.file<- "/media/scratch/Aideens_program/SDDS_ex.txt"
out.file<- "/media/scratch/Aideens_program/SKDP.txt"


##############################################################################################
## ## ### AIDEEN TO CHANGE THE PATH
## genes.file <- "/media/scratch/Aideens_program/aml.txt.csv"
## Project.data<-"/media/UQCCG/Sequencing/Projects/TGCM-AML"
## Sample.name <- "ALL" # "ALL" use "ALL" is want to get all in the /media/UQCCG/Sequencing/Projects/SKDP directory!
## #Sample.name <- c("NSGC-9.3_C4GD4ACXX-8-15","SKDP-55.3_C0J43ACXX-2-16","SKDP-6.1")
## conversion.file <- "/media/UQCCG/Sequencing/Data/Genomes/GeneNameConversion.csv"
## out.file<- "/media/UQCCG/Sequencing/Projects/Aglient_clinical_exomes/out.txt"
##############################################################################################



path.type<-"/"  # "\\" use for windows!!!

## Preset parameters
the.genome<-"hg19" 
genomes.path<-"/mnt/UQCCG/Sequencing/Data/Genomes/hg19/" ## path to genomes files



force.redo.with.new.ranges<-FALSE # force.redo.with.new.ranges<-TRUE # set true is want to recalculate : IF FALSE it will over-write out of date files: see all.QC.column.labels
extend.exons<-0 # all to exons start/end to get additional coverage -200 
##### need to choose one of the below:

######################################## no need to change below ###################
all.QC.column.labels<-c("ID","Sample","Recipe","Capture.Method","Lane","Run","target_file","total_reads","total_mapped_reads","total_dup_reads","unmapped_reads","total.bases",paste("on.target.bases.",extend.exons,"bp",sep=""),"on.target.bases", "percent.on.target.bases","median.max.coverage", "median.mean.coverage","mean.mean.coverage","percent.ccds.gt.1","percent.ccds.gt.5","percent.ccds.gt.10","percent.ccds.gt.15","percent.ccds.gt.30","percent_Duplicated","percent_Unmapped","Description")



##################################################
# STEP 1 - READ IN THE GENES LIST
##################################################

#-----------------------
## Load in the genes
#-----------------------

## Read in the genes list
genes<- {}
genes<-read.delim(genes.file, header=F,sep="\t",skip=0,fill=TRUE,stringsAsFactors=FALSE)

### Make sure no spaces - bad spaces!
genes[,1] <- gsub("\\s+","",genes[,1])
genes <- t(genes)
genes <- toupper(genes)
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

## genes
##   KMT2D:->"MLL2"
## NEURL1->NEURL

print(paste("Missing genes ",paste(genes[!(genes %in% conversion[,"Approved Symbol"])],collpase=", "), "check HUGO gene symbol"))

conversion[1:5,]
grepl(genes[!(genes %in% conversion[,"Approved Symbol"])],conversion[,"Approved Symbol"])
conversion[grep("^NEU$",conversion[,"Approved Symbol"]),]
conversion[grep("^9148$",conversion[,2])]
9148,

### error => need to check that the gene exists
genes.hugo <- conversion[conversion[,"Approved Symbol"] %in% genes,"Entrez Gene ID (mapped data supplied by NCBI)"]



  cds <- {}

  columns(txdb) # what you can get
  cds <- cds(txdb, vals=list(gene_id=genes.hugo), columns=c( "gene_id",  "tx_name", "GENEID","CDSID","EXONRANK","EXONID"))
  cds
 
  #####
   ## seq <- as.vector(seqnames(cds))
   ## seq
  ## start <- as.vector(start(cds))
  ## end <- end(cds)
  ## tx.name<-unlist(lapply(cds$"tx_name",function(x) paste(x,collapse=";")))
  ## gene.id<-unlist(lapply(cds$"GENEID",function(x) paste(x,collapse=";")))
  ## exonRank<-unlist(lapply(cds$"EXONRANK",function(x) paste(x,collapse=";")))
  ## exonID<-unlist(lapply(cds$"EXONID",function(x) paste(x,collapse=";")))
  
  
  rl.from.gr<-as(cds, "RangesList") # data.gr contains
 # rl.from.gr<-reduce(rl.from.gr) # use to collapse regions
 # names(rl.from.gr)
empty.chromo<-lapply(rl.from.gr,length)==0
rl.from.gr<-rl.from.gr[!empty.chromo]


##################################################
# STEP 2 - PROCESS OUR BAM FILE
##################################################


options(show.error.messages = TRUE)


##num.cores<-1 # number of CPU coores to use
##options(cores=num.cores )
## I had to install doMC ## install packages("doMC")
##library(doMC)
##registerDoMC(cores=num.cores)



#############




# where we find the QC dir
QC.root.directory<-paste(Project.data,"QC",sep=path.type)
xx<-try(setwd( QC.root.directory  ),silent=TRUE)
if(inherits(xx, "try-error")){system(paste("mkdir",QC.root.directory,sep=" "));setwd( QC.root.directory  )}

# where we find the QC files
QC.directory<-paste(Project.data,"QC","AlignedQC",sep=path.type)
xx<-try(setwd( QC.directory  ),silent=TRUE)
if(inherits(xx, "try-error")){system(paste("mkdir",QC.directory,sep=" "));setwd( QC.directory);QC.files<-dir(getwd())}else{QC.files<-dir(getwd())}
QC.directory

## go to the bam directory
setwd(QC.directory)

## get all the files
files<-dir(getwd())

## grab only the *.bam files
sample.bams<-files[grepl(".regions.RData$",files)] # sample.bams<-sample.bams.ori[48:132] # sample.bams<-sample.bams.ori[1:47] # sample.bams.ori<-sample.bams

## Now get the bam for your Sample.name that you set above eg "SKDP178.2"
grep("SKDP-36.3",sample.bams)


if(sum(Sample.name!="ALL")==length(Sample.name)){

use<-{}  
wanted<-rep(FALSE,length(sample.bams))
irun<-54
for(irun in 1:length(Sample.name)){

possible<-sample.bams[grep(paste(Sample.name[irun],"_",sep=""),sample.bams,fixed=TRUE)]

if(length(possible)==0){
possible<-sample.bams[grep(Sample.name[irun],sample.bams,fixed=TRUE)]
}
  
best<-possible[grepl("combined",possible)]
length(best)
if(length(best)==0 & length(possible)>0){best<-possible[1]}
print(paste(irun,Sample.name[irun],"From possible:",paste(possible,collapse=", "), "-> using -> ", best))
  
use<-c(use,best)
 }  # grep("21",sample.bams)
sample.bams<-use
}

######################################################3

if (length(sample.bams) == 0) {
  print(paste("Cant find region  file for",Sample.name,sep=" "))
 }else{
print(paste("Doing samples: ",sample.bams))
}


##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
########################## loop over BAMS
## grab the bam header
all.counts<-{}
all.samples<-{}
ib<-2
sample.bams # sample.bams<-sample.bams[c(1,4)]
for(ib in 1:length(sample.bams)){
#sample.bam<-sample.bams[ib]
sample.bam<-paste(QC.directory, sample.bams[ib],sep=path.type)
  print(paste("Doing:",sample.bam))
## the.header<-scanBamHeader(sample.bam) # list with names
## the.bam.files<-names(the.header)
sample.bam.counts<-{}
total.reads<-{}
total.coverage<-{}





####################### Load in the re-calculated coverage information (derived form the bam file previosly #####
print(QC.directory)
pre.loaded.regions <-  sample.bam
pre.loaded.regions

  if( !file.exists(pre.loaded.regions) ) {
    print(paste("File:", pre.loaded.regions,"missing"))
     next
  }

load(pre.loaded.regions)
#cov ### object cov is loaded from above file
rownames(data)

## Sample name extracted from bam file
the.sample<-as.character(data["Sample",])


the.current.bam<-sample.bam

## get run id
a.run<- as.character(data["Run",])
## extract.value.from.DS(the.readgroup["DS",],match.string="SampleRef:") # genome ref
a.lane<-as.character(data["Lane",])# LAne extract.value.from.DS(the.readgroup["DS",ii],match.string="Lane:")

# what capture was used
a.capture<- as.character(data["Capture.Method",]) # capture.tech extract.value.from.DS(the.readgroup["DS",ii],match.string="Description:")
a.recipe<- as.character(data["Recipe",]) # recipe extract.value.from.DS(the.readgroup["DS",ii],match.string="Recipe:") 
a.sample<-the.sample # extract.value.from.DS(the.readgroup["DS",ii],match.string="ParticipantCode:") 
a.sample.ID<-as.character(data["ID",])  # ID extract.value.from.DS(the.readgroup["DS",ii],match.string="SampleID:") 
a.DS<-as.character(data["Description",]) # a.DS<-the.readgroup["DS",ii]

## write.table(paste(the.sample,a.capture,sep=","),file=out.file,col.names=FALSE,row.names=FALSE,sep=",",quote=FALSE,append=FALSE)




######################################################################################################
####### Convert the regions to format to be used with existing coverage code
######################################################################################################

  #targets.file.ori<-targets.file
#rl.from.gr<-rl.from.gr[names(cov)] ## must have the same chromom names in both files
### AIDEEN TO CHANGE THE PATH
setwd(code.dir)
source("auto_loop_over_UQCCG_coverage_with_cov_genes.r")  ## This ones uses a max of 8 cores 

the.counts[1:5,]

########################## WRITE OUTPUT 
 ## "chr,start,end,length,Sum,Max,Mean,bp.gt.1,bp.gt.5,bp.gt.10,bp.gt.15,bp.gt.30,percent.ccds.gt.1,percent.ccds.gt.5,percent.ccds.gt.10,percent.ccds.gt.15,percent.ccds.gt.30"
## paste(colnames(the.counts),collapse=",")
output.labels<-c("Sum","Max","Mean","bp.gt.1","bp.gt.5","bp.gt.10","bp.gt.15","bp.gt.30","percent.ccds.gt.1","percent.ccds.gt.5","percent.ccds.gt.10","percent.ccds.gt.15","percent.ccds.gt.30")
the.sample
a.capture
new.labels<-paste(a.sample,output.labels,sep=":")
colnames(the.counts)[colnames(the.counts) %in% output.labels]<-new.labels

if(is.null(dim(all.counts))){
all.counts<-the.counts
all.samples<-the.sample
}else{
  all.counts<-cbind(all.counts,the.counts[,new.labels])
  all.samples<-c(all.samples,the.sample)

}


}## loop over coverage files (samples)


order.wanted<- expand.labels.to.samples.complex(all.samples,output.labels,paste.after=TRUE,seperator=":")
order.wanted<-c("chr","start","end",order.wanted)

  #####
   chr <- as.vector(seqnames(cds))
   start <- as.vector(start(cds))
  end <- end(cds)
  tx.name<-unlist(lapply(cds$"tx_name",function(x) paste(x,collapse=";")))
  gene.id<-unlist(lapply(cds$"GENEID",function(x) paste(x,collapse=";")))
  exonRank<-unlist(lapply(cds$"EXONRANK",function(x) paste(x,collapse=";")))
   exonID<-unlist(lapply(cds$"EXONID",function(x) paste(x,collapse=";")))



the.loci<-cbind(chr,start,end,gene.id,exonRank,exonID,tx.name)
key.counts<-build.key(all.counts,c("chr","start","end"))
key.loci<-build.key(the.loci,c("chr","start","end"))

posns<-match(key.counts,key.loci)
sum(is.na(posns))
the.loci<-the.loci[posns,]

dim(the.loci)
dim(all.counts)

the.loci[1:5,]
all.counts[1:5,1:5]


Gene.Name<-match.cols.and.collapse.table(the.loci,"gene.id",";",conversion,"Entrez Gene ID (mapped data supplied by NCBI)","Approved Symbol","delimit")
# match.cols.and.collapse.tabletable1,table1.match.col,table.delim,tab2,tab2.match.col,tab2.collapse.cols,gene.collapse)



cbind(the.loci,Gene.Name,all.counts[,order.wanted])[1:5,1:20]

write.table( cbind(the.loci,Gene.Name,all.counts[,order.wanted]),file=out.file,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE,append=FALSE)




