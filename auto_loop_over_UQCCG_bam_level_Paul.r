
## ##################################### overlap for introns/exons for NUMBER OF READS not bp coverage
## estimated library size
## "Based on the Lander-Waterman equation that states: C/X = 1 - exp( -N/X ) where X = number of distinct molecules in library N = number of read pairs C = number of distinct fragments observed in read pairs"###################  #START ###############


#########################
## in case of multiple RG's in a BAM file this assumes that the same target file is used for each RG
## The coverage calculation is perfoemed overe the headers in the BAm file which shoudl be the same as the Bam file name.
## The coverge calculation IS NOT  broken d
##

###############################libraries needed
library(GenomicFeatures)
library(multicore)
library(Rsamtools)
library(BSgenome)
library(foreach)  ## could not get this to work well
library(doMC)


##   paramAll <- ScanBamParam(scanBamFlag(isUnmappedQuery=NA, isDuplicate=NA, isPaired=NA),what=c("rname"))
## all.reads<-countBam(the.bam.files[ii],param=paramAll) # same as  (total.reads+total.unmapped.reads+total.dup.reads)
##############the reads contain all the sequence reads 
#total.reads<-length(the.reads)

  countTotalReadsInBam <- function(seqname, bamFile, human.chromlens, ...){
          param <-ScanBamParam( which=GRanges(seqname, IRanges(1, human.chromlens[seqname])),flag=scanBamFlag(isUnmappedQuery=NA, isDuplicate=NA, isPaired=NA),what=c("qname"))
    count.reads <- countBam(bamFile,param=param)
    count.reads
}

  countUnmappedReadsInBam <- function(seqname, bamFile, human.chromlens, ...){
    param <- ScanBamParam( which=GRanges(seqname, IRanges(1, human.chromlens[seqname])),flag=scanBamFlag(isUnmappedQuery=TRUE),what=c("qname"))
    count.reads <- countBam(bamFile,param=param)
    count.reads
}

    countDupReadsInBam <- function(seqname, bamFile, human.chromlens, ...){
    param <- ScanBamParam( which=GRanges(seqname, IRanges(1, human.chromlens[seqname])),flag=scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=TRUE, isPaired=NA),what=c("rname"))
    count.reads <- countBam(bamFile,param=param)
    count.reads
}

  get.list.elements<-function(x,column){
    x[column]}


extractReadsFromBam <- function(seqname, bamFile, human.chromlens, ...){
    param <- ScanBamParam( which=GRanges(seqname, IRanges(1, human.chromlens[seqname])),flag=scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=FALSE, isPaired=NA))
    the.reads <-readBamGappedAlignments(bamFile,param=param)
    the.reads
}

extractCovFromBam <- function(seqname, bamFile, human.chromlens, ...){
    param <- ScanBamParam( which=GRanges(seqname, IRanges(1, human.chromlens[seqname])),flag=scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=FALSE, isPaired=NA))
    the.reads <-readBamGappedAlignments(bamFile,param=param)
    cov<-coverage(the.reads,width=human.chromlens) 
    cov
}

extractCovFromBamOneChr<- function(seqname, bamFile, human.chromlens, ...){
    param <- ScanBamParam( which=GRanges(seqname, IRanges(1, human.chromlens[seqname])),flag=scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=FALSE, isPaired=NA))
    the.reads <-readBamGappedAlignments(bamFile,param=param)
    cov<-coverage(the.reads,width=human.chromlens) 
    cov[seqname]
}


extractReadsFromBAMwithCigar <- function(file) # Get reads but take cigar into account
{
  ## This ScanBamParam object allows us to load only the necessary
  ## information from the file.
  param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE,
                                         isDuplicate=NA),
                        what=c("rname", "pos", "cigar"))
  
  bam <- scanBam(file, param=param)[[1]]
  ## Note that unmapped reads and reads that are PCR/optical duplicates
  ## have already been filtered out by using the ScanBamParam object above.
  irl <- cigarToIRangesListByRName(bam$cigar, bam$rname, bam$pos)
  irl <- irl[elementLengths(irl) != 0]  # drop empty elements
  irl
}

extractUnmappedReadsFromBAM <- function(file) # Get reads but take cigar into account
{
  ## This ScanBamParam object allows us to load only the necessary
  ## information from the file.
  param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=TRUE), what=c("qname"))
  
  bam <- scanBam(file, param=param)[[1]]
  bam #bam$qname[1:5]
}

extract.value.from.info<-function(info,match.string){ #will be fooled if have GGMAF=  #ok if GMAF not present need [^;]
match.string.general<-paste(match.string,"[a-zA-Z0-9,\\.]*",sep="") #regexec(";GMAF=[a-zA-Z0-9,\\.]*;",indels[has.gmaf,"INFO"]) ";" not needed
match<- regexec(match.string.general,info)
start<-unlist(match)
length<-unlist(lapply(match,function(x)  attr(x,"match.length")))
end<-start+length-1 # minus one to get last position
start<-start+nchar(as.character(match.string))
substr(info,start=start,stop=end)
}

extract.value.from.format<-function(info,match.string){ #will be fooled if have GGMAF=  #ok if GMAF not present need [^;]
match.string.general<-paste(match.string,"[a-zA-Z0-9\\ \\.]*",sep="") #regexec(";GMAF=[a-zA-Z0-9,\\.]*;",indels[has.gmaf,"INFO"]) ";" not needed
match<- regexec(match.string.general,info)
start<-unlist(match)
length<-unlist(lapply(match,function(x)  attr(x,"match.length")))
end<-start+length-1 # minus one to get last position
start<-start+nchar(as.character(match.string))
substr(info,start=start,stop=end)
}

extract.value.from.DS<-function(info,match.string){ #will be fooled if have GGMAF=  #ok if GMAF not present need [^;]
match.string.general<-paste(match.string,"[a-zA-Z0-9,:\\+\\ \\.\\-]*",sep="")  #regexec(";GMAF=[a-zA-Z0-9,\\.]*;",indels[has.gmaf,"INFO"]) ";" not needed
match<- regexec(match.string.general,info)
start<-unlist(match)
length<-unlist(lapply(match,function(x)  attr(x,"match.length")))
end<-start+length-1 # minus one to get last position
start<-start+nchar(as.character(match.string))
substr(info,start=start,stop=end)
}
## the.readgroup["DS",]
## extract.value.from.DS(the.readgroup["DS",],match.string="SampleRef:")
## extract.value.from.DS(the.readgroup["DS",],match.string="Lane:")
## extract.value.from.DS(the.readgroup["DS",],match.string="Description:")
## extract.value.from.DS(the.readgroup["DS",],match.string="Recipe:")
## extract.value.from.DS(the.readgroup["DS",],match.string="ParticipantCode:")
## extract.value.from.DS(the.readgroup["DS",],match.string="SampleID:")


#non.project.directories<-c("Genomes","R-codes","annovar","Sequence_Reads","scripts","bcos_srv"
                         

get.readgroup.info<-function(sample.bam,readgroup.labels=c("ID","PL","LB","DS","SM")){ #returms a matrix of readgroups
  the.header<-scanBamHeader(sample.bam) # list with names is sample.bams is a vertor then it does all of them at once!
  the.bam.files<-names(the.header)
  a.readgroup<-matrix(data=NA,nrow=length(readgroup.labels),ncol=1)
  rownames(a.readgroup)<-readgroup.labels
  readgroup.cols<- {}

  i.ar<-1
  for(ib in 1:length(the.bam.files)){  ## perhaps different BAm files
  ## the.readgroup<-the.header[[the.bam.files[ib]]][["text"]][["@RG"]] # this only provides the FIRST RG
  the.readgroup.list<-the.header[[the.bam.files[ib]]][["text"]][names(the.header[[the.bam.files[ib]]][["text"]])=="@RG"] ## perhaps many RG per BAm file

  for(iread in 1:length(the.readgroup.list)){  # loop over RG
  the.readgroup<-the.readgroup.list[[iread]]

  the.RG.desc<-unlist(lapply(strsplit(the.readgroup,split=":"),function(x) x[1]))

  for(ir in 1:length(the.readgroup)){the.readgroup[ir]<-gsub(paste("^",the.RG.desc[ir],":",sep=""),"",the.readgroup[ir])} #get rid of label from text
  names(the.readgroup)<-the.RG.desc
  
  if(i.ar==1){a.readgroup[readgroup.labels,1]<-the.readgroup[readgroup.labels];readgroup.cols<-the.bam.files[ib]}else{
    a.readgroup<-cbind(a.readgroup[readgroup.labels,],the.readgroup[readgroup.labels])
    readgroup.cols<-c(readgroup.cols,the.bam.files[ib])
  }

  ############### fix the capture technology labels ";" delinates label seperators 
   if(is.na(a.readgroup["DS",i.ar])){a.readgroup["DS",i.ar]<-"NA"}  # so can just match with text later
   a.readgroup["DS",i.ar]<-gsub("TruD;TruX","TruD:TruX",a.readgroup["DS",i.ar]) # special fix
   a.readgroup["DS",i.ar]<-gsub("TruD;NimX3","TruD:NimX3",a.readgroup["DS",i.ar]) # special fix
   a.readgroup["DS",i.ar]<-gsub("KapD;NimX3","KapD:NimX3",a.readgroup["DS",i.ar]) # special fix

a.readgroup["DS",i.ar]<-gsub("NxtD;NxtX","NxtD:NxtX",a.readgroup["DS",i.ar]) # special fix
a.readgroup["DS",i.ar]<-gsub("NxtD;NxtXR","NxtD:NxtXR",a.readgroup["DS",i.ar]) # special fix
a.readgroup["DS",i.ar]<-gsub("NxtD;NxtXE","NxtD:NxtXE",a.readgroup["DS",i.ar]) # special fix


  a.readgroup["DS",i.ar]<-gsub("NEBuD;TruX","NEBuD;TruX",a.readgroup["DS",i.ar]) # special fix
   a.readgroup["DS",i.ar]<-gsub("TruD;Agl1.2","TruD:Agl1.2",a.readgroup["DS",i.ar]) # special fix
   a.readgroup["DS",i.ar]<-gsub("TruD;NimX","TruD:NimX",a.readgroup["DS",i.ar]) # special fix  ;Description:Exome;
a.readgroup["DS",i.ar]<-gsub("TruD;NimX3","TruD:NimX3",a.readgroup["DS",i.ar]) # special fix  ;Description:Exome;
   a.readgroup["DS",i.ar]<-gsub(";Description:Exome;",";Description:TruD:NimX;",a.readgroup["DS",i.ar]) # special fix
   ## "NA" nimblegen v1 - "TruD:NimX" nimblegen v2 - "TruD:TruX" illumina v2 - "TruD;" whole genome
	print(a.readgroup["DS",i.ar])
  if(!grepl("^NA",a.readgroup["DS",i.ar]) & !grepl("TruD:Agl1.2",a.readgroup["DS",i.ar])& !grepl("NEBuD;TruX",a.readgroup["DS",i.ar])  & !grepl("KapD;NimX3",a.readgroup["DS",i.ar])   & !grepl("TruD;NimX3",a.readgroup["DS",i.ar]) & !grepl("TruD:NimX",a.readgroup["DS",i.ar]) & !grepl("NxtD;NxtX",a.readgroup["DS",i.ar]) &  !grepl("NxtD;NxtXR",a.readgroup["DS",i.ar])  & !grepl("NxtD;NxtXE",a.readgroup["DS",i.ar])  & !grepl("TruD:NimX",a.readgroup["DS",i.ar]) & !grepl("TruD;",a.readgroup["DS",i.ar]) ){print("WARNING unknown capture technology in description field - see subroutine get.readgroup.info")}
  
  i.ar<-i.ar+1
    }}
colnames(a.readgroup)<-readgroup.cols
  
 a.readgroup
}

get.genome.info<-function(sample.bam,readgroup.labels=c("AS")){ #returms a matrix of readgroups
  the.header<-scanBamHeader(sample.bam) # list with names is sample.bams is a vertor then it does all of them at once!
  the.bam.files<-names(the.header)
  a.readgroup<-matrix(data=NA,nrow=length(readgroup.labels),ncol=1)
  rownames(a.readgroup)<-readgroup.labels
  for(ib in 1:length(the.bam.files)){
  the.readgroup<-the.header[[the.bam.files[ib]]][["text"]][["@SQ"]] ##only gets the first one
  the.RG.desc<-unlist(lapply(strsplit(the.readgroup,split=":"),function(x) x[1]))
  for(ir in 1:length(the.readgroup)){the.readgroup[ir]<-gsub(paste("^",the.RG.desc[ir],":",sep=""),"",the.readgroup[ir])}
  names(the.readgroup)<-the.RG.desc

    if(ib==1){a.readgroup[readgroup.labels,1]<-the.readgroup[readgroup.labels];readgroup.cols<-the.bam.files[ib]}else{
    a.readgroup<-cbind(a.readgroup,the.readgroup[readgroup.labels])
    readgroup.cols<-c(readgroup.cols,the.bam.files[ib])
  }
  }
  
  ## if(ib==1){a.readgroup[readgroup.labels,1]<-the.readgroup[readgroup.labels];colnames(a.readgroup)[1]<-the.bam.files[ib]}
  ## else{a.readgroup[names(the.readgroup),ib]<-the.readgroup;colnames(a.readgroup)[ib]<-the.bam.files[ib]}
  ##   }
  
  ##  a.readgroup["DS",ib]<-gsub("TruD;TruX","TruD:TruX",a.readgroup["DS",ib]) # special fix
  ##  a.readgroup["DS",ib]<-gsub("TruD;NimX","TruD:NimX",a.readgroup["DS",ib]) # special fix
  ##  ## "NA" nimblegen v1 - "TruD:NimX" nimblegen v2 - "TruD:TruX" illumina v2 - "TruD;" whole genome
  ## if(!is.na(a.readgroup["DS",ib]) & (!grepl("TruD:TruX",a.readgroup["DS",ib]) & !grepl("TruD:NimX",a.readgroup["DS",ib]) & !grepl("TruD;",a.readgroup["DS",ib]) )){print("WARNING unknown capture technology in description field - see subroutine get.readgroup.info")}
 colnames(a.readgroup)<-readgroup.cols
for(ic in 1:dim(a.readgroup)[2]){a.readgroup[is.na(a.readgroup[,ic]),ic]<-"NA"} # replace NA by "NA"
  a.readgroup
}

# sample.bam<-sample.bams[i]
get.targets.file<-function(sample.bam){

the.readgroup<-get.readgroup.info(sample.bam) # colnames is the bam file
the.genome<-get.genome.info(sample.bam)  # colnames is the bam file
the.capture<-extract.value.from.DS(the.readgroup["DS",],match.string="Description:")                                      # colnames is the bam file
names(the.capture)<-colnames(the.readgroup)

a.targets.file <- list()

for(i in 1:dim(the.readgroup)[2]){
  targets.file<-NA
  from.bam<-colnames(the.readgroup)[i]

eep <- "meaow"
print(eep)
################ Define different capture technologyies used  
if(( the.genome["AS",from.bam]=="UCSC_ALL_FULL_CHROMS_HG19.nix") | (the.genome["AS",from.bam]=="NA") ){  # a hg19 genome
eep <- "meaow2"
print(eep)
print(the.capture[from.bam])
    the.ref.genome<-"hg19";the.ref.genome.library<-"BSgenome.Hsapiens.UCSC.hg19"; the.ref.genome.object<-"Hsapiens" # used to get chromsome lengths
    if(the.capture[from.bam]=="TruD:TruX"){targets.file<-"Human_Exome_Targets_illumina_v2_hg19_targets.RData"}
 if(the.capture[from.bam]=="NxtD:NxtX"){targets.file<-"NexteraRapidCapture_ExpandedExome_TargetedRegions_hg19_targets.RData"} ## old one gets replaced by NxtD;NxtXE
if(the.capture[from.bam]=="NxtD:NxtXE"){targets.file<-"NexteraRapidCapture_ExpandedExome_TargetedRegions_hg19_targets.RData"}
if(the.capture[from.bam]=="NxtD:NxtXR"){targets.file<-"NexteraRapidCapture_Exome_TargetedRegions_v1.2_hg19_targets.RData"}
   # if(the.capture[from.bam]=="TruD:NimX"){targets.file<-"Human_Exome_Targets_Nimble_v2_hg19_targets.RData"} #(version=="v2")
  #  if(the.capture[from.bam]=="TruD:NimX3"){targets.file<-"Human_Exome_Targets_Nimble_v3_hg19_targets.RData"} #(version=="v2")TruD;Agl1.2
if(the.capture[from.bam]=="TruD:NimX"){targets.file<-"Human_Exome_Targets_Nimble_v3_hg19_targets.RData"} #(version=="v2")TruD;Agl1.2
if(the.capture[from.bam]=="TruD:NimX3"){targets.file<-"Human_Exome_Targets_Nimble_v3_hg19_targets.RData"} #(version=="v2")TruD;Agl1.2
if(the.capture[from.bam]=="KapD:NimX3"){targets.file<-"Human_Exome_Targets_Nimble_v3_hg19_targets.RData"} #(version=="v2")TruD;Agl1.2
if(the.capture[from.bam]=="NEBuD"){targets.file<-"Human_Exome_Targets_illumina_v2_hg19_targets.RData"} #(version=="v2")TruD;Agl1.2
    if(the.capture[from.bam]=="TruD:Agl1.2"){targets.file<-"Human_Exome_Targets_Nimble_v3_hg19_targets.RData"} #(version=="v2"
    if(the.capture[from.bam]==""){targets.file<-"Human_Exome_Targets_Nimble_v1_hg19_targets.RData"} #(version=="v1")
}

  
if(the.genome["AS",from.bam]=="mm9.nix"){
  the.ref.genome<-"mm9";the.ref.genome.library<-"BSgenome.Mmusculus.UCSC.mm9"; the.ref.genome.object<-"Mmusculus"
 if(the.capture[from.bam]=="TruD:NimX"){ targets.file<-"Mouse_Exome_Targets_Nimble_v100803_targets.RData"}
}
###############

a.target.info<-list(targets.file=targets.file,the.ref.genome=the.ref.genome,the.ref.genome.library=the.ref.genome.library,the.ref.genome.object=the.ref.genome.object)

if(i==1){a.targets.file<-list(a.target.info)}else{a.targets.file<-c(a.targets.file,list(a.target.info))}
}



names(a.targets.file)<-colnames(the.readgroup)
print(a.targets.file)
a.targets.file
}




basesAboveThreshold<-function(x,thresh=0){sum(x>=thresh)}

here <- "EEP"
print(here)
#print(targets.file)
args <- commandArgs(TRUE)
#print(args)
#print(args[1])
#print(args[2])

bam.file.to.do <- args[1]
do.projects.dir <- args[2]
#out.file <- args[2]


#print(the.file)
#the.file <- as.numeric(the.file)

## in the pipeline it is 0-23 so need to and 1 for the R code
#the.file <- the.file + 1
#print(the.file)
#print(out.file)
#q()

code.dir<-"/dmf/uqdi/Core_Services/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/"  ### code dir is where "auto_loop_over_UQCCG_coverage_with_cov.r" is located
## code.dir<-"/home/mmarshall/PBSHOME/PaulCoverageScriptsFeb2013"  ### code dir is where "auto_loop_over_UQCCG_coverage_with_cov.r" is located


## -> default


#################################### OPTIONS AND PATHS NEEDED
options(show.error.messages = TRUE)
#UQCCG.data<-"/mnt/GENO/Working/Pipeline/1359_0045_RSGB/" ### ROOT directory for projects
#UQCCG.data<-"/mnt/GENO/Working/Pipeline/89_0170/" ### ROOT directory for projects
UQCCG.data<-"/media/UQCCG/Sequencing/Projects" ### ROOT directory for projects
the.genome<-"hg19" ## the.genome<-"mm9"

num.cores<-10 # number of CPU coores to use
### genomes.path<-"/mnt/UQCCG/Sequencing/Data/Genomes/hg19" ## path to TARGETS files - the RDATA file that define the capture regions NAMES LIKE: "Human_Exome_Targets_illumina_v2_hg19_targets.RData"
genomes.path<-"/dmf/uqdi/Core_Services/UQCCG/Sequencing/Data/Genomes/hg19" ## path to TARGETS files - the RDATA file that define the capture regions NAMES LIKE: "Human_Exome_Targets_illumina_v2_hg19_targets.RData"

force.redo.with.new.ranges<-FALSE # force.redo.with.new.ranges<-TRUE # set true is want to recalculate EVERYTHING again  : IF FALSE it will over-write out of date files: see all.QC.column.labels
extend.exons<-200 # all to exons start/end to get additional coverage
##### need to choose one of the below:
################################################################33

options(cores=num.cores )
registerDoMC(cores=num.cores) ### Need to register cores for FOREACH method

###################################### EXAMPLES OF HOW TO SELECT SPECIFIC PROJECTS
restrict.analysis<-TRUE
## restrict.to.projects<- c("RSGB")         #AOGC FMHT-N,FTP 
## restrict.to.projects<- c("SKDP")         #AOGC FMHT-N,FTP 
## restrict.to.projects<- c("75_0119_SKDP")  
                                        #AOGC FMHT-N,FTP 
restrict.to.projects<- c("Aglient_clinical_exomes")  
## restrict.to.run<-"EXCLUDE" # "ALL" # "INCLUDE" # "EXCLUDE"
## the.run<-c("D1PVTACXXADADADAS") # only used with "INCLUDE" AND "EXCLUDE"
## the.run<-c("ReCal","DH10","DH09","00","05") # only used with "INCLUDE" AND "EXCLUDE"

restrict.to.run<-"" # "ALL" # "INCLUDE" # "EXCLUDE"
the.run<-""
#############

#restrict.analysis<-TRUE
#restrict.to.projects<- c("SKDP-FMDP","SKDP-FAM-12","SKDP-FAM-16","SKDP-FAM-2","SKDP-FAM-24","SKDP-FAM-36","SKDP-FAM-38","SKDP-FAM-4","SKDP-FAM-40","SKDP-FAM-41","SKDP-FAM-42","SKDP-FAM-44","SKDP-FAM-49","SKDP-FAM-51","SKDP-FAM-61","SKDP-FAM-66","SKDP-FAM-7","SKDP-FAM-70","SKDP-FAM-71","SKDP-FAM-73","SKDP-FAM-76","SKDP-FAM-80", "SKDP-FAM-84","SKDP-FAM-88","SKDP-FAM-90","SKDP-FAM-91","SKDP-FAM-92","SKDP-FAM-99")         #AOGC FMHT-N,FTP 
#restrict.to.run<-"EXCLUDE" # "ALL" # "INCLUDE" # "EXCLUDE"
#the.run<-"D0C3VABXX" # only used with "INCLUDE" AND "EXCLUDE"

#############

#restrict.analysis<-FALSE
#restrict.to.projects<-""
#restrict.to.run<-"" # "ALL" # "INCLUDE" # "EXCLUDE"
#the.run<-"" # only used with "INCLUDE" AND "EXCLUDE"
#############          

#############





############################################################################
#non.exome.projects<-c("BTCK","ALSPAC","AS-WGS","PLD","PLD-whole genome","AML-GENOME","bcos_srv","Mowry","QBI-DISC1")  ## things where I have not defined a capture method





######################################## no need to change below ###################
all.QC.column.labels<-c("ID","Sample","Recipe","Capture.Method","Lane","Run","target_file","total_reads","total_mapped_reads","total_dup_reads","unmapped_reads","total.bases",paste("on.target.bases.",extend.exons,"bp",sep=""),"on.target.bases", "percent.on.target.bases","median.max.coverage", "median.mean.coverage","mean.mean.coverage","percent.ccds.gt.1","percent.ccds.gt.5","percent.ccds.gt.10","percent.ccds.gt.15","percent.ccds.gt.30","percent.ccds.gt.40","percent.ccds.gt.50","percent.ccds.gt.60","percent.ccds.gt.70","percent.ccds.gt.80","percent.ccds.gt.90","percent.ccds.gt.100","percent_Duplicated","percent_Unmapped","Description")

projects<-dir(UQCCG.data)
projects<-projects[!(projects %in% non.exome.projects)]
projects
##

all.data<-{}
count<-0


if(restrict.analysis){
 projects<-projects[projects %in% restrict.to.projects]
}

print(projects)
targets.file.ori<-""
####### ip<-1
for(ip in 1:length(projects)){
print(projects[ip])

do.projects.dir<-paste(UQCCG.data,projects[ip],sep="/")
on.target.stats<-{}  # this will be rebuilt from the individual 

BAM.directory<-paste(do.projects.dir,"BAM",sep="/")
#BAM.directory<-paste(UQCCG.data,projects[ip],"BAM",sep="/")
x<-try(setwd(BAM.directory ))
if(inherits(x, "try-error")){next}

QC.root.directory<-paste(do.projects.dir,"QC",sep="/")
#QC.root.directory<-paste(UQCCG.data,projects[ip],"QC",sep="/")
xx<-try(setwd( QC.root.directory  ),silent=TRUE)
if(inherits(xx, "try-error")){system(paste("mkdir",QC.root.directory,sep=" "));setwd( QC.root.directory  )}


QC.directory<-paste(do.projects.dir,"QC/AlignedQC",sep="/")
#QC.directory<-paste(UQCCG.data,projects[ip],"QC/AlignedQC",sep="/")
xx<-try(setwd( QC.directory  ),silent=TRUE)
if(inherits(xx, "try-error")){system(paste("mkdir",QC.directory,sep=" "));setwd( QC.directory);QC.files<-dir(getwd())}else{QC.files<-dir(getwd())}
  ## QC.files<-QC.files[grepl(".QC$",QC.files)]


setwd(BAM.directory)
#files<-dir(getwd())
#wow <- "get files"
#print(wow)
#print(files)
#sample.bams<-files[grepl(".ReCal.sort.bam$",files)] # sample.bams<-sample.bams.ori[48:132] # sample.bams<-sample.bams.ori[1:47] # sample.bams.ori<-sample.bams

#if((restrict.to.run=="INCLUDE") | (restrict.to.run=="EXCLUDE")){
#  if(restrict.to.run=="INCLUDE"){
#    wanted<-rep(FALSE,length(sample.bams))
#    for(irun in 1:length(the.run)){wanted<-wanted | grepl(the.run[irun],sample.bams)}  # grep("21",sample.bams)
#    sample.bams<-sample.bams[wanted]
#  }
#  if(restrict.to.run=="EXCLUDE"){
#        wanted<-rep(FALSE,length(sample.bams))
#    for(irun in 1:length(the.run)){wanted<-wanted | grepl(the.run[irun],sample.bams)}  # grep("21",sample.bams)
#    sample.bams<-sample.bams[!wanted]
#  }}

#print(sample.bams)
#exit
#k<-1
# i<-1
here <- "here"
print(here)
###i <- the.file OLD WHEN USED $PBS_ARRAY_INDEX
#for(i in 1:length(sample.bams)){

## the bam file is passed in from the pbs file
print(bam.file.to.do)
#print(sample.bams[i])

the.header<-scanBamHeader(bam.file.to.do) # list with names
#the.header<-scanBamHeader(sample.bams[i]) # list with names

the.bam.files<-names(the.header)
print(the.bam.files)
sample.bam.counts<-{}
total.reads<-{}
total.coverage<-{}

the.readgroup<-get.readgroup.info(bam.file.to.do)
the.genome<-get.genome.info(bam.file.to.do)
the.target.files<-get.targets.file(bam.file.to.do)
#the.readgroup<-get.readgroup.info(sample.bams[i])
#the.genome<-get.genome.info(sample.bams[i])
#the.target.files<-get.targets.file(sample.bams[i])



# ii<-1
print(here)
#ii <- the.file
print(here)
for( ii in 1:length(the.bam.files)){  ### loop over bam files
  print(the.bam.files[ii])
print(here)

  ###### check if already done:
sample.bam.counts<-paste(gsub(".bam$","",basename((the.bam.files[ii]))),"QC",sep=".")

  ###################### THIS DOESNT WORK ON CLUSTER ###########
#if(!force.redo.with.new.ranges){
  here <- "here 1"
  print(here)
#if( (sample.bam.counts %in% QC.files)  & (paste(sample.bam.counts,"regions.RData",sep=".") %in% QC.files) ){ ## probably has been done
#  QC.file.chk<-read.delim(paste(QC.directory,sample.bam.counts,sep="/"),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
#  if( ("ID" %in% rownames(QC.file.chk)) & ("Capture.Method"  %in% rownames(QC.file.chk)) &  length(all.QC.column.labels)==length(rownames(QC.file.chk))  ){ # this file has been processed and is ok
on.target.stats<-{}
#       if( is.null(on.target.stats) ){ on.target.stats<-as.matrix(QC.file.chk[all.QC.column.labels,])}else{on.target.stats<-cbind(on.target.stats,QC.file.chk[all.QC.column.labels,])} # complile on.target.stats
#       print(paste(the.bam.files[ii],"Already DONE skipping",sep=" "))
#       next# now force skip of analysis below
 #                  } 
#}else{print(paste(the.bam.files[ii],"NOT DONE running",sep=" "))}

#} # end force redo
##########################

the.sample<-the.readgroup["SM",ii]
  here <- "here 12"
  print(here)
the.current.bam<-the.bam.files[ii]
here <- "here 2"
  print(here)
a.run<- extract.value.from.DS(the.readgroup["DS",ii],match.string="FCID:") # runID
## extract.value.from.DS(the.readgroup["DS",],match.string="SampleRef:") # genome ref
a.lane<-extract.value.from.DS(the.readgroup["DS",ii],match.string="Lane:") # LAne
a.capture<- extract.value.from.DS(the.readgroup["DS",ii],match.string="Description:") # capture.tech
a.recipe<- extract.value.from.DS(the.readgroup["DS",ii],match.string="Recipe:") # recipe
a.sample<-the.sample # extract.value.from.DS(the.readgroup["DS",ii],match.string="ParticipantCode:") #sample  DS does not always exist
a.sample.ID<-extract.value.from.DS(the.readgroup["DS",ii],match.string="SampleID:")  # ID
a.DS<-the.readgroup["DS",ii]
here <- "pah 2"
print(here)
targets.file<-the.target.files[[ii]]$targets.file

if(!exists("targets.file.ori")){targets.file.ori<-targets.file} # for first time through

  ######################################## Load the approproate Genome package = "BSgenome.Hsapiens.UCSC.hg19" typically
  if(is.na(targets.file)){print("ERROR targets file not obtained");next}

if(targets.file != targets.file.ori | !exists("the.chroms")){

  load(paste(genomes.path,targets.file,sep="/"))
 
  print(paste(genomes.path,targets.file,sep="/"))
  package<-the.target.files[[ii]]$the.ref.genome.library
  
 
  library(package,character.only=TRUE)

  the.chroms<-seqlengths(eval(as.name(the.target.files[[ii]]$the.ref.genome.object)))

  genome(data.gr)<-the.target.files[[ii]]$the.ref.genome
  gr<-data.gr

  rl.from.gr<-as(data.gr, "RangesList") # data.gr contains
  rl.from.gr.expanded<-rl.from.gr+extend.exons
  rl.from.gr.expanded<-reduce(rl.from.gr.expanded)
  targets.file.ori<-targets.file
  print("Got new targets file")
}


here <- "pah 3"
print(here)
  
  
the.bam.chrs<-names(the.header[[the.bam.files[ii]]]$targets) ### caus the bam chromosome names migh be differe to yours : chr1 vs chromosme1 etc
the.bam.chr.lengths<-as.integer(the.header[[the.bam.files[ii]]]$targets)
human.chromlens<-the.chroms[the.bam.chrs]


  

##   ########################## NO MULTICORE
## #isPaired=TRUE,isUnmappedQuery=FALSE,hasUnmappedMate=FALSE) what=c("rname", "pos", "cigar")
## #ScanBamParam()
##   # system.time({
## param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=FALSE, isPaired=NA)) # don't use what=c("rname", "pos", "cigar") as it is additionally added to meta cols (two cigars)
## the.reads <-readBamGappedAlignments(the.bam.files[ii],param=param)  ##20 sec for whole exome GappedAlignments of length 22373288 # takes about 150s
## #})
## ##############################################


## ############################ without multicore 
## paramUnmapped <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=TRUE),what=c("qname"))
## unmapped.reads<-countBam(the.bam.files[ii],param=paramUnmapped)

## paramDup <- ScanBamParam(scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=TRUE, isPaired=NA),what=c("rname"))  # in testing duplicates were all mapped
## dup.reads<-countBam(the.bam.files[ii],param=paramDup)
## ##############################################

###### unmapped reads count
count.reads<- mclapply(names(human.chromlens), countTotalReadsInBam, the.bam.files[ii], human.chromlens)
  
records<-sum(unlist(sapply(count.reads,get.list.elements,"records")),na.rm=TRUE)
nucleotides<-sum(unlist(sapply(count.reads,get.list.elements,"nucleotides")),na.rm=TRUE)
file<-unique(unlist(sapply(count.reads,get.list.elements,"file")))
total.reads.have<-data.frame(file=file,records=records,nucleotides=nucleotides)
#####
## > total.reads.have
##                                         file  records nucleotides
## 1 AOGC-01-0003_C0356ABXX-1-01.ReCal.sort.bam 49276580  4976934580
  
###### unmapped reads count
count.reads<- mclapply(names(human.chromlens), countUnmappedReadsInBam, the.bam.files[ii], human.chromlens)
  
records<-sum(unlist(sapply(count.reads,get.list.elements,"records")),na.rm=TRUE)
nucleotides<-sum(unlist(sapply(count.reads,get.list.elements,"nucleotides")),na.rm=TRUE)
file<-unique(unlist(sapply(count.reads,get.list.elements,"file")))
unmapped.reads<-data.frame(file=file,records=records,nucleotides=nucleotides)
#####
  
##### dup read count 
count.reads<- mclapply(names(human.chromlens), countDupReadsInBam, the.bam.files[ii], human.chromlens)

records<-sum(unlist(sapply(count.reads,get.list.elements,"records")),na.rm=TRUE)
nucleotides<-sum(unlist(sapply(count.reads,get.list.elements,"nucleotides")),na.rm=TRUE)
file<-unique(unlist(sapply(count.reads,get.list.elements,"file")))
dup.reads<-data.frame(file=file,records=records,nucleotides=nucleotides)
##### 
  
###################################################
### code chunk number 13: summaryByChromosome
###################################################  NOTE :: SAVE MEMORY  TO DO in chunks put loop over here and coverage Cov is small and the .reads is large!!
                                             #####  Add coverge calculation to extract reads from BAM  !!! This done per chromosome
#system.time({
## the.reads <- mclapply(names(human.chromlens)[18:19], extractReadsFromBam, the.bam.files[ii], human.chromlens) #multicore read per chromosome takes about 57s
## the.reads <- do.call(c,the.reads)
print(paste("Calculate coverage using: ", getOption("cores")," cores",sep=""))
cov<-foreach(the.chr=iter(names(human.chromlens),chunksize=1),.combine='c', .multicombine=TRUE,.inorder=TRUE) %dopar% extractCovFromBamOneChr(the.chr,the.bam.files[ii], human.chromlens)
print(paste("Got coverge for",the.bam.files[ii],sep=" "))
## the.cov<-unlist(lapply(cov.all,sum))
## cov<-cov.all[the.cov!=0]
## cov
## lapply(cov,sum)
############################################## 
## cov.all<- mclapply(names(human.chromlens), extractCovFromBam, the.bam.files[ii], human.chromlens)
##  ## This cov.all is names(human.chromlens) long: each elements is also names(human.chromlens)long ones one of those has data
##  ## I assue here the the order of names(human.chromlens) is preserved
## for(icov in 1:length(cov.all)){
##    cov.all[[icov]]<-
## if(length(names(human.chromlens))>1){
##   cov<-cov.all[[1]]
##   for(icov in 2:length(cov.all)){
##   cov[[names(human.chromlens)[icov]]]<- cov.all[[icov]][[names(human.chromlens)[icov]]]
##                       }}
########################################################################################################
  ## seqname<-names(human.chromlens)[18]
  ##  bamFile<-the.bam.files[ii]

########Multicore run 
## print(paste("START processing",samples.processing,"samples",sep=" "))  

## ###

  
#######################################################################################
###construct the column names



##   paramAll <- ScanBamParam(scanBamFlag(isUnmappedQuery=NA, isDuplicate=NA, isPaired=NA),what=c("rname"))
## all.reads<-countBam(the.bam.files[ii],param=paramAll) # same as  (total.reads+total.unmapped.reads+total.dup.reads)
##############the reads contain all the sequence reads 
#total.reads<-length(the.reads)
  
total.reads<-as.vector(as.numeric(total.reads.have["records"]))
total.unmapped.reads<-as.vector(as.numeric(unmapped.reads["records"]))
total.dup.reads<-as.vector(as.numeric(dup.reads["records"]))
sum.total.reads<-as.vector(as.numeric(total.reads))

total.unmapped.reads/(sum.total.reads)
total.dup.reads/(sum.total.reads-total.unmapped.reads)
total.dup.reads/(sum.total.reads) # this is what picard mark duplicated provides



## system.time({
## cov<-coverage(the.reads,width=human.chromlens) # readBamGappedAlignments use this is have a RangesList and RangedData objects where they are the length of each list
##   rm("the.reads")
## })
# cov<-coverage(the.reads) #extractReadsFromBAM


  ######################################################
  #####################################################
  # Internal loop that processes the cov object

print("about to start -> auto_loop_over_UQCCG_coverage_with_cov")
source(paste(code.dir,"auto_loop_over_UQCCG_coverage_with_cov.r",sep="/"))  ## This ones uses a max of 8 cores 
## source(paste(code.dir,"auto_loop_over_UQCCG_coverage_with_cov_foreach.r",sep="/")) ## should be better is have lots of cores
## foreach version does::3541.770 1160.620 1048.636| 3883.990 1355.240 1138.665 (7|8 cores)   parrallel does::3439.870  209.390  807.275 
  ######################################################
  #####################################################

print("done auto loop over UQCCG coverage with cov")
########################## WRITE OUTPUT
 colnames(data)<-gsub(".QC$",".bam",sample.bam.counts)
print(here)
 write.table(data,file=paste(QC.directory,sample.bam.counts,sep="/"),col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE)
 save(list=c("the.counts","cov","data"),file=paste(QC.directory,paste(sample.bam.counts,"regions.RData",sep="."),sep="/"))   
 rm(the.counts,cov)
print(here)
  
if( is.null(on.target.stats) ) { on.target.stats<-data}else{
 on.target.stats<-cbind(on.target.stats,data)}

#k<-k+1

  ## java -Xmx2g -jar /media/Bioinform-D/Research/Picard/picard-tools-1.61/picard-tools-1.61/CalculateHsMetrics.jar BAIT_INTERVALS=/media/scratch/Genomes/ill_baits.list TARGET_INTERVALS=/media/scratch/Genomes/ill_baits.list INPUT=AOGC-14-2169_D0C59ABXX-3-10.ReCal.sort.bam  OUTPUT=AOGC-14-2169_D0C59ABXX-3-10.ReCal.sort

  ##  java -Xmx2g -jar /media/Bioinform-D/Research/Picard/picard-tools-1.61/picard-tools-1.61/CollectMultipleMetrics.jar   INPUT=D0C4LABXX-8-06.ReCal.sort.bam OUTPUT=D0C4LABXX-8-06.ReCal.sort

  
} # loop over the.bam.files (incase may SM in each BAM file)

#} # loop over sample.bams 
print(here)

rownames(on.target.stats)<-all.QC.column.labels # may not be set if combed run and non-run samples

write.table(on.target.stats,file=paste(QC.directory,"quality.control.summary",sep="/"),col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE,append=TRUE)
#write.table(on.target.stats,file=paste(QC.directory,paste(projects[ip],"quality.control.summary",sep="."),sep="/"),col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE,append=TRUE)



######################Now compile summaries where there are samples with multiple bam files assume these are in the same project file




#} # loop over projects

