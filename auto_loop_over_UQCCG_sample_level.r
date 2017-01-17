
library(GenomicFeatures)
library("doParallel")
##
try(library(multicore),silent=TRUE)
library(Rsamtools)

## num.cores<-5
##  registerDoParallel(cores=num.cores)
host = system("hostname",intern=T)
if(host %in% c("DI-LW-BRN-011") ){
  root.path="/mnt/UQCCG"
}else{
  root.path="/dmf/uqdi/Core_Services/UQCCG"
}

code.dir<-paste(root.path,"/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts",sep="")
setwd(code.dir)
source("annotate_SNPs_subroutines.r")

the.genome<-"hg19" ## the.genome<-"mm9"

UQCCG.data<-paste(root.path,"/Sequencing/Projects",sep="")
genomes.path<-paste(root.path,"/Sequencing/Data/Genomes/hg19",sep="") ## path to directory with capture location Rdata files



options(show.error.messages = TRUE)
possible.capture.methods<-c("TruD:TruX","TruD:NimX","TruD:NimX3","TruD:Agl1.2","NA","NxtD:NxtXE","NxtD:NxtXR")
#UQCCG.data<-"/mnt/GENO/Working/Pipeline/1359_0035_IYMD/IYMD"
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/pipeline/140613_SN7001359_0059_AC4MJJACXX/"
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/pipeline/140312_SN7001359_0051_AH801VADXX/"


## ############
## restrict.analysis<-TRUE
## restrict.to.projects<- c("AMAS")        #AOGC FMHT-N,FTP 
## restrict.to.run<-"" # "ALL" # "INCLUDE" # "EXCLUDE"
## the.run<-"" # only used with "INCLUDE" AND "EXCLUDE"
## ############# 



## ############
restrict.analysis<-FALSE
restrict.to.projects<- ""        #AOGC FMHT-N,FTP 
restrict.to.run<-"ALL" # "ALL" # "INCLUDE" # "EXCLUDE"
the.run<-"" # only used with "INCLUDE" AND "EXCLUDE"
## #############


#############


## path to genomes files
force.redo.with.new.ranges<-FALSE # force.redo.with.new.ranges<-TRUE # set true is want to recalculate
extend.exons<-200 # all to exons start/end to get additional coverage
##### need to choose one of the below:


############################################################################

non.exome.projects<-c("ALSPAC","LKAF","Head and Neck TCGA","AML-PacBio","BTCK","ALSPAC","AS-WGS","PLD","PLD-whole genome","AML-GENOME","bcos_srv","Mowry","QBI-DISC1","directory.structure.txt.gz")



######################################## no need to change below ###################
all.QC.column.labels<-c("ID","Sample","Recipe","Capture.Method","Lane","Run","target_file","total_reads","total_mapped_reads","total_dup_reads","unmapped_reads","total.bases",paste("on.target.bases.",extend.exons,"bp",sep=""),"on.target.bases", "percent.on.target.bases","median.max.coverage", "median.mean.coverage","mean.mean.coverage","percent.ccds.gt.1","percent.ccds.gt.5","percent.ccds.gt.10","percent.ccds.gt.15","percent.ccds.gt.30","percent_Duplicated","percent_Unmapped","Description")

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
####### ip<-30
for(ip in 1:length(projects)){
print(projects[ip])
on.target.stats<-{}  # this will be rebuilt from the individual 

BAM.directory<-paste(UQCCG.data,projects[ip],"BAM",sep="/")
x<-try(setwd(BAM.directory ))
if(inherits(x, "try-error")){next}

QC.root.directory<-paste(UQCCG.data,projects[ip],"QC",sep="/")
xx<-try(setwd( QC.root.directory  ),silent=TRUE)
if(inherits(xx, "try-error")){system(paste("mkdir",QC.root.directory,sep=" "));setwd( QC.root.directory  )}


QC.directory<-paste(UQCCG.data,projects[ip],"QC/AlignedQC",sep="/")
xx<-try(setwd( QC.directory  ),silent=TRUE)
if(inherits(xx, "try-error")){system(paste("mkdir",QC.directory,sep=" "));setwd( QC.directory);QC.files<-dir(getwd())}else{QC.files<-dir(getwd())}
  ## QC.files<-QC.files[grepl(".QC$",QC.files)]


setwd(BAM.directory)
files<-dir(getwd())
sample.bams<-files[grepl("\\.bam$",files)] # sample.bams<-sample.bams.ori[48:132] # sample.bams<-sample.bams.ori[1:47] # sample.bams.ori<-sample.bams

if((restrict.to.run=="INCLUDE") | (restrict.to.run=="EXCLUDE")){
  if(restrict.to.run=="INCLUDE"){
    wanted<-rep(FALSE,length(sample.bams))
    for(irun in 1:length(the.run)){wanted<-wanted | grepl(the.run[irun],sample.bams)}  # grep("21",sample.bams)
    sample.bams<-sample.bams[wanted]
  }
  if(restrict.to.run=="EXCLUDE"){
        wanted<-rep(FALSE,length(sample.bams))
    for(irun in 1:length(the.run)){wanted<-wanted | grepl(the.run[irun],sample.bams)}  # grep("21",sample.bams)
    sample.bams<-sample.bams[!wanted]
  }}

print(sample.bams)


k<-1
all.readgroup<-get.readgroup.info(sample.bams)
all.genome<-get.genome.info(sample.bams)
all.target.files<-get.targets.file(sample.bams)
  

# duplicated.samples<-all.readgroup["SM",][duplicated(all.readgroup["SM",])] ### spoofed by add bams have two samples 
# duplicated.samples<-all.readgroup["SM",] ## can contain mutilple read groups
duplicated.samples<-all.readgroup["SM",][!duplicated(names(all.readgroup["SM",]))] ## ** coverage at bam level so look for non-dup bam names
duplicated.samples<-duplicated.samples[duplicated(duplicated.samples)]

 #  merged BAM files with same SM can make 3 apparent files **coverage done at BAM file level not SM level**
if(length(names(duplicated.samples)) > length(unique(names(duplicated.samples)) ) ){
duplicated.samples<-duplicated.samples[!duplicated(names(duplicated.samples))] # only unique BAm files
duplicated.samples<-duplicated.samples[duplicated(duplicated.samples)] # duplicated samples
}

### added by Mhairi - Jan 2013
duplicated.samples <- unique(duplicated.samples)

  if(length(duplicated.samples)==0){next}
## i <- 1
#print(duplicated.samples)
#length(duplicated.samples)
#q()



for(i in 1:length(duplicated.samples)){
####for(i in 8:length(duplicated.samples)){ 
  
sample.output<-paste(duplicated.samples[i],".combined.ReCal.sort.QC",sep="")   # names of file data is saved to

the.readgroup<-all.readgroup[,all.readgroup["SM",]==duplicated.samples[i] ]
the.bam.files<-colnames(the.readgroup)
bam.roots<-dirname((the.bam.files))
sample.bam.counts<-paste(gsub(".bam$","",basename((the.bam.files))),"QC",sep=".")
cov.objects<-paste(sample.bam.counts,"regions.RData",sep=".")
names(cov.objects)<-the.bam.files

######### check the cov.objects exist ##
cov.objects<-cov.objects[cov.objects %in% QC.files]

#if(length(cov.objects)==0){} # chk to make use can read the coverage object


if(length(cov.objects)==0){print("WARNING These bam files nor preprocess -skipping");next}
print("past next")
  ###################### already run so skip out
if(!force.redo.with.new.ranges){
  
if( (sample.output %in% QC.files) ){ ## probably has been done
  print("in test loop next")
  QC.file.chk<-read.delim(paste(QC.directory,sample.output,sep="/"),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
  present<-rownames(QC.file.chk) %in% all.QC.column.labels 
  QC.file.chk<-subset(QC.file.chk,subset=present)  ### Get rid of extra informtion that might  have been collected
  
  if( ("ID" %in% rownames(QC.file.chk)) & ("Capture.Method"  %in% rownames(QC.file.chk)) &  sum(all.QC.column.labels %in% rownames(QC.file.chk))==length(all.QC.column.labels)   ){ # this file has been processed and is ok
       if( is.null(on.target.stats) ){ on.target.stats<-as.matrix(QC.file.chk[all.QC.column.labels,])}else{on.target.stats<-cbind(on.target.stats,QC.file.chk[all.QC.column.labels,])} # complile on.target.stats
       print(paste(duplicated.samples[i],"DONE",sep=" "))
       next# now force skip of analysis below
                   } 
}else{print(paste(duplicated.samples[i],"NOT combined running",sep=" "))}

} # end force redo
##########################



 ############### load ths coverage objects and check they are up to date

 ############### NOte they could be on different genomes with different chromsomes definitions so addition of cov objects could fail the chr hanes are from the targets file
 ############## so a cov object with all chromsomes is made
############# NOte cov object of different lensthsor with different names can't be simply added together

 
setwd(QC.directory)
old.cov.object<-rep(FALSE,times=length(cov.objects))
 cov.all<-{}
 data.all<-{}


ic<-1
for(ic in 1:length(cov.objects)){
  print(cov.objects[ic])
  load(cov.objects[ic])
  if( !( ("ID" %in% rownames(data)) & ("Capture.Method"  %in% rownames(data)) &  sum(all.QC.column.labels %in% rownames(data))==length(all.QC.column.labels)    ) ){old.cov.object[ic]<-TRUE;next} #the data object denotes it is old
  if(is.null(cov.all)){cov.all<-cov;chrs.known<-names(cov.all)}else{
    chrs.known<-names(cov.all)
    new.chrs<-names(cov)
    common.chrs<-new.chrs[new.chrs %in% chrs.known]
    new.chrs<-new.chrs[!(new.chrs %in% chrs.known)]
    if(length(common.chrs)>0){
      for(icc in 1:length(common.chrs)){ cov.all[common.chrs[icc]]<- cov.all[common.chrs[icc]] +  cov[common.chrs[icc]] } # all together chromosomes
       }
    if(length(new.chrs)>0){cov.all<-c(cov.all,cov[new.chrs]);print(paste("WARNING new chromosomes:",toString(new.chrs),"WITHIN:", toString(cov.objects),"Adding",sep=" ")) }
       } # new cove object

  #############may loop over man opbjects some have old coverage stats
  data<-as.matrix(data[all.QC.column.labels,])
  rownames(data)<-all.QC.column.labels
  if(is.null(data.all)){data.all<-data }else{data.all1<-cbind(data.all,data[all.QC.column.labels,])}
} # loop over coverage



cov<-cov.all
data.all
sum(old.cov.object)
#data.all<-cbind(data.all[all.QC.column.labels,],data[all.QC.column.labels,],stringsAsFactors=FALSE)

cov.objects<-cov.objects[!old.cov.object]
rm(cov.all)
rm(the.counts)
 



if(length(cov.objects)<2){next} # there are not enough RDATA files to proceed


the.readgroup<-the.readgroup[,names(cov.objects)]
the.genome<-all.genome[,names(cov.objects)]
the.target.files<-all.target.files[names(cov.objects)]
                             
the.sample<-duplicated.samples[i]

the.current.bam<-paste(the.bam.files,collapse=";")
  
a.run<- paste(extract.value.from.DS(the.readgroup["DS",],match.string="FCID:"),collapse=";") # runID
## extract.value.from.DS(the.readgroup["DS",],match.string="SampleRef:"),collapse=";") # genome ref
a.lane<-paste(extract.value.from.DS(the.readgroup["DS",],match.string="Lane:"),collapse=";") # LAne
a.capture<- paste(extract.value.from.DS(the.readgroup["DS",],match.string="Description:"),collapse=";") # capture.tech
a.recipe<- paste(extract.value.from.DS(the.readgroup["DS",],match.string="Recipe:"),collapse=";") # recipe
a.sample<-the.sample
a.sample.ID<-paste(extract.value.from.DS(the.readgroup["DS",],match.string="SampleID:"),collapse=";")  # ID
a.DS<-paste(the.readgroup["DS",],collapse=";")


 num.genomes.same<-length(unique(the.genome))
## THIS was TAKE OUT TO MERGE ChMND as differnt bams although both on hg19 had different hg19.nix
if(num.genomes.same>1){print(paste("ERROR for",toString(the.bam.files),"genomes differ can't combine",sep=" "));next}

 ####################### find the lost inclusive target file and determine it's position in 
possible.targets<- extract.value.from.DS(the.readgroup["DS",],match.string="Description:")
num.possible.targets<-length(unique(the.genome))
if(num.possible.targets==1){ii<-1}else{
# possible.capture.methods<-c("TruD:TruX","TruD:NimX","NA") # defined above
  posns<-match( possible.targets,possible.capture.methods) ## Mhairi changed capture.methods -> possible.capture.methods
  ii<-which.min(posns)
}  #ii not contains the best targets option 

  
  ########################################

targets.file<-the.target.files[[ii]]$targets.file
  if(is.na(targets.file)){print("ERROR targets file not obtained");next}

if(targets.file != targets.file.ori){
  load(paste(genomes.path,targets.file,sep="/"))
  package<-the.target.files[[ii]]$the.ref.genome.library
  library(package,character.only=TRUE)
  the.chroms<-seqlengths(eval(as.name(the.target.files[[ii]]$the.ref.genome.object)))

  genome(data.gr)<-the.target.files[[ii]]$the.ref.genome
  gr<-data.gr

  rl.from.gr<-as(data.gr, "RangesList") # data.gr contains
  rl.from.gr.expanded<-rl.from.gr+extend.exons
  rl.from.gr.expanded<-reduce(rl.from.gr.expanded)
  targets.file.ori<-targets.file
}

human.chromlens<-the.chroms[names(cov)]
if(length(human.chromlens)<15){
  print("WARNING fewer than expected chomosomes")
  if(length(human.chromlens)<1){ print("ERROR NO chomosomes");next}
   }

##   paramAll <- ScanBamParam(scanBamFlag(isUnmappedQuery=NA, isDuplicate=NA, isPaired=NA),what=c("rname"))
## all.reads<-countBam(the.bam.files[ii],param=paramAll) # same as  (total.reads+total.unmapped.reads+total.dup.reads)
##############the reads contain all the sequence reads 
total.reads<-sum(as.numeric(as.matrix(data.all["total_mapped_reads",])) , na.rm=TRUE)
total.unmapped.reads<-sum(as.numeric(as.matrix(data.all["unmapped_reads",])) , na.rm=TRUE)
total.dup.reads<-sum(as.numeric(as.matrix(data.all["total_dup_reads",])) , na.rm=TRUE) ## note this is not really true would need to run picard on the 2 files!
sum.total.reads<-sum(as.numeric(as.matrix(data.all["total_reads",])) , na.rm=TRUE)


## change Aug 2013

total.unmapped.reads/(sum.total.reads)
total.dup.reads/(sum.total.reads-total.unmapped.reads)
total.dup.reads/(sum.total.reads) # this is what picard mark duplicated provides


## total.unmapped.reads/(total.reads+total.unmapped.reads+total.dup.reads)
## total.dup.reads/(total.reads)
## total.dup.reads/(total.reads+total.unmapped.reads+total.dup.reads) # this is what picard mark duplicated provides

  ######################################################
  #####################################################
  # Internal loop that processes the cov object
 setwd(code.dir)
source("auto_loop_over_UQCCG_coverage_with_cov.r")
 setwd(QC.directory)
  ######################################################
  #####################################################

########################## WRITE OUTPUT
 
 colnames(data)<-sample.output
 write.table(data,file=paste(QC.directory,sample.output,sep="/"),col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE)
 save(list=c("the.counts","cov","data"),file=paste(QC.directory,paste(sample.output,"regions.RData",sep="."),sep="/"))   
  

if( is.null(on.target.stats) ) { on.target.stats<-data}else{
 on.target.stats<-cbind(on.target.stats,data)}

k<-k+1
  
} # sample loop


if(!is.null(on.target.stats)){ ## only do these lines if got extra data
  
rownames(on.target.stats)<-all.QC.column.labels # may not be set if combed run and non-run samples
xx<-try(write.table(on.target.stats,file=paste(QC.directory,paste(projects[ip],"quality.control.summary",sep="."),sep="/"),col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE), silent=TRUE)
if(inherits(xx, "try-error")){write.table(on.target.stats,file=paste(QC.directory,paste(projects[ip],"quality.control.summary.1",sep="."),sep="/"),col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE)   }


######################Now compile summaries where there are samples with multiple bam files assume these are in the same project file

} # test is.null(on.target.stats)


} # loop over projects





