




options(show.error.messages = TRUE)
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/pipeline/NEX_140613_SN7001359_0059_AC4MJJACXX_140618_SN989_0198_BC4LL7ACXX/"
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/pipeline/140613_SN7001359_0059_AC4MJJACXX/"
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/pipeline/140424_SN989_0195_BC4LV1ACXX/"
UQCCG.data<-"/mnt/UQCCG/Sequencing/Projects"
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/pipeline/141110_SN7001359_0068_BHBE20ADXX/"
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/pipeline/130513_7001408_0052_AC1JRDACXX_QBI/"
combined.on.target.stats.dir<-"/mnt/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC"
combined.output.file.bam<-paste("QC_stat_BAM_",format(Sys.time(), "%a_%b_%d_%Y"),".txt",sep="")
combined.output.file.sample<-paste("QC_stat_SAMPLE_",format(Sys.time(), "%a_%b_%d_%Y"),".txt",sep="")
code.dir<-"/mnt/UQCCG/Sequencing/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
QC.extensions<-c(".ReCal.sort.QC",".sort.QC")

## #############
## restrict.analysis<-TRUE
## restrict.to.projects<- c("SKDP-FMDP","SKDP-FAM-12","SKDP-FAM-16","SKDP-FAM-2","SKDP-FAM-24","SKDP-FAM-36","SKDP-FAM-38","SKDP-FAM-4","SKDP-FAM-40","SKDP-FAM-41","SKDP-FAM-42","SKDP-FAM-44","SKDP-FAM-49","SKDP-FAM-51","SKDP-FAM-61","SKDP-FAM-66","SKDP-FAM-7","SKDP-FAM-70","SKDP-FAM-71","SKDP-FAM-73","SKDP-FAM-76","SKDP-FAM-80", "SKDP-FAM-84","SKDP-FAM-88","SKDP-FAM-90","SKDP-FAM-91","SKDP-FAM-92","SKDP-FAM-99")         #AOGC FMHT-N,FTP 
## restrict.to.run<-"EXCLUDE" # "ALL" # "INCLUDE" # "EXCLUDE"
## the.run<-"D0C3VABXX" # only used with "INCLUDE" AND "EXCLUDE"
## all.QC.column.labels<-c("ID","Sample","Recipe","Capture.Method","Lane","Run","target_file","total_reads","total_mapped_reads","total_dup_reads","unmapped_reads","total.bases",paste("on.target.bases.",extend.exons,"bp",sep=""),"on.target.bases", "percent.on.target.bases","median.max.coverage", "median.mean.coverage","mean.mean.coverage","percent.ccds.gt.1","percent.ccds.gt.5","percent.ccds.gt.10","percent.ccds.gt.15","percent.ccds.gt.30","percent_Duplicated","percent_Unmapped","Description")
## ################

## #############
## restrict.analysis<-FALSE
## restrict.to.projects<- ""        #AOGC FMHT-N,FTP 
## restrict.to.run<-"" # "ALL" # "INCLUDE" # "EXCLUDE"
## the.run<-"" # only used with "INCLUDE" AND "EXCLUDE"
## all.QC.column.labels<-c("ID","Sample","Recipe","Capture.Method","Lane","Run","target_file","total_reads","total_mapped_reads","total_dup_reads","unmapped_reads","total.bases",paste("on.target.bases.",extend.exons,"bp",sep=""),"on.target.bases", "percent.on.target.bases","median.max.coverage", "median.mean.coverage","mean.mean.coverage","percent.ccds.gt.1","percent.ccds.gt.5","percent.ccds.gt.10","percent.ccds.gt.15","percent.ccds.gt.30","percent_Duplicated","percent_Unmapped","Description")
## ################

## restrict.analysis<-TRUE
## restrict.to.projects<- c("2012-04-30_QBI")         #AOGC FMHT-N,FTP 
## restrict.to.run<-"" # "ALL" # "INCLUDE" # "EXCLUDE"
## the.run<-"" # only used with "INCLUDE" AND "EXCLUDE"
## all.QC.column.labels<-c("ID","Sample","Recipe","Capture.Method","Lane","Run","target_file","total_reads","total_mapped_reads","total_dup_reads","unmapped_reads","total.bases",paste("on.target.bases.",extend.exons,"bp",sep=""),"on.target.bases", "percent.on.target.bases","median.max.coverage", "median.mean.coverage","mean.mean.coverage","percent.ccds.gt.1","percent.ccds.gt.5","percent.ccds.gt.10","percent.ccds.gt.15","percent.ccds.gt.30","percent_Duplicated","percent_Unmapped","Description","median.insert")
## #############

#############
restrict.analysis<-FALSE
restrict.to.projects<- ""        #AOGC FMHT-N,FTP 
restrict.to.run<-"" # "ALL" # "INCLUDE" # "EXCLUDE"
the.run<-"" # only used with "INCLUDE" AND "EXCLUDE"
################

# ###################################### EXAMPLES OF HOW TO SELECT SPECIFIC PROJECTS
# restrict.analysis<-TRUE
# restrict.to.projects<- c("MODY")         #AOGC FMHT-N,FTP 
# restrict.to.run<-"INCLUDE" # "ALL" # "INCLUDE" # "EXCLUDE"
# the.run<-"" # c("DH10","DH09","00","05") # only used with "INCLUDE" AND "EXCLUDE"
# #############

###################################### EXAMPLES OF HOW TO SELECT SPECIFIC PROJECTS
# restrict.analysis<-TRUE
# restrict.to.projects<- c("CHMND_NIM_EXOME")         #AOGC FMHT-N,FTP 
# restrict.to.run<-"INCLUDE" # "ALL" # "INCLUDE" # "EXCLUDE"
# the.run<-"" # c("DH10","DH09","00","05") # only used with "INCLUDE" AND "EXCLUDE"
# #############

## restrict.analysis<-TRUE
## restrict.to.projects<- c("AMAS","EDSJ","NSLM")         #AOGC FMHT-N,FTP 
## restrict.to.run<-"INCLUDE" # "ALL" # "INCLUDE" # "EXCLUDE"
## the.run<-"" # c("DH10","DH09","00","05") # only used with "INCLUDE" AND "EXCLUDE"
#############

############################################################################
#non.exome.projects<-c("ALSPAC","AS-WGS","PLD","PLD-whole genome","AML-GENOME","bcos_srv","Mowry","QBI-DISC1") # don't look in these folders;
non.exome.projects<-c("AML-PacBio","BTCK","ALSPAC","AS-WGS","PLD","PLD-whole genome","AML-GENOME","bcos_srv","Mowry","QBI-DISC1")
extend.exons<-200
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
combined.on.target.stats<-{} # container for all data
combined.on.target.stats.columns<-{}
####### ip<-1
for(ip in 1:length(projects)){
print(projects[ip])
on.target.stats<-{}  # this will be rebuilt from the individual 
on.target.stats.columns<-{}
on.target.stats.bam<-{}
QC.root.directory<-paste(UQCCG.data,projects[ip],"QC",sep="/")
xx<-try(setwd( QC.root.directory  ),silent=TRUE)
if(inherits(xx, "try-error")){print(paste(projects[ip],": NOT PROCESSED YET",sep=" "));next}


QC.directory<-paste(UQCCG.data,projects[ip],"QC/AlignedQC",sep="/")
xx<-try(setwd( QC.directory  ),silent=TRUE)
if(inherits(xx, "try-error")){print(paste(projects[ip],": NOT PROCESSED YET",sep=" "));next}else{QC.files<-dir(getwd())}
if(length(QC.files)<1){print("No QC files");next}



possible.QC.files<-rep(FALSE,times=length(QC.files))
for(iq in 1:length(QC.extensions)){possible.QC.files<-possible.QC.files |  grepl(paste(QC.extensions[iq],"$",sep=""),QC.files)}
possible.QC.files<-unique(QC.files[possible.QC.files]) # unique in case cat more than on with the QC.extension
if(length(possible.QC.files)<1){print("No possible QC files");next}

## ii<-1
#print(possible.QC.files)
for(ii in 1:length(possible.QC.files)){
  sample.bam.counts<-possible.QC.files[ii]
  true.extension<-{}
  for(iq in 1:length(QC.extensions)){true.extension<-c(true.extension , grepl(paste(QC.extensions[iq],"$",sep=""),sample.bam.counts)) }
   true.extension<-QC.extensions[ true.extension]
  if(length(true.extension)>1){true.extension<-true.extension[which.max(nchar(true.extension))]} # if more than one extension fits choose longes one

  
  QC.file.chk<-read.delim(paste(QC.directory,sample.bam.counts,sep="/"),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
  present<-rownames(QC.file.chk) %in% all.QC.column.labels 
  QC.file.chk<-subset(QC.file.chk,subset=present)  ### Get rid of extra informtion that might  have been collected
  #print(rownames(QC.file.chk))
  if( ("ID" %in% rownames(QC.file.chk)) & ("Capture.Method"  %in% rownames(QC.file.chk)) &  sum(all.QC.column.labels %in% rownames(QC.file.chk))==length(all.QC.column.labels) ){ # this file has been processed and is ok
       if( is.null(on.target.stats) ){ on.target.stats<-as.matrix(QC.file.chk[all.QC.column.labels,])
                                       on.target.stats.columns<- gsub(paste(true.extension,"$",sep=""),"",sample.bam.counts)
                                       on.target.stats.bam<- colnames(QC.file.chk)
                                     }else{
                                       on.target.stats<-cbind(on.target.stats,QC.file.chk[all.QC.column.labels,])
                                       on.target.stats.columns<-c(on.target.stats.columns,gsub(paste(true.extension,"$",sep=""),"",sample.bam.counts))
                                       on.target.stats.bam<- c(on.target.stats.bam,colnames(QC.file.chk))
                                     } # complile on.target.stats
}else{print(paste(sample.bam.counts,"out of date- excluding"))}

} #possible.QC.files
##########################
if(is.null(on.target.stats)){next} # No valid files QC were found so on.target.stats is null skip below
on.target.stats<-rbind(projects[ip],on.target.stats.bam,as.matrix(on.target.stats))

rownames(on.target.stats)<-c("Project","BAM",all.QC.column.labels) # may not be set if combed run and non-run samples
colnames(on.target.stats)<-on.target.stats.columns
xx<-try(write.table(on.target.stats,file=paste(QC.directory,paste(projects[ip],"quality.control.summary",sep="."),sep="/"),col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE))

       if( is.null(combined.on.target.stats) ){ combined.on.target.stats<-on.target.stats
                                                combined.on.target.stats.columns<-on.target.stats.columns
                                     }else{
                                       combined.on.target.stats<-cbind(combined.on.target.stats,on.target.stats)
                                       combined.on.target.stats.columns<-c(combined.on.target.stats.columns,on.target.stats.columns)
                                     } # complile on.target.stats

######################Now compile summaries where there are samples with multiple bam files assume these are in the same project file

} # loop over projects

#print(combined.on.target.stats)
if(!is.null(combined.on.target.stats)){

rownames(combined.on.target.stats)<-c("Project","BAM",all.QC.column.labels)  # may not be set if combed run and non-run samples

colnames(combined.on.target.stats)<-combined.on.target.stats.columns

combined.on.target.stats<-t(combined.on.target.stats)

combined.on.target.stats[1:2,1:10]
#print("HERE")
############ MHAIRI ADD DATA FROM /UQCCG/Data/QC for all samples summary/HiSeq_Runs.csv to combined.on.target.stats here ########

####################################################################################


the.order<-order(combined.on.target.stats[,"Project"],combined.on.target.stats[,"Sample"],combined.on.target.stats[,"percent.ccds.gt.10"])


#combined.on.target.stats[the.order,c(1,4,22)][1:10,]



getwd()
print(paste(combined.on.target.stats.dir,combined.output.file.bam,sep="/"))
write.table(combined.on.target.stats[the.order,],file=paste(combined.on.target.stats.dir,combined.output.file.bam,sep="/"),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

the.order<-order(combined.on.target.stats[,"Sample"],combined.on.target.stats[,"percent.ccds.gt.10"],decreasing = TRUE) # sort in descreasing coverage
combined.on.target.stats<-combined.on.target.stats[the.order,]
combined.on.target.stats[,c(4,22)][1:5,]
combined.on.target.stats<-combined.on.target.stats[!duplicated(combined.on.target.stats[,"Sample"]),]
print(paste(combined.on.target.stats.dir,combined.output.file.sample,sep="/"))

dim(combined.on.target.stats)
write.table(combined.on.target.stats,file=paste(combined.on.target.stats.dir,combined.output.file.sample,sep="/"),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

}else{print("WARNING: NO output written, no valid QC files were detected")}
