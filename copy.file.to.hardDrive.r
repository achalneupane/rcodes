
coverage.data.file<-"/media/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_BAM_Tue_Aug_06_2013.txt"
cov<-read.delim(coverage.data.file,header=TRUE,sep="\t",skip=0,fill=TRUE,stringsAsFactors=FALSE)
clean.aogc<-read.delim("/media/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/clean_AOGC_SEQ_samples.txt",header=FALSE,sep="\t",skip=0,fill=TRUE,stringsAsFactors=FALSE)
number.needed<-200

start.path<-"/media/UQCCG/Sequencing/Projects/AOGC-NGS/BAM"
end.path<-"/media/Seagate Expansion Drive/AOGC_controls"

cov[1:5,]

table(cov[,"Capture.Method"])
table(cov[,"Project"])
Project

wanted<-cov[,"Capture.Method"]=="TruD:TruX" &   cov[,"Project"]=="AOGC-NGS" & cov[,"Sample"] %in% clean.aogc[,1]

sum(wanted)
dim(cov)

cov<-cov[wanted,]

the.order<-order(as.numeric(cov[,"percent.ccds.gt.10"]),decreasing=TRUE)

cov[the.order,][1:10,]

cov<-cov[the.order,][1:number.needed,]

dim(cov)
cov[1:5,]
i<-1
if(!grepl("/$",end.path)){end.path<-paste(end.path,"/",sep="")}
end.path<-gsub(" ","\\ ",end.path,fixed=TRUE)

the.bam<-paste(cov[,"Sample"],"_",cov[,"ID"],".ReCal.sort.bam",sep="")
the.bai<-paste(cov[,"Sample"],"_",cov[,"ID"],".ReCal.sort.bai",sep="")

###### check
the.files<-dir(start.path)
sum(!(the.bam %in% the.files)) ## zero if all present
sum(!(the.bai %in% the.files)) ## zero if all present

for(i in 1:length(the.bam)){
  print(paste("Doing ",i," of ", length(the.bam)))
  system( paste("cp",paste(start.path,the.bai[i],sep="/"),end.path,sep=" " ) )
  system( paste("cp",paste(start.path,the.bam[i],sep="/"),end.path,sep=" " ) )
}


