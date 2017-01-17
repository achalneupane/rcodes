

#################################### OPTIONS AND PATHS NEEDED
options(show.error.messages = TRUE)

UQCCG.data<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes" ### ROOT directory for projects
the.sep<-"/" # "\\" for windows



directories<-list.dirs(UQCCG.data,recursive =FALSE)
directories
##

i<-3
for (i in 1 : length(directories)){

  sub.directories<-list.dirs(directories[i],recursive =FALSE)
sub.directories

  j<-1
for (j in 1 : length(sub.directories)){

  snp.dir<-paste(sub.directories[j],"SNPs",sep=the.sep)

xx<-try(setwd( snp.dir ),silent=TRUE)
if(inherits(xx, "try-error")){print("ERROR Could not find SNP dir"); next}

  the.files<-dir(getwd())
  vcf.files<-the.files[grep(".vcf$",the.files)]
  if(length(vcf.files)<1){next})
  
  
if(inherits(xx, "try-error")){print("ERROR Could not find SNP dir")}



  chromo1<-try(scan(vcf.files[1],what=character(),n=5000,sep="\n",skip=0,fill=TRUE)) ## find the start of the vcf file
skip.lines<-grep("^#CHROM",chromo1)
if(length(skip.lines)>1){print("ERROR multiple chrom lables found");skip.lines<-skip.lines[1]}
skip.lines<-skip.lines
 
options(show.error.messages = TRUE)
column.labels<-read.delim(vcf.files[1],header=F,nrows=1,skip=(skip.lines-1),sep="\t",fill=TRUE,stringsAsFactors=FALSE)


column.labels


  

} # j loop

} # i loop
