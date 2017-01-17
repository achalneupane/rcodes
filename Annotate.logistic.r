


###########################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
plot.dir<-"/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/ForearmFrLoTrauma"
assoc.file<-"ForearmFrLoTrauma.logistic.covar.assoc.logistic" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
bim.file<-"GBS_wtccc_commom.bim" ### bim file - will extract out snps used
ann.file
#################################################





core.ann<-c("chr","start","end","REF","ALT","TYPE")
code.dir<-"/media/Bioinform-D/Research/AML sequencing"
genome.build<-"hg19"
target<-"temp"



ann<-read.delim("/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/zCall_AOGC.with.evoker_corrected_clean_FINAL.analysis.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

colnames(ann)


write.table(assoc.ann,file="annotated_exome_chip.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)






##################################
setwd(code.dir)
source("ucsc.table.names.r")   # load in the UCSC tables these use the db file names and not their lable-names 
source("annotate_SNPs_subroutines.r")
source("ucsc.table.names.processor.r")
setwd(plot.dir)

library("Biostrings")
library("IRanges")
if(genome.build=="hg19"){library("BSgenome.Hsapiens.UCSC.hg19"); have.chr<-seqnames(Hsapiens);chr.lengths<-seqlengths(Hsapiens)}
if(genome.build=="hg18"){library("BSgenome.Hsapiens.UCSC.hg18"); have.chr<-seqnames(Hsapiens);chr.lengths<-seqlengths(Hsapiens)}
if(genome.build=="mm9"){library("BSgenome.Mmusculus.UCSC.mm9"); have.chr<-seqnames(Mmusculus)}
if(genome.build=="mm10"){library("BSgenome.Mmusculus.UCSC.mm10"); have.chr<-seqnames(Mmusculus)}
######## CAREFUL this code used by read.plink.assoc.files.r


sample.files<-bim.file
names(sample.files)<-"snp"
sample.grs<-paste(names(sample.files),"grs",sep=".")
snp.dir<-plot.dir
setwd(code.dir)
source("read.plink.bim.files.r")
indels<-snp.grs


if(!grepl("^chr",indels[1,"chr"])){indels[,"chr"]<-paste("chr",indels[,"chr"],sep="")}

#######################################################################

assoc<-read.delim(paste(plot.dir,assoc.file,sep="/"),header=T,sep="",fill=TRUE,stringsAsFactors=FALSE)
## bim<-read.delim(paste(plot.dir,bim.file,sep="/"),header=F,sep="",fill=TRUE,stringsAsFactors=FALSE)
## bim[500:550,]
## colnames(fam)<-c("chr","SNP","cm","POS","SEX","status"

assoc[1:5,]
keep<-assoc[,"TEST"]=="ADD"
assoc<-assoc[keep,]



posns<-match(assoc[,"SNP"],indels[,"SNP"])
missing<-is.na(posns)
sum(missing)
indels<-indels[posns[!missing],]

dim(assoc)
dim(indels)

no.chr<-check.chr.and.positions(indels,chr.lengths)

sum(no.chr)
assoc<-assoc[!no.chr,]
indels<-indels[!no.chr,]
dim(assoc)
dim(indels)
setwd(plot.dir)
write.table(assoc,file="final.assoc",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

assoc[1:50,]
##################################


#annotate.GWAS<-function(indels,genome.build,geneanno.DB,filter.DB)

## geneanno.DB<-c("refGene")  # returns 2 extra columns c("refGene","knownGene","ensGene")
## names(geneanno.DB)<-c("refGene") # c("refGene","knownGene","ensGene")

## filter.DB<-c("snp137")  # returns 2  extra columns (DB,score) c("snp132","1000g2011may_all")
## names(filter.DB)<-c("ID")

             
write.table(indels[,core.ann],file=target,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

#table_annovar.pl ex1.human humandb/ -protocol refGene,phastConsElements44way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp135,avsift,ljb_all -operation g,r,r,f,f,f,f,f -nastring NA

system(paste("table_annovar.pl ",target,"  ",anno.DB.location," -buildver ",genome.build," -protocol refGene,snp137 -operation g,f --outfile ",paste(target,"txt",sep=".")," ",sep="" )  )
#system(paste("table_annovar.pl ",target,"  ",anno.DB.location," -buildver ",genome.build," -protocol refGene,snp137 -operation g,f --outfile ",paste(target,"txt",sep=".")," ",sep="" )  )


### get location of header:
## input<-paste(target,"txt",paste(genome.build,"multianno",sep="_"),"txt",sep=".")
## column.labels<-read.table(input,header=F,nrows=1,skip=0,fill=TRUE,stringsAsFactors=FALSE)
## num.vars<-dim(column.labels)[2]
## skip.lines<-1


## ann<-try(scan(input,what=character(num.vars),skip=skip.lines,fill=TRUE))
## num.lines<-length(ann)/(num.vars)
## dim(ann)<-c(num.vars,num.lines)
## ann<-t(ann)
## colnames(ann)<-column.labels 
############### ADD IF HERE IF H


ann<-read.delim(paste(target,"txt",paste(genome.build,"multianno",sep="_"),"txt",sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
ann[ann==""]<-NA  ## Haploviews does not like ""
dim(ann)
dim(assoc)
dim(indels)
dim(ann)
ann.key<-build.key(ann,c("Chr","Start","End","Ref","Alt"))
indel.key<-build.key(indels,c("chr","start","end","REF","ALT"))

posns<-match(indel.key,ann.key)
missing<-is.na(posns)
sum(!missing)
assoc.ann<-cbind(assoc,ann[posns,])
remove.cols<-c("Chr")
assoc.ann<-assoc.ann[,!(colnames(assoc.ann) %in% remove.cols)]


assoc.ann[1:5,]
ann[1:5,]
setwd(plot.dir)
assoc.ann[,"AAChange.refGene"]<-gsub(":",".",assoc.ann[,"AAChange.refGene"])
assoc.ann[,"ExonicFunc.refGene"]<-gsub(" ",".",assoc.ann[,"ExonicFunc.refGene"]) # does not like spaces

write.table(assoc.ann,file="final_annotated.assoc",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

write.table(assoc.ann[1:10,],file="final_annotated2.assoc",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
[,c(1:16)]
/media/UQCCG-Analysis/Myb chipseq
