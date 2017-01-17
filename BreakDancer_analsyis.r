options(width=150)

#####################################################################################################
######################################################start here
#####################################################################################################

library(GenomicFeatures)



## setwd("/media/scratch/AML sequence/AML_Exome_2/Orginal_BAMS_run_2/breakdancer") # full run
## skip.lines<-c(2,3,3)
## file.names<-c("control_old.breakmax","aml_old.breakmax","relapse_old.breakmax")
#########################################################AML WHOLE GENOME  set up read in for sigle files #####
file.dir<-"/media/scratch/AML sequence/AML_Genome"  #Now /media/scratch/AML sequence/AML_Exome_2/Orginal_BAMS_run_2/Final_Bams
bam.dir<-"/media/scratch/AML sequence/AML_Genome" #Now /media/scratch/AML sequence/AML_Exome_2/Orginal_BAMS_run_2/Final_Bams
ann.dir<-"/media/Bioinform-D/Research/annovar/humandb"
sample.names<-c("Normal","Diagnosis","Relapse") ##$$$ NOrmal is AWAYS in the FIRST posn
sample.files<-c("control_breakdancer.txt","aml_breakdancer.txt","relapse_breakdancer.txt") #sample.files<-c("control_old.breakmax","aml_old.breakmax","relapse_old.breakmax")
breakdancer.version<-"new"
skip.lines<-c(13,11,11)
sample.bams<-c("cleaned_genome_hg18_control_sorted.bam","cleaned_genome_hg18_aml_sorted.bam","cleaned_genome_hg18_relapse_sorted.bam")
names(sample.bams)<-sample.names[1:length(sample.bams)] # Normal MUST BE FIRST and in same order as 'sample.names"
sample.grs<-paste(sample.names,".gr",sep="")
## column.labels<-c("chr","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT");
the.library<-list()
z.thresh.for.reads<-2
genome.build<-"hg18"
dbSNP.build<-"130"
#########################################################data.names<-c("control","aml","relapse")




#########################################################AML EXOME set up read in for sigle files #####
file.dir<-"/media/scratch/AML_sequence/AML_Exome_2/Final_Bams"  #Now /media/scratch/AML sequence/AML_Exome_2/Orginal_BAMS_run_2/Final_Bams
bam.dir<-"/media/scratch/AML_sequence/AML_Exome_2/Final_Bams" #Now /media/scratch/AML sequence/AML_Exome_2/Orginal_BAMS_run_2/Final_Bams
ann.dir<-"/media/Bioinform-D/Research/annovar/humandb"
sample.names<-c("Normal","Diagnosis","Relapse") ##$$$ NOrmal is AWAYS in the FIRST posn
sample.files<-c("control_breakmax_2","aml_breakmax_2","relapse_breakmax_2") #sample.files<-c("control_old.breakmax","aml_old.breakmax","relapse_old.breakmax")
breakdancer.version<-"old"
skip.lines<-c(2,3,3)
sample.bams<-c("realign_cleaned_control.sort.bam","realign_cleaned_aml.sort.bam","realign_cleaned_relapse.sort.bam")
names(sample.bams)<-sample.names[1:length(sample.bams)] # Normal MUST BE FIRST and in same order as 'sample.names"
sample.grs<-paste(sample.names,".gr",sep="")
## column.labels<-c("chr","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT");
the.library<-list()
z.thresh.for.reads<-2
genome.build<-"hg19"
dbSNP.build<-"131"
#########################################################data.names<-c("control","aml","relapse")


## setwd("/media/scratch/AML sequence/AML_Exome")
## sample.files<-c("exome_breakmax_control","exome_breakmax_diagnosis","exome_breakmax_relapse")
## skip.lines<-c(1,1,1)

setwd(file.dir)
dir(getwd())
library(Rsamtools)

for(i in 1:length(sample.files)){


  ### the.libaray is all the read groups that are in the bam file (may have different means and stds
lib<-read.delim(sample.files[i],header=F,nrow=skip.lines[i],sep="\t",fill=TRUE,stringsAsFactors=FALSE)
lib<-lib[grep("mean:",lib[,2]),] 
lib.means<-as.numeric(gsub("mean:","",lib[,2]))
lib.stds<-as.numeric(gsub("std:","",lib[,3]))
if(breakdancer.version=="old"){lib.readlens<-as.numeric(gsub("readlen:","",lib[,6]))}else{lib.readlens<-as.numeric(gsub("readlen:","",lib[,5]))}
if(breakdancer.version=="old"){lib.names<-gsub("#","recalibrated_",lib[,1]);lib.names<-gsub(".sort.bam","",lib.names)}else{
  lib.names<-gsub("library:","",lib[,6]) ## new version gets libray
  ## if(sample.files[i]=="control_breakdancer.txt"){lib.names<-paste
  }

names(lib.means)<-lib.names
names(lib.stds)<-lib.names
names(lib.readlens)<-lib.names


#############################
#### the.library contails library names but need the ID's as these label the raeds
bam.header<-scanBamHeader(sample.bams[i]) # sample.bams and samples.files in the same order
bam.read.groups<-( bam.header[[1]]$text[ names(bam.header[[1]]$text)=="@RG" ] )


bam.read.groups<-lapply(bam.read.groups,function(x){
       line<-strsplit(x,split=":")
       unlist(lapply(line,function(y) {value<-y[2]; names(value)<-as.character(y[1]); value}))
     })
bam.read.IDs<-unlist(lapply(bam.read.groups,function(x) x["ID"]))
names(bam.read.IDs)<-unlist(lapply(bam.read.groups,function(x) x["LB"]))

lib.means<-lib.means[names(bam.read.IDs)]
lib.stds<-lib.stds[names(bam.read.IDs)]
lib.readlens<-lib.readlens[names(bam.read.IDs)]

names(lib.means)<-(bam.read.IDs)
names(lib.stds)<-(bam.read.IDs)
names(lib.readlens)<-(bam.read.IDs)
####### did this way incase some ID's have the same libray name (should not happen) 

the.library[[sample.grs[i]]]<-list(means=lib.means,stds=lib.stds,readlens=lib.readlens)
#####################



#### the.library contails libray names but need the ID's as these label the raeds

## data<-eval(as.name(sample.names[i]))
data<-read.delim(sample.files[i],header=T,skip=skip.lines[i],sep="\t",fill=TRUE,stringsAsFactors=FALSE)
start=as.numeric(data[,"Pos1"])
end=as.numeric(data[,"Pos2"])
dim(data[start>end,])  ### some instances where order is wrong

strand<-rep("+",times=length(start))  
strand[start>end]<-"-"
start.use<-start
end.use<-end


################## flip and set strands:
start.use[start>end]<-end[start>end]
end.use[start>end]<-start[start>end]

if(!("Allele_frequency" %in% colnames(data))){data[,"Allele_frequency"]<-NA}

data.gr<-GRanges(seqnames =data[,"X.Chr1"],ranges = IRanges(start=as.numeric(start.use),end=as.numeric(end.use)),
                 chr2=data[,"Chr2"],
                 Orientation1=data[,"Orientation1"],
                 Orientation2=data[,"Orientation2"],
                 Type=data[,"Type"],
                 Size=data[,"Size"],
                 Score=as.numeric(data[,"Score"]),
                 num_Reads=as.numeric(data[,"num_Reads"]),
                 num_Reads_lib=data[,"num_Reads_lib"],
                 Allele_frequency=data[,"Allele_frequency"],
                 strand="+")

assign(sample.grs[i],value=data.gr)
## data.gr
print(paste("Done",sample.grs[i],sep=" "))

}

names(the.library)<-sample.grs

## special case for aml below
## names(the.library[["Normal.gr"]]$means)<-paste(names(the.library[["Normal.gr"]]$means),"_A",sep="")
## names(the.library[["Normal.gr"]]$stds)<-paste(names(the.library[["Normal.gr"]]$stds),"_A",sep="")
## names(the.library[["Normal.gr"]]$readlens)<-paste(names(the.library[["Normal.gr"]]$readlens),"_A",sep="")

Normal.gr
Diagnosis.gr
Relapse.gr

#####################################
####################################

######################################
######################################
########## Quality thresholds

#### how they are distrbuted
aden<-density(values(Normal.gr)[,"Score"])
aden<-density(values(Diagnosis.gr)[,"Score"])
aden<-density(values(Relapse.gr)[,"Score"])
plot(aden,new=FALSE)
abline(v=95,col="blue",lwd=2)
abline(v=90,col="red",lwd=2)
abline(v=80,col="orange",lwd=2)
abline(v=70,col="green",lwd=2)

####################### Note comparions done before quality score consider case where filter  for quality first and then do comparisons
########################################################################################################## AML

################# AML EXOME##
quality.names<-c("Score","num_Reads","Allele_frequency")  #
quality.type<-c("numeric","numeric","numeric")
quality.relapse<-c(99,4,0.05) ### remove the "." to keep the SNPs intact
quality.tumor<-c(99,4,0.05)
quality.normal<-c(80,2,0.01)  ## factor: "test in equal to"  number: "is greater than"

names(quality.relapse)<-quality.names
names(quality.tumor)<-quality.names
names(quality.normal)<-quality.names

################# AML GENOME##
quality.names<-c("Score","num_Reads")  #
quality.type<-c("numeric","numeric")
quality.relapse<-c(99,4) ### remove the "." to keep the SNPs intact
quality.tumor<-c(99,4)
quality.normal<-c(80,2)  ## factor: "test in equal to"  number: "is greater than"

names(quality.relapse)<-quality.names
names(quality.tumor)<-quality.names
names(quality.normal)<-quality.names

##### names must match sample.names
quality.cut<-list()
quality.class<-list()
sample.grs
quality.cut[sample.grs[1]]<-list(quality.normal)
quality.cut[sample.grs[2]]<-list(quality.tumor)
quality.cut[sample.grs[3]]<-list(quality.relapse)

quality.class[sample.grs[1]]<-list(quality.type)
quality.class[sample.grs[2]]<-list(quality.type)
quality.class[sample.grs[3]]<-list(quality.type)

##########################################################################################################

##########################################################################################################

quality.cut
quality.class
sample.grs
quality.cut[sample.grs] ## this must work else loop below will fail
quality.class[sample.grs] ## this must work else loop below will fail

quality.thresh<-list()
k<-0
for(i in 1:length(sample.grs)){
print(sample.grs[i])
data<-eval(as.name(sample.grs[i]))
quality.names<- names(quality.cut[[i]])
quality.type<-quality.class[[i]]
 for(j in 1:length(quality.names)){
 k<-k+1
 if(quality.type[j]=="factor"){quality.thresh[[k]] <- values(data)[,quality.names[j]]==quality.cut[[sample.grs[i]]][quality.names[j]]}
 if(quality.type[j]=="numeric"){quality.thresh[[k]] <- as.numeric(as.character(values(data)[,quality.names[j] ] )) >= as.numeric(quality.cut[[sample.grs[i]]][quality.names[j]]) }

 print(paste("passing filter",quality.names[j],":",sum(quality.thresh[[k]]),sep=" "))
names(quality.thresh)[k]<-paste(sample.grs[i],quality.names[j],sep="::")
}
}
## quality.thresh



##################################### Threshold the normal sample:
names(quality.thresh)

########################################## CHOOSE ONE ####################################################################
#########AML exome
normal.filter <- quality.thresh$"Normal.gr::Score" & quality.thresh$"Normal.gr::num_Reads" & quality.thresh$"Normal.gr::Allele_frequency"
sum(normal.filter)
Normal.filtered.gr<-Normal.gr[normal.filter]

#########AML genome
normal.filter <- quality.thresh$"Normal.gr::Score" & quality.thresh$"Normal.gr::num_Reads"
sum(normal.filter)
Normal.filtered.gr<-Normal.gr[normal.filter]
####################################### make a matrix of cross comparisons based on the filter criteria ####################


########################################
####################################### make a matrix of cross comparisons based on the filter criteria
## sample.names<-c("Normal","Diagnosis","Relapse","All")
sample.names<-c(sample.names,"Normal.filtered")
sample.grs<-paste(sample.names,".gr",sep="")
sample.grs

####################################### make a matrix of cross comparisons
sample.grs.mat<-rep(sample.grs,times=length(sample.grs))
dim(sample.grs.mat)<-c(length(sample.grs),length(sample.grs))
sample.grs.mat<-t(sample.grs.mat)
sample.grs.mat<-apply(sample.grs.mat,2,function(x) paste(sample.grs,x,sep="::"))
sample.grs.mat



comparisons<-list()
k<-0
##everthing but diaginal
for(i in 1:length(sample.grs)){   ## replace sample.grs by 
  for(j in 1:length(sample.grs)){
    if(i==j){next}
    k<-k+1
    target<-sample.grs.mat[i,j]
    print(paste(i,j,target,sep=":"))
    query.subject<-unlist(strsplit(target,split="::"))
    query<-eval(as.name(query.subject[1]))
    subject<-eval(as.name(query.subject[2])) ## assuming that control is listed first
    test<-query %in% subject ##test is the length of query a.screen.cuts[[the.screen]]$low.cut[the.score]
    comparisons[[k]]<-test
    names(comparisons)[k]<-target
  }}
 ## comparisons 
names(comparisons)


######################################
########### Construst the hit classes you want from the comparisons avalibale:
names(comparisons)
names(quality.thresh)

########AML exome
tumor.snp<-!comparisons$"Relapse.gr::Normal.filtered.gr" & quality.thresh$"Relapse.gr::Score" & quality.thresh$"Relapse.gr::num_Reads" & quality.thresh$"Relapse.gr::Allele_frequency"
sum(tumor.snp)

tumor.snp<-!comparisons$"Diagnosis.gr::Normal.filtered.gr" & quality.thresh$"Diagnosis.gr::Score" & quality.thresh$"Diagnosis.gr::num_Reads" & quality.thresh$"Diagnosis.gr::Allele_frequency"
sum(tumor.snp)
##############

########AML genome
tumor.snp<-!comparisons$"Relapse.gr::Normal.filtered.gr" & quality.thresh$"Relapse.gr::Score" & quality.thresh$"Relapse.gr::num_Reads" 
sum(tumor.snp)

tumor.snp<-!comparisons$"Diagnosis.gr::Normal.filtered.gr" & quality.thresh$"Diagnosis.gr::Score" & quality.thresh$"Diagnosis.gr::num_Reads"
sum(tumor.snp)


## relapse.hits<-!comparisons$"relapse::control" & quality.thresh$"relapse" & !comparisons$"relapse::aml" ## aml.hits<-!comparisons$"aml::control" & quality.thresh$"aml" &  !comparisons$"aml::relapse" 
## relapse.aml.hits<-comparisons$"relapse::aml" & !comparisons$"relapse::control" & quality.thresh$"relapse"
## aml.relapse.hits<-comparisons$"aml::relapse" & !comparisons$"aml::control" & quality.thresh$"aml"
## the.comparisons<-c("relapse.hits","aml.hits","relapse.aml.hits","aml.relapse.hits")
## for(i in 1:length(the.comparisons)){
## data<-eval(as.name(the.comparisons[i]))
## sum(is.na(data))
## print(paste("For",the.comparisons[i],"number passing filter: ",sum(data),sep=" "))
## }



################################## CHOOSE ONE
indels<-as.data.frame(Relapse.gr[tumor.snp]) ### convert genomic ranges to a dta frame
target<-"relpase.breakdancer.txt"

indels<-as.data.frame(Diagnosis.gr[tumor.snp]) ### convert genomic ranges to a dta frame
target<-"aml.breakdancer.txt"

indels[indels[,"seqnames"]=="chr13",]#  28608234

indels[indels[,"seqnames"]=="chr13",]


plot(indels[,"Size"],indels[,"Allele_frequency"])
plot(indels[,"Size"],indels[,"Allele_frequency"],xlim=c(0,200),ylim=c(0,0.1))
 chr13  28608234  28608363   130      + chr13       26+26-       26+26-  ITX  -69    99         8  s2_1|3:s1_1|3:s2_2|2             0.08
#########################################################
################################ Dump can can use annovar
chr13  28608234  28608363


###############  consider CTX issues! transchromosomal indels ###
indels[indels[,"Type"]=="CTX",]
indels[as.character(indels[,"seqnames"])!=as.character(indels[,"chr2"]),]
if(length(indels[indels[,"Type"]=="CTX","Type"])>0){
  indels[,"True end"]<-indels[,"end"]
indels[indels[,"Type"]=="CTX","end"]<-indels[indels[,"Type"]=="CTX","start"]+  median(indels[,"width"],na.rm=TRUE)
}
####################################

core.ann<-c("seqnames","start","end")

indel.ann<-as.data.frame(indels)
indel.ann[1:5,]
dim(indel.ann)
indel.ann.ori<-indel.ann

indel.ann[,"seqnames"]<- gsub("chr","",indel.ann[,"seqnames"])
Ref<-rep(0,times=dim(indel.ann[1])[1])
Obs<-rep("-",times=dim(indel.ann[1])[1])
## Obs[indel.ann[,"Type"]=="DEL"]<-"-" ## ..9695176  
## 20003944  	0  	-  	comments: a 342kb deletion encompassing GJB6
## Ref[indel.ann[,"Type"]=="INS"]<-"-"

indel.ann<-cbind(indel.ann[,colnames(indel.ann) %in% core.ann],Ref,Obs,indel.ann[,!(colnames(indel.ann) %in% core.ann)])
col.ann<-colnames(indel.ann)
indel.ann[1:5,]
write.table(indel.ann,file=target,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

genome.build
dbSNP.build
#### skip step one so have miRNA and mnon protein coding regions
ann.dbtype<-"refGene"
system(paste("annotate_variation.pl -geneanno  -buildver ",genome.build," -dbtype refGene ",target," /media/Bioinform-D/Research/annovar/humandb/ --outfile ",paste(target,".step1",sep=""),sep="") )
system(paste("auto_annovar.pl -buildver ",genome.build," -verdbsnp ",dbSNP.build,"  -dbtype refGene --nocheckfile -step 2-3,8-9 ",target," /media/Bioinform-D/Research/annovar/humandb/",sep=""))
system(paste("annotate_variation.pl -geneanno  -buildver ",genome.build," -dbtype ",ann.dbtype," ",paste(target,".step8.varlist",sep="")," /media/Bioinform-D/Research/annovar/humandb/",sep=""))
## using ensGene gets the noncoding transcripts
system(paste("annotate_variation.pl -regionanno   -buildver ",genome.build," -dbtype dgv ",paste(target,".step8.varlist",sep="")," /media/Bioinform-D/Research/annovar/humandb/",sep=""))

## system(paste("annotate_variation.pl -regionanno   -buildver ",genome.build," -dbtype snp130 ",paste(target,".step8.varlist",sep="")," /media/Bioinform-D/Research/annovar/humandb/",sep=""))
######################################
######################################

known.str<-read.delim(paste(target,".step8.varlist.hg19_dgv",sep=""),header=F,skip=0,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
indel.ann.exon<-read.delim(paste(target,".step8.varlist.exonic_variant_function",sep=""),header=F,skip=0,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
indel.ann<-read.delim(paste(target,".step8.varlist.variant_function",sep=""),header=F,skip=0,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
colnames(indel.ann.exon)<-c("line","coding type","exon",col.ann)
colnames(indel.ann)<-c("coding type","gene",col.ann)
colnames(known.str)<-c("DB","Name",col.ann)

rownames(known.str)<-paste(paste("chr",known.str[,"seqnames"],sep=""),paste(known.str[,"start"],end=known.str[,"end"],sep="-"),sep=":")
rownames(indel.ann.exon)<-paste(paste("chr",indel.ann.exon[,"seqnames"],sep=""),paste(indel.ann.exon[,"start"],end=indel.ann.exon[,"end"],sep="-"),sep=":")
rownames(indel.ann)<-paste(paste("chr",indel.ann[,"seqnames"],sep=""),paste(indel.ann[,"start"],end=indel.ann[,"end"],sep="-"),sep=":")

## known.str[,"seqnames"]<-paste("chr",known.str[,"seqnames"],sep="")
## indel.ann.exon[,"seqnames"]<-paste("chr",indel.ann.exon[,"seqnames"],sep="")
## indel.ann[,"seqnames"]<-paste("chr",indel.ann[,"seqnames"],sep="")

colnames(known.str)[colnames(known.str)=="seqnames"]<-"chr"
colnames(indel.ann.exon)[colnames(indel.ann.exon)=="seqnames"]<-"chr"
colnames(indel.ann)[colnames(indel.ann)=="seqnames"]<-"chr"


known.str[1:5,]
indel.ann[1:5,]
indel.ann.exon[1:5,]
indel.ann.exon[,"line"]<- as.integer(gsub("line","",indel.ann.exon[,"line"]))

indel.ann[,"exon"]<-NA
indel.ann[indel.ann.exon[,"line"],"exon"]<-indel.ann.exon[,"exon"]
indel.ann[rownames(known.str),1:5]

known.str[,"Name"]<-gsub("Score=0;","",known.str[,"Name"])

dim(indel.ann)
dim(known.str)

tx.names<-unlist(strsplit(indel.ann[,"exon"],split=":"))
tx.names<-tx.names[grepl("ENST",tx.names)]
gene.names<-gsub("\\(dist=",",",indel.ann[,"gene"])
gene.names<-unlist(strsplit(gene.names,split=","))
gene.names<-gene.names[!grepl("\\)$",gene.names)]
gene.names



indel.ann[indel.ann[,"Type"]=="CTX",]
################## Build the summary file

allele.labels<-c("P-value","Allele Freq","Major Allele Count","Minor Allele Count","Missing Mates","Coverage")
the.sample.order<-names(sample.bams)

sample.class<-as.character( t(matrix(data=rep(the.sample.order,times=length(allele.labels)),nrow=length(the.sample.order),ncol=length(allele.labels),byrow=FALSE)) )
sample.class
sample.class<-paste(sample.class,allele.labels,sep=".")
dim(sample.class)<-c(length(allele.labels),length(the.sample.order))
sample.class<-as.character(t(sample.class))

allele.mat<-matrix(data=NA,ncol=length(sample.class),nrow=dim(indel.ann)[1])
colnames(allele.mat)<-sample.class
allele.mat[1:2,]
indel.ann[1:5,]


summary<-cbind(indel.ann[,1:8],indel.ann[,"exon"],known.str[rownames(indel.ann),"Name"],allele.mat,NA,NA,NA,NA,indel.ann[,9:(dim(indel.ann)[2]-1)],NA,NA,NA,NA)
colnames(summary)<-c(colnames(indel.ann)[1:8],"exon","Known Varient",colnames(allele.mat),"Known Cancer Gene","Known in AML","Known in Haematopoietic Tissue","Description",colnames(indel.ann)[9:(dim(indel.ann)[2]-1)],"COSMIC-Primary Hist","COSMIC-Hist SubType","COSMIC-Primary Site","COSMIC-Combinations")
summary[1:2,]               

########################################################################################################






################################### annotate with gene :                       
library(biomaRt)
listMarts()
## mart <- useMart("ensembl")
## listDatasets(mart)
if(genome.build=="hg18"){mart=useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="may2009.archive.ensembl.org",path="/biomart/martservice",archive=FALSE)}else{
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")}
## listFilters(mart) gene.names<-summary[,"Gene Symbol"]  #### change to hg18 here 
## mart=useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="may2009.archive.ensembl.org",path="/biomart/martservice",archive=FALSE)

gene.names<-unique(gene.names)
gene.names
if(ann.dbtype=="ensGene"){
  ann.refseq<-getBM(attributes = c( "ensembl_gene_id","hgnc_symbol","gene_biotype","entrezgene","description"), filters = "ensembl_gene_id", values=gene.names, mart = mart)}else
{
ann.refseq<-getBM(attributes = c( "ensembl_gene_id","hgnc_symbol","gene_biotype","entrezgene","description"), filters = "hgnc_symbol", values=gene.names, mart = mart)
}

ann.refseq

gene.names<-gsub("\\(dist=",",",summary[,"gene"])
gene.names<-unlist(strsplit(gene.names,split=","))
gene.names<-gene.names[!grepl("\\)$",gene.names)]
gene.names

#grep(paste("^ann.refseq[,"hgnc_symbol"]summary[,"gene"],ann.refseq[,"hgnc_symbol"])


           
posns<-match(summary[,"gene"],ann.refseq[,"hgnc_symbol"])
missing<-is.na(posns)
sum(missing)
posns
                       
summary[!missing,"Description"]<-ann.refseq[posns[!missing],"description"]
summary[1,]

summary[summary[,"chr"]=="13",]
##############################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
################### Alleel frequency via Scan Bam:
options(width=150) ## output length
library(Rsamtools)
library(GenomicRanges)
data<-summary
core.ann<-c( "chr","start","end","strand") #'elementMetadata' cannot use "seqnames", "ranges", "strand", "seqlengths", "start", "end", "width", or "element" as column names
data.gr<-GRanges(seqnames =paste("chr",data[,"chr"],sep=""),ranges = IRanges(start=as.numeric(data[,"start"]),end=as.numeric(data[,"end"])),strand=data[,"strand"])

reserved.gr.names<-c("seqnames", "ranges", "strand", "seqlengths", "isCircular", "start", "end", "width", "element")
names.have<-colnames(data)[!(colnames(data) %in% core.ann)]
names.have[names.have %in% reserved.gr.names]<-paste(names.have[names.have %in% reserved.gr.names],1,sep=".")
colnames(data)[!(colnames(data) %in% core.ann)]<-names.have

remove.cols.from.values<-grepl("COSMIC",colnames(data))
## colnames(data) %in% c("COSMIC")
values(data.gr)<-data[,(!(colnames(data) %in% core.ann)) & !remove.cols.from.values ]
##### grange must have "chr" in seqnames

## data.gr<-GRanges(seqnames =paste("chr",13,sep=""),ranges = IRanges(start=as.numeric(28608234),end=as.numeric(28608363)),strand="+")
which<-  data.gr
which
setwd(bam.dir)

files<-dir(getwd())
files[grepl("bam",files)]
## files[match(bam.file,files)]  ### bam file must exist


## scanBamWhat()
## scanBamFlag()
## params<-ScanBamParam(which=which,flag=scanBamFlag(isPaired=TRUE,isUnmappedQuery=TRUE,hasUnmappedMate=FALSE,isFirstMateRead=NA,isSecondMateRead=NA),what=c("qname","flag","rname","seq","strand","pos","mpos","qwidth","cigar","qual","mapq","isize", "mrnm" ),tag="RG" ) ## this one for missing mates

## params<-ScanBamParam(which=which,flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=NA,isValidVendorRead=TRUE),simpleCigar = FALSE,reverseComplement = FALSE,what=c("qname","flag","rname","seq","strand","pos","mpos","qwidth","cigar","qual","mapq","isize", "mrnm" ),tag="RG" )  ### NOTE isValidVendorRead=FALSE shoudl be TRUE

std.expand<-6 ## to get allele frequency look out to this range
######################## Phred scores DECODED are ascii 0-93
str<-as.character(PhredQuality(0:93))
str<-unlist(strsplit(str,split = character(0)))
score<-0:93
names(score)<-str
## score
## the.qual<-unlist(strsplit(as.character(aln1[[1]]$qual[7]),split = character(0)))
## the.qual

## score[the.qual]  # provides the quality
## library("ChIPsim") ### needed for decodeQuality
## -10*log10(decodeQuality(aln1[[1]]$qual[7], type = c("Sanger")))

#####################################################
for(ii in 1:length(sample.bams)){
  
  expand<-as.integer(max(the.library[[sample.grs[ii]]]$means,na.rm=TRUE) + std.expand*max(the.library[[sample.grs[ii]]]$stds,na.rm=TRUE))
  which.expanded<-which+expand

  params<-ScanBamParam(which=which.expanded,flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=NA,isValidVendorRead=TRUE),simpleCigar = FALSE,reverseComplement = FALSE,what=c("qname","flag","rname","seq","strand","pos","mpos","qwidth","cigar","qual","mapq","isize", "mrnm" ),tag="RG" )
  
aln1 <- scanBam(sample.bams[ii],param=params)
## aln1 lapply(aln1,function(x) length(x$seq))
regions<-names(aln1)
regions<-strsplit(gsub("-",":",regions),split=":")
## regions
print(paste("Doing",sample.bams[ii],sep=" "))
print(paste("Nunmber regions:",length(regions),sep=" "))

base.qual.filter<-15 # bases with a Phred less than 20 are ignored 20-99% 10-90% ok  15 appears to be about what GATK is using
read.qual.filter<-20 # mapQ takes into account the paired read GAYK perhaps about 20

missing.mates.count<-{}
major.allele.count<-{}
minor.allele.count<-{}
alleles.freq<-{}
coverage<-{}


for (i in 1:length(regions)){
   ## loc.wanted<-as.numeric(regions[[i]][3]) # loc.wanted<-101401645
   ## posns<-loc.wanted-aln1[[i]]$pos+1  ##  ## posns<-loc.wanted-aln1[[i]]$mpos+1

   poor.read<- aln1[[i]]$mapq <= read.qual.filter # sort(aln1[[i]]$mapq,decreasing=TRUE)
   ######################################## map to location in read
   ## posns
   ## test<-cigarOpTable(aln1[[i]]$cigar)  ## as.character(aln1[[i]]$seq)  "ILLUMINA-A52086_200000:1:94:3590:14054#0"
   ## print(paste("test no deletions:",i,":",sum(apply(test[,c("I","D","N","H","P")],2,sum))),sep="")

   inserts<- aln1[[i]]$isize
   qname<- aln1[[i]]$qname
   read.groups<- aln1[[i]]$tag$RG
   library.means<-the.library[[sample.grs[[ii]] ]]$means[read.groups]
   library.stds<-the.library[[sample.grs[[ii]] ]]$stds[read.groups]
   
   ## missing.mates<-is.na(inserts)
   ## signs<-inserts/abs(inserts)
   z.score<-(abs(inserts)-library.means)/library.stds
   names(z.score)<-qname
   z.score<-z.score[!poor.read]
   
   ## sort(z.score[!missing.mates])  ## [1:5]
   ## sum(missing.mates)
   sort(z.score[!missing.mates])
   unique.z.score<-z.score[match(unique(qname),names(z.score))] ## the qnames are identical for the mates , some mates may be OUTside the ROI however , missing mates are single entries
   missing.mates<-is.na(unique.z.score)
   
   sum(missing.mates)
   sum( (abs(  unique.z.score) > 2),na.rm=TRUE )
  ## as.numeric(sort(unique.z.score[!missing.mates]))  ## [1:5]
   
   coverage[i]<-length(unique.z.score)
   minor.allele.count[i]<- sum( (abs(  unique.z.score) > z.thresh.for.reads)   | missing.mates )
   major.allele.count[i]<- coverage[i]- minor.allele.count[i]
   alleles.freq[i]<- minor.allele.count[i]/(coverage[i]) ##minor allele frequency
   missing.mates.count[i]<-sum(missing.mates)
    minor.allele.count[i]
    alleles.freq[i]
   coverage[i]
   ## print(toString(paste(unlist(labels(counts)),":",counts," (",mean.quals,")",sep="")))
   
}



## key.gr<-data.gr+expand
## summary.key<-paste(seqnames(key.gr),":",start(key.gr),"-",end(key.gr),sep="")
summary.key<-paste("chr",summary[,"chr"],":",as.numeric(summary[,"start"])-expand,"-",as.numeric(summary[,"end"])+expand,sep="")

posns<-match(summary.key,names(aln1))
missing<-is.na(posns)
sum(missing)
posns
    
summary[!missing,paste(names(sample.bams)[ii],"Allele Freq",sep=".")]<-round(alleles.freq[posns[!missing]],digits=3)
summary[!missing,paste(names(sample.bams)[ii],"Major Allele Count",sep=".")]<-major.allele.count[posns[!missing]]
summary[!missing,paste(names(sample.bams)[ii],"Minor Allele Count",sep=".")]<-minor.allele.count[posns[!missing]]
summary[!missing,paste(names(sample.bams)[ii],"Missing Mates",sep=".")]<-missing.mates.count[posns[!missing]]
summary[!missing,paste(names(sample.bams)[ii],"Coverage",sep=".")]<-coverage[posns[!missing]]

    }

## summary[!missing,"Normal Coverage"]<-coverage[posns[!missing]]
summary[1:10,]
summary[summary[,"chr"]=="13",]

summ<-summary # summary<-summ
#######################fix zero counts and missing alleles
#### Major/MInor Allele can be missing cause of no sequence reads NA
## Or after filtering reading the BAm there are no alleles "o" will result  minor missing cause is homozygous

for( i in 1:length(names(sample.bams))){

### fix cases where th allele frequency is undefiles
 problem.freq<- !is.finite(summary[,paste(names(sample.bams)[i],"Allele Freq",sep=".")])
 summary[problem.freq,paste(names(sample.bams)[i],"Allele Freq",sep=".")]<-0

  problem.freq<- !is.finite(summary[,paste(names(sample.bams)[i],"Coverage",sep=".")])
 summary[problem.freq,paste(names(sample.bams)[i],"Coverage",sep=".")]<-0
 ## summary[problem.freq,]
 
 ###check
 print(paste(paste("Problem Frequencies:",names(sample.bams)[i],sum(problem.freq),sep=" ")))
}


################ check now that all the major and minor alleles are the same acrocc samples

######################################################## Model the Cancer / Normal contrbution
tumor.contamination<-c(0,0.3,0.15)# 30%
names(tumor.contamination)<-names(sample.bams)  ## sample.names<-c("Normal","Diagnosis","Relapse")

tumor.contamination<-c(0,0.3)# 30%
names(tumor.contamination)<-names(sample.bams)  

#### Major/MInor Allele can be missing cause of no sequence reads NA
## Or after filtering reading the BAm there are no alleles "o" will result  minor missing cause is homozygous
for( i in 2:length(names(sample.bams))){


 ##############if there are no reads then just use the REF and ALT to set:
 
  sum(!is.finite(summary[,paste(names(sample.bams)[i],"Allele Freq",sep=".")]))
 
 
 ##30% of reads will be normal subtract this expected amount first
 normal.allele.freq<-summary[,paste(sample.names[1],"Allele Freq",sep=".")]

 
 normal.coverge.in.tumor<-(summary[,paste(sample.names[i],"Coverage",sep=".")]*tumor.contamination[sample.names[i]])
 normal.minor.in.tumor<-normal.coverge.in.tumor*normal.allele.freq
 normal.major.in.tumor<- normal.coverge.in.tumor-normal.minor.in.tumor
 normal.major.in.tumor+ normal.minor.in.tumor
 
 true.sample.cov<-summary[,paste(sample.names[i],"Coverage",sep=".")]- normal.coverge.in.tumor
 true.sample.major.counts<-summary[,paste(names(sample.bams)[i],"Major Allele Count",sep=".")] - normal.major.in.tumor
 true.sample.minor.counts<-summary[,paste(names(sample.bams)[i],"Minor Allele Count",sep=".")] - normal.minor.in.tumor
  true.sample.major.counts+ true.sample.minor.counts

 true.sample.freq<-true.sample.minor.counts/ true.sample.cov
 true.sample.freq[!is.finite(true.sample.freq)]<-0
 
 ## summary[,paste(sample.names[3],"Allele Freq",sep=".")]
 sort(normal.allele.freq)
 normal.allele.freq[normal.allele.freq<0.01]<-0.01
 true.sample.cov[true.sample.cov==0]<-1
## Z=X-np/sqrt(np(1-p)) z sore of the minor allele frequency:
 z= (true.sample.minor.counts- true.sample.cov*normal.allele.freq) / (sqrt(true.sample.cov*normal.allele.freq*(1-normal.allele.freq)))
 p<-signif(2*(pnorm(z, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)),digits=3)        ##this is the 2 tailed
 p
 
 summary[,paste(sample.names[i],"P-value",sep=".")]<-p
 ## summary[alleles.differ,paste(sample.names[i],"P-value",sep=".")]<-0
                                }

 summary[,paste(sample.names[1],"P-value",sep=".")]<-1
 summary[,paste(sample.names[2],"P-value",sep=".")]
summary[,paste(sample.names[3],"P-value",sep=".")]
 ## summary[alleles.differ,]
########################### Z score to P-value ##########################
##  ## erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
##  erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE) ### 1-erf:  2 tailed sigma=1 68% p=0.32; sigma=2 95% p=0.05;  sigma=3 99.7% p=0.003;  = 2*erfc(x)
## z<-0
## erfc(z/sqrt(2))
##  ### this is 2 tailed (eg. higher and lower than 2 sd if x=2 )
## 1-pnorm(x, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE) ### this is 1 tailed (eg. higher  2 sd if x=2 )
## need to abs(z)
## x<-1:6
## signif(2*(pnorm(x, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)),digits=3) ## is the 2 tailed one to USE!

#########################################################################


colnames(summary)
summary[1:2,]
summary[,"gene"]
colnames(indel.ann)[colnames(indel.ann)=="seqnames"]<-"chr"

gene.names<-gsub("\\(dist=",",",indel.ann[,"gene"])
gene.names<-(strsplit(gene.names,split=","))
gene.names<-unlist(lapply(gene.names,function(x) x[1]))
summary[,"Gene Symbol"]<-gene.names
gene.names

################# Annotate with Known cancer gene
setwd(ann.dir)
load("cosmic.known.cancer.genes.RData")
cosmic.known.cancer.genes[1:5,]

posns<-match(summary[,"Gene Symbol"],cosmic.known.cancer.genes[,"gene"])

sort(summary[,"Gene Symbol"])
cosmic.known.cancer.genes[grepl("@",cosmic.known.cancer.genes[,"gene"]),"gene"] ## gene groups not handled
gene.groups<-cosmic.known.cancer.genes[grepl("@",cosmic.known.cancer.genes[,"gene"]),"gene"]
gene.groups<-gsub("@","",gene.groups)
apply(as.matrix(gene.groups),1,function(x) grep(paste("^",x,sep=""),summary[,"Gene Symbol"]))

missing<-is.na(posns)
sum(missing)
posns
cosmic.known.cancer.genes[posns[!missing],]

summary[!missing,"Known Cancer Gene"]<-paste(cosmic.known.cancer.genes[posns[!missing],"gene"],
                                             " (",cosmic.known.cancer.genes[posns[!missing],"cancer.type"],")",
                                             " [",cosmic.known.cancer.genes[posns[!missing],"cancer.class"],"]"
                                             ,sep="")
 summary[,"Known Cancer Gene"]

################# Annotate Mardis AML
setwd("/media/Bioinform-D/Research/annovar/humandb")
load("mardis aml.RData")
aml[1:5,]

posns<-match(summary[,"Gene Symbol"],aml[,"gene"])

missing<-is.na(posns)
sum(missing)
posns
aml[posns[!missing],]

summary[!missing,"Known in AML"]<-paste(aml[posns[!missing],"gene"]," [",aml[posns[!missing],"num.observed"],"]"  ,sep="")
 summary[,"Known in AML"]
################# Annotate with hemopoetic tissue
setwd("/media/Bioinform-D/Research/annovar/humandb")
load("hemo.tissue.RData")
hemo.tissue[1:5,]

posns<-match(summary[,"Gene Symbol"],hemo.tissue[,"gene"])

missing<-is.na(posns)
sum(missing)
posns
hemo.tissue[posns[!missing],]

summary[!missing,"Known in Haematopoietic Tissue"]<-paste(hemo.tissue[posns[!missing],"gene"]," [",hemo.tissue[posns[!missing],"num.observed"],"]"  ,sep="")
 summary[,"Known in Haematopoietic Tissue"]

############################ Filter based on    P-value
the.order<-order(summary[,paste(names(sample.bams)[2],"P-value",sep=".")])
summary<-summary[the.order,]


summ<-summary
min.p.value<-0.005

not.significant<-  (summary[,paste(names(sample.bams)[3],"P-value",sep=".")] > min.p.value)

## summary[not.significant,]
sum(not.significant)
sum(!not.significant)
summary<-summary[!not.significant,]

summary[1:5,]
summary[,1:16]
##########################################################################

#################### COSMIC
setwd(ann.dir)
load("COMIC_genes.RData")


gene.places<-apply(as.matrix(gene.names),1,function(x) grep(paste("^",x,"$",sep=""),cosmic.genes[,"Gene.name"],perl=TRUE)  )
cosmic.ann<-lapply(gene.places,function(x) toString(unique(paste(cosmic.genes[x,"Primary.histology"],cosmic.genes[x,"Histology.subtype"],cosmic.genes[x,"Primary.site"],sep="::"))))
length(cosmic.ann)
cosmic.ann[1:5]
cosmic.ann<-unlist(cosmic.ann)


cosmic.ann.prim<-lapply(gene.places,function(x) toString(unique(cosmic.genes[x,"Primary.histology"])))
cosmic.ann.prim<-unlist(cosmic.ann.prim)

cosmic.ann.hist<-lapply(gene.places,function(x) toString(unique(cosmic.genes[x,"Histology.subtype"])))
cosmic.ann.hist<-unlist(cosmic.ann.hist)

cosmic.ann.site<-lapply(gene.places,function(x) toString(unique(cosmic.genes[x,"Primary.site"])))
cosmic.ann.site<-unlist(cosmic.ann.site)


length(cosmic.ann)
cosmic.ann[1:5]
cosmic.ann.site[1:5]
cosmic.ann.hist[1:5]
cosmic.ann.prim[1:5]


summary[,"COSMIC-Primary Hist"]<-cosmic.ann.prim
summary[,"COSMIC-Hist SubType"]<-cosmic.ann.hist
summary[,"COSMIC-Primary Site"]<-cosmic.ann.site
summary[,"COSMIC-Combinations"]<-cosmic.ann
summary[1:2,]
##########################################

setwd(bam.dir)

summary[,"Known Varient"]<-known.str[rownames(indel.ann),"Name"]
posns<-match(summary[,"Gene Symbol"],summary.aml[,"Gene Symbol"])
missing<-is.na(posns)
summary[!missing,"In Diagnosis"]<-"YES"
summary[is.na(summary[,"In Diagnosis"]),"In Relapse"]<-"NO"

write.table(summary,file="summary.LARGE.INDELs.Replase.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## save.image("large indel relapse.RData")
summary.relapse<-summary
## save(list=c("summary.relapse"),file="summary.relapse.RData") # load("summary.relapse.RData");summary<-summary.relapse


posns<-match(summary[,"Gene Symbol"],summary.relapse[,"Gene Symbol"])
missing<-is.na(posns)
summary[!missing,"In Relapse"]<-"YES"
summary[is.na(summary[,"In Relapse"]),"In Relapse"]<-"NO"

  
write.table(summary,file="summary.LARGE.INDELs.AML.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## save.image("large indel aml.RData")
summary.aml<-summary
## save(list=c("summary.aml"),file="summary.aml.RData")   # load("summary.aml.RData");summary<-summary.aml


