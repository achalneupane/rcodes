
library(GenomicFeatures)
library(HardyWeinberg)
library(Biostrings)

code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir) # load in the UCSC tables these use the db file names and not their lable-names 
source("annotate_SNPs_subroutines.r")
options("width"=250,"max.print"=1000)
core.ann<-c("chr","start","end","REF","ALT","TYPE")



############################################### Seqence data file:
project.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/2013-10-27_AML_with_AOGCControl_NoFailedLane.TGCM-AML.All-maf-filtered.txt"
###########################################################

##################################################
#setwd(analysis.dir)

## grep("Gene.Name",a.indel[,16])
## save(list=c("column.labels"),file="column.labels.RData")
## load("column.labels.RData")
################## fast read ###########
column.labels<-read.delim(project.file,header=F,nrows=1,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
num.vars<-dim(column.labels)[2]
a.indel<-scan(project.file,what=character(num.vars),skip=1,sep="\t",fill=TRUE,na.strings="",quote="\"")
num.lines<-length(a.indel)/(num.vars)
dim(a.indel)<-c(num.vars,num.lines)
a.indel<-t(a.indel)
colnames(a.indel)<-column.labels
########################################
the.chr<-unique(a.indel[,"chr"])
print(paste("Doing Chromosome ",the.chr))

if(sum(!grepl("^chr",the.chr))>0){
a.indel[,"chr"]<-paste("chr",a.indel[,"chr"],sep="")
}

a.indel[1:5,"chr"]
key<-build.key(a.indel,core.ann)




#################################### got some missing gene names still.
## all.genes[grep("GTPBP4",names(all.genes))]
no.gene<-is.na(a.indel[,"Gene.Names"]) | a.indel[,"Gene.Names"]=="NA"
all.genes<-sort(table(a.indel[,"Gene.Names"]),decreasing=TRUE)

ens.to.hgnc<-read.delim("/media/UQCCG/Sequencing/Data/Genomes/hg19/ENSG_to_HGNC.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
ens.to.hgnc[1:5,]
genes<-unique(a.indel[no.gene,"Gene.Embl"])

posns<-match(genes,ens.to.hgnc[,"Ensembl.Gene.ID"])
missing<-is.na(posns)

ens.to.hgnc<-ens.to.hgnc[posns[!missing],]
dups<-duplicated(ens.to.hgnc[,"Ensembl.Gene.ID"])
ens.to.hgnc<-ens.to.hgnc[!dups,]
ens.to.hgnc[1:5,]

posns<-match(a.indel[no.gene,"Gene.Embl"],ens.to.hgnc[,"Ensembl.Gene.ID"])
missing<-is.na(posns)
sum(missing)
a.indel[no.gene,"Gene.Names"]<-ens.to.hgnc[posns,"HGNC.symbol"]

a.indel[no.gene,"Gene.Names"][1:10]
ens.to.hgnc[posns,"HGNC.symbol"][1:10]

no.gene<-is.na(a.indel[,"Gene.Names"]) | a.indel[,"Gene.Names"]=="NA"
a.indel[no.gene,"Gene.Names"]<-a.indel[no.gene,"Gene.Embl"]
a.dash<-a.indel[,"Gene.Names"]=="-"
sum(a.dash)
a.indel[a.dash,][1:5,1:50]
a.indel[a.dash,"Gene.Names"]<-a.indel[a.dash,"Feature.Embl"]

a.dash<-a.indel[,"Gene.Names"]=="-"
sum(a.dash)
a.indel[a.dash,][1:5,1:50]

final.names<-a.indel[a.dash,"knownGene::gene"]
final.names<-gsub("\\(dist=\\d*\\)","",final.names)
final.names<-gsub("\\(dist=NONE\\)","",final.names)
final.names[grepl("^NO_ANNOVAR",final.names)]<-"-"

a.indel[a.dash,"Gene.Names"]<-final.names
a.dash<-a.indel[,"Gene.Names"]=="-"
sum(a.dash)
a.indel[a.dash,c("chr","start")]
final.names<-build.key(a.indel[a.dash,],c("chr","start"))
a.indel[a.dash,"Gene.Names"]<-final.names

unannotated.hits<-a.dash

all.genes<-sort(table(a.indel[,"Gene.Names"]),decreasing=TRUE)
all.genes[1:30]


#################################################################################################################

 #



############################## real in the "sample Sheet" to identify samples :

the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/BAM/TGCM-AML-combine_SampleSheet.csv"

sample.sheet.full<-read.delim(the.sample.sheet,header=T,sep=",",fill=TRUE,stringsAsFactors=FALSE)
sample.sheet.full[1:5,]
colnames(sample.sheet.full)
dim(sample.sheet.full)
aml.samples<-sample.sheet.full[sample.sheet.full[,"AffectionStatus"]==2,"ParticipantCode"]
normal.samples<-sample.sheet.full[sample.sheet.full[,"AffectionStatus"]==1,"ParticipantCode"]
all.possible.samples<-sample.sheet.full[,"ParticipantCode"]


######################### Get all NMP1 mutations by gene NAME:
target<-grepl("NPM1",a.indel[,"Gene.Names"])
sum(target)
a.indel[target,c(1:16,16,30)]
a.indel[target,c(1:6,9,10,11,28,29,30,31,32,33)]

###########################################


################# GET by location/region using Granges
## FLT3
## the.chr<-"chr13"
## begin<-28608214
## end<-28608336


## NPM1 location range
the.chr<-"chr5"
begin<-170837522
end<-170837676

#########


a.indel[1:5,"chr"]
a.indel.gr<-GRanges(seqnames =a.indel[,"chr"],ranges = IRanges(start=as.numeric(a.indel[,"start"]),end=as.numeric(a.indel[,"end"])),strand="+")

indel.gr<-GRanges(seqnames =the.chr,ranges = IRanges(start=as.numeric(begin),end=as.numeric(end)),strand="+")

common<-overlapsAny(a.indel.gr,indel.gr)
sum(common)
a.indel[common ,1:6]


a.indel[common & !is.flat.geno,the.samples.HC.GT][1:2,]
a.indel[common & !is.flat.geno,the.samples.HC.GT][,1]
a.indel[common & !is.flat.geno,the.samples.HC.GQ]

###################################################################


### which samples: have mutstions

the.samples<-aml.samples
all.possible.samples


the.samples.HC<-all.possible.samples[all.possible.samples %in% the.samples]
the.samples.HC.GT<-paste(the.samples.HC,"GT",sep=".")
the.samples.HC.AD<-paste(the.samples.HC,"AD",sep=".")
is.flat.geno<-grepl(":flat$",a.indel[,"TYPE"])
the.samples.HC.GT
the.samples.HC.AD

a.indel[common & !is.flat.geno, the.samples.HC.GT]

hits<-apply(a.indel[common & !is.flat.geno,the.samples.HC.GT],2,function(x) sum(x!="0/0"))
sum(hits)
hits
#x<-a.indel[common & !is.flat.geno,c(the.samples.HC.GT,the.samples.HC.GQ)][1,]


muts.in.cases<-apply(a.indel[common & !is.flat.geno,the.samples.HC.GT],1,function(x) { paste(names(x)[x!="0/0" & !is.na(x)],collapse=",")})
muts.in.cases ## samples with muststions


figure<- match(loci,key)


########################################################
check<-16

quality.cases<-rep("",times=length(loci))


a.indel.sub<-a.indel[common & !is.flat.geno,]
loci<-key[common & !is.flat.geno]
muts.in.cases<-gsub(".GT","",muts.in.cases)
figure<- match(loci,key)

for(check in 1:length(loci)){
print(check)
#check<-"chr11:130066457:130066457:-:A:indel"
# posn<-grep(loci[check],key)
posn<-check

if(muts.in.cases[check]!=""){
#the.gt<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GT",sep=".")
#the.gq<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GQ",sep=".")
the.dp<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"AD",sep=".")

quality.cases[check]<-paste(a.indel.sub[posn,the.dp],collapse=";")

#a.indel[posn,the.gq]
#a.indel[posn,the.gt]
a.indel[posn,the.dp]
}


} # end check
##########################################################################
quality.cases
annotations<-a.indel[common & !is.flat.geno,c(1:6,9,10,11,28,29,30,31,32,33)]
out<-cbind(annotations,muts.in.cases,quality.cases)
out
getwd()
write.table(out,file="NMP1_mutations_in_cases.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)








