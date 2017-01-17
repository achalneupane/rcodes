######## ONLY NEED TO CHOOSE A DIRECTORY AND EXTENSIONS - used tab delimited files 


#source("http://bioconductor.org/biocLite.R")
# biocLite(c("HardyWeinberg"))n
# install.packages("HardyWeinberg")


###############################################
#analysis.dir<-"/media/ga-apps/UQCCG/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/Analysis"
#annotate.dir<-"/media/ga-apps/UQCCG/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/Annotate"

analysis.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-08-21_MND_AOGC_Tulane/Analysis"
annotate.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-08-21_MND_AOGC_Tulane/Annotate"



project.extension<-"_analysis-maf-filtered.txt.small.RData"
project.name<-"2014-08-21_MND_AOGC_Tulane" ## prefix for output file
fam<-c("ALL.ALL_GENOTYPES") #  ALL or  c() ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/BAM/TGCM-AML-combine_SampleSheet.csv"
the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-08-21_MND_AOGC_Tulane/Analysis/MND_AffStatus_using.csv"
bad.samples.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-08-21_MND_AOGC_Tulane/Analysis/bad.samples.csv"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)
remove.cols<-c()

#regions.file<-"/media/scratch2/AOGC-NGS/GFOS/gefos.seq/METHODS/0613-skatmeta-gefos/static/Homo_sapiens.GRCh37.70.protein_coding.genespace_boundaries.5k.split100k.txt"
core.ann<-c("chr","start","end","REF","ALT","TYPE") # out put to annanlsys programs and need foe colun labels
dont.build.summary<-FALSE ##

GATK.SB<-TRUE
maf.threshold.filter.to.use<-c(0.05)
a.label<-"NMD_maf_filtered"
###########################################################################

###############################################
#analysis.dir<-"/media/ga-apps/UQCCG/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/Analysis"
#annotate.dir<-"/media/ga-apps/UQCCG/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/Annotate"

analysis.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/Data_Freeze_ChNMD_AOGC_Tulane/Analysis"
annotate.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/Data_Freeze_ChNMD_AOGC_Tulane/Annotate"

#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/Data_Freeze_ChNMD_AOGC_Tulane/Analysis/Data_Freeze_ChNMD_AOGC_Tulane.chr1.ALL.ALL_GENOTYPES_analysis.txt

project.extension<-"_analysis-maf-filtered.txt.small.RData"
project.name<-"Data_Freeze_ChNMD_AOGC_Tulane" ## prefix for output file
fam<-c("ALL.ALL_GENOTYPES") #  ALL or  c() ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/BAM/TGCM-AML-combine_SampleSheet.csv"
the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-08-21_MND_AOGC_Tulane/Analysis/MND_AffStatus_using.csv"
bad.samples.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-08-21_MND_AOGC_Tulane/Analysis/bad.samples.csv"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)
remove.cols<-c()

#regions.file<-"/media/scratch2/AOGC-NGS/GFOS/gefos.seq/METHODS/0613-skatmeta-gefos/static/Homo_sapiens.GRCh37.70.protein_coding.genespace_boundaries.5k.split100k.txt"
core.ann<-c("chr","start","end","REF","ALT","TYPE") # out put to annanlsys programs and need foe colun labels
dont.build.summary<-FALSE ##

GATK.SB<-TRUE
maf.threshold.filter.to.use<-c(0.01)
a.label<-"NMD_maf_filtered"
###########################################################################


###############################################
#analysis.dir<-"/media/ga-apps/UQCCG/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/Analysis"
#annotate.dir<-"/media/ga-apps/UQCCG/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/Annotate"

analysis.dir<-"/media/old-scratch/media/scratch2/AOGC-NGS/Analysis"
annotate.dir<-"/media/old-scratch/media/scratch2/AOGC-NGS/Annotate"
#/media/old-scratch/media/scratch2/AOGC-NGS/Analysis/AOGC-Genotyping.output.chr9.AOGC_ALL.analysis-maf-filtered.txt.BMD_EFF_STD_HIP.coding.small.RData


project.extension<-".small_final.RData"
project.name<-"AOGC-Genotyping.output." ## prefix for output file
fam<-c(".AOGC_ALL.analysis-maf-filtered.txt.BMD_EFF_STD_HIP.coding",".AOGC_ALL.analysis-maf-filtered.txt.BMD_EFF_STD_HIP.non.coding",".AOGC_ALL.analysis-maf-filtered.txt.BMD_EFF_STD_HIP.all") #  ALL or  c() ""-one project (the prefix of the summary files to collect
fam<-c(".AOGC_ALL.analysis-maf-filtered.txt.BMD_EFF_STD_FN.coding",".AOGC_ALL.analysis-maf-filtered.txt.BMD_EFF_STD_FN.all") 
#fam<-c(".AOGC_ALL.analysis-maf-filtered.txt.BMD_EFF_STD_HIP") 
#the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/BAM/TGCM-AML-combine_SampleSheet.csv"

a.label<-"AOGC_coding_FN"
###########################################################################


###############################################
## /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-01-15_LeoPharma_Dec2014Freeze/Analysis/2015-01-15_LeoPharma_Dec2014Freeze.chrX.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt.SampleProject.coding.somatic.small_final.RData

analysis.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-01-15_LeoPharma_Dec2014Freeze/Analysis"
annotate.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-01-15_LeoPharma_Dec2014Freeze/Annotate"
#/media/old-scratch/media/scratch2/AOGC-NGS/Analysis/AOGC-Genotyping.output.chr9.AOGC_ALL.analysis-maf-filtered.txt.BMD_EFF_STD_HIP.coding.small.RData


project.extension<-"small_final.RData"
project.name<-"2015-01-15_LeoPharma_Dec2014Freeze." ## prefix for output file
#fam<-c(".ALL.ALL_GENOTYPES_analysis-maf-filtered.txt.SampleProject.coding.somatic.",".ALL.ALL_GENOTYPES_analysis-maf-filtered.txt.SampleProject.coding.",".ALL.ALL_GENOTYPES_analysis-maf-filtered.txt.SampleProject.all.somatic.") #  ALL or  c() ""-one project (the prefix of the summary files to collect
fam<-c(".ALL.ALL_GENOTYPES_analysis-maf-filtered.txt.SampleProject.coding.somatic.") 

a.label<-"Leo_SSC"
###########################################################################



options(width=250,max.print=5000)

code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")
source("hwe.r")



library(skatMeta)  ## ridge regression
#library(SKAT) ## skat method
library(GenomicFeatures)
library(HardyWeinberg)
library(Biostrings)




  
#### test fam list
files<-dir(analysis.dir)
the.extension<-paste(project.extension,"$",sep="")
files<-files[grepl(the.extension ,files)]
toString( unique(unlist( mapply(function(x){x[length(x)]}, strsplit(gsub(the.extension,"",files),split=".",fixed=TRUE)   ))) )
toString( files)
files
####


#############################################################################################################
              
######################################### Predefined variables required
##################################################################################
#### assume has format project.chr.fam.extension or chr.project.fam.extension
setwd(analysis.dir)
getwd()

files<-dir(analysis.dir)
the.extension<-paste(project.extension,"$",sep="")
files<-files[grepl(the.extension ,files)]
if(length(fam)==1 & (fam=="ALL" | fam=="All" | fam=="all")){
  fam<-unique(unlist( mapply(function(x){x[length(x)]}, strsplit(gsub(the.extension,"",files),split=".",fixed=TRUE)   )))
}

fam

#
ifam<-1
for(ifam in 1:length(fam)){

the.extension<-paste(fam[ifam],project.extension,"$",sep="")
project.files<-files[grepl(the.extension ,files)]
print(sort(paste("Doing: ",project.files,sep=""))) # project.files<-project.files[1:22]

indels<-{}
the.col<-{}
project.files
#

#meta.results.skat

## meta.results.skatO.all<-meta.results.skatO.all
## meta.results.burden.all<-meta.results.burden.all
## pheno.use.all<-pheno.use.all
## snpinfo.all<-snpinfo.all
## genotypes.all<-genotypes.all
## pass.all<-pass.all
## high.missing.all<-high.missing.all
## annotations.all<-annotations.all
## help.all<-help.all
## key.all<-key.all
## summary.geno.extra.all<-summary.geno.extra.all
project.files
ichr<-1

for(ichr in 1:length(project.files)){ ### loop over chromosomes

print(project.files[ichr])
load(project.files[ichr])
meta.results.burden[1:50,]
 
## meta.results.burden[meta.results.burden[,"gene"] %in% clusters[,1],]
## meta.results.burden[meta.results.burden[,"gene"] %in% clusters[,2],]
cols.annotations<-c(1:6,16,28,7,30,34,37:42,43,14,32,33)
#colnames(a.indel)[cols.annotations]
annotations<-a.indel[,cols.annotations]

## the.samples.use<-pheno[,"SAMPLE"]
## the.samples.use<-paste(the.samples.use,".GT",sep="")
## genotypes<-a.indel[pass,the.samples.use]
## print("start QC")
## genotypes[genotypes=="NA"]<-NA
## genotypes[genotypes=="0/0"]<-0
## genotypes[genotypes=="0/1"]<-1
## genotypes[genotypes=="1/1"]<-2

## ########### prevent any averaging
## dim(genotypes)
## genotypes[is.na(genotypes)]<-0
## dim(genotypes)
## num.col<-dim(genotypes)[2]
## num.row<-dim(genotypes)[1]
## ## genotypes[1:5,1:20]
## genotypes<-as.numeric(as.matrix(genotypes))
## dim(genotypes)<-c(num.row,num.col)
## genotypes<-t(genotypes)


## c("case.control","snpinfo.ori","formula","clusters","pheno.types","ipheno","clusters.wanted","genotypes","p","meta.results.skat","meta.results.skatO","meta.results.burden","pheno","target.pheno.col","snpinfo","fil.genotypes","pass","high.missing.table","a.indel","help","key","summary.geno.extra","full.qual","bad.effect","maf.filter","in.common.hit.gene","on.x.y","unannotated.hits","not.flat.genotype","are.repeats","are.in.repeats","ok.missing","hw.controls.ok.filt","no.genotypes","rare.in.Control","rare.in.Control.filt","in.any.normal","in.any.normal.filt","are.in.repeats.back","are.in.repeats.forward","all.genes")

if(ichr==1){
meta.results.skatO.all<-meta.results.skatO
meta.results.burden.all<-meta.results.burden
meta.results.skat.all<-meta.results.skat
#pheno.use.all<-pheno.use
snpinfo.all<-snpinfo
snpinfo.ori.all<-snpinfo.ori
genotypes.all<-genotypes
fil.genotypes.all<-fil.genotypes
pass.all<-pass
high.missing.all<-high.missing.table
annotations.all<-annotations
help.all<-help
key.all<-key
summary.geno.extra.all<-summary.geno.extra
in.any.normal.all<-in.any.normal
in.any.normal.filt.all<-in.any.normal.filt


}else{
meta.results.skatO.all<-rbind(meta.results.skatO.all,meta.results.skatO)
meta.results.burden.all<-rbind(meta.results.burden.all,meta.results.burden)
meta.results.skat.all<-rbind(meta.results.skat.all,meta.results.skat)

#pheno.use.all<-rbind(pheno.use.all,pheno.use
snpinfo.all<-rbind(snpinfo.all,snpinfo)
genotypes.all<-cbind(genotypes.all,genotypes)
fil.genotypes.all<-rbind(fil.genotypes.all,fil.genotypes)
snpinfo.ori.all<-rbind(snpinfo.ori.all,snpinfo.ori)
pass.all<-c(pass.all,pass)
high.missing.all<-rbind(high.missing.all,high.missing.table)
annotations.all<-rbind(annotations.all,annotations)
help.all<-rbind(help.all,help)
key.all<-c(key.all,key)
in.any.normal.all<-c(in.any.normal.all,in.any.normal)
in.any.normal.filt.all<-c(in.any.normal.filt.all,in.any.normal.filt)
summary.geno.extra.all<-rbind(summary.geno.extra.all,summary.geno.extra)
}



} # loop over projects



#} # loop over Pheno

meta.results.skatO<-meta.results.skatO.all
meta.results.burden<-meta.results.burden.all
meta.results.skat<-meta.results.skat.all
#pheno.use<-pheno.use.all
snpinfo<-snpinfo.all
genotypes<-genotypes.all
fil.genotypes<-fil.genotypes.all
pass<-pass.all
high.missing<-high.missing.all
annotations<-annotations.all
help<-help.all
key<-key.all
in.any.normal<-in.any.normal.all
in.any.normal.filt<-in.any.normal.filt.all

summary.geno.extra<-summary.geno.extra.all

snpinfo.raw<-snpinfo

fam.name<-fam[ifam]
fam.name<-gsub(".ALL.ALL_GENOTYPES_analysis-maf-filtered.txt.SampleProject.","",fam.name)
fam.name<-gsub(".AOGC_ALL.analysis-maf-filtered.txt.","",fam.name)
snap.file<-paste(a.label,fam.name,sep=".")
snap.file
image.file<-paste(snap.file,"Full_COLLECTED","RData",sep=".")
image.file
save.image(file=snap.file)


the.order<-     order(meta.results.burden[,"p"])
sum(is.na(meta.results.burden[,"p"])) ## bad p-values shoudl not happen
meta.results.burden<-  meta.results.burden[the.order,]

meta.results.burden[1:5,]
meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]

the.order<-     order(meta.results.skat[,"p"])
meta.results.skat<-  meta.results.skat[the.order,]
meta.results.skat[1:5,]

the.order<-     order(meta.results.skatO[,"p"])
sum(is.na(meta.results.skatO[,"p"])) ## bad p-values shoudl not happen
meta.results.skatO<-  meta.results.skatO[the.order,]
meta.results.skatO[1:5,]

the.genes<-meta.results.burden[,"gene"]
the.genes<-the.genes[!(the.genes %in% clusters.wanted)]
the.genes<-unique(the.genes)



posns<-match(the.genes,meta.results.burden[,"gene"])
missing<-is.na(posns)
sum(missing)
data<-meta.results.burden[posns,]


the.tests<-c("burden","skat","skatO")
#colnames(data)[colnames(data)=="p"]<-"P_burden"
colnames(data)<-paste(colnames(data),"burden",sep="_")
data[1:5,]


posns<-match(the.genes,meta.results.skatO[,"gene"])
missing<-is.na(posns)
sum(missing)
extra<-meta.results.skatO[posns,c("gene","p","pmin","rho")]
colnames(extra)<-paste(colnames(extra),"skatO",sep="_")
extra[1:5,]
data<-cbind(extra,data)
data[1:5,]

posns<-match(the.genes,meta.results.skat[,"gene"])
missing<-is.na(posns)
sum(missing)
extra<-meta.results.skat[posns,c("gene","p","Qmeta")]
colnames(extra)<-paste(colnames(extra),"skat",sep="_")
extra[1:5,]
data<-cbind(extra,data)
data[1:5,]

order.wanted<-c("gene_burden",paste("p",the.tests,sep="_"))
others<-colnames(data)[!(colnames(data) %in% order.wanted)]
data<-data[,c(order.wanted,others)]
data[1:5,]
colnames(data)[colnames(data)=="gene_burden"]<-"gene"
data[1:5,]
getwd()
write.table(data,file=paste(snap.file,"assoc_summary","txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

} # loop over fam


################## need to recalcalcolate clusters wanted in done per chromosome

################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
##################### RELOAD########################
library(skatMeta)  ## ridge regression
#library(SKAT) ## skat method
library(GenomicFeatures)
library(HardyWeinberg)
library(Biostrings)

## analysis.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis"
## setwd(analysis.dir)
## getwd()

## snap.file<-"coding.0.01.all.geno.all.filters_no.imput"
## snap.file<-"coding.0.01.all.geno.all.filters.NEW"
## snap.file<-"coding.0.001.all.geno.all.filters.NEW"

#load(paste(snap.file,"RData",sep="."))
#options(width=200)
meta.results.burden[1:50,]
meta.results.skat[1:5,]
meta.results.skatO[1:5,]



## target.pheno.col<-fam.name
## formula<-paste(target.pheno.col,"~",paste(covars,collapse="+"),sep="")
print(formula)
formula<-formula(formula)
pheno[1:5,]

###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
## test<-c("MYCBP2","SLC25A24","TMCO3","C13orf35")
## test<-c("MYCBP2","SLC25A24","TMCO3","C13orf35")
## test<-c("chr22:41252435-41252687:ST13")
meta.results.burden[meta.results.burden[,"gene"] %in% clusters[,"SCC_Genes"]  ,]
posns<-match(clusters[,"SCC_Genes"],meta.results.burden[,"gene"])
cbind(posns,clusters[,"SCC_Genes"])

clusters

top<-350
to.unwind.name<-"TOP_350"

to.unwind<-c(meta.results.burden[1:top,"gene"],meta.results.skatO[1:top,"gene"])
## to.unwind<-c("ST14") #, "MCM7", "RNPC3")
## to.unwind<-c("FANC_complex.all") # to.unwind<-meta.results.burden[8,"gene"]
## to.unwind<-c("NPM1")
## to.unwind<-c("FLT3")
## to.unwind<-c("Clinical")
## to.unwind<-c(clusters.wanted[!(clusters.wanted %in% c("Ubin.proteo","lipid_raft","caveolae","Checkpoint_extendedx1","Checkpoint_extendedx2"))])
## to.unwind
## to.unwind.name<-to.unwind


snpinfo.ex<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,]
loci<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,"Name"] # this is IDH1 not IDH1 in cluster # are the snp.names
the.genes<-unique(snpinfo.ex[,"gene"])
the.genes<-the.genes[!(the.genes %in% clusters.wanted)]

the.genes #245 ### if used a cluster name need to do back up to (**)

############repest to clean out cluster names 
to.unwind<-the.genes
snpinfo.ex<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,]
loci<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,"Name"] # this is IDH1 not IDH1 in cluster # are the snp.names
the.genes<-unique(snpinfo.ex[,"gene"])
the.genes<-the.genes[!(the.genes %in% clusters.wanted)]

the.genes

meta.results.skatO[1:5,]
meta.results.burden[1:5,]
the.genes.burden<-meta.results.burden[meta.results.burden[,"gene"] %in% the.genes,]

the.genes.burden
write.table(the.genes.burden,file=paste(to.unwind.name,"conponents:","Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

the.genes.burden<-meta.results.skatO[meta.results.skatO[,"gene"] %in% the.genes,]
the.genes.burden
write.table(the.genes.burden,file=paste(paste(to.unwind.name,collapse="."),"conponents:","SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
#subset<-(rownames(genotypes) %in% loci) # indicated in genotypes




dim(genotypes)
genotypes[1:5,1:5]
genotypes.ex<-genotypes[,loci]


dim(genotypes.ex)
genotypes.ex[is.na(genotypes.ex)]<-0
dim(genotypes.ex)


snpinfo.ex<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,]
dim(snpinfo.ex)
dim(genotypes.ex)
dim(pheno)
snpinfo.ex[1:5,]
# summary.geno.extra[loci,]
#high.missing[loci,]
sum(are.in.repeats[loci])

    
# qual[loci,]
## snpinfo[1:5,]
## qual[1:5,c("FILTER_PASS", "FILTER_100" )]
save(list=c("case.control","snpinfo.ori","formula","clusters","pheno.types","ipheno","clusters.wanted","genotypes","p","meta.results.skat","meta.results.skatO","meta.results.burden","pheno","target.pheno.col","snpinfo","fil.genotypes","pass","high.missing.table","a.indel","help","key","summary.geno.extra","full.qual","bad.effect","maf.filter","in.common.hit.gene","on.x.y","unannotated.hits","not.flat.genotype","are.repeats","are.in.repeats","ok.missing","hw.controls.ok.filt","no.genotypes","rare.in.all","are.in.repeats.back","are.in.repeats.forward"),file=paste(paste(project.files[ichr],".",pheno.types[ipheno],".",snap.file,".small_final.RData",sep="")) )


case.control<- colnames(pheno)[ipheno]
target.pheno.col<-pheno.types[ipheno]
 if(target.pheno.col %in% case.control){
 cohort.seq.ex <- skatCohort(Z=genotypes.ex,formula, SNPInfo = snpinfo.ex, data=pheno,aggregateBy="Name",family=binomial(),verbose=FALSE)

}else{
cohort.seq.ex <- skatCohort(Z=genotypes.ex,formula, SNPInfo = snpinfo.ex, data=pheno,aggregateBy="Name",family=gaussian(),verbose=FALSE) ## genes and clusters
}


meta.results.burden.ex<-burdenMeta(cohort.seq.ex,wts=1,mafRange = c(0,1),SNPInfo = snpinfo.ex,aggregateBy = "Name")
meta.results.burden.ex
pheno[1:5,]



cohort.seq.test <- skatCohort(genotypes.ex, formula, SNPInfo = snpinfo.ex, data=pheno,aggregateBy = "cluster",verbose=FALSE)

## meta.results.burden.test<-burdenMeta(cohort.seq.test,wts=1,mafRange = c(0,1),SNPInfo = snpinfo.ex,aggregateBy = "cluster")
## #meta.results.burden.test

## ## meta.results.skat.ex<-skatMeta(cohort.seq,SNPInfo = snpinfo)
## meta.results.skatO.test<-skatOMeta(cohort.seq.test,burden.wts =1,SNPInfo = snpinfo.ex,aggregateBy="cluster")
## #meta.results.skatO.test

muts.in.cases<-apply(genotypes.ex[pheno[,"cancer"],],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
muts.in.controls<-apply(genotypes.ex[pheno[,"Control"],],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})

## muts.in.cases<-apply(genotypes.ex[pheno[,"low"],],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
## muts.in.controls<-apply(genotypes.ex[pheno[,"high"],],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})

## muts.in.cases<-apply(genotypes.ex[pheno[,"SampleProject"]==1,],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
## muts.in.controls<-apply(genotypes.ex[pheno[,"SampleProject"]==0,],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})

figure<- match(loci,key)


########################################################
check<-16

quality.cases<-rep("",times=length(loci))
quality.controls<-rep("",times=length(loci))
a.indel.sub<-a.indel[figure,]

for(check in 1:length(loci)){
print(check)
#check<-"chr11:130066457:130066457:-:A:indel"
# posn<-grep(loci[check],key)
posn<-check

if(muts.in.cases[check]!=""){
#the.gt<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GT",sep=".")
  the.gq<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GQ",sep=".")
#the.dp<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"DP",sep=".")
#the.ad<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"AD",sep=".") # 4 AOGC


quality.cases[check]<-paste(a.indel.sub[posn,the.gq],collapse=",")
#quality.cases[check]<-paste(a.indel.sub[posn,the.ad],collapse=";") # 4 AOGC

## a.indel[posn,the.gq]
## a.indel[posn,the.gt]
## a.indel[posn,the.dp]
}

if(muts.in.controls[check]!=""){
#the.gt<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"GT",sep=".")
#the.gq<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"GQ",sep=".")
#the.gq<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"DP",sep=".")
the.ad<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"AD",sep=".")

#quality.controls[check]<-paste(a.indel.sub[posn,the.gq],collapse=",")
quality.controls[check]<-paste(a.indel.sub[posn,the.ad],collapse=";")

## a.indel[posn,the.gq]
## a.indel[posn,the.gt]
## a.indel[posn,the.dp]
}

} # end check
##########################################################################

                             
#figure
length(figure)
dim(meta.results.burden.ex)
length(muts.in.cases)
length(muts.in.controls)
#pass[figure]
#help[figure,]
dim(annotations)
dim(help)
dim(summary.geno.extra)
length(figure)

                             
sum(meta.results.burden.ex[,"gene"]!=loci)
## colnames(a.indel)[1:50]

## key[grep("chr17",key)[1:100]]
## grep("chr17:41197708",key)
## key[grep("10088407",key)]
#out<-cbind(meta.results.burden.ex,a.indel[figure,c(1:6,16,28,7,30,34,37:42,43)],summary.geno.extra[figure,],high.missing[figure,],help[figure,])
## out<-cbind(meta.results.burden.ex,a.indel[figure,c(1:6,16,28,7,30,34,37:42,43,14,32,33)],summary.geno.extra[figure,c("GENO.AML","GENO.Control","GENO.AML.filt","GENO.Control.filt")],high.missing[figure,])
summary.geno.extra[figure,]
annotations[figure,]
help[figure,]

dim(meta.results.burden.ex)
#out<-cbind(meta.results.burden.ex,a.indel[figure,c(1:6,16,43,28,7,30,34,37:42)],summary.geno.extra[figure,c("GENO.AML","GENO.Control","GENO.AML.filt","GENO.Control.filt")],help[figure,],muts.in.cases,muts.in.controls)
summary.geno.extra.cols<-colnames(summary.geno.extra)[grepl("^GENO",colnames(summary.geno.extra))]

                             
out<-cbind(meta.results.burden.ex,annotations[figure,],summary.geno.extra[figure,summary.geno.extra.cols],help[figure,],muts.in.cases,quality.cases,muts.in.controls,quality.controls)

#out<-cbind(meta.results.burden.ex,annotations[figure,],muts.in.cases,muts.in.controls)
dim(out)


## table(out[,"refGene::location"])
## table(out[,"Consequence.Embl"])

paste(paste(to.unwind,collapse="."))
paste(to.unwind.name,collapse=".")
  paste(paste(to.unwind.name,collapse="."),"GENOTYPE.conponents:","SkatO","clusters",snap.file,"txt",sep=".")
write.table(out,file=paste(paste(to.unwind.name,collapse="."),"GENOTYPE.conponents:","SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

getwd()

#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################






dim(genotypes)
rownames(genotypes)

####### update affextion status

icc<-2
 if(target.pheno %in% case.control){
   for(icc in 1:length(case.control.classes)){
     recode<-  pheno[,target.pheno] %in% names(case.control.classes)[icc]
     pheno[recode,target.pheno]<-as.numeric(case.control.classes[icc])
   }}
  
pheno[,target.pheno]<-as.numeric(pheno[,target.pheno])
formula

## pheno[1:5,]
## pheno[,target.pheno]



snpinfo.raw[1:5,]
tail(snpinfo.raw)



dups<-duplicated(snpinfo[,"Name"])
sum(dups)
snpinfo[dups,][1:10,]
snpinfo<-snpinfo[!dups,]
tail(snpinfo)


ic<-1
for(ic in 1:length(clusters.wanted)){
   cluster.genes<-clusters[,clusters.wanted[ic]]
   cluster.genes<-cluster.genes[cluster.genes!="" | is.na(cluster.genes)]
   print(paste(clusters.wanted[ic],paste(cluster.genes,collapse=","),sep=": "))
   extra<-snpinfo.raw[snpinfo.raw[,"gene"] %in%  cluster.genes,]
   print(dim(extra))
   extra[,"cluster"]<-clusters.wanted[ic]
   snpinfo<-rbind(snpinfo, extra)

 }
snpinfo.ori<-snpinfo


tail(snpinfo)

meta.results.skatO<-meta.results.skatO.all
meta.results.burden<-meta.results.burden.all

## meta.results.skatO.all<-meta.results.skatO
## meta.results.burden.all<-meta.results.burden


meta.results.skatO<-meta.results.skatO[ !(meta.results.skatO[,"gene"] %in% clusters.wanted ) ,]
meta.results.burden<-meta.results.burden[ !(meta.results.burden[,"gene"] %in% clusters.wanted) ,]


the.order<-     order(meta.results.burden[,"p"])
sum(is.na(meta.results.burden[,"p"])) ## bad p-values shoudl not happen
meta.results.burden<-  meta.results.burden[the.order,]

meta.results.burden[1:50,]
meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]

the.order<-     order(meta.results.skat[,"p"])
meta.results.skat<-  meta.results.skat[the.order,]
meta.results.skat[1:50,]

the.order<-     order(meta.results.skatO[,"p"])
sum(is.na(meta.results.skatO[,"p"])) ## bad p-values shoudl not happen
meta.results.skatO<-  meta.results.skatO[the.order,]
meta.results.skatO[1:50,]

clusters.wanted<-c("crappy_clusteres")
z<-qchisq(meta.results.burden[ !(meta.results.burden[,"gene"] %in% clusters.wanted) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) ## if have no chisq valuse
z0<-qchisq(meta.results.skatO[ !(meta.results.skatO[,"gene"] %in% clusters.wanted ) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) 
z[1:5]

 par(mfrow=c(1,1),font=2,font.lab=2,font.axis=2,mgp=c(3.5,1,0),mar=c(5,5,4,2)+0.1)


meta.results.burden[ !(meta.results.burden[,"gene"] %in% clusters.wanted) ,][1:10,]


p.vals<-meta.results.burden[ !(meta.results.burden[,"gene"] %in% clusters.wanted) ,"p"]
 z0=qnorm(meta.results.skatO[,"p"]/2)
median(z,na.rm=TRUE)/0.456  #1.071491
median(z0,na.rm=TRUE)/0.456  #1.071491
mean(as.numeric(indata[,"CHISQ"]))  #1.381459

z0=qnorm(meta.results.burden[,"p"]/2)
lambda = round(median(z0^2)/.454,3)
lambda


range(as.numeric(indata[,"CHISQ"]))

source("http://bioconductor.org/biocLite.R") 
biocLite("GWASTools")
setRep
qq.Plot(pvals)



## Reads data
## S <- read.table(input,header=F)
## if (stat_type == "Z")
##    z=S[,1]
## if (stat_type == "CHISQ")
##    z=sqrt(S[,1])
## if (stat_type == "PVAL")
##    z0=qnorm(meta.results.skatO[,"p"]/2)
## ## calculates lambda
## lambda = round(median(z0^2)/.454,3)
## lambda
range(z)
z<-z0 # use SKAT0

the.plot<-my.qq.plot(z,dist="chisq",df=1,ncp=0,col="blue",ylab="Observed chi-squared value",xlab="Expected chi-squared value",main="",cex=1,xlim=c(0,22),ylim=c(0,76),cex.lab=2.0,cex.axis=2.0,font.lab=2,font.axis=2,lwd=2,line="robust",plot.it=TRUE) # function defined below

qq<-  qq.data(z,plot.it=FALSE)       ## qq plot used same method as in car library
points(qq$x,qq$y,col="magenta",pch=21)

symbols<-meta.results.burden[,"gene"]

symbols<-meta.results.skatO[,"gene"]
#####annotate curve
selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="red",cex=1,offset=1,atpen='TRUE') ##plate row col symbol
selected.data<-identify(qq$x,qq$y,labels=labels[qq$ord],col="red",cex=1,atpen='TRUE') ## sybmol
selected.data<-identify(qq$x,qq$y,labels=as.character(round(data.in[qq$ord],2)),col="forestgreen",cex=1.25,atpen='TRUE') # observed score
#####


leg.txt<-c(expression(paste("All SNPs ",bold(lambda)==1.07)),expression(paste("remove MHC region 6p (25.5-33.6)Mb  ",lambda==1.06)),"95% confidence intervals")

legend(8.0,6,leg.txt,col=c("magenta","blue","red"),lty=c(-1,-1,2),pch=c(1,1,-1),cex=1.25)

label<-"Burden_coding_0.01_labels"
label<-"Burden_coding_0.01"

savePlot(paste(label,"tiff",sep="."),type="tiff")
savePlot(paste(label,"png",sep="."),type="png")
save.image("qq.AS.paper.final.RData")



data <- data.van
 Mobspval <- sort(data$P)
 Mobspval <- Mobspval[!Mobspval==0]
 o <- -(log10(Mobspval))
 e <- -log10( 1:length(o)/length(o) )

 #Mobsmax <- 3 #trunc(max(Mlogobspval))+1
 #Mexpmax <- trunc(max(Mlogexppval))+1
 #if (is.infinite(Mobsmax)) {Mobsmax <- 3} else {Mobsmax <- Mobsmax}
 #plot(c(0,Mexpmax), c(0,Mexpmax), col="gray", lwd=1, type="l", xlab="Expected -log10 P value", ylab="Observed -log10 P value", xlim=c(0,Mexpmax), ylim=c(0,Mobsmax), las=1, xaxs="i", yaxs="i", bty="l",main=hlabel)
 #plot(c(0,Mexpmax), c(0,Mexpmax), col="gray", lwd=1, type="l", xlab="Expected -log10 P value", ylab="Observed -log10 P value", las=1, xaxs="i", yaxs="i", bty="l",main=hlabel)
 # plot(c(0,Mexpmax), c(0,Mexpmax), col="gray", lwd=1, type="l", xlab="Expected -log10 P value", ylab="Observed -log10 P value", las=1, xaxs="i", yaxs="i", bty="l",main=hlabel)
#points(Mlogexppval,Mlogobspval, pch=23, cex=.4, bg="black")
plot(e,o, pch=23, cex=.4, bg="black",main=hlabel, ylab="Observed -log10 P value",xlab="Expected -log10 P value")
 ########################################################
 # Exclude MHC
 #############
 A <- ((data[,"BP"]>20500000) & (data[,"BP"]<38500000) & (data[,"CHR"]==6))
 data.nomhc <- data[!A,]

 cMobspval <- sort(data.nomhc$P)
 cMlogobspval <- -(log10(cMobspval))
 cMexppval <- c(1:length(cMobspval))
 cMlogexppval <- -(log10( (cMexppval-0.5)/length(cMexppval)))
 points(cMlogexppval, cMlogobspval, pch = 23, cex=.4, bg="red",col="red")

 ########################################################
 # Control SNPs
 ############
 null.pos <- match(tls[,1],data[,"SNP"])
 v <- data[null.pos,]
 v <- v[!is.na(v[,1]),]

 #Second QQ
 tMobspval <- sort(v$P)
 tMlogobspval <- -(log10(tMobspval)) 
 tMexppval <- c(1:length(tMobspval))
 tMlogexppval <- -(log10( (tMexppval-0.5)/length(tMexppval))) 
 points(tMlogexppval, tMlogobspval, pch=23, cex=.4, col="blue",main=hlabel)

 legend.col <- c("black","red","blue")
 legend.words <- c("All SNPs","MHC SNPs excluded","Control SNPs only")
 x.pos <- 3
 y.pos <- 100
 #if(i==1) {x.pos <- 0.5;y.pos <- 10}
 legend(x.pos,y.pos,col=legend.col,legend=legend.words,pch=23)

}
  


#

pheno.use[1:5,]
pheno[1:5,]

tail(pheno.use)
tail(pheno)

posns<-match(pheno.use[,"SAMPLE"],pheno[,"SAMPLE"])
missing<-is.na(posns)
sum(missing) ## MUST BE ZERO
pheno.use[,"AffectionStatus"]<-pheno[posns,"AffectionStatus"]
tail(pheno.use)
tail(pheno)

table(pheno.use[,"AffectionStatus"])
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


meta.results.burden[1:50,]
meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]



test<-snpinfo[snpinfo[,"cluster"]==a.cluster,"cluster"]
meta.results.burden[meta.results.burden[,"gene"] %in% test,]

## formula <- "AffectionStatus ~ 1"
## formula<-formula(formula)
formula
to.unwind<-c("MND_omim")
to.unwind<-c("ALS_JI") # to.unwind<-meta.results.burden.test[1:44,"gene"] # to.unwind<-to.unwind[!(to.unwind %in% clusters.wanted)]
to.unwind<-c(the.genes)
clusters.wanted

to.unwind



snpinfo.ex<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,]
loci<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,"Name"] # this is IDH1 not IDH1 in cluster # are the snp.names
the.genes<-unique(snpinfo.ex[,"gene"])


meta.results.skatO[1:50,]
meta.results.burden[1:50,]
the.genes.burden<-meta.results.burden[meta.results.burden[,"gene"] %in% the.genes,]

the.genes.burden
#write.table(the.genes.burden,file=paste(to.unwind,"conponents:","Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

the.genes.burden<-meta.results.skatO[meta.results.skatO[,"gene"] %in% the.genes,]
the.genes.burden
#write.table(the.genes.burden,file=paste(paste(to.unwind,collapse="."),"conponents:","SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
#subset<-(rownames(genotypes) %in% loci) # indicated in genotypes



sum(subset)
length(subset)
snpinfo[1:5,]


dim(genotypes)
genotypes[1:5,1:5]
genotypes.ex<-genotypes[,loci]


dim(genotypes.ex)
genotypes.ex[is.na(genotypes.ex)]<-0
dim(genotypes.ex)


snpinfo.ex<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,]
dim(snpinfo.ex)
dim(genotypes.ex)
dim(pheno.use)
snpinfo.ex[1:5,]

    # snpinfo.ex[,"cluster"]<-snpinfo.ex[,"gene"]
# qual[loci,]
## snpinfo[1:5,]
## qual[1:5,c("FILTER_PASS", "FILTER_100" )]

 if(target.pheno.col %in% case.control){
 cohort.seq <- skatCohort(Z=genotypes,formula, SNPInfo = snpinfo, data=pheno,aggregateBy="cluster",family=binomial(),verbose=FALSE)

}else{
cohort.seq <- skatCohort(Z=genotypes,formula, SNPInfo = snpinfo, data=pheno,aggregateBy="cluster",family=gaussian(),verbose=FALSE) ## genes and clusters
}


cohort.seq.test <- skatCohort(genotypes.ex, formula, SNPInfo = snpinfo.ex, data=pheno.use,aggregateBy = "cluster",verbose=FALSE)

meta.results.burden.test<-burdenMeta(cohort.seq.test,wts=1,mafRange = c(0,1),SNPInfo = snpinfo.ex,aggregateBy = "cluster")
meta.results.burden.test
## cohort.seq.test$SOD1$maf
### meta.results.burden.test<-burdenMeta(cohort.seq.test,wts=10000,mafRange = c(0,1),SNPInfo = snpinfo.ex,aggregateBy = "cluster")
### meta.results.burden.test
### dbeta(maf,1,25)

## meta.results.skat.ex<-skatMeta(cohort.seq,SNPInfo = snpinfo)
meta.results.skatO.test<-skatOMeta(cohort.seq.test,burden.wts =1,SNPInfo = snpinfo.ex,aggregateBy="cluster")
meta.results.skatO.test




cohort.seq.ex <- skatCohort(genotypes.ex, formula, SNPInfo = snpinfo.ex, data=pheno.use,aggregateBy = "Name",verbose=FALSE)
## meta.results.skat.ex<-skatMeta(cohort.seq,SNPInfo = snpinfo)
meta.results.burden.ex<-burdenMeta(cohort.seq.ex,wts=1,mafRange = c(0,1),SNPInfo = snpinfo.ex,aggregateBy = "Name")
meta.results.burden.ex
pheno[1:5,]



pheno.use[1:5,]

rownames(genotypes.ex)[pheno.use[,"AffectionStatus"]==1]
rownames(genotypes.ex)[pheno.use[,"AffectionStatus"]==0]

muts.in.cases<-apply(genotypes.ex[pheno.use[,"AffectionStatus"]==1,],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
muts.in.controls<-apply(genotypes.ex[pheno.use[,"AffectionStatus"]==0,],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})

snpinfo.ex[1:10,]
muts.in.cases[1:5]
dim(genotypes.ex)
muts.in.cases.gene<-tapply(muts.in.cases,snpinfo.ex[,"gene"],function(x) {length(x[x!=""])})
muts.in.controls.gene<-tapply(muts.in.controls,snpinfo.ex[,"gene"],function(x) {length(x[x!=""])})
               


muts.in.cases.gene[13]
muts.in.controls.gene


muts.in.cases[13]


figure<- match(loci,key)
figure
#pass[figure]
#help[figure,]
dim(annotations)
dim(help)
dim(summary.geno.extra)
length(figure)
meta.results.burden.ex
sum(meta.results.burden.ex[,"gene"]!=loci)
## colnames(a.indel)[1:50]

## key[grep("chr17",key)[1:100]]
## grep("chr17:41197708",key)
## key[grep("10088407",key)]
#out<-cbind(meta.results.burden.ex,a.indel[figure,c(1:6,16,28,7,30,34,37:42,43)],summary.geno.extra[figure,],high.missing[figure,],help[figure,])
## out<-cbind(meta.results.burden.ex,a.indel[figure,c(1:6,16,28,7,30,34,37:42,43,14,32,33)],summary.geno.extra[figure,c("GENO.AML","GENO.Control","GENO.AML.filt","GENO.Control.filt")],high.missing[figure,])
## summary.geno.extra[figure,]
## annotations[figure,]
geno.col<-grep("^GENO",colnames(summary.geno.extra))
out<-cbind(meta.results.burden.ex,muts.in.cases,muts.in.controls,annotations[figure,],summary.geno.extra[figure,][,geno.col] )  # help[figure,])


## meta.results.burden.ex
## out
snap.file<-"prelim_0.001_coding"
out[13,]

if(length(to.unwind)>5){to.unwind<-paste("Top_Hits")}

write.table(out,file=paste(paste(to.unwind,collapse="."),"GENOTYPE.conponents","SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

meta.results.burden.test
the.order<-     order(meta.results.burden.test[,"p"])
meta.results.burden.test<-  meta.results.burden.test[the.order,]

genes.order<-meta.results.burden.test[,"gene"]
muts.in.cases.gene<-muts.in.cases.gene[genes.order]
muts.in.controls.gene<-muts.in.controls.gene[genes.order]

posns<-match(genes.order,meta.results.skatO.test[,"gene"])
meta.results.skatO.test<-meta.results.skatO.test[posns,]
SkatO.p<-meta.results.skatO.test[,"p"]

out.gene<-cbind(meta.results.burden.test,SkatO.p,muts.in.cases.gene,muts.in.controls.gene)
out.gene
write.table(out.gene,file=paste(paste(to.unwind,collapse="."),"conponents.","SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
getwd()

#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################

column.labels<-read.delim(vcf.files,header=F,nrows=1,skip=skip.lines,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")

Coverage at location
file<-"/media/UQCCG/Sequencing/CompleteGenomics/Coverage/Paul.txt"
file<-"/media/UQCCG/Sequencing/CompleteGenomics/Coverage/Russel.txt"


cov1<-read.delim(file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"",as.is=TRUE)
dim(cov1)
colnames(cov1)[1:50]
tail(cov1[,1:5])

wanted<-  cov1[,"Approved.Symbol"]=="FLT3" &  cov1[,"exonRank"]==14
sum(wanted)

cov1<-cov1[wanted,]

dim(cov1)
colnames(cov1)<-gsub("^X","",colnames(cov1))
cov1[,1:10]
#colnames(cov1)
wanted<-  grepl(".Mean$",colnames(cov1))
sum(wanted)
cov1<-cov1[,wanted]
dim(cov1)
dim(cov)
cov<-cbind(cov,cov1)
# cov<-cov1
cov
dim(cov)

coverage[1:5,]

the.samples<-gsub(".Mean","",colnames(cov))
posns<-match(the.samples,coverage[,"Sample"])
missing<-is.na(posns)
the.samples[missing]
dim(t(cov))
dim(coverage[posns,c("Sample","Capture.Method","median.mean.coverage")])
data<-cbind(t(cov),coverage[posns,c("Sample","Capture.Method","median.mean.coverage")])
data[1:5,]
data<-data[!missing,]

remove<-data[,"Sample"] %in% c(51:100)

data<-data[!remove,]

col=rep("blue",times=dim(data)[1])
col[data[,"Capture.Method"]=="TruD:TruX"]<-"red"

the.plot<-my.qq.plot(z,dist="chi
sq",df=1,ncp=0,col="blue",ylab="Observed chi-squared value",xlab="Expected chi-squared value",main="",cex=1,xlim=c(0,22),ylim=c(0,76),cex.lab=2.0,cex.axis=2.0,font.lab=2,font.axis=2,lwd=2,line="robust",plot.it=TRUE)


data[1:5,]
plot(data[,4],data[,1],col=col,xlab="Mean Exome Coverage",ylab="Mean IDT coverage",main="",cex=1,,cex.lab=2.0,cex.axis=2.0)

lm(data[,1]~data[,4])
abline(a=-14.6,b=2.172)

selected.data<-identify(data[,4],data[,1],labels=data[,"Sample"],col="red",cex=1,offset=1,atpen='TRUE')


leg.txt<-c("Nextera","TrueSeq")

legend(70,70,leg.txt,col=c("blue","red"),lty=c(-1,-1,2),pch=c(1,1,-1),cex=1.25)

label<-"exome_coverage"
savePlot(paste(label,"tiff",sep="."),type="tiff")
savePlot(paste(label,"png",sep="."),type="png")
getwd()

r<-seq(from=0.01, to=0.3, by=0.001)
r
#cbind(r,(5+5/r))
log10(500)
plot(r,(5+5/r),xlab="Allelic Ratio",col="blue",ylab="Estimated FLT3-ITD coverage",main="",cex=1,cex.lab=2.0,cex.axis=2.0,type="l")

label<-"ITD.coverage.needed"
savePlot(paste(label,"tiff",sep="."),type="tiff")
savePlot(paste(label,"png",sep="."),type="png")


#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################

my.qq.plot<-function (x, distribution = "chisq",df=1,ncp=0, ylab = deparse(substitute(x)),
    xlab = paste(distribution, "quantiles"), main = NULL, las = par("las"), 
    envelope = 0.95, labels = FALSE, col = palette()[2], lwd = 2, 
    pch = 1, cex = 1, line = c("quartiles", "robust", "none"),xlim=c(0,100),ylim=c(0,20),font.lab=2,font.axis=2,font.main=2,cex.lab=2.5,cex.axis=1.0,plot.it=TRUE, ...){
    result <- NULL
    line <- match.arg(line)
    good <- !is.na(x)
    ord <- order(x[good])
    ord.x <- x[good][ord]
    q.function <- eval(parse(text = paste("q", distribution, 
        sep = "")))
    d.function <- eval(parse(text = paste("d", distribution, 
        sep = "")))
    n <- length(ord.x)
    P <- ppoints(n)
    z <- q.function(P,df=df,ncp=ncp, ...)
    if(plot.it){
    plot(z, ord.x, xlab = xlab, ylab = ylab, main = main, las = las, 
        col = col, pch = pch,cex = cex,xlim=xlim,ylim=ylim,font.lab=font.lab,font.axis=font.axis,font.main=2,cex.lab=cex.lab,cex.axis=cex.axis)}
    if (line == "quartiles") {
        Q.x <- quantile(ord.x, c(0.25, 0.75))
        Q.z <- q.function(c(0.25, 0.75),df=df,ncp=ncp, ...)
        b <- (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
        a <- Q.x[1] - b * Q.z[1]
        if(plot.it){
        abline(a, b, col = "red", lwd = lwd)}
    }
    if (line == "robust") {
        if (!require("MASS")) 
            stop("MASS package not available")
        coef <- coefficients(rlm(ord.x ~ z))
        a <- coef[1]
        b <- coef[2]
          if(plot.it){
        abline(a, b,col="red")}
    }         ###################  Envelope function
    if (line != "none" & envelope != FALSE) {
        zz <- qnorm(1 - (1 - envelope)/2)
        SE <- (b/d.function(z,df=df,ncp=ncp, ...)) * sqrt(P * (1 - P)/n)
        fit.value <- a + b * z
        upper <- fit.value + zz * SE
        lower <- fit.value - zz * SE
          if(plot.it){
        lines(z, upper, lty = 2, lwd = lwd, col = "red")
        lines(z, lower, lty = 2, lwd = lwd, col = "red")}
    }       #####################
    if (labels[1] == TRUE & length(labels) == 1)
        labels <- seq(along = z)
    if (labels[1] != FALSE) {
        selected <- identify(z, ord.x, labels[good][ord])
        result <- seq(along = x)[good][ord][selected]
    }
    if (is.null(result)) 
        invisible(list(result=result,a=a,b=b,x=z,y = ord.x,ord=ord,upper=upper,lower=lower))
    else {sort(result)
           invisible(list(result=result,a=a,b=b,x=z,y = ord.x,ord=ord,upper=upper,lower=lower))}
}



      qq.data<- function (x, plot.it = TRUE, distribution = "chisq", df=1,ncp=0, xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)) , ...)
{
    good <- !is.na(x)
    ord <- order(x[good])
    ord.x <- x[good][ord]
    q.function <- eval(parse(text = paste("q", distribution, 
        sep = "")))
    n <- length(ord.x)
    P <- ppoints(n)
    z <- q.function(P,df=df,ncp=ncp, ...)

    if (plot.it)
        plot(z, ord.x, xlab = xlab, ylab = ylab, ...)
    invisible(list(x = z, y = ord.x, ord=ord))
} ##ord is the order if use identify


######################################### END SECTION

## qq<- qq.data(data.in,distribution="norm",the.mean=the.mean,the.sd=the.sd,plot.it=FALSE)

## my.qq.plot(data.in,distribution="norm",col="blue",xlab="Expected Score",ylab="Observed score",xlim=range(qq$x), ylim=range(data.in),main=paste("Screen:",the.screen,"with 95% confidence intervals for",":",the.score,sep=" "),the.mean=the.mean,the.sd=the.sd,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex=1.5)
########################## USE FUNTIONS BELOW IF YOU REFERNCE FUNCYION IS A NORMAL DISTRIBUTION NOT A CHISQ
my.qq.plot.mean<-function (x, distribution = "norm", ylab = deparse(substitute(x)),
    xlab = paste(distribution, "quantiles"), main = NULL, las = par("las"), 
    envelope = 0.95, labels = FALSE, col = palette()[2], lwd = 2, the.mean=0,the.sd=1,cex.lab=2,
    pch = 1, cex = 1, line = c("quartiles", "robust", "none"),xlim=c(0,100),ylim=c(0,20),font.lab=2,font.axis=2,font.main=2,cex.axis=1,cex.main=1,
    ...)
{
    result <- NULL
    line <- match.arg(line)
    good <- !is.na(x)
    ord <- order(x[good])
    ord.x <- x[good][ord]
    q.function <- eval(parse(text = paste("q", distribution, 
        sep = "")))
    d.function <- eval(parse(text = paste("d", distribution, 
        sep = "")))
    n <- length(ord.x)
    P <- ppoints(n)
    z <- q.function(P, mean=the.mean, sd=the.sd, ...)
    plot(z, ord.x, xlab = xlab, ylab = ylab, main = main, las = las, 
        col = col, pch = pch,cex = cex,xlim=xlim,ylim=ylim,cex.lab=cex.lab,font.lab=font.lab,font.axis=font.axis,font.main=font.main,cex.main=cex.main,cex.axis=cex.axis)
    if (line == "quartiles") {
        Q.x <- quantile(ord.x, c(0.25, 0.75))
        Q.z <- q.function(c(0.25, 0.75), mean=the.mean, sd=the.sd, ...)
        b <- (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
        a <- Q.x[1] - b * Q.z[1]
        abline(a, b, col = "red", lwd = lwd)
    }
    if (line == "robust") {
        if (!require("MASS")) 
            stop("MASS package not available")
        coef <- coefficients(rlm(ord.x ~ z))
        a <- coef[1]
        b <- coef[2]
        abline(a, b)
    }         ###################  Envelope function
    if (line != "none" & envelope != FALSE) {
        zz <- qnorm(1 - (1 - envelope)/2)
        SE <- (b/d.function(z, mean=the.mean, sd=the.sd, ...)) * sqrt(P * (1 - P)/n)
        fit.value <- a + b * z
        upper <- fit.value + zz * SE
        lower <- fit.value - zz * SE
        lines(z, upper, lty = 2, lwd = lwd/2, col = "red")
        lines(z, lower, lty = 2, lwd = lwd/2, col = "red")
    }       #####################
    if (labels[1] == TRUE & length(labels) == 1)
        labels <- seq(along = z)
    if (labels[1] != FALSE) {
        selected <- identify(z, ord.x, labels[good][ord])
        result <- seq(along = x)[good][ord][selected]
    }
    if (is.null(result)) 
        invisible(result)
    else sort(result)
}

