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
library("BSgenome.Hsapiens.UCSC.hg19")
the.chroms<-seqlengths(Hsapiens)


############################################ Nimblegen and illuma capture loci ####################################################




the.sample.sheet

sample.sheet.full<-read.delim(the.sample.sheet,header=T,sep=",",fill=TRUE,stringsAsFactors=FALSE)
sample.sheet.full[1:5,]
colnames(sample.sheet.full)

sample.sheet.full[sample.sheet.full[,"AffectionStatus"]==2,"AffectionStatus"]<- "NMD"
sample.sheet.full[sample.sheet.full[,"AffectionStatus"]==1,"AffectionStatus"]<- "Control"
sample.sheet.full[sample.sheet.full[,"AffectionStatus"]==-9,"AffectionStatus"]<- "AOGC"

sample.sheet.full[,"SampleProject"]<-sample.sheet.full[,"AffectionStatus"] ## SampleProject used to subset groups
ex.controls<-grepl("^SH",sample.sheet.full[,"ParticipantCode"])
sample.sheet.full[ex.controls,"SampleProject"]<-"ex.Control"

##### fix 0 and 9 for missing to NA

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################


clusters<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-08-21_MND_AOGC_Tulane/Analysis/Final_NMD_clusters.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
clusters


clusters.wanted<-colnames(clusters)
ic<-1


gene.aliases<-read.delim("/media/UQCCG/Software/annovar/humandb/Gene_symbol_aliases.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
gene.aliases[1:5,]
gene.aliases<-unwind(gene.aliases, "Aliases",delimit=", ")
gene.aliases[1:5,]
###############  all.genes

all.genes<-gene.aliases[,"Approved.Symbol"] ### all genes is otherwise the names in the indel file
all.genes<-table(all.genes)
ic<-5
recode<-{}
for(ic in 1:length(clusters.wanted)){
   cluster.genes<-clusters[,clusters.wanted[ic]]
   cluster.genes<-cluster.genes[cluster.genes!="" | is.na(cluster.genes)]
   cluster.genes<-unique(cluster.genes)
#   all.genes[1:5]
   missing.name<-!(cluster.genes %in% names(all.genes))
   if(sum( missing.name)>0){
     posns<-match(cluster.genes[missing.name],gene.aliases[, "Aliases"])
     missing<-is.na(posns)
     recode<-cbind(cluster.genes[missing.name][!missing],gene.aliases[posns[!missing], "Approved.Symbol"])
     colnames(recode)<-c("old","new")
 ###### transfer to new gene lists
     posns<-match(clusters[,clusters.wanted[ic]],recode[,"old"])
     missing<-is.na(posns)
    clusters[!missing,clusters.wanted[ic]]<-recode[posns[!missing],"new"] ### redefine the clusters
  }
  }
 #########################################################    
     

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################


############## ASSUME HERE THAT SampleProjuect has classes Control and AML
pheno.types<-c("AffectionStatus")
 names(pheno.types)<-c("AffectionStatus")

case.control<-c("AffectionStatus")
# case.control.classes<-c(0,1,-9)
case.control.classes<-c(0,1,0)
#names(case.control.classes)<-c("Control","NMD","AOGC")
names(case.control.classes)<-c("Control","NMD","AOGC")
case.control.classes
# ib<-1
for(ib in 1:length(case.control)){
  if(!(case.control[ib] %in% colnames(sample.sheet.full))){next}
  sample.sheet.full[(  !(sample.sheet.full[,case.control[ib]] %in% names(case.control.classes))  |  is.na(sample.sheet.full[,case.control[ib]]) | sample.sheet.full[,case.control[ib]]==0 | sample.sheet.full[,case.control[ib]]==9)  ,case.control[ib]]<-NA
}

colnames(sample.sheet.full)[colnames(sample.sheet.full)=="ParticipantCode"]<-"SAMPLE"
#tapply(sample.sheet.full[,"SampleProject"],sample.sheet.full[,"SampleProject"],length)
sample.sheet.full[1:5,]

bad.samples.AOGC.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/exclusions/Samples_To_Exclude_Contaminated.csv"
bad.samples.AOGC<-read.delim(bad.samples.AOGC.file,header=F,sep="\t",fill=TRUE,stringsAsFactors=FALSE)


bad.samples<-read.delim(bad.samples.file,header=F,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

bad.samples<-c(as.character(bad.samples.AOGC[,1]),as.character(bad.samples[,1]))


sample.sheet.full<-sample.sheet.full[!(sample.sheet.full[,"SAMPLE"] %in% bad.samples),]
dim(sample.sheet.full)

control.samples<-{}


all.samples<-sample.sheet.full[,"SAMPLE"]





table(sample.sheet.full[,case.control])
table(sample.sheet.full[,"SampleProject"])
sum(is.na(sample.sheet.full[,case.control])) ## AOGC
length(all.samples)
## o.remove.all<-expand.labels.to.samples(remove.from.all.samples,all.samples)
## to.remove.samples<-unique(to.remove.all)
## remove.cols<-unique(c(remove.cols,to.remove.samples))

  
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
ipheno<-1
for(ipheno in 1:length(pheno.types)){
  print(paste("Doing phenotype:",pheno.types[ipheno]))


target.pheno<-names(pheno.types)[ipheno]
target.pheno.col<-pheno.types[ipheno]

length(sample.sheet.full[,target.pheno.col])
pheno<-sample.sheet.full[!is.na(sample.sheet.full[,target.pheno.col]) ,] ## pheno only contains SAMPLES that have a phenotype
print(dim(pheno))
print(paste("Number Samples:",dim(pheno)[1]))
covars<-c("PCA1","PCA2","PCA3","PCA4")
covars<-c("the.run")
covars<-c("1")
dim(pheno)
pheno[1:5,]


formula<-paste(target.pheno.col,"~",paste(covars,collapse="+"),sep="")
print(formula)
formula<-formula(formula)

#seq.type[1:57,]
posns<-match(pheno[,"SAMPLE"],seq.type[,"Sample"])
missing<-is.na(posns)
sum(missing)


capture<-rep("unknown",times=dim(pheno)[1])  # seq.type[posns,"Capture.Method"]
table(capture)
## the.runs<-seq.type[posns,"Run"]
table(capture) ### 
pheno<-cbind(pheno,capture,stringsAsFactors=FALSE)
pheno[1:5,]
pheno[pheno[,"SampleProject"]=="AOGC","capture"]<-"TruD:TruX"
pheno[pheno[,"SampleProject"]=="NMD","capture"]<-"TruD:NimX3"
pheno[pheno[,"SampleProject"]=="Control","capture"]<-"TruD:NimX3"
pheno[pheno[,"SampleProject"]=="ex.Control","capture"]<-"TruD:NimX"
table(pheno[,"capture"])


the.samples<-paste(pheno[,"SAMPLE"],"GT",sep=".")  ## samples same order as in pheno
print(paste("Number samples: ",length(the.samples),sep=""))



#### assume has format project.chr.fam.extension or chr.project.fam.extension
setwd(analysis.dir)
getwd()

files<-dir(analysis.dir)
the.extension<-paste(project.extension,"$",sep="")
files<-files[grepl(the.extension ,files)]
if(fam=="ALL" | fam=="All" | fam=="all"){
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
ichr<-1

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

for(ichr in 1:length(project.files)){ ### loop over chromosomes

print(project.files[ichr])
load(project.files[ichr])
meta.results.burden[1:50,]
 
meta.results.burden[meta.results.burden[,"gene"] %in% clusters[,1],]
meta.results.burden[meta.results.burden[,"gene"] %in% clusters[,2],]
annotations<-a.indel[,c(1:6,16,28,7,30,34,37:42,43,14,32,33)]

if(ichr==1){
meta.results.skatO.all<-meta.results.skatO
meta.results.burden.all<-meta.results.burden
#pheno.use.all<-pheno.use
snpinfo.all<-snpinfo
genotypes.all<-genotypes
pass.all<-pass
high.missing.all<-high.missing
annotations.all<-annotations
help.all<-help
key.all<-key
summary.geno.extra.all<-summary.geno.extra
}else{
meta.results.skatO.all<-rbind(meta.results.skatO.all,meta.results.skatO)
meta.results.burden.all<-rbind(meta.results.burden.all,meta.results.burden)
#pheno.use.all<-rbind(pheno.use.all,pheno.use
snpinfo.all<-rbind(snpinfo.all,snpinfo)
genotypes.all<-cbind(genotypes.all,genotypes)
pass.all<-cbind(pass.all,pass)
high.missing.all<-rbind(high.missing.all,high.missing)
annotations.all<-rbind(annotations.all,annotations)
help.all<-rbind(help.all,help)
key.all<-c(key.all,key)
summary.geno.extra.all<-rbind(summary.geno.extra.all,summary.geno.extra)
}



} # loop over projects

} # loop over fam

} # loop over Pheno

meta.results.skatO<-meta.results.skatO.all
meta.results.burden<-meta.results.burden.all
#pheno.use<-pheno.use.all
snpinfo<-snpinfo.all
genotypes<-genotypes.all
pass<-pass.all
high.missing<-high.missing.all
annotations<-annotations.all
help<-help.all
key<-key.all
summary.geno.extra<-summary.geno.extra.all

snpinfo.raw<-snpinfo




## save(list=c("meta.results.skat","meta.results.skatO","meta.results.burden","pheno.use","snpinfo","genotypes","pass","high.missing","annotations","help","key","summary.geno.extra"),file=paste(paste(project.files[ichr],".small.RData",sep="")) )
## print(paste("Done: ",project.files[ichr]))
## save.image(file=paste(project.files[ichr],".RData",sep=""))




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

## the.order<-     order(meta.results.skat[,"p"])
## meta.results.skat<-  meta.results.skat[the.order,]
## meta.results.skat[1:50,]

the.order<-     order(meta.results.skatO[,"p"])
sum(is.na(meta.results.skatO[,"p"])) ## bad p-values shoudl not happen
meta.results.skatO<-  meta.results.skatO[the.order,]
meta.results.skatO[1:50,]


z<-qchisq(meta.results.burden[ !(meta.results.burden[,"gene"] %in% clusters.wanted) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) ## if have no chisq valuse
z0<-qchisq(meta.results.skatO[ !(meta.results.skatO[,"gene"] %in% clusters.wanted ) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) 
z[1:5]

 par(mfrow=c(1,1),font=2,font.lab=2,font.axis=2,mgp=c(3.5,1,0),mar=c(5,5,4,2)+0.1)

median(z,na.rm=TRUE)/0.456  #1.071491
median(z0,na.rm=TRUE)/0.456  #1.071491
mean(as.numeric(indata[,"CHISQ"]))  #1.381459

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
QC_stat_BAM_Tue_Nov_18_2014.txt:
coverage<-read.delim("/media/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_BAM_Tue_Nov_18_2014.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")

Coverage at location
file<-"/media/UQCCG/Sequencing/CompleteGenomics/Coverage/Paul.txt"
file<-"/media/UQCCG/Sequencing/CompleteGenomics/Coverage/Russel.txt"
file<-"/media/UQCCG/Sequencing/CompleteGenomics/Coverage/AMAS.txt"

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

