
###### THIS IS the super annotion run USE ALL THE FILTER AND NOVEL DATABASES
UQCCG.data<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes" # path to project data 
project<-"2014-07-15_ChMND_UnifiedGenotyperForIBD" # this is the project directory
project.name<-"2014-07-15_ChMND_UnifiedGenotyperForIBD" # this is the prefix of the SNP and DINDEL files AND the project output NAME

annotate.dir<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-07-15_ChMND_UnifiedGenotyperForIBD/Annotate"
analysis.dir<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-07-15_ChMND_UnifiedGenotyperForIBD/Analysis"

affection.status<-"sample_sheet"
the.sample.sheet<-"Ch_MND_SampleSheet.csv"


skip.annovar.run<-FALSE # if FALSE it just redoes the files annovar summarization IF not already run
update.annovar.annotations<-TRUE # set TRUE will re-read VCF file and force annovar to always run
GATK.SB<-TRUE
genome.build<-"hg19"
dbSNP.build<-"131"
vcf.type="v4" # "annovar" "v4" "v3" "plink_assoc"

bam.extension<-".ReCal.sort.bam"
combined.extension<-""
snp.extension<-".recalibrated.filtered.vcf"  # "_snps.raw.vcf"
indel.extension<-".indelFiltered.vcf"    # "_DINDEL.raw.vcf"  
small.extension<-"_All_VARIANTS.raw.vcf"  # "_All_VARIANTS.raw.vcf"
variant.types<-c("snp") ##MUST have a extension type for each indel defined
names(variant.types)<-c("v4") ### define data type for when reading below
use.variants<-c("snp")
################3

###### THIS IS the super annotion run USE ALL THE FILTER AND NOVEL DATABASES
project.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis" # path to project data 
#project<-"2015-01-15_LeoPharma_Dec2014Freeze" # this is the project directory /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/SNPs/2015-02-11_LeoPharmFinalAnalysis.chrX.output.recalibrated.filtered.vcf
project.name<-"2015-02-11_LeoPharmFinalAnalysis" # this is the prefix of the SNP and DINDEL files AND the project output NAME



##affection.status<-"sample_sheet" This parameter is not used, not sure why here. Mhairi Feb 2015
the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/BAM/LeoPharma_Feb2015.Sample_Sheet_NEW.txt"

skip.annovar.run<-FALSE # if FALSE it just redoes the files annovar summarization IF not already run
update.annovar.annotations<-TRUE # set TRUE will re-read VCF file and force annovar to always run
GATK.SB<-TRUE
genome.build<-"hg19"
dbSNP.build<-"131"
vcf.type="v4" # "annovar" "v4" "v3" "plink_assoc"

bam.extension<-".ReCal.sort.bam"
combined.extension<-""
snp.extension<-".recalibrated.filtered.vcf"  # "_snps.raw.vcf"
indel.extension<-".indelFiltered.vcf"    # "_DINDEL.raw.vcf"  
small.extension<-"_All_VARIANTS.raw.vcf"  # "_All_VARIANTS.raw.vcf"
variant.types<-c("snp") ##MUST have a extension type for each indel defined
names(variant.types)<-c("v4") ### define data type for when reading below
use.variants<-c("snp")
################3


library(GenomicFeatures)
library(Rsamtools)
require(biomaRt)

#library(foreach)  ## could not get this to work well
library(doMC)
num.cores<-5
registerDoMC(cores=num.cores)
if(!exists("read.in.all.vcf")){read.in.all.vcf<-FALSE} # if FALSE will use Max.reads to read in file in chunks
if(!exists("max.reads")){max.reads<-10000 } # 0- get ALL IN ONE GO  :
## NA get all
## 1000 samples - 10000 
## 300 sample 70000 (7 cores) < 14 Gb
## 55 samples 20000 reads < 5Gb # 550,000 > 32Gb memory blown
## 40 samples 100,000 reafs 10Gb
###############10000 -1000 samples about 23Gb RAM' 12hrs
###############ALL   -58 samples about  23Gb 30 mins

if(!exists("project.name.alternative")){project.name.alternative<-project.name}
project.dir<-paste(UQCCG.data,project,sep="/")

snp.dir<-paste(project.dir,"SNPs",sep="/")
bam.dir<-paste(project.dir,"BAM",sep="/") # need the bam fir cause it contains the sample sheet
small.dir<-paste(project.dir,"SNPs",sep="/")
indel.dir<-paste(project.dir,"DINDELs",sep="/")
analysis.dir<-paste(project.dir,"Analysis",sep="/")


final.QC.directory<-"/mnt/UQCCG/Data/QC for all samples summary/Genetic_QC"

############################################################# set up the filter options
if(read.in.all.vcf){max.reads<-0} # max.reads of 0 will cause it to read in the whole file
missing.threshold<-0.80 # less.than   #  missing.threshold<-1.0
inbreeding.threshold<-0.75 # less than must be (1-inbreeding.threshold)*100 PERCENT NOT 0/1 OR MISSING # inbreeding.threshold<-1
hwe.threshold<-1e-4
use.variants<-c("snp") ### THIS MUST BE LOWER CASE SNP ELSE GENOTYPE FILTERING WILL FAIL

QC.dir<-paste(project.dir,"QC",sep="/")
if(!exists("annotate.dir")){annotate.dir<-paste(project.dir,"Annotate",sep="/")}
if(!exists("analysis.dir")){analysis.dir<-annotate.dir}

code.dir<-"/mnt/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
core.ann<-c("chr","start","end","REF","ALT","TYPE") 
setwd(code.dir)
source("annotate_SNPs_subroutines.r")
source("hwe.r")
###############################



################# Set up filtering


#if("SB" %in% info.types){GATK.SB<-TRUE; print("found SB tag")}else{GATK.SB<-FALSE; print("WARNING-No SB tag")}
## if(GATK.SB){  # old GATk has SB and the new has FS


############### just list ALL posibilities are then use indels to get relivant ones
##############$ coula also use start.data as well..


global.quality.labs<-c("QUAL","QD","FS","HRun","SB","FILTER","FILTER","TYPE") ### these become the "good.qual" filter
global.quality.names<-c("QUAL","QD","FS","HRun","SB","FILTER_PASS","FILTER_100","flat")
global.quality.cut<-c(50,0.5,60,5,1,"PASS","99.90to100.00","flat")
global.quality.type<-c("numeric","numeric","numeric","numeric","numeric","factor","factor","factor")
global.quality.dirn<-c("greater","greater","less","less","less","exact","ends_with","not_ends_with")

pass.filter<-c("FILTER_PASS","QUAL","QD","FS","HRun","SB","flat")  ## all these are required to be TRUE  (can have FILTER AND FILTER_100


names(global.quality.cut)<-global.quality.labs
names(global.quality.dirn)<-global.quality.labs
names(global.quality.type)<-global.quality.names

global.quality.cut
global.quality.dirn
global.quality.type

## quality.cut<-global.quality.cut
## quality.type<-global.quality.type
## quality.dirn<-global.quality.dirn

####################################################################################END: set up the filter options

variant.types<-variant.types[variant.types %in% use.variants]

i<- 1
for(i in 1:length(variant.types)){  ### CHOOSE BY project.name * extension
files<-dir( eval(as.name( paste(variant.types[i],"dir",sep=".") )) )
the.extension<-paste(eval(as.name( paste(variant.types[i],"extension",sep=".") )),"$",sep="")
the.extension<-paste(combined.extension,the.extension,sep="")
the.files<-files[grepl(the.extension ,files)]
the.files<-the.files[grepl(paste("^",project.name,sep=""),the.files)]

names(the.files)<-gsub(the.extension,"",the.files)
assign( paste(variant.types[i],"files",sep="."),value=the.files)
}

snp.files
snp.dir

########################### QC analysis only used SNPs
sample.files<-snp.files
names(sample.files)<-rep("snp",times=length(sample.files))

print(paste("Doing:  ",toString(sample.files),sep=""))

a.single.test<-FALSE ## need so will overwrite Aij
# isamp<-14 ; a.single.test<-TRUE
for (isamp in 1:length(sample.files)){

  
print(paste("Doing: ",sample.files[isamp],sep=""))

################## get initial read information
start.data<-prepare.for.Vcf.file.read(sample.files[isamp])
for(i in 1:length(start.data)){assign( names(start.data)[i],value=start.data[[i]])}
###########################################



setwd(eval(as.name(paste(names(sample.files)[isamp],"dir",sep="."))))
con <- file(sample.files[isamp], open="r")  # close(con)
num.lines<-1 # so does while llop at least once
reads<-max.reads  #1.3M lines in snp file 50000 goes out to 24Gb without QC cgeck 
counter<- -1
while (num.lines >0  ){
counter<-counter+1
print(counter)

if(counter==0){
indels<-try(scan(con,what=character(num.vars),skip=(reads*counter)+skip.lines,nlines=reads,sep="\t",fill=TRUE,na.strings="",quote="\""))
}else{
indels<-try(scan(con,what=character(num.vars),nlines=reads,sep="\t",fill=TRUE,na.strings="",quote="\""))
}


#$$$$$$$$$$$$$$



num.lines<-length(indels)/(num.vars)
print(num.lines)
if(num.lines==0){next}
dim(indels)<-c(num.vars,num.lines)
indels<-t(indels)
colnames(indels)<-column.labels

if(dim(indels)[1]<100){next} ### not worth the effort... not enough data


format.posn<-match("FORMAT",colnames(indels))
samples.processing<-length(column.labels)-format.posn # number of sample in the vcf file
samples.order.in.ALL<-column.labels[(format.posn+1):length(column.labels)] # samlple labels

the.samples<-samples.order.in.ALL

if(length(the.samples)>=40){
allele.min<-max(c(10,as.integer(samples.processing/10)))
}
if(length(the.samples)<=40 & length(the.samples)>=15  ){
allele.min<-2
}
if(length(the.samples)<15  ){
allele.min<-1
}
                                        # typically 10 
#allele.min/length(the.samples) # 10/40

####################################### FINISHED Read in data DO PROCESSIng below
###################################################################################################

indels<-  process.GATK.indels(indels,sample.files[isamp],vcf.type,format.types,info.types,info.class,num.cores) ## names(the.sample.file) must be the mutation type

rownames(indels)<-build.key(indels,core.ann,add.chr.label=FALSE)
######## Check not multi ALT alleles if there are then I need to flatten that line- not done here see build.annotation.files.r
####################################################################################################


all.sample.labels<-colnames(indels)
all.sample.labels<-all.sample.labels[grep("\\.GT$",all.sample.labels)]
all.sample.labels<-gsub(".GT$","",all.sample.labels)

############################################33




if(!grepl("^chr",indels[1,"chr"])){
key.indels<-build.key(indels,core.ann,add.chr.label=TRUE)
indels[,"chr"]<-paste("chr",indels[,"chr"],sep="") ## REQUIRE "chr in chromosome labels below
rownames(indels)<-key.indels
}else{key.indels<-build.key(indels,core.ann)      }
################################################################################################################


################################################################################################################
################################################################################################################
################################################################################################################

the.samples<-all.sample.labels
the.samples
print("Get filtered genoypes QC")
all.genotypes<-filtered.genotype.summary(indels,the.samples,prefix="",suffix=".ALL",20,0.2,0.80,0.05,0.95,10,5) # most often used

## while((dim(a.indel)[1] %% num.bits)< 2){num.bits<-num.bits+1} ### go don't get matrix issues
## num.bits
## fil.genotypes<-foreach(a.indel.bit=iter(a.indel,by='row',chunksize=as.integer(dim(a.indel)[1]/num.bits) ), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype(a.indel.bit,gsub(".GT$","",the.samples),prefix="",suffix="",20,0.02,0.98,0.20,0.80,7,2)


## all.genotypes<-filtered.genotype.summary(indels,the.samples,prefix="",suffix=".ALL",30,0.2,0.80,0.05,0.95,30,15) ## agggressive
## het.read.thresh<-20  # depth must exceed this to apply filer on genotypes
## het.low<-0.20
## het.high<-0.80
## het.low.ref<-0.05 #ref calles only
## het.high.alt<-0.95
## depth.het.low<-10
## depth.homo.low<-5






position.filter.full<-position.quality.filter(indels,global.quality.cut,global.quality.type,global.quality.dirn)
pass.filter.have<-pass.filter[pass.filter %in% colnames(position.filter.full)]

## position.filter.full[1:5,]

## names(all.genotypes)
dim(all.genotypes$"genotypes")
## all.genotypes$"genotypes"[1:5,1:20] ## sum(is.na(all.genotypes$"genotypes"[,45]))
## all.genotypes$"summary.geno.group"[45:50,]
## all.genotypes$"summary.depths.group"[1:2,]
#http://www.nature.com/nrg/journal/v11/n11/extref/nrg2865-s1.pdf

## test<-16
## tapply(all.genotypes$"genotypes"[,test],all.genotypes$"genotypes"[,test],length)

print("start QC")
all.genotypes$"genotypes"[all.genotypes$"genotypes"=="NA"]<-NA
all.genotypes$"genotypes"[all.genotypes$"genotypes"=="0/0"]<-0
all.genotypes$"genotypes"[all.genotypes$"genotypes"=="0/1"]<-1
all.genotypes$"genotypes"[all.genotypes$"genotypes"=="1/1"]<-2

p<-as.numeric(all.genotypes$"summary.geno.group"[,"MAF.ALL"])
p[is.na(p)]<-0

hw.p.ok<-getHWE(all.genotypes$"summary.geno.group"[,"GENO.ALL"])

hw.p.ok[45:50]

hw.p.ok<-hw.p.ok >  hwe.threshold
sum(hw.p.ok)
#test[1:500]

## Hg18 pseudoautosomal regions (PARs)
## chrY:1-2709520 and chrY:57443438-57772954 
## chrX:1-2709520 and chrX:154584238-154913754


## Hg19 pseudoautosomal regions (PARs)
## chrY:10001-2649520 and chrY:59034050-59363566 
## chrX:60001-2699520 and chrX:154931044-155260560

#mouse PAR: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC311143/figure/F6/  Fxy to tel:  166350000  End: 166,404,117
if(genome.build=="mm9"){
  parX<-(indels[,"chr"]=="chrX" & (indels[,"start"]> 166350000 & indels[,"start"]< 166650000) )
  print("Pseudo-Autosomal region on mouse mm9 coordinates")
}else{
parX<-(indels[,"chr"]=="chrX" & (indels[,"start"]> 60001 & indels[,"start"]< 2699520) ) | (indels[,"chr"]=="chrX" & (indels[,"start"]> 154931044 & indels[,"start"]< 155260560) )
print("Pseudo-Autosomal region on human hg19 coordinates")
}

on.X<-indels[,"chr"]=="chrX"


## missing.threshold<-0.20 # less.than 
## inbreeding.threshold<-0.95 # less than must be (1-inbreeding.threshold)*100 PERCENT NOT 0/1 OR MISSIN
####END identify pseudo-autosomeal regions and x chromosome

start.duplicated<-duplicated(indels[,"start"]) # don't need cause used TYPE==snp==exact

missing.rate<-(as.numeric(all.genotypes$"summary.geno.group"[,"MISSING.Alleles.ALL"])/2)/length(the.samples)
missing.ok<-missing.rate < missing.threshold  ## genotype is call not NA

inbreeding.rate<-((as.numeric(all.genotypes$"summary.geno.group"[,"MISSING.Alleles.ALL"])/2) + as.numeric(all.genotypes$"summary.geno.group"[,"ALT_HETRO.ALL"]))/length(the.samples)
inbreeding.ok<- (inbreeding.rate < inbreeding.threshold) | on.X


 ## Hardy-Weinberg Disequilibrium
# Given a SNP's frequency, we can calculated the expected number of homozygote wildtypes, heterozygotes, and homozygote mutants:

## hw_test <- function(x)
## {
## 	nid <- length(x)

## 	# calculate the frequency:
## 	p <- sum(x, na.rm=T) / (2*nid)

## 	# What are the expected genotype proportons?
## 	expected <- c(p^2 * nid, 2*p*(1-p) * nid, (1-p)^2 * nid)

## 	# What are the observed genotype proportions?
## 	observed <- c(
## 		sum(x == 2, na.rm=T),
## 		sum(x == 1, na.rm=T),
## 		sum(x == 0, na.rm=T))

## 	# Perform goodness-of-fit chi square test
## 	test <- fisher.test(rbind(observed, expected))
## 	return(test$p.value)
## }
## ## > rbind(expected,observed)
## ##          [,1] [,2] [,3]
## ## expected    1    2    1
## ## observed    2    4    2
##########################################################

REF.length<-nchar(as.character(indels[,"REF"]))
ALT.length<-nchar(as.character(indels[,"ALT"]))

large.indel<-REF.length>1 | ALT.length>1

are.repeats<-identify.repeats(indels,di.run.max=2,homo.run.max=4)

length(large.indel) 
length(are.repeats)
sum(are.repeats)
rownames(indels)[are.repeats][1:20]

#################### in repeats looking  forward

chk.in.repeat<-large.indel & !are.repeats
if(sum(chk.in.repeat)>1){
are.sub.repeat<-indentify.IN.repeat(indels[chk.in.repeat,],looking="forward",bases.about=6,di.run.max=3,homo.run.max=5,genome="BSgenome.Hsapiens.UCSC.hg19")
remove.repeats<-key.indels[chk.in.repeat][are.sub.repeat]
are.in.repeats.forward<- key %in% remove.repeats
}else{
  are.in.repeats.forward<-rep(FALSE,times=length(REF.length))
}
  
#remove.repeats[1:20]
sum(are.in.repeats.forward)
## [1] 6988

###################### in repeats looking back

sum(chk.in.repeat)
chk.in.repeat<-large.indel & !are.repeats & !are.in.repeats.forward
if(sum(chk.in.repeat)>1){
are.sub.repeat<-indentify.IN.repeat(indels[chk.in.repeat,],looking="back",bases.about=6,di.run.max=3,homo.run.max=5,genome="BSgenome.Hsapiens.UCSC.hg19")
remove.repeats<-key[chk.in.repeat][are.sub.repeat]
are.in.repeats.back<- key %in% remove.repeats
}else{
  are.in.repeats.back<-rep(FALSE,times=length(REF.length))
}
  
sum(are.in.repeats.back)
## [1] 3224

are.in.repeats<- are.in.repeats.back | are.in.repeats.forward

length(are.in.repeats)
sum(are.in.repeats)



##########################################################################
##########################################################################
##########################################################################

is.unwound.geno<-grepl("snp:\\d+$",indels[,"TYPE"]) | grepl("indel:\\d+$",indels[,"TYPE"])










##########################################################


maf.fil<-(p >= allele.min/length(the.samples)) & (p <= (1-allele.min/length(the.samples))) # maf.fil<-(p >= 0.01 & p <= 0.1)
maf.fil[is.na(maf.fil)]<- FALSE
sum(maf.fil)
#pass.filter

all.true<-length(pass.filter.have) ### only use snps and higly filtered data in the QC: position.filter.full[1:5,pass.filter]
position.filter.QC<-apply(position.filter.full[,pass.filter.have],1,function(x) sum(x)==all.true) # sum(position.filter)
# position.filter.QC<-  position.filter.QC & maf.fil &  missing.ok & inbreeding.ok
position.filter.QC<-  position.filter.QC & maf.fil &  missing.ok & inbreeding.ok & !is.unwound.geno & !are.in.repeats & !are.repeats & hw.p.ok
## apply(position.filter.full,2,sum)
print(sum(position.filter.QC))



## COMMENTED OUT BY MHAITI - July 2014
#### read plink
# ### cfind out what's in common
# chip.dir <- "/mnt/UQCCG-Analysis/AOGC_exome_chip/working_genotypes/"
# setwd(chip.dir)
# chip.file <- "recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL_chr22"
# chip.indels<- read.plink(chip.file)
# #geno.bim <- read.table("/mnt/UQCCG-Analysis/AOGC_exome_chip/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.bim", quote="\"")
# #chip.indels<-read.plink("file.chr1")
# 
# print("start QC")
# ## Fix the dosages
# chip.indels[chip.indels=="NA"]<-NA
# chip.indels[chip.indels=="0/0"]<-0
# chip.indels[chip.indels=="0/1"]<-1
# chip.indels[chip.indels=="1/1"]<-2
# chip.indels[chip.indels=="GENO"]<-"snp"
# 
# #all.genotypes$"genotypes"[all.genotypes$"genotypes"=="0/0"]<-0
# #all.genotypes$"genotypes"[all.genotypes$"genotypes"=="0/1"]<-1
# #all.genotypes$"genotypes"[all.genotypes$"genotypes"=="1/1"]<-2
# 
# core.ann<-c("chr", "start", "end", "REF", "ALT","TYPE")
# chip.key<-build.key(chip.indels,core.ann) # core.ann<-c("chr", start, end, REF ALT)
# seq.key <- rownames(all.genotypes$"genotypes")
# posns <- seq.key %in% chip.key
# all.genotypes$"genotypes" <- all.genotypes$"genotypes"[posns,]


## match to get common

## truen 
##  & in.common

## Merge IN COMMON genoypes cbind? sse and Genotype -- assume all genotyes are good quaility but check frequency. > 0.05%

## position.filter.QC<-  position.filter.QC[in.common] AND THE OTHER THAT ARE INPUT to genetic.sample***

## all.genotypes$"genotypes" IS A MATRIX

p<-as.numeric(all.genotypes$"summary.geno.group"[,"MAF.ALL"])
p[is.na(p)]<-0

write.plink(indels[position.filter.QC,],paste(project.name,"plink",counter,sep="_"))
## input genoytpe in doage'
## p allel frequency
## positio.filetr 
## the.samples
## on.X

## in future wirte this data using write.plink '
## /media/UQCCG-Analysis/AOGC_exome_chip/meta_analysis_single_point/final.meta.analysis.single.point.significant.xlsx
## sum(missing.ok)
## sum(inbreeding.ok)
## sum(position.filter.QC)
## sum(position.filter.full[,"RPA"]
## position.filter.full[1:5,]
## test<-indels[,"FS"]> 60
## indels[test,][1:4,c("FS","SB")]
##     plot(indels[,"FS"],indels[,"SB"],ylim=c(-4,0.1))
##   sum(as.numeric(indels[,"SB"])>0)
    
if(dim(all.genotypes$"genotypes")[1]<200){num.bits<-1; print(paste("Only",dim(all.genotypes$"genotypes")[1],"genotypes remaining",sep=" ")) }else{num.bits<-num.cores}  # dim(all.genotypes$"genotypes")  all.genotypes$"genotypes"[1:5,]
if(sum(position.filter.QC)!=0){
  
print("do QC check")
Aij<-foreach(genotypes.bit=iter(as.matrix(c(1:(dim(all.genotypes$"genotypes")[1]))),by='row',chunksize=as.integer(dim(all.genotypes$"genotypes")[1]/num.bits)), .combine='+', .multicombine=FALSE,.inorder=FALSE) %dopar% genetic.sample.QC.accumilate(all.genotypes$"genotypes"[as.integer(genotypes.bit),],p[as.integer(genotypes.bit)],position.filter.QC[as.integer(genotypes.bit)],the.samples,on.X[as.integer(genotypes.bit)],parX[as.integer(genotypes.bit)])

}else{
  print("WARNING  0 position left after filtering skipping Aij calculation")
}

if((counter==0 & isamp==1) | a.single.test ){
Aij.total<-Aij
a.single.test<-FALSE
print("first assign")
plink.files.all<-
}else{
Aij.total<-Aij.total+Aij
print("Add Aij")
}

} ## loop over data chunks

close(con)
 } ## loop over isamp-sample files which are all snps

related<-sum.QC.matrix(Aij.total,the.samples) # slight differences cause not corrected (correction small even for 200 genotypes)
#####################################################################################


############# write OUTPUT

xx<-try(setwd( analysis.dir  ),silent=TRUE)
if(inherits(xx, "try-error")){system(paste("mkdir",analysis.dir,sep=" "));setwd( analysis.dir)}
  ## QC.files<-QC.files[grepl(".QC$",QC.files)]
print(project.name.alternative)
setwd(analysis.dir)
write.table(related,file=paste(project.name.alternative,"genetic_QC.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
#save(list=c("Aij.total","related","the.samples"),file="AOGC_all_related.RData")

###############################get the sample sheet and merge where possible

if(exists("the.sample.sheet")){

xx<-try(sample.sheet.full<-read.delim(paste(bam.dir,the.sample.sheet,sep="/"),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE))
if(inherits(xx, "try-error")){
sample.sheet.full<-read.delim(the.sample.sheet,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
}

if(!("ParticipantCode" %in% colnames(sample.sheet.full))){ # somtime it's called "ParticipantCode" or "Participant Code"
  if(sum(grepl("^Participant",colnames(sample.sheet.full) ))==1){colnames(sample.sheet.full)[grepl("^Participant",colnames(sample.sheet.full) )]<-"ParticipantCode"}
}

if(!("SampleProject" %in% colnames(sample.sheet.full))){ # somtime it's called "ParticipantCode" or "Participant Code"
  if(sum(grepl("^Sample Project",colnames(sample.sheet.full) ))==1){colnames(sample.sheet.full)[grepl("^Sample Project",colnames(sample.sheet.full) )]<-"SampleProject"}
}

if(!("SampleProject" %in% colnames(sample.sheet.full))){ # somtime it's called "ParticipantCode" or "Participant Code"
  if(sum(grepl("^Sample.Project",colnames(sample.sheet.full) ))==1){colnames(sample.sheet.full)[grepl("^Sample.Project",colnames(sample.sheet.full) )]<-"SampleProject"}
}
} #exists("the.sample.sheet"


#########################################################################333

sample.sheet.full[1:5,]
related[1:5,]
samples.not.in.samplesheet<- the.samples[!(the.samples %in% sample.sheet.full[,"ParticipantCode"])]
if(length(samples.not.in.samplesheet)>0){print(paste("WARNING Samples missing from SampleSheet: ",toString(samples.not.in.samplesheet),sep=""))}


posns<-match(related[,"sample_A"],sample.sheet.full[,"ParticipantCode"])
diagonal<-related[,"sample_A"]==related[,"sample_B"]
Ident.Samples<-diagonal
related.sheet<-cbind(related,Ident.Samples,sample.sheet.full[posns,])
related.sheet[1:5,]

write.table(related.sheet,file=paste(project.name.alternative,"genetic_sampleSheet_QC.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

   
setwd(final.QC.directory)
out.file.name<-paste("Vcf_Merge.",project.name.alternative,".",format(Sys.time(), "%a_%b_%d_%Y"),".genetic_QC.txt",sep="")
write.table(related.sheet,out.file.name,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

diagonal<-related.sheet[,"sample_A"]==related.sheet[,"sample_B"]
bad.sex<-(related.sheet[,"sex_Predicted"]!=related.sheet[,"Sex"] | related.sheet[,"sex_Predicted"]==9  ) & diagonal & !is.na(related.sheet[,"Sex"])

## related.sheet[bad.sex,]
## tapply(related.sheet[,"sex_Predicted"],related.sheet[,"sex_Predicted"],length)
## tapply(related.sheet[,"Sex"],related.sheet[,"Sex"],length)

## out.file.name<-paste("Vcf_Merge.",project.name.alternative,".DIAGONAL.",format(Sys.time(), "%a_%b_%d_%Y"),".genetic_QC.txt",sep="")
## write.table(related.sheet[diagonal | bad.sex,],out.file.name,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


are.related<-as.numeric(related[,"IBS"])>0.06 & is.finite(as.numeric(related[,"IBS"]))
out.file.name<-paste("Vcf_Merge.",project.name.alternative,".ARE_RELATED.",format(Sys.time(), "%a_%b_%d_%Y"),".genetic_QC.txt",sep="")
write.table(related.sheet[are.related | bad.sex,],out.file.name,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)





###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
