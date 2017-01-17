#library(SKAT) ## skat method
suppressMessages(library(GenomicFeatures))
suppressMessages(library(HardyWeinberg))
suppressMessages(library(Biostrings))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg19"))
#library(doMC)
suppressMessages(library(foreach))
suppressMessages(library(doSNOW))
suppressMessages(library(doParallel))
suppressMessages(library(skatMeta))  ## ridge regression
#library(SKAT) ## skat method
suppressMessages(library(GenomicFeatures))
suppressMessages(library(HardyWeinberg))
suppressMessages(library(Biostrings))

# options(width=250,max.print=5000)

code.dir<-"P:\\PhD\\Exome Capture Data and Mutations\\Raw Sequencing Data\\AML with Controls - FANC Set"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")
source("hwe.r")

analysis.dir<-"P:\\PhD\\Exome Capture Data and Mutations\\Raw Sequencing Data\\AML with Controls - FANC Set\\Analysis"
annotate.dir<-"P:\\PhD\\Exome Capture Data and Mutations\\Raw Sequencing Data\\AML with Controls - FANC Set\\Annotate"

project.extension<-".All-maf-filtered.txt"
project.name<-"2015-01-27_AML_with_AOGCControl" ## prefix for output file
fam<-c("TGCM-AML") #  ALL or  c() ""-one project (the prefix of the summary files to collect

the.sample.sheet<-"P:\\PhD\\Exome Capture Data and Mutations\\Raw Sequencing Data\\AML with Controls - FANC Set\\TGCM-AML-combine_SampleSheet.csv"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)
remove.cols<-c()

#regions.file<-"/media/scratch2/AOGC-NGS/GFOS/gefos.seq/METHODS/0613-skatmeta-gefos/static/Homo_sapiens.GRCh37.70.protein_coding.genespace_boundaries.5k.split100k.txt"
core.ann<-c("chr","start","end","REF","ALT","TYPE") # out put to annanlsys programs and need foe colun labels
dont.build.summary<-FALSE ##

GATK.SB<-TRUE
maf.threshold.filter.to.use<-c(0.05)
a.label<-"CoVarRun.noControl.AML.regions"
###########################################################################


###################################### load ola aogc######################3
## load("/media/Bioinform-D/Research/annovar/humandb/aogc.count.data.RData")
## #print(colnames(indels))
## print(colnames(geno.aogc))
## use.key<-build.key(geno.aogc,core.ann)
## insert.location<-70  ### this is where to add AOGC data INSERTS AFTER THE LOCATION
## ################################ add aogc

## ann<-readRDS("/media/scratch2/AOGC-NGS/Analysis/AOGC_sequnnce_LS/AOGC_sequence_10_LS_ANNOTATION.rds")

geneanno.DB<-c("refGene","knownGene","ensGene") # returns 2 extra columns
names(geneanno.DB)<-c("refGene","knownGene","ensGene")


##################################################### DEFINE A GENE LIST  #####################################################


##################################################### PILOT - GENE LIST #####################################################
## gene.list.file<-"/media/Bioinform-D/Research/exome Chip/AOGC_Gene_Hits_0.05_0.2-sent.txt"
## gene.list<-read.delim(gene.list.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
## gene.list[1:5,]
## colnames(gene.list)[colnames(gene.list)=="chr"]<-"CHR"
## colnames(gene.list)[colnames(gene.list)=="start"]<-"START"
## colnames(gene.list)[colnames(gene.list)=="end"]<-"END"
## colnames(gene.list)[colnames(gene.list)=="value"]<-"GENE_NAME"
## gene.list[1:5,]

######################################################################################################################################
######################################################################################################################################

##################################################### GEFOS - GENE LIST #####################################################

## gene.list.file<-"/media/scratch2/AOGC-NGS/GFOS/gefos.seq/METHODS/0613-skatmeta-gefos/static/Homo_sapiens.GRCh37.70.protein_coding.genespace_boundaries.5k.split100k.txt"
## gene.list<-read.delim(gene.list.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

## #target<-"CHR"
## target<-"GENE_NAME"
## #target<-"GENE_ID"
## sort(tapply(gene.list[,target],gene.list[,target],length),decreasing=TRUE)[1:150]
## order.by<-order(gene.list[,"CHR"],gene.list[,"START"],gene.list[,"END"])
## gene.list<-gene.list[order.by,]
## gene.list[1:5,]


######################################################################################################################################
######################################################################################################################################
############################################SET UP FNCTIONAL  FILTERS #######################################################

####################### MUTATION TYPE DEINITIONS 

possible.mutations<-c("frameshift substitution","nonframeshift substitution","downstream","frameshift deletion","frameshift insertion","intergenic","intronic","ncRNA_exonic","ncRNA_intronic","ncRNA_splicing","ncRNA_UTR3","ncRNA_UTR5","ncRNA_UTR5;ncRNA_UTR3","nonframeshift deletion","nonframeshift insertion","nonsynonymous SNV","splicing","stopgain SNV","stoploss SNV","synonymous SNV","unknown","upstream","upstream;downstream","UTR3","UTR5","UTR5;UTR3")

# select for these
interesting.coding.mutations<-c("frameshift substitution","nonframeshift substitution","nonframeshift deletion","nonframeshift insertion","frameshift deletion","frameshift insertion","nonsynonymous SNV","stopgain SNV","stoploss SNV","splicing")

interesting.mutations.use<-c("frameshift substitution","nonframeshift substitution","nonframeshift deletion","nonframeshift insertion","frameshift deletion","frameshift insertion","nonsynonymous SNV","stopgain SNV","stoploss SNV","splicing","ncRNA_exonic")

wanted.noncoding.subtypes<-c("miRNA","lincRNA") # filter by interesting to prefiler and vep.noncoding so dones get ncRNA intronic ::use gerp.score.threshold.low only these subtypes

interesting.to.prefilter<-c("UTR3","UTR5","UTR5;UTR3","snoRNA","snRNA","antisense","sense_intronic","ncRNA_exonic","ncRNA_splicing") #use gerp.score.threshold

extra.vep.annotations<-c("Uploaded_variation","Gene","Feature","Protein_position","Amino_acids")

#"not_assigned",
vep.types<-c("stop_gained","stop_lost","missense_variant","splice_acceptor_variant","splice_donor_variant","splice_region_variant","initiator_codon_variant","stop_retained_variant","incomplete_terminal_codon_variant","frameshift_variant","inframe_deletion","inframe_insertion","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_exon_variant","NC_stop_gained","NC_stop_lost","NC_splice_acceptor_variant","NC_splice_donor_variant","NC_splice_region_variant","NC_initiator_codon_variant","NC_stop_retained_variant","NC_non_coding_exon_variant","NC_incomplete_terminal_codon_variant","NC_3_prime_UTR_variant","mature_miRNA_variant","NC_5_prime_UTR_variant","TF_binding_site_variant","TFBS_ablation","TFBS_amplification","regulatory_region_variant","intron_variant","NC_intron_variant","synonymous_variant","coding_sequence_variant","NC_synonymous_variant","upstream_gene_variant","downstream_gene_variant","intergenic_variant","NC_intergenic_variant","NMD_transcript_variant","nc_transcript_variant","NC_nc_transcript_variant","feature_truncation","feature_elongation")

## vep.types<-c( "not_assigned","stop_gained","stop_lost","stop_lost,NMD_transcript_variant","stop_gained,splice_region_variant,NMD_transcript_variant","initiator_codon_variant,splice_region_variant","splice_region_variant,3_prime_UTR_variant","stop_gained,NMD_transcript_variant","missense_variant,splice_region_variant","missense_variant","splice_acceptor_variant","splice_acceptor_variant,nc_transcript_variant","splice_region_variant,3_prime_UTR_variant,NMD_transcript_variant","splice_donor_variant,nc_transcript_variant","splice_region_variant,intron_variant,NMD_transcript_variant","splice_donor_variant","splice_region_variant","splice_region_variant,5_prime_UTR_variant","splice_region_variant,synonymous_variant","splice_region_variant,intron_variant,nc_transcript_variant","splice_region_variant,non_coding_exon_variant,nc_transcript_variant","missense_variant,NMD_transcript_variant","splice_region_variant,intron_variant","NMD_transcript_variant","intron_variant,NMD_transcript_variant","mature_miRNA_variant","5_prime_UTR_variant","5_prime_UTR_variant,NMD_transcript_variant","non_coding_exon_variant,nc_transcript_variant","3_prime_UTR_variant,NMD_transcript_variant","non_coding_exon_variant","TF_binding_site_variant","intron_variant,nc_transcript_variant","synonymous_variant,NMD_transcript_variant","3_prime_UTR_variant","regulatory_region_variant","upstream_gene_variant","downstream_gene_variant","intergenic_variant","intron_variant","synonymous_variant")

#"not_assigned",
vep.coding<-c("stop_gained","stop_lost","missense_variant","splice_acceptor_variant","splice_donor_variant","splice_region_variant","initiator_codon_variant","stop_retained_variant","incomplete_terminal_codon_variant","frameshift_variant","inframe_deletion","inframe_insertion")

vep.noncoding<-c("5_prime_UTR_variant","3_prime_UTR_variant","non_coding_exon_variant","NC_stop_gained","NC_stop_lost","NC_splice_acceptor_variant","NC_splice_donor_variant","NC_splice_region_variant","NC_initiator_codon_variant","NC_stop_retained_variant","NC_non_coding_exon_variant","NC_incomplete_terminal_codon_variant","NC_3_prime_UTR_variant","mature_miRNA_variant","NC_5_prime_UTR_variant","TF_binding_site_variant","TFBS_ablation","TFBS_amplification","regulatory_region_variant")              

vep.unwanted<-c("intron_variant","NC_intron_variant","synonymous_variant","coding_sequence_variant","NC_synonymous_variant","upstream_gene_variant","downstream_gene_variant","intergenic_variant","NC_intergenic_variant","NMD_transcript_variant","nc_transcript_variant","NC_nc_transcript_variant","feature_truncation","feature_elongation")

missense.variant<-c("nonsynonymous SNV","missense_variant")

hwe.control.threshold<-1e-8
gerp.score.threshold.high<-2.5 # gerp score >= will be included
gerp.score.threshold.low<-2.0 # gerp score >= will be included
gerp.score.threshold.unknown<-0


#generic.filter.DB
maf.threshold.filter.to.use<-sort(maf.threshold.filter.to.use)

maf.threshold.filter<-maf.threshold.filter.to.use
interesting.mutations<-interesting.mutations.use


if(GATK.SB) {
  global.quality.labs<-c("QUAL","QD","HRun","SB","FILTER","FILTER","PolyPhen.scores","PolyPhen.scores","SIFT.scores","mut.taster::score","phylo::score","PolyPhen.desc","SIFT.desc","GERP::score","GERP::score","GERP::score","MAF.ALL","MAF.HIGH","MAF.LOW","TYPE") ### THESE ARE THE COLUMN LABELS IN THE DATA these become the "good.qual" filter 
  global.quality.names<-c("QUAL","QD","HRun","SB","FILTER_PASS","FILTER_100","PolyPhen.low","PolyPhen.high","SIFT.high","mut.taster.high","phylo.high","PolyPhen.bad","SIFT.bad","GERP.high","GERP.low","GERP.unknown","MAF.ALL","MAF.HIGH","MAF.LOW","flat") ### THESE ARE THE COLUMN LABELS IN the quality.filter TABLE
  #global.quality.cut<-c(50,0.5,5,1,"PASS","TruthSensitivityTranche99.90to100.00",0.1,0.4,0.4,0.4,0.4,"damaging","deleterious",2,2,0.25,0.25,0.25)
  global.quality.cut<-c(50,0.5,5,1,"PASS","TruthSensitivityTranche99.90to100.00",0.1,0.4,0.4,0.4,0.4,"damaging","deleterious",gerp.score.threshold.high,gerp.score.threshold.low,gerp.score.threshold.unknown,0.25,0.25,0.25,"flat")
  global.quality.type<-c("numeric","numeric","numeric","numeric","factor","factor","numeric","numeric","numeric","numeric","numeric","factor","factor","numeric","numeric","numeric","numeric","numeric","numeric","factor")
  global.quality.dirn<-c("greater","greater","less","less","exact","exact","greater","greater","greater","greater","greater","ends_with","exact","greater","greater","exact","greater","greater","greater","ends_with")
  
} else {
  global.quality.labs<-c("QUAL","QD","HRun","FS","FILTER","FILTER") ### these become the "good.qual" filter
  global.quality.names<-c("QUAL","QD","HRun","FS","FILTER_PASS","FILTER_100")
  global.quality.cut<-c(50,0.5,5,60,"PASS","TruthSensitivityTranche99.90to100.00")
  global.quality.type<-c("numeric","numeric","numeric","numeric","factor","factor")
  global.quality.dirn<-c("greater","greater","less","less","exact","exact")
}
names(global.quality.cut)<-global.quality.labs
names(global.quality.dirn)<-global.quality.labs
names(global.quality.type)<-global.quality.names
global.labs<-unique(global.quality.labs)

global.quality.cut
global.quality.dirn
global.quality.type

quality.cut<-global.quality.cut
quality.type<-global.quality.type
quality.dirn<-global.quality.dirn
######################################################################################################################

# coverage information
seq.type.file<-"P:\\PhD\\Exome Capture Data and Mutations\\Raw Sequencing Data\\AML with Controls - FANC Set\\QC_stat_SAMPLE_Tue_Oct_14_2014.txt"
seq.type<-read.delim(seq.type.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

seq.type[1:5,]
## seq.type<-seq.type[seq.type[,"Project"]=="AOGC-NGS",]
nim.samples<-seq.type[seq.type[,"Project"]=="AOGC-NGS" & seq.type[,"Capture.Method"]=="TruD:NimX","Sample"]
ill.samples<-seq.type[seq.type[,"Project"]=="AOGC-NGS" & seq.type[,"Capture.Method"]=="TruD:TruX","Sample"]

nim.samples<-paste(nim.samples,"GT",sep=".")
ill.samples<-paste(ill.samples,"GT",sep=".")


length(nim.samples)
length(ill.samples)

############################################ Nimblegen and illuma capture loci ####################################################

the.chroms<-seqlengths(Hsapiens)


## load("/media/UQCCG/Sequencing/Data/Genomes/hg19/Human_Exome_Targets_illumina_v2_hg19_targets.RData")
## ill.gr<-data.gr
## load("/media/UQCCG/Sequencing/Data/Genomes/hg19/Human_Exome_Targets_Nimble_v2_hg19_targets.RData")
## nim.gr<-data.gr

## genome(ill.gr)<-"hg19"
## genome(nim.gr)<-"hg19"
## names(ill.gr)
## the.chromo<- as.character(unique(seqnames(ill.gr)))
## sum(the.chromo!=as.character(unique(seqnames(nim.gr))))==0 ## must be true else chr in different order


## human.chromlens<-the.chroms[the.chromo]

## ## ill.gr<-ill.gr+200
## ## nim.gr<-nim.gr+200
## ## overlaps<-overlapsAny(ill.gr,nim.gr)
## ## possible.loci.ori<-ill.gr[overlaps]

## ill.gr<-reduce(ill.gr)
## nim.gr<-reduce(nim.gr)
##       cov.ill.gr<-coverage(ill.gr,weight=5,width=human.chromlens) ## Problems exist here is the ir.gerp are not unique get coverage*weight
##       cov.nim.gr<-coverage(nim.gr,weight=5,width=human.chromlens) ## Problems exist here is the ir.gerp are not unique get coverage*weight

##       cov.nim.gr<-cov.nim.gr[names(cov.ill.gr)] # make sure in same order
##       cov.all<-cov.nim.gr+cov.ill.gr

##    x<- slice(cov.all,lower=9,upper=11) ### over is 10 otherwise is 5
##  x[[5]][1:5]


## x
## ########check


## regionViews<-x
## the.counts<-{}
## order.chromos<-names(regionViews)

## # ik<-13
## for (ik in 1:length(order.chromos)){
##   chromo<-order.chromos[ik]
##  ## print( chromo)

##   a.chr<-chromo
##   a.start<-start(regionViews[[chromo]])
##   a.end<-end(regionViews[[chromo]])
##   a.width<-width(regionViews[[chromo]])


##   a.set<-cbind(a.chr,a.start,a.end,a.width)

##   the.counts<-rbind(the.counts,a.set)
## }
## colnames(the.counts) <- c("chr","start","end","length")
## ## }) # system.time
## the.counts[1:5,]

## write.table(the.counts,file="Common_target_loci_between_Nimblegen_v2.and.v3.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## save(list=c("the.counts"),file="Common_target_loci_between_Nimblegen_v2.and.v3.Rdata")
## getwd()

#load("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/Common_target_loci_between_Nimblegen_v2.and.v3.Rdata")
#load("/media/UQCCG/Sequencing/Data/Genomes/hg19/NexteraRapidCapture_Exome_TargetedRegions_hg19_targets.RData")

######################################################################################################################################

## gene.gr<-GRanges(seqnames =gene.list[,"CHR"],ranges = IRanges(start=as.numeric(gene.list[,"START"]),end=as.numeric(gene.list[,"END"])),strand="+")
## length(gene.gr)
## length(reduce(gene.gr))

## location.key<-build.key(gene.list,c("CHR","START","END"))
## unique.locations<-unique(location.key)
## length(unique.locations)
## dim(gene.list)

## posns<-grep("GRIA1.",gene.list[,target],fixed=TRUE)
## gene.list[posns,]
## vep.coding<-c("stop_gained","stop_lost","missense_variant","initiator_codon_variant","stop_retained_variant","incomplete_terminal_codon_variant","frameshift_variant","inframe_deletion","inframe_insertion")
## vep.other.coding<-c("synonymous_variant","coding_sequence_variant")


the.sample.sheet

sample.sheet.full<-read.delim(the.sample.sheet,header=T,sep=",",fill=TRUE,stringsAsFactors=FALSE)
sample.sheet.full[1:5,]
colnames(sample.sheet.full)
dim(sample.sheet.full)


##### fix 0 and 9 for missing to NA

## pheno.types<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_RAD","BMD_EFF_STD_LS","BMD_EFF_STD_FN","EVER_FX_50_EXCL_TRIVIAL")
## names(pheno.types)<-c("HIP","RAD","LS","FN","FX")


######### Check and fix the sample sheet

# coverage information
coverage<-seq.type # read.delim("/media/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_SAMPLE_Tue_Oct_14_2014.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
dim(coverage)
coverage[1:5,]
unique(coverage[,"Project"])
aml.projects<-c("AML-exome","RSGB_AML","AMAS", "TGCM-AML")
aml.samples<-coverage[coverage[,"Project"] %in% aml.projects,"Sample"]
length(aml.samples)
aml.have<-sample.sheet.full[,"ParticipantCode"] %in% aml.samples
table(sample.sheet.full[,"SampleProject"])
table(sample.sheet.full[aml.have,"SampleProject"])
table(sample.sheet.full[!aml.have,"SampleProject"])
a.control<-sample.sheet.full[,"SampleProject"]=="Control"
#sample.sheet.full[aml.have & a.control,]


## cg.samples<-read.delim("/media/UQCCG/Sequencing/CompleteGenomics/917_Data_Delivery_091214.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
## cg.samples[1:5,]
## recodes<-read.delim("/media/UQCCG/Sequencing/CompleteGenomics/RecodingSampleID.csv",header=T,sep=",",fill=TRUE,stringsAsFactors=FALSE) 
## recodes[1:5,]

## analysis.samples<-read.delim("/media/UQCCG/Sequencing/CompleteGenomics/regions_Jonathan/rosetta.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
## analysis.samples[1:5,]


posns<-match(sample.sheet.full[,"ParticipantCode"],coverage[,"Sample"])
missing<-is.na(posns)
sum(missing)
## sample.sheet.full[posns[missing],"ParticipantCode"]
## coverage[,"Sample"]
## sample.sheet.full[missing,]
## coverage<-coverage[posns,]
colnames(coverage)
Capture<-coverage[posns,"Capture.Method"]

cbind(sample.sheet.full,Capture)

############## ASSUME HERE THAT SampleProjuect has classes Control and AML
pheno.types<-c("SampleProject")
names(pheno.types)<-c("SampleProject")

table(sample.sheet.full[,pheno.types[1]])

case.control<-c("SampleProject")
case.control.classes<-c(0,1)
names(case.control.classes)<-c("Control","AML")
case.control.classes
# ib<-1
for(ib in 1:length(case.control)){
  if(!(case.control[ib] %in% colnames(sample.sheet.full))){next}
  sample.sheet.full[(  !(sample.sheet.full[,case.control[ib]] %in% names(case.control.classes))  |  is.na(sample.sheet.full[,case.control[ib]]) | sample.sheet.full[,case.control[ib]]==0 | sample.sheet.full[,case.control[ib]]==9)  ,case.control[ib]]<-NA
}


#tapply(sample.sheet.full[,"SampleProject"],sample.sheet.full[,"SampleProject"],length)
sample.sheet.full[1:5,]
control.samples<-{}

colnames(sample.sheet.full)[colnames(sample.sheet.full)=="ParticipantCode"]<-"SAMPLE"
all.samples<-sample.sheet.full[,"SAMPLE"]


length(all.samples)
## o.remove.all<-expand.labels.to.samples(remove.from.all.samples,all.samples)
## to.remove.samples<-unique(to.remove.all)
## remove.cols<-unique(c(remove.cols,to.remove.samples))

# fam  
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
if(fam=="ALL" | fam=="All" | fam=="all"){
  fam<-unique(unlist( mapply(function(x){x[length(x)]}, strsplit(gsub(the.extension,"",files),split=".",fixed=TRUE)   )))
}

fam
## ifam<-2

## geno.aogc<-read.delim("/media/scratch2/AOGC-NGS/Analysis/AOGC-Genotyping.output.AOGC_ALL.geno.all.txt",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
## setwd("/media/Bioinform-D/Research/annovar/humandb/")
## save(list=c("geno.aogc","use.key"),file="aogc.count.data.RData")

#load("/media/Bioinform-D/Research/annovar/humandb/aogc.count.data.RData")
#print(colnames(indels))
## print(colnames(geno.aogc))
## use.key<-build.key(geno.aogc,core.ann)
## insert.location<-70 
### this is where to add AOGC data
################################ add aogc
## project.files

#
ifam<-1
#for(ifam in 1:length(fam)) {

the.extension<-paste(fam[ifam],project.extension,"$",sep="")
project.files<-files[grepl(the.extension ,files)]
print(sort(paste("Doing: ",project.files,sep=""))) # project.files<-project.files[1:22]

indels<-{}
the.col<-{}
project.files

ichr<-1

#for(ichr in 1:length(project.files)) { ### loop over chromosomes

setwd(analysis.dir)

## grep("Gene.Name",a.indel[,16])
## save(list=c("column.labels"),file="column.labels.RData")
## load("column.labels.RData")
################## fast read ###########
column.labels<-read.delim("P:\\PhD\\Exome Capture Data and Mutations\\Raw Sequencing Data\\AML with Controls - FANC Set\\2013-10-27_AML_with_AOGCControl_NoFailedLane.TGCM-AML.All-maf-filtered.txt",header=F,nrows=1,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
num.vars<-dim(column.labels)[2]
# read in file of massive data - modify
a.indel<-scan("P:\\PhD\\Exome Capture Data and Mutations\\Raw Sequencing Data\\AML with Controls - FANC Set\\2013-10-27_AML_with_AOGCControl_NoFailedLane.TGCM-AML.All-maf-filtered.txt",what=character(num.vars),skip=1,sep="\t",fill=TRUE,na.strings="",quote="\"")
num.lines<-length(a.indel)/(num.vars)
dim(a.indel)<-c(num.vars,num.lines)
a.indel<-t(a.indel)
colnames(a.indel)<-column.labels
########################################
the.chr<-unique(a.indel[,"chr"])
print(paste("Doing Chromosome ",the.chr))

if(sum(!grepl("^chr",the.chr))>0) {
  a.indel[,"chr"]<-paste("chr",a.indel[,"chr"],sep="")
}

a.indel[1:50,"chr"]
key<-build.key(a.indel,core.ann)
rownames(a.indel)<-key
rownames(a.indel)[1:5]
#tapply(gene.list[,"CHR"],gene.list[,"CHR"],length)
all.possible.samples<-gsub(".GT$","",colnames(a.indel)[grep(".GT$",colnames(a.indel))],perl=TRUE)
length(all.possible.samples)
pheno.types

## bad.samples
## colnames(a.indel)[grep("8014",colnames(a.indel))]
## dim(a.indel)

########################### REMOVE BAD SAMPLES HERE
## analysis.samples[1:5,]
## bad.samples<-as.character(c(1:101))
## bad.samples.labels<-expand.labels.to.samples(bad.samples,c("GT","AD","DP","GQ"),paste.after=TRUE)
## a.indel<-a.indel[,colnames(a.indel)[!(colnames(a.indel) %in% bad.samples.labels)]]

## all.possible.samples

## colnames(a.indel)[1:20]
a.indel[1:5,1:10]

col<-grepl("REF",a.indel[,"REF"])
a.indel[col,1:15]
a.indel<-a.indel[!col,]


colnames(a.indel)[1:50]
dim(a.indel)
sort(unique(a.indel[,"chr"]))
sort(unique(a.indel[,"FILTER"]))

test<-"AOGC-08-0287.GT"
wanted<-a.indel[,test]!="0/0" | is.na(a.indel[,test])
sum(wanted)
a.indel[wanted,c("chr","start", "end","Consequence.Embl",test)][1:10,]

sort(table(a.indel[wanted, "Consequence.Embl"]))
sort(table(a.indel[wanted, "TYPE"]))
test<-all.possible.samples
#ALT.Alleles.Control

#################################### got some missing gene names still.
## all.genes[grep("GTPBP4",names(all.genes))]
no.gene<-is.na(a.indel[,"Gene.Names"]) | a.indel[,"Gene.Names"]=="NA"
all.genes<-sort(table(a.indel[,"Gene.Names"]),decreasing=TRUE)

ens.to.hgnc<-read.delim("P:\\PhD\\Exome Capture Data and Mutations\\Raw Sequencing Data\\AML with Controls - FANC Set\\ENSG_to_HGNC.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
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
##    MUC4       TTN     MUC16     MUC5B     MUC12   FAM182B    ZNF717     NBPF1   FAM182A    AHNAK2 
##    1357       540       366       347       326       289       284       283       279       269 
##   MST1L     HYDIN       FLG     SYNE1    NBPF14 POM121L8P       NEB      MUC6 POM121L1P    MST1P2 
##     267       255       249       238       233       213       212       211       204       198 
## PDE4DIP  HLA-DRB1     CTBP2     MUC17     OBSCN  FLJ45445     USH2A     LRP1B      SSPO    DNAH11 
##     196       192       186       184       184       178       175       172       172       167

##  MUC16       MUC4    ANKRD36        TTN   HLA-DRB1      SYNE1      FRG1B     ZNF717        NEB 
##   1176       1171        849        726        712        642        598        576        546 
##   FRG1      MUC12   HLA-DRB5  ANKRD30BL      CSMD1       RYR1      CDC27     DNAH17    FAM182B 
##    487        470        467        450        430        429        398        390        378 
## DNAH11      KMT2C      PKHD1       SSPO ANKRD20A8P      FCGBP      CTBP2      MUC5B      DNAH5 
##    376        364        362        355        354        347        342        341        337 
##  GPR98   MIR1273H      LAMA5 
##    334        329        328 

grep("NOTCH1",names(all.genes))
common.hit.genes<-names(all.genes)[1:10] # common.hit.genes<-names(all.genes)[1:4] - top 10 most ocmmon genes
all.genes["SCN2A"]

###############################################

# tells skat what variants belong to which gene and which genes belong to which pathway
snpinfo.raw<-cbind(key,a.indel[,"Gene.Names"],a.indel[,"Gene.Names"])
snpinfo.raw[1:5,]
tail(snpinfo.raw)

colnames(snpinfo.raw)<-c("Name","gene","cluster")
dim(snpinfo.raw)
snpinfo.raw[1:5,]
dim(snpinfo.raw)
snpinfo<-snpinfo.raw

## other.clusters<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/Final_FANC_clusters.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
## other.clusters<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/other_clusters.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
## other.clusters[1:5,]
## colnames(other.clusters)
## for(ic in 1:dim(other.clusters)[2]){
##    cluster.genes<-clusters[clusters[,"Cluster"]==clusters.wanted[ic],"Genes.assigned.to.group"]
##    cluster.genes<-unlist(strsplit(cluster.genes,split=", "))
##    cluster.genes
##    last.cluster.length<-length(cluster.genes)
##    if(ic==1){clinical.genes<-cluster.genes}
##    if(ic==2){fanc.genes<-cluster.genes}
##   snpinfo[ snpinfo[,"gene"] %in%  cluster.genes  ,"cluster"] <-clusters.wanted[ic]

##  }


#  /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/other_clusters.csv
#clusters<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/Clusters Definitions.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
### CLUSTERS TO CHANGE
clusters<-read.delim("P:\\PhD\\Exome Capture Data and Mutations\\Raw Sequencing Data\\AML with Controls - FANC Set\\Final_FANC_clusters.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
clusters

#clusters.wanted<-c("Clinical","FANC - ACID")
clusters.wanted<-colnames(clusters)
ic<-1
snpinfo[1:5,]

cbind(unique(clusters[,22]),unique(clusters[,22]))
snpinfo[1:5,]

gene.aliases<-read.delim("P:\\PhD\\Exome Capture Data and Mutations\\Raw Sequencing Data\\AML with Controls - FANC Set\\Gene_symbol_aliases.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
gene.aliases[1:5,]
gene.aliases<-unwind(gene.aliases, "Aliases",delimit=", ")
gene.aliases[1:5,]
###############  all.genes

gene.aliases[grepl("SLX1",gene.aliases[,1]),]
gene.aliases[grepl("SLX1",gene.aliases[,2]),]
############FIX MISSING gene names in pathways tests
#Check clusters have genes that are missing from the sequencing list then make sure have the approved gene name


ic<-5
recode<-{}
for(ic in 1:length(clusters.wanted)){
  #  print(ic)
  cluster.genes<-clusters[,clusters.wanted[ic]]
  cluster.genes<-cluster.genes[cluster.genes!="" | is.na(cluster.genes)]
  cluster.genes<-unique(cluster.genes)
  #   all.genes[1:5]
  missing.name<-!(cluster.genes %in% names(all.genes))
  
  if(sum( missing.name)>0) {
    posns<-match(cluster.genes[missing.name],gene.aliases[, "Aliases"])
    missing<-is.na(posns)
    
    if(sum(missing)>0) {
      print(paste("in cluster",clusters.wanted[ic], " missing"))
      print(cluster.genes[missing.name][!missing])
    }
    
    recode<-cbind(cluster.genes[missing.name][!missing],gene.aliases[posns[!missing], "Approved.Symbol"])
    colnames(recode)<-c("old","new")
    ###### transfer to new gene lists
    posns<-match(clusters[,clusters.wanted[ic]],recode[,"old"])
    missing<-is.na(posns)
    clusters[!missing,clusters.wanted[ic]]<-recode[posns[!missing],"new"] ### redefine the clusters
  }
}
#########################################################    

snpinfo<-snpinfo.raw
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

## HAVE
# snpinfo.raw (original from a.indel)
# snpinfo # with extra clusters
# snpinfo.raw a permanent copy of snpinfo

## "FANCM " "MHF1"   "MHF2"   "FAAP24"
clusters[,1]
chk<-apply(clusters,2,function(x){ length(x[x!=""])})
write.table(clusters,file="clusters_as.using.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

############################
a.indel[1:20,1:5]

the.samples<-colnames(a.indel)[grepl(".GT$", colnames(a.indel))]

filtered.genotype(a.indel[1:20,],gsub(".GT$","",the.samples),prefix="",suffix="",20,0.02,0.98,0.20,0.80,7,2)

num.cores<-4
num.bits<-num.cores
##registerDoMC(cores=num.cores)
registeredCores <- makeCluster(num.cores)
registerDoParallel(registeredCores)

colnames(a.indel)[grepl(".GQ$", colnames(a.indel))]
while((dim(a.indel)[1] %% num.bits)< 2) {
  num.bits<-num.bits+1
} ### go don't get matrix issues

# modify if cant use multicore
#fil.genotypes<-foreach(a.indel.bit=iter(a.indel,by='row',chunksize=as.integer(dim(a.indel)[1]/num.bits) ), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype(a.indel.bit,gsub(".GT$","",the.samples),prefix="",suffix="",20,0.02,0.98,0.20,0.80,7,2)
fil.genotypes <- filtered.genotype(a.indel,gsub(".GT$","",the.samples),prefix="",suffix="",20,0.02,0.98,0.20,0.80,7,2)
# 20,0.02,0.98,0.20,0.8,10,5 # for cancers where het may be way out
## length(col)
## fil.genotypes[col,1:5]
## fil.genotypes<-fil.genotypes[!col,]
# 0.02< het(0/1)<0.98
# ref(0/0)< 0.2
# (1/1) > 0.8

dim(fil.genotypes)
colnames(fil.genotypes)[1:5]
rownames(fil.genotypes)[1:5]

dim(fil.genotypes)
dim(a.indel)

tail(rownames(a.indel))

## a.indel<-a.indel[1:342791,]
## fil.genotypes<-fil.genotypes[1:342791,]
############################### do one phenotype #################### 
#
ipheno<-1
# for(ipheno in 1:length(pheno.types)){

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

# what send to skatmeta
formula<-paste(target.pheno.col,"~",paste(covars,collapse="+"),sep="")
print(formula)
formula<-formula(formula)

#seq.type[1:57,]
posns<-match(pheno[,"SAMPLE"],seq.type[,"Sample"])
missing<-is.na(posns)
sum(missing)

capture<-seq.type[posns,"Capture.Method"]

## the.runs<-seq.type[posns,"Run"]

## the.run<-the.runs
## the.run[the.runs=="C1F6GACXX"]<-0
## the.run[the.runs!="C1F6GACXX"]<-1
##  the.run<-as.numeric(the.run)

## names(the.runs)<-the.runs
## unique.runs<-unique(the.runs)
## unique.runs.fill<-c(1:length(unique.runs))
## names(unique.runs.fill)<-unique.runs
## the.run<-unique.runs.fill[the.runs]
## the.run<-as.numeric(the.run)

table(capture) ### all illume here
pheno<-cbind(pheno,capture)
pheno[1:5,]
pheno<-pheno[pheno[,"SAMPLE"] %in% all.possible.samples,]

the.samples<-paste(pheno[,"SAMPLE"],"GT",sep=".")  ## samples same order as in pheno
print(paste("Number samples: ",length(the.samples),sep=""))

table(pheno[,"SampleProject"])

## ###############GEFOS SNP type restrictions

######################## make sure pheno has sample samples as A.indel

##### AFFSTAT_IN_WORDS is mostly correct to the regression #####
pheno[1:5,]

summary.geno.extra<-{}

target<-"AML"
use.samples<-the.samples[pheno[,"SampleProject"]=="AML"]
length(use.samples)
genotypes<-a.indel[,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")
#summary.geno[1:5,]
if(is.null(dim(summary.geno.extra))) {
  summary.geno.extra<-summary.geno
} else {
  summary.geno.extra<-cbind(summary.geno.extra,summary.geno)
}

summary.geno.extra[1:5,]

target<-"Control"
use.samples<-the.samples[pheno[,"SampleProject"]=="Control"]
length(use.samples)
genotypes<-a.indel[,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")
#summary.geno[1:5,]

if(is.null(dim(summary.geno.extra))) {
  summary.geno.extra<-summary.geno
} else {
  summary.geno.extra<-cbind(summary.geno.extra,summary.geno)
}

target<-"AML.filt"
use.samples<-the.samples[pheno[,"SampleProject"]=="AML"]
length(use.samples)
genotypes<-fil.genotypes[,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")
#summary.geno[1:5,]
if(is.null(dim(summary.geno.extra))) {
  summary.geno.extra<-summary.geno
} else {
  summary.geno.extra<-cbind(summary.geno.extra,summary.geno)
}

target<-"Control.filt"
use.samples<-the.samples[pheno[,"SampleProject"]=="Control"]
length(use.samples)
genotypes<-fil.genotypes[,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")
#summary.geno[1:5,]
if(is.null(dim(summary.geno.extra))) {
  summary.geno.extra<-summary.geno
} else {
  summary.geno.extra<-cbind(summary.geno.extra,summary.geno)
}

a.indel[1:5,1:10]
summary.geno.extra[1:5,]
#colnames(summary.geno.extra)
#rownames(summary.geno.extra)<-key
#getHWE(obs_hets, obs_hom1, obs_hom2)

hw.p.control<-getHWE(summary.geno.extra[,"GENO.Control"]) ## used 16 CPUs
hw.p.control.filt<-getHWE(summary.geno.extra[,"GENO.Control.filt"]) ## used 16 CPUs

## class(summary.geno.extra)
hw.p.control[1:5]
length(hw.p.control)
names(hw.p.control)<-key
names(hw.p.control.filt)<-key
######### testing
## hw.p.control[200100:200195]
## hw.p.control[22:25]
## summary.geno.extra[22:25,"GENO.Control"]

## test.loc<-c(24,200129:200143)
## geno<-strsplit(summary.geno.extra[test.loc,"GENO.Control"],split=",")
## AA<-as.numeric(unlist(lapply(geno,function(x) x[1])))
## BB<-as.numeric(unlist(lapply(geno,function(x) x[3])))
## AB<-as.numeric(unlist(lapply(geno,function(x) x[2])))
## x<-cbind(AA,AB,BB)
## x
## HWExactMat(x)
## hw.p.control[test.loc] ## correct in original order

######################## inbreeding and no genotypes

group.maf.thresh<-0.20 # group.maf.thresh<-0.10  if more common that this in Group then discard: Discard Rare but common in this cohort ## not a point level
missing.threshold<-0.20 # missing.threshold<-0.50  60 % of genotypes missing
## missing.threshold.nimblgen<-0.20
## missing.threshold.illumina<-0.20 

## rare.in.group<-( summary.geno.extra[,c("MAF.HIGH","MAF.LOW")]< group.maf.thresh) | ( summary.geno.extra[,c("MAF.HIGH","MAF.LOW")] > (1-group.maf.thresh)) | (is.na( summary.geno.extra[,c("MAF.HIGH","MAF.LOW")]))
rare.in.group.table<-( summary.geno.extra[,c("MAF.AML","MAF.Control","MAF.AML.filt","MAF.Control.filt")]< group.maf.thresh) | ( summary.geno.extra[,c("MAF.AML","MAF.Control","MAF.AML.filt","MAF.Control.filt")] > (1-group.maf.thresh)) | (is.na( summary.geno.extra[,c("MAF.AML","MAF.Control","MAF.AML.filt","MAF.Control.filt")]))
rare.in.group.table[1:5,]
#rare.in.group<-cbind( indels[,c("MAF.HIGH","MAF.LOW")]< group.maf.thresh) , (indels[,c("MAF.HIGH","MAF.LOW")] > (1-group.maf.thresh)) , (is.na(indels[,c("MAF.HIGH","MAF.LOW")])))
rare.in.group<-combine.boolean(rare.in.group.table,c("MAF.AML","MAF.Control"),"OR")
rare.in.group.filt<-combine.boolean(rare.in.group.table,c("MAF.AML.filt","MAF.Control.filt"),"OR")
sum(!rare.in.group)
sum(!rare.in.group.filt)

#no.genotypes<-(indels[,c("MAF.HIGH","MAF.LOW")]== 0)  | (is.na(indels[,c("MAF.HIGH","MAF.LOW")])) # no genotypes in test classes for a mutataion after individaul quality filtering

### checking if BOTH not mono mrphic
no.genotypes.table<-(summary.geno.extra[,c("MAF.AML","MAF.Control","MAF.AML.filt","MAF.Control.filt")]== 0) | summary.geno.extra[,c("MAF.AML","MAF.Control","MAF.AML.filt","MAF.Control.filt")]=="NaN"  | (is.na( summary.geno.extra[,c("MAF.AML","MAF.Control","MAF.AML.filt","MAF.Control.filt")])) # no genotypes in test classes for a mutataion after individaul quality filtering
no.genotypes.table[1:5,]
no.genotypes<-combine.boolean(no.genotypes.table,c("MAF.AML","MAF.Control"),"AND") # was AND
no.genotypes.filt<-combine.boolean(no.genotypes.table,c("MAF.AML.filt","MAF.Control.filt"),"AND") # was AND
sum(no.genotypes)
sum(no.genotypes.filt)

summary.geno.extra[1:5,]
missing.targets<-c("AML","Control","AML.filt","Control.filt")
high.missing<-{}
imt<-1

for(imt in 1:length(missing.targets)) {
  the.missing.alleles<-paste("MISSING.Alleles",missing.targets[imt],sep=".")
  the.total.alleles<-paste("TOTAL.Alleles",missing.targets[imt],sep=".")
  a.missing.test<-as.numeric(summary.geno.extra[,the.missing.alleles])/(as.numeric(summary.geno.extra[,the.total.alleles])+as.numeric(summary.geno.extra[,the.missing.alleles]))
  
  if(is.null(length(high.missing)) | length(high.missing)==0){
    high.missing<-a.missing.test
  } else {
    high.missing<-cbind(high.missing,a.missing.test)
  }
}
colnames(high.missing)<-missing.targets
rownames(high.missing)<-key

## nimblegen.total.missing<-subset(high.missing,select=c("nimblegen"))
## nimblegen.total.missing[1:5,]
## nimblegen.total.missing<-nimblegen.total.missing > missing.threshold.nimblgen
## ## nimblegen.total.missing<-combine.boolean(high.total.missing,c("LOW","HIGH","LOW.pheno","HIGH.pheno"),"OR")
## sum(nimblegen.total.missing)


## illumina.total.missing<-subset(high.missing,select=c("illumina"))
## illumina.total.missing[1:5,]
## illumina.total.missing<-illumina.total.missing > missing.threshold.illumina
## ## nimblegen.total.missing<-combine.boolean(high.total.missing,c("LOW","HIGH","LOW.pheno","HIGH.pheno"),"OR")
## sum(illumina.total.missing)

## sum(high.total.missing | nimblegen.total.missing | illumina.total.missing)
high.missing[1:5,]
#high.missing[places,]
#ok.missing.test[places,]
ok.missing.test<-high.missing <= missing.threshold
ok.missing.test[1:5,]
ok.missing<-combine.boolean(ok.missing.test,c("AML","Control"),"AND") # was AND
ok.missing.filt<-combine.boolean(ok.missing.test,c("AML.filt","Control.filt"),"AND") # was AND
sum(!ok.missing)
sum(!ok.missing.filt)

## ok.missing[places]
## diff<-ok.missing & !ok.missing.filt
## ####testing
## sum(diff)

## grep(TRUE,diff)[1:10]
## target<-"chr10:93396:93396:C:T:snp"
## target<-key[diff][83:84] # target<-key[test]
## high.missing[target,]
## chk<-gsub(".GT$",".AD",the.samples[1:96])
## chk.GT<-the.samples[1:96]
## a.indel[target,chk]
## a.indel[target,chk.GT]
## fil.genotypes[target,chk.GT]
## summary.geno.extra[target,paste("GENO",c("AML","Control","AML.filt","Control.filt"),sep=".")]
## hw.p.control.filt[target]

#test<-pass & !rare.in.group
## > summary.geno.extra[target,paste("GENO",c("AML","Control","AML.filt","Control.filt"),sep=".")]
##                                            GENO.AML   GENO.Control GENO.AML.filt GENO.Control.filt
## chr11:6579106:6579106:C:A:snp:6579106      "9,52,35"  "27,97,76"   "9,48,13"     "26,87,50"       
## chr11:7022531:7022531:A:T:snp:7022531      "49,34,13" "74,107,19"  "49,34,0"     "74,97,2"        
## chr12:11183451:11183451:A:T:snp:11183451   "53,39,4"  "120,67,13"  "53,31,0"     "120,63,0"       
## chr16:84808824:84808824:C:G:snp:84808824   "37,45,14" "70,97,33"   "36,44,3"     "70,96,4"        
## chr17:21216964:21216964:T:C:snp            "24,72,0"  "54,143,0"   "24,64,0"     "39,59,0"        
## chr17:72927123:72927123:G:T:snp:72927123   "42,39,15" "85,95,18"   "42,39,2"     "83,81,11"       
## chr19:13135824:13135824:G:T:snp            "0,41,55"  "5,81,112"   "0,41,48"     "0,30,99"        
## chr19:55285045:55285045:G:T:snp            "2,44,50"  "5,75,118"   "2,41,46"     "1,53,83"        
## chr19:55286796:55286796:G:A:snp            "8,50,38"  "5,88,105"   "8,50,38"     "4,86,103"       
## chr19:55286854:55286854:A:G:snp            "5,53,38"  "9,86,103"   "4,52,37"     "8,75,95"        
## chr1:13183439:13183439:T:C:snp             "52,44,0"  "96,103,0"   "52,44,0"     "96,103,0"       
## chr1:18808526:18808526:A:C:snp:18808526    "38,40,18" "69,87,36"   "37,40,2"     "54,60,20"       
## chr1:248457979:248457979:C:A:snp:248457979 "53,39,4"  "93,92,14"   "53,39,0"     "93,91,0"        
## chr21:15013735:15013735:A:G:snp            "39,39,18" "54,104,42"  "39,39,18"    "54,104,41"      
## chr21:43541342:43541342:C:T:snp:43541342   "38,39,19" "67,106,27"  "38,39,0"     "67,105,4"       
## chr4:187538942:187538942:T:G:snp:187538942 "27,54,15" "73,94,33"   "27,52,0"     "73,87,2"        
## chr6:31237833:31237833:T:C:snp:31237833    "58,28,10" "97,76,27"   "54,28,3"     "86,60,14"       
## chr6:31238259:31238259:G:T:snp:31238259    "56,32,8"  "97,79,24"   "52,30,2"     "80,57,22"       
## chr7:87160618:87160618:A:C:snp:87160618    "28,44,24" "61,87,52"   "28,44,2"     "61,81,26" 

################## get loci in common
## load("/media/Bioinform-D/Research/AML sequencing/Human_Exome_Targets_illumina_v2_hg19_targets.RData")
## ill.gr<-data.gr
## load("/media/Bioinform-D/Research/AML sequencing/Human_Exome_Targets_Nimble_v2_hg19_targets.RData")
## nim.gr<-data.gr
## genome(ill.gr)<-"hg19"
## genome(nim.gr)<-"hg19"
## ill.gr<-ill.gr+200
## nim.gr<-nim.gr+200
## #overlaps<-countOverlaps(ill.gr,nim.gr)
## overlaps<-overlapsAny(ill.gr,nim.gr)
## possible.loci<-ill.gr[overlaps]

## my.loci<-IRanges(start=as.numeric(a.indel[,"start"]),end=as.numeric(a.indel[,"end"]))
## on.chr<-as.vector(seqnames(possible.loci.ori))==a.indel[1,"chr"] # paste("chr",a.indel[1,"chr"],sep="")
## possible.space<-IRanges(start=as.numeric(start(possible.loci.ori))[on.chr],end=as.numeric(end(possible.loci.ori))[on.chr])
## common.loci<-overlapsAny(my.loci,possible.space)
## sum(common.loci)
## overlaps<-findOverlaps(possible.space,my.loci)
## snpinfo<-cbind(snp.names[subjectHits(overlaps)],gene.on.chr[queryHits(overlaps),"GENE_NAME"])
## colnames(snpinfo)<-c("Name","gene")
###############################################################################################################################################
################# HOWEVER ALT ALLELES ALL IS USED AS MINOR ALLELE FREQUNCY  not alternative allele frequencies

maf.lt.all<-a.indel[,colnames(a.indel)[grepl("MAF.lt:",colnames(a.indel))]]
maf.lt.all[1:5,]
#maf.lt.all<-as.logical(maf.lt.all)
as.logical(maf.lt.all[1:5,5])

####################################################################################
##################################### make the POSITION filter matrix QUALITY.FILTER FUNATIOAL group-MAF done here
global.labs[!(global.labs %in% c(colnames(a.indel),colnames(summary.geno.extra))  )]
if(sum( !(global.labs %in% c(colnames(a.indel),colnames(summary.geno.extra)))  )>0) {
  print(paste("WARNING postion filters missing for data",
              toString(global.labs[!(global.labs %in% c(colnames(a.indel))) ])))
}
#grep("ljb_gerp",colnames(a.indel))
##############$$$$$$$$$$$$$$ POSITION FILTER: combines a.indel, filter table and filter.table.pholy in this case.
qual<-position.quality.filter(  cbind( a.indel[,colnames(a.indel) %in% global.labs],summary.geno.extra[,colnames(summary.geno.extra) %in% global.labs] ), global.quality.cut,global.quality.type,global.quality.dirn)
dim(qual)
dim(a.indel)
rownames(qual)<-key
qual[1:5,]
#########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

the.types<-names(tapply(a.indel[,"ensGene::location"],a.indel[,"ensGene::location"],length))
the.types<-unique(unlist(strsplit(the.types,split=";")))
if( sum( !(the.types %in% c("NA","exonic",possible.mutations)))>0  ) {
  print("WARNING ANNOVAR HAS NEW MUTATION TYPES DEFINED- REVIEW line 532")
  print( the.types[!(the.types %in% c("exonic",possible.mutations))] )
}

#  filter.table.pholy[1:5,]
wanted.muts.coding<-test.for.coding.type(a.indel,geneanno.DB,interesting.coding.mutations)
sum(wanted.muts.coding)
####################### mhairi CHECK
wanted.muts.coding.vep<-test.wanted.mutation(a.indel[,"Consequence.Embl"],vep.coding,delimit.by=",")  # filter.table.pholy[,"Consequence"] %in%  vep.coding
## a.indel[wanted.muts.coding.vep,]
sum(wanted.muts.coding.vep)

missense.coding<-test.for.coding.type(a.indel,geneanno.DB,missense.variant)
sum(missense.coding)
####################### mhairi CHECK
missense.coding.vep<-test.wanted.mutation(a.indel[,"Consequence.Embl"],missense.variant,delimit.by=",")  # filter.table.pholy[,"Consequence"] %in%  vep.coding
## a.indel[wanted.muts.coding.vep,]
sum(missense.coding.vep)

is.missense<-missense.coding.vep | missense.coding

qual[1:5,"PolyPhen.low"] # true is above 0.1
a.indel[1:5,"PolyPhen.scores"]
is.benign.missense<-is.missense & !qual[,"PolyPhen.low"]
sum(!is.benign.missense)
dim(a.indel)
## sum(wanted.muts.coding.vep)
##  sum(wanted.muts.coding)
############## noncoding double filtered for biotype ###

interesting.NONcoding.mutations<-interesting.mutations[!(interesting.mutations %in% interesting.coding.mutations)] # "ncRNA_exonic"
wanted.muts.NONcoding<-test.for.coding.type(a.indel,geneanno.DB,interesting.NONcoding.mutations)
sum(wanted.muts.NONcoding)

wanted.interesting.to.prefilter<-test.for.coding.type(a.indel,geneanno.DB,interesting.to.prefilter) #"UTR3"      "UTR5"      "UTR5;UTR3" "snoRNA"
sum(wanted.interesting.to.prefilter)
# unique(a.indel[,"Consequence.Embl"])
####################### mhairi CHECK
wanted.interesting.to.prefilter.vep<-test.wanted.mutation(a.indel[,"Consequence.Embl"],vep.noncoding,delimit.by=",")  #filter.table.pholy[,"Consequence"] %in% vep.noncoding
sum(wanted.interesting.to.prefilter.vep)

wanted.muts.NONcoding.keep<-rep(FALSE,times=dim(a.indel)[1]) # "ncRNA_exonic" and "ncRNA_exonic" 
a.type<-wanted.noncoding.subtypes # list biotypes want to keep
if(!("gene_biotype" %in% colnames(a.indel))) {
  colnames(a.indel)[colnames(a.indel)=="Gene.Biotype"] <- "gene_biotype"
}

for(itype in 1:length(a.type)) {
  the.test<-grepl(a.type[itype],a.indel[,"gene_biotype"])
  the.test[is.na(the.test)] <- FALSE
  wanted.muts.NONcoding.keep<- wanted.muts.NONcoding.keep | the.test
}

#  HERE wanted.muts.NONcoding.keep JUST DENOTES  "miRNA" &  "lincRNA" at this point USE wanted.muts.NONcoding to restrict to exones and splice BELOW
wanted.muts.NONcoding.keep<-wanted.muts.NONcoding.keep & wanted.muts.NONcoding

dim(qual)
colnames(qual)
qual[1:5,]

# things want
bad.coding<-wanted.muts.coding | wanted.muts.coding.vep | (wanted.muts.NONcoding.keep & (qual[,"GERP.low"] | qual[,"GERP.unknown"]) ) | ( (wanted.interesting.to.prefilter | wanted.interesting.to.prefilter.vep) & (qual[,"GERP.high"]) )
sum(bad.coding)

bad.coding<-wanted.muts.coding | wanted.muts.coding.vep 
bad.non.coding<-(wanted.muts.NONcoding.keep & (qual[,"GERP.low"] | qual[,"GERP.unknown"]) ) | ( (wanted.interesting.to.prefilter | wanted.interesting.to.prefilter.vep) & (qual[,"GERP.high"]) )
bad.effect<-bad.coding | bad.non.coding

sum(bad.effect)
sum(bad.non.coding)
sum(bad.coding)

# bad.coding<-test.for.coding.type(a.indel,geneanno.DB,c("stopgain SNV","stoploss SNV","frameshift deletion","frameshift insertion"))
#bad.frame<-test.for.coding.type(geneanno.table,geneanno.DB,c("frameshift deletion","frameshift insertion"))


## basic.qual<-combine.boolean(qual,c("QUAL", "QD", "HRun", "SB"),"AND")
## gatk.qual<-combine.boolean(qual,c("FILTER_PASS", "FILTER_100" ),"OR")
## full.qual<-combine.boolean(cbind(basic.qual,gatk.qual),"all","OR")

## full.qual<-gatk.qual
## sum(full.qual)
full.qual<-qual[,"FILTER_PASS"]
sum(full.qual)
qual[full.qual,c("FILTER_PASS", "FILTER_100")][1:30,]

#any.functional<-combine.boolean(qual,c("PolyPhen.low","SIFT.high","mut.taster.high","phylo.high","PolyPhen.bad","SIFT.bad","GERP.high","ljb_gerp.high"),"OR")
functional<-combine.boolean(qual,c("PolyPhen.low","SIFT.high","PolyPhen.bad","SIFT.bad","GERP.high"),"OR")# include really bad protein changes $$$$ BUT ADD ANY GERP.hiogh might let through a lot of junk
functional.coding<-combine.boolean(qual,c("PolyPhen.low","SIFT.high","PolyPhen.bad","SIFT.bad"),"OR")
sum(functional.coding)
functional<-functional | bad.coding 
sum(functional)

########################  FREQUENCY FILTERS
######### given a 0.01 threshold 6sd would allow 10 alt alleles at 6sd
######## given a 0.005 threshold 6sd would allow 7 alt alleles at 6sd
n<-400 # number of controls
p<-0.001 # maf threshold
## n*p
## ## #np(1-p)# sd =sqrt(var)

## ## ## n<-length(the.samples[pheno[,"SampleProject"]=="Control"])
## ## ## p<-0.01
## ## ## p<-0.005 # maf threshold
## ## sqrt(n*p*(1-p)) #  sd =sqrt(var)

## ## ## #z= 
(13- n*p) / sqrt(n*p*(1-p))      # p=0.01   2=0 7=3.5sd 16=6.03sd
## ## ##                             # p=0.005 2=1 4=3  5=4  6=5 7=6sd
## ##                               # p=0.001 2=4sd 3=6sd  5=4  6=5 7=15sd
## ##                               # p=0.001 2=4sd 3=6sd  5=4  6=5 7=15sd
##                                   #p=0.05  0=10 20=3.2 29=6.16

n<-max(as.integer(summary.geno.extra[,"TOTAL.Alleles.Control"]))
#p<-0.001
p<-0.001 ########### set MAF threshols HERE

sd.thresh<-6
n
p

alt.counts.thresh<-1
while( (alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) <= sd.thresh){alt.counts.thresh<-alt.counts.thresh+1}
alt.counts.thresh

summary.geno.extra[1:5,]
rare.in.controls<-as.numeric(summary.geno.extra[,"ALT.Alleles.Control"])< alt.counts.thresh
rare.in.controls.filt<-as.numeric(summary.geno.extra[,"ALT.Alleles.Control.filt"])< alt.counts.thresh
sum(rare.in.controls)
sum(rare.in.controls.filt)
names(rare.in.controls)<-key
names(rare.in.controls.filt)<-key

length(maf.filter)
length(rare.in.controls)

maf.lt.all[1:5,]
#maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.001"])
#maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.01"])
maf.col<-paste("MAF.lt",p,sep=":")
maf.col
maf.filter<-as.logical(maf.lt.all[,maf.col])
## maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.5"])
sum(maf.filter)
names(maf.filter)<-key
#pass<- rare.in.group & !no.genotypes & !high.missing & common.loci

####################################
################Poly morphic SITE TESTS

REF.length<-nchar(as.character(a.indel[,"REF"]))
ALT.length<-nchar(as.character(a.indel[,"ALT"]))

large.indel<-REF.length>1 | ALT.length>1

are.repeats<-identify.repeats(a.indel,di.run.max=3,homo.run.max=5)

length(large.indel) 
length(are.repeats)
sum(are.repeats)
rownames(a.indel)[are.repeats][1:20]

#################### in repeats looking  forward

chk.in.repeat<-large.indel & !are.repeats
are.sub.repeat<-indentify.IN.repeat(a.indel[chk.in.repeat,],looking="forward",bases.about=6,di.run.max=3,homo.run.max=5,genome="BSgenome.Hsapiens.UCSC.hg19")
remove.repeats<-key[chk.in.repeat][are.sub.repeat]
are.in.repeats.forward<- key %in% remove.repeats

remove.repeats[1:20]
sum(are.in.repeats.forward)
## [1] 6988

###################### in repeats looking back

sum(chk.in.repeat)
chk.in.repeat<-large.indel & !are.repeats & !are.in.repeats.forward

are.sub.repeat<-indentify.IN.repeat(a.indel[chk.in.repeat,],looking="back",bases.about=6,di.run.max=3,homo.run.max=5,genome="BSgenome.Hsapiens.UCSC.hg19")
remove.repeats<-key[chk.in.repeat][are.sub.repeat]
are.in.repeats.back<- key %in% remove.repeats

remove.repeats[1:20]
sum(are.in.repeats.back)
## [1] 3224

are.in.repeats<- are.in.repeats.back | are.in.repeats.forward

length(are.in.repeats)
sum(are.in.repeats)

#are.in.repeats[places]
## test<-a.indel[large.indel,c("chr","start","end","REF","ALT","TYPE")][1:500,]
## save(list=c("test"),file="test.RData")
## are.sub.repeat<-indentify.IN.repeat(test,looking="forward",bases.about=6,di.run.max=3,homo.run.max=5,genome="BSgenome.Hsapiens.UCSC.hg19")
## test[are.sub.repeat,]

##########################################################################
##########################################################################
##########################################################################
not.flat.genotype<-!qual[,"flat"]
sum(not.flat.genotype)

colnames(qual)

is.unwound.geno<-grepl("snp:\\d+$",a.indel[,"TYPE"]) | grepl("indel:\\d+$",a.indel[,"TYPE"])
#a.indel[!not.flat.genotype,"TYPE"]
#grepl("indel:\\d+$",a.indel[places,"TYPE"])
#is.unwound.geno[places]

hwe.control.threshold
hw.p.control.filt[1:50]
hw.p.control[1:50]
hw.controls.ok<-hw.p.control > hwe.control.threshold
hw.controls.ok.filt<-hw.p.control.filt > hwe.control.threshold
sum(hw.controls.ok)
sum(hw.controls.ok.filt)
sum(!hw.controls.ok)
length(hw.controls.ok.filt)
#hw.controls.ok[loci]

in.common.hit.gene <- a.indel[,"Gene.Names"] %in% common.hit.genes
sum(in.common.hit.gene)
the.chr
on.x.y<-a.indel[,"chr"] %in% c("X","Y","23","24","chrX","chrY")
sum(on.x.y)

#table(a.indel[,"TYPE"])
snp.only<-grepl("^snp",a.indel[,"TYPE"])
## bad.coding<-wanted.muts.coding | wanted.muts.coding.vep
##pass<-full.qual & bad.coding & maf.filter & rare.in.group & !no.genotypes  & not.flat.genotype & rare.in.controls
#pass<-full.qual & bad.effect & maf.filter & rare.in.group & !no.genotypes & !in.common.hit.gene  & hw.controls.ok & !on.x.y & !unannotated.hits & not.flat.genotype & !are.repeats

#ok.missing.filt ok.missing snp.only &

## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene  & !on.x.y & !unannotated.hits & not.flat.genotype & !are.repeats & !are.in.repeats &
## ( ok.missing.filt | is.unwound.geno)   & hw.controls.ok.filt & !no.genotypes.filt   & rare.in.controls.filt & rare.in.group

## sum(pass)
## length(pass)

## key[test]

## sum(maf.filter)
## sum(pass)
## length( rare.in.group)
## size.no.extras<-length(pass)

## sum(functional.coding)
##     sum(pass) 
##     sum(full.qual) 
##     sum(bad.effect)
##    sum(bad.coding) 
##     sum(bad.non.coding) 
##     sum(maf.filter) 
##     sum(!rare.in.group)
##   sum(!rare.in.group.filt) 
##     sum(no.genotypes)
##   sum(!no.genotypes.filt) 
##     sum(in.common.hit.gene) 
##     sum(!not.flat.genotype) 
##     sum(!hw.controls.ok)
##    sum(hw.controls.ok.filt)
##     sum(on.x.y) 
##     sum(unannotated.hits) 
##     sum(are.repeats) 

##     sum(summary.geno.extra) 
##      sum(!ok.missing) 
##     sum(are.in.repeats)

## sum(ok.missing.filt)
## sum(hw.p.control.filt)
## sum(rare.in.group.filt)
## sum(no.genotypes.filt)
## sum(rare.in.controls.filt)
## #& not.flat.genotype

## sum(are.repeats | are.in.repeats)

##########################################################################
##########################################################################

######## print out clusters hack
## combine<-clusters[,1]
## for(i in 2:13){

## combine<-c(combine,clusters[,i])

## }

## combine<-unique(combine)

## cbind(combine,combine)

icc<-2
if(target.pheno %in% case.control) {
  
  for(icc in 1:length(case.control.classes)) {
    recode<-  pheno[,target.pheno] %in% names(case.control.classes)[icc]
    pheno[recode,target.pheno]<-as.numeric(case.control.classes[icc])
  }
}

pheno[,target.pheno]<-as.numeric(pheno[,target.pheno])
formula

pheno[1:5,]
pheno[,target.pheno]

#no.pheno.samples<-is.na(pheno[,target.pheno]) # these already removed

## pass<-full.qual & functional & maf.filter & rare.in.group & !no.genotypes  & not.flat.genotype & !(high.total.missing | nimblegen.total.missing | illumina.total.missing)
## sum(pass)

## 16      3639869 3639869 A       C       snp     nonsynonymous SNV       SLX4:NM_032444:exon12:c.T3770G:p.V1257G, (FANCP)
## 17      41197708        41197708        T       G       snp     nonsynonymous SNV       BRCA1:NM_007300:exon24:c.A5642C:p.H1881P,BRCA1
## 3       10088285        10088285        T       G       snp     nonsynonymous SNV       FANCD2:NM_033084:exon15:c.T1156G:p.F386V
## 3       10088343        10088343        A       G       snp     nonsynonymous SNV       FANCD2:NM_033084:exon15:c.A1214G:p.N405S,
## 3       10088412        10088412        G       T       snp     splicing                        FANCD2(NM_033084:exon15:c.1278+5G>T,NM_001018115:exon15:c.1278+5G>T)    
## 3       10089644        10089644        G       A       snp     nonsynonymous SNV       FANCD2:NM_033084:exon16:c.G1322A:p.S441N,
## 3       10089738        10089738        A       G       snp     splicing                        FANCD2(NM_033084:exon16:c.1413+3A>G,NM_001018115:exon16:c.1413+3A>G)    
## 3       10105570        10105570        A       G       snp     nonsynonymous SNV       FANCD2:NM_033084:exon21:c.A1922G:p.H641R,FANCD2:
## 3       10106408        10106408        C       T       snp     splicing                        FANCD2(NM_033084:exon23:c.2022-5C>T,NM_001018115:exon23:c.2022-5C>T)    
## 3       10106532        10106532        C       T       snp     nonsynonymous SNV       FANCD2:NM_033084:exon23:c.C2141T:p.P714L,
## 3       10114944        10114944        A       C       snp     nonsynonymous SNV       FANCD2:NM_033084:exon28:c.A2613C:p.K871N,
## 3       10128939        10128939        G       A       snp     nonsynonymous SNV       FANCD2:NM_033084:exon34:c.G3457A:p.E1153K,
## 3       10088407        10088410        AGTA    -       indel   frameshift deletion     FANCD2:NM_033084:exon15:c.1278_1278del:p.426_426del,

#ENSG00000269323 # junked
ignore.FLT3 <- TRUE
if( ("13" %in% the.chr) & !("chr13:28626716:28626716:C:T:CREST" %in% key) & !ignore.FLT3 ){ 
  # add flt3-ITD in chr213 and not already there
  
  #    /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/extra_FANC.csv
  ## extra<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/extra_FANC.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
  ## if(sum(!grepl("^chr",extra[1,"chr"]))>0){
  ## extra[,"chr"]<-paste("chr",extra[,"chr"],sep="")
  ## }
  ##     key.extra<-build.key(extra,core.ann)
  
  ##    figure<- match(key.extra,key)
  ## pass[figure]
  ## help[figure,]
  
  ## colnames(a.indel)[1:50]
  
  ## key[grep("chr17",key)[1:100]]
  ## grep("chr17:41197708",key)
  ## key[grep("10088407",key)]
  ## out<-cbind(a.indel[figure,c(1:6,16,30,34,37:42)],summary.geno.extra[figure,],help[figure,])
  ## colnames(out)
  ## out[1:5,]
  ## write.table(out,file="extra.summary.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
  
  extra<-a.indel[1,]
  extra[the.samples]<-"0/0"
  itd.pos<-paste(c(6,7,12,17,18,22,23,25,26,30,33,42,46:50,54,64,66,75,76,84,91,93,99),"GT",sep=".")
  itd.pos
  
  #### Make up some dummy genotypes
  extra[itd.pos]<-"0/1"
  extra.geno<-extra[the.samples]
  
  the.samples.AD<-gsub(".GT$",".AD",the.samples)
  extra[the.samples.AD]<-"100,0"
  
  the.samples.DP<-gsub(".GT$",".DP",the.samples)
  extra[the.samples.DP]<-"100"
  
  itd.pos<-paste(c(6,7,12,17,18,22,23,25,26,30,33,42,46:50,54,64,66,75,76,84,91,93,99),"AD",sep=".")
  itd.pos
  
  #### Make up some dummy genotypes
  extra[itd.pos]<-"50,50"
  
  extra[core.ann]<-c("chr13","28626716","28626716","C","T","CREST")
  extra["Gene.Names"]<-"FLT3"
  
  target<-"AML"
  use.samples<-the.samples[pheno[,"SampleProject"]=="1"]
  length(use.samples)
  genotypes.extra<-t(as.matrix(extra[use.samples]))
  dim(genotypes.extra)
  aml.extra<-genotype.summary(as.matrix(genotypes.extra))
  colnames(aml.extra)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")
  aml.extra
  
  target<-"Control"
  use.samples<-the.samples[pheno[,"SampleProject"]=="0"]
  length(use.samples)
  genotypes.extra<-t(as.matrix(extra[use.samples]))
  dim(genotypes.extra)
  control.extra<-genotype.summary(as.matrix(genotypes.extra))
  colnames(control.extra)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")    
  control.extra
  summary.geno.extra.bit<-cbind(aml.extra,control.extra)
  
  
  
  target<-"AML.filt"
  use.samples<-the.samples[pheno[,"SampleProject"]=="1"]
  length(use.samples)
  genotypes.extra<-t(as.matrix(extra[use.samples]))
  dim(genotypes.extra)
  control.extra<-genotype.summary(as.matrix(genotypes.extra))
  colnames(control.extra)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")    
  control.extra
  summary.geno.extra.bit<-cbind(summary.geno.extra.bit,control.extra)
  
  
  
  target<-"Control.filt"
  use.samples<-the.samples[pheno[,"SampleProject"]=="0"]
  length(use.samples)
  genotypes.extra<-t(as.matrix(extra[use.samples]))
  dim(genotypes.extra)
  control.extra<-genotype.summary(as.matrix(genotypes.extra))
  colnames(control.extra)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")    
  control.extra
  summary.geno.extra.bit<-cbind(summary.geno.extra.bit,control.extra)
  
  summary.geno.extra<-rbind( summary.geno.extra,summary.geno.extra.bit)
  
  a.indel<-rbind(a.indel,extra)
  fil.genotypes<-rbind(fil.genotypes,extra.geno)
  key<-build.key(a.indel,core.ann)
  rownames(a.indel)<-key
  
  pass<-c(pass,TRUE)
  full.qual<-c(full.qual,TRUE)
  bad.effect<-c(bad.effect,TRUE)
  maf.filter<-c(maf.filter,TRUE)
  rare.in.group<-c(rare.in.group,TRUE)
  no.genotypes<-c(no.genotypes,FALSE)
  in.common.hit.gene<-c(in.common.hit.gene,FALSE)
  not.flat.genotype<-c(not.flat.genotype,TRUE)
  hw.controls.ok<-c(hw.controls.ok,TRUE)
  hw.controls.ok.filt<-c(hw.controls.ok,TRUE)
  on.x.y<-c(on.x.y,FALSE)
  unannotated.hits<-c(unannotated.hits,FALSE)
  are.repeats<-c(are.repeats,FALSE)
  bad.coding<-c( bad.coding,TRUE)
  bad.non.coding<-c(bad.non.coding,FALSE)
  functional.coding<-c(functional.coding,TRUE)
  
  ok.missing<-c(ok.missing,TRUE)
  are.in.repeats<-c(are.in.repeats,FALSE)
  
  ok.missing.filt<-c(ok.missing.filt,TRUE)
  hw.p.control.filt<-c(hw.p.control.filt,TRUE)
  rare.in.group.filt<-c(rare.in.group.filt,TRUE)
  no.genotypes.filt<-c(no.genotypes.filt,FALSE)
  rare.in.controls.filt<-c(rare.in.controls.filt,TRUE)
  is.unwound.geno<-c(is.unwound.geno,FALSE)
  ## tail(unannotated.hits)
  ##     tail(are.repeats)
  
  ## unannotated.hits<-unannotated.hits[1:342792]
  ## are.repeats<-are.repeats[1:342792]
  
  
  names(hw.controls.ok.filt)<-key
  names(ok.missing.filt.ok)<-key
  names(ok.missing.filt)<-key
  names(hw.p.control.filt)<-key
  names(rare.in.group.filt)<-key
  names(no.genotypes.filt)<-key
  names(rare.in.controls.filt)<-key
  
  names(pass)<-key
  names(full.qual)<-key
  names(bad.effect)<-key
  names(maf.filter)<-key
  names(rare.in.group)<-key
  names(no.genotypes)<-key
  names(in.common.hit.gene)<-key
  names(not.flat.genotype)<-key
  names(hw.controls.ok)<-key
  names(on.x.y)<-key
  names(unannotated.hits)<-key
  names(are.repeats)<-key
  names(bad.coding)<-key
  names(bad.non.coding)<-key
  names(summary.geno.extra)<-key
  names(ok.missing)<-key
  names(are.in.repeats)<-key
  names(functional.coding)<-key
  names(is.unwound.geno)<-key
  
  maf.lt.all<-rbind(maf.lt.all,maf.lt.all[2,]) #  dim(maf.lt.all)
  
  help<-cbind(full.qual,bad.coding,maf.filter,rare.in.group,no.genotypes,in.common.hit.gene ,hw.controls.ok,on.x.y,unannotated.hits,not.flat.genotype,are.repeats,are.in.repeats,ok.missing,ok.missing.filt,is.unwound.geno,(ok.missing.filt | is.unwound.geno) ,hw.p.control.filt,rare.in.group.filt,no.genotypes.filt,rare.in.controls.filt)
} # ftl3 additions

# length is number of variants
length(in.common.hit.gene)
length(bad.coding)
length(key)
length(functional.coding)
length(pass) 
length(full.qual) 
length(bad.effect) 
length(maf.filter) 
length(rare.in.group) 
length(no.genotypes) 
length(in.common.hit.gene) 
length(not.flat.genotype) 
length(hw.controls.ok) 
length(on.x.y) 
length(unannotated.hits) 
length(are.repeats) 
length(bad.coding) 
length(bad.non.coding) 
dim(summary.geno.extra) 
length(ok.missing) 
length(are.in.repeats)
length(ok.missing.filt)
length(hw.p.control.filt)
length(rare.in.group.filt)
length(no.genotypes.filt)
length(rare.in.controls.filt)
length(is.unwound.geno)
length(hw.controls.ok.filt)
length(hw.controls.ok)
dim( maf.lt.all)
tail(hw.controls.ok.filt)
#hw.controls.ok.filt<-hw.controls.ok.filt[1:342792]

################
#CG indels test output

## no.genotypes.table.AML<-(summary.geno.extra[,c("MAF.AML","MAF.AML.filt")]== 0) | summary.geno.extra[,c("MAF.AML","MAF.AML.filt")]=="NaN"  | (is.na( summary.geno.extra[,c("MAF.AML","MAF.AML.filt")])) # no genotypes in test classes for a mutataion after individaul quality filtering
## no.genotypes.table.AML[1:5,]
## no.genotypes.AML<-no.genotypes.table.AML[,"MAF.AML"] # was AND
## no.genotypes.AML.filt<-no.genotypes.table.AML[,"MAF.AML.filt"] #  combine.boolean(no.genotypes.table,c("MAF.AML.filt"),"AND") # was AND
## sum(no.genotypes.AML)
## sum(no.genotypes.filt.AML)

## pass<- full.qual & !no.genotypes.AML  & !in.common.hit.gene  & !are.repeats & !are.in.repeats & !is.unwound.geno & ok.missing & ok.missing.filt  & hw.controls.ok.filt & !no.genotypes.filt 

## sum(pass)
## sum(full.qual)
## sum(full.qual   & !in.common.hit.gene)
## sum(full.qual   & !in.common.hit.gene  & !are.repeats)
## sum(full.qual   & !in.common.hit.gene  & !are.repeats & !are.in.repeats)
## sum(full.qual   & !in.common.hit.gene  & !are.repeats & !are.in.repeats & !is.unwound.geno & hw.controls.ok.filt &  hw.controls.ok & !no.genotypes.filt )
## sum(full.qual   & !in.common.hit.gene  & !are.repeats & !are.in.repeats & !is.unwound.geno & hw.controls.ok.filt &  hw.controls.ok & !no.genotypes.filt & ok.missing & ok.missing.filt )
## core<-full.qual   & !in.common.hit.gene  & !are.repeats & !are.in.repeats & !is.unwound.geno & hw.controls.ok.filt &  hw.controls.ok & !no.genotypes.filt

## summary.geno.extra[core & !ok.missing,][1:5,]
## summary.geno.extra[no.genotypes.AML,][1:5,]
## summary.geno.extra[no.genotypes.AML,][1:5,]

#save.image("aml_HC_indel_working.RData")

length(in.common.hit.gene)
length(bad.coding)
length(key)
length(functional.coding)
length(pass) 
length(full.qual) 
length(bad.coding) 
length(maf.filter) 
length(rare.in.group) 
length(no.genotypes) 
length(in.common.hit.gene) 
length(not.flat.genotype) 
length(hw.controls.ok) 
length(on.x.y) 
length(unannotated.hits) 
length(are.repeats) 
length(bad.coding) 
length(bad.non.coding) 
dim(summary.geno.extra) 
length(ok.missing) 
length(are.in.repeats)
length(ok.missing.filt)
length(hw.p.control.filt)
length(rare.in.group.filt)
length(rare.in.group)
length(no.genotypes.filt)
length(rare.in.controls.filt)
length(is.unwound.geno)
length(hw.controls.ok.filt)
length(hw.controls.ok)
dim( maf.lt.all)
tail(hw.controls.ok.filt)

## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype & !are.repeats & !are.in.repeats &
## ( ok.missing.filt )   & hw.controls.ok.filt & hw.AOGC.ok & hw.AOGC.ok.filt & !no.genotypes.filt   & rare.in.controls.filt & rare.in.group & rare.in.AOGC & rare.in.AOGC.filt &
##  rare.in.ex.controls.filt  & hw.ex.controls.ok.filt & ( ok.missing.filt.ex ) # NMD

#  this is key - what ssends variants to be tested
pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene  & !on.x.y & !unannotated.hits & not.flat.genotype & !are.repeats & !are.in.repeats &
  ok.missing.filt & hw.controls.ok.filt & !no.genotypes.filt & rare.in.controls.filt

& rare.in.group

pass.ori<-pass

#pass<- full.qual & maf.filter & !in.common.hit.gene  & !on.x.y & !unannotated.hits & not.flat.genotype & !are.repeats & !are.in.repeats & ok.missing.filt & hw.controls.ok.filt & !no.genotypes.filt & rare.in.controls.filt & rare.in.group

help<-cbind(full.qual,bad.coding,maf.filter,rare.in.group,no.genotypes,in.common.hit.gene ,hw.controls.ok,on.x.y,unannotated.hits,not.flat.genotype,are.repeats,are.in.repeats,ok.missing,ok.missing.filt,is.unwound.geno,(ok.missing.filt | is.unwound.geno) ,hw.p.control.filt,rare.in.group.filt,no.genotypes.filt,rare.in.controls.filt )

dim(fil.genotypes)
sum(pass)
length(pass)
dim(snpinfo.ori)
snpinfo[1:5,]
a.indel[1:5,1:50]
pass[1:5]
dim(a.indel)

bad.genotypes<-c("chr6:35425714:35425714:-:C:indel","chr11:108126934:108126934:A:T:snp")
bad.genotypes %in% names(pass)
pass[ names(pass) %in% bad.genotypes]<-FALSE


################################# GEFOS FILTERING cause sending all
## snpinfo[grep("chr13:28626716:28626716:C:T:CREST",snpinfo[,"Name"]),]
## snpinfo.ori[grep("chr13:28626716:28626716:C:T:CREST",snpinfo.ori[,"Name"]),]

## table(a.indel[pass,"refGene::gene"])
## table(a.indel[pass,"Gene.name"])

#pass<-pass[the.snps] ### GEOFS 
genotypes<-a.indel[pass,the.samples] ## ordered correctly for phenotypes and have phenotypes
#genotypes<-fil.genotypes[pass,the.samples]

snp.names<-key[pass] ## GEFOS ony name with start

#### snpinfo now A different size than a.indel since added pathways!!!
snpinfo<-snpinfo.ori[snpinfo.ori[,"Name"] %in% snp.names,]
if( sum(!(snp.names %in% snpinfo.ori[,"Name"]))>0){print("WARINING snp.names not in snpinfo- unusual!")}
dim(snpinfo)
length(snp.names)
dim(genotypes)

# 414 639
print("start QC")
#RNPC3
genotypes[genotypes=="NA"]<-NA

genotypes[genotypes=="0/0"]<-0
genotypes[genotypes=="0/1"]<-1
genotypes[genotypes=="1/1"]<-2

########### prevent any averaging
dim(genotypes)
genotypes[is.na(genotypes)]<-0
dim(genotypes)
########### prevent any averaging


########################################## GEFOS MINOR ALLELE TRANSFORMATION
## flip.geno<-gsub("2","3",genotypes[to.flip,])
## #flip.geno[1:15,1:10]
## flip.geno<-gsub("0","2",flip.geno)
## flip.geno<-gsub("3","0",flip.geno)
## genotypes[to.flip,]<-flip.geno
##########################################################################

num.col<-dim(genotypes)[2]
num.row<-dim(genotypes)[1]
## genotypes[1:5,1:20]
genotypes<-as.numeric(as.matrix(genotypes))
dim(genotypes)<-c(num.row,num.col)
genotypes<-t(genotypes) # samples x SNPS

colnames(genotypes)<-snp.names
rownames(genotypes)<-gsub(".GT$","",the.samples)

#################################

dim(genotypes)
dim(pheno)
pheno[1:5,]
snpinfo[1:5,]
genotypes[1:5,1:5]

# cohort.seq.gene<-   skatCohort(Z=genotypes,formula, SNPInfo = snpinfo, data=pheno,family=binomial(),aggregateBy="gene",verbose=FALSE)
## meta.results.burden.gene<-burdenMeta(cohort.seq.gene,wts=1,mafRange = c(0,1),SNPInfo = snpinfo,aggregateBy = "gene")
## meta.results.skatO.gene<-skatOMeta(cohort.seq.gene,burden.wts =1,SNPInfo = snpinfo,aggregateBy="gene")

## the.order.gene<-     order(meta.results.burden.gene[,"p"])
## meta.results.burden.gene<-meta.results.burden.gene[the.order.gene,]
## meta.results.burden.gene[1:50,]

## the.order.gene<-     order(meta.results.skatO.gene[,"p"])
## meta.results.skatO.gene<-meta.results.skatO.gene[the.order.gene,]
## meta.results.skatO.gene[1:50,]

cohort.seq <- skatCohort(Z=genotypes,formula, SNPInfo = snpinfo, data=pheno,aggregateBy="cluster",family=binomial(),verbose=FALSE) ## genes and clusters

meta.results.burden<-burdenMeta(cohort.seq,wts=1,mafRange = c(0,1),SNPInfo = snpinfo,aggregateBy="cluster")

meta.results.skat<-skatMeta(cohort.seq,SNPInfo = snpinfo,aggregateBy="cluster")
meta.results.skatO<-skatOMeta(cohort.seq,burden.wts =1,SNPInfo = snpinfo,aggregateBy="cluster")

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

meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,]
meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]
##         meta.results.burden.gene[meta.results.burden.gene[,"gene"] %in% fanc.genes,]
## snpinfo.ori<-snpinfo 
## meta.results.burden.gene[meta.results.burden.gene[,"gene"] %in% other.clusters[,2],]
## meta.results.burden.gene[meta.results.burden.gene[,"gene"] %in% other.clusters[,3],]

## meta.results.skatO.gene[meta.results.skatO.gene[,"gene"] %in% clinical.genes,]
## meta.results.skatO.gene[meta.results.skatO.gene[,"gene"] %in% fanc.genes,]

## meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,]

setwd(analysis.dir)
getwd()
snap.file<-"coding.0.01.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN"
snap.file<-"coding.0.001.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN"


write.table(meta.results.burden[1:50,],file=paste("Burden","Top100",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,],file=paste("Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


write.table(meta.results.skat[1:50,],file=paste("Skat",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skatO[1:50,],file=paste("SkatO",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skatO[1:50,],file=paste("SkatO",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,],file=paste("SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,],file=paste("Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


colnames(a.indel)[c(1:8,13,16,28,7,30,34,37:42,43,14,32,33,1276,1310)]
#annotations<-a.indel[,c(1:8,13,16,28,7,30,34,37:42,43,14,32,33,1276,1310)]
annotations<-a.indel[,c(1:8,13,16,28,7,30,34,37:42,43,14,32,33)]

save(list=c("maf.col","alt.counts.thresh","formula","p","cohort.seq","meta.results.skat","meta.results.skatO","meta.results.burden","pheno","snpinfo","genotypes","pass","help","high.missing","annotations","key","summary.geno.extra","clusters.wanted"),file=paste(snap.file,"RData",sep="."))
getwd()
#save.image(file=paste("AML_TCGA_image_paper",snap.file,"RData",sep="."))
save(list=c("clusters.wanted"),file="clusters.wanted.RData")
getwd()

## setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-11-04_AML_TCGA_Replication/Analysis")
## load("AML_TCGA_image_paper.coding.0.01.all.geno.all.filters_no.imput_paper_TCGA_REP.RData")

## setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis")
## load("coding.0.001.all.geno.all.filters_no.imput_paper.RData")
#load("AML_HC_image_paper.coding.0.01.all.geno.all.filters_no.imput_paper.RData")

###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
##################### RELOAD########################

## analysis.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis"
## setwd(analysis.dir)
## getwd()

## snap.file<-"coding.0.01.all.geno.all.filters_no.imput"
## snap.file<-"coding.0.01.all.geno.all.filters.NEW"
## snap.file<-"coding.0.001.all.geno.all.filters.NEW"

load(paste(snap.file,"RData",sep="."))
meta.results.burden[1:50,]
meta.results.skat[1:50,]
meta.results.skatO[1:50,]

covars<-"1"
target.pheno.col<-"SampleProject"
formula<-paste(target.pheno.col,"~",paste(covars,collapse="+"),sep="")
print(formula)
formula<-formula(formula)
pheno[1:5,]

###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#to.unwind<-c(meta.results.burden[1:50,"gene"],meta.results.skatO[1:50,"gene"])
to.unwind<-c("FANC_complex.all") # to.unwind<-meta.results.burden[8,"gene"] - give it cluster/gene name and also name in snpinfo file
## to.unwind<-c("NPM1")
## to.unwind<-c("FLT3")
## to.unwind<-c("Clinical")
to.unwind<-c(clusters.wanted[!(clusters.wanted %in% c("Ubin.proteo","lipid_raft","caveolae","Checkpoint_extendedx1","Checkpoint_extendedx2"))])

to.unwind
to.unwind.name<-to.unwind
to.unwind.name<-"ALL_clusters_ALL_mutations"
# to.unwind.name<-"ALL_clusters"
# to.unwind.name<-"ALL_significant"
# to.unwind.name<-"ALL_significant"

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

meta.results.skatO[1:50,]
meta.results.burden[1:50,]
the.genes.burden<-meta.results.burden[meta.results.burden[,"gene"] %in% the.genes,]

the.genes.burden
write.table(the.genes.burden,file=paste(to.unwind.name,"conponents:","Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

the.genes.burden<-meta.results.skatO[meta.results.skatO[,"gene"] %in% the.genes,]
the.genes.burden
write.table(the.genes.burden,file=paste(paste(to.unwind.name,collapse="."),"conponents:","SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
#subset<-(rownames(genotypes) %in% loci) # indicated in genotypes

#sum(subset)
#length(subset)
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
dim(pheno)
snpinfo.ex[1:5,]
# summary.geno.extra[loci,]
#high.missing[loci,]
sum(are.in.repeats[loci])

# qual[loci,]
## snpinfo[1:5,]
## qual[1:5,c("FILTER_PASS", "FILTER_100" )]

cohort.seq.ex <- skatCohort(genotypes.ex, formula, SNPInfo = snpinfo.ex, data=pheno,aggregateBy = "Name",verbose=FALSE)
## meta.results.skat.ex<-skatMeta(cohort.seq,SNPInfo = snpinfo)
meta.results.burden.ex<-burdenMeta(cohort.seq.ex,wts=1,mafRange = c(0,1),SNPInfo = snpinfo.ex,aggregateBy = "Name")
meta.results.burden.ex
pheno[1:5,]

cohort.seq.test <- skatCohort(genotypes.ex, formula, SNPInfo = snpinfo.ex, data=pheno,aggregateBy = "cluster",verbose=FALSE)

meta.results.burden.test<-burdenMeta(cohort.seq.test,wts=1,mafRange = c(0,1),SNPInfo = snpinfo.ex,aggregateBy = "cluster")
#meta.results.burden.test

## meta.results.skat.ex<-skatMeta(cohort.seq,SNPInfo = snpinfo)
meta.results.skatO.test<-skatOMeta(cohort.seq.test,burden.wts =1,SNPInfo = snpinfo.ex,aggregateBy="cluster")
#meta.results.skatO.test

muts.in.cases<-apply(genotypes.ex[pheno[,"SampleProject"]==1,],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
muts.in.controls<-apply(genotypes.ex[pheno[,"SampleProject"]==0,],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})

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
    the.gt<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GT",sep=".")
    #the.gq<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GQ",sep=".")
    #the.dp<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"DP",sep=".")
    
    quality.cases[check]<-paste(a.indel.sub[posn,the.gt],collapse=",")
    
    #a.indel[posn,the.gq]
    a.indel[posn,the.gt]
    #a.indel[posn,the.dp]
  }
  
  if(muts.in.controls[check]!=""){
    the.gt<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"GT",sep=".")
    #the.gq<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"GQ",sep=".")
    #the.gq<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"DP",sep=".")
    quality.controls[check]<-paste(a.indel.sub[posn,the.gt],collapse=",")
    
    #a.indel[posn,the.gq]
    a.indel[posn,the.gt]
    #a.indel[posn,the.dp]
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



out<-cbind(meta.results.burden.ex,annotations[figure,],summary.geno.extra[figure,c("GENO.AML","GENO.Control","GENO.AML.filt","GENO.Control.filt")],help[figure,],muts.in.cases,quality.cases,muts.in.controls,quality.controls)

#out<-cbind(meta.results.burden.ex,annotations[figure,],muts.in.cases,muts.in.controls)
dim(out)


## table(out[,"refGene::location"])
## table(out[,"Consequence.Embl"])

################ PROBLEM IN NAME - TOO LONG FOR WINDOWS!!!

paste(paste(to.unwind,collapse="."))
paste(to.unwind.name,collapse=".")
paste(paste(to.unwind.name,collapse="."),"GENOTYPE.conponents:","SkatO","clusters",snap.file,"txt",sep=".")
#write.table(out,file=paste(paste(to.unwind.name,collapse="."),"GENOTYPE.conponents:","SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(out,file=to.unwind.name,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

getwd()

#######################################################################
