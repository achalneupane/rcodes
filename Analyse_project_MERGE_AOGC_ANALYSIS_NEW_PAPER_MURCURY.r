
######## ONLY NEED TO CHOOSE A DIRECTORY AND EXTENSIONS - used tab delimited files 
###############################################
###### THIS IS the super annotion run USE ALL THE FILTER AND NOVEL DATABASES
##Build AOGC

##################### if have a genotype component
genotype.file.location<-"/media/UQCCG-Analysis/AOGC_exome_chip/working_genotypes"
genotype.file.prefix<-"recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL"


related.to.remove<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.recode.txt",header=F,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
#related.to.remove<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.txt",header=F,fill=TRUE,sep=",",stringsAsFactors=FALSE)
related.to.remove[1:5,]
dim(related.to.remove)
#sum(related.to.remove[,1] %in% fam[,1])

#########################################################
#/media/old-scratch/media/scratch2/AOGC-NGS/Analysis/AOGC-Genotyping.output.chr2.AOGC_ALL.analysis-maf-filtered.txt
annotate.dir<-"/media/old-scratch/media/scratch2/AOGC-NGS/Annotate" # dir(annotate.dir)
analysis.dir<-"/media/old-scratch/media/scratch2/AOGC-NGS/Analysis" # dir(analysis.dir)
project.extension<-".analysis-maf-filtered.txt" ## justX1 the exterion not fam.extension!
project.name<-"AOGC-Genotyping.output"
fam<-c("AOGC_ALL") #  ALL or  c() ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_SAMPLES_PHENOTYPES_Nov.1.2013_RESIDUALS.txt"

the.sample.sheet<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
## ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
contaminated.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/contaminated_AOGC_SEQ_samples.txt"
contaminated<-read.table(contaminated.file,header=F,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)
remove.cols<-c()

regions.file<-"/media/scratch2/AOGC-NGS/GFOS/gefos.seq/METHODS/0613-skatmeta-gefos/static/Homo_sapiens.GRCh37.70.protein_coding.genespace_boundaries.5k.split100k.txt"
core.ann<-c("chr","start","end","REF","ALT","TYPE") # out put to annanlsys programs and need foe colun labels
dont.build.summary<-TRUE ##
GATK.SB<-TRUE
maf.threshold.filter.to.use<-c(0.05)

a.label<-"exome_chip.PASS.GENE.regions"
dont.build.summary<-TRUE


#########################################################
#/media/UQCCG/Sequencing/Projects/BONE-GENOME/Analysis/mergeCHR.chrALL.ACCupdate.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
annotate.dir<-"/media/UQCCG/Sequencing/Projects/BONE-GENOME/Annotate" # dir(annotate.dir)
analysis.dir<-"/media/UQCCG/Sequencing/Projects/BONE-GENOME/Analysis" # dir(analysis.dir)
project.extension<-"_analysis-maf-filtered.txt" ## justX1 the exterion not fam.extension!
project.name<-"mergeCHR."
fam<-c(".ACCupdate.ALL.ALL_GENOTYPES") #  ALL or  c() ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_SAMPLES_PHENOTYPES_Nov.1.2013_RESIDUALS.txt"

the.sample.sheet<-"/media/UQCCG/Sequencing/Projects/BONE-GENOME/Analysis/sample_sheet.csv"
## ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
contaminated.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/contaminated_AOGC_SEQ_samples.txt"
contaminated<-read.table(contaminated.file,header=F,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)
remove.cols<-c()

regions.file<-"/media/scratch2/AOGC-NGS/GFOS/gefos.seq/METHODS/0613-skatmeta-gefos/static/Homo_sapiens.GRCh37.70.protein_coding.genespace_boundaries.5k.split100k.txt"
core.ann<-c("chr","start","end","REF","ALT","TYPE") # out put to annanlsys programs and need foe colun labels
dont.build.summary<-TRUE ##
GATK.SB<-TRUE
maf.threshold.filter.to.use<-c(0.05)

a.label<-"exome_chip.PASS.GENE.regions"
dont.build.summary<-TRUE


###############
library(skatMeta)  ## ridge regression
#library(SKAT) ## skat method
library(GenomicFeatures)
library(HardyWeinberg)
library(Biostrings)

options(width=250,max.print=5000)

code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")
source("hwe.r")
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
##################################################### DEFINE A GENE LIST  #####################################################
##################################################### DEFINE A GENE LIST  #####################################################




##################################################### PILOT - GENE LIST #####################################################
gene.list.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP_lists_used_in_exomeChip/AOGC_Gene_Hits_0.05_0.2-sent.txt"
gene.list<-read.delim(gene.list.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
gene.list[1:5,]
colnames(gene.list)[colnames(gene.list)=="chr"]<-"CHR"
colnames(gene.list)[colnames(gene.list)=="start"]<-"START"
colnames(gene.list)[colnames(gene.list)=="end"]<-"END"
colnames(gene.list)[colnames(gene.list)=="value"]<-"GENE_NAME"
gene.list[1:5,]

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
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

############################################# POPULATION MAF FILTER - PART A
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER

maf.threshold<-0.0  #  MAF threshold for annovar calling zero useful to get back all results !!do not modify!!
maf.threshold.filter.to.use<-c(0.001,0.01,0.05)
maf.threshold.filter.to.use<-sort(as.numeric(maf.threshold.filter.to.use))

filter.cols.novel.use<-c("PopFreqMax","NHBLI_6500_ANNOVAR_ALL","NHBLI_6500_ALL","NHLBI_5400_ALL","NHLBI_5400_EUR","NHLBI_5400_AFR","1000genome","1000genome_asian","1000genome_mine","snp141","snp141_clinical","snp137","CG69","EUR_ASN_AFR_INDEL","AOGC-NGS_ALL","AOGC-NGS_ALL_OLD","Chinese") ##
filter.cols.maf.use<-c("PopFreqMax","NHBLI_6500_ANNOVAR_ALL","NHBLI_6500_ALL","NHBLI_6500_EA","NHBLI_6500_AA","NHLBI_5400_ALL","1000genome","snp141","snp137","snp135")
 maf.threshold.filter<-maf.threshold.filter.to.use

############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER



######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
############################################SET UP FNCTIONAL  FILTERS #######################################################

####################### MUTATION TYPE DEINITIONS 

possible.mutations<-c("frameshift substitution","nonframeshift substitution","downstream","frameshift deletion","frameshift insertion","intergenic","intronic","ncRNA_exonic","ncRNA_intronic","ncRNA_splicing","ncRNA_UTR3","ncRNA_UTR5","ncRNA_UTR5;ncRNA_UTR3","nonframeshift deletion","nonframeshift insertion","nonsynonymous SNV","splicing","stopgain SNV","stoploss SNV","synonymous SNV","unknown","upstream","upstream;downstream","UTR3","UTR5","UTR5;UTR3")

interesting.coding.mutations<-c("frameshift substitution","nonframeshift substitution","nonframeshift deletion","nonframeshift insertion","frameshift deletion","frameshift insertion","nonsynonymous SNV","stopgain SNV","stoploss SNV","splicing")

interesting.mutations.use<-c("frameshift substitution","nonframeshift substitution","nonframeshift deletion","nonframeshift insertion","frameshift deletion","frameshift insertion","nonsynonymous SNV","stopgain SNV","stoploss SNV","splicing","ncRNA_exonic")

wanted.noncoding.subtypes<-c("miRNA","lincRNA") # filter by interesting to prefiler and vep.noncoding so dones get ncRNA intronic ::use gerp.score.threshold.low only these subtypes

interesting.to.prefilter<-c("UTR3","UTR5","UTR5;UTR3","snoRNA","snRNA","antisense","sense_intronic","ncRNA_exonic","ncRNA_splicing") #use gerp.score.threshold

extra.vep.annotations<-c("Uploaded_variation","Gene","Feature","Protein_position","Amino_acids")

vep.types<-c("not_assigned","stop_gained","stop_lost","missense_variant","splice_acceptor_variant","splice_donor_variant","splice_region_variant","initiator_codon_variant","stop_retained_variant","incomplete_terminal_codon_variant","frameshift_variant","inframe_deletion","inframe_insertion","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_exon_variant","NC_stop_gained","NC_stop_lost","NC_splice_acceptor_variant","NC_splice_donor_variant","NC_splice_region_variant","NC_initiator_codon_variant","NC_stop_retained_variant","NC_non_coding_exon_variant","NC_incomplete_terminal_codon_variant","NC_3_prime_UTR_variant","mature_miRNA_variant","NC_5_prime_UTR_variant","TF_binding_site_variant","TFBS_ablation","TFBS_amplification","regulatory_region_variant","intron_variant","NC_intron_variant","synonymous_variant","coding_sequence_variant","NC_synonymous_variant","upstream_gene_variant","downstream_gene_variant","intergenic_variant","NC_intergenic_variant","NMD_transcript_variant","nc_transcript_variant","NC_nc_transcript_variant","feature_truncation","feature_elongation")

## vep.types<-c( "not_assigned","stop_gained","stop_lost","stop_lost,NMD_transcript_variant","stop_gained,splice_region_variant,NMD_transcript_variant","initiator_codon_variant,splice_region_variant","splice_region_variant,3_prime_UTR_variant","stop_gained,NMD_transcript_variant","missense_variant,splice_region_variant","missense_variant","splice_acceptor_variant","splice_acceptor_variant,nc_transcript_variant","splice_region_variant,3_prime_UTR_variant,NMD_transcript_variant","splice_donor_variant,nc_transcript_variant","splice_region_variant,intron_variant,NMD_transcript_variant","splice_donor_variant","splice_region_variant","splice_region_variant,5_prime_UTR_variant","splice_region_variant,synonymous_variant","splice_region_variant,intron_variant,nc_transcript_variant","splice_region_variant,non_coding_exon_variant,nc_transcript_variant","missense_variant,NMD_transcript_variant","splice_region_variant,intron_variant","NMD_transcript_variant","intron_variant,NMD_transcript_variant","mature_miRNA_variant","5_prime_UTR_variant","5_prime_UTR_variant,NMD_transcript_variant","non_coding_exon_variant,nc_transcript_variant","3_prime_UTR_variant,NMD_transcript_variant","non_coding_exon_variant","TF_binding_site_variant","intron_variant,nc_transcript_variant","synonymous_variant,NMD_transcript_variant","3_prime_UTR_variant","regulatory_region_variant","upstream_gene_variant","downstream_gene_variant","intergenic_variant","intron_variant","synonymous_variant")


vep.coding<-c("not_assigned","stop_gained","stop_lost","missense_variant","splice_acceptor_variant","splice_donor_variant","splice_region_variant","initiator_codon_variant","stop_retained_variant","incomplete_terminal_codon_variant","frameshift_variant","inframe_deletion","inframe_insertion")

vep.noncoding<-c("5_prime_UTR_variant","3_prime_UTR_variant","non_coding_exon_variant","NC_stop_gained","NC_stop_lost","NC_splice_acceptor_variant","NC_splice_donor_variant","NC_splice_region_variant","NC_initiator_codon_variant","NC_stop_retained_variant","NC_non_coding_exon_variant","NC_incomplete_terminal_codon_variant","NC_3_prime_UTR_variant","mature_miRNA_variant","NC_5_prime_UTR_variant","TF_binding_site_variant","TFBS_ablation","TFBS_amplification","regulatory_region_variant")              

vep.unwanted<-c("intron_variant","NC_intron_variant","synonymous_variant","coding_sequence_variant","NC_synonymous_variant","upstream_gene_variant","downstream_gene_variant","intergenic_variant","NC_intergenic_variant","NMD_transcript_variant","nc_transcript_variant","NC_nc_transcript_variant","feature_truncation","feature_elongation")



missense.variant<-c("nonsynonymous SNV","missense_variant")

hwe.control.threshold<-1e-8
gerp.score.threshold.high<-2.5 # gerp score >= will be included
gerp.score.threshold.low<-2.0 # gerp score >= will be included
gerp.score.threshold.unknown<-0


#generic.filter.DB

 interesting.mutations<-interesting.mutations.use


if(GATK.SB){
global.quality.labs<-c("QUAL","QD","HRun","SB","FILTER","FILTER","PolyPhen.scores","PolyPhen.scores","SIFT.scores","mut.taster::score","phylo::score","PolyPhen.desc","SIFT.desc","GERP::score","GERP::score","GERP::score","MAF.ALL","MAF.HIGH","MAF.LOW","TYPE") ### THESE ARE THE COLUMN LABELS IN THE DATA these become the "good.qual" filter 
global.quality.names<-c("QUAL","QD","HRun","SB","FILTER_PASS","FILTER_100","PolyPhen.low","PolyPhen.high","SIFT.high","mut.taster.high","phylo.high","PolyPhen.bad","SIFT.bad","GERP.high","GERP.low","GERP.unknown","MAF.ALL","MAF.HIGH","MAF.LOW","flat") ### THESE ARE THE COLUMN LABELS IN the quality.filter TABLE
#global.quality.cut<-c(50,0.5,5,1,"PASS","TruthSensitivityTranche99.90to100.00",0.1,0.4,0.4,0.4,0.4,"damaging","deleterious",2,2,0.25,0.25,0.25)
global.quality.cut<-c(50,0.5,5,1,"PASS","TruthSensitivityTranche99.90to100.00",0.1,0.4,0.4,0.4,0.4,"damaging","deleterious",gerp.score.threshold.high,gerp.score.threshold.low,gerp.score.threshold.unknown,0.25,0.25,0.25,"flat")
global.quality.type<-c("numeric","numeric","numeric","numeric","factor","factor","numeric","numeric","numeric","numeric","numeric","factor","factor","numeric","numeric","numeric","numeric","numeric","numeric","factor")
global.quality.dirn<-c("greater","greater","less","less","exact","exact","greater","greater","greater","greater","greater","ends_with","exact","greater","greater","exact","greater","greater","greater","ends_with")
}else{
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


seq.type.file<-"/media/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_BAM_Fri_Apr_26_2013.txt"
seq.type<-read.delim(seq.type.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

seq.type[1:5,]
seq.type<-seq.type[seq.type[,"Project"]=="AOGC-NGS",]
nim.samples<-seq.type[seq.type[,"Capture.Method"]=="TruD:NimX","Sample"]
ill.samples<-seq.type[seq.type[,"Capture.Method"]=="TruD:TruX","Sample"]

nim.samples<-paste(nim.samples,"GT",sep=".")
ill.samples<-paste(ill.samples,"GT",sep=".")


length(nim.samples)
length(ill.samples)

############################################ Nimblegen and illuma capture loci ####################################################


## library("BSgenome.Hsapiens.UCSC.hg19")
## the.chroms<-seqlengths(Hsapiens)


## load("/media/UQCCG/Sequencing/Data/Genomes/hg19/Human_Exome_Targets_illumina_v2_hg19_targets.RData")
## ill.gr<-data.gr
## load("/media/UQCCG/Sequencing/Data/Genomes/hg19/Human_Exome_Targets_Nimble_v2_hg19_targets.RData")
## nim.gr<-data.gr

## genome(ill.gr)<-"hg19"
## genome(nim.gr)<-"hg19"

## ill.gr<-ill.gr+200
## nim.gr<-nim.gr+200

## overlaps<-overlapsAny(ill.gr,nim.gr)
## possible.loci.ori<-ill.gr[overlaps]

######################################################################################################################################
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

the.sample.sheet

sample.sheet.full<-read.delim(the.sample.sheet,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
sample.sheet.full[1:5,]
colnames(sample.sheet.full)
dim(sample.sheet.full)


##### fix 0 and 9 for missing to NA

## pheno.types<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_RAD","BMD_EFF_STD_LS","BMD_EFF_STD_FN","EVER_FX_50_EXCL_TRIVIAL")
## names(pheno.types)<-c("HIP","RAD","LS","FN","FX")

sample.sheet.full[1:5,pheno.types]
table(sample.sheet.full$BMD_AFFSTAT)

pheno.types<-c("SampleProject") ## vales is column header
names(pheno.types)<-c("HIP") ### name is output columns

## case.control.classes<-c(0,1,0)
## names(case.control.classes)<-c("Control","NMD","AOGC")
## case.control.classes



case.control<-c("AffectionStatus")
case.control.classes<-c(1,2)
names(case.control.classes)<-c("1","2")
case.control.classes
# ib<-1
for(ib in 1:length(case.control)){
  if(!(case.control[ib] %in% colnames(sample.sheet.full))){next}
  sample.sheet.full[(  !(sample.sheet.full[,case.control[ib]] %in% names(case.control.classes))  |  is.na(sample.sheet.full[,case.control[ib]]) | sample.sheet.full[,case.control[ib]]==0 | sample.sheet.full[,case.control[ib]]==9)  ,case.control[ib]]<-NA
}

#tapply(sample.sheet.full[,"SampleProject"],sample.sheet.full[,"SampleProject"],length)

control.samples<-{}
all.samples<-sample.sheet.full[,"PATIENT"]

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

ifam<-1
for(ifam in 1:length(fam)){

the.extension<-paste(fam[ifam],project.extension,"$",sep="")
project.files<-files[grepl(the.extension ,files)]
print(sort(paste("Doing: ",project.files,sep=""))) # project.files<-project.files[1:22]

indels<-{}
the.col<-{}
project.files
# ichr<-13
project.files
if(length(project.files)>1){
 project.files<-project.files[!grepl("chrALL",project.files)]
 

  if(length(project.files)!=24){
    print("########################################### WARNING #################################")
    print("less that 24 chromosomes detected")
    print(fam[ifam])
    print("########################################### WARNING #################################") 
  }
  project.files


ichr<-1

for(ichr in 1:length(project.files)){

setwd(analysis.dir)
################## fast read ###########
column.labels<-read.delim(project.files[ichr],header=F,nrows=1,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="") # quote="\""
num.vars<-dim(column.labels)[2]
a.indel<-scan(project.files[ichr],what=character(num.vars),skip=1,sep="\t",fill=TRUE,na.strings="",quote="")
num.lines<-length(a.indel)/(num.vars)
dim(a.indel)<-c(num.vars,num.lines)
a.indel<-t(a.indel)
colnames(a.indel)<-column.labels
########################################
## chk<-seq(from=1,to=279233293,by=num.vars)
## table(a.indel[chk])


table(a.indel[,"chr"])
a.indel[1:5,1:10]

the.chr<-a.indel[1,"chr"]
print(paste("Doing Chromosome ",the.chr))


if(!grepl("^chr",the.chr)){
a.indel[,"chr"]<-paste("chr",a.indel[,"chr"],sep="")
}
key<-build.key(a.indel,core.ann)
rownames(a.indel)<-key

####### REMOVE BAD SAMPLES
bad.samples<-contaminated[,1]
bad.samples.labels<-expand.labels.to.samples(bad.samples,c("GT","AD","DP"),paste.after=TRUE)

a.indel<-a.indel[,colnames(a.indel)[!(colnames(a.indel) %in% bad.samples.labels)]]
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#################### Read plink file
## if(grepl("^chr",the.chr)){plink.chr<-gsub("^chr","",the.chr)}else{
##   plink.chr<-the.chr}




############################################# POPULATION MAF FILTER - PART B
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER


maf.threshold.filter<-maf.threshold.filter.to.use
filter.cols.maf<-filter.cols.maf.use
filter.cols.novel<-filter.cols.novel.use

#############################get the filters so have minor allele frequencies:

the.combined.DBs<-gsub("::maf$","",colnames(a.indel)[grep("::maf$",colnames(a.indel))]) ## get maf databases
the.combined.DBs<-the.combined.DBs[!(the.combined.DBs %in% c("ID","snp132"))] ## remove artifacts
names.the.combined.DBs<-the.combined.DBs
#the.combined.DBs<-paste(the.combined.DBs,"::maf",sep="")
names(the.combined.DBs)<-names.the.combined.DBs

 present.cols<-filter.cols.maf %in% names(the.combined.DBs)
if(sum(!present.cols)>0){  ### asked to filter on columns that don't exist
  print("WARNING asked to filter using MAF with databaes that have not been annotated")
  print(paste("missing MAF Db:",toString(filter.cols.maf[!present.cols]),"-> ignored",sep=" "))
  filter.cols.maf<-filter.cols.maf[present.cols]
}

 present.cols<-filter.cols.novel %in% names(the.combined.DBs)
if(sum(!present.cols)>0){  ### asked to filter on columns that don't exist
  print("WARNING asked to filter using NOVEL with databaes that have not been annotated")
  print(paste("missing MAF Db:",toString(filter.cols.novel[!present.cols]),"-> ignored",sep=" "))
  filter.cols.novel<-filter.cols.novel[present.cols]
}

print(filter.cols.maf)
print(filter.cols.novel)

maf.target<-the.combined.DBs[names(the.combined.DBs) %in% filter.cols.maf] # the.combined.DBs
maf.target

for(i in 1:length(maf.threshold.filter)){
a.filter<-paste(filter.cols.maf,"::maf-filter::",maf.threshold.filter[i],sep="")
print(a.filter)
 assign(paste("filter.cols.maf",maf.threshold.filter[i],sep="."),value=a.filter)
}

maf.threshold.filter.all<-c(0.0,maf.threshold.filter)
assign(paste("filter.cols.maf",maf.threshold.filter.all[1],sep="."),value=paste(filter.cols.novel,"::maf",sep=""))
print(eval(as.name(paste("filter.cols.maf",maf.threshold.filter.all[1],sep=".")) ) )

#  a.indel[1:5,paste(the.combined.DBs,"::maf",sep="")]


################################################################################
## THE FILTER TABLE IS MADE FROM : the.combined.DBs<-c(filter.DB,generic.filter.DB,function.filter.DB)  only thes DBs can be used for MAF filtereing
target.table<-a.indel[,paste(the.combined.DBs,"::maf",sep="")]
#target.table[1:5,]
extra<-matrix(data=NA,nrow=dim(target.table)[1],ncol=length(maf.target)*length(maf.threshold.filter))

colnames(extra)<-expand.labels.to.samples.complex(maf.target,maf.threshold.filter,paste.after=TRUE,seperator="::maf-filter::")
#extra[1:5,]

target.table<-cbind(target.table,extra)
rm(extra)
####### generate maf.target appreded to the allel frequency
maf.threshold.filter  # filters to apply
maf.target  #<-the.combined.DBs[names(the.combined.DBs) %in% all.filter.cols.maf] # the.combined.DBs

k<-2
i<-2


for(k in 1:length(maf.threshold.filter)){ #maf.threshold.filter.all contains novel "0"
for(i in 1:length(maf.target)){
  an.inversion<-(target.table[, paste(maf.target[i],"maf",sep="::") ] > (1-maf.threshold.filter[k]) )  # 0.999 on a 0.01 filter is also a hit
  an.inversion[is.na(an.inversion)]<-FALSE
#  a.test<-(target.table[,paste(maf.target[i],"maf",sep="::")] < maf.threshold.filter[k]) # normal test for minot allele
  target.table[,paste(maf.target[i],"maf-filter",maf.threshold.filter[k],sep="::")]<-(target.table[,paste(maf.target[i],"maf",sep="::")] < maf.threshold.filter[k]) | an.inversion  ## incase ref and alt are inverted
  target.table[is.na(target.table[,paste(maf.target[i],"maf-filter",maf.threshold.filter[k],sep="::")]),paste(maf.target[i],"maf-filter",maf.threshold.filter[k],sep="::")]<-TRUE # if not found in database then assume has small allele frequency
#  ref.inversion[,k ]<-ref.inversion[,k] &   an.inversion
}
print(paste("Done ",maf.threshold.filter[k],sep=""))
    }

#target.table[1:5,]
###################

## colnames(target.table)

## the.DBs<-the.combined.DBs
## the.DBs
## names(the.DBs)
## for(i in 1:length(the.DBs)){
  
## target.string<-paste("^",the.DBs[i],sep="")
## if(grepl("++",target.string,fixed=TRUE)){target.string<-gsub("++","",target.string,fixed=TRUE)} ##gerp++ casles problems in grep (++)
## colnames(target.table)<-gsub(target.string,names(the.DBs)[i],colnames(target.table))

## ## colnames(target.table)<-gsub(paste("^",the.DBs[i],sep=""),names(the.DBs)[i],colnames(target.table))
## }
## colnames(target.table)

## target.table[1:2,]

## filter.table<-target.table
## ########
# filter.cols.maf.0

  ################################# Get allele frequency table
################# HOWEVER ALT ALLELES ALL IS USED AS MINOR ALLELE FREQUNCY  not alternative allele frequencies
maf.lt.all<-data.frame(key=key,stringsAsFactors=FALSE)

  
  imaf<-1
  for(imaf in 1:length(maf.threshold.filter.all)){
  a.filter.cols.maf<-eval(as.name( paste("filter.cols.maf",maf.threshold.filter.all[imaf],sep=".")  ))
  maf.lt<-rep(TRUE,times=dim(target.table)[1])

  imaff<-1
  for(imaff in 1:length(a.filter.cols.maf)){
    if(maf.threshold.filter.all[imaf]==0){ # different case for novel NOT found (rather than less than)
  #     maf.lt<-maf.lt & !target.table[,a.filter.cols.maf[imaff]]  ### this is correct as using the variable string name derived from filter.cols.novel
        maf.lt<-maf.lt & ( is.na(target.table[,a.filter.cols.maf[imaff]])  | target.table[,a.filter.cols.maf[imaff]]=="NA" )
     }else{
    maf.lt<-maf.lt & as.logical(target.table[,a.filter.cols.maf[imaff]])
  }
  }
# filtered<-maf.lt & wanted.muts.fil
 maf.lt.all<-cbind(maf.lt.all,maf.lt)
}
  
if(dim(maf.lt.all)[2]>1){ colnames(maf.lt.all)<-c("key",paste("MAF.lt:",maf.threshold.filter.all,sep=""))}
maf.lt.all<- maf.lt.all[,colnames(maf.lt.all)!="key"]
maf.lt.all[1:4,]

colnames(a.indel)[colnames(a.indel) == "MAF.lt:0.5"]<-"MAF.lt:0.05"
posns<-match(colnames(maf.lt.all),colnames(a.indel))
missing<-is.na(posns)

## i<-1
## if(i in 1:sum(!missing)){
  
##  a.indel[, posns[!missing][i] ]  <- maf.lt.all[,!missing][i]
## }

dim(a.indel)
maf.lt.all[1:4,]
rm(target.table)



############ FIX ALLELE FREQUENCY PROBLEMS



## if(plink.chr=="X"){plink.chr<-23}
## if(plink.chr=="Y"){plink.chr<-24}
## if(plink.chr=="XY"){plink.chr<-25}
## if(plink.chr=="M"){plink.chr<-26}
## plink.chr
## plink.file<-paste(genotype.file.prefix,"_chr",plink.chr,sep="")
## plink.file
## g.indel<-read.plink(paste(genotype.file.location,plink.file,sep="/"))

## ## test<-read.table("plink.frq",header=T)
## ## ori<-read.table("ori.frq",header=T)

## dim(a.indel)
## dim(g.indel)
## target<-"37442658"
## posn1<-match(target,a.indel[,"end"])
## posn2<-match(target,g.indel[,"start"])
## posn1
## posn2

## ####### REMOVE BAD SAMPLES

## a.indel[posn1,1:10]
## g.indel[posn2,1:10]
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################







#tapply(gene.list[,"CHR"],gene.list[,"CHR"],length)
all.possible.samples<-gsub(".GT$","",colnames(a.indel)[grep(".GT$",colnames(a.indel))],perl=TRUE)
length(all.possible.samples)
pheno.types



 a.indel[1,gsub(".GT",".AD",use.samples)]
fil.genotypes[1,use.samples]

to.fix<-gsub(".GT",".AD",use.samples)
for( io in 1:length(to.fix)){
 bad<-a.indel[,to.fix[io]]=="50,50" & a.indel[,use.samples[io]]=="0/0"
a.indel[bad,to.fix[io]]<-"100,0"
}



#################################### got some missing geen names still.
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
#a.indel[a.dash,][1:5,1:50]
a.indel[a.dash,"Gene.Names"]<-a.indel[a.dash,"Feature.Embl"]

a.dash<-a.indel[,"Gene.Names"]=="-"
sum(a.dash)
#a.indel[a.dash,][1:5,1:50]

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

grep("NOTCH1",names(all.genes))
common.hit.genes<-names(all.genes)[1:30]
common.hit.genes<-all.genes[all.genes>250]
common.hit.genes
## all.genes["SCN2A"]

###############################################

#########################################################################################################################

snpinfo.raw<-cbind(key,a.indel[,"Gene.Names"],a.indel[,"Gene.Names"])
snpinfo.raw[1:5,]
tail(snpinfo.raw)

colnames(snpinfo.raw)<-c("Name","gene","cluster")
dim(snpinfo.raw)
snpinfo.raw[1:5,]
dim(snpinfo.raw)
snpinfo<-snpinfo.raw


clusters<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP_lists_used_in_exomeChip/clusters.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
clusters

#clusters.wanted<-c("Clinical","FANC - ACID")
clusters.wanted<-colnames(clusters)
ic<-1
snpinfo[1:5,]

#cbind(unique(clusters[,22]),unique(clusters[,22]))
snpinfo[1:5,]



gene.aliases<-read.delim("/media/UQCCG/Software/annovar/humandb/Gene_symbol_aliases.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
gene.aliases[1:5,]
gene.aliases<-unwind(gene.aliases, "Aliases",delimit=", ")
gene.aliases[1:5,]
###############  all.genes


## gene.aliases[grepl("SLX1",gene.aliases[,2]),]
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
   if(sum( missing.name)>0){
     posns<-match(cluster.genes[missing.name],gene.aliases[, "Aliases"])
     missing<-is.na(posns)
     if(sum(missing)>0){
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



### HAVE
# snpinfo.raw (original from a.indel)
# snpinfo # with extra clusters
# snpinfo.raw a permanent copy of snpinfo

## "FANCM " "MHF1"   "MHF2"   "FAAP24"
clusters[,1]
chk<-apply(clusters,2,function(x){ length(x[x!=""])})
write.table(clusters,file="clusters_as.using.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

############################

library(doMC)
num.cores<-4
num.bits<-num.cores
registerDoMC(cores=num.cores)

the.samples<-colnames(a.indel)[grepl(".GT$", colnames(a.indel))]

while((dim(a.indel)[1] %% num.bits)< 2){num.bits<-num.bits+1} ### go don't get matrix issues
num.bits
#(dim(a.indel)[1] %% num.bits)
## fil.genotypes<-foreach(a.indel.bit=iter(a.indel,by='row',chunksize=as.integer(dim(a.indel)[1]/num.bits) ), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype(a.indel.bit,gsub(".GT$","",the.samples),prefix="",suffix="",20,0.02,0.98,0.20,0.80,7,2)

fil.genotypes<-foreach(a.indel.bit=iter(a.indel,by='row',chunksize=as.integer(dim(a.indel)[1]/num.bits) ), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype.old(a.indel.bit,gsub(".GT$","",the.samples),prefix="",suffix="",20,0.02,0.98,0.20,0.80,10,5)
# 20,0.02,0.98,0.20,0.8,10,5 # for cancers where het may be way out
#test<-filtered.genotype.old(a.indel[1:5,],gsub(".GT$","",the.samples),prefix="",suffix="",20,0.02,0.98,0.20,0.80,10,5)
# 0.02< het(0/1)<0.98
# ref(0/0)< 0.2
# (1/1) > 0.8
##  a.indel[1,gsub(".GT",".AD",use.samples)]
## fil.genotypes[1,use.samples]




a.indel[1,use.samples]

dim(fil.genotypes)
colnames(fil.genotypes)[1:5]
rownames(fil.genotypes)[1:5]

dim(fil.genotypes)
dim(a.indel)

tail(rownames(a.indel))

## a.indel<-a.indel[1:342791,]
## fil.genotypes<-fil.genotypes[1:342791,]
############################### do one phenotype #################### 


############################### do one phenotype ##################################################################
ipheno<-1
for(ipheno in 1:length(pheno.types)){


##################################################### set formula
print(paste("Doing phenotype",pheno.types[ipheno]))
target.pheno<-names(pheno.types)[ipheno]
target.pheno.col<-pheno.types[ipheno]

if(target.pheno.col %in% case.control){
  covars<-c("AGE_SCAN","PCA1","PCA2","PCA3","PCA4") #AGE_SCAN,PCA1,PCA2,PCA3,PCA4 #covars<-c("1")
}else{
covars<-c("1")
}


formula<-paste(target.pheno.col,"~",paste(covars,collapse="+"),sep="")
print(formula)
formula<-formula(formula)


###############################################################

############################ subset samples with phenotype and covars - assume traits and covars need same phenotypes
sample.sheet.full[1:2,]
colnames(sample.sheet.full)
sample.sheet.full[1:5,]


#if(sum(covars==1) & length(covars)==1){
  got.all.covars<-rep(TRUE,times=dim(sample.sheet.full)[1])
#}else{
#got.all.covars<-apply(sample.sheet.full[,covars],1,function(x) (sum(is.na(as.numeric(x))) ==0 ))
#}

if(!("PATIENT" %in% colnames(sample.sheet.full))){ ## add sample column if using PATIENT
  PATIENT<-sample.sheet.full[,"ParticipantCode"]
  sample.sheet.full<-cbind(PATIENT,sample.sheet.full)
 #  colnames(pheno)[colnames(pheno)=="PATIENT"]<- "SAMPLE"
}

subset.samples<- (sample.sheet.full[,"PATIENT"] %in% all.possible.samples ) & !is.na(sample.sheet.full[,target.pheno.col]) & sample.sheet.full[,target.pheno.col]!="NA"  & got.all.covars
                  
sum(subset.samples)
pheno<-sample.sheet.full[subset.samples ,] ## pheno only contains SAMPLES that have a phenotype
colnames(pheno)
if(!("SAMPLE" %in% colnames(pheno))){ ## add sample column if using PATIENT
  SAMPLE<-pheno[,"PATIENT"]
  pheno<-cbind(SAMPLE,pheno)
 #  colnames(pheno)[colnames(pheno)=="PATIENT"]<- "SAMPLE"
}
colnames(pheno)
print(dim(pheno))
print(paste("Number Samples:",dim(pheno)[1]))

dim(pheno)
pheno[1:5,1:5]



the.samples<-paste(pheno[,"SAMPLE"],"GT",sep=".")  ## samples same order as in pheno
print(paste("Number samples: ",length(the.samples),sep=""))

#seq.type[1:57,]
posns<-match(pheno[,"SAMPLE"],seq.type[,"Sample"])
missing<-is.na(posns)
sum(missing)
pheno[missing,1:5]

capture<-seq.type[posns,"Capture.Method"]


table(capture) ### all illume here
pheno<-cbind(pheno,capture)

## ###############GEFOS SNP type restrictions

## the.snps<-grepl("^snp",a.indel[,"TYPE"]) ## use SNPs only for GEFOS
## the.dups<-duplicated(a.indel[,"start"]) ## GEFOS annotate by just position so remove polymorphic locations (keep the first)
## sum(the.dups)
## a.indel[the.dups,][1:5,1:10]
## the.snps<-the.snps & !the.dups
## a.indel[1:5,c(core.ann,"MAF.ALL")]
## a.indel[the.snps,"MAF.ALL"][1:10]
## ##########################################################################
## ########################################## GEFOS MINOR ALLELE TRANSFORMATION
## to.flip<-as.numeric(a.indel[the.snps,"MAF.ALL"]) >0.5 ##GefOS want to use the minor allele
## sum(to.flip)
## a.indel[the.snps,][to.flip,c(core.ann,"MAF.ALL")][1:5,]
## ################

###### the.samples and pheno in same order but the.samples has .GT extension.
pheno[,pheno.types[ipheno]]

## hist(as.numeric(pheno[,pheno.types[ipheno]]))

## if(length(unique((pheno[,pheno.types[ipheno]])))>4){ 
##   low<-as.numeric(pheno[,pheno.types[ipheno]]) <=0 ## quantative
##   high<-as.numeric(pheno[,pheno.types[ipheno]]) >0
## }else{ ## case.control
##   low<-pheno[,pheno.types[ipheno]]==1
##   high<-pheno[,pheno.types[ipheno]]==2
## }
## sum(low)
## sum(high)

## illumina<-the.samples %in% ill.samples
## nimble<-the.samples %in% nim.samples

all<-rep(TRUE,length(the.samples))

dim(pheno)
length(all)


############## SET up groups to be analysed

pheno<-cbind(pheno,all)

## th project becomes the master copy of the targets below.
the.projects<-c("all")
names(the.projects)<-the.projects
colnames(pheno)

pheno[1:5,]
                       

summary.geno.extra<-{}
####################################################################################
#################################################################################### REGULAR
#### MAY NEED TO ADJUST way use samples in selected based on selection below.
targets<-the.projects  #c("NMD","ex.Control","AOGC")
targets
###### the.samples and pheno in same order but the.samples has .GT extension.
it<-1
for(it in 1:length(targets)){
#use.samples<-the.samples[pheno[,"SampleProject"]==targets[it]]
use.samples<-the.samples[pheno[,targets[it]]]
length(use.samples)
genotypes<-a.indel[,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),names(targets)[it],sep=".")
#summary.geno[1:5,]
if(is.null(dim(summary.geno.extra))){
  summary.geno.extra<-summary.geno
}else{
  summary.geno.extra<-cbind(summary.geno.extra,summary.geno)
}
print(paste("Done: ",targets[it],sep=""))
} ## loop over targets
#################################################################################### FILTERED

targets<-the.projects #c("NMD","ex.Control","AOGC")
targets
names(targets)<-paste(targets,".filt",sep="")
targets
it<-1
for(it in 1:length(targets)){
#use.samples<-the.samples[pheno[,"SampleProject"]==targets[it]]
use.samples<-the.samples[pheno[,targets[it]]]
length(use.samples)
genotypes<-fil.genotypes[,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),names(targets)[it],sep=".")
#summary.geno[1:5,]
if(is.null(dim(summary.geno.extra))){
  summary.geno.extra<-summary.geno
}else{
  summary.geno.extra<-cbind(summary.geno.extra,summary.geno)
}
print(paste("Done: ",targets[it],sep=""))
} ## loop over targets


## #########################################
chk<-grep("^GENO",colnames(summary.geno.extra))
summary.geno.extra[1:5,chk]
#

######################################


a.indel[1:5,1:10]
summary.geno.extra[1:5,]
colnames(summary.geno.extra)
#rownames(summary.geno.extra)<-key
#getHWE(obs_hets, obs_hom1, obs_hom2)
hw.target<-"all"

hw.p.control<-getHWE(summary.geno.extra[,paste("GENO.",hw.target,sep="")]) ## used 16 CPUs
hw.p.control.filt<-getHWE(summary.geno.extra[,paste("GENO.",hw.target,".filt",sep="")]) ## used 16 CPUs


##### class(summary.geno.extra)

length(hw.p.control)
names(hw.p.control)<-key
names(hw.p.control.filt)<-key

hw.p.control[1:5]
hw.p.control.filt[1:5]
######### testing

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

summary.geno.extra[1:5,]
rownames(summary.geno.extra)<-key

######################## inbreeding and no genotypes

group.maf.thresh<-0.20 # group.maf.thresh<-0.10  if more common that this in Group then discard: Discard Rare but common in this cohort
missing.threshold<-0.20 # missing.threshold<-0.50  60 % of genotypes missing
missing.threshold.nimblgen<-0.20
missing.threshold.illumina<-0.20

targets<-the.projects

rare.in.group.table<-( summary.geno.extra[,paste("MAF.",c(the.projects,paste(targets,".filt",sep="")),sep="") ]< group.maf.thresh) | ( summary.geno.extra[,paste("MAF.",c(the.projects,paste(targets,".filt",sep="")),sep="")] > (1-group.maf.thresh)) | (is.na( summary.geno.extra[,paste("MAF.",c(the.projects,paste(targets,".filt",sep="")),sep="")]))
rare.in.group.table[1:5,]
summary.geno.extra[1:5,]

## rare.in.group.table<-( summary.geno.extra[,paste("MAF.",c(the.projects,paste(targets,".filt",sep="")),sep="") ]< group.maf.thresh) | ( summary.geno.extra[,paste("MAF.",c(the.projects,paste(targets,".filt",sep="")),sep="")] > (1-group.maf.thresh)) | (is.na( summary.geno.extra[,paste("MAF.",c(the.projects,paste(targets,".filt",sep="")),sep="")]))
## rare.in.group.table[1:5,]
## summary.geno.extra[1:5,]

#rare.in.group.test<-colnames(rare.in.group)
rare.in.group.test<-paste("MAF.",c("all","all.filt"),sep="")
rare.in.group.test
rare.in.group<-combine.boolean(rare.in.group.table,rare.in.group.test,"OR")
sum(!rare.in.group)
rare.in.group[1:10]
## rare.in.group<-combine.boolean(rare.in.group.table,c("MAF.AML","MAF.Control"),"OR")
## rare.in.group.filt<-combine.boolean(rare.in.group.table,c("MAF.AML.filt","MAF.Control.filt"),"OR")
## sum(!rare.in.group)
## sum(!rare.in.group.filt)

no.genotypes.test<-paste("MAF.",c("all","all.filt"),sep="")
no.genotypes.test
no.genotypes<-(summary.geno.extra[,no.genotypes.test]== 0)  | (is.na( summary.geno.extra[,no.genotypes.test])) # no genotypes in test classes for a mutataion after individaul quality filtering
no.genotypes[1:5,]
no.genotypes<-combine.boolean(no.genotypes,colnames(no.genotypes),"AND")
no.genotypes[1:5]
sum(no.genotypes)

summary.geno.extra[1:5,]
missing.targets<-c("all","all.filt")

high.missing<-{}
imt<-1
for(imt in 1:length(missing.targets)){
    the.missing.alleles<-paste("MISSING.Alleles",missing.targets[imt],sep=".")
    the.total.alleles<-paste("TOTAL.Alleles",missing.targets[imt],sep=".")
    a.missing.test<-as.numeric(summary.geno.extra[,the.missing.alleles])/(as.numeric(summary.geno.extra[,the.total.alleles])+as.numeric(summary.geno.extra[,the.missing.alleles]))
    
     if(is.null(length(high.missing)) | length(high.missing)==0){
       high.missing<-a.missing.test
       }else{
          high.missing<-cbind(high.missing,a.missing.test)
        }
  }
high.missing.table<-high.missing
colnames(high.missing.table)<-missing.targets
rownames(high.missing.table)<-key
                 
## high.missing<- cbind(as.numeric(summary.geno.extra[,"MISSING.Alleles.LOW"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.LOW"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.LOW"])),
##                      as.numeric(summary.geno.extra[,"MISSING.Alleles.HIGH"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.HIGH"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.HIGH"])),
##                      as.numeric( summary.geno.extra[,"MISSING.Alleles.LOW.pheno"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.LOW.pheno"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.LOW.pheno"])),
##                      as.numeric(summary.geno.extra[,"MISSING.Alleles.HIGH.pheno"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.HIGH.pheno"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.HIGH.pheno"])),
##                      as.numeric(summary.geno.extra[,"MISSING.Alleles.nimblegen"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.nimblegen"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.nimblegen"])),
##                      as.numeric(summary.geno.extra[,"MISSING.Alleles.illumina"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.illumina"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.illumina"]))
##                      )

colnames(high.missing.table)
high.missing.table[1:5,]

high.missing.subset<-c("all","all.filt")

high.total.missing<-subset(high.missing.table,select=c("all","all.filt"))
high.total.missing[1:5,]



high.total.missing<-high.total.missing <= missing.threshold
ok.missing<-combine.boolean(high.total.missing,c("all","all.filt"),"AND")
sum(ok.missing)

missing.targets
the.projects
colnames(high.missing.table)

## nimblegen.total.missing<-subset(high.missing.table,select=c("nimble"))
## nimblegen.total.missing[1:5,]
## nimblegen.total.missing<-nimblegen.total.missing > missing.threshold.nimblgen
## ## nimblegen.total.missing<-combine.boolean(high.total.missing,c("LOW","HIGH","LOW.pheno","HIGH.pheno"),"OR")
## sum(nimblegen.total.missing)


## illumina.total.missing<-subset(high.missing.table,select=c("illumina"))
## illumina.total.missing[1:5,]
## illumina.total.missing<-illumina.total.missing > missing.threshold.illumina
## ## nimblegen.total.missing<-combine.boolean(high.total.missing,c("LOW","HIGH","LOW.pheno","HIGH.pheno"),"OR")
## sum(illumina.total.missing)

#sum(high.total.missing | nimblegen.total.missing | illumina.total.missing)
high.missing[1:5,]
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


# maf.lt.all<-a.indel[,colnames(a.indel)[grepl("MAF.lt:",colnames(a.indel))]]
maf.lt.all[1:5,]
#maf.lt.all<-as.logical(maf.lt.all)
as.logical(maf.lt.all[1:5,5])


####################################################################################
##################################### make the POSITION filter matrix QUALITY.FILTER FUNATIOAL group-MAF done here
   global.labs[!(global.labs %in% c(colnames(a.indel),colnames(summary.geno.extra))  )]
 if(sum( !(global.labs %in% c(colnames(a.indel),colnames(summary.geno.extra)))  )>0){print(paste("WARNING postion filters missing for data",
          toString(global.labs[!(global.labs %in% c(colnames(a.indel))) ])))}
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
if( sum( !(the.types %in% c("NA","exonic",possible.mutations)))>0  ){print("WARNING ANNOVAR HAS NEW MUTATION TYPES DEFINED- REVIEW line 532")
                                                                print( the.types[!(the.types %in% c("exonic",possible.mutations))] )   }







#  filter.table.pholy[1:5,]
wanted.muts.coding<-test.for.coding.type(a.indel,geneanno.DB,interesting.coding.mutations)
sum(wanted.muts.coding)
   ####################### mhairi CHECK
wanted.muts.coding.vep<-test.wanted.mutation(a.indel[,"Consequence.Embl"],vep.coding,delimit.by=",")  # filter.table.pholy[,"Consequence"] %in%  vep.coding
  ## a.indel[wanted.muts.coding.vep,]
sum(wanted.muts.coding.vep)
  
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
for(itype in 1:length(a.type)){
  the.test<-grepl(a.type[itype],a.indel[,"gene_biotype"])
  the.test[is.na(the.test)]<-FALSE
  wanted.muts.NONcoding.keep<- wanted.muts.NONcoding.keep | the.test
}

 #  HERE wanted.muts.NONcoding.keep JUST DENOTES  "miRNA" &  "lincRNA" at this point USE wanted.muts.NONcoding to restrict to exones and splice BELOW
  wanted.muts.NONcoding.keep<-wanted.muts.NONcoding.keep & wanted.muts.NONcoding

dim(qual)
colnames(qual)
qual[1:5,]

bad.coding<-wanted.muts.coding | wanted.muts.coding.vep | (wanted.muts.NONcoding.keep & (qual[,"GERP.low"] | qual[,"GERP.unknown"]) ) | ( (wanted.interesting.to.prefilter | wanted.interesting.to.prefilter.vep) & (qual[,"GERP.high"]) )
sum(bad.coding)

bad.coding<-wanted.muts.coding | wanted.muts.coding.vep 
bad.non.coding<-(wanted.muts.NONcoding.keep & (qual[,"GERP.low"] | qual[,"GERP.unknown"]) ) | ( (wanted.interesting.to.prefilter | wanted.interesting.to.prefilter.vep) & (qual[,"GERP.high"]) )
bad.effect<-bad.coding | bad.non.coding


# bad.coding<-test.for.coding.type(a.indel,geneanno.DB,c("stopgain SNV","stoploss SNV","frameshift deletion","frameshift insertion"))
#bad.frame<-test.for.coding.type(geneanno.table,geneanno.DB,c("frameshift deletion","frameshift insertion"))


## basic.qual<-combine.boolean(qual,c("QUAL", "QD", "HRun", "SB"),"AND")
## gatk.qual<-combine.boolean(qual,c("FILTER_PASS", "FILTER_100" ),"OR")
## full.qual<-combine.boolean(cbind(basic.qual,gatk.qual),"all","OR")


## basic.qual<-combine.boolean(qual,c("QUAL", "QD", "HRun", "SB","FILTER_100"),"AND")
## gatk.qual<-qual[,"FILTER_PASS"]
## full.qual<-combine.boolean(cbind(basic.qual,gatk.qual),"all","OR")

full.qual<-qual[,"FILTER_PASS"]
sum(full.qual)

## full.qual<-gatk.qual
## sum(full.qual)
#full.qual<-qual[,"FILTER_PASS"]
sum(full.qual)



#any.functional<-combine.boolean(qual,c("PolyPhen.low","SIFT.high","mut.taster.high","phylo.high","PolyPhen.bad","SIFT.bad","GERP.high","ljb_gerp.high"),"OR")
functional<-combine.boolean(qual,c("PolyPhen.low","SIFT.high","PolyPhen.bad","SIFT.bad","GERP.high"),"OR")
sum(functional)
functional<-functional | bad.coding # include really bad protein changes
sum(functional)


## maf.lt.all[1:5,]
## maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.5"])

#pass<- rare.in.group & !no.genotypes & !high.missing & common.loci

##########################################################################
##########################################################################
##########################################################################

p<-0.01
sd.thresh<-6

p

colnames(summary.geno.extra)
######################################
n<-max(as.integer(summary.geno.extra[,"TOTAL.Alleles.all"]))
n
alt.counts.thresh<-1
while( (alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) <= sd.thresh){alt.counts.thresh<-alt.counts.thresh+1}
alt.counts.thresh

summary.geno.extra[1:5,]
rare.in.all<-as.numeric(summary.geno.extra[,"ALT.Alleles.all"])< alt.counts.thresh
rare.in.all.filt<-as.numeric(summary.geno.extra[,"ALT.Alleles.all.filt"])< alt.counts.thresh
sum(rare.in.all)
sum(rare.in.all.filt)
names(rare.in.all)<-key
names(rare.in.all.filt)<-key
##############################



## n<-max(as.integer(summary.geno.extra[,"TOTAL.Alleles.low"]))
## n
## alt.counts.thresh<-1
## while( (alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) <= sd.thresh){alt.counts.thresh<-alt.counts.thresh+1}
## alt.counts.thresh

## summary.geno.extra[1:5,]
## rare.in.low<-as.numeric(summary.geno.extra[,"ALT.Alleles.low"])< alt.counts.thresh
## rare.in.low.filt<-as.numeric(summary.geno.extra[,"ALT.Alleles.low.filt"])< alt.counts.thresh
## sum(rare.in.low)
## sum(rare.in.low.filt)
## names(rare.in.low)<-key
## names(rare.in.low.filt)<-key
## ##############################

## n<-max(as.integer(summary.geno.extra[,"TOTAL.Alleles.high"]))
## n
## alt.counts.thresh<-1
## while( (alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) <= sd.thresh){alt.counts.thresh<-alt.counts.thresh+1}
## alt.counts.thresh

## summary.geno.extra[1:5,]
## rare.in.high<-as.numeric(summary.geno.extra[,"ALT.Alleles.high"])< alt.counts.thresh
## rare.in.high.filt<-as.numeric(summary.geno.extra[,"ALT.Alleles.high.filt"])< alt.counts.thresh
## sum(rare.in.high)
## sum(rare.in.high.filt)
## names(rare.in.high)<-key
## names(rare.in.high.filt)<-key
## ##############################


#length(maf.filter)
#length(rare.in.controls)

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
rownames(a.indel)[1:20]

#################### in repeats looking  forward

chk.in.repeat<-large.indel & !are.repeats
sum(chk.in.repeat)
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



##########################################################################
##########################################################################
##########################################################################
not.flat.genotype<-!qual[,"flat"]
sum(not.flat.genotype)

colnames(qual)

is.unwound.geno<-grepl("snp:\\d+$",a.indel[,"TYPE"]) | grepl("indel:\\d+$",a.indel[,"TYPE"])
#a.indel[!not.flat.genotype,"TYPE"]



## hwe.control.threshold
## hw.p.control.filt[1:5]
## hw.p.control[1:5]
## hw.p.ex.control.filt[1:5]
## hw.p.ex.control[1:5]





hw.controls.ok<-hw.p.control > hwe.control.threshold
hw.controls.ok.filt<-hw.p.control.filt > hwe.control.threshold




#hw.controls.ok[loci]

in.common.hit.gene <- a.indel[,"Gene.Names"] %in% common.hit.genes
in.common.hit.gene[1:5] 
sum(in.common.hit.gene)
the.chr
on.x.y<-a.indel[,"chr"] %in% c("X","Y","23","24","chrX","chrY")
sum(on.x.y)

#table(a.indel[,"TYPE"])
snp.only<-grepl("^snp",a.indel[,"TYPE"])



##########################################################################
##########################################################################
##########################################################################
icc<-3
 if(target.pheno.col %in% case.control){
   for(icc in 1:length(case.control.classes)){
     recode<-  pheno[,target.pheno.col] %in% names(case.control.classes)[icc]
     pheno[recode,target.pheno.col]<-as.numeric(case.control.classes[icc])
   }}
  
pheno[,target.pheno.col]<-as.numeric(pheno[,target.pheno.col])
formula

pheno[1:5,]
table(pheno[,target.pheno.col])

## pheno.use<-pheno[pheno[,target.pheno] %in% c(0,1),]
## dim(pheno.use)
the.samples.use<-pheno[,"SAMPLE"]
the.samples.use<-paste(the.samples.use,".GT",sep="")




##########################################################################
##########################################################################
##########################################################################
types<-c("all","non.coding","coding")

itypes<-1
for(itypes in 1:length(types)){

snap.file<-types[itypes]

if(types[itypes]=="all"){
pass<- full.qual &  bad.effect & maf.filter   & !in.common.hit.gene  & !unannotated.hits & not.flat.genotype & !are.repeats & !are.in.repeats &
ok.missing & hw.controls.ok.filt & !no.genotypes &  rare.in.all #  & !on.x.y
}

if(types[itypes]=="non.coding"){
pass<- full.qual &  bad.non.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype & !are.repeats  &
ok.missing & hw.controls.ok.filt & !no.genotypes &  rare.in.all #  & !on.x.y
}

if(types[itypes]=="coding"){
pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype & !are.repeats &
ok.missing & hw.controls.ok.filt & !no.genotypes &  rare.in.all #  & !on.x.y
}



sum(pass)

help<-cbind(full.qual,bad.coding,maf.filter,rare.in.group,no.genotypes,in.common.hit.gene ,hw.controls.ok,on.x.y,unannotated.hits,not.flat.genotype,are.repeats,are.in.repeats,ok.missing,ok.missing,is.unwound.geno,is.unwound.geno ,hw.p.control.filt)

## help<-cbind(full.qual,bad.coding,maf.filter,rare.in.group,no.genotypes,in.common.hit.gene ,hw.controls.ok,on.x.y,unannotated.hits,not.flat.genotype,are.repeats,are.in.repeats,ok.missing,ok.missing,is.unwound.geno,is.unwound.geno ,hw.p.control.filt,rare.in.group.filt,no.genotypes.filt,rare.in.controls.filt )
## pass<-full.qual & functional & maf.filter & rare.in.group & !no.genotypes  & not.flat.genotype & !(high.total.missing | nimblegen.total.missing | illumina.total.missing)
## sum(pass)

## sum(pass & not.flat.genotype)
## a.indel[(pass & not.flat.genotype),"TYPE"]


################################# GEFOS FILTERING cause sending all
#pass<-pass[the.snps] ### GEOFS 
genotypes<-a.indel[pass,the.samples.use] ## ordered correctly for phenotypes
snp.names<-key[pass] ## GEFOS ony name with start

#### snpinfo now A different size than a.indel since added pathways!!!
snpinfo<-snpinfo.ori[snpinfo.ori[,"Name"] %in% snp.names,]
if( sum(!(snp.names %in% snpinfo.ori[,"Name"]))>0){print("WARINING snp.names not in snpinfo- unusual!")}
dim(snpinfo)
length(snp.names)
dim(genotypes)


dim(genotypes)

print("start QC")
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


#tapply(as.vector(genotypes),as.vector(genotypes),length)
num.col<-dim(genotypes)[2]
num.row<-dim(genotypes)[1]
## genotypes[1:5,1:20]
genotypes<-as.numeric(as.matrix(genotypes))
dim(genotypes)<-c(num.row,num.col)
genotypes<-t(genotypes) # samples x SNPS

colnames(genotypes)<-snp.names
rownames(genotypes)<-gsub(".GT$","",the.samples.use)



########## make snpInfo
#tapply(gene.list[,"CHR"],gene.list[,"CHR"],length)
## gene.list[gene.list[,"CHR"]=="Y",]
## a.indel[the.snps,1:20][1:12,]

#a.indel[pass,][1:5,1:10]
 if(target.pheno.col %in% case.control){
 cohort.seq <- skatCohort(Z=genotypes,formula, SNPInfo = snpinfo, data=pheno,aggregateBy="cluster",family=binomial(),verbose=FALSE)

}else{
cohort.seq <- skatCohort(Z=genotypes,formula, SNPInfo = snpinfo, data=pheno,aggregateBy="cluster",family=gaussian(),verbose=FALSE) ## genes and clusters
}

meta.results.burden<-burdenMeta(cohort.seq,wts=1,mafRange = c(0,1),SNPInfo = snpinfo,aggregateBy="cluster")
meta.results.skat<-skatMeta(cohort.seq,SNPInfo = snpinfo,aggregateBy="cluster")
meta.results.skatO<-skatOMeta(cohort.seq,burden.wts =1,SNPInfo = snpinfo,aggregateBy="cluster")


#meta.results.skat<-{}
#meta.results.skatO<-{}

the.order<-     order(meta.results.burden[,"p"])
sum(is.na(meta.results.burden[,"p"])) ## bad p-values shoudl not happen
meta.results.burden<-  meta.results.burden[the.order,]

## meta.results.burden[1:50,]
## meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]

the.order<-     order(meta.results.skat[,"p"])
meta.results.skat<-  meta.results.skat[the.order,]
meta.results.skat[1:50,]

the.order<-     order(meta.results.skatO[,"p"])
sum(is.na(meta.results.skatO[,"p"])) ## bad p-values shoudl not happen
meta.results.skatO<-  meta.results.skatO[the.order,]
meta.results.skatO[1:50,]

## meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,]
## meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]
clusters[1:5,]

meta.results.burden[meta.results.burden[,"gene"] %in% clusters[,1],]
meta.results.burden[meta.results.burden[,"gene"] %in% clusters[,2],]


## meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,]



## setwd(analysis.dir) meta.results.skat<-{}
## getwd()
## bad.non.coding
## snap.file<-"coding.0.01.all.geno.all.filters"
## snap.file<-"coding.0.001.all.geno.all.filters_no.imput"

## write.table(meta.results.burden[1:50,],file=paste("Burden","Top50",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,],file=paste("Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

snap.file
write.table(meta.results.burden,file=paste("Burden",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skat,file=paste("Skat",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skatOfile=paste("SkatO",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,],file=paste("SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


#annotations<-a.indel[,c(1:6,16,28,7,30,34,37:42,43,14,32,33)]

## save(list=c("meta.results.skat","meta.results.skatO","meta.results.burden","pheno.use","snpinfo","genotypes","pass","high.missing","annotations","help","key","summary.geno.extra"),file=paste(paste(project.files[ichr],".small.RData",sep="")) )



save(list=c("case.control","snpinfo.ori","formula","clusters","pheno.types","ipheno","clusters.wanted","genotypes","p","meta.results.skat","meta.results.skatO","meta.results.burden","pheno","target.pheno.col","snpinfo","fil.genotypes","pass","high.missing.table","a.indel","help","key","summary.geno.extra","full.qual","bad.effect","maf.filter","in.common.hit.gene","on.x.y","unannotated.hits","not.flat.genotype","are.repeats","are.in.repeats","ok.missing","hw.controls.ok.filt","no.genotypes","rare.in.all","are.in.repeats.back","are.in.repeats.forward"),file=paste(paste(project.files[ichr],".",pheno.types[ipheno],".",snap.file,".small_final.RData",sep="")) )

} # itype

print(paste("Done: ",pheno.types[ipheno],"->",project.files[ichr]))
#save.image(file=paste(project.files[ichr],".",pheno.types[ipheno],".",snap.file,".RData",sep=""))

} # pheno loop

} # loop over projects


} # loop over fam


> getwd()
[1]


## save.image("UQ_mercury_LATEST.RData")
## load("/media/UQCCG/Sequencing/Projects/BONE-GENOME/Analysis/UQ_mercury_LATEST.RData")
## load("/media/UQCCG/Sequencing/Projects/BONE-GENOME/Analysis/Carrie.RData")




## save(list=c("cpheno"),file="pheno.Rdata")

## load("pheno.Rdata")


## pheno[1:5,]
## class(cpheno)
## class(pheno)
## identical(pheno,cpheno)
## sum(pheno[,1]!=cpheno[,1])
## pheno[,1]!=cpheno[,1]

## pheno[20:30,]
## cpheno[20:30,]

## grep("MrOS025",cpheno[,1])
## grep("MrOS025",pheno[,1])

## dim(genotypes)
## grep("MrOS025",rownames(genotypes))
## rownames(genotypes)==pheno[,"SAMPLE"]


## rownames(genotypes)==pheno[,"SAMPLE"] ## must be all true
## posns<-match(rownames(genotypes),pheno[,"SAMPLE"])
## missing<-is.na(posns)
## sum(missing) # check to be sure is 0 else a sample is missing
## pheno<-pheno[posns,]
## rownames(genotypes)==pheno[,"SAMPLE"]

################################## RESTART LOCATION ####################
################################## RESTART LOCATION ####################
################################## RESTART LOCATION ####################
################################## RESTART LOCATION ####################
################################## RESTART LOCATION ####################
################################## RESTART LOCATION ####################


analysis.dir<-"/media/UQCCG/Sequencing/Projects/BONE-GENOME/Analysis"
setwd(analysis.dir)

library(skatMeta)  ### you need this R-library
load("MERCURY_BONE.RData")

traits<-read.delim("Bone_samples2.csv",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
traits[1:5,]
pheno[1:5,]

posns<-match(pheno[,"PATIENT"],traits[,"SAMPLE"])
missing<-is.na(posns)
sum(missing) ## should be zero else missing samples
traits<-traits[posns,] # pheno and traits now in the same order for samples

pheno[,3]!=traits[,2]
pheno[,"SampleProject"]<-as.numeric(traits[,"TRAIT"])

hist(traits[,"TRAIT"])

HIGH<-as.numeric(pheno[,"SampleProject"])>0
LOW<-as.numeric(pheno[,"SampleProject"])<=0

pheno<-cbind(pheno,HIGH,LOW)
################################## RESTART LOCATION ####################



######################  CARRIE YOU NEED TO ADD THE *RESIDUAL BMD* TO COLUMN **SampleProject** in pheno: (it's all NA at the moment)
pheno[1:5,]
## > pheno[1:5,]
##    SAMPLE PATIENT SampleProject FamilyCode ParticipantCode PaternalID MaternalID Sex AffectionStatus  all
## 1 MrOS001 MrOS001       5.16962        ALL         MrOS001          0          0   0               1 TRUE
## 2 MrOS002 MrOS002      -3.90818        ALL         MrOS002          0          0   0               1 TRUE
## 3 MrOS003 MrOS003       3.78133        ALL         MrOS003          0          0   0               1 TRUE
## 4 MrOS004 MrOS004      -3.06563        ALL         MrOS004          0          0   0               1 TRUE
## 5 MrOS005 MrOS005       3.67238        ALL         MrOS005          0          0   0               1 TRUE

covars<-"1"
target.pheno.col<-"SampleProject"
formula<-paste(target.pheno.col,"~",paste(covars,collapse="+"),sep="")
print(formula)
formula<-formula(formula)
#### ADD the residual BMD to the SampleProject column:
formula # this is the formule # skat will use
a.indel[1:5,1:50]
dim(a.indel) # 440130   1379 all the annotated genotypes

############################ this is for the coding mtation with MAF< 0.01

snap.file<-"MrOS.coding_WORKING.0.01" ### used to name output files 


pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype & !are.repeats &
ok.missing & hw.controls.ok.filt & !no.genotypes &  rare.in.all #  & !on.x.y

sum(pass) # 56385

help<-cbind(full.qual,bad.coding,maf.filter,rare.in.group,no.genotypes,in.common.hit.gene ,hw.controls.ok,on.x.y,unannotated.hits,not.flat.genotype,are.repeats,ok.missing,ok.missing,is.unwound.geno,is.unwound.geno ,hw.p.control.filt, rare.in.all) ## summary of QC




################################# GEFOS FILTERING cause sending all

genotypes<-a.indel[pass,the.samples.use] ## ordered correctly for phenotypes
snp.names<-key[pass] ## GEFOS ony name with start

#### snpinfo now A different size than a.indel since added pathways!!!
snpinfo<-snpinfo.ori[snpinfo.ori[,"Name"] %in% snp.names,]
if( sum(!(snp.names %in% snpinfo.ori[,"Name"]))>0){print("WARINING snp.names not in snpinfo- unusual!")}
dim(snpinfo)
length(snp.names)
dim(genotypes)


dim(genotypes)

print("start QC")
genotypes[genotypes=="NA"]<-NA
genotypes[genotypes=="0/0"]<-0
genotypes[genotypes=="0/1"]<-1
genotypes[genotypes=="1/1"]<-2

########### prevent any averaging
dim(genotypes)
genotypes[is.na(genotypes)]<-0
dim(genotypes)
########### prevent any averaging


#tapply(as.vector(genotypes),as.vector(genotypes),length)
num.col<-dim(genotypes)[2]
num.row<-dim(genotypes)[1]
## genotypes[1:5,1:20]
genotypes<-as.numeric(as.matrix(genotypes))
dim(genotypes)<-c(num.row,num.col)
genotypes<-t(genotypes) # samples x SNPS

colnames(genotypes)<-snp.names
rownames(genotypes)<-gsub(".GT$","",the.samples.use)


#a.indel[pass,][1:5,1:10]
###  if(target.pheno.col %in% case.control){
###  cohort.seq <- skatCohort(Z=genotypes,formula, SNPInfo = snpinfo, data=pheno,aggregateBy="cluster",family=binomial(),verbose=FALSE)

### }else{
cohort.seq <- skatCohort(Z=genotypes,formula, SNPInfo = snpinfo, data=pheno,aggregateBy="cluster",family=gaussian(),verbose=TRUE) ## genes and clusters
#}

meta.results.burden<-burdenMeta(cohort.seq,wts=1,mafRange = c(0,1),SNPInfo = snpinfo,aggregateBy="cluster")
meta.results.skat<-skatMeta(cohort.seq,SNPInfo = snpinfo,aggregateBy="cluster")
meta.results.skatO<-skatOMeta(cohort.seq,burden.wts =1,SNPInfo = snpinfo,aggregateBy="cluster")


#meta.results.skat<-{}
#meta.results.skatO<-{}

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

clusters[1:5,]

##### write results-----------------------

write.table(meta.results.burden,file=paste("Burden",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skat,file=paste("Skat",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skatO,file=paste("SkatO",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)




#################### UNWIND GENE CODE #############
#############################################################################
#############################################################################
#############################################################################
##########################################################################################################################################################
#############################################################################
#############################################################################
### unwind the top 50 genes


analysis.dir<-"/media/UQCCG/Sequencing/Projects/BONE-GENOME/Analysis"  ## directory where "UQDI_murcury.RData" and you sample/trits file located
setwd(analysis.dir)
load("UQDI_murcury_LATEST.RData")

################################## RESTART LOCATION ####################



the.top<-50
to.unwind<-unique(c(meta.results.burden[1:the.top,"gene"],meta.results.skatO[1:the.top,"gene"]))
to.unwind.name<-"TOP_50" ## refic for output file name

############ or you can do a single gene or many genes this way
## to.unwind<-c("BCL2L12","SETD8") #"test.set"
## to.unwind.name<-"TEST_GENE_LIST"
##############################################
## /media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/Analyse_project_MERGE_AOGC_ANALYSIS_NEW_PAPER_CHIP.r
## /media/UQCCG/Sequencing/Projects/BONE-GENOME/Analysis/UQ_mercury_LATEST.RData
## /media/UQCCG/Sequencing/Projects/BONE-GENOME/Analysis/TOP_50.GENOTYPE.conponents:.Burden.clusters.MrOS.coding_WORKING.0.01.txt
## /media/UQCCG/Sequencing/Projects/BONE-GENOME/Analysis/meta.analysis/META.STDERR.TRAIT.p_0.01.SUMMARY_UQDI_run.xlsx
## /media/UQCCG/Sequencing/Projects/BONE-GENOME/Analysis/meta.analysis/META.SAMPLE.SIZE.TRAIT.SKATO.WITH_BENIGN_MISSENSE_p_0.01.SUMMARY_UQDIrun.xlsx

snpinfo.ex<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,]
loci<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,"Name"] # this is IDH1 not IDH1 in cluster # are the snp.names
the.genes<-unique(snpinfo.ex[,"gene"])
the.genes<-the.genes[!(the.genes %in% clusters.wanted)]

the.genes #

############repest to clean out cluster names 
to.unwind<-the.genes
snpinfo.ex<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,]
loci<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,"Name"] # this is IDH1 not IDH1 in cluster # are the snp.names
the.genes<-unique(snpinfo.ex[,"gene"])
the.genes<-the.genes[!(the.genes %in% clusters.wanted)]

the.genes
the.genes.burden<-meta.results.burden[meta.results.burden[,"gene"] %in% the.genes,]

the.genes.burden
write.table(the.genes.burden,file=paste(to.unwind.name,"conponents:","Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

the.genes.burden<-meta.results.skatO[meta.results.skatO[,"gene"] %in% the.genes,]
the.genes.burden
write.table(the.genes.burden,file=paste(paste(to.unwind.name,collapse="."),"conponents:","SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


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
dim(pheno)
snpinfo.ex[1:5,]

####################### do single point analysis

cohort.seq.ex <- skatCohort(genotypes.ex, formula, SNPInfo = snpinfo.ex, data=pheno,aggregateBy = "Name",verbose=FALSE)
## meta.results.skat.ex<-skatMeta(cohort.seq,SNPInfo = snpinfo)
meta.results.burden.ex<-burdenMeta(cohort.seq.ex,wts=1,mafRange = c(0,1),SNPInfo = snpinfo.ex,aggregateBy = "Name")
meta.results.burden.ex[1:5,]
pheno[1:5,]


 if(target.pheno.col %in% case.control){
 cohort.seq.test <- skatCohort(genotypes.ex, formula, SNPInfo = snpinfo.ex, data=pheno,aggregateBy = "cluster",verbose=FALSE)
}else{
cohort.seq.test <- skatCohort(genotypes.ex, formula, SNPInfo = snpinfo.ex, data=pheno,aggregateBy = "cluster",family=gaussian(),verbose=FALSE)
}




meta.results.burden.test<-burdenMeta(cohort.seq.test,wts=1,mafRange = c(0,1),SNPInfo = snpinfo.ex,aggregateBy = "cluster")
#meta.results.burden.test


muts.in.cases<-apply(genotypes.ex[pheno[,"HIGH"],],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
muts.in.controls<-apply(genotypes.ex[pheno[,"LOW"],],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})


##controls<- paste(pheno[pheno[,"Control"],"SAMPLE"],".GT",sep="")
## a.indel[figure,controls]
##   table(a.indel[figure,controls][2,])  

## muts.in.cases<-apply(genotypes.ex[pheno[,"SampleProject"]==1,],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
## muts.in.controls<-apply(genotypes.ex[pheno[,"SampleProject"]==0,],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})

figure<- match(loci,key)


########################################################
check<-16

quality.cases<-rep("",times=length(loci))
quality.controls<-rep("",times=length(loci))

depth.cases<-rep("",times=length(loci))
depth.controls<-rep("",times=length(loci))

a.indel.sub<-a.indel[figure,]

check<-1
for(check in 1:length(loci)){
# print(check)
#check<-"chr11:130066457:130066457:-:A:indel"
# posn<-grep(loci[check],key)
posn<-check

if(muts.in.cases[check]!=""){
#the.gt<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GT",sep=".")
the.gq<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GQ",sep=".")
the.ad<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"AD",sep=".")

quality.cases[check]<-paste(a.indel.sub[posn,the.gq],collapse=",")
depth.cases[check]<-paste(a.indel.sub[posn,the.ad],collapse=";")

a.indel[posn,the.gq]
## a.indel[posn,the.gt]
## a.indel[posn,the.dp]
}

if(muts.in.controls[check]!=""){
#the.gt<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"GT",sep=".")
the.gq<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"GQ",sep=".")
the.ad<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"AD",sep=".")

quality.controls[check]<-paste(a.indel.sub[posn,the.gq],collapse=",")
depth.controls[check]<-paste(a.indel.sub[posn,the.ad],collapse=";")

a.indel[posn,the.gq]
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
annotations<-a.indel[,c(1:6,8,11,16,28,7,30,34,35,36,37:42,43,14,32,33)]
dim(annotations)
dim(help)
dim(summary.geno.extra)
length(figure)

                             
sum(meta.results.burden.ex[,"gene"]!=loci)
dim(meta.results.burden.ex)
a.functions<-a.indel[,c("PolyPhen.scores","SIFT.scores","PolyPhen.desc","SIFT.desc")]
#,is.benign.missense[figure]
                             
out<-cbind(meta.results.burden.ex,annotations[figure,],a.functions[figure,],annotations[figure,],summary.geno.extra[figure,colnames(summary.geno.extra)[grep("^GENO",colnames(summary.geno.extra))]],help[figure,],muts.in.cases,quality.cases,depth.cases,muts.in.controls,quality.controls,depth.controls)


dim(out)
out[1:4,]

getwd()
setwd(analysis.dir)

paste(paste(to.unwind,collapse="."))
paste(to.unwind.name,collapse=".")
  paste(paste(to.unwind.name,collapse="."),"GENOTYPE.conponents:","SkatO","clusters",snap.file,"txt",sep=".")
write.table(out,file=paste(paste(to.unwind.name,collapse="."),"GENOTYPE.conponents:","Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


#######################################################################
#######################################################################
#######################################################################





#############################################################################
#############################################################################
#############################################################################
#############
#############################################################################################################################################
#############################################################################
#############################################################################
################ META ANALYSIS


AOGC at HIP:
 "Number Samples: 949"

MRoS
 "Number Samples: 200"

########################################################################## SkatO.MrOS.coding.0.01.txt


## SKAT 0.01  USING SAMPLE sizes for weighting
setwd("/media/UQCCG/Sequencing/Projects/BONE-GENOME/Analysis/meta.analysis")
#files<-c("AOGC_HIP.BMD_EFF_STD_HIP.coding.assoc_summary.txt","SkatO.MrOS_4MAR2015.txt")
files<-c("AOGC_HIP.BMD_EFF_STD_HIP.coding.assoc_summary.txt","SkatO.MrOS.coding_WORKING.0.01.csv")
sizes<-c(949,208) ### add the sample sizes to abobe files (using excell) in a column called nmiss
traits<-c("BONE")
i<-1


system(paste("sed s/chip.TRAIT/",files[1],"/ run_config_for_P_and_BETA_TEMPLATE.txt > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.TRAIT/",files[2],"/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("sed s/chip.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt1",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
paste("sample.size.CONFIG",traits[i],"txt1",sep=".")
system(paste("./metal ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("cp ","META.SAMPLE.SIZE.TRAIT1.tbl","META.SAMPLE.SIZE.TRAIT.SKATO.WITH_BENIGN_MISSENSE_p_0.01.txt", sep=" "))


aogc<-read.delim(files[1],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
MrOs<-read.delim(files[2],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
meta<-read.delim("META.SAMPLE.SIZE.TRAIT.SKATO.WITH_BENIGN_MISSENSE_p_0.01.txt",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)

aogc[1:5,]
MrOs[1:5,]
meta[1:5,]

colnames(aogc)<-paste(colnames(aogc),"AOGC",sep=".")
colnames(MrOs)<-paste(colnames(MrOs),"MROS",sep=".")

posns<-match(meta[,"MarkerName"], aogc[,1])
missing<-is.na(posns)
sum(missing)
sum(!missing)
#meta[missing,"MarkerName"]
aogc<-aogc[posns,]

posns<-match(meta[,"MarkerName"], MrOs[,1])
missing<-is.na(posns)
sum(missing)
#meta[missing,"MarkerName"]
MrOs<-MrOs[posns,]

dim(MrOs)
dim(aogc)

meta<-cbind(meta,aogc,MrOs)
meta[1:5,]
dim(meta)
dim(aogc)
getwd()
write.table(meta,file="META.SAMPLE.SIZE.TRAIT.SKATO.WITH_BENIGN_MISSENSE_p_0.01.SUMMARY_UQDIrun.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

#############################################################################################################



## BURDEN 0.01  USING inverse variance
setwd("/media/UQCCG/Sequencing/Projects/BONE-GENOME/Analysis/meta.analysis")
files<-c("AOGC_HIP.BMD_EFF_STD_HIP.coding.assoc_summary.txt","Burden.MrOS.coding_WORKING.0.01.csv")
#files<-c("AOGC_HIP.BMD_EFF_STD_HIP.coding.assoc_summary.txt","Burden.MrOS_4MAR2015.txt")
sizes<-c(949,200) ### add the sample sizes to abobe files (using excell) in a column called nmiss
traits<-c("BONE")
i<-1

#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-11-04_AML_TCGA_Replication/Analysis/meta_analyis/ALL_clusters.conponents:.Burden.clusters.coding.0.001.all.geno.all.filters_no.imput_paper.txt

system(paste("sed s/chip.TRAIT/",files[1],"/ run_config_for_inverse_varience_TEMPLATE.txt > ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.TRAIT/",files[2],"/ ",paste("STDERR.CONFIG",traits[i],"txt",sep=".")," > ",paste("STDERR.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("sed s/chip.NMISS/","nmiss","/ ",paste("STDERR.CONFIG",traits[i],"txt1",sep=".")," > ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.NMISS/","nmiss","/ ",paste("STDERR.CONFIG",traits[i],"txt",sep=".")," > ",paste("STDERR.CONFIG",traits[i],"txt1",sep="."),sep=""))
paste("STDERR.CONFIG",traits[i],"txt1",sep=".")
system(paste("./metal ",paste("STDERR.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("cp ","META.STDERR.TRAIT1.tbl","META.STDERR.TRAIT.p_0.01.txt", sep=" "))



aogc<-read.delim(files[1],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
MrOs<-read.delim(files[2],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
meta<-read.delim("META.STDERR.TRAIT.p_0.01.txt",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)



aogc[1:5,]
MrOs[1:5,]
meta[1:5,]

colnames(aogc)<-paste(colnames(aogc),"AOGC",sep=".")
colnames(MrOs)<-paste(colnames(MrOs),"MROS",sep=".")

posns<-match(meta[,"MarkerName"], aogc[,1])
missing<-is.na(posns)
sum(missing)
# meta[missing,"MarkerName"]
aogc<-aogc[posns,]

posns<-match(meta[,"MarkerName"], MrOs[,1])
missing<-is.na(posns)
sum(missing)
# meta[missing,"MarkerName"]
MrOs<-MrOs[posns,]

dim(MrOs)
dim(aogc)

meta<-cbind(meta,aogc,MrOs)
meta[1:5,]
getwd()
write.table(meta,file="META.STDERR.TRAIT.p_0.01.SUMMARY_UQDI_run.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

































