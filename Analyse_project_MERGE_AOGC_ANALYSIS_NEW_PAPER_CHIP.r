
######## ONLY NEED TO CHOOSE A DIRECTORY AND EXTENSIONS - used tab delimited files 
###############################################
###### THIS IS the super annotion run USE ALL THE FILTER AND NOVEL DATABASES
##Build AOGC

##################### if have a genotype component
genotype.file.location<-"/media/UQCCG-Analysis/AOGC_exome_chip/working_genotypes"
genotype.file.prefix<-"recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL"

annotate.dir<-"/media/UQCCG-Analysis/AOGC_exome_chip/working_genotypes/Annotate" # dir(annotate.dir)
analysis.dir<-"/media/UQCCG-Analysis/AOGC_exome_chip/working_genotypes/Analysis" # dir(analysis.dir)
project.extension<-".bim" ## justX1 the exterion not fam.extension!
project.name<-genotype.file.prefix
fam<-c("") #


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

############################## Sample Sheet

the.sample.sheet<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
## ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
contaminated.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/contaminated_AOGC_SEQ_samples.txt"
contaminated<-read.table(contaminated.file,header=F,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)
remove.cols<-c()


############################################################

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
## print(colnames(geno.aogc))/media/UQCCG/UQCCG-Projects/Popularion stratification/ALL_eth_650Y_forward_hg19_Best.bed
## use.key<-build.key(geno.aogc,core.ann)
## insert.location<-70  ### this is where to add AOGC data INSERTS AFTER THE LOCATION
## ################################ add aogc

## ann<-readRDS("/media/scratch2/AOGC-NGS/Analysis/AOGC_sequnnce_LS/AOGC_sequence_10_LS_ANNOTATION.rds")


geneanno.DB<-c("refGene","knownGene","ensGene") # returns 2 extra columns
names(geneanno.DB)<-c("refGene","knownGene","ensGene")



## geneanno.DB<-c("refGene") # returns 2 extra columns
## names(geneanno.DB)<-c("refGene")





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

pheno.types<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_FN") ## vales is column header
names(pheno.types)<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_FN") ### name is output columns

## case.control.classes<-c(0,1,0)
## names(case.control.classes)<-c("Control","NMD","AOGC")
## case.control.classes



case.control<-c("BMD_AFFSTAT")
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

# ifam<-1
for(ifam in 1:length(fam)){

the.extension<-paste(fam[ifam],project.extension,"$",sep="")
project.files<-files[grepl(the.extension ,files)]
print(sort(paste("Doing: ",project.files,sep=""))) # project.files<-project.files[1:22]

indels<-{}
the.col<-{}
project.files
# ichr<-13
project.files
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
column.labels<-read.delim(project.files[ichr],header=F,nrows=1,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
num.vars<-dim(column.labels)[2]
a.indel<-scan(project.files[ichr],what=character(num.vars),skip=1,sep="\t",fill=TRUE,na.strings="",quote="\"")
num.lines<-length(a.indel)/(num.vars)
dim(a.indel)<-c(num.vars,num.lines)
a.indel<-t(a.indel)
colnames(a.indel)<-column.labels
########################################


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


ann.snp<-
  load("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/Analysis/assoc.ann_quick_chip.RData")

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
num.cores<-8
num.bits<-num.cores
registerDoMC(cores=num.cores)

the.samples<-colnames(a.indel)[grepl(".GT$", colnames(a.indel))]

while((dim(a.indel)[1] %% num.bits)< 2){num.bits<-num.bits+1} ### go don't get matrix issues
num.bits
#(dim(a.indel)[1] %% num.bits)
## fil.genotypes<-foreach(a.indel.bit=iter(a.indel,by='row',chunksize=as.integer(dim(a.indel)[1]/num.bits) ), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype(a.indel.bit,gsub(".GT$","",the.samples),prefix="",suffix="",20,0.02,0.98,0.20,0.80,7,2)

fil.genotypes<-foreach(a.indel.bit=iter(a.indel,by='row',chunksize=as.integer(dim(a.indel)[1]/num.bits) ), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype(a.indel.bit,gsub(".GT$","",the.samples),prefix="",suffix="",20,0.02,0.98,0.20,0.80,10,5)
# 20,0.02,0.98,0.20,0.8,10,5 # for cancers where het may be way out

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


############################### do one phenotype ##################################################################
# ipheno<-1
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

hist(as.numeric(pheno[,pheno.types[ipheno]]))

if(length(unique((pheno[,pheno.types[ipheno]])))>4){ 
  low<-as.numeric(pheno[,pheno.types[ipheno]]) <=0 ## quantative
  high<-as.numeric(pheno[,pheno.types[ipheno]]) >0
}else{ ## case.control
  low<-pheno[,pheno.types[ipheno]]==1
  high<-pheno[,pheno.types[ipheno]]==2
}
sum(low)
sum(high)

illumina<-the.samples %in% ill.samples
nimble<-the.samples %in% nim.samples
all<-rep(TRUE,length(the.samples))

dim(pheno)
length(all)


############## SET up groups to be analysed

pheno<-cbind(pheno,all,low,high,illumina,nimble)

## th project becomes the master copy of the targets below.
the.projects<-c("all","low","high","illumina","nimble")
names(the.projects)<-the.projects
colnames(pheno)


                       

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
rare.in.group.test<-paste("MAF.",c("low","high","low.filt","high.filt"),sep="")
rare.in.group.test
rare.in.group<-combine.boolean(rare.in.group.table,rare.in.group.test,"OR")
sum(!rare.in.group)
rare.in.group[1:10]
## rare.in.group<-combine.boolean(rare.in.group.table,c("MAF.AML","MAF.Control"),"OR")
## rare.in.group.filt<-combine.boolean(rare.in.group.table,c("MAF.AML.filt","MAF.Control.filt"),"OR")
## sum(!rare.in.group)
## sum(!rare.in.group.filt)

no.genotypes.test<-paste("MAF.",c("low","high","low.filt","high.filt"),sep="")
no.genotypes.test
no.genotypes<-(summary.geno.extra[,no.genotypes.test]== 0)  | (is.na( summary.geno.extra[,no.genotypes.test])) # no genotypes in test classes for a mutataion after individaul quality filtering
no.genotypes[1:5,]
no.genotypes<-combine.boolean(no.genotypes,colnames(no.genotypes),"AND")
no.genotypes[1:5]
sum(no.genotypes)

summary.geno.extra[1:5,]
missing.targets<-c("low","high","low.filt","high.filt","illumina","nimble")

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

high.missing.subset<-c("low","high","low.filt","high.filt")

high.total.missing<-subset(high.missing.table,select=c("low","high","low.filt","high.filt"))
high.total.missing[1:5,]



high.total.missing<-high.total.missing <= missing.threshold
ok.missing<-combine.boolean(high.total.missing,c("low","high","low.filt","high.filt"),"AND")
sum(ok.missing)

missing.targets
the.projects
colnames(high.missing.table)

nimblegen.total.missing<-subset(high.missing.table,select=c("nimble"))
nimblegen.total.missing[1:5,]
nimblegen.total.missing<-nimblegen.total.missing > missing.threshold.nimblgen
## nimblegen.total.missing<-combine.boolean(high.total.missing,c("LOW","HIGH","LOW.pheno","HIGH.pheno"),"OR")
sum(nimblegen.total.missing)


illumina.total.missing<-subset(high.missing.table,select=c("illumina"))
illumina.total.missing[1:5,]
illumina.total.missing<-illumina.total.missing > missing.threshold.illumina
## nimblegen.total.missing<-combine.boolean(high.total.missing,c("LOW","HIGH","LOW.pheno","HIGH.pheno"),"OR")
sum(illumina.total.missing)

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



n<-max(as.integer(summary.geno.extra[,"TOTAL.Alleles.low"]))
n
alt.counts.thresh<-1
while( (alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) <= sd.thresh){alt.counts.thresh<-alt.counts.thresh+1}
alt.counts.thresh

summary.geno.extra[1:5,]
rare.in.low<-as.numeric(summary.geno.extra[,"ALT.Alleles.low"])< alt.counts.thresh
rare.in.low.filt<-as.numeric(summary.geno.extra[,"ALT.Alleles.low.filt"])< alt.counts.thresh
sum(rare.in.low)
sum(rare.in.low.filt)
names(rare.in.low)<-key
names(rare.in.low.filt)<-key
##############################

n<-max(as.integer(summary.geno.extra[,"TOTAL.Alleles.high"]))
n
alt.counts.thresh<-1
while( (alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) <= sd.thresh){alt.counts.thresh<-alt.counts.thresh+1}
alt.counts.thresh

summary.geno.extra[1:5,]
rare.in.high<-as.numeric(summary.geno.extra[,"ALT.Alleles.high"])< alt.counts.thresh
rare.in.high.filt<-as.numeric(summary.geno.extra[,"ALT.Alleles.high.filt"])< alt.counts.thresh
sum(rare.in.high)
sum(rare.in.high.filt)
names(rare.in.high)<-key
names(rare.in.high.filt)<-key
##############################


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
pass<- full.qual &  bad.non.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype & !are.repeats & !are.in.repeats &
ok.missing & hw.controls.ok.filt & !no.genotypes &  rare.in.all #  & !on.x.y
}

if(types[itypes]=="coding"){
pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype & !are.repeats & !are.in.repeats &
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

## the.order<-     order(meta.results.burden[,"p"])
## sum(is.na(meta.results.burden[,"p"])) ## bad p-values shoudl not happen
## meta.results.burden<-  meta.results.burden[the.order,]

## meta.results.burden[1:50,]
## meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]

## the.order<-     order(meta.results.skat[,"p"])
## meta.results.skat<-  meta.results.skat[the.order,]
## meta.results.skat[1:50,]

## the.order<-     order(meta.results.skatO[,"p"])
## sum(is.na(meta.results.skatO[,"p"])) ## bad p-values shoudl not happen
## meta.results.skatO<-  meta.results.skatO[the.order,]
## meta.results.skatO[1:50,]

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


## write.table(meta.results.skat[1:50,],file=paste("Skat",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(meta.results.skatO[1:50,],file=paste("SkatO",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(meta.results.skatO[1:50,],file=paste("SkatO",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
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
[1] "/media/old-scratch/media/scratch2/AOGC-NGS/Analysis"



















###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
         gene           p       beta        se    cmafTotal     cmafUsed nsnpsTotal nsnpsUsed nmiss
630    LRRC27 0.000581132 -0.8617283 0.2504827 0.0205479452 0.0205479452         17        17     0
417     RRP12 0.001092519  0.6965400 0.2132995 0.0252897787 0.0252897787         27        27     0
267   DNAJB12 0.003407281  1.0395384 0.3549870 0.0100105374 0.0100105374          5         5     0
196      MBL2 0.004716617  4.3307620 1.5325919 0.0005268704 0.0005268704          1         1     0
202     BICC1 0.008973762 -0.8112837 0.3104736 0.0131717597 0.0131717597          9         9     0
515   CCDC147 0.009978610  0.7571183 0.2938474 0.0147523709 0.0147523709         16        16     0
374       IDE 0.013858515  1.1983774 0.4869652 0.0052687039 0.0052687039         10        10     0
2     ZMYND11 0.014240787  1.5376796 0.6273346 0.0031612223 0.0031612223          3         3     0
601      LHPP 0.015168353  1.0015424 0.4124401 0.0073761855 0.0073761855          5         5     0
91     DNAJC1 0.015806218 -1.4024205 0.5811071 0.0036880927 0.0036880927          5         5     0
315     MAT1A 0.016315103  1.8434004 0.7675113 0.0021074816 0.0021074816          3         3     0




print(paste("Doing phenotype",pheno.types[ipheno]))
target.pheno<-names(pheno.types)[ipheno]
target.pheno.col<-pheno.types[ipheno]











print("Saving results")
out.file<-paste("AOGC_sequence_", the.chr,"_",target.pheno,a.label,".rds", collapse="T", sep="")
out.file.res.skat<-paste("AOGC_sequence_skat", the.chr,"_",target.pheno,a.label,".rds", collapse="T", sep="")
out.file.res.burden<-paste("AOGC_sequence_burden", the.chr,"_",target.pheno,a.label,".rds", collapse="T", sep="")
out.file.res.skatO<-paste("AOGC_sequence_skatO", the.chr,"_",target.pheno,a.label,".rds", collapse="T", sep="")
out.ann.file<-paste("AOGC_sequence_", the.chr,"_",target.pheno,a.label,"_ANNOTATION.rds", collapse="T", sep="")
out.ann.extra.file<-paste("AOGC_sequence_", the.chr,"_",target.pheno,a.label,"_ANNOTATION_EXTRA.rds", collapse="T", sep="")

saveRDS(cohort.seq, file=out.file)
saveRDS(meta.results.skat, file=out.file.res.skat)
saveRDS(meta.results.burden, file=out.file.res.burden)
saveRDS(meta.results.skatO, file=out.file.res.skatO)
saveRDS(ann,file=out.ann.file)
saveRDS(summary.geno.extra.out,file=out.ann.extra.file)

##################################ADD extras here
#setwd(code.dir)
# source("MERGE_alt_alleles_AOGC.r")
#source("MERGE_alt_alleles_AOGC_as_controls.r")

#/media/old-scratch/media/scratch2/AOGC-NGS/Analysis/to_US/AOGC_HIP.BMD_EFF_STD_HIP.all.assoc_summary.txt###############################



if(dont.build.summary){indels<-{}}
## setwd(code.dir)
## source("ucsc.table.names.r") 

} # ipheno

if(!dont.build.summary){
if(ipheno==1){
if(is.null(dim(indels))){
  indels<-a.indel
  the.col<-colnames(a.indel)
                        }else{
                          
  if(sum(!(the.col %in% colnames(a.indel)))>0  | sum(!(colnames(a.indel) %in% the.col))>0){
    print("error colnames don't match")
    print(the.col[!(the.col %in% colnames(a.indel))])
    print(colnames(a.indel)[!(colnames(a.indel) %in% the.col)])
    next
    } # columns don't match
  
  if(sum(!(colnames(a.indel) %in% the.col))>0 ){
                        a.indel<-a.indel[,the.col]     } # reduce a.indels }
                                                                                    
  indels<-rbind(indels,a.indel[,the.col])
     } ## is null so  not first
} # make for one phenotype.
}

} # ichr

setwd(annotate.dir)
#write.table(indels,file="AOGC_pass_filtered.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


} # ifam loop



getwd()
#####################
######################

pheno.types<-c("BMD_EFF_STD_HIP")
names(pheno.types)<-c("HIP")


setwd(analysis.dir)
## files<-dir(getwd())
## files<-files[files %in% "*burden*"]
## files
## burden.files<-grepl("^AOGC_sequence_skatO",files)

## the.extension<-paste(fam[ifam],project.extension,"$",sep="")
## project.files<-files[grepl(the.extension ,files)]
## print(sort(paste("Doing: ",project.files,sep=""))) # project.files<-project.files[1:22]




indels<-{}
the.col<-{}
project.file


ipheno<-1
for(ipheno in 1:length(pheno.types)){
target.pheno<-names(pheno.types)[ipheno]

iregion<-2 # iregion<-1
region.labels<-c("PILOT.GENE.regions","GENE.regions")
region.labels<-c("FILTERED.PILOT.GENE.regions","FILTERED.GENE.regions")
region.labels<-c("PASS.PILOT.GENE.regions","PASS.GENE.regions")
for(iregion in 1:length(region.labels)){
a.label<-region.labels[iregion]


all.chr<-paste("chr",c(1:22,"X","Y"),sep="")
all.chr<-paste("chr",c(1:22,"X"),sep="")
#all.res
# ichr<-10
all.res<-{}
for(ichr in 1:length(all.chr)){
the.chr<-all.chr[ichr]
print(the.chr)
res.skat<-readRDS(  paste("AOGC_sequence_skat", the.chr,"_",target.pheno,a.label,".rds", collapse="T", sep=""))
res.burden<-readRDS(paste("AOGC_sequence_burden", the.chr,"_",target.pheno,a.label,".rds", collapse="T", sep=""))
res.skatO<-readRDS(paste("AOGC_sequence_skatO", the.chr,"_",target.pheno,a.label,".rds", collapse="T", sep=""))
ann.chr<-readRDS(paste("AOGC_sequence_", the.chr,"_",target.pheno,a.label,"_ANNOTATION.rds", collapse="T", sep=""))
ann.extra.chr<-readRDS(paste("AOGC_sequence_", the.chr,"_",target.pheno,a.label,"_ANNOTATION_EXTRA.rds", collapse="T", sep=""))

res.skatO[1:5,]
res.burden[1:5,]
res.skat[1:5,]
ann.chr[1:5,]


order.by<-order(res.skatO[,"p"],decreasing=FALSE)

res.skatO<-res.skatO[order.by,]


posns<-match(res.skatO[,"gene"],res.burden[,"gene"])
res.burden<-res.burden[posns,]

posns<-match(res.skatO[,"gene"],res.skat[,"gene"])
res.skat<-res.skat[posns,]

res.skatO[1:5,]
res.burden[1:5,]
res.skat[1:5,]

colnames(res.skatO)<-paste(colnames(res.skatO),"skatO",sep=".")
colnames(res.burden)<-paste(colnames(res.burden),"burden",sep=".")
colnames(res.skat)<-paste(colnames(res.skat),"skat",sep=".")


a.res<-cbind(the.chr,res.skatO,res.burden,res.skat)
a.res[1:5,]
if(is.null(dim(all.res))){
  all.res<-a.res
  ann<-ann.chr
  ann.extra<-ann.extra.chr
  
}else{
  all.res<-rbind(all.res,a.res)
  ann<-rbind(ann,ann.chr)
    ann.extra<-rbind(ann.extra,ann.extra.chr)
}

}


order.by<-order(all.res[,"p.skatO"],decreasing=FALSE)

all.res<-all.res[order.by,]
all.res[1:10,]



wanted<-c("gene.skatO","p.skatO","p.burden","p.skat","nmiss.skatO","nsnps.skatO","errflag.skatO","beta.burden","se.burden","cmafTotal.burden","cmafUsed.burden","nsnpsTotal.burden","nsnpsUsed.burden","Qmeta.skat","cmaf.skat")


all.res<-all.res[,wanted]

if(iregion==1){
posns<-match(all.res[,"gene.skatO"],ann[,"Gene.Names"])
ann[posns,][1:5,c("Gene.Names","description","gene_biotype","OMIM (Gene::Status::OMIM::description::disease)" )]
all.res<-cbind(ann[posns,c("Gene.Names","description","gene_biotype","OMIM (Gene::Status::OMIM::description::disease)" )],all.res)
all.res.gene<-all.res
}

if(iregion==2){
posns<-match(all.res[,"gene.skatO"],gene.list[,"GENE_NAME"])
gene.list[posns,][1:5,]
gene.list[1:5,]
all.res<-cbind(all.res,gene.list[posns,])
all.res.pilot<-all.res
}


write.table(all.res.gene,file="pass.Gene.based.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(all.res.pilot,file="pass.Gene.EXON.based.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

all.res[1:15,1:5]
all.res.pilot[1:5,]
all.res.gene[1:5,]

region.cols.wanted<-c("gene.skatO","p.skatO","p.burden","p.skat","nmiss.skatO","nsnps.skatO","nsnpsTotal.burden","errflag.skatO")
region.cols.wanted<-c("gene.skatO","p.skatO","p.burden","p.skat","nsnpsTotal.burden")

posns<-match(all.res.gene[,"Gene.Names"],all.res.pilot[,"hgnc_symbol"])

all.res.pilot.subset<-all.res.pilot[posns,region.cols.wanted]
colnames(all.res.pilot.subset)<-paste(colnames(all.res.pilot.subset),"EXONS",sep=".")


cbind(all.res.gene,all.res.pilot[posns,region.cols.wanted])[1:5,]

cbind(all.res.gene[,1:8],all.res.pilot.subset,all.res.gene[,9:dim(all.res.gene)[2]])[1:5,]

gene.region.compare<-cbind(all.res.gene[,1:8],all.res.pilot.subset,all.res.gene[,9:dim(all.res.gene)[2]])
write.table(gene.region.compare,file="pass.Compare.Gene.EXON.based.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


gene.region.compare[1:5,]
lim = 10
qq <- parse.pvals.qq(gene.region.compare[,"p.skat"],lim=lim)
 
ylab = expression(Observed~~-log[10](italic(p)))
xlab = expression(Expected~~-log[10](italic(p)))
 
plot(qq$e,
     qq$o,
     xlim=c(0,lim),ylim=c(0,lim),
     pch=20,col='deepskyblue',
     xlab=xlab,ylab=ylab, main="Gene Based")
     abline(coef=c(0,1), col=1, lwd=2)

savePlot("PassGeneBased.jpg",type="jpeg")

#write.table(geno.all,file=paste(project.name,fam[ifam],"geno.all.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

#geno.all<-read.delim(paste(project.name,fam[ifam],"geno.all.txt",sep="."),header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)


parse.pvals.qq <- function(pvector,lim=7) {

  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )

  o[o>lim] <- lim

  out <- list(o=o,e=e)

  return(out)
}



write.table(a.indel[wanted,],file="OPRM1.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

library(plyr)

options(stringsAsFactors=FALSE)    
setwd("~/Documents/Data/trial_GWAS/run1/assoc_results")

SNP = read.table("LM-cases-assoc.csv", header=T,  na.strings = "NA")  # skip importing the the #eigenvalue labels 

head(SNP)
str(SNP)

#----------my script------------------------------------
#  the expected dist of P values is NOT normal but an even distribution from 0-1
#  the expected dist of STAT is TDist which seems to == rnorm in this case anyway (see below)
#  dof for t tests is n1+n2 -2 and thus for its distribution too

# SNP = subset(file1, TEST=="ADD")


#  check normality of data
plot(sort(SNP$P[20000:70000]), cex=0.2, col="grey")
lines(sort(seq(0,1, length.out=50000)), col="red", lwd=2, lty=2)

plot(sort(SNP$STAT[40000:70000]), col="grey")
lines(sort(rt(30000, df=99998)), col="red", lwd=4, lty=2)
lines(sort(rt(30000, df=120)), col="green", lwd=4, lty=2)
lines(sort(rt(30000, df=12)), col="cyan", lwd=4, lty=2)

# ---  my script for T STATISTIC  
#  i think negative t scores are as important as positive ones, they just have mean1- mean2 backwards

dof = (length(SNP$STAT) *2) -2
t = sort(-log10(abs(rt(length(SNP$STAT), df=dof))))
# you can't log a negative number. the distribution is symmetrical around 0 anyway, so abs() or chop it off here
test = sort(-log10(abs(SNP$STAT)))
#  test2= qqnorm(SNP$STAT[1:20000], col="forestgreen")

par(bg="white")
jpeg("run2-STAT-myscript.jpg", width=800, height=600)
plot(t, test, main="QQ of run2 (no controls) -log STAT (dof=540k)", col="forestgreen", xlim=c(0,7), ylim=c(0,7), ylab="-log(STAT)", xlab="-log(tDist)")
abline(coef=c(0,1), col=1, lwd=2)
dev.off()


# -------------  my P value QQ

#  using a vector to generate expected (adrian's method) appears better than runif
#runif = -log10(sort(runif(length(SNP$P), min=0, max=1)))    #  other alternative is ppoints
vector = -log10(1:length(SNP$P)/length(SNP$P))
	
test.p = -log10(sort((SNP$P[1:length(SNP$P)]), decreasing=F))

par(bg="white")
jpeg("QQ-run2-myscript.jpg", width=800, height=600)
#plot(runif, test.p, main="QQ of log p", col="forestgreen", xlim=c(0,6), ylim=c(0,6), ylab="-log(P)", xlab="-log(sequence)")
plot (vector, test.p, main="QQ of run2 no controls",col="blue", type="p", cex=1, pch=20)
abline(coef=c(0,1), col=1, lwd=1)
# legend("bottomright", legend=c("runif","vector","seq","rnorm"), col=c("forestgreen","navy","deeppink","purple"), lty=1, title="gradient", cex=0.8, lwd=2) 
dev.off()


# -------------------   broad  institute script modified for STAT


observed <- sort(abs(SNP$STAT), decreasing=T)
lobs <- -(log10(observed))

dof = (length(SNP$STAT) *2) -2
lexp = sort(-log10(abs(rt(length(SNP$STAT), df=dof))))

png("QQ-run2-broad-STAT.png", units="px", width=800, height=600)
plot(c(-5,7), c(-5,7), main="QQ no controls run2",col="red", lwd=3, type="l", xlab="Expected (-logSTAT)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=20, cex=.8, bg="deepskyblue") 
dev.off()


# -------------------  RAW   broad script
	
observed <- sort(SNP$P)

lobs <- -(log10(observed))

expected <- c(1:length(observed)) 

lexp <- -(log10(expected / (length(expected)+1)))
	
plot(-log10(expected/(length(expected)+1)), type="l")
lines(lexp, col="red", add=T)	
	
	
png("QQ_run2-broad-P.png", units="px", width=800, height=600)
plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=20, cex=.8, bg="grey") 
dev.off()




#   adrian's script   --------------------------------

	
parse.pvals.qq <- function(pvector,lim=7) {

  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )

  o[o>lim] <- lim

  out <- list(o=o,e=e)

  return(out)
}

lim = 9
qq <- parse.pvals.qq(SNP[,"P"],lim=lim)
 
ylab = expression(Observed~~-log[10](italic(p)))
xlab = expression(Expected~~-log[10](italic(p)))
 
plot(qq$e,
     qq$o,
     xlim=c(0,lim),ylim=c(0,lim),
     pch=20,col='deepskyblue',
     xlab=xlab,ylab=ylab, main="title")



	abline(coef=c(0,1), col=1, lwd=2)   # my abline
















































########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

### list of snps use this on command line :: 74th is the column with the rs ids...
## head -50 chr10.AOGC-NGS.2013.pointwise.txt | cut -f 74
## /media/ga-apps/UQCCG/Programming/scripts/PerlScripts/GrepMafFiles.pl 35SNPs_or_proxies.csv 74 chr10.AOGC-NGS.2013.pointwise.txt
## the.samples<-sample.sheet.full[,"ParticipantCode"]
## the.samples
## indels<-indels[,paste(the.samples,"GT",sep=".")]
# Asian-AML   Control EXOME-AML 

