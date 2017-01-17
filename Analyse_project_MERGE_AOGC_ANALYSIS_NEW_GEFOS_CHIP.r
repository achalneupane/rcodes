
######## ONLY NEED TO CHOOSE A DIRECTORY AND EXTENSIONS - used tab delimited files 

options(width=250,max.print=5000)

code.dir<-"/media/Bioinform-D/Research/AML sequencing"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")

###############################################
###### THIS IS the super annotion run USE ALL THE FILTER AND NOVEL DATABASES
##Build AOGC

##################### if have a genotype component
genotype.file.location<-"/media/UQCCG-Analysis/AOGC_exome_chip/working_genotypes"
genotype.file.prefix<-"recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL"
#########################################################

annotate.dir<-"/media/scratch2/AOGC-NGS/Annotate" # dir(annotate.dir)
analysis.dir<-"/media/scratch2/AOGC-NGS/Analysis" # dir(analysis.dir)
project.extension<-".analysis.txt" ## justX1 the exterion not fam.extension!
project.name<-"AOGC-Genotyping.output"
fam<-c("AOGC_ALL") #  ALL or  c() ""-one project (the prefix of the summary files to collect
the.sample.sheet<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_SAMPLES_PHENOTYPES_Nov.1.2013_RESIDUALS.txt"

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
library(MultiPhen)
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
gene.list.file<-"/media/Bioinform-D/Research/exome Chip/AOGC_Gene_Hits_0.05_0.2-sent.txt"
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




gerp.score.threshold.high<-2.0 # gerp score >= will be included
gerp.score.threshold.low<-1.5 # gerp score >= will be included
gerp.score.threshold.unknown<-0


#generic.filter.DB
maf.threshold.filter.to.use<-sort(maf.threshold.filter.to.use)

 maf.threshold.filter<-maf.threshold.filter.to.use
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
load("/media/Bioinform-D/Research/AML sequencing/Human_Exome_Targets_illumina_v2_hg19_targets.RData")
ill.gr<-data.gr
load("/media/Bioinform-D/Research/AML sequencing/Human_Exome_Targets_Nimble_v2_hg19_targets.RData")
nim.gr<-data.gr

genome(ill.gr)<-"hg19"
genome(nim.gr)<-"hg19"

ill.gr<-ill.gr+200
nim.gr<-nim.gr+200

overlaps<-overlapsAny(ill.gr,nim.gr)
possible.loci.ori<-ill.gr[overlaps]

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




sample.sheet.full<-read.delim(the.sample.sheet,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
sample.sheet.full[1:5,]
colnames(sample.sheet.full)


##### fix 0 and 9 for missing to NA

## pheno.types<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_RAD","BMD_EFF_STD_LS","BMD_EFF_STD_FN","EVER_FX_50_EXCL_TRIVIAL")
## names(pheno.types)<-c("HIP","RAD","LS","FN","FX")

pheno.types<-c("BMD_EFF_STD_HIP")
 names(pheno.types)<-c("HIP")

case.control<-c("EVER_FX_50_EXCL_TRIVIAL")
# ib<-1
for(ib in 1:length(case.control)){
  if(!(case.control[ib] %in% colnames(sample.sheet.full))){next}
  sample.sheet.full[(is.na(sample.sheet.full[,case.control[ib]]) | sample.sheet.full[,case.control[ib]]==0 | sample.sheet.full[,case.control[ib]]==9)  ,case.control[ib]]<-NA
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

#################### Read plink file
if(grepl("^chr",the.chr)){plink.chr<-gsub("^chr","",the.chr)}else{
  plink.chr<-the.chr}

if(plink.chr=="X"){plink.chr<-23}
if(plink.chr=="Y"){plink.chr<-24}
if(plink.chr=="XY"){plink.chr<-25}
if(plink.chr=="M"){plink.chr<-26}
plink.chr



plink.file<-paste(genotype.file.prefix,"_chr",plink.chr,sep="")


plink.file
g.indel<-read.plink(paste(genotype.file.location,plink.file,sep="/"))

## test<-read.table("plink.frq",header=T)
## ori<-read.table("ori.frq",header=T)

dim(a.indel)
dim(g.indel)
target<-"37442658"
posn1<-match(target,a.indel[,"end"])
posn2<-match(target,g.indel[,"start"])
posn1
posn2

a.indel[posn1,1:10]
g.indel[posn2,1:10]






#tapply(gene.list[,"CHR"],gene.list[,"CHR"],length)
all.possible.samples<-gsub(".GT$","",colnames(a.indel)[grep(".GT$",colnames(a.indel))],perl=TRUE)
length(all.possible.samples)
pheno.types





############################### do one phenotype #################### 
# ipheno<-1
for(ipheno in 1:length(pheno.types)){
  
print(paste("Doing phenotype",pheno.types[ipheno]))


target.pheno<-names(pheno.types)[ipheno]
target.pheno.col<-pheno.types[ipheno]


pheno<-sample.sheet.full[!is.na(sample.sheet.full[,target.pheno.col]) ,] ## pheno only contains SAMPLES that have a phenotype
print(dim(pheno))
print(paste("Number Samples:",dim(pheno)[1]))
covars<-c("PCA1","PCA2","PCA3","PCA4")
dim(pheno)
pheno[1:5,]


formula<-paste(target.pheno.col,"~",paste(covars,collapse="+"),sep="")
print(formula)

the.samples<-paste(pheno[,"SAMPLE"],"GT",sep=".")  ## samples same order as in pheno
print(paste("Number samples: ",length(the.samples),sep=""))


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






######################################## ADDITIONAL QC #################################
##########################################################################
## ipheno<-2
## target.pheno<-names(pheno.types)[ipheno]
## target.pheno.col<-pheno.types[ipheno]
## target.pheno.col<-"BMD_EFF_STD_HIP"
## target.pheno.col

## the.low<-as.numeric(pheno[,target.pheno.col])<0
## the.high<-as.numeric(pheno[,target.pheno.col])>0
## tapply(pheno[the.low,"AFFSTAT_IN_WORDS"],pheno[the.low,"AFFSTAT_IN_WORDS"],length)
## tapply(pheno[the.high,"AFFSTAT_IN_WORDS"],pheno[the.high,"AFFSTAT_IN_WORDS"],length)
## pheno[1:5,"AFFSTAT_IN_WORDS"]
## tapply
## plot(pheno[,target.pheno.col])
## text(c(1:dim(pheno)[1]),pheno[,target.pheno.col],labels=pheno[,"AFFSTAT_IN_WORDS"])

##### AFFSTAT_IN_WORDS is mostly correct to the regression #####
pheno[1:5,]


summary.geno.extra<-{}

target<-"LOW"
use.samples<-the.samples[pheno[,"AFFSTAT_IN_WORDS"]=="low"]
length(use.samples)
genotypes<-a.indel[,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")
#summary.geno[1:5,]
if(is.null(dim(summary.geno.extra))){
  summary.geno.extra<-summary.geno
}else{
  summary.geno.extra<-cbind(summary.geno.extra,summary.geno)
}



target<-"HIGH"
use.samples<-the.samples[pheno[,"AFFSTAT_IN_WORDS"]=="high"]
length(use.samples)
genotypes<-a.indel[,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")
#summary.geno[1:5,]
if(is.null(dim(summary.geno.extra))){
  summary.geno.extra<-summary.geno
}else{
  summary.geno.extra<-cbind(summary.geno.extra,summary.geno)
}

target<-"LOW.pheno"
if(target.pheno.col %in% case.control){
  use.samples<-the.samples[as.numeric(pheno[,target.pheno.col])==1]
                           }else{
                   use.samples<-the.samples[as.numeric(pheno[,target.pheno.col])<=0]
                 }
length(use.samples)
genotypes<-a.indel[,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")
#summary.geno[1:5,]
if(is.null(dim(summary.geno.extra))){
  summary.geno.extra<-summary.geno
}else{
  summary.geno.extra<-cbind(summary.geno.extra,summary.geno)
}

target<-"HIGH.pheno"
if(target.pheno.col %in% case.control){
  use.samples<-the.samples[as.numeric(pheno[,target.pheno.col])==2]
                    }else{
  use.samples<-the.samples[as.numeric(pheno[,target.pheno.col])>0]
      }
length(use.samples)
genotypes<-a.indel[,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")
#summary.geno[1:5,]
if(is.null(dim(summary.geno.extra))){
  summary.geno.extra<-summary.geno
}else{
  summary.geno.extra<-cbind(summary.geno.extra,summary.geno)
}


target<-"nimblegen"
use.samples<-the.samples[the.samples %in% nim.samples]
 length(use.samples)
genotypes<-a.indel[,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")
#summary.geno[1:5,]
if(is.null(dim(summary.geno.extra))){
  summary.geno.extra<-summary.geno
}else{
  summary.geno.extra<-cbind(summary.geno.extra,summary.geno)
}

target<-"illumina"
use.samples<-the.samples[the.samples %in% ill.samples]
 length(use.samples)
genotypes<-a.indel[,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")
#summary.geno[1:5,]
if(is.null(dim(summary.geno.extra))){
  summary.geno.extra<-summary.geno
}else{
  summary.geno.extra<-cbind(summary.geno.extra,summary.geno)
}











 summary.geno.extra[1:5,]
rownames(summary.geno.extra)<-key

   ## class(summary.geno.extra)

######################## inbreeding and no genotypes

group.maf.thresh<-0.20 # group.maf.thresh<-0.10  if more common that this in Group then discard: Discard Rare but common in this cohort
missing.threshold<-0.60 # missing.threshold<-0.50  60 % of genotypes missing
missing.threshold.nimblgen<-0.20
missing.threshold.illumina<-0.20 

## rare.in.group<-( summary.geno.extra[,c("MAF.HIGH","MAF.LOW")]< group.maf.thresh) | ( summary.geno.extra[,c("MAF.HIGH","MAF.LOW")] > (1-group.maf.thresh)) | (is.na( summary.geno.extra[,c("MAF.HIGH","MAF.LOW")]))
rare.in.group<-( summary.geno.extra[,c("MAF.HIGH","MAF.LOW","MAF.HIGH.pheno","MAF.LOW.pheno")]< group.maf.thresh) | ( summary.geno.extra[,c("MAF.HIGH","MAF.LOW","MAF.HIGH.pheno","MAF.LOW.pheno")] > (1-group.maf.thresh)) | (is.na( summary.geno.extra[,c("MAF.HIGH","MAF.LOW","MAF.HIGH.pheno","MAF.LOW.pheno")]))

#rare.in.group<-cbind( indels[,c("MAF.HIGH","MAF.LOW")]< group.maf.thresh) , (indels[,c("MAF.HIGH","MAF.LOW")] > (1-group.maf.thresh)) , (is.na(indels[,c("MAF.HIGH","MAF.LOW")])))
rare.in.group<-combine.boolean(rare.in.group,colnames(rare.in.group),"OR")
sum(!rare.in.group)



#no.genotypes<-(indels[,c("MAF.HIGH","MAF.LOW")]== 0)  | (is.na(indels[,c("MAF.HIGH","MAF.LOW")])) # no genotypes in test classes for a mutataion after individaul quality filtering
no.genotypes<-(summary.geno.extra[,c("MAF.HIGH","MAF.LOW","MAF.HIGH.pheno","MAF.LOW.pheno")]== 0)  | (is.na( summary.geno.extra[,c("MAF.HIGH","MAF.LOW","MAF.HIGH.pheno","MAF.LOW.pheno")])) # no genotypes in test classes for a mutataion after individaul quality filtering
no.genotypes[1:5,]
no.genotypes<-combine.boolean(no.genotypes,colnames(no.genotypes),"AND")
no.genotypes[1:5]
sum(no.genotypes)

summary.geno.extra[1:5,]
missing.targets<-c("LOW","HIGH","LOW.pheno","HIGH.pheno","nimblegen","illumina")
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
colnames(high.missing)<-missing.targets
rownames(high.missing)<-key
                 
## high.missing<- cbind(as.numeric(summary.geno.extra[,"MISSING.Alleles.LOW"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.LOW"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.LOW"])),
##                      as.numeric(summary.geno.extra[,"MISSING.Alleles.HIGH"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.HIGH"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.HIGH"])),
##                      as.numeric( summary.geno.extra[,"MISSING.Alleles.LOW.pheno"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.LOW.pheno"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.LOW.pheno"])),
##                      as.numeric(summary.geno.extra[,"MISSING.Alleles.HIGH.pheno"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.HIGH.pheno"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.HIGH.pheno"])),
##                      as.numeric(summary.geno.extra[,"MISSING.Alleles.nimblegen"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.nimblegen"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.nimblegen"])),
##                      as.numeric(summary.geno.extra[,"MISSING.Alleles.illumina"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.illumina"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.illumina"]))

                     
##                      )
colnames(high.missing)
high.total.missing<-subset(high.missing,select=c("LOW","HIGH","LOW.pheno","HIGH.pheno"))
high.total.missing[1:5,]
high.total.missing<-high.total.missing > missing.threshold
high.total.missing<-combine.boolean(high.total.missing,c("LOW","HIGH","LOW.pheno","HIGH.pheno"),"OR")
sum(high.total.missing)


nimblegen.total.missing<-subset(high.missing,select=c("nimblegen"))
nimblegen.total.missing[1:5,]
nimblegen.total.missing<-nimblegen.total.missing > missing.threshold.nimblgen
## nimblegen.total.missing<-combine.boolean(high.total.missing,c("LOW","HIGH","LOW.pheno","HIGH.pheno"),"OR")
sum(nimblegen.total.missing)


illumina.total.missing<-subset(high.missing,select=c("illumina"))
illumina.total.missing[1:5,]
illumina.total.missing<-illumina.total.missing > missing.threshold.illumina
## nimblegen.total.missing<-combine.boolean(high.total.missing,c("LOW","HIGH","LOW.pheno","HIGH.pheno"),"OR")
sum(illumina.total.missing)

sum(high.total.missing | nimblegen.total.missing | illumina.total.missing)
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

my.loci<-IRanges(start=as.numeric(a.indel[,"start"]),end=as.numeric(a.indel[,"end"]))
on.chr<-as.vector(seqnames(possible.loci.ori))==a.indel[1,"chr"] # paste("chr",a.indel[1,"chr"],sep="")
possible.space<-IRanges(start=as.numeric(start(possible.loci.ori))[on.chr],end=as.numeric(end(possible.loci.ori))[on.chr])
common.loci<-overlapsAny(my.loci,possible.space)
sum(common.loci)
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
 if(sum( !(global.labs %in% c(colnames(a.indel),colnames(summary.geno.extra)))  )>0){print(paste("WARNING postion filters missing for data",
          toString(global.labs[!(global.labs %in% c(colnames(a.indel),colnames(filter.table),colnames(filter.table.pholy))) ])))}
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



# bad.coding<-test.for.coding.type(a.indel,geneanno.DB,c("stopgain SNV","stoploss SNV","frameshift deletion","frameshift insertion"))
#bad.frame<-test.for.coding.type(geneanno.table,geneanno.DB,c("frameshift deletion","frameshift insertion"))


basic.qual<-combine.boolean(qual,c("QUAL", "QD", "HRun", "SB"),"AND")
gatk.qual<-combine.boolean(qual,c("FILTER_PASS", "FILTER_100" ),"OR")
full.qual<-combine.boolean(cbind(basic.qual,gatk.qual),"all","OR")


basic.qual<-combine.boolean(qual,c("QUAL", "QD", "HRun", "SB","FILTER_100"),"AND")
gatk.qual<-qual[,"FILTER_PASS"]
full.qual<-combine.boolean(cbind(basic.qual,gatk.qual),"all","OR")
sum(full.qual)

## full.qual<-gatk.qual
## sum(full.qual)
#full.qual<-qual[,"FILTER_PASS"]
sum(full.qual)
qual[full.qual,c("FILTER_PASS", "FILTER_100")][1:2,]


#any.functional<-combine.boolean(qual,c("PolyPhen.low","SIFT.high","mut.taster.high","phylo.high","PolyPhen.bad","SIFT.bad","GERP.high","ljb_gerp.high"),"OR")
functional<-combine.boolean(qual,c("PolyPhen.low","SIFT.high","PolyPhen.bad","SIFT.bad","GERP.high"),"OR")
sum(functional)
functional<-functional | bad.coding # include really bad protein changes
sum(functional)


maf.lt.all[1:5,]
maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.5"])

#pass<- rare.in.group & !no.genotypes & !high.missing & common.loci

##########################################################################
##########################################################################
##########################################################################
not.flat.genotype<-!qual[,"flat"]
sum(not.flat.genotype)
a.indel[!not.flat.genotype,"TYPE"]

pass<-full.qual & functional & maf.filter & rare.in.group & !no.genotypes  & not.flat.genotype
sum(pass)

## pass<-full.qual & functional & maf.filter & rare.in.group & !no.genotypes  & not.flat.genotype & !(high.total.missing | nimblegen.total.missing | illumina.total.missing)
## sum(pass)

## sum(pass & not.flat.genotype)
## a.indel[(pass & not.flat.genotype),"TYPE"]
##########################################################################
##########################################################################
##########################################################################
## filtered.genotype.summary()
## full.genotype.summary()

################################# GEFOS FILTERING cause sending all
#pass<-pass[the.snps] ### GEOFS 
genotypes<-a.indel[pass,the.samples] ## ordered correctly for phenotypes
snp.names<-key[pass] ## GEFOS ony name with start


dim(genotypes)

print("start QC")
genotypes[genotypes=="NA"]<-NA
genotypes[genotypes=="0/0"]<-0
genotypes[genotypes=="0/1"]<-1
genotypes[genotypes=="1/1"]<-2

## genotypes[to.flip,][1:15,1:10]
## genotypes[!to.flip,][1:5,1:10]


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


########## make snpInfo
#tapply(gene.list[,"CHR"],gene.list[,"CHR"],length)
## gene.list[gene.list[,"CHR"]=="Y",]
## a.indel[the.snps,1:20][1:12,]

#a.indel[pass,][1:5,1:10]
gene.list[gene.list[,"CHR"]==the.chr,][1:5,]

the.chr<-a.indel[1,"chr"]
my.snps<-IRanges(start=as.numeric(a.indel[pass,"start"]),end=as.numeric(a.indel[pass,"end"]))
# a.indel[the.snps,][,1]

gene.on.chr<-gene.list[gene.list[,"CHR"]==the.chr,]



gene.space<-IRanges(start=as.numeric(gene.on.chr[,"START"]),end=as.numeric(gene.on.chr[,"END"]))
gene.space
my.snps

overlaps<-findOverlaps(gene.space,my.snps)


snpinfo<-cbind(snp.names[subjectHits(overlaps)],gene.on.chr[queryHits(overlaps),"GENE_NAME"])
## snpinfo<-cbind(snp.names,a.indel[pass,"Gene.Names"])
snpinfo[1:5,]



colnames(snpinfo)<-c("Name","gene")
dim(snpinfo)
 snpinfo[1:5,]
dim(snpinfo)

colnames(genotypes)<-snp.names
rownames(genotypes)<-gsub(".GT$","",the.samples)

####### clean up missing
## overlaps
## snpinfo[1:5,]
 summary.geno.extra[1:5,]
 summary.geno.extra[snpinfo[,"Name"],][1:5,]


missing.snps<-!(snp.names %in% snpinfo[,"Name"])
 sum(missing.snps)
sum(!missing.snps)
## length(missing.snps)
## snp.names[missing.snps][1:5]


genotypes<-genotypes[,!missing.snps]


## snpinfo[1:5,]
## genotypes[1:5,1:20]
## gene.on.chr[1:5,]
## dim(genotypes)
## dim(pheno)
## dim(snpinfo)
## formula



#match(rownames(genotypes),pheno[,"SAMPLE"])
## ann<-toString(colnames(a.indel)[1:70])


#key<-build.key(a.indel[the.snps,],core.ann)
#key<-build.key(a.indel[the.snps,],"start")
posns<-match(snpinfo[,"Name"],snp.names)
missing<-is.na(posns)
print(sum(missing))

ann<-a.indel[pass,][posns,c(core.ann,"chr","start","end","REF","ALT","TYPE","refGene::location","refGene::type","refGene::gene","knownGene::location","knownGene::type","knownGene::gene","ensGene::location","ensGene::type","ensGene::gene","Gene.Names","description","gene_biotype","OMIM (Gene::Status::OMIM::description::disease)","Consequence.Embl","Uploaded_variation.Embl","Gene.Embl","Feature.Embl","Protein_position.Embl","Amino_acids.Embl","FILTER","wanted.muts","wanted.muts.coding","MAF.lt:0","MAF.lt:0.001","MAF.lt:0.005","MAF.lt:0.01","MAF.lt:0.025","MAF.lt:0.5","gerp.scores","MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO","MAF.ALL","ALT.Alleles.ALL","REF.Alleles.ALL","TOTAL.Alleles.ALL","MISSING.Alleles.ALL","ALT_HOMO.ALL","ALT_HETRO.ALL","GENO.ALL","MAF.OTHER","ALT.Alleles.OTHER","REF.Alleles.OTHER","TOTAL.Alleles.OTHER","MISSING.Alleles.OTHER","ALT_HOMO.OTHER","ALT_HETRO.OTHER","GENO.OTHER","Hetero.ALT.reads","Hetero.REF.reads","Hetero.Read.Balance","PolyPhen.desc","PolyPhen.scores","SIFT.desc","SIFT.scores","Hetero.ALT.reads.ALL","Hetero.REF.reads.ALL")]


summary.geno.extra.out<-summary.geno.extra[pass,][posns,]
summary.geno.extra.out[1:5,]

high.missing.out<-high.missing[pass,][posns,]
high.missing.out[1:5,]


dim(summary.geno.extra.out)
dim(ann)

## ann[1:5,]
## snpinfo[1:5,]

print(paste("Running skatCohort With Formula ", formula, sep=""))
cohort.seq <- skatCohort(genotypes, formula, SNPInfo = snpinfo, data=pheno,verbose=FALSE)


## skatMeta(..., SNPInfo=NULL, wts = function(maf){dbeta(maf,1,25)}, method = "saddlepoint", 
##              snpNames = "Name", aggregateBy = "gene", mafRange = c(0,0.5), verbose=FALSE)


meta.results.skat<-skatMeta(cohort.seq,SNPInfo = snpinfo)
meta.results.burden<-burdenMeta(cohort.seq,wts=1,mafRange = c(0,1),SNPInfo = snpinfo)
## meta.results.skatO<-skatOMeta(cohort.seq,burden.wts =1,SNPInfo = snpinfo)
meta.results.skatO<-meta.results.burden


###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
## test<-c("MYCBP2","SLC25A24","TMCO3","C13orf35")
## test<-c("MYCBP2","SLC25A24","TMCO3","C13orf35")
## test<-c("chr22:41252435-41252687:ST13")
## test<-c("chr4:125631155-125631932:ANKRD50","chr4:88416133-88416279:SPARCL1","chr4:120085366-120085549:MYOZ2")

## meta.results.skat[meta.results.skat[,"gene"] %in% test,]
## meta.results.burden[meta.results.burden[,"gene"] %in% test,]

## loci<-snpinfo[snpinfo[,"gene"] %in% test,"Name"]
## #summary.geno.extra.out[loci,]
## high.missing.out[loci,]
## qual[loci,]
## a.indel[loci,extra]
## snpinfo[1:5,]
## qual[1:5,c("FILTER_PASS", "FILTER_100" )]

## tapply(a.indel[pass,"Consequence.Embl"],a.indel[pass,"Consequence.Embl"],length)
## ##  dbeta(x, shape1, shape2, ncp = 0, log = FALSE)
## ## shape1, shape2: positive parameters of the Beta distribution.
## ##  ncp: non-centrality parameter.
## extra<-c("MAF.lt:0.5","Consequence.Embl","wanted.muts","wanted.muts.coding",unique(global.quality.labs)[unique(global.quality.labs) %in% colnames(a.indel)],"Hetero.ALT.reads","Hetero.REF.reads","Hetero.Read.Balance","culprit")
## ## Warning message:
## ## In check_format_skat(Z, SNPInfo, nullmodel, aggregateBy, snpNames) :
## ##   Some missing genotypes - will be imputed to average dose
## exclude.samples<-gsub(".GT$","",nim.samples)

## dim(genotypes)
## dim(pheno)
## exclude<-rownames(genotypes) %in% exclude.samples
## sum(exclude)
## genotypes.ex<-genotypes[!exclude,]
## pheno.ex<-pheno[!exclude,]
## dim(genotypes.ex)
## dim(pheno.ex)

## cohort.seq.ex <- skatCohort(genotypes.ex, formula, SNPInfo = snpinfo, data=pheno.ex,verbose=FALSE)
## meta.results.skat.ex<-skatMeta(cohort.seq,SNPInfo = snpinfo)
## meta.results.burden.ex<-burdenMeta(cohort.seq.ex,wts=1,mafRange = c(0,1),SNPInfo = snpinfo)
## #meta.results.skatO.ex<-skatOMeta(cohort.seq.ex,burden.wts =1,SNPInfo = snpinfo)


## meta.results.skat.ex[meta.results.skat.ex[,"gene"] %in% test,]
## meta.results.burden.ex[meta.results.burden.ex[,"gene"] %in% test,]



###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$













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


###############################



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

