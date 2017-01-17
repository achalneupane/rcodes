
#### red in per patiant shared regions....exampls
data<-cbind(c("chr1","chr1","chr1"),c(10,50,50),c(30,70,100))
data2<-cbind(c("chr1","chr1","chr1"),c(10,60,10),c(30,70,40))



Views on a 249250621-length Rle subject

views:
    start end width
[1]    10  30    21 [10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10]
[2]    60  70    11 [10 10 10 10 10 10 10 10 10 10 10]

overlapping<-cbind(start(x[["chr1"]]),end(x[["chr1"]]))
colnames(overlapping)<-c("starts","ends")
overlapping

write.table(sample.sheet.full,file="sampe.sheet.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


cyto<-read.delim("/media/TRI-T-DRIVE-tpleo/uqdi/Core_Services/UQCCG/Sequencing/CompleteGenomics/MergedAML_Cytogenetics_etc.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
cg.done<-read.delim("/media/UQCCG/Sequencing/CompleteGenomics/917_Data_Delivery_091214.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
cyto[1:5,]
cg.done[1:5,]


posns<-match(cg.done[,"Cust..Sample.ID"],cyto[,"ID"])
posns
missing<-is.na(posns)
dim(cg.done)
sum(!missing)

cg.done[missing,"Cust..Sample.ID"]


cbind(cg.done[missing,"Cust..Sample.ID"],cg.done[missing,"Cust..Sample.ID"])

cyto.have<-
  cyto[posns[!missing],]

compare<-cbind(cg.done[!missing,],cyto[posns[!missing],])
getwd()
setwd("/media/UQCCG/Sequencing/CompleteGenomics")
write.table(compare,file="cg.cyto2.csv",col.names=TRUE,row.names=FALSE,sep="\t")



#######################################

tcga<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014/2014-11-04_AML_TCGA_Replication/Analysis/Phenotypes_TCGA.csv",header=T,sep="\t",fill=TRUE,as.is=TRUE,stringsAsFactors=FALSE)
tcga[1:5,1:11]
colnames(tcga)

tail(sample.sheet.full)
cases<-sample.sheet.full[sample.sheet.full$AffectionStatus==2,"ParticipantCode"]

cases.short<-unlist(lapply(strsplit(cases,split="-"),function(x) x[3]))
cases<-cbind(cases,cases.short)
cases[1:5,]







posns<-match(cases[,"cases.short"],tcga[,"TCGA.Patient.ID"])
posns
missing<-is.na(posns)
sum(missing)
cases[missing,] # TCGA-AB-2852-03A-01W-0726-08

extra<-tcga[posns,c("TCGA.Patient.ID","Sex","Race","Age","FAB","X.BM.Blast")]
dim(extra)
dim(caeses)

cases[1:5,]
extra[1:5,]
cases<-cbind(cases,extra)

cases[is.na(cases[,"TCGA.Patient.ID"]),]
cases[is.na(cases[,"TCGA.Patient.ID"]),"TCGA.Patient.ID"]<-2852
cases[is.na(cases[,"TCGA.Patient.ID"]),"Race"]<-"unknown"


tail(sample.sheet.full)
posns<-match(sample.sheet.full[,"ParticipantCode"],cases[,"cases"])
posns
missing<-is.na(posns)
sum(missing)
add.extra<-extra[posns,]
sample.sheet.full<-cbind(sample.sheet.full,add.extra)

a.sheet<-cbind("HIP","ALL",all.possible.samples,0,0,0,1)
a.sheet[1:5,]
getwd()

write.table(a.sheet,file="sample_sheet.csv",col.names=TRUE,row.names=FALSE,sep=",")


######## ONLY NEED TO CHOOSE A DIRECTORY AND EXTENSIONS - used tab delimited files 
#source("http://bioconductor.org/biocLite.R")
# biocLite(c("HardyWeinberg"))n
# install.packages("HardyWeinberg")
###############################################
#analysis.dir<-"/media/ga-apps/UQCCG/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/Analysis"
#annotate.dir<-"/media/ga-apps/UQCCG/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/Annotate"

analysis.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis"
annotate.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/Annotate"

project.extension<-".All-maf-filtered.txt"
project.name<-"2013-02-27_AML_with_AOGCControl" ## prefix for output file
fam<-c("TGCM-AML") #  ALL or  c() ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/BAM/TGCM-AML-combine_SampleSheet.csv"
the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/BAM/TGCM-AML-combine_SampleSheet.csv"

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

###############################################
#analysis.dir<-"/media/ga-apps/UQCCG/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/Analysis"

#annotate.dir<-"/media/ga-apps/UQCCG/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/Annotate"

analysis.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-06-24_AML_RSGB_AOGC_withHaplotypeCaller/Analysis"
annotate.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-10-27_AML_with_AOGCControl_NoFailedLane/Annotate"

project.extension<-".wanted.All-maf-filtered.txt"
project.name<-"2014-06-24_AML_RSGB_AOGC_withHaplotypeCaller" ## prefix for output file
fam<-c("ALL") #  ALL or  c() ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/BAM/TGCM-AML-combine_SampleSheet.csv"
the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-06-24_AML_RSGB_AOGC_withHaplotypeCaller/BAM/TGCM-AML_RSGB_PILOT_SampleSheet_controls.csv"

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

###############################################  NEW HC u  with all data
#analysis.dir<-"/media/ga-apps/UQCCG/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/Analysis"

#annotate.dir<-"/media/ga-apps/UQCCG/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/Annotate"
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-10-09_AML_CompleteGenomics_HC/Analysis/2014-10-09_AML_CompleteGenomics_HC.chrALL..ALL.ALL_GENOTYPES_.analysis-maf-filtered.txt
analysis.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014/2014-11-04_AML_TCGA_Replication/Analysis"
annotate.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014/2014-11-04_AML_TCGA_Replication/Annotate"

## /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-10-09_AML_CompleteGenomics_HC/Analysis/2014-10-09_AML_CompleteGenomics_HC.chr1.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
##/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-10-09_AML_CompleteGenomics_HC/Analysis/2014-10-09_AML_CompleteGenomics_HC.chrALL.Indel.ALL.ALL_GENOTYPES_analysis.txt

## project.extension<-".analysis-maf-filtered.txt"
## project.name<-"2014-10-09_AML_CompleteGenomics_HC." ## prefix for output file
## fam<-c("chrALL_GENOTYPES") #  ALL or  c() ""-one project (the prefix of the summary files to collect
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-11-04_AML_TCGA_Replication/Analysis/2014-11-04_AML_TCGA_Replication.chrALL.ACC_good_qual.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt

project.extension<-"_analysis-maf-filtered.txt"
project.name<-"2014-10-09_AML_CompleteGenomics_HC." ## prefix for output file
fam<-c("chrALL.ACC_good_qual.ALL.ALL_GENOTYPES") #  ALL or  c() ""-one project (the prefix of the summary files to collect

#the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/BAM/TGCM-AML-combine_SampleSheet.csv"
the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014/2014-11-04_AML_TCGA_Replication/Analysis/Full.sample_sheet.wP.csv"

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



contaminated.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/contaminated_AOGC_SEQ_samples.txt"
contaminated<-read.table(contaminated.file,header=F,fill=TRUE,sep="\t",stringsAsFactors=FALSE)



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


############################################# POPULATION MAF FILTER - PART A
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER

maf.threshold<-0.0  #  MAF threshold for annovar calling zero useful to get back all results !!do not modify!!
maf.threshold.filter.to.use<-c(0.001,0.01,0.05)
maf.threshold.filter.to.use<-sort(as.numeric(maf.threshold.filter.to.use))

filter.cols.novel.use<-c("NHBLI_6500_ANNOVAR_ALL","NHBLI_6500_ALL","NHLBI_5400_ALL","NHLBI_5400_EUR","NHLBI_5400_AFR","1000genome","1000genome_asian","1000genome_mine","snp141","snp141_clinical","snp137","CG69","EUR_ASN_AFR_INDEL","AOGC-NGS_ALL","AOGC-NGS_ALL_OLD","Chinese") ##
filter.cols.maf.use<-c("PopFreqMax","NHBLI_6500_ANNOVAR_ALL","NHBLI_6500_ALL","NHBLI_6500_EA","NHBLI_6500_AA","NHLBI_5400_ALL","1000genome","snp141","snp137","snp135")
 maf.threshold.filter<-maf.threshold.filter.to.use

############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER



##################################################### DEFINE A GENE LIST  #####################################################
##################################################### DEFINE A GENE LIST  #####################################################
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


seq.type.file<-"/media/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_SAMPLE_Tue_Oct_14_2014.txt"
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

library("BSgenome.Hsapiens.UCSC.hg19")
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
sample.sheet.full[1:10,]
colnames(sample.sheet.full)
dim(sample.sheet.full)

## names<-colnames(sample.sheet.full.rep)[ colnames(sample.sheet.full.rep) %in% colnames(sample.sheet.full)]


## names
## sample.sheet.full.1<-sample.sheet.full[,names]
## sample.sheet.full.2<-sample.sheet.full.rep[,names]
## sample.sheet.full<-rbind(sample.sheet.full.2,sample.sheet.full.1)

##### fix 0 and 9 for missing to NA

## pheno.types<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_RAD","BMD_EFF_STD_LS","BMD_EFF_STD_FN","EVER_FX_50_EXCL_TRIVIAL")
## names(pheno.types)<-c("HIP","RAD","LS","FN","FX")


######### Check and fix the sample sheet

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

################################################################### AOGC STATS
the.sample.sheet.aogc<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"


sum(sample.sheet.full[,"ParticipantCode"] %in% contaminated[,1])


sample.sheet.full.aogc<-read.delim(the.sample.sheet.aogc,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
sample.sheet.full.aogc[1:5,]
posns<-match(sample.sheet.full[,"ParticipantCode"],sample.sheet.full.aogc[,"PATIENT"])
missing<-is.na(posns)
sum(missing)
sample.sheet.full.aogc<-sample.sheet.full.aogc[posns[!missing],]
dim(sample.sheet.full.aogc)

table(sample.sheet.full.aogc$CENTRE) ## center 15 is saliva (NZ cohort)
mean((sample.sheet.full.aogc$AGE_SCAN),na.rm=TRUE) 
range((sample.sheet.full.aogc$AGE_SCAN),na.rm=TRUE)

################################################################## GET COVERAGE STATS

## 203    AOGC-02-0395           38.32074
## 222    AOGC-02-0441           47.39313
## 249    AOGC-02-0512           36.96861
## 250    AOGC-02-0513           26.59086
## 270    AOGC-03-0017           49.88467
## 284    AOGC-03-0046           49.70396
## 286    AOGC-03-0049           49.17130
## 407    AOGC-08-0167           48.66813
## 408    AOGC-08-0169           44.18966

## bad.coverage<-the.coverage[,target]<50

posns<-match(sample.sheet.full[,"ParticipantCode"],coverage[,"Sample"])
missing<-is.na(posns)
sum(missing)
sample.sheet.full[missing,"ParticipantCode"]
## ## sample.sheet.full[posns[missing],"ParticipantCode"]
## ## coverage[,"Sample"]
## ## sample.sheet.full[missing,]
## ## coverage<-coverage[posns,]
## colnames(coverage)
the.coverage<-coverage[posns,]

the.coverage<-cbind(sample.sheet.full,the.coverage)
the.coverage[a.control,][1:5,]

target<-"percent.ccds.gt.10"
## target<-"total_reads"
 
     
mean(the.coverage[a.control & !bad.cov ,target])
median(the.coverage[a.control,target])
range(the.coverage[a.control,target])
## sd(the.coverage[a.control,target])
## table(the.coverage[a.control,"Capture.Method"])
bad.cov<-as.numeric(the.coverage[,target])<70

mean(the.coverage[!a.control,target])
median(the.coverage[!a.control,target])
range(the.coverage[!a.control,target])
## sd(the.coverage[!a.control,target])
## table(the.coverage[!a.control,"Capture.Method"])

the.coverage[as.numeric(the.coverage[a.control,target])<70,c("Sample",target)]

bad.coverge.sample<-the.coverage[as.numeric(the.coverage[a.control,target])<70,"Sample"]
bad.coverge.sample<-bad.coverge.sample[!is.na(bad.coverge.sample)]
## mean(the.coverage[!a.control,c("percent.ccds.gt.10")])

## the.coverage[a.control,c("Sample","percent.ccds.gt.10")]
## the.coverage[!a.control,c("Sample","percent.ccds.gt.10")]


########################################## coverage stats









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


bad.samples<-sample.sheet.full[sample.sheet.full[,"SampleProject"]=="Asian-AML","ParticipantCode"] # "1"  "9"  "63" "74" "83" "84" "99"

tail(sample.sheet.full)
test<-((sample.sheet.full[,"Race"]!="W" | is.na(sample.sheet.full[,"Race"])) & sample.sheet.full[,"AffectionStatus"]==2) 

sample.sheet.full[test,]

bad.samples<-sample.sheet.full[test,"ParticipantCode"] # "1"  "9"  "63" "74" "83" "84" "99"

contaminated.aogc<-sample.sheet.full[sample.sheet.full[,"ParticipantCode"] %in% contaminated[,1],"ParticipantCode"]


bad.samples<-c(bad.samples,contaminated.aogc)


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
for(ifam in 1:length(fam)){

the.extension<-paste(fam[ifam],project.extension,"$",sep="")
project.files<-files[grepl(the.extension ,files)]
print(sort(paste("Doing: ",project.files,sep=""))) # project.files<-project.files[1:22]

indels<-{}
the.col<-{}
project.files
#
ichr<-1


for(ichr in 1:length(project.files)){ ### loop over chromosomes

setwd(analysis.dir)

## grep("Gene.Name",a.indel[,16])
## save(list=c("column.labels"),file="column.labels.RData")
## load("column.labels.RData")
################## fast read ###########
column.labels<-read.delim(project.files[ichr],header=F,nrows=1,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
num.vars<-dim(column.labels)[2]
a.indel<-scan(project.files[ichr],what=character(num.vars),skip=1,sep="\t",fill=TRUE,na.strings="",quote="\"")
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

a.indel[1:50,"chr"]
key<-build.key(a.indel,core.ann)
rownames(a.indel)<-key
rownames(a.indel)[1:5]





########################### REMOVE BAD SAMPLES HERE
## analysis.samples[1:5,]
bad.samples
bad.samples.labels<-expand.labels.to.samples(bad.samples,c("GT","AD","DP","GQ"),paste.after=TRUE)
a.indel<-a.indel[,colnames(a.indel)[!(colnames(a.indel) %in% bad.samples.labels)]]

#colnames(a.indel)[(colnames(a.indel) %in% bad.samples.labels)]
bad.samples
#tapply(gene.list[,"CHR"],gene.list[,"CHR"],length)
all.possible.samples<-gsub(".GT$","",colnames(a.indel)[grep(".GT$",colnames(a.indel))],perl=TRUE)
length(all.possible.samples)
pheno.types





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




## all.possible.samples

## colnames(a.indel)[1:20]
a.indel[1:5,1:10]

## col<-grepl("REF",a.indel[,"REF"])
## sum(col)
## a.indel[col,1:15]
## a.indel<-a.indel[!col,]


## colnames(a.indel)[1:50]
## dim(a.indel)
## sort(unique(a.indel[,"chr"]))
## sort(unique(a.indel[,"FILTER"]))

## test<-"AOGC-08-0287.GT"
## wanted<-a.indel[,test]!="0/0" | is.na(a.indel[,test])
## sum(wanted)
## a.indel[wanted,c("chr","start", "end","Consequence.Embl",test)][1:10,]

## sort(table(a.indel[wanted, "Consequence.Embl"]))
## sort(table(a.indel[wanted, "TYPE"]))
## test<-all.possible.samples
## ALT.Alleles.Control

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
common.hit.genes<-names(all.genes)[all.genes>250] # common.hit.genes<-names(all.genes)[1:4]
all.genes["SCN2A"]

###############################################

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
clusters<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/Final_FANC_clusters.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
clusters

#clusters.wanted<-c("Clinical","FANC - ACID")
clusters.wanted<-colnames(clusters)
ic<-1
snpinfo[1:5,]

cbind(unique(clusters[,22]),unique(clusters[,22]))
snpinfo[1:5,]



gene.aliases<-read.delim("/media/UQCCG/Software/annovar/humandb/Gene_symbol_aliases.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
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
a.indel[1:20,1:5]


#filtered.genotype(a.indel[1:20,],gsub(".GT$","",the.samples),prefix="",suffix="",20,0.02,0.98,0.20,0.80,7,2)



library(doMC)
num.cores<-6
num.bits<-num.cores
registerDoMC(cores=num.cores)

the.samples<-colnames(a.indel)[grepl(".GT$", colnames(a.indel))]

colnames(a.indel)[grepl(".GQ$", colnames(a.indel))]
while((dim(a.indel)[1] %% num.bits)< 2){num.bits<-num.bits+1} ### go don't get matrix issues
fil.genotypes<-foreach(a.indel.bit=iter(a.indel,by='row',chunksize=as.integer(dim(a.indel)[1]/num.bits) ), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype.old(a.indel.bit,gsub(".GT$","",the.samples),prefix="",suffix="",20,0.02,0.98,0.20,0.80,7,2)
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

    ## AML Control 
    ##  89     197 
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
if(is.null(dim(summary.geno.extra))){
  summary.geno.extra<-summary.geno
}else{
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
if(is.null(dim(summary.geno.extra))){
  summary.geno.extra<-summary.geno
}else{
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
if(is.null(dim(summary.geno.extra))){
  summary.geno.extra<-summary.geno
}else{
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
if(is.null(dim(summary.geno.extra))){
  summary.geno.extra<-summary.geno
}else{
  summary.geno.extra<-cbind(summary.geno.extra,summary.geno)
}


#########################################
#

######################################


a.indel[1:5,1:10]
summary.geno.extra[1:5,]
#colnames(summary.geno.extra)
#rownames(summary.geno.extra)<-key
#getHWE(obs_hets, obs_hom1, obs_hom2)

hw.p.control<-getHWE(summary.geno.extra[,"GENO.Control"]) ## used 16 CPUs
hw.p.control.filt<-getHWE(summary.geno.extra[,"GENO.Control.filt"]) ## used 16 CPUs


###
 
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
#ok.missing.test[places,]
ok.missing.test<-high.missing <= missing.threshold
ok.missing.test[1:5,]
ok.missing<-combine.boolean(ok.missing.test,c("AML","Control"),"AND") # was AND
ok.missing.filt<-combine.boolean(ok.missing.test,c("AML.filt","Control.filt"),"AND") # was AND
sum(!ok.missing)
sum(!ok.missing.filt)
ok.missing[1:5]
ok.missing.filt<-combine.boolean(ok.missing.test,c("AML.filt","Control.filt"),"AND") 

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


#maf.lt.all<-a.indel[,colnames(a.indel)[grepl("MAF.lt:",colnames(a.indel))]]
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



missense.coding<-test.for.coding.type(a.indel,geneanno.DB,missense.variant)
sum(missense.coding)
   ####################### mhairi CHECK
missense.coding.vep<-test.wanted.mutation(a.indel[,"Consequence.Embl"],missense.variant,delimit.by=",")  # filter.table.pholy[,"Consequence"] %in%  vep.coding
  ## a.indel[wanted.muts.coding.vep,]
sum(missense.coding.vep)


is.missense<-missense.coding.vep | missense.coding

qual[1:50,"PolyPhen.low"] # true is above 0.1
a.indel[1:50,"PolyPhen.scores"]
a.indel[1:50,c("PolyPhen.desc","SIFT.desc")]
qual[1:5,"SIFT.high"] 
a.indel[1:5,"SIFT.scores"]
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
if(!("gene_biotype" %in% colnames(a.indel))){colnames(a.indel)[colnames(a.indel)=="Gene.Biotype"]<- "gene_biotype"}

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




#########################  FREQUENCY FILTERS
######### given a 0.01 threshold 6sd would allow 10 alt alleles at 6sd
######## given a 0.005 threshold 6sd would allow 7 alt alleles at 6sd
  n<-476 # number of controls
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
p<-0.01
#p<-0.01 ########### set MAF threshols HEREX1

sd.thresh<-4
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
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
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
 if(target.pheno %in% case.control){
   for(icc in 1:length(case.control.classes)){
     recode<-  pheno[,target.pheno] %in% names(case.control.classes)[icc]
     pheno[recode,target.pheno]<-as.numeric(case.control.classes[icc])
   }}
  
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
ignore.FLT3<-TRUE
  if( ("13" %in% the.chr) & !("chr13:28626716:28626716:C:T:CREST" %in% key) & !ignore.FLT3 ){ # add flt3-ITD in chr213 and not already there

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

           help<-cbind(full.qual,bad.coding,maf.filter,rare.in.group,no.genotypes,in.common.hit.gene ,hw.controls.ok,on.x.y,unannotated.hits,not.flat.genotype,are.repeats,are.in.repeats,ok.missing,ok.missing.filt,is.unwound.geno,(ok.missing.filt | is.unwound.geno) ,hw.p.control.filt,rare.in.group.filt,no.genotypes.filt,rare.in.controls.filt
                       )

  } # ftl3 additions


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



######################################################################################################
n<-max(as.integer(summary.geno.extra[,"TOTAL.Alleles.Control"]))

p<-0.01
#p<-0.01 ########### set MAF threshols HEREX1

sd.thresh<-4
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
######################################################################################################




maf.lt.all[1:5,]
#maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.001"])
#maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.01"])


#maf.aogc.total<-as.logical(a.indel[,"MAF.lt:0.01"])
maf.aogc.total<-as.logical(a.indel[,maf.col])
## maf.aogc.total<-as.logical(a.indel[,"MAF.lt:0.005"])

maf.aogc.total[1:5]


is.missense<-missense.coding.vep | missense.coding

qual[1:50,"PolyPhen.low"] # true is above 0.1
a.indel[1:50,"PolyPhen.scores"]
a.indel[1:50,c("PolyPhen.desc","SIFT.desc")]
qual[1:5,"SIFT.high"] 
a.indel[1:5,"SIFT.scores"]

#predict.benign<-a.indel[,"SIFT.desc"]=="tolerated" & a.indel[,"PolyPhen.desc"]=="benign"

predict.benign<-qual[,"PolyPhen.low"]


sum((predict.benign))

is.benign.missense<-is.missense & !predict.benign
sum(is.benign.missense)

dim(a.indel)
table(a.indel[,"PolyPhen.desc"])
table(a.indel[,"SIFT.desc"])

pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene  & !on.x.y & !unannotated.hits & not.flat.genotype  & ok.missing.filt & hw.controls.ok.filt & !no.genotypes.filt  & rare.in.controls.filt & rare.in.controls  & !are.repeats & !are.in.repeats & maf.aogc.total #  & !is.benign.missense #  #  & !is.benign.missense  # & maf.aogc.total & rare.in.group

#& !are.repeats & !are.in.repeats

bad.genotypes<-c("chr6:35425714:35425714:-:C:indel","chr11:108126934:108126934:A:T:snp")
bad.genotypes %in% names(pass)
pass[ names(pass) %in% bad.genotypes]<-FALSE

pass.ori<-pass


sum(pass)
sum(pass.ori)

## summary.geno.extra[pass & !pass.3 ,c("GENO.AML","GENO.Control","GENO.AML.filt","GENO.Control.filt")][1:5,]
## maf.aogc.total[pass & !pass.3][1:5]
#pass<- full.qual             & maf.filter   & !in.common.hit.gene  & !on.x.y & !unannotated.hits & not.flat.genotype & !are.repeats & !are.in.repeats & ok.missing.filt & hw.controls.ok.filt & !no.genotypes.filt & rare.in.controls.filt & rare.in.group



help<-cbind(full.qual,bad.coding,maf.filter,rare.in.group,no.genotypes,in.common.hit.gene ,hw.controls.ok,on.x.y,unannotated.hits,not.flat.genotype,are.repeats,are.in.repeats,ok.missing,ok.missing.filt,is.unwound.geno,(ok.missing.filt | is.unwound.geno) ,hw.p.control.filt,rare.in.group.filt,no.genotypes.filt,rare.in.controls.filt )

 dim(fil.genotypes)
sum(pass)
length(pass)
dim(snpinfo.ori)
snpinfo[1:5,]
a.indel[1:5,1:50]
pass[1:5]
dim(a.indel)




################################# GEFOS FILTERING cause sending all
## snpinfo[grep("chr13:28626716:28626716:C:T:CREST",snpinfo[,"Name"]),]
## snpinfo.ori[grep("chr13:28626716:28626716:C:T:CREST",snpinfo.ori[,"Name"]),]

## table(a.indel[pass,"refGene::gene"])
## table(a.indel[pass,"Gene.name"])
;

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
formula

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

meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]
meta.results.burden[1:50,]


the.order<-     order(meta.results.skat[,"p"])
meta.results.skat<-  meta.results.skat[the.order,]
meta.results.skat[1:50,]

the.order<-     order(meta.results.skatO[,"p"])
sum(is.na(meta.results.skatO[,"p"])) ## bad p-values shoudl not happen
meta.results.skatO<-  meta.results.skatO[the.order,]
meta.results.skatO[1:50,]

meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,]

##         meta.results.burden.gene[meta.results.burden.gene[,"gene"] %in% fanc.genes,]
## snpinfo.ori<-snpinfo 
## meta.results.burden.gene[meta.results.burden.gene[,"gene"] %in% other.clusters[,2],]
## meta.results.burden.gene[meta.results.burden.gene[,"gene"] %in% other.clusters[,3],]

        ## meta.results.skatO.gene[meta.results.skatO.gene[,"gene"] %in% clinical.genes,]
        ## meta.results.skatO.gene[meta.results.skatO.gene[,"gene"] %in% fanc.genes,]




## meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,]



setwd(analysis.dir)
getwd()

snap.file<-"coding.0.01.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN_wSTRAT"
snap.file<-"coding.0.001.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN_wSTRAT"

snap.file<-"coding.0.01.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN_wSTRAT_no_benign"
snap.file<-"coding.0.001.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN_wSTRAT_no_benign"

snap.file<-"coding.0.01.all.geno.all.filters_no.imput_paper_wSTRAT"
snap.file<-"coding.0.001.all.geno.all.filters_no.imput_paper_wSTRAT"

snap.file<-"coding.0.01.all.geno.all.filters_no.imput_paper_wSTRAT_no_benign"
snap.file<-"coding.0.001.all.geno.all.filters_no.imput_paper_wSTRAT_no_benign"


snap.file<-"coding.0.001.all.geno.all.filters_no.imput_paper_wSTRAT"
snap.file<-"coding.0.001.all.geno.all.filters_no.imput_paper_wSTRAT_rare"

snap.file<-"coding.0.01.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN"
snap.file<-"coding.0.001.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN"

snap.file<-"coding.0.01.all.geno.all.filters_no.imput_paper_no_benign"
snap.file<-"coding.0.001.all.geno.all.filters_no.imput_paper"

snap.file<-"coding.0.01.all.geno.all.filters"
snap.file<-"coding.0.001.all.geno.all.filters_no.imput"

snap.file<-"coding.0.01.all.geno.all.filters_no.imput_HC_indels"

write.table(meta.results.burden[1:50,],file=paste("Burden","Top50",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,],file=paste("Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


write.table(meta.results.skat[1:50,],file=paste("Skat",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skatO[1:50,],file=paste("SkatO",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skatO[1:50,],file=paste("SkatO",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skatO[meta.results.skat[,"gene"] %in% clusters.wanted,],file=paste("Skat","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,],file=paste("SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,],file=paste("Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


colnames(a.indel)[c(1:8,13,16,28,7,30,34,37:42,43,14,32,33)] # 1276,1310
annotations<-a.indel[,c(1:8,13,16,28,7,30,34,37:42,43,14,32,33)]
      
save(list=c("case.control","snpinfo.ori","formula","clusters","pheno.types","ipheno","clusters.wanted","genotypes","p","meta.results.skat","meta.results.skatO","meta.results.burden","pheno","target.pheno.col","snpinfo","fil.genotypes","pass","high.missing.table","a.indel","help","key","summary.geno.extra","full.qual","bad.effect","maf.filter","in.common.hit.gene","on.x.y","unannotated.hits","not.flat.genotype","are.repeats","are.in.repeats","ok.missing","hw.controls.ok.filt","no.genotypes","rare.in.Control","rare.in.Control.filt","in.any.normal","in.any.normal.filt","are.in.repeats.back","are.in.repeats.forward","all.genes"),file=paste(paste(project.files[ichr],".",pheno.types[ipheno],".",snap.file,".small_final.RData",sep="")) )
getwd()

save.image(file=paste(snap.file,"RData",sep="."))

## save(list=c("clusters.wanted"),file="clusters.wanted.RData")
## getwd()
## load(paste(snap.file,"RData",sep="."))


## setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-11-04_AML_TCGA_Replication/Analysis")
## load("AML_TCGA_image_paper.coding.0.01.all.geno.all.filters_no.imput_paper_TCGA_REP.RData")

## setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis")
load("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/coding.0.001.all.geno.all.filters_no.imput_paper_wSTRAT_no_benign.RData")
load("AML_HC_image_paper.coding.0.01.all.geno.all.filters_no.imput_paper.RData")
load("AML_AOGC_image.coding.0.001.all.geno.all.filters.NEW.RData")
load("coding.0.001.all.geno.all.filters_no.imput_paper_wSTRAT_no_benign.RData")
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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

load(paste(snap.file,"RData",sep="."))
options(width=200)
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
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
## test<-c("MYCBP2","SLC25A24","TMCO3","C13orf35")
## test<-c("MYCBP2","SLC25A24","TMCO3","C13orf35")
## test<-c("chr22:41252435-41252687:ST13")
## test<-c("SETD8")
test<-c("SEC61A1","ST14","GPANK1","EEF1A2")
test<-c("C19orf40")
test<-c("FANCP")
test<-c("IDH1")
test<-c("clinical") # not sig after coverage filtering at 0.01 0.1997688 from 0.000142
test<-clinical.genes
test<-fanc.genes

snpinfo[1:5,]
a.cluster<-"random.745"
test<-snpinfo[snpinfo[,"cluster"]==a.cluster,"gene"]
test

 test<-c("LOC100268168") #,"NOTCH1")
snpinfo[1:5,]
## meta.results.skat[meta.results.skat[,"gene"] %in% test,]


meta.results.burden[1:20,]
meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]



test<-snpinfo[snpinfo[,"cluster"]==a.cluster,"cluster"]
meta.results.burden[meta.results.burden[,"gene"] %in% test,]


to.unwind<-c(meta.results.burden[1:50,"gene"],meta.results.skatO[1:50,"gene"])
#to.unwind<-c("FANCD2_minimal_mono_ubi") #, "MCM7", "RNPC3")
to.unwind<-c("FANC_complex.all") # to.unwind<-meta.results.burden[8,"gene"]
#to.unwind<-c("BLM.Complex_AND_Checkpoint")
#to.unwind<-c("NPM1")
## to.unwind<-c("FLT3")
## to.unwind<-c("Clinical")
to.unwind<-c("DDX41","TET2", "GATA2", "ASXL1", "NOTCH1", "IDH1", "JAK2","MET")

to.unwind %in% meta.results.burden[,"gene"]
to.unwind %in% meta.results.skatO[,"gene"]


dim( meta.results.burden)
to.unwind<-c(clusters.wanted[!(clusters.wanted %in% c("Ubin.proteo","lipid_raft","caveolae","Checkpoint_extendedx1","Checkpoint_extendedx2"))])

to.unwind
#to.unwind.name<-to.unwind
to.unwind.name<-"ALL_clusters_ALL_mutations"
# to.unwind.name<-"Collaboration_Genes"
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

# the.genes.burden
write.table(the.genes.burden,file=paste(to.unwind.name,"conponents:","Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

the.genes.burden<-meta.results.skatO[meta.results.skatO[,"gene"] %in% the.genes,]
#the.genes.burden
write.table(the.genes.burden,file=paste(paste(to.unwind.name,collapse="."),"conponents:","SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
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
#meta.results.burden.ex
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
#the.gt<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GT",sep=".")
#the.gq<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GQ",sep=".")
the.gq<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"DP",sep=".")

quality.cases[check]<-paste(a.indel.sub[posn,the.gq],collapse=",")

## a.indel[posn,the.gq]
## a.indel[posn,the.gt]
## a.indel[posn,the.dp]
}

if(muts.in.controls[check]!=""){
#the.gt<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"GT",sep=".")
#the.gq<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"GQ",sep=".")
the.gq<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"DP",sep=".")
quality.controls[check]<-paste(a.indel.sub[posn,the.gq],collapse=",")

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
## summary.geno.extra[figure,]
## annotations[figure,]
## help[figure,]

dim(meta.results.burden.ex)
#out<-cbind(meta.results.burden.ex,a.indel[figure,c(1:6,16,43,28,7,30,34,37:42)],summary.geno.extra[figure,c("GENO.AML","GENO.Control","GENO.AML.filt","GENO.Control.filt")],help[figure,],muts.in.cases,muts.in.controls)

maf.old<-annotations[,"MAF.lt:0.01" ]
maf.new<-maf.lt.all[,"MAF.lt:0.01"]
maf.aogc<-a.indel[,"AOGC-NGS_ALL::maf"]
a.functions<-a.indel[,c("PolyPhen.scores","SIFT.scores","PolyPhen.desc","SIFT.desc")]
                             
out<-cbind(meta.results.burden.ex,maf.new[figure],maf.old[figure],maf.aogc[figure],a.functions[figure,],is.benign.missense[figure],annotations[figure,],summary.geno.extra[figure,c("GENO.AML","GENO.Control","GENO.AML.filt","GENO.Control.filt")],help[figure,],muts.in.cases,quality.cases,muts.in.controls,quality.controls)

#out<-cbind(meta.results.burden.ex,annotations[figure,],muts.in.cases,muts.in.controls)

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



meta.results.burden[1:50,]
meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]

the.order<-     order(meta.results.skat[,"p"])
meta.results.skat<-  meta.results.skat[the.order,]
meta.results.skat[1:50,]

the.order<-     order(meta.results.skatO[,"p"])
sum(is.na(meta.results.skatO[,"p"])) ## bad p-values shoudl not happen
meta.results.skatO<-  meta.results.skatO[the.order,]
meta.results.skatO[1:50,]


dim(clusters)
snpinfo.sub<-snpinfo[snpinfo[,"cluster"] %in% clusters.wanted,]
genes.cl<-unique(snpinfo.sub[,"gene"])
genes.cl<-genes.cl[!(genes.cl %in% clusters.wanted)]
genes.cl
clusters.wanted

genes.and.clusters<-c(genes.cl,clusters.wanted)
meta.results.burden[1:5,]

#################################### just want to plot

subset<-meta.results.burden[ !(meta.results.burden[,"gene"] %in% genes.and.clusters) ,]
subset[1:5,]

dups<-duplicated(subset[,"gene"])
sum(dups)

z<-qchisq(subset[ ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) ## if have no chisq valuse
#z0<-qchisq(meta.results.skatO[ !(meta.results.skatO[,"gene"] %in% genes.and.clusters ) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) 
z[1:5]
subset[1:10,]

p.val<-as.numeric(subset[ ,"p"])
par(mfrow=c(1,1),font=2,font.lab=2,font.axis=2,mgp=c(3.5,1,0),mar=c(5,5,4,2)+0.1)


##  z<-rchisq(length(p.val), df=1, ncp = 0) ## null test


median(z,na.rm=TRUE)/0.456  #1.071491
median(z0,na.rm=TRUE)/0.456  #1.071491



################## p-values
z0=qnorm(p.val/2)
lambda = round(median(z0^2)/0.454,3)
lambda



## source("http://bioconductor.org/biocLite.R") 
## biocLite("GWASTools")
## setRep
## qq.Plot(pvals)



## Reads data
## S <- read.table(input,header=F)
## if (stat_type == "Z")
##    z=S[,1]
## if (stat_type == "CHISQ")
##    z=sqrt(S[,1])
## if (stat_type == "PVAL")
##    z0=qnorm(meta.results.skatO[,"p"]/2)
## ## calculates lambda
lambda = round(median(z0^2)/.454,3)
## lambda
range(z)


the.plot<-my.qq.plot(z,dist="chisq",df=1,ncp=0,col="blue",ylab="Observed chi-squared value",xlab="Expected chi-squared value",main="",cex=1,xlim=c(0,22),ylim=c(0,80),cex.lab=2.0,cex.axis=2.0,font.lab=2,font.axis=2,lwd=2,line="robust",plot.it=TRUE) # function defined below


z.all<-qchisq(meta.results.burden[ !(meta.results.burden[,"gene"] %in% clusters.wanted) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)
range(z.all)
qq<-  qq.data(z.all,plot.it=FALSE)       ## qq plot used same method as in car library
points(qq$x,qq$y,col="magenta",pch=21)

symbols<-meta.results.burden[!(meta.results.burden[,"gene"] %in% clusters.wanted),"gene"]

#symbols<-meta.results.skatO[,"gene"]
#####annotate curve
selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="red",cex=1,offset=1,atpen='TRUE') ##plate row col symbol
selected.data<-identify(qq$x,qq$y,labels=labels[qq$ord],col="red",cex=1,atpen='TRUE') ## sybmol
selected.data<-identify(qq$x,qq$y,labels=as.character(round(data.in[qq$ord],2)),col="forestgreen",cex=1.25,atpen='TRUE') # observed score
#####


leg.txt<-c("All Genes","Remove Clinical Genes")

legend(2,60,leg.txt,col=c("magenta","blue"),lty=c(-1,-1,2),pch=c(1,1,-1),cex=2.25)

label<-"AML_paper_0.01_all_filters_NO_STRAT"


savePlot(paste(label,"tiff",sep="."),type="tiff")
savePlot(paste(label,"png",sep="."),type="png")
save.image("qq.AS.paper.final.RData")



o <- -(log10(p.val))
e <- -log10( 1:length(o)/length(o) )

plot(e,o, pch=23, cex=.4, bg="black",main=hlabel, ylab="Observed -log10 P value",xlab="Expected -log10 P value")
abline(coef=c(0,1), col=1, lwd=2)



## qq <- parse.pvals.qq(gene.region.compare[,"p.skatO"],lim=lim)
## qq <- parse.pvals.qq(all.res.gene[,"p.skatO"],lim=lim)

qq <- parse.pvals.qq(gene.region.compare[,"p.skatO"],lim=lim)
qq <- parse.pvals.qq(all.res.gene[,"p.skatO"],lim=lim)
 
ylab = expression(Observed~~-log[10](italic(p)))
xlab = expression(Expected~~-log[10](italic(p)))
 
plot(qq$e,
     qq$o,
     xlim=c(0,lim),ylim=c(0,lim),
     pch=20,col='deepskyblue',
     xlab=xlab,ylab=ylab, main="Gene Based")
     abline(coef=c(0,1), col=1, lwd=2)

savePlot("FilteredGeneBased.jpg",type="jpeg")

#write.table(geno.all,file=paste(project.name,fam[ifam],"geno.all.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

#geno.all<-read.delim(paste(project.name,fam[ifam],"geno.all.txt",sep="."),header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)


parse.pvals.qq <- function(pvector,lim=7) {

  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )

  o[o>lim] <- lim

  out <- list(o=o,e=e)

  return(out)
}






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




















############################################### META ANALYSIS
############################################### META ANALYSIS
############################################### META ANALYSIS
############################################### META ANALYSIS
############################################### META ANALYSIS
############################################### META ANALYSIS
  1   2  TCGA
476 102

  1   2  DISCOVERY
197  89




## setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-11-04_AML_TCGA_Replication/Analysis/meta_analyis/")
## files<-c("Burden.clusters.coding.0.01.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN.txt","Burden.clusters.coding.0.01.all.geno.all.filters_no.imput_paper.txt")
## sizes<-c(150,96)

setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014/2014-11-04_AML_TCGA_Replication/Analysis/meta_analysis/")
files<-c("Burden.clusters.coding.0.01.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN_wSTRAT_no_benign.txt","Burden.clusters.coding.0.01.all.geno.all.filters_no.imput_paper_wSTRAT_no_benign.txt")
sizes<-c(102,89)

traits<-c("AML")
i<-1


#### NEED TP HAND MODIFY FILES TO CHNAGES NMISS TO NUMBERS ABOVE


#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-11-04_AML_TCGA_Replication/Analysis/meta_analyis/ALL_clusters.conponents:.Burden.clusters.coding.0.001.all.geno.all.filters_no.imput_paper.txt


#setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014/2014-11-04_AML_TCGA_Replication/Analysis/meta_analyis")

system(paste("sed s/chip.TRAIT/",files[1],"/ run_config_for_P_and_BETA_TEMPLATE.txt > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.TRAIT/",files[2],"/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("sed s/chip.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt1",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
paste("sample.size.CONFIG",traits[i],"txt1",sep=".")
system(paste("./metal ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("cp ","META.SAMPLE.SIZE.TRAIT1.tbl","META.SAMPLE.SIZE.TRAIT.p_0.01.txt", sep=" "))



system(paste("sed s/chip.TRAIT/",files[1],"/ run_config_for_inverse_varience_TEMPLATE.txt > ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.TRAIT/",files[2],"/ ",paste("STDERR.CONFIG",traits[i],"txt",sep=".")," > ",paste("STDERR.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("sed s/chip.NMISS/","nmiss","/ ",paste("STDERR.CONFIG",traits[i],"txt1",sep=".")," > ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.NMISS/","nmiss","/ ",paste("STDERR.CONFIG",traits[i],"txt",sep=".")," > ",paste("STDERR.CONFIG",traits[i],"txt1",sep="."),sep=""))
paste("STDERR.CONFIG",traits[i],"txt1",sep=".")
system(paste("./metal ",paste("STDERR.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("cp ","META.STDERR.TRAIT1.tbl","META.STDERR.TRAIT.p_0.01.txt", sep=" "))



tcga<-read.delim(files[1],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
dis<-read.delim(files[2],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
meta<-read.delim("META.STDERR.TRAIT.p_0.01.txt",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)



tcga[1:5,]
dis[1:5,]
meta[1:5,]

colnames(tcga)<-paste(colnames(tcga),"TCGA",sep=".")
colnames(dis)<-paste(colnames(dis),"DIS",sep=".")

posns<-match(meta[,"MarkerName"], tcga[,1])
missing<-is.na(posns)
sum(missing)
meta[missing,"MarkerName"]
tcga<-tcga[posns,]

posns<-match(meta[,"MarkerName"], dis[,1])
missing<-is.na(posns)
sum(missing)
meta[missing,"MarkerName"]
dis<-dis[posns,]

dim(dis)
dim(tcga)

meta<-cbind(meta,tcga,dis)
meta[1:5,]
getwd()
write.table(meta,file="META.STDERR.TRAIT.p_0.01.SUMMARY.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)






)
bad.genotypes %in% names(pass)
pass[ names(pass) %in% bad.genotypes]<-FALSE

pass.ori<-pass


sum(pass)
sum(pass.ori)

## summary.geno.extra[pass & !pass.3 ,c("GENO.AML","GENO.Control","GENO.AML.filt","GENO.Control.filt")][1:5,]
## maf.aogc.total[pass & !pass.3][1:5]
#pass<- full.qual             & maf.filter   & !in.common.hit.gene  & !on.x.y & !unannotated.hits & not.flat.genotype & !are.repeats & !are.in.repeats & ok.missing.filt & hw.controls.ok.filt & !no.genotypes.filt & rare.in.controls.filt & rare.in.group



help<-cbind(full.qual,bad.coding,maf.filter,rare.in.group,no.genotypes,in.common.hit.gene ,hw.controls.ok,on.x.y,unannotated.hits,not.flat.genotype,are.repeats,are.in.repeats,ok.missing,ok.missing.filt,is.unwound.geno,(ok.missing.filt | is.unwound.geno) ,hw.p.control.filt,rare.in.group.filt,no.genotypes.filt,rare.in.controls.filt )

 dim(fil.genotypes)
sum(pass)
length(pass)
dim(snpinfo.ori)
snpinfo[1:5,]
a.indel[1:5,1:50]
pass[1:5]
dim(a.indel)




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
formula

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

meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]
meta.results.burden[1:50,]


the.order<-     order(meta.results.skat[,"p"])
meta.results.skat<-  meta.results.skat[the.order,]
meta.results.skat[1:50,]

the.order<-     order(meta.results.skatO[,"p"])
sum(is.na(meta.results.skatO[,"p"])) ## bad p-values shoudl not happen
meta.results.skatO<-  meta.results.skatO[the.order,]
meta.results.skatO[1:50,]

meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,]

##         meta.results.burden.gene[meta.results.burden.gene[,"gene"] %in% fanc.genes,]
## snpinfo.ori<-snpinfo 
## meta.results.burden.gene[meta.results.burden.gene[,"gene"] %in% other.clusters[,2],]
## meta.results.burden.gene[meta.results.burden.gene[,"gene"] %in% other.clusters[,3],]

        ## meta.results.skatO.gene[meta.results.skatO.gene[,"gene"] %in% clinical.genes,]
        ## meta.results.skatO.gene[meta.results.skatO.gene[,"gene"] %in% fanc.genes,]




## meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,]



setwd(analysis.dir)
getwd()

snap.file<-"coding.0.01.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN_wSTRAT"
snap.file<-"coding.0.001.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN_wSTRAT"

snap.file<-"coding.0.01.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN_wSTRAT_no_benign"
snap.file<-"coding.0.001.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN_wSTRAT_no_benign"

snap.file<-"coding.0.01.all.geno.all.filters_no.imput_paper_wSTRAT"
snap.file<-"coding.0.001.all.geno.all.filters_no.imput_paper_wSTRAT"

snap.file<-"coding.0.01.all.geno.all.filters_no.imput_paper_wSTRAT_no_benign"
snap.file<-"coding.0.001.all.geno.all.filters_no.imput_paper_wSTRAT_no_benign"


snap.file<-"coding.0.001.all.geno.all.filters_no.imput_paper_wSTRAT"
snap.file<-"coding.0.001.all.geno.all.filters_no.imput_paper_wSTRAT_rare"

snap.file<-"coding.0.01.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN"
snap.file<-"coding.0.001.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN"

snap.file<-"coding.0.01.all.geno.all.filters_no.imput_paper_no_benign"
snap.file<-"coding.0.001.all.geno.all.filters_no.imput_paper"

snap.file<-"coding.0.01.all.geno.all.filters"
snap.file<-"coding.0.001.all.geno.all.filters_no.imput"

snap.file<-"coding.0.01.all.geno.all.filters_no.imput_HC_indels"

write.table(meta.results.burden[1:50,],file=paste("Burden","Top50",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,],file=paste("Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


write.table(meta.results.skat[1:50,],file=paste("Skat",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skatO[1:50,],file=paste("SkatO",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skatO[1:50,],file=paste("SkatO",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,],file=paste("SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,],file=paste("Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


colnames(a.indel)[c(1:8,13,16,28,7,30,34,37:42,43,14,32,33)] # 1276,1310
annotations<-a.indel[,c(1:8,13,16,28,7,30,34,37:42,43,14,32,33)]
      
save(list=c("case.control","snpinfo.ori","formula","clusters","pheno.types","ipheno","clusters.wanted","genotypes","p","meta.results.skat","meta.results.skatO","meta.results.burden","pheno","target.pheno.col","snpinfo","fil.genotypes","pass","high.missing.table","a.indel","help","key","summary.geno.extra","full.qual","bad.effect","maf.filter","in.common.hit.gene","on.x.y","unannotated.hits","not.flat.genotype","are.repeats","are.in.repeats","ok.missing","hw.controls.ok.filt","no.genotypes","rare.in.Control","rare.in.Control.filt","in.any.normal","in.any.normal.filt","are.in.repeats.back","are.in.repeats.forward","all.genes"),file=paste(paste(project.files[ichr],".",pheno.types[ipheno],".",snap.file,".small_final.RData",sep="")) )
getwd()

save.image(file=paste(snap.file,"RData",sep="."))

## save(list=c("clusters.wanted"),file="clusters.wanted.RData")
## getwd()
## load(paste(snap.file,"RData",sep="."))


## setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-11-04_AML_TCGA_Replication/Analysis")
## load("AML_TCGA_image_paper.coding.0.01.all.geno.all.filters_no.imput_paper_TCGA_REP.RData")

## setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis")
load("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/coding.0.001.all.geno.all.filters_no.imput_paper_wSTRAT_no_benign.RData")
load("AML_HC_image_paper.coding.0.01.all.geno.all.filters_no.imput_paper.RData")
load("AML_AOGC_image.coding.0.001.all.geno.all.filters.NEW.RData")
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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

load(paste(snap.file,"RData",sep="."))
options(width=200)
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
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
## test<-c("MYCBP2","SLC25A24","TMCO3","C13orf35")
## test<-c("MYCBP2","SLC25A24","TMCO3","C13orf35")
## test<-c("chr22:41252435-41252687:ST13")
## test<-c("SETD8")
test<-c("SEC61A1","ST14","GPANK1","EEF1A2")
test<-c("C19orf40")
test<-c("FANCP")
test<-c("IDH1")
test<-c("clinical") # not sig after coverage filtering at 0.01 0.1997688 from 0.000142
test<-clinical.genes
test<-fanc.genes

snpinfo[1:5,]
a.cluster<-"random.745"
test<-snpinfo[snpinfo[,"cluster"]==a.cluster,"gene"]
test

 test<-c("LOC100268168") #,"NOTCH1")
snpinfo[1:5,]
## meta.results.skat[meta.results.skat[,"gene"] %in% test,]


meta.results.burden[1:20,]
meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]



test<-snpinfo[snpinfo[,"cluster"]==a.cluster,"cluster"]
meta.results.burden[meta.results.burden[,"gene"] %in% test,]


to.unwind<-c(meta.results.burden[1:50,"gene"],meta.results.skatO[1:50,"gene"])
#to.unwind<-c("FANCD2_minimal_mono_ubi") #, "MCM7", "RNPC3")
to.unwind<-c("FANC_complex.all") # to.unwind<-meta.results.burden[8,"gene"]
#to.unwind<-c("BLM.Complex_AND_Checkpoint")
#to.unwind<-c("NPM1")
## to.unwind<-c("FLT3")
## to.unwind<-c("DDX41","TET2", "GATA2", "ASXL1", "NOTCH1", "IDH1", "JAK2")
to.unwind<-c(clusters.wanted[!(clusters.wanted %in% c("Ubin.proteo","lipid_raft","caveolae","Checkpoint_extendedx1","Checkpoint_extendedx2"))])

to.unwind
#to.unwind.name<-to.unwind
to.unwind.name<-"ALL_clusters_ALL_mutations"
to.unwind.name<-"Collaboration_genes"
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

#the.genes.burden
write.table(the.genes.burden,file=paste(to.unwind.name,"conponents:","Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

the.genes.burden<-meta.results.skatO[meta.results.skatO[,"gene"] %in% the.genes,]
#the.genes.burden
write.table(the.genes.burden,file=paste(paste(to.unwind.name,collapse="."),"conponents:","SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
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
#meta.results.burden.ex
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
#the.gt<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GT",sep=".")
#the.gq<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GQ",sep=".")
the.gq<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"DP",sep=".")

quality.cases[check]<-paste(a.indel.sub[posn,the.gq],collapse=",")

## a.indel[posn,the.gq]
## a.indel[posn,the.gt]
## a.indel[posn,the.dp]
}

if(muts.in.controls[check]!=""){
#the.gt<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"GT",sep=".")
#the.gq<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"GQ",sep=".")
the.gq<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"DP",sep=".")
quality.controls[check]<-paste(a.indel.sub[posn,the.gq],collapse=",")

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
## summary.geno.extra[figure,]
## annotations[figure,]
## help[figure,]

dim(meta.results.burden.ex)
#out<-cbind(meta.results.burden.ex,a.indel[figure,c(1:6,16,43,28,7,30,34,37:42)],summary.geno.extra[figure,c("GENO.AML","GENO.Control","GENO.AML.filt","GENO.Control.filt")],help[figure,],muts.in.cases,muts.in.controls)

maf.old<-annotations[,"MAF.lt:0.01" ]
maf.new<-maf.lt.all[,"MAF.lt:0.01"]
maf.aogc<-a.indel[,"AOGC-NGS_ALL::maf"]
a.functions<-a.indel[,c("PolyPhen.scores","SIFT.scores","PolyPhen.desc","SIFT.desc")]
                             
out<-cbind(meta.results.burden.ex,maf.new[figure],maf.old[figure],maf.aogc[figure],a.functions[figure,],is.benign.missense[figure],annotations[figure,],summary.geno.extra[figure,c("GENO.AML","GENO.Control","GENO.AML.filt","GENO.Control.filt")],help[figure,],muts.in.cases,quality.cases,muts.in.controls,quality.controls)

#out<-cbind(meta.results.burden.ex,annotations[figure,],muts.in.cases,muts.in.controls)

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



meta.results.burden[1:50,]
meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]

the.order<-     order(meta.results.skat[,"p"])
meta.results.skat<-  meta.results.skat[the.order,]
meta.results.skat[1:50,]

the.order<-     order(meta.results.skatO[,"p"])
sum(is.na(meta.results.skatO[,"p"])) ## bad p-values shoudl not happen
meta.results.skatO<-  meta.results.skatO[the.order,]
meta.results.skatO[1:50,]


dim(clusters)
snpinfo.sub<-snpinfo[snpinfo[,"cluster"] %in% clusters.wanted,]
genes.cl<-unique(snpinfo.sub[,"gene"])
genes.cl<-genes.cl[!(genes.cl %in% clusters.wanted)]
genes.cl
clusters.wanted

genes.and.clusters<-c(genes.cl,clusters.wanted)
meta.results.burden[1:5,]

#################################### just want to plot

subset<-meta.results.burden[ !(meta.results.burden[,"gene"] %in% genes.and.clusters) ,]
subset[1:5,]

dups<-duplicated(subset[,"gene"])
sum(dups)

z<-qchisq(subset[ ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) ## if have no chisq valuse
#z0<-qchisq(meta.results.skatO[ !(meta.results.skatO[,"gene"] %in% genes.and.clusters ) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) 
z[1:5]
subset[1:10,]

p.val<-as.numeric(subset[ ,"p"])
par(mfrow=c(1,1),font=2,font.lab=2,font.axis=2,mgp=c(3.5,1,0),mar=c(5,5,4,2)+0.1)


##  z<-rchisq(length(p.val), df=1, ncp = 0) ## null test


median(z,na.rm=TRUE)/0.456  #1.071491
median(z0,na.rm=TRUE)/0.456  #1.071491



################## p-values
z0=qnorm(p.val/2)
lambda = round(median(z0^2)/0.454,3)
lambda



## source("http://bioconductor.org/biocLite.R") 
## biocLite("GWASTools")
## setRep
## qq.Plot(pvals)



## Reads data
## S <- read.table(input,header=F)
## if (stat_type == "Z")
##    z=S[,1]
## if (stat_type == "CHISQ")
##    z=sqrt(S[,1])
## if (stat_type == "PVAL")
##    z0=qnorm(meta.results.skatO[,"p"]/2)
## ## calculates lambda
lambda = round(median(z0^2)/.454,3)
## lambda
range(z)


the.plot<-my.qq.plot(z,dist="chisq",df=1,ncp=0,col="blue",ylab="Observed chi-squared value",xlab="Expected chi-squared value",main="",cex=1,xlim=c(0,22),ylim=c(0,80),cex.lab=2.0,cex.axis=2.0,font.lab=2,font.axis=2,lwd=2,line="robust",plot.it=TRUE) # function defined below


z.all<-qchisq(meta.results.burden[ !(meta.results.burden[,"gene"] %in% clusters.wanted) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)
range(z.all)
qq<-  qq.data(z.all,plot.it=FALSE)       ## qq plot used same method as in car library
points(qq$x,qq$y,col="magenta",pch=21)

symbols<-meta.results.burden[!(meta.results.burden[,"gene"] %in% clusters.wanted),"gene"]

#symbols<-meta.results.skatO[,"gene"]
#####annotate curve
selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="red",cex=1,offset=1,atpen='TRUE') ##plate row col symbol
selected.data<-identify(qq$x,qq$y,labels=labels[qq$ord],col="red",cex=1,atpen='TRUE') ## sybmol
selected.data<-identify(qq$x,qq$y,labels=as.character(round(data.in[qq$ord],2)),col="forestgreen",cex=1.25,atpen='TRUE') # observed score
#####


leg.txt<-c("All Genes","Remove Clinical Genes")

legend(2,60,leg.txt,col=c("magenta","blue"),lty=c(-1,-1,2),pch=c(1,1,-1),cex=2.25)

label<-"AML_paper_0.01_all_filters_NO_STRAT"


savePlot(paste(label,"tiff",sep="."),type="tiff")
savePlot(paste(label,"png",sep="."),type="png")
save.image("qq.AS.paper.final.RData")



o <- -(log10(p.val))
e <- -log10( 1:length(o)/length(o) )

plot(e,o, pch=23, cex=.4, bg="black",main=hlabel, ylab="Observed -log10 P value",xlab="Expected -log10 P value")
abline(coef=c(0,1), col=1, lwd=2)



## qq <- parse.pvals.qq(gene.region.compare[,"p.skatO"],lim=lim)
## qq <- parse.pvals.qq(all.res.gene[,"p.skatO"],lim=lim)

qq <- parse.pvals.qq(gene.region.compare[,"p.skatO"],lim=lim)
qq <- parse.pvals.qq(all.res.gene[,"p.skatO"],lim=lim)
 
ylab = expression(Observed~~-log[10](italic(p)))
xlab = expression(Expected~~-log[10](italic(p)))
 
plot(qq$e,
     qq$o,
     xlim=c(0,lim),ylim=c(0,lim),
     pch=20,col='deepskyblue',
     xlab=xlab,ylab=ylab, main="Gene Based")
     abline(coef=c(0,1), col=1, lwd=2)

savePlot("FilteredGeneBased.jpg",type="jpeg")

#write.table(geno.all,file=paste(project.name,fam[ifam],"geno.all.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

#geno.all<-read.delim(paste(project.name,fam[ifam],"geno.all.txt",sep="."),header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)


parse.pvals.qq <- function(pvector,lim=7) {

  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )

  o[o>lim] <- lim

  out <- list(o=o,e=e)

  return(out)
}






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




















############################################### META ANALYSIS
############################################### META ANALYSIS
############################################### META ANALYSIS
############################################### META ANALYSIS
############################################### META ANALYSIS
############################################### META ANALYSIS
  1   2  TCGA
476 102

  1   2  DISCOVERY
197  89




## setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-11-04_AML_TCGA_Replication/Analysis/meta_analyis/")
## files<-c("Burden.clusters.coding.0.01.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN.txt","Burden.clusters.coding.0.01.all.geno.all.filters_no.imput_paper.txt")
## sizes<-c(150,96)

setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014/2014-11-04_AML_TCGA_Replication/Analysis/meta_analysis/")
files<-c("Burden.clusters.coding.0.01.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN_wSTRAT_no_benign.txt","Burden.clusters.coding.0.01.all.geno.all.filters_no.imput_paper_wSTRAT_no_benign.txt")
sizes<-c(102,89)

traits<-c("AML")
i<-1


#### NEED TP HAND MODIFY FILES TO CHNAGES NMISS TO NUMBERS ABOVE


#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-11-04_AML_TCGA_Replication/Analysis/meta_analyis/ALL_clusters.conponents:.Burden.clusters.coding.0.001.all.geno.all.filters_no.imput_paper.txt


#setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014/2014-11-04_AML_TCGA_Replication/Analysis/meta_analyis")

system(paste("sed s/chip.TRAIT/",files[1],"/ run_config_for_P_and_BETA_TEMPLATE.txt > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.TRAIT/",files[2],"/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("sed s/chip.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt1",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
paste("sample.size.CONFIG",traits[i],"txt1",sep=".")
system(paste("./metal ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("cp ","META.SAMPLE.SIZE.TRAIT1.tbl","META.SAMPLE.SIZE.TRAIT.p_0.01.txt", sep=" "))



system(paste("sed s/chip.TRAIT/",files[1],"/ run_config_for_inverse_varience_TEMPLATE.txt > ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.TRAIT/",files[2],"/ ",paste("STDERR.CONFIG",traits[i],"txt",sep=".")," > ",paste("STDERR.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("sed s/chip.NMISS/","nmiss","/ ",paste("STDERR.CONFIG",traits[i],"txt1",sep=".")," > ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.NMISS/","nmiss","/ ",paste("STDERR.CONFIG",traits[i],"txt",sep=".")," > ",paste("STDERR.CONFIG",traits[i],"txt1",sep="."),sep=""))
paste("STDERR.CONFIG",traits[i],"txt1",sep=".")
system(paste("./metal ",paste("STDERR.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("cp ","META.STDERR.TRAIT1.tbl","META.STDERR.TRAIT.p_0.01.txt", sep=" "))



tcga<-read.delim(files[1],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
dis<-read.delim(files[2],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
meta<-read.delim("META.STDERR.TRAIT.p_0.01.txt",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)



tcga[1:5,]
dis[1:5,]
meta[1:5,]

colnames(tcga)<-paste(colnames(tcga),"TCGA",sep=".")
colnames(dis)<-paste(colnames(dis),"DIS",sep=".")

posns<-match(meta[,"MarkerName"], tcga[,1])
missing<-is.na(posns)
sum(missing)
meta[missing,"MarkerName"]
tcga<-tcga[posns,]

posns<-match(meta[,"MarkerName"], dis[,1])
missing<-is.na(posns)
sum(missing)
meta[missing,"MarkerName"]
dis<-dis[posns,]

dim(dis)
dim(tcga)

meta<-cbind(meta,tcga,dis)
meta[1:5,]
getwd()
write.table(meta,file="META.STDERR.TRAIT.p_0.01.SUMMARY.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)






























## files<-c("Burden.clusters.coding.0.001.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN.txt","Burden.clusters.coding.0.001.all.geno.all.filters_no.imput_paper.txt")
## sizes<-c(150,96)
## traits<-c("AML")


setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014/2014-11-04_AML_TCGA_Replication/Analysis/meta_analysis/")
files<-c("Burden.clusters.coding.0.001.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN_wSTRAT_no_benign.txt","Burden.clusters.coding.0.001.all.geno.all.filters_no.imput_paper_wSTRAT_no_benign.txt")
sizes<-c(102,89)


i<-1

system(paste("sed s/chip.TRAIT/",files[1],"/ run_config_for_P_and_BETA_TEMPLATE.txt > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.TRAIT/",files[2],"/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("sed s/chip.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt1",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
paste("sample.size.CONFIG",traits[i],"txt1",sep=".")
system(paste("./metal ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("cp ","META.SAMPLE.SIZE.TRAIT1.tbl","META.SAMPLE.SIZE.TRAIT.p_0.001.txt", sep=" "))



system(paste("sed s/chip.TRAIT/",files[1],"/ run_config_for_inverse_varience_TEMPLATE.txt > ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.TRAIT/",files[2],"/ ",paste("STDERR.CONFIG",traits[i],"txt",sep=".")," > ",paste("STDERR.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("sed s/chip.NMISS/","nmiss","/ ",paste("STDERR.CONFIG",traits[i],"txt1",sep=".")," > ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.NMISS/","nmiss","/ ",paste("STDERR.CONFIG",traits[i],"txt",sep=".")," > ",paste("STDERR.CONFIG",traits[i],"txt1",sep="."),sep=""))
paste("STDERR.CONFIG",traits[i],"txt1",sep=".")
system(paste("./metal ",paste("STDERR.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("cp ","META.STDERR.TRAIT1.tbl","META.STDERR.TRAIT.p_0.001.txt", sep=" "))



tcga<-read.delim(files[1],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
dis<-read.delim(files[2],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
meta<-read.delim("META.STDERR.TRAIT.p_0.001.txt",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)

tcga[1:5,]
dis[1:5,]
meta[1:5,]

colnames(tcga)<-paste(colnames(tcga),"TCGA",sep=".")
colnames(dis)<-paste(colnames(dis),"DIS",sep=".")

posns<-match(meta[,"MarkerName"], tcga[,1])
missing<-is.na(posns)
sum(missing)
meta[missing,"MarkerName"]
tcga<-tcga[posns,]

posns<-match(meta[,"MarkerName"], dis[,1])
missing<-is.na(posns)
sum(missing)
meta[missing,"MarkerName"]
dis<-dis[posns,]

dim(dis)
dim(tcga)

meta<-cbind(meta,tcga,dis)
meta[1:5,]
getwd()
write.table(meta,file="META.STDERR.TRAIT.p_0.001.SUMMARY.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



##########################################################################

setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014/2014-11-04_AML_TCGA_Replication/Analysis/meta_analysis/")
files<-c("SkatO.clusters.coding.0.001.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN_wSTRAT_no_benign.txt","SkatO.clusters.coding.0.001.all.geno.all.filters_no.imput_paper_wSTRAT_no_benign.txt")
sizes<-c(102,89)


i<-1

system(paste("sed s/chip.TRAIT/",files[1],"/ run_config_for_P_and_BETA_TEMPLATE.txt > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.TRAIT/",files[2],"/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("sed s/chip.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt1",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
paste("sample.size.CONFIG",traits[i],"txt1",sep=".")
system(paste("./metal ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("cp ","META.SAMPLE.SIZE.TRAIT1.tbl","META.SAMPLE.SIZE.TRAIT.SKATO.p_0.001.txt", sep=" "))


tcga<-read.delim(files[1],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
dis<-read.delim(files[2],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
meta<-read.delim("META.SAMPLE.SIZE.TRAIT.SKATO.p_0.001.txt",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)

tcga[1:5,]
dis[1:5,]
meta[1:5,]

colnames(tcga)<-paste(colnames(tcga),"TCGA",sep=".")
colnames(dis)<-paste(colnames(dis),"DIS",sep=".")

posns<-match(meta[,"MarkerName"], tcga[,1])
missing<-is.na(posns)
sum(missing)
meta[missing,"MarkerName"]
tcga<-tcga[posns,]

posns<-match(meta[,"MarkerName"], dis[,1])
missing<-is.na(posns)
sum(missing)
meta[missing,"MarkerName"]
dis<-dis[posns,]

dim(dis)
dim(tcga)

meta<-cbind(meta,tcga,dis)
meta[1:5,]
getwd()
write.table(meta,file="META.SAMPLE.SIZE.TRAIT.SKATO.p_0.001.SUMMARY.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



##########################################################################

setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014/2014-11-04_AML_TCGA_Replication/Analysis/meta_analysis/")
files<-c("SkatO.clusters.coding.0.01.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN_wSTRAT_no_benign.txt","SkatO.clusters.coding.0.01.all.geno.all.filters_no.imput_paper_wSTRAT_no_benign.txt")
sizes<-c(102,89)


i<-1

system(paste("sed s/chip.TRAIT/",files[1],"/ run_config_for_P_and_BETA_TEMPLATE.txt > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.TRAIT/",files[2],"/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("sed s/chip.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt1",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
paste("sample.size.CONFIG",traits[i],"txt1",sep=".")
system(paste("./metal ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("cp ","META.SAMPLE.SIZE.TRAIT1.tbl","META.SAMPLE.SIZE.TRAIT.SKATO.p_0.01.txt", sep=" "))


tcga<-read.delim(files[1],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
dis<-read.delim(files[2],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
meta<-read.delim("META.SAMPLE.SIZE.TRAIT.SKATO.p_0.01.txt",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)

tcga[1:5,]
dis[1:5,]
meta[1:5,]

colnames(tcga)<-paste(colnames(tcga),"TCGA",sep=".")
colnames(dis)<-paste(colnames(dis),"DIS",sep=".")

posns<-match(meta[,"MarkerName"], tcga[,1])
missing<-is.na(posns)
sum(missing)
meta[missing,"MarkerName"]
tcga<-tcga[posns,]

posns<-match(meta[,"MarkerName"], dis[,1])
missing<-is.na(posns)
sum(missing)
meta[missing,"MarkerName"]
dis<-dis[posns,]

dim(dis)
dim(tcga)

meta<-cbind(meta,tcga,dis)
meta[1:5,]
getwd()
write.table(meta,file="META.SAMPLE.SIZE.TRAIT.SKATO.p_0.01.SUMMARY.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)




################## NO BENIGN
################## NO BENIGN
################## NO BENIGN
################## NO BENIGN
################## NO BENIGN
################## NO BENIGN
################## NO BENIGN
################## NO BENIGN


  1   2  TCGA
476 102

  1   2  DISCOVERY
197  89




Burden 0.001 stderr & size stderror reported


setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014/2014-11-04_AML_TCGA_Replication/Analysis/meta_analysis/")
files<-c("Burden.clusters.coding.0.001.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN_wSTRAT.txt","Burden.clusters.coding.0.001.all.geno.all.filters_no.imput_paper_wSTRAT.txt")
sizes<-c(102,89)
traits<-c("AML")

i<-1

system(paste("sed s/chip.TRAIT/",files[1],"/ run_config_for_P_and_BETA_TEMPLATE.txt > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.TRAIT/",files[2],"/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("sed s/chip.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt1",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
paste("sample.size.CONFIG",traits[i],"txt1",sep=".")
system(paste("./metal ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("cp ","META.SAMPLE.SIZE.TRAIT1.tbl","META.SAMPLE.SIZE.TRAIT.WITH_BENIGN_MISSENSE_p_0.001.txt", sep=" "))



system(paste("sed s/chip.TRAIT/",files[1],"/ run_config_for_inverse_varience_TEMPLATE.txt > ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.TRAIT/",files[2],"/ ",paste("STDERR.CONFIG",traits[i],"txt",sep=".")," > ",paste("STDERR.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("sed s/chip.NMISS/","nmiss","/ ",paste("STDERR.CONFIG",traits[i],"txt1",sep=".")," > ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.NMISS/","nmiss","/ ",paste("STDERR.CONFIG",traits[i],"txt",sep=".")," > ",paste("STDERR.CONFIG",traits[i],"txt1",sep="."),sep=""))
paste("STDERR.CONFIG",traits[i],"txt1",sep=".")
system(paste("./metal ",paste("STDERR.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("cp ","META.STDERR.TRAIT1.tbl","META.STDERR.TRAIT.WITH_BENIGN_MISSENSE_p_0.001.txt", sep=" "))



tcga<-read.delim(files[1],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
dis<-read.delim(files[2],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
meta<-read.delim("META.STDERR.TRAIT.WITH_BENIGN_MISSENSE_p_0.001.txt",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)

tcga[1:5,]
dis[1:5,]
meta[1:5,]

colnames(tcga)<-paste(colnames(tcga),"TCGA",sep=".")
colnames(dis)<-paste(colnames(dis),"DIS",sep=".")

posns<-match(meta[,"MarkerName"], tcga[,1])
missing<-is.na(posns)
sum(missing)
meta[missing,"MarkerName"]
tcga<-tcga[posns,]

posns<-match(meta[,"MarkerName"], dis[,1])
missing<-is.na(posns)
sum(missing)
meta[missing,"MarkerName"]
dis<-dis[posns,]

dim(dis)
dim(tcga)

meta<-cbind(meta,tcga,dis)
meta[1:5,]
getwd()
write.table(meta,file="META.STDERR.TRAIT.WITH_BENIGN_MISSENSE_p_0.001.SUMMARY.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

########################
Burden 0.01 stderr & size , std error reported


setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014/2014-11-04_AML_TCGA_Replication/Analysis/meta_analysis/")
files<-c("Burden.clusters.coding.0.01.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN_wSTRAT.txt","Burden.clusters.coding.0.01.all.geno.all.filters_no.imput_paper_wSTRAT.txt")
sizes<-c(102,89)

traits<-c("AML")
i<-1


system(paste("sed s/chip.TRAIT/",files[1],"/ run_config_for_P_and_BETA_TEMPLATE.txt > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.TRAIT/",files[2],"/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("sed s/chip.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt1",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
paste("sample.size.CONFIG",traits[i],"txt1",sep=".")
system(paste("./metal ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("cp ","META.SAMPLE.SIZE.TRAIT1.tbl","META.SAMPLE.SIZE.TRAIT.WITH_BENIGN_MISSENSE_p_0.01.txt", sep=" "))



system(paste("sed s/chip.TRAIT/",files[1],"/ run_config_for_inverse_varience_TEMPLATE.txt > ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.TRAIT/",files[2],"/ ",paste("STDERR.CONFIG",traits[i],"txt",sep=".")," > ",paste("STDERR.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("sed s/chip.NMISS/","nmiss","/ ",paste("STDERR.CONFIG",traits[i],"txt1",sep=".")," > ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.NMISS/","nmiss","/ ",paste("STDERR.CONFIG",traits[i],"txt",sep=".")," > ",paste("STDERR.CONFIG",traits[i],"txt1",sep="."),sep=""))
paste("STDERR.CONFIG",traits[i],"txt1",sep=".")
system(paste("./metal ",paste("STDERR.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("cp ","META.STDERR.TRAIT1.tbl","META.STDERR.TRAIT.WITH_BENIGN_MISSENSE_p_0.01.txt", sep=" "))



tcga<-read.delim(files[1],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
dis<-read.delim(files[2],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
meta<-read.delim("META.STDERR.TRAIT.WITH_BENIGN_MISSENSE_p_0.01.txt",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)



tcga[1:5,]
dis[1:5,]
meta[1:5,]

colnames(tcga)<-paste(colnames(tcga),"TCGA",sep=".")
colnames(dis)<-paste(colnames(dis),"DIS",sep=".")

posns<-match(meta[,"MarkerName"], tcga[,1])
missing<-is.na(posns)
sum(missing)
meta[missing,"MarkerName"]
tcga<-tcga[posns,]

posns<-match(meta[,"MarkerName"], dis[,1])
missing<-is.na(posns)
sum(missing)
meta[missing,"MarkerName"]
dis<-dis[posns,]

dim(dis)
dim(tcga)

meta<-cbind(meta,tcga,dis)
meta[1:5,]
getwd()
write.table(meta,file="META.STDERR.TRAIT.WITH_BENIGN_MISSENSE_p_0.01.SUMMARY.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)






##########################################################################


SKAT 0.001  size , size reported

setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014/2014-11-04_AML_TCGA_Replication/Analysis/meta_analysis/")
files<-c("SkatO.clusters.coding.0.001.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN_wSTRAT.txt","SkatO.clusters.coding.0.001.all.geno.all.filters_no.imput_paper_wSTRAT.txt")
sizes<-c(102,89)


i<-1

system(paste("sed s/chip.TRAIT/",files[1],"/ run_config_for_P_and_BETA_TEMPLATE.txt > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.TRAIT/",files[2],"/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("sed s/chip.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt1",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
paste("sample.size.CONFIG",traits[i],"txt1",sep=".")
system(paste("./metal ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("cp ","META.SAMPLE.SIZE.TRAIT1.tbl","META.SAMPLE.SIZE.TRAIT.SKATO.WITH_BENIGN_MISSENSE_p_0.001.txt", sep=" "))


tcga<-read.delim(files[1],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
dis<-read.delim(files[2],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
meta<-read.delim("META.SAMPLE.SIZE.TRAIT.SKATO.WITH_BENIGN_MISSENSE_p_0.001.txt",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)

tcga[1:5,]
dis[1:5,]
meta[1:5,]

colnames(tcga)<-paste(colnames(tcga),"TCGA",sep=".")
colnames(dis)<-paste(colnames(dis),"DIS",sep=".")

posns<-match(meta[,"MarkerName"], tcga[,1])
missing<-is.na(posns)
sum(missing)
meta[missing,"MarkerName"]
tcga<-tcga[posns,]

posns<-match(meta[,"MarkerName"], dis[,1])
missing<-is.na(posns)
sum(missing)
meta[missing,"MarkerName"]
dis<-dis[posns,]

dim(dis)
dim(tcga)

meta<-cbind(meta,tcga,dis)
meta[1:5,]
getwd()
write.table(meta,file="META.SAMPLE.SIZE.TRAIT.SKATO.WITH_BENIGN_MISSENSE_p_0.001.SUMMARY.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



##########################################################################

SKAT 0.01  size , size reported

setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014/2014-11-04_AML_TCGA_Replication/Analysis/meta_analysis/")
files<-c("SkatO.clusters.coding.0.01.all.geno.all.filters_no.imput_paper_TCGA_REP_CLEAN_wSTRAT.txt","SkatO.clusters.coding.0.01.all.geno.all.filters_no.imput_paper_wSTRAT.txt")
sizes<-c(102,89)


i<-1

system(paste("sed s/chip.TRAIT/",files[1],"/ run_config_for_P_and_BETA_TEMPLATE.txt > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.TRAIT/",files[2],"/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("sed s/chip.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt1",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("sed s/seq.NMISS/","nmiss","/ ",paste("sample.size.CONFIG",traits[i],"txt",sep=".")," > ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
paste("sample.size.CONFIG",traits[i],"txt1",sep=".")
system(paste("./metal ",paste("sample.size.CONFIG",traits[i],"txt1",sep="."),sep=""))
system(paste("cp ","META.SAMPLE.SIZE.TRAIT1.tbl","META.SAMPLE.SIZE.TRAIT.SKATO.WITH_BENIGN_MISSENSE_p_0.01.txt", sep=" "))


tcga<-read.delim(files[1],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
dis<-read.delim(files[2],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
meta<-read.delim("META.SAMPLE.SIZE.TRAIT.SKATO.WITH_BENIGN_MISSENSE_p_0.01.txt",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)

tcga[1:5,]
dis[1:5,]
meta[1:5,]

colnames(tcga)<-paste(colnames(tcga),"TCGA",sep=".")
colnames(dis)<-paste(colnames(dis),"DIS",sep=".")

posns<-match(meta[,"MarkerName"], tcga[,1])
missing<-is.na(posns)
sum(missing)
meta[missing,"MarkerName"]
tcga<-tcga[posns,]

posns<-match(meta[,"MarkerName"], dis[,1])
missing<-is.na(posns)
sum(missing)
meta[missing,"MarkerName"]
dis<-dis[posns,]

dim(dis)
dim(tcga)

meta<-cbind(meta,tcga,dis)
meta[1:5,]
getwd()
write.table(meta,file="META.SAMPLE.SIZE.TRAIT.SKATO.WITH_BENIGN_MISSENSE_p_0.01.SUMMARY.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)






######################################################
######################################################
######################################################
######################################################
######################################################
######################################################









write.table(bad.samples,file="excluded_samples_TCGA.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



write.table(bad.samples,file="excluded_samples_DISCOVERY.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)




system(paste("sed s/TRAIT/common.",traits[i],"/ run_config_for_P_and_BETA_TEMPLATE.txt > ",paste("sample.size.CONFIG.common",traits[i],"txt",sep="."),sep=""))
system(paste("./metal ",paste("sample.size.CONFIG.common",traits[i],"txt",sep="."),sep=""))


system(paste("sed s/TRAIT/",traits[i],"/ run_config_for_inverse_varience_TEMPLATE.txt > ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("./metal ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))

system(paste("sed s/TRAIT/common.",traits[i],"/ run_config_for_inverse_varience_TEMPLATE.txt > ",paste("STDERR.CONFIG.common",traits[i],"txt",sep="."),sep=""))
system(paste("./metal ",paste("STDERR.CONFIG.common",traits[i],"txt",sep="."),sep=""))











####################################
chk<-out[,"refGene::location"]=="synonymous SNV"
chk2<-out[,"Consequence.Embl"]=="missense_variant"
sum(chk & chk2)
sort(table(out[chk,"Consequence.Embl"]))
out[chk & chk2,1:25][1:5,]

out
out.idh1<-out #"chr2:209113112:209113112:C:T:snp"

write.table(out,file="fanc-acid.summary.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(out,file="fanc-acid.4 missing.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
 write.table(out,file="fanc-acid.0.001-2.protein.summary.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
 write.table(out,file="clinical.0.001.protein.summary.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
out[1:5,]

  meta.results.burden.ex[  meta.results.burden.ex[,"gene"] %in% remove.repeats,]
## #meta.results.skatO.ex<-skatOMeta(cohort.seq.ex,burden.wts =1,SNPInfo = snpinfo)

grep(TRUE,diff)[1:10]
target<-c("chr1:2116899:2116899:A:G:snp","chr1:2121211:2121211:C:T:snp","chr15:89803883:89803883:T:C:snp","chr14:45652935:45652935:A:G:snp")
target<-"chr19:33465099:33465099:C:T:snp"
target<-"chr11:130059563:130059563:G:A:snp"
target<-grep("chr2:209113113:209113113:G",key)
target
#target<-key[diff][83:84] # target<-key[test]
target<-key[target]
target<-key[figure]

pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene  & !on.x.y & !unannotated.hits & not.flat.genotype & !are.repeats &
( ok.missing.filt | is.unwound.geno)   & hw.controls.ok.filt & !no.genotypes.filt   & !are.in.repeats & rare.in.controls.filt & rare.in.group

a.indel[target,1:40]

out<-cbind(a.indel[target,1:40],summary.geno.extra[target,c("GENO.AML","GENO.Control","GENO.AML.filt","GENO.Control.filt")])
help[target,]
pass[target]
ok.missing.filt[target]
maf.filter[target] 
high.missing[target,]
chk<-gsub(".GT$",".AD",the.samples[1:96])
chk.GT<-the.samples[1:96]

chk<-gsub(".GT$",".AD",the.samples)
chk.GT<-the.samples


a.indel[target,chk]
a.indel[target,chk.GT]
fil.genotypes[target,chk.GT]
summary.geno.extra[target,paste("GENO",c("AML","Control","AML.filt","Control.filt"),sep=".")]
hw.p.control.filt[target]

                                                AML Control  AML.filt Control.filt
chr2:209113113:209113113:G:A:snp:209113113:flat   0       0 0.5000000         0.29
chr2:209113113:209113113:G:A:snp:209113113        0       0 0.5104167         0.29
chr2:209113113:209113113:G:T:snp:209113113        0       0 0.5312500         0.29
## ## tapply(a.indel[pass,"Consequence.Embl"],a.indel[pass,"Consequence.Embl"],length)
## ## ##  dbeta(x, shape1, shape2, ncp = 0, log = FALSE)
## ## ## shape1, shape2: positive parameters of the Beta distribution.
## ## ##  ncp: non-centrality parameter.

extra.out<-c("FILTER","TYPE","MAF.lt:0.5","Consequence.Embl","wanted.muts","wanted.muts.coding",unique(global.quality.labs)[unique(global.quality.labs) %in% colnames(a.indel)],"Hetero.ALT.reads","Hetero.REF.reads","Hetero.Read.Balance","culprit")


        gene                         p        beta        se  cmafTotal   cmafUsed nsnpsTotal nsnpsUsed nmiss
1042     WT1 0.00007100568909052773730  0.94444308 0.2377219 0.13764770 0.13764770          9         9     0
1807     MLL 0.48467366086356300503013 -0.32851222 0.4701054 0.03378378 0.03378378         16        16     0
2187    KRAS 0.07033093737069703865750  1.59454406 0.8810739 0.01013514 0.01013514          5         5     0
3033    FLT3 0.00000000000005425888317  2.38660803 0.3173157 0.06925676 0.06925676         10        10     0
4420    IDH2 0.00000001305450061629625  2.32224967 0.4084620 0.04729730 0.04729730          5         5     0
5478    TP53 0.00004001465115838068794  1.99534188 0.4857925 0.02367167 0.02367167         11        11     0
11171 DNMT3A 0.00000000000019795761198  2.24740471 0.3057624 0.07794105 0.07794105         30        30     0
13567    KIT 0.01521559788883776638546  1.32404414 0.5455013 0.02370054 0.02370054         10        10     0
13787   TET2 0.00000047126510932962245  1.12239340 0.2228007 0.12718477 0.12718477         52        52     0
14788   NPM1 0.00000000000000003337779  1.69694129 0.2012002 0.10979730 0.10979730          6         6     0
17452   JAK2 0.84085058870232365357822  0.04984451 0.2482227 0.15168108 0.15168108         12        12     0
>



5697            TTC19 0.00000000000000000001958964 -1.8425116 0.19887755 0.29245955 0.29245955          6         6     0
15152            NPM1 0.00000000000000000960301848  0.7070687 0.08242233 0.26520270 0.26520270         15        15     0
13394            PFN2 0.00000000000000026646051843  1.9890686 0.24293636 0.12668919 0.12668919          6         5     0
3103             FLT3 0.00000000000005425888316660  2.3866080 0.31731574 0.06925676 0.06925676         10        10     0
11441          DNMT3A 0.00000000000019795761198317  2.2474047 0.30576235 0.07794105 0.07794105         30        30     0
15162    LOC100268168 0.00000000002476538050455725  1.2574656 0.18839142 0.32075104 0.32075104          7         7     0
2875             KSR2 0.00000000003062675194293715 -0.9080621 0.13668377 0.41430453 0.41430453         14        14     0
4910             GGA2 0.00000000010816520700968617  0.4888360 0.07572894 1.01858108 1.01858108         15        15     0
11254         PACSIN2 0.00000000021489113624525686  2.2058416 0.34736005 0.06624158 0.06624158          9         9     0
9451             CD1A 0.00000000127822168012818693  1.9147310 0.31543594 0.09819093 0.09819093          6         6     0
5429  ENSG00000269323 0.00000000390495998358280713 -1.3197635 0.22413815 0.13682432 0.13682432          1         1     0
4536             IDH2 0.00000001305450061629624858  2.3222497 0.40846204 0.04729730 0.04729730          5         5     0
16893         ATXN7L1 0.00000002337016933282044505  0.5516432 0.09877230 0.46168151 0.46168151         20        20     0
721         LOC728407 0.00000013493389842786308576  1.5188057 0.28808826 0.49166667 0.49166667          1         1     0
4800        LOC440335 0.00000017189409453387654263  1.6801551 0.32141208 0.08992448 0.08992448          2         2     0
10638          EEF1A2 0.00000030640515080703841849  1.2747271 0.24899620 0.18978953 0.18978953          2         2     0
209              PARG 0.00000070845983409388752649  1.5324555 0.30902322 0.47445532 0.47445532          3         3     0
3520           STXBP6 0.00000077629917146060417535  1.7819896 0.36063660 0.06460484 0.06460484          5         5     0
14512           NIPBL 0.00000105409135833204367882  1.4049795 0.28783122 0.09060888 0.09060888         16        16     0
17242           KMT2C 0.00000196638747493748888254 -0.9783022 0.20566185 0.26116071 0.26116071          2         2     0
2222             LDHB 0.00000202648349448258055131 -0.6908779 0.14542455 0.37440685 0.37440685          4         4     0
13836            KLF3 0.00000366262309136746226783  1.3769195 0.29741251 0.10642464 0.10642464          6         6     0
3025  ENSR00000430498 0.00000889885831862108461743  1.1566910 0.26037916 0.15593220 0.15593220          1         1     0
2712             FGD6 0.00000989160569885710713636  0.3892227 0.08806881 0.38446519 0.38446519         16        16     0
4123            INO80 0.00001022417375355785278319  1.1280452 0.25565475 0.11323372 0.11323372         17        17     0
12297            IDH1 0.00001125814375019362615672  1.5378205 0.35018271 0.03885135 0.03885135          7         7     0
16833         TSC22D4 0.00001337573547972480456471  1.8395756 0.42251542 0.04641653 0.04641653          5         5     0
540             WBP1L 0.00001442643020075763023451 -1.1895626 0.27426537 0.09799027 0.09799027          8         8     0
16583            OGDH 0.00001894088326779406772924  1.0058254 0.23516949 0.16048446 0.16048446         10        10     0
2062             CD27 0.00002418586493523308856993 -1.6657034 0.39450463 0.05574324 0.05574324          2         2     0
6907             CNN2 0.00003705952582898055575793 -0.9727442 0.23581183 0.09631838 0.09631838          7         7     0
5616             TP53 0.00004001465115838068794391  1.9953419 0.48579254 0.02367167 0.02367167         11        11     0
3191              ESD 0.00004143513725315481078517  0.8959166 0.21855207 0.20439189 0.20439189          4         4     0
12767          CTNNB1 0.00004377859525646177034103  0.7841088 0.19187427 0.25210448 0.25210448         11        11     0
8156             ZNF8 0.00004839096663123216608563 -1.1702579 0.28800923 0.10472973 0.10472973          6         6     0
18290           PTPN3 0.00004851169760004435047018 -0.7725056 0.19014662 0.47316903 0.47316903          8         8     0
6710            GAREM 0.00005049050955042321724500  1.0778286 0.26591091 0.54560811 0.54560811         10        10     0
10023          HNRNPU 0.00005088402190605824181686  0.6647065 0.16406313 0.29325553 0.29325553         10        10     0
13664          CCDC58 0.00005148194684455684465340  1.4406924 0.35583237 0.07094595 0.07094595          3         3     0
18168            SHC3 0.00005324718275289896197771 -0.8552526 0.21164920 0.17059658 0.17059658          4         4     0
1063              WT1 0.00007100568909052773730432  0.9444431 0.23772190 0.13764770 0.13764770          9         9     0
14396       WDFY3-AS1 0.00007448685883887118689903  1.3467120 0.33995181 0.07993197 0.07993197          1         1     0
4010            NEDD8 0.00007814162348186583608497  0.8274176 0.20947125 0.21428571 0.21428571          3         3     0
16229           TULP4 0.00008986033811572452945776  0.5485538 0.14006383 0.43248993 0.43248993         23        23     0
1972  ENSR00000558810 0.00009000456176326984095593  1.8173676 0.46407959 0.03938356 0.03938356          1         1     0
15987           PRDM1 0.00009393036056275626418353  0.3916096 0.10026469 0.82994055 0.82994055         18        18     0
15365           ATXN1 0.00009444772816356851680240  0.1766385 0.04524054 2.69360977 2.69360977         41        41     0
6679          C18orf8 0.00011118710945817823006372  3.1580161 0.81712596 0.01182432 0.01182432          6         6     0
12614            EMC3 0.00011295421737071749813779  0.9412806 0.24379606 0.28210032 0.28210032          4         4     0
14906         SLC23A1 0.00011398259503656368992960  3.1532978 0.81718759 0.01186484 0.01186484          8         7     0



          gene          p       beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
1016     FANCF 0.75904796  0.2957806 0.9642988 0.009504010 0.009504010          2         2     0
3495     FANCM 0.39824998  0.3771907 0.4465111 0.042229730 0.042229730         16        16     0
4404     FANCI 0.66334497 -0.2261604 0.5195511 0.030405405 0.030405405         12        12     0
5265     FANCA 0.40425225 -0.3841600 0.4605953 0.035472973 0.035472973         17        17     0
6362  C17orf70 0.05315640  1.2173415 0.6295569 0.020461234 0.020461234         10        10     0
6385    STRA13 0.30527330 -0.9827900 0.9586413 0.005669328 0.005669328          2         2     0
7280  C19orf40 0.08503474 -1.2392754 0.7195919 0.013183442 0.013183442          2         2     0
8073   C1orf86 0.57068855  0.2946712 0.5196677 0.030775238 0.030775238          7         7     0
11359    FANCL 0.17286745  0.6546883 0.4803121 0.032094595 0.032094595          9         9     0
12317   FANCD2 0.39288470 -0.5871702 0.6872329 0.016891892 0.016891892          7         7     0
15277    FANCE 0.48768622 -1.4850169 2.1398136 0.001689189 0.001689189          1         1     0
17575    FANCG 0.32554696 -1.4900680 1.5156478 0.003378378 0.003378378          2         2     0
17737    FANCC 0.65107059  0.1894990 0.4189902 0.041252909 0.041252909          6         6     0

################################### protein 0.01
        gene                          p       beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
13492   NPM1 0.000000000000000004666874  1.0474066 0.1209304 0.179054054 0.179054054          8         8     0
2743    FLT3 0.000000000000158669842946  2.4365276 0.3301671 0.065878378 0.065878378          8         8     0
10150 DNMT3A 0.000000000000626740759891  2.2837670 0.3174308 0.069471967 0.069471967         25        25     0
3983    IDH2 0.000000013054500616296249  2.3222497 0.4084620 0.047297297 0.047297297          5         5     0
12494   TET2 0.000016213366527191842669  1.2089637 0.2804029 0.070957437 0.070957437         40        40     0
4953    TP53 0.000108295608479553532970  1.9264148 0.4976258 0.021982481 0.021982481         10        10     0
12295    KIT 0.070330937370696955390770  1.5945441 0.8810739 0.010135135 0.010135135          5         5     0
936      WT1 0.070330937370697080290860  1.5945441 0.8810739 0.010135135 0.010135135          5         5     0
1989    KRAS 0.202970214166854234782988  1.5782177 1.2396312 0.005067568 0.005067568          2         2     0
15772   JAK2 0.506363538396016887865869 -0.4805368 0.7231439 0.015202703 0.015202703          8         8     0
1633     MLL 0.732588605446840124280072 -0.2056990 0.6020149 0.018581081 0.018581081         11        11     0

################################### protein 0.01
10317    FANCL 0.0004122055  1.6604799 0.4701054 0.033783784 0.033783784          8         8     0
5744  C17orf70 0.0531564049  1.2173415 0.6295569 0.020461234 0.020461234         10        10     0
5766    STRA13 0.3052732963 -0.9827900 0.9586413 0.005669328 0.005669328          2         2     0
15879    FANCG 0.3255469610 -1.4900680 1.5156478 0.003378378 0.003378378          2         2     0
11180   FANCD2 0.3928846996 -0.5871702 0.6872329 0.016891892 0.016891892          7         7     0
3157     FANCM 0.3982499824  0.3771907 0.4465111 0.042229730 0.042229730         16        16     0
4761     FANCA 0.4042522470 -0.3841600 0.4605953 0.035472973 0.035472973         17        17     0
916      FANCF 0.4498418852  0.8126484 1.0753889 0.006756757 0.006756757          1         1     0
6601  C19orf40 0.4876862217 -1.4850169 2.1398136 0.001689189 0.001689189          1         1     0
13812    FANCE 0.4876862217 -1.4850169 2.1398136 0.001689189 0.001689189          1         1     0
16019    FANCC 0.6244209596 -0.2771562 0.5660935 0.025337838 0.025337838          4         4     0
7354   C1orf86 0.6734397267  0.2314917 0.5492980 0.027254067 0.027254067          5         5     0
3970     FANCI 0.7840547028 -0.1462422 0.5336548 0.028716216 0.028716216         11        11     0


################################### protein 0.001
     gene                          p      beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
11649   NPM1 0.000000000000000004666874 1.0474066 0.1209304 0.179054054 0.179054054          8         8     0
2369    FLT3 0.000000000000000063908810 3.2306785 0.3865474 0.054054054 0.054054054          6         6     0
8721  DNMT3A 0.000000000000552230211637 2.3671831 0.3282382 0.062650013 0.062650013         21        21     0
3420    IDH2 0.000000000023130646061980 3.3067633 0.4946721 0.033783784 0.033783784          2         2     0
10766   TET2 0.000009553518252604718656 1.3067709 0.2951794 0.062511491 0.062511491         36        36     0
4265    TP53 0.000748277501256778697532 1.7680941 0.5244624 0.018604103 0.018604103          9         9     0
800      WT1 0.003655478822618108813297 3.1255708 1.0753889 0.006756757 0.006756757          4         4     0
10603    KIT 0.011978874707657898024404 3.1149033 1.2396312 0.005067568 0.005067568          3         3     0
1717    KRAS 0.202970214166854234782988 1.5782177 1.2396312 0.005067568 0.005067568          2         2     0
13607   JAK2 0.715431886461324495485314 0.3512669 0.9635083 0.008445946 0.008445946          5         5     0
1401     MLL 0.781879288747886591615099 0.1890884 0.6829481 0.013513514 0.013513514          8         8     0

     gene           p        beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
8874     FANCL 0.002711277  1.54750263 0.5160566 0.027027027 0.027027027          5         5     0
3408     FANCI 0.281888405  0.82393519 0.7656772 0.013513514 0.013513514          8         8     0
4950    STRA13 0.305273296 -0.98279000 0.9586413 0.005669328 0.005669328          2         2     0
13700    FANCG 0.325546961 -1.49006803 1.5156478 0.003378378 0.003378378          2         2     0
6349   C1orf86 0.401499884  0.49081068 0.5850328 0.023781175 0.023781175          3         3     0
2725     FANCM 0.434345506  0.56533746 0.7231439 0.015202703 0.015202703          8         8     0
11926    FANCE 0.487686222 -1.48501695 2.1398136 0.001689189 0.001689189          1         1     0
14150    FANCC 0.487686222 -1.48501695 2.1398136 0.001689189 0.001689189          1         1     0
4094     FANCA 0.648916787 -0.34858796 0.7656772 0.013513514 0.013513514          7         7     0
9624    FANCD2 0.973273037  0.04153204 1.2396312 0.005067568 0.005067568          3         3     0
4930  C17orf70 0.976555202  0.03643125 1.2396673 0.005084863 0.005084863          3         3     0

          gene                                       p      beta         se cmafTotal  cmafUsed nsnpsTotal nsnpsUsed nmiss
546    Clinical 0.0000000000000000000000000000004888122 0.8565784 0.07393617 0.4495089 0.4495089        104       104     0
908 FANC - ACID 0.0867467126344465128129357367470220197 0.4301199 0.25111881 0.1156164 0.1156164         43        43     0
> 
                gene                          p       beta        se  cmafTotal   cmafUsed nsnpsTotal nsnpsUsed nmiss
11649            NPM1 0.000000000000000004666874  1.0474066 0.1209304 0.17905405 0.17905405          8         8     0
2369             FLT3 0.000000000000000063908810  3.2306785 0.3865474 0.05405405 0.05405405          6         6     0
8721           DNMT3A 0.000000000000552230211637  2.3671831 0.3282382 0.06265001 0.06265001         21        21     0
3420             IDH2 0.000000000023130646061980  3.3067633 0.4946721 0.03378378 0.03378378          2         2     0
11550    LOC100268168 0.000000000023533911748733  1.2857857 0.1924188 0.31399428 0.31399428          3         3     0
4127  ENSG00000269323 0.000000003904959983582807 -1.3197635 0.2241381 0.13682432 0.13682432          1         1     0
12234         SERINC1 0.000000094211826775019504  1.8572919 0.3479672 0.12785741 0.12785741          2         2     0
544         LOC728407 0.000000134933898427863086  1.5188057 0.2880883 0.49166667 0.49166667          1         1     0
3636        LOC440335 0.000000171894094533876543  1.6801551 0.3214121 0.08992448 0.08992448          2         2     0
9362             IDH1 0.000000331598514265690453  3.2136150 0.6295600 0.02027027 0.02027027          3         3     0
4276             PER1 0.000000692890021435440133  1.7087936 0.3442826 0.08116412 0.08116412          6         6     0
8107           EEF1A2 0.000001137775237446866202  1.2327350 0.2533272 0.18641115 0.18641115          1         1     0
13094           KMT2C 0.000001966387474937488883 -0.9783022 0.2056618 0.26116071 0.26116071          2         2     0
10766            TET2 0.000009553518252604718656  1.3067709 0.2951794 0.06251149 0.06251149         36        36     0
8299             PCNT 0.000011442460687260778087  1.1341969 0.2584800 0.13554180 0.13554180         12        12     0
1952             LRP1 0.000034034260333051306034  0.8104806 0.1955489 0.22651456 0.22651456          9         9     0
7999            UBE2C 0.000041498170035095508933  1.9445022 0.4743873 0.03914591 0.03914591          1         1     0
10813        KIAA1109 0.000047255810233395838932  2.2193480 0.5454551 0.02364865 0.02364865         13        13     0
8207            RUNX1 0.000047326938703021625054  1.5047527 0.3698590 0.03909016 0.03909016         15        15     0
7302             DDR2 0.000058857795100619874381 -1.6531321 0.4114988 0.05084746 0.05084746          2         2     0
11153            BDP1 0.000061989348153295098726  0.9370378 0.2339607 0.20975925 0.20975925          9         9     0
7023             NRAS 0.000074914710439908353693  2.5831005 0.6522797 0.01525424 0.01525424          8         8     0
3412             KIF7 0.000100898997892848101368  1.2854867 0.3305936 0.09599800 0.09599800          8         8     0
7082         ADAMTSL4 0.000101601001243716289693  1.2695473 0.3266358 0.08385783 0.08385783          7         7     0
7934           CEP250 0.000111187109458178230064  3.1580161 0.8171260 0.01182432 0.01182432          7         7     0
9519           ANKMY1 0.000111187109458178230064  3.1580161 0.8171260 0.01182432 0.01182432          7         7     0
10803           USP53 0.000111187109458178230064  3.1580161 0.8171260 0.01182432 0.01182432          7         7     0
1707             LDHB 0.000116389200518012228438 -0.9555849 0.2479714 0.21114865 0.21114865          2         2     0
6797            CDCP2 0.000122802383916957911060  0.9933946 0.2586653 0.35780486 0.35780486          5         5     0
7359            ASTN1 0.000133839078173962654824  1.9595385 0.5130642 0.02372881 0.02372881         13        13     0
9623             EMC3 0.000182299524791116322091  0.9236984 0.2468224 0.28040541 0.28040541          3         3     0
2504            FARP1 0.000191424749198074471484  1.3362401 0.3582344 0.05514385 0.05514385          5         5     0
12470          TNRC18 0.000234083870525570765389  0.5140810 0.1397311 0.36964260 0.36964260         13        13     0
8218             TTC3 0.000247927027443518377750  1.1483377 0.3133774 0.09534780 0.09534780         10        10     0
14094          NOTCH1 0.000249298429951354749094  1.2113617 0.3307039 0.06834617 0.06834617         11        11     0
7977          GDAP1L1 0.000249478395604726900427 -1.2256653 0.3346257 0.07263514 0.07263514          1         1     0
11848          GPANK1 0.000314055267712844826938  1.2264781 0.3403649 0.07506803 0.07506803          5         5     0
12779            MCM7 0.000354372117981734016494  3.1471264 0.8810739 0.01013514 0.01013514          5         5     0
2037         C12orf50 0.000354372117981734287544  3.1471264 0.8810739 0.01013514 0.01013514          5         5     0
5693           SLC7A9 0.000354372117981735100696  3.1471264 0.8810739 0.01013514 0.01013514          6         6     0
10541           ARAP2 0.000354372117981735100696  3.1471264 0.8810739 0.01013514 0.01013514          6         6     0
11347         SLC23A1 0.000356907628275577887972  3.1455490 0.8810928 0.01014663 0.01014663          6         6     0
217            SPOCK2 0.000358679910400708907257 -0.9409022 0.2636499 0.12325957 0.12325957          3         3     0
9460           SPATA3 0.000366044773861016385728  3.1398917 0.8811439 0.01018702 0.01018702          5         5     0
3974              FUK 0.000368341828536694387754  1.7763957 0.4987375 0.02676248 0.02676248          5         5     0
13467          FER1L6 0.000461698838858031076170  2.1220873 0.6059563 0.02198244 0.02198244         12        12     0
5289             TLE6 0.000473051272974335946641 -1.4972329 0.4283230 0.04250552 0.04250552          3         3     0
3516            METRN 0.000488996172629530792726  0.6917143 0.1983858 0.24723979 0.24723979          5         5     0
5347            LONP1 0.000598160234910257077902  2.2889379 0.6668530 0.15473912 0.15473912          8         8     0
445            TCF7L2 0.000617325306994607092732  1.3700787 0.4001531 0.04918521 0.04918521          4         4     0
>


> meta.results.skatO.gene[1:50,]
                 gene                         p               pmin rho       cmaf nmiss nsnps errflag
13492            NPM1 0.00000000000000006231578  0.000000000000000   0 0.17905405     0     8       3
10150          DNMT3A 0.00000000063524399003710 -0.000000011481598   0 0.06947197     0    25       3
2743             FLT3 0.00000000071406168256950  0.000000000000000   0 0.06587838     0     8       3
3983             IDH2 0.00000000185436090349299 -0.000000001299192   0 0.04729730     0     5       3
4790  ENSG00000269323 0.00000000390495929759913  0.000000003904959   0 0.13682432     0     1       0
10898            IDH1 0.00000003199432113946245  0.000000336874034   1 0.02027027     0     3       0
14171         SERINC1 0.00000005535758817594879  0.000000090207034   1 0.12785741     0     2       0
4222        LOC440335 0.00000007760763769242383  0.000000169958773   1 0.08992448     0     2       0
641         LOC728407 0.00000013493379736382203  0.000000134933797   0 0.49166667     0     1       0
9443           EEF1A2 0.00000113777580490786562  0.000001137775805   0 0.18641115     0     1       0
1745             ST14 0.00000611062893059399721  0.000000774375733   0 0.09877338     0     6       0
8506             DDR2 0.00000635028110935179648  0.000041332534274   1 0.05593220     0     3       0
15171           KMT2C 0.00001624214541258868269  0.000001964479786   1 0.26116071     0     2       0
12494            TET2 0.00001784037973307827301  0.000016213028713   1 0.07095744     0    40       0
7888            CDCP2 0.00002905436713756232202  0.000004173148585   0 0.35950554     0     6       0
8073            RNPC3 0.00003485771459188366722  0.000064092765796   0 0.26271186     0     7       0
9318            UBE2C 0.00004149812741518308111  0.000041498127415   0 0.03914591     0     1       0
9554            RUNX1 0.00006414818616932162719  0.000047326710961   1 0.03909016     0    15       0
4967             PER1 0.00007598954113991147950  0.000076821033159   1 0.12576632     0    17       0
7113           ZNF880 0.00009357194667626277389  0.000083158238226   0 0.23370192     0     7       0
1977             LDHB 0.00011677481253043246702  0.000116389126280   1 0.21114865     0     2       0
8160             NRAS 0.00014215324891966822228  0.000074914501136   1 0.01525424     0     8       0
4953             TP53 0.00014976007670256414337  0.000108295491227   1 0.02198248     0    10       0
16350          NOTCH1 0.00015788371991295785304  0.000064771648038   0 0.09234495     0    17       0
9284          GDAP1L1 0.00021973013925981576796  0.000184702254448   1 0.07432432     0     2       0
6141             TLE6 0.00024465429727435966190  0.000189855323963   1 0.04757308     0     5       0
12108          FGFRL1 0.00030240431511864730848  0.000266572893705   0 0.26207047     0     5       0
11178            EMC3 0.00037495772538978234785  0.000182299540965   1 0.28040541     0     3       0
12936            BDP1 0.00039352602231099734166  0.000180463974624   1 0.24016465     0    18       0
14809            MCM7 0.00048124778478639620993  0.000111186970999   1 0.01182432     0     6       0
7996          COL24A1 0.00055122382437280048259  0.000286734556616   0 0.15423729     0    14       0
2245            NABP2 0.00067773351061196007369  0.000677733510612   0 0.11718750     0     1       0
8243            ANXA9 0.00069925444700109397334  0.000362147710139   1 0.01016949     0     4       0
265            SPOCK2 0.00070550090794654462655  0.000358679931356   1 0.12325957     0     3       0
10317           FANCL 0.00078542350021463982861  0.000412205725316   1 0.03378378     0     8       0
11930          LRRIQ4 0.00078647783891075081532  0.000731331535248   1 0.03716216     0     4       0
2271             LRP1 0.00088440249666421849421  0.000468516330787   1 0.27077381     0    23       0
2040         ADAMTS20 0.00099548907401215600520  0.000629621878342   1 0.15371622     0    13       0
9661             PCNT 0.00103788866818594455697  0.000523955889346   1 0.19473237     0    30       0
12996        ATP6AP1L 0.00109459889658832312809  0.000238452879667   1 0.01520270     0     4       0
11140            CHL1 0.00111354688852394233603  0.000776550554834   1 0.06418919     0     6       0
2874           COMMD6 0.00116585873464020574157  0.001165858734640   0 0.05323194     0     1       0
7037             VRK3 0.00119730526807206960442  0.001073935077417   1 0.01858108     0     7       0
2627         ATP6V0A2 0.00120585642429437236468  0.001079285181628   1 0.01689189     0     5       0
8234         ADAMTSL4 0.00124712487726952921813  0.000650456467526   1 0.10090030     0    15       0
4083            METRN 0.00126149029609892693836  0.000621033308648   1 0.24892898     0     6       0
7073            KLK14 0.00132600922668462458537  0.001141255779499   0 0.02612879     0     3       0
14285         PLEKHG1 0.00144443572977550477310  0.001022375924098   1 0.03209459     0    12       0
1764  ENSG00000254418 0.00144593322772500133995  0.001445933227725   0 0.20992366     0     1       0
13167         SLC23A1 0.00147425012502366554849  0.000356907836456   1 0.01014663     0     6       0




############################with repeats controled
        gene                          p      beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
11384   NPM1 0.000000000000000004666874 1.0474066 0.1209304 0.179054054 0.179054054          8         8     0
2326    FLT3 0.000000000000000063908810 3.2306785 0.3865474 0.054054054 0.054054054          6         6     0
8506  DNMT3A 0.000000000000552230211637 2.3671831 0.3282382 0.062650013 0.062650013         21        21     0
3331    IDH2 0.000000000023130646061980 3.3067633 0.4946721 0.033783784 0.033783784          2         2     0
10515   TET2 0.000009553518252604718656 1.3067709 0.2951794 0.062511491 0.062511491         36        36     0
4151    TP53 0.000748277501256778697532 1.7680941 0.5244624 0.018604103 0.018604103          9         9     0
786      WT1 0.003655478822618108813297 3.1255708 1.0753889 0.006756757 0.006756757          4         4     0
10354    KIT 0.011978874707657898024404 3.1149033 1.2396312 0.005067568 0.005067568          3         3     0
1686    KRAS 0.202970214166854234782988 1.5782177 1.2396312 0.005067568 0.005067568          2         2     0
13287   JAK2 0.715431886461324495485314 0.3512669 0.9635083 0.008445946 0.008445946          5         5     0
1372     MLL 0.781879288747886591615099 0.1890884 0.6829481 0.013513514 0.013513514          8         8     0
>         meta.results.burden.gene[meta.results.burden.gene[,"gene"] %in% fanc.genes,]
          gene          p        beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
8657     FANCL 0.06550967  1.41019676 0.7656772 0.013513514 0.013513514          4         4     0
3319     FANCI 0.28188840  0.82393519 0.7656772 0.013513514 0.013513514          8         8     0
4822    STRA13 0.30527330 -0.98279000 0.9586413 0.005669328 0.005669328          2         2     0
13379    FANCG 0.32554696 -1.49006803 1.5156478 0.003378378 0.003378378          2         2     0
6183   C1orf86 0.40149988  0.49081068 0.5850328 0.023781175 0.023781175          3         3     0
2667     FANCM 0.43434551  0.56533746 0.7231439 0.015202703 0.015202703          8         8     0
11657    FANCE 0.48768622 -1.48501695 2.1398136 0.001689189 0.001689189          1         1     0
13820    FANCC 0.48768622 -1.48501695 2.1398136 0.001689189 0.001689189          1         1     0
3988     FANCA 0.64891679 -0.34858796 0.7656772 0.013513514 0.013513514          7         7     0
9390    FANCD2 0.97327304  0.04153204 1.2396312 0.005067568 0.005067568          3         3     0
4802  C17orf70 0.97655520  0.03643125 1.2396673 0.005084863 0.005084863          3         3     0

> meta.results.burden.gene[1:50,]
                 gene                          p       beta        se  cmafTotal   cmafUsed nsnpsTotal nsnpsUsed nmiss
11384            NPM1 0.000000000000000004666874  1.0474066 0.1209304 0.17905405 0.17905405          8         8     0
2326             FLT3 0.000000000000000063908810  3.2306785 0.3865474 0.05405405 0.05405405          6         6     0
8506           DNMT3A 0.000000000000552230211637  2.3671831 0.3282382 0.06265001 0.06265001         21        21     0
3331             IDH2 0.000000000023130646061980  3.3067633 0.4946721 0.03378378 0.03378378          2         2     0
4018  ENSG00000269323 0.000000003904959983582807 -1.3197635 0.2241381 0.13682432 0.13682432          1         1     0
3539        LOC440335 0.000000171894094533876543  1.6801551 0.3214121 0.08992448 0.08992448          2         2     0
9136             IDH1 0.000000331598514265690453  3.2136150 0.6295600 0.02027027 0.02027027          3         3     0
4162             PER1 0.000000692890021435440133  1.7087936 0.3442826 0.08116412 0.08116412          6         6     0
7910           EEF1A2 0.000001137775237446866202  1.2327350 0.2533272 0.18641115 0.18641115          1         1     0
12186          TNRC18 0.000001552073217241876052  1.1277451 0.2347312 0.19070665 0.19070665         10        10     0
10515            TET2 0.000009553518252604718656  1.3067709 0.2951794 0.06251149 0.06251149         36        36     0
8096             PCNT 0.000011442460687260778087  1.1341969 0.2584800 0.13554180 0.13554180         12        12     0
8006            RUNX1 0.000011520978735816074680  1.8369797 0.4187840 0.02877979 0.02877979         13        13     0
1916             LRP1 0.000034034260333051306034  0.8104806 0.1955489 0.22651456 0.22651456          9         9     0
7803            UBE2C 0.000041498170035095508933  1.9445022 0.4743873 0.03914591 0.03914591          1         1     0
10562        KIAA1109 0.000047255810233395838932  2.2193480 0.5454551 0.02364865 0.02364865         13        13     0
7126             DDR2 0.000058857795100619874381 -1.6531321 0.4114988 0.05084746 0.05084746          2         2     0
10896            BDP1 0.000061989348153295098726  0.9370378 0.2339607 0.20975925 0.20975925          9         9     0
6850             NRAS 0.000074914710439908353693  2.5831005 0.6522797 0.01525424 0.01525424          8         8     0
6908         ADAMTSL4 0.000101601001243716289693  1.2695473 0.3266358 0.08385783 0.08385783          7         7     0
7739           CEP250 0.000111187109458178230064  3.1580161 0.8171260 0.01182432 0.01182432          7         7     0
9289           ANKMY1 0.000111187109458178230064  3.1580161 0.8171260 0.01182432 0.01182432          7         7     0
10552           USP53 0.000111187109458178230064  3.1580161 0.8171260 0.01182432 0.01182432          7         7     0
7183            ASTN1 0.000133839078173962654824  1.9595385 0.5130642 0.02372881 0.02372881         13        13     0
9389             EMC3 0.000182299524791116322091  0.9236984 0.2468224 0.28040541 0.28040541          3         3     0
2458            FARP1 0.000191424749198074471484  1.3362401 0.3582344 0.05514385 0.05514385          5         5     0
4229           MYO15A 0.000233148869783894401551  1.3709464 0.3725304 0.04229875 0.04229875         21        21     0
8017             TTC3 0.000247927027443518377750  1.1483377 0.3133774 0.09534780 0.09534780         10        10     0
13765          NOTCH1 0.000249298429951354749094  1.2113617 0.3307039 0.06834617 0.06834617         11        11     0
7781          GDAP1L1 0.000249478395604726900427 -1.2256653 0.3346257 0.07263514 0.07263514          1         1     0
1185         ARHGEF17 0.000253292668381546249353  2.6464736 0.7232961 0.01531949 0.01531949          8         8     0
11579          GPANK1 0.000314055267712844826938  1.2264781 0.3403649 0.07506803 0.07506803          5         5     0
12486            MCM7 0.000354372117981734016494  3.1471264 0.8810739 0.01013514 0.01013514          5         5     0
1999         C12orf50 0.000354372117981734287544  3.1471264 0.8810739 0.01013514 0.01013514          5         5     0
5545           SLC7A9 0.000354372117981735100696  3.1471264 0.8810739 0.01013514 0.01013514          6         6     0
10294           ARAP2 0.000354372117981735100696  3.1471264 0.8810739 0.01013514 0.01013514          6         6     0
11087         SLC23A1 0.000356907628275577887972  3.1455490 0.8810928 0.01014663 0.01014663          6         6     0
3426            METRN 0.000365365220460703169283  0.7069607 0.1983664 0.24496706 0.24496706          4         4     0
9233           SPATA3 0.000366044773861016385728  3.1398917 0.8811439 0.01018702 0.01018702          5         5     0
13155          FER1L6 0.000461698838858031076170  2.1220873 0.6059563 0.02198244 0.02198244         12        12     0
5149             TLE6 0.000473051272974335946641 -1.4972329 0.4283230 0.04250552 0.04250552          3         3     0
7943             PIGU 0.000617045303329858394782  0.8522241 0.2488965 0.26689189 0.26689189          2         2     0
436            TCF7L2 0.000617325306994607092732  1.3700787 0.4001531 0.04918521 0.04918521          4         4     0
421            SORCS3 0.000743219662802089187266  2.5827199 0.7656772 0.01351351 0.01351351          6         6     0
422            SORCS1 0.000743219662802090921989  2.5827199 0.7656772 0.01351351 0.01351351          8         8     0
4151             TP53 0.000748277501256778697532  1.7680941 0.5244624 0.01860410 0.01860410          9         9     0
7329             CD34 0.000782871117631803257325  1.2453682 0.3707798 0.05753230 0.05753230          4         4     0
1728         ADAMTS20 0.000814007272223931277640  0.8533616 0.2548879 0.15033784 0.15033784         11        11     0
3179           IGDCC4 0.000840322856677447847795  0.6763688 0.2025566 0.19018193 0.19018193          7         7     0
6721          COL24A1 0.000841237701234497006059  0.9069638 0.2716390 0.13389831 0.13389831          8         8     0



                 gene                         p              pmin rho       cmaf nmiss nsnps errflag
11384            NPM1 0.00000000000000006231578 0.000000000000000   0 0.17905405     0     8       3
2326             FLT3 0.00000000004159119762795 0.000000000000000   0 0.05405405     0     6       3
8506           DNMT3A 0.00000000065664519662610 0.000000000000000   1 0.06265001     0    21       3
4018  ENSG00000269323 0.00000000390495929759913 0.000000003904959   0 0.13682432     0     1       0
9136             IDH1 0.00000003199432113946245 0.000000336874034   1 0.02027027     0     3       0
3331             IDH2 0.00000004493313973139144 0.000000000000000   0 0.03378378     0     2       3
3539        LOC440335 0.00000007760763769242383 0.000000169958773   1 0.08992448     0     2       0
1916             LRP1 0.00000055844612457180626 0.000034033979921   1 0.22651456     0     9       0
7126             DDR2 0.00000112303401150393009 0.000058857561437   1 0.05084746     0     2       0
7910           EEF1A2 0.00000113777580490786562 0.000001137775805   0 0.18641115     0     1       0
12186          TNRC18 0.00000254387052620701859 0.000001552393842   1 0.19070665     0    10       0
1462             ST14 0.00000798384789661489150 0.000000469149787   0 0.09194148     0     4       0
10515            TET2 0.00000973448554996648209 0.000009551470734   1 0.06251149     0    36       0
4162             PER1 0.00000994446631615928415 0.000000697164039   1 0.08116412     0     6       0
6778            RNPC3 0.00003268132259280769904 0.000067241379935   0 0.25762712     0     4       0
7803            UBE2C 0.00004149812741518308111 0.000041498127415   0 0.03914591     0     1       0
6721          COL24A1 0.00004360811649046598205 0.000096512602278   0 0.13389831     0     8       0
10896            BDP1 0.00005667937916215882738 0.000061989094361   1 0.20975925     0     9       0
13765          NOTCH1 0.00006777611311859203130 0.000007700593850   0 0.06834617     0    11       0
5967           ZNF880 0.00009357194667626277389 0.000083158238226   0 0.23370192     0     7       0
7183            ASTN1 0.00009768630291826020423 0.000133838990002   1 0.02372881     0    13       0
6908         ADAMTSL4 0.00010755872574991883677 0.000101600820104   1 0.08385783     0     7       0
10562        KIAA1109 0.00010841578313424904347 0.000047255556232   1 0.02364865     0    13       0
6850             NRAS 0.00014215324891966822228 0.000074914501136   1 0.01525424     0     8       0
2458            FARP1 0.00018633779647004084152 0.000191424694611   1 0.05514385     0     5       0
7781          GDAP1L1 0.00024947839359046791513 0.000249478393590   0 0.07263514     0     1       0
10197          FGFRL1 0.00025960576293948740313 0.000241706127955   0 0.03828477     0     3       0
9389             EMC3 0.00037495772538978234785 0.000182299540965   1 0.28040541     0     3       0
4229           MYO15A 0.00046589230545476220972 0.000233148713961   1 0.04229875     0    21       0
7248             ASPM 0.00057904110521865776243 0.000424932659092   0 0.20086407     0    13       0
3426            METRN 0.00071442981815005197899 0.000365365399523   1 0.24496706     0     4       0
7739           CEP250 0.00074783202647760822421 0.000111186970999   1 0.01182432     0     7       0
9289           ANKMY1 0.00074783202647760822421 0.000111186970999   1 0.01182432     0     7       0
10552           USP53 0.00074783202647760822421 0.000111186970999   1 0.01182432     0     7       0
5149             TLE6 0.00075323556471991437965 0.000473051548984   1 0.04250552     0     3       0
3047              MGA 0.00078088636041159006224 0.000394817033164   0 0.06250000     0    14       0
8017             TTC3 0.00088817956132411686205 0.000247927107401   1 0.09534780     0    10       0
7329             CD34 0.00091171936144897422916 0.000782871711179   1 0.05753230     0     4       0
9357             CHL1 0.00091824710097400277873 0.000859782542944   0 0.05743243     0     3       0
1185         ARHGEF17 0.00094721919192385626868 0.000253292746374   1 0.01531949     0     8       0
10033          LRRIQ4 0.00108328441973140317603 0.000988094126919   1 0.03547297     0     3       0
2435           COMMD6 0.00116585873464020574157 0.001165858734640   0 0.05323194     0     1       0
11043           SEPT8 0.00120710039548990841381 0.001162522989611   1 0.07311776     0     2       0
7943             PIGU 0.00121815113374191207572 0.000617045957877   1 0.26689189     0     2       0
12486            MCM7 0.00125280761917496356699 0.000354372250721   1 0.01013514     0     5       0
13155          FER1L6 0.00127857000269924439735 0.000461699042854   1 0.02198244     0    12       0
12316          CAMK2B 0.00129517724533686996377 0.000999241709910   1 0.01302213     0     4       0
1728         ADAMTS20 0.00135014904903250606075 0.000814007843442   1 0.15033784     0    11       0
9233           SPATA3 0.00143334790237578148242 0.000366044954399   1 0.01018702     0     5       0
1999         C12orf50 0.00146533556131319894766 0.000354372250721   1 0.01013514     0     5       0


######################### non coding 0.01 repeats
                  gene                        p       beta         se  cmafTotal   cmafUsed nsnpsTotal nsnpsUsed nmiss
8586             PFN2 0.0000000000000002664605  1.9890686 0.24293636 0.12668919 0.12668919          6         5     0
5932             CD1A 0.0000000012442374240917  2.4832436 0.40880228 0.04810997 0.04810997          1         1     0
7133          PACSIN2 0.0000000045295542039248  2.8369042 0.48381665 0.03583618 0.03583618          1         1     0
8920             KLF3 0.0000020035235062399152  1.4215841 0.29908765 0.10472973 0.10472973          5         5     0
1717             FGD6 0.0000024571346728112104  0.4630778 0.09828347 0.33885708 0.33885708          8         8     0
10814            OGDH 0.0000065260074012094771  1.3335337 0.29577625 0.10810811 0.10810811          2         2     0
5042             ZNF8 0.0000086678253943043828 -1.3333883 0.29977335 0.09290541 0.09290541          1         1     0
9370            NIPBL 0.0000100892675686425552  1.8434779 0.41752502 0.04919837 0.04919837         10        10     0
357             WBP1L 0.0000106962172988840060 -1.4606236 0.33176402 0.07432432 0.07432432          5         5     0
6777  ENSR00000403778 0.0000152165703799475389  0.9719161 0.22469295 0.19661017 0.19661017          1         1     0
4556             CNN2 0.0000175535283945161972 -1.6720611 0.38940081 0.05743243 0.05743243          1         1     0
2704            INO80 0.0000226149987306957645  1.1821126 0.27897297 0.09797297 0.09797297         10        10     0
10553           TULP4 0.0000465778583229100321  0.5968424 0.14656630 0.39695946 0.39695946          6         6     0
6139            KDM5B 0.0000535444630009138749  1.6502905 0.40852910 0.04745763 0.04745763          3         3     0
9298        WDFY3-AS1 0.0000744868588388711869  1.3467120 0.33995181 0.07993197 0.07993197          1         1     0
1203  ENSR00000558810 0.0000900045617632698410  1.8173676 0.46407959 0.03938356 0.03938356          1         1     0
10958         TSC22D4 0.0000946051296121490500  1.7679752 0.45285933 0.03956129 0.03956129          2         2     0
10874           PHTF2 0.0000970987856811571758  1.3762771 0.35309713 0.06761385 0.06761385          6         6     0
10638          OSTCP1 0.0001193190917018499715 -0.9609236 0.24975136 0.17398649 0.17398649          1         1     0
336              TLX1 0.0001249485584720845467  1.3695816 0.35701441 0.05513308 0.05513308          1         1     0
8897       KCNIP4-IT1 0.0001255141308119878877  2.4144325 0.62956210 0.02043293 0.02043293          3         3     0
7404            ALMS1 0.0001488084439207132338  0.9730538 0.25653602 0.19584726 0.19584726          2         2     0
7085          FAM118A 0.0001670385597882903701 -1.1808086 0.31368889 0.06250000 0.06250000          2         2     0
680               WT1 0.0001858114946950466728  0.9312925 0.24917096 0.12500000 0.12500000          3         3     0
1373          SLCO1C1 0.0002124346118168629117  1.1325269 0.30577835 0.08360045 0.08360045          2         2     0
1179             FLI1 0.0002446207337969641903  1.4469747 0.39450463 0.05574324 0.05574324          7         7     0
5365           RNF19B 0.0002470716685476660322  0.8887876 0.24248857 0.23898305 0.23898305          5         5     0
9784             SGCD 0.0002471438600087611610  1.0068381 0.27470197 0.12842466 0.12842466         11        11     0
4614            PTPRS 0.0002477966203335093837  1.0306785 0.28125829 0.09094817 0.09094817          5         5     0
6426             IL20 0.0002776610166622033656 -1.2748651 0.35069081 0.06440678 0.06440678          1         1     0
7616             GPD2 0.0003396946912176378457  1.0984161 0.30656483 0.08569442 0.08569442          6         6     0
11108         CREB3L2 0.0003457397710789674047  0.6916436 0.19328438 0.15878378 0.15878378          9         9     0
167           HNRNPH3 0.0003509574215695936269 -1.0435742 0.29195332 0.07939189 0.07939189          6         6     0
2078           CAB39L 0.0003543721179817351007  3.1471264 0.88107388 0.01013514 0.01013514          5         5     0
9796             EBF1 0.0003685039301610075761  0.4282638 0.12024244 0.82650677 0.82650677         25        25     0
8365            ROBO1 0.0003765031993305783392  0.5064511 0.14242042 0.36995557 0.36995557         13        13     0
5593            TTLL7 0.0003805656309857651833  1.3538053 0.38100946 0.05202352 0.05202352          4         4     0
10642 ENSR00001233890 0.0004362132705023125719 -0.5808398 0.16514519 0.56925676 0.56925676          3         3     0
12111         CACNA1B 0.0005414104587663745886  0.7545655 0.21812128 0.18169838 0.18169838          6         6     0
3494             GSE1 0.0005479380695629277763  0.8455973 0.24466411 0.15046024 0.15046024         12        12     0
9548             PJA2 0.0006180456170429991705  1.0213026 0.29831519 0.09206880 0.09206880          8         8     0
3015  ENSR00000405933 0.0006615683494434497275  1.0246404 0.30091964 0.08781362 0.08781362          1         1     0
11303           CHMP7 0.0007242222624701354214  2.4298699 0.71884512 0.01196707 0.01196707          4         4     0
5027           ZNF628 0.0007687302761574642132  0.7849972 0.23336514 0.23742507 0.23742507          2         2     0
6133            PTPN7 0.0008007722223019250740  1.4241425 0.42479665 0.04830215 0.04830215          3         3     0
8095          ZNF385D 0.0009115871559858098249  2.5409567 0.76616044 0.01391300 0.01391300          4         4     0
208            DNAJC9 0.0009439504589197324430  2.5347086 0.76653249 0.01459854 0.01459854          1         1     0
7634             TBR1 0.0009510971873220742251 -1.1039775 0.33407228 0.08277027 0.08277027          2         2     0
4287            VAMP2 0.0009902905940165425587 -0.5612425 0.17042100 0.31418919 0.31418919          2         2     0
6998             MYH9 0.0010931689914459362476 -0.7305933 0.22373907 0.20101351 0.20101351          4         4     0


######################### non coding 0.01 repeats
        gene            p        beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
680      WT1 0.0001858115  0.93129252 0.2491710 0.125000000 0.125000000          3         3     0
9095    TET2 0.0299625606  1.15136274 0.5304389 0.025718858 0.025718858         11        11     0
377     SMC3 0.0405430821  3.10430839 1.5156478 0.003378378 0.003378378          1         1     0
2005    FLT3 0.0405430821  3.10430839 1.5156478 0.003378378 0.003378378          2         2     0
8973     KIT 0.0666101830  1.40459975 0.7657436 0.013565400 0.013565400          5         5     0
3684    TP53 0.1482276942  3.09378531 2.1398136 0.001689189 0.001689189          1         1     0
4776   CEBPA 0.1482276942  3.09378531 2.1398136 0.001689189 0.001689189          1         1     0
9824    NPM1 0.1841520121  1.27961512 0.9635083 0.008445946 0.008445946          4         4     0
7199  DNMT3A 0.1855514020  1.27559874 0.9635491 0.008469086 0.008469086          5         5     0
1384    KRAS 0.2029702142  1.57821767 1.2396312 0.005067568 0.005067568          3         3     0
11147   EZH2 0.2261878604 -1.50027397 1.2396528 0.005084746 0.005084746          1         1     0
11843 HNRNPK 0.2967280545  0.27443305 0.2629990 0.120013546 0.120013546         17        17     0
6847   U2AF1 0.3255469610 -1.49006803 1.5156478 0.003378378 0.003378378          2         2     0
1813  PTPN11 0.3389335569  0.43321603 0.4530247 0.025464240 0.025464240          3         3     0
11560  RAD21 0.4876862217 -1.48501695 2.1398136 0.001689189 0.001689189          1         1     0
1112     MLL 0.5063635384 -0.48053684 0.7231439 0.015202703 0.015202703          5         5     0
6097   FAM5C 0.5966678409  0.80208760 1.5156653 0.003389831 0.003389831          2         2     0
11697   JAK2 0.6527915370  0.11713622 0.2603678 0.136478376 0.136478376          4         4     0
5740    NRAS 0.7193820883  0.34618391 0.9635364 0.008474576 0.008474576          4         4     0
6822   RUNX1 0.9240786289  0.02617505 0.2746671 0.128419293 0.128419293         10        10     0



######################### non coding filtered on coverage
               gene                        p       beta        se  cmafTotal   cmafUsed nsnpsTotal nsnpsUsed nmiss
6806            PFN2 0.0000000000000001114682  2.0241415 0.2441112 0.12530921 0.12530921          6         5     0
3952          ZNF628 0.0000000000763594977791  2.2952454 0.3527021 0.17957746 0.17957746          1         1     0
4071         PLEKHG5 0.0000000003445125198729  2.1833791 0.3478216 0.14858757 0.14858757          4         2     0
2860            TOX3 0.0000000014051499963940  1.2112405 0.2000434 0.40459856 0.40459856          6         6     0
5667         PACSIN2 0.0000000054936038001572  2.8224972 0.4840100 0.03620690 0.03620690          1         1     0
8943            OGDH 0.0000001819885167122744  1.6599143 0.3181831 0.09363958 0.09363958          1         1     0
3232             HLF 0.0000001910567376793365  1.4806432 0.2843107 0.24700530 0.24700530          9         9     0
1367            FGD6 0.0000004676907084027239  0.6008102 0.1192293 0.27386780 0.27386780          7         7     0
3422 ENSR00001346374 0.0000005302364056717515  1.8926496 0.3773950 0.09950249 0.09950249          1         1     0
7113            KLF3 0.0000025701135128080424  1.4229483 0.3025948 0.10135135 0.10135135          4         4     0
2823       RAB11FIP3 0.0000036411347416446038  1.1422161 0.2466520 0.35106264 0.35106264          4         4     0
9719   MIR4668::UGCG 0.0000059064803024314693  1.1747033 0.2593330 0.25735667 0.25735667          2         2     0
2125          CTAGE5 0.0000093042918947091136  2.2217441 0.5012123 0.04587156 0.04587156          1         1     0
8917           ELFN1 0.0000105178969030470281  1.0023946 0.2274942 0.25821881 0.25821881          7         7     0
5385 ENSR00000403778 0.0000107203678988820426  1.0161314 0.2308283 0.18439716 0.18439716          1         1     0
4950 ENSR00000163508 0.0000110620078604507429  1.4339182 0.3262390 0.15422886 0.15422886          1         1     0
264            WBP1L 0.0000127148251343299580 -1.6776411 0.3843432 0.05915065 0.05915065          3         3     0
7413       WDFY3-AS1 0.0000128087409190540472  2.0501932 0.4698671 0.05066079 0.05066079          1         1     0
5651 ENSR00001041998 0.0000153456918746836024  1.3283770 0.3072339 0.19902816 0.19902816          3         3     0
3104         SMARCE1 0.0000317375752644696637  1.1276416 0.2710269 0.12575235 0.12575235          5         5     0
9173           STMN2 0.0000391947788728819934  1.3895546 0.3379123 0.09031778 0.09031778          5         4     0
2398            CHD2 0.0000581461650955096936  1.5328199 0.3812785 0.14252285 0.14252285          8         8     0
6017            GPD2 0.0000629135851035580237  2.3211096 0.5800443 0.05082363 0.05082363          3         3     0
2204           INO80 0.0000638779349663875942  1.3169498 0.3294015 0.08261883 0.08261883          6         5     0
7480           NIPBL 0.0000719068498351913536  1.7726203 0.4465166 0.04223546 0.04223546          9         8     0
8391           TULP4 0.0000731871771328772568  0.5854931 0.1476401 0.38934029 0.38934029          5         5     0
4650            CD1A 0.0000796868558910846452  2.4086200 0.6104975 0.03437500 0.03437500          1         1     0
1602            TESC 0.0000850582623921027536  1.0979859 0.2794083 0.16770021 0.16770021          3         3     0
4227          RNF19B 0.0000942900632103326456  0.9478244 0.2427310 0.23392144 0.23392144          3         3     0
7787            SGCD 0.0001056345600476162592  1.1412240 0.2943369 0.10732379 0.10732379          7         7     0
921  ENSR00000558810 0.0001077384744300534071  2.5975897 0.6707839 0.03900709 0.03900709          1         1     0
4492            ST7L 0.0001105057034761962271  3.1594783 0.8171871 0.01188747 0.01188747          4         4     0
2794            GSE1 0.0001262546609200849513  0.9839160 0.2566524 0.13690552 0.13690552         12        11     0
526              WT1 0.0001274740341488569799  0.9579426 0.2500315 0.12331654 0.12331654          2         2     0
3274           DDX42 0.0001319223475810524021  1.0344252 0.2705906 0.15066531 0.15066531          6         5     0
4805         ADIPOR1 0.0001638439756260418653  1.5537428 0.4122329 0.05239422 0.05239422          2         2     0
8424          TRIM15 0.0001679680242171599633  0.8816094 0.2342911 0.18179003 0.18179003          2         2     0
1018           CCND2 0.0001908784036524938197  0.6188038 0.1658639 0.38348364 0.38348364          8         7     0
6428           SATB1 0.0002040133432748618642  0.8688103 0.2339288 0.18704653 0.18704653         12        11     0
6937        KIAA0226 0.0002356787522383769649  0.8563891 0.2328828 0.33055865 0.33055865          7         7     0
3886           CCDC8 0.0002466912647625096595  1.3877475 0.3785794 0.10700968 0.10700968          3         3     0
8799 ENSR00000632639 0.0002480111822073928583  0.9662919 0.2637040 0.15920398 0.15920398          1         1     0
2424 ENSR00000405933 0.0002491497564689616859  2.0992192 0.5730667 0.04248366 0.04248366          1         1     0
2661          LPCAT2 0.0002760752901902625281  1.9104187 0.5253060 0.04461848 0.04461848          2         2     0
7992         C6orf62 0.0003376895230738654837  0.4669573 0.1302703 0.45821546 0.45821546         10        10     0
5505            GNAZ 0.0003459948047349956221  1.7477164 0.4884372 0.11190389 0.11190389          3         3     0
7098      KCNIP4-IT1 0.0003543721179817342875  3.1471264 0.8810739 0.01013514 0.01013514          2         1     0
3314            GRB2 0.0003691597007200729291  0.5912133 0.1660151 0.33588687 0.33588687          5         5     0
4802           KDM5B 0.0003720231489334260562  2.1581550 0.6063640 0.02608516 0.02608516          2         2     0
7034 ENSR00001378190 0.0003914465133046934474  2.5799237 0.7276031 0.02112676 0.02112676          1         1     0


######################### non coding filtered on coverage
      gene           p        beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
526     WT1 0.000127474  0.95794260 0.2500315 0.123316537 0.123316537          2         2     0
7239   TET2 0.028779893  1.43628883 0.6568930 0.019178441 0.019178441          9         8     0
7818   NPM1 0.040543082  3.10430839 1.5156478 0.003378378 0.003378378          2         2     0
1632   FLT3 0.040868247  3.09934770 1.5156741 0.003389869 0.003389869          2         2     0
1431 PTPN11 0.047016610  1.23430223 0.6214545 0.014604318 0.014604318          2         2     0
9480 HNRNPK 0.084284744  0.48249708 0.2794900 0.103299198 0.103299198         15        13     0
3801  CEBPA 0.148227694  3.09378531 2.1398136 0.001689189 0.001689189          1         1     0
2958   TP53 0.149552034  3.08375427 2.1398383 0.001700680 0.001700680          1         1     0
4767  FAM5C 0.151590107  3.06844828 2.1398759 0.001718213 0.001718213          2         1     0
1090   KRAS 0.202916790  1.57844402 1.2396624 0.005084785 0.005084785          3         3     0
843     MLL 0.203449805  1.57655863 1.2396432 0.005073294 0.005073294          3         3     0
7155    KIT 0.205670314  1.56938906 1.2400739 0.005279519 0.005279519          3         3     0
5715 DNMT3A 0.457143001  0.79972571 1.0755369 0.006834352 0.006834352          4         4     0
9242  RAD21 0.478735186 -1.51579254 2.1399273 0.001742160 0.001742160          1         1     0
5424  RUNX1 0.730590306  0.09851200 0.2860890 0.113612033 0.113612033          5         5     0
9677   JAK2 0.857395075 -0.06226895 0.3465333 0.094664396 0.094664396          2         2     0
4507   NRAS 0.976355499  0.03674339 1.2397257 0.005119613 0.005119613          3         3     0

################# coding 
        gene                          p        beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
13492   NPM1 0.000000000000000004666874  1.04740664 0.1209304 0.179054054 0.179054054          8         8     0
2743    FLT3 0.000000000000158669842946  2.43652759 0.3301671 0.065878378 0.065878378          8         8     0
10150 DNMT3A 0.000000000000626740759891  2.28376703 0.3174308 0.069471967 0.069471967         25        25     0
3983    IDH2 0.000000013054500616296249  2.32224967 0.4084620 0.047297297 0.047297297          5         5     0
12494   TET2 0.000016213366527191842669  1.20896367 0.2804029 0.070957437 0.070957437         40        40     0
9554   RUNX1 0.000047326938703021625054  1.50475272 0.3698590 0.039090165 0.039090165         15        15     0
8160    NRAS 0.000074914710439908353693  2.58310051 0.6522797 0.015254237 0.015254237          8         8     0
4953    TP53 0.000108295608479553532970  1.92641482 0.4976258 0.021982481 0.021982481         10        10     0
6607   CEBPA 0.055405577783932999369476  1.83600014 0.9584055 0.005287955 0.005287955          3         3     0
2522  PTPN11 0.067090467507398429680698  1.96910959 1.0753889 0.006756757 0.006756757          4         4     0
9606   U2AF1 0.067090467507398443558486  1.96910959 1.0753889 0.006756757 0.006756757          4         4     0
12295    KIT 0.070330937370696955390770  1.59454406 0.8810739 0.010135135 0.010135135          5         5     0
936      WT1 0.070330937370697080290860  1.59454406 0.8810739 0.010135135 0.010135135          5         5     0
15587  RAD21 0.142886245022605579135799  0.34028191 0.2322542 0.198270728 0.198270728          6         6     0
1989    KRAS 0.202970214166854234782988  1.57821767 1.2396312 0.005067568 0.005067568          2         2     0
517     SMC3 0.325546960976486610128688 -1.49006803 1.5156478 0.003378378 0.003378378          2         2     0
15966 HNRNPK 0.422434173461554707262877  0.61053173 0.7610666 0.010135135 0.010135135          6         6     0
15772   JAK2 0.506363538396016887865869 -0.48053684 0.7231439 0.015202703 0.015202703          8         8     0
8638   FAM5C 0.596667840900447998819800  0.80208760 1.5156653 0.003389831 0.003389831          2         2     0
1633     MLL 0.732588605446840124280072 -0.20569895 0.6020149 0.018581081 0.018581081         11        11     0


##### protein 5%
                 gene                          p       beta         se  cmafTotal   cmafUsed nsnpsTotal nsnpsUsed nmiss
14059            NPM1 0.000000000000000004666874  1.0474066 0.12093043 0.17905405 0.17905405          8         8     0
10574          DNMT3A 0.000000000000626740759891  2.2837670 0.31743081 0.06947197 0.06947197         25        25     0
2869             FLT3 0.000000000042940967530786  1.8258515 0.27691446 0.08952703 0.08952703         10        10     0
11366           IKZF2 0.000000000162553782222710  1.6070721 0.25137595 0.18243243 0.18243243          4         4     0
4964  ENSG00000269323 0.000000003904959983582807 -1.3197635 0.22413815 0.13682432 0.13682432          1         1     0
4141             IDH2 0.000000013054500616296249  2.3222497 0.40846204 0.04729730 0.04729730          5         5     0
9834           EEF1A2 0.000001137775237446866202  1.2327350 0.25332720 0.18641115 0.18641115          1         1     0
755        MRGPRG-AS1 0.000001620994244998847405  1.6952647 0.35349606 0.05434106 0.05434106          6         6     0
9955            RUNX1 0.000011520978735816074680  1.8369797 0.41878401 0.02877979 0.02877979         13        13     0
14247          GABBR1 0.000030668322866356776714  0.5362051 0.12863428 0.26351351 0.26351351          6         6     0
9705            UBE2C 0.000041498170035095508933  1.9445022 0.47438728 0.03914591 0.03914591          1         1     0
5248           FAM83G 0.000074104900624649342042  1.3362992 0.33721882 0.06233780 0.06233780          9         8     0
8480             NRAS 0.000074914710439908353693  2.5831005 0.65227971 0.01525424 0.01525424          8         8     0
5133             TP53 0.000108295608479553532970  1.9264148 0.49762578 0.02198248 0.02198248         10        10     0
13549        ATP6AP1L 0.000120345718131096676338  2.4209233 0.62956005 0.02027027 0.02027027          5         5     0
6916            WDR62 0.000132231614563317053735  0.6136045 0.16053427 0.28490021 0.28490021         13        13     0
3534  ENSG00000268657 0.000142026995963789453041  1.5318377 0.40262660 0.04974230 0.04974230          3         3     0
892            SWAP70 0.000166810470471911846159  1.8081867 0.48031207 0.03209459 0.03209459          4         4     0
2283             RARG 0.000168895995722435617124  0.6830158 0.18158060 0.16816741 0.16816741          6         6     0
16308        KIAA0196 0.000171163465494147894359  1.6485204 0.43865027 0.04391892 0.04391892          6         6     0
9594             PIGU 0.000173165348652006353828  0.9121882 0.24290983 0.27364865 0.27364865          3         3     0
4240            METRN 0.000181467358188834620080  0.6762444 0.18064448 0.26185895 0.26185895          6         6     0
11633            EMC3 0.000182299524791116322091  0.9236984 0.24682237 0.28040541 0.28040541          3         3     0
3949          SLC24A1 0.000185839035721821097761  0.4930339 0.13191447 0.37579767 0.37579767         12        12     0
8859             DDR2 0.000185920484066668837775 -1.3474965 0.36054218 0.06440678 0.06440678          5         4     0
6372             TLE6 0.000189855349907728799532 -1.5249470 0.40859788 0.04757308 0.04757308          5         5     0
16779           ABCA1 0.000210268124718437780015  0.4606855 0.12429628 0.40328369 0.40328369         22        22     0
16684            ROR2 0.000263981755791403032742  1.5034248 0.41208918 0.04243136 0.04243136         10        10     0
690            IFITM3 0.000277277470632597870196  0.6201126 0.17056429 0.21096161 0.21096161          5         5     0
15248            OGDH 0.000325766738098032249127  1.2447853 0.34636073 0.07095744 0.07095744          8         8     0
6284             THEG 0.000332875478878449752862  0.5912790 0.16478094 0.20949989 0.20949989         10        10     0
9515        LOC284788 0.000353468521122107771307  0.8876798 0.24846958 0.13344595 0.13344595          1         1     0
213            PCDH15 0.000359745451234783932253  0.2528272 0.07086007 1.01300444 1.01300444         35        35     0
14094           F13A1 0.000439290968136652447158 -0.5372435 0.15283094 0.31250573 0.31250573         10        10     0
4275           IGFALS 0.000550888288549952298383  0.3747050 0.10846213 0.39657298 0.39657298          5         5     0
14844         SLC18B1 0.000568241175026907826065  1.0389074 0.30145265 0.08108108 0.08108108          4         4     0
14134           ATXN1 0.000592694281715230300231  0.1935701 0.05635334 1.62772567 1.62772567         25        25     0
15388         ZNF804B 0.000657888939569081744249  0.2151153 0.06314752 0.61486486 0.61486486         19        19     0
9668          GDAP1L1 0.000692928251073881293701 -1.0309938 0.30391621 0.09459459 0.09459459          4         4     0
12037           FOXP1 0.000716162631546294215640 -1.2240824 0.36179965 0.05912162 0.05912162          7         7     0
13741            SRA1 0.000719218043990806207497  2.4310968 0.71880295 0.01195183 0.01195183          4         4     0
2294             ATF7 0.000731331122414082554842 -1.5988321 0.47336956 0.03716216 0.03716216          3         3     0
5363           UNC45B 0.000859460792559572591721 -0.6919260 0.20760483 0.23986486 0.23986486         12        12     0
13339           CDH12 0.000882773914947601189181 -0.4779997 0.14374005 0.41216216 0.41216216          9         9     0
6841           SLC7A9 0.000920033854839213418274  1.7001351 0.51303102 0.02364865 0.02364865          8         8     0
7247              BAX 0.000984494814669964338552  0.4477650 0.13589553 0.31189039 0.31189039          7         7     0
4381          SEC14L5 0.000988804896881187618540  1.6903010 0.51319328 0.02386440 0.02386440         11        11     0
16375          ZNF707 0.001036536316145084081139 -0.8043018 0.24518276 0.14020843 0.14020843          3         3     0
8465            LRIG2 0.001051771831947953842953  0.8401587 0.25643530 0.16272928 0.16272928          6         6     0
7303             VRK3 0.001073934611807498456029  1.9688328 0.60201492 0.01858108 0.01858108          7         7     0
> 



##############nw data missing done before 
> > > >                  gene                          p      beta        se  cmafTotal   cmafUsed nsnpsTotal nsnpsUsed nmiss
13234            NPM1 0.000000000000000005198732 1.0461338 0.1209554 0.17966102 0.17966102          8         8     0
2711             FLT3 0.000000000000049318060989 2.4659912 0.3273275 0.06768195 0.06768195          9         9     0
9943           DNMT3A 0.000000000000199583042110 2.3521909 0.3200663 0.06791091 0.06791091         24        24     0
3912             IDH2 0.000000013605859243950595 2.3197880 0.4085373 0.04740084 0.04740084          5         5     0
11529         SEC61A1 0.000000124700250433274801 2.1204570 0.4011103 0.06864493 0.06864493          5         5     0
1733             ST14 0.000000209348865691834268 1.9585871 0.3773161 0.05861313 0.05861313          5         5     0
704        MRGPRG-AS1 0.000000738935917701869801 1.7696535 0.3574451 0.05286128 0.05286128          6         6     0
2149             CSAD 0.000002029191699452423640 1.5548682 0.3273066 0.14261090 0.14261090          9         9     0
10401           MYO7B 0.000020847494809388570432 1.3513460 0.3175442 0.09384357 0.09384357         18        18     0
12256            TET2 0.000027671120912236711522 0.7256282 0.1731053 0.14518446 0.14518446         44        44     0
9361            RUNX1 0.000032710565942251596070 1.7704006 0.4262199 0.02740480 0.02740480         12        12     0
14781 ENSG00000257743 0.000042388014398671702702 1.9023615 0.4646635 0.02702703 0.02702703          5         5     0
6425           ZNF254 0.000057467362506624575702 1.6948940 0.4213035 0.04401626 0.04401626          8         8     0
11046            GLB1 0.000064323090770531730193 1.8035849 0.4513066 0.03961539 0.03961539          2         2     0
8207             LMNA 0.000070478854723002212204 1.8877066 0.4749349 0.03979004 0.03979004          4         4     0
2301             HELB 0.000078136748618861764687 1.8325983 0.4639437 0.03900693 0.03900693         10        10     0
1146           VPS37C 0.000084173871405202670506 1.8652174 0.4743451 0.04192747 0.04192747          5         5     0
8066            ANXA9 0.000114012786876461089873 3.1531366 0.8171595 0.01186441 0.01186441          5         5     0
4843             TP53 0.000116093623281911869239 1.9184083 0.4977408 0.02216952 0.02216952         10        10     0
12745        ATP6AP1L 0.000120345718131096676338 2.4209233 0.6295600 0.02027027 0.02027027          5         5     0
7710           ACOT11 0.000124673505871855801967 0.7948071 0.2071564 0.09713898 0.09713898         11        11     0
9968           ATRAID 0.000124818098998251202992 1.6603143 0.4327719 0.03716216 0.03716216          6         6     0
15704            ROR2 0.000160064636162976807118 1.5819575 0.4190714 0.04171823 0.04171823         10        10     0
4184            ABCC6 0.000160453098764581165905 1.1363830 0.3010839 0.06992748 0.06992748         10        10     0
14576          GIGYF1 0.000164915830796788964903 1.4362089 0.3812139 0.03723240 0.03723240          7         7     0
11997            TLR1 0.000169399407266112797272 1.8925524 0.5032370 0.02878576 0.02878576         10        10     0
15350        KIAA0196 0.000171163465494147894359 1.6485204 0.4386503 0.04391892 0.04391892          6         6     0
1752  ENSG00000268844 0.000197390549831760232519 1.9793946 0.5317618 0.05172414 0.05172414          1         1     0
7066           ZNF304 0.000228376293438168156448 1.9157186 0.5198180 0.03078272 0.03078272          2         2     0
6254             MRI1 0.000253049384171016161771 1.8811752 0.5141012 0.02547580 0.02547580          4         4     0
15169           PREX2 0.000277942757790603802402 1.1343359 0.3120564 0.05912162 0.05912162         18        18     0
11177           IP6K2 0.000281234864335923795062 1.8003365 0.4956876 0.03527690 0.03527690          4         4     0
7145            ACAP3 0.000302869423169560743076 1.1385046 0.3151279 0.07570491 0.07570491          5         5     0
12313        KIAA1109 0.000315524229182963754630 1.5873223 0.4406524 0.03558464 0.03558464         19        19     0
4604           ATP2C2 0.000315617699817664996292 1.1454037 0.3179794 0.08658198 0.08658198         19        19     0
9600          KREMEN1 0.000320169949832403081917 2.2510153 0.6255585 0.01691503 0.01691503          7         7     0
1783            PRMT8 0.000321607141685174727306 2.3789902 0.6613369 0.02583449 0.02583449          3         3     0
14805         TAS2R60 0.000362685106858874357211 2.3406188 0.6563994 0.01858108 0.01858108          1         1     0
162            RASSF4 0.000407743552660958842333 2.2096821 0.6250833 0.01751409 0.01751409          5         5     0
786             TRIM3 0.000429704007948382665609 2.3123523 0.6567065 0.01894891 0.01894891          4         4     0
15002          ENTPD4 0.000441509267988472264930 1.2170605 0.3463519 0.04095219 0.04095219         10        10     0
10342            EDAR 0.000448892920076970067863 1.5830456 0.4510700 0.02197675 0.02197675          4         4     0
11362          CRYBG3 0.000479381684615500815432 2.1165303 0.6061048 0.02211152 0.02211152         10        10     0
5904            CTDP1 0.000517015944726954227578 1.0401717 0.2996074 0.05947742 0.05947742         10        10     0
210           RHOBTB1 0.000535126076548453040796 2.0145107 0.5818034 0.02121967 0.02121967          4         4     0
10858          KLHL30 0.000562384633953766409700 2.0860473 0.6048030 0.02474246 0.02474246          8         8     0
7288           TMEM82 0.000568470616489318426694 0.8521315 0.2472650 0.06345112 0.06345112          9         9     0
2483            RAD9B 0.000571678983766462224352 1.8387324 0.5337849 0.02888874 0.02888874          4         4     0
16080           ABCA2 0.000577529196555780135330 1.4510557 0.4215791 0.05815759 0.05815759         12        12     0
15410           PYCRL 0.000609260226812426944516 3.0247655 0.8825098 0.01109613 0.01109613          5         5     0
> meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]
           gene                                p      beta        se cmafTotal  cmafUsed nsnpsTotal nsnpsUsed nmiss
439    Clinical 0.000000000000000000000005269345 1.2476189 0.1234700 0.3602883 0.3602883        142       142     0
670 FANC - ACID 0.529508831193668427772536233533 0.1219868 0.1940129 0.2019703 0.2019703         78        78     0
>         meta.results.burden.gene[meta.results.burden.gene[,"gene"] %in% clinical.genes,]
        gene                          p        beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
13234   NPM1 0.000000000000000005198732  1.04613378 0.1209554 0.179661017 0.179661017          8         8     0
2711    FLT3 0.000000000000049318060989  2.46599120 0.3273275 0.067681945 0.067681945          9         9     0
9943  DNMT3A 0.000000000000199583042110  2.35219094 0.3200663 0.067910908 0.067910908         24        24     0
3912    IDH2 0.000000013605859243950595  2.31978798 0.4085373 0.047400835 0.047400835          5         5     0
12256   TET2 0.000027671120912236711522  0.72562817 0.1731053 0.145184460 0.145184460         44        44     0
9361   RUNX1 0.000032710565942251596070  1.77040060 0.4262199 0.027404805 0.027404805         12        12     0
4843    TP53 0.000116093623281911869239  1.91840825 0.4977408 0.022169522 0.022169522         10        10     0
7980    NRAS 0.003710317430548133472296  3.12063001 1.0754138 0.006779661 0.006779661          4         4     0
12057    KIT 0.066285869450683806980429  1.11727855 0.6083786 0.028400244 0.028400244          6         6     0
928      WT1 0.070330937370697080290860  1.59454406 0.8810739 0.010135135 0.010135135          5         5     0
10687   IDH1 0.115793254911813789376218  0.88463808 0.5625040 0.022127090 0.022127090          2         2     0
1971    KRAS 0.202970214166854234782988  1.57821767 1.2396312 0.005067568 0.005067568          2         2     0
507     SMC3 0.318927922043224509884851 -1.51072866 1.5157862 0.003442433 0.003442433          2         2     0
15679 HNRNPK 0.422434173461554707262877  0.61053173 0.7610666 0.010135135 0.010135135          6         6     0
1620     MLL 0.547357757609294282019619  0.31896817 0.5300903 0.025337838 0.025337838         13        13     0
9415   U2AF1 0.594362502534595771308545  0.80712018 1.5156478 0.003378378 0.003378378          2         2     0
2501  PTPN11 0.595513595715084664838912  0.80460856 1.5156609 0.003384104 0.003384104          2         2     0
8471   FAM5C 0.596667840900447998819800  0.80208760 1.5156653 0.003389831 0.003389831          2         2     0
15485   JAK2 0.751870165974434678801686 -0.18493752 0.5849226 0.023648649 0.023648649          9         9     0
15311  RAD21 0.787526013436009808543758  0.15766029 0.5849569 0.023677673 0.023677673          5         5     0
14821   EZH2 0.962014617123815418686661  0.04196169 0.8810739 0.010135135 0.010135135          5         5     0



> > > >                  gene                          p      beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
12804            NPM1 0.000000000000000005198732 1.0461338 0.1209554 0.179661017 0.179661017          8         8     0
9610           DNMT3A 0.000000000000199583042110 2.3521909 0.3200663 0.067910908 0.067910908         24        24     0
3787             IDH2 0.000000013605859243950595 2.3197880 0.4085373 0.047400835 0.047400835          5         5     0
11155         SEC61A1 0.000000073049858231976606 2.5227418 0.4686060 0.051323584 0.051323584          3         3     0
2080             CSAD 0.000000200290463446029697 1.8664681 0.3590005 0.132360067 0.132360067          7         7     0
1675             ST14 0.000008025302852671984085 1.5505244 0.3473002 0.075642916 0.075642916          5         5     0
11855            TET2 0.000012074718120871907201 1.2367108 0.2825970 0.070860136 0.070860136         39        39     0
6289          U2AF1L4 0.000015314117990915600986 1.8070436 0.4178985 0.063863900 0.063863900          4         4     0
10051           MYO7B 0.000020847494809388570432 1.3513460 0.3175442 0.093843572 0.093843572         18        18     0
9052            RUNX1 0.000032710565942251596070 1.7704006 0.4262199 0.027404805 0.027404805         12        12     0
7642            RNPC3 0.000058432551854491600416 1.6933303 0.4213257 0.044085210 0.044085210          4         4     0
15545           ABCA2 0.000091213341086191364023 1.7045241 0.4356222 0.053713149 0.053713149         11        11     0
14058            MCM7 0.000111187109458177864145 3.1580161 0.8171260 0.011824324 0.011824324          6         6     0
4699             TP53 0.000116093623281911869239 1.9184083 0.4977408 0.022169522 0.022169522         10        10     0
1695  ENSG00000268844 0.000197390549831760232519 1.9793946 0.5317618 0.051724138 0.051724138          1         1     0
12331        ATP6AP1L 0.000238452901613662066353 2.6570861 0.7231439 0.015202703 0.015202703          4         4     0
4055            ABCC6 0.000281867651053124270918 1.5647591 0.4308947 0.049657213 0.049657213          7         7     0
11907        KIAA1109 0.000315524229182963754630 1.5873223 0.4406524 0.035584641 0.035584641         19        19     0
1727            PRMT8 0.000321607141685174727306 2.3789902 0.6613369 0.025834488 0.025834488          3         3     0
6903         C19orf33 0.000333359051844815310571 1.4745121 0.4109686 0.043071161 0.043071161          1         1     0
7806            ANXA9 0.000362147663061860734714 3.1422261 0.8811048 0.010169492 0.010169492          4         4     0
9446           CYP2D6 0.000385127214111089725794 2.2221026 0.6259319 0.017335119 0.017335119          9         9     0
12496         SLC23A1 0.000428671337539374803453 3.1037873 0.8813139 0.010658217 0.010658217          6         6     0
10995          CRYBG3 0.000479381684615500815432 2.1165303 0.6061048 0.022111523 0.022111523         10        10     0
9025             MRAP 0.000667145971961554166106 1.5671786 0.4605642 0.055685526 0.055685526          2         2     0
959             NR1H3 0.000702887484659180214637 2.4355644 0.7187841 0.011905512 0.011905512          7         7     0
3100             ESR2 0.000743219662802089295686 2.5827199 0.7656772 0.013513514 0.013513514          7         7     0
10310            NRP2 0.000765713948857463851810 2.5767097 0.7657611 0.013571927 0.013571927          8         8     0
8336            CENPF 0.000963201086109111061923 1.5206316 0.4606494 0.035593220 0.035593220         16        16     0
13565         PLEKHG1 0.001022374605895130881411 1.5774870 0.4803121 0.032094595 0.032094595         12        12     0
6653             VRK3 0.001073934611807498456029 1.9688328 0.6020149 0.018581081 0.018581081          7         7     0
2520         ATP6V0A2 0.001079284197161915715910 2.2465641 0.6872329 0.016891892 0.016891892          5         5     0
12508            SRA1 0.001103427165890638605994 2.3479072 0.7196133 0.012931967 0.012931967          4         4     0
202           RHOBTB1 0.001133542127159656255284 3.1363116 0.9635083 0.008445946 0.008445946          3         3     0
7064           TMEM82 0.001152063413686501332917 0.8161676 0.2510907 0.058181408 0.058181408          6         6     0
8137            ASTN1 0.001152957272663993407158 1.3527126 0.4161848 0.037449589 0.037449589         17        17     0
8801          L3MBTL1 0.001158574963517203845059 3.1305133 0.9635648 0.008480538 0.008480538          5         5     0
12629           NDST1 0.001236116432906644337086 1.8172006 0.5625289 0.022173643 0.022173643          7         7     0
10488          COL6A3 0.001289237272950148746262 0.8789480 0.2731039 0.100100072 0.100100072         33        33     0
9456          PACSIN2 0.001379644627650715655101 1.6620740 0.5195806 0.030434114 0.030434114          8         8     0
9933           FER1L5 0.001447720681670354924905 1.7490471 0.5491576 0.027084758 0.027084758         16        16     0
8522          TMEM234 0.002092515051375646445431 1.5683870 0.5097499 0.037848606 0.037848606          1         1     0
11530            WFS1 0.002102272226138989126565 1.0002333 0.3252378 0.066281244 0.066281244         20        20     0
14498          ENTPD4 0.002181773594185488270719 1.1121746 0.3629453 0.034195435 0.034195435          9         9     0
7615            ABCA4 0.002191838480875605600640 0.9061983 0.2958603 0.076288481 0.076288481         29        29     0
1024            OR5M3 0.002283666905670040866982 2.4928257 0.8171528 0.011841620 0.011841620          5         5     0
4434          CNTNAP4 0.002291739615674785561505 2.4919049 0.8171349 0.011830050 0.011830050          5         5     0
9635           ATRAID 0.002306314748613086298284 2.4903213 0.8171260 0.011824324 0.011824324          5         5     0
1736           GALNT8 0.002306314748613088900370 2.4903213 0.8171260 0.011824324 0.011824324          7         7     0
9707            CDKL4 0.002311738543013478613952 2.4898404 0.8171575 0.011847229 0.011847229          5         5     0
>         meta.results.burden.gene[meta.results.burden.gene[,"gene"] %in% clinical.genes,]
        gene                          p        beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
12804   NPM1 0.000000000000000005198732  1.04613378 0.1209554 0.179661017 0.179661017          8         8     0
9610  DNMT3A 0.000000000000199583042110  2.35219094 0.3200663 0.067910908 0.067910908         24        24     0
3787    IDH2 0.000000013605859243950595  2.31978798 0.4085373 0.047400835 0.047400835          5         5     0
11855   TET2 0.000012074718120871907201  1.23671084 0.2825970 0.070860136 0.070860136         39        39     0
9052   RUNX1 0.000032710565942251596070  1.77040060 0.4262199 0.027404805 0.027404805         12        12     0
4699    TP53 0.000116093623281911869239  1.91840825 0.4977408 0.022169522 0.022169522         10        10     0
7724    NRAS 0.003710317430548133472296  3.12063001 1.0754138 0.006779661 0.006779661          4         4     0
10330   IDH1 0.003760195379923354525725  3.11642276 1.0755154 0.006872852 0.006872852          1         1     0
2627    FLT3 0.007293311413177344139369  1.50846696 0.5622008 0.022073837 0.022073837          7         7     0
894      WT1 0.070330937370697080290860  1.59454406 0.8810739 0.010135135 0.010135135          5         5     0
11666    KIT 0.073757958808056370281214  1.57580826 0.8812694 0.010265529 0.010265529          5         5     0
1907    KRAS 0.202970214166854234782988  1.57821767 1.2396312 0.005067568 0.005067568          2         2     0
494     SMC3 0.318927922043224509884851 -1.51072866 1.5157862 0.003442433 0.003442433          2         2     0
15155 HNRNPK 0.422434173461554707262877  0.61053173 0.7610666 0.010135135 0.010135135          6         6     0
14968   JAK2 0.506363538396016887865869 -0.48053684 0.7231439 0.015202703 0.015202703          8         8     0
9100   U2AF1 0.594362502534595771308545  0.80712018 1.5156478 0.003378378 0.003378378          2         2     0
2418  PTPN11 0.595513595715084664838912  0.80460856 1.5156609 0.003384104 0.003384104          2         2     0
8193   FAM5C 0.596667840900447998819800  0.80208760 1.5156653 0.003389831 0.003389831          2         2     0
1566     MLL 0.732588605446840124280072 -0.20569895 0.6020149 0.018581081 0.018581081         11        11     0
14795  RAD21 0.961099484015240812517789  0.04297574 0.8811217 0.010164159 0.010164159          4         4     0
14320   EZH2 0.962014617123815418686661  0.04196169 0.8810739 0.010135135 0.010135135          5         5     0






####### filtered genotypes: pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene  & !on.x.y & !unannotated.hits & not.flat.genotype & !are.repeats &( ok.missing.filt | is.unwound.geno)   & hw.controls.ok.filt & !no.genotypes.filt   & !are.in.repeats & rare.in.controls.filt & rare.in.group

> > >                  gene                          p      beta        se  cmafTotal   cmafUsed nsnpsTotal nsnpsUsed nmiss
13351            NPM1 0.000000000000000005198732 1.0461338 0.1209554 0.17966102 0.17966102          8         8     0
10030          DNMT3A 0.000000000000199583042110 2.3521909 0.3200663 0.06791091 0.06791091         24        24     0
13594          GPANK1 0.000000011804785007255658 2.0650544 0.3621300 0.10794559 0.10794559          6         6     0
3944             IDH2 0.000000013605859243950595 2.3197880 0.4085373 0.04740084 0.04740084          5         5     0
11630         SEC61A1 0.000000124700250433274801 2.1204570 0.4011103 0.06864493 0.06864493          5         5     0
711        MRGPRG-AS1 0.000000738935917701869801 1.7696535 0.3574451 0.05286128 0.05286128          6         6     0
2166             CSAD 0.000002029191699452423640 1.5548682 0.3273066 0.14261090 0.14261090          9         9     0
7002           ZNF880 0.000004229308603069661974 1.8574131 0.4038045 0.11837622 0.11837622          8         8     0
1749             ST14 0.000006837954464606166843 1.5324597 0.3406472 0.07902129 0.07902129          6         6     0
13536           HLA-A 0.000011041783483456846593 0.5465112 0.1243287 0.28932879 0.28932879         19        19     0
6549          U2AF1L4 0.000011740494742377153553 1.7422707 0.3975652 0.06555309 0.06555309          5         5     0
2318             HELB 0.000025771180707801876196 1.9108266 0.4541005 0.04395743 0.04395743         11        11     0
6481           ZNF254 0.000027691079140785789095 0.9712906 0.2317193 0.08286761 0.08286761          9         9     0
9443            RUNX1 0.000032710565942251596070 1.7704006 0.4262199 0.02740480 0.02740480         12        12     0
14915 ENSG00000257743 0.000042388014398671702702 1.9023615 0.4646635 0.02702703 0.02702703          5         5     0
7967            RNPC3 0.000058432551854491600416 1.6933303 0.4213257 0.04408521 0.04408521          4         4     0
11141            GLB1 0.000064323090770531730193 1.8035849 0.4513066 0.03961539 0.03961539          2         2     0
8281             LMNA 0.000070478854723002212204 1.8877066 0.4749349 0.03979004 0.03979004          4         4     0
10488           MYO7B 0.000111318148685442098816 1.1076194 0.2866141 0.11425173 0.11425173         19        19     0
8139            ANXA9 0.000114012786876461089873 3.1531366 0.8171595 0.01186441 0.01186441          5         5     0
4887             TP53 0.000116093623281911869239 1.9184083 0.4977408 0.02216952 0.02216952         10        10     0
12855        ATP6AP1L 0.000120345718131096676338 2.4209233 0.6295600 0.02027027 0.02027027          5         5     0
10055          ATRAID 0.000124818098998251202992 1.6603143 0.4327719 0.03716216 0.03716216          6         6     0
10776            IDH1 0.000147522836223088713590 1.7863450 0.4706850 0.04067390 0.04067390          4         4     0
15843            ROR2 0.000160064636162976807118 1.5819575 0.4190714 0.04171823 0.04171823         10        10     0
12102            TLR1 0.000169399407266112797272 1.8925524 0.5032370 0.02878576 0.02878576         10        10     0
15486        KIAA0196 0.000171163465494147894359 1.6485204 0.4386503 0.04391892 0.04391892          6         6     0
12361            TET2 0.000185530298357044874527 0.5103515 0.1365326 0.17558987 0.17558987         45        45     0
14474            OGDH 0.000194889275992655281160 1.3032912 0.3498251 0.06953254 0.06953254          7         7     0
1768  ENSG00000268844 0.000197390549831760232519 1.9793946 0.5317618 0.05172414 0.05172414          1         1     0
7129           ZNF304 0.000228376293438168156448 1.9157186 0.5198180 0.03078272 0.03078272          2         2     0
6309             MRI1 0.000253049384171016161771 1.8811752 0.5141012 0.02547580 0.02547580          4         4     0
7210            ACAP3 0.000302869423169560743076 1.1385046 0.3151279 0.07570491 0.07570491          5         5     0
12419        KIAA1109 0.000315524229182963754630 1.5873223 0.4406524 0.03558464 0.03558464         19        19     0
4647           ATP2C2 0.000315617699817664996292 1.1454037 0.3179794 0.08658198 0.08658198         19        19     0
9686          KREMEN1 0.000320169949832403081917 2.2510153 0.6255585 0.01691503 0.01691503          7         7     0
1799            PRMT8 0.000321607141685174727306 2.3789902 0.6613369 0.02583449 0.02583449          3         3     0
2103          BCDIN3D 0.000355552591566583072857 1.3655719 0.3824006 0.06073695 0.06073695          3         3     0
14939         TAS2R60 0.000362685106858874357211 2.3406188 0.6563994 0.01858108 0.01858108          1         1     0
949             ABTB2 0.000383803694537056662526 1.3745016 0.3870770 0.06083512 0.06083512          6         6     0
794             TRIM3 0.000429704007948382665609 2.3123523 0.6567065 0.01894891 0.01894891          4         4     0
15137          ENTPD4 0.000441509267988472264930 1.2170605 0.3463519 0.04095219 0.04095219         10        10     0
10429            EDAR 0.000448892920076970067863 1.5830456 0.4510700 0.02197675 0.02197675          4         4     0
7781           ACOT11 0.000503425798686950938705 0.7018158 0.2017334 0.11438036 0.11438036         12        12     0
5956            CTDP1 0.000517015944726954227578 1.0401717 0.2996074 0.05947742 0.05947742         10        10     0
213           RHOBTB1 0.000535126076548453040796 2.0145107 0.5818034 0.02121967 0.02121967          4         4     0
10948          KLHL30 0.000562384633953766409700 2.0860473 0.6048030 0.02474246 0.02474246          8         8     0
2501            RAD9B 0.000571678983766462224352 1.8387324 0.5337849 0.02888874 0.02888874          4         4     0
16222           ABCA2 0.000577529196555780135330 1.4510557 0.4215791 0.05815759 0.05815759         12        12     0
14711             ZAN 0.000596752484487304423084 0.5176847 0.1507928 0.25183902 0.25183902         39        39     0
>         gene                          p        beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
13351   NPM1 0.000000000000000005198732  1.04613378 0.1209554 0.179661017 0.179661017          8         8     0
10030 DNMT3A 0.000000000000199583042110  2.35219094 0.3200663 0.067910908 0.067910908         24        24     0
3944    IDH2 0.000000013605859243950595  2.31978798 0.4085373 0.047400835 0.047400835          5         5     0
9443   RUNX1 0.000032710565942251596070  1.77040060 0.4262199 0.027404805 0.027404805         12        12     0
4887    TP53 0.000116093623281911869239  1.91840825 0.4977408 0.022169522 0.022169522         10        10     0
10776   IDH1 0.000147522836223088713590  1.78634504 0.4706850 0.040673902 0.040673902          4         4     0
12361   TET2 0.000185530298357044874527  0.51035154 0.1365326 0.175589866 0.175589866         45        45     0
8053    NRAS 0.007276758615167416496816  2.66958354 0.9946636 1.840112994 1.840112994          8         8     0
2731    FLT3 0.016007804546955266278285  0.93099866 0.3865090 0.045722486 0.045722486          9         9     0
12162    KIT 0.066285869450683806980429  1.11727855 0.6083786 0.028400244 0.028400244          6         6     0
936      WT1 0.070330937370697080290860  1.59454406 0.8810739 0.010135135 0.010135135          5         5     0
1987    KRAS 0.202970214166854234782988  1.57821767 1.2396312 0.005067568 0.005067568          2         2     0
512     SMC3 0.318927922043224509884851 -1.51072866 1.5157862 0.003442433 0.003442433          2         2     0
15818 HNRNPK 0.422434173461554707262877  0.61053173 0.7610666 0.010135135 0.010135135          6         6     0
1633     MLL 0.547357757609294282019619  0.31896817 0.5300903 0.025337838 0.025337838         13        13     0
9498   U2AF1 0.594362502582130414197081  0.80712018 1.5156478 1.003378378 1.003378378          4         4     0
8547   FAM5C 0.596667840900447998819800  0.80208760 1.5156653 0.003389831 0.003389831          2         2     0
2519  PTPN11 0.635146991466079557930868  0.64281004 1.3547247 0.753384104 0.753384104          4         4     0
15623   JAK2 0.751870165974434678801686 -0.18493752 0.5849226 0.023648649 0.023648649          9         9     0
15447  RAD21 0.787526013436009808543758  0.15766029 0.5849569 0.023677673 0.023677673          5         5     0
14955   EZH2 0.962014617123815418686661  0.04196169 0.8810739 0.010135135 0.010135135          5         5     0
>           gene           p         beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
10190    FANCL 0.009713958  1.630338920 0.6304851 0.021112253 0.021112253          7         7     0
5673  C17orf70 0.116968166  0.787652229 0.5024487 0.030952513 0.030952513         10        10     0
13691    FANCE 0.123163401 -1.059459207 0.6872329 0.016891892 0.016891892          2         2     0
11051   FANCD2 0.139117951 -0.768517923 0.5195906 0.030451214 0.030451214          8         8     0
6498  C19orf40 0.158692635 -1.515847751 1.0754643 0.006825939 0.006825939          1         1     0
916      FANCF 0.451435530  0.544661693 0.7232993 0.015410959 0.015410959          2         2     0
15873    FANCC 0.624420960 -0.277156188 0.5660935 0.025337838 0.025337838          4         4     0
3136     FANCM 0.656014862  0.092626079 0.2079513 0.165998622 0.165998622         23        23     0
4703     FANCA 0.750213021 -0.131985079 0.4145801 0.045684388 0.045684388         18        18     0
15727    FANCG 0.945737438  0.042848200 0.6295600 0.020270270 0.020270270          4         4     0
3930     FANCI 0.948995022 -0.030280911 0.4733696 0.037162162 0.037162162         12        12     0
7237   C1orf86 0.989379441 -0.008377649 0.6293642 0.021303642 0.021303642          2         2     0
>


########################################## BEST

####### full genotypes: pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene  & !on.x.y & !unannotated.hits & not.flat.genotype & !are.repeats &( ok.missing.filt | is.unwound.geno)   & hw.controls.ok.filt & !no.genotypes.filt   & !are.in.repeats & rare.in.controls.filt & rare.in.group
## from regular genotypes maf 0.01-10
52    Clinical 0.00000000000000000000000000002083576 0.68524895 0.06086059 0.7760321 0.7760321        188       188     0
707 FANC - ACID 0.77103190867222715088047380049829371 0.03695912 0.12699620 0.4411895 0.4411895         93        93     0

> > > >
  gene                          p       beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
12729     NPM1 0.000000000000000004666874  1.0474066 0.1209304 0.179054054 0.179054054          8         8     0
2711    FLT3 0.000000000000049318060989  2.46599120 0.3273275 0.067681945 0.067681945          9         9     0
9540    DNMT3A 0.000000000000185717732325  2.3545462 0.3199679 0.067717581 0.067717581         24        24     0
3762      IDH2 0.000000013054500616296249  2.3222497 0.4084620 0.047297297 0.047297297          5         5     0
10256     IDH1 0.000000331598514265690453  3.2136150 0.6295600 0.020270270 0.020270270          3         3     0
1660      ST14 0.000000396027727362394940  1.9358463 0.3817592 0.055818388 0.055818388          4         4     0
11779     TET2 0.000019808388425572099386  1.2116587 0.2839577 0.067567568 0.067567568         38        38     0
8986     RUNX1 0.000029515227742021809366  1.7799352 0.4261098 0.027067462 0.027067462         12        12     0
7576     RNPC3 0.000058475485235818226389  1.6931988 0.4213111 0.044067797 0.044067797          4         4     0
7658      NRAS 0.000074914710439908353693  2.5831005 0.6522797 0.015254237 0.015254237          8         8     0
4664      TP53 0.000108295608479553532970  1.9264148 0.4976258 0.021982481 0.021982481         10        10     0
13979     MCM7 0.000111187109458177864145  3.1580161 0.8171260 0.011824324 0.011824324          6         6     0
9377    CYP2D6 0.000114075565379125875849  2.3233540 0.6021360 0.018726133 0.018726133         10        10     0
12256 ATP6AP1L 0.000238452901613662066353  2.6570861 0.7231439 0.015202703 0.015202703          4         4     0
5756      TLE6 0.000258761091211884709361 -1.5154747 0.4148088 0.045872403 0.045872403          4         4     0
12948   GPANK1 0.000264697418450174023780  1.2990350 0.3561338 0.071655059 0.071655059          4         4     0
11831 KIAA1109 0.000308685940444626853371  1.5895272 0.4405686 0.035472973 0.035472973         19        19     0
7740     ANXA9 0.000362147663061860734714  3.1422261 0.8811048 0.010169492 0.010169492          4         4     0
10922   CRYBG3 0.000457354068794540470633  2.1235191 0.6059294 0.021959459 0.021959459         10        10     0
14521    PRKDC 0.000481131631454768124781  0.8484258 0.2430291 0.134599946 0.134599946         22        22     0
949      NR1H3 0.000687113422372042552747  2.4398835 0.7187395 0.011858915 0.011858915          7         7     0
3078      ESR2 0.000743219662802089295686  2.5827199 0.7656772 0.013513514 0.013513514          7         7     0
10236     NRP2 0.000752044732056882016952  2.5803454 0.7657113 0.013536653 0.013536653          8         8     0
8271     CENPF 0.000963201086109111061923  1.5206316 0.4606494 0.035593220 0.035593220         16        16     0
13488  PLEKHG1 0.001022374605895130881411  1.5774870 0.4803121 0.032094595 0.032094595         12        12     0
6595      VRK3 0.001073934611807498456029  1.9688328 0.6020149 0.018581081 0.018581081          7         7     0
2498  ATP6V0A2 0.001079284197161915715910  2.2465641 0.6872329 0.016891892 0.016891892          5         5     0
2590     RNF17 0.001099120584976548842954  1.6422868 0.5031758 0.028716216 0.028716216         11        11     0
201    RHOBTB1 0.001133542127159656255284  3.1363116 0.9635083 0.008445946 0.008445946          3         3     0
8736   L3MBTL1 0.001137640815658285363410  3.1353567 0.9635184 0.008451672 0.008451672          5         5     0
12421  SLC23A1 0.001141781308313696103804  3.1343953 0.9635286 0.008457437 0.008457437          5         5     0
6010  C19orf57 0.001141783739387275311422  3.1343871 0.9635262 0.008457398 0.008457398          4         4     0
12553    NDST1 0.001143436272736687679685  1.8291797 0.5623700 0.021959459 0.021959459          7         7     0
6488     PPP5C 0.001171382463498273291194  3.1275459 0.9635794 0.008497715 0.008497715          3         3     0
10413   COL6A3 0.001203607223687304774232  0.8839309 0.2729847 0.099679379 0.099679379         33        33     0
9387   PACSIN2 0.001366863775572994561119  1.6633733 0.5195511 0.030405405 0.030405405          8         8     0
4257     MMP15 0.001452027290110410352730 -1.4203910 0.4460881 0.042384569 0.042384569          5         5     0
8240      CD34 0.001864112754271283686758  1.1082620 0.3562319 0.059227213 0.059227213          5         5     0
8071     ASTN1 0.002088028437452790157108  1.2603555 0.4095498 0.038988816 0.038988816         17        17     0
12433     SRA1 0.002125094123910511306752  2.3383450 0.7611382 0.010227693 0.010227693          3         3     0
7549     ABCA4 0.002185309251448841645626  0.9064358 0.2958516 0.076271186 0.076271186         29        29     0
1014     OR5M3 0.002283666905670040866982  2.4928257 0.8171528 0.011841620 0.011841620          5         5     0
9635     CDKL4 0.002306314748613079793071  2.4903213 0.8171260 0.011824324 0.011824324          5         5     0
14671   NIPAL2 0.002306314748613079793071  2.4903213 0.8171260 0.011824324 0.011824324          3         3     0
6186    SLC7A9 0.002306314748613086298284  2.4903213 0.8171260 0.011824324 0.011824324          7         7     0
9565    ATRAID 0.002306314748613086298284  2.4903213 0.8171260 0.011824324 0.011824324          5         5     0
4403   CNTNAP4 0.002306314748613087165646  2.4903213 0.8171260 0.011824324 0.011824324          5         5     0
1720    GALNT8 0.002306314748613088900370  2.4903213 0.8171260 0.011824324 0.011824324          7         7     0
6702     TARM1 0.002306314748613088900370  2.4903213 0.8171260 0.011824324 0.011824324          4         4     0
6311    SAMD4B 0.002333289395956713706964  2.4875738 0.8171619 0.011847464 0.011847464          4         4     0
9294   TMPRSS6 0.002339903226242971787802  2.4868960 0.8171679 0.011853033 0.011853033          7         7     0
>         gene                          p        beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
12729   NPM1 0.000000000000000004666874  1.04740664 0.1209304 0.179054054 0.179054054          8         8     0
9540  DNMT3A 0.000000000000185717732325  2.35454622 0.3199679 0.067717581 0.067717581         24        24     0
3762    IDH2 0.000000013054500616296249  2.32224967 0.4084620 0.047297297 0.047297297          5         5     0
10256   IDH1 0.000000331598514265690453  3.21361502 0.6295600 0.020270270 0.020270270          3         3     0
11779   TET2 0.000019808388425572099386  1.21165871 0.2839577 0.067567568 0.067567568         38        38     0
8986   RUNX1 0.000029515227742021809366  1.77993519 0.4261098 0.027067462 0.027067462         12        12     0
7658    NRAS 0.000074914710439908353693  2.58310051 0.6522797 0.015254237 0.015254237          8         8     0
4664    TP53 0.000108295608479553532970  1.92641482 0.4976258 0.021982481 0.021982481         10        10     0
2606    FLT3 0.007139709251128845308998  1.51291969 0.5623700 0.021959459 0.021959459          7         7     0
2397  PTPN11 0.067090467507398429680698  1.96910959 1.0753889 0.006756757 0.006756757          4         4     0
9034   U2AF1 0.067090467507398443558486  1.96910959 1.0753889 0.006756757 0.006756757          4         4     0
11592    KIT 0.070330937370696955390770  1.59454406 0.8810739 0.010135135 0.010135135          5         5     0
886      WT1 0.070330937370697080290860  1.59454406 0.8810739 0.010135135 0.010135135          5         5     0
1891    KRAS 0.202970214166854234782988  1.57821767 1.2396312 0.005067568 0.005067568          2         2     0
492     SMC3 0.325546960976486610128688 -1.49006803 1.5156478 0.003378378 0.003378378          2         2     0
15071 HNRNPK 0.422434173461554707262877  0.61053173 0.7610666 0.010135135 0.010135135          6         6     0
14884   JAK2 0.506363538396016887865869 -0.48053684 0.7231439 0.015202703 0.015202703          8         8     0
8127   FAM5C 0.596667840900447998819800  0.80208760 1.5156653 0.003389831 0.003389831          2         2     0
1551     MLL 0.732588605446840124280072 -0.20569895 0.6020149 0.018581081 0.018581081         11        11     0
14240   EZH2 0.962014617123815418686661  0.04196169 0.8810739 0.010135135 0.010135135          5         5     0
14713  RAD21 0.962014617123815418686661  0.04196169 0.8810739 0.010135135 0.010135135          4         4     0
>           gene           p       beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
9698     FANCL 0.009701225  1.6282316 0.6295600 0.020270270 0.020270270          7         7     0
5409  C17orf70 0.136874385  1.0759111 0.7232861 0.015330052 0.015330052          7         7     0
14987    FANCG 0.325546961 -1.4900680 1.5156478 0.003378378 0.003378378          2         2     0
10524   FANCD2 0.392884700 -0.5871702 0.6872329 0.016891892 0.016891892          7         7     0
2993     FANCM 0.398249982  0.3771907 0.4465111 0.042229730 0.042229730         16        16     0
4490     FANCA 0.404252247 -0.3841600 0.4605953 0.035472973 0.035472973         17        17     0
866      FANCF 0.449841885  0.8126484 1.0753889 0.006756757 0.006756757          1         1     0
13037    FANCE 0.487686222 -1.4850169 2.1398136 0.001689189 0.001689189          1         1     0
15122    FANCC 0.624420960 -0.2771562 0.5660935 0.025337838 0.025337838          4         4     0
8451   C1orf86 0.703563876 -0.2497855 0.6564424 0.018644068 0.018644068          1         1     0
3749     FANCI 0.784054703 -0.1462422 0.5336548 0.028716216 0.028716216         11        11     0
> 

####### full genotypes: pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene  & !on.x.y & !unannotated.hits & not.flat.genotype & !are.repeats &( ok.missing.filt | is.unwound.geno)   & hw.controls.ok.filt & !no.genotypes.filt   & !are.in.repeats & rare.in.controls.filt & rare.in.group
## from regular genotypes maf 0.001-2
            gene                                           p      beta         se  cmafTotal   cmafUsed nsnpsTotal nsnpsUsed nmiss
371     Clinical 0.00000000000000000000000000000000004925784 0.9007790 0.07294287 0.50365483 0.50365483        147       147     0
1248 FANC - ACID 0.43818185292591255164396102372847963124514 0.2559758 0.33017759 0.06589567 0.06589567         36        36     0


10597     NPM1 0.000000000000000004666874 1.047407 0.1209304 0.179054054 0.179054054          8         8     0
7891    DNMT3A 0.000000000000552230211637 2.367183 0.3282382 0.062650013 0.062650013         21        21     0
3110      IDH2 0.000000000023130646061980 3.306763 0.4946721 0.033783784 0.033783784          2         2     0
8489      IDH1 0.000000331598514265690453 3.213615 0.6295600 0.020270270 0.020270270          3         3     0
5516    ZNF880 0.000001625172471612332086 1.456413 0.3037235 0.085255050 0.085255050          5         5     0
9793      TET2 0.000011451729355498889545 1.314288 0.2995343 0.059121622 0.059121622         34        34     0
9836  KIAA1109 0.000012012177824843884501 2.461700 0.5623700 0.021959459 0.021959459         12        12     0
7428     RUNX1 0.000029515227742021809366 1.779935 0.4261098 0.027067462 0.027067462         12        12     0
6332      NRAS 0.000074914710439908353693 2.583101 0.6522797 0.015254237 0.015254237          8         8     0
8621    COL6A3 0.000075092562780992657442 1.914991 0.4836391 0.035484464 0.035484464         19        19     0
3936    MYO15A 0.000084047277731249094360 1.627472 0.4138458 0.033841318 0.033841318         19        19     0
7184    CEP250 0.000111187109458178230064 3.158016 0.8171260 0.011824324 0.011824324          7         7     0
8638    ANKMY1 0.000111187109458178230064 3.158016 0.8171260 0.011824324 0.011824324          7         7     0
9827     USP53 0.000111187109458178230064 3.158016 0.8171260 0.011824324 0.011824324          7         7     0
6654     ASTN1 0.000133839078173962654824 1.959539 0.5130642 0.023728814 0.023728814         13        13     0
11632     MCM7 0.000354372117981734016494 3.147126 0.8810739 0.010135135 0.010135135          5         5     0
1867  C12orf50 0.000354372117981734287544 3.147126 0.8810739 0.010135135 0.010135135          5         5     0
5126    SLC7A9 0.000354372117981735100696 3.147126 0.8810739 0.010135135 0.010135135          6         6     0
9582     ARAP2 0.000354372117981735100696 3.147126 0.8810739 0.010135135 0.010135135          6         6     0
12254   FER1L6 0.000362685106858877338767 2.340619 0.6563994 0.018581081 0.018581081         11        11     0
8583    SPATA3 0.000366044773861016385728 3.139892 0.8811439 0.010187021 0.010187021          5         5     0
1014    PCNXL3 0.000367371992882505613368 3.139089 0.8811535 0.010192866 0.010192866          6         6     0
1646      MLL2 0.000742263683766061733646 1.769209 0.5244475 0.018581081 0.018581081         11        11     0
403     SORCS3 0.000743219662802089187266 2.582720 0.7656772 0.013513514 0.013513514          6         6     0
404     SORCS1 0.000743219662802090921989 2.582720 0.7656772 0.013513514 0.013513514          8         8     0
1825      HELB 0.000743219662802090921989 2.582720 0.7656772 0.013513514 0.013513514          5         5     0
3863      TP53 0.000748277501256778697532 1.768094 0.5244624 0.018604103 0.018604103          9         9     0
6796      CD34 0.000782871117631803257325 1.245368 0.3707798 0.057532298 0.057532298          4         4     0
11994   ENTPD4 0.000919245595398182320859 2.161655 0.6522515 0.015202703 0.015202703          5         5     0
1282   DSCAML1 0.000928345275044910908484 2.159931 0.6522738 0.015225607 0.015225607          9         9     0
4411    MYO15B 0.001079143986350689348591 1.360112 0.4160590 0.037202362 0.037202362         19        19     0
7983     THADA 0.001079284197161923739006 2.246564 0.6872329 0.016891892 0.016891892         10        10     0
4423      EVPL 0.001130514823533310337159 2.237922 0.6873530 0.016995823 0.016995823         10        10     0
67      MLLT10 0.001133542127159656688964 3.136312 0.9635083 0.008445946 0.008445946          5         5     0
328     CRTAC1 0.001133542127159656688964 3.136312 0.9635083 0.008445946 0.008445946          5         5     0
751      PAMR1 0.001133542127159656688964 3.136312 0.9635083 0.008445946 0.008445946          5         5     0
5456      VRK3 0.001133542127159656688964 3.136312 0.9635083 0.008445946 0.008445946          5         5     0
8143     SMYD1 0.001133542127159656688964 3.136312 0.9635083 0.008445946 0.008445946          4         4     0
8528     TTLL4 0.001133542127159656688964 3.136312 0.9635083 0.008445946 0.008445946          5         5     0
9059    CRYBG3 0.001133542127159656688964 3.136312 0.9635083 0.008445946 0.008445946          5         5     0
9808       CFI 0.001133542127159656688964 3.136312 0.9635083 0.008445946 0.008445946          5         5     0
12106      TOX 0.001133542127159656688964 3.136312 0.9635083 0.008445946 0.008445946          5         5     0
10885   UNC5CL 0.001137640815658282978165 3.135357 0.9635184 0.008451672 0.008451672          5         5     0
10323  SLC23A1 0.001141781308313696103804 3.134395 0.9635286 0.008457437 0.008457437          5         5     0
2971    IGDCC4 0.001145964219784313762457 3.133427 0.9635388 0.008463241 0.008463241          5         5     0
8472      NRP2 0.001150190176437126780840 3.132453 0.9635491 0.008469086 0.008469086          5         5     0
6073    RAD54L 0.001154481159148085288660 3.131391 0.9635364 0.008474576 0.008474576          5         5     0
6279    CELSR2 0.001154595264666994063746 2.233908 0.6873848 0.017042346 0.017042346         10        10     0
6732      PKP1 0.001158678519010327128552 3.130429 0.9635465 0.008480341 0.008480341          5         5     0
4758   TMPRSS9 0.001162843638144572117055 3.129539 0.9635751 0.008486381 0.008486381          5         5     0
>         gene                          p       beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
10597   NPM1 0.000000000000000004666874  1.0474066 0.1209304 0.179054054 0.179054054          8         8     0
7891  DNMT3A 0.000000000000552230211637  2.3671831 0.3282382 0.062650013 0.062650013         21        21     0
3110    IDH2 0.000000000023130646061980  3.3067633 0.4946721 0.033783784 0.033783784          2         2     0
8489    IDH1 0.000000331598514265690453  3.2136150 0.6295600 0.020270270 0.020270270          3         3     0
9793    TET2 0.000011451729355498889545  1.3142885 0.2995343 0.059121622 0.059121622         34        34     0
7428   RUNX1 0.000029515227742021809366  1.7799352 0.4261098 0.027067462 0.027067462         12        12     0
6332    NRAS 0.000074914710439908353693  2.5831005 0.6522797 0.015254237 0.015254237          8         8     0
3863    TP53 0.000748277501256778697532  1.7680941 0.5244624 0.018604103 0.018604103          9         9     0
2164    FLT3 0.002032749790892977338752  2.3481990 0.7610666 0.010135135 0.010135135          5         5     0
736      WT1 0.003655478822618108813297  3.1255708 1.0753889 0.006756757 0.006756757          4         4     0
9639     KIT 0.011978874707657898024404  3.1149033 1.2396312 0.005067568 0.005067568          3         3     0
12229  RAD21 0.040543082123378960945903  3.1043084 1.5156478 0.003378378 0.003378378          2         2     0
1980  PTPN11 0.067090467507398429680698  1.9691096 1.0753889 0.006756757 0.006756757          4         4     0
7470   U2AF1 0.067090467507398443558486  1.9691096 1.0753889 0.006756757 0.006756757          4         4     0
1573    KRAS 0.202970214166854234782988  1.5782177 1.2396312 0.005067568 0.005067568          2         2     0
12530 HNRNPK 0.422434173461554707262877  0.6105317 0.7610666 0.010135135 0.010135135          6         6     0
11839   EZH2 0.449841885210977288078738  0.8126484 1.0753889 0.006756757 0.006756757          4         4     0
408     SMC3 0.487686221714785539393944 -1.4850169 2.1398136 0.001689189 0.001689189          1         1     0
6702   FAM5C 0.596667840900447998819800  0.8020876 1.5156653 0.003389831 0.003389831          2         2     0
12373   JAK2 0.715431886461324495485314  0.3512669 0.9635083 0.008445946 0.008445946          5         5     0
1292     MLL 0.781879288747886591615099  0.1890884 0.6829481 0.013513514 0.013513514          8         8     0
>           gene          p        beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
8032     FANCL 0.06709047  1.96910959 1.0753889 0.006756757 0.006756757          3         3     0
3098     FANCI 0.28188840  0.82393519 0.7656772 0.013513514 0.013513514          8         8     0
12462    FANCG 0.32554696 -1.49006803 1.5156478 0.003378378 0.003378378          2         2     0
2482     FANCM 0.43434551  0.56533746 0.7231439 0.015202703 0.015202703          8         8     0
10852    FANCE 0.48768622 -1.48501695 2.1398136 0.001689189 0.001689189          1         1     0
12879    FANCC 0.48768622 -1.48501695 2.1398136 0.001689189 0.001689189          1         1     0
3719     FANCA 0.64891679 -0.34858796 0.7656772 0.013513514 0.013513514          7         7     0
8725    FANCD2 0.97327304  0.04153204 1.2396312 0.005067568 0.005067568          3         3     0
4472  C17orf70 0.97655520  0.03643125 1.2396673 0.005084863 0.005084863          3         3     0



          gene          p        beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
1233       ATM 0.01923871  1.32515302 0.5660935 0.025337838 0.025337838         15        15     0
3941     TOP3A 0.02192947  2.20796334 0.9635083 0.008445946 0.008445946          5         5     0
8032     FANCL 0.06709047  1.96910959 1.0753889 0.006756757 0.006756757          3         3     0
2791      FAN1 0.15753607  1.15493162 0.8171260 0.011824324 0.011824324          7         7     0
3098     FANCI 0.28188840  0.82393519 0.7656772 0.013513514 0.013513514          8         8     0
3119       BLM 0.32554696 -1.49006803 1.5156478 0.003378378 0.003378378          2         2     0
12462    FANCG 0.32554696 -1.49006803 1.5156478 0.003378378 0.003378378          2         2     0
2482     FANCM 0.43434551  0.56533746 0.7231439 0.015202703 0.015202703          8         8     0
10852    FANCE 0.48768622 -1.48501695 2.1398136 0.001689189 0.001689189          1         1     0
12879    FANCC 0.48768622 -1.48501695 2.1398136 0.001689189 0.001689189          1         1     0
1189    MRE11A 0.54921451 -0.57708133 0.9635083 0.008445946 0.008445946          5         5     0
3064       FAH 0.59551360  0.80460856 1.5156609 0.003384104 0.003384104          2         2     0
3719     FANCA 0.64891679 -0.34858796 0.7656772 0.013513514 0.013513514          7         7     0
8725    FANCD2 0.97327304  0.04153204 1.2396312 0.005067568 0.005067568          3         3     0
4472  C17orf70 0.97655520  0.03643125 1.2396673 0.005084863 0.005084863          3         3     0
> sum(meta.results.burden[,"gene"] %in% fanc.cluster)
Error in meta.results.burden[, "gene"] %in% fanc.cluster : 
  error in evaluating the argument 'table' in selecting a method for function '%in%': Error: object 'fanc.cluster' not found
> meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]
           gene                                           p      beta         se cmafTotal  cmafUsed nsnpsTotal nsnpsUsed nmiss
269    Clinical 0.00000000000000000000000000000000004925784 0.9007790 0.07294287 0.5036548 0.5036548        147       147     0
425 FANC - ACID 0.01778860926305391235158204210620169760659 0.5807776 0.24505496 0.1267122 0.1267122         72        72     0


### real genos 0.01 no filter tests
## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene  & !on.x.y & !unannotated.hits & not.flat.genotype & !are.repeats & ( ok.missing | is.unwound.geno) & hw.controls.ok & !no.genotypes & !are.in.repeats & rare.in.controls & rare.in.group
> > > >           gene                          p      beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
13158     NPM1 0.000000000000000004666874 1.0474066 0.1209304 0.179054054 0.179054054          8         8     0
9885    DNMT3A 0.000000000000626740759891 2.2837670 0.3174308 0.069471967 0.069471967         25        25     0
3871      IDH2 0.000000013054500616296249 2.3222497 0.4084620 0.047297297 0.047297297          5         5     0
10623     IDH1 0.000000331598514265690453 3.2136150 0.6295600 0.020270270 0.020270270          3         3     0
1711      ST14 0.000000398444428453536362 1.8902375 0.3728500 0.059389817 0.059389817          5         5     0
9307     RUNX1 0.000011520978735816074680 1.8369797 0.4187840 0.028779791 0.028779791         13        13     0
12182     TET2 0.000016213366527191842669 1.2089637 0.2804029 0.070957437 0.070957437         40        40     0
7867     RNPC3 0.000017560009932915971728 1.6933256 0.3943606 0.047457627 0.047457627          6         6     0
9081     UBE2C 0.000041498170035095508933 1.9445022 0.4743873 0.039145907 0.039145907          1         1     0
7951      NRAS 0.000074914710439908353693 2.5831005 0.6522797 0.015254237 0.015254237          8         8     0
4813      TP53 0.000108295608479553532970 1.9264148 0.4976258 0.021982481 0.021982481         10        10     0
14448     MCM7 0.000111187109458177864145 3.1580161 0.8171260 0.011824324 0.011824324          6         6     0
15951   NOTCH1 0.000142414483804481471372 1.0743503 0.2824312 0.092344953 0.092344953         17        17     0
12676 ATP6AP1L 0.000238452901613662066353 2.6570861 0.7231439 0.015202703 0.015202703          4         4     0
12845  SLC23A1 0.000356907628275577887972 3.1455490 0.8810928 0.010146626 0.010146626          6         6     0
8034     ANXA9 0.000362147663061860734714 3.1422261 0.8811048 0.010169492 0.010169492          4         4     0
516     TCF7L2 0.000617325306994607092732 1.3700787 0.4001531 0.049185215 0.049185215          4         4     0
12235 KIAA1109 0.000664814726269654884296 1.4730227 0.4327719 0.037162162 0.037162162         20        20     0
984      NR1H3 0.000687113422372042552747 2.4398835 0.7187395 0.011858915 0.011858915          7         7     0
12857     SRA1 0.000719218043990806207497 2.4310968 0.7188029 0.011951831 0.011951831          4         4     0
3162      ESR2 0.000743219662802089295686 2.5827199 0.7656772 0.013513514 0.013513514          7         7     0
1772    GALNT8 0.000743219662802090921989 2.5827199 0.7656772 0.013513514 0.013513514          8         8     0
10603     NRP2 0.000752044732056882016952 2.5803454 0.7657113 0.013536653 0.013536653          8         8     0
8573     CENPF 0.000963201086109111061923 1.5206316 0.4606494 0.035593220 0.035593220         16        16     0
4099   SEC14L5 0.000988804896881187618540 1.6903010 0.5131933 0.023864397 0.023864397         11        11     0
14243   CAMK2B 0.000999240785028693507056 2.3601313 0.7172037 0.013022127 0.013022127          4         4     0
13941  PLEKHG1 0.001022374605895130881411 1.5774870 0.4803121 0.032094595 0.032094595         12        12     0
7180     NPHP4 0.001059862888942325089003 1.2595164 0.3846869 0.059876832 0.059876832         24        24     0
6846      VRK3 0.001073934611807498456029 1.9688328 0.6020149 0.018581081 0.018581081          7         7     0
2572  ATP6V0A2 0.001079284197161915715910 2.2465641 0.6872329 0.016891892 0.016891892          5         5     0
205    RHOBTB1 0.001133542127159656255284 3.1363116 0.9635083 0.008445946 0.008445946          3         3     0
9042   L3MBTL1 0.001137640815658285363410 3.1353567 0.9635184 0.008451672 0.008451672          5         5     0
12978    NDST1 0.001143436272736687679685 1.8291797 0.5623700 0.021959459 0.021959459          7         7     0
9771      MLC1 0.001145907487241422167620 3.1334404 0.9635387 0.008463163 0.008463163          5         5     0
6729     PPP5C 0.001171382463498273291194 3.1275459 0.9635794 0.008497715 0.008497715          3         3     0
11808   FGFRL1 0.001196865201944575943585 1.5303958 0.4723989 0.038284767 0.038284767          3         3     0
9716    CYP2D6 0.001222570793866454858662 1.8189709 0.5625286 0.022162559 0.022162559         11        11     0
12299  TMEM154 0.001234952964252972055018 3.1137951 0.9638198 0.008622428 0.008622428          4         4     0
272     DNAJC9 0.001325403525163939633461 2.2092866 0.6881609 0.017976919 0.017976919          2         2     0
9726   PACSIN2 0.001366863775572994561119 1.6633733 0.5195511 0.030405405 0.030405405          8         8     0
11302   CRYBG3 0.001406269910620657143455 1.8678690 0.5849226 0.023648649 0.023648649         11        11     0
2272      HELB 0.001419332295452793235926 1.7519940 0.5490956 0.027027027 0.027027027         10        10     0
10784   COL6A3 0.001723090950845458397406 0.8514486 0.2716617 0.101368568 0.101368568         34        34     0
14895   ENTPD4 0.002006604281164279298538 1.1211255 0.3629115 0.033789510 0.033789510          9         9     0
4753     RPAIN 0.002059928109132580942991 2.3454000 0.7611344 0.010199924 0.010199924          6         6     0
8369     ASTN1 0.002088028437452790157108 1.2603555 0.4095498 0.038988816 0.038988816         17        17     0
1049     OR5M3 0.002283666905670040866982 2.4928257 0.8171528 0.011841620 0.011841620          5         5     0
9983     CDKL4 0.002306314748613079793071 2.4903213 0.8171260 0.011824324 0.011824324          5         5     0
15158   NIPAL2 0.002306314748613079793071 2.4903213 0.8171260 0.011824324 0.011824324          3         3     0
6415    SLC7A9 0.002306314748613086298284 2.4903213 0.8171260 0.011824324 0.011824324          7         7     0
>         gene                          p        beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
13158   NPM1 0.000000000000000004666874  1.04740664 0.1209304 0.179054054 0.179054054          8         8     0
9885  DNMT3A 0.000000000000626740759891  2.28376703 0.3174308 0.069471967 0.069471967         25        25     0
3871    IDH2 0.000000013054500616296249  2.32224967 0.4084620 0.047297297 0.047297297          5         5     0
10623   IDH1 0.000000331598514265690453  3.21361502 0.6295600 0.020270270 0.020270270          3         3     0
9307   RUNX1 0.000011520978735816074680  1.83697974 0.4187840 0.028779791 0.028779791         13        13     0
12182   TET2 0.000016213366527191842669  1.20896367 0.2804029 0.070957437 0.070957437         40        40     0
7951    NRAS 0.000074914710439908353693  2.58310051 0.6522797 0.015254237 0.015254237          8         8     0
4813    TP53 0.000108295608479553532970  1.92641482 0.4976258 0.021982481 0.021982481         10        10     0
2678    FLT3 0.007139709251128845308998  1.51291969 0.5623700 0.021959459 0.021959459          7         7     0
6423   CEBPA 0.055405577783932999369476  1.83600014 0.9584055 0.005287955 0.005287955          3         3     0
2468  PTPN11 0.067090467507398429680698  1.96910959 1.0753889 0.006756757 0.006756757          4         4     0
9357   U2AF1 0.067090467507398443558486  1.96910959 1.0753889 0.006756757 0.006756757          4         4     0
11989    KIT 0.070330937370696955390770  1.59454406 0.8810739 0.010135135 0.010135135          5         5     0
914      WT1 0.070330937370697080290860  1.59454406 0.8810739 0.010135135 0.010135135          5         5     0
1946    KRAS 0.202970214166854234782988  1.57821767 1.2396312 0.005067568 0.005067568          2         2     0
505     SMC3 0.325546960976486610128688 -1.49006803 1.5156478 0.003378378 0.003378378          2         2     0
15570 HNRNPK 0.422434173461554707262877  0.61053173 0.7610666 0.010135135 0.010135135          6         6     0
15377   JAK2 0.506363538396016887865869 -0.48053684 0.7231439 0.015202703 0.015202703          8         8     0
8426   FAM5C 0.596667840900447998819800  0.80208760 1.5156653 0.003389831 0.003389831          2         2     0
1599     MLL 0.732588605446840124280072 -0.20569895 0.6020149 0.018581081 0.018581081         11        11     0
14714   EZH2 0.962014617123815418686661  0.04196169 0.8810739 0.010135135 0.010135135          5         5     0
15201  RAD21 0.962014617123815418686661  0.04196169 0.8810739 0.010135135 0.010135135          4         4     0
>           gene           p       beta        se   cmafTotal    cmafUsed nsnpsTotal nsnpsUsed nmiss
10050    FANCL 0.009701225  1.6282316 0.6295600 0.020270270 0.020270270          7         7     0
5589  C17orf70 0.053156405  1.2173415 0.6295569 0.020461234 0.020461234         10        10     0
5611    STRA13 0.305273296 -0.9827900 0.9586413 0.005669328 0.005669328          2         2     0
15483    FANCG 0.325546961 -1.4900680 1.5156478 0.003378378 0.003378378          2         2     0
10894   FANCD2 0.392884700 -0.5871702 0.6872329 0.016891892 0.016891892          7         7     0
3075     FANCM 0.398249982  0.3771907 0.4465111 0.042229730 0.042229730         16        16     0
4627     FANCA 0.404252247 -0.3841600 0.4605953 0.035472973 0.035472973         17        17     0
894      FANCF 0.449841885  0.8126484 1.0753889 0.006756757 0.006756757          1         1     0
6417  C19orf40 0.487686222 -1.4850169 2.1398136 0.001689189 0.001689189          1         1     0
13473    FANCE 0.487686222 -1.4850169 2.1398136 0.001689189 0.001689189          1         1     0
15622    FANCC 0.624420960 -0.2771562 0.5660935 0.025337838 0.025337838          4         4     0
7155   C1orf86 0.673439727  0.2314917 0.5492980 0.027254067 0.027254067          5         5     0
3858     FANCI 0.784054703 -0.1462422 0.5336548 0.028716216 0.028716216         11        11     0
## the.an<-
##   target<-"1"
## use.samples<-the.samples[pheno[,"SampleProject"]==target]
## the.an<-gsub(".GT",".AD",use.samples)
## use.samples
## snpinfo[1:5,]
## genos<-a.indel[loci,use.samples]
## a.indel[loci,the.an]
## colnames(genos)[genos[2,]!="0/0"]
## snpinfo[1:5,]
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


## "51.GT"  "52.GT"  "53.GT"  "54.GT"  "55.GT"  "56.GT"  "57.GT"  "58.GT"  "59.GT"  "60.GT"  "61.GT"  "62.GT"  "63.GT"  "64.GT"  "65.GT"  "66.GT"  "67.GT"  "68.GT"  "69.GT"  "70.GT"  "71.GT"  "72.GT"  "73.GT"  "74.GT"  "75.GT"  "76.GT"  "77.GT" "78.GT"  "79.GT"  "80.GT"  "81.GT"  "82.GT"  "83.GT"  "84.GT"  "85.GT"  "86.GT"  "87.GT"  "88.GT"  "89.GT"  "91.GT"  "92.GT"  "93.GT"  "94.GT"  "96.GT"  "97.GT"  "98.GT"  "99.GT"  "100.GT

http://www.hindawi.com/journals/ijcb/2012/161837/
http://www.genome.jp/dbget-bin/get_linkdb?-t+8+path:hsa00020
KEGG GENES

hsa:1431             CS; citrate synthase (EC:2.3.3.1); K01647 citrate synthase [EC:2.3.3.1] 
hsa:1737             DLAT, DLTA, PDC-E2, PDCE2; dihydrolipoamide S-acetyltransferase (EC:2.3.1.12); K00627 pyruvate dehyd 
hsa:1738             DLD, DLDD, DLDH, E3, GCSL, LAD, PHE3; dihydrolipoamide dehydrogenase (EC:1.8.1.4); K00382 dihydrolip 
hsa:1743             DLST, DLTS; dihydrolipoamide S-succinyltransferase (E2 component of 2-oxo-glutarate complex) (EC:2.3 
hsa:2271             FH, HLRCC, LRCC, MCL, MCUL1; fumarate hydratase (EC:4.2.1.2); K01679 fumarate hydratase, class II [E 
hsa:3417             IDH1, HEL-216, HEL-S-26, IDCD, IDH, IDP, IDPC, PICD; isocitrate dehydrogenase 1 (NADP+), soluble (EC 
hsa:3418             IDH2, D2HGA2, ICD-M, IDH, IDHM, IDP, IDPM, mNADP-IDH; isocitrate dehydrogenase 2 (NADP+), mitochondr 
hsa:3419             IDH3A; isocitrate dehydrogenase 3 (NAD+) alpha (EC:1.1.1.41); K00030 isocitrate dehydrogenase (NAD+) 
hsa:3420             IDH3B, H-IDHB, RP46; isocitrate dehydrogenase 3 (NAD+) beta (EC:1.1.1.41); K00030 isocitrate dehydro 
hsa:3421             IDH3G, H-IDHG; isocitrate dehydrogenase 3 (NAD+) gamma (EC:1.1.1.41); K00030 isocitrate dehydrogenas 
hsa:4190             MDH1, HEL-S-32, MDH-s, MDHA, MGC:1375, MOR2; malate dehydrogenase 1, NAD (soluble) (EC:1.1.1.37 1.1. 
hsa:4191             MDH2, M-MDH, MDH, MGC:3559, MOR1; malate dehydrogenase 2, NAD (mitochondrial) (EC:1.1.1.37); K00026  
hsa:47               ACLY, ACL, ATPCL, CLATP; ATP citrate lyase (EC:2.3.3.8); K01648 ATP citrate (pro-S)-lyase [EC:2.3.3. 
hsa:48               ACO1, ACONS, HEL60, IREB1, IREBP, IREBP1, IRP1; aconitase 1, soluble (EC:4.2.1.3); K01681 aconitate  
hsa:4967             OGDH, AKGDH, E1k, OGDC; oxoglutarate (alpha-ketoglutarate) dehydrogenase (lipoamide) (EC:1.2.4.2); K 
hsa:50               ACO2, ACONM, ICRD; aconitase 2, mitochondrial (EC:4.2.1.3); K01681 aconitate hydratase [EC:4.2.1.3] 
hsa:5091             PC, PCB; pyruvate carboxylase (EC:6.4.1.1); K01958 pyruvate carboxylase [EC:6.4.1.1] 
hsa:5105             PCK1, PEPCK-C, PEPCK1, PEPCKC; phosphoenolpyruvate carboxykinase 1 (soluble) (EC:4.1.1.32); K01596 p 
hsa:5106             PCK2, PEPCK, PEPCK-M, PEPCK2; phosphoenolpyruvate carboxykinase 2 (mitochondrial) (EC:4.1.1.32); K01 
hsa:5160             PDHA1, PDHA, PDHCE1A, PHE1A; pyruvate dehydrogenase (lipoamide) alpha 1 (EC:1.2.4.1); K00161 pyruvat 
hsa:5161             PDHA2, PDHAL; pyruvate dehydrogenase (lipoamide) alpha 2 (EC:1.2.4.1); K00161 pyruvate dehydrogenase 
hsa:5162             PDHB, PDHBD, PDHE1-B, PHE1B; pyruvate dehydrogenase (lipoamide) beta (EC:1.2.4.1); K00162 pyruvate d 
hsa:55753            OGDHL; oxoglutarate dehydrogenase-like (EC:1.2.4.-); K00164 2-oxoglutarate dehydrogenase E1 componen 
hsa:6389             SDHA, CMD1GG, FP, PGL5, SDH1, SDH2, SDHF; succinate dehydrogenase complex, subunit A, flavoprotein ( 
hsa:6390             SDHB, CWS2, IP, PGL4, SDH, SDH1, SDH2, SDHIP; succinate dehydrogenase complex, subunit B, iron sulfu 
hsa:6391             SDHC, CYB560, CYBL, PGL3, QPS1, SDH3; succinate dehydrogenase complex, subunit C, integral membrane  
hsa:6392             SDHD, CBT1, CII-4, CWS3, PGL, PGL1, QPs3, SDH4, cybS; succinate dehydrogenase complex, subunit D, in 
hsa:8801             SUCLG2, GBETA; succinate-CoA ligase, GDP-forming, beta subunit (EC:6.2.1.4); K01900 succinyl-CoA syn 
hsa:8802             SUCLG1, GALPHA, MTDPS9, SUCLA1; succinate-CoA ligase, alpha subunit (EC:6.2.1.4 6.2.1.5); K01899 suc 
hsa:8803             SUCLA2, A-BETA, MTDPS5, SCS-betaA; succinate-CoA ligase, ADP-forming, beta subunit (EC:6.2.1.5); K01 


 http://www.genome.jp/dbget-bin/get_linkdb?-t+genes+path:hsa04666
 Fc gamma R-mediated
ID                   Definition
----------------------------------------------------------------------------------------------------
hsa:10000            AKT3, MPPH, PKB-GAMMA, PKBG, PRKBG, RAC-PK-gamma, RAC-gamma, STK-2; v-akt murine thymoma viral oncog 
hsa:100137049        PLA2G4B, HsT16992, cPLA2-beta; phospholipase A2, group IVB (cytosolic) (EC:3.1.1.4); K16342 cytosoli 
hsa:10092            ARPC5, ARC16, dJ127C7.3, p16-Arc; actin related protein 2/3 complex, subunit 5, 16kDa; K05754 actin  
hsa:10093            ARPC4, ARC20, P20-ARC; actin related protein 2/3 complex, subunit 4, 20kDa; K05755 actin related pro 
hsa:10094            ARPC3, ARC21, p21-Arc; actin related protein 2/3 complex, subunit 3, 21kDa; K05756 actin related pro 
hsa:10095            ARPC1B, ARC41, p40-ARC, p41-ARC; actin related protein 2/3 complex, subunit 1B, 41kDa; K05757 actin  
hsa:10109            ARPC2, ARC34, PNAS-139, p34-Arc; actin related protein 2/3 complex, subunit 2, 34kDa; K05758 actin r 
hsa:10163            WASF2, IMD2, SCAR2, WASF4, WAVE2, dJ393P12.2; WAS protein family, member 2; K05748 WAS protein famil 
hsa:10451            VAV3; vav 3 guanine nucleotide exchange factor; K05730 guanine nucleotide exchange factor VAV 
hsa:10552            ARPC1A, Arc40, HEL-68, SOP2Hs, SOP2L; actin related protein 2/3 complex, subunit 1A, 41kDa; K05757 a 
hsa:1072             CFL1, CFL, HEL-S-15; cofilin 1 (non-muscle); K05765 cofilin 
hsa:1073             CFL2, NEM7; cofilin 2 (muscle); K05765 cofilin 
hsa:10810            WASF3, Brush-1, SCAR3, WAVE3; WAS protein family, member 3; K06083 WAS protein family, member 3 
hsa:123745           PLA2G4E; phospholipase A2, group IVE (EC:3.1.1.4); K16342 cytosolic phospholipase A2 [EC:3.1.1.4] 
hsa:1398             CRK, CRKII, p38; v-crk avian sarcoma virus CT10 oncogene homolog; K04438 proto-oncogene C-crk 
hsa:1399             CRKL; v-crk avian sarcoma virus CT10 oncogene homolog-like; K04438 proto-oncogene C-crk 
hsa:1785             DNM2, CMT2M, CMTDI1, CMTDIB, DI-CMTB, DYN2, DYNII, LCCS5; dynamin 2 (EC:3.6.5.5); K01528 dynamin GTP 
hsa:1794             DOCK2; dedicator of cytokinesis 2; K12367 dedicator of cytokinesis protein 2 
hsa:207              AKT1, AKT, CWS6, PKB, PKB-ALPHA, PRKBA, RAC, RAC-ALPHA; v-akt murine thymoma viral oncogene homolog  
hsa:208              AKT2, HIHGHH, PKBB, PKBBETA, PRKBB, RAC-BETA; v-akt murine thymoma viral oncogene homolog 2 (EC:2.7. 
hsa:2209             FCGR1A, CD64, CD64A, FCRI, IGFR1; Fc fragment of IgG, high affinity Ia, receptor (CD64); K06498 high 
hsa:2212             FCGR2A, CD32, CD32A, CDw32, FCG2, FCGR2, FCGR2A1, FcGR, IGFR2; Fc fragment of IgG, low affinity IIa, 
hsa:2213             FCGR2B, CD32, CD32B, FCG2, FCGR2, IGFR2; Fc fragment of IgG, low affinity IIb, receptor (CD32); K125 
hsa:2214             FCGR3A, CD16, CD16A, FCG3, FCGR3, FCGRIII, FCR-10, FCRIII, FCRIIIA, IGFR3, IMD20; Fc fragment of IgG 
hsa:23396            PIP5K1C, LCCS3, PIP5K-GAMMA, PIP5Kgamma, PIPKIg_v4; phosphatidylinositol-4-phosphate 5-kinase, type  
hsa:23533            PIK3R5, F730038I15Rik, FOAP-2, P101-PI3K, p101; phosphoinositide-3-kinase, regulatory subunit 5; K02 
hsa:255189           PLA2G4F, PLA2G4FZ; phospholipase A2, group IVF (EC:3.1.1.4); K16342 cytosolic phospholipase A2 [EC:3 
hsa:27040            LAT, LAT1, pp36; linker for activation of T cells; K07362 linker for activation of T cells 
hsa:273              AMPH, AMPH1; amphiphysin; K12562 amphiphysin 
hsa:283748           PLA2G4D, cPLA2delta; phospholipase A2, group IVD (cytosolic) (EC:3.1.1.4); K16342 cytosolic phosphol 
hsa:2934             GSN, ADF, AGEL; gelsolin; K05768 gelsolin 
hsa:3055             HCK, JTK9, p59Hck, p61Hck; hemopoietic cell kinase (EC:2.7.10.2); K08893 hemopoietic cell kinase [EC 
hsa:3635             INPP5D, SHIP, SHIP-1, SHIP1, SIP-145, hp51CN, p150Ship; inositol polyphosphate-5-phosphatase, 145kDa 
hsa:3636             INPPL1, OPSMD, SHIP2; inositol polyphosphate phosphatase-like 1 (EC:3.1.3.86); K15909 phosphatidylin 
hsa:382              ARF6; ADP-ribosylation factor 6; K07941 ADP-ribosylation factor 6 
hsa:3984             LIMK1, LIMK, LIMK-1; LIM domain kinase 1 (EC:2.7.11.1); K05743 LIM domain kinase 1 [EC:2.7.11.1] 
hsa:3985             LIMK2; LIM domain kinase 2 (EC:2.7.11.1); K05744 LIM domain kinase 2 [EC:2.7.11.1] 
hsa:4067             LYN, JTK8, p53Lyn, p56Lyn; v-yes-1 Yamaguchi sarcoma viral related oncogene homolog (EC:2.7.10.2); K 
hsa:4082             MARCKS, 80K-L, MACS, PKCSL, PRKCSL; myristoylated alanine-rich protein kinase C substrate; K12561 my 
hsa:4651             MYO10; myosin X; K12559 myosin X 
hsa:5058             PAK1, PAKalpha; p21 protein (Cdc42/Rac)-activated kinase 1 (EC:2.7.11.1); K04409 p21-activated kinas 
hsa:50807            ASAP1, AMAP1, CENTB4, DDEF1, PAG2, PAP, ZG14P; ArfGAP with SH3 domain, ankyrin repeat and PH domain  
hsa:5290             PIK3CA, CLOVE, CWS5, MCAP, MCM, MCMTC, PI3K, p110-alpha; phosphatidylinositol-4,5-bisphosphate 3-kin 
hsa:5291             PIK3CB, P110BETA, PI3K, PI3KBETA, PIK3C1; phosphatidylinositol-4,5-bisphosphate 3-kinase, catalytic  
hsa:5293             PIK3CD, APDS, IMD14, P110DELTA, PI3K, p110D; phosphatidylinositol-4,5-bisphosphate 3-kinase, catalyt 
hsa:5294             PIK3CG, PI3CG, PI3K, PI3Kgamma, PIK3, p110gamma, p120-PI3K; phosphatidylinositol-4,5-bisphosphate 3- 
hsa:5295             PIK3R1, AGM7, GRB1, p85, p85-ALPHA; phosphoinositide-3-kinase, regulatory subunit 1 (alpha); K02649  
hsa:5296             PIK3R2, MPPH, P85B, p85, p85-BETA; phosphoinositide-3-kinase, regulatory subunit 2 (beta); K02649 ph 
hsa:5321             PLA2G4A, PLA2G4, cPLA2-alpha; phospholipase A2, group IVA (cytosolic, calcium-dependent) (EC:3.1.1.5 
hsa:5335             PLCG1, NCKAP3, PLC-II, PLC1, PLC148, PLCgamma1; phospholipase C, gamma 1 (EC:3.1.4.11); K01116 phosp 
hsa:5336             PLCG2, APLAID, FCAS3, PLC-IV, PLC-gamma-2; phospholipase C, gamma 2 (phosphatidylinositol-specific)  
hsa:5337             PLD1; phospholipase D1, phosphatidylcholine-specific (EC:3.1.4.4); K01115 phospholipase D1/2 [EC:3.1 
hsa:5338             PLD2; phospholipase D2 (EC:3.1.4.4); K01115 phospholipase D1/2 [EC:3.1.4.4] 
hsa:55616            ASAP3, ACAP4, CENTB6, DDEFL1, UPLC1; ArfGAP with SH3 domain, ankyrin repeat and PH domain 3; K12488  
hsa:5578             PRKCA, AAG6, PKC-alpha, PKCA, PRKACA; protein kinase C, alpha (EC:2.7.11.13); K02677 classical prote 
hsa:5579             PRKCB, PKC-beta, PKCB, PRKCB1, PRKCB2; protein kinase C, beta (EC:2.7.11.13); K02677 classical prote 
hsa:5580             PRKCD, CVID9, MAY1, PKCD, nPKC-delta; protein kinase C, delta (EC:2.7.10.2 2.7.11.13); K06068 novel  
hsa:5581             PRKCE, PKCE, nPKC-epsilon; protein kinase C, epsilon (EC:2.7.11.13); K18050 novel protein kinase C e 
hsa:5582             PRKCG, PKC-gamma, PKCC, PKCG, SCA14; protein kinase C, gamma (EC:2.7.11.13); K02677 classical protei 
hsa:5594             MAPK1, ERK, ERK-2, ERK2, ERT1, MAPK2, P42MAPK, PRKM1, PRKM2, p38, p40, p41, p41mapk, p42-MAPK; mitog 
hsa:5595             MAPK3, ERK-1, ERK1, ERT2, HS44KDAP, HUMKER1A, P44ERK1, P44MAPK, PRKM3, p44-ERK1, p44-MAPK; mitogen-a 
hsa:5604             MAP2K1, CFC3, MAPKK1, MEK1, MKK1, PRKMK1; mitogen-activated protein kinase kinase 1 (EC:2.7.12.2); K 
hsa:56848            SPHK2, SK_2, SK-2, SPK_2, SPK-2; sphingosine kinase 2 (EC:2.7.1.91); K04718 sphingosine kinase [EC:2 
hsa:5788             PTPRC, B220, CD45, CD45R, GP180, L-CA, LCA, LY5, T200; protein tyrosine phosphatase, receptor type,  
hsa:5879             RAC1, Rac-1, TC-25, p21-Rac1; ras-related C3 botulinum toxin substrate 1 (rho family, small GTP bind 
hsa:5880             RAC2, EN-7, Gx, HSPC022, p21-Rac2; ras-related C3 botulinum toxin substrate 2 (rho family, small GTP 
hsa:5894             RAF1, CRAF, NS5, Raf-1, c-Raf; v-raf-1 murine leukemia viral oncogene homolog 1 (EC:2.7.11.1); K0436 
hsa:6198             RPS6KB1, PS6K, S6K, S6K-beta-1, S6K1, STK14A, p70_S6KA, p70(S6K)-alpha, p70-S6K, p70-alpha; ribosoma 
hsa:6199             RPS6KB2, KLS, P70-beta, P70-beta-1, P70-beta-2, S6K-beta2, S6K2, SRK, STK14B, p70(S6K)-beta, p70S6Kb 
hsa:65108            MARCKSL1, F52, MACMARCKS, MLP, MLP1, MRP; MARCKS-like 1; K13536 MARCKS-related protein 
hsa:653361           NCF1, NCF1A, NOXO2, SH3PXD1A, p47phox; neutrophil cytosolic factor 1; K08011 neutrophil cytosolic fa 
hsa:6850             SYK, p72-Syk; spleen tyrosine kinase (EC:2.7.10.2); K05855 spleen tyrosine kinase [EC:2.7.10.2] 
hsa:7408             VASP; vasodilator-stimulated phosphoprotein; K06274 vasodilator-stimulated phosphoprotein 
hsa:7409             VAV1, VAV; vav 1 guanine nucleotide exchange factor; K05730 guanine nucleotide exchange factor VAV 
hsa:7410             VAV2, VAV-2; vav 2 guanine nucleotide exchange factor; K05730 guanine nucleotide exchange factor VAV 
hsa:7454             WAS, IMD2, SCNX, THC, THC1, WASP; Wiskott-Aldrich syndrome; K05747 Wiskott-Aldrich syndrome protein 
hsa:81873            ARPC5L, ARC16-2; actin related protein 2/3 complex, subunit 5-like; K05754 actin related protein 2/3 
hsa:8394             PIP5K1A; phosphatidylinositol-4-phosphate 5-kinase, type I, alpha (EC:2.7.1.68); K00889 1-phosphatid 
hsa:8395             PIP5K1B, MSS4, STM7; phosphatidylinositol-4-phosphate 5-kinase, type I, beta (EC:2.7.1.68); K00889 1 
hsa:8398             PLA2G6, CaI-PLA2, GVI, INAD1, IPLA2-VIA, NBIA2, NBIA2A, NBIA2B, PARK14, PLA2, PNPLA9, iPLA2, iPLA2be 
hsa:8503             PIK3R3, p55, p55-GAMMA; phosphoinositide-3-kinase, regulatory subunit 3 (gamma); K02649 phosphoinosi 
hsa:85477            SCIN; scinderin; K05768 gelsolin 
hsa:8611             PPAP2A, LLP1a, LPP1, PAP-2a, PAP2; phosphatidic acid phosphatase type 2A (EC:3.1.3.4); K01080 phosph 
hsa:8612             PPAP2C, LPP2, PAP-2c, PAP2-g; phosphatidic acid phosphatase type 2C (EC:3.1.3.4); K01080 phosphatida 
hsa:8613             PPAP2B, Dri42, LPP3, PAP2B, VCIP; phosphatidic acid phosphatase type 2B (EC:3.1.3.4); K01080 phospha 
hsa:8853             ASAP2, AMAP2, CENTB3, DDEF2, PAG3, PAP, Pap-alpha, SHAG1; ArfGAP with SH3 domain, ankyrin repeat and 
hsa:8877             SPHK1, SPHK; sphingosine kinase 1 (EC:2.7.1.91); K04718 sphingosine kinase [EC:2.7.1.91] 
hsa:8936             WASF1, SCAR1, WAVE, WAVE1; WAS protein family, member 1; K05753 WAS protein family, member 1 
hsa:8976             WASL, N-WASP, NWASP; Wiskott-Aldrich syndrome-like; K05747 Wiskott-Aldrich syndrome protein 
hsa:9846             GAB2; GRB2-associated binding protein 2; K08091 growth factor receptor bound protein 2-associated pr 
hsa:998              CDC42, CDC42Hs, G25K; cell division cycle 42; K04393 cell division control protein 42 


http://pid.nci.nih.gov/MoleculePage?molid=503236


 http://pid.nci.nih.gov/search/advanced_landing.shtml?what=graphic&svg=&jpg=true&xml=&biopax=&complex_uses=on&family_uses=on&degree=1&molecule=&pathway=FANCONI&macro_process=&source_id=5&evidence_code=NIL&evidence_code=IAE&evidence_code=IC&evidence_code=IDA&evidence_code=IFC&evidence_code=IGI&evidence_code=IMP&evidence_code=IOS&evidence_code=IPI&evidence_code=RCA&evidence_code=RGE&evidence_code=TAS&output-format=graphic&Submit=Go

 colnames(clusters)[1:13]

                                                                                                    for
combine<-rbind(clusters[,1],)

                                                                                                    
 Extra fancoi tsoff
FANCD2	FANCD2/FANCI/BRCA2/PALB2,FANCG/BRCA2/FANCD2/XRCC3,FANCD2/FANCI,FANCD2/FANCI/FAN1,FANCD2/FANCI/H2AX,FANCD2/FANCI,FANCD2/FANCI
FANCJ	FANCJ/BLM/TOP3A/BLAP75
FANCM	FA core complex/FANCM/FAAP24/MHF1/MHF2/BLM/TOP3A/BLAP75,FANCM/FAAP24/MHF1/MHF2
H2AX	FANCD2/FANCI/H2AX
PALB2	FANCD2/FANCI/BRCA2/PALB2
FANCG	FA core complex/FANCM/FAAP24/MHF1/MHF2/BLM/TOP3A/BLAP75,FANCA/FANCB/FANCC/FANCF/FANCG/FANCL/UBE2T/HES1/FAAP100,FANCG/BRCA2/FANCD2/XRCC3,FA core complex,FA core complex
FAAP24	FA core complex/FANCM/FAAP24/MHF1/MHF2/BLM/TOP3A/BLAP75,FANCM/FAAP24/MHF1/MHF2
beta TrCP1-2	
BRCA2	FANCD2/FANCI/BRCA2/PALB2,FANCG/BRCA2/FANCD2/XRCC3
FAN1	FANCD2/FANCI/FAN1
BRCA1	
XRCC3	FANCG/BRCA2/FANCD2/XRCC3
FANCE	FA core complex/FANCM/FAAP24/MHF1/MHF2/BLM/TOP3A/BLAP75,FA core complex,FA core complex
CHK1	



                                                                                                    
BLM
BLAP75
TOP3A
                                                                                                    
BRCA2
RMI1
BRCA1
FAH
ATM
TOP3A
MRE11A
BLM
XRCC3
FAN1

Molecules	Uses in complexes (duplicate names indicate differences in component properties)
FANCD2	FANCD2/FANCI/BRCA2/PALB2,FANCG/BRCA2/FANCD2/XRCC3,FANCD2/FANCI,FANCD2/FANCI/FAN1,FANCD2/FANCI/H2AX,FANCD2/FANCI,FANCD2/FANCI
FANCJ	FANCJ/BLM/TOP3A/BLAP75
FANCM	FA core complex/FANCM/FAAP24/MHF1/MHF2/BLM/TOP3A/BLAP75,FANCM/FAAP24/MHF1/MHF2
H2AX	FANCD2/FANCI/H2AX
PALB2	FANCD2/FANCI/BRCA2/PALB2
FANCG	FA core complex/FANCM/FAAP24/MHF1/MHF2/BLM/TOP3A/BLAP75,FANCA/FANCB/FANCC/FANCF/FANCG/FANCL/UBE2T/HES1/FAAP100,FANCG/BRCA2/FANCD2/XRCC3,FA core complex,FA core complex
FAAP24	FA core complex/FANCM/FAAP24/MHF1/MHF2/BLM/TOP3A/BLAP75,FANCM/FAAP24/MHF1/MHF2
beta TrCP1-2	
BRCA2	FANCD2/FANCI/BRCA2/PALB2,FANCG/BRCA2/FANCD2/XRCC3
FAN1	FANCD2/FANCI/FAN1
BRCA1	
XRCC3	FANCG/BRCA2/FANCD2/XRCC3
FANCE	FA core complex/FANCM/FAAP24/MHF1/MHF2/BLM/TOP3A/BLAP75,FA core complex,FA core complex
CHK1	

                                                                                                    

ubinquitin mediated proteolysis

K02207               UBE2R, UBC3, CDC34; ubiquitin-conjugating enzyme E2 R [EC:6.3.2.19] 
K03094               SKP1, CBF3D; S-phase kinase-associated protein 1 
K03175               TRAF6; TNF receptor-associated factor 6 
K03178               UBE1, UBA1; ubiquitin-activating enzyme E1 [EC:6.3.2.19] 
K03347               CUL1, CDC53; cullin 1 
K03348               APC1; anaphase-promoting complex subunit 1 
K03349               APC2; anaphase-promoting complex subunit 2 
K03350               APC3, CDC27; anaphase-promoting complex subunit 3 
K03351               APC4; anaphase-promoting complex subunit 4 
K03352               APC5; anaphase-promoting complex subunit 5 
K03353               APC6, CDC16; anaphase-promoting complex subunit 6 
K03354               APC7; anaphase-promoting complex subunit 7 
K03355               APC8, CDC23; anaphase-promoting complex subunit 8 
K03356               APC9; anaphase-promoting complex subunit 9 
K03357               APC10, DOC1; anaphase-promoting complex subunit 10 
K03358               APC11; anaphase-promoting complex subunit 11 
K03359               APC12, CDC26; anaphase-promoting complex subunit 12 
K03360               GRR1; F-box and leucine-rich repeat protein GRR1 
K03361               CDC4; F-box and WD-40 domain protein CDC4 
K03362               FBXW1_11, BTRC, beta-TRCP; F-box and WD-40 domain protein 1/11 
K03363               CDC20; cell division cycle 20, cofactor of APC complex 
K03364               CDH1; cell division cycle 20-like protein 1, cofactor of APC complex 
K03868               RBX1, ROC1; RING-box protein 1 
K03869               CUL3; cullin 3 
K03870               CUL2; cullin 2 
K03871               VHL; von Hippel-Lindau disease tumor supressor 
K03872               TCEB1; transcription elongation factor B, polypeptide 1 
K03873               TCEB2; transcription elongation factor B, polypeptide 2 
K03875               SKP2, FBXL1; F-box and leucine-rich repeat protein 1 (S-phase kinase-associated protein 2) 
K03876               MEL26; BTB domain ubiquitin protein ligase cofactor 
K04416               MAP3K1, MEKK1; mitogen-activated protein kinase kinase kinase 1 [EC:2.7.11.25] 
K04506               SIAH1; E3 ubiquitin-protein ligase SIAH1 [EC:6.3.2.19] 
K04552               UBE2L3, UBCH7; ubiquitin-conjugating enzyme E2 L3 [EC:6.3.2.19] 
K04553               UBE2L6, UBCH8; ubiquitin-conjugating enzyme E2 L6 [EC:6.3.2.19] 
K04554               UBE2J2, NCUBE2, UBC6; ubiquitin-conjugating enzyme E2 J2 [EC:6.3.2.19] 
K04555               UBE2G2, UBC7; ubiquitin-conjugating enzyme E2 G2 [EC:6.3.2.19] 
K04556               PARK2; parkin [EC:6.3.2.19] 
K04649               HIP2, UBC1; ubiquitin-conjugating enzyme (huntingtin interacting protein 2) [EC:6.3.2.19] 
K04678               SMURF; E3 ubiquitin ligase SMURF1/2 [EC:6.3.2.19] 
K04694               SOCS1, JAB; suppressor of cytokine signaling 1 
K04696               SOCS3, CIS3; suppressor of cytokine signaling 3 
K04706               PIAS1; E3 SUMO-protein ligase PIAS1 [EC:6.3.2.-] 
K04707               CBL; E3 ubiquitin-protein ligase CBL [EC:6.3.2.19] 
K04725               XIAP, BIRC4; E3 ubiquitin-protein ligase XIAP 
K05630               WWP2, AIP2; atrophin-1 interacting protein 2 (WW domain containing E3 ubiquitin protein ligase 2) [E 
K05632               AIP4, ITCH; atrophin-1 interacting protein 4 [EC:6.3.2.19] 
K05633               AIP5, WWP1; atrophin-1 interacting protein 5 (WW domain containing E3 ubiquitin protein ligase 1) [E 
K06643               MDM2; E3 ubiquitin-protein ligase Mdm2 [EC:6.3.2.19] 
K06688               UBE2C, UBC11; ubiquitin-conjugating enzyme E2 C [EC:6.3.2.19] 
K06689               UBE2D_E, UBC4, UBC5; ubiquitin-conjugating enzyme E2 D/E [EC:6.3.2.19] 
K07868               RHOBTB1_2; Rho-related BTB domain-containing protein 1/2 
K08285               TRIM18, MID1; midline 1 [EC:6.3.2.19] 
K09561               STUB1, CHIP; STIP1 homology and U-box containing protein 1 [EC:6.3.2.19] 
K10054               PML, TRIM19; probable transcription factor PML 
K10099               FBXO2, NFB42; F-box protein 2 
K10140               DDB2; DNA damage-binding protein 2 
K10143               RFWD2, COP1; E3 ubiquitin-protein ligase RFWD2 [EC:6.3.2.19] 
K10144               RCHY1, PIRH2; RING finger and CHY zinc finger domain-containing protein 1 [EC:6.3.2.19] 
K10259               MET30; F-box and WD-40 domain protein MET30 
K10260               FBXW7, SEL10; F-box and WD-40 domain protein 7 
K10264               FBXW8; F-box and WD-40 domain protein 8 
K10291               FBXO4; F-box protein 4 
K10447               KLHL9_13; kelch-like protein 9/13 
K10456               KLHL19, KEAP1, INRF2; kelch-like protein 19 
K10570               ERCC8, CKN1, CSA; DNA excision repair protein ERCC-8 
K10571               DET1; de-etiolated-1 
K10573               UBE2A, UBC2, RAD6A; ubiquitin-conjugating enzyme E2 A [EC:6.3.2.19] 
K10574               UBE2B, RAD6B; ubiquitin-conjugating enzyme E2 B [EC:6.3.2.19] 
K10575               UBE2G1, UBC7; ubiquitin-conjugating enzyme E2 G1 [EC:6.3.2.19] 
K10576               UBE2H, UBC8; ubiquitin-conjugating enzyme E2 H [EC:6.3.2.19] 
K10577               UBE2I, UBC9; ubiquitin-conjugating enzyme E2 I [EC:6.3.2.19] 
K10578               UBE2J1, NCUBE1, UBC6; ubiquitin-conjugating enzyme E2 J1 [EC:6.3.2.19] 
K10579               UBE2M, UBC12; ubiquitin-conjugating enzyme E2 M [EC:6.3.2.19] 
K10580               UBE2N, BLU, UBC13; ubiquitin-conjugating enzyme E2 N [EC:6.3.2.19] 
K10581               UBE2O; ubiquitin-conjugating enzyme E2 O [EC:6.3.2.19] 
K10582               UBE2Q; ubiquitin-conjugating enzyme E2 Q [EC:6.3.2.19] 
K10583               UBE2S, E2EPF; ubiquitin-conjugating enzyme E2 S [EC:6.3.2.19] 
K10584               UBE2U; ubiquitin-conjugating enzyme E2 U [EC:6.3.2.19] 
K10585               UBE2Z; ubiquitin-conjugating enzyme E2 Z [EC:6.3.2.19] 
K10586               BIRC6, BRUCE; baculoviral IAP repeat-containing protein 6 (apollon) [EC:6.3.2.19] 
K10587               UBE3A, E6AP; ubiquitin-protein ligase E3 A [EC:6.3.2.19] 
K10588               UBE3B; ubiquitin-protein ligase E3 B [EC:6.3.2.19] 
K10589               UBE3C; ubiquitin-protein ligase E3 C [EC:6.3.2.19] 
K10590               TRIP12; E3 ubiquitin-protein ligase TRIP12 [EC:6.3.2.19] 
K10591               NEDD4, RSP5; E3 ubiquitin-protein ligase NEDD4 [EC:6.3.2.19] 
K10592               HUWE1, MULE, ARF-BP1; E3 ubiquitin-protein ligase HUWE1 [EC:6.3.2.19] 
K10593               EDD1, UBR5; E3 ubiquitin-protein ligase EDD1 [EC:6.3.2.19] 
K10594               HERC1; E3 ubiquitin-protein ligase HERC1 [EC:6.3.2.19] 
K10595               HERC2; E3 ubiquitin-protein ligase HERC2 [EC:6.3.2.19] 
K10596               UBE4A; ubiquitin conjugation factor E4 A [EC:6.3.2.19] 
K10597               UBE4B, UFD2; ubiquitin conjugation factor E4 B [EC:6.3.2.19] 
K10598               PPIL2, CYC4, CHP60; peptidyl-prolyl cis-trans isomerase-like 2 [EC:5.2.1.8] 
K10599               PRPF19, PRP19; pre-mRNA-processing factor 19 [EC:6.3.2.19] 
K10600               UBOX5, UIP5; U-box domain-containing protein 5 
K10601               SYVN1, HRD1; E3 ubiquitin-protein ligase synoviolin [EC:6.3.2.19] 
K10602               NHLRC1; E3 ubiquitin-protein ligase NHLRC1 [EC:6.3.2.19] 
K10603               AIRE; autoimmune regulator 
K10604               MGRN1; E3 ubiquitin-protein ligase MGRN1 [EC:6.3.2.19] 
K10605               BRCA1; breast cancer type 1 susceptibility protein 
K10606               FANCL, PHF9; E3 ubiquitin-protein ligase FANCL [EC:6.3.2.19] 
K10607               TRIM32, HT2A; tripartite motif-containing protein 32 [EC:6.3.2.19] 
K10608               TRIM37, MUL; tripartite motif-containing protein 37 [EC:6.3.2.19] 
K10609               CUL4; cullin 4 
K10610               DDB1; DNA damage-binding protein 1 
K10611               RBX2, ROC2, RNF7; RING-box protein 2 
K10612               CUL5; cullin 5 
K10613               CUL7; cullin 7 
K10614               HERC3; E3 ubiquitin-protein ligase HERC3 [EC:6.3.2.19] 
K10615               HERC4; E3 ubiquitin-protein ligase HERC4 [EC:6.3.2.19] 
K10684               UBLE1A, SAE1; ubiquitin-like 1-activating enzyme E1 A [EC:6.3.2.19] 
K10685               UBLE1B, SAE2, UBA2; ubiquitin-like 1-activating enzyme E1 B [EC:6.3.2.19] 
K10686               UBE1C, UBA3; ubiquitin-activating enzyme E1 C [EC:6.3.2.19] 
K10687               UBE2F; ubiquitin-conjugating enzyme E2 F [EC:6.3.2.19] 
K10688               UBE2W, UBC16; ubiquitin-conjugating enzyme E2 W [EC:6.3.2.19] 
K10698               UBE1L, UBA7; ubiquitin-activating enzyme E1-like [EC:6.3.2.19] 
K10699               UBE1L2, UBA6; ubiquitin-activating enzyme E1-like protein 2 [EC:6.3.2.19] 
K12416               SWM1, APC13; anaphase-promoting complex subunit SWM1 
K12456               APC13; anaphase-promoting complex subunit 13 
K13305               NEDD4L; E3 ubiquitin-protein ligase NEDD4-like [EC:6.3.2.19] 
K16060               BIRC2_3; baculoviral IAP repeat-containing protein 2/3 
K16061               BIRC7_8; baculoviral IAP repeat-containing protein 7/8 
K16063               PIAS2; E3 SUMO-protein ligase PIAS2 [EC:6.3.2.-] 
K16064               PIAS3; E3 SUMO-protein ligase PIAS3 [EC:6.3.2.-] 
K16065               PIAS4; E3 SUMO-protein ligase PIAS4 [EC:6.3.2.-] 

                                                                                                    

myotubularin related protein 15 (1017 aa) 


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
write.table(indels,file="RunCovar_AML_pass_noControl_Coding_filtered.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


} # ifam loop


q()

getwd()
#####################
######################

pheno.types<-c("BMD_EFF_STD_HIP")
names(pheno.types)<-c("HIP")



## files<-dir(getwd())
## files<-files[files %in% "*burden*"]
## files
## burden.files<-grepl("^AOGC_sequence_skatO",files)

## the.extension<-paste(fam[ifam],project.extension,"$",sep="")
## project.files<-files[grepl(the.extension ,files)]
## print(sort(paste("Doing: ",project.files,sep=""))) # project.files<-project.files[1:22]




#indels<-{}
the.col<-{}
project.file
setwd(analysis.dir)

ipheno<-1
for(ipheno in 1:length(pheno.types)){
target.pheno<-names(pheno.types)[ipheno]

iregion<-1 # iregion<-1
region.labels<-c("PILOT.GENE.regions","GENE.regions")
region.labels<-c("FILTERED.PILOT.GENE.regions","FILTERED.GENE.regions")

region.labels<-c("AML.regions")
region.labels<-c("Coding.noControl.AML.regions")

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

}  ## loop over chr


order.by<-order(all.res[,"p.skatO"],decreasing=FALSE)

all.res<-all.res[order.by,]
all.res[1:45,]



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


write.table(all.res.gene,file="AML.coding.noControl.filtered.Gene.based.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(all.res.pilot,file="filtered.Gene.EXON.based.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

all.res[1:15,1:5]
all.res.pilot[1:5,]
all.res.gene[1:6,]

region.cols.wanted<-c("gene.skatO","p.skatO","p.burden","p.skat","nmiss.skatO","nsnps.skatO","nsnpsTotal.burden","errflag.skatO")

posns<-match(all.res.gene[,"Gene.Names"],all.res.pilot[,"hgnc_symbol"])

all.res.pilot.subset<-all.res.pilot[posns,region.cols.wanted]
colnames(all.res.pilot.subset)<-paste(colnames(all.res.pilot.subset),"EXONS",sep=".")


cbind(all.res.gene,all.res.pilot[posns,region.cols.wanted])[1:5,]

cbind(all.res.gene[,1:8],all.res.pilot.subset,all.res.gene[,9:dim(all.res.gene)[2]])[1:5,]

gene.region.compare<-cbind(all.res.gene[,1:8],all.res.pilot.subset,all.res.gene[,9:dim(all.res.gene)[2]])
write.table(gene.region.compare,file="filtered.Compare.Gene.EXON.based.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

lim = 28
qq <- parse.pvals.qq(gene.region.compare[,"p.skatO"],lim=lim)
qq <- parse.pvals.qq(all.res.gene[,"p.skatO"],lim=lim)
 
ylab = expression(Observed~~-log[10](italic(p)))
xlab = expression(Expected~~-log[10](italic(p)))
 
plot(qq$e,
     qq$o,
     xlim=c(0,lim),ylim=c(0,lim),
     pch=20,col='deepskyblue',
     xlab=xlab,ylab=ylab, main="Gene Based")
     abline(coef=c(0,1), col=1, lwd=2)

savePlot("FilteredGeneBased.jpg",type="jpeg")

#write.table(geno.all,file=paste(project.name,fam[ifam],"geno.all.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

#geno.all<-read.delim(paste(project.name,fam[ifam],"geno.all.txt",sep="."),header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)


parse.pvals.qq <- function(pvector,lim=7) {

  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )

  o[o>lim] <- lim

  out <- list(o=o,e=e)

  return(out)
}



## test<-c("MYCBP2","SLC25A24","TMCO3","C13orf35")
## test<-c("MYCBP2","SLC25A24","TMCO3","C13orf35")
## test<-c("chr22:41252435-41252687:ST13")
rownames(all.res.gene)<-1:dim(all.res.gene)[1]
test<-c("FLT3","NPM1","DNMT3A","IDH2","IDH1","TET2","RUNX1","TP53","NRAS","CEBPA","WT1","PTPN11","KIT","U2AF1","KRAS","SMC1A","SMC3","PHF6","STAG2","RAD21","FAM5C","EZH2","NHRNPK")

all.res.gene[all.res.gene[,"gene.skatO"] %in% test,]
meta.results.burden[meta.results.burden[,"gene"] %in% test,]

## loci<-snpinfo[snpinfo[,"gene"] %in% test,"Name"]
## #summary.geno.extra.out[loci,]
## high.missing.out[loci,]
## qual[loci,]
## a.indel[loci,extra]
## snpinfo[1:5,]
## qual[1:5,c("FILTER_PASS", "FILTER_100" )]








































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

