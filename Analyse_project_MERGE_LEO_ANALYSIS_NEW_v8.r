## the.sample.sheet

## sample.sheet.full<-read.delim(the.sample.sheet,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
## sample.sheet.full[1:5,]
## colnames(sample.sheet.full)
## dim(sample.sheet.full)

## table(sample.sheet.full[,7])
## AK<- grepl("AK1$",sample.sheet.full[,3]) | grepl("AK2$",sample.sheet.full[,3]) | grepl("AK$",sample.sheet.full[,3])
## sample.sheet.full[AK,]

## table(sample.sheet.full[,7])
## PD<- grepl("PD1$",sample.sheet.full[,3]) | grepl("PD2$",sample.sheet.full[,3]) | grepl("PD$",sample.sheet.full[,3])
## sample.sheet.full[PD,]

## table(sample.sheet.full[,7])
## SCC<- grepl("SCC1$",sample.sheet.full[,3]) | grepl("SCC2$",sample.sheet.full[,3]) | grepl("SCC$",sample.sheet.full[,3])
## sample.sheet.full[SCC,]

## Normal<-sample.sheet.full[,7] ==2 & !(AK | SCC | PD)

## sample.sheet.full[Normal,]

## sample.sheet.full[AK,"SampleProject"]<-"AK"
## sample.sheet.full[PD,"SampleProject"]<-"PD"
## sample.sheet.full[SCC,"SampleProject"]<-"SCC"
## sample.sheet.full[Normal,"SampleProject"]<-"Normal"

## table(sample.sheet.full[,"SampleProject"])
## setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-01-15_LeoPharma_Dec2014Freeze/BAM")
## write.table(sample.sheet.full,file="LeoPharma_Dec2014.Sample_Sheet_NEW.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## ######## ONLY NEED TO CHOOSE A DIRECTORY AND EXTENSIONS - used tab delimited files 
###############################################
###### THIS IS the super annotion run USE ALL THE FILTER AND NOVEL DATABASES
##Build AOGC

##################### if have a genotype component
## genotype.file.location<-"/media/UQCCG-Analysis/AOGC_exome_chip/working_genotypes"
## genotype.file.prefix<-"recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL"


## related.to.remove<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.recode.txt",header=F,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
## #related.to.remove<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.txt",header=F,fill=TRUE,sep=",",stringsAsFactors=FALSE)
## related.to.remove[1:5,]
## dim(related.to.remove)
#sum(related.to.remove[,1] %in% fam[,1])

qc<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Vcf_Merge.GQ-20.Min.ALT-41_2015-02-11_LeoPharmFinalAnalysis.ARE_RELATED.Thu_Mar_12_2015.genetic_QC.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

qc[1:5,]
related.controls<-qc[,"AffectionStatus"]==1 & qc[,"IBS"] < 0.9 &  qc[,"IBS"] > 0.2
sum(related.controls)
qc[related.controls,][1:20,]
bad.controls<-unique(qc[related.controls,"sample_B"])

extra<-unique(qc[grepl("^SKDP-200",qc[,"sample_B"]),"sample_B"])

extra.cases<-c("COLL-FAM-1.201","COLL-FAM-1.3","LPH-001-10_AK1","LPH-001-14_SCC","LPH-001-13_AK1","LPH-001-25_AK2")
extra.cases %in% qc[,"sample_B"]

bad.controls<-unique(c(bad.controls,extra,extra.cases))
length(bad.controls)
 write.table(bad.controls,file="LeoPharma_contaimated_march_31.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
#########################################################

#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-01-15_LeoPharma_Dec2014Freeze/Analysis/2015-01-15_LeoPharma_Dec2014Freeze.chrALL.ACC.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-01-15_LeoPharma_Dec2014Freeze/Analysis/2015-01-15_LeoPharma_Dec2014Freeze.chr1.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
annotate.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-01-15_LeoPharma_Dec2014Freeze/Annotate" # dir(annotate.dir)
analysis.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-01-15_LeoPharma_Dec2014Freeze/Analysis" # dir(analysis.dir)
project.extension<-"analysis-maf-filtered.txt" ## justX1 the exterion not fam.extension!
project.name<-"2015-01-15_LeoPharma_Dec2014Freeze"
#fam<-c(".ALL.ALL_GENOTYPES_") ## FOR all chromosomes
#fam<-c(".chrALL.ACC.ALL.ALL_GENOTYPES_") ## for merged chromosomes
#fam<-c(".chrALL.ACC_good_qual.ALL.ALL_GENOTYPES_") # PASS and wanted 
fam<-c(".chrALL.ACC_wanted.ALL.ALL_GENOTYPES_")   # wanted ALL or  c() ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_SAMPLES_PHENOTYPES_Nov.1.2013_RESIDUALS.txt"
run.per.chromsome<- FALSE ## set true if doing per chromosome
the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-01-15_LeoPharma_Dec2014Freeze/BAM/LeoPharma_Dec2014.Sample_Sheet_NEW.txt"
## ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
contaminated.file<-c()
contaminated<-c()

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)
remove.cols<-c()

core.ann<-c("chr","start","end","REF","ALT","TYPE") # out put to annanlsys programs and need foe colun labels
dont.build.summary<-TRUE ##
GATK.SB<-TRUE
maf.threshold.filter.to.use<-c(0.05)

a.label<-"coding_0.001"
dont.build.summary<-TRUE

#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-01-15_LeoPharma_Dec2014Freeze/Analysis/2015-01-15_LeoPharma_Dec2014Freeze.chrALL.ACC.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/2015-02-11_LeoPharmFinalAnalysis.chrALL.ACC_wanted.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
annotate.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Annotate" # dir(annotate.dir)
analysis.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis" # dir(analysis.dir)
project.extension<-"analysis-maf-filtered.txt" ## justX1 the exterion not fam.extension!
project.name<-"2015-02-11_LeoPharmFinalAnalysis"
#fam<-c(".ALL.ALL_GENOTYPES_") ## FOR all chromosomes
#fam<-c(".chrALL.ACC.ALL.ALL_GENOTYPES_") ## for merged chromosomes
#fam<-c(".chrALL.ACC_good_qual.ALL.ALL_GENOTYPES_") # PASS and wanted 
fam<-c(".chrALL.ACC_wanted.ALL.ALL_GENOTYPES_")   # wanted ALL or  c() ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_SAMPLES_PHENOTYPES_Nov.1.2013_RESIDUALS.txt"
run.per.chromsome<- TRUE ## set true if doing per chromosome
the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/BAM/LeoPharma_Feb2015.Sample_Sheet_NEW.txt"
## ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
contaminated.file<-c("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/LeoPharma_contaimated_march_31.txt") # COLL-FAM-1.201,COLL-FAM-1.3
#contaminated<-c()
contaminated<-read.table(contaminated.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)
remove.cols<-c()

core.ann<-c("chr","start","end","REF","ALT","TYPE") # out put to annanlsys programs and need foe colun labels
dont.build.summary<-TRUE ##
GATK.SB<-TRUE
maf.threshold.filter.to.use<-c(0.025)

a.label<-"coding_0.001"
dont.build.summary<-TRUE

###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-03_LeoPharma_NovoAlign/Analysis/2015-06-03_LeoPharma_NovoAlign.chrALL.ACC_0.01.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-03_LeoPharma_NovoAlign/Analysis/2015-06-03_LeoPharma_NovoAlign.BEST.chrALL.ACC_0.025.ALL.ALL_GENOTYPES_analysis.txt

annotate.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-03_LeoPharma_NovoAlign/Annotate" # dir(annotate.dir)
analysis.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-03_LeoPharma_NovoAlign/Analysis" # dir(analysis.dir)
#project.extension<-"analysis-maf-filtered.txt" ## justX1 the exterion not fam.extension!
project.extension<-"analysis.txt" ## 
project.name<-"2015-02-11_LeoPharmFinalAnalysis"
#fam<-c(".ALL.ALL_GENOTYPES_") ## FOR all chromosomes
#fam<-c(".chrALL.ACC.ALL.ALL_GENOTYPES_") ## for merged chromosomes
#fam<-c(".chrALL.ACC_good_qual.ALL.ALL_GENOTYPES_") # PASS and wanted 
fam<-c(".BEST.chrALL.ACC_0.025.ALL.ALL_GENOTYPES_")   # wanted ALL or  c() ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_SAMPLES_PHENOTYPES_Nov.1.2013_RESIDUALS.txt"
run.per.chromsome<- FALSE ## set true if doing per chromosome
the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-03_LeoPharma_NovoAlign/BAM/LeoPharma_June2015.Sample_Sheet_NEW.csv"
## ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
contaminated.file<-""  # COLL-FAM-1.201,COLL-FAM-1.3
contaminated<-c()
#contaminated<-read.table(contaminated.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)
remove.cols<-c()

core.ann<-c("chr","start","end","REF","ALT","TYPE") # out put to annanlsys programs and need foe colun labels
dont.build.summary<-TRUE ##
GATK.SB<-TRUE
maf.threshold.filter.to.use<-c(0.01)

a.label<-"coding_0.001"
dont.build.summary<-TRUE








################################## Other input files needed - path required in file names 
gene.symbol.file.for.clusters<-"/media/UQCCG/Software/annovar/humandb/Gene_symbol_aliases.txt"
cluster.definition.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-01-15_LeoPharma_Dec2014Freeze/Analysis/SCC_clusters.csv" ## tab delimited with header
ensembleID.to.HGNC.map.file<-"/media/UQCCG/Sequencing/Data/Genomes/hg19/ENSG_to_HGNC.txt" # tab delimited with header
coverage.file<-"/media/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_SAMPLE_coverage_LEO_2015.csv" # tab delimited with header
###########################
#/media/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_SAMPLE_Fri_Feb_06_2015.txt  latest
core.ann<-c("chr","start","end","REF","ALT","TYPE") # out put to annanlsys programs and need foe colun labels
dont.build.summary<-TRUE ##
GATK.SB<-TRUE


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





##################################################### GEFOS - GENE LIST #####################################################
mafs<-colnames(a.indel)[grepl("maf",colnames(a.indel))]
############################################# POPULATION MAF FILTER - PART A
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER

maf.threshold<-0.0  #  MAF threshold for annovar calling zero useful to get back all results !!do not modify!!
maf.threshold.filter.to.use<-c(0.001,0.005,0.01)
maf.threshold.filter.to.use<-sort(as.numeric(maf.threshold.filter.to.use))

filter.cols.novel.use<-c("NHBLI_6500_ANNOVAR_ALL","NHBLI_6500_ALL","NHLBI_5400_ALL","NHLBI_5400_EUR","NHLBI_5400_AFR","1000genome","1000genome_asian","1000genome_mine","snp141","snp141_clinical","snp137","CG69","EUR_ASN_AFR_INDEL","AOGC-NGS_ALL","AOGC-NGS_ALL_OLD","Chinese") ##
filter.cols.maf.use<-c("PopFreqMax","NHBLI_6500_ANNOVAR_ALL","NHBLI_6500_ALL","NHBLI_6500_EA","NHBLI_6500_AA","NHLBI_5400_ALL","1000genome","snp141","snp137","snp135","AOGC-NGS_ALL")
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

possible.mutations<-c("frameshift substitution","nonframeshift substitution","downstream","frameshift deletion","frameshift insertion","intergenic","intronic","ncRNA_exonic","ncRNA_intronic","ncRNA_splicing","ncRNA_UTR3","ncRNA_UTR5","ncRNA_UTR5;ncRNA_UTR3","nonframeshift deletion","nonframeshift insertion","nonsynonymous SNV","splicing","stopgain SNV","stoploss SNV","synonymous SNV","unknown","upstream","upstream;downstream","UTR3","UTR5","UTR5;UTR3","stopgain","stoploss")

interesting.coding.mutations<-c("frameshift substitution","nonframeshift substitution","nonframeshift deletion","nonframeshift insertion","frameshift deletion","frameshift insertion","nonsynonymous SNV","stopgain SNV","stoploss SNV","splicing","stopgain","stoploss")

interesting.mutations.use<-c("frameshift substitution","nonframeshift substitution","nonframeshift deletion","nonframeshift insertion","frameshift deletion","frameshift insertion","nonsynonymous SNV","stopgain SNV","stoploss SNV","splicing","ncRNA_exonic","stopgain","stoploss")

wanted.noncoding.subtypes<-c("miRNA","lincRNA") # filter by interesting to prefiler and vep.noncoding so dones get ncRNA intronic ::use gerp.score.threshold.low only these subtypes

interesting.to.prefilter<-c("UTR3","UTR5","UTR5;UTR3","snoRNA","snRNA","antisense","sense_intronic","ncRNA_exonic","ncRNA_splicing") #use gerp.score.threshold

extra.vep.annotations<-c("Uploaded_variation","Gene","Feature","Protein_position","Amino_acids")

vep.types<-c("not_assigned","stop_gained","stop_lost","missense_variant","splice_acceptor_variant","splice_donor_variant","splice_region_variant","initiator_codon_variant","stop_retained_variant","incomplete_terminal_codon_variant","frameshift_variant","inframe_deletion","inframe_insertion","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_exon_variant","NC_stop_gained","NC_stop_lost","NC_splice_acceptor_variant","NC_splice_donor_variant","NC_splice_region_variant","NC_initiator_codon_variant","NC_stop_retained_variant","NC_non_coding_exon_variant","NC_incomplete_terminal_codon_variant","NC_3_prime_UTR_variant","mature_miRNA_variant","NC_5_prime_UTR_variant","TF_binding_site_variant","TFBS_ablation","TFBS_amplification","regulatory_region_variant","intron_variant","NC_intron_variant","synonymous_variant","coding_sequence_variant","NC_synonymous_variant","upstream_gene_variant","downstream_gene_variant","intergenic_variant","NC_intergenic_variant","NMD_transcript_variant","nc_transcript_variant","NC_nc_transcript_variant","feature_truncation","feature_elongation")

## vep.types<-c( "not_assigned","stop_gained","stop_lost","stop_lost,NMD_transcript_variant","stop_gained,splice_region_variant,NMD_transcript_variant","initiator_codon_variant,splice_region_variant","splice_region_variant,3_prime_UTR_variant","stop_gained,NMD_transcript_variant","missense_variant,splice_region_variant","missense_variant","splice_acceptor_variant","splice_acceptor_variant,nc_transcript_variant","splice_region_variant,3_prime_UTR_variant,NMD_transcript_variant","splice_donor_variant,nc_transcript_variant","splice_region_variant,intron_variant,NMD_transcript_variant","splice_donor_variant","splice_region_variant","splice_region_variant,5_prime_UTR_variant","splice_region_variant,synonymous_variant","splice_region_variant,intron_variant,nc_transcript_variant","splice_region_variant,non_coding_exon_variant,nc_transcript_variant","missense_variant,NMD_transcript_variant","splice_region_variant,intron_variant","NMD_transcript_variant","intron_variant,NMD_transcript_variant","mature_miRNA_variant","5_prime_UTR_variant","5_prime_UTR_variant,NMD_transcript_variant","non_coding_exon_variant,nc_transcript_variant","3_prime_UTR_variant,NMD_transcript_variant","non_coding_exon_variant","TF_binding_site_variant","intron_variant,nc_transcript_variant","synonymous_variant,NMD_transcript_variant","3_prime_UTR_variant","regulatory_region_variant","upstream_gene_variant","downstream_gene_variant","intergenic_variant","intron_variant","synonymous_variant")


vep.coding<-c("not_assigned","stop_gained","stop_lost","missense_variant","splice_acceptor_variant","splice_donor_variant","splice_region_variant","initiator_codon_variant","stop_retained_variant","incomplete_terminal_codon_variant","frameshift_variant","inframe_deletion","inframe_insertion")

vep.noncoding<-c("5_prime_UTR_variant","3_prime_UTR_variant","non_coding_exon_variant","NC_stop_gained","NC_stop_lost","NC_splice_acceptor_variant","NC_splice_donor_variant","NC_splice_region_variant","NC_initiator_codon_variant","NC_stop_retained_variant","NC_non_coding_exon_variant","NC_incomplete_terminal_codon_variant","NC_3_prime_UTR_variant","mature_miRNA_variant","NC_5_prime_UTR_variant","TF_binding_site_variant","TFBS_ablation","TFBS_amplification","regulatory_region_variant")              

vep.unwanted<-c("intron_variant","NC_intron_variant","synonymous_variant","coding_sequence_variant","NC_synonymous_variant","upstream_gene_variant","downstream_gene_variant","intergenic_variant","NC_intergenic_variant","NMD_transcript_variant","nc_transcript_variant","NC_nc_transcript_variant","feature_truncation","feature_elongation")



missense.variant<-c("nonsynonymous SNV","missense_variant")
synonymous.variant<-c("synonymous SNV","synonymous_variant")

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

#"/media/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_SAMPLE_Wed_Feb_04_2015.txt"  # leo coverage
seq.type.file<-coverage.file  ## Troels-covergae prior to top up
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

## ########check

the.sample.sheet

sample.sheet.full<-read.delim(the.sample.sheet,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
sample.sheet.full[1:5,]
colnames(sample.sheet.full)
dim(sample.sheet.full)
table(sample.sheet.full[,"SampleProject"])

##### fix 0 and 9 for missing to NA

## pheno.types<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_RAD","BMD_EFF_STD_LS","BMD_EFF_STD_FN","EVER_FX_50_EXCL_TRIVIAL")
## names(pheno.types)<-c("HIP","RAD","LS","FN","FX")
     ## AK Control  Normal      PD     SCC 
     ## 50     291      24      25      24



## extra<-all.possible.samples[!(all.possible.samples %in% sample.sheet.full[,"SAMPLE"])]
## sample.sheet.full[,"SAMPLE"][!( sample.sheet.full[,"SAMPLE"] %in% all.possible.samples)]
## all.possible.samples[grepl("\\d$",all.possible.samples,perl=TRUE) & grepl("^LPH",all.possible.samples,perl=TRUE) & !grepl("AK1$",all.possible.samples,perl=TRUE) & !grepl("AK2$",all.possible.samples,perl=TRUE) ]

## all.possible.samples[grep("SCC$",all.possible.samples,perl=TRUE)]

## all.possible.samples[grep("AK1",all.possible.samples,perl=TRUE)]

## table(sample.sheet.full$SampleProject)

pheno.types<-c("SampleProject") ## vales is column header
names(pheno.types)<-c("SampleProject") ### name is output columns

## case.control.classes<-c(0,1,0)
## names(case.control.classes)<-c("Control","NMD","AOGC")
## case.control.classes

     ## AK Control  Normal      PD     SCC 
     ## 50     292      25      25      24


  ## AK Control  Normal      PD     SCC 
  ##    50     253      24      25      24 

case.control<-c("SampleProject")
case.control.classes<-c(1,0,0,0,1)
names(case.control.classes)<-c("AK","Control","Normal","PD","SCC" )
case.control.classes
# ib<-1
for(ib in 1:length(case.control)){
  if(!(case.control[ib] %in% colnames(sample.sheet.full))){next}
  sample.sheet.full[(  !(sample.sheet.full[,case.control[ib]] %in% names(case.control.classes))  |  is.na(sample.sheet.full[,case.control[ib]]) | sample.sheet.full[,case.control[ib]]==0 | sample.sheet.full[,case.control[ib]]==9)  ,case.control[ib]]<-NA
}

table(sample.sheet.full$SampleProject)
#tapply(sample.sheet.full[,"SampleProject"],sample.sheet.full[,"SampleProject"],length)
colnames(sample.sheet.full)[colnames(sample.sheet.full)=="ParticipantCode"]<-"SAMPLE"
control.samples<-{}
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


if(run.per.chromsome){
 project.files<-project.files[!grepl("chrALL",project.files)]

  if(length(project.files)!=24){
    print("########################################### WARNING #################################")
    print("less that 24 chromosomes detected")
    print(fam[ifam])
    print("########################################### WARNING #################################") 
  }
}
  project.files


ichr<-1

for(ichr in 1:length(project.files)){

setwd(analysis.dir)
################## fast read ###########  load("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/LEO_Pharma_wanted_with_synon_ALL_0.01_muts_Feb23")
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
if(length(contaminated)>0){
bad.samples<-contaminated[,1]
bad.samples.labels<-expand.labels.to.samples(bad.samples,c("GT","AD","DP","GQ"),paste.after=TRUE)

if(length(bad.samples.labels)>1){
a.indel<-a.indel[,colnames(a.indel)[!(colnames(a.indel) %in% bad.samples.labels)]]
}
}
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

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
sum(missing)
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

ens.to.hgnc<-read.delim(ensembleID.to.HGNC.map.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
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


##               MUC4                TTN              MUC16              CTBP2            FAM182B               MUC6              MUC5B              NBPF1               RYR2              CDC27              MUC12              MST1L             ZNF717 
##               1010                528                355                308                232                228                217                214                196                195                195                171                170 
##            ANKRD36              DNAH5            FAM182A              HYDIN              NBPF8           HLA-DRB1              CSMD3                NEB              LRP1B              SYNE1            NBPF25P              CSMD1               MUC2 
##                160                156                154                154                149                148                147                147                144                144                142                140                138 
## PTGER4P2-CDK2AP2P2              USH2A             MST1P2              DLEC1 
##                137                136                133                132

grep("NOTCH1",names(all.genes))
## common.hit.genes<-names(all.genes)[1:30]
## common.hit.genes<-all.genes[all.genes>400]
common.hit.genes<-{}
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


clusters<-read.delim(cluster.definition.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
clusters

#clusters.wanted<-c("Clinical","FANC - ACID")
clusters.wanted<-colnames(clusters)
ic<-1
snpinfo[1:5,]

#cbind(unique(clusters[,22]),unique(clusters[,22]))
snpinfo[1:5,]



gene.aliases<-read.delim(gene.symbol.file.for.clusters,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
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
#save(list=c("a.indel"),file="indels.RData")

num.cores<-2
library(doMC)
registerDoMC(cores=num.cores)
## library("doParallel")
## registerDoParallel(cores=num.cores)
num.bits<-num.cores
#



the.samples<-colnames(a.indel)[grepl(".GT$", colnames(a.indel))]

while((dim(a.indel)[1] %% num.bits)< 2){num.bits<-num.bits+1} ### go don't get matrix issues
num.bits
#(dim(a.indel)[1] %% num.bits)
## fil.genotypes<-foreach(a.indel.bit=iter(a.indel,by='row',chunksize=as.integer(dim(a.indel)[1]/num.bits) ), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype(a.indel.bit,gsub(".GT$","",the.samples),prefix="",suffix="",20,0.02,0.98,0.20,0.80,7,2)

fil.genotypes<-foreach(a.indel.bit=iter(a.indel,by='row',chunksize=as.integer(dim(a.indel)[1]/num.bits) ), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype(a.indel.bit,gsub(".GT$","",the.samples),prefix="",suffix="",20,0.02,0.98,0.20,0.80,7,2)
# 20,0.02,0.98,0.20,0.8,10,5 # for cancers where het may be way out

# 0.02< het(0/1)<0.98
# ref(0/0)< 0.2
# (1/1) > 0.8
## rownames(fil.genotypes)<-key
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
# for(ipheno in 1:length(pheno.types)){  ## short circut look for now


##################################################### set formula
print(paste("Doing phenotype:",pheno.types[ipheno]))
target.pheno<-names(pheno.types)[ipheno]
target.pheno.col<-pheno.types[ipheno]

if(target.pheno.col %in% case.control){
  covars<-c("1") # c("AGE_SCAN","PCA1","PCA2","PCA3","PCA4") #AGE_SCAN,PCA1,PCA2,PCA3,PCA4 #covars<-c("1")
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
all.possible.samples[(!(all.possible.samples %in% sample.sheet.full[,"SAMPLE"] ) )]
subset.samples<- (sample.sheet.full[,"SAMPLE"] %in% all.possible.samples ) & !is.na(sample.sheet.full[,target.pheno.col]) & sample.sheet.full[,target.pheno.col]!="NA"  & got.all.covars
                  
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
dim(sample.sheet.full)
length(all.possible.samples)

########## make sure sample sheet has no extras

pheno[1:5,]



the.samples<-paste(pheno[,"SAMPLE"],"GT",sep=".")  ## samples same order as in pheno
print(paste("Number samples: ",length(the.samples),sep=""))

#seq.type[1:57,]
posns<-match(pheno[,"SAMPLE"],seq.type[,"Sample"])
missing<-is.na(posns)
sum(missing)
pheno[missing,1:5]

if(!exists("pheno.ori")){
pheno.ori<-pheno
}
#pheno<-pheno.ori

### Troels will update but all are Nextera
## capture<-seq.type[posns,"Capture.Method"]
## table(capture) ### all illume here
## pheno<-cbind(pheno,capture)


################################################################################
###### PLACE to set others phenotype classes here ( invasive sub types) etc
## pheno[,pheno.types[ipheno]]

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
## all<-rep(TRUE,length(the.samples))

## dim(pheno)
## length(all)  pheno<-pheno.ori
table(pheno.ori$SampleProject)

     ## AK Control  Normal      PD     SCC  # latestes with QC
     ## 47     133      24      25      23
    ## AK Control  Normal      PD     SCC 
    ##  50     133      25      25      24

cancer<-rep(FALSE,times=dim(pheno)[1])
cancer[pheno$SampleProject %in% c("AK","SCC")]<-TRUE

normal<-rep(FALSE,times=dim(pheno)[1])
normal[pheno$SampleProject %in% c("Normal")]<-TRUE

Control<-rep(FALSE,times=dim(pheno)[1])
Control[pheno$SampleProject %in% c("Control")]<-TRUE

PD<-rep(FALSE,times=dim(pheno)[1])
PD[pheno$SampleProject %in% c("PD")]<-TRUE

SCC<-rep(FALSE,times=dim(pheno)[1])
SCC[grepl("_SCC$",pheno$SAMPLE)]<-TRUE

AK<-rep(FALSE,times=dim(pheno)[1])
AK[grepl("_AK1$",pheno$SAMPLE) | grepl("_AK2$",pheno$SAMPLE)]<-TRUE


the.SCC<-pheno$SAMPLE[grepl("_SCC$",pheno$SAMPLE)]
the.AK<-pheno$SAMPLE[grepl("_AK1$",pheno$SAMPLE) | grepl("_AK2$",pheno$SAMPLE)]

sum(PD)
sum(normal)
####"LPH-001-27 Blood ok
## "LPH-001-27_PD"
## normal[pheno$SAMPLE %in% c("LPH-001-27_PD")]<-TRUE
## PD[pheno$SAMPLE %in% c("LPH-001-27_PD")]<-FALSE


pheno<-cbind(pheno,cancer,Control,normal,PD,SCC,AK)

     ## AK Control  Normal      PD     SCC 
     ## 47     133      24      25      23 
#pheno<-cbind(pheno,SCC,AK)
############################################################## pheno<-pheno.ori the.samples<-paste(pheno[,"SAMPLE"],"GT",sep=".") 



## 27 now ok
## normal[pheno$SAMPLE %in% c("LPH-001-27_PD")]<-TRUE
## PD[pheno$SAMPLE %in% c("LPH-001-27_PD")]<-FALSE
## pheno[pheno[,"SAMPLE"] %in%  c("LPH-001-27_AK1","LPH-001-27_AK2","LPH-001-27_PD","LPH-001-27_SCC"),"cancer"]<-FALSE
## pheno[pheno[,"SAMPLE"] %in%  c("LPH-001-27_AK1","LPH-001-27_AK2","LPH-001-27_PD","LPH-001-27_SCC"),"normal"]<-FALSE
## pheno[pheno[,"SAMPLE"] %in%  c("LPH-001-27_AK1","LPH-001-27_AK2","LPH-001-27_PD","LPH-001-27_SCC"),"PD"]<-FALSE


      
############## SET up groups to be analysed
table(pheno$SampleProject)

pheno[1:5,]

## th project becomes the master copy of the targets below.
the.projects<-c("cancer","Control","normal","PD","SCC","AK")
names(the.projects)<-the.projects
colnames(pheno)

pheno[pheno[,the.projects[1]],"SAMPLE"]
pheno[pheno[,the.projects[2]],"SAMPLE"]
pheno[pheno[,the.projects[3]],"SAMPLE"]
pheno[pheno[,the.projects[4]],"SAMPLE"] 


for (ir in 1: length(the.projects)){
  print(paste(the.projects[ir],"Num. samples:",sum(pheno[,the.projects[ir]])))
      }
## [1] "cancer Num. samples: 67"
## [1] "Control Num. samples: 133"
## [1] "normal Num. samples: 24"
## [1] "PD Num. samples: 24"

##NOVO
## [1] "cancer Num. samples: 74"
## [1] "Control Num. samples: 133"
## [1] "normal Num. samples: 25"
## [1] "PD Num. samples: 25"
## [1] "SCC Num. samples: 24"
## [1] "AK Num. samples: 50"

##  > pheno[1:5,]
##   SampleProject FamilyCode SAMPLE PaternalID MaternalID Sex AffectionStatus cancer Control normal    PD
## 1             0        ALL    860          0          0  -9               1  FALSE    TRUE  FALSE FALSE
## 2             0        ALL    861          0          0  -9               1  FALSE    TRUE  FALSE FALSE
## 3             0        ALL    862          0          0  -9               1  FALSE    TRUE  FALSE FALSE
## 4             0        ALL   2529          0          0  -9               1  FALSE    TRUE  FALSE FALSE
## 5             0        ALL   2530          0          0  -9               1  FALSE    TRUE  FALSE FALSE                      
# summary.geno.extra.ori <-summary.geno.extra
summary.geno.extra<-{}
####################################################################################
#################################################################################### REGULAR
#### MAY NEED TO ADJUST way use samples in selected based on selection below.
targets<-the.projects  #c("NMD","ex.Control","AOGC")
targets
###### the.samples and pheno in same order but the.samples has .GT extension.
it<-1

the.samples<-paste(pheno[,"SAMPLE"],"GT",sep=".")
for(it in 1:length(targets)){
#use.samples<-the.samples[pheno[,"SampleProject"]==targets[it]]
use.samples<-the.samples[pheno[,targets[it]]]
print(targets[it])
print(use.samples)
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
print(targets[it])
print(use.samples)
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

#summary.geno.extra.ori<-summary.geno.extra


## ######################################### use.sample<-the.AK

## chk<-grep("^GENO",colnames(summary.geno))
## summary.geno["chr19:50169131:50169131:C:T:snp",chk]
## summary.geno.extra.ori<-summary.geno.extra
a.snp<-"chr15:40675107:40675107:C:T:snp" # KNSTRN
a.snp<-"chr19:50169131:50169131:C:T:snp" # BCL@L12
a.snp<-"chr6:31940123:31940123:G:A:snp" #STK19

chk<-grep("^GENO",colnames(summary.geno.extra))
summary.geno.extra[a.snp,chk]

## chk<-grep("^GENO",colnames(summary.geno.extra.use))
## summary.geno.extra.use[a.snp,chk]

chk<-grep("^GENO",colnames(summary.geno.extra.ori))
summary.geno.extra.ori[a.snp,chk]


###################################### GET HARDY Weinbery equilibrium of conrols
###################################### GET HARDY Weinbery equilibrium of conrols
###################################### GET HARDY Weinbery equilibrium of conrols
###################################### GET HARDY Weinbery equilibrium of conrols


a.indel[1:5,1:10]
summary.geno.extra[1:5,]
colnames(summary.geno.extra)
rownames(summary.geno.extra)<-key
#getHWE(obs_hets, obs_hom1, obs_hom2)
hw.target<-"Control"  ## what to calculate HW with

hw.p.control<-getHWE(summary.geno.extra[,paste("GENO.",hw.target,sep="")]) ## used 16 CPUs
hw.p.control.filt<-getHWE(summary.geno.extra[,paste("GENO.",hw.target,".filt",sep="")]) ## used 16 CPUs


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
rare.in.group.test<-paste("MAF.",c("cancer","Control"),sep="")
rare.in.group.test
rare.in.group<-combine.boolean(rare.in.group.table,rare.in.group.test,"OR")
sum(!rare.in.group)
summary.geno.extra[!rare.in.group,paste("GENO.",c("cancer","Control"),sep="")][1:5,]

rare.in.group[1:10]
## rare.in.group<-combine.boolean(rare.in.group.table,c("MAF.AML","MAF.Control"),"OR")
## rare.in.group.filt<-combine.boolean(rare.in.group.table,c("MAF.AML.filt","MAF.Control.filt"),"OR")
## sum(!rare.in.group)
## sum(!rare.in.group.filt)

no.genotypes.test<-paste("MAF.",c("cancer.filt","Control.filt"),sep="")
no.genotypes.test
no.genotypes<-(summary.geno.extra[,no.genotypes.test]== 0)  | (is.na( summary.geno.extra[,no.genotypes.test])) # no genotypes in test classes for a mutataion after individaul quality filtering
no.genotypes[1:5,]
no.genotypes.filt<-combine.boolean(no.genotypes,colnames(no.genotypes),"AND")
no.genotypes.filt[1:5]

no.genotypes.test<-paste("MAF.",c("cancer","Control"),sep="")
no.genotypes.test
no.genotypes<-(summary.geno.extra[,no.genotypes.test]== 0)  | (is.na( summary.geno.extra[,no.genotypes.test])) # no genotypes in test classes for a mutataion after individaul quality filtering
no.genotypes[1:5,]
no.genotypes<-combine.boolean(no.genotypes,colnames(no.genotypes),"AND")
no.genotypes[1:5]

no.genotypes.cancer.test<-paste("MAF.",c("cancer"),sep="")
no.genotypes.cancer<-(summary.geno.extra[,no.genotypes.cancer.test]== 0)

sum(no.genotypes)

summary.geno.extra[1:5,]
missing.targets<-c("cancer.filt","Control.filt")

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

high.missing.subset<-c("cancer.filt","Control.filt")

high.total.missing<-subset(high.missing.table,select=c("cancer.filt","Control.filt"))
high.total.missing[1:5,]



high.total.missing<-high.total.missing <= missing.threshold
ok.missing<-combine.boolean(high.total.missing,c("cancer.filt","Control.filt"),"AND")
sum(ok.missing)

missing.targets
the.projects
colnames(high.missing.table)


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
#maf.lt.all[1:5,]
#maf.lt.all<-as.logical(maf.lt.all)
#as.logical(maf.lt.all[1:5,1])


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


###################### look for missense
missense.coding<-test.for.coding.type(a.indel,geneanno.DB,missense.variant)
sum(missense.coding)
   ####################### mhairi CHECK
missense.coding.vep<-test.wanted.mutation(a.indel[,"Consequence.Embl"],missense.variant,delimit.by=",")  # filter.table.pholy[,"Consequence"] %in%  vep.coding
  ## a.indel[wanted.muts.coding.vep,]
sum(missense.coding.vep)

is.missense<-missense.coding.vep | missense.coding

## qual[1:50,"PolyPhen.low"] # true is above 0.1
## a.indel[1:50,"PolyPhen.scores"]
## a.indel[1:50,c("PolyPhen.desc","SIFT.desc")]
## qual[1:5,"SIFT.high"] 
## a.indel[1:5,"SIFT.scores"
synonymous<-test.for.coding.type(a.indel,geneanno.DB,synonymous.variant)
sum(synonymous)



is.benign.missense<-is.missense & !qual[,"PolyPhen.low"] & a.indel[,"PolyPhen.desc"] !="unknown" # chr9:21971016:21971016:G:A:snp

sum(!is.benign.missense) # is.benign.missense["chr9:21971016:21971016:G:A:snp"]
dim(a.indel)
###################### look for missense
###################### look for missense
###################### look for missense


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
  the.test<-grepl(a.type[itype],a.indel[,"Gene.Biotype"])
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
########################################################################## /media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/Gerp.test.r

p<-0.001
sd.thresh<-6

p

colnames(summary.geno.extra)
######################################
n<-max(as.integer(summary.geno.extra[,"TOTAL.Alleles.Control"]))
n
alt.counts.thresh<-1
while( (alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) <= sd.thresh){alt.counts.thresh<-alt.counts.thresh+1}
alt.counts.thresh

summary.geno.extra[1:5,]
rare.in.Control<-as.numeric(summary.geno.extra[,"ALT.Alleles.Control"])< alt.counts.thresh
rare.in.Control.filt<-as.numeric(summary.geno.extra[,"ALT.Alleles.Control.filt"])< alt.counts.thresh
sum(rare.in.Control)
sum(rare.in.Control.filt)
names(rare.in.Control)<-key
names(rare.in.Control.filt)<-key

summary.geno.extra[!rare.in.Control,][1:5,1:20]
##############################



#################################################################


#length(maf.filter)
#length(rare.in.controls)

#maf.lt.all[1:5,]
#maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.001"])
#maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.01"])
maf.col<-paste("MAF.lt",p,sep=":")
maf.col
if(exists("maf.lt.all")){
maf.filter<-as.logical(maf.lt.all[,maf.col])
}else{
maf.filter<-as.logical(a.indel[,maf.col])
}
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



#########################  Choose Somatic

in.any.normal<-as.numeric(summary.geno.extra[,"ALT.Alleles.normal"]) > 0
in.any.normal.filt<-as.numeric(summary.geno.extra[,"ALT.Alleles.normal.filt"]) > 0

summary.geno.extra[!in.any.normal.filt,grepl("normal",colnames(summary.geno.extra))][1:5,]
sum(!in.any.normal)






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



## pheno.use<-pheno[pheno[,target.pheno] %in% c(0,1),]
## dim(pheno.use) pheno.ori.ori<-pheno.ori
if(!exists("pheno.ori") | !identical(colnames(pheno),colnames(pheno.ori))){
pheno.ori<-pheno
}
pheno<-pheno.ori[!pheno.ori[,"PD"],]
the.samples.use<-pheno[,"SAMPLE"]
the.samples.use<-paste(the.samples.use,".GT",sep="")
table(pheno[, "SampleProject"])

sum(pheno[,"PD"]) ## pehno used to do asssocation tests - excluded PD
sum(pheno.ori[,"PD"])
sum(pheno[,"SCC"])


############################################################################################
##########################################################################################
## look for mutations mutation is germiline controls is not signifcatly different to those called in cancer

threshold=1
        
dim(pheno.ori)
length(the.samples)
pheno.ori[1:5,]

normals<-pheno.ori[pheno.ori[,"normal"],"SAMPLE"]
#normals<-normals[normals!="LPH-001-27_PD"]
normals
normal.alt.counts<-alt.reads.reference.calls(a.indel,normals,threshold=1)





Control<-pheno.ori[pheno.ori[,"Control"],"SAMPLE"]
Control
Control.alt.counts<-alt.reads.reference.calls(a.indel,Control,threshold=1)
#Control.alt.counts.true<-alt.reads.Non.reference.calls(a.indel,Control,threshold=1)
Control.alt.counts[1:5,]

##  normal.allele.freq<-Control.alt.counts[test,"Read.Balance"]/100
## tumor.sample.major.counts<-normal.alt.counts[test,"REF.reads"]
## tumor.sample.minor.counts<-normal.alt.counts[test,"ALT.reads"]
test<-c("chr19:2939268:2939289:ACCACCCTTACCCAAGGAGGCA:-:indel","chr13:21729956:21729956:A:G:snp","chr2:198363406:198363406:C:T:snp","chr22:20760282:20760282:A:C:snp")
test<-c("chr15:40675107:40675107:C:T:snp","chr19:50169131:50169131:C:T:snp","chr17:7578263:7578263:G:A:snp","chr17:7574003:7574003:G:A:snp")


## poissonian.model.position(Control.alt.counts[test,"Read.Balance"]/100,normal.alt.counts[test,"REF.reads"], normal.alt.counts[test,"ALT.reads"],av.tumor.contamination=0)





poss.model<-poissonian.model.position(Control.alt.counts[,"Read.Balance"]/100,normal.alt.counts[,"REF.reads"], normal.alt.counts[,"ALT.reads"],1)
## poss.model[test,]
poss.model[1:5,]
dim(poss.model)
dim(a.indel)
sum(rownames(Control.alt.counts)!=rownames(normal.alt.counts))
sum(rownames(a.indel)!=rownames(poss.model))

sum(poss.model[pass,"P-value"]<1e-5,na.rm=TRUE)
poss.model[1:5,]  ### used in conjection with additional snp filtering
 # 2*(pnorm(abs(c(1:8)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)) #2*(pnorm(abs(seq(from=3,to=5,by=0.1)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
# 2*(pnorm(abs(4.4), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)) ##4 sd greater than mean

pass.possonian.control.model<- poss.model[,"P-value"]> 1e-5  | is.na( poss.model[,"P-value"])

length(pass.possonian.control.model)
#length(pass)
sum(!pass.possonian.control.model) #2172

prev.contamin<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/containamination loci.csv",header=T,sep="\t",fill=TRUE,skip=0,stringsAsFactors=FALSE)
prev.contamin[1:5,1:5]

pass.possonian.control.model[1:5]
names(pass.possonian.control.model)[1:5]
prev.contamin[1:5,"gene"]

grep("chr19:50169131:50169131:C:T:snp",names(pass.possonian.control.model))

posns<-match(prev.contamin[,"gene"],names(pass.possonian.control.model))
posns[1:20]
chk<-cbind(prev.contamin[,1:5],pass.possonian.control.model[posns])
chk[1:60,]
getwd()
write.table(chk,file="contamination_BWA_Novo_Check.txt",sep="\t",row.names=F,col.names=T)

##########################################################################
##########################################################################




##########################################################################
##########################################################################
##########################################################################
## extra filtering
files<-dir("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Snp_read_filter")
files<-files[grepl("ALL.ALL_GENOTYPES_analysis-maf-filtered.readfilter.summary.txt$",files)]
files
i<-1
setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Snp_read_filter")
for (i in 1:length(files)){
a.filt<-read.delim(files[i],header=T,sep="\t",fill=TRUE,skip=5,stringsAsFactors=FALSE)
if(i==1){
  filt<-a.filt
}else{
  filt<-rbind(filt,a.filt)
}

}
setwd(analysis.dir)

dim(a.indel)
dim(filt)
a.indel[1:5,1:10]
table(filt[,"chr"])

filt.key<-build.key(filt,core.ann,add.chr.label=TRUE)
#rownames(filt)<-filt.key
filt.key[1:5]
rownames(a.indel)[1:5]

posns<-match(rownames(a.indel),filt.key)
missing<-(is.na(posns))
sum(missing)
rownames(a.indel)[missing][1:10]



filt<-filt[posns[!missing],]
filt.key<-filt.key[posns[!missing]]

sum(missing[pass])

## filt<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Snp_read_filter/LEO_Pharma_WANTED_ALL_noGQ_filter_0.01_muts_Mar06.filter.pass.somatic.readfilter.txt",header=T,sep="\t",fill=TRUE,skip=5,stringsAsFactors=FALSE)
## filt[1:5,1:20]

## filt.key<-build.key(filt,core.ann)
## rownames(filt)<-filt.key

## chr6:32063957:32063957:C:G:snp
## grep("32063957",filt[,"start"])


## filt[173713,]

pass.filt<-strsplit(filt[,"FILTER_SUMMARY"],split=";")
pass.filt[173713]
pass.filt<-unlist(lapply(pass.filt,function(x) {sum(as.logical(x[1:3]))==0}))
pass.filt[173713]
sum(!pass.filt)
filt.key[173713]
snp.fail.filt<-filt.key[!pass.filt]

bad.qual.locations<-rownames(a.indel) %in% snp.fail.filt
sum(bad.qual.locations)

filt["chr15:40675107:40675107:C:T:snp",]
"chr6:32063957:32063957:C:G:snp" %in% snp.fail.filt

### don't sue the 4 filter picking up somatic calls in cases 
##                                   chr    start      end REF ALT TYPE   GENE         FILTER_SUMMARY                     SUMMARY_CALLED     SUMMARY_NOT_CALLED
## chr15:40675107:40675107:C:T:snp chr15 40675107 40675107   C   T  snp KNSTRN FALSE;FALSE;FALSE;TRUE 686;113;4;27;86;0;0.04;0.24;0.76;0 42339;182;5;23;159;1;8
       
## other.bad<-c("DNAH5","CSMD3","CSMD1","PCLO","INTS1","TYRO3","TTN","TESK1","MUC17","OR4A15","OR6C1","BCL2L11")

## other.bad %in% bad.genes

## bad.genes<-c(bad.genes,other.bad)

## save.image("APRIL_1st_no_containiation_snpFilters_weights.RData")

## save.image("APRIL_1st_no_containiation_snpFilters_weights_possonian_rescue.RData")
##########################################################################
##########################################################################
##########################################################################

##### for Maris idea of using coding only I will additioanlly add that the the SNP must be in the CDS region of the SNP:

library(GenomicFeatures)
library(Rsamtools)
library(BSgenome)
library(Homo.sapiens)
## transcripts(Homo.sapiens, columns=c("TXNAME","SYMBOL"))
## cds(Homo.sapiens, columns=c("TXNAME","SYMBOL"))

## source("http://bioconductor.org/biocLite.R")
## biocLite("TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts")
## biocLite("mirbase.db")

## ## library(TxDb.Hsapiens.UCSC.hg19.knownGene)
## ##  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
## library(mirbase.db)
## library("TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts")

## /media/scratch/annovar/annovar_latest


## ls(TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts)
## ## microRNAs(TxDb.Hsapiens.UCSC.hg19.knownGene)
##      ## End(Not run)
## mir<-mirbase
## x<-mirbaseCHRLOC
## class(x)
## x
## mapped_keys <- mappedkeys(x)
## # Convert to a list
## xx <- as.list(x[mapped_keys])
## if(length(xx) > 0) {
## # Get the CHRLOC for the first five entries
## xx[1:5]
## }

## names(xx[1])
## lapply(xx,length )

## the.chr<-lapply(xx,function(x){ names(x)[[1]] })


## TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscript

## getwd()
## setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis")

## write.table(a.indel[,core.ann],file="Annovar.update.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
## system("/media/scratch/annovar/annovar_latest/annotate_variation.pl -regionanno -build hg19  -dbtype wgRna  Annovar.update.txt  /media/UQCCG/Software/annovar/humandb/hg19/")
## system("/media/scratch/annovar/annovar_latest/annotate_variation.pl -regionanno -build hg19  -dbtype targetScanS  Annovar.update.txt  /media/UQCCG/Software/annovar/humandb/hg19/")
## #system("/media/scratch/annovar/annovar_latest/annotate_variation.pl -regionanno -build hg19  -dbtype gerp++elem  Annovar.update.txt  /media/UQCCG/Software/annovar/humandb/hg19/")
## system("/media/scratch/annovar/annovar_latest/annotate_variation.pl -regionanno -build hg19  -dbtype gerp++gt2  Annovar.update.txt  /media/UQCCG/Software/annovar/humandb/hg19/")
## files.ann<-c("wgRna","targetScanS")# ,"gerp++elem")

## iff<-3
## all.data<-{}
## for (iff in 1:length(files)){
## the.data<-try(read.delim(paste("Annovar.update.txt.hg19",files.ann[iff],sep="_"),header=F,skip=0,sep="\t",fill=TRUE,stringsAsFactors=FALSE,colClasses = "character"),silent=TRUE)
## the.data[1:5,]
## the.data[,2]<-gsub("Name=","",the.data[,2])

## if(iff==2){
##   targetsc<-strsplit(the.data[,2],split=";")
##   score<-unlist(lapply(targetsc,function(x) x[1]))
##   miRNA.binding<-unlist(lapply(targetsc,function(x) x[2]))
##   the.data[,2]<-miRNA.binding
##   the.data[,1]<-score
## }

## key.ann<-build.key(the.data,colnames(the.data)[((dim(the.data)[2])-5):(dim(the.data)[2])])
## key.ann[1:5]

## posns<-match(key,key.ann)
## missing<-is.na(key)
## sum(missing)

## colnames(the.data)[2]<-files.ann[iff]
## colnames(the.data)[1]<-paste("Score",files.ann[iff],sep="_")

## if(iff==1){
##   all.data<-the.data[posns,c(1,2)]
## }else{
##   all.data<-cbind(all.data,the.data[posns,c(1,2)])
## }
## }

## dim(all.data)
## all.data[1:50,]
## all.data<-cbind(key,all.data)

## ita<-2
## target.add<-c("wgRna","targetScanS")
## extra.all<-{}
## for(ita in 1:length(target.add)){
## wgRNA<-all.data[!is.na(all.data[,target.add[ita]]),]

## extra<-matrix(data=NA,nrow=dim(wgRNA)[1],ncol=dim(snpinfo.ori)[2])
## colnames(extra)<-colnames(snpinfo.ori)
## extra<-as.matrix(extra)
## extra[1:5,]
## snpinfo.ori[1:5,]
## wgRNA[1:5,]
## extra[,"cluster"]<-target.add[ita]
## extra[,"gene"]<-as.character(wgRNA[,target.add[ita]])
## extra[,"Name"]<-as.character(wgRNA[,1])

## extra[1:5,]
## dim(extra)

## if(iff==1){
##   extra.all<-extra
## }else{
##   extra.all<-rbind(extra.all,extra)
## }
## }

## ### om;y run next 3 lines once!
## ## snpinfo.ori[226937:226952,"gene"]<-extra[,"gene"]
## ## interesting<-snpinfo.ori[,"cluster"] %in% c("targetScanS","wgRna")
## ## snpinfo.ori[interesting,"cluster"]<-snpinfo.ori[interesting,"gene"]
## extra.miRNA<-extra.all
## snpinfo.ori[1:5,]
## extra.all[1:5,]
## dim(extra.all)
## tail(extra.all)
## snpinfo.ori[(226952-500):226952,]
## #
## if( sum(extra.all[,"gene"] %in% snpinfo.ori[,"gene"])==0 ){snpinfo.ori<-rbind(snpinfo.ori,extra.all)} #only do once
## sum(is.na(snpinfo.ori[,"gene"]))



## ## colnames(a.indel)[1:50]
## ## sum(as.logical(a.indel[, "wanted.muts.coding"]))

## write.table(a.indel[inspect,1:30],file="TTN_anotation_complications.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)




all.coding<-cds(Homo.sapiens, columns=c("CDSID"))
## all.coding[1:5,]
all.coding<-reduce(all.coding) ## so all locations are unque

head(start(all.coding))

a.indel.gr<-GRanges(seqnames =a.indel[,"chr"],ranges = IRanges(start=as.numeric(a.indel[,"start"]),end=as.numeric(a.indel[,"end"])),strand="*")
a.indel.gr


dim(a.indel)
sum(as.logical(a.indel[,"wanted.muts"]))
sum(as.logical(a.indel[,"wanted.muts.coding"]))


#the.coding.snps <- findOverlaps(a.indel.gr,all.coding,ignore.strand = TRUE)
the.coding.snps <- overlapsAny(a.indel.gr,all.coding,ignore.strand = TRUE)
sum(the.coding.snps)
length(the.coding.snps)

############################
## the.coding.snps  is now a boolean if SNP add this to the selection process:
weights[1:5,]

# FOR gene modeling # 

## pass_coding<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  & 
##   ok.missing & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal & pass.possonian.control.model & !bad.qual.locations & the.coding.snps & !no.genotypes.cancer

## pass_non.coding<- full.qual &  !bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  &
##   ok.missing & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal & pass.possonian.control.model & !bad.qual.locations & the.coding.snps & !no.genotypes.cancer

#nocoalign pass.possonian.control.model & !bad.qual.locations not required
pass_coding<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  & 
  ok.missing & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal  & the.coding.snps & !no.genotypes.cancer

pass_non.coding<- full.qual &  !bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  &
  ok.missing & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal  & the.coding.snps & !no.genotypes.cancer

sum(pass_coding) # 47150 _> 31805 ->33366-> 33666
sum(pass_non.coding) # 21842 _>14990 -> 15675->15675 (recovered

## sum(pass_coding & no.genotypes.cancer)
## sum(pass_coding & !no.genotypes.cancer)
## summary.geno.extra[pass_coding & no.genotypes.cancer,][1:5,]

sort(table(a.indel[pass_coding,"Consequence.Embl"]),decreasing=TRUE)[1:20]

##                                             missense_variant                                                   stop_gained                       missense_variant,NMD_transcript_variant 
##                                                         40683                                                          2083                                                          1723 
##                 non_coding_exon_variant,nc_transcript_variant                        missense_variant,splice_region_variant                      splice_region_variant,synonymous_variant 
##                                                          1066                                                          1052                                                           502 
##                         frameshift_variant,feature_truncation                                              inframe_deletion                         frameshift_variant,feature_elongation 
##                                                           356                                                           225                                                           141 
##                            stop_gained,NMD_transcript_variant                                                intron_variant                                             inframe_insertion 
##                                                           114                                                            87                                                            74 
##                                           3_prime_UTR_variant                          splice_region_variant,intron_variant                             stop_gained,splice_region_variant 
##                                                            72                                                            65                                                            55 
## missense_variant,splice_region_variant,NMD_transcript_variant                                           5_prime_UTR_variant                                       initiator_codon_variant 
##                                                            52                                                            41                                                            39 
##                    3_prime_UTR_variant,NMD_transcript_variant                                     regulatory_region_variant 
##                                                            36

##                                                 missense_variant                                                      stop_gained                          missense_variant,NMD_transcript_variant 
##                                                            27928                                                             1782                                                             1135 
##                           missense_variant,splice_region_variant                    non_coding_exon_variant,nc_transcript_variant                         splice_region_variant,synonymous_variant 
##                                                              775                                                              696                                                              355 
##                               stop_gained,NMD_transcript_variant                                stop_gained,splice_region_variant                                                   intron_variant 
##                                                               82                                                               52                                                               51 
##                             splice_region_variant,intron_variant                                              3_prime_UTR_variant                            frameshift_variant,feature_truncation 
##                                                               48                                                               47                                                               33 
##    missense_variant,splice_region_variant,NMD_transcript_variant                                              5_prime_UTR_variant                       3_prime_UTR_variant,NMD_transcript_variant 
##                                                               33                                                               32                                                               28 
##                                          initiator_codon_variant                                            upstream_gene_variant                                               synonymous_variant 
##                                                               27                                                               26                                                               25 
## splice_region_variant,3_prime_UTR_variant,NMD_transcript_variant  splice_region_variant,synonymous_variant,NMD_transcript_variant 
                                                              19                                                               17


sort(table(a.indel[pass_non.coding,"Consequence.Embl"]),decreasing=TRUE)[1:20]



#### Troels I was concerned about the 5' and 3' entries in here even though they are small in number: these two lines check those ... they appear to be ok and are coging in thother transcripts
## Consequence.Embl is NOT a reliable list for  SNPs which are not traditionally damaging: (cause I report intron BEFORE synonymous is SNP is both!
#### one of refGene::location   knownGene::location ensGene::location   Consequence.Embl 
chk<-pass_coding & (a.indel[,"Consequence.Embl"] %in% c("5_prime_UTR_variant" ,"3_prime_UTR_variant,NMD_transcript_variant"))
a.indel[chk,c("chr","start","end","refGene::location","knownGene::location","ensGene::location","Consequence.Embl")]
#### THIS check idicated we are correct..
#############################################################

sort(table(a.indel[pass_coding,"Consequence.Embl"]),decreasing=TRUE)[1:20]

sort(table(a.indel[pass_non.coding,"refGene::location"]))  ### is if soe just REFSEQ
sort(table(a.indel[pass_non.coding,"knownGene::location"]))
sort(table(a.indel[pass_non.coding,"ensGene::location"]))

### like wise these look ok ( Consequence.Embl ) recodes the most 
chk<-pass_non.coding & (a.indel[,"Consequence.Embl"] %in% c("5_prime_UTR_variant" ,"intron_variant"))  
a.indel[chk,c("chr","start","end","refGene::location","knownGene::location","ensGene::location","Consequence.Embl","Gene.Names")][1:50,]



###################### SO I also do it though postions and use use the annova annotations.... lest check they are the same
synonymous<-test.for.coding.type(a.indel,geneanno.DB,synonymous.variant)
sum(synonymous) #17386

missense<-test.for.coding.type(a.indel,geneanno.DB,missense.variant)


## pass_non.coding.alternative<- full.qual &  !bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  &
##   ok.missing & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal & pass.possonian.control.model & !bad.qual.locations & synonymous

pass_non.coding.alternative<- full.qual &  !bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  &
  ok.missing & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal  & synonymous & !no.genotypes.cancer

sum(pass_non.coding.alternative)
sum(pass_non.coding)


(30916+ 1246+ 900)/15868


sort(table(a.indel[missense,"Consequence.Embl"]),decreasing=TRUE)[1:20]
sort(table(a.indel[synonymous,"Consequence.Embl"]),decreasing=TRUE)[1:20]

sort(table(a.indel[pass_non.coding.alternative,"Consequence.Embl"]),decreasing=TRUE)[1:20]
chk<-a.indel[,"Consequence.Embl"]=="non_coding_exon_variant,nc_transcript_variant"

ensGene::location


synon<-a.indel[,"ensGene::location"]=="synonymous SNV"
miss<-a.indel[,"ensGene::location"]=="nonsynonymous SNV"
sum(is.na(synon))

pass_coding<-full.qual  & not.flat.genotype & !in.any.normal  & !no.genotypes.cancer  & maf.filter  & !in.common.hit.gene   & !unannotated.hits  & !no.genotypes.filt & hw.controls.ok.filt & !in.any.normal.filt &  rare.in.Control.filt  & the.coding.snps &  bad.coding
pass_non.coding<- full.qual  & not.flat.genotype & !in.any.normal  & !no.genotypes.cancer  & maf.filter  & !in.common.hit.gene   & !unannotated.hits & !no.genotypes.filt & hw.controls.ok.filt & !in.any.normal.filt &  rare.in.Control.filt  & the.coding.snps &  !bad.coding

miss_coding<-full.qual  & not.flat.genotype & !in.any.normal  & !no.genotypes.cancer  & maf.filter  & !in.common.hit.gene   & !unannotated.hits  & !no.genotypes.filt & hw.controls.ok.filt & !in.any.normal.filt &  rare.in.Control.filt  & the.coding.snps &  miss

synon_coding<- full.qual  & not.flat.genotype & !in.any.normal  & !no.genotypes.cancer  & maf.filter  & !in.common.hit.gene   & !unannotated.hits & !no.genotypes.filt & hw.controls.ok.filt & !in.any.normal.filt &  rare.in.Control.filt  & the.coding.snps &  synon 


sum(miss_coding | synon_coding) # 50495

sum(miss_coding) #33071
sum(synon_coding) # 17424
sum(miss_coding)/sum(synon_coding) #1.898014

pass<- full.qual  & not.flat.genotype & !in.any.normal  & !no.genotypes.cancer  & maf.filter  & !in.common.hit.gene   & !unannotated.hits & !no.genotypes.filt & hw.controls.ok.filt & !in.any.normal.filt &  rare.in.Control.filt  & the.coding.snps & (miss_coding | synon_coding)

sum(pass) # 50495


miss.or.synon<-miss_coding | synon_coding
setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-03_LeoPharma_NovoAlign/Analysis")

       save(list=c("miss_coding","synon_coding","synonymous","no.genotypes.cancer","cellularity","cancer","PDs","the.samples.use","filt","snp.fail.filt","use.wieght","weights","core.ann","case.control","snpinfo.ori","pheno.types","ipheno","clusters.wanted","p","meta.results.skat","meta.results.skatO","meta.results.burden","pheno","pheno.ori","target.pheno.col","pass","a.indel","full.qual","bad.coding","bad.effect","maf.filter","in.common.hit.gene","on.x.y","unannotated.hits","not.flat.genotype","are.repeats","are.in.repeats","ok.missing","hw.controls.ok.filt","no.genotypes","rare.in.Control","rare.in.Control.filt","in.any.normal","in.any.normal.filt","are.in.repeats.back","are.in.repeats.forward","all.genes","contaminated","is.benign.missense","pass.possonian.control.model","bad.qual.locations"),file="Felicity_MugSIG_DATA.RData" )


sum(pass_coding) # 47150 _> 31805 ->33366-> 33666
sum(pass_non.coding) # 21842 _>14990 -> 15675->15675 (recovered
sum(pass_coding)/sum(pass_non.coding)





sort(table(a.indel[missense,"ensGene::location"]),decreasing=TRUE)[1:20]
sort(table(a.indel[synonymous,"ensGene::location"]),decreasing=TRUE)[1:20]



target<-missense & chk

target<-synonymous & chk & missense
sum(target)
a.indel[target,c(core.ann,"refGene::location","knownGene::location", "ensGene::location", "Consequence.Embl")][1:20,]

sum(chk)
sum(chk & synonymous)



################################ test specific genes 
## chk<- a.indel[,"Gene.Names"]=="MUC16" & pass_non.coding
## sum(chk)

## a.indel[chk,c(core.ann,"refGene::location","knownGene::location",
##               "ensGene::location", "Consequence.Embl","MAF.lt:0","MAF.lt:0.001","MAF.lt:0.005","MAF.lt:0.01","MAF.lt:0.025" )]
## a.chk<-c("chr19:8993452:8993452:A:G:snp","chr19:8995964:8995964:G:T:snp")
## a.chk %in% rownames(a.indel)
## mafs<-colnames(a.indel)[grepl("maf",colnames(a.indel))]
## a.indel[chk,c(core.ann,"refGene::location","knownGene::location",
##               "ensGene::location", "Consequence.Embl","MAF.lt:0","MAF.lt:0.001","MAF.lt:0.005","MAF.lt:0.01","MAF.lt:0.025",mafs )]
## maf.lt.all[chk,]
## maf.filter[chk]
## a.chk %in% rownames(a.indel)
## a.indel[a.chk,1:50]
## ## > sum(pass_non.coding.alternative)
## ## [1] 22989
## ## > sum(pass_non.coding)
## ## [1] 22704

sum(pass_non.coding.alternative != pass_non.coding)
#285 different
chk<-pass_non.coding.alternative != pass_non.coding
cbind(synonymous[chk],a.indel[chk,c("chr","start","end","refGene::location","knownGene::location","ensGene::location","Consequence.Embl")])[1:50,]
sort(table(a.indel[chk,"ensGene::location"]))
########## TRoles instecking these looks like using CDS lost a few 285 SYNs probaly due to diferences in CDS vs transcripts 
sort(table(a.indel[chk,"Gene.Names"]))
######## there only appear to be 1 or two in a few genes I think either approach is ok  ...

rownames(a.indel)[1:50]
C.T.snps<- grepl(":C:T:",rownames(a.indel))
G.A.snps<- grepl(":G:A:",rownames(a.indel))
sum(C.T.snps)

modif<-C.T.snps | G.A.snps




## target<-modif & pass_coding & a.indel[ ,"Gene.Names"] =="TTN"
## sort(tapply(a.indel[target ,"Consequence.Embl"],a.indel[target,"Consequence.Embl"],length),decreasing=TRUE)[1:10]

## target<-modif & pass_non.coding & a.indel[ ,"Gene.Names"] =="TTN"
## sort(tapply(a.indel[target ,"Consequence.Embl"],a.indel[target,"Consequence.Embl"],length),decreasing=TRUE)[1:10]

## inspect<-a.indel[ ,"Gene.Names"] =="TTN" & (a.indel[ ,"Consequence.Embl"] %in% c("intron_variant,nc_transcript_variant","intron_variant","non_coding_exon_variant,nc_transcript_variant"))
## a.indel[inspect,1:30]
## write.table(a.indel[inspect,1:30],file="TTN_anotation_complications.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


################################################



target<-modif & pass_coding 

pass_coding_count<-as.data.frame(tapply(a.indel[target ,"Gene.Names"],a.indel[target,"Gene.Names"],length))
pass_coding_count$Gene<-rownames(pass_coding_count)
colnames(pass_coding_count)<-c("NonSynonymous","Gene")

target.non<-modif & pass_non.coding 
pass_non.coding_count<-as.data.frame(tapply(a.indel[target.non,"Gene.Names"],a.indel[target.non,"Gene.Names"],length))
pass_non.coding_count$Gene<-rownames(pass_non.coding_count)
colnames(pass_non.coding_count)<-c("Synonymous","Gene")
df <- merge(pass_non.coding_count,pass_coding_count,by="Gene")
#weights<-df[df$Synonymous > 5 & df$NonSynonymous > 5,]
weights<-df
getwd()
write.table(weights,file="new_weights_final_MARIA_w_zeros_GENO.recovered.txt",sep="\t",row.names=F,col.names=T)
write.table(weights,file="new_weights_final_MARIA_w_zeros_GENO.NOTrecovered.txt",sep="\t",row.names=F,col.names=T)
weights[1:5,]






## load("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/APRIL_1st_no_containiation_snpFilters_weights.RData")
library(robust)
library(robustbase)
library(nls2)
library(akima)
## weights<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Felicity_filter/gene_synonymous_nonsynon_counts_above5.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)



## weights<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Felicity_filter/new_weights_final.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
## weights<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Felicity_filter/new_weights_final_MARIA.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
## weights<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Felicity_filter/new_weights_final_MARIA_w_zeros.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

#weights<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Felicity_filter/gene_synonymous_nonsynon_counts_all.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

weights<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-03_LeoPharma_NovoAlign/Analysis/new_weights_final_MARIA_w_zeros_GENO.NOTrecovered.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

weights<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-03_LeoPharma_NovoAlign/Analysis/new_weights_final_MARIA_w_zeros_GENO.recovered.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

weights[1:5,]

 weights[weights$Gene %in% c("TP53","TTN","BCL2L12","KNSTRN"),]
##    plot(
##   x,
##   y,
##   col=cols,
##   pch=20,
##   cex=3.0,
##   xlim=x.range,ylim=y.range,
##   ylab=expression(-log[10](italic(p))),
##   xlab="chromosome 6 position (Mb)",
##   bty='l',xaxt='n',cex.lab=3.0,cex.axis=3.0,font=2,font.lab=2)

## axis(1,at=seq(29e6,34e6,by=0.5e6),labels=seq(29,34,by=0.5),
##      cex.axis=3.0,font=2,font.lab=1)


# snps<-identify(x,y,label=t$ID[toplot],cex=1.5)
# t[toplot,][snps,]


par(mar=c(5.5,5.5,2,2),mgp=c(3,1,0))
            
## original.model<-use.wieght
#the.model.no.geno.recovery<-the.model

plot(as.numeric(weights[,2]) ,as.numeric(weights[,3]),main="For C>T motif SNPs",xlab="Synonymous Counts",ylab="Non-synonymousCounts",cex.lab=2.0,cex.axis=2.0,font=2,font.lab=2)
identify(as.numeric(weights[,2]) ,as.numeric(weights[,3]),labels=weights[,1],atpen=TRUE,font=2,cex=1.5)

## the.model<-try(lmrob(NonSynonymous~Synonymous,data=weights,silent=FALSE)) # ,subset=c(1:dim(test)[1])[!red])
##                            abline(coef(the.model),lty=10,col="green",lwd=2)

## abline(original.model,lty=10,col="purple",lwd=2)
## abline(use.wieght,lty=10,col="red",lwd=2)


the.model<-try(lm(NonSynonymous~Synonymous,data=weights,subset=c(1:dim(weights)[1])[weights$Synonymous >=5 & weights$NonSynonymous >=0 ],silent=FALSE) )
               abline(coef(the.model),lty=10,col="red",lwd=2)
use.wieght<-coef(the.model)
 #abline(a=0,b=coef(the.model),lty=10,col="black",lwd=2)


the.model<-try(lm(NonSynonymous~Synonymous,data=weights,subset=c(1:dim(weights)[1])[weights$Synonymous >=10 & weights$NonSynonymous >=10],silent=FALSE) )
abline(coef(the.model),lty=10,col="purple",lwd=2)



## the.model<-try(lmrob(NonSynonymous~Synonymous,data=weights,subset=c(1:dim(weights)[1])[!(weights$Gene %in% c("TTN"))],silent=FALSE) ) #
## abline(coef(the.model),lty=10,col="green",lwd=2)

## the.model<-try(lmrob(NonSynonymous~Synonymous,data=weights,subset=c(1:dim(weights)[1])[weights$Synonymous >=10],silent=FALSE) )
##                abline(coef(the.model),lty=10,col="red",lwd=2)

#use.wieght<-coef(the.model)
                ## abline(c(0,coef(the.model)[2]),lty=10,col="cyan",lwd=2)

## the.model<-try(lmrob(NonSynonymous~Synonymous-1,data=weights,subset=c(1:dim(weights)[1])[weights$Synonymous >=10],silent=FALSE) )
##                abline(c(0,coef(the.model)),lty=10,col="cyan",lwd=2)


the.model<-try(lm(NonSynonymous~Synonymous,data=weights,subset=c(1:dim(weights)[1])[(weights$Synonymous >=5 & weights$NonSynonymous >=0)  & !(weights$Gene %in% c("TTN","MUC16"))],silent=FALSE) )
 abline(coef(the.model),lty=10,col="blue",lwd=2)

the.model<-try(lm(NonSynonymous~Synonymous,data=weights,subset=c(1:dim(weights)[1])[(weights$Synonymous >=5 & weights$NonSynonymous >=0) & !(weights$Gene %in% c("TTN"))],silent=FALSE) )
 abline(coef(the.model),lty=10,col="forestgreen",lwd=2) 




leg.txt<-c("Fit No. Synon >5","Fit No. Synon/nonSynon >10","Fit No. Synon >5 no TTN","Fit No. Synon >5 no TTN,MUC16")
legend(0,170,leg.txt,col=c("red","purple","forestgreen","blue"),lty=c(2,2,2,2),pch=c(-1,-1,-1,-1),cex=2.0)


  ## abline(coef(the.model),lty=10,col="red",lwd=2)

savePlot(filename=paste("UV_wieght_model_final_recovery.png",sep=''),type="png")
savePlot(filename=paste("UV_wieght_model_final_recovery.jpeg",sep=''),type="jpeg")
savePlot(filename=paste("UV_wieght_model_final_recovery.tiff",sep=''),type="tiff")
#savePlot(filename=paste(fig.prefix,"aa.bmp",sep=''),type="bmp")
dev.print(svg,paste("UV_wieght_model_final_recovery.svg",sep=''))

savePlot(filename=paste("UV_wieght_model_final_No_recovery.png",sep=''),type="png")
savePlot(filename=paste("UV_wieght_model_final_No_recovery.jpeg",sep=''),type="jpeg")
savePlot(filename=paste("UV_wieght_model_final_No_recovery.tiff",sep=''),type="tiff")
#savePlot(filename=paste(fig.prefix,"aa.bmp",sep=''),type="bmp")
dev.print(svg,paste("UV_wieght_model_final_No_recovery.svg",sep=''))
getwd()

#weights<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Felicity_filter/gene_synonymous_nonsynon_counts_all.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
use.wieght
weights[1:5,]
min.val<-(1-use.wieght[1])/use.wieght[2]
min.val
the.weights<-weights$Synonymous*use.wieght[2] + use.wieght[1]
names(the.weights)<-weights$Gene ### use to weight the SNPS

##################### just checks that no negative weights
get.min.alleles<-cbind(weights,the.weights)
get.min.alleles[1:5,]

positive.weights<-get.min.alleles[,"the.weights"]>=0

get.min.alleles<-get.min.alleles[positive.weights,]
order.by<-order(get.min.alleles[,"the.weights"])
get.min.alleles[order.by,][1:20,]

min(get.min.alleles[positive.weights,"Synonymous"])
min(get.min.alleles[positive.weights,"NonSynonymous"])
#too.few.mutations<-weights$Synonymous<=5
plot(sort(the.weights))

max(the.weights)
min(the.weights)
##################################


the.weights<-1/the.weights


sum(weights$Synonymous<=5) # 6959-> 4191
#weights[weights$Synonymous<5,]

too.few.mutations<-(weights$Synonymous <=5) # can go to 3
sum(too.few.mutations) #5160
sum(!too.few.mutations) #212
#sum(weights[weights$Synonymous<5,"Gene"] %in% tight.weight.sig.genes)

#unweighted<-weights[too.few.mutations,"Gene"][(weights[too.few.mutations,"Gene"]%in% tight.weight.sig.genes)]
#the.weights[unweighted]

                           
## > the.weights[unweighted]
##       AFF2       ANO5    BCL2L12      BMPER   C12orf42    CD163L1      CDH10    COL11A2    COL19A1     COL1A2    COL24A1     COL4A1     COL4A2     COL4A4     CTNNA2      DCLK1       DGKK      FCRL5 
##  -9.596140  -6.714160  -9.596140 -12.478120 -12.478120  -3.832179  -6.714160  -3.832179 -12.478120 -12.478120  -9.596140  -6.714160  -9.596140  -3.832179 -12.478120 -12.478120 -12.478120  -6.714160 
##       FMN2       FSCB      GRIA2       HFM1      ITIH6      KCNQ3      LPPR4     LRRIQ1     LRRTM4      MYO16      NLRP4     OR4A15      PAPPA     PAPPA2      PLCB1      PREX2      PRKD1      PTPRC 
##  -3.832179  -9.596140  -6.714160  -9.596140  -6.714160  -6.714160 -12.478120  -9.596140  -6.714160 -12.478120  -3.832179  -9.596140  -3.832179  -3.832179  -6.714160 -12.478120 -12.478120  -3.832179 
##      RIMS1      SCN2A     SHANK2      TRHDE      TRPA1     UNC13C      ZBED9    ZNF804A 
##  -3.832179  -3.832179  -3.832179  -6.714160  -6.714160  -6.714160  -3.832179 -12.478120 
weights[!too.few.mutations,][1:10,]
weights[too.few.mutations,][1:10,]

max(the.weights[!too.few.mutations])
min(the.weights[!too.few.mutations])


dim(weights)
length(the.weights)

the.weights[1:5]
weights[1:5,]

sum(weights$Gene !=names(the.weights))

the.weights<-the.weights[!too.few.mutations]  ## useed in assoctaion test
weights<-weights[!too.few.mutations,]

max(the.weights)
min(the.weights)

## sort(the.weights)
plot(the.weights)
if("the.weights" %in% colnames(weights)){
  weights[,"the.weights"]<-the.weights
}else{
weights<-cbind(weights,the.weights)
}

sort(the.weights,decreasing=TRUE)[1:10]
sort(the.weights,decreasing=FALSE)[1:10]
weights[weights$Gene %in% c("TP53","BCL2L12","KNSTRN"),]

##         TTN       MUC16       DNAH5        MYH2        RYR2      AHNAK2        FAT3      ABCA13        MUC4       OBSCN        RYR1       CSMD3        PLEC        PCLO      ZNF469      DNAH10 
## 0.004585467 0.010008095 0.022142761 0.022142761 0.022142761 0.025382305 0.025382305 0.027385597 0.027385597 0.027385597 0.027385597 0.029732204 0.029732204 0.032518650 0.032518650 0.035881386 
##       HYDIN       MUC5B       USH2A       WDFY4       XIRP2      COL6A3       CSMD1        CUBN      DNAH17        GRM1       HERC2       KMT2D      MYO18B        PKD1       SYNE1       DNAH7 
## 0.035881386 0.035881386 0.035881386 0.035881386 0.035881386 0.040019812 0.040019812 0.040019812 0.040019812 0.040019812 0.040019812 0.040019812 0.040019812 0.040019812 0.040019812 0.045237318 
##         FN1       LRP1B       MYO5B        RELN        RYR3        TNXB       TRRAP        ANK2         DMD         EYS         FLG       IGFN1        LRP1        LRP2        MYH4        MYH6 
## 0.045237318 0.045237318 0.045237318 0.045237318 0.045237318 0.045237318 0.045237318 0.052019223 0.052019223 0.052019223 0.052019223 0.052019223 0.052019223 0.052019223 0.052019223 0.052019223 
##        MYH7       PRDM9       WDR87        ACAN     CACNA1G     CACNA1H     CACNA1S     CNTNAP4       GPR98       GREB1       MYO7B      NOTCH1        RTL1         VWF        APOB     CCDC168 
## 0.052019223 0.052019223 0.052019223 0.061193209 0.061193209 0.061193209 0.061193209 0.061193209 0.061193209 0.061193209 0.061193209 0.061193209 0.061193209 0.061193209 0.074295849 0.074295849 
##       DNAH1       DNAH3       DOCK8     DYNC1H1        FAT4       LAMA1       MACF1      RNF213       SALL1       SCN5A   SPATA31D1       ABCC6       ABCC8      ATP10A          C3       DCHS2 
## 0.074295849 0.074295849 0.074295849 0.074295849 0.074295849 0.074295849 0.074295849 0.074295849 0.074295849 0.074295849 0.074295849 0.094538314 0.094538314 0.094538314 0.094538314 0.094538314 
##       DNAH2       DNAH8       DNAH9       DOCK2        FASN       FRAS1        HRNR       HSPG2        MGAM        MYH1        MYLK     PKHD1L1      PRUNE2         RP1        SSPO       TTC28 
## 0.094538314 0.094538314 0.094538314 0.094538314 0.094538314 0.094538314 0.094538314 0.094538314 0.094538314 0.094538314 0.094538314 0.094538314 0.094538314 0.094538314 0.094538314 0.094538314 
##        VCAN       ZFHX3       ZFHX4    C20orf26      COL7A1         DCC      GRIN2B   KIAA1549L      KIF26B       LAMA5       MUC12        MYPN      NLRP12       NLRP5        PCNT      PIK3CG 
## 0.094538314 0.094538314 0.094538314 0.129941992 0.129941992 0.129941992 0.129941992 0.129941992 0.129941992 0.129941992 0.129941992 0.129941992 0.129941992 0.129941992 0.129941992 0.129941992 
##       SCN1A     SHROOM4        SPTB       TECTA          TG        TP53      TRANK1         ZAN       ZZEF1       ABCA3       ABCC9        ANK3     ANKRD11         CAD        CHD7     CNTNAP5 
## 0.129941992 0.129941992 0.129941992 0.129941992 0.129941992 0.129941992 0.129941992 0.129941992 0.129941992 0.207737759 0.207737759 0.207737759 0.207737759 0.207737759 0.207737759 0.207737759
##########################################################################
##########################################################################
##########################################################################

############################# rescue somatics with few calls with possion model

if(!exists("a.indel.ori")){ ### don't do revovery is alrready done!!
cellularity<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Sequenza/NSLM.cellularity.summary.use.csv",header=T,sep="\t",fill=TRUE,skip=0,stringsAsFactors=FALSE)

cellularity[1:5,]


PDs<-pheno.ori[pheno.ori[,"PD"],"SAMPLE"]
## PDs<-c(PDs,"LPH-001-27_PD")
PDs
#PDs.alt.counts<-alt.reads.reference.calls(a.indel,PDs,threshold=1)


cancer<-pheno.ori[pheno.ori[,"cancer"],"SAMPLE"]
cancer
#cancer.alt.counts<-alt.reads.reference.calls(a.indel,cancer,threshold=1)
#cancer.alt.counts.true<-alt.reads.Non.reference.calls(a.indel,cancer,threshold=1)
Controls<-pheno.ori[pheno.ori[,"Control"],"SAMPLE"]



## normal.alt.counts[test,]
## PDs.alt.counts[test,]
## cancer.alt.counts[c("chr13:21729956:21729956:A:G:snp","chr2:198363406:198363406:C:T:snp","chr22:20760282:20760282:A:C:snp"),]
## cancer.alt.counts.true[c("chr13:21729956:21729956:A:G:snp","chr2:198363406:198363406:C:T:snp","chr22:20760282:20760282:A:C:snp"),]



## normal.alt.counts[c("chr15:40675107:40675107:C:T:snp","chr19:50169131:50169131:C:T:snp","chr17:7578263:7578263:G:A:snp","chr17:7574003:7574003:G:A:snp"),]
## PDs.alt.counts[c("chr15:40675107:40675107:C:T:snp","chr19:50169131:50169131:C:T:snp","chr17:7578263:7578263:G:A:snp","chr17:7574003:7574003:G:A:snp"),]
## cancer.alt.counts[c("chr15:40675107:40675107:C:T:snp","chr19:50169131:50169131:C:T:snp","chr17:7578263:7578263:G:A:snp","chr17:7574003:7574003:G:A:snp"),]
## cancer.alt.counts.true[c("chr15:40675107:40675107:C:T:snp","chr19:50169131:50169131:C:T:snp","chr17:7578263:7578263:G:A:snp","chr17:7574003:7574003:G:A:snp"),]


test<-c("chr13:21729956:21729956:A:G:snp","chr2:198363406:198363406:C:T:snp","chr22:20760282:20760282:A:C:snp")
test<-c("chr19:50169131:50169131:C:T:snp","chr15:40675107:40675107:C:T:snp","chr6:31940123:31940123:G:A:snp","chr19:46520242:46520242:C:T:snp","chr12:103762692:103762692:G:A:snp","chr17:7578212:7578212:G:A:snp")


test<-c("chr12:11286221:11286221:T:A:snp","chr6:132029865:132029865:G:A:snp","chr7:82784807:82784807:C:A:snp","chr7:137206693:137206693:G:A:snp")

 test<-c("chr5:138268367:138268367:T:A:snp") #,"chr12:11286221:11286221:T:A:snp")
a.indel[test,paste(normals,".AD",sep="")]
a.indel[test,paste(controls,".GT",sep="")]
a.indel[test,paste(normals,".GT",sep="")]
a.indel[test,paste(PDs,".AD",sep="")]
a.indel[test,paste(PDs,".GT",sep="")]
a.indel[test,paste(cancer,".AD",sep="")]
a.indel[test,paste(cancer,".GT",sep="")]
#a.indel.ori[test,paste(cancer,".GT",sep="")]

#test<-c("chr15:40675107:40675107:C:T:snp") 
#cbind(a.indel[test,paste(cancer,".GT",sep="")],a.indel.ori[test,paste(cancer,".GT",sep="")])
Control.alt.counts[1:5,]
geno.p1<-genotype.p.values(a.indel.ori[test,],c(PDs),Control.alt.counts[test,"Read.Balance"]/100,cellularity)

geno.p1<-genotype.p.values(a.indel.ori[test,],c(cancer),Control.alt.counts[test,"Read.Balance"]/100,cellularity)

geno.p2<-genotype.p.values(a.indel[test,],c(PDs),normal.alt.counts[test,"Read.Balance"]/100,cellularity)


2*(pnorm(abs(c(4)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
p.threshold=0.0026 # z=3:   2*(pnorm(abs(c(1:8)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
p.threshold=2*(pnorm(abs(c(6)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
p.threshold
found.genotype<-  geno.p1<=p.threshold
apply(found.genotype,1,sum)

apply(found.genotype,1,sum)


rbind(a.indel.ori[test,paste0(c(PDs),".GT")][2,],
      a.indel[test,paste0(c(PDs),".GT")][2,],
a.indel[test,paste0(c(PDs),".AD")][2,],
      signif(geno.p1[2,],2),
      found.genotype[2,])


trial<-rbind(a.indel.ori[test,paste0(c(cancer),".GT")][2,],
      a.indel[test,paste0(c(cancer),".GT")][2,],
a.indel[test,paste0(c(cancer),".AD")][2,],
                signif(geno.p1[2,],2),
      found.genotype[2,])

rownames(trial)<-c("original","corrected","AD","p-value","revovered")

trial

order.by<-order(as.numeric(trial["p-value",]),decreasing=FALSE)

trial[,order.by]

geno.p[found.genotype]<-"0/1" # p.threshold<- 0.0026
geno.p[!found.genotype]<-"0/0"
colnames(geno.p)<-paste(colnames(geno.p),".GT",sep="")
colnames(genotypes)
geno.p[test,]
genotypes[test,colnames(geno.p)]



use.samples<-paste(c(cancer),".GT",sep="")
genotypes<-a.indel[test,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),"CANCER",sep=".")
summary.geno



use.samples<-paste(c(PDs),".GT",sep="")
genotypes<-a.indel[test,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),"PD",sep=".")
summary.geno

summary.geno.extra[test,grepl("GENO",colnames(summary.geno.extra))]
summary.geno.extra.ori[test,grepl("GENO",colnames(summary.geno.extra.ori))]

snp.only<-grepl("^snp",a.indel[,"TYPE"]) ### the p.values codes uses "AD" so does not work for SNPs

has.one.geno<-as.numeric(summary.geno.extra[,"ALT.Alleles.cancer"])>0

#geno.p<-genotype.p.values(a.indel[snp.only & has.one.geno ,] ,c(cancer,PDs),Control.alt.counts[snp.only & has.one.geno,"Read.Balance"]/100,cellularity) ## COULD USE NORMAL HERE
geno.p<-genotype.p.values(a.indel.ori[snp.only & has.one.geno ,] ,c(cancer,PDs),AD.extension="AD",Control.alt.counts[snp.only & has.one.geno,"Read.Balance"]/100,cellularity)
                           ## Control.alt.counts[1:5,]
                           ##  normal.alt.counts[1:5,]
#p.threshold=0.0026 # z=3:   2*(pnorm(abs(c(1:8)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))

p.threshold=2*(pnorm(abs(c(6)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
p.threshold
found.genotype<-  geno.p <= p.threshold
geno.p[found.genotype]<-"0/1" # p.threshold<- 0.0026
geno.p[!found.genotype]<-"0/0"
colnames(geno.p)<-paste(colnames(geno.p),".GT",sep="")


################################################################check normaml revovery
normal<-pheno.ori[pheno.ori[,"normal"],"SAMPLE"]
normal

####### set sample as pure
no.cellularity<-cbind(normal,1)
colnames(no.cellularity)<-c("Sample","Cellularity")
###########################################################
geno.p.n<-genotype.p.values(a.indel.ori[snp.only,] ,c(normal),AD.extension="AD",Control.alt.counts[snp.only,"Read.Balance"]/100,no.cellularity)

p.threshold=2*(pnorm(abs(c(6)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
p.threshold
found.genotype.n<-  geno.p.n <= p.threshold
geno.p.n[found.genotype.n]<-"0/1" # p.threshold<- 0.0026
geno.p.n[!found.genotype.n]<-"0/0"
colnames(geno.p.n)<-paste(colnames(geno.p.n),".GT",sep="")


posns<-match(rownames(a.indel),rownames(geno.p.n))
missing<-is.na(posns)


use.samples<-paste(normal,".GT",sep="")
print(use.samples)
length(use.samples)
genotypes<-geno.p.n[posns,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),"normal",sep=".")
summary.geno.normals.rec<-summary.geno

summary.geno.normals.rec[1:5,]

in.any.normal.rec<-as.numeric(summary.geno.normals.rec[,"ALT.Alleles.normal"]) > 0

sum(!in.any.normal)
sum(!in.any.normal.rec)

sum(in.any.normal.rec & !in.any.normal)

save(list=c("summary.geno.normals.rec","in.any.normal.rec","geno.p.n"),file="Normal.geno.recovered.RData")
setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-03_LeoPharma_NovoAlign/Analysis")

summary.geno.extra[!in.any.normal.filt,grepl("normal",colnames(summary.geno.extra))][1:5,]
sum(!in.any.normal)
# /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-03_LeoPharma_NovoAlign/Analysis/Normal.geno.recovered.RData
## colnames(genotypes)
## colnames(geno.p)

## test<-c("chr13:21729956:21729956:A:G:snp","chr2:198363406:198363406:C:T:snp","chr22:20760282:20760282:A:C:snp")
## test<-c("chr19:50169131:50169131:C:T:snp") # ,"chr15:40675107:40675107:C:T:snp","chr19:50169131:50169131:C:T:snp","chr17:7578263:7578263:G:A:snp","chr17:7574003:7574003:G:A:snp")
## test<-c("chr12:11286221:11286221:T:A:snp","chr6:132029865:132029865:G:A:snp","chr7:82784807:82784807:C:A:snp","chr7:137206693:137206693:G:A:snp")

## geno.p[test,]
## genotypes[test,colnames(geno.p)]

## dim(genotypes)
## dim(geno.p)
# summary.geno.extra.ori<- summary.geno.extra
#  a.indel.ori<-a.indel Only do this ONCE

sum(!( colnames(geno.p) %in% colnames(a.indel))) # must ge zero
to.transfer<-colnames(geno.p)[colnames(geno.p) %in% colnames(a.indel)]

posns<-match(rownames(a.indel),rownames(geno.p))
missing<-is.na(posns)
sum(!missing)
dim(geno.p)

a.indel[!missing,to.transfer]<-geno.p[posns[!missing],to.transfer]


## chk<-c("chr15:40675107:40675107:C:T:snp","chr19:50169131:50169131:C:T:snp") # BCL2L12
## loc<-key %in% chk
## targets<-the.projects  #c("NMD","ex.Control","AOGC")
## targets
## ###### the.samples and pheno in same order but the.samples has .GT extension.
## ## it<-1
## ## for(it in 1:length(targets)){
## ## use.samples<-the.samples[pheno[,targets[it]]]
## use.samples<-the.samples[pheno.ori[,"cancer"]]
## print(targets[it])
## print(use.samples)
## length(use.samples)
## genotypes<-a.indel.ori[loc,use.samples]
## dim(genotypes)
## summary.geno<-genotype.summary(as.matrix(genotypes))
## colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),names(targets)[it],sep=".")
## summary.geno
## }

## summary.geno.extra.ori[loc,c("MAF.cancer","GENO.cancer")]
## summary.geno.extra[loc,c("MAF.cancer","GENO.cancer")]

#summary.geno[1:5,]
###################################################
## a.indel is now the recovered genotypes
## a.indel.ori is the original


## "chr15:40675107:40675107:C:T:snp" %in% snp.fail.filt
## sum((snp.fail.filt %in% names(pass)))
## pass[ names(pass) %in% snp.fail.filt]<-FALSE

## pass.possonian.control.model<- poss.model[,"P-value"]> 1e-5  | is.na( poss.model[,"P-value"])

## length(pass.possonian.control.model)
## length(pass)

## ## sum(pass & !pass.possonian.control.model)
## ## snp.fail<-names(pass) %in% snp.fail.filt
## ## sum(snp.fail & !pass.possonian.control.model)
## ## sum(snp.fail)
## ## sum(!pass.possonian.control.model)

## pass[!pass.possonian.control.model]<-FALSE


## types<-c("all.somatic","coding.somatic","coding")
## types<-c("coding.somatic","coding.somatic.with.Indels","coding.somatic.with.Indels.noBenign","coding.somatic.with.Indels.noBenign.wRepeats","synonymous.with.Indels.noBenign.wRepeats","coding.somatic.with.repeats")
types<-c("coding.somatic.with.Indels.noBenign","Single Point","sliding.window")


itypes<-1

       
## for(itypes in 1:length(types)){

snap.file<-types[itypes]
snap.file
#snap.file<-paste(snap.file,"QC_no27.FINAL_6sd",sep="")
snap.file<-paste(snap.file,".withGQ.filter.recovery_weights_FULLQC_TIGHT_sd6_rare_recNOR",sep="")
snap.file<-paste(snap.file,".withGQ.filter.recovery_weights_FULLQC_TIGHT_sd6_rare",sep="")
snap.file<-paste(snap.file,".withGQ.filter.recovery_weights_TIGHT",sep="")
## snap.file<-paste(snap.file,".withGQ.filter.No_recovery_weights_TIGHT",sep="")
## snap.file<-paste(snap.file,".withGQ.filter.weights.lm.TTN_TIGHT_FULL_QUAL",sep="")
snap.file
## paste(snap.file,".withGQ.filter",sep="")

## paste(snap.file,".noGQ.filter",sep="")

## if(types[itypes]=="all.somatic"){
## pass<- full.qual &  bad.effect & maf.filter   & !in.common.hit.gene  & !unannotated.hits & not.flat.genotype & !are.repeats & !are.in.repeats &
## ok.missing & hw.controls.ok.filt & !no.genotypes &  rare.in.Control.filt & !in.any.normal.filt # & !on.x.y
## }


## if(types[itypes]=="all"){
## pass<- full.qual &  bad.effect & maf.filter   & !in.common.hit.gene  & !unannotated.hits & not.flat.genotype & !are.repeats & !are.in.repeats &
## ok.missing & hw.controls.ok.filt & !no.genotypes &  rare.in.Control.filt # & !on.x.y
## }

## if(types[itypes]=="non.coding"){
##   pass<- full.qual &  bad.non.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype & !are.repeats & !are.in.repeats &
## ok.missing.filt & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt #  & !on.x.y
## }

## if(types[itypes]=="coding.somatic.with.repeats"){
## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  &
## ok.missing & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal   &  snp.only # & !on.x.y
## }


## if(types[itypes]=="coding.somatic"){
  
## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype & !are.repeats & !are.in.repeats &
## ok.missing & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal   &  snp.only # & !on.x.y

## }

## if(types[itypes]=="coding.somatic.with.Indels"){
  
## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype & !are.repeats & !are.in.repeats &
## ok.missing & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal  # &!is.benign.missense  &  snp.only # & !on.x.y #  59753

## }
##  & !no.genotypes.cancer

## pass<- full.qual &  maf.filter  & not.flat.genotype  & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal & pass.possonian.control.model & !bad.qual.locations & !no.genotypes.cancer #  90263
## sum(pass) #133812

## pass<- full.qual  & not.flat.genotype & !in.any.normal  & !no.genotypes.cancer #  90263
## sum(pass) #133812


## pass_coding<-full.qual  & not.flat.genotype & !in.any.normal  & !no.genotypes.cancer &  bad.coding 
## pass_non.coding<- full.qual  & not.flat.genotype & !in.any.normal  & !no.genotypes.cancer &  !bad.coding 

## sum(pass_coding) # 47150 _> 31805 ->33366-> 33666
## sum(pass_non.coding) # 21842 _>14990 -> 15675->15675 (recovered
## sum(pass_coding)/sum(pass_non.coding)


## pass_coding<-full.qual  & not.flat.genotype & !in.any.normal  & !no.genotypes.cancer &  bad.coding & maf.filter  & !in.common.hit.gene   & !unannotated.hits  & !no.genotypes.filt & hw.controls.ok.filt & !in.any.normal.filt &  rare.in.Control.filt  & the.coding.snps 
## pass_non.coding<- full.qual  & not.flat.genotype & !in.any.normal  & !no.genotypes.cancer &  !bad.coding & maf.filter  & !in.common.hit.gene   & !unannotated.hits & !no.genotypes.filt & hw.controls.ok.filt & !in.any.normal.filt &  rare.in.Control.filt  & the.coding.snps 

## sum(pass_coding) # 47150 _> 31805 ->33366-> 33666
## sum(pass_non.coding) # 21842 _>14990 -> 15675->15675 (recovered
## sum(pass_coding)/sum(pass_non.coding)

## /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-03_LeoPharma_NovoAlign/Analysis/a.indel_for_felicity_filtering_quick.txt
## pass_coding<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  & 
##   ok.missing & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal  & the.coding.snps & !no.genotypes.cancer

## pass_non.coding<- full.qual &  !bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  &
##   ok.missing & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal  & the.coding.snps & !no.genotypes.cancer

## sum(pass_coding) # 47150 _> 31805 ->33366-> 33666
## sum(pass_non.coding) # 21842 _>14990 -> 15675->15675 (recovered /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/2015-03-16_AllAMLandLungHC.2.BEST.chrALL.ACC_SUBSET.ALL.ALL_GENOTYPES_analysis.txt


## wanted.genes<-unique(c(wanted.gene,c("NUDT18","BAD","ANKRD30BL,MIR663B","MIR663B","CDH11","PRKCG","NLRP4","NFKBIE","TNXB","ATN1","DEFB127","SPEF2","ZNF319","HYDIN","XIRP2","MMP16","STIP1")))
## sum(pass)
## pass<-pass & a.indel[,"Gene.Names"] %in% wanted.genes
## getwd()
## setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-03_LeoPharma_NovoAlign/Analysis")
## write.table(a.indel[pass,],file="a.indel_for_felicity_filtering_quick.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


pass<- full.qual &  maf.filter  & not.flat.genotype  & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal & pass.possonian.control.model & !bad.qual.locations
sum(pass)

sum(pass & in.any.normal.rec)

pass<- full.qual &  maf.filter  & not.flat.genotype  & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal & pass.possonian.control.model & !bad.qual.locations &  bad.coding

sum(pass)
sum(pass & in.any.normal.rec)

pass<- full.qual &  maf.filter  & not.flat.genotype  & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal & pass.possonian.control.model & !bad.qual.locations &  !bad.coding

sum(pass)
sum(pass & in.any.normal.rec)

## ## a.indel[!not.flat.genotype,core.ann][1:20,]

##  write.table(meta.results.burden,file=paste("Burden","ALL",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


if(types[itypes]=="Single Point"){
  
pass<- full.qual &  maf.filter  & not.flat.genotype  & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal & pass.possonian.control.model & !bad.qual.locations  #  90263

}

if(types[itypes]=="coding.somatic.with.Indels.noBenign"){
  
## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  &
## ok.missing & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal & !is.benign.missense & pass.possonian.control.model & !bad.qual.locations  #  43074

pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  &
ok.missing & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal & !is.benign.missense  & pass.possonian.control.model & !bad.qual.locations & !in.any.normal.rec

}

if(types[itypes]=="sliding.window"){
  
pass<- full.qual &  maf.filter  & not.flat.genotype  & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal & pass.possonian.control.model & !bad.qual.locations  #  90263
load("sliding_window_window.all.RData")
window.all[1:5,]
colnames(window.all)[1]<-"Name"

posns<-match(window.all[,"Name"],snpinfo.ori[,"Name"])
missing<-is.na(posns)
sum(missing)
window.all[!missing,"gene"]<-snpinfo.ori[posns[!missing],"gene"]
window.all[1:5,]
#snpinfo.ori.ori<-snpinfo.ori
rownames(window.all)<-window.all[,"Name"]
snpinfo.ori<-window.all


}


####
#snpinfo.ori.genes<-snpinfo.ori



## if(types[itypes]=="coding.somatic.with.Indels.noBenign.wRepeats"){
  
## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  &
## ok.missing & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal & !is.benign.missense & pass.possonian.control.model & !bad.qual.locations # &  &  snp.only # & !on.x.y #  42004

sum(pass)




## }

## if(types[itypes]=="synonymous.with.Indels.noBenign.wRepeats"){
  
## pass<- full.qual &  synonymous & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  &
## ok.missing & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !is.benign.missense # &  &  snp.only # & !on.x.y #  42004

## }



## if(types[itypes]=="coding"){
## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype & !are.repeats & !are.in.repeats &
## ok.missing & hw.controls.ok.filt & !no.genotypes &  rare.in.Control.filt #  & !on.x.y
## }

## pass<- full.qual & maf.filter     & !unannotated.hits & not.flat.genotype &
## ok.missing & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal # use for MutSig

sum(pass)
#45953

## elp<-cbind(pass,full.qual,bad.coding,maf.filter,rare.in.group,no.genotypes,in.common.hit.gene ,hw.controls.ok,on.x.y,unannotated.hits,not.flat.genotype,are.repeats,are.in.repeats,ok.missing,ok.missing,is.unwound.geno,is.unwound.geno ,hw.p.control.filt,is.benign.missense, pass.possonian.control.model,bad.qual.locations)


#

## pass<-a.indel[,"Gene.Names"] %in% c("TP53","NOTCH1","NOTCH2","ATM","ACD","ASIP","BAP1","CASP8","CCND1","CDK4","MC1R","MITF","MTAP","MX2","OCA2","PARP1","PLA2G6","POT1","SLC45A2","TERF2IP","TERT","TYR","TYRP1","VDR","BCL2L12","KNSTRN","ISX","CDKN2A","BCL2L11","STK19","FJX1","TRHDE")
## sum(pass)


## interesting.gene<-c("BCL2L12","CCDC61","STK19","KNSTRN","TRHDE","FREM2","EBNA1BP2","PHACTR3","DCLK1","LRRIQ1","PHACTR3","CSMD3")

## interesting.gene<-c("CDKN2A")
## pass<-  (a.indel[,"Gene.Names"] %in% interesting.gene)


## pass<- full.qual & (a.indel[,"Gene.Names"] %in% interesting.gene)  & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  &
## ok.missing & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal # & !is.benign.missense 
## sum(pass)
## snap.file<-"GENE_SYNON_TEST"



## help<-cbind(full.qual,bad.coding,maf.filter,rare.in.group,no.genotypes,in.common.hit.gene ,hw.controls.ok,on.x.y,unannotated.hits,not.flat.genotype,are.repeats,are.in.repeats,ok.missing,ok.missing,is.unwound.geno,is.unwound.geno ,hw.p.control.filt,rare.in.group.filt,no.genotypes.filt,rare.in.controls.filt )
## pass<-full.qual & functional & maf.filter & rare.in.group & !no.genotypes  & not.flat.genotype & !(high.total.missing | nimblegen.total.missing | illumina.total.missing)
## sum(pass)

sum(pass)

sum(as.numeric(summary.geno.extra[pass,"ALT.Alleles.normal"]))

hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.normal"])!=0
sum(hit & pass)
test<-pass & hit
## summary.geno.extra[test,]
## in.any.normal[test]
## summary.geno.extra[1:5,1:40]


## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  &
## ok.missing & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal & !is.benign.missense  & pass.possonian.control.model & !bad.qual.locations 

colnames(summary.geno.extra)
hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.cancer"])!=0
sum(hit & pass)  # 27582

hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.SCC"])!=0
sum(hit & pass)  #  8224 length(the.SCC)
sum(hit & pass)/length(the.SCC) #  373.8182
# > length(the.SCC)  [1] 22

hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.AK"])!=0
sum(hit & pass)  #  10912
sum(hit & pass)/length(the.AK)  # 242.4889
## > length(the.AK)[1] 45

hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.PD"])!=0
sum(hit & pass)  # 1853
sum(hit & pass)/25  #77.20833

test<-pass & hit
summary.geno.extra[test,][1:5,]
in.any.normal[test]
all.possible.samples

the.samples<-c(the.samples.use,pheno.ori[pheno.ori[,"PD"],"SAMPLE"])
genotypes<-a.indel[pass,paste(the.samples,".GT",sep="")] ## ordered correctly for phenotypes

genotypes[genotypes=="NA"]<-NA
genotypes[genotypes=="0/0"]<-0
genotypes[genotypes=="0/1"]<-1
genotypes[genotypes=="1/1"]<-2


dim(genotypes)
genotypes[is.na(genotypes)]<-0
dim(genotypes)
dim(genotypes)
colnames(genotypes)


#geno.per.sample<-apply(genotypes,2,function(x) sum(as.numeric(x)))

geno.per.sample<-apply(genotypes,2,function(x) sum(x[x!=0]))

names(geno.per.sample)<-gsub(".GT$","",names(geno.per.sample))
geno.per.sample
AKs<-grepl("_AK1$",names(geno.per.sample)) | grepl("_AK2$",names(geno.per.sample))
PDs<-grepl("_PD$",names(geno.per.sample))
SCCs<-grepl("_SCC$",names(geno.per.sample))


names(geno.per.sample)[AKs]
names(geno.per.sample)[PDs]
names(geno.per.sample)[SCCs]

sum(geno.per.sample[AKs])
sum(geno.per.sample[PDs])
sum(geno.per.sample[SCCs])


median(geno.per.sample[AKs])
median(geno.per.sample[PDs])
median(geno.per.sample[SCCs])

mean(geno.per.sample[AKs])
mean(geno.per.sample[PDs])
mean(geno.per.sample[SCCs])

sd(geno.per.sample[AKs])
sd(geno.per.sample[PDs])
sd(geno.per.sample[SCCs])

sort(geno.per.sample[AKs])
sort(geno.per.sample[PDs])
sort(geno.per.sample[SCCs])



(geno.per.sample)[AKs]
(geno.per.sample)[PDs]
(geno.per.sample)[SCCs]

t.test(geno.per.sample[SCCs],geno.per.sample[AKs])
hist(c(geno.per.sample[SCCs],geno.per.sample[AKs]))


counts<-geno.per.sample[SCCs | PDs | AKs]
counts
posns<-match(names(counts),cellularity[,"Sample"])

counts<-cbind(names(counts),names(counts),counts,cellularity[posns,"Cellularity"])
colnames(counts)<-c("Sample","Type","Counts","Cellularity")

coverage.file<-"/media/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_SAMPLE_coverage_LEO_2015.csv"
seq.type.file<-coverage.file  ## Troels-covergae prior to top up
seq.type<-read.delim(seq.type.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

seq.type[1:5,]
posns<-match(counts[,"Sample"],seq.type[,"Sample"])
counts<-cbind(counts,seq.type[posns,"percent.ccds.gt.30"])
colnames(counts)<-c(c("Sample","Type","Counts","Cellularity"),"Coverage")

counts[1:5,]

pheno1<-rep(NA,times=dim(counts)[1])

AKs<-grepl("_AK1$",counts[,"Sample"]) | grepl("_AK2$",counts[,"Sample"])
PDs<-grepl("_PD$",counts[,"Sample"])
SCCs<-grepl("_SCC$",counts[,"Sample"])
pheno1[AKs]<-"AK"
pheno1[SCCs]<-"SCC"
pheno1[PDs]<-"PD"
counts[,"Type"]<-pheno1
counts[1:5,]
## > counts[1:5,]
##                Sample           Type  Counts Cellularity Coverage       
## LPH-001-10_AK2 "LPH-001-10_AK2" "AK"  "247"  "0.36"      "85.2949501683"
## LPH-001-10_PD  "LPH-001-10_PD"  "PD"  "64"   "0.78"      "95.8668793316"
## LPH-001-10_SCC "LPH-001-10_SCC" "SCC" "1021" "0.37"      "89.7898658435"
## LPH-001-12_AK1 "LPH-001-12_AK1" "AK"  "171"  "0.33"      "90.6810583461"
## LPH-001-12_AK2 "LPH-001-12_AK2" "AK"  "425"  "0.28"      "91.1308513939"
                    

form<-"Type~Coverage+Cellularity+Counts"
data.test<-counts[counts[,"Type"]!="PD",]
data.test[data.test[,"Type"]=="SCC","Type"]<-1
data.test[data.test[,"Type"]=="AK","Type"]<-0

data.test<-as.data.frame(data.test)
data.test[1:5,]


form<-"Type~Counts"
form = as.formula(form)
form
fit1 = glm(as.numeric(as.character(data.test[,"Type"]))~ as.numeric(as.character(data.test[,"Counts"])) + as.numeric(as.character(data.test[,"Coverage"])) + as.numeric(as.character(data.test[,"Cellularity"])) ,family=binomial("logit"))
summary(fit1)

summary(fit1)
############## after revoery novo
## > hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.cancer"])!=0
## > sum(hit & pass)  # 27582
## [1] 30419
## > hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.SCC"])!=0
## sum(hit & pass)  #  8224 length(the.SCC)
## > sum(hit & pass)/length(the.SCC) #  373.8182
## [1] 9811
## > [1] 408.7917
## > hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.AK"])!=0
## sum(hit & pass)  #  10912
## sum(hit & pass)/length(the.AK)  # 242.4889
## > hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.PD"])!=0
## [1] 22783
## > [1] 455.66
## > sum(hit & pass)  # 1853
## sum(hit & pass)/25  #77.20833
## > test<-pass & hit
## [1] 1867
## > [1] 74.68

##### for ori no recovery novo
## > hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.cancer"])!=0
## > sum(hit & pass)  # 27582
## [1] 30387
## > hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.SCC"])!=0
## > sum(hit & pass)  #  8224 length(the.SCC)
## [1] 8826
## > sum(hit & pass)/length(the.SCC) #  373.8182
## [1] 367.75
## > hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.AK"])!=0
## > sum(hit & pass)  #  10912
## sum(hit & pass)/length(the.AK)  # 242.4889
## [1] 21947
## > [1] 438.94
## > hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.PD"])!=0
## sum(hit & pass)  # 1853
## > [1] 298
## > sum(hit & pass)/25  #77.20833
## [1] 11.92


## after geno recovery
## > hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.cancer"])!=0
## > sum(hit & pass)  # 27582
## [1] 30361
## > hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.SCC"])!=0
## sum(hit & pass)  #  8224 length(the.SCC)
## > [1] 9878
## > sum(hit & pass)/length(the.SCC) #  373.8182
## hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.AK"])!=0
## sum(hit & pass)  #  10912
## [1] 411.5833
## > > [1] 22769
## > sum(hit & pass)/length(the.AK)  # 242.4889
## [1] 455.38
## > hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.PD"])!=0
## > sum(hit & pass)  # 1853
## [1] 1891
## > sum(hit & pass)/25  #77.20833
## [1] 75.64

## after geno recovery full filtering
##  hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.cancer"])!=0
## > sum(hit & pass)  # 27582
## [1] 28768
## > hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.SCC"])!=0
## > sum(hit & pass)  #  8224 length(the.SCC)
## [1] 9311
## > sum(hit & pass)/length(the.SCC) #  373.8182
## [1] 387.9583
## > hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.AK"])!=0
## > sum(hit & pass)  #  10912
## [1] 21575
## > sum(hit & pass)/length(the.AK)  # 242.4889
## [1] 431.5
## > hit<-as.numeric(summary.geno.extra[,"ALT.Alleles.PD"])!=0
## > sum(hit & pass)  # 1853
## [1] 1745
## > sum(hit & pass)/25  #77.20833
## [1] 69.8



## a.indel[test,][,c("Gene.Names")]
## help[test,][1:3,]


chk<-"chr15:40675107:40675107:C:T:snp" # KNSTRN
chk<-"chr19:50169131:50169131:C:T:snp" # BCL2L12
chk<-"chr1:198723533:198723533:C:T:snp"
loc<-grep(chk,key)
#a.indel[,c(1:6,16,28,7,30,34,37:42,43,14,32,33)]
a.indel[loc,c(1:6,16,28,7,30,34,37:42,43,14,32,33)]
help[loc,]
pass[loc]
help[loc,]
summary.geno.extra[loc,]
summary.geno.extra.ori[loc,]
in.any.normal[loc]

## target.sample<-pheno[pheno[,"cancer"],"SAMPLE"]
## target.sample<-colnames(a.indel)[grep("_PD.GT",colnames(a.indel))]
## target.sample<-paste(target.sample,".GT",sep="")
## a.indel[loc,target.sample]

## inspect<-"LPH-001-9."
## inspect<-grepl(inspect,colnames(a.indel),fixed=TRUE)
## a.indel[loc,inspect]

## inspect<-"LPH-001-9"
## inspect<-grepl(inspect,colnames(a.indel),fixed=TRUE)
## a.indel[loc,inspect]


## sum(pass & not.flat.genotype)
## a.indel[(pass & not.flat.genotype),"TYPE"]

## ############ EXTRA FILTERING FROM FELICITY
## "chr15:40675107:40675107:C:T:snp" %in% snp.fail.filt
## sum((snp.fail.filt %in% names(pass)))
## pass[ names(pass) %in% snp.fail.filt]<-FALSE

## pass.possonian.control.model<- poss.model[,"P-value"]> 1e-5  | is.na( poss.model[,"P-value"])

## length(pass.possonian.control.model)
## length(pass)

## ## sum(pass & !pass.possonian.control.model)
## ## snp.fail<-names(pass) %in% snp.fail.filt
## ## sum(snp.fail & !pass.possonian.control.model)
## ## sum(snp.fail)
## ## sum(!pass.possonian.control.model)

## pass[!pass.possonian.control.model]<-FALSE



#################### FORCE GOOD CALLS

#pass[ names(pass) %in% "KNSTRN"]

## interesting<-snpinfo.ori[,"Name"] %in% extra.miRNA[,"Name"]
## #snpinfo.ori[interesting,"cluster"]<-snpinfo.ori[interesting,"gene"]
## sum(interesting)
## snpinfo.ori[interesting,]
## snpinfo.ori[interesting,"Name"]

## good.genotypes<-c("chr19:50169131:50169131:C:T:snp","chr19:50169104:50169104:C:T:snp","chr19:50169132:50169132:C:T:snp",
##                   "chr9:21971017:21971017:G:A:snp","chr9:21971016:21971016:G:A:snp","chr9:21971000:21971000:C:A:snp","chr9:21971056:21971056:C:T:snp","chr9:21971096:21971096:C:A:snp","chr9:21971099:21971099:G:A:snp","chr9:21971116:21971116:G:A:snp","chr9:21971015:21971015:C:A:snp","chr9:21994381:21994381:G:A:snp","chr9:21971040:21971040:C:T:snp","chr9:21994138:21994138:C:T:snp",snpinfo.ori[interesting,"Name"]) 

good.genotypes<-c("chr19:50169131:50169131:C:T:snp","chr19:50169104:50169104:C:T:snp","chr19:50169132:50169132:C:T:snp",  # BCL2L12
                  "chr9:21971017:21971017:G:A:snp","chr9:21971016:21971016:G:A:snp","chr9:21971000:21971000:C:A:snp","chr9:21971056:21971056:C:T:snp","chr9:21971096:21971096:C:A:snp","chr9:21971099:21971099:G:A:snp","chr9:21971116:21971116:G:A:snp","chr9:21971015:21971015:C:A:snp","chr9:21994381:21994381:G:A:snp","chr9:21971040:21971040:C:T:snp","chr9:21994138:21994138:C:T:snp", #CDKN2A
                  "chr7:1544063:1544063:G:A:snp","chr7:1544064:1544064:G:A:snp" # INTS1

                  ) ## chr9 are CDKN2A low coverge exone 1 and 2
good.genotypes %in% names(pass)
good.genotypes %in% snp.fail.filt
pass[ names(pass) %in% good.genotypes]<-TRUE

############## add weights
#the.weights ## weights for gene with more than 5 syn mutations

help<-cbind( pass,full.qual,bad.coding,maf.filter,in.common.hit.gene,unannotated.hits,not.flat.genotype,ok.missing,hw.controls.ok.filt,no.genotypes.filt,rare.in.Control.filt,in.any.normal.filt,in.any.normal,is.benign.missense,pass.possonian.control.model,bad.qual.locations )
#####################

snpinfo.ori[1:5,]
gene.weights<-rep(1,times=dim(snpinfo.ori)[1])
names(gene.weights)<-snpinfo.ori[,"gene"]
gene.weights[1:10]
posns<-match(names(gene.weights),names(the.weights))
missing<-is.na(posns)
sum(missing)
#cbind( names(gene.weights)[!missing],gene.weights[!missing],the.weights[posns[!missing]])[1:100,]
gene.weights[!missing]<-as.numeric(the.weights[posns[!missing]])
gene.weights["KNSTRN"]

C.T.snps<- grepl(":C:T:",snpinfo.ori[,"Name"])
G.A.snps<- grepl(":G:A:",snpinfo.ori[,"Name"])
sum(C.T.snps)

gene.weights[!(C.T.snps | G.A.snps)]<-1 ### don't weight non-UV SNPS

sum(grepl(":C:T:",snpinfo.ori[,"Name"]))

## > sum(grepl(":A:T:",snpinfo.ori[,"Name"]))
## [1] 4003
## > sum(grepl(":G:T:",snpinfo.ori[,"Name"]))
## [1] 8527
## > sum(grepl(":T:C:",snpinfo.ori[,"Name"]))
## [1] 14657
## > sum(grepl(":G:C:",snpinfo.ori[,"Name"]))
## [1] 7360
## > sum(grepl(":G:A:",snpinfo.ori[,"Name"]))
## [1] 60564
## > sum(grepl(":A:C:",snpinfo.ori[,"Name"]))
## [1] 4555
## > sum(grepl(":T:A:",snpinfo.ori[,"Name"]))
## [1] 4194
## > sum(grepl(":C:A:",snpinfo.ori[,"Name"]))
## [1] 20401
## > sum(grepl(":C:T:",snpinfo.ori[,"Name"]))
## [1] 60671
  
gene.weights[1:5]
snpinfo.ori[snpinfo.ori[,"gene"]=="KNSTRN",] 
#snpinfo[snpinfo[,"gene"]=="KNSTRN",] 
dim(snpinfo.ori)
length(gene.weights)

pass[unique(snpinfo.ori[snpinfo.ori[,"gene"]=="KNSTRN","Name"])]

sum(pass)



############### remove c->T uv
## a.C.T.snps<- grepl(":C:T:",key)
## a.G.A.snps<- grepl(":G:A:",key)
## uv.motif<-a.C.T.snps | a.G.A.snps
## sum(uv.motif)
## sum(pass & !uv.motif)
## sum(pass)
## pass[uv.motif]<-FALSE
################################# GEFOS FILTERING cause sending all
#pass<-pass[the.snps] ### GEOFS
## the.samples.use.no27<-the.samples.use[!(the.samples.use %in% c("LPH-001-27_AK1.GT","LPH-001-27_AK2.GT","LPH-001-27_PD.GT","LPH-001-27_SCC.GT"))]
## pheno<-pheno[pheno[,"SAMPLE"] %in% gsub(".GT$","",the.samples.use.no27),]
## dim(pheno)


############################ SCC vs AK

## #pheno.use<-pheno

## SCC<-rep(FALSE,times=dim(pheno)[1])
## SCC[grepl("_SCC$",pheno$SAMPLE)]<-TRUE

## AK<-rep(FALSE,times=dim(pheno)[1])
## AK[grepl("_AK1$",pheno$SAMPLE) | grepl("_AK2$",pheno$SAMPLE)]<-TRUE


## the.SCC<-pheno$SAMPLE[grepl("_SCC$",pheno$SAMPLE)]
## the.AK<-pheno$SAMPLE[grepl("_AK1$",pheno$SAMPLE) | grepl("_AK2$",pheno$SAMPLE)]


## sum(SCC)
## sum(AK)
## sum(pheno$PD)
## sum(pheno$SCC)
## sum(pheno$AK)
## table(pheno$AffectionStatus)
## pheno[1:5,]
## pheno$AffectionStatus<-0
## pheno[pheno$SCC,"AffectionStatus"]<-1
## pheno[pheno$AK,"AffectionStatus"]<-2
## the.samples.use<-paste(c(the.AK,the.SCC),".GT",sep="")
## pheno<-pheno[pheno[,"SAMPLE"] %in% gsub(".GT$","",the.samples.use),]
## the.samples.use.no27<-paste( pheno[,"SAMPLE"],".GT",sep="")

## the.samples.use.no27 %in% colnames(a.indel)
#############################

sum(pass)

genotypes<-a.indel[pass,the.samples.use] ## ordered correctly for phenotypes
snp.names<-key[pass] ## GEFOS ony name with start

#### snpinfo now A different size than a.indel since added pathways!!!  snpinfo[snpinfo[,"gene"]=="KNSTRN",]
snpinfo<-snpinfo.ori[snpinfo.ori[,"Name"] %in% snp.names,]
gene.weights.subset<-gene.weights[snpinfo.ori[,"Name"] %in% snp.names] # weight in same order as snpinfo.ori
snpinfo<-cbind(snpinfo,gene.weights.subset)
snpinfo[1:5,]
sum(is.na(as.numeric(snpinfo[,"gene.weights.subset"])))

###################################################



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

### single pont tests: cohort.seq.single<-cohort.seq meta.results.burden.single<-meta.results.burden
### single pont tests:
## snpinfo.ori[1:5,]
## cohort.seq <- skatCohort(Z=genotypes,formula, SNPInfo = snpinfo, data=pheno,aggregateBy="gene",family=binomial(),verbose=FALSE)
## meta.results.burden<-singlesnpMeta(cohort.seq,SNPInfo = snpinfo,aggregateBy="gene")
### single pont tests:
### single pont tests:
### single pont tests:

 if(target.pheno.col %in% case.control){
 cohort.seq <- skatCohort(Z=genotypes,formula, SNPInfo = snpinfo, data=pheno,aggregateBy="cluster",family=binomial(),verbose=FALSE)

}else{
cohort.seq <- skatCohort(Z=genotypes,formula, SNPInfo = snpinfo, data=pheno,aggregateBy="cluster",family=gaussian(),verbose=FALSE) ## genes and clusters
}

meta.results.burden<-burdenMeta(cohort.seq,wts="gene.weights.subset",mafRange = c(0,1),SNPInfo = snpinfo,aggregateBy="cluster")


# meta.results.burden<-burdenMeta(cohort.seq,wts=1,mafRange = c(0,1),SNPInfo = snpinfo,aggregateBy="cluster")
#meta.results.skat<-skatMeta(cohort.seq,SNPInfo = snpinfo,aggregateBy="cluster")
#meta.results.skatO<-skatOMeta(cohort.seq,burden.wts =1,SNPInfo = snpinfo,aggregateBy="cluster",method = "integration")
## print("start SKAT)")
## meta.results.skatO<-skatOMeta(cohort.seq,burden.wts ="gene.weights.subset",SNPInfo = snpinfo,aggregateBy="cluster",method = "integration")


## meta.results.skatO<-{}

the.order<-     order(meta.results.burden[,"p"])
sum(is.na(meta.results.burden[,"p"])) ## bad p-values shoudl not happen
meta.results.burden<-  meta.results.burden[the.order,]
meta.results.burden[1:50,]
## ## meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,] meta.results.burden.3sd<-meta.results.burden

meta.results.skat<-{}
meta.results.skatO<-{}
## ## the.order<-     order(meta.results.skat[,"p"])
## ## meta.results.skat<-  meta.results.skat[the.order,]
## ## meta.results.skat[1:50,]

the.order<-     order(meta.results.skatO[,"p"])
sum(is.na(meta.results.skatO[,"p"])) ## bad p-values shoudl not happen
meta.results.skatO<-  meta.results.skatO[the.order,]
meta.results.skatO[1:50,]

## ## meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,]
## ## meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]
## clusters[1:5,]

## meta.results.burden[meta.results.burden[,"gene"] %in% clusters[,1],]
## meta.results.burden[meta.results.burden[,"gene"] %in% clusters[,2],]


## meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,]
if(types[itypes]=="sliding.window"){
  posns<-match(meta.results.burden[,"gene"],snpinfo[,"cluster"])
  missing<-is.na(posns)
  the.gene<-snpinfo[posns,"gene"]
  meta.results.burden<-cbind(the.gene,meta.results.burden)
}
                           

## setwd(analysis.dir) meta.results.skat<-{}
## getwd()
## bad.non.coding
## snap.file<-"coding.0.01.all.geno.all.filters"
## snap.file<-"coding.0.001.all.geno.all.filters_no.imput"

 write.table(meta.results.burden,file=paste("Burden","ALL",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
 write.table(meta.results.skatO,file=paste("SKATO","ALL",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,],file=paste("Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


## write.table(meta.results.skat[1:50,],file=paste("Skat",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(meta.results.skatO[1:50,],file=paste("SkatO",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(meta.results.skatO[1:50,],file=paste("SkatO",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,],file=paste("SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


 annotations<-a.indel[,c(1:6,16,28,7,30,34,37:42,43,14,32,33)]

## save(list=c("meta.results.skat","meta.results.skatO","meta.results.burden","pheno.use","snpinfo","genotypes","pass","high.missing","annotations","help","key","summary.geno.extra"),file=paste(paste(project.files[ichr],".small.RData",sep="")) )




getwd()
save(list=c("synonymous","no.genotypes.cancer","cellularity","cancer","PDs","Control.alt.counts", "normal.alt.counts","the.samples.use","gene.weights","gene.weights.subset","filt","snp.fail.filt","use.wieght","weights","core.ann","case.control","snpinfo.ori","formula","clusters","pheno.types","ipheno","clusters.wanted","p","meta.results.skat","meta.results.skatO","meta.results.burden","pheno","pheno.ori","target.pheno.col","snpinfo","pass","high.missing.table","a.indel","a.indel.ori","help","key","summary.geno.extra","summary.geno.extra.ori","full.qual","bad.coding","bad.effect","maf.filter","in.common.hit.gene","on.x.y","unannotated.hits","not.flat.genotype","are.repeats","are.in.repeats","ok.missing","hw.controls.ok.filt","no.genotypes","rare.in.Control","rare.in.Control.filt","in.any.normal","in.any.normal.filt","are.in.repeats.back","are.in.repeats.forward","all.genes","contaminated","is.benign.missense","pass.possonian.control.model","bad.qual.locations"),file=paste(snap.file,".small_final.RData",sep="") )

} # itype

print(paste("Done: ",pheno.types[ipheno],"->",project.files[ichr]))
save.image(file=paste(project.files[ichr],".",pheno.types[ipheno],".",snap.file,".RData",sep=""))

} # pheno loop

} # loop over projects


} # loop over fam


 getwd()
[1] "/media/old-scratch/media/scratch2/AOGC-NGS/Analysis"
load("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/ALL_with_synon_Mar11_2015.RData")
load("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-03_LeoPharma_NovoAlign/Analysis/coding.somatic.with.Indels.noBenign.withGQ.filter.recovery_weights_FULLQC_TIGHT.small_final.RData")
colnames(a.indel)[1:50]

getwd()

setwd("/media/scratch/SCC_LOCAL")

save.image("FINAL_with_wights_lm_TTN_GENO_recover_UNK_JUN_16_novo.RData")
load("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-03_LeoPharma_NovoAlign/Analysis/FINAL_with_wights_lm_TTN_GENO_recover_UNK_JUN_16_novo.RData")

## save.image("FINAL_with_wights_6sd_rescure_JUN19.RData")
## save.image("FINAL_with_wights_6sd_resure_UNK_JUN.RData") 
## save.image("FINAL_with_wights_6sd_resure_UNK_27_APR.RData") ## @7th aprol has Maria weights  -15.36010     2.88198
## save.image("FINAL_with_wights_6sd_resure_UNK_14_MARR.RData")

###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
######################## YOU can re start the caulstion here / modify pass etc from here
## keep<-as.logical(a.indel[,"wanted.muts"])
hits<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/TOP_350.GENOTYPE.conponents-.Burden.clusters.coding.somatic.with.Indels.noBenign.mismatchfilter.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
hits[1:5,1:20]
hits<-unique(hits[,"Gene.Names"])

keep<-as.logical(a.indel[,"wanted.muts"]) | synonymous
sum(keep)

a.indel<-a.indel[keep,]

save(list=c("a.indel"),file="LEO_Pharma_wanted_with_synon_ALL_0.01_muts_QC_MAR31_a.indel.RData" )



## setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis")
## load("LEO_Pharma_wanted_with_synon_ALL_0.01_muts_Feb23")
## #save.image(file=paste("LEO_Pharma_WANTED_ALL_noGQ_filter_0.01_muts_Mar06",".RData",sep=""))

## setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-01-15_LeoPharma_Dec2014Freeze/Analysis")
## load("LEO_Pharma_WANTED_muts_Feb05.RData") ### by default the last run "coding"
## load("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/LEO_Pharma_ALL_0.01_muts_Feb23.RData")




## pass.coding.somatic.with.Indels<-pass.coding.somatic.with.Indels[keep]
## pass.coding.somatic.with.Indels.noBenign<-pass.coding.somatic.with.Indels.noBenign[keep]
## pass.coding.somatic.with.Indels.noBenign.wRepeats<-pass.coding.somatic.with.Indels.noBenign.wRepeats[keep]

## summary.geno.extra.keep<-summary.geno.extra[keep,]
## a.indel.keep<-a.indel[keep,]
## "a.indel.keep","summary.geno.extra.keep"
## save(list=c("pass.coding.somatic.with.Indels","pass.coding.somatic.with.Indels.noBenign","pass.coding.somatic.with.Indels.noBenign.wRepeats"),file="LEO_Pharma_WANTED_SUBSET2_0.01_muts_Feb23" )
##  rm("a.indel.keep")
##  rm("summary.geno.extra.keep")


###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
##################### RELOAD########################
## load("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/FINAL_with_wights_6sd_resure_UNK_27_APR.RData")

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
snap.file<-"coding_no_CT"


options(width=200)
the.order<-     order(meta.results.burden[,"p"])
sum(is.na(meta.results.burden[,"p"])) ## bad p-values shoudl not happen
meta.results.burden<-  meta.results.burden[the.order,]

meta.results.burden[1:50,]
## meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]

covars
target.pheno.col
#formula<-paste(target.pheno.col,"~",paste(covars,collapse="+"),sep="")
print(formula)
#formula<-formula(formula)
pheno[1:5,]

###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
## BCL2L12  KNSTRN     ISX   STK19  CDKN2A 
##      24     389     594     615     616 
#### this part breaks down the gene or pathway data 
## test<-c("MYCBP2","SLC25A24","TMCO3","C13orf35")
## test<-c("MYCBP2","SLC25A24","TMCO3","C13orf35")
## test<-c("chr22:41252435-41252687:ST13")
## test<-c("SETD8")
net<-c("BCL2L12","KNSTRN","ISX","CDKN2A","BCL2L11","STK19","FJX1","TTN","MUC4")
file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-01-15_LeoPharma_Dec2014Freeze/Analysis/String/INTS1_fill_network.txt"
file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-01-15_LeoPharma_Dec2014Freeze/Analysis/String/spring.ori.2.more.txt"

net<-read.delim(file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
dim(net)
net[1:5,]
net<-c(net[,1],net[,2])

sort(table(net),decreasing=TRUE)
net<-unique(net)
posns<-match(net,meta.results.burden[,"gene"])
names(posns)<-net
sort(posns)

posns<-match(net,meta.results.skatO[,"gene"])
names(posns)<-net
sort(posns)

meta.results.burden[meta.results.burden[,"gene"] %in% net,]

check<-grepl("^MED",meta.results.burden[,"gene"])
check<-grepl("^MED",meta.results.burden[,"gene"])


meta.results.burden[,"gene"][check]

posns<-grep("KNSTRN",a.indel[,16],as.is=TRUE)

a.indel[posns,1:20]


meta.results.burden[1:50,]
meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]
clusters.wanted<-gsub("^ ","",clusters.wanted)
########
########
the.top<-550
to.unwind<-c(meta.results.burden[1:the.top,"gene"])# ,meta.results.skatO[1:the.top,"gene"])



## to.unwind<-c("FGFR3") #, "MCM7", "RNPC3")
## interesting.gene<-c("BCL2L12","CCDC61","STK19","KNSTRN","TRHDE","FREM2","EBNA1BP2","PHACTR3","DCLK1","LRRIQ1","PHACTR3","CSMD3")
## ## to.unwind<-c("FANC_complex.all") # to.unwind<-meta.results.burden[8,"gene"]
##  to.unwind<-c("BCL2L12","CCDC61","STK19","KNSTRN","TRHDE","FREM2","EBNA1BP2","PHACTR3","DCLK1","LRRIQ1","PHACTR3","CSMD3") #"test.set"

## to.unwind<-c("TP53","NOTCH1","NOTCH2","FAT4","STK19","ISX","TRHDE","ARHGAP35","PREX1","KL","PIK3CA","KNSTRN","BCL2L12")
## to.unwind<-c("TP53","NOTCH1","NOTCH2","ATM","ACD","ASIP","BAP1","CASP8","CCND1","CDK4","MC1R","MITF","MTAP","MX2","OCA2","PARP1","PLA2G6","POT1","SLC45A2","TERF2IP","TERT","TYR","TYRP1","VDR","BCL2L12","KNSTRN","ISX","CDKN2A","BCL2L11","STK19","FJX1","TRHDE")

#to.unwind<-c("CDKN2A") #,"TCEA1","POLR2A","CTDP1")
## to.unwind<-c("Clinical")
## to.unwind<-c("Ubin.proteo","lipid_raft","caveolae","Citric")
#to.unwind<-c(clusters.wanted[!(clusters.wanted %in% c("Ubin.proteo","lipid_raft","caveolae","Checkpoint_extendedx1","Checkpoint_extendedx2"))])
#grep(to.unwind,meta.results.burden[,"gene"])
#to.unwind
# to.unwind.name<-to.unwind[1]
to.unwind.name<-"TOP_550_LENGTH_CONTROL"
#match(net,meta.results.burden[,"gene"])
# to.unwind.name<-"SYNON_test"
# to.unwind.name<-"Pathways"
# to.unwind.name<-"ALL_significant"
# to.unwind.name<-"ALL_significant"

snpinfo.ex<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,]
loci<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,"Name"] # this is IDH1 not IDH1 in cluster # are the snp.names
loci<-unique(loci)
the.genes<-unique(snpinfo.ex[,"cluster"])
the.genes<-unique(snpinfo.ex[,"gene"])
the.genes<-the.genes[!(the.genes %in% clusters.wanted)]

the.genes #245 ### if used a cluster name need to do back up to (**)

############repest to clean out cluster names 

the.genes.burden<-meta.results.burden[meta.results.burden[,"gene"] %in% the.genes,]

#the.genes.burden
write.table(the.genes.burden,file=paste(to.unwind.name,"conponents:","Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

the.genes.burden<-meta.results.skatO[meta.results.skatO[,"gene"] %in% the.genes,]
#the.genes.burden
write.table(the.genes.burden,file=paste(paste(to.unwind.name,collapse="."),"conponents:","SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



########### single point only
## length(loci)
## meta.results.burden[1:5,]
## loci<-meta.results.burden[1:550,"Name"]
###############


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


########### single point only
## snpinfo.ex<-snpinfo[snpinfo[,"Name"] %in% loci,]
## dim(snpinfo.ex)
## meta.results.burden.ex<-meta.results.burden[1:550,]
             ###############
# summary.geno.extra[loci,]
#high.missing[loci,]
#sum(are.in.repeats[loci])

    
# qual[loci,]
## snpinfo[1:5,]
## qual[1:5,c("FILTER_PASS", "FILTER_100" )]

cohort.seq.ex <- skatCohort(genotypes.ex, formula, SNPInfo = snpinfo.ex, data=pheno,aggregateBy = "Name",verbose=FALSE)
## meta.results.skat.ex<-skatMeta(cohort.seq,SNPInfo = snpinfo)
meta.results.burden.ex<-burdenMeta(cohort.seq.ex,wts=1,mafRange = c(0,1),SNPInfo = snpinfo.ex,aggregateBy = "Name")
meta.results.burden.ex[1:5,]
pheno[1:5,]
dim(meta.results.burden.ex)
dim(genotypes.ex)

cohort.seq.test <- skatCohort(genotypes.ex, formula, SNPInfo = snpinfo.ex, data=pheno,aggregateBy = "cluster",verbose=FALSE)

meta.results.burden.test<-burdenMeta(cohort.seq.test,wts=1,mafRange = c(0,1),SNPInfo = snpinfo.ex,aggregateBy = "cluster")
#meta.results.burden.test

## meta.results.skat.ex<-skatMeta(cohort.seq,SNPInfo = snpinfo)
#meta.results.skatO.test<-skatOMeta(cohort.seq.test,burden.wts =1,SNPInfo = snpinfo.ex,aggregateBy="cluster")
#meta.results.skatO.test
figure<- match(loci,key)

#genotypes.PD<-a.indel[figure, c("LPH-001-27_PD.GT",paste(pheno.ori[pheno.ori[,"PD"],"SAMPLE"],".GT",sep="")) ]
genotypes.PD<-a.indel[figure, c(paste(pheno.ori[pheno.ori[,"PD"],"SAMPLE"],".GT",sep="")) ]
genotypes.PD<-t(genotypes.PD)

genotypes.PD[genotypes.PD=="NA"]<-NA
genotypes.PD[genotypes.PD=="0/0"]<-0
genotypes.PD[genotypes.PD=="0/1"]<-1
genotypes.PD[genotypes.PD=="1/1"]<-2

rownames(genotypes.PD)<-gsub(".GT","",rownames(genotypes.PD))

dim(genotypes.ex)
dim(genotypes.PD)

options(max.print=200)
muts.in.PD<-apply(genotypes.PD,2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
muts.in.cases<-apply(genotypes.ex[pheno[,"cancer"],],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
muts.in.controls<-apply(genotypes.ex[pheno[,"Control"],],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})


 controls<- paste(pheno[pheno[,"Control"],"SAMPLE"],".GT",sep="")
## a.indel[figure,controls]
##   table(a.indel[figure,controls][2,])  

## muts.in.cases<-apply(genotypes.ex[pheno[,"SampleProject"]==1,],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
## muts.in.controls<-apply(genotypes.ex[pheno[,"SampleProject"]==0,],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})




########################################################
check<-16

quality.cases<-rep("",times=length(loci))
quality.controls<-rep("",times=length(loci))
quality.PD<-rep("",times=length(loci))

depth.cases<-rep("",times=length(loci))
depth.controls<-rep("",times=length(loci))
depth.PD<-rep("",times=length(loci))

a.indel.sub<-a.indel[figure,]
dim(a.indel.sub)
check<-1
for(check in 1:length(loci)){
# print(check)
#check<-"chr11:130066457:130066457:-:A:indel"
# posn<-grep(loci[check],key)
posn<-check


if(muts.in.PD[check]!=""){
#the.gt<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GT",sep=".")
the.gq<-paste(unlist(strsplit(muts.in.PD[check],split=",")),"GQ",sep=".")
the.ad<-paste(unlist(strsplit(muts.in.PD[check],split=",")),"AD",sep=".")

quality.PD[check]<-paste(a.indel.sub[posn,the.gq],collapse=",")
depth.PD[check]<-paste(a.indel.sub[posn,the.ad],collapse=";")

a.indel[posn,the.gq]
## a.indel[posn,the.gt]
## a.indel[posn,the.dp]
}




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
dim(a.indel)
dim(poss.model)
length(quality.cases)
length(figure)
dim(meta.results.burden.ex)
                             
#sum(meta.results.burden.ex[,"gene"]!=loci)
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
a.functions<-a.indel[,c("PolyPhen.scores","SIFT.scores","PolyPhen.desc","SIFT.desc")]

                             
out<-cbind(meta.results.burden.ex,poss.model[figure,],a.functions[figure,],annotations[figure,],is.benign.missense[figure],annotations[figure,],summary.geno.extra[figure,colnames(summary.geno.extra)[grep("^GENO",colnames(summary.geno.extra))]],help[figure,],muts.in.cases,quality.cases,depth.cases,muts.in.PD,quality.PD,depth.PD,muts.in.controls,quality.controls,depth.controls)

#all.data[figure,]
#out<-cbind(meta.results.burden.ex,annotations[figure,],muts.in.cases,muts.in.controls)
dim(out)
out[,1:13]



## table(out[,"refGene::location"])
## table(out[,"Consequence.Embl"])
getwd()
setwd(analysis.dir)
paste(paste(to.unwind,collapse="."))
paste(to.unwind.name,collapse=".")
  paste(paste(to.unwind.name,collapse="."),"GENOTYPE.conponents.","SkatO","clusters",snap.file,"txt",sep=".")

order.by<-order(out[,"p"],decreasing=FALSE)

out[order.by,][1:10,1:10]
write.table(out[order.by,],file=paste(paste(to.unwind.name,collapse="."),"GENOTYPE.conponents.","Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

getwd()
setwd(analysis.dir)

save(list=c("pheno.ori","pheno","a.indel","summary.geno.extra","use.wieght","pass.possonian.control.model","fil.genotypes"),file="No_revovery_model.RData")


#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################





viral<-("CDK4","CDK6","CDKN2A","CCND1","MDM2","KAT2B","EP300","RBPJ","CHD4","CDKN1B","RB1","ATM","CREB1")
hepB<-c("PRKCG","ATM","CDKN1A")


dim(genotypes)
genotypes[1:5,1:10]
sort(apply(genotypes,1,function(x) sum(x!=0,na.rm=TRUE)))
hist(sort(apply(genotypes,1,function(x) sum(x!=0,na.rm=TRUE))))














x<-out[,"muts.in.cases"]
x<-paste(x,collpase=",")
x<-gsub(" ,","",x)
x<-unlist(strsplit(x,split=","))
x<-x[x!=""]
x<-unique(x)

tapply(out[,"muts.in.cases"],out[,"Gene.Names"],function(x){
 x<-paste(x,collpase=",")
x<-gsub(" ,","",x)
x<-unlist(strsplit(x,split=","))
x<-x[x!=""]
x<-unique(x)
 x<-length(x)
}
       )


file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/paper/fanc snps for paul.txt"

filter<-read.delim(file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

dim(filter)
filter[1:5,]



filter[1:5]
################ Q-Qplot of data
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
## filt<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/TOP_350.GENOTYPE.conponents-.Burden.clusters.coding.somatic.with.Indels.noBenign.filter.update.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
## filt[1:5,]
## write.table(filt[1:350,],file="string.txt",,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


## bad<-filt[,"STRAND_BIAS"]=="true" |  filt[,"MISMATCH"]=="true"  | filt[,"READ_END"]=="true"
## tapply(filt[bad,"Gene.Names"],bad,length)
## bad.genes<-filt[bad,"Gene.Names"]
## bad.genes
## table(bad.genes)
## bad

       
## interesting<-c("TP53","NOTCH1","NOTCH2","FAT4","STK19","ISX","TRHDE","ARHGAP35","PREX1","KL","PIK3CA","KNSTRN","BCL2L12")

## interesting[interesting %in% bad.genes]
## bad.genes<-bad.genes[!(bad.genes %in% interesting)] # checked this is ok as were just single point mutations


## other.bad<-c("DNAH5","CSMD3","CSMD1","PCLO","INTS1","TYRO3","TTN","TESK1","MUC17","OR4A15","OR6C1","BCL2L11")

## other.bad %in% bad.genes

## bad.genes<-c(bad.genes,other.bad)

## meta.results.burden[1:50,]
## meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]

## the.order<-     order(meta.results.skat[,"p"])
## meta.results.skat<-  meta.results.skat[the.order,]
## meta.results.skat[1:50,]

## the.order<-     order(meta.results.skatO[,"p"])
## sum(is.na(meta.results.skatO[,"p"])) ## bad p-values shoudl not happen
## meta.results.skatO<-  meta.results.skatO[the.order,]
## meta.results.skatO[1:50,]

       

## dim(clusters)
## snpinfo.sub<-snpinfo[snpinfo[,"cluster"] %in% clusters.wanted,]
## genes.cl<-unique(snpinfo.sub[,"gene"])
## genes.cl<-genes.cl[!(genes.cl %in% clusters.wanted)]
## genes.cl
## clusters.wanted

## genes.and.clusters<-c(genes.cl,clusters.wanted)
## meta.results.burden[1:5,]

#################################### just want to plot

## subset<-meta.results.burden[ !(meta.results.burden[,"gene"] %in% genes.and.clusters) ,]
/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Burden.ALL.coding.somatic.with.Indels.noBenign.wRepeatsno27.final_6sd.xlsx
/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/TOP_550_LENGTH_CONTROL.GENOTYPE.conponents:.Burden.no27.final_6sd.txt
/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/TOP_550_LENGTH_CONTROL.GENOTYPE.conponents:.Burden.no27.final_6sd.txt
file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/TOP_350_LENGTH_CONTROL.GENOTYPE.conponents_.Burden.clusters.coding.somatic.with.Indels.noBenign.wRepeats.igv.csv"
file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Burden.ALL.coding.somatic.with.Indels.noBenign.wRepeats.withGQ.filter.QCed_weights.txt"

file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Burden.ALL.coding.somatic.with.Indels.noBenign.wRepeatsQC_no27.FINAL_6sd_NO_wieights.txt"
file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Burden.ALL.coding.somatic.with.Indels.noBenign.wRepeatsno27.final_6sd.txt"

file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Burden.ALL.coding.somatic.with.Indels.noBenign.wRepeats.withGQ.filter.QCed_weights_TIGHT.txt"


meta.results.burden<-read.delim(file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)




## posns<-match(meta.results.burden[,"gene"],rownames(poss.model))
## meta.results.burden<-cbind(poss.model[posns,],meta.results.burden)
## meta.results.burden[1:5,1:10]
## write.table(meta.results.burden,file="containamination loci.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)




clusters.wanted<-c("SCC_Genes"," BRCA.Genes"," RAD51.Paralogues"," BLM.Complex"," Checkpoint.Proteins"," citric"," lipid_raft"," caveolae"," Checkpoint_extendedx1.noP53"," Checkpoint_extendedx2.noP53"," Citric.no.IDH"," BLM.Complex_AND_Checkpoint"," FANCD2_minimal_mono_ubi"," MLH_cluster")

dups<-duplicated(meta.results.burden[,"gene"])
meta.results.burden[grep("RYR",meta.results.burden[,"gene"]),]




sum(dups)

file
meta.results.burden[1:5,]
       dim(meta.results.burden)
## bad.genes<-c(bad.genes,clusters.wanted)
subset<-meta.results.burden[ !(meta.results.burden[,"gene"] %in%  clusters.wanted) ,]
       subset[1:10,]
subset<-subset[!is.na(as.numeric(subset[ ,"p"])),]

## subset[1:10,]
## no.weight.sig.genes<-subset[as.numeric(subset[ ,"p"])<1.2e-6,"gene"]
## tight.weight.sig.genes<-subset[as.numeric(subset[ ,"p"])<1.2e-6,"gene"]

symbols<-subset[,"gene"]

## subset[1:10,]
## posns<-match(interesting,subset[,"gene"])
## names(posns)<-interesting
## posns
## write.table(subset[1:350,],file="string.txt",,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## dups<-duplicated(subset[,"gene"])
## sum(dups)

## write.table(subset[1:350,],file="string.txt",,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
z1<-qchisq(subset[ ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)
z<-qchisq(subset[ ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) ## if have no chisq valuse
#z0<-qchisq(meta.results.skatO[ !(meta.results.skatO[,"gene"] %in% genes.and.clusters ) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) 
z[1:5]
subset[1:10,]

p.value<-as.numeric(subset[ ,"p"])
par(mfrow=c(1,1),font=2,font.lab=2,font.axis=2,mgp=c(3.5,1,0),mar=c(5,5.5,4,2)+0.1)


##  z<-rchisq(length(p.val), df=1, ncp = 0) ## null test


median(z,na.rm=TRUE)/0.456  #1.071491
#median(z0,na.rm=TRUE)/0.456  #1.071491


sum(is.na(p.val))
################## p-values
z0=qnorm(p.value[!is.na(p.value)]/2)
lambda = round(median(z0^2)/0.454,3)
lambda


plot(x=c(1:10))
## source("http://bioconductor.org/biocLite.R") 
## biocLite("GWASTools")
## setRep
## qq.Plot(pvals)


help[1:5,]
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
signif(range(z),digits=3)

z[!is.na(z)]
the.plot<-my.qq.plot(z,dist="chisq",df=1,ncp=0,col="blue",ylab="Observed chi-squared value",xlab="Expected chi-squared value",main="",cex=1.5,xlim=c(5,18),ylim=c(20,100),cex.lab=3.0,cex.axis=2.0,font.lab=2,font.axis=2,lwd=3,line="robust",plot.it=TRUE) # function defined below


## z.all<-qchisq(meta.results.burden[ !(meta.results.burden[,"gene"] %in% clusters.wanted) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)
## range(z.all)
qq<-  qq.data(z,plot.it=FALSE)       ## qq plot used same method as in car library
#points(qq$x,qq$y,col="magenta",pch=21)

symbols<-subset[,"gene"]


interesting<-c("TP53","NOTCH1","NOTCH2","FAT1","DGKI","COL19A1") # ,"COL19A1"
interesting<-c("TP53","NOTCH1","NOTCH2","FAT1") # talk
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="magenta",pch=20,cex=2)
text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="magenta",pos=4,offset=1,cex=2.25)  #,pch=25)

interesting<-c("STK19","KNSTRN","BCL2L12","CCDC61","EBNA1BP2","PHACTR3","CDKN2A")
interesting<-c("STK19","KNSTRN","CDKN2A","DGKI")
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="forestgreen",pch=20,cex=2)
text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="forestgreen",pos=2,offset=1,cex=2.25)  #,pch=25)



interesting<-c("C12orf42")
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="forestgreen",pch=20,cex=2)
text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="forestgreen",pos=4,offset=2,cex=1.3)  #,pch=25)

interesting<-c("TP53","NOTCH1","NOTCH2","CDKN2A","KNSTRN")
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="red",pch=25,cex=2,lwd=2)
#text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="forestgreen",pos=2,offset=1,cex=1.3)


interesting<-c("TRHDE","PAPPA")
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col=c("forestgreen","magenta"),pch=20,cex=2)
#points(qq$x[show],qq$y[show],col=c("cyan"),pch=20,cex=2)
#text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col=c("forestgreen","magenta"),pos=4,offset=1,cex=1.3)  #,pch=25)
selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="forestgreen",cex=1.3,offset=1,atpen='TRUE')
selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="magenta",cex=1.3,offset=1,atpen='TRUE')





abline(h=qchisq(1.6e-6,df=1,ncp=0,lower.tail=FALSE,log.p=FALSE),col="green",lw=2)

leg.txt<-c(expression(paste("Genes")),"95% confidence intervals","Genes with multiple mutations","Genes with focal mutations")
legend(6,80,leg.txt,col=c("blue","red","magenta","forestgreen"),lty=c(-1,2,-1,-1),pch=c(1,-1,19,19),cex=2.0)

setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Figures")
fig.prefix<-"Q-Q_Final_6sd_tight_weights"
fig.prefix<-"Q-Q_Final_conference"


savePlot(filename=paste(fig.prefix,".png",sep=''),type="png")
savePlot(filename=paste(fig.prefix,".jpeg",sep=''),type="jpeg")
savePlot(filename=paste(fig.prefix,".tiff",sep=''),type="tiff")
savePlot(filename=paste(fig.prefix,".bmp",sep=''),type="bmp")
dev.print(svg,paste(fig.prefix,".svg",sep=''))



# FREM2 -check why it has moved about
# TNXB has pseudo genes and may be dodgy not caused by rescue but in redo with 27
# UNC80 has moved about chr2:210704125:210704125:C:T:snp


#####annotate curve
selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="red",cex=1,offset=1,atpen='TRUE') ##plate row col symbol

selected.data<-identify(qq$x,qq$y,labels=qq$ord,col="red",cex=1,offset=1,atpen='TRUE')


selected.data<-identify(qq$x,qq$y,labels=labels[qq$ord],col="red",cex=1,atpen='TRUE') ## sybmol
selected.data<-identify(qq$x,qq$y,labels=as.character(round(data.in[qq$ord],2)),col="forestgreen",cex=1.25,atpen='TRUE') # observed score
#####


############################################ WIGHT MODEL EFFEECTS
file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Burden.ALL.coding.somatic.with.Indels.noBenign.wRepeatsQC_no27.FINAL_6sd_NO_wieights.txt"
file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Burden.ALL.coding.somatic.with.Indels.noBenign.wRepeats.withGQ.filter.QCed_weights_TIGHT.txt"

prev.study.genes<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/KNSTRN_PUB_Resequenced_genes.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

prev.study.genes<-prev.study.genes[,1]

meta.results.burden<-read.delim(file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
clusters.wanted<-c("SCC_Genes"," BRCA.Genes"," RAD51.Paralogues"," BLM.Complex"," Checkpoint.Proteins"," citric"," lipid_raft"," caveolae"," Checkpoint_extendedx1.noP53"," Checkpoint_extendedx2.noP53"," Citric.no.IDH"," BLM.Complex_AND_Checkpoint"," FANCD2_minimal_mono_ubi"," MLH_cluster","Citric.no.IDH")

subset<-meta.results.burden[ !(meta.results.burden[,"gene"] %in%  clusters.wanted) ,]
subset[1:10,]
subset<-subset[!is.na(as.numeric(subset[ ,"p"])),]

## subset[1:10,]
## no.weight.sig.genes<-subset[as.numeric(subset[ ,"p"])<1.2e-6,"gene"]
## tight.weight.sig.genes<-subset[as.numeric(subset[ ,"p"])<1.2e-6,"gene"]
## no.weight.sig.genes<-subset[as.numeric(subset[ ,"p"])<1.2e-6,"gene"]
## tight.weight.sig.genes<-subset[as.numeric(subset[ ,"p"])<1.2e-6,"gene"]

change<-no.weight.sig.genes[!(no.weight.sig.genes %in% tight.weight.sig.genes)]

 ## [1] "TTN"           "RYR2"          "DMD"           "CCDC168"       "MUC16"         "PCLO"          "ABCA13"        "DNAH5"         "SYNE1"         "RELN"          "CSMD1"         "CUBN"         
## [13] "MGAM"          "FAT4"          "DOCK2"         "ZFHX4"         "LRP2"           "PKHD1L1"       "SPHKAP"        "EYS"           "RYR3"          "SCN1A"         "GPR98"        
## [25] "HYDIN"         "XIRP2"         "CNTNAP4"       "CNTNAP5"       "USH2A"         "DNAH8"         "APOB"          "VCAN"          "ANK2"          "FAT3"          "RYR1"          "MUC5B"        
## [37] "SCN9A"         "TLR4"          "GRIN2A"        "DNAH10"              "DCC"           "IGFN1"         "SCN7A"         "ABCC9"         "FREM2"
length(tight.weight.sig.genes)/length(no.weight.sig.genes) #0.65  about 35% of genes removed 

length(tight.weight.sig.genes)
length(no.weight.sig.genes)

diff<-tight.weight.sig.genes[!(tight.weight.sig.genes %in% no.weight.sig.genes)] ## No new genes are discovered 

z<-qchisq(subset[ ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)



the.plot<-my.qq.plot(z,dist="chisq",df=1,ncp=0,col="blue",ylab="Observed chi-squared value",xlab="Expected chi-squared value",main="",xlim=c(-1,18),ylim=c(-10,100),cex=1,cex.lab=2.0,cex.axis=2.0,font.lab=2,font.axis=2,lwd=2,line="robust",plot.it=TRUE) # function defined below



qq<-  qq.data(z,plot.it=FALSE)       ## qq plot used same method as in car library
#points(qq$x,qq$y,col="magenta",pch=21)

symbols<-subset[,"gene"]


interesting<-no.weight.sig.genes # ,"COL19A1"
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="magenta",pch=20,cex=2)
#text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="magenta",pos=4,offset=1,cex=1.3)  #,pch=25)

c("TTN","DNAH10","DNAH5","MUC4","OBSCN","MYH2", "APOB","MUC16","RYR2","AHNAK2","FAT3","RYR1","PLEC","PCLO","ZNF469")

long.genes<-c("MUC16","RYR2","PCLO","APOB", "ZNF521", "FAT4", "MUC5B" )
interesting<-long.genes
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="black",pch=22,cex=1.3)
text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col=c("red"),pos=2,offset=1.5,cex=1.3)  #,pch=25)

long.genes<-c("TTN","DNAH10","MUC4","OBSCN","FAT3")
interesting<-long.genes
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="black",pch=22,cex=1.3)
text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col=c("red"),pos=4,offset=1,cex=1.3)  #,pch=25)


## interesting<-prev.study.genes
## show<-symbols[qq$ord] %in% interesting
## points(qq$x[show],qq$y[show],col="black",pch=20,cex=1.5)
## #text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col=c("red"),pos=1,offset=1,cex=1.3)  #,pch=25)

## prev.study.genes[prev.study.genes %in% tight.weight.sig.genes]
##  ## [1] "KNSTRN" "CDKN2A" "COL1A2" "COL4A1" "COL4A2" "CSMD3"  "DCLK1"  "FAT1"  
##  ## [9] "FLNA"   "NOTCH1" "NOTCH2" "TP53"

## interesting<-prev.study.genes[prev.study.genes %in% tight.weight.sig.genes]

## "KNSTRN" "CDKN2A" "COL1A2" "COL4A1" "COL4A2" "CSMD3"  "DCLK1"  "FAT1"  
##  [9] "FLNA"   "NOTCH1" "NOTCH2" "TP53"  
## show<-symbols[qq$ord] %in% interesting
## #points(qq$x[show],qq$y[show],col="black",pch=20,cex=2.5)
## text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col=c("red"),pos=2,offset=1,cex=1.3) 

selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="magenta",cex=1.3,offset=1,atpen='TRUE')

abline(h=qchisq(1.6e-6,df=1,ncp=0,lower.tail=FALSE,log.p=FALSE),col="green",lw=2)

leg.txt<-c(expression(paste("Genes")),"95% confidence intervals","Significant genes with No UV weights")
legend(0,80,leg.txt,col=c("blue","red","magenta"),lty=c(-1,2,-1),pch=c(1,-1,19),cex=2.0)

setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Figures")
fig.prefix<-"Q-Q_Final_effect_of_weights"


savePlot(filename=paste(fig.prefix,".png",sep=''),type="png")
savePlot(filename=paste(fig.prefix,".jpeg",sep=''),type="jpeg")
savePlot(filename=paste(fig.prefix,".tiff",sep=''),type="tiff")
savePlot(filename=paste(fig.prefix,".bmp",sep=''),type="bmp")
dev.print(svg,paste(fig.prefix,".svg",sep=''))


############################## prev .stidy


the.plot<-my.qq.plot(z,dist="chisq",df=1,ncp=0,col="blue",ylab="Observed chi-squared value",xlab="Expected chi-squared value",main="",xlim=c(-1,18),ylim=c(-10,100),cex=1,cex.lab=2.0,cex.axis=2.0,font.lab=2,font.axis=2,lwd=2,line="robust",plot.it=TRUE) # function defined below



qq<-  qq.data(z,plot.it=FALSE)       ## qq plot used same method as in car library
#points(qq$x,qq$y,col="magenta",pch=21)

symbols<-subset[,"gene"]


interesting<-no.weight.sig.genes # ,"COL19A1"
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="magenta",pch=20,cex=2)
#text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="magenta",pos=4,offset=1,cex=1.3)  #,pch=25)

c("TTN","DNAH10","DNAH5","MUC4","OBSCN","MYH2", "APOB","MUC16","RYR2","AHNAK2","FAT3","RYR1","PLEC","PCLO","ZNF469")

## long.genes<-c("MUC16","RYR2","PCLO","APOB", "ZNF521", "FAT4", "MUC5B" )
## interesting<-long.genes
## show<-symbols[qq$ord] %in% interesting
## points(qq$x[show],qq$y[show],col="black",pch=22,cex=1.3)
## #text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col=c("red"),pos=2,offset=1.5,cex=1.3)  #,pch=25)

## long.genes<-c("TTN","DNAH10","MUC4","OBSCN","FAT3")
## interesting<-long.genes
## show<-symbols[qq$ord] %in% interesting
## points(qq$x[show],qq$y[show],col="black",pch=22,cex=1.3)
## #text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col=c("red"),pos=1,offset=1,cex=1.3)  #,pch=25)


interesting<-prev.study.genes
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="cyan",pch=20,cex=1.0)
#text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col=c("red"),pos=1,offset=1,cex=1.3)  #,pch=25)

prev.study.genes[prev.study.genes %in% tight.weight.sig.genes]
 ## [1] "KNSTRN" "CDKN2A" "COL1A2" "COL4A1" "COL4A2" "CSMD3"  "DCLK1"  "FAT1"  
 ## [9] "FLNA"   "NOTCH1" "NOTCH2" "TP53"

prev.study.genes[prev.study.genes %in% change]

interesting<-prev.study.genes[prev.study.genes %in% tight.weight.sig.genes]

"KNSTRN" "CDKN2A" "COL1A2" "COL4A1" "COL4A2" "CSMD3"  "DCLK1"  "FAT1"  
 [9] "FLNA"   "NOTCH1" "NOTCH2" "TP53"

interesting<-prev.study.genes[prev.study.genes %in% change]
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="black",pch=22,cex=1.5)
text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col=c("black"),pos=2,offset=1,cex=1.3)


interesting<-prev.study.genes[prev.study.genes %in% tight.weight.sig.genes]

interesting<-c("COL1A2","COL4A1","COL4A2","DCLK1","FLNA","NOTCH1","NOTCH2","TP53")


show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="black",pch=22,cex=2.0)
text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col=c("forestgreen"),pos=2,offset=1,cex=1.3) 

interesting<-c("KNSTRN","CDKN2A","CSMD3","FAT1")
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="black",pch=22,cex=2.0)
text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col=c("forestgreen"),pos=4,offset=1,cex=1.3) 


#selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="magenta",cex=1.3,offset=1,atpen='TRUE')

abline(h=qchisq(1.6e-6,df=1,ncp=0,lower.tail=FALSE,log.p=FALSE),col="green",lw=2)

leg.txt<-c(expression(paste("Genes")),"95% confidence intervals","Replicated Genes by Lee et.al.")
legend(0,80,leg.txt,col=c("blue","red","cyan"),lty=c(-1,2,-1),pch=c(1,-1,19),cex=2.0)

setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Figures")
fig.prefix<-"Q-Q_Final_lee_replication_genes"


savePlot(filename=paste(fig.prefix,".png",sep=''),type="png")
savePlot(filename=paste(fig.prefix,".jpeg",sep=''),type="jpeg")
savePlot(filename=paste(fig.prefix,".tiff",sep=''),type="tiff")
savePlot(filename=paste(fig.prefix,".bmp",sep=''),type="bmp")
dev.print(svg,paste(fig.prefix,".svg",sep=''))







################################## with blinded gene names

file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Burden.ALL.coding.somatic.with.Indels.noBenign.wRepeats.withGQ.filter.QCed_weights_TIGHT.txt"


meta.results.burden<-read.delim(file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)




## posns<-match(meta.results.burden[,"gene"],rownames(poss.model))
## meta.results.burden<-cbind(poss.model[posns,],meta.results.burden)
## meta.results.burden[1:5,1:10]
## write.table(meta.results.burden,file="containamination loci.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)




clusters.wanted<-c("SCC_Genes"," BRCA.Genes"," RAD51.Paralogues"," BLM.Complex"," Checkpoint.Proteins"," citric"," lipid_raft"," caveolae"," Checkpoint_extendedx1.noP53"," Checkpoint_extendedx2.noP53"," Citric.no.IDH"," BLM.Complex_AND_Checkpoint"," FANCD2_minimal_mono_ubi"," MLH_cluster")

dups<-duplicated(meta.results.burden[,"gene"])
meta.results.burden[grep("RYR",meta.results.burden[,"gene"]),]




sum(dups)

file
meta.results.burden[1:5,]
       dim(meta.results.burden)
## bad.genes<-c(bad.genes,clusters.wanted)
subset<-meta.results.burden[ !(meta.results.burden[,"gene"] %in%  clusters.wanted) ,]
       subset[1:10,]
subset<-subset[!is.na(as.numeric(subset[ ,"p"])),]

## subset[1:10,]
## no.weight.sig.genes<-subset[as.numeric(subset[ ,"p"])<1.2e-6,"gene"]
## tight.weight.sig.genes<-subset[as.numeric(subset[ ,"p"])<1.2e-6,"gene"]

symbols<-subset[,"gene"]

## subset[1:10,]
## posns<-match(interesting,subset[,"gene"])
## names(posns)<-interesting
## posns
## write.table(subset[1:350,],file="string.txt",,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## dups<-duplicated(subset[,"gene"])
## sum(dups)

## write.table(subset[1:350,],file="string.txt",,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
z1<-qchisq(subset[ ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)
z<-qchisq(subset[ ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) ## if have no chisq valuse
#z0<-qchisq(meta.results.skatO[ !(meta.results.skatO[,"gene"] %in% genes.and.clusters ) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) 
z[1:5]
subset[1:10,]

p.value<-as.numeric(subset[ ,"p"])
par(mfrow=c(1,1),font=2,font.lab=2,font.axis=2,mgp=c(3.5,1,0),mar=c(5,5,4,2)+0.1)


##  z<-rchisq(length(p.val), df=1, ncp = 0) ## null test


median(z,na.rm=TRUE)/0.456  #1.071491
#median(z0,na.rm=TRUE)/0.456  #1.071491


sum(is.na(p.val))
################## p-values
z0=qnorm(p.value[!is.na(p.value)]/2)
lambda = round(median(z0^2)/0.454,3)
lambda


plot(x=c(1:10))
## source("http://bioconductor.org/biocLite.R") 
## biocLite("GWASTools")
## setRep
## qq.Plot(pvals)


help[1:5,]
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
signif(range(z),digits=3)

z[!is.na(z)]
the.plot<-my.qq.plot(z,dist="chisq",df=1,ncp=0,col="blue",ylab="Observed chi-squared value",xlab="Expected chi-squared value",main="",cex=1,xlim=c(5,18),ylim=c(20,100),cex.lab=2.0,cex.axis=2.0,font.lab=2,font.axis=2,lwd=2,line="robust",plot.it=TRUE) # function defined below


## z.all<-qchisq(meta.results.burden[ !(meta.results.burden[,"gene"] %in% clusters.wanted) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)
## range(z.all)
qq<-  qq.data(z,plot.it=FALSE)       ## qq plot used same method as in car library
#points(qq$x,qq$y,col="magenta",pch=21)

symbols<-subset[,"gene"]


interesting<-c("TP53","NOTCH1","NOTCH2","CDKN2A")
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="red",pch=25,cex=2,lwd=2)
text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="forestgreen",pos=2,offset=1,cex=1.3)

interesting<-c("KNSTRN")
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="red",pch=25,cex=2,lwd=2)
text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="magenta",pos=4,offset=1,cex=1.3)


interesting<-c("STK19","DCLK1")
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col=c("magenta","forestgreen"),pch=20,cex=2)
#points(qq$x[show],qq$y[show],col=c("cyan"),pch=20,cex=2)
text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col=c("magenta","forestgreen"),pos=2,offset=1,cex=1.3)  #,pch=25)


## selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="forestgreen",cex=1.3,offset=1,atpen='TRUE')
## selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="magenta",cex=1.3,offset=1,atpen='TRUE')





abline(h=qchisq(1.6e-6,df=1,ncp=0,lower.tail=FALSE,log.p=FALSE),col="green",lw=2)

leg.txt<-c(expression(paste("Genes")),"95% confidence intervals","Genes with multiple mutations","Genes with focal mutations")
legend(6,80,leg.txt,col=c("blue","red","magenta","forestgreen"),lty=c(-1,2,-1,-1),pch=c(1,-1,19,19),cex=2.0)

setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Figures")
fig.prefix<-"Q-Q_Final_6sd_tight_weights_NO_NAMES"


savePlot(filename=paste(fig.prefix,".png",sep=''),type="png")
savePlot(filename=paste(fig.prefix,".jpeg",sep=''),type="jpeg")
savePlot(filename=paste(fig.prefix,".tiff",sep=''),type="tiff")
savePlot(filename=paste(fig.prefix,".bmp",sep=''),type="bmp")
dev.print(svg,paste(fig.prefix,".svg",sep=''))


##############################VIRAL GENES


## the.plot<-my.qq.plot(z,dist="chisq",df=1,ncp=0,col="blue",ylab="Observed chi-squared value",xlab="Expected chi-squared value",main="",cex=1,xlim=c(5,18),ylim=c(20,100),cex.lab=2.0,cex.axis=2.0,font.lab=2,font.axis=2,lwd=2,line="robust",plot.it=TRUE) # function defined below

the.plot<-my.qq.plot(z,dist="chisq",df=1,ncp=0,col="blue",ylab="Observed chi-squared value",xlab="Expected chi-squared value",main="",xlim=c(-1,18),ylim=c(-10,100),cex=1,cex.lab=2.0,cex.axis=2.0,font.lab=2,font.axis=2,lwd=2,line="robust",plot.it=TRUE) # function defined below

## z.all<-qchisq(meta.results.burden[ !(meta.results.burden[,"gene"] %in% clusters.wanted) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)
## range(z.all)
qq<-  qq.data(z,plot.it=FALSE)       ## qq plot used same method as in car library
#points(qq$x,qq$y,col="magenta",pch=21)

symbols<-subset[,"gene"]

viral<-c("CDK4","CDK6","CDKN2A","CCND1","MDM2","KAT2B","EP300","RBPJ","CHD4","CDKN1B","RB1","ATM","CREB1","PRKCG","ATM","CDKN1A")
hepB<-c("PRKCG","ATM","CDKN1A")

interesting<-viral
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="red",pch=25,cex=2,lwd=2)


interesting<-c("CHD4","CDKN2A")
show<-symbols[qq$ord] %in% interesting
#points(qq$x[show],qq$y[show],col="red",pch=25,cex=2,lwd=2)
text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="magenta",pos=2,offset=1,cex=1.3)

interesting<-c("PRKCG")
show<-symbols[qq$ord] %in% interesting
#points(qq$x[show],qq$y[show],col="red",pch=25,cex=2,lwd=2)
text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="magenta",pos=4,offset=1,cex=1.3)

interesting<-c("EBNA1BP2")
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="red",pch=20,cex=2,lwd=2)
text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="forestgreen",pos=2,offset=1.5,cex=1.3)



## selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="forestgreen",cex=1.3,offset=1,atpen='TRUE')
## selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="magenta",cex=1.3,offset=1,atpen='TRUE')
## interesting<-c("TP53","NOTCH1","NOTCH2","CDKN2A","KNSTRN")
## show<-symbols[qq$ord] %in% interesting
## points(qq$x[show],qq$y[show],col="red",pch=25,cex=2,lwd=2)
## text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="forestgreen",pos=2,offset=1,cex=1.3)



interesting<-c("STK19","DCLK1")




abline(h=qchisq(1.6e-6,df=1,ncp=0,lower.tail=FALSE,log.p=FALSE),col="green",lw=2)

leg.txt<-c(expression(paste("Genes")),"95% confidence intervals","KEGG:viral carinogenesis genes","EBV associated")
legend(0,80,leg.txt,col=c("blue","red","magenta","forestgreen"),lty=c(-1,2,-1,-1),pch=c(1,-1,19,19),cex=2.0)

setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Figures")
fig.prefix<-"Q-Q_Final_6sd_tight_weights_VIRAL"


savePlot(filename=paste(fig.prefi  x,".png",sep=''),type="png")
savePlot(filename=paste(fig.prefix,".jpeg",sep=''),type="jpeg")
savePlot(filename=paste(fig.prefix,".tiff",sep=''),type="tiff")
savePlot(filename=paste(fig.prefix,".bmp",sep=''),type="bmp")
dev.print(svg,paste(fig.prefix,".svg",sep=''))



#################################### Single point

file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Burden.ALL.coding.somatic.with.Indels.noBenign.wRepeats.withGQ.filter.QCed_weights_TIGHT.txt"
meta.results.burden.ori<-read.delim(file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Final_results/Burden.ALL.Single Point.withGQ.filter.QCed_weights_TIGHT.txt"
meta.results.burden<-read.delim(file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

## posns<-match(meta.results.burden[,"gene"],rownames(poss.model))
## meta.results.burden<-cbind(poss.model[posns,],meta.results.burden)
## meta.results.burden[1:5,1:10]
## write.table(meta.results.burden,file="containamination loci.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)




clusters.wanted<-c("SCC_Genes"," BRCA.Genes"," RAD51.Paralogues"," BLM.Complex"," Checkpoint.Proteins"," citric"," lipid_raft"," caveolae"," Checkpoint_extendedx1.noP53"," Checkpoint_extendedx2.noP53"," Citric.no.IDH"," BLM.Complex_AND_Checkpoint"," FANCD2_minimal_mono_ubi"," MLH_cluster")

bad<-c( "ATN1","NUDT18","OR5D13","MUC4")


meta.results.burden[1:5,]
dups<-duplicated(meta.results.burden[,"Name"])
sum(dups)


meta.results.burden[dups,][1:5,]
meta.results.burden[meta.results.burden[,"Name"]==meta.results.burden[dups,"Name"][1],]
meta.results.burden<-meta.results.burden[!dups,]

       dim(meta.results.burden)
## bad.genes<-c(bad.genes,clusters.wanted)
subset<-meta.results.burden[ !(meta.results.burden[,"gene"] %in%  clusters.wanted) & !(meta.results.burden[,"gene"] %in%  bad)  ,]
 subset[1:10,]
subset<-subset[!is.na(as.numeric(subset[ ,"p"])),]

## subset[1:10,]
## no.weight.sig.genes<-subset[as.numeric(subset[ ,"p"])<1.2e-6,"gene"]
## tight.weight.sig.genes<-subset[as.numeric(subset[ ,"p"])<1.2e-6,"gene"]
subset[1:10,]
symbols<-subset[,"gene"]

## subset[1:10,]
## posns<-match(interesting,subset[,"gene"])
## names(posns)<-interesting
## posns
## write.table(subset[1:350,],file="string.txt",,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## dups<-duplicated(subset[,"gene"])
## sum(dups)

## write.table(subset[1:350,],file="string.txt",,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
z1<-qchisq(subset[ ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)
z<-qchisq(subset[ ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) ## if have no chisq valuse
#z0<-qchisq(meta.results.skatO[ !(meta.results.skatO[,"gene"] %in% genes.and.clusters ) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) 
z[1:5]
subset[1:10,]

p.value<-as.numeric(subset[ ,"p"])
par(mfrow=c(1,1),font=2,font.lab=2,font.axis=2,mgp=c(3.5,1,0),mar=c(5,5,4,2)+0.1)


##  z<-rchisq(length(p.val), df=1, ncp = 0) ## null test


median(z,na.rm=TRUE)/0.456  #1.071491
#median(z0,na.rm=TRUE)/0.456  #1.071491


sum(is.na(p.val))
################## p-values
z0=qnorm(p.value[!is.na(p.value)]/2)
lambda = round(median(z0^2)/0.454,3)
lambda


plot(x=c(1:10))
## source("http://bioconductor.org/biocLite.R") 
## biocLite("GWASTools")
## setRep
## qq.Plot(pvals)


help[1:5,]
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
signif(range(z),digits=3)

the.plot<-my.qq.plot(z,dist="chisq",df=1,ncp=0,col="blue",ylab="Observed chi-squared value",xlab="Expected chi-squared value",main="",cex=1,xlim=c(7,22),ylim=c(10,150),cex.lab=2.0,cex.axis=2.0,font.lab=2,font.axis=2,lwd=2,line="robust",plot.it=TRUE) # function defined below


## z.all<-qchisq(meta.results.burden[ !(meta.results.burden[,"gene"] %in% clusters.wanted) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)
## range(z.all)
qq<-  qq.data(z,plot.it=FALSE)       ## qq plot used same method as in car library
#points(qq$x,qq$y,col="magenta",pch=21)

symbols<-subset[,"gene"]

known<-meta.results.burden.ori[meta.results.burden.ori[,"p"]< 1.6e-6 & !is.na(meta.results.burden.ori[,"p"]),"gene"]
have<-meta.results.burden[meta.results.burden[,"p"]< 1.6e-5 & !is.na(meta.results.burden[,"p"]),"gene"]


interesting<-known # ,"COL19A1"
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="black",pch=20,cex=1)
#text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="magenta",pos=4,offset=1,cex=1.3)  #,pch=25)



interesting<-have[ !(have %in% known) & !(have %in% bad) ]
interesting

"INTS1"
"DEFB127"
"NFKBIE"
"SPEF2"
"TNXB"
"ZNF319" 
"HYDIN"
"XIRP2"
"TNXB"
"MMP16"
"TNXB"

 ## [1] "INTS1"   "INTS1"   "DEFB127" "NFKBIE"  "SPEF2"   "TNXB"    "ZNF319" 
 ## [8] "HYDIN"   "XIRP2"   "TNXB"    "MMP16"   "TNXB"

[1] "INTS1"    "INTS1"    "DEFB127"  "NFKBIE"   "SPEF2"    "TNXB"    
 [7] "ZNF319"   "HYDIN"    "XIRP2"    "TNXB"     "MMP16"    "TNXB"    
[13] "STIP1"    "BAD"      "B4GALNT1" "PHACTR3"  "XIRP2"    "FREM2"   
[19] "PNN"      "KCNH5"    "LHCGR"    "CNTNAP5"  "ZNF467"   "DCAF13"  

to.color<-unique(interesting)
names(to.color)<-rainbow(length(to.color))
cols<-rep("blue",times=length(symbols))
names(cols)<-symbols
cols<-names(to.color)[match(names(cols),to.color)]
cols[is.na(cols)]<-"blue"

show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col=cols[qq$ord][show],pch=20,cex=2)



#text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="forestgreen",pos=2,offset=1,cex=1.0)  #,pch=25)

## qq$ord[c(87908,87945,87970,87979,87983,87987,87990)]
## selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="forestgreen",cex=1.3)
## show<-selected.data$ind
## #points(qq$x[show],qq$y[show],col="forestgreen",pch=20,cex=2)
## text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="forestgreen",pos=4,offset=1.5,cex=1.2)  #,pch=25)

## ## interesting<-c("TP53","NOTCH1","NOTCH2","CDKN2A","KNSTRN")
## ## show<-symbols[qq$ord] %in% interesting
## ## points(qq$x[show],qq$y[show],col="red",pch=25,cex=2,lwd=2)
## #text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="forestgreen",pos=2,offset=1,cex=1.3)

## symbols[qq$ord[selected.data$ind]]
## ## interesting<-c("TRHDE","PAPPA")
## ## show<-symbols[qq$ord] %in% interesting
## ## points(qq$x[show],qq$y[show],col=c("forestgreen","magenta"),pch=20,cex=2)
## #points(qq$x[show],qq$y[show],col=c("cyan"),pch=20,cex=2)
## #text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col=c("forestgreen","magenta"),pos=4,offset=1,cex=1.3)  #,pch=25)
## selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="forestgreen",cex=1.3,offset=2,pos=4,atpen='TRUE')
## selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="forestgreen",cex=1.3,offset=2,pos=1,atpen='TRUE')
## selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="forestgreen",cex=1.3,pos=2,offset=1,atpen='TRUE')
## selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="magenta",cex=1.3,offset=1,atpen='TRUE')


abline(h=qchisq(1.6e-6,df=1,ncp=0,lower.tail=FALSE,log.p=FALSE),col="green",lw=2)

leg.txt<-c(expression(paste("Genes")),"95% confidence intervals","Previously identified using Gene tests",to.color)
legend(7,150, bty="n",leg.txt,col=c("blue","red","black",names(to.color)),lty=c(-1,2,-1,rep(-1,times=length(to.color))),pch=c(1,-1,20,rep(19,times=length(to.color))),cex=1.5)

setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Figures")
fig.prefix<-"Q-Q_Final_6sd_tight_weights_single_point"


savePlot(filename=paste(fig.prefix,".png",sep=''),type="png")
savePlot(filename=paste(fig.prefix,".jpeg",sep=''),type="jpeg")
savePlot(filename=paste(fig.prefix,".tiff",sep=''),type="tiff")
savePlot(filename=paste(fig.prefix,".bmp",sep=''),type="bmp")
dev.print(svg,paste(fig.prefix,".svg",sep=''))



# FREM2 -check why it has moved about
# TNXB has pseudo genes and may be dodgy not caused by rescue but in redo with 27
# UNC80 has moved about chr2:210704125:210704125:C:T:snp


#####annotate curve




####################################SLIDINg window

file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Burden.ALL.coding.somatic.with.Indels.noBenign.wRepeats.withGQ.filter.QCed_weights_TIGHT.txt"
meta.results.burden.ori<-read.delim(file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Final_results/Burden.ALL.sliding.window.withGQ.filter.QCed_weights_TIGHT.txt"
meta.results.burden<-read.delim(file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

## posns<-match(meta.results.burden[,"gene"],rownames(poss.model))
## meta.results.burden<-cbind(poss.model[posns,],meta.results.burden)
## meta.results.burden[1:5,1:10]
## write.table(meta.results.burden,file="containamination loci.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)




clusters.wanted<-c("SCC_Genes"," BRCA.Genes"," RAD51.Paralogues"," BLM.Complex"," Checkpoint.Proteins"," citric"," lipid_raft"," caveolae"," Checkpoint_extendedx1.noP53"," Checkpoint_extendedx2.noP53"," Citric.no.IDH"," BLM.Complex_AND_Checkpoint"," FANCD2_minimal_mono_ubi"," MLH_cluster")

bad<-c( "ATN1","NUDT18","OR5D13","MUC4","OR9G4","OR8H3","ANKRD30BL,MIR663B")
significant<-meta.results.burden.ori[1:200,"gene"]

meta.results.burden[1:5,]
dups<-duplicated(meta.results.burden[,"gene"])
sum(dups)


meta.results.burden[dups,][1:5,]
meta.results.burden[meta.results.burden[,"gene"]==meta.results.burden[dups,"gene"][1],]
#meta.results.burden<-meta.results.burden[!dups,]

       dim(meta.results.burden)
## bad.genes<-c(bad.genes,clusters.wanted)
subset<-meta.results.burden[ !(meta.results.burden[,"gene"] %in%  clusters.wanted) & !(meta.results.burden[,"gene"] %in%  bad) & !(meta.results.burden[,"gene"] %in%  significant) ,]
subset<-meta.results.burden[ !(meta.results.burden[,"the.gene"] %in%  clusters.wanted) & !(meta.results.burden[,"the.gene"] %in%  bad)  ,]
 subset[1:10,]
subset<-subset[!is.na(as.numeric(subset[ ,"p"])),]

## subset[1:10,]
## no.weight.sig.genes<-subset[as.numeric(subset[ ,"p"])<1.2e-6,"gene"]
## tight.weight.sig.genes<-subset[as.numeric(subset[ ,"p"])<1.2e-6,"gene"]
subset[1:10,]
symbols<-subset[,"gene"]

## subset[1:10,]
## posns<-match(interesting,subset[,"gene"])
## names(posns)<-interesting
## posns
## write.table(subset[1:350,],file="string.txt",,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## dups<-duplicated(subset[,"gene"])
## sum(dups)

## write.table(subset[1:350,],file="string.txt",,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
z1<-qchisq(subset[ ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)
z<-qchisq(subset[ ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) ## if have no chisq valuse
#z0<-qchisq(meta.results.skatO[ !(meta.results.skatO[,"gene"] %in% genes.and.clusters ) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) 
z[1:5]
subset[1:10,]

p.value<-as.numeric(subset[ ,"p"])
par(mfrow=c(1,1),font=2,font.lab=2,font.axis=2,mgp=c(3.5,1,0),mar=c(5,5,4,2)+0.1)


##  z<-rchisq(length(p.val), df=1, ncp = 0) ## null test


median(z,na.rm=TRUE)/0.456  #1.071491
#median(z0,na.rm=TRUE)/0.456  #1.071491


sum(is.na(p.val))
################## p-values
z0=qnorm(p.value[!is.na(p.value)]/2)
lambda = round(median(z0^2)/0.454,3)
lambda


plot(x=c(1:10))
## source("http://bioconductor.org/biocLite.R") 
## biocLite("GWASTools")
## setRep
## qq.Plot(pvals)


help[1:5,]
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
signif(range(z),digits=3)

the.plot<-my.qq.plot(z,dist="chisq",df=1,ncp=0,col="blue",ylab="Observed chi-squared value",xlab="Expected chi-squared value",main="",cex=1,xlim=c(10,18),ylim=c(10,155),cex.lab=2.0,cex.axis=2.0,font.lab=2,font.axis=2,lwd=2,line="robust",plot.it=TRUE) # function defined below


## z.all<-qchisq(meta.results.burden[ !(meta.results.burden[,"gene"] %in% clusters.wanted) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)
## range(z.all)
qq<-  qq.data(z,plot.it=FALSE)       ## qq plot used same method as in car library
#points(qq$x,qq$y,col="magenta",pch=21)

symbols<-subset[,"the.gene"]

known<-meta.results.burden.ori[meta.results.burden.ori[,"p"]< 1.6e-6 & !is.na(meta.results.burden.ori[,"p"]),"gene"]
have<-meta.results.burden[meta.results.burden[,"p"]< 1.6e-6 & !is.na(meta.results.burden[,"p"]),"the.gene"]


interesting<-known # ,"COL19A1"
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="black",pch=20,cex=1)
#text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="magenta",pos=4,offset=1,cex=1.3)  #,pch=25)



interesting<-have[ !(have %in% known) & !(have %in% bad) ]
unique(interesting)

"INTS1"
"NFKBIE"
"ANKRD30BL"
"PRB3"
"TNXB"     
"BAD"
"CDH11"
"CHST4"
"XIRP2"



## INTS1
## DEFB127
## NFKBIE
## SPEF2
## ANKRD30BL
## PRB3
## TNXB     
## BAD
## CDH11
## CHST4
## XIRP2
## ZNF319 
## INTS1
## HYDIN
## MMP16


## [13] "PRB3"      "DOCK2"

##  [1] "INTS1"     "NFKBIE"    "ANKRD30BL" "ANKRD30BL" "PRB3"      "TNXB"     
##  [7] "BAD"       "CDH11"     "CHST4"     "OR9G4"     "OR8H3"     "XIRP2"    
## [13] "PRB3"      "DOCK2"

to.color<-unique(interesting)
names(to.color)<-rainbow(length(to.color))
cols<-rep("blue",times=length(symbols))
names(cols)<-symbols
cols<-names(to.color)[match(names(cols),to.color)]
cols[is.na(cols)]<-"blue"

show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col=cols[qq$ord][show],pch=20,cex=2)



#text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="forestgreen",pos=2,offset=1,cex=1.0)  #,pch=25)

## qq$ord[c(87908,87945,87970,87979,87983,87987,87990)]
## selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="forestgreen",cex=1.3)
## show<-selected.data$ind
## #points(qq$x[show],qq$y[show],col="forestgreen",pch=20,cex=2)
## text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="forestgreen",pos=4,offset=1.5,cex=1.2)  #,pch=25)

## ## interesting<-c("TP53","NOTCH1","NOTCH2","CDKN2A","KNSTRN")
## ## show<-symbols[qq$ord] %in% interesting
## ## points(qq$x[show],qq$y[show],col="red",pch=25,cex=2,lwd=2)
## #text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="forestgreen",pos=2,offset=1,cex=1.3)

## symbols[qq$ord[selected.data$ind]]
## ## interesting<-c("TRHDE","PAPPA")
## ## show<-symbols[qq$ord] %in% interesting
## ## points(qq$x[show],qq$y[show],col=c("forestgreen","magenta"),pch=20,cex=2)
## #points(qq$x[show],qq$y[show],col=c("cyan"),pch=20,cex=2)
## #text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col=c("forestgreen","magenta"),pos=4,offset=1,cex=1.3)  #,pch=25)
#selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="forestgreen",cex=1.3,offset=2,pos=4,atpen='TRUE')
## selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="forestgreen",cex=1.3,offset=2,pos=1,atpen='TRUE')
## selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="forestgreen",cex=1.3,pos=2,offset=1,atpen='TRUE')
## selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="magenta",cex=1.3,offset=1,atpen='TRUE')


abline(h=qchisq(1.6e-6,df=1,ncp=0,lower.tail=FALSE,log.p=FALSE),col="green",lw=2)

leg.txt<-c(expression(paste("Genes")),"95% confidence intervals","Previously identified using Gene tests",to.color)
legend(10,150, bty="n",leg.txt,col=c("blue","red","black",names(to.color)),lty=c(-1,2,-1,rep(-1,times=length(to.color))),pch=c(1,-1,20,rep(19,times=length(to.color))),cex=1.5)

setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Figures")
fig.prefix<-"Q-Q_Final_6sd_tight_weights_sliding_window"


savePlot(filename=paste(fig.prefix,".png",sep=''),type="png")
savePlot(filename=paste(fig.prefix,".jpeg",sep=''),type="jpeg")
savePlot(filename=paste(fig.prefix,".tiff",sep=''),type="tiff")
savePlot(filename=paste(fig.prefix,".bmp",sep=''),type="bmp")
dev.print(svg,paste(fig.prefix,".svg",sep=''))



# FREM2 -check why it has moved about
# TNXB has pseudo genes and may be dodgy not caused by rescue but in redo with 27
# UNC80 has moved about chr2:210704125:210704125:C:T:snp


#####annotate curve





















############################################# BARPLOTS
############################################# BARPLOTS
############################################# BARPLOTS
############################################# BARPLOTS
############################################# BARPLOTS
############################################# BARPLOTS
############################################# BARPLOTS


file <- "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Figures/leopharma_mutation_matrix_with_counts.csv"
counts<-read.delim(file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

freq<-counts[,c("GENE","tumour_percent","tumour_count","pd_percent","pd_count")]
freq[1:5,]
counts[1:5,]
counts$GENE
###### supplement
slide.hits<-c( "INTS1","NFKBIE","PRB3","TNXB","BAD","CHST4","XIRP2","DOCK2")

file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/TOP_550_LENGTH_CONTROL.GENOTYPE.conponents..Burden.clusters.coding.somatic.with.Indels.noBenign.ORI.withGQ.filter.QCed_weights_TIGHT.txt"
file <- "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Final_results/TOP_550_LENGTH_CONTROL.GENOTYPE.conponents..Burden.clusters.sliding.window.withGQ.filter.QCed_weights_TIGHT.txt"
components<-read.delim(file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

i<-3
for( i in 1:length(slide.hits)){
  do<-slide.hits[i]
  cancers<-components[components[,"Gene.Names"]==do,"muts.in.cases"]
  cancers<-unique(unlist(strsplit(cancers,split=",")))
  cancers<-cancers[cancers!=""]

   pd<-components[components[,"Gene.Names"]==do,"muts.in.PD"]
  pd<-unique(unlist(strsplit(pd,split=",")))
  pd<-pd[pd!=""]
  extra<-c(paste(do,"*",sep=""),100*length(cancers)/70,length(cancers),100*length(pd)/25,length(pd))
 freq<-rbind(freq,extra)
  
}
  
freq

target<-"tumour_percent"



freq<-freq[order(freq[,target],decreasing=TRUE),]
col<-rep("black",times=dim(freq)[1])
labels<-rep("",times=dim(freq)[1])
names(col)<-freq[,"GENE"]
names(labels)<-freq[,"GENE"]
col[names(col) %in% c("TP53","NOTCH1","NOTCH2","KNSTRN","PRKCG") ]<-"magenta"
col[names(col) %in% c("STK19","KNSTRN","BCL2L12","CCDC61","EBNA1BP2","PHACTR3","CDKN2A","UNC80","C12orf42","TRHDE","DCLK1",paste(slide.hits,"*",sep=""))]<-"forestgreen"

labels[names(labels) %in% c("TP53","NOTCH1","NOTCH2","KNSTRN","PRKCG","CDKN2A","EBNA1BP2","STK19")]<-names(labels)[names(labels) %in% c("TP53","NOTCH1","NOTCH2","KNSTRN","PRKCG","CDKN2A","EBNA1BP2","STK19")]

col[ c("TP53","NOTCH1","NOTCH2","KNSTRN","STK19","CDKN2A")]    
bp<-barplot(as.numeric(freq[,target]),col=col,main=("Percent SCC or AK samples"),axisnames=TRUE,font.main = 4,lwd=2,cex.axis=2 )
#text(bp, par("usr")[3] - 0.025, srt = 45, adj = 1, labels = as.character(labels), xpd = TRUE, font = 2,cex=1.5)

text(bp, par("usr")[3] - 0.025, srt = 45, adj = 1, labels = as.character(names(labels)), xpd = TRUE, font = 2.5)

#fig.prefix<-"barplot_SCCAK_noNames"
fig.prefix<-"barplot_SCCAK_Names"

savePlot(filename=paste(fig.prefix,".png",sep=''),type="png")
savePlot(filename=paste(fig.prefix,".jpeg",sep=''),type="jpeg")
savePlot(filename=paste(fig.prefix,".tiff",sep=''),type="tiff")
#savePlot(filename=paste(fig.prefix,".bmp",sep=''),type="bmp")
dev.print(svg,paste(fig.prefix,".svg",sep=''))
############################################# BARPLOTS
############################################# BARPLOTS
############################################# BARPLOTS

target<-"pd_percent"

bp<-barplot(as.numeric(freq[,target]),col=col,main=("Percent PD samples"),axisnames=TRUE,font.main = 4,lwd=2,cex.axis=2 )
#text(bp, par("usr")[3] - 0.025, srt = 45, adj = 1, labels = as.character(labels), xpd = TRUE, font = 2,cex=1.5)

text(bp, par("usr")[3] - 0.025, srt = 45, adj = 1, labels = as.character(names(labels)), xpd = TRUE, font = 2)

fig.prefix<-"barplot_PD_NoNames"
fig.prefix<-"barplot_PD_Names"

savePlot(filename=paste(fig.prefix,".png",sep=''),type="png")
savePlot(filename=paste(fig.prefix,".jpeg",sep=''),type="jpeg")
savePlot(filename=paste(fig.prefix,".tiff",sep=''),type="tiff")
#savePlot(filename=paste(fig.prefix,".bmp",sep=''),type="bmp")
dev.print(svg,paste(fig.prefix,".svg",sep=''))






 require(grDevices) # for colours
     tN <- table(Ni <- stats::rpois(100, lambda = 5))
     r <- barplot(tN, col = rainbow(20))
     #- type = "h" plotting *is* 'bar'plot
     lines(r, tN, type = "h", col = "red", lwd = 2)
     
     barplot(tN, space = 1.5, axisnames = FALSE,
             sub = "barplot(..., space= 1.5, axisnames = FALSE)")
     
     barplot(VADeaths, plot = FALSE)
     barplot(VADeaths, plot = FALSE, beside = TRUE)
     
     mp <- barplot(VADeaths) # default
     tot <- colMeans(VADeaths)
     text(mp, tot + 3, format(tot), xpd = TRUE, col = "blue")
     barplot(VADeaths, beside = TRUE,
             col = c("lightblue", "mistyrose", "lightcyan",
                     "lavender", "cornsilk"),
             legend = rownames(VADeaths), ylim = c(0, 100))
     title(main = "Death Rates in Virginia", font.main = 4)
     
     hh <- t(VADeaths)[, 5:1]
     mybarcol <- "gray20"
     mp <- barplot(hh, beside = TRUE,
             col = c("lightblue", "mistyrose",
                     "lightcyan", "lavender"),
             legend = colnames(VADeaths), ylim = c(0,100),
             main = "Death Rates in Virginia", font.main = 4,
             sub = "Faked upper 2*sigma error bars", col.sub = mybarcol,
             cex.names = 1.5)
     segments(mp, hh, mp, hh + 2*sqrt(1000*hh/100), col = mybarcol, lwd = 1.5)
     stopifnot(dim(mp) == dim(hh))  # corresponding matrices
     mtext(side = 1, at = colMeans(mp), line = -2,
           text = paste("Mean", formatC(colMeans(hh))), col = "red")
     
     # Bar shading example
     barplot(VADeaths, angle = 15+10*1:5, density = 20, col = "black",
             legend = rownames(VADeaths))
     title(main = list("Death Rates in Virginia", font = 4))
     
     # border :
     barplot(VADeaths, border = "dark blue") 
     
     
     # log scales (not much sense here):
     barplot(tN, col = heat.colors(12), log = "y")
     barplot(tN, col = gray.colors(20), log = "xy")
     
     # args.legend
     barplot(height = cbind(x = c(465, 91) / 465 * 100,
                            y = c(840, 200) / 840 * 100,
                            z = c(37, 17) / 37 * 100),
             beside = FALSE,
             width = c(465, 840, 37),
             col = c(1, 2),
             legend.text = c("A", "B"),
             args.legend = list(x = "topleft"))
     


###Tumour
library(ggplot2)
data<-counts
data.tumour<-data[,c("GENE","tumour_percent")]
data.tumour<-data.tumour[order(-data.tumour$tumour_percent),] 
genes<-as.character(data.tumour$GENE)
colours<-rep("black",length(genes))
interesting.green<-c("TRHDE","STK19","KNSTRN","BCL2L12","CCDC61","EBNA1BP2")
interesting.magenta<-c("TP53","NOTCH1","NOTCH2","FAT1","DGKI","COL19A1","PAPPA") 
interesting.red<-c("CHR19_CNV")
colours
posns<-match(interesting.green,genes)
colours[posns]<-"forestgreen"
posns<-match(interesting.magenta,genes)
colours[posns]<-"magenta"
posns<-match(interesting.red,genes)
colours[posns]<-"red"
data.tumour <- transform(data.tumour, GENE = factor(GENE, levels = data.tumour$GENE))
bp<-ggplot(data=data.tumour, aes(x=GENE, y=tumour_percent, fill=GENE)) + geom_bar(stat="identity")  + xlab("Gene") + ylab("Percent SCC or AK samples") 
bp<-bp + theme(legend.position="none",axis.text.x = element_text(angle = 90, hjust = 1))+ scale_fill_manual(values=colours)
bp<-bp + theme(text = element_text(size=10),axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))
plot(bp)
ggsave(filename = "scc_ak_counts_by_gene.png",plot=bp,width=12.9,height=7.42)

###PD
data.pd<-data[,c("GENE","pd_percent")]
#data.pd<-data.pd[order(-data.tumour$pd_count),] 
data.pd <- transform(data.pd, GENE = factor(GENE, levels = data.tumour$GENE))
png("pd_counts_by_gene.png")
bp<-ggplot(data=data.pd, aes(x=GENE, y=pd_percent, fill=GENE)) + geom_bar(stat="identity")  + xlab("Gene") + ylab("Percent PD samples") 
bp<-bp + theme(legend.position="none",axis.text.x = element_text(angle = 90, hjust = 1))+ scale_fill_manual(values=colours)
bp<-bp + theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))
ggsave(filename = "pd_counts_by_gene.png",plot=bp,width=12.9,height=7.42)

data.pd<-data.pd[order(-data.pd$pd_percent),] 
genes<-as.character(data.pd$GENE)
colours<-rep("black",length(genes))
posns<-match(interesting.green,genes)
colours[posns]<-"forestgreen"
posns<-match(interesting.magenta,genes)
colours[posns]<-"magenta"
posns<-match(interesting.red,genes)
colours[posns]<-"red"
data.pd <- transform(data.pd, GENE = factor(GENE, levels = data.pd$GENE))
bp<-ggplot(data=data.pd, aes(x=GENE, y=pd_percent, fill=GENE)) + geom_bar(stat="identity")  + xlab("Gene") + ylab("Percent PD samples") 
bp<-bp + theme(legend.position="none",axis.text.x = element_text(angle = 90, hjust = 1))+ scale_fill_manual(values=colours)
bp<-bp + theme( axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))
ggsave(filename = "pd_counts_by_gene2.png",plot=bp,width=12.9,height=7.42)












o <- -(log10(p.value))
e <- -log10( 1:length(o)/length(o) )
       sum(as.numeric(is.na(o)))

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






      qq.data.mean<- function (x, plot.it = TRUE, distribution = "norm", df=1, the.mean=0, the.sd=1,  xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)) , ...)
{
    good <- !is.na(x)
    ord <- order(x[good])
    ord.x <- x[good][ord]
    q.function <- eval(parse(text = paste("q", distribution, 
        sep = "")))
    n <- length(ord.x)
    P <- ppoints(n)
    z <- q.function(P, mean=the.mean, sd=the.sd, ...)

    if (plot.it)
        plot(z, ord.x, xlab = xlab, ylab = ylab, ...)
    invisible(list(x = z, y = ord.x, ord=ord))
}











############################## replication selection Pac-Bio pac-bio



data<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Gene_expression/imm_cancer_compare.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
data[1:5,]

sum(data[,1] %in% data[,2])
data[data[,1] %in% data[,2],]

sd.thresh<-6

p<-0.05
n=230
alt.counts.thresh=0

z<-(alt.counts.thresh- n*p) / sqrt(n*p*(1-p))
z
2*(pnorm(abs(z), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
 # 2*(pnorm(abs(c(1:8)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)) #2*(pnorm(abs(seq(from=3,to=5,by=0.1)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
# 2*(pnorm(abs(4.4), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)) ##4 sd greater than mean

pass.possonian.control.model<- poss.model[,"P-value"]> 1e-5  | is.na( poss.model[,"P-value"])


colnames(summary.geno.extra)
######################################
n<-max(as.integer(summary.geno.extra[,"TOTAL.Alleles.Control"]))
n
alt.counts.thresh<-1
while( (alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) <= sd.thresh){alt.counts.thresh<-alt.counts.thresh+1}
alt.counts.thresh






file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Final_results/TOP_550_LENGTH_CONTROL.GENOTYPE.conponents..Burden.clusters.coding.somatic.with.Indels.noBenign.ORI.NO.GENO.REVOVERY.csv"

file <- "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/SCC_recurrent_point_mutation_tights_replication.csv"

rep<-read.delim(file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)


colnames(rep)
rep[1:5,]
#rep.do<-rep[rep[,"REP.2"]!="",]
rep.do<-rep[rep[,"REP.1"]!="",]

## rep.do<-rep
## rep.do<-rep.do[rep.do[,"Gene.Names"] %in% c("EBNA1BP2","TP53"),] ## for EB samlps
dim(rep.do)
rep.do[,"Gene.Names"]
unique(rep.do[,"Gene.Names"])
    
rep.do.samples<-c(rep.do[,"muts.in.cases"],rep.do[,"muts.in.PD"])
rep.do.samples<-rep.do.samples[rep.do.samples!=""]
rep.do.samples<-unlist(apply(as.matrix(rep.do.samples),1,function(x){ strsplit(x,split=",")}) )
chk<-sort(table(rep.do.samples ),decreasing=TRUE)

chk

rep.do.samples



sort(table(rep.do.samples ),decreasing=TRUE)[1:48]


LPH-001-21_AK1 LPH-001-25_AK1  LPH-001-4_AK1 LPH-001-14_AK2 LPH-001-18_AK2 
             9              8              8              7              7 
  LPH-001-1_PD LPH-001-21_AK2  LPH-001-3_SCC  LPH-001-6_AK1 LPH-001-14_AK1 
             7              7              7              7              6 
LPH-001-15_AK2 LPH-001-16_AK1 LPH-001-18_AK1  LPH-001-1_AK1  LPH-001-1_SCC 
             6              6              6              6              6 
 LPH-001-25_PD  LPH-001-3_AK1  LPH-001-3_AK2  LPH-001-4_AK2  LPH-001-8_AK1 
             6              6              6              6              6 
 LPH-001-12_PD LPH-001-24_AK1 LPH-001-24_SCC LPH-001-25_SCC LPH-001-26_SCC 
             5              5              5              5              5 
  LPH-001-3_PD  LPH-001-4_SCC  LPH-001-6_AK2   LPH-001-6_PD   LPH-001-7_PD 
             5              5              5              5              5 
 LPH-001-9_AK2  LPH-001-10_PD LPH-001-12_AK2 LPH-001-12_SCC LPH-001-16_AK2 
             5              4              4              4              4 
 LPH-001-16_PD LPH-001-16_SCC LPH-001-17_AK1  LPH-001-1_AK2 LPH-001-21_SCC 
             4              4              4              4              4 
 LPH-001-28_PD   LPH-001-9_PD LPH-001-12_AK1 LPH-001-13_SCC LPH-001-17_SCC 
             4              4              3              3              3


##LPH-001-20_AK1  LPH-001-21_PD  LPH-001-22_PD #LPH-001-22_SCC LPH-001-24_AK2
##              3              3              3              3              3 
##  LPH-001-24_PD LPH-001-26_AK1 LPH-001-28_SCC   LPH-001-5_PD  LPH-001-5_SCC 
##              3              3              3              3              3 
##  LPH-001-6_SCC  LPH-001-8_AK2   LPH-001-8_PD  LPH-001-9_AK1 LPH-001-10_AK2 
##              3              3              3              3              2 
## LPH-001-13_AK2  LPH-001-14_PD  LPH-001-15_PD LPH-001-17_AK2 LPH-001-22_AK1 
##              2              2              2              2              2 
## LPH-001-22_AK2 LPH-001-26_AK2 LPH-001-30_SCC  LPH-001-5_AK1  LPH-001-5_AK2 
##              2              2              2              2              2 
##  LPH-001-13_PD  LPH-001-17_PD LPH-001-18_SCC LPH-001-20_AK2  LPH-001-27_PD 
##              1              1              1              1              1 
## LPH-001-28_AK1 LPH-001-28_AK2 LPH-001-30_AK2  LPH-001-30_PD   LPH-001-4_PD 
##              1              1              1              1              1 
##  LPH-001-7_AK1  LPH-001-7_SCC  LPH-001-8_SCC 
##              1              1              1 
test<-"LPH-001-27_SCC"
test<-"LPH-001-17_SCC" # "EBNA1BP2" "COL4A4"  
test<-"LPH-001-20_AK1" # "TRHDE"    "C12orf42"
test<-"LPH-001-17_AK1" # "EBNA1BP2" "BCL2L12" 
test<-"LPH-001-16_SCC" # "EBNA1BP2" "TP53"    
test<-"LPH-001-15_AK1" #"COL4A4" "DGKI"  
test<-"LPH-001-12_AK2" #  "BCL2L12" "DGKI" 
test<-"LPH-001-12_AK1" # "KNSTRN" "TP53"

test<-"LPH-001-24_PD" # "KNSTRN"  "BCL2L12"
test<-"LPH-001-28_PD" #"C12orf42" "TP53" 
test<-"LPH-001-9_PD" # "EBNA1BP2" "BCL2L12" 


chk[1:48]
#to.chk<-chk # EB test
to.chk<-chk[chk==4]
i<-1
collect<-{}
for (i in 1:length(to.chk)){
test<-names(to.chk)[i]
if(!grepl("_PD",test)){
test.in<-grepl(test,c(rep.do[,"muts.in.cases"]))
print(test)
print(rep.do[test.in,"Gene.Names"])
collect<-c(collect, paste(test,length(unique(rep.do[test.in,"Gene.Names"])),sep=":"))
print("---------------------------------------")
}else{
print("PD--------------")
test.in<-grepl(test,rep.do[,"muts.in.PD"])
print(test)
print(rep.do[test.in,"Gene.Names"])
collect<-c(collect, paste(test,length(unique(rep.do[test.in,"Gene.Names"])),sep=":"))
print("---------------------------------------")
}
 }

print(rep.do[test.in,"Gene.Names"])
    
collect

[1] "LPH-001-12_AK1"
[1] "KNSTRN" "TP53"   "INTS1" 
[1] "---------------------------------------"
[1] "LPH-001-13_SCC"
[1] "CCDC61"  "CCDC61"  "BCL2L12"
[1] "---------------------------------------"
[1] "LPH-001-17_SCC"
[1] "EBNA1BP2" "INTS1"    "INTS1"   
[1] "---------------------------------------"
[1] "LPH-001-20_AK1"
[1] "TRHDE"    "C12orf42" "INTS1"   
[1] "---------------------------------------"
[1] "PD--------------"
[1] "LPH-001-21_PD"
[1] "DCLK1" "TP53"  "STK19"
[1] "---------------------------------------"
[1] "PD--------------"
[1] "LPH-001-22_PD"
[1] "BCL2L12" "UNC80"   "INTS1"  
[1] "---------------------------------------"
[1] "LPH-001-22_SCC"
[1] "KNSTRN" "TP53"   "INTS1" 
[1] "---------------------------------------"
[1] "LPH-001-24_AK2"
[1] "C12orf42" "CCDC61"   "INTS1"   
[1] "---------------------------------------"
[1] "PD--------------"
[1] "LPH-001-24_PD"
[1] "KNSTRN"  "BCL2L12" "INTS1"  
[1] "---------------------------------------"
[1] "LPH-001-26_AK1"
[1] "TP53"   "CCDC61" "INTS1" 
[1] "---------------------------------------"
[1] "LPH-001-28_SCC"
[1] "CCDC61" "INTS1"  "INTS1" 
[1] "---------------------------------------"
[1] "PD--------------"
[1] "LPH-001-5_PD"
[1] "STK19" "INTS1" "INTS1"
[1] "---------------------------------------"
[1] "LPH-001-5_SCC"
[1] "BCL2L12" "INTS1"   "INTS1"  
[1] "---------------------------------------"
[1] "LPH-001-6_SCC"
[1] "KNSTRN" "INTS1"  "INTS1" 
[1] "---------------------------------------"
[1] "LPH-001-8_AK2"
[1] "CCDC61" "INTS1"  "INTS1" 
[1] "---------------------------------------"
[1] "PD--------------"
[1] "LPH-001-8_PD"
[1] "C12orf42" "TP53"     "CCDC61"  
[1] "---------------------------------------"
[1] "LPH-001-9_AK1"
[1] "EBNA1BP2" "INTS1"    "INTS1"   
[1] "---------------------------------------"
> n

gsub(":2$","",collect[grepl(":2",collect)])
EB and TP53 ORI
[1] "LPH-001-21_AK2" "LPH-001-17_SCC" "LPH-001-6_AK1"  "LPH-001-15_AK2"
[5] "LPH-001-21_SCC" "LPH-001-3_SCC"  "LPH-001-4_AK1"


both width recovered
 [1] "LPH-001-21_AK1" "LPH-001-3_AK1"  "LPH-001-3_PD"   "LPH-001-15_AK2"
 [5] "LPH-001-21_AK2" "LPH-001-24_SCC" "LPH-001-6_PD"   "LPH-001-8_AK1" 
 [9] "LPH-001-9_AK1"  "LPH-001-16_SCC" "LPH-001-17_AK1" "LPH-001-17_SCC"
[13] "LPH-001-21_SCC" "LPH-001-22_AK2" "LPH-001-4_SCC"  "LPH-001-9_PD"




gsub(":1$","",collect[grepl(":1",collect)])

TP53 only  ORI no EB only found
 [1] "LPH-001-6_SCC"  "LPH-001-17_AK2" "LPH-001-1_AK2"  "LPH-001-26_AK1"
 [5] "LPH-001-30_AK1" "LPH-001-8_AK1"  "LPH-001-10_SCC" "LPH-001-16_AK1"
 [9] "LPH-001-1_AK1"  "LPH-001-22_AK2" "LPH-001-24_AK2" "LPH-001-25_SCC"
[13] "LPH-001-26_SCC" "LPH-001-28_AK1" "LPH-001-28_SCC" "LPH-001-30_AK2"
[17] "LPH-001-5_SCC"  "LPH-001-9_AK1"  "LPH-001-10_AK2" "LPH-001-12_AK1"
[21] "LPH-001-12_AK2" "LPH-001-13_AK2" "LPH-001-14_AK2" "LPH-001-15_AK1"
[25] "LPH-001-15_PD"  "LPH-001-15_SCC" "LPH-001-16_AK2" "LPH-001-17_AK1"
[29] "LPH-001-18_AK2" "LPH-001-18_SCC" "LPH-001-1_PD"   "LPH-001-20_AK2"
[33] "LPH-001-26_AK2" "LPH-001-28_PD"  "LPH-001-3_AK1"  "LPH-001-4_SCC" 
[37] "LPH-001-5_AK1"  "LPH-001-5_AK2"  "LPH-001-6_PD"   "LPH-001-7_AK2"

# recoved EN but no TP53
"LPH-001-4_AK1"   "LPH-001-7_PD" "LPH-001-6_AK1" "LPH-001-4_AK2"



"LPH-001-15_AK1" -> "LPH-001-9_PD"
"LPH-001-12_AK1"-> "LPH-001-28_PD"

LPH-001-12_AK1
LPH-001-13_SCC
LPH-001-17_SCC


Want to add:
"LPH-001-21_PD"
"LPH-001-22_PD"
"LPH-001-24_PD"
"LPH-001-8_PD"


REaplace
"LPH-001-12_AK1"
"LPH-001-13_SCC"
"LPH-001-17_SCC" 
"LPH-001-12_AK2"


wanted<-names(sort(table(rep.do.samples ),decreasing=TRUE))[1:45]
wanted[wanted=="LPH-001-12_AK1"]<-"LPH-001-21_PD"
wanted[wanted=="LPH-001-13_SCC"]<-"LPH-001-22_PD"
wanted[wanted=="LPH-001-17_SCC"]<-"LPH-001-24_PD"
wanted[wanted=="LPH-001-12_AK2"]<-"LPH-001-8_PD"

wanted
select<-sort(table(rep.do.samples ),decreasing=TRUE)[wanted]
select


## LPH-001-21_AK1 LPH-001-21_AK2  LPH-001-4_AK1 LPH-001-14_AK2 LPH-001-18_AK1 LPH-001-18_AK2  LPH-001-3_SCC  LPH-001-6_AK1   LPH-001-6_PD LPH-001-14_AK1 LPH-001-15_AK2 LPH-001-16_AK1 LPH-001-25_AK1 
##              6              6              6              5              5              5              5              5              5              4              4              4              4 
##  LPH-001-3_AK1  LPH-001-3_AK2  LPH-001-4_AK2  LPH-001-4_SCC   LPH-001-7_PD  LPH-001-8_AK1   LPH-001-8_PD  LPH-001-12_PD LPH-001-12_SCC LPH-001-13_SCC LPH-001-16_AK2  LPH-001-1_AK1   LPH-001-1_PD 
##              4              4              4              4              4              4              4              3              3              3              3              3              3 
##  LPH-001-1_SCC  LPH-001-21_PD LPH-001-21_SCC LPH-001-24_AK1 LPH-001-24_SCC  LPH-001-25_PD LPH-001-26_SCC   LPH-001-3_PD  LPH-001-9_AK2  LPH-001-10_PD  LPH-001-28_PD LPH-001-12_AK2  LPH-001-14_PD 
##              3              3              3              3              3              3              3              3              3              2              2              2              2 
##   LPH-001-9_PD  LPH-001-16_PD LPH-001-16_SCC LPH-001-17_AK1 LPH-001-17_SCC LPH-001-20_AK1 
##              2              2              2              2              2              2 
> 

select.blood<-gsub("-","_",names(select))
select.blood<-  unlist(apply(as.matrix(select.blood),1,function(x){      unlist(strsplit(x,split="_"))[3] }))
 sort(table(select.blood ),decreasing=TRUE)

select.samples<-select
names(select.samples)<-select.blood
select.samples

sort(tapply(select.samples,names(select.samples),sum))

10 20 22 13 26  9 17  7 15 25  8  1 12 18 24  6 14 16  4  3 21 
 2  2  2  3  3  3  4  4  6  7  8  9 10 10 10 10 11 11 14 16 18


> sort(tapply(select.samples,names(select.samples),sum))
22 10 17 28 26  7 15 12  8  9 14 18 24  6 16 25  4  1 21  3 
 3  4  4  4  5  5  6  9  9  9 13 13 13 17 18 19 19 23 23 24 

run

21,3,1 bloods
select.use<-names(select)
length(select.use)

i<-1
j<-1
rep.do[,"muts.in.cases"]
gene.hits<-{}

 # for(j in 1:dim(rep.do)[1]){
    
    sum.test<-0
    for(i in 1:length(select.use)){
   test1<-grepl(select.use[i],rep.do[,"muts.in.cases"])
   test2<-grepl(select.use[i],rep.do[,"muts.in.PD"])
   sum.test<-test1 | test2 ## genes that have hits in that sample
 
   print(sum.test)
   found<-rep.do[sum.test,c("gene","REP.1","Gene.Names")]
   dups<-duplicated(found[,"REP.1"])
   found<-found[!dups,]
   print(found)
    if(i==1){gene.hits<-found}else{
      gene.hits<-rbind(gene.hits,found)
   
       }
 }
gene.hits

table(gene.hits[,"Gene.Names"])

 BCL2L12 C12orf42   CCDC61    DCLK1 EBNA1BP2    INTS1   KNSTRN    STK19 
      25       17       13        9       21       42       14       13 
    TP53    TRHDE    UNC80 
      19       12       12 


   rep.do[,c("Gene.Names","muts.in.cases","muts.in.PD")]
  
c(rep.do[,"muts.in.cases"],rep.do[,"muts.in.PD"])
select.use

 [1] "LPH-001-21_AK1" "LPH-001-21_AK2" "LPH-001-4_AK1"  "LPH-001-14_AK2" "LPH-001-18_AK1" "LPH-001-18_AK2" "LPH-001-3_SCC"  "LPH-001-6_AK1"  "LPH-001-6_PD"   "LPH-001-14_AK1" "LPH-001-15_AK2"
[12] "LPH-001-16_AK1" "LPH-001-25_AK1" "LPH-001-3_AK1"  "LPH-001-3_AK2"  "LPH-001-4_AK2"  "LPH-001-4_SCC"  "LPH-001-7_PD"   "LPH-001-8_AK1"  "LPH-001-8_PD"   "LPH-001-12_PD"  "LPH-001-12_SCC"
[23] "LPH-001-13_SCC" "LPH-001-16_AK2" "LPH-001-1_AK1"  "LPH-001-1_PD"   "LPH-001-1_SCC"  "LPH-001-21_PD"  "LPH-001-21_SCC" "LPH-001-24_AK1" "LPH-001-24_SCC" "LPH-001-25_PD"  "LPH-001-26_SCC"
[34] "LPH-001-3_PD"   "LPH-001-9_AK2"  "LPH-001-10_PD"  "LPH-001-28_PD"  "LPH-001-12_AK2" "LPH-001-14_PD"  "LPH-001-9_PD"   "LPH-001-16_PD"  "LPH-001-16_SCC" "LPH-001-17_AK1" "LPH-001-17_SCC"
[45] "LPH-001-20_AK1" "LPH-001-18_SCC" "LPH-001-13_SCC"

"LPH-001-21"
"LPH-001-3"
"LPH-001-1"

selected<-c("LPH-001-21_AK1","LPH-001-21_AK2","LPH-001-4_AK1","LPH-001-14_AK2","LPH-001-18_AK1","LPH-001-18_AK2","LPH-001-3_SCC","LPH-001-6_AK1","LPH-001-6_PD","LPH-001-14_AK1","LPH-001-15_AK2",
"LPH-001-16_AK1","LPH-001-25_AK1","LPH-001-3_AK1","LPH-001-3_AK2","LPH-001-4_AK2","LPH-001-4_SCC","LPH-001-7_PD","LPH-001-8_AK1","LPH-001-8_PD","LPH-001-12_PD","LPH-001-12_SCC",
"LPH-001-13_SCC","LPH-001-16_AK2","LPH-001-1_AK1","LPH-001-1_PD","LPH-001-1_SCC","LPH-001-21_PD","LPH-001-21_SCC","LPH-001-24_AK1","LPH-001-24_SCC","LPH-001-25_PD","LPH-001-26_SCC",
"LPH-001-3_PD","LPH-001-9_AK2","LPH-001-10_PD","LPH-001-28_PD","LPH-001-12_AK2","LPH-001-14_PD","LPH-001-9_PD","LPH-001-16_PD","LPH-001-16_SCC","LPH-001-17_AK1","LPH-001-17_SCC",
"LPH-001-20_AK1","LPH-001-18_SCC","LPH-001-10_SCC","LPH-001-12_SCC")
length(unique(selected))

are.SCC<-
  sum(grepl("SCC$",selected))
  sum(grepl("AK2$",selected) | grepl("AK1$",selected) )
  sum(grepl("PD$",selected))

new<-rep.do.samples # [1:117]

new<-new[grepl("^LPH",new)]
new<-unique(new[!(new %in% selected)])
rep[1:5,]

rest<-pheno.ori[pheno.ori[,"SCC"],"SAMPLE"]

rest<-unique(rest[!(rest %in% selected)])

a.gene<-"INTS1"

gene.burden<-rep[rep[,"Gene.Names"]==a.gene,c("gene","Gene.Names","muts.in.cases","muts.in.PD")]
gene.burden
sum(duplicated(selected))

all<-paste(c(as.character(gene.burden[, "muts.in.cases"]),as.character(gene.burden[, "muts.in.PD"])),collapse=",")
all<-unique(unlist(strsplit(all,split=",")))

all



selected[selected %in% all]
selected[!(selected %in% all)]

new[new %in% all]
new[!(new %in% all)]

rest[rest %in% all]
rest[!(rest %in% all)]

get in some INTS1 "-" samples 

"LPH-001-17_SCC"->"LPH-001-15_SCC"
"LPH-001-16_SCC"->"LPH-001-27_SCC"
"LPH-001-9_PD"->"LPH-001-20_PD"
"LPH-001-20_AK1"->"LPH-001-20_AK2"
"LPH-001-12_AK2"->"LPH-001-30_AK1"
"LPH-001-9_AK2"->"LPH-001-7_AK2"
"LPH-001-24_AK1"->"LPH-001-15_AK1" 
"LPH-001-1_AK1"-> "LPH-001-28_AK1"

selected<-c("LPH-001-21_AK1","LPH-001-21_AK2","LPH-001-4_AK1","LPH-001-14_AK2","LPH-001-18_AK1","LPH-001-18_AK2","LPH-001-3_SCC","LPH-001-6_AK1","LPH-001-6_PD","LPH-001-14_AK1","LPH-001-15_AK2",
"LPH-001-16_AK1","LPH-001-25_AK1","LPH-001-3_AK1","LPH-001-3_AK2","LPH-001-4_AK2","LPH-001-4_SCC","LPH-001-7_PD","LPH-001-8_AK1","LPH-001-8_PD","LPH-001-12_PD","LPH-001-12_SCC",
"LPH-001-13_SCC","LPH-001-16_AK2","LPH-001-1_AK1","LPH-001-1_PD","LPH-001-1_SCC","LPH-001-21_PD","LPH-001-21_SCC","LPH-001-24_AK1","LPH-001-24_SCC","LPH-001-25_PD","LPH-001-26_SCC",
"LPH-001-3_PD","LPH-001-9_AK2","LPH-001-10_PD","LPH-001-28_PD","LPH-001-12_AK2","LPH-001-14_PD","LPH-001-9_PD","LPH-001-16_PD","LPH-001-16_SCC","LPH-001-17_AK1","LPH-001-17_SCC",
"LPH-001-20_AK1","LPH-001-18_SCC","LPH-001-10_SCC","LPH-001-12_SCC")

selected.exp<-c("LPH-001-21_AK1","LPH-001-21_AK2","LPH-001-4_AK1","LPH-001-14_AK2","LPH-001-18_AK1","LPH-001-18_AK2","LPH-001-3_SCC","LPH-001-6_AK1","LPH-001-6_PD","LPH-001-14_AK1","LPH-001-15_AK2",
"LPH-001-16_AK1","LPH-001-25_AK1","LPH-001-3_AK1","LPH-001-3_AK2","LPH-001-4_AK2","LPH-001-4_SCC","LPH-001-7_PD","LPH-001-8_AK1","LPH-001-8_PD","LPH-001-12_PD","LPH-001-12_SCC",
"LPH-001-13_SCC","LPH-001-16_AK2","LPH-001-28_AK1","LPH-001-1_PD","LPH-001-1_SCC","LPH-001-21_PD","LPH-001-21_SCC","LPH-001-15_AK1","LPH-001-24_SCC","LPH-001-25_PD","LPH-001-26_SCC",
"LPH-001-3_PD","LPH-001-7_AK2","LPH-001-10_PD","LPH-001-28_PD","LPH-001-30_AK1","LPH-001-14_PD","LPH-001-20_PD","LPH-001-16_PD","LPH-001-27_SCC","LPH-001-17_AK1","LPH-001-15_SCC",
"LPH-001-20_AK2","LPH-001-18_SCC","LPH-001-10_SCC","LPH-001-9_SCC")
length(unique(selected.exp))
selected.exp[duplicated(selected.exp)]
setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Gene_expression")
write.table(selected.exp,file="nano_string_slected_samples.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

selected.exp[selected.exp %in% all]
selected.exp[!(selected.exp %in% all)]

pheno.ori[pheno.ori[,"SCC"],"SAMPLE"] %in% selected





####################################### build a sliging window

x<-select.blood[1]

Sliding window analysis

pass<- full.qual &  maf.filter  & not.flat.genotype  & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal & pass.possonian.control.model & !bad.qual.locations  #  90263

snpinfo<-snpinfo.ori[pass,]
sort(table(snpinfo[,"gene"]))
/media/UQCCG/Sequencing/Data/Genomes/Exome Cature annotations/Nextera/nexterarapidcapture_exome_targetedregions_v1.2.bed
getwd()


save(list=c("snpinfo.ori","pass","summary.geno.extra","annotations"),file="core.RData")

library(skatMeta)  ## ridge regression
#library(SKAT) ## skat method
library(GenomicFeatures)
library(HardyWeinberg)
library(Biostrings)

NexteraRapidCapture_Exome_TargetedRegions_hg19_targets
setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis")
options(max.print=2000)
load("core.RData")
load("/media/UQCCG/Sequencing/Data/Genomes/NexteraRapidCapture_Exome_TargetedRegions_v1.2_hg19_targets.RData")

data.gr
data.gr<-data.gr+200
data.gr<-reduce(data.gr)
caploci<-as.data.frame(data.gr)

the.chrs<-table(caploci[,"seqnames"])
caploci[1:5,]
i<-1
annotations[1:5,]
length(pass)
dim(annotations)
sum(pass)
annoations.pass<-annotations[pass,]

window.all<-{}
for( i in 1:length(the.chrs)){
annoations.pass.chr<-annoations.pass[annoations.pass[,"chr"]==names(the.chrs)[i],]
print(dim(annoations.pass.chr))
loci = IRanges(start=as.numeric(annoations.pass.chr[,"start"]),end=as.numeric(annoations.pass.chr[,"end"]))
window<-loci+150
#window<-reduce(window)
# a.hist<-hist(width(window))
window

counts<-countOverlaps(window,loci)
length(counts)
window<-window[counts>1,]

#overlap<-findOverlaps(loci,window)
overlap<-findOverlaps(window,loci)

(overlap)[1:20]
queryHits(overlap)[1:10]
length(queryHits(overlap))
queryHits(overlap)

test.common<-tapply(subjectHits(overlap),queryHits(overlap),function(x) toString(x))
test.common[1:5]

dups<-duplicated(test.common)
test.common<-test.common[!dups]
test.common[1:5]

use.windows<-queryHits(overlap) %in% names(test.common)
use.windows[1:20]

gene<-paste(names(the.chrs)[i],queryHits(overlap)[use.windows],sep="_")
cluster<-gene
Names<-rownames(annoations.pass.chr)[subjectHits(overlap)[use.windows]]

a.window<-cbind(Names,gene,cluster)
if(i==1){window.all<-a.window}else{window.all<-rbind(window.all,a.window)}

} # i loop

window.all[1:5,]
dim(window.all)
length(window)
table(counts)
  annoations.pass
getwd()
colnames(window.all)[1]<-"Name"
save(list=c("window.all"),file="sliding_window_window.all.RData")
  
## cap.chr<-  caploci[annoations.pass,]
## the.sum<-cumsum(as.numeric(caploci[1:5,"width"]))
## step=500
## cap.ranges<-cbind(the.sum[-1*length(the.sum)],the.sum[-1])
## colnames(cap.ranges)<-c("start","end")
## ranges = IRanges(start=as.numeric(cap.ranges[,"start"]),end=as.numeric(cap.ranges[,"end"]))
## data.gr<-GRanges(seqnames =data[,"chr"],ranges = IRanges(start=as.numeric(data[,"start"]),end=as.numeric(data[,"end"])),strand=data[,"strand"])
## slide1.start<-seq(from=the.sum[1],to=tail(the.sum,n=1),by=step)
## slide1.end<-seq(from=(the.sum[1]+step),to=tail(the.sum,n=1),by=step)
## slide1.start<-slide1.start[-1*length(slide1.start)]
## slide<-cbind(slide1.start,slide1.end)
## colnames(slide)<-c("start","end")
## slide = IRanges(start=as.numeric(slide[,"start"]),end=as.numeric(slide[,"end"]))
## hist(width(data.gr))


/media/UQCCG/Sequencing/CompleteGenomics/917_Data_Delivery_2015.xlsx

/media/UQCCG/Sequencing/CompleteGenomics/Chort_descriptions/Sequencing comparisons-ver 9.xlsx
