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

########################## MAKE contaminated FILE 

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
run.per.chromsome<- FALSE ## set true if doing per chromosome
the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/BAM/LeoPharma_Feb2015.Sample_Sheet_NEW.txt"
## ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
contaminated.file<-c("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/LeoPharma_contaminated_march_31.txt") # COLL-FAM-1.201,COLL-FAM-1.3
#contaminated<-c()
contaminated<-read.table(contaminated.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)
remove.cols<-c()

core.ann<-c("chr","start","end","REF","ALT","TYPE") # out put to annanlsys programs and need foe colun labels
dont.build.summary<-TRUE ##
GATK.SB<-TRUE
maf.threshold.filter.to.use<-c(0.05)

a.label<-"coding_0.001"
dont.build.summary<-TRUE



################################## Other input files needed - path required in file names 
gene.symbol.file.for.clusters<-"/media/UQCCG/Software/annovar/humandb/Gene_symbol_aliases.txt"
cluster.definition.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-01-15_LeoPharma_Dec2014Freeze/Analysis/SCC_clusters.csv" ## tab delimited with header
ensembleID.to.HGNC.map.file<-"/media/UQCCG/Sequencing/Data/Genomes/hg19/ENSG_to_HGNC.txt" # tab delimited with header
coverage.file<-"/media/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_SAMPLE_coverage_LEO_2015.csv" # tab delimited with header
###########################
#/media/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_SAMPLE_Fri_Feb_06_2015.txt  latest



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

############################################# POPULATION MAF FILTER - PART A
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER

maf.threshold<-0.0  #  MAF threshold for annovar calling zero useful to get back all results !!do not modify!!
maf.threshold.filter.to.use<-c(0.001,0.005,0.01)
maf.threshold.filter.to.use<-sort(as.numeric(maf.threshold.filter.to.use))

filter.cols.novel.use<-c("NHBLI_6500_ANNOVAR_ALL","NHBLI_6500_ALL","NHLBI_5400_ALL","NHLBI_5400_EUR","NHLBI_5400_AFR","1000genome","1000genome_asian","1000genome_mine","snp141","snp141_clinical","snp137","CG69","EUR_ASN_AFR_INDEL","AOGC-NGS_ALL","AOGC-NGS_ALL_OLD","Chinese") ##
filter.cols.maf.use<-c("PopFreqMax","NHBLI_6500_ANNOVAR_ALL","NHBLI_6500_ALL","NHBLI_6500_EA","NHBLI_6500_AA","NHLBI_5400_ALL","1000genome","snp141","snp137","snp135","PopFreqMax","AOGC-NGS_ALL")
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
bad.samples<-contaminated[,1]
bad.samples.labels<-expand.labels.to.samples(bad.samples,c("GT","AD","DP","GQ"),paste.after=TRUE)

if(length(bad.samples.labels)>1){
a.indel<-a.indel[,colnames(a.indel)[!(colnames(a.indel) %in% bad.samples.labels)]]
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

grep("NOTCH1",names(all.genes))
common.hit.genes<-names(all.genes)[1:30]
common.hit.genes<-all.genes[all.genes>400]
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

library(doMC)
num.cores<-4
num.bits<-num.cores
registerDoMC(cores=num.cores)

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
pheno[1:5,]



the.samples<-paste(pheno[,"SAMPLE"],"GT",sep=".")  ## samples same order as in pheno
print(paste("Number samples: ",length(the.samples),sep=""))

#seq.type[1:57,]
posns<-match(pheno[,"SAMPLE"],seq.type[,"Sample"])
missing<-is.na(posns)
sum(missing)
pheno[missing,1:5]




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
## length(all)
table(pheno$SampleProject)

     ## AK Control  Normal      PD     SCC  # latestes with QC
     ## 47     133      24      25      23

cancer<-rep(FALSE,times=dim(pheno)[1])
cancer[pheno$SampleProject %in% c("AK","SCC")]<-TRUE

normal<-rep(FALSE,times=dim(pheno)[1])
normal[pheno$SampleProject %in% c("Normal")]<-TRUE

Control<-rep(FALSE,times=dim(pheno)[1])
Control[pheno$SampleProject %in% c("Control")]<-TRUE

PD<-rep(FALSE,times=dim(pheno)[1])
PD[pheno$SampleProject %in% c("PD")]<-TRUE


sum(PD)
sum(normal)
####"LPH-001-27 Blood is missing so will use PD
## "LPH-001-27_PD"
normal[pheno$SAMPLE %in% c("LPH-001-27_PD")]<-TRUE
PD[pheno$SAMPLE %in% c("LPH-001-27_PD")]<-FALSE


pheno<-cbind(pheno,cancer,Control,normal,PD)

     ## AK Control  Normal      PD     SCC 
     ## 47     133      24      25      23 

##############################################################

############## SET up groups to be analysed
table(pheno$SampleProject)

pheno[1:5,]

## th project becomes the master copy of the targets below.
the.projects<-c("cancer","Control","normal","PD")
names(the.projects)<-the.projects
colnames(pheno)


for (ir in 1: length(the.projects)){
  print(paste(the.projects[ir],"Num. samples:",sum(pheno[,the.projects[ir]])))
      }


##  > pheno[1:5,]
##   SampleProject FamilyCode SAMPLE PaternalID MaternalID Sex AffectionStatus cancer Control normal    PD
## 1             0        ALL    860          0          0  -9               1  FALSE    TRUE  FALSE FALSE
## 2             0        ALL    861          0          0  -9               1  FALSE    TRUE  FALSE FALSE
## 3             0        ALL    862          0          0  -9               1  FALSE    TRUE  FALSE FALSE
## 4             0        ALL   2529          0          0  -9               1  FALSE    TRUE  FALSE FALSE
## 5             0        ALL   2530          0          0  -9               1  FALSE    TRUE  FALSE FALSE                      

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



## #########################################
chk<-grep("^GENO",colnames(summary.geno.extra))
summary.geno.extra["chr19:50169104:50169104:C:T:snp",chk]
#

######################################


a.indel[1:5,1:10]
summary.geno.extra[1:5,]
colnames(summary.geno.extra)
#rownames(summary.geno.extra)<-key
#getHWE(obs_hets, obs_hom1, obs_hom2)
hw.target<-"Control"  ## what to calculate HW with

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



is.benign.missense<-is.missense & !qual[,"PolyPhen.low"]

sum(!is.benign.missense)
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
##########################################################################

p<-0.001
sd.thresh<-2

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
## dim(pheno.use)
if(!exists("pheno.ori")){
pheno.ori<-pheno
}
pheno<-pheno.ori[!pheno.ori[,"PD"],]
the.samples.use<-pheno[,"SAMPLE"]
the.samples.use<-paste(the.samples.use,".GT",sep="")
table(pheno[, "SampleProject"])

sum(pheno[,"PD"])


##########################################################################
##########################################################################
##########################################################################

####START HERE FELIICTY / TROELS
library(skatMeta)  ## ridge regression
#library(SKAT) ## skat method
library(GenomicFeatures)
library(HardyWeinberg)


code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")
source("hwe.r")



load("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/APRIL_1st_no_containiation_snpFilters_weights.RData")

### NOTE the above CONTAINS the filter and weights data already so YOU CAN jump to LINE 1724 (define "pass") starigh away unless you want to modify weights or additional filtering
##########################################################################
##########################################################################
##########################################################################
## extra filtering


filt<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Snp_read_filter/LEO_Pharma_WANTED_ALL_noGQ_filter_0.01_muts_Mar06.filter.pass.somatic.readfilter.txt",header=T,sep="\t",fill=TRUE,skip=5,stringsAsFactors=FALSE)
filt[1:5,1:20]

filt.key<-build.key(filt,core.ann)
rownames(filt)<-filt.key

pass.filt<-strsplit(filt[,"FILTER_SUMMARY"],split=";")
pass.filt[1:5]
pass.filt<-unlist(lapply(pass.filt,function(x) {sum(as.logical(x[1:3]))==0}))
pass.filt[1:5]
sum(!pass.filt)

snp.fail.filt<-rownames(filt)[!pass.filt]

snp.fail.filt[1:50] # list of bad snps 

filt["chr15:40675107:40675107:C:T:snp",1:10]
"chr15:40675107:40675107:C:T:snp" %in% snp.fail.filt

### don't sue the 4 filter picking up somatic calls in cases 
##                                   chr    start      end REF ALT TYPE   GENE         FILTER_SUMMARY                     SUMMARY_CALLED     SUMMARY_NOT_CALLED
## chr15:40675107:40675107:C:T:snp chr15 40675107 40675107   C   T  snp KNSTRN FALSE;FALSE;FALSE;TRUE 686;113;4;27;86;0;0.04;0.24;0.76;0 42339;182;5;23;159;1;8
       
## other.bad<-c("DNAH5","CSMD3","CSMD1","PCLO","INTS1","TYRO3","TTN","TESK1","MUC17","OR4A15","OR6C1","BCL2L11")

## other.bad %in% bad.genes

## bad.genes<-c(bad.genes,other.bad)


##########################################################################
##########################################################################
##########################################################################
library(robust)
library(robustbase)
library(nls2)
library(akima)
weights<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Felicity_filter/gene_synonymous_nonsynon_counts_above5.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
#weights<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Felicity_filter/gene_synonymous_nonsynon_counts_all.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

weights[1:5,]

 
               

#the.model

plot(as.numeric(weights[,2]) ,as.numeric(weights[,3]),main="For C>T motif SNPs",xlab="Synonymous",ylab="NonSynonymous",cex.lab=1.5)
identify(as.numeric(weights[,2]) ,as.numeric(weights[,3]),labels=weights[,1])

the.model<-try(lmrob(NonSynonymous~Synonymous,data=weights,silent=FALSE)) # ,subset=c(1:dim(test)[1])[!red])
abline(coef(the.model),lty=10,col="red",lwd=2)

the.model<-try(lmrob(NonSynonymous~Synonymous,data=weights,subset=c(1:dim(weights)[1])[!(weights$Gene %in% c("TTN"))],silent=FALSE) ) #
abline(coef(the.model),lty=10,col="green",lwd=2)

the.model<-try(lmrob(NonSynonymous~Synonymous,data=weights,subset=c(1:dim(weights)[1])[weights$Synonymous >=10],silent=FALSE) )
               abline(coef(the.model),lty=10,col="red",lwd=2)

use.wieght<-coef(the.model) ## use THIS weighting
                abline(c(0,coef(the.model)[2]),lty=10,col="cyan",lwd=2)

the.model<-try(lmrob(NonSynonymous~Synonymous-1,data=weights,subset=c(1:dim(weights)[1])[weights$Synonymous >=10],silent=FALSE) )
               abline(c(0,coef(the.model)),lty=10,col="purple",lwd=2)


the.model<-try(lmrob(NonSynonymous~Synonymous,data=weights,subset=c(1:dim(weights)[1])[(weights$Synonymous >=10) & !(weights$Gene %in% c("TTN"))],silent=FALSE) )
 abline(coef(the.model),lty=10,col="blue",lwd=2) 
               
  abline(coef(the.model),lty=10,col="red",lwd=2)

savePlot(filename=paste("UV_wieght_model.png",sep=''),type="png")
savePlot(filename=paste("UV_wieght_model.jpeg",sep=''),type="jpeg")
savePlot(filename=paste("UV_wieght_model.tiff",sep=''),type="tiff")
#savePlot(filename=paste(fig.prefix,"aa.bmp",sep=''),type="bmp")
dev.print(svg,paste("UV_wieght_model.svg",sep=''))
getwd()

#weights<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Felicity_filter/gene_synonymous_nonsynon_counts_all.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
use.wieght
weights[1:5,]
min.val<-(1-use.wieght[1])/use.wieght[2]
min.val # check not getting negative weights
the.weights<-weights$Synonymous*use.wieght[2] + use.wieght[1]
names(the.weights)<-weights$Gene
max(the.weights)
min(the.weights)
the.weights<-1/the.weights

sort(the.weights)
plot(the.weights)
if("the.weights" %in% colnames(weights)){
  weights[,"the.weights"]<-the.weights
}else{
weights<-cbind(weights,the.weights)
}
##########################################################################
##########################################################################
##########################################################################
types<-c("all.somatic","coding.somatic","coding")
types<-c("coding.somatic","coding.somatic.with.Indels","coding.somatic.with.Indels.noBenign","coding.somatic.with.Indels.noBenign.wRepeats","synonymous.with.Indels.noBenign.wRepeats","coding.somatic.with.repeats")

itypes<-4 #"coding.somatic.with.Indels.noBenign.wRepeats"

       
## for(itypes in 1:length(types)){

## snap.file<-types[itypes]
snap.file
paste(snap.file,".withGQ.filter.QCed",sep="")

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

## if(types[itypes]=="coding.somatic.with.Indels.noBenign"){
  
## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype & !are.repeats & !are.in.repeats & ok.missing & hw.controls.ok.filt & !no.genotypes &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal & !is.benign.missense # &  &  snp.only # & !on.x.y #  42004

## }

## if(types[itypes]=="coding.somatic.with.Indels.noBenign.wRepeats"){
  
pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  &
ok.missing & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control.filt & !in.any.normal.filt & !in.any.normal & !is.benign.missense # &  &  snp.only # & !on.x.y #  42004


## pass sets which snps will be used
## *.filt derived form "filtered" genothes (coverage and GQ)
## full.qual= filter=pass
## bad.coding<-protein changing
## af filetr < 0.001 poulation filter appled to controls
## sum(in.common.hit.gene ) no used here
## unannotated.hits  no gene name
## not.flat.genotype not the flat multi allle snps so don't double count

## ok missing < 20% missing
## hw.controls.ok.filt HW atfer coverge and GQ filter


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

help<-cbind(full.qual,bad.coding,maf.filter,rare.in.group,no.genotypes,in.common.hit.gene ,hw.controls.ok,on.x.y,unannotated.hits,not.flat.genotype,are.repeats,are.in.repeats,ok.missing,ok.missing,is.unwound.geno,is.unwound.geno ,hw.p.control.filt,is.benign.missense)

#

## pass<-a.indel[,"Gene.Names"] %in% c("TP53","NOTCH1","NOTCH2","ATM","ACD","ASIP","BAP1","CASP8","CCND1","CDK4","MC1R","MITF","MTAP","MX2","OCA2","PARP1","PLA2G6","POT1","SLC45A2","TERF2IP","TERT","TYR","TYRP1","VDR","BCL2L12","KNSTRN","ISX","CDKN2A","BCL2L11","STK19","FJX1","TRHDE")
## sum(pass)



## help<-cbind(full.qual,bad.coding,maf.filter,rare.in.group,no.genotypes,in.common.hit.gene ,hw.controls.ok,on.x.y,unannotated.hits,not.flat.genotype,are.repeats,are.in.repeats,ok.missing,ok.missing,is.unwound.geno,is.unwound.geno ,hw.p.control.filt,rare.in.group.filt,no.genotypes.filt,rare.in.controls.filt )
## pass<-full.qual & functional & maf.filter & rare.in.group & !no.genotypes  & not.flat.genotype & !(high.total.missing | nimblegen.total.missing | illumina.total.missing)
## sum(pass)
## loc<-grep("chr19:50169131:50169131:C:T:snp",key)
## pass[loc]
## help[loc,]
## summary.geno.extra[loc,]
## in.any.normal[loc]

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

############ EXTRA FILTERING FROM FELICITY
"chr15:40675107:40675107:C:T:snp" %in% snp.fail.filt
sum((snp.fail.filt %in% names(pass)))
pass[ names(pass) %in% snp.fail.filt]<-FALSE
#################### FORCE GOOD CALLS

#pass[ names(pass) %in% "KNSTRN"]

good.genotypes<-c("chr19:50169131:50169131:C:T:snp","chr19:50169104:50169104:C:T:snp","chr19:50169132:50169132:C:T:snp")
good.genotypes %in% names(pass)
good.genotypes %in% snp.fail.filt
pass[ names(pass) %in% good.genotypes]<-TRUE

############## add weights
the.weights ## weights for gene with more than 5 syn mutations
snpinfo.ori[1:5,]
gene.weights<-rep(1,times=dim(snpinfo.ori)[1])
names(gene.weights)<-snpinfo.ori[,"gene"]
gene.weights[1:10]
posns<-match(names(gene.weights),names(the.weights))
missing<-is.na(posns)
sum(missing)
cbind( names(gene.weights)[!missing],gene.weights[!missing],the.weights[posns[!missing]])[1:100,]
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

dim(snpinfo.ori)
length(gene.weights)


pass[unique(snpinfo.ori[snpinfo.ori[,"gene"]=="KNSTRN","Name"])]

sum(pass)
################################# GEFOS FILTERING cause sending all
#pass<-pass[the.snps] ### GEOFS 
genotypes<-a.indel[pass,the.samples.use] ## ordered correctly for phenotypes
snp.names<-key[pass] ## GEFOS ony name with start

#### snpinfo now A different size than a.indel since added pathways!!!  snpinfo[snpinfo[,"gene"]=="KNSTRN",]
snpinfo<-snpinfo.ori[snpinfo.ori[,"Name"] %in% snp.names,]
gene.weights.subset<-gene.weights[snpinfo.ori[,"Name"] %in% snp.names] # weight in same order as snpinfo.ori
snpinfo<-cbind(snpinfo,gene.weights.subset)
snpinfo[1:5,]


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
 cohort.seq <- skatCohort(Z=genotypes,formula, SNPInfo = snpinfo, data=pheno,aggregateBy="cluster",family=binomial(),verbose=FALSE) # for Leo case control

}else{
cohort.seq <- skatCohort(Z=genotypes,formula, SNPInfo = snpinfo, data=pheno,aggregateBy="cluster",family=gaussian(),verbose=FALSE) ## genes and clusters
}

meta.results.burden<-burdenMeta(cohort.seq,wts="gene.weights.subset",mafRange = c(0,1),SNPInfo = snpinfo,aggregateBy="cluster")
# meta.results.burden<-burdenMeta(cohort.seq,wts=1,mafRange = c(0,1),SNPInfo = snpinfo,aggregateBy="cluster") OLD VESION
#meta.results.skat<-skatMeta(cohort.seq,SNPInfo = snpinfo,aggregateBy="cluster")
meta.results.skatO<-skatOMeta(cohort.seq,burden.wts =1,SNPInfo = snpinfo,aggregateBy="cluster",method = "integration")



## meta.results.skatO<-{}

the.order<-     order(meta.results.burden[,"p"])
sum(is.na(meta.results.burden[,"p"])) ## bad p-values shoudl not happen
meta.results.burden<-  meta.results.burden[the.order,]
meta.results.burden[1:50,]
## ## meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]

meta.results.skat<-{}
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


# annotations<-a.indel[,c(1:6,16,28,7,30,34,37:42,43,14,32,33)]

## save(list=c("meta.results.skat","meta.results.skatO","meta.results.burden","pheno.use","snpinfo","genotypes","pass","high.missing","annotations","help","key","summary.geno.extra"),file=paste(paste(project.files[ichr],".small.RData",sep="")) )



save(list=c("gene.weights","gene.weights.subset","filt","snp.fail.filt","use.wieght","weights","core.ann","case.control","snpinfo.ori","formula","clusters","pheno.types","ipheno","clusters.wanted","genotypes","p","meta.results.skat","meta.results.skatO","meta.results.burden","pheno","pheno.ori","target.pheno.col","snpinfo","fil.genotypes","pass","high.missing.table","a.indel","help","key","summary.geno.extra","full.qual","bad.effect","maf.filter","in.common.hit.gene","on.x.y","unannotated.hits","not.flat.genotype","are.repeats","are.in.repeats","ok.missing","hw.controls.ok.filt","no.genotypes","rare.in.Control","rare.in.Control.filt","in.any.normal","in.any.normal.filt","are.in.repeats.back","are.in.repeats.forward","all.genes","contaminated"),file=paste(paste(project.files[ichr],".",pheno.types[ipheno],".",snap.file,".small_final.RData",sep="")) )

} # itype

print(paste("Done: ",pheno.types[ipheno],"->",project.files[ichr]))
save.image(file=paste(project.files[ichr],".",pheno.types[ipheno],".",snap.file,".RData",sep=""))

} # pheno loop

} # loop over projects


} # loop over fam


 getwd()
[1] "/media/old-scratch/media/scratch2/AOGC-NGS/Analysis"
load("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/ALL_with_synon_Mar11_2015.RData")

colnames(a.indel)[1:50]
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
snap.file<-"synonymous.coding.0.001.snpOLNY"


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

## net<-c("BCL2L12","KNSTRN","ISX","CDKN2A","BCL2L11","STK19","FJX1","TTN","MUC4")
## file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-01-15_LeoPharma_Dec2014Freeze/Analysis/String/INTS1_fill_network.txt"
## file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-01-15_LeoPharma_Dec2014Freeze/Analysis/String/spring.ori.2.more.txt"

## net<-read.delim(file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
## dim(net)
## net[1:5,]
## net<-c(net[,1],net[,2])

## sort(table(net),decreasing=TRUE)
## net<-unique(net)
## posns<-match(net,meta.results.burden[,"gene"])
## names(posns)<-net
## sort(posns)

## posns<-match(net,meta.results.skatO[,"gene"])
## names(posns)<-net
## sort(posns)

## meta.results.burden[meta.results.burden[,"gene"] %in% net,]

## check<-grepl("^MED",meta.results.burden[,"gene"])
## check<-grepl("^MED",meta.results.burden[,"gene"])


## meta.results.burden[,"gene"][check]

## posns<-grep("KNSTRN",a.indel[,16],as.is=TRUE)

## a.indel[posns,1:20]


## meta.results.burden[1:350,]
## meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]

####################### summaruze at VARIANT LEVEL set by unwind
########
the.top<-350
to.unwind<-c(meta.results.burden[1:the.top,"gene"],meta.results.skatO[1:the.top,"gene"])


posn.get<-rep(FALSE,times=dim(meta.results.burden)[1])
posn.get[1:the.top] <- TRUE
sum(posn.get)
wanted<-posn.get & meta.results.burden[,"nsnpsUsed"]<100
sum(wanted)

to.unwind<-c(meta.results.burden[wanted ,"gene"])
net<-
## to.unwind<-c("BCL2L12") #, "MCM7", "RNPC3")
## to.unwind<-c("FANC_complex.all") # to.unwind<-meta.results.burden[8,"gene"]
 to.unwind<-c("BCL2L12","KNSTRN","ISX","CDKN2A","BCL2L11","STK19","FJX1","TRHDE") #"test.set"
to.unwind<-c("TP53","NOTCH1","NOTCH2","FAT4","STK19","ISX","TRHDE","ARHGAP35","PREX1","KL","PIK3CA","KNSTRN","BCL2L12")
to.unwind<-c("TP53","NOTCH1","NOTCH2","ATM","ACD","ASIP","BAP1","CASP8","CCND1","CDK4","MC1R","MITF","MTAP","MX2","OCA2","PARP1","PLA2G6","POT1","SLC45A2","TERF2IP","TERT","TYR","TYRP1","VDR","BCL2L12","KNSTRN","ISX","CDKN2A","BCL2L11","STK19","FJX1","TRHDE")

#to.unwind<-c("KNSTRN") #,"TCEA1","POLR2A","CTDP1")
## to.unwind<-c("Clinical")
#  to.unwind<-c(clusters.wanted[!(clusters.wanted %in% c("Ubin.proteo","lipid_raft","caveolae","Checkpoint_extendedx1","Checkpoint_extendedx2"))])
#grep(to.unwind,meta.results.burden[,"gene"])
#to.unwind
# to.unwind.name<-to.unwind[1]
to.unwind.name<-"TOP_350"
#match(net,meta.results.burden[,"gene"])
# to.unwind.name<-"Aideen.set"
# to.unwind.name<-"Hot.set.PASS"
# to.unwind.name<-"ALL_significant"
# to.unwind.name<-"ALL_significant"

to.unwind.name
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
the.genes.burden<-meta.results.burden[meta.results.burden[,"gene"] %in% the.genes,]

the.genes.burden
write.table(the.genes.burden,file=paste(to.unwind.name,"conponents:","Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

the.genes.burden<-meta.results.skatO[meta.results.skatO[,"gene"] %in% the.genes,]
the.genes.burden
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
#sum(are.in.repeats[loci])

    
# qual[loci,]
## snpinfo[1:5,]
## qual[1:5,c("FILTER_PASS", "FILTER_100" )]

cohort.seq.ex <- skatCohort(genotypes.ex, formula, SNPInfo = snpinfo.ex, data=pheno,aggregateBy = "Name",verbose=FALSE)
## meta.results.skat.ex<-skatMeta(cohort.seq,SNPInfo = snpinfo)
meta.results.burden.ex<-burdenMeta(cohort.seq.ex,wts=1,mafRange = c(0,1),SNPInfo = snpinfo.ex,aggregateBy = "Name")
meta.results.burden.ex[1:5,]
pheno[1:5,]



cohort.seq.test <- skatCohort(genotypes.ex, formula, SNPInfo = snpinfo.ex, data=pheno,aggregateBy = "cluster",verbose=FALSE)

meta.results.burden.test<-burdenMeta(cohort.seq.test,wts=1,mafRange = c(0,1),SNPInfo = snpinfo.ex,aggregateBy = "cluster")
#meta.results.burden.test

## meta.results.skat.ex<-skatMeta(cohort.seq,SNPInfo = snpinfo)
#meta.results.skatO.test<-skatOMeta(cohort.seq.test,burden.wts =1,SNPInfo = snpinfo.ex,aggregateBy="cluster")
#meta.results.skatO.test



muts.in.cases<-apply(genotypes.ex[pheno[,"cancer"],],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
muts.in.controls<-apply(genotypes.ex[pheno[,"Control"],],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})


 controls<- paste(pheno[pheno[,"Control"],"SAMPLE"],".GT",sep="")
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

                             
out<-cbind(meta.results.burden.ex,annotations[figure,],a.functions[figure,],is.benign.missense[figure],annotations[figure,],summary.geno.extra[figure,colnames(summary.geno.extra)[grep("^GENO",colnames(summary.geno.extra))]],help[figure,],muts.in.cases,quality.cases,depth.cases,muts.in.controls,quality.controls,depth.controls)

#out<-cbind(meta.results.burden.ex,annotations[figure,],muts.in.cases,muts.in.controls)
dim(out)
out[1:4,]

## table(out[,"refGene::location"])
## table(out[,"Consequence.Embl"])

paste(paste(to.unwind,collapse="."))
paste(to.unwind.name,collapse=".")
  paste(paste(to.unwind.name,collapse="."),"GENOTYPE.conponents:","SkatO","clusters",snap.file,"txt",sep=".")
write.table(out,file=paste(paste(to.unwind.name,collapse="."),"GENOTYPE.conponents:","Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

getwd()
setwd(analysis.dir)
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################






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







################ Q-Qplot of data
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
filt<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/TOP_350.GENOTYPE.conponents-.Burden.clusters.coding.somatic.with.Indels.noBenign.filter.update.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
filt[1:5,]
write.table(filt[1:350,],file="string.txt",,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


do<-filt[filt[,"Gene.Names"] %in% interesting,]

do[1:5,]

i<-1
ckk<-interesting[i]
tapply(

bad<-filt[,"STRAND_BIAS"]=="true" |  filt[,"MISMATCH"]=="true"  | filt[,"READ_END"]=="true"
tapply(filt[bad,"Gene.Names"],bad,length)
bad.genes<-filt[bad,"Gene.Names"]
bad.genes
table(bad.genes)
bad

interesting<-c("TP53","NOTCH1","NOTCH2","FAT4","STK19","ISX","TRHDE","ARHGAP35","PREX1","KL","PIK3CA","KNSTRN","BCL2L12")

interesting[interesting %in% bad.genes]
bad.genes<-bad.genes[!(bad.genes %in% interesting)] # checked this is ok as were just single point mutations


other.bad<-c("DNAH5","CSMD3","CSMD1","PCLO","INTS1","TYRO3","TTN","TESK1","MUC17","OR4A15","OR6C1","BCL2L11")

other.bad %in% bad.genes

bad.genes<-c(bad.genes,other.bad)

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

## subset<-meta.results.burden[ !(meta.results.burden[,"gene"] %in% genes.and.clusters) ,]


## bad.genes<-c(bad.genes,clusters.wanted)
subset<-meta.results.burden[ !(meta.results.burden[,"gene"] %in% bad.genes) ,]
subset<-subset[!is.na(as.numeric(subset[ ,"p"])),]

subset[1:50,]
posns<-match(interesting,subset[,"gene"])
names(posns)<-interesting
posns
write.table(subset[1:350,],file="string.txt",,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
dups<-duplicated(subset[,"gene"])
sum(dups)

write.table(subset[1:350,],file="string.txt",,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

z<-qchisq(subset[ ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) ## if have no chisq valuse
#z0<-qchisq(meta.results.skatO[ !(meta.results.skatO[,"gene"] %in% genes.and.clusters ) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) 
z[1:5]
subset[1:10,]

p.val<-as.numeric(subset[ ,"p"])
par(mfrow=c(1,1),font=2,font.lab=2,font.axis=2,mgp=c(3.5,1,0),mar=c(5,5,4,2)+0.1)


##  z<-rchisq(length(p.val), df=1, ncp = 0) ## null test


median(z,na.rm=TRUE)/0.456  #1.071491
median(z0,na.rm=TRUE)/0.456  #1.071491


sum(is.na(p.val))
################## p-values
z0=qnorm(p.val[!is.na(p.val)]/2)
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


the.plot<-my.qq.plot(z[!is.na(z)],dist="chisq",df=1,ncp=0,col="blue",ylab="Observed chi-squared value",xlab="Expected chi-squared value",main="",cex=1,xlim=c(0,22),ylim=c(0,90),cex.lab=2.0,cex.axis=2.0,font.lab=2,font.axis=2,lwd=2,line="robust",plot.it=TRUE) # function defined below


## z.all<-qchisq(meta.results.burden[ !(meta.results.burden[,"gene"] %in% clusters.wanted) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)
## range(z.all)
qq<-  qq.data(z,plot.it=FALSE)       ## qq plot used same method as in car library
#points(qq$x,qq$y,col="magenta",pch=21)

symbols<-subset[,"gene"]
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="magenta",pch=25)

#symbols<-meta.results.skatO[,"gene"]
#####annotate curve
selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="red",cex=1,offset=1,atpen='TRUE') ##plate row col symbol

selected.data<-identify(qq$x,qq$y,labels=qq$ord,col="red",cex=1,offset=1,atpen='TRUE')


selected.data<-identify(qq$x,qq$y,labels=labels[qq$ord],col="red",cex=1,atpen='TRUE') ## sybmol
selected.data<-identify(qq$x,qq$y,labels=as.character(round(data.in[qq$ord],2)),col="forestgreen",cex=1.25,atpen='TRUE') # observed score
#####


## leg.txt<-c("All Genes","Remove Clinical Genes")

## legend(2,60,leg.txt,col=c("magenta","blue"),lty=c(-1,-1,2),pch=c(1,1,-1),cex=2.25)

## label<-"AML_paper_0.01_all_filters_NO_STRAT"

leg.txt<-c(expression(paste("All Genes ")),"95% confidence intervals")
legend(2,80,leg.txt,col=c("blue","red"),lty=c(-1,2),pch=c(1,-1),cex=2.25)

setwd( analysis.dir)
label<-"SCC1"
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




















