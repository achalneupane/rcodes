

# 
# ###### THIS IS the super annotion run USE ALL THE FILTER AND NOVEL DATABASES
# UQCCG.data<-"/mnt/ga-apps/UQCCG/Data/Sequence_Genotypes" # path to project data 
# project<-"2013-10-27_AML_with_AOGCControl_NoFailedLane" # this is the project directory
# project.name<-"2013-10-27_AML_with_AOGCControl_NoFailedLane" # this is the prefix of the SNP and DINDEL files AND the project output NAME
# 
# annotate.dir<-"/mnt/ga-apps/UQCCG/Data/Sequence_Genotypes/2013-10-27_AML_with_AOGCControl_NoFailedLane/Annotate"
# analysis.dir<-"/mnt/ga-apps/UQCCG/Data/Sequence_Genotypes/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis"
# 
# affection.status<-"sample_sheet"
# the.sample.sheet<-"TGCM-AML_SampleSheet.csv"
# restrict.family<-c()
# 
# filter.cols.novel<-c("NHBLI_6500_ANNOVAR_ALL","NHBLI_6500_ALL","NHLBI_5400_ALL","NHLBI_5400_EUR","NHLBI_5400_AFR","1000genome","1000genome_asian","1000genome_mine","snp137","snp137_clinical","snp135","CG69","EUR_ASN_AFR_INDEL","AOGC-NGS_ALL","AOGC-NGS_ALL_OLD","Chinese") ##
# filter.cols.maf<-c("NHBLI_6500_ANNOVAR_ALL","NHBLI_6500_ALL","NHBLI_6500_EA","NHBLI_6500_AA","NHLBI_5400_ALL","1000genome","snp137","AOGC-NGS_ALL")
# 
# GATK.SB<-TRUE
# genome.build<-"hg19"
# Skip.Analysis<-FALSE
# 
# 
# alt.allele.other.use<-TRUE
# alt.allele.other.use.threshold<-0
# make.nonDefined.samples.as.OTHER<-TRUE
# include.homozygotes<-TRUE ## used for cancer
# ###############################################

args <- commandArgs(TRUE)
ichr <- args[1]
ichr <- as.numeric(ichr)
#ichr <- 10
# ichr <- 22
#the.chr.index<-c(1:24)
#for(ichr in 1:14){
 

##### THIS IS the super annotion run USE ALL THE FILTER AND NOVEL DATABASES
 UQCCG.data<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes"
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/genotyping"

 project<-"2014-03-19_MND_MODY_LKAF_Nimbelgen" # this is the project directory
#project.name<-"20131115_GroupCallSKDP_MODY_PCC_ChMND" # this is the prefix of the SNP and DINDEL files AND the project output NAME
 project.name <- project
 #annotate.dir<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/20131115_GroupCallSKDP_MODY_PCC_ChMND/Annotate"
# analysis.dir<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/20131115_GroupCallSKDP_MODY_PCC_ChMND/Analysis"
# 
annotate.dir <-paste(UQCCG.data,project,"Annotate",sep="/")
analysis.dir<-paste(UQCCG.data,project,"Analysis",sep="/")

 affection.status<-"sample_sheet"
 the.sample.sheet<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-03-19_MND_MODY_LKAF_Nimbelgen/BAM/Ch_MND_SampleSheet.csv"
 restrict.family<-c("Ch_MND_F")



# # UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/genotyping//"
# UQCCG.data<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/"
#  project<-"2014-01-06_Dec2013SkdpGenotyped" # this is the project directory
#  project.name<-"2014-01-06_Dec2013SkdpGenotyped" # this is the prefix of the SNP and DINDEL files AND the project output NAME
# # 
#  annotate.dir<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-01-06_Dec2013SkdpGenotyped/Annotate"
#  analysis.dir<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-01-06_Dec2013SkdpGenotyped/Analysis"
# # 
#  affection.status<-"sample_sheet"
#  the.sample.sheet<-"SKDP_SampleSheet.csv"
#  restrict.family<-c()


# ##### THIS IS the super annotion run USE ALL THE FILTER AND NOVEL DATABASES
# UQCCG.data<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes//"
# project<-"2013-11-17_GroupCallIYMD_Lung_AOGC_Controls" # this is the project directory
# project.name<-"20131117_GroupCallIYMD_Lung_AOGC_Controls" # this is the prefix of the SNP and DINDEL files AND the project output NAME
# 
# annotate.dir<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes//2013-11-17_GroupCallIYMD_Lung_AOGC_Controls/Annotate"
# analysis.dir<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes//2013-11-17_GroupCallIYMD_Lung_AOGC_Controls/Analysis"
# 
# affection.status<-"sample_sheet"
# the.sample.sheet<-"IYMD_SampleSheet_Project.csv"
# restrict.family<-c()




code.dir<-"/mnt/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/"



filter.cols.novel<-c("NHBLI_6500_ANNOVAR_ALL","NHBLI_6500_ALL","NHLBI_5400_ALL","NHLBI_5400_EUR","NHLBI_5400_AFR","1000genome","1000genome_asian","1000genome_mine","snp137","snp137_clinical","snp135","CG69","EUR_ASN_AFR_INDEL","AOGC-NGS_ALL","AOGC-NGS_ALL_OLD","Chinese") ##
filter.cols.maf<-c("NHBLI_6500_ANNOVAR_ALL","NHBLI_6500_ALL","NHBLI_6500_EA","NHBLI_6500_AA","NHLBI_5400_ALL","1000genome","snp137","AOGC-NGS_ALL")

GATK.SB<-TRUE
genome.build<-"hg19"
Skip.Analysis<-FALSE


alt.allele.other.use<-TRUE
alt.allele.other.use.threshold<-0
#control.ONLY.in.other<-TRUE
include.control.in.other<-TRUE
include.homozygotes<-TRUE ## used for cancer
###############################################



#############################################################################################################
############################################# BEGIN #######################################################
options(stringsASFactors=FALSE)
options(width=250,max.print=5000)

library(GenomicFeatures)
library(Rsamtools)
#require(biomaRt)
library(DBI)
library(RMySQL)
library(multicore)
library(foreach)  ## could not get this to work well
library(doMC)
if(!exists("num.cores")){num.cores<-12}
registerDoMC(cores=num.cores)

db.details<-list(user="gerp2",pass="GeRpUQ8!6", dbname="gerp", host="ga-apps.di.uq.edu.au")

minor.allele.DB<-"snp137" # This database only has minor allele frquecies and not alternative allele frequencies
primary.snp.annotation.DB<-"snp137"
secondary.snp.annotation.DB<-"snp132"
######################################################### Set up the basics for each run #####
project.extension<-".table.RData"
project.dir<-paste(UQCCG.data,project,sep="/")
bam.dir<-paste(project.dir,"BAM",sep="/")
polyphen.dir<-paste(project.dir,"PolyPhen",sep="/")

QC.dir<-paste(project.dir,"QC",sep="/")
if(!exists("annotate.dir")){annotate.dir<-"/mnt/Bioinform-D/Research/annotate"}

if(!exists("analysis.dir")){analysis.dir<-annotate.dir}
if(!exists("Skip.Analysis")){Skip.Analysis<-FALSE}
if(!exists("restrict.family")){restrict.family<-{}}
if(!exists("force.redo")){force.redo<-TRUE}

if(!exists("alt.allele.other.use")){alt.allele.other.use<-TRUE} ## if true ALL other samples not in current project/family are "other"
################ CAN RUN IN THREE MODES
if(!exists("include.control.in.other")){include.control.in.other<-FALSE} # any sample NOT in current project/family AND has is labelled   CONTROL or control on SampleProject
if(!exists("control.ONLY.in.other")){control.ONLY.in.other<-FALSE} # ONLY samples listed as CONTROL or control in SampleProject are used as other
if(!exists("make.nonDefined.samples.as.OTHER")){make.nonDefined.samples.as.OTHER<-FALSE}  ## other samples not listed in the sample sheet are used at other
### if ALL of above are FALSE and alt.allele.other.use then no alt allele threshold is applied to gene hit counts. BUT alt alleles.OTHER is still reported and is ALL samples not in current project/family
## This is the "usual" run more
#######################
if(!exists("alt.allele.other.use.threshold") & alt.allele.other.use ){alt.allele.other.use.threshold<--9}  # if alt.allele.other.use.threshold not set ,set to max maf ->thesshold= no.other.samples*2 *maf

if(!exists("use.only.GATK.pass")){use.only.GATK.pass<-FALSE} # if true in gene/sample summaries only use TRUE
if(!exists("include.homozygotes")){include.homozygotes<-FALSE}  ## CANCER STUDIES in gene hit and sample hit summay count 0/1 and 1/1 (both count as 1)



## anno.DB.location.core<-"/mnt/Bioinform-D/Research/annovar/humandb"
## anno.DB.location<-paste(anno.DB.location.core,genome.build,sep="/")
#anno.DB.location<-"/mnt/Bioinform-D/Research/annovar/humandb/hg19" ## location of the annotation database

gerp.element.threshold<- 1 # any hit is significant
gerp.score.threshold<-2.5 # gerp score >= will be included
gerp.score.threshold.low<-2 # gerp score >= will be included

force.functional.update<-FALSE # will re-real functional data normally set to FALSE
force.GERP.update<-FALSE  # will reget get scores normally set to FALSE
update.annotations.with.vep<-FALSE ## TRUE normally undate gene names and geneanno.table with filter.table.poly... wanted muts traken care of below
update.gene.names<-FALSE # TRUE normally
skip.gene.analysis<-FALSE # TRUE normally



splice.threshold<-5 # default is 2 in annovar. I have seen  3 and 5 are used in dbSNP
maf.threshold<-0.0  #  MAF threshold for annovar calling zero useful to get back all results !!do not modify!!
maf.threshold.filter.to.use<-c(0.001,0.005,0.01,0.025) #c(0.001,0.01,0.025,0.05) MAF threshold for boolean filtering operations NOVEL always done
# maf.threshold.filter.to.use<-c(0.001,0.005) # cancer sequence
#maf.threshold.filter<-0.05  maf.threshold.filter.to.use<-c(0.5,1)
core.ann<-c("chr","start","end","REF","ALT","TYPE") # out put to annanlsys programs and need foe colun labels

###############

if(alt.allele.other.use.threshold==-9){alt.allele.other.use.threshold<-max(maf.threshold.filter.to.use)}


#PASS  VQSRTrancheINDEL99.00to99.90 VQSRTrancheINDEL99.90to100.00    VQSRTrancheSNP99.00to99.90   VQSRTrancheSNP99.90to100.00 
global.quality.labs<-c("QUAL","QD","FS","HRun","SB","FILTER","FILTER","TYPE") ### these become the "good.qual" filter
global.quality.names<-c("QUAL","QD","FS","HRun","SB","FILTER_PASS","FILTER_100","flat")
global.quality<-c(50,0.5,60,5,1,"PASS","99.90to100.00","flat")
global.quality.type<-c("numeric","numeric","numeric","numeric","numeric","factor","factor","factor")
global.quality.dirn<-c("greater","greater","less","less","less","exact","ends_with","ends_with")

secondary.filter<-c("QUAL","QD","FS","HRun","SB") # names in global.quailty.filter that are TRUE : list all wanted ones not present in data are ignored

indiv.quality.labs<-c("GT","GT","GT","DP","DP","DP")  ### column ingenomic range object
indiv.quality.names<-c("GT_A","GT_AB","GT_B","DP_High","DP_Low","DP_Thresh") ## anme that will be use for that filter
indiv.quality<-c("0/0","0/1","1/1",7,5,5)  ## the filter to apply : factors "==" used numeric ">" used (geather than)
indiv.quality.type<-c("factor","factor","factor","numeric","numeric","numeric")
indiv.quality.dirn<-c("exact","exact","exact","greater","greater","less")


####################################


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


#generic.filter.DB
maf.threshold.filter.to.use<-sort(maf.threshold.filter.to.use) # cause when subset the final data use the largest MAf liste

setwd(code.dir)
## source("ucsc.table.names.r") 
                                        # load in the UCSC tables these use the db file names and not their lable-names 
source("annotate_SNPs_subroutines.r")

##################################################################################                   
######################################### Predefined variables required
##################################################################################


setwd(annotate.dir)
files<-dir(annotate.dir)
the.extension<-paste(project.extension,"$",sep="")
project.files<-files[grepl(paste("^",project.name,sep="") ,files) & grepl(the.extension ,files)]
if(length(project.files)==0){project.files<-files[grepl(the.extension ,files)] } ### use if no project files


print(paste("Doing: ",project.files,sep="")) # project.files<-project.files[c(-15)]
 #ichr<-21
print(project.files)
print(length(project.files))
print("Chr")
print(ichr)
print("Project files...")
print(project.files[ichr])
#q()
#for(ichr in 1:length(project.files)){

  setwd(annotate.dir)
  print(paste("Loading: ",project.files[ichr],sep=""))
  load(project.files[ichr])
  get<-try(load(paste(gsub(".txt$","",target),".GeneLists.RData",sep="")),silent=TRUE)
 ## setwd(annotate.dir)
 ## load(paste(project.name,".table.RData",sep=""))  # load data from the annotation run

  print(paste("Got: ",target,sep=""))
  print(paste("From project: ",project.name,sep=""))
  


##########################
  ## indel.ori<-indels
  ## load("indels.RData")
  ## tapply(indel.ori[,"FILTER"],indel.ori[,"FILTER"],length)
  ## key<-build.key(indels,core.ann)
  ## key.ori<-build.key(indel.ori,core.ann)
  ## sum(key!=key.ori)

########################################### SET AFFECTION STATUS #############################
########################################### SET AFFECTION STATUS #############################
########################################### SET AFFECTION STATUS #############################

if(update.gene.names){
ann.order<-colnames(geneanno.table)[grepl("::gene",colnames(geneanno.table))]
gene.names<-get.gene.names(geneanno.table,ann.order)
}

  
## attributes(gene.names)
## gene.names$"ensGene::gene"
  
##### use these when calling pheno analysis  needed
## setwd(annotate.dir)
##  load(paste(project.name,".table.RData",sep=""))
                                        # load data from the annotation run

 maf.threshold.filter<-maf.threshold.filter.to.use
 interesting.mutations<-interesting.mutations.use

if(!grepl("^chr",indels[1,"chr"])){
key.indels<-build.key(indels,core.ann,add.chr.label=TRUE)
}else{key.indels<-build.key(indels,core.ann)      }
## indels<-as.data.frame(All.grs) ### convert genomic ranges to a dta frame
## colnames(indels)[colnames(indels)=="seqnames"]<-"chr"
## key.indels<-paste(indels[,"chr"],indels[,"start"],indels[,"end"],indels[,"REF"],indels[,"ALT"],indels[,"TYPE"],sep=":")
## key<-paste("chr",rownames(geneanno.table),sep="")

if(dim(geneanno.table)[1]!=length(key.indels)){
  print("ERROR ERROR genenano and indels matrics different sizes")
  
  if(!grepl("^chr",rownames(geneanno.table)[1])){
    key<-paste("chr",rownames(geneanno.table),sep="")
  }else{key<-rownames(geneanno.table)     }

posns<-match(key.indels,key)
missing<-is.na(posns)
sum(missing)
indels<-indels[!missing,]
#All.grs<-All.grs[!missing,]
key.indels<-build.key(indels,core.ann,add.chr.label=TRUE)
sum(key!=key.indels) # should be zero
  } #dim(geneanno.table)[1]!=length(key.indels)

###################filters to apply
print(filter.cols.maf)
print(filter.cols.novel)
###################filters to apply



#############################get the filters so have minor allele frequencies:



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


maf.target<-the.combined.DBs[names(the.combined.DBs) %in% filter.cols.maf] # the.combined.DBs

  

for(i in 1:length(maf.threshold.filter)){
a.filter<-paste(filter.cols.maf,"::maf-filter::",maf.threshold.filter[i],sep="")
 assign(paste("filter.cols.maf",maf.threshold.filter[i],sep="."),value=a.filter)
}

maf.threshold.filter.all<-c(0.0,maf.threshold.filter)
assign(paste("filter.cols.maf",maf.threshold.filter.all[1],sep="."),value=filter.cols.novel)


  


################################################################################
## THE FILTER TABLE IS MADE FROM : the.combined.DBs<-c(filter.DB,generic.filter.DB,function.filter.DB)  only thes DBs can be used for MAF filtereing
target.table<-filter.table
#target.table[1:2,]

####### generate maf.target appreded to the allel frequency
maf.threshold.filter  # filters to apply

maf.target  #<-the.combined.DBs[names(the.combined.DBs) %in% all.filter.cols.maf] # the.combined.DBs

## ref.inversion<-matrix(data=TRUE,nrow=dim(target.table)[1],ncol=length(maf.threshold.filter))
## colnames(ref.inversion)<-maf.threshold.filter

### ADD COLUMNS LIKE: "hg19_snp137_maf.txt::maf-filter::0.05" that contain BOOLEAN to the FILTER TABLE "< filter" | > (1-filter) PASS:
#the.gen<-expand.labels.to.samples(c("GT"),samples.order.in.ALL) ### use this to predefine space
## the.new.filter.labels<-expand.labels.to.samples(maf.threshold.filter,maf.target)
## filter.table.extra<-matrix(data=TRUE,nrow=dim(target.table)[1],ncol=length(the.new.filter.labels))
## target.table<-cbind(target.table,filter.table.extra)



        
  

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

#rm(filter.table.extra)
#if(sum(ref.inversion[,1])>0){print("WARNING ref and ALT inversions noted")}
############ FIX ALLELE FREQUENCY PROBLEMS
#################################### ENS MAF 
############# Convert to database "names"

colnames(target.table)

the.DBs<-the.combined.DBs
the.DBs
colnames(the.DBs)
for(i in 1:length(the.DBs)){
  
target.string<-paste("^",the.DBs[i],sep="")
if(grepl("++",target.string,fixed=TRUE)){target.string<-gsub("++","",target.string,fixed=TRUE)} ##gerp++ casles problems in grep (++)
colnames(target.table)<-gsub(target.string,names(the.DBs)[i],colnames(target.table))

## colnames(target.table)<-gsub(paste("^",the.DBs[i],sep=""),names(the.DBs)[i],colnames(target.table))
}
colnames(target.table)

#target.table[1:2,]

filter.table<-target.table
########


  ################################# Get allele frequency table
################# HOWEVER ALT ALLELES ALL IS USED AS MINOR ALLELE FREQUNCY  not alternative allele frequencies
maf.lt.all<-data.frame(key=key.indels,stringsAsFactors=FALSE)

  
  
  for(imaf in 1:length(maf.threshold.filter.all)){
  a.filter.cols.maf<-eval(as.name( paste("filter.cols.maf",maf.threshold.filter.all[imaf],sep=".")  ))
  maf.lt<-rep(TRUE,times=dim(filter.table)[1])
  for(imaff in 1:length(a.filter.cols.maf)){
    if(maf.threshold.filter.all[imaf]==0){ # different case for novel NOT found (rather than less than)
       maf.lt<-maf.lt & !filter.table[,a.filter.cols.maf[imaff]]
     }else{
    maf.lt<-maf.lt & filter.table[,a.filter.cols.maf[imaff]]
  }
  }
# filtered<-maf.lt & wanted.muts.fil
 maf.lt.all<-cbind(maf.lt.all,maf.lt)
}
  

  if(dim(maf.lt.all)[2]>1){ colnames(maf.lt.all)<-c("key",paste("MAF.lt:",maf.threshold.filter.all,sep=""))}

 maf.lt.all<- maf.lt.all[,colnames(maf.lt.all)!="key"]



############ FIX ALLELE FREQUENCY PROBLEMS
#################################### ENS MAF 
############# Convert to database "names"
## i<-3
## sum(target.table[,paste(maf.target[minor.allele.DB],"maf",sep="::")]> 0.6 & !is.na(target.table[,paste(maf.target[minor.allele.DB],"maf",sep="::")]))

######################################## FIX MINOR ALLELE PROBLEMS IN ALT ALLELES summary.geno.group #####################33
################# MAFs > 0.05 or max(are included when doing filtering using ---- maf.threshold.filter[k]) | (1-maf.threshold.filter[k]) ) --- SO FILTERING ROBUST TO REF FLIPS
################# HOWEVER ALT ALLELES ALL IS USED AS MINOR ALLELE FREQUNCY  not alternative allele frequencies
target.table<-filter.table

################################not sure what to do about allel flips yet..
 skip.flip<-TRUE
  if(!skip.flip){
## not.mafs<-sum(as.numeric(target.table[,paste(minor.allele.DB,"maf",sep="::")])> 0.5 & is.finite(as.numeric(target.table[,paste(minor.allele.DB,"maf",sep="::")])) & !is.na(as.numeric(target.table[,paste(minor.allele.DB,"maf",sep="::")]))  )
## if(not.mafs > 0){
##   print(paste("WARNING minor.allele.DB: ",minor.allele.DB," contains ",not.mafs ," snps with MAF>0.5",sep=""))
## }

##### local mag > 0.05 and not NAN , not too may missing 1/2=0.5, is a otherwise rare SNP not an non infomrtive with ALT.Alleles would become zero on flip
min.counts.missing<- max((0.25*length(all.sample.labels)*2),10)
freq.flip<-(as.numeric(summary.geno.group[,"MAF.ALL"]) > 0.5 & is.finite(as.numeric(summary.geno.group[,"MAF.ALL"])) & !is.na(as.numeric(summary.geno.group[,"MAF.ALL"])) &
                as.numeric(summary.geno.group[,"MISSING.Alleles.ALL"]) < min.counts.missing & maf.lt.all[, dim(maf.lt.all)[2]] & as.numeric(summary.geno.group[,"REF.Alleles.ALL"]) !=0  ) # MAF > 0.5 and not NaN
summary.geno.group<-ref.inversion.for.summary(freq.flip,summary.geno.group)  ### make ALT allele the minor allele

}
  freq.flip<-rep(FALSE,times=dim(indels)[1])

### some of these are just sampling error exclude
#  sum(freq.flip)
## cbind(indels[,c("FILTER","ID")],summary.geno.group,maf.lt.all,filter.table)[freq.flip,][7,]
#cbind(indels[,c("FILTER","ID")],summary.geno.group2,maf.lt.all)[freq.flip,][1:5,]
##  common.snps<-maf.lt.all[, dim(maf.lt.all)[2]]
  
## common.snps<-as.numeric(target.table[,paste(minor.allele.DB,"maf",sep="::")])>= max(maf.threshold.filter) & is.finite(as.numeric(target.table[,paste(minor.allele.DB,"maf",sep="::")])) & !is.na(as.numeric(target.table[,paste(minor.allele.DB,"maf",sep="::")]))  

## ref.inverions that  are common snps (ones with maf > max(maf.threshold.filter) will be ingored and the ALT alleles left as is (the will become invisable using ALT alleles all filtering
## sum(ref.inversion & common.snps)
#ref.inversion[common.snps]<-FALSE

#summary.geno.group<-ref.inversion.for.summary(freq.flip,summary.geno.group)  ### make ALT allele the minor allele



##############testing#################
## sum(ref.inversion)
## sum(common.snps)
## cbind(summary.geno.group,target.table[,paste(minor.allele.DB,"maf",sep="::")],ref.inversion,common.snps)[ref.inversion,][1:10,]
## cbind(summary.geno.group.ori,target.table[,paste(minor.allele.DB,"maf",sep="::")],ref.inversion,common.snps)[ref.inversion,][1:10,]
## summary.geno.group[28:32,]
## summary.depths.group[1:2,]
## summary.het.indiv[1:2,] summary.geno.group.ori<-summary.geno.group

######################## Construct table of gene lists whene gene names were matched
######################## Here OMIM table shoudl always be present to gene.table is never empty

gene.lists<-c("AML_OTHER_Mutations","AML_MITO_Mutations","AML_LEY_RELAPSE_Mutations","AML_FRANC_Mutations","AML_ENU_Mutations","AML_DINDEL_Relapse","AML_SNP_Mutations","AML_DINDEL_Mutations","AML_CLINICAL_Mutations","AML_ASH_Mutations","AML_ALL_Mutations","omim","ALS","skeletome","mouse.defect","sewell.cycling","lung.genes","Dequeant.cycling","ingenuity.bone.genes","HypoMag","ProtonPump")

#  gene.lists<-{}
  
current.objects<-ls()
gene.lists<-gene.lists[gene.lists %in% current.objects]
gene.table<-{}
  if(length(gene.lists)!=0){
for(it in 1:length(gene.lists)){
  a.temp<-eval(as.name(gene.lists[it]))
  if(is.null(dim(gene.table)) ){gene.table<-a.temp}else{ if( dim(a.temp)[1]==dim(gene.table)[1] ){gene.table<-cbind(gene.table,a.temp)} }
                                  }
rm("a.temp")
if(dim(gene.table)[1]!=dim(indels)[1]){print("ERROR gene.table DIFFERENT number of rows to indels-> making empty genetable");gene.table<-as.matrix(key.indels) ; colnames(gene.table)<-"key"}
dim(gene.table)
}else{ gene.table<-as.matrix(key.indels) ; colnames(gene.table)<-"key"}
  
colnames(gene.table)
#####################################################
## test<-""
## chk<-cbind(indels[,c(core.ann,"ID")],gene.table,geneanno.table)


  ###########
## output summary of the whole lot here
##toString(ls())
## ###########
############# UPDATE PHOLYPHEN

  
functional.data.file<-paste(gsub(".txt$","",target),".FunctionalPredictions.RData",sep="")
if(functional.data.file %in% files & !force.functional.update){
    setwd(annotate.dir) # or where ever this file is saved
    print(paste("Load Functional Data:",functional.data.file))
    load(functional.data.file)
  }else{
    setwd(code.dir)
    print("REDOing functional annotation")
    source("pholyphen_read_and_Update.r")
    save(list=c("filter.table.pholy","pholy.data","sift.data","regulation.data"),file=paste(gsub(".txt$","",target),".FunctionalPredictions.RData",sep=""))
     }
  print(paste("GOT Functional Data:",functional.data.file))
 ## dim(filter.table.pholy)
# filter.table.pholy[1:5,]
## print(tapply(filter.table.pholy[,"PolyPhen.desc"],filter.table.pholy[,"PolyPhen.desc"],length))
## print(tapply(filter.table.pholy[,"SIFT.desc"],filter.table.pholy[,"SIFT.desc"],length))
## print(tapply(filter.table.pholy[,"regulation.desc"],filter.table.pholy[,"regulation.desc"],length))
## print(sort(tapply(filter.table.pholy[,"Consequence"],filter.table.pholy[,"Consequence"],length)))
## print(sort(tapply(filter.table.pholy[,"Feature"],filter.table.pholy[,"Feature"],length)))
## test<-"inframe_deletion"
## test<-grep(test,poly[,"Consequence"])
## length(test)
## poly[test,][1:20,]


  ##################################### merge VEP and annovar ENSEMBLE annotations
## clean.unique.combine<-function(x,y){
##   z<-unique(c(x,y))
##   z[!is.na(z) & !grepl("^PLACE_",z)] #z[z!="-" & !is.na(z) & !grepl("^PLACE_",z)]
## }

  ### Merge TWO LISTs columnwise using mapply:
  ## sort(tapply(filter.table.pholy[,"Gene"],filter.table.pholy[,"Gene"],length),decreasing=TRUE)[1:10]

  ### this may help in gene counts esp in MT which are not annotated with annovar.
  if(update.annotations.with.vep){
    vep.has.no.value<-is.na(filter.table.pholy[,"Gene"]) | filter.table.pholy[,"Gene"]=="-"  #intergenic_variant in vep Consequence not labelled 
    vep.has.no.feature<-is.na(filter.table.pholy[,"regulation.feature"]) | filter.table.pholy[,"regulation.feature"]=="-"

    
    ## sum(vep.has.no.value)
    ## sum(vep.has.no.feature)
    ## sum(vep.has.no.value & vep.has.no.feature)
    ## sum(vep.has.no.value & !vep.has.no.feature)
    ## filter.table.pholy[(vep.has.no.value & !vep.has.no.feature),][1:20,] regulation.feature
    ## filter.table.pholy[(!vep.has.no.feature),][1:200,]
    
    new.values.gene<-filter.table.pholy[(!vep.has.no.value | !vep.has.no.feature) ,"Gene"]
    new.values.gene[new.values.gene=="-"]<-NA
    new.values.feature<-filter.table.pholy[(!vep.has.no.value | !vep.has.no.feature) ,"regulation.feature"]
    new.values<-mapply(c, new.values.gene, new.values.feature, SIMPLIFY=FALSE)
    ## new.values[(vep.has.no.value & !vep.has.no.feature)][1:20]

     
  ## testing  
  ##   sum(!vep.has.no.value)
  ## attributes(gene.names)
  ## t1<-gene.names$"ensGene::gene"[!vep.has.no.value][95:100]
  ## t2<-as.list(filter.table.pholy[!vep.has.no.value,"Gene"][95:100])
  ## t3<-mapply(c, gene.names$"ensGene::gene"[180:200], as.list(filter.table.pholy[,"Gene"][180:200]), SIMPLIFY=FALSE) 
  ## vep.annotations<-mapply(clean.unique.combine, gene.names$"ensGene::gene"[!vep.has.no.value], as.list(filter.table.pholy[!vep.has.no.value,"Gene"]), SIMPLIFY=FALSE)
  ##  vep.annotations<-mapply(clean.unique.combine, gene.names$"ensGene::gene"[(!vep.has.no.value | !vep.has.no.feature)], new.values, SIMPLIFY=FALSE)                                                                                                                                
  ## gene.names$"ensGene::gene"[(!vep.has.no.value | !vep.has.no.feature)]<-vep.annotations  # gene.names$"ensGene::gene"[(!vep.has.no.value | !vep.has.no.feature)][1:10]

               
   vep.annotations<-mapply(clean.unique.combine, gene.names$"ensGene::gene"[grep(TRUE,(!vep.has.no.value | !vep.has.no.feature))], new.values, SIMPLIFY=FALSE)
    
  gene.names$"ensGene::gene"[grep(TRUE,(!vep.has.no.value | !vep.has.no.feature))]<-vep.annotations 


      ## fix missing annovar annottaion for ensGene
    missing.annovar.anno<-is.na(geneanno.table[!vep.has.no.value,"ensGene::gene"])
    geneanno.table[!vep.has.no.value,][missing.annovar.anno,"ensGene::gene"] <- filter.table.pholy[!vep.has.no.value,][missing.annovar.anno,"Gene"]

    
  ##### Uncomment below if also want to change the geneanno.table  
  ## test<-unlist(collapse.gene.names(vep.annotations,1,delimit.by=","))
  ## geneanno.table[(!vep.has.no.value | !vep.has.no.feature),"ensGene::gene"]<-paste(geneanno.table[(!vep.has.no.value | !vep.has.no.feature),"ensGene::gene"],test,sep=",")

#################################
    ## update GENE ANNOTATIONS
##     require("biomaRt")
## if(genome.build=="hg18"){mart<-useMart("ensembl_mart_51",host="may2009.archive.ensembl.org",dataset="hsapiens_gene_ensembl",archive=TRUE)
##     }else if(genome.build=="hg19"){
##       mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
##     }else if(genome.build=="mm9"){
##       mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl");print("here:")
##  }
## genes<-unique(unlist(gene.names[["ensGene::gene"]]))
## genes<-genes[!grepl("NONE",genes) | !grepl("^PLACE_",genes) | !grepl("^NO_ANNOVAR_",genes) ]  #### NONE(dist=NONE exists test<-!grepl("NONE",genes) & !grepl("^PLACE_",genes)


    
## if(genome.build=="hg19"){
## gene.ann.wanted<-c("ensembl_gene_id","hgnc_symbol","gene_biotype","description")
## }

## if(genome.build=="hg18"){
## gene.ann.wanted<-c("ensembl_gene_id","hgnc_symbol","biotype","description")
## }
  
## if(genome.build=="mm9"){
## gene.ann.wanted<-c("ensembl_gene_id","mgi_symbol","gene_biotype","description")
## }
## if(length(genes)>0){ 
## gene.desc<-getBM(attributes = gene.ann.wanted, filters = "ensembl_gene_id", values=genes, mart = mart)
## }else{gene.desc<-matrix(data=NA,nrow=1,ncol=length(gene.ann.wanted)); colnames(gene.desc)<-gene.ann.wanted}

    


    
## ## gene.desc[1,]
## ######### get rid of dupicate entries with resplce to the fist column  dim(gene.desc)[1]==length(unique(gene.desc[,1]))
##   dups<-duplicated(gene.desc[,1])
##   gene.desc<-gene.desc[!dups,]
## ####################### match genes to another column #########################################
## gene.desc.table<-match.cols.and.collapse(gene.names,"ensGene::gene",gene.desc,"ensembl_gene_id",gene.ann.wanted,"delimit")
## gene.desc.table[1:5,]
########################################################################################

load("/mnt/UQCCG/Software/annovar/humandb/Biomart.gene.desc.RData")
  dups<-duplicated(gene.desc[,1])
  gene.desc<-gene.desc[!dups,]
####################### match genes to another column #########################################

gene.desc.table<-match.cols.and.collapse(gene.names,"ensGene::gene",gene.desc,"Ensembl.Gene.ID",colnames(gene.desc),"delimit")

## the above is much faster then below when doing delimited concatination compated to tapply
## list.element.lengths<-unlist(lapply(gene.names[["ensGene::gene"]],length))        
## indel.index<-rep(1:length(list.element.lengths),times=list.element.lengths)
## flat.gene.list<-unlist(the.genes[[annotation.labels[i]]])  
## gene.desc.table<-match.cols.and.collapse.fast(list.element.lengths,indel.index,flat.gene.list,gene.desc,"ensembl_gene_id",gene.ann.wanted,"delimit")

print("Gene desc table")
#gene.desc.table[1:5,]


    
  } ## end update.annotations.with.vep
## ####################
## get rid of maf > threshold here as not wanted

################## filter.table.pholy[wanted.regulation,][1:20,]
 if(exists("poly") & class(poly)!="function"){rm("poly")} ## poly is also a function

########################### GERP UPDATE

  setwd(annotate.dir)
  gerp.data.file<-paste(gsub(".txt$","",target),".GerpPredictions.RData",sep="")
  if(force.GERP.update | !(gerp.data.file %in% files) ){
    gerp.scores<-get.GERP.MULT(db.details,indels[,c("chr","start","end","REF")])
    save(list=c("gerp.scores"),file=paste(gsub(".txt$","",target),".GerpPredictions.RData",sep=""))
  }else{
    print(paste("Load Gerp:",gerp.data.file))
    load(gerp.data.file)
    gerp.keys<-build.key(indels,c("chr","start","end","REF"))  
    if( sum(names(gerp.scores)!=gerp.keys)>0 | length(gerp.scores)!=length(gerp.keys) ){print("Gerp not matching indels recalculating")
                                                                                       gerp.scores<-get.GERP.MULT(db.details,indels[,c("chr","start","end","REF")])
                                                                                      save(list=c("gerp.scores"),file=paste(gsub(".txt$","",target),".GerpPredictions.RData",sep=""))
                                                                                      } ## gerp scores out of alignment    
  }
 print(paste("GOT Gerp:",gerp.data.file))
## missing<-is.na(filter.table[,"GERP::score"])
## sum(missing)
## sum(filter.table[!missing,"GERP::score"] !=gerp.scores[!missing])
## gerp.element.threshold<- 1 # any hit is significant
## gerp.score.threshold<-3 # gerp score >= will be included

  if("GERP::score" %in% colnames(filter.table)){
  filter.table[,"GERP::score"] <-gerp.scores
    }else{
  filter.table<-cbind(filter.table,gerp.scores)
  colnames(filter.table)[colnames(filter.table)=="gerp.scores"]<-"GERP::score"
    }
  
  gerp.hit<- gerp.scores >= gerp.score.threshold & !is.na(gerp.scores)
  gerp.hit.low<- (gerp.scores >= gerp.score.threshold.low & !is.na(gerp.scores)) |  (gerp.scores ==0 & !is.na(gerp.scores))
  gerp.unknown<-gerp.scores ==0
  
  ## sum(gerp.hit)
  ## sum(gerp.hit.low)
  ## sum(gerp.unknown)
##   sum(indel.type & !is.na(Gerp.scores))
############ indels are are a problem casue they have not been annotated would need a Run length encode object that has max GERP over the Indel span
############ Here have use the gerp elemets that are annoted with the regionanno.DB  (constrained regions) now sure if that beeter than doing the above?
  ## cbind(indels[,core.ann],Gerp.scores,regionanno.table[,"gerpelem::chrom, chromStart, chromEnd, score"])[!is.na(gerp.scores),][1:20,]
## this indicates that a treshold on the Gerp elemnts might not make any sense.
## think the first method is prefered as some indels are only 2-3 bp and are much smaller that the gerp elemnet length which I think is too course
# about 3/4 for gerp element hits are not gerp hits also  tends to support that gerp elments are too score and use scored directly
#low.ALT.counts<-get.alleles.per.region(indels,data.gr,genos.low[,"ALT"],pass) ## sometink like thisrequired over elements

### Retire the below
## gerp.elements<-strsplit(regionanno.table[,"gerpelem::chrom, chromStart, chromEnd, score"],split=":")
## gerp.elements<-as.numeric(unlist(lapply(gerp.elements,function(x) x[length(x)])))
## gerp.element.hit<-gerp.elements >= gerp.element.threshold & !is.na(gerp.elements)
## sum(gerp.element.hit) # about 3/4 for gerp element hits are not gerp hits think I should not use these


##  vep.coding<-c("not_assigned","stop_gained","stop_lost","missense_variant","splice_acceptor_variant","splice_donor_variant","splice_region_variant","initiator_codon_variant","stop_retained_variant")
## vep.noncoding<-c("mature_miRNA_variant","non_coding_exon_variant","TF_binding_site_variant","5_prime_UTR_variant","3_prime_UTR_variant","regulatory_region_variant")              
## vep.unwanted<-c("nc_transcript_variant","synonymous_variant","upstream_gene_variant","downstream_gene_variant","intergenic_variant","intron_variant","NMD_transcript_variant")



  
################################### complete wanted and types of mutaions up front ########

the.types<-names(tapply(geneanno.table[,"ensGene::location"],geneanno.table[,"ensGene::location"],length))
the.types<-unique(unlist(strsplit(the.types,split=";")))
if( sum( !(the.types %in% c("NA","exonic",possible.mutations)))>0  ){print("WARNING ANNOVAR HAS NEW MUTATION TYPES DEFINED- REVIEW line 532")
                                                                print( the.types[!(the.types %in% c("exonic",possible.mutations))] )   }

#  filter.table.pholy[1:5,]
wanted.muts.coding<-test.for.coding.type(geneanno.table,geneanno.DB,interesting.coding.mutations)

   ####################### mhairi CHECK
wanted.muts.coding.vep<-test.wanted.mutation(filter.table.pholy[,"Consequence"],vep.coding,delimit.by=",")  # filter.table.pholy[,"Consequence"] %in%  vep.coding
  ## filter.table.pholy[wanted.muts.coding.vep,]

  
  ## sum(wanted.muts.coding.vep)
  ##  sum(wanted.muts.coding)
############## noncoding double filtered for biotype ###

interesting.NONcoding.mutations<-interesting.mutations[!(interesting.mutations %in% interesting.coding.mutations)] # "ncRNA_exonic"
wanted.muts.NONcoding<-test.for.coding.type(geneanno.table,geneanno.DB,interesting.NONcoding.mutations)

wanted.interesting.to.prefilter<-test.for.coding.type(geneanno.table,geneanno.DB,interesting.to.prefilter) #"UTR3"      "UTR5"      "UTR5;UTR3" "snoRNA"


  ####################### mhairi CHECK
wanted.interesting.to.prefilter.vep<-test.wanted.mutation(filter.table.pholy[,"Consequence"],vep.noncoding,delimit.by=",")  #filter.table.pholy[,"Consequence"] %in% vep.noncoding

  ## filter.table.pholy[wanted.interesting.to.prefilter.vep,]
## sum(wanted.interesting.to.prefilter.vep)  
## sum(wanted.interesting.to.prefilter) # 124532

## wanted.regulation<-!is.na(filter.table.pholy[,"regulation.desc"])
## sum(wanted.regulation)
## sum(wanted.regulation & !wanted.interesting.to.prefilter.vep)
## filter.table.pholy[wanted.regulation & !wanted.interesting.to.prefilter.vep,]

indel.widths<-as.integer(indels[,"end"])-as.integer(indels[,"start"]) ## most are zero as start at same location
if(sum( indel.widths<0 | is.na(indel.widths))>0){print("ERROR indels width negative OR NA")}

  
###############################################################  
## Filter interesting.NONcoding.mutations for Biotypes to exclude all but wanted.noncoding.subtypes this removes Snow RNA's
print("Biotypes")
## print(toString(names(tapply(gene.desc.table[wanted.muts.NONcoding,"gene_biotype"],gene.desc.table[wanted.muts.NONcoding,"gene_biotype"],length))))
 ## "3prime_overlapping_ncrna, antisense, antisense::antisense, antisense::miRNA, antisense::pseudogene, lincRNA, lincRNA::antisense, lincRNA::lincRNA, lincRNA::processed_transcript, lincRNA::pseudogene, lincRNA::rRNA, lincRNA::snRNA, miRNA, misc_RNA, processed_transcript, processed_transcript::3prime_overlapping_ncrna, processed_transcript::antisense, processed_transcript::pseudogene, protein_coding, protein_coding::antisense, protein_coding::lincRNA, protein_coding::protein_coding, protein_coding::protein_coding::protein_coding, protein_coding::pseudogene, protein_coding::snRNA, pseudogene, pseudogene::antisense, pseudogene::lincRNA, pseudogene::protein_coding, pseudogene::pseudogene, pseudogene::snRNA, rRNA, sense_intronic, sense_overlapping, sense_overlapping::sense_intronic, snoRNA, snRNA, snRNA::lincRNA, snRNA::pseudogene"

wanted.muts.NONcoding.keep<-rep(FALSE,times=dim(geneanno.table)[1]) # "ncRNA_exonic" and "ncRNA_exonic" 
a.type<-wanted.noncoding.subtypes # list biotypes want to keep
#itype <- 1
for(itype in 1:length(a.type)){
  print(colnames(gene.desc.table))
  the.test<-grepl(a.type[itype],gene.desc.table[,"Gene.Biotype"])
  the.test[is.na(the.test)]<-FALSE
  wanted.muts.NONcoding.keep<- wanted.muts.NONcoding.keep | the.test
}
print("done miRNA")
 #  HERE wanted.muts.NONcoding.keep JUST DENOTES  "miRNA" &  "lincRNA" at this point USE wanted.muts.NONcoding to restrict to exones and splice BELOW
  wanted.muts.NONcoding.keep<-wanted.muts.NONcoding.keep & wanted.muts.NONcoding

  

  
##   sum(gerp.hit)
##   sum(gerp.hit.low)
##   sum(gerp.unknown)
## sum(wanted.muts.coding)
## sum(wanted.muts.coding.vep)
## sum(wanted.muts.NONcoding.keep)
## sum((wanted.interesting.to.prefilter | wanted.interesting.to.prefilter.vep))
## sum(  (wanted.interesting.to.prefilter | wanted.interesting.to.prefilter.vep) & gerp.hit )

## wanted.muts.print<-wanted.muts.coding | wanted.muts.coding.vep | (wanted.muts.NONcoding.keep ) | (  (wanted.interesting.to.prefilter | wanted.interesting.to.prefilter.vep)  )
wanted.muts<-wanted.muts.coding | wanted.muts.coding.vep | (wanted.muts.NONcoding.keep & (gerp.hit.low | gerp.unknown) ) | ( (wanted.interesting.to.prefilter | wanted.interesting.to.prefilter.vep) & (gerp.hit) )
  
#  | gerp.hit | gerp.unknown


### used for MS Whole genome run
## wanted.muts<-wanted.muts.coding | wanted.muts.coding.vep | (wanted.muts.NONcoding.keep & gerp.hit.low) | (  (wanted.interesting.to.prefilter | wanted.interesting.to.prefilter.vep) & gerp.hit ) 



wanted.muts.coding <- wanted.muts.coding | wanted.muts.coding.vep

  
  sum(wanted.muts)
  prefilter<-(wanted.muts |  (gerp.hit.low | gerp.unknown)) & maf.lt.all[,dim(maf.lt.all)[2]] #### what will be printed out
  sum(prefilter) # 206839 (maf filter) goes to 88021 for 50 exomes wanted | gerp.hit=2 | gerp.unknown
#
# maf.lt.all[,dim(maf.lt.all)[2]]


  ########## below will remove those that do not pass from analsysis
Prefilter.all.data<-FALSE

if(Prefilter.all.data){
  prefilter<-wanted.muts & maf.lt.all[,dim(maf.lt.all)[2]] ### wanted and Rare
  sum(prefilter)
 ## current.objects<-ls()
 ##  io<-11
 ##  test<-current.objects
  indels.size<-dim(indels)[1]
  for(io in 1:length(current.objects)){

    ## if(class(eval(as.name(current.objects[io]))) %in% c("IRanges","Rle","RleViews","array")){
    ##   current.objects[io]
    ##   print(class(eval(as.name(current.objects[io]))))
    ##   print(paste(io,current.objects[io],class(eval(as.name(current.objects[io]))),sep="-> "))
    ## }
    test[io]<-  class(eval(as.name(current.objects[io])))[1]
 
    if(class(eval(as.name(current.objects[io])))[1] %in% c("matrix","data.frame","array")){
       if(dim(eval(as.name(current.objects[io])))[1]==indels.size){
         print(paste(io,current.objects[io],class(eval(as.name(current.objects[io])))[1],sep="-> "))
         temp<-eval(as.name(current.objects[io]))
         assign(current.objects[io],value=temp[prefilter,])
    }}

      if(class(eval(as.name(current.objects[io])))[1] %in% c("numeric","logical","integer","character")){
       if(length(eval(as.name(current.objects[io])))==indels.size){
         print(paste(io,current.objects[io],class(eval(as.name(current.objects[io]))),sep="-> "))
         temp<-eval(as.name(current.objects[io]))
         assign(current.objects[io],value=temp[prefilter])
    }}

  }
}


  

#### use gene.desc.table here instead? may have missing va
Gene.Names<-collapse.gene.names(gene.names,"refGene::gene",delimit.by="::")
Gene.Names[grepl("^PLACE_",Gene.Names)]<-NA
#####################################################


  ## quality.names[!(quality.names %in% colnames(indels))]
  ## colnam  - 19 >SKDP-36.4    SKDP-80.1       SKDP-80.2       SKDP-80.202     SKDP-80.3       SKDP-80.4       UK10K_CIL5165136
###################################################################################################
###################################################  RESTART RESTART RESTART  ################################################################
###################################################  RESTART RESTART RESTART  ################################################################
###################################################  RESTART RESTART RESTART  ################################################################
###################################################  RESTART RESTART RESTART  ################################################################
###################################################  RESTART RESTART RESTART  ################################################################
###################################################  RESTART RESTART RESTART  ################################################################
###################################################  RESTART RESTART RESTART  ################################################################
###################################################  RESTART RESTART RESTART  ################################################################
if(length(all.sample.labels)>0){ # all.sample.labels==0 means it is a plink.assoc or bim file annotation run

## if(sum(affection.status!="sample_sheet")>0){
## posns<-match(all.sample.labels,affection.status)
## missing<-is.na(posns)
## posns
## sample.labels<-affection.status[posns[!missing]]
## names(sample.labels)<-names(affection.status)[posns[!missing]]
## family.name.prefix <- project
## }else{
  
if(exists("sample.labels")){rm("sample.labels")}



sheet.dir<-dirname(the.sample.sheet)
the.sample.sheet<-basename(the.sample.sheet)
if(sheet.dir=="."){sheet.dir<-bam.dir} # IF no path specified the assume it is BAM dir

fam.format<-FALSE
if(grepl(".fam$",the.sample.sheet)){fam.format<-TRUE}
the.sample.sheet

if(fam.format){
  sample.sheet.full<-read.delim(paste(sheet.dir,the.sample.sheet,sep="/"),header=F,sep="\t",fill=TRUE,stringsAsFactors=FALSE) # read a fam file
  sample.sheet.full<-cbind(project.name,sample.sheet.full)
  colnames(sample.sheet.full)<-c("SampleProject","FamilyCode","ParticipantCode","PaternalID","MaternalID","Sex","AffectionStatus")
 
}else{
sample.sheet.full<-read.delim(paste(sheet.dir,the.sample.sheet,sep="/"),header=T,sep=",",fill=TRUE,stringsAsFactors=FALSE)
}



if(!("ParticipantCode" %in% colnames(sample.sheet.full))){ # somtime it's called "ParticipantCode" or "Participant Code"
  print("here")
  if(sum(grepl("^Participant",colnames(sample.sheet.full) ))==1){colnames(sample.sheet.full)[grepl("^Participant",colnames(sample.sheet.full) )]<-"ParticipantCode"}
}

if(!("SampleProject" %in% colnames(sample.sheet.full))){ # somtime it's called "ParticipantCode" or "Participant Code"
  print("Fam")
  if(sum(grepl("^Sample Project",colnames(sample.sheet.full) ))==1){colnames(sample.sheet.full)[grepl("^Sample Project",colnames(sample.sheet.full) )]<-"SampleProject"}
}

if(!("SampleProject" %in% colnames(sample.sheet.full))){ # somtime it's called "ParticipantCode" or "Participant Code"
  print("or")
  if(sum(grepl("^Sample.Project",colnames(sample.sheet.full) ))==1){colnames(sample.sheet.full)[grepl("^Sample.Project",colnames(sample.sheet.full) )]<-"SampleProject"}
}

## } # got the sample sheet


if(sum(!(sample.sheet.full[,"ParticipantCode"] %in% all.sample.labels  ))!=0){
  print(paste("ERROR: REMOVING: samples not annotated by sample sheet: ", toString(sample.sheet.full[,"ParticipantCode"][!(sample.sheet.full[,"ParticipantCode"] %in% all.sample.labels  )])
              ))
  ## print("Reducing sample sheet")
  ## all.sample.labels[!(all.sample.labels %in% sample.sheet.full[,"ParticipantCode"])]
  
}
sample.sheet.full<-sample.sheet.full[(sample.sheet.full[,"ParticipantCode"] %in% all.sample.labels  ),]


if(sum(sample.sheet.full[,"ParticipantCode"] %in% all.sample.labels)!=length(sample.sheet.full[,"ParticipantCode"])){print("WARNING: sample sheet contains samples not in sequence data ")}


family.name.prefix.vector<-sample.sheet.full[,"SampleProject"]
family.name.prefix.vector<-unique(family.name.prefix.vector)
family.name.prefix.vector       #family.name.prefix.vector.ori !=family.name.prefix.vector  family.name.prefix.vector<-family.name.prefix.vector[-7]
family.name.prefix.vector<-family.name.prefix.vector[family.name.prefix.vector!=""]
#sample.sheet.full[,"ParticipantCode"]
maf.threshold.filter.all

if(project=="2013-02-21_GroupCallSKDP_Run75_119_WDR60_ContaminatedPlates"){family.name.prefix.vector<-c("SKDP-MELO","SKDP-FAM-57","SKDP-FAM-77","SKDP-FOP")}



if(exists("restrict.family") & !is.null(restrict.family)){print(paste("Restrcting analysis to :",restrict.family,sep=""));family.name.prefix.vector<-family.name.prefix.vector[family.name.prefix.vector %in% restrict.family]} ##restrict to just one family
if(exists("exclude.family") & !is.null(restrict.family)){print(paste("Excluding from  analysis:",restrict.family,sep=""));family.name.prefix.vector<-family.name.prefix.vector[!(family.name.prefix.vector %in% exclude.family)]} 

family.name.prefix.vector
# ifam<-1  family.name.prefix.vector<-family.name.prefix.vector[c(8,14)]
for(ifam in 1:length(family.name.prefix.vector)){

family.name.prefix<-family.name.prefix.vector[ifam]
print(paste("Doing",family.name.prefix))

################ If !force.redo then will skip analysis that are already complete:

########################################################################################
########################################################################################
########################################################################################

if(exists("family.name.prefix")){
  #file.out.name<-paste(project,family.name.prefix,sep=".")
  file.out.name<-paste(gsub(".txt$","",target),family.name.prefix,sep=".")
}else{file.out.name<-gsub(".txt$","",target)}

paste(file.out.name,"analysis.txt",sep=".")

if(!force.redo){
  
if(file.exists( paste(analysis.dir,paste(file.out.name,"analysis-maf-filtered.txt",sep="."),sep="/") ) ){
  print(paste( paste(file.out.name,"analysis-maf-filtered.txt",sep="."), "-> exists - SKIPPING"))
  next
}
}

########################################################################################
########################################################################################
########################################################################################


wanted.samples<-{}
wanted.samples<-c(wanted.samples,grep(paste("^",family.name.prefix,"$",sep=""),sample.sheet.full[,"SampleProject"]))

wanted.samples<-unique(wanted.samples)
sample.sheet<-sample.sheet.full[wanted.samples,]
sample.labels<-sample.sheet[,"ParticipantCode"]
affection.status<-sample.sheet[,"AffectionStatus"]
### WARNING WARNING CHANGE
affection.status[is.na(affection.status)]<-2  ### no affection statues then set as affected !!
###
affection.status<-gsub("2","Aff",affection.status);affection.status<-gsub("1","UnAff",affection.status)
affection.status<-gsub("0","unknown",affection.status);affection.status<-gsub("9","unknown",affection.status)
sample.labels<-gsub("#",".",sample.labels)
#sample.labels<-gsub("-",".",sample.labels)
names(sample.labels)<-affection.status

duplicate.samples<-duplicated(sample.labels)
if(sum(duplicate.samples)>0){ # just in case sample sheet as duplicate enties for the same sample
  print("WARNING SAMPLE sheet has duplicate entries for the sample- for this project ID")
  sample.labels<-sample.labels[!duplicate.samples] 
}





##############WARNING THIS IS WRONG FOR 2 AFFECTED IN COMPOUND HER AND OTHER BEDS######
if(exists("the.mother") | exists("the.father")){rm("the.mother","the.father")}
   if(sum(c("PaternalID","MaternalID") %in% colnames(sample.sheet))==2 ){ # maternal and paternal id found
   
     the.mother<-unique(sample.sheet[,"MaternalID"])
     the.father<-unique(sample.sheet[,"PaternalID"])

     the.mother<-gsub("#",".",the.mother);the.mother<-gsub("-",".",the.mother)
     the.father<-gsub("#",".",the.father);the.father<-gsub("-",".",the.father)

     the.mother<-the.mother[!(the.mother=="x" | the.mother=="" | the.mother=="0" | is.na(the.mother))]
     the.father<-the.father[!(the.father=="x" | the.father==""| the.father=="0" | is.na(the.father))]

     
     if(sum(the.mother %in% the.father)>0){print("Parents mixed up in PED mother=father")}
     if( length(the.mother)>1 | length(the.father)>1 | sum(c(the.mother,the.father) %in% sample.labels)<2){print("To many parents do not model Compound Het");rm("the.mother","the.father")}
}

 
print(sample.labels)
## controls<-c("BalbC","C3H","C57BL6","DBA")
## sample.labels<- controls

#####################



######################

the.affected<-sample.labels[names(sample.labels)=="Aff"]
the.Unaffected<-sample.labels[names(sample.labels)=="UnAff"]
the.unknown<-sample.labels[names(sample.labels)=="unknown"]

num.affected<-length(sample.labels[names(sample.labels)=="Aff"])
num.Unaffected<-length(sample.labels[names(sample.labels)=="UnAff"])
num.unknown<-length(sample.labels[names(sample.labels)=="unknown"])


##############get other alt alleles.cases

########################### Genotypes summary for just the 
## the.samples<-samples.order.in.ALL[samples.order.in.ALL %in% all.bm.samples] all.sample.labels

#paste(the.samples,"GT",sep=".")[!(paste(the.samples,"GT",sep=".") %in% colnames(indels))]
the.samples<-c(the.affected,the.Unaffected,the.unknown)

genotypes<-indels[,paste(the.samples,"GT",sep=".")]
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO")
rownames(summary.geno)<-key.indels
#summary.geno<-ref.inversion.for.summary(freq.flip,summary.geno)


################# DEFINE OTHER samples

## if(!exists("alt.allele.other.use")){alt.allele.other.use<-FALSE} ## in alt all counts use a specifc all.alleles threhold !caseful with "include.control.in.other<-TRUE if sample in 2 peds!
## ################ CAN RUN IN THREE MODES
## if(!exists("include.control.in.other")){include.control.in.other<-FALSE} # any sample NOT in current family AND item in sampleproject as CONTROL or control set to other group
## if(!exists("control.ONLY.in.other")){control.ONLY.in.other<-FALSE} # item in sampleproject as CONTROL or control set to other group other fams NOT USED
## if(!exists("make.nonDefined.samples.as.OTHER")){make.nonDefined.samples.as.OTHER<-FALSE}  ## other group called samples not listed in the sample sheet become   CONTROLS automatically
## #######################

the.samples.other<-all.sample.labels[!(all.sample.labels %in% the.samples)]
if(include.control.in.other){
  the.samples.other<-unique(c(the.samples.other,sample.sheet.full[sample.sheet.full[,"SampleProject"]=="CONTROL","ParticipantCode"],
                              sample.sheet.full[sample.sheet.full[,"SampleProject"]=="control","ParticipantCode"],
                                sample.sheet.full[sample.sheet.full[,"SampleProject"]=="Control","ParticipantCode"],
                              all.sample.labels[!(all.sample.labels %in% sample.sheet.full[,"ParticipantCode"])]
                              ))
}

if(control.ONLY.in.other){
  the.samples.other<-unique(c(sample.sheet.full[sample.sheet.full[,"SampleProject"]=="CONTROL","ParticipantCode"],sample.sheet.full[sample.sheet.full[,"SampleProject"]=="control","ParticipantCode"]))
  }

if(make.nonDefined.samples.as.OTHER){
the.samples.other<-unique(all.sample.labels[!(all.sample.labels %in% sample.sheet.full[,"ParticipantCode"])])
  }

################################



if(length(the.samples.other)==0){summary.geno.other<-matrix(data=NA,nrow=dim(indels)[1],ncol=dim(summary.geno)[2]) }else{
genotypes<-indels[,paste(the.samples.other,"GT",sep=".")]
summary.geno.other<-genotype.summary(as.matrix(genotypes))
} ## ens else for not other groups
colnames(summary.geno.other)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),"OTHER",sep=".")
rownames(summary.geno.other)<-key.indels
summary.geno.other<-ref.inversion.for.summary(freq.flip,summary.geno.other)

#  cbind(indels[,c("FILTER","ID")],summary.geno,maf.lt.all)[freq.flip,][1:5,]
genotypes<-indels[,paste(the.samples,"GT",sep=".")]
allele.depths<-indels[,paste(the.samples,"AD",sep=".")]
allele.depths[genotypes!="0/1"] <- NA
summary.depths<-apply(as.matrix(allele.depths),1,allele.summary.with.apply)
summary.depths<-t(summary.depths)
colnames(summary.depths)<-c("Hetero.ALT.reads","Hetero.REF.reads","Hetero.Read.Balance")
rownames(summary.depths)<-key.indels

## summary.depths.info<-to.info.format(summary.depths)
## summary.geno.info<-to.info.format(summary.geno)
## summary.filter.info<-to.info.format(indels[,"FILTER"],labels=c("FILTER"))
## summary.geno.info[1:5]
## summary.depths.info[1:5]
## summary.filter.info[1:5]

## full.info.all<-paste(summary.depths.info,summary.geno.info,summary.filter.info,indels[,"INFO"],sep=";")
## the.af.all<-summary.geno[,"MAF"]

###########################





########################################### SET AFFECTION STATUS #############################
########################################### SET AFFECTION STATUS #############################

if(!Skip.Analysis){
 
######################################################################################################################
######################################################################################################################
### Do filtering
######################################################################################################################
######################################################################################################################


############ require DP information in the analysis for each genotype GET from AN
if( sum(paste(sample.labels,"DP",sep=".") %in% colnames(indels))==length(sample.labels)){
  allele.depths<-as.matrix(indels[,paste(sample.labels,"DP",sep=".")])
  print("DP information present")
   }else{
     print("DP information NOT present")
  allele.depths.group<-indels[,paste(sample.labels,"AD",sep=".")]
  allele.depths<-apply(as.matrix(allele.depths.group),2,allele.DP.from.AN)
  colnames(allele.depths)<-gsub(".AD",".DP",colnames(allele.depths))
  indels<-cbind(indels,allele.depths)
}
###### set up the coverage and quality filter

setwd(code.dir)
source("run.filtering.r")
setwd(annotate.dir)

dim(quality.thresh)
######################################################################################################################
######################################################################################################################



###############################################################################################################################################################
####################################################### Finished quality calculation ###################################################################
###############################################################################################################################################################
######################define good quality (passing all global characterists)
#c("QUAL","QD","HRun","SB","FILTER_PASS","FILTER_100")



gatk.filter<-c("FILTER_PASS")
good.qual.filter<-apply( as.matrix(quality.thresh[,gatk.filter]) ,1,sum)
good.qual.filter<-good.qual.filter==length(gatk.filter)
########################################################3
tapply(indels[,"FILTER"],indels[,"FILTER"],length)

if(!use.only.GATK.pass){
traditional<-secondary.filter[secondary.filter %in% colnames(quality.thresh)] # c("QUAL","QD","HRun","SB")
good.qual<- apply( quality.thresh[,traditional] ,1,sum)
good.qual<-good.qual==length(traditional)
good.qual<-good.qual | good.qual.filter
}else{
good.qual<-good.qual.filter
}
sum(good.qual)
sum(good.qual.filter)

## good.qual<-good.qual.filter # this just filter =pass
########### model dominance - same indel in all affected not in unaffected


global.qual<-global.quality.labs

## dominant.alleles<-c("GT_AB","GT_B")
## dominant.coverage<-c("DP_High","DP_Low")


dominant.alleles<-c("GT_AB","GT_B")
dominant.coverage<-c("GT_AB","DP_Thresh") # use DP_thresh when what to include 1/1 only when coverage is LOW 

ref.alleles<-c("GT_A")
ref.corverage<-c("DP_Low")

if(num.affected==0){dom.aff.GT=c();dom.aff.GT.DP=c()}else{
dom.aff.GT<-expand.labels.to.samples(dominant.alleles,sample.labels[names(sample.labels)=="Aff"])
dom.aff.GT.DP<-expand.labels.to.samples(dominant.coverage,sample.labels[names(sample.labels)=="Aff"]) ##GT.DP cause match genotype to DP 0/1 high coverage than 1/1
}

if(num.Unaffected==0){dom.Unaff.GT=c();dom.Unaff.GT.DP=c()}else{
dom.Unaff.GT<-expand.labels.to.samples(ref.alleles,sample.labels[names(sample.labels)=="UnAff"]) ###must use direct not alternative
dom.Unaff.GT.DP<-expand.labels.to.samples(ref.corverage,sample.labels[names(sample.labels)=="UnAff"])
}


#### careful direct not alternative cause faulein coverge will make it appear to be the other genotope 0/0 (zero count = 1/1 Dp=2)
if(num.affected==0){dom.aff.ALT.counts<-rep(0,times=dim(indels)[1])}else{
dom.aff.ALT.counts<-  apply(as.matrix(quality.thresh[,dom.aff.GT] & quality.thresh[,dom.aff.GT.DP]) ,1,sum)
}# number affected with  0/1 and 1/1 for each indels matrix & matrix works fine (per element)
#dom.aff.ALT.counts.fil[1:5]

if(num.Unaffected==0){dom.Unaff.REF.counts<-rep(0,times=dim(indels)[1])}else{
dom.Unaff.REF.counts<-apply( as.matrix(quality.thresh[,dom.Unaff.GT] & quality.thresh[,dom.Unaff.GT.DP]),1,sum) # number affected with  0/1 and 1/1 for each indels
}
                      

###### Fit dominate profile no coverage or quality
## dom.aff.ALT.counts<-apply(as.matrix(quality.thresh[,dom.aff.GT]),1,sum) # number affected with  0/1 and 1/1 for each indels # need as .matrix incase there is just one sample-> becomes a vector 
## dom.Unaff.REF.counts<-apply(as.matrix(quality.thresh[,c(dom.Unaff.GT)]),1,sum) # number unaffected with  0/0 for each indels

#tapply(dom.aff.ALT.counts,dom.aff.ALT.counts,length)
#tapply(dom.Unaff.REF.counts,dom.Unaff.REF.counts,length)

dom.profile<-(dom.aff.ALT.counts==num.affected) & (dom.Unaff.REF.counts==num.Unaffected)  ### dom profile

##################################################################################################################################################################

###### Fit dominate profile WITH coverage

dominant.alleles<-c("GT_AB","GT_B")
dominant.coverage<-c("DP_High","DP_Thresh") # use DP_thresh when what to include 1/1 only when coverage is LOW 

ref.alleles<-c("GT_A")
ref.corverage<-c("DP_Low")

if(num.affected==0){dom.aff.GT=c();dom.aff.GT.DP=c()}else{
dom.aff.GT<-expand.labels.to.samples(dominant.alleles,sample.labels[names(sample.labels)=="Aff"])
dom.aff.GT.DP<-expand.labels.to.samples(dominant.coverage,sample.labels[names(sample.labels)=="Aff"]) ##GT.DP cause match genotype to DP 0/1 high coverage than 1/1
}

if(num.Unaffected==0){dom.Unaff.GT=c();dom.Unaff.GT.DP=c()}else{
dom.Unaff.GT<-expand.labels.to.samples(ref.alleles,sample.labels[names(sample.labels)=="UnAff"]) ###must use direct not alternative
dom.Unaff.GT.DP<-expand.labels.to.samples(ref.corverage,sample.labels[names(sample.labels)=="UnAff"])
}



#### careful direct not alternative cause faulein coverge will make it appear to be the other genotope 0/0 (zero count = 1/1 Dp=2)
if(num.affected==0){dom.aff.ALT.counts.fil<-rep(0,times=dim(quality.thresh)[1])}else{
dom.aff.ALT.counts.fil<-  apply( as.matrix(quality.thresh[,dom.aff.GT] & quality.thresh[,dom.aff.GT.DP]) ,1,sum)
}# number affected with  0/1 and 1/1 for each indels matrix & matrix works fine (per element)
#dom.aff.ALT.counts.fil[1:5]

if(num.Unaffected==0){dom.Unaff.REF.counts.fil<-rep(0,times=dim(quality.thresh)[1])}else{
dom.Unaff.REF.counts.fil<-apply( as.matrix(quality.thresh[,dom.Unaff.GT] & quality.thresh[,dom.Unaff.GT.DP]),1,sum) # number affected with  0/1 and 1/1 for each indels
}
                      
## tapply(dom.aff.ALT.counts.fil,dom.aff.ALT.counts.fil,length)
## tapply(dom.Unaff.REF.counts.fil,dom.Unaff.REF.counts.fil,length)

dom.profile.fil<-(dom.aff.ALT.counts.fil==num.affected) & (dom.Unaff.REF.counts.fil==num.Unaffected)  ### dom profile with coverage

#dom.profile.fil.qual<-dom.profile.fil & good.qual  ### dom profile with coverage and quality


############################### IN GROUP filter options #################################
## summary.geno.group[1:5,]
## summary.depths.group[1:5,]


#ok.in.group<-summary.geno.group[,"ALT.Alleles.ALL"] <= (num.affected+3) # this is what I typically use
ok.in.group<-as.numeric(summary.geno.group[,"ALT.Alleles.ALL"]) <= (num.affected)*10 # in case is fully recessive
on.xychromo<-indels[,"chr"]=="chrX" | indels[,"chr"]=="chrY"

sum(ok.in.group)
sum(on.xychromo)
extra.counts.for.sex<-as.numeric(summary.geno.group[,"ALT.Alleles.ALL"]) <= (num.affected+3)*10
sum(on.xychromo & extra.counts.for.sex)
ok.in.group<-ok.in.group | (on.xychromo & extra.counts.for.sex)

############### other filter for gene counts:
if(alt.allele.other.use){
  alt.alleles.other.test<-(as.numeric(summary.geno.other[,"ALT.Alleles.OTHER"]) <= alt.allele.other.use.threshold*length(the.samples.other)) | is.na(summary.geno.other[,"ALT.Alleles.OTHER"])
 # sum(alt.alleles.other.test)
  ok.in.group<-ok.in.group &  alt.alleles.other.test
}
sum(ok.in.group)
not.flat.genotype<-!quality.thresh[,"flat"]
#################################################3
## ok.balance.in.group<-((summary.depths.group[,"Hetero.Read.Balance.ALL"] >= 20) & (summary.depths.group[,"Hetero.Read.Balance.ALL"] <= 80)) | is.na(summary.depths.group[,"Hetero.Read.Balance.ALL"]) | (summary.depths.group[,"Hetero.ALT.reads.ALL"] <= num.affected*5) | grepl("^indel",indels[,"TYPE"])

## sum(grepl("^indel",indels[,"TYPE"]))
## sum(grepl("^snp",indels[,"TYPE"]))
## ok.balance.in.group<-((summary.depths.group[,"Hetero.Read.Balance.ALL"] >= 20) & (summary.depths.group[,"Hetero.Read.Balance.ALL"] <= 80)) | is.na(summary.depths.group[,"Hetero.Read.Balance.ALL"]) | (summary.depths.group[,"Hetero.ALT.reads.ALL"] <= num.affected*10) | indels[,"TYPE"]=="indel"
## test<-is.na(summary.depths.group[,"Hetero.Read.Balance.ALL"])
## sum(ok.balance.in.group)
#summary.depths.group[ok.balance.in.group,][1:10,]


######FIT PHENOTYPE WITH QUALITY AND MAF
wanted.muts.fil<-wanted.muts & good.qual ## boolean if want that indel and it was good coverage
#######################FIT Phenotype (recessive requires the parents)


####################### FIT COMPOUND HET POSSIBLES###################################
key.hetero<-build.key(indels,core.ann)
Compound.Hetero<-data.frame(key=key.hetero,stringsAsFactors=FALSE)
Denovo.DOM<-data.frame(key=key.hetero,stringsAsFactors=FALSE)

if(exists("the.mother") & exists("the.father")){

ref.the.mother<-expand.labels.to.samples(c("GT_A"),the.mother)
ref.the.father<-expand.labels.to.samples(c("GT_A"),the.father)

compound.het.FM<-quality.thresh[,ref.the.mother] & !quality.thresh[,ref.the.father]
compound.het.MF<-!quality.thresh[,ref.the.mother] & quality.thresh[,ref.the.father]

Denovo.DOM.test<-quality.thresh[,ref.the.mother] & quality.thresh[,ref.the.father]
Compound.Hetero.test<-compound.het.FM | compound.het.MF  # indels[Compound.Hetero,the.gen][1:50,37:40]

Compound.Hetero<-cbind(Compound.Hetero,Compound.Hetero.test)
Denovo.DOM<-cbind(Denovo.DOM,Denovo.DOM.test)
}
####################### FIT COMPOUND HET POSSIBLES###################################


the.affected<-sample.labels[names(sample.labels)=="Aff"]
the.Unaffected<-sample.labels[names(sample.labels)=="UnAff"]
the.unknown<-sample.labels[names(sample.labels)=="unknown"]

num.affected<-length(sample.labels[names(sample.labels)=="Aff"])
num.Unaffected<-length(sample.labels[names(sample.labels)=="UnAff"])
num.unknown<-length(sample.labels[names(sample.labels)=="unknown"])

#"ok.in.group" CONTAINS alt.allele.other.use
## use.the.globals<-c("wanted.muts","good.qual","ok.in.group","ok.balance.in.group")
use.the.globals<-c("wanted.muts","good.qual","ok.in.group","not.flat.genotype")
#----------FMHT in R stop here
##  a.filter.cols.maf %in% colnames(filter.table)
## paste(sample.labels,"GT",sep=".") %in% colnames(indels)
## colnames(indels)[grep("\\.GT$",colnames(indels))]
## Error in `colnames<-`(`*tmp*`, value = c("2009-048", "2009-165", "2008-308",  : 
##   length of 'dimnames' [2] not equal to array extent
## In addition: There were 50 or more warnings (use warnings() to see the first 50)
                                 
################# NOW TO PHENO NOVEL ( with 0% MAF )
########################################################################################
####################################################################################################################
####################################################################################################################
#
#a.filter.cols.maf[!(a.filter.cols.maf %in% colnames(filter.table))]
# imaf<-1

for(imaf in 1:length(maf.threshold.filter.all)){
  if(maf.threshold.filter.all[imaf]!=0){next} # novel done above need different filter actions
  a.filter.cols.maf<-eval(as.name( paste("filter.cols.maf",maf.threshold.filter.all[imaf],sep=".")  ))
  print(maf.threshold.filter.all[imaf])


  ################ Counting alternate allele counts  in 2 classes  - so ro ref to refernce alleles###### include.homozygotes<-TRUE
if(include.homozygotes){
the.attributes<-c("GT_AB","GT_B")
the.attributes.coverage<-c("GT_AB","GT_B") # use GT_AB so there is no threshold
}else{
the.attributes<-c("GT_AB","GT_B")
the.attributes.coverage<-c("GT_AB","DP_Thresh") # use GT_AB so there is no threshold inclde GT_B is few reads
}
  
  if(num.affected==0){dom.aff.GT=c()}else{
dom.aff.GT<-expand.labels.to.samples(the.attributes,sample.labels[names(sample.labels)=="Aff"])
}
if(num.Unaffected==0){dom.Unaff.GT=c()}else{
dom.Unaff.GT<-expand.labels.to.samples(the.attributes,sample.labels[names(sample.labels)=="UnAff"]) ###must use direct not alternative 
}
if(num.affected==0){dom.aff.GT.DP=c()}else{
dom.aff.GT.DP<-expand.labels.to.samples(the.attributes.coverage,sample.labels[names(sample.labels)=="Aff"]) ##GT.DP cause match genotype to DP 0/1 high coverage than 1/1
}
if(num.Unaffected==0){dom.Unaff.GT.DP=c()}else{
dom.Unaff.GT.DP<-expand.labels.to.samples(the.attributes.coverage,sample.labels[names(sample.labels)=="UnAff"])
}
#######################
  

if(length(dom.aff.GT.DP)==0 | length(dom.aff.GT)==0 ){input<-as.matrix(quality.thresh[,dom.aff.GT])}else{input<-as.matrix(quality.thresh[,dom.aff.GT] & quality.thresh[,dom.aff.GT.DP])} #can't multiply a null matrix
pheno.qual.aff<-CountOverGene.withFilterAndGlobal(input,the.affected,the.attributes,gene.names,filter.table,a.filter.cols.maf,"OR",filter.action="EXCLUDE",global.vars=use.the.globals )


if(length(dom.Unaff.GT.DP)==0 | length(dom.Unaff.GT)==0 ){input<-as.matrix(quality.thresh[,dom.Unaff.GT])}else{input<-as.matrix(quality.thresh[,dom.Unaff.GT] & quality.thresh[,dom.Unaff.GT.DP])} #can't multiply a null matrix                       
pheno.qual.Unaff<-CountOverGene.withFilterAndGlobal(input,the.Unaffected,the.attributes,gene.names,filter.table,a.filter.cols.maf,"OR",filter.action="EXCLUDE",global.vars=use.the.globals )
  

##$$
pheno.extension<-paste(":",maf.threshold.filter.all[imaf],".maf.qual",sep="")
the.list<-pheno.qual.aff
the.list.Unaff<-pheno.qual.Unaff
##$$

## names(the.list.Unaff)
## the.list$hit[1:10]
## the.list$hit.count[1:5]
## dim(the.list$gene.counts)
## the.list.Unaff$gene.counts[1:5,]
## # CHECK

colnames(the.list$gene.counts)<-paste(colnames(the.list$gene.counts),pheno.extension,sep="")
colnames(the.list.Unaff$gene.counts)<-paste(colnames(the.list.Unaff$gene.counts),pheno.extension,sep="")
assign(paste("phenotype.gene.hit",pheno.extension,sep=""),value=the.list$hit)
assign(paste("phenotype.gene.hit.count",pheno.extension,sep=""),value=the.list$hit.count)

assign(paste("phenotype.gene.hit.Unaff",pheno.extension,sep=""),value=the.list.Unaff$hit)
assign(paste("phenotype.gene.hit.count.Unaff",pheno.extension,sep=""),value=the.list.Unaff$hit.count)
  
assign(paste("gene.counts",pheno.extension,sep=""),value=the.list$gene.counts)
assign(paste("gene.counts.Unaff",pheno.extension,sep=""),value=the.list.Unaff$gene.counts)
####################################################################################################################
####################################################################################################################
## wanted<-geneanno.table[,"refGene::gene"]=="MAFB"
## the.list$hit[wanted]
## the.list$hit.count[wanted]
## the.list$gene.counts[wanted,]
## the.list.Unaff$gene.counts[1:5,]

##   the.gen<-expand.labels.to.samples(c("GT"),all.sample.labels)
##   indels[wanted,c(core.ann,"FILTER",the.gen)]

##  indels[wanted,1:30]

######FIT PHENO NOVEL  WITH COVERAGE QUALITY 
########################################################################################
####################################################################################################################
####################################################################################################################
  ################ Counting alternate allele counts  in 2 classes  - so ro ref to refernce alleles#####
  
if(include.homozygotes){
the.attributes<-c("GT_AB","GT_B")
the.attributes.coverage<-c("DP_High","DP_Low") 
}else{  
the.attributes<-c("GT_AB","GT_B")
the.attributes.coverage<-c("DP_High","DP_Thresh") # use GT_AB so there is no threshold
}


  
if(num.affected==0){dom.aff.GT=c()}else{
dom.aff.GT<-expand.labels.to.samples(the.attributes,sample.labels[names(sample.labels)=="Aff"])
}
if(num.Unaffected==0){dom.Unaff.GT=c()}else{
dom.Unaff.GT<-expand.labels.to.samples(the.attributes,sample.labels[names(sample.labels)=="UnAff"]) ###must use direct not alternative 
}
if(num.affected==0){dom.aff.GT.DP=c()}else{
dom.aff.GT.DP<-expand.labels.to.samples(the.attributes.coverage,sample.labels[names(sample.labels)=="Aff"]) ##GT.DP cause match genotype to DP 0/1 high coverage than 1/1
}
if(num.Unaffected==0){dom.Unaff.GT.DP=c()}else{
dom.Unaff.GT.DP<-expand.labels.to.samples(the.attributes.coverage,sample.labels[names(sample.labels)=="UnAff"])
}

dom.aff.GT
dom.Unaff.GT
dom.aff.GT.DP
dom.Unaff.GT.DP
quality.thresh[1:5,dom.aff.GT]
quality.thresh[1:5,dom.Unaff.GT]
#######################
  
if(length(dom.aff.GT.DP)==0 | length(dom.aff.GT)==0 ){input<-as.matrix(quality.thresh[,dom.aff.GT])}else{input<-as.matrix(quality.thresh[,dom.aff.GT] & quality.thresh[,dom.aff.GT.DP])} #can't multiply a null matrix

pheno.qual.cov<-CountOverGene.withFilterAndGlobal(input,the.affected,the.attributes,gene.names,filter.table,a.filter.cols.maf,"OR",filter.action="EXCLUDE",global.vars=use.the.globals )

if(length(dom.Unaff.GT.DP)==0 | length(dom.Unaff.GT)==0 ){input<-as.matrix(quality.thresh[,dom.Unaff.GT])}else{input<-as.matrix(quality.thresh[,dom.Unaff.GT] & quality.thresh[,dom.Unaff.GT.DP])} #can't multiply a null matrix                       
pheno.qual.cov.Unaff<-CountOverGene.withFilterAndGlobal(input,the.Unaffected,the.attributes,gene.names,filter.table,a.filter.cols.maf,"OR",filter.action="EXCLUDE",global.vars=use.the.globals )


##$$
pheno.extension<-paste(":",maf.threshold.filter.all[imaf],".maf.qual.cov",sep="")
the.list<-pheno.qual.cov
the.list.Unaff<-pheno.qual.cov.Unaff
##$$


colnames(the.list$gene.counts)<-paste(colnames(the.list$gene.counts),pheno.extension,sep="")
colnames(the.list.Unaff$gene.counts)<-paste(colnames(the.list.Unaff$gene.counts),pheno.extension,sep="")
assign(paste("phenotype.gene.hit",pheno.extension,sep=""),value=the.list$hit)
assign(paste("phenotype.gene.hit.count",pheno.extension,sep=""),value=the.list$hit.count)
  
assign(paste("phenotype.gene.hit.Unaff",pheno.extension,sep=""),value=the.list.Unaff$hit)
assign(paste("phenotype.gene.hit.count.Unaff",pheno.extension,sep=""),value=the.list.Unaff$hit.count)
  
assign(paste("gene.counts",pheno.extension,sep=""),value=the.list$gene.counts)
assign(paste("gene.counts.Unaff",pheno.extension,sep=""),value=the.list.Unaff$gene.counts)
####################################################################################################################
####################################################################################################################
} # loop over maf that only does maf=0



####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
################# NOW TO PHENO with MAF



for(imaf in 1:length(maf.threshold.filter.all)){
  if(maf.threshold.filter.all[imaf]==0){next} # novel done above need different filter actions
  a.filter.cols.maf<-eval(as.name( paste("filter.cols.maf",maf.threshold.filter.all[imaf],sep=".")  ))
  print(maf.threshold.filter.all[imaf])
########################################################################################
  ######FIT PHENOTYPE WITH   MAF, QUALITY 
####################################################################################################################
####################################################################################################################

  ################ Counting alternate allele counts  in 2 classes  - so ro ref to refernce alleles######

if(include.homozygotes){
the.attributes<-c("GT_AB","GT_B")
the.attributes.coverage<-c("GT_AB","GT_B") # use GT_AB so there is no threshold
}else{
the.attributes<-c("GT_AB","GT_B")
the.attributes.coverage<-c("GT_AB","DP_Thresh") # use GT_AB so there is no threshold inclde GT_B is few reads
}
   # use GT_AB so there is no threshold
  
  if(num.affected==0){dom.aff.GT=c()}else{
dom.aff.GT<-expand.labels.to.samples(the.attributes,sample.labels[names(sample.labels)=="Aff"])
}
if(num.Unaffected==0){dom.Unaff.GT=c()}else{
dom.Unaff.GT<-expand.labels.to.samples(the.attributes,sample.labels[names(sample.labels)=="UnAff"]) ###must use direct not alternative 
}
if(num.affected==0){dom.aff.GT.DP=c()}else{
dom.aff.GT.DP<-expand.labels.to.samples(the.attributes.coverage,sample.labels[names(sample.labels)=="Aff"]) ##GT.DP cause match genotype to DP 0/1 high coverage than 1/1
}
if(num.Unaffected==0){dom.Unaff.GT.DP=c()}else{
dom.Unaff.GT.DP<-expand.labels.to.samples(the.attributes.coverage,sample.labels[names(sample.labels)=="UnAff"])
}
#######################


if(length(dom.aff.GT.DP)==0 | length(dom.aff.GT)==0 ){input<-as.matrix(quality.thresh[,dom.aff.GT])}else{input<-as.matrix(quality.thresh[,dom.aff.GT] & quality.thresh[,dom.aff.GT.DP])} # can't muliply to NULLS
pheno.maf.qual<-CountOverGene.withFilterAndGlobal(input,the.affected,the.attributes,gene.names,filter.table,a.filter.cols.maf,"AND",filter.action="INCLUDE",global.vars=use.the.globals )
  
if(length(dom.Unaff.GT.DP)==0 | length(dom.Unaff.GT)==0 ){input<-as.matrix(quality.thresh[,dom.Unaff.GT])}else{input<-as.matrix(quality.thresh[,dom.Unaff.GT] & quality.thresh[,dom.Unaff.GT.DP])} #can't multiply a null matrix
pheno.maf.qual.Unaff<-CountOverGene.withFilterAndGlobal(input,the.Unaffected,the.attributes,gene.names,filter.table,a.filter.cols.maf,"AND",filter.action="INCLUDE",global.vars=use.the.globals )




##$$
pheno.extension<-paste(":",maf.threshold.filter.all[imaf],".maf.qual",sep="")
the.list<-pheno.maf.qual
the.list.Unaff<-pheno.maf.qual.Unaff
##$$

colnames(the.list$gene.counts)<-paste(colnames(the.list$gene.counts),pheno.extension,sep="")
colnames(the.list.Unaff$gene.counts)<-paste(colnames(the.list.Unaff$gene.counts),pheno.extension,sep="")
assign(paste("phenotype.gene.hit",pheno.extension,sep=""),value=the.list$hit)
assign(paste("phenotype.gene.hit.count",pheno.extension,sep=""),value=the.list$hit.count)

assign(paste("phenotype.gene.hit.Unaff",pheno.extension,sep=""),value=the.list.Unaff$hit)
assign(paste("phenotype.gene.hit.count.Unaff",pheno.extension,sep=""),value=the.list.Unaff$hit.count)
  
assign(paste("gene.counts",pheno.extension,sep=""),value=the.list$gene.counts)
assign(paste("gene.counts.Unaff",pheno.extension,sep=""),value=the.list.Unaff$gene.counts)
####################################################################################################################
####################################################################################################################

######FIT PHENOTYPE WITH   MAF, COVERAGE, QUALITY 
########################################################################################
####################################################################################################################
####################################################################################################################

  ################ Counting alternate allele counts  in 2 classes  - so ro ref to refernce alleles######

if(include.homozygotes){
the.attributes<-c("GT_AB","GT_B")
the.attributes.coverage<-c("DP_High","DP_Low") 
}else{  
the.attributes<-c("GT_AB","GT_B")
the.attributes.coverage<-c("DP_High","DP_Thresh") # use GT_AB so there is no threshold
}
 # use GT_AB so there is no threshold
  
if(num.affected==0){dom.aff.GT=c()}else{
dom.aff.GT<-expand.labels.to.samples(the.attributes,sample.labels[names(sample.labels)=="Aff"])
}
if(num.Unaffected==0){dom.Unaff.GT=c()}else{
dom.Unaff.GT<-expand.labels.to.samples(the.attributes,sample.labels[names(sample.labels)=="UnAff"]) ###must use direct not alternative 
}
if(num.affected==0){dom.aff.GT.DP=c()}else{
dom.aff.GT.DP<-expand.labels.to.samples(the.attributes.coverage,sample.labels[names(sample.labels)=="Aff"]) ##GT.DP cause match genotype to DP 0/1 high coverage than 1/1
}
if(num.Unaffected==0){dom.Unaff.GT.DP=c()}else{
dom.Unaff.GT.DP<-expand.labels.to.samples(the.attributes.coverage,sample.labels[names(sample.labels)=="UnAff"])
}

dom.aff.GT
dom.Unaff.GT
dom.aff.GT.DP
dom.Unaff.GT.DP

#######################

## dom.aff.alt.counts.fil<-  apply( (quality.thresh[,dom.aff.GT] & quality.thresh[,dom.aff.GT.DP]) ,1,sum) # number affected with  0/1 and 1/1 for each indels matrix & matrix works fine (per element)
  ### coverage given by (quality.thresh[,dom.aff.GT]
if(length(dom.aff.GT.DP)==0 | length(dom.aff.GT)==0 ){input<-as.matrix(quality.thresh[,dom.aff.GT])}else{input<-as.matrix(quality.thresh[,dom.aff.GT] & quality.thresh[,dom.aff.GT.DP])} # can't muliply to NULLS

pheno.maf.qual.cov<-CountOverGene.withFilterAndGlobal(input,the.affected,the.attributes,gene.names,filter.table,a.filter.cols.maf,"AND",filter.action="INCLUDE",global.vars=use.the.globals )

if(length(dom.Unaff.GT.DP)==0 | length(dom.Unaff.GT)==0 ){input<-as.matrix(quality.thresh[,dom.Unaff.GT])}else{input<-as.matrix(quality.thresh[,dom.Unaff.GT] & quality.thresh[,dom.Unaff.GT.DP])} #can't multiply a null matrix

pheno.maf.qual.cov.Unaff<-CountOverGene.withFilterAndGlobal(input,the.Unaffected,the.attributes,gene.names,filter.table,a.filter.cols.maf,"AND",filter.action="INCLUDE",global.vars=use.the.globals )

##$$
pheno.extension<-paste(":",maf.threshold.filter.all[imaf],".maf.qual.cov",sep="")
the.list<-pheno.maf.qual.cov
the.list.Unaff<-pheno.maf.qual.cov.Unaff
##$$

colnames(the.list$gene.counts)<-paste(colnames(the.list$gene.counts),pheno.extension,sep="")
colnames(the.list.Unaff$gene.counts)<-paste(colnames(the.list.Unaff$gene.counts),pheno.extension,sep="")
assign(paste("phenotype.gene.hit",pheno.extension,sep=""),value=the.list$hit)
assign(paste("phenotype.gene.hit.count",pheno.extension,sep=""),value=the.list$hit.count)

assign(paste("phenotype.gene.hit.Unaff",pheno.extension,sep=""),value=the.list.Unaff$hit)
assign(paste("phenotype.gene.hit.count.Unaff",pheno.extension,sep=""),value=the.list.Unaff$hit.count)
  
assign(paste("gene.counts",pheno.extension,sep=""),value=the.list$gene.counts)
assign(paste("gene.counts.Unaff",pheno.extension,sep=""),value=the.list.Unaff$gene.counts)
####################################################################################################################
####################################################################################################################

} # loop over MAF cases

## pheno.qual$hit[1:10]
## pheno.qual$hit.count[1:10]
## pheno.qual$gene.counts[1:5,]
## ### phenotype::maf:0 qual gene identified and counts
## pheno.qual.cov$hit[1:10]
## pheno.qual.cov$hit.count[1:10]
## pheno.qual.cov$gene.counts[1:5,]
## ### phenotype::maf:5% qual gene identified and counts 
## pheno.maf.qual$hit[1:10]
## pheno.maf.qual$hit.count[1:10]
## pheno.maf.qual$gene.counts[1:5,] # good gene hits per individual
## ### phenotype::maf:5% qual coverage      gene identified and counts 
## pheno.maf.qual.cov$hit[1:10]
## pheno.maf.qual.cov$hit.count[1:10]
## pheno.maf.qual.cov$gene.counts[1:5,] # good gene hits per individual  gene.counts.maf.qual,gene.counts.maf.qual.cov,gene.counts.qual,gene.counts.qual.cov,
####################################################################################################################
####################################################################################################################




## dim(indels)
## length(key.indels)
## indels[1:2,]

## quality.thresh[1:2,]
## geneanno.table[1:2,]
## gene.desc.table[1:2,]
## regionanno.table[1:2,]
## filter.table[1:2,]

## wanted.muts[1:5] ## interesting mutations interesting.mutations
## wanted.muts.fil[1:5] ## interesting mutations  quality
## good.qual[1:5] ## QUAL and QC


## dom.aff.ALT.counts[1:5]
## dom.Unaff.REF.counts[1:5]
## dom.profile[1:5]   ## all affected 0/1 or 1/1
## dom.profile.fil[1:5] ## includes coverage high and low on varients
## dom.profile.fil.qual[1:5] ## includes coverage high and low on varients and quality

### phenotype::maf:0 qual gene identified and counts # phenotype.gene.hit.qual,phenotype.gene.hit.count.qual,phenotype.gene.hit.qual.cov,phenotype.gene.hit.count.qual.cov,


ls()[grepl("^pheno",ls())]
### phenotype::maf:0 qual gene identified and counts # phenotype.gene.hit.qual,phenotype.gene.hit.count.qual,phenotype.gene.hit.qual.cov,phenotype.gene.hit.count.qual.cov,



################# combine data for different MAF runs (including novel)
## ls()[grepl("gene.counts",ls())]
## grep("20:39317279:39317279:G:A:snp",key)

## filter.table[100746:100748,]
##  maf.lt.all[100746:100748,]

## Gene.Names<-collapse.gene.names(gene.names,"refGene::gene",delimit.by="::")

if(class(Gene.Names)=="list"){print("Error Gene.Names a list: L1338");Gene.Names<-unlist(lapply(Gene.Names,function(x) paste(x,collapse="::")))}

### key<-paste(indels[,"chr"],":",indels[,"start"],":",indels[,"end"],":",indels[,"REF"],":",indels[,"ALT"],":",indels[,"TYPE"],sep="")
key<-build.key(indels,core.ann)
## maf.lt.all<-data.frame(key=key,stringsAsFactors=FALSE)
gene.hit.count.fil.sample.all<-data.frame(key=key,stringsAsFactors=FALSE)
gene.hit.count.sample.all<-data.frame(key=key,stringsAsFactors=FALSE)
gene.hit.count.fil.all<-data.frame(key=key,stringsAsFactors=FALSE)
gene.hit.count.fil.aff.all<-data.frame(key=key,stringsAsFactors=FALSE)
gene.hit.count.fil.Unaff.all<-data.frame(key=key,stringsAsFactors=FALSE)

phenotype.gene.hit.count.maf.qual.all<-data.frame(key=key,stringsAsFactors=FALSE)
phenotype.gene.hit.count.maf.qual.cov.all<-data.frame(key=key,stringsAsFactors=FALSE)
phenotype.gene.hit.maf.qual.all<-data.frame(key=key,stringsAsFactors=FALSE)
phenotype.gene.hit.maf.qual.cov.all<-data.frame(key=key,stringsAsFactors=FALSE)

phenotype.gene.hit.count.maf.qual.all.Unaff<-data.frame(key=key,stringsAsFactors=FALSE)
phenotype.gene.hit.count.maf.qual.cov.all.Unaff<-data.frame(key=key,stringsAsFactors=FALSE)
phenotype.gene.hit.maf.qual.all.Unaff<-data.frame(key=key,stringsAsFactors=FALSE)
phenotype.gene.hit.maf.qual.cov.all.Unaff<-data.frame(key=key,stringsAsFactors=FALSE)


print("Do gene hit summary")

print(maf.threshold.filter.all)
for(imaf in 1:length(maf.threshold.filter.all)){
  
##   a.filter.cols.maf<-eval(as.name( paste("filter.cols.maf",maf.threshold.filter.all[imaf],sep=".")  ))
##   maf.lt<-rep(TRUE,times=dim(filter.table)[1])
  
##   for(i in 1:length(a.filter.cols.maf)){
##     if(maf.threshold.filter.all[imaf]==0){ # different case for novel NOT found (rather than less than)
##        maf.lt<-maf.lt & !filter.table[,a.filter.cols.maf[i]]
##      }else{
##     maf.lt<-maf.lt & filter.table[,a.filter.cols.maf[i]]
##   }
##   }
## }
  
## # filtered<-maf.lt & wanted.muts.fil
##  maf.lt.all<-cbind(maf.lt.all,maf.lt)


gene.counts.maf.qual<-eval(as.name(paste("gene.counts:",maf.threshold.filter.all[imaf],".maf.qual",sep="")))
gene.counts.Unaff.maf.qual<-eval(as.name(paste("gene.counts.Unaff:",maf.threshold.filter.all[imaf],".maf.qual",sep="")))
gene.hit.count.sample<-cbind(gene.counts.maf.qual,gene.counts.Unaff.maf.qual)
gene.hit.count.sample[1:5,]
gene.hit.count.sample<-gene.hit.count.sample[,paste(sample.labels[(names(sample.labels) %in% c("Aff","UnAff"))],":",maf.threshold.filter.all[imaf],".maf.qual",sep="")] # exclude unknows as have no counts
if(is.null(dim(gene.hit.count.sample))){ gene.hit.count.sample<-as.matrix(gene.hit.count.sample);colnames(gene.hit.count.sample)<-paste(sample.labels[(names(sample.labels) %in% c("Aff","UnAff"))],":",maf.threshold.filter.all[imaf],".maf.qual",sep="")} ##one affected and no unaffect the gene hit count is an array
                   
gene.counts.maf.qual.cov<-eval(as.name(paste("gene.counts:",maf.threshold.filter.all[imaf],".maf.qual.cov",sep="")))
gene.counts.Unaff.maf.qual.cov<-eval(as.name(paste("gene.counts.Unaff:",maf.threshold.filter.all[imaf],".maf.qual.cov",sep="")))
gene.hit.count.fil.sample<-cbind(gene.counts.maf.qual.cov,gene.counts.Unaff.maf.qual.cov)
gene.hit.count.fil.sample[1:5,]
gene.hit.count.fil.sample<-gene.hit.count.fil.sample[,paste(sample.labels[(names(sample.labels) %in% c("Aff","UnAff"))],":",maf.threshold.filter.all[imaf],".maf.qual.cov",sep="")] # clean up incase on aff or unaffected
if(is.null(dim(gene.hit.count.fil.sample))){ gene.hit.count.fil.sample<-as.matrix(gene.hit.count.fil.sample);colnames(gene.hit.count.fil.sample)<-paste(sample.labels[(names(sample.labels) %in% c("Aff","UnAff"))],":",maf.threshold.filter.all[imaf],".maf.qual.cov",sep="")} ##
# the.affected
gene.hit.count.fil<-apply(gene.hit.count.fil.sample,1,function(x) sum(x,na.rm=TRUE)) ## get total counts across all samples 
gene.hit.count.fil.all<-cbind(gene.hit.count.fil.all,gene.hit.count.fil)
  
if(length(the.affected)>0){gene.hit.count.fil.aff<-apply(as.matrix(gene.hit.count.fil.sample[,paste(the.affected,":",maf.threshold.filter.all[imaf],".maf.qual.cov",sep="")]),1,function(x) sum(x,na.rm=TRUE))
                           gene.hit.count.fil.aff.all<-cbind(gene.hit.count.fil.aff.all,gene.hit.count.fil.aff) }
if(length(the.Unaffected)>0){gene.hit.count.fil.Unaff<-apply(as.matrix(gene.hit.count.fil.sample[,paste(the.Unaffected,":",maf.threshold.filter.all[imaf],".maf.qual.cov",sep="")]),1,function(x) sum(x,na.rm=TRUE))
                            gene.hit.count.fil.Unaff.all<-cbind(gene.hit.count.fil.Unaff.all,gene.hit.count.fil.Unaff)}

 gene.hit.count.fil.sample.all<-cbind(gene.hit.count.fil.sample.all,gene.hit.count.fil.sample) # fil=qual and coveregae
 gene.hit.count.sample.all<-cbind(gene.hit.count.sample.all,gene.hit.count.sample) # fil=qual only
  
 phenotype.gene.hit.count.maf.qual<-eval(as.name(paste("phenotype.gene.hit.count:",maf.threshold.filter.all[imaf],".maf.qual",sep="")))
 phenotype.gene.hit.count.maf.qual.cov<-eval(as.name(paste("phenotype.gene.hit.count:",maf.threshold.filter.all[imaf],".maf.qual.cov",sep="")))
 phenotype.gene.hit.maf.qual<-eval(as.name(paste("phenotype.gene.hit:",maf.threshold.filter.all[imaf],".maf.qual",sep="")))
 phenotype.gene.hit.maf.qual.cov<-eval(as.name(paste("phenotype.gene.hit:",maf.threshold.filter.all[imaf],".maf.qual.cov",sep="")))


 phenotype.gene.hit.count.maf.qual.Unaff<-eval(as.name(paste("phenotype.gene.hit.count.Unaff:",maf.threshold.filter.all[imaf],".maf.qual",sep="")))
 phenotype.gene.hit.count.maf.qual.cov.Unaff<-eval(as.name(paste("phenotype.gene.hit.count.Unaff:",maf.threshold.filter.all[imaf],".maf.qual.cov",sep="")))
 phenotype.gene.hit.maf.qual.Unaff<-eval(as.name(paste("phenotype.gene.hit.Unaff:",maf.threshold.filter.all[imaf],".maf.qual",sep="")))
 phenotype.gene.hit.maf.qual.cov.Unaff<-eval(as.name(paste("phenotype.gene.hit.Unaff:",maf.threshold.filter.all[imaf],".maf.qual.cov",sep="")))


 phenotype.gene.hit.count.maf.qual.all<-cbind(phenotype.gene.hit.count.maf.qual.all,phenotype.gene.hit.count.maf.qual)
 phenotype.gene.hit.count.maf.qual.cov.all<-cbind(phenotype.gene.hit.count.maf.qual.cov.all,phenotype.gene.hit.count.maf.qual.cov)
 phenotype.gene.hit.maf.qual.all<-cbind(phenotype.gene.hit.maf.qual.all,phenotype.gene.hit.maf.qual)
 phenotype.gene.hit.maf.qual.cov.all<-cbind(phenotype.gene.hit.maf.qual.cov.all,phenotype.gene.hit.maf.qual.cov)

 phenotype.gene.hit.count.maf.qual.all.Unaff<-cbind(phenotype.gene.hit.count.maf.qual.all.Unaff,phenotype.gene.hit.count.maf.qual.Unaff)
 phenotype.gene.hit.count.maf.qual.cov.all.Unaff<-cbind(phenotype.gene.hit.count.maf.qual.cov.all.Unaff,phenotype.gene.hit.count.maf.qual.cov.Unaff)
 phenotype.gene.hit.maf.qual.all.Unaff<-cbind(phenotype.gene.hit.maf.qual.all.Unaff,phenotype.gene.hit.maf.qual.Unaff)
 phenotype.gene.hit.maf.qual.cov.all.Unaff<-cbind(phenotype.gene.hit.maf.qual.cov.all.Unaff,phenotype.gene.hit.maf.qual.cov.Unaff)


} #loop over imaf


############################################# make pretty column names 
## if(dim(maf.lt.all)[2]>1)                  {   colnames(maf.lt.all)<-c("key",paste("MAF.lt:",maf.threshold.filter.all,sep=""))}
if(dim(gene.hit.count.fil.all)[2]>1)       {  colnames(gene.hit.count.fil.all)<-c("key",paste("Gene.Hit.Count.QC.MAF:",maf.threshold.filter.all,sep=""))}
if(dim(gene.hit.count.fil.aff.all)[2]>1)  {  colnames(gene.hit.count.fil.aff.all)<-c("key",paste("Aff.Gene.Hit.Count.QC.MAF:",maf.threshold.filter.all,sep=""))}
if(dim(gene.hit.count.fil.Unaff.all)[2]>1){  colnames(gene.hit.count.fil.Unaff.all)<-c("key",paste("Unff.Gene.Hit.Count.QC.MAF:",maf.threshold.filter.all,sep=""))}

if(dim(phenotype.gene.hit.count.maf.qual.all)[2]>1) {  colnames(phenotype.gene.hit.count.maf.qual.all)<-c("key",paste("Shared.Gene.Aff.Count.Q.MAF:",maf.threshold.filter.all,sep=""))}
if(dim(phenotype.gene.hit.count.maf.qual.cov.all)[2]>1) {  colnames(phenotype.gene.hit.count.maf.qual.cov.all)<-c("key",paste("Shared.Gene.Aff.Count.QC.MAF:",maf.threshold.filter.all,sep=""))}
if(dim(phenotype.gene.hit.maf.qual.all)[2]>1)   {  colnames(phenotype.gene.hit.maf.qual.all)<-c("key",paste("Shared.Gene.Aff.Hit.Q.MAF:",maf.threshold.filter.all,sep=""))}
if(dim(phenotype.gene.hit.maf.qual.cov.all)[2]>1) {  colnames(phenotype.gene.hit.maf.qual.cov.all)<-c("key",paste("Shared.Gene.Aff.Hit.QC.MAF:",maf.threshold.filter.all,sep=""))}

if(dim(phenotype.gene.hit.count.maf.qual.all.Unaff)[2]>1) {  colnames(phenotype.gene.hit.count.maf.qual.all.Unaff)<-c("key",paste("Shared.Gene.UnAff.Count.Q.MAF:",maf.threshold.filter.all,sep=""))}
if(dim(phenotype.gene.hit.count.maf.qual.cov.all.Unaff)[2]>1) {  colnames(phenotype.gene.hit.count.maf.qual.cov.all.Unaff)<-c("key",paste("Shared.Gene.UnAff.Count.QC.MAF:",maf.threshold.filter.all,sep=""))}
if(dim(phenotype.gene.hit.maf.qual.all.Unaff)[2]>1)   {  colnames(phenotype.gene.hit.maf.qual.all.Unaff)<-c("key",paste("Shared.Gene.UnAff.Hit.Q.MAF:",maf.threshold.filter.all,sep=""))}
if(dim(phenotype.gene.hit.maf.qual.cov.all.Unaff)[2]>1) {  colnames(phenotype.gene.hit.maf.qual.cov.all.Unaff)<-c("key",paste("Shared.Gene.UnAff.Hit.QC.MAF:",maf.threshold.filter.all,sep=""))}

if(dim(gene.hit.count.fil.sample.all)[2]>1) {  colnames(gene.hit.count.fil.sample.all)<-gsub(".maf.qual.cov","(MAF).QC",colnames(gene.hit.count.fil.sample.all))}
if(dim(gene.hit.count.sample.all)[2]>1)  {  colnames(gene.hit.count.sample.all)<-gsub(".maf.qual","(MAF).Q",colnames(gene.hit.count.sample.all))}

## only "dom.profile"          "dom.profile.fil"      "dom.profile.fil.qual" exist but I want better names in the summary file
#dom.profile.Q<-dom.profile.qual 
dom.profile.C<-dom.profile.fil ## includes coverage high and low on varients
#dom.profile.QC<-dom.profile.fil.qual
#dom.Unaff.REF.counts
dom.aff.ALT.counts.C<-dom.aff.ALT.counts.fil
dom.Unaff.REF.counts.C<-dom.Unaff.REF.counts.fil

################### clean up pointless "exonic" lables in the gene annotations
for(i in 1:length(geneanno.DB)){
  geneanno.table[,paste(geneanno.DB[i],"::location",sep="")]<-gsub(";exonic","",geneanno.table[,paste(geneanno.DB[i],"::location",sep="")])
#  print(tapply(geneanno.table[,paste(geneanno.DB[i],"::location",sep="")],geneanno.table[,paste(geneanno.DB[i],"::location",sep="")],length))
}
#sum(wanted.muts)



##############################################################################################

## maf.lt.all[1:5,]
## gene.hit.count.fil.all[1:5,]
## gene.hit.count.fil.aff.all[1:5,]
## gene.hit.count.fil.Unaff.all[1:5,]
## phenotype.gene.hit.count.maf.qual.all[1:5,]
## phenotype.gene.hit.count.maf.qual.cov.all[1:5,]
## phenotype.gene.hit.maf.qual.all[1:5,]
## phenotype.gene.hit.maf.qual.cov.all[1:5,]
## gene.hit.count.fil.sample.all[1:5,]
## gene.hit.count.sample.all[1:5,]


#######################################################################################################
################################################ SHip.Analysis - Do not run analysis


} # end Skip.Analysis


###############################################

#########################  define order and columns wanted
#indels[1:5,core.ann]
wanted.qualities<-c("QUAL","FILTER","SB","FS","HRun")


annotation.in.filter.table<-c("Consequence",extra.vep.annotations)

colnames(filter.table.pholy)<-gsub(".Embl","",colnames(filter.table.pholy))
filter.table.annotation<-colnames(filter.table.pholy) %in% annotation.in.filter.table
colnames(filter.table.pholy)[filter.table.annotation]<-paste(colnames(filter.table.pholy)[filter.table.annotation],".Embl",sep="")


indels.order.GT.AN<-c(expand.labels.to.samples(c(the.affected,the.Unaffected,the.unknown),c("GT","AD"),paste.after=TRUE))
indels.order.SNP.QUAL<-wanted.qualities[wanted.qualities %in% colnames(indels)]


indels.order.1st<-c(indels.order.GT.AN,indels.order.SNP.QUAL)
filter.ID.cols<-colnames(indels) %in% c("ID","FILTER")


indels.order.2nd<-expand.labels.to.samples(c(the.affected,the.Unaffected,the.unknown),c("DP","GQ","PL"),paste.after=TRUE) ##"DP","GQ","PL" rarely used removed in unwanted.cols
indels.order.2nd<-indels.order.2nd[indels.order.2nd %in% colnames(indels)]
indels.order.other<-colnames(indels)[!(colnames(indels) %in% c(core.ann,indels.order.1st,indels.order.2nd))]


## indels[1:5,indels.order.1st]
print ("eep")
gene.desc.table.labels.wanted<-c("Description","Gene.Biotype") # gene.desc.table.labels.wanted<-c("description")
print ("eep gene.desc.table")
print(gene.desc.table)
print("gene.desc.table.labels.wanted")
print(gene.desc.table.labels.wanted)
gene.desc.table.labels.other<-colnames(gene.desc.table)[!( colnames(gene.desc.table) %in% gene.desc.table.labels.wanted )]

#filter.cols.maf.0

unwanted.columns<-c("refGene" ,"knownGene", "ensGene","gene.table",
                    "snp132::found","1000genome","snp132","CG69","EUR_ASN_AFR_INDEL","ALL_EXON_SNP","1000genomeV2.txt","RNSH_FMHT_maf.txt","ljb_pp2","avsift","ljb_mt","ljb_phylop",
                    colnames(filter.table)[grepl("::maf-filter",colnames(filter.table))],
                     colnames(regionanno.table)[grepl("::HIT",colnames(regionanno.table))],filter.cols.maf.0,
                    expand.labels.to.samples(c(the.affected,the.Unaffected,the.unknown),c("GQ","PL"),paste.after=TRUE),
                      names(the.combined.DBs) ## these are boolean filter.table columns 
                    )



## posn<-grep("50623682", indels[,"start"])

##      indels[posn,]
# skeletome ,mouse.defect,sewell.cycling,Dequeant.cycling,ingenuity.bone.genes,omim

## objects.to.delete.for.space<-c("poly","a.poly","quality.thresh")
## rm(list=c(objects.to.delete.for.space))

## dim(indels)
## geneanno.table
## Gene.Names
## gene.desc.table[,gene.desc.table.labels.wanted]


if(genome.build=="hg19" |  genome.build=="hg18"){
if(!Skip.Analysis){
summary<-cbind(indels[,c(core.ann)],geneanno.table,Gene.Names,gene.desc.table[,gene.desc.table.labels.wanted],gene.table,filter.table.pholy[,filter.table.annotation],indels[,filter.ID.cols],wanted.muts,wanted.muts.coding,good.qual,wanted.muts.fil,ok.in.group,maf.lt.all,Denovo.DOM,Compound.Hetero,gerp.scores,summary.geno,summary.geno.group,summary.geno.other,summary.depths,indels[,indels.order.GT.AN],indels[,indels.order.SNP.QUAL],filter.table.pholy[,c("PolyPhen.desc","PolyPhen.scores","SIFT.desc","SIFT.scores")],dom.profile,dom.aff.ALT.counts,dom.Unaff.REF.counts,dom.aff.ALT.counts.C,dom.Unaff.REF.counts.C,gene.hit.count.fil.all,gene.hit.count.fil.aff.all,gene.hit.count.fil.Unaff.all   ,phenotype.gene.hit.maf.qual.all,phenotype.gene.hit.maf.qual.all.Unaff, phenotype.gene.hit.count.maf.qual.all,phenotype.gene.hit.count.maf.qual.all.Unaff,phenotype.gene.hit.maf.qual.cov.all, phenotype.gene.hit.maf.qual.cov.all.Unaff ,gene.hit.count.sample.all,summary.depths.group,gene.desc.table[,gene.desc.table.labels.other],regionanno.table,indels[,indels.order.2nd],filter.table,indels[,indels.order.other])

}else{ ## Skip.Analysis so no gene information
  
summary<-cbind(indels[,c(core.ann)],geneanno.table,Gene.Names,gene.desc.table[,gene.desc.table.labels.wanted],gene.table,filter.table.pholy[,filter.table.annotation],indels[,filter.ID.cols],wanted.muts,wanted.muts.coding ,maf.lt.all,gerp.scores,summary.geno,summary.geno.group,summary.geno.other,summary.depths,indels[,indels.order.GT.AN],indels[,indels.order.SNP.QUAL],filter.table.pholy[,c("PolyPhen.desc","PolyPhen.scores","SIFT.desc","SIFT.scores")],summary.depths.group,gene.desc.table[,gene.desc.table.labels.other],regionanno.table,indels[,indels.order.2nd],filter.table,indels[,indels.order.other])

}  ## Skip.Analysis so no gene information


## summary<-cbind(indels[prefilter,c(core.ann)],geneanno.table[prefilter,],Gene.Names[prefilter],gene.desc.table[prefilter,gene.desc.table.labels.wanted],gene.table[prefilter,],filter.table.pholy[prefilter,filter.table.annotation],indels[prefilter,filter.ID.cols],wanted.muts[prefilter],wanted.muts.coding[prefilter],good.qual[prefilter],wanted.muts.fil[prefilter],ok.in.group[prefilter],maf.lt.all[prefilter,],indels[prefilter,"FILTER"],Denovo.DOM[prefilter,],Compound.Hetero[prefilter,],gerp.scores[prefilter],summary.geno[prefilter,],summary.geno.group[prefilter,],summary.depths[prefilter,],indels[prefilter,indels.order.GT.AN],indels[prefilter,indels.order.SNP.QUAL],filter.table.pholy[prefilter,c("PolyPhen.desc","PolyPhen.scores","SIFT.desc","SIFT.scores")],dom.profile[prefilter],dom.aff.ALT.counts[prefilter],dom.Unaff.REF.counts[prefilter],dom.aff.ALT.counts.C[prefilter],dom.Unaff.REF.counts.C[prefilter],gene.hit.count.fil.all[prefilter,],gene.hit.count.fil.aff.all[prefilter,],gene.hit.count.fil.Unaff.all[prefilter,],phenotype.gene.hit.maf.qual.all[prefilter,],phenotype.gene.hit.count.maf.qual.all[prefilter,],phenotype.gene.hit.maf.qual.cov.all[prefilter,],gene.hit.count.sample.all[prefilter,],summary.depths.group[prefilter,],gene.desc.table[prefilter,gene.desc.table.labels.other],regionanno.table[prefilter,],indels[prefilter,indels.order.2nd],filter.table[prefilter,],indels[prefilter,indels.order.other])
## colnames(gene.hit.count.fil.sample.all) # gene hit counts PER SAMPLE with quality, coverage,MAF  no longer outpt


}else{
summary<-cbind(indels[,c(core.ann,"ID")],geneanno.table,Gene.Names,gene.desc.table[,gene.desc.table.labels.wanted],gene.table,wanted.muts,wanted.muts.coding,good.qual,wanted.muts.fil,maf.lt.all,summary.geno,summary.geno.group,summary.depths,indels[,indels.order.GT.AN],indels[,indels.order.SNP.QUAL],filter.table.pholy,dom.profile,dom.aff.ALT.counts,dom.Unaff.REF.counts,dom.aff.ALT.counts.C,dom.Unaff.REF.counts.C,   gene.hit.count.fil.all,gene.hit.count.fil.aff.all,gene.hit.count.fil.Unaff.all   ,phenotype.gene.hit.maf.qual.all,phenotype.gene.hit.count.maf.qual.all,phenotype.gene.hit.maf.qual.cov.all,phenotype.gene.hit.count.maf.qual.cov.all ,gene.hit.count.sample.all,gene.hit.count.fil.sample.all ,summary.depths.group,gene.desc.table[,gene.desc.table.labels.other],regionanno.table,indels[,indels.order.2nd],filter.table,indels[,indels.order.other])

}

other.sample.info<-{}
if(exists("family.name.prefix")){
gatk.labels<-c("GT","AD","DP","GQ","PL")
gatk.labels<-paste(gatk.labels,"$",sep="")
  for(i in 1:length(gatk.labels)){
other.sample.info<-c(other.sample.info, indels.order.other[ grepl(gatk.labels[i],indels.order.other)])
}}
  
dim(summary)
rm("quality.thresh")
unwanted.columns.full<-c(other.sample.info,unwanted.columns,colnames(summary)[grepl("^key",colnames(summary))]) # remove the unwanted "key columns"

summary<-summary[,colnames(summary)[ !(colnames(summary) %in% unwanted.columns.full)]]

ref.genotype.labels<-expand.labels.to.samples(c("GT"),sample.labels)
ref.genotype.labels


## ref.genotype.labels<-expand.labels.to.samples(c("GT_A"),sample.labels)
## ref.count<-apply(as.matrix(quality.thresh[,ref.genotype.labels]),1,sum) ### All NA's will also get through here quality.thresh[1:5,ref.genotype.labels]
## ref.only<- ref.count>=length(sample.labels)
## sum(ref.only)
## dim(summary)
## summary[summary[,"start"]=="103143402",c("chr","start","end",indels.order.1st)] # a.test<-grep("103143402",indels[,"start"]) ; 

ref.only<-summary.geno[,"ALT.Alleles"]=="0" & summary.geno[,"MISSING.Alleles"]=="0" # ref.only[a.test]
summary<-summary[!ref.only,]

#summary[1:5,]
if(exists("family.name.prefix")){
  #file.out.name<-paste(project,family.name.prefix,sep=".")
  file.out.name<-paste(gsub(".txt$","",target),family.name.prefix,sep=".")
}else{file.out.name<-gsub(".txt$","",target)}

paste(file.out.name,"analysis.txt",sep=".")
paste(file.out.name,"analysis-maf-filtered.txt",sep=".")
getwd()

xx<-try(setwd( analysis.dir  ),silent=TRUE)
if(inherits(xx, "try-error")){system(paste("mkdir",analysis.dir,sep=" "))
                              setwd( analysis.dir  )}
                              

print("Write table summary")
print(file.out.name)
write.table(summary,file=paste(file.out.name,"analysis.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

#  length(maf.lt.all[!ref.only,dim(maf.lt.all)[2]])  # length( wanted.muts[!ref.only])
colnames(maf.lt.all)
## summary.reduced<-summary[ (maf.lt.all[!ref.only,dim(maf.lt.all)[2]]  & wanted.muts[!ref.only]) ,] ### reduce by the LAST maf.lt which is the largest
## summary.reduced<-summary[maf.lt.all[!ref.only,dim(maf.lt.all)[2]],] ### reduce by the LAST maf.lt which is the largest
summary.reduced<-summary[prefilter[!ref.only],] ### reduce by the LAST maf.lt which is the largest
dim(summary.reduced)
print("Print summary reduced")
print(file.out.name)
write.table(summary.reduced,file=paste(file.out.name,"analysis-maf-filtered.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
print("PROJECT ANNOTATION FINISHED")
#summary.reduced[1:5,]

##### SANITY CHECK that there is not a memory failure

# imem<-1
for(imem in 1:length(ref.genotype.labels)){
  if(sum(summary[,ref.genotype.labels[imem]] != indels[!ref.only,ref.genotype.labels[imem]],na.rm=TRUE)!=0){print("ERROR ERROR STRANGE MEMEORY ERROR IN RUN 1")}
}

for(imem in 1:length(ref.genotype.labels)){
  if(sum(summary.reduced[,ref.genotype.labels[imem]] != summary[prefilter[!ref.only],ref.genotype.labels[imem]],na.rm=TRUE)!=0){print("ERROR ERROR STRANGE MEMEORY ERROR IN RUN 2")}
}


print("#################################")
print(file.out.name)
print("#################################")
rm(summary)
rm(summary.reduced)
#save.image(paste(file.out.name,"-full.RData",sep=""))


} #loop over faimilies ifam



} # length(all.sample.labels)>0 that is SKIP analysis and just output summary

print("PROJECT PIPELINE FINISHED......")

if(length(all.sample.labels)<=0){
 
if(class(Gene.Names)=="list"){print("Error Gene.Names a list: L1338");Gene.Names<-unlist(lapply(Gene.Names,function(x) paste(x,collapse="::")))}
annotation.in.filter.table<-c("Consequence",extra.vep.annotations)

colnames(filter.table.pholy)<-gsub(".Embl","",colnames(filter.table.pholy))
filter.table.annotation<-colnames(filter.table.pholy) %in% annotation.in.filter.table
colnames(filter.table.pholy)[filter.table.annotation]<-paste(colnames(filter.table.pholy)[filter.table.annotation],".Embl",sep="")

gene.desc.table.labels.wanted<-c("Description","Gene.Biotype") # gene.desc.table.labels.wanted<-c("description")
gene.desc.table.labels.other<-colnames(gene.desc.table)[!( colnames(gene.desc.table) %in% gene.desc.table.labels.wanted )]

unwanted.columns<-c("refGene" ,"knownGene", "ensGene","gene.table",
                    "snp132::found","1000genome","snp132","CG69","EUR_ASN_AFR_INDEL","ALL_EXON_SNP","1000genomeV2.txt","RNSH_FMHT_maf.txt","ljb_pp2","avsift","ljb_mt","ljb_phylop",
                    colnames(filter.table)[grepl("::maf-filter",colnames(filter.table))],
                     colnames(regionanno.table)[grepl("::HIT",colnames(regionanno.table))],filter.cols.maf.0
                    )
## ##########################################################################################################################################################
## ############## below use only for annotation output
  filter.table.pholy[1:5,]
  gene.desc.table.labels.wanted<-c("Description","Gene.Biotype") # gene.desc.table.labels.wanted<-c("description")
  gene.desc.table.labels.other<-colnames(gene.desc.table)[!( colnames(gene.desc.table) %in% gene.desc.table.labels.wanted )]

  summary<-cbind(indels,geneanno.table,Gene.Names,gene.desc.table[,gene.desc.table.labels.wanted],gene.table,filter.table.pholy[,filter.table.annotation],wanted.muts,wanted.muts.coding,maf.lt.all,gerp.scores,filter.table.pholy[,c("PolyPhen.desc","PolyPhen.scores","SIFT.desc","SIFT.scores")],gene.desc.table[,gene.desc.table.labels.other],regionanno.table,filter.table)


  dim(summary)
  unwanted.columns<-c("junk","shit","SHIT","JUNK","crap")
  unwanted.columns.full<-c(unwanted.columns,colnames(summary)[grepl("^key",colnames(summary))]) # remove the unwanted "key columns"
  summary<-summary[,colnames(summary)[ !(colnames(summary) %in% unwanted.columns.full)]]
  if(exists("family.name.prefix")){
                                        #file.out.name<-paste(project,family.name.prefix,sep=".")
    file.out.name<-paste(gsub(".txt$","",target),family.name.prefix,sep=".")
  }else{file.out.name<-gsub(".txt$","",target)}
  paste(file.out.name,"analysis.txt",sep=".")
  getwd()
  xx<-try(setwd( analysis.dir  ))
  if(inherits(xx, "try-error")){system(paste("mkdir",paste("'",analysis.dir,"'",sep=""),sep=" "));setwd( analysis.dir  )}

  write.table(summary,file=paste(file.out.name,"analysis.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
############## below use only for annotation output
## ##########################################################################################################################################################
} # length(all.sample.labels)<=0
  

#} # ichr loop over projects








