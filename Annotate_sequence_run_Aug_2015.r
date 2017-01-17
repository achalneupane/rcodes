
#options(error=recover)
#/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/Annotate_sequence_run_Aug_2014.r
###### THIS IS the super annotion run USE ALL THE FILTER AND NOVEL DATABASES

# AML
## SKDP - 1 
#UQCCG.data<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/" # path to project data 
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/genotyping/2014-07-18_SKDP_MODY/" # path to project data 
#UQCCG.data<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/" # path to project data 
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/genotyping/" # path to project data 
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/HaplotypeCaller/" # path to project data 
# project<-"2014-07-21_LeoPharma_1st_NexteraSamples" # this is the project directory
#project<-"2014-06-27_All_AOGC_HC" # this is the project directory
#project<-"2014-05-13_SKDP_LGCA" # this is the project directory
# project.name<-"20131115_GroupCallSKDP_MODY_PCC_ChMND" # this is the prefix of the SNP and DINDEL files AND the project output NAME
# # 
# annotate.dir<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes//20131115_GroupCallSKDP_MODY_PCC_ChMND/Annotate"
# analysis.dir<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes//20131115_GroupCallSKDP_MODY_PCC_ChMND/Analysis"

## UQCCG.data<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes"
## project<-"2014-11-07_AML_AOGC_Replication"
## ichr<-10


#args <- commandArgs(TRUE)

args.temp <- commandArgs(TRUE)
if(length(args.temp)>0){args<-args.temp}
print(paste("arguemensts are",args))


ichr <- args[1]
ichr <- as.numeric(ichr)

UQCCG.data <- args[2]
project <- args[3]

if(length(args)>=4){
  max.reads <- as.numeric(args[4])
}
 
print(paste(" UQGGC.data-> ",UQCCG.data, ", project-> ",project,", ichr-> ",ichr,sep=" "))
#q()


project.name <- project
#project.name<-"2014-03-05_AML_RSGB_AOGC_Control" # this is the prefix of the SNP and DINDEL files AND the project output NAME
# 
 annotate.dir <-paste(UQCCG.data,project,"Annotate2",sep="/")
 analysis.dir<-paste(UQCCG.data,project,"Analysis",sep="/")
print(annotate.dir)





skip.annovar.run<-FALSE # if FALSE it just redoes the files annovar summarization IF not already run
update.annovar.annotations<-FALSE # set TRUE will re-read VCF file and force annovar to always run or  FALSE will use the precalcuated Annovar annotation
force.VEP.read<-TRUE ## set to FALSE is default, set to TRUE and update.annovar.annotations set to FALSE will not run annovar but will read in vcf file
GATK.SB<-TRUE
genome.build<-"hg19"
dbSNP.build<-"131"
vcf.type="v4" # "annovar" "v4" "v3" "plink_assoc"

bam.extension<-".ReCal.sort.bam"
combined.extension<-""
snp.extension<-".output.recalibrated.filtered.vcf"  # "_snps.raw.vcf"
indel.extension<-".indelFiltered.vcf"    # "_DINDEL.raw.vcf"  
small.extension<-"_All_VARIANTS.raw.vcf"  # "_All_VARIANTS.raw.vcf"
variant.types<-c("snp","indel","small") ##MUST have a extension type for each indel defined
names(variant.types)<-c("v4","v4","v4") ### define data type for when reading below
use.variants<-c("snp","indel")

#args <- commandArgs(TRUE)
#ichr <- args[1]
#ichr <- as.numeric(ichr)
###############################################

######################################################################################################################
######################################################################################################################
######################################################################################################################
############################################# BEGIN #######################################################
######################################################### Set up the basics for each run #####

## force.functional.update<-TRUE # will re-real functional data normally set to TRUE id have functional.RData file already set to FALSE (will ALWAYS chech is pholyPhen run is TRUE)
## force.GERP.update<-TRUE  # will  get scores normally set to TRUE
## update.annotations.with.vep<-TRUE ##undate gene names and geneanno.table with filter.table.poly... wanted muts traken care of below
## skip.gene.list.matches<-FALSE ### usually false match against lists of gene names


force.functional.update<-TRUE # will re-read functional data normally set to TRUE.  If have functional.RData file already set to FALSE (will ALWAYS chech is pholyPhen run is TRUE)
force.GERP.update<-TRUE  # will  get scores normally set to TRUE # force.GERP.update<-FALSE
update.annotations.with.vep<-TRUE ##undate gene names and geneanno.table with filter.table.poly... wanted muts traken care of below
skip.gene.list.matches<-FALSE ### usually false match against lists of gene names


options(stringsASFactors=FALSE)
options(width=250,max.print=2000)
options("scipen" = 15) # so 150000 does not go to 1.5e5 ## use format() too
library(GenomicFeatures)
library(Rsamtools)
#require(biomaRt)
library(multicore)
#library(parallel)
library(DBI)
library(RMySQL)
num.cores<-4 # 7 usually
## if(!exists("num.cores")){num.cores<-6}
## registerDoMC(cores=num.cores)


###### read in piecewise  AND filter out columns ## default read all and no filtering ok up to @300 samples####
if(!exists("max.reads")){max.reads<-0} # 0- get ALL IN ONE GO  :
if(!exists("remove.extensions")){remove.extensions<-{}} # remove.extensions<-c(".GQ",".DP",".PL") in above typical
if(!exists("remove.cols")){remove.cols<-{}}
if(!exists("force.VEP.read")){if(update.annovar.annotations){force.VEP.read<-TRUE}else{force.VEP.read<-FALSE} }# is !exist force.VEP.read : if update.annovar.annotations=TRUE force.VEP.rea
# !update.annovar.annotations & !force.VEP.read # update.annovar.annotations==FALSE force.VEP.read==TRUE forces vep read but does not run annovar -dangerous)



###################################################

## NA get all
## 1000 samples - 10000 
## 300 sample 70000 (7 cores) < 14 Gb
## 55 samples 20000 reads < 5Gb # 550,000 > 32Gb memory blown
## 40 samples 100,000 reafs 10Gb
###############10000 -1000 samples about 23Gb RAM' 12hrs
###############ALL   -58 samples about  23Gb 30 mins

project.dir<-paste(UQCCG.data,project,sep="/")

bam.dir<-paste(project.dir,"BAM",sep="/")
xx<-try(setwd( bam.dir ),silent=TRUE)
if(inherits(xx, "try-error")){bam.dir<-project.dir}


files<-dir(project.dir)
if(!("bam.dir") %in% files ){snp.dir<-paste(project.dir,"SNPs",sep="/")}else if(!grepl(project,snp.dir)){snp.dir<-paste(project.dir,"SNPs",sep="/")} # | !grepl(project,snp.dir) incase just doing another peoject


if(!exists("snp.dir")){snp.dir<-paste(project.dir,"SNPs",sep="/")}else if(!grepl(project,snp.dir)){snp.dir<-paste(project.dir,"SNPs",sep="/")} # | !grepl(project,snp.dir) incase just doing another peoject
small.dir<-paste(project.dir,"SNPs",sep="/")
indel.dir<-paste(project.dir,"DINDELs",sep="/")
polyphen.dir<-paste(project.dir,"PolyPhen2",sep="/")

#polyphen.dir2<-paste(project.dir,"test",sep="/")


xx<-try(setwd( snp.dir ),silent=TRUE)
if(inherits(xx, "try-error")){snp.dir<-project.dir;print("ERROR Could not find SNP dir")}
xx<-try(setwd( small.dir ),silent=TRUE)
if(inherits(xx, "try-error")){small.dir<-project.dir}
xx<-try(setwd( indel.dir ),silent=TRUE)
if(inherits(xx, "try-error")){indel.dir<-project.dir;print("ERROR Could not find SNP dir")}
xx<-try(setwd( polyphen.dir ),silent=TRUE)
if(inherits(xx, "try-error")){
  xx<-try(system(paste("mkdir",polyphen.dir,sep=" ")),silent=TRUE)
  xx<-try(paste("mkdir",paste("'",polyphen.dir,"'",sep=""),sep=" "),silent=TRUE)
  if(inherits(xx, "try-error")){
  system(polyphen.dir<-project.dir) }
       }

QC.dir<-paste(project.dir,"QC",sep="/")
if(!exists("annotate.dir")){annotate.dir<-"/media/Bioinform-D/Research/annotate"}
#analysis.dir<-"/media/Bioinform-D/Research/annotate/final"

xx<-try(setwd( annotate.dir ),silent=TRUE)
system(paste("mkdir",paste("'",annotate.dir,"'",sep=""),sep=" "));setwd(annotate.dir  )

xx<-try(setwd( annotate.dir  ),silent=TRUE)
if(inherits(xx, "try-error")){
  xx<-try(paste("mkdir",paste("'",annotate.dir,"'",sep=""),sep=" "),silent=TRUE)
  xx<-try(setwd(annotate.dir ),silent=TRUE)
  if(inherits(xx, "try-error")){
  system(annotate.dir<-project.dir) }
       }


code.dir<-"/mnt/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/"
#code.dir<-"/home/mmarshall/PBSHOME/PaulGenotypingPipeline"
db.details<-list(user="gerp2",pass="GeRpUQ8!6", dbname="gerp", host="ga-apps.di.uq.edu.au")





variant.types<-variant.types[variant.types %in% use.variants]
splice.threshold<-5 # default is 2 in annovar. I have seen  3 and 5 are used in dbSNP
maf.threshold<-0.0  #  MAF threshold for annovar calling zero useful to get back all results !!do not modify!!
#af.threshold.filter<-c(0.001,0.01,0.025,0.05) # MAF threshold for boolean filtering operations NOVEL always done
#maf.threshold.filter<-0.05
core.ann<-c("chr","start","end","REF","ALT","TYPE") # out put to annals programs and need foe colun labels


extra.vep.annotations<-c("Uploaded_variation","Gene","Feature","Protein_position","Amino_acids")



vep.types<-c("not_assigned","stop_gained","stop_lost","missense_variant","splice_acceptor_variant","splice_donor_variant","splice_region_variant","initiator_codon_variant","stop_retained_variant","incomplete_terminal_codon_variant","frameshift_variant","inframe_deletion","inframe_insertion","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_exon_variant","NC_stop_gained","NC_stop_lost","NC_splice_acceptor_variant","NC_splice_donor_variant","NC_splice_region_variant","NC_initiator_codon_variant","NC_stop_retained_variant","NC_non_coding_exon_variant","NC_incomplete_terminal_codon_variant","NC_3_prime_UTR_variant","mature_miRNA_variant","NC_5_prime_UTR_variant","TF_binding_site_variant","TFBS_ablation","TFBS_amplification","regulatory_region_variant","intron_variant","NC_intron_variant","synonymous_variant","coding_sequence_variant","NC_synonymous_variant","upstream_gene_variant","downstream_gene_variant","intergenic_variant","NC_intergenic_variant","NMD_transcript_variant","nc_transcript_variant","NC_nc_transcript_variant","feature_truncation","feature_elongation")

#interesting.mutations<-c("nonsynonymous SNV","splicing","frameshift substitution","stopgain SNV","stoploss SNV") #3 older annovar
## possible.mutations<-c("downstream","frameshift substitution","intergenic","intronic","ncRNA","nonframeshift substitution","nonsynonymous SNV","splicing","stopgain SNV","stoploss SNV","synonymous SNV","unknown","upstream","upstream;downstream","UTR3","UTR5","UTR5;UTR3") ## older annovar

interesting.mutations<-c("frameshift substitution","nonframeshift substitution","nonframeshift deletion","nonframeshift insertion","frameshift deletion","frameshift insertion","nonsynonymous SNV","stopgain SNV","stoploss SNV","splicing","ncRNA_exonic")

## "ncRNA_UTR3","ncRNA_UTR5","ncRNA_UTR5;ncRNA_UTR3","UTR3","UTR5","UTR5;UTR3" ## other UTR's I could use

possible.mutations<-c("frameshift substitution","nonframeshift substitution","downstream","frameshift deletion","frameshift insertion","intergenic","intronic","ncRNA_exonic","ncRNA_intronic","ncRNA_splicing","ncRNA_UTR3","ncRNA_UTR5","ncRNA_UTR5;ncRNA_UTR3","nonframeshift deletion","nonframeshift insertion","nonsynonymous SNV","splicing","stopgain SNV","stoploss SNV","synonymous SNV","unknown","upstream","upstream;downstream","UTR3","UTR5","UTR5;UTR3")


###################### define annotation databases ################


if(genome.build=="hg18"){

anno.DB.location.core<-"/mnt/UQCCG/Software/annovar/humandb"
anno.DB.location<-paste(anno.DB.location.core,genome.build,sep="/")

  
regionanno.DB<-c("gwasCatalog","segdup","snp130","tfbs","dgv","mce28way","mirnatarget","mirna","rnaGene","cpgIslandExt")  ##use colsWanted -all so must have colun labels defined c("gwasCatalog","omimGene","segdup","snp132","tfbsConsSites","dgv")
names(regionanno.DB)<-c("gwasCatalog","segdup","snp130","tfbs","dgv","conserved_28way","mirnatarget","mirna","rnaGene","cpgIslandExt")  ##use colsWanted -all so must have colun labels defined c("gwasCatalog","omimGene","segdup","snp132","tfbsConsSites","dgv")

regionanno.DB.score<-c("wgEncodeRegDnaseClustered","wgEncodeRegTfbsClustered") # region annotation the require a score REQUIRE -scorecolumn 5 (not a core annovar database)
regionanno.DB.score.col<-c("5","5") ## number of the column to use 
names(regionanno.DB.score)<-c("wgEncodeRegDnaseClustered","wgEncodeRegTfbsClustered") 


geneanno.DB<-c("refGene","knownGene","ensGene") # returns 2 extra columns
names(geneanno.DB)<-c("refGene","knownGene","ensGene")

filter.DB<-c("snp132","1000g2010jul_ceu","1000g2010jul_yri","1000g2010jul_jptchb")  # returns 2  extra columns (DB,score) c("snp132","1000g2010nov_all")
names(filter.DB)<-c("snp132","1000genome_CEU","1000genome_YRI","1000genome_JPTCHB") # these will be the names in the summary table

generic.filter.DB<-c("hg18_ljb_gerp++.txt","hg18_esp5400_all.txt","hg18_esp5400_ea.txt","hg18_esp5400_aa.txt","hg18_cg69.txt")    # returns 2  extra columns (DB,score)

names(generic.filter.DB)<-c("GERP","NHLBI_5400_all","NHLBI_5400_EUR","NHLBI_5400_AFR","CG69")

function.filter.DB<-c("ljb_pp2","avsift","ljb_mt","ljb_phylop") # returns 2  extra columns (DB,score)
names(function.filter.DB)<-c("pholyphen","sift","mut.taster","phylo")

the.combined.DBs<-c(filter.DB,generic.filter.DB,function.filter.DB) # used when processing the  have two columns the  key : eg :  DB, MAF,chr,start,end, REF,ALT,score
length(the.combined.DBs)


###### CHOOSE WHAT TO USE IN THE ANALYSIS

  #  the **names** of all DB that have all frequency in the 6th column to use for maf filtering


all.filter.cols.novel<-c("snp135","NHLBI_5400_all","snp132","1000genome","snp137","snp132","CG69","EUR_ASN_AFR_INDEL","ALL_EXON_SNP","1000genomeV2","RNSH_FMHT","MCTO","AOGC-NGS","AOGC-NGS-All")
all.filter.cols.maf<-c("snp135","NHLBI_5400_all","1000genome","snp134","snp137","CG69","EUR_ASN_AFR_INDEL","ALL_EXON_SNP","1000genomeV2","RNSH_FMHT","MCTO","AOGC-NGS","AOGC-NGS-All")                         

}

if(genome.build=="hg19"){

  anno.DB.location.core<-"/mnt/UQCCG/Software/annovar/humandb"
anno.DB.location<-paste(anno.DB.location.core,genome.build,sep="/")

  #### fix problem with grep++elem in variable names - just remove ++ from the DB name and sue that dbnam in the ucsc.tables.names.r

regionanno.DB<-c("jaxQtlAsIs","gwasCatalog","segdup","tfbsConsSites","dgv","mirna")  ##use colsWanted -all so must have colun labels defined c("jaxQtlAsIs","gerp++elem","Gerp","gwasCatalog","omimGene","segdup","tfbsConsSites","dgv","mirna") 
names(regionanno.DB)<-c("jaxQtlAsIs","gwasCatalog","segdup","tfbsConsSites","dgv","mirna")  ##use colsWanted -all so must have c("jaxQtlAsIs","gerpelem","Gerp","gwasCatalog","omimGene","segdup","tfbsConsSites","dgv","mirna")


regionanno.DB.score<-c("wgEncodeRegDnaseClustered","wgEncodeRegTfbsClustered") # region annotation the require a score REQUIRE -scorecolumn 5 (not a core annovar database)
regionanno.DB.score.col<-c("5","5") ## number of the column to use 
names(regionanno.DB.score)<-c("wgEncodeRegDnaseClustered","wgEncodeRegTfbsClustered") 

geneanno.DB<-c("refGene","knownGene","ensGene") # returns 2 extra columns
names(geneanno.DB)<-c("refGene","knownGene","ensGene")


## these ones have a score  column
filter.DB<-c("snp138","1000g2012apr_all","1000g2012apr_asn")  # returns 2  extra columns (DB,score) c("snp132","1000g2011may_all")
names(filter.DB)<-c("ID","1000genome","1000genome_asian") # these will be the names in the summary table


## generic.filter.DB<-c("hg19_1000_genomes_feb_maf.txt","hg19_snp135_maf.txt","hg19_snp135_pathalogical_maf.txt","hg19_snp135_clinical_maf.txt","hg19_Gerp_Annovar_2.0.txt","hg19_NHBLI_5400_ALL_maf.txt","hg19_NHBLI_5400_EA_maf.txt","hg19_NHBLI_5400_AA_maf.txt","hg19_NHBLI_5400_ALL_QCFail_maf.txt","hg19_NHBLI_5400_EA_QCFail_maf.txt","hg19_NHBLI_5400_AA_QCFail_maf.txt","hg19_snp134_maf.txt","hg19_cg69.txt","hg19_EUR_ASN_AFR_DINDEL_maf.txt","hg19_ALL_EXON_SNP_maf.txt","ALL.wgs.phase1_maf.txt","RNSH_FMHT_maf.txt","MCTO_maf.txt","AOGC-NGS_ALL_maf.txt","AOGC-NGS_ALL_QCFail_maf.txt")

## names(generic.filter.DB)<-c("1000genome_mine","snp135","snp135_pathalogical","snp135_clinical","GERP","NHLBI_5400_ALL","NHLBI_5400_EUR","NHLBI_5400_AFR","NHLBI_5400_ALL_QCFail","NHLBI_5400_EUR_QCFail","NHLBI_5400_AFR_QCFail","snp134","CG69","EUR_ASN_AFR_INDEL","ALL_EXON_SNP","1000genomeV2","RNSH_FMHT","MCTO","AOGC-NGS_ALL","AOGC-NGS_ALL_QCFail")

generic.filter.DB<-c("hg19_popfreq_max.txt","hg19_esp6500si_all.txt","hg19_NHBLI_6500_ALL_maf.txt","hg19_esp6500_ea.txt","hg19_esp6500_aa.txt","hg19_snp141_pubmed_maf.txt","hg19_snp141_omim_maf.txt","hg19_snp141_clinical_maf.txt","hg19_snp141_maf.txt","hg19_snp141_papu_maf.txt",
                     "hg19_1000_genomes_feb_maf.txt","hg19_snp137_maf.txt","hg19_snp137_pathalogical_maf.txt","hg19_snp137_clinical_maf.txt",
                     "hg19_NHBLI_5400_ALL_maf.txt","hg19_NHBLI_5400_EA_maf.txt","hg19_NHBLI_5400_AA_maf.txt","hg19_NHBLI_5400_ALL_QCFail_maf.txt","hg19_NHBLI_5400_EA_QCFail_maf.txt","hg19_NHBLI_5400_AA_QCFail_maf.txt",
                     "hg19_cg69.txt","hg19_EUR_ASN_AFR_DINDEL_maf.txt",
                     "AOGC-NGS_ALL_maf.txt","AOGC-NGS_ALL_QCFail_maf.txt","hg19_Chinese_maf.txt")    # returns 2  extra columns (DB,score)

names(generic.filter.DB)<-c("PopFreqMax","NHBLI_6500_ANNOVAR_ALL","NHBLI_6500_ALL","NHBLI_6500_EA","NHBLI_6500_AA","snp141_pubmed","snp141_omim","snp141_clinical","snp141","snp141papu",
                            "1000genome_mine","snp137","snp137_pathalogical","snp137_clinical",
                            "NHLBI_5400_ALL","NHLBI_5400_EUR","NHLBI_5400_AFR","NHLBI_5400_ALL_QCFail","NHLBI_5400_EUR_QCFail","NHLBI_5400_AFR_QCFail",
                            "CG69","EUR_ASN_AFR_INDEL",
                            "AOGC-NGS_ALL","AOGC-NGS_ALL_QCFail","Chinese")

## 137->141
## 135->137
#### thies ones don't have a score column

function.filter.DB<-c("cadd","ljb23_sift", "ljb23_pp2hdiv", "ljb23_pp2hvar", "ljb23_lrt", "ljb23_mt", "ljb23_ma", "ljb23_fathmm", "ljb23_metasvm", "ljb23_metalr", "ljb23_gerp++", "ljb23_phylop", "ljb23_siphy","cosmic68")
names(function.filter.DB)<-c("CADD","sift","pholyphen_GWAS","pholyphen","LRT","mut.taster","mut.accessor","FATHMM","MetaSVM","MetaLR","gerp.ann","PhyloP","SiPhy","cosmic68") 
  

## function.filter.DB<-c("ljb_pp2","avsift","ljb_mt","ljb_phylop","cosmic68") # returns 2  extra columns (DB,score) c("gerp++gt2","ljb_gerp++")
## names(function.filter.DB)<-c("pholyphen","sift","mut.taster","phylo","cosmic") 
                   #                  c("GERP","ljb_gerp")

## Score (dbtype) 	# variants in LJB23 build hg19 	Categorical Prediction
## SIFT (sift) 	77593284 	D: Deleterious (sift<=0.05); T: tolerated (sift>0.05)
## PolyPhen 2 HDIV (pp2_hdiv) 	72533732 	D: Probably damaging (>=0.957), P: possibly damaging (0.453<=pp2_hdiv<=0.956); B: benign (pp2_hdiv<=0.452)
## PolyPhen 2 HVar (pp2_hvar) 	72533732 	D: Probably damaging (>=0.909), P: possibly damaging (0.447<=pp2_hdiv<=0.909); B: benign (pp2_hdiv<=0.446)
## LRT (lrt) 	68069321 	D: Deleterious; N: Neutral; U: Unknown
## MutationTaster (mt) 	88473874 	A" ("disease_causing_automatic"); "D" ("disease_causing"); "N" ("polymorphism"); "P" ("polymorphism_automatic"
## MutationAssessor (ma) 	74631375 	H: high; M: medium; L: low; N: neutral. H/M means functional and L/N means non-functional
## FATHMM (fathmm) 	70274896 	D: Deleterious; T: Tolerated
## MetaSVM (metasvm) 	82098217 	D: Deleterious; T: Tolerated
## MetaLR (metalr) 	82098217 	D: Deleterious; T: Tolerated
## GERP++ (gerp++) 	89076718 	higher scores are more deleterious
## PhyloP (phylop) 	89553090 	higher scores are more deleterious
## SiPhy (siphy) 	88269630 	higher scores are more deleterious 

  

the.combined.DBs<-c(filter.DB,generic.filter.DB,function.filter.DB) # used when processing the  have two columns the  key : eg :  DB, MAF,chr,start,end, REF,ALT,score
length(the.combined.DBs) ## ONLY THE COMBINED DBs: the.combined.DBs CAN BE USED FOR MAF FILTERING


###### CHOOSE WHAT TO USE IN THE ANALYSIS

  #  the **names** of all DB that have all frequency in the 6th column to use for maf filtering

## filter.cols.novel<-c("1000genome","1000genome_asian","1000genome_mine","1000genomeV2","snp137","snp137_clinical","snp135","NHLBI_5400_ALL","snp134","CG69","EUR_ASN_AFR_INDEL","ALL_EXON_SNP","1000genomeV2","AOGC-NGS_ALL") ## f just if found or not"RNSH_FMHT","MCTO"
## filter.cols.maf<-c("1000genome_mine","snp135","NHLBI_5400_ALL","1000genome","snp134","1000genomeV2","AOGC-NGS_ALL")

## all.filter.cols.novel<-c("NHBLI_6500_ALL","1000genome_asian","1000genome_mine","snp135","NHLBI_5400_ALL","snp132","1000genome","snp137","CG69","EUR_ASN_AFR_INDEL","ALL_EXON_SNP","1000genomeV2","MCTO","AOGC-NGS_ALL","AOGC-NGS_ALL_OLD")
## all.filter.cols.maf<-c("1000genome_asian","1000genome_mine","snp135","NHBLI_6500_ALL","NHLBI_5400_ALL","1000genome","snp137","CG69","EUR_ASN_AFR_INDEL","ALL_EXON_SNP","1000genomeV2","AOGC-NGS_ALL","AOGC-NGS_ALL_OLD")    

#all.filter.cols.maf<-c("snp135","NHLBI_5400_ALL","1000genome","snp134","snp132","CG69","EUR_ASN_AFR_INDEL","ALL_EXON_SNP","1000genomeV2","RNSH_FMHT","MCTO","AOGC-NGS","AOGC-NGS-ALL")  
}

if(genome.build=="mm9"){

anno.DB.location.core<-"/mnt/UQCCG/Software/annovar/mousedb"
anno.DB.location<-paste(anno.DB.location.core,genome.build,sep="/")
  
regionanno.DB<-c("jaxPhenotype","segdup","miRNA","mce30way")  ##use colsWanted -all so must have colun labels defined c("gwasCatalog","omimGene","segdup","snp132","tfbsConsSites","dgv")
names(regionanno.DB)<-c("MGI_phenotype","segdup","miRNA","mce30way") # phastConsElements30way

geneanno.DB<-c("refGene","knownGene","ensGene") # returns 2 extra columns
names(geneanno.DB)<-c("refGene","knownGene","ensGene")

filter.DB<-c("snp128")  # returns 2  extra columns (DB,score) c("snp132","1000g2010nov_all")
names(filter.DB)<-c("snp128") # these will be the names in the summary table

generic.filter.DB<-c("mm9_snp128_maf.txt","Harwell-ENU_maf.txt","TCGM-ENU_maf.txt","mm9_snpSANGER_maf.txt","mm9_snpSANGER-ALL_maf.txt","ENU18Flt3ITD21_maf.txt","ENU18Flt3ITD24_maf.txt","ENU20Flt3ITD37_maf.txt","Harwell-Control_maf.txt","Harwell-Control-ALL_maf.txt","TCGM-ENU-ALL_maf.txt","Harwell-Bone-ALL_maf.txt","Harwell-Bone_maf.txt","Harwell-Kidney-ALL_maf.txt","Harwell-Kidney_maf.txt")    # returns 2  extra columns (DB,score) ,"Harwell-ENU_maf.txt"
names(generic.filter.DB)<-c("SNP128","Harwell-ENU","TCGM-ENU","snpSANGER","snpSANGER-ALL","ENU18Flt3ITD21","ENU18Flt3ITD24","ENU20Flt3ITD37","Harwell-Control","Harwell-Control-ALL","TCGM-ENU-ALL","Harwell-Bone-ALL","Harwell-Bone","Harwell-Kidney-ALL","Harwell-Kidney") # ,"Harwell-ENU"

function.filter.DB<-c() # returns 2  extra columns (DB,score) These are not avalaibel from mouse from annovar 
names(function.filter.DB)<-c()

the.combined.DBs<-c(filter.DB,generic.filter.DB) # used when processing the  have two columns the  key : eg :  DB, MAF,chr,start,end, REF,ALT,score
length(the.combined.DBs)

## all.filter.cols.novel<-c("snp128","SNP128","TCGM-ENU","snpSANGER","snpSANGER-ALL","Harwell-Control","Harwell-ENU","ENU18Flt3ITD21","ENU18Flt3ITD24","ENU20Flt3ITD37","Harwell-Control-ALL","TCGM-ENU-ALL","Harwell-Bone-ALL","Harwell-Bone","Harwell-Kidney-ALL","Harwell-Kidney") # ALL POSSIBLE NOVEL FILTER DATASETS this is used in the super annotation run
## all.filter.cols.maf<-c("snp128","SNP128","TCGM-ENU","snpSANGER","snpSANGER-ALL","Harwell-Control","Harwell-ENU","ENU18Flt3ITD21","ENU18Flt3ITD24","ENU20Flt3ITD37","Harwell-Control-ALL","TCGM-ENU-ALL","Harwell-Bone-ALL","Harwell-Bone","Harwell-Kidney-ALL","Harwell-Kidney") # ALL POSSIBLE MAF FILTER DATASETS this is used in the super annotation run


###### CHOOSE WHAT TO USE IN A SPECIFIC  ANALYSIS




}

############################################################
############################################################

##############################INPORTANT TECHNICALITY
## ANNOVAR fails if length(REF)> length(ALT)
### ANNOVAR DB does not need to be ordered
## the indel definations  used in GATK ARE DIFFERENT to annovar in format .EG
## 1	98619	98622	TGAC	-	56	685  # dindel vcf using convert2annovar.pl
## 1	98618	98622	TTGAC	T                    # dindel processed by me AND used in my DINDEL definations
## So I process dbsnp132, DINDEL so the filter will work fine for those ut annovar formated files will be wrong.
## I will use my formatting but will but in the future will need to modify 1) the was DINDELS are represented in GATK AND then I can use convert2annovar.pl
## I could convert my vcf files using convert2annovar but I will need to change it to get the columns I want

####effect is
## 1	98619	98622	TGAC	-	56	685 # is not detected in a test BUT
## generic	0.0295571	1	98618	98622	TTGAC	T	0.0295571	56	AF=0.0295571;NS=179;DP=685;HP=2;NFS=44;NRS=31 # is found
## that is it is just matching REF and ALT in the DB and the input file to see if they are the same.

##########################################################  functions required at beginning
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
generic.filter.DB
#maf.threshold.filter<-sort(maf.threshold.filter) # cause when subset the final data use the largest MAf liste

setwd(code.dir)
source("ucsc.table.names.r")   # load in the UCSC tables these use the db file names and not their lable-names 
source("annotate_SNPs_subroutines.r")
source("ucsc.table.names.processor.r")
##################################################################################                   
######################################### Predefined variables required
##################################################################################




## if(!skip.annovar.run){  ### already did annovar and an RData file for the project exists
## ########################################################################################################################
########################################################################################################################
########################################################################################################################
################################################## generate the files lists for variant types
files<-dir(bam.dir)
bam.files<-files[grepl("bam$",files)]
names(bam.files)<-gsub(bam.extension,"",bam.files)

i<- 2
for(i in 1:length(variant.types)){  ### CHOOSE BY project.name * extension
  eep <- "variant type"
  print(eep)
  print(variant.types[i])
files<-dir( eval(as.name( paste(variant.types[i],"dir",sep=".") )) )
the.extension<-paste(eval(as.name( paste(variant.types[i],"extension",sep=".") )),"$",sep="")
the.extension<-paste(combined.extension,the.extension,sep="")
the.files<-files[grepl(the.extension ,files)]

  if(project.name=="2015-03-16_AllAMLandLung" & variant.types[i]=="indel"){
 the.files<-the.files[grepl(paste("^","AML",sep=""),the.files) & grepl(paste("^","AML.chr",sep=""),the.files) ]
 names(the.files)<-gsub(the.extension,"",gsub("AML","2015-03-16_AllAMLandLung",the.files))

 
  }else{
the.files<-the.files[grepl(paste("^",project.name,sep=""),the.files) | grepl(paste("^chr\\S+",project.name,sep=""),the.files) ]
names(the.files)<-gsub(the.extension,"",the.files)
}


assign( paste(variant.types[i],"files",sep="."),value=the.files)
}

snp.files
indel.files
 length(snp.files)
 length(indel.files)
## names(indel.files)[!(names(indel.files) %in% names(snp.files))]



#################################
## IF HAVE DIFFERENT CHROMOSOMES
#################################
## loop over project.name_extra:  project_chr1.snps  project_chr1.dindels
## uses the "SNP" to capture the prefic : so always assumes there are SNPs!


######## remove chrGL files for now ... as no SNPs are defined don't know about genes either?
uncharacterised.chr<-grepl("chrGL",snp.files)
snp.files<-snp.files[!uncharacterised.chr]

if(exists("indel.files")){ ## just to annotaion of bims then no indels file
uncharacterised.chr<-grepl("chrGL",indel.files)
indel.files<-indel.files[!uncharacterised.chr]
}


#print("here got indels")
## snp.files<-snp.files[c(-1,-2,-10)]
## snp.files<-snp.files[-1*c(1:8)]
## snp.files<-snp.files[-1*c(1:3)]
## snp.files<-snp.files[-1]
## snp.files<-snp.files[-1]
## snp.files<-snp.files[-1*c(1:2)]
## snp.files
## snp.files<-snp.files[c(-1,-2,-3,-4,-5,-6,-7,-8,-10)]
#ichr<-24
#ichr <- 20
#for(ichr in 1:length(snp.files)){

indels<-{}
geneanno.table<-{}
gene.desc.table<-{}
regionanno.table<-{}
filter.table<-{}



## for(ichr in 2:9){
print(paste("doing chromosome SNP file: ",snp.files[ichr],sep=""))
print(paste("doing chromosomes INDEL file: ",indel.files[ichr],sep=""))
  
the.extension<-paste(eval(as.name( paste(variant.types[variant.types=="snp"],"extension",sep=".") )),"$",sep="")
the.extension<-paste(combined.extension,the.extension,sep="")
summary.file<-gsub(the.extension,"",snp.files[ichr])


##################### GET the project name and the variant imput file names ################
#summary.file<-paste(project.name,combined.extension,sep="")
print(summary.file)

if(length(snp.files)>1){
target<-summary.file
target<-paste(target,"txt",sep=".")
target
}else{
target<-paste(project.name,"txt",sep=".") # "fam26.txt"
target
}

sample.files<-{}

## if(project=="TGCM-AML"){  
## for(i in 1:length(variant.types)){
## the.files<-eval(as.name(paste(variant.types[i],"files",sep=".")))
## if(i==1){
## sample.files[i]<-the.files[grepl(summary.file,names(the.files))]
## }
## if(i==2){
## sample.files[i]<-the.files[grepl(gsub(".output","",names(summary.file)),names(the.files))]
## }  
## }
## names(sample.files)<-variant.types
## print(sample.files)

## }else{

  
for(i in 1:length(variant.types)){
the.files<-eval(as.name(paste(variant.types[i],"files",sep=".")))
if(length(the.files)!=0){
sample.files[i]<-the.files[grepl(paste(summary.file,"$",sep=""),names(the.files))]
names(sample.files)[i]<-variant.types[i]
} # no files then none of that variant class available
}

print(sample.files)
#} ## "TGCM-AML"

############# sample.files now looks like this:
## > sample.files
##                       snp                     indel 
##   "FMHT_All_snps.raw.vcf" "FMHT_All_DINDEL.raw.vcf"
####################################################################/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/DINDELs/AML.chr14.indelFiltered.vcf#
######################## READ IN THE DATA 

sample.grs<-paste(names(sample.files),"grs",sep=".")


if(!skip.annovar.run){  ### skip annovar

  
if(!update.annovar.annotations & !force.VEP.read){ # reload previous data if already read in and saved - use to re-run annovar for updates

print(annotate.dir)
  setwd(annotate.dir)
  print(paste("Loading: ",paste(gsub(".txt$","",target),"_indel_files.RData",sep="")))
  load(paste(gsub(".txt$","",target),"_indel_files.RData",sep=""))     
  if(!grepl("^chr",indels[1,"chr"])){
    key.indels<-build.key(indels,core.ann,add.chr.label=TRUE)
    indels[,"chr"]<-gsub("chr","",indels[,"chr"])
  }else{key.indels<-build.key(indels,core.ann)      }
  files.in.annotate<-dir(annotate.dir) ## needs there for pholyphen etc

}else{

print(paste("Reading in vcf genotypes:-> ",UQCCG.data, ", project-> ",project,", ichr-> ",ichr," Files: ",toString(sample.files),sep=" "))
setwd(code.dir)
if(vcf.type=="v4" | vcf.type=="v3")source("read.sample.files.MultiAllele.piecewise.r")
#if(vcf.type=="mach_assoc"){source("read.assoc.files.r")}
if(vcf.type=="bim"){source("read.plink.bim.files.r")}
if(vcf.type=="plink_assoc"){source("read.plink.assoc.files.r")}
if(vcf.type=="legend"){source("read.legend.files.r")}
if(vcf.type=="annovar"){source("read.annovar.files.r")}
#########################################################
######################################################## FINISHED READING IN DATA
#########################################################

######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################

if(exists("snp.grs")){print(dim(snp.grs))}
if(exists("indel.grs")){print(dim(indel.grs))}



############### Merge all indel types AND check the the same samples are included in each ######

################################# NOW concatinate the individual genotype calls 
indels<-{}
#### make sure all columns have the same metadata columns ### Collect all metacols
meta.cols<-{}
i<-1

for(i in 1:length(sample.grs)){
meta.cols.temp<-colnames( eval(as.name( sample.grs[i]) ) )
#meta.cols.temp<-colnames(values( eval(as.name( sample.grs[i]) )) )  # if sample.gr is a Grange object
if(exists(sample.grs[i])){print( paste(sample.grs[i],"rows:", dim( eval(as.name( sample.grs[i]) ))[1],sep=" ") )}

if(i==1){meta.cols<-meta.cols.temp; the.samples<-gsub(".GT$","",meta.cols.temp[grepl(".GT",meta.cols.temp)]) }else{
  meta.cols<-c(meta.cols,meta.cols.temp[!(meta.cols.temp %in% meta.cols) ])}

the.samples.next<-gsub(".GT$","",meta.cols.temp[grepl(".GT",meta.cols.temp)])
if( (length(the.samples.next) != length(the.samples)) | sum(the.samples %in%  the.samples.next) !=length(the.samples)){
  print("----------------------------------------------------------------------")
  print("ERROR ERROR samples do not match")
  print(the.samples)
  print(the.samples.next)
}
}


###############
##########combine all Granges order of elementMetavalues is important
for(i in 1:length(sample.grs)){
a.grs<-eval(as.name( sample.grs[i]) )

 if(exists(sample.grs[i])){ rm(list=as.character(sample.grs[i]) )  }


## the.meta.cols<-colnames(eval(as.name( sample.grs[i]) ) )
the.meta.cols<-colnames(a.grs)
## the.meta.cols<-colnames(values( a.grs) ) # if sample.gr is a Grange object
missing.meta.cols<-meta.cols[!(meta.cols %in% the.meta.cols)]

if(length(missing.meta.cols)>0 | (sum(meta.cols[1:length(the.meta.cols)] != the.meta.cols)!=0) ){  
  extra.blank<-matrix(data=NA,nrow=dim(a.grs)[1],ncol=length(missing.meta.cols))
  colnames(extra.blank)<-missing.meta.cols
  a.grs<-cbind(a.grs,extra.blank)[,meta.cols]
  ## values(a.grs)<-cbind(values(a.grs),DataFrame(extra.blank))[,meta.cols]  # if sample.gr is a Grange object NEED DataFrame call!
}

if(i==1){indels<-a.grs}else{indels<-rbind(indels,a.grs)}
## if(i==1){All.grs<-a.grs}else{All.grs<-c(All.grs,a.grs)} # if sample.gr is a Grange object
}

rm("a.grs")

#All.grs[1:3,]
## length(All.grs)
## indels<-as.data.frame(All.grs) ### convert genomic ranges to a dta frame # if sample.gr is a Grange object
## setwd(annotate.dir)
## save(list=c("snp.grs","indel.grs","project.name"),file=paste(gsub(".txt$","",target),"_tmp_indel_files.RData",sep=""))
## setwd("/media/ga-apps/UQCCG/Data/Sequence_Genotypes/2012-11-12_GroupCallQIMR/DINDELs")
## load(paste(gsub(".txt$","",target),"_tmp_indel_files.RData",sep=""))

dim(indels)
for(i in 1:length(sample.grs)){
  if(exists(sample.grs[i])){ rm(list=as.character(sample.grs[i]) )  }
}


colnames(indels)[colnames(indels)=="seqnames"]<-"chr"
indels[,"chr"]<-gsub("chr","",indels[,"chr"])
key.indels<-build.key(indels,core.ann) #sep is always ":"
#key.indels<-paste(indels[,"chr"],indels[,"start"],indels[,"end"],indels[,"REF"],indels[,"ALT"],indels[,"TYPE"],sep=":")

if(length(key.indels)!=length(unique(key.indels))){
  print("WARNING duplicate entries found") ## Hapmap / 1000 Genomes an even within dbSNP can have duplicates key above shoudl work but if not
  posns<-match(unique(key.indels),key.indels)
  missing<-is.na(posns)
  print((key.indels[missing]))
  indels<-indels[posns[!missing],]
  ## All.grs<-All.grs[posns[!missing],] # if sample.gr is a Grange object
#  indels<-as.data.frame(All.grs) ### convert genomic ranges to a dta frame
#  colnames(indels)[colnames(indels)=="seqnames"]<-"chr"
#  indels[,"chr"]<-gsub("chr","",indels[,"chr"])
  key.indels<-build.key(indels,core.ann)
}

col.ann<-colnames(indels)

######################################################################################################################
#####################################################################################################################

## key.indels[1:5]
## "1" "1666447"  "1666451"  "+"    "rs112262693" "GGTTT" "GTTTT" "22893.34" "PASS"                         "indel:1666447"  "GTTTT"
## TGAACTGTGCACTCATCCATTTTCTG  [-/GTTT/TTTT]  TTTTTTGTTTGTTTTTATTATTTTT
## abSNP: 1666447:1666448
#####################################################################################################################
######################################################################################################################
#################### write all data 


new.maf.DB<-paste(gsub(".txt$","",target),"_maf.txt",sep="")
new.maf.DB.ALL<-paste(gsub(".txt$","",target),"-ALL_maf.txt",sep="")
new.maf.DB
new.maf.DB.ALL
#indels[1:5,core.ann]
#indels[1:5,]
dim(indels)




setwd(annotate.dir) ## Added by Mhairi Oct 2014
write.table(indels[,core.ann],file=target,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE) ## used for annovar
## save(list=c("target","All.grs","snp.grs","indel.grs"),file=paste(project.name,"_indel_files.RData",sep="")) # if sample.gr is a Grange object
save(list=c("target","indels","project.name"),file=paste(gsub(".txt$","",target),"_indel_files.RData",sep="")) # if sample.gr is a Grange object
## load( paste(gsub(".txt$","",target),"_indel_files.RData",sep="") )

if("AF" %in% colnames(indels)){
filter.vals<-(tapply(indels[,"FILTER"],indels[,"FILTER"],length))
if("PASS" %in% names(filter.vals)){
  passing<-indels[,"FILTER"]=="PASS"  
  write.table(indels[passing,c("chr","start","end","REF","ALT","AF","FILTER")],file=new.maf.DB,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
  write.table(indels[,c("chr","start","end","REF","ALT","AF","FILTER")],file=new.maf.DB.ALL,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
  }else{
  write.table(indels[,c("chr","start","end","REF","ALT","AF","GENO")],file=new.maf.DB,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
}}

} # end else update.annovar.annotations

  ## write.table(control.DB,file="Harwell-Control_maf.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

## "annotate_variation.pl -batchsize 20m -regionanno -colsWanted 1,2,3,4  -buildver hg19 -dbtype Gerp Eu_Locus.txt  /media/Bioinform-D/Research/annovar/humandb/hg19 --outfile Eu_Locus.txt.Gerp "

######################################################################################################################
######################################################################################################################
######################################################################################################################
################################## RUN ANNOVAR        ######################################## update.annovar.annotations<-FALSE
setwd(annotate.dir)

files.in.annotate<-dir(annotate.dir)
  
#update.annovar.annotations<-FALSE
############### region_annotation
## regionanno.DB<-c("segdup","dgv","tfbsConsSites","gwasCatalog","omimGene","snp131")
print("regionanno.DB")
print(regionanno.DB)
for(i in 1:length(regionanno.DB)){
if( (paste(target,regionanno.DB[i],"log",sep=".") %in% files.in.annotate ) & !update.annovar.annotations){next} # perhaps have added batabases but dones want to update
get.cols<-toString(eval(as.name(paste(regionanno.DB[i],"labels.posns",sep="."))))
get.cols<-gsub(" ","",get.cols)
print(get.cols)
system(
   paste("annotate_variation.pl -regionanno -colsWanted ",get.cols,"  -buildver ",genome.build," -dbtype ",regionanno.DB[i]," ",target,"  ",anno.DB.location," --outfile ",paste(target,regionanno.DB[i],sep=".")," ",sep="" )
      )
}


if(exists("regionanno.DB.score")){
for(i in 1:length(regionanno.DB.score)){
if( (paste(target,regionanno.DB.score[i],"log",sep=".") %in% files.in.annotate ) & !update.annovar.annotations){next}
get.cols<-toString(eval(as.name(paste(regionanno.DB.score[i],"labels.posns",sep="."))))
get.cols<-gsub(" ","",get.cols)
print(get.cols)
system(
   paste("annotate_variation.pl -regionanno -colsWanted ",get.cols,"  -buildver ",genome.build," -dbtype ",regionanno.DB.score[i]," ",target,"  ",anno.DB.location," -scorecolumn ", regionanno.DB.score.col[i], " --outfile ",paste(target,regionanno.DB.score[i],sep=".")," ",sep="" )
      )
}
}
## annotate_variation.pl -batchsize 20m -regionanno -colsWanted 5  -buildver hg19 -dbtype omimGene check.txt  /media/Bioinform-D/Research/annovar/humandb/ --outfile check.txt.omimGene
############################################################

################### filter known
### WARNING snp132 only returns the snp name AND NOT the maf so it will only annotate in this mode:

## filter.DB<-c("1000g2010nov_all","snp132") ## here dbsnp132  is just annoting not doing filtering  ## 1000g2010nov_all==hg19_ALL.sites.2010_11.txt
for(i in 1:length(filter.DB)){
if( (paste(target,filter.DB[i],"log",sep=".") %in% files.in.annotate ) & !update.annovar.annotations){next} # perhaps have added batabases but dones want to update
system(
    paste("annotate_variation.pl -batchsize 50m -filter --score_threshold ",maf.threshold," -buildver ",genome.build," -dbtype ",filter.DB[i]," ",target,"  ",anno.DB.location," --outfile ",paste(target,filter.DB[i],sep=".")," ",sep="" )
      )
}
############################################################

##FMHT stopped here
################ filter generic with maf
## generic.filter.DB<-c("hg19_snp132_maf.txt","hg19_cg46.txt")  ## "dropped" have maf>= maf.threshold    "filtered have maf< maf.threshold  $$OR$$  if (maf==0)  ## "dropped" are FOUND    "filtered"  are NotFOUND
for(i in 1:length(generic.filter.DB)){
#  if(is.null(generic.filter.DB[i])){next}
#  print(i)
# for(i in 12:13){
if( (paste(target,generic.filter.DB[i],"log",sep=".") %in% files.in.annotate ) & !update.annovar.annotations){next} # perhaps have added batabases but dones want to update
system(
    paste("annotate_variation.pl -batchsize 1000m -filter --score_threshold ",maf.threshold," -buildver ",genome.build," -dbtype  generic -genericdbfile ",generic.filter.DB[i]," ",target,"  ",anno.DB.location," --outfile ",paste(target,generic.filter.DB[i],sep=".")," ",sep="" )
      )
}

############################################################

################ filter for function   ## "dropped" are FOUND    "filtered"  are NotFOUND
## function.filter.DB<-c("ljb_pp2","ljb_mt","ljb_phylop","avsift")
for(i in 1:length(function.filter.DB)){
if(is.null(function.filter.DB[i])){next}
if( (paste(target,function.filter.DB[i],"log",sep=".") %in% files.in.annotate ) & !update.annovar.annotations){next} # perhaps have added batabases but dones want to update
system(
    paste("annotate_variation.pl -batchsize 100m -filter  -buildver ",genome.build," -dbtype ",function.filter.DB[i]," ",target,"  ",anno.DB.location," --outfile ",paste(target,function.filter.DB[i],sep=".")," ",sep="" )
      )
}

## annotate_variation.pl -batchsize 20m -filter --separate  -buildver hg19  -dbtype ljb_pp2 check.txt  /media/Bioinform-D/Research/annovar/humandb/ --outfile test.txt.hg19_ljb_pp2.txt
############################################################

################ gene annotation ensGene knownGene refGene avaible (gene is the same as RefGene)
## geneanno.DB<-c("refGene","knownGene")
for(i in 1:length(geneanno.DB)){
if( (paste(target,geneanno.DB[i],"log",sep=".") %in% files.in.annotate ) & !update.annovar.annotations){next} # perhaps have added batabases but dones want to update
system(
    paste("annotate_variation.pl -batchsize 50m -geneanno --splicing_threshold ",splice.threshold," -buildver ",genome.build," -dbtype ",geneanno.DB[i]," ",target,"  ",anno.DB.location," --outfile ",paste(target,geneanno.DB[i],sep=".")," ",sep="" )
      )
}

if(grepl("^chrMT",target) | grepl(".chrMT.",target)){  ## mitochrondrial using 1000 genomes FASTA - which we use replaces the ensGene one written above for mitochromdria (only works if doing per chromosome
 print("running mitochrondria")
  if( (paste(target,geneanno.DB[i],"log",sep=".") %in% files.in.annotate ) & !update.annovar.annotations){next} # perhaps have added batabases but dones want to update
system(
    paste("annotate_variation.pl -batchsize 50m -geneanno --splicing_threshold ",splice.threshold," -buildver ",genome.build, "-dbtype MT_GRCh37_ensGene ",target,"  ",anno.DB.location," --outfile ",paste(target,"ensGene",sep=".")," ",sep="" )
      )
}
#annotate_variation.pl -buildver hg19 -dbtype MT_GRCh37_ensGene chrMT.BTCK.2.txt /media/Bioinform-D/Research/annovar/humandb/hg19
print("got past annovar")
############################################################
#  update.annovar.annotations<-TRUE
############################################

## /usr/local/bin/annotate_variation.pl -batchsize 100m -filter --score_threshold 0 -buildver hg19 -dbtype generic -genericdbfile hg19_Gerp_Annovar_2.0.txt QBI_All.txt /media/Bioinform-D/Research/annovar/humandb/hg19 --outfile QBI_All.txt.hg19_Gerp_Annovar_2.0.txt
##########################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################

##################################### process all the files that have been written #################################
###### all the regionanno                                       ; files need to be parsed and broken up
###### all the  filter.DB generic.filter.DB function.filter.DB   have just the score only need to parse the "_dropped" with maf set to 0 "filtered contains novel cases"
###### geneanno  have 2 files => .exonic_variant_function (has protein change  **  line1	nonsynonymous SNV	OR4F5:NM_001005484:exon1:c.A421G:p.T141A,
######                           .variant_function   (has location and gene  **    exonic	OR4F5




} ### skip annovar


##################### some lines sin case want to do some developement work
## setwd(annotate.dir)
## load(paste(target,"_indel_files.RData",sep=""))
##   if(!grepl("^chr",indels[1,"chr"])){
##     key.indels<-build.key(indels,core.ann,add.chr.label=TRUE)
##     indels[,"chr"]<-gsub("chr","",indels[,"chr"])
##   }else{key.indels<-build.key(indels,core.ann)      }





{ # this is a geneal bracket that matches with second last line below.

  
if(exists("regionanno.DB.score")){
  regionanno.DB<-c(regionanno.DB,regionanno.DB.score[!(regionanno.DB.score %in% regionanno.DB)]) # may already have sone one step of loop)    
}

setwd(annotate.dir)

if(!grepl("^chr",indels[1,"chr"])){  ### just in case restarted from a strange place
key.indels<-build.key(indels,core.ann,add.chr.label=TRUE)
}else{key.indels<-build.key(indels,core.ann)      }


options(show.error.messages = TRUE)
the.files<-dir(annotate.dir)
the.files<-the.files[!grepl("\\.log$",the.files) & grepl(paste("^",target,sep=""),the.files)]
the.files # will contain that will be process and some extra
## indels<-read.delim(target,header=F,skip=0,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
## regionanno.DB<-c("segdup","dgv","tfbsConsSites","gwasCatalog","omimGene","snp131")

setwd(annotate.dir)
ann.table<-data.frame(key=key.indels,stringsAsFactors=FALSE)
ann.order<-{}

# wait here for SKDP
############################################## regionanno.DB
for(i in 1:length(regionanno.DB)){
the.DB<-regionanno.DB[i]
if(grepl("+",the.DB,fixed=TRUE)){the.DB<-gsub("+","\\+",the.DB,fixed=TRUE)}

data.file<-the.files[grepl(paste("^",target,".",the.DB,".",genome.build,sep=""),the.files)]
####### regionanno do not contain dropped or filtered  could have annotated with same dataset I filtered with

data.file<-data.file[ !(grepl("_dropped$",data.file) | grepl("_filtered$",data.file))]

if(grepl("\\+",the.DB,fixed=TRUE)){the.DB<-gsub("\\+","",the.DB,fixed=TRUE)}

print(paste(the.DB,data.file,sep=" -> "))
the.labels<-eval(as.name(paste(the.DB,"labels",sep=".")))
the.labels.wanted<-eval(as.name(paste(the.DB,"labels.wanted",sep=".")))
  
print("Getting data.file")
the.data<-try(read.delim(data.file,header=F,skip=0,sep="\t",fill=TRUE,stringsAsFactors=FALSE,colClasses = "character"),silent=TRUE)
#inherits(the.data, "try-error")

if(!inherits(the.data, "try-error")){

colnames(the.data)<-c("DB","ann",core.ann)
#the.key<-paste(the.data[,"chr"],the.data[,"start"],the.data[,"end"],the.data[,"REF"],the.data[,"ALT"],the.data[,"TYPE"],sep=":")

the.key<-build.key(the.data,core.ann,add.chr.label=TRUE) #sep is always ":"
if(length(the.key)!=length(unique(the.key))){
  print("WARNING duplicate entries found") ## Hapmap / 1000 Genomes an even within dbSNP can have duplicates key above shoudl work but if not
  posns<-match(unique(the.key),the.key)
  missing<-is.na(posns)
  the.data<-the.data[posns[!missing],]
  the.key<-the.key[posns[!missing]]
}

rownames(the.data)<-the.key

the.data[1:5,]
length(the.key)>0

posns<-match(key.indels,the.key)
missing<-is.na(posns) # 
sum(missing)

     }else{missing<-rep(TRUE,times=length(key.indels))}## no data

the.ann.label<-paste(the.DB,toString(the.labels.wanted),sep="::")
inherits(the.data, "try-error")

if(!inherits(the.data, "try-error")){
  
   ann.table[!missing,the.ann.label]<-the.data[posns[!missing],"ann"]
   
      }else {ann.table[!missing,the.ann.label]<-as.character(NA)}

ann.table[,the.DB]<-!missing
ann.table[1:5,]
ann.order<-c(ann.order,the.DB)


  }
ann.table[1,]
regionanno.table<-ann.table
regionanno.ann.order<-ann.order
rownames(regionanno.table)<-regionanno.table[,"key"]
regionanno.table<-regionanno.table[,-1]
#regionanno.table[1,]
regionanno.table<-regionanno.table[,c(ann.order,colnames(regionanno.table)[!(colnames(regionanno.table) %in% ann.order)])]
#regionanno.table[1:5,]
print(paste("regionanno.table id size:",dim(regionanno.table),sep=""))
## regionanno.table[1:50,"omimGene::name"]


########################################################
########################################################
########################### FILTER DB ##################


  
######################## filter Dbs these have a _dropped ####
ann.table<-data.frame(key=key.indels,stringsAsFactors=FALSE)
ann.order<-{}

#the.DBs<-c(filter.DB,generic.filter.DB,function.filter.DB)
the.DBs<-the.combined.DBs
for(i in 1:length(the.DBs)){
the.DB<-the.DBs[i]

if(grepl("+",the.DB,fixed=TRUE)){the.DB<-gsub("+","\\+",the.DB,fixed=TRUE)}

target.string<-paste("^",target,".",the.DB,"\\S*_dropped$",sep="")
#if(grepl("++",target.string,fixed=TRUE)){target.string<-gsub("++","\\++",target.string,fixed=TRUE)} ##gerp++ casles problems in grep (++)

data.file<-the.files[grepl(target.string,the.files)]

if(grepl("\\+",the.DB,fixed=TRUE)){the.DB<-gsub("\\+","",the.DB,fixed=TRUE)}
print(paste(the.DB,data.file,sep=" -> "))

the.labels<-eval(as.name(paste(the.DB,"labels",sep=".")))
the.labels.wanted<-eval(as.name(paste(the.DB,"labels.wanted",sep=".")))
  
print("also getting data file")
the.data<-try(read.delim(data.file,header=F,skip=0,sep="\t",fill=TRUE,stringsAsFactors=FALSE,colClasses = "character"),silent=TRUE) 
if(!inherits(the.data, "try-error")){

colnames(the.data)<-the.labels
#the.key<-paste(the.data[,"chr"],the.data[,"start"],the.data[,"end"],the.data[,"REF"],the.data[,"ALT"],the.data[,"TYPE"],sep=":")

the.key<-build.key(the.data,core.ann,add.chr.label=TRUE) #sep is always ":"
if(length(the.key)!=length(unique(the.key))){
  print("WARNING duplicate entries found") ## Hapmap / 1000 Genomes an even within dbSNP can have duplicates key above shoudl work but if not
  posns<-match(unique(the.key),the.key)

  missing<-is.na(posns)
  the.data<-the.data[posns[!missing],]
  the.key<-the.key[posns[!missing]]
}

rownames(the.data)<-the.key

the.data[1:5,]
length(the.key)>0

posns<-match(key.indels,the.key)
missing<-is.na(posns) # 
sum(missing)

}else{print("WARNING could not find DATABASE OR NO MATCHES");missing<-rep(TRUE,times=length(key.indels))}## no data

 the.ann.label<-paste(the.DB,toString(the.labels.wanted),sep="::")

if(!inherits(the.data, "try-error")){
ann.table[!missing,the.ann.label]<-the.data[posns[!missing],the.labels.wanted]
}else{ann.table[!missing,the.ann.label]<-as.character(NA)}
ann.table[,the.DB]<-!missing
ann.table
ann.order<-c(ann.order,the.DB)
}
#ann.table[10,]
filter.table<-ann.table
filter.ann.order<-ann.order
rownames(filter.table)<-filter.table[,"key"]
filter.table<-filter.table[,-1]
filter.table<-filter.table[,c(ann.order,colnames(filter.table)[!(colnames(filter.table) %in% ann.order)])]

#filter.table[1:2,]

########################################################
########################################################
######################## Do geneanno.DB #### There are 2 databases here
## ann.table<-data.frame(key=key.indels,location=rep(NA,times=length(key.indels)),type=rep(NA,times=length(key.indels)),gene=rep(NA,times=length(key.indels))stringsAsFactors=FALSE)

ann.table<-data.frame(key=key.indels)
core.table<-data.frame(location=rep(NA,times=length(key.indels)),type=rep(NA,times=length(key.indels)),gene=rep(NA,times=length(key.indels)),stringsAsFactors=FALSE)
ann.order<-{}
the.DBs<-c(geneanno.DB)
i<- 3
for(i in 1:length(the.DBs)){
  
the.DB<-the.DBs[i]
print(the.DB)
  
ann.table.titles<-colnames(ann.table)
ann.table<-cbind(ann.table,core.table)
colnames(ann.table)<-c(ann.table.titles,paste(the.DB,colnames(core.table),sep="::"))

data.file.variant_function<-the.files[grepl(paste("^",target,".",the.DB,"\\S*\\.variant_function$",sep=""),the.files)]
data.file.exonic_variant_function<-the.files[grepl(paste("^",target,".",the.DB,"\\S*\\.exonic_variant_function$",sep=""),the.files)]

print("Getting data.file.variant_function")
the.data.variant_function<-try(read.delim(data.file.variant_function,header=F,skip=0,sep="\t",fill=TRUE,stringsAsFactors=FALSE,colClasses = "character"),silent=TRUE)
print("data.file.exonic_variant_function")
print(data.file.exonic_variant_function)
the.data.exonic_variant_function<-try(read.delim(data.file.exonic_variant_function,header=F,skip=0,sep="\t",fill=TRUE,stringsAsFactors=FALSE,colClasses = "character"),silent=TRUE)

####do first one
print("yip")
the.data<-the.data.exonic_variant_function
dim(the.data)
the.labels<-eval(as.name(paste(the.DB,"labels.exonic_variant_function",sep=".")))
the.labels.wanted<-eval(as.name(paste(the.DB,"labels.wanted.exonic_variant_function",sep=".")))
print(the.labels.wanted)
print("yip 2")
if(!inherits(the.data, "try-error")){
colnames(the.data)<-the.labels
#the.key<-paste(the.data[,"chr"],the.data[,"start"],the.data[,"end"],the.data[,"REF"],the.data[,"ALT"],the.data[,"TYPE"],sep=":")
the.key<-build.key(the.data,core.ann,add.chr.label=TRUE) #sep is always ":"
if(length(the.key)!=length(unique(the.key))){
  print("WARNING duplicate entries found") ## Hapmap / 1000 Genomes an even within dbSNP can have duplicates key above shoudl work but if not
  posns<-match(unique(the.key),the.key)
  missing<-is.na(posns)
  the.data<-the.data[posns[!missing],]
  the.key<-the.key[posns[!missing]]
}
print("wowser")
rownames(the.data)<-the.key
print("pah")
posns<-match(key.indels,the.key)
missing<-is.na(posns) # 
## sum(missing)
ann.table[!missing,paste(the.DB,"location",sep="::")]<-the.data[posns[!missing],"location"]
ann.table[!missing,paste(the.DB,"type",sep="::")]<-the.data[posns[!missing],"type"]
print("pah 111")
}## no data does nothing

## ann.table
print("hwew")
the.data<-the.data.variant_function
the.labels<-eval(as.name(paste(the.DB,"labels.variant_function",sep=".")))
the.labels.wanted<-eval(as.name(paste(the.DB,"labels.wanted.variant_function",sep=".")))
print("wagga")
if(!inherits(the.data, "try-error")){
colnames(the.data)<-the.labels
#the.key<-paste(the.data[,"chr"],the.data[,"start"],the.data[,"end"],the.data[,"REF"],the.data[,"ALT"],the.data[,"TYPE"],sep=":")
the.key<-build.key(the.data,core.ann,add.chr.label=TRUE) #sep is always ":"
if(length(the.key)!=length(unique(the.key))){
  print("WARNING duplicate entries found") ## Hapmap / 1000 Genomes an even within dbSNP can have duplicates key above shoudl work but if not
  posns<-match(unique(the.key),the.key)
  missing<-is.na(posns)
  the.data<-the.data[posns[!missing],]
  the.key<-the.key[posns[!missing]]
}

rownames(the.data)<-the.key

print("MEaow")
####FIX the GENE names ##############################
posns<-match(key.indels,the.key)
missing<-is.na(posns) #
ann.table[!missing,paste(the.DB,"gene",sep="::")]<-the.data[posns[!missing],"type"]

#####################################################################################


############## give annotation to non-exonic hits (not already assigned ##################

already.assigned<-!is.na(ann.table[,paste(the.DB,"location",sep="::")]) ## Exclude those already assigned
#sum(already.assigned)
#test<-ann.table[already.assigned,]

posns<-match(key.indels[!already.assigned],the.key) # only check those not already assigned to exons
missing<-is.na(posns) # 
## sum(missing)
#length(ann.table[!already.assigned,][!missing,paste(the.DB,"location",sep="::")]) ; ann.table[!already.assigned,][!missing,][481100,]
#length(the.data[posns[!missing],"location"])  the.data[posns[!missing],][481100,]
print("Gosh")
ann.table[!already.assigned,][!missing,paste(the.DB,"location",sep="::")] <-the.data[posns[!missing],"location"] # location in class of the indel
ann.table[!already.assigned,][!missing,paste(the.DB,"type",sep="::")] <-the.data[posns[!missing],"type"] 
#########################################################################

########################### CLEAN UP THE GENE NAMES  ############
## sum(is.na(ann.table[,paste(the.DB,"location",sep="::")])
## test<-is.na(ann.table[,paste(the.DB,"location",sep="::")])
## tapply(ann.table[,paste(the.DB,"location",sep="::")],ann.table[,paste(the.DB,"location",sep="::")],length)
    
genes.involved<-ann.table[,paste(the.DB,"gene",sep="::")]
genes.involved<-strsplit(genes.involved,split=";") # "PLOD1;PLOD1" or  "CDK11A,CDK11B;CDK11A,CDK11B" etc
multiple.genes<-unlist(lapply(genes.involved,length))
## tapply(multiple.genes,multiple.genes,length)
multiple.genes<-multiple.genes>1    
    
genes.involved<-lapply(genes.involved[multiple.genes],unique)
genes.involved<-collapse.gene.names(genes.involved,"null",delimit.by=";") # single list list name is not defined
ann.table[multiple.genes,paste(the.DB,"gene",sep="::")]<-genes.involved
print("MEaow 222")

#################################################April 4th 2013 decided that I do not need this part #########################
## #the.data[posns[!missing],"type"][1:100]
## ############extra exonic annotions 
## ########## SOME locations might need to be added too synonymous and Splicing...for example the exon annotions only include one possibility
## posns<-match(key.indels[already.assigned],the.key) # only check those not already assigned to exons
## missing<-is.na(posns) #
## # ann.table[already.assigned,][missing,]

## exonic.labels<-ann.table[already.assigned,][!missing,paste(the.DB,"location",sep="::")] ## NEW
## #exonic.labels<-ann.table[posns[!missing],paste(the.DB,"location",sep="::")] ## OLD march2013  Exclude those already assigned length(ann.table[posns[!missing],paste(the.DB,"location",sep="::")])
## #tapply(exonic.labels,exonic.labels,length)
## # sum(is.na(exonic.labels))
## #grep("chr17:48264412:48264412:A:-:indel",ann.table[posns[!missing],"key"]) # ann.table[posns[!missing],][86764,]
## gene.exonic.labels<-the.data[posns[!missing],"location"] # gene.exonic.labels<-the.data[posns[!missing],"location"] 
## #tapply(gene.exonic.labels,gene.exonic.labels,length)
## full.labels<-paste(exonic.labels,gene.exonic.labels,sep=";")


##  # tapply(full.labels,full.labels,length)  # synonymous SNV;exonic;splicing  are ones that I might have missed before
##   ## test<-grep("synonymous SNV;exonic;splicing",full.labels)[1:500]
##   ##  ann.table[already.assigned,][test,]
##   ##  the.data[posns[!missing],][test,]

## ## #this was in original code potential problem is wanted mutation appears after an unwanted mutations miRNA_exonic before synonoymous etc 
## ## ## exonic;splicing	KLHL17;KLHL17 AND   exonic;splicing	FAAH;FAAH  etc not of interest since are captured as nonsynon. and synon. in the exon part
## not.true.splicing<-gene.exonic.labels=="exonic;splicing" & !(grepl("exon",the.data[posns[!missing],"type"])) # does not contain "exon" in "type" 
## # tapply(full.labels,full.labels,length)  # synonymous SNV;exonic;splicing  are ones that I might have missed before
## # tapply(full.labels[!not.true.splicing],full.labels[!not.true.splicing],length)
## ann.table[already.assigned,][!not.true.splicing,paste(the.DB,"location",sep="::")]<-full.labels[!not.true.splicing]


## ### now fix those gene names ( also accounted for in gene.names
## not.true.splicing.gene<-ann.table[already.assigned,][not.true.splicing,paste(the.DB,"gene",sep="::")]
## not.true.splicing.gene<-strsplit(not.true.splicing.gene,split=";") # "PLOD1;PLOD1" or  "CDK11A,CDK11B;CDK11A,CDK11B" etc
## not.true.splicing.gene<-lapply(not.true.splicing.gene,unique)
## not.true.splicing.gene<-collapse.gene.names(not.true.splicing.gene,"null",delimit.by=";") # single list list name is not defined 

## ann.table[already.assigned,][not.true.splicing,paste(the.DB,"gene",sep="::")] <- not.true.splicing.gene  # fix the gene names
## print(paste("TEST for: ",target,"- sum(not.true.splicing) ",sum(not.true.splicing)," length(not.true.splicing.gene): ",length(not.true.splicing.gene),sep=""))
##     ################
## these ones do not annoted the exon component
## cause a probelm cause stuffs up the gene name CDC2L1,CDC2L2;CDC2L1(uc001agw.1:exon6:c.222+1T>C,uc001agy.1:exon5:c.324+1T>C,uc001agx.1:exon5:c.324+1T>C),CDC2L2
###################################

}   ## no data does nothing
## ann.table
print("wanted..")
wanted.muts<-test.wanted.mutation(ann.table[,paste(the.DB,"location",sep="::")],interesting.mutations,delimit.by=";")
print("wanted muts")
ann.table[,the.DB]<-wanted.muts
print("ann")
#tapply(ann.table[,paste(the.DB,"location",sep="::")],ann.table[,paste(the.DB,"location",sep="::")],length)
#tapply(wanted.muts,ann.table[,paste(the.DB,"location",sep="::")],sum)

ann.order<-c(ann.order,the.DB)

print("ann order")

} # loop over refGene and knownGene

  



  print("done refGene * knownGene")
ann.table[1:5,]
geneanno.table<-ann.table
geneanno.ann.order<-ann.order
rownames(geneanno.table)<-geneanno.table[,"key"]
geneanno.table<-geneanno.table[,-1]
geneanno.table[1:5,]
geneanno.table<-geneanno.table[,c(colnames(geneanno.table)[!(colnames(geneanno.table) %in% ann.order)],ann.order)]

## test<-grep("135103311",indels[,"start"])
## indels[test,1:10]
## geneanno.table[test,]


####################################################
## missing annovar annotations  identified  remove and handled below
####################################################

## missing<-rep(FALSE,times=dim(geneanno.table)[1])
## for(i in 1:length(geneanno.DB)){
## a.test<-is.na(geneanno.table[,paste(geneanno.DB[i],"location",sep="::")])
## missing<- missing | a.test
## }
## sum(missing)
## geneanno.table[missing,]
## if(sum(missing)>0){

##   print("Error missing annotations- removing -DANGEROUS")
##   print(key.indels[missing])
##   write.table(indels[missing,],file=paste(gsub(".txt$","",target),"_MISSING_ANNOTATIONS.txt",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



  
##   geneanno.table<-geneanno.table[!missing,]
##   regionanno.table<-regionanno.table[!missing,]
##   filter.table<-filter.table[!missing,]
##    indels<-indels[!missing,]               
##    key.indels<-key.indels[!missing]
##  }

## geneanno.table[1,]
## filter.table[1,]
## regionanno.table[1,]
## ##########################################################3


######################################## CLEANUP TIME

################# GET THE GENE NAMES
## gene.names must be for each line, genes can be unordered

##snp132 in filter.table is filter to just in snp132 or noy  snp132::maf ->snp132::rs
target.table<-geneanno.table
target.table[1:5,]
ann.order<-geneanno.ann.order
print("do boolean")
## raname the boolean filter to DB::HIT
for(i in 1:length(ann.order)){
colnames(target.table)<-gsub(paste("^",ann.order[i],"$",sep=""),paste(ann.order[i],"HIT",sep="::"),colnames(target.table))
}
gene.names<-get.gene.names(target.table,paste(ann.order,"gene",sep="::")) # give table and column labels: return cleaned up gene names (returned as a list)
names(gene.names)
gene.names[["knownGene::gene"]][1:5]
length(gene.names[[3]])

############# UPDATE PHOLYPHEN

print("do Polyphen")
functional.data.file<-paste(gsub(".txt$","",target),".FunctionalPredictions.RData",sep="")
print("for func")
if(functional.data.file %in% files.in.annotate & !force.functional.update){ ## usually alyways run so will update poly-genic positions
    print("true ")
    setwd(annotate.dir) # or where ever this file is saved
    print(paste("loading:",functional.data.file))
    load(functional.data.file)
  }else{
    setwd(code.dir)   
      print("source poly_read_update")
    source("pholyphen_read_and_Update.r")
     setwd(annotate.dir) 
    save(list=c("filter.table.pholy","pholy.data","sift.data","regulation.data"),file=paste(gsub(".txt$","",target),".FunctionalPredictions.RData",sep=""))
     }
print("done polyphen")
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

#rm(poly.data)
#rm(sift.data)
#rm(regulation.data)
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
    ## sum((!vep.has.no.value | !vep.has.no.feature))
    ## sum(is.na(vep.has.no.value))
    ## sum(is.na(vep.has.no.feature))

 ##    save(list=c("vep.has.no.value","vep.has.no.feature","gene.names","filter.table.pholy","new.values"),file="fucked.Rdata")
  ##   sum(!vep.has.no.value)
  ## attributes(gene.names)
  ## t1<-gene.names$"ensGene::gene"[!vep.has.no.value][95:100]
  ## t2<-as.list(filter.table.pholy[!vep.has.no.value,"Gene"][95:100])
  ## t3<-mapply(c, gene.names$"ensGene::gene"[180:200], as.list(filter.table.pholy[,"Gene"][180:200]), SIMPLIFY=FALSE)
  ## vep.annotations<-mapply(clean.unique.combine, gene.names$"ensGene::gene"[!vep.has.no.value], as.list(filter.table.pholy[!vep.has.no.value,"Gene"]), SIMPLIFY=FALSE)
    
  if(sum((!vep.has.no.value | !vep.has.no.feature)) != length(new.values)){print("ERROR ERROR dimension failure in vep.annotation.update")}
   vep.annotations<-mapply(clean.unique.combine, gene.names$"ensGene::gene"[(!vep.has.no.value | !vep.has.no.feature)], new.values, SIMPLIFY=FALSE)
    
   gene.names$"ensGene::gene"[(!vep.has.no.value | !vep.has.no.feature)]<-vep.annotations


  ## fix missing annovar annottaion for ensGene
    missing.annovar.anno<-is.na(geneanno.table[!vep.has.no.value,"ensGene::gene"])
    geneanno.table[!vep.has.no.value,][missing.annovar.anno,"ensGene::gene"] <- filter.table.pholy[!vep.has.no.value,][missing.annovar.anno,"Gene"]
    

  ## gene.names$"ensGene::gene"[missing]
  ## vep.annotations<-mapply(clean.unique.combine, gene.names$"ensGene::gene"[grep(TRUE,(!vep.has.no.value | !vep.has.no.feature))], new.values, SIMPLIFY=FALSE)  
  ## gene.names$"ensGene::gene"[grep(TRUE,(!vep.has.no.value | !vep.has.no.feature))]<-vep.annotations  # gene.names$"ensGene::gene"[(!vep.has.no.value | !vep.has.no.feature)][1:10]
    
  } ## end update.annotations.with.vep
print("end update.annotations.with.vep")

  ##### Uncomment below if also want to change the geneanno.table  
  ## test<-unlist(collapse.gene.names(vep.annotations,1,delimit.by=","))
  ## geneanno.table[(!vep.has.no.value | !vep.has.no.feature),"ensGene::gene"]<-paste(geneanno.table[(!vep.has.no.value | !vep.has.no.feature),"ensGene::gene"],test,sep=",")

#################################

#################################################################################
###geneanno cleanup  gene.names contains the best guess gene names
##################################################################################

####################################################
## missing annovar annotations  identified  remove and handled below
####################################################

for(i in 1:length(geneanno.DB)){
no.annovar<-is.na(geneanno.table[,paste(geneanno.DB[i],"gene",sep="::")])
## sum(no.annovar)

geneanno.table[no.annovar,paste(geneanno.DB[i],"gene",sep="::")]<-paste("NO_ANNOVAR_",key.indels[no.annovar],sep="")
geneanno.table[no.annovar,paste(geneanno.DB[i],"type",sep="::")]<-paste("NO_ANNOVAR_",key.indels[no.annovar],sep="")
geneanno.table[no.annovar,paste(geneanno.DB[i],"location",sep="::")]<-"unknown"

}
print("done length(geneanno.DB)")
geneanno.table[1,]
filter.table[1,]
regionanno.table[1,]
##########################################################3









#################### GET ADDITIONAL ANNOTATIONS WITH BIOMART
##get gene information

#gene.desc <- read.csv("/home/mmarshall/HPC/PaulGenotypingPipeline/martquery_1120061743_599.txt")
#save(list=c("gene.desc"),file="/home/mmarshall/HPC/PaulGenotypingPipeline/gene.desc.RData")
load("/mnt/UQCCG/Software/annovar/humandb/Biomart.gene.desc.RData")
print("loaded biomart")
#require("biomaRt")
## listMarts()  ; listMarts(host="july2009.archive.ensembl.org",archive=TRUE)
## ensembl=useMart("ensembl")  
## listDatasets(ensembl)  mart<-useMart("ensembl_mart_51",host="may2009.archive.ensembl.org",archive=TRUE); listDatasets(mart)
##  listFilters(mart)       mart<-useMart("ensembl_mart_51",host="may2009.archive.ensembl.org",dataset="hsapiens_gene_ensembl",archive=TRUE)

  
#if(!exists("mart")){



# if(genome.build=="hg18"){mart<-useMart("ensembl_mart_51",host="may2009.archive.ensembl.org",dataset="hsapiens_gene_ensembl",archive=TRUE)
#     }else if(genome.build=="hg19"){
#       mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
#     }else if(genome.build=="mm9"){
#       mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl");print("here:")
#  }
#} # is alread have mart then skip
  
#mart



#genes<-unique(unlist(gene.names[["refGene::gene"]]))

#genes<-unique(unlist(gene.names[["ensGene::gene"]]))
#genes<-genes[!grepl("NONE",genes) & !grepl("^PLACE",genes)]  #### NONE(dist=NONE exists)

#genes[1:10] # !grepl("^PLACE",genes)[1:10]
#USE the two beow for annotion with REFGENE
#gene.ann.wanted<-c("hgnc_symbol","gene_biotype","description")
#gene.desc<-getBM(attributes = gene.ann.wanted, filters = "hgnc_symbol", values=genes, mart = mart)


#USE the two beow for annotion with ENSGENE
# if(genome.build=="hg19"){
# gene.ann.wanted<-c("ensembl_gene_id","hgnc_symbol","gene_biotype","description")
# }
# 
# if(genome.build=="hg18"){
# gene.ann.wanted<-c("ensembl_gene_id","hgnc_symbol","biotype","description")
# }
#   
# if(genome.build=="mm9"){
# gene.ann.wanted<-c("ensembl_gene_id","mgi_symbol","gene_biotype","description")
# }
# if(length(genes)>0){ 
# gene.desc<-getBM(attributes = gene.ann.wanted, filters = "ensembl_gene_id", values=genes, mart = mart)
# }else{gene.desc<-matrix(data=NA,nrow=1,ncol=length(gene.ann.wanted)); colnames(gene.desc)<-gene.ann.wanted}
# #FIX THIS SO DOES THE CHROMOSOME REGION TO GET MiRNA and LOC names name use uscs ids if run with annovar
# ## gene.desc<-getBM( filters = "embl", values=genes, mart = mart)
# ## colnames(gene.desc)<- c("hgnc_symbol","gene_biotype","description")dbass3_name
# gene.desc[1,]

######### get rid of dupicate entries with resplce to the fist column  dim(gene.desc)[1]==length(unique(gene.desc[,1]))
  dups<-duplicated(gene.desc[,1])
  gene.desc<-gene.desc[!dups,]
print("done dup entries")
####################### match genes to another column #########################################

gene.desc.table<-match.cols.and.collapse(gene.names,"ensGene::gene",gene.desc,"Ensembl.Gene.ID",colnames(gene.desc),"delimit")
print("done match genes")
## the above is much faster then below when doing delimited concatination compated to tapply
## list.element.lengths<-unlist(lapply(gene.names[["ensGene::gene"]],length))        
## indel.index<-rep(1:length(list.element.lengths),times=list.element.lengths)
## flat.gene.list<-unlist(the.genes[[annotation.labels[i]]])  
## gene.desc.table<-match.cols.and.collapse.fast(list.element.lengths,indel.index,flat.gene.list,gene.desc,"ensembl_gene_id",gene.ann.wanted,"delimit")


gene.desc.table[1:5,]
########################################################################################


###############
###filter regionanno
############
##snp132 in filter.table is filter to just in snp132 or noy  snp132::maf ->snp132::rs
target.table<-regionanno.table
ann.order<-regionanno.ann.order
print("At  line 1501 get gerp data")
## raname the boolean filter to DB::HIT
for(i in 1:length(ann.order)){
colnames(target.table)<-gsub(paste("^",ann.order[i],"$",sep=""),paste(ann.order[i],"HIT",sep="::"),colnames(target.table))
}
target.table[1,]



  setwd(annotate.dir)
  gerp.data.file<-paste(gsub(".txt$","",target),".GerpPredictions.RData",sep="")
  if(force.GERP.update | !(gerp.data.file %in% files.in.annotate) ){
    if(exists("db.details")){db.details<-list(user="gerp2",pass="GeRpUQ8!6", dbname="gerp", host="di-mysql01.di.uq.edu.au")} #host="ga-apps.di.uq.edu.au"
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
  
rm(gerp.scores)
print("got gerp scores")
## if("omimGene" %in% regionanno.DB){  ## this for the OLD version of omim 
## omim.desc<-try(read.delim(paste(anno.DB.location,"/",genome.build,"_","omimMorbidMap.txt",sep="") ,header=F,skip=0,sep="\t",fill=TRUE,stringsAsFactors=FALSE,colClasses = "character"),silent=TRUE) 
## colnames(omim.desc)<-omim.desc.labels
## omim.desc[1:5,]
## ###  
## if(sum("omimGene::name" %in% colnames(target.table))>0){
## omim.names<-target.table[,"omimGene::name"]
## omim.names<-gsub("Name=","",omim.names)
## omim.names<-(strsplit(omim.names,split=","))

## omim.desc.table<-match.cols.and.collapse(list(omim.names),1,omim.desc,"omim::name",colnames(omim.desc),"delimit") # list is expected to be a list of lists where as omin is just a list
## ## omim.desc.table[omim.desc.table[,"omim::name"]=="NA","omim::name"]<-NA   ## some "NA" rather than NA
## ## omim.desc.table[is.na(omim.desc.table[,"omim::name"]),"omim::name"]

## target.table[,"omimGene::name"]<-apply(omim.desc.table,1,function(x) {if(is.na(x[3])){return("NA")}else{paste(x,collapse="::")}})
## colnames(target.table)[colnames(target.table)=="omimGene::name"]<-paste("omimGene",paste(colnames(omim.desc.table),collapse=","),sep="::")

## }} #"omimGene::name" does not exist so do nothing


target.table[1:5,]
colnames(target.table)
regionanno.table<-target.table
print("done regionano.table")
######################################################################################################################
####################################### ANNOVAR FAILURE ##############################################
##### finished all annotations have tables



geneanno.table[1:2,]
gene.desc.table[1:2,]
regionanno.table[1:2,]
filter.table[1:2,]

# filter all samples at ones takes a lot of memory case cause memory overflows! "quality.cut","quality.type","quality.dirn","quality.thresh",  


######################################################################################################################
######################################################################################################################

############# UPDATE PHOLYPHEN
## print("Do functional annotation")
## functional.data.file<-paste(gsub(".txt$","",target),".FunctionalPredictions.RData",sep="")
## if(functional.data.file %in% files & !redo.pholyphen ){
##     setwd(annotate.dir) # or where ever this file is saved
##     load(functional.data.file)
##   }else{
##     setwd(code.dir)
##     source("pholyphen_read_and_Update.r")
##   }
## print("Finish functional annotation")

############################################################################
###########################################################################

print("Getting omim file")
a.bone.file<-read.delim(paste(anno.DB.location.core,"/","omim.txt",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1,]
the.colnames<-colnames(a.bone.file)
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
a.bone.file<-apply(a.bone.file,1,function(x) if(x[1]=="NA" | is.na(x[1])){NA}else{paste(x,collapse="::")})
a.bone.file[1:5]
omim<-as.data.frame(a.bone.file)
colnames(omim)<-paste("OMIM (",paste(the.colnames,collapse="::"),")",sep="")
omim[1:5,]

################################################# GET OTHER ANOTATIONS WANTED anno.DB.location.core
#if(genome.build=="hg19" |  genome.build=="hg18"){
if(!skip.gene.list.matches){
  print("Do Gene List matches")
#if(!exists("omim")){

#}

#if(!exists("ingenuity.bone.genes")){
  print("Getting bone list")
a.bone.file<-read.delim(paste(anno.DB.location.core,"/","bone.function.list.ingenuity.txt",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit")
ingenuity.bone.genes<-as.data.frame(a.bone.file[,2])
colnames(ingenuity.bone.genes)<-c("ingenuity.bone.genes")
#}

#if(!exists("Dequeant.cycling")){
  print("Dequeant.cycling.txt")
a.bone.file<-read.delim(paste(anno.DB.location.core,"/","Dequeant.cycling.txt",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,"Symbol",c("Symbol","Entrez.Gene.Name"),"delimit")
Dequeant.cycling<-as.data.frame(a.bone.file[,2])
colnames(Dequeant.cycling)<-c("Dequeant.cycling")
#}

  print("mouse.defect")
#if(!exists("mouse.defect")){
a.bone.file<-read.delim(paste(anno.DB.location.core,"/","mouse.defect.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
mouse.defect<-as.data.frame(a.bone.file[,1])
colnames(mouse.defect)<-c("mouse.defect")
#}

#if(!exists("ProtonPump")){
  print("ProtonPump")
a.bone.file<-read.delim(paste(anno.DB.location.core,"/","proton.all.ingenuity.txt",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
ProtonPump<-as.data.frame(a.bone.file[,1])
colnames(ProtonPump)<-c("ProtonPump")
#}

#if(!exists("HypoMag")){
  print("HypoMag")
a.bone.file<-read.delim(paste(anno.DB.location.core,"/","hypomagnesaemia.txt",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
HypoMag<-as.data.frame(a.bone.file[,1])
colnames(HypoMag)<-c("Hypomagnesaemia") #  HypoMag,ProtonPump
#}

  print("sewell.cycling")
#if(!exists("sewell.cycling")){
a.bone.file<-read.delim(paste(anno.DB.location.core,"/","sewell.cycling.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
sewell.cycling<-as.data.frame(a.bone.file[,1])
colnames(sewell.cycling)<-c("sewell.cycling")
#}
  
print("skeletome")
#if(!exists("skeletome")){
a.bone.file<-read.delim(paste(anno.DB.location.core,"/","skeletome.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
a.bone.file<-apply(a.bone.file,1,function(x) if(x[1]=="NA" | is.na(x[1])){NA}else{paste(x,collapse="::")})
skeletome<-as.data.frame(a.bone.file)
colnames(skeletome)<-c("skeletome")
#}
  
  print("LungGenes")
  #if(!exists("skeletome")){
  a.bone.file<-read.delim("/mnt/UQCCG/Sequencing/Projects/IYMD-Lung/TPCH_DRgenesList.csv",header=T,sep=",",fill=TRUE,stringsAsFactors=FALSE)
  a.bone.file[1:5,]
  a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
  a.bone.file<-apply(a.bone.file,1,function(x) if(x[1]=="NA" | is.na(x[1])){NA}else{paste(x,collapse="::")})
  lung.genes<-as.data.frame(a.bone.file)
  colnames(lung.genes)<-c("lung.genes")
  #}

  print("GeneList_NeurodegenAll.csv")
#if(!exists("skeletome")){
a.bone.file<-read.delim(paste(anno.DB.location.core,"/","GeneList_NeurodegenAll.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
a.bone.file<-apply(a.bone.file,1,function(x) if(x[1]=="NA" | is.na(x[1])){NA}else{paste(x,collapse="::")})
Neurodegen<-as.data.frame(a.bone.file)
colnames(Neurodegen)<-c("Neurodegen")
#}

  print("MND_omim.csv")
#if(!exists("skeletome")){
a.bone.file<-read.delim(paste(anno.DB.location.core,"/","MND_omim.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
a.bone.file<-apply(a.bone.file,1,function(x) if(x[1]=="NA" | is.na(x[1])){NA}else{paste(x,collapse="::")})
MND_omim<-as.data.frame(a.bone.file)
colnames(MND_omim)<-c("MND_omim")
#}

  print("skeletome")
#if(!exists("skeletome")){
a.bone.file<-read.delim(paste(anno.DB.location.core,"/","ALS_ji.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
a.bone.file<-apply(a.bone.file,1,function(x) if(x[1]=="NA" | is.na(x[1])){NA}else{paste(x,collapse="::")})
ALS_ji<-as.data.frame(a.bone.file)
colnames(ALS_ji)<-c("ALS_ji")
#}

  
#  /media/UQCCG/Software/annovar/humandb/GeneList_NeurodegenAll.csv
#/media/UQCCG/Software/annovar/humandb/MND_omim.csv
#  /media/UQCCG/Software/annovar/humandb/ALS_ji.csv

  ### TOOK OUT BY Mhairi - April 2014 as broken!! 
# print("ALS Genes")
# a.bone.file<-read.delim(paste(anno.DB.location.core,"/","ALS.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
# a.bone.file[1:5,]
# a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[3],colnames(a.bone.file),"delimit") # macth column must in be collapse list
# a.bone.file<-apply(a.bone.file,1,function(x) if(x[1]=="NA" | is.na(x[1])){NA}else{paste(x,collapse="::")})
# ALS<-as.data.frame(a.bone.file)
# colnames(ALS)<-c("ALS Genes")
  
# skeletome ,mouse.defect,sewell.cycling,Dequeant.cycling,ingenuity.bone.genes

#if(!exists("pot.genes")){
  print("Pot genes")
pot.genes<-read.delim(paste(anno.DB.location.core,"/","potassium.all.ingenuity.txt",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
pot.genes[1:5,]
pot.genes<-match.cols.and.collapse(gene.names,"refGene::gene",pot.genes,"Potassium.Gene",colnames(pot.genes),"delimit")
pot.genes<-as.data.frame(pot.genes[,1])
colnames(pot.genes)<-c("Potassium")

  print("Got past gene lists")

  ## Removed "ALS" as was broken above
   save(list=c("Neurodegen","MND_omim","ALS_ji","omim","ingenuity.bone.genes","Dequeant.cycling","mouse.defect","sewell.cycling","skeletome","lung.genes","pot.genes","HypoMag","ProtonPump"),file=paste(gsub(".txt$","",target),".GeneLists.RData",sep=""))
  #save(list=c("ALS_ji.csv","MND_omim.csv","ALS_ji.csv","omim","ingenuity.bone.genes","Dequeant.cycling","mouse.defect","sewell.cycling","skeletome","lung.genes","pot.genes","HypoMag","ProtonPump"),file=paste(gsub(".txt$","",target),".GeneLists.RData",sep=""))
  
  #########################################################

if(grepl("AML",project)){
  
a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_ALL_Mutations.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
AML_ALL_Mutations<-as.data.frame(a.bone.file[,1]) #
colnames(AML_ALL_Mutations)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)

a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_ASH_Mutations.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
AML_ASH_Mutations<-as.data.frame(a.bone.file[,1])
colnames(AML_ASH_Mutations)<-colnames(a.bone.file)[1] # unique(AML_ASH_Mutations)

a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_CLINICAL_Mutations.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
AML_CLINICAL_Mutations<-as.data.frame(a.bone.file[,1])
colnames(AML_CLINICAL_Mutations)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)

a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_DINDEL_Mutations.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
AML_DINDEL_Mutations<-as.data.frame(a.bone.file[,1])
colnames(AML_DINDEL_Mutations)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)

a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_SNP_Mutations.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
AML_SNP_Mutations<-as.data.frame(a.bone.file[,1])
colnames(AML_SNP_Mutations)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)

a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_DINDEL_Relapse.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
AML_DINDEL_Relapse<-as.data.frame(a.bone.file[,1])
colnames(AML_DINDEL_Relapse)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)

a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_ENU_Mutations.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
AML_ENU_Mutations<-as.data.frame(a.bone.file[,1])
colnames(AML_ENU_Mutations)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)

a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_FRANC_Mutations.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
AML_FRANC_Mutations<-as.data.frame(a.bone.file[,1])
colnames(AML_FRANC_Mutations)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)

a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_LEY_RELAPSE_Mutations.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
AML_LEY_RELAPSE_Mutations<-as.data.frame(a.bone.file[,1])
colnames(AML_LEY_RELAPSE_Mutations)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)

a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_MITO_Mutations.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
AML_MITO_Mutations<-as.data.frame(a.bone.file[,1]) # unique(as.data.frame(a.bone.file[,1]))
colnames(AML_MITO_Mutations)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)

a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_OTHER_Mutations.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
AML_OTHER_Mutations<-as.data.frame(a.bone.file[,1])
colnames(AML_OTHER_Mutations)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)

   save(list=c("Neurodegen","MND_omim","ALS_ji","omim","ingenuity.bone.genes","Dequeant.cycling","mouse.defect","sewell.cycling","lung.genes","skeletome","pot.genes","HypoMag","ProtonPump","AML_OTHER_Mutations","AML_MITO_Mutations","AML_LEY_RELAPSE_Mutations","AML_FRANC_Mutations","AML_ENU_Mutations","AML_DINDEL_Relapse","AML_SNP_Mutations","AML_DINDEL_Mutations","AML_CLINICAL_Mutations","AML_ASH_Mutations","AML_ALL_Mutations"),file=paste(gsub(".txt$","",target),".GeneLists.RData",sep=""))
} ## project AML


} #skip.gene.list.matches
## a.bone.file[ (a.bone.file[,1]!="NA" & !is.na(a.bone.file[,1])),]

## ls()
 ## collapse the gene.names list using refGene

  
## indels<-as.data.frame(All.grs) ### convert genomic ranges to a dta frame
## colnames(indels)[colnames(indels)=="seqnames"]<-"chr"
## indels[,"chr"]<-gsub("chr","",indels[,"chr"])
## key.indels<-paste(indels[,"chr"],indels[,"start"],indels[,"end"],indels[,"REF"],indels[,"ALT"],indels[,"TYPE"],sep=":")

print(" build key..")
if(!grepl("^chr",indels[1,"chr"])){
key.indels<-build.key(indels,core.ann,add.chr.label=TRUE)
}else{key.indels<-build.key(indels,core.ann)      }

print("sample lables")
all.sample.labels<-colnames(indels)
all.sample.labels<-all.sample.labels[grep(".GT$",all.sample.labels)]
all.sample.labels<-gsub(".GT$","",all.sample.labels)
all.sample.labels ### these are ALL the sample labels!
#all.sample.labels<-all.sample.labels[all.sample.labels!="2011-171"] ### one for 2012-07-19_GroupCallGAII

print("do annovar samp labels")
if(length(all.sample.labels)!=0 ){ # if aannovar than these are already calculated and encoded in annnovar
########################### ALL ###############
the.samples<-all.sample.labels
the.samples

genotypes<-as.matrix(indels[,paste(the.samples,"GT",sep=".")])
summary.geno.group<-genotype.summary(genotypes)
colnames(summary.geno.group)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),"ALL",sep=".")
rownames(summary.geno.group)<-key.indels

if( sum(paste(the.samples,"AD",sep=".") %in% colnames(indels))==length(the.samples)){ ### no AD in indels (annovar format maybe)

allele.depths.group<-indels[,paste(the.samples,"AD",sep=".")]
summary.het.indiv<-apply(as.matrix(allele.depths.group),2,allele.summary.individuals)
colnames(summary.het.indiv)<-gsub("$.AD",".HET",colnames(summary.het.indiv))

if( sum(paste(the.samples,"DP",sep=".") %in% colnames(indels))==length(the.samples)){
  allele.depths<-as.matrix(indels[,paste(the.samples,"DP",sep=".")])
   }else{
  allele.depths<-apply(as.matrix(allele.depths.group),2,allele.DP.from.AN)
  colnames(allele.depths)<-gsub(".AD",".DP",colnames(allele.depths))
}

allele.depths.group[genotypes!="0/1"] <- NA
summary.depths.group<-apply(allele.depths.group,1,allele.summary.with.apply)
summary.depths.group<-t(summary.depths.group)
colnames(summary.depths.group)<-paste(c("Hetero.ALT.reads","Hetero.REF.reads","Hetero.Read.Balance"),"ALL",sep=".")
rownames(summary.depths.group)<-key.indels

}else{
  summary.depths.group<-{}
  summary.het.indiv<-{}
  allele.depths<-{}
} ### no DP in indels (annovar format maybe)

## chk<-15437
## chk<-353579
## summary.geno.group[chk,]
## tapply(genotypes[chk,],genotypes[chk,],length)
## summary.depths.group[1:2,]
## summary.het.indiv[1:2,]

summary.geno.group[1:2,]
summary.depths.group[1:2,]
summary.het.indiv[1:2,]
genotypes[1:5,]
allele.depths[1:5,]
###################################################################
}else{ ### no sample labels plink fomrt fam file
summary.geno.group<-{}
summary.depths.group<-{}
summary.het.indiv<-{}
genotypes<-{}
allele.depths<-{}
}

print("ANNOTATION PIPELINE FINISHED......")
#print(all.sample.labels)
 save(list=c("target","project.name","geneanno.DB","indels","key.indels","geneanno.table","gene.desc.table","regionanno.table","filter.table","gene.names","the.combined.DBs","all.sample.labels","interesting.mutations" , "omim","summary.geno.group","summary.depths.group","summary.het.indiv"),file=paste(gsub(".txt$","",target),".table.RData",sep=""))

## GERP scores Functional scores  and gene lists saved above
## load(paste(gsub(".txt$","",target),".table.RData",sep=""))




## save(list=c("project.name","geneanno.DB","All.grs","indels","key.indels","geneanno.table","gene.desc.table","regionanno.table","filter.table","gene.names","maf.threshold.filter","the.combined.DBs","all.filter.cols.novel","all.filter.cols.maf","all.sample.labels" , "omim","ingenuity.bone.genes","Dequeant.cycling","mouse.defect","sewell.cycling","skeletome","pot.genes","interesting.mutations"),file=paste(project.name,".table.RData",sep=""))

} # not a loop bracket just to group output
#} #loop over ichr

 


