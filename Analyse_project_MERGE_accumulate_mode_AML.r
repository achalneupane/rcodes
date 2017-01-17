# 

###############################################

#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/genotyping"
UQCCG.data<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/"
project<-"2013-10-27_AML_with_AOGCControl_NoFailedLane"
project.name <- project
annotate.dir <-paste(UQCCG.data,project,"Annotate",sep="/")
analysis.dir<-paste(UQCCG.data,project,"Analysis",sep="/")

#annotate.dir<-"~/PBSHOME/HiSeq/genotyping/20131115_GroupCallSKDP_MODY_PCC_ChMND//Annotate" # dir(annotate.dir)
#analysis.dir<-"~/PBSHOME/HiSeq/genotyping/20131115_GroupCallSKDP_MODY_PCC_ChMND//Analysis" # dir(analysis.dir)

#project.extension<-".analysis-maf-filtered.txt" # use this one for per-chromosome tests
project.extension<-"analysis.txt" # use this one for compete listing


#".analysis-maf-filtered.txt"   ".analysis-maf-filtered.txt.geno.all.txt"  ## just the exterion not fam.extension! ".analysis-maf-filtered.txt.geno.all.txt" #

fam<-c("TGCM-AML") #  ALL or  all ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-03-19_MND_MODY_LKAF_Nimbelgen/BAM/140206_SN775_0138_AC3751ACXX-SampleSheet_MODY_SJ.csv"
the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-10-27_AML_with_AOGCControl_NoFailedLane/BAM/TGCM-AML_SampleSheet_controls.csv"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)

# 300 samples is too many and exceeds memory limit
# remove.from.controls<-c("Q","AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q") # expand.labels.to.samples(remove.from.controls,control.samples)
# remove.from.all.samples<-c("AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q","Q.MAF:0.001","Q.MAF:0.005","Q.MAF:0") #expand.labels.to.samples(remove.from.all.samples,all.samples)

remove.cols<-c(     ) ### may cause error reload after read in a.indels
#dont.build.summary<-FALSE # FALSE then get a combined result
dont.build.summary<-TRUE # FALSE then get a combined result
accumulate<-TRUE
GQ.thresh<-0 
# # ###############################################


###############################################

#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/genotyping"
UQCCG.data<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/"
project<-"2014-06-24_AML_RSGB_AOGC_withHaplotypeCaller"
project.name <- project
annotate.dir <-paste(UQCCG.data,project,"Annotate",sep="/")
analysis.dir<-paste(UQCCG.data,project,"Analysis",sep="/")


project.extension<-"analysis-maf-filtered.txt" # use this one for compete listing


#".analysis-maf-filtered.txt"   ".analysis-maf-filtered.txt.geno.all.txt"  ## just the exterion not fam.extension! ".analysis-maf-filtered.txt.geno.all.txt" #

fam<-c("ALL") #  ALL or  all ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-03-19_MND_MODY_LKAF_Nimbelgen/BAM/140206_SN775_0138_AC3751ACXX-SampleSheet_MODY_SJ.csv"
the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-06-24_AML_RSGB_AOGC_withHaplotypeCaller/BAM/TGCM-AML_RSGB_PILOT_SampleSheet_controls.csv"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)

# 300 samples is too many and exceeds memory limit
# remove.from.controls<-c("Q","AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q") # expand.labels.to.samples(remove.from.controls,control.samples)
# remove.from.all.samples<-c("AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q","Q.MAF:0.001","Q.MAF:0.005","Q.MAF:0") #expand.labels.to.samples(remove.from.all.samples,all.samples)

remove.cols<-c(     ) ### may cause error reload after read in a.indels
dont.build.summary<-FALSE # FALSE then get a combined result
# # ###############################################

###############################################
## /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-10-09_AML_CompleteGenomics_HC/Analysis/2014-10-09_AML_CompleteGenomics_HC.chr1.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
## /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-10-09_AML_CompleteGenomics_HC/Analysis/2014-10-09_AML_CompleteGenomics_HC.chr2.ALL.ALL_GENOTYPES_analysis.txt
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/genotyping"
UQCCG.data<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/"
project<-"2014-10-09_AML_CompleteGenomics_HC"
project.name <- project
annotate.dir <-paste(UQCCG.data,project,"Annotate",sep="/")
analysis.dir<-paste(UQCCG.data,project,"Analysis",sep="/")


project.extension<-"analysis.txt" # use this one for compete listing

#project.extension<-"analysis-maf-filtered.txt"


#".analysis-maf-filtered.txt"   ".analysis-maf-filtered.txt.geno.all.txt"  ## just the exterion not fam.extension! ".analysis-maf-filtered.txt.geno.all.txt" #

fam<-c(".ALL.ALL_GENOTYPES_") #  ALL or  all ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-03-19_MND_MODY_LKAF_Nimbelgen/BAM/140206_SN775_0138_AC3751ACXX-SampleSheet_MODY_SJ.csv"
the.sample.sheet<-"no_file" # Use no_file to just use all in the genotyping run "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-06-24_AML_RSGB_AOGC_withHaplotypeCaller/BAM/TGCM-AML_RSGB_PILOT_SampleSheet_controls.csv"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)

# 300 samples is too many and exceeds memory limit
# remove.from.controls<-c("Q","AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q") # expand.labels.to.samples(remove.from.controls,control.samples)
# remove.from.all.samples<-c("AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q","Q.MAF:0.001","Q.MAF:0.005","Q.MAF:0") #expand.labels.to.samples(remove.from.all.samples,all.samples)

remove.cols<-c(     ) ### may cause error reload after read in a.indels
dont.build.summary<-TRUE # FALSE then get a combined result
accumulate<-TRUE
GQ.thresh<-0 ### this will only write FILTER=PASS
# # ###############################################


###############################################
## /media/UQCCG/Sequencing/Data/Sequence_Genotypes/Data_Freeze_ChNMD_AOGC_Tulane/Analysis/Data_Freeze_ChNMD_AOGC_Tulane.chr1.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
## /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-10-09_AML_CompleteGenomics_HC/Analysis/2014-10-09_AML_CompleteGenomics_HC.chr2.ALL.ALL_GENOTYPES_analysis.txt
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/genotyping"
UQCCG.data<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/"
project<-"Data_Freeze_ChNMD_AOGC_Tulane"
project.name <- project
annotate.dir <-paste(UQCCG.data,project,"Annotate",sep="/")
analysis.dir<-paste(UQCCG.data,project,"Analysis",sep="/")


#project.extension<-"analysis.txt" # use this one for compete listing

project.extension<-"analysis-maf-filtered.txt"


#".analysis-maf-filtered.txt"   ".analysis-maf-filtered.txt.geno.all.txt"  ## just the exterion not fam.extension! ".analysis-maf-filtered.txt.geno.all.txt" #

fam<-c(".ALL.ALL_GENOTYPES_") #  ALL or  all ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-03-19_MND_MODY_LKAF_Nimbelgen/BAM/140206_SN775_0138_AC3751ACXX-SampleSheet_MODY_SJ.csv"
the.sample.sheet<-"no_file" # Use no_file to just use all in the genotyping run "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-06-24_AML_RSGB_AOGC_withHaplotypeCaller/BAM/TGCM-AML_RSGB_PILOT_SampleSheet_controls.csv"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)

# 300 samples is too many and exceeds memory limit
# remove.from.controls<-c("Q","AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q") # expand.labels.to.samples(remove.from.controls,control.samples)
# remove.from.all.samples<-c("AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q","Q.MAF:0.001","Q.MAF:0.005","Q.MAF:0") #expand.labels.to.samples(remove.from.all.samples,all.samples)

remove.cols<-c(     ) ### may cause error reload after read in a.indels
dont.build.summary<-TRUE # FALSE then get a combined result
accumulate<-TRUE
GQ.thresh<-0 ### this will only write FILTER=PASS
# # ###############################################



###############################################

##/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/2013-10-27_AML_with_AOGCControl_NoFailedLane.chrY.ALL.ALL_GENOTYPES_WITH_SYNON_analysis-maf-filtered.txt
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/genotyping"
UQCCG.data<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/"
project<-"2013-10-27_AML_with_AOGCControl_NoFailedLane"
project.name <- project
annotate.dir <-paste(UQCCG.data,project,"Annotate",sep="/")
analysis.dir<-paste(UQCCG.data,project,"Analysis",sep="/")


#project.extension<-"analysis.txt" # use this one for compete listing

project.extension<-"analysis-maf-filtered.txt"
#project.extension<-"analysis.txt"

#".analysis-maf-filtered.txt"   ".analysis-maf-filtered.txt.geno.all.txt"  ## just the exterion not fam.extension! ".analysis-maf-filtered.txt.geno.all.txt" #

fam<-c(".ALL.ALL_GENOTYPES_WITH_SYNON_") #  ALL or  all ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-03-19_MND_MODY_LKAF_Nimbelgen/BAM/140206_SN775_0138_AC3751ACXX-SampleSheet_MODY_SJ.csv"
the.sample.sheet<-"no_file" # Use no_file to just use all in the genotyping run "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-06-24_AML_RSGB_AOGC_withHaplotypeCaller/BAM/TGCM-AML_RSGB_PILOT_SampleSheet_controls.csv"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)

# 300 samples is too many and exceeds memory limit
# remove.from.controls<-c("Q","AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q") # expand.labels.to.samples(remove.from.controls,control.samples)
# remove.from.all.samples<-c("AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q","Q.MAF:0.001","Q.MAF:0.005","Q.MAF:0") #expand.labels.to.samples(remove.from.all.samples,all.samples)

remove.cols<-c(     ) ### may cause error reload after read in a.indels
dont.build.summary<-TRUE # FALSE then get a combined result
accumulate<-TRUE
GQ.thresh<-0 ### this will only write FILTER=PASS
# # ###############################################


###############################################
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-11-04_AML_TCGA_Replication
##/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/2013-10-27_AML_with_AOGCControl_NoFailedLane.chrY.ALL.ALL_GENOTYPES_WITH_SYNON_analysis-maf-filtered.txt
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/genotyping"
UQCCG.data<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014/"
project<-"2014-11-04_AML_TCGA_Replication"
project.name <- project
annotate.dir <-paste(UQCCG.data,project,"Annotate",sep="/")
analysis.dir<-paste(UQCCG.data,project,"Analysis",sep="/")


#project.extension<-"analysis.txt" # use this one for compete listing

project.extension<-"analysis-maf-filtered.txt"
#project.extension<-"analysis.txt"
#".analysis-maf-filtered.txt"   ".analysis-maf-filtered.txt.geno.all.txt"  ## just the exterion not fam.extension! ".analysis-maf-filtered.txt.geno.all.txt" #
fam<-c(".ALL.ALL_GENOTYPES_") #  ALL or  all ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-03-19_MND_MODY_LKAF_Nimbelgen/BAM/140206_SN775_0138_AC3751ACXX-SampleSheet_MODY_SJ.csv"
the.sample.sheet<-"no_file" # Use no_file to just use all in the genotyping run "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-06-24_AML_RSGB_AOGC_withHaplotypeCaller/BAM/TGCM-AML_RSGB_PILOT_SampleSheet_controls.csv"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)

# 300 samples is too many and exceeds memory limit
# remove.from.controls<-c("Q","AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q") # expand.labels.to.samples(remove.from.controls,control.samples)
# remove.from.all.samples<-c("AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q","Q.MAF:0.001","Q.MAF:0.005","Q.MAF:0") #expand.labels.to.samples(remove.from.all.samples,all.samples)

remove.cols<-c(     ) ### may cause error reload after read in a.indels
dont.build.summary<-TRUE # FALSE then get a combined result
accumulate<-TRUE
GQ.thresh<-0 ### this will only write FILTER=PASS
# # ###############################################

###############################################
 #/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-01-15_LeoPharma_Dec2014Freeze/Analysis/2015-01-15_LeoPharma_Dec2014Freeze.chr19.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/genotyping"
UQCCG.data<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/"
project<-"2015-01-15_LeoPharma_Dec2014Freeze"
project.name <- project
annotate.dir <-paste(UQCCG.data,project,"Annotate",sep="/")
analysis.dir<-paste(UQCCG.data,project,"Analysis",sep="/")


#project.extension<-"analysis.txt" # use this one for compete listing
project.extension<-"variants.analysis.txt" ; fam<-c(".large.")
#project.extension<-"analysis-maf-filtered.txt"
#project.extension<-"analysis.txt"
#".analysis-maf-filtered.txt"   ".analysis-maf-filtered.txt.geno.all.txt"  ## just the exterion not fam.extension! ".analysis-maf-filtered.txt.geno.all.txt" #
fam<-c(".ALL.ALL_GENOTYPES_") #  ALL or  all ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-03-19_MND_MODY_LKAF_Nimbelgen/BAM/140206_SN775_0138_AC3751ACXX-SampleSheet_MODY_SJ.csv"
the.sample.sheet<-"no_file" # Use no_file to just use all in the genotyping run "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-06-24_AML_RSGB_AOGC_withHaplotypeCaller/BAM/TGCM-AML_RSGB_PILOT_SampleSheet_controls.csv"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)

# 300 samples is too many and exceeds memory limit
# remove.from.controls<-c("Q","AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q") # expand.labels.to.samples(remove.from.controls,control.samples)
# remove.from.all.samples<-c("AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q","Q.MAF:0.001","Q.MAF:0.005","Q.MAF:0") #expand.labels.to.samples(remove.from.all.samples,all.samples)

remove.cols<-c(     ) ### may cause error reload after read in a.indels
dont.build.summary<-TRUE # FALSE then get a combined result
accumulate<-TRUE
GQ.thresh<-0 ### this will only write FILTER=PASS
# # ###############################################




###############################################
# /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/2015-02-11_LeoPharmFinalAnalysis.chr1.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/genotyping"
UQCCG.data<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/"
project<-"2015-02-11_LeoPharmFinalAnalysis"
project.name <- project
annotate.dir <-paste(UQCCG.data,project,"Annotate",sep="/")
analysis.dir<-paste(UQCCG.data,project,"Analysis",sep="/")


#project.extension<-"analysis.txt" # use this one for compete listing
#project.extension<-"variants.analysis.txt" ; fam<-c(".large.")
project.extension<-"analysis-maf-filtered.txt"
#project.extension<-"analysis.txt"
#".analysis-maf-filtered.txt"   ".analysis-maf-filtered.txt.geno.all.txt"  ## just the exterion not fam.extension! ".analysis-maf-filtered.txt.geno.all.txt" #
fam<-c(".ALL.ALL_GENOTYPES_") #  ALL or  all ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-03-19_MND_MODY_LKAF_Nimbelgen/BAM/140206_SN775_0138_AC3751ACXX-SampleSheet_MODY_SJ.csv"
the.sample.sheet<-"no_file" # Use no_file to just use all in the genotyping run "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-06-24_AML_RSGB_AOGC_withHaplotypeCaller/BAM/TGCM-AML_RSGB_PILOT_SampleSheet_controls.csv"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)

# 300 samples is too many and exceeds memory limit
# remove.from.controls<-c("Q","AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q") # expand.labels.to.samples(remove.from.controls,control.samples)
# remove.from.all.samples<-c("AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q","Q.MAF:0.001","Q.MAF:0.005","Q.MAF:0") #expand.labels.to.samples(remove.from.all.samples,all.samples)

remove.cols<-c(     ) ### may cause error reload after read in a.indels
dont.build.summary<-TRUE # FALSE then get a combined result
accumulate<-TRUE
GQ.thresh<-0 ### this will only write FILTER=PASS
# # ###############################################


###############################################
# /media/UQCCG/Sequencing/Projects/BONE-GENOME/Analysis/mergeCHR2update.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/genotyping"
UQCCG.data<-"/media/UQCCG/Sequencing/Projects/"
project<-"BONE-GENOME"
project.name <- "mergeCHR"
annotate.dir <-paste(UQCCG.data,project,"Annotate",sep="/")
analysis.dir<-paste(UQCCG.data,project,"Analysis",sep="/")


#project.extension<-"analysis.txt" # use this one for compete listing
#project.extension<-"variants.analysis.txt" ; fam<-c(".large.")
project.extension<-"analysis-maf-filtered.txt"
#project.extension<-"analysis.txt"
#".analysis-maf-filtered.txt"   ".analysis-maf-filtered.txt.geno.all.txt"  ## just the exterion not fam.extension! ".analysis-maf-filtered.txt.geno.all.txt" #
fam<-c("update.ALL.ALL_GENOTYPES_") #  ALL or  all ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-03-19_MND_MODY_LKAF_Nimbelgen/BAM/140206_SN775_0138_AC3751ACXX-SampleSheet_MODY_SJ.csv"
the.sample.sheet<-"no_file" # Use no_file to just use all in the genotyping run "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-06-24_AML_RSGB_AOGC_withHaplotypeCaller/BAM/TGCM-AML_RSGB_PILOT_SampleSheet_controls.csv"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)

# 300 samples is too many and exceeds memory limit
# remove.from.controls<-c("Q","AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q") # expand.labels.to.samples(remove.from.controls,control.samples)
# remove.from.all.samples<-c("AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q","Q.MAF:0.001","Q.MAF:0.005","Q.MAF:0") #expand.labels.to.samples(remove.from.all.samples,all.samples)

remove.cols<-c(     ) ### may cause error reload after read in a.indels
dont.build.summary<-FALSE # FALSE then get a combined result
accumulate<-FALSE
GQ.thresh<-0 ### this will only write FILTER=PASS
# # ###############################################



###############################################
# /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-03_LeoPharma_NovoAlign/Analysis/2015-06-03_LeoPharma_NovoAlign.chr1.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/genotyping"
UQCCG.data<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/"
project<-"2015-06-03_LeoPharma_NovoAlign"
project.name <- project
annotate.dir <-paste(UQCCG.data,project,"Annotate",sep="/")
analysis.dir<-paste(UQCCG.data,project,"Analysis",sep="/")


#project.extension<-"analysis.txt" # use this one for compete listing
#project.extension<-"variants.analysis.txt" ; fam<-c(".large.")
#project.extension<-"analysis-maf-filtered.txt"
project.extension<-"analysis.txt"
#".analysis-maf-filtered.txt"   ".analysis-maf-filtered.txt.geno.all.txt"  ## just the exterion not fam.extension! ".analysis-maf-filtered.txt.geno.all.txt" #
fam<-c(".ALL.ALL_GENOTYPES_") #  ALL or  all ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-03-19_MND_MODY_LKAF_Nimbelgen/BAM/140206_SN775_0138_AC3751ACXX-SampleSheet_MODY_SJ.csv"
the.sample.sheet<-"no_file" # Use no_file to just use all in the genotyping run "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-06-24_AML_RSGB_AOGC_withHaplotypeCaller/BAM/TGCM-AML_RSGB_PILOT_SampleSheet_controls.csv"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)

# 300 samples is too many and exceeds memory limit
# remove.from.controls<-c("Q","AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q") # expand.labels.to.samples(remove.from.controls,control.samples)
# remove.from.all.samples<-c("AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q","Q.MAF:0.001","Q.MAF:0.005","Q.MAF:0") #expand.labels.to.samples(remove.from.all.samples,all.samples)

remove.cols<-c(     ) ### may cause error reload after read in a.indels
dont.build.summary<-TRUE # FALSE then get a combined result
accumulate<-TRUE
GQ.thresh<-0 ### this will only write FILTER=PASS
# # ###############################################


###############################################
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/2015-03-16_AllAMLandLung.chr2.ALL.ALL_GENOTYPES_analysis.txt
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/genotyping"
UQCCG.data<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/"
project<-"2015-03-16_AllAMLandLung"
project.name <- project
annotate.dir <-paste(UQCCG.data,project,"Annotate",sep="/")
analysis.dir<-paste(UQCCG.data,project,"Analysis",sep="/")


#project.extension<-"analysis.txt" # use this one for compete listing
#project.extension<-"variants.analysis.txt" ; fam<-c(".large.")
project.extension<-"analysis-maf-filtered.txt"
#project.extension<-"analysis.txt"
#".analysis-maf-filtered.txt"   ".analysis-maf-filtered.txt.geno.all.txt"  ## just the exterion not fam.extension! ".analysis-maf-filtered.txt.geno.all.txt" #
fam<-c(".ALL.ALL_GENOTYPES_") #  ALL or  all ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-03-19_MND_MODY_LKAF_Nimbelgen/BAM/140206_SN775_0138_AC3751ACXX-SampleSheet_MODY_SJ.csv"
the.sample.sheet<-"no_file" # Use no_file to just use all in the genotyping run "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-06-24_AML_RSGB_AOGC_withHaplotypeCaller/BAM/TGCM-AML_RSGB_PILOT_SampleSheet_controls.csv"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)

related.or.bad.file<- "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/related.or.bad.txt"
related.or.bad<-read.table(related.or.bad.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)


remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)

# 300 samples is too many and exceeds memory limit
# remove.from.controls<-c("Q","AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q") # expand.labels.to.samples(remove.from.controls,control.samples)
# remove.from.all.samples<-c("AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q","Q.MAF:0.001","Q.MAF:0.005","Q.MAF:0") #expand.labels.to.samples(remove.from.all.samples,all.samples)

remove.cols<-c(     ) ### may cause error reload after read in a.indels
dont.build.summary<-TRUE # FALSE then get a combined result
accumulate<-TRUE
GQ.thresh<-0 ### this will only write FILTER=PASS
# # ###############################################








###############################################
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/SVAnalysis/2015-03-16_AllAMLandLung.chr2.large.variants.analysis.rare.txt
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/genotyping"
UQCCG.data<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/"
project<-"2015-03-16_AllAMLandLung"
project.name <- project
annotate.dir <-paste(UQCCG.data,project,"Annotate",sep="/")
analysis.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/SVAnalysis"


#project.extension<-"analysis.txt" # use this one for compete listing
#project.extension<-"variants.analysis.txt" ; fam<-c(".large.")
#project.extension<-"analysis-maf-filtered.txt"
project.extension<-".large.variants.analysis.rare.txt"
#".analysis-maf-filtered.txt"   ".analysis-maf-filtered.txt.geno.all.txt"  ## just the exterion not fam.extension! ".analysis-maf-filtered.txt.geno.all.txt" #
fam<-c("") #  ALL or  all ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-03-19_MND_MODY_LKAF_Nimbelgen/BAM/140206_SN775_0138_AC3751ACXX-SampleSheet_MODY_SJ.csv"
the.sample.sheet<-"no_file" # Use no_file to just use all in the genotyping run "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-06-24_AML_RSGB_AOGC_withHaplotypeCaller/BAM/TGCM-AML_RSGB_PILOT_SampleSheet_controls.csv"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)

related.or.bad.file<- "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/related.or.bad.txt"
related.or.bad<-read.table(related.or.bad.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)


remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)

# 300 samples is too many and exceeds memory limit
# remove.from.controls<-c("Q","AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q") # expand.labels.to.samples(remove.from.controls,control.samples)
# remove.from.all.samples<-c("AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q","Q.MAF:0.001","Q.MAF:0.005","Q.MAF:0") #expand.labels.to.samples(remove.from.all.samples,all.samples)

remove.cols<-c(     ) ### may cause error reload after read in a.indels
dont.build.summary<-TRUE # FALSE then get a combined result
accumulate<-TRUE
GQ.thresh<-0 ### this will only write FILTER=PASS
# # ###############################################



###############################################
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/2015-03-16_AllAMLandLung.chr2.ALL.ALL_GENOTYPES_analysis.txt
#UQCCG.data<-"/home/mmarshall/PBSHOME/HiSeq/genotyping"
UQCCG.data<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/"
project<-"2015-08-14_AML_mixedAligners"
project.name <- project
annotate.dir <-paste(UQCCG.data,project,"Annotate",sep="/")
analysis.dir<-paste(UQCCG.data,project,"Analysis",sep="/")


#project.extension<-"analysis.txt" # use this one for compete listing
#project.extension<-"variants.analysis.txt" ; fam<-c(".large.")
#project.extension<-"analysis-maf-filtered.txt"
project.extension<-"analysis.txt"
#".analysis-maf-filtered.txt"   ".analysis-maf-filtered.txt.geno.all.txt"  ## just the exterion not fam.extension! ".analysis-maf-filtered.txt.geno.all.txt" #
fam<-c(".ALL.ALL_GENOTYPES_") #  ALL or  all ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-03-19_MND_MODY_LKAF_Nimbelgen/BAM/140206_SN775_0138_AC3751ACXX-SampleSheet_MODY_SJ.csv"
the.sample.sheet<-"no_file" # Use no_file to just use all in the genotyping run "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-06-24_AML_RSGB_AOGC_withHaplotypeCaller/BAM/TGCM-AML_RSGB_PILOT_SampleSheet_controls.csv"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)

related.or.bad.file<- "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/related.or.bad.txt"
related.or.bad<-read.table(related.or.bad.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)


remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)

# 300 samples is too many and exceeds memory limit
# remove.from.controls<-c("Q","AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q") # expand.labels.to.samples(remove.from.controls,control.samples)
# remove.from.all.samples<-c("AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q","Q.MAF:0.001","Q.MAF:0.005","Q.MAF:0") #expand.labels.to.samples(remove.from.all.samples,all.samples)

remove.cols<-c(     ) ### may cause error reload after read in a.indels
dont.build.summary<-TRUE # FALSE then get a combined result
accumulate<-TRUE
GQ.thresh<-0 ### this will only write FILTER=PASS
# # ###############################################



############################################# BEGIN #######################################################
options(stringsASFactors=FALSE)
options(width=250,max.print=4000)
genome.build<-"hg19"
if(!exists("dont.build.summary")){dont.build.summary<-FALSE} ## if TRUE does not make summary indel file
if(!exists("set.unlabelled.to.control")){set.unlabelled.to.control<-TRUE} # set samples not listed in the samplesheet to controls 
######################################################### Set up the basics for each run #####
########################### Set up for annovar #################################     
maf.threshold<-0.0 
if(!exists("GQ.thresh")){GQ.thresh<-0}
anno.DB.location.core<-"/mnt/UQCCG/Software/annovar/humandb"
anno.DB.location<-paste(anno.DB.location.core,genome.build,sep="/")

######################### Info the run
generic.filter.DB<-c("hg19_sequnome.txt")    # returns 2  extra columns (DB,score)
names(generic.filter.DB)<-c("AML-sequnome")
update.annovar.annotations<-TRUE
core.ann<-c("chr","start","end","REF","ALT","TYPE")
##################################
code.dir<-"/mnt/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/"
setwd(code.dir)
source("ucsc.table.names.r")   # load in the UCSC tables these use the db file names and not their lable-names 
source("annotate_SNPs_subroutines.r")
source("ucsc.table.names.processor.r")





#### test fam list
files<-dir(analysis.dir)
the.extension<-paste(project.extension,"$",sep="")
files<-files[grepl(the.extension ,files)]
toString( unique(unlist( mapply(function(x){x[length(x)]}, strsplit(gsub(the.extension,"",files),split=".",fixed=TRUE)   ))) )
toString( files)
####


#############################################################################################################

######################################### Predefined variables required
##################################################################################


#### assume has format project.chr.fam.extension or chr.project.fam.extension
setwd(analysis.dir)

files<-dir(analysis.dir)
the.extension<-paste(project.extension,"$",sep="")
files<-files[grepl(the.extension ,files)]
if(fam=="ALL" | fam=="All" | fam=="all" ){
  fam<-unique(unlist( mapply(function(x){x[length(x)]}, strsplit(gsub(the.extension,"",files),split=".",fixed=TRUE)   )))
}

fam
## ifam<-2

## geno.aogc<-read.delim("/media/scratch2/AOGC-NGS/Analysis/AOGC-Genotyping.output.AOGC_ALL.geno.all.txt",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
## setwd("/media/Bioinform-D/Research/annovar/humandb/")
## save(list=c("geno.aogc","use.key"),file="aogc.count.data.RData")


## load("/mnt/UQCCG/Software/annovar/humandb/aogc.count.data.RData")
## #print(colnames(indels))
## print(colnames(geno.aogc))
## use.key<-build.key(geno.aogc,core.ann)
## insert.location<-70  ### this is where to add AOGC data INSERTS AFTER THE LOCATION
################################ add aogc


#
ifam<-1
for(ifam in 1:length(fam)){
  
  the.extension<-paste(fam[ifam],project.extension,"$",sep="")
  project.files<-files[grepl(the.extension ,files)]
  print(sort(paste("Doing: ",project.files,sep=""))) # project.files<-project.files[1:22]
  
  indels<-{}
  the.col<-{}
  
  
  data.summary<-{} # used if doing genotype counts
  project.files
  # ichr<-1




 
  project.files<-project.files[!grepl("chrALL",project.files)] ## remove old combined runs 


  
  if(length(project.files)!=24){
    print("########################################### WARNING #################################")
    print("less that 24 chromosomes detected")
    print(fam[ifam])
    print("########################################### WARNING #################################") 
  }
  project.files  #  project.files <- project.files[-23]
   ichr<-20
  for(ichr in 1:length(project.files)){
    print(project.files[ichr])
    
    
    
    ################## fast read ###########
    column.labels<-read.delim(project.files[ichr],header=F,nrows=1,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
    num.vars<-dim(column.labels)[2]
    a.indel<-scan(project.files[ichr],what=character(num.vars),skip=1,sep="\t",fill=TRUE,na.strings="",quote="") # quote="\""
    num.lines<-length(a.indel)/(num.vars)
    dim(a.indel)<-c(num.vars,num.lines)
    a.indel<-t(a.indel)
    colnames(a.indel)<-column.labels

    print(table(a.indel[,"chr"]))
    ########################################
    #a.indel<-read.delim(project.files[ichr],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
    all.possible.samples<-gsub(".GT$","",colnames(a.indel)[grep(".GT$",colnames(a.indel))],perl=TRUE)
    print(length( all.possible.samples))
    #if(set.unlabelled.to control){ grep("PL$",column.labels)

  ##   sum(all.possible.samples %in% sample.sheet.full[,"ParticipantCode"])
  ## all.possible.samples[ !(all.possible.samples %in% sample.sheet.full[,"ParticipantCode"])]
## the.chr<-a.indel[1,"chr"]
## print(paste("Doing Chromosome ",the.chr))


## if(!grepl("^chr",the.chr)){
## a.indel[,"chr"]<-paste("chr",a.indel[,"chr"],sep="")
## }
## key<-build.key(a.indel,core.ann)
## rownames(a.indel)<-key

##     mafs<-colnames(a.indel)[grepl("maf",colnames(a.indel))]

## a.indel["chr19:8995964:8995964:G:T:snp",c(core.ann,"refGene::location","knownGene::location",
##               "ensGene::location", "Consequence.Embl","MAF.lt:0","MAF.lt:0.001","MAF.lt:0.005","MAF.lt:0.01","MAF.lt:0.025",mafs )]

    
## a.chk<-c("chr19:8993452:8993452:A:G:snp","chr19:8995964:8995964:G:T:snp")
## a.chk %in% rownames(a.indel)
    
    ## print(dim(a.indel))
    ## a.indel[1:2,1:50]

#    table(a.indel[,"FILTER"])
#    table(a.indel[,"FILTER"])
#    table(a.indel[,"wanted.muts"])



######################### subset samples #################
########################################### samlpe subsetting using sample sheet

    ### no_file just use all samples - all assumed as affected
if(the.sample.sheet=="no.file" | the.sample.sheet=="no_file"){
  all.samples<-all.possible.samples
  sample.sheet.full<-cbind(fam[ifam],fam[ifam], all.samples,".",".","-9",2)
colnames(sample.sheet.full)<-c("SampleProject","FamilyCode","ParticipantCode","PaternalID","MaternalID","Sex","AffectionStatus")
  control.samples<-{} # assume no controls  control.samples<-c()
  

    
  }else{
sample.sheet.full<-read.delim(the.sample.sheet,header=T,sep=",",fill=TRUE,stringsAsFactors=FALSE)
sample.sheet.full[1:5,]
dim(sample.sheet.full)
colnames(sample.sheet.full)
tapply(sample.sheet.full[,"SampleProject"],sample.sheet.full[,"SampleProject"],length)

all.samples<-unique(sample.sheet.full[,"ParticipantCode"])
control.samples<-unique(sample.sheet.full[sample.sheet.full[,"SampleProject"]=="Control" | sample.sheet.full[,"SampleProject"]=="control"  | sample.sheet.full[,"SampleProject"]=="CONTROL" , "ParticipantCode"])
}

sample.sheet.full[1:5,]

#######################identify cols to remove


    
if(!is.null(remove.from.controls) & !is.null(control.samples)){
to.remove.control<-expand.labels.to.samples.complex(remove.from.controls,control.samples,paste.after=FALSE,seperator=".")
to.remove.control.2<-expand.labels.to.samples.complex(remove.from.controls,control.samples,paste.after=FALSE,seperator=":") ## get "S02-F21-P01:0.005(MAF).Q"
to.remove.control<-unique(c(to.remove.control,to.remove.control.2))
}else{
  to.remove.control<-{}
}

if(!is.null(remove.from.all.samples)){
to.remove.all<-expand.labels.to.samples(remove.from.all.samples,all.samples)
}else{
  to.remove.all<-{}
}

to.remove.samples<-unique(c(to.remove.control,to.remove.all))
remove.cols<-unique(c(remove.cols,to.remove.samples))


if(!is.null(remove.cols)){
a.indel<-a.indel[,colnames(a.indel)[!(colnames(a.indel) %in% remove.cols)]]
  }
   ##  pass.qual<-a.indel[,"FILTER"]=="PASS" & as.logical(a.indel[,"wanted.muts"]) 
   ##  print(sum(pass.qual))
   ## a.indel<-a.indel[pass.qual,c(1:51,2062:2131)]
    
    ### below used a one-off prefilter
    ## pass.qual<-a.indel[,"FILTER"]=="PASS" & as.logical(a.indel[,"wanted.muts.coding"]) & as.logical(a.indel[,"MAF.lt:0.001"])
    ## print(sum(pass.qual))
    ## a.indel<-a.indel[pass.qual,c(1:51,2062:2131)]
    
    ##################################ADD extras here
    #setwd(code.dir)
    # source("MERGE_alt_alleles_AOGC.r")
    # source("MERGE_alt_alleles_AOGC_as_controls.r")
    # source("MERGE_genotype.counts.r") ## counts tyes of mutations 
    
   # all.lung.genes <- cbind(a.indel[,"refGene::gene"],a.indel[,"refGene::gene"])
  #  colnames(all.lung.genes) <- c("LungGenes","LungGenes")
    
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
    #print(project.files[ichr])
    #print(dim(a.indel))
    
   # HERE MHAIRI 
    
   

    
  
    
## "FANCM " is "FANCM" (no white space) in one column
## "CHEK1 " is "CHEK1" (no white space)  error only in FA_BRCA_HRR
## DNMT3a  is  DNMT3A
## FANCD1   is BRCA2
## CTIP is RBBP8
## "SLX1"  is  "SLX1A"  I think as there is no human gene called "SLX1" and is the same as GIYD1
## GIYD1  is also "SLX1A" I think!
## UAF1 and WDR48 are the same gene - you use both WDR48
## WDR48 and ATP5E and CTIP you won't see as 
    
 wanted.genes<-unique(c("DDX41","TET2", "GATA2", "ASXL1", "NOTCH1", "IDH1", "JAK2","WT1","MLL","KRAS","FLT3","IDH2","IDH1","TP53","KIT","NPM1","JAK2","DNMT3A","TET2","RUNX1","NRAS","CEBPA","PTPN11","U2AF1","SMC1A","SMC3","PHF6","STAG2","RAD21","FAM5C","EZH2","HNRNPK","FANCA","FANCB","FANCC","FANCD1","FANCD2","FANCE","FANCF","FANCG","FANCI","BRIP1","FANCL","FANCM","PALB2","RAD51C","SLX4","ERCC4","APITD1","STRA13","C1orf86","C19orf40","C17orf70","SLX1","MUS81","ERCC1","FAN1","EME1","EME2","MRE11A","NBN1","RAD50","FAND1","BRCA1","BARD1","RAD51","RAD51B","RAD51D","XRCC2","XRCC3","RMI1","RMI2","BLM","TOP3A","RPA1","RPA2","RPA3","ATM","ATR","ATRIP","CHECK1","RAD9A","RAD17","CS","DLAT","DLD","DLST","FH","IDH1","IDH2","IDH3A","IDH3B","IDH3G","MDH1","MDH2","ACLY","ACO1","OGDH","ACO2","PC","PCK1","PCK2","PDHA1","PDHA2","PDHB","OGDHL","SDHA","SDHB","SDHC","SDHD","SUCLG2","SUCLG1","SUCLA2","EGFR","ERBB2","ERBB3","ERBB4","MET","CHEK2","FANCD1","MLH1","MLH3","NBN","RBBP8","BRAC2","SPNS1","PKN1","ARTN","FAH","FANCA","ANAPC2","NOS3","HES1","C21orf33","NDUFA1","NDUFA2","NDUFA3","NDUFA4","NDUFA4L2","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFA10","NDUFA11","NDUFA12","NDUFA13","NDUFAB1","NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10","NDUFB11","NDUFC1","NDUFC2","NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS5","NDUFS6","NDUFS7","NDUFS8","NDUFV1","NDUFV2","NDUFV3","SDHA","SDHB","SDHC","SDHD","CYC1","UQCR10","UQCR11","UQCRB","UQCRC1","UQCRC2","UQCRFS1","UQCRH","UQCRQ","COX4I1","COX4I2","COX5A","COX5B","COX6A1","COX6A2","COX6B1","COX6B2","COX6C","COX7A1","COX7A2","COX7B","COX7B2","COX7C","COX8A","COX8C","ATP5A1","ATP5B","ATP5C1","ATP5D","ATP5E","ATP5F1","ATP5G1","ATP5G2","ATP5G3","ATP5H","ATP5I","ATP5J","ATP5J2","ATP5L","ATP5L2","ATP5O","ATPIF1","IDH1","IDH2","TET2","FLT3","DNMT3A","NRAS","KRAS","PTPN11","BRAF","EGFR","ERBB2","HER2","ERBB3","ERBB4","MET","ATM","ATR","ATRIP","BARD1","BLM","BRCA1","CHEK1","CHEK2","EME1","EME2","ERCC1","C17ORF70","C1orf86","C19orf40","FAN1","FANCB","FANCC","BRCA2","FANCD2","FANCE","FANCF","FANCG","FANCI","BRIP1","FANCL","FANCM","PALB2","RAD51C","SLX4","ERCC4","APITD1","STRA13","MLH1","MLH3","MRE11A","MUS81","NBN","RAD9A","RAD17","RAD50","RAD51","RAD51B","RAD51D","CTIP","RMI1","RMI2","RPA1","RPA2","RPA3","SLX1A","TOP3A","UAF1","UBA52","UBC","UBE2T","USP1","XRCC2","XRCC3","WDR48"))

    ############################ GQ.thresh<-20
    if(accumulate){
   #   best.quality<-as.numeric(a.indel[,"GQ_MEAN"])>=GQ.thresh & a.indel[,"FILTER"]=="PASS"
     # best.quality<- a.indel[, "FILTER"]=="PASS" #  &  a.indel[,"wanted.muts"]=="TRUE" #  & a.indel[,"MAF.lt:0.005"]=="TRUE"
        best.quality<-a.indel[,"MAF.lt:0.025"]=="TRUE" | a.indel[,"Gene.Names"] %in% wanted.genes
        best.quality2<-a.indel[,"Gene.Names"] %in%  wanted.genes
        sum(best.quality2)

        if(exists("wanted.genes")){
            poly.gene.site<-grep("::",a.indel[,"Gene.Names"])

            length(poly.gene.site)
            ipoly<-1

            if(length(poly.gene.site)>0){
                also.get<-{}
                for( ipoly in 1 :length(poly.gene.site)){
                    some.genes<-unlist(strsplit(a.indel[poly.gene.site[ipoly],"Gene.Names"],split="::") )
                    found<-sum(some.genes %in% wanted.genes)>0
                    if(found){
                        also.get<-c(also.get,poly.gene.site[ipoly])
                    } # found
                } # has a poly site loop
                best.quality2[also.get]<-TRUE
                best.quality[also.get]<-TRUE 
            } # has a poly

        } ## using a gene list
        


      
      if(ichr==1){
        max.maf<-max(as.numeric(gsub("^MAF.lt:","",colnames(a.indel)[grepl("^MAF.lt:",colnames(a.indel))])))
        
        max.maf.col<-paste("MAF.lt:",max.maf,sep="")
      }

 #     an.indel<-grepl("^indel",a.indel[,"TYPE"])
 #sum(an.indel)
 ##      dim(a.indel[an.indel,])
    
       if(ichr==1){
 #      write.table(a.indel,file=paste(project.name,"HC.2.BEST.chrALL.ACC",fam[ifam],project.extension,sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
#       write.table(a.indel[best.quality,],file=paste(project.name,"HC.2.BEST.chrALL.ACC_0.025",fam[ifam],project.extension,sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
       write.table(a.indel[best.quality2,],file=paste(project.name,"HC.ALL.BEST.chrALL.ACC_SUBSET2",fam[ifam],project.extension,sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
       ## write.table(a.indel[an.indel,],file=paste(project.name,".chrALL.Indel",fam[ifam],project.extension,sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
       ## write.table(a.indel[best.quality & an.indel,],file=paste(project.name,".chrALL.Indel_good_qual",fam[ifam],project.extension,sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

       
      }else{
 #      write.table(a.indel,file=paste(project.name,"HC.2.BEST.chrALL.ACC",fam[ifam],project.extension,sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
 #      write.table(a.indel[best.quality,],file=paste(project.name,"HC.2.BEST.chrALL.ACC_0.025",fam[ifam],project.extension,sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
       write.table(a.indel[best.quality2,],file=paste(project.name,"HC.ALL.BEST.chrALL.ACC_SUBSET2",fam[ifam],project.extension,sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
# write.table(summary.reduced,file=paste(file.out.name,"ALL_GENOTYPES_analysis-maf-filtered.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE) ## write in geno dump
       ## write.table(a.indel[an.indel,],file=paste(project.name,".chrALL.Indel",fam[ifam],project.extension,sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
       ## write.table(a.indel[best.quality & an.indel,],file=paste(project.name,".chrALL.Indel_good_qual",fam[ifam],project.extension,sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=TRUE)

       
             }

    } # accumulate
    ## setwd(code.dir)
    ## source("ucsc.table.names.r")
    
  

     rm(a.indel)
    
    if(dont.build.summary){indels<-{}}
    ## setwd(code.dir)
    ## source("ucsc.table.names.r")



    
    } # ichr

  
} # ifam
  

  ################## write per family data

#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/2015-08-14_AML_mixedAlignersHC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
###############################################   END   ######################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################


# 
# 
# rm(indels)
# sum(indels[,"wanted.muts"])
# dim(indels)
# 
# colnames(indels)
# 
# new.counts.Q<-indels[,c(1:6,33,54:368)]
# new.counts.Q[1:5,]
# colnames(new.counts.Q)
# setwd("/media/Bioinform-D/Research/AML-with Controls")
# setwd(analysis.dir)
# 
# write.table(new.counts.Q,file=paste(project.name,fam[ifam],"new.counts.all.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# 
# #write.table(geno.all,file=paste(project.name,fam[ifam],"geno.all.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# 
# #geno.all<-read.delim(paste(project.name,fam[ifam],"geno.all.txt",sep="."),header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
# 
# ## write.table(indels,file="aogc_specificity_analyis.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# 
# 
# ### list of snps use this on command line :: 74th is the column with the rs ids...
# ## head -50 chr10.AOGC-NGS.2013.pointwise.txt | cut -f 74
# ## /media/ga-apps/UQCCG/Programming/scripts/PerlScripts/GrepMafFiles.pl 35SNPs_or_proxies.csv 74 chr10.AOGC-NGS.2013.pointwise.txt
# ## the.samples<-sample.sheet.full[,"ParticipantCode"]
# ## the.samples
# ## indels<-indels[,paste(the.samples,"GT",sep=".")]
# # Asian-AML   Control EXOME-AML 
# 
# geno.all<-{}
# 
# target<-"AML"
# the.samples<-sample.sheet.full[sample.sheet.full[,"SampleProject"]==target,"ParticipantCode"]
# the.samples
# genotypes<-indels[,paste(the.samples,"GT",sep=".")]
# summary.geno<-genotype.summary(as.matrix(genotypes))
# colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")
# summary.geno[1:5,]
# 
# geno.all<-summary.geno
# 
# target<-"Control"
# the.samples<-sample.sheet.full[sample.sheet.full[,"SampleProject"]==target,"ParticipantCode"]
# the.samples
# genotypes<-indels[,paste(the.samples,"GT",sep=".")]
# summary.geno<-genotype.summary(as.matrix(genotypes))
# colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")
# summary.geno[1:5,]
# 
# geno.all<-cbind(geno.all,summary.geno)
# 
# target<-"Asian-AML"
# the.samples<-sample.sheet.full[sample.sheet.full[,"SampleProject"]==target,"ParticipantCode"]
# the.samples
# genotypes<-indels[,paste(the.samples,"GT",sep=".")]
# summary.geno<-genotype.summary(as.matrix(genotypes))
# colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")
# summary.geno[1:5,]
# 
# geno.all<-cbind(geno.all,summary.geno)
# 
# geno.all[1:5,]
# 
# #insert about line after col 48:
# setwd("/media/Bioinform-D/Research/AML-with Controls")
# setwd(analysis.dir)
# write.table(cbind(indels[,1:6],geno.all),file=paste(project.name,fam[ifam],"geno.all.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
#   
# ## remove.cols<-c("ALS.Genes","mouse.defect","sewell.cycling","Dequeant.cycling","ingenuity.bone.genes","Hypomagnesaemia","ProtonPump"
# ##                colnames(indels)[grepl(".AD$",colnames(indels))]
# ##                colnames(indels)[grepl(".Q$",colnames(indels))]
# ##                )
# 
# 
# 
# 
# ann.order<-c("refGene","knownGene","ensGene") # indels[1:5,paste(ann.order,"gene",sep="::")]
# gene.names<-get.gene.names(indels,paste(ann.order,"gene",sep="::"))
# anno.DB.location.core<-"/media/Bioinform-D/Research/annovar/humandb"
# 
# 
# 
# require("biomaRt")
# 
# if(genome.build=="hg18"){mart<-useMart("ensembl_mart_51",host="may2009.archive.ensembl.org",dataset="hsapiens_gene_ensembl",archive=TRUE)
#     }else if(genome.build=="hg19"){
#       mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
#     }else if(genome.build=="mm9"){
#       mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl");print("here:")
#  }
# 
# genes<-unique(unlist(gene.names[["ensGene::gene"]]))
# genes<-genes[!grepl("NONE",genes)]  #### NONE(dist=NONE exists)
# 
# 
# if(genome.build=="hg19"){
# gene.ann.wanted<-c("ensembl_gene_id","hgnc_symbol","gene_biotype","description")
# }
# if(genome.build=="hg18"){
# gene.ann.wanted<-c("ensembl_gene_id","hgnc_symbol","biotype","description")
# }
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
# ######### get rid of dupicate entries with resplce to the fist column  dim(gene.desc)[1]==length(unique(gene.desc[,1]))
#   dups<-duplicated(gene.desc[,1])
#   gene.desc<-gene.desc[!dups,]
# ####################### match genes to another column #########################################
# 
# gene.desc.table<-match.cols.and.collapse(gene.names,"ensGene::gene",gene.desc,"ensembl_gene_id",gene.ann.wanted,"delimit")
# 
# gene.desc.table[1:5,]
# 
# 
# indels[1:50,"refGene::gene"]
# gene.desc.table[1:50,"hgnc_symbol"]
# sum(is.na(indels[,"refGene::gene"]))
# 
# ################# add to refGene to get better gene list matches 
# test.is.missing<-grepl("NONE",indels[,"refGene::gene"]) | grepl("NO_ANNOVAR_",indels[,"refGene::gene"]) | is.na(indels[,"refGene::gene"]) | indels[,"refGene::gene"]=="NA"
# sum(test.is.missing)
# 
# ## gene.names.hgnc<-get.gene.names(gene.desc.table,"hgnc_symbol")
# 
# ## sum(test.is.missing)
# ## if(sum(test.is.missing)){
# ##   locations<-
# ##   mapply(clean.unique.combine, gene.names$"ensGene::gene"[grep(TRUE,test.is.missing)], gene.names.hgnc, SIMPLIFY=FALSE)
# ## indels[test.is.missing,"refGene::gene"]<-gene.desc.table[test.is.missing,"hgnc_symbol"]
# ## gene.names<-get.gene.names(indels,paste(ann.order,"gene",sep="::"))
# ## }
# 
# ## gene.desc.table[grep("::",gene.desc.table[,"hgnc_symbol"]),"hgnc_symbol"]
# ## indels[grep("::",gene.desc.table[,"hgnc_symbol"]),"refGene::gene"]
# ## gene.names[[1]][1:20]
# ## gene.names.hgnc[[1]][1:20]
# 
# ## clean.unique.combine<-function(x,y){
# ##   z<-unique(c(x,y))
# ##   z[!is.na(z) & !grepl("^PLACE_",z)] #z[z!="-" & !is.na(z) & !grepl("^PLACE_",z)]
# ## }
# 
# 
# 
# if(project.name=="2013-02-27_AML_with_AOGCControl"){
# 
# 
# a.bone.file<-read.delim(paste(anno.DB.location.core,"/", "AML_publications.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
# a.bone.file[1:5,]
# a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
# AML_publications<-as.data.frame(a.bone.file[,2]) #
# colnames(AML_publications)<-colnames(a.bone.file)[2] # unique(AML_ALL_Mutations)
# ## colnames(AML_publications)<-"AML_publications"
# 
# 
# a.bone.file<-read.delim(paste(anno.DB.location.core,"/", "AML_MDS_publications.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
# a.bone.file[1:5,]
# a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
# AML_MDS_publications<-as.data.frame(a.bone.file[,2]) #
# colnames(AML_MDS_publications)<-colnames(a.bone.file)[2] # unique(AML_ALL_Mutations)
# ## colnames(AML_MDS_publications)<-"AML_MDS_publications" # unique(AML_ALL_Mutations)
# 
# a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_Validated.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
# a.bone.file[1:5,]
# a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
# AML_Validated<-as.data.frame(a.bone.file[,1]) #
# colnames(AML_Validated)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)
# 
# a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_SNP_ALL_Confident.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
# a.bone.file[1:5,]
# a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
# AML_SNP_ALL_Confident<-as.data.frame(a.bone.file[,1]) #
# colnames(AML_SNP_ALL_Confident)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)
# 
# a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_SNP_Diagnosis.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
# a.bone.file[1:5,]
# a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
# AML_SNP_Diagnosis<-as.data.frame(a.bone.file[,1]) #
# colnames(AML_SNP_Diagnosis)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)
# 
# a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_SNP_Relapse.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
# a.bone.file[1:5,]
# a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
# AML_SNP_Relapse<-as.data.frame(a.bone.file[,1]) #
# colnames(AML_SNP_Relapse)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)
# 
# a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_DINDEL_Diagnosis.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
# a.bone.file[1:5,]
# a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
# AML_DINDEL_Diagnosis<-as.data.frame(a.bone.file[,1]) #
# colnames(AML_DINDEL_Diagnosis)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)
# 
# a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_DINDEL_Relapse.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
# a.bone.file[1:5,]
# a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
# AML_DINDEL_Relapse<-as.data.frame(a.bone.file[,1]) #
# colnames(AML_DINDEL_Relapse)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)
# 
#  
#   
# a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_ALL_Mutations.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
# a.bone.file[1:5,]
# a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
# AML_ALL_Mutations<-as.data.frame(a.bone.file[,1]) #
# colnames(AML_ALL_Mutations)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)
# 
# a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_ASH_Mutations.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
# a.bone.file[1:5,]
# a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
# AML_ASH_Mutations<-as.data.frame(a.bone.file[,1])
# colnames(AML_ASH_Mutations)<-colnames(a.bone.file)[1] # unique(AML_ASH_Mutations)
# 
# a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_CLINICAL_Mutations.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
# a.bone.file[1:5,]
# a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
# AML_CLINICAL_Mutations<-as.data.frame(a.bone.file[,1])
# colnames(AML_CLINICAL_Mutations)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)
# 
# 
# a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_SNP_Mutations.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
# a.bone.file[1:5,]
# a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
# AML_SNP_Mutations<-as.data.frame(a.bone.file[,1])
# colnames(AML_SNP_Mutations)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)
# 
# a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_DINDEL_Relapse.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
# a.bone.file[1:5,]
# a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
# AML_DINDEL_Relapse<-as.data.frame(a.bone.file[,1])
# colnames(AML_DINDEL_Relapse)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)
# 
# a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_ENU_Mutations.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
# a.bone.file[1:5,]
# a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
# AML_ENU_Mutations<-as.data.frame(a.bone.file[,1])
# colnames(AML_ENU_Mutations)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)
# 
# a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_FRANC_Mutations.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
# a.bone.file[1:5,]
# a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
# AML_FRANC_Mutations<-as.data.frame(a.bone.file[,1])
# colnames(AML_FRANC_Mutations)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)
# 
# a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_LEY_RELAPSE_Mutations.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
# a.bone.file[1:5,]
# a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
# AML_LEY_RELAPSE_Mutations<-as.data.frame(a.bone.file[,1])
# colnames(AML_LEY_RELAPSE_Mutations)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)
# 
# a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_MITO_Mutations.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
# a.bone.file[1:5,]
# a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
# AML_MITO_Mutations<-as.data.frame(a.bone.file[,1]) # unique(as.data.frame(a.bone.file[,1]))
# colnames(AML_MITO_Mutations)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)
# 
# a.bone.file<-read.delim(paste(anno.DB.location.core,"/","AML_OTHER_Mutations.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
# a.bone.file[1:5,]
# a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
# AML_OTHER_Mutations<-as.data.frame(a.bone.file[,1])
# colnames(AML_OTHER_Mutations)<-colnames(a.bone.file)[1] # unique(AML_ALL_Mutations)
# 
# 
# } ## project AML
# 
# 
# 
# 
# 
# 
# 
# ######################## Here OMIM table shoudl always be present to gene.table is never empty
# ## target<-gsub(".TGCM-AML.analysis-maf-filtered.txt","",project.files[ichr])
# ## get<-try(load(paste(gsub(".txt$","",target),".GeneLists.RData",sep="")),silent=TRUE)   
# 
# gene.lists<-c("AML_Validated","AML_SNP_ALL_Confident","AML_SNP_Diagnosis","AML_SNP_Relapse","AML_DINDEL_Diagnosis","AML_DINDEL_Relapse","AML_CLINICAL_Mutations","AML_FRANC_Mutations","AML_ENU_Mutations","AML_publications","AML_MDS_publications","AML_SNP_Mutations","AML_ALL_Mutations", "AML_OTHER_Mutations","AML_MITO_Mutations","AML_LEY_RELAPSE_Mutations","AML_ASH_Mutations")
# current.objects<-ls()
# gene.lists<-gene.lists[gene.lists %in% current.objects]
# gene.table<-{}
# for(it in 1:length(gene.lists)){
#  # print(it)
#   if(!exists(gene.lists[it])){print("missing")}
#   a.temp<-eval(as.name(gene.lists[it]))
#    print(paste(gene.lists[it],"-->",colnames(a.temp),sep=""))
#   if(is.null(dim(gene.table)) ){gene.table<-a.temp}else{ if( dim(a.temp)[1]==dim(indels)[1] ){gene.table<-cbind(gene.table,a.temp)}else{print("dimension mismatch")} }
#                                   }
# rm("a.temp")
# if(dim(gene.table)[1]!=dim(indels)[1]){print("ERROR gene.table DIFFERENT number of rows to indels");gene.table<-""}
# dim(gene.table)
# colnames(gene.table)
# gene.table[1:5,]
# #####################################################
# 
# 
# 
# 
# 
# ####################################### RUN ANNOVAR ONCE ############################################
# 
# indels[,"chr"]<-gsub("chr","",indels[,"chr"])
# key.indels<-build.key(indels,core.ann)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# genome.build<-"hg19"
# maf.threshold<-0.0 
# 
# anno.DB.location.core<-"/media/Bioinform-D/Research/annovar/humandb"
# anno.DB.location<-paste(anno.DB.location.core,genome.build,sep="/")
# 
# ######################### Info the run
# generic.filter.DB<-c("hg19_sequnome.txt")    # returns 2  extra columns (DB,score)
# names(generic.filter.DB)<-c("AML-sequnome")
# 
# update.annovar.annotations<-TRUE
# 
# 
# setwd(code.dir)
# source("ucsc.table.names.r")   # load in the UCSC tables these use the db file names and not their lable-names 
# source("annotate_SNPs_subroutines.r")
# source("ucsc.table.names.processor.r")
# #####################
# 
# setwd(annotate.dir)
# files.in.annotate<-dir(annotate.dir)
# 
# target<-"extra_match.txt"
# write.table(indels[,core.ann],file=target,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE) ## used for annovar
# 
# 
# for(i in 1:length(generic.filter.DB)){
# if( (paste(target,generic.filter.DB[i],"log",sep=".") %in% files.in.annotate ) & !update.annovar.annotations){next} # perhaps have added batabases but dones want to update
# system(
#     paste("annotate_variation.pl  -filter  -otherinfo  --score_threshold ",maf.threshold," -buildver ",genome.build," -dbtype  generic -genericdbfile ",generic.filter.DB[i]," ",target,"  ",anno.DB.location," --outfile ",paste(target,generic.filter.DB[i],sep=".")," ",sep="" )
#       )
# }
# 
# 
# #annotate_variation.pl  -filter get.cols 4,5 --score_threshold 0 -buildver hg19 -dbtype  generic -genericdbfile hg19_sequnome.txt extra_match.txt  /media/Bioinform-D/Research/annovar/humandb/hg19 --outfile extra_match.txt.hg19_sequnome.txt2 
# ##############################################################################################################################################
# ##############################################################################################################################################
# ##############################################################################################################################################
# ##############################################################################################################################################
# ##############################################################################################################################################
# ##############################################################################################################################################
# ##############################################################################################################################################
# ##############################################################################################################################################
# 
# ######################## filter Dbs these have a _dropped ####
# 
# if(!grepl("^chr",indels[1,"chr"])){  ### just in case restarted from a strange place
# key.indels<-build.key(indels,core.ann,add.chr.label=TRUE)
# }else{key.indels<-build.key(indels,core.ann)      }
# 
# 
# ann.table<-data.frame(key=key.indels,stringsAsFactors=FALSE)
# ann.order<-{}
# 
# the.combined.DBs<-c(generic.filter.DB)  # c(filter.DB,generic.filter.DB,function.filter.DB)
# the.DBs<-the.combined.DBs
# 
# the.files<-dir(annotate.dir)
# the.files<-the.files[!grepl("\\.log$",the.files) & grepl(paste("^",target,sep=""),the.files)]
# the.files
# 
# for(i in 1:length(the.DBs)){
# the.DB<-the.DBs[i]
# 
# if(grepl("+",the.DB,fixed=TRUE)){the.DB<-gsub("+","\\+",the.DB,fixed=TRUE)}
# 
# target.string<-paste("^",target,".",the.DB,"\\S*_dropped$",sep="")
# #if(grepl("++",target.string,fixed=TRUE)){target.string<-gsub("++","\\++",target.string,fixed=TRUE)} ##gerp++ casles problems in grep (++)
# 
# data.file<-the.files[grepl(target.string,the.files)]
# 
# if(grepl("\\+",the.DB,fixed=TRUE)){the.DB<-gsub("\\+","",the.DB,fixed=TRUE)}
# print(paste(the.DB,data.file,sep=" -> "))
# 
# the.labels<-eval(as.name(paste(the.DB,"labels",sep=".")))
# the.labels.wanted<-eval(as.name(paste(the.DB,"labels.wanted",sep=".")))
#   
# the.data<-try(read.delim(data.file,header=F,skip=0,sep="\t",fill=TRUE,stringsAsFactors=FALSE,colClasses = "character"),silent=TRUE) 
# if(!inherits(the.data, "try-error")){
# 
# colnames(the.data)<-the.labels
# #the.key<-paste(the.data[,"chr"],the.data[,"start"],the.data[,"end"],the.data[,"REF"],the.data[,"ALT"],the.data[,"TYPE"],sep=":")
# 
# the.key<-build.key(the.data,core.ann,add.chr.label=TRUE) #sep is always ":"
# if(length(the.key)!=length(unique(the.key))){
#   print("WARNING duplicate entries found") ## Hapmap / 1000 Genomes an even within dbSNP can have duplicates key above shoudl work but if not
#   posns<-match(unique(the.key),the.key)
# 
#   missing<-is.na(posns)
#   the.data<-the.data[posns[!missing],]
#   the.key<-the.key[posns[!missing]]
# }
# 
# rownames(the.data)<-the.key
# 
# the.data[1:5,]
# length(the.key)>0
# 
# posns<-match(key.indels,the.key)
# missing<-is.na(posns) # 
# sum(missing)
# 
# }else{print("WARNING could not find DATABASE OR NO MATCHES");missing<-rep(TRUE,times=length(key.indels))}## no data
# 
#  the.ann.label<-paste(the.DB,toString(the.labels.wanted),sep="::")
# 
# if(!inherits(the.data, "try-error")){
# ann.table[!missing,the.ann.label]<-the.data[posns[!missing],the.labels.wanted]
# }else{ann.table[!missing,the.ann.label]<-as.character(NA)}
# ann.table[,the.DB]<-!missing
# #ann.table
# ann.order<-c(ann.order,the.DB)
# }
# #############################3
# rownames(ann.table)<-ann.table[,"key"]
# ann.table<-ann.table[,-1]
# 
# ann.table[1:5,]
# sum(!is.na(ann.table[,1]))
# sum(ann.table[,2])
#  ann.table[!is.na(ann.table[,1]),] # check have WT1
# 
# 
# 
# setwd(analysis.dir)
# 
# write.table(ann.table,file=paste(project.name,fam[ifam],"ann.table.all.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# 
# 
# 
# #########################################################
# #Now create final merged outout file
# colnames(indels)
# 
# indels<-indels[,colnames(indels)[!(colnames(indels) %in% remove.cols)]]
# tapply(indels[,"chr"],indels[,"chr"],length)
# 
# colnames(indels)
# 
# ########## specyail for aml
# sum(geno.all[,"ALT.Alleles.Asian-AML"] != indels[,"ALT.Alleles.Asian-AML"])
# colnames(new.counts.Q)
# colnames(ann.table)
# colnames(geno.all)
# indels[,"ok.in.group"]<-new.counts.Q[,"ok.in.group"]
# insert.location<-53
# colnames(indels)
# test<-cbind(indels[1:5,1:81],ann.table[1:5,],indels[1:5,84:288],new.counts.Q[1:5],indels[1:5,595:dim(indels)[2]])
# 
# 
# indels<-cbind(indels[,1:40],geno.all[,7:30],indels[,65:81],ann.table,indels[,84:288],new.counts.Q,indels[,595:dim(indels)[2]])
# ######
# 
# insert.location<-40
# indels<-cbind(indels[,1:insert.location],geno.all,gene.table,ann.table,indels[,(insert.location+1):dim(indels)[2]])
# 
# 
# sum(as.logical(indels[,"wanted.muts"]))
# sum(as.logical(indels[,"wanted.muts"]) | as.logical(ann.table[,2]))
# indels[,"wanted.muts"]<-as.logical(indels[,"wanted.muts"]) | as.logical(ann.table[,2])
# indels[,"wanted.muts.coding"]<-as.logical(indels[,"wanted.muts.coding"]) | as.logical(ann.table[,2])
# 
# maf.cols<-colnames(indels)[grep("^MAF.lt",colnames(indels))]
# #imaf<-1 
# for(imaf in 1:length(maf.cols)){
# indels[,maf.cols[imaf]]<-as.logical(indels[,maf.cols[imaf]]) | as.logical(ann.table[,2])
# }
# #sum(indels[1:5,maf.cols[length(maf.cols)]]!=(indels[1:5,maf.cols[length(maf.cols)]] | indels[1:5,"hg19_sequnome.txt"]))
# 
# ref.genotype.labels<-colnames(indels)[grepl(".GT",colnames(indels))]
# ref.genotype.labels
# ref.counts<-apply(indels[,ref.genotype.labels],1,function(x){ sum( x=="0/0"| is.na(x))})
# dim(indels)
# length(ref.counts)
# ref.counts<-ref.counts==length(ref.genotype.labels)
# sum(!ref.counts)
# 
# #indels[!ref.counts,][1,ref.genotype.labels]
#                                         # ref.only[a.test]
# #indels<-indels[!ref.counts,]
# 
# 
# getwd()
# setwd(analysis.dir)
# analysis.dir
# dim(indels)
# write.table(indels,file=paste(project.name,fam[ifam],"ALL_with_ref.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# 
# indels<-indels[!ref.counts,]
# 
# write.table(indels,file=paste(project.name,fam[ifam],"ALL.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# 
# indels<-indels[as.logical(indels[,maf.cols[length(maf.cols)]]),]
# 
# dim(indels)
# write.table(indels,file=paste(project.name,fam[ifam],"ALL-maf-filtered.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# 
# 
# indels<-indels[as.logical(indels[,"wanted.muts"]),]  ### wanted and maf< max.maf
# dim(indels)
# write.table(indels,file=paste(project.name,fam[ifam],"ALL-wanted.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# 
# 
# 
# write.table(indels[,colnames(geno.all)],file=paste(project.name,fam[ifam],"All-wanted_geno.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# 
# 
# ###############################################
# sum(!indels[,"wanted.muts"])
# indels[!indels[,"wanted.muts"],1:5]
# 
# setwd("/media/ga-apps/UQCCG/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/Analysis")
# 
# setwd(annotate.dir)
# 
# load("update.RData")
# setwd(analysis.dir)
# indels<-read.delim("2013-02-27_AML_with_AOGCControl.TGCM-AML.All.txt",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
# 
# indels<-read.delim("AOGC-Genotyping.output.AOGCHBMFam.wanted.All-maf-filtered.all.txt",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
# # "AOGC-Genotyping.output.AOGCHBMFam.ALL-wanted.txt"
# #use<-read.delim("2013-02-27_AML_with_AOGCControl.TGCM-AML.All-wanted.txt",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
# 
# geno.all<-read.delim(paste(project.name,fam[ifam],"geno.all.txt",sep="."),header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
# ## write.table(new.counts.Q,file=paste(project.name,fam[ifam],"new.counts.all.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# 
# new.counts.Q<-read.delim(paste(project.name,fam[ifam],"new.counts.all.txt",sep="."),header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
# ann.table<-read.delim(paste(project.name,fam[ifam],"ann.table.all.txt",sep="."),header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
# 
# dim(new.counts.Q)
# dim(ann.table)
# 
# dim(gene.table) 
# dim(geno.all)
# dim(indels)
# 
# colnames(geno.all)
# geno.all[1000,1:25]
# indels[1000000,1:105]
# sum(geno.all[,"start"]!=indels[,"start"])
# indels[,colnames(geno.all)]<-geno.all
# 
# core.ann
# use.key<-build.key(geno.all,core.ann)
# indels.key<-build.key(indels,core.ann)
# 
# posns<-match(indels.key,use.key)
# missing<-is.na(posns)
# sum(missing) ## should be zero
# 
# ## indels[missing,1:7]
# ## grep("9047242",geno.all[,"start"])
# 
# geno.all.match<-geno.all[posns[!missing],]
# dim(geno.all.match)
# dim(indels)
# dim(gene.table)
# 
# colnames(gene.table)
# colnames(indels)
# colnames(geno.all.match)
# 
# geno.all.match[1000,]
# indels[1000,1:72]
# insert.location<-37
# indels<-cbind(indels[,1:insert.location],geno.all.match[,7:dim(geno.all.match)[2]],gene.table,indels[,(insert.location+1):dim(indels)[2]])
# 
# write.table(indels,file="AOGC-Genotyping.output.AOGCHBMFam.wanted.All-maf-filtered.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# 
# 
# #write.table(indels,file="2013-02-27_AML_with_AOGCControl.TGCM-AML.All-wanted.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# ##############################################################################################################################################
# ##############################################################################################################################################
# ##############################################################################################################################################
# ##############################################################################################################################################
# ##############################################################################################################################################
# ##############################################################################################################################################
# ##############################################################################################################################################
# ##############################################################################################################################################
# ##############################################################################################################################################
# 
# 
# 
# ################ collapse a large list of publications into one list
# 
# 
# setwd("/media/ga-apps/UQCCG/Projects/TGCM-AML/Analysis/AML_Prelim_analysis")
# file<-"Publication_gene_lists.csv"
# file<-"Publication_MDS_gene_lists.csv"
# 
# 
# setwd("/media/Bioinform-D/Research/annovar/humandb")
# file<-"bone_density.raw.csv"
# 
# data<-read.delim(file,header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
# dim(data)
# data[1:50,]
# 
# data.all<-{}
# genes.all<-{}
# for(i in 1:dim(data)[2]){
#   missing<- data[,i]=="" | is.na(data[,i])
#   
#   genes<-data[,i]
#   genes[missing]<-"missing"
# 
#   genes.all<-c(genes.all,genes)
#   data.all<-c(data.all,rep(colnames(data)[i],times=length(data[,i]) ))
#   
#  
# }
# 
# 
# 
# junk<-grepl("^missing",genes.all)
# data.all<-data.all[!junk]
# genes.all<-genes.all[!junk]
# 
# 
# #data.all<-cbind(genes,data.all)
# 
# unique.genes<-unique(genes.all)
# anno<-{}
# for(i in 1:length(unique.genes)){
#   posns<-grep(unique.genes[i],genes.all)
#   a.list<-paste(data.all[posns],collapse=";")
#   a.list<-paste(length(posns),"-",a.list,sep="")
#   anno<-c(anno,a.list)
# }
# 
# data<-cbind(unique.genes,anno)
# colnames(data)<-c("Gene","MDS_Papers")
# 
# write.table(data,file="AML_MDS_publications.csv",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# write.table(data,file="AML_publications.csv",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# 
# 
# 
# output<-cbind(unique.genes,unique.genes)
# 
# 
# colnames(output)<-c("BoneDensity_Sequencing","BoneDensity_Sequencing")
# colnames(output)<-c("BoneDensity_GWAS","BoneDensity_GWAS")
# output[1:5,]
# 
# write.table(output,file="BoneDensity_Sequencing.csv",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# write.table(output,file="BoneDensity_GWAS.csv",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# ###############################################
# 
# 
# 
# 
# setwd("/media/scratch2/AOGC-NGS/Analysis")
# file<-"AOGC-NGS.all.CHR.analysis-maf-filtered.txt" # "AOGC-NGS.wanted.CHR.analysis-maf-filtered.txt"
# 
# aogc<-read.delim(file,header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
# 
# ## colnames(aogc)
# ## tapply(aogc[,"chr"],aogc[,"chr"],length)
# ## get<-aogc[,"chr"] %in% c("1","2","3","4","5","6","7","8","9")
# ## sum(get)
# ## sum(!get)
# ## write.table(aogc[get,],file="AOGC-NGS.chr1to9.analysis.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# ## write.table(aogc[!get,],file="AOGC-NGS.chr9toY.analysis.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# 
# sum(grepl(".GT",colnames(aogc)))
# aogc<-aogc[aogc[,"MAF.lt.0"],]
# counts<-tapply(aogc[,"ALT.Alleles.ALL"],aogc[,"Gene.Names"],sum)
# 
# counts.all<-tapply(aogc[,"TOTAL.Alleles.ALL"],aogc[,"Gene.Names"],sum)
# 
# counts.all<-counts.all[names(counts)]
# counts<-counts/counts.all
# counts<-counts*96
# counts
# 
# counts<-cbind(names(counts),counts)
# colnames(counts)<-c("Gene","count")
# setwd("/media/Bioinform-D/Research/annovar/humandb")
# write.table(counts,file="aogc.to.96.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# 
# 
# ### list of snps use this on command line :: 74th is the column with the rs ids...
# ## head -50 chr10.AOGC-NGS.2013.pointwise.txt | cut -f 74
# ## /media/ga-apps/UQCCG/Programming/scripts/PerlScripts/GrepMafFiles.pl 35SNPs_or_proxies.csv 74 chr10.AOGC-NGS.2013.pointwise.txt  
# 
# remove.cols<-c("ALS.Genes","mouse.defect","sewell.cycling","Dequeant.cycling","ingenuity.bone.genes","Hypomagnesaemia","ProtonPump"
#                colnames(indels)[grepl(".AD$",colnames(indels))]
#                colnames(indels)[grepl(".Q$",colnames(indels))]
#                )
# 
# 
# 
# 
# ann.order<-c("refGene","knownGene","ensGene") # indels[1:5,paste(ann.order,"gene",sep="::")]
# gene.names<-get.gene.names(indels,paste(ann.order,"gene",sep="::"))
# 
# 
#        aogc[1:10,c("Gene.Names", "ALT.Alleles.ALL","MAF.lt.0")]"TOTAL.Alleles.ALL"
