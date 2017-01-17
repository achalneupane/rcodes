
ALDH1B1:c.1132C>T (p.Gln378Ter) chr9

ALDH2:c.1510G>A (p.Glu504Lys) chr12 






coverage.file<-"/media/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_SAMPLE_Mon_Jul_20_2015.txt"


tcga<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/New_sample_set/Full.sample_sheet.wP.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

good.samples<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/New_sample_set/pheno.ori.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

aogc<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/New_sample_set/AOGC_novoaligned_with_gvcf.tsv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

coverage.file<-"/media/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_SAMPLE_Mon_Jul_20_2015.txt"
seq.type.file<-coverage.file  ## Troels-covergae prior to top up
seq.type<-read.delim(seq.type.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE,check.names=FALSE)


contaminated.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/contaminated_AOGC_SEQ_samples.txt"
contaminated<-read.table(contaminated.file,header=F,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

seq.type[1:5,]
good.samples[1:5,]
tcga[1:5,]
aogc[1:5,]



posns<-match(aogc[,"SM"],seq.type[,"Sample"])
missing<-is.na(posns)
sum(missing)# 0
seq.type<-seq.type[posns,] # just get aogc

order.by<-order(seq.type[,"percent.ccds.gt.10"],decreasing=TRUE)
seq.type<-seq.type[order.by,]


hist(seq.type[,"percent.ccds.gt.10"])
bad.coverage<-seq.type[,"percent.ccds.gt.10"]<75
sum(!bad.coverage)
sum(bad.coverage)
## low.coverage<-bad.coverage & seq.type[,"Project"]=="AOGC"
## sum(low.coverage)
low.coverage<-seq.type[bad.coverage,"Sample"]


sum(tcga[,"ParticipantCode"] %in% low.coverage) #131
sum(good.samples[,"SAMPLE"] %in% low.coverage) # 15


tcga.samples<-tcga[,"ParticipantCode"] 
low.coverage
tcga.samples

exclude.samples<-unique(c(low.coverage,tcga.samples))
)
seq.type.sub<-seq.type[!(seq.type[,"Sample"] %in% exclude.samples),]
seq.type.sub[1:5,]
dim(seq.type.sub)
tail(seq.type.sub)

pheno.ori[1:5,]

new.controls<-seq.type[,"Sample"]
dim(pheno.ori)
table(pheno.ori[,"Project"])

is.AOGC<-pheno.ori[,"Project"] =="AOGC-NGS"
sum(is.na(is.AOGC))
pheno.ori[is.na(is.AOGC),]

is.AOGC[is.na(is.AOGC)] <- TRUE
use<-pheno.ori[ !is.AOGC,]
table(use[,"Project"])
table(use[,"SampleProject"])
table(use[,"AffectionStatus"])

use.case<-use[,c("SAMPLE","Project", "PaternalID","MaternalID","Sex","AffectionStatus")]
use.control<-cbind(seq.type.sub[,c("Sample","Project")], 0,0,2,1)
colnames(use.control)<-c("SAMPLE","Project", "PaternalID","MaternalID","Sex","AffectionStatus")
use.control[1:5,]


the.samples<-rbind(use.case,use.control)



seq.bam.file<-"/media/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_BAM_Mon_Jul_20_2015.txt" ## Troels-covergae prior to top up
seq.bam<-read.delim(seq.bam.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE,check.names=FALSE)


posns<-match(the.samples[,"SAMPLE"],seq.bam[,"Sample"])
missing<-is.na(posns)
sum(!missing)# 0


keep<-seq.bam[,"Sample"] %in% the.samples[,"SAMPLE"]
sum(keep)

seq.bam[keep,c("Project","BAM","ID","Sample")][1:5,]

the.samples.bam<-seq.bam[keep,c("Project","BAM","ID","Sample")]
write(the.samples.bam

seq.type<-seq.type[posns,] # just get aogc
write.table(the.samples.bam,file="new_AML_CALL_SAMPLE_BAMS.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(the.samples,file="new_AML_CALL_SAMPLEs.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
getwd()

################# build the sample sheet
################# build the sample sheet
################# build the sample sheet
################# build the sample sheet
################# build the sample sheet
################# build the sample sheet
the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/sample_sheet.for.analysis.csv"

the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/2015-08-14_AML_mixedAligners.chr2.ALL.Sample_Sheet_NEW.txt"
sample.sheet.full<-read.delim(the.sample.sheet,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

#all.possible.samples<-unique(sample.sheet.full[,"ParticipantCode"])


#the.sample.sheet<-"/media/UQCCG/Sequencing/CompleteGenomics/Chort_descriptions/Sequencing comparisons-ver 4.csv"
#sample.sheet.new.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/2015-03-16_AllAMLandLung.ALL.Sample_Sheet_NEW.csv"
sample.sheet.new.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/2015-08-14_AML_mixedAligners.chr2.ALL.Sample_Sheet_NEW.txt"

sample.sheet.new<-read.delim(sample.sheet.new.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
colnames(sample.sheet.new)
#"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/BAM/TGCM-AML-combine_SampleSheet.csv"
qc<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/Vcf_Merge.GQ-20.Min.ALT-14_2015-03-16_AllAMLandLung.ARE_RELATED.Fri_Jul_03_2015.genetic_QC.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

the.sample.sheet<-"/media/UQCCG/Sequencing/CompleteGenomics/Chort_descriptions/Sequencing comparisons-ver 7.csv"
sample.sheet.full<-read.delim(the.sample.sheet,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
#sample.sheet.full<-read.delim(the.sample.sheet,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

colnames(sample.sheet.full)[colnames(sample.sheet.full)=="Sequence.Number"]<-"ParticipantCode"
sample.sheet.full[1:5,1:10]
sample.sheet.new[1:10,]
colnames(sample.sheet.full)
dim(sample.sheet.full)
sample.sheet.full[1:5,1:10]
qc[1:5,]

sample.sheet.full[1:5,1:10]
dim(sample.sheet.full)
found<-sample.sheet.full[,"ParticipantCode"] %in% all.possible.samples 
sample.sheet.full[!found,"ParticipantCode"] # [1] there are annotated but not sequenced:  "95"           "AMLM12035D-T" "AMLM12008CAE" "AMLM12024L-G" "AMLM12008AM"  "1406601-PQ" 

found<- all.possible.samples %in% sample.sheet.full[,"ParticipantCode"]
sum(!found)
all.possible.samples[!found] 

all.possible.samples[grep("PGB",all.possible.samples)]
#"AMLM12030PGB"    "AMLM12PAH030PGB"
all.possible.samples[grep("BJD",all.possible.samples)]

all.possible.samples[grep("1409101",all.possible.samples)]

## > all.possible.samples[grep("PGB",all.possible.samples)]
## [1] "AMLM12030PGB"    "AMLM12PAH030PGB"  ## PAH is crap coverage
## > all.possible.samples[grep("BJD",all.possible.samples)]
## [1] "AMLM12038BJD"    "AMLM12PAH038BJD" ## PAH is crap coverage  place in contaminated
## ParticipantCode:AMAS-25.3-Diagnostic;LibraryPlateLocation:NLS-20140812-EDAM/SJ/JH-AMAS-LeoPharma-D11;OldPCode:1420901;
## AMAS-18.3- is 25.3

## 13468(SDDS)	AMLM12035D-T are sibs
## AMLM12035D-T	G11F
## 13380	AMAS-18.3-DiagnosticMouseBlood
## 13474	AMLM12035D-T
## AMAS-18.3-DiagnosticMouseBlood	AMAS-25.3-Diagnostic
## AMAS-18.3-DiagnosticMouseSpleenLiver	AMAS-25.3-Diagnostic
## 13380	AMAS-18.3-DiagnosticMouseSpleenLiver
## AMAS-18.3-DiagnosticMouseBlood	AMLM12033MD
## 13381	AMAS-18.3-DiagnosticMouseSpleenLiver
## 37.2	AMAS-18.3-DiagnosticMouseBlood
## AMAS-18.3-DiagnosticMouseSpleenLiver	AMLM12004VDP
## 13379	AMAS-18.3-DiagnosticMouseSpleenLiver

### chk 25.3 is not rea;ted to 18.3 in any way
# in RSH <- AML
## RG	ID:H1178ADXX-1-05	PL:illumina	LB:H1178ADXX-1-05	DS:FCID:H1178ADXX;Lane:1;SampleID:H1178ADXX-1-05;SampleRef:Human;IndexSeq:GGACTCCT;Description:NxtD;NxtXR;Control:N;Recipe:101+9+101;Operator:SharonSong;JessicaHarris;SampleProject:RSGB;ParticipantCode:AMLM12RMH026J-N;LibraryPlateLocation:NLA-20130927-RSGB-SS-E04;	SM:AMLM12RMH026J-N
## @RG	ID:H1178ADXX-1-11	PL:illumina	LB:H1178ADXX-1-11	DS:FCID:H1178ADXX;Lane:1;SampleID:H1178ADXX-1-11;SampleRef:Human;IndexSeq:AAGAGGCA;Description:NxtD;NxtXR;Control:N;Recipe:101+9+101;Operator:SharonSong;JessicaHarris;SampleProject:RSGB;ParticipantCode:AMLM12RMH026J-N;LibraryPlateLocation:NLA-20130927-RSGB-SS-E05;	SM:AMLM12RMH026J-N
## @RG	ID:H1178ADXX-2-05	PL:illumina	LB:H1178ADXX-2-05	DS:FCID:H1178ADXX;Lane:2;SampleID:H1178ADXX-2-05;SampleRef:Human;IndexSeq:GGACTCCT;Description:NxtD;NxtXR;Control:N;Recipe:101+9+101;Operator:SharonSong;JessicaHarris;SampleProject:RSGB;ParticipantCode:AMLM12RMH026J-N;LibraryPlateLocation:NLA-20130927-RSGB-SS-E04;	SM:AMLM12RMH026J-N
## @RG	ID:H1178ADXX-2-11	PL:illumina	LB:H1178ADXX-2-11	DS:FCID:H1178ADXX;Lane:2;SampleID:H1178ADXX-2-11;SampleRef:Human;IndexSeq:AAGAGGCA;Description:NxtD;NxtXR;Control:N;Recipe:101+9+101;Operator:SharonSong;JessicaHarris;SampleProject:RSGB;ParticipantCode:AMLM12RMH026J-N;LibraryPlateLocation:NLA-20130927-RSGB-SS-E05;	SM:AMLM12RMH026J-

############# set affection statues in rought samples sheet
annotated.AML<-sample.sheet.new[,"ParticipantCode"] %in%  sample.sheet.full[,"ParticipantCode"]
sum(annotated.AML)
sample.sheet.new[annotated.AML,"AffectionStatus"]<-2
sample.sheet.new[!annotated.AML,"AffectionStatus"]<-1
#######################################################################################

posns<-match(sample.sheet.new[,"ParticipantCode"], sample.sheet.full[,"ParticipantCode"])
missing<-is.na(posns)
sum(missing)

sample.sheet.full[!missing,"ParticipantCode"]
SampleProject<-rep(NA,times=dim(sample.sheet.new)[1])

sample.sheet<-cbind(sample.sheet.new,sample.sheet.full[posns,])


table(sample.sheet.new[missing,"AffectionStatus"]) #461 all 2
a.case<-sample.sheet.new[,"AffectionStatus"]==2
a.control<-sample.sheet.new[,"AffectionStatus"]==1
a.unknown<-sample.sheet.new[,"AffectionStatus"]==9
#sample.sheet.new[a.case,]
sample.sheet.new[missing & a.case,"ParticipantCode"]
#    "AMLM12PAH030PGB" "AMLM12PAH038BJD"


sample.sheet.new[missing & a.control,"ParticipantCode"]
sample.sheet.new[missing & a.unknown,"ParticipantCode"]
 # "0413E1210023"  "0413E1211353"  "0413E14100212" "0413E14113512" "0413E2100212"  "0413E2113512" these are lung samples

asians<-c("AMLM12PAH030PGB","AMLM12PAH038BJD","NSGC-23.2","NSGC-23.3","NSGC-23.4","13380","13381","13455","13456","13457","13379","13474","37.2","37.3","63","74","83","84","9","99","AMAS-25.3-Diagnostic","AMLM12034H-F","AMLM12038BJD","AMLM12004VDP","AMLM12005R-G","AMLM12033MD","AMAS-18.3-Diagnostic","AMAS-18.3-DiagnosticMouseBlood","AMAS-18.3-DiagnosticMouseSpleenLiver","MODY_250.3","SKDP-200.3083","SKDP-200.3065","SKDP-200.3024","SKDP-200.3025","SKDP-200.3023","SKDP-200.7095","SKDP-200.3045","SKDP-200.305","SKDP-200.3084","SKDP-200.7036")
sample.sheet[1:5,1:10]

sample.sheet[a.case,"ParticipantCode"]


sample.sheet[a.case,"SampleProject"]<-"AML"
sample.sheet[a.control,"SampleProject"]<-"Control"
sample.sheet[a.unknown,"SampleProject"]<-"Control"
sample.sheet[a.case  & (sample.sheet[,"ParticipantCode"] %in% asians) ,"SampleProject"]<-"Asian-AML"
sample.sheet[a.control  & (sample.sheet[,"ParticipantCode"] %in% asians) ,"SampleProject"]<-"Asian-Control"
sum(is.na(sample.sheet[,"SampleProject"]))

colnames(sample.sheet)

table(sample.sheet[,"SampleProject"])
          ## AML     Asian-AML Asian-Control       Control 
          ## 175            17            23           363

## mixed aligeners
          ## AML     Asian-AML Asian-Control       Control 
          ## 173            15            13           448 

table(sample.sheet[, "Status" ])
               ## 1 Blood Culture         1 Spleen+Liver Culture              2m POST ALLOGRAFT              8m POST ALLOGRAFT                      Diagnosis                        Pre-AML                        Relapse Relapse (Chloroma-Skin-biopsy) 
               ##               1                              1                              1                              1                            171                              1                              8                              1 
               ##       Remission 
               ##               5 
not.diagnosis<-sample.sheet[,"Status"]!="Diagnosis" | is.na(sample.sheet[,"Status"])

sum(not.diagnosis)
sample.sheet[(not.diagnosis & a.case),c("SampleProject","ParticipantCode","Status")]

sample.sheet[(not.diagnosis & a.case),"SampleProject"]<-paste(sample.sheet[(not.diagnosis & a.case),"SampleProject"],"NotDiagnosis",sep="-")
sample.sheet[(not.diagnosis & a.case),c("SampleProject","ParticipantCode","Status")]
colnames(sample.sheet)

a.child<-as.numeric(sample.sheet[,"Age.at.DX"])<16 | is.na(sample.sheet[,"Age.at.DX"])

sample.sheet[(a.child & a.case),c("SampleProject","ParticipantCode","Age.at.DX","Status")]
sample.sheet[(a.child & a.case),"SampleProject"]<-paste(sample.sheet[(a.child & a.case),"SampleProject"],"Child",sep="-")
sample.sheet[(a.child & a.case),c("SampleProject","ParticipantCode","Age.at.DX","Sequence")]
table(sample.sheet[,"SampleProject"])

                         ## AML                    AML-Child       AML-NotDiagnosis-Child                    Asian-AML              Asian-AML-Child Asian-AML-NotDiagnosis-Child                Asian-Control                      Control 
                         ## 136                           22                           17                           11                            2                            4                           23                          363     
#16                           11                            2                            4                           23
qc[1:5,]
are.related<-(qc[,"sample_A"] != qc[,"sample_B"] ) & as.numeric(qc[,"IBS"])>=0.2
related<-unique(qc[are.related,"sample_A"])

a.true.case<-sample.sheet[,"SampleProject"]=="AML"
a.true.control<-sample.sheet[,"SampleProject"]=="Control"
are.related<-(sample.sheet[,"ParticipantCode"] %in% related)

sample.sheet[(are.related & a.true.case),c("SampleProject","ParticipantCode","Status")]
## 107           AML    AMLM12030PGB Diagnosis # realed to AMLM12PAH030PGB
## 104           AML    AMLM12035D-T Diagnosis
## 105           AML    AMLM12038BJD Diagnosis # AMLM12038BJD	AMLM12PAH038BJD	1.19

##  "95"           "AMLM12035D-T" "AMLM12008CAE" "AMLM12024L-G" "AMLM12008AM"  "1406601-PQ" 

bad<-c("AMLM12PAH030PGB","AMLM12PAH038BJD","AMLM12035D-T","1406601-PQ")
 
are.related<-(qc[,"sample_A"] != qc[,"sample_B"] ) & as.numeric(qc[,"IBS"])>=0.2 & ( !(qc[,"sample_A"] %in% bad) | !(qc[,"sample_B"] %in% bad) )
related<-unique(qc[are.related,"sample_A"])
a.true.case<-sample.sheet[,"SampleProject"]=="AML"
a.true.control<-sample.sheet[,"SampleProject"]=="Control"
are.related<-(sample.sheet[,"ParticipantCode"] %in% related)

sample.sheet[(are.related & a.true.case),c("SampleProject","ParticipantCode")]
##     SampleProject ParticipantCode  Sequence
## 107           AML    AMLM12030PGB Diagnosis
## 104           AML    AMLM12035D-T Diagnosis


sample.sheet[(are.related & a.true.control),c("SampleProject","ParticipantCode","Status")]


a.lung<-grepl("^0413",sample.sheet[,"ParticipantCode"])
a.indg<-grepl("^SKDP.200",sample.sheet[,"ParticipantCode"])
sample.sheet[a.indg,"ParticipantCode"]
# "SKDP-200.3023" "SKDP-200.3024" "SKDP-200.3025" "SKDP-200.3045" "SKDP-200.305"  "SKDP-200.3065" "SKDP-200.3083" "SKDP-200.3084" "SKDP-200.7036" "SKDP-200.7093" "SKDP-200.7095" "SKDP-200.7096"
sample.sheet[a.lung,"ParticipantCode"]
# "0413E1210023"  "0413E1211353"  "0413E14100212" "0413E14113512" "0413E2100212"  "0413E2113512" 

contaminted<-unique(c(sample.sheet[a.lung,"ParticipantCode"],sample.sheet[a.indg,"ParticipantCode"],sample.sheet[(are.related & a.true.control),c("ParticipantCode")]))
length(contaminted) #  148
contaminted<-c(contaminted,bad)

ok<-!(sample.sheet[,"ParticipantCode"] %in% contaminted)

t(table(sample.sheet[ok,"SampleProject"]))

  ##      AML AML-Child AML-NotDiagnosis-Child Asian-AML Asian-AML-Child Asian-AML-NotDiagnosis-Child Asian-Control Control
  ## [1,] 135        22                     16        11               2                            2            13     244

  ##      AML AML-Child AML-Child-Child-Child AML-NotDiagnosis-Child-Child-Child Asian-AML Asian-AML-Child Asian-AML-Child-Child-Child Asian-AML-NotDiagnosis-Child-Child-Child Asian-Control Control
  ## [1,] 132         3                    22                                 16        10               1                           2                                        2            13     447

AML:135
AML-Child:22
AML-NotDiagnosis-Child:1
Asian-AML:11
Asian-AML-Child:2
Asian-AML-NotDiagnosis-Child :13
Asian-Control Control:22
Control:244
Sequence.Number

#posns<-match(sample.sheet.full[,"ParticipantCode"],qc[,"sample_A"])
posns<-match(sample.sheet[,"ParticipantCode"],qc[,"sample_A"])
missing<-is.na(posns)
sum(missing)

sample.sheet[1:5,]
qc[1:5,]
sample.sheet[!missing,"Sex"]<-qc[posns[!missing],"sex_Predicted"]
table(sample.sheet[,"Sex"])
sample.sheet[ sample.sheet[,"Sex"]==-9, "ParticipantCode"]
sample.sheet[ sample.sheet[,"Sex"]==-9, "Sex"]<-2 ### all these were AOGC
sample.sheet[ sample.sheet[,"Sex"]==9, ] #  AMAS-2.3-Relapse2   AMLM12036T-S 
sample.sheet[ sample.sheet[,"Sex"]==9, "Sex"]<-1 # both are male
##   1   2 
## 161 488 

are.related<-(qc[,"sample_A"] != qc[,"sample_B"] ) & as.numeric(qc[,"IBS"])>=0.2
related<-unique(qc[are.related,"sample_A"])
sample.sheet[1:5,]


getwd()
setwd( "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis")
setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis")
write.table(sample.sheet,file="sample_sheet.for.analysis.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(contaminted,file="related.or.bad.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
#########################################################


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
contaminated.file<-c("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/AML_contaminated.txt") # COLL-FAM-1.201,COLL-FAM-1.3
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
#maf.threshold.filter.to.use<-c(0.01)
maf.threshold.filter.to.use<-c(0.001,0.01,0.025)
a.label<-"coding_0.001"
dont.build.summary<-TRUE

sample.types<-c("AK","Control","Normal","PD","SCC" )



###############################################  NEW HC u  with all data
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/SVAnalysis/2015-03-16_AllAMLandLungLI.BEST.chrALL.ACC.large.variants.analysis.rare.txt
## /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/2015-03-16_AllAMLandLung.BEST.chrALL.ACC_SUBSET.ALL.ALL_GENOTYPES_analysis.txt
## /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/2015-03-16_AllAMLandLung.BEST.chrALL.ACC_0.025.ALL.ALL_GENOTYPES_analysis.txt

analysis.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/SVAnalysis"
annotate.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Annotate"

## /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-10-09_AML_CompleteGenomics_HC/Analysis/2014-10-09_AML_CompleteGenomics_HC.chr1.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
##/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-10-09_AML_CompleteGenomics_HC/Analysis/2014-10-09_AML_CompleteGenomics_HC.chrALL.Indel.ALL.ALL_GENOTYPES_analysis.txt

## project.extension<-".analysis-maf-filtered.txt"
## project.name<-"2014-10-09_AML_CompleteGenomics_HC." ## prefix for output file
## fam<-c("chrALL_GENOTYPES") #  ALL or  c() ""-one project (the prefix of the summary files to collect
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-11-04_AML_TCGA_Replication/Analysis/2014-11-04_AML_TCGA_Replication.chrALL.ACC_good_qual.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt

project.extension<-".large.variants.analysis.rare.txt"
project.name<-"2015-03-16_AllAMLandLung." ## prefix for output file
fam<-c("chrALL.ACC") #  ALL or  c() ""-one project (the prefix of the summary files to collect
#fam<-c("BEST.chrALL.ACC_0.025.ALL.ALL_GENOTYPES")
run.per.chromsome<-FALSE # TRue is doing per chromosome

#the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/BAM/TGCM-AML-combine_SampleSheet.csv"
#the.sample.sheet<-"/media/UQCCG/Sequencing/CompleteGenomics/Chort_descriptions/Sequencing comparisons-ver 6.csv"
the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/sample_sheet.for.analysis.txt"
related.or.bad.file<- "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/related.or.bad.txt"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)
remove.cols<-c()

#regions.file<-"/media/scratch2/AOGC-NGS/GFOS/gefos.seq/METHODS/0613-skatmeta-gefos/static/Homo_sapiens.GRCh37.70.protein_coding.genespace_boundaries.5k.split100k.txt"
core.ann<-c("chr","start","end","REF","ALT","TYPE") # out put to annanlsys programs and need foe colun labels
dont.build.summary<-FALSE ##

GATK.SB<-TRUE


a.label<-"CoVarRun.noControl.AML.regions"
maf.threshold.filter.to.use<-c(0.001,0.005,0.01,0.025)
sample.types<-c("Control","AML","AML-Child","AML-NotDiagnosis-Child","Asian-AML","Asian-AML-Child","Asian-AML-NotDiagnosis-Child","Asian-Control")
hw.target<-"Control"
#high.missing.subset<-c("AML.filt","Control.filt")
cancer.group<-"AML"  ## used for no.genotypes.cancer etc
missing.targets<-c("AML.filt","Control.filt")


SnpStats.file.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Snp_read_filter"
SnpStats.file.suffix<-"HC.2.BEST.chrALL.ACC_SUBSET.ALL.ALL_GENOTYPES_analysis.txt_snpstats.stats.txt$"
#SnpStats.bam.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/snpstats_input.txt"
SnpStats.bam.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/snpstats_input_aligner.txt"
SnpStats.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/snp_stats/2015-03-16_AllAMLandLungHC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt_snpstats.read.filter.txt"
SnpStats.stats.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/snp_stats/2015-03-16_AllAMLandLungHC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt_snpstats.stats.txt"

###########################################################################









###############################################  NEW HC u  with all data
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis
## /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/2015-03-16_AllAMLandLung.BEST.chrALL.ACC_SUBSET.ALL.ALL_GENOTYPES_analysis.txt
## /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/2015-03-16_AllAMLandLung.BEST.chrALL.ACC_0.025.ALL.ALL_GENOTYPES_analysis.txt

analysis.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis"
annotate.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Annotate"

## /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-10-09_AML_CompleteGenomics_HC/Analysis/2014-10-09_AML_CompleteGenomics_HC.chr1.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
##/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-10-09_AML_CompleteGenomics_HC/Analysis/2014-10-09_AML_CompleteGenomics_HC.chrALL.Indel.ALL.ALL_GENOTYPES_analysis.txt

## project.extension<-".analysis-maf-filtered.txt"
## project.name<-"2014-10-09_AML_CompleteGenomics_HC." ## prefix for output file
## fam<-c("chrALL_GENOTYPES") #  ALL or  c() ""-one project (the prefix of the summary files to collect
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-11-04_AML_TCGA_Replication/Analysis/2014-11-04_AML_TCGA_Replication.chrALL.ACC_good_qual.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/2015-03-16_AllAMLandLungHC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
project.extension<-"_analysis-maf-filtered.txt"
project.name<-"2015-03-16_AllAMLandLung." ## prefix for output file
fam<-c("HC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES") #  ALL or  c() ""-one project (the prefix of the summary files to collect
#fam<-c("BEST.chrALL.ACC_0.025.ALL.ALL_GENOTYPES")
run.per.chromsome<-FALSE # TRue is doing per chromosome

#the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/BAM/TGCM-AML-combine_SampleSheet.csv"
#the.sample.sheet<-"/media/UQCCG/Sequencing/CompleteGenomics/Chort_descriptions/Sequencing comparisons-ver 6.csv"
the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/sample_sheet.for.analysis.txt"
related.or.bad.file<- "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/related.or.bad.txt"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)
remove.cols<-c()

#regions.file<-"/media/scratch2/AOGC-NGS/GFOS/gefos.seq/METHODS/0613-skatmeta-gefos/static/Homo_sapiens.GRCh37.70.protein_coding.genespace_boundaries.5k.split100k.txt"
core.ann<-c("chr","start","end","REF","ALT","TYPE") # out put to annanlsys programs and need foe colun labels
dont.build.summary<-FALSE ##

GATK.SB<-TRUE


a.label<-"CoVarRun.noControl.AML.regions"
maf.threshold.filter.to.use<-c(0.001,0.005,0.01,0.025)
sample.types<-c("Control","AML","AML-Child","AML-NotDiagnosis-Child","Asian-AML","Asian-AML-Child","Asian-AML-NotDiagnosis-Child","Asian-Control")
hw.target<-"Control"
high.missing.subset<-c("AML.filt","Control.filt")
cancer.group<-"AML"  ## used for no.genotypes.cancer etc
missing.targets<-c("AML.filt","Control.filt")
targets.analysis<-c("AML","Control")  ## used to do comaprion of rare in group ect


SnpStats.file.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Snp_read_filter"
SnpStats.file.suffix<-"HC.2.BEST.chrALL.ACC_SUBSET.ALL.ALL_GENOTYPES_analysis.txt_snpstats.stats.txt$"
#SnpStats.bam.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/snpstats_input.txt"
SnpStats.bam.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/snpstats_input_aligner.txt"
SnpStats.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/snp_stats/2015-03-16_AllAMLandLungHC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt_snpstats.read.filter.txt"
SnpStats.stats.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/snp_stats/2015-03-16_AllAMLandLungHC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt_snpstats.stats.txt"
###########################################################################






###############################################  mixed aligens with AOGC
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis
## /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/2015-03-16_AllAMLandLung.BEST.chrALL.ACC_SUBSET.ALL.ALL_GENOTYPES_analysis.txt
## /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/2015-08-14_AML_mixedAlignersHC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
##/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/2015-08-14_AML_mixedAlignersHC.ALL.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis.txt


analysis.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis"
annotate.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Annotate"

## /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-10-09_AML_CompleteGenomics_HC/Analysis/2014-10-09_AML_CompleteGenomics_HC.chr1.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
##/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-10-09_AML_CompleteGenomics_HC/Analysis/2014-10-09_AML_CompleteGenomics_HC.chrALL.Indel.ALL.ALL_GENOTYPES_analysis.txt

## project.extension<-".analysis-maf-filtered.txt"
## project.name<-"2014-10-09_AML_CompleteGenomics_HC." ## prefix for output file
## fam<-c("chrALL_GENOTYPES") #  ALL or  c() ""-one project (the prefix of the summary files to collect
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-11-04_AML_TCGA_Replication/Analysis/2014-11-04_AML_TCGA_Replication.chrALL.ACC_good_qual.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/2015-03-16_AllAMLandLungHC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt

## project.extension<-"_analysis-maf-filtered.txt"
## project.name<-"2015-08-14_AML_mixedAligners." ## prefix for output file
## fam<-c("HC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES") #  ALL or  c() ""-one project (the prefix of the summary files to collect


project.extension<-"_analysis.txt"
project.name<-"2015-08-14_AML_mixedAligners." ## prefix for output file
fam<-c("HC.ALL.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES") #  ALL or  c() ""-one project (the prefix of the summary files to collect


#fam<-c("BEST.chrALL.ACC_0.025.ALL.ALL_GENOTYPES")
run.per.chromsome<-FALSE # TRue is doing per chromosome

#the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/BAM/TGCM-AML-combine_SampleSheet.csv"
#the.sample.sheet<-"/media/UQCCG/Sequencing/CompleteGenomics/Chort_descriptions/Sequencing comparisons-ver 6.csv"
the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/sample_sheet.for.analysis.txt"
related.or.bad.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/related.or.bad.txt"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)
remove.cols<-c()

#regions.file<-"/media/scratch2/AOGC-NGS/GFOS/gefos.seq/METHODS/0613-skatmeta-gefos/static/Homo_sapiens.GRCh37.70.protein_coding.genespace_boundaries.5k.split100k.txt"
core.ann<-c("chr","start","end","REF","ALT","TYPE") # out put to annanlsys programs and need foe colun labels
dont.build.summary<-FALSE ##

GATK.SB<-TRUE


a.label<-"CoVarRun.noControl.AML.regions"
maf.threshold.filter.to.use<-c(0.001,0.005,0.01,0.025,0.05)

sample.types<-c("Control","AML","AML-Child","AML-NotDiagnosis-Child","Asian-AML","Asian-AML-Child","Asian-AML-NotDiagnosis-Child","Asian-Control") ## these are the disease type
case.control.classes<-c(0,1,9,9,9,9,9,9) ### Must be the same length as above


hw.target<-"Control"
high.missing.subset<-c("AML.filt","Control.filt")
cancer.group<-"AML"  ## used for no.genotypes.cancer etc
missing.targets<-c("AML.filt","Control.filt")
targets.analysis<-c("AML","Control")  ## used to do comaprion of rare in group ect


## SnpStats.file.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Snp_read_filter"
## SnpStats.file.suffix<-"HC.2.BEST.chrALL.ACC_SUBSET.ALL.ALL_GENOTYPES_analysis.txt_snpstats.stats.txt$"
## #SnpStats.bam.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/snpstats_input.txt"
## SnpStats.bam.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/snpstats_input_aligner.txt" ##*** May need updating
## SnpStats.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/snpstats/2015-08-14_AML_mixedAlignersHC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt_snpstats.read.filter.txt"
## SnpStats.stats.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/snpstats/2015-08-14_AML_mixedAlignersHC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt_snpstats.stats.txt"
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/snpstats/2015-08-14_AML_mixedAlignersHC.ALL.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis.txt_snpstats.stats.txt


SnpStats.file.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Snp_read_filter"
SnpStats.file.suffix<-"HC.ALL.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis.plus.EXTRAS.txt_snpstats.stats.txt$"
#SnpStats.bam.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/snpstats_input.txt"
SnpStats.bam.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/snpstats_input_aligner.txt" ##*** May need updating
SnpStats.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/snpstats/2015-08-14_AML_mixedAlignersHC.ALL.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis.plus.EXTRAS.txt_snpstats.read.filter.txt"
SnpStats.stats.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/snpstats/2015-08-14_AML_mixedAlignersHC.ALL.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis.plus.EXTRAS.txt_snpstats.stats.txt"
###########################################################################






################################## Other input files needed - path required in file names
num.cores<-5 # use for filtering don't make too large else will exceed memory

gene.symbol.file.for.clusters<-"/media/UQCCG/Software/annovar/humandb/Gene_symbol_aliases.txt"

cluster.definition.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/clusters/Clusters for Skat analysis.txt"
#cluster.definition.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/Final_FANC_clusters.csv" ## tab delimited with header
ensembleID.to.HGNC.map.file<-"/media/UQCCG/Sequencing/Data/Genomes/hg19/ENSG_to_HGNC.txt" # tab delimited with header
coverage.file<-"/media/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_SAMPLE_Mon_Jul_20_2015.txt" # tab delimited with header


## #/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/snp_stats/2015-03-16_AllAMLandLungHC.2.BEST.chrALL.ACC_SUBSET.ALL.ALL_GENOTYPES_analysis.txt_snpstats.stats.txt
## SnpStats.file.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Snp_read_filter"
## SnpStats.file.suffix<-"HC.2.BEST.chrALL.ACC_SUBSET.ALL.ALL_GENOTYPES_analysis.txt_snpstats.stats.txt$"
## #SnpStats.bam.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/snpstats_input.txt"
## SnpStats.bam.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/snpstats_input_aligner.txt"
## SnpStats.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/snp_stats/2015-03-16_AllAMLandLungHC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt_snpstats.read.filter.txt"
## SnpStats.stats.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/snp_stats/2015-03-16_AllAMLandLungHC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt_snpstats.stats.txt"

BWA.bad<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/containamination loci.csv",header=T,sep="\t",fill=TRUE,skip=0,stringsAsFactors=FALSE)
 

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
source("Sliding_window_generation.r")
source("hwe.r")
setwd(analysis.dir)
##########################################################################
########################### START ADDITIONAL FILTERING " : bad.qual.locations############################
##########################################################################
##########################################################################
## extra filtering
#SNP_STATS:version1
#FILTER_SUMMARY:alt_mismatch5_proportion_gt_0.8;Strand_bias_proportion_gt_0.9;Read_end3_proportion_gt_0.8;High_not_called_alts_gte_3
#SUMMARY_CALLED:ref_count;alt_count;mismatch_alt_count;fstrand_alt_count;rstrand_alt_count;read_end_alt_count;prop_mismatch;prop_fstrand;prop_rstrand;prop_read_end
#SUMMARY_NOT_CALLED:ref_count;alt_count;mismatch_alt_count;fstrand_alt_count;rstrand_alt_count;read_end_alt_count;high_alt_count
#SAMPLE.RF:GATK_called_snp;ref_count;alt_count;mismatch_alt_count;fstrand_alt_count;rstrand_alt_count;read_end_alt_count


if(is.null(SnpStats.bam.file) | SnpStats.bam.file==""){
files<-dir(SnpStats.file.dir)
files<-files[grepl(SnpStats.file.suffix,files)]
files
i<-1
setwd(SnpStats.file.dir)
for (i in 1:length(files)){
a.filt<-read.delim(files[i],header=T,sep="\t",fill=TRUE,skip=5,stringsAsFactors=FALSE)
if(i==1){
  filt<-a.filt
}else{
  filt<-rbind(filt,a.filt)
}

}

}else{
  print(paste0("Reading SNPstats file",SnpStats.file))
  filt<-read.table(SnpStats.file,header=T,skip=5,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
  filt<-filt[,1:11]
}

## filt[1:5,"chr"]
## table(filt[,"chr"]
      )

## filt[ filt[,"chr"]=="chr10" ,][1:20,]
do.chr<-gsub("^chr","",filt[,"chr"])
table(do.chr)
filt[,"chr"]<-do.chr

filt.key<-build.key(filt,core.ann,add.chr.label=TRUE)
tail(filt.key)

setwd(analysis.dir)

colnames(filt)
pass.filt<-strsplit(filt[,"FILTER_SUMMARY"],split=";")

pass.filt<-unlist(lapply(pass.filt,function(x) {sum(as.logical(x[1:3]))==0}))
#pass.filt[173713]
sum(!pass.filt)
#filt.key[173713]
snp.fail.filt<-filt.key[!pass.filt]
snp.fail.filt[1:5]

#filt<-filt[!pass.filt,]

if(exists("SnpStats.stats.file")){
a.indel.stats<-read.table(SnpStats.stats.file,header=T,skip=5,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
a.indel.stats[1:5,1:12]
dim(a.indel.stats)

do.chr<-gsub("^chr","",a.indel.stats[,"chr"])
table(do.chr)
a.indel.stats[,"chr"]<-do.chr


key.stats<-build.key(filt,core.ann,add.chr.label=TRUE)
key.stats[1:5]
}
## all.possible.samples.stats<-gsub(".FAD$","",colnames(a.indel.stats)[grep(".FAD$",colnames(a.indel.stats))],perl=TRUE)
## all.possible.samples.stats
## length(all.possible.samples)
## length(all.possible.samples.stats)
## all.possible.samples[!(all.possible.samples %in% all.possible.samples.stats)]
##########################################################################
########################### END ADDITIONAL FILTERING ############################
##########################################################################
a.indel.stats[1:5,1:15]


contaminated.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/contaminated_AOGC_SEQ_samples.txt"
contaminated<-read.table(contaminated.file,header=F,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

#related.or.bad.file<- "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/related.or.bad.txt"
related.or.bad<-read.table(related.or.bad.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

if(SnpStats.bam.file!="" | exists(SnpStats.bam.file)) {
SnpStats.bam<-read.delim(SnpStats.bam.file,header=F,sep="\t",fill=TRUE,skip=0,stringsAsFactors=FALSE,check.names=FALSE)
# SnpStats.bam[1:5,]
}




###########################
#/media/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_SAMPLE_Fri_Feb_06_2015.txt  latest
core.ann<-c("chr","start","end","REF","ALT","TYPE") # out put to annanlsys programs and need foe colun labels
dont.build.summary<-TRUE ##
GATK.SB<-TRUE
#max(BWA.bad[,"P.value"])




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
#mafs<-colnames(a.indel)[grepl("maf",colnames(a.indel))]
############################################# POPULATION MAF FILTER - PART A
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER

maf.threshold<-0.0  #  MAF threshold for annovar calling zero useful to get back all results !!do not modify!!

maf.threshold.filter.to.use<-sort(as.numeric(maf.threshold.filter.to.use))

filter.cols.novel.use<-c("PopFreqMax","NHBLI_6500_ANNOVAR_ALL","NHBLI_6500_ALL","NHLBI_5400_ALL","NHLBI_5400_EUR","NHLBI_5400_AFR","1000genome","1000genome_asian","1000genome_mine","snp141","snp141_clinical","snp137","CG69","EUR_ASN_AFR_INDEL","AOGC-NGS_ALL","AOGC-NGS_ALL_OLD","Chinese") ##
#filter.cols.maf.use<-c("NHBLI_6500_ANNOVAR_ALL","NHBLI_6500_ALL","NHBLI_6500_EA","NHBLI_6500_AA","NHLBI_5400_ALL","1000genome","snp141","snp137","snp135","AOGC-NGS_ALL")
filter.cols.maf.use<-c("NHBLI_6500_ANNOVAR_ALL","NHBLI_6500_ALL","NHBLI_6500_EA","NHBLI_6500_AA","NHLBI_5400_ALL","1000genome","snp141","snp137","snp135")  # don't use AOGC if they are controls
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

interesting.to.prefilter<-c("UTR3","UTR5","UTR5;UTR3","snoRNA","snRNA","antisense","sense_intronic","ncRNA_exonic","ncRNA_splicing","intronic") #use gerp.score.threshold

extra.vep.annotations<-c("Uploaded_variation","Gene","Feature","Protein_position","Amino_acids")

vep.types<-c("not_assigned","stop_gained","stop_lost","missense_variant","splice_acceptor_variant","splice_donor_variant","splice_region_variant","initiator_codon_variant","stop_retained_variant","incomplete_terminal_codon_variant","frameshift_variant","inframe_deletion","inframe_insertion","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_exon_variant","NC_stop_gained","NC_stop_lost","NC_splice_acceptor_variant","NC_splice_donor_variant","NC_splice_region_variant","NC_initiator_codon_variant","NC_stop_retained_variant","NC_non_coding_exon_variant","NC_incomplete_terminal_codon_variant","NC_3_prime_UTR_variant","mature_miRNA_variant","NC_5_prime_UTR_variant","TF_binding_site_variant","TFBS_ablation","TFBS_amplification","regulatory_region_variant","intron_variant","NC_intron_variant","synonymous_variant","coding_sequence_variant","NC_synonymous_variant","upstream_gene_variant","downstream_gene_variant","intergenic_variant","NC_intergenic_variant","NMD_transcript_variant","nc_transcript_variant","NC_nc_transcript_variant","feature_truncation","feature_elongation")

## vep.types<-c( "not_assigned","stop_gained","stop_lost","stop_lost,NMD_transcript_variant","stop_gained,splice_region_variant,NMD_transcript_variant","initiator_codon_variant,splice_region_variant","splice_region_variant,3_prime_UTR_variant","stop_gained,NMD_transcript_variant","missense_variant,splice_region_variant","missense_variant","splice_acceptor_variant","splice_acceptor_variant,nc_transcript_variant","splice_region_variant,3_prime_UTR_variant,NMD_transcript_variant","splice_donor_variant,nc_transcript_variant","splice_region_variant,intron_variant,NMD_transcript_variant","splice_donor_variant","splice_region_variant","splice_region_variant,5_prime_UTR_variant","splice_region_variant,synonymous_variant","splice_region_variant,intron_variant,nc_transcript_variant","splice_region_variant,non_coding_exon_variant,nc_transcript_variant","missense_variant,NMD_transcript_variant","splice_region_variant,intron_variant","NMD_transcript_variant","intron_variant,NMD_transcript_variant","mature_miRNA_variant","5_prime_UTR_variant","5_prime_UTR_variant,NMD_transcript_variant","non_coding_exon_variant,nc_transcript_variant","3_prime_UTR_variant,NMD_transcript_variant","non_coding_exon_variant","TF_binding_site_variant","intron_variant,nc_transcript_variant","synonymous_variant,NMD_transcript_variant","3_prime_UTR_variant","regulatory_region_variant","upstream_gene_variant","downstream_gene_variant","intergenic_variant","intron_variant","synonymous_variant")


vep.coding<-c("not_assigned","stop_gained","stop_lost","missense_variant","splice_acceptor_variant","splice_donor_variant","splice_region_variant","initiator_codon_variant","stop_retained_variant","incomplete_terminal_codon_variant","frameshift_variant","inframe_deletion","inframe_insertion")

vep.noncoding<-c("5_prime_UTR_variant","3_prime_UTR_variant","non_coding_exon_variant","NC_stop_gained","NC_stop_lost","NC_splice_acceptor_variant","NC_splice_donor_variant","NC_splice_region_variant","NC_initiator_codon_variant","NC_stop_retained_variant","NC_non_coding_exon_variant","NC_incomplete_terminal_codon_variant","NC_3_prime_UTR_variant","mature_miRNA_variant","NC_5_prime_UTR_variant","TF_binding_site_variant","TFBS_ablation","TFBS_amplification","regulatory_region_variant","intron_variant")              

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
seq.type<-read.delim(seq.type.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
coverage<-seq.type # read.delim("/media/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_SAMPLE_Tue_Oct_14_2014.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

#to.fix<-(coverage[,"Project"]=="SDDS")
#coverage[to.fix,"Sample"]<-gsub("^0+","", coverage[to.fix,"Sample"])  ## SDDS project has proceeding zeros


seq.type[1:5,]
seq.type.sub<-seq.type[seq.type[,"Project"]=="AOGC-NGS",]
nim.samples<-seq.type.sub[seq.type.sub[,"Capture.Method"]=="TruD:NimX","Sample"]
ill.samples<-seq.type.sub[seq.type.sub[,"Capture.Method"]=="TruD:TruX","Sample"]

nim.samples<-paste(nim.samples,"GT",sep=".")
ill.samples<-paste(ill.samples,"GT",sep=".")


length(nim.samples)
length(ill.samples)

############################################ Nimblegen and illuma capture loci ####################################################

## ########check

the.sample.sheet

sample.sheet.full<-read.delim(the.sample.sheet,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) #,check.names=FALSE )
sample.sheet.full[1:5,1:10]
colnames(sample.sheet.full)
dim(sample.sheet.full)
table(sample.sheet.full[,"SampleProject"])
table(sample.sheet.full[,"Sex"])
##### fix 0 and 9 for missing to NA
 sample.sheet.full[sample.sheet.full[,"SampleProject"]=="AML", c("SampleProject","SAMPLE","Complex")]
 sample.sheet.full[sample.sheet.full[,"SampleProject"]=="AML", c("SampleProject","SAMPLE","Complex")]

table(sample.sheet.full[sample.sheet.full[,"SampleProject"]=="AML", c("Complex")])

pca<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/SNPs/pca_aml.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)


pca[1:5,]
sample.sheet.full[1:5,]

posns<-match(sample.sheet.full[,"ParticipantCode"],pca[,"Sample"])
missing<-is.na(posns)
sum(missing)
sample.sheet.full[missing,"ParticipantCode"]

extra<-pca[posns,]
sample.sheet.full<-cbind(sample.sheet.full,extra)
sample.sheet.full[1:5,]

select<-as.numeric(sample.sheet.full[,"PCA.1"])< -0.005
sum(select)

table(sample.sheet.full[select,"SampleProject"])
table(sample.sheet.full[!select,"SampleProject"])

sample.sheet.full[!select,c("SampleProject","Sample")]
#write.table(sample.sheet.full,file="sample_sheet.for.analysis_wPCA.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
### > 
##                    AML              AML-Child AML-NotDiagnosis-Child                Control 
##                    131                     17                     13                    436 
## > table(sample.sheet.full[!select,"SampleProject"])

##                    Asian-AML              Asian-AML-Child Asian-AML-NotDiagnosis-Child                Asian-Control 
##                           15                            7                            5                           25 
################################################# coverage check
################################################# coverage check
######### Check and fix the sample sheet


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

table(coverage[, "Capture.Method"])

#SnpStats.bam[1:5,]

if(exists("coverage")){

posns<-match(sample.sheet.full[,"ParticipantCode"], coverage[,"Sample"])
missing<-is.na(posns)
sum(missing)
seq.project<-coverage[posns, c("Project","Capture.Method")]
capture<-coverage[posns, "Capture.Method"]
table(capture)
capture<-strsplit(capture,split=";")
lengths<-unlist(lapply(capture,length))
capture[lengths==0]<-"unknown"

capture<-unlist(lapply(capture,function(x) { paste(unique(x),collapse=";")})) #x[(is.na(x) | x=="")]<-"NONE";
table(capture)
#cbind(capture,seq.project)
 
sample.sheet.full<-cbind(sample.sheet.full,capture,seq.project)
}
table(sample.sheet.full[,"capture"])
table(sample.sheet.full[,"capture"])


if(exists("SnpStats.bam")){
posns<-match(sample.sheet.full[,"ParticipantCode"],SnpStats.bam[,2])
missing<-is.na(posns)
sum(missing)
Aligner<-SnpStats.bam[posns,3]
sample.sheet.full[is.na(Aligner),"ParticipantCode"] ### these wree all aogc
Aligner[is.na(Aligner)]<-"novoalign"
sample.sheet.full<-cbind(sample.sheet.full,Aligner)
}
table(sample.sheet.full[,"Aligner"])
table(sample.sheet.full[,"Aligner"])
colnames(sample.sheet.full)
          ## bwa     novoalign novoalign;bwa 
          ##  81           557            11
################################################# Build snpstats file
################################################# coverage snpstats file





## a.control<-grepl("Control",sample.sheet.full[,"SampleProject"])
## a.case<-grepl("AML",sample.sheet.full[,"SampleProject"])
## AFF_STATUS<-rep(NA,times=dim(sample.sheet.full)[1])
## sum(a.case & a.control) # should be 0
## AFF_STATUS[a.case]<-2
## AFF_STATUS[a.control]<-1

## #sample.sheet.full[aml.have & a.control,]

## posns<-match(sample.sheet.full[,"ParticipantCode"], coverage[,"Sample"])
## missing<-is.na(posns)
## sum(missing)

## sample.sheet.full.coverage<-coverage[posns,]
## sample.sheet.full.coverage[1:5,]

## SnpStats.input.file<-cbind(paste("/mnt/UQCCG/Sequencing/Projects/",sample.sheet.full.coverage[,"Project"],"/BAM/",sample.sheet.full.coverage[,"BAM"],sep=""),sample.sheet.full.coverage[,"ID"],AFF_STATUS,sample.sheet.full[,"ParticipantCode"])
## SnpStats.input.file[1:5,]

## colnames(SnpStats.input.file)<-c("BAM","SAMPLE_ID","AFF_STATUS","ParticipantCode")

## all.bams<-SnpStats.input.file[,"SAMPLE_ID"]


## mult.bams<-strsplit(all.bams,split=";")
## num.mult.bams<-unlist(lapply(mult.bams,length))
## first.bam<-unlist(lapply(mult.bams,function(x) x[1]))
## #has.mult.bams
## SnpStats.input.file[,"SAMPLE_ID"]<-first.bam
## SnpStats.input.file<-cbind(SnpStats.input.file,num.mult.bams,all.bams)
## SnpStats.input.file[1:5,]


## getwd()
## write.table(SnpStats.input.file,file="SnpStats.input.file.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## ## /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/2015-03-16_AllAMLandLung.chr13.ALL.ALL_GENOTYPES_analysis.txt



##  [1] "860"  "861"  "862"  "2529" "2530" "2531" "2567" "2688" "2689" "2851" "2924" "2944" "3198" "3199" "3200" "3204" "3255" "3256" "3257" "3258" "3259" "3342" "3343" "3344" "3345" "2568" "2569" "2669" "2670" "2671" "3040" "3041" "3169" "3170" "3211"
## [36] "3212" "3282" "3287" "3288" "3289" "3291" "3293" "3294" "3295" "3309"


################################################# coverage check
################################################# coverage check



########################################################### define bad samples
#bad.samples<-sample.sheet.full[ (sample.sheet.full[,"SampleProject"]!="AML" & sample.sheet.full[,"SampleProject"]!="Control"  )   ,"ParticipantCode"] # "1"  "9"  "63" "74" "83" "84" "99"
#table(sample.sheet.full[ !(sample.sheet.full[,"ParticipantCode"] %in% bad.samples) ,"SampleProject"])

    ## AML Control 
    ## 136     363
## tail(sample.sheet.full)
## test<-((sample.sheet.full[,"Race"]!="W" | is.na(sample.sheet.full[,"Race"])) & sample.sheet.full[,"AffectionStatus"]==2) 
## sample.sheet.full[test,]
## bad.samples<-sample.sheet.full[test,"ParticipantCode"] # "1"  "9"  "63" "74" "83" "84" "99"

contaminated.aogc<-sample.sheet.full[sample.sheet.full[,"ParticipantCode"] %in% contaminated[,1],"ParticipantCode"]
related.or.bad[,1]
bad.samples<-unique(c(contaminated.aogc,related.or.bad[,1]))

table(sample.sheet.full[ !(sample.sheet.full[,"ParticipantCode"] %in% bad.samples) ,"SampleProject"])
################################################################################









pheno.types<-c("SampleProject") ## vales is column header
names(pheno.types)<-c("SampleProject") ### name is output columns

## case.control.classes<-c(0,1,0)
## names(case.control.classes)<-c("Control","NMD","AOGC")
## case.control.classes


case.control<-c("SampleProject")
sample.types
# case.control.classes<-c(0,1,9,9,9,9,9,9)
names(case.control.classes)<-sample.types
case.control.classes
# ib<-1
for(ib in 1:length(case.control)){
  if(!(case.control[ib] %in% colnames(sample.sheet.full))){next}
  sample.sheet.full[(  !(sample.sheet.full[,case.control[ib]] %in% names(case.control.classes))  |  is.na(sample.sheet.full[,case.control[ib]]) | sample.sheet.full[,case.control[ib]]==0 | sample.sheet.full[,case.control[ib]]==9)  ,case.control[ib]]<-NA
}



                         ## AML                    AML-Child       AML-NotDiagnosis-Child                    Asian-AML              Asian-AML-Child Asian-AML-NotDiagnosis-Child                Asian-Control                      Control 
                         ## 136                           22                           17                           11                            2                            4                           23                          363

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


## load("2015-03-16_AllAMLandLung.BEST.chrALL.ACC_0.025.ALL.ALL_GENOTYPES_analysis.txt_SUBSET.RData")
## column.labels<-read.delim(project.files[ichr],header=F,nrows=1,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
## length(column.labels)
## dim(a.indel)
## colnames(a.indel)<-column.labels
## load("2015-03-16_AllAMLandLung.BEST.chrALL.ACC_0.025.ALL.ALL_GENOTYPES_analysis.txt.RData")
## column.labels<-colnames(a.indel)
## a.indel<-as.matrix(a.indel)
## wanted.genes<-c("DDX41","TET2", "GATA2", "ASXL1", "NOTCH1", "IDH1", "JAK2","WT1","MLL","KRAS","FLT3","IDH2","IDH1","TP53","KIT","NPM1","JAK2","DNMT3A","TET2","RUNX1","NRAS","CEBPA","PTPN11","U2AF1","SMC1A","SMC3","PHF6","STAG2","RAD21","FAM5C","EZH2","HNRNPK","FANCA","FANCB","FANCC","FANCD1","FANCD2","FANCE","FANCF","FANCG","FANCI","BRIP1","FANCL","FANCM","PALB2","RAD51C","SLX4","ERCC4","APITD1","STRA13","C1orf86","C19orf40","C17orf70","SLX1","MUS81","ERCC1","FAN1","EME1","EME2","MRE11A","NBN1","RAD50","FAND1","BRCA1","BARD1","RAD51","RAD51B","RAD51D","XRCC2","XRCC3","RMI1","RMI2","BLM","TOP3A","RPA1","RPA2","RPA3","ATM","ATR","ATRIP","CHECK1","RAD9A","RAD17","CS","DLAT","DLD","DLST","FH","IDH1","IDH2","IDH3A","IDH3B","IDH3G","MDH1","MDH2","ACLY","ACO1","OGDH","ACO2","PC","PCK1","PCK2","PDHA1","PDHA2","PDHB","OGDHL","SDHA","SDHB","SDHC","SDHD","SUCLG2","SUCLG1","SUCLA2")
## wanted.genes2<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/Gene_Lists/mam_new/ALL_genes.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="")
## wanted.genes<-unique(c(wanted.genes,wanted.genes2[,1]))
## wanted<-test[,"Gene.Names"] %in% wanted.genes
## sum(wanted)
## a.indel<-test[wanted,]
## save(list=c("a.indel"),file="2015-03-16_AllAMLandLung.BEST.chrALL.ACC_0.025.ALL.ALL_GENOTYPES_analysis.txt_SUBSET.RData")

all.possible.samples<-gsub(".GT$","",colnames(a.indel)[grep(".GT$",colnames(a.indel))],perl=TRUE)
all.possible.samples.stats<-gsub(".FAD$","",colnames(a.indel.stats)[grep(".FAD$",colnames(a.indel.stats))],perl=TRUE)
length(all.possible.samples)
length(all.possible.samples.stats)

print(paste0("Samples missing in genotyping vs Snpstats:",all.possible.samples[!(all.possible.samples %in% all.possible.samples.stats)]))
#print(all.possible.samples[!(all.possible.samples %in% all.possible.samples.stats)]) #

pheno.types


the.chr<-a.indel[1,"chr"]
print(paste("Doing Chromosome ",the.chr))
print(table(a.indel[,"chr"]))

if(!grepl("^chr",the.chr)){
a.indel[,"chr"]<-paste("chr",a.indel[,"chr"],sep="")
}
key<-build.key(a.indel,core.ann)
rownames(a.indel)<-key
length(key)
key[1:5]
key.stats[1:5]
####### REMOVE BAD SAMPLES

## key.stats<-build.key(a.indel.stats,core.ann[1:5],add.chr=TRUE)
## key<-build.key(a.indel,core.ann[1:5])
## length(key)
## length(key.stats)  ##/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/snpstats/2015-08-14_AML_mixedAlignersHC.ALL.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis.txt_snpstats.stats.txt
## sum(key!=key.stats)
## cbind(key,key.stats)[key!=key.stats,][1:10,]


posns<-match(key,key.stats)
missing<-is.na(posns)
sum(missing)

#  write.table(a.indel[missing,],file=paste(project.name,"EXTRAS_NEEDED.txt",sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=TRUE)

## an.indel<-grepl("^indel",a.indel[,"TYPE"])
## a.indel[missing & !an.indel ,core.ann][1:10,]



a.indel.stats<-a.indel.stats[posns,]

rownames(a.indel.stats)<-key

dim(a.indel)
dim(a.indel.stats)





if(length(bad.samples)>0){

bad.samples.labels<-expand.labels.to.samples(bad.samples,c("GT","AD","DP","GQ"),paste.after=TRUE)
if(length(bad.samples.labels)>1){
a.indel<-a.indel[,colnames(a.indel)[!(colnames(a.indel) %in% bad.samples.labels)]]
}

bad.samples.labels<-expand.labels.to.samples(bad.samples,c("FAD","TAD","DUP"),paste.after=TRUE)
if(length(bad.samples.labels)>1){
a.indel.stats<-a.indel.stats[,colnames(a.indel.stats)[!(colnames(a.indel.stats) %in% bad.samples.labels)]]
}

} # have bad samples so remove


#########################################################################################################################
#########################################################################################################################
#########################################################################################################################


########################################################################################################################
####################################### EXTRA snp filtering ##################################################
########################alignemnt issues files
BWA.bad[1:5,]
BWA.bad[1:5,c(core.ann,"gene","P.value")]

## getwd()
## write.table(BWA.bad[,c(core.ann,"gene","P.value")],file="Discordane_positions_BWA_NOVOALIGN.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
#sum(poss.model[pass,"P-value"]<1e-5,na.rm=TRUE)
#poss.model[1:5,]  ### used in conjection with additional snp filtering
 # 2*(pnorm(abs(c(1:8)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)) #2*(pnorm(abs(seq(from=3,to=5,by=0.1)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
# 2*(pnorm(abs(6), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)) ##4 sd greater than mean

posns<-match(key,BWA.bad[,"gene"])
poss.model<-BWA.bad[posns, c("gene","P.value")]
pass.possonian.control.model<- poss.model[,"P.value"]> 2*(pnorm(abs(6), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))  | is.na( poss.model[,"P.value"])

length(pass.possonian.control.model)  ## these are known BWA problem sites
#length(pass)
sum(!pass.possonian.control.model)



########################alignement issues files from filt array

bad.qual.locations<-key %in% snp.fail.filt
sum(bad.qual.locations)
length(bad.qual.locations)
#length(pass)


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


#grep("NOTCH1",names(all.genes))

common.hit.genes<-{}
## all.genes["SCN2A"]

###############################################

#########################################################################################################################

snpinfo.raw<-cbind(key,a.indel[,"Gene.Names"],a.indel[,"Gene.Names"])
snpinfo.raw[1:5,]
tail(snpinfo.raw)
colnames(snpinfo.raw)<-c("Name","gene","cluster")

            poly.gene.site<-grep("::",snpinfo.raw[,"gene"])

            length(poly.gene.site) #### if 
            ipoly<-2
            all.extra<-{}
            if(length(poly.gene.site)>0){
               
                for( ipoly in 1 :length(poly.gene.site)){
                    gene<-unlist(strsplit(snpinfo.raw[poly.gene.site[ipoly],"gene"],split="::") )
                    if(is.null(dim(all.extra))){
                    all.extra<-cbind(snpinfo.raw[poly.gene.site[ipoly],"Name"],gene,gene)
                    colnames(all.extra)<-colnames(snpinfo.raw)
                     }else{
                    extra<-cbind(snpinfo.raw[poly.gene.site[ipoly],"Name"],gene,gene)
                    colnames(extra)<-colnames(snpinfo.raw)
                    all.extra<-rbind(all.extra,extra)
                }
                       
 
            } }# has a poly

# all.extra[1:5,]
        

snpinfo.raw<-rbind(snpinfo.raw,all.extra)


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


ic<-16
recode<-{}
recode.record<-{}
for(ic in 1:length(clusters.wanted)){
    print(ic)
   cluster.genes<-clusters[,clusters.wanted[ic]]
   cluster.genes<-cluster.genes[cluster.genes!="" | is.na(cluster.genes)]
   cluster.genes<-unique(cluster.genes)
#   all.genes[1:5]
   missing.name<-!(cluster.genes %in% names(all.genes))
   if(sum( missing.name)>0){
     posns<-match(cluster.genes[missing.name],gene.aliases[, "Aliases"])
     missing<-is.na(posns)
     

     
     if(sum(missing)>0){
       print("--------------------------------------------------------")
       print(paste("in cluster",clusters.wanted[ic], "PERMANENTLY missing"))
       print(cluster.genes[missing.name][missing])
       print("--------------------------------------------------------")


       
     }
     recode<-cbind(cluster.genes[missing.name][!missing],gene.aliases[posns[!missing], "Approved.Symbol"])

     print(paste("recoded"))
     print(recode)

     
     colnames(recode)<-c("old","new")
 ###### transfer to new gene lists
     posns<-match(clusters[,clusters.wanted[ic]],recode[,"old"])
     missing<-is.na(posns)

     if(is.null(dim(recode.record))){recode.record<-recode}else{recode.record<-rbind(recode.record,recode)}

    clusters[!missing,clusters.wanted[ic]]<-recode[posns[!missing],"new"] ### redefine the clusters
  }

   
  }
 #########################################################    
     



snpinfo<-snpinfo.raw
ic<-5
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

snpinfo.ori[1:20,]

### HAVE
# snpinfo.raw (original from a.indel)
# snpinfo # with extra clusters
# snpinfo.raw a permanent copy of snpinfo

## "FANCM " "MHF1"   "MHF2"   "FAAP24"
clusters[,1]
chk<-apply(clusters,2,function(x){ length(x[x!=""])})
write.table(clusters,file="clusters_as.using.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

recode.record<-recode.record[!duplicated(recode.record[,"old"]),]
write.table(recode.record,file="Recoded_genes.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
############################
#save(list=c("a.indel"),file="indels.RData")


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
 covars<-c("PCA.1","PCA.2","PCA.3","PCA.4") #  covars<-c("1") # c("AGE_SCAN","PCA1","PCA2","PCA3","PCA4") #AGE_SCAN,PCA1,PCA2,PCA3,PCA4 #covars<-c("1")
}else{
covars<-c("1")
}


formula<-paste(target.pheno.col,"~",paste(covars,collapse="+"),sep="")
print(formula)
formula<-formula(formula)


###############################################################

############################ subset samples with phenotype and covars - assume traits and covars need same phenotypes
# sample.sheet.full<-sample.sheet.full[,1:31]
sample.sheet.full[1:2,]
colnames(sample.sheet.full)
sample.sheet.full[1:5,]
table(sample.sheet.full[,"SampleProject" ])

if ( !("SAMPLE" %in% colnames(sample.sheet.full))){
  SAMPLE<- sample.sheet.full[,"ParticipantCode"]
  sample.sheet.full<-cbind(sample.sheet.full,SAMPLE)
}

#if(sum(covars==1) & length(covars)==1){
 got.all.covars<-rep(TRUE,times=dim(sample.sheet.full)[1])
#}else{
#got.all.covars<-apply(sample.sheet.full[,covars],1,function(x) (sum(is.na(as.numeric(x))) ==0 ))
#}
all.possible.samples[(!(all.possible.samples %in% sample.sheet.full[,"SAMPLE"] ) )]
subset.samples<- (sample.sheet.full[,"SAMPLE"] %in% all.possible.samples ) & !is.na(sample.sheet.full[,target.pheno.col]) & sample.sheet.full[,target.pheno.col]!="NA"  & got.all.covars

    
sum(subset.samples)

pheno.ori<-sample.sheet.full[subset.samples ,] ## pheno.ori only contains SAMPLES that have a phenotype
colnames(pheno.ori)
if(!("SAMPLE" %in% colnames(pheno.ori))){ ## add sample column if using PATIENT
  SAMPLE<-pheno.ori[,"PATIENT"]
  pheno.ori<-cbind(SAMPLE,pheno.ori)
 #  colnames(pheno.ori)[colnames(pheno.ori)=="PATIENT"]<- "SAMPLE"
}
colnames(pheno.ori)
print(dim(pheno.ori))
print(paste("Number Samples:",dim(pheno.ori)[1]))

dim(pheno.ori)
dim(sample.sheet.full)
length(all.possible.samples)

########## make sure sample sheet has no extras



pheno.ori[1:5,]



the.samples<-paste(pheno.ori[,"SAMPLE"],"GT",sep=".")  ## samples same order as in pheno.ori
print(paste("Number samples: ",length(the.samples),sep=""))

#seq.type[1:57,]
posns<-match(pheno.ori[,"SAMPLE"],seq.type[,"Sample"])
missing<-is.na(posns)
sum(missing)
pheno.ori[missing,1:5]

if(!exists("pheno.ori")){
pheno.ori<-pheno
}
#pheno<-pheno.ori

#grep("AMLM12005R-G",seq.type[,"Sample"])

### Troels will update but all are Nextera
## capture<-seq.type[posns,"Capture.Method"]
## table(capture) ### all illume here
## pheno<-cbind(pheno,capture)
##       SampleProject FamilyCode SAMPLE PaternalID MaternalID
## NA.8        Control        ALL   3170          0          0
## NA.10       Control        ALL   3212          0          0
## NA.15       Control        ALL   3291          0          0
## NA.17       Control        ALL   3294          0          0
## NA.18       Control        ALL   3295          0          0
## NA.19       Control        ALL   3309          0          0

################################################################################

table(pheno.ori$SampleProject)

                         ## AML                    AML-Child       AML-NotDiagnosis-Child                    Asian-AML              Asian-AML-Child Asian-AML-NotDiagnosis-Child                Asian-Control                      Control 
                         ## 131                           17                           13                           15                            7                            5                           25                          429 

table(pheno.ori$Project)

    ## AMAS AOGC-NGS     MODY RSGB_AML     SDDS     SKDP TGCM-AML 
    ##   42      323       20       50       86       25       96
    
table(pheno.ori$capture)

## sum(!(pheno.ori.use[,"SAMPLE"]==pheno.ori[,"SAMPLE"]))
## the.projects<-c("cancer","Control","normal","PD","SCC","AK")

## cancer<-rep(FALSE,times=dim(pheno.ori)[1])
## cancer[pheno$SampleProject %in% c("AK","SCC")]<-TRUE

## normal<-rep(FALSE,times=dim(pheno.ori)[1])
## normal[pheno.ori$SampleProject %in% c("Normal")]<-TRUE

## Control<-rep(FALSE,times=dim(pheno.ori)[1])
## Control[pheno.ori$SampleProject %in% c("Control")]<-TRUE

## PD<-rep(FALSE,times=dim(pheno.ori)[1])
## PD[pheno.ori$SampleProject %in% c("PD")]<-TRUE

## SCC<-rep(FALSE,times=dim(pheno.ori)[1])
## SCC[grepl("_SCC$",pheno.ori$SAMPLE)]<-TRUE

## AK<-rep(FALSE,times=dim(pheno.ori)[1])
## AK[grepl("_AK1$",pheno.ori$SAMPLE) | grepl("_AK2$",pheno.ori$SAMPLE)]<-TRUE


## the.SCC<-pheno.ori$SAMPLE[grepl("_SCC$",pheno.ori$SAMPLE)]
## the.AK<-pheno.ori$SAMPLE[grepl("_AK1$",pheno.ori$SAMPLE) | grepl("_AK2$",pheno.ori$SAMPLE)]

## sum(PD)
## sum(normal)
## ####"LPH-001-27 Blood ok
## ## "LPH-001-27_PD"
## ## normal[pheno.ori$SAMPLE %in% c("LPH-001-27_PD")]<-TRUE
## ## PD[pheno.ori$SAMPLE %in% c("LPH-001-27_PD")]<-FALSE


## pheno.ori<-cbind(pheno.ori,cancer,Control,normal,PD,SCC,AK)

## pheno.ori[pheno.ori[,the.projects[1]],"SAMPLE"]
## pheno.ori[pheno.ori[,the.projects[2]],"SAMPLE"]
## pheno.ori[pheno.ori[,the.projects[3]],"SAMPLE"]
## pheno.ori[pheno.ori[,the.projects[4]],"SAMPLE"]

##NOVO
## [1] "cancer Num. samples: 74"
## [1] "Control Num. samples: 133"
## [1] "normal Num. samples: 25"
## [1] "PD Num. samples: 25"
## [1] "SCC Num. samples: 24"
## [1] "AK Num. samples: 50"

     ## AK Control  Normal      PD     SCC 
     ## 47     133      24      25      23 
#pheno.ori<-cbind(pheno.ori,SCC,AK)
############################################################## 
############## SET up groups to be analysed

########### special case if just using AOGC for QC
AOGC<-pheno.ori$Project=="AOGC-NGS" & !is.na(pheno.ori$Project)
pheno.ori$SAMPLE[AOGC]
#########################################################

the.projects<-sample.types
table(pheno.ori$SampleProject)
ist<-1

pheno.names.ori<-colnames(pheno.ori)
for(ist in 1:length(sample.types)){

a.type<-rep(FALSE,times=dim(pheno.ori)[1])
a.type[pheno.ori$SampleProject %in% sample.types[ist]]<-TRUE
assign(sample.types[ist],value=a.type)
pheno.ori<-cbind(pheno.ori,a.type)

     }
colnames(pheno.ori)<-c(pheno.names.ori,sample.types)
colnames(pheno.ori)
########### special case if just using AOGC for QC
AOGC<-pheno.ori$Project=="AOGC-NGS" & !is.na(pheno.ori$Project)
pheno.ori$SAMPLE[AOGC]
pheno.ori<-cbind(pheno.ori,AOGC)
    sum(AOGC)


sample.types.full<-c(sample.types,"AOGC")



#################
table(pheno.ori$capture)
nextera.case<-pheno.ori$capture=="NxtD:NxtXR" & pheno.ori$AML
sum(nextera.case)
pheno.ori<-cbind(pheno.ori,nextera.case)
sample.types.full<-c(sample.types.full,"nextera.case")
#################

table(pheno.ori$capture)
nextera.control<-pheno.ori$capture=="NxtD:NxtXR" & pheno.ori$Control
sum(nextera.control)
pheno.ori<-cbind(pheno.ori,nextera.control)
sample.types.full<-c(sample.types.full,"nextera.control")


#################
table(pheno.ori$capture)
trueSeq.case<-pheno.ori$capture=="TruD:TruX" & pheno.ori$AML
sum(trueSeq.case)
pheno.ori<-cbind(pheno.ori,trueSeq.case)
sample.types.full<-c(sample.types.full,"trueSeq.case")
#################

table(pheno.ori$capture)
trueSeq.control<-pheno.ori$capture=="TruD:TruX" & pheno.ori$Control
sum(trueSeq.control)
pheno.ori<-cbind(pheno.ori,trueSeq.control)
sample.types.full<-c(sample.types.full,"trueSeq.control")


####################3

###############################################################################
table(pheno.ori$Aligner)
bwa.case<-(pheno.ori$Aligner=="bwa" | pheno.ori$Aligner=="novoalign;bwa" ) & pheno.ori$AML
sum(bwa.case)
pheno.ori<-cbind(pheno.ori,bwa.case)
sample.types.full<-c(sample.types.full,"bwa.case")


## table(pheno.ori$Aligner)
## bwa.control<-(pheno.ori$Aligner=="bwa" | pheno.ori$Aligner=="novoalign;bwa" ) & pheno.ori$Control
## sum(bwa.control) #0
## pheno.ori<-cbind(pheno.ori,bwa.case)
## sample.types.full<-c(sample.types,"bwa.case")
###############################################################################

###############################################################################
table(pheno.ori$Aligner)
novoalign.case<-pheno.ori$Aligner=="novoalign"  & pheno.ori$AML
sum(novoalign.case)
pheno.ori<-cbind(pheno.ori,novoalign.case)
sample.types.full<-c(sample.types.full,"novoalign.case")
###############################################################################

## table(pheno.ori$AffectionStatus)
## sum(is.na((pheno.ori$AffectionStatus)))

PD<-pheno.ori$AffectionStatus==2 & !pheno.ori$AML
pheno.ori$SAMPLE[PD]
pheno.ori$SampleProject[PD]
pheno.ori<-cbind(pheno.ori,PD)
    sum(PD)
sample.types.full<-c(sample.types.full,"PD")



#################
table(pheno.ori$capture)
nextera.PD<-pheno.ori$capture=="NxtD:NxtXR" & pheno.ori$PD
sum(nextera.PD)
pheno.ori<-cbind(pheno.ori,nextera.PD)
sample.types.full<-c(sample.types.full,"nextera.PD")
#################

#################
table(pheno.ori$capture)
trueSeq.PD<-pheno.ori$capture=="TruD:TruX" & pheno.ori$PD
sum(trueSeq.PD)
pheno.ori<-cbind(pheno.ori,trueSeq.PD)
sample.types.full<-c(sample.types.full,"trueSeq.PD")
#################

table(pheno.ori$Aligner)
bwa.PD<-(pheno.ori$Aligner=="bwa" | pheno.ori$Aligner=="novoalign;bwa" ) & pheno.ori$PD
sum(bwa.PD)
pheno.ori<-cbind(pheno.ori,bwa.PD)
sample.types.full<-c(sample.types.full,"bwa.PD")


###############################################################################
table(pheno.ori$Aligner)
novoalign.PD<-pheno.ori$Aligner=="novoalign"  & pheno.ori$PD
sum(novoalign.PD)
pheno.ori<-cbind(pheno.ori,novoalign.PD)
sample.types.full<-c(sample.types.full,"novoalign.PD")



####################3
#########################################################
##############################################################
##############################################################
##############################################################
sample.types.full
the.projects.ori<-the.projects


names(the.projects)<-the.projects
colnames(pheno.ori)

the.projects<-sample.types.full
names(the.projects)<-the.projects
the.projects

    

for (ir in 1: length(the.projects)){
  print(paste(the.projects[ir],"Num. samples:",sum(pheno.ori[,the.projects[ir]])))
      }
## [1] "Control Num. samples: 429"
## [1] "AML Num. samples: 131"
## [1] "AML-Child Num. samples: 17"
## [1] "AML-NotDiagnosis-Child Num. samples: 13"
## [1] "Asian-AML Num. samples: 15"
## [1] "Asian-AML-Child Num. samples: 7"
## [1] "Asian-AML-NotDiagnosis-Child Num. samples: 5"
## [1] "Asian-Control Num. samples: 25"
## [1] "AOGC Num. samples: 323"
## [1] "nextera.case Num. samples: 42"
## [1] "nextera.control Num. samples: 104"
## [1] "trueSeq.case Num. samples: 89"
## [1] "trueSeq.control Num. samples: 321"
## [1] "bwa.case Num. samples: 42"
## [1] "novoalign.case Num. samples: 89"
## [1] "PD Num. samples: 57"
## [1] "PD Num. samples: 57"
## [1] "nextera.PD Num. samples: 50"
## [1] "trueSeq.PD Num. samples: 7"
## [1] "bwa.PD Num. samples: 50"
## [1] "novoalign.PD Num. samples: 7"

 pheno.ori[1:5,]
     
# summary.geno.extra.ori <-summary.geno.extra

summary.geno.extra<-{}
####################################################################################
#################################################################################### REGULAR
#### MAY NEED TO ADJUST way use samples in selected based on selection below.
targets<-the.projects  #c("NMD","ex.Control","AOGC")
#targets<-sample.types # so include AOGC
targets



###### the.samples and pheno.ori in same order but the.samples has .GT extension.
it<-5

the.samples<-paste(pheno.ori[,"SAMPLE"],"GT",sep=".")
for(it in 1:length(targets)){
#use.samples<-the.samples[pheno.ori[,"SampleProject"]==targets[it]]
use.samples<-the.samples[pheno.ori[,targets[it]]]
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

#targets<-the.projects #c("NMD","ex.Control","AOGC")
targets
names(targets)<-paste(targets,".filt",sep="")
targets
it<-1
for(it in 1:length(targets)){
#use.samples<-the.samples[pheno.ori[,"SampleProject"]==targets[it]]
use.samples<-the.samples[pheno.ori[,targets[it]]]
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
summary.geno.extra[1:5,]
colnames(summary.geno.extra)[grepl("GENO",colnames(summary.geno.extra))]







genotypes<-a.indel[,paste(all.possible.samples,"GT",sep=".")] #### use a.indel.ori otherwise
a.indel.stats<-cbind(a.indel.stats,genotypes)

tail(colnames(a.indel.stats))
## a.indel<-a.indel.ori
## summary.geno.extra<-summary.geno.extra.ori


## a.indel.use<-a.indel
## summary.geno.extra.use <-summary.geno.extra

######################################################################################################
######################################################################################################
##########################################################################
##########################################################################
########################### START GENOTYPE REVOVERY  MODEL ############################
##########################################################################
######################################################################################################
######################################################################################################

use.genotype.recovery<-FALSE
force.recovery.model<-FALSE
if(use.genotype.recovery & (!exists("a.indel.ori") | force.recovery.model)){ ### don't do revovery is alrready done!!
## cellularity<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Sequenza/NSLM.cellularity.summary.use.csv",header=T,sep="\t",fill=TRUE,skip=0,stringsAsFactors=FALSE)
a.indel.ori<-a.indel
summary.geno.extra.ori<-summary.geno.extra


## alt.counts.thresh.ori<-alt.counts.thresh
## rare.in.Control.ori<-rare.in.Control
## rare.in.Control.filt.ori<-rare.in.Control.filt
## has.one.geno.ori<-has.one.geno
## rare.in.group.ori<-rare.in.group
## no.genotypes.filt.ori<-no.genotypes.filt
## no.genotypes.ori<-no.genotypes
## ok.missing.ori<-ok.missing


######################################################################################################
######################################################################################################
##################### Only revover rare genotypes likely to be somatic
n<-max(as.integer(summary.geno.extra.ori[,"TOTAL.Alleles.Control"]))

#p<-0.001
p<-0.01 ########### set MAF threshols HEREX1
sd.thresh<-6
n
p

alt.counts.thresh<-1
while( (alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) <= sd.thresh){alt.counts.thresh<-alt.counts.thresh+1}
alt.counts.thresh


## alt.counts.thresh<-1
## while( (alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) <= sd.thresh){alt.counts.thresh<-alt.counts.thresh+1}
## alt.counts.thresh

n<-as.numeric(summary.geno.extra.ori[,"TOTAL.Alleles.Control"])
alt.counts.thresh  <- round((sd.thresh*sqrt(n*p*(1-p)) + n*p))
alt.counts.thresh.ori<-alt.counts.thresh


alt.counts.thresh[1:50]
summary.geno.extra.ori[1:5,grepl("GENO",colnames(summary.geno.extra))]
## rare.in.Control<-as.numeric(summary.geno.extra[,"ALT.Alleles.Control"])< alt.counts.thresh
## rare.in.Control.filt <-as.numeric(summary.geno.extra[,"ALT.Alleles.Control.filt"])< alt.counts.thresh

rare.in.Control<-as.numeric(summary.geno.extra.ori[,"ALT.Alleles.Control"])<= alt.counts.thresh
rare.in.Control.filt <-as.numeric(summary.geno.extra.ori[,"ALT.Alleles.Control.filt"])<= alt.counts.thresh



sum(rare.in.Control)
sum(rare.in.Control.filt )
names(rare.in.Control)<-key
names(rare.in.Control.filt )<-key


length(maf.filter)
length(rare.in.Control)

maf.lt.all[1:5,]
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


## test<-"chr5:68667282:68667282:T:G:snp:68667282"
## summary.geno.extra[test,"ALT.Alleles.Control"]
## rare.in.Control[test]
######################################################################################################




## a.indel<-a.indel.ori
## summary.geno.extra<-summary.geno.extra.ori


#colnames(pheno.ori)

#pheno[pheno[,"AML"],c("SAMPLE","Blast.Count..", "Control","AML","SampleProject", "AffectionStatus"  )]
cellularity<- pheno.ori[,c("SAMPLE","Blast.Count..")]
cellularity[,"Blast.Count.."]<-cellularity[,"Blast.Count.."]/100
cellularity[is.na(cellularity[,"Blast.Count.."]) & pheno.ori[,"Control"] ,"Blast.Count.."]<-1 ## unnoen then make pure
# cellularity[pheno.ori[,"AML"],]
# cellularity[pheno.ori[,"Control"],]
# sample.types

PDs<-pheno.ori[pheno.ori[,"AML-Child"] | pheno.ori[,"Asian-AML-Child"] | pheno.ori[,"Asian-AML"]  | pheno.ori[,"AML-NotDiagnosis-Child"] | pheno.ori[, "Asian-AML-NotDiagnosis-Child"] | pheno.ori[,"Asian-Control"],"SAMPLE"]
## PDs<-c(PDs,"LPH-001-27_PD")   table(pheno.ori[pheno.ori[,"AffectionStatus"]==2,"SAMPLE"] )
PDs
#PDs.alt.counts<-alt.reads.reference.calls(a.indel,PDs,threshold=1)


cancer<-pheno.ori[pheno.ori[,"AML"],"SAMPLE"]
cancer
#cancer.alt.counts<-alt.reads.reference.calls(a.indel,cancer,threshold=1)
#cancer.alt.counts.true<-alt.reads.Non.reference.calls(a.indel,cancer,threshold=1)
Controls<-pheno.ori[pheno.ori[,"Control"],"SAMPLE"]

#indels,the.samples,AD.extension="AD",threshold=1,prefix="",suffix=""

a.indel.stats[1:5,1:20]
#Control.alt.counts<-alt.reads.reference.calls(a.indel.ori,Controls,AD.extension="AD",threshold=1)


Control.alt.counts<-alt.reads.reference.calls(a.indel.stats,Controls,AD.extension="FAD",threshold=1,prefix="",suffix="") ## needs a sample.GT column



## ##################################testing
## #Control.alt.counts.true<-alt.reads.Non.reference.calls(a.indel,Control,threshold=1)
## #Control.alt.counts.ori<-Control.alt.counts
## Control.alt.counts[1:50,]

## a.indel.ori[grepl("chr2:209113113:209113113",key),1:10]
## a.indel.ori[grepl("chr2:209113113:209113113",key),paste(cancer,".AD",sep="")]
## a.indel.ori[grepl("chr2:209113113:209113113",key),paste(cancer,".AD",sep="")]

## ## test<-c("chr12:11286221:11286221:T:A:snp","chr6:132029865:132029865:G:A:snp","chr7:82784807:82784807:C:A:snp","chr7:137206693:137206693:G:A:snp")

##  test<-c("chr15:90631934:90631934:C:T:snp","chr2:209113113:209113113:G:A:snp:209113113","chr15:90631838:90631838:C:T:snp") #,"chr12:11286221:11286221:T:A:snp")

##  test<-c("chr2:209113113:209113113:G:A:snp:209113113")
## test<-c("chr15:90631934:90631934:C:T:snp","chr15:90631838:90631838:C:T:snp")
## ## a.indel[test,paste(normals,".AD",sep="")]
## ## a.indel[test,paste(controls,".GT",sep="")]
## ## a.indel[test,paste(normals,".GT",sep="")]
## ## a.indel[test,paste(PDs,".AD",sep="")]
## ## a.indel[test,paste(PDs,".GT",sep="")]
## a.indel[test,paste(cancer,".AD",sep="")]
## a.indel[test,paste(cancer,".GT",sep="")]
## #a.indel.ori[test,paste(cancer,".GT",sep="")]

## #test<-c("chr15:40675107:40675107:C:T:snp") 
## #cbind(a.indel[test,paste(cancer,".GT",sep="")],a.indel.ori[test,paste(cancer,".GT",sep="")])
## Control.alt.counts[1:5,]
## geno.p1<-genotype.p.values(a.indel.ori[test,],c(PDs),AD.extension="AD",Control.alt.counts[test,"Read.Balance"]/100,cellularity)

## geno.p0<-genotype.p.values(a.indel.ori[test,],c(cancer),AD.extension="AD",Control.alt.counts[test,"Read.Balance"]/100,cellularity)
## geno.p1<-genotype.p.values(a.indel.stats[test,],c(cancer),AD.extension="FAD",Control.alt.counts[test,"Read.Balance"]/100,cellularity)
## #geno.p2<-genotype.p.values(a.indel[test,],c(PDs),normal.alt.counts[test,"Read.Balance"]/100,cellularity)


## 2*(pnorm(abs(c(4)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
## p.threshold=0.0026 # z=3:   2*(pnorm(abs(c(1:8)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
## p.threshold=2*(pnorm(abs(c(4)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
## p.threshold
## found.genotype<-  geno.p1<=p.threshold
## apply(found.genotype,1,sum)
## dim(geno.p1)
## apply(found.genotype,1,sum)


## rbind(a.indel.ori[test,paste0(c(PDs),".GT")][2,],
##       a.indel[test,paste0(c(PDs),".GT")][2,],
## a.indel[test,paste0(c(PDs),".AD")][2,],
##       signif(geno.p1[2,],2),
##       found.genotype[2,])


## geno.p[found.genotype]<-"0/1" # p.threshold<- 0.0026
## geno.p[!found.genotype]<-"0/0"
## colnames(geno.p)<-paste(colnames(geno.p),".GT",sep="")
## colnames(genotypes)
## geno.p[test,]
## genotypes[test,colnames(geno.p)]



## use.samples<-paste(c(cancer),".GT",sep="")
## genotypes<-a.indel[test,use.samples]
## dim(genotypes)
## summary.geno<-genotype.summary(as.matrix(genotypes))
## colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),"CANCER",sep=".")
## summary.geno



## ## use.samples<-paste(c(PDs),".GT",sep="")
## ## genotypes<-a.indel[test,use.samples]
## ## dim(genotypes)
## ## summary.geno<-genotype.summary(as.matrix(genotypes))
## ## colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),"PD",sep=".")
## ## summary.geno

## #summary.geno.extra[test,grepl("GENO",colnames(summary.geno.extra))]
## summary.geno.extra.ori[test,grepl("GENO",colnames(summary.geno.extra.ori))]


## ############## finish testing
## ##################################testing



# maf.filter
# rare.in.Control 
snp.only<-grepl("^snp",a.indel.ori[,"TYPE"]) ### the p.values codes uses "AD" so does not work for SNPs
is.flat<-grepl("flat$",a.indel.ori[,"TYPE"])
has.one.geno<-as.numeric(summary.geno.extra.ori[,"ALT.Alleles.AML"])>0 | as.numeric(summary.geno.extra.ori[,"ALT.Alleles.PD"])>0
alt.count.thresh<-1  # use  p[true.sample.minor.counts <= alt.count.thresh]<-1 default more than one needed
alt.count.thresh.include<-2 # ignore above if  more that this number of samples have true calls (so if have 2 samples can recover samples with one call)
has.one.geno.ori<-has.one.geno
to.recover<-snp.only & has.one.geno & !is.flat # & maf.filter & rare.in.Control  ## maf.filter & rare.in.Control probably not required
recover.samples<-c(cancer,PDs,Controls)
sum(to.recover)
AD.lower.tail <-FALSE





#geno.p<-genotype.p.values(a.indel.stats[to.recover ,],recover.samples,AD.extension="FAD",Control.alt.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh,AD.lower.tail) ## COULD USE NORMAL HERE
geno.p<-genotype.p.values.row(a.indel.stats[to.recover ,],recover.samples,AD.extension="FAD",Control.alt.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh,alt.count.thresh.include,AD.lower.tail) ## COULD USE NORMAL HERE
                           ## Control.alt.counts[1:5,]
                           ##  normal.alt.counts[1:5,]
#p.threshold=0.0026 # z=3:   2*(pnorm(abs(c(1:8)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))



## a.indel[test,)1:15]
## p.threshold.z.thresh<-2
## samples<-colnames(geno.p)[geno.p[test,] <=2*(pnorm(abs(c(p.threshold.z.thresh)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))]


## samples


## test<-c("chr2:209113113:209113113:G:A:snp:209113113") # IDN
## test<-c("chrX:123195620:123195620:G:T:snp")#  $STAK2
## test<-c("chr9:5073770:5073770:G:T:snp")


## ############################## check have for that mutation 
## ## posn<-grep(test,key)
## ## has.one.geno[posn]
## ## to.recover[posn]
## ## has.one.geno[posn]
## ## geno.p[test,]
## ## got<-rownames(a.indel.stats[to.recover ,])
## ## got[1:5]
## ## grep(test,got)

## samples<-c("AMAS-25.3-Diagnostic","AMAS-18.3-Diagnostic","87","AMLM12003H-J","AMLM12003H-K","AMLM12015WPS","AMLM12021D-M","AMLM12027NM","AMLM12028DAK","AMLM12040A-J","AMLM12PAH037A-B")

## number<-10 ## check top 10
## first20<-sort(geno.p[test,recover.samples])[1:number]
## samples<-names(first20)

## geno.p[test,samples]

## a.test<-cbind(1:number,geno.p[test,samples])
## colnames(a.test)<-c("Rank","P-value")
## a.test
## loc<-1
## trial<-cbind(a.indel[test,paste0(c(samples),".GT")],
##              a.indel.stats[test,paste0(c(samples),".GT")],
##       a.indel[test,paste0(c(samples),".GQ")],
## a.indel.stats[test,paste0(c(samples),".FAD")],signif(geno.p[test,samples],2))
## trial
## chk<-grep("^GENO",colnames(summary.geno))
## t(trial)

## range<-1:length(to.recover)
## recover.samples<-cancer
## a.indel.stats[to.recover ,paste(recover.samples,"GT",sep=".")][range,1:10] 
## a.indel.use[to.recover ,paste(recover.samples,"GT",sep=".")][range,1:10] =="0/0"
## geno.p[80:100,1:10]
## a.indel.stats[to.recover ,paste(recover.samples,"FAD",sep=".")][range,1:20] 
## a.indel.stats[to.recover ,paste(recover.samples,"DUP",sep=".")][range,1:10] 

## range.c<-1:length(recover.samples) #30:50
## a.indel.stats[to.recover ,paste(recover.samples,"GT",sep=".")][,range.c] 
## a.indel.use[to.recover ,paste(recover.samples,"GT",sep=".")][,range.c] =="0/0"
## geno.p[,range.c]
## a.indel.stats[to.recover ,paste(recover.samples,"FAD",sep=".")][,range.c] 
## a.indel.stats[to.recover ,paste(recover.samples,"DUP",sep=".")][,range.c]  

# dbinom(4,100,0.001) ## getting 4 alleles in 100 at 0.001 
p.threshold.z.thresh<-4
p.threshold=2*(pnorm(abs(c(p.threshold.z.thresh)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
p.threshold #  z=4 0.00006334248   z=6 0.000000001973175
found.genotype<-  geno.p <= p.threshold
geno.p[found.genotype]<-"0/1" # p.threshold<- 0.0026
geno.p[!found.genotype]<-"0/0"
colnames(geno.p)<-paste(colnames(geno.p),".GT",sep="")


# geno.p[test,paste(samples,".GT",sep="")] 
############# I don't want to correct GATK calls that are already 0/1 to 0/0:


genotypes<-a.indel.stats[to.recover ,paste(recover.samples,"GT",sep=".")] 
ref.call<-genotypes =="0/0"
geno.p[!ref.call]<-genotypes[!ref.call]


sum(!( colnames(geno.p) %in% colnames(a.indel))) # must ge zero
to.transfer<-colnames(geno.p)[colnames(geno.p) %in% colnames(a.indel)]

posns<-match(rownames(a.indel),rownames(geno.p))
missing<-is.na(posns)
sum(missing)
sum(!missing)
dim(geno.p)

a.indel[!missing,to.transfer]<-geno.p[posns[!missing],to.transfer]



## ########################### a.indel.stats contains non-recovered genotypes add the new ones now: would affest lib and aligner models slightly?
############################## problem is that bad aligner and lib position shoudl probbaly NOT be recovered so think this soudl remain as is
## ## genotypes<-a.indel[,paste(all.possible.samples,"GT",sep=".")] #### use a.indel.ori otherwise
## ## a.indel.stats<-cbind(a.indel.stats,genotypes)
## ## tail(colnames(a.indel.stats))

## posns<-match(rownames(a.indel.stats),rownames(geno.p))
## missing<-is.na(posns)
## sum(missing)
## sum(!missing)
## dim(geno.p)

## sum(!( colnames(geno.p) %in% colnames(a.indel.stats))) # must ge zero
## a.indel.stats[!missing,to.transfer]<-geno.p[posns[!missing],to.transfer]
## ########################################################




summary.geno.extra<-{}
####################################################################################
#################################################################################### REGULAR
#### MAY NEED TO ADJUST way use samples in selected based on selection below.
targets<-the.projects  #c("NMD","ex.Control","AOGC")
targets

###### the.samples and pheno in same order but the.samples has .GT extension.
it<-1

the.samples<-paste(pheno.ori[,"SAMPLE"],"GT",sep=".")
for(it in 1:length(targets)){
#use.samples<-the.samples[pheno.ori[,"SampleProject"]==targets[it]]
use.samples<-the.samples[pheno.ori[,targets[it]]]
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
#use.samples<-the.samples[pheno.ori[,"SampleProject"]==targets[it]]
use.samples<-the.samples[pheno.ori[,targets[it]]]
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



a.indel[1:5,1:10]
summary.geno.extra[1:5,]
colnames(summary.geno.extra)


rownames(summary.geno.extra)<-key

## ######################################### Check a few
## use.sample<-the.AK

## chk<-grep("^GENO",colnames(summary.geno))
## summary.geno["chr19:50169131:50169131:C:T:snp",chk]
## summary.geno.extra.ori<-summary.geno.extra
## a.snp<-"chr15:40675107:40675107:C:T:snp" # KNSTRN
## a.snp<-"chr19:50169131:50169131:C:T:snp" # BCL@L12
## a.snp<-"chr6:31940123:31940123:G:A:snp" #STK19

## chk<-grep("^GENO",colnames(summary.geno.extra))
## summary.geno.extra[a.snp,chk]

## ## chk<-grep("^GENO",colnames(summary.geno.extra.use))
## ## summary.geno.extra.use[a.snp,chk]

## chk<-grep("^GENO",colnames(summary.geno.extra.ori))
## summary.geno.extra.ori[a.snp,chk]


################################################# redo filtered data


while((dim(a.indel)[1] %% num.bits)< 2){num.bits<-num.bits+1} ### go don't get matrix issues
num.bits
#(dim(a.indel)[1] %% num.bits)
## fil.genotypes<-foreach(a.indel.bit=iter(a.indel,by='row',chunksize=as.integer(dim(a.indel)[1]/num.bits) ), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype(a.indel.bit,gsub(".GT$","",the.samples),prefix="",suffix="",20,0.02,0.98,0.20,0.80,7,2)

fil.genotypes<-foreach(a.indel.bit=iter(a.indel,by='row',chunksize=as.integer(dim(a.indel)[1]/num.bits) ), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype(a.indel.bit,gsub(".GT$","",the.samples),prefix="",suffix="",20,0.02,0.98,0.20,0.80,7,2)
# 20,0.02,0.98,0.20,0.8,10,5 # for cancers where het may be way out

# 0.02< het(0/1)<0.98
# ref(0/0)< 0.2
# (1/1) > 0.8


rownames(fil.genotypes)<-key
dim(fil.genotypes)
## colnames(fil.genotypes)[1:5]
## rownames(fil.genotypes)[1:5]

dim(fil.genotypes)
dim(a.indel)

tail(rownames(a.indel))




} #end genotype recovery

## a.indel.use<-a.indel
## summary.geno.extra.use <-summary.geno.extra

############################################################################################
############################ END reCovery #############################
############################################################################################
##########################################################################################



############################################################################################
############################ START BWA NOVOALIGN FIX #############################
############################################################################################
##########################################################################################
#SNP_STATS:version1
#FILTER_SUMMARY:alt_mismatch5_proportion_gt_0.8;Strand_bias_proportion_gt_0.9;Read_end3_proportion_gt_0.8;High_not_called_alts_gte_3
#SUMMARY_CALLED:ref_count;alt_count;mismatch_alt_count;fstrand_alt_count;rstrand_alt_count;read_end_alt_count;prop_mismatch;prop_fstrand;prop_rstrand;prop_read_end
#SUMMARY_NOT_CALLED:ref_count;alt_count;mismatch_alt_count;fstrand_alt_count;rstrand_alt_count;read_end_alt_count;high_alt_count
#SAMPLE.RF:GATK_called_snp;ref_count;alt_count;mismatch_alt_count;fstrand_alt_count;rstrand_alt_count;read_end_alt_count

## SDHA that align pooly with Nextera
## chr5:251680:251680:G:T:snp
## chr5:251672:251672:C:A:snp
## chr5:251674:251674:T:G:snp
## chr12:56676244:56676244:C:T:snp

## no support on ther strand for data
## chr5:251680:251680:G:T:snp
##
## test.samples<-c("AMLM12003H-J","AMLM12006CC","AMLM12007FF","AMLM12011J-G","AMLM12014DJG","AMLM12015WPS","AMLM12021D-M","AMLM12021M-R","AMLM12022N-A","AMLM12025GWR","AMLM12026MJC","AMLM12027NM","AMLM12028DAK","AMLM12030PGB","AMLM12031CDF","AMLM12036A-S","AMLM12036T-S","AMLM12038J-H","AMLM12040A-J","AMLM12PAH016K-B","AMLM12RMH026J-N")

## test.samples<-c("AMLM12001KP","AMLM12002K-B","AMLM12003H-J","AMLM12003H-K","AMLM12011J-G","AMLM12014N-R","AMLM12015WPS","AMLM12016D-F","AMLM12018AES","AMLM12018W-M","AMLM12019K-S","AMLM12019S-P","AMLM12020M-B","AMLM12021D-M","AMLM12021KS","AMLM12021M-R","AMLM12022N-A","AMLM12025GWR","AMLM12026MJC","AMLM12027MLM","AMLM12027NM","AMLM12028DAK","AMLM12029P-L","AMLM12030PGB","AMLM12031D-M","AMLM12032NRO","AMLM12036A-S","AMLM12036R-L","AMLM12036T-S","AMLM12037L-T","AMLM12037SAT","AMLM12038J-H","AMLM12039M-A","AMLM12040A-J","AMLM12PAH037A-B","AMLM12RMH026J-N")
## test.samples<-c("AMLM12001KP","AMLM12002K-B","AMLM12003H-J","100","35")
## test.samples<-c("AMLM12011J-G")

## test.loc<-c("chr5:251680:251680:G:T:snp","chr12:56676244:56676244:C:T:snp")
## get.samples<-expand.labels.to.samples(test.samples,c("FAD","TAD","DUP"),paste.after=TRUE)

## get.samples<-expand.labels.to.samples(test.samples,c("FAD"),paste.after=TRUE)
## get.samples<-expand.labels.to.samples(cancer,c("DUP"),paste.after=TRUE)

## rownames(filt)<-key
## filt[test.loc,1:12]

## grep(test.loc,key)

## a.indel.stats[test.loc,get.samples]            
## a.indel.stats[1:5,1:20]





#look for mutations mutation is germiline controls is not signifcatly different to those called in cancer



## table(pheno.ori[,"Aligner"])
## table(pheno.ori[pheno.ori[,"Aligner"]=="novoalign;bwa","SampleProject"])
## table(pheno.ori[pheno.ori[,"Aligner"]=="bwa","SampleProject"])
## table(pheno.ori[pheno.ori[,"Aligner"]=="novoalign","SampleProject"])

tapply(pheno.ori[,"Aligner"],pheno.ori[,"SampleProject"],length)

bwa<-pheno.ori[pheno.ori[,"Aligner"] %in% c("bwa","novoalign;bwa")  ,"SAMPLE"] ## Leo Pharma
novoalign<-pheno.ori[ !(pheno.ori[,"Aligner"] %in% c("bwa","novoalign;bwa")) & pheno.ori[,"AffectionStatus"]!=1  ,"SAMPLE"]
#normals<-normals[normals!="LPH-001-27_PD"]

bwa
novoalign
#normal.alt.counts<-alt.reads.reference.calls(a.indel,novoalign,threshold=1)

novo.alt.counts<-alt.reads.reference.calls(a.indel.stats,novoalign,AD.extension="FAD",threshold=1)
bwa.alt.counts<-alt.reads.reference.calls(a.indel.stats,bwa,AD.extension="FAD",threshold=1) ## needs a sample.GT column


##  normal.allele.freq<-Control.alt.counts[test,"Read.Balance"]/100
## tumor.sample.major.counts<-normal.alt.counts[test,"REF.reads"]
## tumor.sample.minor.counts<-normal.alt.counts[test,"ALT.reads"]
## test<-c("chr19:2939268:2939289:ACCACCCTTACCCAAGGAGGCA:-:indel","chr13:21729956:21729956:A:G:snp","chr2:198363406:198363406:C:T:snp","chr22:20760282:20760282:A:C:snp")
#test<-c("chr5:251680:251680:G:T:snp","chr12:56676244:56676244:C:T:snp")


## poissonian.model.position(Control.alt.counts[test,"Read.Balance"]/100,normal.alt.counts[test,"REF.reads"], normal.alt.counts[test,"ALT.reads"],av.tumor.contamination=0)
## poss.model.ori<-poss.model
## pass.possonian.aligner.model.ori<-pass.possonian.aligner.model

poss.model<-poissonian.model.position(bwa.alt.counts[,"Read.Balance"]/100,novo.alt.counts[,"REF.reads"], novo.alt.counts[,"ALT.reads"],1,lower.tail=FALSE) ## novo as reference negtive z more in bwa
## poss.model[test,]
poss.model[1:5,]
dim(poss.model)
dim(a.indel)
sum(rownames(bwa.alt.counts)!=rownames(novo.alt.counts))
sum(rownames(a.indel)!=rownames(poss.model))

#sum(poss.model[pass,"P-value"]<1e-5,na.rm=TRUE)
#poss.model[test,]  ### used in conjection with additional snp filtering
 # 2*(pnorm(abs(c(1:8)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)) #2*(pnorm(abs(seq(from=3,to=5,by=0.1)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
# 2*(pnorm(abs(4.4), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)) ##4 sd greater than mean
#signif(2*(pnorm(-6, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)),digits=3)


pass.possonian.aligner.model<- as.numeric(poss.model[,"Z-score"]) >= -8  | is.na( poss.model[,"P-value"])

#poss.model[!pass.possonian.aligner.model,][1:10,]
#test %in% rownames(poss.model)[!pass.possonian.aligner.model]

length(pass.possonian.aligner.model)
#length(pass)
sum(!pass.possonian.aligner.model) #60 with subset 2
#sum(!pass.possonian.aligner.model.ori)

#write.table(chk,file="contamination_BWA_Novo_Check.txt",sep="\t",row.names=F,col.names=T)


################################# library test
################################# library test
################################# library test
################################# library test
################################# library test
################################# library test
################################# library test
################################# library test

tapply(pheno.ori[,"capture"],pheno.ori[,"capture"],length)




## lib1<-unique(c("AMLM12001KP","AMLM12003H-J","AMLM12003H-K","AMLM12006CC","AMLM12007FF","AMLM12011J-G","AMLM12012SB","AMLM12014DJG","AMLM12014N-R","AMLM12015WPS","AMLM12018AES","AMLM12018W-M","AMLM12019S-P","AMLM12020M-B","AMLM12021KS","AMLM12021M-R","AMLM12025GWR","AMLM12026MJC","AMLM12027MLM","AMLM12030PGB","AMLM12031CDF","AMLM12032NRO","AMLM12036A-S","AMLM12036R-L","AMLM12037SAT","AMLM12038J-H","AMLM12039M-A","AMLM12003H-J","AMLM12003H-K","AMLM12014DJG","AMLM12014N-R","AMLM12014PJV","AMLM12015WPS","AMLM12016D-F","AMLM12018AES","AMLM12021D-M","AMLM12025GWR","AMLM12026MJC","AMLM12027NM","AMLM12028DAK","AMLM12030PGB","AMLM12031CDF","AMLM12036A-S","AMLM12038J-H","AMLM12040A-J","AMLM12PAH037A-B","AMLM12001KP","AMLM12003H-J","AMLM12003H-K","AMLM12007FF","AMLM12011J-G","AMLM12012SB","AMLM12015WPS","AMLM12020M-B","AMLM12021KS","AMLM12021M-R","AMLM12025GWR","AMLM12026MJC","AMLM12027MLM","AMLM12030PGB","AMLM12032NRO","AMLM12001KP","AMLM12003H-J","AMLM12011J-G","AMLM12020M-B","AMLM12021M-R","AMLM12028DAK","AMLM12030PGB","AMLM12031CDF","AMLM12032NRO","AMLM12037SAT","AMLM12038J-H","AMLM12040A-J","AMLM12PAH037A-B","AMLM12RMH026J-N","AMLM12001KP","AMLM12002K-B","AMLM12003H-K","AMLM12011J-G","AMLM12015WPS","AMLM12020M-B","AMLM12021D-M","AMLM12030PGB","AMLM12036R-L","AMLM12036T-S","AMLM12038J-H","AMLM12039M-A","AMLM12040A-J","AMAS-1.3-Diagnostic","AMAS-12.3-Diagnostic","AMAS-13.3-Diagnostic","AMAS-15.3-Diagnostic","AMAS-16.3-Diagnostic","AMAS-17.3-Relapse","AMAS-18.3-Diagnostic","AMAS-20.3-Diagnostic","AMAS-22.3-Diagnostic","AMAS-24.3-2ndAMLDiagnostic","AMAS-4.3-PreAML","AMAS-4.3-Relapse","AMAS-4.3-Remission","AMAS-5.3-Diagnostic","AMAS-5.3-Remission","AMAS-6.3-Diagnostic","AMAS-7.3-Diagnostic","AMAS-8.3-Diagnostic","AMAS-9.3-Diagnostic","AMLM12034H-F","AMLM12038BJD","AMAS-2.3-Remission","AMAS-24.3-2ndAMLDiagnostic","AMAS-3.3-Diagnostic","AMAS-3.3-Relapse","AMAS-3.3-Remission","AMAS-5.3-8PostAllograft","AMAS-5.3-Diagnostic","AMAS-9.3-Diagnostic","AMLM12005R-G","AMAS-1.3-Diagnostic","AMAS-12.3-Diagnostic","AMAS-13.3-Diagnostic","AMAS-15.3-Diagnostic","AMAS-16.3-Diagnostic","AMAS-17.3-Relapse","AMAS-5.3-Remission","AMAS-6.3-Diagnostic","AMAS-7.3-Diagnostic","AMAS-8.3-Diagnostic","AMLM12034H-F","AMAS-10.3-Diagnostic","AMAS-11.3-Diagnostic","AMAS-15.3-Diagnostic","AMAS-16.3-Diagnostic","AMAS-17.3-Relapse","AMAS-2.3-Remission","AMAS-23.3-Diagnostic","AMAS-4.3-Diagnostic","AMAS-4.3-Remission","AMAS-6.3-Diagnostic","AMAS-9.3-Diagnostic","AMAS-1.3-Diagnostic","AMAS-1.3-Remission","AMAS-11.3-Diagnostic","AMAS-16.3-Diagnostic","AMAS-18.3-Diagnostic","AMAS-2.3-Relapse","AMAS-2.3-Remission","AMAS-21.3-Diagnostic","AMAS-23.3-Diagnostic","AMAS-4.3-Diagnostic","AMAS-4.3-Relapse","AMAS-4.3-Remission","AMAS-5.3-8PostAllograft","AMAS-5.3-Diagnostic","AMAS-5.3-PostAllograft","AMAS-5.3-Remission","AMAS-7.3-Diagnostic","AMAS-8.3-Diagnostic"))



lib1<-pheno.ori[grepl("NxtXR",pheno.ori[,"capture"]) & pheno.ori[,"AffectionStatus"]!=1   ,"SAMPLE"] ## Leo Pharma
are.lib1<-pheno.ori[,"SAMPLE"] %in% lib1

lib2<-pheno.ori[ !(are.lib1) & pheno.ori[,"AffectionStatus"]!=1  ,"SAMPLE"]
are.lib2<-pheno.ori[,"SAMPLE"] %in% lib2


table(pheno.ori[are.lib1,"capture"])
table(pheno.ori[are.lib1,"SampleProject"])
table(pheno.ori[are.lib2,"capture"])


lib1.alt.counts<-alt.reads.reference.calls(a.indel.stats,lib1,AD.extension="FAD",threshold=1)
lib2.alt.counts<-alt.reads.reference.calls(a.indel.stats,lib2,AD.extension="FAD",threshold=1) ## needs a sample.GT column
#test<-c("chr7:44663912:44663912:A:T:snp","chrX:123195619:123195619:A:T:snp","chr7:44663913:44663913:C:T:snp","chrX:123200021:123200021:A:T:snp","chr5:251785:251785:G:A:snp")
## poss.model.lib.ori<-poss.model.lib
## pass.possonian.lib.model.ori<-pass.possonian.lib.model

    
poss.model.lib<-poissonian.model.position(lib2.alt.counts[,"Read.Balance"]/100,lib1.alt.counts[,"REF.reads"], lib1.alt.counts[,"ALT.reads"],1,lower.tail=FALSE) ## lib1 had not ref reads if significant


#poss.model<-poissonian.model.position(lib1.alt.counts[,"Read.Balance"]/100,lib2.alt.counts[,"REF.reads"], lib2.alt.counts[,"ALT.reads"],1,lower.tail=FALSE) ## lib2 as reference negtive z more in bwa

#poss.model.lib[test,]


dim(poss.model.lib)
dim(a.indel.stats)

pass.possonian.lib.model<- as.numeric(poss.model.lib[,"Z-score"]) <= 6  | is.na( poss.model.lib[,"P-value"])

#poss.model[!pass.possonian.aligner.model,][1:10,]
#test %in% rownames(poss.model)[!pass.possonian.aligner.model]
length(pass.possonian.lib.model)
#length(pass)
sum(!pass.possonian.lib.model) #739
#sum(!pass.possonian.lib.model.ori)








##########################################################################
##########################################################################
############################ END BWA NOVOALIGN FIX #############################
##########################################################################
##########################################################################




################################################################################################################
################################################################################################################
################################ VARIANT dependednt filters ######################################
################################################################################################################
################################################################################################################
################################################################################################################


## a.indel.use<-a.indel
## summary.geno.extra.use <-summary.geno.extra

## a.indel<-a.indel.ori
## summary.geno.extra<-summary.geno.extra.ori

## a.indel<-a.indel.use
## summary.geno.extra <-summary.geno.extra.use

###################################### GET HARDY Weinbery equilibrium of conrols
###################################### GET HARDY Weinbery equilibrium of conrols
###################################### GET HARDY Weinbery equilibrium of conrols
###################################### GET HARDY Weinbery equilibrium of conrols

colnames(summary.geno.extra)
#getHWE(obs_hets, obs_hom1, obs_hom2)
#hw.target<-"Control"  ## what to calculate HW with

hw.p.control<-getHWE(summary.geno.extra[,paste("GENO.",hw.target,sep="")]) ## used 16 CPUs
hw.p.control.filt<-getHWE(summary.geno.extra[,paste("GENO.",hw.target,".filt",sep="")]) ## used 16 CPUs


length(hw.p.control)
names(hw.p.control)<-key
names(hw.p.control.filt)<-key

hw.p.control[1:5]
hw.p.control.filt[1:5]

hw.controls.ok<-hw.p.control > hwe.control.threshold
hw.controls.ok.filt<-hw.p.control.filt > hwe.control.threshold

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

rare.in.group.table<-( summary.geno.extra[,paste("MAF.",c(targets.analysis,paste(targets.analysis,".filt",sep="")),sep="") ]< group.maf.thresh) | ( summary.geno.extra[,paste("MAF.",c(targets.analysis,paste(targets.analysis,".filt",sep="")),sep="")] > (1-group.maf.thresh)) | (is.na( summary.geno.extra[,paste("MAF.",c(targets.analysis,paste(targets.analysis,".filt",sep="")),sep="")]))
rare.in.group.table[1:5,]
# summary.geno.extra[1:5,]

## rare.in.group.table<-( summary.geno.extra[,paste("MAF.",c(the.projects,paste(targets,".filt",sep="")),sep="") ]< group.maf.thresh) | ( summary.geno.extra[,paste("MAF.",c(the.projects,paste(targets,".filt",sep="")),sep="")] > (1-group.maf.thresh)) | (is.na( summary.geno.extra[,paste("MAF.",c(the.projects,paste(targets,".filt",sep="")),sep="")]))
## rare.in.group.table[1:5,]
## summary.geno.extra[1:5,]

#rare.in.group.test<-colnames(rare.in.group)
rare.in.group.test<-paste("MAF.",targets.analysis,sep="")
rare.in.group.test
rare.in.group<-combine.boolean(rare.in.group.table,rare.in.group.test,"OR")
#sum(!rare.in.group)
#summary.geno.extra[!rare.in.group,paste("GENO.",targets.analysis,sep="")][1:5,]



rare.in.group[1:10]
## rare.in.group<-combine.boolean(rare.in.group.table,c("MAF.AML","MAF.Control"),"OR")
## rare.in.group.filt<-combine.boolean(rare.in.group.table,c("MAF.AML.filt","MAF.Control.filt"),"OR")
## sum(!rare.in.group)
## sum(!rare.in.group.filt)

no.genotypes.test<-paste("MAF.",targets.analysis,".filt",sep="")
no.genotypes.test
no.genotypes<-(summary.geno.extra[,no.genotypes.test]== 0)  | (is.na( summary.geno.extra[,no.genotypes.test])) # no genotypes in test classes for a mutataion after individaul quality filtering
no.genotypes[1:5,]
summary.geno.extra[1:5,no.genotypes.test]
no.genotypes.filt<-combine.boolean(no.genotypes,colnames(no.genotypes),"AND")
no.genotypes.filt[1:5]

no.genotypes.test<-paste("MAF.",targets.analysis,sep="")
no.genotypes.test
no.genotypes<-(summary.geno.extra[,no.genotypes.test]== 0)  | (is.na( summary.geno.extra[,no.genotypes.test])) # no genotypes in test classes for a mutataion after individaul quality filtering
no.genotypes[1:5,]
no.genotypes<-combine.boolean(no.genotypes,colnames(no.genotypes),"AND")
no.genotypes[1:5]

colnames(summary.geno.extra)[grepl("^MAF",colnames(summary.geno.extra))]

summary.geno.extra["chr7:150700484:150700484:G:A:snp",no.genotypes.test] # [,]

no.genotypes.extra.test<-paste("MAF.",cancer.group,sep="")
#no.genotypes.cancer<-(as.numeric(summary.geno.extra[,no.genotypes.cancer.test])== 0  & is.finite(as.numeric(summary.geno.extra[,no.genotypes.cancer.test]) ) )
no.genotypes.cancer<-(summary.geno.extra[,no.genotypes.extra.test]== 0)               
sum(no.genotypes.cancer)

no.genotypes.extra.test<-paste("MAF.","PD",sep="")
no.genotypes.PD<-(summary.geno.extra[,no.genotypes.extra.test]== 0)               
sum(no.genotypes.PD)


no.genotypes.extra.test<-paste("MAF.","PD",".filt",sep="")
no.genotypes.PD.filt<-(summary.geno.extra[,no.genotypes.extra.test]== 0)               
sum(no.genotypes.PD.filt)




sum(no.genotypes)
    length(no.genotypes)

#summary.geno.extra[1:5,]

print(paste0("High missing using: ",missing.targets))
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

   summary.geno.extra[1:5,c("GENO.AML.filt","GENO.Control","MISSING.Alleles.AML.filt","TOTAL.Alleles.AML.filt","MISSING.Alleles.Control.filt","TOTAL.Alleles.Control.filt")]
summary.geno.extra.ori[1:5,grep("AML$",colnames(summary.geno.extra.ori))]  
## high.missing<- cbind(as.numeric(summary.geno.extra[,"MISSING.Alleles.LOW"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.LOW"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.LOW"])),
##                      as.numeric(summary.geno.extra[,"MISSING.Alleles.HIGH"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.HIGH"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.HIGH"])),
##                      as.numeric( summary.geno.extra[,"MISSING.Alleles.LOW.pheno"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.LOW.pheno"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.LOW.pheno"])),
##                      as.numeric(summary.geno.extra[,"MISSING.Alleles.HIGH.pheno"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.HIGH.pheno"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.HIGH.pheno"])),
##                      as.numeric(summary.geno.extra[,"MISSING.Alleles.nimblegen"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.nimblegen"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.nimblegen"])),
##                      as.numeric(summary.geno.extra[,"MISSING.Alleles.illumina"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.illumina"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.illumina"]))
##                      )

colnames(high.missing.table)
high.missing.table[1:5,]

very.high.missing.controls.thresh<-0.9
very.high.missing.controls<- high.missing.table[,"Control.filt"]> very.high.missing.controls.thresh
sum(very.high.missing.controls)

high.total.missing<-subset(high.missing.table,select=missing.targets)
high.total.missing[1:5,]

high.total.missing<-high.total.missing <= missing.threshold
ok.missing<-combine.boolean(high.total.missing,missing.targets,"AND")
sum(ok.missing)

missing.targets
the.projects
colnames(high.missing.table)





################################################################################################################
################################################################################################################
################################ END VARIANT dependednt filters ######################################
################################################################################################################
################################################################################################################
################################################################################################################





################################################################################################################
################################################################################################################
################################ START POSITION dependednt filters ######################################
################################################################################################################
################################################################################################################
################################################################################################################



#sum(high.total.missing | nimblegen.total.missing | illumina.total.missing)
#high.missing[1:5,]
################## get loci in common

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


## test<-c("chr15:89848747:89848747:C:A:snp")
## test.pos<-grep(test,key)
## wanted.muts.NONcoding.keep[test.pos]
## wanted.interesting.to.prefilter[test.pos]
## wanted.interesting.to.prefilter.vep[test.pos]
## qual[test.pos,"GERP.low"]
## qual[test.pos,"GERP.unknown"]
## qual[test.pos,"GERP.high"]

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

## chk<-"chr21:36252969:36252969:A:G:snp"
## a.indel[chk,1:50]

## maf.lt.all[1:5,]
## maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.5"])

#pass<- rare.in.group & !no.genotypes & !high.missing & common.loci

##########################################################################
##########################################################################



#################################################################


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

#chk.in.repeat<-large.indel & !are.repeats
chk.in.repeat<- !are.repeats

are.sub.repeat<-indentify.IN.repeat(a.indel[chk.in.repeat,],looking="forward",bases.about=6,di.run.max=3,homo.run.max=5,genome="BSgenome.Hsapiens.UCSC.hg19")
remove.repeats<-key[chk.in.repeat][are.sub.repeat]
are.in.repeats.forward<- key %in% remove.repeats

remove.repeats[1:20]
sum(are.in.repeats.forward)
## [1] 6988

###################### in repeats looking back are.repeats[789]

sum(chk.in.repeat)
#chk.in.repeat<-large.indel & !are.repeats & !are.in.repeats.forward
chk.in.repeat<- !are.repeats & !are.in.repeats.forward

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



#hw.controls.ok[loci]

in.common.hit.gene <- a.indel[,"Gene.Names"] %in% common.hit.genes
in.common.hit.gene[1:5] 
sum(in.common.hit.gene)
the.chr
on.x.y<-a.indel[,"chr"] %in% c("X","Y","23","24","chrX","chrY")
sum(on.x.y)

#table(a.indel[,"TYPE"])
snp.only<-grepl("^snp",a.indel[,"TYPE"])



########################################################################## /media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/Gerp.test.r

## in.any.normal<-as.numeric(summary.geno.extra[,"ALT.Alleles.normal"]) > 0
## in.any.normal.filt<-as.numeric(summary.geno.extra[,"ALT.Alleles.normal.filt"]) > 0

## summary.geno.extra[!in.any.normal.filt,grepl("normal",colnames(summary.geno.extra))][1:5,]
## sum(!in.any.normal)

####################################### nextera bias
n.nextera.case<-max(as.integer(summary.geno.extra[,"ALT.Alleles.nextera.case"]))
n.case<-max(as.integer(summary.geno.extra[,"ALT.Alleles.AML"]))

p.nextera.default<-n.nextera.case/n.case

p.cohort<-as.integer(summary.geno.extra[,"TOTAL.Alleles.nextera.case"])/as.integer(summary.geno.extra[,"TOTAL.Alleles.AML"])
p.cohort[!is.finite(p.cohort)]<-p.nextera.default
names(p.cohort)<-rownames(summary.geno.extra)
test.nextera<-dbinom(as.integer(summary.geno.extra[,"ALT.Alleles.nextera.case"]), as.integer(summary.geno.extra[,"ALT.Alleles.AML"]), p.cohort, log=FALSE)
#test.nextera[as.integer(summary.geno.extra[,"ALT.Alleles.AML"])==0]<-1 ## case where ther are zero trials
names(test.nextera)<-rownames(summary.geno.extra)
nextera.bias<-(as.numeric(test.nextera) < 0.001)
#######################################################

####################################### trueSeq bias
n.trueSeq.case<-max(as.integer(summary.geno.extra[,"ALT.Alleles.trueSeq.case"]))
n.case<-max(as.integer(summary.geno.extra[,"ALT.Alleles.AML"]))

p.trueSeq.default<-n.trueSeq.case/n.case

p.cohort<-as.integer(summary.geno.extra[,"TOTAL.Alleles.trueSeq.case"])/as.integer(summary.geno.extra[,"TOTAL.Alleles.AML"])
p.cohort[!is.finite(p.cohort)]<-p.trueSeq.default
names(p.cohort)<-rownames(summary.geno.extra)
test.trueSeq<-dbinom(as.integer(summary.geno.extra[,"ALT.Alleles.trueSeq.case"]), as.integer(summary.geno.extra[,"ALT.Alleles.AML"]), p.cohort, log=FALSE)
#test.trueSeq[as.integer(summary.geno.extra[,"ALT.Alleles.AML"])==0]<-1 ## case where ther are zero trials
names(test.trueSeq)<-rownames(summary.geno.extra)
trueSeq.bias<-(as.numeric(test.trueSeq) < 0.001)
#######################################################


####################################### bwa bias
n.bwa.case<-max(as.integer(summary.geno.extra[,"ALT.Alleles.bwa.case"]))
n.case<-max(as.integer(summary.geno.extra[,"ALT.Alleles.AML"]))

p.bwa.default<-n.bwa.case/n.case

p.cohort<-as.integer(summary.geno.extra[,"TOTAL.Alleles.bwa.case"])/as.integer(summary.geno.extra[,"TOTAL.Alleles.AML"])
p.cohort[!is.finite(p.cohort)]<-p.bwa.default
names(p.cohort)<-rownames(summary.geno.extra)
test.bwa<-dbinom(as.integer(summary.geno.extra[,"ALT.Alleles.bwa.case"]), as.integer(summary.geno.extra[,"ALT.Alleles.AML"]), p.cohort, log=FALSE)
#test.bwa[as.integer(summary.geno.extra[,"ALT.Alleles.AML"])==0]<-1 ## case where ther are zero trials
names(test.bwa)<-rownames(summary.geno.extra)
bwa.bias<-(as.numeric(test.bwa) < 0.001)
#######################################################

####################################### novoalign bias
n.novoalign.case<-max(as.integer(summary.geno.extra[,"ALT.Alleles.novoalign.case"]))
n.case<-max(as.integer(summary.geno.extra[,"ALT.Alleles.AML"]))

p.novoalign.default<-n.novoalign.case/n.case

p.cohort<-as.integer(summary.geno.extra[,"TOTAL.Alleles.novoalign.case"])/as.integer(summary.geno.extra[,"TOTAL.Alleles.AML"])
p.cohort[!is.finite(p.cohort)]<-p.novoalign.default
names(p.cohort)<-rownames(summary.geno.extra)
test.novoalign<-dbinom(as.integer(summary.geno.extra[,"ALT.Alleles.novoalign.case"]), as.integer(summary.geno.extra[,"ALT.Alleles.AML"]), p.cohort, log=FALSE)
#test.novoalign[as.integer(summary.geno.extra[,"ALT.Alleles.AML"])==0]<-1 ## case where ther are zero trials
names(test.novoalign)<-rownames(summary.geno.extra)
novoalign.bias<-(as.numeric(test.novoalign) < 0.001)
#######################################################




############################################### repeat for the PD cohort
####################################### nextera bias
n.nextera.PD<-max(as.integer(summary.geno.extra[,"ALT.Alleles.nextera.PD"]))
n.PD<-max(as.integer(summary.geno.extra[,"ALT.Alleles.PD"]))

p.nextera.default<-n.nextera.PD/n.PD

p.cohort<-as.integer(summary.geno.extra[,"TOTAL.Alleles.nextera.PD"])/as.integer(summary.geno.extra[,"TOTAL.Alleles.PD"])
p.cohort[!is.finite(p.cohort)]<-p.nextera.default
names(p.cohort)<-rownames(summary.geno.extra)
test.PD.nextera<-dbinom(as.integer(summary.geno.extra[,"ALT.Alleles.nextera.PD"]), as.integer(summary.geno.extra[,"ALT.Alleles.PD"]), p.cohort, log=FALSE)
#test.PD.nextera[as.integer(summary.geno.extra[,"ALT.Alleles.PD"])==0]<-1 ## case where ther are zero trials
names(test.PD.nextera)<-rownames(summary.geno.extra)
nextera.PD.bias<-(as.numeric(test.PD.nextera) < 0.001)
#######################################################

####################################### trueSeq bias
n.trueSeq.PD<-max(as.integer(summary.geno.extra[,"ALT.Alleles.trueSeq.PD"]))
n.PD<-max(as.integer(summary.geno.extra[,"ALT.Alleles.PD"]))

p.trueSeq.default<-n.trueSeq.PD/n.PD

p.cohort<-as.integer(summary.geno.extra[,"TOTAL.Alleles.trueSeq.PD"])/as.integer(summary.geno.extra[,"TOTAL.Alleles.PD"])
p.cohort[!is.finite(p.cohort)]<-p.trueSeq.default
names(p.cohort)<-rownames(summary.geno.extra)
test.PD.trueSeq<-dbinom(as.integer(summary.geno.extra[,"ALT.Alleles.trueSeq.PD"]), as.integer(summary.geno.extra[,"ALT.Alleles.PD"]), p.cohort, log=FALSE)
#test.PD.trueSeq[as.integer(summary.geno.extra[,"ALT.Alleles.PD"])==0]<-1 ## case where ther are zero trials
names(test.PD.trueSeq)<-rownames(summary.geno.extra)
trueSeq.PD.bias<-(as.numeric(test.PD.trueSeq) < 0.001)
#######################################################


####################################### bwa bias
n.bwa.PD<-max(as.integer(summary.geno.extra[,"ALT.Alleles.bwa.PD"]))
n.PD<-max(as.integer(summary.geno.extra[,"ALT.Alleles.PD"]))

p.bwa.default<-n.bwa.PD/n.PD

p.cohort<-as.integer(summary.geno.extra[,"TOTAL.Alleles.bwa.PD"])/as.integer(summary.geno.extra[,"TOTAL.Alleles.PD"])
p.cohort[!is.finite(p.cohort)]<-p.bwa.default
names(p.cohort)<-rownames(summary.geno.extra)
test.PD.bwa<-dbinom(as.integer(summary.geno.extra[,"ALT.Alleles.bwa.PD"]), as.integer(summary.geno.extra[,"ALT.Alleles.PD"]), p.cohort, log=FALSE)
#test.PD.bwa[as.integer(summary.geno.extra[,"ALT.Alleles.PD"])==0]<-1 ## case where ther are zero trials
names(test.PD.bwa)<-rownames(summary.geno.extra)
bwa.PD.bias<-(as.numeric(test.PD.bwa) < 0.001)
#######################################################

####################################### novoalign bias
n.novoalign.PD<-max(as.integer(summary.geno.extra[,"ALT.Alleles.novoalign.PD"]))
n.PD<-max(as.integer(summary.geno.extra[,"ALT.Alleles.PD"]))

p.novoalign.default<-n.novoalign.PD/n.PD

p.cohort<-as.integer(summary.geno.extra[,"TOTAL.Alleles.novoalign.PD"])/as.integer(summary.geno.extra[,"TOTAL.Alleles.PD"])
p.cohort[!is.finite(p.cohort)]<-p.novoalign.default
names(p.cohort)<-rownames(summary.geno.extra)
test.PD.novoalign<-dbinom(as.integer(summary.geno.extra[,"ALT.Alleles.novoalign.PD"]), as.integer(summary.geno.extra[,"ALT.Alleles.PD"]), p.cohort, log=FALSE)
#test.PD.novoalign[as.integer(summary.geno.extra[,"ALT.Alleles.PD"])==0]<-1 ## case where ther are zero trials
names(test.PD.novoalign)<-rownames(summary.geno.extra)
novoalign.PD.bias<-(as.numeric(test.PD.novoalign) < 0.001)
#######################################################


## indels<-a.indel
use.samples.GQ<-paste(pheno.ori[,"SAMPLE"],"GT",sep=".")
use.samples.GQ<-use.samples.GQ[pheno.ori[,"AML"]]
the.QG.AML<-genotype.GQ.summary(a.indel,use.samples.GQ,"AML")

####
use.samples.GQ<-paste(pheno.ori[,"SAMPLE"],"GT",sep=".")
use.samples.GQ<-use.samples.GQ[pheno.ori[,"PD"]]
the.QG.PD<-genotype.GQ.summary(a.indel,use.samples.GQ,"PD")

use.samples.GQ<-paste(pheno.ori[,"SAMPLE"],"GT",sep=".")
use.samples.GQ<-use.samples.GQ[pheno.ori[,"Control"]]
the.QG.Controls<-genotype.GQ.summary(a.indel,use.samples.GQ,"Control")


the.QG.Controls[1:5,]
the.QG.AML[1:5,]
the.QG.PD[1:5,]
dim(the.QG.AML)
dim(the.QG.PD)
dim(a.indel)


GQ.AML.pass<-the.QG.AML[,"the.mean.AML"] > 50 | is.na(the.QG.AML[,"the.mean.AML"])
GQ.Control.pass<-the.QG.Controls[,"the.mean.Control"] > 50 | is.na(the.QG.Controls[,"the.mean.Control"])


##########################################################################
##########################################################################
##########################################################################
icc<-3
 if(target.pheno.col %in% case.control){
   for(icc in 1:length(case.control.classes)){
     recode<-  pheno.ori[,target.pheno.col] %in% names(case.control.classes)[icc]
     pheno.ori[recode,target.pheno.col]<-as.numeric(case.control.classes[icc])
   }}
  
pheno.ori[,target.pheno.col]<-as.numeric(pheno.ori[,target.pheno.col])
formula



## pheno.use<-pheno[pheno[,target.pheno] %in% c(0,1),]
## dim(pheno.use) pheno.ori.ori<-pheno.ori
if(!exists("pheno.ori") | !identical(colnames(pheno),colnames(pheno.ori))){
pheno.ori<-pheno
}

pheno<-pheno.ori[pheno.ori[, "SampleProject"] %in% c(0,1) ,] #SCC dim(pheno.ori)


    table(pheno.ori[, "Project"])
pheno<-pheno.ori[pheno.ori[, "SampleProject"] %in% c(0,1) & !(pheno.ori[, "Project"] %in% c("MODY","SDDS","SKDP"))  & !is.na((pheno.ori[, "Project"])),] #SCC dim(pheno.ori)



######################################### complex karyotype
## colnames(pheno.ori)
## table(sample.sheet.full[sample.sheet.full[,"SampleProject"]=="AML", c("Complex")])
## pheno<-pheno.ori[pheno.ori[, "SampleProject"] %in% c(1),]
## dim(pheno)

## (table(pheno$Complex))
##  ##   A  C  N 
##  ## 2 35 23 71

##  (table(pheno$SampleProject))

## ## C vs N or A
## complex<-pheno[,"Complex"] %in% c("C")
## not.complex<-pheno[,"Complex"] %in% c("A","N")
## exclude<-pheno[,"Complex"]==""
## 0   1 
## 106  23

## ## ## C or A v N
## complex<-pheno[,"Complex"] %in% c("C","A")
## not.complex<-pheno[,"Complex"] %in% c("N")
## exclude<-pheno[,"Complex"]==""
##  0  1 
## 71 58 


## ## ## C vs N
## complex<-pheno[,"Complex"] %in% c("C")
## not.complex<-pheno[,"Complex"] %in% c("N")
## exclude<-pheno[,"Complex"] %in% c("A") | pheno[,"Complex"]==""
##  0  1 
## 71 23

## ## ## C va A dry well
## ## ## complex<-pheno[,"Complex"] %in% c("C",)
## ## ## not.complex<-pheno[,"Complex"] %in% c("A")
## ## ## exclude<-pheno[,"Complex"] %in% c("N") | pheno[,"Complex"]==""

## ## ## A vs N
## complex<-pheno[,"Complex"] %in% c("A")
## not.complex<-pheno[,"Complex"] %in% c("N")
## exclude<-pheno[,"Complex"] %in% c("C") | pheno[,"Complex"]==""
## ##  0  1 
## ## 71 35 

## sum(complex)
## sum(not.complex)
## pheno[complex, "SampleProject"]<-1
## pheno[not.complex, "SampleProject"]<-0
## pheno<-pheno[!exclude,]
##  (table(pheno$SampleProject))


## ########catagorical
## ## C<-pheno[,"Complex"] %in% c("C")
## ## N<-pheno[,"Complex"] %in% c("N")
## ## A<-pheno[,"Complex"] %in% c("A")

## ## pheno[C, "SampleProject"]<-0
## ## pheno[A, "SampleProject"]<-1
## ## pheno[N, "SampleProject"]<-2

## ## exclude<-pheno[,"Complex"]==""

 #################### END COMPLEX KAROTYPE

#pheno<-pheno.ori
## pheno[!(pheno[, "Project"] %in% c("MODY","SDDS","SKDP"),"SAMPLE"]
##  pheno.ori(pheno.ori[, "Project"] %in% c("AOGC-NGS")) ,c("SAMPLE","SampleProject")]
## pheno[pheno[, "SampleProject"] %in% c(0),c("SAMPLE","SampleProject","Project")]

the.samples.use<-pheno[,"SAMPLE"]
the.samples.use<-paste(the.samples.use,".GT",sep="")






## heno[pheno[, "capture"]=="NA",]

 table(pheno[, "SampleProject"])
##   0   1 
## 323 135 
## table(pheno[, "capture"])
## KapD:NimX3           NxtD:NxtX            NxtD:NxtXR    TruD:NimX     TruD:NimX3;NxtD:NxtXR     TruD:TruX  Unknown
##     1                  0                     161            3                     1                183         6

table(pheno[, "Project"])
## AOGC-NGS RSGB_AML TGCM-AML 
##      323       45       90 
    



## ##                                  [,1]   [,2]      [,3]                             [,4] [,5]              
## ## chr7:44663912:44663912:A:T:snp   "107"  "OGDH"    "0.0000000000000237397842208046" NA   "Russell=26"      
## ## chr8:117868531:117868531:T:A:snp "965"  "RAD21"   "0.00000260458212690608"         NA   "Russell=9"       
## ## chrX:47003866:47003866:T:C:snp   "334"  "NDUFB11" "0.00000121806308157221"         NA   "Russell=11"      
## ## chrX:123184970:123184970:G:T:snp "226"  "STAG2"   "0.000000000128483683375683"     NA   "Russell=19"      
## ## chrX:123184971:123184971:C:T:snp "868"  "STAG2"   "0.00298023223876953"            NA   "Russell=5"       
## ## chrX:123195620:123195620:G:T:snp "1624" "STAG2"   "0.00152118575741536"            NA   "Russell=7 ;Tom=1"
## ## chrX:123204993:123204993:A:T:snp "1083" "STAG2"   "0.0195261498146793"             NA   "Russell=3 ;Tom=1"

## test<-"chr12:125396186:125396186:C:A:snp"
## max(test.nextera)
## test.nextera[test]
## as.integer(summary.geno.extra[test,"ALT.Alleles.trueSeq.case"])
## as.integer(summary.geno.extra[test,"ALT.Alleles.AML"])
## as.integer(summary.geno.extra[test,"ALT.Alleles.nextera.case"])
## as.integer(summary.geno.extra[test,"ALT.Alleles.AML"])
## p.cohort[test]

## colnames(pheno.ori)
## pheno.ori[pheno.ori$"trueSeq.PD","SAMPLE"] 

## as.integer(summary.geno.extra[test,"ALT.Alleles.trueSeq.PD"])
## as.integer(summary.geno.extra[test,"ALT.Alleles.PD"])
## as.integer(summary.geno.extra[test,"ALT.Alleles.nextera.PD"])
## as.integer(summary.geno.extra[test,"ALT.Alleles.PD"])
## p.cohort[test]
## dbinom(57,57,0.87,log=FALSE)


## pbinom(1,6,0.5, lower.tail = TRUE, log.p = FALSE)
## qbinom(0.0005,2,0.5, lower.tail = FALSE,log = FALSE)
## dbinom(2,2,0.3,log=FALSE)


############################ compare normal and 
## p=0.1
## n=10
## x<-seq(from=0,to=5,by=0.5)
## y<-x*sqrt(n*p*(1-p)) + n*p
## pbinom(y, n, p, lower.tail = FALSE, log.p = TRUE)
## pnorm(x, mean = 0, sd = 1, lower.tail = FALSE, log.p = TRUE)
## dbinom(as.integer(y), n, p, log = TRUE)
## ############################




############################## ADD WEIGHT MODEL HERE
############################## ADD WEIGHT MODEL HERE
############################## ADD WEIGHT MODEL HERE
############################## ADD WEIGHT MODEL HERE
## identical(snpinfo.ori,snpinfo.ori.keep)
## identical(snpinfo.ori,snpinfo.sliding)
######################################################################################################
n<-max(as.integer(summary.geno.extra[,"TOTAL.Alleles.Control"]))
#n<-max(as.integer(summary.geno.extra[,"TOTAL.Alleles.AOGC"]))

    
p<-0.001
# p<-0.01 ########### set MAF threshols HEREX1
# p<-0.005 ########### set MAF threshols HEREX1
# p<-0.05 ########### set MAF threshols HEREX1

do.MAFS<-c(0.001,0.005,0.01)
ido.mafs<-1
for (ido.mafs in 1:length(do.MAFS)){
    
######################################################################################################
n<-max(as.integer(summary.geno.extra[,"TOTAL.Alleles.Control"]))
#n<-max(as.integer(summary.geno.extra[,"TOTAL.Alleles.AOGC"]))

    
p<-do.MAFS[ido.mafs]

    
sd.thresh<-4
n
p

alt.counts.thresh<-1
while( (alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) <= sd.thresh){alt.counts.thresh<-alt.counts.thresh+1}
alt.counts.thresh

## alt.counts.thresh<-1
## while( (alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) <= sd.thresh){alt.counts.thresh<-alt.counts.thresh+1}
## alt.counts.thresh

n<-as.numeric(summary.geno.extra[,"TOTAL.Alleles.Control"])
# n<-as.numeric(summary.geno.extra[,"TOTAL.Alleles.AOGC"])
#  n  
alt.counts.thresh  <- round((sd.thresh*sqrt(n*p*(1-p)) + n*p))

alt.counts.thresh[1:50]
summary.geno.extra[1:5,]


rare.in.Control<-as.numeric(summary.geno.extra[,"ALT.Alleles.Control"])<= alt.counts.thresh
rare.in.Control.filt <-as.numeric(summary.geno.extra[,"ALT.Alleles.Control.filt"])<= alt.counts.thresh

## rare.in.Control<-as.numeric(summary.geno.extra[,"ALT.Alleles.AOGC"])<= alt.counts.thresh
## rare.in.Control.filt <-as.numeric(summary.geno.extra[,"ALT.Alleles.AOGC.filt"])<= alt.counts.thresh
    
alt.counts.thresh.4.rare.in.Controls<-alt.counts.thresh

    
sum(rare.in.Control)
sum(rare.in.Control.filt )
names(rare.in.Control)<-key
names(rare.in.Control.filt )<-key


################################################ AOGC controls ################
if(sum( c("TOTAL.Alleles.AOGC","ALT.Alleles.AOGC","ALT.Alleles.AOGC.filt") %in%  colnames(summary.geno.extra))==3  ){
n<-as.numeric(summary.geno.extra[,"TOTAL.Alleles.AOGC"])
#  n  
alt.counts.thresh  <- round((sd.thresh*sqrt(n*p*(1-p)) + n*p))

alt.counts.thresh[1:50]
summary.geno.extra[1:5,]


## rare.in.Control<-as.numeric(summary.geno.extra[,"ALT.Alleles.Control"])<= alt.counts.thresh
## rare.in.Control.filt <-as.numeric(summary.geno.extra[,"ALT.Alleles.Control.filt"])<= alt.counts.thresh

rare.in.AOGC<-as.numeric(summary.geno.extra[,"ALT.Alleles.AOGC"])<= alt.counts.thresh
rare.in.AOGC.filt <-as.numeric(summary.geno.extra[,"ALT.Alleles.AOGC.filt"])<= alt.counts.thresh    
#alt.counts.thresh.4.rare.in.Controls<-alt.counts.thresh

    
sum(rare.in.AOGC)
sum(rare.in.AOGC.filt )
names(rare.in.AOGC)<-key
names(rare.in.AOGC.filt )<-key
}
################################################################################

length(maf.filter)
length(rare.in.Control)

maf.lt.all[1:5,]
maf.lt.all[1:5,]
#maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.001"])
#maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.01"])
maf.col<-paste("MAF.lt",p,sep=":")
maf.col
 maf.filter<-as.logical(maf.lt.all[,maf.col])
## maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.5"])
sum(maf.filter)
names(maf.filter)<-key


colnames(pheno)
table(pheno$Aligner)
table(pheno$Capture.Method)





## test<-"chr5:68667282:68667282:T:G:snp:68667282"
## summary.geno.extra[test,"ALT.Alleles.Control"]
## rare.in.Control[test]
######################################################################################################






types<-c("coding","sliding.window","Single Point","bad.effect.with.Indels.noBenign")


itypes<-1

       
## for(itypes in 1:length(types)){

snap.file<-types[itypes]
snap.file

#snap.file<-paste(snap.file,".AOGC.",p,".FINAL.PCA.EXTRA_norecovery_GENOfilt",sep="")
#snap.file<-paste(snap.file,".AOGC.",p,".NEWwBIASwQGwMissing",sep="")
# snap.file<-paste(snap.file,".AOGC.",p,".KARYO_CorA.vs.N.CADD10.NEWwBIASwQGwmissing",sep="")
#snap.file<-paste("NEWwBIASwQG.",snap.file,".",p,".",project.files[ichr],sep="")

## snap.file<-paste(snap.file,".AOGC.",p,"CADD10.NEWwBIASwQGw",sep="")
snap.file
#snap.file<-paste(snap.file,".AOGC.",p,".KARYO_CvsAorN.CADD10.NEWwBIASwQG",sep="")
snap.file<-paste(snap.file,".AOGC.",p,".KARYO_AvsN.NEWwBIASwQG",sep="")
## paste(snap.file,".withGQ.filter",sep="")

## paste(snap.file,".noGQ.filter",sep="")
## a.indel.use<-a.indel
## summary.geno.extra.use<-summary.geno.extra

## a.indel<-a.indel.ori
## summary.geno.extra<-summary.geno.extra.ori

##  a.indel<-a.indel.use
## summary.geno.extra<-summary.geno.extra.use


if(types[itypes]=="Single Point"){
  
pass<- full.qual &  maf.filter  & not.flat.genotype  & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.control.filt & !in.any.normal.filt  & pass.possonian.control.model & !bad.qual.locations  #  90263

}

if(types[itypes]=="coding"){
  



## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene  & !on.x.y & !unannotated.hits & not.flat.genotype  & ok.missing & hw.controls.ok.filt & !no.genotypes.filt  & rare.in.Control.filt & rare.in.Control  & !are.repeats & !are.in.repeats  # previous

########## with no missing use don't use the .flt data  rare.in.Control.filt this one is location specific  FIX run ## no changes with repeats filters
## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control & pass.possonian.control.model & !bad.qual.locations   & pass.possonian.aligner.model

  
                                        # & !very.high.missing.controls # & !is.benign.missense #  # & ok.missing # & !are.repeats & !are.in.repeats #  & ok.missing.filt

## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene  & !no.genotypes   & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  # all dump

## cadd.missing<-is.na(cadd)
## sum(pass & cadd.missing) # 0 so all pass have a CADD (pass ontained using CADD
## cadd.10<- cadd>10 & !is.na(cadd)
## cbind(cadd,cadd.10)[1:100,]


    ##################### use for karytpe but not needed has no genotype location are not important
a.summary.geno<-genotype.summary(a.indel[,the.samples.use])
colnames(a.summary.geno)<-c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO")
test.no.geno<-as.numeric(a.summary.geno[,"MAF"])==0 # & is.finite(as.numeric(a.summary.geno[,"MAF"])
test.no.geno[is.na(test.no.geno)]<-TRUE        # test.no.geno[1802]                    
 

pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene & !no.genotypes  & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model & !no.genotypes.filt &  !novoalign.bias & !bwa.bias & !trueSeq.bias & !nextera.bias & GQ.Control.pass & GQ.AML.pass  #  & cadd.10 # &  ok.missing  #  & ok.missing  ## CADD and missing for SA
#######################

print("geno pass")
print(sum(pass & !test.no.geno))

#pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene & !no.genotypes  & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model & !no.genotypes.filt &  !novoalign.bias & !bwa.bias & !trueSeq.bias & !nextera.bias & GQ.Control.pass & GQ.AML.pass #  & cadd.10 # &  ok.missing  #  & ok.missing  ## CADD and missing for SA paper


## loci.1<-a.indel[,"Gene.Names"]=="TP53" & pass
## sum(loci.1)


## loci.2<-a.indel[,"Gene.Names"]=="TP53"  & pass
## sum(loci.2)
## loci<-loci.2 & !loci.1
## sum(loci)

## a.summary.geno[loci,]

## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene & !no.genotypes  & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model & !no.genotypes.filt & !novoalign.bias & !bwa.bias & !trueSeq.bias & !nextera.bias  & GQ.Control.pass & GQ.AML.pass ### one I'm using

## a.indel[,"chr"]=="chr10" &
##& ok.missing & !very.high.missing.controls  & rare.in.Control.filt 


pass.all.cohorts<- full.qual &  bad.coding  & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model & !no.genotypes.filt & !novoalign.bias & !bwa.bias & !trueSeq.bias & !nextera.bias #

## sum(pass.all.cohorts)
## sum(pass.all.cohorts & !pass )
## ## sum(pass & !rare.in.AOGC.filt) 
## ## sum(pass & !rare.in.Control.filt)

## ## #######final


 sum(pass.all.cohorts)

 sum(pass)
## sum(pass.all.cohorts | pass)

## ## pass.0.001.use<-pass #1091 # 984 814 #  
## pass.PD.cohorts.0.001.use<-pass.all.cohorts
## pass.all.cohorts.0.001.use<-pass.all.cohorts | pass
## ## maf.filter.0.001<- maf.filter
## ## rare.in.Control.0.001<-rare.in.Control
## ## rare.in.Control.filt.0.001<-rare.in.Control.filt
## ## alt.counts.thresh.4.rare.in.Controls.0.001<-alt.counts.thresh.4.rare.in.Controls

## pass.0.005.use<-pass #1360 # 1218 # 1087
## pass.all.cohorts.0.005.use<-pass.all.cohorts
## pass.PD.cohorts.0.005.use<-pass.all.cohorts
## pass.all.cohorts.0.005.use<-pass.all.cohorts | pass
## maf.filter.0.005<- maf.filter
## rare.in.Control.0.005<-rare.in.Control
## rare.in.Control.filt.0.005<-rare.in.Control.filt
## alt.counts.thresh.4.rare.in.Controls.0.005<-alt.counts.thresh.4.rare.in.Controls

## ## pass.0.01.use<-pass #1438 #1279 w bias  #1143 with no missing
## pass.PD.cohorts.0.01.use<-pass.all.cohorts
## pass.all.cohorts.0.01.use<-pass.all.cohorts | pass

## maf.filter.0.01<- maf.filter
## rare.in.Control.0.01<-rare.in.Control
## rare.in.Control.filt.0.01<-rare.in.Control.filt
## alt.counts.thresh.4.rare.in.Controls.0.01<-alt.counts.thresh.4.rare.in.Controls



 sum(pass)


## sum(pass.0.001.use & ! pass.0.01.use) #0


                                        # & !no.genotyp
#!very.high.missing.controls["chr2:212295788:212295788:C:T:snp:212295788"]
# pass<- full.qual &  bad.coding    & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt & !no.genotypes.filt  & pass.possonian.control.model & !bad.qual.locations   & pass.possonian.aligner.model & !are.repeats & !are.in.repeats ### filter for no maf
## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene  & !on.x.y & !unannotated.hits & not.flat.genotype  & ok.missing.filt & hw.controls.ok.filt & !no.genotypes.filt  & rare.in.controls.filt & rare.in.controls  & !are.repeats & !are.in.repeats # old AML


}

if(types[itypes]=="ALL.with.Indels.noBenign"){
  
## pass<- full.qual &  bad.effect & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control & pass.possonian.control.model & !bad.qual.locations   & pass.possonian.aligner.model

}




if(types[itypes]=="sliding.window"){
  
pass<- full.qual & maf.filter   & !in.common.hit.gene & !no.genotypes  & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model & !no.genotypes.filt

sum(pass) #6611 at 0.01

snpinfo.sliding<-sliding.window.generation(subset(a.indel,subset=pass)) ## use subset
snpinfo.sliding[1:5,]

using.sliding.window<-TRUE
if(!exists("snpinfo.ori.keep")){
snpinfo.ori.keep<-snpinfo.ori
}
snpinfo.ori<-snpinfo.sliding

}


##### full dump #  pass<-rep(TRUE,times=dim(a.indel)[1]) ; names(pass)<-key

# print(sum(pass))


## pass<-a.indel[,"Gene.Names"] %in% c("TP53","NOTCH1","NOTCH2","ATM","ACD","ASIP","BAP1","CASP8","CCND1","CDK4","MC1R","MITF","MTAP","MX2","OCA2","PARP1","PLA2G6","POT1","SLC45A2","TERF2IP","TERT","TYR","TYRP1","VDR","BCL2L12","KNSTRN","ISX","CDKN2A","BCL2L11","STK19","FJX1","TRHDE")
## sum(pass)


## good.genotypes<-c("chr9:5073770:5073770:G:T:snp","chr11:108224608:108224608:-:T:indel")
## summary.geno.extra[good.genotypes ,grepl("^GENO",colnames(summary.geno.extra))]
## ,"chr19:50169131:50169131:C:T:snp","chr19:50169104:50169104:C:T:snp","chr19:50169132:50169132:C:T:snp",  # BCL2L12
##                   "chr9:21971017:21971017:G:A:snp","chr9:21971016:21971016:G:A:snp","chr9:21971000:21971000:C:A:snp","chr9:21971056:21971056:C:T:snp","chr9:21971096:21971096:C:A:snp","chr9:21971099:21971099:G:A:snp","chr9:21971116:21971116:G:A:snp","chr9:21971015:21971015:C:A:snp","chr9:21994381:21994381:G:A:snp","chr9:21971040:21971040:C:T:snp","chr9:21994138:21994138:C:T:snp", #CDKN2A
##                   "chr7:1544063:1544063:G:A:snp","chr7:1544064:1544064:G:A:snp" # INTS1

#                  ) ## chr9 are CDKN2A low coverge exone 1 and 2

#pass[good.genotypes]


good.genotypes %in% names(pass)
good.genotypes %in% snp.fail.filt
#pass[ names(pass) %in% good.genotypes]<-TRUE





## bad.genotypes<-c("chr11:108121426:108121426:A:T:snp") # chr2:209113113:209113113
##  help[ rownames(help) %in% bad.genotypes]
## pass[ names(pass) %in% bad.genotypes]<-FALSE

############## add weights
#the.weights ## weights for gene with more than 5 syn mutations

## help<-cbind( pass,full.qual,bad.coding,bad.effect,maf.filter,in.common.hit.gene,unannotated.hits,not.flat.genotype,ok.missing,very.high.missing.controls,hw.controls.ok.filt,no.genotypes.filt,no.genotypes,rare.in.Control,rare.in.Control.filt,is.benign.missense,pass.possonian.control.model,are.repeats,are.in.repeats,pass.possonian.aligner.model,pass.possonian.lib.model ,bad.qual.locations,very.high.missing.controls )

    
help<-cbind( pass,full.qual,bad.coding,bad.effect,in.common.hit.gene,unannotated.hits,not.flat.genotype,ok.missing,very.high.missing.controls,hw.controls.ok.filt,no.genotypes,no.genotypes.filt,is.benign.missense,pass.possonian.control.model,are.repeats,are.in.repeats,pass.possonian.aligner.model,pass.possonian.lib.model,bad.qual.locations )
#####################

## snpinfo.ori[1:5,]
## help[good.genotypes,]


##################### check for missing data

## drop<-read.table("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/validated/dropped variants for Paul to Check.csv",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
## drop<-read.table("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/validated/dropped FANC variants for Paul to Check.csv",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

## drop[1:5,]

## bad.genotypes<-build.key(drop,core.ann,add.chr=TRUE)

## are.drops<-key %in% bad.genotypes
## sum(are.drops)
## dim(drop)

##  drop[!(bad.genotypes %in% key),]


## help[are.drops ,]

## a.indel[ pass,c(core.ann,"FILTER")]

## help[are.drops & ! pass,]

## a.indel[are.drops & ! pass,c(core.ann,"FILTER")]

## write.table(help[are.drops & ! pass,],"truth_table_missing_FANCS.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


## summary.geno.extra[are.drops & ! pass ,no.genotypes.test]
## summary.geno.extra[are.drops & ! pass ,grepl("^GENO",colnames(summary.geno.extra))]
## no.genotypes[are.drops & ! pass]
## no.genotypes.filt[are.drops & ! pass]

## samples<-drop[,"WES.ID"]
## samples<-unlist(strsplit(samples,split=", "))


## samples


## test<-c("chr5:68692376:68692376:-:AA:indel:68692375")
## test<c("chr1:161326566:161326566:A:T:snp")
## help[test,]
## summary.geno.extra[test,grepl("^GENO",colnames(summary.geno.extra))]
    


## test<-c("chr2:209113113:209113113:G:A:snp:209113113") # IDN
## test<-c("chrX:123195620:123195620:G:T:snp")#  $STAK2
## samples<-c("87","AMLM12003H-J","AMLM12003H-K","AMLM12015WPS","AMLM12021D-M","AMLM12027NM","AMLM12028DAK","AMLM12040A-J","AMLM12PAH037A-B")

## chk<-rownames(a.indel)[are.drops & ! pass]

## test<-"chr3:142215368:142215369:GG:-:indel" #   intronic
## a.indel[test,]
## help[test,]

## loc<-1
## trial<-rbind(a.indel[test,paste0(c(samples),".GT")],
##              a.indel.use[test,paste0(c(samples),".GT")],
##       a.indel[test,paste0(c(samples),".GQ")],
## a.indel.stats[test,paste0(c(samples),".FAD")],signif(geno.p[test,samples],2))
## trial
## chk<-grep("^GENO",colnames(summary.geno))
## trial

## range<-1:length(to.recover)
## recover.samples<-cancer
## a.indel.stats[to.recover ,paste(recover.samples,"GT",sep=".")][range,1:10] 
## a.indel.use[to.recover ,paste(recover.samples,"GT",sep=".")][range,1:10] =="0/0"
## geno.p[80:100,1:10]
## a.indel.stats[to.recover ,paste(recover.samples,"FAD",sep=".")][range,1:20] 
## a.indel.stats[to.recover ,paste(recover.samples,"DUP",sep=".")][range,1:10] 




## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene & !no.genotypes  & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model & !no.genotypes.filt



## bad.genotypes
## length(key)
## dim(a.indel)


## a.indel[are.drops & ! pass,c(core.ann,"FILTER")]

## key.short<-build.key(a.indel,c("chr","start"))
## key.short[1:5]

## bad.genotypes.key<-build.key(drop,c("chr","start"),add.chr=TRUE)
##  drop[!(bad.genotypes %in% key),]


## write.table(meta.results.burden,file=paste("Burden","ALL",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## #
                                        # ## >  drop[!(bad.genotypes %in% key),]
## ##    chr     start       end      REF ALT  TYPE   refGene..location refGene..type refGene..gene
## ## 8    7  99057799  99057799        T   C   snp   nonsynonymous SNV                      ATP5J2
## ## 28   3 142215369 142215370 GAAGTAAC   - indel frameshift deletion                         ATR
## ## 30   3 142215369 142215370 GAAGTAAC   - indel frameshift deletion                         ATR
## #chr11:108224608:108224608:-:T:indel # "VQSRTrancheINDEL99.00to99.90"

## #ATP5J2,ATP5J2-PTCD1

############################################## end check for missing

## length(key)
############################################
############################################
############################################
############### APPLY WEIGHTS ##################
############################################
############################################

############################################
############################################
############################################
############### END APPLY WEIGHTS ##################
############################################
############################################

#############################

sum(pass)
#  genotypes<-fil.genotypes[pass,the.samples.use] 
genotypes<-a.indel[pass,the.samples.use] ## ordered correctly for phenotypes
snp.names<-key[pass] ## GEFOS ony name with start

#### snpinfo now A different size than a.indel since added pathways!!!  snpinfo[snpinfo[,"gene"]=="KNSTRN",]
snpinfo<-snpinfo.ori[snpinfo.ori[,"Name"] %in% snp.names,]

if(!exists("gene.weights")){
  gene.weights.subset<-1
}else{
gene.weights.subset<-gene.weights[snpinfo.ori[,"Name"] %in% snp.names] # weight in same order as snpinfo.ori
}
snpinfo<-cbind(snpinfo,gene.weights.subset)
snpinfo[1:5,]
sum(is.na(as.numeric(snpinfo[,"gene.weights.subset"])))

###################################################/media/scratch/software/matlab/network.lic



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

meta.results.skat<-skatMeta(cohort.seq,SNPInfo = snpinfo,aggregateBy="cluster")
meta.results.skatO<-skatOMeta(cohort.seq,burden.wts =1,SNPInfo = snpinfo,aggregateBy="cluster",method = "saddlepoint") #method = "integration"
## print("start SKAT)")
## meta.results.skatO<-skatOMeta(cohort.seq,burden.wts ="gene.weights.subset",SNPInfo = snpinfo,aggregateBy="cluster",method = "integration")
if(types[itypes]=="sliding.window"){
  posns<-match(meta.results.burden[,"gene"],snpinfo[,"cluster"])
  missing<-is.na(posns)
  the.gene<-snpinfo[posns,"gene"]
  meta.results.burden<-cbind(the.gene,meta.results.burden)

    posns<-match(meta.results.skatO[,"gene"],snpinfo[,"cluster"])
  missing<-is.na(posns)
  the.gene<-snpinfo[posns,"gene"]
  meta.results.skatO<-cbind(the.gene,meta.results.skatO)
}


## meta.results.skatO<-{}

the.order<-     order(meta.results.burden[,"p"])
sum(is.na(meta.results.burden[,"p"])) ## bad p-values shoudl not happen
meta.results.burden<-  meta.results.burden[the.order,]
meta.results.burden[1:15,]
meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,] #meta.results.burden.3sd<-meta.results.burden
# meta.results.burden[meta.results.burden[,"gene"] %in% c("JAK2","OGDHL","NDUFB"),]




#meta.results.skatO<-{}
the.order<-     order(meta.results.skat[,"p"])
meta.results.skat<-  meta.results.skat[the.order,]
#meta.results.skat[1:50,]



the.order<-     order(meta.results.skatO[,"p"])
sum(is.na(meta.results.skatO[,"p"])) ## bad p-values shoudl not happen
meta.results.skatO<-  meta.results.skatO[the.order,]
meta.results.skatO[1:10,]
meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,]
## ## meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]
## clusters[1:5,]

## meta.results.burden[meta.results.burden[,"gene"] %in% clusters[,1],]
## meta.results.burden[meta.results.burden[,"gene"] %in% clusters[,2],]


## meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,]
                        

setwd(analysis.dir) #meta.results.skat<-{}
## getwd()
## bad.non.coding
## snap.file<-"coding.0.01.all.geno.all.filters"
## snap.file<-"coding.0.001.all.geno.all.filters_no.imput"
print(paste("Burden","ALL",snap.file,"txt",sep="."))
 write.table(meta.results.burden,file=paste("Burden","ALL",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
 write.table(meta.results.skatO,file=paste("SKATO","ALL",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
 write.table(meta.results.skat,file=paste("Skat","ALL",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

 write.table(meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,],file=paste("Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,],file=paste("SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skat[meta.results.skat[,"gene"] %in% clusters.wanted,],file=paste("Skat","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

## write.table(meta.results.skat[1:50,],file=paste("Skat",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(meta.results.skatO[1:50,],file=paste("SkatO",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(meta.results.skatO[1:50,],file=paste("SkatO",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,],file=paste("SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
##  write.table(pheno.ori,file="pheno.ori.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

 annotations<-a.indel[,c(1:6,16,28,7,30,34,37:42,43,14,32,33)]

## save(list=c("meta.results.skat","meta.results.skatO","meta.results.burden","pheno.use","snpinfo","genotypes","pass","high.missing","annotations","help","key","summary.geno.extra"),file=paste(paste(project.files[ichr],".small.RData",sep="")) )

} #MAFS

# test total snps      0.001/0.005/0.1
# AML vs Control       984/1218/1279
# AML vs Control CADD  814/977/1022
# C vs A or N CADD      323/405/437 #  Samples: 23(C) vs 106(AorN)
# C vs A or N           375/494/537      #  Samples: 23(C) vs 106(AorN)

# C or A vs N CADD     323/405/437 # Samples:  58(C or A) 71(N)
# C or A vs N          375/494/537      # Samples:  58(C or A) 71(N)

# C vs N CADD         235/297/322  # Samples:  23(C) 71(N) 
# C vs N              273/365/400  # Samples:  23(C) 71(N)  

# A vs N CADD         257/331/360 # Samples:  35(A) 71(N)


getwd()
save(list=c("synonymous","no.genotypes.cancer","cellularity","cancer","PDs","Control.alt.counts", "normal.alt.counts","the.samples.use","gene.weights","gene.weights.subset","filt","snp.fail.filt","use.wieght","weights","core.ann","case.control","snpinfo.ori","formula","clusters","pheno.types","ipheno","clusters.wanted","p","meta.results.skat","meta.results.skatO","meta.results.burden","pheno","pheno.ori","target.pheno.col","snpinfo","pass","high.missing.table","a.indel","a.indel.ori","help","key","summary.geno.extra","summary.geno.extra.ori","full.qual","bad.coding","bad.effect","maf.filter","in.common.hit.gene","on.x.y","unannotated.hits","not.flat.genotype","are.repeats","are.in.repeats","ok.missing","hw.controls.ok.filt","no.genotypes","rare.in.Control","rare.in.Control.filt","in.any.normal","in.any.normal.filt","are.in.repeats.back","are.in.repeats.forward","all.genes","contaminated","is.benign.missense","pass.possonian.control.model","pass.possonian.aligner.model","bad.qual.locations","filt"),file=paste(snap.file,".small_final.RData",sep="") )

} # itype

print(paste("Done: ",pheno.types[ipheno],"->",project.files[ichr]))
save.image(file=paste(project.files[ichr],".",pheno.types[ipheno],".",snap.file,".RData",sep=""))

} # pheno loop

} # loop over projects


} # loop over fam


 getwd()
## [1] "/media/old-scratch/media/scratch2/AOGC-NGS/Analysis"
## load("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/ALL_with_synon_Mar11_2015.RData")
## load("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-03_LeoPharma_NovoAlign/Analysis/coding.somatic.with.Indels.noBenign.withGQ.filter.recovery_weights_FULLQC_TIGHT.small_final.RData")
## colnames(a.indel)[1:50]

## getwd()

## setwd("/media/scratch/SCC_LOCAL")

## save.image("FINAL_with_wights_lm_TTN_GENO_recover_UNK_JUN_16_novo.RData")

save.image("SUBSET_HC3_FINAL.RData") ## a.indel IS recovered
save.image("SUBSET_AOGC_FINAL.RData")
save.image("FINAL_with_CADD_noRecovery.RData")
## save.image("FINAL_with_wights_6sd_rescure_JUN19.RData")
## save.image("FINAL_with_wights_6sd_resure_UNK_JUN.RData") 
## save.image("FINAL_with_wights_6sd_resure_UNK_27_APR.RData") ## @7th aprol has Maria weights  -15.36010     2.88198
## save.image("FINAL_with_wights_6sd_resure_UNK_14_MARR.RData")

###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
> meta.results.burden
  gene            p    beta        se  cmafTotal   cmafUsed nsnpsTotal nsnpsUsed nmiss
1 TP53 0.0001037268 2.61468 0.6735907 0.04263566 0.04263566          9         8     0
> meta.results.skatO
  gene             p         pmin rho       cmaf nmiss nsnps errflag
1 TP53 0.00006543338 0.0001117029   1 0.04263566     0     8       0
> meta.results.skat
  gene           p    Qmeta       cmaf nmiss nsnps
1 TP53 0.001286095 3088.434 0.04263566     0     8



> meta.results.burden
  gene            p    beta        se  cmafTotal   cmafUsed nsnpsTotal nsnpsUsed nmiss
1 TP53 0.0001037268 2.61468 0.6735907 0.04263566 0.04263566          8         8     0
> meta.results.skatO
  gene            p         pmin rho       cmaf nmiss nsnps errflag
1 TP53 0.0003460207 0.0001117029   1 0.04263566     0     8       0
> meta.results.skat
  gene           p    Qmeta       cmaf nmiss nsnps
1 TP53 0.001286095 3088.434 0.04263566     0     8


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

########
the.top<-1:dim(meta.results.burden)[1]
the.top<-1:500
to.unwind<-c(meta.results.burden[the.top,"gene"])# ,meta.results.skatO[1:the.top,"gene"])


clusters.wanted.subset<-c("FANC_complex.all", "Clinical" ,"BRCA.Genes","RAD51.Paralogues","BLM.Complex","Checkpoint.Proteins","citric","Citric_final","BLM.Complex_AND_Checkpoint","FANCD2_minimal_mono_ubi","MLH_cluster","Richard", "C1Alpha","C1Delta","C1Beta","C2","C3","C4","C5","MRC","MRC_IDH","C1","NDUFA","NDUFB","NDUFS","NDUFCV")
## to.unwind<-c("FGFR3") #, "MCM7", "RNPC3")
## interesting.gene<-c("BCL2L12","CCDC61","STK19","KNSTRN","TRHDE","FREM2","EBNA1BP2","PHACTR3","DCLK1","LRRIQ1","PHACTR3","CSMD3")
## ## to.unwind<-c("FANC_complex.all") # to.unwind<-meta.results.burden[8,"gene"]
##  to.unwind<-c("BCL2L12","CCDC61","STK19","KNSTRN","TRHDE","FREM2","EBNA1BP2","PHACTR3","DCLK1","LRRIQ1","PHACTR3","CSMD3") #"test.set"

## to.unwind<-c("TP53","NOTCH1","NOTCH2","FAT4","STK19","ISX","TRHDE","ARHGAP35","PREX1","KL","PIK3CA","KNSTRN","BCL2L12")
## to.unwind<-c("TP53","NOTCH1","NOTCH2","ATM","ACD","ASIP","BAP1","CASP8","CCND1","CDK4","MC1R","MITF","MTAP","MX2","OCA2","PARP1","PLA2G6","POT1","SLC45A2","TERF2IP","TERT","TYR","TYRP1","VDR","BCL2L12","KNSTRN","ISX","CDKN2A","BCL2L11","STK19","FJX1","TRHDE")

#to.unwind<-c("CDKN2A") #,"TCEA1","POLR2A","CTDP1")
## to.unwind<-c("Clinical")

## to.unwind<-c("FANCD2_minimal_mono_ubi")
## to.unwind<-c("BLM.Complex")
## to.unwind<-c("BLM.Complex_AND_Checkpoint")


to.unwind<-c(clusters.wanted.subset)
#to.unwind<-c("FANC_complex.all","Citric_final","citric","Clinical")


#sum(!(the.genes %in% wanted.genes))
## to.unwind<-c("Citric_final")
## to.unwind<-c("citric")

## to.unwind<-c("Ubin.proteo","lipid_raft","caveolae","Citric")
#to.unwind<-c(clusters.wanted[!(clusters.wanted %in% c("Ubin.proteo","lipid_raft","caveolae","Checkpoint_extendedx1","Checkpoint_extendedx2"))])
#grep(to.unwind,meta.results.burden[,"gene"])
#to.unwind
to.unwind.name<-to.unwind[1]
# to.unwind.name<-"EVERYTHING"
# to.unwind.name<-"TOP500"
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

sort(the.genes) #245 ### if used a cluster name need to do back up to (**) the.genes<-c(the.genes,"STAG2")

############repest to clean out cluster names 

the.genes.burden<-meta.results.burden[meta.results.burden[,"gene"] %in% the.genes,]

the.genes.burden
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
PDs<-pheno.ori[pheno.ori[,"AML-Child"] | pheno.ori[,"Asian-AML-Child"] | pheno.ori[,"Asian-AML"]  | pheno.ori[,"AML-NotDiagnosis-Child"] | pheno.ori[, "Asian-AML-NotDiagnosis-Child"],"SAMPLE"]
genotypes.PD<-a.indel[figure, c(paste(PDs,".GT",sep="")) ]
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
muts.in.cases<-apply(genotypes.ex[pheno[,"AML"],],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
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
depth.fad.cases<-rep("",times=length(loci))
dup.cases<-rep("",times=length(loci))

depth.controls<-rep("",times=length(loci))
depth.fad.controls<-rep("",times=length(loci))
dup.controls<-rep("",times=length(loci))

depth.PD<-rep("",times=length(loci))
depth.fad.PD<-rep("",times=length(loci))
dup.PD<-rep("",times=length(loci))

a.indel.sub<-a.indel[figure,]
a.indel.stats.sub<-a.indel.stats[figure,]


somatic.matrix.desc.full.sub<-somatic.matrix.desc.full[figure,]
somatic.matrix.p.full.sub<-somatic.matrix.p.full[figure,]

somatic.cases<-rep("",times=length(loci))
somatic.PD<-rep("",times=length(loci))

somatic.p.cases<-rep("",times=length(loci))
somatic.p.PD<-rep("",times=length(loci))



# a.indel.stats.sub[1:5,1:20]
# dim(a.indel.sub)

check<-1
for(check in 1:length(loci)){
# print(check)
#check<-"chr11:130066457:130066457:-:A:indel"
# posn<-grep(loci[check],key)
posn<-check


if(muts.in.PD[check]!=""){
#the.gt<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GT",sep=".")
the.gq<-paste(unlist(strsplit(muts.in.PD[check],split=",")),"GQ",sep=".")
quality.PD[check]<-paste(a.indel.sub[posn,the.gq],collapse=",")

the.ad<-paste(unlist(strsplit(muts.in.PD[check],split=",")),"AD",sep=".")
depth.PD[check]<-paste(a.indel.sub[posn,the.ad],collapse=";")

the.fad<-paste(unlist(strsplit(muts.in.PD[check],split=",")),"FAD",sep=".")
depth.fad.PD[check]<-paste(a.indel.stats.sub[posn,the.fad],collapse=";")

the.dup<-paste(unlist(strsplit(muts.in.PD[check],split=",")),"DUP",sep=".")
dup.PD[check]<-paste(a.indel.stats.sub[posn,the.dup],collapse=";")

the.ad.soma<-paste(unlist(strsplit(muts.in.PD[check],split=",")),"GT",sep=".")
somatic.PD[check]<-paste(somatic.matrix.desc.full.sub[posn,the.ad.soma],collapse=";")

the.ad.soma<-paste(unlist(strsplit(muts.in.PD[check],split=",")),"GT",sep=".")
somatic.p.PD[check]<-paste(signif(somatic.matrix.p.full.sub[posn,the.ad.soma],digits=4),collapse=";")


a.indel[posn,the.gq]
## a.indel[posn,the.gt]
## a.indel[posn,the.dp]
}




if(muts.in.cases[check]!=""){
#the.gt<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GT",sep=".")
the.gq<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GQ",sep=".")
quality.cases[check]<-paste(a.indel.sub[posn,the.gq],collapse=",")

the.ad<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"AD",sep=".")
depth.cases[check]<-paste(a.indel.sub[posn,the.ad],collapse=";")

the.fad<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"FAD",sep=".")
depth.fad.cases[check]<-paste(a.indel.stats.sub[posn,the.fad],collapse=";")

the.dup<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"DUP",sep=".")
dup.cases[check]<-paste(a.indel.stats.sub[posn,the.dup],collapse=";")


the.ad.soma<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GT",sep=".")
somatic.cases[check]<-paste(somatic.matrix.desc.full.sub[posn,the.ad.soma],collapse=";")

the.ad.soma<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GT",sep=".")
somatic.p.cases[check]<-paste(signif(somatic.matrix.p.full.sub[posn,the.ad.soma],digits=4),collapse=";")



a.indel[posn,the.gq]
## a.indel[posn,the.gt]
## a.indel[posn,the.dp]
}

if(muts.in.controls[check]!=""){
#the.gt<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"GT",sep=".")
the.gq<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"GQ",sep=".")
quality.controls[check]<-paste(a.indel.sub[posn,the.gq],collapse=",")

the.ad<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"AD",sep=".")
depth.controls[check]<-paste(a.indel.sub[posn,the.ad],collapse=";")

the.fad<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"FAD",sep=".")
depth.fad.controls[check]<-paste(a.indel.stats.sub[posn,the.fad],collapse=";")

the.dup<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"DUP",sep=".")
dup.controls[check]<-paste(a.indel.stats.sub[posn,the.dup],collapse=";")

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
muts.in.cases[1:10]
x<-muts.in.cases[6]

capture.counts.cases<-apply(as.matrix(muts.in.cases),1,function(x){
if(x!=""){
    x<-unlist(strsplit(x,split=","))
    x<-pheno.ori[pheno.ori[,"SAMPLE"] %in% x,"capture"]
    x<-table(x)
    x<-x[x!=0]
    x<-paste(paste(names(x),x,sep="="),collapse=" ;")
}else{x<-""}
x
}
)

Aligner.counts.cases<-apply(as.matrix(muts.in.cases),1,function(x){
if(x!=""){
    x<-unlist(strsplit(x,split=","))
    x<-pheno.ori[pheno.ori[,"SAMPLE"] %in% x,"Aligner"]
    x<-table(x)
    x<-x[x!=0]
    x<-paste(paste(names(x),x,sep="="),collapse=" ;")
}else{x<-""}
x
}
)


source.counts.cases<-apply(as.matrix(muts.in.cases),1,function(x){
if(x!=""){
    x<-unlist(strsplit(x,split=","))
    x<-pheno.ori[pheno.ori[,"SAMPLE"] %in% x,"sample.Source"]
    x<-table(x)
    x<-x[x!=0]
    x<-paste(paste(names(x),x,sep="="),collapse=" ;")
}else{x<-""}
x
}
)

########################
capture.counts.PD<-apply(as.matrix(muts.in.PD),1,function(x){
if(x!=""){
    x<-unlist(strsplit(x,split=","))
    x<-pheno.ori[pheno.ori[,"SAMPLE"] %in% x,"capture"]
    x<-table(x)
    x<-x[x!=0]
    x<-paste(paste(names(x),x,sep="="),collapse=" ;")
}else{x<-""}
x
}
)

Aligner.counts.PD<-apply(as.matrix(muts.in.PD),1,function(x){
if(x!=""){
    x<-unlist(strsplit(x,split=","))
    x<-pheno.ori[pheno.ori[,"SAMPLE"] %in% x,"Aligner"]
    x<-table(x)
    x<-x[x!=0]
    x<-paste(paste(names(x),x,sep="="),collapse=" ;")
}else{x<-""}
x
}
)


source.counts.PD<-apply(as.matrix(muts.in.PD),1,function(x){
if(x!=""){
    x<-unlist(strsplit(x,split=","))
    x<-pheno.ori[pheno.ori[,"SAMPLE"] %in% x,"sample.Source"]
    x<-table(x)
    x<-x[x!=0]
    x<-paste(paste(names(x),x,sep="="),collapse=" ;")
}else{x<-""}
x
}
)


extra.lib.info.cases<-cbind(capture.counts.cases,Aligner.counts.cases,source.counts.cases)
extra.lib.info.PD<-cbind(capture.counts.PD,Aligner.counts.PD,source.counts.PD)
################################ end counting in 



toString(colnames(a.indel)[c(1:6,8,11,16,28,7,30,34,35,36,37:42,43,14,32,33)])
colnames(a.indel)[1:60]

 ann.cols<-c("chr","start","end","REF","ALT","TYPE","refGene::type","knownGene::type","ensGene::type","Gene.Names","Genes.mentioned.at.ASH","refGene::location","knownGene::location","ensGene::location","OMIM (Gene::Status::OMIM::description::disease)","Consequence.Embl","Uploaded_variation.Embl","Gene.Embl","Feature.Embl", "Protein_position.Embl", "Amino_acids.Embl" , "ensGene::type","ID::maf","FILTER")# ,"rs.id")


annotations<-a.indel[,ann.cols]
dim(annotations)
dim(help)
dim(summary.geno.extra)
dim(a.indel)
dim(poss.model)
length(quality.cases)
length(figure)
dim(meta.results.burden.ex)
gerp.scores<-a.indel[,"gerp.scores"]
                             
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
if(!exists("summary.geno.extra.ori")){summary.geno.extra.ori<-summary.geno.extra}
if(!exists("summary.geno.extra.ori")){summary.geno.extra.ori<-summary.geno.extra}
## if(!exists("pass.old")){pass.old<-pass}
## if(!exists("pass.new")){pass.new<-pass}
#out<-cbind(meta.results.burden.ex,a.indel[figure,c(1:6,16,43,28,7,30,34,37:42)],summary.geno.extra[figure,c("GENO.AML","GENO.Control","GENO.AML.filt","GENO.Control.filt")],help[figure,],muts.in.cases,muts.in.controls)
a.functions<-a.indel[,c("PolyPhen.scores","SIFT.scores","PolyPhen.desc","SIFT.desc")]


posns<-match(key,filt.key)
missing<-is.na(posns)
sum(missing)
filt.sub<-filt[posns,]



meta.results.burden.ex[1:5,]
snpinfo.sliding[1:5,]
if(types[itypes]=="sliding.window"){ ### att the cluster name:
  
posns<-match(meta.results.burden.ex[,"gene"],snpinfo.sliding[,"Name"])
missing<-is.na(posns)
sum(missing)
the.window<-snpinfo.sliding[posns,"cluster"]
meta.results.burden.ex<-cbind(the.window,meta.results.burden.ex)

if(exists("snpinfo.ori.keep")){
snpinfo.ori<-snpinfo.ori.keep
rm(snpinfo.ori.keep)
}
using.sliding.window<-FALSE
}




all.GQ<-cbind(the.QG.AML,the.QG.PD,the.QG.Controls)
colnames(all.GQ)<-paste(colnames(all.GQ),"GQ",sep=".")

abundance<-cbind(test.nextera,test.trueSeq,test.bwa,test.novoalign,test.PD.nextera,test.PD.trueSeq,test.PD.bwa,test.PD.novoalign)
abundance[1:5,]
colnames(abundance)<-paste(colnames(abundance),"Pval",sep=".")
colnames(abundance)<-gsub("^test","EnRiched",colnames(abundance))

truth.table<-cbind(GQ.AML.pass,GQ.Control.pass,novoalign.bias,bwa.bias,trueSeq.bias,nextera.bias)


# filt.sub[figure,c("FILTER_SUMMARY","SUMMARY_CALLED","SUMMARY_NOT_CALLED")]
                             
## out<-cbind(meta.results.burden.ex,a.functions[figure,],gerp.scores[figure],annotations[figure,],maf.lt.all[figure,],is.benign.missense[figure],annotations[figure,],summary.geno.extra[figure,colnames(summary.geno.extra)[grep("^GENO",colnames(summary.geno.extra))]], filt.sub[figure,c("FILTER_SUMMARY","SUMMARY_CALLED","SUMMARY_NOT_CALLED")],pass.lose.filters[figure],pass.0.01.use[figure],pass.0.001.use[figure],pass.0.01.bad.loc[figure],pass.0.001.bad.loc[figure],pass.0.01.old[figure],pass.0.001.old[figure],help[figure,],high.missing.table[figure,],validated.posn[figure],poss.model[figure,],poss.model.lib[figure,],muts.in.cases,somatic.cases,somatic.p.cases,quality.cases,depth.fad.cases,depth.cases,dup.cases,muts.in.PD,somatic.PD,somatic.p.PD,quality.PD,depth.fad.PD,depth.PD,muts.in.controls,quality.controls,depth.fad.controls,depth.controls,dup.controls,summary.geno.extra.ori[figure,colnames(summary.geno.extra.ori)[grep("^GENO",colnames(summary.geno.extra.ori))]])

## out<-cbind(meta.results.burden.ex,a.functions[figure,],gerp.scores[figure],annotations[figure,],maf.lt.all[figure,],is.benign.missense[figure],annotations[figure,],summary.geno.extra[figure,colnames(summary.geno.extra)[grep("^GENO",colnames(summary.geno.extra))]], filt.sub[figure,c("FILTER_SUMMARY","SUMMARY_CALLED","SUMMARY_NOT_CALLED")],pass[figure],pass.0.01.use[figure],pass.0.005.use[figure],pass.0.001.use[figure],help[figure,],alt.counts.thresh.4.rare.in.Controls,high.missing.table[figure,],poss.model[figure,],poss.model.lib[figure,],muts.in.cases,somatic.cases,somatic.p.cases,quality.cases,depth.fad.cases,depth.cases,dup.cases,muts.in.PD,somatic.PD,somatic.p.PD,quality.PD,depth.fad.PD,depth.PD,muts.in.controls,quality.controls,depth.fad.controls,depth.controls,dup.controls,summary.geno.extra.ori[figure,colnames(summary.geno.extra.ori)[grep("^GENO",colnames(summary.geno.extra.ori))]]) ### use for AMP with one
enum<-1:dim(meta.results.burden.ex)[1]

out<-cbind(enum,meta.results.burden.ex,a.functions[figure,],gerp.scores[figure],annotations[figure,],maf.lt.all[figure,],is.benign.missense[figure],annotations[figure,],summary.geno.extra[figure,colnames(summary.geno.extra)[grep("^GENO",colnames(summary.geno.extra))]], filt.sub[figure,c("FILTER_SUMMARY","SUMMARY_CALLED","SUMMARY_NOT_CALLED")],pass.0.001.use[figure],pass.PD.cohorts.0.001.use[figure],pass.all.cohorts.0.001.use[figure],maf.filter.0.001[figure],rare.in.Control.0.001[figure],rare.in.Control.filt.0.001[figure],alt.counts.thresh.4.rare.in.Controls.0.001[figure],pass.0.005.use[figure],pass.PD.cohorts.0.005.use[figure],pass.all.cohorts.0.005.use[figure],maf.filter.0.005[figure],rare.in.Control.0.005[figure],rare.in.Control.filt.0.005[figure],alt.counts.thresh.4.rare.in.Controls.0.005[figure],pass.0.01.use[figure],pass.PD.cohorts.0.01.use[figure],pass.all.cohorts.0.01.use[figure],maf.filter.0.01[figure],rare.in.Control.0.01[figure],rare.in.Control.filt.0.01[figure],all.GQ[figure,],abundance[figure,],truth.table[figure,],alt.counts.thresh.4.rare.in.Controls.0.01[figure],help[figure,],high.missing.table[figure,],poss.model[figure,],poss.model.lib[figure,],validated.posn.FANC[figure],validated.posn.CITRIC[figure],muts.in.cases,somatic.cases,somatic.p.cases,quality.cases,depth.fad.cases,depth.cases,dup.cases,extra.lib.info.cases,muts.in.PD,somatic.PD,somatic.p.PD,quality.PD,depth.fad.PD,depth.PD,extra.lib.info.PD,muts.in.controls,quality.controls,depth.fad.controls,depth.controls,dup.controls,summary.geno.extra.ori[figure,colnames(summary.geno.extra.ori)[grep("^GENO",colnames(summary.geno.extra.ori))]]) ### use for out



#all.data[figure,]
#out<-cbind(meta.results.burden.ex,annotations[figure,],muts.in.cases,muts.in.controls)
dim(out)
#out[,1:13]
# help["chr7:150700484:150700484:G:A:snp",]


## table(out[,"refGene::location"])
## table(out[,"Consequence.Embl"]) # to.unwind.name<-"IDH"
getwd()
setwd(analysis.dir)
paste(paste(to.unwind,collapse="."))
paste(to.unwind.name,collapse=".")
paste(paste(to.unwind.name,collapse="."),p,"GENOTYPE.conponents.","SkatO","clusters",snap.file,"txt",sep=".")

order.by<-order(out[,"p"],decreasing=FALSE)
#enum<-1:dim(meta.results.burden.ex)[1]
out[order.by,][1:10,1:10]
setwd(analysis.dir)


write.table(out[order.by,],file=paste(paste(to.unwind.name,collapse="."),p,"GENOTYPE.conponents.","Burden","clusters_FINAL",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

getwd()

#"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/paper/Discovery_paper_final/EVERYTHING.GENOTYPE.conponents..Burden.clusters_FINAL.coding.somatic.with.Indels.AOGC.ALL.FINAL.PCA.txt"
file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/paper/Discovery_paper_final/EVERYTHING.GENOTYPE.FINAL.PCA.EXTRA.txt"
out<-read.table(file,header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

out[1:5,1:10]
rownames(out)<-out$gene

colnames(out)
# out<-out[order.by,]
## a.indel.use<-a.indel
## summary.geno.extra.use<-summary.geno.extra
## save.image("HC3.recovered_AML_CONTROL.RData")

#save(list=c("pheno.ori","pheno","a.indel","summary.geno.extra","use.wieght","pass.possonian.control.model","fil.genotypes"),file="No_revovery_model.RData")
out<-out.ori
test.logic<-(out[,"Gene.Names"] %in% genes.tested) & as.logical(out[,"pass.0.005.use"])  & !novoalign.bias & !bwa.bias & !trueSeq.bias & !nextera.bias  & GQ.pass  & as.logical(out[,"ok.missing"]) 


all.GQ<-cbind(the.QG.AML,the.QG.PD,the.QG.Controls)
colnames(all.GQ)<-paste(colnames(all.GQ),"GQ",sep=".")

abundance<-cbind(test.nextera,test.trueSeq,test.bwa,test.novoalign,test.PD.nextera,test.PD.trueSeq,test.PD.bwa,test.PD.novoalign)
abundance[1:5,]
colnames(abundance)<-paste(colnames(abundance),"Pval",sep=".")
colnames(abundance)<-gsub("^test","EnRiched",colnames(abundance))

the.QG.AML[1:5,]
the.QG.Controls[1:5,]

GQ.AML.pass<-the.QG.AML[,"the.mean.AML"] > 50 | is.na(the.QG.AML[,"the.mean.AML"])
GQ.Control.pass<-the.QG.Controls[,"the.mean.Control"] > 50 | is.na(the.QG.Controls[,"the.mean.Control"])

sum(GQ.Control.pass)
sum(GQ.AML.pass)
sum(GQ.AML.pass & GQ.Control.pass)

truth.table<-cbind(GQ.AML.pass,GQ.Control.pass,novoalign.bias,bwa.bias,trueSeq.bias,nextera.bias)
truth.table[1:5,]

out.new<-cbind(out[,1:104],all.GQ,abundance,truth.table,out[,105:159])
order.by<-order(out.new[,"p"],decreasing=FALSE)
#enum<-1:dim(meta.results.burden.ex)[1]
out[order.by,][1:10,1:10]
setwd(analysis.dir)


write.table(out.new,file="EVERYTHING_with_extra_filters.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

all.GQ[1:5,]
## hell<-cbind(novoalign.PD.bias,bwa.PD.bias,trueSeq.PD.bias,nextera.PD.bias,novoalign.bias,bwa.bias,trueSeq.bias,nextera.bias,the.QG.Controls,the.QG.AML,the.QG.PD)
## hell<-cbind(the.QG.AML,novoalign.PD.bias,bwa.PD.bias,trueSeq.PD.bias,nextera.PD.bias,novoalign.bias,bwa.bias,trueSeq.bias,nextera.bias)
hell<-cbind(test.logic,GQ.pass,the.QG.AML,the.QG.Controls,the.QG.PD,test.nextera,test.trueSeq,test.bwa,test.novoalign,test.PD.nextera,test.PD.trueSeq,test.PD.bwa,test.PD.novoalign,the.QG.Controls,the.QG.PD)




colnames(out)<-gsub("\\[figure\\]","",colnames(out))


file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/paper/AdelaidePassVariants_FANC_0.001.txt"
file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/paper/AdelaidePassVariants_MRC_0.005.txt"
file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/paper/pass.0.005.adel.csv"

filt<-read.table(file,header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
filt[1:5,]
dim(filt)

posns<-match(rownames(a.indel),filt[,"gene"])
missing<-is.na(posns)
sum(missing)
cadd<-filt[posns,"CADD"]
length(cadd)
cadd[1:100]
sum(is.na(cadd))
sum(!is.na(cadd))

dim(filt)
table( filt[,"AdelaidePASS"]) ## 584 for FANC
a.pass<-filt[as.logical(filt[,"AdelaidePASS"]),"UniqueID"]
length(a.pass)
dim(out)

pass.adel<-rownames(out) %in% a.pass
sum(pass.adel)

colnames(out)

genes.tested<-unique(as.character(out[pass.adel,"Gene.Names"]))
genes.tested
length(genes.tested)
length(unique(as.character(out[,"Gene.Names"]))) # 170/284
test<-"pass.all.cohorts.0.001.use"
test<-"pass.0.005.use"

sum( as.logical(out[,test ]) & pass.adel)
sum(  pass.adel & !as.logical(out[,test ]) )
the.QG.AML[1:5,]

GQ.pass<-the.QG.AML[,"the.mean.AML"] > 50 | is.na(the.QG.AML[,"the.mean.AML"])
sum(GQ.pass)


sum(rownames(out)!=names(GQ.pass))
## posns<-match(rownames(out),names(GQ.pass))

170/284 genes 
1310 vars at 0.005 in those genes originaly
3 removed because of capture or aligner bias
45 removed as mean GQ is < 50
152 removed by ok.missing  (64 re-induced by your selection)

651 selected
524 rejected


## missing<-is.na(posns)
## sum(missing)
## GQ.pass<-GQ.pass[posns]
sum( (out[,"Gene.Names"] %in% genes.tested) & as.logical(out[,"pass.0.005.use"]) )
test.logic<-(out[,"Gene.Names"] %in% genes.tested) & as.logical(out[,"pass.0.005.use"])  & !novoalign.bias & !bwa.bias & !trueSeq.bias & !nextera.bias  & GQ.pass  & as.logical(out[,"ok.missing"]) 
sum(test.logic)
sum(pass.adel)

sum(test.logic & pass.adel) # 649 of 1262 (1111 with ok.missing) 1310 ##45 removed 
sum(  pass.adel & !test.logic) 

discord<-pass.adel & !test.logic
sum(discord) # 64

discord<-test.logic & !pass.adel
concord<-test.logic & pass.adel
sum(concord)
sum(discord)

colnames(out)
out.sub<-out[discord,]



 
## gsub(", ",'\\",\\"',
## paste(colnames(out)[91:159],collapse=",")

##      "GENO.Control"                                    "GENO.AML"                                        "GENO.AML-Child"                                 
##  [73] "GENO.AML-NotDiagnosis-Child"                     "GENO.Asian-AML"                                  "GENO.Asian-AML-Child"                            "GENO.Asian-AML-NotDiagnosis-Child"              
##  [77] "GENO.Asian-Control"                              "GENO.AOGC"                                       "GENO.Control.filt"                               "GENO.AML.filt"                                  
##  [81] "GENO.AML-Child.filt"                             "GENO.AML-NotDiagnosis-Child.filt"                "GENO.Asian-AML.filt"                             "GENO.Asian-AML-Child.filt"                      
##  [85] "GENO.Asian-AML-NotDiagnosis-Child.filt"          "GENO.Asian-Control.filt"                         "GENO.AOGC.filt"

extra.cols<-c("pass.0.001.use","pass.all.cohorts.0.001.use","maf.filter.0.001","rare.in.Control.0.001","rare.in.Control.filt.0.001","alt.counts.thresh.4.rare.in.Controls.0.001","pass.0.005.use","pass.all.cohorts.0.005.use","maf.filter.0.005","rare.in.Control.0.005","rare.in.Control.filt.0.005","alt.counts.thresh.4.rare.in.Controls.0.005","pass.0.01.use","pass.all.cohorts.0.01.use","maf.filter.0.01","rare.in.Control.0.01","rare.in.Control.filt.0.01","alt.counts.thresh.4.rare.in.Controls.0.01","pass","full.qual","bad.coding","bad.effect","in.common.hit.gene","unannotated.hits","not.flat.genotype","ok.missing","very.high.missing.controls","hw.controls.ok.filt","no.genotypes","no.genotypes.filt","is.benign.missense","pass.possonian.control.model","are.repeats","are.in.repeats","pass.possonian.aligner.model","pass.possonian.lib.model","bad.qual.locations","AML.filt","Control.filt","Z-score","P-value","Z-score","P-value","validated.posn.FANC","validated.posn.CITRIC","muts.in.cases","somatic.cases","somatic.p.cases","quality.cases","depth.fad.cases","depth.cases","dup.cases","capture.counts.cases","Aligner.counts.cases","source.counts.cases","muts.in.PD","somatic.PD","somatic.p.PD","quality.PD","depth.fad.PD","depth.PD","capture.counts.PD","Aligner.counts.PD","source.counts.PD","muts.in.controls","quality.controls","depth.fad.controls","depth.controls","dup.controls")

sum(out.sub[,"ok.missing"]) ##64 that fail ok.missing  
cols<-c("p","Gene.Names","ID::maf","GENO.AML","GENO.AML.filt","GENO.AOGC","GENO.AOGC.filt","pass.0.005.use","ok.missing","pass.all.cohorts.0.001.use","muts.in.cases","quality.cases")

cols<-c("gene","p","Gene.Names","Consequence.Embl","gerp.scores","ID::maf","GENO.AML","GENO.AML.filt","GENO.AOGC","GENO.AOGC.filt", "GENO.Control", "GENO.Control.filt", "pass.0.005.use","ok.missing","pass.all.cohorts.0.001.use","muts.in.cases","quality.cases") 

cols %in% colnames(out)
extra.cols %in% colnames(out)



## hell<-cbind(novoalign.PD.bias,bwa.PD.bias,trueSeq.PD.bias,nextera.PD.bias,novoalign.bias,bwa.bias,trueSeq.bias,nextera.bias,the.QG.Controls,the.QG.AML,the.QG.PD)
## hell<-cbind(the.QG.AML,novoalign.PD.bias,bwa.PD.bias,trueSeq.PD.bias,nextera.PD.bias,novoalign.bias,bwa.bias,trueSeq.bias,nextera.bias)
hell<-cbind(test.logic,GQ.pass,the.QG.AML,test.nextera,test.trueSeq,test.bwa,test.novoalign,test.PD.nextera,test.PD.trueSeq,test.PD.bwa,test.PD.novoalign,the.QG.Controls,the.QG.PD)
sum(rownames(the.QG.AML) != names(test.nextera))
sum(rownames(hell) != rownames(out))

## posns<-match(rownames(out),rownames(hell))
## missing<-is.na(posns)
## sum(missing)
## hell<-hell[posns,]

hell[1:5,]

cbind(out[discord,cols],high.missing.table[discord,],hell[discord,])[1:20,]

hist(high.missing.table[discord,"AML.filt"],main="64 vars: missing rate in AML cohort (threhold=0.2)")
#hist(high.missing.table[concord,"AML.filt"],main="64 vars: missing rate in AML cohort (threhold=0.2)")
#savePlot("MRC.0.005.ok.missing.FALSE.png",type="png")

the.rejected<-cbind(pass.adel[discord],high.missing.table[discord,],hell[discord,],out[discord,cols[!(cols %in% extra.cols)]],out[discord,extra.cols])
the.concord<-cbind(pass.adel[concord],high.missing.table[concord,],hell[concord,],out[concord,cols[!(cols %in% extra.cols)]],out[concord,extra.cols])


table(summary.geno.extra[concord,"ALT.Alleles.AML"])
table(summary.geno.extra[discord,"ALT.Alleles.AML"])
chk<-summary.geno.extra[discord,"ALT.Alleles.AML"]==19
sum(chk)
the.rejected[chk,]

the.concord["chr2:25457162:25457162:A:G:snp",]
the.concord["chr15:90630398:90630398:C:T:snp",]

 rejected[1:3,]    
write.table(the.rejected,file="rejected.MRC.0.005.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(the.concord,file="concordant.MRC.0.005.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


filt.mrc<-filt
dim(filt.mrc)


filt<-rbind(filt,filt.mrc)
dim(filt)
filt[1:5,]

table(filt$AdelaidePASS)

length(unique(filt[,"gene"]))
sub<-out[as.logical(out[,"pass.0.01.use[figure]"]),]
dim(sub)
colnames(out)

table(out$TYPE)

1438

a.key<-build.key(out,core.ann,add.chr.label=FALSE)
posns<-match(filt[,"gene"],a.key)
missing<-is.na(posns) # 
sum(missing)
sub<-out[posns,]

cbind(filt,sub[c("gene","pass.0.01.use[figure]","CAD")])[1:10,]


write.table(a.indel[,core.ann],file=target,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE) 


############################################################33
test<-a.indel[,c("chr","start","REF","ALT")]
colnames(test)<-c("CHROM", "POS", "REF", "ALT")
test[1:5,]
test[,1]<-gsub("^chr","",test[,1])
write.table(test,file="test.vcf",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

write.table(a.indel[,core.ann],file=target,col.names="FALSE",row.names=FALSE,sep="\t",quote=FALSE)


target<-"update.annovar"
annovar.excutable.path<-"/media/UQCCG/Software/annovar"
genome.build<-"hg19"
  anno.DB.location.core<-"/media/UQCCG/Software/annovar/humandb"
anno.DB.location<-paste(anno.DB.location.core,genome.build,sep="/")
getwd()


function.filter.DB<-c("cadd")
names(function.filter.DB)<-c("CADD")


i<-1
for(i in 1:length(function.filter.DB)){

system(
    
    paste(annovar.excutable.path,"/","annotate_variation.pl -batchsize 100m -filter  -buildver ",genome.build," -dbtype ",function.filter.DB[i]," ",target,"  ",anno.DB.location," --outfile ",paste(target,function.filter.DB[i],sep=".")," ",sep="" )


      )
}
#annotate_variation.pl example/ex1.avinput humandb/ -filter -dbtype cadd -buildver hg19 -out ex1 

data.file<-paste(target,function.filter.DB[i],paste(genome.build,function.filter.DB[i],"dropped",sep="_"),sep=".")
the.data<-try(read.delim(data.file,header=F,skip=0,sep="\t",fill=TRUE,stringsAsFactors=FALSE,colClasses = "character"),silent=TRUE)
colnames(the.data)<-c("DB","ann",core.ann)
the.data[1:5,]


the.key<-build.key(the.data,core.ann,add.chr.label=FALSE)
a.key<-build.key(out,core.ann,add.chr.label=FALSE)

posns<-match(a.key,the.key)
missing<-is.na(posns) # 
sum(missing)
CAD<-the.data[posns,"ann"]
out<-cbind(out,CAD)

#######################################################################
dim(a.indel)
length(CAD)



update.annovar.cadd.hg19_cadd_dropped 


#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################


#GUess SOMATIC VS GERMLINE

Somatic in cancer chort
1) Is rescued from 0/0  ### only do IF have done genotype recovery
2) If have a 0/1 in controls and is difference in alleles 
3) check a few known somatic locations 


1):
cellularity<- pheno.ori[,c("SAMPLE","Blast.Count..")]
cellularity[,"Blast.Count.."]<-cellularity[,"Blast.Count.."]/100
cellularity[is.na(cellularity[,"Blast.Count.."]) & pheno.ori[,"Control"] ,"Blast.Count.."]<-1 ## unnoen then make pure
# cellularity[pheno.ori[,"AML"],]
# cellularity[pheno.ori[,"Control"],]
# sample.types

PDs<-pheno.ori[pheno.ori[,"AML-Child"] | pheno.ori[,"Asian-AML-Child"] | pheno.ori[,"Asian-AML"]  | pheno.ori[,"AML-NotDiagnosis-Child"] | pheno.ori[, "Asian-AML-NotDiagnosis-Child"] | pheno.ori[,"Asian-Control"],"SAMPLE"]
## PDs<-c(PDs,"LPH-001-27_PD")   table(pheno.ori[pheno.ori[,"AffectionStatus"]==2,"SAMPLE"] )
PDs
#PDs.alt.counts<-alt.reads.reference.calls(a.indel,PDs,threshold=1)


cancer<-pheno.ori[pheno.ori[,"AML"],"SAMPLE"]
cancer
Controls<-pheno.ori[pheno.ori[,"Control"],"SAMPLE"]

a.indel.stats[1:5,1:20]


Control.alt.counts<-alt.reads.reference.calls(a.indel.stats,Controls,AD.extension="FAD",threshold=1) ## needs a sample.GT column
Control.alt.counts[1:5,]
## snp.only<-grepl("^snp",a.indel.ori[,"TYPE"]) ### the p.values codes uses "AD" so does not work for SNPs
## is.flat<-grepl("flat$",a.indel.ori[,"TYPE"])
## has.one.geno<-as.numeric(summary.geno.extra.ori[,"ALT.Alleles.AML"])>0
## alt.count.thresh<-1  # use  p[true.sample.minor.counts <= alt.count.thresh]<-1 default more than one needed
## has.one.geno.ori<-has.one.geno
## to.recover<-snp.only & has.one.geno & !is.flat  & maf.filter & rare.in.Control  ## maf.filter & rare.in.Control probably not required
recover.samples<-c(cancer,PDs)
## sum(to.recover)
sum(pass)
to.recover<-pass


alt.count.thresh<-1  # use  p[true.sample.minor.counts <= alt.count.thresh]<-1 default more than one needed
alt.count.thresh.include<-2 # ignore above mor that this number of samples have true calls
AD.lower.tail <-FALSE


# geno.p<-genotype.p.values(a.indel.stats[to.recover ,] ,c(cancer,PDs),AD.extension="FAD",Control.alt.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh) ## COULD USE NORMAL HERE
# geno.p<-genotype.p.values(a.indel.stats[to.recover ,],recover.samples,AD.extension="FAD",Control.alt.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh) ## COULD USE NORMAL

geno.p<-genotype.p.values.row(a.indel.stats[to.recover ,],recover.samples,AD.extension="FAD",Control.alt.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh,alt.count.thresh.include,AD.lower.tail) ## COULD USE 
                                        #geno.p<-genotype.p.values(a.indel.stats[to.recover ,],recover.samples,AD.extension="FAD",Control.alt.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh,AD.lower.tail) ## COULD USE NORMAL HERE


p.threshold.z.thresh<-6
p.threshold=2*(pnorm(abs(c(p.threshold.z.thresh)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
p.threshold
found.genotype<-  geno.p <= p.threshold
geno.p[found.genotype]<-"0/1" # p.threshold<- 0.0026
geno.p[!found.genotype]<-"0/0"
colnames(geno.p)<-paste(colnames(geno.p),".GT",sep="")


############# I don't want to correct GATK calls that are already 0/1 to 0/0:


genotypes<-a.indel.stats[to.recover ,paste(recover.samples,"GT",sep=".")] 
ref.call<-genotypes =="0/0"
geno.p[!ref.call]<-genotypes[!ref.call]


found.genotype[!ref.call]<-FALSE ## no recovery of 0/1 or 1/1
dim(found.genotype)

#############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
2)
# summary.geno.extra.ori<-summary.geno.extra


genotypes<-a.indel.stats[to.recover ,paste(recover.samples,"GT",sep=".")] 
ref.call<-genotypes =="0/0"
geno.p[!ref.call]<-genotypes[!ref.call]


Control.alt.Non.REF.counts<-alt.reads.Non.reference.calls(a.indel.stats,Controls,AD.extension="FAD",threshold=1) ## needs a sample.GT column

Control.alt.Non.REF.counts[1:30,]
hist(Control.alt.Non.REF.counts[,"Read.Balance"])
dim(Control.alt.Non.REF.counts)

het50<-Control.alt.Non.REF.counts[,"Read.Balance"]==0
########################################### ASSUME IS 50%
Control.alt.Non.REF.counts[ Control.alt.Non.REF.counts[,"Read.Balance"]==0,"Read.Balance"]<-50.00
########################################### ASSUME IS 50%

has.an.alt <- as.numeric(summary.geno.extra.ori[to.recover , "ALT.Alleles.Control" ])>0

snp.only<-grepl("^snp",a.indel.ori[,"TYPE"]) ### the p.values codes uses "AD" so does not work for SNPs
is.flat<-grepl("flat$",a.indel.ori[,"TYPE"])
has.one.geno<-as.numeric(summary.geno.extra.ori[,"ALT.Alleles.AML"])>0
alt.count.thresh<-1  # use  p[true.sample.minor.counts <= alt.count.thresh]<-1 default more than one needed
has.one.geno.ori<-has.one.geno

to.recover<-snp.only & has.one.geno & !is.flat   ## maf.filter & rare.in.Control probably not required

########################################### ASSUME IS 50%
to.recover<-snp.only &  !is.flat   ## maf.filter & rare.in.Control probably not required
########################################### ASSUME IS 50%

recover.samples<-c(cancer,PDs)
sum(to.recover)



alt.count.thresh<-1  # use  p[true.sample.minor.counts <= alt.count.thresh]<-1 default more than one needed
alt.count.thresh.include<-2 # ignore above more that this number of samples have true calls
AD.lower.tail <-"BOTH"

#geno.p.Non.REF<-genotype.p.values(a.indel.stats[to.recover ,],recover.samples,AD.extension="FAD",Control.alt.Non.REF.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh,AD.lower.tail) ## COULD USE NORMAL
geno.p.Non.REF<-genotype.p.values.row(a.indel.stats[to.recover ,],recover.samples,AD.extension="FAD",Control.alt.Non.REF.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh,alt.count.thresh.include,AD.lower.tail) ##
#geno.p<-genotype.p.values(a.indel.stats[to.recover ,],recover.samples,AD.extension="FAD",Control.alt.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh,AD.lower.tail) ## COULD USE NORMAL HERE

geno.p.Non.REF[1:5,1:5]

p.threshold.z.thresh<-6
p.threshold=2*(pnorm((c(p.threshold.z.thresh)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
p.threshold
found.genotype.Non.REF<-  geno.p.Non.REF <= p.threshold

## sort(apply(found.genotype.Non.REF[has.an.alt,],1,sum),decreasing=TRUE)[1:20]
## sort(apply(found.genotype.Non.REF[!has.an.alt,],1,sum),decreasing=TRUE)[1:20]

genotypes<-a.indel.stats[to.recover ,paste(recover.samples,"GT",sep=".")] 
ref.call<-genotypes =="0/0"
#geno.p[!ref.call]<-genotypes[!ref.call]


geno.p.Non.REF[!has.an.alt,]<-1
found.genotype.Non.REF[!has.an.alt,]<-FALSE  #### don't do ones that were no 0/1 calles in controls

############# I don't want to correct GATK calls that are already 0/1 to 0/0:
found.genotype.Non.REF[ref.call]<-FALSE ## no recovery of 0/0
geno.p.Non.REF[ref.call]<-1


##  found.genotype.Non.REF[1:5,1:5]
##              found.genotype.Non.REF[test,paste0(c(samples),".GT")]
##          a.indel.stats[test,paste0(c(samples),".FAD")]

## summary.geno.extra[test,c("GENO.AML","GENO.Control")]
#################################################################
#somatic.matrix<-found.genotype | found.genotype.Non.REF
somatic.matrix<-found.genotype.Non.REF


somatic.matrix.desc<-matrix(data="GERM",nrow=dim(somatic.matrix)[1],ncol=dim(somatic.matrix)[2])
somatic.matrix.desc[somatic.matrix]<-"SOMA"

sort(apply(somatic.matrix.desc,1,function(x){sum(x=="SOMA")}),decreasing=TRUE)[1:20]

colnames(somatic.matrix.desc)<-paste(recover.samples,".GT",sep="")
rownames(somatic.matrix.desc)<-rownames(geno.p.Non.REF)

somatic.matrix.desc.full<-matrix(data="GERM",nrow=dim(a.indel)[1],ncol=length(all.possible.samples))
colnames(somatic.matrix.desc.full)<-paste(all.possible.samples,".GT",sep="")
rownames(somatic.matrix.desc.full)<-rownames(a.indel)

sum(!( colnames(somatic.matrix.desc) %in% colnames(a.indel))) # must ge zero
to.transfer<-colnames(somatic.matrix.desc)[colnames(somatic.matrix.desc) %in% colnames(somatic.matrix.desc.full)]

posns<-match(rownames(somatic.matrix.desc.full),rownames(somatic.matrix.desc))
missing<-is.na(posns)
sum(missing)
sum(!missing)
dim(somatic.matrix.desc.full)

somatic.matrix.desc.full[!missing,to.transfer]<-somatic.matrix.desc[posns[!missing],to.transfer]


#################################################################
colnames(geno.p.Non.REF)<-paste(colnames(geno.p.Non.REF),"GT",sep=".")
colnames(found.genotype.Non.REF)<-paste(colnames(found.genotype.Non.REF),"GT",sep=".")
found.genotype.Non.REF[1:5,1:5]

somatic.matrix.p.full<-matrix(data=1,nrow=dim(a.indel)[1],ncol=length(all.possible.samples))
colnames(somatic.matrix.p.full)<-paste(all.possible.samples,".GT",sep="")
rownames(somatic.matrix.p.full)<-rownames(a.indel)

to.transfer<-colnames(geno.p.Non.REF)[colnames(geno.p.Non.REF) %in% colnames(somatic.matrix.p.full)]

posns<-match(rownames(somatic.matrix.p.full),rownames(geno.p.Non.REF))
missing<-is.na(posns)
sum(missing)
sum(!missing)
dim(somatic.matrix.p.full)

somatic.matrix.p.full[!missing,to.transfer]<-geno.p.Non.REF[posns[!missing],to.transfer]


##########################################################
##########################################################
##########################################################


test<-c("chr7:44663912:44663912:A:T:snp")
samples<-c("AMLM12001KP","AMLM12003H-J","AMLM12003H-K","AMLM12006CC")
geno.p.Non.REF[1:5,1:5]
geno.p.Non.REF[test,paste0(c(samples),".GT")]
geno.p.Non.REF[test,]

trial<-rbind(a.indel[test,paste0(c(samples),".GT")],
             somatic.matrix.desc.full[test,paste0(c(samples),".GT")],
             signif( somatic.matrix.p.full[test,paste0(c(samples),".GT")],digits=4),
             found.genotype.Non.REF[test,paste0(c(samples),".GT")],
         a.indel.stats[test,paste0(c(samples),".FAD")]
             )
trial
summary.geno.extra[test,c("GENO.AML","GENO.Control")]
############use this in above

max(posns[!missing])
dim(somatic.matrix.desc)
##########################################################



chr15:90631934:90631934:C:T:snp
chr2:209113113:209113113:G:A:snp:209113113




p.threshold.z.thresh<-2
samples<-colnames(geno.p)[geno.p[test,] <=2*(pnorm(abs(c(p.threshold.z.thresh)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))]

samples<-c("21","31","66","68","81","89","AMLM12037SAT")
test<-c("chr14:45645955:45645955:A:C:snp") #,
test<-c("chr8:90994994:90994994:G:A:snp")# ,
test<-c("chr15:90631934:90631934:C:T:snp")
test<-c("chr2:209113113:209113113:G:A:snp:209113113")


summary.geno.extra[test,"GENO.Control"]
loc<-1

trial<-rbind(a.indel[test,paste0(c(samples),".GT")],
             a.indel.ori[test,paste0(c(samples),".GT")],
         a.indel.stats[test,paste0(c(samples),".FAD")],
                          found.genotype.Non.REF[test,samples],
               found.genotype[test,samples],
             signif(geno.p.Non.REF[test,samples],digits=4)
             )
trial


range<-1:length(to.recover)
recover.samples<-cancer
a.indel.stats[to.recover ,paste(recover.samples,"GT",sep=".")][range,1:10] 
a.indel.use[to.recover ,paste(recover.samples,"GT",sep=".")][range,1:10] =="0/0"
geno.p[80:100,1:10]
a.indel.stats[to.recover ,paste(recover.samples,"FAD",sep=".")][range,1:20] 
a.indel.stats[to.recover ,paste(recover.samples,"DUP",sep=".")][range,1:10] 











###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###############$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
######################## YOU can re start the caulstion here / modify pass etc from here
## keep<-as.logical(a.indel[,"wanted.muts"])
a.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/sanger_tested/FA-BRCA-HRR Network validated variants.csv"
a.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/sanger_tested/Missing variants 20150730.csv"



## extra filtering
#SNP_STATS:version1
#FILTER_SUMMARY:alt_mismatch5_proportion_gt_0.8;Strand_bias_proportion_gt_0.9;Read_end3_proportion_gt_0.8;High_not_called_alts_gte_3
#SUMMARY_CALLED:ref_count;alt_count;mismatch_alt_count;fstrand_alt_count;rstrand_alt_count;read_end_alt_count;prop_mismatch;prop_fstrand;prop_rstrand;prop_read_end
#SUMMARY_NOT_CALLED:ref_count;alt_count;mismatch_alt_count;fstrand_alt_count;rstrand_alt_count;read_end_alt_count;high_alt_count
#SAMPLE.RF:GATK_called_snp;ref_count;alt_count;mismatch_alt_count;fstrand_alt_count;rstrand_alt_count;read_end_alt_count

  print(paste0("Reading SNPstats file",SnpStats.file))
  filt<-read.table(SnpStats.file,header=T,skip=5,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
  filt<-filt[,1:11]

filt[1:5,]
filt.key<-build.key(filt,core.ann,add.chr.label=TRUE)


colnames(filt)
pass.filt<-strsplit(filt[,"FILTER_SUMMARY"],split=";")

pass.filt<-unlist(lapply(pass.filt,function(x) {sum(as.logical(x[1:3]))==0}))
#pass.filt[173713]
sum(!pass.filt)
#filt.key[173713]
snp.fail.filt<-filt.key[!pass.filt]
snp.fail.filt[1:5]


bad.qual.locations<-key %in% snp.fail.filt
sum(bad.qual.locations)
length(bad.qual.locations)


dim(filt)
length(pass.filt)
pa

posns<-match(key,filt.key)
missing<-is.na(posns)
sum(missing)
filt.sub<-filt[posns,]
setwd(analysis.dir)
pass.filt.sub<-pass.filt[posns]

pass.t<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control & pass.possonian.control.model   & pass.possonian.aligner.model # & !no.genotypes.filt  & !bad.qual.locations
sum(pass.t)

###############################################

miss.loc<-read.delim(a.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
key[1:5]

key.short<-build.key(a.indel,c("chr","start","REF","ALT"))
key.miss<-build.key(miss.loc,c("Chr","Position","Ref.Base","Mut.Base"),add.chr<-TRUE)

key.miss[1:5]
key.short[1:5]
are.missing<-!(key.miss %in% key.short)
sum(are.missing)
sum(!are.missing)



missing.keys.loc<-( key.short %in% key.miss )
validated.posn<-( key.short %in% key.miss ) 
help[validated.posn,]


validated.posn.FANC<-validated.posn
validated.posn.CITRIC<-validated.posn

sum(validated.posn.FANC)
sum(validated.posn.CITRIC)

validated.posn.FANC,validated.posn.CITRIC

key.miss<-c("chrX:123195619:123195619:A:T:snp","chr5:251785:251785:G:A:snp")
missing.keys.loc<-( key %in% key.miss )


help[missing.keys.loc,]



cbind(filt.sub[missing.keys.loc,],pass.filt.sub[missing.keys.loc],pass[missing.keys.loc])
pass.filt.sub[missing.keys.loc]


pass.filt[chk,]
miss.loc[are.missing,]



##   WES.ID Chr Position         Ref.Base Mut.Base   Gene Ref.Transcript.ID Type.of.Mutation Consequence.of.Mutation     DNA.Change Protein.Change Mutation.Status Allele.Load.... CADD.Score
## 1      0   3 10140453                G        A FANCD2         NM_033084              SNV                Missense      c.4235G>A       p.S1412N         Somatic              44      19.08
## 5      4   9 86617335 ATCAGAGAATAGCATT        -   RMI1         NM_024945         Deletion              Frameshift c.1419_1434del   p.473_478del        Germline              39      23.50




#"FANCD1"->"BRAC2"
##    WES.ID Chr Position         Ref.Base Mut.Base   Gene Ref.Transcript.ID Type.of.Mutation Consequence.of.Mutation     DNA.Change Protein.Change Reported.in.FA.mut.DB Reported.in.COSMIC
## 15     32  13 32914891              AAA        - FANCD1         NM_000059         Deletion          Non-Frameshift c.6399_6401del p.2133_2134del                                         
## 17     30   2 58388670              ATA        -  FANCL      NM_001114636         Deletion          Non-Frameshift c.1022_1024del   p.341_342del                     N                  N
## 31      4   9 86617335 ATCAGAGAATAGCATT        -   RMI1         NM_024945         Deletion              Frameshift c.1419_1434del   p.473_478del                     -                  N

miss.genes<-unique(miss.loc[are.missing,"Gene"])

miss.genes[!(miss.genes %in% unique(a.indel[,"Gene.Names"]))]

 ann.cols<-c("chr","start","end","REF","ALT","TYPE","refGene::type","knownGene::type","Gene.Names","Genes.mentioned.at.ASH","refGene::location","OMIM (Gene::Status::OMIM::description::disease)","Consequence.Embl","Uploaded_variation.Embl","Gene.Embl","Feature.Embl", "Protein_position.Embl", "Amino_acids.Embl" , "ensGene::type","ID::maf","FILTER")


 ann.cols<-c("chr","start","end","REF","ALT","TYPE","Gene.Embl","Feature.Embl", "Protein_position.Embl", "Amino_acids.Embl" ,"ID::maf","FILTER")

chk<-a.indel[,"Gene.Names"]=="NPM1" & a.indel[,"ensGene::location"]=="frameshift insertion"

cbind(pass[chk],a.indel[chk,ann.cols])
help[chk,]

chr5:170837544:170837544:-:CTGC:indel
chr5:170837543:170837543:-:TCTG:indel
chr5:170837545:170837545:-:TGCT:indel:170837545
chr5:170819932:170819934:TGC:-:indel

summary.geno.extra[chk,c("GENO.Control","GENO.AML","GENO.PD")]
summary.geno.extra[chk,grepl("^GENO",colnames(summary.geno.extra))]


hits<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/TOP_350.GENOTYPE.conponents-.Burden.clusters.coding.somatic.with.Indels.noBenign.mismatchfilter.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
hits[1:5,1:20]
hits<-unique(hits[,"Gene.Names"])

keep<-as.logical(a.indel[,"wanted.muts"]) | synonymous
sum(keep)

a.indel<-a.indel[keep,]

save(list=c("a.indel"),file="LEO_Pharma_wanted_with_synon_ALL_0.01_muts_QC_MAR31_a.indel.RData" )








######################## TEST SOME LOCATIONS
options(max.print=2000)
options(width=2000)
 ann.cols<-c("chr","start","end","REF","ALT","TYPE","refGene::type","Gene.Embl","Feature.Embl", "Protein_position.Embl", "Amino_acids.Embl" ,"ID::maf","FILTER")


a.indel[a.indel[,"Gene.Names"]=="JAK2",ann.cols]
test<-c("chr9:5073770:5073770:G:T:snp")
annotations[test,]

filter.cols.maf
a.indel[chk,paste(filter.cols.maf,"maf",sep="::")]

                             NHBLI_6500_ANNOVAR_ALL::maf NHBLI_6500_ALL::maf NHBLI_6500_EA::maf NHBLI_6500_AA::maf NHLBI_5400_ALL::maf 1000genome::maf snp141::maf snp137::maf AOGC-NGS_ALL::maf
chr9:5073770:5073770:G:T:snp "0.000231"                  "0.000231"          "0.000233"         "0.000227"         "0.000279"          "0"             "0"         "0"         "0.00201816"     
chr9:5069166:5069166:A:G:snp "0.000077"                  "7.7e-05"           "0.000116"         "0.000000"         "9.3e-05"           "NA"            "0"         "0"         "NA"             
> JAK2 V617F in AOGC cohort

summary.geno.extra[test,]

chk<-c("chr9:5073770:5073770:G:T:snp","chr9:5069166:5069166:A:G:snp")

cbind(a.indel[chk,ann.cols],pass[chk])
cbind(pass[chk],a.indel[chk,ann.cols],help[chk,])
cbind(maf.lt.all[chk,],help[chk,])


chk<-a.indel[,"Gene.Names"]=="NPM1" & a.indel[,"ensGene::location"]=="frameshift insertion"









































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


interesting<-c("TP53","NOTCH1","NOTCH2","FAT1","DGKI","COL19A1") # ,"COL19A1"
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="magenta",pch=20,cex=2)
text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="magenta",pos=4,offset=1,cex=1.3)  #,pch=25)

interesting<-c("STK19","KNSTRN","BCL2L12","CCDC61","EBNA1BP2","PHACTR3","CDKN2A")
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="forestgreen",pch=20,cex=2)
text(qq$x[show],qq$y[show],labels=symbols[qq$ord[show]],col="forestgreen",pos=2,offset=1,cex=1.3)  #,pch=25)

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


savePlot(filename=paste(fig.prefix,".png",sep=''),type="png")
savePlot(filename=paste(fig.prefix,".jpeg",sep=''),type="jpeg")
savePlot(filename=paste(fig.prefix,".tiff",sep=''),type="tiff")
savePlot(filename=paste(fig.prefix,".bmp",sep=''),type="bmp")
dev.print(svg,paste(fig.prefix,".svg",sep=''))



# FREM2 -check why it has moved about
# TNXB has pseudo genes and may be dodgy not caused by rescue but in redo with 27
# UNC80 has moved about chr2:210704125:210704125:C:T:snp
/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/sample_sheet.for.analysis.txt
/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/2015-03-16_AllAMLandLung.BEST.chrALL.ACC_SUBSET.ALL.ALL_GENOTYPES_analysis.txt
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

###### supplement
slide.hits<-c( "INTS1","NFKBIE","PRB3","TNXB","BAD","CHST4","XIRP2","DOCK2")


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









############################## replication selection
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


##########################st jude single sample dump

sample<-"AMAS-13.3-Diagnostic"
pheno.ori[grep(sample,pheno.ori[,"SAMPLE"]),]
sample.bam
SnpStats.bam[grep(sample,SnpStats.bam[,2]),]



to.unwind<-c("Clinical")


to.unwind.name<-to.unwind[1]
#

snpinfo.ex<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,]
loci<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,"Name"] # this is IDH1 not IDH1 in cluster # are the snp.names
loci<-unique(loci)
the.genes<-unique(snpinfo.ex[,"cluster"])
the.genes<-unique(snpinfo.ex[,"gene"])
the.genes<-the.genes[!(the.genes %in% clusters.wanted)]

sort(the.genes)



sample.labels<-colnames(a.indel)[grep(sample,colnames(a.indel))]
sample.labels
annotations[1:5,]



sum(full.qual)
full.qual.test<-full.qual |   ( grepl("^indel",a.indel[,"TYPE"]) & (as.numeric(a.indel[,"AMAS-13.3-Diagnostic.GQ"]) >=90) & !is.na(as.numeric(a.indel[,"AMAS-13.3-Diagnostic.GQ"])) )

sum(full.qual.test)

has.geno<-a.indel[,"AMAS-13.3-Diagnostic.GT"]!="0/0" & a.indel[,"AMAS-13.3-Diagnostic.GT"]!="NA"

wanted.genes<-a.indel[,"Gene.Names"] %in% the.genes
sum(wanted.genes)
sum(has.geno)

sum(wanted.genes & has.geno)

figure<-full.qual.test & has.geno  & maf.filter    & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control & pass.possonian.control.model & !bad.qual.locations   & pass.possonian.aligner.model

sum(figure)

out<-cbind(annotations[figure,],gerp.scores[figure],bad.effect[figure],bad.coding[figure],a.indel[figure,sample.labels],summary.geno.extra[figure,colnames(summary.geno.extra)[grep("^GENO",colnames(summary.geno.extra))]])
out


getwd()
setwd(analysis.dir)
out[order.by,][1:10,1:10]
setwd(analysis.dir)
write.table(out,file=paste(to.unwind.name,"genes",sample,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

figure<-  has.geno  & maf.filter    & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt & !no.genotypes.filt &  rare.in.Control & pass.possonian.control.model & !bad.qual.locations   & pass.possonian.aligner.model

sum(figure)

gatk.filters<-c("FILTER", colnames(a.indel)[1913:1936])

out<-cbind(annotations[figure,],gerp.scores[figure],bad.effect[figure],bad.coding[figure],a.indel[figure,sample.labels],summary.geno.extra[figure,colnames(summary.geno.extra)[grep("^GENO",colnames(summary.geno.extra))]],a.indel[figure,gatk.filters])
a.indel[figure,c(core.ann,"FILTER")]

write.table(out,file=paste(to.unwind.name,"genes",sample,"ALL","txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



colnames(a.indel)[c(1:6,8,11,16,28,7,30,34,35,36,37:42,43,14,32,33)]
annotations<-a.indel[,c(1:6,8,11,16,28,7,30,34,35,36,37:42,43,14,32,33)]
gerp.scores<-a.indel[,"gerp.scores"]


figure<-  has.geno & rare.in.Control 
sum(figure)
dim(a.indel)

colnames(a.indel)[1930:1963]

SVTYPE<-a.indel[,"SVTYPE"]

gatk.filters<-c("FILTER", colnames(a.indel)[1930:1963])

out<-cbind(annotations[figure,],SVTYPE[figure],gerp.scores[figure],bad.effect[figure],bad.coding[figure],a.indel[figure,sample.labels],summary.geno.extra[figure,colnames(summary.geno.extra)[grep("^GENO",colnames(summary.geno.extra))]],a.indel[figure,gatk.filters])
a.indel[figure,c(core.ann,"FILTER")]



write.table(out,file=paste(to.unwind.name,"SV_CNV",sample,"ALL","txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
