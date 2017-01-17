


# project<-"2013-11-17_GroupCallIYMD_Lung_AOGC_Controls" # this is the project directory
# project.name<-"20131117_GroupCallIYMD_Lung_AOGC_Controls" # this is the prefix of the SNP

UQCCG.data<-"/media/UQCCG/Sequencing/Projects"
project<-"BONE-GENOME" # this is the project directory
chrs.to.run<-c(21:24)

iannotate <-7

for(iannotate in 1:length(chrs.to.run)){
  args<-c(chrs.to.run[iannotate],UQCCG.data,project)
  source("/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/Annotate_sequence_run_Aug_2014_BONE.r")

    }




###############################################################





## getwd()
## write.table(t(cbind(column.labels,column.labels)),file="Bone_samples.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



## UQCCG.data<-"/media/UQCCG/Sequencing/Projects"
## project<-"BONE-GENOME" # this is the project directory
## chrs.to.run<-c(1:24)

## iannotate <-1

## for(iannotate in 1:length(chrs.to.run)){
##   args<-c(chrs.to.run[iannotate],UQCCG.data,project,"/media/UQCCG/Sequencing/Projects/BONE-GENOME/MarkersMrOS.csv")
##   source("/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/Analyse_project_new_GENO_DUMP.r")

##     }



UQCCG.data<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/"
project<-"2015-03-16_AllAMLandLung" # this is the project directory
chrs.to.run<-c(5,11)
chrs.to.run<-c(1:24)
iannotate <-1

for(iannotate in 1:length(chrs.to.run)){
  args<-c(chrs.to.run[iannotate],UQCCG.data,project)
  source("/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/Annotate_sequence_run_Aug_2015.r")

    }



UQCCG.data<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/"
project<-"2015-06-04_MODY_NSGC_DISH_PCC" # this is the project directory
chrs.to.run<-c(5,11)
chrs.to.run<-c(1:24)
iannotate <-1

for(iannotate in 1:length(chrs.to.run)){
  args<-c(chrs.to.run[iannotate],UQCCG.data,project)
  source("/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/Annotate_sequence_run_2015.r")

    }


 ##                                   2015-06-04_MODY_NSGC_DISH_PCC.chr10                                    2015-06-04_MODY_NSGC_DISH_PCC.chr11                                    2015-06-04_MODY_NSGC_DISH_PCC.chr12 
## "2015-06-04_MODY_NSGC_DISH_PCC.chr10.output.recalibrated.filtered.vcf" "2015-06-04_MODY_NSGC_DISH_PCC.chr11.output.recalibrated.filtered.vcf" "2015-06-04_MODY_NSGC_DISH_PCC.chr12.output.recalibrated.filtered.vcf" 
##                                    2015-06-04_MODY_NSGC_DISH_PCC.chr13                                    2015-06-04_MODY_NSGC_DISH_PCC.chr14                                    2015-06-04_MODY_NSGC_DISH_PCC.chr15 
## "2015-06-04_MODY_NSGC_DISH_PCC.chr13.output.recalibrated.filtered.vcf" "2015-06-04_MODY_NSGC_DISH_PCC.chr14.output.recalibrated.filtered.vcf" "2015-06-04_MODY_NSGC_DISH_PCC.chr15.output.recalibrated.filtered.vcf" 
##                                    2015-06-04_MODY_NSGC_DISH_PCC.chr16                                    2015-06-04_MODY_NSGC_DISH_PCC.chr17                                    2015-06-04_MODY_NSGC_DISH_PCC.chr18 
## "2015-06-04_MODY_NSGC_DISH_PCC.chr16.output.recalibrated.filtered.vcf" "2015-06-04_MODY_NSGC_DISH_PCC.chr17.output.recalibrated.filtered.vcf" "2015-06-04_MODY_NSGC_DISH_PCC.chr18.output.recalibrated.filtered.vcf" 
##                                    2015-06-04_MODY_NSGC_DISH_PCC.chr19                                     2015-06-04_MODY_NSGC_DISH_PCC.chr1                                    2015-06-04_MODY_NSGC_DISH_PCC.chr20 
## "2015-06-04_MODY_NSGC_DISH_PCC.chr19.output.recalibrated.filtered.vcf"  "2015-06-04_MODY_NSGC_DISH_PCC.chr1.output.recalibrated.filtered.vcf" "2015-06-04_MODY_NSGC_DISH_PCC.chr20.output.recalibrated.filtered.vcf" 
##                                    2015-06-04_MODY_NSGC_DISH_PCC.chr21                                    2015-06-04_MODY_NSGC_DISH_PCC.chr22                                     2015-06-04_MODY_NSGC_DISH_PCC.chr2 
## "2015-06-04_MODY_NSGC_DISH_PCC.chr21.output.recalibrated.filtered.vcf" "2015-06-04_MODY_NSGC_DISH_PCC.chr22.output.recalibrated.filtered.vcf"  "2015-06-04_MODY_NSGC_DISH_PCC.chr2.output.recalibrated.filtered.vcf" 
##                                     2015-06-04_MODY_NSGC_DISH_PCC.chr3                                     2015-06-04_MODY_NSGC_DISH_PCC.chr4                                     2015-06-04_MODY_NSGC_DISH_PCC.chr5 
##  "2015-06-04_MODY_NSGC_DISH_PCC.chr3.output.recalibrated.filtered.vcf"  "2015-06-04_MODY_NSGC_DISH_PCC.chr4.output.recalibrated.filtered.vcf"  "2015-06-04_MODY_NSGC_DISH_PCC.chr5.output.recalibrated.filtered.vcf" 
##                                     2015-06-04_MODY_NSGC_DISH_PCC.chr6                                     2015-06-04_MODY_NSGC_DISH_PCC.chr7                                     2015-06-04_MODY_NSGC_DISH_PCC.chr8 
##  "2015-06-04_MODY_NSGC_DISH_PCC.chr6.output.recalibrated.filtered.vcf"  "2015-06-04_MODY_NSGC_DISH_PCC.chr7.output.recalibrated.filtered.vcf"  "2015-06-04_MODY_NSGC_DISH_PCC.chr8.output.recalibrated.filtered.vcf" 
##                                     2015-06-04_MODY_NSGC_DISH_PCC.chr9                                     2015-06-04_MODY_NSGC_DISH_PCC.chrX                                     2015-06-04_MODY_NSGC_DISH_PCC.chrY 
##  "2015-06-04_MODY_NSGC_DISH_PCC.chr9.output.recalibrated.filtered.vcf"  "2015-06-04_MODY_NSGC_DISH_PCC.chrX.output.recalibrated.filtered.vcf"  "2015-06-04_MODY_NSGC_DISH_PCC.chrY.output.recalibrated.filtered.vcf" 
