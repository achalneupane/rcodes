
UQCCG.data<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015"
project<-"2015-07-14_MODY_NSGC_DISH_PCC"
the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-07-14_MODY_NSGC_DISH_PCC/BAM/2015-07-14_MODY_NSGC_DISH_PCC-SampleSheet.csv"
chrs.to.run<-c(1:22,"X","Y")

iannotate <-7

for(iannotate in 1:length(chrs.to.run)){
  args<-c(chrs.to.run[iannotate],UQCCG.data,project)
  source("/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/Annotate_sequence_run_Aug_2014.r")

    }


for(iannotate in 1:length(chrs.to.run)){
  args<-c(chrs.to.run[iannotate],UQCCG.data,project,the.sample.sheet)
  quit()
  source("/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/Analyse_project_new_2015.r")

    }








UQCCG.data<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015"
project<-"2015-07-14_MODY_NSGC_DISH_PCC"
the.sample.sheet<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-07-14_MODY_NSGC_DISH_PCC/BAM/2015-07-14_MODY_NSGC_DISH_PCC-SampleSheet.csv"


for(iannotate in 1:length(chrs.to.run)){
  args<-c(UQCCG.data,project,the.sample.sheet)
  quit()
  source("/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/Analyse_project_MERGE.r")

    }


