




#setwd("/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd("/dmf/uqdi/Core_Services/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts")
load("GERP.test.RData")
########################################### GERP  or just use annovar for gerp> 2
############### ADD IF HERE IF H
library(multicore)
library(DBI)
library(RMySQL)
library(IRanges)

source("annotate_SNPs_subroutines.r")
num.cores<-7 # number of threads to run

#db.details<-list(user="gerp2",pass="GeRpUQ8!6", dbname="gerp", host="ga-apps.di.uq.edu.au")
db.details<-list(
  host="di-mysql01.di.uq.edu.au",
  dbname="gerp",
  user="gerp2",
  pass="GeRpUQ8!6"
)
force.GERP.update<-TRUE
files.in.annotate<-dir(getwd())
target<-basename(bim.file)
target
indels[1:5,]
indels[1:5,c("chr","start","end","REF")]

  gerp.data.file<-paste(gsub(".txt$","",target),".GerpPredictions.RData",sep="")


    if(exists("db.details")){db.details<-list(user="gerp2",pass="GeRpUQ8!6", dbname="gerp", host="ga-apps.di.uq.edu.au")}
    gerp.scores<-get.GERP.MULT(db.details,indels[,c("chr","start","end","REF")]) ### gerp subroutine in annotate.snp.subroutines


   length(gerp.scores)
           sum(gerp.scores!=0) ## most should NOT BE ZERO!!!!
           gerp.scores[1:60]

#    save(list=c("gerp.scores"),file=paste(gsub(".txt$","",target),".GerpPredictions.RData",sep=""))


#########################
