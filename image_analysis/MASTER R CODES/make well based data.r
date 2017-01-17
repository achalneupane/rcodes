

################################ this just concatinates the well data into a single file
############################### use Add well based data to sort out missing info etc
############################### for cellomics based data

setwd("/media/scratch/Data/Cellomics/millian")
#list well files
files.wells<-dir(getwd())
files.wells<-files.wells[grep("well",files.wells)]
files.wells

###########################
#concatinate well data
data.well<-{}
for(i in 1:length(files.wells)){
  data<-read.delim(files.wells[i],header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
  data.well<-rbind(data.well,data)
}

dim(data.well)
data.well[1:5,]
 write.table(data.well,"well_based_data.txt",row.names=FALSE,sep="\t")


########################################################
data.well<-{}
unique(data.well[,"BarCode"])
setwd("/media/Bioinform-D/Research/Cellomics/millian/Fawzi Plate reader raw data/siFB1/siFB1 Resazurin Data")
setwd("/media/Bioinform-D/Research/Cellomics/millian/Fawzi Plate reader raw data/siFB2/Resazurin")
setwd("/media/Bioinform-D/Research/Cellomics/millian/Fawzi Plate reader raw data/siFB3/Resazurin")

data.well<-{}
setwd("/media/Bioinform-D/Research/Cellomics/millian/Fawzi Plate reader raw data/siFB2/AK")
setwd("/media/Bioinform-D/Research/Cellomics/millian/Fawzi Plate reader raw data/siFB3/AK")

files<-dir(getwd())
files<-files[-grep(".xls",files)]
## files<-files[-grep(".xlsx",files)]
files
barcodes<-files

files<-files[grep("SI",files)]
barcodes<-gsub("SI fb","si FB",barcodes)

barcodes<-gsub("SI ","si ",barcodes)

files
barcodes
i<-1
data[1:5,]
for( i in 1:length(files)){
   data<-read.delim(files[i],header=F,sep="",skip=1,fill=TRUE,stringsAsFactors=FALSE)
   data[,1]<-gsub(":","",data[,1])
   print(data[1:3,])
   row<-unlist(lapply(strsplit(data[,1],split=""),function(x) x[1]))
   col<-as.numeric(unlist(lapply(strsplit(data[,1],split="\\D"),function(x) x[2])))
  data.well<-rbind(data.well,cbind(barcodes[i],row,col,data[,2]))
 }


dim(data.well)

data.well[1:5,]
setwd("/media/scratch/Data/Cellomics/millian")



colnames(data.well)<-c("BarCode","Row","Col","resazurin")
 write.table(data.well,"resazurin.txt",row.names=FALSE,sep="\t")

colnames(data.well)<-c("BarCode","Row","Col","ak_leakage")
 write.table(data.well,"ak_leakage.txt",row.names=FALSE,sep="\t")
