fam.file<-"AllVertFX_OP.fam" # fam file to change affection 
ped.file<-"ALL_VERT_FX_OP.ped" ## use affection status in this fam file
out.file<-"AllVertFX_OPupdatedPED.fam" ## 

fam.old.file<-paste(fam.file,"OLD",sep="_") ## keep a copy of the original fam file 
working.dir<-"C:/Users/eduncan/Desktop/Phenotypes_AOGC_exome_chip_analysis"
setwd(working.dir)

 

fam <- read.table(fam.file, quote="\"")
write.table(fam,file=fam.old.file,quote=FALSE,col.names=FALSE,row.names=FALSE)
ped <- read.delim(ped.file, header=F)

fam[1:5,]
ped[1:5,]

posns<-match(fam[,1],ped[,1])
missing<-is.na(posns)
sum(missing)

length(posns)
dim(fam)
posns[1:5]

 

ped.reorder<-ped[posns,]
# above same length and in the same order
fam[,6]<-ped.reorder[,6] # change affection status
fam[,5]<-ped.reorder[,5] # change gender status

fam[1:5,]
ped.reorder[1:5,]

write.table(fam,file=fam.file,quote=FALSE,col.names=FALSE,row.names=FALSE)

getwd()

########################################################################

 

