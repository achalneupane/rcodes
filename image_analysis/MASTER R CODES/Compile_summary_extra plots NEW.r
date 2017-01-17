setwd("/media/Bioinform-D/Research/Cellomics/Hits Plate 5")
ann<-read.delim("Hits 5 annotation.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
ann<-read.delim("Hits fo Hits 2B annotation.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
colnames(ann)[1]<-"plate"
colnames(ann)[4]<-"Symbol"
ann[1:5,]


ann<-read.delim("Ca Gene set with controls.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
######## run through below and annotaion are ok
ann.new<-{}
plate1.1<-ann[ann[,"plate"]==1,]
dim(plate1.1)
plate1.2<-plate1.1
plate1.3<-plate1.1
plate1.1[,"plate"]<-1.1
plate1.2[,"plate"]<-1.2
plate1.3[,"plate"]<-1.3
ann.new<-rbind(ann.new,plate1.1,plate1.2,plate1.3)

plate1.1<-ann[ann[,"plate"]==2,]
dim(plate1.1)
plate1.2<-plate1.1
plate1.3<-plate1.1
plate1.1[,"plate"]<-2.1
plate1.2[,"plate"]<-2.2
plate1.3[,"plate"]<-2.3
ann.new<-rbind(ann.new,plate1.1,plate1.2,plate1.3)

plate1.1<-ann[ann[,"plate"]==3,]
dim(plate1.1)
plate1.2<-plate1.1
plate1.3<-plate1.1
plate1.1[,"plate"]<-3.1
plate1.2[,"plate"]<-3.2
plate1.3[,"plate"]<-3.3
ann.new<-rbind(ann.new,plate1.1,plate1.2,plate1.3)
ann.old<-ann
ann <- ann.new
save(list=c("ann"),file="Ca_ann.RData")



ann<-read.delim("siCA3 plate reader data.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)


colnames(ann)[c(1,4,5)]<-c("plate","row","col")
ann[,"plate"]<-gsub(">","",ann[,"plate"])
### run tests below
save(list=c("ann"),file="Leo_ann.RData")

setwd("/media/Bioinform-D/Research/Cellomics/Ca screen/Latest")
write.table(ann,"Ca_annotations2.txt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
ann<-read.delim("Ca_annotations2.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
load("Ca_ann2.RData")
ann<-read.delim("Helga_annotations.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

ann3<-ann
ann3[,"plate"]<-3
ann4<-ann
ann4[,"plate"]<-4
ann5<-ann
ann5[,"plate"]<-5
ann<-rbind(ann3,ann4,ann5)
write.table(ann,"Helga_annotations2.txt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

##############################repealted plates
files<-dir(getwd())
files<-files[grepl("SUMMARY",files)]
files
ann<-read.delim("LizHitsPlateMap.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
ann[1:5,]
ann<-cbind("1a",ann,stringsAsFactors=FALSE)
colnames(ann)<-c("plate","row","col","Symbol")
######### go below and sort this data
files<-strsplit(files,split="_")
plates<-unlist(lapply(files,function(x) x[4]))
TREAT<-strsplit(plates,split="",perl=TRUE)
TREAT<-unlist(lapply(TREAT,function(x) x[2]))

ann.new<-{}
for(i in 1:length(plates)){
  temp<-cbind(plates[i],ann[,c("row","col","Symbol")],TREAT[i],stringsAsFactors=FALSE)
  ann.new<-rbind(ann.new,temp)
}
colnames(ann.new)<-c("plate","row","col","Symbol","TREAT")
ann<-ann.new

###########################


colnames(ann)<-c("row","col","Symbol")
ann<-cbind(plate=1,ann)

setwd("/media/Bioinform-D/Research/Cellomics/Jane")
ann<-read.delim("jane.ann3.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
ann<-ann[,1:10]
colnames(ann)[c(2,7,8,9)]<-c("plate","row","col","Symbol")
tapply(ann[,"plate"],ann[,"plate"],length)
tapply(ann[,"row"],ann[,"row"],length)
tapply(ann[,"col"],ann[,"col"],length)
tapply(ann[,"BarCode"],ann[,"BarCode"],length)

barcodes<-unique(ann[,"BarCode"])
nrows
## ann.plate<-ann
a.plate<-matrix("MOCK",nrow=col.type*row.type,ncol=dim(ann)[2])
colnames(a.plate)<-colnames(ann)
a.plate[,"row"]<-ann.plate[,"row"]
a.plate[,"col"]<-ann.plate[,"col"]
a.plate[( (a.plate[,"row"] %in% c("A","B","C")) & (a.plate[,"col"] %in% c("6"))),"Symbol"]<-"CCND1"
a.plate[( (a.plate[,"row"] %in% c("D","E","F")) & (a.plate[,"col"] %in% c("6"))),"Symbol"]<-"HSPCN111-N"
a.plate[( (a.plate[,"row"] %in% c("G","H")) & (a.plate[,"col"] %in% c("6"))),"Symbol"]<-"MOCK"

a.plate[( (a.plate[,"row"] %in% c("A","B","C")) & (a.plate[,"col"] %in% c("12"))),"Symbol"]<-"HSPCN111-N"
a.plate[( (a.plate[,"row"] %in% c("D","E","F")) & (a.plate[,"col"] %in% c("12"))),"Symbol"]<-"CCND1"
a.plate[( (a.plate[,"row"] %in% c("G","H")) & (a.plate[,"col"] %in% c("12"))),"Symbol"]<-"UNTRANS"
a.plate

tapply(a.plate[,"Symbol"],a.plate[,"Symbol"],length)
a.plate[,c("row","col","Symbol")]

load("working.RData")


split<-strsplit(barcodes,split="\\.")
plate<-unlist(lapply(split,function(x) x[3]))
barcodes<-sort(barcodes)
ann.all<-{}
for(i in 1:length(barcodes)){
ann.all<-rbind(ann.all,cbind(barcodes[i],plate[i],a.plate))
}
colnames(ann.all)[1:2]<-c("BarCode","plate")
ann.all<-ann.all[,-10]
ann.all<-as.data.frame(ann.all,stringsAsFactors=FALSE)
ann<-as.data.frame(ann,stringsAsFactors=FALSE)

ann.old<-ann
ann<-ann.all
colnames(ann)[1:2]<-c("plate","core.plate")  ## make the Barcode the plate nname
colnames(ann.old)[3]<-c("plate")

have<-paste(ann.old[,"plate"],ann.old[,"row"],ann.old[,"col"],sep="::")
got<-paste(ann[,"plate"],ann[,"row"],ann[,"col"],sep="::")

ann[ann[,"plate"]=="MUS.a.1",]
posns<-match(have,got)
missing<-is.na(posns)
sum(missing)
ann[posns,c("Accessions","GeneID","MusIDT","M.position","Symbol" ,"TREAT") ]<-ann.old[,c("Accessions","GeneID","MusIDT","M.position","Symbol" ,"TREAT")]


 well.type<-96
row.type<-8
col.type<-12

> tapply(ann[,"plate"],ann[,"plate"],length)
 1  2  3  4  5  6  7 
80 80 80 80 80 49 70 
> tapply(ann[,"row"],ann[,"row"],length)
 a  A  b  B  c  C  d  D  e  E  f  F  g  G  h  H 
 3 64  3 63  3 63  3 63  3 63  3 63  3 58  3 58 
> tapply(ann[,"col"],ann[,"col"],length)
 1  2  3  4  5  7  8  9 10 11 
54 54 54 54 54 49 56 48 48 48 
> plates

ann[ann[,"row"]=="a","row"]<-"A"
ann[ann[,"row"]=="b","row"]<-"B"
ann[ann[,"row"]=="c","row"]<-"C"
ann[ann[,"row"]=="d","row"]<-"D"
ann[ann[,"row"]=="e","row"]<-"E"
ann[ann[,"row"]=="f","row"]<-"F"
ann[ann[,"row"]=="g","row"]<-"G"
ann[ann[,"row"]=="h","row"]<-"H"
tapply(ann[,"plate"],ann[,"plate"],length)
tapply(ann[,"row"],ann[,"row"],length)
tapply(ann[,"col"],ann[,"col"],length)



posns<-strsplit(ann[,"X384.Absolute.Well.ID"],split="",)






plate<-unlist(lapply(posns,function(x) x[1]))
row<-unlist(lapply(posns,function(x) x[2]))
posns<-strsplit(ann[,"X384.Absolute.Well.ID"],split="\\D",)
col<-as.numeric(unlist(lapply(posns,function(x) x[2])))

ann<-cbind(plate,row,col,ann)
cont<-read.delim("millian.controls.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
colnames(cont)[c(2,3,4)]<-c("plate","row","col")
 tapply(cont[,"Row"],cont[,"Row"],length)
 tapply(cont[,"siRNA.Plate"],cont[,"siRNA.Plate"],length)
 tapply(cont[,"Row"],cont[,"Row"],length)
cont.have<-paste(cont[,"plate"],cont[,"row"],cont[,"col"],sep="")
length(cont.have)
length(unique(cont.have))
posns<-match(unique(cont.have),cont.have)
sum(is.na(posns))
cont<-cont[posns,]

dim(ann)
colnames(ann)[7]<-"Symbol"
ann[1,]
ann.cont<-ann[1:dim(cont)[1],]
ann.cont<-apply(ann.cont,2,function(x) x<-rep(NA,time=length(x)))
ann.cont[,"plate"]<-cont[,"plate"]
ann.cont[,"row"]<-cont[,"row"]
ann.cont[,"col"]<-cont[,"col"]
ann.cont[,"symbol"]<-cont[,"Gene.Name"]

ann<-rbind(ann,ann.cont)
ann[ann[,"row"]=="I",1:4]
ann<-ann[-1040,]


ann.ori<-ann
files<-dir(getwd())
files<-files[4:27]
files<-gsub("plate_","",files)
files<-gsub("_SUMMARY.RData","",files)
files
the.plates<-strsplit(files,split="_")
the.plates<-lapply(the.plates,function(x) x[2])
the.plates<-unlist(the.plates)

the.plates<-strsplit(the.plates,split="\\.")
reps<-unlist(lapply(the.plates,function(x) x[2]))
plates<-unlist(lapply(the.plates,function(x) x[1]))


> plates
 [1] "1" "1" "1" "1" "2" "2" "2" "2" "3" "3" "3" "3" "1" "1" "1" "2" "2" "2" "3"
[20] "3" "3" "1" "1" "1" "2" "2" "2" "3" "3" "3"
> files
 [1] "si FB1_1.1" "si FB1_1.2" "si FB1_1.3" "si FB1_1.4" "si FB1_2.1"
 [6] "si FB1_2.2" "si FB1_2.3" "si FB1_2.4" "si FB1_3.1" "si FB1_3.2"
[11] "si FB1_3.3" "si FB1_3.4" "si FB2_1.1" "si FB2_1.2" "si FB2_1.3"
[16] "si FB2_2.1" "si FB2_2.2" "si FB2_2.3" "si FB2_3.1" "si FB2_3.2"
[21] "si FB2_3.3" "si FB3_1.1" "si FB3_1.2" "si FB3_1.3" "si FB3_2.1"
[26] "si FB3_2.2" "si FB3_2.3" "si FB3_3.1" "si FB3_3.2" "si FB3_3.3"
> reps
 [1] "1" "2" "3" "4" "1" "2" "3" "4" "1" "2" "3" "4" "1" "2" "3" "1" "2" "3" "1"
[20] "2" "3" "1" "2" "3" "1" "2" "3" "1" "2" "3"
>

plate.locs<-ann[,1:length(unique(ann[,"plate"]))]
for(i in 1:dim(plate.locs)[2]){ plate.locs[,i]<-ann[,"plate"]==i }
plate.locs[,1]<-ann[,"plate"]==1
plate.locs[,2]<-ann[,"plate"]==2
plate.locs[,3]<-ann[,"plate"]==3

ann.total<-{}
for(i in 1:length(reps)){
   ann.total<-rbind(ann.total,cbind(files[i],ann[ann[,"plate"]==plates[i],]))
 }
ann<-ann.total
colnames(ann)[1:2]<-c("plate","core.plate")

save(list=c("ann"),file="Jane_ann.RData")
load("Jane_ann.RData")
colnames(ann)[1:2]<-c("plate","core.plate")


expt.a<-grepl("MUS.a",ann[,"plate"])
ann<-ann[expt.a,]

load("Jane_validation_ann.RData") 
test<-read.delim("mtoSymbol2.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
ann[1:5,]
ann[,"M.position"]<-ann[,"Symbol"]
ann[,"Accessions"]<-NA
ann[,"GeneID"]<-NA
ann[1:5,]
test[1:50,]
for(i in 1:dim(test)[1]){

 posns<-grep(paste("^",test[i,"M.position"],"$",sep=""),ann[,"M.position"])
 missing<-is.na(posns)
 print(i)
print( posns[!missing])
 if(length(posns)>0){
ann[posns[!missing],"Symbol"]<-test[i,"Symbol"]
ann[posns[!missing],"Accessions"]<-test[i,"Accessions"]
ann[posns[!missing],"GeneID"]<-test[i,"GeneID"]
}
}
 
ann[posns[!missing],]
 
test<-read.delim("MUS JANE 3 plate map.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
colnames(test)[c(3,6,7,8)]<-c("Symbol","plate","row","col")

test[,"plate"]<-paste("MUS.a.",test[,"plate"],sep="")

tapply(ann[,"plate"],ann[,"plate"],length)
tapply(ann[,"row"],ann[,"row"],length)
tapply(ann[,"col"],ann[,"col"],length)

tapply(test[,"plate"],test[,"plate"],length)
tapply(test[,"row"],test[,"row"],length)
tapply(test[,"col"],test[,"col"],length)


given<-paste(test[,"plate"],test[,"row"],test[,"col"],sep="::")
got<-paste(ann[,"plate"],ann[,"row"],ann[,"col"],sep="::")


posns<-match(got,given)
missing<-is.na(posns)
sum(missing)

ann[1:15,core.ann]
test[posns[1:15],core.ann]
test<-test[posns,]

ann[1:15,core.ann]
test[1:15,core.ann]

diff<-ann[,"M.position"] != test[,"M.position"]

a.control<-ann[,"Symbol"] %in% c("UNTRANSDUCED","CCND1","HSPCN111-N","MOCK")

diff<- diff & !a.control
sum(diff)

cbind(ann[diff,core.ann],test[diff,core.ann])

tapply(ann[diff,"plate"],ann[diff,"plate"],length)
tapply(ann[diff,"row"],ann[diff,"row"],length)
tapply(ann[diff,"col"],ann[diff,"col"],length)

tapply(test[,"plate"],test[,"plate"],length)
tapply(test[,"row"],test[,"row"],length)
tapply(test[diff,"col"],test[diff,"col"],length)
sort(tapply(test[diff,"Symbol"],test[diff,"Symbol"],length))


core.ann<-c(core.ann,"M.position","Accessions")
ann[posns,c("Accessions","GeneID","MusIDT","M.position","Symbol" ,"TREAT") ]<-ann.old[,c("Accessions","GeneID","MusIDT","M.position","Symbol" ,"TREAT")]

        B3GALNT1          FASTKD2 gene:12057:::---          GPATCH4             IRS1            MTP18            SFRS2             SHQ1 
               1                1                1                1                1                1                1                1 
         STARD10            WDR75     UNTRANSDUCED            CCND1       HSPCN111-N             MOCK 
               1                1               16               48               48               62 

test[diff,core.ann]
      plate row col           Symbol M.position M.position.1  Accessions
305 MUS.a.4   B   2            WDR75       M4B2         M4B2    BC040567
306 MUS.a.4   B   3             IRS1       M4B3         M4B3    BC053895
317 MUS.a.4   C   2            SFRS2       M4C2         M4C2    BC066958
343 MUS.a.4   E   4             SHQ1       M4E4         M4E4    BC047879
338 MUS.a.4   E  10            MTP18      M4E10        M4E10    BC001608
350 MUS.a.4   F  10          STARD10      M4F10        M4F10    BC007919
368 MUS.a.4   G   5          GPATCH4       M4G5         M4G5    BC056904
373 MUS.a.4   H   1         B3GALNT1       M4H1         M4H1    BC047618
660 MUS.a.7   G   9          FASTKD2       M7G9         M7G9    BC001544
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
## plate checks

well.type<-384
row.type<-16
col.type<-24

well.type<-96
row.type<-8
col.type<-12

row.index<-c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")

tapply(ann[,"plate"],ann[,"plate"],length)
tapply(ann[,"row"],ann[,"row"],length)
tapply(ann[,"col"],ann[,"col"],length)
tapply(ann[,"Symbol"],ann[,"Symbol"],length)

plate.ids<-unique(ann[,"plate"])
plate.ids<-sort(plate.ids)
for(i in 1:length(plate.ids)){
 print(plate.ids[i])
  print(tapply( ann[ann[,"plate"]==plate.ids[i],"row"],ann[ann[,"plate"]==plate.ids[i],"row"],length)) }
############### test plate nummbering
## plate.ids<-unique(ann[,"plate"])
plate.ids.t<-rep(plate.ids,well.type)
dim(plate.ids.t)<-c(length(plate.ids),well.type)
plate.ids.t<-t(plate.ids.t)
plate.ids.t<-as.vector(plate.ids.t)

sum(plate.ids.t!=ann[,"plate"])  # must be zero
####################################

############### test row nummbering
row.ids<-row.index[1:row.type]
row.ids.t<-rep(row.ids,col.type)
dim(row.ids.t)<-c(length(row.ids),col.type)
row.ids.t<-t(row.ids.t)
row.ids.t<-as.vector(row.ids.t)
row.ids.t<-rep(row.ids.t,length(plate.ids))

sum(row.ids.t!=ann[,"row"]) # must be zero
test<-row.ids.t!=ann[,"row"]
ann[test,1:4]

###################################
 ############### test row nummbering
col.ids<-1:col.type
col.ids.t<-rep(col.ids,row.type)
col.ids.t<-rep(col.ids.t,length(plate.ids))

sum(col.ids.t!=ann[,"col"]) # must be zero
test<-col.ids.t!=ann[,"col"]
ann[test,1:4]

############################################################
##################################### FIX row and column ordering
############################################################
order.wanted<-paste(plate.ids.t,row.ids.t,col.ids.t,sep="::")
order.have<-paste(ann[,"plate"],ann[,"row"],ann[,"col"],sep="::")
posns<-match(order.wanted,order.have)
sum(is.na(posns)) # must be zero

ann<-ann[posns,]
################## now go back to above and run checks again

############################################################
#
################## now go back to above and run checks again
target<-"GeneID"
sum(ann1[,target] !=ann[,target])

to.fix<-ann1[,target] !=ann[,target]

cbind(ann1[to.fix,core.ann],ann[to.fix,core.ann])

ann1[to.fix,target]<-ann[to.fix,target]

colnames(ann)[8]<-"Symbol"
ann[,"plate"]<-"2B"



###################################
setwd( "/media/Bioinform-D/Research/Cellomics/Joseph screen")
setwd( "/media/Bioinform-D/Research/Cellomics/Kiril screen")
save(list=c("ann"),file="annPlate5.RData")
save(list=c("ann"),file="annPlate2B.RData")
save(list=c("ann"),file="Ca3_ann.RData")
save(list=c("ann"),file="Jane_validation_ann.RData")
############################################### Plate 2B
setwd("/media/Bioinform-D/Research/Cellomics/Hits Plate 5")
load("annPlate2B.RData")

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S") # only the numerical ones


the.screen<-"Plate2B"
files<- paste("plate_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Plate2B_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Plate2B_summary_NOTGREEN_DNA.txt"
#######################################################

############################################### Plate 5
setwd("/media/Bioinform-D/Research/Cellomics/Hits Plate 5")
load("annPlate5.RData")

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S")


the.screen<-"Plate5"
files<- paste("plate_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Plate5_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Plate5_summary_NOTGREEN_DNA.txt"
#######################################################


############################################### Leo  ### latest DNA and flexible annotation
setwd( "/media/Bioinform-D/Research/Cellomics/Leo screen")
load("Leo_ann.RData")
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann
sum.vars<-c(core.ann,"Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S",colnames(ann)[!place.core.ann])

the.screen<-"Leo"
files<- paste("Plates_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Leo_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Leo_summary_NOTGREEN_DNA.txt"
#######################################################



############################################### inhibitor
setwd( "/media/Bioinform-D/Research/Cellomics/inhibitor")
load("ann3000.RData")
ann<-ann[1:96,]  # fake annotation one plate
sum.vars<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S","Description","Horfeome.position", "MGC.Accession","Pilot.Position","Gene.Name.pilotPos","Kiril.position")


the.screen<-"Inhibitor"
files<- paste("plate_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Inhibitor_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Inhibitor_summary_NOTGREEN_DNA.txt"
#######################################################


############################################### Ca screen
setwd( "/media/Bioinform-D/Research/Cellomics/Ca screen")
load("Ca_ann.RData")


sum.vars<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S","Entrez_ID","Cell_loc","Description","Gene.ID","Accession.Number","GI.Number","Library.Type")

the.screen<-"Ca siRNA"
files<- paste("Ca_siCa_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Ca_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Ca_summary_NOTGREEN_DNA.txt"
well.type<-384
row.type<-16
col.type<-24
#######################################################


############################################### Ca screen 2
setwd( "/media/Bioinform-D/Research/Cellomics/Ca screen/Ca screen 2")
#load("Ca_ann2.RData") # same as Ca_ann2 except 
map.file<-"cellomics.to.expt.map.csv" # "plate" must be same "plate in annoation AND IN THE SAME ORDER
annotations.file<-"Ca_annotations2.txt"

sum.vars<-c("plate","row","col","Symbol","%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S","Entrez_ID","Cell_loc","Description","Gene.ID","Accession.Number","GI.Number","Library.Type")

the.screen<-"Ca siRNA2"
#files<- paste("Ca_siCa_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="") # defined below
field.output.file<-"Ca_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Ca_summary_NOTGREEN_DNA.txt"
well.type<-384
row.type<-16
col.type<-24

#######################################################


############################################### Ca screen 2 " same as above but the very latest
setwd( "/media/Bioinform-D/Research/Cellomics/Ca screen/Ca screen 2")
load("Ca2_ann.RData")  ### "plate" must match "plate" in file  cellomics to expt map.csv
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","G1","G2","S","Gr-S+G2/G1","Gr-S+G2+ >4N/G1","Gr-<2N","Gr->4N","Gr-G1","Gr-G2","Gr-S","NtGr-S+G2/G1","NtGr-S+G2+ >4N/G1","NtGr-<2N","NtGr->4N","NtGr-G1","NtGr-G2","NtGr-S" )


the.screen<-"Ca siRNA2"
files<- paste("plate_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Ca2_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Ca2_summary_NOTGREEN_DNA.txt"
well.type<-384
row.type<-16
col.type<-24

#######################################################

###############################################Ca screen 3  ### latest DNA and flexible annotation
setwd("/media/Bioinform-D/Research/Cellomics/Ca screen/Ca screen 3")
load("Ca3_ann.RData")  ### "plate" must match "plate" in file  cellomics to expt map.csv
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","G1","G2","S","Gr-S+G2/G1","Gr-S+G2+ >4N/G1","Gr-<2N","Gr->4N","Gr-G1","Gr-G2","Gr-S","NtGr-S+G2/G1","NtGr-S+G2+ >4N/G1","NtGr-<2N","NtGr->4N","NtGr-G1","NtGr-G2","NtGr-S"
             )

the.screen<-"Ca siRNA3"
files<- paste("plate_siCA3_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Ca3_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Ca3_summary_NOTGREEN_DNA.txt"
well.type<-384
row.type<-16
col.type<-24
#######################################################


############################################### Leo  ### latest DNA and flexible annotation
setwd( "/media/Bioinform-D/Research/Cellomics/Leo screen")
load("Leo_ann.RData")
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann
sum.vars<-c(core.ann,"Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S",colnames(ann)[!place.core.ann])

the.screen<-"Leo"
files<- paste("Plates_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Leo_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Leo_summary_NOTGREEN_DNA.txt"
#######################################################

############################################### Helga
setwd( "/media/Bioinform-D/Research/Cellomics/Helga")
#load("Ca_ann2.RData") # same as Ca_ann2 except 
map.file<-"cellomics.to.expt.map.csv" # "plate" must be same "plate in annoation AND IN THE SAME ORDER
annotations.file<-"Helga_annotations2.txt"

sum.vars<-c("plate","row","col","Symbol","%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S","Entrez_ID","Cell_loc","Description","Gene.ID","Accession.Number","GI.Number","Library.Type")

the.screen<-"Helga"
#files<- paste("Ca_siCa_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="") # defined below
field.output.file<-"Helga_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Helga_summary_NOTGREEN_DNA.txt"
well.type<-96
row.type<-8
col.type<-12

#######################################################
############################################### Kiril  ### latest DNA and flexible annotation
setwd("/media/Bioinform-D/Research/Cellomics/Kiril screen")
load("ann2800.RData")  ### "plate" must match "plate" in file  cellomics to expt map.csv
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field")

the.screen<-"kiril"
files<- paste("Kiril_Plate",unique(ann[,"plate"]),"_SUMMARY1",".RData",sep="")
field.output.file<-"Kiril_field_summary_NOTGREEN_DNA_2010.txt"
well.output.file<-"Kiril_summary_NOTGREEN_DNA_2010.txt"
well.type<-96
row.type<-8
col.type<-12
#######################################################


############################################### Kiril  ### latest DNA and flexible annotation
setwd("/media/Bioinform-D/Research/Cellomics/Joseph screen")
load("ann1500.RData")  ### "plate" must match "plate" in file  cellomics to expt map.csv
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field")

the.screen<-"Joseph"
files<- paste("Joseph_Plate_",unique(ann[,"plate"]),"_SUMMARY1",".RData",sep="")
field.output.file<-"Joseph_field_summary_NOTGREEN_DNA_2010.txt"
well.output.file<-"Joesph_summary_NOTGREEN_DNA_2010.txt"
well.type<-96
row.type<-8
col.type<-12
#######################################################


 # annotation for KIRIL screen
# annotation for kiril screen
#sum.vars<-c("plate","row","col","Symbol","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","%R","R","RG_Low","RG_Mid","RG_High","G","Num-cells","score","RG/G","RnG/nG","RnG","nR","nG","gt4N","Total","Trans_%","<cells_field>","sd_cells_field","Description","Horfeome.position", "MGC.Accession","Pilot.Position","Gene.Name.pilotPos","Kiril.position")
load("ann3000.RData")
sum.vars<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","Description","Horfeome.position", "MGC.Accession","Pilot.Position","Gene.Name.pilotPos","Kiril.position")
  #kiril


files<- paste("Kiril_Plate",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Kiril_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Kiril_summary_NOTGREEN_DNA.txt"


load("ann1500.RData")
sum.vars<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","Description","GenBank_Accn","dbEST_Id" ,"MGC.Accession")



 #joesph
files<- paste("Joseph_plate_",unique(ann[,"plate"]),"_SUMMARY1",".RData",sep="")
field.output.file<-"Joseph_field_summary_NOTGREEN_DNA1.txt"
well.output.file<-"Joseph_summary_NOTGREEN_DNA1.txt"

setwd( "/media/Bioinform-D/Research/Cellomics/Hugo screen")
load("ann.RData")
sum.vars<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","Description","GenBank_Accn","dbEST_Id" ,"MGC.Accession")
 #joesph
files<- paste("Hugo_hugoA_00",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Hugo_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Hugo_summary_NOTGREEN_DNA.txt"
the.screen<-"Hugo"

# summ<-c(rep.int("0",dim(ann)[1]*length(sum.vars)))



############################################### millian  ### latest DNA and flexible annotation
setwd("/media/Bioinform-D/Research/Cellomics/millian")
load("siFB_ann.RData")  ### "plate" must match "plate" in file  cellomics to expt map.csv
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","G1","G2","S","Gr-S+G2/G1","Gr-S+G2+ >4N/G1","Gr-<2N","Gr->4N","Gr-G1","Gr-G2","Gr-S","NtGr-S+G2/G1","NtGr-S+G2+ >4N/G1","NtGr-<2N","NtGr->4N","NtGr-G1","NtGr-G2","NtGr-S"
             )

the.screen<-"siFB"
files<- paste("plate_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"siFB_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"siFB_summary_NOTGREEN_DNA.txt"
well.type<-384
row.type<-16
col.type<-24
#######################################################

############################################### Kiril 2010  ### latest DNA and flexible annotation
setwd("/media/Bioinform-D/Research/Cellomics/Kiril screen 2010")
load("ann2800.RData")  ### "plate" must match "plate" in file  cellomics to expt map.csv
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","G1","G2","S","Gr-S+G2/G1","Gr-S+G2+ >4N/G1","Gr-<2N","Gr->4N","Gr-G1","Gr-G2","Gr-S","NtGr-S+G2/G1","NtGr-S+G2+ >4N/G1","NtGr-<2N","NtGr->4N","NtGr-G1","NtGr-G2","NtGr-S"
             )

the.screen<-"Kiril"
files<- paste("Kiril_Plate",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Kiril_2010_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Kiril_2010_summary_NOTGREEN_DNA.txt"
well.type<-96
row.type<-8
col.type<-12
#######################################################
############################################### JOseph 2010  ### latest DNA and flexible annotation
setwd("/media/Bioinform-D/Research/Cellomics/Joseph screen 2010")
load("ann1500.RData")  ### "plate" must match "plate" in file  cellomics to expt map.csv
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","G1","G2","S","Gr-S+G2/G1","Gr-S+G2+ >4N/G1","Gr-<2N","Gr->4N","Gr-G1","Gr-G2","Gr-S","NtGr-S+G2/G1","NtGr-S+G2+ >4N/G1","NtGr-<2N","NtGr->4N","NtGr-G1","NtGr-G2","NtGr-S"
             )

the.screen<-"Joseph"
files<- paste("plate_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Joseph_2010_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Joseph_2010_summary_NOTGREEN_DNA.txt"
well.type<-96
row.type<-8
col.type<-12
######################################################


############################################### Jane  ### latest DNA and flexible annotation
setwd("/media/Bioinform-D/Research/Cellomics/Jane")
load("Jane_ann.RData")  ### "plate" must match "plate" in file  cellomics to expt map.csv
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","G1","G2","S","Gr-S+G2/G1","Gr-S+G2+ >4N/G1","Gr-<2N","Gr->4N","Gr-G1","Gr-G2","Gr-S","NtGr-S+G2/G1","NtGr-S+G2+ >4N/G1","NtGr-<2N","NtGr->4N","NtGr-G1","NtGr-G2","NtGr-S"
             )

the.screen<-"Jane"
files<- paste("plate_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Jane_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Jane_summary_NOTGREEN_DNA.txt"
well.type<-96
row.type<-8
col.type<-12
#######################################################

############################################### Jane VALIDATION  ### latest DNA and flexible annotation
setwd("/media/Bioinform-D/Research/Cellomics/Jane/Validation")
load("Jane_validation_ann.RData")  ### "plate" must match "plate" in file  cellomics to expt map.csv
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","G1","G2","S","Gr-S+G2/G1","Gr-S+G2+ >4N/G1","Gr-<2N","Gr->4N","Gr-G1","Gr-G2","Gr-S","NtGr-S+G2/G1","NtGr-S+G2+ >4N/G1","NtGr-<2N","NtGr->4N","NtGr-G1","NtGr-G2","NtGr-S"
             )

the.screen<-"Jane_validation"
files<- paste("plate_LIZ_HITS_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Jane_validation_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Jane_validation_summary_NOTGREEN_DNA.txt"
well.type<-96
row.type<-8
col.type<-12
#######################################################

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
################################################################## START
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

########################## the map and key file strategy is new:
map<-read.delim(map.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
ann<-read.delim(annotations.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
### cleck annotation file comtains all plates in map file but map file may conatin more entries than the annotation file
(sum(unique(ann[,"plate"]) %in% map[,"Plate"]) ==  length(unique(ann[,"plate"]))) # MUST BE TRUE
#(sum(unique(ann[,"plate"]) %in% map[,"Plate"]) ==  length(unique(ann[,"plate"]))) & (length(unique(ann[,"plate"]))== length(map[,"Plate"])) # MUST BE TRUE



#########################for flexible annotations IF a map file exits:
map<-read.delim(map.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
posns<-match(unique(ann[,"plate"]),map[,"Plate"])
posns<-posns[!is.na(posns)]
map<-map[posns,]
files<- paste("plate_",map[,"BarCode"],"_SUMMARY",".RData",sep="")
files
##########################################################################

###################### for flexible annotation START: begin immediately files and ann loaded already
files
unique(ann[,"plate"])

row.index<-c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
### get the samples in the column names
####################### NOT THIS FAILS IF ANN NOT IN ROW WISE ORDER!
summ={}
cell.counts={}
counter<-0
n<-4 #number of significant figures to use

for(i in 1:length(files)) {
RG.P.field<-{}# dummy cause this variable not awyas available!
load(files[i])
print(i)

#(x[ii]-n[ii])/sd.expected[ii]
# get prolif probablitity from well est and then just use green cells:
#prolif.prob<-(red.c/(all.c-big.c))
prolif.prob<-(redAndNotGreen/(notGreen.c))
percent.red<-(red.c/(all.c-big.c))

stand.dev<-sqrt(green.c*prolif.prob*(1-prolif.prob))    #sigma=sqrt(var)=sqrt( np(1-p) )
mean=green.c*prolif.prob  # mean=n*p
all.prob<- (redAndGreen-mean)/stand.dev # z=x-m/s

stand.devHigh<-sqrt(greenHigh*prolif.prob*(1-prolif.prob))    #sigma=sqrt(var)=sqrt( np(1-p) )
meanHigh=greenHigh*prolif.prob
all.probHigh<- (redAndGreenHigh-meanHigh)/stand.devHigh

stand.devMid<-sqrt(greenMid*prolif.prob*(1-prolif.prob))    #sigma=sqrt(var)=sqrt( np(1-p) )
meanMid=greenMid*prolif.prob
all.probMid<- (redAndGreenMid-meanMid)/stand.devMid

stand.devLow<-sqrt(greenLow*prolif.prob*(1-prolif.prob))    #sigma=sqrt(var)=sqrt( np(1-p) )
meanLow=greenLow*prolif.prob
all.probLow<- (redAndGreenLow-meanLow)/stand.devLow


#### REMEMBER arrays fill columnwise first in R
test.same<- c(rep(FALSE,length(cells.P.field)) )
row.labels<-rep(rownames(red.c),dim(red.c)[2])
col.labels<-rep(colnames(red.c),dim(red.c)[1])
dim(col.labels)<-dim(t(cells.P.field))
col.labels<-as.integer(t(col.labels))

### test that there are the same number of fields for all quantities (can be different as filed with too few green or red cells are skipped sometimes
for(j in 1:length(cells.P.field)) {
test.same[j]<- (length(as.numeric(unlist(strsplit(cells.P.field[[j]],",")))) == length(as.numeric(unlist(strsplit(R.P.field[[j]],",")))) ) & ( length(as.numeric(unlist(strsplit(cells.P.field[[j]],",")))) == length(as.numeric(unlist(strsplit(G.P.field[[j]],","))))  )
}


### need as vectore aer strsplit etc can return a matrix or a vector from unlist
cells.f<-as.vector(unlist(lapply(cells.P.field[test.same],function(x) { as.numeric(unlist(strsplit(x,",")))  }  ) ))
r.cells.f<-as.vector(unlist(lapply(R.P.field[test.same],function(x) { as.numeric(unlist(strsplit(x,",")))  }  ) ))
g.cells.f<-as.vector(unlist(lapply(G.P.field[test.same],function(x) { as.numeric(unlist(strsplit(x,",")))  }  ) ))

if(!is.null(dim(RG.P.field)[1])){
  rg.cells.f<-as.vector(unlist(lapply(RG.P.field[test.same],function(x) { as.numeric(unlist(strsplit(x,",")))  }  ) ))
}else{
 rg.cells.f <-rep(0,times=length(r.cells.f)) ###not correct but will do for now only effect field data
}

if(!exists("RnG.P.field")){rng.cells.f <-r.cells.f-rg.cells.f}else{
if(!is.null(dim(RnG.P.field)[1])){
  rng.cells.f<-as.vector(unlist(lapply(RnG.P.field[test.same],function(x) { as.numeric(unlist(strsplit(x,",")))  }  ) ))
}else{
 rng.cells.f <-rep(0,times=length(r.cells.f))
}}

if(!exists("nG.P.field")){ng.cells.f<-cells.f-g.cells.f}else{
if(!is.null(dim(nG.P.field)[1])){
  ng.cells.f<-as.vector(unlist(lapply(nG.P.field[test.same],function(x) { as.numeric(unlist(strsplit(x,",")))  }  ) ))
}else{
 ng.cells.f<-rep(0,times=length(r.cells.f))
}}
    #these values don't always exist

fields<-as.vector(unlist(lapply(cells.P.field[test.same],function(x) { length(unlist(strsplit(x,",")))  }  ) ))
trans.well<-as.vector(unlist(apply(cbind(fields,(green.c[test.same]/(green.c[test.same]+notGreen.c[test.same]))),1,function(x) {rep(x[2],times=x[1])})))

# extimate based on well average
percent.red.well<- as.vector(unlist(apply(cbind(fields,(red.c[test.same]/(all.c[test.same]-big.c[test.same]))),1,function(x) {rep(x[2],times=x[1])})))  ### WARNING maybe should use the not green estimate
prolif.well<- as.vector(unlist(apply(cbind(fields,(redAndNotGreen[test.same]/(notGreen.c[test.same]))),1,function(x) {rep(x[2],times=x[1])}))) #use this
#prolif.well<- as.vector(unlist(apply(cbind(fields,(redAndGreen[test.same]/(green.c[test.same]))),1,function(x) {rep(x[2],times=x[1])}))) #test

####### ann sorted row-wise so must fix:
ann.plate<-as.matrix(ann[(1+counter*well.type):(well.type+counter*well.type),"Symbol"])
ann.plate.label<-as.matrix(ann[(1+counter*well.type):(well.type+counter*well.type),"plate"])

dim(ann.plate)=dim(t(red.c))
ann.plate<-t(ann.plate)

dim(ann.plate.label)=dim(t(red.c))
ann.plate.label<-t(ann.plate.label)

ann.well<- as.vector(as.vector(unlist(apply(cbind(fields,ann.plate[test.same]),1,function(x) {rep(x[2],times=x[1])}))))
ann.plate.well<- as.vector(as.vector(unlist(apply(cbind(fields,ann.plate.label[test.same]),1,function(x) {rep(x[2],times=x[1])}))))

trans.p.field<-g.cells.f/cells.f

if(!is.null(dim(RG.P.field)[1])){
 prolif.p.field<-rng.cells.f/ng.cells.f
 percent.red.field<-r.cells.f/cells.f
}else{
 prolif.p.field<-r.cells.f/cells.f  # not strickly correct but will do for now
 percent.red.field<-r.cells.f/cells.f
}

row.labels<-row.labels[test.same]
col.labels<-col.labels[test.same]

well.labels<-as.vector(unlist(apply(cbind(fields,paste(row.labels,col.labels,sep="")),1,function(x) {rep(x[2],times=x[1])})))
field.labels <- as.vector(unlist(apply(as.matrix(fields),1,function(x) {c(1:x)})))

              
temporary<-t(rbind(cells.f,r.cells.f,g.cells.f,rg.cells.f,ng.cells.f,rng.cells.f,ann.plate.well,well.labels,field.labels,trans.p.field,trans.well,prolif.p.field,prolif.well,percent.red.field,percent.red.well,ann.well))
cell.counts<-rbind( cell.counts,temporary)


controls<-temporary[temporary[,"ann.well"]=="PLV101" | temporary[,"ann.well"]=="MOCK",]
control.prolif<-median(as.numeric(controls[,"prolif.well"]),na.rm=TRUE)

stand.dev.red.control<-sqrt(red.c*control.prolif*(1-control.prolif))    #sigma=sqrt(var)=sqrt( np(1-p) )
mean.red.control=(all.c-big.c)*control.prolif
prob.red.control<- (red.c-mean.red.control)/stand.dev.red.control




cells.f.mean<-unlist(lapply(cells.P.field,function(x) { mean(as.numeric(unlist(strsplit(x,","))) ) }  ) )
cells.f.sd<-unlist(lapply(cells.P.field,function(x) { sd(as.numeric(unlist(strsplit(x,","))) ) }  )   )
dim(cells.f.mean)<-dim(red.c)
dim(cells.f.sd)<-dim(red.c)


## core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S")


## vals<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","Description","Horfeome.position", "MGC.Accession","Pilot.Position","Gene.Name.pilotPos","Kiril.position")


##flex annotation prior to May 7th 2009
## vals<-cbind(round(as.numeric(t(100*percent.red)),n),round(as.numeric(t(100*green.c/(green.c+notGreen.c))),n),round(as.numeric(t(100*prolif.prob)),n),round(as.numeric(t(all.prob)),n),round(as.numeric(t(all.probLow)),n),round(as.numeric(t(all.probMid)),n),round(as.numeric(t(all.probHigh)),n),as.numeric(t(redAndGreen)),as.integer(t(mean)),as.integer(t(stand.dev)),as.numeric(t(redAndGreenLow)),as.integer(t(meanLow)),as.numeric(t(redAndGreenMid)),as.integer(t(meanMid)),as.numeric(t(redAndGreenHigh)), as.integer(t(meanHigh)),as.numeric(t(red.c)),as.numeric(t(redAndNotGreen)),as.numeric(t(green.c)), as.numeric(t(notGreen.c)),as.numeric(t(greenLow)),as.numeric(t(greenMid)), as.numeric(t(greenHigh)), round(as.numeric(t(DNA.G1andG2$adenG/DNA.G1andG2$adenNG)),n), round(as.numeric(t(DNA.aboveG1$adenG/DNA.aboveG1$adenNG)),n),  round(as.numeric(t(DNA.lt2$adenG/DNA.lt2$adenNG)),n),  round(as.numeric(t(DNA.gt4$adenG/DNA.gt4$adenNG)),n), round(as.numeric(t(DNA.S$adenG/DNA.S$adenNG)),n), round(as.numeric(t(DNA.G2/DNA.G1)),n),  round(as.numeric(t(DNA.G1)),n),round(as.numeric(t(DNA.fit.success$adenG)),n), round(as.numeric(t(DNA.fit.success$adenNG)),n), round(as.numeric(t(all.c-big.c)),n),round(as.numeric(t(score)),n),round(as.numeric(t(redAndGreen/green.c)),n), round(as.numeric(t(redAndNotGreen/notGreen.c)),n), as.numeric(t(notRed.c)),as.numeric(t(big.c)),as.numeric(t(all.c)),round(as.numeric(t(cells.f.mean)),n),round(as.numeric(t(cells.f.sd)),n),round(as.numeric(t(DNA.G1andG2$aden)),n),round(as.numeric(t(DNA.aboveG1$aden)),n),round(as.numeric(t(DNA.lt2$aden)),n), round(as.numeric(t(DNA.gt4$aden)),n), round(as.numeric(t(DNA.S$aden)),n) )

## core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","G1","G2","S","Gr-S+G2/G1","Gr-S+G2+ >4N/G1","Gr-<2N","Gr->4N","Gr-G1","Gr-G2","Gr-S","NtGr-S+G2/G1","NtGr-S+G2+ >4N/G1","NtGr-<2N","NtGr->4N","NtGr-G1","NtGr-G2","NtGr-S")

vals<-cbind(round(as.numeric(t(100*percent.red)),n),round(as.numeric(t(100*green.c/(green.c+notGreen.c))),n),round(as.numeric(t(100*prolif.prob)),n),round(as.numeric(t(all.prob)),n),round(as.numeric(t(all.probLow)),n),round(as.numeric(t(all.probMid)),n),round(as.numeric(t(all.probHigh)),n),as.numeric(t(redAndGreen)),as.integer(t(mean)),as.integer(t(stand.dev)),as.numeric(t(redAndGreenLow)),as.integer(t(meanLow)),as.numeric(t(redAndGreenMid)),as.integer(t(meanMid)),as.numeric(t(redAndGreenHigh)), as.integer(t(meanHigh)),as.numeric(t(red.c)),as.numeric(t(redAndNotGreen)),as.numeric(t(green.c)), as.numeric(t(notGreen.c)),as.numeric(t(greenLow)),as.numeric(t(greenMid)), as.numeric(t(greenHigh)), round(as.numeric(t(DNA.G1andG2$adenG/DNA.G1andG2$adenNG)),n), round(as.numeric(t(DNA.aboveG1$adenG/DNA.aboveG1$adenNG)),n),  round(as.numeric(t(DNA.lt2$adenG/DNA.lt2$adenNG)),n),  round(as.numeric(t(DNA.gt4$adenG/DNA.gt4$adenNG)),n), round(as.numeric(t(DNA.S$adenG/DNA.S$adenNG)),n), round(as.numeric(t(DNA.G2/DNA.G1)),n),  round(as.numeric(t(DNA.G1)),n),round(as.numeric(t(DNA.fit.success$adenG)),n), round(as.numeric(t(DNA.fit.success$adenNG)),n), round(as.numeric(t(all.c-big.c)),n),round(as.numeric(t(score)),n),round(as.numeric(t(redAndGreen/green.c)),n), round(as.numeric(t(redAndNotGreen/notGreen.c)),n), as.numeric(t(notRed.c)),as.numeric(t(big.c)),as.numeric(t(all.c)),round(as.numeric(t(cells.f.mean)),n),round(as.numeric(t(cells.f.sd)),n),round(as.numeric(t(DNA.G1andG2$aden)),n),round(as.numeric(t(DNA.aboveG1$aden)),n),round(as.numeric(t(DNA.lt2$aden)),n), round(as.numeric(t(DNA.gt4$aden)),n), round(as.numeric(t(DNA.inG1$aden)),n), round(as.numeric(t(DNA.inG2$aden)),n),round(as.numeric(t(DNA.S$aden)),n),round(as.numeric(t(DNA.G1andG2$adenG)),n),round(as.numeric(t(DNA.aboveG1$adenG)),n),round(as.numeric(t(DNA.lt2$adenG)),n), round(as.numeric(t(DNA.gt4$adenG)),n), round(as.numeric(t(DNA.inG1$adenG)),n), round(as.numeric(t(DNA.inG2$adenG)),n),round(as.numeric(t(DNA.S$adenG)),n),round(as.numeric(t(DNA.G1andG2$adenNG)),n),round(as.numeric(t(DNA.aboveG1$adenNG)),n),round(as.numeric(t(DNA.lt2$adenNG)),n), round(as.numeric(t(DNA.gt4$adenNG)),n), round(as.numeric(t(DNA.inG1$adenNG)),n), round(as.numeric(t(DNA.inG2$adenNG)),n),round(as.numeric(t(DNA.S$adenNG)),n)
  ) ## for 2010 flexible annotations

## vals<-cbind(round(as.numeric(t(100*percent.red)),n),round(as.numeric(t(100*green.c/(green.c+notGreen.c))),n),round(as.numeric(t(100*prolif.prob)),n),round(as.numeric(t(all.prob)),n),round(as.numeric(t(all.probLow)),n),round(as.numeric(t(all.probMid)),n),round(as.numeric(t(all.probHigh)),n),as.numeric(t(redAndGreen)),as.integer(t(mean)),as.integer(t(stand.dev)),as.numeric(t(redAndGreenLow)),as.integer(t(meanLow)),as.numeric(t(redAndGreenMid)),as.integer(t(meanMid)),as.numeric(t(redAndGreenHigh)), as.integer(t(meanHigh)),as.numeric(t(red.c)),as.numeric(t(redAndNotGreen)),as.numeric(t(green.c)), as.numeric(t(notGreen.c)),as.numeric(t(greenLow)),as.numeric(t(greenMid)), as.numeric(t(greenHigh)), round(as.numeric(t(DNA.G1andG2$adenG/DNA.G1andG2$adenNG)),n), round(as.numeric(t(DNA.aboveG1$adenG/DNA.aboveG1$adenNG)),n),  round(as.numeric(t(DNA.lt2$adenG/DNA.lt2$adenNG)),n),  round(as.numeric(t(DNA.gt4$adenG/DNA.gt4$adenNG)),n), round(as.numeric(t(DNA.S$adenG/DNA.S$adenNG)),n), round(as.numeric(t(DNA.G2/DNA.G1)),n),  round(as.numeric(t(DNA.G1)),n),round(as.numeric(t(DNA.fit.success$adenG)),n), round(as.numeric(t(DNA.fit.success$adenNG)),n), round(as.numeric(t(all.c-big.c)),n),round(as.numeric(t(score)),n),round(as.numeric(t(redAndGreen/green.c)),n), round(as.numeric(t(redAndNotGreen/notGreen.c)),n), as.numeric(t(notRed.c)),as.numeric(t(big.c)),as.numeric(t(all.c)),round(as.numeric(t(cells.f.mean)),n),round(as.numeric(t(cells.f.sd)),n)) # for Kirl and joseph redo

colnames(vals)<-core.vars

#vals.ann<-cbind(ann[(1+counter*96):(96+counter*96),1:4],vals,ann[(1+counter*96):(96+counter*96),5:8] ) #Hugo
#vals.ann<-cbind(ann[(1+counter*96):(96+counter*96),1:4],vals,ann[(1+counter*96):(96+counter*96),5:8] ) #Joesph
#vals.ann<-cbind(ann[(1+counter*96):(96+counter*96),1:4],vals,ann[(1+counter*96):(96+counter*96),5:10] ) #Kiril
core.ann<-c("plate", "row", "col", "Symbol")
col<-colnames(ann) %in% core.ann

vals.ann<-cbind(ann[(1+counter*well.type):(well.type+counter*well.type),col],vals,ann[(1+counter*well.type):(well.type+counter*well.type),!col]) #General
## vals.ann<-cbind(ann[(1+counter*well.type):(well.type+counter*well.type),core.ann],vals,ann[(1+counter*well.type):(well.type+counter*well.type),colnames(ann)[!place.core.ann]]) #General

counter<-counter+1
summ<-rbind(summ,vals.ann)
                          }

summ[1:5,]
#colnames(summ)<-sum.vars[1:ncol(summ)] # no longer needed
write.table(summ,well.output.file,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)




#################################################################################
######################################END SUMMARY GENERATION ##################
#################################################################################
#################################################################################
#################################################################################
#################################################################################


colnames(cell.counts)<-c("total","R","G","RG","nG","RnG","plate","well","field","trans.field","trans.well","prolif.field","prolif.well","percent.red.field","percent.red.well","symbol")
cell.counts[1:5,]

write.table(cell.counts,paste(gsub(".txt","",field.output.file),".RAW.txt",sep=""),row.names=FALSE,sep="\t")

cell.counts<-as.data.frame(cell.counts,stringsAsFactors = FALSE )
 ## num.cols<-c(1,2,3,4,6,7,8,9,10)
num.cols<-c("total","R","G","RG","nG","RnG","field","trans.field","trans.well","prolif.field","prolif.well","percent.red.field","percent.red.well")
 for(ic in 1:length(num.cols)){
 cell.counts[,num.cols[ic]]<-as.numeric(cell.counts[,num.cols[ic]]) }

cell.counts<-cell.counts[is.finite(cell.counts[,"total"]),]
stand.val.R<-sqrt(cell.counts[,"total"]*cell.counts[,"prolif.well"]*(1-cell.counts[,"prolif.well"]))    #sigma=sqrt(var)=sqrt( np(1-p) )
mean.val.R=cell.counts[,"total"]*cell.counts[,"prolif.well"]
prob.prolif<- (cell.counts[,"R"]-mean.val.R)/stand.val.R
stand.val.G<-sqrt(cell.counts[,"total"]*cell.counts[,"trans.well"]*(1-cell.counts[,"trans.well"]))    #sigma=sqrt(var)=sqrt( np(1-p) )
mean.val.G=cell.counts[,"total"]*cell.counts[,"trans.well"]
prob.trans<- (cell.counts[,"G"]-mean.val.G)/stand.val.G
cells<-cbind(cell.counts,prob.prolif,mean.val.R,stand.val.R,prob.trans,mean.val.G,stand.val.G)
colnames(cells)<-c(colnames(cell.counts), "prolif z-score","prolif-expected","prolif-sd","transducted z-score","transduced-expected","transduced-sd")
cells[1:5,]
write.table(cells,field.output.file,row.names=FALSE,sep="\t")



###################################################################################################
###################################################################################################
################################# restart from written files
#######################
summ<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
cells<-read.delim(field.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
col.array<-rainbow(length(unique(cells[,"plate"])))
colnames(cells)<-c("total","R","G","plate","well","field","transduced-field","transduced-well","prolif-field","prolif-well","symbol","prolif z-score","prolif-expected","prolif-sd","transducted z-score","transduced-expected","transduced-sd")

colnames(cells)<-c("total","R","G","plate","well","field","transduced-field","transduced-well","prolif-field","prolif-well","percent.red.field","percent.red.well","symbol", "prolif z-score","prolif-expected","prolif-sd","transducted z-score","transduced-expected","transduced-sd")
colnames(summ)<-sum.vars # sum.vars found at top of file

the.screen<-"Joseph"
the.screen<-"Kiril"
the.screen<-"Hugo"
the.screen<-"Ca siRNA"

sort(tapply(cells[,"symbol"],cells[,"symbol"],length)) ## look for controls
###################################################################################################
###################################################################################################
################## transduction AND proliferation index ALL PLATES:
the.plates<-unique(cells[,"plate"])
col.array<-rainbow(length(the.plates))

color = colorRampPalette(c("blue","green","red"),space="Lab",interpolate="linear")
col.array<-color(length(the.plates))

names(col.array)<-the.plates
## names(col.array)<-unlist(lapply(strsplit(files,split="_"),function(x) x[3])) # if plate ids in 3rd column
y.range.low<- range(cells[is.finite(cells[,"prolif z-score"]),"prolif z-score"])[1]
y.range.high<- range(cells[is.finite(cells[,"prolif z-score"]),"prolif z-score"])[2]
y.range.high<-20  

 sample<-rep(TRUE,times=dim(cells)[1])
  plot(cells[sample,"total"],cells[sample,"transducted z-score"],col=col.array[cells[sample,"plate"]],pch=19,cex=1.0,ylim=c(y.range.low,y.range.high),xlab="Cells per field",ylab="Transduction z-score (bases on well estimates)",font.lab=2,cex.lab=1.5,main=paste(the.screen,"Screen :  (PLV101-square)",sep=" "))
  abline(h=0)
#sample<-cells[,"symbol"]=="MOCK"
#points(cells[sample,"total"],cells[sample,"transducted z-score"],pch=23,cex=1.6,col="black",lwd=2)
sort(tapply(cells[sample,"symbol"],cells[sample,"symbol"],length))
the.control<-"HSPCN111-N"

the.control<-"No Lipid"
the.control2<-"PLV411"

the.control<-"plv411"
the.control2<-"MOCK"

sample<-cells[,"symbol"]==the.control

points(cells[sample,"total"],cells[sample,"prolif z-score"],pch=22,cex=1.2,col="black",lwd=1)

sample<-cells[,"symbol"]==the.control2
points(cells[sample,"total"],cells[sample,"prolif z-score"],pch=23,cex=1.2,col="red",lwd=1)
legend(y.range.low,y.range.high,legend=c(names(col.array),the.control,the.control2),col=c(col.array,"black","red"),pch=c(rep(19,times=length(col.array)),22,23),cex=1.2)

savePlot(paste(the.screen," transduction z _F.jpeg",sep=""),type="jpeg")

sample<-rep(TRUE,times=dim(cells)[1])
## sample<-cells[,"plate"]<=3
## sample<-cells[,"plate"]>3 & cells[,"plate"]<=9

  plot(cells[sample,"total"],cells[sample,"prolif z-score"],col=col.array[cells[sample,"plate"]],pch=19,cex=0.85,ylim=c(y.range.low,y.range.high),xlab="Cells per field",ylab="Proliferative z-score (bases on well estimates)",font.lab=2,cex.lab=1.5,main=paste(the.screen,"Screen : (MOCK-diamond) (control-square)",sep=" ") )
  abline(h=0)


sample<-cells[,"symbol"]==the.control

points(cells[sample,"total"],cells[sample,"prolif z-score"],pch=22,cex=1.2,col="black",lwd=1)

sample<-cells[,"symbol"]==the.control2
points(cells[sample,"total"],cells[sample,"prolif z-score"],pch=23,cex=1.2,col="red",lwd=1)

# legend(y.range.low,y.range.high,legend=names(col.array),col=col.array,pch=19,cex=1.5)
#legend(1000,y.range.high,legend=names(col.array),col=col.array,pch=19,cex=1.0)
sample<-cells[,"plate"]<=3
& cells[,"plate"]<=6
points(cells[sample,"total"],cells[sample,"prolif z-score"],pch=21,cex=1.2,col="black",lwd=1)

legend(y.range.low,y.range.high,legend=c(names(col.array),the.control,the.control2),col=c(col.array,"black","red"),pch=c(rep(19,times=length(col.array)),22,23),cex=1.2)

savePlot(paste(the.screen," proliferation z _F.jpeg",sep=""),type="jpeg")


###########################################
# for Ca screen only
plate.series<-3
per.r<-cbind(plate.series,summ[summ[,"plate"]==paste(plate.series,1,sep="."),c("row","col","Symbol","%R")],summ[summ[,"plate"]==paste(plate.series,2,sep="."),c("%R")],summ[summ[,"plate"]==paste(plate.series,3,sep="."),c("%R")],NA,NA)
colnames(per.r)[c(1,5:9)]<-c("plate series","Rep:1 %R","Rep:2 %R","Rep:3 %R","Mean %R","SD %R")
per.r[,"Mean %R"]<-apply(per.r[,c("Rep:1 %R","Rep:2 %R","Rep:3 %R")],1,function(x) mean(x,na.rm=TRUE))
per.r[,"SD %R"]<-apply(per.r[,c("Rep:1 %R","Rep:2 %R","Rep:3 %R")],1,function(x) sd(x,na.rm=TRUE))
write.table(per.r,paste("plates",plate.series,"NORMALIZED percentR.txt",sep="_"),row.names=FALSE,sep="\t")



###################################################################################################
###################################################################################################
################################# NORMALIZE  ##############################################
###################################################################################################
############################################### Ca screen
setwd( "/media/Bioinform-D/Research/Cellomics/Ca screen")
sum.vars<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S","Entrez_ID","Cell_loc","Description","Gene.ID","Accession.Number","GI.Number","Library.Type")
well.output.file<-"Ca_summary_NOTGREEN_DNA.txt"
well.type<-384
row.type<-16
col.type<-24
the.screen<-"Ca siRNA"
plates.all<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized
colnames(plates.all)<-sum.vars
#######################################################
###############################################################
#####################################################################################
###START#########################
 setwd("/media/Bioinform-D/Research/Cellomics/Hugo screen")
sum.vars<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","Description","GenBank_Accn","dbEST_Id" ,"MGC.Accession")
well.output.file<-"Hugo_summary_NOTGREEN_DNA.txt"
the.screen="Hugo"
plates.all<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized
colnames(plates.all)<-sum.vars
############################################
############################################### Leo  ### latest DNA and flexible annotation
setwd( "/media/Bioinform-D/Research/Cellomics/Leo screen")
load("Leo_ann.RData") # neded for annotation file 
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann
sum.vars<-c(core.ann,"Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S",colnames(ann)[!place.core.ann])
well.output.file<-"Leo_summary_NOTGREEN_DNA.txt"
the.screen<-"Leo"
plates.all<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized
colnames(plates.all)<-sum.vars
#######################################################
###########################################################
setwd( "/media/Bioinform-D/Research/Cellomics/Kiril screen")
sum.vars<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","Description","Horfeome.position", "MGC.Accession","Pilot.Position","Gene.Name.pilotPos","Kiril.position")
well.output.file<-"Kiril_summary_NOTGREEN_DNA.txt"
the.screen="Kiril"
plates.all<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized
colnames(plates.all)<-sum.vars

###START#########################
 setwd("/media/Bioinform-D/Research/Cellomics/Joseph screen")
sum.vars<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","Description","GenBank_Accn","dbEST_Id" ,"MGC.Accession")
well.output.file<-"Joseph_summary_NOTGREEN_DNA.txt"
the.screen="Joseph"
plates.all<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized
colnames(plates.all)<-sum.vars
############################################

###START#########################
 setwd("/media/Bioinform-D/Research/Cellomics/Hugo screen")
sum.vars<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","Description","GenBank_Accn","dbEST_Id" ,"MGC.Accession")
well.output.file<-"Hugo_summary_NOTGREEN_DNA.txt"
the.screen="Hugo"
plates.all<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized
colnames(plates.all)<-sum.vars
############################################


library(robust)
library(robustbase)
library(nls2)
library(akima)

############# read in the reults one screen or many screens #############


colnames(plates.all)<-gsub("%","percent",colnames(plates.all))
colnames(plates.all)<-gsub("/","_",colnames(plates.all))

test.scores<-c("Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","percentR" )
zero.center<-c( TRUE,           TRUE,         TRUE,          TRUE,         FALSE)

test.scores<-c("S+G2_G1","S","percentR" )
zero.center<-c( FALSE,FALSE,         FALSE)

##########################################################################
################### NORMALIZATION LOOP
the.plates<-unique(plates.all[,"plate"])
for(itest in 1:length(test.scores)){
the.score<-test.scores[itest]
data.in<- (plates.all[is.finite(plates.all[,the.score]),the.score])

jpeg( paste(the.screen,"_raw_","S+G1_G2",".jpeg",sep=""), width = 1280, height = 1024, units = "px",quality = 100 )

ahist<-hist(data.in,col="pink",breaks=100,main=paste(the.screen," UnNormalised for ",the.score,sep=""),xlab=the.score,cex.lab=1.5,cex.main=2.0)
aden<-density(data.in,n=1024,adjust=1)
apeak<-max(ahist$counts)/max(aden$y)

the.test<-shapiro.test(data.in)
text(ahist$breaks[floor((length(ahist$breaks)/2)+5)],max(ahist$counts)*0.9, paste("Shapiro Test:",signif(the.test$p.value,3),sep=" "),cex=1.5)
lines(aden$x,apeak*aden$y,col="green",lwd=2)
dev.off()
#savePlot(paste(the.screen,the.score,".jpeg",sep=""),type="jpeg")


the.plates<-unique(plates.all[,"plate"])
the.medians<-tapply(plates.all[,the.score],plates.all[,"plate"], function(x) median(x,na.rm=TRUE))
the.fit.medians<-the.medians



###################### filter the results
filter<-is.finite(plates.all[,the.score])

RG.expected.min<-15 # typically 50 set to 15 for Leo screen
cut.low<- -20 # for Not Green well estimates
cut.high<- 10 # for  NOt Green  well estimates

########################################### SET positive controls - exclude from normalization  ############################
pos.control<-plates.all[,"Symbol"]=="CCNE1" | plates.all[,"Symbol"]=="101-CCNE1" | plates.all[,"Symbol"]=="411-CCNE1" | plates.all[,"Symbol"]=="101-CYCD2"
#pos.control<-plates.all[,"Symbol"]=="GFP" | plates.all[,"Symbol"]=="PLK1" | plates.all[,"Symbol"]=="ATP2C1" | plates.all[,"Symbol"]=="empty"
#################################################################################################################

#GFP no good at all
#PLK1 should kill all but is variable
#ATP2C1 increase proliferation ?
#empty

########################################### exclude low R and G weels from normalizationn  #####################
RG.small<-plates.all[,"RG_expected"] < RG.expected.min & !is.na(plates.all[,"RG_expected"])
## RG.small<-plates.all[,"RG_High"] <20
## low.val<-plates.all[,the.score] < cut.low 
## high.val<- plates.all[,the.score] > cut.high 
## filter<-filter & !RG.small & !low.val & !high.val & !pos.control
filter<-filter & !RG.small & !pos.control
################################################################################################################


########################################### SET the control wells for normalization  ############################
filter.plv<-(plates.all[,"Symbol"]=="PLV101" | plates.all[,"Symbol"]=="plv101") & !RG.small # ARVEC screen
##filter.plv<-(plates.all[,"Symbol"]=="NT") &  !RG.small # CA screenor use NT or ATP2C1 use NT
#################################################################################################################

the.plv.medians<-tapply(plates.all[filter.plv,the.score],plates.all[filter.plv,"plate"], function(x) median(x,na.rm=TRUE))
the.filtered.medians<-tapply(plates.all[filter,the.score],plates.all[filter,"plate"], function(x) median(x,na.rm=TRUE))
#the.filtered.medians<-tapply(plates.all[filter & !low.val & !high.val,the.score],plates.all[filter & !low.val & !high.val,"plate"], function(x) median(x,na.rm=TRUE))

number<-tapply(plates.all[filter,"plate"],plates.all[filter,"plate"],length)
the.number.plv.medians<-tapply(plates.all[filter.plv,the.score],plates.all[filter.plv,"plate"], length)

jpeg( paste(the.screen,"_BOXPLOT_PLATE_MEDIANS_",the.score,".jpeg",sep=""), width = 1280, height = 1024, units = "px",quality = 100 )
#to.box<-plates.all[filter.plv,c("plate",the.score)]
to.box<-plates.all[filter,c("plate",the.score)]
colnames(to.box)<-c("x","y")
#the.boxplot<-boxplot(y~x,data=to.box,ylab=the.score,xlab="Plate",main="For PLV101s with RG >50")
the.boxplot<-boxplot(y~x,data=to.box,ylab=the.score,xlab="Plate",main="Boxplot for PLATE MEDIANS",cex.lab=1.5,cex.axis=1.5,cex.main=2,col=rainbow(length(the.plates)))
text(1,max(to.box$y,na.rm=TRUE),labels=toString(paste(the.boxplot$names,"(",the.boxplot$n,")",sep="")),pos=4,cex=1.0,cex.lab=1.5,font=1.5,col="red")
dev.off()
#savePlot( paste(the.screen,"_BOXPLOT_PLATE_MEDIANS_",the.score,".jpeg",sep=""),type="jpeg")

jpeg( paste(the.screen,"_BOXPLOT_CONTROL_MEDIANS_",the.score,".jpeg",sep=""), width = 1280, height = 1024, units = "px",quality = 100 )
#to.box<-plates.all[filter.plv,c("plate",the.score)]
to.box<-plates.all[filter.plv,c("plate",the.score)]
colnames(to.box)<-c("x","y")
the.boxplot<-boxplot(y~x,data=to.box,ylab=the.score,xlab="Plate",main=paste("For PLV101s with RG > ",RG.expected.min,sep=""),col=rainbow(length(the.boxplot$names)))
#the.boxplot<-boxplot(y~x,data=to.box,ylab=the.score,xlab="Plate",main="Boxplot for PLATE CONTROLS :ATP2C1",cex.lab=1.5,cex.axis=1.5,cex.main=2,col=rainbow(length(the.plates)))
text(1,max(to.box$y),labels=toString(paste(the.boxplot$names,"(",the.boxplot$n,")",sep="")),pos=4,cex=1.0,cex.lab=1.5,font=1.5,col="red")
dev.off()
#savePlot(paste(the.screen,"_BOXPLOT_CONTROL_MEDIANS_",the.score,".jpeg",sep=""),type="jpeg")

################### if the plv.means runs out for the plate use the filtered medians
if(length(the.filtered.medians) != length(the.plates)){print("ERROR NORMALIZATION FAILURE - no genes on plate")}
the.plv.medians<-the.plv.medians[as.character(the.plates)]
the.filtered.medians<-the.filtered.medians[as.character(the.plates)]
missing<-is.na(the.plv.medians)
the.plv.medians[missing]<-the.filtered.medians[missing]
names(the.plv.medians)[missing]<-names(the.filtered.medians)[missing]

#x11( width = 11, height = 11)
jpeg( paste(the.screen,"_raw_ALL_WELLS_",the.score,".jpeg",sep=""), width = 1280, height = 1024, units = "px",quality = 100 )

nf<-layout(matrix(c(1:((ceiling(length(the.plates)/3))*3)),(ceiling(length(the.plates)/3)),3,byrow=TRUE),heights=c(rep(1,times=ceiling(length(the.plates)/3)) ),widths=c(rep(1,times=3) ) )
layout.show(nf)
par(mar=c(3,3.5,2.1,2.1),mgp=c(2,1,0)) #c(bottom, left, top, right)

for(i in 1:length(the.plates)){
data.plate<-(plates.all[plates.all[,"plate"]==the.plates[i] & filter, the.score])
ahist.plate<-hist(data.plate,col=col.array[i],breaks=100,main=paste("Plate=",the.plates[i],sep=""),xlab=the.score,cex.main=1.5,cex.lab=1.5,cex.axis=2)
data.plate<-data.plate[!is.na(data.plate)]
aden<-density(data.plate,n=1024,adjust=1)
apeak<-max(ahist.plate$counts)/max(aden$y)
lines(aden$x,apeak*aden$y,col="green",lwd=2)
the.fit.medians[names(the.fit.medians)==the.plates[i]]<-coef(ltsreg(data.plate~1 ))
#the.fit.medians1<-the.fit.medians[the.plates[i]]
 #the.model<- ltsreg(data ~ plate , data=data.plate)


#abline(v=the.medians[names(the.medians)==the.plates[i]],col="blue",lwd=3)
#abline(v=the.fit.medians[names(the.fit.medians)==the.plates[i]],col="orange2",lwd=3)
abline(v=the.filtered.medians[names(the.filtered.medians)==the.plates[i]],col="black",lwd=4)
abline(v=the.plv.medians[names(the.plv.medians)==the.plates[i]],col="cyan",lwd=3,lty=2)
text(max(ahist.plate$mid),max(ahist.plate$counts)-0.5,paste("median= ",round(the.filtered.medians[names(the.filtered.medians)==the.plates[i]],2),"(",number[names(number)==the.plates[i]],")",sep=""),cex=1.5,pos=2,font=1.5)
}

#savePlot(paste(the.screen,"_raw_ALL_WELLS_",the.score,".jpeg",sep=""),type="jpeg")
dev.off()

#x11()
jpeg(paste(the.screen,"_NORMALIZED_",the.score,".jpeg",sep=""), width = 1280, height = 1024, units = "px",quality = 100 )

############### meadian normalize data
#centers.all<-the.fit.medians
centers.all<-the.filtered.medians
#centers.all<-the.plv.medians
for(i in 1:length(the.plates)){
plates.all[plates.all[,"plate"]==the.plates[i],the.score] <-plates.all[plates.all[,"plate"]==the.plates[i],the.score]-
  centers.all[names(centers.all)==the.plates[i]]
}

#### modification so don't have neagiatve values
if(!zero.center[itest]){
  min.value<-abs((min(plates.all[,the.score],na.rm=TRUE)))
 plates.all[,the.score]<- plates.all[,the.score]+min.value}

#sum(is.na(plates.all[,the.score]))


data.in<- (plates.all[filter,the.score])

#3sum(is.na(data.in))
#sum(is.na(filter))
#data.in<-plates.all.ori[filter,the.score]-coef(ltsreg(plates.all.ori[filter,the.score]~1 ))

ahist<-hist(data.in,col="pink",breaks=100,main=paste(the.screen," Normalised for ",the.score,sep=""),xlab=the.score,cex.lab=1.5,cex.main=2.0)
aden<-density(data.in,n=1024,adjust=1)
apeak<-max(ahist$counts)/max(aden$y)
#apeak<-30/max(aden$y)
lines(aden$x,apeak*aden$y,col="green",lwd=2)

the.mean<-median(data.in,na.rm=TRUE)
the.sd<-sd(data.in,na.rm=TRUE)
####################### fit normal to data
 to.fit<-data.frame(y=aden$y*apeak, x=aden$x)
A <- max(to.fit$y)/dnorm(the.mean,mean=the.mean,sd=the.sd, log=FALSE)
  the.fit<-nls(y~A*dnorm(x,mean=the.mean,sd=the.sd, log=FALSE)
               ,data=to.fit
               ,start=list(A=A,the.mean=the.mean,the.sd=the.sd)
                ,trace=TRUE
               ,control=list(maxiter=1000, minFactor=1/4048,tol=1e-4, warnOnly = TRUE))

coef(the.fit)
curve(coef(the.fit)["A"]*dnorm(x,mean=coef(the.fit)["the.mean"],sd=coef(the.fit)["the.sd"], log=FALSE),add=TRUE, col="red",lwd=2,lty="dashed")

wings<-(data.in > (coef(the.fit)["the.mean"]+coef(the.fit)["the.sd"]*3) ) | (data.in < (coef(the.fit)["the.mean"]-coef(the.fit)["the.sd"]*3))
the.test<-shapiro.test(data.in[!wings])
text(coef(the.fit)["the.mean"]+coef(the.fit)["the.sd"]*2,max(to.fit$y),paste("Shapiro Test:",signif(the.test$p.value,3),sep=" "),cex=1.5,pos=4)
if(!zero.center[itest]){
text(coef(the.fit)["the.mean"]+coef(the.fit)["the.sd"]*2,max(to.fit$y),paste("Distribution centered to",min.value,sep=" "),cex=1.5,pos=1)}
dev.off()

}  ### end loop over normalization

#normalized<-plates.all[,c(colnames(plates.all)[1:4],test.scores,colnames(plates.all)[27:36],colnames(plates.all)[45:dim(plates.all)[2]])]
normalized<-plates.all
#write.table(normalized,paste(the.screen,"NORMALIZED","txt",sep="."),row.names=FALSE,sep="\t")
write.table(normalized,paste(the.screen,"NORMALIZED","txt",sep="."),row.names=FALSE,sep="\t")




#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
################################# restart from written files  Q_Q PLOTS
#######################
plates.all<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized
colnames(plates.all)<-sum.vars

normalized<-read.delim(paste(the.screen,"NORMALIZED","txt",sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
colnames(normalized)<-sum.vars
confirmed<-read.delim("confirmed.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
colnames(confirmed)[8:dim(confirmed)[2]]<-sum.vars[5:length(sum.vars)]

#######################

plates.all<-normalized
plates.all<-joseph
plates.all<-joseph.norm
a.score<-1:4
names(a.score)<-c("Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green")
a.cut<-c(9.0,8.9,8.0,8.0)
names(a.cut)<-c("Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green")

par(mar=c(4.5,5.5,2.5,1.1),mgp=c(2.8,1,0))  #c(bottom, left, top, right) use 4 for last one
nf<-layout(matrix(c(1,2,3,4),4,byrow=TRUE),heights=c(1,1,1,1))
nf<-layout(matrix(c(1,2,3,4),2,2,byrow=TRUE),heights=c(1,1,1,1))
layout.show(nf)

the.score<-"Z-Low Green"
the.score<-"Z-Mid Green"
the.score<-"Z-High Green" 
the.score<-"Z-all Green"
the.score<-"percentR"
#par(mar=c(3,5.5,2.1,5.1),mgp=c(2.5,1,0))

labels<-as.character(paste(plates.all[is.finite(plates.all[,the.score]),"plate"],plates.all[is.finite(plates.all[,the.score]),"row"],plates.all[is.finite(plates.all[,the.score]),"col"],"-",plates.all[is.finite(plates.all[,the.score]),"Symbol"],sep=""))
xlabels<-as.character(paste(plates.all[is.finite(plates.all[,the.score]),"plate"],plates.all[is.finite(plates.all[,the.score]),"row"],plates.all[is.finite(plates.all[,the.score]),"col"],sep=""))
symbols<-plates.all[is.finite(plates.all[,the.score]),"Symbol"]
RGs<-plates.all[is.finite(plates.all[,the.score]),"RG_expected"]
data.in<- plates.all[is.finite(plates.all[,the.score]),the.score]

the.mean<-mean(data.in)
#the.mean<-min.value

#my.qq.plot(data.in,distribution="norm",col="blue",xlab="Expected Score",ylab="Observed score",xlim=c(-15+the.mean,15+the.mean), ylim=c(-15+the.mean,20+the.mean),main=paste("Normalised ",the.screen,the.score,"  With 95% confidence intervals",sep=" "),the.mean=the.mean,the.sd=sd(data.in),cex.lab=1.5,cex.axis=1.5,cex.main=1.5)

#my.qq.plot(data.in,distribution="norm",col="blue",xlab="Expected Score",ylab="Observed score",xlim=c(5+the.mean,15+the.mean), ylim=c(0+the.mean,20+the.mean),main=paste("Normalised ",the.screen,the.score,"  With 95% confidence intervals",sep=" "),the.mean=the.mean,the.sd=sd(data.in),cex.lab=1.5,cex.axis=1.5,cex.main=1.5)

my.qq.plot(data.in,distribution="norm",col="blue",xlab="Expected Score",ylab="Observed score",xlim=c(6,20), ylim=c(6,35),main=paste("Normalised Screen A with 95% confidence intervals",":",the.score,sep=" "),the.mean=the.mean,the.sd=sd(data.in),cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex=1.5)

my.qq.plot(data.in,distribution="norm",col="blue",xlab="Expected Score",ylab="Observed score",xlim=c(-25+the.mean,35+the.mean), ylim=c(-15+the.mean,35+the.mean),main=paste("Normalised  with 95% confidence intervals",":",the.score,sep=" "),the.mean=the.mean,the.sd=sd(data.in),cex.lab=1.5,cex.axis=1.5,cex.main=2.0,cex=1.5)

#my.qq.plot(data.in,distribution="norm",col="blue",xlab="",ylab="Observed score",xlim=c(5,14.5), ylim=c(5,20),main=paste(the.screen,":",the.score,sep=" "),the.mean=the.mean,the.sd=sd(data.in),cex.lab=1.5,cex.axis=1.50,cex.main=1.5,cex=2.0)

qq<- qq.data(data.in,distribution="norm",the.mean=the.mean,the.sd=sd(data.in),plot.it=FALSE)
#posns<-grep("NT",labels[qq$ord]) # CA SCREEN
posns<-grep("PLV101",labels[qq$ord])
lowRG<-RGs[qq$ord][posns]<RG.expected.min # set during normalization typicall7 50
points(qq$x[posns[!lowRG]],qq$y[posns[!lowRG]],col="cyan",pch=19,cex=1.5)
#points(qq$x[posns[lowRG]],qq$y[posns[lowRG]],col="cyan",pch=21,lwd=2)

#posns<-c(1:length(symbols))[symbols[qq$ord] %in% c("CCNE1","101-CCNE1","411-CCNE1","101-CYCD2")]
posns<-c(1:length(symbols))[symbols[qq$ord] %in% c("CCNE1","101-CCNE1","411-CCNE1","101-CYCD2","1H7-CyclD2")]
#posns<-c(1:length(symbols))[symbols[qq$ord] %in% c("PLK1")] # CA SCREEN
lowRG<-RGs[qq$ord][posns]<29
points(qq$x[posns[!lowRG]],qq$y[posns[!lowRG]],col="red",pch=19,cex=1.5)
pos.control<-plates.all[,"Symbol"]=="CCNE1" | plates.all[,"Symbol"]=="101-CCNE1" | plates.all[,"Symbol"]=="411-CCNE1" | plates.all[,"Symbol"]=="101-CYCD2"
#posns<-c(1:length(symbols))[symbols[qq$ord] %in% c("ATP2C1")] ###### CA SCREEN
lowRG<-RGs[qq$ord][posns]<29
points(qq$x[posns[!lowRG]],qq$y[posns[!lowRG]],col="green",pch=19,cex=1.5)
#GFP no good at all
#PLK1 should kill all but is variable
#ATP2C1 increase proliferation ?
#empty



subset<-conf.screen.norm[conf.norm.hits.posns[[a.score[the.score]]],c("plate","Validated","row","col",the.score)]
subset<-subset[subset[,"Validated"]=="N",c("plate","row","col",the.score)]
subset.xlabels<-as.character(paste(subset[,"plate"],subset[,"row"],subset[,"col"],sep=""))
posns<-c(1:length(xlabels))[xlabels[qq$ord] %in%  subset.xlabels]
points(qq$x[posns],qq$y[posns],col="green",pch=19,cex=1.5)
abline(h=a.cut[the.score],col="forestgreen",lwd=2)

subset<-conf.screen.norm[conf.norm.hits.posns[[a.score[the.score]]],c("plate","Validated","row","col",the.score)]
subset<-subset[subset[,"plate"]<=15,c("plate","Validated","row","col",the.score)]
subset<-subset[subset[,"Validated"]=="Y",c("plate","row","col",the.score)]
subset.xlabels<-as.character(paste(subset[,"plate"],subset[,"row"],subset[,"col"],sep=""))
posns<-c(1:length(xlabels))[xlabels[qq$ord] %in%  subset.xlabels]
points(qq$x[posns],qq$y[posns],col="forestgreen",pch=19,cex=1.5)




savePlot(filename=paste("Q-Q -joseph with confirmation2",the.score,"jpeg",sep="."),type="jpeg")
savePlot(filename=paste("Q-Q -joseph with confirmation2",the.score,"tiff",sep="."),type="tiff")

##### Ca screen
legend(20,80,legend=c("NT","PLK1","ATP2C1"),col=c("cyan","red","green"),pch=20,cex=2)
text(60,20,labels="Mean:50, Stand. Dev:10",cex=1.5)

legend(4.75,22,legend=c("PLV101","Positive Control","Confirmed","Not Confirmed"),col=c("cyan","red","forestgreen","green"),pch=19,cex=2.0,bty="n")

legend(6,16,legend=c("PLV101","Positive Control"),col=c("cyan","red"),pch=19,cex=2.0)
legend(-10,20,legend=c("PLV101","Positive Control"),col=c("cyan","red"),pch=19,cex=2.0)
abline(h=8.9,col="forestgreen",lwd=2)
text(-1,9.5,labels=c("Cut at Z-Low=8.9"),pos=4,cex=1.5,col="forestgreen")

############################# DRAW levels

selected.data<-identify(qq$x,qq$y,labels=as.character(round(data.in[qq$ord],3)),col="forestgreen",font=2,cex=1.5,offset=1,atpen='TRUE',plot=TRUE)
savePlot(filename=paste(the.screen,"Normalized Q-Q Plot for",the.score,"jpeg",sep="."),type="jpeg")

savePlot(filename=paste("Q-Q -joseph LOW plv1o1 and controls",the.score,"jpeg",sep="."),type="jpeg")
savePlot(filename=paste("Q-Q -joseph LOW plv1o1 and controls ZOOM2",the.score,"jpeg",sep="."),type="jpeg")

selected.data<-identify(qq$x,qq$y,labels=labels[qq$ord],col="red",cex=1,atpen='TRUE')



## subset<- kiril.norm[kiril.norm.hits.posns$z.low,c("plate","row","col","Z-Low Green")]
## subset<-subset[subset[,"plate"]<=15,c("plate","row","col","Z-Low Green")]
## subset<-subset[subset[,"Z-Low Green"]>10,c("plate","row","col","Z-Low Green")]
## subset.xlabels<-as.character(paste(subset[,"plate"],subset[,"row"],subset[,"col"],sep=""))
## posns<-c(1:length(xlabels))[xlabels[qq$ord] %in%  subset.xlabels]
## points(qq$x,qq$y,col="magenta",pch=21)

cuts<-c(0.9,0.7,0.7,0.7,27) # q-q cuts for JOseph UNnormalized
cuts<-c(9.0,8.9,8.0,8.0,27) # q-q cuts for Joseph Normalized

cuts<-c(6.0,6.7,7.5,13,21) # q-q cuts for Kiril UNnormalized
cuts<-c(10.0,10.0,12.0,12.0,27) # q-q cuts for Kiril Normalized



###############################################################
#####################################################################################
###########################################################
setwd( "/media/Bioinform-D/Research/Cellomics/Kiril screen")
sum.vars<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","Description","Horfeome.position", "MGC.Accession","Pilot.Position","Gene.Name.pilotPos","Kiril.position")
well.output.file<-"Kiril_summary_NOTGREEN_DNA.txt"
kiril<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized
the.screen="Kiril"
kiril.norm<-read.delim(paste(the.screen,"NORMALIZED","txt",sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
kiril.norm.plv<-read.delim(paste(the.screen,"NORMALIZED_BY_PLV","txt",sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
colnames(kiril.norm.plv)<-sum.vars
colnames(kiril.norm)<-sum.vars
colnames(kiril)<-sum.vars

###START#########################
 setwd("/media/Bioinform-D/Research/Cellomics/Joseph screen")
sum.vars<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","Description","GenBank_Accn","dbEST_Id" ,"MGC.Accession")
well.output.file<-"Joseph_summary_NOTGREEN_DNA.txt"
joseph<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized
the.screen="Joseph"
joseph.norm<-read.delim(paste(the.screen,"NORMALIZED","txt",sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
colnames(joseph.norm)<-sum.vars
colnames(joseph)<-sum.vars
############################################

###START#########################
 setwd("/media/Bioinform-D/Research/Cellomics/Hugo screen")
sum.vars<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","Description","GenBank_Accn","dbEST_Id" ,"MGC.Accession")
well.output.file<-"Hugo_summary_NOTGREEN_DNA.txt"
hugo<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized
the.screen="Hugo"
hugo.norm<-read.delim(paste(the.screen,"NORMALIZED","txt",sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
colnames(hugo.norm)<-sum.vars
colnames(hugo)<-sum.vars
############################################




RG.label<-c("RG_expected","RG_Low_expected","RG_Mid_expected","RG_High_expected")
RG.cut<-c(30,25,25,25)
            

screen2<-joseph.norm
screen2<-confirmed
screen2.cuts<-c(9.0,8.9,8.0,8.0) # q-q cuts for Joseph Normalized
screen1<-kiril.norm
screen1.cuts<-c(9.0,8.9,8.0,8.0) #q-q cuts for Kiril Normalized
screen1.cuts<-c(10,10,10,10) #q-q cuts for Kiril Normalized and by PLV ok

screen2<-joseph
screen2.cuts<-c(0.9,0.7,0.6,1.0) # q-q cuts for JOseph UNnormalized
screen1<-kiril
screen1.cuts<-c(5.0,6.0,6.5,9) # q-q cuts for Kiril UNnormalized

screen2<-hugo.norm
screen2.cuts<-c(4.4,4.4,4.3,4.0) # q-q cuts for Hugo normalized
screen1<-joseph.norm
screen1.cuts<-c(9.0,8.9,8.0,8.0) # q-q cuts for Joseph Normalized        


common<-intersect(screen2[,"Symbol"],screen1[,"Symbol"])
length(common)

##################### subset
posns<-match(screen2[,"Symbol"],common) ### posns of  joseph in common length of common
#posns<-posns[!is.na(posns)]
screen2<-screen2[!is.na(posns),]
dim(screen2)
########################
##################### subset
posns<-match(screen1[,"Symbol"],common) ### length common
#posns<-posns[!is.na(posns)]
#screen1<-screen1[posns,]
screen1<-screen1[!is.na(posns),]
dim(screen1)
########################
##  dim(screen2)
## [1] 1728   48
## [1] 2784   50

length(common)






symbols1<-{}
symbols2<-{}
screen1.f.counts<-list(z.all=NA,z.low=NA,z.mid=NA,z.high=NA)
screen2.f.counts<-list(z.all=NA,z.low=NA,z.mid=NA,z.high=NA)
screen1.f.posns<-list(z.all=c(),z.low=c(),z.mid=c(),z.high=c())
screen2.f.posns<-list(z.all=c(),z.low=c(),z.mid=c(),z.high=c())
test.scores<-c("Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green")
for(itest in 1:length(test.scores)){
  the.score<-test.scores[itest]
############### crap filter
crap<-screen2[,RG.label[itest]] <= RG.cut[itest] | is.na(screen2[,the.score])
screen2.f<-screen2[!crap,]

crap<-screen1[,RG.label[itest]] <= RG.cut[itest] | is.na(screen1[,the.score])
screen1.f<-screen1[!crap,]

############ reorder
## the.order<-order(screen2.f[,the.score],decreasing=TRUE)
## screen2.f<-screen2.f[the.order,]
## the.order<-order(screen1.f[,the.score],decreasing=TRUE)
## screen1.f<-screen1.f[the.order,]
            
screen2.cuts[itest]
screen1.cuts[itest]
symbols2<-c(symbols2,screen2.f[screen2.f[,the.score] >= screen2.cuts[itest],"Symbol"])
symbols1<-c(symbols1,screen1.f[screen1.f[,the.score] >= screen1.cuts[itest],"Symbol"])

screen2.f.counts[itest]<-sum(screen2.f[,the.score] >= screen2.cuts[itest])
screen1.f.counts[itest]<-sum(screen1.f[,the.score] >= screen1.cuts[itest])
  
screen2.f.posns[[itest]]<-rownames(screen2.f)[screen2.f[,the.score] >= screen2.cuts[itest]] # positions in unfiltered
screen1.f.posns[[itest]]<-rownames(screen1.f)[screen1.f[,the.score] >= screen1.cuts[itest]]
  
}

length(unique(symbols1))
length(unique(symbols2))
the.common<-intersect(symbols1,symbols2)
length(the.common)

unique(symbols1)
unique(symbols2)
unique(the.common)
hugo.norm.hits<-symbols2


############# enclosed are done WITHOUT taking the intersection #########
conf.norm.hits<-symbols2
conf.norm.hits.posns<-screen2.f.posns
conf.screen.norm<-screen2

hugo.norm.hits<-symbols2
hugo.norm.hits.posns<-screen2.f.posns
hugo.screen.norm<-screen2

joseph.norm.hits<-symbols2
joseph.norm.hits.posns<-screen2.f.posns
joseph.screen.norm<-screen2

kiril.norm.hits<-symbols1
kiril.norm.hits.posns<-screen1.f.posns
kiril.screen.norm<-screen1
########################################################################

kiril.joseph.unnorm.all.commom<-the.common
kiril.joseph.norm.all.commom<-the.common  ### based on wider z-score cut offs


kiril.joseph.norm.all.commom.extended<-the.common
hugo.joseph.norm.all.commom<-the.common


setdiff(kiril.joseph.norm.all.commom,kiril.joseph.unnorm.all.commom)
setdiff(kiril.joseph.unnorm.all.commom,kiril.joseph.norm.all.commom)
intersect(kiril.joseph.norm.all.commom,kiril.joseph.unnorm.all.commom)

intersect(kiril.joseph.norm.all.commom,hugo.joseph.norm.all.commom)
setdiff(kiril.joseph.norm.all.commom,hugo.joseph.norm.all.commom)
setdiff(hugo.joseph.norm.all.commom,kiril.joseph.norm.all.commom)
save.image("HITS for Kiril and joseph.RData")
############ Unnormalized:


hugo.screen.hits<-hugo.screen.norm[unique(unlist(hugo.norm.hits.posns)),]
joseph.screen.hits<-joseph.screen.norm[unique(unlist(joseph.norm.hits.posns)),]
kiril.screen.hits<-kiril.screen.norm[unique(unlist(kiril.norm.hits.posns)),]

write.table(joseph.screen.hits,"joseph_z-score_HITS.txt",row.names=FALSE,sep="\t")
write.table(kiril.screen.hits,"kiril_z-score_HITS.txt",row.names=FALSE,sep="\t")
write.table(hugo.screen.hits,"hugo_z-score_HITS.txt",row.names=FALSE,sep="\t")



 ##########################
my.qq.plot<-function (x, distribution = "norm", ylab = deparse(substitute(x)),
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






      qq.data<- function (x, plot.it = TRUE, distribution = "norm", df=1, the.mean=0, the.sd=1,  xlab = deparse(substitute(x)),
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


#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
##############################################################
#############################################################
#############################################################
write.table(joseph.screen.hits,"joseph_z-score_HITS.txt",row.names=FALSE,sep="\t")
write.table(kiril.screen.hits,"kiril_z-score_HITS.txt",row.names=FALSE,sep="\t")
write.table(hugo.screen.hits,"hugo_z-score_HITS.txt",row.names=FALSE,sep="\t")
##########################################################
 setwd("/media/Bioinform-D/Research/Cellomics/Joseph screen")
setwd("/media/Bioinform-D/Research/Cellomics/Kiril screen")
file<-"joseph_z-score_HITS.txt"
hits<-read.delim(file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized
colnames(hits)<-sum.vars
ann.new<-read.delim("Kiril normalized annotations.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
conf<-ann.new[ann.new[,"plate"] %in% c(30,31,32,33,34),] # this in now the comformation plates
tapply(conf[,"plate"],conf[,"plate"],length)
## 30 31 32 33 34 
## 96 96 96 96 94  34E3  34F1, 34F3, missing
#tpos.joseph<-paste("P",ann[,"plate"],ann[,"row"],ann[,"col"],sep="")
tpos.hits<-paste("P",hits[,"plate"],hits[,"row"],hits[,"col"],sep="")
########## plate 15 is column 12 from pilot plates 1-14
conf[conf[,"KPposition"]=="P14G12","KPposition"]<-"P15G11"
conf[conf[,"KPposition"]=="P12A12","KPposition"]<-"P15A9"
#in.common<-match(conf[,"KPposition"],tpos.joseph)
in.common<-match(conf[,"KPposition"],tpos.hits)
in.common<-unique(in.common)
in.common<-in.common[!is.na(in.common)]  ## 57 in replication
 hit.conf<-hits[in.common,]
write.table(hit.conf,"joseph_z-score_HITS_inCONF.txt",row.names=FALSE,sep="\t")

load("ann1500.RData")

save.image("working.compile.RData")

redd

slope<-{}
inter<-{}
plates.all<-kiril.norm
to.box<-plates.all[filter,c("plate","Trans_%","%R")]

to.box<-plates.all[filter,c("plate","Trans_%","RG/G")]
to.box[,"RG/G"]<-to.box[,"RG/G"]*100


colnames(to.box)<-c("plate","x","y")
the.plates<-unique(to.box[,"plate"])
col.array<-rainbow(length(the.plates))
plot(to.box[,"x"],to.box[,"y"],pch=19,col=col.array[to.box[,"plate"]],cex=0.75,xlab="% Transduction",ylab=" % R")
tapply(to.box[,"plate"],to.box[,"plate"],length)


for(i in 1:length(the.plates)){
  print(i)

  to.box.plate<-to.box[to.box[,"plate"]==the.plates[i],c("x","y")]
  #  plot(to.box.plate[,"x"],to.box.plate[,"y"],pch=19,col=col.array[i],cex=0.75,xlab="% Transduction",ylab=" % R")
  the.model<-lmrob(y~x,data=to.box.plate)
  slope[i]<-coef(the.model)[2]
  inter[i]<-coef(the.model)[1]
  abline(coef(the.model),lty=10,col=col.array[i],lwd=2)
}
savePlot("Transduction vs proliferation of Green.jpeg",type="jpeg")
x11()

#print(as.character(the.model$converged))
lm.red.slope[i,j]<-coef(the.model)[2]
lm.red.inter[i,j]<-coef(the.model)[1]
lm.red.coverg[i,j]<-the.model$converged



boxplot(y~x,data=to.box,ylab="%R",xlab="Plate",main="For filtered with RG >50")
savePlot("Boxplot %R for all plates.jpeg",type="jpeg")

plates.all.f<-plates.all[filter,]

the.order<-order(plates.all.f[,"Trans_%"])
col.array<-rainbow(length(unique(plates.all.f[,"plate"])))

plates.all.f<-plates.all.f[the.order,]
plot(plates.all.f[,"Trans_%"],plates.all.f[,"%R"],pch=19,col=col.array[plates.all.f[,"plate"]],cex=0.75)

points(plates.all.f[,"Trans_%"],plates.all.f[,"RG/G"]*100,pch=21,col=col.array[plates.all.f[,"plate"]],cex=1.1)
legend(70,32,legend=c("% prolif notGreen","%proif Green",labels(col.array)),col<-col.array,cex=0.85)

plot(plates.all.f[,"Trans_%"],plates.all.f[,"RG/G"]*100,pch=19,col=col.array[plates.all.f[,"plate"]])

col.array<-rainbow(length(unique(cells[,"plate"])))
 sample<-rep(TRUE,times=dim(cells)[1])
  plot(cells[sample,"total"],cells[sample,"transducted z-score"],col=col.array[cells[sample,"plate"]],pch=19,cex=1.0,ylim=c(-20,20),xlab="Cells per field",ylab="Transduction z-score (bases on well estimates)",font.lab=2,cex.lab=1.5,main=paste(the.screen,"Screen :  (PLV101-square)",sep=" "))
  abline(h=0)


x11()
to.box<-plates.all[filter.plv,c("plate","%R")]
colnames(to.box)<-c("x","y")
boxplot(y~x,data=to.box,ylab="%R",xlab="Plate",main="For PLV101s with RG >50")
savePlot("Boxplot %R for PLV101 all plates.jpeg",type="jpeg")

to.box<-plates.all[filter.plv,c("plate",the.score)]
colnames(to.box)<-c("x","y")
boxplot(y~x,data=to.box,ylab=the.score,xlab="Plate",main="For PLV101s with RG >50")
savePlot(paste("Boxplot ",the.score," for PLV101 all plates.jpeg",sep=""),type="jpeg")

to.box<-plates.all[filter.plv,c("plate","RG/G")]
colnames(to.box)<-c("x","y")
to.box[,"y"]<-to.box[,"y"]*100
boxplot(y~x,data=to.box,ylab="%RG/G",xlab="Plate",main="For PLV101s filtered with RG >50")
savePlot("Boxplot %RG for PLV101 all plates.jpeg",type="jpeg")

to.box<-plates.all[filter.plv,c("plate","Total")]
colnames(to.box)<-c("x","y")
to.box[,"y"]<-to.box[,"y"]*100
boxplot(y~x,data=to.box,ylab="Total",xlab="Plate",main="For PLV101s filtered with RG >50")
savePlot("Boxplot Total cells for PLV101 all plates.jpeg",type="jpeg")


to.box<-plates.all[filter,c("plate","%R")]
colnames(to.box)<-c("x","y")
boxplot(y~x,data=to.box,ylab="%R",xlab="Plate",main="For filtered with RG >50")
savePlot("Boxplot %R for all plates.jpeg",type="jpeg")

to.box<-plates.all[filter,c("plate",the.score)]
colnames(to.box)<-c("x","y")
boxplot(y~x,data=to.box,ylab=the.score,xlab="Plate",main="For filtered with RG >50")
savePlot(paste("Boxplot ",the.score," for all plates.jpeg",sep=""),type="jpeg")

to.box<-plates.all[filter,c("plate","RG/G")]
colnames(to.box)<-c("x","y")
to.box[,"y"]<-to.box[,"y"]*100
boxplot(y~x,data=to.box,ylab="%RG/G",xlab="Plate",main="For filtered with RG >50")
savePlot("Boxplot %RG for all plates.jpeg",type="jpeg")

to.box<-plates.all[filter,c("plate","Total")]
colnames(to.box)<-c("x","y")
to.box[,"y"]<-to.box[,"y"]*100
boxplot(y~x,data=to.box,ylab="Total",xlab="Plate",main="For filtered with RG >50")
savePlot("Boxplot Total cells for all plates.jpeg",type="jpeg")



#,notch=TRUE)
boxplot(y~x,data=plates.all[filter.plv,])
#data.plate<-(plates.all[filter ,c("plate",the.score)])
plates.all[(plates.all[,"Symbol"]=="PLV101" | plates.all[,"Symbol"]=="plv101"),]
plates.all[(plates.all[,"Symbol"]=="PLV101" | plates.all[,"Symbol"]=="plv101") & plates.all[,"plate"]==19 ,1:24]
plates.all[(plates.all[,"Symbol"]=="PLV101" | plates.all[,"Symbol"]=="plv101") & plates.all[,"plate"]==25 ,1:24]
