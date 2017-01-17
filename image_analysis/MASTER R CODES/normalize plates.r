############################################### millian  ### latest DNA and flexible annotation
setwd("/media/Bioinform-D/Research/Cellomics/Kiril screen")
load("ann2800.RData")  ### "plate" must match "plate" in file  cellomics to expt map.csv
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","G1","G2","S","Gr-S+G2/G1","Gr-S+G2+ >4N/G1","Gr-<2N","Gr->4N","Gr-G1","Gr-G2","Gr-S","NtGr-S+G2/G1","NtGr-S+G2+ >4N/G1","NtGr-<2N","NtGr->4N","NtGr-G1","NtGr-G2","NtGr-S"
             )

the.screen<-"Kiril"
files<- paste("plate_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Kiril_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Kiril.NOTGREEN3.Wells.txt"
well.type<-96
row.type<-8
col.type<-2
exported.well.data<-"well_based_data.txt"
#######################################################


############################################### Plate 2B
setwd("/media/Bioinform-D/Research/Cellomics/Hits Plate 5")
load("annPlate2B.RData")

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S")



the.screen<-"Plate2B"
files<- paste("plate_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Plate2B_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Plate2B_summary_NOTGREEN_DNA.txt"
well.type<-96
row.type<-8
col.type<-12
#######################################################

############################################### Plate 5
setwd("/media/Bioinform-D/Research/Cellomics/Hits Plate 5")
load("annPlate5.RData")

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S")



the.screen<-"Plate5"
files<- paste("plate_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Plate5_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Plate5_summary_NOTGREEN_DNA.txt"
well.type<-96
row.type<-8
col.type<-12
#######################################################



############################################### Helga
setwd( "/media/Bioinform-D/Research/Cellomics/Helga")
#load("Ca_ann2.RData") # same as Ca_ann2 except 
map.file<-"cellomics.to.expt.map.csv" # "plate" must be same "plate in annoation AND IN THE SAME ORDER
annotations.file<-"Helga_annotations2.txt"

sum.vars<-c("plate","row","col","Symbol","%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S")

the.screen<-"Helga"
#files<- paste("Ca_siCa_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="") # defined below
field.output.file<-"Helga_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Helga_summary_NOTGREEN_DNA.txt"
well.type<-96
row.type<-8
col.type<-12

#######################################################
#########################################
########################################
Use this file to combine well based and cell based data typicall done after normalization


############################################### Ca screen 2 " same as above but the very latest
setwd( "/media/Bioinform-D/Research/Cellomics/Ca screen/Ca screen 2")
load("Ca2_ann.RData")  ### renames the SUMMARY files using the map file
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","G1","G2","S","Gr-S+G2/G1","Gr-S+G2+ >4N/G1","Gr-<2N","Gr->4N","Gr-G1","Gr-G2","Gr-S","NtGr-S+G2/G1","NtGr-S+G2+ >4N/G1","NtGr-<2N","NtGr->4N","NtGr-G1","NtGr-G2","NtGr-S" )


the.screen<-"Ca siRNA2"
files<- paste("plate_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Ca2_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Ca siRNA2.NOTGREEN3.Wells.txt"
well.type<-384
row.type<-16
col.type<-24

map.file<-"cellomics.to.expt.map.csv" #
exported.well.data<-"well_based_data.txt"

#######################################################

###############################################Ca screen 3  ### latest DNA and flexible annotation
setwd("/media/Bioinform-D/Research/Cellomics/Ca screen/Ca screen 3")
load("Ca3_ann.RData")
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann


core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S")

the.screen<-"Ca siRNA3"
files<- paste("plate_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Ca3_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Ca siRNA3.NOTGREEN3.Wells.txt" # name now after adding well data
well.type<-384
row.type<-16
col.type<-24

map.file<-"cellomics.to.expt.map.csv" #
exported.well.data<-"well_based_data.txt"
#######################################################
###START#########################
 setwd("/media/Bioinform-D/Research/Cellomics/Joseph screen/joseph redo 2010")
sum.vars<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","Description","GenBank_Accn","dbEST_Id" ,"MGC.Accession")
well.output.file<-"Joseph_summary_NOTGREEN_DNA.txt"
the.screen="Joseph"
plates.all<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized
colnames(plates.all)<-sum.vars
############################################


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
well.output.file<-"siFB.NOTGREEN3.Wells.txt"
well.type<-384
row.type<-16
col.type<-24
exported.well.data<-"well_based_data.txt"
#######################################################

############################################### JOseph 2010  ### latest DNA and flexible annotation
setwd("/media/Bioinform-D/Research/Cellomics/Joseph screen 2010")
load("ann1500.RData")  ### "plate" must match "plate" in file  cellomics to expt map.csv
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","G1","G2","S","Gr-S+G2/G1","Gr-S+G2+ >4N/G1","Gr-<2N","Gr->4N","Gr-G1","Gr-G2","Gr-S","NtGr-S+G2/G1","NtGr-S+G2+ >4N/G1","NtGr-<2N","NtGr->4N","NtGr-G1","NtGr-G2","NtGr-S"
             )

the.screen<-"Joseph_2010"
files<- paste("plate_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Joseph_2010_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Joseph_2010_summary_NOTGREEN_DNA.txt"
well.type<-96
row.type<-8
col.type<-12

normalized.file<-paste(the.screen,"NORMALIZED","txt",sep=".")
######################################################

############################################### Jane  ### latest DNA and flexible annotation
setwd("/media/Bioinform-D/Research/Cellomics/Jane")
load("Jane_ann.RData")  ### "plate" must match "plate" in file  cellomics to expt map.csv
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","G1","G2","S","Gr-S+G2/G1","Gr-S+G2+ >4N/G1","Gr-<2N","Gr->4N","Gr-G1","Gr-G2","Gr-S","NtGr-S+G2/G1","NtGr-S+G2+ >4N/G1","NtGr-<2N","NtGr->4N","NtGr-G1","NtGr-G2","NtGr-S")

the.screen<-"Jane"
files<- paste("plate_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Jane_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Jane_summary_NOTGREEN_DNA.txt"
well.type<-96
row.type<-8
col.type<-12

normalized.file<-paste(the.screen,"NORMALIZED","txt",sep=".")
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

normalized.file<-paste(the.screen,"NORMALIZED","txt",sep=".")
#######################################################


####################################START HERE
library(robust)
library(robustbase)
library(nls2)
library(akima)

############# read in the reults one screen or many screens #############
plates.all<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)


## colnames(plates.all)[1:length(sum.vars)]<-sum.vars # old way
#################################################################
#########flexible annotion way
col<-colnames(plates.all) %in% core.ann
labels<-c(colnames(plates.all)[col],core.vars,colnames(ann)[!place.core.ann])## 
colnames(plates.all)[1:length(labels)]<-labels  ## takes account of possibility of added well based data

##########one off fixed
## plates.all[,"plate"]<-gsub("siCA3_","",plates.all[,"plate"])


############################################ START HERE IS HAVE PLATES_ALL 

############################# fix the column lables so no / or % so files names are ok
colnames(plates.all)
labels<-gsub("%","percent_",colnames(plates.all))
labels<-gsub("/","_over_",labels)
colnames(plates.all)<-labels
colnames(plates.all)
plates.all.ori<-plates.all
###########################set up normalization parameters:
sort(tapply(plates.all[,"Symbol"],plates.all[,"Symbol"],length)) ## look for controls

if(the.screen=="Jane" | the.screen=="Jane_validation"){
  plates.all<-cbind(plates.all,plates.all[,"RG_over_G"]); colnames(plates.all)[dim(plates.all)[2]]<-"RG_over_G_Z"
}
                       
############ 2 color overexpression choose what to normalize:
test.scores<-c( "G1 posn","Num-cells","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","Relative S+G2+>4N","Gr-S+G2+ >4N_over_G1","RG_over_G","RG_over_G_Z","Gr-G1"  )
zero.center<-c(FALSE, FALSE, TRUE, TRUE,  TRUE,   TRUE, FALSE, FALSE, FALSE,TRUE,FALSE)
use.mean<-rep(FALSE,times=length(test.scores)); names(use.mean)<-test.scores
use.mean["RG_over_G_Z"] <-TRUE
make.Z.score<-rep(FALSE,times=length(test.scores)); names(make.Z.score)<-test.scores
make.Z.score["RG_over_G_Z"] <-TRUE      
use.filter.RG.expected.min<-TRUE ;  RG.expected.min<-10 # typically 50 10 for Jane
use.filter.G.expected.min<-TRUE ;  G.expected.min<-50 # typically 50
use.medians=TRUE # usually TRUE (if false uses norm.control.list to normalize)
if(the.screen=="Jane_validation"){
 use.medians=FALSE
}
#  # 
## use.plv=TRUE # if True it overrides use med

##### for JANE###
pos.control.list<-c("CCND1","MOCK","HSPCN111-N","UNTRANSDUCED","no DNA","PLV411") # controls that would affect mediann center normalization
norm.control.list<-c("HSPCN111-N") # controls that could be used for  normalization
######

##### for JANE VALIDATION###
pos.control.list<-c("CCND1","C-MYC","HSPCN","UTMOCK","MOCK") # controls that would affect mediann center normalization
norm.control.list<-c("plv411") # controls that could be used for  normalization
######

#### for other
pos.control.list<-c("CCNE1","101-CCNE1","411-CCNE1","101-CYCD2","not transduced","MOCK") # controls that would affect mediann center normalization
norm.control.list<-c("PLV101","plv101","PLV101-2","CPNE3","PSTPIP2") # controls that could be used for  normalization

###################################
######for single color siRNA
## test.scores<-c("G1.posn")
use.medians=FALSE # usually TRUE (if false uses norm.control.list to normalize)
test.scores<-c("G1 posn","S+G2_over_G1","percent_R", "ak_leakage","resazurin","Num-cells", "S+G2+ >4N_over_G1"  )
zero.center<-c(FALSE, FALSE,FALSE, FALSE,FALSE,FALSE,FALSE)
use.mean<-rep(FALSE,times=length(test.scores)); names(use.mean)<-test.scores
## use.mean["RG_over_G_Z"] <-TRUE
make.Z.score<-rep(FALSE,times=length(test.scores)); names(make.Z.score)<-test.scores
## make.Z.score["RG_over_G_Z"] <-TRUE  
use.filter.RG.expected.min<-FALSE ;  RG.expected.min<-15 # typically 50
use.filter.G.expected.min<-FALSE ;  G.expected.min<-50 # typically 50


pos.control.list<-c("LAMA","GFP","PLK1","ATP2C1","empty","No Cells","No Lipid","mediaonly","cellonly","old control")
norm.control.list<-c("NT") # controls that could be used for  normalization
use.filter.RG.expected.min<-FALSE # not relevant for one color screen
##############



length(zero.center)==length(test.scores) # must be TRUE
plates.all[1:5,test.scores] # test have names correct


##########################################################################
################### NORMALIZATION LOOP
the.plates<-as.character(unique(plates.all[,"plate"])) #### modified Oct 26th
col.array<-rainbow(length(the.plates))



for(itest in 1:length(test.scores)){
the.score<-test.scores[itest]
the.score
data.in<- (plates.all[is.finite(plates.all[,the.score]),the.score])

jpeg( paste(the.screen,"_raw_",the.score,".jpeg",sep=""), width = 1280, height = 1024, units = "px",quality = 100 )

ahist<-hist(data.in,col="pink",breaks=100,main=paste(the.screen," UnNormalised for ",toString(the.score),sep=""),xlab=the.score,cex.lab=1.5,cex.main=2.0)
aden<-density(data.in,n=1024,adjust=1)
apeak<-max(ahist$counts)/max(aden$y)

if(length(data.in)>5000){the.test<-shapiro.test(data.in[1:4999])}else{the.test<-shapiro.test(data.in)}
text(ahist$breaks[floor((length(ahist$breaks)/2)+5)],max(ahist$counts)*0.9, paste("Shapiro Test:",signif(the.test$p.value,3),sep=" "),cex=1.5)
lines(aden$x,apeak*aden$y,col="green",lwd=2)
dev.off()
#savePlot(paste(the.screen,the.score,".jpeg",sep=""),type="jpeg")

the.medians<-tapply(plates.all[,the.score],plates.all[,"plate"], function(x) median(x,na.rm=TRUE))
the.fit.medians<-the.medians



###################### filter the results
filter<-is.finite(plates.all[,the.score])


########################################### SET positive controls - exclude from normalization  ############################
pos.control<-plates.all[,"Symbol"] %in% pos.control.list  # identify positive controls
#################################################################################################################

########################################### exclude low R and G weels from normalizationn  #####################

if(use.filter.RG.expected.min){RG.small<-plates.all[,"RG_expected"] < RG.expected.min & !is.na(plates.all[,"RG_expected"]) ; filter<-filter & !RG.small}
if(use.filter.G.expected.min){G.small<-plates.all[,"G"] < G.expected.min & !is.na(plates.all[,"G"]) ; filter<-filter & !G.small}

filter<-filter &  !pos.control
################################################################################################################

########################################### SET the control wells for normalization  ############################
filter.plv<-(plates.all[,"Symbol"] %in% norm.control.list) & is.finite(plates.all[,the.score])
if(use.filter.RG.expected.min){filter.plv<-filter.plv & !RG.small}#  NT
#################################################################################################################


the.plv.medians<-tapply(plates.all[filter.plv,the.score],plates.all[filter.plv,"plate"], function(x) median(x,na.rm=TRUE))
the.filtered.medians<-tapply(plates.all[filter,the.score],plates.all[filter,"plate"], function(x) median(x,na.rm=TRUE))

the.plv.means<-tapply(plates.all[filter.plv,the.score],plates.all[filter.plv,"plate"], function(x) mean(x,na.rm=TRUE))
the.filtered.means<-tapply(plates.all[filter,the.score],plates.all[filter,"plate"], function(x) mean(x,na.rm=TRUE))

the.plv.sd<-tapply(plates.all[filter.plv,the.score],plates.all[filter.plv,"plate"], function(x) sd(x,na.rm=TRUE))
the.filtered.sd<-tapply(plates.all[filter,the.score],plates.all[filter,"plate"], function(x) sd(x,na.rm=TRUE))

the.plv.medians<-the.plv.medians[the.plates]
the.filtered.medians<-the.filtered.medians[the.plates]
names(the.plv.medians)[is.na(the.plv.medians)] <-the.plates[is.na(the.plv.medians)] 
names(the.filtered.medians)[is.na(the.filtered.medians)] <-the.plates[is.na(the.filtered.medians)] 

the.plv.means<-the.plv.means[the.plates]
the.filtered.means<-the.filtered.means[the.plates]
names(the.plv.means)[is.na(the.plv.means)] <-the.plates[is.na(the.plv.means)] 
names(the.filtered.means)[is.na(the.filtered.means)] <-the.plates[is.na(the.filtered.means)] 

the.plv.sd<-the.plv.sd[the.plates]
the.filtered.sd<-the.filtered.sd[the.plates]
names(the.plv.sd)[is.na(the.plv.sd)] <-the.plates[is.na(the.plv.sd)] 
names(the.filtered.sd)[is.na(the.filtered.sd)] <-the.plates[is.na(the.filtered.sd)] 


if(use.mean[the.score]){
the.plv.medians<-the.plv.means
the.filtered.medians<-the.filtered.means
}

   
#the.filtered.medians<-tapply(plates.all[filter & !low.val & !high.val,the.score],plates.all[filter & !low.val & !high.val,"plate"], function(x) median(x,na.rm=TRUE))

number<-tapply(plates.all[filter,"plate"],plates.all[filter,"plate"],length)
the.number.plv.medians<-tapply(plates.all[filter.plv,the.score],plates.all[filter.plv,"plate"], length)
print(paste("Number if plv controls for median",the.number.plv.medians,sep=" "))
print(paste("Number if cells for median",number,sep=" "))
      
jpeg( paste(the.screen,"_BOXPLOT_PLATE_MEDIANS_",the.score,".jpeg",sep=""), width = 1280, height = 1024, units = "px",quality = 100 )
#to.box<-plates.all[filter.plv,c("plate",the.score)]
to.box<-plates.all[filter,c("plate",the.score)]
all.data.median<-median(to.box[,2],na.rm=TRUE)
all.data.sd<-sd(to.box[,2],na.rm=TRUE)
all.data.min<-min(to.box[,2])

the.box.max<-all.data.median + 6*all.data.sd
the.box.min<-min(to.box[,2])
colnames(to.box)<-c("x","y")
#the.boxplot<-boxplot(y~x,data=to.box,ylab=the.score,xlab="Plate",main="For PLV101s with RG >50")
the.boxplot<-boxplot(y~x,data=to.box,ylab=the.score,xlab="Plate",main=paste(the.screen," BOXPLOT PLATE MEDIANS-NORMALIZED ",the.score,".jpeg",sep=" "),cex.lab=1.5,cex.axis=1.5,cex.main=2,col=rainbow(length(the.plates)),ylim=c(the.box.min,the.box.max))
text(0,the.box.max,labels=toString(paste(the.boxplot$names,"(",the.boxplot$n,")",sep="")),pos=4,cex=0.9,cex.lab=1.5,font=1.5,col="red")

dev.off()
#savePlot( paste(the.screen,"_BOXPLOT_PLATE_MEDIANS_",the.score,".jpeg",sep=""),type="jpeg")

jpeg( paste(the.screen,"_BOXPLOT_CONTROL_MEDIANS_",the.score,".jpeg",sep=""), width = 1280, height = 1024, units = "px",quality = 100 )
#to.box<-plates.all[filter.plv,c("plate",the.score)]
to.box<-plates.all[filter.plv,c("plate",the.score)]
colnames(to.box)<-c("x","y")
the.boxplot<-boxplot(y~x,data=to.box,ylab=the.score,xlab="Plate",main=paste(the.score," For Controls : ",toString(norm.control.list),sep=""),col=rainbow(length(the.boxplot$names)),ylim=c(the.box.min,the.box.max))
#the.boxplot<-boxplot(y~x,data=to.box,ylab=the.score,xlab="Plate",main="Boxplot for PLATE CONTROLS :ATP2C1",cex.lab=1.5,cex.axis=1.5,cex.main=2,col=rainbow(length(the.plates)))
text(1,max(to.box$y),labels=toString(paste(the.boxplot$names,"(",the.boxplot$n,")",sep="")),pos=4,cex=0.9,cex.lab=1.5,font=1.5,col="red")
dev.off()
#savePlot(paste(the.screen,"_BOXPLOT_CONTROL_MEDIANS_",the.score,".jpeg",sep=""),type="jpeg")

################### if the plv.means runs out for the plate use the filtered medians
if(length(the.filtered.medians) != length(the.plates)){print("ERROR NORMALIZATION FAILURE - no genes on plate")}


## the.plv.medians<-the.plv.medians[as.character(the.plates)]
## the.filtered.medians<-the.filtered.medians[as.character(the.plates)]
## missing<-is.na(the.plv.medians)
## the.plv.medians[missing]<-the.filtered.medians[missing]
## names(the.plv.medians)[missing]<-names(the.filtered.medians)[missing]

#x11( width = 11, height = 11)
jpeg( paste(the.screen,"_raw_ALL_WELLS_",the.score,".jpeg",sep=""), width = 1280, height = 1024, units = "px",quality = 100 )

nf<-layout(matrix(c(1:((ceiling(length(the.plates)/3))*3)),(ceiling(length(the.plates)/3)),3,byrow=TRUE),heights=c(rep(1,times=ceiling(length(the.plates)/3)) ),widths=c(rep(1,times=3) ) )
layout.show(nf)
par(mar=c(3,3.5,2.1,2.1),mgp=c(2,1,0)) #c(bottom, left, top, right)

for(i in 1:length(the.plates)){
data.plate<-(plates.all[plates.all[,"plate"]==the.plates[i] & filter, the.score])
data.plate<- data.plate[data.plate< (all.data.median + 8*all.data.sd)]

if(length(data.plate)<1){
data.plate<-c(all.data.median,all.data.median,all.data.median)
ahist.plate<-hist(data.plate,col="white",border="white",breaks=100,main=paste("Plate=",the.plates[i],sep=""),xlab=the.score,cex.main=1.5,cex.lab=1.5,cex.axis=2,xlim=c(all.data.min,(all.data.median + 4*all.data.sd)))
text(all.data.median,2,"NO DATA TO PLOT",col="red",cex=2.5,pos=4,font=1.5)
}else{
ahist.plate<-hist(data.plate,col=col.array[i],breaks=100,main=paste("Plate=",the.plates[i],sep=""),xlab=the.score,cex.main=1.5,cex.lab=1.5,cex.axis=2,xlim=c(all.data.min,(all.data.median + 4*all.data.sd)))
}


data.plate<-data.plate[!is.na(data.plate)]
aden<-density(data.plate,n=1024,adjust=1)
apeak<-max(ahist.plate$counts)/max(aden$y)
lines(aden$x,apeak*aden$y,col="green",lwd=2)
the.fit.medians[names(the.fit.medians)==the.plates[i]]<-coef(ltsreg(data.plate~1 ))
#the.fit.medians1<-the.fit.medians[the.plates[i]]
 #the.model<- ltsreg(data ~ plate , data=data.plate)


#abline(v=the.medians[names(the.medians)==the.plates[i]],col="blue",lwd=3)
#abline(v=the.fit.medians[names(the.fit.medians)==the.plates[i]],col="orange2",lwd=3)
abline(v=the.filtered.medians[names(the.filtered.medians)==the.plates[i]],col="magenta",lwd=4)
abline(v=the.plv.medians[names(the.plv.medians)==the.plates[i]],col="cyan",lwd=3,lty=2)
text(max(ahist.plate$mid),max(ahist.plate$counts)-0.5,paste("median= ",round(the.filtered.medians[names(the.filtered.medians)==the.plates[i]],2),"(",number[names(number)==the.plates[i]],")",sep=""),cex=1.5,pos=2,font=1.5)
}

#savePlot(paste(the.screen,"_raw_ALL_WELLS_",the.score,".jpeg",sep=""),type="jpeg")
dev.off()

#x11()
jpeg(paste(the.screen,"_NORMALIZED_",the.score,".jpeg",sep=""), width = 1280, height = 1024, units = "px",quality = 100 )

############### meadian normalize data### cound extend extend this to do some be controls some by median
if(use.medians){
centers.all<-the.filtered.medians
sd.all<-the.filtered.sd
}else{
### centers.all<-the.fit.medians
  centers.all<-the.plv.medians
  sd.all<-the.plv.sd
}

#centers.all<-the.plv.medians
if(make.Z.score[the.score]){
  for(i in 1:length(the.plates)){
  plates.all[plates.all[,"plate"]==the.plates[i],the.score] <-(plates.all[plates.all[,"plate"]==the.plates[i],the.score]-
  centers.all[names(centers.all)==the.plates[i]])/sd.all[names(sd.all)==the.plates[i]]
  }
}else{
    for(i in 1:length(the.plates)){
    plates.all[plates.all[,"plate"]==the.plates[i],the.score] <-plates.all[plates.all[,"plate"]==the.plates[i],the.score]-
    centers.all[names(centers.all)==the.plates[i]]
    }
}

plates.all[plates.all[,"plate"]==the.plates[i],the.score][1:100]
test<-plates.all[plates.all[,"plate"]==the.plates[i],the.score]-centers.all[names(centers.all)==the.plates[i]]
plates.all[plates.all[,"plate"]==the.plates[i],the.score][1:5]
#### modification so don't have neagiatve values
if(!zero.center[itest]){
  min.value<-abs((min(plates.all[,the.score],na.rm=TRUE)))
 plates.all[,the.score]<- plates.all[,the.score]+min.value
}

data.in<- (plates.all[filter,the.score])
the.mean<-median(data.in[filter],na.rm=TRUE)
the.sd<-sd(data.in[filter],na.rm=TRUE)


data.in<-data.in[data.in <= (the.mean+ 3*the.sd)]

ahist<-hist(data.in,col="pink",breaks=25,main=paste(the.screen," Normalised for ",the.score,sep=""),xlab=the.score,cex.lab=1.5,cex.main=2.0)




aden<-density(data.in,n=1024,adjust=1)
apeak<-max(ahist$counts)/max(aden$y)
#apeak<-30/max(aden$y)
lines(aden$x,apeak*aden$y,col="green",lwd=2)


####################### fit normal to data


 to.fit<-data.frame(y=aden$y*apeak, x=aden$x)
A <- max(to.fit$y)/dnorm(the.mean,mean=the.mean,sd=the.sd, log=FALSE)
  try(the.fit<-nls(y~A*dnorm(x,mean=the.mean,sd=the.sd, log=FALSE)
               ,data=to.fit
               ,start=list(A=A,the.mean=the.mean,the.sd=the.sd)
                ,lower=list(A/10,the.mean/1.5,the.sd/3.5)
                ,upper=list(A*10,the.mean*1.5,the.sd*3.5)
                ,trace=TRUE
               ,control=list(maxiter=1000, minFactor=1/4048,tol=1e-4,algorithm="port", warnOnly = TRUE)),silent=TRUE)


if(!exists("the.fit")){
  text(the.mean+the.sd,max(to.fit$y),"Failed to fit normal",cex=1.5,pos=4)
 }else{
coef(the.fit)
curve(coef(the.fit)["A"]*dnorm(x,mean=coef(the.fit)["the.mean"],sd=coef(the.fit)["the.sd"], log=FALSE),add=TRUE, col="red",lwd=2,lty="dashed")

wings<-(data.in > (coef(the.fit)["the.mean"]+coef(the.fit)["the.sd"]*3) ) | (data.in < (coef(the.fit)["the.mean"]-coef(the.fit)["the.sd"]*3))

if(sum(!wings)>15){
if(length(data.in[!wings])>5000){the.test<-shapiro.test(data.in[!wings][1:4999])}else{the.test<-shapiro.test(data.in[!wings])}
}

text(coef(the.fit)["the.mean"]+coef(the.fit)["the.sd"]*2,max(to.fit$y),paste("Shapiro Test:",signif(the.test$p.value,3),sep=" "),cex=1.5,pos=4)
if(!zero.center[itest]){
text(coef(the.fit)["the.mean"]+coef(the.fit)["the.sd"]*2,max(to.fit$y),paste("Distribution centered to",min.value,sep=" "),cex=1.5,pos=1)}
}


#### for publication
## par(mar=c(5,6.0,2.5,5.1),mgp=c(3.5,1,0))
## ahist<-hist(data.in,col="pink",breaks=25,xlab="z-score",font=2,font.lab=2,cex.axis=2,cex.lab=3,cex.main=2.0,xlim=c(-15,15),main="")
## curve(coef(the.fit)["A"]*dnorm(x,mean=coef(the.fit)["the.mean"],sd=coef(the.fit)["the.sd"], log=FALSE),add=TRUE, col="red",lwd=3,lty="solid")
## text(-16,135,paste("Shapiro Test p-value=",signif(the.test$p.value,2),sep=" "),cex=2,pos=4,font=2)
## savePlot("normalized data.tiff",type="tiff")
## savePlot("normalized data.jpeg",type="jpeg")
dev.off()

}  ### end loop over normalization

#normalized<-plates.all[,c(colnames(plates.all)[1:4],test.scores,colnames(plates.all)[27:36],colnames(plates.all)[45:dim(plates.all)[2]])]
normalized<-plates.all
#write.table(normalized,paste(the.screen,"NORMALIZED","txt",sep="."),row.names=FALSE,sep="\t")
write.table(normalized,paste(the.screen,"NORMALIZED","txt",sep="."),row.names=FALSE,sep="\t")


###### construct inter screen ratios
if(the.screen=="Jane"){
  tapply(normalized[,"plate"],normalized["plate"],length)
  mus.a.expt<-grepl("MUS.a",normalized[,"plate"])
  mus.b.expt<-grepl("MUS.b",normalized[,"plate"])
  mus.c.expt<-grepl("MUS.c",normalized[,"plate"])

 sum( normalized[mus.a.expt,"row"] !=  normalized[mus.c.expt,"row"])
sum( normalized[mus.a.expt,"col"] !=  normalized[mus.b.expt,"col"])

  
  mus.a.nums<-normalized[mus.a.expt,"RG_over_G"]
  mus.a.nums<-c(mus.a.nums,mus.a.nums,mus.a.nums)
  ratios[mus.a.expt]
  
  ratios[mus.a.expt]<-normalized[,"RG_over_G"]
  
 normalized<-cbind(normalized,normalized[,"RG_over_G"]/mus.a.nums)
  colnames(normalized)[dim(normalized)[2]]<-"RG_over_G_over_EtoH"
  normalized[,"RG_over_G_over_EtoH"]
}
###### construct inter screen ratios
if(the.screen=="Jane_validation"){
  tapply(normalized[,"plate"],normalized["plate"],length)
  mus.1a.expt<-grepl("1a",normalized[,"plate"])
  mus.2a.expt<-grepl("2a",normalized[,"plate"])
  mus.3a.expt<-grepl("3a",normalized[,"plate"])

  mus.1b.expt<-grepl("1b",normalized[,"plate"])
  mus.2b.expt<-grepl("2b",normalized[,"plate"])
  mus.3b.expt<-grepl("3b",normalized[,"plate"])

  mus.1c.expt<-grepl("1c",normalized[,"plate"])
  mus.2c.expt<-grepl("2c",normalized[,"plate"])
  mus.3c.expt<-grepl("3c",normalized[,"plate"])

  mus.4a.expt<-grepl("4a",normalized[,"plate"])
  mus.4b.expt<-grepl("4b",normalized[,"plate"])
  mus.4c.expt<-grepl("4c",normalized[,"plate"])

 sum( normalized[mus.1a.expt,"row"] !=  normalized[mus.1c.expt,"row"])
sum( normalized[mus.1a.expt,"col"] !=  normalized[mus.1b.expt,"col"])
  
  sum( normalized[mus.2a.expt,"row"] !=  normalized[mus.2c.expt,"row"])
sum( normalized[mus.2a.expt,"col"] !=  normalized[mus.2b.expt,"col"])
  
  sum( normalized[mus.3a.expt,"row"] !=  normalized[mus.3c.expt,"row"])
sum( normalized[mus.3a.expt,"col"] !=  normalized[mus.3b.expt,"col"])
  
  sum( normalized[mus.4a.expt,"row"] !=  normalized[mus.4c.expt,"row"])
sum( normalized[mus.4a.expt,"col"] !=  normalized[mus.4b.expt,"col"])

  > unique(normalized[,"plate"])
 [1] "1a" "1b" "1c" "2a" "2b" "2c" "3a" "3b" "3c" "4a" "4b" "4c"
  
  mus.1a.nums<-normalized[mus.1a.expt,"RG_over_G"]
  mus.2a.nums<-normalized[mus.2a.expt,"RG_over_G"]
  mus.3a.nums<-normalized[mus.3a.expt,"RG_over_G"]
  mus.4a.nums<-normalized[mus.4a.expt,"RG_over_G"]
  
  mus.a.nums<-c(mus.1a.nums,mus.1a.nums,mus.1a.nums,mus.2a.nums,mus.2a.nums,mus.2a.nums,mus.3a.nums,mus.3a.nums,mus.3a.nums,mus.4a.nums,mus.4a.nums,mus.4a.nums)
  no.greens<-is.na(mus.a.nums) | !is.finite(mus.a.nums) | mus.a.nums==0
  ## ratios[mus.a.expt]
  length(mus.a.nums)
  dim(normalized)
  ## ratios[mus.a.expt]<-normalized[,"RG_over_G"]
  
 normalized<-cbind(normalized,normalized[,"RG_over_G"]/mus.a.nums)
  colnames(normalized)[dim(normalized)[2]]<-"RG_over_G_over_EtoH"
  normalized[no.greens,"RG_over_G_over_EtoH"]<-1
   normalized[,"RG_over_G_over_EtoH"]
   normalized[1:50,"RG_over_G_over_EtoH"]
   normalized[ is.na(normalized[,"RG_over_G_over_EtoH"]),"RG_over_G_over_EtoH"]<-1
}

#####################################
normalized.ori.ca3<-normalized
normalized.ori.ca2<-normalized
normalized <- plates.all.ori # for plotting unnormalized
normalized <- normalized.ori.ca2
#############################
setwd("/media/Bioinform-D/Research/Cellomics/Jane/Validation/Green 32 run")
normalized<-read.delim(paste(the.screen,"NORMALIZED","txt",sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized


col<-colnames(normalized) %in% core.ann
labels<-c(colnames(normalized)[col],core.vars,colnames(ann)[!place.core.ann])## 
colnames(normalized)[1:length(labels)]<-labels
plates.all<-normalized


normalized[1:5,]


var.x<-"Z-all Green"
var.y<-"RG_over_G_Z"

var.x<-"Z.all.Green"
var.y<-"RnG_over_nG"

var.x<-"Z.all.Green"
var.y<-"Trans_percent_"

var.x<-"Z.all.Green"
var.y<-"G"


pos.control.list<-c("NT","LAMA","GFP","PLK1","ATP2C1","empty","No Cells","No Lipid")
norm.control.list<-c("NT") ## NT for Ca2  No Lipid

pos.control.list<-c("CCNE1","101-CCNE1","411-CCNE1","101-CYCD2","not transduced","MOCK") # controls that would affect mediann center normalization
norm.control.list<-c("PLV101","plv101","PLV101-2","CPNE3","PSTPIP2")

pos.control.list<-c("CCND1","MOCK","HSPCN111-N","UNTRANSDUCED","no DNA","PLV411") # controls that would affect mediann center normalization
norm.control.list<-c("HSPCN111-N")

pos.control.list<-c("CCND1","C-MYC") # controls that would affect mediann center normalization
norm.control.list<-c("plv411") # controls that could be used for  normalization

normalized<-normalized[normalized[,"RG"]>49,]
normalized<-normalized[grepl("c",normalized[,"plate"]),]

library(robust)
library(robustbase)
the.model<-lmrob(RG_over_G_Z~ Z-all.Green ,data=normalized)
summary(the.model)

plot(normalized[,var.x],normalized[,var.y],ylab=var.y,xlab=var.x,main=paste("Normalized data used -",the.screen,sep=""),cex.lab=1.5,cex.main=1.5)

## plot(normalized[,"Z-all Green"],normalized[,"RnG_over_nG"],ylab="Red NOt Green/ (Not Green)",xlab="Z-all Green",main="Normalized data used - Joseph Screen")
hit.cut<-5
abline(v=hit.cut,col="red",lwd=2)

abline(coef(the.model),lty=10,col="red",lwd=2)

text(-10,15,"R^2=0.7 ; p-value: < 2e-16",cex=1.5,font=2)


pos.control<-normalized[,"Symbol"] %in% pos.control.list
points(normalized[pos.control,var.x],normalized[pos.control,var.y],col="red",pch=19,cex=2.0 )
normalized[pos.control,c("plate","row","col","Symbol",var.x,var.y)]


norm.control<-normalized[,"Symbol"] %in% norm.control.list
points(normalized[norm.control,var.x],normalized[norm.control,var.y],col="blue",pch=19,cex=2.0 )

hit.cut<-8.9
hits<-normalized[,"Z-all Green"]>8.9
points(normalized[hits,var.x],normalized[hits,var.y],col="orange",pch=19,cex=1.0 )

mock.control.list<-c("UTMOCK")
mock.control<-normalized[,"Symbol"] %in% mock.control.list
points(normalized[mock.control,var.x],normalized[mock.control,var.y],col="green",pch=20,cex=2.0 )

leg.txt<-c("CCNE1","PLV101","Z-score Hit")


legend(0.3,0.1,leg.txt,col=c("red","blue","orange"),pch=19,cex=2.0)

legend(-12,0.1,leg.txt,col=c("red","blue"),pch=19,cex=2.0)

abline(a=0,b=1,lwd=2)

text(11,0.01," <-Hit Cutoff",col="red",cex=2.0,pos=1)

selected.data<-identify(normalized[,var.x],normalized[,var.y],labels=paste(normalized[,"Symbol"],"",normalized[,"plate"],normalized[,"row"],normalized[,"col"],sep=":"),col="forestgreen",font=2,cex=1.0,offset=1,atpen='TRUE',plot=TRUE)

savePlot("green32.jpeg",type="jpeg")
savePlot("green128 treatc.jpeg",type="jpeg")
savePlot("Duka.plot.jpeg",type="jpeg")
savePlot("Duka.plot.ZvsTransduction.jpeg",type="jpeg")
savePlot("Duka.plot.ZvsNumGreen.jpeg",type="jpeg")
savePlot("Proliferation in transduced vs untransduced normalized.jpeg",type="jpeg")
savePlot("Proliferation in transduced vs untransduced Unnormalized.jpeg",type="jpeg")
lines(y~x)



var.x<-"Z.all.Green"
var.y<-"RG_over_G_Z"

var.x<-"Z.all.Green"
var.y<-"RnG_over_nG"

var.x<-"Z.all.Green"
var.y<-"Trans_percent_"

var.x<-"Z.all.Green"
var.y<-"G"
pos.control.list<-c("CCND1","MOCK","HSPCN111-N","UNTRANSDUCED","no DNA","PLV411") # controls that would affect mediann center normalization
norm.control.list<-c("HSPCN111-N")




plot(normalized[,var.x],normalized[,var.y],ylab=var.y,xlab=var.x,main="Z-score well based VS Z score plate based",cex.lab=1.5,cex.main=1.5)
sort(tapply(normalized[,"Symbol"],normalized[,"Symbol"],length))

pos.control.list<-c("CCNE1","101-CCNE1","411-CCNE1","101-CYCD2","not transduced","MOCK") # controls that would affect mediann center normalization
norm.control.list<-c("PLV101","plv101","PLV101-2","CPNE3","PSTPIP2")


pos.control.list<-c("CCNE1","101-CCNE1","411-CCNE1","101-CYCD2") # controls that would affect mediann center normalization
pos.control<-normalized[,"Symbol"] %in% pos.control.list
points(normalized[pos.control,var.x],normalized[pos.control,var.y],col="red",pch=19,cex=1.5 )
normalized[pos.control,]

norm.control.list<-c("HSPCN111-N")
norm.control.list<-c("PLV101","plv101","PLV101-2")
norm.control<-normalized[,"Symbol"] %in% norm.control.list
points(normalized[norm.control,var.x],normalized[norm.control,var.y],col="blue",pch=19,cex=0.8 )

hit.cut<-9.0
hits<-normalized[,"Z.all.Green"]>hit.cut
points(normalized[hits,var.x],normalized[hits,var.y],col="orange",pch=19,cex=1.0 )

mock.control.list<-c("MOCK")
mock.control<-normalized[,"Symbol"] %in% mock.control.list
points(normalized[mock.control,var.x],normalized[mock.control,var.y],col="green",pch=20,cex=0.7)

leg.txt<-c("CCND1","HSPCN111-N","MOCK")


legend(-20,25.1,leg.txt,col=c("red","blue","green"),pch=19,cex=2.0)

legend(-12,0.1,leg.txt,col=c("red","blue"),pch=19,cex=2.0)

abline(a=0,b=1,lwd=2)

text(11,0.01," <-Hit Cutoff",col="red",cex=2.0,pos=1)

selected.data<-identify(normalized[,var.x],normalized[,var.y],labels=paste(normalized[,"Symbol"],"",normalized[,"plate"],sep=":"),col="forestgreen",font=2,cex=0.8,offset=4,atpen='TRUE',plot=TRUE)

selected.data<-identify(normalized[,var.x],normalized[,var.y],labels=paste(normalized[,"Symbol"],sep=":"),col="forestgreen",font=2,cex=0.8,offset=4,atpen='TRUE',plot=TRUE)

savePlot("Duka.plot.jpeg",type="jpeg")
savePlot("Duka.plot.ZvsTransduction.jpeg",type="jpeg")
savePlot("Duka.plot.ZvsNumGreen.jpeg",type="jpeg")
savePlot("Proliferation in transduced vs untransduced normalized.jpeg",type="jpeg")
savePlot("Proliferation in transduced vs untransduced Unnormalized.jpeg",type="jpeg")
savePlot("Z score well vs Z score plate.jpeg",type="jpeg")

library(robust)
library(robustbase)
the.model<-lmrob(RG_over_G_Z~Z.all.Green,data=normalized)
summary(the.model)

plot((test[,xaxis]),(test[,yaxis]),pch=20,cex=0.1,xlab="Total Intensity of DAPI",ylab="EDU Intensity",main=paste(file.plate.name," Well:",row.label.letter,col.label.number,sep=" "),font.lab=2,cex.lab=1.5,cex.main=2.0,ylim=c(0,red.yrange.on.plot),xlim=c(0,max.ObjectTotalIntenCh1))

abline(coef(the.model),lty=10,col="red",lwd=2)

text(-5,5,"R^2=0.7 ; p-value: < 2e-16",cex=1.5,font=2)

     the.model2<-lm(RG_over_G_Z~Z.all.Green,data=normalized)


sd(normalized[(grepl("MUS.a",normalized[,"plate"]) & normalized[,"Symbol"]=="HSPCN111-N"),"RG_over_G"])
mean(normalized[(grepl("MUS.a",normalized[,"plate"]) & normalized[,"Symbol"]=="HSPCN111-N"),"RG_over_G"])

4*(0.25)^2/(0.16)^2
mean(normalized[(grepl("MUS.a",normalized[,"plate"]) & normalized[,"Symbol"]=="PLV411"),"RG_over_G"])
sd(normalized[(grepl("MUS.a",normalized[,"plate"]) & normalized[,"Symbol"]=="PLV411"),"RG_over_G"])


normalized[(grepl("MUS.a",normalized[,"plate"]) & normalized[,"Symbol"]=="HSPCN111-N"),1:4]

plates.all[(grepl("MUS.a",normalized[,"plate"]) & normalized[,"Symbol"]=="HSPCN111-N"),1:4]

tapply(plates.all[(grepl("MUS.a",plates.all[,"plate"]) & plates.all[,"Symbol"]=="HSPCN111-N"),"RG_over_G"],plates.all[(grepl("MUS.a",plates.all[,"plate"]) & plates.all[,"Symbol"]=="HSPCN111-N"),"plate"],mean)

tapply(plates.all[(grepl("MUS.a",plates.all[,"plate"]) & plates.all[,"Symbol"]=="HSPCN111-N"),"RG_over_G"],plates.all[(grepl("MUS.a",plates.all[,"plate"]) & plates.all[,"Symbol"]=="HSPCN111-N"),"plate"],sd)

sd(tapply(plates.all[(grepl("MUS.a",plates.all[,"plate"]) & plates.all[,"Symbol"]=="HSPCN111-N"),"RG_over_G"],plates.all[(grepl("MUS.a",plates.all[,"plate"]) & plates.all[,"Symbol"]=="HSPCN111-N"),"plate"],mean))

all<-1:1000
all2<-1:1000
for (i in 1:1000){
all[i]<-sd(rnorm(mean=0.23,sd=0.027,n=16))
all2[i]<-mean(rnorm(mean=0.23,sd=0.027,n=16))
}
mean(all)
max(all)
min(all)
3*sd(all)/0.027*100


max(all)/0.027
0.027/min(all)

0.027/max(all)
min(all)/0.027


0.23/max(all2)
min(all2)/0.23

max(all2)/0.23
0.23/min(all2)
