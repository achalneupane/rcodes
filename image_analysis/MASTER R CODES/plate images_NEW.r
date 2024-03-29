




target<-"X.R"
plot.new()
temp2<-matrix(0:10)  #temp2<-matrix(1:10,4,5) paste(expression(R^2),"for","rs999999",sep=" ")
image(temp2, col = colorss,xlab=bquote(bold(.(target))),main=expression(bold(Scale)),axes=FALSE,cex=1,cex.lab=2.5,cex.main=2,lwd=3)
### NOTE USE APPLICATION -> CHARACTER MAP TO GET UNICODE ID
axis(1,at=seq(0,1,0.1),labels=as.character(seq(0,1,0.1)),cex.axis=1.25,font=2 )
box()
the.colors<-colors()
image(x=1:65,y=1:10,matrix(1:650,65,10), col = the.colors[1:650])
text(1,10,labels="x")
axis(1,at=seq(0,1,0.1),labels=as.character(seq(0,1,0.1)),cex.axis=1.25,font=2 )


the.symbols<-matrix(c(1:650),nrow=65,ncol=10,byrow=TRUE)
the.symbols<-t(the.symbols)[,ncol(t(the.symbols)):1]
the.symbols<-as.character(the.symbols)
text(c(1:65)/0.1,c(1:10)/0.1,labels=the.symbols,cex.axis=0.5,font=2 ,las=2)



well.output.file<-"Kiril_summary_NOTGREEN_DNA.txt"
            
###START#########################
setwd( "/media/Bioinform-D/Research/Cellomics/Kiril screen")
sum.vars<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","Description","Horfeome.position", "MGC.Accession","Pilot.Position","Gene.Name.pilotPos","Kiril.position")

well.output.file<-"Kiril_summary_NOTGREEN_DNA.txt"
well.type<-96
row.type<-8
col.type<-12

plates.all<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized

the.screen="Kiril"
normalized<-read.delim(paste(the.screen,"NORMALIZED","txt",sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
plates.all<-normalized


###START#########################
setwd( "/media/Bioinform-D/Research/Cellomics/Joseph screen")
sum.vars<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","Description","GenBank_Accn","dbEST_Id" ,"MGC.Accession")

well.output.file<-"Joseph_summary_NOTGREEN_DNA.txt"
well.type<-96
row.type<-8
col.type<-12

plates.all<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized

the.screen="Joseph"
normalized<-read.delim(paste(the.screen,"NORMALIZED","txt",sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

############################################

############################################
############################################### Inhibitor screen
setwd( "/media/Bioinform-D/Research/Cellomics/inhibitor")

sum.vars<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","Description","Horfeome.position", "MGC.Accession","Pilot.Position","Gene.Name.pilotPos","Kiril.position")

well.output.file<-"Inhibitor_summary_NOTGREEN_DNA.txt"
well.type<-96
row.type<-8
col.type<-12
the.screen<-"Inhibitor"

#normalized<-read.delim(paste(the.screen,"NORMALIZED","txt",sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
normalized<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

#######################################################
############################################### Ca screen
setwd( "/media/Bioinform-D/Research/Cellomics/Ca screen")
sum.vars<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S","Entrez_ID","Cell_loc","Description","Gene.ID","Accession.Number","GI.Number","Library.Type")


field.output.file<-"Ca_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Ca_summary_NOTGREEN_DNA.txt"

normalized.file<-paste(the.screen,"NORMALIZED2","txt",sep=".")
normalized.file<-paste(the.screen,"NORMALIZED2","Wells","txt",sep=".")
well.type<-384
row.type<-16
col.type<-24
the.screen<-"Ca siRNA"

## normalized<-read.delim(paste(the.screen,"NORMALIZED","txt",sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)



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
well.type<-96
row.type<-8
col.type<-12
#######################################################



############################################### Ca screen latest
setwd( "/media/Bioinform-D/Research/Cellomics/Ca screen/Latest")
#load("Ca_ann2.RData") # same as Ca_ann2 except 
map.file<-"cellomics.to.expt.map.csv" # "plate" must be same "plate in annoation AND IN THE SAME ORDER
annotations.file<-"Ca_annotations2.txt"

sum.vars<-c("plate","row","col","Symbol","%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S","Entrez_ID","Cell_loc","Description","Gene.ID","Accession.Number","GI.Number","Library.Type")

the.screen<-"Ca siRNA2"
#files<- paste("Ca_siCa_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="") # defined below
field.output.file<-"Ca siRNA2_summary_NOTGREEN_DNA.txt"
well.output.file<-"Ca siRNA2.NOTGREEN3.Wells.txt"
well.type<-384
row.type<-16
col.type<-24

exported.well.data<-"Well_based_data.txt"
normalized.file<-paste(the.screen,"NORMALIZED","txt",sep=".")

#######################################################

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
normalized.file<-paste(the.screen,"NORMALIZED","txt",sep=".")
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
normalized.file<-paste(the.screen,"NORMALIZED","txt",sep=".")
#######################################################

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

normalized.file<-paste(the.screen,"NORMALIZED","txt",sep=".")


targets<-c("S+G2_over_G1","percent_R", "ak_leakage","resazurin","Num-cells")

##"S+G2_over_G1" > 40 look crappy
load("a.screen.cuts.RData")
breakup<-TRUE
#######################################################

targets<-c("S+G2_over_G1","percent_R","resazurin","Num-cells")

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

targets<-c("Num-cells","Z-all Green","Relative S+G2+_gt_4N","Gr-S+G2+ _gt_4N_over_G1","RG_over_G","RG_over_G_Z"  )
breakup<-TRUE

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

targets<-c("Num-cells","Z-all Green","Relative S+G2+_gt_4N","Gr-S+G2+ _gt_4N_over_G1","RG_over_G","RG_over_G_Z"  )
breakup<-TRUE
#######################################################



if(breakup){  ### only use breakup if multiple screens are in the same file

normalized<-read.delim(normalized.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
the.order<-order(normalized[,"plate"])
normalized<-normalized[the.order,]


col<-colnames(normalized) %in% core.ann
labels<-c(colnames(normalized)[col],core.vars,colnames(ann)[!place.core.ann])## 
colnames(normalized)[1:length(labels)]<-labels
##############################################################################################################
##############################################################################################################
##############################################################################################################
############################ NEED TO RESTART FROM HERE EVERY TIME AS PLATES.ALL IS MODIFIED ##############################

plates.all<-normalized
labels<-gsub("%","percent_",colnames(plates.all))
labels<-gsub("/","_over_",labels)
labels<-gsub(">","_gt_",labels)
colnames(plates.all)<-labels
colnames(plates.all)<-gsub("/","_",colnames(plates.all))

norm.all<-plates.all



########################RESTART
norm.all[1:5,]
dim(norm.all)
norm<-norm.all


## #####
## the.cells<-"HACK"
## the.screen<-"si FB1"
## wanted.screen<-grepl(the.screen,norm.all[,"plate"])
## sum(wanted.screen)
## norm<-norm.all[wanted.screen,]
## unique(norm[,"plate"])
## low.cut<-a.screen.cuts[[the.screen]]$low.cut
## high.cut<-a.screen.cuts[[the.screen]]$high.cut
## dim(norm)
## #####
## #####
## the.cells<-"C33A"
## the.screen<-"si FB2"
## wanted.screen<-grepl(the.screen,norm.all[,"plate"])
## sum(wanted.screen)
## norm<-norm.all[wanted.screen,]
## unique(norm[,"plate"])
## low.cut<-a.screen.cuts[[the.screen]]$low.cut
## high.cut<-a.screen.cuts[[the.screen]]$high.cut
## dim(norm)
## #####
## #####
## the.cells<-"CASKI"
## the.screen<-"si FB3"
## wanted.screen<-grepl(the.screen,norm.all[,"plate"])
## sum(wanted.screen)
## norm<-norm.all[wanted.screen,]
## unique(norm[,"plate"])
## low.cut<-a.screen.cuts[[the.screen]]$low.cut
## high.cut<-a.screen.cuts[[the.screen]]$high.cut
## dim(norm)
## #####


## the.screen<-"si FB1"
the.cells<-"EtOH"
the.screen<-"a"
wanted.screen<-grepl(the.screen,norm[,"plate"])
sum(wanted.screen)
norm<-norm[wanted.screen,]
unique(norm[,"plate"])
#####

#####
## the.cells<-"C33A"
## the.screen<-"si FB2"
the.cells<-"ICI"
the.screen<-"b"
wanted.screen<-grepl(the.screen,norm[,"plate"])
sum(wanted.screen)
norm<-norm[wanted.screen,]
unique(norm[,"plate"])
#####

#####
the.cells<-"EX"
the.screen<-"c"
wanted.screen<-grepl(the.screen,norm[,"plate"])
sum(wanted.screen)
norm<-norm[wanted.screen,]
unique(norm[,"plate"])
#####




dim(norm)
norm[1:5,]
}

plot(1:100) ## to resize 
#################################################################################
#normalized.file<-paste(the.screen,"NORMALIZED","txt",sep=".")
## targets<-c("Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","%R","Num-cells","MEAN_Area","MEAN_LWR","MEAN_P2A","MEAN_Total_I","SD_LWR","SD_P2A")
## targets<-c("%R","S+G2_G1","S","ak_leakage","resazurin","Cell_count","MEAN_Area","MEAN_LWR","MEAN_P2A","MEAN_Avg_Inten","MEAN_Total_Inten","SD_LWR","SD_P2A")
## targets<-c("G1 posn","S+G2_over_G1","percent_R", "ak_leakage","resazurin","Num-cells"  )

on.sheet<-4

for(itarget in 1:length(targets)){
target<-targets[itarget]

###### for old annotation
## plates.all<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized
## colnames(plates.all)<-sum.vars
## normalized<-read.delim(normalized.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
## colnames(normalized)[1:length(sum.vars)]<-sum.vars


##### flexible annotation
if(breakup){plates.all<-norm}else{
normalized<-read.delim(normalized.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
the.order<-order(normalized[,"plate"])
normalized<-normalized[the.order,]


col<-colnames(normalized) %in% core.ann
labels<-c(colnames(normalized)[col],core.vars,colnames(ann)[!place.core.ann])## 
colnames(normalized)[1:length(labels)]<-labels
##############################################################################################################
##############################################################################################################
##############################################################################################################
############################ NEED TO RESTART FROM HERE EVERY TIME AS PLATES.ALL IS MODIFIED ##############################

plates.all<-normalized
labels<-gsub("%","percent_",colnames(plates.all))
labels<-gsub("/","_over_",labels)
labels<-gsub(">","_gt_",labels)
colnames(plates.all)<-labels
colnames(plates.all)<-gsub("/","_",colnames(plates.all))
}
## on.sheet<-8
## nf<-layout(matrix(c(1,1,2:9),5,2,byrow=TRUE),heights=c(0.35,1,1,1,1))


## nf<-layout(matrix(c(1,2:4),4,1,byrow=TRUE),heights=c(0.35,1,1,1,1)) ##this for on.sheet=3 only
nf<-layout(matrix(c(1,2:(on.sheet+1)),on.sheet+1,1,byrow=TRUE),heights=c(0.35,rep(1,times=on.sheet)))
layout.show(nf) # for use in mon-floating scale bar
#temp2<-matrix(1:10,4,5) paste(expression(R^2),"for","rs999999",sep=" ")



#plates.all<-plates.all[plates.all[,"plate"]==1,]
##############
##############
#color = colorRampPalette(c("green","red","blue"),space="rgb",interpolate="linear") # Lab/rgb spline/linear
## color = colorRampPalette(c("forestGreen","white","red"),space="rgb",interpolate="linear")
## color = colorRampPalette(c("forestGreen","yellow","white","pink","purple3"),space="rgb",interpolate="linear")
## color = colorRampPalette(c("yellow","white","red"),space="rgb",interpolate="linear")

## color = colorRampPalette(c("forestGreen","yellowgreen","white","orange","red3"),space="rgb",interpolate="linear")
## color = colorRampPalette(c("blue","lightseagreen","white","orange","purple3"),space="rgb",interpolate="linear")
## color = colorRampPalette(c("turquoise4","yellowgreen","white","orange","purple3"),space="rgb",interpolate="linear")
color = colorRampPalette(c("darkgreen","yellowgreen","white","orange","purple3"),space="rgb",interpolate="linear")
##color = colorRampPalette(c("darkgreen","gold1","white","orange","purple3"),space="rgb",interpolate="linear")
colorss<-color(11) #  color = colorRampPalette(c("blue","red"))

## columns.alpha<-c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
## names(columns.alpha)<-as.character(c(1:length(columns.alpha)))

columns.alpha<-c(1:row.type)
names(columns.alpha)<-c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")[1:length(columns.alpha)]
#plates.all<-read.delim(the.file,header=T,sep="\t",fill=TRUE)
#fields<-read.delim("Joseph_field_summary_NOTGREEN.txt",header=T,sep="\t",fill=TRUE)
the.expt<-the.screen
#colnames(plates.all)[c(11,26,27)]<-c("%R","Trans %","Mean_cells_field")



max.val<-max(plates.all[is.finite(plates.all[,target]),target])
min.val<-min(plates.all[is.finite(plates.all[,target]),target])

the.median<-median(plates.all[is.finite(plates.all[,target]),target])
the.sd<-sd(plates.all[is.finite(plates.all[,target]),target])

max.val2<-the.median+ 4*the.sd
min.val2<-the.median- 4*the.sd

if(max.val > max.val2){max.val<-max.val2}
if(min.val <  min.val2){min.val<-min.val2}

plates.all[(is.infinite(plates.all[,target]) | plates.all[,target] > max.val ) & !is.na(plates.all[,target])  ,target] <-max.val
plates.all[(is.infinite(plates.all[,target]) | plates.all[,target] < min.val ) & !is.na(plates.all[,target])  ,target] <-min.val

## plates.all[(is.infinite(plates.all[,target]) | plates.all[,target] < min.val ),target][408:413]
## plates.all[,target][408:413]
## plates.all[(is.infinite(plates.all[,target])),target][408:413]
## is.infinite(plates.all[,target])[408:413]
## (plates.all[,target] < min.val)[408:413]
## (is.infinite(plates.all[,target]) | plates.all[,target] < min.val )[408:413]

plates.all[plates.all[,target]==min.val & !is.na(plates.all[,target]),]
plates.all[plates.all[,target]==max.val & !is.na(plates.all[,target]),]
plates.all[is.infinite(plates.all[,target]),]
breaks<-(max.val-min.val)/length(colorss)

if(min.val < 0){
col.index<-floor((plates.all[,target]+abs(min.val))/breaks)
}else{
  col.index<-floor((plates.all[,target]-abs(min.val))/breaks)
}
col.index[grep(length(colorss),col.index)]<-col.index[grep(length(colorss),col.index)]-1 # change highest value
col.index<-col.index+1 # index not starting at zero

col=colorss[col.index]
col[is.na(col)]<-"gray40"  ##q
plates.all<-cbind(plates.all,col)
colnames(plates.all)[dim(plates.all)[2]]<-"colors"
            
par.opts.start<-par(no.readonly=TRUE)

  ############# DO BIOMART QUERYIES FIRST
layout.show(nf)
par(mar=c(2.0,12.5,1.5,12.5),mgp=c(0.55,0.55,0)) #c(bottom, left, top, right)
image(matrix(0:10), col = colorss,xlab=bquote(bold(.(target))),main=expression(bold(Scale)),axes=FALSE,cex=1,cex.lab=1.75,cex.main=1.75,lwd=3)
text(seq(0,1,0.1),0.1,labels=as.character(round(seq(min.val+breaks/2,max.val,breaks),2)),cex=2.0,font=2 )
box()

plate.list<-unique(plates.all[,"plate"])

par(mar=c(1.5,2.0,1.0,1),mgp=c(0.5,0.35,0)) #c(bottom, left, top, right)
   
plates.all[1:5,]
target

for(i in 1:length(plate.list)){

### NOTE USE APPLICATION -> CHARACTER MAP TO GET UNICODE ID
a.plate<-plates.all[plates.all[,"plate"]==plate.list[i],target]
the.symbols<-plates.all[plates.all[,"plate"]==plate.list[i],"Symbol"]
rows<-plates.all[plates.all[,"plate"]==plate.list[i],"row"]
cols<-plates.all[plates.all[,"plate"]==plate.list[i],"col"]
the.colors<-as.character(plates.all[plates.all[,"plate"]==plate.list[i],"colors"])
index<-paste(rows,cols,sep="")
a.plate<-matrix(a.plate,nrow=length(unique(rows)),ncol=length(unique(cols)),byrow=TRUE)
the.colors<-matrix(the.colors,nrow=length(unique(rows)),ncol=length(unique(cols)),byrow=TRUE)
the.symbols<-matrix(the.symbols,nrow=length(unique(rows)),ncol=length(unique(cols)),byrow=TRUE)
### images start from bottom row fills next row up same column row-wise 
the.colors<-t(the.colors)[,ncol(t(the.colors)):1]
the.colors<-as.character(the.colors)

the.symbols<-t(the.symbols)[,ncol(t(the.symbols)):1]
the.symbols<-as.character(the.symbols)

#plot.new()
image(x=1:ncol(a.plate),y=1:nrow(a.plate),matrix(1:well.type,col.type,row.type), col = the.colors,xlab="",ylab="",main=bquote(bold("Plate :"~.(plate.list[i]))),axes=FALSE,cex=1,cex.lab=2.5,cex.main=1.25,lwd=3,font=2)
### NOTE USE APPLICATION -> CHARACTER MAP TO GET UNICODE ID
axis(1,at=seq(1,ncol(a.plate),1),labels=as.character(seq(1,ncol(a.plate),1)),cex.axis=1.0,font=2 )
axis(2,at=seq(1,nrow(a.plate),1),labels=as.character(names(columns.alpha[seq(nrow(a.plate),1,-1)])),cex.axis=1.0,font=2 ,las=2)
text(cols,as.numeric((columns.alpha[rows])),labels=as.character(the.symbols),cex.axis=1.0,font=2 ,las=2)
#text(3,1,labels=as.character(the.symbols[1]))

box()
if(floor(i/on.sheet)==ceiling(i/on.sheet)){
 savePlot(paste("Plate Images-",the.expt,"-",target,"-",i-on.sheet,"-",i,".tiff",sep=""),type="tiff")
 savePlot(paste("Plate Images-",the.expt,"-",target,"-",i-on.sheet,"-",i,".jpeg",sep=""),type="jpeg") 
layout.show(nf)
par(mar=c(2.0,12.5,1.5,12.5),mgp=c(0.55,0.55,0)) #c(bottom, left, top, right)
image(matrix(0:10), col = colorss,xlab=bquote(bold(.(target))),main=expression(bold(Scale)),axes=FALSE,cex=1,cex.lab=1.5,cex.main=1.5,lwd=3)
text(seq(0,1,0.1),0.1,labels=as.character(round(seq(min.val+breaks/2,max.val,breaks),2)),cex=2.0,font=2)
box()
par(mar=c(1.5,1.5,1.0,1),mgp=c(0.5,0.35,0))
                            }
}

} # over targets
    
 #savePlot(paste("Plate Images-",the.expt,"-",target,"-",i-on.sheet,"-",i,".jpeg",sep=""),type="jpeg") 




