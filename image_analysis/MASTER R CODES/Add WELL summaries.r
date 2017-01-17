
#########################################
########################################
#Use this file to combine well based and cell based data typically done BEFORE normalization
############################################### Ca screen latest
setwd( "/media/Bioinform-D/Research/Cellomics/Ca screen/Latest")
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

exported.well.data<-"Well_based_data.txt"

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
well.output.file<-"Ca2_summary_NOTGREEN_DNA.txt"
well.type<-384
row.type<-16
col.type<-24

map.file<-"cellomics.to.expt.map.csv" #
exported.well.data<-"well_based_data.txt"

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

map.file<-"cellomics.to.expt.map.csv" #
exported.well.data<-"well_based_data.txt"
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
well.output.file<-"siFB_summary_NOTGREEN_DNA.txt"
well.type<-384
row.type<-16
col.type<-24
exported.well.data<-"well_based_data.txt"
#######################################################

############# specal cases when ak and rezaurin in barcode row col format:





plates.all<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized
normalized<-plates.all


###for flexaible annotion use this one
colnames(plates.all)<-c(core.ann,core.vars,colnames(ann)[!place.core.ann])## 

 ## sum.vars<-c(core.ann,"Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S")

         
###### if normalized exists use this one:
## normalized<-read.delim(normalized.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Normalized



map<-read.delim(map.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)


exported<-read.delim(exported.well.data,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
exported[1:5,]

########## rename export columns #NOt plate in data aligned to Barcode by default without a map file
exported<-exported[!is.na(exported[,"Row"]),]

tapply(normalized[,"plate"],normalized[,"plate"],length)
tapply(exported[,"BarCode"],exported[,"BarCode"],length)
map

#UPDs<-tapply(exported[,"UPD"],exported[,"UPD"],length)


###############if NOT map file and Barcodes(well) === plate(cell data) then use this only:
if(!exists("map")){
  colnames(exported)[colnames(exported)=="BarCode"]<-"Plate"
}else{  #### not jump tp $$$$

BarCode<-unique(exported[,"BarCode"])

posns<- match(BarCode,map[,"BarCode"])
missing<-is.na(posns)
sum(missing) # should be zero

################## should map the plate to the barcode ### not implemeted yet
BarCode<-map[,"BarCode"]

posns<-match(exported[,"BarCode"],BarCode)
missing<-is.na(posns)
sum(missing) ## should be zero
posns<-posns[!missing]
exported<-cbind(Plate=map[posns,"Plate"],exported) ## add plates names  to we
}



######################  $$$
exported[1:5,]
exported[,"Col"] <- exported[,"Col"]+1
exported[,"Row"] <- exported[,"Row"]+1
row.index<-c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
exported[,"Row"] <- row.index[exported[,"Row"]]
########remap
tapply(normalized[,"plate"],normalized[,"plate"],length)
tapply(exported[,"Plate"],exported[,"Plate"],length) ## these must match
unique(exported[,"Plate"])
unique(normalized[,"plate"])


order.wanted<-paste(normalized[,"plate"],normalized[,"row"],normalized[,"col"])

order.have<-paste(exported[,"Plate"] ,exported[,"Row"] ,exported[,"Col"] ) ## use this if have a map file
##order.have<-paste(exported[,"BarCode"] ,exported[,"Row"] ,exported[,"Col"] )

posns<-match(order.wanted,order.have)
present<-!is.na(posns)
sum(!present) # if zero  some exported data missing which can happen
## posns<-posns[present]

normalized<-cbind(normalized,exported[posns,])
## ori.size<-ncol(normalized)+1
## normalized<-cbind(normalized,matrix(data=NA,nrow=nrow(normalized),ncol=ncol(exported)))
## colnames(normalized)[ori.size:ncol(normalized)]<-colnames(exported)
## normalized[present,colnames(exported)]<-exported[posns,]

##checks
sum(normalized[,"row"] != normalized[,"Row"] & !is.na(normalized[,"Row"])) # must be zero
sum(normalized[,"col"] != normalized[,"Col"] & !is.na(normalized[,"Col"])) # must be zero
sum(normalized[,"plate"] != normalized[,"Plate"] & !is.na(normalized[,"Plate"])) 

sum(normalized[,"plate"] != normalized[,"BarCode"] & !is.na(normalized[,"BarCode"])) # MAYBE be zero is Barcode same as plate name
###### if writing normalized
write.table(normalized,paste(the.screen,"NOTGREEN2","Wells","txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t")

normalized[1:5,]
######################################## stop here for adding well data now add plate raeder data:
normalized.ori<-normalized
normalized<-normalized.ori


########################### Add plate reader data ##########################################
#############################################################################################
#############################################################################################
#############################################################################################

#################### barcode row col format
well.output.file<-paste(the.screen,"NOTGREEN2","Wells","txt",sep=".")
plates.all<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized
normalized<-plates.all


################# Choose one
exported.well.data<-"resazurin.txt"
exported.well.data<-"ak_leakage.txt"


## IF have a map file## 
map<-read.delim(map.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
######

exported<-read.delim(exported.well.data,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
exported[1:5,]


###############if NOT map file and Barcodes(well) === plate(cell data) then use this only:
if(!exists("map")){
  colnames(exported)[colnames(exported)=="BarCode"]<-"Plate"
}else{  #### not jump tp $$$$

BarCode<-unique(exported[,"BarCode"])

posns<- match(BarCode,map[,"BarCode"])
missing<-is.na(posns)
sum(missing) # should be zero

################## should map the plate to the barcode ### not implemeted yet
BarCode<-map[,"BarCode"]

posns<-match(exported[,"BarCode"],BarCode)
missing<-is.na(posns)
sum(missing) ## should be zero
posns<-posns[!missing]
exported<-cbind(Plate=map[posns,"Plate"],exported) ## add plates names  to we
}



######################  $$$
exported[1:5,]
## exported[,"Col"] <- exported[,"Col"]+1
## exported[,"Row"] <- exported[,"Row"]+1
## row.index<-c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
## exported[,"Row"] <- row.index[exported[,"Row"]]
########remap
tapply(normalized[,"plate"],normalized[,"plate"],length)
tapply(exported[,"Plate"],exported[,"Plate"],length) ## these must match
unique(exported[,"Plate"])
unique(normalized[,"plate"])

order.wanted<-paste(normalized[,"plate"],normalized[,"row"],normalized[,"col"])

order.have<-paste(exported[,"Plate"] ,exported[,"Row"] ,exported[,"Col"] ) ## use this if have a map file
##order.have<-paste(exported[,"BarCode"] ,exported[,"Row"] ,exported[,"Col"] )

posns<-match(order.wanted,order.have)
present<-!is.na(posns)
sum(!present) # if zero  some exported data missing which can happen
## posns<-posns[present]

normalized<-cbind(normalized,exported[posns,])
## ori.size<-ncol(normalized)+1
## normalized<-cbind(normalized,matrix(data=NA,nrow=nrow(normalized),ncol=ncol(exported)))
## colnames(normalized)[ori.size:ncol(normalized)]<-colnames(exported)
## normalized[present,colnames(exported)]<-exported[posns,]

##checks
sum(normalized[,"row"] != normalized[,"Row"] & !is.na(normalized[,"Row"])) # must be zero
sum(normalized[,"col"] != normalized[,"Col"] & !is.na(normalized[,"Col"])) # must be zero
sum(normalized[,"plate"] != normalized[,"Plate"] & !is.na(normalized[,"Plate"])) 

sum(normalized[,"plate"] != normalized[,"BarCode"] & !is.na(normalized[,"BarCode"])) # MAYBE be zero is Barcode same as plate name
###### if writing normalized
normalized[1:5,]
### go up and do next

write.table(normalized,paste(the.screen,"NOTGREEN3","Wells","txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t")


##########################################################################################################################################


###################### Grid formaT
######################below if is in grid formaT
###################assumes plate names are in 1st column and start with a ">" and no missing cols or rows ########
setwd("/media/Bioinform-D/Research/Cellomics/Ca screen/Ca screen 2")
row.type<-16
col.type<-24
files<-c("ak_leakage.txt","resazurin.txt")


############### test row nummbering
############### test row nummbering
unique(new.data[,"plate"])
unique(normalized[,"plate"])

row.ids<-row.index[1:row.type]
row.ids.t<-rep(row.ids,col.type)
dim(row.ids.t)<-c(length(row.ids),col.type)
row.ids.t<-t(row.ids.t)
row.ids.t<-as.vector(row.ids.t)
col.ids.t<-rep(1:col.type,times=row.type)

for(ifile in 1:length(files)){

label<-gsub(".txt","",files[ifile])


reader<-read.delim(files[ifile],header=F,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
apply(reader,1,function(x) sum(is.na(x)))
plate.starts<-grep("^>",reader[,1])
plate.names<-gsub(">","",reader[plate.starts,1])
plate.names
new.data<-{}
for(i in 1:length(plate.starts)){
 start<- (plate.starts[i]+1)
 data<-reader[start:(start+row.type-1),]
 data<-as.numeric(t(as.matrix(data)))
 data<-cbind(plate=plate.names[i],row=row.ids.t,col=col.ids.t,the.data=data)
new.data<-rbind(new.data,data)
}
colnames(new.data)[4]<-label
new.data[,"plate"]<-gsub("siCA2_","",new.data[,"plate"]) # casue max has a "siCA2_" ask him to remove



order.wanted<-paste(normalized[,"plate"],normalized[,"row"],normalized[,"col"])
order.have<-paste(new.data[,"plate"] ,new.data[,"row"] ,new.data[,"col"] )
posns<-match(order.wanted,order.have)
present<-!is.na(posns)
posns<-posns[present]

######## remap and export
## rownames(norm.ori)<-paste(norm.ori[,"plate"],norm.ori[,"row"],norm.ori[,"col"])
## rownames(normalized)<-paste(normalized[,"plate"],normalized[,"row"],normalized[,"col"])
## rownames(exported)<-paste(exported[,"Plate"] ,exported[,"Row"] ,exported[,"Col"] )

ori.size<-ncol(normalized)+1
normalized<-cbind(normalized,matrix(data=NA,nrow=nrow(normalized),ncol=1))
colnames(normalized)[ori.size:ncol(normalized)]<-label
normalized[present,label]<-new.data[posns,label]
} #loop over files


########################################################################
###########special case ofe ca3
new.data<-read.delim("siCA3 plate reader data.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
new.data<-new.data[,c( "Plate.ID", "row", "col", "AK.Raw.Data", "Raw.Data", "Symbol")]
colnames(new.data)<-c( "plate", "row", "col", "ak_leakage", "resazurin", "Symbol")
new.data[,"plate"]<-gsub(">siCA3_","",new.data[,"plate"])
unique(new.data[,"plate"])
unique(normalized[,"plate"])
####

order.wanted<-paste(normalized[,"plate"],normalized[,"row"],normalized[,"col"])
order.have<-paste(new.data[,"plate"] ,new.data[,"row"] ,new.data[,"col"] )
posns<-match(order.wanted,order.have)
 sum(is.na(posns)
present<-!is.na(posns)
posns<-posns[present]

label<-"ak_leakage"
     label<-"resazurin"
 ori.size<-ncol(normalized)+1
normalized<-cbind(normalized,matrix(data=NA,nrow=nrow(normalized),ncol=1))
colnames(normalized)[ori.size:ncol(normalized)]<-label
normalized[present,label]<-new.data[posns,label]
##############################################################

write.table(normalized,paste(the.screen,"NOTGREEN3","Wells","txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t")
