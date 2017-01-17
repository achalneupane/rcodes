
#########################################
########################################
Use this file to combine well based and cell based data typicall done after normalization


############################################### Ca screen
setwd( "/media/Bioinform-D/Research/Cellomics/Ca screen")
load("Ca_ann.RData")

sum.vars<-c("plate","row","col","Symbol","Trans_%","%R","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S","Entrez_ID","Cell_loc","Description","Gene.ID","Accession.Number","GI.Number","Library.Type")

the.screen<-"Ca siRNA"

############ dump excell files to tab delimed files 
files<- paste("Ca_siCa_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Ca_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Ca_summary_NOTGREEN_DNA.txt"
normalized.file<-paste(the.screen,"NORMALIZED2","txt",sep=".")
well.type<-384
row.type<-16
col.type<-24

map.file<-"cellomics.to.expt.map.csv"
exported.well.data<-"Well_based_data.csv"
##################### reda in file#####################

plates.all<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized
normalized<-plates.all

###### if normalized exists use this one:
normalized<-read.delim(normalized.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Normalized



map<-read.delim(map.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) 
exported<-read.delim(exported.well.data,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)


########## rename export columns # remove clank lines:
exported<-exported[!is.na(exported[,"Row"]),]

tapply(exported[,"UPD"],exported[,"UPD"],length)
#UPDs<-tapply(exported[,"UPD"],exported[,"UPD"],length)
UPDs<-unique(exported[,"UPD"])
posns<- match(exported[,"UPD"],UPDs)
exported<-cbind(Plate=map[posns,"Plate"],exported)

exported[,"Col"] <- exported[,"Col"]+1
exported[,"Row"] <- exported[,"Row"]+1
row.index<-c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
exported[,"Row"] <- row.index[exported[,"Row"]]
########remap
## order.wanted<-paste(ann[,"plate"],ann[,"row"],ann[,"col"])
order.wanted<-paste(normalized[,"plate"],normalized[,"row"],normalized[,"col"])
rownames(ann)<-order.wanted
order.have<-paste(exported[,"Plate"] ,exported[,"Row"] ,exported[,"Col"] )
posns<-match(order.wanted,order.have)
present<-!is.na(posns)
posns<-posns[present]
test<-exported[posns,]
######## remap and export
## rownames(norm.ori)<-paste(norm.ori[,"plate"],norm.ori[,"row"],norm.ori[,"col"])
## rownames(normalized)<-paste(normalized[,"plate"],normalized[,"row"],normalized[,"col"])
## rownames(exported)<-paste(exported[,"Plate"] ,exported[,"Row"] ,exported[,"Col"] )

ori.size<-ncol(normalized)+1
normalized<-cbind(normalized,matrix(data=NA,nrow=nrow(normalized),ncol=ncol(exported)))
colnames(normalized)[ori.size:ncol(normalized)]<-colnames(exported)
normalized[present,colnames(exported)]<-exported[posns,]
## normalized[rownames(exported),colnames(exported)]<-exported # goes with rownames method

## normalized[1:100,c(1:3,57,60:63)]
## sum(is.na(normalized[,57]))
## normalized[is.na(normalized[,57]),c(1:3,57,60:63)]

###### if writing normalized
write.table(normalized,paste(the.screen,"NotGreen2","Wells","txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t")



###################assumes plate names are in 1st column and start with a ">" and no missing cols or rows ########
setwd("/media/scratch/Data/Cellomics/Ca-screen-latest")
row.type<-16
col.type<-24
file<-"resazurin.txt"
file<-"ak_leakage.txt"
############### test row nummbering
############### test row nummbering
row.ids<-row.index[1:row.type]
row.ids.t<-rep(row.ids,col.type)
dim(row.ids.t)<-c(length(row.ids),col.type)
row.ids.t<-t(row.ids.t)
row.ids.t<-as.vector(row.ids.t)
row.ids.t<-rep(row.ids.t,length(plate.ids))
col.ids.t<-rep(1:row.type,times=col.type)
label<-gsub(".txt","",file)


reader<-read.delim(file,header=F,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
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


 order.wanted<-paste(normalized[,"plate"],normalized[,"row"],normalized[,"col"])
rownames(ann)<-order.wanted
order.have<-paste(exported[,"Plate"] ,exported[,"Row"] ,exported[,"Col"] )
posns<-match(order.wanted,order.have)
present<-!is.na(posns)
posns<-posns[present]
test<-exported[posns,]
######## remap and export
## rownames(norm.ori)<-paste(norm.ori[,"plate"],norm.ori[,"row"],norm.ori[,"col"])
## rownames(normalized)<-paste(normalized[,"plate"],normalized[,"row"],normalized[,"col"])
## rownames(exported)<-paste(exported[,"Plate"] ,exported[,"Row"] ,exported[,"Col"] )

ori.size<-ncol(normalized)+1
normalized<-cbind(normalized,matrix(data=NA,nrow=nrow(normalized),ncol=ncol(exported)))
colnames(normalized)[ori.size:ncol(normalized)]<-colnames(exported)
normalized[present,colnames(exported)]<-exported[posns,]
