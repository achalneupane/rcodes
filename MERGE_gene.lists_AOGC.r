
current.dir<-getwd()

ann.order<-c("refGene","knownGene","ensGene") # a.indel[1:5,paste(ann.order,"gene",sep="::")]
gene.names<-get.gene.names(a.indel,paste(ann.order,"gene",sep="::"))
anno.DB.location.core<-"/media/Bioinform-D/Research/annovar/humandb"






a.bone.file<-read.delim(paste(anno.DB.location.core,"/", "BoneDensity_GWAS.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
BoneDensity_GWAS<-as.data.frame(a.bone.file[,2]) #
colnames(BoneDensity_GWAS)<-colnames(a.bone.file)[2] # unique(AML_ALL_Mutations)
## colnames(AML_publications)<-"AML_publications"


a.bone.file<-read.delim(paste(anno.DB.location.core,"/", "BoneDensity_Sequencing.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
BoneDensity_Sequencing<-as.data.frame(a.bone.file[,2]) #
colnames(BoneDensity_Sequencing)<-colnames(a.bone.file)[2] # unique(AML_ALL_Mutations)
a.bone.file[!is.na(a.bone.file[,2]),]
## colnames(AML_MDS_publications)<-"AML_MDS_publications" # unique(AML_ALL_Mutations)


a.bone.file<-read.delim(paste(anno.DB.location.core,"/","omim.txt",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1,]
the.colnames<-colnames(a.bone.file)
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
a.bone.file<-apply(a.bone.file,1,function(x) if(x[1]=="NA" | is.na(x[1])){NA}else{paste(x,collapse="::")})
a.bone.file[1:5]
omim<-as.data.frame(a.bone.file)
colnames(omim)<-paste("OMIM (",paste(the.colnames,collapse="::"),")",sep="")
omim[1:5,]


#if(!exists("ingenuity.bone.genes")){
a.bone.file<-read.delim(paste(anno.DB.location.core,"/","bone.function.list.ingenuity.txt",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit")
ingenuity.bone.genes<-as.data.frame(a.bone.file[,2])
colnames(ingenuity.bone.genes)<-c("ingenuity.bone.genes")
#}

#if(!exists("Dequeant.cycling")){
a.bone.file<-read.delim(paste(anno.DB.location.core,"/","Dequeant.cycling.txt",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,"Symbol",c("Symbol","Entrez.Gene.Name"),"delimit")
Dequeant.cycling<-as.data.frame(a.bone.file[,2])
colnames(Dequeant.cycling)<-c("Dequeant.cycling")
#}

#if(!exists("mouse.defect")){
a.bone.file<-read.delim(paste(anno.DB.location.core,"/","mouse.defect.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
mouse.defect<-as.data.frame(a.bone.file[,1])
colnames(mouse.defect)<-c("mouse.defect")
#}

#if(!exists("ProtonPump")){
a.bone.file<-read.delim(paste(anno.DB.location.core,"/","proton.all.ingenuity.txt",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
ProtonPump<-as.data.frame(a.bone.file[,1])
colnames(ProtonPump)<-c("ProtonPump")
#}


#if(!exists("sewell.cycling")){
a.bone.file<-read.delim(paste(anno.DB.location.core,"/","sewell.cycling.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
sewell.cycling<-as.data.frame(a.bone.file[,1])
colnames(sewell.cycling)<-c("sewell.cycling")
#}
  

#if(!exists("skeletome")){
a.bone.file<-read.delim(paste(anno.DB.location.core,"/","skeletome.csv",sep=""),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
a.bone.file[1:5,]
a.bone.file<-match.cols.and.collapse(gene.names,"refGene::gene",a.bone.file,colnames(a.bone.file)[1],colnames(a.bone.file),"delimit") # macth column must in be collapse list
a.bone.file<-apply(a.bone.file,1,function(x) if(x[1]=="NA" | is.na(x[1])){NA}else{paste(x,collapse="::")})
skeletome<-as.data.frame(a.bone.file)
colnames(skeletome)<-c("skeletome")
#}



  

   save(list=c("BoneDensity_GWAS","BoneDensity_Sequencing","omim","ingenuity.bone.genes","Dequeant.cycling","mouse.defect","sewell.cycling","skeletome"),file=paste(gsub(".txt$","",target),".GeneLists.RData",sep=""))





######################## Here OMIM table shoudl always be present to gene.table is never empty
## target<-gsub(".TGCM-AML.analysis-maf-filtered.txt","",project.files[ichr])
## get<-try(load(paste(gsub(".txt$","",target),".GeneLists.RData",sep="")),silent=TRUE)   
gene.lists<-c("omim","ALS","skeletome","mouse.defect","sewell.cycling","Dequeant.cycling","ingenuity.bone.genes","BoneDensity_GWAS","BoneDensity_Sequencing")

  
## gene.lists<-c("BoneDensity_GWAS","BoneDensity_Sequencing","AML_Validated","AML_SNP_ALL_Confident","AML_SNP_Diagnosis","AML_SNP_Relapse","AML_DINDEL_Diagnosis","AML_DINDEL_Relapse","AML_CLINICAL_Mutations","AML_FRANC_Mutations","AML_ENU_Mutations","AML_publications","AML_MDS_publications","AML_SNP_Mutations","AML_ALL_Mutations", "AML_OTHER_Mutations","AML_MITO_Mutations","AML_LEY_RELAPSE_Mutations","AML_ASH_Mutations")

current.objects<-ls()
gene.lists<-gene.lists[gene.lists %in% current.objects]
gene.table<-{}
for(it in 1:length(gene.lists)){
  print(it)
  if(!exists(gene.lists[it])){print("missing")}
  a.temp<-eval(as.name(gene.lists[it]))
   print(paste(gene.lists[it],"-->",colnames(a.temp),sep=""))
  if(is.null(dim(gene.table)) ){gene.table<-a.temp}else{ if( dim(a.temp)[1]==dim(gene.table)[1] ){gene.table<-cbind(gene.table,a.temp)}else{print("dimension mismatch")} }
                                  }
rm("a.temp")
if(dim(gene.table)[1]!=dim(a.indel)[1]){print("ERROR gene.table DIFFERENT number of rows to a.indel");gene.table<-""}
dim(gene.table)
colnames(gene.table)
gene.table[1:5,]
#####################################################
setwd(current.dir)
