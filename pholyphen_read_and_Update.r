
setwd(annotate.dir)
print("WOWSER in pholy...")
##########################################IF no old Phophen then will make a bland filetr table .Pholy
if(grepl("^chr",indels[1,"chr"])){
key.pholy<-paste(indels[,"chr"],":",indels[,"start"],".",indels[,"REF"],indels[,"ALT"],sep="")}else{key.pholy<-paste("chr",indels[,"chr"],":",indels[,"start"],".",indels[,"REF"],indels[,"ALT"],sep="")}
print("a")
filter.table.pholy<-data.frame(key=key.pholy,stringsAsFactors=FALSE)
polyphen.dir # <-"/media/ga-apps/UQCCG/Data/Sequence_Genotypes/2012-01-27_AOGC/PolyPhen"
polyphen.types<-c("Mendelian","GWAS") # polyphen.types<-c("Mendelian","Complex")
polyphen.files<-try(dir(polyphen.dir),silent = TRUE)
polyphen.files
polyphen.types
#if( sum(grepl(polyphen.types[1],polyphen.files)) == 0  & sum(grepl(polyphen.types[2],polyphen.files)) == 0 ){ polyphen.types<-project} # old type with only one predict type
#if( (!grepl(paste("^",polyphen.types[1],sep=""),polyphen.files)) | (!grepl(paste("^",polyphen.types[2],sep=""),polyphen.files)) ){ polyphen.types<-project} # old type with only one predict type
print("b")
i <- 1
for( i in 1:length(polyphen.types)){
#  the.poly.file<-polyphen.files[grepl(paste("^",polyphen.types[i],sep=""),polyphen.files)]
  the.poly.file<-polyphen.files[grepl(polyphen.types[i],polyphen.files)]
  if(length(the.poly.file)==0){
    #the.poly.file<-"dummy.txt"}
#x<-try(pholy<-read.table(paste(polyphen.dir,the.poly.file,sep="/"),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE,comment.char = ""),silent = TRUE)  

#if(inherits(x, "try-error")){
a.filter.table.pholy<-matrix(data=NA,nrow=length(key.indels),ncol=length(c("prediction","pph2_prob")))
colnames(a.filter.table.pholy)<-paste(polyphen.types[i],c("prediction","pph2_prob"),sep=".") # c("snp_id","pp2::Predict","pp2::score")
#
filter.table.pholy<-cbind(filter.table.pholy,a.filter.table.pholy) 
}else{
## pholy[1:5, "pph2_prob"]
## pholy[1:5, "prediction"]

pholy[,1]<-gsub("\\s","",pholy[, 1])
pholy[,"prediction"]<-gsub("\\s","",pholy[, "prediction"])
pholy[,"prediction"]<-gsub("damaging$"," damaging",pholy[, "prediction"])

## tapply(pholy[, "prediction"],pholy[, "prediction"],length)

posns<-match(key.pholy,pholy[,1])
missing<-is.na(posns)
sum(!missing)
if(i==1){filter.table.pholy[!missing,]<-pholy[posns[!missing], c("snp_id")] } #c("snp_id","prediction","pph2_class","pph2_prob","pph2_FPR","pph2_TPR","pph2_FDR"

a.filter.table.pholy<-pholy[posns, c("prediction","pph2_prob")]
colnames(a.filter.table.pholy)<-paste(polyphen.types[i],colnames(a.filter.table.pholy),sep=".")
  
filter.table.pholy<-cbind(filter.table.pholy,a.filter.table.pholy) #c("snp_id","prediction","pph2_class","pph2_prob","pph2_FPR","pph2_TPR","pph2_FDR"


# filter.table.pholy[1:20,]
}
}
print("c")
###################################### a.filter.table.pholy[1:5,]
print(filter.table.pholy[1:5,])
## print(tapply(filter.table.pholy[,2],filter.table.pholy[,2],length))
## print(tapply(filter.table.pholy[,4],filter.table.pholy[,4],length))
#put data together

## ######################PholyPhen Check
## interesting.coding.mutations.pholy<-c("frameshift substitution","nonframeshift substitution","nonframeshift deletion","nonframeshift insertion","nonsynonymous SNV","stopgain SNV","stoploss SNV","splicing")
## interesting.strict.coding.mutations<-interesting.coding.mutations.pholy[!(interesting.coding.mutations.pholy %in% c("splicing"))]
## wanted.muts.STRICT.coding<-test.for.coding.type(geneanno.table,geneanno.DB,interesting.strict.coding.mutations)

## poly.scores.cols<-paste(polyphen.types,c("pph2_prob"),sep=".")
## missing.poly<-apply(filter.table.pholy[,poly.scores.cols],1,function(x,num.cols){sum(is.na(x))==num.cols},length(poly.scores.cols)) # apply extra arguemnets lists after
## print(paste("Missing PholyPhen for :",sum( missing.poly & wanted.muts.STRICT.coding)," of :",sum(wanted.muts.STRICT.coding)))
## #######################################

  ###################################




  pholy.files<-dir(polyphen.dir) # location where would be if precalculated already

# if( grepl("_predictions.txt$",pholy.files)){

## if(project=="TGCM-AML"){
##   the.current.chr<-gsub("txt$","",target)
##   the.current.chr<-gsub("TGCM_AML.output.","",the.current.chr)
##     vep.fileS<-pholy.files[(grepl(the.current.chr,pholy.files) & grepl("redictions.txt$",pholy.files))]
##     setwd(polyphen.dir)
##  }else{ 

 if ( sum(  ( grepl(gsub("txt$","",target),pholy.files,fixed=TRUE) & grepl("redictions.txt$",pholy.files) ) | ( grepl(gsub(".txt$","_",target),pholy.files,fixed=TRUE) & grepl("redictions.txt$",pholy.files)  )    ) >0 ) { ## carefil chr1 vs chr11
    vep.fileS<-pholy.files[ ( grepl(gsub("txt$","",target),pholy.files,fixed=TRUE) & grepl("redictions.txt$",pholy.files) ) | ( grepl(gsub(".txt$","_",target),pholy.files,fixed=TRUE) & grepl("redictions.txt$",pholy.files)  )] ## target.txt -> target. fixed so get "." chr11 not confused with chr1 ALSO use _predictions instead of .predisctions so do 1_ and 10_ as well
    setwd(polyphen.dir)
  }else{

  pholy.format<-to.pholy.format(indels)
  dim(pholy.format)
 
  pholy.output.file<-paste(gsub(".txt$","",target),"_to_pholyPhen.txt",sep="")
  vep.file<-paste(gsub(".txt$","",target),"_predictions.txt",sep="")
  
  if( !(vep.file %in% pholy.files) ){

  xx<-try(setwd( polyphen.dir ),silent=TRUE)
  if(inherits(xx, "try-error")){system(paste("mkdir",polyphen.dir,sep=" "));setwd( polyphen.dir);QC.files<-dir(getwd())}else{QC.files<-dir(getwd())}

  
  setwd(polyphen.dir)
  write.table(pholy.format,file=pholy.output.file,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=FALSE)

## perl /media/scratch2/PholyPhen/variant_effect_predictor/variant_effect_predictor.pl -i test.txt -sift b -polyphen b -numbers -regulatory -domains -protein  -cache -dir  /media/scratch2/PholyPhen/.vep -o test.prodictions.txt
### Needed to copy variant_effect_predictor.pl to /usr/local/bin and use chmod 777 for this to work. Useing "perl" first caused problems. 
  system(paste("perl ",VEP.excutable.path,"/","variant_effect_predictor.pl -i ",pholy.output.file," -sift b -polyphen b -numbers -regulatory -domains -protein -force_overwrite -offline -cache -dir ",  VEP.excutable.path ," -o ",vep.file,sep="" ))
#VEP.excutable.path<-"/dmf/uqdi/Core_Services/UQCCG/Software/VEP" # variant_effect_predictor.pl
   pholy.files<-dir(polyphen.dir)
   vep.fileS<-pholy.files[(grepl(gsub(".txt$","",target),pholy.files) & grepl("redictions.txt$",pholy.files))]
#  cbind(filter.table.pholy[missing.poly & wanted.muts.STRICT.coding,],geneanno.table[missing.poly & wanted.muts.STRICT.coding,])[1:10,] -fork 4
}
} # esle catch is was already run

#} #project=="TGCM-AML"
print("d")
################## read in vep files  ivep<-1
for(ivep in 1:length(vep.fileS)){
  vep.file<-vep.fileS[ivep]

the.head<-try(scan(vep.file,what=character(),n=40,sep="\n",skip=0,fill=TRUE)) ## find the start of the vcf file
to.skip.lines<-grep("^#Uploaded_variation",the.head)

if(length(to.skip.lines)==0){print("ERROR pholyphen file missing header- NO predictions provided");next} # not assigned so makewill then trad all in next line 

options(show.error.messages = TRUE)
column.labels<-read.delim(vep.file,header=F,nrows=1,skip=(to.skip.lines-1),sep="\t",fill=TRUE,stringsAsFactors=FALSE)
num.vars<-dim(column.labels)[2]
column.labels<-gsub("^#","",column.labels)

### sometimes #Uploaded_variation" dones not preceed the data and there is an extra line
to.skip.lines<-match(TRUE,!grepl("^#",the.head))-1 # first line not starting with "#" , go one back since need the "skip" number
if(is.na(to.skip.lines)){to.skip.lines<-length(the.head)}
  
con.p <- file(vep.file, open="r")  # close(con)
num.lines<-1 # so does while llop at least once
reads.p<-0  #1.3M lines in snp file 50000 goes out to 24Gb without QC cgeck 
counter.p<- -1
#while (num.lines >0  ){
counter.p<-counter.p+1
print(counter.p)

if(counter.p==0){
a.poly<-try(scan(con.p,what=character(num.vars),skip=(reads.p*counter.p)+to.skip.lines,nlines=reads.p,sep="\t",fill=TRUE,na.strings=""))
}else{
a.poly<-try(scan(con.p,what=character(num.vars),nlines=reads,sep="\t",fill=TRUE,na.strings=""))
}
close(con.p)

num.lines<-length(a.poly)/(num.vars)
if(num.lines!=0){ 
print(num.lines)
#if(num.lines==0){next}
dim(a.poly)<-c(num.vars,num.lines)
a.poly<-t(a.poly)
colnames(a.poly)<-column.labels

if(ivep==1){
  poly<-a.poly
}else{
  poly<-rbind(poly,a.poly)
   }
 } # 
}  # loop over VEP files
  
print("e")
## unique(poly[,"Consequence"])
## dim(poly)

###################  GET  missing indels annotations from my unpacking methods
key.indels.poly<-build.key(indels,c("chr","start","ALT"))
print("key")
print(dim(poly))
print(poly[1:10,])
key.pholy<-build.pholy.key(poly)
print("pholy")
posns<-match(key.indels.poly,key.pholy)
print("posns")
missing<-is.na(posns)
print("missing")
print("f")
## #####################################################if pholyphen fucked use this to skip
## pholy.format<-to.pholy.format(indels[missing,]) 
## pholy.output.file<-paste(gsub(".txt$","",target),"_missing_to_pholyPhen.txt",sep="")
## vep.file<-paste(gsub(".txt$","",target),"_missing_predictions.txt",sep="")
## setwd(polyphen.dir)
## write.table(pholy.format,file=pholy.output.file,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=FALSE)
## missing<-FALSE
## #######################################################
print("in pholy...")
if(sum(missing)>0){
  print("moo pholy")
pholy.format<-to.pholy.format(subset(indels,missing)) # use sunset cause if pass only one element will throw an error (matrix-> vector)
pholy.output.file<-paste(gsub(".txt$","",target),"_missing_to_pholyPhen.txt",sep="")
vep.file<-paste(gsub(".txt$","",target),"_missing_predictions.txt",sep="")
setwd(polyphen.dir)
write.table(pholy.format,file=pholy.output.file,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=FALSE)


## perl /media/scratch2/PholyPhen/variant_effect_predictor/variant_effect_predictor.pl -i test.txt -sift b -polyphen b -numbers -regulatory -domains -protein  -cache -dir  /media/scratch2/PholyPhen/.vep -o test.prodictions.txt
### Needed to copy variant_effect_predictor.pl to /usr/local/bin and use chmod 777 for this to work. Useing "perl" first caused problems.

  system(paste("perl /mnt/UQCCG/Software/VEP/variant_effect_predictor.pl -i ",pholy.output.file," -sift b -polyphen b -numbers -regulatory -domains -protein -force_overwrite -offline -cache -dir  /mnt/UQCCG/Software/VEP/ -o ",vep.file,sep="" ))
 
  
the.head<-try(scan(vep.file,what=character(),n=1000,sep="\n",skip=0,fill=TRUE)) ## find the start of the vcf file
to.skip.lines<-grep("^#Uploaded_variation",the.head)
if(length(to.skip.lines)==0){print("ERROR pholyphen file missing header- NO predictions provided -MISSING INDELS");}else{ # not assigned so makewill then trad all in next line
  
options(show.error.messages = TRUE)
column.labels<-read.delim(vep.file,header=F,nrows=1,skip=(to.skip.lines-1),sep="\t",fill=TRUE,stringsAsFactors=FALSE)
num.vars<-length(column.labels)
column.labels<-gsub("^#","",column.labels)
### sometimes #Uploaded_variation" dones not preceed the data and there is an extra line
to.skip.lines<-grep(TRUE,!grepl("^#",the.head))-1 # first line not starting with "#" , go one back since need the "skip" number
if(length(to.skip.lines)==0){to.skip.lines<-length(the.head)}

con.p <- file(vep.file, open="r")  # close(con)
num.lines<-1 # so does while llop at least once
reads.p<-0  #1.3M lines in snp file 50000 goes out to 24Gb without QC cgeck 
counter.p<- -1
#while (num.lines >0  ){
counter.p<-counter.p+1
print(counter.p)
if(counter.p==0){
a.poly<-try(scan(con.p,what=character(num.vars),skip=(reads.p*counter.p)+to.skip.lines,nlines=reads.p,sep="\t",fill=TRUE,na.strings=""))
}else{
a.poly<-try(scan(con.p,what=character(num.vars),nlines=reads,sep="\t",fill=TRUE,na.strings=""))
}
close(con.p)
num.lines<-length(a.poly)/(num.vars)
if(num.lines!=0){
print(num.lines)
#if(num.lines==0){next}
dim(a.poly)<-c(num.vars,num.lines)
a.poly<-t(a.poly)
colnames(a.poly)<-column.labels
if(!exists("poly")){poly<-{}} # got no data

 if(is.null(dim(poly))){
 poly<-a.poly
 }else{
 poly<-rbind(poly,a.poly)
  }
 } #lines to process
} # skip cause pholyphen has no header


} # missing indels in VEP - mainly due to my mucking about in unwinding indels
####################################  Finished update up all VEP I can get!!!

## has.sift<-grepl("SIFT=",poly[,"Extra"])
## sum(has.pholy)
## sum(has.sift)
## poly[has.pholy,][1:10,]

########################use code below for new vep types:
## types<-sort(tapply(poly[,"Consequence"],poly[,"Consequence"],length))
## write.table(cbind(names(types),types),file="/media/Bioinform-D/Research/AML sequencing/function.types.txt",quote=TRUE,row.names=FALSE,sep="\t")

## vep.types<-read.delim("/media/Bioinform-D/Research/AML sequencing/function.types.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
## vep.names<-vep.types[,"Variant"]
## vep.types<-vep.names[!duplicated(vep.names)]
## toString(vep.names)
## ## the use replace to fix

## vep.types<-c("not_assigned","stop_gained","stop_lost","missense_variant","splice_acceptor_variant","splice_donor_variant","splice_region_variant","initiator_codon_variant","stop_retained_variant","incomplete_terminal_codon_variant","frameshift_variant","inframe_deletion","inframe_insertion","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_exon_variant","NC_stop_gained","NC_stop_lost","NC_splice_acceptor_variant","NC_splice_donor_variant","NC_splice_region_variant","NC_initiator_codon_variant","NC_stop_retained_variant","NC_non_coding_exon_variant","NC_incomplete_terminal_codon_variant","NC_3_prime_UTR_variant","mature_miRNA_variant","NC_5_prime_UTR_variant","TF_binding_site_variant","TFBS_ablation","TFBS_amplification","regulatory_region_variant","intron_variant","NC_intron_variant","synonymous_variant","coding_sequence_variant","NC_synonymous_variant","upstream_gene_variant","downstream_gene_variant","intergenic_variant","NC_intergenic_variant","NMD_transcript_variant","nc_transcript_variant","NC_nc_transcript_variant","feature_truncation","feature_elongation")

## Consequence<-get.Consequence.from.VEP(indels,vep.types,poly,"Consequence",c("Uploaded_variation","Gene","Feature"))

print("skipped")

 if(!is.null(dim(poly))){ ## no function data so just jump out  c("Uploaded_variation","Gene","Feature","Protein_position","Amino_acids")

   
Consequence<-get.Consequence.from.VEP(indels,vep.types,poly,"Consequence",extra.vep.annotations) ## vep types is ordered list of all possible mutation worst at front

pholy.data<-get.function.from.VEP(indels,poly,"PolyPhen=")
sift.data<-get.function.from.VEP(indels,poly,"SIFT=")

## tapply(poly[,"Feature_type"],poly[,"Feature_type"],length)
regulation.data<-get.regulation.from.VEP(indels,poly,c("MotifFeature","RegulatoryFeature"))

sift.data[,"SIFT.scores"]<-1.0-as.numeric(sift.data[,"SIFT.scores"])
pholy.calibration<-functional.score.to.desc.calibrate(pholy.data,"PolyPhen.scores","PolyPhen.desc")
sift.calibration<-functional.score.to.desc.calibrate(sift.data,"SIFT.scores","SIFT.desc")
##################################################################### finish get pholyPhen


############################## Combine functional data Update where necessaey - not needed after fis pholyphen
############################## WARNING ONLY SCORES ARE UPDATED. DESCRIPTions ARE NOT
#filter.table.pholy<-filter.table.pholy.ori

filter.table.pholy<-cbind(filter.table.pholy,pholy.data,sift.data,regulation.data,Consequence)
#filter.table.pholy[110:120,]

###################### Merge OLD/NEW PholyPhen
  poly.scores.cols<-paste(polyphen.types,c("pph2_prob"),sep=".")
  missing.poly<-apply(filter.table.pholy[,poly.scores.cols],1,function(x,num.cols){sum(is.na(x))==num.cols},length(poly.scores.cols))

  new.missing.poly<-is.na(filter.table.pholy[,"PolyPhen.scores"])
  transfer<-!missing.poly & new.missing.poly
  # filter.table.pholy[transfer,]
# sum(transfer)
# sum(missing.poly)
# sum(new.missing.poly)
# sum(!missing.poly & new.missing.poly) # present in old pholyphen and missing in new
  filter.table.pholy[transfer,"PolyPhen.scores"]<-filter.table.pholy[transfer,"GWAS.pph2_prob"]
  filter.table.pholy[transfer,"PolyPhen.desc"]<-filter.table.pholy[transfer,"GWAS.prediction"]
  filter.table.pholy[transfer,"PolyPhen.desc"]<-gsub("possibly damaging","possibly_damaging",filter.table.pholy[transfer,"PolyPhen.desc"])
  filter.table.pholy[transfer,"PolyPhen.desc"]<-gsub("probably damaging","probably_damaging",filter.table.pholy[transfer,"PolyPhen.desc"])
###################


########################### MERGE ANNOVAR DATA 
# filter.table[,c("GERP::score","ljb_gerp::score","pholyphen::score","sift::score","mut.taster::score","phylo::score")]
# filter.table.pholy[1:5,]
###################### Merge SIFT  
# names(function.filter.DB)
#### filter.table[,c("GERP::score","ljb_gerp::score","pholyphen::score","sift::score","mut.taster::score","phylo::score")]
SKIP.merge<-FALSE
  if("sift" %in% names(function.filter.DB)){ filter.sift<-function.filter.DB[names(function.filter.DB)=="sift"];SKIP.merge<-FALSE}else{SKIP.merge<-TRUE} #

if(!SKIP.merge){
  missing.poly<-is.na(filter.table[,filter.sift]) # from ANNOVAR and in filter table
  new.missing.poly<-is.na(filter.table.pholy[,"SIFT.scores"])
  transfer<-!missing.poly & new.missing.poly
  # sum(transfer)
  # cbind(filter.table.pholy[transfer,],filter.table[transfer,"sift::score"])
  filter.table.pholy[transfer,"SIFT.scores"]<-filter.table[transfer,filter.sift]
#if(sum(transfer)!=0 & dim(sift.calibration)[1]!=0) { 
  if(sum(transfer)!=0 & !is.null(dim(sift.calibration))) {  #### Changed by Mhairi June 2014
  filter.table.pholy[transfer,"SIFT.desc"]<-apply.calibration(filter.table[transfer,filter.sift],sift.calibration)
                      }
} # skip merge
###################### Merge PholyPhen
#### filter.table[,c("GERP::score","ljb_gerp::score","pholyphen::score","sift::score","mut.taster::score","phylo::score")]
#  if("pholyphen::score" %in% colnames(filter.table)){filter.pholy<-"pholyphen::score" }else{filter.pholy<-"ljb_pp2::score"}
SKIP.merge<-FALSE
  if("pholyphen" %in% names(function.filter.DB)){ filter.pholy<-function.filter.DB[names(function.filter.DB)=="pholyphen"];SKIP.merge<-FALSE}else{SKIP.merge<-TRUE} #
if(!SKIP.merge){
  missing.poly<-is.na(filter.table[,filter.pholy]) # from ANNOVAR and in filter table
  new.missing.poly<-is.na(filter.table.pholy[,"PolyPhen.scores"])
  transfer<-!missing.poly &   new.missing.poly
  # sum(transfer)
  # cbind(filter.table.pholy[transfer,],filter.table[transfer,"pholyphen::score"])
  filter.table.pholy[transfer,"PolyPhen.scores"]<-filter.table[transfer,filter.pholy]
 ### if(sum(transfer)!=0 & dim(sift.calibration)[1]!=0){ 
    if(sum(transfer)!=0 & !is.null(dim(pholy.calibration))){   #### Changed by Mhairi June 2014
  filter.table.pholy[transfer,"PolyPhen.desc"]<-apply.calibration(filter.table[transfer,filter.pholy],pholy.calibration)
                       }
  } # skip merge

} #no functional data so jump out
##############################
 print("done")
 setwd(annotate.dir) ### casue later assume this this is the working dir

## new.missing<-is.na(test[,"PolyPhen.scores"])
## sum(!missing.poly & new.missing)
## a.test<-!missing.poly & new.missing
## sum(a.test) # 92 for AOGC chr 10 ### these appear to be transcript misses in Ensemble 
## ## a.test<-missing.poly & !new.missing # so ne method working better but missing soem
## ## sum(a.test) # 965 for AOGC chr 10 
##########################################################################################################
