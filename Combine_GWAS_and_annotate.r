
## shellfish.py --file <.geno and .map prefix> --numpcs 10 --evecs <.evecs
## file> --snpload --out <out file prefix>

## could try admixture as well
## http://www.genetics.ucla.edu/software/admixture/
###########################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-TRUE
plot.dir<-"/media/scratch2/Chinese_MND_GWAS/Katie"
assoc.file<-"/media/scratch2/Chinese_MND_GWAS/Katie/eigenstrat_res_chisq.txt" # chisq file
bim.file<-"/media/scratch2/Chinese_MND_GWAS/Katie/eigenstrat.snp.txt"  ### bim file

bim.file.cols<-c("SNP","chr","AF","POS","A1","A2")
assoc.file.cols<-c("Chisq","Chisq-eigenstrat")
bim.has.header<-F
assoc.has.header<-T
#################################################


###########################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-TRUE
plot.dir<-"/media/scratch2/GBS"
assoc.file1<-"/media/scratch2/GBS/assoc_res_bcformat.txt" # chisq file
assoc.file2<-"/media/scratch2/GBS/final_annotated.assoc" 

bim.file<-"/media/scratch2/GBS/GBS_hg19.bim"  ### bim file that is the superset of all genotypes

bim.file.cols<-c("chr","SNP","AF","POS","A1","A2")
assoc.file.cols1<-c("Chisq","Chisq-eigenstrat")
assoc.file.cols2<-c("Chisq","Chisq-eigenstrat")

bim.has.header<-FALSE
assoc.has.header1<-TRUE
assoc.has.header2<-TRUE
#################################################






setwd(plot.dir)
bim<-read.table(bim.file,header=bim.has.header,fill=TRUE,stringsAsFactors=FALSE)
assoc1<-read.table(assoc.file1,header=assoc.has.header1,fill=TRUE,stringsAsFactors=FALSE)
assoc2<-read.table(assoc.file2,header=assoc.has.header2,fill=TRUE,stringsAsFactors=FALSE)

if(!bim.has.header){
colnames(bim)<-bim.file.cols
}

if(!assoc.has.header1<){
colnames(assoc1)<-bim.file.cols1
}

if(!assoc.has.header2<){
colnames(assoc1)<-bim.file.cols2
}

dim(assoc1)
dim(assoc2)
dim(bim)
assoc1[1:2,]
assoc2[1:2,]
bim[1:2,]

## olnames(assoc)<-assoc.file.cols
## colnames(bim)<-bim.file.cols
## assoc[1:5,]
## bim[1:5,]


## ############################################## combine data anyway using largest data file
## p.value<-pchisq(as.numeric(assoc[,"Chisq"]),1,lower.tail=FALSE)
## p.value.eigenstrat <- pchisq(as.numeric(assoc[,"Chisq-eigenstrat"]),1,lower.tail=FALSE)
## data<-cbind(bim,assoc,p.value,p.value.eigenstrat)
## ##############################################

assoc1[,"SNP"]<-gsub("^RS","rs",assoc1[,"SNP"])
not.an.rs<- !grep("rs",assoc1[,"SNP"])
sum(not.an.rs) ## if not zero do something



######### checck bim is a superset
missing<-!(assoc1[,"SNP"] %in% bim[,"SNP"])
sum(missing) # if not zero are missing some
##### check to see if those missing are worth looking at;
assoc1[missing,][1:5,]
the.order<-order(as.numeric(assoc1[missing,"EIGENSTRAT_P"]))
assoc1[missing,][the.order,] ## all chr24 or not significant

### remove 
 assoc1<-assoc1[!missing,]           

missing<-!(assoc2[,"SNP"] %in% bim[,"SNP"])
sum(missing)

            
bim[1:5,]
assoc1[1:5,]
assoc2[1:5,]
dim(bim)
dim(assoc1)
dim(assoc2)

colnames(assoc1)<-paste(colnames(assoc1),"Matt",sep="_")
colnames(assoc2)<-paste(colnames(assoc2),"WTCCC",sep="_")


posns<-match(bim[,"SNP"],assoc1[,"SNP_Matt"])
missing<-is.na(posns)
sum(missing)


data<-cbind(bim,assoc1[posns,c("SNP_Matt","EIGENSTRAT_P_Matt")])

posns<-match(bim[,"SNP"],assoc2[,"SNP_WTCCC"])
missing<-is.na(posns)
sum(missing)


data<-cbind(data,assoc2[posns,c(2,4:12)])

both.missing<-is.na(data[,"SNP_Matt"]) & is.na(data[,"SNP_WTCCC"])
both.missing[1:10]

sum(both.missing)

data<-data[!both.missing,]

data[1:5,]         

write.table(data,file="final.combined.assoc",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


############### STOP give "final.eigen.assoc" to katie


##################################
core.ann<-c("chr","start","end","REF","ALT","TYPE")
code.dir<-"/media/Bioinform-D/Research/AML sequencing"
genome.build<-"hg19"
target<-"temp"

data[1:5,]
indels<-cbind(data[,c("chr","POS","POS","A1","A2")],"snp")
indels[1:5,]
colnames(indels)<-core.ann

annotate.GWAS<-function(indels,genome.build,geneanno.DB,filter.DB)

geneanno.DB<-c("refGene")  # returns 2 extra columns c("refGene","knownGene","ensGene")
names(geneanno.DB)<-c("refGene") # c("refGene","knownGene","ensGene")

filter.DB<-c("snp137")  # returns 2  extra columns (DB,score) c("snp132","1000g2011may_all")
names(filter.DB)<-c("ID")

setwd(code.dir)
source("ucsc.table.names.r")   # load in the UCSC tables these use the db file names and not their lable-names 
source("annotate_SNPs_subroutines.r")
source("ucsc.table.names.processor.r")
setwd(code.dir)
source("get.forward.strand.allele.r")

anno.DB.location.core<-"/media/Bioinform-D/Research/annovar/humandb"
anno.DB.location<-paste(anno.DB.location.core,genome.build,sep="/")

setwd(plot.dir)
write.table(indels[,core.ann],file=target,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

write.table(indels[,core.ann],file=target,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

#table_annovar.pl ex1.human humandb/ -protocol refGene,phastConsElements44way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp135,avsift,ljb_all -operation g,r,r,f,f,f,f,f -nastring NA

system(paste("table_annovar.pl ",target,"  ",anno.DB.location," -buildver ",genome.build," -protocol refGene,snp137 -operation g,f --outfile ",paste(target,"txt",sep=".")," ",sep="" )  )
#system(paste("table_annovar.pl ",target,"  ",anno.DB.location," -buildver ",genome.build," -protocol refGene,snp137 -operation g,f --outfile ",paste(target,"txt",sep=".")," ",sep="" )  )


### get location of header:
## input<-paste(target,"txt",paste(genome.build,"multianno",sep="_"),"txt",sep=".")
## column.labels<-read.table(input,header=F,nrows=1,skip=0,fill=TRUE,stringsAsFactors=FALSE)
## num.vars<-dim(column.labels)[2]
## skip.lines<-1


## ann<-try(scan(input,what=character(num.vars),skip=skip.lines,fill=TRUE))
## num.lines<-length(ann)/(num.vars)
## dim(ann)<-c(num.vars,num.lines)
## ann<-t(ann)
## colnames(ann)<-column.labels 
############### ADD IF HERE IF H


ann<-read.delim(paste(target,"txt",paste(genome.build,"multianno",sep="_"),"txt",sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
ann[ann==""]<-NA  ## Haploviews does not like ""
dim(ann)
dim(assoc)
dim(indels)
dim(ann)
ann.key<-build.key(ann,c("Chr","Start","End","Ref","Alt"))
indel.key<-build.key(indels,c("chr","start","end","REF","ALT"))

posns<-match(indel.key,ann.key)
missing<-is.na(posns)
sum(missing)
assoc.ann<-cbind(data,ann[posns,])
remove.cols<-c("Chr","Start","End","Ref","Alt","SNP_WTCCC","A1_WTCCC","TEST_WTCCC","NMISS_WTCCC")
assoc.ann<-assoc.ann[,!(colnames(assoc.ann) %in% remove.cols)]


assoc.ann[1:5,]
ann[1:5,]
setwd(plot.dir)
assoc.ann[,"AAChange.refGene"]<-gsub(":",".",assoc.ann[,"AAChange.refGene"])
assoc.ann[,"ExonicFunc.refGene"]<-gsub(" ",".",assoc.ann[,"ExonicFunc.refGene"]) # does not like spaces



write.table(assoc.ann,file="final_combined_annotated.assoc",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

write.table(assoc.ann[1:10,],file="final_annotated2.assoc",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
[,c(1:16)]




##############
assoc.ann<-assoc.ann[,c(1:8,14,9:13,15:19)]











