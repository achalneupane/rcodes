
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
assoc.file<-"/media/scratch2/GBS/GBS.assoc.assoc.logistic" # chisq file
bim.file<-"/media/scratch2/GBS/GBS_wtccc_commom.bim"  ### bim file

bim.file.cols<-c("SNP","chr","AF","POS","A1","A2")
assoc.file.cols<-c("Chisq","Chisq-eigenstrat")
bim.has.header<-F
assoc.has.header<-T
#################################################






setwd(plot.dir)
bim<-read.table(bim.file,header=bim.has.header,fill=TRUE,stringsAsFactors=FALSE)
assoc<-read.table(assoc.file,header=assoc.has.header,fill=TRUE,stringsAsFactors=FALSE)

dim(assoc)
dim(bim)
assoc[1:2,]
bim[1:2,]

colnames(assoc)<-assoc.file.cols
colnames(bim)<-bim.file.cols
assoc[1:5,]
bim[1:5,]

p.value<-pchisq(as.numeric(assoc[,"Chisq"]),1,lower.tail=FALSE)
p.value.eigenstrat <- pchisq(as.numeric(assoc[,"Chisq-eigenstrat"]),1,lower.tail=FALSE)

data<-cbind(bim,assoc,p.value,p.value.eigenstrat)

write.table(data,file="final.eigen.assoc",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


############### STOP give "final.eigen.assoc" to katie


##################################
core.ann<-c("chr","start","end","REF","ALT","TYPE")
code.dir<-"/media/Bioinform-D/Research/AML sequencing"
genome.build<-"hg19"
target<-"temp"

data[1:5,]
indels<-cbind(data[,c("chr","POS","POS","A1","A2")],"snp")
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

setwd(
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
sum(!missing)
assoc.ann<-cbind(assoc,ann[posns,])
remove.cols<-c("Chr")
assoc.ann<-assoc.ann[,!(colnames(assoc.ann) %in% remove.cols)]


assoc.ann[1:5,]
ann[1:5,]
setwd(plot.dir)
assoc.ann[,"AAChange.refGene"]<-gsub(":",".",assoc.ann[,"AAChange.refGene"])
assoc.ann[,"ExonicFunc.refGene"]<-gsub(" ",".",assoc.ann[,"ExonicFunc.refGene"]) # does not like spaces

write.table(assoc.ann,file="final_annotated.assoc",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

write.table(assoc.ann[1:10,],file="final_annotated2.assoc",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
[,c(1:16)]











## filter.DB<-c("1000g2010nov_all","snp132") ## here dbsnp132  is just annoting not doing filtering  ## 1000g2010nov_all==hg19_ALL.sites.2010_11.txt
for(i in 1:length(filter.DB)){
if( (paste(target,filter.DB[i],"log",sep=".") %in% files.in.annotate ) & !update.annovar.annotations){next} # perhaps have added batabases but dones want to update
system(
    paste("annotate_variation.pl -batchsize 50m -filter --score_threshold ",maf.threshold," -buildver ",genome.build," -dbtype ",filter.DB[i]," ",target,"  ",anno.DB.location," --outfile ",paste(target,filter.DB[i],sep=".")," ",sep="" )
      )
}
############################################################


for(i in 1:length(geneanno.DB)){
if( (paste(target,geneanno.DB[i],"log",sep=".") %in% files.in.annotate ) & !update.annovar.annotations){next} # perhaps have added batabases but dones want to update
system(
    paste("annotate_variation.pl -batchsize 50m -geneanno --splicing_threshold ",splice.threshold," -buildver ",genome.build," -dbtype ",geneanno.DB[i]," ",target,"  ",anno.DB.location," --outfile ",paste(target,geneanno.DB[i],sep=".")," ",sep="" )
      )
}

if(grepl("^chrMT",target) | grepl(".chrMT.",target)){  ## mitochrondrial using 1000 genomes FASTA - which we use replaces the ensGene one written above for mitochromdria (only works if doing per chromosome

  if( (paste(target,geneanno.DB[i],"log",sep=".") %in% files.in.annotate ) & !update.annovar.annotations){next} # perhaps have added batabases but dones want to update
system(
    paste("annotate_variation.pl -batchsize 50m -geneanno --splicing_threshold ",splice.threshold," -buildver ",genome.build, "-dbtype MT_GRCh37_ensGene ",target,"  ",anno.DB.location," --outfile ",paste(target,"ensGene",sep=".")," ",sep="" )
      )
}


