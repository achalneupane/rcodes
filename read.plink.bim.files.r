
#chrom	SNP	position	ALLELES	FREQ1	RSQR	EFFECT1	STDERR	WALDCHISQ	-log10(Pval)	origin
#SNP POS Al1 Al2 CHR   	MARKER          ALLELES  FREQ1    RSQR   EFFECT1  OR      STDERR  WALDCHISQ PVALUE     LRCHISQ LRPVAL


# i<-1
for(i in 1:length(sample.files)){
  
################ LARGE FILES ########
print(sample.files[i])
## setwd(eval(as.name(paste(names(sample.files)[i],"dir",sep=".")))) ## different directory for each indel type
setwd(snp.dir)
######BELOW PROCESSING  this for snp for indel varient types in vcf4.0 or vcf 3.0 format

  
### get location of header:
chromo1<-try(scan(sample.files[i],what=character(),n=200,sep="\n",skip=0,fill=TRUE)) ## find the start of the vcf file
skip.lines<-grep("CHR",chromo1)
if(length(skip.lines)>1){print("ERROR multiple chrom lables found");skip.lines<-skip.lines[1]}
if(length(skip.lines)==0){print("NO column labels assume is a bim file");skip.lines<-0}



if(skip.lines==0){column.labels<-c("chr","SNP","junk","POS","REF","ALT");num.vars<-length(column.labels)}else{
options(show.error.messages = TRUE)
column.labels<-read.table(sample.files[i],header=F,nrows=1,skip=(skip.lines-1),fill=TRUE,stringsAsFactors=FALSE)
num.vars<-dim(column.labels)[2]
}

indels<-try(scan(sample.files[i],what=character(num.vars),skip=skip.lines,fill=TRUE))
num.lines<-length(indels)/(num.vars)
dim(indels)<-c(num.vars,num.lines)
indels<-t(indels)

############### ADD IF HERE IF HAVE DIFFERENT ASSOC FORMATS column.labels<-c("chr","start","end","REF","ALT","maf","geno")


column.labels[column.labels=="chrom"]<-"chr"
column.labels[column.labels=="CHR"]<-"chr"
column.labels[column.labels=="position"]<-"POS"
column.labels[column.labels=="BP"]<-"POS"
column.labels[column.labels=="BP"]<-"POS"
column.labels[column.labels=="Al1"]<-"A1"
column.labels[column.labels=="Al2"]<-"A2"
colnames(indels)<-column.labels 


####################################### FINISHED Read in data

######## in plink assoc files need to redo  to get alleles correct
if(sum(c("REF","ALT") %in% column.labels)!=2){ ## see if have REF and ALT
if(sum(c("A1","A2") %in% column.labels)!=2){
 has.comma<-grepl(",",indels[,"ALLELES"])
 if(sum(has.comma)!=dim(indels)[1]){ print("ERROR alleles not formated as expected");next}
 alleles<-strsplit(indels[,"ALLELES"],split=",")
 REF<-unlist(lapply(alleles,function(x) x[1]))
 ALT<-unlist(lapply(alleles,function(x) x[2]))           
                                   }else{
                                     REF<-indels[,"A1"]
                                     ALT<-indels[,"A2"]
}

indels<-cbind(indels,REF=REF,ALT=ALT)
}
#################### FINISHED processing individual filters
print(indels[1:5,])
############################### MAKE the final data file

########################### PUT ALLELES IN ANNOVAR FORMAT from convert2annovar.pl line 1083########
ref.length<-nchar(as.character(indels[,"REF"]))
alt.length<-nchar(as.character(indels[,"ALT"]))
is.snp<-(ref.length==1 & alt.length==1)

# sum(!is.snp)

POS.end<-as.numeric(indels[,"POS"])
del<-ref.length > alt.length
ins<-(ref.length <= alt.length) & !is.snp
## indels[del,][1:5,]
## POS.end[del][1:5]
### deletion or block substitution
head<-substr(as.character(indels[del,"REF"]),1,alt.length[del])
head.is.mut<-(head==as.character(indels[del,"ALT"]))
indels[del,"REF"][head.is.mut]<-substr(as.character(indels[del,"REF"][head.is.mut]),(alt.length[del][head.is.mut]+1),ref.length[del][head.is.mut])
indels[del,"ALT"][head.is.mut]<-"-"
indels[del,"POS"][head.is.mut]<-as.numeric(indels[del,"POS"][head.is.mut]) + nchar(as.character(head[head.is.mut]))
POS.end[del]<-POS.end[del]+ref.length[del]-1  # same for both head is mut and not head is mut

## indels
## POS.end
### insertion or block substitution
head<-substr(as.character(indels[ins,"ALT"]),1,ref.length[ins])
head.is.ref<-(head==as.character(indels[ins,"REF"]))
indels[ins,"ALT"][head.is.ref]<-substr(as.character(indels[ins,"ALT"][head.is.ref]),(ref.length[ins][head.is.ref]+1),alt.length[ins][head.is.ref])
indels[ins,"REF"][head.is.ref]<-"-"
indels[ins,"POS"][head.is.ref]<-as.numeric(indels[ins,"POS"][head.is.ref]) + ref.length[ins][head.is.ref]-1
POS.end[ins]<-POS.end[ins]+ref.length[ins]-1
########################################################################################################
#$$$  all.alleles  data<-cbind(indels,"SNP","+")  ; colnames(data)<-c(colnames(indels),"TYPE","strand")

info<-indels[,!(colnames(indels) %in% c("chr","POS"))]

indels<-cbind(indels[,"chr"],indels[,"POS"],POS.end,"+",names(sample.files)[i],info)

colnames(indels)<-c("chr","start","end","strand","TYPE",colnames(info))
indels[1:5,]
gr.core.ann<-c( "chr","start","end","strand") #'elementMetadata' cannot use "seqnames", "ranges", "strand", "seqlengths", "start", "end", "width", or "element" as column names


##########################
setwd(code.dir)
source("get.forward.strand.allele.r")

setwd(snp.dir)
assign(sample.grs[i],value=indels)

## data.gr<-GRanges(seqnames =indels[,"chr"],ranges = IRanges(start=as.numeric(indels[,"start"]),end=as.numeric(indels[,"end"])),strand=indels[,"strand"])
## values(data.gr)<-indels[,!(colnames(data) %in% gr.core.ann)]
## assign(sample.grs[i],value=data.gr) sample.grs[i]<-"snp"
## data.gr
print(paste("Done",sample.grs[i],sep=" "))
}
