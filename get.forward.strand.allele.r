


library("Biostrings")
library("IRanges")

if(!exists("genome.build")){genome.build<-"hg19"}
if(genome.build=="hg19"){library("BSgenome.Hsapiens.UCSC.hg19"); have.chr<-seqnames(Hsapiens);chr.lengths<-seqlengths(Hsapiens)}
if(genome.build=="hg18"){library("BSgenome.Hsapiens.UCSC.hg18"); have.chr<-seqnames(Hsapiens);chr.lengths<-seqlengths(Hsapiens)}
if(genome.build=="mm9"){library("BSgenome.Mmusculus.UCSC.mm9"); have.chr<-seqnames(Mmusculus)}
if(genome.build=="mm10"){library("BSgenome.Mmusculus.UCSC.mm10"); have.chr<-seqnames(Mmusculus)}
######## CAREFUL this code used by read.plink.assoc.files.r
print(paste("Genome build is:",genome.build))
print("inside forward strand")
 print(dim(indels))
print("inside forward strand")

if(!grepl("^chr",indels[1,"chr"])){indels[,"chr"]<-paste("chr",indels[,"chr"],sep="")}

if(!("SNP" %in% colnames(indels))){
  print("WARNING columns SNP not in indels change col  ID->SNP")
    if("ID" %in% colnames(indels)){ colnames(indels)[colnames(indels)=="ID"]<-"SNP"}else{
      print("FAIL columns SNP or ID  not in indels USING KEY")
      SNP<-build.key(indels,core.ann)
      indels<-cbind(indels,SNP)
      }
}


if(!exists("sample.files")){sample.files<-c("data");i<-1}
if(!exists("snp.dir")){snp.dir<-getwd()}

## if(!grepl("^chr",indels[1,"chr"])){
## key.indels<-build.key(all,core.ann,add.chr.label=TRUE)
##   indels[,"chr"]<-paste("chr",indels[,"chr"],sep="")
## }else{key.indels<-build.key(indels,core.ann)      }
## indels[1:5,]


## names(types)[grep("frameshift_variant",names(types))]
## types<-sort(tapply(poly[,"Consequence"],poly[,"Consequence"],length))
## to.get<-"Consequence"
## order.by<-vep.types
##  to.get.other<-c("Uploaded_variation","Gene","Feature","Protein_position","Amino_acids")
## to.get.other<-{}




indels<-as.matrix(indels)  ## else get factor issues

#no.chr<-check.chr.and.positions(indels,chr.lengths)

indels[,"chr"]<-gsub("^\\s+","",indels[,"chr"])
indels[,"chr"]<-gsub("\\s+ $","",indels[,"chr"])

indels[,"start"]<-gsub("^\\s+","",indels[,"start"])
indels[,"start"]<-gsub("\\s+ $","",indels[,"start"])

indels[,"end"]<-gsub("^\\s+","",indels[,"end"])
indels[,"end"]<-gsub("\\s+ $","",indels[,"end"])

 positions<-chr.lengths[indels[,"chr"]]

#as.integer(indels[,"start"])

indels[indels[,"chr"]=="chr23","chr"]<-"chrX"  ## no Y observed
indels[indels[,"chr"]=="chr24","chr"]<-"chrY"
indels[indels[,"chr"]=="chr25","chr"]<-"chrXY"
indels[indels[,"chr"]=="chr26","chr"]<-"chrM"

## tapply(indels[,"chr"],indels[,"chr"],length)

no.chr<-check.chr.and.positions(indels,chr.lengths)
sum(no.chr)

if(sum(no.chr)>1){
print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
exclude.file<-paste("to.exclude",sample.files[i],"txt",sep=".")
print(paste("Missing chromosomes/off end of chromosomes  so won't be tested write to ",exclude.file,sep=""))
print(indels[no.chr,])
setwd(snp.dir)
if("SNP" %in% colnames(indels)){
write.table(indels[no.chr,"SNP"],file=exclude.file,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
}else{
write.table(indels[no.chr,"ID"],file=exclude.file,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
}
setwd(code.dir)
print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
}
## the.chrom[1:5]
## tapply(the.chrom,the.chrom,length)


    ##  X    X chromosome                    -> 23
    ##  Y    Y chromosome                    -> 24
    ##  XY   Pseudo-autosomal region of X    -> 25
    ##  MT   Mitochondrial                   -> 26
# Hsapiens


the.chrom<-as.character(indels[!no.chr,"chr"])
starts<-as.integer(as.character(indels[!no.chr,"start"]))
ends<-as.integer(as.character(indels[!no.chr,"end"]))
all.genomic.have<-getSeq(Hsapiens, the.chrom, starts, ends,as.character=TRUE) ## get forward stand alleles


all.genomic<-indels[,"REF"]
all.genomic[!no.chr]<-all.genomic.have #
all.genomic.ori<-all.genomic
## all.genomic[1:5]
## indels[,"REF.base"]<-
## alleles<-c("A","C","G","T")

ref.ori<-indels[,"REF"]
alt.ori<-indels[,"ALT"]

to.flip<-rep(0,times=length(ref.ori))
##TRY 
## (1) SWAP (+1 to to.flip)
## (2) Complement
## (3) SWAP  (-1 to to.flip)
## (4) Complent back so  no-fix's are back to original
## (5) assume what is left if forward strand alleles that are insertions
## to.flip !=0 need flipping
########### (1) ##########################################################################################
ref.diff<-indels[,"REF"]!=all.genomic & indels[,"REF"]!="-"  ## if a insertion can't tell if it match the genome
sum(ref.diff)

if(sum(ref.diff)!=0){
to.flip[ref.diff]<-to.flip[ref.diff]+1
## sum(ref.diff) #4601
print(paste(sum(ref.diff)," STEP1: strand differences or minor allele not REF"))
## cbind(indels,all.genomic)[ref.diff,][1:10,]

### swap REF AND ALT to see is just minor allele / major allele swap
the.ref<-indels[ref.diff,"REF"]
indels[ref.diff,"REF"]<-indels[ref.diff,"ALT"]
indels[ref.diff,"ALT"]<-the.ref
###### swap 0/1 --- keep track of times swapped

############ for indels must chenge start and end and the get new genomic sequence
an.indel<- (nchar(as.character(indels[ref.diff,"REF"])) >1) | (nchar(as.character(indels[ref.diff,"ALT"])) > 1) | indels[ref.diff,"REF"]=="-" | indels[ref.diff,"ALT"]=="-"
sum(an.indel)

if(sum(an.indel)>0){
indels[ref.diff,][an.indel,c("chr","start","end","REF","ALT")]<-redo.start.end.annovar(indels[ref.diff,][an.indel,c("chr","start","end","REF","ALT")])


no.chr<-check.chr.and.positions(indels[ref.diff,][an.indel,c("chr","start","end","REF","ALT")],chr.lengths)

the.chrom<-as.character(indels[ref.diff,"chr"][an.indel][!no.chr])
starts<-as.integer(as.character(indels[ref.diff,"start"][an.indel][!no.chr]))
ends<-as.integer(as.character(indels[ref.diff,"end"][an.indel][!no.chr]))
all.genomic.indels<-getSeq(Hsapiens, the.chrom, starts, ends,as.character=TRUE)
all.genomic[ref.diff][an.indel][!no.chr]<-all.genomic.indels
}
#cbind(all.genomic[!no.chr][ref.diff][an.indel],all.genomic.indels)
#### TOP /BOT diff
####################### issue here is that for INSERTIONS turned to DELETIONS may my chance match the reference mase cause insert/deletion is at a homopolymer run!
####################### so ca't tell IF it was more frequent or A1/ A2 or insetion or deletion


} # something to do in step 1

########### (2) ##########################################################################################

## ref.diff<-indels[,"REF"]!=all.genomic
## an.indel<- (nchar(as.character(indels[ref.diff,"REF"])) >1) | (nchar(as.character(indels[ref.diff,"ALT"])) > 1) | indels[ref.diff,"REF"]=="-" | indels[ref.diff,"ALT"]=="-"
## sum(ref.diff[an.indel]) # 6
## forward<-cbind(indels,all.genomic)[ref.diff,]
## save(list=c("forward"),file="forward.RData")
## cbind(forward[an.indel,],all.genomic.ori[ref.diff][an.indel],indels.ori[ref.diff,][an.indel,])

## cbind(forward,all.genomic.ori[ref.diff],indels[ref.diff,])[1:10,]
########################################################################################################
## IUPAC_CODE_MAP
   ##   A      C      G      T      M      R      W      S      Y      K      V      H      D      B      N 
   ## "A"    "C"    "G"    "T"   "AC"   "AG"   "AT"   "CG"   "CT"   "GT"  "ACG"  "ACT"  "AGT"  "CGT" "ACGT"

ref.diff<-indels[,"REF"]!=all.genomic
sum(ref.diff)


print(paste(sum(ref.diff),":STEP 2 Strand differences use complment DANGER with polymorphic alleles ASSUMING ONE IS THE REF ALLELE to flip written"))
print(paste("Flip file called:",paste("to.flip",sample.files[i],"txt",sep=".")," put in :",snp.dir))
 setwd(snp.dir)      
write.table(indels[ref.diff,"SNP"],file=paste("to.flip",sample.files[i],"txt",sep="."),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
setwd(code.dir)
## load("/media/scratch2/AOGC-NGS/ExomeChip/polymorphic.snps.RData") # save(list=c("poly.morphic.snps","all.data"),file="polymorphic.snps.RData")
## dim(all.data)
## cbind(forward,all.genomic.ori[ref.diff],indels[ref.diff,])[1:10,]
## posns<-match(indels[,"SNP"],poly.morphic.snps)
## poly.morphic.snps<-poly.morphic.snps[posns]
## cbind(forward,all.genomic.ori[ref.diff],indels[ref.diff,],poly.morphic.snps[ref.diff])[1:10,]
## test<- ref.diff & !is.na(poly.morphic.snps)
## sum(test)
## posns<-match(indels[,"SNP"],all.data[,"SNP"])
## all.data<-all.data[posns,]
## cbind(indels,all.genomic.ori,poly.morphic.snps,all.data)[test,]

if(sum(ref.diff)!=0){ ## something to do
  
######### this bad codes does not work for inserions or deletions 
#bad.codes<-(!(indels[ref.diff,"REF"] %in% names(IUPAC_CODE_MAP)) | !(indels[ref.diff,"ALT"] %in% names(IUPAC_CODE_MAP))) ## not a IUPAC id
## bad.codes2<-(nchar(indels[ref.diff,"REF"])==1 | nchar(indels[ref.diff,"REF"])==1 | indels[ref.diff,"REF"]!="-" | indels[ref.diff,"ALT"]!="-" ) ## An insert or del
## ## Can't handle indrts and del cause I have lost confortance with the real base, so only do Snps till this is fixed
## ##  [1,] "chr10" "61846"  "61846"  "+"    "snp" "chr10:61846"                        "C" "CTTT"          "-" "AAA"          "C"        
## ##  [2,] "chr10" "72170"  "72170"  "+"    "snp" "rs201275645"                        "C" "CT"            "-" "T"            "C"        
## ##  [3,] "chr10" "79178"  "79178"  "+"    "snp" "rs201215709"                        "A" "AT"            "-" "T"            "A"        
## ##  [4,] "chr10" "79179"  "79179"  "+"    "snp" "rs199708591"                        "T" "TG"            "-" "G"            "T"        
## ##  [5,] "chr10" "95429"  "95429"  "+"    "snp" "rs35442274"                         "C" "CA"            "-" "A"            "C"


## bad.codes<-bad.codes1 & !bad.codes2

## bad.codes.ref<-indels[ref.diff,"REF"]=="-"
## bad.codes.alt<-indels[ref.diff,"ALT"]=="-"
## bad.codes<-bad.codes.ref | bad.codes.alt
## sum(bad.codes)
## cbind(indels,all.genomic)[ref.diff,][bad.codes,]


#ref.diff[ref.diff] <- !bad.codes ## don't do ones that have bad codes


to.convert.ref<-as.character(indels[ref.diff,"REF"])
to.convert.alt<-as.character(indels[ref.diff,"ALT"])

one.base<-nchar(to.convert.ref)==1 &  nchar(to.convert.alt)==1
sum(one.base)
sum(ref.diff)

  
to.convert.ref[one.base]<-flip.one.base(to.convert.ref[one.base])       # handles "-"          
to.convert.alt[one.base]<-flip.one.base(to.convert.alt[one.base])      # handles "-"      

#cbind(forward[an.indel,],to.convert.ref[an.indel],to.convert.alt[an.indel],all.genomic.ori[ref.diff][an.indel],indels.ori[ref.diff,][an.indel,])

# icon<-1
if(length(to.convert.ref[!one.base])!=0){
for(icon in 1:length(to.convert.ref[!one.base])){
to.convert.ref[!one.base][icon]<-as.character(reverseComplement(DNAString(to.convert.ref[!one.base][icon])))
to.convert.alt[!one.base][icon]<-as.character(reverseComplement(DNAString(to.convert.alt[!one.base][icon])))
}
} ## more than one multibase



indels[ref.diff,"REF"]<-to.convert.ref
indels[ref.diff,"ALT"]<-to.convert.alt


} # somthing to do in step 2

########### (3) ##########################################################################################
#Migh have amajor / mino problem still 
ref.diff<-indels[,"REF"]!=all.genomic
sum(ref.diff)

if(sum(ref.diff)!=0){ ## something to do


to.flip[ref.diff]<-to.flip[ref.diff]-1

## sum(ref.diff) #4601
print(paste(sum(ref.diff),":STEP 3  minor allele not REF after complmenet swap these"))
## cbind(indels,all.genomic)[ref.diff,][1:10,]

### swap REF AND ALT to see is just minor allele / major allele swap
the.ref<-indels[ref.diff,"REF"]
indels[ref.diff,"REF"]<-indels[ref.diff,"ALT"]
indels[ref.diff,"ALT"]<-the.ref
###### swap 0/1 --- keep track of times swapped

############ for indels must chenge start and end and the get new genomic sequence
an.indel<- (nchar(as.character(indels[ref.diff,"REF"])) >1) | (nchar(as.character(indels[ref.diff,"ALT"])) > 1) | indels[ref.diff,"REF"]=="-" | indels[ref.diff,"ALT"]=="-"
sum(an.indel)

if(sum(an.indel)>0){
indels[ref.diff,][an.indel,c("chr","start","end","REF","ALT")]<-redo.start.end.annovar(indels[ref.diff,][an.indel,c("chr","start","end","REF","ALT")])
no.chr<-check.chr.and.positions(indels[ref.diff,][an.indel,c("chr","start","end","REF","ALT")],chr.lengths)
the.chrom<-as.character(indels[ref.diff,"chr"][an.indel][!no.chr])
starts<-as.integer(as.character(indels[ref.diff,"start"][an.indel][!no.chr]))
ends<-as.integer(as.character(indels[ref.diff,"end"][an.indel][!no.chr]))
all.genomic.indels<-getSeq(Hsapiens, the.chrom, starts, ends,as.character=TRUE)
all.genomic[ref.diff][an.indel][!no.chr]<-all.genomic.indels
}

#cbind(all.genomic[!no.chr][ref.diff][an.indel],all.genomic.indels)
#### TOP /BOT diff



} # something to do in step 3


## ref.diff<-indels[,"REF"]!=all.genomic
## sum(ref.diff)

## an.indel<- (nchar(as.character(indels[ref.diff,"REF"])) >1) | (nchar(as.character(indels[ref.diff,"ALT"])) > 1) | indels[ref.diff,"REF"]=="-" | indels[ref.diff,"ALT"]=="-"
## sum(ref.diff[an.indel]) # 21
## forward.remain<-cbind(indels,all.genomic)[ref.diff,]
## save(list=c("forward.remain"),file="forward.remain.RData")
## cbind(forward.remain[an.indel,],all.genomic.ori[ref.diff][an.indel],indels.ori[ref.diff,][an.indel,])
## cbind(forward.remain,all.genomic.ori[ref.diff],indels.ori[ref.diff,],indels.ori[ref.diff,])

ref.diff<-indels[,"REF"]!=all.genomic
sum(ref.diff)


########### (4) ##########################################################################################
if(sum(ref.diff)!=0){ ## something to do
  print(paste(sum(ref.diff),":STEP 4 cant be fixed after MINOR/COMPLEMENT/MINOR COMPLEMENT complement back to original state "))

print(paste("Un-fixable snp file called:",paste("un.fixable",sample.files[i],"txt",sep=".")," put in :",snp.dir))
setwd(snp.dir)      
write.table(indels[ref.diff,"SNP"],file=paste("un.fixable",sample.files[i],"txt",sep="."),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
setwd(code.dir)
  
######### this bad codes does not work for inserions or deletions 
#bad.codes<-(!(indels[ref.diff,"REF"] %in% names(IUPAC_CODE_MAP)) | !(indels[ref.diff,"ALT"] %in% names(IUPAC_CODE_MAP))) ## not a IUPAC id
## bad.codes2<-(nchar(indels[ref.diff,"REF"])==1 | nchar(indels[ref.diff,"REF"])==1 | indels[ref.diff,"REF"]!="-" | indels[ref.diff,"ALT"]!="-" ) ## An insert or del
## ## Can't handle indrts and del cause I have lost confortance with the real base, so only do Snps till this is fixed
## ##  [1,] "chr10" "61846"  "61846"  "+"    "snp" "chr10:61846"                        "C" "CTTT"          "-" "AAA"          "C"        
## ##  [2,] "chr10" "72170"  "72170"  "+"    "snp" "rs201275645"                        "C" "CT"            "-" "T"            "C"        
## bad.codes.ref<-indels[ref.diff,"REF"]=="-"
## bad.codes.alt<-indels[ref.diff,"ALT"]=="-"
## bad.codes<-bad.codes.ref | bad.codes.alt
## sum(bad.codes)
## cbind(indels,all.genomic)[ref.diff,][bad.codes,]


#ref.diff[ref.diff] <- !bad.codes ## don't do ones that have bad codes


to.convert.ref<-indels[ref.diff,"REF"]
to.convert.alt<-indels[ref.diff,"ALT"]

one.base<-nchar(to.convert.ref)==1 &  nchar(to.convert.alt)==1
sum(one.base)
sum(ref.diff)

  
to.convert.ref[one.base]<-flip.one.base(to.convert.ref[one.base])       # handles "-"          
to.convert.alt[one.base]<-flip.one.base(to.convert.alt[one.base])      # handles "-"      

## cbind(forward.remain[an.indel,],to.convert.ref[an.indel],to.convert.alt[an.indel],all.genomic.ori[ref.diff][an.indel],indels.ori[ref.diff,][an.indel,])

# icon<-1
if(length(to.convert.ref[!one.base])!=0){
for(icon in 1:length(to.convert.ref[!one.base])){
to.convert.ref[!one.base][icon]<-as.character(reverseComplement(DNAString(to.convert.ref[!one.base][icon])))
to.convert.alt[!one.base][icon]<-as.character(reverseComplement(DNAString(to.convert.alt[!one.base][icon])))
}
} ## more than one multibase



indels[ref.diff,"REF"]<-to.convert.ref
indels[ref.diff,"ALT"]<-to.convert.alt


} # somthing to do in step 4



ref.diff<-indels[,"REF"]!=all.genomic
sum(ref.diff)

########### (5) ##########################################################################################
if(sum(ref.diff)!=0){ ## something to do
#cbind(forward.remain,all.genomic.ori[ref.diff],indels[ref.diff,],all.genomic[ref.diff],indels.ori[ref.diff,])


## exm-IND15-41876706 out by one base +1 CT genomic wanted? 
## exm-IND16-87694958
## think both above are just denovo insertions

#### OF ones I have dome before like above they are either polymorphic OR the are native insertions and so am trying to
#### match the new ALT alleles with the reference which will never match. So in these cases make sure any indels are in fact INSERTIONS
### the insertion IS assumed now to be on the forward stand, it be be different which would change the protein....
### swap REF AND ALT to see is just minor allele / major allele swap

is.del<-indels[,"ALT"]=="-"
ref.diff<-ref.diff & is.del
sum(ref.diff)

################ make remaining fails insertions
if(sum(ref.diff)!=0){
to.flip[ref.diff]<-to.flip[ref.diff]+1
### swap REF AND ALT to see is just minor allele / major allele swap
the.ref<-indels[ref.diff,"REF"]
indels[ref.diff,"REF"]<-indels[ref.diff,"ALT"]
indels[ref.diff,"ALT"]<-the.ref
###### swap 0/1 --- keep track of times swapped

############ for indels must chenge start and end and the get new genomic sequence
an.indel<- (nchar(as.character(indels[ref.diff,"REF"])) >1) | (nchar(as.character(indels[ref.diff,"ALT"])) > 1) | indels[ref.diff,"REF"]=="-" | indels[ref.diff,"ALT"]=="-"
sum(an.indel)
if(sum(an.indel)>0){
indels[ref.diff,][an.indel,c("chr","start","end","REF","ALT")]<-redo.start.end.annovar(indels[ref.diff,][an.indel,c("chr","start","end","REF","ALT")])
}
############ 5 make reamining fails insertions 


}
} # finish step 5
ref.diff<-indels[,"REF"]!=all.genomic
sum(ref.diff)
to.flip<-to.flip!=0
#sum(to.flip)

the.samples<-colnames(indels)[grepl(".GT$",colnames(indels))]
if(length(the.samples)>0 & sum(to.flip)>0){
    geno.to.flip<-indels[to.flip,the.samples]
    geno.to.flip[geno.to.flip=="0/0"]<-"9/9"
    geno.to.flip[geno.to.flip=="1/1"]<-"0/0"
    geno.to.flip[geno.to.flip=="9/9"]<-"1/1"
    indels[to.flip,the.samples]<-geno.to.flip
  }


rm(ref.ori)
rm(alt.ori)
rm(all.genomic)
rm(all.genomic.have)

## out.list<-list(column.labels,skip.lines,num.vars,info.types,info.class,info.description,format.types,format.class,format.description)
## names(out.list)<-c("column.labels","skip.lines","num.vars","info.types","info.class","info.description","format.types","format.class","format.description")
## out.list

## ################## Unwind using code like
## start.data<-prepare.for.Vcf.file.read(sample.files[isamp])
## for(i in 1:length(start.data)){assign( names(start.data)[i],value=start.data[[i]])}
## ###########################################

##      SNP                                                                   junk REF  ALT       chr     start       end         strand TYPE  SNP                                                                   junk REF  ALT 
## [1,] "UQDIexomechipWithAOGC_chr3:195510046:195510046:A:T-1_B_R_2111118932" "0"  "A"  "T"  "G"  "chr3"  "195510046" "195510046" "+"    "snp" "UQDIexomechipWithAOGC_chr3:195510046:195510046:A:T-1_B_R_2111118932" "0"  "A"  "T" 
## [2,] "UQDIexomechipWithAOGC_chr6:31238909:31238909:A:T-1_T_R_2111118934"   "0"  "T"  "A"  "G"  "chr6"  "31238909"  "31238909"  "+"    "snp" "UQDIexomechipWithAOGC_chr6:31238909:31238909:A:T-1_T_R_2111118934"   "0"  "T"  "A" 
## [3,] "UQDIexomechipWithAOGC_chr12:32794600:32794600:C:G-1_T_R_2111118516"  "0"  "C"  "G"  "A"  "chr12" "32794600"  "32794600"  "+"    "snp" "UQDIexomechipWithAOGC_chr12:32794600:32794600:C:G-1_T_R_2111118516"  "0"  "C"  "G" 
## [4,] "exm-IND15-41876706"                                                  "0"  "AG" "-"  "TC" "chr15" "44089415"  "44089416"  "+"    "snp" "exm-IND15-41876706"                                                  "0"  "AG" "-" 
## [5,] "exm-IND16-87694958"                                                  "0"  "-"  "AC" "C"  "chr16" "89167458"  "89167458"  "+"    "snp" "exm-IND16-87694958"                                                  "0"  "-"  "AC"
## [6,] "exm2216440"                                                          "0"  "A"  "T"  "C"  "chrM"  "14000"     "14000"     "+"    "snp" "exm2216440"                                                          "0"  "A"  "T" 
###### swap 0/1 --- keep track of times swapped
## cbind(indels,all.genomic)[ref.diff,][1:50,]  # 6
# indels
##        chr     start       end REF ALT origin MAF Lead.SNP Comment pval all.genomic
## 48431 chr3 195506099 195506099   G   C   rare 0.5      NaN   0,1,0   NA           T
## 5002  chr3 195508021 195508021   T   A   rare 0.5      NaN   0,1,0   NA           C
## 5173  chr3 195510046 195510046   T   A   rare 0.5      NaN   0,1,0   NA           G
## 5430  chr3 195514414 195514414   T   A   rare NaN      0.5   0,0,0   NA           G
## 7896  chr6  31238909  31238909   T   A   rare 0.5      NaN   0,1,0   NA           G

#}else{print("NO strand differences")} #skip doing the TOP/BOT conversion
