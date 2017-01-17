


#### scans through sub- directories making up a new bed file

### Assume
## snp list files have extension "_snp.list"
## evoker makes scores files of the same name but append ".score"
## beb file save under the same name as original bed with a prepend of snp list used.
## Assume gentain was the the only save bed file
## maybe "0" are bad genotype
## there is a common fam file in each sub-directory
## Need to extract out the SNPs from the list from the approprate bed file



ROOT.dir<-"/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/cluster_viz"
restrict.analysis<-FALSE # restrict.to.projects
ignore.dirs<-c()
 ## numof of lists processing
core.lists<-c("above.high","above.low","special","below.high","below.low","low")
core.beds.roots<-c("AOGC_gentrain","AOGC_opticall") # orginal prefix  used by Evoker add .chr.bed for real name to exclde orginal from analysis
#core.fams<-c("AOGC_gentrain.fam","AOGC_opticall.fam") # orginal ones used by Evoker
the.bim.root<-"AOGC_gentrain" # .chr.bim added later
the.fam<-"AOGC_gentrain.fam" # same fam no matter what the chr
bed.user.name.delim<-"_" ## file name "processor"_"list.name"_"chromosome"


projects<-list.dirs(ROOT.dir)

projects<-projects[!(projects %in% ignore.dirs)]
projects<-projects[!(projects %in% ROOT.dir)] ## remove the root
projects
##

all.data<-{}
count<-0


if(restrict.analysis){
  projects<-projects[projects %in% restrict.to.projects]
}

print(projects)
targets.file.ori<-""
all.beds<-{}
all.scores<-{}
all.snp.lists<-{}
####### ip<-1
for(ip in 1:length(projects)){
print(projects[ip])
on.target.stats<-{}  # this will be rebuilt from the individual

the.chr<-basename(projects[ip])
#the.dir<-dirname(projects[ip])

setwd(projects[ip])
files<-list.files(getwd())
files
snp.lists<-files[grep("_snp.list$",files)]
names(snp.lists)<-gsub("_snp.list$","",snp.lists)
names(snp.lists)<-gsub(paste(".",the.chr,sep=""),"",names(snp.lists))

snp.lists
scores<-files[grep("snp.list.scores$",files)]
names(scores)<-gsub("_snp.list.scores$","",scores)
names(scores)<-gsub(paste(".",the.chr,sep=""),"",names(scores))
scores

core.beds<-c(paste(core.beds.roots,".",gsub("chr","",the.chr),".bed",sep="") , paste(core.lists,".bed",sep="")) # original AND those created by below analsys
beds<-files[grep(".bed$",files)]
beds<-beds[!(beds %in% core.beds) & !grepl("^all.lists",beds) ] # make some "all.list.bed"s in below too!
names(beds)<-gsub("_snp.list.scores$","",beds)
names(beds)<-unlist(lapply(names(beds), function(x){unlist(strsplit(x,split=bed.user.name.delim))[2]}))

print(paste("Chromosome",the.chr,"---------------------------------------------------------"))
if( (length(scores)!=length(snp.lists)) | (sum(!(names(scores) %in% names(snp.lists)))!=0)){ # missing a score
  all.names<-unique(c(names(scores), names(snp.lists)))
   print(paste("Chromosome",the.chr))
   print(paste("missing snp list:", all.names[!(all.names %in% names(snp.lists))])) 
   print(paste("missing scores:", all.names[!(all.names %in% names(scores))]))
}
  
if( (length(beds)!=length(snp.lists)) | (sum(!(names(beds) %in% names(snp.lists)))!=0)){ # missing a bed or incorrebly names
  all.names<-unique(c(names(beds), names(snp.lists)))
   print(paste("Chromosome",the.chr))
   print(paste("missing beds:", all.names[!(all.names %in% names(beds))]))
   print(paste("Additional beds:", names(beds)[!(names(beds) %in% core.lists)]))

}

assign(paste("beds",the.chr,sep="."),value=beds)
assign(paste("scores",the.chr,sep="."),value=scores)
assign(paste("snp.lists",the.chr,sep="."),value=snp.lists)

all.beds<-c(all.beds,paste("beds",the.chr,sep="."))
all.scores<-c(all.scores,paste("scores",the.chr,sep="."))
all.snp.lists<-c(all.snp.lists,paste("snp.lists",the.chr,sep="."))

} ## loop over ip "prjects==chr

assign(data.names[i],value=data)
}
##############



############ for each snp.list extract out just tghose genotypes
beds
scores
snp.lists
print(projects)

# ip<-22

maybe.snps<-{}
final.plink.files<-{}

for(ip in 1:length(projects)){
print(projects[ip])
setwd(projects[ip])


the.chr<-basename(projects[ip])
the.chr.int<-gsub("^chr","",the.chr)
the.bim<-paste(the.bim.root,the.chr.int,"bim",sep=".")

beds<-eval(as.name(paste("beds",the.chr,sep=".")))
scores<-eval(as.name(paste("scores",the.chr,sep=".")))
snp.lists<-eval(as.name(paste("snp.lists",the.chr,sep=".")))

beds
scores
snp.lists
core.lists
getwd()
# ic<-1
maybe.snps.chr<-{}
plink.files.made<-{}
for(ic in 1:length(core.lists)){
  the.list<-core.lists[ic]
  num.snps.tested<-eval(system(paste("wc","-l",snp.lists[the.list],sep=" "),intern=TRUE)) ## if no snps plinkwill crash
  if(grepl("^0",num.snps.tested)){next}
  
  system(paste("plink","--bed",beds[the.list],"--bim",the.bim,"--fam",the.fam,"--extract",snp.lists[the.list],"--make-bed","--out",the.list,"--noweb",sep=" "))

  plink.files.made<-c(plink.files.made,paste(the.list,c("bed","bim","fam"),sep=".",collapse=" "))
  
  a.score<- read.table(scores[the.list],header=FALSE,fill=TRUE,stringsAsFactors=FALSE)
  maybe.snps.chr<-c(maybe.snps.chr,a.score[a.score[,2]==0,1])
}
maybe.snps<-c(maybe.snps,maybe.snps.chr)
write.table(maybe.snps.chr,file=paste("maybe.snps",the.chr,sep="."),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(plink.files.made,file=paste("plink.files.made",the.chr,sep="."),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

system( paste("plink","--merge-list",paste("plink.files.made",the.chr,sep="."),"--make-bed","--out","all.lists","--noweb",sep=" ") )
system( paste("plink","--bfile","all.lists","--exclude",paste("maybe.snps",the.chr,sep="."),"--make-bed","--out",paste("all.lists",the.chr,sep="."),"--noweb",sep=" ") )
final.plink.files<-c(final.plink.files,paste(paste(projects[ip],"/",paste("all.lists",the.chr,sep="."),sep=""),c("bed","bim","fam"),sep=".",collapse=" "))
}

final.plink.files
setwd(ROOT.dir)
write.table(final.plink.files,file="final.plink.files",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(maybe.snps,file="bad.checked.snps",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
system( paste("plink","--merge-list","final.plink.files","--make-bed","--out","all.good.evoker","--noweb",sep=" ") )

#all.good.evoker ## contains the final good checked snps
## Merging 23 samples, final sample contains 8015 individuals and 1404 markers
## Before frequency and genotyping pruning, there are 1404 SNPs
## 8015 founders and 0 non-founders found
## Total genotyping rate in remaining individuals is 0.991828
## 0 SNPs failed missingness test ( GENO > 1 )
## 0 SNPs failed frequency test ( MAF < 0 )
## After frequency and genotyping pruning, there are 1404 SNPs
## After filtering, 0 cases, 0 controls and 8015 missing
## After filtering, 0 males, 0 females, and 8015 of unspecified sex
## Writing pedigree information to [ all.good.evoker.fam ] 
## Writing map (extended format) information to [ all.good.evoker.bim ] 
## Writing genotype bitfile to [ all.good.evoker.bed ] 
## Using (default) SNP-major mode

1)## Move the data to all.good.evoker and bad snps to
"/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/cluster_viz" -> /media/scratch2/AOGC-NGS/ExomeChip


2) remove bad a new good snps from original data :

cut -f 2 all.good.evoker.bim > good.evoker.snps

cat good.evoker.snps  bad.checked.snps > evoker.checked.snps

plink --bfile  Exome_Plus_AOGC_Gentrain_TOP --exclude evoker.checked.snps --remove excluded.samples.txt  --make-bed --out Exome_Plus_AOGC_Gentrain_TOP.no.evoker.clean

plink --bfile Exome_Plus_AOGC_Gentrain_TOP.no.evoker.clean --bmerge all.good.evoker.bed all.good.evoker.bim all.good.evoker.fam --make-bed --out Exome_Plus_AOGC_Gentrain_TOP_corrected_clean

## 8015 founders and 0 non-founders found
## 28998 heterozygous haploid genotypes; set to missing
## Writing list of heterozygous haploid genotypes to [ Exome_Plus_AOGC_Gentrain_TOP_corrected_clean.hh ]
## 401 SNPs with no founder genotypes observed

#15 Oct decided to use zCall
plink --bfile  zCall_AOGC --exclude evoker.checked.snps --remove excluded.samples.txt  --make-bed --out  zCall_AOGC_TOP.no.evoker.clean

plink --bfile all.good.evoker --remove excluded.samples.txt  --make-bed --out all.good.evoker.clean

plink --bfile zCall_AOGC_TOP.no.evoker.clean --bmerge all.good.evoker.clean.bed all.good.evoker.clean.bim all.good.evoker.clean.fam --make-bed --out zCall_AOGC_TOP.with.evoker_corrected_clean


Before frequency and genotyping pruning, there are 266163 SNPs
7665 founders and 0 non-founders found
29945 heterozygous haploid genotypes; set to missing
Writing list of heterozygous haploid genotypes to [ zCall_AOGC_TOP.with.evoker_corrected_clean.hh ]
404 SNPs with no founder genotypes observed

After frequency and genotyping pruning, there are 266163 SNPs
After filtering, 0 cases, 0 controls and 7665 missing
After filtering, 0 males, 0 females, and 7665 of unspecified sex


######################## MOVE data to the forward strand:

setwd("/media/scratch2/AOGC-NGS/ExomeChip")

# bim<-read.table("Exome_Plus_AOGC_Gentrain_TOP_corrected_clean.bim",header=F,fill=TRUE,skip=0,stringsAsFactors=FALSE)
bim<-read.table("zCall_AOGC_TOP.with.evoker_corrected_clean.bim",header=F,fill=TRUE,skip=0,stringsAsFactors=FALSE)

dim(bim)
bim[1:5,]
colnames(bim)<-c("chr","SNP","cm","POS","A1","A2")


manifest<-read.table("HumanExome-12v1-1_A.csv",header=T,fill=TRUE,skip=7,stringsAsFactors=FALSE,sep=",")
dim(manifest)
## manifest<-manifest[,c("IlmnID","Name","IlmnStrand","SNP","Chr","MapInfo","SourceStrand")]
## manifest[1:5,]
## unique(manifest[,"SourceStrand"])

orient<-strsplit(manifest[,"IlmnID"],split="_")
ill.st<-unlist(lapply(orient,function(x) x[2]))
hg19.st<-unlist(lapply(orient,function(x) x[3]))
unique(ill.st)
unique(hg19.st)

controls<-is.na(hg19.st) ## all were controls 
sum(controls)

manifest<-manifest[!controls,]
 ## bim file has versions in snp names
#ver2 ver 3

vers<-grepl("_ver",bim[,"SNP"])
original.vers<-bim[vers,"SNP"]
sum(vers)
## test<-bim[vers,][1:20,2]
## sub("_ver\\d","",test,perl=TRUE)
bim[vers,"SNP"]<-sub("_ver\\d","",bim[vers,"SNP"],perl=TRUE)
bim[vers,][1:5,]

posns<-match(bim[,"SNP"],manifest[,"Name"])
missing<-is.na(posns)
sum(missing)


all.data<-cbind(bim,manifest[posns,])
dim(bim)
dim(all.data)

length(posns)
dim(all.data)

all.data[1:2,]

missing.in.manifest<-missing
uqdi<-grepl("^UQDI",all.data[,"SNP"])
sum(uqdi)
sum(missing.in.manifest)
sum(missing.in.manifest & !uqdi) # 0 all are uqdi

all.data[uqdi,][1:5,]

orient<-strsplit(all.data[,"IlmnID"],split="_")
orient[[1]]
ill.st<-unlist(lapply(orient,function(x) x[2]))
hg19.st<-unlist(lapply(orient,function(x) x[3]))
unique(ill.st) # "M" "P" "T" "B"
unique(hg19.st) # "R" "F"
sum(is.na(hg19.st))
sum(is.na(ill.st))


###need to do uqdi seperately
orient<-strsplit(all.data[uqdi,"SNP"],split="_")
orient[[1]]
ill.st.uqdi<-unlist(lapply(orient,function(x) x[3]))
hg19.st.uqdi<-unlist(lapply(orient,function(x) x[4]))
unique(ill.st.uqdi) # "M" "P" "T" "B"
unique(hg19.st.uqdi) # "R" "F"
sum(is.na(hg19.st.uqdi))
sum(is.na(ill.st.uqdi))

ill.st[uqdi]<-ill.st.uqdi
hg19.st[uqdi]<-hg19.st.uqdi
sum(is.na(hg19.st.uqdi))
sum(is.na(ill.st.uqdi))

all.data<-cbind(all.data,ill.st,hg19.st)

all.data[1:5,]
unique(all.data[,"ill.st"]) # "M" "P" "T" "B"
unique(all.data[,"hg19.st"])




man.uqdi<-read.table("/media/Bioinform-D/Research/exome Chip/Final targets/UQDIexomechipWithAOGC_By ExistingDesign - WGGT Final_103877.score.csv",header=T,sep=",",fill=TRUE,skip=20,stringsAsFactors=FALSE)

man.uqdi[1:5,]
all.data[uqdi,][1:5,]

keep.cols<- c("Genome_Build_Version","Chromosome","Coordinate","Final_Score","Sequence","Sequence_Orientation","Validation_Class","Assay_Design_Id","Ilmn_Id","Normalization_Bin","Bead_Types.Assay","Design_Date","Assay_Type")

unique(man.uqdi[,"Sequence_Orientation"]) ## all forward

posns<-match(all.data[,"SNP"],man.uqdi[,"Locus_Name"])
missing<-is.na(posns)
sum(!missing)

uqdi.manifest.details<-man.uqdi[posns,keep.cols]

all.data<-cbind(all.data,uqdi.manifest.details)
all.data[uqdi,][1:5,]
man.uqdi[1,]

##### make both the forward and top the same for uqdi so won't otherwise flip these
sum(is.na(all.data[uqdi,"Sequence"]))
all.data[uqdi,"SourceSeq"]<-all.data[uqdi,"Sequence"]
all.data[uqdi,"TopGenomicSeq"]<-all.data[uqdi,"Sequence"]
###



############################################ check is Source seq and top Se qre different

indel.alleles<-all.data[,"TopGenomicSeq"]
indel.alleles<-gsub("[","::",indel.alleles,fixed=TRUE)
indel.alleles<-gsub("]","::",indel.alleles,fixed=TRUE)
indel.alleles<-gsub("/","::",indel.alleles,fixed=TRUE)

orient<-strsplit(indel.alleles,split="::")
orient[[1]]
a1.hg19.top<-unlist(lapply(orient,function(x) x[2]))
a2.hg19.top<-unlist(lapply(orient,function(x) x[3]))
unique(a1.hg19) ## a1 is always the deletion
unique(a2.hg19) ## contains alleles *the insertion*


indel.alleles<-all.data[,"SourceSeq"]
indel.alleles<-gsub("[","::",indel.alleles,fixed=TRUE)
indel.alleles<-gsub("]","::",indel.alleles,fixed=TRUE)
indel.alleles<-gsub("/","::",indel.alleles,fixed=TRUE)

orient<-strsplit(indel.alleles,split="::")
orient[[1]]
a1.hg19.frd<-unlist(lapply(orient,function(x) x[2]))
a2.hg19.frd<-unlist(lapply(orient,function(x) x[3]))
unique(a1.hg19) ## a1 is always the deletion
unique(a2.hg19) ## contains alleles *the insertion*



to.flip.al<- a1.hg19.frd != a1.hg19.top |  a2.hg19.frd != a2.hg19.top 


### IF ON TOP already
# B & F flip
# B & R ok 
# T & F ok
# T & R flip

# M & F flip
# M & R ok 
# P & F ok
# P & R flip
to.flip<- (all.data[,"ill.st"]=="B" & all.data[,"hg19.st"]=="F") | (all.data[,"ill.st"]=="T" & all.data[,"hg19.st"]=="R") |
             (all.data[,"ill.st"]=="M" & all.data[,"hg19.st"]=="F") |  (all.data[,"ill.st"]=="P" & all.data[,"hg19.st"]=="R")
length(to.flip)
sum(to.flip)

sum(!to.flip & to.flip.al) #775
#Some of these DO need flipping SO the T_F like annotations are a problem
## $$$ Not that to.flip.al does not always identify ones that need flipping cause to C/G A/T (same problem as plink flipping)snps eg

##        chr    SNP cm    POS A1 A2                  IlmnID   Name IlmnStrand   SNP AddressA_ID                                   AlleleA_ProbeSeq AddressB_ID
## 84583    1 exm158  0 878744  G  C exm158-0_B_F_1921504479 exm158        BOT [G/C]    32734962 TGGACCGTGGATGACGTCTGCAGCTTCGTGGGGGGCCTGTCTGGCTGTGG    47757125
## 151215   1 exm355  0 889403  T  A exm355-0_T_R_2058862089 exm355        TOP [A/T]    88736127 CCGTCCAGCAGCCCGCTCTGGGGGAAGCTTCGTGTGGACATCAAGGCTTA    31637906

##                                          AlleleB_ProbeSeq GenomeBuild Chr MapInfo  Ploidy      Species    Source SourceVersion SourceStrand
## 84583  TGGACCGTGGATGACGTCTGCAGCTTCGTGGGGGGCCTGTCTGGCTGTGC        37.1   1  878744 diploid Homo sapiens ExomeSNPs             0          BOT
## 151215 CCGTCCAGCAGCCCGCTCTGGGGGAAGCTTCGTGTGGACATCAAGGCTTT        37.1   1  889403 diploid Homo sapiens ExomeSNPs             0          BOT
##                                                                                                                            SourceSeq
## 84583  ACGTCACCAAGTGGACCGTGGATGACGTCTGCAGCTTCGTGGGGGGCCTGTCTGGCTGTG[C/G]AGAGTACACTCGGGTAAGGGGGGGCCCCAGTTCCTGGGGCGGGGCTGGAGCTGGCTGGCA
## 151215 GATACAGCCAGGCCCCCGTGCCCTCCCCACCAGAATAGCACCTGTATGGCCGAGCCCAGG[A/T]AAGCCTTGATGTCCACACGAAGCTTCCCCCAGAGCGGGCTGCTGGACGGCTGCAGCATCC
##                                                                                                                        TopGenomicSeq BeadSetID Exp_Clusters RefStrand ill.st hg19.st
## 84583  TGCCAGCCAGCTCCAGCCCCGCCCCAGGAACTGGGGCCCCCCCTTACCCGAGTGTACTCT[C/G]CACAGCCAGACAGGCCCCCCACGAAGCTGCAGACGTCATCCACGGTCCACTTGGTGACGT
##151215 GGATGCTGCAGCCGTCCAGCAGCCCGCTCTGGGGGAAGCTTCGTGTGGACATCAAGGCTT[A/T]CCTGGGCTCGGCCATACAGGTGCTATTCTGGTGGGGAGGGCACGGGGGCCTGGCTGTATC
## test<-to.flip & !to.flip.al & !uqdi
## sum(test)
## all.data[test,][1:5,]
## tapply(all.data[test,"SourceStrand"],all.data[test,"SourceStrand"],length)


to.flip<-to.flip | to.flip.al
sum(to.flip) # 132642



cbind(original.vers,bim[vers,"SNP"])[1:30,]
corrected.vers<-bim[vers,"SNP"]
bim[vers,"SNP"]<-original.vers ##  so can flip with existing bi file
all.data[vers,"SNP"]<-as.character(original.vers) ##  so can flip with existing bi file
sum(to.flip & uqdi) # 11746
all.data[vers,][1:5,]
#############################
all.data[1:5,]
#now correct alleles for indels
unique(all.data[,"A1"]) # "A" "G" "C" "T" "0" "D" "I"
unique(all.data[,"A2"])

a.indel<- (all.data[,"A1"] %in% c("I","D")) | (all.data[,"A2"] %in% c("I","D"))

sum(a.indel) # 136
sum(a.indel & to.flip) # 1 (use TOP alleles so match SNPs and flip later) !
cbind(all.data[a.indel,],to.flip[a.indel],to.flip.al[a.indel])[1:13,]

################ GET INDEL ALLELES ###################
## exm-rs3834129 is the forward AGTAAG 4820
##                      TOP is       CTTACT

indel.alleles<-all.data[a.indel,"TopGenomicSeq"]
indel.alleles<-gsub("[","::",indel.alleles,fixed=TRUE)
indel.alleles<-gsub("]","::",indel.alleles,fixed=TRUE)
indel.alleles<-gsub("/","::",indel.alleles,fixed=TRUE)

orient<-strsplit(indel.alleles,split="::")
orient[[1]]
a1.hg19<-unlist(lapply(orient,function(x) x[2]))
a2.hg19<-unlist(lapply(orient,function(x) x[3]))
unique(a1.hg19) ## a1 is always the deletion
unique(a2.hg19) ## contains alleles *the insertion allele*

#### all these indels are insertions but the look
#### all are on the forward strand except one.

all.data[a.indel,][1:5,]


###### get a1 and a2 ordering which is based on the frequency
##### whe you swap indels for start and end " end can be affected"
unique(all.data[a.indel,"A1"]) #  "D" "I" "0"
unique(all.data[a.indel,"A2"]) # "I" "D"
a2.is.deletion<-all.data[a.indel,"A2"]=="D"
sum(a2.is.deletion) #43
all.data[a.indel,][a2.is.deletion,][1:5,]

a1.hg19[a2.is.deletion]<-a2.hg19[a2.is.deletion]
a2.hg19[a2.is.deletion]<-"-"

#a1.hg19 a2.hg19 are top strand alleles- no indels needed flipping so this is safe
sum(a.indel & to.flip) # 0 (all on Top start no flips needed !

A1.base<-all.data[,"A1"]
A2.base<-all.data[,"A2"]

A1.base[a.indel]<-a1.hg19
A2.base[a.indel]<-a2.hg19



unique(A1.base) #contains - and 0
unique(A2.base) #contains - and 0


#A1.base and A2.base CONTAIN THE TOP STAND ALLELES

all.data<-cbind(all.data,A1.base,A2.base)

all.data[a.indel,][1:5,]
test<-is.na(a1.hg19)
sum(test) #0
all.data[a.indel,][a2.is.deletion,][1:5,]

################# fix monomorphc calls ### use the forward stand and then flip as here have mono & to.flip uqdi dont have the forward stad alleles
unique(all.data[,"A1"]) #  "D" "I" "0"
unique(all.data[,"A2"]) # "I" "D"

sum(all.data[,"A1"]=="0") ## monomorphic call _> 426
sum(all.data[,"A2"]=="0") # 399 both alleles 0 in -> 402 in zcall

mono<-all.data[,"A1"]=="0" | all.data[,"A2"]=="0" # sum(all.data[,"A2"]=="0")
sum(mono) # 426 in zcall
all.data[mono,][1:5,]
sum(mono)
sum(mono & to.flip) #3271 -> zcall 232
sum(mono & to.flip & a.indel) # 0 no nomomorphic indels

all.data[1:5,]
indel.alleles<-all.data[mono,"SourceSeq"]
sum(is.na(indel.alleles)) #586 -> zcall 253
indel.alleles[is.na(indel.alleles)]<-all.data[mono,][is.na(indel.alleles),c("Sequence")]

sum(is.na(indel.alleles)) #0

indel.alleles<-gsub("[","::",indel.alleles,fixed=TRUE)
indel.alleles<-gsub("]","::",indel.alleles,fixed=TRUE)
indel.alleles<-gsub("/","::",indel.alleles,fixed=TRUE)

orient<-strsplit(indel.alleles,split="::")
orient[[1]]
a1.hg19<-unlist(lapply(orient,function(x) x[2]))
a2.hg19<-unlist(lapply(orient,function(x) x[3]))
unique(a1.hg19) ## a1 is always the deletion
unique(a2.hg19) ## contains alleles *the insertion*

###these I think are all forward strand
length(mono)
length(to.flip)
dim(all.data)
length(a1.hg19)
sum(mono)


all.data[a.indel,][1:5,]

A1.base.f.m<-all.data[,"A1"]
A2.base.f.m<-all.data[,"A2"]

A1.base.f.m[mono]<-a1.hg19  
A2.base.f.m[mono]<-a2.hg19

### above are foward strand so flip to TOP
cbind(all.data,A1.base.f.m,A2.base.f.m)[mono & to.flip,][1:5,]

#### these alles are derived from the forward strand data ( cause only had that for uqdi) so pinto TOP - no indels so just use the one base

flip.one.base<-function(to.flip){
  ### flip one base only to.flip is just an vector of A,T,G,C's
  As<-to.flip=="A"
  Ts<-to.flip=="T"
  Gs<-to.flip=="G"
  Cs<-to.flip=="C"
 # Dels<-to.flip=="-"
  
   to.flip[As]<-"T"
    to.flip[Ts]<-"A"
    to.flip[Gs]<-"C"
    to.flip[Cs]<-"G"
 #  to.flip[Dels]<-"-"
  to.flip
}

A1.base.f.m[mono & to.flip]<-flip.one.base(A1.base.f.m[mono & to.flip])
A2.base.f.m[mono & to.flip]<-flip.one.base(A2.base.f.m[mono & to.flip])
cbind(all.data,A1.base.f.m,A2.base.f.m)[mono & to.flip,][1:5,] ### looks good
sum(mono & to.flip & a.indel) #0 so flip one base is ok

### check have at lease one base
ok<-(A1.base.f.m[mono] == all.data[mono,"A1.base"] | A1.base.f.m[mono] == all.data[mono,"A2.base"] |
     A2.base.f.m[mono] == all.data[mono,"A1.base"] | A2.base.f.m[mono] == all.data[mono,"A2.base"] |
     (all.data[mono,"A1.base"]=="0" & all.data[mono,"A2.base"]=="0")
      )

sum(!ok) # a1 match or a2 match or both 0 
cbind(all.data,A1.base.f.m,A2.base.f.m)[mono,][ ok,][1:5,]




#1.base.f.m and A1.base.f.m now on TOP strand  ### reorder so have correct order
do.nothing<-(all.data[mono,"A1.base"]=="0" & all.data[mono,"A2.base"]=="0")
a1.to.a2<-(all.data[mono,"A1.base"]=="0" & A2.base.f.m[mono] == all.data[mono,"A2.base"] )  & !do.nothing
sum(a1.to.a2) #16
cbind(all.data,A1.base.f.m,A2.base.f.m)[mono,][a1.to.a2,][1:5,]
all.data[mono,][a1.to.a2,"A1.base"]<-A1.base.f.m[mono][a1.to.a2]

a1.to.a1<-(all.data[mono,"A1.base"]=="0" & A1.base.f.m[mono] == all.data[mono,"A2.base"] )  & !do.nothing
cbind(all.data,A1.base.f.m,A2.base.f.m)[mono,][a1.to.a1,][1:5,]
all.data[mono,][a1.to.a1,"A1.base"]<-A2.base.f.m[mono][a1.to.a1]


sum(all.data[mono,"A2.base"]=="0" & !do.nothing) ## no a2 0 and a1 not
cbind(all.data,A1.base.f.m,A2.base.f.m)[mono,][do.nothing,][1:5,]


########### pUT missing mono alleles in TOP back into A1.base
all.data[mono,][do.nothing,"A1.base"]<-A1.base.f.m[mono][do.nothing]

all.data<-as.matrix(all.data)
all.data[mono,][do.nothing,"A2.base"]<-A2.base.f.m[mono][do.nothing]

test<-all.data[,"A1.base"]=="0"
sum(test)
all.data[1:5,]


#check no "0"
unique(all.data[,"A1.base"]) #  
unique(all.data[,"A2.base"]) # 


all.data[uqdi,][1:5,]

## > getwd()
## [1] "/media/scratch2/AOGC-NGS/ExomeChip"

gsub("^\\s+","","  fg3456fg")


all.data[,"POS"]<-gsub("^\\s+","",all.data[,"POS"])
sum(grepl("^ ",all.data[,"POS"]))

all.data[,"chr"]<-gsub("^\\s+","",all.data[,"chr"])
sum(grepl("^ ",all.data[,"chr"]))
sum(grepl("^ ",all.data[,"A1.base"]))

write.table(all.data,file="SNP_annotation_summary_technical_zCall.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(all.data[to.flip,"SNP"],file="flip_top_to_forward_zCall.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)


ann<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP_annotation_summary_technical_zCall_final.txt",header=T,fill=TRUE,skip=0,stringsAsFactors=FALSE)


all.data[1:5,]
write.table(all.data[,c("chr","SNP","cm","POS","A1.base","A2.base")],file="zCall_AOGC_TOP.no.evoker.clean_Alleles.bim",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

sum(to.flip)
sum(grepl("ver",all.data[,c("SNP")]))

############ test was kept in order
bim<-read.table("zCall_AOGC_TOP.with.evoker_corrected_clean.bim",header=F,fill=TRUE,skip=0,stringsAsFactors=FALSE)
dim(bim)
dim(all.data)
bim[1:5,]
colnames(bim)<-c("chr","SNP","cm","POS","A1","A2")
test<-all.data[,"SNP"]!=bim[,"SNP"]
sum(test)
cbind(all.data[,"SNP"],bim[,"SNP"])[test,][1:10,]

all.data[a.indel & to.flip,]



ins<-all.data[,"A2.base"]=="-"
sum(ins & to.flip)
test<-(ins & to.flip)
all.data[test,]

test<-"exm-rs3834129"
loc<-grep(test,all.data[,"SNP"])
loc

cbind(all.data,to.flip)[loc:(loc+1),]


#### at this point just had modify 
zCall_AOGC_TOP.no.evoker.clean_Alleles.bim


#############################################################

plink --bfile zCall_AOGC_TOP.with.evoker_corrected_clean --flip  flip_top_to_forward_zCall.txt --make-bed --out zCall_AOGC.with.evoker_corrected_clean

plink --bed zCall_AOGC_TOP.with.evoker_corrected_clean.bed --bim zCall_AOGC_TOP.no.evoker.clean_Alleles.bim --fam zCall_AOGC_TOP.with.evoker_corrected_clean.fam --flip  flip_top_to_forward_zCall.txt --make-bed --out zCall_AOGC.with.evoker_corrected_clean_Alleles

## 7665 founders and 0 non-founders found
## 29945 heterozygous haploid genotypes; set to missing
## Writing list of heterozygous haploid genotypes to [ zCall_AOGC.with.evoker_corrected_clean_Alleles.hh ]
## 404 SNPs with no founder genotypes observed
## Warning, MAF set to 0 for these SNPs (see --nonfounders)
## Writing list of these SNPs to [ zCall_AOGC.with.evoker_corrected_clean_Alleles.nof ]
## Total genotyping rate in remaining individuals is nan
## 0 SNPs failed missingness test ( GENO > 1 )
## 0 SNPs failed frequency test ( MAF < 0 )
## After frequency and genotyping pruning, there are 266163 SNPs
## After filtering, 0 cases, 0 controls and 7665 missing
## After filtering, 0 males, 0 females, and 7665 of unspecified sex
## Writing pedigree information to [ zCall_AOGC.with.evoker_corrected_clean_Alleles.fam ] 
## Writing map (extended format) information to [ zCall_AOGC.with.evoker_corrected_clean_Alleles.bim ] 
## Writing genotype bitfile to [ zCall_AOGC.with.evoker_corrected_clean_Alleles.bed ] 

# 71   2 exm-IND2-165344412 I and D are swapped


############ plink does NOT flip multi-allel bases!

## Options in effect:
## 	--bed zCall_AOGC_TOP.with.evoker_corrected_clean.bed
## 	--bim zCall_AOGC_TOP.no.evoker.clean_Alleles.bim
## 	--fam zCall_AOGC_TOP.with.evoker_corrected_clean.fam
## 	--flip flip_top_to_forward_zCall.txt
## 	--make-bed
## 	--out zCall_AOGC.with.evoker_corrected_clean_Alleles


## pleo@bioinform01:/media/scratch2/AOGC-NGS/ExomeChip$ grep exm-rs3834129 zCall_AOGC.with.evoker_corrected_clean_Alleles.bim
## 2	exm-rs3834129	0	202097532	-	CTTACT
## pleo@bioinform01:/media/scratch2/AOGC-NGS/ExomeChip$ grep exm-rs3834129 flip_top_to_forward_zCall.txt
## exm-rs3834129
## pleo@bioinform01:/media/scratch2/AOGC-NGS/ExomeChip$  grep exm-rs3834129 zCall_AOGC_TOP.no.evoker.clean_Alleles.bim
## 2	exm-rs3834129	0	202097532	-	CTTACT


##        SourceSeq                                                                                                                           
## 4820   "TTATGAATGAGCCGAGGAAGGCACTGAGACGTTAAGTAACTTGCCCAAGGTCACGCAGCT[-/AGTAAG]TGGCAGAGCAAGAATTACTATGGCTTTATAAGCCTAGGAAAAAGTCTGAAAGAATCAAAA"
##        TopGenomicSeq                                                                                                                        
## 4820   "TTTTGATTCTTTCAGACTTTTTCCTAGGCTTATAAAGCCATAGTAATTCTTGCTCTGCCA[-/CTTACT]AGCTGCGTGACCTTGGGCAAGTTACTTAACGTCTCAGTGCCTTCCTCGGCTCATTCATAA

## #### MANUALLY CHANGED BIM FILE TO FIX THIS OR NEED TO DO REVERSE COMPLEMENT IN THIS CODE #######

### RUN my annoation code and 1474 STILL WOONG with 6 not fixable
rm(forward)
load("/media/Bioinform-D/Research/AML sequencing/forward.RData")
dim(forward)

rm(forward.remain)
load("/media/Bioinform-D/Research/AML sequencing/forward.remain.RData")
dim(forward.remain)

ignore<-forward[,"SNP"] %in% forward.remain[,"SNP"]
sum(ignore)
forward<-forward[!ignore,]
dim(forward)
forward[1:6,]

posns<-match(all.data[,"SNP"],forward[,"SNP"])
missing<-is.na(posns)
forward<-forward[posns,]
more.flips<-!missing

multi.base<- (nchar(as.character(all.data[,"A1.forward"])) >1) | (nchar(as.character(all.data[,"A2.forward"])) > 1)
sum(multi.base)
sum(multi.base & more.flips) #0 would need to fix manually
sum(to.flip & more.flips)  # 1114 flipped in error
sum(to.flip.al & more.flips) #393 from my strand info

to.flip<-to.flip & more.flips
all.data[,"to.flip"]<-to.flip

write.table(all.data[more.flips,"SNP"],file="MORE_flip_top_to_forward_zCall.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)


plink --bfile zCall_AOGC.with.evoker_corrected_clean_Alleles --flip  MORE_flip_top_to_forward_zCall.txt --make-bed --out zCall_AOGC.with.evoker_corrected_clean_FINAL






bim.FLIPPED<-read.table("zCall_AOGC.with.evoker_corrected_clean_FINAL.bim",header=F,fill=TRUE,skip=0,stringsAsFactors=FALSE)

colnames(bim.FLIPPED)<-c("chr","SNP","cm","POS","A1","A2")
bim.FLIPPED[1:5,]
unique(bim.FLIPPED[,"chr"])

colnames(all.data)[colnames(all.data)=="A1.base"]<-"A1.TOP"
colnames(all.data)[colnames(all.data)=="A2.base"]<-"A2.TOP"

test<-all.data[,"SNP"]!=bim.FLIPPED[,"SNP"]
sum(test)
cbind(all.data[,"SNP"],bim.FLIPPED[,"SNP"])[test,][1:10,]

#orde moved about

posns<-match(all.data[,"SNP"],bim.FLIPPED[,"SNP"])
missing<-is.na(posns)
sum(missing) #0
bim.FLIPPED<-bim.FLIPPED[posns,]
test<-all.data[,"SNP"]!=bim.FLIPPED[,"SNP"]
sum(test)
A1.forward<-bim.FLIPPED[,"A1"]
A2.forward<-bim.FLIPPED[,"A2"]

if("A1.forward" %in% colnames(all.data)){
  all.data[,"A1.forward"]<-A1.forward
  all.data[,"A2.forward"]<-A2.forward
}else{

all.data<-cbind(all.data,A1.forward,A2.forward)
all.data<-cbind(all.data,to.flip)

}

#### 344 poly morphic  survive hwe Qc for imputation $$ check with UCSC and these look real)


all.data[1:5,]


#### HAs final flips with 
write.table(all.data,file="SNP_annotation_summary_technical_zCall_final.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
save(list=c("all.data"),file="SNP_annotation_summary_technical_zCall_final.RData")
#load("SNP_annotation_summary_technical_zCall_final.RData")
#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
 
## Annotated with /media/Bioinform-D/Research/AML sequencing/Annotate_sequence_run.r
## Collated with /media/Bioinform-D/Research/AML sequencing/Analyse_project_new.r
## Output in /media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/

setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip")

ann<-read.table("/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/Analysis/zCall_AOGC.with.evoker_corrected_clean_FINAL.analysis.txt",header=F,fill=TRUE,skip=0,stringsAsFactors=FALSE)

convert<-read.table("/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/exomeChip.to.rs.snp137.link.txt",header=T,fill=TRUE,skip=0,stringsAsFactors=FALSE)
convert[1:5,]
colnames(convert)[colnames(convert)=="ID..maf"]<-"rs"
missing<-is.na(convert[,"rs"])
sum(missing) # 16208
convert<-convert[!missing,]
dups<-duplicated(convert[,"rs"])
sum(dups)

dup.snps<-convert[dups,"rs"]

convert[convert[,"rs"] %in% dup.snps,]
pholymorphic.rs<-convert[convert[,"rs"] %in% dup.snps,]

pholymorphic.rs[,1:10]

write.table(pholymorphic.rs,file="exomeChip.to.rs.snp137.PolymorphicExceptions.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
convert<-convert[!(convert[,"rs"] %in% dup.snps),] # leave the dups as is


write.table(convert,file="exomeChip.to.rs.snp137.link.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(convert[,c("SNP","rs")],file="exomeChip.to.rs.snp137.convert.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

plink --bfile recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL --update-map exomeChip.to.rs.snp137.convert.txt --update-name --make-bed --out recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL_withRS



#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
 plink --bfile HBM_final_to_impute --flip to.flip.HBM_final_to_impute.bim.txt --make-bed --out HBM_final_to_impute_Best


### get.coomon.snp.in.bims.r to get common.snps.txt (polymorphic SNps werels
not recoded
plink --bfile ALL_eth_650Y_forward_hg19_Best --extract common_snps.txt  --make-bed --out ALL_eth_Best_Common
plink --bfile AOGC_merge_common._Best --extract common_snps.txt  --make-bed --out AOGC_merge_Best_Common
plink --bfile recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL_withRS --remove related.to.remove.recode.txt    --extract common_snps.txt  --make-bed --out AOGC_exomeChip_withRS

#### get rid of common ones from AOGC merge and exomechip

cut -d " "  -f 1,2   AOGC_exomeChip_withRS.fam > aogc.exome.samples
plink --bfile AOGC_merge_Best_Common --remove aogc.exome.samples --make-bed --out AOGC_merge_Best_Common_Unique
## 957 samples remail
 
plink --bfile  AOGC_exomeChip_withRS --bmerge AOGC_merge_Best_Common_Unique.bed AOGC_merge_Best_Common_Unique.bim AOGC_merge_Best_Common_Unique.fam --make-bed --out exomeChip.GWAS 
## error use to rs7824 which is ply morphic -remove
plink --bfile  AOGC_exomeChip_withRS --bmerge ALL_eth_Best_Common.bed ALL_eth_Best_Common.bim ALL_eth_Best_Common.fam --make-bed --out exomeChip.650Y
## error use to rs7824 which is ply morphic -remove
plink --bfile  AOGC_merge_Best_Common_Unique --bmerge ALL_eth_Best_Common.bed ALL_eth_Best_Common.bim ALL_eth_Best_Common.fam --make-bed --out GWAS.650Y

plink --bfile  AOGC_exomeChip_withRS --exclude exomeChip.GWAS.missnp --make-bed --out AOGC_exomeChip_withRS.1
plink --bfile  AOGC_merge_Best_Common_Unique --exclude exomeChip.GWAS.missnp --make-bed --out AOGC_merge_Best_Common_Unique.1
plink --bfile  ALL_eth_Best_Common --exclude exomeChip.GWAS.missnp --make-bed --out ALL_eth_Best_Common.1

#use recoded data to mege with AOGC

plink --bfile  AOGC_exomeChip_withRS.1 --bmerge AOGC_merge_Best_Common_Unique.1.bed AOGC_merge_Best_Common_Unique.1.bim AOGC_merge_Best_Common_Unique.1.fam --make-bed --out exomeChip.GWAS
plink --bfile  exomeChip.GWAS --bmerge ALL_eth_Best_Common.1.bed ALL_eth_Best_Common.1.bim ALL_eth_Best_Common.1.fam --make-bed --out exomeChip.GWAS.650Y


## ####################### Clean up for PCA 


plink --bfile  exomeChip.GWAS.650Y --allow-no-sex --maf 0.05 --max-maf 0.95 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out exomeChip.GWAS.650Y.f
./shellfish.py --pca --numpcs 10 --maxprocs 2 --file exomeChip.GWAS.650Y.f --out exomeChip.GWAS.650Y.f.1 

plink --bfile exomeChip.GWAS.650Y.f --keep exomeChip.GWAS.650Y.f.1.evecs_keep.txt --allow-no-sex --maf 0.05 --max-maf 0.95 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out exomeChip.GWAS.650Y.f.1

plink --bfile exomeChip.GWAS.650Y.f.1  --exclude Price2008_hg19_paul.txt  --range --make-bed --out exomeChip.GWAS.650Y.f.1.ld  --allow-no-sex

./shellfish.py --pca --numpcs 10 --maxprocs 7 --file exomeChip.GWAS.650Y.f.1.ld --out exomeChip.GWAS.650Y.f.2.ld



./shellfish.py --pca --numpcs 10 --maxprocs 7 --file exomeChip.GWAS.650Y.f.1 --out good
./shellfish.py --pca --numpcs 10 --maxprocs 7 --file exomeChip.GWAS.650Y.f.1.ld --out bad

plink --bfile exomeChip.GWAS.650Y.f.1.ld --extract snps.trial.txt  --make-bed --out snp.trial

./shellfish.py --pca --numpcs 10 --maxprocs 7 --file snp.trial --out snp.trial

##$$$$$ /media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/exomeChip.GWAS.650Y.f.2.evecs_pca_eigenvectors.txt ARE THE EIGENVECTORS

##%%%%%%%%%%%%%%%%%5 THE above is the last working solution
##%%%%%%%%%%%%%%%%%5
##%%%%%%%%%%%%%%%%%5
##%%%%%%%%%%%%%%%%%5
##%%%%%%%%%%%%%%%%%5
##%%%%%%%%%%%%%%%%%5
##%%%%%%%%%%%%%%%%%5
########################## So no matter what I do I can't get PCAS to work right
########################## tried remove snps with phenotype and genotype assoctionas
########################## GO back to rull dataset and leave in stratification outlyers
########################## test associations afterwarsds

 ############### some related samples were renamed reidentify related samples
key<-read.table("/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/pheno_gwas_id_missmatches.csv",header=T,fill=TRUE,skip=0,stringsAsFactors=FALSE)
rel<-read.table("/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/IBD/related.to.remove.txt",header=F,fill=TRUE,skip=0,stringsAsFactors=FALSE)

key[1:5,]
rel[1:5,]

posns<-match(rel[,1],key[,"PhenoID"])
missing<-is.na(posns)
rel[!missing,1]<-key[posns[!missing],"GWASID"]
rel[!missing,2]<-key[posns[!missing],"GWASID"]
rel[10,]
key[43,]

write.table(rel,file="/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/IBD/related.to.remove.recode.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
snps.with.probs.or.uqdi.txt

snps<-read.table("snps.with.probs.or.uqdi.txt",header=T,fill=TRUE,stringsAsFactors=FALSE)
all.snps<-unique(snps[,1])
length(all.snps)

write.table(all.snps,file="snps.with.probs.or.uqdi.or.targets.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)


exomeChip.GWAS.HBM.650Y.ld.f.b
snps.any.probs.txt
plink --bfile exomeChip.GWAS.HBM.650Y.ld.f.b --remove related.to.remove.txt  --make-bed --out AOGC.pca.c
plink --bfile exomeChip.GWAS.HBM.650Y.ld.f.b --remove related.to.remove.txt  --exclude snps.any.probs.txt --allow-no-sex  --make-bed --out AOGC.pca
./shellfish.py --pca --numpcs 10 --maxprocs 7 --file AOGC.pca --out AOGC.pca

##OOpps forgor HBM samples 7389 snps in common use c("exomeChip.GWAS.650Y.f.1.ld.bim","HBM_final_to_impute_Best.bim")

## plink --bfile HBM_final_to_impute_Best --extract common_snps_HBM.txt  --make-bed --out HBM_Common
## plink --bfile HBM_Common --remove aogc.exome.samples --make-bed --out HBM_Common_Unique
## plink --bfile HBM_Common_Unique --keep to.keep.HBM.txt --make-bed --out HBM_Common_Unique.c

plink --bfile recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL_withRS  --remove related.to.remove.recode.txt --extract common_snps_Chip.GWAS.AOGC.650.HBM.txt  --make-bed --out AOGC_exomeChip_withRS
cut -d " "  -f 1,2   AOGC_exomeChip_withRS.fam > aogc.exome.samples
plink --bfile HBM_final_to_impute_Best --remove aogc.exome.samples --make-bed --out HBM_Best_Unique
plink --bfile AOGC_merge_common._Best --remove aogc.exome.samples --make-bed --out AOGC_merge_common-Unique




plink --bfile ALL_eth_650Y_forward_hg19_Best --extract common_snps_Chip.GWAS.AOGC.650.HBM.txt  --make-bed --out ALL_eth_Best_Common
plink --bfile AOGC_merge_common-Unique --extract common_snps_Chip.GWAS.AOGC.650.HBM.txt  --make-bed --out AOGC_merge_Best_Common
plink --bfile HBM_Best_Unique --extract common_snps_Chip.GWAS.AOGC.650.HBM.txt  --make-bed --out HBM_Best_Unique.noExomeC_Common
#rs7824
plink --bfile HBM_Best_Unique.noExomeC_Common --merge-list files.to.merge.txt --allow-no-sex  --make-bed --out exomeChip.GWAS.HBM.650Y

plink --bfile  exomeChip.GWAS.HBM.650Y  --exclude Price2008_hg19_paul.txt  --range --make-bed --out  exomeChip.GWAS.HBM.650Y.ld  --allow-no-sex
plink --bfile exomeChip.GWAS.HBM.650Y.ld --allow-no-sex --maf 0.05 --max-maf 0.95 --geno 0.01 --hwe 0.00001  --make-bed --out exomeChip.GWAS.HBM.650Y.ld.f
#### 5 samples has high missing rates 
plink --bfile exomeChip.GWAS.HBM.650Y.ld.f --allow-no-sex  --mind 0.02 --make-bed --out  exomeChip.GWAS.HBM.650Y.ld.f.b ## 5 removed here 


plink --bfile exomeChip.GWAS.HBM.650Y.ld --allow-no-sex --maf 0.01 --max-maf 0.99 --geno 0.001 --hwe 0.0000001 
### check those samples in the original data: and they DO NOT have high MISSING RATES- strange they appear to be missing common genotypes ... killing thes 5 samples
## FID	IID	MISS_PHENO	N_MISS	N_GENO 	F_MISS
## hi_s_216	hi_s_216	N	3514	7024	0.5003
## AOGC-03-0039	AOGC-03-0039	N	3492	7024	0.4972
## AOGC-02-0385	AOGC-02-0385	N	1829	7024	0.2604
## lo_s_302	lo_s_302	N	431	7024	0.06136
## AOGC-02-0214	AOGC-02-0214	N	209	7024	0.02976


#### do first pass PCA 

./shellfish.py --pca --numpcs 10 --maxprocs 8 --file exomeChip.GWAS.HBM.650Y.ld.f.b --out exomeChip.GWAS.HBM.650Y.ld.f.b.1
####################### Clean up for PCA 


plink --bfile exomeChip.GWAS.HBM.650Y.ld.f.b --keep exomeChip.GWAS.HBM.650Y.ld.f.b.1.evecs_keep.txt --allow-no-sex --maf 0.05 --max-maf 0.95 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out exomeChip.GWAS.HBM.650Y.ld.f.b.1


./shellfish.py --pca --numpcs 10 --maxprocs 8 --file exomeChip.GWAS.HBM.650Y.ld.f.b.1 --out exomeChip.GWAS.HBM.650Y.ld.f.b.2
######### BAD!

           -9        Adygei        Basque       Bedouin Biaka_Pygmies 
         7592            16            24            11             1 
         Case       Control         Druze        French       Italian 
          781           444            46            29            13 
     Orcadian   Palestinian       Russian     Sardinian        Tuscan 
           16            13            25            28             8
####### keep spikes in place cut at 0.005

plink --bfile exomeChip.GWAS.HBM.650Y.ld.f.b --keep exomeChip.GWAS.HBM.650Y.ld.f.b.1.evecs_wSpikes.keep.txt --allow-no-sex --maf 0.05 --max-maf 0.95 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out exomeChip.GWAS.HBM.650Y.ld.f.b.1s

./shellfish.py --pca --numpcs 10 --maxprocs 8 --file exomeChip.GWAS.HBM.650Y.ld.f.b.1s --out exomeChip.GWAS.HBM.650Y.ld.f.b.2s

exomeChip.GWAS.HBM.650Y.ld.f.b.1.evecs_keep.txt


plink --bfile exomeChip.GWAS.HBM.650Y.ld.f.b --keep exomeChip.GWAS.HBM.650Y.ld.f.b.1.evecs_wSpikes.keep.01.txt --exclude snps.with.probs.txt --allow-no-sex --maf 0.05 --max-maf 0.95 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out exomeChip.GWAS.HBM.650Y.ld.f.b.1s2



./shellfish.py --pca --numpcs 10 --maxprocs 8 --file exomeChip.GWAS.HBM.650Y.ld.f.b.1s2 --out exomeChip.GWAS.HBM.650Y.ld.f.b.2s2
## fails

plink --bfile exomeChip.GWAS.HBM.650Y.ld.f.b.1s2 --exclude snps.with.probs.or.uqdi.or.targets.txt  --make-bed --allow-no-sex --out exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf
./shellfish.py --pca --numpcs 10 --maxprocs 8 --file exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf --out exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf1

./shellfish.py --file exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf --numpcs 10 --evecs exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf1.evecs --snpload --out exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf

load<-read.table("/media/Bioinform-D/Research/shellfish/shellfish/exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.snpload.map",header=F,fill=TRUE,stringsAsFactors=FALSE)
colnames(load)<-c("chr","SNP","cm","POS","a1","a2","maf",paste("PCA",1:10,sep=""))
load[1:5,]

## order.by<-order(abs(load[,"PCA1"]),decreasing=TRUE)
## load[order.by,][1:10,P]


hist(as.numeric(load[,"PCA1"]),breaks=200)

remove<-abs(as.numeric(load[,"PCA1"]))>5
sum(remove)
high.loads<-load[remove,"SNP"]

write.table(high.loads,file="/media/Bioinform-D/Research/shellfish/shellfish/snps.with.loads.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

plink --bfile exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf --exclude snps.with.loads.txt  --make-bed --allow-no-sex --out exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.nl
./shellfish.py --pca --numpcs 10 --maxprocs 8 --file exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.nl --out exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.nl

./shellfish.py --file exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.nl --numpcs 10 --evecs exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.nl.evecs --snpload --out exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.nl


load<-read.table("/media/Bioinform-D/Research/shellfish/shellfish/exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.nl.snpload.map",header=F,fill=TRUE,stringsAsFactors=FALSE)
colnames(load)<-c("chr","SNP","cm","POS","a1","a2","maf",paste("PCA",1:10,sep=""))
load[1:5,]
hist(as.numeric(load[,"PCA2"]),breaks=200)
hist(as.numeric(load[,"PCA2"]),breaks=200,xlim=c(-10,10))
remove<-abs(as.numeric(load[,"PCA2"]))>2.5
sum(remove)
high.loads1<-load[remove,"SNP"]
length(high.loads1)

hist(as.numeric(load[,"PCA3"]),breaks=200)
hist(as.numeric(load[,"PCA3"]),breaks=200,xlim=c(-10,10))
remove<-abs(as.numeric(load[,"PCA3"]))>3
sum(remove)
high.loads2<-load[remove,"SNP"]
length(high.loads2)

high.loads<-unique(c(high.loads1,high.loads2))
length(high.loads)

write.table(high.loads,file="/media/Bioinform-D/Research/shellfish/shellfish/snps.with.loads23.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

plink --bfile exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.nl --exclude snps.with.loads23.txt  --make-bed --allow-no-sex --out exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.nl23
./shellfish.py --pca --numpcs 10 --maxprocs 8 --file exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.nl23 --out exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.nl23

############ yes FUCKING YES!!!!!!!!!!
#/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/PCA_calc/exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.nl23.evecs_pca_eigenvectors.txt
## http://www.sciencedirect.com/science/article/pii/S0002929713002103 loco-ld looks goo too
http://www.cs.tau.ac.il/~heran/cozygene/software.shtml

plink --bfile exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.nl23 --keep exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.nl23.evecs_keep_6SD.txt --make-bed --allow-no-sex --out exomeChip.FINAL
./shellfish.py --pca --numpcs 10 --maxprocs 8 --file  exomeChip.FINAL --out  exomeChip.FINAL

./shellfish.py --file exomeChip.FINAL --numpcs 10 --evecs exomeChip.FINAL.evecs --snpload --out exomeChip.FINAL

plink --file data --indep 50 5 1.5 --out ldtrimset

load<-read.table("/media/Bioinform-D/Research/shellfish/shellfish/exomeChip.FINAL.snpload.map",header=F,fill=TRUE,stringsAsFactors=FALSE)
colnames(load)<-c("chr","SNP","cm","POS","a1","a2","maf",paste("PCA",1:10,sep=""))
load[1:5,]
hist(as.numeric(load[,"PCA1"]),breaks=200)
hist(as.numeric(load[,"PCA1"]),breaks=200,xlim=c(-10,10))
remove<-abs(as.numeric(load[,"PCA1"]))>3
sum(remove)
high.loads<-load[remove,"SNP"]
length(high.loads)

write.table(high.loads,file="/media/Bioinform-D/Research/shellfish/shellfish/snps.with.loads4.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

plink --bfile  exomeChip.FINAL --exclude snps.with.loads4.txt  --make-bed --allow-no-sex --out exomeChip.FINAL.1
./shellfish.py --pca --numpcs 10 --maxprocs 8 --file  exomeChip.FINAL.1 --out  exomeChip.FINAL.1

./shellfish.py --file exomeChip.FINAL.1 --numpcs 10 --evecs exomeChip.FINAL.1.evecs --snpload --out exomeChip.FINAL.1

load<-read.table("/media/Bioinform-D/Research/shellfish/shellfish/exomeChip.FINAL.1.snpload.map",header=F,fill=TRUE,stringsAsFactors=FALSE)
colnames(load)<-c("chr","SNP","cm","POS","a1","a2","maf",paste("PCA",1:10,sep=""))
load[1:5,]
hist(as.numeric(load[,"PCA2"]),breaks=200)
hist(as.numeric(load[,"PCA2"]),breaks=200,xlim=c(-10,10))
remove<-abs(as.numeric(load[,"PCA2"]))>5
sum(remove)
high.loads1<-load[remove,"SNP"]
length(high.loads1)

hist(as.numeric(load[,"PCA1"]),breaks=200)
hist(as.numeric(load[,"PCA1"]),breaks=200,xlim=c(-10,10))
remove<-abs(as.numeric(load[,"PCA1"]))>5
sum(remove)
high.loads2<-load[remove,"SNP"]
length(high.loads2)

high.loads<-unique(c(high.loads1,high.loads2))
length(high.loads)

write.table(high.loads,file="/media/Bioinform-D/Research/shellfish/shellfish/snps.with.loads5.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

plink --bfile  exomeChip.FINAL.1 --exclude snps.with.loads5.txt  --make-bed --allow-no-sex --out exomeChip.FINAL.2

plink --bfile  exomeChip.FINAL.2 --exclude ldtrimset.prune.out  --make-bed --allow-no-sex --out exomeChip.FINAL.2.p

./shellfish.py --pca --numpcs 10 --maxprocs 8 --file  exomeChip.FINAL.2.p --out  exomeChip.FINAL.2.p

############### best solution founds: exomeChip.FINAL.2.p
./shellfish.py --pca --numpcs 10 --maxprocs 8 --file exomesChip.rs.f.ld.c --out exomesChip.rs.f.ld.c.1

cut -f 2 exomeChip.FINAL.2.p.bim > working.strat.snps.txt
########### test thes on all ethnicityes

plink --bfile exomeChip.GWAS.HBM.650Y.ld.f.b --extract working.strat.snps.txt --make-bed --out exomeChip.GWAS.HBM.650Y.strat.FINAL --allow-no-sex

/media/Bioinform-D/Research/shellfish/shellfish/


./shellfish.py --pca --numpcs 10 --maxprocs 8 --file exomeChip.GWAS.HBM.650Y.strat.FINAL --out exomeChip.GWAS.HBM.650Y.strat.FINAL

################## so these snps work CHECK ARE SELF-CONSISTANT :

plink --bfile exomeChip.GWAS.HBM.650Y.strat.FINAL --keep exomeChip.GWAS.HBM.650Y.strat.FINAL.evecs_keep_6SD.txt --make-bed --out exomeChip.FINAL.BEST

./shellfish.py --pca --numpcs 10 --maxprocs 2 --file exomeChip.FINAL.BEST --out exomeChip.FINAL.BEST






















## plink --bfile exomeChip.FINAL.1 --indep 50 5 1.5 --out ldtrimset

## plink --bfile exomeChip.FINAL.1  --exclude ldtrimset.prune.out --make.bed --out exomeChip.FINAL.1.tr

## exomeChip.GWAS.HBM.650Y.ld.f.b.1.evecs_wSpikes.01.keep.txt
## ############### major problem need to resolve  use tricks  --out exomeChip.GWAS.HBM.650Y.ld.f.b.2

## plink --bfile recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL_withRS  --remove related.to.remove.recode.txt  --maf 0.05 --max-maf 0.95 --geno 0.05 --hwe 0.0000001  --make-bed --out exomesChip.rs

## plink --bfile  exomesChip.rs --allow-no-sex  --mind 0.02 --make-bed --out  exomesChip.rs.f

## plink --bfile  exomesChip.rs.f  --exclude Price2008_hg19_paul.txt  --range --make-bed --out exomesChip.rs.f.ld  --allow-no-sex

## plink --bfile exomesChip.rs.f.ld --exclude snps.with.probs.or.uqdi.or.targets.txt --allow-no-sex   --make-bed --out exomesChip.rs.f.ld.c


## ############### best solution founds:
## ./shellfish.py --pca --numpcs 10 --maxprocs 8 --file exomesChip.rs.f.ld.c --out exomesChip.rs.f.ld.c.1

## cut -f 2 exomesChip.rs.f.ld.c.bim > working.strat.snps.txt
## ########### test thes on all ethnicityes

## plink --bfile exomeChip.GWAS.HBM.650Y.ld.f.b --keep working.strat.snps.txt --make-bed --out exomeChip.GWAS.HBM.650Y.strat.FINAL --allow-no-sex




## all.snps<-{}

## setwd("/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/PCA_calc")
## ped.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1.fam"

## ped<-read.table(ped.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
## hbm.ped<-read.table("/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/PCA_calc/HBM_Best_Unique.noExomeC_Common.fam",header=F,fill=TRUE,stringsAsFactors=FALSE)
## aogc.ped<-read.table("/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/PCA_calc/AOGC_merge_Best_Common.fam",header=F,fill=TRUE,stringsAsFactors=FALSE)
## echip.ped<-read.table("/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/PCA_calc/AOGC_exomeChip_withRS.fam",header=F,fill=TRUE,stringsAsFactors=FALSE)


## #################choose 1
## cases<-ped[,1] %in% hbm.ped[,1]

## cases<-ped[,1] %in% aogc.ped[,1]

## cases<-ped[,1] %in% echip.ped[,1]
## ########################

## controls<-!cases
## sum(cases)
## sum(controls)

## ped[cases,6]<-2
## ped[controls,6]<-1
## ped[1:5,]

## write.table(ped,file="HBM.test.fam",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

## write.table(ped,file="aogc.test.fam",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

## write.table(ped,file="echip.test.fam",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)


## #################choose 1
## plink --bed exomeChip.GWAS.HBM.650Y.ld.f.b.1.bed --bim exomeChip.GWAS.HBM.650Y.ld.f.b.1.bim --fam HBM.test.fam --allow-no-sex --assoc

## plink --bed exomeChip.GWAS.HBM.650Y.ld.f.b.1.bed --bim exomeChip.GWAS.HBM.650Y.ld.f.b.1.bim --fam aogc.test.fam --allow-no-sex --assoc

## plink --bed exomeChip.GWAS.HBM.650Y.ld.f.b.1.bed --bim exomeChip.GWAS.HBM.650Y.ld.f.b.1.bim --fam echip.test.fam --allow-no-sex --assoc
## ##########################


## assoc<-read.table("plink.assoc",header=T,fill=TRUE,stringsAsFactors=FALSE)
## assoc[1:5,]
## hist(-log10(assoc[,"P"])) # lo
## remove<-assoc[,"P"] < 0.01
## sum(remove)
## snps<-assoc[remove,"SNP"]
## all.snps<-c(all.snps,snps)
## all.snps

## all.snps<-unique(all.snps)
## length(all.snps)

## write.table(all.snps,file="snps.with.chip.bias.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

## plink --bfile exomeChip.GWAS.HBM.650Y.ld.f.b --exclude snps.with.chip.bias.txt  --allow-no-sex --maf 0.05 --max-maf 0.95 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB

## ./shellfish.py --pca --numpcs 10 --maxprocs 8 --file exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB --out exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB

## plink --bfile exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB --keep exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.evecs_keep.txt --maf 0.05 --max-maf 0.95 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out  exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2 




## ./shellfish.py --pca --numpcs 10 --maxprocs 8 --file exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2 --out exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.3


## ########### phenotypes??
## #pheno.file<-"/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_SAMPLES_PHENOTYPES.txt" # in /media/Bioinform-D/Research/GWAS extreme regression/AOGC_HBM/
## setwd("/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/PCA_calc")
## ped.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1.fam"

## ped<-read.table(ped.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
## pheno.file<-"/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_SAMPLES_PHENOTYPES_Nov.1.2013.txt"
## pheno<-read.delim(pheno.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
## pheno[1:5,]
## posns<-match(ped[,1],pheno[,"PATIENT"])
## missing<-is.na(posns)
## ped[missing,]

## write.table(ped[missing,],file="Genotyped.but.no.phenotypes.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## no.pheno<- ped[missing,1]

## pheno.have<-pheno[posns,]
## dim(pheno.have)
## pheno.have[40:50,]

## aff<-pheno.have[,"TOT_HIP_Z"]
## aff[1:50]
## table(aff)
## miss<-pheno.have[is.na(aff), "PATIENT"]
## miss<-unique(miss,no.pheno)

## ped[,6]<-aff
## write.table(cbind(miss,miss),file="exclude.samples.from.pheno.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(ped,file="pheno.fam",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)


## plink --bed exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.bed --bim exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.bim --fam pheno.fam --exclude exclude.samples.from.pheno.txt  --allow-no-sex --assoc

## ### STILL NOT WORKING
## #####################################
## aff<-pheno.have[,"BMD_AFFSTAT"]
## aff[1:50]
## length(aff)

## aff[is.na(aff)]<--9
## table(aff)

## ped[,6]<-aff
## #write.table(cbind(miss,miss),file="exclude.samples.from.pheno.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(ped,file="pheno_BMD_AFF.fam",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
## ########################################

## ####################################
## aff<-pheno.have[,"HEIGHT"]
## aff[1:50]
## length(aff)

## aff[is.na(aff)]<--9
## table(aff)

## ped[,6]<-aff
## #write.table(cbind(miss,miss),file="exclude.samples.from.pheno.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(ped,file="pheno_HEIGHT.fam",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
## ########################################

## ####################################
## aff<-pheno.have[,"WEIGHT"]
## aff[1:50]
## length(aff)

## aff[is.na(aff)]<--9
## table(aff)

## ped[,6]<-aff
## #write.table(cbind(miss,miss),file="exclude.samples.from.pheno.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(ped,file="pheno_WEIGHT.fam",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
## ########################################

## plink --bed exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.nPB.bed --bim exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.nPB.bim --fam pheno_BMD_AFF.fam --allow-no-sex --assoc

## plink --bed exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.nPB.bed --bim exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.nPB.bim  --fam pheno_HEIGHT.fam --allow-no-sex --assoc

## plink --bed exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.nPB.bed --bim exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.nPB.bim  --fam pheno_WEIGHT.fam --allow-no-sex --assoc


## assoc<-read.table("/media/Bioinform-D/Research/shellfish/shellfish/plink.qassoc",header=T,fill=TRUE,stringsAsFactors=FALSE)

## assoc<-read.table("plink.assoc",header=T,fill=TRUE,stringsAsFactors=FALSE)
## assoc<-read.table("plink.qassoc",header=T,fill=TRUE,stringsAsFactors=FALSE)

## assoc[1:5,]
## hist(-log10(assoc[,"P"])) # lo
## remove<-assoc[,"P"] < 0.01
## sum(remove)
## snps<-assoc[remove,"SNP"]
## all.snps<-c(snps)
## all.snps
## assoc[order(assoc[,"P"]),][1:50,]

## all.snps<-unique(all.snps)
## length(all.snps)

## write.table(all.snps,file="snps.with.pheno.BMD_AFF.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(all.snps,file="snps.with.pheno.HEIGHT.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(all.snps,file="snps.with.pheno.WEIGHT.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

## all.snps<-read.table("snps.all.pheno.txt",header=F,fill=TRUE,stringsAsFactors=FALSE)
## all.snps<-read.table("snps.with.probs.txt",header=F,fill=TRUE,stringsAsFactors=FALSE)
## all.snps<-as.vector(all.snps[,1])
## all.snps<-unique(all.snps)
## length(all.snps)
## write.table(all.snps,file="snps.any.probs.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)


## plink --bfile exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2  --exclude snps.with.pheno.bias.txt --allow-no-sex --maf 0.05 --max-maf 0.95 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.nPB

## ./shellfish.py --pca --numpcs 10 --maxprocs 8 --file exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.nPB --out exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.nPB.1

## plink --bfile exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.nPB --remove all.HBM.txt --make-bed --out sanity
## ./shellfish.py --pca --numpcs 10 --maxprocs 8 --file sanity --out sanity.1
## ############## test to see if the set of SNP can do population stratification:

## plink --bfile exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.nPB --exclude snps.any.pheno.txt  --allow-no-sex --maf 0.05 --max-maf 0.95 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.nPB2

## ./shellfish.py --pca --numpcs 10 --maxprocs 8 --file exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.nPB2 --out exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.nPB.2

## plink --bfile exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.nPB2  --exclude Price2008_hg19_paul.txt  --range --make-bed  --allow-no-sex

## plink --bfile exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB --exclude snps.with.pheno.bias.txt --allow-no-sex --maf 0.05 --max-maf 0.95 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.nPB

## ./shellfish.py --pca --numpcs 10 --maxprocs 8 --file exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.nPB --out exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.nPB



## ## plink --bfile exomeChip.GWAS.650Y.f.1.ld --bmerge  HBM_Common_Unique.bed   HBM_Common_Unique.bim   HBM_Common_Unique.fam --make-bed --out exomeChip.GWAS.HBM.650Y.f.1.ld

## ## plink --bfile exomeChip.GWAS.HBM.650Y.f.1.ld --allow-no-sex --maf 0.05 --max-maf 0.95 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out exomeChip.GWAS.HBM.650Y.f.1.ld.c

## ## ./shellfish.py --pca --numpcs 10 --maxprocs 7 --file exomeChip.GWAS.HBM.650Y.f.1.ld.c --out exomeChip.GWAS.HBM.650Y.f.2.ld.c



##  HBM_Common
## /media/Bioinform-D/Research/shellfish/shellfish/exomeChip.GWAS.650Y.f.1.ld.bed
## /media/Bioinform-D/Research/shellfish/shellfish/exomeChip.GWAS.650Y.f.1.ld.fam
## /media/Bioinform-D/Research/shellfish/shellfish/exomeChip.GWAS.650Y.f.1.ld.bim

























dim(all.data)
key<-paste(all.data[,"chr"],all.data[,"POS"],sep=":")
same<-tapply(key,key,length)
same[1:5]
sum(same>1)# 900
dups<-unique(names(same)[same>1])
length(dups)
dups[1:5]

duplicates<-key %in% dups
sum(duplicates) #1802

the.dups<-all.data[duplicates,] #900 but some are REAL duplications

key<-paste(the.dups[,"chr"],the.dups[,"POS"],the.dups[,"A1.forward"],the.dups[,"A2.forward"],sep=":")
same<-tapply(key,key,length)


same[1:5]
sum(same>1) #789
dups<-unique(names(same)[same>1])
length(dups)
dups[1:5]

the.dups[1:5,]

true.duplicates<-key %in% dups
sum(true.duplicates) #1600

illumina.doubles<-the.dups[true.duplicates,]
illumina.doubles[1:6,]

write.table(illumina.doubles,file="illumina.replicate.snps.on.exomeChipv1.1.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

poly.snps<-the.dups[!true.duplicates,]
poly.snps[1:6,]
dim(poly.snps)
write.table(poly.snps,file="polymorphic.snps.on.exomeChipv1.1.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


poly.morphic.snps<-poly.snps[,"SNP"]
save(list=c("poly.morphic.snps","all.data"),file="polymorphic.snps.RData")
getwd()


real.duplicates<-

missing<-is.na(posns)
sum(!missing)

load("/media/Bioinform-D/Research/AML sequencing/forward.remain.RData")
forward<-forward.remain
forward[1:5,]
posns<-match(forward[,"SNP"],all.data[,"SNP"])
sum(is.na(posns))

test<-cbind(forward,all.data[posns,])

test[1:5,]
test<-"exm-rs3890745"

loc<-grep(test,all.data[,"SNP"])
loc

cbind(all.data,to.flip)[loc:(loc+1),]


cbind(all.data,to.flip)[(loc-2):(loc+2),]




ls()ls()



## beds<-beds[core.lists]
## scores
## snp.lists

## phyper(4, 416, (48000-416), 485, lower.tail = TRUE, log.p = FALSE) 

## choose(6,0)*choose(43,6)/choose(49,6)
## choose(6,6)*choose(43,0)/choose(49,6)

## exp(lchoose(6,6)+lchoose(43,0)-lchoose(49,6))

## choose(485,4)*choose((48000-485),(485-4))/choose(48000,485)
## (lchoose(485,4)+lchoose((48000-485),(485-4)))-lchoose(48000,485)

## exp(lchoose(485,4)+lchoose((48000-485),(485-4))-lchoose(48000,485))

## exp(lchoose(485,1)+lchoose((48000-485),(485-1))-lchoose(48000,485))
## ls()

## n<-7
## exp(lchoose(485,n)+lchoose((48000-485),(485-n))-lchoose(48000,485))



#################### add to final annotations 



code.dir<-"/media/Bioinform-D/Research/AML sequencing"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")

ann<-read.delim("/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/zCall_AOGC.with.evoker_corrected_clean_FINAL.analysis.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
load("/media/UQCCG-Analysis/AOGC_exome_chip/SNP_lists_used_in_exomeChip/uqdi.design.RData")
common.all<-cbind(common.all,common[,c("F_A","F_U","P","OR","origin")],stringsAsFactors=FALSE)

core.ann<-c("chr","start")
the.key<-build.key(ann,core.ann)
design<-rep(NA,times=length(the.key))
data.names<-c("hla","pop.strat2","pop.strat","OA.gws","gfos.gws","gfos.meta","height.gws","rare","rare2","common.all")

i<-1

for(i in 1:length(data.names)){
data<-eval(as.name(data.names[i]))
print("------------------------------------------")
print(data.names[i])
print(dim(data))
print(data[1,c(core.ann,"origin")])
print("------------------------------------------")
key<-build.key(data,core.ann)
posns<-match(the.key,key)
missing<-is.na(posns)
print(sum(!missing))

## data[posns[!missing],c(core.ann,"origin")][1:5,]
## ann[!missing,][1:5,1:10]

design[!missing]<-paste(design[!missing],data[posns[!missing],"origin"],sep=";")


#assign(data.names[i],value=data)
}
sum(!is.na(design))
sum(grepl("UQDI",ann[,"SNP"]))

missing<-is.na(design)
done<-grepl("^UQDI",ann[,"SNP"])

sum(missing & done)
design[missing & done][1:10]
ann[missing & done,][1:5,]


design<-gsub("^\\;+","",design)
design<-gsub("NA;","",design)

tapply(design,design,length)
> tapply(design,design,length)
##               common             GFOS.gws      GFOS.gws;common 
##                 2704                  171                    1 
##   GFOS.gws;GFOS.meta            GFOS.meta     GFOS.meta;common 
##                   14                 4096                    2 
##           Height.gws                   HL         HL;GFOS.meta 
##                  150                   79                    1 
##               OA.gws            pop.strat           pop.strat2 
##                   22                   15                    4 
## pop.strat2;pop.strat                 rare                rare2 
##                   25                13917                 1961 
##          rare;common 
##                  102

ann[1:5,1:15]
ann<-cbind(ann[,1:5],design,ann[,6:180])
getwd()

setwd("/media/UQCCG-Analysis/AOGC_exome_chip/Analysis")
write.table(ann,file="annoatation_extensive_exome_chip.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
#ann<-read.delim("/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/annoatation_extensive_exome_chip.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
colnames(ann)

ann.short<-ann[,c(1:12,22,23,27:31,34,39,40:54),]
write.table(ann.short,file="annoatation_short_exome_chip.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

#ann.short<-read.delim("/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/annoatation_short_exome_chip.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

logistic.files<-c(
                  "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/Hip/totalhipBMDassocinexcel.txt",
                  "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/LS/LSBMDassocinexcel.txt",
                  "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/FN/FemNeckBMDassocinexcel.txt",
                  "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/AllFr_50_OP/AllFr_50_OP.logistic.covar.assoc.logistic",
                  "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/ForearmFrLoTrauma/ForearmFrLoTrauma.logistic.covar.assoc.logistic",
                  "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/everFx18nottrivia/everFx18nottrivia.LOGISTIC.COVAR.assoc.logistic",
                  "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/EverFx/everFx_final_analysis.logistic.covar.assoc.logistic",
                  "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/AllVertFx_OP/AllVertFX_OP.logistic.covar.assoc.logistic",
                  "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/AllNonVertFX_50_OP/AllNonVertFX_50_OP.logistic.covar.assoc.logistic",
                  "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/EMcC_BP_GRP/EMcC_BP_GRP.logistic.covar.assoc.logistic",
                  "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/AllHipFr_50_OP/AllHipFr_50_OP.logistic.covar.assoc.logistic"
                  )


## ann.short[ann.short[,"design"]!="",design]
## tapply(ann.short[,"design"],ann.short[,"design"],length)
## for(i in 1:dim(ann.short)[2]){
  
## ann.short[,i]<-gsub(":",".",ann.short[,i])
## ann.short[,i]<-gsub(" ",".",ann.short[,i])
## ann.short[,i]<-gsub("(",".",ann.short[,i],fixed=TRUE)
## ann.short[,i]<-gsub(")",".",ann.short[,i],fixed=TRUE)
## #ann.short[,i]<-gsub("",".",ann.short[,i])
## ann.short[is.na(ann.short[,i]),i]<-"."
## }


colnames(ann.short)[colnames(ann.short)=="chr"]<-"chr.ann"
## ann.short[is.na(ann.short[,"design"]),"design"]<-NA
#gsub("",".",test[300:303,"hgnc_symbol"])

## ann
#setwd("/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/AllFr_50_OP")

wanted.cols<-c("design","refGene..location","refGene..type","refGene..gene","Gene.Names","description","skeletome","mouse.defect","sewell.cycling","Dequeant.cycling","ingenuity.bone.genes","Consequence.Embl","Amino_acids.Embl","wanted.muts","wanted.muts.coding","MAF.lt.0","MAF.lt.0.001","MAF.lt.0.005","MAF.lt.0.01","MAF.lt.0.025","MAF.lt.0.5","gerp.scores","PolyPhen.desc","PolyPhen.scores","SIFT.desc","SIFT.scores","ensembl_gene_id")

#wanted<-c("CHR","SNP","POS","A1","TEST","NMISS","OR","SE","L95","U95","STAT","P","design")

cols.to.halpview.fix<-c("design","refGene..location","refGene..type","refGene..gene","Gene.Names","description","skeletome","mouse.defect","sewell.cycling","Dequeant.cycling","ingenuity.bone.genes","Consequence.Embl","Amino_acids.Embl","wanted.muts","wanted.muts.coding","MAF.lt.0","MAF.lt.0.001","MAF.lt.0.005","MAF.lt.0.01","MAF.lt.0.025","MAF.lt.0.5","gerp.scores","PolyPhen.desc","PolyPhen.scores","SIFT.desc","SIFT.scores","ensembl_gene_id")



#wanted<-c("CHR","SNP","POS","A1","TEST","NMISS","OR","SE","L95","U95","STAT","P","design","refGene..location","refGene..type","refGene..gene","Gene.Names","skeletome","mouse.defect","sewell.cycling","Dequeant.cycling","ingenuity.bone.genes","Consequence.Embl","Amino_acids.Embl","wanted.muts","wanted.muts.coding","MAF.lt.0","MAF.lt.0.001","MAF.lt.0.005","MAF.lt.0.01","MAF.lt.0.025","MAF.lt.0.5","gerp.scores","PolyPhen.desc","PolyPhen.scores","SIFT.desc","SIFT.scores","ensembl_gene_id","hgnc_symbol")

i<-2
#for(i in 1:length(logistic.files)){

  for(i in 1:length(logistic.files)){
print(logistic.files[i])
logistic<-read.table(logistic.files[i],header=T,fill=TRUE,stringsAsFactors=FALSE)

colnames(logistic)[colnames(logistic)=="BP"]<-"POS"

if("TEST" %in% colnames(logistic)){
wanted<-grepl("ADD",logistic[,"TEST"])
logistic<-logistic[wanted,]
}

posns<-match(logistic[,"SNP"],ann.short[,"SNP"])
test<-cbind(logistic,ann.short[posns,]) # ,stringsAsFactors=FALSE

dim(test)

#ifh<-1
for(ifh in 1:length(cols.to.halpview.fix)){
  
test[,cols.to.halpview.fix[ifh]]<-gsub(":",".",test[,cols.to.halpview.fix[ifh]])
test[,cols.to.halpview.fix[ifh]]<-gsub(" ",".",test[,cols.to.halpview.fix[ifh]])
test[,cols.to.halpview.fix[ifh]]<-gsub("(",".",test[,cols.to.halpview.fix[ifh]],fixed=TRUE)
test[,cols.to.halpview.fix[ifh]]<-gsub(")",".",test[,cols.to.halpview.fix[ifh]],fixed=TRUE)
#test[,i]<-gsub("",".",test[,i])
test[is.na(test[,cols.to.halpview.fix[ifh]]),cols.to.halpview.fix[ifh]]<-"."
}


write.table(test[,c(colnames(logistic),wanted.cols)],file=paste(logistic.files[i],"ann",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

}


######################################################### single point for exome chip
######################################################### single point for exome chip
######################################################### single point for exome chip
#fam.template.file<-"recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL"

fam.template.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.fam"
#fam.template.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL_chr1.fam"
#fam.template.file<-"/media/UQCCG/UQCCG-Projects/PAUL_LEO/AOGC exome chip core genotyping/keep.txt"
annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"

#recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL
#the.chr<-basename(fam.template.file)
## the.fam.dir<-dirname(fam.template.file)
## files<-dir(the.fam.dir)  
## projects<-gsub(".bed","",files[grepl(".bed$",files)])  # files[grepl(".bed$",files)]
setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/")
projects<-"recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL"

fam<-read.table(fam.template.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

related.to.remove<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.recode.txt",header=F,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
related.to.remove<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.txt",header=F,fill=TRUE,sep=",",stringsAsFactors=FALSE)
related.to.remove[1:5,]

sum(related.to.remove[,1] %in% fam[,1])


traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")

#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
names(traits)<-rep("logistic",times=length(traits))
names(traits)[1:3]<-rep("assoc",times=3)
traits
traits %in% colnames(ann)


fam[1:5,]

i<-1
the.fam.dir<-dirname(fam.template.file)
setwd(the.fam.dir)


for (i in 4:length(traits)){

  
fam[,6]<--9
posns<-match(fam[,1],ann[,"PATIENT"])
missing<-is.na(posns)
print(sum(missing))

no.pheno<-ann[posns,traits[i]]==3 | ann[posns,traits[i]]==-9 | is.na(ann[posns,traits[i]])  # ann[posns,traits[i]]

print(sum(no.pheno))


missing<-missing | no.pheno
#cbind(fam[no.pheno,],ann[posns[no.pheno],c("PATIENT",traits[i])])
print(sum(missing))

if(names(traits)[i]=="logistic"){
 using.covars<-c("AGE_SCAN","PCA1","PCA2","PCA3","PCA4")
  covars<- ann[posns,c("PATIENT","PATIENT", using.covars)]
 colnames(covars)[1:2]<-c("FID","IID")
 covars[1:5,]
  missing.covars<-is.na(covars[,using.covars])
   missing.covars<-apply(missing.covars,1,sum,na.rm=TRUE)
 missing.covars<-  missing.covars>0
missing<-missing | missing.covars
write.table(covars[!missing,],file=paste("covars",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
}


has.missing<-sum(missing)>0
has.missing
write.table(fam[missing,c(1,1)],file=paste("remove_from",traits[i],sep="."),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

fam[!missing,6]<-ann[posns[!missing],traits[i]]
the.fam<-paste(traits[i],"fam",sep=".")

write.table(fam,file=the.fam,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)



ip<-1
projects
for(ip in 1:length(projects)){
print(projects[ip])
#setwd(projects[ip])

the.bed<-paste(projects[ip],"bed",sep=".")
the.bim<-paste(projects[ip],"bim",sep=".")

#  system( paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam,"--remove",paste("remove_from",traits[i],sep="."),"--hardy","--allow-no-sex","--out",paste(traits[i],projects[ip],sep="."),"--noweb",sep=" ") )


if(names(traits)[i]=="assoc"){
  system( paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam,"--remove",paste("remove_from",traits[i],sep="."),"--assoc","--allow-no-sex","--out",paste(traits[i],projects[ip],sep="."),"--noweb",sep=" ") )
}

if(names(traits)[i]=="logistic"){
  system( paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam,"--remove",paste("remove_from",traits[i],sep="."),"--covar",paste("covars",traits[i],sep="."),"--covar-name AGE_SCAN,PCA1,PCA2,PCA3,PCA4","--logistic --ci 0.95 ","--allow-no-sex","--out",paste(traits[i],projects[ip],sep="."),"--noweb",sep=" ") )
}
} # ip loop over chromosomes
} # i loop over traits


######################################################################################
######################################################################################
######################################################################################
######################################################################################
############################ AOGC VCF conversion

options(width=250,max.print=5000)

code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")


annotation.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)


#traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
## traits<-c("TOT_HIP_GCM","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")

## traits<-c("BMD_EFF_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")


traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")

names(traits)<-rep("logistic",times=length(traits))
names(traits)[1:3]<-rep("assoc",times=3)
traits
traits %in% colnames(ann)




#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")


input.dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2012-06-18_All_AOGC/Analysis"
files<-dir(input.dir)  
project.files<-files[ grepl("^AOGC-Genotyping.output.",files) & grepl(".AOGC_ALL.analysis.txt$",files)]
project.files
#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2012-06-18_All_AOGC/Analysis/AOGC-Genotyping.output.chr8.AOGC_ALL.analysis.txt

ichr<-4
i<-5 # traits

library(doMC)
num.cores<-7
registerDoMC(cores=num.cores)


for(ichr in 2:length(project.files)){
################## fast read ###########
  setwd(input.dir)
column.labels<-read.delim(project.files[ichr],header=F,nrows=1,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
num.vars<-dim(column.labels)[2]
a.indel<-scan(project.files[ichr],what=character(num.vars),skip=1,sep="\t",fill=TRUE,na.strings="",quote="\"")
num.lines<-length(a.indel)/(num.vars)
dim(a.indel)<-c(num.vars,num.lines)
a.indel<-t(a.indel)
colnames(a.indel)<-column.labels
########################################
the.chr<-a.indel[1,"chr"]
print(paste("Doing Chromosome ",the.chr))

if(!grepl("^chr",the.chr)){
a.indel[,"chr"]<-paste("chr",a.indel[,"chr"],sep="")
}
core.ann<-c("chr","start","end","REF","ALT","TYPE")
key<-build.key(a.indel,core.ann)

rownames(a.indel)<-key
plink.chr<-the.chr


if(plink.chr=="X"){plink.chr<-23}
if(plink.chr=="Y"){plink.chr<-24}
if(plink.chr=="XY"){plink.chr<-25}
if(plink.chr=="M"){plink.chr<-26}
plink.chr

the.samples<-colnames(a.indel)[grepl(".GT$",colnames(a.indel))]
length(the.samples)
the.samples.fam<-gsub(".GT$","",the.samples)


####
a.fam<-cbind(the.samples.fam,the.samples.fam,0,0,9,-9)
colnames(a.fam)<-c("PATIENT","SAMPLE","mother","father","sex","pheno")

setwd("/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink")
the.fam<-paste(traits[i],"fam",sep=".")
a.fam<-read.table(the.fam,header=F,fill=TRUE,stringsAsFactors=FALSE)
colnames(a.fam)<-c("PATIENT","SAMPLE","mother","father","sex","pheno")
### colnames(fam)<-c("PATIENT","SAMPLE","mother","father","sex","pheno")


if(dim(a.indel)[1]<200){
num.bits<-1; print(paste("Only",dim(a.indel)[1],"genotypes remaining",sep=" ")) }else{num.bits<-num.cores}

#chunksize=as.integer(dim(a.indel)[1]/num.bits) )

fil.genotypes<-foreach(a.indel.bit=iter(a.indel,by='row',chunksize=5000 ), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype(a.indel.bit,the.samples.fam,prefix="",suffix="",20,0.2,0.80,0.05,0.95,10,5)

dim(fil.genotypes)

a.indel<-cbind(a.indel[,c("chr","start","end","REF","ALT","TYPE")],fil.genotypes)

root<-paste("AOGC_to_vcf_filtered",plink.chr,sep="_")
write.plink(a.indel,root)

} # short_cut loop to write filtered genotypes





################## below used for testing vcf vs plink
################## below used for testing vcf vs plink
################## below used for testing vcf vs plink
################## below used for testing vcf vs plink
################## below used for testing vcf vs plink
################## below used for testing vcf vs plink

use.covar<-read.table("covars.EVER_FX",header=T)


use.fam<-a.fam
table(use.fam[,6])


ann[1:5,]
length(the.samples.fam)
c("PCA1","PCA2","PCA3","PCA4")
posns<-match(a.fam[,"PATIENT"],ann[,"PATIENT"])
length(posns)
missing<-is.na(posns)
sum(missing)
a.ann<-ann[posns,]
dim(a.ann)
a.fam[,6]<-a.ann[,traits[i]]
a.fam[a.fam[,6]==3,6]<--9
a.fam[is.na(a.ann[,"PCA1"]),6]<--9
a.fam[,6]

table(use.fam[,6])
table(a.fam[,6])

sum(use.fam[,1]!=a.fam[,1])
sum(use.fam[,6]!=a.fam[,6])

a.covar<-a.ann[,colnames(use.covar)]
a.covar[1:5,]
sum(use.covar[,6]!=a.covar[,6])
## for EVER_FX fam file is ok



"/media/U5CG-nalysis/AOGC_exome_chip/AOGC_vcf_to_plink/AOGC_to_vcf_filtered_chr13"
a.indel<-read.plink.raw("/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink/AOGC_to_vcf_filtered_13") #A1->ALT POS->start ALT is 1 REF 0

a.indel<-read.plink.raw("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2012-06-18_All_AOGC/PassSNPsConvertedToPlink/AOGC-Genotyping.PASS.ValidSamples.chr5")


dim(a.indel)
a.indel[1:15,1:10]
rownames(g.indel)<-build.key(g.indel,c("chr","start","REF","ALT"))

## the.bed<-paste(projects[ip],"bed",sep=".")
## the.bim<-paste(projects[ip],"bim",sep=".")
the.samples<-colnames(a.indel)[grepl(".GT$",colnames(a.indel))]
length(the.samples)
the.samples.fam<-gsub(".GT$","",the.samples)


####
a.fam<-cbind(the.samples.fam,the.samples.fam,0,0,9,-9)
colnames(a.fam)<-c("PATIENT","SAMPLE","mother","father","sex","pheno")
##   colnames(test) <- sample.geno.order
##   rownames(test) <- indels[,"SNP"]
num.samples<-dim(a.fam)[1]
#
# res <- logscan(genotypes, a.fam,"pheno",num.samples)


i<-1
for (i in 1:length(traits)){


fam<-a.fam 
fam[,6]<--9
posns<-match(fam[,1],ann[,"PATIENT"])
missing<-is.na(posns)
print(sum(missing))

no.pheno<-ann[posns,traits[i]]==3 | ann[posns,traits[i]]==-9 | is.na(ann[posns,traits[i]]) | ann[posns,"GENDER"]!=2   # ann[posns,traits[i]]

print(sum(no.pheno))


missing<-missing | no.pheno
cbind(fam[no.pheno,],ann[posns[no.pheno],c("PATIENT",traits[i])])
print(sum(missing))

if(names(traits)[i]=="logistic"){
 using.covars<-c("AGE_SCAN","PCA1","PCA2","PCA3","PCA4") #  using.covars<-c("CENTRE","WEIGHT","AGE_SCAN","PCA1","PCA2","PCA3","PCA4")
  covars<- ann[posns,c("PATIENT","PATIENT", using.covars)]
 colnames(covars)[1:2]<-c("FID","IID")
 covars[1:5,]
  missing.covars<-is.na(covars[,using.covars])
   missing.covars<-apply(missing.covars,1,sum,na.rm=TRUE)
 missing.covars<-  missing.covars>0
missing<-missing | missing.covars
#write.table(covars[!missing,],file=paste("covars",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
}


has.missing<-sum(missing)>0
has.missing
#write.table(fam[missing,c(1,1)],file=paste("remove_from",traits[i],sep="."),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

fam[!missing,6]<-ann[posns[!missing],traits[i]]
the.fam<-paste(traits[i],"fam",sep=".")

fam[!missing,][1:15,]
ann[posns[!missing],c("PATIENT",traits[i])][1:5,]

fam<-fam[!missing,]
dim(fam)

fam[1:5,]
the.samples.use<-paste(fam[,"PATIENT"],".GT",sep="")
genotypes<-a.indel[,the.samples.use]

genotypes[1:3,1:5]

print("start genotype conversion QC")

grep("14358303",g.indel[,"start"])


genotypes<-filtered.genotype(a.indel[19617:19620,],fam[,1],prefix="",suffix="",30,0.2,0.80,0.05,0.95,15,7)


genotypes<-fil.genotypes[19617:19620,the.samples.use]


a.indel[1:6,1:10]
colnames(a.indel)<-a.indel[,"SNP"]
genotypes<-a.indel[5415:5420,the.samples.use]

genotypes<-g.indel[19617:19620,the.samples.use]

genotypes<-a.indel[41458:41460,the.samples.use]
genotypes<-g.indel[41458:41460,the.samples.use]
#all.genotypes<-filtered.genotype.summary(a.indel[41458:41460,],fam[,1],prefix="",suffix="",20,0.2,0.80,0.05,0.95,10,5)
genotypes<-filtered.genotype(a.indel[41458:41460,],fam[,1],prefix="",suffix="",20,0.2,0.80,0.05,0.95,10,5)
genotypes<-genotypes[,the.samples.use]

#a.indel[a.indel=="NA"]<-NA
genotypes[genotypes=="0/0"]<-0
genotypes[genotypes=="0/1"]<-1
genotypes[genotypes=="1/1"]<-2

genotypes[1:3,1:15]
dim(genotypes)
dim(fam)
num.samples<-length(the.samples.use)

dim(genotypes)
dim(fam)
num.samples
apply(genotypes,1,function(x) sum(as.numeric(x),na.rm=TRUE))
logistic.run(genotypes, fam,"pheno",num.samples)




require(doMC)
num.cores<-16
registerDoMC(cores=num.cores)


if(dim(a.indel)[1]<200){
num.bits<-1; print(paste("Only",dim(a.indel)[1],"genotypes remaining",sep=" ")) }else{num.bits<-num.cores}  # dim(genotypes)  genotypes[1:5,]

fil.genotypes<-foreach(a.indel.bit=iter(a.indel,by='row',chunksize=as.integer(dim(a.indel)[1]/num.bits)), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype.summary(a.indel.bit,fam[,1],prefix="",suffix="",20,0.2,0.80,0.05,0.95,10,5)



##              V1           V2 V3 V4 V5         V6
## 1  AOGC-01-0004 AOGC-01-0004  0  0 -9 -1.1744368
## 2  AOGC-02-0001 AOGC-02-0001  0  0 -9 -1.3493510
## 3  AOGC-02-0004 AOGC-02-0004  0  0 -9 -1.9811935
## 4  AOGC-02-0006 AOGC-02-0006  0  0 -9 -0.2589538
## 5  AOGC-02-0007 AOGC-02-0007  0  0 -9 -1.2621613
## 6  AOGC-02-0009 AOGC-02-0009  0  0 -9 -1.0031078
## 7  AOGC-02-0010 AOGC-02-0010  0  0 -9 -0.7473245
## 8  AOGC-02-0012 AOGC-02-0012  0  0 -9 -1.0981900
## 9  AOGC-02-0013 AOGC-02-0013  0  0 -9 -1.3716957
## 10 AOGC-02-0015 AOGC-02-0015  0  0 -9 -2.3591726
## 11 AOGC-02-0016 AOGC-02-0016  0  0 -9 -1.1362144
## 12 AOGC-02-0018 AOGC-02-0018  0  0 -9 -1.3437141
## 13 AOGC-02-0019 AOGC-02-0019  0  0 -9 -1.1903162
## 14 AOGC-02-0022 AOGC-02-0022  0  0 -9 -1.2017207
## 15 AOGC-02-0023 AOGC-02-0023  0  0 -9 -1.4214780
##   rownames(test) <- indels[,"SNP"]
num.samples<-dim(fam)[1]
#
res <- logistic.run(genotypes[41458:41460,], fam,"pheno",num.samples)
res <- logistic.run(genotypes, fam,"pheno",num.samples)

write.table(fam,file=the.fam,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)



ip<-8
projects

genotypes

   UQDIexomechipWithAOGC_chr5:14358303:14358303:A:G-1_B_R_2098209433 

grep("41858623",a.indel[,"start"])
grep("60148436",g.indel[,"start"])
28436472

exm826163_ver4
[2,]    "10:60148436:C:G" "0"                   "916" "-0.141094604875547" "0.0501590271836513" "-2.81294540181064" "0.00501410692669864"
[2,] "chr10:60148436:60148436:C:G:snp" "0"                   "916" "-0.141094604875547" "0.0501590271836513" "-2.81294540181064" "0.00501410692669864" Filtered
[2,] "chr10:60148436:60148436:C:G:snp" "0.000514933058702369" "971" "-0.974365644383531" "1.52536201288705"  "-0.638776655083571" "0.523119187980993" Unfiltered
9.855e-09:0.9994

a.indel[41459,1:50]



res


g



## ## bim match
## grep 41858623:41858623 BMD_EFF_STD_HIP.AOGC_vcf_final_chr17.qassoc
##  CHR                                                                                           SNP         BP    NMISS       BETA         SE         R2        T            P 
##   17                                                                                  17:41858623:41858623:A:G:snp   41858623      969     0.1737     0.1232   0.002051     1.41        0.159 


## grep 41858623 BMD_EFF_STD_HIP.recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.qassoc
##   17     UQDIexomechipWithAOGC_chr17:41858623:41858623:G:A-1_T_F_2098253894   41858623     6916     0.1397     0.0279   0.003613    5.007    5.673e-07

## ## cluster looks ok
## 0.02995 fro seq using logistic- seq center as covar
##  0.4957 using older residuals
## same with non-standard (except for BETA)
## Checked two residual files are ok


#### why did this not replicate????















########################################################## single point for AOGC
########################################################## single point for AOGC
########################################################## single point for AOGC
########################################################## single point for AOGC
########################################################## single point for AOGC
########################################################## single point for AOGC

## Imputed genotypes against exome

## /media/UQCCG-Analysis/AOGC_exome_chip/AOGC_ImputedExomeOnly/

## /media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/ForearmFrLoTrauma/ForearmFrLoTrauma.logistic.covar.assoc.logistic.ann

## /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2012-06-18_All_AOGC/PassSNPsConvertedToPlink/cd 

## prior to EMCC_BP_GRP.AOGC_to_vcf_filtered rooted


fam.template.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink/AOGC_vcf_final_chr1.fam"
project.prefix<-"AOGC_vcf_final_chr"

fam.template.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink/AOGC_to_vcf_filtered_10.fam"
project.prefix<-"AOGC_to_vcf_filtered_"

annotation.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
contaminated.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/contaminated_AOGC_SEQ_samples.txt"

## annotation.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_SAMPLES_PHENOTYPES_Nov.1.2013_RESIDUALS.txt"
## annotation.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/pheno.sequence_HIP_LS_FN_FX.txt"


#sum(ann[,"BMD_EFF_STD_HIP"] !=  ann1[,"BMD_EFF_STD_HIP"])

#the.chr<-basename(fam.template.file)
the.fam.dir<-dirname(fam.template.file)


files<-dir(the.fam.dir)  
projects<-gsub(".bed","",files[grepl(".bed$",files)])  # files[grepl(".bed$",files)] ### these are sequencing ones

projects<-projects[grepl(paste("^",project.prefix,sep=""),projects)]
projects

fam<-read.table(fam.template.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
contaminated<-read.table(contaminated.file,header=F,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

#traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
## traits<-c("TOT_HIP_GCM","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")

## traits<-c("BMD_EFF_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")


traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")

names(traits)<-rep("logistic",times=length(traits))
names(traits)[1:3]<-rep("assoc",times=3)
traits
traits %in% colnames(ann)




#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")



fam[1:5,]

i<-1
setwd(the.fam.dir)

#for (i in 1:length(traits)){
## for (i in 1:12){
## for (i in 1:3){
## for (i in 11:17){
for (i in 18:length(traits)){
## print(logistic.files[i])
## logistic<-read.table(logistic.files[i],header=T,fill=TRUE,stringsAsFactors=FALSE)

## colnames(logistic)[colnames(logistic)=="BP"]<-"POS"

## if("ADD" %in% colnames(logistic)){
## wanted<-grepl("ADD",logistic[,"TEST"])
## logistic<-logistic[wanted,]
## }

## posns<-match(logistic[,"SNP"],ann.short[,"SNP"])
## test<-cbind(logistic,ann.short[posns,]) # ,stringsAsFactors=FALSE

## dim(test)
  
fam[,6]<--9
posns<-match(fam[,1],ann[,"PATIENT"])
missing<-is.na(posns)
print(sum(missing))

no.pheno<-ann[posns,traits[i]]==3 | ann[posns,traits[i]]==-9 | is.na(ann[posns,traits[i]])  # ann[posns,traits[i]]
bad.sample<-fam[,1] %in% contaminated[,1] # contamited in sequecing
print(sum(no.pheno))
sum(bad.sample)

missing<-missing | no.pheno | bad.sample
cbind(fam[no.pheno,],ann[posns[no.pheno],c("PATIENT",traits[i])])
print(sum(missing))

if(names(traits)[i]=="logistic"){
 using.covars<-c("AGE_SCAN","PCA1","PCA2","PCA3","PCA4") #  using.covars<-c("CENTRE","WEIGHT","AGE_SCAN","PCA1","PCA2","PCA3","PCA4")
  covars<- ann[posns,c("PATIENT","PATIENT", using.covars)]
 colnames(covars)[1:2]<-c("FID","IID")
 covars[1:5,]
  missing.covars<-is.na(covars[,using.covars])
   missing.covars<-apply(missing.covars,1,sum,na.rm=TRUE)
 missing.covars<-  missing.covars>0
missing<-missing | missing.covars
write.table(covars[!missing,],file=paste("covars",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
}


has.missing<-sum(missing)>0
has.missing
write.table(fam[missing,c(1,1)],file=paste("remove_from",traits[i],sep="."),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

fam[!missing,6]<-ann[posns[!missing],traits[i]]
the.fam<-paste(traits[i],"fam",sep=".")

fam[!missing,][1:15,]
ann[posns[!missing],c("PATIENT",traits[i])][1:5,]

##              V1           V2 V3 V4 V5         V6
## 1  AOGC-01-0004 AOGC-01-0004  0  0 -9 -1.1744368
## 2  AOGC-02-0001 AOGC-02-0001  0  0 -9 -1.3493510
## 3  AOGC-02-0004 AOGC-02-0004  0  0 -9 -1.9811935
## 4  AOGC-02-0006 AOGC-02-0006  0  0 -9 -0.2589538
## 5  AOGC-02-0007 AOGC-02-0007  0  0 -9 -1.2621613
## 6  AOGC-02-0009 AOGC-02-0009  0  0 -9 -1.0031078
## 7  AOGC-02-0010 AOGC-02-0010  0  0 -9 -0.7473245
## 8  AOGC-02-0012 AOGC-02-0012  0  0 -9 -1.0981900
## 9  AOGC-02-0013 AOGC-02-0013  0  0 -9 -1.3716957
## 10 AOGC-02-0015 AOGC-02-0015  0  0 -9 -2.3591726
## 11 AOGC-02-0016 AOGC-02-0016  0  0 -9 -1.1362144
## 12 AOGC-02-0018 AOGC-02-0018  0  0 -9 -1.3437141
## 13 AOGC-02-0019 AOGC-02-0019  0  0 -9 -1.1903162
## 14 AOGC-02-0022 AOGC-02-0022  0  0 -9 -1.2017207
## 15 AOGC-02-0023 AOGC-02-0023  0  0 -9 -1.4214780



write.table(fam,file=the.fam,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)



ip<-1
projects
#for(ip in 11:11){
for(ip in 1:length(projects)){
print(projects[ip])
#setwd(projects[ip])

the.bed<-paste(projects[ip],"bed",sep=".")
the.bim<-paste(projects[ip],"bim",sep=".")


if(names(traits)[i]=="assoc"){
  system( paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam,"--remove",paste("remove_from",traits[i],sep="."),"--assoc","--allow-no-sex","--out",paste(traits[i],projects[ip],sep="."),"--noweb",sep=" ") )
}

  ## system( paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam,"--remove",paste("remove_from",traits[i],sep="."),"--covar",paste("covars",traits[i],sep="."),"--covar-name CENTRE,WEIGHT,AGE_SCAN,PCA1,PCA2,PCA3,PCA4","--linear --ci 0.95 --vif 500","--allow-no-sex","--out",paste(traits[i],projects[ip],sep="."),"--noweb",sep=" ") )



if(names(traits)[i]=="logistic"){
  system( paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam,"--remove",paste("remove_from",traits[i],sep="."),"--covar",paste("covars",traits[i],sep="."),"--covar-name AGE_SCAN,PCA1,PCA2,PCA3,PCA4","--logistic --ci 0.95 ","--allow-no-sex","--out",paste(traits[i],projects[ip],sep="."),"--noweb",sep=" ") )
}
} # ip loop over chromosomes
} # i loop over traits







#/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/clean_AOGC_SEQ_samples.txt.csv

################ hardy for filtered
########################################################## single point for AOGC



fam.template.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink/AOGC_vcf_final_chr1.fam"
project.prefix<-"AOGC_vcf_final_chr"

fam.template.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink/AOGC_to_vcf_filtered_10.fam"
project.prefix<-"AOGC_to_vcf_filtered_"

annotation.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/clean_AOGC_SEQ_samples.txt.csv

## annotation.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_SAMPLES_PHENOTYPES_Nov.1.2013_RESIDUALS.txt"
## annotation.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/pheno.sequence_HIP_LS_FN_FX.txt"


#sum(ann[,"BMD_EFF_STD_HIP"] !=  ann1[,"BMD_EFF_STD_HIP"])

#the.chr<-basename(fam.template.file)
the.fam.dir<-dirname(fam.template.file)


files<-dir(the.fam.dir)  
projects<-gsub(".bed","",files[grepl(".bed$",files)])  # files[grepl(".bed$",files)] ### these are sequencing ones

projects<-projects[grepl(paste("^",project.prefix,sep=""),projects)] # 
projects

fam<-read.table(fam.template.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)


#traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
## traits<-c("TOT_HIP_GCM","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")

## traits<-c("BMD_EFF_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")


traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")

names(traits)<-rep("logistic",times=length(traits))
names(traits)[1:3]<-rep("assoc",times=3)
traits
traits %in% colnames(ann)


#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")




i<-1
setwd(the.fam.dir)

for (i in 1:length(traits)){

  
ip<-1
projects
#for(ip in 11:11){
for(ip in 1:length(projects)){
print(projects[ip])
#setwd(projects[ip])

the.bed<-paste(projects[ip],"bed",sep=".")
the.bim<-paste(projects[ip],"bim",sep=".")
the.fam<-paste(traits[i],"fam",sep=".")


  system( paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam,"--remove",paste("remove_from",traits[i],sep="."),"--hardy","--allow-no-sex","--out",paste(traits[i],projects[ip],sep="."),"--noweb",sep=" ") )




} # ip loop over chromosomes
} # i loop over traits


















##################################### AOGC beta rescale
#res <- logscan(genotypes, a.fam,"pheno",num.samples)
# logistic.run(genotypes[41458:41460,], fam,"pheno",num.samples)
#genotypes, fam,"pheno",num.samples)
logistic.run<- function(gdata, pdata,pcolumn,num.samp) {
 ## gdata <- genotypes
 ## pdata <- fam
 ## pcolumn <- "pheno"
 ## num.samp <- num.samples

  
  ## print(dim(gdata))
  ## print(rownames(gdata))
  n <- dim(gdata)[1]
##   if(is.null(n)){ # only one element
##     gdata<-as.matrix(gdata)
##     if(dim(gdata)[1]!=1){gdata<-t(gdata)}
## }
#  n <- names(gdata)
  markers <- length(n) - 1

  all.missing<-apply(gdata,1,function(x) sum(is.na(as.numeric(x))))
  all.missing<-grep(dim(gdata)[2],all.missing,fixed=TRUE)
  do.rows<-1:dim(gdata)[1]
  if(length(all.missing)>0){do.rows<-do.rows[!(do.rows %in% all.missing)]}
  ## all.missing[1:5]
  ## sum(missing)


  
  results <- matrix(data=NA,nrow=nrow(gdata),ncol=7)
  #results <- matrix(nrow=ncol(gdata),ncol=5)
  colnames(results) <- c("SNP","REFfreq","ngeno","Estimate", "Std Error", "Z", "pval" )
#  i<-3
  for( i in 1:length(do.rows)) {
 #   print(i)
    
 #   for( i in 2:length(n)) {  
   # pdata[,"missing"]<-as.numeric(gdata[i,])
    ## for aogc & hbm use both
  #  fit <- lm(gdata[i,] ~ pdata[,"PCA1"]+ pdata[,"PCA2"] + pdata[,"PCA3"] + pdata[,"PCA4"]+pdata[,pcolumn])
   fit <- lm(as.numeric(pdata[,pcolumn])~as.numeric(gdata[do.rows[i],]) ) # same results as plink
 #  fit <- lm(as.numeric(gdata[do.rows[i],]) ~ as.numeric(pdata[,pcolumn]) )
    a <- summary(fit)
    
    ### get coeff for pheno1 only
#     Estimate <- a$coefficients[2,1]
#     Std_error <- a$coefficients[2,2]
#     T <- a$coefficients[2,3]
#     pval <- a$coefficients[6,1:4]
      the.geno<-as.numeric(gdata[do.rows[i],])
      the.geno<-the.geno[!is.na(the.geno)]
      p= sum(the.geno)/(2*length(the.geno))
  # p= 1- sum(as.numeric(gdata[do.rows[i],]),na.rm=TRUE)/(2*num.samp)

    if( dim(a$coefficients)[1]==1){ # no slope (monomorphic typically)  so get P if intercept - bad! as dim(a$coefficients)[1]==1 NOT 2 (slope and intercept)
      results[do.rows[i],] <- c(rownames(gdata)[do.rows[i]],p,length(the.geno),rep(NA,times=4))
      }else{
    # results[do.rows[i],] <- c(rownames(gdata)[do.rows[i]],p,length(the.geno),a$coefficients[dim(a$coefficients)[1],1:4])    #a$coefficients[c(Estimate, Std_error, T, pval )
     results[do.rows[i],] <- c(rownames(gdata)[do.rows[i]],p,length(the.geno),a$coefficients["as.numeric(gdata[do.rows[i], ])",1:4]) ## just in case ordr changes
   
     }
  
  }

  if(length(all.missing)>0){results[all.missing,"SNP"]<-rownames(gdata)[all.missing]}
  #rownames(results) <- seq(from = 1, to = markers)
  #colnames(results) <- c("Estimate", "Std Error", "Z", "pval" )
 
  return(results)
  
}   






logscan<- function(gdata, pdata,pcolumn,num.samp) {
 ## gdata <- genotypes[c(41,54,76,77,81),]# [1,]
 ## pdata <- a.fam
 ## pcolumn <- "pheno"
 ## num.samp <- num.samples

  
  ## print(dim(gdata))
  ## print(rownames(gdata))
  n <- dim(gdata)[1]
##   if(is.null(n)){ # only one element
##     gdata<-as.matrix(gdata)
##     if(dim(gdata)[1]!=1){gdata<-t(gdata)}
## }
#  n <- names(gdata)
  markers <- length(n) - 1

  all.missing<-apply(gdata,1,function(x) sum(is.na(x)))
  all.missing<-grep(dim(gdata)[2],all.missing,fixed=TRUE)
  do.rows<-1:dim(gdata)[1]
  if(length(all.missing)>0){do.rows<-do.rows[!(do.rows %in% all.missing)]}
  ## all.missing[1:5]
  ## sum(missing)


  
  results <- matrix(data=NA,nrow=nrow(gdata),ncol=7)
  #results <- matrix(nrow=ncol(gdata),ncol=5)
  colnames(results) <- c("SNP","REFfreq","ngeno","Estimate", "Std Error", "Z", "pval" )
#  i<- 1
  for( i in 1:length(do.rows)) {
 #   print(i)
    

  #  plot(as.numeric(gdata[do.rows[i],]), as.numeric(pdata[,pcolumn])) 
  #  fit <- lm(gdata[i,] ~ pdata[,"PCA1"]+ pdata[,"PCA2"] + pdata[,"PCA3"] + pdata[,"PCA4"]+pdata[,pcolumn])
 #  fit <- lm(as.numeric(pdata[,pcolumn])~as.numeric(gdata[do.rows[i],]) ) # same results as plink
   fit <- lm(as.numeric(gdata[do.rows[i],]) ~ as.numeric(pdata[,pcolumn]) )
    a <- summary(fit)
    
    ### get coeff for pheno1 only
#     Estimate <- a$coefficients[2,1]
#     Std_error <- a$coefficients[2,2]
#     T <- a$coefficients[2,3]
#     pval <- a$coefficients[6,1:4]
      the.geno<-as.numeric(gdata[do.rows[i],])
      the.geno<-the.geno[!is.na(the.geno)]
      p= sum(the.geno)/(2*length(the.geno))
    
  # p= 1- sum(as.numeric(gdata[do.rows[i],]),na.rm=TRUE)/(2*num.samp)
     if( dim(a$coefficients)[1]==1){ # no slope (monomorphic typically)  so get P if intercept - bad! as dim(a$coefficients)[1]==1 NOT 2 (slope and intercept)
      results[do.rows[i],] <- c(rownames(gdata)[do.rows[i]],p,length(the.geno),rep(NA,times=4))
      }else{
    # results[do.rows[i],] <- c(rownames(gdata)[do.rows[i]],p,length(the.geno),a$coefficients[dim(a$coefficients)[1],1:4])    #a$coefficients[c(Estimate, Std_error, T, pval )
     results[do.rows[i],] <- c(rownames(gdata)[do.rows[i]],p,length(the.geno),a$coefficients["as.numeric(pdata[, pcolumn])",1:4]) ## just in case ordr changes
   
     }
    results[do.rows[i],] <- c(rownames(gdata)[do.rows[i]],p,length(the.geno),a$coefficients[dim(a$coefficients)[1],1:4])    #a$coefficients[c(Estimate, Std_error, T, pval )
  
  }

  if(length(all.missing)>0){results[all.missing,"SNP"]<-rownames(gdata)[all.missing]}
  #rownames(results) <- seq(from = 1, to = markers)
  #colnames(results) <- c("Estimate", "Std Error", "Z", "pval" )
 
  return(results)
  
}   



code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir)
source("annotate_SNPs_subroutines.r") # has same version as the above two









annotation.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

genotype.file.location<-"/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink"
project.prefix<-"AOGC_vcf_final_chr"

genotype.file.location<-"/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink"
project.prefix<-"AOGC_to_vcf_filtered_"

the.seq.dir<-"/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink" # dirname(fam.template.file)
files<-dir(the.seq.dir)  
projects<-gsub(".bed","",files[grepl(".bed$",files)])  # files[grepl(".bed$",files)]
projects


projects<-projects[grepl(paste("^",project.prefix,sep=""),projects)] # 
projects

traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")

names(traits)<-rep("logistic",times=length(traits))
names(traits)[1:3]<-rep("assoc",times=3)
traits
traits %in% colnames(ann)


#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
names(traits)<-rep(".assoc.logistic",times=length(traits))
names(traits)[1:3]<-rep(".qassoc",times=3)
traits


library(doMC)
num.cores<-7
registerDoMC(cores=num.cores)

i<-2

setwd(genotype.file.location)
  

ip<-11
projects
#for(ip in 11:11){
for(ip in 1:length(projects)){
print(projects[ip])
#setwd(projects[ip])


a.indel<-read.plink.raw(paste(genotype.file.location,projects[ip],sep="/")) #A1->ALT POS->start ALT is 1 REF 0
dim(a.indel)

###  write.plink(a.indel[132080:132090,] ,"test")

print("start genotype conversion QC")
#a.indel[a.indel=="NA"]<-NA
a.indel[a.indel=="0/0"]<-0
a.indel[a.indel=="0/1"]<-1
a.indel[a.indel=="1/1"]<-2

#all.possible.samples<-gsub(".GT$","",colnames(a.indel)[grep(".GT$",colnames(a.indel))],perl=TRUE)

for (i in 1:3){

the.fam<-paste(traits[i],"fam",sep=".")
a.fam<-read.table(the.fam,header=F,fill=TRUE,stringsAsFactors=FALSE)
colnames(a.fam)<-c("PATIENT","SAMPLE","mother","father","sex","pheno")
a.fam<-a.fam[a.fam[,"pheno"]!=-9,]
a.fam[1:5,]
the.samples<-paste(a.fam[,"PATIENT"],".GT",sep="")
genotypes<-a.indel[,the.samples]
rownames(genotypes)<-a.indel[,"SNP"]



## the.bed<-paste(projects[ip],"bed",sep=".")
## the.bim<-paste(projects[ip],"bim",sep=".")

##   colnames(test) <- sample.geno.order
##   rownames(test) <- indels[,"SNP"]
num.samples<-dim(a.fam)[1]
#
# res <- logscan(genotypes, a.fam,"pheno",num.samples)
#foreach(genotypes.bit=iter(genotypes,by='row',chunksize=as.integer(dim(genotypes)[1]/num.bits)), .combine='cbind', .multicombine=TRUE,.inorder=TRUE) %dopar%  print(dim(genotypes.bit))

if(dim(genotypes)[1]<200){
num.bits<-1; print(paste("Only",dim(genotypes)[1],"genotypes remaining",sep=" ")) }else{num.bits<-num.cores}  # dim(genotypes)  genotypes[1:5,]

res<-foreach(genotypes.bit=iter(genotypes,by='row',chunksize=as.integer(dim(genotypes)[1]/num.bits)), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% logscan(genotypes.bit,a.fam,"pheno",num.samples)

write.table(res,file=paste("rescale",traits[i],projects[ip],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


} # i loop over traits
} # ip loop over chromosomes


######################### RESCALE
######################### RESCALE
######################### RESCALE
######################### RESCALE
######################### RESCALE
######################### RESCALE
######################### RESCALE
######################### RESCALE
######################### RESCALE










annotation.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)







the.seq.dir<-"/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink" # dirname(fam.template.file)
files<-dir(the.seq.dir)  
projects<-gsub(".bed","",files[grepl(".bed$",files)])  # files[grepl(".bed$",files)]


project.prefix<-"AOGC_vcf_final_chr"
project.prefix<-"AOGC_to_vcf_filtered_"

projects<-projects[grepl(paste("^",project.prefix,sep=""),projects)] # 
projects


traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")

names(traits)<-rep("logistic",times=length(traits))
names(traits)[1:3]<-rep("assoc",times=3)
traits
traits %in% colnames(ann)


#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
names(traits)<-rep(".assoc.logistic",times=length(traits))
names(traits)[1:3]<-rep(".qassoc",times=3)
traits

i<-1
for (i in 1:3){
  
ip<-1
logistic.seq<-{}
for(ip in 1:length(projects)){

the.fam<-paste(traits[i],"fam",sep=".")
a.fam<-read.table(the.fam,header=F,fill=TRUE,stringsAsFactors=FALSE)
colnames(a.fam)<-c("PATIENT","SAMPLE","mother","father","sex","pheno")
a.fam<-a.fam[a.fam[,"pheno"]!=-9,]
a.fam[1:5,]

the.samples<-a.fam[,"PATIENT"]
#the.samples<-paste(a.fam[,"PATIENT"],".GT",sep="")
print(length(the.samples))

logistic.file.seq<-paste(traits[i],projects[ip],sep=".")
print(logistic.file.seq)
setwd(the.seq.dir)
a.logistic.seq<-read.table(paste(logistic.file.seq,names(traits)[i],sep=""),header=T,fill=TRUE,stringsAsFactors=FALSE)

colnames(a.logistic.seq)[colnames(a.logistic.seq)=="BP"]<-"POS"

if("TEST" %in% colnames(a.logistic.seq)){
wanted<-grepl("ADD",a.logistic.seq[,"TEST"])
a.logistic.seq<-a.logistic.seq[wanted,]
}


a.logistic.seq[1:5,]
traits[i]
missing<-is.na(ann[,traits[i]])
gos<-grepl("^GOS",ann[,"PATIENT"])
## sum(gos)
## ann[gos,c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN")][1:5,]
extreme<-ann[,"PATIENT"] %in% a.fam[,"PATIENT"]

var<-sd(ann[!missing & !extreme,traits[i]])^2


rescale.file<-paste("rescale",traits[i],projects[ip],sep=".")
res<-read.delim(file=rescale.file,header=T,fill=TRUE,stringsAsFactors=FALSE)

posns<-match(a.logistic.seq[,"SNP"],res[,"SNP"])
missing<-is.na(posns)
sum(missing)
res<-res[posns,]

res[1:5,]
a.logistic.seq[1:5,]

p<-as.numeric(res[,"REFfreq"])
p.scale<-2*p*(1-p)
BETA.TRUE <- (as.numeric(res[,"Estimate"])*var)/p.scale
SE.TRUE <- (as.numeric(res[,"Std.Error"])*var)/p.scale
#BETA <- as.numeric(a.logistic.seq[,"BETA"] )
#cbind(a.logistic.seq[,c("SNP","BETA","SE","P")],res[,c("SNP","REFfreq","Estimate","Std.Error")],p.scale,var,BETA.TRUE,SE.TRUE)[1:5,]


a.logistic.seq[,"BETA"]<-BETA.TRUE
a.logistic.seq[,"SE"]<-SE.TRUE

recaled.logistic.file.seq<-paste("NewBETA",traits[i],projects[ip],sep=".")
print(recaled.logistic.file.seq)
setwd(the.seq.dir)


write.table(a.logistic.seq,file=recaled.logistic.file.seq,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

}}



################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis

################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis


################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis

################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
#### all these in


annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)



the.seq.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/AOGC_vcf_to_plink" # dirname(fam.template.file)
files<-dir(the.seq.dir)  
projects<-gsub(".bed","",files[grepl(".bed$",files)])  # files[grepl(".bed$",files)]


project.prefix.old<-"AOGC_vcf_final_chr" ## old run moved to /media/UQCCG-Analysis/AOGC_exome_chip/meta_analysis_single_point/non-filtered-seq
output.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/notFiltered"

project.prefix<-"AOGC_to_vcf_filtered_" ## outpt used same name as above
output.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/filtered"

projects.old<-projects[grepl(paste("^",project.prefix.old,sep=""),projects)]
projects<-projects[grepl(paste("^",project.prefix,sep=""),projects)] # 
projects
projects.old

the.chr.old<-gsub(project.prefix.old,"",projects.old)
the.chr<-gsub(project.prefix,"",projects)

sum(the.chr.old !=the.chr) # should be ZERO else out of order




########### projects, pojects.old ref.cohort all in the same order
projects

traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")

names(traits)<-rep("logistic",times=length(traits))
names(traits)[1:3]<-rep("assoc",times=3)
traits
traits %in% colnames(ann)


#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
names(traits)<-rep(".assoc.logistic",times=length(traits))
names(traits)[1:3]<-rep(".qassoc",times=3)
traits

# Corresdonding in same order for exome chip
 logistic.project.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/"
 logistic.project<-"recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL"
logistic.files<-paste(logistic.project.dir,traits,".",logistic.project,names(traits),sep="")
logistic.files
## logistic.files<-c(
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/Hip/totalhipBMDassocinexcel.txt",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/LS/LSBMDassocinexcel.txt",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/FN/FemNeckBMDassocinexcel.txt",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/EverFx/everFx_final_analysis.logistic.covar.assoc.logistic",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/everFx18nottrivia/everFx18nottrivia.LOGISTIC.COVAR.assoc.logistic",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/AllFr_50_OP/AllFr_50_OP.logistic.covar.assoc.logistic",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/AllVertFx_OP/AllVertFX_OP.logistic.covar.assoc.logistic",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/AllNonVertFX_50_OP/AllNonVertFX_50_OP.logistic.covar.assoc.logistic",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/ForearmFrLoTrauma/ForearmFrLoTrauma.logistic.covar.assoc.logistic",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/AllHipFr_50_OP/AllHipFr_50_OP.logistic.covar.assoc.logistic",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/EMcC_BP_GRP/EMcC_BP_GRP.logistic.covar.assoc.logistic"
##                   )

#annotation_short_exome_chip.txt

ann.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/annoatation_short_exome_chip.txt"
ann.short<-read.delim(ann.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

## below commented out was to get the rs and design score into ann.short
## ann.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/annoatation_extensive_exome_chip.txt"
## ann.ex<-read.delim(ann.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
## ann.ex<-ann.ex[,c(1:6,25,27,7,73,107)]
## colnames(ann.ex)[11]<-"rs.id"
## colnames(ann.ex)[7]<-"OMIM"
## colnames(ann.ex)
## colnames(ann.short)

## /media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP_lists_used_in_exomeChip/UQDIexomechipWithAOGC_By ExistingDesign - WGGT Final_103877.score.csv
## read.delim(ann.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

## design.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP_lists_used_in_exomeChip/UQDIexomechipWithAOGC_By ExistingDesign - WGGT Final_103877.score.csv"


## design<-read.delim(design.file,header=T,sep=",",skip=20,fill=TRUE,stringsAsFactors=FALSE)
## colnames(design)
## dim(design) #25879
## design[1:5,]
## design[1:5,c("Locus_Name","Final_Score","Normalization_Bin","Bead_Types.Assay")] 
## ann.ex[1:5,]
## ann.short[1:5,]

## ## posns<-match(ann.ex[,"SNP"],design[,"Locus_Name"])
## ## missing<-is.na(posns)
## ## sum(!missing) #23212
## ## design<-design[posns,]
## ## ann.ex<-cbind(ann.ex,design[,c("Locus_Name","Final_Score","Normalization_Bin","Bead_Types.Assay")])
## ann.uqdi<-ann.short[grepl("^UQDI",ann.short[,"SNP"]),]
## dim(ann.uqdi) #23121
## posns<-match(design[,"Locus_Name"],ann.uqdi[,"SNP"])
## missing<-is.na(posns)
## sum(missing) #23212

## design[,"Locus_Name"]
## write.table(design[missing,],file="/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/uqdi_missing-design.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE) 
## ann.ex[posns,c("rs.id","Final_Score","Normalization_Bin","Bead_Types.Assay","OMIM")][1:5,]


## ann.short<-cbind(ann.short,ann.ex[posns,c("rs.id","Final_Score","Normalization_Bin","Bead_Types.Assay","OMIM")]

##  table(ann.short[,"Normalization_Bin"])
##  table(design[,"Normalization_Bin"])

##  write.table(ann.short,file="/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/annoatation_short_exome_chip.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE) 
                 

########################## get processed bim files that are all on the forward strand in annovar format"
code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")

######

bim.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.bim"
bim.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.bim"
sample.files<-basename(bim.file)
names(sample.files)<-"GENO"
snp.dir<-dirname(bim.file)
sample.grs<-paste(names(sample.files),"grs",sep=".")
sample.grs
setwd(code.dir)
source("read.plink.bim.files.r")
bim.chip<-GENO.grs
dim(bim.chip)
bim.chip[1:5,]


# use
# pleo@bioinform01:/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink$ cat *.bim >  AOGC_vcf_final_chrALL.bim
bim.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink/AOGC_vcf_final_chrALL.bim"
bim.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink/AOGC_to_vcf_filtered_chrALL.bim"



sample.files<-basename(bim.file)
names(sample.files)<-"GENO"
snp.dir<-dirname(bim.file)
sample.grs<-paste(names(sample.files),"grs",sep=".")
sample.grs
sample.files
setwd(code.dir)
source("read.plink.bim.files.r")
bim.seq<-GENO.grs
dim(bim.seq)
bim.seq[1:5,]


##############################these bim files have the minor allele as A1 

bim.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.bim"
bim.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.bim"
bim<-read.table(bim.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
colnames(bim)<-c("chr","SNP","cm","start","ALT","REF")
bim.chip<-bim

## dim(bim.chip)
## dim(bim.chip1)
## diff<-(bim.chip[,"ALT"]!=bim.chip1[,"A1"])
## bim.chip[1:5,]
## bim.chip1[1:5,]
## sum(diff)
## cbind(bim.chip,bim.chip1)[diff,][1:10,]


bim.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink/AOGC_to_vcf_filtered_chrALL.bim"
bim<-read.table(bim.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
colnames(bim)<-c("chr","SNP","cm","start","ALT","REF") # A1 is labelled ALT here
bim.seq<-bim

bim.chip[1:5,]
bim.seq[1:5,]

#########

# save(list=c("bim.seq","bim.chip"),file="/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_A1_A2.RData")
# save(list=c("bim.seq","bim.chip"),file="/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_from_filtered_VCF.RData")
#save(list=c("bim.seq","bim.chip"),file="/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_fromVCF.RData")
##      chr  SNP                                           junk POS         REF ALT
## [1,] "12" "12:9009999:9009999:NA:-:indel"               "0"  "9009999"   NA  "-"
## [2,] "14" "14:102547769:102547769:NA:-:indel:102547769" "0"  "102547769" NA  "-"
## [3,] "16" "16:70011947:70011947:NA:-:indel:70011947"    "0"  "70011947"  NA  "-"
## [4,] "17" "17:38153714:38153714:NA:-:indel:38153714"    "0"  "38153714"  NA  "-"
## [5,] "20" "20:62606067:62606067:NA:-:indel"             "0"  "62606067"  NA  "-"
## [6,] "2"  "2:39210719:39210719:NA:-:indel:39210719"     "0"  "39210719"  NA  "-"

## bim.chip<-read.delim(bim.chip.file,header=F,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
i<-2
projects

for (i in 1:length(traits)){

#for (i in 1:12){
#  bim.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.bim"
#load("/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_fromVCF.RData")
load("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_A1_A2.RData")
print(logistic.files[i])
logistic<-read.table(logistic.files[i],header=T,fill=TRUE,stringsAsFactors=FALSE)

colnames(logistic)[colnames(logistic)=="BP"]<-"POS"

if("TEST" %in% colnames(logistic)){
wanted<-grepl("ADD",logistic[,"TEST"])
logistic<-logistic[wanted,]
}
dim(logistic)
dim(bim.chip)


posns<-match(logistic[,"SNP"],bim.chip[,"SNP"])
missing<-is.na(posns)
sum(missing)
bim.chip<-bim.chip[posns,]

#core.ann<-c("chr","start","end","REF","ALT")
core.ann<-c("chr","start","REF","ALT")
key.logistic<-build.key(bim.chip,core.ann)
dups<-duplicated(key.logistic)
sum(dups) #804
## dups<-unique(key.logistic[dups])
## dups<-key.logistic %in% dups
## logistic[dups,][1:50,] #### dups all have the same valuse so can remove any one 
bim.chip<-bim.chip[!dups,]
logistic<-logistic[!dups,]
key.logistic<-key.logistic[!dups]


rownames(logistic)<-key.logistic
TYPE<-rep("GENO",times=length(key.logistic))
logistic<-cbind(logistic,TYPE,bim.chip[, core.ann],stringsAsFactors=FALSE)

## diff<-(logistic[,"A1"]!=bim.chip[,"ALT"])
## sum(diff)
## logistic[diff,][1:20,]

ok<-logistic[,"CHR"]==2 & ( as.numeric(logistic[,"start"]) > 119000000 &  as.numeric(logistic[,"start"]) < 119729378)
sum(ok)
out<-logistic[ok,]
posns<-match(out[,"SNP"],ann.short[,"SNP"])
missing<-is.na(posns)
sum(missing)
out<-cbind(out,ann.short[posns,])
out
write.table(out,file="AOGC-LS-replication.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

## test<-key.logistic %in% c("chr10:101117384:101117384:G:A", "chr10:102307870:102307870:A:G")
## logistic[test,]
## bim.chip[test,]

## dim(test)

ip<-1 # i<-5
logistic.seq<-{}


###################read in all the sequence data
for(ip in 1:length(projects)){

setwd(the.seq.dir)
  
if(i %in% 1:3){ # use rescalled data
#  NewBETA
logistic.file.seq<-paste("NewBETA",traits[i],projects[ip],sep=".")
print(logistic.file.seq)
setwd(the.seq.dir)
a.logistic.seq<-read.table(logistic.file.seq,header=T,fill=TRUE,stringsAsFactors=FALSE)
}else{
logistic.file.seq<-paste(traits[i],projects[ip],sep=".")
print(logistic.file.seq)
setwd(the.seq.dir)
a.logistic.seq<-read.table(paste(logistic.file.seq,names(traits)[i],sep=""),header=T,fill=TRUE,stringsAsFactors=FALSE)
}


colnames(a.logistic.seq)[colnames(a.logistic.seq)=="BP"]<-"POS"

if("TEST" %in% colnames(a.logistic.seq)){
wanted<-grepl("ADD",a.logistic.seq[,"TEST"])
a.logistic.seq<-a.logistic.seq[wanted,]
}



##############


if(ip==1){logistic.seq<-a.logistic.seq}else{logistic.seq<-rbind(logistic.seq,a.logistic.seq)}


}## loop over ip
#setwd(projects[ip])

dim(logistic.seq)
dim(bim.seq)
#colnames(logistic.seq)[colnames(logistic.seq)=="BP"]<-"POS"

logistic.seq[1:5,]
logistic[1:5,]
bim.seq[1:5,]

posns<-match(logistic.seq[,"SNP"],bim.seq[,"SNP"]) # for annotation
missing<-is.na(posns)
sum(missing)
logistic.seq<-logistic.seq[!missing,]
bim.seq<-bim.seq[posns[!missing],]



dim(logistic.seq)
dim(bim.seq)
#logistic.seq[missing,][1:10,]
#bim.seq[1:5,]
################ remove flat genotypes  ### actaully let use there and reject the others using dup below and flat wlays first
## is.flat<-grepl(":flat$",bim.seq[,"SNP"]) # is.flat<-grepl(":flat$",logistic.seq[,"SNP"])
## sum(is.flat)
## ## bim.seq[is.flat,"SNP"][1:15]
## bim.seq<-bim.seq[!is.flat,]
## logistic.seq<-logistic.seq[!is.flat,]

core.ann<-c("chr","start","REF","ALT")
key.logistic.seq<-build.key(bim.seq,core.ann)
dups<-duplicated(key.logistic.seq)
sum(dups)
if(sum(dups)>0){
  
## dups<-unique(key.logistic.seq[dups])
## dups<-key.logistic.seq %in% dups
#logistic.seq[dups,][1:50,] #### dups all all unwound element not caught with end 
bim.seq<-bim.seq[!dups,]
logistic.seq<-logistic.seq[!dups,]
key.logistic.seq<-key.logistic.seq[!dups]
# key.logistic.seq[dups][1:10]
}



rownames(logistic.seq)<-key.logistic.seq
TYPE<-rep("SEQ",times=length(key.logistic.seq))
logistic.seq<-cbind(logistic.seq,TYPE,bim.seq[, core.ann],stringsAsFactors=FALSE)



logistic.seq[1:5,]
logistic[1:5,]

###########################
if("A1" %in% colnames(logistic)){ ### make sure ALT is A1 so effect direction is consistant - dont need to change sign of dirn then

ref.diff<-logistic[,"A1"]!=logistic[,"ALT"]
sum(ref.diff)
logistic[ref.diff,][1:10,]
the.ref<-logistic[ref.diff,"REF"]
logistic[ref.diff,"REF"]<-logistic[ref.diff,"ALT"]
logistic[ref.diff,"ALT"]<-the.ref
ref.diff<-logistic[,"A1"]!=logistic[,"ALT"]
sum(ref.diff)



sum(is.na(logistic.seq[,"A1"]))
remove<-is.na(logistic.seq[,"A1"]) | is.na(logistic.seq[,"ALT"]) |  is.na(logistic.seq[,"REF"]) # crap results
logistic.seq<-logistic.seq[!remove,]

ref.diff<-logistic.seq[,"A1"]!=logistic.seq[,"ALT"] 
sum(ref.diff)
logistic.seq[ref.diff,][1:10,]
the.ref<-logistic.seq[ref.diff,"REF"]
logistic.seq[ref.diff,"REF"]<-logistic.seq[ref.diff,"ALT"]
logistic.seq[ref.diff,"ALT"]<-the.ref
ref.diff<-logistic.seq[,"A1"]!=logistic.seq[,"ALT"]
sum(ref.diff)

# redo.start.end.annovar 

}

### using minor allele now so should be fineetal will sort out differences







logistic.seq[1:5,]
logistic[1:5,]
################### match up snp ids:
posns<-match(rownames(logistic.seq),rownames(logistic))
missing<-is.na(posns)
sum(!missing)
## uqdi<-grepl("^UQDI",logistic[,"SNP"])
## sum(uqdi)
## sum(missing)
## test<-missing & uqdi
## sum(test)
##  cbind(logistic.seq[!missing,],logistic[posns[!missing],])[1:5,]
logistic.seq[!missing,"SNP"]<-logistic[posns[!missing],"SNP"]

## logistic[1:5,]
## logistic.seq[1:5,]
 sum(is.na(logistic[,"chr"]))
 sum(is.na(logistic.seq[,"chr"]))
#order.by<-order(as.numeric(logistic[,"P"]),decreasing=FALSE)


#logistic<-logistic[order.by,]

if(!("BETA" %in% colnames(logistic))){
BETA=log(as.numeric(logistic[,"OR"]))
logistic<-cbind(logistic,BETA,stringsAsFactors=FALSE)


## exp(BETA +- 1.96SE)=95% CI
## so
## se<--(log(as.numeric(logistic[,"L95"]))-BETA)/1.96
## logistic<-cbind(logistic,BETA,se,stringsAsFactors=FALSE)
## and se = SE in plink so !! SE Is the standard error of BETA!! no need to mess with
}

if(!("BETA" %in% colnames(logistic.seq))){
BETA=log(as.numeric(logistic.seq[,"OR"]))
logistic.seq<-cbind(logistic.seq,BETA,stringsAsFactors=FALSE)

#logistic.seq[,"SE"]<-log(as.numeric(logistic.seq[,"SE"]))
}
print("PAST BETA")
logistic[posns[!missing],][1:5,]
logistic.seq[!missing,][1:5,]

setwd(output.dir)

write.table(logistic,file=paste("chip",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(logistic.seq,file=paste("seq",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

write.table(logistic[posns[!missing],],file=paste("chip.common",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(logistic.seq[!missing,],file=paste("seq.common",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


##print(doing merge)
common.cols<-colnames(logistic)[colnames(logistic) %in% colnames(logistic.seq)]
the.merge<-rbind(logistic[,common.cols],logistic.seq[,common.cols])
the.merge[1:5,]
order.by<-order(the.merge[,"CHR"],the.merge[,"start"])
the.merge<-the.merge[order.by,]
the.merge[1:5,]
write.table(the.merge,file=paste("merge",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
##

system(paste("sed s/TRAIT/",traits[i],"/ run_config_for_P_and_BETA_TEMPLATE.txt > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("./metal ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))

system(paste("sed s/TRAIT/common.",traits[i],"/ run_config_for_P_and_BETA_TEMPLATE.txt > ",paste("sample.size.CONFIG.common",traits[i],"txt",sep="."),sep=""))
system(paste("./metal ",paste("sample.size.CONFIG.common",traits[i],"txt",sep="."),sep=""))


system(paste("sed s/TRAIT/",traits[i],"/ run_config_for_inverse_varience_TEMPLATE.txt > ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))
system(paste("./metal ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))

system(paste("sed s/TRAIT/common.",traits[i],"/ run_config_for_inverse_varience_TEMPLATE.txt > ",paste("STDERR.CONFIG.common",traits[i],"txt",sep="."),sep=""))
system(paste("./metal ",paste("STDERR.CONFIG.common",traits[i],"txt",sep="."),sep=""))


} ## loop over traits







####################
###################




#########################  ANNOTATE META ANALYSIS/media/UQCCG-Analysis
#########################  ANNOTATE META ANALYSIS/media/UQCCG-Analysis
#########################  ANNOTATE META ANALYSIS/media/UQCCG-Analysis
#########################  ANNOTATE META ANALYSIS/media/UQCCG-Analysis
#########################  ANNOTATE META ANALYSIS/media/UQCCG-Analysis




#########################  ANNOTATE META ANALYSIS/media/UQCCG-Analysis
#########################  ANNOTATE META ANALYSIS/media/UQCCG-Analysis
#########################  ANNOTATE META ANALYSIS/media/UQCCG-Analysis
#########################  ANNOTATE META ANALYSIS/media/UQCCG-Analysis
#########################  ANNOTATE META ANALYSIS/media/UQCCG-Analysis



#traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")


ROOT.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/filtered"
fam.template.file<-"recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.fam"
unfiltered.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/notFiltered/"

code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")

annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")
#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
names(traits)<-rep(".assoc.logistic",times=length(traits))
names(traits)[1:3]<-rep(".qassoc",times=3)
traits




## hwe.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/zCall_AOGC.with.evoker_corrected_clean_FINAL.hwe"
## hwe<-read.table(hwe.file,header=T,fill=TRUE,stringsAsFactors=FALSE)
## wanted<-grepl("ALL",hwe[,"TEST"])
## hwe.ori<-hwe[wanted,]
## hwe.ori[1:5,]

## frq.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/zCall_AOGC.with.evoker_corrected_clean_FINAL.frq"
## frq.ori<-read.table(frq.file,header=T,fill=TRUE,stringsAsFactors=FALSE)
## frq.ori[1:5,]

## ann.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/annoatation_short_exome_chip.txt"
## ann.short<-read.delim(ann.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

## posns<-match(ann.short[,"SNP"],hwe.ori[,"SNP"])
## missing<-is.na(posns)
## sum(missing)
## hwe.ori<-hwe.ori[posns,]

## posns<-match(ann.short[,"SNP"],frq.ori[,"SNP"])
## missing<-is.na(posns)
## sum(missing)
## frq.ori<-frq.ori[posns,]

## ann.short.ori<-cbind(ann.short,hwe.ori,frq.ori,stringsAsFactors=FALSE)
## ann.short[1:5,]

the.chip.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/"
the.seq.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/AOGC_vcf_to_plink/" # dirname(fam.template.file)
files<-dir(the.seq.dir)  
projects<-gsub(".bed","",files[grepl(".bed$",files)])  # files[grepl(".bed$",files)]
project.prefix.old<-"AOGC_vcf_final_chr"
project.prefix<-"AOGC_to_vcf_filtered_"
projects.old<-projects[grepl(paste("^",project.prefix.old,sep=""),projects)]
projects<-projects[grepl(paste("^",project.prefix,sep=""),projects)] #

projects
projects.old

the.chr.old<-gsub(project.prefix.old,"",projects.old)
the.chr<-gsub(project.prefix,"",projects)

sum(the.chr.old !=the.chr)
the.chr<-the.chr[!(the.chr %in% c(23,24))]


### to use this would need to add reat the the single point sequening results to get the key right from the SNP name and then align those.
## ref.cohort<-paste("/media/UQCCG/Sequencing/Projects/AOGC-NGS/Australian Reference Genome Cohort/AOGC_Controls/AOGC-Genotyping.output.chr",the.chr,".AOGC_ALL.controls.txt",sep="")

## ref.cohort
## ip<-1
## for(ip in 1:length(ref.cohort)){
##   print(ref.cohort[ip])
## a.ref<- read.table(ref.cohort[ip],header=T,fill=TRUE,stringsAsFactors=FALSE,sep="\t",quote="")
## if(ip==1){ref.ann<-a.ref}else{ref.ann<-rbind(ref.ann,a.ref)}
## }
  
## ref.ann.key<-build.key(ref.ann,c("chr","start","end","REF","ALT","TYPE"))
## ref.ann.key[1:20]

########### need to load in seq from one of the runs below
## seq.key<-seq[,c("chr","start","start","REF","ALT","SNP")]
## colnames(seq.key)<-c("chr","start","end","REF","ALT","SNP")

## seq.key[,c("chr","start","end","REF","ALT")]<-redo.start.end.annovar(seq.key[,c("chr","start","end","REF","ALT")]) ## put in annovar format
## seq.key<-get.forward.strand.allele(seq.key)
## dim(seq.key)
## seq.key[1:5,]
## seq.key[,"chr"]<-gsub("^chr","",seq.key[,"chr"])
## seq.key.use<-build.key(seq.key,c("chr","start","end","REF","ALT"))
## ref.ann.key.use<-build.key(ref.ann,c("chr","start","end","REF","ALT"))

## seq.key.use[1:5]
## ref.ann.key.use[1:5]
## posns<-match(seq.key.use,ref.ann.key.use)
## missing<-is.na(posns)
## sum(missing)
## save(list=c("ref.ann.key","ref.ann.key.use","ref.ann","seq.key","seq.key.use"),file="/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/ref.ann.RData")

###### coomom ones miss chip and sequence data:
## media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/filtered$ grep exm-rs10048146 META.STDERR.common.BMD_EFF_STD_HIP1.tbl
## media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/filtered$ grep exm-rs10048146 META.STDERR.BMD_EFF_STD_HIP1.tbl
## exm-rs10048146	a	g	0.0531	0.0193	0.006038	+?


# Corresdonding in same order for exome chip
 logistic.project.dir<-ROOT.dir
 logistic.project<-gsub(".fam$","",fam.template.file)
 if(!grepl("/$",logistic.project.dir)){logistic.project.dir<-paste(logistic.project.dir,"/",sep="")}

#logistic.files<-paste(logistic.project.dir,traits,".",logistic.project,names(traits),sep="")
hwe.files.chip<-paste(the.chip.dir,traits,".",logistic.project,".hwe",sep="")
logistic.files<-paste(logistic.project.dir,"chip.",traits,sep="")
logistic.seq.files<-paste(logistic.project.dir,"seq.",traits,sep="")
logistic.seq.files.unfilt<-paste(unfiltered.dir,"seq.",traits,sep="")

bim.seq.files<-paste(projects,".bim",sep="")


#"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/filtered/"
meta.dir<-logistic.project.dir
meta.stderr.files<-paste(meta.dir,"META.STDERR.",traits,"1.tbl",sep="")
meta.size.files<-paste(meta.dir,"META.SAMPLE.SIZE.",traits,"1.tbl",sep="")


load("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_A1_A2_withKEY.RData") # contains "bim.chip" "bim.seq" "chipSeqkey" key using posn and ref/alt
load("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/ref.ann.RData")
bim.chip.ori<-bim.chip
bim.seq.ori<-bim.seq
ref.ann.key.ori<-build.key(ref.ann,c("chr","start","end","ALT","REF","TYPE"))
ref.ann.ori<-ref.ann

ann.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/annoatation_short_exome_chip.txt"
ann.short.ori<-read.delim(ann.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

gefos<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP_lists_used_in_exomeChip/GEFOS SNPs in final file.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
dim(gefos)


                                 
  i<-1

 all.meta<-{}

for (i in 1:length(traits)){
print(i)
print(traits[i])
  ## for (i in 12:17){
ann.short<-ann.short.ori
ref.ann.key<-ref.ann.key.ori
ref.ann<-ref.ann.ori
bim.seq<-bim.seq.ori
setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/AOGC_vcf_to_plink")
ip<-1

print(logistic.files[i])
chip<-read.table(logistic.files[i],header=T,fill=TRUE,stringsAsFactors=FALSE)
seq<-read.table(logistic.seq.files[i],header=T,fill=TRUE,stringsAsFactors=FALSE)
seq.unfilt<-read.table(logistic.seq.files.unfilt[i],header=T,fill=TRUE,stringsAsFactors=FALSE)
hwe.chip<-read.table(hwe.files.chip[i],header=T,fill=TRUE,stringsAsFactors=FALSE)
stderr<-read.table(meta.stderr.files[i],header=T,fill=TRUE,stringsAsFactors=FALSE)
sample.size<-read.table(meta.size.files[i],header=T,fill=TRUE,stringsAsFactors=FALSE)



ip<-1
for(ip in 1:length(projects)){
setwd(the.seq.dir)
############ read the hwe data
print(paste((traits)[i],projects[ip],"hwe",sep="."))

a.bim<-read.table(bim.seq.files[ip],header=F,fill=TRUE,stringsAsFactors=FALSE)

a.hwe<- read.table(paste((traits)[i],projects[ip],"hwe",sep="."),header=T,fill=TRUE,stringsAsFactors=FALSE)
a.hwe.old<- read.table(paste((traits)[i],projects.old[ip],"hwe",sep="."),header=T,fill=TRUE,stringsAsFactors=FALSE)

a.hwe<-as.matrix(a.hwe) ## colClasses="character"
a.hwe.old<-as.matrix(a.hwe.old)  
##############

if( sum( c("AFF","UNAFF") %in% a.hwe[,"TEST"] )==2 ){ # case/control unwind AFF and UNAFF
  aff<-a.hwe[,"TEST"]=="AFF"
  unaff<-a.hwe[,"TEST"]=="UNAFF"
  GENO.aff<-a.hwe[aff,"GENO"]
  GENO.unaff<-a.hwe[unaff,"GENO"]
 sum(aff)
 sum(unaff)
 if(sum(a.hwe[unaff,"SNP"]!=a.hwe[unaff,"SNP"])>0){print("Error in hwe alignments")}
  a.hwe<-cbind(a.hwe[unaff,c("CHR","SNP","A1","A2","P")],GENO.aff,GENO.unaff,stringsAsFactors=FALSE)
  aff<-a.hwe.old[,"TEST"]=="AFF"
  unaff<-a.hwe.old[,"TEST"]=="UNAFF"
  GENO.aff<-a.hwe.old[aff,"GENO"]
  GENO.unaff<-a.hwe.old[unaff,"GENO"]
 sum(aff)
 sum(unaff)
 if(sum(a.hwe.old[unaff,"SNP"]!=a.hwe.old[unaff,"SNP"])>0){print("Error in hwe alignments")}
  a.hwe.old<-cbind(a.hwe.old[unaff,c("CHR","SNP","A1","A2","P")],GENO.aff,GENO.unaff,stringsAsFactors=FALSE)
}

if(ip==1){hwe<-a.hwe}else{hwe<-rbind(hwe,a.hwe)}
if(ip==1){hwe.old<-a.hwe.old}else{hwe.old<-rbind(hwe.old,a.hwe.old)}
if(ip==1){bim.seq<-a.bim}else{bim.seq<-rbind(bim.seq,a.bim)}
#dim(a.logistic.seq)
} ## loop over ip


bim.seq[1:5,]
colnames(bim.seq)<-c("chr","SNP","junk","start","A1","A2")


if( sum( c("AFF","UNAFF") %in% hwe.chip[,"TEST"] )==2 ){ # case/control unwind AFF and UNAFF
  aff<-hwe.chip[,"TEST"]=="AFF"
  unaff<-hwe.chip[,"TEST"]=="UNAFF"
  GENO.aff<-hwe.chip[aff,"GENO"]
  GENO.unaff<-hwe.chip[unaff,"GENO"]
 sum(aff)
 sum(unaff)
 if(sum(hwe.chip[unaff,"SNP"]!=hwe.chip[unaff,"SNP"])>0){print("Error in hwe alignments")}
  hwe.chip<-cbind(hwe.chip[unaff,c("CHR","SNP","A1","A2","P")],GENO.aff,GENO.unaff,stringsAsFactors=FALSE)
}


## chip[1:5,]
## seq[1:5,]
## seq.unfilt[1:5,]
## sample.size[1:5,]
## stderr[1:5,]
## hwe.chip[1:5,]
## hwe[1:5,]
## hwe.old[1:5,]


dim(bim.seq)
dim(hwe)
dim(hwe.old)
dim(bim.seq)
dim(seq)
dim(stderr)
dim(sample.size)
dim(ref.ann)
dim(hwe)
dim(hwe.old)
#####################################################
### sort out hwe for chips ########
print("sort out hwe")

## sum(is.na(hwe.old[,"A1"]))
## sum(hwe.old[,"A1"]=="0")

## hwe.old.ori<-hwe.old #  hwe.old<-hwe.old.ori
## hwe.ori<-hwe #  hwe<-hwe.ori


hwe.old[hwe.old[,"A1"]=="0" & !is.na(hwe.old[,"A1"]),"A1"]<-"-"
hwe.old[hwe.old[,"A2"]=="0" & !is.na(hwe.old[,"A2"]),"A2"]<-"-"

hwe[hwe[,"A1"]=="0" & !is.na(hwe[,"A1"]),"A1"]<-"-"
hwe[hwe[,"A2"]=="0" & !is.na(hwe[,"A2"]),"A2"]<-"-"


#seq has exome chip ids use SNP names cause have flat genotypes and other problems
############## Sort out hwe for sequencing
hwe.key<-hwe[,"SNP"]
hwe.old.key<-hwe.old[,"SNP"]

### need to remove "end" which is variable based on allele length
hwe.key<-strsplit(hwe.key,split=":")
hwe.key<-unlist(lapply(hwe.key,function(x)  paste(c(x[1:2],x[4],x[5],x[6:length(x)]),collapse=":")))

hwe.old.key<-strsplit(hwe.old.key,split=":")
hwe.old.key<-unlist(lapply(hwe.old.key,function(x)  paste(c(x[1:2],x[4],x[5],x[6:length(x)]),collapse=":")))

hwe.key[1:5]
hwe.old.key[1:5]

sum(is.na(hwe.key))
sum(is.na(hwe.old.key))    
## test<-cbind(hwe.old.key,hwe.key)
## test[1:5,]

## test[grep("10:8106233",test[,1]),]

##      hwe.old.key                         hwe.key                            
## [1,] "10:8106233:-:T:indel:8106231:flat" "10:8106233:T:-:indel:8106231:flat"
## [2,] "10:8106233:T:-:indel:8106231"      "10:8106233:T:-:indel:8106231"     
## [3,] "10:8106233:T:-:indel:8106231"      "10:8106233:-:T:indel:8106231"

## hwe.old.key[grep("10:8106233:T:-:indel:8106231",hwe.old.key)]
## hwe.key[grep("10:8106233:T:-:indel:8106231",hwe.key)]
## hwe.old[grep("10:8106233:T:-:indel:8106231",hwe.old.key),]

posns<-match(hwe.old.key,hwe.key)
missing<-(is.na(posns))
sum(missing) #24279 HIP


########## hwe old has allele differences 
hwe.old.key.f<-strsplit(hwe.old[missing,"SNP"],split=":")
hwe.old.key.f[1:5]
hwe.old.key.f<-unlist(lapply(hwe.old.key.f,function(x)  paste(c(x[1:2],x[5],x[4],x[6:length(x)]),collapse=":")))
hwe.old.key[missing]<-hwe.old.key.f

hwe.old.flip<-missing


posns<-match(hwe.old.key,hwe.key)
missing<-is.na(posns)
sum(missing) #0

posns<-match(hwe.key,hwe.old.key)
missing<-is.na(posns)
sum(missing) #476

hwe.key.f<-strsplit(hwe[missing,"SNP"],split=":")
hwe.key.f[1:5]
hwe.key.f<-unlist(lapply(hwe.key.f,function(x)  paste(c(x[1:2],x[5],x[4],x[6:length(x)]),collapse=":")))
hwe.key[missing]<-hwe.key.f

hwe.flip<-missing

posns<-match(hwe.old.key,hwe.key)
missing<-is.na(posns)
print(sum(missing)) #0

posns<-match(hwe.key,hwe.old.key)
missing<-is.na(posns)
print(sum(missing)) #0


########### in general hwe shoudl be flipped in most cases

to.mark.cols<-colnames(hwe)[grep("^GENO",colnames(hwe))]
for(icol in 1:length(to.mark.cols)){
  hwe[hwe.flip,to.mark.cols[icol]]<-paste(hwe[hwe.flip,to.mark.cols[icol]],":MF",sep="")
  hwe.old[hwe.old.flip,to.mark.cols[icol]]<-paste(hwe.old[hwe.old.flip,to.mark.cols[icol]],":MF",sep="")
}
                                 
posns<-match(hwe.key,hwe.old.key)
missing<-is.na(posns)
sum(missing) #0
hwe.old<-hwe.old[posns,]
hwe.old.key<-hwe.old.key[posns]

hwe.key.ori<-hwe.key
hwe.old.key.ori<-hwe.old.key

sum(is.na(hwe.key))
sum(is.na(hwe.old.key))

sum(is.na( hwe.old[,"SNP"]))
sum(is.na( hwe[,"SNP"]))
## hwe.old.key[grep("10:93576",hwe.old.key)]


##  1 BMD_EFF_STD_HIP UQDIexomechipWithAOGC_chr10:93416:93416:G:A-1_T_F_2098216216
## 2 BMD_EFF_STD_HIP                                       10:93576:93576:A:C:snp
## 3 BMD_EFF_STD_HIP UQDIexomechipWithAOGC_chr10:93694:93694:T:C-1_T_R_2098216233
## 4 BMD_EFF_STD_HIP                                       10:94018:94018:C:T:snp
## 5 BMD_EFF_STD_HIP                                       10:94066:94066:T:C:snp
## hwe.key      hwe.key.ori
## 1 UQDIexomechipWithAOGC_chr10:93416:93416:G:A-1_T_F_2098216216 10:93416:A:G:snp
## 2                                                         <NA>             <NA>
## 3 UQDIexomechipWithAOGC_chr10:93694:93694:T:C-1_T_R_2098216233 10:93694:C:T:snp
## 4                                                         <NA>             <NA>
## 5                                                         <NA>             <NA>

##                                                             snps       rs.id
## 1 UQDIexomechipWithAOGC_chr10:93416:93416:G:A-1_T_F_2098216216 rs200558688
## 2                                       10:93576:93576:A:C:snp        <NA>
## 3 UQDIexomechipWithAOGC_chr10:93694:93694:T:C-1_T_R_2098216233 rs144539776
## 4                                       10:94018:94018:C:T:snp        <NA>
## 5                                       10:94066:94066:T:C:snp        <NA



## > hwe.old[grep("10:8106233:T:-:indel:8106231",hwe.old.key),]
##      CHR  SNP                                         TEST      A1  A2 
## [1,] "10" "10:8106233:8106233:-:T:indel:8106231:flat" "ALL(QT)" "-" "T"
## [2,] "10" "10:8106233:8106233:T:-:indel:8106231"      "ALL(QT)" "T" "-"
## [3,] "10" "10:8106233:8106233:T:-:indel:8106231"      "ALL(QT)" "T" "-"
##      GENO          O.HET.     E.HET.    
## [1,] "173/104/297" "0.181200" "0.476700"
## [2,] "21/39/514"   "0.067940" "0.131200"
## [3,] "235/85/254"  "0.148100" "0.499500



##########hwe and hwe.old now in same order
#############################################################
########### Get in same order as seq
print("order hwe to seq")

dim(bim.seq)
dim(hwe)
dim(hwe.old)
dim(stderr)
dim(seq)
dim(ann.short)

## sum(grepl("flat$",seq[,"SNP"]))

## ann.short[1:5,] ### use ann.short to match
## ref.ann[1:2,]
## hwe.key[1:5]
###############################
## aligen hwe to seq
## seq has flat
## seq has exome ID
## match all I can then match but position and alleles only

## hwe.old.ori.2<-hwe.old #  hwe.old<-hwe.old.ori.2
## hwe.ori.2<-hwe #  hwe<-hwe.ori.2
## hwe.key.ori<-hwe.key # hwe.key<-hwe.key.ori
## hwe.old.key.ori<-hwe.old.key # hwe.old.key<-hwe.old.key.ori

snps<-unique(c(seq$SNP,chip$SNP)) ## all possible markers
snps.chip<-unique(chip$SNP) ## ones on the chip
length(snps)


sum(is.na(seq[,c("SNP")]))
sum(is.na(chip[,c("SNP")]))
sum(is.na(snps))

snps.posn<-rbind(seq[,c("SNP","chr","start","REF","ALT")],chip[,c("SNP","chr","start","REF","ALT")])
posns<-match(snps,snps.posn[,"SNP"])
missing<-(is.na(posns))
sum(missing) # should be zero
snps.posn<-snps.posn[posns,]

sum(snps !=snps.posn[,"SNP"])
sum(is.na(snps.posn[,"SNP"]))



# cbind(snps,snps.posn)[1:5,]

#seq.key<-seq[,"SNP"]
seq.key<-snps


### need to remove end which is variable based on allele length
a.snp.name<-grepl("^rs",seq.key) | grepl("^UQDI",seq.key) | grepl("^exm",seq.key)
sum(a.snp.name)
 ##### do no ruin SNP ids                
seq.key.temp<-strsplit(seq.key[!a.snp.name],split=":")
seq.key.temp<-unlist(lapply(seq.key.temp,function(x) paste(c(x[1:2],x[4],x[5],x[6:length(x)]),collapse=":") ))
seq.key[!a.snp.name]<-seq.key.temp                     
seq.key[1:5]
hwe.key[1:5]
seq.key.temp[1:5]
seq.key[!a.snp.name][1:5]

#####################

posns<-match(hwe.key,seq.key)
missing<-(is.na(posns))
sum(missing)
sum(!missing)


#### make the keys for hwe the SNP names 
hwe.key[!missing]<-snps.posn[posns[!missing],"SNP"]
hwe.old.key[!missing]<-snps.posn[posns[!missing],"SNP"]

#### seq.key and hwe.key do not have "end" position
###

########################################
hwe.key.s<-strsplit(hwe.key[missing],split=":")
hwe.key.s<-unlist(lapply(hwe.key.s,function(x)  paste(c(x[1],x[2],x[3],x[4]),collapse=":")))


### need to remove end which is variable based on allele length
seq.key.s<-build.key(snps.posn,c("chr","start","ALT","REF"))

## hwe.key.ori<-hwe.key
## hwe.old.key.ori<-hwe.key.old

seq.key.s[1:5]
hwe.key.s[1:5]

posns.s<-match(hwe.key.s,seq.key.s)
missing.s<-is.na(posns.s)
sum(missing.s) #0

## sum(is.na(snps.posn[,c("SNP")]))
## sum(is.na(seq[,c("SNP")]))
## sum(is.na(chip[,c("SNP")]))
## sum(is.na(snps))
## sum(is.na(seq.key))
## sum(is.na(hwe.key))
## sum(is.na(hwe.key.s))
## sum(is.na(seq.key.s))
## length(hwe.key[missing][!missing.s])
## length(snps.posn[posns.s[!missing.s],"SNP"])
## hwe.key[missing][!missing.s][1000:1010]
## seq.key.s[posns.s[!missing.s]][1000:1010]
## snps.posn[posns.s[!missing.s],"SNP"][1000:1010]
## seq.key[posns.s[!missing.s]][1000:1010]
## hwe.old.key.t<-hwe.old.key #  hwe.old<-hwe.old.ori.2
## hwe.key.t<-hwe.key #

hwe.key[missing][!missing.s]<-snps.posn[posns.s[!missing.s],"SNP"]
hwe.old.key[missing][!missing.s]<-snps.posn[posns.s[!missing.s],"SNP"]


## sum(is.na(snps.posn[,c("SNP")]))
## sum(is.na(seq[,c("SNP")]))
## sum(is.na(chip[,c("SNP")]))
## sum(is.na(snps))
## sum(is.na(seq.key))
## sum(is.na(hwe.key))
## sum(is.na(hwe.key.s))
## sum(is.na(seq.key.s))
## sum(is.na(hwe.key[missing][!missing.s]))
## sum(is.na(seq.key.s[posns.s[!missing.s]]))
## sum(is.na(snps.posn[posns.s[!missing.s],"SNP"]))

## grep("10:92627",snps.posn[posns.s[!missing.s],"SNP"])
## grep("10:92627",snps.posn[,"SNP"])
## grep("10:92627",hwe.key[!missing])
## hwe.key[missing][!missing.s][1:5]
## ## hwe.ori.2<-hwe #
## hwe.key[1:5]
## seq.key[posns[!missing]][1:5]
## snps.posn[posns[!missing],][1:5,]
## sum(snps !=snps.posn[,"SNP"])
## sum(is.na(snps))
## sum(is.na(snps.posn[,"SNP"]))


posns<-match(hwe.key,snps.posn[,"SNP"])
missing<-is.na(posns)
print(sum(missing))



##  hwe.key[missing][1:10]
 ## [1] "10:92627:CTGC:-:indel"         "10:92756:A:G:snp"             
 ## [3] "10:92788:G:A:snp"              "10:92855:T:C:snp"             
 ## [5] "10:92864:G:A:snp"              "10:92927:A:T:snp"             
 ## [7] "10:92941:AGAACACAGTAA:-:indel" "10:92969:C:-:indel"           
 ## [9] "10:92981:T:C:snp"              "10:92990:T:G:snp" 


## if(sum(missing>0)){
## hwe.key.s<-strsplit(hwe.key[missing],split=":")
## hwe.key.s<-unlist(lapply(hwe.key.s,function(x)  paste(c(x[1],x[2],x[4],x[3]),collapse=":")))
## hwe.key.s[1:5]
## posns.s<-match(hwe.key.s,seq.key.s)
## missing.s<-is.na(posns.s)
## sum(missing.s)


## to.mark.cols<-colnames(hwe)[grep("^GENO",colnames(hwe))]
## for(icol in 1:length(to.mark.cols)){
##   hwe[hwe.flip,to.mark.cols[icol]]<-paste(hwe[hwe.flip,to.mark.cols[icol]],":MF",sep="")
##   hwe.old[hwe.old.flip,to.mark.cols[icol]]<-paste(hwe.old[hwe.old.flip,to.mark.cols[icol]],":MF",sep="")
## }

## }
## cbind(hwe.key[missing],hwe.key.s[missing],seq[posns[missing],])[1:10,]

############################### align ref.ann in the same was as hwe.key
#####################
print("sort out ref.ann")
sum(snps !=snps.posn[,"SNP"])

posns<-match(ref.ann.key,snps)
missing<-(is.na(posns))
sum(missing)
sum(!missing)

ref.ann.key[1:5]
seq.key[1:5]
#### make the keys for hwe the SNP names 
ref.ann.key[!missing]<-snps.posn[posns[!missing],"SNP"]


#### seq.key and hwe.key do not have "end" position
########################################
ref.ann.key.s<-strsplit(ref.ann.key[missing],split=":")
ref.ann.key.s<-unlist(lapply(ref.ann.key.s,function(x)  paste(c(x[1],x[2],x[4],x[5]),collapse=":")))

seq.key.s[1:5]
ref.ann.key.s[1:5]

posns.s<-match(ref.ann.key.s,seq.key.s)
missing.s<-is.na(posns.s)
sum(missing.s) #0
ref.ann.key.s[1:10]

ref.ann.key[missing][!missing.s]<-snps.posn[posns.s[!missing.s],"SNP"]

###################################### make second pass flipping alleles

posns<-match(ref.ann.key,snps)
missing<-(is.na(posns))
sum(missing)
sum(!missing)

ref.ann.key.s<-strsplit(ref.ann.key[missing],split=":")
ref.ann.key.s<-unlist(lapply(ref.ann.key.s,function(x)  paste(c(x[1],x[2],x[5],x[4]),collapse=":")))

seq.key.s[1:5]
ref.ann.key.s[1:5]

posns.s<-match(ref.ann.key.s,seq.key.s)
missing.s<-is.na(posns.s)
sum(missing.s) #0
ref.ann.key.s[1:10]
ref.ann.key[missing][!missing.s]<-snps.posn[posns.s[!missing.s],"SNP"]

ref.ann.key[1:5]
seq.key[1:5]

posns<-match(ref.ann.key,snps)
missing<-(is.na(posns))
sum(missing)
###########################ref.ann



##########hwe and hwe.old now in same order
#############################################################
########### Get in same order as seq
## sample.size[1:5,]
## stderr[1:5,]
## hwe.chip[1:5,]
## hwe[1:5,]
## hwe.old[1:5,]
## seq[1:5,]
## chip[1:5,]
## ann.short[1:5,]

## sum(snps !=snps.posn[,"SNP"])
## hwe.key #(has snps
## hwe.old.key # (has snps
## ref.ann.key # (has snps
#####################################################################################################################

posns<-match(ann.short[,"SNP"],snps) ## make chip the default id name before.
missing<-is.na(posns)
sum(missing) ##47 all markers in snps
### these are exome chip id that I have as UQDI

ann.short.key<-build.key(ann.short[missing,],c("chr","start","REF","ALT"))
#chip.key<-build.key(chip,c("chr","start","REF","ALT"))
posns.s<-match(ann.short.key,seq.key.s) ## seq.key.s, snp.key, snps, in same order
missing.s<-is.na(posns.s)
sum(missing.s) #30
#ann.short[missing,"SNP"][!missing.s] <-snps.posn[posns.s[!missing.s],"SNP"]
#ann.short.key[!missing.s]

########## these 17 need to be flipped
ann.short.key.f<-strsplit(ann.short.key[missing.s],split=":")
ann.short.key.f[1:5]
ann.short.key.f<-unlist(lapply(ann.short.key.f,function(x)  paste(c(x[1],x[2],x[4],x[3]),collapse=":")))
ann.short.key[missing.s]<-ann.short.key.f

### re match again
#posns.s<-match(ann.short.key,chip.key) ## make chip the default id name before.
posns.s<-match(ann.short.key,seq.key.s) ## make chip the default id name before.
missing.s<-is.na(posns.s)
sum(missing.s) # 0

ann.short[missing,"SNP"][!missing.s] <-snps.posn[posns.s[!missing.s],"SNP"]
################## ann.short now corrected 
posns<-match(ann.short[,"SNP"],snps) ## make chip the default id name before.
missing<-is.na(posns)
sum(missing)



## getwd() "/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/AOGC_vcf_to_plink"
## write.table(cbind(ann.short[missing,][!missing.s,1:10],chip[posns.s[!missing.s],]),file="exome.UQDI.dup.snps.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



############ use snps to aligne everthing against

print("match up all arrays")

chip[1:5,]
seq[1:5,]
seq.unfilt[1:5,]
sample.size[1:5,]
stderr[1:5,]
hwe.chip[1:5,]
hwe[1:5,]
hwe.old[1:5,]
seq.unfilt[1:5,]

dim(bim.seq)
dim(hwe)
dim(hwe.old)
dim(bim.seq)
dim(seq)
length(snps)
                                             

posns<-match(snps,seq[,"SNP"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
seq<-seq[posns,]

posns<-match(snps,seq.unfilt[,"SNP"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
seq.unfilt<-seq.unfilt[posns,]

posns<-match(snps,ann.short[,"SNP"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(missing)
ann.short<-ann.short[posns,]

posns<-match(snps,chip[,"SNP"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
chip<-chip[posns,]

posns<-match(snps,stderr[,"MarkerName"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
stderr<-stderr[posns,]

posns<-match(snps,sample.size[,"MarkerName"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
sample.size<-sample.size[posns,]


posns<-match(snps,hwe.chip[,"SNP"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
hwe.chip<-hwe.chip[posns,]


posns<-match(snps,hwe.key) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
hwe<-hwe[posns,]
hwe.old<-hwe.old[posns,]
hwe.key<-hwe.key[posns]
hwe.key.ori<-hwe.key.ori[posns]

posns<-match(snps,ref.ann.key) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
ref.ann<-ref.ann[posns,]




dim(seq.unfilt)
dim(seq)
dim(chip)
dim(stderr)
dim(sample.size)
dim(ann.short)
dim(hwe)
dim(hwe.old)
dim(hwe.chip)
dim(ref.ann)
dim(ann.short)
length(hwe.key)
length(hwe.key.ori)

P.value.StdErr<-stderr[,"P.value"]
P.value.Size<-sample.size[,"P.value"]
Weight<-sample.size[,"Weight"]

P.Seq.UnFilt<-seq.unfilt[,"P"]
P.Seq.Filt<-seq[,"P"]
#P.Seq<-rep(NA,times=length(P.Seq.Filt))
P.Chip<-chip[,"P"]

BETA.chip<-chip[,"BETA"]
BETA.seq<-seq[,"BETA"]

SE.chip<-chip[,"SE"]
SE.seq<-seq[,"SE"]

print("combine and print")
#### hwe by trait and pheno so can modify hwe.ori<-hwe ; hwe.old.ori<-hwe.old # hwe<-hwe.ori ; hwe.old<-hwe.old.ori
hwe<-hwe[,(grepl("^GENO",colnames(hwe)) |  colnames(hwe)=="P")]
hwe.old<-hwe.old[,(grepl("^GENO",colnames(hwe.old)) |  colnames(hwe.old)=="P")]
hwe.chip<-hwe.chip[,(grepl("^GENO",colnames(hwe.chip)) |  colnames(hwe.chip)=="P")]

colnames(hwe)[colnames(hwe)=="P"]<-"P.hardy"
colnames(hwe.old)[colnames(hwe.old)=="P"]<-"P.hardy"
colnames(hwe.chip)[colnames(hwe.chip)=="P"]<-"P.hardy"

colnames(hwe)<-paste(colnames(hwe),"Seq.Filt",sep=".")
colnames(hwe.old)<-paste(colnames(hwe.old),"Seq.UnFilt",sep=".")
colnames(hwe.chip)<-paste(colnames(hwe.chip),"Chip",sep=".")

all.hwe<-cbind(hwe,hwe.old,hwe.chip)

all.hwe[1:5,]
all.hwe<-as.matrix(all.hwe)

extra.cols<-c("GENO.Chip","GENO.Seq.Filt","GENO.Seq.UnFilt","GENO.aff.Chip","GENO.aff.Seq.Filt","GENO.unaff.Chip","GENO.unaff.Seq.Filt","GENO.aff.Seq","GENO.unaff.Seq","P.hardy.Chip","P.hardy.Seq.Filt","P.hardy.Seq.UnFilt")
extra<-matrix(data=NA,nrow=dim(all.hwe)[1],ncol=length(extra.cols))
colnames(extra)<-extra.cols
extra[,extra.cols[extra.cols %in% colnames(all.hwe)]]<-all.hwe[,extra.cols[extra.cols %in% colnames(all.hwe)]]
                                        # reorganise HWE
extra[1:5,]

dim(extra)
length(BETA.chip)


wanted.cols<-c("rs.id","design","refGene..type","Gene.Names","description","gerp.scores","PolyPhen.desc","PolyPhen.scores","SIFT.desc","SIFT.scores","skeletome","mouse.defect","sewell.cycling","Dequeant.cycling","ingenuity.bone.genes","Consequence.Embl","Amino_acids.Embl","MAF.lt.0.001","MAF.lt.0.5","Final_Score")


# wanted.cols %in% colnames(ann.short)

meta<-cbind(traits[i],stderr[,c("MarkerName","Allele1","Allele2","Effect","StdErr","Direction")],P.value.StdErr,P.value.Size,Weight,P.Chip,P.Seq.Filt,P.Seq.UnFilt,extra,hwe.key,hwe.key.ori,snps,ann.short[,wanted.cols],chip[,c("chr","start","REF","ALT","SNP")],BETA.chip,BETA.seq,SE.chip,SE.seq,ref.ann[,c("Hetero.ALT.reads","Hetero.REF.reads","Hetero.Read.Balance","FILTER","gerp.scores")],stringsAsFactors=FALSE)

## dim(meta)
## colnames(meta)


getwd() # "/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/AOGC_vcf_to_plink"
setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point")

out.file<-paste("META.ALL.",traits[i],".txt",sep="")
out.file.sig<-paste("META.SIG.",traits[i],".txt",sep="")
out.file.chip<-paste("META.ChipCentric.",traits[i],".txt",sep="")
out.file.gefos<-paste("META.GEFOS.",traits[i],".txt",sep="")

write.table(meta,file=out.file,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


posns<-match(snps,snps.chip)
missing<-is.na(posns)
sum(missing) ## shoudl be zero
sum(!missing)

write.table(meta[!missing,],file=out.file.chip,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


posns<-match(snps,gefos[,1])
missing<-is.na(posns)
sum(missing) ## shoudl be zero
sum(!missing)
write.table(meta[!missing,],file=out.file.gefos,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


#######################


sig.seq<-as.numeric(as.character(meta[,"P.Seq.Filt"]))< 0.001 & !is.na(as.numeric(as.character(meta[,"P.Seq.Filt"])))
sig.seq.unfilt<-as.numeric(as.character(meta[,"P.Seq.UnFilt"]))< 0.001 & !is.na(as.numeric(as.character(meta[,"P.Seq.UnFilt"])))
sig.chip<-as.numeric(as.character(meta[,"P.Chip"]))< 0.001 & !is.na(as.numeric(as.character(meta[,"P.Chip"])))
sig.sdterr<-as.numeric(as.character(meta[,"P.value.StdErr"]))< 0.0001 & !is.na(as.numeric(as.character(meta[,"P.value.StdErr"])))
sig.size<-as.numeric(as.character(meta[,"P.value.Size"]))< 0.0001 & !is.na(as.numeric(as.character(meta[,"P.value.Size"])))
is.UQDI<-grepl("^UQDI",snps)

possible.hits<-c("SRSF9","HELZ2","INTS6","ZKSCAN8","VWC2","SIRT1","HRNR","PCDH7","TCHH","AHNAK2","PGM5","HPX","PAPPA2","ARHGAP35","RERG","CFDP1","SRPX2","GATC::SRSF9","CNTNAP4")
# possible.hits %in% meta[,"Gene.Names"]
is.possible<-  meta[,"Gene.Names"] %in% possible.hits

## sum(is.possible)
## sum(sig.seq)
## sum(sig.seq.unfilt)
## sum(sig.chip)
## sum(sig.size)
## sum(sig.sdterr)
## sum(is.UQDI)

is.sig<-sig.seq | sig.seq.unfilt | sig.chip | sig.sdterr | sig.size |  is.UQDI | is.possible
## sum(is.sig)
write.table(meta[is.sig,],file=out.file.sig,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

############ put everything in the order of to.check


if(is.null(dim(all.meta))){
  all.meta<-meta[is.sig,]
}else{
  all.meta<-rbind(all.meta,meta[is.sig,])
}



}  ## loop over traits

getwd()
save(list=c("all.meta"),file="all.meta.sig.RData")
write.table(all.meta,file="META.all.sig.combined.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



######################################################################################################  FINAL REDO OF CLUSTER CHECKS

############################## compile merges list of all genotypes to check
############################## compile merges list of all genotypes to check
############################## compile merges list of all genotypes to check
############################## compile merges list of all genotypes to check
############################## compile merges list of all genotypes to check
############################## compile merges list of all genotypes to check
options(stringsAsFactors=FALSE)



already.checked<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/non-filtered-seq/META.BOTH.common.all_p0.0005.postED review.txt"
chked<-read.delim(already.checked,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)


################################compile list of checked genotypes
chked<-unique(chked[,"SNP"])
Plan<-rep("post.ED.1",times=length(chked))
chk.all<-cbind(chked,Plan)
chk.all[1:5,]

chked<-read.delim("/media/UQCCG/UQCCG-Projects/PAUL_LEO/AOGC exome chip core genotyping/set1-SNP.listre-checket/all.lists_filtered-clean.f.bim",header=F,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
chked[1:5,]
chked<-unique(chked[,2])
Plan<-rep("checked.ED.2",times=length(chked))
chked<-cbind(chked,Plan)
chk.all<-rbind(chk.all,chked)

chked<-read.delim("/media/UQCCG/UQCCG-Projects/PAUL_LEO/AOGC exome chip core genotyping/set2 -2014-06-23-HWE/all.lists_filtered-clean.f.bim",header=F,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
chked[1:5,]
chked<-unique(chked[,2])
Plan<-rep("HWE.3",times=length(chked))
chked<-cbind(chked,Plan)
chked[1:5,]
chk.all<-rbind(chk.all,chked)

chked<-read.delim("/media/UQCCG/UQCCG-Projects/PAUL_LEO/AOGC exome chip core genotyping/set3- 2014-06-23-Mhairi/all.lists_filtered-clean.f.bim",header=F,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
chked[1:5,]
chked<-unique(chked[,2])
Plan<-rep("Mh.rescore.4",times=length(chked))
chked<-cbind(chked,Plan)
chked[1:5,]
chk.all<-rbind(chk.all,chked)


#good.evoker<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/cluster_viz/good.evoker.snps",header=F,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
#good.evoker<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/cluster_viz/good.evoker.snps",header=F,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
chked<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/cluster_viz/good.evoker.snps",header=F,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
chked[1:5,]
chked<-unique(chked[,1])
Plan<-rep("good.evoker",times=length(chked))
chked<-cbind(chked,Plan)
chked[1:5,]
chk.all<-rbind(chk.all,chked)


#bad.evoker<-
chked<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/cluster_viz/bad.checked.snps",header=F,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
chked[1:5,]
chked<-unique(chked[,1])
Plan<-rep("bad.evoker",times=length(chked))
chked<-cbind(chked,Plan)
chked[1:5,]
chk.all<-rbind(chk.all,chked)

sum(duplicated(chk.all[,1]))
dups<-chk.all[duplicated(chk.all[,1]),1]
dups<-chk.all[chk.all[,1] %in% dups,]
dups<-tapply(dups[,2],dups[,1],function(x) paste(x,collapse=":"))

dups<-cbind(names(dups),dups)
colnames(dups)<-colnames(chk.all)

chk.all<-chk.all[!(chk.all[,1] %in% dups[,1]),]
chk.all<-rbind(chk.all,dups)

hits<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/punative.hit.list.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
hits[1:5,]
table(hits$Hit)
hits<-hits[hits$Hit !="",]
table(hits$Hit)
chk.all[1:5,]
dim(chk.all)

traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")
#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
names(traits)<-rep(".assoc.logistic",times=length(traits))
names(traits)[1:3]<-rep(".qassoc",times=3)
traits
i<-5



## out.file<-paste("META.ALL.",traits[i],".txt",sep="")
## out.file.sig<-paste("META.SIG.",traits[i],".txt",sep="")
## out.file.chip<-paste("META.ChipCentric.",traits[i],".txt",sep="")
## out.file.gefos<-paste("META.GEFOS.",traits[i],".txt",sep="")
## choose from above

in.files<-paste("META.SIG.",traits,".txt",sep="")
out.file<-paste("META.SIG.SUMMARY.txt")
options(scipen=4)

in.files<-paste("META.ChipCentric.",traits,".txt",sep="")
out.file<-paste("META.ChipCentric.SUMMARY.txt")
options(scipen=4)

in.files<-paste("META.GEFOS.",traits,".txt",sep="")
out.file<-paste("META.GEFOS.SUMMARY.txt")
options(scipen=4)


## seq.key.new<-build.key(ref.ann,c("chr","start","end","REF","ALT","TYPE"))
## seq.key.s<-build.key(ref.ann,c("chr","start","REF","ALT"))
## sum(!(seq.key[,"SNP"] %in% seq.key.s))

i<-1
setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point")
for (i in 1:length(traits)){
print(i)

in.file<-in.files[i]
a.sig<-read.delim(in.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
if(dim(a.sig)[1]==0){next}

colnames(a.sig)
colnames(a.sig)[1]<-"traits"
## chked[1:5,]
a.sig[1:5,]

########################### aligned checked all thiese are exome chip
## colnames(chked)
colnames(a.sig)
chk.all[1:5,]
posns<-match(chk.all[,"chked"],a.sig[,"snps"])
#posns<-match(chk.all[,"chked"],a.sig[,"MarkerName"])
sum(is.na(posns)) ## not all checked will be signifivant , SO EXPECT SOME MISSING (except with all)

posns<-match(a.sig[,"snps"],chk.all[,"chked"])
sum(is.na(posns))

chked<-chk.all[posns,"Plan"]

hits[1:5,]
posns<-match(a.sig[,"snps"],hits[,"MarkerName"])
Hit<-hits[posns,"Hit"]
table(Hit)
############################################################



a.sig<-cbind(chked,Hit,a.sig)

if(i==1){
  all.sig<-a.sig
}else{
    all.sig<-rbind(all.sig,a.sig)
  }

}  ### loop over traits




dim(all.sig)
all.sig[1:3,]

sum(is.na(all.sig$snps))
table(all.sig$traits)
order.by<-order(all.sig[,"snps"])

all.sig<-all.sig[order.by,]

lowest.p<-tapply(all.sig[,"P.value.StdErr"],all.sig[,"SNP"],function(x) min(as.numeric(x,na.rm=TRUE)))

#plan<-tapply(chked[,"Plan"],chked[,"SNP"],function(x) length(unique(x)))
## sum(plan!=1)  ## only one plan per snp
## plan[plan!=1]
## plan<-tapply(chked[,"Plan"],chked[,"SNP"],function(x) unique(x))
## unique(plan)
## plan[plan==""]<-"no change"
## plan<-cbind(names(plan),plan)
## colnames(plan)<-c("SNP","Plan")



length(lowest.p)
length(unique(all.sig[,"SNP"]))

lowest.p[1:40]
order.by<-order(as.numeric(lowest.p),decreasing=FALSE)
lowest.p<-lowest.p[order.by]

posns<-match(all.sig[,"SNP"],names(lowest.p))
length(posns)
posns[1:20]
dim(all.sig)
lowp<-cbind(names(lowest.p),lowest.p)
colnames(lowp)<-c("SNP","P")

## plan<-cbind(names(plan),plan)
## colnames(plan)<-c("SNP","Plan")
## plan[1:5,]

all.sig.f<-merge(all.sig,lowp,by="SNP",all.x=TRUE)
#all.sig.final<-merge(all.sig.f,plan,by="SNP",all.x=TRUE)
all.sig.final<-all.sig.f
all.sig.final[1:5,]

order.by<-order(as.numeric(as.character(all.sig.final[,"P"])),decreasing=FALSE)
all.sig.final[order.by,"P"][1:5]
all.sig.final[order.by,][1:5,]
all.sig.final<-all.sig.final[order.by,]

write.table(all.sig.final,file=out.file,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
getwd()

#/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/META.SIG.SUMMARY.txt
########################################################################################################################################


have<-unique(all.sig.final$traits)
have
## subset1 %in% have
## subset2 %in% have
## subset3 %in% have

## subset4<-have[!(have %in% c(subset1,subset2, subset3))]
## toString(subset4)

#/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/META.SIG.SUMMARY.QT.txt

subset<-c("BMD_AFFSTAT","BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN")
name<-"QT"
write.table(all.sig.final[(all.sig.final[,"traits"] %in% subset),],file=gsub(".txt",paste(".",name,".txt",sep=""),out.file),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


subset<-c("HIP_FR_50_OP","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","FR_50_OP","NEW_FR_50_OP_vs_no.adult.fx","NEW_FR_50_OP_vs_never_fx","EVER_FX")
name<-"FX"
write.table(all.sig.final[(all.sig.final[,"traits"] %in% subset),],file=gsub(".txt",paste(".",name,".txt",sep=""),out.file),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

subset<-c("NONVERT_OP_FX_50","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","VERT_FX_OP","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX")
name<-"VERT_FX"
write.table(all.sig.final[(all.sig.final[,"traits"] %in% subset),],file=gsub(".txt",paste(".",name,".txt",sep=""),out.file),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

subset<-c("NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","FR_18_NOT_TRIVIA","FOREARM_FR_LOTRAUMA","EMCC_BP_GRP","HIP_FR_50_EVER","RECODED_EMCC_BP_GRP")
name<-"FX_other"
write.table(all.sig.final[(all.sig.final[,"traits"] %in% subset),],file=gsub(".txt",paste(".",name,".txt",sep=""),out.file),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)







## (If this proves too big we could put LS by itself)

 

## 2.       HIP_FR_50_EVER

## By itself

 

## 3.       FOREARM_FR_LOTRAUMA

## NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX

 

## 4.       FR_50_OP

## NEW_FR_50_OP_vs_never_fx


## NEW_FR_50_OP_vs_no.adult.fx



 

## 5.
## NONVERT_OP_FX_50

## NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX

## NEW_NONVERT_OP_FX_50_VS_NEVER_FX N

 

## 6.        

## VERT_FX_OP


## NEW_VERT_FX_OP_VS_NEVER_FX


## NEW_VERT_FX_OP_VS_NO_ADULT_FX
















chipSeqkey[1:5,]
all.sig.final[1:3,]

unique(all.sig.final[,"SNP"])[1:50]
dim(all.sig.final)
table(all.sig.final[,"Direction"])
seq.only<-grepl("^\\?",all.sig.final[,"Direction"]) & is.na(all.sig.final[,"GENO"]) & !grepl("^e",all.sig.final[,"SNP"]) & !grepl("^UQ",all.sig.final[,"SNP"])
seq.only[test]

## grepl("^\\?",all.sig.final[test,"SNP"])
## is.na(all.sig.final[test,"GENO"])
## !grepl("^e",all.sig.final[test,"SNP"])
## !grepl("^UQ",all.sig.final[test,"SNP"])

sum(!seq.only)
table(all.sig.final[!seq.only,"Direction"])

sig<-as.numeric(as.character(all.sig.final[,"P"]))<0.00001 | as.numeric(as.character(all.sig.final[,"P.Chip"]))<0.00001 | is.na(as.numeric(as.character(all.sig.final[,"P.Chip"]))) | is.na(as.numeric(as.character(all.sig.final[,"P"])))
sum(!sig)
all.sig.final[!sig,][1:5,]
all.sig.final[!sig,][500:505,]

table(all.sig.final[,"Plan"])
reviewed<- !is.na(all.sig.final[,"Plan"])
sum(reviewed)


test<-grep("7:82595154:82595154:T:A:snp",all.sig.final[,"SNP"])

pase
sum((!seq.only & sig))
sum((!seq.only & sig) | reviewed)


##     wanted<-(!seq.only & !sig) & reviewed
## wanted[test]

sig[test]
!seq.only[test]
wanted<-(!seq.only & sig) | reviewed




## test<-grep("exm1346556",all.sig.final[,"SNP"])

##    all.sig.final[test,]  


              
    
unique(all.sig.final[,"Plan"])
all.sig.final[,"P.hardy.Seq.Filt"]<-as.numeric(as.character(all.sig.final[,"P.hardy.Seq.Filt"]))
all.sig.final[,"P.hardy.Seq"]<-as.numeric(as.character(all.sig.final[,"P.hardy.Seq"]))
all.sig.final[,"P.value.StdErr"]<-as.numeric(as.character(all.sig.final[,"P.value.StdErr"]))  

write.table(all.sig.final[wanted,],file="final.meta.analysis.single.point.significant.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
getwd()

              

all.sig.final[1:5,]
uqdi<-grepl("^UQDI",all.sig.final[,"SNP"])


              



              
write.table(all.sig.final[uqdi,],file="final.meta.analysis.single.point.uqdi.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
getwd()
##########################             
##########################

"MarkerName"
"MarkerName"
              

data<-read.delim("/media/UQCCG-Analysis/AOGC_exome_chip/meta_analysis_single_point/final.meta.analysis.single.point.significant.reorg.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

good.evoker<-read.delim("/media/UQCCG-Analysis/AOGC_exome_chip/cluster_viz/good.evoker.snps",header=F,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
bad.evoker<-read.delim("/media/UQCCG-Analysis/AOGC_exome_chip/cluster_viz/bad.checked.snps",header=F,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

data[1:5,]
good.evoker[1:5,]

posns<-match(data[,"SNP"],good.evoker[,1])
missing<-is.na(posns)
sum(!missing)
posns[!missing][1:50]

cbind(good.evoker[posns[!missing],1][1:5],data[!missing,][1:5,1:3])
dim(data)
length(posns)

checked<-!is.na(data[!missing,"Plan"])

data[!missing,][!checked,"Plan"]<-"evoker.good"
data[!missing,][checked,"Plan"]<-paste("evoker.good",data[!missing,][checked,"Plan"],sep="::")

posns<-match(data[,"SNP"],bad.evoker[,1])
missing<-is.na(posns)
sum(!missing) # 0 so all evoker have been removed

write.table(data,file="final.meta.analysis.single.point.significant.reorg.evoker.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
table(data[,"Plan"])

#/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/final.meta.analysis.single.point.significant.reorg.evoker.txt
########################
######################
#gene list for last pass

data<-read.delim("/media/UQCCG-Analysis/AOGC_exome_chip/meta_analysis_single_point/final.meta.analysis.single.point.significant.reorg.evoker.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)


dim(data)
colnames(data)

sum(is.na(data[,"P_min"]))/22

signif<-as.numeric(data[,"P_min"]) < 5e-5 | is.na(data[,"P_min"])
sum(signif)/22
designed<-!is.na(data[,"design"])
sum(designed)/22
no.plan<-is.na(data[,"Plan"])
sum(signif & designed & no.plan)

data<-data[signif & designed & no.plan,]
dim(data)
dups<-duplicated(data[,"SNP"])
sum(dups)
data<-data[!dups,]
dim(data)
write.table(data,file="final.meta.analysis.single.point.significant.reorg.evoker.to.REEVOKER.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

sum(is.na(data[,"chr"])
data[is.na(data[,"chr"]),]
data[is.na(data[,"chr"]),"chr"]<-6
data[is.na(data[,"chr"]),"start"]<-31234114


    
    
    table(data[,"chr"])
setwd("/media/UQCCG-Analysis/AOGC_exome_chip/meta_analysis_single_point")
the.chr<-names(table(data[,"chr"]))
for(i in 1:length(the.chr)){
  out.file<-paste("meta.snp.test_chr",the.chr[i],".txt",sep="")
  print(length((data[data[,"chr"]==the.chr[i],"SNP"])))
  write.table(data[data[,"chr"]==the.chr[i],"SNP"],file=out.file,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
}
  
getwd()











    

data[!missing,][checked,"Plan"]
data[checked,"Plan"]
grep("evoker.good",data[,"Plan"])[1:10]

best.col.order<-colnames(data)

/home/pleo/Dropbox/AOGC paper/final.meta.analysis.single.point.significant.xlsx

/media/UQCCG-Analysis/AOGC_exome_chip/meta_analysis_single_point/final.meta.analysis.single.point.significant.xlsx


uqdi<-grepl("^UQDI",data[,"MarkerName"])
no.geno<-grepl("^0/0",data[,"GENO"]) | grepl("^0/1",data[,"GENO"])
signif<-as.numeric(data[,"P_min"]) < 1e-6 & !is.na(data[,"P_min"])

sum(signif)
sum(no.geno)
sum(uqdi)
sum(!no.geno & uqdi & signif)
data[!no.geno & uqdi & signif,c("MarkerName","GENO","P_min","design")][1:200,]

test<-data[!no.geno & uqdi & signif,]
dups<-duplicated(test[,"MarkerName"])
test<-test[!dups,]
dim(test)
test

    

              
load("/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_fromVCF.RData")

all.sig.final[1:5,]


bim.seq[1:5,]
bim.chip[1:5,]

posns<-match(chip[,"SNP"],seq[,"SNP"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
common<-chip[!missing,"SNP"]


#Collect stuff that is
## 1)significant ORIGINALLY or in meta analysis
## 2)in common no matter what

## posns<-match(chip[,"SNP"],seq[,"SNP"]) ## make chip the default id name before.
## missing<-is.na(posns)
## sum(!missing)
## common<-chip[!missing,"SNP"]

## hist(as.numeric(seq[,"NMISS"]))
## hist(as.numeric(chip[,"NMISS"]))

num.samples<-max(as.numeric(seq[,"NMISS"]),na.rm=TRUE)
min.miss<-num.samples-0.2*num.samples
min(as.numeric(seq[,"NMISS"],na.rm=TRUE))
sig.seq<-as.numeric(as.character(seq[,"P"]))< 0.00001 & !is.na(as.numeric(as.character(seq[,"P"]))) & as.numeric(seq[,"NMISS"]) > min.miss
sum(sig.seq)
sig.seq<-seq[sig.seq,"SNP"]

num.samples<-max(as.numeric(chip[,"NMISS"]),na.rm=TRUE)
min.miss<-num.samples-0.2*num.samples
min(as.numeric(chip[,"NMISS"],na.rm=TRUE))
sig.chip<-as.numeric(chip[,"P"])< 0.0001 & !is.na(as.numeric(chip[,"P"])) & as.numeric(chip[,"NMISS"]) > min.miss
sum(sig.chip)
sig.chip<-chip[sig.chip,"SNP"]


to.check<-unique(c(common,sig.seq,sig.chip))

################# in exome chip but not sequencing
UQDIexomechipWithAOGC_chr13:51939613:51939613:C:T-1_T_R_2098221723: HIP
exm634087

## chr13	51936239	51939981	3743	+	chr13:51936240-51939981:INTS6	5119	43	0	5852	31	0	12	-1.967301147	0.049148503	ENSG00000102786	INTS6	protein_coding	integrator complex subunit 6 [Source:HGNC Symbol;Acc:14879]	NA	NA

############################################################

pleo@DI-LW-BRN-011:/media/UQCCG-Analysis/AOGC_exome_chip/meta_analysis_single_point/notFiltered$ head -1  META.BOTH.common.RESCALED.SIGNIF2.EVER_FX1.tblALL > EVER_FX.headpleo@DI-LW-BRN-011:/media/UQCCG-Analysis/AOGC_exome_chip/meta_analysis_single_point/notFiltered$ grep INST6 META.BOTH.common.RESCALED.SIGNIF2.EVER_FX1.tblALL > INTS6pleo@DI-LW-BRN-011:/media/UQCCG-Analysis/AOGC_exome_chip/meta_analysis_single_point/notFiltered$ wc -l INTS6

/media/old-scratch/media/Bioinform-D/Research/exome Chip$ grep INTS6 AOGC_Gene_Hits_0.05_0.2-sent.txt 

0 INTS6
pleo@DI-LW-BRN-011:/media/UQCCG-Analysis/AOGC_exome_chip/meta_analysis_single_point/notFiltered$ grep INTS6 META.BOTH.common.RESCALED.SIGNIF2.EVER_FX1.tblALL > INTS6
pleo@DI-LW-BRN-011:/media/UQCCG-Analysis/AOGC_exome_chip/meta_analysis_single_point/notFiltered$ wc -l INTS6
16 INTS6
pleo@DI-LW-BRN-011:/media/UQCCG-Analysis/AOGC_exome_chip/meta_analysis_single_point/notFiltered$ cat EVER_FX.head INTS6 > INTS6.txt
pleo@DI-LW-BRN-011:/media/UQCCG-Analysis/AOGC_exome_chip/meta_analysis_single_point/notFiltered$ 


    
setwd("/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink")
EVER_FX.AOGC_to_vcf_filtered_13.assoc.logistic

a.logistic.seq<-read.table("/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink/EVER_FX.AOGC_to_vcf_filtered_13.assoc.logistic",header=T,fill=TRUE,stringsAsFactors=FALSE)
a.logistic.seq<-read.table("/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink/BMD_AFFSTAT.AOGC_to_vcf_filtered_13.assoc.logistic",header=T,fill=TRUE,stringsAsFactors=FALSE)
a.logistic.seq[1:5,]

colnames(a.logistic.seq)[colnames(a.logistic.seq)=="BP"]<-"POS"
a.logistic.seq[grep("51939613",a.logistic.seq[,"POS"]),]

if("TEST" %in% colnames(a.logistic.seq)){
wanted<-grepl("ADD",a.logistic.seq[,"TEST"])
a.logistic.seq<-a.logistic.seq[wanted,]
}
a.logistic.seq[1:5,]

grep("51939613",a.logistic.seq[,"POS"])

a.logistic.seq[17410:17464,]

ref.ann[1:5,]

a.seq[1:5,]

    

FOR DESIGN INFO

## check region associated from origin design
## /media/old-scratch/media/Bioinform-D/Research/exome Chip
## grep SFRS9 AOGC_Gene_Hits_0.05_0.2-sent.txt
/media/old-scratch/media/Bioinform-D/Research/exome Chip/Final targets/For_Exomechip_FromPaul.csv

## check for pseudo genes: (origin design probes)
## /media/old-scratch/media/Bioinform-D/Research/exome Chip/Final targets
## grep chr12:120899502:120899502 WGGTScoreResults.v6.csv


### chedk to see if other loci in sequencing:
/media/UQCCG-Analysis/AOGC_exome_chip/meta_analysis_single_point/notFiltered
 head -1   META.BOTH.common.RESCALED.SIGNIF2.BMD_EFF_STD_HIP1.tblALL > head
grep SFRS9 META.BOTH.common.RESCALED.SIGNIF2.BMD_EFF_STD_HIP1.tblALL

  }



UQDIexomechipWithAOGC_chr12:120779978:120779978:C:A-1_T_F_2098221188

## UQDIexomechipWithAOGC_chr6:43255471:43255471:C:T-1_T_R_2098211709
## UQDIexomechipWithAOGC_chr5:14358303:14358303:A:G-1_B_R_2098209433

## pleo@bioinform01:/media/UQCCG-Analysis/AOGC_exome_chip/working_genotypes$ grep 14358303 exomeChip_working_chr5.bim
## 5	UQDIexomechipWithAOGC_chr5:14358303:14358303:A:G-1_B_R_2098209433	0	14358303	G	A



## pleo@bioinform01:/media/UQCCG-Analysis/AOGC_exome_chip/working_genotypes$ head -1 BMD_EFF_STD_HIP.recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.qassoc

##    5      UQDIexomechipWithAOGC_chr5:14358303:14358303:A:G-1_B_R_2098209433   14358303     6916     0.9427     0.2147   0.002781    4.391    1.146e-05 
##  CHR                                                                    SNP         BP    NMISS       BETA         SE         R2        T            P


## pleo@bioinform01:/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink$ grep 14358303:14358303 BMD_EFF_STD_HIP.AOGC_vcf_final_chr5.qassoc
##    5                                                                   5:14358303:14358303:G:A:snp   14358303      911    -0.1379      0.188  0.0005914  -0.7334       0.4635 
##  CHR                                                                                           SNP         BP    NMISS       BETA         SE         R2        T            P 

## pleo@bioinform01:/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink$ grep 14358303:14358303 AOGC_vcf_final_chr5.bim
## 5 5:14358303:14358303:G:A:snp 0 14358303 G A


## #####################
## UQDIexomechipWithAOGC_chr17:41858623:41858623:G:A-1_T_F_2098253894
## pleo@bioinform01:/media/UQCCG-Analysis/AOGC_exome_chip/AOGC_vcf_to_plink$ grep 41858623:41858623 AOGC_vcf_final_chr17.bim
## 17 17:41858623:41858623:A:G:snp 0 41858623 A G

## pleo@bioinform01:/media/UQCCG-Analysis/AOGC_exome_chip/working_genotypes$ grep 41858623 recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL_chr17.bim
## 17	UQDIexomechipWithAOGC_chr17:41858623:41858623:G:A-1_T_F_2098253894	0	41858623	A	G

## ## bim match
## grep 41858623:41858623 BMD_EFF_STD_HIP.AOGC_vcf_final_chr17.qassoc
##  CHR                                                                                           SNP         BP    NMISS       BETA         SE         R2        T            P 
##   17                                                                                  17:41858623:41858623:A:G:snp   41858623      969     0.1737     0.1232   0.002051     1.41        0.159 


## grep 41858623 BMD_EFF_STD_HIP.recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.qassoc
##   17     UQDIexomechipWithAOGC_chr17:41858623:41858623:G:A-1_T_F_2098253894   41858623     6916     0.1397     0.0279   0.003613    5.007    5.673e-07

## ## cluster looks ok
## 0.02995 fro seq using logistic- seq center as covar
##  0.4957 using older residuals
## same with non-standard (except for BETA)
## Checked two residual files are ok


#### why did this not replicate????
    a.fam.ori<-a.fam
    
annotation.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/pheno.sequence_HIP_LS_FN_FX_RAD_noHeight.txt" # used in GEFOS
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)    
    
fam.file1<-"/media/UQCCG/GWAS/GEFOS_uk10k/AOGC_HBM_Genotypes_ImputedUk10K/HBM/chr1/HBM_1kg_1.fam" 
fam.file2<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/AOGC_1kg_1.fam"
fam1<-read.table(fam.file1,header=F,fill=TRUE,stringsAsFactors=FALSE)
fam2<-read.table(fam.file2,header=F,fill=TRUE,stringsAsFactors=FALSE)
a.fam<-rbind(fam1,fam2)
    dim(a.fam)
    

annotation.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
contaminated.file<-"/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/contaminated_AOGC_SEQ_samples.txt"
contaminated<-read.table(contaminated.file,header=F,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

    


    


    forarm<-read.delim("/media/old-scratch/media/Bioinform-D/Research/GWAS extreme regression/AOGC_HBM/Forearm BMD phenotype file from BC.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) #
forarm[,1]
forarm[!is.na(forarm[,"DISTFOREARM_BMD"]),"RTOT_BMD"]<-forarm[!is.na(forarm[,"DISTFOREARM_BMD"]),"DISTFOREARM_BMD"]
forarm[1:5,]

 pheno<-ann   
dim(pheno)
pheno[1:2,]

posns<-match(pheno[,"PATIENT"],forarm[,"PATIENT"])
missing<-is.na(posns)
sum(missing)
sum(!missing)
RTOT_BMD_USED<-forarm[posns,"RTOT_BMD"]
pheno<-cbind(pheno,RTOT_BMD_USED)
dim(pheno)
pheno[1:2,]
ann<-pheno

    
##                PATIENT      VISIT COHORT CENTRE GENDER
## 107   AOGC-02-0105  1/08/2008      1      2      2
## 246   AOGC-02-0250 12/11/2008      1      2      2
## 399   AOGC-02-0405 26/05/2009      1      2      2
## 943   AOGC-08-0170 18/08/2009      1      8      2
## 963   AOGC-08-0194  1/09/2009      1      8      2
## 1035  AOGC-08-0279 30/09/2009      1      8      2
## 1039  AOGC-08-0284  2/10/2009      1      8      2
## 1040  AOGC-08-0285  1/10/2009      1      8      2
## 1046  AOGC-08-0291  6/10/2009      1      8      2
## 3638  AOGC-14-0064 25/03/2009      1     14      2
## 3641  AOGC-14-0081 25/03/2009      1     14      2
## 3642  AOGC-14-0093 25/03/2009      1     14      2
## 3646  AOGC-14-0154 25/03/2009      1     14      2
## 3760  AOGC-14-1181 25/03/2009      1     14      2
## 4113  AOGC-14-4035 25/03/2009      1     14      2
## 9013  AOGC-15-0016 20/10/2009      1     15      2
## 1018  AOGC-08-0258 29/09/2009      1      8      2
## 880   AOGC-08-0083  9/07/2009      1      8      2
## 15158  S01-F08-P01 13/01/2004      4     21      2
## 15164  S01-F14-P01 22/12/2008      4     21      2
## 15218  S02-F17-P01  4/10/2006      4     22      2
## 15243  S04-F14-P01 13/09/2007      4     24      2
## 15511  S09-F83-P01 12/11/2008      4     29      2
## 15523  S11-F02-P01 20/06/2005      4     31      2
## 15555  S13-F08-P01 12/08/2008X1      4     33      2
## 15610  S14-F15-P01 14/11/2003      4     34      2

## contaminated<- a.ann[!(a.ann[,"PATIENT"] %in% ann2[,"PATIENT"]),"PATIENT"]   
##   write.table(cbind(contaminated,contaminated),file="/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/contaminated_AOGC_SEQ_samples.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)  
colnames(ann)
#target<-"TOT_HIP_GCM"
target<-"FN_GCM2"
target<-"LS_GCM"
target<-"PATIENT_RAD"
target<-"WEIGHT"
target<-"AGE_SCAN"
target<-"HEIGHT"
target<-"BMI"
target<-"RTOT_BMD_USED"

fam<-a.fam

    
posns<-match(fam[,1],ann[,"PATIENT"])
missing<-is.na(posns)
sum(missing)

a.ann<-ann[posns[!missing],]
wanted<-!is.na(a.ann[,"PCA1"]) & !is.na(a.ann[,target]) & !(a.ann[,"PATIENT"] %in% contaminated[,1])
sum(wanted)


#ann[,"PATIENT_FN"]<-as.numeric(ann[,"predicted_FN"])+ as.numeric(ann[,"value.resid_FN"])


sum(!is.na(a.ann[wanted,target]))
mean(a.ann[wanted,target],na.rm=TRUE)
sd(a.ann[wanted,target],na.rm=TRUE)
median(a.ann[wanted,target],na.rm=TRUE)
min(a.ann[wanted,target],na.rm=TRUE)
max(a.ann[wanted,target],na.rm=TRUE)

a.ann[,"BMI"]<-as.numeric(a.ann[,"WEIGHT"])/((as.numeric(a.ann[,"HEIGHT"])/100)^2)

###########################
 # check the review of remainng snps



    









    ###########################
  RESTART WITH NEW SNPS
  RESTART WITH NEW SNPS
  RESTART WITH NEW SNPS
  RESTART WITH NEW SNPS
    RESTART WITH NEW SNPS

#######################################################################
######################################################################
##### reevoker and clecked a lot of new clustre need to recompile these and
##### RERUN THE SINGLE POINTS TESTS
####  REDO THE MEtAL ANaNLYSIS AND REGENERALE THE ANNOTATION /RESULT FILE
    ########
# ONES TO CHECK ARE WRITTEN TO FILES lIKE   /media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/meta.snp.test_chr23.txt
# THE CHECKED RESULTS ARE NOW IN /media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP list re-checked ED chr1-19/meta.snp.test_chr23.txt.scores
# scores file will exist but if none changed bed will not be present


    #
########################
######################
#gene list for last pass
#### must up to data data:

    
ann.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/annoatation_short_exome_chip.txt"
ann.short<-read.delim(ann.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
ann.short[1:5,]

data<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/final.meta.analysis.single.point.significant.reorg.evoker.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

#### THIS bit takes SNPS, SCORE and LIST files to get new bed/bim/ fam fo list those types

    

    

#------------------------------------------------------------------
# loccation of rescored files    
ROOT.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/2014-06-23_SNPs_HWE.gt.10-13andP.lt.0.001" # location of recorces beds

##IMPORTANT  The bed/list/score file must contain the chromosome name prefix chromsome sufix
## names of the list/score/bed files must be the chromosome
### set up the list names
check.flag<-"Mhairi.ED.Final"
list.file.prefix<-"hwe_andp."
list.file.sufix<-".txt"
### set up the score file names    
score.file.prefix<-list.file.prefix
score.file.sufix<-paste(list.file.sufix,".scores",sep="")
# set up to get the bed.file.names
bed.file.prefix<-"hwe_andp."
bed.file.sufix<-".bed"    
#---------------------------------------------------------------------------------    

#------------------------------------------------------------------
# loccation of rescored files    
ROOT.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/2014-06-23_MhairiRescore_InteresingSNPs" # location of recorces beds

##IMPORTANT  The bed/list/score file must contain the chromosome name prefix chromsome sufix
## names of the list/score/bed files must be the chromosome
### set up the list names

check.flag<-"Mhairi.ED"
list.file.prefix<-"Mhairi.redo."
list.file.sufix<-".txt"
### set up the score file names    
score.file.prefix<-list.file.prefix
score.file.sufix<-paste(list.file.sufix,".scores",sep="")
# set up to get the bed.file.names
bed.file.prefix<-"AOGC_gentrain."
bed.file.sufix<-".bed"    
#---------------------------------------------------------------------------------

    
    
#------------------------------------------------------------------
# loccation of rescored files    
ROOT.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP.list.re-checked.ED.chr1-19" # location of recorces beds

##IMPORTANT  The bed/list/score file must contain the chromosome name prefix chromsome sufix
## names of the list/score/bed files must be the chromosome
### set up the list names
check.flag<-"ED.modified"
    
list.file.prefix<-"meta.snp.test_"
list.file.sufix<-".txt"
### set up the score file names    
score.file.prefix<-list.file.prefix
score.file.sufix<-paste(list.file.sufix,".scores",sep="")
# set up to get the bed.file.names
bed.file.prefix<-"meta.snp.test_"
bed.file.sufix<-".bed"    
#---------------------------------------------------------------------------------    

    

############ location of fam and bim files to use GENERIC
the.bim.root<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Cluster_vis_clean/Cluster_viz/" # .chr.bim added later
the.fam<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Cluster_vis_clean/Cluster_viz/chr1/AOGC_gentrain.fam"


    

    
    
    
#core.lists<-c("meta.snp.test_chr")

    
restrict.analysis<-FALSE # restrict.to.projects
ignore.dirs<-c()
 ## numof of lists processing


    
 # need the score files
## core.beds.roots<-c("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/cluster_viz/AOGC_genmeta.snp.test_chr3_rescoredED.1.bed") # just used gentarian to fix
## the.bim.root<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP.list.re-checked.ED.chr1-19/AOGC_gentrain" # .chr.bim added later
## the.fam<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP.list.re-checked.ED.chr1-19/AOGC_opticall" # same fam no matter what the chr

## core.beds.roots<-c("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/cluster_viz/AOGC_genmeta.snp.test_chr3_rescoredED.1.bed") # just used gentarian to fix
## the.bim.root<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP.list.re-checked.ED.chr1-19/AOGC_gentrain" # .chr.bim added later
## the.fam<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP.list.re-checked.ED.chr1-19/AOGC_opticall"

    
   


## projects<-list.dirs(ROOT.dir)
## projects<-projects[!(projects %in% ignore.dirs)]
## projects<-projects[!(projects %in% ROOT.dir)] ## remove the root
projects<-(ROOT.dir)
##

all.data<-{}
count<-0


    
print(projects)
targets.file.ori<-""
all.beds<-{}
all.scores<-{}
all.snp.lists<-{}
#######
    ip<-1
    projects # this is directory where looking Will get all the chromomes there
for(ip in 1:length(projects)){
print(projects[ip])
on.target.stats<-{}  # this will be rebuilt from the individual

#.chr<-basename(projects[ip])
#the.dir<-dirname(projects[ip])

setwd(projects[ip])
files<-list.files(getwd())
files
snp.lists<-files[grepl(paste("^",list.file.prefix,sep=""),files) & grepl(paste(list.file.sufix,"$",sep=""),files)]
names.snp.lists<-gsub(paste(list.file.sufix,"$",sep=""),"",snp.lists)
names.snp.lists<-gsub(paste("^",list.file.prefix,sep=""),"",names.snp.lists)
names(snp.lists)<-names.snp.lists
snp.lists

## scores<-files[grepl("^meta.snp.test_chr",files) & grepl(".txt.scores$",files)] # files[grep("snp.list.scores$",files)]
## names(scores)<-gsub(".txt.scores$","",scores)

## scores<-files[grepl(paste("^",score.file.prefix,sep=""),files) & grepl(paste(score.file.sufix,"$",sep=""),files)]
## names(snp.lists)<-gsub(".txt$","",snp.lists)

scores<-files[grepl(paste("^",score.file.prefix,sep=""),files) & grepl(paste(score.file.sufix,"$",sep=""),files)]
names.scores<-gsub(paste(score.file.sufix,"$",sep=""),"",scores)
names.scores<-gsub(paste("^",score.file.prefix,sep=""),"",names.scores)
names(scores)<-names.scores
scores



beds<-files[grepl(paste("^",bed.file.prefix,sep=""),files) & grepl(paste(bed.file.sufix,"$",sep=""),files)]
names.beds<-gsub(paste(bed.file.sufix,"$",sep=""),"",beds)
names.beds<-gsub(paste("^",bed.file.prefix,sep=""),"",names.beds)
names(beds)<-names.beds
beds





if( (length(scores)!=length(snp.lists)) | (sum(!(names(scores) %in% names(snp.lists)))!=0)){ # missing a score
  all.names<-unique(c(names(scores), names(snp.lists)))
   print(paste("Chromosome",the.chr))
   print(paste("missing snp list:", all.names[!(all.names %in% names(snp.lists))])) 
   print(paste("missing scores:", all.names[!(all.names %in% names(scores))]))
}
  
if( (length(beds)!=length(snp.lists)) | (sum(!(names(beds) %in% names(snp.lists)))!=0)){ # missing a bed or incorrebly names
  all.names<-unique(c(names(beds), names(snp.lists)))
   print(paste("Chromosome",the.chr))
   print(paste("missing beds:", all.names[!(all.names %in% names(beds))]))
   print(paste("Additional beds:", names(beds)[!(names(beds) %in% core.lists)]))

}

## assign(paste("beds",the.chr,sep="."),value=beds)
## assign(paste("scores",the.chr,sep="."),value=scores)
## assign(paste("snp.lists",the.chr,sep="."),value=snp.lists)

## all.beds<-c(all.beds,paste("beds",the.chr,sep="."))
## all.scores<-c(all.scores,paste("scores",the.chr,sep="."))
## all.snp.lists<-c(all.snp.lists,paste("snp.lists",the.chr,sep="."))

} ## loop over ip which is a loop over directories at the moment


##############







## cd "/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/rechecked evoker data"
## The Bed save by evoker have ALL snps so need to extract the one you want!

    #######################
    ######################
    ########################


# At this stage have the bed 


    #########################
    #######################
    #######################



############ for each snp.list extract out just tghose genotypes
beds
scores
snp.lists
print(projects)

#
    ip<-1

maybe.snps<-{}
final.plink.files<-{}

for(ip in 1:length(projects)){
print(projects[ip])
setwd(projects[ip])


## the.chr<-basename(projects[ip])
## the.chr.int<-gsub("^chr","",the.chr)
## the.bim<-paste(the.bim.root,the.chr.int,"bim",sep=".")

## beds<-eval(as.name(paste("beds",the.chr,sep=".")))
## scores<-eval(as.name(paste("scores",the.chr,sep=".")))
## snp.lists<-eval(as.name(paste("snp.lists",the.chr,sep=".")))

beds
scores
snp.lists

scores<-scores[names(beds)]
snp.lists<-snp.lists[names(beds)]
sum(names(beds)!= names(scores))
sum(names(beds)!= names(snp.lists))


beds


## the.bim.root<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Cluster_vis_clean/Cluster_viz/" # .chr.bim added later
## the.fam<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Cluster_vis_clean/Cluster_viz/chr1/AOGC_gentrain.fam"



# ic<-1
maybe.snps.chr<-{}
plink.files.made<-{}
all.scores<-{}
ic<-1
for(ic in 1:length(beds)){
#  for(ic in 9:9){

#the.chr<-gsub("^meta.snp.test_","",names(beds)[ic])

    the.chr<-names(beds)[ic]

  the.bim<- paste(the.bim.root,the.chr,"/AOGC_gentrain.",gsub("^chr","",the.chr),".bim",sep="")
 # the.fam<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP.list.re-checked.ED.chr1-19/chr19/AOGC_gentrain.fam" # defined above
  the.fam
  
  the.list<-snp.lists[ic]
  num.snps.tested<-eval(system(paste("wc","-l", the.list,sep=" "),intern=TRUE)) ## if no snps plinkwill crash
  if(grepl("^0",num.snps.tested)){next}
  
  system(paste("plink","--bed",beds[ic],"--bim",the.bim,"--fam",the.fam,"--extract",snp.lists[ic],"--make-bed","--out",paste("extracted",names(beds)[ic],sep="."),"--noweb",sep=" "))

  plink.files.made<-c(plink.files.made,paste(paste("extracted",names(beds)[ic],sep="."),c("bed","bim","fam"),sep=".",collapse=" "))
  
  a.score<- read.table(scores[ic],header=FALSE,fill=TRUE,stringsAsFactors=FALSE)
  maybe.snps.chr<-c(maybe.snps.chr,a.score[a.score[,2]==0,1])
  if(ic==1){all.scores<-a.score}else{all.scores<-rbind(all.scores,a.score)}
} ## loop over ic
} ## loop over ip

plink.files.made
  write.table(plink.files.made[2:length(plink.files.made)],file=paste("plink.files.made","txt",sep="."),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

first<-unlist(strsplit(plink.files.made[1],split=" "))[1]
first<-gsub(".bed$","",first)
print(first)
paste("plink","--bfile",first,"--merge-list",paste("plink.files.made","txt",sep="."),"--make-bed","--out","all.lists","--noweb",sep=" ")

  system( paste("plink","--bfile",first,"--merge-list",paste("plink.files.made","txt",sep="."),"--make-bed","--out","all.lists","--noweb",sep=" ") )

ic<-1
all.scores<-{}
for(ic in 1:length(scores)){
  a.score<- read.table(scores[ic],header=FALSE,fill=TRUE,stringsAsFactors=FALSE)
  if(ic==1){all.scores<-a.score}else{all.scores<-rbind(all.scores,a.score)}
}
  dim(all.scores)

  all.scores[all.scores[,2]==-1,2]<-check.flag
  all.scores[all.scores[,2]==1,2]<-paste(check.flag,"modified")


  getwd()
  
write.table(all.scores,file="checked.snps.r2.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

## modified.snp.calls<-  all.scores[all.scores[,2]=="ckecked.modified.r2",1]

##   modified.snp.calls[1:5]
##   dim(all.scores) #1079
##   length( modified.snp.calls) # 71 of 1079
  

## write.table(modified.snp.calls,file="modified.snp.calls.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
  

## write.table(all.scores,file="checked.snps.r2.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)



  
## system( paste("plink","--bfile","all.lists","--extract","modified.snp.calls.txt","--make-bed","--out","all.lists_modified_calls","--noweb",sep=" ") )

  
#### now recode this set


### check the sample names are same not recoded
##### between realed.to.remove and PCA only one bad sample was one sneaking through  AOGC-02-0398 on June 24 decided to address this by making full exclude list
## related.to.remove<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.recode.txt",header=F,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
## related.to.remove[1:5,]
## bad.chips<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Cluster_vis_clean/Cluster_viz/excluded.samples.txt",header=F,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

## dim(bad.chips)
## dim(related.to.remove)
## bad.chips<-bad.chips[!(bad.chips[,1] %in% related.to.remove[,1]),] # bad chips not caught by related
## dim(bad.chips)

## related.to.remove<-rbind(related.to.remove,bad.chips)
## dups<-duplicated(related.to.remove[,1])
## sum(dups)
## related.to.remove[dups,]
## related.to.remove[related.to.remove[,1] %in% related.to.remove[dups,1],] 
## related.to.remove<-related.to.remove[!dups,]

## write.table(related.to.remove,file="/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.recode.txt",col.names=FALSE,row.names=FALSE,sep=",",quote=FALSE)

#### remove bad samples
  system( paste("plink","--bfile","all.lists","--remove","/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.txt","--make-bed","--out","all.lists_filtered","--noweb",sep=" ") )
#### now recode this set

################# recode GWAS to phenotype IDS:
  #working with
fam.template.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL_chr1.fam"
fam<-read.table(fam.template.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
exclude<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.recode.txt",header=F,fill=TRUE,stringsAsFactors=FALSE)

related.to.remove<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.txt",header=F,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
related.to.remove[1:5,]

annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)


# new checked genotypes in:"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP.list.re-checked.ED.chr1-19"
fam.new<-read.table("all.lists_filtered.fam",header=F,fill=TRUE,stringsAsFactors=FALSE)
recodes<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/pheno_gwas_id_missmatches.csv",,header=T,fill=TRUE,stringsAsFactors=FALSE)
recodes[1:5,]


dim(fam)
dim(fam.new)


sum(!(fam[,1] %in% fam.new[,1]))
fam[!(fam[,1] %in% fam.new[,1]),][1:10,]

recodes[1:5,]
recodes[,2] %in% fam[,1] # fam uses GWASID
recodes[,2] %in% exclude[,1] # exclude used GWASID
recodes[,2] %in% ann[,"PATIENT"]  # ann used GWASID

sum(fam.new[,1] %in% related.to.remove[,1]) ### cause were nor recoded 


recodes[,1] %in% fam.new[,1] # fam.new (original data  uses PhenoID
recodes[1:5,]
# fam.new (original data  uses PhenoID

### convert from 
posns<-match(fam.new[,1],recodes[,"PhenoID"])
missing<-is.na(posns)
sum(!missing)
fam.new[!missing,1]<-recodes[posns[!missing],"GWASID"]
fam.new[!missing,2]<-recodes[posns[!missing],"GWASID"]
sum(duplicated(fam.new[,1]))

## replace fam file with recoded
#setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP.list.re-checked.ED.chr1-19")
write.table(fam.new,file="all.lists_filtered.fam",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

system( paste("plink","--bfile","all.lists_filtered","--remove"," /media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.recode.txt","--make-bed","--out","all.lists_filtered-clean","--noweb",sep=" ") )

getwd()
#### samples are all recoded

   

#/media/UQCCG-Analysis/AOGC_exome_chip/Phenotypes/clean_AOGC_SEQ_samples.txt.csv

################ hardy for filtered
########################################################## single point for AOGC
## Writing pedigree information to [ all.lists_modified_calls_clean.fam ] 
## Writing map (extended format) information to [ all.lists_modified_calls_clean.bim ] 
## Writing genotype bitfile to [ all.lists_modified_calls_clean.bed ] 
## Using (default) SNP-major mode

## Analysis finished: Wed May 14 18:57:38 2014

######################################################### single point for exome chip

###### get these new SNPs onto the same strand as the orhin ones
# 1) extract from the origin set in "/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes"
# 2) change name
# 3) try to merge
# 4) flip
# 5) re-check
getwd()
fam.template.file<-"all.lists_filtered-clean.fam"

system("cut -f 2 all.lists_filtered-clean.bim  > extracted.snps")



system("plink --bfile /media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL --extract extracted.snps --make-bed --out extracted")
system("plink --bfile extracted --bmerge all.lists_filtered-clean.bed all.lists_filtered-clean.bim all.lists_filtered-clean.fam")

system("plink --bfile extracted --extract plink.missnp --freq --out extracted")
system("plink --bfile all.lists_filtered-clean --extract plink.missnp --freq --out redone")
### need this cause original data was in TOP strand
test<-system("grep Error plink.log | wc -l",intern=TRUE)
    test
if(as.numeric(test)>0){
  print("Flipping required")
  system( paste("plink","--bfile","all.lists_filtered-clean","--flip plink.missnp --make-bed","--out","all.lists_filtered-clean.f","--noweb",sep=" ") )
  fam.template.file<-"all.lists_filtered-clean.f.fam"
  system("plink --bfile extracted --bmerge all.lists_filtered-clean.f.bed all.lists_filtered-clean.f.bim all.lists_filtered-clean.f.fam")
  test<-system("grep Error plink.log | wc -l",intern=TRUE)
  if(as.numeric(test)==0){print("Flip fixed all issues")}else{print("ERROR inconistances remaon see plonk .log")}
}


## Below is just some testing to check now hits behave
## system( paste("plink","--bfile","hits","--remove","/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/IBD/related.to.remove.txt","--make-bed","--out","hits-clean","--noweb",sep=" ") )
## system( paste("plink","--bfile","hits-clean","--bmerge all.lists_filtered-clean.bed all.lists_filtered-clean.bim all.lists_filtered-clean.fam","--make-bed","--out","all-test","--noweb",sep=" ") )
## /media/UQCCG/UQCCG-Projects/AOGC_exome_chip/2014-06-23_MhairiRescore_InteresingSNPs/all-test-merge.missnp
## system( paste("plink","--bfile","all.lists_filtered-clean","--flip all-test-merge.missnp --make-bed","--out","all.lists_filtered-clean.f","--noweb",sep=" ") )
## system( paste("plink","--bfile","hits-clean","--bmerge all.lists_filtered-clean.f.bed all.lists_filtered-clean.f.bim all.lists_filtered-clean.f.fam","--make-bed","--out","all-test","--noweb",sep=" ") )


## system( paste("plink","--bfile","hits-clean","--bmerge all.lists_filtered-clean.bed all.lists_filtered-clean.bim all.lists_filtered-clean.fam","--make-bed","--out","all-test","--noweb",sep=" ") )

## fam.template.file<-"all-test.fam" #### USE to test just the current hits

############### original data
## fam.template.file<-"recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.fam"
## annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
## projects<-"recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL"
## ROOT.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes"

## ROOT.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/2014-06-23_MhairiRescore_InteresingSNPs"


setwd(ROOT.dir)
fam.template.file
# fam.template.file set above in cases where flipping data is required
fam<-read.table(fam.template.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
projects<-gsub(".fam$","",fam.template.file)

annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
## annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/test.ori"
## annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/test.ori.1.pca"
## annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/test.ori.2pca"
## annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/test.ori.no.age.sq"
## annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/test.ori.center.lin"
## annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/test.ori.gender"
#annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_SAMPLES_PHENOTYPES.txt"
#annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_SAMPLES_PHENOTYPES_Nov.1.2013_RESIDUALS.txt"
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
ann[1:3,]

## ann[1:100,c("BMD_AFFSTAT","AFFSTAT_IN_WORDS")]
## table(ann[,"AFFSTAT_IN_WORDS"])
## table(ann[,"BMD_AFFSTAT"])
## are.na<-is.na(ann[,"BMD_EFF_STD_HIP"])
## ann[,"BMD_EFF_STD_HIP"]
## hi<-ann[,"TOT_HIP_Z"] >= 1.5 & !is.na(ann[,"TOT_HIP_Z"]) & !is.na(ann[,"PCA1"])
## low<-ann[,"TOT_HIP_Z"] <= -1.5 & !is.na(ann[,"TOT_HIP_Z"]) & !is.na(ann[,"PCA1"])
## ## hi<-ann[,"BMD_EFF_STD_HIP"] >= 1.5 & !is.na(ann[,"BMD_EFF_STD_HIP"])
## ## low<-ann[,"BMD_EFF_STD_HIP"] <= -1.5 & !is.na(ann[,"BMD_EFF_STD_HIP"])

## sum(hi)
## sum(low)
## sum(hi & is.na(ann[,"BMD_AFFSTAT"]))
## sum(low & is.na(ann[,"BMD_AFFSTAT"]))
## miss.hi<-(hi & is.na(ann[,"BMD_AFFSTAT"]))
## miss.low<-(low & is.na(ann[,"BMD_AFFSTAT"]))
## sum(miss.hi)
## sum(miss.low)
## ann[miss.hi,"BMD_AFFSTAT"]<-1
## write.table(ann,file="/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS_2.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


    
#hist(ann[posns[!missing],traits[i]])

traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")



#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
names(traits)<-rep("logistic",times=length(traits))
names(traits)[1:3]<-rep("assoc",times=3)

#traits
traits %in% colnames(ann)


fam[1:5,]

i<-1
ip<-1
 # set on line 4600


for (i in 1:length(traits)){


fam[,6]<--9
posns<-match(fam[,1],ann[,"PATIENT"])
missing<-is.na(posns)
print(sum(!missing))

no.pheno<-ann[posns,traits[i]]==3 | ann[posns,traits[i]]==-9 | is.na(ann[posns,traits[i]]) | ann[posns,"GENDER"]!=2  # ann[posns,traits[i]]

#print(sum(no.pheno))
#table(ann[posns[!missing],traits[i]])

missing<-missing | no.pheno
#cbind(fam[no.pheno,],ann[posns[no.pheno],c("PATIENT",traits[i])])
print(sum(!missing))

if(names(traits)[i]=="logistic"){
# using.covars<-c("AGE_SCAN","AGEXAGE","PCA1","PCA2","PCA3","PCA4")
  using.covars<-c("AGE_SCAN","PCA1","PCA2","PCA3","PCA4")

 
  covars<- ann[posns,c("PATIENT","PATIENT", using.covars)]
 colnames(covars)[1:2]<-c("FID","IID")
 covars[1:5,]
  missing.covars<-is.na(covars[,using.covars])
   missing.covars<-apply(missing.covars,1,sum,na.rm=TRUE)
 missing.covars<-  missing.covars>0
missing<-missing | missing.covars
write.table(covars[!missing,],file=paste("covars",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
}


has.missing<-sum(missing)>0
has.missing
write.table(fam[missing,c(1,1)],file=paste("remove_from",traits[i],sep="."),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

#chk<-fam[missing,c(1,1)]

fam[!missing,6]<-ann[posns[!missing],traits[i]] ### put in the traints
the.fam<-paste(traits[i],"fam",sep=".")

write.table(fam,file=the.fam,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
dim(fam)


ip<-1
projects
for(ip in 1:length(projects)){
print(projects[ip])
print(traits[i])
#setwd(projects[ip])

the.bed<-paste(projects[ip],"bed",sep=".")
the.bim<-paste(projects[ip],"bim",sep=".")


if(names(traits)[i]=="assoc"){
  system( paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam,"--remove",paste("remove_from",traits[i],sep="."),"--assoc","--allow-no-sex","--out",paste(traits[i],projects[ip],sep="."),"--noweb",sep=" ") )
}

if(names(traits)[i]=="logistic"){
  system( paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam,"--remove",paste("remove_from",traits[i],sep="."),"--covar",paste("covars",traits[i],sep="."),"--covar-name",toString(using.covars),"--logistic --ci 0.95 ","--allow-no-sex","--out",paste(traits[i],projects[ip],sep="."),"--noweb",sep=" ") ) # --vif 500 
}


} # ip loop over chromosomes
} # i loop over traits



    ########### below lines ofr testing only
##   system( paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam,"--remove",paste("remove_from",traits[i],sep="."),"--covar",paste("covars",traits[i],sep="."),"--covar-name",toString(using.covars),"--extract hits.txt  --linear --hide-covar --ci 0.95 ","--allow-no-sex","--out",paste(traits[i],projects[ip],sep="."),"--noweb",sep=" ") )

   ## system( paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam,"--remove",paste("remove_from",traits[i],sep="."),"--covar",paste("covars",traits[i],sep="."),"--covar-name",toString(using.covars),"--extract hits.txt  --logistic --hide-covar --ci 0.95 ","--allow-no-sex","--out",paste(traits[i],projects[ip],sep="."),"--noweb",sep=" ") )

## system( paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam,"--remove",paste("remove_from",traits[i],sep="."),"--extract hits.txt --assoc","--allow-no-sex","--out", paste(traits[i],projects[ip],sep="."),"--noweb",sep=" ") )

###### below is a test showing we are robust

## ############## ORIGINAL
##  CHR                                                                   SNP         BP    NMISS       BETA         SE         R2        T            P 
##    1   UQDIexomechipWithAOGC_chr1:152081443:152081443:C:T-1_T_R_2110821056  152081443     6905     0.9065     0.2946    0.00137    3.077     0.002098 
##    1   UQDIexomechipWithAOGC_chr1:152193920:152193920:T:G-1_B_F_2098201016  152193920     6905     -1.704     0.4417   0.002152   -3.859     0.000115 
##    4     UQDIexomechipWithAOGC_chr4:31144670:31144670:G:A-1_B_R_2098208305   31144670     6904     0.2852      0.159   0.000466    1.794      0.07289 
##    7     UQDIexomechipWithAOGC_chr7:49952040:49952040:C:G-1_B_F_2110822923   49952040     6904     0.3206     0.1206   0.001023    2.658      0.00787 
##    9     UQDIexomechipWithAOGC_chr9:71145183:71145183:C:T-1_B_F_2098215432   71145183     6904     -1.038     0.2945   0.001796   -3.524    0.0004286 
##   10    UQDIexomechipWithAOGC_chr10:69677269:69677269:T:G-1_B_F_2110819177   69677269     6904     0.3429     0.3953   0.000109   0.8675       0.3857 
##   12  UQDIexomechipWithAOGC_chr12:120899502:120899502:A:C-1_B_R_2098221197  120899502     6904       1.71     0.3948   0.002712    4.332    1.495e-05 
##   14  UQDIexomechipWithAOGC_chr14:105419224:105419224:G:A-1_B_R_2098223279  105419224     6904   0.004469     0.0553  9.463e-07  0.08082       0.9356 
##   15    UQDIexomechipWithAOGC_chr15:98508942:98508942:G:A-1_B_R_2098251903   98508942     6904    -0.1095    0.04434  0.0008836   -2.471      0.01351 
##   16    UQDIexomechipWithAOGC_chr16:76555177:76555177:A:T-1_T_R_2098225575   76555177     6904    -0.2122    0.07514   0.001154   -2.824     0.004757


## ##original + gender (44 males) 
##  CHR                                                                   SNP         BP    NMISS       BETA         SE         R2        T            P 
##    1   UQDIexomechipWithAOGC_chr1:152081443:152081443:C:T-1_T_R_2110821056  152081443     6905     0.9068     0.2944   0.001372     3.08     0.002079 
##    1   UQDIexomechipWithAOGC_chr1:152193920:152193920:T:G-1_B_F_2098201016  152193920     6905     -1.705     0.4414   0.002156   -3.862    0.0001134 
##    4     UQDIexomechipWithAOGC_chr4:31144670:31144670:G:A-1_B_R_2098208305   31144670     6904     0.2851     0.1589  0.0004661    1.794      0.07286 
##    7     UQDIexomechipWithAOGC_chr7:49952040:49952040:C:G-1_B_F_2110822923   49952040     6904     0.3202     0.1206   0.001021    2.656     0.007917 
##    9     UQDIexomechipWithAOGC_chr9:71145183:71145183:C:T-1_B_F_2098215432   71145183     6904     -1.037     0.2943   0.001796   -3.524    0.0004278 
##   10    UQDIexomechipWithAOGC_chr10:69677269:69677269:T:G-1_B_F_2110819177   69677269     6904     0.3425     0.3951  0.0001089    0.867       0.3859 
##   12  UQDIexomechipWithAOGC_chr12:120899502:120899502:A:C-1_B_R_2098221197  120899502     6904      1.709     0.3946    0.00271     4.33     1.51e-05 
##   14  UQDIexomechipWithAOGC_chr14:105419224:105419224:G:A-1_B_R_2098223279  105419224     6904   0.004654    0.05527  1.027e-06   0.0842       0.9329 
##   15    UQDIexomechipWithAOGC_chr15:98508942:98508942:G:A-1_B_R_2098251903   98508942     6904    -0.1093    0.04431  0.0008803   -2.466      0.01369 
##   16    UQDIexomechipWithAOGC_chr16:76555177:76555177:A:T-1_T_R_2098225575   76555177     6904    -0.2122     0.0751   0.001155   -2.826     0.004733

## extreme.regression_aogc_hbm_testing_post_analysis.r
## #plink --linear no AGExAGE cause is colinear NOTE the  (sd(residulas)  is 0.13478 so 1/sd is 7.419 for HIP so these BETA * 7.419 gve us our numbers 
##  CHR                                                                   SNP         BP   A1       TEST    NMISS       BETA       SE      L95      U95         STAT            P 
##    1   UQDIexomechipWithAOGC_chr1:152081443:152081443:C:T-1_T_R_2110821056  152081443    T        ADD     6905     0.1216  0.03954  0.04406    0.199        3.074     0.002118
##    1   UQDIexomechipWithAOGC_chr1:152193920:152193920:T:G-1_B_F_2098201016  152193920    G        ADD     6905    -0.2297  0.05931   -0.346  -0.1135       -3.874    0.0001082
##    4     UQDIexomechipWithAOGC_chr4:31144670:31144670:G:A-1_B_R_2098208305   31144670    A        ADD     6904     0.0375  0.02135 -0.004343  0.07934        1.757      0.07904
##    7     UQDIexomechipWithAOGC_chr7:49952040:49952040:C:G-1_B_F_2110822923   49952040    G        ADD     6904    0.04246  0.01619  0.01073  0.07418        2.623     0.008736
##    9     UQDIexomechipWithAOGC_chr9:71145183:71145183:C:T-1_B_F_2098215432   71145183    T        ADD     6904    -0.1389  0.03951  -0.2164 -0.06152       -3.517    0.0004392
##   10    UQDIexomechipWithAOGC_chr10:69677269:69677269:T:G-1_B_F_2110819177   69677269    G        ADD     6904    0.04621  0.05303 -0.05774   0.1502       0.8712       0.3837
##   12  UQDIexomechipWithAOGC_chr12:120899502:120899502:A:C-1_B_R_2098221197  120899502    C        ADD     6904     0.2255  0.05297   0.1217   0.3293        4.258    2.092e-05
##   14  UQDIexomechipWithAOGC_chr14:105419224:105419224:G:A-1_B_R_2098223279  105419224    A        ADD     6904  0.0009415 0.007427 -0.01362   0.0155       0.1268       0.8991
##   15    UQDIexomechipWithAOGC_chr15:98508942:98508942:G:A-1_B_R_2098251903   98508942    A        ADD     6904   -0.01438  0.00595 -0.02604 -0.00272       -2.417      0.01567
##   16    UQDIexomechipWithAOGC_chr16:76555177:76555177:A:T-1_T_R_2098225575   76555177    T        ADD     6904   -0.02743  0.01009  -0.0472 -0.00766       -2.719     0.006556

## #### originla but center as linaer in lme
##  CHR                                                                   SNP         BP    NMISS       BETA         SE         R2        T            P 
##    1   UQDIexomechipWithAOGC_chr1:152081443:152081443:C:T-1_T_R_2110821056  152081443     6905     0.9063     0.2943   0.001372    3.079     0.002083 
##    1   UQDIexomechipWithAOGC_chr1:152193920:152193920:T:G-1_B_F_2098201016  152193920     6905     -1.702     0.4412    0.00215   -3.857     0.000116 
##    4     UQDIexomechipWithAOGC_chr4:31144670:31144670:G:A-1_B_R_2098208305   31144670     6904     0.2849     0.1589  0.0004659    1.794      0.07292 
##    7     UQDIexomechipWithAOGC_chr7:49952040:49952040:C:G-1_B_F_2110822923   49952040     6904     0.3202     0.1205   0.001022    2.657       0.0079 
##    9     UQDIexomechipWithAOGC_chr9:71145183:71145183:C:T-1_B_F_2098215432   71145183     6904     -1.039     0.2942   0.001803   -3.531    0.0004175 
##   10    UQDIexomechipWithAOGC_chr10:69677269:69677269:T:G-1_B_F_2110819177   69677269     6904     0.3418     0.3949  0.0001086   0.8656       0.3867 
##   12  UQDIexomechipWithAOGC_chr12:120899502:120899502:A:C-1_B_R_2098221197  120899502     6904      1.709     0.3944   0.002714    4.334    1.486e-05 
##   14  UQDIexomechipWithAOGC_chr14:105419224:105419224:G:A-1_B_R_2098223279  105419224     6904   0.003615    0.05525  6.204e-07  0.06544       0.9478 
##   15    UQDIexomechipWithAOGC_chr15:98508942:98508942:G:A-1_B_R_2098251903   98508942     6904    -0.1085    0.04429  0.0008693    -2.45      0.01429 
##   16    UQDIexomechipWithAOGC_chr16:76555177:76555177:A:T-1_T_R_2098225575   76555177     6904    -0.2133    0.07507   0.001168   -2.841     0.004514

## ####- original same model based on LME
##    1   UQDIexomechipWithAOGC_chr1:152081443:152081443:C:T-1_T_R_2110821056  152081443     6905     0.9066     0.2943   0.001373    3.081     0.002073 
##    1   UQDIexomechipWithAOGC_chr1:152193920:152193920:T:G-1_B_F_2098201016  152193920     6905     -1.703     0.4412   0.002154    -3.86    0.0001142 
##    4     UQDIexomechipWithAOGC_chr4:31144670:31144670:G:A-1_B_R_2098208305   31144670     6904     0.2853     0.1589  0.0004672    1.796       0.0725 
##    7     UQDIexomechipWithAOGC_chr7:49952040:49952040:C:G-1_B_F_2110822923   49952040     6904       0.32     0.1205   0.001021    2.656     0.007932 
##    9     UQDIexomechipWithAOGC_chr9:71145183:71145183:C:T-1_B_F_2098215432   71145183     6904     -1.037     0.2942   0.001799   -3.527    0.0004237 
##   10    UQDIexomechipWithAOGC_chr10:69677269:69677269:T:G-1_B_F_2110819177   69677269     6904     0.3434     0.3949  0.0001096   0.8697       0.3845 
##   12  UQDIexomechipWithAOGC_chr12:120899502:120899502:A:C-1_B_R_2098221197  120899502     6904      1.708     0.3944    0.00271    4.331    1.505e-05 
##   14  UQDIexomechipWithAOGC_chr14:105419224:105419224:G:A-1_B_R_2098223279  105419224     6904   0.004235    0.05524  8.513e-07  0.07665       0.9389 
##   15    UQDIexomechipWithAOGC_chr15:98508942:98508942:G:A-1_B_R_2098251903   98508942     6904    -0.1092    0.04429  0.0008792   -2.464      0.01375 
##   16    UQDIexomechipWithAOGC_chr16:76555177:76555177:A:T-1_T_R_2098225575   76555177     6904     -0.212    0.07506   0.001154   -2.824     0.004751 
## ####################### redo hardy for extra

## #### lme above but with 2PCA
##  CHR                                                                   SNP         BP    NMISS       BETA         SE         R2        T            P 
##    1   UQDIexomechipWithAOGC_chr1:152081443:152081443:C:T-1_T_R_2110821056  152081443     6905     0.9069     0.2942   0.001374    3.082     0.002062 
##    1   UQDIexomechipWithAOGC_chr1:152193920:152193920:T:G-1_B_F_2098201016  152193920     6905     -1.709     0.4411    0.00217   -3.875    0.0001076 
##    4     UQDIexomechipWithAOGC_chr4:31144670:31144670:G:A-1_B_R_2098208305   31144670     6904     0.2807     0.1588  0.0004524    1.767      0.07721 
##    7     UQDIexomechipWithAOGC_chr7:49952040:49952040:C:G-1_B_F_2110822923   49952040     6904     0.3183     0.1205    0.00101    2.642     0.008265 
##    9     UQDIexomechipWithAOGC_chr9:71145183:71145183:C:T-1_B_F_2098215432   71145183     6904     -1.039     0.2941   0.001807   -3.534    0.0004114 
##   10    UQDIexomechipWithAOGC_chr10:69677269:69677269:T:G-1_B_F_2110819177   69677269     6904      0.336     0.3948  0.0001049    0.851       0.3948 
##   12  UQDIexomechipWithAOGC_chr12:120899502:120899502:A:C-1_B_R_2098221197  120899502     6904      1.714     0.3943    0.00273    4.346    1.404e-05 
##   14  UQDIexomechipWithAOGC_chr14:105419224:105419224:G:A-1_B_R_2098223279  105419224     6904   0.003646    0.05523  6.313e-07  0.06601       0.9474 
##   15    UQDIexomechipWithAOGC_chr15:98508942:98508942:G:A-1_B_R_2098251903   98508942     6904    -0.1082    0.04428   0.000865   -2.444      0.01453 
##   16    UQDIexomechipWithAOGC_chr16:76555177:76555177:A:T-1_T_R_2098225575   76555177     6904    -0.2116    0.07505   0.001151    -2.82     0.004813

## #### lme above but with 1PCA
##  CHR                                                                   SNP         BP    NMISS       BETA         SE         R2        T            P 
##    1   UQDIexomechipWithAOGC_chr1:152081443:152081443:C:T-1_T_R_2110821056  152081443     6905     0.9068     0.2942   0.001374    3.082     0.002064 
##    1   UQDIexomechipWithAOGC_chr1:152193920:152193920:T:G-1_B_F_2098201016  152193920     6905     -1.709     0.4411    0.00217   -3.875    0.0001077 
##    4     UQDIexomechipWithAOGC_chr4:31144670:31144670:G:A-1_B_R_2098208305   31144670     6904     0.2807     0.1588  0.0004522    1.767      0.07725 
##    7     UQDIexomechipWithAOGC_chr7:49952040:49952040:C:G-1_B_F_2110822923   49952040     6904     0.3184     0.1205   0.001011    2.643     0.008243 
##    9     UQDIexomechipWithAOGC_chr9:71145183:71145183:C:T-1_B_F_2098215432   71145183     6904     -1.039     0.2941   0.001807   -3.534    0.0004116 
##   10    UQDIexomechipWithAOGC_chr10:69677269:69677269:T:G-1_B_F_2110819177   69677269     6904     0.3358     0.3948  0.0001048   0.8506        0.395 
##   12  UQDIexomechipWithAOGC_chr12:120899502:120899502:A:C-1_B_R_2098221197  120899502     6904      1.714     0.3943    0.00273    4.347    1.401e-05 
##   14  UQDIexomechipWithAOGC_chr14:105419224:105419224:G:A-1_B_R_2098223279  105419224     6904   0.003631    0.05523  6.261e-07  0.06574       0.9476 
##   15    UQDIexomechipWithAOGC_chr15:98508942:98508942:G:A-1_B_R_2098251903   98508942     6904    -0.1083    0.04428  0.0008657   -2.445      0.01449 
##   16    UQDIexomechipWithAOGC_chr16:76555177:76555177:A:T-1_T_R_2098225575   76555177     6904    -0.2116    0.07505    0.00115   -2.819     0.004828 

## #### no age.sq
##  CHR                                                                   SNP         BP    NMISS       BETA         SE         R2        T            P 
##    1   UQDIexomechipWithAOGC_chr1:152081443:152081443:C:T-1_T_R_2110821056  152081443     6905     0.9063     0.2943   0.001372    3.079     0.002083 
##    1   UQDIexomechipWithAOGC_chr1:152193920:152193920:T:G-1_B_F_2098201016  152193920     6905     -1.702     0.4412    0.00215   -3.857     0.000116 
##    4     UQDIexomechipWithAOGC_chr4:31144670:31144670:G:A-1_B_R_2098208305   31144670     6904     0.2849     0.1589  0.0004659    1.794      0.07292 
##    7     UQDIexomechipWithAOGC_chr7:49952040:49952040:C:G-1_B_F_2110822923   49952040     6904     0.3202     0.1205   0.001022    2.657       0.0079 
##    9     UQDIexomechipWithAOGC_chr9:71145183:71145183:C:T-1_B_F_2098215432   71145183     6904     -1.039     0.2942   0.001803   -3.531    0.0004175 
##   10    UQDIexomechipWithAOGC_chr10:69677269:69677269:T:G-1_B_F_2110819177   69677269     6904     0.3418     0.3949  0.0001086   0.8656       0.3867 
##   12  UQDIexomechipWithAOGC_chr12:120899502:120899502:A:C-1_B_R_2098221197  120899502     6904      1.709     0.3944   0.002714    4.334    1.486e-05 
##   14  UQDIexomechipWithAOGC_chr14:105419224:105419224:G:A-1_B_R_2098223279  105419224     6904   0.003615    0.05525  6.204e-07  0.06544       0.9478 
##   15    UQDIexomechipWithAOGC_chr15:98508942:98508942:G:A-1_B_R_2098251903   98508942     6904    -0.1085    0.04429  0.0008693    -2.45      0.01429 
##   16    UQDIexomechipWithAOGC_chr16:76555177:76555177:A:T-1_T_R_2098225575   76555177     6904    -0.2133    0.07507   0.001168   -2.841     0.004514 




########### Get Hardy for new tests

fam.template.file # set above
projects # set above

traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")

#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
names(traits)<-rep("logistic",times=length(traits))
names(traits)[1:3]<-rep("assoc",times=3)
traits
traits %in% colnames(ann)


fam[1:5,]


i<-1
for (i in 1:length(traits)){

  ### alredy got remove file from above
## write.table(fam[missing,c(1,1)],file=paste("remove_from",traits[i],sep="."),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

## fam[!missing,6]<-ann[posns[!missing],traits[i]]
the.fam<-paste(traits[i],"fam",sep=".")
### already wrote fam
## write.table(fam,file=the.fam,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)


ip<-1
projects
for(ip in 1:length(projects)){
print(projects[ip])
#setwd(projects[ip])

the.bed<-paste(projects[ip],"bed",sep=".")
the.bim<-paste(projects[ip],"bim",sep=".")

  system( paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam,"--remove",paste("remove_from",traits[i],sep="."),"--hardy","--allow-no-sex","--out",paste(traits[i],projects[ip],sep="."),"--noweb",sep=" ") )

} # ip loop over projects / chromosomes
} # i loop over traits

    

################################ prepare for META analysis  $$$ of extra checks $$$
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
################################ prepare for META analysis
#### all these in
#output.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/filtered"
output.dir<-ROOT.dir  

annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)



the.seq.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/AOGC_vcf_to_plink" # dirname(fam.template.file)
files<-dir(the.seq.dir)  
projects<-gsub(".bed","",files[grepl(".bed$",files)])  # files[grepl(".bed$",files)]


## project.prefix<-"AOGC_vcf_final_chr" ## old run moved to /media/UQCCG-Analysis/AOGC_exome_chip/meta_analysis_single_point/non-filtered-seq
## output.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/notFiltered"

## project.prefix<-"AOGC_to_vcf_filtered_" ## outpt used same name as above
## output.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/filtered"

project.prefix<-"AOGC_to_vcf_filtered_" ## outpt used same name as above
project.prefix.old<-"AOGC_vcf_final_chr"

    
projects.old<-projects[grepl(paste("^",project.prefix.old,sep=""),projects)]
projects<-projects[grepl(paste("^",project.prefix,sep=""),projects)] #

projects
projects.old

the.chr.old<-gsub(project.prefix.old,"",projects.old)
the.chr<-gsub(project.prefix,"",projects)

sum(the.chr.old !=the.chr) # should be ZERO else out of order




########### projects, pojects.old ref.cohort all in the same order
projects

traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")

names(traits)<-rep("logistic",times=length(traits))
names(traits)[1:3]<-rep("assoc",times=3)
traits
traits %in% colnames(ann)


#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
names(traits)<-rep(".assoc.logistic",times=length(traits))
names(traits)[1:3]<-rep(".qassoc",times=3)
traits

# Corresdonding in same order for exome chip
logistic.project.dir<-ROOT.dir
if(!grepl("/$",logistic.project.dir)){logistic.project.dir<-paste(logistic.project.dir,"/",sep="")}
logistic.project<-gsub(".fam$","",fam.template.file)
logistic.files<-paste(logistic.project.dir,traits,".",logistic.project,names(traits),sep="")
logistic.files
## logistic.files<-c(
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/Hip/totalhipBMDassocinexcel.txt",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/LS/LSBMDassocinexcel.txt",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/FN/FemNeckBMDassocinexcel.txt",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/EverFx/everFx_final_analysis.logistic.covar.assoc.logistic",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/everFx18nottrivia/everFx18nottrivia.LOGISTIC.COVAR.assoc.logistic",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/AllFr_50_OP/AllFr_50_OP.logistic.covar.assoc.logistic",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/AllVertFx_OP/AllVertFX_OP.logistic.covar.assoc.logistic",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/AllNonVertFX_50_OP/AllNonVertFX_50_OP.logistic.covar.assoc.logistic",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/ForearmFrLoTrauma/ForearmFrLoTrauma.logistic.covar.assoc.logistic",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/AllHipFr_50_OP/AllHipFr_50_OP.logistic.covar.assoc.logistic",
## "/media/UQCCG-Analysis/AOGC_exome_chip/for Paul for logistic regression analysis on bigger computer/EMcC_BP_GRP/EMcC_BP_GRP.logistic.covar.assoc.logistic"
##                   )

#annotation_short_exome_chip.txt

ann.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/annoatation_short_exome_chip.txt"
ann.short<-read.delim(ann.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)


    
########################## get processed bim files that are all on the forward strand in annovar format"
code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")

######



##############################these bim files have the minor allele as A1

## fam.template.file
## # fam.template.file set above in cases where flipping data is required
## fam<-read.table(fam.template.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
## projects<-gsub(".fam$","",fam.template.file)


logistic.project.dir<-ROOT.dir
setwd(logistic.project.dir)
bim.file<-gsub(".fam$",".bim",fam.template.file)
bim<-read.table(bim.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
colnames(bim)<-c("chr","SNP","cm","start","ALT","REF")
bim.chip<-bim
bim.chip.ori<-bim.chip
## dim(bim.chip)
## dim(bim.chip1)
## diff<-(bim.chip[,"ALT"]!=bim.chip1[,"A1"])
## bim.chip[1:5,]
## bim.chip1[1:5,]
## sum(diff)
## cbind(bim.chip,bim.chip1)[diff,][1:10,]


bim.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/AOGC_vcf_to_plink/AOGC_to_vcf_filtered_chrALL.bim"
bim<-read.table(bim.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
colnames(bim)<-c("chr","SNP","cm","start","ALT","REF") # A1 is labelled ALT here
bim.seq.ori<-bim
bim.seq<-bim.seq.ori


    
bim.chip[1:5,]
bim.seq[1:5,]
dim(bim.seq)
  dim(bim.chip)
#########

# save(list=c("bim.seq","bim.chip"),file="/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_A1_A2.RData")
# save(list=c("bim.seq","bim.chip"),file="/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_from_filtered_VCF.RData")
#save(list=c("bim.seq","bim.chip"),file="/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_fromVCF.RData")
##      chr  SNP                                           junk POS         REF ALT
## [1,] "12" "12:9009999:9009999:NA:-:indel"               "0"  "9009999"   NA  "-"
## [2,] "14" "14:102547769:102547769:NA:-:indel:102547769" "0"  "102547769" NA  "-"
## [3,] "16" "16:70011947:70011947:NA:-:indel:70011947"    "0"  "70011947"  NA  "-"
## [4,] "17" "17:38153714:38153714:NA:-:indel:38153714"    "0"  "38153714"  NA  "-"
## [5,] "20" "20:62606067:62606067:NA:-:indel"             "0"  "62606067"  NA  "-"
## [6,] "2"  "2:39210719:39210719:NA:-:indel:39210719"     "0"  "39210719"  NA  "-"

## bim.chip<-read.delim(bim.chip.file,header=F,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
i<-1
projects

    
for (i in 1:length(traits)){

#for (i in 1:12){
  
#load("/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_fromVCF.RData")
load("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_A1_A2.RData")
print(logistic.files[i])
logistic<-read.table(logistic.files[i],header=T,fill=TRUE,stringsAsFactors=FALSE)

colnames(logistic)[colnames(logistic)=="BP"]<-"POS"

if("TEST" %in% colnames(logistic)){
wanted<-grepl("ADD",logistic[,"TEST"])
logistic<-logistic[wanted,]
}

bim.chip<-bim.chip.ori
dim(logistic)
dim(bim.chip)


posns<-match(logistic[,"SNP"],bim.chip[,"SNP"])
missing<-is.na(posns)
sum(missing)
bim.chip<-bim.chip[posns,]

#core.ann<-c("chr","start","end","REF","ALT")
core.ann<-c("chr","start","REF","ALT")
key.logistic<-build.key(bim.chip,core.ann)
dups<-duplicated(key.logistic)
sum(dups) #804
## dups<-unique(key.logistic[dups])
## dups<-key.logistic %in% dups
## logistic[dups,][1:50,] #### dups all have the same valuse so can remove any one 
bim.chip<-bim.chip[!dups,]
logistic<-logistic[!dups,]
key.logistic<-key.logistic[!dups]


rownames(logistic)<-key.logistic
TYPE<-rep("GENO",times=length(key.logistic))
logistic<-cbind(logistic,TYPE,bim.chip[, core.ann],stringsAsFactors=FALSE)

## diff<-(logistic[,"A1"]!=bim.chip[,"ALT"])
## sum(diff)
## logistic[diff,][1:20,]

## ok<-logistic[,"CHR"]==2 & ( as.numeric(logistic[,"start"]) > 119000000 &  as.numeric(logistic[,"start"]) < 119729378)
## sum(ok)
## out<-logistic[ok,]
## posns<-match(out[,"SNP"],ann.short[,"SNP"])
## missing<-is.na(posns)
## sum(missing)
## out<-cbind(out,ann.short[posns,])
## out
## write.table(out,file="AOGC-LS-replication.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

## test<-key.logistic %in% c("chr10:101117384:101117384:G:A", "chr10:102307870:102307870:A:G")
## logistic[test,]
## bim.chip[test,]

## dim(test)

ip<-1 # i<-5
logistic.seq<-{}
bim.seq<-bim.seq.ori

###################read in all the sequence data NEED TO DO IN LOOP CASE ONE FOR EACH TEST!
for(ip in 1:length(projects)){

setwd(the.seq.dir)
  
if(i %in% 1:3){ # use rescalled data
#  NewBETA
logistic.file.seq<-paste("NewBETA",traits[i],projects[ip],sep=".")
print(logistic.file.seq)
setwd(the.seq.dir)
## a.logistic.seq<-read.table(logistic.file.seq,header=T,fill=TRUE,stringsAsFactors=FALSE) # old<-a.logistic.seq

column.labels<-read.table(logistic.file.seq,header=F,nrows=1,stringsAsFactors=FALSE)
num.vars<-dim(column.labels)[2]
a.logistic.seq<-scan(logistic.file.seq,what=character(num.vars),skip=1,fill=TRUE)
num.lines<-length(a.logistic.seq)/(num.vars)
dim(a.logistic.seq)<-c(num.vars,num.lines)
a.logistic.seq<-t(a.logistic.seq)
colnames(a.logistic.seq)<-column.labels


}else{
logistic.file.seq<- paste(traits[i],projects[ip],sep=".")
logistic.file.seq<-paste(logistic.file.seq,names(traits)[i],sep="") #
print(logistic.file.seq)
setwd(the.seq.dir)
#a.logistic.seq<-read.table(logistic.file.seq,header=T,fill=TRUE,stringsAsFactors=FALSE)

column.labels<-read.table(logistic.file.seq,header=F,nrows=1,fill=TRUE,stringsAsFactors=FALSE)
num.vars<-dim(column.labels)[2]
a.logistic.seq<-scan(logistic.file.seq,what=character(num.vars),skip=1,fill=TRUE)
num.lines<-length(a.logistic.seq)/(num.vars)
dim(a.logistic.seq)<-c(num.vars,num.lines)
a.logistic.seq<-t(a.logistic.seq)
colnames(a.logistic.seq)<-column.labels
}


colnames(a.logistic.seq)[colnames(a.logistic.seq)=="BP"]<-"POS"

if("TEST" %in% colnames(a.logistic.seq)){
wanted<-grepl("ADD",a.logistic.seq[,"TEST"])
a.logistic.seq<-a.logistic.seq[wanted,]
}



##############


if(ip==1){logistic.seq<-a.logistic.seq}else{logistic.seq<-rbind(logistic.seq,a.logistic.seq)}


}## loop over ip
#setwd(projects[ip])
#hold<-logistic.seq # logistic.seq<-hold
dim(logistic.seq)
dim(bim.seq)
#colnames(logistic.seq)[colnames(logistic.seq)=="BP"]<-"POS"

logistic.seq[1:5,]
logistic[1:5,]
bim.seq[1:5,]

posns<-match(logistic.seq[,"SNP"],bim.seq[,"SNP"]) # for annotation
missing<-is.na(posns)
sum(missing)
logistic.seq<-logistic.seq[!missing,]
bim.seq<-bim.seq[posns[!missing],]



dim(logistic.seq)
dim(bim.seq)
#logistic.seq[missing,][1:10,]
#bim.seq[1:5,]
################ remove flat genotypes  ### actaully let use there and reject the others using dup below and flat wlays first
## is.flat<-grepl(":flat$",bim.seq[,"SNP"]) # is.flat<-grepl(":flat$",logistic.seq[,"SNP"])
## sum(is.flat)
## ## bim.seq[is.flat,"SNP"][1:15]
## bim.seq<-bim.seq[!is.flat,]
## logistic.seq<-logistic.seq[!is.flat,]

core.ann<-c("chr","start","REF","ALT")

key.logistic.seq<-build.key(bim.seq,core.ann)
dups<-duplicated(key.logistic.seq)
sum(dups)
if(sum(dups)>0){
  
## dups<-unique(key.logistic.seq[dups])
## dups<-key.logistic.seq %in% dups
#logistic.seq[dups,][1:50,] #### dups all all unwound element not caught with end 
bim.seq<-bim.seq[!dups,]
logistic.seq<-logistic.seq[!dups,]
key.logistic.seq<-key.logistic.seq[!dups]
# key.logistic.seq[dups][1:10]
}



rownames(logistic.seq)<-key.logistic.seq
TYPE<-rep("SEQ",times=length(key.logistic.seq))
logistic.seq<-cbind(logistic.seq,TYPE,bim.seq[, core.ann],stringsAsFactors=FALSE)



logistic.seq[1:5,]
logistic[1:5,]


##################### do BETA here in case need to flip direction later


if(!("BETA" %in% colnames(logistic))){
BETA=log(as.numeric(logistic[,"OR"]))
logistic<-cbind(logistic,BETA,stringsAsFactors=FALSE)

## exp(BETA +- 1.96SE)=95% CI
## so
## se<--(log(as.numeric(logistic[,"L95"]))-BETA)/1.96
## logistic<-cbind(logistic,BETA,se,stringsAsFactors=FALSE)
## and se = SE in plink so !! SE Is the standard error of BETA!! no need to mess with
}

if(!("BETA" %in% colnames(logistic.seq))){
BETA=log(as.numeric(logistic.seq[,"OR"]))
logistic.seq<-cbind(logistic.seq,BETA,stringsAsFactors=FALSE)

#logistic.seq[,"SE"]<-log(as.numeric(logistic.seq[,"SE"]))
}
print("PAST BETA")



###########################
if("A1" %in% colnames(logistic)){ ### make sure ALT is A1 so effect direction is consistant - dont need to change sign of dirn then

ref.diff<-logistic[,"A1"]!=logistic[,"ALT"] & logistic[,"A1"] !=0
sum(ref.diff)

if(sum(ref.diff)>0){
logistic[ref.diff,][1:10,]
the.ref<-logistic[ref.diff,"REF"]
logistic[ref.diff,"REF"]<-logistic[ref.diff,"ALT"]
logistic[ref.diff,"ALT"]<-the.ref
logistic[ref.diff,"BETA"]<--1*logistic[ref.diff,"BETA"]

#key.logistic.ref.diff<-build.key(logistic[ref.diff,],core.ann)
#rownames(logistic)[ref.diff]<-key.logistic.ref.diff

ref.diff<-logistic[,"A1"]!=logistic[,"ALT"]
sum(ref.diff)

}



sum(is.na(logistic.seq[,"A1"]))
remove<-is.na(logistic.seq[,"A1"]) | is.na(logistic.seq[,"ALT"]) |  is.na(logistic.seq[,"REF"]) # crap results
sum(remove)
logistic.seq<-logistic.seq[!remove,]

ref.diff<-logistic.seq[,"A1"]!=logistic.seq[,"ALT"]  & logistic.seq[,"A1"] !=0
sum(ref.diff)

#logistic.seq[ref.diff,][1:10,]

if(sum(ref.diff)>0){
logistic.seq[ref.diff,][1:10,]
the.ref<-logistic.seq[ref.diff,"REF"]
logistic.seq[ref.diff,"REF"]<-logistic.seq[ref.diff,"ALT"]
logistic.seq[ref.diff,"ALT"]<-the.ref
logistic.seq[ref.diff,"BETA"]<--1*logistic.seq[ref.diff,"BETA"]

#key.logistic.seq.ref.diff<-build.key(logistic.seq[ref.diff,],core.ann)
#dups<-duplicated(key.logistic.seq.ref.diff)
#print(sum(dups))
#made.dups<-key.logistic.seq.ref.diff %in% rownames(logistic.seq)[!ref.diff]
#rownames(logistic.seq)[ref.diff]<-key.logistic.seq.ref.diff

ref.diff<-logistic.seq[,"A1"]!=logistic.seq[,"ALT"]
sum(ref.diff)
}
# redo.start.end.annovar 

}

### using minor allele now so should be fineetal will sort out differences
## key.logistic.seq.ref.diff[made.dups]

## target<-grep("10:7804196:-:A",rownames(logistic.seq))
## logistic.seq[(target-1):(target+5),]
## 10:70241990:A:-   10                10:70241990:70241990:-:A:indel 70241990  -
## 10:70241990:-:A   10                10:70241990:70241990:A:-:indel 70241990  -


############## need to reset the rownames now case the flip above will have changed sum
############## need to reset the rownames now case the flip above will have changed sum
  
core.ann<-c("chr","start","REF","ALT")
key.logistic<-build.key(logistic,core.ann)
dups<-duplicated(key.logistic)
sum(dups) #804
## dups<-unique(key.logistic[dups])
## dups<-key.logistic %in% dups
## logistic[dups,][1:10,] #### dups all have the same valuse so can remove any one 
bim.chip<-bim.chip[!dups,]
logistic<-logistic[!dups,]
key.logistic<-key.logistic[!dups]
rownames(logistic)<-key.logistic

key.logistic.seq<-build.key(logistic.seq,core.ann)
dups<-duplicated(key.logistic.seq)
sum(dups)
if(sum(dups)>0){ 
## dups<-unique(key.logistic.seq[dups])
## dups<-key.logistic.seq %in% dups
#logistic.seq[dups,][1:10,] #### dups all all unwound element not caught with end 
bim.seq<-bim.seq[!dups,]
logistic.seq<-logistic.seq[!dups,]
key.logistic.seq<-key.logistic.seq[!dups]
# key.logistic.seq[dups][1:10]
}

rownames(logistic.seq)<-key.logistic.seq
  


############################ change this so w.r.t smaller logistic file? hold1<-logistic hold2<-logistic.seq  logistic<-hold1
logistic.seq[1:5,]
logistic[1:5,]
################### match up snp ids:
posns<-match(rownames(logistic),rownames(logistic.seq))
missing<-is.na(posns)
sum(!missing)

## uqdi<-grepl("^UQDI",logistic[,"SNP"])
## sum(uqdi)
## sum(missing)
## test<-missing & uqdi
## sum(test)

## test<-missing & !ok.miss
## sum(test)
## logistic[test,"SNP"]
## logistic.seq[grep("1:152193954:G:C",rownames(logistic.seq)),]
## logistic.seq[posns[test],]

sum(missing)
if(sum(missing)>0){ ## if in seq then flip and try again
ref.diff<-missing
the.ref<-logistic[ref.diff,"REF"]
logistic[ref.diff,"REF"]<-logistic[ref.diff,"ALT"]
logistic[ref.diff,"ALT"]<-the.ref
logistic[ref.diff,"BETA"]<--1*logistic[ref.diff,"BETA"]

#key.logistic.seq.ref.diff<-build.key(logistic.seq[ref.diff,],core.ann)
#dups<-duplicated(key.logistic.seq.ref.diff)
#print(sum(dups))
#made.dups<-key.logistic.seq.ref.diff %in% rownames(logistic.seq)[!ref.diff]
#rownames(logistic.seq)[ref.diff]<-key.logistic.seq.ref.diff


key.logistic.ref.diff<-build.key(logistic[ref.diff,],core.ann) ## flipped value

made.dups<-key.logistic.ref.diff %in% rownames(logistic)[!ref.diff]
make.dup.value.ori<-rownames(logistic)[ref.diff][made.dups] ## original value

#grep("1:152080604:G:C",rownames(logistic))

### make.dups are flips that are already present tis casue rownames assigmnet to fail
if( sum(made.dups)>0 ){
  a.made.dups<-rownames(logistic) %in% make.dup.value.ori # location not to be flipped
  ref.diff<- ref.diff & !a.made.dups
  key.logistic.ref.diff<- key.logistic.ref.diff[!made.dups]

### swap back that duplicate
the.ref<-logistic[a.made.dups,"REF"]
logistic[a.made.dups,"REF"]<-logistic[a.made.dups,"ALT"]
logistic[a.made.dups,"ALT"]<-the.ref
logistic[a.made.dups,"BETA"]<--1*logistic[a.made.dups,"BETA"]
  
}

  
rownames(logistic)[ref.diff]<-key.logistic.ref.diff


#### now swap back the ones that were not matched - these are non-eomes SNPS
posns<-match(rownames(logistic),rownames(logistic.seq))
missing<-is.na(posns)
sum(missing) # these missing are exomechip intergenic typically - swap THEM back to original

if( sum(made.dups)>0 ){ # don't swap back these as Ihave swapped in the above
  missing<-missing & ! a.made.dups
}


if(sum(missing)>0){ ## if in seq then flip and try again
ref.diff<-missing
the.ref<-logistic[ref.diff,"REF"]
logistic[ref.diff,"REF"]<-logistic[ref.diff,"ALT"]
logistic[ref.diff,"ALT"]<-the.ref
logistic[ref.diff,"BETA"]<--1*logistic[ref.diff,"BETA"]

## key.logistic.ref.diff<-build.key(logistic[ref.diff,],core.ann)
## rownames(logistic)[ref.diff]<-key.logistic.ref.diff
}

} ###

posns<-match(rownames(logistic),rownames(logistic.seq))
missing<-is.na(posns)
sum(!missing)



# cbind(logistic.seq[posns[!missing],],logistic[!missing,])[1:5,]

logistic.seq[posns[!missing],"SNP"]<-logistic[!missing,"SNP"]

## logistic[1:5,]
## logistic.seq[1:5,]
 ## sum(is.na(logistic[,"chr"]))
 ## sum(is.na(logistic.seq[,"chr"]))
#order.by<-order(as.numeric(logistic[,"P"]),decreasing=FALSE)


#logistic<-logistic[order.by,]

## logistic[!missing,][1:2,]
## logistic.seq[posns[!missing],][1:2,]

setwd(output.dir)

### below is commented as just doing common
## write.table(logistic,file=paste("chip",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## write.table(logistic.seq,file=paste("seq",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

#write.table(logistic[!missing,],file=paste("chip.common.extra",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

write.table(logistic,file=paste("chip.common.extra",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE) # keep all chip
write.table(logistic.seq[posns[!missing],],file=paste("seq.common.extra",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


##print(doing merge)
## common.cols<-colnames(logistic)[colnames(logistic) %in% colnames(logistic.seq)]
## the.merge<-rbind(logistic[,common.cols],logistic.seq[,common.cols])
## the.merge[1:5,]
## order.by<-order(the.merge[,"CHR"],the.merge[,"start"])
## the.merge<-the.merge[order.by,]
## the.merge[1:5,]
## write.table(the.merge,file=paste("merge.extra",traits[i],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
##

## system(paste("sed s/TRAIT/",traits[i],"/ run_config_for_P_and_BETA_TEMPLATE.txt > ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))
## system(paste("./metal ",paste("sample.size.CONFIG",traits[i],"txt",sep="."),sep=""))

system(paste("cp /media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/run_config_for_P_and_BETA_TEMPLATE.txt", output.dir,sep=" "))
system(paste("cp /media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/run_config_for_inverse_varience_TEMPLATE.txt", output.dir,sep=" "))
system(paste("cp /media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/metal", output.dir,sep=" "))

system(paste("sed s/TRAIT/common.extra.",traits[i],"/ run_config_for_P_and_BETA_TEMPLATE.txt > ",paste("sample.size.CONFIG.common.extra",traits[i],"txt",sep="."),sep=""))
system(paste("./metal ",paste("sample.size.CONFIG.common.extra",traits[i],"txt",sep="."),sep=""))


## system(paste("sed s/TRAIT/",traits[i],"/ run_config_for_inverse_varience_TEMPLATE.txt > ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))
## system(paste("./metal ",paste("STDERR.CONFIG",traits[i],"txt",sep="."),sep=""))

system(paste("sed s/TRAIT/common.extra.",traits[i],"/ run_config_for_inverse_varience_TEMPLATE.txt > ",paste("STDERR.CONFIG.common.extra",traits[i],"txt",sep="."),sep=""))
system(paste("./metal ",paste("STDERR.CONFIG.common.extra",traits[i],"txt",sep="."),sep=""))


} ## loop over traits


    
  ############# some testing

##   data<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/final.meta.analysis.single.point.significant.reorg.evoker.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

##   hwe.extra<-  read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP.list.re-checked.ED.chr1-19/all.lists_modified_calls_clean.hwe",header=T,fill=TRUE,stringsAsFactors=FALSE)
##  hwe.extra<- hwe.extra[hwe.extra[,"TEST"]=="ALL",]
##   hwe.extra[1:5,]

## target<-"UQDIexomechipWithAOGC_chr1:28422829:28422829:C:A-1_T_F_2098199532"

## target<-"UQDIexomechipWithAOGC_chr11:111386239:111386239:A:C-1_B_R_2098219387"

## target<-"UQDIexomechipWithAOGC_chr10:114045896:114045896:A:C-1_T_F_2098217256"

## posn.chip<-grep(target,logistic[,"SNP"])
## posn.hwe<-grep(target,hwe.extra[,"SNP"])
## posn.seq<-grep(target,logistic.seq[,"SNP"])
## data.posn<-data[,"SNP"]==target & data[,"traits"] == traits[i] & !is.na(data[,"traits"])  & !is.na(data[,"SNP"]) 

##    hwe.extra[posn.hwe,]
## logistic.seq[posn.seq,]
##   logistic[posn.chip,]
## data[data.posn,]
##   ############# some testing




########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$





######################
finished reprocessing fo redone samples to meta analysiis stage

#------------------------------------------------------------------
# loccation of rescored files    
ROOT.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/2014-06-23_SNPs_HWE.gt.10-13andP.lt.0.001" # location of recorces beds

##IMPORTANT  The bed/list/score file must contain the chromosome name prefix chromsome sufix
## names of the list/score/bed files must be the chromosome
### set up the list names
check.flag<-"Mhairi.ED.Final"
list.file.prefix<-"hwe_andp."
list.file.sufix<-".txt"
### set up the score file names    
score.file.prefix<-list.file.prefix
score.file.sufix<-paste(list.file.sufix,".scores",sep="")
# set up to get the bed.file.names
bed.file.prefix<-"hwe_andp."
bed.file.sufix<-".bed"    
#---------------------------------------------------------------------------------    

#------------------------------------------------------------------
# loccation of rescored files    
ROOT.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/2014-06-23_MhairiRescore_InteresingSNPs" # location of recorces beds

##IMPORTANT  The bed/list/score file must contain the chromosome name prefix chromsome sufix
## names of the list/score/bed files must be the chromosome
### set up the list names

check.flag<-"Mhairi.ED"
list.file.prefix<-"Mhairi.redo."
list.file.sufix<-".txt"
### set up the score file names    
score.file.prefix<-list.file.prefix
score.file.sufix<-paste(list.file.sufix,".scores",sep="")
# set up to get the bed.file.names
bed.file.prefix<-"AOGC_gentrain."
bed.file.sufix<-".bed"    
#---------------------------------------------------------------------------------

    
    
#------------------------------------------------------------------
# loccation of rescored files  #1  
ROOT.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP.list.re-checked.ED.chr1-19" # location of recorces beds

##IMPORTANT  The bed/list/score file must contain the chromosome name prefix chromsome sufix
## names of the list/score/bed files must be the chromosome
### set up the list names
check.flag<-"ED.modified"
    
list.file.prefix<-"meta.snp.test_"
list.file.sufix<-".txt"
### set up the score file names    
score.file.prefix<-list.file.prefix
score.file.sufix<-paste(list.file.sufix,".scores",sep="")
# set up to get the bed.file.names
bed.file.prefix<-"meta.snp.test_"
bed.file.sufix<-".bed"    
#---------------------------------------------------------------------------------    

    

############ location of fam and bim files to use GENERIC
the.bim.root<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Cluster_vis_clean/Cluster_viz/" # .chr.bim added later
the.fam<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Cluster_vis_clean/Cluster_viz/chr1/AOGC_gentrain.fam"

 fam.template.file<-"all.lists_filtered-clean.f.fam" ## has been flipped to match
    

    






##  NEED TO GET BEFORE FILTERED
##
## data<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/final.meta.analysis.single.point.Filtered2.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)











  ############ insert into previos data chip P, meta P x2 hwe for chip


########################## get processed bim files that are all on the forward strand in annovar format"
code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")

######


annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
  
traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")

names(traits)<-rep("logistic",times=length(traits))
names(traits)[1:3]<-rep("assoc",times=3)
traits
traits %in% colnames(ann)


#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
names(traits)<-rep(".assoc.logistic",times=length(traits))
names(traits)[1:3]<-rep(".qassoc",times=3)
traits


## fam.template.file
## # fam.template.file set above in cases where flipping data is required
## fam<-read.table(fam.template.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
## projects<-gsub(".fam$","",fam.template.file)


    
the.seq.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/AOGC_vcf_to_plink" # dirname(fam.template.file)
files<-dir(the.seq.dir)  
projects<-gsub(".bed","",files[grepl(".bed$",files)])  # files[grepl(".bed$",files)]


project.prefix<-"AOGC_to_vcf_filtered_" ## outpt used same name as above
project.prefix.old<-"AOGC_vcf_final_chr"

projects.old<-projects[grepl(paste("^",project.prefix.old,sep=""),projects)]
projects<-projects[grepl(paste("^",project.prefix,sep=""),projects)] #

projects
projects.old

the.chr.old<-gsub(project.prefix.old,"",projects.old)
the.chr<-gsub(project.prefix,"",projects)

sum(the.chr.old !=the.chr)


# Corresdonding in same order for exome chip
 logistic.project.dir<-ROOT.dir
 logistic.project<-gsub(".fam$","",fam.template.file)
 if(!grepl("/$",logistic.project.dir)){logistic.project.dir<-paste(logistic.project.dir,"/",sep="")}

#logistic.files<-paste(logistic.project.dir,traits,".",logistic.project,names(traits),sep="")
hwe.files.chip<-paste(logistic.project.dir,traits,".",logistic.project,".hwe",sep="")
logistic.files<-paste(logistic.project.dir,"chip.common.extra.",traits,sep="")
logistic.seq.files<-paste(logistic.project.dir,"seq.common.extra.",traits,sep="")

#"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/filtered/"
meta.dir<-logistic.project.dir
meta.stderr.files<-paste(meta.dir,"META.STDERR.common.extra.",traits,"1.tbl",sep="")
meta.size.files<-paste(meta.dir,"META.SAMPLE.SIZE.common.extra.",traits,"1.tbl",sep="")

#data<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/final.meta.analysis.single.point.significant.reorg.evoker.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

#write.table(all.sig.final,file="final.meta.analysis.single.point.Filtered2.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
getwd()
setwd(logistic.project.dir)



ann.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/annoatation_short_exome_chip.txt"
ann.short.ori<-read.delim(ann.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)






############
########### want some extra info in data: RS.id, p.hardy.chip Stderr.chip  Stderr.seq 

##   colnames(data)[2]
## colnames(data)[colnames(data)=="traits.i."]<-"traits"
## data[1:20,1:10]
  
  i<-1

 all.meta<-{}
for (i in 1:length(traits)){


print(logistic.files[i])
chip<-read.table(logistic.files[i],header=T,fill=TRUE,stringsAsFactors=FALSE)
seq<-read.table(logistic.seq.files[i],header=T,fill=TRUE,stringsAsFactors=FALSE)
hwe.chip<-read.table(hwe.files.chip[i],header=T,fill=TRUE,stringsAsFactors=FALSE)
stderr<-read.table(meta.stderr.files[i],header=T,fill=TRUE,stringsAsFactors=FALSE)
sample.size<-read.table(meta.size.files[i],header=T,fill=TRUE,stringsAsFactors=FALSE)

ip<-1
for(ip in 1:length(projects)){
setwd(the.seq.dir)
############ read the hwe data
print(paste((traits)[i],projects[ip],"hwe",sep="."))

a.hwe<- read.table(paste((traits)[i],projects[ip],"hwe",sep="."),header=T,fill=TRUE,stringsAsFactors=FALSE)
a.hwe.old<- read.table(paste((traits)[i],projects.old[ip],"hwe",sep="."),header=T,fill=TRUE,stringsAsFactors=FALSE)

a.hwe<-as.matrix(a.hwe) ## colClasses="character"
a.hwe.old<-as.matrix(a.hwe.old)  
##############

if( sum( c("AFF","UNAFF") %in% a.hwe[,"TEST"] )==2 ){ # case/control unwind AFF and UNAFF
  aff<-a.hwe[,"TEST"]=="AFF"
  unaff<-a.hwe[,"TEST"]=="UNAFF"
  GENO.aff<-a.hwe[aff,"GENO"]
  GENO.unaff<-a.hwe[unaff,"GENO"]
 sum(aff)
 sum(unaff)
 if(sum(a.hwe[unaff,"SNP"]!=a.hwe[unaff,"SNP"])>0){print("Error in hwe alignments")}
  a.hwe<-cbind(a.hwe[unaff,c("CHR","SNP","A1","A2","P")],GENO.aff,GENO.unaff,stringsAsFactors=FALSE)
  aff<-a.hwe.old[,"TEST"]=="AFF"
  unaff<-a.hwe.old[,"TEST"]=="UNAFF"
  GENO.aff<-a.hwe.old[aff,"GENO"]
  GENO.unaff<-a.hwe.old[unaff,"GENO"]
 sum(aff)
 sum(unaff)
 if(sum(a.hwe.old[unaff,"SNP"]!=a.hwe.old[unaff,"SNP"])>0){print("Error in hwe alignments")}
  a.hwe.old<-cbind(a.hwe.old[unaff,c("CHR","SNP","A1","A2","P")],GENO.aff,GENO.unaff,stringsAsFactors=FALSE)
}

if(ip==1){hwe<-a.hwe}else{hwe<-rbind(hwe,a.hwe)}
if(ip==1){hwe.old<-a.hwe.old}else{hwe.old<-rbind(hwe.old,a.hwe.old)}

#dim(a.logistic.seq)
} ## loop over ip


if( sum( c("AFF","UNAFF") %in% hwe.chip[,"TEST"] )==2 ){ # case/control unwind AFF and UNAFF
  aff<-hwe.chip[,"TEST"]=="AFF"
  unaff<-hwe.chip[,"TEST"]=="UNAFF"
  GENO.aff<-hwe.chip[aff,"GENO"]
  GENO.unaff<-hwe.chip[unaff,"GENO"]
 sum(aff)
 sum(unaff)
 if(sum(hwe.chip[unaff,"SNP"]!=hwe.chip[unaff,"SNP"])>0){print("Error in hwe alignments")}
  hwe.chip<-cbind(hwe.chip[unaff,c("CHR","SNP","A1","A2","P")],GENO.aff,GENO.unaff,stringsAsFactors=FALSE)
}


setwd(logistic.project.dir)



chip[1:5,]
seq[1:5,]
#seq.unfilt[1:5,]
sample.size[1:5,]
stderr[1:5,]
hwe.chip[1:5,]
hwe[1:5,]
hwe.old[1:5,]

core.ann<-c("chr","start","REF","ALT")
seq.key<-build.key(seq,core.ann)


############## Sort out hwe for sequencing
posns<-match(hwe.old[,"SNP"],hwe[,"SNP"])
hwe<-hwe[posns,]

grep("^UQDI",hwe.old[,"SNP"])

hwe.old.key<-strsplit(hwe.old[,"SNP"],split=":")
start<-unlist(lapply(hwe.old.key,function(x) x[2]))
start[1:10]
hwe.old<-cbind(hwe.old,start)

hwe.key<-build.key(hwe.old,c("CHR","start","A2","A1"))

hwe.key[1:10]
posns<-match(seq.key,hwe.key)
missing<-is.na(posns)
hwe<-hwe[posns,]
hwe.old<-hwe.old[posns,]
hwe.old[,"SNP"]<-seq[,"SNP"]
hwe[,"SNP"]<-seq[,"SNP"]
###############################

dim(chip)
dim(seq)
dim(stderr)
## dim(logistic.seq)


## logistic.seq[1:5,]
## dim(logistic.seq)
## dim(seq)

######### alignen all to chip

posns<-match(chip[,"SNP"],seq[,"SNP"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
seq<-seq[posns,]

posns<-match(chip[,"SNP"],ann.short.ori[,"SNP"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
ann.short<-ann.short.ori[posns,]

posns<-match(chip[,"SNP"],stderr[,"MarkerName"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
stderr<-stderr[posns,]

posns<-match(chip[,"SNP"],sample.size[,"MarkerName"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
sample.size<-sample.size[posns,]

posns<-match(chip[,"SNP"],hwe[,"SNP"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
hwe<-hwe[posns,]

posns<-match(chip[,"SNP"],hwe.old[,"SNP"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
hwe.old<-hwe.old[posns,]

posns<-match(chip[,"SNP"],hwe.chip[,"SNP"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
hwe.chip<-hwe.chip[posns,]


ann.short[1:5,]


dim(seq)
dim(chip)
dim(stderr)
dim(sample.size)
dim(ann.short)
dim(hwe)
dim(hwe.old)
dim(hwe.chip)


P.value.StdErr<-stderr[,"P.value"]
P.value.Size<-sample.size[,"P.value"]
Weight<-sample.size[,"Weight"]

P.Seq.Filt<-seq[,"P"]
P.Seq<-rep(NA,times=length(P.Seq.Filt))
P.Chip<-chip[,"P"]

BETA.chip<-chip[,"BETA"]
BETA.seq<-seq[,"BETA"]

SE.chip<-chip[,"SE"]
SE.seq<-seq[,"SE"]

#### hwe by trait and pheno so can modify hwe.ori<-hwe ; hwe.old.ori<-hwe.old # hwe<-hwe.ori ; hwe.old<-hwe.old.ori
hwe<-hwe[,(grepl("^GENO",colnames(hwe)) |  colnames(hwe)=="P")]
hwe.old<-hwe.old[,(grepl("^GENO",colnames(hwe.old)) |  colnames(hwe.old)=="P")]
hwe.chip<-hwe.chip[,(grepl("^GENO",colnames(hwe.chip)) |  colnames(hwe.chip)=="P")]

colnames(hwe)[colnames(hwe)=="P"]<-"P.hardy"
colnames(hwe.old)[colnames(hwe.old)=="P"]<-"P.hardy"
colnames(hwe.chip)[colnames(hwe.chip)=="P"]<-"P.hardy"

colnames(hwe)<-paste(colnames(hwe),"Seq.Filt",sep=".")
colnames(hwe.old)<-paste(colnames(hwe.old),"Seq",sep=".")
colnames(hwe.chip)<-paste(colnames(hwe.chip),"Chip",sep=".")

all.hwe<-cbind(hwe,hwe.old,hwe.chip)

all.hwe[1:5,]
all.hwe<-as.matrix(all.hwe)

extra.cols<-c("GENO.Chip","GENO.Seq.Filt","GENO.Seq","GENO.aff.Chip","GENO.aff.Seq.Filt","GENO.unaff.Chip","GENO.unaff.Seq.Filt","GENO.aff.Seq","GENO.unaff.Seq","P.hardy.Chip","P.hardy.Seq.Filt","P.hardy.Seq")
extra<-matrix(data=NA,nrow=dim(all.hwe)[1],ncol=length(extra.cols))
colnames(extra)<-extra.cols

extra[,extra.cols[extra.cols %in% colnames(all.hwe)]]<-all.hwe[,extra.cols[extra.cols %in% colnames(all.hwe)]]
                                        # reorganise HWE
extra[1:5,]

dim(extra)
length(BETA.chip)


wanted.cols<-c("rs.id","design","refGene..type","Gene.Names","description","gerp.scores","PolyPhen.desc","PolyPhen.scores","SIFT.desc","SIFT.scores","skeletome","mouse.defect","sewell.cycling","Dequeant.cycling","ingenuity.bone.genes","Consequence.Embl","Amino_acids.Embl","MAF.lt.0.001","MAF.lt.0.5","Final_Score")


#wanted<-c("CHR","SNP","POS","A1","TEST","NMISS","OR","SE","L95","U95","STAT","P","design")

## cols.to.halpview.fix<-c("GENO","design","refGene..type","Gene.Names","description","gerp.scores","PolyPhen.desc","PolyPhen.scores","SIFT.desc","SIFT.scores","skeletome","mouse.defect","sewell.cycling","Dequeant.cycling","ingenuity.bone.genes","Consequence.Embl","Amino_acids.Embl","MAF","MAF.lt.0.001","MAF.lt.0.5")

wanted.cols %in% colnames(ann.short)

meta<-cbind(traits[i],stderr[,c("MarkerName","Allele1","Allele2","Effect","StdErr","Direction")],P.value.StdErr,P.value.Size,Weight,P.Chip,P.Seq.Filt,P.Seq,extra,ann.short[,wanted.cols],chip[,c("chr","start","REF","ALT","SNP")],BETA.chip,BETA.seq,SE.chip,SE.seq,stringsAsFactors=FALSE)

dim(meta)
colnames(meta)

if(is.null(dim(all.meta))){
  all.meta<-meta
}else{
  all.meta<-rbind(all.meta,meta)
}


} ## loop over traits

###################################

##combine with original data # R:2 is old ED # R:3 is latest

    
    getwd()
    dim(all.meta)

    save(list="all.meta",file="all.meta.RData")
colnames(all.meta)

colnames(all.meta)[colnames(all.meta) == "traits[i]"]<-"traits"
colnames(all.meta)[colnames(all.meta) == "GENO.Chip"]<-"GENO"

ann.short<-ann.short.ori
dim(ann.short)

nyway the recoded Plink files are in:

 

/mnt/UQCCG/UQCCG-Projects/PAUL_LEO/AOGC exome chip core
genotyping/Set1_New_Plink_Files/Set1_recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.*

/mnt/UQCCG/UQCCG-Projects/PAUL_LEO/AOGC exome chip core
genotyping/Set2_New_Plink_Files/Set2_recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.*

/mnt/UQCCG/UQCCG-Projects/PAUL_LEO/AOGC exome chip core
genotyping/Set3_New_Plink_Files/
Set3_recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.*



## code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
## setwd(code.dir)
source("/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/annotate_SNPs_subroutines.r")



######## do one at a time Must be in correct order
    #------------------------------------------------------------------
# loccation of rescored files  #1  
ROOT.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP.list.re-checked.ED.chr1-19" # location of recorces beds
check.flag<-"ED.modified"
#24794    54  1078 snps
#---------------------------------------------------------------------------------      

#------------------------------------------------------------------
# loccation of rescored files  #2   
ROOT.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/2014-06-23_SNPs_HWE.gt.10-13andP.lt.0.001" # location of recorces beds
check.flag<-"Mhairi.ED.Final"
# 13800    54   600 snps
#---------------------------------------------------------------------------------    
#------------------------------------------------------------------
# loccation of rescored files  #3  
ROOT.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/2014-06-23_MhairiRescore_InteresingSNPs" # location of recorces beds
check.flag<-"Mhairi.ED"
#/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/2014-06-23_MhairiRescore_InteresingSNPs" # location of recorces beds
# 2162   54 94 snps
#
#---------------------------------------------------------------------------------


## /media/UQCCG/UQCCG-Projects/AOGC_exome_chip/IBD/related.to.remove.recode.txt # m=284
## /media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.recode.txt # n=373 contains

## pleo@DI-LW-BRN-011:/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/exclusions$ grep AOGC-02-0398 *
## Samples_To_Exclude_Contaminated.csv:AOGC-02-0398	AOGC-02-0398
## Samples_To_Exclude.txt:AOGC-02-0398

#AOGC-02-0398

/media/UQCCG/UQCCG-Projects/PAUL_LEO/AOGC exome chip core genotyping

fam.new<-read.table("/media/UQCCG/UQCCG-Projects/PAUL_LEO/AOGC exome chip core genotyping/all.lists_filtered-clean.f.fam",header=F,fill=TRUE,stringsAsFactors=FALSE)
fam<-read.table("/media/UQCCG/UQCCG-Projects/PAUL_LEO/AOGC exome chip core genotyping/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.fam",header=F,fill=TRUE,stringsAsFactors=FALSE)
fam<-read.table("/media/UQCCG/UQCCG-Projects/PAUL_LEO/AOGC exome chip core genotyping/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL2.fam",header=F,fill=TRUE,stringsAsFactors=FALSE)

related.to.remove<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.recode.txt",header=F,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
related.to.remove[1:5,]

/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/IBD/related.to.remove.recode.txt

recodes<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/pheno_gwas_id_missmatches.csv",,header=T,fill=TRUE,stringsAsFactors=FALSE)
recodes[1:5,]


fam<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/2014-06-23_MhairiRescore_InteresingSNPs/BMD_EFF_STD_HIP.fam",header=F,fill=TRUE,stringsAsFactors=FALSE)
remove<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/2014-06-23_MhairiRescore_InteresingSNPs/remove_from.BMD_EFF_STD_HIP",header=F,fill=TRUE,stringsAsFactors=FALSE)

fam<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/BMD_EFF_STD_HIP.fam",header=F,fill=TRUE,stringsAsFactors=FALSE)


fam<-fam[!(fam[,1] %in% remove[,1]),]
dim(fam)
dim(fam.new)
recodes[1:5,]

sum(fam[,1] %in% recodes[,1] )
#fam[test,1] %in% related.to.remove[,1]
#ann[,"PATIENT"]
sum((fam[,1] %in% related.to.remove[,1]))

fam[(fam[,1] %in% related.to.remove[,1]),]


fam[!(fam[,1] %in% fam.new[,1]),][1:10,]

recodes[1:5,]
recodes[,2] %in% fam.new[,2] # fam uses GWASID
recodes[,2] %in% fam[,1] # fam uses GWASID
recodes[,2] %in% exclude[,1] # exclude used GWASID
recodes[,2] %in% ann[,"PATIENT"]  # ann used GWASID


recodes[,1] %in% fam.new[,1] # fam.new (original data  uses PhenoID
recodes[1:5,]
# fam.new (original data  uses PhenoID

### convert from 
posns<-match(fam.new[,1],recodes[,"PhenoID"])
missing<-is.na(posns)
sum(!missing)
fam.new[!missing,1]<-recodes[posns[!missing],"GWASID"]
fam.new[!missing,2]<-recodes[posns[!missing],"GWASID"]
sum(duplicated(fam.new[,1]))

/media/UQCCG/UQCCG-Projects/PAUL_LEO/AOGC exome chip core genotyping/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.bed

setwd("/media/UQCCG/UQCCG-Projects/PAUL_LEO/AOGC exome chip core genotyping/")
system( paste("plink","--bfile","recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL","--remove","/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.recode.txt","--make-bed","--out","recode_sampls *FINAL*leNames_zCall_AOGC.with.evoker_corrected_clean_FINAL2","--noweb",sep=" ") )


fam<-read.table("/media/UQCCG/UQCCG-Projects/PAUL_LEO/AOGC exome chip core genotyping/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL2.fam",header=F,fill=TRUE,stringsAsFactors=FALSE)
test<-!(fam[,1] %in% fam.new[,1])
sum(test) 0 so ok

write.table(fam[,c(1,2)],file="keep.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

plink --bfile recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL --remove /media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.recode.txt --make-bed --out recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL2 --noweb

keep.txt from recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL2
##1 set 2 set and 3 set run through

wc -l all.lists_filtered-clean.f.bim
plink --bfile all.lists_filtered-clean.f --keep keep.txt --out all.lists_filtered-clean.f2  --make-bed --allow-no-sex --noweb
cut -f 2 all.lists_filtered-clean.f2.bim > snps.txt
wc -l snps.txt
plink --bfile recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL2 --exclude snps.txt --make-bed --out recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL3  --allow-no-sex --noweb

plink -bfile recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL3  --bmerge all.lists_filtered-clean.f2.bed all.lists_filtered-clean.f2.bim all.lists_filtered-clean.f2.fam --make-bed --out recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL2 --allow-no-sex --noweb


266163 variants and 7619 people pass filters and QC.
266163 variants and 7619 people pass filters and QC
266163 variants and 7619 people pass filters and QC.

recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL2
move to
/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes
and renamed
recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL
(old data in first <- run <- data)


plink --bfile all.lists_filtered-clean.f --keep keep.txt --out all.lists_filtered-clean.f2  --make-bed --allow-no-sex --noweb

plink --bed all.lists_filtered-clean.f2.bed --bim all.lists_filtered-clean.f2.bim --fam BMD_EFF_STD_HIP.fam --remove remove_from.BMD_EFF_STD_HIP --assoc --allow-no-sex --out test.no.bad --noweb


"plink --bed all.lists_filtered-clean.f.bed --bim all.lists_filtered-clean.f.bim --fam BMD_EFF_STD_HIP.fam --remove remove_from.BMD_EFF_STD_HIP --assoc --allow-no-sex --out BMD_EFF_STD_HIP.all.lists_filtered-clean.f --noweb"


## Warning: Variants 'exm596' and 'exm2253593' have the same position.
## Warning: Variants 'exm2277033' and 'exm1771279' have the same position.
## Warning: Variants 'exm5910' and 'exm2250954_ver2' have the same position.
## 899 more same-position warnings: see log file.



## replace fam file with recoded
#setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP.list.re-checked.ED.chr1-19")
write.table(fam.new,file="all.lists_filtered.fam",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP.list.re-checked.ED.chr1-19/all.lists_filtered-clean.f.bim

ROOT.dirs<-c(
"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP.list.re-checked.ED.chr1-19",
"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/2014-06-23_SNPs_HWE.gt.10-13andP.lt.0.001",
"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/2014-06-23_MhairiRescore_InteresingSNPs"
             )

check.flags<-c(
"ED.modified",
"Mhairi.ED.Final",
"Mhairi.ED"   )           

i<-3
meta.full<-{}
### do above in order
for(i in 1:length(ROOT.dirs)){
ROOT.dir<-ROOT.dirs[i]
check.flag<-check.flags[i]

setwd(ROOT.dir)
load("all.meta.RData")
colnames(all.meta)
length(unique(all.meta$SNP))

colnames(all.meta)[colnames(all.meta) == "traits[i]"]<-"traits"
colnames(all.meta)[colnames(all.meta) == "GENO.Chip"]<-"GENO"
colnames(all.meta)
Plan<-rep(check.flag,times=dim(all.meta)[1])
all.meta<-cbind(all.meta,Plan)

print(dim(all.meta))
sum(is.na(all.meta[,"SNP"]))
sum(all.meta[,"SNP"]!=all.meta[,"MarkerName"] & !is.na(all.meta[,"MarkerName"]))
test<-is.na(all.meta[,"SNP"])



if(is.null(meta.full)){meta.full<-all.meta}else{

  dup.snps<-meta.full[,"SNP"] %in% all.meta[,"SNP"]
  print(sum(dup.snps))
  key.all.meta<-build.key(all.meta,c("SNP","traits"))
  key.meta.full<-build.key(meta.full,c("SNP","traits"))

  posns<-match(key.meta.full,key.all.meta)
  missing<-is.na(posns)
 sum(!missing)
  common<-cbind(meta.full[!missing,c("SNP","traits","P.Chip")],all.meta[posns[!missing],c("SNP","traits","P.Chip")])
  print(dim(common))
  if(dim(common)[1]!=0){print("warning problem")}
  print(dim(meta.full))
    meta.full<-rbind(meta.full,all.meta)
  print(dim(meta.full))
  
  
                         }
}


#meta.full.ori<-meta.full
write.table(meta.full,file="/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/rechecked-redone-genotypes.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)  ########## these are all the resone genotypes
########################### ALL to existing data ##########

#curent.working.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/first cut short list ori.txt"
#curent.working.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/first cut short list ori.Mhairi.ED.txt"
curent.working.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/rechecked-redone-genotypes.txt"

all.meta<-read.delim(curent.working.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
#data<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/final.meta.analysis.single.point.significant.reorg.evoker.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)


### tring to work out
data<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/final.meta.analysis.single.point.Filtered2.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) # orginal one but does not have "plan" see line 3900
###
 already.checked<-"/media/UQCCG-Analysis/AOGC_exome_chip/meta_analysis_single_point/non-filtered-seq/META.BOTH.common.all_p0.0005.postED review.txt"
 chked<-read.delim(already.checked,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)



## data<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/final.meta.analysis.Chip.centric.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ## does not contain all traints
## old.plans<-data[,c("SNP","traits","Plan")]

colnames(data)[colnames(data)=="traits.i."]<-"traits"
table(data$traits)
table(all.meta$traits)
dim(all.meta)
dim(data)


colnames(all.meta)
colnames(data)

colnames(all.meta)[ !(colnames(all.meta) %in% colnames(data))]
colnames(data)[ !(colnames(data) %in% colnames(all.meta))]
colnames(data)

if(sum(!extra.cols %in% colnames(extra)) > 0){ ## must be zero! else below assignmnet will fail
extra.cols<-intersect(colnames(data),colnames(all.meta))
}else{
  extra.cols<-colnames(all.meta)
}

key.data<-build.key(data,c("traits","MarkerName"))
key.meta<-build.key(all.meta,c("traits","MarkerName"))

posns<-match(key.meta,key.data)
missing<-is.na(posns)
sum(missing)

    
length(key.meta)
length(missing)
length(posns)
dim(all.meta)
length(posns[!missing])
sum(!missing)

colnames(all.meta)
colnames(data)
## all.meta[missing,][1:15,1:5]
## grep("UQDIexomechipWithAOGC_chr7:141619286:141619286:T:G-1_B_F_2098213710",all.meta[,"MarkerName"])
all.meta<-as.matrix(all.meta)

# options(width=200)
data[posns[!missing],extra.cols[extra.cols %in% colnames(data)]][1:15,c("traits","MarkerName","P.Chip","P.Seq.Filt")]
all.meta[!missing, extra.cols[extra.cols %in% colnames(all.meta)] ] [1:15,c("traits","MarkerName","P.Chip","P.Seq.Filt")]


data[posns[!missing],extra.cols[extra.cols %in% colnames(data)]]<-all.meta[!missing, extra.cols[extra.cols %in% colnames(all.meta)] ]

data[posns[!missing],extra.cols[extra.cols %in% colnames(data)]][1:15,c("traits","MarkerName","P.Chip","P.Seq.Filt")]
    
data[posns[!missing],"Plan"]<-paste(data[posns[!missing],"Plan"],check.flag,sep=":")
data[posns[!missing],"Plan"]<-gsub("^NA:","",data[posns[!missing],"Plan"])


## data2[posns[!missing],extra.cols[extra.cols %in% colnames(data)]][15:35,c("traits","MarkerName","P.Chip","P.Seq.Filt")]
## all.meta[!missing, extra.cols[extra.cols %in% colnames(all.meta)] ] [15:35,c("traits","MarkerName","P.Chip","P.Seq.Filt")]

extra.cols.add<-colnames(data)
extra<-matrix(data=NA,nrow=sum(missing),ncol=length(extra.cols.add))
colnames(extra)<-extra.cols.add
dim(data)
dim(extra)
sum(!missing)

extra.cols<-colnames(all.meta)

if(sum(!extra.cols %in% colnames(extra)) > 0){ ## must be zero! else below assignmnet will fail
extra.cols<-insect(colnames(data),colnames(all.meta))
}else{
  extra.cols<-colnames(all.meta)
}

    
colnames(data)
colnames(extra)
dim(all.meta[missing, extra.cols[extra.cols %in% colnames(all.meta)] ])
dim(extra)
#all.meta<-as.matrix(all.meta)

extra[,extra.cols[extra.cols %in% colnames(extra)]]<-all.meta[missing, extra.cols[extra.cols %in% colnames(all.meta)] ]

extra[,"Plan"]<-check.flag

data<-rbind(data,extra)

table(data[,"Hit"])
table(data[,"Plan"])
getwd()
out.file<-paste(gsub(".txt$","",curent.working.file),".",check.flag,".txt",sep="")
    out.file

colnames(data)
data[is.na(data[,"P.Chip"]),"P.Chip"]<-99  # [1:5]
data[is.na(data[,"P.Seq.Filt"]),"P.Seq.Filt"] <-99
data[is.na(data[,"P.Seq"]),"P.Seq"]<-99



    sum(is.na(data[,"P_min"]))
sum(is.na(data[,"P.value.StdErr"]))


no.stderr<-is.na(data[,"P.value.StdErr"]) | is.na(as.numeric(data[,"P.value.StdErr"]))
    sum(no.stderr)
fix<-pmin(as.numeric(data[no.stderr,c("P.Chip")]),as.numeric(data[no.stderr,"P.Seq.Filt"]),as.numeric(data[no.stderr,"P.Seq"]))
fix
data[no.stderr,"P.value.StdErr"]<-fix
test<-is.na(as.numeric(data[,"P.value.StdErr"]))
    sum(test)

data[test,][1:5,]

   

########### get ID sorted out marker- snp exectp where one is NA:
test<-data[,"MarkerName"] != data[,"SNP"] & !is.na(data[,"MarkerName"]) & !is.na(data[,"SNP"])    
sum(test) # 0
test<-is.na(data[,"SNP"])
data[test,"SNP"] <- data[test,"MarkerName"]
    sum(test)
test<-is.na(data[,"MarkerName"])
       sum(test)
data[test,"MarkerName"] <- data[test,"SNP"]
 

######### set p_min correctly
test<-is.na(data[,"P_min"]) | is.na(as.numeric(data[,"P_min"]))
sum(test)

data[test,][1:5,]
    
lowest.p<-tapply(data[,"P.value.StdErr"],data[,"SNP"],function(x) min(as.numeric(x,na.rm=TRUE)))


lowest.p[1:40]
order.by<-order(as.numeric(lowest.p),decreasing=FALSE)
lowest.p<-lowest.p[order.by]

posns<-match(data[,"SNP"],names(lowest.p))
missing<-is.na(posns)
sum(missing)
length(posns)

posns[1:20]
dim(data)
cbind(data[,"P_min"],lowest.p[posns])[test,][1:50,]


data[grep("exm1346556",data[,"SNP"]),c("traits","SNP","P.value.StdErr")]

data[,"P_min"]<-lowest.p[posns]
#############################################
    
out.file    
write.table(data,file=out.file,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

















data<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/final.meta.analysis.single.point.Filtered2.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)


gefos<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP_lists_used_in_exomeChip/GEFOS SNPs in final file.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

dim(gefos)

gefos[1:5,]

data[1:5,1:10]

gefos[,1] %in%  data[,"MarkerName"]

dim( data)

markers<-unique(data[,"MarkerName"])
length(markers)





#########################  ANNOTATE META ANALYSIS/media/UQCCG-Analysis
#########################  ANNOTATE META ANALYSIS/media/UQCCG-Analysis
#########################  ANNOTATE META ANALYSIS/media/UQCCG-Analysis
#########################  ANNOTATE META ANALYSIS/media/UQCCG-Analysis
#########################  ANNOTATE META ANALYSIS/media/UQCCG-Analysis
########### this so go get CHIP focussed data ######


#traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")

code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")

traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")
#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
names(traits)<-rep(".assoc.logistic",times=length(traits))
names(traits)[1:3]<-rep(".qassoc",times=3)
traits




hwe.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/zCall_AOGC.with.evoker_corrected_clean_FINAL.hwe"
hwe<-read.table(hwe.file,header=T,fill=TRUE,stringsAsFactors=FALSE)
wanted<-grepl("ALL",hwe[,"TEST"])
hwe.ori<-hwe[wanted,]
hwe.ori[1:5,]

frq.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/zCall_AOGC.with.evoker_corrected_clean_FINAL.frq"
frq.ori<-read.table(frq.file,header=T,fill=TRUE,stringsAsFactors=FALSE)
frq.ori[1:5,]

ann.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/annoatation_short_exome_chip.txt"
ann.short<-read.delim(ann.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

posns<-match(ann.short[,"SNP"],hwe.ori[,"SNP"])
missing<-is.na(posns)
sum(missing)
hwe.ori<-hwe.ori[posns,]

posns<-match(ann.short[,"SNP"],frq.ori[,"SNP"])
missing<-is.na(posns)
sum(missing)
frq.ori<-frq.ori[posns,]

ann.short.ori<-cbind(ann.short,hwe.ori,frq.ori,stringsAsFactors=FALSE)
ann.short[1:5,]


the.seq.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/AOGC_vcf_to_plink" # dirname(fam.template.file)
files<-dir(the.seq.dir)  
projects<-gsub(".bed","",files[grepl(".bed$",files)])  # files[grepl(".bed$",files)]
project.prefix.old<-"AOGC_vcf_final_chr"
project.prefix<-"AOGC_to_vcf_filtered_"
projects.old<-projects[grepl(paste("^",project.prefix.old,sep=""),projects)]
projects<-projects[grepl(paste("^",project.prefix,sep=""),projects)] #

projects
projects.old

the.chr.old<-gsub(project.prefix.old,"",projects.old)
the.chr<-gsub(project.prefix,"",projects)

sum(the.chr.old !=the.chr)
the.chr<-the.chr[!(the.chr %in% c(23,24))]


### to use this would need to add reat the the single point sequening results to get the key right from the SNP name and then align those.
## ref.cohort<-paste("/media/UQCCG/Sequencing/Projects/AOGC-NGS/Australian Reference Genome Cohort/AOGC_Controls/AOGC-Genotyping.output.chr",the.chr,".AOGC_ALL.controls.txt",sep="")

## ref.cohort
## ip<-1
## for(ip in 1:length(ref.cohort)){
##   print(ref.cohort[ip])
## a.ref<- read.table(ref.cohort[ip],header=T,fill=TRUE,stringsAsFactors=FALSE,sep="\t",quote="")
## if(ip==1){ref.ann<-a.ref}else{ref.ann<-rbind(ref.ann,a.ref)}
## }
  
## ref.ann.key<-build.key(ref.ann,c("chr","start","end","REF","ALT","TYPE"))
## ref.ann.key[1:20]

########### need to load in seq from one of the runs below
## seq.key<-seq[,c("chr","start","start","REF","ALT","SNP")]
## colnames(seq.key)<-c("chr","start","end","REF","ALT","SNP")

## seq.key[,c("chr","start","end","REF","ALT")]<-redo.start.end.annovar(seq.key[,c("chr","start","end","REF","ALT")]) ## put in annovar format
## seq.key<-get.forward.strand.allele(seq.key)
## dim(seq.key)
## seq.key[1:5,]
## seq.key[,"chr"]<-gsub("^chr","",seq.key[,"chr"])
## seq.key.use<-build.key(seq.key,c("chr","start","end","REF","ALT"))
## ref.ann.key.use<-build.key(ref.ann,c("chr","start","end","REF","ALT"))

## seq.key.use[1:5]
## ref.ann.key.use[1:5]
## posns<-match(seq.key.use,ref.ann.key.use)
## missing<-is.na(posns)
## sum(missing)
## save(list=c("ref.ann.key","ref.ann.key.use","ref.ann","seq.key","seq.key.use"),file="/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/ref.ann.RData")

setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point")

filtered.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/filtered"
unfiltered.dir<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/notFiltered"

#load("/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward.RData")
#load("/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_fromVCF.RData")
## load("/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_A1_A2.RData") # contains "bim.chip" "bim.seq"

## bim.seq[1:5,]
## bim.chip[1:5,]
## chipkey<-build.key(bim.chip,c("chr","start","REF","ALT"))
## seqkey<-build.key(bim.seq,c("chr","start","REF","ALT"))
## posns<-match(chipkey,seqkey)
## missing<-is.na(posns)
## sum(!missing)
## chipSeqkey<-cbind(bim.chip[!missing,],bim.seq[posns[!missing],"SNP"])
## chipSeqkey[1:5,]
## colnames(chipSeqkey)[7]<-"marker"
## save(list=c("bim.chip","bim.seq","chipSeqkey"),file="/media/UQCCG-Analysis/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_A1_A2_withKEY.RData")

load("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/bim_chip_seq_on_forward_A1_A2_withKEY.RData") # contains "bim.chip" "bim.seq" "chipSeqkey" key using posn and ref/alt

load("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/ref.ann.RData")
bim.chip.ori<-bim.chip
bim.seq.ori<-bim.seq
ref.ann.key.ori<-ref.ann.key
                                 
i<-1
for (i in 10:length(traits)){
#  for (i in 4:9){
#    for (i in 2:3){
#for (i in 1:length(traits)){

  ## for (i in 12:17){
ann.short<-ann.short.ori
ref.ann.key<-ref.ann.key.ori
bim.seq<-bim.seq.ori
setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/AOGC_vcf_to_plink")
ip<-1

for(ip in 1:length(projects)){


if(i %in% 1:3){ # use rescalled data
#  NewBETA
logistic.file.seq<-paste("NewBETA",traits[i],projects[ip],sep=".")
print(logistic.file.seq)
setwd(the.seq.dir)
a.logistic.seq<-read.table(logistic.file.seq,header=T,fill=TRUE,stringsAsFactors=FALSE)
}else{
logistic.file.seq<-paste(traits[i],projects[ip],sep=".")
print(logistic.file.seq)
setwd(the.seq.dir)
a.logistic.seq<-read.table(paste(logistic.file.seq,names(traits)[i],sep=""),header=T,fill=TRUE,stringsAsFactors=FALSE)
}

colnames(a.logistic.seq)[colnames(a.logistic.seq)=="BP"]<-"POS"

if("TEST" %in% colnames(a.logistic.seq)){
wanted<-grepl("ADD",a.logistic.seq[,"TEST"])
a.logistic.seq<-a.logistic.seq[wanted,]
}

if(ip==1){logistic.seq<-a.logistic.seq}else{logistic.seq<-rbind(logistic.seq,a.logistic.seq)}

############ read the hwe data
print(paste((traits)[i],projects[ip],"hwe",sep="."))

a.hwe<- read.table(paste((traits)[i],projects[ip],"hwe",sep="."),header=T,fill=TRUE,stringsAsFactors=FALSE)
a.hwe.old<- read.table(paste((traits)[i],projects.old[ip],"hwe",sep="."),header=T,fill=TRUE,stringsAsFactors=FALSE)

a.hwe<-as.matrix(a.hwe) ## colClasses="character"
a.hwe.old<-as.matrix(a.hwe.old)  
##############

if( sum( c("AFF","UNAFF") %in% a.hwe[,"TEST"] )==2 ){ # case/control unwind AFF and UNAFF
  aff<-a.hwe[,"TEST"]=="AFF"
  unaff<-a.hwe[,"TEST"]=="UNAFF"
  GENO.aff<-a.hwe[aff,"GENO"]
  GENO.unaff<-a.hwe[unaff,"GENO"]
 sum(aff)
 sum(unaff)
 if(sum(a.hwe[unaff,"SNP"]!=a.hwe[unaff,"SNP"])>0){print("Error in hwe alignments")}
  a.hwe<-cbind(a.hwe[unaff,c("CHR","SNP","A1","A2","P")],GENO.aff,GENO.unaff,stringsAsFactors=FALSE)
  aff<-a.hwe.old[,"TEST"]=="AFF"
  unaff<-a.hwe.old[,"TEST"]=="UNAFF"
  GENO.aff<-a.hwe.old[aff,"GENO"]
  GENO.unaff<-a.hwe.old[unaff,"GENO"]
 sum(aff)
 sum(unaff)
 if(sum(a.hwe.old[unaff,"SNP"]!=a.hwe.old[unaff,"SNP"])>0){print("Error in hwe alignments")}
  a.hwe.old<-cbind(a.hwe.old[unaff,c("CHR","SNP","A1","A2","P")],GENO.aff,GENO.unaff,stringsAsFactors=FALSE)
}




if(ip==1){hwe<-a.hwe}else{hwe<-rbind(hwe,a.hwe)}
if(ip==1){hwe.old<-a.hwe.old}else{hwe.old<-rbind(hwe.old,a.hwe.old)}

#dim(a.logistic.seq)
} ## loop over ip






dim(logistic.seq)
dim(bim.seq)
#colnames(logistic.seq)[colnames(logistic.seq)=="BP"]<-"POS"
#gefos[,1] %in%  data[,"MarkerName"]
## logistic.seq[1:5,]
## logistic[1:5,]
## bim.seq[1:5,]

## posns<-match(logistic.seq[,"SNP"],bim.seq[,"SNP"]) # for annotation
## missing<-is.na(posns)
## sum(missing)
## logistic.seq<-logistic.seq[!missing,]
## bim.seq<-bim.seq[posns[!missing],]
###################################################



out.file<-paste("META.BOTH.ALLCHIP.RESCALED.",traits[i],"1.tbl",sep="")
out.file.sig<-paste("META.BOTH.ALLCHIP.RESCALED.SIGNIF2.",traits[i],"1.tbl",sep="")
out.file.ann<-paste("META.BOTH.ALLCHIP.RESCALED.",traits[i],"1.tbl.ann",sep="")


  
setwd(filtered.dir)
## sample.size.file<-paste("META.SAMPLE.SIZE.common.",traits[i],"1.tbl",sep="")
## stderr.file<-paste("META.STDERR.common.",traits[i],"1.tbl",sep="")

sample.size.file<-paste("META.SAMPLE.SIZE.",traits[i],"1.tbl",sep="")
stderr.file<-paste("META.STDERR.",traits[i],"1.tbl",sep="")

chip<-read.delim(paste("chip",traits[i],sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
seq<-read.delim(paste("seq",traits[i],sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

sample.size<-read.delim(sample.size.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
stderr<-read.delim(stderr.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)


setwd(unfiltered.dir)
### chip.unfilt<-read.delim(paste("chip",traits[i],sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) # chip  should be the same
seq.unfilt<-read.delim(paste("seq",traits[i],sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)




chip[1:5,]
seq[1:5,]
seq.unfilt[1:5,]
sample.size[1:5,]
stderr[1:5,]
hwe[1:5,]
hwe.old[1:5,]


dim(seq)
dim(stderr)
dim(logistic.seq)
length(ref.ann.key.ori)

logistic.seq[1:5,]
dim(logistic.seq)
dim(seq)

posns<-match(chip[,"SNP"],seq[,"SNP"]) ## make chip the default id name before.
missing<-is.na(posns)
sum(!missing)
common<-chip[!missing,"SNP"]


#Collect stuff that is
## 1)significant ORIGINALLY or in meta analysis
## 2)in common no matter what

## posns<-match(chip[,"SNP"],seq[,"SNP"]) ## make chip the default id name before.
## missing<-is.na(posns)
## sum(!missing)
## common<-chip[!missing,"SNP"]

## hist(as.numeric(seq[,"NMISS"]))
## hist(as.numeric(chip[,"NMISS"]))

num.samples<-max(as.numeric(seq[,"NMISS"]),na.rm=TRUE)
min.miss<-num.samples-0.2*num.samples
min(as.numeric(seq[,"NMISS"],na.rm=TRUE))
sig.seq<-as.numeric(as.character(seq[,"P"]))< 0.00001 & !is.na(as.numeric(as.character(seq[,"P"]))) & as.numeric(seq[,"NMISS"]) > min.miss
sum(sig.seq)
sig.seq<-seq[sig.seq,"SNP"]

num.samples<-max(as.numeric(chip[,"NMISS"]),na.rm=TRUE)
min.miss<-num.samples-0.2*num.samples
min(as.numeric(chip[,"NMISS"],na.rm=TRUE))
sig.chip<-as.numeric(chip[,"P"])< 0.0001 & !is.na(as.numeric(chip[,"P"])) & as.numeric(chip[,"NMISS"]) > min.miss
sum(sig.chip)
sig.chip<-chip[sig.chip,"SNP"]


to.check<-unique(c(common,sig.seq,sig.chip))

############ put everything in the order of to.check




seq.key[1:5,]
ref.ann.key.use[1:5]
seq.key.use[1:5]
names(seq.key.use)[1:5]
names(ref.ann.key.use)[1:5]

posns<-match(chipSeqkey[,"marker"],hwe[,"SNP"])
missing<-is.na(posns)
hwe[posns[!missing],"SNP"]<-chipSeqkey[!missing,"SNP"]

posns<-match(chipSeqkey[,"marker"],hwe.old[,"SNP"])
missing<-is.na(posns)
hwe.old[posns[!missing],"SNP"]<-chipSeqkey[!missing,"SNP"]

## chipSeqkey[!missing,"SNP"][1:10]
## hwe[posns[!missing],"SNP"][1:10]

## grep("7:141619286",hwe[,"SNP"]) # 1494718
## grep("7:141619286",ref.ann.key.use)
## seq.key.use[grep("7:141619286",ref.ann.key.use)]
## seq.key[grep("7:141619286",ref.ann.key.use),]



posns<-match(chip[,"SNP"],seq.key[,"SNP"])
sum(is.na(posns))
targets<-seq.key.use[posns]
posns<-match(targets,ref.ann.key.use)

## ref.ann.key.ori[1:50]
## ## get meta in same order
## posns<-match(seq[,"SNP"],ref.ann.key.ori)
## missing<-is.na(posns)
## sum(!missing)


                                             

posns<-match(chip[,"SNP"],sample.size[,"MarkerName"])
missing<-is.na(posns)
sum(missing) # sample.size[missing,]
sample.size<-sample.size[posns,]

posns<-match(chip[,"SNP"],stderr[,"MarkerName"])
missing<-is.na(posns)
sum(!missing)
stderr<-stderr[posns,]

posns<-match(chip[,"SNP"],seq[,"SNP"])
missing<-is.na(posns)
sum(missing)
seq<-seq[posns,]

posns<-match(chip[,"SNP"],seq.unfilt[,"SNP"])
missing<-is.na(posns)
sum(missing)
#stderr[missing,][1:20,]
seq.unfilt<-seq.unfilt[posns,]

## posns<-match(stderr[,"MarkerName"],chip.unfilt[,"SNP"])
## missing<-is.na(posns)
## sum(missing)
## #stderr[missing,][1:20,]
## chip.unfilt<-chip.unfilt[posns,]

posns<-match(chip[,"SNP"],ann.short[,"SNP"])
missing<-is.na(posns)
sum(!missing)
ann.short<-ann.short[posns,]

posns<-match(chip[,"SNP"],hwe[,"SNP"])
missing<-is.na(posns)
sum(!missing)
hwe<-hwe[posns,]

posns<-match(chip[,"SNP"],hwe.old[,"SNP"])
missing<-is.na(posns)
sum(!missing)
hwe.old<-hwe.old[posns,]

dim(seq)
dim(chip)
dim(stderr)
dim(sample.size)
dim(ann.short)
dim(hwe)
dim(hwe.old)

P.value.StdErr<-stderr[,"P.value"]
P.value.Size<-sample.size[,"P.value"]
Weight<-sample.size[,"Weight"]

P.Seq.Filt<-seq[,"P"]
P.Seq<-seq.unfilt[,"P"]
P.Chip<-chip[,"P"]
BETA.Seq.Filt<-seq[,"BETA"]
BETA.Seq<-seq.unfilt[,"BETA"]
BETA.Chip<-chip[,"BETA"]
#### hwe by trait and pheno so can modify hwe.ori<-hwe ; hwe.old.ori<-hwe.old # hwe<-hwe.ori ; hwe.old<-hwe.old.ori
hwe<-hwe[,(grepl("^GENO",colnames(hwe)) |  colnames(hwe)=="P")]
hwe.old<-hwe.old[,(grepl("^GENO",colnames(hwe.old)) |  colnames(hwe.old)=="P")]

colnames(hwe)[colnames(hwe)=="P"]<-"P.hardy"
colnames(hwe.old)[colnames(hwe.old)=="P"]<-"P.hardy"

colnames(hwe)<-paste(colnames(hwe),"Seq.Filt",sep=".")
colnames(hwe.old)<-paste(colnames(hwe.old),"Seq",sep=".")


all.hwe<-cbind(hwe,hwe.old)

all.hwe[1:5,]


extra.cols<-c("GENO.Seq.Filt","GENO.Seq","GENO.aff.Seq.Filt","GENO.unaff.Seq.Filt","GENO.aff.Seq","GENO.unaff.Seq","P.hardy.Seq.Filt","P.hardy.Seq")
extra<-matrix(data=NA,nrow=dim(all.hwe)[1],ncol=length(extra.cols))
colnames(extra)<-extra.cols

extra[,extra.cols[extra.cols %in% colnames(all.hwe)]]<-all.hwe[,extra.cols[extra.cols %in% colnames(all.hwe)]]
                                        # reorganise HWE



colnames(chip)[colnames(chip)=="STAT"]<-"T"


wanted.cols<-c("rs.id","GENO","design","refGene..type","Gene.Names","description","gerp.scores","PolyPhen.desc","PolyPhen.scores","SIFT.desc","SIFT.scores","skeletome","mouse.defect","sewell.cycling","Dequeant.cycling","ingenuity.bone.genes","Consequence.Embl","Amino_acids.Embl","MAF","MAF.lt.0.001","MAF.lt.0.5")


#wanted<-c("CHR","SNP","POS","A1","TEST","NMISS","OR","SE","L95","U95","STAT","P","design")

cols.to.halpview.fix<-c("GENO","design","refGene..type","Gene.Names","description","gerp.scores","PolyPhen.desc","PolyPhen.scores","SIFT.desc","SIFT.scores","skeletome","mouse.defect","sewell.cycling","Dequeant.cycling","ingenuity.bone.genes","Consequence.Embl","Amino_acids.Embl","MAF","MAF.lt.0.001","MAF.lt.0.5")



meta<-cbind(chip[,c("SNP","POS","REF","ALT","NMISS","BETA","SE","T","P")],stderr[,c("MarkerName","Allele1","Allele2","Effect","StdErr","Direction")],P.value.StdErr,P.value.Size,Weight,P.Chip,P.Seq.Filt,P.Seq,extra,BETA.Chip,BETA.Seq.Filt,BETA.Seq,ann.short[,wanted.cols],seq[,c("chr","start","REF","ALT","SNP")],stringsAsFactors=FALSE)

write.table(meta,file=paste(out.file.sig,"ALL",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

## gefos[,1] %in%  chip[,"SNP"]
## gefos[,1] %in%  meta[,"SNP"]
## gefos[!(gefos[,1] %in%  meta[,"MarkerName"]),]
##chip[1:5,c("SNP","POS","REF","ALT","NMISS","BETA","SE","T","P")]
## posns<-match(to.check,meta[,"MarkerName"])
## missing<-is.na(posns)
## sum(missing)


## meta<-meta[posns[!missing],]

## setwd("/media/UQCCG-Analysis/AOGC_exome_chip/meta_analysis_single_point")
## write.table(meta,file=out.file.sig,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
}





######################################################################################################  FINAL REDO OF CLUSTER CHECKS

############################## compile merges list of all genotypes to check
############################## compile merges list of all genotypes to check
############################## compile merges list of all genotypes to check
############################## compile merges list of all genotypes to check
############################## compile merges list of all genotypes to check
############################## compile merges list of all genotypes to check

load("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/ref.ann.RData")


ann.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Analysis/annoatation_short_exome_chip.txt"
ann.short.ori<-read.delim(ann.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)



already.checked<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/non-filtered-seq/META.BOTH.common.all_p0.0005.postED review.txt"
chked<-read.delim(already.checked,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)


chked[1:5,]


traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")
#traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
names(traits)<-rep(".assoc.logistic",times=length(traits))
names(traits)[1:3]<-rep(".qassoc",times=3)
traits
i<-1


#c("ref.ann.key","ref.ann","seq.key","seq.key.use")
all.sig[1:5,]
seq.key[1:5,]
## rownames(seq.key)<- seq.key.use
## rownames(ref.ann)<- ref.ann.key.use

#out.file.sig<-paste("META.BOTH.common.RESCALED.SIGNIF.",traits[i],"1.tbl",sep="")
#setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point")
setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/notFiltered")
for (i in 1:length(traits)){
print(i)

paste(out.file.sig,"ALL",sep="")
out.file.sig<-paste("META.BOTH.ALLCHIP.RESCALED.SIGNIF2.",traits[i],"1.tblALL",sep="")
a.sig<-read.delim(out.file.sig,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
if(dim(a.sig)[1]==0){next}
## a.sig[1:5,]

## chked[1:5,]

## colnames(chked)
## colnames(a.sig)





posns<-match(a.sig[,"MarkerName"],seq.key[,"SNP"])
sum(is.na(posns))
targets<-seq.key.use[posns]
posns<-match(targets,ref.ann.key.use)
extra<-ref.ann[posns,c( "chr","start","end","Hetero.ALT.reads","Hetero.REF.reads","Hetero.Read.Balance")]
extra<-ref.ann[posns,c( "Hetero.ALT.reads","Hetero.REF.reads","Hetero.Read.Balance","Gene.Names","description","Amino_acids.Embl","FILTER","gerp.scores")]
 ## if(sum(c("GENO.aff.Seq.Filt","GENO.unaff.Seq.Filt","GENO.aff.Seq","GENO.unaff.Seq") %in% colnames(a.sig))==4){
 ##   a.sig<-

a.sig<-cbind(traits[i],a.sig,extra,stringsAsFactors=F)

if(i==1){
  all.sig<-a.sig
}else{
    all.sig<-rbind(all.sig,a.sig)
  }

} ###  loop over traits




options(scipen=4)

dim(all.sig)
all.sig[1:5,]

##########################################################################
######### set p_min correctly
test<-is.na(data[,"P_min"]) | is.na(as.numeric(data[,"P_min"]))
sum(test)

data[test,][1:5,]
min(as.numeric(all.sig[all.sig[,"SNP"]=="exm1000038","P"]),na.rm=TRUE)
lowest.p<-tapply(all.sig[,"P"],all.sig[,"SNP"],function(x) min(as.numeric(x),na.rm=TRUE))
lowest.p[!is.finite(lowest.p)]<-NA
lowest.p["exm1000038"]
lowest.p[1:40]
order.by<-order(as.numeric(lowest.p),decreasing=FALSE)
lowest.p<-lowest.p[order.by]

posns<-match(all.sig[,"SNP"],names(lowest.p))
missing<-is.na(posns)
sum(missing)
length(posns)
P_min<-lowest.p[posns]
P_min[1:20]
dim(all.sig)
cbind(all.sig[,"P"],P_min)[1:50,]

all.sig<-cbind(all.sig[,1:16],P_min,all.sig[,17:dim(all.sig)[2]])
######### set p_min correctly
data<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/first cut short list ori.Mhairi.ED.ED.modified.Mhairi.ED.Final.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
 plan<-tapply(data[,"Plan"],data[,"SNP"],function(x) paste(unique(x[!is.na(x) | x !=""]),collapse=";") )
plan[1:5]

plan[plan=="NA"]<-NA
plan<-gsub("^NA;","",plan)
plan<-gsub(";NA$","",plan)

posns<-match(all.sig[,"SNP"],names(plan))
missing<-is.na(posns)
sum(!missing)
length(posns)
Plan<-plan[posns]

Plan[1:50]
length(Plan)
cbind(all.sig,Plan)[1:50,]
all.sig<-all.sig[,-1*grep("Plan",colnames(all.sig))]


all.sig<-cbind(all.sig[,1:16],Plan,all.sig[,17:dim(all.sig)[2]])
all.sig[1:2,]

setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point")
write.table(all.sig,file="final.meta.analysis.Chip.centric.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(all.sig,file="final.meta.analysis.Chip.centric_final.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


gefos<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/SNP_lists_used_in_exomeChip/GEFOS SNPs in final file.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

gefos[,1] %in% all.sig[,"SNP"]

gefos.data<- all.sig[all.sig[,"SNP"] %in% gefos[,1] ,]

dim(gefos.data)

write.table(gefos.data,file="final.meta.analysis.Chip.centric.GEFOS.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


########################################################################## SPLIT UP BIG FILE
table(data$Plan)
wanted<-data$Plan=="evoker.good::altered" & !is.na(data$Plan)

test<-"UQDIexomechipWithAOGC_chr1:152081443:152081443:C:T-1_T_R_2110821056"
test<-"UQDIexomechipWithAOGC_chr1:12144552:12144552:A:C-1_T_F_2098199048"


wanted<-data$MarkerName==test & data$traits=="BMD_EFF_STD_HIP" & !is.na(data$MarkerName)


sum(wanted)
data[wanted,]
 data[wanted,c("MarkerName","traits","Plan","P.Chip")]

wanted<-all.sig$MarkerName==test & all.sig$traits=="BMD_EFF_STD_HIP" & !is.na(all.sig$MarkerName)
sum(wanted)
all.sig[wanted,]

wanted<-data$MarkerName==test & data$traits=="BMD_EFF_STD_HIP" & !is.na(data$MarkerName)



wanted.2<-all.sig$SNP=="UQDIexomechipWithAOGC_chr6:42574146:42574146:G:T-1_T_R_2098241709"
sum(wanted.2)
all.sig[wanted.2,]

1   UQDIexomechipWithAOGC_chr1:152081443:152081443:C:T-1_T_R_2110821056  152081443     6905     0.9065     0.2946    0.00137    3.077     0.002098 



plink --bed /media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.bed --bim /media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.bim --fam BMD_EFF_STD_HIP.fam --remove remove_from.BMD_EFF_STD_HIP --assoc --allow-no-sex --out BMD_EFF_STD_HIP./media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL --noweb


all.sig<-read.delim("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/final.meta.analysis.Chip.centric.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point")
dim(all.sig)

all.sig[1:15,1:5]

colnames(all.sig)[1]<-"traits"

table(all.sig[,"traits"])

subset<-c("BMD_EFF_STD_FN","BMD_EFF_STD_HIP","BMD_EFF_STD_LS")

subset<-c("BMD_EFF_STD_FN","BMD_EFF_STD_HIP")

subset<-c("BMD_EFF_STD_HIP")

wanted<-all.sig[,"traits"] %in% subset

getwd()

write.table(all.sig[wanted,],file="final.meta.analysis.Chip.centric.GEFOS_FN_LS_HIP.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(all.sig[wanted,],file="final.meta.analysis.Chip.centric.GEFOS_FN_HIP.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(all.sig[wanted,],file="final.meta.analysis.Chip.centric.GEFOS_HIP.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)








pleo@DI-LW-BRN-011:/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes$ grep exm1327856 BMD_EFF_STD_HIP.recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.qassoc
  17                                                             exm1327856   42225547     6886    0.07347    0.01636    0.00292     4.49    7.234e-06 
pleo@DI-LW-BRN-011:/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes$ grep exm1327856 recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.bim
17	exm1327856	0	42225547	C	A
pleo@DI-LW-BRN-011:/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes$ grep 42225547 recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.bim
17	exm1327856	0	42225547	C	A
pleo@DI-LW-BRN-011:/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes$ grep exm1327856 BMD_EFF_STD_HIP.recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.hwe
  17                                                             exm1327856  ALL(QT)    C    A        608/2819/3459   0.4094   0.4143       0.3228
pleo@DI-LW-BRN-011:/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes$ head -1 BMD_EFF_STD_HIP.recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.hwe
 CHR                                                                    SNP     TEST   A1   A2                 GENO   O(HET)   E(HET)            P 
pleo@DI-LW-BRN-011:/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes$ 



pleo@DI-LW-BRN-011:/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/AOGC_vcf_to_plink$ grep 42225547 BMD_EFF_STD_HIP.AOGC_vcf_final_chr17.qassoc
  17                                                                                  17:42225547:42225547:C:A:snp   42225547      971    0.04076    0.07717  0.0002878   0.5281       0.5975 
pleo@DI-LW-BRN-011:/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/AOGC_vcf_to_plink$ head -1 BMD_EFF_STD_HIP.AOGC_vcf_final_chr17.qassoc
 CHR                                                                                                           SNP         BP    NMISS       BETA         SE         R2        T            P 
pleo@DI-LW-BRN-011:/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/AOGC_vcf_to_plink$ grep 42225547 BMD_EFF_STD_HIP.AOGC_to_vcf_filtered_17.qassoc
  17                                                                                  17:42225547:42225547:C:A:snp   42225547      941    0.05705    0.07859  0.0005609    0.726        0.468 
pleo@DI-LW-BRN-011:/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/AOGC_vcf_to_plink$ grep 42225547 BMD_EFF_STD_HIP.AOGC_to_vcf_filtered_17.hwe
  17                                                                                  17:42225547:42225547:C:A:snp  ALL(QT)    C    A           77/402/462   0.4272   0.4163       0.4808
pleo@DI-LW-BRN-011:/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/AOGC_vcf_to_plink$ grep 42225547 BMD_AFFSTAT.AOGC_vcf_final_chr17.hwe 
  17                                                                                  17:42225547:42225547:C:A:snp      ALL    C    A           79/436/469   0.4431   0.4215       0.1128
  17                                                                                  17:42225547:42225547:C:A:snp      AFF    C    A           40/220/234   0.4453   0.4229        0.287
  17                                                                                  17:42225547:42225547:C:A:snp    UNAFF    C    A           39/216/235   0.4408     0.42       0.3326


######################################################
## fam.template.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.fam"
## #fam.template.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL_chr1.fam"
## #fam.template.file<-"/media/UQCCG/UQCCG-Projects/PAUL_LEO/AOGC exome chip core genotyping/keep.txt"
## annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"

## #recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL
## #the.chr<-basename(fam.template.file)
## ## the.fam.dir<-dirname(fam.template.file)
## ## files<-dir(the.fam.dir)  
## ## projects<-gsub(".bed","",files[grepl(".bed$",files)])  # files[grepl(".bed$",files)]
## setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/")
## projects<-"recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL"

## fam<-read.table(fam.template.file,header=F,fill=TRUE,stringsAsFactors=FALSE)
## ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

## related.to.remove<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.recode.txt",header=F,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
## related.to.remove<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.txt",header=F,fill=TRUE,sep=",",stringsAsFactors=FALSE)
## related.to.remove[1:5,]
## dim(related.to.remove)
## sum(related.to.remove[,1] %in% fam[,1]) # 0


## traits<-c("BMD_EFF_STD_HIP","BMD_EFF_STD_LS","BMD_EFF_STD_FN","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP","NEW_FR_18_NOT_TRIVIA_VS_NEVER_FX","NEW_FR_50_OP_vs_never_fx","NEW_FR_50_OP_vs_no.adult.fx","NEW_VERT_FX_OP_VS_NEVER_FX","NEW_VERT_FX_OP_VS_NO_ADULT_FX","NEW_NONVERT_OP_FX_50_VS_NEVER_FX","NEW_NONVERT_OP_FX_50_VS_NO_ADULT_FX","HIP_FR_50_EVER","new_FOREARM_FR_LOTRAUMA_VS_NEVER_FX","NEW_FOREARM_FR_LOTRAUMA_VS_NO_ADULT_FX","RECODED_EMCC_BP_GRP")

## #traits<-c("TOT_HIP_GCM","LS_GCM","FN_GCM2","BMD_AFFSTAT","EVER_FX","FR_18_NOT_TRIVIA","FR_50_OP","VERT_FX_OP","NONVERT_OP_FX_50","FOREARM_FR_LOTRAUMA","HIP_FR_50_OP","EMCC_BP_GRP")
## names(traits)<-rep("logistic",times=length(traits))
## names(traits)[1:3]<-rep("assoc",times=3)
## traits
## traits %in% colnames(ann)



code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")


media/UQCCG/GWAS/1000G_UK10K_common_genes/genome.pos.txt
imputable.bim<-read.table("/media/UQCCG/GWAS/AOGC_UK10K_1000G_imputed/AOGC_10kmerge_ALL.bim",header=F,fill=TRUE,stringsAsFactors=FALSE)
colnames(imputable.bim)<-c("chr","SNP","cm","pos","A1","A2")

imputable.r2<-read.table("/media/UQCCG/GWAS/AOGC_UK10K_1000G_imputed/AOGC_10kmerge_ALL_info.txt",header=T,fill=TRUE,stringsAsFactors=FALSE)
imputable.r2[1:5,]
imputable.bim[1:5,]

dups<-duplicated(imputable.bim[,"SNP"])
sum(dups)
imputable.bim<-imputable.bim[!dups,]

dups<-duplicated(imputable.r2[,"rs_id"])
sum(dups)
imputable.rs<-imputable.r2[!dups,]


posns<-match(imputable.bim[,"SNP"],imputable.r2[,"rs_id"])
missing<-is.na(posns)
sum(missing)

imputable.r2<-imputable.r2[posns,]

save(list=c("imputable.r2","imputable.bim"),file="/media/UQCCG/GWAS/AOGC_UK10K_1000G_imputed/imputable.RData")

load("/media/UQCCG/GWAS/AOGC_UK10K_1000G_imputed/imputable.RData")


final<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/Emma_2015_analysis/Yippee! finished!! from ED/replicate_AOGC_sequencing_final_unfiltered.csv",header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

final[1:5,]


load("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/all.meta.sig.RData")
colnames(all.meta)

#.snps<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/Unique SNPs to check LD.csv",header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

ld.snps<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/Emma_2015_analysis/Yippee! finished!! from ED/EDs final list for both!!!.csv",header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

ld.snps[1:5,]
ld.snps<-ld.snps[,1]
ld.snps.1<-unique(ld.snps)





ld.snps<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/Emma_2015_analysis/Yippee! finished!! from ED/SNPs with MAF_above_5%_check imputability.no duplicates.rsID included_FRACTURE.csv",header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)


ld.snps[1:5,]
ld.snps<-ld.snps[,1]
ld.snps<-unique(ld.snps)
ld.snps.1 %in% ld.snps
ld.snps<-unique(c(ld.snps,ld.snps.1)) ### come for testing


plan<-read.table("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/Emma_2015_analysis/Yippee! finished!! from ED/The final list_BMD SNPs with articulated plan.csv",header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
test<-unique(plan[,"SNP"])



length(test)
length(ld.snps)
test %in% ld.snps
sum(!test %in% ld.snps)
test<-test[(!test %in% ld.snps)]
test<-test[test!=""]
length(test)
colnames(all.meta)

freq[1:5,]
freq[freq[,"SNP"]%in% test,]

ld.snps<-c(ld.snps,unique(test))

colnames(plan)



plan[plan[,"SNP"] %in% test, c("refGene..type"  ,"rs.id","Comment.in.December","CONSIDER.REPLICATION","Summary.comment")]
all.meta[1:5,]

posns<-match(ld.snps,all.meta[, "MarkerName"])
missing<-is.na(posns)
sum(missing)


replicate<-cbind(all.meta[posns[!missing],c( "MarkerName", "rs.id","chr","start","REF","ALT","refGene..type","Gene.Names","description","gerp.scores","PolyPhen.desc","PolyPhen.scores","GENO.Chip","P.Chip")])   
replicate[1:5,]
## ld.snps<-c(a.data[,"SNP"],"HLA_B_1501","HLA_C_03","HLA_C_0303", "HLA_C_0702","HLA_B_44","HLA_B_4402","HLA_A_03","HLA_A_0301","HLA_DQB1_0603","HLA_DQA1_0103","HLA_DQB1_0602","HLA_DQA1_0102", "HLA_DRB1_1501","HLA_DRB1_0401", "HLA_DRB1_1301")

# files<-paste(paste( "LD.with2",ld.snps,sep="."),"ld",sep=".") # gies with cases only below 
files<-paste(paste( "LD.with",ld.snps,sep="."),"ld",sep=".")
#setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/") ## AA regression
files

## data.file<- "/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/META.SIG.SUMMARY.txt"
## data<-read.table(data.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

###################################
chk.snps<-c("rs2523765","rs190921232")
all.meta[all.meta[,"rs.id"] %in% chk.snps,1:29]
chk.snps<-unique(all.meta[all.meta[,"rs.id"] %in% chk.snps,"MarkerName"])
chk.snps %in% ld.snps


#########################






annotation.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/AOGC_HBM_ALL_PHENOTYPES_RESIDUALS_UPDATED FX OPTIONS.txt"
ann<-read.table(annotation.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
outlyers.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/related.to.remove.recode.txt"


fam.ori.file<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.fam"
the.bed<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.bed"
the.bim<-"/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.bim"
#the.dosage<-"/media/scratch/cervical_cancer/snp2hla/genotyped_common_snp/association/data/cervical_cancer_common_snp2hla_r2.dosage"
the.fam<-fam.ori.file
  to.remove<-outlyers.file
 setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/Emma_2015_analysis/Yippee! finished!! from ED")

i<-1

system(      paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam," --remove ", outlyers.file,"--allow-no-sex --freq" ,"--out","FREQ" ,"--noweb",sep=" ") ) 
freq<-read.table("FREQ.frq",header=T,fill=TRUE,stringsAsFactors=FALSE)
freq[1:5,]

setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/Emma_2015_analysis/Yippee! finished!! from ED/LD_files")

for(i in 1:length(ld.snps)){

## system(      paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam,"--filter-cases  --remove ", outlyers.file," --r2 dprime --ld-snp ",ld.snps[i],"--allow-no-sex --ld-window-r2 0  --ld-window-kb 500000  --ld-window 99999 ","--hide-covar","--out",paste( "LD.with",ld.snps[i],sep=".") ,"--noweb",sep=" ") ) # no ld with dosage


system(      paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam, "--r2 --ld-snp ",ld.snps[i],"--allow-no-sex --ld-window-r2 0  --ld-window-kb 500000  --ld-window 99999 ","--hide-covar","--out",paste( "LD.with",ld.snps[i],sep=".") ,"--noweb",sep=" ") ) # no ld with dosage


}

## plink --bed /media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.bed --bim /media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.bim --fam /media/UQCCG/UQCCG-Projects/AOGC_exome_chip/working_genotypes/recode_sampleNames_zCall_AOGC.with.evoker_corrected_clean_FINAL.fam --r2 --ld-snp  exm1131025 --allow-no-sex --ld-window-r2 0  --ld-window-kb 500000  --ld-window 99999  --out LD.with.exm1131025 --noweb


extra<-"exm-rs2239806"
system(      paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam, "--r2 --ld-snp ",extra,"--allow-no-sex --ld-window-r2 0  --ld-window-kb 500000  --ld-window 99999 ","--hide-covar","--out",paste( "LD.with",extra,sep=".") ,"--noweb",sep=" ") )

####################### choose to use R2 or Dprime
#use.r2.DP.col<-"DP" ;use.r2.DP.col.name="Dprime"

use.r2.DP.col<-"R2" ;use.r2.DP.col.name="R2"
r2.thresh<-0.7


############################# RESTART HERE IF HAVE FILES ALREADY

files<-paste( "LD.with",ld.snps,"ld",sep=".")

i<-1
print(length(files))
for (i in 1:length(files)){
  print(length(files)-i)
 data<-read.table(files[i],header=T,fill=TRUE,stringsAsFactors=FALSE)
 all.data<-data
 all.data

colnames(all.data)
order.by<-order(as.numeric(all.data[,use.r2.DP.col]),decreasing=TRUE)
all.data<-all.data[order.by,]
#print(all.data[all.data[,"SNP_B"] %in% ld.snps,])
all.data[all.data[,"SNP_B"] %in% ld.snps,]
  if(i==1){
    ld.data<-all.data[all.data[,"SNP_B"] %in% ld.snps,]
  }else{
    ld.data<-rbind(ld.data,all.data[all.data[,"SNP_B"] %in% ld.snps,])
  }

}

 freq[1:5,] 
ld.data[1:5,]

posns<-match(ld.data[,"SNP_A"],freq[,"SNP"])

MAF<-freq[posns,"MAF"]
ld.data<-cbind(ld.data,MAF)

ld.data[1:15,]

## table(data[,"traits"])
## posns<-match(ld.data[,"SNP_A"],data[,"SNP"])

## P<-data[posns,"P.Chip"]
## ld.data<-cbind(ld.data,MAF)

#### remove identical
ld<-ld.data[ld.data[,"SNP_A"]!=ld.data[,"SNP_B"],]
ld[1:15,]
dim(ld)
max.ld<-tapply(ld[,"R2"],ld[,"SNP_A"],max)
max.ld[1:5]

sort(max.ld)
snps.in.ld.ori<-max.ld[max.ld>=r2.thresh]

snps.not.in.ld<-max.ld[max.ld<r2.thresh]

length(snps.in.ld.ori)
length(snps.not.in.ld)
length(max.ld)
sort(snps.in.ld.ori)[1:10]
targets<-"exm-rs2596560"
targets<-to.remove
ld[ld[,"SNP_A"] %in% targets,]
ld[ld[,"SNP_A"] %in% chk.snps,]
#ld.data["AA_DRB1_11_32660115_PV" ,"AA_DRB1_96_32657590_QY","AA_DRB1_96_32657590_HE")
###########################################################################################################


################# recursively remove the sample that appear more than once with the lowerest call rate till all only singles
snps.in.ld<-snps.in.ld.ori

length(snps.in.ld)
ld[1:5,]
ld.trim<-ld[ ( ld[,"SNP_A"] %in% names(snps.in.ld.ori) | ld[,"SNP_B"] %in% names(snps.in.ld.ori))  &  as.numeric(as.character(ld[,"R2"]))>=r2.thresh ,]
#ld.trim<-ld[ld[,"SNP_B"] %in% names(snps.in.ld.ori),]
ld.trim[1:5,]
dim(ld.trim)
to.keep.snps<-{}
to.exclude.snps<-{}
while((dim(ld.trim)[1]!=0)){


order.by<-order(as.numeric(ld.trim[,"MAF"]),decreasing=TRUE)
ld.trim<-ld.trim[order.by,]
ld.trim[1:5,]  
to.remove<-ld.trim[1,"SNP_A"]
to.remove
to.keep.snps<-c(to.keep.snps,to.remove)
to.keep.snps

targets<-to.remove
to.exclude<-ld.trim[ld.trim[,"SNP_A"] %in% targets  ,]
to.exclude
to.exclude<-to.exclude[as.numeric(to.exclude[,"R2"])>=r2.thresh,]
to.exclude
#### find out which column ha the most entries
to.exclude.snps<-c(to.exclude.snps,unique(to.exclude[,"SNP_B"]))
to.exclude.snps
done<- ld.trim[,"SNP_A"] %in% c(to.keep.snps,to.exclude.snps)
print("--------------------------------")
print(dim(ld.trim))
print(sum(done))
ld.trim<-ld.trim[!done,] ## remove new target and the SNps in LD with it from consideration
print(dim(ld.trim))

                       }

to.keep.snps
length(to.keep.snps)
#"UQDIexomechipWithAOGC_chr8:120008371:120008371:T:C-1_B_F_2098194396" "UQDIexomechipWithAOGC_chr7:96116813:96116813:T:C-1_B_F_2098189130"   "exm-rs718314"                                                        "exm-rs1812175" 

length(to.exclude.snps)
to.exclude.snps

# "exm-rs11995824"                                                     "exm-rs6993813"                                                      "UQDIexomechipWithAOGC_chr7:96119481:96119481:G:A-1_T_F_2098189136"  "UQDIexomechipWithAOGC_chr12:26470850:26470850:T:A-1_T_F_2098196022" "exm-rs7689420

replicate[1:5,]
freq[1:5,]


posns<-match(replicate[, "MarkerName"], freq[,"SNP"])
missing<-is.na(posns)
sum(missing)

MAF<-freq[posns,"MAF"]
LD.reject<-replicate[, "MarkerName"] %in% to.exclude.snps
LD.reject

imputable.bim[1:5,]

key<-build.key(replicate,c("chr","start"))
key.imput<-build.key(imputable.bim,c("chr","pos"))


key.imput[1:5]
key[1:5]

posns<-match(key,key.imput)
missing<-is.na(posns)
sum(missing)
sum(!missing)
## can.impute<-unique(replicate[!missing], "MarkerName"])
## can.impute

imputable.r2[1:5,]

Can.impute<-as.numeric(imputable.r2[posns,"info"])>=0.8

replicate<-cbind(replicate,MAF,LD.reject,Can.impute,imputable.bim[posns,],imputable.r2[posns,]) # replicate<-cbind(replicate,MAF,LD.reject,imputable.bim[posns,],imputable.r2[posns,])
replicate[1:5,]




no.result<-is.na(replicate[,"SNP"])
sum(no.result)
posns<-match(replicate[no.result, "rs.id"],imputable.bim[,"SNP"])
missing<-is.na(posns)
sum(missing)
sum(!missing)

replicate.extra<-cbind(replicate[no.result,],MAF[no.result],LD.reject[no.result],imputable.bim[posns,],imputable.r2[posns,])
colnames(replicate.extra)<-colnames(replicate)
replicate.extra

final[1:5,1:5]
final[,"MarkerName"] %in% replicate[,"MarkerName"]

replicate[,"MarkerName"] %in% final[,"MarkerName"]

setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/Emma_2015_analysis")
write.table(replicate,file="replicate_FX_snps_not_in_0.7_LD.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)


setwd("/media/UQCCG/UQCCG-Projects/AOGC_exome_chip/meta_analysis_single_point/Emma_2015_analysis/Yippee! finished!! from ED")

write.table(replicate,file="replicate_AOGC_sequencing.ALL.extra.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

write.table(replicate,file="replicate_AOGC_sequencing.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(replicate,file="replicate_AOGC_sequencing.high.MAF.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

replicate.1[(replicate.1[,"MarkerName"] %in% replicate[,"MarkerName"]),1:5]

# replicate.1<-replicate

######### now remove half of the singles chosing ones with the worst call rate. 
######### related now only contains unique samples
related[1:8,]

ld[1:8,]

posns1<-match(related[,"ID1"],ld[,"IID"])
posns2<-match(related[,"ID2"],ld[,"IID"])


#to.swap<-as.numeric(related[,"Call rate ID1"]) > as.numeric(related[,"Call rate ID2"])
to.swap<-as.numeric(ld[posns1,"call.rate"]) > as.numeric(ld[posns2,"call.rate"]) 

#cbind(related,ld[posns1,"call.rate"],as.numeric(ld [posns2,"call.rate"]),to.swap)[1:10,]

related[to.swap,"ID1"]<-related[to.swap,"ID1"]

to.remove.samples<-c(to.remove.samples,related[,"ID1"])
length(to.remove.samples)

write.table(cbind(to.remove.samples,to.remove.samples),file="related.to.remove.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
length(to.remove.samples)
################################################# Check frequency for calling algorithms
