


# assume assoc has columns
if(!exists("target.pval.col")){
target.pval.col<-"P"
}

if(!exists("the.SNP.col.name")){
the.SNP.col.name<- "SNP"
}

if(!exists("the.chr.col.name")){
the.chr.col.name<- "CHR"
}

if(!exists("the.pos.col.name")){
the.pos.col.name<- "BP"
}



if(!exists("work.dir")){
    work.dir<-dirname(file.assoc)
}
if(!exists("plot.dir")){
    plot.dir<-dirname(file.assoc)
}

if(!exists("rates.dir")){
    if(genome.build=="hg18"){
    rates.dir<-paste(path.to.UQCCG,"/GWAS/Recombination_Rates/Recombination Rates_hg18",sep="")
    rate.file.prefix<-file<-"genetic_map_chr"
    rate.file.suffix<-"_b36.txt"
    rate.file.columns<-c("position","rate","map")
} #hg18
           if(genome.build=="hg19"){
    rates.dir<-paste(path.to.UQCCG,"/GWAS/Recombination_Rates/genetic_map_HapMapII_GRCh37",sep="")
    rate.file.prefix<-file<-"genetic_map_GRCh37_chr"
    rate.file.suffix<-".txt"
    rate.file.columns<-c("chr","position","rate","map")
} # hg19
}


####################################### read in data

    
###########################################

setwd(work.dir)
imput<-read.delim(file.assoc,header=T,skip=0,sep="",fill=TRUE,stringsAsFactors=FALSE)
imput[1:5,]
  ### needed for point types later

colnames(imput)[colnames(imput)==the.SNP.col.name]<-"SNP"
imput[,"SNP"]<-gsub("^RS","rs",as.character(imput[,"SNP"]))
rownames(imput)<-as.character(imput[,"SNP"])

colnames(imput)[colnames(imput)==target.pval.col]<-"Pval"
colnames(imput)[colnames(imput)==the.pos.col.name]<-"position"
colnames(imput)[colnames(imput)==the.chr.col.name]<-"chrom"

imput[is.na(imput[,"Pval"]),"Pval"]<-1 # for conditional regression
imput[,"Pval"] <- -log10(imput[,"Pval"] )
imput[1:5,]

if(exists("bim.file")){
snps<-read.delim(bim.file.genotyped,header=F,skip=0,sep="",fill=TRUE,stringsAsFactors=FALSE)
snps<-snps[,2]
snps<-gsub("^RS","rs",snps)
}else{
    snps<-as.character(imput[,"SNP"])
}


imput.ori<-imput

#snps<-rownames(imput)[imput[,"Genotype.all"]=="Genotyped"] ## genotyped SNPS
#discovery<-discovery.ori
############### subset imput to interested range
the.range<-(imput.ori[,"chrom"]==the.chr) & (as.numeric(imput.ori[,"position"])  >= low.cut) & (as.numeric(imput.ori[,"position"]) <= high.cut)
imput<-imput[the.range,]
dim(imput)
if(is.null(dim(imput))){print("ERROR NO genotypes left in desited range")}
############ restrict range #######
imput<-imput[!is.na(imput[,"Pval"]),]
max(imput[,"position"])
min(imput[,"position"])

if(!exists("the.snp")){the.snp<-imput[which.max(imput[,"Pval"]),"SNP"]}
if(the.snp==""){the.snp<-imput[which.max(imput[,"Pval"]),"SNP"]}

#imput[imput[,"Pval"]==max(imput[,"Pval"]),]

ymax.range<-max(imput[,"Pval"])+1.0
pch.array<-rep(23,times=dim(imput)[1])
posns<-match(rownames(imput),snps)
genotyped<-!is.na(posns)
pch.array[genotyped]<-19


####################################### get ld information

if(got.ld){
##################################################### NO LD imformaton
setwd(work.dir)

extra<-try(read.delim(file.plink.ld,header=T,skip=0,sep="",fill=TRUE,stringsAsFactors=FALSE),silent=TRUE)
if(inherits(extra, "try-error")){extra<-{}}else{
extra[extra[,"SNP_B"]==the.snp,]
dim(extra)
}
 ## LD form 1000 genomes but does not contain positio info
extra2<-try(read.delim(file.1000G.ld,header=T,skip=0,sep="",fill=TRUE,stringsAsFactors=FALSE),silent=TRUE)

if(!inherits(extra2, "try-error")){
  
extra[1:5,]
extra2[1:5,]

# tapply(extra2[,"SNP_A"],extra2[,"SNP_A"],length)
## posns<-match(extra2[,"SNP_B"],extra[,"SNP_B"])
## cbind(extra[posns,],extra2)[1:5,]

found<-extra2[,"SNP_B"] %in% extra[,"SNP_B"]
sum(!found)

## extra4<-merge(extra,extra2,by.x="SNP_B",by.y="SNP_B",all=FALSE,all.y=TRUE,all.x=TRUE) # merges columnwise like cbin common stuff
extra<-rbind(extra,extra2[!found,])

}
}else{ #no.ld go get sfrom scratch assimung a bed/bim/fam is in the working directory

    the.bim<-bim.file
    
    if(!exists("the.bed")){
        the.bed<-gsub(".bim$",".bed",bim.file)
    }

        if(!exists("the.fam")){
        the.fam<-gsub(".bim$",".fam",bim.file)
    }

    
system(      paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam, "--r2 --ld-snp ",the.snp,"--allow-no-sex --ld-window-r2 0  --ld-window-kb 500000  --ld-window 99999 ","--hide-covar","--out",paste( "LD.with",the.snp,sep=".") ,"--noweb",sep=" ") )

    
extra<-read.table(paste( "LD.with",the.snp,"ld",sep="."),header=T,fill=TRUE,stringsAsFactors=FALSE)
extra[1:5,]

} # no.ld

 ################################################### got ld in array extra
    
       ################################################### get rates

  #################### read recombination file downloaded from   http://ftp.hapmap.org/recombination/latest/rates/

                        
  setwd(rates.dir)
options(show.error.messages = TRUE)
file<-paste(rate.file.prefix,the.chr, rate.file.suffix,sep="")
  file
recomb<-try(scan(file,what=character(length(rate.file.columns)),skip=1,sep="\t",fill=TRUE))
num.lines<-length(recomb)/length(rate.file.columns)
dim(recomb)<-c(length(rate.file.columns),num.lines)
recomb<-t(recomb)
colnames(recomb)<- rate.file.columns
recomb[1:5,]
the.range<-(as.numeric(recomb[,"position"])  >= low.cut) & (as.numeric(recomb[,"position"]) <= high.cut)
recomb<-recomb[the.range,]

recomb.low<-0

 setwd(plot.dir)




###### imput, and extra may have different number of SNPs as the coveresion to the impuated data to a ped file may have caused extra missing snps
##### get only the common snps
common.snps<-intersect(rownames(imput),extra[,"SNP_B"])
length(common.snps)
dim(extra)
dim(imput)
#reorder common so in accenting order
## posns<-match(common.snps,extra[,"SNP_B"])
## extra<-extra[posns[!is.na(posns)],]
## the.order<-order(extra[,"BP_B"])
## extra<-extra[the.order,]

imput[,"R2"]<-1.1
posns<-match(rownames(imput),extra[,"SNP_B"],)
missing<-is.na(posns)
sum(missing)
lost<-imput[setdiff(1:dim(imput)[1],posns[!is.na(posns)]),]
imput[!missing,"R2"]<-as.numeric(extra[posns[!missing],"R2"])

## dim(lost)
## sum(!missing)
## imput[!missing,][1:5,]

## the.order<-order(imput[!missing,][,"Pval"],decreasing=TRUE)
## imput![missing,][the.order,][1:5,]

imput[missing,"R2"]<-0
#imput<-imput[!missing,]

# must be zero
##################### fix most associated SNP:
imput[the.snp,] ## must be present else lost in intersection with extra
if(imput[the.snp,"R2"]!=1){imput[the.snp,"R2"]<-1; print("WARNING")} # should be one anyway!

###########################################################################################################
##########################################################################################################





############################################################################### ONCE ###################
############ for archive version get host name from ensemble -> BioMart-> archive


if(genome.build=="hg18"){mart=useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="may2009.archive.ensembl.org",path="/biomart/martservice",archive=FALSE)
    } else if(genome.build=="hg19"){
      mart <-useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
      print("using hg19")
    } else if(genome.build=="mm9"){
      mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl");print("here:")
     } else if(genome.build=="hg20"){
      mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")                                 
     } else {print("ERROR genome build not found cant call BIOMART")}# is alread have mart then skip
  
## mart



## listMarts(host="may2009.archive.ensembl.org",path="/biomart/martservice",archive=FALSE)
##   #http://www.ensembl.org/info/website/archives/assembly.html  may 2009 last NCBI 36 - reference 
## listAttributes(mart)[1:200,]
 
#mart<-  useMart("ensembl",dataset="hsapiens_gene_ensembl")"
####################################################################################################

## mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")


a.filter<-c( "chromosome_name", "start" , "end")
fil.vals<-list(as.character(the.chr), low.cut, high.cut)
fil.vals
### below for latest version
if(genome.build=="hg18"){
    exons<-getBM(attributes = c("ensembl_gene_id","external_gene_id","5_utr_start","5_utr_end","ensembl_exon_id","exon_chrom_start","exon_chrom_end","3_utr_start","3_utr_end","rank","strand","gene_biotype","chromosome_name","start_position","end_position"), filters = a.filter, values=fil.vals, mart = mart) ## works for old mart
}

if(genome.build=="hg19"){                        
 exons<-getBM(attributes = c("ensembl_gene_id","5_utr_start","5_utr_end","ensembl_exon_id","exon_chrom_start","exon_chrom_end","3_utr_start","3_utr_end","rank","strand","gene_biotype","chromosome_name","start_position","end_position"), filters = a.filter, values=fil.vals, mart = mart)
}

## fil = listFilters(mart)
## fil[grep("exon",fil[,1]),]

unique.genes<-unique(exons[,1])
unique.genes

colors.ori<-c("darkblue","green","red","purple","forestgreen","gold3","magenta","orange","black","cyan","salmon","rosybrown4","plum4","orchid2","orangered3","olivedrab2","indianred","grey63","brown","aquamarine3","seagreen2")

colors.ori.expanded<-rep(colors.ori,times=length(unique.genes))# longer than it need to be
col.array<-colors.ori.expanded[1:length(unique.genes)]
names(col.array)<-unique.genes
unique(exons[,"gene_biotype"])
keep<-exons[,"gene_biotype"] %in% c("protein_coding","miRNA")
exons<-exons[keep,]
## plus<-getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","5_utr_start","3_utr_start","ensembl_exon_id","exon_chrom_start","exon_chrom_end","rank","strand","gene_biotype"), filters = a.filter, values=fil.vals, mart = mart)


fil.vals<-list(as.character(the.chr),  low.cut, high.cut)
ann<-getBM(attributes = c( "ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position","strand","hgnc_symbol","gene_biotype"), filters = a.filter, values=fil.vals, mart = mart)
ann
keep<-ann[,"gene_biotype"] %in% c("protein_coding","miRNA")
if(sum(keep)>1){                                                   
ann<-ann[keep,]
}
# ann

the.order<-order(exons[,"ensembl_gene_id"],exons[,"rank"])
exons<-exons[the.order,]

