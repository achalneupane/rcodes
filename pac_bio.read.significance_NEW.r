
#############################################
## 1) Read in data at SNP location
## 2) estimate backgroud substitution error rare for "non SNP" alleles for a given SNP
## 3) For each well apply binomial approximation and obtain a p-value for null hypothesis  for SNP bases (2 per well) 
## "LPH-001-26_SCC_lbc33_TP53_Ex6_F_chr17_7577528_7577548.txt" no snp name
## [1] "LPH-001-26_SCC_lbc33_TP53_Ex6_F_chr17_7577528_7577548.txt"
## [1] "ERROR file conatins NO SNP name"
## [1] "LPH-001-26_SCC_lbc33_TP53_Ex6_F_chr17_7577528_7577548.txt"
## [1] "LPH-001-26_SCC_lbc33_TP53_Ex6_F_chr17_7577529_7577549.txt"
## [1] "ERROR file conatins NO SNP name"
## [1] "LPH-001-26_SCC_lbc33_TP53_Ex6_F_chr17_7577529_7577549.txt"
## [1] "LPH-001-26_SCC_lbc33_TRHDE_Ex13_F_chr12_73012761_73012781.txt"
## [1] "No DATA"

##################################################################################################################
## 1) Read in data at SNP location
#args <- commandArgs(TRUE)


############# directory of pileups 
#dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-03_LeoPharma_NovoAlign/Felicity_analysis/pileup_zmw/"
dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-03_LeoPharma_NovoAlign/Felicity_analysis/Pacbio_pileup/pileup_zmw_filter/"
dir<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-03_LeoPharma_NovoAlign/Felicity_analysis/Pacbio_pileup_version2/pileup_zmw_filter"
args<-dir(dir)
args<-args[!grepl(".Pac.bio.call.txt",args)]
args<-args[grepl("^LPH",args)]
#args<-args[!grepl("_SCC_",args) & !grepl("_PD_",args) & !grepl("AK1_",args) & !grepl("AK2_",args) ]


args[1:40]
## args<-"LPH-001-26_SCC_lbc33_DCLK1_Ex2_F_chr13_36700087_36700107.txt"
## args<-"LPH-001-14_AK1_lbc58_BCL2L12_Ex1_F_chr19_50169121_50169141.txt"
## args<-args[grep("BCL2L12",args)]


ifiles<-1
args[1:100]
ifiles<-53
length(args)
#args[35:40]


new.cells<-c("m151216_051615_42229_c100935702550000001823204005251625","m151216_093800_42229_c100935702550000001823204005251626")
do.cells<-new.cells


ref.posn<-4 ## posn of REF in (chr:start:end:REF:ALT:type) string
alt.posn<-5 ## posn of REF in (chr:start:end:REF:ALT:type) string
sub.read.threshold<-10 # don't consider celles with fewer sub-reads
min.alt.reads.dups<- 5  # in a deteroduplex require a more than this number of reads for the ALT allele

filter.cells<-FALSE
all.calls<-{}
for(ifiles in 1:length(args)){  


file<-args[ifiles]
print(file)
#file<-"pileup_zmw/LPH-001-3_SCC_lbc56_EBNA1BP2_Ex1_F_chr1_43637970_43637990.txt"
pile.up.file<-paste(dir,file,sep="/")

data<-read.delim(pile.up.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

if(dim(data)[1] ==0){
    print("No DATA")
    next
}

setwd(dirname(pile.up.file))
basename<-basename(pile.up.file)
basename<-gsub(".txt","",basename)

dim(data)
#data[1:5,]

##############



if(filter.cells){
    keep<-rep(FALSE,times=dim(data)[1])
    for(icell in 1:length(do.cells)){
        a.keep<-grepl(do.cells[icell],data[,"read_name"])
        print(length(a.keep))
        keep<-keep | a.keep
    }
    data<-data[keep,]
} # filter.cells
    

##### get alleles for this location from "snp" column in format (chr:start:end:REF:ALT:type)



if(sum(is.na(data[,"snp"]))==dim(data)[1]){
    print("ERROR file conatins NO SNP name")
    print(file)
    next
}

alleles<-strsplit(data[,"snp"],split=":")
REFs<-unlist(lapply(alleles,function(x) {x[ref.posn] } ))
ALTs<-unlist(lapply(alleles,function(x) {x[alt.posn] } ))

#### tridy UP AND CHECKS
### check have a ref and an alt always
non.snp.locations<-is.na(REFs) | is.na(ALTs)
print(paste("non.snp.locations:", sum(non.snp.locations)))
#print(data[bad.locations,])
data<-cbind(REFs,ALTs,data,stringsAsFactors=FALSE)
data[is.na(data[,"REFs"]),"REFs"]<-data[is.na(data[,"REFs"]),"ref_base"]

## posn of REF in (chr:start:end:REF:ALT:type) string

data.other<-data[non.snp.locations,] ## this pile up data beside the SNP location (can use to get backgroud error rate)
data<-data[!non.snp.locations,] ## this pile up data AT the SNP location

#data[1:5,]



###################### Get
#############################################
## 2) estimate backgroud substitution error rare form "non SNP" alleles for a given SNP across all wells
##        estimate backgroud substitution error rare form "non REF" alleles for a given SNP (assumes no SNP beside that location)
########### LOOP over SNPs (this could be over all samples, but can just use an individals)


################ get extimate of the backgroud error rate for a snps 
bases<-c("A","C","G","T")
base.cols<-paste(bases,"total",sep="_")  ### use these to get base count these exclude deletions.
the.snps<-unique(data[,"snp"])

pval.columns<-c("error.rate","total.bases",paste(base.cols,"pval",sep="."))
snp.error<-matrix(data=NA,nrow=length(the.snps),ncol=length(pval.columns))
colnames(snp.error)<-pval.columns
rownames(snp.error)<-the.snps


#### for( isample= 1:length(samples)) {  ### is have data for all samples i one file
i<-1
for (i in 1:length(the.snps)){
    data.sub <-data[data[,"snp"]==the.snps[i],]
    allele.counts<-apply(data.sub[,base.cols],2,function(x) sum(x,na.rm=TRUE))
    total.counts<-sum(allele.counts)
   ##  A_total C_total G_total T_total 
   ## 3123      75    6060      74
   #data.sub[1:5,]
    ### get backgroup subsitution error rate by excluing bases when snps is expected.
    snp.bases<-unique(c(as.character(data.sub[,"REFs"]),as.character(data.sub[,"ALTs"])))
    background.bases<-base.cols[!(base.cols %in% paste(snp.bases,"total",sep="_"))]
    background.error<-sum(allele.counts[background.bases])/total.counts
    snp.error[the.snps[i],"error.rate"]<-background.error ### pval.data same size as data keep error rate for that snp
    snp.error[the.snps[i],"total.bases"]<-total.counts
}


snp.error  ## error rate for SNp based on all wells at SNp location (excluding teh ref and alt alleles

#data[1:5,]
dim(data)


#data.other[1:5,]
################ get extimate of the backgroud error rate for a snps 



dim(data.other)[1]/length(unique(data.other[,"read_name"]))  ## this is the number of non snp locations 

#### for( isample= 1:length(samples)) {  ### is have data for all samples i one file
#have.reads<-data.other[,"totalCount"] > sub.read.threshold
#sum(have.reads)
#sum(!have.reads)
#unique.cells<-unique(data.other[have.reads,"read_name"])
unique.cells<-unique(data.other[,"read_name"]) 

pval.columns<-c("error.rate","total.bases","background.bases",paste(base.cols,"pval",sep="."))
pval.extra.other<-matrix(data=NA,nrow=length(unique.cells),ncol=length(pval.columns))
colnames(pval.extra.other)<-pval.columns
rownames(pval.extra.other)<-unique.cells


#pval.extra.other[1:5,]
ic<-2
for (ic in 1:length(unique.cells)){
the.cell<-unique.cells[ic]
the.cell
data.cell<-data.other[data.other[,"read_name"] == the.cell,]
# data.cell<-data.cell[data.cell[,"totalCount"] > sub.read.threshold,]
# sub.read.threshold
# if(dim(data.cell)[1]==0){next} # not enough sub-reads in that well- skip
dim(data.cell)
#data.cell[1:15,]


######## loop over all base types
    allele.counts.cell<-apply(data.cell[,base.cols],2,function(x) sum(x,na.rm=TRUE))
    total.counts.cell<-sum(allele.counts.cell)

#data.cell[1:15,c("REFs",base.cols)]

ref.base.posns<-paste(data.cell[,"REFs"],"total",sep="_")
for (i in 1:length(ref.base.posns)){
    data.cell[i,ref.base.posns[i]]<-0
}


    allele.counts.cell.nr<-apply(data.cell[,base.cols],2,function(x) sum(x,na.rm=TRUE))
    allele.counts.cell.nr
    background.bases.cell<-sum(allele.counts.cell.nr)

  
    background.error<-background.bases.cell/total.counts.cell
    if(total.counts.cell<sub.read.threshold){background.error<-0} # does get erroe rate form well with little coverage
    pval.extra.other[the.cell,"error.rate"]<-background.error
    pval.extra.other[the.cell,"total.bases"]<-total.counts.cell
    pval.extra.other[the.cell,"background.bases"]<-background.bases.cell ### pval.data same size as data keep error rate for that snp
} # loop over cells
#pval.extra.other[1:5,]
#hist(pval.extra.other[,"error.rate"])

#pval.extra.other[is.na(pval.extra.other[,"error.rate"]),][1:10,]


max(pval.extra.other[,"error.rate"],na.rm=FALSE)
mean(pval.extra.other[,"error.rate"],na.rm=FALSE)
median(pval.extra.other[,"error.rate"],na.rm=FALSE)

pval.extra.other.sum<-sum(pval.extra.other[,"background.bases"])/sum(pval.extra.other[,"total.bases"])
pval.extra.sum<-snp.error[,"error.rate"]
max.error<-max(pval.extra.other[,"error.rate"],na.rm=TRUE)
sd.error<-sd(pval.extra.other[,"error.rate"],na.rm=TRUE)
pval.extra.sum
pval.extra.other.sum
max.error
sd.error


###################################################################
############################ build array to hold p-values #########

unique.cells<-unique(data[,"read_name"])
base.p.labels<-paste(base.cols,"pval",sep=".")
pval.columns<-c("error.rate","total.bases","background.bases",base.p.labels)
pval.extra<-matrix(data=NA,nrow=length(unique.cells),ncol=length(pval.columns))
colnames(pval.extra)<-pval.columns
rownames(pval.extra)<-unique.cells


sum(rownames(pval.extra.other) %in% rownames(pval.extra))
sum(!(rownames(pval.extra.other) %in% rownames(pval.extra)))


posns<-match(rownames(pval.extra),rownames(pval.extra.other))
missing<-is.na(posns)
sum(missing)
pval.extra.other<-pval.extra.other[posns,]

posns<-match(rownames(pval.extra),data[,"read_name"])
missing<-is.na(posns)
sum(missing)
data<-data[posns,]

######################data and other roor matrics in the same order
## pval.extra.other[1:5,]
## pval.extra[1:5,]

print(paste("Error rate comparison form SNP location or, SNP region",pval.extra.sum,pval.extra.other.sum ))

##################################################################################################################
##################################################################################################################
backgroud.error<-max(pval.extra.sum,pval.extra.other.sum)


allele.counts<-apply(data[,base.cols],1,function(x) sum(x,na.rm=TRUE))
if(!is.finite(backgroud.error)){backgroud.error<-0.02} # fore cases when this are no reads this just fille the location
if(backgroud.error==0){backgroud.error<-0.02}

pval.extra[,"error.rate"]<-backgroud.error
pval.extra[,"total.bases"]<-allele.counts
## 3) For each well apply binomial approximation and obtain a p-value for null hypothesis  each well ... do for all bases though only need ref and alt
##

#data[1:5,]
#pval.extra[1:5,]


i<-3
for (i in 1:length(base.cols)){

#p.base<-dbinom(data[,base.cols[i]], pval.extra[,"total.bases"], pval.extra[,"error.rate"], log= FALSE) # original
p.base<-pbinom(data[,base.cols[i]], pval.extra[,"total.bases"], pval.extra[,"error.rate"], lower.tail = FALSE, log= FALSE) 
p.base[data[,base.cols[i]]==0]<-1
#pbinom(measure, n, p, lower.tail = TRUE, log.p = FALSE) + pbinom(measure, n, p, lower.tail = FALSE, log.p = FALSE)

 pval.extra[,paste(base.cols[i],"pval",sep=".")]<-signif(p.base,digits=7)
}

# pval.extra[1:5,]

#data[1:5,]

pval.thresh<-2*(pnorm(abs(c(3)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)) # 3sd
# pval.thresh<-1e-6 # original
pval.thresh

boolean.test<-subset(pval.extra,select=base.p.labels)
boolean.test<-boolean.test<=pval.thresh
colnames(boolean.test)<-gsub(".pval$",".count",colnames(boolean.test))
#boolean.test[1:5,]
alleles.called<-apply(boolean.test,1,sum)
#alleles.called[1:5]

cols<-c("REFs","ALTs", "chr", "position","snp","ref_base","A_forward","A_reverse","A_total","C_forward","C_reverse","C_total","G_forward","G_reverse","G_total","T_forward","T_reverse","T_total","del_forward","del_reverse","del_total","total_forward","total_reverse","totalCount")

dim(data)
#pval.extra[1:3,]
#data[1:3,]
#boolean.test[1:3,]
summary<-cbind(pval.extra,boolean.test,alleles.called,data[,cols])

#summary<-cbind(pval.extra,boolean.test,alleles.called,data[,cols])
length(colnames(summary))
length(unique(colnames(summary)))

#summary[1:5,]
table(summary[,"alleles.called"])
test<-summary[,"alleles.called"]==0
#summary[test,][1:10,]
#data[test,]

##################need test here to select heteroduplex calls
reject<-summary[,"alleles.called"]==0 | summary[,"alleles.called"]>2 | summary[,"total.bases"] < sub.read.threshold
print(paste("bases called"))
print(table(summary[!reject,"alleles.called"]))
print("-----------------------------------")
passing<-summary[!reject,]
test<-rownames(passing)[passing[,"alleles.called"]==2]
#passing[test,][1:10,]

##################################### choose which to call ##############  sub.read.threshold<-30


# "m150921_155657_42229_c100905102550000001823198804291661_s1_p0/69500"
# "LPH-001-1_lbc96_BCL2L12_Ex1_F_chr19_50169121_50169141.txt"


chk.duplex<-subset(summary,subset=(rownames(summary) %in% test ),select=c(cols,"error.rate"))

number<-dim(chk.duplex)[1]
number<-min(number,5)
print(file)
print(paste0("Duplexs found :",dim(chk.duplex)[1]))
#print(chk.duplex[1:number,])
columns<-c("single.error","duplex.error","balance.error","strand.error","min.alt.reads.error","other.ref.too.large","error.rate.duplex","allowed.errors", "alts.total","alts.reads","alts.dels","refs.total","refs.reads","refs.dels", "refs.total.other","refs.reads.other","refs.dels.other")
      
dirns<-c("forward","reverse")

if(dim(chk.duplex)[1]>0){

unique.cells<-rownames(chk.duplex)
duplex.checks<-matrix(data=NA,nrow=length(unique.cells),ncol=length(columns))
colnames(duplex.checks)<-columns
rownames(duplex.checks)<-unique.cells

dim(duplex.checks)

# grep("m151216_051615_42229_c100935702550000001823204005251625_s1_p0/125158",rownames(chk.duplex))
#grep("m151216_051615_42229_c100935702550000001823204005251625_s1_p0/125158",rownames(chk.duplex))

#i<-25
#chk.duplex[i,]
for( i in 1:dim(chk.duplex)[1]){
    alt.locations<-paste(chk.duplex[i,"ALTs"],dirns,sep="_")
    the.order<-order(chk.duplex[i,alt.locations],decreasing=TRUE)
    alt.loc<-dirns[the.order[1]]

    alt.locations<-paste(chk.duplex[i,"REFs"],dirns,sep="_")
    the.order<-order(chk.duplex[i,alt.locations],decreasing=TRUE)
    ref.loc<-dirns[the.order[1]]

    if( alt.loc==ref.loc){strand.error<-TRUE}else{strand.error<-FALSE}
    
    alts.reads<-chk.duplex[i,paste(chk.duplex[i,"ALTs"],alt.loc,sep="_")]
    alts.dels<-chk.duplex[i,paste("del",alt.loc,sep="_")]
    alts.total<-alts.reads+alts.dels

    refs.reads<-chk.duplex[i,paste(chk.duplex[i,"REFs"],ref.loc,sep="_")]
    refs.dels<-chk.duplex[i,paste("del",ref.loc,sep="_")]
    refs.total<-refs.reads+refs.dels

    refs.reads.other<-chk.duplex[i,paste(chk.duplex[i,"REFs"],dirns[dirns!=ref.loc],sep="_")]
    refs.dels.other<-chk.duplex[i,paste("del",dirns[dirns!=ref.loc],sep="_")]
    refs.total.other<-refs.reads.other+refs.dels.other
    

    cbind(alts.total,alts.reads,alts.dels,refs.total,refs.reads,refs.dels)
    
    error.rate.duplex<-chk.duplex[i,"error.rate"]
   # error.rate.duplex<-abs(chk.duplex[i,"totalCount"]-(alts.total+refs.total))/chk.duplex[i,"totalCount"] ### this is the error rate for basex excluding the main heteroduplex alleles
   # if(error.rate.duplex<chk.duplex[i,"error.rate"]){error.rate.duplex<-chk.duplex[i,"error.rate"]} # error rate lass than previosly calculted rate
    
    mean.error<- refs.total*error.rate.duplex
    sd.error<-sqrt(refs.total*error.rate.duplex*(1-error.rate.duplex))
    allowed.errors<-ceiling(mean.error+3*sd.error)+3 ## plus 2 incase loop startes on wrong base

    single.error<-FALSE
    balance.error<-FALSE
    duplex.error<-FALSE
    min.alt.reads.error<- FALSE
   other.ref.too.large<- FALSE
 
   duplex.checks[i,]<-c(single.error,duplex.error,balance.error,strand.error,min.alt.reads.error,other.ref.too.large,error.rate.duplex,allowed.errors,alts.total,alts.reads,alts.dels,refs.total,refs.reads,refs.dels, refs.total.other,refs.reads.other,refs.dels.other)
#        columns<-c("duplex.error","balance.error","strand.error","error.rate.duplex","min.alt.reads.error","allowed.errors","alts.total","alts.reads","alts.dels","refs.total","refs.reads","refs.dels")                                        
}
######################redo using the median error rate - same can be inflated


hist(duplex.checks[,"error.rate.duplex"])

median(duplex.checks[,"error.rate.duplex"])

# duplex.checks[1:5,] mean.error.rate.duplex<-chk.duplex[i,"error.rate"]
mean.error.rate.duplex<-max(mean(duplex.checks[,"error.rate.duplex"],trim=0.1),median(duplex.checks[,"error.rate.duplex"]))
print(paste0("duplex error rate ",mean.error.rate.duplex))
mean.error.rate.duplex<-min(mean.error.rate.duplex,0.05)
print(paste0("Chosen duplex error rate ",mean.error.rate.duplex))
#sd.error.rate.duplex<-sd(duplex.checks[,"error.rate.duplex"])
i<-4
#duplex.checks[i,]
for( i in 1:dim(duplex.checks)[1]){
    refs.total<-duplex.checks[i,"refs.total"]
    alts.total<-duplex.checks[i,"alts.total"]
    alts.reads<-duplex.checks[i,"alts.reads"]
    
    mean.error<- refs.total*mean.error.rate.duplex
    sd.error<-sqrt(refs.total*mean.error.rate.duplex*(1-mean.error.rate.duplex))
    allowed.errors<-as.integer(mean.error+3*sd.error)+3 ## plus 2 incase loop startes on wrong base

    refs.total.other<-duplex.checks[i,"refs.total.other"]
    refs.reads.other<-duplex.checks[i,"refs.reads.other"]
    

    diff.hetero<-abs(refs.total-alts.total)
    balance.error<-diff.hetero > allowed.errors
    min.alt.reads.error<-alts.reads<=min.alt.reads.dups
    other.ref.too.large<- refs.reads.other > allowed.errors
    
    
    duplex.error<-balance.error | strand.error | min.alt.reads.error |  other.ref.too.large

    duplex.checks[i,"duplex.error"]<-duplex.error
    duplex.checks[i,"balance.error"]<-balance.error
    duplex.checks[i,"allowed.errors"]<-allowed.errors
    duplex.checks[i,"min.alt.reads.error"]<-min.alt.reads.error
    duplex.checks[i,"other.ref.too.large"]<- other.ref.too.large
                                                
}

  sum( duplex.checks[,"duplex.error"]) 
#duplex.checks.ori<-duplex.checks
#summary[1:5,]
#duplex.checks[1:5,]

posns<-match(rownames(summary),rownames(duplex.checks))
if(length(posns)>2){  ### just in case only one entry
duplex.checks<-duplex.checks[posns,]
rownames(duplex.checks)<-rownames(summary)
}

}else{ # no duplexs make a dummy  duplex.checks same size as summary
unique.cells<-rownames(summary)
duplex.checks<-matrix(data=NA,nrow=length(unique.cells),ncol=length(columns))
colnames(duplex.checks)<-columns
rownames(duplex.checks)<-unique.cells
}  # duplex >0



#grep("m151215_115908_42229_c100935702550000001823204005251621_s1_p0/131474",rownames(duplex.checks))
#duplex.checks[754,]
#duplex.checks[25,]

#summary[,"single.error"]
duplex.checks[is.na(duplex.checks[,"duplex.error"]),"duplex.error"]<-0 ## Na are not errors
duplex.checks[is.na(duplex.checks[,"balance.error"]),"balance.error"]<-0
duplex.checks[is.na(duplex.checks[,"strand.error"]),"strand.error"]<-0
duplex.checks[is.na(duplex.checks[,"min.alt.reads.error"]),"min.alt.reads.error"]<-0
duplex.checks[is.na(duplex.checks[,"other.ref.too.large"]),"min.alt.reads.error"]<-0
# summary[1:5,]
#duplex.checks[1:5,]

colnames(duplex.checks) %in% colnames(summary) 
summary<-cbind(summary,duplex.checks)

colnames(summary)
#################### final check that balance in ok for single allele calls
#chk.singles<-summary[,summary[,"alleles.called"]==1
the.pval.cols<-paste(base.cols,"pval",sep=".")

chk.single<-subset(summary,subset=(summary[,"alleles.called"]==1),select=c(cols,"error.rate",the.pval.cols,"single.error"))
       
#chk.single[i,]
#i<-1
if(dim(chk.single)[1]>0 ){
for( i in 1:dim(chk.single)[1]){
    
    called.allele<-names(sort(chk.single[i,the.pval.cols],decreasing=FALSE))[1]
    called.allele<-gsub("_total.pval","",called.allele)
    
    alt.locations<-paste(called.allele,dirns,sep="_")
    the.order<-order(chk.single[i,alt.locations],decreasing=TRUE)
    ref.loc<-dirns[the.order[1]]

    refs.reads<-chk.single[i,paste(called.allele,ref.loc,sep="_")]
    refs.dels<-chk.single[i,paste("del",ref.loc,sep="_")]
    refs.total<-refs.reads+refs.dels

    refs.reads.other<-chk.single[i,paste(called.allele,dirns[dirns!=ref.loc],sep="_")]
    refs.dels.other<-chk.single[i,paste("del",dirns[dirns!=ref.loc],sep="_")]
    refs.total.other<-refs.reads.other+refs.dels.other
    

    cbind(alts.total,alts.reads,alts.dels,refs.total,refs.reads,refs.dels)
    
    error.rate<-chk.single[i,"error.rate"]
   # error.rate.duplex<-abs(chk.single[i,"totalCount"]-(alts.total+refs.total))/chk.single[i,"totalCount"] ### this is the error rate for basex excluding the main heteroduplex alleles
   # if(error.rate.duplex<chk.single[i,"error.rate"]){error.rate.duplex<-chk.single[i,"error.rate"]} # error rate lass than previosly calculted rate
    
    mean.error<- refs.total*error.rate.duplex
    sd.error<-sqrt(refs.total*error.rate.duplex*(1-error.rate.duplex))
    allowed.errors<-ceiling(mean.error+3*sd.error)+3 ## plus 2 incase loop startes on wrong base
    single.error<-abs(refs.total-refs.total.other)> allowed.errors
    chk.single[,"single.error"]<-single.error
    
#        columns<-c("duplex.error","balance.error","strand.error","error.rate.duplex","min.alt.reads.error","allowed.errors","alts.total","alts.reads","alts.dels","refs.total","refs.reads","refs.dels")                                        
}


posns<-match(rownames(summary),rownames(chk.single))
summary[,"single.error"]<-chk.single[posns,"single.error"]
summary[is.na(summary[,"single.error"]),"single.error"]<-FALSE
    
} ## no singles to analyse
    






##################################################################################################################################
##################################################################################################################################
################################################### MAKe calles ########################################################
##################################################################################################################################
##################################################################################################################################

reject<-summary[,"alleles.called"]==0 | summary[,"alleles.called"]>2 | (summary[,"total.bases"] < sub.read.threshold) | (summary[,"duplex.error"]==1) | summary[,"single.error"]
called.bases<-paste(bases,"total.count",sep="_")

reject.no.call.boolean<-(summary[,"alleles.called"]==0)
reject.multi.allele.boolean<-(summary[,"alleles.called"]>2)
reject.coverage.boolean<-((summary[,"total.bases"] < sub.read.threshold))
reject.duplex.error.boolean<-((summary[,"duplex.error"]==1))
reject.duplex.min.alt.reads.boolean<-((summary[,"min.alt.reads.error"]==1))
reject.other.ref.too.large.boolean<-sum((summary[,"other.ref.too.large"]==1))
reject.single.error.boolean<-sum(summary[,"single.error"])

sum(!reject & summary[,"single.error"])
reject.no.call<-sum(summary[,"alleles.called"]==0)
reject.multi.allele<-sum(summary[,"alleles.called"]>2)
reject.coverage<-sum((summary[,"total.bases"] < sub.read.threshold))
reject.duplex.error<-sum((summary[,"duplex.error"]==1))
reject.duplex.min.alt.reads<-sum((summary[,"min.alt.reads.error"]==1))
reject.other.ref.too.large<-sum((summary[,"other.ref.too.large"]==1))
reject.single.error<-sum(summary[,"single.error"])

print(basename)
print("Alleles Called")
print(table(summary[!reject,"alleles.called"]))


Well.ID<-rownames(summary)
write.table(cbind(Well.ID,reject,reject.no.call.boolean,reject.multi.allele.boolean,reject.coverage.boolean,reject.duplex.error.boolean,reject.duplex.min.alt.reads.boolean,reject.other.ref.too.large.boolean,reject.single.error.boolean,mean.error.rate.duplex,summary),file=paste(dir,"/significance/",basename,".summary.txt",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
getwd()

called.bases<-paste(bases,"total.count",sep="_")


SNP<-unique(summary[,"snp"])
REF<-unique(summary[,"REFs"])
ALT<-unique(summary[,"ALTs"])

ALT.col<-paste(ALT,"total.count",sep="_")
REF.col<-paste(REF,"total.count",sep="_")

chr<-unique(summary[,"chr"])
position<-unique(summary[,"position"])
rejected<-sum(reject)
non.rejected<-sum(!reject)

if(sum(!reject)>0){
the.counts<-apply(summary[!reject,called.bases],2,sum)
the.counts<-the.counts[called.bases]
num.duplex<-sum(summary[!reject,"alleles.called"]==2)
ALT.counts<-sum(summary[!reject,ALT.col])
REF.counts<-sum(summary[!reject,REF.col])
MAF<-ALT.counts/(ALT.counts+REF.counts)

}else{
the.counts<-rep(0,times=length(called.bases))
names(the.counts)<-called.bases
num.duplex<-0
ALT.counts<-0
REF.counts<-0
MAF<-0
}

a.call<-c(basename,SNP,chr,position,REF,ALT,REF.counts,ALT.counts,MAF,the.counts,num.duplex,rejected,non.rejected,reject.no.call,reject.multi.allele,reject.coverage,reject.duplex.error,reject.duplex.min.alt.reads,reject.other.ref.too.large,reject.single.error,mean.error.rate.duplex)



if(ifiles==1){
    all.calls<-a.call
}else{
    all.calls<-rbind(all.calls,a.call)
  }

#print(all.calls)
print("-----------------------------------")
print("-----------------------------------")
print("-----------------------------------")

} # loop over files

#150/12000= 0.013
grep("LPH-001-21_lbc94_TP53_Ex6_F_chr17_7577528_7577548",args)
summary[test,]


#all.calls.ori<-all.calls

dim(all.calls)
all.calls[1:5,]
#colnames(all.calls)<-c("basename","SNP","chr","position","REF","ALT","REF.counts","ALT.counts","MAF",called.bases,"num.duplex","rejected","non.rejected","reject.no.call","reject.multi.allele","reject.coverage","reject.duplex.error","reject.duplex.min.alt.reads","mean.error.rate.duplex")

colnames(all.calls)<-c("basename","SNP","chr","position","REF","ALT","REF.counts","ALT.counts","MAF",called.bases,"num.duplex","rejected","non.rejected","reject.no.call","reject.multi.allele","reject.coverage","reject.duplex.error","reject.duplex.min.alt.reads","reject.other.ref.too.large","reject.single.error","mean.error.rate.duplex")
getwd()
write.table(all.calls,file="Summary_BEST_trimmed.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
















## ################################### DEVELOPEMENT WORK
## ##Bias coin example
## refs.total<-6
## error.rate.duplex<-0.3

## dbinom(as.integer(0),   refs.total, error.rate.duplex, log= FALSE)  # probability of 0 heads
## dbinom(as.integer(6),   refs.total, error.rate.duplex, log= FALSE)  # probability of 6 heads

## mean.error<- refs.total*error.rate.duplex
## sd.error<-sqrt(refs.total*error.rate.duplex*(1-error.rate.duplex))
    
##         dbinom(as.integer(mean.error+sd.error),   refs.total, error.rate.duplex, log= FALSE)
##         dbinom(as.integer(mean.error+2*sd.error),   refs.total, error.rate.duplex, log= FALSE)
##         dbinom(as.integer(mean.error+3*sd.error),   refs.total, error.rate.duplex, log= FALSE)
##         dbinom(as.integer(mean.error+4*sd.error),   refs.total, error.rate.duplex, log= FALSE)
##         dbinom(as.integer(mean.error+5*sd.error),   refs.total, error.rate.duplex, log= FALSE)
##         dbinom(as.integer(mean.error+6*sd.error),   refs.total, error.rate.duplex, log= FALSE)


    
## p=0.5
## n=100
## x<-seq(from=0,to=6,by=1)  ## z scores
## y<-x*sqrt(n*p*(1-p)) + n*p  ## actual values
## x
## y
## pbinom(y, n, p, lower.tail = FALSE, log.p = TRUE)
## pnorm(x, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
## pnorm(x, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)


## dnorm(as.integer(y), mean = n*p , sd = sqrt(n*p*(1-p)), log = FALSE) ## normal distribution value
## dbinom(as.integer(y), n, p, log = FALSE) ## normal distribution value

## ########### function values change with variance
## dnorm(0, mean = 0 , sd = sqrt(0.2), log = FALSE)
## dnorm(0, mean = 0 , sd = sqrt(1), log = FALSE)
## dnorm(0, mean = 0 , sd = sqrt(5), log = FALSE)

## ########################### use cumulative probability
## pnorm(0, mean = 0 , sd = sqrt(1), lower.tail = TRUE, log = FALSE)
## pnorm(0, mean = 0 , sd = sqrt(5), lower.tail = TRUE, log = FALSE)

## pnorm(2, mean = 0 , sd = 1, lower.tail = TRUE, log = FALSE)


## (2*pnorm(x, mean = 0 , sd = 1, lower.tail = FALSE, log = FALSE) )-1 #### this is the error function
## -(2*pnorm(x, mean = 0 , sd = 1, lower.tail = FALSE, log = FALSE) )+1 #### this is the error function

## qnorm(0.5, mean = 0 , sd = 1, lower.tail = TRUE, log = FALSE)
## qnorm(0.9, mean = 0 , sd = 1, lower.tail = FALSE, log = FALSE)
## pbinom(q, size, prob)

## p=0.00003
## n=40
## pbinom(22, n, p, lower.tail = TRUE, log.p = FALSE) + pbinom(22, n, p, lower.tail = FALSE, log.p = FALSE)
## ########## is one 


## dbinom- gives the value

## ##############3
## -(2*pnorm(x, mean = 0 , sd = 3, lower.tail = FALSE, log = FALSE) )+1 #### this is the error function
## pnorm(x, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)

## pnorm(2, mean = 0 , sd = sqrt(5), lower.tail = TRUE, log = FALSE)


## dnorm(0, mean = 0 , sd = sqrt(0.5), log = FALSE)

## pnorm(-1, mean = -2 , sd = sqrt(0.5), lower.tail = TRUE, log = FALSE)


## dnorm(0, mean = 0 , sd = sqrt(0.2), log = FALSE)


## p=0.001
## n=50
## n*p
## measure<-2
## sd<-sqrt(n*p*(1-p))

## x*sqrt(n*p*(1-p)) + n*p
## z<-(measure-n*p)/sd
## z
## -(2*pnorm(z, mean = 0 , sd = 1, lower.tail = FALSE, log = FALSE) )+1 #### this is the error function
## (2*pnorm(z, mean = 0 , sd = 1, lower.tail = FALSE, log = FALSE) ) #### this is the complementary error function
## pnorm(z, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)


## #pnbinom(measure, n, p, lower.tail = FALSE, log.p = FALSE)
## pbinom(measure, n, p, lower.tail = FALSE, log.p = FALSE)

## pbinom(measure, n, p, lower.tail = TRUE, log.p = FALSE) + pbinom(measure, n, p, lower.tail = FALSE, log.p = FALSE)


## pbinom(1000, n, p, lower.tail = FALSE, log.p = FALSE) ### need more reads than expected

## k<-seq(from=0,to=100,by=1)
## dbinom(as.integer(mean.error+6*sd.error),   n, p, log= FALSE)
##  plot (k, dbinom(k, n, p, log = FALSE), type = "l")
## # plot (k, dnbinom(k, n, p, log = FALSE), type = "l")
## dbinom(k, n, p, log = TRUE)
##  plot (k, pbinom(k, n, p, lower.tail = TRUE, log.p = FALSE), type = "l")
## ############################
## ## pbinom(1,6,0.5, lower.tail = TRUE, log.p = FALSE)
## ## qbinom(0.0005,2,0.5, lower.tail = FALSE,log = FALSE)
## ## dbinom(2,2,0.3,log=FALSE)


## ## p=0.015
## ## n=50
## ## alt.counts.thresh<-2

## ## z<-(alt.counts.thresh- n*p) / sqrt(n*p*(1-p))
## ## z

## ## 2*(pnorm(abs(c(z)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))

## ## qbinom(c(0,0.5,0.05,0.005,0.0005), n, p, lower.tail = FALSE, log.p = FALSE)

## ## pbinom(abs(c(alt.counts.thresh)), n, p, lower.tail = FALSE, log.p = FALSE)

## ## rbinom(abs(10), n, p)

## ## require(graphics)
## ## # Compute P(45 < X < 55) for X Binomial(100,0.5)
## ## sum(dbinom(46:54, 100, 0.5))

## ## ## Using "log = TRUE" for an extended range :
## ## n <- 2000
## ## k <- seq(0, n, by = 20)
## ## plot (k, dbinom(k, n, pi/10, log = TRUE), type = "l", ylab = "log density",
## ##       main = "dbinom(*, log=TRUE) is better than  log(dbinom(*))")
## ## lines(k, log(dbinom(k, n, pi/10)), col = "red", lwd = 2)
## ## ## extreme points are omitted since dbinom gives 0.
## ## mtext("dbinom(k, log=TRUE)", adj = 0)
## ## mtext("extended range", adj = 0, line = -1, font = 4)
## ## mtext("log(dbinom(k))", col = "red", adj = 1)




## ## pbinom(q, size, prob, lower.tail = TRUE, log.p = FALSE)
         
## ## ## p.threshold=0.0026 # z=3:   2*(pnorm(abs(c(1:8)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
## ## ## p.threshold=2*(pnorm(abs(c(2)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
## ## ## p.threshold
