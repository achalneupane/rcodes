
#############################################
## 1) Read in data at SNP location
## 2) estimate backgroud substitution error rare for "non SNP" alleles for a given SNP
## 3) For each well apply binomial approximation and obtain a p-value for null hypothesis  for SNP bases (2 per well) 


##################################################################################################################
## 1) Read in data at SNP location
pile.up.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Pac_bio_replication/LPH-001-1_AK1_lbc62_STK19_SNP.csv"
data<-read.delim(pile.up.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

setwd(dirname(pile.up.file))

dim(data)
data[1:5,]

##### get alleles for this location from "snp" column in format (chr:start:end:REF:ALT:type)
ref.posn<-4
alt.posn<-5

alleles<-strsplit(data[,"snp"],split=":")
REFs<-unlist(lapply(alleles,function(x) {x[4] } ))
ALTs<-unlist(lapply(alleles,function(x) {x[5] } ))

#### tridy UP AND CHECKS
### check have a ref and an alt always
bad.locations<-is.na(REFs) | is.na(ALTs)
print("Number of bad locations:", sum(bad.locations))
print(data[bad.locations,])
data<-cbind(REFs,ALTs,data)

data<-data[!bad.locations,]

data[1:5,]

##################################################################################################################
data.sub[1:5,]

###################### Get
#############################################
## 2) estimate backgroud substitution error rare form "non SNP" alleles for a given SNP

########### LOOP over SNPs (this could be over all samples, but can just use an individals)


################ get extimate of the backgroud error rate for a snps 
bases<-c("A","C","G","T")
base.cols<-paste(bases,"total",sep="_")
the.snps<-unique(data[,"snp"])

pval.columns<-c("error.rate","total.bases",paste(base.cols,"pval",sep="."))
pval.extra<-matrix(data=NA,nrow=dim(data)[1],ncol=length(pval.columns))
colnames(pval.extra)<-pval.columns

#### for( isample= 1:length(samples)) {  ### is have data for all samples i one file
i<-1
for (i in 1:length(the.snps)){
    data.sub <-data[data[,"snp"]==the.snps[i],]
    allele.counts<-apply(data.sub[,base.cols],2,function(x) sum(x,na.rm=TRUE))
    total.counts<-sum(allele.counts)
   ##  A_total C_total G_total T_total 
   ## 3123      75    6060      74
   data.sub[1:5,]
    ### get backgroup subsitution error rate by excluing bases when snps is expected.
    snp.bases<-unique(c(as.character(data.sub[,"REFs"]),as.character(data.sub[,"ALTs"])))
    background.bases<-base.cols[!(base.cols %in% paste(snp.bases,"total",sep="_"))]
    background.error<-sum(allele.counts[background.bases])/total.counts
    pval.extra[data[,"snp"]==the.snps[i],"error.rate"]<-background.error ### pval.data same size as data keep error rate for that snp
}

##################################################################################################################
##################################################################################################################
## 3) For each well apply binomial approximation and obtain a p-value for null hypothesis  each well ... do for all bases though only need ref and alt
##

data[1:5,]
pval.extra[1:5,]

base.counts<-apply(data[,base.cols],1,function(x) sum(x,na.rm=TRUE))
base.counts[1:5]
data[1:5,base.cols]
pval.extra[,"total.bases"]<-base.counts ## total bases counted in a well

pval.extra[1:5,]
## > pval.extra[1:5,]
##      error.rate p.ref p.alt
##      error.rate total.bases A_total.pval C_total.pval G_total.pval T_total.pval
## [1,] 0.01596657          13           NA           NA           NA           NA
## [2,] 0.01596657         119           NA           NA           NA           NA
## [3,] 0.01596657         121           NA           NA           NA           NA
## [4,] 0.01596657          58           NA           NA           NA           NA
## [5,] 0.01596657          10           NA           NA           NA           NA

########### get refernce p_values looping over all possible reference bases
i<-3
for (i in 1:length(base.cols)){

## p.base<-pbinom(data[,base.cols[i]], pval.extra[,"total.bases"], pval.extra[,"error.rate"], lower.tail = FALSE, log.p = FALSE) ### this is wrong correted Nov10 2015 use below

p.base<-dbinom(data[,base.cols[i]], pval.extra[,"total.bases"], pval.extra[,"error.rate"], lower.tail = FALSE, log.p = FALSE)


 pval.extra[,paste(base.cols[i],"pval",sep=".")]<-signif(p.base,digits=7)
}

pval.extra[1:5,]

data[1:5,]
 


write.table(cbind(data,pval.extra),file="Pac.bio.call.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
getwd()

150/12000= 0.013


## n=6
## p=0.3
## alt.counts.thresh=3
## (alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) 
## z<-(alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) 
## (pnorm((alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) , mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
## dbinom(alt.counts.thresh,n,p,log=FALSE)


## p=0.015
## n=50
## alt.counts.thresh<-2

## z<-(alt.counts.thresh- n*p) / sqrt(n*p*(1-p))
## z

## 2*(pnorm(abs(c(z)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))

## qbinom(c(0,0.5,0.05,0.005,0.0005), n, p, lower.tail = FALSE, log.p = FALSE)

## pbinom(abs(c(alt.counts.thresh)), n, p, lower.tail = FALSE, log.p = FALSE) 

## rbinom(abs(10), n, p)

## require(graphics)
## # Compute P(45 < X < 55) for X Binomial(100,0.5)
## sum(dbinom(46:54, 100, 0.5))

## ## Using "log = TRUE" for an extended range :
## n <- 2000
## k <- seq(0, n, by = 20)
## plot (k, dbinom(k, n, pi/10, log = TRUE), type = "l", ylab = "log density",
##       main = "dbinom(*, log=TRUE) is better than  log(dbinom(*))")
## lines(k, log(dbinom(k, n, pi/10)), col = "red", lwd = 2)
## ## extreme points are omitted since dbinom gives 0.
## mtext("dbinom(k, log=TRUE)", adj = 0)
## mtext("extended range", adj = 0, line = -1, font = 4)
## mtext("log(dbinom(k))", col = "red", adj = 1)




## pbinom(q, size, prob, lower.tail = TRUE, log.p = FALSE)
         
## ## p.threshold=0.0026 # z=3:   2*(pnorm(abs(c(1:8)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
## ## p.threshold=2*(pnorm(abs(c(2)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
## ## p.threshold
