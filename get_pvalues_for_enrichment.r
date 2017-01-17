
save(list=c("pheno.ori","a.indel.ori","Control.alt.counts","cellularity","summary.geno.extra.ori"),file="p-value.model.test.data.RData")






############# z scored based on a binomial distribution
#README
# http://en.wikipedia.org/wiki/Binomial_distribution
#  Mean in n*p
# variance np(1-p) sd =sqrt(var)
# Models use n as the coverage or reads at a given postion
## Z=X-np/sqrt(np(1-p)) :z-score of the minor allele frequency: (np is mean in Normal, X is measurement in tumor , the standard deviation is sqrt(np(1-p))
# and p as the frequency of the alternative allele in the control/Normal sample.

# the programs cacluales a z-score (the number of stand deviations the allele frequnecy in the tumor sample is away from the control group)
# the z-score is then converged to a p-value.
# a count for the control/tumor mox ins included which extimates the number of alternative alleles might be due to the control sample along

# the cellulaity is used to estimate the mix of tumor to normal tissue
# alt.reads.reference.calls is used to help estimate the alternative allele frequency in control or normal samples.

p<-2/66686
p<-0.001
n=32

(0-n*p)/sqrt(n*p*(1-p))

sd.thresh<-2

p

colnames(summary.geno.extra)
######################################
n<-max(as.integer(summary.geno.extra[,"TOTAL.Alleles.Control"]))
n
alt.counts.thresh<-1
while( (alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) <= sd.thresh){alt.counts.thresh<-alt.counts.thresh+1}
alt.counts.thresh

############################# rescue somatics with few calls with possion model
## cellularity<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Sequenza/NSLM.cellularity.summary.use.csv",header=T,sep="\t",fill=TRUE,skip=0,stringsAsFactors=FALSE)

## cellularity[1:5,]




code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")

working.directory<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/Figures" # directory where   "p-value.model.test.data.RData" located  
setwd(working.directory)

load("p-value.model.test.data.RData")
a.indel<-a.indel.ori
summary.geno.extra<-summary.geno.extra.ori

normals<-pheno.ori[pheno.ori[,"normal"],"SAMPLE"]
normals<-normals[normals!="LPH-001-27_PD"]
normals

Control<-pheno.ori[pheno.ori[,"Control"],"SAMPLE"]
Control
Control.alt.counts<-alt.reads.reference.calls(a.indel,Control,threshold=1)


PDs<-pheno.ori[pheno.ori[,"PD"],"SAMPLE"]
PDs<-c(PDs,"LPH-001-27_PD")
PDs
#PDs.alt.counts<-alt.reads.reference.calls(a.indel,PDs,threshold=1)


cancer<-pheno.ori[pheno.ori[,"cancer"],"SAMPLE"]
cancer

PDs<-pheno.ori[pheno.ori[,"PD"],"SAMPLE"]
PDs<-c(PDs,"LPH-001-27_PD")
PDs
#cancer.alt.counts<-alt.reads.reference.calls(a.indel,cancer,threshold=1)
#cancer.alt.counts.true<-alt.reads.Non.reference.calls(a.indel,cancer,threshold=1)


test<-c("chr13:21729956:21729956:A:G:snp","chr2:198363406:198363406:C:T:snp","chr22:20760282:20760282:A:C:snp")
test<-c("chr19:50169131:50169131:C:T:snp","chr15:40675107:40675107:C:T:snp","chr19:50169131:50169131:C:T:snp","chr17:7578263:7578263:G:A:snp","chr17:7574003:7574003:G:A:snp")
test<-c("chr12:11286221:11286221:T:A:snp","chr6:132029865:132029865:G:A:snp","chr7:82784807:82784807:C:A:snp","chr7:137206693:137206693:G:A:snp")

 test<-c("chr5:138268367:138268367:T:A:snp","chr12:11286221:11286221:T:A:snp")
a.indel[test,paste(normals,".AD",sep="")]
a.indel[test,paste(normals,".GT",sep="")]
a.indel[test,paste(PDs,".AD",sep="")]
a.indel[test,paste(cancer,".AD",sep="")]
a.indel[test,paste(cancer,".GT",sep="")]
##         normals
## normals.AD<-gsub(".GT",".AD",normals)

## a.indel["chr13:21729956:21729956:A:G:snp",normals]
## a.indel["chr13:21729956:21729956:A:G:snp",normals.AD]

geno.p<-genotype.p.values(a.indel[test,],cancer,Control.alt.counts[test,"Read.Balance"]/100,cellularity)

p.threshold=0.0026 # z=3:   2*(pnorm(abs(c(1:8)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))

found.genotype<-  geno.p<=p.threshold
geno.p[found.genotype]<-"0/1" # p.threshold<- 0.0026
geno.p[!found.genotype]<-"0/0"
colnames(geno.p)<-paste(colnames(geno.p),".GT",sep="")




snp.only<-grepl("^snp",a.indel[,"TYPE"]) ### the p.values codes uses "AD" so does not work for SNPs

has.one.geno<-as.numeric(summary.geno.extra[,"ALT.Alleles.cancer"])>0

geno.p<-genotype.p.values(a.indel[snp.only & has.one.geno ,] ,c(cancer,PDs),Control.alt.counts[snp.only & has.one.geno,"Read.Balance"]/100,cellularity) ## COULD USE NORMAL HERE


                           ## Control.alt.counts[1:5,]
                           ##  normal.alt.counts[1:5,]
#p.threshold=0.0026 # z=3:   2*(pnorm(abs(c(1:8)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
p.threshold=2*(pnorm(abs(c(6)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
p.threshold
found.genotype<-  geno.p <= p.threshold
geno.p[found.genotype]<-"0/1" # p.threshold<- 0.0026
geno.p[!found.genotype]<-"0/0"
colnames(geno.p)<-paste(colnames(geno.p),".GT",sep="")

## colnames(genotypes)
## colnames(geno.p)

## test<-c("chr13:21729956:21729956:A:G:snp","chr2:198363406:198363406:C:T:snp","chr22:20760282:20760282:A:C:snp")
## test<-c("chr19:50169131:50169131:C:T:snp") # ,"chr15:40675107:40675107:C:T:snp","chr19:50169131:50169131:C:T:snp","chr17:7578263:7578263:G:A:snp","chr17:7574003:7574003:G:A:snp")
## test<-c("chr12:11286221:11286221:T:A:snp","chr6:132029865:132029865:G:A:snp","chr7:82784807:82784807:C:A:snp","chr7:137206693:137206693:G:A:snp")

## geno.p[test,]
## genotypes[test,colnames(geno.p)]

## dim(genotypes)
## dim(geno.p)

#  a.indel.ori<-a.indel Only do this ONCE

sum(!( colnames(geno.p) %in% colnames(a.indel))) # must ge zero
to.transfer<-colnames(geno.p)[colnames(geno.p) %in% colnames(a.indel)]

posns<-match(rownames(a.indel),rownames(geno.p))
missing<-is.na(posns)
sum(!missing)

a.indel[!missing,to.transfer]<-geno.p[posns[!missing],to.transfer]


chk<-c("chr15:40675107:40675107:C:T:snp","chr19:50169131:50169131:C:T:snp") # BCL2L12
loc<-key %in% chk
targets<-the.projects  #c("NMD","ex.Control","AOGC")
targets
###### the.samples and pheno in same order but the.samples has .GT extension.
## it<-1
## for(it in 1:length(targets)){
## use.samples<-the.samples[pheno[,targets[it]]]
use.samples<-the.samples[pheno.ori[,"cancer"]]
print(targets[it])
print(use.samples)
length(use.samples)
genotypes<-a.indel.ori[loc,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),names(targets)[it],sep=".")
summary.geno
}

summary.geno.extra.ori[loc,c("MAF.cancer","GENO.cancer")]
summary.geno.extra[loc,c("MAF.cancer","GENO.cancer")]

#summary.geno[1:5,]
###################################################
