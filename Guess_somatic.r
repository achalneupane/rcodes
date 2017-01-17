




cellularity<- pheno.ori[,c("SAMPLE","Blast.Count..")]
cellularity[,"Blast.Count.."]<-cellularity[,"Blast.Count.."]/100
cellularity[is.na(cellularity[,"Blast.Count.."]) & pheno.ori[,"Control"] ,"Blast.Count.."]<-1 ## unnoen then make pure
# cellularity[pheno.ori[,"AML"],]
# cellularity[pheno.ori[,"Control"],]
# sample.types

PDs<-pheno.ori[pheno.ori[,"AML-Child"] | pheno.ori[,"Asian-AML-Child"] | pheno.ori[,"Asian-AML"]  | pheno.ori[,"AML-NotDiagnosis-Child"] | pheno.ori[, "Asian-AML-NotDiagnosis-Child"] | pheno.ori[,"Asian-Control"],"SAMPLE"]
## PDs<-c(PDs,"LPH-001-27_PD")   table(pheno.ori[pheno.ori[,"AffectionStatus"]==2,"SAMPLE"] )
PDs
#PDs.alt.counts<-alt.reads.reference.calls(a.indel,PDs,threshold=1)


cancer<-pheno.ori[pheno.ori[,"AML"],"SAMPLE"]
cancer
Controls<-pheno.ori[pheno.ori[,"Control"],"SAMPLE"]

a.indel.stats[1:5,1:20]


Control.alt.counts<-alt.reads.reference.calls(a.indel.stats,Controls,AD.extension="FAD",threshold=1) ## needs a sample.GT column

## snp.only<-grepl("^snp",a.indel.ori[,"TYPE"]) ### the p.values codes uses "AD" so does not work for SNPs
## is.flat<-grepl("flat$",a.indel.ori[,"TYPE"])
## has.one.geno<-as.numeric(summary.geno.extra.ori[,"ALT.Alleles.AML"])>0
## alt.count.thresh<-1  # use  p[true.sample.minor.counts <= alt.count.thresh]<-1 default more than one needed
## has.one.geno.ori<-has.one.geno
## to.recover<-snp.only & has.one.geno & !is.flat  & maf.filter & rare.in.Control  ## maf.filter & rare.in.Control probably not required
recover.samples<-c(cancer,PDs)
## sum(to.recover)
sum(pass)
to.recover<-pass


alt.count.thresh<-1  # use  p[true.sample.minor.counts <= alt.count.thresh]<-1 default more than one needed
alt.count.thresh.include<-2 # ignore above mor that this number of samples have true calls
AD.lower.tail <-FALSE


# geno.p<-genotype.p.values(a.indel.stats[to.recover ,] ,c(cancer,PDs),AD.extension="FAD",Control.alt.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh) ## COULD USE NORMAL HERE
# geno.p<-genotype.p.values(a.indel.stats[to.recover ,],recover.samples,AD.extension="FAD",Control.alt.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh) ## COULD USE NORMAL

geno.p<-genotype.p.values.row(a.indel.stats[to.recover ,],recover.samples,AD.extension="FAD",Control.alt.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh,alt.count.thresh.include,AD.lower.tail) ## COULD USE 
                                        #geno.p<-genotype.p.values(a.indel.stats[to.recover ,],recover.samples,AD.extension="FAD",Control.alt.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh,AD.lower.tail) ## COULD USE NORMAL HERE


p.threshold.z.thresh<-6
p.threshold=2*(pnorm(abs(c(p.threshold.z.thresh)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
p.threshold
found.genotype<-  geno.p <= p.threshold
geno.p[found.genotype]<-"0/1" # p.threshold<- 0.0026
geno.p[!found.genotype]<-"0/0"
colnames(geno.p)<-paste(colnames(geno.p),".GT",sep="")


############# I don't want to correct GATK calls that are already 0/1 to 0/0:


genotypes<-a.indel.stats[to.recover ,paste(recover.samples,"GT",sep=".")] 
ref.call<-genotypes =="0/0"
geno.p[!ref.call]<-genotypes[!ref.call]


found.genotype[!ref.call]<-FALSE ## no recovery of 0/1 or 1/1
dim(found.genotype)










genotypes<-a.indel.stats[to.recover ,paste(recover.samples,"GT",sep=".")] 
ref.call<-genotypes =="0/0"
geno.p[!ref.call]<-genotypes[!ref.call]


Control.alt.Non.REF.counts<-alt.reads.Non.reference.calls(a.indel.stats,Controls,AD.extension="FAD",threshold=1) ## needs a sample.GT column

has.an.alt <- as.numeric(summary.geno.extra.ori[to.recover , "ALT.Alleles.Control" ])>0

snp.only<-grepl("^snp",a.indel.ori[,"TYPE"]) ### the p.values codes uses "AD" so does not work for SNPs
is.flat<-grepl("flat$",a.indel.ori[,"TYPE"])
has.one.geno<-as.numeric(summary.geno.extra.ori[,"ALT.Alleles.AML"])>0
alt.count.thresh<-1  # use  p[true.sample.minor.counts <= alt.count.thresh]<-1 default more than one needed
has.one.geno.ori<-has.one.geno

to.recover<-snp.only & has.one.geno & !is.flat   ## maf.filter & rare.in.Control probably not required
recover.samples<-c(cancer,PDs)
sum(to.recover)
AD.lower.tail <-"BOTH"


alt.count.thresh<-1  # use  p[true.sample.minor.counts <= alt.count.thresh]<-1 default more than one needed
alt.count.thresh.include<-2 # ignore above more that this number of samples have true calls
AD.lower.tail <-"BOTH"

#geno.p.Non.REF<-genotype.p.values(a.indel.stats[to.recover ,],recover.samples,AD.extension="FAD",Control.alt.Non.REF.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh,AD.lower.tail) ## COULD USE NORMAL
geno.p.Non.REF<-genotype.p.values.row(a.indel.stats[to.recover ,],recover.samples,AD.extension="FAD",Control.alt.Non.REF.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh,alt.count.thresh.include,AD.lower.tail) ##
#geno.p<-genotype.p.values(a.indel.stats[to.recover ,],recover.samples,AD.extension="FAD",Control.alt.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh,AD.lower.tail) ## COULD USE NORMAL HERE

geno.p.Non.REF[1:5,1:5]

p.threshold.z.thresh<-6
p.threshold=2*(pnorm((c(p.threshold.z.thresh)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
p.threshold
found.genotype.Non.REF<-  geno.p.Non.REF <= p.threshold

## sort(apply(found.genotype.Non.REF[has.an.alt,],1,sum),decreasing=TRUE)[1:20]
## sort(apply(found.genotype.Non.REF[!has.an.alt,],1,sum),decreasing=TRUE)[1:20]

genotypes<-a.indel.stats[to.recover ,paste(recover.samples,"GT",sep=".")] 
ref.call<-genotypes =="0/0"
#geno.p[!ref.call]<-genotypes[!ref.call]


geno.p.Non.REF[!has.an.alt,]<-1
found.genotype.Non.REF[!has.an.alt,]<-FALSE  #### don't do ones that were no 0/1 calles in controls

############# I don't want to correct GATK calls that are already 0/1 to 0/0:
found.genotype.Non.REF[ref.call]<-FALSE ## no recovery of 0/0
geno.p.Non.REF[ref.call]<-1


##  found.genotype.Non.REF[1:5,1:5]
##              found.genotype.Non.REF[test,paste0(c(samples),".GT")]
##          a.indel.stats[test,paste0(c(samples),".FAD")]

## summary.geno.extra[test,c("GENO.AML","GENO.Control")]
#################################################################
#somatic.matrix<-found.genotype | found.genotype.Non.REF
somatic.matrix<-found.genotype.Non.REF


somatic.matrix.desc<-matrix(data="GERM",nrow=dim(somatic.matrix)[1],ncol=dim(somatic.matrix)[2])
somatic.matrix.desc[somatic.matrix]<-"SOMA"

sort(apply(somatic.matrix.desc,1,function(x){sum(x=="SOMA")}),decreasing=TRUE)[1:20]

colnames(somatic.matrix.desc)<-paste(recover.samples,".GT",sep="")
rownames(somatic.matrix.desc)<-rownames(geno.p.Non.REF)

somatic.matrix.desc.full<-matrix(data="GERM",nrow=dim(a.indel)[1],ncol=length(all.possible.samples))
colnames(somatic.matrix.desc.full)<-paste(all.possible.samples,".GT",sep="")
rownames(somatic.matrix.desc.full)<-rownames(a.indel)

sum(!( colnames(somatic.matrix.desc) %in% colnames(a.indel))) # must ge zero
to.transfer<-colnames(somatic.matrix.desc)[colnames(somatic.matrix.desc) %in% colnames(somatic.matrix.desc.full)]

posns<-match(rownames(somatic.matrix.desc.full),rownames(somatic.matrix.desc))
missing<-is.na(posns)
sum(missing)
sum(!missing)
dim(somatic.matrix.desc.full)

somatic.matrix.desc.full[!missing,to.transfer]<-somatic.matrix.desc[posns[!missing],to.transfer]


#################################################################
colnames(geno.p.Non.REF)<-paste(colnames(geno.p.Non.REF),"GT",sep=".")
colnames(found.genotype.Non.REF)<-paste(colnames(found.genotype.Non.REF),"GT",sep=".")
found.genotype.Non.REF[1:5,1:5]

somatic.matrix.p.full<-matrix(data=1,nrow=dim(a.indel)[1],ncol=length(all.possible.samples))
colnames(somatic.matrix.p.full)<-paste(all.possible.samples,".GT",sep="")
rownames(somatic.matrix.p.full)<-rownames(a.indel)

to.transfer<-colnames(geno.p.Non.REF)[colnames(geno.p.Non.REF) %in% colnames(somatic.matrix.p.full)]

posns<-match(rownames(somatic.matrix.p.full),rownames(geno.p.Non.REF))
missing<-is.na(posns)
sum(missing)
sum(!missing)
dim(somatic.matrix.p.full)

somatic.matrix.p.full[!missing,to.transfer]<-geno.p.Non.REF[posns[!missing],to.transfer]

