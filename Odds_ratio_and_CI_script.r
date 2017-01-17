#Need bim/bed/fam files from wence the results came from:

#Requirements
#1. Bim,Bed and Fam files
#2. FaST_LMM results file
#3. Freq file - plink --noweb --allow-no-sex --bfile  --freq --out
#4. Assoc file

#Odds Ratios

#getwd()
setwd("/media/scratch/uveitis/OR_CI_Paul_21q22_ERAP1_analysis_17Dec2012/ASwithUveitis")
options(width=300)
list=ls()
rm(list=ls())
options(stringsAsFactors=FALSE)

#freq file
f <- read.table(file="Freq.frq",header=TRUE,stringsAsFactors=F)
head(f)
nrow(f)

#fast lmm result file
p <- read.table(file="Result_fileWithUveitis.txt",header=F,stringsAsFactors=F)
colnames(p) <- c("SNP","Chromosome","GeneticDistance","Position","Pvalue","Qvalue","N","NullLogLike","AltLogLike","SNPWeight","SNPWeightSE","WaldStat","NullLogDelta","NullGeneticVar","NullResidualVar","NullBias")
head(p)
nrow(p)
colnames(p)
p <- p[,c(1,5,8,9,10,11)]
head(p)

#assoc file
a <- read.table(file="Assoc.assoc",header=TRUE,stringsAsFactors=F)
head(a)
nrow(f)

#fam file
fam <- read.table("ASwithUveitisCasesandControls.fam")
head(fam)

#the money
new <- f[,c(2,5)]
colnames(new) <- c("SNP","MAF_ALL")
head(new)
new$CA.fr <- a[,5]
new$CO.fr <- a[,6]
head(new)

Nsamples <- nrow(fam)
Ncases <- nrow(fam[fam$V6==2,])
Ncontrols <- nrow(fam[fam$V6==1,])
v <- Ncases/Nsamples
head(new)
head(p)

new1 <- merge(new,p,by.x="SNP",by.y="SNP")
head(new1)

new1$BETA <-as.numeric(new1$SNPWeight)/(sqrt(2*as.numeric(new1$MAF_ALL)*(1-as.numeric(new1$MAF_ALL))))
new1$OR <- 1 + new1$BETA*new1$MAF_ALL*(1-new1$MAF_ALL)/(v*(1-v)*new1$CO.fr*(1-new1$CA.fr))
head(new1)

#Confidence Intervals
new1$CorrSE <- as.numeric(new1$SNPWeightSE)*(sqrt(Nsamples))
new1$LL <- new1$OR - (1.96*new1$CorrSE)
new1$UL <- new1$OR + (1.96*new1$CorrSE)

head(new1)
new1[grep(96150086,new1[,"SNP"]),]
new1[grep(39388614,new1[,"SNP"]),]


S1 <- new1[grep(96150086,new1[,"SNP"]),]
S2 <- new1[grep(39388614,new1[,"SNP"]),]

S3 <- rbind(S1,S2)
colnames(S3) <- colnames(S1)
write.table(S3,"ASwithUveitisRESULT.txt",col.names=T,row.names=T,quote=F)
