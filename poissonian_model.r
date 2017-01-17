


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





sample.bams<-c("realign_cleaned_control.sort.bam","realign_cleaned_aml.sort.bam","realign_cleaned_relapse.sort.bam")
sample.names<-c("Normal","Diagnosis","Relapse")
names(sample.bams)<-sample.names[1:length(sample.bams)] # Normal MUST BE FIRST and in same order as 'sample.names"


tumor.contamination<-c(0,0.3,0.15)# 30% tumor at Diagnosis 85% at Relapse.
names(tumor.contamination)<-names(sample.bams)  ## sample.names<-c("Normal","Diagnosis","Relapse")


##### summary data of allele counts per SNP or mutation events vs indel .
summary<-read.delim("summary.FINAL_analysis.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
colnames(summary)[1:41]

#### loop over tumor classes
for( i in 2:length(names(sample.bams))){
 
 ##30% of reads will be normal subtract this expected amount first
 normal.allele.freq<-summary[,paste(sample.names[1],"Allele.Freq",sep=".")]

#"Normal.Allele.Freq" column is data (minor.allele.count/(total.allel.count))
# see the test below
 test<-summary[,paste(sample.names[1],"Minor.Allele.Count",sep=".")]/(summary[,paste(sample.names[1],"Major.Allele.Count",sep=".")]+summary[,paste(sample.names[1],"Minor.Allele.Count",sep=".")])

 test[1:10]
normal.allele.freq[1:10]


 ##### account for expect number of alleles from the Normal tissue in the tumor assuming tumor.contamination set above 
 
 normal.coverge.in.tumor<-(summary[,paste(sample.names[i],"Coverage",sep=".")]*tumor.contamination[sample.names[i]])
 normal.minor.in.tumor<-normal.coverge.in.tumor*normal.allele.freq
 normal.major.in.tumor<- normal.coverge.in.tumor-normal.minor.in.tumor ## count of major or genomic reference allele (forward strand genomic reference)
# normal.major.in.tumor+ normal.minor.in.tumor

 ### remove the reads from "Normal" tissue within the tumor sample:
 
 true.sample.cov<-summary[,paste(sample.names[i],"Coverage",sep=".")]- normal.coverge.in.tumor
 true.sample.major.counts<-summary[,paste(names(sample.bams)[i],"Major.Allele.Count",sep=".")] - normal.major.in.tumor
 true.sample.minor.counts<-summary[,paste(names(sample.bams)[i],"Minor.Allele.Count",sep=".")] - normal.minor.in.tumor
# true.sample.major.counts+ true.sample.minor.counts

 true.sample.freq<-true.sample.minor.counts/ true.sample.cov
 true.sample.freq[!is.finite(true.sample.freq)]<-0
 
 ## summary[,paste(sample.names[3],"Allele Freq",sep=".")]
 sort(normal.allele.freq)
 normal.allele.freq[normal.allele.freq<0.01]<-0.005 ### set a minimum low allele freqency so maths does not go crazy
 true.sample.cov[true.sample.cov==0]<-1  ## case when normal has coverge but tumor does not!
## Z=X-np/sqrt(np(1-p)) z-score of the minor allele frequency: (np is mean in Normal, X is measurement in tumor , the standard deviation is sqrt(np(1-p))
## normal.allele.freq<- 0.878/100
## true.sample.cov<- 1227 
##  true.sample.minor.counts<-   8
 p=2/66686
 n=
 z= (true.sample.minor.counts- true.sample.cov*normal.allele.freq) / (sqrt(true.sample.cov*normal.allele.freq*(1-normal.allele.freq)))
 p<-signif(2*(pnorm(z, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)),digits=3)        ##this is the 2 tailed
 z
 p
 true.sample.cov*normal.allele.freq

 ##### check this is what I have already calculated
 p[1:10]
summary[,paste(sample.names[i],"P.value",sep=".")][1:10]

##### $$$$ SO these  where the sorted in the columns 
 summary[,paste(sample.names[i],"P.value.NEW",sep=".")]<-p
 ###

                                }

0 summary[,paste(sample.names[1],"P.value.NEW",sep=".")]<-1 ### normal is the reference
 summary[,paste(sample.names[2],"P.value.NEW",sep=".")]
summary[,paste(sample.names[3],"P.value.NEW",sep=".")]
