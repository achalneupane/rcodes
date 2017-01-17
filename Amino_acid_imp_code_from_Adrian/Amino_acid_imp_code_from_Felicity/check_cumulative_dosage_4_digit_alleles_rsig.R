host = system("hostname",intern=T)
host
if(host == "irving") {
  this.dir="/home/fsnewell/Desktop/cervical_cancer_hla"  
} else {
  this.dir="/home/tfnewell/PBSHOME/Projects/Cervical_Cancer/data/merged/merge_info795/association/cervical_cancer_hla_genotyped"
}
setwd(this.dir) 

fam.file = "data/cervical_cancer_forward_qc_snp2hla_merged_hg19.fam"
col1 = read.table("data/col1.dosage")
dose.file = "data/cervical_cancer_forward_qc_snp2hla_merged_hg19.dosage"


loci = c("HLA_A","HLA_C","HLA_B","HLA_DRB1",
         "HLA_DQA1","HLA_DQB1","HLA_DPA1","HLA_DPB1")

samples.out = list()

r2.file="data/cervical_cancer_forward_qc_snp2hla_merged_hg19.r2"
r2.data <- as.matrix(read.delim(r2.file, sep="\t"))
#Get min r squared

r2.data.min<-cbind(r2.data[,1],apply(r2.data, 1, min))
#Get rid of NA
posns.all <- !is.na(r2.data.min[,2]) 
r2.data.min=r2.data.min[posns.all,]
#Get values above 0.5
posns.wanted=r2.data.min[,2] > 0.5

r2.wanted <- r2.data.min[posns.wanted,]

fam = read.table(fam.file,header=F)

for( i in 1:length(loci)) {
  cat(i,"\n")
  locus = loci[i]
  rows = grep(locus,col1$V1)
  tmp.file <- tempfile()
  cmd <- paste("sed -n ",min(rows),",",max(rows),
               "p ",dose.file," > ",tmp.file,sep='')
  system(cmd,intern=FALSE)
  data = read.table(tmp.file,header=F,sep='\t',
    colClasses=c(rep("character",3),rep("numeric",nrow(fam))))
  allele.names = col1[rows,]
  allele.names.char<-as.character(allele.names)
  ids.codes = sapply(strsplit(allele.names.char,split="_"),function(x) return(x[3]))

  are.4.digit = nchar(ids.codes) == 4
  data2 = data[are.4.digit,]
  ##keep those with high R2
  
  data2.wanted<- data2[,1] %in% r2.wanted[,1]
  data3<-data2[data2.wanted,]

  sums = colSums(data3[,4:ncol(data3)])
  samples.out[[i]] = fam[which(sums > 2.5) , ]
}


cat("Done\n")

cat("length ", length(samples.out), "\n")

all = do.call('rbind',samples.out)
all.ids = unique(all$V1)

dose = read.table(dose.file,header=F)
dose.r2 = dose[,1] %in% r2.wanted[,1]

dosage.keep = dose[dose.r2,]

write.table(dosage.keep,"test.r2.sig.dosage",
            sep="\t", col.names=F, row.names=F, quote=F)


write.table(data.frame(V1=all.ids,V2=all.ids),
            file='data/dosage_samples_outliers_rsig.txt',
            row.names=F,col.names=F,quote=F)
