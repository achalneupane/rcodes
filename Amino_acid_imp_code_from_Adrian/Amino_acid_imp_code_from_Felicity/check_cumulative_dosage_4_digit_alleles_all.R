host = system("hostname",intern=T)

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
  

  sums = colSums(data2[,4:ncol(data2)])
  samples.out[[i]] = fam[which(sums > 2.5) , ]
}


cat("Done\n")

cat("length ", length(samples.out), "\n")

all = do.call('rbind',samples.out)
all.ids = unique(all$V1)



write.table(data.frame(V1=all.ids,V2=all.ids),
            file='data/dosage_samples_outliers_all.txt',
            row.names=F,col.names=F,quote=F)
