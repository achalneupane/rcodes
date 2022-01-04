#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(tidyr)

setwd(args[1])
    
# read the file with variants that have one allele as zero (missing)
data <- read.table(paste(args[2], sep=""),head=F,stringsAsFactors=F)
colnames(data) <- c("Chr","SNP","cm","Position","allele1","allele2")
num_line <- nrow(data)

line=paste("read your SNP file that needs allele update with ",num_line," lines.", sep="")
write(line,file="R_report_allele_update.txt",append=TRUE)


# read the manifest file that will be used to update the missing alele
array_manifest <- read.delim(paste(args[3], sep=""), head=T,,stringsAsFactors=T,sep="," , skip=7)
#colnames(array_manifest) <- c("IlmnID","Name","IlmnStrand","SNP","AddressA_ID","AlleleA_ProbeSeq","AddressB_ID","AlleleB_ProbeSeq","GenomeBuild","Chr","MapInfo","Ploidy","Species","Source","SourceVersion","SourceStrand","SourceSeq","TopGenomicSeq","BeadSetID","Exp_Clusters","Intensity_Only","RefStrand")

num_lines <- nrow(array_manifest)

line=paste("read the manifest file with ",num_lines," lines.",sep="")
write(line,file="R_report_allele_update.txt",append=TRUE)

array_manifest_sub1 <- array_manifest[array_manifest$Name %in% data$SNP,]

array_manifest_sub <- array_manifest_sub1[,c(2,4,10,11)]
array_manifest_sub$SNP <- gsub("\\[","", array_manifest_sub$SNP)
array_manifest_sub$SNP <- gsub("\\]","", array_manifest_sub$SNP)

manifest_edit <- separate(data = array_manifest_sub, col = SNP, into = c("alleleA_manifest", "allelB_manifest"), sep = "/")
colnames(manifest_edit) <- c("Name","alleleA_manifest","allelB_manifest","Chr_manifest","MapInfo_manifest")

merge1 <- merge(data,manifest_edit, by.x = "SNP", by.y = "Name",all.x=T)
num_lines <- nrow(merge1)
line=paste("Merged the SNP file with manifest file. The merged file has ", num_lines, " variants.",sep="")
write(line,file="R_report_allele_update.txt",append=TRUE)

# remove variants that do not have manifest allele info. if it is complete set no need to subset
merge1_sub <- merge1[!is.na(merge1$alleleA_manifest),]

num_lines <- nrow(merge1_sub)
line=paste("The merged file is filtered for variants that do not have manifest allele info. This resulted in ",num_lines," variants to be resolved.",sep="")
write(line,file="R_report_allele_update.txt",append=TRUE)


var_notin_array <- data[!data$SNP %in% array_manifest_sub$Name,]
write.table (var_notin_array[,2], "variants_not_in_manifest_to_remove.txt", col.names=F, row.names=F, quote=F, sep="\t")
num_lines <- nrow(var_notin_array)
line=paste(num_lines, " variants do not exist in manifest file.", sep="")
write(line,file="R_report_allele_update.txt",append=TRUE)



RC <- function(x){
  if (x == "A"){y <- "T"}
  if (x == "T"){y <- "A"}
  if (x == "G"){y <- "C"}
  if (x == "C"){y <- "G"}
  if (x == "I"){y <- "I"}
  if (x == "D"){y <- "D"}
  return(y)
}



#### allele update

#allele1 = col 5
#allele2 = col 6
#alleleA (manifest)= col 7
#alleleB (manifest) = col 8

### only update missing alleles (one of the alleles = 0, do not touch those that have alphabetic Allele1 and 2)
### NOTE: based on Excel data, it seems that the missing allele is always the first allele and I did not find any cases that allele 2 was zero. However the loop below is written for both

#variants with missing allele1
aa <- sum (merge1_sub[,5] == 0)

line=paste("num of variants with missing allele1: ",aa, sep="")
write(line,file="R_report_allele_update.txt",append=TRUE)


#variants with missing allele2
#sum (merge1_sub[,6] == 0)
#0

update_data <- merge1_sub

for (i in 1:nrow(merge1_sub)){
  # if allele1 is missing and allele 2 is present
  if(merge1_sub[i,5] == 0 & merge1_sub[i,6] !=0 & !is.na(merge1_sub[i,7]) & !is.na(merge1_sub[i,8])){ # if allele 1 is missing(0) but allele 2 is present, and alleles A and B are not missing proceed:
    if (merge1_sub[i,6] %in%  c( merge1_sub[i,7], merge1_sub[i,8]) ) # if allele2 is one of the AlleleA or AlleleB from manifest file
    {update_data [i,5] <- c(merge1_sub[i,7], merge1_sub[i,8])[!(c(merge1_sub[i,7], merge1_sub[i,8]) %in% merge1_sub[i,6])] } # Allele 1 (missing) should be the other allele (from alleles A and B)
    if (!merge1_sub[i,6] %in%  c( merge1_sub[i,7], merge1_sub[i,8]) & merge1_sub[i,6] %in%  c( RC(merge1_sub[i,7]), RC(merge1_sub[i,8]))) #if allele2 is NOT one of the AlleleA or AlleleB from manifest file, but present in alleles A and B reverse complements (RC), proceed:
    {update_data [i,5] <- RC(c(merge1_sub[i,7], merge1_sub[i,8])[!(c(merge1_sub[i,7], merge1_sub[i,8]) %in% RC(merge1_sub[i,6]))]) } # Allele 1 (missing) should be the RC of the other allele (from Allele A and B)
    if (!merge1_sub[i,6] %in%  c( merge1_sub[i,7], merge1_sub[i,8]) & !merge1_sub[i,6] %in%  c( RC(merge1_sub[i,7]), RC(merge1_sub[i,8]))) # if allele2 is not one of AlleleA or B or their RC, it means there is an issue and the alleles do not match. report this as a " No_match" character in allele to be investigated or removed
    {update_data [i,5] <- "No_match" }
  }
  
  else {next}
  
}

line="The alleles are updated"
write(line,file="R_report_allele_update.txt",append=TRUE)


#variants with missing allele1 in updated_data
bb <- sum (update_data[,5] == 0)

line=paste("num of variants with missing allele1 in updated file: ",bb, sep="")
write(line,file="R_report_allele_update.txt",append=TRUE)


cc <- sum (update_data[,5] == "No_match")
line=paste("num of variants with 'Not Matching' allele1 in updated file: ",cc, sep="")
write(line,file="R_report_allele_update.txt",append=TRUE)


old_codes <- merge1_sub[,c(1,5,6)]
allele_change <- merge(old_codes,update_data,by= "SNP",all.x=T)
allele_change1 <- allele_change[,c(1,2,3,7,8)]
write.table(allele_change1,"allele_change.txt", col.names=F, row.names=F, quote=F, sep="\t")

