#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(tidyr)

setwd(args[1])

data <- read.table(args[2],head=F,stringsAsFactors=F)
colnames(data) <- c("Chr","SNP","cm","Position","allele1","allele2")
num_line <- nrow(data)
line=paste("Read your SNP file that needs InDel allele update with ",num_line," lines.", sep="")
write(line,file="R_report_InDel_update.txt",append=TRUE)


# read the manifest file that will be used to update the missing alele
array_manifest <- read.delim(args[3], head=T,stringsAsFactors=F,sep="," , skip=7)
array_manifest_sub1 <- array_manifest[array_manifest$Name %in% data$SNP,]
num_lines <- nrow(array_manifest)
line=paste("Read the manifest file with ",num_lines," lines.",sep="")
write(line,file="R_report_InDel_update.txt",append=TRUE)

array_manifest_sub <- array_manifest_sub1[,c(2,10,11,18)]
array_manifest_sub$TopGenomicSeq <- gsub("^.*\\[","", array_manifest_sub$TopGenomicSeq)
array_manifest_sub$TopGenomicSeq <- gsub("\\].*$","", array_manifest_sub$TopGenomicSeq)

var_notin_array <- data[!data$SNP %in% array_manifest_sub$Name,]
write.table (var_notin_array[,2], "InDels_not_in_manifest_to_remove.txt", col.names=F, row.names=F, quote=F, sep="\t")
num_lines <- nrow(var_notin_array)
line=paste(num_lines," variants not in manifest file. The list can be found in InDels_not_in_manifest_to_remove.txt ",sep="")
write(line,file="R_report_InDel_update.txt",append=TRUE)

manifest_edit <- separate(data = array_manifest_sub, col = TopGenomicSeq, into = c("TopGenomic_alleleA_manifest", "TopGenomic_alleleB_manifest"), sep = "/")
colnames(manifest_edit) <- c("Name","Chr_manifest","MapInfo_manifest","TopGenomic_alleleA_manifest","TopGenomic_alleleB_manifest")


merge1 <- merge(data,manifest_edit, by.x = "SNP", by.y = "Name",all.x=T)
merge1_sub <- merge1[!is.na(merge1$TopGenomic_alleleA_manifest),]
num_lines <- nrow(merge1)
line=paste("Merged the SNP file with manifest file. The merged file has ", num_lines, " variants.",sep="")
write(line,file="R_report_InDel_update.txt",append=TRUE)


line="updating the InDel alleles...."
write(line,file="R_report_InDel_update.txt",append=TRUE)

update_indel <- merge1_sub
for (i in 1:nrow(merge1_sub)){
  # if allele1 is missing and allele 2 is present
  if(merge1_sub[i,5] == "I" & merge1_sub[i,6] == "D" & !is.na(merge1_sub[i,9]) & !is.na(merge1_sub[i,10])){ # if allele 1 is an insertion and allele 2 is a delition, and alleles A and B are not missing proceed:
    update_indel[i,6] <- "-" # assign missing to allele2
    update_indel[i,5] <- c(merge1_sub[i,9],merge1_sub[i,10])[(c(merge1_sub[i,9], merge1_sub[i,10]) != "-")] # allele 1 (insertion) whill the manifest allele that is not a dash
  }
  if (merge1_sub[i,5] == "D" & merge1_sub[i,6] == "I" & !is.na(merge1_sub[i,9]) & !is.na(merge1_sub[i,10])){
    update_indel[i,5] <- "-" # assign missing to allele1
    update_indel[i,6] <- c(merge1_sub[i,9],merge1_sub[i,10])[(c(merge1_sub[i,9], merge1_sub[i,10]) != "-")]
  }
  else {next} 
}


aa <- sum (update_indel[,5] == "I" | update_indel[,5] == "D" )
bb <- sum (update_indel[,6] == "I" | update_indel[,6] == "D" )

line=paste("Number of variants with allele1 codes as D/I after updating the InDel alleles: ", aa, sep="")
write(line,file="R_report_InDel_update.txt",append=TRUE)

line=paste("Number of variants with allele2 codes as D/I after updating the InDel alleles: ", bb, sep="")
write(line,file="R_report_InDel_update.txt",append=TRUE)


old_codes <- merge1_sub[,c(1,5,6)]
allele_change <- merge(old_codes,update_indel,by= "SNP",all.x=T)
allele_change1 <- allele_change[,c(1,2,3,7,8)]
write.table(allele_change1,"allele_change_Indels.txt", col.names=F, row.names=F, quote=F, sep="\t")












































