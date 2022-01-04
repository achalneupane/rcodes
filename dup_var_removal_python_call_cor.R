#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(spgs)
library(ggplot2)

setwd(args[1])
data <- read.table(args[2],head=F,stringsAsFactors=F)
colnames(data) <- c("Chr","SNP","cm","Position","allele1","allele2", "dummy_ID")
num_line <- nrow(data)
line=paste("read variant file with dups with ",num_line," lines.", sep="")
write(line,file="R_report_dup_var_removal.txt",append=TRUE)


miss <- read.table(args[3], head=T,stringsAsFactors=F)
data_miss <- merge(data,miss, by.x ="SNP", by.y ="SNP")
num_line <- nrow(miss)
line=paste("read call rate (missingness) file with ",num_line," lines.", sep="")
write(line,file="R_report_dup_var_removal.txt",append=TRUE)
#order the file and write the SNP names to be used as index for VCF sorting
data_miss_order <- data_miss[order(data_miss$Chr, data_miss$dummy_ID,data_miss$allele1, data_miss$allele2),]
write.table(data_miss_order[,1], "Index_for_vcf_sorting.txt", col.names=F,row.names=F, quote=F, sep="\t")

line=("Calculating call correlation for pairwise duplicates...")
write(line,file="R_report_dup_var_removal.txt",append=TRUE)

##run python sript to get the call cors
#order the vcf file withour header based on the index generated anove
system("awk 'NR==FNR{o[FNR]=$1; next} {t[$3]=$0} END{for(x=1; x<=FNR; x++){y=o[x]; print t[y]}}' Index_for_vcf_sorting.txt duplicate_SNPs_noheader.vcf > vcf_index_sort.vcf")
#execute the call corr script in python3
system("cat vcf_index_sort.vcf | python3 call_corr_calc.py call_cor_pairwise.txt")

call_cor <- read.table("call_cor_pairwise.txt", head=F,stringsAsFactors=F )
colnames(call_cor)<- c("SNP","call_cor")

combo_order3 <- cbind(data_miss_order,call_cor$call_cor)


#### remove the duplicates ( variants with same pos and same alleles, and n% call corelation if no missing data in pair)
line=("start the loop to identify duplicate variants based on pos and allele...")
write(line,file="R_report_dup_var_removal.txt",append=TRUE)

# the changes are made in no_dups data frame
no_dups <- combo_order3

for(i in 2:nrow(combo_order3)){
  if (combo_order3[i,7] != combo_order3[(i-1),7]) # if the dummyIDs (chr_pos) for two consecutive rows are not the same, move to next row
  {next}
  if ((combo_order3[i,7] == combo_order3[(i-1),7]) & (combo_order3[i,5] == 0 | combo_order3[i,6] ==0 |combo_order3[(i-1),5] == 0 |combo_order3[(i-1),6] ==0)) #if the dummyIDs for two consecutive rows are the same, but one of the alleles in the two duplicat markers is zero, move to next row because we don't have enough info with missing allele to judge
  {next}
  if (combo_order3[i,7] == combo_order3[(i-1),7]) { # if the dummyIDs of two consecutive rows are the same and we know that none of the alleles are zero (previous conditional has filtered), then:
    if ( sum(c(combo_order3[i,5],combo_order3[i,6]) %in% c(combo_order3[(i-1),5],combo_order3[(i-1),6]),na.rm = TRUE) == 2 | sum(c(combo_order3[i,5],combo_order3[i,6]) %in% c(reverseComplement(combo_order3[(i-1),5],case="upper"),reverseComplement(combo_order3[(i-1),6],case="upper"),na.rm = TRUE) == 2)) # if the two alleles in row i exist in the alleles of row i-1 (duplicate vars) or exist in the RC of alleles of row i-1, then proceed. SO, one of the statements seperated by pipe has to be true. for each statement to be true, both alleles should exist among the alleles of the dup var, meaning the sum should yield 2 True values becasue sum() only outputs the sum of True values.
    { if (combo_order3[i,11] == combo_order3[(i-1),11] & combo_order3[(i),12] >= as.numeric(args[4])) # if none of the duplicates that have the same alleles (or RC), have any missng data, check the correlation between the calls. if they are more than 99% correlated, they are identical alleles, remove the ith row
    {no_dups[i,] <- NA}
      if (combo_order3[i,11] == combo_order3[(i-1),11] & combo_order3[(i),12] < as.numeric(args[4]) ) # if none of the duplicates that have the same alleles (or RC), have any missng data, check the correlation between the calls. if they are more than 99% correlated, they are identical alleles, remove the ith row
      {no_dups[i,] <- NA
       no_dups[(i-1),] <- NA}
      if (combo_order3[i,11] != combo_order3[(i-1),11]) # if the missing value is different between two duplicates:
      {aa <- which.max(c(combo_order3[i,11],combo_order3[(i-1),11])) # find which variant has the highest missing rate (lowwest call rate), and remove that row
      if (aa==1){no_dups[i,] <- NA} 
      if (aa==2){no_dups[(i-1),] <- NA}
      }
    }
    else {next} 
  }
  
}
#Num  f variants removed after the first loop:
mynum <- sum(is.na(no_dups$SNP))
line=paste("Variants removed after the first loop: ",mynum, sep="")
write(line,file="R_report_dup_var_removal.txt",append=TRUE)

# remove the rows with all NAs
no_dups2 <- no_dups %>% filter_all(any_vars(!is.na(.)))
bb <- nrow(no_dups2)
line=paste("Variants remained after the first loop: ",bb, sep="")
write(line,file="R_report_dup_var_removal.txt",append=TRUE)


############second loop
line=paste("The second loop for duplicate removal is started to remove duplicate left from triplicates and higher orders: ",bb, sep="")
write(line,file="R_report_dup_var_removal.txt",append=TRUE)

no_dups3 <- no_dups2

for (i in 2:nrow(no_dups2)){
  if (no_dups2[i,7] != no_dups2[(i-1),7]) # if the dummyIDs (chr_pos) for two consecutive rows are not the same, move to next row
  {next}
  if ((no_dups2[i,7] == no_dups2[(i-1),7]) & (no_dups2[i,5] == 0 | no_dups2[i,6] ==0 |no_dups2[(i-1),5] == 0 |no_dups2[(i-1),6] ==0)) #if the dummyIDs for two consecutive rows are the same, but one of the alleles in the two duplicat markers is zero, move to next row because we don't have enough info with missing allele to judge
  {next}
  if (no_dups2[i,7] == no_dups2[(i-1),7]) { # if the dummyIDs of two consecutive rows are the same and we know that none of the alleles are zero (previous conditional has filtered), then:
    if ( sum(c(no_dups2[i,5],no_dups2[i,6]) %in% c(no_dups2[(i-1),5],no_dups2[(i-1),6]),na.rm = TRUE) == 2 | sum(c(no_dups2[i,5],no_dups2[i,6]) %in% c(reverseComplement(no_dups2[(i-1),5],case="upper"),reverseComplement(no_dups2[(i-1),6],case="upper"),na.rm = TRUE) == 2)) # if the two alleles in row i exist in the alleles of row i-1 (duplicate vars) or exist in the RC of alleles of row i-1, then proceed. SO, one of the statements seperated by pipe has to be true. for each statement to be true, both alleles should exist among the alleles of the dup var, meaning the sum should yield 2 True values becasue sum() only outputs the sum of True values.
    { if (no_dups2[i,11] == no_dups2[(i-1),11] & no_dups2[(i),12] >= as.numeric(args[4])) # if none of the duplicates that have the same alleles (or RC), have any missng data, check the correlation between the calls. if they are more than 99% correlated, they are identical alleles, remove the ith row
    {no_dups3[i,] <- NA}
      if (no_dups2[i,11] == no_dups2[(i-1),11] & no_dups2[(i),12] < as.numeric(args[4])) # if none of the duplicates that have the same alleles (or RC), have any missng data, check the correlation between the calls. if they are more than 99% correlated, they are identical alleles, remove the ith row
      {no_dups3[i,] <- NA
      no_dups3[(i-1),] <- NA}
      if (no_dups2[i,11] != no_dups2[(i-1),11]) # if the missing value is different between two duplicates:
      {aa <- which.max(c(no_dups2[i,11],no_dups2[(i-1),11])) # find which variant has the highest missing rate (lowwest call rate), and remove that row
      if (aa==1){no_dups3[i,] <- NA} 
      if (aa==2){no_dups3[(i-1),] <- NA}
      }
    }
    else {next} 
  }
}


ee <- sum(is.na(no_dups3$SNP))
line=paste("Variants removed after the second loop: ", ee,sep="\t")
write(line,file="R_report_dup_var_removal.txt",append=TRUE)

# remove the rows with all NAs
no_dups4 <- no_dups3 %>% filter_all(any_vars(!is.na(.)))

# get the list of removed dups
lis_removed_dups <- combo_order3$SNP[!combo_order3$SNP %in% no_dups4$SNP]

write.table(lis_removed_dups,"dups_to_be_excluded.txt", row.names=F, col.names=F, quote=F, sep="\t")

line=("List of removed dup variants is written to file 'dups_to_be_excluded.txt'.")
write(line,file="R_report_dup_var_removal.txt",append=TRUE)

###check the higher order reps after removing variants with one missing allele
no_missing_allele <- no_dups4[no_dups4$allele1 != 0,]

no_miss_higher_order <- as.data.frame(table(no_missing_allele$dummy_ID))
colnames(no_miss_higher_order) <- c("dummyID","count")

no_miss_higher_order_IDs <- no_miss_higher_order[no_miss_higher_order$count > 2, ]
ff <- nrow(no_miss_higher_order_IDs)
line=paste("Number of triplicates and quadruplates after dup removal: ",ff,sep="")
write(line,file="R_report_dup_var_removal.txt",append=TRUE)


higher_order_rep_extract <- no_missing_allele[no_missing_allele$dummy_ID %in% no_miss_higher_order_IDs$dummyID,]
higher_order_rep_extract$call_cor <- NULL


higher_order_rep_extract_order <- higher_order_rep_extract[order(higher_order_rep_extract$dummy_ID ,higher_order_rep_extract$allele1, higher_order_rep_extract$allele2),]


write.table(higher_order_rep_extract_order, "Variants_with_triplicate_and_quadruplicate_after_dup_removal.txt", row.names=F, quote=F, sep="\t")

line="The list of variants with triplicates and quadruplicates is written in file Variants_with_triplicate_and_quadruplicate_after_dup_removal.txt. This file is expeted to be empty because plink does not consider multiallelic as dups, so we should not have multiallelic variants with the same posotion remaining after this script."
write(line,file="R_report_dup_var_removal.txt",append=TRUE)




######## plot the call rate for duplicates

#plot call rate for only dups (remove higher orders)
table_reps <- as.data.frame(table(combo_order3$dummy_ID),stringsAsFactors=F)
higher_order <- table_reps[table_reps$Freq > 2,]
colnames(higher_order) <- c("SNP","Freq")
sum(higher_order$Freq)

gg <- sum(higher_order$Freq)
line=paste("Number of triplicats and quadruplates excluded from plot data: ", gg, sep="")
write(line,file="R_report_dup_var_removal.txt",append=TRUE)

plot_data <- combo_order3[!combo_order3$dummy_ID %in% higher_order$SNP,]
write.table(plot_data,"duplicate_SNP_for_call_cor_plot.txt", col.names=F, row.names=F, quote=F, sep="\t")

system("cat duplicate_SNP_for_call_cor_plot.txt | python3 call_cor_plot_subset.py dup_subset_for_hist.txt")

hist_prep<- read.delim("dup_subset_for_hist.txt", head=F,stringsAsFactors=F )
colnames(hist_prep) <- colnames(plot_data)
colnames(hist_prep)[12] <- "call_cor"


pdf("hist_pairwise_call_cor.pdf")
ggplot(hist_prep, aes(x=call_cor)) + geom_histogram (bins=nrow(hist_prep))+ 
  scale_x_continuous(breaks=seq(0,100,5)) +
  theme(text = element_text(size=10),axis.text.x = element_text(angle=90, hjust=1)) 
dev.off()

# get a detailed frq table
p <- hist(hist_prep$call_cor, breaks=nrow(hist_prep)) 
myhist <- cbind.data.frame(p$mids,p$counts)
write.table(myhist,"freq_table_pairwise_call_cor.txt", col.names=T, row.names=F, quote=F, sep="\t")

line="A histogram of call correlation between pairwise duplicates with the same call rate is generated in the file hist_pairwise_call_cor.pdf. The frequency table is stored in freq_table_pairwise_call_cor.txt"
write(line,file="R_report_dup_var_removal.txt",append=TRUE)


