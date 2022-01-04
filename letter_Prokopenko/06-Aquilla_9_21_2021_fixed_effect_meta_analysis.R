
### Quick Manhattan
setwd("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/fixed_effect_meta-analysis")
library (qqman)
library(data.table)
LOGISTIC<- fread("METAANALYSIS_inputfile1-2_1.tbl",head=T)
LOGISTIC$CHR <- as.numeric(as.character(sapply(strsplit(LOGISTIC$MarkerName,":"), `[`, 1)))
LOGISTIC$BP <- as.numeric(as.character(sapply(strsplit(LOGISTIC$MarkerName,":"), `[`, 2)))

colnames(LOGISTIC) <- c("SNP","Allele1" , "Allele2", "Weight", "Zscore", "P", "Direction", "CHR", "BP")
sort(unique(LOGISTIC$CHR))
LOGISTIC <- LOGISTIC[!is.na(LOGISTIC$CHR),]


jpeg("QQ_AQUILLA_fixed_effect_meta_analysis.jpg", units="mm", width=190, height=142, res=1000)
qq(LOGISTIC$P)
dev.off()

## Read annotation file1
anno_file1 <- fread("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean3-snpeff-dbnsfp-ExAC.0.3.GRCh38.vcf.tsv", header = T, sep = "\t")

## Read Anotation file2
anno_file2 <- fread("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/FBAT/FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1_with_STATUS_nonADSP_post_QC2-snpeff-dbnsfp-ExAC.0.3.GRCh38.vcf.tsv", header = T, sep = "\t")

anno_file <- rbind.data.frame(anno_file1, anno_file2)
anno_file <- anno_file[-which(duplicated(anno_file)), ]

## Add gene annotation
anno_file$key <- sapply(strsplit(anno_file$ID,";"), `[`, 1)
head(anno_file)
## Annotations split by | will give: 2 consequences, 3 gene, region type (codng/non coding), 5 gene region
anno_file$consequence <- sapply(strsplit(anno_file$INFO,"\\|"), `[`, 2)
anno_file$gene <- sapply(strsplit(anno_file$INFO,"\\|"), `[`, 3)
anno_file$type <- sapply(strsplit(anno_file$INFO,"\\|"), `[`, 4)
anno_file$region <- sapply(strsplit(anno_file$INFO,"\\|"), `[`, 5)
## No drop INFO field
anno_file <- as.data.frame(anno_file)
anno_file <- anno_file[,-grep("INFO",colnames(anno_file))]
# anno_file[, -which(names(anno_file) %in% "INFO")]


# https://stackoverflow.com/questions/13774773/check-whether-values-in-one-data-frame-column-exist-in-a-second-data-frame
LOGISTIC <- cbind(LOGISTIC,anno_file[match(LOGISTIC$SNP, anno_file$key),])
LOGISTIC_ANNO <- LOGISTIC


library(lattice)
source("https://raw.githubusercontent.com/achalneupane/rcodes/master/qqrunif_with_lambda.r")
my.pvalues <- LOGISTIC$P


# This is a helper functoin to calculate LAMBDA
inflation <- function(pvalues) {
  # chisq <- qchisq(1-my.pvalues, 1, lower.tail = F)
  chisq <- qchisq(1 - pvalues, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  return(lambda)
}

# library("GenABEL")
# estlambda(my.pvalues, method="median")


LAMBDA <- round(inflation(my.pvalues), digits = 3)
LAMBDA


PValueFILE <- "METAANALYSIS_inputfile1-2_1.tbl"


### MAX MAF 0.01
# jpeg(paste0("QQ_plot_", PValueFILE, ".jpeg"), height = 20, width = 15, units='cm', res = 300)
jpeg(paste0("QQ_plot_", PValueFILE, ".jpeg"), units="mm", width=190, height=142, res=1000)
qqunif.plot(my.pvalues, LAMBDA= LAMBDA)
dev.off()




# LOGISTIC
LOGISTIC$CHR <- sapply(strsplit(LOGISTIC$SNP,":"), `[`, 1)
LOGISTIC$BP <- sapply(strsplit(LOGISTIC$SNP,":"), `[`, 2)

LOGISTIC$CHR[LOGISTIC$CHR == "X"] <- 23
LOGISTIC$CHR[LOGISTIC$CHR == "Y"] <- 24
LOGISTIC$CHR <- as.numeric(LOGISTIC$CHR)
LOGISTIC$BP <- as.numeric(LOGISTIC$BP)

LOGISTIC$SNP <- paste(LOGISTIC$SNP, LOGISTIC$gene, sep = "_")

# jpeg(paste0("Manhattan_", PValueFILE, ".jpeg"), height = 20, width = 30, units='cm', res = 300, pointsize = 14)
jpeg(paste0("Manhattan_", PValueFILE, ".jpeg"), height = 15, width = 30, units='cm', res = 300, pointsize = 19)
# LOGISTIC <- LOGISTIC[!grepl("23|24",LOGISTIC$CHR ),]
# LOGISTIC.t <- LOGISTIC
# LOGISTIC.t <- LOGISTIC.t[!is.na(LOGISTIC.t$CHR),]
# LOGISTIC.t$SNP <- as.character(gsub(".*_","",LOGISTIC.t$SNP))
# ###

# manhattan(LOGISTIC, main = "", ylim=c(0,22), cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = -log10(1e-06), genomewideline = -log10(1e-08), annotatePval = -log10(1e-08), chrlabs = as.character(c(1:23)))
manhattan(LOGISTIC, main = "", ylim=c(0,10), col = c("blue4", "orange3"), suggestiveline = -log10(1e-06), genomewideline = -log10(1e-08), annotateTop = FALSE, annotatePval = 1e-05, chrlabs = as.character(1:22))
dev.off()


write.table(LOGISTIC_ANNO, "Fixed_effect_Meta_analysis_results.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)


#####################################################

#####################################################


































## Match with Tanzi replicates
Discovery_replicated_Tanzi <- read.table("Discovery_replicated_Tanzi.csv", sep =",", header = TRUE)


## Replicated genes from Discovery
Discovery_replicated_Tanzi$KEY_VAR <- paste(Discovery_replicated_Tanzi$Chromosome, Discovery_replicated_Tanzi$Position, sep =":")
Discovery_replicated_Tanzi$KEY <- paste(Discovery_replicated_Tanzi$Chromosome, Discovery_replicated_Tanzi$Position, Discovery_replicated_Tanzi$Effect_allele, Discovery_replicated_Tanzi$Other_allele, sep =":")
LOGISTIC_ANNO$KEY_VAR <- paste(LOGISTIC_ANNO$CHR, LOGISTIC_ANNO$BP, sep =":")


sum(Discovery_replicated_Tanzi$KEY_VAR %in% LOGISTIC_ANNO$KEY_VAR)
# 18
sum(Discovery_replicated_Tanzi$KEY %in% LOGISTIC_ANNO$key)
# 0

Discovery_replicated_Tanzi$FASE_MATCH_by_site_P <- LOGISTIC_ANNO$P[match(Discovery_replicated_Tanzi$KEY_VAR, LOGISTIC_ANNO$KEY_VAR)]

## Match variants from Single-variant-analysis 
require('gtools')
library(stringr)

Discovery_replicated_Tanzi$FASE_MATCH_by_nearest_gene_P <- ""
Discovery_replicated_Tanzi$FASE_MATCH_by_nearest_gene_P_sorted <- ""

Discovery_replicated_Tanzi$Nearest_protein_coding_gene <- as.character(Discovery_replicated_Tanzi$Nearest_protein_coding_gene)
for(i in 1:nrow(Discovery_replicated_Tanzi)){
MATCH_GENE <- Discovery_replicated_Tanzi$Nearest_protein_coding_gene[i]
if (sum(grepl(paste0("^",MATCH_GENE,"$"), LOGISTIC_ANNO$gene)) < 1){
  next
}
print(paste0("DOING gene: ", MATCH_GENE, ", ROW num ", i))
wanted.var.p.values <- paste0(LOGISTIC_ANNO$SNP[grepl(paste0("^",MATCH_GENE,"$"), LOGISTIC_ANNO$gene)], ", P=", LOGISTIC_ANNO$P [grepl(paste0("^",MATCH_GENE,"$"), LOGISTIC_ANNO$gene)])
# No need to show vars with p= NA
wanted.var.p.values <- wanted.var.p.values [! grepl("P=NA", wanted.var.p.values)]
# sort them by p-values
wanted.var.p.values <- mixedsort(wanted.var.p.values)

Discovery_replicated_Tanzi$FASE_MATCH_by_nearest_gene_P[i] <- paste(wanted.var.p.values, collapse = ";")
wanted.var.p.values.sorted <- wanted.var.p.values[mixedorder(gsub('.*P=(.*)','\\1',wanted.var.p.values))]
Discovery_replicated_Tanzi$FASE_MATCH_by_nearest_gene_P_sorted[i] <- paste(wanted.var.p.values.sorted, collapse = ";")
}

# NOTE: Genes replicated in Prokopenko et al is indicated by REPLICATED column.
# Genes replicated between Prokopenko et al. and FASe data is represented by
# column FASE_MATCH_by_nearest_gene_P

write.table(Discovery_replicated_Tanzi, "Single_Variant_Discovery_replicated_Prokopenko_replicated_also_in_FASe.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)





###########################################################################################
########################## MATCH VARIANTS FROM SPATIAL clustering #########################
###########################################################################################




Spatial_replicated_Tanzi <- read.table("Spatial_replicated_Tanzi.csv", sep =",", header = TRUE)

singlevar_Tanzi$Nearest.protein.coding.gene <- as.character(singlevar_Tanzi$Nearest.protein.coding.gene)

# sum(LOGISTIC_ANNO$gene %in% singlevar_Tanzi$Nearest.protein.coding.gene)


LOGISTIC_ANNO <- cbind(LOGISTIC_ANNO, singlevar_Tanzi[match(LOGISTIC_ANNO$gene, singlevar_Tanzi$Nearest.protein.coding.gene), c("Nearest.protein.coding.gene", "VEP.consequence")])
REPLICATED_singled_var <- LOGISTIC_ANNO[!is.na(LOGISTIC_ANNO$Nearest.protein.coding.gene),]

colnames(REPLICATED_singled_var) <- c("Marker", "Allele", "afreq", "fam_size", "S-E(S)", "Var(S)", "Z", "P", "CHROM", "POS", "ID", "REF", "ALT", "key", "consequence", "gene", "type", "region", "TANZI_Nearest.protein.coding.gene", "TANZI_VEP.consequence")
write.table(REPLICATED_singled_var, "/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/FBAT/Replicated_single_variant_genes.csv", sep ="\t", col.names = T, quote = F, row.names = FALSE)

## below suggestive significance
singlevar_Tanzi$key <- paste(singlevar_Tanzi$Chromosome, singlevar_Tanzi$Position, sep =":")
REPLICATED_singled_var_below_suggestive <- REPLICATED_singled_var[REPLICATED_singled_var$P < 0.05,]
write.table(REPLICATED_singled_var_below_suggestive, "/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/FBAT/REPLICATED_singled_var_below_suggestive.csv", sep ="\t", col.names = T, quote = F, row.names = FALSE)

SUGGESTIVE_SINGLE_VAR <- {}
for (i in 1:nrow(singlevar_Tanzi)){
  SUGGESTIVE_SINGLE_VAR_tmp <-  REPLICATED_singled_var_below_suggestive[grepl(singlevar_Tanzi$key[i], REPLICATED_singled_var_below_suggestive$key),]
  SUGGESTIVE_SINGLE_VAR <- rbind.data.frame(SUGGESTIVE_SINGLE_VAR, SUGGESTIVE_SINGLE_VAR_tmp)
}


library("GenomicRanges")
q=GRanges(seqnames=LOGISTIC_ANNO$`#CHROM`,
          ranges=IRanges(start = LOGISTIC_ANNO$POS, end = LOGISTIC_ANNO$POS)
)

q

ALL_overlappingVAR <- {}
for (i in 1:nrow(singlevar_Tanzi)){
  gr=GRanges(seqnames=singlevar_Tanzi$Chromosome[i],
             ranges=IRanges(start = singlevar_Tanzi$Position[i]-100000, end = singlevar_Tanzi$Position[i]+100000))
  overlappingVAR <- subsetByOverlaps(q, gr)
  overlappingVAR <- as.data.frame(overlappingVAR)
  if(nrow(overlappingVAR) != 0){
    overlappingVAR$TANZI_var <- paste(singlevar_Tanzi$Chromosome[i], singlevar_Tanzi$Position[i], sep =":")
    ALL_overlappingVAR <- rbind.data.frame(ALL_overlappingVAR, overlappingVAR)
  }
}


ALL_overlappingVAR$key <- paste(ALL_overlappingVAR$seqnames, ALL_overlappingVAR$start, sep = ":")
hundreadKB_SV <- {}
for (i in 1:nrow(ALL_overlappingVAR)){
  hundreadKB_SV_tmp <- LOGISTIC_ANNO [grepl(ALL_overlappingVAR$key[i], LOGISTIC_ANNO$key),]  
  hundreadKB_SV_tmp$TANZI_var <- ALL_overlappingVAR$TANZI_var[i]
  hundreadKB_SV <- rbind.data.frame(hundreadKB_SV, hundreadKB_SV_tmp)
}


write.table(hundreadKB_SV, "/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/FBAT/Replicated_single_variant_genes_within_100KB.csv", sep ="\t", col.names = T, quote = F, row.names = FALSE)

## below suggestive significance
hundreadKB_SV$key2 <- paste(hundreadKB_SV$`#CHROM`, hundreadKB_SV$POS, sep =":")
hundreadKB_SV$P <- as.numeric(hundreadKB_SV$P)
hundreadKB_SV_below_suggestive <- hundreadKB_SV[hundreadKB_SV$P < 0.05,]

write.table(hundreadKB_SV_below_suggestive, "/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/FBAT/hundreadKB_SV_below_suggestive.csv", sep ="\t", col.names = T, quote = F, row.names = FALSE)

SUGGESTIVE_SINGLE_VAR_100KB <- {}
for (i in 1:nrow(singlevar_Tanzi)){
  SUGGESTIVE_SINGLE_VAR_100KB_tmp <-  hundreadKB_SV_below_suggestive[grepl(singlevar_Tanzi$key[i], hundreadKB_SV_below_suggestive$key2),]
  SUGGESTIVE_SINGLE_VAR_100KB <- rbind.data.frame(SUGGESTIVE_SINGLE_VAR_100KB, SUGGESTIVE_SINGLE_VAR_100KB_tmp)
}












########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
## Start demographic table
demographic_1590 <- PHENO.FINAL
table(table(as.character(demographic_1590$FID)))
# > table(table(as.character(demographic_1803$FID)))

table(table(as.character(demographic_1590$FID)))
# > table(table(as.character(demographic_1803$FID)))
# 1   2   3   4   5   6   7   8   9 
# 365  11  95 128  59  28  15   3   3 

########################################################################################################
####################### Age-of-onset stratified by STATUS, ethnicity and sex ###########################
########################################################################################################

covars <- demographic_1590
colnames(covars)[colnames(covars) == "AAO"] <- "AGE_AT_ONSET"
colnames(covars)[colnames(covars) == "ALA"] <- "AGE_LAST_VISIT"
covars$ETHNICITY <- "NHW"

## Set binary values to APOE 
covars$APOE4ANY <- covars$APOE
sum(grepl("22|23|33", covars$APOE4ANY))
covars$APOE4ANY[grepl("22|23|33|32", covars$APOE4ANY)] <- 0
covars$APOE4ANY[grepl("24|34|44|42|43", covars$APOE4ANY)] <- 1
table(covars$APOE4ANY)
table(covars$STATUS)


library(dplyr)
df <- covars %>%
  mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999, 0), "unknown"))) %>%
  group_by(STATUS, ETHNICITY) %>%
  summarise('<=65' = sum(as.numeric(as.character(AGE_AT_ONSET)) <= 65 & as.numeric(as.character(AGE_AT_ONSET)) > 0 ,
                         na.rm = TRUE),
            '<=70'= sum(as.numeric(as.character(AGE_AT_ONSET)) <= 70  & as.numeric(as.character(AGE_AT_ONSET)) > 0,
                        na.rm = TRUE),
            '<=75'= sum(as.numeric(as.character(AGE_AT_ONSET)) <= 75  & as.numeric(as.character(AGE_AT_ONSET)) > 0,
                        na.rm = TRUE))

df <- as.data.frame(df)
df
## require(tidyverse)
# df <- df %>% pivot_wider(Ethnicity, names_from = STATUS, values_from = c(`<65`,`<70`,`<75`))
library(reshape2)
# summarize values
# df <- melt(df, id.vars = 1:2, measure.vars = 3:6)
df <- reshape2::melt(df, id.vars = 1:2, measure.vars = 3:5)
df <- reshape2::dcast(df, ETHNICITY ~ STATUS + variable)
# colnames(df) <- gsub(".*_", "", colnames(df))
df
###############################################################
################# age at last assessment ######################
###############################################################

library(dplyr)
df2 <- covars %>%
  mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999), NA))) %>%
  group_by(STATUS, ETHNICITY) %>%
  summarise('>=70'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 70,
                        na.rm = TRUE),
            '>=80'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 80,
                        na.rm = TRUE),
            '>=85'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 85,
                        na.rm = TRUE))


df2 <- as.data.frame(df2)
## require(tidyverse)
# df <- df %>% pivot_wider(Ethnicity, names_from = STATUS, values_from = c(`<65`,`<70`,`<75`))
library(reshape2)
# summarize values
# df2 <- melt(df2, id.vars = 1:2, measure.vars = 3:6)
df2 <- reshape2::melt(df2, id.vars = 1:2, measure.vars = 3:5)
df2 <- reshape2::dcast(df2, ETHNICITY ~ STATUS + variable)
# colnames(df) <- gsub(".*_", "", colnames(df))
df2

Age_stratified <- cbind(df,df2)
Age_stratified
write.csv(Age_stratified, "stratified_classification.csv", quote = FALSE)



FIX_weird_NUMS <- function(x){
  x <- as.numeric(as.character(x))
  x [x %in% c(-9, 888, 999, 0)] <- NA
  return(x)
}

covars[grepl("AGE",colnames(covars))] <- sapply(covars[grepl("AGE",colnames(covars))], FIX_weird_NUMS)


# N Control Missing age
MISSING <- sum(is.na(covars$AGE_LAST_VISIT[covars$STATUS == 1]))
MISSING
# N Cases Missing age
MISSING <- sum(is.na(covars$AGE_AT_ONSET[covars$STATUS == 2]))
MISSING





table(APOE=covars$APOE4ANY, STATUS=covars$STATUS)

## Status: control = 1, case = 2, MCI = 3
## Sex: 1 = male, 2 = female, -9 = unknown
df <- as.data.frame(table(covars$SEX, covars$ETHNICITY)[2:3,])
df <- rbind(df,colSums(df))
df
# % Percent of females:
(df[2,]/df[3,])*100

## Total
table(covars$ETHNICITY)

## % APOE4
# df <- table(covars$APOE4ANY, covars$ETHNICITY)[2:3,]
# df <- rbind(df,colSums(df))
# df
# (df[2,]/df[3,])*100
(table( covars[, c("APOE4ANY", "ETHNICITY") ] )[3,]/ as.vector(aggregate(covars$IID, list(covars$ETHNICITY), length)[2])) * 100

# Number of cases Controls by Ethnicity
table(covars$STATUS, covars$ETHNICITY)

# % APOE4+ for each ethnicity
# table(covars$APOE4ANY, covars$STATUS, covars$ETHNICITY)
# table(STATUS=NHW$status==1, APOE=NHW$apoe4any==1 )

# % CONTROLS APOE4+
(table( covars[ covars$STATUS == "1" , c("APOE4ANY", "ETHNICITY") ] )[3,]/table( covars[ covars$STATUS == "1" , c("ETHNICITY") ] )) * 100

# % CASES APOE4+
(table( covars[ covars$STATUS == "2" , c("APOE4ANY", "ETHNICITY") ] )[3,]/table( covars[ covars$STATUS == "2" , c("ETHNICITY") ] )) * 100

# % MCI APOE4+
(table( covars[ covars$STATUS == "3" , c("APOE4ANY", "ETHNICITY") ] )[3,]/table( covars[ covars$STATUS == "3" , c("ETHNICITY") ] )) * 100
################################### 10/27/2020 Vicky only wants age of onset for Cases and MCI and Age at last assessment for CO
# CO > 75, 80 and 85

attr(covars, "spec") <- NULL
class(covars) <- setdiff(class(covars), "spec_tbl_df")

FIX_weird_NUMS <- function(x){
x <- as.numeric(as.character(x))
x [x %in% c(-9, 888, 999, 0)] <- NA
return(x)
}

covars[grepl("AGE",colnames(covars))] <- sapply(covars[grepl("AGE",colnames(covars))], FIX_weird_NUMS)

# Average age (Cases==2)
MEAN_CA <- mean(as.numeric(as.character(covars[covars$STATUS == 2,"AGE_AT_ONSET"])), na.rm= T)
SD_CA <- sd(as.numeric(as.character(covars[covars$STATUS == 2,"AGE_AT_ONSET"])), na.rm= T)
RANGE_CA <- range(as.numeric(as.character(covars[covars$STATUS == 2,"AGE_AT_ONSET"])), na.rm= T)


# Average age (Controls)
MEAN_CO <- mean(as.numeric(as.character(covars[covars$STATUS == 1,"AGE_AT_ONSET"])), na.rm= T)
SD_CO <- sd(as.numeric(as.character(covars[covars$STATUS == 1,"AGE_AT_ONSET"])), na.rm= T)
RANGE_CO <- range(as.numeric(as.character(covars[covars$STATUS == 1,"AGE_AT_ONSET"])), na.rm= T)


# End demographic table
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################

# For Gene-based analysis, run rscript SKAT-unrelated-loop-OPTIMAL-cohort-covariate.R
########################################################################
