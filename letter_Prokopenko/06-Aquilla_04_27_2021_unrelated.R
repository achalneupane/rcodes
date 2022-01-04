## Remove preceding zeros from ADNI samples before matching the phenotype data. For example, ADNI_0002 should be ADNI_2


setwd("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/")

FAM <- read.table("AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1.fam", sep =" ")
dim(FAM)
FAM$V2 <- as.character(FAM$V2)



## cd fernandezv@virtual-workstation1:/gscmnt/gc2645/wgs/WXS_Aquilla/03-FINAL/VCFs/JOINT_CALL_7983$ 
## scp ID_list_brian.list achal@fenix.psych.wucon.wustl.edu://40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/
# PROJECTinfo <- read.delim("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/ID_list_brian.list", header = F, sep = " ", stringsAsFactors = FALSE)
temp <- "/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/ID_list_brian.list"
PROJECTinfo <- read.table(text = gsub(",", "\t", readLines(temp)))

PROJECTinfo$V2 <- as.character(PROJECTinfo$V2)
PROJECTinfo$PR <- gsub("/.*","", sapply(strsplit(PROJECTinfo$V2,split='^', fixed=TRUE), "[", 3))
PROJECTinfo <- PROJECTinfo[c(1,3)]
colnames(PROJECTinfo) <- c("SM", "PR_AN")


sum(FAM$V2 %in% PROJECTinfo$SM)
# > FAM$V2[!FAM$V2 %in% PROJECTinfo$SM]
# [1] "010-0009-0155433" "B6-48789"         "B6-48790"




FAM <- cbind(FAM, PROJECTinfo[match(FAM$V2, PROJECTinfo$SM), ])
sum(FAM$V2 == FAM$SM)


head(FAM)





PHENO <- read.csv("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/Aquilla_9515_phenotype-UPDATED.csv", header = TRUE, sep = ",", check.names = F)
PHENO$`Aquilla vcfID` <- as.character(PHENO$`Aquilla vcfID`)


dim(PHENO)
# View(PHENO)
FAM$IDs.with.no.Rep <- gsub("_R1$|_R2$", "", FAM$V2)




PHENO$ID_PR_PHENO <- paste(PHENO$`Aquilla vcfID`, PHENO$`Seq round / Project`, sep = ":")
PHENO$ID_PR_PHENO[grepl ("ADNI", PHENO$`Seq round / Project`)] <- paste0(paste(PHENO$`IID (Best ID)`, PHENO$`Seq round / Project`, sep = ":")[grepl ("ADNI", PHENO$`Seq round / Project`)], "_WGS")
FAM$ID_PR_FAM <- paste(FAM$IDs.with.no.Rep, FAM$PR_AN, sep =":")
sum(FAM$ID_PR_FAM %in% PHENO$ID_PR_PHENO)
sum(!FAM$ID_PR_FAM %in% PHENO$ID_PR_PHENO)


FAM_ADNI <- FAM[grepl("ADNI", FAM$V2),]
FAM_NON_ADNI <- FAM[!grepl("ADNI", FAM$V2),]
sum(!FAM_NON_ADNI$ID_PR_FAM %in% PHENO$ID_PR_PHENO)


MATCHED_non_ADNI <- FAM_NON_ADNI[FAM_NON_ADNI$ID_PR_FAM %in% PHENO$ID_PR_PHENO,]
NOTMATCHED_non_ADNI <- FAM_NON_ADNI[!FAM_NON_ADNI$ID_PR_FAM %in% PHENO$ID_PR_PHENO,]



FAM_ADNI$ID_PR_FAM <- gsub("_0+", "_", FAM_ADNI$V2)

MATCHED_ADNI <- FAM_ADNI[FAM_ADNI$ID_PR_FAM %in% PHENO$`IID (Best ID)`,]
NOTMATCHED_ADNI <- FAM_ADNI[!FAM_ADNI$ID_PR_FAM %in% PHENO$`IID (Best ID)`,]

MATCHED_non_ADNI <- cbind(MATCHED_non_ADNI, PHENO[match(MATCHED_non_ADNI$ID_PR_FAM, PHENO$ID_PR_PHENO),])
MATCHED_ADNI <- cbind(MATCHED_ADNI, PHENO[match(MATCHED_ADNI$ID_PR_FAM, PHENO$`IID (Best ID)`),])

head(MATCHED_non_ADNI)
head(MATCHED_ADNI)
MATCHED <- rbind.data.frame(MATCHED_ADNI, MATCHED_non_ADNI)

NOTMATCHED_non_ADNI

library(plyr)
MATCHED <- rbind.fill(MATCHED, NOTMATCHED_non_ADNI)

FAM_FINAL <- cbind.data.frame(FID=as.character(MATCHED$V1), IID=MATCHED$V2,
                              SEX=MATCHED$Sex, STATUS=MATCHED$STATUS,
                              APOE=MATCHED$APOE, AAO=MATCHED$AAO,
                              ALA=MATCHED$ALA, AGE=MATCHED$AGE, ADCO = MATCHED$ADCO,
                              WXS_TYPE= MATCHED$WXS,
                              PR= as.character(MATCHED$PR_AN))


## Fix the AGE information (If CA[2] use AAO )
FAM_FINAL$AGE [grepl("2", FAM_FINAL$STATUS)] <- FAM_FINAL$AAO [grepl("2", FAM_FINAL$STATUS)]
FAM_FINAL$AGE [grepl("1", FAM_FINAL$STATUS)] <- FAM_FINAL$ALA [grepl("1", FAM_FINAL$STATUS)]
FAM_FINAL$AGE [grepl("-9", FAM_FINAL$STATUS)] <- FAM_FINAL$ALA [grepl("-9", FAM_FINAL$STATUS)]



FAM_FINAL$WXS <- NA
FAM_FINAL$PR <- as.character(FAM_FINAL$PR)
unique(FAM_FINAL$PR)

FAM_FINAL$WXS [grepl("Genentech_WGS", FAM_FINAL$PR)] <- "2"
FAM_FINAL$WXS [grepl("Genentech_WES", FAM_FINAL$PR)] <- "1"
FAM_FINAL$WXS [grepl("Otogenetics_WES", FAM_FINAL$PR)] <- "1"
FAM_FINAL$WXS [grepl("TGI_WES", FAM_FINAL$PR)] <- "1"
FAM_FINAL$WXS [grepl("MGI_FASeEOAD_201605", FAM_FINAL$PR)] <- "1"
FAM_FINAL$WXS [grepl("Broad_WGS", FAM_FINAL$PR)] <- "2"
FAM_FINAL$WXS [grepl("Macrogen_WGS", FAM_FINAL$PR)] <- "2"
FAM_FINAL$WXS [grepl("LOAD_WES", FAM_FINAL$PR)] <- "1"
FAM_FINAL$WXS [grepl("phs000572_201508", FAM_FINAL$PR)] <- "1"
FAM_FINAL$WXS [grepl("MAPT_A152T", FAM_FINAL$PR)] <- "2"
FAM_FINAL$WXS [grepl("phs000572_201707", FAM_FINAL$PR)] <- "1"
FAM_FINAL$WXS [grepl("phs000572_201802", FAM_FINAL$PR)] <- "1"
FAM_FINAL$WXS [grepl("phs000572_201612", FAM_FINAL$PR)] <- "2"
FAM_FINAL$WXS [grepl("201907_USUHS_gDNA_SHERMAN", FAM_FINAL$PR)] <- "2"
FAM_FINAL$WXS [grepl("ADNI_WGS", FAM_FINAL$PR)] <- "2"
FAM_FINAL$WXS [grepl("MGI_Imaging_201612", FAM_FINAL$PR)] <- "1"
FAM_FINAL$WXS [grepl("MGI_Gregg_201704", FAM_FINAL$PR)] <- "1"
FAM_FINAL$WXS [grepl("MGI_DIANEXR_201706", FAM_FINAL$PR)] <- "2"
FAM_FINAL$WXS [grepl("CACHE_WGS", FAM_FINAL$PR)] <- "2"
FAM_FINAL$WXS [grepl("201909_MGI_gDNA_LINDSEY", FAM_FINAL$PR)] <- "2"
FAM_FINAL$WXS [grepl("Mayo_FTLD-TDP_201901", FAM_FINAL$PR)] <- "2"
FAM_FINAL$WXS [grepl("201812_MGI_DIANWGS_REDCLOUD", FAM_FINAL$PR)] <- "2"
FAM_FINAL$WXS [grepl("201907_USUHS_gDNA_SHERMA", FAM_FINAL$PR)] <- "2"
FAM_FINAL$WXS [grepl("MGI_DIANEXR_201902", FAM_FINAL$PR)] <- "2"

FAM_FINAL$WXS_TYPE <- as.character(FAM_FINAL$WXS_TYPE)
FAM_FINAL$WXS_TYPE [grepl("2", FAM_FINAL$WXS)] <- "WGS"
FAM_FINAL$WXS_TYPE [grepl("1", FAM_FINAL$WXS)] <- "WES"
unique(FAM_FINAL$WXS_TYPE)


df <- FAM_FINAL
tmpDF <- data.frame(matrix("", ncol = length(unique(df$PR)), nrow = nrow(df)))
colnames(tmpDF) <- unique(df$PR)
df <- cbind(df, tmpDF)

new.cols <- unique(as.character(df$PR))
df[, new.cols] <- "1"
for (i in 1:nrow(df)){
  df[i,grep(df$PR[i], colnames(df))] <- "2"
}


## Bias due to batch effect (complete phenotype is FASe_3894_Phenoscope.txt)
setwd("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files")
write.table(df, "FASe_2441_Phenoscope.txt", sep =" ", col.names = T, quote = F, row.names = FALSE)


IBD <- fread("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-IBD.genome", header = T, sep = " ")

# IBD_PROBLEM <- IBD[(IBD$Z0 < 0.5 &
#                       IBD$Z0 > 0.10 &
#                       IBD$Z1 < 0.9 &
#                       IBD$Z1 > 0.6),] 
# # IBD_PROBLEM <- IBD[(IBD$Z0 < 0.3), ]


IBD$COLOR <- with(IBD, ifelse((Z1 + Z0) < 0.98, "Expected",
                              ifelse(0.1 < Z0 & Z0 < 0.33 & 0.6 < Z1 & Z1 < 0.93, "Unexpected", "Expected")))

sum(grepl ("Unexpected", IBD$COLOR) )

IBD_PROBLEM <- IBD[grepl ("Unexpected", IBD$COLOR),]


IBD_PROBLEM <- cbind(IBD_PROBLEM, IID1_PR = FAM_FINAL$PR[match(IBD_PROBLEM$IID1, FAM_FINAL$IID)], IID2_PR = FAM_FINAL$PR[match(IBD_PROBLEM$IID2, FAM_FINAL$IID)])
write.table(IBD_PROBLEM, "Problem_IBDs.csv", sep =",", col.names = T, quote = F, row.names = FALSE)


Problem_samples_PR <- cbind.data.frame(as.vector(rbind(IBD_PROBLEM$IID1, IBD_PROBLEM$IID2)), as.vector(rbind(IBD_PROBLEM$IID1_PR, IBD_PROBLEM$IID2_PR)))
colnames(Problem_samples_PR) <- c("SM", "PR")
library(dplyr)
Problem_samples_PR <- distinct(Problem_samples_PR,SM, .keep_all= TRUE)


# IBD$COLOR <- with(IBD, ifelse((Z1 + Z0) < 0.95, "Other",
#                               ifelse(0.1 < Z0 & Z0 < 0.5 & 0.5 < Z1 & Z1 < 0.9, "Unexpected", "Expected")))
# ggplot(IBD, aes(Z0, Z1)) +
#   geom_point(aes(color = COLOR)) +
#   scale_color_manual(values = c(Expected="blue", Unexpected="red", Expected="blue"))



P <- ggplot(IBD, aes(Z0, Z1)) +
  geom_point(aes(color = COLOR)) +
  scale_color_manual(values = c(Expected="black", Unexpected="red"))  

ggsave("Replication dataset-2441-IBD_Problem.jpg", plot = P, device = NULL, scale = 1, width = 16, height = 9, dpi = 300, limitsize = TRUE)

PROBLEM_SAMPLES <- unique(c(IBD_PROBLEM$IID1, IBD_PROBLEM$IID2))
write.table(PROBLEM_SAMPLES, "IBD_Problem_Samples.csv", sep =",", col.names = T, quote = F, row.names = FALSE)

## Problem samples that are several times repeated
# IID1	n	PR	
# 62033	2010	MGI_FASeEOAD_201605	
# 64069	393	Broad_WGS	
# MAP_63762	380	201909_MGI_gDNA_LINDSEY	
# MAP_65789	137	201909_MGI_gDNA_LINDSEY

## duplicates and problem SM list
duplicates <- IBD[(IBD$Z0 < 0.25 &
                     IBD$Z1 < 0.25), ]

FOUR_PROBLEM_SM <- c("62033", "64069", "MAP_63762", "MAP_65789")
# samples with discordant sex
DIS_SEX <- c("MAP_22392", "MAP_62727", "MAP_63762", "MAP_64559")

PROB_SMs <- unique(c(FOUR_PROBLEM_SM, DIS_SEX, duplicates$IID1))

PHENO <- read.table("FASe_2441_Phenoscope.txt", header = T)

PROB_SMs <- PHENO[match(PROB_SMs, PHENO$IID), c("FID", "IID")]
write.table(PROB_SMs, "Duplicates_and_problem_Samples.txt", col.names = T, quote = F, row.names = FALSE)


## Read the heterogygosity
HETR <- read.delim("AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS_QC_het.het", header = T, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
head(HETR)
dim(HETR)
# Num_Sites - Num Homozygous/ Num_Sites
HETR$HETEROZYGOSITY_RATE <- (HETR$`N(NM)` - HETR$`O(HOM)`)/HETR$`N(NM)`
# HETR$HETROZYGOSITY_RATE
HETR$SM <- as.numeric(rownames(HETR))
P <- ggplot(HETR, aes(x=SM, y=HETEROZYGOSITY_RATE)) + 
  geom_point() + theme(text = element_text(size=16))
ggsave("Replication dataset-2441-Heterozygosity_rate.jpg", plot = P, device = NULL, scale = 1, width = 8, height = 6, dpi = 300, limitsize = TRUE)

breaks <- c(0,0.05,0.1,0.2)
tags <- c("[0-0.05]","[0.05-0.1]", "[0.1-0.2]")
# bucketing HETR bins
HET_GROUPS <- cut(HETR$HETEROZYGOSITY_RATE, 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=tags)
# inspect HET GROUPS
summary(HET_GROUPS)

## Check in this file to see if the problem IDs have similar issues in Aquilla call:
## /100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-a/01-Aquilla-preQC/03-PLINK-QC-files_all_Aquilla/Aquilla_7983_WXS_SNPS-INDELS-geno0.05-hwe-mind0.1-WXSmissingCLEAN-sexupdate-IBD.genome

library(data.table)
IBD_Aquilla <- fread("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-a/01-Aquilla-preQC/03-PLINK-QC-files_all_Aquilla/Aquilla_7983_WXS_SNPS-INDELS-geno0.05-hwe-mind0.1-WXSmissingCLEAN-sexupdate-IBD.genome", header = T, sep = "\t")
IBD_Aquilla <- IBD_Aquilla[,c(2,4,7,8,9,10)]

IBD_Aquilla$KEY <- paste(IBD_Aquilla$IID1, IBD_Aquilla$IID2, sep = ":")

IBD_PROBLEM <- IBD[grepl ("Unexpected", IBD$COLOR),]
IBD_PROBLEM$KEY1 <- paste(IBD_PROBLEM$IID1, IBD_PROBLEM$IID2, sep = ":")
IBD_PROBLEM$KEY2 <- paste(IBD_PROBLEM$IID2, IBD_PROBLEM$IID1, sep = ":")

IBD_Aquilla$MATCH1 <- ifelse(IBD_Aquilla$KEY %in% IBD_PROBLEM$KEY1, "YES", "NO")
IBD_Aquilla$MATCH2 <- ifelse(IBD_Aquilla$KEY %in% IBD_PROBLEM$KEY2, "YES", "NO")

IBD_Aquilla$MATCH <- paste0(IBD_Aquilla$MATCH1, IBD_Aquilla$MATCH2)
IBD_Aquilla$TYPE <- ifelse(grepl("YES", IBD_Aquilla$MATCH), "MATCH", "NOMATCH")
head(IBD_Aquilla)
IBD_Aquilla_small <- IBD_Aquilla[,c("Z0", "Z1", "PI_HAT", "TYPE")]

library(ggplot2)
library(plyr)
P <- ggplot(IBD_Aquilla_small, aes(Z0, Z1)) +
  geom_point(aes(color = TYPE)) +
  geom_point(data = subset(IBD_Aquilla_small, TYPE == 'MATCH'),  aes(x = Z0, y = Z1, color = TYPE)) +
  scale_color_manual(values = c(NOMATCH="black", MATCH="red"))  


ggsave("Problem_IBD_pairs_matched_in_full_Aquilla.jpg", plot = P, device = NULL, scale = 1, width = 16, height = 9, dpi = 300, limitsize = TRUE)



## Check heterozygosity of all Aquilla
## Read the heterogygosity
HETR_Aquilla <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-a/01-Aquilla-preQC/03-PLINK-QC-files_all_Aquilla/Aquilla_7983_WXS_SNPS-INDELS-geno0.05-hwe-mind0.1-WXSmissingCLEAN-sexupdate_QC_het.het", header = T, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
head(HETR_Aquilla)
dim(HETR_Aquilla)

# Num_Sites - Num Homozygous/ Num_Sites
HETR_Aquilla$HETEROZYGOSITY_RATE <- (HETR_Aquilla$`N(NM)` - HETR_Aquilla$`O(HOM)`)/HETR_Aquilla$`N(NM)`
HETR_Aquilla$SM <- as.numeric(rownames(HETR_Aquilla))

HETR_Aquilla$Problem_samples <- ifelse(HETR_Aquilla$IID %in% PROBLEM_SAMPLES, "YES", "NO")

HETR_Aquilla$Problem_samples <- factor(HETR_Aquilla$Problem_samples, levels = c("YES", "NO"))

library(plyr)
P <- ggplot(HETR_Aquilla) +
  geom_point(aes(x = SM, y = HETEROZYGOSITY_RATE)) +
  geom_point(data = subset(HETR_Aquilla, Problem_samples == 'YES'),  aes(x = SM, y = HETEROZYGOSITY_RATE, color = Problem_samples)) +
  scale_color_manual(values = c(NO="black", YES="blue")) 

ggsave("Aquilla-7983-Heterozygosity_rate.jpg", plot = P, device = NULL, scale = 1, width = 8, height = 6, dpi = 300, limitsize = TRUE)



##################### Vicky's code to plot heterozygosity
## replication dataset
HET<-read.table("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS_QC_het.het", head=T)
library(ggplot2)
# first plot Observed vs Expected Heterozigosity
P <- ggplot(HET, aes(x=O.HOM.,y=E.HOM.)) + geom_point() + ggtitle("Replication_set_2441 Heterozygosity")
ggsave("AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS_QC_het.jpg", plot=P, device = NULL, scale = 1, width = 10, height = 10, dpi = 300, limitsize = TRUE)
# then plot the p-values to see "normality"
P <- ggplot(HET,aes(x=F)) + geom_histogram(binwidth=0.005) + ggtitle("Replication_set_2441 Fvalues")
ggsave("AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS_QC_het.fpvals.jpg", plot=P, device = NULL, scale = 1, width = 10, height = 10, dpi = 300, limitsize = TRUE)


## Aquilla
HET<-read.table("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-a/01-Aquilla-preQC/03-PLINK-QC-files_all_Aquilla/Aquilla_7983_WXS_SNPS-INDELS-geno0.05-hwe-mind0.1-WXSmissingCLEAN-sexupdate_QC_het.het", head=T)
library(ggplot2)
# first plot Observed vs Expected Heterozigosity
P <- ggplot(HET, aes(x=O.HOM.,y=E.HOM.)) + geom_point() + ggtitle("Aquilla_7882 Heterozygosity")
ggsave("Aquilla_7983_WXS_SNPS-INDELS-geno0.05-hwe-mind0.1-WXSmissingCLEAN-sexupdate.het.jpg", plot=P, device = NULL, scale = 1, width = 10, height = 10, dpi = 300, limitsize = TRUE)
# then plot the p-values to see "normality"
P <- ggplot(HET,aes(x=F)) + geom_histogram(binwidth=0.005) + ggtitle("Aquilla_7882 Heterozygosity Fvalues")
ggsave("Aquilla_7983_WXS_SNPS-INDELS-geno0.05-hwe-mind0.1-WXSmissingCLEAN-sexupdate.het.fpvals.jpg", plot=P, device = NULL, scale = 1, width = 10, height = 10, dpi = 300, limitsize = TRUE)

####################

## Replot IBD after removing sex discordant and problem samples

# After removing 234 Duplicates_and_problem_Samples.txt
IBD <- fread("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean1-IBD.genome", header = T, sep = " ")
P <- ggplot(IBD, aes(Z0, Z1)) +
  geom_point()
ggsave("AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean1-IBD.jpg", plot = P, device = NULL, scale = 1, width = 16, height = 9, dpi = 300, limitsize = TRUE)

## Now, remove one of the two related samples
sum(IBD$Z0 < 0.38)

IBD_related2 <- IBD[IBD$Z0 < 0.39,]
IBD_UNRELATED <- IBD[!IBD$Z0 < 0.39,]

# ggplot(IBD_UNRELATED, aes(Z0, Z1)) +
#   geom_point()

Remove_relatives <- IBD_related2[,c("FID1", "IID1")]
# > dim(Remove_relatives)
# [1] 126   2
colnames(Remove_relatives) <- c("FID", "IID")


## Also remove any ADSPs
FAM <- read.table("AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean1.fam", header = FALSE, check.names = FALSE)
head(FAM)
FAM <- cbind.data.frame(FAM, PHENO[match(FAM$V2, PHENO$IID),])
dim(FAM)

ADSP <- FAM[grepl("phs", FAM$PR),]
ADSP <- ADSP[,1:2]
dim(ADSP)
colnames(ADSP) <- c("FID", "IID")

samples_in_discovery <- read.table("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/FBAT/FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1_with_STATUS_nonADSP_post_QC2-geno-0.02-maxmaf-0.01.fam", header = FALSE, sep = " ")
samples_in_discovery <- samples_in_discovery[1:2]
colnames(samples_in_discovery) <- c("FID", "IID")
samples_in_discovery <- FAM[na.omit(match(samples_in_discovery$IID, FAM$V2)), c("FID", "IID")]

## list of samples to remove: 1 relatives, 2 ADSP, 3 samples common in discovery set
Remove_relatives_and_ADSPs <- rbind.data.frame(Remove_relatives, ADSP, samples_in_discovery)
Remove_relatives_and_ADSPs <- distinct(Remove_relatives_and_ADSPs)


write.table(Remove_relatives_and_ADSPs, "Relatives_and_problem_Samples_samples_in_discovery.txt", col.names = T, sep = "\t", quote = F, row.names = FALSE)


# Now plot the IBD
IBD <- fread("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean2-IBD.genome", header = T, sep = " ")
P <- ggplot(IBD, aes(Z0, Z1)) +
  geom_point()
ggsave("AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean2-IBD.jpg", plot = P, device = NULL, scale = 1, width = 16, height = 9, dpi = 300, limitsize = TRUE)

## HAPMAP PCA
PCA <- read.table("AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean2-HAPMAP-MERGED3-for_PCA.eigenvec", header =T, stringsAsFactors=FALSE)
HAPMAP.ethnicty <- read.table("relationships_w_pops_121708.txt", header = T )
head(HAPMAP.ethnicty)

PCA$COHORT <- "FASe"
PCA$COHORT <- HAPMAP.ethnicty$population[match(PCA$IID, HAPMAP.ethnicty$IID)]
PCA <- PCA[c(1:4,23)]
PCA$COHORT <- as.character(PCA$COHORT)
PCA$COHORT[is.na(PCA$COHORT)] <- "FASe"
write.table(PCA, "FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic_post_QC_with_STATUS-hwe-geno0.05-mind0.1_post_QC2-HAPMAP-MERGED3-for_PCA.eigenvec-PC1-PC2-COHORT.txt", sep ="\t", col.names = T, quote = F)


library(ggplot2)
PCAs<- read.table("AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean2-HAPMAP-MERGED3-for_PCA.eigenvec-PC1-PC2-COHORT.txt", header=T)
# ##plotting:
# ggplot(PCAs, aes(x=PC1, y=PC2, color=COHORT)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Replication 1780") +
#   scale_color_manual(values = c('green', 'black', 'red', "blue"))
# ggsave("Replication_1780-PCs-COHORT.jpg", plot = last_plot(), device = NULL, scale = 1, width = 8, height = 4, dpi = 600, limitsize = TRUE)

## Plot PCA with ordered CEU and Replication
library(ggplot2)
library(plyr)
ggplot(PCAs) +
  geom_point(aes(x=PC1, y=PC2, color=COHORT)) +
  geom_point(data = subset(PCAs, COHORT == 'CEU'),  aes(x = PC1, y = PC2, color = COHORT)) +
  scale_color_manual(values = c('green', 'black', 'red', "blue"))
ggsave("Replication_1780-PCs-COHORT.jpg", plot = last_plot(), device = NULL, scale = 1, width = 8, height = 4, dpi = 600, limitsize = TRUE)

PCAs$COHORT <- as.character(PCAs$COHORT)
REPLICATION <- PCAs[grepl ("Replication", PCAs$COHORT),]
sum(PCAs$PC1 > 0.00 & PCAs$PC2 < 0.125)

REPLICATION <- REPLICATION[REPLICATION$PC1 > 0.00 & REPLICATION$PC2 < 0.125,]
ggplot(REPLICATION) +
  geom_point(aes(x=PC1, y=PC2, color=COHORT)) +
  geom_point(data = subset(REPLICATION, COHORT == 'CEU'),  aes(x = PC1, y = PC2, color = COHORT)) +
  scale_color_manual(values = c('green', 'black', 'red', "blue"))


## Extract 1590 NHW samples
write.table(REPLICATION[1:2], "REPLICATION_Data_1590_NHW.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)

REPLICATION.PHENO <- read.table("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/Replication_1590_Phenoscope.txt", header = T)


## Plot PCA of 1590 NHW with HAPMAP 
PCA <- read.table("AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean3-HAPMAP-MERGED3-for_PCA.eigenvec", header =T, stringsAsFactors=FALSE)
HAPMAP.ethnicty <- read.table("relationships_w_pops_121708.txt", header = T )
head(HAPMAP.ethnicty)

PCA$COHORT <- "Replication"
PCA$COHORT <- HAPMAP.ethnicty$population[match(PCA$IID, HAPMAP.ethnicty$IID)]
PCA <- PCA[c(1:4,23)]
PCA$COHORT <- as.character(PCA$COHORT)
PCA$COHORT[is.na(PCA$COHORT)] <- "Replication"
write.table(PCA, "AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean3-HAPMAP-MERGED3-for_PCA.eigenvec-PC1-PC2-COHORT.txt", sep ="\t", col.names = T, quote = F)





library(ggplot2)
PCAs<- read.table("AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean3-HAPMAP-MERGED3-for_PCA.eigenvec-PC1-PC2-COHORT.txt", header=T)
##plotting:
ggplot(PCAs) +
  geom_point(aes(x=PC1, y=PC2, color=COHORT)) + xlab("PC1") + ylab("PC2") + ggtitle("Replication 1590") +
  geom_point(data = subset(PCAs, COHORT == 'CEU'),  aes(x = PC1, y = PC2, color = COHORT)) +
  scale_color_manual(values = c('green', 'black', 'red', "blue"))
ggsave("Replication_1590_NHW_with_HAPMAP.jpg", plot = last_plot(), device = NULL, scale = 1, width = 8, height = 4, dpi = 600, limitsize = TRUE)



## PCA of 1590 without HAPMAP
PCAs<- read.table("AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean3-PCAS.eigenvec", header=T)
##plotting:
library(ggplot2)
ggplot(PCAs, aes(x=PC1, y=PC2)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Replication 1590")
# ggsave("Replication_NHW_1590-PCs-COHORT_without_HAPMAP.jpg", plot = last_plot(), device = NULL, scale = 1, width = 6, height = 4, dpi = 600, limitsize = TRUE)


## add PCs to PHENO
PHENO <- read.table("FASe_2441_Phenoscope.txt", header = T, check.names = FALSE)

head(PCAs)
colnames(PCAs)[1:2] <- paste0("PC_", colnames(PCAs)[1:2])
PHENO.FINAL <- cbind.data.frame(PHENO[match(PCAs$PC_IID, PHENO$IID),], PCAs[,c(3:12)])

write.table(PHENO.FINAL, "Replication_1590_Phenoscope.txt", sep =" ", col.names = T, quote = F, row.names = FALSE)

write.table(PHENO.FINAL[,1:2], "Replication_1590_samples.txt", sep =" ", col.names = T, quote = F, row.names = FALSE)







### Quick Manhattan
setwd("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files")
library (qqman)
library(data.table)
LOGISTIC<- fread("Manhattan_AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean3.assoc.logistic",head=T)

jpeg("QQ_Manhattan_AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean3.assoc.logistic.jpg", units="mm", width=190, height=142, res=1000)
qq(LOGISTIC$P)
dev.off()

library (qqman)
LOGISTIC<- fread("Manhattan_AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean3.assoc.logistic",head=T)
jpeg("Manhattan_Manhattan_AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean3.assoc.logistic.jpg", units="mm", width=190, height=142, res=1000)
manhattan(LOGISTIC, main = "", ylim=c(0,20), col = c("blue4", "orange3"), suggestiveline = -log10(1e-06), genomewideline = -log10(1e-08), annotateTop = FALSE, annotatePval = 1e-06, chrlabs = as.character(1:24)) 
dev.off()


##################################################################



anno_file <- fread("AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean3-snpeff-dbnsfp-ExAC.0.3.GRCh38.vcf.tsv", header = T, sep = "\t")

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


# jpeg(paste0("QQ_plot_", PValueFILE, ".jpeg"), height = 20, width = 15, units='cm', res = 300)
jpeg(paste0("POSTER_QQ_plot_", PValueFILE, ".jpeg"), height = 8, width = 8, units='cm', res = 300, pointsize = 22)
qqunif.plot(my.pvalues, LAMBDA= LAMBDA)
dev.off()




# LOGISTIC
# LOGISTIC$CHR <- sapply(strsplit(LOGISTIC$Marker,":"), `[`, 1)
# LOGISTIC$BP <- sapply(strsplit(LOGISTIC$Marker,":"), `[`, 2)
# LOGISTIC$SNP <- LOGISTIC$Marker

LOGISTIC$CHR[LOGISTIC$CHR == "X"] <- 23
LOGISTIC$CHR[LOGISTIC$CHR == "Y"] <- 24
LOGISTIC$CHR <- as.numeric(LOGISTIC$CHR)
LOGISTIC$BP <- as.numeric(LOGISTIC$BP)

# jpeg(paste0("Manhattan_", PValueFILE, ".jpeg"), height = 20, width = 30, units='cm', res = 300, pointsize = 14)
jpeg(paste0("POSTER_Manhattan_", PValueFILE, ".jpeg"), height = 15, width = 30, units='cm', res = 300, pointsize = 19)
# LOGISTIC <- LOGISTIC[!grepl("23|24",LOGISTIC$CHR ),]
# LOGISTIC.t <- LOGISTIC
# LOGISTIC.t <- LOGISTIC.t[!is.na(LOGISTIC.t$CHR),]
# LOGISTIC.t$SNP <- as.character(gsub(".*_","",LOGISTIC.t$SNP))
# ###

# manhattan(LOGISTIC, main = "", ylim=c(0,22), cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = -log10(1e-06), genomewideline = -log10(1e-08), annotatePval = -log10(1e-08), chrlabs = as.character(c(1:23)))
manhattan(LOGISTIC, main = "", ylim=c(0,10), col = c("blue4", "orange3"), suggestiveline = -log10(1e-06), genomewideline = -log10(1e-08), annotateTop = FALSE, annotatePval = 1e-06, chrlabs = as.character(1:23))
dev.off()

########################################## MAF <= 0.01
library (qqman)
LOGISTIC<-fread("Manhattan_LOGISTIC_AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean3-geno-0.02-maxmaf-0.01.assoc.logistic",head=T)

my.pvalues <- LOGISTIC$P[!is.na(LOGISTIC$P)]

jpeg("QQ_Manhattan_LOGISTIC_AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean3-geno-0.02-maxmaf-0.01.assoc.logistic2.jpg", units="mm", width=190, height=142, res=1000)
qqunif.plot(my.pvalues, LAMBDA= LAMBDA)
dev.off()



anno_file <- fread("AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean3-snpeff-dbnsfp-ExAC.0.3.GRCh38.vcf.tsv", header = T, sep = "\t")

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


# LOGISTIC
# LOGISTIC$CHR <- sapply(strsplit(LOGISTIC$Marker,":"), `[`, 1)
# LOGISTIC$BP <- sapply(strsplit(LOGISTIC$Marker,":"), `[`, 2)
# LOGISTIC$SNP <- LOGISTIC$Marker

LOGISTIC$CHR[LOGISTIC$CHR == "X"] <- 23
LOGISTIC$CHR[LOGISTIC$CHR == "Y"] <- 24
LOGISTIC$CHR <- as.numeric(LOGISTIC$CHR)
LOGISTIC$BP <- as.numeric(LOGISTIC$BP)

jpeg("Manhattan_Manhattan_LOGISTIC_AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean3-geno-0.02-maxmaf-0.01.assoc.logistic.jpg", height = 15, width = 30, units='cm', res = 300, pointsize = 19)
# manhattan(LOGISTIC, main = "", ylim=c(0,22), cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = -log10(1e-06), genomewideline = -log10(1e-08), annotatePval = -log10(1e-08), chrlabs = as.character(c(1:23)))
manhattan(LOGISTIC, main = "", ylim=c(0,10), col = c("blue4", "orange3"), suggestiveline = -log10(1e-06), genomewideline = -log10(1e-08), annotateTop = FALSE, annotatePval = 1e-06, chrlabs = as.character(1:23))
dev.off()


##############################














































#########################END #######################################
write.table(LOGISTIC, "Firth-Fallback_replication_study_results.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)



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
