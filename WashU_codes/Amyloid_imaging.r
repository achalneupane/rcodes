###########################################################################################
############################       ANALYSIS      ##########################################
###########################################################################################

setwd("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/06-Amyloid_Imaging-Achal")
AMYLOID <- read.table("/80/GWAS/amyloid_imaging_GWAS/amyloid_7_cohorts_pheno_complete.txt", header = TRUE, colClasses=c('character'))
head(AMYLOID)

AMYLOID$sample_id <- as.character(AMYLOID$sample_id)
AMYLOID$CLEANED_ID <- as.character(AMYLOID$sample_id)



WXS_FAM <- read.delim("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/06-Amyloid_Imaging-Achal/Bloomfield-WES-WGS-project.csv", header = T, colClasses=c('character'), sep = ",")
head(WXS_FAM)
dim(WXS_FAM)


WXS_FAM$Bloomfield.gvcf.id..SM...9810. <- as.character(WXS_FAM$Bloomfield.gvcf.id..SM...9810.)
WXS_FAM$CLEANED_ID <- WXS_FAM$Bloomfield.gvcf.id..SM...9810.
sum(WXS_FAM$CLEANED_ID %in% AMYLOID$sample_id)
# 330

## CREATE a FOUND and NOT FOUND column
AMYLOID$FOUND <- ifelse (AMYLOID$CLEANED_ID %in% WXS_FAM$CLEANED_ID, "YES", "NO")

# Check How many overlap
table(AMYLOID$FOUND )
# NO  YES 
# 5997  330

## MAP samples
sum(grepl("^MAP_", WXS_FAM$CLEANED_ID))
# 2227

###################
## CLEAN MAP IDs ##
###################
## GWAS MAP IDs are missing MAP_ prefix, so fix that
AMYLOID$CLEANED_ID[(grepl("^MAP_", AMYLOID$IID))] <- paste0("MAP_", AMYLOID$CLEANED_ID[(grepl("^MAP_", AMYLOID$IID))])

## RUN FOUND not Found
AMYLOID$FOUND <- ifelse (AMYLOID$CLEANED_ID %in% WXS_FAM$CLEANED_ID, "YES", "NO")
# Check How many overlap
table(AMYLOID$FOUND)
# NO  YES 
# 5352  975 

####################
## CLEAN ADNI IDs ##
####################
## GWAS ADNI IDs are missing ADNI_ prefix, so fix that
AMYLOID$CLEANED_ID[(grepl("ADNI_CL", AMYLOID$study))] <- sapply(strsplit(AMYLOID$IID[(grepl("ADNI_CL", AMYLOID$study))], split="\\^"), "[", 1)

## RUN FOUND not Found
AMYLOID$FOUND <- ifelse (AMYLOID$CLEANED_ID %in% WXS_FAM$CLEANED, "YES", "NO")
# Check How many overlap
table(AMYLOID$FOUND)
# NO  YES 
# 4660 1667

####################
## CLEAN DIAN IDs ##
####################
AMYLOID$CLEANED_ID[grepl("DIAN", AMYLOID$IID)] <- gsub("DIAN_", "", sapply(strsplit(AMYLOID$IID[grepl("DIAN", AMYLOID$IID)], split = "\\^"), "[", 1))
## RUN FOUND not Found
AMYLOID$FOUND <- ifelse (AMYLOID$CLEANED_ID %in% WXS_FAM$CLEANED, "YES", "NO")
# Check How many overlap
table(AMYLOID$FOUND)
# NO  YES 
# 4660 1667




#########################################################################
##                 Do the same but for WXS data                        ##
#########################################################################
# Object AMYLOID is from above
head(WXS_FAM)

sum(WXS_FAM$CLEANED_ID %in% AMYLOID$CLEANED_ID)
# 1667

## CREATE a FOUND and NOT FOUND column
WXS_FAM$FOUND <- ifelse (WXS_FAM$CLEANED_ID %in% AMYLOID$CLEANED_ID, "YES", "NO")

# Check How many overlap
table(WXS_FAM$FOUND)
# NO  YES 
# 8143 1667 

## MAP samples
sum(grepl("^MAP_", WXS_FAM$CLEANED))
# 2227

##############################################
## CLEAN duplicate .WGS or .WES samples IDs ##
##############################################
WXS_FAM$CLEANED_ID[(grepl(".WGS|.WES", WXS_FAM$CLEANED_ID))] <- gsub("\\..*","",WXS_FAM$CLEANED_ID[(grepl(".WGS|.WES", WXS_FAM$CLEANED_ID))])

## RUN FOUND not Found
WXS_FAM$FOUND <- ifelse (WXS_FAM$CLEANED_ID %in% AMYLOID$CLEANED_ID, "YES", "NO")
# Check How many overlap
table(WXS_FAM$FOUND)
# NO  YES 
# 7701 2109



# Check again to see which ones are not found in Muhammed's AMYLOID data
## RUN FOUND not Found
AMYLOID$FOUND <- ifelse (AMYLOID$CLEANED_ID %in% WXS_FAM$CLEANED_ID, "YES", "NO")
# Check How many overlap
table(AMYLOID$FOUND)
# NO  YES 
# 4442 1885

####################
## CLEAN DIAN IDs ##
####################

WXS_FAM$CLEANED_ID[grepl("[[:digit:]]+-+[[:digit:]]+-+[[:digit:]]", WXS_FAM$BEST_ID..9810.) & grepl("DIAN", WXS_FAM$Seq_project)
                   & !grepl("EXR",  WXS_FAM$BEST_ID..9810.)] <- sapply(strsplit(WXS_FAM$CLEANED_ID[grepl("[[:digit:]]+-+[[:digit:]]+-+[[:digit:]]", WXS_FAM$BEST_ID..9810.) & grepl("DIAN", WXS_FAM$Seq_project)
                                                                                                   & !grepl("EXR",  WXS_FAM$BEST_ID..9810.)], split="-"), "[", 3)

## RUN FOUND not Found
WXS_FAM$FOUND <- ifelse (WXS_FAM$CLEANED_ID %in% AMYLOID$CLEANED_ID, "YES", "NO")
# Check How many overlap
table(WXS_FAM$FOUND)
# NO  YES 
# 7480 2330 

ALL_SEQ_FOUND <- WXS_FAM

WXS_FAM <- WXS_FAM[grepl("YES", WXS_FAM$FOUND),]
WXS_FAM_B4_QC <- WXS_FAM

# Run again on AMYLOID to see which ones are not found
## RUN FOUND not Found
AMYLOID$FOUND <- ifelse (AMYLOID$CLEANED_ID %in% WXS_FAM$CLEANED_ID, "YES", "NO")
# Check How many overlap
table(AMYLOID$FOUND)
# NO  YES 
# 4222 2105


ALL_AMYLOID <- AMYLOID

## We do not have A4_SUVR yet
AMYLOID <- AMYLOID[!grepl("A4_SUVR", AMYLOID$study),]


# Keep only those IDs that overlap with the latest FAM post IBD
BLOOMFIELD_FAM <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/03-PLINK-QC-files/Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2.fam", header = F )
BLOOMFIELD_FAM$V2 <- as.character(BLOOMFIELD_FAM$V2)

## Exclude
WXS_FAM$Bloomfield.gvcf.id..SM...9810.[!(WXS_FAM$Bloomfield.gvcf.id..SM...9810. %in% BLOOMFIELD_FAM$V2)]

WXS_FAM <- WXS_FAM[(WXS_FAM$Bloomfield.gvcf.id..SM...9810. %in% BLOOMFIELD_FAM$V2),]
dim(WXS_FAM)
# 1999
# ADD FID from FAM file
WXS_FAM <- cbind.data.frame(setNames(BLOOMFIELD_FAM[match(WXS_FAM$Bloomfield.gvcf.id..SM...9810., BLOOMFIELD_FAM$V2),1:2], c("FID", "IID")), WXS_FAM)
sum(WXS_FAM$IID == WXS_FAM$Bloomfield.gvcf.id..SM...9810.)
# 1999


write.table(WXS_FAM[1:2], "/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/06-Amyloid_Imaging-Achal/amyloid_imaging_sample_list.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)

#############################################################################
################################ IBD ########################################
#############################################################################
## IBD
library(data.table)
IBD <- fread("Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_Amyloid_Imaging-IBD.genome")

library(ggplot2)
p <- ggplot(IBD, aes(x=Z0, y=Z1))+ geom_point() + ggtitle("Bloomfield-Amyloid_Imaging - 1999")
p + annotate("rect", xmin = 0.6, xmax = 1.02, ymin= -0.01, ymax= 0.4, 
             fill=NA, colour="red") +
  annotate("text", x=0.6, y=0.25, label="Unrelated", size=4, color = "black")
ggsave("Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_Amyloid_Imaging-IBD.jpg", plot = last_plot(), device = NULL, scale = 1, width = 8, height = 5, dpi = 300, limitsize = TRUE)




## IDs around 1 and 0.5 are possibly the related samples, so removing those samples as well
related <- IBD [IBD$Z1 > 0.4,]

related$nIID1 <- with(transform(related, n = 1),  ave(n, IID1, FUN = length))
related$nIID2 <- with(transform(related, n = 1),  ave(n, IID2, FUN = length))

# Between nIID1 and nIID2, find those with max repitions
related[, TO_REMOVE_FID := cbind(FID1, FID2)[cbind(.I, max.col(.SD, "first"))], .SDcols = nIID1:nIID2]
related[, TO_REMOVE_IID := cbind(IID1, IID2)[cbind(.I, max.col(.SD, "first"))], .SDcols = nIID1:nIID2]

write.table(cbind.data.frame(related$TO_REMOVE_FID, related$TO_REMOVE_FID), "/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/06-Amyloid_Imaging-Achal/remove_IDs-AN.list", sep ="\t", col.names = F, quote = F, row.names = FALSE)



## After removing related samples, perform IBD again. Also, do PCA with HAPMAP and plot that PCA
library(ggplot2)
IBD<-read.table("Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_Amyloid_Imaging-CLEAN1-IBD.genome", head=T)
ggplot(IBD, aes(x=Z0, y=Z1))+ geom_point() + ggtitle("Bloomfield-Amyloid_Imaging - 1989")
ggsave("Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_Amyloid_Imaging-CLEAN1-IBD.jpg", plot = last_plot(), device = NULL, scale = 1, width = 8, height = 5, dpi = 300, limitsize = TRUE)



# Merge HAPMAP ethnicity
PCA <- read.table("Bloomfield_9810-hwe-geno0.05-mind0.1-WXSm-missing-projects-include-good-IDS-V2_Amyloid_Imaging-CLEAN1_no_CHR-HAPMAP-MERGED3-for_PCA.eigenvec", header =T, stringsAsFactors=FALSE)
dim(PCA)
HAPMAP.ethnicty <- read.table("relationships_w_pops_121708.txt", header = T )
head(HAPMAP.ethnicty)


PCA$COHORT <- "Amyloid_Imaging"
PCA$COHORT <- HAPMAP.ethnicty$population[match(PCA$IID, HAPMAP.ethnicty$IID)]
PCA <- PCA[c("FID", "IID", c(paste0("PC", 1:10), "COHORT"))]
PCA$COHORT <- as.character(PCA$COHORT)
PCA$COHORT[is.na(PCA$COHORT)] <- "Amyloid_Imaging"
# write.table(PCA, "Bloomfield_Amyloid_Imaging_1989-round1.txt", sep ="\t", col.names = T, quote = F)



## Generate a new file that has IID, PC1,PC2, and a new column COHORT 
p <- ggplot(PCA, aes(x=PC1, y=PC2, color=COHORT)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Bloomfield-Amyloid-Imaging 1989") +
  scale_color_manual(values = c('green', 'black', 'red', "blue"))


## Select NHW samples only
p + annotate("rect", xmin=-0.012, xmax=0.025, ymin=-0.02, ymax= 0.055, 
             fill=NA, colour="red") +
  annotate("text", x=0.012, y=0.025, label="NHW", size=4, color = "red")

ggsave("Bloomfield-Amyloid-Imaging-PCs-COHORT.jpg", plot = last_plot(), device = NULL, scale = 1, width = 12, height = 8, dpi = 600, limitsize = TRUE)

## Select NHW only
sum(PCA$PC1 > -0.012 & PCA$PC1 < 0.025 & PCA$PC2 > -0.02 & PCA$PC2 < 0.055)

SELECTED_PC_SAMPLES <- PCA[(PCA$PC1 > -0.012 & PCA$PC1 < 0.025 & PCA$PC2 > -0.02 & PCA$PC2 < 0.055),]
table(SELECTED_PC_SAMPLES$COHORT)
# Amyloid_Imaging             CEU 
# 1882             165

## NHW samples
SELECTED_PC_SAMPLES <- SELECTED_PC_SAMPLES[grepl("Amyloid_Imaging", SELECTED_PC_SAMPLES$COHORT),]

#############################################################################
############################## GET PHENOTYPE ################################
#############################################################################
PHENO <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/01-Bloomfield-preQC/03-PLINK-QC-files/Bloomfield_8751_metadata_20211203.csv", header = T, sep = ",", stringsAsFactors = T)
dim(PHENO)
# 8751
head(PHENO)

PHENO_ALL <- PHENO

# Have pheno for all samples?
sum(WXS_FAM$Bloomfield.gvcf.id..SM...9810. %in% PHENO$Bloomfield.gvcf.id..SM...9810.)
# 1999
PHENO <- cbind(WXS_FAM[1:2], PHENO[match(WXS_FAM$Bloomfield.gvcf.id..SM...9810., PHENO$Bloomfield.gvcf.id..SM...9810.),])
PHENO <- PHENO[colnames(PHENO)[!grepl("^X", colnames(PHENO))]]



# KEEP only NHW PHENO
PHENO <- PHENO[(PHENO$IID %in% SELECTED_PC_SAMPLES$IID),]

colnames(PHENO) [grepl("STATUS", colnames(PHENO))] <- "STATUS"
colnames(PHENO) [grepl("SEX", colnames(PHENO))] <- "SEX"


## Pheno SEX
PHENO_SEX <- cbind.data.frame(FID=PHENO$FID, IID=PHENO$IID, SEX=PHENO$SEX)
## Recode SEX 
PHENO_SEX$SEX[PHENO_SEX$SEX == -9] <- 0
write.table(PHENO_SEX, "/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/06-Amyloid_Imaging-Achal/amyloid_imaging_PHENO_SEX.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)
## STATUS
STATUS <- cbind.data.frame(FID=PHENO$FID, IID=PHENO$IID, SEX=PHENO$STATUS)
write.table(STATUS, "/40/AD/AD_Seq_Data/05.-Analyses/07-Bloomfield_202109/06-Amyloid_Imaging-Achal/amyloid_imaging_PHENO_STATUS.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)

############## Carlos asked me how many WU ADNI and UPITT
# RECd_from_Muhammad
table(ALL_AMYLOID$study)
# A4_SUVR      ADNI_CL ADNIDOD_SUVR ADRC_AV45_CL     ADRC_PIB     DIAN_PIB         HABS        UPitt 
# 3244         1154          171          557          358          223          263          357
sum(grepl("MAP", AMYLOID$IID))
# 915
sum(grepl("UPitt", AMYLOID$study))
# 357
sum(grepl("ADNI_CL", AMYLOID$study))
# 1154

# ALL SEQ in JOint CALL (JOINT_VCF)
# MAP
sum(grepl("MAP", ALL_SEQ_FOUND$Bloomfield.gvcf.id..SM...9810.))
# 2228
## Site 27
sum(grepl("^27_", ALL_SEQ_FOUND$Bloomfield.gvcf.id..SM...9810.))
# 240
sum(grepl("ADNI", ALL_SEQ_FOUND$Bloomfield.gvcf.id..SM...9810.))
# 841 
sum(grepl("UPitt", ALL_SEQ_FOUND$Seq_project))
# 700



# PHENO (WXS_PHENO)
# MAP
sum(grepl("MAP", PHENO_ALL$Bloomfield.gvcf.id..SM...9810.))
# 1672
## Site 27
sum(grepl("^27_", PHENO_ALL$Bloomfield.gvcf.id..SM...9810.))
# 218
sum(grepl("ADNI", PHENO_ALL$Bloomfield.gvcf.id..SM...9810.))
# 818
sum(grepl("UPitt", PHENO_ALL$Seq_project))
# 680



## Amyloid imaging in Joint call Before QC (Amyloid WXS (overlapping, Pre-QC))
sum(grepl("MAP", WXS_FAM_B4_QC$CLEANED_ID))
# 1087
sum(grepl("^27_", WXS_FAM_B4_QC$CLEANED_ID))
# 0
sum(grepl("ADNI", WXS_FAM_B4_QC$CLEANED_ID))
# 692
sum(grepl("UPitt", WXS_FAM_B4_QC$Seq_project))
# 330


# POST QC SEQ (Amyloid WXS (overlapping; POST-IBD on WXS))
WXS_FAM
sum(grepl("^MAP", WXS_FAM$CLEANED_ID))
# 779
sum(grepl("ADNI", WXS_FAM$CLEANED_ID))
# 684
sum(grepl("UPitt", WXS_FAM$Seq_project))
# 320


tt <- WXS_FAM[!grepl("^MAP|ADNI", WXS_FAM$CLEANED_ID),]
tt <- tt [!grepl("UPitt", tt$Seq_project),]
dim(tt)


## WXS POST-QC (Amyloid WXS (overlapping; POST-PCA-IBD))
sum(grepl("^MAP", PHENO$IID))
# 694
sum(grepl("ADNI", WXS_FAM$IID))
# 684
sum(grepl("UPitt", WXS_FAM$Seq_project))
# 320







