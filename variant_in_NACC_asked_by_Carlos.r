# Email from Carlos on: Wed 11/10/2021 2:30 PM
# Subject: Early onset locus?
# Can you run a quick analyses for that SNP In NACC? - Carlos
# rs143080277 (chr2:105749599:T:C)
# https://gnomad.broadinstitute.org/variant/2-105749599-T-C?dataset=gnomad_r3

# Function to change column names
change_names <- function(x){
  colnames(x) [grepl("AGE_AT_LAST_VISIT|AGE_AT_LAST_ASSESSMENT|AGE_AT_EXAM|AGE_AT_LAST|age_last", colnames(x), ignore.case = TRUE)] <- "AGE_LAST_VISIT"
  colnames(x) [grepl("AGE_AT_DEATH", colnames(x), ignore.case = TRUE)] <- "AGE_AT_DEATH"
  colnames(x) [grepl("ONSET|AAO", colnames(x), ignore.case = TRUE)] <- "AGE_AT_ONSET"
  colnames(x) [grepl("GENDER|SEX|sex", colnames(x), ignore.case = TRUE)] <- "SEX"
  colnames(x)[grepl("final_CC_status", colnames(x), ignore.case = TRUE)] <- "Final_CC_status"
  colnames(x)[grepl("Ethnicity", colnames(x), ignore.case = TRUE)] <- "ETHNICITY"
  x}

# ADGC NHW covariate
ADGC_NHW_covar <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/NHW_covariate.tsv", header = T)
sum(grepl("NACC", ADGC_NHW_covar$IID))
# extract NACC list from covariate files
ADGC_covar_NACC <- ADGC_NHW_covar[grepl("NACC", ADGC_NHW_covar$IID),]
ADGC_covar_NACC$CLEANED_ID  <- paste0("NACC", str_extract(ADGC_covar_NACC$IID, "(?<=NACC)[0-9]*"))
ADGC_covar_NACC <- ADGC_covar_NACC[!duplicated(ADGC_covar_NACC$CLEANED_ID),]

ADGC_covar_NACC <- ADGC_covar_NACC[,c("IID", "AGE_AT_ONSET", "AGE_LAST_VISIT", "SEX", "STATUS", "ETHNICITY", "COHORT", "CLEANED_ID")]


# recode STATUS
table(ADGC_covar_NACC$STATUS)
ADGC_covar_NACC$STATUS <- as.character(ADGC_covar_NACC$STATUS)
ADGC_covar_NACC$STATUS[ADGC_covar_NACC$STATUS=="CA"] <- 2
ADGC_covar_NACC$STATUS[ADGC_covar_NACC$STATUS=="CO"] <- 1
ADGC_covar_NACC$STATUS[ADGC_covar_NACC$STATUS=="MCI"] <- 3
ADGC_covar_NACC$STATUS[ADGC_covar_NACC$STATUS=="unknown"] <- -9
ADGC_covar_NACC_ALL <- ADGC_covar_NACC
# Select CA and CO by age; CA <= 65; CO > 70
ADGC_covar_NACC <- ADGC_covar_NACC[which((ADGC_covar_NACC$STATUS == "2" & ADGC_covar_NACC$AGE_AT_ONSET <= 65)| (ADGC_covar_NACC$STATUS == "1" & ADGC_covar_NACC$AGE_LAST_VISIT > 70)),]
sum(grepl("NACC", ADGC_covar_NACC$IID))



# NACC PHENO with selected STATUS from Fengxian
NACC <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/WASHU-GWAS/NACC_unique_core_pheno_20210927.csv", sep = ",", header = T)
# Selecting only the AD samples
NACC$STATUS <- as.character(NACC$final_CC_status)
NACC <- NACC[grepl("^CA$|^CO$|^Neuro_AD$|^Neuro_AD_DLB$|^Neuro_CO$|^Neuro_PreSymptomatic_AD$|OT\\(CO\\)", NACC$STATUS),]
table(NACC$STATUS)

# Recode CA/CO
NACC$STATUS[grepl("^CO$|^Neuro_CO$|OT\\(CO\\)", NACC$STATUS)] <- 1
NACC$STATUS[grepl("^CA$|^Neuro_AD$|^Neuro_AD_DLB$|^Neuro_PreSymptomatic_AD$", NACC$STATUS)] <- 2
table(NACC$STATUS)

NACC <- change_names(NACC)
NACC$SEX <- as.character(NACC$SEX)
# Recode SEX
NACC$SEX[NACC$SEX == "Female"] <- 2
NACC$SEX[NACC$SEX == "Male"] <- 1

NACC_ALL <- NACC
# Select CA and CO by age; CA <= 65; CO > 70
NACC <- NACC[which((NACC$STATUS == "2" & NACC$AGE_AT_ONSET <= 65)| (NACC$STATUS == "1" & NACC$AGE_LAST_VISIT > 70)),]
# Get clean NACC IDs
NACC$IID <- NACC$NACC_ID
NACC$CLEANED_ID  <- paste0("NACC", str_extract(NACC$NACC_ID, "(?<=NACC)[0-9]*"))
NACC <- NACC[!duplicated(NACC$CLEANED_ID),]


# Find NACC from Fengxian that are not included in ADGC PHENO
NACC <- NACC[!NACC$CLEANED_ID %in% ADGC_covar_NACC$CLEANED_ID,]

# COHORT (STUDY) is unknown in FENGXIAN's PHENO  
NACC$COHORT <- NA

NACC <- NACC[,c("IID", "AGE_AT_ONSET", "AGE_LAST_VISIT", "SEX", "STATUS", "ETHNICITY", "COHORT", "CLEANED_ID")]
NACC <- rbind.data.frame(ADGC_covar_NACC, NACC)

# How many of these are in AGDC_NHW .FAM?
ADGC_NHW_FAM <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ADGC_NHW_Cohort.fam", header = F)
NACC_ADGC_NHW_FAM <- NACC_ADGC_NHW_FAM[grepl("NACC", ADGC_NHW_FAM$V2),]
NACC_ADGC_NHW_FAM$CLEANED_ID <- paste0("NACC", str_extract(NACC_ADGC_NHW_FAM$V2, "(?<=NACC)[0-9]*"))
NACC_ADGC_NHW_FAM_ALL <- NACC_ADGC_NHW_FAM
sum(NACC$CLEANED_ID %in% NACC_ADGC_NHW_FAM$CLEANED_ID)
# 7169
NACC <- NACC[NACC$CLEANED_ID %in% NACC_ADGC_NHW_FAM$CLEANED_ID,]
# Get FID from FAM file in Pheno file
NACC$IID <- NACC$CLEANED_ID
NACC$FID <- as.character(NACC_ADGC_NHW_FAM$V1[match(NACC$IID, NACC_ADGC_NHW_FAM$V2)])
# relocate FID to the first position
NACC <- NACC %>%
  relocate(FID)

NACC$COHORT <- as.character(NACC$COHORT)

# Find STUDY information for samples from Fengxian in ADGC
sum(is.na(NACC$COHORT))
# 2050 samples

# All FAM files
fileName <- "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/01-EOAD-preQC/NACC_analysis/myfinal.list"

## Loop over a file connection
conn <- file(fileName,open="r")
linn <-readLines(conn)
df1 <- {}
for (i in 1:length(linn)){
print(linn[i])
COHORT <- str_match(linn[i], "ADGC_NHW-\\s*(.*?)\\s*_geno") [,2] 
df.tmp <- read.table(linn[i], header = F, colClasses='character')  
df.tmp$COHORT <- COHORT
df1 <- rbind.data.frame(df1,df.tmp)
}

close(conn)

df1 <- df1[grepl("NACC", df1$V2),]
df1$CLEANED_ID <- paste0("NACC", str_extract(df1$V2, "(?<=NACC)[0-9]*"))
df1 <- df1[!duplicated(df1$CLEANED_ID),]

# Get the STUDY information for the samples
sum(NACC$CLEANED_ID %in% df1$CLEANED_ID)
# 7169
NACC$STUDY <- df1$COHORT[match(NACC$CLEANED_ID, df1$CLEANED_ID)]

table(NACC$STATUS)
# 1 (CO)    2 (CA) 
# 4963 2206 

# Model the phenotype file
table(NACC$STUDY)
# ADC1 ADC10 ADC11 ADC12  ADC2  ADC3  ADC4  ADC5  ADC6  ADC7  ADC8  ADC9 
# 675  1096   558   386   259   433   475   644   604   804   462   773
## ADC10 has the max num of samples

NACC$STUDY <- as.factor(NACC$STUDY)
CLEAN <- function(x){ colnames(x) <- gsub("STUDY", "", colnames(x)); x } 
NACC <- cbind(NACC, CLEAN(as.data.frame(model.matrix(~STUDY, data=NACC))))

# merge PCs
PCA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/01-EOAD-preQC/NACC_analysis/ADGC_NHW_Cohort-PCAS.eigenvec", header = T)
NACC <- cbind.data.frame(NACC,PCA[match(NACC$IID, PCA$IID),3:12])


write.table(NACC, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/01-EOAD-preQC/phenofiles/NHW_NACC_final_phenotype_11_16_2021.tsv", sep ="\t", col.names = T, quote = F, row.names = FALSE)
write.table(NACC[1:2], "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/01-EOAD-preQC/phenofiles/NACC_NHW.IDlist", sep ="\t", col.names = FALSE, quote = F, row.names = FALSE)
# After analysis
df <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/01-EOAD-preQC/NACC_analysis/ADGC_NHW_Cohort_rs143080277.ped", stringsAsFactors = FALSE, header = T, colClasses = c("character"))
# Label headers 
cols <- paste0(rep("rs143080277", each = 2), c("_A1", "_A2"))
cols <- c("FID", "IID", "PAT", "MAT", "SEX", "STATUS", cols)
colnames(df) <- cols

rownames(df) <- df$IID
df <- df[cols[-c(1:6)]]

# ## If A1 == A2 ---> Homozygous; If A1 != A2 ---> Heterozygous (carrier);
# odds <- seq(1, ncol(df), by = 2)
# odds
# df <- do.call(cbind, Map(function(z, nm) setNames(cbind(z, Status = +(z[[1]] != z[[2]])), c(names(z), nm)),
#                          split.default(df, rep(odds, each = 2)),
#                          paste0(gsub("_A1","", colnames(df)[odds]), "_Status"), USE.NAMES = FALSE))


# Carriers: If A1 == A2 & A1 ==1 ---> Homozygous_Nonreference; If A1 != A2 ---> Heterozygous
homozygous_nonreference <- df[(df$rs143080277_A1 == df$rs143080277_A2) & df$rs143080277_A1 ==1,]
heterozygous <- df[(df$rs143080277_A1 != df$rs143080277_A2),]
homozygous_reference <- df[(df$rs143080277_A1 == df$rs143080277_A2) & df$rs143080277_A1 ==2,]

HOMO_NONREF <- nrow(homozygous_nonreference)
HETERO <- nrow(heterozygous)
HOMO_REF <- nrow(homozygous_reference)

# # Get Heterozygosity Status
# df <- df[grepl("Status", colnames(df))]
# ALL_CARRIERS <- sum(df ==1)
# 
# 
# # get Carrier count
# df.TF <- df ==1
# if (sum(df.TF) == 0){
#   next
# }

# cc <- colSums(df.TF)[colSums(df.TF) > 0]
# variants <- as.data.frame(as.vector(cc))

# # Get the list of carriers
# Carrier.samples.tmp <- row.names(df.TF)[df.TF[,names(cc[n])]]
# Carrier.samples <- paste(names(cc[n]), paste(Carrier.samples.tmp, sep="", collapse=","), sep = "\t" )

Carrier.samples <- c(rownames(homozygous_nonreference), rownames(heterozygous))
NCarriers <- length(Carrier.samples)

# Age range of carriers
# AAL (age range 55-100)
tt <- NACC[match(Carrier.samples, NACC$CLEANED_ID),]

table(tt$STATUS)
# 1  2 
# 58 35 

# AAL (age range 71-100, 58 of 93 were controls)
aggregate(tt$AGE_LAST_VISIT, by=list(tt$STATUS), range)
# AAO (age range 52-65, 35 of 93 were cases)
aggregate(tt$AGE_AT_ONSET, by=list(tt$STATUS), range)

# Total samples 7169 of which 93 were carriers (90 heterozygous, 3 homozygous nonreference). Out of 93 carriers, 35 are cases with AAO range 52-65. AAL range for 58 CO carriers was 71-100. 

###########################
## Analysis with CO > 80 ##
###########################

# Plot PCA
library(ggplot2)
ggplot(NACC, aes(x=PC1, y=PC2)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("ADGC_NACC")
ggsave("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/01-EOAD-preQC/NACC_analysis/1769_NACC_PCA_without_HAPMAP.jpg", plot = last_plot(), device = NULL, scale = 1, width = 16, height = 9, dpi = 300, limitsize = TRUE)




##############################
## Analysis II with CO > 80 ##
##############################

# Age threshold for CO > 80 
NACC_CO80 <- NACC[which((NACC$STATUS == "2" & NACC$AGE_AT_ONSET <= 65)| (NACC$STATUS == "1" & NACC$AGE_LAST_VISIT > 80)),]
dim(NACC_CO80)
# [1] 3981   32

table(NACC_CO80$STATUS)
# 1    2 
# 1775 2206 

write.table(NACC_CO80, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/01-EOAD-preQC/phenofiles/NHW_NACC_CO80_final_phenotype_11_16_2021.tsv", sep ="\t", col.names = T, quote = F, row.names = FALSE)
write.table(NACC_CO80[1:2], "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/01-EOAD-preQC/phenofiles/NACC_NHW_CO80.IDlist", sep ="\t", col.names = FALSE, quote = F, row.names = FALSE)

####################

# Carrier counts
df <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/01-EOAD-preQC/NACC_analysis/ADGC_NHW_NACC_CO80_Cohort_rs143080277.ped", stringsAsFactors = FALSE, header = T, colClasses = c("character"))
# Label headers 
cols <- paste0(rep("rs143080277", each = 2), c("_A1", "_A2"))
cols <- c("FID", "IID", "PAT", "MAT", "SEX", "STATUS", cols)
colnames(df) <- cols

rownames(df) <- df$IID
df <- df[cols[-c(1:6)]]

# Carriers: If A1 == A2 & A1 ==1 ---> Homozygous_Nonreference; If A1 != A2 ---> Heterozygous
homozygous_nonreference <- df[(df$rs143080277_A1 == df$rs143080277_A2) & df$rs143080277_A1 ==1,]
heterozygous <- df[(df$rs143080277_A1 != df$rs143080277_A2),]
homozygous_reference <- df[(df$rs143080277_A1 == df$rs143080277_A2) & df$rs143080277_A1 ==2,]

HOMO_NONREF <- nrow(homozygous_nonreference)
# 2
HETERO <- nrow(heterozygous)
# 50
HOMO_REF <- nrow(homozygous_reference)
# 3928


Carrier.samples <- c(rownames(homozygous_nonreference), rownames(heterozygous))
NCarriers <- length(Carrier.samples)
# 52

# Age range of carriers
tt <- NACC[match(Carrier.samples, NACC$CLEANED_ID),]

table(tt$STATUS)
# 1  2 
# 17 35 

# AAL (age range 81-100, 17 of 52 were controls)
aggregate(tt$AGE_LAST_VISIT, by=list(tt$STATUS), range)
# AAO (age range 52-65, 35 of 52 were cases)
aggregate(tt$AGE_AT_ONSET, by=list(tt$STATUS), range)

# Total samples 3981 of which 52 were carriers (50 heterozygous, 2 homozygous nonreference). Out of 52 carriers, 35 are cases with AAO range 52-65. AAL range for 17 CO carriers was 81-100. 


###########################################################################################
# ALL NACC in ADGC NHW covariate file
dim(ADGC_covar_NACC_ALL)
# 17012
# ALL NACC in ADGC covariate file and are also in ADGC FAM file
sum(ADGC_covar_NACC_ALL$CLEANED_ID %in% NACC_ADGC_NHW_FAM_ALL$CLEANED_ID)
# 17012
ADGC_covar_NACC_ALL <- ADGC_covar_NACC_ALL[ADGC_covar_NACC_ALL$CLEANED_ID %in% NACC_ADGC_NHW_FAM_ALL$CLEANED_ID,]

# ALL AD NACC from Fengxian
dim(NACC_ALL)
# 28319
NACC_ALL$CLEANED_ID  <- paste0("NACC", str_extract(NACC_ALL$NACC_ID, "(?<=NACC)[0-9]*"))

# ALL AD NACC from Fengxian and are also in ADGC NHW FAM file
sum(NACC_ALL$CLEANED_ID %in% NACC_ADGC_NHW_FAM_ALL$CLEANED_ID)
# 11905
NACC_ALL <- NACC_ALL[NACC_ALL$CLEANED_ID %in% NACC_ADGC_NHW_FAM_ALL$CLEANED_ID,]


## Add STUDY info to Fengxian's NACC samples
NACC_ALL$COHORT <- df1$COHORT[match(NACC_ALL$CLEANED_ID, df1$CLEANED_ID)]
table(NACC$COHORT)
# NACC_ALL_CACO <- NACC_ALL[which((NACC_ALL$STATUS == "2" & NACC_ALL$AGE_AT_ONSET <= 65)| (NACC_ALL$STATUS == "1" & NACC_ALL$AGE_LAST_VISIT > 70)),]
# ADGC_covar_NACC_ALL_CACO <- ADGC_covar_NACC_ALL[which((ADGC_covar_NACC_ALL$STATUS == "2" & ADGC_covar_NACC_ALL$AGE_AT_ONSET <= 65)| (ADGC_covar_NACC_ALL$STATUS == "1" & ADGC_covar_NACC_ALL$AGE_LAST_VISIT > 70)),]
# NACC_NHW_ALL_PHENO_CACO <- rbind.data.frame(NACC_ALL_CACO[,c("CLEANED_ID", "AGE_AT_ONSET", "AGE_LAST_VISIT", "SEX", "STATUS", "ETHNICITY", "COHORT")], ADGC_covar_NACC_ALL_CACO[,c("CLEANED_ID", "AGE_AT_ONSET", "AGE_LAST_VISIT", "SEX", "STATUS", "ETHNICITY", "COHORT")])
# NACC_NHW_ALL_PHENO_CACO <- NACC_NHW_ALL_PHENO_CACO[!duplicated(NACC_NHW_ALL_PHENO_CACO$CLEANED_ID),]
# dim(NACC_NHW_ALL_PHENO_CACO)
# # 7169

NACC_NHW_ALL_PHENO <- rbind.data.frame(NACC_ALL[,c("CLEANED_ID", "AGE_AT_ONSET", "AGE_LAST_VISIT", "SEX", "STATUS", "ETHNICITY", "COHORT")], ADGC_covar_NACC_ALL[,c("CLEANED_ID", "AGE_AT_ONSET", "AGE_LAST_VISIT", "SEX", "STATUS", "ETHNICITY", "COHORT")])
NACC_NHW_ALL_PHENO <- NACC_NHW_ALL_PHENO[!duplicated(NACC_NHW_ALL_PHENO$CLEANED_ID),]

dim(NACC_NHW_ALL_PHENO)
# 17947

table(NACC_NHW_ALL_PHENO$COHORT)
# ADC1 ADC10 ADC11 ADC12  ADC2  ADC4  ADC5  ADC6  ADC7  ADC8  ADC9 
# 2208  1295   455   366   271   181   271   429   213   159   194 

##########################
##########################
##### With all NHW #######
##########################
##########################
##################################################### Count carriers
## CO > 70

df <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/01-EOAD-preQC/NHW_analysis/analysis_for_early_onset_locus_rs143080277/ADGC_NHW_Cohort_rs143080277_CO70.ped", stringsAsFactors = FALSE, header = T, colClasses = c("character"))
# Label headers 
cols <- paste0(rep("rs143080277", each = 2), c("_A1", "_A2"))
cols <- c("FID", "IID", "PAT", "MAT", "SEX", "STATUS", cols)
colnames(df) <- cols

df$KEY <- paste(df$FID, df$IID, sep = ":")
rownames(df) <- df$KEY
df <- df[cols[-c(1:6, 9)]]

# ## If A1 == A2 ---> Homozygous; If A1 != A2 ---> Heterozygous (carrier);
# odds <- seq(1, ncol(df), by = 2)
# odds
# df <- do.call(cbind, Map(function(z, nm) setNames(cbind(z, Status = +(z[[1]] != z[[2]])), c(names(z), nm)),
#                          split.default(df, rep(odds, each = 2)),
#                          paste0(gsub("_A1","", colnames(df)[odds]), "_Status"), USE.NAMES = FALSE))


# Carriers: If A1 == A2 & A1 ==1 ---> Homozygous_Nonreference; If A1 != A2 ---> Heterozygous
homozygous_nonreference <- df[(df$rs143080277_A1 == df$rs143080277_A2) & df$rs143080277_A1 ==1,]
heterozygous <- df[(df$rs143080277_A1 != df$rs143080277_A2),]
homozygous_reference <- df[(df$rs143080277_A1 == df$rs143080277_A2) & df$rs143080277_A1 ==2,]

HOMO_NONREF <- nrow(homozygous_nonreference)
HETERO <- nrow(heterozygous)
HOMO_REF <- nrow(homozygous_reference)

# # Get Heterozygosity Status
# df <- df[grepl("Status", colnames(df))]
# ALL_CARRIERS <- sum(df ==1)
# 
# 
# # get Carrier count
# df.TF <- df ==1
# if (sum(df.TF) == 0){
#   next
# }

# cc <- colSums(df.TF)[colSums(df.TF) > 0]
# variants <- as.data.frame(as.vector(cc))

# # Get the list of carriers
# Carrier.samples.tmp <- row.names(df.TF)[df.TF[,names(cc[n])]]
# Carrier.samples <- paste(names(cc[n]), paste(Carrier.samples.tmp, sep="", collapse=","), sep = "\t" )

Carrier.samples <- c(rownames(homozygous_nonreference), rownames(heterozygous))
NCarriers <- length(Carrier.samples)
# 181

# Age range of carriers
tt <- combined_NHW_CO70[match(Carrier.samples, combined_NHW_CO70$KEY2),]
dim(tt)
# 181 77
table(tt$STATUS)
# 1   2 
# 135  46 

# AAL (age range 71-100, 135 of 181 were controls)
aggregate(tt$AGE_LAST_VISIT, by=list(tt$STATUS), range)
# AAO (age range 52-65, 46 of 181 were cases)
aggregate(tt$AGE_AT_ONSET, by=list(tt$STATUS), range)
###################################################################################

## CO > 80

df <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/01-EOAD-preQC/NHW_analysis/analysis_for_early_onset_locus_rs143080277/ADGC_NHW_Cohort_rs143080277_CO80.ped", stringsAsFactors = FALSE, header = T, colClasses = c("character"))
# Label headers 
cols <- paste0(rep("rs143080277", each = 2), c("_A1", "_A2"))
cols <- c("FID", "IID", "PAT", "MAT", "SEX", "STATUS", cols)
colnames(df) <- cols

df$KEY <- paste(df$FID, df$IID, sep = ":")
rownames(df) <- df$KEY
df <- df[cols[-c(1:6, 9)]]

# ## If A1 == A2 ---> Homozygous; If A1 != A2 ---> Heterozygous (carrier);
# odds <- seq(1, ncol(df), by = 2)
# odds
# df <- do.call(cbind, Map(function(z, nm) setNames(cbind(z, Status = +(z[[1]] != z[[2]])), c(names(z), nm)),
#                          split.default(df, rep(odds, each = 2)),
#                          paste0(gsub("_A1","", colnames(df)[odds]), "_Status"), USE.NAMES = FALSE))


# Carriers: If A1 == A2 & A1 ==1 ---> Homozygous_Nonreference; If A1 != A2 ---> Heterozygous
homozygous_nonreference <- df[(df$rs143080277_A1 == df$rs143080277_A2) & df$rs143080277_A1 ==1,]
heterozygous <- df[(df$rs143080277_A1 != df$rs143080277_A2),]
homozygous_reference <- df[(df$rs143080277_A1 == df$rs143080277_A2) & df$rs143080277_A1 ==2,]

HOMO_NONREF <- nrow(homozygous_nonreference)
HETERO <- nrow(heterozygous)
HOMO_REF <- nrow(homozygous_reference)

# # Get Heterozygosity Status
# df <- df[grepl("Status", colnames(df))]
# ALL_CARRIERS <- sum(df ==1)
# 
# 
# # get Carrier count
# df.TF <- df ==1
# if (sum(df.TF) == 0){
#   next
# }

# cc <- colSums(df.TF)[colSums(df.TF) > 0]
# variants <- as.data.frame(as.vector(cc))

# # Get the list of carriers
# Carrier.samples.tmp <- row.names(df.TF)[df.TF[,names(cc[n])]]
# Carrier.samples <- paste(names(cc[n]), paste(Carrier.samples.tmp, sep="", collapse=","), sep = "\t" )

Carrier.samples <- c(rownames(homozygous_nonreference), rownames(heterozygous))
NCarriers <- length(Carrier.samples)
# 94

# Age range of carriers
tt <- combined_NHW_CO70[match(Carrier.samples, combined_NHW_CO70$KEY2),]
dim(tt)
# 94 77
table(tt$STATUS)
# 1  2 
# 48 46

# AAL (age range 81-100, 48 of 94 were controls)
aggregate(tt$AGE_LAST_VISIT, by=list(tt$STATUS), range)
# AAO (age range 52-65, 46 of 94 were cases)
aggregate(tt$AGE_AT_ONSET, by=list(tt$STATUS), range)






# Now create a column for match and mismatch 
ADGC_NHW_covar$CLEANED_ID <- as.character(ADGC_NHW_covar$IID)
ADGC_NHW_covar$FOUND <- ifelse(ADGC_NHW_covar$CLEANED_ID %in% ADGC_NHW_FAM$V2, "YES", "NO")
table(ADGC_NHW_covar$FOUND)
# NO   YES 
# 20971 29273

# First clean NACC IDs in NHW covar
ADGC_NHW_covar$CLEANED_ID[grepl("NACC", ADGC_NHW_covar$IID)] <- paste0("NACC", str_extract(ADGC_NHW_covar$IID, "(?<=NACC)[0-9]*")) [grepl("NACC", ADGC_NHW_covar$IID)]

# CHECK MATCHES !!
ADGC_NHW_covar$FOUND <- ifelse(ADGC_NHW_covar$CLEANED_ID %in% ADGC_NHW_FAM$V2, "YES", "NO")
table(ADGC_NHW_covar$FOUND)
