setwd("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101")
wantID <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/FASe_ID_list_round1.txt", header = FALSE, sep =",", check.names = FALSE)
head(wantID)
colnames(wantID) <- c("Original_ID", "VCF_ID", "PR")
dim(wantID)

wantID[] <- lapply(wantID, function(x) as.vector(as.character(x)))



wantID[grepl("2_18_4",wantID$Original_ID),]


LOOKUP <- read.csv("/40/AD/AD_Seq_Data/03.-phenotype/2019-04-Aurora-phenotype/20191210-Aurora_phenotye_UPDATED.csv", header = TRUE)
head(LOOKUP)
dim(LOOKUP)

LOOKUP[] <- lapply(LOOKUP, function(x) as.vector(as.character(x)))

LOOKUP[grepl("0_62802_45", LOOKUP$IID),]
LOOKUP[grepl("10R_R54_10", LOOKUP$IID),]


wanted.final <- wantID
wanted.final$IID <- "EMPTY"
wanted.final$FID <- "EMPTY"
wanted.final$Sex <- "EMPTY"
wanted.final$STATUS <- "EMPTY"
wanted.final$SEXCHECK <- "EMPTY"
wanted.final$phenoSCOPE <- "EMPTY"
wanted.final$phenoProject <- "EMPTY"

# i=274
# i =247
warning()
for (i in 1:length(wantID$Original_ID)){
print(paste0("DOING: ", i))  
varTosearch <- paste0(wantID$Original_ID[i],"$", collapse = "$")
# print(grep(varTosearch, LOOKUP$IID))
# tmpLOOKUP <- c(wantID[i,],LOOKUP[grep(varTosearch, LOOKUP$IID), c("IID", "FID", "Sex", "STATUS", "F.SEXCHECK.", "SCOPE", "Project")])
tmpLOOKUP <- c(LOOKUP[grep(varTosearch, LOOKUP$IID), c("IID", "FID", "Sex", "STATUS", "F.SEXCHECK.", "SCOPE", "Project")])
tmpLOOKUP

## Remove suplicates from different projects

if(is.na(tmpLOOKUP$SCOPE) || length(tmpLOOKUP$SCOPE)==0){
  tmpLOOKUP$SCOPE <- "NONE" 
}

if(is.na(tmpLOOKUP$Project) || length(tmpLOOKUP$Project)==0){
  tmpLOOKUP$Project <- "NONE" 
}

indextograb <- c(grep(wanted.final[i,"PR"], tmpLOOKUP$Project),grep(wanted.final[i,"PR"], tmpLOOKUP$SCOPE))
if (length(indextograb) == 0){ 
print("BAD Integer!!")
next  
}

# tmpLOOKUP$IID <- tmpLOOKUP$IID[grepl(varTosearch, tmpLOOKUP$IID)][indextograb]
tmpLOOKUP$IID <- tmpLOOKUP$IID[indextograb]
tmpLOOKUP$FID <- tmpLOOKUP$FID[indextograb]
tmpLOOKUP$Sex <- tmpLOOKUP$Sex[indextograb]
tmpLOOKUP$STATUS <- tmpLOOKUP$STATUS[indextograb]
tmpLOOKUP$F.SEXCHECK. <- tmpLOOKUP$F.SEXCHECK.[indextograb]


tmpLOOKUP$SCOPE <- tmpLOOKUP$SCOPE[grepl(wanted.final[i,"PR"], tmpLOOKUP$SCOPE)]
tmpLOOKUP$Project <- tmpLOOKUP$Project[grepl(wanted.final[i,"PR"], tmpLOOKUP$Project)]



if((wanted.final[i,"PR"]==tmpLOOKUP$SCOPE || wanted.final[i,"PR"]==tmpLOOKUP$Project) && wanted.final[i,"IID"]=="EMPTY"){
# if(wanted.final[i,"PR"]==tmpLOOKUP$SCOPE || wanted.final[i,"PR"]==tmpLOOKUP$PR){
if(length(tmpLOOKUP$IID)>0){
wanted.final$IID[i] <- tmpLOOKUP$IID
}
if(length(tmpLOOKUP$FID)>0){
wanted.final$FID[i] <- tmpLOOKUP$FID
}
if(length(tmpLOOKUP$Sex)>0){
wanted.final$Sex[i] <- tmpLOOKUP$Sex
}
if(length(tmpLOOKUP$STATUS)>0){
wanted.final$STATUS[i] <- tmpLOOKUP$STATUS
}
if(length(tmpLOOKUP$F.SEXCHECK.)>0){
wanted.final$SEXCHECK[i] <- tmpLOOKUP$F.SEXCHECK.
}
if(length(tmpLOOKUP$SCOPE)>0){
wanted.final$phenoSCOPE[i] <- tmpLOOKUP$SCOPE
}
if(length(tmpLOOKUP$Project)>0){
wanted.final$phenoProject[i] <- tmpLOOKUP$Project
}
}
}

# Counts of bad ones
sum(grepl("EMPTY",wanted.final$IID))
GOOD1 <-  wanted.final[!grepl("EMPTY",wanted.final$IID),]
dim(GOOD1)
BAD1 <- wanted.final[grepl("EMPTY",wanted.final$IID),]
dim(BAD1)


# different projects in BAD1
unique(BAD1$PR)
# [1] "Macrogen_WGS"              "LOAD_WES"                  "Genentech_WGS"             "Genentech_WES"             "Broad_WGS"                
# [6] "TGI_WES"                   "MGI_FASeEOAD_201605"       "Otogenetics_WES"           "201907_USUHS_gDNA_SHERMAN" "phs000572_201508"         
# [11] "phs000572_201612"

# MacrogenWGS
Macrogen_WGS_BAD1 <- BAD1[grepl("Macrogen_WGS", BAD1$PR),]
# first find exact match
Macrogen_WGS_GOOD2 <- cbind(Macrogen_WGS_BAD1[c("Original_ID", "VCF_ID", "PR")], LOOKUP[which(LOOKUP$IID %in% Macrogen_WGS_BAD1$Original_ID),c("IID", "FID", "Sex", "STATUS", "F.SEXCHECK.", "SCOPE", "Project")])
Macrogen_WGS_BAD2 <- Macrogen_WGS_BAD1[which(! Macrogen_WGS_BAD1$Original_ID %in% Macrogen_WGS_GOOD2$Original_ID),]


# LOAD_WES
LOAD_WES_BAD1 <- BAD1[grepl("LOAD_WES", BAD1$PR),]
LOAD_WES_LOOKUP <- LOOKUP[which(LOOKUP$IID %in% LOAD_WES_BAD1$Original_ID),c("IID", "FID", "Sex", "STATUS", "F.SEXCHECK.", "SCOPE", "Project")]
LOAD_WES_LOOKUP <- LOAD_WES_LOOKUP[grepl("WES", LOAD_WES_LOOKUP$SCOPE),]
indexed <- match(LOAD_WES_LOOKUP$IID, LOAD_WES_BAD1$Original_ID)
LOAD_WES_GOOD2 <- cbind(LOAD_WES_BAD1[indexed, c("Original_ID", "VCF_ID", "PR")], LOAD_WES_LOOKUP)
LOAD_WES_BAD2 <- LOAD_WES_BAD1[which(! LOAD_WES_BAD1$Original_ID %in% LOAD_WES_GOOD2$Original_ID),]

# Genentech_WGS
Genentech_WGS_BAD1 <- BAD1[grepl("Genentech_WGS", BAD1$PR),]
Genentech_WGS_LOOKUP <- LOOKUP[which(LOOKUP$IID %in% Genentech_WGS_BAD1$Original_ID),c("IID", "FID", "Sex", "STATUS", "F.SEXCHECK.", "SCOPE", "Project")]
Genentech_WGS_LOOKUP <- Genentech_WGS_LOOKUP[grepl("WGS", Genentech_WGS_LOOKUP$SCOPE),]
indexed <- match(Genentech_WGS_LOOKUP$IID, Genentech_WGS_BAD1$Original_ID)
Genentech_WGS_GOOD2 <- cbind(Genentech_WGS_BAD1[indexed, c("Original_ID", "VCF_ID", "PR")], Genentech_WGS_LOOKUP)
Genentech_WGS_BAD2 <- Genentech_WGS_BAD1[which(! Genentech_WGS_BAD1$Original_ID %in% Genentech_WGS_GOOD2$Original_ID),]

# Genentech_WES
Genentech_WES_BAD1 <- BAD1[grepl("Genentech_WES", BAD1$PR),]
Genentech_WES_LOOKUP <- LOOKUP[which(LOOKUP$IID %in% Genentech_WES_BAD1$Original_ID),c("IID", "FID", "Sex", "STATUS", "F.SEXCHECK.", "SCOPE", "Project")]
Genentech_WES_LOOKUP <- Genentech_WES_LOOKUP[grepl("WES", Genentech_WES_LOOKUP$SCOPE),]
# first find exact match
indexed <- match(Genentech_WES_LOOKUP$IID, Genentech_WES_BAD1$Original_ID)
Genentech_WES_GOOD2 <- cbind(Genentech_WES_BAD1[indexed, c("Original_ID", "VCF_ID", "PR")], Genentech_WES_LOOKUP)
Genentech_WES_BAD2 <- Genentech_WES_BAD1[which(! Genentech_WES_BAD1$Original_ID %in% Genentech_WES_GOOD2$Original_ID),]

# Broad_WGS
Broad_WGS_BAD1 <- BAD1[grepl("Broad_WGS", BAD1$PR),]
Broad_WGS_LOOKUP <- LOOKUP[which(LOOKUP$IID %in% Broad_WGS_BAD1$Original_ID),c("IID", "FID", "Sex", "STATUS", "F.SEXCHECK.", "SCOPE", "Project")]
Broad_WGS_LOOKUP <- Broad_WGS_LOOKUP[grepl("WGS", Broad_WGS_LOOKUP$SCOPE),]
# first find exact match
indexed <- match(Broad_WGS_LOOKUP$IID, Broad_WGS_BAD1$Original_ID)
Broad_WGS_GOOD2 <- cbind(Broad_WGS_BAD1[indexed, c("Original_ID", "VCF_ID", "PR")], Broad_WGS_LOOKUP)
Broad_WGS_BAD2 <- Broad_WGS_BAD1[which(! Broad_WGS_BAD1$Original_ID %in% Broad_WGS_GOOD2$Original_ID),]

# TGI_WES
TGI_WES_BAD1 <- BAD1[grepl("TGI_WES", BAD1$PR),]
TGI_WES_LOOKUP <- LOOKUP[which(LOOKUP$IID %in% TGI_WES_BAD1$Original_ID),c("IID", "FID", "Sex", "STATUS", "F.SEXCHECK.", "SCOPE", "Project")]
TGI_WES_LOOKUP <- TGI_WES_LOOKUP[grepl("WES", TGI_WES_LOOKUP$SCOPE),]
# first find exact match
indexed <- match(TGI_WES_LOOKUP$IID, TGI_WES_BAD1$Original_ID)
TGI_WES_GOOD2 <- cbind(TGI_WES_BAD1[indexed, c("Original_ID", "VCF_ID", "PR")], TGI_WES_LOOKUP)
TGI_WES_BAD2 <- TGI_WES_BAD1[which(! TGI_WES_BAD1$Original_ID %in% TGI_WES_GOOD2$Original_ID),]

# MGI_FASeEOAD_201605
MGI_FASeEOAD_201605_BAD1 <- BAD1[grepl("MGI_FASeEOAD_201605", BAD1$PR),]
MGI_FASeEOAD_201605_LOOKUP <- LOOKUP[which(LOOKUP$IID %in% MGI_FASeEOAD_201605_BAD1$Original_ID),c("IID", "FID", "Sex", "STATUS", "F.SEXCHECK.", "SCOPE", "Project")]
MGI_FASeEOAD_201605_LOOKUP <- MGI_FASeEOAD_201605_LOOKUP[grepl("WES", MGI_FASeEOAD_201605_LOOKUP$SCOPE),]
# first find exact match
indexed <- match(MGI_FASeEOAD_201605_LOOKUP$IID, MGI_FASeEOAD_201605_BAD1$Original_ID)
MGI_FASeEOAD_201605_GOOD2 <- cbind(MGI_FASeEOAD_201605_BAD1[indexed, c("Original_ID", "VCF_ID", "PR")], MGI_FASeEOAD_201605_LOOKUP)
MGI_FASeEOAD_201605_BAD2 <- MGI_FASeEOAD_201605_BAD1[which(! MGI_FASeEOAD_201605_BAD1$Original_ID %in% MGI_FASeEOAD_201605_GOOD2$Original_ID),]

# Otogenetics_WES
Otogenetics_WES_BAD1 <- BAD1[grepl("Otogenetics_WES", BAD1$PR),]
Otogenetics_WES_LOOKUP <- LOOKUP[which(LOOKUP$IID %in% Otogenetics_WES_BAD1$Original_ID),c("IID", "FID", "Sex", "STATUS", "F.SEXCHECK.", "SCOPE", "Project")]
Otogenetics_WES_LOOKUP <- Otogenetics_WES_LOOKUP[grepl("WES", Otogenetics_WES_LOOKUP$SCOPE),]
# first find exact match
indexed <- match(Otogenetics_WES_LOOKUP$IID, Otogenetics_WES_BAD1$Original_ID)
Otogenetics_WES_GOOD2 <- cbind(Otogenetics_WES_BAD1[indexed, c("Original_ID", "VCF_ID", "PR")], Otogenetics_WES_LOOKUP)
Otogenetics_WES_BAD2 <- Otogenetics_WES_BAD1[which(! Otogenetics_WES_BAD1$Original_ID %in% Otogenetics_WES_GOOD2$Original_ID),]

# phs000572_201508
phs000572_201508_BAD1 <- BAD1[grepl("phs000572_201508", BAD1$PR),]
phs000572_201508_LOOKUP <- LOOKUP[which(LOOKUP$IID %in% phs000572_201508_BAD1$Original_ID),c("IID", "FID", "Sex", "STATUS", "F.SEXCHECK.", "SCOPE", "Project")]
phs000572_201508_LOOKUP <- phs000572_201508_LOOKUP[grepl("WES", phs000572_201508_LOOKUP$SCOPE),]
# first find exact match
indexed <- match(phs000572_201508_LOOKUP$IID, phs000572_201508_BAD1$Original_ID)
phs000572_201508_GOOD2 <- cbind(phs000572_201508_BAD1[indexed, c("Original_ID", "VCF_ID", "PR")], phs000572_201508_LOOKUP)
phs000572_201508_BAD2 <- phs000572_201508_BAD1[which(! phs000572_201508_BAD1$Original_ID %in% phs000572_201508_GOOD2$Original_ID),]


# phs000572_201612
phs000572_201612_BAD1 <- BAD1[grepl("phs000572_201612", BAD1$PR),]
phs000572_201612_LOOKUP <- LOOKUP[which(LOOKUP$IID %in% phs000572_201612_BAD1$Original_ID),c("IID", "FID", "Sex", "STATUS", "F.SEXCHECK.", "SCOPE", "Project")]
phs000572_201612_LOOKUP <- phs000572_201612_LOOKUP[grepl("WES", phs000572_201612_LOOKUP$SCOPE),]
# first find exact match
indexed <- match(phs000572_201612_LOOKUP$IID, phs000572_201612_BAD1$Original_ID)
phs000572_201612_GOOD2 <- cbind(phs000572_201612_BAD1[indexed, c("Original_ID", "VCF_ID", "PR")], phs000572_201612_LOOKUP)
phs000572_201612_BAD2 <- phs000572_201612_BAD1[which(! phs000572_201612_BAD1$Original_ID %in% phs000572_201612_GOOD2$Original_ID),]



GOOD2 <- rbind(Macrogen_WGS_GOOD2, LOAD_WES_GOOD2, Genentech_WGS_GOOD2, Genentech_WES_GOOD2, Broad_WGS_GOOD2, TGI_WES_GOOD2, MGI_FASeEOAD_201605_GOOD2, Otogenetics_WES_GOOD2, phs000572_201508_GOOD2, phs000572_201612_GOOD2)
dim(GOOD2)
BAD2 <- rbind(Macrogen_WGS_BAD2, LOAD_WES_BAD2, Genentech_WGS_BAD2, Genentech_WES_BAD2, Broad_WGS_BAD2, TGI_WES_BAD2, MGI_FASeEOAD_201605_BAD2, Otogenetics_WES_BAD2, phs000572_201508_BAD2, phs000572_201612_BAD2)
BAD2 <- rbind(BAD2, BAD1[grepl("gDNA_SHERMAN", BAD1$PR),])
dim(BAD2)



ADSP <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/10026-subject_info-AN.csv", header = TRUE, sep =",", check.names = FALSE)
# GOOD3 <- merge(BAD2[c("Original_ID", "VCF_ID", "PR")], ADSP, by.x = "Original_ID", by.y = "Subject")
ADSP <- ADSP[na.omit(match(BAD2$Original_ID, ADSP$Subject)),]
tmpBAD <- BAD2[which(BAD2$Original_ID%in%ADSP$Subject),]
GOOD3 <- cbind(tmpBAD[which(tmpBAD$Original_ID%in%ADSP$Subject), c("Original_ID", "VCF_ID", "PR")], ADSP)
BAD3 <- BAD2[which(! BAD2$Original_ID %in% GOOD3$Subject),]

# SHERMAN
SHERMAN <- read.table("/40/AD/AD_Seq_Data/01.-RawData/201907_USUHS_gDNA_SHERMAN/fastq.sampleid.key.csv", header = TRUE, sep = ",")

BAD3$SHERMAN_BARCODE <- "EMPTY"
# check all samples are included
# sum(wantID$Original_ID %in% c(GOOD1$Original_ID, GOOD2$Original_ID, GOOD3$Original_ID,BAD3$Original_ID))

BAD3$SHERMAN_BARCODE <- SHERMAN$DepositID[match(BAD3$Original_ID, SHERMAN$TAGCID)]


ADSP[grepl ("A-LOAD-LD000988", ADSP$Subject),]



write.table(GOOD1, file ="pheno_GOOD_LIST1.csv", quote = FALSE, sep = ",", col.names = TRUE )
write.table(GOOD2, file ="pheno_GOOD_LIST2.csv", quote = FALSE, sep = ",", col.names = TRUE )
write.table(GOOD3, file ="pheno_GOOD_LIST3.csv", quote = FALSE, sep = ",", col.names = TRUE )
write.table(BAD3, file ="Aquilla_pheno_BAD_LIST.csv", quote = FALSE, sep = ",", col.names = TRUE )



## After finalizing
final.df <- read.csv("Final_Pheno_FASe_3894.csv", header = TRUE, sep =",", check.names = FALSE)


dim(final.df)
head(final.df)

length(unique(as.character(final.df$VCF_ID)))
# wantID$VCF_ID[which(!wantID$VCF_ID %in% final.df$VCF_ID)]
# # "7_108_10026_R1" "10R_R8_42_R1"   "25_10_183_R1"   "25_19_383_R1"   "10J_58_1_R1" 
# wantID[grepl("7_108_10026_R1", wantID$VCF_ID),]
# wantID[grepl("10R_R8_42_R1", wantID$VCF_ID),]
# wantID[grepl("10R_R8_42_R1", wantID$VCF_ID),]
# wantID[grepl("25_10_183_R1", wantID$VCF_ID),]
# wantID[grepl("25_19_383_R1", wantID$VCF_ID),]
# wantID[grepl("10J_58_1_R1", wantID$VCF_ID),]

# Family size
# We need to remove the 24 size family
table(table(final.df$FID))
# 1    2    3    4    5    6    7    8    9   10   24 
# 2142   16   98  136   74   37   24    2    8    1    1 

sort(table(final.df$FID), decreasing = T)
final.df[grepl("203" , final.df$FID),]


## Bias due to batch effect
setwd("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files3")
bias.df <- read.csv("Final_Pheno_FASe_3894.pheno", header = TRUE, sep ="\t", check.names = FALSE)
head(bias.df)
new.cols <- unique(as.character(bias.df$Project))
bias.df[, new.cols] <- "1"
head(bias.df)
df <- bias.df
for (i in 1:nrow(df)){
  df[i,grep(df$Project[i], colnames(df))] <- "2"
}

# remove IDs that got dropped in filtered .Fam file

filtered.fam <- read.csv("FASe_3894_WXS_SNPS_INDELS-hwe-geno0.05-mind0.1-WXSm.fam", header = FALSE, sep =" ", check.names = FALSE)
dim(filtered.fam)

df <- df[which(df$IID %in% filtered.fam$V2),]

# key.fam <- paste(filtered.fam$V1, filtered.fam$V2, sep = ":")
# key.df <- paste(df$FID, df$IID, sep = ":")
# sum(sort(key.fam) == sort (key.df))
# 
# metadata.df <- read.csv("Final_Pheno_FASe_3894.txt", header = TRUE, sep ="\t", check.names = FALSE)
# metadata.df <-  metadata.df[ which(metadata.df$IID %in% df$IID),]
# key.meta.df <- paste(metadata.df$FID, metadata.df$IID, sep = ":")
# sum(sort(key.meta.df) == sort (key.df))

write.table(df, file ="FASe_3894_Phenoscope.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

# Percentage by projects
diff.df <- data.frame(sort((table(wantID$PR)/3894)*100, decreasing = T))


####################recode sex and FID
sex_recode <- read.table("SEX_recode.txt", header = T)
dim(sex_recode)
PHENO <- read.table("Final_Pheno_FASe_3894.pheno", header = T)
head(PHENO)
dim(PHENO)

PHENO[grep("62033", PHENO$IID),]

match(PHENO$IID, sex_recode$IID)

setwd("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files2/test_original/")

FAM <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files2/test_original/FASe_3894_WXS_SNPS_INDELS.fam", header = F)
PhenoKey <- paste(PHENO$IID,PHENO$IID, sep = ":")
FAMkey <- paste(FAM$V1,FAM$V2, sep = ":")
grep("62033",FAM$V1)

# https://zzz.bwh.harvard.edu/plink/dataman.shtml
# plink --bfile mydata --update-ids recoded.txt --make-bed --out mydata2
# changes ID codes for individuals specified in recoded.txt, which should be in the format of four columnds per row: old FID, old IID, new FID, new IID, e.g.
NewFam <- cbind.data.frame(oldFID=FAM$V1,oldIID=FAM$V2)

head(PHENO)
PHENO$KeyID <- PHENO$IID
recodeIID <- merge(NewFam, PHENO, by.x = "oldIID", by.y = "KeyID")
head(recodeIID)

recodeIID[grep("62033",recodeIID$IID),]

dim(recodeIID)


recodeFID <- cbind.data.frame(oldFID=recodeIID$oldFID, oldIID=recodeIID$oldFID, newFID=recodeIID$FID, newIID=recodeIID$IID)

recodeSex <- cbind.data.frame(FID=recodeIID$FID, IID=recodeIID$IID, SEX=recodeIID$SEX)

write.table(recodeFID, file ="FID_IID_recode.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(recodeSex, file ="SEX_recode.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

################
# sex_recode$PhenoSex <- final.df$SEX [match(sex_recode$IID, final.df$VCF_ID)]
# sum(sex_recode$PhenoSex == sex_recode$SEX)
# sex_recode[grep ("62547",sex_recode$IID),]
# discordant.sex[grep("62547",discordant.sex$IID),]
# test_sex <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files2/original_bfile/FASe_3894_WXS_SNPS_INDELS_sex.sexcheck", header = TRUE)
# 
# dim(test_sex)
# varTestSex <- test_sex$SNPSEX[match(discordant.sex$IID, test_sex$IID)]
# sum(discordant.sex$PEDSEX == varTestSex)








# # Check Sex
# check.sex <- read.csv("FASe_3894_WXS_SNPS_INDELS.fam", header = FALSE, sep =" ", check.names = FALSE)
# head(check.sex)
# colnames(check.sex) <- c("FID", "IID", "V3", "V4", "Sex", "Status")
# # sorted.sexID <- sort(check.sex$IID) 
# # sorted.final.dfID <- sort(final.df$VCF_ID) 
# # sum(sorted.sexID == sorted.final.dfID)
# 
# fam.Sex <- paste(check.sex$IID, check.sex$Sex, sep =":")
# final.Sex <- paste(final.df$IID, final.df$SEX, sep =":")
# sum(final.Sex %in% fam.Sex)
# 
# 
# ## check which project does discordant sex belong to
# # discordant.sex <- read.delim("discordant_sex_FASe_3894_WXS_SNPS_INDELS-hwe-geno0.05-mind0.1-WXSm.tsv", header = FALSE, sep ="\t", check.names = FALSE)
# discordant.sex <- read.table("discordant_sex_FASe_3894_WXS_SNPS_INDELS-hwe-geno0.05-mind0.1-WXSm.csv", header = TRUE, sep =",", check.names = FALSE)
# 
# 
# dim(discordant.sex)
# head(discordant.sex)
# 
# colnames(discordant.sex) <- c("FID", "IID", "PEDSEX", "SNPSEX", "STATUS", "INBREED_COEFF" )
# discordant.sex$PEDSEX
# table(final.df$Project[which(final.df$VCF_ID %in% discordant.sex$IID)])
# 
# # as.character(discordant.sex$IID) [1:10]
# # as.character(final.df$VCF_ID)[1:10]
# 
# PhenoSex <-final.df$SEX[match(discordant.sex$IID, final.df$VCF_ID)]
# sum(discordant.sex$PEDSEX %in% PhenoSex)
# dim(discordant.sex)
# 
# discordant.sex$Project <- final.df$Project [match(discordant.sex$IID, final.df$VCF_ID)]
# discordant.sex[grep("8_64042_14", discordant.sex$IID),]
# 
# write.table(discordant.sex, file ="discordant_sex_FASe_3894_WXS_SNPS_INDELS-hwe-geno0.05-mind0.1-WXSm_V2.csv", quote = FALSE, sep = ",", col.names = TRUE, row.names = FALSE)



#####################################################################################
#########################IBD Filter for duplicates ##################################
#####################################################################################

library(ggplot2)
library(plotly)
library(ggrepel)
# setwd("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files2/test_original/")
setwd("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files3/")
# IBD<- fread("FASe_3894_WXS_SNPS_INDELS-hwe-geno0.05-mind0.1-IBD.genome")
IBD<- fread("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic_post_QC_with_STATUS-hwe-geno0.05-mind0.1-IBD.genome")
dim(IBD)
head(IBD)
# ## parent-offspring : 316
# sum(IBD$Z0 < 0.25 &
#       IBD$Z1 < 0.75) 

# remove non-Europeans
# PCAs<- read.table("FASe_3894_WXS_SNPS_INDELS-hwe-geno0.05-mind0.1-PCAS.eigenvec", header=T)
PCAs<- read.table("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic_post_QC_with_STATUS-hwe-geno0.05-mind0.1-PCAS.eigenvec", header=T)

# With HAPMAP
PCAs<- read.table("FASe_3894_WXS_SNPS_INDELS-hwe-geno0.05-mind0.1-HAPMAP-MERGED3-for_PCA.eigenvec-PC1-PC2-COHORT.txt", header=T)
p <- ggplot(PCAs, aes(x=PC1, y=PC2, color=COHORT)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Aquilla 3891") +
  scale_color_manual(values = c('green', 'black', 'red', "blue"))

ggplotly(p)

# PCAs<- read.table("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic_post_QC_with_STATUS-hwe-geno0.05-mind0.1-HAPMAP-MERGED3-for_PCA.eigenvec", header=T)
# PCAs<- read.table("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic_post_QC_with_STATUS-hwe-geno0.05-mind0.1_post_QC-HAPMAP-MERGED3-for_PCA.eigenvec", header=T)
# p <- ggplot(PCAs, aes(x=PC1, y=PC2)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Aquilla_FASe 3891")
# p


# # Any unwanted points, select interactively
# ggplotly(p)

# PCA_NHW <- PCAs[(PCAs$PC1 > -0.02 &
#                    PCAs$PC1 < 0.025 &
#                    PCAs$PC2 < 0.065 &
#                    PCAs$PC2 > -0.025),] 

# PCA_NHW <- PCAs[(PCAs$PC1 > -0.1 &
#                    PCAs$PC1 < 0.05 &
#                    PCAs$PC2 < 0.045 &
#                    PCAs$PC2 > -0.1),] 

PCA_NHW <- PCAs[(PCAs$PC1 > -0.001 &
                   PCAs$PC1 < 0.02 &
                   PCAs$PC2 < 0.010 &
                   PCAs$PC2 > -0.01),]



p <- ggplot(PCA_NHW, aes(x=PC1, y=PC2, color=COHORT)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Aquilla 3891") +
  scale_color_manual(values = c('green', 'black', 'red', "blue"))

ggplotly(p)

PCA_NHW <- PCA_NHW[grepl ("FASe", PCA_NHW$COHORT ),]

p<- ggplot(PCA_NHW, aes(x=PC1, y=PC2)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Aquilla_FASe 3891") +
geom_text_repel(aes(label=ifelse((PC1< (-0.014)) & (PC2  < 0.025),as.character(IID),''))) 
p

# by FID
p<- ggplot(PCA_NHW, aes(x=PC1, y=PC2)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Aquilla_FASe 3891") +
  geom_text_repel(aes(label=ifelse((PC1< (-0.014)) & (PC2  < 0.025),as.character(FID),''))) 
p


# # Do it again removing more specific samples
# removedFID <- as.character(PCA_NHW$FID[(PCA_NHW$PC1< (-0.014)) & (PCA_NHW$PC2  < 0.025)])
# removedIID <- as.character(PCA_NHW$IID[(PCA_NHW$PC1< (-0.014)) & (PCA_NHW$PC2  < 0.025)])
# Do it again removing more specific samples
removedFID <- as.character(PCA_NHW$FID[(PCA_NHW$PC1< (-0.014)) & (PCA_NHW$PC2  < 0.025)])
removedIID <- as.character(PCA_NHW$IID[(PCA_NHW$PC1< (-0.014)) & (PCA_NHW$PC2  < 0.025)])  

sum(PCA_NHW$FID %in% removedFID)
sum(PCA_NHW$IID %in% removedIID)


PCA_NHW <- PCA_NHW[which(!PCA_NHW$FID %in% removedFID),]

# Now Plot again
p<- ggplot(PCA_NHW, aes(x=PC1, y=PC2)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Aquilla_FASe 3891") +
  geom_text_repel(aes(label=ifelse((PC1< (-0.014)) & (PC2  < 0.025),as.character(IID),''))) 
p

wanted.famID <- unique(as.character(PCA_NHW$FID))

wantedPCA_NHW <- {}
for (i in 1:length(wanted.famID)){
tmp <- PCAs[which(wanted.famID[i]== PCAs$FID),]
wantedPCA_NHW <- rbind(wantedPCA_NHW, tmp)
}
dim(wantedPCA_NHW)
# 3747

p <- ggplot(wantedPCA_NHW, aes(x=PC1, y=PC2), label = as.factor(IID)) + geom_point() + xlab("PC1") + ylab("PC2") 
p
ggplotly(p)


write.csv(wantedPCA_NHW, "Aquilla_FASe_NHW_final.csv", quote = FALSE, row.names = FALSE)


# check again 
PCAs<- read.table("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic_post_QC_with_STATUS-hwe-geno0.05-mind0.1_post_QC-PCAS.eigenvec", header=T)
p <- ggplot(PCAs, aes(x=PC1, y=PC2), label = as.factor(IID)) + geom_point() + xlab("PC1") + ylab("PC2") 
p
ggplotly(p)


View(PCAs)

which(!PCAs$IID %in% PCA_NHW$IID )

FILE="FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic_post_QC_with_STATUS-hwe-geno0.05-mind0.1_post_QC-HAPMAP-MERGED3-for_PCA.eigenvec"
PCAs<- read.table("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic_post_QC_with_STATUS-hwe-geno0.05-mind0.1_post_QC-HAPMAP-MERGED3-for_PCA.eigenvec-PC1-PC2-COHORT-After-removing-NHW.txt", header=T)
##plotting:

# This sample also seems to be non-european
# 17_18_3
newPCAs <- PCAs[!grepl("17_18_3", PCAs$IID),]
p <- ggplot(newPCAs, aes(x=PC1, y=PC2, color=COHORT)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Aquilla 3891") +
  scale_color_manual(values = c('green', 'black', 'red', "blue"))  

ggplotly(p)

# Remove more points

library(ggplot2)
#read data:
# PCAs<-read.table("aurora_7234_WXS_SNPS_INDELS-hwe-geno0.05-mind0.1-PCAS.eigenvec", header=T)
PCAs<- read.table("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic_post_QC_with_STATUS-hwe-geno0.05-mind0.1_post_QC-PCAS.eigenvec", header=T)
##plotting:
p <- ggplot(PCAs, aes(x=PC1, y=PC2)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Aquilla_FASe 3891")
ggplotly(p)


PCAs<- read.table("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic_post_QC_with_STATUS-hwe-geno0.05-mind0.1_post_QC-HAPMAP-MERGED3-for_PCA.eigenvec-PC1-PC2-COHORT-After-removing-NHW.txt", header=T)
##plotting:
ggplot(PCAs, aes(x=PC1, y=PC2, color=COHORT)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Aquilla 3891") +
  scale_color_manual(values = c('green', 'black', 'red', "blue"))  

# Do it again removing more specific samples
removedFID <- as.character(PCAs$FID[(PCAs$PC1 > (-0.05)) & (PCAs$PC2  < 0.04)])
removedIID <- as.character(PCA_NHW$IID[(PCA_NHW$PC1< (-0.014)) & (PCA_NHW$PC2  < 0.025)])  

sum(PCA_NHW$FID %in% removedFID)
sum(PCA_NHW$IID %in% removedIID)
PCA_NHW <- PCA_NHW[which(!PCA_NHW$FID %in% removedFID),]

#####################################

# df[rowSums(apply(df[, 1:2], 2, function(y) y %in% x)) == 2,]



parent_offsping <- IBD[(IBD$Z0 < 0.25 & IBD$Z1 > 0.75), ]


parent_offsping$Relationship <- "parent-offspring"
dim(parent_offsping)
# 985
# 105

sibPairs <- IBD[(IBD$Z0 < 0.5 &
                   IBD$Z0 > 0.10 &
                   IBD$Z1 < 0.75 &
                   IBD$Z1 > 0.25),] 


sibPairs$Relationship <- "sib-pairs"
dim(sibPairs)
# 5556
# 2360

# ## possible duplicates : 201
# sum(IBD$Z0 < 0.25 &
#       IBD$Z1 < 0.25) 

duplicates <- IBD[(IBD$Z0 < 0.25 &
                     IBD$Z1 < 0.25), ]

duplicates$Relationship <- "duplicates"
dim(duplicates)
# 159
# 116
relatives.ALL <- rbind(parent_offsping, sibPairs, duplicates)

# # Add pheno info
# relatives.ALL <- cbind(relatives.ALL, final.df[match(relatives.ALL$IID1, final.df$VCF_ID), c("VCF_ID", "SEX", "STATUS")], final.df[match(relatives.ALL$IID2, final.df$VCF_ID), c("VCF_ID", "SEX", "STATUS")])

# Extract only NHW samples from relatives.ALL pairs

# all(apply(df[, 1:2], 2, function(y) all(x %in% y)))
# x <- c("A", "B", "D")
# df <- structure(c("A", "B", "C", "D", "B", "A", "D", "E", "3", "4", 
#                   "5", "6"), .Dim = 4:3, .Dimnames = list(NULL, c("A1", "B1", "C1"
#                   )))
# df[rowSums(apply(df[, 1:2], 2, function(y) y %in% x)) == 2,]
# tt <- relatives.ALL[rowSums(apply(relatives.ALL[, c("IID1", "IID2")], 2, function(y) y %in% wantedPCA_NHW$IID)) == 2,]

#### Label the ethnicity first
relatives.ALL$IID1_Ethnicity <- ifelse(relatives.ALL$IID1 %in% wantedPCA_NHW$IID, "NHW", "Other")
relatives.ALL$IID2_Ethnicity <- ifelse(relatives.ALL$IID2 %in% wantedPCA_NHW$IID, "NHW", "Other")

## Test
grep("62038", wantedPCA_NHW$IID) # NHW
grep("25_20_401", wantedPCA_NHW$IID) #Other
grep("MAP_62038", wantedPCA_NHW$FID) # NHW
grep("25_20_401", wantedPCA_NHW$FID) #Other

write.csv(relatives.ALL, "Related_IBD_Pairs.csv", quote = FALSE, row.names = FALSE)


Non.duplicate.relatives <- relatives.ALL[!grepl("duplicates",relatives.ALL$Relationship),]


################################
# Now get the final list of samples to keep in the analysis
# remove.df <- read.table("Samples_to_removed_after_IBD_FASe_3891.csv", header = F)
# remove.df <- read.table("Samples_to_remove_after_IBD_FASe.csv", header = F)
remove.df <- read.table("duplicates_to_remove.txt", header = F)
removeSamples <- unique(as.character(remove.df$V1))
length(removeSamples)

FinalWantedSamples <- wantedPCA_NHW[!wantedPCA_NHW$IID %in% removeSamples,1:2]


# Remove samples with discordant sex
discordant.sex <- c("0_62519_39", "10R_R8_36", "27_164_85644", "26_FQX_FQX77608")
FinalWantedSamples <- FinalWantedSamples[-match(discordant.sex, FinalWantedSamples$IID),]
dim(FinalWantedSamples)
table(table(as.character(FinalWantedSamples$FID)))


# write.table(FinalWantedSamples, "Samples_to_keep_for_Aquilla_FASe_3606.txt", sep="\t", quote = FALSE, row.names = FALSE)
write.table(FinalWantedSamples, "Samples_to_keep_for_Aquilla_FASe_1608.txt", sep="\t", quote = FALSE, row.names = FALSE)
table(table(as.character(PCA_NHW$FID)))
table(table(as.character(FinalWantedSamples$FID)))

PHENO <- read.table ("FASe_3894_Phenoscope.txt", header = TRUE)
head(PHENO)
PHENO$KeyID <- PHENO$IID

## Family within PHENO 3891; families with 7 families with household size 9 were removed after PC and IBD
table(table(PHENO$FID))
# 1    2    3    4    5    6    7    8    9   10   24 
# 2141   16   98  136   74   39   22    2    8    1    1 

# Families with household 9
which(table(PHENO$FID)==9)
# 0_3103  0_3251 10R_R78  27_131   4_162   4_553   4_558   4_649 
# 46      75     137     315     353     361     363     381

FINAL_PHENO_AFTER_QC <- merge(FinalWantedSamples, PHENO, by.x = "IID", by.y = "KeyID")

FINAL_PHENO_AFTER_QC <- FINAL_PHENO_AFTER_QC[-(1:2)]
colnames(FINAL_PHENO_AFTER_QC) [1:2] <- c("FID", "IID")

write.table(FINAL_PHENO_AFTER_QC, "After_QC_Pheno_Aquilla_FASe_3606.txt", sep = "\t", quote = FALSE, row.names = FALSE)

####################################################### After QC 
FINAL_PHENO_AFTER_QC <- read.table("After_QC_Pheno_Aquilla_FASe_3606.txt", header = T)

Sam1 <- relatives.ALL$IID1[grep("Other",relatives.ALL$IID1_Ethnicity)]
Sam2 <- relatives.ALL$IID2[grep("Other",relatives.ALL$IID2_Ethnicity)]
SAMPLES <- c(Sam1, Sam2)
# This be zero as there should not be any non-NHW
sum(as.character(FINAL_PHENO_AFTER_QC$IID) %in% SAMPLES)

length(FINAL_PHENO_AFTER_QC$FID)

table(table(as.character(FINAL_PHENO_AFTER_QC$FID)))
# > table(table(FINAL_PHENO_AFTER_QC$FID))
# 0    1    2    3    4    5    6    7    8    9   24 
# 100 2080   14   96  131   64   30   16    4    2    1
sum(table(table(FINAL_PHENO_AFTER_QC$FID)))
# 

# Vicky suggested to remove family 24 members (FID="LD0241F")
FINAL_PHENO_AFTER_QC <- FINAL_PHENO_AFTER_QC[!grepl("^203$", FINAL_PHENO_AFTER_QC$FID),]

# remove individuals from household size 1 (less than 2)
families_to_include <- names(which(table(FINAL_PHENO_AFTER_QC$FID) >1))
FINAL_PHENO_AFTER_QC <- FINAL_PHENO_AFTER_QC[FINAL_PHENO_AFTER_QC$FID %in% families_to_include,]

# Sex
table(FINAL_PHENO_AFTER_QC$SEX )
# 1   2 
# 584 918 

# CACO
table(FINAL_PHENO_AFTER_QC$STATUS)
# -9    1    2 
# 147  286 1069 

table(table(as.character(FINAL_PHENO_AFTER_QC$FID)))

head(FINAL_PHENO_AFTER_QC)
dim(FINAL_PHENO_AFTER_QC)
write.table(FINAL_PHENO_AFTER_QC, "After_QC_Pheno_Demographic_Aquilla_FASe_1502.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Final list of samples to include in the analysis (excluding families from household 1)
head(FINAL_PHENO_AFTER_QC[1:2])
write.table(FINAL_PHENO_AFTER_QC[1:2], "After_QC_Final_list_Aquilla_FASe_1502.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# # Phenotype for rvtest
phenoPed <- FINAL_PHENO_AFTER_QC[1:2]
phenoPed$fatid <- "0"
phenoPed$matid <- "0"
phenoPed <- cbind(phenoPed, SEX=FINAL_PHENO_AFTER_QC$SEX, STATUS=FINAL_PHENO_AFTER_QC$STATUS)

write.table(phenoPed, "phenotype.ped", sep = "\t", quote = FALSE, row.names = FALSE)






# Now from this cluster of relatives.ALL, select samples who is younger/older
# (depending on the analysis) and whether are CA or CO based on the proportion
# of cases/controls in the whole cohort. Then after selecting the samples we
# want to include, we can can exclude the other set.

## PCA with 1502 Fase data only (without HAPMAP)
PCAs<- read.table("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1_post_QC2-PCAS.eigenvec", header=T)
dim(PCAs)


# Merge PCs to phenotype data
pheno <- read.table("phenotype1.ped" , header = T,  comment.char = '&', check.names = F)
dim(pheno)
head(pheno)
PhenoPC <- cbind(pheno, PCAs[match(pheno$IND_ID, PCAs$IID), c("IID", "PC1", "PC2", "PC3", "PC4")])
head(PhenoPC)
colnames(PhenoPC)[1:dim(pheno)[2]] <- colnames(pheno)
write.table(PhenoPC, "phenotype1_with_PC1-PC4_with_1502_samples_and_no_hapmap.ped", sep =" ", col.names = T, quote = F, row.names = F)


######################## Epacts results and plots

# https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
library(lattice)
source("https://raw.githubusercontent.com/achalneupane/rcodes/master/qqrunif.r")
library (qqman)
# LOGISTIC<-read.table("EPACTs_Complete_stats_single_variant_min_maf_0.txt",  comment.char = '&', header = TRUE, sep ="\t", check.names = F)


## Min MAF 0
# LOGISTIC <- fread("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_min-maf_0_EPACTs_Complete_stats_single_variant.txt", header = TRUE, sep ="\t", check.names = F)
LOGISTIC <- fread("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_min-maf_0_EPACTs_Complete_stats_single_variant.txt", header = TRUE, sep ="\t", check.names = F)

## MAX MAF 0.01
LOGISTIC <- fread("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_max-maf_0.01_EPACTs_Complete_stats_single_variant.txt", header = TRUE, sep ="\t", check.names = F)
## MAX MAF 0.001
# LOGISTIC <- fread("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_max-maf_0.001_EPACTs_Complete_stats_single_variant.txt", header = TRUE, sep ="\t", check.names = F)
## MMSKAT MAX MAF 0.01 
LOGISTIC <- fread("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr.gene.mmskat.epacts", header = TRUE, sep ="\t", check.names = F)


LOGISTIC <- as.data.frame(LOGISTIC)
# LOGISTIC <- LOGISTIC[LOGISTIC$MAF <= 0.01,]




head(LOGISTIC)
LOGISTIC <- LOGISTIC[!is.na(LOGISTIC$PVALUE),]
LOGISTIC$PVALUE <- as.numeric(LOGISTIC$PVALUE)


# my.pvalues <- LOGISTIC$PVALUE[!is.na(LOGISTIC$log10)]
# # Pvalues below MAF1%
# LOGISTIC_MAF_1Percent <- LOGISTIC[LOGISTIC$MAF > 0.005 & LOGISTIC$MAF <= 0.01,]
# my.pvalues_within_MAF1percent <- LOGISTIC_MAF_1Percent$PVALUE[!is.na(LOGISTIC_MAF_1Percent$PVALUE)]

my.pvalues <- LOGISTIC$PVALUE

# ## library(lattice);
# qqmath(~-log10(my.pvalues),
#        distribution=function(x){-log10(qunif(1-x))}
# );

# #Make plot
# #Calculate expectations
# exp.pvalues<-(rank(LOGISTIC$PVALUE, ties.method="first")+.5)/(length(my.pvalues)+1)
# plot(-log10(exp.pvalues), -log10(LOGISTIC$PVALUE), asp=1)
# abline(0,1)

# qqunif.plot(my.pvalues) # raw p-values

### MAX MAF 0.01
jpeg("QQ_plot_EPACTs_Complete_stats_single_variant_max_maf_0.01.jpeg", height = 20, width = 15, units='cm', res = 300)
### MIN MAF 0
# jpeg("QQ_plot_EPACTs_Complete_stats_single_variant_min_maf_0.jpeg", height = 20, width = 15, units='cm', res = 300)
### MMSKAT MAX MAF 0.01
jpeg("QQ_plot_EPACTs_Complete_stats_gene-wise_mmskat_max_maf_0.01.jpeg", height = 20, width = 15, units='cm', res = 300)

qqunif.plot(LOGISTIC$PVALUE)
dev.off()

# # qqunif.plot(my.pvalues) # raw p-values
# tiff("QQ_plot_EPACTs_Complete_stats_single_variant_min_maf_0_Pvalues_within_MAF1percent.jpeg", height = 20, width = 15, units='cm', compression = "lzw", res = 300)
# qqunif.plot(my.pvalues_within_MAF1percent)
# dev.off()



library (qqman) 

sum(my.pvalues < `^`(10,-8))




## Min maf=0
# TOP.genes <- read.table("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_min-maf_0.epacts.top5000",  comment.char = '&', header = TRUE, sep ="\t", check.names = F)
### Max MAF 0.01
# TOP.genes <- read.table("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_max-maf_0.01.epacts.top5000",  comment.char = '&', header = TRUE, sep ="\t", check.names = F)
### MMSKAT MAX MAF 0.01
TOP.genes <- read.table("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr.gene.mmskat.epacts.top5000",  comment.char = '&', header = TRUE, sep ="\t", check.names = F)

min(TOP.genes$PVALUE)
# 6.911e-11
head(TOP.genes, 20)
colnames(TOP.genes) [1] <- "CHR"

TOP.genes$CHRposition <- paste(TOP.genes$CHR, TOP.genes$BEG, sep = ":")

sum(TOP.genes$PVALUE < `^`(10,-6))
Significant <- TOP.genes[TOP.genes$PVALUE < `^`(10,-6),]

# Read Annotated variants

### Min MAF 0
# annotatedGenes <- read.table("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_min-maf_0.epacts_EPACTS_variants_to_search_result.txt", header = F, sep = "\t")
### Max MAF 0.01
# annotatedGenes <- read.table("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_max-maf_0.01.epacts_EPACTS_variants_to_search_result.txt", header = F, sep = "\t")
### MMSkat MAX MAF 0.01
annotatedGenes <- read.table("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr.gene.mmskat.epacts_EPACTS_variants_to_search_result.txt", header = F, sep = "\t")

head(annotatedGenes)
annotatedGenes <- annotatedGenes[c(3,8)]

Features <- cbind(as.character(annotatedGenes$V3), do.call(rbind, lapply(strsplit(as.character(annotatedGenes$V8),split='|',fixed=TRUE), `[`, c(2:4, 11))))
colnames(Features) <- c("Variant","Var_Type", "Effect", "Gene", "Pchange" )
Features <- cbind.data.frame(Features, varID=gsub(";.*","",Features[,1]))
head(Features)
# grep("77092549:A:G", Features$varID)





colnames(LOGISTIC) [1] <- "CHR" 
colnames(LOGISTIC) [2] <- "BP" 
colnames(LOGISTIC)[colnames(LOGISTIC)=="PVALUE"] <- "P" 
LOGISTIC$CHR <- as.character(LOGISTIC$CHR)
LOGISTIC$CHR[grepl("X", LOGISTIC$CHR) ] <- 23
LOGISTIC$CHR[grepl("Y", LOGISTIC$CHR) ] <- 24
LOGISTIC$CHR <- as.numeric(LOGISTIC$CHR)
LOGISTIC$SNP <- as.character(vapply(strsplit(LOGISTIC$MARKER_ID,"_"), `[`, 3, FUN.VALUE=character(1)))

############# Only for MM SKAT (Max MAF 0.01)
head(LOGISTIC)
# LOGISTIC$SNP <- as.character(vapply(strsplit(LOGISTIC$MARKER_ID,"-"), `[`, 1, FUN.VALUE=character(1)))
# LOGISTIC$SNP <- as.character(vapply(strsplit(LOGISTIC$SNP,"_"), `[`, 2, FUN.VALUE=character(1)))
LOGISTIC$SNP <- LOGISTIC$MARKER_ID
#############


LOGISTIC$SNP <- as.character(LOGISTIC$SNP)

TopFeatures <- merge(LOGISTIC, Features, by.x = "SNP", by.y = "varID")



## Min MAF 0
# write.table(TopFeatures, "Annotated_TOP_500_GENES_FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_min-maf_0.epacts.top500.csv", sep = ",", quote = FALSE, row.names = FALSE)
### MAX MAF 0.01
# write.table(TopFeatures, "Annotated_TOP_500_GENES_FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_max-maf_0.01.epacts.top500.csv", sep = ",", quote = FALSE, row.names = FALSE)
### MMSkat MAX MAF 0.01
write.table(TopFeatures, "Annotated_TOP_500_GENES_FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_MMSkat_max-maf_0.01.epacts.top500.csv", sep = ",", quote = FALSE, row.names = FALSE)



LOGISTIC.t <- merge(LOGISTIC, Features, by.x = "SNP", by.y = "varID", all.x = TRUE)

LOGISTIC.t$SNP[!is.na(LOGISTIC.t$Gene)] <- as.character(LOGISTIC.t$Gene[which(!is.na(LOGISTIC.t$Gene))] )


# na.omit(match (LOGISTIC$SNP,Features[,"varID"]))
# TopFeatures <- cbind(LOGISTIC[which(LOGISTIC$SNP %in% Features[,"varID"]),], Features)

# write.table(TopFeatures, "Top_genes_with_annotation_min_maf_0_from_epacts.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# 
# 
# Features[,"varID"] <- as.character(Features[,"varID"])
# wantedFeatures <- LOGISTIC$SNP[which(LOGISTIC$SNP %in% Features[,"varID"])] 
# wantedGene <- Features[,"Gene"] [which(Features[,"varID"] %in% wantedFeatures)]
# LOGISTIC$SNP[which(LOGISTIC$SNP %in% Features[,"varID"])]  <- as.character(wantedGene)


# # I will remove CHRX and Y completely
LOGISTIC.t <- LOGISTIC.t[!grepl("23|24",LOGISTIC.t$CHR ),]

# test <- LOGISTIC.t[order(LOGISTIC.t$P),]
# test <- test[1:100,c("SNP", "P", "CHR", "BP")]
# manhattan(test, main = "", ylim=c(0,15), col = c("blue4", "orange3"), suggestiveline = -log10(1e-06), genomewideline = -log10(1e-08), annotateTop = TRUE, annotatePval = -log10(1e-08), chrlabs = as.character(sort(unique(test$CHR)))) 


# tiff("Manhattan_plot_EPACTs_Complete_stats_single_variant_max_maf_0.01.tiff", height = 20, width = 30, units='cm', compression = "lzw", res = 300)

### Min MAF 0
# jpeg("Manhattan_plot_EPACTs_Complete_stats_single_variant_min_maf_0.jpeg", height = 20, width = 30, units='cm', res = 300)
### Max MAF 0.01
# jpeg("Manhattan_plot_EPACTs_Complete_stats_single_variant_max_maf_0.01.jpeg", height = 20, width = 30, units='cm', res = 300)
### MMSkat Max MAF 0.01
jpeg("Manhattan_plot_EPACTs_Complete_stats_single_variant_MMSkat_max_maf_0.01.jpeg", height = 20, width = 30, units='cm', res = 300)
LOGISTIC <- LOGISTIC[!grepl("23|24",LOGISTIC$CHR ),]
LOGISTIC.t <- LOGISTIC
LOGISTIC.t <- LOGISTIC.t[!is.na(LOGISTIC.t$CHR),]
LOGISTIC.t$SNP <- as.character(gsub(".*_","",LOGISTIC.t$SNP))
###

# manhattan(LOGISTIC, main = "", ylim=c(0,15), cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = -log10(1e-06), genomewideline = -log10(1e-08), annotatePval = -log10(1e-08), chrlabs = c(1:22, "X", "Y")) 
manhattan(LOGISTIC.t, main = "", ylim=c(0,15), col = c("blue4", "orange3"), suggestiveline = -log10(1e-06), genomewideline = -log10(1e-08), annotateTop = TRUE, annotatePval = -log10(1e-06), chrlabs = as.character(1:22)) 
dev.off()


# source("https://raw.githubusercontent.com/achalneupane/rcodes/master/manhattanplot.r")
# # https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_Manhattan_Plots_in_R
# manhattan.plot(LOGISTIC$CHR, LOGISTIC$BP, LOGISTIC$P)

################################################################################################
######################## Epacts results and plots (with SEX, PC1, PC2, PC3) ####################
################################################################################################
# https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
setwd("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files2")
library(lattice)
source("https://raw.githubusercontent.com/achalneupane/rcodes/master/qqrunif_with_lambda.r")
library (qqman)

# This is a helper functoin to calculate LAMBDA
inflation <- function(pvalues) {
  chisq <- qchisq(1 - pvalues, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  return(lambda)
}


# LOGISTIC<-read.table("EPACTs_Complete_stats_single_variant_min_maf_0.txt",  comment.char = '&', header = TRUE, sep ="\t", check.names = F)


## 1. Min MAF 0
LOGISTIC <- fread("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_PC3_with_no_HAPMAP_min-maf_0_EPACTs_Complete_stats_single_variant.txt", header = TRUE, sep ="\t", check.names = F)
## MIN MAF > 0.01
LOGISTIC <- LOGISTIC[LOGISTIC$MAF > 0.01,]

## 2.  MAX MAF 0.01
LOGISTIC <- fread("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_PC3_with_no_HAPMAP_max-maf_0.01_EPACTs_Complete_stats_single_variant.txt", header = TRUE, sep ="\t", check.names = F)


## 3. MAX MAF 0.01 mmSKAT
## MMSKAT MAX MAF 0.01 
LOGISTIC <- fread("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_SEX_PC1_PC2_PC3_without_HAPMAP.gene.mmskat_EPACTs_Complete_stats_single_variant.txt", header = TRUE, sep ="\t", check.names = F)


LOGISTIC <- as.data.frame(LOGISTIC)
# LOGISTIC <- LOGISTIC[LOGISTIC$MAF <= 0.01,]

head(LOGISTIC)
LOGISTIC <- LOGISTIC[!is.na(LOGISTIC$PVALUE),]
LOGISTIC$PVALUE <- as.numeric(LOGISTIC$PVALUE)

my.pvalues <- LOGISTIC$PVALUE



LAMBDA <- round(inflation(my.pvalues), digits = 3)
LAMBDA

### MIN MAF 0
jpeg("QQ_plot_EPACTs_Complete_stats_single_variant_min_maf_0_SEX_PC1_PC2_PC3.jpeg", height = 20, width = 15, units='cm', res = 300)
jpeg("QQ_plot_EPACTs_Complete_stats_single_variant_min_maf_0.01_SEX_PC1_PC2_PC3.jpeg", height = 20, width = 15, units='cm', res = 300)
### MAX MAF 0.01
jpeg("QQ_plot_EPACTs_Complete_stats_single_variant_max_maf_0.01_SEX_PC1_PC2_PC3.jpeg", height = 20, width = 15, units='cm', res = 300)
### MMSKAT MAX MAF 0.01
jpeg("QQ_plot_EPACTs_Complete_stats_gene-wise_mmskat_max_maf_0.01_SEX_PC1_PC2_PC3.jpeg", height = 20, width = 15, units='cm', res = 300)

qqunif.plot(LOGISTIC$PVALUE, LAMBDA= LAMBDA)
dev.off()

# # qqunif.plot(my.pvalues) # raw p-values
# tiff("QQ_plot_EPACTs_Complete_stats_single_variant_min_maf_0_Pvalues_within_MAF1percent.jpeg", height = 20, width = 15, units='cm', compression = "lzw", res = 300)
# qqunif.plot(my.pvalues_within_MAF1percent)
# dev.off()



library (qqman) 

sum(my.pvalues < `^`(10,-8))




## Min maf=0
TOP.genes <- read.table("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_PC3_with_no_HAPMAP_min-maf_0.epacts.top5000",  comment.char = '&', header = TRUE, sep ="\t", check.names = F)
### Max MAF 0.01
TOP.genes <- read.table("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_PC3_with_no_HAPMAP_max-maf_0.01.epacts.top5000",  comment.char = '&', header = TRUE, sep ="\t", check.names = F)
### MMSKAT MAX MAF 0.01
TOP.genes <- read.table("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_SEX_PC1_PC2_PC3_without_HAPMAP.gene.mmskat.epacts.top5000",  comment.char = '&', header = TRUE, sep ="\t", check.names = F)

min(TOP.genes$PVALUE)
# 6.911e-11
head(TOP.genes, 20)
colnames(TOP.genes) [1] <- "CHR"

TOP.genes$CHRposition <- paste(TOP.genes$CHR, TOP.genes$BEG, sep = ":")

sum(TOP.genes$PVALUE < `^`(10,-6))
Significant <- TOP.genes[TOP.genes$PVALUE < `^`(10,-6),]

# Read Annotated variants

### Min MAF 0
annotatedGenes <- read.table("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_PC3_with_no_HAPMAP_min-maf_0.epacts_EPACTS_variants_to_search_result.txt", header = F, sep = "\t")
### Max MAF 0.01
annotatedGenes <- read.table("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_PC3_with_no_HAPMAP_max-maf_0.01.epacts_EPACTS_variants_to_search_result.txt", header = F, sep = "\t")
### MMSkat MAX MAF 0.01
annotatedGenes <- read.table("FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_SEX_PC1_PC2_PC3_without_HAPMAP.gene.mmskat.epacts_EPACTS_variants_to_search_result.txt", header = F, sep = "\t")


head(annotatedGenes)
annotatedGenes <- annotatedGenes[c(3,8)]

Features <- cbind(as.character(annotatedGenes$V3), do.call(rbind, lapply(strsplit(as.character(annotatedGenes$V8),split='|',fixed=TRUE), `[`, c(2:4, 11))))
colnames(Features) <- c("Variant","Var_Type", "Effect", "Gene", "Pchange" )
Features <- cbind.data.frame(Features, varID=gsub(";.*","",Features[,1]))
head(Features)
# grep("77092549:A:G", Features$varID)





colnames(LOGISTIC) [1] <- "CHR" 
colnames(LOGISTIC) [2] <- "BP" 
colnames(LOGISTIC)[colnames(LOGISTIC)=="PVALUE"] <- "P" 
LOGISTIC$CHR <- as.character(LOGISTIC$CHR)
LOGISTIC$CHR[grepl("X", LOGISTIC$CHR) ] <- 23
LOGISTIC$CHR[grepl("Y", LOGISTIC$CHR) ] <- 24
LOGISTIC$CHR <- as.numeric(LOGISTIC$CHR)
LOGISTIC$SNP <- as.character(vapply(strsplit(LOGISTIC$MARKER_ID,"_"), `[`, 3, FUN.VALUE=character(1)))
LOGISTIC$Variant_full <- LOGISTIC$SNP
LOGISTIC$Variant_full <- LOGISTIC$SNP
# ############# Only for MM SKAT (Max MAF 0.01)
head(LOGISTIC)
# LOGISTIC$SNP <- as.character(vapply(strsplit(LOGISTIC$MARKER_ID,"-"), `[`, 1, FUN.VALUE=character(1)))
# LOGISTIC$SNP <- as.character(vapply(strsplit(LOGISTIC$SNP,"_"), `[`, 2, FUN.VALUE=character(1)))
LOGISTIC$SNP <- LOGISTIC$MARKER_ID
#############


LOGISTIC$SNP <- as.character(LOGISTIC$SNP)

TopFeatures <- merge(LOGISTIC, Features, by.x = "SNP", by.y = "varID")



## Min MAF 0
write.table(TopFeatures, "Annotated_TOP_500_GENES_FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_PC3_without_hapmap_min-maf_0.epacts.top500.csv", sep = ",", quote = FALSE, row.names = FALSE)
### MAX MAF 0.01
write.table(TopFeatures, "Annotated_TOP_500_GENES_FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_PC3_without_hapmap_max-maf_0.01.epacts.top500.csv", sep = ",", quote = FALSE, row.names = FALSE)
### MMSkat MAX MAF 0.01
# write.table(TopFeatures, "Annotated_TOP_500_GENES_FASe_1502_from_Aquilla_after_QC-biallelic_FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_with_no_chr_epacts_single_variant_cov_SEX-PC1-PC2_PC3_without_hapmap_MMSkat_max-maf_0.01.epacts.top500.csv", sep = ",", quote = FALSE, row.names = FALSE)



LOGISTIC.t <- merge(LOGISTIC, Features, by.x = "SNP", by.y = "varID", all.x = TRUE)

LOGISTIC.t$SNP[!is.na(LOGISTIC.t$Gene)] <- as.character(LOGISTIC.t$Gene[which(!is.na(LOGISTIC.t$Gene))] )
# # I will remove CHRX and Y completely
LOGISTIC.t <- LOGISTIC.t[!grepl("23|24",LOGISTIC.t$CHR ),]

# # test <- LOGISTIC.t[order(LOGISTIC.t$P),]
# # test <- test[1:100,c("SNP", "P", "CHR", "BP")]
# # manhattan(test, main = "", ylim=c(0,15), col = c("blue4", "orange3"), suggestiveline = -log10(1e-06), genomewideline = -log10(1e-08), annotateTop = TRUE, annotatePval = -log10(1e-08), chrlabs = as.character(sort(unique(test$CHR)))) 

############# Only for MM SKAT (Max MAF 0.01)
LOGISTIC <- LOGISTIC[!grepl("23|24",LOGISTIC$CHR ),]
LOGISTIC.t <- LOGISTIC
LOGISTIC.t <- LOGISTIC.t[!is.na(LOGISTIC.t$CHR),]
LOGISTIC.t$SNP <- as.character(gsub(".*_","",LOGISTIC.t$SNP))
############

Genes.Maf_LT_0.01_EPACTS <- LOGISTIC.t[LOGISTIC.t$P  < `^`(10,-6),]
Genes.Maf_LT_0.01_EPACTS$SNP
write.table(Genes.Maf_GT_0.01_EPACTS, "TopGenes_with_pValue_below_1e-6_EPACTS_rare_variants.csv", sep =",", row.names = FALSE, col.names = T, quote = F)


LOGISTIC.t$SNP[LOGISTIC.t$P  > `^`(10,-8)] <- " "

#### Only for MAF > 0.01
TopGeneOrdered <- LOGISTIC.t[(order(LOGISTIC.t$P)),]
head(TopGeneOrdered)
TopGeneOrdered$OR <- exp(TopGeneOrdered$BETA)
TopGeneTable <- TopGeneOrdered[1:20, c("Variant_full", "MAF", "Var_Type", "OR", "P")]
# TopGeneTable <- TopGeneOrdered[1:20, c("SNP", "MAF", "BETA", "Variant", "Var_Type", "P", "Effect", "Variant_full")]
TopGeneTable$P <- format(TopGeneTable$P, scientific=TRUE)
write.table(TopGeneTable, "TopGene_with_pValue_common_variants_MAF_GT_0.01.csv", sep =",", row.names = FALSE, col.names = T, quote = F)
####


# tiff("Manhattan_plot_EPACTs_Complete_stats_single_variant_max_maf_0.01.tiff", height = 20, width = 30, units='cm', compression = "lzw", res = 300)

### Min MAF 0
jpeg("Manhattan_plot_EPACTs_Complete_stats_single_variant_min_maf_0_SEX_PC1_PC2_PC3_without_HAPMAP.jpeg", height = 20, width = 30, units='cm', res = 300)
jpeg("Manhattan_plot_EPACTs_Complete_stats_single_variant_min_maf_0.01_SEX_PC1_PC2_PC3_without_HAPMAP.jpeg", height = 20, width = 30, units='cm', res = 300)
### Max MAF 0.01
jpeg("Manhattan_plot_EPACTs_Complete_stats_single_variant_max_maf_0.01_SEX_PC1_PC2_PC3_without_HAPMAP.jpeg", height = 20, width = 30, units='cm', res = 300)
### MMSkat Max MAF 0.01
jpeg("Manhattan_plot_EPACTs_Complete_stats_single_variant_MMSkat_max_maf_0.01_SEX_PC1_PC2_PC3_without_HAPMAP.jpeg", height = 20, width = 30, units='cm', res = 300)



# manhattan(LOGISTIC, main = "", ylim=c(0,15), cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = -log10(1e-06), genomewideline = -log10(1e-08), annotatePval = -log10(1e-08), chrlabs = c(1:22, "X", "Y")) 
manhattan(LOGISTIC.t, main = "", ylim=c(0,15), col = c("blue4", "orange3"), suggestiveline = -log10(1e-06), genomewideline = -log10(1e-08), annotateTop = FALSE, annotatePval = -log10(1e-06), chrlabs = as.character(1:22)) 
dev.off()


################################################################################
################################## FarVat ######################################
################################################################################
setwd("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files2/FarVat/")

# https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
library(lattice)
source("https://raw.githubusercontent.com/achalneupane/rcodes/master/qqrunif_with_lambda.r")
library (qqman)
# LOGISTIC<-read.table("EPACTs_Complete_stats_single_variant_min_maf_0.txt",  comment.char = '&', header = TRUE, sep ="\t", check.names = F)


## 1.  MAX MAF 0.01
LOGISTIC <- fread("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1_post_QC2_with_STATUS.BLUP.RESULTS_MAF1_percent.THEO.gene.res", header = TRUE, sep ="\t", check.names = F)
## Without SEX
# LOGISTIC <- fread("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1_post_QC2_with_STATUS.BLUP.RESULTS_MAF1_percent.THEO.gene.res", header = TRUE, sep ="\t", check.names = F)



## 2. CADD > 20
LOGISTIC <- fread("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1_post_QC2_with_STATUS.BLUP.RESULTS_CADD_above_20.THEO.gene.res", header = TRUE, sep ="\t", check.names = F)


LOGISTIC <- as.data.frame(LOGISTIC)
# LOGISTIC <- LOGISTIC[LOGISTIC$MAF <= 0.01,]

head(LOGISTIC)
LOGISTIC$BP <- LOGISTIC$START
## with SKAT_O pvalues
LOGISTIC$PVALUE <- LOGISTIC$P_SKATO
## With CLP p-values
# LOGISTIC$PVALUE <- LOGISTIC$P_CLP

LOGISTIC$SNP <- LOGISTIC$GENE

LOGISTIC <- LOGISTIC[!is.na(LOGISTIC$PVALUE),]
LOGISTIC$PVALUE <- as.numeric(LOGISTIC$PVALUE)

my.pvalues <- LOGISTIC$PVALUE

LAMBDA <- round(inflation(my.pvalues), digits = 3)
LAMBDA

### MAX MAF 0.01
jpeg("QQ_plot_FarVat_gene-based_max_maf_0.01_SEX_PC1_PC2_PC3.jpeg", height = 20, width = 15, units='cm', res = 300)
# jpeg("QQ_plot_FarVat_gene-based_max_maf_0.01_SEX_PC1_PC2_PC3_P_CLP.jpeg", height = 20, width = 15, units='cm', res = 300)
### CADD > 20
jpeg("QQ_plot_FarVat_gene-based_CADD_above_20_SEX_PC1_PC2_PC3.jpeg", height = 20, width = 15, units='cm', res = 300)
# jpeg("QQ_plot_FarVat_gene-based_CADD_above_20_SEX_PC1_PC2_PC3_P_CLP.jpeg", height = 20, width = 15, units='cm', res = 300)

qqunif.plot(my.pvalues, LAMBDA = LAMBDA)
dev.off()



sum(my.pvalues < `^`(10,-6))




LOGISTIC$CHR[grepl("X", LOGISTIC$CHR) ] <- 23
LOGISTIC$CHR[grepl("Y", LOGISTIC$CHR) ] <- 24
LOGISTIC$CHR <- as.numeric(LOGISTIC$CHR)

LOGISTIC$SNP <- as.character(LOGISTIC$SNP)

LOGISTIC$VarPos <- paste(LOGISTIC$CHR,LOGISTIC$START, sep = ":")


# # I will remove CHRX and Y completely
LOGISTIC.t <- LOGISTIC[!grepl("23|24",LOGISTIC$CHR ),]
## SKAT P value
LOGISTIC.t$P <- LOGISTIC.t$P_SKATO
# ## CLP P values
# LOGISTIC.t$P <- LOGISTIC.t$P_CLP

## MAF 1%
# Genes
common.genes_with_Epacts_rare_variantTest.and.farvat.maf0.01 <- LOGISTIC.t[na.omit(match(Genes.Maf_LT_0.01_EPACTS$Gene , LOGISTIC.t$GENE)),]
# Variants
# common.Variants_with_Epacts_rare_variantTest.and.farvat.maf0.01 <- LOGISTIC.t[na.omit(match(gsub("::","",gsub("[a-zA-Z ]", "", Genes.Maf_LT_0.01_EPACTS$Variant_full)) , LOGISTIC.t$VarPos)),]
write.table(common.genes_with_Epacts_rare_variantTest.and.farvat.maf0.01, "Genes_significant_at_p_10-6_in_EPACTs_rareVar_test_also_common_in_farvat_test_maf_0.01.csv", sep =",", row.names = FALSE, col.names = T, quote = F)

## CADD > 20
common.genes_with_Epacts_rare_variantTest.and.farvat.CADD20 <- LOGISTIC.t[na.omit(match(Genes.Maf_LT_0.01_EPACTS$Gene , LOGISTIC.t$GENE)),]
write.table(common.genes_with_Epacts_rare_variantTest.and.farvat.CADD20, "Genes_significant_at_p_10-6_in_EPACTs_rareVar_test_also_common_in_farvat_test_CADD20.csv", sep =",", row.names = FALSE, col.names = T, quote = F)


LOGISTIC.t$SNP[LOGISTIC.t$P  > `^`(10,-4)] <- " "

### Max MAF 0.01
jpeg("Manhattan_plot_FarVat_gene-based_max_maf_0.01_SEX_PC1_PC2_PC3_without_HAPMAP.jpeg", height = 20, width = 30, units='cm', res = 300)
# jpeg("Manhattan_plot_FarVat_gene-based_max_maf_0.01_SEX_PC1_PC2_PC3_without_HAPMAP_CLP_P_values.jpeg", height = 20, width = 30, units='cm', res = 300)
### CADD > 20
jpeg("Manhattan_plot_FarVat_gene-based_CADD_above_20_SEX_PC1_PC2_PC3_without_HAPMAP.jpeg", height = 20, width = 30, units='cm', res = 300)
# jpeg("Manhattan_plot_FarVat_gene-based_CADD_above_20_SEX_PC1_PC2_PC3_without_HAPMAP_CLP_P_values.jpeg", height = 20, width = 30, units='cm', res = 300)



# manhattan(LOGISTIC, main = "", ylim=c(0,15), cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = -log10(1e-06), genomewideline = -log10(1e-08), annotatePval = -log10(1e-08), chrlabs = c(1:22, "X", "Y")) 
manhattan(LOGISTIC.t, main = "", ylim=c(0,15), col = c("blue4", "orange3"), suggestiveline = -log10(1e-06), genomewideline = -log10(1e-08), annotateTop = T, annotatePval = -log10(1e-04), chrlabs = as.character(1:22)) 
dev.off()





















################################
# HAPMAP
setwd("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files4")
PCA <- read.table("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1_post_QC_with_STATUS-HAPMAP-MERGED3-for_PCA.eigenvec", header =T, stringsAsFactors=FALSE)
HAPMAP.ethnicty <- read.table("relationships_w_pops_121708.txt", header = T )
head(HAPMAP.ethnicty)

PCA$COHORT <- "FASe"
PCA$COHORT <- HAPMAP.ethnicty$population[match(PCA$IID, HAPMAP.ethnicty$IID)]
PCA <- PCA[c(1:4,23)]
PCA$COHORT <- as.character(PCA$COHORT)
PCA$COHORT[is.na(PCA$COHORT)] <- "FASe"
write.table(PCA, "FASe_3894_WXS_SNPS_INDELS-hwe-geno0.05-mind0.1-HAPMAP-MERGED3-for_PCA.eigenvec-PC1-PC2-COHORT.txt", sep ="\t", col.names = T, quote = F)
####################
# After removing NHW
####################
# PCA <- read.table("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1_post_QC-HAPMAP-MERGED3-for_PCA.eigenvec", header =T, stringsAsFactors=FALSE)
PCA <- read.table("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic_post_QC_with_STATUS-hwe-geno0.05-mind0.1_post_QC-HAPMAP-MERGED3-for_PCA.eigenvec",  header =T, stringsAsFactors=FALSE)
HAPMAP.ethnicty <- read.table("relationships_w_pops_121708.txt", header = T )
head(HAPMAP.ethnicty)

PCA$COHORT <- "FASe"
PCA$COHORT <- HAPMAP.ethnicty$population[match(PCA$IID, HAPMAP.ethnicty$IID)]
PCA <- PCA[c(1:4,23)]
PCA$COHORT <- as.character(PCA$COHORT)
PCA$COHORT[is.na(PCA$COHORT)] <- "FASe"
write.table(PCA, "FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic_post_QC_with_STATUS-hwe-geno0.05-mind0.1_post_QC-HAPMAP-MERGED3-for_PCA.eigenvec-PC1-PC2-COHORT-After-removing-NHW.txt", sep ="\t", col.names = T, quote = F)



##########PCA





library(ggplot2)
#read data:
# PCAs<-read.table("aurora_7234_WXS_SNPS_INDELS-hwe-geno0.05-mind0.1-HAPMAP-MERGED3-for_PCA.eigenvec-PC1-PC2-COHORT.txt", header=T)
# Ethnicity map ftp://ftp.ncbi.nlm.nih.gov/hapmap/phase_3/relationships_w_pops_121708.txt
# PCAs<- read.table("FASe_3894_WXS_SNPS_INDELS-hwe-geno0.05-mind0.1-HAPMAP-MERGED3-for_PCA.eigenvec", header=T)
### Check R code to create this file for PCA
PCAs<- read.table("FASe_3894_WXS_SNPS_INDELS-hwe-geno0.05-mind0.1-HAPMAP-MERGED3-for_PCA.eigenvec-PC1-PC2-COHORT-After-removing-NHW.txt", header=T)
##plotting:
p <- ggplot(PCAs, aes(x=PC1, y=PC2, color=COHORT)) + geom_point() + xlab("PC1") + ylab("PC2") + ggtitle("Aquilla 3891") +
  scale_color_manual(values = c('green', 'black', 'red', "blue"))  
ggsave("Aquilla-FASe-PCs-COHORT_After_Removing_NHW_with_HAPMAP.jpg", p, plot = last_plot(), device = NULL, scale = 1, width = 16, height = 9, dpi = 300, limitsize = TRUE)


######################### 
df <- read.table("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files2/FarVat/phenotype1_with_PC1-PC4_with_1502_samples_and_no_hapmap.ped", header = T)
dim(df)
head(df)
WXS_Percent <- bias.df[bias.df$IID %in% df$IID,]
(table(as.character(WXS_Percent$WXS))/sum(table(as.character(WXS_Percent$WXS))))*100
# 1 (WES)       2 (WGS)
# 72.70306 27.29694

aurora_Age <- read.table("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/00-Aurora-IDlists_Age_list.csv", header = T, sep = ",")
dim(aurora_Age)
head(aurora_Age)

sum(aurora_Age$IID %in% df$IID)

df2 <- merge(aurora_Age, df, by = "IID", all =TRUE)
dim(df2)

df3 <- df2[which(df2$IID %in% df$IID),]
dim(df3)
sum(is.na(df3$AAO))
# df.matched <- df3[!is.na(df3$AAO) & df3$STATUS.x ,]

# View(df.matched)
Not.Matched <- df [(which(!df3$IID %in% df.matched$IID_confirm)),]

# AAO for Case(2); ALA for Control (1)
df.not.matched1 <- df3[is.na(df3$AAO) & df3$STATUS.y==2,]
df.not.matched2 <- df3[is.na(df3$ALA) & df3$STATUS.y==1,]

age_not_available <- rbind(df.not.matched1, df.not.matched2)

unavailable <- (unique(c(as.character(Not.Matched$IID), as.character(age_not_available$IID))))

Age.Not.Found <- df[df$IID %in% unavailable,]

colnames(Age.Not.Found)[grep("confirm", colnames(Age.Not.Found))] <- "Recoded_IID"

# From SHERMAN object from Recode_WU_IDs.r
Age.Not.Found$Recoded_IID <- as.character(Age.Not.Found$Recoded_IID)
keep.id <- Age.Not.Found$Recoded_IID

Age.Not.Found$Recoded_IID <- SHERMAN$recoded_ID[match(as.character(Age.Not.Found$IID), as.character(SHERMAN$ID) )]
Age.Not.Found$Recoded_IID <- as.character(Age.Not.Found$Recoded_IID)
Age.Not.Found$Recoded_IID[is.na(Age.Not.Found$Recoded_IID)] <- keep.id [is.na(Age.Not.Found$Recoded_IID)]


Age.not.found2 <- df3[df3$AAO==0 & df3$ALA==0,]
IIDs.for.age.not.found <- unique(as.character(Age.not.found2$IID))
IIDs.for.age.not.found.df <- df[match(IIDs.for.age.not.found, df$IID_confirm),]

head(IIDs.for.age.not.found.df)
head(Age.Not.Found)
colnames(IIDs.for.age.not.found.df)[7] <- "Recoded_IID"
Age.Not.Found <- rbind(Age.Not.Found, IIDs.for.age.not.found.df)

write.table(Age.Not.Found, "FASe_353_samples_with_no_age.csv", sep =",", col.names = T, quote = F, row.names = FALSE)



# Confirm case control status

df <- read.table("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files2/FarVat/phenotype1_with_PC1-PC4_with_1502_samples_and_no_hapmap.ped", header = T)
df_silva <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/03-silva_201703/04-FASe-complete/silva_FASE_complete_1290ID.PHENO", header = T)


sum(df$IID%in% df_silva$vcfID)
df2 <- cbind(df,df_silva[match(df$IID, df_silva$vcfID),])
# > sum(grepl ("-9", df2$STATUS))
# [1] 147      

mismatch <- df2[which(df2$STATUS != df2$PHENOTYPE),]
sum(grepl("-9", mismatch$STATUS))

table(mismatch$PHENOTYPE)
# 1   2 
# 125  33
table(mismatch$STATUS)
# -9  2 
# 87 71

table(table(final.df$FID))



########################################

df <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/20200518-201909_MGI_gDNA_LINDSEY/Sex_OK_List_Kristy.txt", header = TRUE)
head(df)
dim(df)


coverage <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/20200518-201909_MGI_gDNA_LINDSEY/20200518-201909_MGI_gDNA_LINDSEY.data.transfer.629.csv", header = TRUE, sep =",")
dim(coverage)
# View(coverage)

# coverage$MEAN_COVERAGE
coverage$DNA <- as.character(coverage$DNA)
df$IID <- as.character(df$IID)
df$Mean_Coverage <- "NULL"


foo = function(x, n, i){
  do.call(c, lapply(x, function(X)
    paste(unlist(strsplit(X, "-"))[(n+1):(i)], collapse = "-")))
}

#USAGE
coverage$IID <- foo(x = coverage$DNA, n = 1, i = 2)

df$Mean_Coverage <-  coverage$MEAN_COVERAGE [match(df$IID, coverage$IID)]

df <- df[!is.na(df$Mean_Coverage),]

df <- df[!grepl("MAP_62589|MAP_64341|MAP_68128|MAP_63386|MAP_64559|MAP_64359|MAP_22392|MAP_23392|MAP_68128|MAP_64559", df$IID),]

write.table(df, "/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/20200518-201909_MGI_gDNA_LINDSEY/20200518-201909_MGI_gDNA_LINDSEY_Sex_OK_List_with_coverage.csv", sep =",", col.names = T, quote = F, row.names = FALSE)


################ All Aquilla
df <- read.csv("/40/AD/AD_Seq_Data/03.-phenotype/2021-01-Aquilla-phenotype/Joint_Call_Round1/Aquilla_9576_phenotype.csv", header = TRUE, sep = ",", check.names = F)
head(df)

# Case=2; control=1
table(df$STATUS)

aquilla <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files_all_Aquilla/AQUILLA_WXS_SNPS_INDELS_picard_biallelic.fam",  header = FALSE, sep = " ", check.names = F)
head(aquilla)
aquilla$V2 %in% df$`Aquilla Renamed ID (7983)`

replicates.aquilla <- aquilla[grep ("_R1$|_R2$", aquilla$V2),]
dim(replicates.aquilla)
sum(aquilla$V2 %in% df$`Aquilla vcfID`)

colnames(df)[1] <- "AquillaVCFID" 
aquilla_merged <- merge(aquilla, df, by.x = "V2", by.y = "AquillaVCFID")
dim(aquilla_merged)

######################################################## 03-05-2021; After Vicky gave me the list for FASe with correct status
# Here, I am conitnuing from "Method 2" on the shell script 01-Aquilla_derived_from_Aurora-preQC_3_05_21.sh

final.df <- read.csv("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/Working-FASe2_from_Vicky.csv", header = TRUE, sep =",", check.names = FALSE)
table(table(final.df$FID))

as.character(final.df$`Aquilla vcfID`[duplicated(final.df$`Aquilla vcfID`)])

colnames(final.df)
final.df <- final.df[,c("Aquilla vcfID", "Seq round / Project", "WXS", "FID", "IID (Best ID)", "PID", "MID", "Sex", "STATUS", "APOE", "APOE4", "AAO", "ALA", "AGE", "ADCO")]
colnames(final.df) <- paste0("Vicky_", c("vcfID", "Project", "WXS", "FID", "IID (Best ID)", "PID", "MID", "Sex", "STATUS", "APOE", "APOE4", "AAO", "ALA", "AGE", "ADCO"))
final.df$Vicky_vcfID <- as.character(final.df$Vicky_vcfID)

all.fam <- read.csv("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files2/FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic.fam", header = FALSE, sep =" ", check.names = FALSE)
dim(all.fam)
head(all.fam)
colnames(all.fam) <- c("FID_FAM", "IID_FAM", "MID_FAM", "PID_FAM", "SEX_FAM", "STATUS_FAM")

sum(final.df$Vicky_vcfID %in% all.fam$IID)

R1_R2 <- all.fam[grep("_R1$|_R2$", all.fam$IID),]
NON_R1_R2 <- all.fam[!grepl("_R1$|_R2$", all.fam$IID),]


WANT1 <- {}
for (i in 1:nrow(NON_R1_R2)){
  wanted.row <- final.df[grep(paste0("^",gsub("_R1|_R2", "", R1_R2$IID[i]),"$"), final.df$Vicky_vcfID),c(1,4,6:15)]
  if(length(wanted.row$Vicky_vcfID)!=0){
    tmpWANT <- cbind(R1_R2[i,], wanted.row)
    WANT1 <- rbind(WANT1, tmpWANT)
    # final.df$vcfID[grep(paste0(paste0("^",all.fam$V2[i],"$"),"|",paste0("^",all.fam$V2[i],"_R1$")), final.df$vcfID)]  
  }
}




WANT2 <- {}
for (i in 1:nrow(NON_R1_R2)){
  wanted.row <- final.df[grep(paste0("^",NON_R1_R2$IID[i],"$"), final.df$Vicky_vcfID),c(1,4,6:15)]
  if(length(wanted.row$Vicky_vcfID)!=0){
  tmpWANT <- cbind(NON_R1_R2[i,], wanted.row)
  WANT2 <- rbind(WANT2, tmpWANT)
  # final.df$vcfID[grep(paste0(paste0("^",all.fam$V2[i],"$"),"|",paste0("^",all.fam$V2[i],"_R1$")), final.df$vcfID)]  
}
}


WANT <- rbind(WANT1, WANT2)

dim(WANT)
WANT <- WANT[-which(duplicated(WANT)), ]
dim(WANT)

wantID <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/FASe_ID_list_round1.txt", header = FALSE, sep =",", check.names = FALSE)
colnames(wantID) <- c("Achal_IID", "Achal_VCFID", "Achal_PR")

WANT.final <- cbind(WANT, wantID[match(WANT$IID_FAM, wantID$Achal_VCFID),])


setwd("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files3")

PHENO <- read.table("FASe_3894_Phenoscope.txt", header = T)
dim(PHENO)
WANT.FINAL <- cbind(WANT.final, PHENO[match(WANT.final$IID_FAM, PHENO$IID),-c(3:4)])

colnames(WANT.FINAL)

WANT.FINAL <- WANT.FINAL[,c( "Vicky_FID", "IID_FAM", "FID_FAM", "Vicky_PID", "Vicky_MID", "Vicky_Sex", "Vicky_STATUS", 
               "Vicky_APOE", "Vicky_APOE4", "Vicky_AAO", "Vicky_ALA", "Vicky_AGE", "Vicky_ADCO", "Project",
               "WXS", "Macrogen_WGS", "MAPT_A152T", "LOAD_WES", "Genentech_WGS", "Genentech_WES", "Broad_WGS",
               "TGI_WES", "MGI_FASeEOAD_201605", "Otogenetics_WES", "phs000572_201508", "phs000572_201707",
               "phs000572_201802", "phs000572_201612", "X201907_USUHS_gDNA_SHERMAN")]

colnames(WANT.FINAL) <- gsub("Vicky_", "", colnames(WANT.FINAL))

table(table(WANT.FINAL$FID))
sort(table(WANT.FINAL$FID), decreasing = T) [1:200]

# Recode APOE (See Kristy's email, Dated 3/5/2021: "Here is the chart from our APOE protocol")

table(WANT.FINAL$APOE)
WANT.FINAL$FID <-  as.character(WANT.FINAL$FID)
colnames(WANT.FINAL)[grep("X201907_USUHS_gDNA_SHERMAN", colnames(WANT.FINAL))] <- "201907_USUHS_gDNA_SHERMAN"

PHENO<- WANT.FINAL[-3]
colnames(PHENO)[2] <- "IID"

# CA=1280; CO=341 (Vicky)
write.table(PHENO, "PHENO.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)

write.table(WANT.FINAL[c(3,2)], "Samples_to_keep_for_Aquilla_FASe_1784.txt", sep ="\t", col.names = F, quote = F, row.names = FALSE)

NewIID <- WANT.FINAL[c(1,2)]

OldIID <- read.table("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files3/FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic_post_QC.fam", header = F, sep = " ")
OldIID <- OldIID[1:2]
FID_IID_recode <- cbind(OldIID, NewIID[match(OldIID$V2, NewIID$IID_FAM),])
colnames(FID_IID_recode) <- c("oldFID", "oldIID", "newFID", "newIID")

write.table(FID_IID_recode, "FID_IID_recode_vicky.txt", sep =" ", col.names = T, quote = F, row.names = FALSE)

# Update SEX
SEX <- WANT.FINAL[c(1:2,6)]
write.table(SEX, "SEX_recode_vicky.txt", sep =" ", col.names = T, quote = F, row.names = FALSE)

# Update STATUS
STATUS <- WANT.FINAL[c(1:2,7)]
write.table(STATUS, "update_pheno.txt", sep =" ", col.names = T, quote = F, row.names = FALSE)


########## After final QC

FINAL.DF.1608 <- PHENO[match(FinalWantedSamples$IID, PHENO$IID),]
write.table(FINAL.DF.1608[1:2], "Extract_1608_after_QC_vicky.txt", sep =" ", col.names = T, quote = F, row.names = FALSE)

## This is for replication study
rep.study <- FINAL.DF.1608[!grepl("phs" , FINAL.DF.1608$Project),]
table(table(rep.study$FID))


df <- read.table("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic_post_QC_with_STATUS-hwe-geno0.05-mind0.1_post_QC2.fam", header = F)
dim(df)

df.merged <- merge(FINAL.DF.1608, df , by.x = "IID", by.y = "V2")
sum(df.merged$STATUS==df.merged$V6)
as.character(df.merged$Sex)==as.character(df.merged$V5)

### PCA without HAPMAP
PCA_final <- read.table("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic_post_QC_with_STATUS-hwe-geno0.05-mind0.1_post_QC2-PCAS.eigenvec",  header =T, stringsAsFactors=FALSE)
PCA_final <- PCA_final[2:6]
FINAL.DF.1608_with_PCA <- merge(FINAL.DF.1608, PCA_final, by.x = "IID", by.y = "IID")
FINAL.DF.1608_with_PCA <- cbind(FINAL.DF.1608_with_PCA[2],FINAL.DF.1608_with_PCA[-2]) 
write.table(FINAL.DF.1608_with_PCA, "Pheno_1608_after_QC_vicky.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)

## After APOE genotype 
dim(FINAL.DF.1608_with_PCA)
# View(FINAL.DF.1608_with_PCA)
# plink --bfile $BFILE --extract APOE_SNPs.list --recode --out APOE_PLINK_OUT

APOEvar1 <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files3/APOE_PLINK_OUT1.ped", header =FALSE, stringsAsFactors=FALSE, check.names = F)
colnames(APOEvar1) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE", "REF", "ALT")
APOEvar1$GENO <- paste(APOEvar1$REF, APOEvar1$ALT, sep ="/")

APOEvar2 <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files3/APOE_PLINK_OUT2.ped", header =FALSE, stringsAsFactors=FALSE, check.names = F)
colnames(APOEvar2) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE", "REF", "ALT")
APOEvar2$GENO <- paste(APOEvar2$REF, APOEvar2$ALT, sep ="/")
APOEVAR <-cbind(FINAL.DF.1608_with_PCA, rs429358 = APOEvar1[match(FINAL.DF.1608_with_PCA$IID, APOEvar1$IID), "GENO"])
APOEVAR <-cbind(APOEVAR, rs7412 = APOEvar2[match(APOEVAR$IID, APOEvar2$IID), "GENO"])

APOEVAR$rs429358 <- as.character(APOEVAR$rs429358)
APOEVAR$rs7412 <- as.character(APOEVAR$rs7412)

APOEVAR$rs429358[APOEVAR$rs429358 == "T/T"] <- 2
APOEVAR$rs429358[APOEVAR$rs429358 == "C/T"] <- 1
APOEVAR$rs429358[APOEVAR$rs429358 == "C/C"] <- 0
APOEVAR$rs429358[APOEVAR$rs429358 == "0/0"] <- -9

APOEVAR$rs7412[APOEVAR$rs7412 == "T/T"] <- 2
APOEVAR$rs7412[APOEVAR$rs7412 == "C/T"] <- 1
APOEVAR$rs7412[APOEVAR$rs7412 == "C/C"] <- 0
APOEVAR$rs7412[APOEVAR$rs7412 == "T/C"] <- 1
APOEVAR$rs7412[APOEVAR$rs7412 == "0/0"] <- -9

# recode PHENO APOE
APOEVAR$APOE[is.na(APOEVAR$APOE)] <- -9
APOEVAR$rs429358PHENO [APOEVAR$APOE == 22]  <- 2
APOEVAR$rs429358PHENO [APOEVAR$APOE == 23 | APOEVAR$APOE == 32]  <- 2
APOEVAR$rs429358PHENO [APOEVAR$APOE == 33]  <- 2
APOEVAR$rs429358PHENO [APOEVAR$APOE == 24 | APOEVAR$APOE == 42]  <- 1
APOEVAR$rs429358PHENO [APOEVAR$APOE == 34 | APOEVAR$APOE == 43]  <- 1
APOEVAR$rs429358PHENO [APOEVAR$APOE == 44]  <- 0
APOEVAR$rs429358PHENO [APOEVAR$APOE == -9]  <- -9


APOEVAR$rs7412PHENO [APOEVAR$APOE == 22]  <- 2
APOEVAR$rs7412PHENO [APOEVAR$APOE == 23 | APOEVAR$APOE == 32]  <- 1
APOEVAR$rs7412PHENO [APOEVAR$APOE == 33]  <- 0
APOEVAR$rs7412PHENO [APOEVAR$APOE == 24 | APOEVAR$APOE == 42]  <- 1
APOEVAR$rs7412PHENO [APOEVAR$APOE == 34 | APOEVAR$APOE == 43]  <- 0
APOEVAR$rs7412PHENO [APOEVAR$APOE == 44]  <- 0
APOEVAR$rs7412PHENO [APOEVAR$APOE == -9]  <- -9

sum(APOEVAR$rs429358PHENO == APOEVAR$rs429358)
# [1] 1529
sum(APOEVAR$rs7412PHENO == APOEVAR$rs7412)
# [1] 1491

table(rs429358PHENO=APOEVAR$rs429358PHENO, rs429358=APOEVAR$rs429358)
FREQrs429358 <- as.data.frame(table(rs429358PHENO=APOEVAR$rs429358PHENO, rs429358=APOEVAR$rs429358, STATUS=APOEVAR$STATUS))
table(rs7412PHENO=APOEVAR$rs7412PHENO, rs7412=APOEVAR$rs7412)
FREQrs7412 <- as.data.frame(table(rs7412PHENO=APOEVAR$rs7412PHENO, rs7412=APOEVAR$rs7412, STATUS=APOEVAR$STATUS))
# table(STATUS=APOEVAR$STATUS, rs429358=APOEVAR$rs429358)
# COUNTS_rs429358 <- as.data.frame(table(STATUS=APOEVAR$STATUS, rs429358=APOEVAR$rs429358))
# table(STATUS=APOEVAR$STATUS, rs429358=APOEVAR$rs429358)
# COUNTS_rs429358 <- as.data.frame(table(STATUS=APOEVAR$STATUS, rs429358=APOEVAR$rs429358))

# # recode PHENO APOE
# APOEVAR$APOE[is.na(APOEVAR$APOE)] <- -9
# APOEVAR$rs429358PHENO [APOEVAR$APOE == 22]  <- 2
# APOEVAR$rs429358PHENO [APOEVAR$APOE == 23]  <- 2
# APOEVAR$rs429358PHENO [APOEVAR$APOE == 32]  <- 1
# APOEVAR$rs429358PHENO [APOEVAR$APOE == 33]  <- 2
# APOEVAR$rs429358PHENO [APOEVAR$APOE == 24]  <- 1
# APOEVAR$rs429358PHENO [APOEVAR$APOE == 42]  <- 1
# APOEVAR$rs429358PHENO [APOEVAR$APOE == 34]  <- 1
# APOEVAR$rs429358PHENO [APOEVAR$APOE == 43]  <- 0
# APOEVAR$rs429358PHENO [APOEVAR$APOE == 44]  <- 0
# APOEVAR$rs429358PHENO [APOEVAR$APOE == -9]  <- -9
# 
# 
# APOEVAR$APOE[is.na(APOEVAR$APOE)] <- -9
# APOEVAR$rs7412PHENO [APOEVAR$APOE == 22]  <- 2
# APOEVAR$rs7412PHENO [APOEVAR$APOE == 23]  <- 1
# APOEVAR$rs7412PHENO [APOEVAR$APOE == 32]  <- 2
# APOEVAR$rs7412PHENO [APOEVAR$APOE == 33]  <- 0
# APOEVAR$rs7412PHENO [APOEVAR$APOE == 24]  <- 1
# APOEVAR$rs7412PHENO [APOEVAR$APOE == 42]  <- 1
# APOEVAR$rs7412PHENO [APOEVAR$APOE == 34]  <- 0
# APOEVAR$rs7412PHENO [APOEVAR$APOE == 43]  <- 1
# APOEVAR$rs7412PHENO [APOEVAR$APOE == 44]  <- 0
# APOEVAR$rs7412PHENO [APOEVAR$APOE == -9]  <- -9
# 
# # sum(APOEVAR$rs429358PHENO == APOEVAR$rs429358)
# # [1] 745
# # sum(APOEVAR$rs7412PHENO == APOEVAR$rs7412)
# # [1] 720







# APOEVAR$rs429358[as.character(APOEVAR$rs429358) == "0/0"] <- -9
# APOEVAR$rs7412[as.character(APOEVAR$rs7412) == "0/0"] <- -9

# APOEVAR$APOE_CODE <- paste(APOEVAR$rs429358, APOEVAR$rs7412, sep = "_")

# APOEVAR$APOE_GENO <- -9
# APOEVAR$APOE[is.na(APOEVAR$APOE)] <- -9
# 
# APOEVAR$APOE_GENO [APOEVAR$APOE_CODE == "T/T_T/T"] <- 22
# APOEVAR$APOE_GENO [APOEVAR$APOE_CODE == "T/T_C/T"] <- 23
# APOEVAR$APOE_GENO [APOEVAR$APOE_CODE == "T/T_C/C"] <- 33
# APOEVAR$APOE_GENO [APOEVAR$APOE_CODE == "C/T_C/T"] <- 24
# APOEVAR$APOE_GENO [APOEVAR$APOE_CODE == "C/T_C/C"] <- 34
# APOEVAR$APOE_GENO [APOEVAR$APOE_CODE == "C/C_C/C"] <- 44

# sum(APOEVAR$APOE == APOEVAR$APOE_GENO)
# # tt <- APOEVAR[APOEVAR$APOE != APOEVAR$APOE_GENO,]

colnames(APOEVAR)[1:6] <- c("##FAM_ID", "IND_ID", "FAT_ID", "MOT_ID", "SEX", "STATUS")

## Missing values are not accpetable in Covariate file, so NA = -9
APOEVAR$SEX[is.na(APOEVAR$SEX)] <- -9

table(STATUS=APOEVAR$STATUS, APOE=APOEVAR$rs429358)

# APOEVAR$APOEmatch <- ifelse (APOEVAR$APOE == APOEVAR$APOE_GENO, "Match", "Mis-match")
write.table(APOEVAR, "/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files3/epacts/Pheno_1608_after_QC_vicky_with_APOE_Genotype_check.ped", sep =" ", col.names = T, quote = F, row.names = FALSE)
write.table(APOEVAR[c(1:6, 29:32, 33,34)], "/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files3/epacts/Pheno_1608_after_QC_vicky_with_APOE_Genotype_check_for_epacts.ped", sep =" ", col.names = T, quote = F, row.names = FALSE)

####
setwd("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files4")
PCA <- read.table("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1_post_QC_with_STATUS_post_QC2-HAPMAP-MERGED3-for_PCA.eigenvec", header =T, stringsAsFactors=FALSE)
HAPMAP.ethnicty <- read.table("relationships_w_pops_121708.txt", header = T )
head(HAPMAP.ethnicty)

PCA$COHORT <- "FASe"
PCA$COHORT <- HAPMAP.ethnicty$population[match(PCA$IID, HAPMAP.ethnicty$IID)]
PCA <- PCA[c(1:4,23)]
PCA$COHORT <- as.character(PCA$COHORT)
PCA$COHORT[is.na(PCA$COHORT)] <- "FASe"
write.table(PCA, "FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic_post_QC_with_STATUS-hwe-geno0.05-mind0.1_post_QC2-HAPMAP-MERGED3-for_PCA.eigenvec-PC1-PC2-COHORT.txt", sep ="\t", col.names = T, quote = F)

###############################


# ## Plink QC
# # https://meyer-lab-cshl.github.io/plinkQC/articles/plinkQC.html
# install.packages("plinkQC")
# package.dir <- find.package('plinkQC')
# install.packages("https://cran.r-project.org/src/contrib/plinkQC_0.3.3.tar.gz", repos=NULL, type="source", ask=FALSE)
# library("plink")
# 
# library(devtools)
# install_github("meyer-lab-cshl/plinkQC")
# 
# getOption("repos")
# options(repos = c(CRAN = "https://cran.rstudio.org"))
# install.packages("plinkQC", dependencies=TRUE, repos='https://cran.r-project.org/src/contrib/plinkQC_0.3.3.tar.gz')



