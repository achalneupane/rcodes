# convert plink BED to PED
echo "19:44908684:T:C" > APOE_SNPs1.list
echo "19:44908822:C:T" > APOE_SNPs2.list
BFILE="AQUILLA_WXS_SNPS_INDELS_picard_biallelic"
plink --bfile $BFILE --extract APOE_SNPs1.list --recode --out APOE_PLINK_OUT1
plink --bfile $BFILE --extract APOE_SNPs2.list --recode --out APOE_PLINK_OUT2


# plink --bfile AQUILLA_WXS_SNPS_INDELS_picard_biallelic --extract snps.txt --recodeA --out APOE_PLINK_OUT
# snps.txt:
# 19:44908684:T:C
# 19:44908822:C:T


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