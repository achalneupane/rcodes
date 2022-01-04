setwd("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/05-FTD-analysis/")

## How many variants there are in each gene
ANNO <- read.table("ALL_FTD_gene_variants.tsv", header = F, stringsAsFactors = F)
library(dplyr) 
unique_VARS <- distinct_at(ANNO, vars(V3, V9))
table(unique_VARS$V9)
# LDHA    LDHB   SARM1 SLC16A1 SLC16A3 SLC16A4 SLC16A7 
# 76      54     330      47     200     112      66 


## How many non-synonymous
ANNO <- read.table("FTD_gene_variants.tsv", header = F, stringsAsFactors = F)
# Table for all Non-synonymous
TABLE <- cbind.data.frame(ANNO$V3, ANNO$V9, ANNO$V7, ANNO$V14, ANNO$V17, ANNO$V15)
colnames(TABLE) <- c("CHR:BP:REF:ALT", "Gene", "EFFECT", "p.Change", "MAF", "CADD")
TABLE <- distinct_at(TABLE, vars(`CHR:BP:REF:ALT`, p.Change), .keep_all = T)

write.table(TABLE, "Non-synonymous_variants_by_p.change.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)




# Get unique variants counts
ANNO <- ANNO %>% distinct(V3, .keep_all = TRUE)
# These are the unique counts of variants
table(ANNO$V9)
# LDHA    LDHB   SARM1 SLC16A1 SLC16A3 SLC16A4 SLC16A7 
# 15      16      26      19      37      37      37

## How many non-synonymous with MAF 1%
MAF <- ANNO[ANNO$V17 <= 0.01,]
table(MAF$V9)
# LDHA    LDHB   SARM1 SLC16A1 SLC16A3 SLC16A4 SLC16A7 
# 6       8      11       8      11      15      16 
sum(table(MAF$V9))
# 75

## How many non-synonymous with CADD 20
CADD <- ANNO[ANNO$V15 >= 20,]
table(CADD$V9)
# LDHA    LDHB   SARM1 SLC16A1 SLC16A3 SLC16A4 SLC16A7 
# 15      12      23      14      30      23      23
sum(table(CADD$V9))
# 140




TYPE <- c("LDHA", "LDHB", "SARM1", "SLC16A1", "SLC16A3", "SLC16A4", "SLC16A7")
ALL_CARRIERS <- {}
Carrier.samples <- {}
ALL_VARIANTS <- list()
for (i in 1:length(TYPE)){
df <- read.table(paste0(TYPE[i],"_variants_in_Aquilla.ped"), stringsAsFactors = FALSE, header = T, colClasses = c("character"))
# df <- read.table("FTD_variants_in_Aquilla.ped", stringsAsFactors = FALSE, header = T, colClasses = c("character"))

vars <- read.table(paste0(TYPE[i], "_SNPs.txt"))
vars <- as.character(t(vars))


# Label our headers with FTD variants 
cols <- paste0(rep(vars, each = 2), c("_A1", "_A2"))
cols <- c("FID", "IID", "PAT", "MAT", "SEX", "STATUS", cols)
colnames(df) <- cols

# subset the data with FTD individuals only
FTDpheno <- read.delim("FTD_phenotype.txt", sep = "\t")

# Only 30 FTD samples available
sum(df$IID %in% FTDpheno$Aquilla.vcfID)
# One sample failed QC: 64405
FTDpheno$Aquilla.vcfID[!FTDpheno$Aquilla.vcfID %in% df$IID]

df <- df [df$IID %in% FTDpheno$Aquilla.vcfID,]
rownames(df) <- df$IID
df <- df[cols[-c(1:6)]]

## If A1 == A2 ---> Homozygous; If A1 != A2 ---> Heterozygous (carrier);

# vars <- unique(sapply(strsplit(colnames(df),"_"), `[`, 1))
# for (j in 1:length(vars)){
#   df[, paste(vars[j], "Status", sep = "_")] <- as.numeric(df[,paste0(vars[j], "_A1")]!=df[,paste0(vars[j], "_A2")])
# }

odds <- seq(1, ncol(df), by = 2)
odds
df <- do.call(cbind, Map(function(z, nm) setNames(cbind(z, Status = +(z[[1]] != z[[2]])), c(names(z), nm)),
                   split.default(df, rep(odds, each = 2)),
                   paste0(gsub("_A1","", colnames(df)[odds]), "_Status"), USE.NAMES = FALSE))

# Get Heterozygosity Status
df <- df[grepl("Status", colnames(df))]
carriers <- sum(df ==1)
names(carriers) <- TYPE[i]
ALL_CARRIERS <- c(ALL_CARRIERS, carriers)

# Also get variants
df.TF <- df ==1
if (sum(df.TF) == 0){
  next
}
cc <- colSums(df.TF)[colSums(df.TF) > 0]
variants <- as.data.frame(as.vector(cc))

# Get the list of carriers
for (n in 1:length(cc)){
Carrier.samples.tmp <- row.names(df.TF)[df.TF[,names(cc[n])]]
Carrier.samples.tmp <- paste(TYPE[i], names(cc[n]), paste(Carrier.samples.tmp, sep="", collapse=","), sep = "\t" )
Carrier.samples <- c(Carrier.samples, Carrier.samples.tmp)
}

colnames(variants) <- "Carriers"
variants$Variant <- gsub("_Status", "", paste(names(cc), TYPE[i], sep = "_"))
ALL_VARIANTS <- as.data.frame(cbind(ALL_VARIANTS, as.matrix(variants)))
}




ALL_CARRIERS
# LDHA    LDHB   SARM1 SLC16A1 SLC16A3 SLC16A4 SLC16A7 
# 0       0       1      18       0       0      12
# Variants and # Carriers
ALL_VARIANTS
# chr17:28396267:C:T_SARM1 chr1:112913924:A:T_SLC16A1 chr12:59779575:A:T_SLC16A7
#      1                         18                              12




# Summary table 1
FTDpheno <- FTDpheno[!FTDpheno$Aquilla.vcfID %in% "64405",]
table(FTDpheno$ADCO)
table(FTDpheno$ADCO, FTDpheno$Sex)
# Females (1 male, 2 female)
table(FTDpheno$ADCO, FTDpheno$Sex)[,2] / (table(FTDpheno$ADCO, FTDpheno$Sex)[,1] + table(FTDpheno$ADCO, FTDpheno$Sex)[,2]) * 100
## % APOE +
FTDpheno$APOE4any <- ifelse(grepl("24|42|34|43|44", FTDpheno$APOE), "+", "-")
table(FTDpheno$ADCO, FTDpheno$APOE4any)[,2] / (table(FTDpheno$ADCO, FTDpheno$APOE4any)[,1] + table(FTDpheno$ADCO, FTDpheno$APOE4any)[,2]) * 100
## Average AAO
aggregate(FTDpheno$AAO~FTDpheno$ADCO, FTDpheno, mean, na.rm =TRUE, na.action = na.pass)

#############################
# Read IBD file
IBD <- fread("Aquilla_7983_WXS_SNPS-INDELS-geno0.05-hwe-mind0.1-WXSmissingCLEAN-sexupdate_cleaned1-IBD-MAF05.genome", select=c("IID1", "IID2", "Z0", "Z1", "Z2", "PI_HAT"))
# IBD <- fread("Aquilla_7983_WXS_SNPS-INDELS-geno0.05-hwe-mind0.1-WXSmissingCLEAN-sexupdate-IBD.genome",  select=c("IID1", "IID2", "Z0", "Z1", "Z2", "PI_HAT"))

parent_offsping <- IBD[(IBD$Z0 < 0.125 & IBD$Z1 > 0.75), ]

parent_offsping$Relationship <- "parent-offspring"

sibPairs <- IBD[(IBD$Z0 < 0.5 &
                   IBD$Z0 > 0.10 &
                   IBD$Z1 < 0.75 &
                   IBD$Z1 > 0.25),] 
sibPairs$Relationship <- "sib-pairs"

duplicates <- IBD[(IBD$Z0 < 0.25 &
                     IBD$Z1 < 0.25), ]
duplicates$Relationship <- "duplicates"

relatives.ALL <- rbind(parent_offsping, sibPairs, duplicates)

FTD_relatives <- relatives.ALL[relatives.ALL$IID1 %in% as.character(FTDpheno$Aquilla.vcfID) & relatives.ALL$IID2 %in% as.character(FTDpheno$Aquilla.vcfID),]

###############################

Carrier.samples
Carrier.samples <- read.table(text = Carrier.samples, sep = "\t", header = FALSE)
colnames(Carrier.samples) <- c("GENE", "variant", "carriers")

TABLE2 <- {}
for (i in 1:nrow(Carrier.samples)){
# how many of the 18 were male? how many were female?
carriers <- as.character(Carrier.samples[i, "carriers"])
carriers <- unlist(strsplit(carriers, split=","))

MALE <- sum(FTDpheno[FTDpheno$Aquilla.vcfID %in% carriers, "Sex"] == 1)
FEMALE <- sum(FTDpheno[FTDpheno$Aquilla.vcfID %in% carriers, "Sex"] == 2)

# age ranges? 
AGE.RANGE <- paste(range(FTDpheno[FTDpheno$Aquilla.vcfID %in% carriers, "AAO"] [!FTDpheno[FTDpheno$Aquilla.vcfID %in% carriers, "AAO"]==0], na.rm = T), collapse = "-")
# APOE genotypes?
APOE <- paste(FTDpheno[FTDpheno$Aquilla.vcfID %in% carriers, "APOE"], collapse = ",")

## are they related? 
FTD_relatives$KEY <- paste(FTD_relatives$IID1,FTD_relatives$IID2, sep = ":")
if (length(carriers) < 2){
RELATED.PAIRS <- ""
} else {
LOOKUP <- combn(carriers, 2, FUN=paste0, collapse = ":")
# RELATED.PAIRS <- paste(FTD_relatives[FTD_relatives$KEY %in% LOOKUP, c("KEY", "PI_HAT")], collapse = ";")
RELATED.PAIRS <- paste(paste0(FTD_relatives$KEY[FTD_relatives$KEY %in% LOOKUP], " (", FTD_relatives$PI_HAT[FTD_relatives$KEY %in% LOOKUP],")" ), collapse = ", ")
}

TABLE2.tmp <- cbind.data.frame(Carrier.samples[i,], APOE= APOE, MALE= MALE, FEMALE = FEMALE, AGE_RANGE = AGE.RANGE, RELATED_PAIRS =  RELATED.PAIRS)
TABLE2 <- rbind.data.frame(TABLE2, TABLE2.tmp)
}

write.table(TABLE2, "Table_of_carriers_AGE_SEX_APOE.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)

##########################
## APOE4 +
sum(grepl("24|42|34|43|44", unlist(strsplit(as.character(TABLE2[1,"APOE"]), split=","))))/ length(unlist(strsplit(as.character(TABLE2[1,"APOE"]), split=","))) *100
sum(grepl("24|42|34|43|44", unlist(strsplit(as.character(TABLE2[2,"APOE"]), split=","))))/ length(unlist(strsplit(as.character(TABLE2[2,"APOE"]), split=","))) * 100
sum(grepl("24|42|34|43|44", unlist(strsplit(as.character(TABLE2[3,"APOE"]), split=","))))/ length(unlist(strsplit(as.character(TABLE2[3,"APOE"]), split=","))) * 100


# Total carriers
length(unique(unlist(strsplit(as.character(TABLE2[,"carriers"]), split=","))))

# Related individuals
length(unique(unlist(strsplit(as.character(TABLE2[1,"RELATED_PAIRS"]), split=","))))
length(unique(unlist(strsplit(as.character(TABLE2[2,"RELATED_PAIRS"]), split=","))))
length(unique(unlist(strsplit(as.character(TABLE2[3,"RELATED_PAIRS"]), split=","))))

########################################################################
## Check variants in genome AD
# https://gnomad.broadinstitute.org/gene/ENSG00000137942?dataset=gnomad_r2_1
# # grep chr17:28396267:C:T Aquilla_7983_WXS_SNPS-INDELS-geno0.05-hwe-mind0.1-WXSmissingCLEAN-sexupdate-annot-snpeff-dbnsfp-FIELDS-simple-HIGH_MODERATE.txt
# 17     28396267        chr17:28396267:C:T      C       T       T       missense_variant        MODERATE        SARM1   ENSG00000004139 transcript      ENST00000585482.6       c.2156C>T       p.Ala719Val     13.6    .       1.977E-4        1.978E-4
# 17      28396267        chr17:28396267:C:T      C       T       T       missense_variant        MODERATE        SARM1   ENSG00000004139 transcript      ENST00000578128.5       c.752C>T        p.Ala251Val     13.6    .       1.977E-4        1.978E-4

## Then search ENSG00000004139 GenomeAD
## Check variants in genome AD
# https://gnomad.broadinstitute.org/gene/ENSG00000004139?dataset=gnomad_r3
# see which one is canonical transcript (which in our case is ENST00000585482.6), then get the p.change for the canonical transcript from the annotated file

########################################################################
## Common ones
cc <- unlist(strsplit(as.character(TABLE2[1,"carriers"]), split=","))
cc <- gsub(" ","", cc)
SARM1 <- unique(unlist(strsplit(sapply(strsplit(cc,"\\("), `[`, 1), ":")))

cc <- unlist(strsplit(as.character(TABLE2[2,"carriers"]), split=","))
cc <- gsub(" ","", cc)
SLC16A1 <- unique(unlist(strsplit(sapply(strsplit(cc,"\\("), `[`, 1), ":")))

cc <- unlist(strsplit(as.character(TABLE2[3,"carriers"]), split=","))
cc <- gsub(" ","", cc)
SLC16A4 <- unique(unlist(strsplit(sapply(strsplit(cc,"\\("), `[`, 1), ":")))

cc <- unlist(strsplit(as.character(TABLE2[4,"carriers"]), split=","))
cc <- gsub(" ","", cc)
SLC16A7 <- unique(unlist(strsplit(sapply(strsplit(cc,"\\("), `[`, 1), ":")))
