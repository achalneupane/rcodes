#############################
## Single variant analysis ##
#############################
# Family (FBAT)
setwd("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files")

NEUP_FAMILIAL_FBAT_ANNO <- read.delim("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/FBAT/FBAT_rare_variant_analysis_results.csv", header = T, sep = "\t", stringsAsFactors = FALSE)
sum(NEUP_FAMILIAL_FBAT_ANNO$P < 0.05, na.rm = T)
colnames(NEUP_FAMILIAL_FBAT_ANNO) <- paste0("NEUP_FAMILIAL_", colnames(NEUP_FAMILIAL_FBAT_ANNO))
# sort
NEUP_FAMILIAL_FBAT_ANNO <- NEUP_FAMILIAL_FBAT_ANNO[order(NEUP_FAMILIAL_FBAT_ANNO$NEUP_FAMILIAL_P),]
NEUP_FAMILIAL_FBAT_ANNO <- NEUP_FAMILIAL_FBAT_ANNO[NEUP_FAMILIAL_FBAT_ANNO$NEUP_FAMILIAL_P < 0.05,]
# write.table(NEUP_FAMILIAL_FBAT_ANNO, "Single_variant_analysis_familial_data_results_p_0.05.csv", sep =",", col.names = T, quote = F, row.names = FALSE)
NEUP_FAMILIAL_FBAT_ANNO$ID <- do.call(paste, c(read.table(text = NEUP_FAMILIAL_FBAT_ANNO$NEUP_FAMILIAL_Marker, sep = ":")[1:2], sep = ":"))
NEUP_FAMILIAL_FBAT_ANNO$geneNAME <- NEUP_FAMILIAL_FBAT_ANNO$NEUP_FAMILIAL_gene


# Replication 
NEUP_UNRELATED_ANNO <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/Firth-Fallback_replication_study_results.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
sum(NEUP_UNRELATED_ANNO$P < 0.05, na.rm = T)
colnames(NEUP_UNRELATED_ANNO) <- paste0("NEUP_UNRELATED_", colnames(NEUP_UNRELATED_ANNO))
NEUP_UNRELATED_ANNO <- NEUP_UNRELATED_ANNO[!is.na(NEUP_UNRELATED_ANNO$NEUP_UNRELATED_P),]
NEUP_UNRELATED_ANNO <- NEUP_UNRELATED_ANNO[order(NEUP_UNRELATED_ANNO$NEUP_UNRELATED_P),]
NEUP_UNRELATED_ANNO <- NEUP_UNRELATED_ANNO[NEUP_UNRELATED_ANNO$NEUP_UNRELATED_P < 0.05,]
# write.table(NEUP_UNRELATED_ANNO, "Single_variant_analysis_unrelated_data_results_p_0.05.csv", sep =",", col.names = T, quote = F, row.names = FALSE)
NEUP_UNRELATED_ANNO$ID <- do.call(paste, c(read.table(text = NEUP_UNRELATED_ANNO$NEUP_UNRELATED_SNP, sep = ":")[1:2], sep = ":"))
NEUP_UNRELATED_ANNO$geneNAME <- NEUP_UNRELATED_ANNO$NEUP_UNRELATED_gene

# Meta-analysis of two datasets
Fixed_Effect_Meta_analysis_result <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/fixed_effect_meta-analysis/Fixed_effect_Meta_analysis_results.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
colnames(Fixed_Effect_Meta_analysis_result) <- paste0("NEUP_META_", colnames(Fixed_Effect_Meta_analysis_result))
Fixed_Effect_Meta_analysis_result <- Fixed_Effect_Meta_analysis_result[order(Fixed_Effect_Meta_analysis_result$NEUP_META_P),]
Fixed_Effect_Meta_analysis_result <- Fixed_Effect_Meta_analysis_result[Fixed_Effect_Meta_analysis_result$NEUP_META_P < 0.05,]
# write.table(Fixed_Effect_Meta_analysis_result, "Single_variant_META_analysis_results_p_0.05.csv", sep =",", col.names = T, quote = F, row.names = FALSE)
Fixed_Effect_Meta_analysis_result$ID <- do.call(paste, c(read.table(text = Fixed_Effect_Meta_analysis_result$NEUP_META_SNP, sep = ":")[1:2], sep = ":"))
Fixed_Effect_Meta_analysis_result$geneNAME <- Fixed_Effect_Meta_analysis_result$NEUP_META_gene


## First filtering paradigm: P< 0.0005 on familial data and P< 0.05 on unrelated data
singlevar_Prokopenko_p0.0005 <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/Rare_variants_showing_association_at_P_5e-04_in_the_NIMH_NIA_ADSP_AD_families.csv", sep =",", header = TRUE)
singlevar_Prokopenko_p0.0005$Prokopenko_Nearest_protein.coding.gene <- as.character(singlevar_Prokopenko_p0.0005$Prokopenko_Nearest_protein.coding.gene)
singlevar_Prokopenko_p0.0005$key <- paste(singlevar_Prokopenko_p0.0005$Prokopenko_Chromosome, singlevar_Prokopenko_p0.0005$Prokopenko_Position, sep = ":")
######################### Comparison ##########################
singlevar_Prokopenko_p0.05 <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/Rare_variants_showing_association_in_Prokopenko_paper_P0.05.csv", sep =",", header = TRUE)
singlevar_Prokopenko_p0.05$Prokopenko_Overlapping_GREAT_associated_genes <- as.character(singlevar_Prokopenko_p0.05$Prokopenko_Overlapping_GREAT_associated_genes)
singlevar_Prokopenko_p0.05$key <- paste(singlevar_Prokopenko_p0.05$Prokopenko_Chromosome, singlevar_Prokopenko_p0.05$Prokopenko_Position, sep = ":")

PROKOPENKO_GENES <- c("FNBP1L", "SEL1L", "LINC00298", "ID2", "C15orf41", "PRKCH", "C2CD3", "KIF2A", "APC", "LHX9", "NALCN", "CTNNA2", "SYTL3", "CLSTN2")

# Prokopenko filter paradigm 1
singlevar_Prokopenko_p0.0005 <- singlevar_Prokopenko_p0.0005[singlevar_Prokopenko_p0.0005$Prokopenko_Nearest_protein.coding.gene %in% PROKOPENKO_GENES,]
singlevar_Prokopenko_p0.0005 <- cbind.data.frame(KEY = singlevar_Prokopenko_p0.0005$key, rsID = singlevar_Prokopenko_p0.0005$Prokopenko_rsID, geneNAME = singlevar_Prokopenko_p0.0005$Prokopenko_Nearest_protein.coding.gene, VEP_GREAT = singlevar_Prokopenko_p0.0005$Prokopenko_VEP_gene_symbol, DiscoveryP = singlevar_Prokopenko_p0.0005$Prokopenko_Discovery_P.value, ReplicationP = singlevar_Prokopenko_p0.0005$Prokopenko_Replication_P.value, Effect_Direction = singlevar_Prokopenko_p0.0005$Prokopenko_Effect.direction, MetaP = singlevar_Prokopenko_p0.0005$Prokopenko_META_P.value, Effect_allele = singlevar_Prokopenko_p0.0005$Prokopenko_Effect_allele, Other_allele = singlevar_Prokopenko_p0.0005$Prokopenko_Other_allele)
singlevar_Prokopenko_p0.0005$Filter <- "paradigm_discoveryP_0.0005"

# Prokopenko filter paradigm 2
singlevar_Prokopenko_p0.05 <- singlevar_Prokopenko_p0.05[singlevar_Prokopenko_p0.05$Prokopenko_Nearest_protein_coding.gene %in% PROKOPENKO_GENES,]
singlevar_Prokopenko_p0.05 <- cbind.data.frame(KEY = singlevar_Prokopenko_p0.05$key, rsID = singlevar_Prokopenko_p0.05$Prokopenko_rsid, geneNAME = singlevar_Prokopenko_p0.05$Prokopenko_Nearest_protein_coding.gene, VEP_GREAT = singlevar_Prokopenko_p0.05$Prokopenko_Overlapping_GREAT_associated_genes, DiscoveryP = singlevar_Prokopenko_p0.05$Prokopenko_DISC_P_value, ReplicationP = singlevar_Prokopenko_p0.05$Prokopenko_REP_P_value, Effect_Direction = singlevar_Prokopenko_p0.05$Prokopenko_Effect_direction, MetaP = singlevar_Prokopenko_p0.05$Prokopenko_META_P_value, Effect_allele = singlevar_Prokopenko_p0.05$Prokopenko_Effect_allele, Other_allele = singlevar_Prokopenko_p0.05$Prokopenko_Other_allele)
singlevar_Prokopenko_p0.05$Filter <- "paradigm_discoveryP_0.05"

PROKOPENKO_SINGLE_VARIANT <- rbind.data.frame(singlevar_Prokopenko_p0.0005, singlevar_Prokopenko_p0.05)
PROKOPENKO_SINGLE_VARIANT <- PROKOPENKO_SINGLE_VARIANT[(PROKOPENKO_SINGLE_VARIANT$DiscoveryP < 0.0005 & PROKOPENKO_SINGLE_VARIANT$ReplicationP < 0.05) | (PROKOPENKO_SINGLE_VARIANT$DiscoveryP < 0.05 & PROKOPENKO_SINGLE_VARIANT$MetaP < 0.0005) ,]
PROKOPENKO_SINGLE_VARIANT$ID <- as.character(PROKOPENKO_SINGLE_VARIANT$KEY)
PROKOPENKO_SINGLE_VARIANT$geneNAME <- as.character(PROKOPENKO_SINGLE_VARIANT$geneNAME)

colnames(PROKOPENKO_SINGLE_VARIANT) [!(grepl("geneNAME|ID", colnames(PROKOPENKO_SINGLE_VARIANT)))] <- paste0(colnames(PROKOPENKO_SINGLE_VARIANT) [!(grepl("geneNAME|ID", colnames(PROKOPENKO_SINGLE_VARIANT)))], "_PROKOPENKO")


## Familial data
SINGLE_VARIANT_ANALYSIS_matched_by_variant_ID_FAMILIAL <- Reduce(function(x,y) merge(x,y,by="ID",all.x= TRUE) ,list(PROKOPENKO_SINGLE_VARIANT,NEUP_FAMILIAL_FBAT_ANNO))
SINGLE_VARIANT_ANALYSIS_matched_by_geneNAME_FAMILIAL <- Reduce(function(x,y) merge(x,y,by="geneNAME",all.x= TRUE) ,list(PROKOPENKO_SINGLE_VARIANT,NEUP_FAMILIAL_FBAT_ANNO))


# write.table(SINGLE_VARIANT_ANALYSIS_matched_by_variant_ID_FAMILIAL, "SINGLE_VARIANT_ANALYSIS_matched_by_variant_ID_FAMILIAL_all.csv", sep =",", col.names = T, quote = F, row.names = FALSE)
# write.table(SINGLE_VARIANT_ANALYSIS_matched_by_geneNAME_FAMILIAL, "SINGLE_VARIANT_ANALYSIS_matched_by_geneNAME_FAMILIAL_all.csv", sep =",", col.names = T, quote = F, row.names = FALSE)

write.table(SINGLE_VARIANT_ANALYSIS_matched_by_variant_ID_FAMILIAL, "SINGLE_VARIANT_ANALYSIS_matched_by_variant_ID_FAMILIAL.csv", sep =",", col.names = T, quote = F, row.names = FALSE)
write.table(SINGLE_VARIANT_ANALYSIS_matched_by_geneNAME_FAMILIAL, "SINGLE_VARIANT_ANALYSIS_matched_by_geneNAME_FAMILIAL.csv", sep =",", col.names = T, quote = F, row.names = FALSE)


# Unrelated data
SINGLE_VARIANT_ANALYSIS_matched_by_variant_ID_UNRELATED <- Reduce(function(x,y) merge(x,y,by="ID",all.x= TRUE) ,list(PROKOPENKO_SINGLE_VARIANT,NEUP_UNRELATED_ANNO))
SINGLE_VARIANT_ANALYSIS_matched_by_geneNAME_UNRELATED <- Reduce(function(x,y) merge(x,y,by="geneNAME",all.x= TRUE) ,list(PROKOPENKO_SINGLE_VARIANT,NEUP_UNRELATED_ANNO))

# write.table(SINGLE_VARIANT_ANALYSIS_matched_by_variant_ID_UNRELATED, "SINGLE_VARIANT_ANALYSIS_matched_by_variant_ID_UNRELATED_all.csv", sep =",", col.names = T, quote = F, row.names = FALSE)
# write.table(SINGLE_VARIANT_ANALYSIS_matched_by_geneNAME_UNRELATED, "SINGLE_VARIANT_ANALYSIS_matched_by_geneNAME_UNRELATED_all.csv", sep =",", col.names = T, quote = F, row.names = FALSE)

write.table(SINGLE_VARIANT_ANALYSIS_matched_by_variant_ID_UNRELATED, "SINGLE_VARIANT_ANALYSIS_matched_by_variant_ID_UNRELATED.csv", sep =",", col.names = T, quote = F, row.names = FALSE)
write.table(SINGLE_VARIANT_ANALYSIS_matched_by_geneNAME_UNRELATED, "SINGLE_VARIANT_ANALYSIS_matched_by_geneNAME_UNRELATED.csv", sep =",", col.names = T, quote = F, row.names = FALSE)

## Meta analysis
SINGLE_VARIANT_ANALYSIS_matched_by_variant_ID_META <- Reduce(function(x,y) merge(x,y,by="ID",all.x= TRUE) ,list(PROKOPENKO_SINGLE_VARIANT,Fixed_Effect_Meta_analysis_result))
SINGLE_VARIANT_ANALYSIS_matched_by_geneNAME_META <- Reduce(function(x,y) merge(x,y,by="geneNAME",all.x= TRUE) ,list(PROKOPENKO_SINGLE_VARIANT,Fixed_Effect_Meta_analysis_result))

# write.table(SINGLE_VARIANT_ANALYSIS_matched_by_variant_ID_META, "SINGLE_VARIANT_ANALYSIS_matched_by_variant_ID_META_all.csv", sep =",", col.names = T, quote = F, row.names = FALSE)
# write.table(SINGLE_VARIANT_ANALYSIS_matched_by_geneNAME_META, "SINGLE_VARIANT_ANALYSIS_matched_by_geneNAME_META_all.csv", sep =",", col.names = T, quote = F, row.names = FALSE)

write.table(SINGLE_VARIANT_ANALYSIS_matched_by_variant_ID_META, "SINGLE_VARIANT_ANALYSIS_matched_by_variant_ID_META.csv", sep =",", col.names = T, quote = F, row.names = FALSE)
write.table(SINGLE_VARIANT_ANALYSIS_matched_by_geneNAME_META, "SINGLE_VARIANT_ANALYSIS_matched_by_geneNAME_META.csv", sep =",", col.names = T, quote = F, row.names = FALSE)



####################################################################################
####################################################################################
####################################################################################

#########################
## gene-based analysis ##
#########################

# Now compare these with Prokopenko's result
PROKOPENKO_SPATIAL <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/Tanzi_spatial_clustering.csv", header = T, sep = ",")
PROKOPENKO_SPATIAL$PROKOPENKO_GENE <- as.character(PROKOPENKO_SPATIAL$PROKOPENKO_GENE)
PROKOPENKO_SPATIAL$geneNAME <- as.character(PROKOPENKO_SPATIAL$PROKOPENKO_GENE)

PROKOPENKO_SPATIAL <- PROKOPENKO_SPATIAL[PROKOPENKO_SPATIAL$geneNAME %in% PROKOPENKO_GENES[-c(1:5)],]

## Familial dataset
## gene-based analysis by MAF
NEUP_FAMILIAL_SKAT_C_MAF <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/02-gene-based/geneset-MAF1PEXAC_SKAT_SKAT_C-GENOME.ASSOC", header = T, sep = "\t", stringsAsFactors = FALSE)
sum(NEUP_FAMILIAL_SKAT_C_MAF$SKAT_C < 0.05, na.rm = T)
NEUP_FAMILIAL_SKAT_C_MAF <- NEUP_FAMILIAL_SKAT_C_MAF[NEUP_FAMILIAL_SKAT_C_MAF$SKAT_C < 0.05,]
colnames(NEUP_FAMILIAL_SKAT_C_MAF) <- paste0("NEUP_SKATC_MAF_FAMILIAL_", colnames(NEUP_FAMILIAL_SKAT_C_MAF))
NEUP_FAMILIAL_SKAT_C_MAF$geneNAME <- NEUP_FAMILIAL_SKAT_C_MAF$NEUP_SKATC_MAF_FAMILIAL_GENE

# Sort by gene name and P value and merge with Prokopenko's gene list
NEUP_FAMILIAL_SKAT_C_MAF <- NEUP_FAMILIAL_SKAT_C_MAF[order(NEUP_FAMILIAL_SKAT_C_MAF$geneNAME, NEUP_FAMILIAL_SKAT_C_MAF$NEUP_SKATC_MAF_FAMILIAL_SKAT_C),]
SPATIAL_CLUSTERING_ANALYSIS_with_MAF <- Reduce(function(x,y) merge(x,y,by="geneNAME",all.x= TRUE) ,list(PROKOPENKO_SPATIAL, NEUP_FAMILIAL_SKAT_C_MAF))


## gene-based analysis by CADD
NEUP_FAMILIAL_SKAT_C_CADD <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/02-gene-based/geneset-CADD20_SKAT_SKAT_C-GENOME.ASSOC", header = T, sep = "\t", stringsAsFactors = FALSE)
sum(NEUP_FAMILIAL_SKAT_C_CADD$SKAT_C < 0.05, na.rm = T)
NEUP_FAMILIAL_SKAT_C_CADD <- NEUP_FAMILIAL_SKAT_C_CADD[NEUP_FAMILIAL_SKAT_C_CADD$SKAT_C < 0.05,]
colnames(NEUP_FAMILIAL_SKAT_C_CADD) <- paste0("NEUP_SKAT_CADD_FAMILIAL_", colnames(NEUP_FAMILIAL_SKAT_C_CADD))
NEUP_FAMILIAL_SKAT_C_CADD$geneNAME  <- NEUP_FAMILIAL_SKAT_C_CADD$NEUP_SKAT_CADD_FAMILIAL_GENE

# Sort by gene name and P value and merge with table above
NEUP_FAMILIAL_SKAT_C_CADD <- NEUP_FAMILIAL_SKAT_C_CADD[order(NEUP_FAMILIAL_SKAT_C_CADD$geneNAME, NEUP_FAMILIAL_SKAT_C_CADD$NEUP_SKAT_CADD_FAMILIAL_SKAT_C),]

SPATIAL_CLUSTERING_ANALYSIS_with_MAF_CADD <- Reduce(function(x,y) merge(x,y,by="geneNAME",all.x= TRUE) ,list(SPATIAL_CLUSTERING_ANALYSIS_with_MAF, NEUP_FAMILIAL_SKAT_C_CADD))





## Unrelated dataset
## gene-based analysis by MAF
NEUP_REPLICATION_SKAT_C_MAF <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/02-gene-based/geneset-MAF1PEXAC_SKAT_SKAT_C-GENOME.ASSOC", header = T, sep = "\t", stringsAsFactors = FALSE)
sum(NEUP_REPLICATION_SKAT_C_MAF$SKAT_C < 0.05, na.rm = T)
NEUP_REPLICATION_SKAT_C_MAF <- NEUP_REPLICATION_SKAT_C_MAF[NEUP_REPLICATION_SKAT_C_MAF$SKAT_C < 0.05,]
colnames(NEUP_REPLICATION_SKAT_C_MAF) <- paste0("NEUP_SKATC_MAF_UNRELATED_", colnames(NEUP_REPLICATION_SKAT_C_MAF))
NEUP_REPLICATION_SKAT_C_MAF$geneNAME <- NEUP_REPLICATION_SKAT_C_MAF$NEUP_SKATC_MAF_UNRELATED_GENE

# Sort by gene name and P value and merge with Prokopenko's gene list
NEUP_REPLICATION_SKAT_C_MAF <- NEUP_REPLICATION_SKAT_C_MAF[order(NEUP_REPLICATION_SKAT_C_MAF$geneNAME, NEUP_REPLICATION_SKAT_C_MAF$NEUP_SKATC_MAF_UNRELATED_SKAT_C),]
SPATIAL_CLUSTERING_ANALYSIS_with_MAF <- Reduce(function(x,y) merge(x,y,by="geneNAME",all.x= TRUE) ,list(SPATIAL_CLUSTERING_ANALYSIS_with_MAF_CADD, NEUP_REPLICATION_SKAT_C_MAF))




## gene-based analysis by CADD
NEUP_REPLICATION_SKAT_C_CADD <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/02-gene-based/geneset-CADD20_SKAT_SKAT_C-GENOME.ASSOC", header = T, sep = "\t", stringsAsFactors = FALSE)
sum(NEUP_REPLICATION_SKAT_C_CADD$SKAT_C < 0.05, na.rm = T)
NEUP_REPLICATION_SKAT_C_CADD <- NEUP_REPLICATION_SKAT_C_CADD[NEUP_REPLICATION_SKAT_C_CADD$SKAT_C < 0.05,]
colnames(NEUP_REPLICATION_SKAT_C_CADD) <- paste0("NEUP_SKAT_CADD_UNRELATED_", colnames(NEUP_REPLICATION_SKAT_C_CADD))
NEUP_REPLICATION_SKAT_C_CADD$geneNAME  <- NEUP_REPLICATION_SKAT_C_CADD$NEUP_SKAT_CADD_UNRELATED_GENE

# Sort by gene name and P value and merge with table above
NEUP_REPLICATION_SKAT_C_CADD <- NEUP_REPLICATION_SKAT_C_CADD[order(NEUP_REPLICATION_SKAT_C_CADD$geneNAME, NEUP_REPLICATION_SKAT_C_CADD$NEUP_SKAT_CADD_UNRELATED_SKAT_C),]

SPATIAL_CLUSTERING_ANALYSIS_with_MAF_CADD <- Reduce(function(x,y) merge(x,y,by="geneNAME",all.x= TRUE) ,list(SPATIAL_CLUSTERING_ANALYSIS_with_MAF, NEUP_REPLICATION_SKAT_C_CADD))

write.table(SPATIAL_CLUSTERING_ANALYSIS_with_MAF_CADD, "GENE_BASED_analysis_MAFF_CADD20_PROKOPENKO_FAMILIAL_Unrelated_data_results_p_0.05.csv", sep =",", col.names = T, quote = F, row.names = FALSE)



################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################
########################
#### Analysis set 2  ###
########################

#################################
#### Single-variant analysis  ###
#################################


# Family (FBAT)
setwd("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files")

NEUP_FAMILIAL_FBAT_ANNO <- read.delim("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/FBAT/FBAT_rare_variant_analysis_results.csv", header = T, sep = "\t", stringsAsFactors = FALSE)
sum(NEUP_FAMILIAL_FBAT_ANNO$P < 0.05, na.rm = T)
colnames(NEUP_FAMILIAL_FBAT_ANNO) <- paste0("NEUP_FAMILIAL_", colnames(NEUP_FAMILIAL_FBAT_ANNO))
# sort
NEUP_FAMILIAL_FBAT_ANNO <- NEUP_FAMILIAL_FBAT_ANNO[order(NEUP_FAMILIAL_FBAT_ANNO$NEUP_FAMILIAL_P),]
NEUP_FAMILIAL_FBAT_ANNO <- NEUP_FAMILIAL_FBAT_ANNO[NEUP_FAMILIAL_FBAT_ANNO$NEUP_FAMILIAL_P < 0.05,]
# write.table(NEUP_FAMILIAL_FBAT_ANNO, "Single_variant_analysis_familial_data_results_p_0.05.csv", sep =",", col.names = T, quote = F, row.names = FALSE)
NEUP_FAMILIAL_FBAT_ANNO$ID <- do.call(paste, c(read.table(text = NEUP_FAMILIAL_FBAT_ANNO$NEUP_FAMILIAL_Marker, sep = ":")[1:2], sep = ":"))
NEUP_FAMILIAL_FBAT_ANNO$geneNAME <- NEUP_FAMILIAL_FBAT_ANNO$NEUP_FAMILIAL_gene


# Replication 
NEUP_UNRELATED_ANNO <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/Firth-Fallback_replication_study_results.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
sum(NEUP_UNRELATED_ANNO$P < 0.05, na.rm = T)
colnames(NEUP_UNRELATED_ANNO) <- paste0("NEUP_UNRELATED_", colnames(NEUP_UNRELATED_ANNO))
NEUP_UNRELATED_ANNO <- NEUP_UNRELATED_ANNO[!is.na(NEUP_UNRELATED_ANNO$NEUP_UNRELATED_P),]
NEUP_UNRELATED_ANNO <- NEUP_UNRELATED_ANNO[order(NEUP_UNRELATED_ANNO$NEUP_UNRELATED_P),]
NEUP_UNRELATED_ANNO <- NEUP_UNRELATED_ANNO[NEUP_UNRELATED_ANNO$NEUP_UNRELATED_P < 0.05,]
# write.table(NEUP_UNRELATED_ANNO, "Single_variant_analysis_unrelated_data_results_p_0.05.csv", sep =",", col.names = T, quote = F, row.names = FALSE)
NEUP_UNRELATED_ANNO$ID <- do.call(paste, c(read.table(text = NEUP_UNRELATED_ANNO$NEUP_UNRELATED_SNP, sep = ":")[1:2], sep = ":"))
NEUP_UNRELATED_ANNO$geneNAME <- NEUP_UNRELATED_ANNO$NEUP_UNRELATED_gene

# Meta-analysis of two datasets
Fixed_Effect_Meta_analysis_result <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/fixed_effect_meta-analysis/Fixed_effect_Meta_analysis_results.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
colnames(Fixed_Effect_Meta_analysis_result) <- paste0("NEUP_META_", colnames(Fixed_Effect_Meta_analysis_result))
Fixed_Effect_Meta_analysis_result <- Fixed_Effect_Meta_analysis_result[order(Fixed_Effect_Meta_analysis_result$NEUP_META_P),]
Fixed_Effect_Meta_analysis_result <- Fixed_Effect_Meta_analysis_result[Fixed_Effect_Meta_analysis_result$NEUP_META_P < 0.05,]
# write.table(Fixed_Effect_Meta_analysis_result, "Single_variant_META_analysis_results_p_0.05.csv", sep =",", col.names = T, quote = F, row.names = FALSE)
Fixed_Effect_Meta_analysis_result$ID <- do.call(paste, c(read.table(text = Fixed_Effect_Meta_analysis_result$NEUP_META_SNP, sep = ":")[1:2], sep = ":"))
Fixed_Effect_Meta_analysis_result$geneNAME <- Fixed_Effect_Meta_analysis_result$NEUP_META_gene


PROKOPENKO_GENES <- c("FNBP1L", "SEL1L", "LINC00298", "ID2", "C15orf41", "PRKCH", "C2CD3", "KIF2A", "APC", "LHX9", "NALCN", "CTNNA2", "SYTL3", "CLSTN2")


NEUP_FAMILIAL_FBAT_ANNO <- NEUP_FAMILIAL_FBAT_ANNO [NEUP_FAMILIAL_FBAT_ANNO$geneNAME %in% PROKOPENKO_GENES,]
NEUP_UNRELATED_ANNO <- NEUP_UNRELATED_ANNO [NEUP_UNRELATED_ANNO$geneNAME %in% PROKOPENKO_GENES,]
Fixed_Effect_Meta_analysis_result <- Fixed_Effect_Meta_analysis_result [Fixed_Effect_Meta_analysis_result$geneNAME %in% PROKOPENKO_GENES,]


SINGLE_VAR_MATCHED_BY_VARIANT <- Reduce(function(x,y) merge(x,y,by="ID",all = TRUE) ,list(NEUP_FAMILIAL_FBAT_ANNO,NEUP_UNRELATED_ANNO, Fixed_Effect_Meta_analysis_result))

write.table(SINGLE_VAR_MATCHED_BY_VARIANT, "Single_variant_analysis_results_for13_genes_p_0.05.csv", sep =",", col.names = T, quote = F, row.names = FALSE)


#########################
## gene-based analysis ##
#########################

# # Now compare these with Prokopenko's result

## Familial dataset
## gene-based analysis by MAF
NEUP_FAMILIAL_SKAT_C_MAF <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/02-gene-based/geneset-MAF1PEXAC_SKAT_SKAT_C-GENOME.ASSOC", header = T, sep = "\t", stringsAsFactors = FALSE)
sum(NEUP_FAMILIAL_SKAT_C_MAF$SKAT_C < 0.05, na.rm = T)
NEUP_FAMILIAL_SKAT_C_MAF <- NEUP_FAMILIAL_SKAT_C_MAF[NEUP_FAMILIAL_SKAT_C_MAF$SKAT_C < 0.05,]
colnames(NEUP_FAMILIAL_SKAT_C_MAF) <- paste0("NEUP_SKATC_MAF_FAMILIAL_", colnames(NEUP_FAMILIAL_SKAT_C_MAF))
NEUP_FAMILIAL_SKAT_C_MAF$geneNAME <- NEUP_FAMILIAL_SKAT_C_MAF$NEUP_SKATC_MAF_FAMILIAL_GENE

# Sort by gene name and P value and merge with Prokopenko's gene list
NEUP_FAMILIAL_SKAT_C_MAF <- NEUP_FAMILIAL_SKAT_C_MAF[order(NEUP_FAMILIAL_SKAT_C_MAF$geneNAME, NEUP_FAMILIAL_SKAT_C_MAF$NEUP_SKATC_MAF_FAMILIAL_SKAT_C),]
# extract Prokopenko's genes
NEUP_FAMILIAL_SKAT_C_MAF <- NEUP_FAMILIAL_SKAT_C_MAF[NEUP_FAMILIAL_SKAT_C_MAF$geneNAME %in% PROKOPENKO_GENES,]


## gene-based analysis by CADD
NEUP_FAMILIAL_SKAT_C_CADD <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/02-gene-based/geneset-CADD20_SKAT_SKAT_C-GENOME.ASSOC", header = T, sep = "\t", stringsAsFactors = FALSE)
sum(NEUP_FAMILIAL_SKAT_C_CADD$SKAT_C < 0.05, na.rm = T)
NEUP_FAMILIAL_SKAT_C_CADD <- NEUP_FAMILIAL_SKAT_C_CADD[NEUP_FAMILIAL_SKAT_C_CADD$SKAT_C < 0.05,]
colnames(NEUP_FAMILIAL_SKAT_C_CADD) <- paste0("NEUP_SKAT_CADD_FAMILIAL_", colnames(NEUP_FAMILIAL_SKAT_C_CADD))
NEUP_FAMILIAL_SKAT_C_CADD$geneNAME  <- NEUP_FAMILIAL_SKAT_C_CADD$NEUP_SKAT_CADD_FAMILIAL_GENE

# Sort by gene name and P value and merge with table above
NEUP_FAMILIAL_SKAT_C_CADD <- NEUP_FAMILIAL_SKAT_C_CADD[order(NEUP_FAMILIAL_SKAT_C_CADD$geneNAME, NEUP_FAMILIAL_SKAT_C_CADD$NEUP_SKAT_CADD_FAMILIAL_SKAT_C),]

# extract Prokopenko's genes
NEUP_FAMILIAL_SKAT_C_CADD <- NEUP_FAMILIAL_SKAT_C_CADD[NEUP_FAMILIAL_SKAT_C_CADD$geneNAME %in% PROKOPENKO_GENES,]





## Unrelated dataset
## gene-based analysis by MAF
NEUP_REPLICATION_SKAT_C_MAF <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/02-gene-based/geneset-MAF1PEXAC_SKAT_SKAT_C-GENOME.ASSOC", header = T, sep = "\t", stringsAsFactors = FALSE)
sum(NEUP_REPLICATION_SKAT_C_MAF$SKAT_C < 0.05, na.rm = T)
NEUP_REPLICATION_SKAT_C_MAF <- NEUP_REPLICATION_SKAT_C_MAF[NEUP_REPLICATION_SKAT_C_MAF$SKAT_C < 0.05,]
colnames(NEUP_REPLICATION_SKAT_C_MAF) <- paste0("NEUP_SKATC_MAF_UNRELATED_", colnames(NEUP_REPLICATION_SKAT_C_MAF))
NEUP_REPLICATION_SKAT_C_MAF$geneNAME <- NEUP_REPLICATION_SKAT_C_MAF$NEUP_SKATC_MAF_UNRELATED_GENE

# Sort by gene name and P value and merge with Prokopenko's gene list
NEUP_REPLICATION_SKAT_C_MAF <- NEUP_REPLICATION_SKAT_C_MAF[order(NEUP_REPLICATION_SKAT_C_MAF$geneNAME, NEUP_REPLICATION_SKAT_C_MAF$NEUP_SKATC_MAF_UNRELATED_SKAT_C),]
# extract Prokopenko's genes
NEUP_REPLICATION_SKAT_C_MAF <- NEUP_REPLICATION_SKAT_C_MAF[NEUP_REPLICATION_SKAT_C_MAF$geneNAME %in% PROKOPENKO_GENES,]





## gene-based analysis by CADD
NEUP_REPLICATION_SKAT_C_CADD <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/02-gene-based/geneset-CADD20_SKAT_SKAT_C-GENOME.ASSOC", header = T, sep = "\t", stringsAsFactors = FALSE)
sum(NEUP_REPLICATION_SKAT_C_CADD$SKAT_C < 0.05, na.rm = T)
NEUP_REPLICATION_SKAT_C_CADD <- NEUP_REPLICATION_SKAT_C_CADD[NEUP_REPLICATION_SKAT_C_CADD$SKAT_C < 0.05,]
colnames(NEUP_REPLICATION_SKAT_C_CADD) <- paste0("NEUP_SKAT_CADD_UNRELATED_", colnames(NEUP_REPLICATION_SKAT_C_CADD))
NEUP_REPLICATION_SKAT_C_CADD$geneNAME  <- NEUP_REPLICATION_SKAT_C_CADD$NEUP_SKAT_CADD_UNRELATED_GENE

# Sort by gene name and P value and merge with table above
NEUP_REPLICATION_SKAT_C_CADD <- NEUP_REPLICATION_SKAT_C_CADD[order(NEUP_REPLICATION_SKAT_C_CADD$geneNAME, NEUP_REPLICATION_SKAT_C_CADD$NEUP_SKAT_CADD_UNRELATED_SKAT_C),]
# extract Prokopenko's genes
NEUP_REPLICATION_SKAT_C_CADD <- NEUP_REPLICATION_SKAT_C_CADD[NEUP_REPLICATION_SKAT_C_CADD$geneNAME %in% PROKOPENKO_GENES,]


## Nothing was found from this analysis ata P 0.05

## For the variants found, I am adding P.change, transcript ID and Effect 
variants_familial <- c("13:101095770:G:A")
variants_unrelated <- c("1:93544273:C:T", "14:81479553:C:A", "11:74034499:G:A", "11:74049368:G:T", "11:74054742:G:A", "11:74054767:T:C", "11:74078100:A:T", "11:74078287:C:T", "11:74109213:C:T", "11:74118187:C:T", "6:158708519:C:T")

# familial
anno_file.familial <- fread("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/02-gene-based/FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1_with_STATUS_nonADSP_post_QC2-snpeff-dbnsfp-FIELDS-simple.txt")
anno_file.familial <- anno_file.familial[anno_file.familial$ID %in% variants_familial,]
# CHROM       POS               ID REF ALT ANN[*].ALLELE  ANN[*].EFFECT ANN[*].IMPACT ANN[*].GENE   ANN[*].GENEID ANN[*].FEATURE ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P dbNSFP_CADD_phred
# 1:    13 101095770 13:101095770:G:A   G   A             A intron_variant      MODIFIER       NALCN ENSG00000102452     transcript  ENST00000251127  c.3163-90C>T             .                 .
# dbNSFP_1000Gp3_AF dbNSFP_ExAC_AF dbNSFP_ExAC_Adj_AF

# unrelated
anno_file.unrelated <- fread("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/02-gene-based/AQUILLA_Brian_2445_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1-WXSm-SCOPEm_with_STATUS-clean3-snpeff-dbnsfp-FIELDS-simple.txt", header = T, sep = "\t")
anno_file.unrelated <- anno_file.unrelated[anno_file.unrelated$ID %in% variants_unrelated,]


## Then look for canonical transcripts: 
## Vicky's comment
# Identify which is the canonical transcript;e.g.https://gnomad.broadinstitute.org/gene/ENSG00000137942?dataset=gnomad_r2_1for FNBP1L the canonical trascsript is "ENST00000271234.7"
# [11:04 AM] Fernandez Hernandez, Victoria
# then, when you report your results, specially the nonsynoimous, refer to the canonical transcript per default or mention if that nonsynoimous effect falls in a different transcript


