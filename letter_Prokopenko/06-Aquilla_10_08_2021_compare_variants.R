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




########################
## Spatial clustering ##
########################





setwd("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/")
# Family (Spatial)
NEUP_DISCOVERY_SPATIAL <- read.delim("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/FBAT/Spatial_clustering_analysis_results.csv", header = T, sep = "\t", stringsAsFactors = FALSE)
sum(NEUP_DISCOVERY_SPATIAL$P < 0.01, na.rm = T) # 529
sum(NEUP_DISCOVERY_SPATIAL$P < 0.05, na.rm = T) #4145

colnames(NEUP_DISCOVERY_SPATIAL) <- paste0("NEUP_SPATIAL_FAMILIAL_", colnames(NEUP_DISCOVERY_SPATIAL))
# sort by P-value
NEUP_DISCOVERY_SPATIAL <- NEUP_DISCOVERY_SPATIAL[order(NEUP_DISCOVERY_SPATIAL$NEUP_SPATIAL_FAMILIAL_P),]
sum(NEUP_DISCOVERY_SPATIAL$NEUP_SPATIAL_FAMILIAL_P < 0.05)
NEUP_DISCOVERY_SPATIAL <- NEUP_DISCOVERY_SPATIAL[NEUP_DISCOVERY_SPATIAL$NEUP_SPATIAL_FAMILIAL_P < 0.05,]
# write.table(NEUP_DISCOVERY_SPATIAL, "Spatial_clustering_analysis_familial_data_results_p_0.05.csv", sep =",", col.names = T, quote = F, row.names = FALSE)

# Now, select regions with P< 0.05
NEUP_FAMILIAL_Cluster.4145.P0.05 <- NEUP_DISCOVERY_SPATIAL[NEUP_DISCOVERY_SPATIAL$NEUP_SPATIAL_FAMILIAL_P < 0.05,]

NEUP_FAMILIAL_Cluster.4145.P0.05$geneNAME <- NEUP_FAMILIAL_Cluster.4145.P0.05$NEUP_SPATIAL_FAMILIAL_gene



# Now compare these with Prokopenko's result
PROKOPENKO_SPATIAL <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/Tanzi_spatial_clustering.csv", header = T, sep = ",")
PROKOPENKO_SPATIAL$PROKOPENKO_GENE <- as.character(PROKOPENKO_SPATIAL$PROKOPENKO_GENE)
PROKOPENKO_SPATIAL$geneNAME <- as.character(PROKOPENKO_SPATIAL$PROKOPENKO_GENE)

# # MERGED_SPATIAL_CLUSTER_ANALYSIS <- Reduce(function(x,y) merge(x,y,by="geneNAME",all= TRUE) ,list(SPATIAL_CLUSTER_RESULT_FAMILIAL_UNRELATED.4145,PROKOPENKO_SPATIAL))
# # PROKOPENKO_SPATIAL_PARADIGM1 <- MERGED_SPATIAL_CLUSTER_ANALYSIS[(MERGED_SPATIAL_CLUSTER_ANALYSIS$PROKOPENKO_DISC_P < 0.0005 & MERGED_SPATIAL_CLUSTER_ANALYSIS$PROKOPENKO_META_P < 0.0005),]
# # PROKOPENKO_SPATIAL_PARADIGM2 <- PROKOPENKO_SPATIAL[PROKOPENKO_SPATIAL$PROKOPENKO_META_P < 0.00005 & PROKOPENKO_SPATIAL$PROKOPENKO_DISC_P < 0.05,]
# 
# PROKOPENKO_SPATIAL_PARADIGM1 <- PROKOPENKO_SPATIAL[(PROKOPENKO_SPATIAL$PROKOPENKO_DISC_P < 0.0005 & PROKOPENKO_SPATIAL$PROKOPENKO_META_P < 0.0005 & PROKOPENKO_SPATIAL$PROKOPENKO_REP_P < 0.05),]
# PROKOPENKO_SPATIAL_PARADIGM1$Filtering <- "Paradigm1"
# # colnames(PROKOPENKO_SPATIAL_PARADIGM1)[colnames(PROKOPENKO_SPATIAL_PARADIGM1) != "geneNAME"] <- paste0(colnames(PROKOPENKO_SPATIAL_PARADIGM1)[colnames(PROKOPENKO_SPATIAL_PARADIGM1) != "geneNAME"], ".PARADIGM1")
# PROKOPENKO_SPATIAL_PARADIGM2 <- PROKOPENKO_SPATIAL[(PROKOPENKO_SPATIAL$PROKOPENKO_META_P < 0.00005 & PROKOPENKO_SPATIAL$PROKOPENKO_DISC_P < 0.05 & PROKOPENKO_SPATIAL$PROKOPENKO_REP_P < 0.05),]
# # colnames(PROKOPENKO_SPATIAL_PARADIGM2)[colnames(PROKOPENKO_SPATIAL_PARADIGM2) != "geneNAME"] <- paste0(colnames(PROKOPENKO_SPATIAL_PARADIGM2)[colnames(PROKOPENKO_SPATIAL_PARADIGM2) != "geneNAME"], ".PARADIGM2")
# PROKOPENKO_SPATIAL_PARADIGM2$Filtering <- "Paradigm2"
# PROKOPENKO_SPATIAL <- rbind.data.frame(PROKOPENKO_SPATIAL_PARADIGM1, PROKOPENKO_SPATIAL_PARADIGM2)

PROKOPENKO_SPATIAL <- PROKOPENKO_SPATIAL[PROKOPENKO_SPATIAL$geneNAME %in% PROKOPENKO_GENES[-c(1:5)],]


## Familial
SPATIAL_CLUSTERING_ANALYSIS <- Reduce(function(x,y) merge(x,y,by="geneNAME",all.x= TRUE) ,list(PROKOPENKO_SPATIAL, NEUP_FAMILIAL_Cluster.4145.P0.05))
# Sort by P value and gene name
SPATIAL_CLUSTERING_ANALYSIS <- SPATIAL_CLUSTERING_ANALYSIS[order(SPATIAL_CLUSTERING_ANALYSIS$NEUP_SPATIAL_FAMILIAL_P, SPATIAL_CLUSTERING_ANALYSIS$geneNAME),]

## Same regions in UNrelated datasets
NEUP_UNRELATED_SPATIAL_CLUSTERING_PROKOPENKO_GENES <-  read.delim("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/02-gene-based/geneset-PROKOPENKO_SKAT_SKAT_C-GENOME.ASSOC", header = T, sep = "\t", stringsAsFactors = FALSE)
NEUP_DISCOVERY_SPATIAL$CLUSTER_ID <- gsub(":", "_", NEUP_DISCOVERY_SPATIAL$NEUP_SPATIAL_FAMILIAL_ClID)
NEUP_UNRELATED_SPATIAL_CLUSTERING_PROKOPENKO_GENES$geneNAME <- NEUP_DISCOVERY_SPATIAL$NEUP_SPATIAL_FAMILIAL_gene [match(NEUP_UNRELATED_SPATIAL_CLUSTERING_PROKOPENKO_GENES$GENE, NEUP_DISCOVERY_SPATIAL$CLUSTER_ID)]

SPATIAL_CLUSTERING_ANALYSIS <- Reduce(function(x,y) merge(x,y,by="geneNAME",all.x= TRUE) ,list(SPATIAL_CLUSTERING_ANALYSIS, NEUP_UNRELATED_SPATIAL_CLUSTERING_PROKOPENKO_GENES))


## gene-based analysis by MAF
NEUP_REPLICATION_SKAT_C_MAF <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/02-gene-based/geneset-MAF1PEXAC_SKAT_SKAT_C-GENOME.ASSOC", header = T, sep = "\t", stringsAsFactors = FALSE)
sum(NEUP_REPLICATION_SKAT_C_MAF$SKAT_C < 0.05, na.rm = T)
NEUP_REPLICATION_SKAT_C_MAF <- NEUP_REPLICATION_SKAT_C_MAF[NEUP_REPLICATION_SKAT_C_MAF$SKAT_C < 0.05,]
colnames(NEUP_REPLICATION_SKAT_C_MAF) <- paste0("NEUP_SKATC_MAF_UNRELATED_", colnames(NEUP_REPLICATION_SKAT_C_MAF))
NEUP_REPLICATION_SKAT_C_MAF$geneNAME <- NEUP_REPLICATION_SKAT_C_MAF$NEUP_SKATC_MAF_UNRELATED_GENE

# Sort by gene name and P value 
NEUP_REPLICATION_SKAT_C_MAF <- NEUP_REPLICATION_SKAT_C_MAF[order(NEUP_REPLICATION_SKAT_C_MAF$geneNAME, NEUP_REPLICATION_SKAT_C_MAF$NEUP_SKATC_MAF_UNRELATED_SKAT_C),]
SPATIAL_CLUSTERING_ANALYSIS_with_MAF <- Reduce(function(x,y) merge(x,y,by="geneNAME",all.x= TRUE) ,list(SPATIAL_CLUSTERING_ANALYSIS, NEUP_REPLICATION_SKAT_C_MAF))




## gene-based analysis by CADD
NEUP_REPLICATION_SKAT_C_CADD <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/02-gene-based/geneset-CADD20_SKAT_SKAT_C-GENOME.ASSOC", header = T, sep = "\t", stringsAsFactors = FALSE)
sum(NEUP_REPLICATION_SKAT_C_CADD$SKAT_C < 0.05, na.rm = T)
NEUP_REPLICATION_SKAT_C_CADD <- NEUP_REPLICATION_SKAT_C_CADD[NEUP_REPLICATION_SKAT_C_CADD$SKAT_C < 0.05,]
colnames(NEUP_REPLICATION_SKAT_C_CADD) <- paste0("NEUP_SKAT_CADD_UNRELATED_", colnames(NEUP_REPLICATION_SKAT_C_CADD))
NEUP_REPLICATION_SKAT_C_CADD$geneNAME  <- NEUP_REPLICATION_SKAT_C_CADD$NEUP_SKAT_CADD_UNRELATED_GENE

# Sort by gene name and P value 
NEUP_REPLICATION_SKAT_C_CADD <- NEUP_REPLICATION_SKAT_C_CADD[order(NEUP_REPLICATION_SKAT_C_CADD$geneNAME, NEUP_REPLICATION_SKAT_C_CADD$NEUP_SKAT_CADD_UNRELATED_SKAT_C),]

SPATIAL_CLUSTERING_ANALYSIS_with_MAF_CADD <- Reduce(function(x,y) merge(x,y,by="geneNAME",all.x= TRUE) ,list(SPATIAL_CLUSTERING_ANALYSIS_with_MAF, NEUP_REPLICATION_SKAT_C_CADD))








########################
## END !!!!!!!!!!!!!!!!!
########################



## Results from gene-based analysis on 4145 regions of unrelated dataset
NEUP_UNRELATED_SKAT_C_4145.CLUSTER <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/02-gene-based/geneset-SPATIALCLUSTER_SKAT_SKAT_C-GENOME.ASSOC", header = T, sep = "\t", stringsAsFactors = FALSE)
colnames(NEUP_UNRELATED_SKAT_C_4145.CLUSTER) <- paste0("NEUP_SPATIAL_UNRELATED_", colnames(NEUP_UNRELATED_SKAT_C_4145.CLUSTER))
NEUP_UNRELATED_SKAT_C_4145.CLUSTER$NEUP_SPATIAL_UNRELATED_CLUSTER_ID <- gsub("_", ":", NEUP_UNRELATED_SKAT_C_4145.CLUSTER$NEUP_SPATIAL_UNRELATED_GENE)

# Merge the results from FAMILIAL and UNRELATED datasets by matching cluster ID
SPATIAL_CLUSTER_RESULT_FAMILIAL_UNRELATED.4145 <- cbind(NEUP_FAMILIAL_Cluster.4145.P0.05, NEUP_UNRELATED_SKAT_C_4145.CLUSTER[match(NEUP_FAMILIAL_Cluster.4145.P0.05$NEUP_SPATIAL_FAMILIAL_ClID, NEUP_UNRELATED_SKAT_C_4145.CLUSTER$NEUP_SPATIAL_UNRELATED_CLUSTER_ID ),])
SPATIAL_CLUSTER_RESULT_FAMILIAL_UNRELATED.4145$geneNAME <- SPATIAL_CLUSTER_RESULT_FAMILIAL_UNRELATED.4145$NEUP_SPATIAL_FAMILIAL_gene

# Perform Fisher's combined probability test
# i=1
SPATIAL_CLUSTER_RESULT_FAMILIAL_UNRELATED.4145$NEUP_SPATIAL_META_P <- ""
for(i in 1:nrow(SPATIAL_CLUSTER_RESULT_FAMILIAL_UNRELATED.4145)){
if (is.na(SPATIAL_CLUSTER_RESULT_FAMILIAL_UNRELATED.4145$NEUP_SPATIAL_UNRELATED_SKAT_C[i])){
    next # skipping NAs
  }
  SPATIAL_CLUSTER_RESULT_FAMILIAL_UNRELATED.4145$NEUP_SPATIAL_META_P[i] <- sumlog(as.vector(c(SPATIAL_CLUSTER_RESULT_FAMILIAL_UNRELATED.4145$NEUP_SPATIAL_FAMILIAL_P[i], SPATIAL_CLUSTER_RESULT_FAMILIAL_UNRELATED.4145$NEUP_SPATIAL_UNRELATED_SKAT_C[i])))$p
}
SPATIAL_CLUSTER_RESULT_FAMILIAL_UNRELATED.4145$NEUP_SPATIAL_META_P <- as.numeric(SPATIAL_CLUSTER_RESULT_FAMILIAL_UNRELATED.4145$NEUP_SPATIAL_META_P)

# sum(SPATIAL_CLUSTER_RESULT_FAMILIAL_UNRELATED.4145$NEUP_SPATIAL_FAMILIAL_P < 0.0005)
# sum(SPATIAL_CLUSTER_RESULT_FAMILIAL_UNRELATED.4145$NEUP_SPATIAL_META_P < 0.0005)

# Now compare these with Prokopenko's result
PROKOPENKO_SPATIAL <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/Tanzi_spatial_clustering.csv", header = T, sep = ",")
PROKOPENKO_SPATIAL$PROKOPENKO_GENE <- as.character(PROKOPENKO_SPATIAL$PROKOPENKO_GENE)
PROKOPENKO_SPATIAL$geneNAME <- as.character(PROKOPENKO_SPATIAL$PROKOPENKO_GENE)

# MERGED_SPATIAL_CLUSTER_ANALYSIS <- Reduce(function(x,y) merge(x,y,by="geneNAME",all= TRUE) ,list(SPATIAL_CLUSTER_RESULT_FAMILIAL_UNRELATED.4145,PROKOPENKO_SPATIAL))
# PROKOPENKO_SPATIAL_PARADIGM1 <- MERGED_SPATIAL_CLUSTER_ANALYSIS[(MERGED_SPATIAL_CLUSTER_ANALYSIS$PROKOPENKO_DISC_P < 0.0005 & MERGED_SPATIAL_CLUSTER_ANALYSIS$PROKOPENKO_META_P < 0.0005),]
# PROKOPENKO_SPATIAL_PARADIGM2 <- PROKOPENKO_SPATIAL[PROKOPENKO_SPATIAL$PROKOPENKO_META_P < 0.00005 & PROKOPENKO_SPATIAL$PROKOPENKO_DISC_P < 0.05,]

PROKOPENKO_SPATIAL_PARADIGM1 <- PROKOPENKO_SPATIAL[(PROKOPENKO_SPATIAL$PROKOPENKO_DISC_P < 0.0005 & PROKOPENKO_SPATIAL$PROKOPENKO_META_P < 0.0005 & PROKOPENKO_SPATIAL$PROKOPENKO_REP_P < 0.05),]
colnames(PROKOPENKO_SPATIAL_PARADIGM1)[colnames(PROKOPENKO_SPATIAL_PARADIGM1) != "geneNAME"] <- paste0(colnames(PROKOPENKO_SPATIAL_PARADIGM1)[colnames(PROKOPENKO_SPATIAL_PARADIGM1) != "geneNAME"], ".PARADIGM1")
PROKOPENKO_SPATIAL_PARADIGM2 <- PROKOPENKO_SPATIAL[(PROKOPENKO_SPATIAL$PROKOPENKO_META_P < 0.00005 & PROKOPENKO_SPATIAL$PROKOPENKO_DISC_P < 0.05 & PROKOPENKO_SPATIAL$PROKOPENKO_REP_P < 0.05),]
colnames(PROKOPENKO_SPATIAL_PARADIGM2)[colnames(PROKOPENKO_SPATIAL_PARADIGM2) != "geneNAME"] <- paste0(colnames(PROKOPENKO_SPATIAL_PARADIGM2)[colnames(PROKOPENKO_SPATIAL_PARADIGM2) != "geneNAME"], ".PARADIGM2")

# MERGED_SPATIAL_CLUSTER_ANALYSIS <- Reduce(function(x,y) merge(x,y,by="geneNAME",all= TRUE) ,list(SPATIAL_CLUSTER_RESULT_FAMILIAL_UNRELATED.4145,PROKOPENKO_SPATIAL_PARADIGM1, PROKOPENKO_SPATIAL_PARADIGM2))
MERGED_SPATIAL_CLUSTER_ANALYSIS_PROKOPENKO <- Reduce(function(x,y) merge(x,y,by="geneNAME",all= TRUE) ,list(PROKOPENKO_SPATIAL_PARADIGM1, PROKOPENKO_SPATIAL_PARADIGM2))
MERGED_SPATIAL_CLUSTER_ANALYSIS_by_geneNAME <- Reduce(function(x,y) merge(x,y,by="geneNAME",all.x = TRUE) ,list(MERGED_SPATIAL_CLUSTER_ANALYSIS_PROKOPENKO, SPATIAL_CLUSTER_RESULT_FAMILIAL_UNRELATED.4145))


###########################################################################
###########################################################################
########################## END !!!!!!!!!!!!!###############################








# Replication (gene-based SKAT RC)
# MAF
NEUP_REPLICATION_SKAT_C_MAF <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/02-gene-based/geneset-MAF1PEXAC_SKAT_SKAT_C-GENOME.ASSOC", header = T, sep = "\t", stringsAsFactors = FALSE)
sum(NEUP_REPLICATION_SKAT_C_MAF$SKAT_C < 0.0005, na.rm = T)
colnames(NEUP_REPLICATION_SKAT_C_MAF) <- paste0("NEUP_SKATC_MAF_UNRELATED_", colnames(NEUP_REPLICATION_SKAT_C_MAF))
# sort by P-value
NEUP_REPLICATION_SKAT_C_MAF <- NEUP_REPLICATION_SKAT_C_MAF[order(NEUP_REPLICATION_SKAT_C_MAF$NEUP_SKATC_MAF_UNRELATED_SKAT_C),]
NEUP_REPLICATION_SKAT_C_MAF <- NEUP_REPLICATION_SKAT_C_MAF[NEUP_REPLICATION_SKAT_C_MAF$NEUP_SKATC_MAF_UNRELATED_SKAT_C < 0.01,]
write.table(NEUP_REPLICATION_SKAT_C_MAF, "GENE_BASED_analysis_MAF0.01_Unrelated_data_results_p_0.01.csv", sep =",", col.names = T, quote = F, row.names = FALSE)


# CADD
NEUP_REPLICATION_SKAT_C_CADD <- read.delim("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/02-gene-based/geneset-CADD20_SKAT_SKAT_C-GENOME.ASSOC", header = T, sep = "\t", stringsAsFactors = FALSE)
sum(NEUP_REPLICATION_SKAT_C_CADD$SKAT_C < 0.01, na.rm = T)
colnames(NEUP_REPLICATION_SKAT_C_CADD) <- paste0("NEUP_SKAT_CADD_UNRELATED_", colnames(NEUP_REPLICATION_SKAT_C_CADD))
NEUP_REPLICATION_SKAT_C_CADD <- NEUP_REPLICATION_SKAT_C_CADD[order(NEUP_REPLICATION_SKAT_C_CADD$NEUP_SKAT_CADD_UNRELATED_GENE, NEUP_REPLICATION_SKAT_C_CADD$NEUP_SKAT_CADD_UNRELATED_SKAT_C),]
NEUP_REPLICATION_SKAT_C_CADD$NEUP_SKATC_CADD_UNRELATED_GENE_NAME <-  NEUP_REPLICATION_SKAT_C_CADD$NEUP_SKAT_CADD_UNRELATED_GENE
# sort by P-value
NEUP_REPLICATION_SKAT_C_CADD <- NEUP_REPLICATION_SKAT_C_CADD[order(NEUP_REPLICATION_SKAT_C_CADD$NEUP_SKAT_CADD_UNRELATED_SKAT_C),]
NEUP_REPLICATION_SKAT_C_CADD <- NEUP_REPLICATION_SKAT_C_CADD[NEUP_REPLICATION_SKAT_C_CADD$NEUP_SKAT_CADD_UNRELATED_SKAT_C < 0.01,]
write.table(NEUP_REPLICATION_SKAT_C_CADD, "GENE_BASED_analysis_CADD20_Unrelated_data_results_p_0.01.csv", sep =",", col.names = T, quote = F, row.names = FALSE)


######################### Comparison ##########################
sum(NEUP_DISCOVERY_SPATIAL$NEUP_SPATIAL_DISCgene %in% NEUP_REPLICATION_SKAT_C_MAF$NEUP_SKATC_MAF_REPGENE)
sum(NEUP_DISCOVERY_SPATIAL$NEUP_SPATIAL_DISCgene %in% NEUP_REPLICATION_SKAT_C_CADD$NEUP_SKAT_CADD_REPGENE)


(m1 <- merge(NEUP_DISCOVERY_SPATIAL, NEUP_REPLICATION_SKAT_C_MAF, by.x = "NEUP_SPATIAL_DISCgene", by.y = "NEUP_SKATC_MAF_REPGENE", all.x = TRUE))

# m2 is the dataframe that has spatial clustering analysis results of familial data
# merged with the (top p value) corresponding genes in CADD and MAF
(m2 <- merge(m1, NEUP_REPLICATION_SKAT_C_CADD, by.x = "NEUP_SPATIAL_DISCgene", by.y = "NEUP_SKAT_CADD_REPGENE", all.x = TRUE))


# Now compare these with Prokopenko's result
PROKOPENKO_SPATIAL <- read.table("https://raw.githubusercontent.com/achalneupane/data/master/Tanzi_spatial_clustering.csv", header = T, sep = ",")
PROKOPENKO_SPATIAL$PROKOPENKO_GENE <- as.character(PROKOPENKO_SPATIAL$PROKOPENKO_GENE)

sum(m2$NEUP_SPATIAL_DISCgene %in% PROKOPENKO_SPATIAL$PROKOPENKO_GENE)

unique(PROKOPENKO_SPATIAL$PROKOPENKO_GENE) 

sum(unique(PROKOPENKO_SPATIAL$PROKOPENKO_GENE) %in% m2$NEUP_SPATIAL_DISCgene)






###########################################################################################
########################## MATCH VARIANTS FROM SPATIAL clustering #########################
###########################################################################################




Spatial_replicated_Tanzi <- read.table("Spatial_replicated_Tanzi.csv", sep =",", header = TRUE)

singlevar_Tanzi$Nearest.protein.coding.gene <- as.character(singlevar_Tanzi$Prokopenko_Nearest_protein.coding.gene)

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
  gr=GRanges(seqnames=singlevar_Tanzi$Prokopenko_Chromosome[i],
             ranges=IRanges(start = singlevar_Tanzi$Prokopenko_Position[i]-100000, end = singlevar_Tanzi$Prokopenko_Position[i]+100000))
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