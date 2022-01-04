## Perform gene-based analysis by spatial clusters
# Family (Spatial)
#cd /100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/02-gene-based/chr1
# mkdir geneset-SPATIALCLUSTER
NEUP_DISCOVERY_SPATIAL <- read.delim("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/FBAT/Spatial_clustering_analysis_results.csv", header = T, sep = "\t", stringsAsFactors = FALSE)
colnames(NEUP_DISCOVERY_SPATIAL) <- paste0("NEUP_SPATIAL_FAMILIAL_", colnames(NEUP_DISCOVERY_SPATIAL))

PROKOPENKO_GENES <- c("FNBP1L", "SEL1L", "LINC00298", "ID2", "PRKCH" , "C15orf41" , "C2CD3", "KIF2A", "APC", "LHX9", "NALCN", "CTNNA2", "SYTL3", "CLSTN2")


NEUP_DISCOVERY_SPATIAL_geneset <- NEUP_DISCOVERY_SPATIAL[NEUP_DISCOVERY_SPATIAL$NEUP_SPATIAL_FAMILIAL_gene %in% PROKOPENKO_GENES,]


NEUP_DISCOVERY_SPATIAL_geneset$genesetID <- gsub(":","_", NEUP_DISCOVERY_SPATIAL_geneset$NEUP_SPATIAL_FAMILIAL_ClID)


NEUP_DISCOVERY_SPATIAL_geneset$NEUP_SPATIAL_FAMILIAL_CLUSTER

## First, we expand rows for each cluster to create genesets with variants included by that cluster
require(data.table)
TT <- as.data.table(NEUP_DISCOVERY_SPATIAL_geneset)
TT <- TT[ , list( NEUP_SPATIAL_FAMILIAL_CLUSTER = unlist( strsplit( NEUP_SPATIAL_FAMILIAL_CLUSTER , "-" ) ) ) , by = list(genesetID, NEUP_SPATIAL_FAMILIAL_ClID, NEUP_SPATIAL_FAMILIAL_consequence, NEUP_SPATIAL_FAMILIAL_gene, NEUP_SPATIAL_FAMILIAL_CHR, NEUP_SPATIAL_FAMILIAL_type)]
GENESET <- unique(TT$genesetID)

# Now, create geneset files
PATH <- "/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/02-gene-based/"
for (i in 1:length(GENESET)){
  print(paste0("Doing i: ", i))
  df <- as.data.frame(TT[grepl(GENESET[i], TT$genesetID), c("genesetID", "NEUP_SPATIAL_FAMILIAL_CHR", "NEUP_SPATIAL_FAMILIAL_CLUSTER")])
  CHR <- paste0("chr", df$NEUP_SPATIAL_FAMILIAL_CHR[1])
  CHR <- gsub("chr23", "chrX", CHR)
  CHR <- gsub("chr24", "chrY", CHR)
  GENESET_PATH <- paste0(PATH, CHR, "/geneset-PROKOPENKO")
  if (!dir.exists(GENESET_PATH)) {dir.create(GENESET_PATH)}
  GENESET_FILE <- paste(df$genesetID[1], "gene", sep =".")
  write.table(df$NEUP_SPATIAL_FAMILIAL_CLUSTER, paste0(GENESET_PATH,"/",GENESET_FILE), col.names = F, row.names = F, quote = F)
}


################################################# Only on those with P < 0.05

## Perform gene-based analysis by spatial clusters
# Family (Spatial)
#cd /100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/02-gene-based/chr1
# mkdir geneset-SPATIALCLUSTER
NEUP_DISCOVERY_SPATIAL <- read.delim("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/FBAT/Spatial_clustering_analysis_results.csv", header = T, sep = "\t", stringsAsFactors = FALSE)
sum(NEUP_DISCOVERY_SPATIAL$P < 0.05, na.rm = T)
colnames(NEUP_DISCOVERY_SPATIAL) <- paste0("NEUP_SPATIAL_FAMILIAL_", colnames(NEUP_DISCOVERY_SPATIAL))

NEUP_DISCOVERY_SPATIAL_geneset <- NEUP_DISCOVERY_SPATIAL[NEUP_DISCOVERY_SPATIAL$NEUP_SPATIAL_FAMILIAL_P < 0.05,]


NEUP_DISCOVERY_SPATIAL_geneset$genesetID <- gsub(":","_", NEUP_DISCOVERY_SPATIAL_geneset$NEUP_SPATIAL_FAMILIAL_ClID)


NEUP_DISCOVERY_SPATIAL_geneset$NEUP_SPATIAL_FAMILIAL_CLUSTER

## First, we expand rows for each cluster to create genesets with variants included by that cluster
require(data.table)
TT <- as.data.table(NEUP_DISCOVERY_SPATIAL_geneset)
TT <- TT[ , list( NEUP_SPATIAL_FAMILIAL_CLUSTER = unlist( strsplit( NEUP_SPATIAL_FAMILIAL_CLUSTER , "-" ) ) ) , by = list(genesetID, NEUP_SPATIAL_FAMILIAL_ClID, NEUP_SPATIAL_FAMILIAL_consequence, NEUP_SPATIAL_FAMILIAL_gene, NEUP_SPATIAL_FAMILIAL_CHR, NEUP_SPATIAL_FAMILIAL_type)]
GENESET <- unique(TT$genesetID)

# Now, create geneset files
PATH <- "/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/06-Aquilla_202101-b/03-plink-QC-files/02-gene-based/"
for (i in 1:length(GENESET)){
print(paste0("Doing i: ", i))
df <- as.data.frame(TT[grepl(GENESET[i], TT$genesetID), c("genesetID", "NEUP_SPATIAL_FAMILIAL_CHR", "NEUP_SPATIAL_FAMILIAL_CLUSTER")])
CHR <- paste0("chr", df$NEUP_SPATIAL_FAMILIAL_CHR[1])
CHR <- gsub("chr23", "chrX", CHR)
CHR <- gsub("chr24", "chrY", CHR)
GENESET_PATH <- paste0(PATH, CHR, "/geneset-SPATIALCLUSTER")
  if (!dir.exists(GENESET_PATH)) {dir.create(GENESET_PATH)}
GENESET_FILE <- paste(df$genesetID[1], "gene", sep =".")
write.table(df$NEUP_SPATIAL_FAMILIAL_CLUSTER, paste0(GENESET_PATH,"/",GENESET_FILE), col.names = F, row.names = F, quote = F)
}
