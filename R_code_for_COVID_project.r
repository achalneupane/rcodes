MAP_ID <- read.table("https://raw.githubusercontent.com/achalneupane/data/master/MasterID_NACC_MAP_20210618.csv", header = T, sep = ",", stringsAsFactors = FALSE)
head(MAP_ID)
# View(MAP_ID)
MAP_ID$NACC_ID <- as.character(MAP_ID$NACC_ID)
MAP_ID$NACC_ID <- gsub("\\s", "", MAP_ID$NACC_ID)
colnames(MAP_ID)[4] <- "NACCID"

NACC <- read.delim("https://raw.githubusercontent.com/achalneupane/data/master/NACCphenotypes_COVID_CC.csv", header = T, sep = ",", stringsAsFactors = FALSE)
# NACC <- read.table("/home/achal/Aquilla_7983_merged_vcfs/NACCphenotypes_COVID_CC_WASHU.csv", header = T, sep = ",", stringsAsFactors = FALSE)
head(NACC)

# View(NACC)

NACC$NACCID <- as.character(NACC$NACCID)
NACC$NACCID <- gsub("\\s", "", NACC$NACCID)

# MAP_NACC <- merge(x=NACC, y=MAP_ID, by = "NACCID")
NACC$ADCLocation <- as.character(NACC$ADCLocation)

MAP_NACC <- cbind(NACC,MAP_ID[match(NACC$NACCID, MAP_ID$NACCID),])

dim(MAP_NACC)
# View(MAP_NACC)



Aquilla<- read.table("/home/achal/Aquilla_7983_merged_vcfs/ID_LIST_8835_edited.txt", header = FALSE, sep = ",")
dim(Aquilla)
# View(Aquilla)
Aquilla <- cbind.data.frame(SM=Aquilla$V1, PR=Aquilla$V3)

EOAD <- read.table("/home/achal/Aquilla_7983_merged_vcfs/202004_USUHS_EOAD-WGS_gDNA_EOLUS_SM_DNA_PR_RGBASE_R1_R2_FULLSM_metafile.csv", header = FALSE, sep = ",")
dim(EOAD)
# View(EOAD)
EOAD <- cbind.data.frame(SM=EOAD$V1, PR=EOAD$V3)


FUS <- read.table("/home/achal/Aquilla_7983_merged_vcfs/202103_ADSP_FUS-familial_WGS_UNNAMED-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc3.csv", header = FALSE, sep = ",")
dim(FUS)
# View(FUS)
FUS <- cbind.data.frame(SM=gsub("-BL-NCR","", FUS$V2), PR= FUS$V4)

SITE27 <- read.table("/home/achal/Aquilla_7983_merged_vcfs/202104_ADSP_site27-sync-n303_WGS_UNNAMED-MGIID-sm-dna-pr-rgbase-cram-smdir-cramloc.csv", header = FALSE, sep = ",")
dim(SITE27)
# View(SITE27)

SITE27 <- cbind.data.frame(SM=gsub("-BL-NCR","", SITE27$V2), PR=SITE27$V4)

ALL_SM <- rbind.data.frame(Aquilla, EOAD, FUS, SITE27)

# View(ALL_SM)

ALL_SM$SM <- as.character(ALL_SM$SM)
ALL_SM$SM <- gsub("\\s", "", ALL_SM$SM)
sum(MAP_NACC$NACCID %in% ALL_SM$SM)

ALL_SM$MAP_NUM <- gsub("MAP_", "", ALL_SM$SM)

MAP_NACC$WXS_available <- ifelse(MAP_NACC$MAP_ID %in%  ALL_SM$MAP_NUM, "yes", "no")
ALL_SM$PR <- as.character(ALL_SM$PR)

MAP_NACC$PR <-""

MAP_NACC$PR <- ALL_SM$PR[match(MAP_NACC$MAP_ID,  ALL_SM$MAP_NUM)]

# MAP_NACC$ADCLocation <- as.character(MAP_NACC$ADCLocation)
MAP_NACC$WXS <- MAP_NACC$PR
table(MAP_NACC$PR)

table(MAP_NACC$WXS_available)


MAP_NACC$WXS [MAP_NACC$PR=="Genentech_WES"] <- "WES"
MAP_NACC$WXS [MAP_NACC$PR=="201909_MGI_gDNA_LINDSEY"] <- "WGS"
MAP_NACC$WXS [MAP_NACC$PR=="202004_USUHS_EOAD-WGS_gDNA_EOLUS"] <- "WGS"
MAP_NACC$WXS [MAP_NACC$PR=="Broad_WGS"] <- "WGS"
MAP_NACC$WXS [MAP_NACC$PR=="Genentech_WGS"] <- "WGS"
MAP_NACC$WXS [MAP_NACC$PR=="MGI_FASeEOAD_201605"] <- "WGS"
MAP_NACC$WXS [MAP_NACC$PR=="MGI_Gregg_201704"] <- "WES"
MAP_NACC$WXS [MAP_NACC$PR=="MGI_Imaging_201612"] <- "WES"
MAP_NACC$WXS [MAP_NACC$PR=="phs000572_201508"] <- "WES"


write.table(MAP_NACC, file ="/home/achal/Aquilla_7983_merged_vcfs/NACC_Phenotypes_COVID_AN.csv", sep = ",", col.names = TRUE, row.names = FALSE)
