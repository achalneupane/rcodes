# # Veryify variants data
ADFTDPDmutDB_patho <- read.csv("https://raw.githubusercontent.com/achalneupane/data/master/ADFTDPD_patho.csv", header = TRUE,  stringsAsFactors=FALSE, fileEncoding="latin1", check.names = FALSE)
ADFTDPDmutDB_possible <- read.csv("https://raw.githubusercontent.com/achalneupane/data/master/ADFTDPD-mutDB_Possible.csv", header = TRUE,  stringsAsFactors=FALSE, fileEncoding="latin1", check.names = FALSE)
pathoYes <- read.csv("https://raw.githubusercontent.com/achalneupane/data/master/pathoYes.csv", header = TRUE, check.names = FALSE)
pathoYes_Sum <- read.csv("https://raw.githubusercontent.com/achalneupane/data/master/pathoYes_sum.csv", header = TRUE, check.names = FALSE)


## exome

getwd()
setwd("/home/achal/carlos_task/AGRF/exome/")
# read two times the vcf file, first for the columns names, second for the data
tmp_vcf<-readLines("21310.VEP.tsv")
tmp_vcf_data<-read.table("21310.VEP.tsv", stringsAsFactors = FALSE)

# filter for the columns names 
# tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))] 
tmp_vcf<-tmp_vcf[-(grep("#Uploaded_variation",tmp_vcf)+1):-(length(tmp_vcf))]

vcf_names<-unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
names(tmp_vcf_data)<-vcf_names


colnames(tmp_vcf_data)

dim(tmp_vcf_data)

# tmp.vep <- head(tmp_vcf_data, 100)



sum(grepl ("non_coding", tmp_vcf_data$Consequence))

# final.vep <- tmp_vcf_data[!grepl("non_synonymous", tmp_vcf_data$Consequence),]
final.vep <- tmp_vcf_data
dim(final.vep)
sum(grepl("coding", final.vep$Consequence))
# coding non-coding both
# final.coding <- final.vep[grepl("coding", final.vep$Consequence),]


# mendelian <- c("APP", "PSEN1", "PSEN2", "GRN", "MAPT", "TARDBP", "VCP",  "C9ORF72", "LRRK2", "PARK2", "PINK1")

mendelian <- c("APP", "PSEN1", "PSEN2", "GRN", "MAPT", "TARDBP", "VCP", "LRRK2", "PARK2", "PINK1", 
               "PRNP", "CHMP2B", "FUS", "TBK1", "PARK7", "SNCA", "UCHL1", "ATP13A2", "GIGYF2", "HTRA2", 
               "PLA2G6", "FBXO7", "VPS35", "EIF4G1", "DNAJC16")
# APOE

sum(final.vep$SYMBOL %in% mendelian)

mendelian.vars <- final.vep[final.vep$SYMBOL %in% mendelian, ]





## Types
no.downstream <- mendelian.vars[!grepl("downstream_", mendelian.vars$Consequence),]
dim(no.downstream)
no.downstream_upstream <- no.downstream[!grepl("upstream_", no.downstream$Consequence), ]
dim(no.downstream_upstream)
no.downstream_upstream_intron <- no.downstream_upstream[!grepl("intron_", no.downstream_upstream$Consequence), ]
dim(no.downstream_upstream_intron)
no.downstream_upstream_intron_UTR <- no.downstream_upstream_intron[!grepl("prime_UTR", no.downstream_upstream_intron$Consequence), ]
dim(no.downstream_upstream_intron_UTR)
no.downstream_upstream_intron_UTR_coding_only <- no.downstream_upstream_intron_UTR[!grepl("non_coding", no.downstream_upstream_intron_UTR$Consequence), ]
dim(no.downstream_upstream_intron_UTR_coding_only) # 181 variants


# # mendelian.vars.non_syn <- mendelian.vars[mendelian.vars$Consequence != "synonymous_variant", ]
# # mendelian.vars.non_syn_coding <- mendelian.vars.non_syn[mendelian.vars.non_syn$Consequence != "intron_variant", ]
# 
# exclude <-  c("downstream_gene_variant", "3_prime_UTR_variant", "5_prime_UTR_variant",
#               "upstream_gene_variant", "splice_region_variant,intron_variant",
#               "splice_region_variant,synonymous_variant", "splice_region_variant,5_prime_UTR_variant", 
#               "intron_variant,non_coding_transcript_variant", "non_coding_transcript_exon_variant")
# 
# mendelian_codin_non_synonymous <- mendelian.vars.non_syn_coding[(!mendelian.vars.non_syn_coding$Consequence %in% exclude),]


# varID <- gsub("_|\\/", ":", mendelian.vars$`#Uploaded_variation`)
varID <- paste(gsub("chr", "", no.downstream_upstream_intron_UTR_coding_only$Location), no.downstream_upstream_intron_UTR_coding_only$USED_REF, no.downstream_upstream_intron_UTR_coding_only$Allele, sep = ":")
no.downstream_upstream_intron_UTR_coding_only <- cbind(varID = varID, no.downstream_upstream_intron_UTR_coding_only)

# # PathoYes_SUM
# pathoYes_Sum$GRCh37 <- as.character(pathoYes_Sum$GRCh37)
# n <- 2
# pat <- paste0('^([^:]+(?::[^:]+){',n-1,'}).*')
# pathoSumLocation <- sub(pat, '\\1', pathoYes_Sum$GRCh37)
# 
# 
# pathoYes_Sum$Location <- pathoSumLocation
# # mapped locations
# match(as.character(pathoYes_Sum$Location), gsub("chr", "", no.downstream_upstream_intron_UTR_coding_only$Location))
# sum(no.downstream_upstream_intron_UTR_coding_only$varID %in% as.character(pathoYes_Sum$GRCh37))

## patho_Yes



pathoYes$GRCh37 <- as.character(pathoYes$GRCh37)
n <- 2
pat <- paste0('^([^:]+(?::[^:]+){',n-1,'}).*')
pathoYesLocation <- sub(pat, '\\1', pathoYes$GRCh37)


pathoYes$Location <- pathoYesLocation

# match by locations
match(pathoYes$Location, gsub("chr", "", no.downstream_upstream_intron_UTR_coding_only$Location)) 
sum(!is.na(match(pathoYes$Location, gsub("chr", "", no.downstream_upstream_intron_UTR_coding_only$Location))))

#match by RsID
match(as.character(pathoYes$`rs#`), gsub("rs", "", no.downstream_upstream_intron_UTR_coding_only$`#Uploaded_variation`))
sum(!is.na(match(as.character(pathoYes$`rs#`), gsub("rs", "", no.downstream_upstream_intron_UTR_coding_only$`#Uploaded_variation`))))

write.table(mendelian_codin_non_synonymous, "Exome_AD_FTD_PD_non_synonymous_coding_variants.tsv", col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mendelian.vars, "Exome_AD_FTD_PD_all_variants.tsv", col.names = TRUE, quote = FALSE, sep = "\t")

# # all *synonymous* variants
# sum(grepl("synonymous", tmp_vcf_data$Consequence))
# all_syno <- tmp_vcf_data[grepl("synonymous", tmp_vcf_data$Consequence), ]
# # all.syno.tmp <- head(all_syno, 100)
# 
# sum(grepl("nonsynonymous", all_syno$Consequence))
# 
# 
# # protein coding variants
# protein_coding <- tmp_vcf_data[grepl("protein_coding", tmp_vcf_data$BIOTYPE), ]
# 
# protein_coding[protein_coding$SYMBOL %in% mendelian, ]
# 
# # length(unique(protein_coding$SYMBOL))
# 
# tmp_vcf_data


































#####genome


setwd("/home/achal/carlos_task/AGRF/genome/")


## Two large to read
############
# # read two times the vcf file, first for the columns names, second for the data
# tmp_vcfG <- readLines("21310.VEP.tsv")
# tmp_vcf_dataG <- read.table("21310.VEP.tsv", stringsAsFactors = FALSE)
# 
# # filter for the columns names 
# # tmp_vcf <- tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))] 
# tmp_vcfG <- tmp_vcf[-(grep("#Uploaded_variation",tmp_vcf)+1):-(length(tmp_vcf))]
############

# reading by lines and patterns instead:
library(data.table)

# paste(mendelian,collapse='\\|')
# tmp_vcf_dataG <- fread(cmd=paste("grep", mendelian, "21310.VEP.tsv"))
tmp_vcf_dataG <- fread(cmd=paste("grep", "coding", "21310.VEP.tsv"))
# tmp_vcf_dataG <- fread(cmd=paste("grep", "coding_sequence_variant", "21310.VEP.tsv"))
mycol <- fread(cmd=paste("grep", "\\#Uploaded_variation", "21310.VEP.tsv"))

colnames(tmp_vcf_dataG) <- colnames(mycol)

dim(tmp_vcf_dataG)

# final.vepG <- tmp_vcf_dataG
# mendelian.varsG <- final.vepG[final.vepG$SYMBOL %in% mendelian, ]
mendelian.varsG <- tmp_vcf_dataG[tmp_vcf_dataG$SYMBOL %in% mendelian, ]

dim(mendelian.varsG)



## Types
no.downstreamG <- mendelian.varsG[!grepl("downstream_", mendelian.varsG$Consequence),]
dim(no.downstreamG)
no.downstream_upstreamG <- no.downstreamG[!grepl("upstream_", no.downstreamG$Consequence), ]
dim(no.downstream_upstreamG)
no.downstream_upstream_intronG <- no.downstream_upstreamG[!grepl("intron_", no.downstream_upstreamG$Consequence), ]
dim(no.downstream_upstream_intronG)
no.downstream_upstream_intron_UTR.G <- no.downstream_upstream_intronG[!grepl("prime_UTR", no.downstream_upstream_intronG$Consequence), ]
dim(no.downstream_upstream_intron_UTR.G)
no.downstream_upstream_intron_UTR_coding_onlyG <- no.downstream_upstream_intron_UTR.G[!grepl("non_coding", no.downstream_upstream_intron_UTR.G$Consequence), ]
dim(no.downstream_upstream_intron_UTR_coding_onlyG) # 181 variants




varID.G <- paste(gsub("chr", "", no.downstream_upstream_intron_UTR_coding_onlyG$Location), no.downstream_upstream_intron_UTR_coding_onlyG$USED_REF, no.downstream_upstream_intron_UTR_coding_onlyG$Allele, sep = ":")
no.downstream_upstream_intron_UTR_coding_onlyG <- cbind(varID.G = varID.G, no.downstream_upstream_intron_UTR_coding_onlyG)



## match vars from WES and WGS
match(no.downstream_upstream_intron_UTR_coding_only$varID, no.downstream_upstream_intron_UTR_coding_onlyG$varID.G)

# match by locations
match(pathoYes$Location, gsub("chr", "", no.downstream_upstream_intron_UTR_coding_onlyG$Location)) 
sum(!is.na(match(pathoYes$Location, gsub("chr", "", no.downstream_upstream_intron_UTR_coding_only$Location))))

# match by rsID
match(as.character(pathoYes$`rs#`), gsub("rs", "", no.downstream_upstream_intron_UTR_coding_onlyG$`#Uploaded_variation`))
sum(!is.na(match(as.character(pathoYes$`rs#`), gsub("rs", "", no.downstream_upstream_intron_UTR_coding_only$`#Uploaded_variation`))))


#### combine WES WGS
no.downstream_upstream_intron_UTR_coding_only$rsID <- no.downstream_upstream_intron_UTR_coding_only$`#Uploaded_variation`
WES <- no.downstream_upstream_intron_UTR_coding_only[, c("varID", "Location", "SYMBOL", "rsID", "HGVSp", "HGVSc", "Consequence", "CDS_position", "Protein_position", "Amino_acids", "IMPACT")]
WES$data <- "WES"
WES$varID <- as.character(paste0("chr", WES$varID))

no.downstream_upstream_intron_UTR_coding_onlyG$rsID <- no.downstream_upstream_intron_UTR_coding_onlyG$`#Uploaded_variation`
WGS <- no.downstream_upstream_intron_UTR_coding_onlyG[, c("varID.G", "Location", "SYMBOL", "rsID", "HGVSp", "HGVSc", "Consequence", "CDS_position", "Protein_position", "Amino_acids", "IMPACT")]
WGS$data <- "WGS"



chr <- sapply(strsplit(WGS$Location,":"), `[`, 1)
POS <- sapply(strsplit(WGS$Location,":"), `[`, 2)
positions <- !grepl("-", POS)

POS[positions] <- paste(POS[positions], POS[positions], sep="-")
build38 <- paste(chr, POS, sep = ":")

build37converted <- read.table("https://raw.githubusercontent.com/achalneupane/data/master/build37", header = FALSE, stringsAsFactors = FALSE)
build37converted <- as.matrix(build37converted)


WGS$Location[positions] <- substr(build37converted[positions], 1, regexpr("-", build37converted[positions])-1)

WGS$varID.G <- paste(WGS$Location, sub("^([^:]+:){2}", "", WGS$varID.G), sep = ":")

WESyes <- WES$varID %in% WGS$varID.G
WGSYes <- WGS$varID.G %in% WES$varID 


ALLvars <- rbind(WES,WGS, use.names = FALSE)
# ALLvars <- ALLvars[!grepl("APOE", ALLvars$SYMBOL), ]
ALLvars$Both_WGS_and_WES <- ifelse(c(WESyes, WGSYes), "Yes", "No") 

ALLvars$rsID[grepl("chr", ALLvars$rsID)] <- NA
View(ALLvars)



#######################
# # get cordinates from build 38 for all
ALLvars$Location
chr <- sapply(strsplit(ALLvars$Location,":"), `[`, 1)
POS <- sapply(strsplit(ALLvars$Location,":"), `[`, 2)
positions <- !grepl("-", POS)

POS[positions] <- paste(POS[positions], POS[positions], sep="-")
build37 <- paste(chr, POS, sep = ":")

build38 <- read.table("https://raw.githubusercontent.com/achalneupane/data/master/build_38", header = FALSE)
build38 <- as.matrix(build38)
ALLvars$Location_build38 <- ALLvars$Location
ALLvars$Location_build38[positions] <- substr(build38[positions], 1, regexpr("-", build38[positions])-1)
# WGS$varID.G <- paste(WGS$Location, sub("^([^:]+:){2}", "", WGS$varID.G), sep = ":")

ALLvars$p.change <- ifelse(grepl("%3D", ALLvars$HGVSp), NA, ALLvars$HGVSp) 
ALLvars$p.change

ALLvars$p.change <- sub("^([^:]+:){1}", "", ALLvars$p.change)

ALLvars$c.change <- sub("^([^:]+:){1}", "", ALLvars$HGVSc)

ALLvars <- ALLvars[, c("varID", "Location", "SYMBOL", "c.change", "p.change", "rsID", "Consequence", "data", "Both_WGS_and_WES", "CDS_position", "Protein_position", "Amino_acids", "Location_build38")]
dim(ALLvars)
ALLvars <- unique(ALLvars)
# ALLvars$Both_WGS_and_WES <- ifelse(grepl("Present", ALLvars$Both_WGS_and_WES), "Yes", "No")



##### match by locations 
# pathoYes
match(pathoYes$Location, gsub("chr", "", ALLvars$Location)) 
sum(!is.na(match(pathoYes$Location, gsub("chr", "", ALLvars$Location))))
# pathoYes_Sum
match(pathoYes_Sum$Location, gsub("chr", "", ALLvars$Location)) 
sum(!is.na(match(pathoYes_Sum$Location, gsub("chr", "", ALLvars$Location))))

# ###ADFTDP_possible
# map with gene and codon
key1 <- paste(ADFTDPDmutDB_possible$GENE, ADFTDPDmutDB_possible$g.position, sep =  ":")
key2 <- paste(ALLvars$SYMBOL, gsub("c.", "g.", ALLvars$c.change), sep = ":")
sum(!is.na(match(key1, key2)))

# map with gene and amino acid change
key1 <- paste(ADFTDPDmutDB$GENE, ADFTDPDmutDB$p.CHANGE, sep =  ":")
key2 <- paste(ALLvars$SYMBOL, ALLvars$p.change, sep = ":")
sum(!is.na(match(key1, key2)))

# ###ADFTDPD_patho
# match location
match(ADFTDPDmutDB_patho$`CHR:POS:REF:ALT(GRCh37)`, gsub("chr", "", ALLvars$Location)) 
sum(!is.na(match(ADFTDPDmutDB_patho$`CHR:POS:REF:ALT(GRCh37)`, gsub("chr", "", ALLvars$Location))))
# map with gene and codon
key1 <- paste(ADFTDPDmutDB_patho$GENE, ADFTDPDmutDB_patho$cCHANGe, sep =  ":")
key2 <- paste(ALLvars$SYMBOL, ALLvars$c.change, sep = ":")
sum(!is.na(match(key1, key2)))

# map with gene and amino acid change
key1 <- paste(ADFTDPDmutDB_patho$GENE, ADFTDPDmutDB_patho$p.CHANGE, sep =  ":")
key2 <- paste(ALLvars$SYMBOL, ALLvars$p.change, sep = ":")
sum(!is.na(match(key2, key1)))



####### ######### match by RsID
dtRSiD <-  gsub("rs", "", ALLvars$rsID)
dtRSiD <- dtRSiD[!is.na(dtRSiD) ]
match(as.character(pathoYes$`rs#`), dtRSiD)
sum(!is.na(match(as.character(pathoYes$`rs#`), dtRSiD)))


FINAL <- ALLvars[,-c(10:11)]
library(dplyr)
FINAL <- distinct(FINAL)
dim(FINAL)

write.table(FINAL, "WES_WGS_AD_FTD_PD_all_variants.tsv", col.names = TRUE, quote = FALSE, sep = "\t")
 
# #######################
# # aggregate(ALLvars$SYMBOL, by = list( ALLvars$Both_WGS_and_WES, ALLvars$data), FUN=table)
# # 
# # aggregate(ALLvars$SYMBOL, by = list( ALLvars$Both_WGS_and_WES), FUN=table)
# # 
# # 
# # table(ALLvars$SYMBOL, group_by(ALLvars$data)
# # table(ALLvars$SYMBOL, ALLvars$data, ALLvars$Both_WGS_and_WES)
# # sum(WES$Location %in% WGS$Location)
# 
# 
# #######################################################################################
# # tmp.vep <- head(tmp_vcf_data, 100)
# 
# 
# 
# # sum(grepl ("non_coding", tmp_vcf_dataG$Consequence))
# # 
# # final.vep <- tmp_vcf_dataG[!grepl("non_synonymous", tmp_vcf_dataG$Consequence),]
# # dim(final.vep)
# # sum(grepl("coding", final.vep$Consequence))
# # final.coding <- final.vep[!grepl("non_coding", final.vep$Consequence),]
# # 
# # 
# # 
# # sum(final.vep$SYMBOL %in% mendelian)
# # 
# # mendelian.vars <- final.vep[final.vep$SYMBOL %in% mendelian, ]
# # dim(mendelian.vars)
# # # non synonymous 
# # mendelian.vars.non_syn <- mendelian.vars[mendelian.vars$Consequence != "synonymous_variant", ]
# # dim(mendelian.vars.non_syn)
# # mendelian.vars.non_syn_coding <- mendelian.vars.non_syn[mendelian.vars.non_syn$Consequence != "intron_variant", ]
# # dim(mendelian.vars.non_syn_coding)
# 
# 
# exclude <-  c("downstream_gene_variant", "3_prime_UTR_variant", "5_prime_UTR_variant",
#               "upstream_gene_variant", "splice_region_variant,intron_variant",
#               "splice_region_variant,synonymous_variant", "splice_region_variant,5_prime_UTR_variant")
# 
# # mendelian_codin_non_synonymous <- mendelian.vars.non_syn_coding[(!mendelian.vars.non_syn_coding$Consequence %in% exclude),]
# 
# include <- c("missense_variants", "stop_gained", "coding_sequence_variant", "inframe_insertion", "inframe_deletion", "frameshift_variant")
# mendelian_codin_non_synonymous_include <- tmp_vcf_dataG[(tmp_vcf_dataG$Consequence %in% include),]
# mendelian_codin_non_synonymous <- tmp_vcf_dataG[(!tmp_vcf_dataG$Consequence %in% exclude),]
# 
# 
# 
# write.table(mendelian_codin_non_synonymous, "WGS_AD_FTD_PD_non_synonymous_coding_variants.tsv", col.names = TRUE, quote = FALSE, sep = "\t")
# write.table(mendelian.vars, "WGS_AD_FTD_PD_all_variants.tsv", col.names = TRUE, quote = FALSE, sep = "\t")
# 

######################################################################################

df <- read.table("//fenix.psych.wucon.wustl.edu/achal/carlos_task/AGRF/genome/WES_WGS_AD_FTD_PD_all_variants.tsv", header = TRUE)
head(df)

###




















#####################################################################################################
# mendelian genes
mend.table <- read.table(text = "Disease	Gene	N	OR	P	N	OR	P	N	OR	P	N	OR	P
AD	APP	19	1.522	0.340	93	0.958	0.396	15	3.106	0.103	31	1.122	0.688
AD	PSEN1	10	3.380	0.122	51	1.656	0.990	10	3.380	0.122	23	1.520	0.880
AD	PSEN2	24	1.011	1.000	74	0.897	0.123	11	1.471	0.763	18	0.810	0.414
AD	PRNP	3	1.259	0.762	28	1.001	0.534	2	0.838	1.000	8	5.677	0.991
FTD	CHMP2B	7	2.238	0.059	25	0.878	0.243	4	NA	0.130	10	0.202	0.026
FTD	FUS	12	0.552	0.159	61	0.976	0.463	9	0.416	0.315	18	0.648	0.246
FTD	GRN	26	0.656	0.119	94	1.044	0.679	19	1.444	0.494	12	4.056	0.991
FTD	MAPT	33	1.101	0.791	97	0.905	0.178	23	0.639	0.299	37	1.064	0.635
FTD	TARDBP	4	1.118	1.000	19	0.867	0.197	2	NA	0.503	8	0.810	0.517
FTD	TBK1	18	0.956	1.000	77	1.207	0.902	13	2.820	0.160	22	1.419	0.843
FTD	VCP	15	0.727	0.459	41	1.071	0.701	12	0.595	0.398	9	6.489	0.995
PD	LRRK2	54	1.142	0.591	267	0.994	0.488	38	0.705	0.318	43	1.461	0.910
PD	PARK2	20	0.831	0.577	82	0.971	0.407	9	1.048	1.000	31	1.284	0.804
PD	PARK7	8	1.471	0.763	30	0.916	0.389	6	1.680	0.694	18	1.274	0.768
PD	PINK1	29	0.935	0.866	91	0.893	0.212	23	0.835	0.673	27	1.105	0.671
PD	SNCA	3	1.399	0.734	11	0.866	0.283	2	NA	0.503	-	-	-
PD	UCHL1	6	1.682	0.521	31	1.045	0.621	5	3.365	0.384	11	2.162	0.933
PD	ATP13A2	51	0.750	0.149	195	0.837	0.006	32	1.018	1.000	79	0.897	0.357
PD	GIGYF2	30	0.769	0.402	165	1.040	0.692	22	0.623	0.377	70	1.146	0.752
PD	HTRA2	11	1.231	0.595	50	1.196	0.909	6	0.166	0.098	25	0.636	0.176
PD	PLA2G6	40	0.976	1.000	136	0.999	0.510	26	0.977	1.000	17	0.911	0.518
PD	FBXO7	20	1.452	0.322	71	1.030	0.611	17	1.200	0.809	42	1.035	0.603
PD	VPS35	7	6.771	0.045	46	0.913	0.361	5	3.365	0.384	30	0.868	0.421
PD	EIF4G1	40	1.713	0.010	204	1.016	0.606	26	1.146	0.843	100	0.954	0.447
PD	DNAJC16	29	1.404	0.186	143	0.918	0.227	20	1.567	0.374	69	0.833	0.262
ALS	SOD1	1	0.000	0.208	8	0.810	0.327	-	-	-	1	NA	1.000
ALS	OPTN	14	0.835	0.684	68	1.092	0.750	12	0.836	0.779	30	1.060	0.631
ALS	UBQLN2	5	1.688	0.346	-	-	-	-	-	-	-	-	-
ALS	PFN1	1	0.000	0.456	18	1.128	0.716	1	0.000	0.456	3	1.621	0.831
ALS	SQSTM1	30	0.758	0.271	99	0.945	0.295	22	0.694	0.399	24	0.957	0.537
TOTAL AD	56	1.174	0.439	246	0.988	0.443	38	2.315	0.028	80	1.283	0.885
TOTAL FTD	115	0.877	0.378	414	0.980	0.367	82	1.057	0.907	116	1.092	0.711
TOTAL PD	348	1.049	0.688	1522	0.976	0.276	237	0.897	0.492	562	1.013	0.573
TOTAL ALS	51	0.823	0.367	193	0.980	0.410	35	0.698	0.307	58	1.070	0.649", header = TRUE)

mend.table$Gene <- as.character(mend.table$Gene)

mendelian <- mend.table$Gene[!mend.table$Disease %in% c("ALS", "TOTAL")]
mendelian <- c(mendelian, "C9ORF72")







