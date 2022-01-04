setwd("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/")
# Read all Covar files
TARCC3_AA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/TARCC3_AA/Covariates/tarcc3_AA_covar_06112020.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(TARCC3_AA)
head(TARCC3_AA)

NIALOAD_NCRAD_AA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/NIALOAD_NCRAD_AA/Covariates/niaload.covar.2017Dec15.txt", header = T, sep = " ", stringsAsFactors = FALSE)
dim(NIALOAD_NCRAD_AA)
colnames(NIALOAD_NCRAD_AA) <- c("FID", "IID", "missing",  "SEX", "age", "APOE", "STATUS", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "p10", "samp_order")
head(NIALOAD_NCRAD_AA)

MIRAGE600_AA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/MIRAGE600_AA/Covariates/mir600.covar.2017Dec19.phen.csv", header = T, sep = ",", stringsAsFactors = FALSE)
dim(MIRAGE600_AA)
head(MIRAGE600_AA)
MIRAGE600_AA$FID <- MIRAGE600_AA$id
colnames(MIRAGE600_AA) <- c("IID", "STATUS", "SEX", "age", "APOE", "pc1", "pc2", "pc3", "FID")
head(MIRAGE600_AA)

MIRAGE300_AA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/MIRAGE300_AA/Covariates/mir300.covar.2017Dec19.phen.csv", header = T, sep = ",", stringsAsFactors = FALSE)
dim(MIRAGE300_AA)
head(MIRAGE300_AA)
MIRAGE300_AA$FID <- MIRAGE300_AA$id
colnames(MIRAGE300_AA) <- c("IID", "STATUS", "SEX", "age", "APOE", "pc1", "pc2", "pc3", "FID")
head(MIRAGE300_AA)

JHU_AA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/JHU_AA/Covariates/jhu.covar.2017Dec15.txt", header = T, sep = " ", stringsAsFactors = FALSE)
dim(JHU_AA)
head(JHU_AA)
colnames(JHU_AA) <- c("FID", "IID", "missing", "SEX", "age", "APOE", "STATUS", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "p10", "samp_order")
head(JHU_AA)

INDIANAPOLIS_AA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/INDIANAPOLIS_AA/Covariates/indianapolis.covar.2017Dec15.txt", header = T, sep = " ", stringsAsFactors = FALSE)
dim(INDIANAPOLIS_AA)
head(INDIANAPOLIS_AA)
colnames(INDIANAPOLIS_AA) <- c("FID", "IID", "missing", "SEX", "age", "APOE", "STATUS", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "p10", "samp_order")
head(INDIANAPOLIS_AA)

CHOP_AA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/CHOP_AA/Covariates/adgcchop.covar.2017Dec15.txt", header = T, sep = " ", stringsAsFactors = FALSE)
dim(CHOP_AA)
head(CHOP_AA)
colnames(CHOP_AA) <- c("FID", "IID", "missing", "SEX", "age", "APOE", "STATUS", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "p10", "samp_order")
head(CHOP_AA)

CHAP_AA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/CHAP_AA/Covariates/chap.covar.2017Dec15.txt", header = T, sep = " ", stringsAsFactors = FALSE)
dim(CHAP_AA)
head(CHAP_AA)
colnames(CHAP_AA) <- c("FID", "IID", "missing", "SEX", "age", "APOE", "STATUS", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "p10", "samp_order")
head(CHAP_AA)

ADC9_AA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/ADC9_AA/Covariates/adc9_AA_covariate_15Nov2019.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC9_AA)
head(ADC9_AA)
colnames(ADC9_AA) <- c("FID", "IID", "cohort", "old_fid", "old_iid", "merged_id", "omit", "SEX",
                       "STATUS", "RACE", "HISPANIC", "Age_at_Onset", "Age_at_LastExam", "Age_at_Death",
                       "aaoaae", "aaoaae2", "apoe_1", "apoe_2", "apoe4any", "apoe4dose", "pc1", "pc2",
                        "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10", "gen_order")
head(ADC9_AA)

ADC8_AA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/ADC8_AA/Covariates/adc8_aa_covar_19Aug2019.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC8_AA)
head(ADC8_AA)
colnames(ADC8_AA) <- c("FID", "IID", "cohort", "old_fid", "old_iid", "merged_id", "omit",
                       "SEX", "STATUS", "RACE", "Ethnicity", "Age_at_Onset", "Age_at_LastExam",
                        "Age_at_Death", "APOE_A", "APOE_B", "aaoaae", "aaoaae2", "apoe_1", "apoe_2",
                        "apoe4any", "apoe4dose", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8",
                        "pc9", "pc10", "gen_order")

wanted.cols <-  c("FID", "IID", "cohort", "old_fid", "old_iid", "merged_id", "omit",
                       "SEX", "STATUS", "RACE", "Age_at_Onset", "Age_at_LastExam", 
                       "Age_at_Death", "aaoaae", "aaoaae2", "apoe_1", "apoe_2",
                       "apoe4any", "apoe4dose", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8",
                       "pc9", "pc10", "gen_order")  

ADC8_AA <- ADC8_AA[wanted.cols]

head(ADC8_AA)

ADC3_AA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/ADC3_AA/Covariates/adc3.covar.2017Dec15.txt", header = T, sep = " ", stringsAsFactors = FALSE)
dim(ADC3_AA)
head(ADC3_AA)
colnames(ADC3_AA) <- c("FID", "IID", "missing", "SEX", "age", "APOE", "STATUS", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "p10", "samp_order")
head(ADC3_AA)

ADC12_AA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/ADC12_AA/Covariates/adc12_aframr.covariate_16Dec_2020.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC12_AA)
head(ADC12_AA)
colnames(ADC12_AA) <- c("FID", "IID", "cohort", "old_fid", "old_iid", "merged_id", "omit", "SEX", "STATUS",
                        "Age_at_Onset", "Age_at_Death", "Age_at_LastExam", "aaoaae_orig", "aaoaae", "aaoaae2",
                        "apoe_1", "apoe_2", "apoe4any", "apoe4dose", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6",
                        "pc7", "pc8", "pc9", "pc10", "gen_order", "omit_reason")  
head(ADC12_AA)

ADC11_AA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/ADC11_AA/Covariates/adc11_african_american.covariate_01Oct_2020.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC11_AA)
head(ADC11_AA)
colnames(ADC11_AA) <- c("FID", "IID", "cohort", "old_fid", "old_iid", "merged_id", "omit",
                        "SEX", "STATUS", "Age_at_Onset", "Age_at_Death", "Age_at_LastExam",
                        "aaoaae_orig", "aaoaae", "aaoaae2", "apoe_1", "apoe_2", "apoe4any",
                        "apoe4dose", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8",
                        "pc9", "pc10", "gen_order")  
head(ADC11_AA)


ADC10_AA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/ADC10_AA/Covariates/adc10_AA_covariate_15Nov2019.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC10_AA)
head(ADC10_AA)
colnames(ADC10_AA) <- c("FID", "IID", "cohort", "old_fid", "old_iid", "merged_id", "omit",
                        "SEX", "STATUS", "RACE", "HISPANIC", "Age_at_Onset", "Age_at_LastExam",
                        "Age_at_Death", "aaoaae", "aaoaae2", "apoe_1", "apoe_2", "apoe4any", 
                        "apoe4dose", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8",
                        "pc9", "pc10", "gen_order")  
head(ADC10_AA)


ADC1_2_AA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/ADC1-2_AA/Covariates/adc12.covar.2017Dec15.txt", header = T, sep = " ", stringsAsFactors = FALSE)
dim(ADC1_2_AA)
head(ADC1_2_AA)
colnames(ADC1_2_AA) <- c("FID", "IID", "missing", "SEX", "age", "APOE", "STATUS", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "p10", "samp_order")
head(ADC1_2_AA)

ACT_AA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/ACT_AA/Covariates/act.covar.2017Dec15.txt", header = T, sep = " ", stringsAsFactors = FALSE)
dim(ACT_AA)
head(ACT_AA)
colnames(ACT_AA) <- c("FID", "IID", "missing", "SEX", "age", "APOE", "STATUS", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "p10", "samp_order")
head(ACT_AA)

ACT3_AA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/ACT3_AA/Covariates/ACT3_afr_amr.covariate_28Jan_2021.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ACT3_AA)
head(ACT3_AA)
colnames(ACT3_AA) <-  c("FID", "IID", "cohort", "old_fid", "old_iid", "merged_id", "omit",
                      "SEX", "STATUS", "Age_at_Onset", "Age_at_Death", "Age_at_LastExam", 
                      "aaoaae_orig", "aaoaae", "aaoaae2", "apoe_1", "apoe_2", "apoe4any", 
                      "apoe4dose", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8",
                      "pc9", "pc10", "gen_order", "omit_reason")  

head(ACT3_AA)

all.dfs <- scan(text = "TARCC3_AA NIALOAD_NCRAD_AA MIRAGE600_AA MIRAGE300_AA JHU_AA INDIANAPOLIS_AA CHOP_AA CHAP_AA ADC9_AA ADC8_AA ADC3_AA ADC12_AA ADC11_AA ADC10_AA ADC1_2_AA ACT_AA ACT3_AA",  what="")
all.dfs
for(i in seq_along(all.dfs)){
assign(all.dfs[i], `[[<-`(get(all.dfs[i]), 'cohort', value=all.dfs[i])) #creating a new column
}

AFRICAN <- Reduce(function(x, y) merge(x, y, all=TRUE), list(TARCC3_AA, NIALOAD_NCRAD_AA, MIRAGE600_AA, MIRAGE300_AA, JHU_AA, INDIANAPOLIS_AA, CHOP_AA, CHAP_AA, ADC9_AA, ADC8_AA, ADC3_AA, ADC12_AA, ADC11_AA, ADC10_AA, ADC1_2_AA, ACT_AA, ACT3_AA))
wanted.cols <- c("FID", "IID", "SEX", "STATUS", "cohort",  "age", "Age_at_Onset", "Age_at_LastExam", "Age_at_Death", "old_fid", "old_iid", "merged_id", "omit", "aaoaae", "aaoaae2", "apoe_1", "apoe_2", "apoe4any", "apoe4dose", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "p10", "gen_order", "aaoaae_orig", "omit_reason", "missing", "samp_order", "APOE", "RACE", "HISPANIC", "DX", "age_baseline")
# wanted.cols <- c("FID", "IID", "SEX", "STATUS", "cohort",  "age", "Age_at_Onset", "Age_at_LastExam", "Age_at_Death", "old_fid", "old_iid", "merged_id", "omit", "aaoaae", "aaoaae2", "apoe_1", "apoe_2", "apoe4any", "apoe4dose", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "p10", "gen_order", "aaoaae_orig", "omit_reason", "missing", "samp_order", "APOE", "RACE", "HISPANIC", "DX", "age_baseline", "APOE_A", "APOE_B")
colnames(AFRICAN)
AFRICAN <- AFRICAN[wanted.cols]

View(AFRICAN)
write.table(AFRICAN, "African_American_covariate.tsv", sep ="\t", col.names = T, quote = F, row.names = FALSE)
############## Hispanic

TARCC3_Hispanic <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Hispanic/TARCC3_Hispanic/Covariates/tarcc3_hispanic_covar_06112020.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(TARCC3_Hispanic)
head(TARCC3_Hispanic)
## NO STATUS for TARCC3_Hispanic

ADC9_Hispanic <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Hispanic/ADC9_Hispanic/Covariates/adc9_hispanic_covariate_15Nov2019.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC9_Hispanic)
head(ADC9_Hispanic)

ADC8_Hispanic <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Hispanic/ADC8_Hispanic/Covariates/adc8_hispanic_covar_19Aug2019.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC8_Hispanic)
head(ADC8_Hispanic)

ADC12_Hispanic <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Hispanic/ADC12_Hispanic/Covariates/adc12_hispanic.covariate_16Dec_2020.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC12_Hispanic)
head(ADC12_Hispanic)
colnames(ADC12_Hispanic)[colnames(ADC12_Hispanic)=="sex"] <- "SEX"
colnames(ADC12_Hispanic) [grep("onset|death|visit", tolower(colnames(ADC12_Hispanic)))] <- c("Age_at_Onset", "Age_at_Death", "Age_at_LastExam")
head(ADC12_Hispanic)

ADC11_Hispanic <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Hispanic/ADC11_Hispanic/Covariates/adc11_hispanic.covariate_01Oct_2020.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC11_Hispanic)
head(ADC11_Hispanic)
colnames(ADC11_Hispanic)[colnames(ADC11_Hispanic)=="sex"] <- "SEX"
colnames(ADC11_Hispanic) [grep("onset|death|visit", tolower(colnames(ADC11_Hispanic)))] <- c("Age_at_Onset", "Age_at_Death", "Age_at_LastExam")
head(ADC11_Hispanic)


ADC10_Hispanic <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Hispanic/ADC10_Hispanic/Covariates/adc10_hispanic_covariate_15Nov2019.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC10_Hispanic)
head(ADC10_Hispanic)

ACT3_Hispanic <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Hispanic/ACT3_Hispanic/Covariates/ACT3_hispanic.covariate_28Jan_2021.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ACT3_Hispanic)
head(ACT3_Hispanic)
colnames(ACT3_Hispanic)[colnames(ACT3_Hispanic)=="sex"] <- "SEX"
colnames(ACT3_Hispanic) [grep("onset|death|visit", tolower(colnames(ACT3_Hispanic)))] <- c("Age_at_Onset", "Age_at_Death", "Age_at_LastExam")
head(ACT3_Hispanic)

all.dfs <- scan(text = "TARCC3_Hispanic ADC9_Hispanic ADC8_Hispanic ADC12_Hispanic ADC11_Hispanic ADC10_Hispanic ACT3_Hispanic",  what="")
all.dfs
for(i in seq_along(all.dfs)){
  assign(all.dfs[i], `[[<-`(get(all.dfs[i]), 'cohort', value=all.dfs[i])) #creating a new column
}


HISPANIC <- Reduce(function(x, y) merge(x, y, all=TRUE), list(TARCC3_Hispanic, ADC9_Hispanic, ADC8_Hispanic, ADC12_Hispanic, ADC11_Hispanic, ADC10_Hispanic, ACT3_Hispanic))
View(HISPANIC)
write.table(HISPANIC, "HISPANIC_covariate.tsv", sep ="\t", col.names = T, quote = F, row.names = FALSE)
#######Asian

JPN <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Asian/JPN/Covariates/jpn2012_asian.covariate_20Jan_2021.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(JPN)
head(JPN)
colnames(JPN)[colnames(JPN)=="sex"] <- "SEX"
colnames(JPN) [grep("onset|death|visit", tolower(colnames(JPN)))] <- c("Age_at_Onset", "Age_at_Death", "Age_at_LastExam")
head(JPN)

ASA_JPN <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Asian/ASA-JPN/Covariates/asajpn_asian.covariate_2Feb_2021.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ASA_JPN)
head(ASA_JPN)
colnames(ASA_JPN)[colnames(ASA_JPN)=="sex"] <- "SEX"
colnames(ASA_JPN) [grep("onset|death|visit", tolower(colnames(ASA_JPN)))] <- c("Age_at_Onset", "Age_at_Death", "Age_at_LastExam")
head(ASA_JPN)

ADC12_Asian <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Asian/ADC12_Asian/Covariates/adc12_asian.covariate_16Dec_2020.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC12_Asian)
head(ADC12_Asian)
colnames(ADC12_Asian)[colnames(ADC12_Asian)=="sex"] <- "SEX"
colnames(ADC12_Asian) [grep("onset|death|visit", tolower(colnames(ADC12_Asian)))] <- c("Age_at_Onset", "Age_at_Death", "Age_at_LastExam")
head(ADC12_Asian)

ADC11_Asian <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Asian/ADC11_Asian/Covariates/adc11_asian.covariate_01Oct_2020.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC11_Asian)
head(ADC11_Asian)
colnames(ADC11_Asian)[colnames(ADC11_Asian)=="sex"] <- "SEX"
colnames(ADC11_Asian) [grep("onset|death|visit", tolower(colnames(ADC11_Asian)))] <- c("Age_at_Onset", "Age_at_Death", "Age_at_LastExam")
head(ADC11_Asian)


ADC10_Asian <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Asian/ADC10_Asian/Covariates/adc10_asian_covariate_11Aug2020.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC10_Asian)
head(ADC10_Asian)

ACT3_Asian <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Asian/ACT3_Asian/Covariates/ACT3_asian.covariate_28Jan_2021.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ACT3_Asian)
head(ACT3_Asian)
colnames(ACT3_Asian)[colnames(ACT3_Asian)=="sex"] <- "SEX"
colnames(ACT3_Asian) [grep("onset|death|visit", tolower(colnames(ACT3_Asian)))] <- c("Age_at_Onset", "Age_at_Death", "Age_at_LastExam")
head(ACT3_Asian)

all.dfs <- scan(text = "JPN ASA_JPN ADC12_Asian ADC11_Asian ADC10_Asian ACT3_Asian",  what="")
all.dfs
for(i in seq_along(all.dfs)){
  assign(all.dfs[i], `[[<-`(get(all.dfs[i]), 'cohort', value=all.dfs[i])) #creating a new column
}


ASIAN <- Reduce(function(x, y) merge(x, y, all=TRUE), list(JPN, ASA_JPN, ADC12_Asian, ADC11_Asian, ADC10_Asian, ACT3_Asian))

write.table(ASIAN, "ASIAN_covariate.tsv", sep ="\t", col.names = T, quote = F, row.names = FALSE)
########## NHW

ACT2 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ACT2/Covariates/act2.covar.2018Nov05.txt", header = T, sep = " ", stringsAsFactors = FALSE)
# ACT2 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ACT2/Covariates/act2.covar.for_snptest.txt", header = T, sep = " ")
# ACT2 <- ACT2[-1,]
dim(ACT2)
colnames(ACT2) [grep("onset", tolower(colnames(ACT2)))] <- "age_onset"
colnames(ACT2) [grep("lastexam", tolower(colnames(ACT2)))] <- "age_last_visit"
colnames(ACT2) [grep("death", tolower(colnames(ACT2)))] <- "age_death"
colnames(ACT2) [grep("SEX", colnames(ACT2))] <- "sex"
colnames(ACT2) [grep("apoe1", colnames(ACT2))] <- "apoe_1"
colnames(ACT2) [grep("race", colnames(ACT2))] <- "RACE"
colnames(ACT2) [grep("hispanic", colnames(ACT2))] <- "HISPANIC"
head(ACT2)

ACT1 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ACT1/Covariates/act.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ACT1)
head(ACT1)

WHICAP <- read.csv("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/WHICAP/Covariates/whicap.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
# WHICAP <- read.csv("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/WHICAP/Covariates/whicap.covar.for_snptest.txt", header = T, sep = " ")
# WHICAP <- WHICAP[-1,]
dim(WHICAP)
head(WHICAP)

WASHU2 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/WASHU2/Covariates/washu2.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(WASHU2)
head(WASHU2)

WASHU1 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/WASHU1/Covariates/washu.covar.2014Mar17_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(WASHU1)
head(WASHU1)

UPITT <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/UPITT/Covariates/upitt.covar.2014Mar17_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(UPITT)
head(UPITT)

UMVUTARC2 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/UMVUTARC2/Covariates/mtv.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(UMVUTARC2)
head(UMVUTARC2)

UMVUMSSM <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/UMVUMSSM/Covariates/umvumssm_combined.covar.2013Dec06.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(UMVUMSSM)
head(UMVUMSSM)

UKS <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/UKS/Covariates/uks.covar.2014Jun27_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(UKS)
head(UKS)

TGEN2 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/TGEN2/Covariates/tgen2.covar.2014Mar17_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(TGEN2)
head(TGEN2)

# /40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/TARCC3_NHW/**********MISSING_COVAR 
TARCC3 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/TARCC3/Covariates/tarcc3_caucasian_covar_06112020.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(TARCC3)
colnames(TARCC3) [grep("onset", tolower(colnames(TARCC3)))] <- "age_onset"
colnames(TARCC3) [grep("lastexam", tolower(colnames(TARCC3)))] <- "age_last_visit"
colnames(TARCC3) [grep("death", tolower(colnames(TARCC3)))] <- "age_death"
colnames(TARCC3) [grep("SEX", colnames(TARCC3))] <- "sex"

head(TARCC3)

TARCC1 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/TARCC1/Covariates/tarc1.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(TARCC1)
head(TARCC1)

ROSMAP2 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ROSMAP2/Covariates/rosmap2.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
# ROSMAP2 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ROSMAP2/Covariates/rosmap2.covar.for_snptest_AN.txt", header = T, sep = "\t")
# ROSMAP2 <- ROSMAP2[-1,]
dim(ROSMAP2)
head(ROSMAP2)

ROSMAP1 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ROSMAP1/Covariates/rosmap.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
# ROSMAP1 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ROSMAP1/Covariates/rosmap1.covar.for_snptest_AN.txt", header = T, sep = "\t")
# ROSMAP1 <- ROSMAP1[-1,]
dim(ROSMAP1)
head(ROSMAP1)

RMayo <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/RMayo/Covariates/rmayo.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(RMayo)
head(RMayo)

OHSU <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/OHSU/Covariates/ohsu.covar.2014Mar17_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(OHSU)
head(OHSU)

NBB <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/NBB/Covariates/nbb.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(NBB)
head(NBB)
colnames(NBB)[1:2] <- c("FID", "IID")
head(NBB)

MIRAGE <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/MIRAGE/Covariates/mirage.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(MIRAGE)
head(MIRAGE)

MAYO <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/MAYO/Covariates/mayo.covar.2014Mar17_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(MAYO)
head(MAYO)


LOAD <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/LOAD/Covariates/load.covar.2014Mar17_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(LOAD)
head(LOAD)

GSK <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/GSK/Covariates/gsk.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(GSK)
head(GSK)

EAS <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/EAS/Covariates/eas.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
# EAS <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/EAS/Covariates/eas.covar.for_snptest_AN.txt", header = T, sep = "\t")
# EAS <- EAS[-1,]
dim(EAS)
head(EAS)

CHAP2 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/CHAP2/Covariates/chap.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(CHAP2)
head(CHAP2)

BIOCARD <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/BIOCARD/Covariates/biocard.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(BIOCARD)
head(BIOCARD)


ADNI <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ADNI/Covariates/adni.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADNI)
head(ADNI)

ADC9 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ADC9/Covariates/adc9_caucasian_covariate_15Nov2019v2.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC9)
colnames(ADC9)[grep("Race", colnames(ADC9))] <- "RACE"
head(ADC9)
ADC9 <- ADC9[!grepl("APOE_A|APOEB",colnames(ADC9))]
head(ADC9)

ADC8 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ADC8/Covariates/adc8_caucasian_covar_11March2019.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC8)
head(ADC8)

ADC7 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ADC7/Covariates/adc7.covar.2015May05_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC7)
head(ADC7)

ADC6 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ADC6/Covariates/adc6.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC6)
head(ADC6)

ADC5 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ADC5/Covariates/adc5.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC5)
head(ADC5)

ADC4 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ADC4/Covariates/adc4.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC4)
head(ADC4)

ADC3 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ADC3/Covariates/adc3.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC3)
head(ADC3)

ADC2 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ADC2/Covariates/adc2.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC2)
head(ADC2)

ADC12 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ADC12/Covariates/adc12_nonhisp_caucasian.covariate_16Dec_2020.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC12)
head(ADC12)

ADC11 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ADC11/Covariates/adc11_nonhispanic_caucasian.covariate_01Oct_2020.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC11)
head(ADC11)

ADC10 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ADC10/Covariates/adc10_caucasian_covariate_11Aug2020.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC10)
colnames(ADC10) [grep("onset", tolower(colnames(ADC10)))] <- "age_onset"
colnames(ADC10) [grep("lastexam", tolower(colnames(ADC10)))] <- "age_last_visit"
colnames(ADC10) [grep("death", tolower(colnames(ADC10)))] <- "age_death"
colnames(ADC10) [grep("SEX", colnames(ADC10))] <- "sex"
colnames(ADC10)[grep("race", colnames(ADC10))] <- "RACE"
head(ADC10)

ADC1 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ADC1/Covariates/adc1.covar.2013Dec06_AN.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ADC1)
head(ADC1)

ACT3 <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ACT3/Covariates/ACT3_nonhisp_caucasian.covariate_28Jan_2021.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
dim(ACT3)
head(ACT3)

# Check the number of samples in NHW covariates and in VCF 
alldf <- c("ACT1", "ACT2", "ACT3", "ADC1", "ADC10", "ADC11", "ADC12",
           "ADC2", "ADC3", "ADC4", "ADC5", "ADC6", "ADC7", "ADC8", "ADC9", "ADNI",
           "BIOCARD", "CHAP2", "EAS", "GSK", "LOAD", "MAYO", "MIRAGE", "NBB", "OHSU",
           "RMayo", "ROSMAP1", "ROSMAP2", "TARCC1", "TARCC3", "TGEN2", "UKS", "UMVUMSSM",
           "UMVUTARC2", "UPITT", "WASHU1", "WASHU2", "WHICAP")

for (i in 1:length(alldf)){
  print(paste(alldf[i], (nrow(eval(parse(text=alldf[i])))), sep =":::"  ))
}



all.dfs <- scan(text = "ACT2 ACT1 WHICAP WASHU2 WASHU1 UPITT UMVUTARC2 UMVUMSSM
            UKS TGEN2 TARCC3 TARCC1 ROSMAP2 ROSMAP1 RMayo OHSU NBB 
            MIRAGE MAYO LOAD GSK EAS CHAP2 BIOCARD ADNI ADC9 ADC8 
            ADC7 ADC6 ADC5 ADC4 ADC3 ADC2 ADC12 ADC11 ADC10 ADC1 ACT3",  what="")
all.dfs
for(i in seq_along(all.dfs)){
  assign(all.dfs[i], `[[<-`(get(all.dfs[i]), 'cohort', value=all.dfs[i])) #creating a new column
}



NHW <- Reduce(function(x, y) merge(x, y, all=TRUE), list(ACT2, ACT1, WHICAP, WASHU2, WASHU1, UPITT, UMVUTARC2, UMVUMSSM,
                                                         UKS, TGEN2, TARCC3, TARCC1, ROSMAP2, ROSMAP1, RMayo, OHSU, NBB, 
                                                         MIRAGE, MAYO, LOAD, GSK, EAS, CHAP2, BIOCARD, ADNI, ADC9, ADC8, 
                                                         ADC7, ADC6, ADC5, ADC4, ADC3, ADC2, ADC12, ADC11, ADC10, ADC1, ACT3))


major.columns <- c("FID", "IID", "sex", "status", "age_onset", "age_death", "age_last_visit", "aaoaae", "aaoaae2", "apoe_1", "apoe_2",
                    "apoe4any", "apoe4dose", "omit", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10", "cohort",
                    "merged_id", "gen_order", "old_fid", "old_iid", "aaoaae_orig", "omit_reason", "HISPANIC", "apoe",
                    "Sum_ID", "DX", "RACE", "age_baseline")

NHW <- NHW[major.columns]

write.table(NHW, "NHW_covariate.tsv", sep ="\t", col.names = T, quote = F, row.names = FALSE)


View(as.data.frame(unique(ASIAN$cohort)))

#################################################
# AFRICAN$ETNICITY <- "AFRICAN"
lst1 <- list(AFRICAN=AFRICAN, HISPANIC=HISPANIC, ASIAN=ASIAN, NHW=NHW)

# add cohortID to each dataframe
f <- function (data, ETHNICITY){
  data$ETHNICITY <- ETHNICITY
  data
}

lst1 <- Map(f, lst1 , names(lst1))

lapply (lst1, head)
lapply (lst1, dim)
lapply(lst1, colnames)



# make all column values upper case
lst2 <- lapply(lst1, function(x){
  colnames(x) <- toupper(colnames(x))
  x})


lapply (lst2, head)
lapply (lst2, dim)
lapply(lst2, colnames)

# check common columns to merge
common_columns <- Reduce(intersect, lapply(lst2, colnames))
common_columns

# all columns
all_variables <- unique(unlist(lapply(lst2, colnames)))

# columns that differ between tables
diff_columns <- setdiff(all_variables, common_columns)
diff_columns



lst2 <- lapply(lst2, function(x){
  colnames(x) [grep("LAST", colnames(x))] <- "AGE_LAST_VISIT"
  colnames(x) [grep("DEATH", colnames(x))] <- "AGE_AT_DEATH"
  colnames(x) [grep("ONSET", colnames(x))] <- "AGE_AT_ONSET"
  x})


common_columns <- Reduce(intersect, lapply(lst2, colnames))
common_columns

# all columns
all_variables <- unique(unlist(lapply(lst2, colnames)))

# columns that differ between tables
diff_columns <- setdiff(all_variables, common_columns)
diff_columns

covars <- Reduce(function(...) merge(..., all=TRUE), lst2)
# covars <- Reduce(function(x, y) merge(x, y, all = TRUE), lst2)
TOPMEDcovars <- covars


for(i in 1:nrow(covars)) {
  if(covars$STATUS[i]==1 & !is.na(covars$STATUS[i])) {
    covars$STATUS[i] = "CO"
  } else if (covars$STATUS[i]==2 & !is.na(covars$STATUS[i])) {
    covars$STATUS[i] = "CA"
  } else if (covars$STATUS[i]==3 & !is.na(covars$STATUS[i])) {
    covars$STATUS[i] = "MCI"
  } else {
    covars$STATUS[i] = "unknown"
  } 
}


# covars$STATUS[covars$STATUS==1] <- "CO"
# covars$STATUS[covars$STATUS==2] <- "CA"
# covars$STATUS[covars$STATUS==3] <- "MCI"
# covars$STATUS[covars$STATUS==-9] <- "unknown"
# covars$STATUS[is.na(covars$STATUS)] <- "unknown"
# # covars1 <- covars %>% 
# #   mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999, 0), "unknown")))

################################### 10/27/2020 Vicky only wants age of onset for Cases and MCI and Age at last assessment for CO
# CO > 75, 80 and 85

########################################################################################################
####################### Age-of-onset stratified by STATUS, ethnicity and sex ###########################
########################################################################################################

# table(covars[covars$AGE_AT_ONSET <65 & covars$AGE_AT_ONSET > 0, c("STATUS", "SEX")])
library(dplyr)
df <- covars %>%
  mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999, 0), "unknown"))) %>%
  group_by(STATUS, ETHNICITY) %>%
  summarise('<=65' = sum(as.numeric(as.character(AGE_AT_ONSET)) <= 65 & as.numeric(as.character(AGE_AT_ONSET)) > 0 ,
                        na.rm = TRUE),
            '<=70'= sum(as.numeric(as.character(AGE_AT_ONSET)) <= 70  & as.numeric(as.character(AGE_AT_ONSET)) > 0,
                       na.rm = TRUE),
            '<=75'= sum(as.numeric(as.character(AGE_AT_ONSET)) <= 75  & as.numeric(as.character(AGE_AT_ONSET)) > 0,
                       na.rm = TRUE))

# df <- covars %>% 
#   mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999, 0), "unknown"))) %>%
#   group_by(STATUS, ETHNICITY) %>%
#   summarise('<=65' = sum(as.numeric(as.character(AGE_AT_ONSET)) <= 65 & as.numeric(as.character(AGE_AT_ONSET)) > 0 , 
#                         na.rm = TRUE), 
#             '<=70'= sum(as.numeric(as.character(AGE_AT_ONSET)) <= 70  & as.numeric(as.character(AGE_AT_ONSET)) > 0, 
#                        na.rm = TRUE),
#             '<=75'= sum(as.numeric(as.character(AGE_AT_ONSET)) <= 75  & as.numeric(as.character(AGE_AT_ONSET)) > 0, 
#                        na.rm = TRUE),
#             'MISSING'= sum(is.na(AGE_AT_ONSET)))

df <- as.data.frame(df)
df
## require(tidyverse)
# df <- df %>% pivot_wider(Ethnicity, names_from = STATUS, values_from = c(`<65`,`<70`,`<75`))
library(reshape2)
# summarize values
# df <- melt(df, id.vars = 1:2, measure.vars = 3:6)
df <- reshape2::melt(df, id.vars = 1:2, measure.vars = 3:5)
df <- reshape2::dcast(df, ETHNICITY ~ STATUS + variable)
# colnames(df) <- gsub(".*_", "", colnames(df))
df
###############################################################
################# age at last assessment ######################
###############################################################

library(dplyr)
df2 <- covars %>%
  mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999), NA))) %>%
  group_by(STATUS, ETHNICITY) %>%
  summarise('>=70'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 70,
                       na.rm = TRUE),
            '>=80'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 80,
                       na.rm = TRUE),
            '>=85'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 85,
                       na.rm = TRUE))


# df2 <- covars %>% 
#   mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999), NA))) %>%
#   group_by(STATUS, ETHNICITY) %>%
#   summarise('>=70'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 70, 
#                        na.rm = TRUE),
#             '>=80'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 80, 
#                        na.rm = TRUE),
#             '>=85'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 85, 
#                        na.rm = TRUE),
#             'MISSING'= sum(is.na(AGE_LAST_VISIT)))

df2 <- as.data.frame(df2)
## require(tidyverse)
# df <- df %>% pivot_wider(Ethnicity, names_from = STATUS, values_from = c(`<65`,`<70`,`<75`))
library(reshape2)
# summarize values
# df2 <- melt(df2, id.vars = 1:2, measure.vars = 3:6)
df2 <- reshape2::melt(df2, id.vars = 1:2, measure.vars = 3:5)
df2 <- reshape2::dcast(df2, ETHNICITY ~ STATUS + variable)
# colnames(df) <- gsub(".*_", "", colnames(df))
df2

Age_stratified <- cbind(df,df2)
Age_stratified
write.csv(Age_stratified, "stratified_classification.csv", quote = FALSE)


# N Control Missing age
MISSING <- covars %>% 
  mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999), NA))) %>%
  group_by(STATUS, ETHNICITY) %>%
  summarise('MISSING'= sum(is.na(AGE_LAST_VISIT)))

MISSING <- reshape2::melt(MISSING, id.vars = 1:2, measure.vars = 3)
MISSING <- reshape2::dcast(MISSING, ETHNICITY ~ STATUS + variable)
MISSING

# N Cases Missing age
MISSING <- covars %>% 
  mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999), NA))) %>%
  group_by(STATUS, ETHNICITY) %>%
  summarise('MISSING'= sum(is.na(AGE_AT_ONSET)))

MISSING <- reshape2::melt(MISSING, id.vars = 1:2, measure.vars = 3)
MISSING <- reshape2::dcast(MISSING, ETHNICITY ~ STATUS + variable)
MISSING

table(APOE4=covars$APOE4ANY, STATUS=covars$STATUS)

## Status: control = 1, case = 2, MCI = 3
## Sex: 1 = male, 2 = female, -9 = unknown
df <- table(covars$SEX, covars$ETHNICITY)[2:3,]
df <- rbind(df,colSums(df))
df
# % Percent of females:
(df[2,]/df[3,])*100

## % APOE4
# df <- table(covars$APOE4ANY, covars$ETHNICITY)[2:3,]
# df <- rbind(df,colSums(df))
# df
# (df[2,]/df[3,])*100
(table( covars[, c("APOE4ANY", "ETHNICITY") ] )[3,]/ as.vector(aggregate(covars$IID, list(covars$ETHNICITY), length)[2])) * 100

# Number of cases Controls by Ethnicity
table(covars$STATUS, covars$ETHNICITY)

# % APOE4+ for each ethnicity
# table(covars$APOE4ANY, covars$STATUS, covars$ETHNICITY)
# table(STATUS=NHW$status==1, APOE=NHW$apoe4any==1 )

# % CONTROLS APOE4+
(table( covars[ covars$STATUS == "CO" , c("APOE4ANY", "ETHNICITY") ] )[3,]/table( covars[ covars$STATUS == "CO" , c("ETHNICITY") ] )) * 100

# % CASES APOE4+
(table( covars[ covars$STATUS == "CA" , c("APOE4ANY", "ETHNICITY") ] )[3,]/table( covars[ covars$STATUS == "CA" , c("ETHNICITY") ] )) * 100

# % MCI APOE4+
(table( covars[ covars$STATUS == "MCI" , c("APOE4ANY", "ETHNICITY") ] )[2,]/table( covars[ covars$STATUS == "MCI" , c("ETHNICITY") ] )) * 100


## APOE4
cbind.data.frame(table( covars[ , c("STATUS","APOE4ANY", "ETHNICITY") ] ))


######################################################################
######################ADGC Demographics by Ethnicity##################
######################################################################

covars <- Reduce(function(...) merge(..., all=TRUE), lst2)
# covars <- Reduce(function(x, y) merge(x, y, all = TRUE), lst2)
TOPMEDcovars <- covars


for(i in 1:nrow(covars)) {
  if(covars$STATUS[i]==1 & !is.na(covars$STATUS[i])) {
    covars$STATUS[i] = "CO"
  } else if (covars$STATUS[i]==2 & !is.na(covars$STATUS[i])) {
    covars$STATUS[i] = "CA"
  } else if (covars$STATUS[i]==3 & !is.na(covars$STATUS[i])) {
    covars$STATUS[i] = "MCI"
  } else {
    covars$STATUS[i] = "unknown"
  } 
}


# covars$STATUS[covars$STATUS==1] <- "CO"
# covars$STATUS[covars$STATUS==2] <- "CA"
# covars$STATUS[covars$STATUS==3] <- "MCI"
# covars$STATUS[covars$STATUS==-9] <- "unknown"
# covars$STATUS[is.na(covars$STATUS)] <- "unknown"
# # covars1 <- covars %>% 
# #   mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999, 0), "unknown")))

################################### 10/27/2020 Vicky only wants age of onset for Cases and MCI and Age at last assessment for CO
# CO > 75, 80 and 85

########################################################################################################
####################### Age-of-onset stratified by STATUS, ethnicity and sex ###########################
########################################################################################################


FIX_weird_NUMS <- function(x){
  x <- as.numeric(as.character(x))
  x [x %in% c(-9, 888, 999, 0)] <- NA
  return(x)
}

covars[grepl("AGE",colnames(covars))] <- sapply(covars[grepl("AGE",colnames(covars))], FIX_weird_NUMS)



# Number of CA <= 65 (Cases==2)
CA.65 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == "CA","AGE_AT_ONSET"])))) <= 65)
CA.70 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == "CA","AGE_AT_ONSET"])))) <= 70)




MEAN_CA <- length(as.numeric(as.character(covars[covars$STATUS == "CA","AGE_AT_ONSET"])), na.rm= T)
SD_CA <- sd(as.numeric(as.character(covars[covars$STATUS == 2,"AGE_AT_ONSET"])), na.rm= T)
RANGE_CA <- range(as.numeric(as.character(covars[covars$STATUS == 2,"AGE_AT_ONSET"])), na.rm= T)




# table(covars[covars$AGE_AT_ONSET <65 & covars$AGE_AT_ONSET > 0, c("STATUS", "SEX")])
library(dplyr)
df <- covars %>%
  mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999, 0), "unknown"))) %>%
  group_by(STATUS, ETHNICITY) %>%
  summarise('<=65' = sum(as.numeric(as.character(AGE_AT_ONSET)) <= 65 & as.numeric(as.character(AGE_AT_ONSET)) > 0 ,
                         na.rm = TRUE),
            '<=70'= sum(as.numeric(as.character(AGE_AT_ONSET)) <= 70  & as.numeric(as.character(AGE_AT_ONSET)) > 0,
                        na.rm = TRUE),
            '<=75'= sum(as.numeric(as.character(AGE_AT_ONSET)) <= 75  & as.numeric(as.character(AGE_AT_ONSET)) > 0,
                        na.rm = TRUE))

# df <- covars %>% 
#   mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999, 0), "unknown"))) %>%
#   group_by(STATUS, ETHNICITY) %>%
#   summarise('<=65' = sum(as.numeric(as.character(AGE_AT_ONSET)) <= 65 & as.numeric(as.character(AGE_AT_ONSET)) > 0 , 
#                         na.rm = TRUE), 
#             '<=70'= sum(as.numeric(as.character(AGE_AT_ONSET)) <= 70  & as.numeric(as.character(AGE_AT_ONSET)) > 0, 
#                        na.rm = TRUE),
#             '<=75'= sum(as.numeric(as.character(AGE_AT_ONSET)) <= 75  & as.numeric(as.character(AGE_AT_ONSET)) > 0, 
#                        na.rm = TRUE),
#             'MISSING'= sum(is.na(AGE_AT_ONSET)))

df <- as.data.frame(df)
df
## require(tidyverse)
# df <- df %>% pivot_wider(Ethnicity, names_from = STATUS, values_from = c(`<65`,`<70`,`<75`))
library(reshape2)
# summarize values
# df <- melt(df, id.vars = 1:2, measure.vars = 3:6)
df <- reshape2::melt(df, id.vars = 1:2, measure.vars = 3:5)
df <- reshape2::dcast(df, ETHNICITY ~ STATUS + variable)
# colnames(df) <- gsub(".*_", "", colnames(df))
df
###############################################################
################# age at last assessment ######################
###############################################################

library(dplyr)
df2 <- covars %>%
  mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999), NA))) %>%
  group_by(STATUS, ETHNICITY) %>%
  summarise('>=70'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 70,
                        na.rm = TRUE),
            '>=80'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 80,
                        na.rm = TRUE),
            '>=85'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 85,
                        na.rm = TRUE))


# df2 <- covars %>% 
#   mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999), NA))) %>%
#   group_by(STATUS, ETHNICITY) %>%
#   summarise('>=70'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 70, 
#                        na.rm = TRUE),
#             '>=80'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 80, 
#                        na.rm = TRUE),
#             '>=85'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 85, 
#                        na.rm = TRUE),
#             'MISSING'= sum(is.na(AGE_LAST_VISIT)))

df2 <- as.data.frame(df2)
## require(tidyverse)
# df <- df %>% pivot_wider(Ethnicity, names_from = STATUS, values_from = c(`<65`,`<70`,`<75`))
library(reshape2)
# summarize values
# df2 <- melt(df2, id.vars = 1:2, measure.vars = 3:6)
df2 <- reshape2::melt(df2, id.vars = 1:2, measure.vars = 3:5)
df2 <- reshape2::dcast(df2, ETHNICITY ~ STATUS + variable)
# colnames(df) <- gsub(".*_", "", colnames(df))
df2

Age_stratified <- cbind(df,df2)
Age_stratified
write.csv(Age_stratified, "stratified_classification.csv", quote = FALSE)


# N Control Missing age
MISSING <- covars %>% 
  mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999), NA))) %>%
  group_by(STATUS, ETHNICITY) %>%
  summarise('MISSING'= sum(is.na(AGE_LAST_VISIT)))

MISSING <- reshape2::melt(MISSING, id.vars = 1:2, measure.vars = 3)
MISSING <- reshape2::dcast(MISSING, ETHNICITY ~ STATUS + variable)
MISSING

# N Cases Missing age
MISSING <- covars %>% 
  mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999), NA))) %>%
  group_by(STATUS, ETHNICITY) %>%
  summarise('MISSING'= sum(is.na(AGE_AT_ONSET)))

MISSING <- reshape2::melt(MISSING, id.vars = 1:2, measure.vars = 3)
MISSING <- reshape2::dcast(MISSING, ETHNICITY ~ STATUS + variable)
MISSING

table(APOE4=covars$APOE4ANY, STATUS=covars$STATUS)

## Status: control = 1, case = 2, MCI = 3
## Sex: 1 = male, 2 = female, -9 = unknown
df <- table(covars$SEX, covars$ETHNICITY)[2:3,]
df <- rbind(df,colSums(df))
df
# % Percent of females:
(df[2,]/df[3,])*100

## % APOE4
# df <- table(covars$APOE4ANY, covars$ETHNICITY)[2:3,]
# df <- rbind(df,colSums(df))
# df
# (df[2,]/df[3,])*100
(table( covars[, c("APOE4ANY", "ETHNICITY") ] )[3,]/ as.vector(aggregate(covars$IID, list(covars$ETHNICITY), length)[2])) * 100

# Number of cases Controls by Ethnicity
table(covars$STATUS, covars$ETHNICITY)

# % APOE4+ for each ethnicity
# table(covars$APOE4ANY, covars$STATUS, covars$ETHNICITY)
# table(STATUS=NHW$status==1, APOE=NHW$apoe4any==1 )

# % CONTROLS APOE4+
(table( covars[ covars$STATUS == "CO" , c("APOE4ANY", "ETHNICITY") ] )[3,]/table( covars[ covars$STATUS == "CO" , c("ETHNICITY") ] )) * 100

# % CASES APOE4+
(table( covars[ covars$STATUS == "CA" , c("APOE4ANY", "ETHNICITY") ] )[3,]/table( covars[ covars$STATUS == "CA" , c("ETHNICITY") ] )) * 100

# % MCI APOE4+
(table( covars[ covars$STATUS == "MCI" , c("APOE4ANY", "ETHNICITY") ] )[2,]/table( covars[ covars$STATUS == "MCI" , c("ETHNICITY") ] )) * 100


## APOE4
cbind.data.frame(table( covars[ , c("STATUS","APOE4ANY", "ETHNICITY") ] ))



######################################################################
##################### GWAS Phenotype data ############################
######################################################################
## Case=2; control=1


change_names <- function(x){
  colnames(x) [grepl("AGE_AT_LAST_VISIT|AGE_AT_LAST_ASSESSMENT|AGE_AT_EXAM|AGE_AT_LAST", colnames(x), ignore.case = TRUE)] <- "AGE_LAST_VISIT"
  colnames(x) [grepl("AGE_AT_DEATH", colnames(x), ignore.case = TRUE)] <- "AGE_AT_DEATH"
  colnames(x) [grepl("ONSET|AAO", colnames(x), ignore.case = TRUE)] <- "AGE_AT_ONSET"
  colnames(x) [grepl("GENDER|SEX", colnames(x), ignore.case = TRUE)] <- "SEX"
  x}


## ADNI Core
ADNI_Core <- read.table("/home/achal/ADGC_GWAS/EOAD_GWAS/ADNI_Core_phenotype_April2021.csv", sep = ",", header = T)
# View(ADNI_Core)
ADNI_Core$STATUS_AN <- as.character(ADNI_Core$Final_CC_status)
as.data.frame(table(ADNI_Core$STATUS_AN))
# > as.data.frame(table(ADNI_Core$STATUS_AN))
# Var1 Freq
# 1        -9  234
# 2 Clinic_CA 2055 CA
# 3 Clinic_CO 1283 CO
# 4 Clinic_OT  127 -9
# 5  Neuro_AD   46 CA
# 6  Neuro_OT   18 -9
# 7    OT(CO)   79 CO


sum(grepl("Clinic_CA|Clinic_CO|Neuro_AD|OT\\(CO\\)", ADNI_Core$STATUS_AN))
ADNI_Core <- ADNI_Core[grepl("Clinic_CA|Clinic_CO|Neuro_AD|OT\\(CO\\)", ADNI_Core$STATUS_AN),]
ADNI_Core$STATUS_AN[grepl("Clinic_CO|OT\\(CO\\)", ADNI_Core$STATUS_AN)] <- 1
ADNI_Core$STATUS_AN[grepl("Clinic_CA|Neuro_AD", ADNI_Core$STATUS_AN)] <- 2

ADNI_Core <- change_names(ADNI_Core)

###############
## PD_seq_835
PD_seq <- read.table("/home/achal/ADGC_GWAS/EOAD_GWAS/PD_seq_835samples_pheno_20170919.csv", sep = ",", header = T)
# View(PD_seq)
as.data.frame(table(PD_seq$case_control_status))
# Case  Control Excluded 
# 562      272        1 
# Don't use cases from PD_seq
PD_seq <- PD_seq[grepl("Control", PD_seq$case_control_status),]
PD_seq$STATUS_AN <- as.character(PD_seq$case_control_status)
table(PD_seq$STATUS_AN)
PD_seq$STATUS_AN[grepl("Control", PD_seq$STATUS_AN)] <- 1
# PD_seq$STATUS_AN[grepl("Case", PD_seq$STATUS_AN)] <- 2
table(PD_seq$STATUS_AN)
PD_seq <- change_names(PD_seq)

###############
## NIALOAD
NIALOAD <- read.table("/home/achal/ADGC_GWAS/EOAD_GWAS/NIALOAD_uniq_core_phenotype_April2021.csv", sep = ",", header = T)
# View(NIALOAD)
as.data.frame(table(NIALOAD$Final_CC_status))
# > as.data.frame(table(NIALOAD$Final_CC_status))
# Var1 Freq
# 1                  9498
# 2               AD 2942 CA
# 3           AD_MCI  353 CA
# 4        AD_Unsure  410 -9
# 5     ADAD_carrier    2 X
# 6      ADAD_family   18 X
# 7  ADAD_noncarrier    4 CO
# 8         C9ORF72+    2 X
# 9               CA  149 CA 
# 10              CO 2420 CO
# 11             DLB   13 X
# 12              DS   14 X
# 13             FTD    7 X
# 14  LRRK2_positive   15 X
# 15             MCI   25 ?
# 16        Neuro_AD   24 CA
# 17    Neuro_AD_DLB    1 CA 
# 18        Neuro_CO    1 CO
# 19       Neuro_DLB    1 -9
# 20        Neuro_OT    1 -9
# 21 Neuro_Presyntom    1 CA
# 22              OT 2528 -9
# 23          OT(CO)  575 CO
# 24              PD   19 X
# 25     PGRN_family    1 -9
NIALOAD$STATUS_AN <- as.character(NIALOAD$Final_CC_status)
sum(grepl("^AD$|AD_MCI|ADAD_noncarrier|^CA$|^CO$|Neuro_AD|Neuro_AD_DLB|Neuro_CO|Neuro_Presyntom|OT\\(CO\\)", NIALOAD$STATUS_AN))
NIALOAD <- NIALOAD[grepl("^AD$|AD_MCI|ADAD_noncarrier|^CA$|^CO$|Neuro_AD|Neuro_AD_DLB|Neuro_CO|Neuro_Presyntom|OT\\(CO\\)", NIALOAD$STATUS_AN),]
table(NIALOAD$STATUS_AN)

NIALOAD$STATUS_AN[grepl("ADAD_noncarrier|^CO$|Neuro_CO|OT\\(CO\\)", NIALOAD$STATUS_AN)] <- 1
NIALOAD$STATUS_AN[grepl("^AD$|AD_MCI|^CA$|Neuro_AD|Neuro_AD_DLB|Neuro_Presyntom", NIALOAD$STATUS_AN)] <- 2
table(NIALOAD$STATUS_AN)
NIALOAD <- change_names(NIALOAD)

###############
## NACC
NACC <- read.table("/home/achal/ADGC_GWAS/EOAD_GWAS/NACC_uniq_core_phenotype_April2021.csv", sep = ",", header = T)
# View(NACC)
as.data.frame(table(NACC$final_CC_status))
# > as.data.frame(table(NACC$final_CC_status))
# Var1  Freq
# 1                              7
# 2             ADAD_carrier    25
# 3              ADAD_family     6
# 4          ADAD_noncarrier     1
# 5                 C9ORF72+     4
# 6                       CA 11060
# 7                       CO  9397
# 8                      DLB   830
# 9                       DS     3
# 10                     FTD  1414
# 11            MAPT_carrier     2
# 12             MAPT_family     7
# 13                Neuro_AD  1805
# 14            Neuro_AD_DLB     3
# 15            Neuro_AD_FTD     4
# 16                Neuro_CO   210
# 17               Neuro_DLB     2
# 18                Neuro_OT   935
# 19 Neuro_PreSymptomatic_AD   104
# 20                      OT 10405
# 21                  OT(CO)  2593
# 22                      PD     9
# 23             PGRN_family    10


NACC$STATUS_AN <- as.character(NACC$final_CC_status)
sum(grepl("ADAD_carrier|ADAD_family|ADAD_noncarrier|^CA$|^CO$|Neuro_AD|Neuro_AD_DLB|Neuro_CO|Neuro_PreSymptomatic_AD|OT\\(CO\\)", NACC$STATUS_AN))
NACC <- NACC[grepl("ADAD_carrier|ADAD_family|ADAD_noncarrier|^CA$|^CO$|Neuro_AD|Neuro_AD_DLB|Neuro_CO|Neuro_PreSymptomatic_AD|OT\\(CO\\)", NACC$STATUS_AN),]
table(NACC$STATUS_AN)

NACC$STATUS_AN[grepl("ADAD_noncarrier|^CO$|Neuro_CO|OT\\(CO\\)", NACC$STATUS_AN)] <- 1
NACC$STATUS_AN[grepl("ADAD_carrier|ADAD_family|^CA$|Neuro_AD|Neuro_AD_DLB|Neuro_PreSymptomatic_AD", NACC$STATUS_AN)] <- 2
table(NACC$STATUS_AN)
NACC <- change_names(NACC)

###############
## MAP_Core
MAP_Core <- read.table("/home/achal/ADGC_GWAS/EOAD_GWAS/MAP_Core_pheno_April2021.csv", sep = ",", header = T)
# View(MAP_Core)
as.data.frame(table(MAP_Core$Final_CC_status))
# > as.data.frame(table(MAP_Core$Final_CC_status))
# Var1 Freq
# 1                   374
# 2     ADAD_carrier    7
# 3      ADAD_family   64
# 4  ADAD_noncarrier    4
# 5         C9ORF72+   15
# 6               CA 1753
# 7          Clin_AD    4
# 8         Clin_DLB    1
# 9         Clin_FTD    1
# 10 Clin_Presyntoma    1
# 11              CO 1594
# 12             DLB   99
# 13             FTD   79
# 14 LRRRK2_positive   31
# 15    MAPT_carrier    3
# 16     MAPT_family   29
# 17        Neuro_AD  311
# 18    Neuro_AD_DLB   16
# 19    Neuro_AD_FTD    6
# 20     Neuro_AD_PD    1
# 21       Neuro_ALS    1
# 22        Neuro_CO   27
# 23       Neuro_DLB   12
# 24       Neuro_FTD    4
# 25        Neuro_OT   33
# 26     Neuro_Parry    3
# 27        Neuro_PD    2
# 28 Neuro_Presyntom   17
# 29       Neuro_PSP    1
# 30              OT  664
# 31          OT(CO)  159
# 32              PD   29
# 33    PGRN_carrier    6
# 34     PGRN_family  208
# 35 PGRN_noncarrier    3



MAP_Core$STATUS_AN <- as.character(MAP_Core$Final_CC_status)
sum(grepl("ADAD_carrier|ADAD_family|ADAD_noncarrier|^CA$|^CO$|Neuro_AD|Neuro_AD_DLB|Neuro_CO|Neuro_PreSymptomatic_AD|OT\\(CO\\)", MAP_Core$STATUS_AN))
MAP_Core <- MAP_Core[grepl("ADAD_carrier|ADAD_family|ADAD_noncarrier|^CA$|^CO$|Neuro_AD|Neuro_AD_DLB|Neuro_CO|Neuro_PreSymptomatic_AD|OT\\(CO\\)", MAP_Core$STATUS_AN),]
table(MAP_Core$STATUS_AN)

MAP_Core$STATUS_AN[grepl("ADAD_noncarrier|^CO$|Neuro_CO|OT\\(CO\\)", MAP_Core$STATUS_AN)] <- 1
MAP_Core$STATUS_AN[grepl("ADAD_carrier|ADAD_family|^CA$|Neuro_AD|Neuro_AD_DLB|Neuro_PreSymptomatic_AD", MAP_Core$STATUS_AN)] <- 2
table(MAP_Core$STATUS_AN)
MAP_Core <- change_names(MAP_Core)

###############
## PPMI_Screening
PPMI_Screening <- read.table("/home/achal/ADGC_GWAS/EOAD_GWAS/PPMI_Screening_Demographics.csv", sep = ",", header = T)
# View(PPMI_Screening)
as.data.frame(table(PPMI_Screening$F_STATUS))
# > as.data.frame(table(PPMI_Screening$F_STATUS))
# Var1 Freq
# 1    S   21
# 2    V 2209

# No age at onset ?

###############
## PPMI_Current
PPMI_Current <- read.table("/home/achal/ADGC_GWAS/EOAD_GWAS/PPMI_Current_Biospecimen_Analysis_Results.csv", sep = ",", header = T)
# View(PPMI_Current)
as.data.frame(table(PPMI_Current$DIAGNOSIS))
# > as.data.frame(table(PPMI_Current$DIAGNOSIS))
# Var1   Freq
# 1          Control  70806 CO
# 2   Genetic Cohort  37063 CO
# 3 Genetic Registry  19393
# 4               PD 153674
# 5        Prodromal   5796
# 6            SWEDD  10696

# No age at onset


########################################################################################################
####################################### Merge all of them ##############################################
########################################################################################################
ADNI_Core <- ADNI_Core[c("rid", "AGE_AT_ONSET",  "AGE_LAST_VISIT", "SEX", "STATUS_AN", "race", "APOE", "flag_APOE4", "Ethnicity", "education", "cdr_last_assessment")]
ADNI_Core <- ADNI_Core[1:7]
colnames(ADNI_Core)[1] <- "ID"
ADNI_Core$COHORT <- "ADNI_Core"

PD_seq <- PD_seq[c("PD_ID", "AGE_AT_ONSET", "AGE_LAST_VISIT", "SEX", "STATUS_AN", "zadock_Joint_VCF_ID", "case_control_status", "tau", "ptau", "ab", "major_project", "sub_project")]
PD_seq <- PD_seq[1:5]
colnames(PD_seq)[1] <- "ID"
PD_seq$COHORT <- "PD_seq"
PD_seq$race <- "Unknown"
PD_seq$APOE <- NA

NIALOAD <- NIALOAD[c("NIALOAD_ID","AGE_AT_ONSET", "AGE_LAST_VISIT", "SEX", "STATUS_AN", "race", "APOE", "AGE_AT_DEATH", "Site_ID", "family_ID", "Father_id", "Mother_id", "Final_CC_status", "YOB", "YOD", "Ethnicity", "Education", "cdr_last", "New_NCRAD_ID", "map_id")]
NIALOAD <- NIALOAD[1:7]
colnames(NIALOAD)[1] <- "ID"
NIALOAD$COHORT <- "NIALOAD"

NACC <- NACC[c("NACC_ID", "AGE_AT_ONSET", "AGE_LAST_VISIT", "SEX", "STATUS_AN", "race", "APOE", "AGE_AT_DEATH", "final_CC_status", "Year_of_Birth", "Year_of_death", "Ethnicity", "Education", "cdr_last", "MAP_ID")]
NACC <- NACC[1:7]
colnames(NACC)[1] <- "ID"
NACC$COHORT <- "NACC"

MAP_Core <- MAP_Core[c("MAP_ID", "AGE_AT_ONSET", "AGE_LAST_VISIT", "SEX", "STATUS_AN", "race", "APOE", "Final_CC_status", "DOB", "flag_APOE4", "Ethnicity", "NIALOAD_ID", "education", "cdr_last_assessment")]
MAP_Core <- MAP_Core[1:7]
colnames(MAP_Core)[1] <- "ID"
MAP_Core$COHORT <- "MAP_Core"


covars <- rbind.data.frame(ADNI_Core, PD_seq, NIALOAD, NACC, MAP_Core)
colnames(covars)[colnames(covars)=="race"] <- "ETHNICITY"
covars$ETHNICITY <- as.character(covars$ETHNICITY)
table(covars$ETHNICITY)
covars$ETHNICITY[grepl("American Indian", covars$ETHNICITY, ignore.case = T)] <- "AMERICAN_INDIAN"
covars$ETHNICITY[grepl("Asian", covars$ETHNICITY, ignore.case = T)] <- "ASIAN"
covars$ETHNICITY[grepl("BLACK", covars$ETHNICITY, ignore.case = T)] <- "BLACK"
covars$ETHNICITY[grepl("Missing|Unknown|Other|More than", covars$ETHNICITY, ignore.case = T)] <- "UNKNOWN"
covars$ETHNICITY[covars$ETHNICITY ==""] <- "UNKNOWN"
covars$ETHNICITY[covars$ETHNICITY =="White"] <- "NHW"

covars$SEX <- as.character(covars$SEX)
covars$SEX[covars$SEX == "Male"] <- 1
covars$SEX[covars$SEX == "Female"] <- 2
covars$SEX[covars$SEX == "Unknown"] <- -9
table(covars$SEX)
colnames(covars)[colnames(covars)=="STATUS_AN"] <- "STATUS"


## reformat for CACO on Spreadsheet
spreadsheets <- cbind.data.frame(FID = covars$ID, IID = covars$ID, Cohort = covars$ETHNICITY, Sex = covars$SEX, CACO = covars$STATUS, Status= covars$STATUS, COHORT = covars$COHORT)
ADNI_Core <- spreadsheets[spreadsheets$COHORT == "ADNI_Core",]
PD_seq <- spreadsheets[spreadsheets$COHORT == "PD_seq",]
NIALOAD <- spreadsheets[spreadsheets$COHORT == "NIALOAD",]
NACC <- spreadsheets[spreadsheets$COHORT == "NACC",]
MAP_Core <- spreadsheets[spreadsheets$COHORT == "MAP_Core",]


setwd("/home/achal/ADGC_GWAS/EOAD_GWAS/CA_CO_status/")
write.csv(ADNI_Core, "ADNI_Core.csv", quote = FALSE, row.names = FALSE)
write.csv(PD_seq, "PD_seq.csv", quote = FALSE, row.names = FALSE)
write.csv(NIALOAD, "NIALOAD.csv", quote = FALSE, row.names = FALSE)
write.csv(NACC, "NACC.csv", quote = FALSE, row.names = FALSE)
write.csv(MAP_Core, "MAP_Core.csv", quote = FALSE, row.names = FALSE)

########################################################################################################
####################### Age-of-onset stratified by STATUS, ethnicity and sex ###########################
########################################################################################################

# table(covars[covars$AGE_AT_ONSET <65 & covars$AGE_AT_ONSET > 0, c("STATUS", "SEX")])

library(dplyr)
df <- covars %>%
  mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999, 0), "unknown"))) %>%
  group_by(STATUS, ETHNICITY) %>%
  summarise('<=65' = sum(as.numeric(as.character(AGE_AT_ONSET)) <= 65 & as.numeric(as.character(AGE_AT_ONSET)) > 0 ,
                         na.rm = TRUE),
            '<=70'= sum(as.numeric(as.character(AGE_AT_ONSET)) <= 70  & as.numeric(as.character(AGE_AT_ONSET)) > 0,
                        na.rm = TRUE),
            '<=75'= sum(as.numeric(as.character(AGE_AT_ONSET)) <= 75  & as.numeric(as.character(AGE_AT_ONSET)) > 0,
                        na.rm = TRUE))

# df <- covars %>% 
#   mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999, 0), "unknown"))) %>%
#   group_by(STATUS, ETHNICITY) %>%
#   summarise('<=65' = sum(as.numeric(as.character(AGE_AT_ONSET)) <= 65 & as.numeric(as.character(AGE_AT_ONSET)) > 0 , 
#                         na.rm = TRUE), 
#             '<=70'= sum(as.numeric(as.character(AGE_AT_ONSET)) <= 70  & as.numeric(as.character(AGE_AT_ONSET)) > 0, 
#                        na.rm = TRUE),
#             '<=75'= sum(as.numeric(as.character(AGE_AT_ONSET)) <= 75  & as.numeric(as.character(AGE_AT_ONSET)) > 0, 
#                        na.rm = TRUE),
#             'MISSING'= sum(is.na(AGE_AT_ONSET)))

df <- as.data.frame(df)
df
## require(tidyverse)
# df <- df %>% pivot_wider(Ethnicity, names_from = STATUS, values_from = c(`<65`,`<70`,`<75`))
library(reshape2)
# summarize values
# df <- melt(df, id.vars = 1:2, measure.vars = 3:6)
df <- melt(df, id.vars = 1:2, measure.vars = 3:5)
df <- dcast(df, ETHNICITY ~ STATUS + variable)
# colnames(df) <- gsub(".*_", "", colnames(df))
df
###############################################################
################# age at last assessment ######################
###############################################################

library(dplyr)
df2 <- covars %>%
  mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999), NA))) %>%
  group_by(STATUS, ETHNICITY) %>%
  summarise('>=70'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 70,
                        na.rm = TRUE),
            '>=80'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 80,
                        na.rm = TRUE),
            '>=85'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 85,
                        na.rm = TRUE))


# df2 <- covars %>% 
#   mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999), NA))) %>%
#   group_by(STATUS, ETHNICITY) %>%
#   summarise('>=70'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 70, 
#                        na.rm = TRUE),
#             '>=80'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 80, 
#                        na.rm = TRUE),
#             '>=85'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 85, 
#                        na.rm = TRUE),
#             'MISSING'= sum(is.na(AGE_LAST_VISIT)))

df2 <- as.data.frame(df2)
## require(tidyverse)
# df <- df %>% pivot_wider(Ethnicity, names_from = STATUS, values_from = c(`<65`,`<70`,`<75`))
library(reshape2)
# summarize values
# df2 <- melt(df2, id.vars = 1:2, measure.vars = 3:6)
df2 <- reshape2::melt(df2, id.vars = 1:2, measure.vars = 3:5)
df2 <- reshape2::dcast(df2, ETHNICITY ~ STATUS + variable)
# colnames(df) <- gsub(".*_", "", colnames(df))
df2

Age_stratified <- cbind(df,df2)
Age_stratified
write.csv(Age_stratified, "stratified_classification.csv", quote = FALSE)


# N Control Missing age
MISSING <- covars %>% 
  mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999), NA))) %>%
  group_by(STATUS, ETHNICITY) %>%
  summarise('MISSING'= sum(is.na(AGE_LAST_VISIT)))

MISSING <- reshape2::melt(MISSING, id.vars = 1:2, measure.vars = 3)
MISSING <- reshape2::dcast(MISSING, ETHNICITY ~ STATUS + variable)
MISSING

# N Cases Missing age
MISSING <- covars %>% 
  mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999), NA))) %>%
  group_by(STATUS, ETHNICITY) %>%
  summarise('MISSING'= sum(is.na(AGE_AT_ONSET)))

MISSING <- reshape2::melt(MISSING, id.vars = 1:2, measure.vars = 3)
MISSING <- reshape2::dcast(MISSING, ETHNICITY ~ STATUS + variable)
MISSING


# Recode APOE4ANY
covars$APOE4ANY[is.na(covars$APOE)] <- -9


covars$APOE4ANY[grepl("22", covars$APOE)] <- 0
covars$APOE4ANY[grepl("23|32", covars$APOE)] <- 0
covars$APOE4ANY[grepl("33", covars$APOE)] <- 0
covars$APOE4ANY[grepl("24|42", covars$APOE)] <- 1
covars$APOE4ANY[grepl("34|43", covars$APOE)] <- 1
covars$APOE4ANY[grepl("44", covars$APOE)] <- 1


table(APOE4ANY=covars$APOE4ANY, STATUS=covars$STATUS)

## Status: control = 1, case = 2, MCI = 3
## Sex: 1 = male, 2 = female, -9 = unknown
df <- table(covars$SEX, covars$ETHNICITY)[2:3,]
df <- rbind(df,colSums(df))
df
# % Percent of females:
(df[2,]/df[3,])*100

## Total
table(covars$ETHNICITY)
 
## % APOE4
# df <- table(covars$APOE4ANY, covars$ETHNICITY)[2:3,]
# df <- rbind(df,colSums(df))
# df
# (df[2,]/df[3,])*100
(table( covars[, c("APOE4ANY", "ETHNICITY") ] )[3,]/ as.vector(aggregate(covars$ID, list(covars$ETHNICITY), length)[2])) * 100

# Number of cases Controls by Ethnicity
table(covars$STATUS, covars$ETHNICITY)

# % APOE4+ for each ethnicity
# table(covars$APOE4ANY, covars$STATUS, covars$ETHNICITY)
# table(STATUS=NHW$status==1, APOE=NHW$apoe4any==1 )

# % CONTROLS APOE4+
(table( covars[ covars$STATUS == "1" , c("APOE4ANY", "ETHNICITY") ] )[3,]/table( covars[ covars$STATUS == "1" , c("ETHNICITY") ] )) * 100

# % CASES APOE4+
(table( covars[ covars$STATUS == "2" , c("APOE4ANY", "ETHNICITY") ] )[3,]/table( covars[ covars$STATUS == "2" , c("ETHNICITY") ] )) * 100

# % MCI APOE4+
(table( covars[ covars$STATUS == "3" , c("APOE4ANY", "ETHNICITY") ] )[3,]/table( covars[ covars$STATUS == "3" , c("ETHNICITY") ] )) * 100
################################### 10/27/2020 Vicky only wants age of onset for Cases and MCI and Age at last assessment for CO
# CO > 75, 80 and 85

# ######################################################################################################
# ####################### Age-of-onset stratified by STATUS, Project and sex ###########################
# ######################################################################################################
# 
# # table(covars[covars$AGE_AT_ONSET <65 & covars$AGE_AT_ONSET > 0, c("STATUS", "SEX")])
# 
# library(dplyr)
# df <- covars %>%
#   mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999, 0), "unknown"))) %>%
#   group_by(STATUS, COHORT) %>%
#   summarise('<=65' = sum(as.numeric(as.character(AGE_AT_ONSET)) <= 65 & as.numeric(as.character(AGE_AT_ONSET)) > 0 ,
#                          na.rm = TRUE),
#             '<=70'= sum(as.numeric(as.character(AGE_AT_ONSET)) <= 70  & as.numeric(as.character(AGE_AT_ONSET)) > 0,
#                         na.rm = TRUE),
#             '<=75'= sum(as.numeric(as.character(AGE_AT_ONSET)) <= 75  & as.numeric(as.character(AGE_AT_ONSET)) > 0,
#                         na.rm = TRUE))
# 
# # df <- covars %>% 
# #   mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999, 0), "unknown"))) %>%
# #   group_by(STATUS, ETHNICITY) %>%
# #   summarise('<=65' = sum(as.numeric(as.character(AGE_AT_ONSET)) <= 65 & as.numeric(as.character(AGE_AT_ONSET)) > 0 , 
# #                         na.rm = TRUE), 
# #             '<=70'= sum(as.numeric(as.character(AGE_AT_ONSET)) <= 70  & as.numeric(as.character(AGE_AT_ONSET)) > 0, 
# #                        na.rm = TRUE),
# #             '<=75'= sum(as.numeric(as.character(AGE_AT_ONSET)) <= 75  & as.numeric(as.character(AGE_AT_ONSET)) > 0, 
# #                        na.rm = TRUE),
# #             'MISSING'= sum(is.na(AGE_AT_ONSET)))
# 
# df <- as.data.frame(df)
# df
# ## require(tidyverse)
# # df <- df %>% pivot_wider(Ethnicity, names_from = STATUS, values_from = c(`<65`,`<70`,`<75`))
# library(reshape2)
# # summarize values
# # df <- melt(df, id.vars = 1:2, measure.vars = 3:6)
# df <- melt(df, id.vars = 1:2, measure.vars = 3:5)
# df <- dcast(df, COHORT ~ STATUS + variable)
# # colnames(df) <- gsub(".*_", "", colnames(df))
# df
# ###############################################################
# ################# age at last assessment ######################
# ###############################################################
# 
# library(dplyr)
# df2 <- covars %>%
#   mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999), NA))) %>%
#   group_by(STATUS, COHORT) %>%
#   summarise('>=70'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 70,
#                         na.rm = TRUE),
#             '>=80'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 80,
#                         na.rm = TRUE),
#             '>=85'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 85,
#                         na.rm = TRUE))
# 
# 
# # df2 <- covars %>% 
# #   mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999), NA))) %>%
# #   group_by(STATUS, ETHNICITY) %>%
# #   summarise('>=70'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 70, 
# #                        na.rm = TRUE),
# #             '>=80'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 80, 
# #                        na.rm = TRUE),
# #             '>=85'= sum(as.numeric(as.character(AGE_LAST_VISIT)) >= 85, 
# #                        na.rm = TRUE),
# #             'MISSING'= sum(is.na(AGE_LAST_VISIT)))
# 
# df2 <- as.data.frame(df2)
# ## require(tidyverse)
# # df <- df %>% pivot_wider(Ethnicity, names_from = STATUS, values_from = c(`<65`,`<70`,`<75`))
# library(reshape2)
# # summarize values
# # df2 <- melt(df2, id.vars = 1:2, measure.vars = 3:6)
# df2 <- melt(df2, id.vars = 1:2, measure.vars = 3:5)
# df2 <- dcast(df2, COHORT ~ STATUS + variable)
# # colnames(df) <- gsub(".*_", "", colnames(df))
# df2
# 
# Age_stratified <- cbind(df,df2)
# Age_stratified
# write.csv(Age_stratified, "stratified_classification_GWAS.csv", quote = FALSE)
# 
# 
# # N Control Missing age
# MISSING <- covars %>% 
#   mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999), NA))) %>%
#   group_by(STATUS, ETHNICITY) %>%
#   summarise('MISSING'= sum(is.na(AGE_LAST_VISIT)))
# 
# MISSING <- reshape2::melt(MISSING, id.vars = 1:2, measure.vars = 3)
# MISSING <- reshape2::dcast(MISSING, ETHNICITY ~ STATUS + variable)
# MISSING
# 
# # N Cases Missing age
# MISSING <- covars %>% 
#   mutate(across(starts_with('AGE'), ~replace(., . %in% c(-9, 888, 999), NA))) %>%
#   group_by(STATUS, ETHNICITY) %>%
#   summarise('MISSING'= sum(is.na(AGE_AT_ONSET)))
# 
# MISSING <- reshape2::melt(MISSING, id.vars = 1:2, measure.vars = 3)
# MISSING <- reshape2::dcast(MISSING, ETHNICITY ~ STATUS + variable)
# MISSING
# 
# table(covars$APOE4ANY, covars$STATUS)
# 
# ## Status: control = 1, case = 2, MCI = 3
# ## Sex: 1 = male, 2 = female, -9 = unknown
# df <- table(covars$SEX, covars$COHORT)[2:3,]
# df <- rbind(df,colSums(df))
# df
# # % Percent of females:
# (df[2,]/df[3,])*100
# 
# ## % APOE4
# # df <- table(covars$APOE4ANY, covars$ETHNICITY)[2:3,]
# # df <- rbind(df,colSums(df))
# # df
# # (df[2,]/df[3,])*100
# (table( covars[, c("APOE4ANY", "ETHNICITY") ] )[3,]/ as.vector(aggregate(covars$IID, list(covars$ETHNICITY), length)[2])) * 100
# 
# # Number of cases Controls by Ethnicity
# table(covars$STATUS, covars$COHORT)
# 
# # % APOE4+ for each ethnicity
# # table(covars$APOE4ANY, covars$STATUS, covars$ETHNICITY)
# # table(STATUS=NHW$status==1, APOE=NHW$apoe4any==1 )
# 
# # % CONTROLS APOE4+
# (table( covars[ covars$STATUS == "CO" , c("APOE4ANY", "ETHNICITY") ] )[3,]/table( covars[ covars$STATUS == "CO" , c("ETHNICITY") ] )) * 100
# 
# # % CASES APOE4+
# (table( covars[ covars$STATUS == "CA" , c("APOE4ANY", "ETHNICITY") ] )[3,]/table( covars[ covars$STATUS == "CA" , c("ETHNICITY") ] )) * 100
# 
# # % MCI APOE4+
# (table( covars[ covars$STATUS == "MCI" , c("APOE4ANY", "ETHNICITY") ] )[3,]/table( covars[ covars$STATUS == "MCI" , c("ETHNICITY") ] )) * 100

