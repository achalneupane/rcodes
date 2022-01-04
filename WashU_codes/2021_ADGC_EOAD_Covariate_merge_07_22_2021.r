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

# Fix APOE4any column; AFRICAN cohorts do not have APOE4any, so I will add values in APOE4any based on APOE1 and APOE2
covars$APOE1_2 <- paste0(covars$APOE_1, covars$APOE_1)

covars$APOE4ANY [is.na(covars$APOE4ANY)] <- covars$APOE1_2[is.na(covars$APOE4ANY)]

HISPANIC <- covars[grepl ("HISPANIC", covars$ETHNICITY),]
AFRICAN <- covars[grepl ("AFRICAN", covars$ETHNICITY),]
NHW <- covars[grepl ("NHW", covars$ETHNICITY),]
ASIAN <- covars[grepl ("ASIAN", covars$ETHNICITY),]

write.table(HISPANIC, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/HISPANIC_covariate.tsv", sep ="\t", col.names = T, quote = F, row.names = FALSE)
write.table(AFRICAN, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/African_American_covariate.tsv", sep ="\t", col.names = T, quote = F, row.names = FALSE)
write.table(NHW, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/NHW_covariate.tsv", sep ="\t", col.names = T, quote = F, row.names = FALSE)
write.table(ASIAN, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ASIAN_covariate.tsv", sep ="\t", col.names = T, quote = F, row.names = FALSE)


######################################################################
######################ADGC Demographics by Ethnicity##################
######################################################################
# Replacing any weird numbers with NAs
FIX_weird_NUMS <- function(x){
  x <- as.numeric(as.character(x))
  x [x %in% c(-9, 888, 999, 0)] <- NA
  return(x)
}

covars[grepl("AGE",colnames(covars))] <- sapply(covars[grepl("AGE",colnames(covars))], FIX_weird_NUMS)

# # Remove duplicates
# library(tidyverse)
# df <- covars %>% group_by_at(vars(IID)) %>% filter(n()>1) %>% ungroup()
## Note ~ 5000 IID seem to have duplicates, but are not duplicates. If you check the STATUS, SEX and mergedID they do not match. 

Get_STATs <- function(covars){
  
TOTAL = nrow(covars)
# CA,CO, MCI
N.controls <- sum(covars$STATUS=="CO")
N.cases <- sum(covars$STATUS=="CA")
N.mci <- sum(covars$STATUS=="MCI")

# Number of CA <= 65 and 70 (Cases==2)
CA.65 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == "CA","AGE_AT_ONSET"])))) <= 65)
CA.70 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == "CA","AGE_AT_ONSET"])))) <= 70)
# Number of CO > 70 (Cases==2)
CO.70 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == "CO","AGE_LAST_VISIT"])))) > 70)
# Number of MCI <= 65 and 70 (Cases==2)
MCI.65 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == "MCI","AGE_AT_ONSET"])))) <= 65)
MCI.70 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == "MCI","AGE_AT_ONSET"])))) <= 70)
# Number of Others
N.OTHERS <- sum(covars$STATUS == "unknown")

# Percent Female
PERC.FEMALE <- (sum(covars$SEX == 2, na.rm = T)/(sum(covars$SEX == 1, na.rm = T)+ sum(covars$SEX == 2, na.rm = T))) *100

# MISSING AGES CA and CO
N.CA.missing.age <- sum(is.na(covars [covars$STATUS == "CA", "AGE_AT_ONSET"]))
N.CO.missing.age <- sum(is.na(covars [covars$STATUS == "CO", "AGE_LAST_VISIT"]))


# Percent APOE4
# PERC.APOE <- (table( covars[, c("APOE4ANY") ] )[3]/(table( covars[, c("APOE4ANY") ] )[2] + table( covars[, c("APOE4ANY") ] )[3]))*100 
POS <- sum(covars[, c("APOE4ANY")] == 1, na.rm = T)
NEG <- sum(covars[, c("APOE4ANY")] == 0, na.rm = T)
UNK <- sum(covars[, c("APOE4ANY")] == -9, na.rm = T)
PERC.APOE <- (POS/ (POS + NEG)) *100

CA.APOE.PERC <- sum(covars[ covars$STATUS == "CA" , c("APOE4ANY")] == 1, na.rm = T)/ (N.cases) *100
CO.APOE.PERC <- sum(covars[ covars$STATUS == "CO" , c("APOE4ANY")] == 1, na.rm = T)/ (N.controls) *100
MCI.APOE.PERC <- sum(covars[ covars$STATUS == "MCI" , c("APOE4ANY")] == 1, na.rm = T)/ (N.mci) *100

# For Ethnicity
STATS <- cbind(ETHNICITY = unique(covars$ETHNICITY), TOTAL = TOTAL, '% FEMALE' = round(PERC.FEMALE, 2), '% APOE' = round(PERC.APOE, 2), 'N CONTROLS (1)' = N.controls, 'N CASES (2)' = N.cases, 'N MCI (3)' = N.mci, 'N CONTROLS > 70 yo' = CO.70,
               'N CONTROLS missing age' = N.CO.missing.age, 'N CASES ≤ 65 yo' = CA.65, 'N CASES ≤ 70 yo' = CA.70, 'N CASES missing age' = N.CA.missing.age, 'N MCI (3) ≤ 65 yo' = MCI.65, 'N MCI (3) ≤ 70 yo' = MCI.70, 'N OTHERS (-9)' = N.OTHERS,
               '% CONTROLS (1) APOE4+' = round(CO.APOE.PERC, 2), '% CASES (2) APOE4+' = round(CA.APOE.PERC, 2), '% MCI (3) APOE4+' = round(MCI.APOE.PERC, 2))

# # ## For STUDY
# STATS <- cbind(ETHNICITY = unique(covars$ETHNICITY), STUDY = unique(covars$COHORT), TOTAL = TOTAL, '% FEMALE' = round(PERC.FEMALE, 2), '% APOE' = round(PERC.APOE, 2), 'N CONTROLS (1)' = N.controls, 'N CASES (2)' = N.cases, 'N MCI (3)' = N.mci, 'N CONTROLS > 70 yo' = CO.70,
#                'N CONTROLS missing age' = N.CO.missing.age, 'N CASES ≤ 65 yo' = CA.65, 'N CASES ≤ 70 yo' = CA.70, 'N CASES missing age' = N.CA.missing.age, 'N MCI (3) ≤ 65 yo' = MCI.65, 'N MCI (3) ≤ 70 yo' = MCI.70, 'N OTHERS (-9)' = N.OTHERS,
#                '% CONTROLS (1) APOE4+' = round(CO.APOE.PERC, 2), '% CASES (2) APOE4+' = round(CA.APOE.PERC, 2), '% MCI (3) APOE4+' = round(MCI.APOE.PERC, 2))


return(STATS)
}

# By Ethnicity
ADGC_GWAS_by_ETHNICITY <- t(do.call(rbind, by(covars, covars$ETHNICITY, FUN = Get_STATs)))

write.table(ADGC_GWAS_by_ETHNICITY, "ADGC_GWAS_by_ETHNICITY.txt", sep ="\t", col.names = T, quote = F, row.names = T)
# BY study
ADGC_GWAS_by_STUDY <- do.call(rbind, by(covars, covars$COHORT, FUN = Get_STATs))

write.table(ADGC_GWAS_by_STUDY, "ADGC_GWAS_by_STUDY.txt", sep ="\t", col.names = T, quote = F, row.names = F)


######################################################################
##################### GWAS Phenotype data ############################
######################################################################
## Case=2; control=1


change_names <- function(x){
  colnames(x) [grepl("AGE_AT_LAST_VISIT|AGE_AT_LAST_ASSESSMENT|AGE_AT_EXAM|AGE_AT_LAST|age_last", colnames(x), ignore.case = TRUE)] <- "AGE_LAST_VISIT"
  colnames(x) [grepl("AGE_AT_DEATH", colnames(x), ignore.case = TRUE)] <- "AGE_AT_DEATH"
  colnames(x) [grepl("ONSET|AAO", colnames(x), ignore.case = TRUE)] <- "AGE_AT_ONSET"
  colnames(x) [grepl("GENDER|SEX|sex", colnames(x), ignore.case = TRUE)] <- "SEX"
  colnames(x)[grepl("final_CC_status", colnames(x), ignore.case = TRUE)] <- "Final_CC_status"
  colnames(x)[grepl("Ethnicity", colnames(x), ignore.case = TRUE)] <- "ETHNICITY"
  x}


RECODE_CACO_SEX <- function(covars){
  covars$STATUS <- as.character(covars$STATUS)
  covars$STATUS[covars$STATUS=="CA"] <- 2
  covars$STATUS[covars$STATUS=="CO"] <- 1
  covars$STATUS[covars$STATUS=="MCI"] <- 3
  covars$STATUS[covars$STATUS=="unknown"] <- -9
  covars$SEX[covars$SEX == "Female"] <- 2
  covars$SEX[covars$SEX == "Male"] <- 1
  return(covars)
}



# Note : After our ADGC meeting on 11/10/2021, Carlos came and asked me to remove any ADAD samples

# Here, I am reading these phenotype files
#################
### ADNI Core ###
#################
ADNI_Core <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/WASHU-GWAS/ADNI_unique_core_pheno_20211102.csv", sep = ",", header = T)
colnames(ADNI_Core)[colnames(ADNI_Core) == "final_CC_status"] <- "Final_CC_status"
ADNI_Core$STATUS_AN <- as.character(ADNI_Core$Final_CC_status)
as.data.frame(table(ADNI_Core$STATUS_AN))
# > as.data.frame(table(ADNI_Core$STATUS_AN))
# Var1 Freq
# 1            236
# 2 Clinic_CA 2122
# 3 Clinic_CO 1318
# 4 Clinic_OT  125
# 5  Neuro_AD   61
# 6  Neuro_OT   19
# 7    OT(CO)   95

# Selecting only the AD samples
sum(grepl("^Clinic_CA$|^Clinic_CO$|^Neuro_AD$|OT\\(CO\\)", ADNI_Core$STATUS_AN))
# 3596
ADNI_Core <- ADNI_Core[grepl("^Clinic_CA$|^Clinic_CO$|^Neuro_AD$|OT\\(CO\\)", ADNI_Core$STATUS_AN),]
ADNI_Core$STATUS_AN[grepl("^Clinic_CO$|OT\\(CO\\)", ADNI_Core$STATUS_AN)] <- 1
ADNI_Core$STATUS_AN[grepl("^Clinic_CA$|^Neuro_AD$", ADNI_Core$STATUS_AN)] <- 2

ADNI_Core <- change_names(ADNI_Core)

###############
### NIALOAD ###
###############
NIALOAD <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/WASHU-GWAS/NIALOAD_uniq_core_phenotype_April2021.csv", sep = ",", header = T)
as.data.frame(table(NIALOAD$Final_CC_status))
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

# Selecting only the AD samples
NIALOAD$STATUS_AN <- as.character(NIALOAD$Final_CC_status)
sum(grepl("^AD$|AD_MCI|^CA$|^CO$|^Neuro_AD$|^Neuro_AD_DLB$|^Neuro_CO$|^Neuro_Presyntom$|OT\\(CO\\)", NIALOAD$STATUS_AN))
# 6466
NIALOAD <- NIALOAD[grepl("^AD$|AD_MCI|^CA$|^CO$|^Neuro_AD$|^Neuro_AD_DLB$|^Neuro_CO$|^Neuro_Presyntom$|OT\\(CO\\)", NIALOAD$STATUS_AN),]
table(NIALOAD$STATUS_AN)

NIALOAD$STATUS_AN[grepl("^CO$|^Neuro_CO$|OT\\(CO\\)", NIALOAD$STATUS_AN)] <- 1
NIALOAD$STATUS_AN[grepl("^AD$|^AD_MCI$|^CA$|^Neuro_AD$|^Neuro_AD_DLB$|^Neuro_Presyntom$", NIALOAD$STATUS_AN)] <- 2
table(NIALOAD$STATUS_AN)
NIALOAD <- change_names(NIALOAD)


############
### NACC ###
############
# NACC <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/WASHU-GWAS/NACC_uniq_core_phenotype_April2021.csv", sep = ",", header = T)
NACC <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/WASHU-GWAS/NACC_unique_core_pheno_20210927.csv", sep = ",", header = T)
NACC_ALL <- NACC
dim(NACC_ALL)
# 43999
# View(NACC)
as.data.frame(table(NACC$final_CC_status))
# Var1  Freq
# 1             ADAD_carrier    25
# 2              ADAD_family     6
# 3          ADAD_noncarrier     1
# 4                 C9ORF72+     4
# 5                       CA 12449
# 6                       CO 10739
# 7                      DLB   978
# 8                       DS     3
# 9                      FTD  1641
# 10            MAPT_carrier     2
# 11             MAPT_family     3
# 12                Neuro_AD  1785
# 13            Neuro_AD_DLB     3
# 14            Neuro_AD_FTD     4
# 15                Neuro_CO   210
# 16               Neuro_DLB     2
# 17                Neuro_OT   957
# 18 Neuro_PreSymptomatic_AD   139
# 19                      OT 12035
# 20                  OT(CO)  2994
# 21                      PD     9
# 22             PGRN_family    10


NACC$STATUS_AN <- as.character(NACC$final_CC_status)
sum(grepl("^CA$|^CO$|^Neuro_AD$|^Neuro_AD_DLB$|^Neuro_CO$|^Neuro_PreSymptomatic_AD$|OT\\(CO\\)", NACC$STATUS_AN))
# 28319
# Selecting only the AD samples
NACC <- NACC[grepl("^CA$|^CO$|^Neuro_AD$|^Neuro_AD_DLB$|^Neuro_CO$|^Neuro_PreSymptomatic_AD$|OT\\(CO\\)", NACC$STATUS_AN),]
table(NACC$STATUS_AN)

NACC$STATUS_AN[grepl("^CO$|^Neuro_CO$|OT\\(CO\\)", NACC$STATUS_AN)] <- 1
NACC$STATUS_AN[grepl("^CA$|^Neuro_AD$|^Neuro_AD_DLB$|^Neuro_PreSymptomatic_AD$", NACC$STATUS_AN)] <- 2
table(NACC$STATUS_AN)
NACC <- change_names(NACC)
colnames(NACC)[colnames(NACC) == "final_CC_status"] <- "Final_CC_status"

################
### MAP_Core ###
################
MAP_Core <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/WASHU-GWAS/MAP_Core_pheno_April2021.csv", sep = ",", header = T)
# Per Vicky on 12/17/2021
# so NHW = White + not Hispanic or LAtino
# latino = White + Hispanic or LAtino

as.data.frame(table(MAP_Core$Final_CC_status))
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
sum(grepl("^CA$|^CO$|^Neuro_AD$|^Neuro_AD_DLB$|^Neuro_CO$|^Neuro_PreSymptomatic_AD$|OT\\(CO\\)", MAP_Core$STATUS_AN))
# 3860
# Selecting only the AD samples
MAP_Core <- MAP_Core[grepl("^CA$|^CO$|^Neuro_AD$|^Neuro_AD_DLB$|^Neuro_CO$|^Neuro_PreSymptomatic_AD$|OT\\(CO\\)", MAP_Core$STATUS_AN),]
table(MAP_Core$STATUS_AN)

MAP_Core$STATUS_AN[grepl("^CO$|^Neuro_CO$|OT\\(CO\\)", MAP_Core$STATUS_AN)] <- 1
MAP_Core$STATUS_AN[grepl("Clin_AD|^CA$|^Neuro_AD$|^Neuro_AD_DLB$|^Neuro_PreSymptomatic_AD$", MAP_Core$STATUS_AN)] <- 2
table(MAP_Core$STATUS_AN)
MAP_Core <- change_names(MAP_Core)

##################
### PD_seq_835 ###
##################
PD_seq <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/WASHU-GWAS/PD_seq_835samples_pheno_20170919.csv", sep = ",", header = T)
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
colnames(PD_seq)[colnames(PD_seq) == "case_control_status"] <- "Final_CC_status"

######################
### PPMI_Screening ###
######################
PPMI_Screening <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/WASHU-GWAS/PPMI_Screening_Demographics.csv", sep = ",", header = T)
as.data.frame(table(PPMI_Screening$F_STATUS))
# > as.data.frame(table(PPMI_Screening$F_STATUS))
# Var1 Freq
# 1    S   21
# 2    V 2209

# No age at onset


####################
### PPMI_Current ###
####################
PPMI_Current <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/WASHU-GWAS/PPMI_Current_Biospecimen_Analysis_Results.csv", sep = ",", header = T)
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
PPMI_Current$Final_CC_status <- PPMI_Current$DIAGNOSIS

########################################################################################################
####################################### Merge all of them ##############################################
########################################################################################################
ADNI_Core <- ADNI_Core[c("rid", "AGE_AT_ONSET",  "AGE_LAST_VISIT", "SEX", "STATUS_AN", "race", "APOE",  "Final_CC_status", "flag_APOE4", "Ethnic", "Education", "cdr_last")]
ADNI_Core <- ADNI_Core[1:8]
colnames(ADNI_Core)[1] <- "ID"
ADNI_Core$COHORT <- "ADNI_Core"

PD_seq <- PD_seq[c("PD_ID", "AGE_AT_ONSET", "AGE_LAST_VISIT", "SEX", "STATUS_AN", "Final_CC_status", "zadock_Joint_VCF_ID", "tau", "ptau", "ab", "major_project", "sub_project")]
PD_seq <- PD_seq[1:6]
colnames(PD_seq)[1] <- "ID"
PD_seq$COHORT <- "PD_seq"
PD_seq$race <- "Unknown"
PD_seq$APOE <- NA

NIALOAD <- NIALOAD[c("NIALOAD_ID","AGE_AT_ONSET", "AGE_LAST_VISIT", "SEX", "STATUS_AN", "race", "APOE", "Final_CC_status", "AGE_AT_DEATH", "Site_ID", "family_ID", "Father_id", "Mother_id", "YOB", "YOD", "ETHNICITY", "Education", "cdr_last", "New_NCRAD_ID", "map_id")]
NIALOAD <- NIALOAD[1:8]
colnames(NIALOAD)[1] <- "ID"
NIALOAD$COHORT <- "NIALOAD"

NACC <- NACC[c("NACC_ID", "AGE_AT_ONSET", "AGE_LAST_VISIT", "SEX", "STATUS_AN", "race", "APOE", "Final_CC_status", "AGE_AT_DEATH", "Year_of_Birth", "Year_of_death", "ETHNICITY", "Education", "cdr_last", "MAP_ID")]
NACC <- NACC[1:8]
colnames(NACC)[1] <- "ID"
NACC$COHORT <- "NACC"

MAP_Core <- MAP_Core[c("MAP_ID", "AGE_AT_ONSET", "AGE_LAST_VISIT", "SEX", "STATUS_AN", "race", "APOE", "Final_CC_status", "DOB", "flag_APOE4", "ETHNICITY", "NIALOAD_ID", "education", "cdr_last_assessment")]
MAP_Core <- MAP_Core[1:8]
colnames(MAP_Core)[1] <- "ID"
MAP_Core$COHORT <- "MAP_Core"


covars_WASHU <- rbind.data.frame(ADNI_Core, PD_seq, NIALOAD, NACC, MAP_Core)
colnames(covars_WASHU)[colnames(covars_WASHU)=="race"] <- "ETHNICITY"
covars_WASHU$ETHNICITY <- as.character(covars_WASHU$ETHNICITY)
table(covars_WASHU$ETHNICITY)

# ## WashU Cases (<= 65 yo)
# CA.65 <- covars_WASHU[covars_WASHU$STATUS_AN == 2 & covars_WASHU$AGE_AT_ONSET <=65,]
# sum(is.na(CA.65$ID))
# CA.65 <- CA.65[!is.na(CA.65$ID),]
# # Writing this as per Vicky's suggestion on 10/14/2021
# write.csv(CA.65, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/WASHU-GWAS/WASHU_All_cases_65_or_younger.csv", quote = FALSE, row.names = FALSE)

covars_WASHU$ETHNICITY[grepl("American Indian", covars_WASHU$ETHNICITY, ignore.case = T)] <- "AMERICAN_INDIAN"
covars_WASHU$ETHNICITY[grepl("Asian", covars_WASHU$ETHNICITY, ignore.case = T)] <- "ASIAN"
covars_WASHU$ETHNICITY[grepl("BLACK", covars_WASHU$ETHNICITY, ignore.case = T)] <- "BLACK"
covars_WASHU$ETHNICITY[grepl("Missing|Unknown|Other|More than", covars_WASHU$ETHNICITY, ignore.case = T)] <- "UNKNOWN"
covars_WASHU$ETHNICITY[covars_WASHU$ETHNICITY ==""] <- "UNKNOWN"
covars_WASHU$ETHNICITY[covars_WASHU$ETHNICITY =="White"] <- "NHW"

covars_WASHU$SEX <- as.character(covars_WASHU$SEX)
covars_WASHU$SEX[covars_WASHU$SEX == "Male"] <- 1
covars_WASHU$SEX[covars_WASHU$SEX == "Female"] <- 2
covars_WASHU$SEX[covars_WASHU$SEX == "Unknown"] <- -9
table(covars_WASHU$SEX)
colnames(covars_WASHU)[colnames(covars_WASHU)=="STATUS_AN"] <- "STATUS"


# ## reformat for CACO on Spreadsheet
# spreadsheets <- cbind.data.frame(FID = covars_WASHU$ID, IID = covars_WASHU$ID, Cohort = covars_WASHU$ETHNICITY, Sex = covars_WASHU$SEX, CACO = covars_WASHU$STATUS, Status= covars_WASHU$STATUS, COHORT = covars_WASHU$COHORT)
# ADNI_Core <- spreadsheets[spreadsheets$COHORT == "ADNI_Core",]
# PD_seq <- spreadsheets[spreadsheets$COHORT == "PD_seq",]
# NIALOAD <- spreadsheets[spreadsheets$COHORT == "NIALOAD",]
# NACC <- spreadsheets[spreadsheets$COHORT == "NACC",]
# MAP_Core <- spreadsheets[spreadsheets$COHORT == "MAP_Core",]
# 
# 
# setwd("/home/achal/ADGC_GWAS/EOAD_GWAS/CA_CO_status/")
# write.csv(ADNI_Core, "ADNI_Core.csv", quote = FALSE, row.names = FALSE)
# write.csv(PD_seq, "PD_seq.csv", quote = FALSE, row.names = FALSE)
# write.csv(NIALOAD, "NIALOAD.csv", quote = FALSE, row.names = FALSE)
# write.csv(NACC, "NACC.csv", quote = FALSE, row.names = FALSE)
# write.csv(MAP_Core, "MAP_Core.csv", quote = FALSE, row.names = FALSE)


######################################################################
################WASHU_GWAS Demographics by Ethnicity##################
######################################################################

# Recode APOE4ANY
covars_WASHU$APOE4ANY[is.na(covars_WASHU$APOE)] <- -9


covars_WASHU$APOE4ANY[grepl("22", covars_WASHU$APOE)] <- 0
covars_WASHU$APOE4ANY[grepl("23|32", covars_WASHU$APOE)] <- 0
covars_WASHU$APOE4ANY[grepl("33", covars_WASHU$APOE)] <- 0
covars_WASHU$APOE4ANY[grepl("24|42", covars_WASHU$APOE)] <- 1
covars_WASHU$APOE4ANY[grepl("34|43", covars_WASHU$APOE)] <- 1
covars_WASHU$APOE4ANY[grepl("44", covars_WASHU$APOE)] <- 1


## recode STATUS
covars_WASHU$STATUS[grepl("1", covars_WASHU$STATUS)] <- "CO"
covars_WASHU$STATUS[grepl("2", covars_WASHU$STATUS)] <- "CA"
covars_WASHU$STATUS[grepl("3", covars_WASHU$STATUS)] <- "MCI"
covars_WASHU$STATUS[grepl("-9", covars_WASHU$STATUS)] <- "unknown"

colnames(covars_WASHU)[grepl ("ID",colnames(covars_WASHU))] <- "IID"

## Save WASHU_Covars
all_covars_WASHU <- covars_WASHU

write.table(covars_WASHU, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/WASHU-GWAS/ALL_WASHU_COVARIATES_SELECTED_BY_AD_STATUS.txt", sep ="\t", col.names = T, quote = F, row.names = T)

## Check WASHU duplicates
library(tidyverse)
df <- covars_WASHU %>% group_by_at(vars(IID)) %>% filter(n()>1) %>% ungroup()
## Note: there are not really duplicates, only duplicate IDs


# ## Remove any overlapping samples
# covars_WASHU_overlapping <- covars_WASHU[covars_WASHU$IID %in% covars$IID,]
# ## ADNI_Core Is not duplicate; removing 12915 NACC samples
# covars_WASHU_overlapping <- covars_WASHU_overlapping[!grepl ("ADNI", covars_WASHU_overlapping$COHORT),]
# covars_WASHU <- covars_WASHU[!covars_WASHU$IID %in% covars_WASHU_overlapping$IID,]


# Get_STATs(covars_WASHU)


Get_STATs <- function(covars){
  
  TOTAL = nrow(covars)
  # CA,CO, MCI
  N.controls <- sum(covars$STATUS=="CO")
  N.cases <- sum(covars$STATUS=="CA")
  N.mci <- sum(covars$STATUS=="MCI")
  
  # Number of CA <= 65 and 70 (Cases==2)
  CA.65 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == "CA","AGE_AT_ONSET"])))) <= 65)
  CA.70 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == "CA","AGE_AT_ONSET"])))) <= 70)
  # Number of CO > 70 (Cases==2)
  CO.70 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == "CO","AGE_LAST_VISIT"])))) > 70)
  CO.80 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == "CO","AGE_LAST_VISIT"])))) > 80)
  # Number of MCI <= 65 and 70 (Cases==2)
  MCI.65 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == "MCI","AGE_AT_ONSET"])))) <= 65)
  MCI.70 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == "MCI","AGE_AT_ONSET"])))) <= 70)
  # Number of Others
  N.OTHERS <- sum(covars$STATUS == "unknown")
  
  # Percent Female
  PERC.FEMALE <- (sum(covars$SEX == 2, na.rm = T)/(sum(covars$SEX == 1, na.rm = T)+ sum(covars$SEX == 2, na.rm = T))) *100
  
  # MISSING AGES CA and CO
  N.CA.missing.age <- sum(is.na(covars [covars$STATUS == "CA", "AGE_AT_ONSET"]))
  N.CO.missing.age <- sum(is.na(covars [covars$STATUS == "CO", "AGE_LAST_VISIT"]))
  
  
  # Percent APOE4
  # PERC.APOE <- (table( covars[, c("APOE4ANY") ] )[3]/(table( covars[, c("APOE4ANY") ] )[2] + table( covars[, c("APOE4ANY") ] )[3]))*100 
  POS <- sum(covars[, c("APOE4ANY")] == 1, na.rm = T)
  NEG <- sum(covars[, c("APOE4ANY")] == 0, na.rm = T)
  UNK <- sum(covars[, c("APOE4ANY")] == -9, na.rm = T)
  PERC.APOE <- (POS/ (POS + NEG)) *100
  
  
  
  CA.APOE.PERC <- (sum(covars[ covars$STATUS == "CA" , c("APOE4ANY")] == 1, na.rm = T)/ (N.cases)) *100
  CO.APOE.PERC <- (sum(covars[ covars$STATUS == "CO" , c("APOE4ANY")] == 1, na.rm = T)/ (N.controls)) *100
  MCI.APOE.PERC <- (sum(covars[ covars$STATUS == "MCI" , c("APOE4ANY")] == 1, na.rm = T)/ (N.mci)) *100
  
  STATS <- cbind(ETHNICITY = unique(covars$ETHNICITY), TOTAL = TOTAL, '% FEMALE' = round(PERC.FEMALE, 2), '% APOE' = round(PERC.APOE, 2), 'N CONTROLS (1)' = N.controls, 'N CASES (2)' = N.cases, 'N MCI (3)' = N.mci, 'N CONTROLS > 70 yo' = CO.70, 'N CONTROLS > 80 yo' = CO.80,
                 'N CONTROLS missing age' = N.CO.missing.age, 'N CASES ≤ 65 yo' = CA.65, 'N CASES ≤ 70 yo' = CA.70, 'N CASES missing age' = N.CA.missing.age, 'N MCI (3) ≤ 65 yo' = MCI.65, 'N MCI (3) ≤ 70 yo' = MCI.70, 'N OTHERS (-9)' = N.OTHERS, 
                 '% CONTROLS (1) APOE4+' = round(CO.APOE.PERC, 2), '% CASES (2) APOE4+' = round(CA.APOE.PERC, 2), '% MCI (3) APOE4+' = round(MCI.APOE.PERC, 2))
  
  # ## For STUDY
  # STATS <- cbind(STUDY = unique(covars$COHORT), TOTAL = TOTAL, '% FEMALE' = round(PERC.FEMALE, 2), '% APOE' = round(PERC.APOE, 2), 'N CONTROLS (1)' = N.controls, 'N CASES (2)' = N.cases, 'N MCI (3)' = N.mci, 'N CONTROLS > 70 yo' = CO.70,
  #                'N CONTROLS missing age' = N.CO.missing.age, 'N CASES ≤ 65 yo' = CA.65, 'N CASES ≤ 70 yo' = CA.70, 'N CASES missing age' = N.CA.missing.age, 'N MCI (3) ≤ 65 yo' = MCI.65, 'N MCI (3) ≤ 70 yo' = MCI.70, 'N OTHERS (-9)' = N.OTHERS, 
  #                '% CONTROLS (1) APOE4+' = round(CO.APOE.PERC, 2), '% CASES (2) APOE4+' = round(CA.APOE.PERC, 2), '% MCI (3) APOE4+' = round(MCI.APOE.PERC, 2))
  
  
  return(STATS)
}


# WASHU_GWAS_by_ETHNICITY <- do.call(rbind, by(covars_WASHU, covars_WASHU$ETHNICITY, FUN = Get_STATs))
# write.table(t(WASHU_GWAS_by_ETHNICITY), "WASHU_GWAS_by_ETHNICITY.txt", sep ="\t", col.names = T, quote = F, row.names = T)

# WASHU_GWAS_by_STUDY <- do.call(rbind, by(covars, covars$COHORT, FUN = Get_STATs))




#######################################################################################
# ## check if these are in our latest GWAs plink file
# # GWAS.fam <- read.table("/40/AD/GWAS_data/Combined_Plink/Full_Imputation/2021_Aug/All_cohorts_merged_sep_2021_clean.fam", stringsAsFactors = FALSE)
# GWAS.fam <- read.table("/40/AD/GWAS_data/Combined_Plink/Full_Imputation/2021_Aug/All_cohorts_merged_Nov_2021.fam", stringsAsFactors = FALSE)
# GWAS.fam$KEY1 <-   sapply(strsplit(GWAS.fam$V2,"\\^"), `[`, 1)

## Use ID matrix instead
ADGC_NHW.fam <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ADGC_NHW_Cohort.fam", stringsAsFactors = FALSE)
ADGC_Hispanic.fam <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Hispanic/ADGC_Hispanic_Cohort.fam", stringsAsFactors = FALSE)
ADGC_Asian.fam <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Asian/ADGC_Asian_Cohort.fam", stringsAsFactors = FALSE)
ADGC_AA.fam <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/ADGC_AA_Cohort.fam", stringsAsFactors = FALSE)

ADGC.fam <- rbind.data.frame(ADGC_NHW.fam, ADGC_Hispanic.fam, ADGC_AA.fam, ADGC_Asian.fam)


# # Cohorts
# table(CA.65$COHORT)

#############
## All MAP ##
#############

MAP_Core <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/WASHU-GWAS/MAP_Core_pheno_April2021.csv", sep = ",", header = T)
# Per Vicky on 12/17/2021
# so NHW = White + not Hispanic or LAtino
# latino = White + Hispanic or LAtino

#############
as.data.frame(table(MAP_Core$Final_CC_status))
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
sum(grepl("^CA$|^CO$|^Neuro_AD$|^Neuro_AD_DLB$|^Neuro_CO$|^Neuro_PreSymptomatic_AD$|OT\\(CO\\)", MAP_Core$STATUS_AN))
# 3860
# Selecting only the AD samples
MAP_Core <- MAP_Core[grepl("^CA$|^CO$|^Neuro_AD$|^Neuro_AD_DLB$|^Neuro_CO$|^Neuro_PreSymptomatic_AD$|OT\\(CO\\)", MAP_Core$STATUS_AN),]
table(MAP_Core$STATUS_AN)

MAP_Core$STATUS_AN[grepl("^CO$|^Neuro_CO$|OT\\(CO\\)", MAP_Core$STATUS_AN)] <- 1
MAP_Core$STATUS_AN[grepl("Clin_AD|^CA$|^Neuro_AD$|^Neuro_AD_DLB$|^Neuro_PreSymptomatic_AD$", MAP_Core$STATUS_AN)] <- 2
table(MAP_Core$STATUS_AN)
MAP_Core <- change_names(MAP_Core)

colnames(MAP_Core)[colnames(MAP_Core)== "MAP_ID"] <- "IID"


MAP_Core$APOE4ANY[grepl("22|23|33|32", MAP_Core$APOE)] <- 0
MAP_Core$APOE4ANY[grepl("24|34|44|42|43", MAP_Core$APOE)] <- 1


# Per Vicky on 12/17/2021
# so NHW = White + not Hispanic or LAtino
# latino = White + Hispanic or LAtino
as.data.frame(table(MAP_Core$ETHNICITY))
as.data.frame(table(MAP_Core$race))


#### ############
## Samira and Carlos suggested to read from IDMatrix instead of FAM file
GWAS.fam <- read.delim("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/WashU_samples/ID_matrix_hg38_Nov2021.csv", sep = ",")
dim(GWAS.fam)


# Extract all MAP samples from the latest GWAS file I got from Samira
## How many of WashU GWAS.fam have MAP
GWAS.fam.MAP <- GWAS.fam[!is.na(GWAS.fam$MAP),]
nrow(GWAS.fam.MAP)
# 6870
# There are 6174 MAP samples in WashU GWAS.fam
GWAS.fam.MAP$KEY1 <- paste0("MAP_", gsub("MAP_", "", GWAS.fam.MAP$MAP))
length(GWAS.fam.MAP$KEY1 [!duplicated(GWAS.fam.MAP$KEY1 )])
# 4870 unique MAP samples


# MAP_Core <- covars_WASHU[grepl("MAP_",covars_WASHU$COHORT),]

MAP_Core$KEY1 <- paste0("MAP_", gsub("MAP_", "", MAP_Core$IID))
# # of samples in phenotype 
nrow(MAP_Core)
# 3860

# MAP_Core is the phenotype file I received from Fengxian
# 3446 of 3860 WASHU samples are available in WashU's GWAS fam file
sum(MAP_Core$KEY1 %in% GWAS.fam.MAP$KEY1)
# 3446 

# Selecting only those that are in GWAS FAM file
MAP_Core <- MAP_Core[MAP_Core$KEY1 %in% GWAS.fam.MAP$KEY1,]
MAP_Core$Final_CC_status <- as.character(MAP_Core$Final_CC_status) 
as.data.frame(table(MAP_Core$Final_CC_status))
# Var1 Freq
# 1           CA 1473
# 2           CO 1466
# 3     Neuro_AD  310
# 4 Neuro_AD_DLB   15
# 5     Neuro_CO   27
# 6       OT(CO)  155


MAP_Core$RACE_NEW <- as.character(MAP_Core$race)
as.data.frame(table(as.character(MAP_Core$RACE_NEW)))
# Var1 Freq
# 1                                               4
# 2                           AFRICAN_AMERICAN  501
# 3             American Indian/Alaskan Native   10
# 4                                      Asian   11
# 5                         More than One Race    3
# 6 Native Hawaiian or Other Pacific islanders    2
# 7                                      White 2829
MAP_Core$RACE_NEW[MAP_Core$RACE_NEW == "Black/African American"] <- "African_American"
MAP_Core$RACE_NEW[MAP_Core$RACE_NEW == "American Indian/Alaskan Native"] <- "AI_AN"
MAP_Core$RACE_NEW[MAP_Core$RACE_NEW == "More than One Race"] <- "Mixed_race"
MAP_Core$RACE_NEW[MAP_Core$RACE_NEW == "Native Hawaiian or Other Pacific islanders"] <- "NH_PI"

MAP_Core$ETHNICITY_NEW <- as.character(MAP_Core$ETHNICITY)
as.data.frame(table(MAP_Core$ETHNICITY_NEW))
# Var1 Freq
# 1     Hispanic or Latino   26
# 2 Not Hispanic or Latino 1275
# 3                Unknown 2059

MAP_Core$ETHNICITY_NEW[MAP_Core$ETHNICITY_NEW == "Hispanic or Latino"] <- "Hispanic_or_Latino"
MAP_Core$ETHNICITY_NEW[MAP_Core$ETHNICITY_NEW == "Not Hispanic or Latino"] <- "Not_Hispanic_or_Latino"

MAP_Core$Race_Ethnicity <- paste(MAP_Core$RACE_NEW, MAP_Core$ETHNICITY_NEW, sep = ":")

as.data.frame(table(MAP_Core$Race_Ethnicity))
# Var1 Freq
# 1                      :Hispanic_or_Latino    1
# 2                  :Not_Hispanic_or_Latino    1
# 3                                 :Unknown    2
# 4      African_American:Hispanic_or_Latino    2
# 5  African_American:Not_Hispanic_or_Latino  193
# 6                 African_American:Unknown  306
# 7             AI_AN:Not_Hispanic_or_Latino    6
# 8                            AI_AN:Unknown    4
# 9             Asian:Not_Hispanic_or_Latino    7
# 10                           Asian:Unknown    4
# 11       Mixed_race:Not_Hispanic_or_Latino    2
# 12                      Mixed_race:Unknown    1
# 13                NH_PI:Hispanic_or_Latino    1
# 14                           NH_PI:Unknown    1
# 15                White:Hispanic_or_Latino   22
# 16            White:Not_Hispanic_or_Latino 1066
# 17                           White:Unknown 1741


## I will also run PCA to check whether White:Unknown are NHW. First I will
# extract these samples from GWAS data. 
sum(GWAS.fam.MAP$KEY1 %in% MAP_Core$KEY1)
# 5138
nrow(MAP_Core)
# 3446
GWAS_DATA <- GWAS.fam.MAP[GWAS.fam.MAP$KEY1 %in% MAP_Core$KEY1,]
GWAS_DATA$IID <- as.character(GWAS_DATA$IID)
GWAS_DATA$Array <- as.character(sapply(strsplit(GWAS_DATA$IID, "\\^"), `[`, 3))
as.data.frame(table(GWAS_DATA$Array))

dim(GWAS_DATA)
# [1] 5138   68

# GWAS_DATA_UNIQUE <- GWAS_DATA[!duplicated(GWAS_DATA$KEY1),]
# 
# write.table(GWAS_DATA_UNIQUE[1:2], "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/WashU_samples/WashU_AD_samples_for_PCA.txt", col.names = F, row.names = F, quote = F, sep = "\t")


write.table(GWAS_DATA[1:2], "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/WashU_samples/WashU_AD_samples.txt", col.names = F, row.names = F, quote = F, sep = "\t")



#########################################################################################



















# Since there are duplicate samples from different arrays, I will check the 
# frequency of their arrays. We want to minimize the number of different arrays
# used to increase the statistical power. I checked with Muhammad to decide
# which ones are the same arrays based on this file 
# (hg38_Array_Grouping_09082021.xlsx): 
# https://wustl.box.com/s/ks5rgedey1itbjhgjiuzp3muc7minuwm


sum(GWAS.fam$KEY1 %in% MAP_Core$KEY1)
# 4680
nrow(MAP_Core)
# 3360
GWAS_DATA <- GWAS.fam[GWAS.fam$KEY1 %in% MAP_Core$KEY1,]
GWAS_DATA$Array <- sapply(strsplit(GWAS_DATA$V2, "\\^"), `[`, 3)
as.data.frame(table(GWAS_DATA$Array))
# Var1 Freq
# 1                 2010_660W  717
# 2                2011OmniEx  226
# 3                2012OmniEx  155
# 4          2013_660k_NACCR1   31
# 5          2013_660k_NACCR2    5
# 6  2013Genentech_Omni1-Quad   49
# 7                2013OmniEx  138
# 8         2013OmniEx_NACCR3  122
# 9         2013OmniEx_NACCR4  116
# 10        2013OmniEx_NACCR5   58
# 11        2013OmniEx_NACCR6  114
# 12            2014BioBankUK   54
# 13               2014CoreEx  345
# 14               2015CoreEx  373
# 15      2015OmniExEx_NACCR7   94
# 16       2016_Human1M-Duov3  129
# 17      2016OmniExEx_NACCR8  147
# 18         201711_CoreExome   91
# 19              2017NeuroX2  564
# 20             201812_GSAv2  262
# 21         2018GSAv1_NACCR9  112
# 22    2019GSAMD24v1_NACCR10  222
# 23    202103_GSAv3_SNOWMASS  521
# 24                     660k   26
# 25                   OmniEx    9

GWAS_DATA$Array[grepl("^2013_660W$|^2010_660W$|^2013_660k_NACCR1$|^660k$|^2013_660k_NACCR2$", GWAS_DATA$Array)] <- "660"
GWAS_DATA$Array[grepl("^ADNIA4_GSA$|^201812_GSAv2$|^ADNI3_GSAv2$|^2019GSAMD24v1_NACCR10$|^2018GSAv1_NACCR9$|^202103_GSAv3_SNOWMASS$", GWAS_DATA$Array)] <- "GSA"
GWAS_DATA$Array[grepl("^ADNI2GO_OmniEx$|^ADNI2_OmniEx$|^2019_OmniExEx$|^2012OmniEx$|2016OmniExEx$|^ADNIDOD_OmniEx$|^2017OmniExEx$|^2011OmniEx$|^2013OmniEx_NACCR3$|^2013OmniEx_NACCR6$|^2013OmniEx_NACCR4$|^2016OmniExEx_NACCR8$|^2015OmniExEx_NACCR7$|^OmniEx$|^2014_Omni5Ex$|^2013OmniEx_NACCR5$|^2013OmniEx$", GWAS_DATA$Array)] <- "OmniEx"
GWAS_DATA$Array[grepl("^201711_CoreExome$|^2016_CoreEx$|^2014CoreEx$|^2015CoreEx$|^2014_CoreEx$|^2017_CoreEx$", GWAS_DATA$Array)] <- "CoreEx"

###############
## To Samira ##
###############
# I have asked Samira to test her Shiny tool for sample extraction based on array freq and missingness
IDMatrix <- read.delim("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/WashU_samples/ID_matrix_hg38_Nov2021.csv", sep = ",")
dim(IDMatrix)

MAP_Core$KEY1 %in% IDMatrix

MAP_Core$KEY <-  gsub("MAP_", "", MAP_Core$KEY1)
WashU_IDMATRIX <- IDMatrix[IDMatrix$MAP %in% MAP_Core$KEY,]
dim(WashU_IDMATRIX)
write.table(WashU_IDMATRIX[1:2], "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/WashU_samples/WashU_AD_samples.txt", col.names = F, row.names = F, quote = F, sep = "\t")
# --maf 0.01 --geno 0.02 --hwe 1e-6
# Prefix="WashU_MAP_for_ADGC"


###############

as.data.frame(table(GWAS_DATA$Array))
# Var1 Freq
# 1 2013Genentech_Omni1-Quad   49
# 2            2014BioBankUK   54
# 3       2016_Human1M-Duov3  129
# 4              2017NeuroX2  564
# 5                      660  779
# 6                   CoreEx  809
# 7                      GSA 1117
# 8                   OmniEx 1179

# How many times each array is repeated
GWAS_DATA$nArray <- with(transform(GWAS_DATA, n = 1),  ave(n, Array, FUN = length))

# Now select the unique samples
library(dplyr)
tt <- GWAS_DATA %>% group_by(KEY1) %>% filter(nArray==max(nArray))

tt[duplicated(tt$KEY1),]




MAP_Core_TABLE <- t(do.call(rbind, by(MAP_Core, MAP_Core$ETHNICITY, FUN = Get_STATs))) 


##########
## NACC ##
##########
# GWAS.NACC <- GWAS.fam[(grepl("NACC", GWAS.fam$V2)),]
# GWAS.NACC$KEY1 <- gsub("UW_|UPENN_", "", GWAS.NACC$KEY1)
# CA.65.NACC <- CA.65[grepl("NACC", CA.65$COHORT),]

# ## 23/ 3621
# sum(CA.65.NACC$ID %in% GWAS.NACC$KEY1)
# CA.65.NACC.NOTMATCHED <-  CA.65.NACC$ID[(!CA.65.NACC$ID %in% GWAS.NACC$KEY1)]
# 
# ## covars.NACC
# 
# sum(CA.65.NACC.NOTMATCHED %in% covars$IID)
# ADGC.fam.NACC <- ADGC.fam[grepl("NACC", ADGC.fam$V2),]
# ADGC.fam.NACC$KEY <- paste0("NACC", str_extract(ADGC.fam.NACC$V2, "(?<=NACC)[0-9]*"))
# 
# # Check how many unmatched now match in ADGC data
# # 2036/ 3621
# sum(CA.65.NACC.NOTMATCHED %in% ADGC.fam.NACC$KEY)

NACC_WASHU <- covars_WASHU[grepl("NACC", covars_WASHU$COHORT),]
dim(NACC_WASHU)
# 28319    10

# Check in ADGC fam if they are all available
ADGC_NACC.fam <- ADGC.fam [grepl("NACC", ADGC.fam$V2),]
dim(ADGC_NACC.fam)
# 22617     7

# clean NACC IDs as some have suffix/prefixes
ADGC_NACC.fam$KEY1 <- paste0("NACC", str_extract(ADGC_NACC.fam$V2, "(?<=NACC)[0-9]*"))

sum(NACC_WASHU$IID %in% ADGC_NACC.fam$KEY1)
# 14649

NACC_WASHU_FOUND_IN_ADGC.fam <- NACC_WASHU[NACC_WASHU$IID %in% ADGC_NACC.fam$KEY1,]

NACC_WASHU_NOTFOUND_IN_ADGC.fam <- NACC_WASHU$IID[!NACC_WASHU$IID %in% ADGC_NACC.fam$KEY1]
length(NACC_WASHU_NOTFOUND_IN_ADGC.fam)
# 13670

## check these in WashU GWAS data if NACC_WASHU_NOTFOUND_IN_ADGC.fam are there
GWAS.fam.NACC <- GWAS.fam[grepl("NACC", GWAS.fam$V2),]
GWAS.fam.NACC$KEY1 <- paste0("NACC", str_extract(GWAS.fam.NACC$V2, "(?<=NACC)[0-9]*"))
# Get those with FULL NACC IDs
sum(grepl("NACC[0-9]", GWAS.fam.NACC$KEY1))
# 171
GWAS.fam.NACC <- GWAS.fam.NACC[grepl("NACC[0-9]", GWAS.fam.NACC$KEY1),]
sum(NACC_WASHU_NOTFOUND_IN_ADGC.fam %in% GWAS.fam.NACC)
# 0 # No match!!
head(NACC_ALL)

# GET map ID of NACC_WASHU
NACC_ALL_with_MAP_ID <- NACC_ALL[!is.na(NACC_ALL$MAP_ID),]
dim(NACC_ALL_with_MAP_ID)
# 1814 # samples are MAP with NACC IDs
NACC_ALL_with_MAP_ID$MAP_ID <- paste0("MAP_", NACC_ALL_with_MAP_ID$MAP_ID)

NACC_WASHU$MAP_ID <- NACC_ALL_with_MAP_ID$MAP_ID[match(NACC_WASHU$IID, NACC_ALL_with_MAP_ID$NACC_ID)]

# How many of WashU PHENO have NACC MAP_IDs
sum(!is.na(NACC_WASHU$MAP_ID))
# 1527

# Check if NACC found in ADGC have MAP IDs
sum(NACC_WASHU_FOUND_IN_ADGC.fam$IID %in% NACC_ALL_with_MAP_ID$NACC_ID)
# 1184


# Check if NACC not found in ADGC have MAP IDs
sum(NACC_WASHU_NOTFOUND_IN_ADGC.fam %in% NACC_ALL_with_MAP_ID$NACC_ID)
# 344

# How many of WashU GWAS have MAP
GWAS.fam.MAP <- GWAS.fam[grepl("^MAP_", GWAS.fam$V2),]
dim(GWAS.fam.MAP)
# There are 6157 MAP samples in WashU GWAS.fam, but only 4626 of them are unique IDs
GWAS.fam.MAP$KEY1 <- paste0("MAP_", str_extract(GWAS.fam.MAP$V2, "(?<=MAP_)[0-9]*"))
length(GWAS.fam.MAP$KEY1 [!duplicated(GWAS.fam.MAP$KEY1 )])
# 4626
# Now check if any of GWAS MAP also have NACC IDs
sum(GWAS.fam.MAP$KEY1 [!duplicated(GWAS.fam.MAP$KEY1 )] %in% NACC_ALL_with_MAP_ID$MAP_ID)
# 1747

## There are only 171 NACC samples in latest GWAS data (from fam file Samira pointed out)
sum(GWAS.fam.NACC$KEY1 %in% NACC_ALL_with_MAP_ID$NACC_ID)
# 0

NACC_WASHU_FOUND_IN_ADGC.fam$Final_CC_status <- as.character(NACC_WASHU_FOUND_IN_ADGC.fam$Final_CC_status)

# Var1 Freq
# 1                      CA 6210
# 2                      CO 6395
# 3                Neuro_AD 1078
# 4            Neuro_AD_DLB    2
# 5                Neuro_CO  130
# 6 Neuro_PreSymptomatic_AD  107
# 7                  OT(CO)  727


NACC_Core_TABLE <- t(do.call(rbind, by(NACC_WASHU_FOUND_IN_ADGC.fam, NACC_WASHU_FOUND_IN_ADGC.fam$ETHNICITY, FUN = Get_STATs)))


# ## Check for DUber's CSF NACC
dim(ADGC_NACC.fam)
# [1] 22617     7
CSF_NACC <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/WASHU-GWAS/DUBER_cruchaga06302021csf_v2.csv", header =T, sep = ",")
head(CSF_NACC)
dim(CSF_NACC)
# 1715   24
CSF_NACC$NACCID <- as.character(CSF_NACC$NACCID)
CSF_NACC$ADGC_GWAS <- ifelse(CSF_NACC$NACCID %in% ADGC_NACC.fam$KEY1, "YES", "NO")
# Only 1220 are available in ADGC
# write.csv(CSF_NACC, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/WASHU-GWAS/DUBER_cruchaga06302021csf_v2_AN.csv", quote = FALSE, row.names = FALSE)

# Get the FID and IID for these
tt <- ADGC_NACC.fam[na.omit(match(CSF_NACC$NACCID, ADGC_NACC.fam$V2)),1:2]

tt <- tt[!duplicated(tt$V2),]
write.table(tt, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/WASHU-GWAS/extract_dubers_NACC.list", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# BFILE="ADGC_NHW_Cohort" # 768 people
# BFILE="ADGC_AA_Cohort" # 27 people
# BFILE="ADGC_Hispanic_Cohort" # 18 people
# BFILE="ADGC_Asian_Cohort" # 17 people
# plink1.9 --bfile ${BFILE} --keep /40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/WASHU-GWAS/extract_dubers_NACC.list --make-bed --keep-allele-order --out /40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/NACC_samples/CSF_NACC_4_DUBER/${BFILE}_CSF_NACC_DUBER
# plink1.9 --bfile ADGC_NHW_Cohort_CSF_NACC_DUBER --merge-list mylist.list --keep-allele-order --make-bed --out ADGC_CSF_NACC_DUBER







# #############
# ## NIALOAD ##
# #############
# GWAS.NIALOAD <- GWAS.fam[(grepl("NIALOAD", GWAS.fam$V2)),]
# GWAS.NIALOAD$KEY1 <- gsub("NIALOAD_", "", GWAS.NIALOAD$KEY1)
# 
# CA.65.NIALOAD <- CA.65[grepl("NIALOAD", CA.65$COHORT),]
# 
# ## 457/ 608
# sum(CA.65.NIALOAD$ID %in% GWAS.NIALOAD$KEY1)
# 
# 
# ##########
# ## ADNI ##
# ##########
# GWAS.ADNI <- GWAS.fam[(grepl("ADNI", GWAS.fam$V2)),]
# GWAS.ADNI$KEY1 <- gsub("UW_|UPENN_", "", GWAS.ADNI$KEY1)
# 
# CA.65.ADNI <- CA.65[grepl("ADNI", CA.65$COHORT),]
# CA.65.ADNI$CLEANED.ID <- paste0("ADNI_", CA.65.ADNI$ID)
# 
# ## 148/ 333
# sum(CA.65.ADNI$CLEANED.ID %in% GWAS.ADNI$KEY1)
# 
# 
# # Total matched samples
# 330 + 23 + 457 + 333
# # 1143


#####################################################################################################################
# Check how many of ADGC fam are available in covars
dim(covars)
dim(ADGC.fam)
sum(covars$IID %in% ADGC.fam$V2)
# 43874

# Get those that were not matched
NHW_NOT_FOUND_IN_COVAR <- ADGC_NHW.fam[!ADGC_NHW.fam$V2 %in% covars$IID,]
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
########################################## CLEANING OF PHENOTYPE FILE ###############################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

##########################
##########################
########## NHW ###########
##########################
##########################
# ADGC NHW covariate
ADGC_NHW_covar <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/NHW_covariate.tsv", header = T)
## ADGC FAM
# ADGC_NHW_FAM <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/ADGC_NHW_Cohort.fam", header = F)
## Identify FAM STUDY
# All FAM files
fileName <- "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_NHW/fam.list"

## Loop over a file connection
conn <- file(fileName,open="r")
linn <-readLines(conn)
df1 <- {}
for (i in 1:length(linn)){
  print(linn[i])
  COHORT <- str_match(linn[i], "ADGC_NHW-\\s*(.*?)\\s*_geno") [,2] 
  df.tmp <- read.table(linn[i], header = F, colClasses = 'character')  
  df.tmp$COHORT <- COHORT
  df1 <- rbind(df1,df.tmp)
  # Sys.sleep(2)
}

close(conn)


df1$V2 <- as.character(df1$V2)
# df1 <- df1[!duplicated(df1$V2),]
df1$KEY <- paste(df1$COHORT, df1$V2, sep = ":")

ADGC_NHW_FAM <- df1

# sum(ADGC_NHW_FAM$V2 %in% df1$V2)




head(ADGC_NHW_covar)
ADGC_NHW_covar$MERGED_ID <- as.character(ADGC_NHW_covar$MERGED_ID)


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
# NO   YES 
# 18230 32014 
# Now copy all MERGED_ID to CLEANED_ID to rows FOUND == "NO", and see if we have any matches
ADGC_NHW_covar$CLEANED_ID[ADGC_NHW_covar$FOUND == "NO"] <- as.character(ADGC_NHW_covar$MERGED_ID[ADGC_NHW_covar$FOUND == "NO"])

# CHECK MATCHES !!
ADGC_NHW_covar$FOUND <- ifelse(ADGC_NHW_covar$CLEANED_ID %in% ADGC_NHW_FAM$V2, "YES", "NO")
table(ADGC_NHW_covar$FOUND)
# NO   YES 
# 10335 39909 


# Now clean MAYO
library(strex)
ADGC_NHW_covar$CLEANED_ID[grepl("MAYO", ADGC_NHW_covar$COHORT)] <- paste0("0_", str_after_nth(ADGC_NHW_covar$CLEANED_ID[grepl("MAYO", ADGC_NHW_covar$COHORT)], "_", 3))

# CHECK MATCHES !!
ADGC_NHW_covar$FOUND <- ifelse(ADGC_NHW_covar$CLEANED_ID %in% ADGC_NHW_FAM$V2, "YES", "NO")
table(ADGC_NHW_covar$FOUND)
# NO   YES 
# 8332 41912 


# Next clean ACT_ACT_ IDs from MERGED_ID column of covar file and Replace in CLEANED_ID
ADGC_NHW_covar$CLEANED_ID[grepl("ACT_ACT+[[:digit:]]+_ACT+[[:digit:]]+$", ADGC_NHW_covar$MERGED_ID)] <- paste0("ACT", str_extract(ADGC_NHW_covar$MERGED_ID, "(?<=ACT_ACT)[0-9]*")) [grepl("ACT_ACT+[[:digit:]]+_ACT+[[:digit:]]+$", ADGC_NHW_covar$MERGED_ID)]

# CHECK MATCHES !!
ADGC_NHW_covar$FOUND <- ifelse(ADGC_NHW_covar$CLEANED_ID %in% ADGC_NHW_FAM$V2, "YES", "NO")
table(ADGC_NHW_covar$FOUND)
# NO   YES 
# 8173 42071 

# Now clean IDs for UMVUMSSM
ADGC_NHW_covar$CLEANED_ID[grepl("UMVUMSSM", ADGC_NHW_covar$COHORT)] <- gsub("_","-", ADGC_NHW_covar$SUM_ID[grepl("UMVUMSSM", ADGC_NHW_covar$COHORT)])

ADGC_NHW_covar$FOUND <- ifelse(ADGC_NHW_covar$CLEANED_ID %in% ADGC_NHW_FAM$V2, "YES", "NO")
table(ADGC_NHW_covar$FOUND)
# NO   YES 
# 6515 43729 


# Now Clean MIRAGE by extracting strings between x and y delimiters
XYDelim = function(x, n, i){
  do.call(c, lapply(x, function(X)
    paste(unlist(strsplit(X, "_"))[(n+1):(i)], collapse = "_")))
}

x <- as.character(ADGC_NHW_covar$MERGED_ID[grepl("MIRAGE", ADGC_NHW_covar$COHORT)])
ADGC_NHW_covar$CLEANED_ID[grepl("MIRAGE", ADGC_NHW_covar$COHORT)] <- paste(XYDelim(x = x, n = 1, i = 2), XYDelim(x = x, n = 1, i = 3), sep = "_")

# CHECK MATCHES !!
ADGC_NHW_covar$FOUND <- ifelse(ADGC_NHW_covar$CLEANED_ID %in% ADGC_NHW_FAM$V2, "YES", "NO")
table(ADGC_NHW_covar$FOUND)
# NO   YES 
# 5456 44788 


# Clean ADNI IDs by adding leading zero if length less than 7
ADGC_NHW_covar$CLEANED_ID [ADGC_NHW_covar$COHORT == "ADNI"] <- str_pad(as.character(ADGC_NHW_covar$IID [ADGC_NHW_covar$COHORT == "ADNI"]), 7, pad = "0")

# CHECK MATCHES !!
ADGC_NHW_covar$FOUND <- ifelse(ADGC_NHW_covar$CLEANED_ID %in% ADGC_NHW_FAM$V2, "YES", "NO")
table(ADGC_NHW_covar$FOUND)
# NO   YES 
# 4980 45264 

# Clean UPITT
ADGC_NHW_covar$CLEANED_ID [grepl("UPITT", ADGC_NHW_covar$COHORT)] <- gsub("UPITT_", "", ADGC_NHW_covar$CLEANED_ID [grepl("UPITT", ADGC_NHW_covar$COHORT)])

# CHECK MATCHES !!
ADGC_NHW_covar$FOUND <- ifelse(ADGC_NHW_covar$CLEANED_ID %in% ADGC_NHW_FAM$V2, "YES", "NO")
table(ADGC_NHW_covar$FOUND)
# NO   YES 
# 2761 47483 

#######################
### Clean UMVUTARC2 ###
#######################
as.character(ADGC_NHW_covar$MERGED_ID [grepl("UMVUTARC2", ADGC_NHW_covar$COHORT)]) 
# First TARC
as.character(ADGC_NHW_covar$MERGED_ID [grepl("UMVUTARC2", ADGC_NHW_covar$COHORT) & grepl("TARC", ADGC_NHW_covar$MERGED_ID)])
ADGC_NHW_covar$CLEANED_ID [grepl("UMVUTARC2", ADGC_NHW_covar$COHORT) & grepl("TARC", ADGC_NHW_covar$MERGED_ID)] <- paste0("TARC2_0_TARC2_0-", sapply(strsplit(as.character(ADGC_NHW_covar$MERGED_ID [grepl("UMVUTARC2", ADGC_NHW_covar$COHORT) & grepl("TARC", ADGC_NHW_covar$MERGED_ID)]), split="_"), "[", 3))
# CHECK MATCHES !!
ADGC_NHW_covar$FOUND <- ifelse(ADGC_NHW_covar$CLEANED_ID %in% ADGC_NHW_FAM$V2, "YES", "NO")
table(ADGC_NHW_covar$FOUND)
# NO   YES 
# 2463 47781 

# Then IHG
as.character(ADGC_NHW_covar$MERGED_ID [grepl("UMVUTARC2", ADGC_NHW_covar$COHORT) & grepl("IHG", ADGC_NHW_covar$MERGED_ID)])
ADGC_NHW_covar$CLEANED_ID [grepl("UMVUTARC2", ADGC_NHW_covar$COHORT) & grepl("IHG", ADGC_NHW_covar$MERGED_ID)] <- paste(sapply(strsplit(gsub("_","-", gsub("MTV_", "", as.character(ADGC_NHW_covar$MERGED_ID [grepl("UMVUTARC2", ADGC_NHW_covar$COHORT) & grepl("IHG", ADGC_NHW_covar$MERGED_ID)]))), split="-"), "[", 1), gsub("_","-", gsub("MTV_", "", as.character(ADGC_NHW_covar$MERGED_ID [grepl("UMVUTARC2", ADGC_NHW_covar$COHORT) & grepl("IHG", ADGC_NHW_covar$MERGED_ID)]))), sep = "_")

# CHECK MATCHES !!
ADGC_NHW_covar$FOUND <- ifelse(ADGC_NHW_covar$CLEANED_ID %in% ADGC_NHW_FAM$V2, "YES", "NO")
table(ADGC_NHW_covar$FOUND)
# NO   YES 
# 2340 47904 

# Then clean VAN
as.character(ADGC_NHW_covar$MERGED_ID [grepl("UMVUTARC2", ADGC_NHW_covar$COHORT) & grepl("VAN", ADGC_NHW_covar$MERGED_ID)])
ADGC_NHW_covar$CLEANED_ID [grepl("UMVUTARC2", ADGC_NHW_covar$COHORT) & grepl("VAN", ADGC_NHW_covar$MERGED_ID)] <- paste(sapply(strsplit(gsub("_","-", gsub("MTV_", "", as.character(ADGC_NHW_covar$MERGED_ID [grepl("UMVUTARC2", ADGC_NHW_covar$COHORT) & grepl("VAN", ADGC_NHW_covar$MERGED_ID)]))), split="-"), "[", 1), gsub("_","-", gsub("MTV_", "", as.character(ADGC_NHW_covar$MERGED_ID [grepl("UMVUTARC2", ADGC_NHW_covar$COHORT) & grepl("VAN", ADGC_NHW_covar$MERGED_ID)]))), sep = "_")

# CHECK MATCHES !!
ADGC_NHW_covar$FOUND <- ifelse(ADGC_NHW_covar$CLEANED_ID %in% ADGC_NHW_FAM$V2, "YES", "NO")
table(ADGC_NHW_covar$FOUND)
# NO   YES 
# 2266 47978 


# Then clean DUK
as.character(ADGC_NHW_covar$MERGED_ID [grepl("UMVUTARC2", ADGC_NHW_covar$COHORT) & grepl("DUK", ADGC_NHW_covar$MERGED_ID)])
ADGC_NHW_covar$CLEANED_ID [grepl("UMVUTARC2", ADGC_NHW_covar$COHORT) & grepl("DUK", ADGC_NHW_covar$MERGED_ID)] <- paste(sapply(strsplit(gsub("_","-", gsub("MTV_", "", as.character(ADGC_NHW_covar$MERGED_ID [grepl("UMVUTARC2", ADGC_NHW_covar$COHORT) & grepl("DUK", ADGC_NHW_covar$MERGED_ID)]))), split="-"), "[", 1), gsub("_","-", gsub("MTV_", "", as.character(ADGC_NHW_covar$MERGED_ID [grepl("UMVUTARC2", ADGC_NHW_covar$COHORT) & grepl("DUK", ADGC_NHW_covar$MERGED_ID)]))), sep = "_")

# CHECK MATCHES !!
ADGC_NHW_covar$FOUND <- ifelse(ADGC_NHW_covar$CLEANED_ID %in% ADGC_NHW_FAM$V2, "YES", "NO")
table(ADGC_NHW_covar$FOUND)
# NO   YES 
# 2263 47981 

# Clean WashU 
ADGC_NHW_covar$CLEANED_ID[grepl("WASHU1", ADGC_NHW_covar$COHORT)] <- gsub("_", "-", gsub("WASHU_","", ADGC_NHW_covar$MERGED_ID[grepl("WASHU1", ADGC_NHW_covar$COHORT)]))

ADGC_NHW_covar$FOUND <- ifelse(ADGC_NHW_covar$CLEANED_ID %in% ADGC_NHW_FAM$V2, "YES", "NO")
table(ADGC_NHW_covar$FOUND)
# NO   YES 
# 1593 48651 




# It looks like most of the samples that were not found are from ADC3 and GSKs
ADGC_NHW_covar$CLEANED_ID [ADGC_NHW_covar$COHORT == "ADC3" & ADGC_NHW_covar$FOUND == "NO"]
NOTFOUND <- ADGC_NHW_covar[ADGC_NHW_covar$FOUND == "NO",]

# Now clean GSK
ADGC_NHW_covar$CLEANED_ID [grepl("GSK", ADGC_NHW_covar$COHORT)] <- ADGC_NHW_covar$MERGED_ID [grepl("GSK", ADGC_NHW_covar$COHORT)]

ADGC_NHW_covar$FOUND <- ifelse(ADGC_NHW_covar$CLEANED_ID %in% ADGC_NHW_FAM$V2, "YES", "NO")
table(ADGC_NHW_covar$FOUND) # The number did not change from before because it was matching the partial IDs, so fixed that here
# NO   YES 
# 1593 48651 

NOTFOUND <- ADGC_NHW_covar[ADGC_NHW_covar$FOUND == "NO",]



# We still do not have pheno information for 1529 NHW samples (based on FAM file) and these are mostly NACC samples
nrow(ADGC_NHW_FAM) - sum(ADGC_NHW_FAM$V2 %in% ADGC_NHW_covar$CLEANED_ID)
ADGC_NHW_FAM$V2[!ADGC_NHW_FAM$V2 %in% ADGC_NHW_covar$CLEANED_ID]


# Let's check if any of these are present in FENGXIAN's Phenofile
FENGXIAN_COVAR <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/WASHU-GWAS/ALL_WASHU_COVARIATES_SELECTED_BY_AD_STATUS.txt", header = T, sep = "\t", stringsAsFactors = F)
FENGXIAN_COVAR_ALL <- FENGXIAN_COVAR
dim(FENGXIAN_COVAR)
head(FENGXIAN_COVAR)
head(ADGC_NHW_covar)
# available in Fengxian's pheno file
sum(ADGC_NHW_FAM$V2[!ADGC_NHW_FAM$V2 %in% ADGC_NHW_covar$CLEANED_ID] %in% FENGXIAN_COVAR$IID) # 935 of those are present in Fengxian's Phenofile
FENGXIAN_SAMPLES <- ADGC_NHW_FAM$V2[!ADGC_NHW_FAM$V2 %in% ADGC_NHW_covar$CLEANED_ID][ADGC_NHW_FAM$V2[!ADGC_NHW_FAM$V2 %in% ADGC_NHW_covar$CLEANED_ID] %in% FENGXIAN_COVAR$IID]

# Merge the samples Pheno from FENGXIAN
FENGXIAN_COVAR <- FENGXIAN_COVAR[FENGXIAN_COVAR$IID %in% FENGXIAN_SAMPLES,]
# Get STUDY info for Fengxian's Pheno
FENGXIAN_COVAR$COHORT <- ADGC_NHW_FAM$COHORT[match(FENGXIAN_COVAR$IID, ADGC_NHW_FAM$V2)]
FENGXIAN_COVAR$KEY <- paste(FENGXIAN_COVAR$COHORT, FENGXIAN_COVAR$IID, sep = ":")


# Missing samples
MISSING <- ADGC_NHW_FAM$V2[!ADGC_NHW_FAM$V2 %in% ADGC_NHW_covar$CLEANED_ID][!ADGC_NHW_FAM$V2[!ADGC_NHW_FAM$V2 %in% ADGC_NHW_covar$CLEANED_ID] %in% FENGXIAN_COVAR$IID]
length(MISSING)
# 594

# Next, there are duplicate sample IDs from two different STUDIES, so I will create a column with look up KEYs: 
# duplicates <- ADGC_NHW_FAM[duplicated(ADGC_NHW_FAM$V2),]
ADGC_NHW_covar$KEY <- paste(ADGC_NHW_covar$COHORT, ADGC_NHW_covar$CLEANED_ID, sep = ":")

# Drop Final_CC_status
FENGXIAN_COVAR <- FENGXIAN_COVAR[!grepl("Final_CC_status", colnames(FENGXIAN_COVAR))]

# Merge ADGC_NHW_covar and FENGXIAN_COVAR
head(ADGC_NHW_covar)
head(FENGXIAN_COVAR)
FENGXIAN_COVAR$FOUND <- "YES"



library(plyr)
combined <- rbind.fill(ADGC_NHW_covar, FENGXIAN_COVAR)






# 594 are missing
nrow(ADGC_NHW_FAM) - sum(ADGC_NHW_FAM$KEY %in% combined$KEY)
# [1] 594



# After running the analysis for the first time, I noticed that some of the samples do not have the right AAL or the latest AAL in ADGC Phenotype file mostly from NACC; AAO is also missing in some cases. I will get pheno information from Fengxian's phenotype file for those samples
dim(FENGXIAN_COVAR_ALL)

FENGXIAN_COVAR_NACC <- FENGXIAN_COVAR_ALL[FENGXIAN_COVAR_ALL$COHORT == "NACC",]
FENGXIAN_COVAR_NACC <- cbind.data.frame(FENGXIAN_COVAR_NACC, combined[match(FENGXIAN_COVAR_NACC$IID, combined$CLEANED_ID),])


FENGXIAN_COVAR_NACC <- FENGXIAN_COVAR_NACC[!is.na(FENGXIAN_COVAR_NACC$CLEANED_ID),]
dim(FENGXIAN_COVAR_NACC)
# 10969    58
FENGXIAN_COVAR_NACC <- FENGXIAN_COVAR_NACC[1:10]

# FILL ADGC Phenotype information from Fengxian's phenotype file
cc <- combined[na.omit(match(FENGXIAN_COVAR_NACC$IID,combined$CLEANED_ID)),]
# combined[na.omit(match(FENGXIAN_COVAR_NACC$IID,combined$CLEANED_ID)), c("SEX", "STATUS", "AGE_AT_ONSET", "AGE_LAST_VISIT", "APOE4ANY")]
combined$SEX[na.omit(match(FENGXIAN_COVAR_NACC$IID,combined$CLEANED_ID))] <- FENGXIAN_COVAR_NACC$SEX
combined$STATUS[na.omit(match(FENGXIAN_COVAR_NACC$IID,combined$CLEANED_ID))] <- FENGXIAN_COVAR_NACC$STATUS
combined$AGE_AT_ONSET[na.omit(match(FENGXIAN_COVAR_NACC$IID,combined$CLEANED_ID))] <- FENGXIAN_COVAR_NACC$AGE_AT_ONSET
combined$AGE_LAST_VISIT[na.omit(match(FENGXIAN_COVAR_NACC$IID,combined$CLEANED_ID))] <- FENGXIAN_COVAR_NACC$AGE_LAST_VISIT
combined$APOE4ANY[na.omit(match(FENGXIAN_COVAR_NACC$IID,combined$CLEANED_ID))] <- FENGXIAN_COVAR_NACC$APOE4ANY




# Now get the IID and FID for these from FAM file
combined$original_FID <- as.character(combined$FID)
combined$original_IID <- as.character(combined$IID)
combined$FID <- as.character(ADGC_NHW_FAM$V1[match(combined$KEY, ADGC_NHW_FAM$KEY)])
combined$IID <- as.character(ADGC_NHW_FAM$V2[match(combined$KEY, ADGC_NHW_FAM$KEY)])


###########################
## Recode SEX and STATUS ##
###########################
# recode STATUS
table(combined$STATUS)
combined <- RECODE_CACO_SEX(combined)

# # NACC PHENO with selected STATUS from Fengxian
# NACC <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/WASHU-GWAS/NACC_unique_core_pheno_20210927.csv", sep = ",", header = T)
# # Selecting only the AD samples
# NACC$STATUS <- as.character(NACC$final_CC_status)
# NACC <- NACC[grepl("^CA$|^CO$|^Neuro_AD$|^Neuro_AD_DLB$|^Neuro_CO$|^Neuro_PreSymptomatic_AD$|OT\\(CO\\)", NACC$STATUS),]
# table(NACC$STATUS)
# 
# # Recode CA/CO
# NACC$STATUS[grepl("^CO$|^Neuro_CO$|OT\\(CO\\)", NACC$STATUS)] <- 1
# NACC$STATUS[grepl("^CA$|^Neuro_AD$|^Neuro_AD_DLB$|^Neuro_PreSymptomatic_AD$", NACC$STATUS)] <- 2
# table(NACC$STATUS)
# 
# NACC <- change_names(NACC)
# NACC$SEX <- as.character(NACC$SEX)
# # Recode SEX
# NACC$SEX[NACC$SEX == "Female"] <- 2
# NACC$SEX[NACC$SEX == "Male"] <- 1

head(ADGC_NHW_FAM)

combined_NHW <- combined[grepl("YES",combined$FOUND),]
combined_NHW$ETHNICITY <- "NHW"

write.table(combined_NHW, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/cleaned_phenotypes/CLEANED_PHENO_ADGC_NHW_49586.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)


################################################################
########################## ADD PCA #############################
################################################################
# PCA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/01-EOAD-preQC/NHW_analysis/analysis_for_early_onset_locus_rs143080277/ADGC_NHW_Cohort-PCAS-ALL.eigenvec", header = T)
# PCA$FID <- as.character(PCA$FID)
# PCA$IID <- as.character(PCA$IID)
# PCA$PC_KEY <- paste(PCA$FID, PCA$IID, sep = ":")
# colnames(PCA)[1:2] <- c(paste0("PC_", colnames(PCA)[1:2] ))
# 
# combined_NHW$KEY2 <- paste(combined_NHW$FID, combined_NHW$IID, sep = ":")
# # Keep only the required columns
# combined_NHW <- combined_NHW[c("FID", "IID", "SEX", "STATUS", "COHORT", "AGE_AT_ONSET", "AGE_LAST_VISIT", "AGE_AT_DEATH", "APOE4ANY", "RACE", "ETHNICITY", "CLEANED_ID", "FOUND", "KEY", "KEY2")]
# 
# # Model the STUDY
# combined_NHW$STUDY <- as.factor(combined_NHW$COHORT)
# CLEAN <- function(x){ colnames(x) <- gsub("STUDY", "", colnames(x)); x } 
# combined_NHW <- cbind(combined_NHW, CLEAN(as.data.frame(model.matrix(~STUDY, data=combined_NHW))))
# 
# 
# 
# # Now merge PCAs
# combined_NHW <- cbind(combined_NHW, PCA[match(combined_NHW$KEY2, PCA$PC_KEY),])
# 
# # Age threshold for CA <=65  & CO 70 
# combined_NHW_CO70 <- combined_NHW[which((combined_NHW$STATUS == "2" & combined_NHW$AGE_AT_ONSET <= 65)| (combined_NHW$STATUS == "1" & combined_NHW$AGE_LAST_VISIT > 70)),]
# table(combined_NHW_CO70$STATUS)
# # 1     2 
# # 14125  3202
# 
# write.table(combined_NHW_CO70, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/cleaned_phenotypes/CLEANED_PHENO_ADGC_NHW_CO70.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)
# write.table(combined_NHW_CO70[1:2], "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/cleaned_phenotypes/CLEANED_PHENO_ADGC_NHW_CO70_IDs.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)
# 
# # Age threshold for CA <=65  & CO 80 
# combined_NHW_CO80 <- combined_NHW[which((combined_NHW$STATUS == "2" & combined_NHW$AGE_AT_ONSET <= 65)| (combined_NHW$STATUS == "1" & combined_NHW$AGE_LAST_VISIT > 80)),]
# table(combined_NHW_CO80$STATUS)
# # 1    2 
# # 5782 3202
# write.table(combined_NHW_CO80, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/cleaned_phenotypes/CLEANED_PHENO_ADGC_NHW_CO80.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)
# write.table(combined_NHW_CO80[1:2], "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/cleaned_phenotypes/CLEANED_PHENO_ADGC_NHW_CO80_IDs.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)
# 
# 
# ##################################################### For report
# tt <- as.data.frame(table(ADGC_NHW_FAM$COHORT))
# paste0(tt[,1]," (",  tt[,2], "),")
# names(table(combined_NHW_CO70$COHORT))
# 
# dim(combined_NHW_CO70)
# # [1] 17327    77
# dim(combined_NHW_CO80)
# # [1] 8984   77
# table(combined_NHW_CO70$STATUS)
# # 1     2 
# # 14125  3202 
# table(combined_NHW_CO80$STATUS)
# # 1    2 
# # 5782 3202 



##########################
##########################
########### AA ###########
##########################
##########################
# ADGC AA covariate
ADGC_AA_covar <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/African_American_covariate.tsv", header = T)
ADGC_AA_covar$COHORT <- as.character(ADGC_AA_covar$COHORT)
ADGC_AA_covar$COHORT[ADGC_AA_covar$COHORT =="ADC1_2_AA"] <- "ADC1-2_AA"

## Identify FAM STUDY
# All FAM files
fileName <- "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_AA/fam.list"

## Loop over a file connection
conn <- file(fileName,open="r")
linn <-readLines(conn)
df1 <- {}
for (i in 1:length(linn)){
  print(linn[i])
  COHORT <- str_match(linn[i], "ADGC_AA-\\s*(.*?)\\s*_geno") [,2] 
  df.tmp <- read.table(linn[i], header = F, colClasses = 'character')  
  df.tmp$COHORT <- COHORT
  df1 <- rbind(df1,df.tmp)
  # Sys.sleep(2)
}

close(conn)


df1$V2 <- as.character(df1$V2)
# df1 <- df1[!duplicated(df1$V2),]
df1$KEY <- paste(df1$COHORT, df1$V2, sep = ":")

ADGC_AA_FAM <- df1

head(ADGC_AA_covar)
ADGC_AA_covar$MERGED_ID <- as.character(ADGC_AA_covar$MERGED_ID)


# Now create a column for match and mismatch 
ADGC_AA_covar$CLEANED_ID <- as.character(ADGC_AA_covar$IID)
# CHECK MATCHES !!
ADGC_AA_covar$FOUND <- ifelse(ADGC_AA_covar$CLEANED_ID %in% ADGC_AA_FAM$V2, "YES", "NO")
table(ADGC_AA_covar$FOUND)
# NO  YES 
# 1180 7383

################
## Clean NACC ##
################
# IDs in AA covar : 09AD14028_NACC367443

ADGC_AA_covar$CLEANED_ID[grepl("[[:digit:]]+AD+[[:digit:]]+_NACC+[[:digit:]]", ADGC_AA_covar$CLEANED_ID)] <- paste0("0_",ADGC_AA_covar$CLEANED_ID[grepl("[[:digit:]]+AD+[[:digit:]]+_NACC+[[:digit:]]", ADGC_AA_covar$CLEANED_ID)])
# CHECK MATCHES !!
ADGC_AA_covar$FOUND <- ifelse(ADGC_AA_covar$CLEANED_ID %in% ADGC_AA_FAM$V2, "YES", "NO")
table(ADGC_AA_covar$FOUND)
# NO  YES 
# 1102 7461

# Clean NACC229788_09AD13003
ADGC_AA_covar$CLEANED_ID[grepl("NACC+[[:digit:]]+_+[[:digit:]]+AD+[[:digit:]]", ADGC_AA_covar$CLEANED_ID)] <- paste0("0_",ADGC_AA_covar$CLEANED_ID[grepl("NACC+[[:digit:]]+_+[[:digit:]]+AD+[[:digit:]]", ADGC_AA_covar$CLEANED_ID)])
# CHECK MATCHES !!
ADGC_AA_covar$FOUND <- ifelse(ADGC_AA_covar$CLEANED_ID %in% ADGC_AA_FAM$V2, "YES", "NO")
table(ADGC_AA_covar$FOUND)
# NO  YES 
# 1041 7522



########################
## Clean MIRAGE600_AA ##
########################

ADGC_AA_covar$CLEANED_ID[ADGC_AA_covar$COHORT == "MIRAGE600_AA"] <- sapply(strsplit(ADGC_AA_covar$CLEANED_ID[ADGC_AA_covar$COHORT == "MIRAGE600_AA"] ,"_"), `[`, 2)
# CHECK MATCHES !!
ADGC_AA_covar$FOUND <- ifelse(ADGC_AA_covar$CLEANED_ID %in% ADGC_AA_FAM$V2, "YES", "NO")
table(ADGC_AA_covar$FOUND)
# NO  YES 
# 614 7949


##################
## Clean JHU_AA ##
##################
foo <- ADGC_AA_covar$CLEANED_ID[ADGC_AA_covar$COHORT == "JHU_AA"]
# get indices of IDs in FAM with grep
INDEX <- NULL
for (i in 1:length(foo)){
 INDEX.tmp <-  grep(foo[i], ADGC_AA_FAM$V2)
 INDEX <- c(INDEX, INDEX.tmp) 
}

# SANITY check
ADGC_AA_covar$CLEANED_ID[ADGC_AA_covar$COHORT == "JHU_AA"] == paste(sapply(strsplit(ADGC_AA_FAM$V2[INDEX], "_"), `[`, 2), sapply(strsplit(ADGC_AA_FAM$V2[INDEX], "_"), `[`, 3), sep = "_")

# Now IDs from FAM to phenotype
ADGC_AA_covar$CLEANED_ID[ADGC_AA_covar$COHORT == "JHU_AA"] <-   ADGC_AA_FAM$V2[INDEX]
# CHECK MATCHES !!
ADGC_AA_covar$FOUND <- ifelse(ADGC_AA_covar$CLEANED_ID %in% ADGC_AA_FAM$V2, "YES", "NO")
table(ADGC_AA_covar$FOUND)
# NO  YES 
# 116 8447


########################
## Clean MIRAGE600_AA ##
########################

ADGC_AA_covar$CLEANED_ID[ADGC_AA_covar$COHORT == "MIRAGE300_AA"] <- sapply(strsplit(ADGC_AA_covar$CLEANED_ID[ADGC_AA_covar$COHORT == "MIRAGE300_AA"] ,"_"), `[`, 2)
# CHECK MATCHES !!
ADGC_AA_covar$FOUND <- ifelse(ADGC_AA_covar$CLEANED_ID %in% ADGC_AA_FAM$V2, "YES", "NO")
table(ADGC_AA_covar$FOUND)
# YES 
# 8563

# Now that I was able to match all PHENO IDs with GENO IDs, I will create a LOOKUP column to MATCH FID from FAM file
ADGC_AA_covar$KEY <- paste(ADGC_AA_covar$COHORT, ADGC_AA_covar$CLEANED_ID, sep = ":")
# Let's see if this KEY has any duplicates
ADGC_AA_covar$KEY[duplicated(ADGC_AA_covar$KEY)]

# Since MIRAGE samples are duplicated, I will clean them again. Before that, let's copy the FID from FAM files for the unique samples
ADGC_AA_covar$FID <- as.character(ADGC_AA_covar$FID)
ADGC_AA_covar$IID <- as.character(ADGC_AA_covar$IID)
ADGC_AA_covar$original_FID <- ADGC_AA_covar$FID
ADGC_AA_covar$original_IID <- ADGC_AA_covar$IID
ADGC_AA_covar$FID <- ADGC_AA_FAM$V1[match(ADGC_AA_covar$KEY, ADGC_AA_FAM$KEY)]

# Now CLEAN MIRAGE IDs and FIDs again
ADGC_AA_covar$FID[grepl("MIRAGE", ADGC_AA_covar$COHORT)] <- sapply(strsplit(ADGC_AA_covar$original_FID[grepl("MIRAGE", ADGC_AA_covar$COHORT)],"_"), `[`, 1)
ADGC_AA_covar$CLEANED_ID[grepl("MIRAGE", ADGC_AA_covar$COHORT)] <- sapply(strsplit(ADGC_AA_covar$original_FID[grepl("MIRAGE", ADGC_AA_covar$COHORT)],"_"), `[`, 2)


# Now cleaned ID becomes IID
ADGC_AA_covar$IID <- ADGC_AA_covar$CLEANED_ID

# Now create KEY2 with FID and IID to cross check the samples in PHENO and GENO datasets
ADGC_AA_covar$KEY2 <- paste(ADGC_AA_covar$FID, ADGC_AA_covar$IID, sep = ":")
ADGC_AA_FAM$KEY2 <- paste(ADGC_AA_FAM$V1, ADGC_AA_FAM$V2, sep = ":")

# CHECK if KEY2 have any duplicates
sum(duplicated(ADGC_AA_covar$KEY2))
# 0

sum(ADGC_AA_covar$KEY2 %in% ADGC_AA_FAM$KEY2)
# 8563

## Recode STATUS and SEX
ADGC_AA_covar <- RECODE_CACO_SEX(ADGC_AA_covar)


write.table(ADGC_AA_covar, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/cleaned_phenotypes/CLEANED_PHENO_ADGC_AA_8563.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)

# ###########################
# ## Recode SEX and STATUS ##
# ###########################
# # recode STATUS
# table(combined$STATUS)
# combined$STATUS <- as.character(combined$STATUS)
# combined$STATUS[combined$STATUS=="CA"] <- 2
# combined$STATUS[combined$STATUS=="CO"] <- 1
# combined$STATUS[combined$STATUS=="MCI"] <- 3
# combined$STATUS[combined$STATUS=="unknown"] <- -9
# combined$SEX[combined$SEX == "Female"] <- 2
# combined$SEX[combined$SEX == "Male"] <- 1
# 
# 
# combined_NHW <- combined[grepl("YES",combined$FOUND),]
# 
# write.table(combined_NHW, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/cleaned_phenotypes/CLEANED_PHENO_ADGC_NHW_49586.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)
# 
# 
# 
# ######################## Merge with PCA
# PCA <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/01-EOAD-preQC/NHW_analysis/analysis_for_early_onset_locus_rs143080277/ADGC_NHW_Cohort-PCAS-ALL.eigenvec", header = T)
# PCA$FID <- as.character(PCA$FID)
# PCA$IID <- as.character(PCA$IID)
# PCA$PC_KEY <- paste(PCA$FID, PCA$IID, sep = ":")
# colnames(PCA)[1:2] <- c(paste0("PC_", colnames(PCA)[1:2] ))
# 
# combined_NHW$KEY2 <- paste(combined_NHW$FID, combined_NHW$IID, sep = ":")
# # Keep only the required columns
# combined_NHW <- combined_NHW[c("FID", "IID", "SEX", "STATUS", "COHORT", "AGE_AT_ONSET", "AGE_LAST_VISIT", "AGE_AT_DEATH", "APOE4ANY", "RACE", "ETHNICITY", "CLEANED_ID", "FOUND", "KEY", "KEY2")]
# 
# # Model the STUDY
# combined_NHW$STUDY <- as.factor(combined_NHW$COHORT)
# CLEAN <- function(x){ colnames(x) <- gsub("STUDY", "", colnames(x)); x } 
# combined_NHW <- cbind(combined_NHW, CLEAN(as.data.frame(model.matrix(~STUDY, data=combined_NHW))))
# 
# 
# 
# # Now merge PCAs
# combined_NHW <- cbind(combined_NHW, PCA[match(combined_NHW$KEY2, PCA$PC_KEY),])
# 
# # Age threshold for CA <=65  & CO 70 
# combined_NHW_CO70 <- combined_NHW[which((combined_NHW$STATUS == "2" & combined_NHW$AGE_AT_ONSET <= 65)| (combined_NHW$STATUS == "1" & combined_NHW$AGE_LAST_VISIT > 70)),]
# table(combined_NHW_CO70$STATUS)
# # 1     2 
# # 14125  3202
# 
# write.table(combined_NHW_CO70, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/cleaned_phenotypes/CLEANED_PHENO_ADGC_NHW_CO70.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)
# write.table(combined_NHW_CO70[1:2], "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/cleaned_phenotypes/CLEANED_PHENO_ADGC_NHW_CO70_IDs.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)
# 
# # Age threshold for CA <=65  & CO 80 
# combined_NHW_CO80 <- combined_NHW[which((combined_NHW$STATUS == "2" & combined_NHW$AGE_AT_ONSET <= 65)| (combined_NHW$STATUS == "1" & combined_NHW$AGE_LAST_VISIT > 80)),]
# table(combined_NHW_CO80$STATUS)
# # 1    2 
# # 5782 3202
# write.table(combined_NHW_CO80, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/cleaned_phenotypes/CLEANED_PHENO_ADGC_NHW_CO80.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)
# write.table(combined_NHW_CO80[1:2], "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/cleaned_phenotypes/CLEANED_PHENO_ADGC_NHW_CO80_IDs.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)
# 

##########################
##########################
######### Asian ##########
##########################
##########################

# ADGC AA covariate
ADGC_Asian_covar <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ASIAN_covariate.tsv", header = T)
ADGC_Asian_covar$COHORT  <- as.character(ADGC_Asian_covar$COHORT )
ADGC_Asian_covar$COHORT [ADGC_Asian_covar$COHORT == "ASA_JPN"] <- "ASA-JPN"

ADGC_Asian_covar$COHORT <- as.character(ADGC_Asian_covar$COHORT)

## Identify FAM STUDY
# All FAM files
fileName <- "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Asian/fam.list"

## Loop over a file connection
conn <- file(fileName,open="r")
linn <-readLines(conn)
df1 <- {}
for (i in 1:length(linn)){
  print(linn[i])
  COHORT <- str_match(linn[i], "ADGC_Asian-\\s*(.*?)\\s*_geno") [,2] 
  df.tmp <- read.table(linn[i], header = F, colClasses = 'character')  
  df.tmp$COHORT <- COHORT
  df1 <- rbind(df1,df.tmp)
  # Sys.sleep(2)
}

close(conn)


df1$V2 <- as.character(df1$V2)
# df1 <- df1[!duplicated(df1$V2),]
df1$KEY <- paste(df1$COHORT, df1$V2, sep = ":")

ADGC_Asian_FAM <- df1

head(ADGC_Asian_covar)
ADGC_Asian_covar$MERGED_ID <- as.character(ADGC_Asian_covar$MERGED_ID)


# Now create a column for match and mismatch 
ADGC_Asian_covar$CLEANED_ID <- as.character(ADGC_Asian_covar$IID)
# CHECK MATCHES !!
ADGC_Asian_covar$FOUND <- ifelse(ADGC_Asian_covar$CLEANED_ID %in% ADGC_Asian_FAM$V2, "YES", "NO")
table(ADGC_Asian_covar$FOUND)
# NO  YES 
# 1753 2989 

NOTFOUND <- ADGC_Asian_covar[ADGC_Asian_covar$FOUND == "NO",]

###############
## Clean JPN ##
###############
# IDs AF_002_080108 with 0_AF_002_080108

ADGC_Asian_covar$CLEANED_ID[ADGC_Asian_covar$COHORT == "JPN"] <- paste0("0_", as.character(ADGC_Asian_covar$IID[ADGC_Asian_covar$COHORT == "JPN"]))

ADGC_Asian_covar$FOUND <- ifelse(ADGC_Asian_covar$CLEANED_ID %in% ADGC_Asian_FAM$V2, "YES", "NO")
table(ADGC_Asian_covar$FOUND)
# YES 
# 4742 

# Check for duplicate IDs
ADGC_Asian_covar$CLEANED_ID[duplicated(ADGC_Asian_covar$CLEANED_ID)]

# Create KEY by COHORT
ADGC_Asian_covar$KEY <- paste(ADGC_Asian_covar$COHORT, ADGC_Asian_covar$CLEANED_ID, sep = ":")
sum(duplicated(ADGC_Asian_covar$KEY))
# 0

# Now extract FID from FAM
ADGC_Asian_covar$original_FID <- ADGC_Asian_covar$FID
ADGC_Asian_covar$original_IID <- ADGC_Asian_covar$IID

ADGC_Asian_covar$FID <- ADGC_Asian_FAM$V1[match(ADGC_Asian_covar$KEY, ADGC_Asian_FAM$KEY)]
ADGC_Asian_covar$IID  <- ADGC_Asian_covar$CLEANED_ID

# Finally, check if KEY by FAM and IID match in both PHENO and GENO datasets
ADGC_Asian_covar$KEY2 <- paste(ADGC_Asian_covar$FID, ADGC_Asian_covar$IID, sep = ":")
ADGC_Asian_FAM$KEY2 <- paste(ADGC_Asian_FAM$V1, ADGC_Asian_FAM$V2, sep = ":")
sum(!ADGC_Asian_covar$KEY2 %in% ADGC_Asian_FAM$KEY2)
# 0 

## recode STATUS and SEX
ADGC_Asian_covar <- RECODE_CACO_SEX(ADGC_Asian_covar)

write.table(ADGC_Asian_covar, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/cleaned_phenotypes/CLEANED_PHENO_ADGC_Asian_4742.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)


##########################
##########################
####### Hispanic #########
##########################
##########################

# ADGC Hispanic covariate
ADGC_Hispanic_covar <- read.table("/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/HISPANIC_covariate.tsv", header = T)
ADGC_Hispanic_covar$COHORT  <- as.character(ADGC_Hispanic_covar$COHORT )

ADGC_Hispanic_covar$COHORT <- as.character(ADGC_Hispanic_covar$COHORT)

## Identify FAM STUDY
# All FAM files
fileName <- "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/ADGC_Hispanic/fam.list"

## Loop over a file connection
conn <- file(fileName,open="r")
linn <-readLines(conn)
df1 <- {}
for (i in 1:length(linn)){
  print(linn[i])
  COHORT <- str_match(linn[i], "ADGC_Hispanic-\\s*(.*?)\\s*_geno") [,2] 
  df.tmp <- read.table(linn[i], header = F, colClasses = 'character')  
  df.tmp$COHORT <- COHORT
  df1 <- rbind(df1,df.tmp)
  # Sys.sleep(2)
}

close(conn)


df1$V2 <- as.character(df1$V2)
# df1 <- df1[!duplicated(df1$V2),]
df1$KEY <- paste(df1$COHORT, df1$V2, sep = ":")

ADGC_Hispanic_FAM <- df1

head(ADGC_Hispanic_covar)
ADGC_Hispanic_covar$MERGED_ID <- as.character(ADGC_Hispanic_covar$MERGED_ID)


# Now create a column for match and mismatch 
ADGC_Hispanic_covar$CLEANED_ID <- as.character(ADGC_Hispanic_covar$IID)
# CHECK MATCHES !!
ADGC_Hispanic_covar$FOUND <- ifelse(ADGC_Hispanic_covar$CLEANED_ID %in% ADGC_Hispanic_FAM$V2, "YES", "NO")
table(ADGC_Hispanic_covar$FOUND)
# NO  YES 
# 975 1317

NOTFOUND <- ADGC_Hispanic_covar[ADGC_Hispanic_covar$FOUND == "NO",]

#####################
## TARCC3_Hispanic ##
#####################

# IDs 082297Sm with paste(FID, IID, sep = "_")
ADGC_Hispanic_covar$FID <- as.character(ADGC_Hispanic_covar$FID)
ADGC_Hispanic_covar$IID <- as.character(ADGC_Hispanic_covar$IID)

ADGC_Hispanic_covar$CLEANED_ID[ADGC_Hispanic_covar$COHORT == "TARCC3_Hispanic"] <- paste(ADGC_Hispanic_covar$FID, ADGC_Hispanic_covar$IID, sep = "_")[ADGC_Hispanic_covar$COHORT == "TARCC3_Hispanic"]

ADGC_Hispanic_covar$FOUND <- ifelse(ADGC_Hispanic_covar$CLEANED_ID %in% ADGC_Hispanic_FAM$V2, "YES", "NO")
table(ADGC_Hispanic_covar$FOUND)
# YES 
# 2292

# Check for duplicate IDs
length(ADGC_Hispanic_covar$CLEANED_ID[duplicated(ADGC_Hispanic_covar$CLEANED_ID)])
# 40 IIDs are duplicated

# Create KEY by COHORT
ADGC_Hispanic_covar$KEY <- paste(ADGC_Hispanic_covar$COHORT, ADGC_Hispanic_covar$CLEANED_ID, sep = ":")
sum(duplicated(ADGC_Hispanic_covar$KEY))
# 0

# Now extract FID from FAM
ADGC_Hispanic_covar$original_FID <- ADGC_Hispanic_covar$FID
ADGC_Hispanic_covar$original_IID <- ADGC_Hispanic_covar$IID

ADGC_Hispanic_covar$FID <- ADGC_Hispanic_FAM$V1[match(ADGC_Hispanic_covar$KEY, ADGC_Hispanic_FAM$KEY)]
ADGC_Hispanic_covar$IID  <- ADGC_Hispanic_covar$CLEANED_ID

# Finally, check if KEY by FAM and IID match in both PHENO and GENO datasets
ADGC_Hispanic_covar$KEY2 <- paste(ADGC_Hispanic_covar$FID, ADGC_Hispanic_covar$IID, sep = ":")
ADGC_Hispanic_FAM$KEY2 <- paste(ADGC_Hispanic_FAM$V1, ADGC_Hispanic_FAM$V2, sep = ":")
sum(!ADGC_Hispanic_covar$KEY2 %in% ADGC_Hispanic_FAM$KEY2)
# 0 

## recode STATUS and SEX
ADGC_Hispanic_covar <- RECODE_CACO_SEX(ADGC_Hispanic_covar)

write.table(ADGC_Hispanic_covar, "/40/AD/GWAS_data/Source_Plink/2021_ADGC_EOAD/all_covariates/cleaned_phenotypes/CLEANED_PHENO_ADGC_Hispanic_2292.txt", sep ="\t", col.names = T, quote = F, row.names = FALSE)


############################################################################################
####################### Get Demographic table for cleaned covar data #######################
############################################################################################

Get_STATs <- function(covars){
  covars$ETHNICITY <- as.character(covars$ETHNICITY)
  TOTAL = nrow(covars)
  # CA,CO, MCI
  N.controls <- sum(covars$STATUS==1)
  N.cases <- sum(covars$STATUS==2)
  N.mci <- sum(covars$STATUS==3)
  
  # Number of CA <= 65 and 70 (Cases==2)
  CA.65 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == 2,"AGE_AT_ONSET"])))) <= 65)
  CA.70 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == 2,"AGE_AT_ONSET"])))) <= 70)
  # Number of CO > 70 (Cases==2)
  CO.70 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == 1,"AGE_LAST_VISIT"])))) > 70)
  CO.80 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == 1,"AGE_LAST_VISIT"])))) > 80)
  # Number of MCI <= 65 and 70 (Cases==2)
  MCI.65 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == 3,"AGE_AT_ONSET"])))) <= 65)
  MCI.70 <- sum(as.vector(na.omit(as.numeric(as.character(covars[covars$STATUS == 3,"AGE_AT_ONSET"])))) <= 70)
  # Number of Others
  N.OTHERS <- sum(covars$STATUS == -9)
  
  # Percent Female
  PERC.FEMALE <- (sum(covars$SEX == 2, na.rm = T)/(sum(covars$SEX == 1, na.rm = T)+ sum(covars$SEX == 2, na.rm = T))) *100
  
  # MISSING AGES CA and CO
  N.CA.missing.age <- sum(is.na(covars [covars$STATUS == 2, "AGE_AT_ONSET"]))
  N.CO.missing.age <- sum(is.na(covars [covars$STATUS == 1, "AGE_LAST_VISIT"]))
  
  
  # Percent APOE4
  # PERC.APOE <- (table( covars[, c("APOE4ANY") ] )[3]/(table( covars[, c("APOE4ANY") ] )[2] + table( covars[, c("APOE4ANY") ] )[3]))*100 
  POS <- sum(covars[, c("APOE4ANY")] == 1, na.rm = T)
  NEG <- sum(covars[, c("APOE4ANY")] == 0, na.rm = T)
  UNK <- sum(covars[, c("APOE4ANY")] == -9, na.rm = T)
  PERC.APOE <- (POS/ (POS + NEG)) *100
  
  CA.APOE.PERC <- sum(covars[ covars$STATUS == 2 , c("APOE4ANY")] == 1, na.rm = T)/ (N.cases) *100
  CO.APOE.PERC <- sum(covars[ covars$STATUS == 1 , c("APOE4ANY")] == 1, na.rm = T)/ (N.controls) *100
  MCI.APOE.PERC <- sum(covars[ covars$STATUS == 3 , c("APOE4ANY")] == 1, na.rm = T)/ (N.mci) *100
  
  # For Ethnicity
  STATS <- cbind(ETHNICITY = unique(covars$ETHNICITY), TOTAL = TOTAL, '% FEMALE' = round(PERC.FEMALE, 2), '% APOE' = round(PERC.APOE, 2), 'N CONTROLS (1)' = N.controls, 'N CASES (2)' = N.cases,
                 'N MCI (3)' = N.mci, 'N CONTROLS > 70 yo' = CO.70, 'N CONTROLS > 80 yo' = CO.80, 'N CONTROLS missing age' = N.CO.missing.age, 'N CASES ≤ 65 yo' = CA.65, 'N CASES ≤ 70 yo' = CA.70,
                 'N CASES missing age' = N.CA.missing.age, 'N MCI (3) ≤ 65 yo' = MCI.65, 'N MCI (3) ≤ 70 yo' = MCI.70, 'N OTHERS (-9)' = N.OTHERS, '% CONTROLS (1) APOE4+' = round(CO.APOE.PERC, 2),
                 '% CASES (2) APOE4+' = round(CA.APOE.PERC, 2), '% MCI (3) APOE4+' = round(MCI.APOE.PERC, 2))
  
  return(STATS)
}





View(as.data.frame(t(Get_STATs(ADGC_AA_covar))))
View(as.data.frame(t(Get_STATs(ADGC_Asian_covar))))
View(as.data.frame(t(Get_STATs(ADGC_Hispanic_covar))))
View(as.data.frame(t(Get_STATs(combined_NHW))))

unname(cbind.data.frame(as.data.frame(t(Get_STATs(ADGC_AA_covar))), as.data.frame(t(Get_STATs(ADGC_Asian_covar))), as.data.frame(t(Get_STATs(ADGC_Hispanic_covar))), as.data.frame(t(Get_STATs(combined_NHW)))))

####################################################################################################
###################################### HapMap All Ethnicities ######################################
####################################################################################################














#############################################################################################################################################################################################################################
#############################################################################################################################################################################################################################
#############################################################################################################################################################################################################################
#############################################################################################################################################################################################################################
#############################################################################################################################################################################################################################



