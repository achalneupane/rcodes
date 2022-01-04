setwd("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/06-Aquilla_202101-a/01-Aquilla-preQC/03-PLINK-QC-files4/ADMIXTURE")
df <- read.table("cv_error.txt", header = T)
library(ggplot2)
ggplot(data=df, aes(x=K, y=Error)) + geom_line()

ggsave("ADMIXTURE-CV-FASE.jpg", plot = last_plot(), device = NULL, scale = 1, width = 16, height = 9, dpi = 300, limitsize = TRUE)


FASE_Ethnicity_from_Aurora <- read.table("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/06-Aquilla_202101-a/01-Aquilla-preQC/03-PLINK-QC-files4/GENESIS/FASe_3894_Ethnicity_from_Aurora.txt", header = T)

FAM <- read.table("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1.fam", header = F)

FAM <- cbind(FAM,FASE_Ethnicity_from_Aurora[match(FAM$V2, FASE_Ethnicity_from_Aurora$IID),])
sum(as.character(FAM$V2) == as.character(FAM$IID))

# plot the Q estimates
tbl=read.table("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1.4.Q")

tbl <- cbind(IID= FAM$IID, Ethnicity=FAM$Population,tbl)
rownames(tbl) <- tbl$IID
library(tidyr)

plot_data <- tbl %>% 
  mutate(id = row_number())%>% 
  gather('pop', 'prob', V1:V4) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

plot_data$Ethnicity <- as.character(plot_data$Ethnicity)
plot_data$Ethnicity[(is.na(plot_data$Ethnicity))] <- "-9"

## With facects
p <- ggplot(plot_data, aes(IID, prob, fill = pop)) +
  geom_col() +
  facet_grid(~likely_assignment, scales = 'free', space = 'free')

p

ggsave("ADMIXTURE-population-k4.jpg", plot = p, device = NULL, scale = 1, width = 16, height = 9, dpi = 300, limitsize = TRUE)

p <- ggplot(plot_data, aes(Ethnicity, prob, fill = pop)) +
  geom_col() +
  facet_grid(~likely_assignment, scales = 'free', space = 'free')

p
ggsave("ADMIXTURE-population--Aurora-assignment-k4.jpg", plot = p, device = NULL, scale = 1, width = 16, height = 9, dpi = 300, limitsize = TRUE)

## Theme classic
ggplot(plot_data, aes(IID, prob, fill = pop)) +
  geom_col() +
  theme_classic()



# TANZI <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/FBAT/FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1_with_STATUS_nonADSP_post_QC2-geno-0.02-maxmaf-0.01.fam")
