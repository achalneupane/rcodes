setwd("/40/AD/GWAS_data/Source_Plink/2020_ADGC_EOAD/Achal_ADGC_EOAD-analyses")

#install.packages("qqman",repos="http://cran.cnr.berkeley.edu/",lib="~" ) # location of installation can be changed but has to correspond with the library location 
library("qqman") 
results_log <- read.table("ADGC_Asian-logistic-OR-CI.assoc.logistic", head=TRUE)
jpeg("ADGC_Asian_QQ-Plot_logistic.jpeg")
qq(results_log$P, main = "ADGC_Asian Q-Q plot of GWAS p-values : log")
dev.off()

#results_as <- read.table("assoc_results.assoc", head=TRUE)
#jpeg("ADGC_AA_QQ-Plot_assoc.jpeg")
#qq(results_as$P, main = "ADGC_AA Q-Q plot of GWAS p-values : log")
#dev.off()

