table <- data.frame()

for( i in 1:length(snps) ) {

   snp <- snps[i]
   snp.data <- gen[,c("IID","PHEN",snp)]
   data <- merge(b27,snp.data,by=1)
   data <- na.omit(data)

   form0 <- as.formula(paste("PHEN ~ B27 + ",snp,sep=''))
   form1 <- as.formula(paste("PHEN ~ B27*",snp,sep=''))

   fit0 <- glm(form0,data=data,family=binomial(logit))
   fit1 <- glm(form1,data=data,family=binomial(logit))
   anv <- anova(fit0,fit1,test='Chisq')

   pval.int <- summary(fit1)$coefficients[4,4]
   pval.anv <- anv$"Pr(>Chi)"[2]

   snp.id <- gsub("_.$","",snp)

   line <- data.frame(SNP=snp.id,
                      PVAL_INT=pval.int,
                      PVAL_ANV=pval.anv)

   table <- rbind(table,line)

}
