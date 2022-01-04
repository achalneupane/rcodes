setwd("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files_all_Aquilla/variant_check")
PHENO <- read.csv("/40/AD/AD_Seq_Data/03.-phenotype/2021-01-Aquilla-phenotype/Joint_Call_Round1/Aquilla_9576_phenotype.csv")
dim(PHENO)
# Case=2; control=1
table(PHENO$STATUS)


tt1 <- read.table("Aquilla_variant_grep_X_135974614.vcf")
tt2 <- read.table("Aquilla_variant_grep_X_135974666.vcf")
tt3 <- read.table("Aquilla_variant_grep_X_135974676.vcf")
tt4 <- read.table("Aquilla_variant_grep_X_135974680.vcf")
tt5 <- read.table("Aquilla_variant_grep_X_135974739.vcf")
tt6 <- read.table("Aquilla_variant_grep_X_135974747.vcf")
tt7 <- read.table("Aquilla_variant_grep_X_135974755.vcf")
tt8 <- read.table("Aquilla_variant_grep_X_135974835.vcf")
tt9 <- read.table("Aquilla_variant_grep_X_135974867.vcf")
tt10 <- read.table("Aquilla_variant_grep_X_135974800.vcf")
tt11 <- read.table("Aquilla_variant_grep_X_135974834.vcf")


VARIANTS <- c("X:135974614:A:G", "X:135974666:C:T", "X:135974676:A:G", "X:135974680:G:A", "X:135974739:C:T",
              "X:135974747:C:T", "X:135974755:G:A", "X:135974835:T:G", "X:135974867:G:A", "X:135974800:C:A", "X:135974834:T:C")

VAR.ALL <- list()
for(i in 1:11){ 
tt <- eval(parse(text=paste0("tt", i)))
dim(tt)
head(tt)

tt <- tt[-(1:8),]
colnames(tt) <- c("ID", "variant")
head(tt)
tt$STATUS <- PHENO$STATUS[match(tt$ID, PHENO$Aquilla.vcfID)]
tt$variant <- as.character(tt$variant)

tt <- tt[grepl ("0/0|0/1|1/1", tt$variant ),]

# pp <- tt[c(1:5),]
# pp <- rbind(pp, tt[grepl("4_596_8|5_26202_105|A-ADC-AD002860|64368|S0085", tt$ID),])
# pp <- as.matrix(pp)
# pp <- as.data.frame(pp)

pp <- table(tt[c("STATUS", "variant")])
pp <- as.data.frame(cbind(pp, Variant= VARIANTS[i]))
pp$STATUS <- rownames(as.data.frame(cbind(pp, Variant= VARIANTS[i])))
VAR.ALL[[i]] <- pp
}

FINAL <- Reduce(function(x, y) merge(x, y, all=TRUE), VAR.ALL)
FINAL <- FINAL[c(2,3,1,4,5)]

# Sort the dataframe
FINAL <- FINAL[with(FINAL, order(Variant, STATUS)), ]

write.table(FINAL, "/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files_all_Aquilla/variant_check/Slc9_Aquilla_Achal.csv", sep =",", col.names = T, quote = F, row.names = FALSE)
