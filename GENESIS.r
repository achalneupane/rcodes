# https://bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/pcair.html

cd /100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files4/GENESIS
# https://www.biostars.org/p/104412/

ln -s  ../FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1.fam
ln -s ../FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1.bed
ln -s ../FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1.bim



# Overview
# This vignette provides a description of how to use the r Biocpkg("GENESIS")
# package to run genetic association tests on array (SNP) data. r
# Biocpkg("GENESIS") uses mixed models for genetic association testing, as
# PC-AiR PCs can be used as fixed effect covariates to adjust for population
# stratification, and a kinship matrix (or genetic relationship matrix)
# estimated from PC-Relate can be used to account for phenotype correlation due
# to genetic similarity among samples.

# Data
# Preparing Scan Annotation Data
# The fitNullModel function in the r Biocpkg("GENESIS") package reads sample
# data from either a standard data.frame class object or a
# ScanAnnotationDataFrame class object as created by the r Biocpkg("GWASTools")
# package. This object must contain all of the outcome and covariate data for
# all samples to be included in the mixed model analysis. Additionally, this
# object must include a variable called "scanID" which contains a unique
# identifier for each sample in the analysis. While a standard data.frame can be
# used, we recommend using a ScanAnnotationDataFrame object, as it can be paired
# with the genotype data (see below) to ensure matching of sample phenotype and
# genotype data. Through the use of r Biocpkg("GWASTools"), a
# ScanAnnotationDataFrame class object can easily be created from a data.frame
# class object. Example R code for creating a ScanAnnotationDataFrame object is
# presented below. Much more detail can be found in the r Biocpkg("GWASTools")
# package reference manual.

library(GENESIS)
library(GWASTools)

# # file path to GDS file
# gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
# # read in GDS data
# HapMap_geno <- GdsGenotypeReader(filename = gdsfile)
# # create a GenotypeData class object
# HapMap_genoData <- GenotypeData(HapMap_geno)
# # load saved matrix of KING-robust estimates
# data("HapMap_ASW_MXL_KINGmat")
# 
# # run PC-AiR
# mypcair <- pcair(HapMap_genoData, 
#                  kinobj = HapMap_ASW_MXL_KINGmat, 
#                  divobj = HapMap_ASW_MXL_KINGmat,
#                  verbose = FALSE)
# mypcs <- mypcair$vectors[,1,drop=FALSE]
# 
# # create a GenotypeBlockIterator object
# HapMap_genoData <- GenotypeBlockIterator(HapMap_genoData) 
# # run PC-Relate
# mypcrel <- pcrelate(HapMap_genoData, pcs = mypcs,
#                     training.set = mypcair$unrels,
#                     verbose = FALSE)
# 
# # generate a phenotype
# set.seed(4)
# pheno <- 0.2*mypcs + rnorm(mypcair$nsamp, mean = 0, sd = 1)
# # mypcair contains PCs from a previous PC-AiR analysis
# # pheno is a vector of Phenotype values
# 
# # make a data.frame
# mydat <- data.frame(scanID = mypcair$sample.id, pc1 = mypcair$vectors[,1], 
#                     pheno = pheno)
# head(mydat)
# 
# # make ScanAnnotationDataFrame
# scanAnnot <- ScanAnnotationDataFrame(mydat)
# scanAnnot
# Reading in Genotype Data
# The assocTestSingle function in the r Biocpkg("GENESIS") package reads genotype data from a GenotypeData class object as created by the r Biocpkg("GWASTools") package. Through the use of r Biocpkg("GWASTools"), a GenotypeData class object can easily be created from:
#   
#   an R matrix of SNP genotype data
# a GDS file
# PLINK files
# Example R code for creating a GenotypeData object is presented below. Much more detail can be found in the r Biocpkg("GWASTools") package reference manual.
# 
# R Matrix
# geno <- MatrixGenotypeReader(genotype = genotype, snpID = snpID, 
#                              chromosome = chromosome, position = position, 
#                              scanID = scanID)
# genoData <- GenotypeData(geno)
# genotype is a matrix of genotype values coded as 0 / 1 / 2, where rows index SNPs and columns index samples
# snpID is an integer vector of unique SNP IDs
# chromosome is an integer vector specifying the chromosome of each SNP
# position is an integer vector specifying the position of each SNP
# scanID is a vector of unique individual IDs
# GDS files
# geno <- GdsGenotypeReader(filename = "genotype.gds")
# genoData <- GenotypeData(geno)
# filename is the file path to the GDS object
# PLINK files
# The r Biocpkg("SNPRelate") package provides the snpgdsBED2GDS function to convert binary PLINK files into a GDS file.

snpgdsBED2GDS(bed.fn = "FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1.bed", 
              bim.fn = "FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1.bim", 
              fam.fn = "FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1.fam", 
              out.gdsfn = "genotype.gds")
# bed.fn is the file path to the PLINK .bed file
# bim.fn is the file path to the PLINK .bim file
# fam.fn is the file path to the PLINK .fam file
# out.gdsfn is the file path for the output GDS file
# Once the PLINK files have been converted to a GDS file, then a GenotypeData object can be created as described above.

# HapMap Data
# To demonstrate association testing with the r Biocpkg("GENESIS") package, we
# analyze SNP data from the Mexican Americans in Los Angeles, California (MXL)
# and African American individuals in the southwestern USA (ASW) population
# samples of HapMap 3. Mexican Americans and African Americans have a diverse
# ancestral background, and familial relatives are present in these data.
# Genotype data at a subset of 20K autosomal SNPs for 173 individuals are
# provided as a GDS file.

# read in GDS data

gdsfile <- "/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files4/GENESIS/genotype.gds"
# gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")

# LD pruning
library(SNPRelate)
# read in GDS data
gds <- snpgdsOpen(gdsfile)
snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6, 
                          ld.threshold=sqrt(0.1), verbose=FALSE)
pruned <- unlist(snpset, use.names=FALSE)
length(pruned)
## [1] 3826
head(pruned)
## [1]  6  7 15 17 22 31
snpgdsClose(gds)



## Create Kinship matrix
library(SNPRelate)
gds <- snpgdsOpen(gdsfile)
king <- snpgdsIBDKING(gds)
KINGmat <- kingToMatrix(king)
snpgdsClose(gds)



library(GWASTools)
HapMap_geno <- GdsGenotypeReader(filename = gdsfile)
# create a GenotypeData class object
HapMap_genoData <- GenotypeData(HapMap_geno)
HapMap_genoData


part <- pcairPartition(kinobj = KINGmat, divobj = KINGmat)
head(part$unrels); head(part$rels)

# run PC-AiR on pruned SNPs
mypcair <- pcair(HapMap_genoData, kinobj = KINGmat, divobj = KINGmat,
                 snp.include = pruned)


# The default is to plot PC values as black dots and blue pluses for individuals
# in the “unrelated subset” and “related subsets” respectively. The plotting
# colors and characters, as well as other standard plotting parameters, can be
# changed with the standard input to the plot function.

# plot top 2 PCs
jpeg("PCAIR_FASe_3982_PC1_PC2.jpeg", width = 1000, height = 1000, units = "px", pointsize = 16, quality = 75, bg = "white")
plot(mypcair)
dev.off()





# plot PCs 3 and 4
jpeg("PCAIR_FASe_3982_PC3_PC4.jpeg", width = 1000, height = 1000, units = "px", pointsize = 16, quality = 75, bg = "white")
plot(mypcair, vx = 3, vy = 4)
dev.off()




# Relatedness Estimation Adjusted for Principal Components (PC-Relate)
# 4.1Running PC-Relate
# PC-Relate uses the ancestry representative principal components (PCs)
# calculated from PC-AiR to adjust for the population structure and ancestry of
# individuals in the sample and provide accurate estimates of recent genetic
# relatedness due to family structure. The pcrelate function performs the
# PC-Relate analysis.

# The training.set input of the pcrelate function allows for the specification
# of which samples are used to estimate the ancestry adjustment at each SNP. The
# adjustment tends to perform best when close relatives are excluded from
# training.set, so the individuals in the “unrelated subset” from the PC-AiR
# analysis are typically a good choice. However, when an “unrelated subset” is
# not available, the adjustment still works well when estimated using all
# samples (training.set = NULL).

# In order to run PC-Relate, we first need to create an iterator object to read
# SNPs in blocks. We create the iterator such that only pruned SNPs are returned
# in each block.

# run PC-Relate
HapMap_genoData <- GenotypeBlockIterator(HapMap_genoData, snpInclude=pruned)
mypcrelate <- pcrelate(HapMap_genoData, pcs = mypcair$vectors[,1:2], 
                       training.set = mypcair$unrels)
# genoData is a GenotypeIterator class object

# pcs is a matrix whose columns are the PCs used for the ancestry adjustment

# training.set is a vector of individual IDs specifying which samples are used
# to esimate the ancestry adjustment at each SNP

# If estimates of IBD sharing probabilities are not desired, then setting the
# input ibd.probs = FALSE will speed up the computation.

# 4.2Output from PC-Relate
# The pcrelate function will return an object of class pcrelate, which is a list
# of two data.frames: kinBtwn with pairwise kinship values, and kinSelf with
# inbreeding coefficients.

jpeg("PC_RELATE_FASe_3982_PC3_PC4.jpeg", width = 1000, height = 1000, units = "px", pointsize = 16, quality = 75, bg = "white")
plot(mypcrelate$kinBtwn$k0, mypcrelate$kinBtwn$kin, xlab="k0", ylab="kinship")
dev.off()

# A function is provided for making a genetic relationship matrix (GRM). Using a
# threshold for kinship will create a sparse matrix by setting kinship for pairs
# less than the threshold to 0. This can be useful to reduce memory usage for
# very large sample sizes.

iids <- as.character(getScanID(HapMap_genoData))
pcrelateToMatrix(mypcrelate, sample.include = iids[1:5], thresh = 2^(-11/2), scaleKin = 2)



# https://rdrr.io/bioc/GENESIS/f/vignettes/assoc_test.Rmd

library(GENESIS)
library(Biobase)
library(SeqVarTools)


df <- read.table("/40/AD/AD_Seq_Data/03.-phenotype/2021-01-Aquilla-phenotype/Joint_Call_Round1/Aquilla_9576_phenotype.csv", sep = ",", header = TRUE)
head(df)
df$Aquilla.vcfID <- as.character(df$Aquilla.vcfID)

df1 <- read.table("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files4/GENESIS/FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1.fam", sep = " ", header = FALSE)
dim(df1)
head(df1)
colnames(df1) <- c("FID", "IID", "MID", "PID", "SEX", "STATUS")
df1$IID <- as.character(df1$IID)

sum(df1$IID %in% df$Aquilla.vcfID)
# sum(grepl("_R1$", df1$IID))

df2 <- cbind(df1,df[match(df1$IID, df$Aquilla.vcfID), c("Aquilla.vcfID", "Sex", "STATUS")])


df2_nonNA <- df2[!is.na(df2$Aquilla.vcfID),]
df2_R1R2 <- df2[grepl ("_R1$|_R2$", as.character(df2$IID)),]


df2_R1R2 <- cbind(df2_R1R2[1:6], df[match(gsub("_R1|_R2","", df2_R1R2$IID), df$Aquilla.vcfID), c("Aquilla.vcfID", "Sex", "STATUS")])

head(df2_nonNA)
head(df2_R1R2)

df2 <- rbind(df2_nonNA, df2_R1R2)
df2 <- df2[c(2,8,9)]
df2$Population <- "Empty"
head(df2)
df2 <- df2[c("IID", "Population", "Sex", "STATUS")]

annot <- df2
rownames(annot) <- NULL

annot$IID <- as.character(annot$IID)
annot$Sex <- as.character(annot$Sex)
annot$STATUS <- as.character(annot$STATUS)
annot$Population <- as.character(annot$Population)

# # reorder rows
annot <- annot[match(mypcair$sample.id, annot$IID),]

# Should be 70 of these
sum(grepl("_R1$|_R2$", annot$IID))


# make a data.frame
mydat <- data.frame(scanID = mypcair$sample.id, pc1 = mypcair$vectors[,1:3], 
                    pheno = annot)


head(mydat)

# make ScanAnnotationDataFrame
scanAnnot <- ScanAnnotationDataFrame(mydat)
scanAnnot


# https://bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/pcair.html#:~:text=This%20vignette%20provides%20a%20description,sharing%20probabilities%2C%20and%20inbreeding%20coefficients.








##################################### END ##########################################
####################################################################################





library(GENESIS)
library(Biobase)
library(SeqVarTools)


df <- read.table("/40/AD/AD_Seq_Data/03.-phenotype/2021-01-Aquilla-phenotype/Joint_Call_Round1/Aquilla_9576_phenotype.csv", sep = ",", header = TRUE)
head(df)
df$Aquilla.vcfID <- as.character(df$Aquilla.vcfID)

df1 <- read.table("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files4/GENESIS/FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic-hwe-geno0.05-mind0.1.fam", sep = " ", header = FALSE)
dim(df1)
head(df1)
colnames(df1) <- c("FID", "IID", "MID", "PID", "SEX", "STATUS")
df1$IID <- as.character(df1$IID)

sum(df1$IID %in% df$Aquilla.vcfID)
# sum(grepl("_R1$", df1$IID))

df2 <- cbind(df1,df[match(df1$IID, df$Aquilla.vcfID), c("Aquilla.vcfID", "Sex", "STATUS")])


df2_nonNA <- df2[!is.na(df2$Aquilla.vcfID),]
df2_R1R2 <- df2[grepl ("_R1$|_R2$", as.character(df2$IID)),]


df2_R1R2 <- cbind(df2_R1R2[1:6], df[match(gsub("_R1|_R2","", df2_R1R2$IID), df$Aquilla.vcfID), c("Aquilla.vcfID", "Sex", "STATUS")])

head(df2_nonNA)
head(df2_R1R2)

df2 <- rbind(df2_nonNA, df2_R1R2)
df2 <- df2[c(2,8,9)]
df2$Population <- "Empty"
head(df2)
df2 <- df2[c("IID", "Population", "Sex", "STATUS")]


# data(sample_annotation_1KG)
# annot <- sample_annotation_1KG
# head(annot)
##   sample.id Population sex
## 1   HG00110        GBR   F
## 2   HG00116        GBR   M
## 3   HG00120        GBR   F
## 4   HG00128        GBR   F
## 5   HG00136        GBR   M
## 6   HG00137        GBR   F
# simulate some phenotype data
# set.seed(4)
# annot$outcome <- rnorm(nrow(annot))


annot <- df2
rownames(annot) <- NULL

annot$IID <- as.character(annot$IID)
annot$Sex <- as.character(annot$Sex)
annot$STATUS <- as.character(annot$STATUS)
annot$Population <- as.character(annot$Population)

# # reorder rows
annot <- annot[match(seqGetData(gds, "sample.id"), annot$IID),]

sum(grepl("_R1$|_R2$", annot$IID))




library(SNPRelate)
# read in GDS data
gds <- snpgdsOpen(HapMap_genoData)
snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6, 
                          ld.threshold=sqrt(0.1), verbose=FALSE)
pruned <- unlist(snpset, use.names=FALSE)
length(pruned)







# Reading in the GRM from PC-Relate
# A mixed model for genetic association testing typically includes a genetic
# relationship matrix (GRM) to account for genetic similarity among sample
# individuals. If we are using kinship coefficient estimates from PC-Relate to
# construct this GRM, then the function pcrelateToMatrix should be used to
# provide the matrix in the appropriate format for fitNullModel.








# mypcrel contains Kinship Estimates from a previous PC-Relate analysis
myGRM <- pcrelateToMatrix(mypcrel)
myGRM[1:5,1:5]
Note that both the row and column names of this matrix are the same scanIDs as used in the scan annotation data.

Mixed Model Association Testing
There are two steps to performing genetic association testing with r Biocpkg("GENESIS"). First, the null model (i.e. the model with no SNP genotype term) is fit using the fitNullModel function. Second, the output of the null model fit is used in conjunction with the genotype data to quickly run SNP-phenotype association tests using the assocTestSingle function. There is a computational advantage to splitting these two steps into two function calls; the null model only needs to be fit once, and SNP association tests can be paralelized by chromosome or some other partitioning to speed up analyses (details below).

Fit the Null Model
The first step for association testing with r Biocpkg("GENESIS") is to fit the mixed model under the null hypothesis that each SNP has no effect. This null model contains all of the covariates, including ancestry representative PCs, as well as any random effects, such as a polygenic effect due to genetic relatedness, but it does not include any SNP genotype terms as fixed effects.

Using the fitNullModel function, random effects in the null model are specified via their covariance structures. This allows for the inclusion of a polygenic random effect using a kinship matrix or genetic relationship matrix (GRM).

Quantitative Phenotypes
A linear mixed model (LMM) should be fit when analyzing a quantitative phenotype. The example R code below fits a basic null mixed model.

# fit the null mixed model
nullmod <- fitNullModel(scanAnnot, outcome = "pheno", covars = "pc1", 
                        cov.mat = myGRM, family = gaussian)
the first argument is the class ScanAnnotationDataFrame or data.frame object containing the sample data
outcome specifies the name of the outcome variable in scanAnnot
covars specifies the names of the covariates in scanAnnot
cov.mat specifies the covariance structures for the random effects included in the model
family should be gaussian for a quantitative phenotype, specifying a linear mixed model
The Average Information REML (AIREML) procedure is used to estimate the variance components of the random effects. When verbose = TRUE, the variance component estimates, the log-likelihood, and the residual sum of squares in each iteration are printed to the R console (shown above). In this example, Sigma^2_A is the variance component for the random effect specified in cov.mat, and Sigma^2_E is the residual variance component.

Multiple Fixed Effect Covariates
The model can be fit with multiple fixed effect covariates by setting covars equal to vector of covariate names. For example, if we wanted to include the variables "pc1", "pc2", "sex", and "age" all as covariates in the model:
  
  nullmod <- fitNullModel(scanAnnot, outcome = "pheno", 
                          covars = c("pc1","pc2","sex","age"), 
                          cov.mat = myGRM, family = gaussian)
Multiple Random Effects
The model also can be fit with multiple random effects. This is done by setting cov.mat equal to a list of matrices. For example, if we wanted to include a polygenic random effect with covariance structure given by the matrix "myGRM" and a household random effect with covariance structure specified by the matrix "H":
  
  nullmod <- fitNullModel(scanAnnot, outcome = "pheno", covars = "pc1",
                          cov.mat = list("GRM" = myGRM, "House" = H), 
                          family = gaussian)
The names of the matrices in cov.mat determine the names of the variance component parameters. Therefore, in this example, the output printed to the R console will include Sigma^2_GRM for the random effect specified by "myGRM", Sigma^2_House for the random effect specified by "H", and Sigma^2_E for the residual variance component.

Note: the row and column names of each matrix used to specify the covariance structure of a random effect in the mixed model must be the unique scanIDs for each sample in the analysis.

Heterogeneous Residual Variances
LMMs are typically fit under an assumption of constant (homogeneous) residual variance for all observations. However, for some outcomes, there may be evidence that different groups of observations have different residual variances, in which case the assumption of homoscedasticity is violated. group.var can be used in order to fit separate (heterogeneous) residual variance components by some grouping variable. For example, if we have a categorical variable "study" in our scanAnnot, then we can estimate a different residual variance component for each unique value of "study" by using the following code:
  
  nullmod <- fitNullModel(scanAnnot, outcome = "pheno", covars = "pc1", 
                          cov.mat = myGRM, family = gaussian, 
                          group.var = "study")
In this example, the residual variance component Sigma^2_E is replaced with group specific residual variance components Sigma^2_study1, Sigma^2_study2, ..., where "study1", "study2", ... are the unique values of the "study" variable.

Binary Phentoypes
Ideally, a generalized linear mixed model (GLMM) would be fit for a binary phenotype; however, fitting a GLMM is much more computationally demanding than fitting an LMM. To provide a compuationally efficient approach to fitting such a model, fitNullModel uses the penalized quasi-likelihood (PQL) approximation to the GLMM (Breslow and Clayton). The implementation of this procedure in r Biocpkg("GENESIS") is the same as in GMMAT (Chen et al.), and more details can be found in that manuscript. If our outcome variable, "pheno", were binary, then the same R code could be used to fit the null model, but with family = binomial.

nullmod <- fitNullModel(scanAnnot, outcome = "pheno", covars = "pc1", 
                        cov.mat = myGRM, family = binomial)
Multiple fixed effect covariates and multiple random effects can be specified for binary phenotypes in the same way as they are for quantitative phenotypes. group.var does not apply here.

Run SNP-Phenotype Association Tests
The second step for association testing with r Biocpkg("GENESIS") is to use the fitted null model to test the SNPs in the GenotypeData object for association with the specified outcome variable. This is done with the assocTestSingle function. The use of assocTestSingle for running association tests with a quantitative or binary phenotype is identical.

Before we can run an association test on a GenotypeData object, we much first decide how many SNPs we want to read at a time. We do this by creating a GenotypeBlockIterator object that defines blocks of SNPs. The default setting is to read 10,000 SNPs in each block, but this may be changed with the snpBlock argument.

genoIterator <- GenotypeBlockIterator(HapMap_genoData, snpBlock=5000)
The example R code below runs the association analyses using the null model we fit using fitNullModel in the previous section.

assoc <- assocTestSingle(genoIterator, null.model = nullmod)
genoData is a GenotypeData class object
null.model is the output from fitNullModel
By default, the function will perform association tests at all SNPs in the genoData object. However, for computational reasons it may be practical to parallelize this step, partitioning SNPs by chromosome or some other pre-selected grouping. If we only want to test a pre-specified set of SNPs, this can be done by passing a vector of snpID values to the snpInclude argument when we create the iterator.

# mysnps is a vector of snpID values for the SNPs we want to test
genoIterator <- GenotypeBlockIterator(HapMap_genoData, snpInclude=mysnps)
assoc <- assocTestSingle(genoIterator, null.model = nullmod)
Output
The Null Model
The fitNullModel function will return a list with a large amount of data. Some of the more useful output for the user includes:
  
  varComp: the variance component estimates for the random effects
fixef: a data.frame with point estimates, standard errors, test statistics, and p-values for each of the fixed effect covariates
fitted.values: the fitted values from the model
resid.marginal and resid.conditional: the marginal and conditional residuals from the model
There are also metrics assessing model fit such as the log-likelihood (logLik), restricted log-likelihood (logLikR), and the Akaike information criterion (AIC). Additionally, there are some objects such as the working outcome vector (workingY) and the Cholesky decomposition of the inverse of the estimated phenotype covariance matrix (cholSigmaInv) that are used by the assocTestSingle function for association testing. Further details describing all of the output can be found with the command help(fitNullModel).

The Association Tests
The assocTestSingle function will return a data.frame with summary information from the association test for each SNP. Each row corresponds to a different SNP.

head(assoc)
variant.id: the unique snp ID
chr: the chromosome
pos: the position
n.obs: the number of samples analyzed at that SNP
freq: the frequency of the tested ("A") allele
MAC: the minor allele count
Score: the value of the score function
Score.SE: the estimated standard error of the score
Score.Stat: the score Z test statistic
Score.pval: the p-value based on the score test statistic
Est: an approximation of the effect size estimate (beta) for that SNP
Est.SE: an approximation of the standard error of the effect size estimate
PVE: an approximation of the proportion of phenotype variance explained
Further details describing all of the output can be found with the command help(assocTestSingle).

Heritability Estimation
It is often of interest to estimate the proportion of the total phenotype variability explained by the entire set of genotyped SNPs avaialable; this provides an estimate of the narrow sense heritability of the trait. One method for estimating heritability is to use the variance component estimates from the null mixed model. r Biocpkg("GENESIS") includes the varCompCI function for computing the proportion of variance explained by each random effect along with 95% confidence intervals.

varCompCI(nullmod, prop = TRUE)
close(genoIterator)
the first argument is the output from fitNullModel
prop is a logical indicator of whether the point estimates and confidence intervals should be returned as the proportion of total variability explained (TRUE) or on the orginal scale (FALSE)
When additional random effects are included in the model (e.g. a shared household effect), varCompCI will also return the proportion of variability explained by each of these components.

Note: varCompCI can not compute proportions of variance explained when heterogeneous residual variances are used in the null model (i.e. group.var is used in fitNullModel). Confidence intervals can still be computed for the variance component estimates on the original scale by setting prop = FALSE.

Note: variance component estimates are not interpretable for binary phenotypes when fit using the PQL method implemented in fitNullModel; proportions of variance explained should not be calculated for these models.



##################################################################################











# # Overview
# # This vignette provides a description of how to use the GENESIS package to
# # analyze sequence data. We demonstrate the use of mixed models for genetic
# # association testing, as PC-AiR PCs can be used as fixed effect covariates to
# # adjust for population stratification, and a kinship matrix (or genetic
# # relationship matrix) estimated from PC-Relate can be used to account for
# # phenotype correlation due to genetic similarity among samples. To illustrate
# # the methods, we use a small subset of data from 1000 Genomes Phase 3.
# 
# # 2Convert VCF to GDS
# # The first step is to convert a VCF file into the GDS file format used by
# # GENESIS. We use the SeqArray package, which defines the extended GDS format
# # used to capture all data in a VCF file. If the VCF files are split by
# # chromosome, they can be combined into a single GDS file.
# 
library(SeqArray)
# vcffile <- system.file("extdata", "1KG",
#                        paste0("1KG_phase3_subset_chr", 1:22, ".vcf.gz"),
#                        package="GENESIS")
# gdsfile <- tempfile()
# seqVCF2GDS(vcffile, gdsfile, verbose=FALSE)
# gds <- seqOpen(gdsfile)
# gds


# BiocManager::install("GENESIS")
library(GENESIS)
library(SeqArray)
# vcffile <- system.file("extdata", "1KG", 
                       # paste0("1KG_phase3_subset_chr", 1:22, ".vcf.gz"), 
                       # package="GENESIS")

ln -s /40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files2/FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic.vcf.gz
ln -s /40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files2/FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic.vcf.gz.tbi

# vcffile <- ("FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic.vcf.gz")

vcffile <- c("../FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic_with_no_chr.vcf.gz")
# gdsfile <- tempfile()
## # convert, save in "tmp.gds" with the default lzma compression algorithm
seqVCF2GDS(vcffile, "/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files4/GENESIS/gdsfile", verbose=TRUE, parallel = 10L)
# seqVCF2GDS(vcffile,"/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files3/GENESIS/gdsfile", storage.option="ZIP_RA", parallel = 10L)
gds <- seqOpen("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files4/GENESIS/gdsfile")
gds


# Create a SeqVarData object
# Next, we combine the GDS file with information about the samples, which we
# store in an AnnotatedDataFrame (defined in the Biobase package). An
# AnnotatedDataFrame combines a data.frame with metadata describing each column.
# A SeqVarData object (defined in the SeqVarTools package), contains both an
# open GDS file and an AnnotatedDataFrame describing the samples. The sample.id
# column in the AnnotatedDataFrame must match the sample.id node in the GDS
# file.

library(GENESIS)
library(Biobase)
library(SeqVarTools)


df <- read.table("/40/AD/AD_Seq_Data/03.-phenotype/2021-01-Aquilla-phenotype/Joint_Call_Round1/Aquilla_9576_phenotype.csv", sep = ",", header = TRUE)
head(df)
df$Aquilla.vcfID <- as.character(df$Aquilla.vcfID)

df1 <- read.table("/100/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files4/FASe_3894_from_AQUILLA_WXS_SNPS_INDELS_picard_biallelic.fam", sep = " ", header = FALSE)
dim(df1)
head(df1)
colnames(df1) <- c("FID", "IID", "MID", "PID", "SEX", "STATUS")
df1$IID <- as.character(df1$IID)

sum(df1$IID %in% df$Aquilla.vcfID)
# sum(grepl("_R1$", df1$IID))

df2 <- cbind(df1,df[match(df1$IID, df$Aquilla.vcfID), c("Aquilla.vcfID", "Sex", "STATUS")])


df2_nonNA <- df2[!is.na(df2$Aquilla.vcfID),]
df2_R1R2 <- df2[grepl ("_R1$|_R2$", as.character(df2$IID)),]


df2_R1R2 <- cbind(df2_R1R2[1:6], df[match(gsub("_R1|_R2","", df2_R1R2$IID), df$Aquilla.vcfID), c("Aquilla.vcfID", "Sex", "STATUS")])

head(df2_nonNA)
head(df2_R1R2)

df2 <- rbind(df2_nonNA, df2_R1R2)
df2 <- df2[c(2,8,9)]
df2$Population <- "Empty"
head(df2)
df2 <- df2[c("IID", "Population", "Sex", "STATUS")]


# data(sample_annotation_1KG)
# annot <- sample_annotation_1KG
# head(annot)
##   sample.id Population sex
## 1   HG00110        GBR   F
## 2   HG00116        GBR   M
## 3   HG00120        GBR   F
## 4   HG00128        GBR   F
## 5   HG00136        GBR   M
## 6   HG00137        GBR   F
# simulate some phenotype data
# set.seed(4)
# annot$outcome <- rnorm(nrow(annot))


annot <- df2
rownames(annot) <- NULL

annot$IID <- as.character(annot$IID)
annot$Sex <- as.character(annot$Sex)
annot$STATUS <- as.character(annot$STATUS)
annot$Population <- as.character(annot$Population)

# # reorder rows
annot <- annot[match(seqGetData(gds, "sample.id"), annot$IID),]

sum(grepl("_R1$|_R2$", annot$IID))


colnames(annot)[1] <- "sample.id" 
# metadata <- data.frame(labelDescription=c("sample id", 
#                                           "1000 genomes population", 
#                                           "sex", 
#                                           "simulated phenotype"),
#                        row.names=names(annot))

metadata <- data.frame(labelDescription=c("sample id", 
                                          "1000 genomes population", 
                                          "sex", 
                                          "simulated phenotype"),
                       row.names=names(annot))


annot <- AnnotatedDataFrame(annot, metadata)
sum(is.na(annot$sample.id))

all.equal(annot$sample.id, seqGetData(gds, "sample.id"))
## [1] TRUE
seqData <- SeqVarData(gds, sampleData=annot)


# 3Population structure and relatedness
# PC-AiR and PC-Relate are described in detail in a separate vignette. Here, we
# demonstrate their use to provide adjustment for population structure and
# relatedness in a mixed model.

# 3.1KING
# Step 1 is to get initial estimates of kinship using KING, which is robust to
# population structure but not admixture. The KING algorithm is available in
# SNPRelate. We select a subset of variants for this calculation with LD
# pruning.

library(SNPRelate)

# # LD pruning to get variant set
snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6,
                          ld.threshold=sqrt(0.1), num.thread=10L, verbose=TRUE)

# Working space: 3,894 samples, 1,990,134 SNVs
# using 10 (CPU) cores
# sliding window: 10,000,000 basepairs, Inf SNPs
# |LD| threshold: 0.316228
# method: correlation
# Chromosome 1: 10.43%, 36,245/347,491
# Chromosome 2: 12.34%, 32,235/261,227
# Chromosome 3: 12.48%, 23,461/188,055
# Chromosome 4: 14.82%, 19,589/132,150
# Chromosome 5: 10.46%, 16,666/159,289
# Chromosome 6: 13.70%, 20,452/149,309
# Chromosome 7: 12.19%, 20,762/170,274
# Chromosome 8: 13.80%, 16,629/120,500
# Chromosome 9: 11.52%, 16,732/145,246
# Chromosome 10: 13.00%, 17,823/137,129
# Chromosome 11: 10.18%, 20,475/201,213
# Chromosome 12: 8.59%, 16,396/190,844
# Chromosome 13: 15.70%, 9,340/59,497
# Chromosome 14: 11.27%, 11,907/105,670
# Chromosome 15: 9.30%, 11,444/123,001
# Chromosome 16: 7.74%, 12,637/163,278
# Chromosome 17: 8.58%, 16,199/188,891
# Chromosome 18: 15.64%, 8,389/53,622
# Chromosome 19: 5.21%, 11,747/225,294
# Chromosome 20: 10.14%, 9,090/89,669
# Chromosome 21: 11.71%, 4,432/37,856
# Chromosome 22: 8.12%, 6,824/84,045
# 359,474 markers are selected in total.
# Warning message:
#   In snpgdsLDpruning(gds, method = "corr", slide.max.bp = 1e+07, ld.threshold = sqrt(0.1),  :
#                        The current version of 'snpgdsLDpruning()' does not support multi-threading.
# > pruned <- unlist(snpset, use.names=FALSE)
# > length(pruned)
# [1] 359474


# snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=50000, 
#                           ld.threshold=sqrt(0.1), num.thread=10L, verbose=TRUE)
# 
# 
# # > library(SeqVarTools)
# # > library(SNPRelate)
# # SNPRelate -- supported by Streaming SIMD Extensions 2 (SSE2)
# # > snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=5+                           ld.threshold=sqrt(0.1), num.threadSNV pruning based on LD:
# #                               Excluding 79,971 SNVs on non-autosomes
# #                             Calculating allele counts/frequencies ...
# #                             [..................................................]  0%, ETC:[==================================================] 100%, completed in 2s
# #                             Excluding 1,343,416 SNVs (monomorphic: TRUE, MAF: NaN, missing rate: NaN)
# #                             Working space: 3,894 samples, 1,990,134 SNVs
# #                             using 10 (CPU) cores
# #                             sliding window: 50,000 basepairs, Inf SNPs
# #                             |LD| threshold: 0.316228
# #                             method: correlation
# #                             Chromosome 1: 35.99%, 125,064/347,491
# #                             Chromosome 2: 36.35%, 94,964/261,227
# #                             Chromosome 3: 37.33%, 70,208/188,055
# #                             Chromosome 4: 37.78%, 49,923/132,150
# #                             Chromosome 5: 35.17%, 56,021/159,289
# #                             Chromosome 6: 38.12%, 56,916/149,309
# #                             Chromosome 7: 36.07%, 61,426/170,274
# #                             Chromosome 8: 36.34%, 43,784/120,500
# #                             Chromosome 9: 36.82%, 53,481/145,246
# #                             Chromosome 10: 37.80%, 51,829/137,129
# #                             Chromosome 11: 35.60%, 71,641/201,213
# #                             Chromosome 12: 34.73%, 66,284/190,844
# #                             Chromosome 13: 38.24%, 22,754/59,497
# #                             Chromosome 14: 35.39%, 37,397/105,670
# #                             Chromosome 15: 32.53%, 40,009/123,001
# #                             Chromosome 16: 33.19%, 54,194/163,278
# #                             Chromosome 17: 35.21%, 66,515/188,891
# #                             Chromosome 18: 39.22%, 21,032/53,622
# #                             Chromosome 19: 33.42%, 75,290/225,294
# #                             Chromosome 20: 35.73%, 32,040/89,669
# #                             Chromosome 21: 34.23%, 12,957/37,856
# #                             Chromosome 22: 35.12%, 29,515/84,045
# #                             1,193,244 markers are selected in total.
# #                             Warning message:
# #                               In snpgdsLDpruning(gds, method = "corr", slide.max.bp = 50000, ld.threshold = sqrt(0.1),  :
# #                                                    The current version of 'snpgdsLDpruning()' does not support multi-threading.
# #                                                  > pruned <- unlist(snpset, use.names=FALSE)
# #                                                  > length(pruned)
# #                                                  [1] 1193244
                                                  


pruned <- unlist(snpset, use.names=FALSE)

king <- snpgdsIBDKING(gds, snp.id=pruned, num.thread=10L, verbose=TRUE)
kingMat <- king$kinship
dimnames(kingMat) <- list(king$sample.id, king$sample.id)


# 3.2PC-AiR
# The next step is PC-AiR, in which we select a set of unrelated samples that is
# maximally informative about all ancestries in the sample, use this unrelated
# set for Principal Component Analysis (PCA), then project the relatives onto
# the PCs. We use a kinship threshold of degree 3 (unrelated is less than first
# cousins). In this example, we use the KING estimates for both kinship (kinobj)
# and ancestry divergence (divobj). KING kinship estimates are negative for
# samples with different ancestry.

pcs <- pcair(seqData, 
             kinobj=kingMat, kin.thresh=2^(-9/2),
             divobj=kingMat, div.thresh=-2^(-9/2),
             snp.include=pruned, num.cores = 10L, verbose = TRUE)


pcs <- pcair(gds,
             kinobj=kingMat, kin.thresh=2^(-9/2),
             divobj=kingMat, div.thresh=-2^(-9/2),
             snp.include=pruned, num.cores = 10L, verbose = TRUE)



# > pcs <- pcair(gds,
#                +              kinobj=kingMat, kin.thresh=2^(-9/2),
#                +              divobj=kingMat, div.thresh=-2^(-9/2),
#                +              snp.include=pruned)
# Identifying relatives for each sample using kinship threshold 0.0441941738241592
# Identifying pairs of divergent samples using divergence threshold -0.0441941738241592
# Partitioning samples into unrelated and related sets...
# ...1000 samples added to related.set...
# Unrelated Set: 2833 Samples
# Related Set: 1062 Samples
# Performing PCA on the Unrelated Set...
# Principal Component Analysis (PCA) on genotypes:
#   Calculating allele counts/frequencies ...
# [==================================================] 100%, completed in 24s
# Excluding 55,065 SNVs (monomorphic: TRUE, MAF: NaN, missing rate: NaN)
# Working space: 2,341 samples, 312,639 SNVs
# using 1 (CPU) core
# CPU capabilities: Double-Precision SSE2
# Wed Apr 21 14:57:42 2021    (internal increment: 1940)
# [==================================================] 100%, completed in 7.6m
# Wed Apr 21 15:05:18 2021    Begin (eigenvalues and eigenvectors)
# Wed Apr 21 15:05:21 2021    Done.
# Predicting PC Values for the Related Set...
# SNP loading:
#   Working space: 2341 samples, 312639 SNPs
# using 1 (CPU) core
# using the top 32 eigenvectors
# Wed Apr 21 15:05:21 2021    (internal increment: 15548)
# [==================================================] 100%, completed in 51s
# Wed Apr 21 15:06:12 2021    Done.
# Error in .InitFile(gdsobj, sample.id = sample.id, snp.id = loadobj$snp.id) :
#   Some of sample.id do not exist!
#   


# We need to determine which PCs are ancestry informative. To do this we need
# population information for the 1000 Genomes samples. We make a parallel
# coordinates plot, color-coding by 1000 Genomes population.

library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(GGally)

pc.df <- as.data.frame(pcs$vectors)
names(pc.df) <- paste0("PC", 1:ncol(pcs$vectors))
pc.df$sample.id <- row.names(pcs$vectors)
pc.df <- left_join(pc.df, pData(annot), by="sample.id")

pop.cols <- setNames(brewer.pal(12, "Paired"),
                     c("ACB", "ASW", "CEU", "GBR", "CHB", "JPT", 
                       "CLM", "MXL", "LWK", "YRI", "GIH", "PUR"))

ggplot(pc.df, aes(PC1, PC2, color=Population)) + geom_point() +
  scale_color_manual(values=pop.cols)


ggparcoord(pc.df, columns=1:10, groupColumn="Population", scale="uniminmax") +
  scale_color_manual(values=pop.cols) +
  xlab("PC") + ylab("")



# PC-Relate
# The first 2 PCs separate populations, so we use them to compute kinship
# estimates adjusting for ancestry. The pcrelate method requires creating a
# SeqVarBlockIterator object, which should iterate over the pruned SNPs only.

seqSetFilter(seqData, variant.id=pruned)
## # of selected variants: 10,823
iterator <- SeqVarBlockIterator(seqData, variantBlock=20000, verbose=FALSE)
pcrel <- pcrelate(iterator, pcs=pcs$vectors[,1:2], training.set=pcs$unrels)
seqResetFilter(seqData, verbose=FALSE)
kinship <- pcrel$kinBtwn

ggplot(kinship, aes(k0, kin)) +
  geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
  geom_point(alpha=0.5) +
  ylab("kinship estimate") +
  theme_bw()


# To improve our estimates for PCs and kinship, we could run another iteration
# of PC-AiR and PC-Relate, this time using the PC-Relate kinship estimates as
# the kinobj argument to pcair. The KING matrix is still used for ancestry
# divergence. We could then use those new PCs to calculate revised kinship
# estimates.

# 4Association tests
# 4.1Null model
# The first step for association testing is to fit the model under the null
# hypothesis that each SNP has no effect. This null model contains all of the
# covariates, including ancestry representative PCs, as well as any random
# effects, such as a polygenic effect due to genetic relatedness, but it does
# not include any SNP genotype terms as fixed effects.

# The type of model fit depends on the arguments to fitNullModel. Including a
# cov.mat argument will result in a mixed model, while omitting this argument
# will run a standard linear model. A logistic model is specified with
# family="binomial". In the case of a logistic model and a covariance matrix,
# fitNullModel will use the GMMAT algorithm. Including a group.var argument will
# allow heteroscedastic variance (for linear models or linear mixed models
# only).

# add PCs to sample annotation in SeqVarData object
annot <- AnnotatedDataFrame(pc.df)
sampleData(seqData) <- annot

# covariance matrix from pcrelate output
grm <- pcrelateToMatrix(pcrel, scaleKin=2)

# fit the null model
nullmod <- fitNullModel(seqData, outcome="outcome", 
                        covars=c("sex", "Population", paste0("PC", 1:2)),
                        cov.mat=grm, verbose=FALSE)
# 4.2Single variant tests
# To run a test using the null model, we first create an iterator object
# specifying how we want variants to be selected. (See the documentation for the
# SeqVarIterator class in SeqVarTools for more details.) For single-variant
# tests (GWAS), it is common to use a block iterator that reads variants in
# blocks (default is 10,000 variants per block).

# For example purposes, we restrict our analysis to chromosome 1. The
# seqSetFilter function can be used to restrict the set of variants tested in
# other ways (e.g., variants that pass a quality filter).

# select chromosome 1
seqSetFilterChrom(seqData, include=1)
## # of selected variants: 1,120
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
assoc <- assocTestSingle(iterator, nullmod, verbose=FALSE)
head(assoc)
##   variant.id chr     pos allele.index n.obs  freq MAC       Score Score.SE
## 1          1   1  828740            1   100 0.035   7   0.4554795 2.888007
## 2          2   1  913272            1   100 0.010   2  -0.2775530 1.403110
## 3          3   1 1171878            1   100 0.005   1   0.2227947 1.015806
## 4          4   1 1242288            1   100 0.025   5  -1.2206682 2.163235
## 5          5   1 1378837            1   100 0.670  66 -11.2824035 6.080219
## 6          6   1 1403820            1   100 0.015   3   2.3221454 1.692363
##   Score.Stat Score.pval         Est    Est.SE          PVE
## 1  0.1577141  0.8746821  0.05461001 0.3462595 0.0002926286
## 2 -0.1978127  0.8431916 -0.14098160 0.7127024 0.0004603457
## 3  0.2193280  0.8263945  0.21591530 0.9844400 0.0005659318
## 4 -0.5642791  0.5725642 -0.26084970 0.4622707 0.0037459651
## 5 -1.8555915  0.0635118 -0.30518497 0.1644678 0.0405079734
## 6  1.3721322  0.1700223  0.81077907 0.5908899 0.0221496917

# The default test is a Score test, but the Wald test is also available for
# continuous outcomes.

# If there are multiallelic variants, each alternate allele is tested
# separately. The allele.index column in the output differentiates between
# different alternate alleles for the same variant.

# We make a QQ plot to examine the results.

qqPlot <- function(pval) {
  pval <- pval[!is.na(pval)]
  n <- length(pval)
  x <- 1:n
  dat <- data.frame(obs=sort(pval),
                    exp=x/n,
                    upper=qbeta(0.025, x, rev(x)),
                    lower=qbeta(0.975, x, rev(x)))
  
  ggplot(dat, aes(-log10(exp), -log10(obs))) +
    geom_line(aes(-log10(exp), -log10(upper)), color="gray") +
    geom_line(aes(-log10(exp), -log10(lower)), color="gray") +
    geom_point() +
    geom_abline(intercept=0, slope=1, color="red") +
    xlab(expression(paste(-log[10], "(expected P)"))) +
    ylab(expression(paste(-log[10], "(observed P)"))) +
    theme_bw()
}    

qqPlot(assoc$Score.pval)

library(qqman)
LOGISTIC <- assoc
colnames(LOGISTIC)[1] <- "SNP"
colnames(LOGISTIC)[2] <- "CHR"
colnames(LOGISTIC)[3] <- "BP"
colnames(LOGISTIC)[10] <- "P"

LOGISTIC$CHR <- as.numeric(LOGISTIC$CHR)

manhattan(LOGISTIC, main = "", ylim=c(0,15), col = c("blue4", "orange3"), suggestiveline = -log10(1e-06), genomewideline = -log10(1e-08), annotateTop = TRUE, annotatePval = -log10(1e-06), chrlabs = as.character(1)) 


# 4.3Aggregate tests
# We can aggregate rare variants for association testing to decrease multiple
# testing burden and increase statistical power. We can create functionally
# agnostic units using a SeqVarWindowIterator. This iterator type generates a
# sliding window over the genome, with user-specified width and step size. We
# can also create units with specific start and end points or containing
# specific variants, using a SeqVarRangeIterator or a SeqVarListIterator.

# In this example, we illustrate defining ranges based on known genes. We run a
# burden test, setting a maximum alternate allele frequency to exclude common
# variants.

library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# return the variants on chromosome 1 as a GRanges object
seqSetFilterChrom(seqData, include=1)
## # of selected variants: 1,120
gr <- granges(gds)

# find variants that overlap with each gene
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gr <- renameSeqlevels(gr, paste0("chr", seqlevels(gr)))
ts <- transcriptsByOverlaps(txdb, gr, columns="GENEID")
# simplistic example - define genes as overlapping transcripts
genes <- reduce(ts)
genes <- renameSeqlevels(genes, sub("chr", "", seqlevels(genes)))

# create an iterator where each successive unit is a different gene
iterator <- SeqVarRangeIterator(seqData, variantRanges=genes, verbose=FALSE)

# do a burden test on the rare variants in each gene
assoc <- assocTestAggregate(iterator, nullmod, AF.max=0.05, test="Burden",
                            verbose=FALSE)
# The output of an aggregate test is a list with two elements: 1) a data.frame
# with the test results for each aggregate unit, and 2) a list of data.frames
# containing the variants in each aggregate unit.

head(assoc$results)
##   n.site n.alt n.sample.alt      Score Score.SE Score.Stat Score.pval
## 1      1     3            3  2.3221454 1.692363  1.3721322 0.17002228
## 2      1     5            5  4.0027945 2.246167  1.7820556 0.07474017
## 3      1     1            1 -0.6661171 1.043845 -0.6381379 0.52338391
## 4      1     9            9 -1.8242825 2.581275 -0.7067371 0.47972988
## 5      0     0            0         NA       NA         NA         NA
## 6      1     1            1 -0.4031965 1.023646 -0.3938825 0.69366776
##          Est    Est.SE         PVE
## 1  0.8107791 0.5908899 0.022149692
## 2  0.7933763 0.4452029 0.037360979
## 3 -0.6113339 0.9579965 0.004790765
## 4 -0.2737938 0.3874055 0.005876132
## 5         NA        NA          NA
## 6 -0.3847838 0.9768997 0.001825195
head(assoc$variantInfo)
## [[1]]
##   variant.id chr     pos allele.index n.obs  freq MAC weight
## 1          6   1 1403820            1   100 0.015   3      1
## 
## [[2]]
##   variant.id chr     pos allele.index n.obs  freq MAC weight
## 1          7   1 1421285            1   100 0.025   5      1
## 
## [[3]]
##   variant.id chr     pos allele.index n.obs  freq MAC weight
## 1         12   1 2023475            1   100 0.005   1      1
## 
## [[4]]
##   variant.id chr     pos allele.index n.obs  freq MAC weight
## 1         21   1 3254100            1   100 0.045   9      1
## 
## [[5]]
## [1] variant.id   chr          pos          allele.index n.obs       
## [6] freq         MAC          weight      
## <0 rows> (or 0-length row.names)
## 
## [[6]]
##   variant.id chr     pos allele.index n.obs  freq MAC weight
## 1         24   1 3818550            1   100 0.005   1      1