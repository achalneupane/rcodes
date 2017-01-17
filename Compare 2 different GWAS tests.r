rm(list=ls())
options(stringsAsFactors=FALSE)
options(width=400)
list.files(path=getwd())

read.in <- commandArgs(trailingOnly=TRUE)

snp   <- read.in[1] #this is the SNP id
snpid <- read.in[1] #this is an alternate it the snp id that sometimes I write in my scripts
chr   <- read.in[2] #chromosome

raw.data <- "/lustre/home/tprobins/uveitis/dataVersion5/ImputingImmunochipDec2013/postQC_Imputed_uveitis_chromosomes"  #data folder where files are

wf <- "/lustre/home/tprobins/uveitis/dataVersion5/OR_comparison/WriteFolder" #write folder
gtl <- "/lustre/home/tprobins/uveitis/dataVersion5/OR_comparison/Gtool_logs" #log folder
adsam <- "/lustre/home/tprobins/uveitis/dataVersion5/OR_comparison/AdjustedSams" #Sam files with eigenvectors included

exp2.top <- read.table("/lustre/home/tprobins/uveitis/dataVersion5/OR_comparison/supportFiles/SetSet2_results_ANNOVAR_annotated.txt",header=TRUE) #one set of snps I want to check
exp2.top.snps <- exp2.top[,1:2] #this just extracts the snp ids and chromosomes form the file

exp3.top <- read.table("/lustre/home/tprobins/uveitis/dataVersion5/OR_comparison/supportFiles/SetSet3_results_ANNOVAR_annotated.txt",header=TRUE) #another set of snps I want to check
exp3.top.snps <- exp3.top[,1:2] #this just extracts the snp ids and chromosomes form the file

as.snps <- read.table("/lustre/home/tprobins/uveitis/dataVersion5/OR_comparison/supportFiles/AS_snp_list.txt") #another set of snps I want to check
colnames(as.snps) <- c("rsid","Chromosome") #this just extracts the snp ids and chromosomes form the file

top.snps <- rbind(exp2.top.snps,exp3.top.snps,as.snps[,c(2,1)])
colnames(top.snps) <- c("Chromosome","SNPid") #top snps

#the above files obviously aren't necessary for an 'args' type script but it wasn't always an 'args' type script

pheno2 <- read.table("/lustre/home/tprobins/uveitis/dataVersion5/OR_comparison/supportFiles/Set2.August.UveitisExperiment.plink.pheno") #phenotype file for set 1 (which I call 2): ID1, ID2, and phenotype as columns
pheno3 <- read.table("/lustre/home/tprobins/uveitis/dataVersion5/OR_comparison/supportFiles/Set3.August.UveitisExperiment.plink.pheno")#phenotype file for set 2 (which I call 3): ID1, ID2, and phenotype as columns

people2.method1 <- "/lustre/home/tprobins/uveitis/dataVersion5/OR_comparison/supportFiles/Set2.August.UveitisExperiment.plink.extract" #plink type extract file with 2 columns
people3.method1 <- "/lustre/home/tprobins/uveitis/dataVersion5/OR_comparison/supportFiles/Set3.August.UveitisExperiment.plink.extract" #plink type extract file with 2 columns

people2.method2 <- "/lustre/home/tprobins/uveitis/dataVersion5/OR_comparison/supportFiles/ControlExtractFile1.plink.extract" #plink type extract file with 2 columns
people3.method2 <- "/lustre/home/tprobins/uveitis/dataVersion5/OR_comparison/supportFiles/ControlExtractFile2.plink.extract" #plink type extract file with 2 columns

###########################################################################################

#add HLA-B27 dosage
b <- read.table("/lustre/home/tprobins/uveitis/dataVersion5/OR_comparison/supportFiles/rs116488202_genotype_uveitis_people.txt",header=TRUE) #this is a read in of a covariate (HLA_b27 in this case)
b$C_dose <- (b$CC*2)+b$CT #this just converts it from 3 column genotype dosage for a single number genotype dosage
b <- b[,-c(2,3,4,5)]
head(b)

#Read in eigens
eigens <- read.table("/lustre/home/tprobins/uveitis/dataVersion5/OR_comparison/supportFiles/Labelled_10_Uveitis_PCAs_withoutASregions_excluded.txt") #this is a read in of a covariate (eigens in this case)

#################

cat(snpid,file=paste(wf,"/SNpextractfile_",snp,".txt",sep="")) #put the snp to extract in a file so gtool can find it

#ExtractSNPInformation
######################
#Exp2 - method 1

system.line <- paste("/lustre/home/tprobins/programs/gtool/gtool -S --g ",raw.data,"/Imputed_Chr",chr,".haps --s ",adsam,"/Imputed_Chr",chr,".sample.adj --sample_id ",people2.method1, " --inclusion ",wf,"/SNpextractfile_",snp,".txt --og ",wf,"/m1_SNP_exp2_",snp,".gen --os ",wf,"/m1_SNP_exp2_",snp,".sample --log ",gtl,"/Exp2_method1_",snp,".log",sep="")
system(system.line)


###########################################################################################################################

 #Exp2 - method 1
 #Converting .gen triple dosage format to a file for conditioning
 x <- read.table(paste(wf,"/m1_SNP_exp2_",snp,".gen",sep=""))
 file.start <- x[1,1:5]
 file.main <- x[1,6:(ncol(x))]
 file.matrix <- matrix(file.main,ncol=3,byrow=T)
 vec1 <- as.numeric(file.matrix[,1])
 vec2 <- as.numeric(file.matrix[,2])
 vec3 <- as.numeric(file.matrix[,3])
 filedf <- cbind(vec1,vec2,vec3)
 filedf <- data.frame(filedf)
 colnames(filedf) <- c("II","IA","AA")

 #For conditioning
 filedf$cond <- (filedf$II*2)+filedf$IA
 sam <- read.table(paste(wf,"/m1_SNP_exp2_",snp,".sample",sep=""),header=TRUE)
 sam <- sam[-1,]
 outfile <- paste(wf,"/m1_Keep_safe_exp2_",snp,".txt",sep="")
 write.table(cbind(sam[,1:2],filedf$cond),paste(outfile),col.names=F,row.names=F,quote=F)
 data1.m1 <- read.table(paste(wf,"/m1_Keep_safe_exp2_",snp,".txt",sep=""))
#############################################################################################################################
#Exp2 - method 2
system.line <- paste("/lustre/home/tprobins/programs/gtool/gtool -S --g ",wf,"/m1_SNP_exp2_",snp,".gen --s ",wf,"/m1_SNP_exp2_",snp,".sample --sample_id ",people2.method2, " --og ",wf,"/m2_SNP_exp2_",snp,".gen --os ",wf,"/m2_SNP_exp2_",snp,".sample --log ",gtl,"/Exp2_method2_",snp,".log",sep="")
system(system.line)

###########################################################################################################################
#Converting .gen triple dosage format to a file for conditioning
 x <- read.table(paste(wf,"/m2_SNP_exp2_",snp,".gen",sep=""))
 file.start <- x[1,1:5]
 file.main <- x[1,6:(ncol(x))]
 file.matrix <- matrix(file.main,ncol=3,byrow=T)
 vec1 <- as.numeric(file.matrix[,1])
 vec2 <- as.numeric(file.matrix[,2])
 vec3 <- as.numeric(file.matrix[,3])
 filedf <- cbind(vec1,vec2,vec3)
 filedf <- data.frame(filedf)
 colnames(filedf) <- c("II","IA","AA")

 #For conditioning
 filedf$cond <- (filedf$II*2)+filedf$IA
 sam <- read.table(paste(wf,"/m2_SNP_exp2_",snp,".sample",sep=""),header=TRUE)
 sam <- sam[-1,]
 outfile <- paste(wf,"/m2_Keep_safe_exp2_",snp,".txt",sep="")
 write.table(cbind(sam[,1:2],filedf$cond),paste(outfile),col.names=F,row.names=F,quote=F)
 data1.m2 <- read.table(paste(wf,"/m2_Keep_safe_exp2_",snp,".txt",sep=""))


############################################################################################################################################################
#Exp3 - method 1

system.line <- paste("/lustre/home/tprobins/programs/gtool/gtool -S --g ",raw.data,"/Imputed_Chr",chr,".haps --s ",adsam,"/Imputed_Chr",chr,".sample.adj --inclusion ",wf,"/SNpextractfile_",snp,".txt --sample_id ",people3.method1, " --og ",wf,"/m1_SNP_exp3_",snp,".gen --os ",wf,"/m1_SNP_exp3_",snp,".sample --log ",gtl,"/Exp3_method1_",snp,".log",sep="")
system(system.line)


###########################################################################################################################
#Converting .gen triple dosage format to a file for conditioning
 x <- read.table(paste(wf,"/m1_SNP_exp3_",snp,".gen",sep=""))
 file.start <- x[1,1:5]
 file.main <- x[1,6:(ncol(x))]
 file.matrix <- matrix(file.main,ncol=3,byrow=T)
 vec1 <- as.numeric(file.matrix[,1])
 vec2 <- as.numeric(file.matrix[,2])
 vec3 <- as.numeric(file.matrix[,3])
 filedf <- cbind(vec1,vec2,vec3)
 filedf <- data.frame(filedf)
 colnames(filedf) <- c("II","IA","AA")

 #For conditioning
 filedf$cond <- (filedf$II*2)+filedf$IA
 sam <- read.table(paste(wf,"/m1_SNP_exp3_",snp,".sample",sep=""),header=TRUE)
 sam <- sam[-1,]
 outfile <- paste(wf,"/m1_Keep_safe_exp3_",snp,".txt",sep="")
 write.table(cbind(sam[,1:2],filedf$cond),paste(outfile),col.names=F,row.names=F,quote=F)
 data2.m1 <- read.table(paste(wf,"/m1_Keep_safe_exp3_",snp,".txt",sep=""))

#########################################################################################################################
#Exp3 - method 2

system.line <- paste("/lustre/home/tprobins/programs/gtool/gtool -S --g ",wf,"/m1_SNP_exp3_",snp,".gen --s ",wf,"/m1_SNP_exp3_",snp,".sample --sample_id ",people3.method2, " --og ",wf,"/m2_SNP_exp3_",snp,".gen --os ",wf,"/m2_SNP_exp3_",snp,".sample --log ",gtl,"/Exp3_method2_",snp,".log",sep="")
system(system.line)

###########################################################################################################################
#Converting .gen triple dosage format to a file for conditioning
 x <- read.table(paste(wf,"/m2_SNP_exp3_",snp,".gen",sep=""))
 file.start <- x[1,1:5]
 file.main <- x[1,6:(ncol(x))]
 file.matrix <- matrix(file.main,ncol=3,byrow=T)
 vec1 <- as.numeric(file.matrix[,1])
 vec2 <- as.numeric(file.matrix[,2])
 vec3 <- as.numeric(file.matrix[,3])
 filedf <- cbind(vec1,vec2,vec3)
 filedf <- data.frame(filedf)
 colnames(filedf) <- c("II","IA","AA")

 #For conditioning
 filedf$cond <- (filedf$II*2)+filedf$IA
 sam <- read.table(paste(wf,"/m2_SNP_exp3_",snp,".sample",sep=""),header=TRUE)
 sam <- sam[-1,]
 outfile <- paste(wf,"/m2_Keep_safe_exp3_",snp,".txt",sep="")
 write.table(cbind(sam[,1:2],filedf$cond),paste(outfile),col.names=F,row.names=F,quote=F)
 data2.m2 <- read.table(paste(wf,"/m2_Keep_safe_exp3_",snp,".txt",sep=""))

#############################################################################################################################################################################################
###########################################################################################################################################################################################

################################################################
#add HLA-B27 dosage
data1.m1.b27 <- merge(data1.m1,b,by.x=1,by.y="ID_1")
data1.m2.b27 <- merge(data1.m2,b,by.x=1,by.y="ID_1")

data2.m1.b27 <- merge(data2.m1,b,by.x=1,by.y="ID_1")
data2.m2.b27 <- merge(data2.m2,b,by.x=1,by.y="ID_1")

###############################################################
#Merge Eigenvectors
d1.m1 <- merge(data1.m1.b27,eigens,by="V1")
d1.m2 <- merge(data1.m2.b27,eigens,by="V1")

d2.m1 <- merge(data2.m1.b27,eigens,by="V1")
d2.m2 <- merge(data2.m2.b27,eigens,by="V1")

d1.m1 <- d1.m1[,-5]
d1.m2 <- d1.m2[,-5]
d2.m1 <- d2.m1[,-5]
d2.m2 <- d2.m2[,-5]

colnames(d1.m1) <- c("ID1","ID2","snp","B27dose","PCA1","PCA2","PCA3","PCA4","PCA5","PCA6","PCA7","PCA8","PCA9","PCA10")
colnames(d1.m2) <- c("ID1","ID2","snp","B27dose","PCA1","PCA2","PCA3","PCA4","PCA5","PCA6","PCA7","PCA8","PCA9","PCA10")
colnames(d2.m1) <- c("ID1","ID2","snp","B27dose","PCA1","PCA2","PCA3","PCA4","PCA5","PCA6","PCA7","PCA8","PCA9","PCA10")
colnames(d2.m2) <- c("ID1","ID2","snp","B27dose","PCA1","PCA2","PCA3","PCA4","PCA5","PCA6","PCA7","PCA8","PCA9","PCA10")

################################################################
#Add phenotypes

head(pheno2)
nrow(pheno2)
nrow(d1.m1)

d1.m1.p <- merge(d1.m1,pheno2,by.x="ID1",by.y="V1")
d1.m1 <- d1.m1.p[,-15]
d1.m2.p <- merge(d1.m2,pheno2,by.x="ID1",by.y="V1")
d1.m2 <- d1.m2.p[,-15]
d2.m1.p <- merge(d2.m1,pheno3,by.x="ID1",by.y="V1")
d2.m1 <- d2.m1.p[,-15]
d2.m2.p <- merge(d2.m2,pheno3,by.x="ID1",by.y="V1")
d2.m2 <- d2.m2.p[,-15]

colnames(d1.m1)[15] <- "Pheno"
colnames(d1.m2)[15] <- "Pheno"
colnames(d2.m1)[15] <- "Pheno"
colnames(d2.m2)[15] <- "Pheno"

#Make logistic regression phenotypes
d1.m1$Pheno <- d1.m1$Pheno -1
d1.m2$Pheno <- d1.m2$Pheno -1
d2.m1$Pheno <- d2.m1$Pheno -1
d2.m2$Pheno <- d2.m2$Pheno -1

###############################################################
#Model 1 Pheno ~ SNP + Eigens
 model1.1 <- glm(Pheno ~ snp + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10, data=d1.m1) #Set 1
 model1.2 <- glm(Pheno ~ snp + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10, data=d2.m1) #Set 2

#Model 2 Pheno ~ SNP + B27 Eigens
 model2.1 <- glm(Pheno ~ snp + B27dose + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7+ PCA8 + PCA9 + PCA10, data=d1.m1) #Set 1
 model2.2 <- glm(Pheno ~ snp + B27dose + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7+ PCA8 + PCA9 + PCA10, data=d2.m1) #Set 2


#############################################################################################################################################               
#Comparison Method 

#(1) Stack one dataset on top of the other
#(2) Add a dummy variable "Z" that indicates the first group (Z= 0) or the second group (Z=1)
#(3) Include an interaction term also SNP*Z
#(4) Fit your logistic regression: Uveitis ~ SNP + B27 + Z + Z*SNP
#(5) P value for interaction term should be the test you want

d1.m2$Z <- 0
d2.m2$Z <- 1

method2 <- rbind(d1.m2,d2.m2)

model.alpha <- glm(Pheno ~ snp + Z + snp*Z + B27dose + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10,data = method2)

m2.pval <- coef(summary(model.alpha))["snp:Z","Pr(>|t|)"]

line.info <- c(snp,beta1,beta2,se1,se2,df,cof1,cof2,pval1,pval2,pval,m2.pval)

#write.table(line.info,paste("fInterum_",i,".txt",sep=""),col.names=F,row.names=F,quote=F)

#colnames() <- c("SNP","beta1","beta2","se1","se2","df","Coef1a","Coef1b","Coef2a","Coef2b","Pval1","Pval2","Pval_mthod1","Pval_mthod2")

write.table(line.info,paste("/lustre/home/tprobins/uveitis/dataVersion5/OR_comparison/Results/fInterum_",snp,".txt",sep=""),col.names=F,row.names=F,quote=F)

#colnames(master) <- c("SNP","beta1","beta2","se1","se2","df","Coef1a","Coef1b","Coef2a","Coef2b","Pval1","Pval2","Pval_mthod1","Pval_mthod2")
#write.table(master,"fFinal_OR_comparison_script.txt",col.names=T,row.names=F,quote=F)


