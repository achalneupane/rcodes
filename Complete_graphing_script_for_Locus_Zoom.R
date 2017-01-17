rm(list=ls())
options(stringsAsFactors=FALSE)
options(width=300)

#INPUTS
################################################################################################################
setwd("/media/scratch/UVEITIS_GWAS/Checking_loci_for_Paul")                                                    #
x <-read.table("ASAU_and_uveitis_4eigen.assoc.logistic_rsCorrected",header=T)                                  #
num <-67478546                                                                                                 #
loci.num <-2                                                                                                   #  
cs <- 1                                                                                                        #
wS <- 300000     #window size                                                                                  #  
################################################################################################################

head(x)
x <- x[x[,"CHR"]==cs,]
x <- x[x[,"TEST"]=="ADD",]
r <- order(x[,"P"])
x<-x[r,]
head(x)

below <- num-wS
above <- num+wS
x1<-x[x["BP"]>below,]
x2<-x1[x1[,"BP"]<above,]
nrow(x2)
write.table(x2,paste("Chr",cs,"LOCUS",loci.num,"_",(num-wS/10000000),"_",(num+wS/10000000),"_preResultfile.txt",sep=""),col.names=F,row.names=F,quote=F)

file1 <- x2


file1 <- file1[,c(2,3,12)]
head(file1)
#SNP <- gsub("RS","rs",file1[,"SNP"])
#file1 <- cbind(SNP,file1[,2:3])
head(file1)
write.table(file1[1,1],file="top_snp.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
topSNP <- as.character(file1[1,1])
extract <- cbind(cs,below,above,"DP")
write.table(extract,paste("CHR",cs,"locus_",loci.num,"plink_extract_file.txt",sep=""),col.names=F,row.names=F,quote=F)


cmd <- paste("plink --noweb --allow-no-sex --bfile ASAU_and_uveitis_qc_rsCorrected --extract CHR",cs,"locus_",loci.num,"plink_extract_file.txt --range --r2 --ld-snp ", topSNP," --ld-window 99999 --ld-window-kb 1000 --ld-window-r2 0 --out ",as.character(file1[1,1]),"_ldfile",sep="" )
system(cmd)

#sort on position
Y <- order(file1[,2])
file1 <- file1[Y,]
file1$POS <- paste("chr",cs,":",file1$BP,sep="")
head(file1)

#remove old position column
file2 <- file1[,c(1,3,4)]
head(file2)

X <- complete.cases(file2)
file2 <- file2[X,]

head(file2) #96011038
tail(file2) #96493512
colnames(file2) <- c("SNP","Pval","POS")

#setwd("/home/probinson/CAST_to_LIX1_cond_2nd_corr/LZ_format_data")

write.table(file2,paste("CHR",cs,"LOCUS",loci.num,"_",below,"_",above,"_metal_file.txt",sep=""),col.names=T,row.names=FALSE,quote=FALSE)

#########################################

file <- read.table(paste(topSNP,"_ldfile.ld",sep=""),header=T)

head(file)

file$snp1 <- paste("chr",cs,":",file$BP_B,sep="")

file$snp2 <- paste("chr",cs,":",file$BP_A,sep="")

zero <- rep(0,nrow(file))

zero.df <- as.data.frame(zero,ncol=1)

head(zero.df)

LD.file <- cbind(file$snp1,file$snp2,zero,file$R2)

colnames(LD.file) <- c("snp1","snp2","dprime","rsquare")

final <- c(LD.file[1,2],LD.file[1,2],"0","1")

LD.file <- rbind(LD.file,final)

tail(LD.file)

head(LD.file)

write.table(LD.file,paste(topSNP,"_corr_format_",below,"_",above,".txt",sep=""),row.names=FALSE,quote=FALSE)

######################################################################################################

#Graphing

cmd2 <- paste("/home/probinson/software/Locus_Zoom/locuszoom/bin/locuszoom --metal CHR",cs,"LOCUS",loci.num,"_",below,"_",above,"_metal_file.txt --ld ",topSNP,"_corr_format_",below,"_",above,".txt --chr ",cs," --start ", below," --end ", above, " --pvalcol Pval --markercol POS --delim space --snpset NULL title=\"Locus ",loci.num," Chromosome ", cs, " Uveitis immunochip \"",sep="")

system(cmd2)


