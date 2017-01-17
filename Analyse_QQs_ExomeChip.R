rm(list=ls())
options(stringsAsFactors=FALSE)
options(width=350)
list.files(path=getwd())

#source("/media/scratch/Exome_chip/GT_Folder/QQ_plots/script.R")

par(mfrow=c(2,5))

tls <- read.table("/media/scratch/Exome_chip/Analysis_AS/Controlsnps/Exome_chip_control_snps_renamed.txt",header=F)
bad.snps <- read.table("/media/scratch/Exome_chip/GT_Folder/QQ_plots/Bad_Exome_chip_snps_renamed.txt")
data.table <- c()

y <- 2   #selected eigenvector number

for (i in 1:7) {
#     i <- 1   i <-2
  res <- read.table(paste("/media/scratch/Exome_chip/GT_Folder/Logistic_regression_by_incremental_eigenvectors/ExomeChip_incremental_",y,".assoc.logistic",sep=""),header=TRUE)
  res <- res[res[,"TEST"]=="ADD",]
  res <- res[!(is.na(res[,"P"])),]
  res <- res[!(res[,"SNP"] %in% bad.snps[,1]),]
   


  data.van <- res
  

  if (i==1) {   fexcl    <- read.table("/media/scratch/Exome_chip/GT_Folder/CreatingFreqExclusions/Exclusions_0.01_exomechip_genotypes.txt")
                data.van <- data.van[!(data.van[,"SNP"] %in% fexcl[,1]),]
                hlabel    <- "less than 0.01 excluded & 2 eigenvectors" }
  if (i==2) {   fexcl    <- read.table("/media/scratch/Exome_chip/GT_Folder/CreatingFreqExclusions/Exclusions_0.05_exomechip_genotypes.txt")
                data.van <- data.van[!(data.van[,"SNP"] %in% fexcl[,1]),]
                hlabel    <- "less than 0.05 excluded & 2 eigenvectors" }
  if (i==3) {   fexcl    <- read.table("/media/scratch/Exome_chip/GT_Folder/CreatingFreqExclusions/Exclusions_0.005_exomechip_genotypes.txt")
                data.van <- data.van[!(data.van[,"SNP"] %in% fexcl[,1]),]
                hlabel    <- "less than 0.005 excluded & 2 eigenvectors" }
  if (i==4)  {  fexcl    <-read.table("/media/scratch/Exome_chip/GT_Folder/CreatingFreqExclusions/Exclusions_0.0016_exomechip_genotypes.txt")
                data.van <- data.van[!(data.van[,"SNP"] %in% fexcl[,1]),]
                hlabel    <- "less than 0.0016 excluded & 2 eigenvectors" }
  if (i==5)  {  fexcl    <-read.table("/media/scratch/Exome_chip/GT_Folder/CreatingFreqExclusions/Exclusions_0.10_exomechip_genotypes.txt")
                data.van <- data.van[!(data.van[,"SNP"] %in% fexcl[,1]),]
                hlabel    <- "less than 0.10 excluded & 2 eigenvectors" }
     if (i==6)  {  fexcl    <-read.table("/media/scratch/Exome_chip/GT_Folder/CreatingFreqExclusions/Exclusions_0.20_exomechip_genotypes.txt")
                data.van <- data.van[!(data.van[,"SNP"] %in% fexcl[,1]),]
                hlabel    <- "less than 0.20 excluded & 2 eigenvectors" }
     if (i==7)  {  fexcl    <-read.table("/media/scratch/Exome_chip/GT_Folder/CreatingFreqExclusions/Exclusions_0.30_exomechip_genotypes.txt")
                data.van <- data.van[!(data.van[,"SNP"] %in% fexcl[,1]),]
                hlabel    <- "less than 0.30 excluded & 2 eigenvectors" }
   
#   qq(data.van[,"P"],title="test")
 print(paste("Number of rows in ",hlabel," data.van =",nrow(data.van)))
 ################
 #No exclusions QQ
 data <- data.van
 Mobspval <- sort(data$P)
 Mobspval <- Mobspval[!Mobspval==0]
 o <- -(log10(Mobspval))
 e <- -log10( 1:length(o)/length(o) )

 #Mobsmax <- 3 #trunc(max(Mlogobspval))+1
 #Mexpmax <- trunc(max(Mlogexppval))+1
 #if (is.infinite(Mobsmax)) {Mobsmax <- 3} else {Mobsmax <- Mobsmax}
 #plot(c(0,Mexpmax), c(0,Mexpmax), col="gray", lwd=1, type="l", xlab="Expected -log10 P value", ylab="Observed -log10 P value", xlim=c(0,Mexpmax), ylim=c(0,Mobsmax), las=1, xaxs="i", yaxs="i", bty="l",main=hlabel)
 #plot(c(0,Mexpmax), c(0,Mexpmax), col="gray", lwd=1, type="l", xlab="Expected -log10 P value", ylab="Observed -log10 P value", las=1, xaxs="i", yaxs="i", bty="l",main=hlabel)
 # plot(c(0,Mexpmax), c(0,Mexpmax), col="gray", lwd=1, type="l", xlab="Expected -log10 P value", ylab="Observed -log10 P value", las=1, xaxs="i", yaxs="i", bty="l",main=hlabel)
#points(Mlogexppval,Mlogobspval, pch=23, cex=.4, bg="black")
plot(e,o, pch=23, cex=.4, bg="black",main=hlabel, ylab="Observed -log10 P value",xlab="Expected -log10 P value")
 ########################################################
 # Exclude MHC
 #############
 A <- ((data[,"BP"]>20500000) & (data[,"BP"]<38500000) & (data[,"CHR"]==6))
 data.nomhc <- data[!A,]

 cMobspval <- sort(data.nomhc$P)
 cMlogobspval <- -(log10(cMobspval))
 cMexppval <- c(1:length(cMobspval))
 cMlogexppval <- -(log10( (cMexppval-0.5)/length(cMexppval)))
 points(cMlogexppval, cMlogobspval, pch = 23, cex=.4, bg="red",col="red")

 ########################################################
 # Control SNPs
 ############
 null.pos <- match(tls[,1],data[,"SNP"])
 v <- data[null.pos,]
 v <- v[!is.na(v[,1]),]

 #Second QQ
 tMobspval <- sort(v$P)
 tMlogobspval <- -(log10(tMobspval)) 
 tMexppval <- c(1:length(tMobspval))
 tMlogexppval <- -(log10( (tMexppval-0.5)/length(tMexppval))) 
 points(tMlogexppval, tMlogobspval, pch=23, cex=.4, col="blue",main=hlabel)

 legend.col <- c("black","red","blue")
 legend.words <- c("All SNPs","MHC SNPs excluded","Control SNPs only")
 x.pos <- 3
 y.pos <- 100
 #if(i==1) {x.pos <- 0.5;y.pos <- 10}
 legend(x.pos,y.pos,col=legend.col,legend=legend.words,pch=23)

}
  
###############################################
savePlot("/media/scratch/Exome_chip/GT_Folder/QQ_plots/ExomeChip_quad_QQ.jpeg",type="jpeg")
 
#source("http://people.virginia.edu/~sdt5z/0STABLE/qqman.r")
ls()

gg.qq(data$P,title="test",spartan=F)

head(res)

plot(e,o,xlim=c(0,7),ylim=c(0,7))
points(tMlogexppval, tMlogobspval, pch=23, cex=.4, col="blue",main=hlabel)
points(cMlogexppval, cMlogobspval, pch = 23, cex=.4, bg="red",col="red")

plot(e,o,xlim=c(0,10),ylim=c(0,10))
points(tMlogexppval, tMlogobspval, pch=23, cex=.4, col="blue",main=hlabel)
points(cMlogexppval, cMlogobspval, pch = 23, cex=.4, bg="red",col="red")

plot(e,o,xlim=c(0,1),ylim=c(0,1))
points(tMlogexppval, tMlogobspval, pch=23, cex=.4, col="blue",main=hlabel)
points(cMlogexppval, cMlogobspval, pch = 23, cex=.4, bg="red",col="red")

plot(o,e,xlim=c(0,50),ylim=c(0,50))
points(tMlogobspval,tMlogexppval, pch=23, cex=.4, col="blue",main=hlabel)
points(cMlogobspval,cMlogexppval, pch = 23, cex=.4, bg="red",col="red")

plot(o,e,xlim=c(0,10),ylim=c(0,10))
points(tMlogobspval,tMlogexppval, pch=23, cex=.4, col="blue",main=hlabel)
points(cMlogobspval,cMlogexppval, pch = 23, cex=.4, bg="red",col="red")

hist(o,breaks=50)
head(data.van)
nrow(data.van)
r <- order(data.van[,"P"])
d <- data.van[r,]
