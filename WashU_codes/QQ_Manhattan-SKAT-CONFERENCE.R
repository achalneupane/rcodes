#!/usr/bin/Rscript --vanilla --slave
# 2016-12-14 - vifehe
# Script to generate QQ plots and manhattanplots in a loop given a set of gene-sets (MAF, CADD ..) and a set of SKAT tests (SAKT, SKAT-O Skat Binary)
# Sets and tests need to match what is said in the input table file

	#2nd- Load the qqman package
	library("qqman")
	library(tools)
	args<-commandArgs(TRUE)
	SET<-args[1]
	TEST<-args[2]
	SKAT<-args[3]
	
	#3rd- Read tables
	graphics<-read.table(paste0(SET,"_",TEST,"_",SKAT,"-GENOME.ASSOC"), header=T)
	colnames(graphics)<-c("CHR","SNP","BP","A","F_A","F_U","A2","CHISQ","P","OR")
	graphics$CHR <- gsub("X", 23, graphics$CHR)
	graphics$CHR <- gsub("Y", 24, graphics$CHR)
	graphics$CHR <- as.numeric(graphics$CHR)
	#3rd- Perform qq-plots
	#tiff(paste0("QQplot-SKAT-",SET,"_",TEST,"_",SKAT,".tif"), units="mm", width=190, height=142, compression="lzw", res=1000)
	jpeg(paste0("QQplot-SKAT-",SET,"_",TEST,"_",SKAT,".jpg"), units="mm", width=190, height=142, res=1000)
	qq(graphics$P, main =(paste0("QQplot SKAT" ,SET,"_",TEST,"_",SKAT)))
	dev.off ()
    #4th- Perform Manhattan Plots
	#tiff(paste0("Manhattan-SKAT-",SET,"_",TEST,"_",SKAT,".tif"), units="mm", width=190, height=142, compression="lzw", res=1000)
	jpeg(paste0("Manhattan-SKAT-",SET,"_",TEST,"_",SKAT,".jpg"), units="mm", width=190, height=142, res=1000)
	manhattan(graphics, main=(paste0("Manhattan Plot of SKAT ",SET,"_",TEST,"_",SKAT)), cex=0.5, cex.axis=1, col=c("aquamarine4", "darkgoldenrod4", "blue", "blueviolet", "chartreuse4", "firebrick2", "deepskyblue", "hotpink", "yellow2", "darkolivegreen", "darkorange", "darkviolet"),suggestiveline=-log10(1e-04), genomewideline= -log10(1e-06), ylim=c(0,8))
	dev.off ()
	#5th- Perform Manhattan Plots 2 colors
	#tiff(paste0("Manhattan-SKAT-",SET,"_",TEST,"_",SKAT,"-2colors-AAIC.tif"), units="mm", width=190, height=142, compression="lzw", res=1000)
	jpeg(paste0("Manhattan-SKAT-",SET,"_",TEST,"_",SKAT,"-2colors-AAIC.jpg"), units="mm", width=190, height=142, res=1000)
	manhattan(graphics, main=(paste0("Manhattan Plot of SKAT ",SET,"_",TEST,"_",SKAT)), cex=0.5, cex.axis=1, col=c("#E69F00", "blueviolet"),suggestiveline=-log10(1e-04), genomewideline= -log10(1e-06), ylim=c(0,14))
	dev.off ()
