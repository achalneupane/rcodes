
## Confidence Interval 95 for the odd ratios is : exp[log(OR) ± 1.96 * standard error

## So for the SNP rs11249215 in AS_AU group: 
## log(OR) = 1.220
## > Standard error (STDERR) = 
## > Upper 95% Limit = exp(1.220 + 1.96 * 0.054) = 3.7653
## > Lower 95% Limit = exp(1.220 - 1.96 * 0.054) = 3.0470
## > 

## ### first must get position as the VCF is in hg19 coords while the GWAS may be in hg18 ... match done via SNP name
## zgrep rs4748516 /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr10.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz

## vcftools --gzvcf /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr10.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz --chr 10 --from-bp 18043635  --to-bp 20043635 --out imp --recode
## vcftools --vcf imp.recode.vcf  --plink-tped --out imp.plink
## plink --tped imp.plink.tped --tfam imp.plink.tfam --make-bed --out imp 
## plink --bfile imp --r2 --ld-snp rs4748516 --ld-window-r2 0 --ld-window-kb 100000  --ld-window 99999 --out chr10B27.1K


                                                  

################################# START DO ONCE ONLY #########
## setwd(plot.dir)
## dir(getwd())

install.packages("/media/UQCCG/GWAS/Recombination_Rates/lodplot_working.tar", repos = NULL, type="source")



require(lodplot)
library(biomaRt)   





#############################################chr2:166,113,288-166,565,815  rs6670304
genome.build<-"hg19"
path.to.UQCCG<- "/media/UQCCG" ### mofify for your machine
file.assoc<-"/media/UQCCG/GWAS/1000G_vcf_ped_gene_expression/association/hapmap3_1000G_BCL2L12_1Mbp_African_filtered.assoc.linear" ## asscoc data from eigenstrat
file.assoc<-"/media/UQCCG/GWAS/1000G_vcf_ped_gene_expression/association/hapmap3_1000G_BCL2L12_1Mbp_African_filtered_MAF0.1.assoc.linear"
 file.assoc<-"BCL2L12_African_cond.rs16981277.assoc.linear"
#file.assoc<-"BCL2L12_African_cond.rs16981277,chr19.5022.assoc.linear"

plot.name<-"African_rs2304206"

target.pval.col<-"P" ## col nae  in assoc with a the p-value

bim.file.genotyped<-"/media/UQCCG/GWAS/1000G_vcf_ped_gene_expression/association/hapmap3_1000G_BCL2L12_1Mbp_African.bim"
#file.mach2dat <- "All.out" ## mach2dat output
got.ld<-FALSE # if false will recalculate
bim.file<-"/media/UQCCG/GWAS/1000G_vcf_ped_gene_expression/association/hapmap3_1000G_BCL2L12_1Mbp_African.bim" # used to get ld


the.chr<-19 ; the.snp="rs2304206" ; low.cut<-49850000   ; high.cut<-50475000  ; left.side=FALSE; ymax.range=6; recomb.high=40 ; target.gene<-"BCL2L12" # the.snp<-"rs16981277"

label.snps<-c("rs16981277","rs11880051", "rs58705717","rs2060263","rs58705717","chr19:50228271:I","rs2304206","rs2304204")
## plot.name<-"African_zoom"
##  the.chr<-19 ; the.snp<-"rs16981277" ; low.cut<-50100000   ; high.cut<-50250000  ; left.side=FALSE; ymax.range=6; recomb.high=40 ; target.gene<-"BCL2L12"
###############################################
## system(      paste("plink","--bed",the.bed,"--bim",the.bim,"--fam",the.fam,"--filter-cases ", "--r2 --ld-snp ",ld.snps[i],"--allow-no-sex --ld-window-r2 0  --ld-window-kb 500000  --ld-window 99999 ","--hide-covar","--out",paste( "LD.with",ld.snps[i],sep=".") ,"--noweb",sep=" ") ) # no ld with dosage


#####################################################get up paths for plot and plot.data

## > imput[1:5,]
##                  chrom              SNP position A1 TEST NMISS       BETA     STAT       Pval
## rs8111874           19        rs8111874 49168942  A  ADD   326 -0.0072190 -0.27920 0.10773839
## chr19:49169216:D    19 chr19:49169216:D 49169216  Y  ADD   326 -0.0009947 -0.06019 0.02136305
## rs67824186          19       rs67824186 49170068  A  ADD   326 -0.0052480 -0.29350 0.11390427
## rs66514709          19       rs66514709 49170468  A  ADD   326 -0.0201900 -0.94480 0.46155195
## rs2353015           19        rs2353015 49171035  G  ADD   326 -0.0050640 -0.14530 0.05325306
with.conservation<-FALSE
the.gene <- 20 # tag used to do adjustments
source(paste0(path.to.UQCCG,"/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/annotate_SNPs_subroutines.r"))

if(with.conservation){
if( left.side==TRUE ){
nf<-layout(matrix(c(1,1,2,2,3,2,2,2,4,4,5,5),6,2,byrow=TRUE),heights=c(1.0,0.25,1.0,4.35,1.5,2),widths=c(1.4,2))
} else {
nf<-layout(matrix(c(1,1,2,2,2,3,2,2,4,4,5,5),6,2,byrow=TRUE),heights=c(1.0,0.25,1.0,4.35,1.5,2),widths=c(2,1.4))
}
}else{

 if( left.side==TRUE ){ 
nf<-layout(matrix(c(1,1,2,2,3,2,2,2,4,4),5,2,byrow=TRUE),heights=c(1.0,0.25,1.0,4.35,2),widths=c(1.4,2))
}else {
nf<-layout(matrix(c(1,1,2,2,2,3,2,2,4,4),5,2,byrow=TRUE),heights=c(1.0,0.25,1.0,4.35,2),widths=c(2,1.4))
}
}
  
layout.show(nf)


source(paste0(path.to.UQCCG,"/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/get_data_4_locus_data.r"))
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#########################################################DO PLOT DO PLOT

source(paste0(path.to.UQCCG,"/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/make_locusZomm_plot.r"))

######################################################### end plot

fig.prefix<-paste0(plot.name,"_ld.",the.snp)
print(fig.prefix)
getwd() 
savePlot(filename=paste(fig.prefix,".png",sep=''),type="png")
savePlot(filename=paste(fig.prefix,".jpeg",sep=''),type="jpeg")
savePlot(filename=paste(fig.prefix,".tiff",sep=''),type="tiff")
savePlot(filename=paste(fig.prefix,".bmp",sep=''),type="bmp")
dev.print(svg,paste(fig.prefix,".svg",sep=''))

dev.off


######################################################### end plot
######################################################### end plot
######################################################### end plot

the.snp


the.plot<-plot(imput[,"position"],imput[,"Pval"],pch=pch.array,ylab=expression(bold(-log[10](P[val]))),xlab="",col=colorss[round(imput[,"R2"]*10,0)+1],main="",xlim=c(low.cut,high.cut),ylim=c(0,ymax.range),axes=FALSE,cex=2.0,cex.lab=2.5,font=2,font.lab=2,lwd=2)



loc<-identify(imput[,"position"],imput[,"Pval"],imput[,"SNP"],append=TRUE)
imput[loc,"SNP"]

"rs16981277" "rs58705717","rs2060263","chr19:50228271:I","rs4287688")

"rs2060263"  "rs58705717"

condition<-c("rs16981277") #,"chr19:50228271:I")

write.table(condition,file="condition.list.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE) 

system("
 plink --bfile hapmap3_1000G_BCL2L12_1Mbp --condition-list condition.list.txt --covar African_PCA.evecs --covar-name PC1,PC2 --geno 0.05 --keep African.fam --linear hide-covar  --out BCL2L12_African_cond.rs16981277

")







################### RESULTS
africa merged : "rs16981277" "rs58705717","chr19:50228271:I","rs4287688") #,"rs2060263"
