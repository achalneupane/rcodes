



Confidence Interval for the odd ratios is : exp[log(OR) ± 1.96 * standard error

So for the SNP rs11249215 in AS_AU group: 
log(OR) = 1.220
> Standard error (STDERR) = 
> Upper 95% Limit = exp(1.220 + 1.96 * 0.054) = 3.7653
> Lower 95% Limit = exp(1.220 - 1.96 * 0.054) = 3.0470
> 
>

OR<-1.178                                                
STDERR<-0.062
LOGOR<-log(OR)
  LOGOR                                              
exp(LOGOR + 1.96 * STDERR) 
exp(LOGOR - 1.96 * STDERR)

                                                
http://hgdownload.cse.ucsc.edu/goldenPath/hg18/snp129Mask/
 ### Recombination Rates  http://ftp.hapmap.org/recombination/latest/
 ### see /scratech/Recombination Rates
  http://cran.ms.unmelb.edu.au/bin/linux/ubuntu

##########install lodplot
install.packages("/media/Bioinform-D/Research/lodplot/lodplot_1.1.tar.gz",lib="/home/pleo/R_latest/library")
############################ to get R^2 data
## to get R^2 for a SNP:
## go to HAPMAP browser
## http://hapmap.ncbi.nlm.nih.gov/cgi-perl/gbrowse/hapmap24_B36/
##   enter region
## in reports and analysis download snp genotype data:
##   Configure : CEU
## load into haploview as HAPMPA download
## the run tagger with force include on the target snp AND get r^2 threshold to zro

 ### OR secion 14 of plink LD measurements

### see assoc_PEX.r in the PEX directory for details  <<<---- ####
## convert_mach.pl 6p22.mach1.out.mlgeno 6p22.mach1.out.mlinfo 6p22.mach1.out.mlqc -legend /media/Bioinform-D/Research/SNP\ imputation/ceu.100/chr6.CEU.legend.txt -ped  CEU.hapmap.fam  -prefix imp.6p22
## convert_mach.pl 6p22.mach1.out.mlgeno 6p22.mach1.out.mlinfo 6p22.mach1.out.mlqc -legend /media/Bioinform-D/Research/SNP\ imputation/ceu.100/genotypes_chr6_CEU_r22_nr.b36_fwd_legend.txt  -ped  CEU.hapmap.fam  -prefix imp.6p22                                                                                                                                 

## convert_mach.pl mach1.out.mlgeno mach1.out.mlinfo mach1.out.mlqc -legend /media/Bioinform-D/Research/SNP\ imputation/ceu.100/genotypes_chr6_CEU_r22_nr.b36_fwd_legend.txt -ped CEU.hapmap.fam   -prefix imp


########################################method used to get LD form HAPMAP:
## got genotypes from hapmap:
## /media/scratch/HapMap Genotypes
## subsetted out CEU other done also for other eth  and uploaded to BCGene
## CEU.hapmap.csv is from CEU.hapmap.fam which was generated from above but fam.id has been replaced by indiv.id since mach assume both are the same??

## Do mach imputation on BC gene using fastest methof and HapMap II data (Marker Map HapMap Rel #22), selected get additional files: (dataset is CEU_Hapmap table ds101538)

##  download  mach1.out.mlgeno mach1.out.mlinfo mach1.out.mlqc to this directory and then run below: (replacing old files)

convert_mach.pl mach1.out.mlgeno mach1.out.mlinfo mach1.out.mlqc -legend /media/Bioinform-D/Research/SNP\ imputation/ceu.legend/genotypes_chr4_CEU_r22_nr.b36_fwd_legend.txt -ped CEU.hapmap.csv  -prefix imp

convert_mach.pl mach1.out.mlgeno mach1.out.mlinfo mach1.out.mlqc -legend /media/Bioinform-D/Research/SNP\ imputation/ceu.legend/genotypes_chr2_CEU_r22_nr.b36_fwd_legend.txt -ped CEU.hapmap.csv  -prefix imp

convert_mach.pl mach1.out.mlgeno mach1.out.mlinfo mach1.out.mlqc -legend /media/Bioinform-D/Research/SNP\ imputation/ceu.legend/genotypes_chr6_CEU_r22_nr.b36_fwd_legend.txt -ped CEU.hapmap.csv  -prefix imp

convert_mach.pl mach1.out.mlgeno mach1.out.mlinfo mach1.out.mlqc -legend /media/Bioinform-D/Research/SNP\ imputation/ceu.legend/genotypes_chr1_CEU_r22_nr.b36_fwd_legend.txt -ped CEU.hapmap.csv  -prefix imp

convert_mach.pl mach1.out.mlgeno mach1.out.mlinfo mach1.out.mlqc -legend /media/Bioinform-D/Research/SNP\ imputation/ceu.legend/genotypes_chr11_CEU_r22_nr.b36_fwd_legend.txt -ped CEU.hapmap.csv  -prefix imp

convert_mach.pl mach1.out.mlgeno mach1.out.mlinfo mach1.out.mlqc -legend /media/Bioinform-D/Research/SNP\ imputation/ceu.legend/genotypes_chr6_CEU_r22_nr.b36_fwd_legend.txt -ped CEU.hapmap.csv  -prefix imp

 convert_mach.pl mach1.out.mlgeno mach1.out.mlinfo mach1.out.mlqc -legend /media/Bioinform-D/Research/SNP\ imputation/ceu.legend/genotypes_chr16_CEU_r22_nr.b36_fwd_legend.txt -ped CEU.hapmap.csv  -prefix imp

################# interploated p-values use mach2dat with eigenvectors as covariates. assumed to be additiative
## to get ld run in /media/scratch/HapMap Genotypes
## plink --bfile imp --r2 --ld-snp rs6710518 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out LD_GALNT3
## plink --bfile CEU.hapmap --r2 --ld-snp rs13204965 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out LD_6p22


/media/ga-apps/UQCCG-Analysis/DORITHs_WORK/DORITHs_WORK/Imputed/Association/AUAS_and_uveitis_vs_AS



## plink --file imp --r2 --ld-snp rs2849576 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out LD_ISBP
## plink --file imp --r2 --ld-snp rs1863196 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out LD_GALNT3
## plink --file imp --r2 --ld-snp rs13204965 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out LD_6p22
## plink --file imp --r2 --ld-snp rs7550034 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out LD_TGFBR3
## plink --file imp --r2 --ld-snp rs1152620 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out LD_LTBP3
## plink --file imp --r2 --ld-snp rs9466056 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out LD_SOX4

tabix -f chr6.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz
tabix chr6.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz 6:162325371-162621118 > genotypes.vcf

vcftools --gzvcf chr6.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz --geno-r2


vcftools --gzvcf chr6.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz  --plink-tped chr6.phase1_release_v3.20101123
vcftools --vcf imp.recode.vcf  --plink-tped --out imp.plink
plink --tped imp.plink.tped --tfam imp.plink.tfam --make-bed --out imp 

############# use the 1000 genomes vcf  file.
"/media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3"
##################################################################################################################################################################
vcftools --gzvcf /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr6.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz --chr 6 --from-bp 162324371 --to-bp 162625118 --out imp --recode
vcftools --gzvcf imp  --plink-tped
plink --bfile imp --r2 --ld-snp rs2849576 --ld-window-r2 0 --ld-window-kb 100000  --ld-window 99999 --out LD_PARK2.1K
plink --bfile Total_cases_samples_snps_GOOD_QC_affection --r2 --ld-snp rs2849576 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out LD_PARK2

vcftools --gzvcf /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr6.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz --chr 4 --from-bp 171618086 --to-bp 171818086 --out imp --recode
vcftools --gzvcf imp  --plink-tped
plink --bfile imp --r2 --ld-snp rs403730 --ld-window-r2 0 --ld-window-kb 100000  --ld-window 99999 --out chr4

zgrep  rs16847563 /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz
##################################################################################################################################################################
### first must get position as the VCF is in hg19 coords while the GWAS may be in hg18 ... match done via SNP name
zgrep  rs16847563 /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz

vcftools --gzvcf /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz --chr 1 --from-bp 173008854 --to-bp 175808854 --out imp --recode
vcftools --vcf imp.recode.vcf  --plink-tped --out imp.plink
plink --tped imp.plink.tped --tfam imp.plink.tfam --make-bed --out imp 
plink --bfile imp --r2 --ld-snp rs16847563 --ld-window-r2 0 --ld-window-kb 100000  --ld-window 99999 --out RC3H1.1K

## from real dataset
plink --bfile Total_cases_samples_snps_GOOD_QC_affection --r2 --ld-snp rs16847563 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out RC3H1
##################################################################################################################################################################

### first must get position as the VCF is in hg19 coords while the GWAS may be in hg18 ... match done via SNP name
zgrep  rs403730 /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr4.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz

vcftools --gzvcf /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr4.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz --chr 4 --from-bp 170481511 --to-bp 172481511 --out imp --recode
vcftools --vcf imp.recode.vcf  --plink-tped --out imp.plink
plink --tped imp.plink.tped --tfam imp.plink.tfam --make-bed --out imp 
plink --bfile imp --r2 --ld-snp rs403730 --ld-window-r2 0 --ld-window-kb 100000  --ld-window 99999 --out chr4.1K

## from real dataset
plink --bfile Total_cases_samples_snps_GOOD_QC_affection --r2 --ld-snp rs403730 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out ch4
##################################################################################################################################################################
172355309

### first must get position as the VCF is in hg19 coords while the GWAS may be in hg18 ... match done via SNP name
zgrep rs6670304 /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz

vcftools --gzvcf /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz --chr 1 --from-bp 75023350  --to-bp 77023350 --out imp --recode
vcftools --vcf imp.recode.vcf  --plink-tped --out imp.plink
plink --tped imp.plink.tped --tfam imp.plink.tfam --make-bed --out imp 
plink --bfile imp --r2 --ld-snp rs6670304 --ld-window-r2 0 --ld-window-kb 100000  --ld-window 99999 --out slc.1K

## from real dataset
plink --bfile /media/ga-apps/UQCCG-Analysis/DORITHs_WORK/DORITHs_WORK/Imputed/Association/AUAS_and_uveitis_vs_AS/Total_cases_samples_snps_GOOD_QC_affection --r2 --ld-snp rs6670304 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out slc
##################################################################################################################################################################

### first must get position as the VCF is in hg19 coords while the GWAS may be in hg18 ... match done via SNP name
zgrep rs2648465 /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr3.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz

vcftools --gzvcf /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr3.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz --chr 3 --from-bp 1 --to-bp 700000 --out imp --recode
vcftools --vcf imp.recode.vcf  --plink-tped --out imp.plink
plink --tped imp.plink.tped --tfam imp.plink.tfam --make-bed --out imp 
plink --bfile imp --r2 --ld-snp rs2648465 --ld-window-r2 0 --ld-window-kb 100000  --ld-window 99999 --out chl1.1K

## from real dataset
plink --bfile /media/ga-apps/UQCCG-Analysis/DORITHs_WORK/DORITHs_WORK/Imputed/Association/AUAS_and_uveitis_vs_AS/Total_cases_samples_snps_GOOD_QC_affection --r2 --ld-snp rs2648465 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out chl1
##################################################################################################################################################################

### first must get position as the VCF is in hg19 coords while the GWAS may be in hg18 ... match done via SNP name
zgrep rs1482071 /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr9.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz

vcftools --gzvcf /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr9.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz --chr 9 --from-bp 28594380 --to-bp 31594380 --out imp --recode
vcftools --vcf imp.recode.vcf  --plink-tped --out imp.plink
plink --tped imp.plink.tped --tfam imp.plink.tfam --make-bed --out imp 
plink --bfile imp --r2 --ld-snp rs1482071 --ld-window-r2 0 --ld-window-kb 100000  --ld-window 99999 --out chr9.1K

## from real dataset
plink --bfile /media/ga-apps/UQCCG-Analysis/DORITHs_WORK/DORITHs_WORK/Imputed/Association/AUAS_and_uveitis_vs_AS/Total_cases_samples_snps_GOOD_QC_affection --r2 --ld-snp rs1482071 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out chr9
##################################################################################################################################################################

##################################################################################################################################################################

### first must get position as the VCF is in hg19 coords while the GWAS may be in hg18 ... match done via SNP name
zgrep rs1537146 /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr9.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz

vcftools --gzvcf /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr9.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz --chr 9 --from-bp 3859303	 --to-bp 5859303 --out imp --recode
vcftools --vcf imp.recode.vcf  --plink-tped --out imp.plink
plink --tped imp.plink.tped --tfam imp.plink.tfam --make-bed --out imp 
plink --bfile imp --r2 --ld-snp rs1537146 --ld-window-r2 0 --ld-window-kb 100000  --ld-window 99999 --out rcl1l.1K

## from real dataset
plink --bfile /media/ga-apps/UQCCG-Analysis/DORITHs_WORK/DORITHs_WORK/Imputed/Association/AUAS_and_uveitis_vs_AS/Total_cases_samples_snps_GOOD_QC_affection --r2 --ld-snp rs1537146 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out rcl11
##################################################################################################################################################################
##  --ld-window-kb 100000  --ld-window 99999

### first must get position as the VCF is in hg19 coords while the GWAS may be in hg18 ... match done via SNP name
zgrep rs6590820 /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr11.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz

vcftools --gzvcf /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr11.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz --chr 11 --from-bp 99644491 --to-bp 101644491 --out imp --recode
vcftools --vcf imp.recode.vcf  --plink-tped --out imp.plink
plink --tped imp.plink.tped --tfam imp.plink.tfam --make-bed --out imp 
plink --bfile imp --r2 --ld-snp rs6590820 --ld-window-r2 0 --ld-window-kb 100000  --ld-window 99999 --out ARHGAP42.1K

## from real dataset
plink --bfile /media/ga-apps/UQCCG-Analysis/DORITHs_WORK/DORITHs_WORK/Imputed/Association/AUAS_and_uveitis_vs_AS/Total_cases_samples_snps_GOOD_QC_affection --r2 --ld-snp rs6590820 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out ARHGAP42
##################################################################################################################################################################

### first must get position as the VCF is in hg19 coords while the GWAS may be in hg18 ... match done via SNP name
zgrep rs6590820 /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr11.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz

vcftools --gzvcf /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr11.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz --chr 11 --from-bp 99644491 --to-bp 101644491 --out imp --recode
vcftools --vcf imp.recode.vcf  --plink-tped --out imp.plink
plink --tped imp.plink.tped --tfam imp.plink.tfam --make-bed --out imp 
plink --bfile imp --r2 --ld-snp rs504372 --ld-window-r2 0 --ld-window-kb 100000  --ld-window 99999 --out PGR.1K

## from real dataset
plink --bfile /media/ga-apps/UQCCG-Analysis/DORITHs_WORK/DORITHs_WORK/Imputed/Association/AUAS_and_uveitis_vs_AS/Total_cases_samples_snps_GOOD_QC_affection --r2 --ld-snp rs504372 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out PGR
##################################################################################################################################################################




### first must get position as the VCF is in hg19 coords while the GWAS may be in hg18 ... match done via SNP name
zgrep rs7187327 /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr16.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz

vcftools --gzvcf /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr16.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz --chr 16 --from-bp 73176712 --to-bp 75176712 --out imp --recode
vcftools --vcf imp.recode.vcf  --plink-tped --out imp.plink
plink --tped imp.plink.tped --tfam imp.plink.tfam --make-bed --out imp 
plink --bfile imp --r2 --ld-snp rs7187327 --ld-window-r2 0 --ld-window-kb 100000  --ld-window 99999 --out chr16.1K

## from real dataset
plink --bfile /media/ga-apps/UQCCG-Analysis/DORITHs_WORK/DORITHs_WORK/Imputed/Association/AUAS_and_uveitis_vs_AS/Total_cases_samples_snps_GOOD_QC_affection --r2 --ld-snp rs7187327 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out chr16
##################################################################################################################################################################

### first must get position as the VCF is in hg19 coords while the GWAS may be in hg18 ... match done via SNP name
zgrep rs8052043 /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr16.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz

vcftools --gzvcf /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr16.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz --chr 16 --from-bp 81000000 --to-bp 83000000 --out imp --recode
vcftools --vcf imp.recode.vcf  --plink-tped --out imp.plink
plink --tped imp.plink.tped --tfam imp.plink.tfam --make-bed --out imp 
plink --bfile imp --r2 --ld-snp rs8052043 --ld-window-r2 0 --ld-window-kb 100000  --ld-window 99999 --out CDH13.1K

## from real dataset
plink --bfile /media/ga-apps/UQCCG-Analysis/DORITHs_WORK/DORITHs_WORK/Imputed/Association/AUAS_and_uveitis_vs_AS/Total_cases_samples_snps_GOOD_QC_affection --r2 --ld-snp rs8052043 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out CDH13
##################################################################################################################################################################



### first must get position as the VCF is in hg19 coords while the GWAS may be in hg18 ... match done via SNP name
zgrep rs2836760 /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz

vcftools --gzvcf /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz --chr 21 --from-bp 39300052 --to-bp 41300052 --out imp --recode
vcftools --vcf imp.recode.vcf  --plink-tped --out imp.plink
plink --tped imp.plink.tped --tfam imp.plink.tfam --make-bed --out imp 
plink --bfile imp --r2 --ld-snp rs2836760 --ld-window-r2 0 --ld-window-kb 100000  --ld-window 99999 --out chr21.1K

## from real dataset
plink --bfile /media/ga-apps/UQCCG-Analysis/DORITHs_WORK/DORITHs_WORK/Imputed/Association/AUAS_and_uveitis_vs_AS/Total_cases_samples_snps_GOOD_QC_affection --r2 --ld-snp rs2836760 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out chr21
##################################################################################################################################################################


### first must get position as the VCF is in hg19 coords while the GWAS may be in hg18 ... match done via SNP name
zgrep rs9444597 /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr6.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz

vcftools --gzvcf /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr6.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz --chr 6 --from-bp 87987617  --to-bp 89987617 --out imp --recode
vcftools --vcf imp.recode.vcf  --plink-tped --out imp.plink
plink --tped imp.plink.tped --tfam imp.plink.tfam --make-bed --out imp 
plink --bfile imp --r2 --ld-snp rs9444597--ld-window-r2 0 --ld-window-kb 100000  --ld-window 99999 --out chr6B27.1K

## from real dataset
plink --bfile /media/ga-apps/UQCCG-Analysis/DORITHs_WORK/DORITHs_WORK/Imputed/Association/AUAS_and_uveitis_vs_AS/Total_cases_samples_snps_GOOD_QC_affection --r2 --ld-snp rs9444597 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out chr6B27
##################################################################################################################################################################


### first must get position as the VCF is in hg19 coords while the GWAS may be in hg18 ... match done via SNP name
zgrep rs4748516 /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr10.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz

vcftools --gzvcf /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr10.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz --chr 10 --from-bp 18043635  --to-bp 20043635 --out imp --recode
vcftools --vcf imp.recode.vcf  --plink-tped --out imp.plink
plink --tped imp.plink.tped --tfam imp.plink.tfam --make-bed --out imp 
plink --bfile imp --r2 --ld-snp rs4748516 --ld-window-r2 0 --ld-window-kb 100000  --ld-window 99999 --out chr10B27.1K

## from real dataset
plink --bfile /media/ga-apps/UQCCG-Analysis/DORITHs_WORK/DORITHs_WORK/Imputed/Association/AUAS_and_uveitis_vs_AS/Total_cases_samples_snps_GOOD_QC_affection --r2 --ld-snp rs4748516 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out chr10B27
##################################################################################################################################################################

### first must get position as the VCF is in hg19 coords while the GWAS may be in hg18 ... match done via SNP name
zgrep rs4748516 /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr10.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz

vcftools --gzvcf /media/ga-apps/UQCCG/GWAS/Phase1_1000Genome_v3/chr10.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz --chr 10 --from-bp 18043635  --to-bp 20043635 --out imp --recode
vcftools --vcf imp.recode.vcf  --plink-tped --out imp.plink
plink --tped imp.plink.tped --tfam imp.plink.tfam --make-bed --out imp 
plink --bfile imp --r2 --ld-snp rs4748516 --ld-window-r2 0 --ld-window-kb 100000  --ld-window 99999 --out chr10B27.1K

## from real dataset
plink --bfile /media/scratch2/cervical/LD-trim/after_QC_filter --r2 --ld-snp rs3134943 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out notch4 --noweb
##################################################################################################################################################################
                                                
## plink --tped plinkformat.tped --tfam plinkformat.tfam --make-bed --out ~/delete
## plink --bfile delete --ld rs961253 rs2423279
## plink --bfile delete --r2 --ld-snp-list list.txt --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --out ld_results

## low.cut<-162325371 ; high.cut<-162621118
# /media/ga-apps/UQCCG-Analysis/DORITHs_WORK/DORITHs_WORK/Imputed/Association/AUAS_and_uveitis_vs_AS
plink --bfile Total_cases_samples_snps_GOOD_QC_affection --r2 --ld-snp rs2849576 --ld-window-r2 0  --ld-window-kb 100000  --ld-window 99999 --out LD_PARK2
grep "rs16847563"  /media/ga-apps/UQCCG/GWAS/1000G_v2_LD/EUR/v2.20101123.EUR.chr6.xt > LD_PARK2.1K  ## but this is DPrime need allele frequencies

source("http://bioconductor.org/biocLite.R")
biocLite()

########## turn model to assoc using TREND
setwd("/media/scratch/china/New Imputation/Files for Figures")
model<-read.delim("plink.model",header=T,skip=0,sep="",fill=TRUE,stringsAsFactors=FALSE)
trend.posns<-grep("TREND",model[,"TEST"])
model<-model[trend.posns,]
write.table(model,"FINAL_TEST.assoc",row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

## GALNT3 – peak SNP rs1863196

## 6q22  - peak SNP rs13204965

 

## For the supplementary paper –

## TGFBR3 - peak SNP rs7550034

## LTPB3 - peak SNP rs1152620X1

## SOX4 - peak SNP rs9466056

## CLCN7 – peak SNP rs13336428




rs4748516 --ld-window-r2 0 --ld-window-kb 100000  --ld-window 99999 --out chr10B27.1K
#############################################chr2:166,113,288-166,565,815  rs6670304
genome.build<-"hg18"
plot.dir<-"/media/scratch2/cervical"
work.dir<-"/media/scratch2/cervical"
rates.dir<-"/media/scratch/Recombination Rates" ##genetic_map_chr1_b36.txt
file.assoc<-"cervical_best.txt" ## asscoc data from eigenstrat
target.pval.col<-"all.Pval"
#file.mach2dat <- "All.out" ## mach2dat output
got.ld<-TRUE
file.plink.ld <- "notch4.ld" ### above plink format get ld status
file.1000G.ld <- ""
with.conservation<-FALSE
the.gene <- 20 ; the.chr<-6 ; the.snp<-"rs3134943" ; low.cut<-30620500    ; high.cut<-34620500   ; left.side=FALSE; ymax.range=21; recomb.high=40
###########################################

################################# START DO ONCE ONLY #########
## setwd(plot.dir)
## dir(getwd())

## discovery<-read.delim(file.assoc,header=T,skip=0,sep="",fill=TRUE,stringsAsFactors=FALSE)
## colnames(discovery)<-c("chrom","SNP","position","A1","test","NMISS","UNKNOWN1","UNKNOWN2","P")
## discovery[,"SNP"]<-gsub("RS","rs",discovery[,"SNP"])
## snps<-as.character(discovery[,"SNP"])
## rownames(discovery)<-snps  ### needed for point types later
## discovery.ori<-discovery




setwd(work.dir)
imput<-read.delim(file.assoc,header=T,skip=0,sep="",fill=TRUE,stringsAsFactors=FALSE)

rownames(imput)<-as.character(imput[,"SNP"])  ### needed for point types later
colnames(imput)[colnames(imput)==target.pval.col]<-"Pval"
colnames(imput)[colnames(imput)=="Pos"]<-"position"
colnames(imput)[colnames(imput)=="POS"]<-"position"
colnames(imput)[colnames(imput)=="CHR"]<-"chrom"
imput[,"Pval"] <- -log10(imput[,"Pval"] )

#snps<-read.delim("TASC.bim",header=F,skip=0,sep="",fill=TRUE,stringsAsFactors=FALSE)
  snps<-imput[,"SNP"]                                   

## if(the.gene==16){snps<-read.delim("uveitis.bim",header=F,skip=0,sep="",fill=TRUE,stringsAsFactors=FALSE)}

## snps<-snps[,2]
## snps<-gsub("RS","rs",snps)
#snps<-rownames(imput)[imput[,"Genotype.all"]=="Genotyped"] ## genotyped SNPS


imput[15000,]
imput.ori<-imput
#save(list=c("imput.ori","snps","discovery.ori"),file="plot_info.RData")

#############################################################################################

## setwd(plot.dir)
## load("plot_info.RData")


####################################START HERE 

imput<-imput.ori
#discovery<-discovery.ori
############### subset imput to interested range
the.range<-(imput.ori[,"chrom"]==the.chr) & (as.numeric(imput.ori[,"position"])  >= low.cut) & (as.numeric(imput.ori[,"position"]) <= high.cut)
imput<-imput[the.range,]
dim(imput)
############ restrict range #######
imput<-imput[!is.na(imput[,"Pval"]),]
max(imput[,"position"])
min(imput[,"position"])
max(imput[,"Pval"])
imput[imput[,"Pval"]==max(imput[,"Pval"]),]

pch.array<-rep(23,times=dim(imput)[1])
posns<-match(rownames(imput),snps)
genotyped<-!is.na(posns)
pch.array[genotyped]<-19

if(got.ld){
##################################################### NO LD imformaton
setwd(work.dir)

extra<-try(read.delim(file.plink.ld,header=T,skip=0,sep="",fill=TRUE,stringsAsFactors=FALSE),silent=TRUE)
if(inherits(extra, "try-error")){extra<-{}}else{
extra[extra[,"SNP_B"]==the.snp,]
dim(extra)
}
 ## LD form 1000 genomes but does not contain positio info
extra2<-try(read.delim(file.1000G.ld,header=T,skip=0,sep="",fill=TRUE,stringsAsFactors=FALSE),silent=TRUE)

if(!inherits(extra2, "try-error")){
  
extra[1:5,]
extra2[1:5,]

# tapply(extra2[,"SNP_A"],extra2[,"SNP_A"],length)
## posns<-match(extra2[,"SNP_B"],extra[,"SNP_B"])
## cbind(extra[posns,],extra2)[1:5,]

found<-extra2[,"SNP_B"] %in% extra[,"SNP_B"]
sum(!found)

## extra4<-merge(extra,extra2,by.x="SNP_B",by.y="SNP_B",all=FALSE,all.y=TRUE,all.x=TRUE) # merges columnwise like cbin common stuff
extra<-rbind(extra,extra2[!found,])

}
}

  


## r={D}/{\sqrt{p_1 p_2 q_1 q_2} ### need the frequencies to do this

## file.1000G.ld<-"LD_PARK2.1K"  ## LD form 1000 genomes but does not contain positio info
## extra2<-try(read.delim(file.1000G.ld,header=F,skip=0,sep="",fill=TRUE,stringsAsFactors=FALSE))
## if(inherits(xx, "try-error")){extra<-{}}else{
## colnames(extra2)<-c("M1","M2","DPRIME","DELTASQ","COUPLING")
##max(extra2[,"DPRIME"])

## correct.target<-extra2[,"M2"]==the.snp
## temp<-extra2[!correct.target,"M2"]
## extra2[!correct.target,"M2"]<-extra2[!correct.target,"M1"]
## extra2[!correct.target,"M1"]<-temp

## found<-extra2[,"M1"] %in% extra[,"SNP_B"]
## sum(!found)

## posns<-match(extra2[,"M1"],extra[,"SNP_B"])
## cbind(extra[posns,],extra2)

## extra2<-extra2[!found,]


## extra.bit<-matrix(data=NA,nrow=dim(extra)[1],ncol=dim(extra)[2])
## colnames(extra.bit)<-colnames(extra)
## extra.bit[,"SNP_A"]<-extra2[,"M_2"]
## extra.bit[,"SNP_B"]<-extra2[,"M_1"]
## extra.bit[,"R2"]<-extra2[,"DPRIME"]         
##             }



###### imput, and extra may have different number of SNPs as the coveresion to the impuated data to a ped file may have caused extra missing snps
##### get only the common snps
common.snps<-intersect(rownames(imput),extra[,"SNP_B"])
length(common.snps)
dim(extra)
dim(imput)
#reorder common so in accenting order
## posns<-match(common.snps,extra[,"SNP_B"])
## extra<-extra[posns[!is.na(posns)],]
## the.order<-order(extra[,"BP_B"])
## extra<-extra[the.order,]

imput[,"R2"]<-1.1
posns<-match(rownames(imput),extra[,"SNP_B"],)
missing<-is.na(posns)
sum(missing)
lost<-imput[setdiff(1:dim(imput)[1],posns[!is.na(posns)]),]
imput[!missing,"R2"]<-as.numeric(extra[posns[!missing],"R2"])

dim(lost)
sum(!missing)
imput[!missing,][1:5,]

the.order<-order(imput[missing,][,"Pval"],decreasing=TRUE)
imput[missing,][the.order,][1:5,]

imput[missing,"R2"]<-0
#imput<-imput[!missing,]

# must be zero
##################### fix most associated SNP:
imput[the.snp,] ## must be present else lost in intersection with extra
if(imput[the.snp,"R2"]!=1){imput[the.snp,"R2"]<-1; print("WARNING")} # should be one anyway!

###########################################################################################################
##########################################################################################################



#colorss<-bluered(max(data)*4)
color = colorRampPalette(c("blue","green","red"),space="Lab",interpolate="linear")
colorss<-color(11)
colorss<-c(colorss,"black")                                       #  color = colorRampPalette(c("blue","red"))
#colorss<-colorss[c(-2,-3,-4,-5,-6)]


############################################################################### ONCE ###################
############ for archive version get host name from ensemble -> BioMart-> archive
library(biomaRt)

## if(genome.build=="hg18"){mart<-useMart("ensembl_mart_51",host="may2009.archive.ensembl.org",dataset="hsapiens_gene_ensembl",archive=TRUE)
##     }else if(genome.build=="hg19"){
##       mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
##     }else if(genome.build=="mm9"){
##       mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl");print("here:")
##  }
## #} # is alread have mart then skip
  
## mart



## listMarts(host="may2009.archive.ensembl.org",path="/biomart/martservice",archive=FALSE)
##   #http://www.ensembl.org/info/website/archives/assembly.html  may 2009 last NCBI 36 - reference 
mart=useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="may2009.archive.ensembl.org",path="/biomart/martservice",archive=FALSE)
 
#mart<-  useMart("ensembl",dataset="hsapiens_gene_ensembl")"
####################################################################################################

colors.ori<-c("darkblue","green","red","purple","forestgreen","gold3","magenta","orange","black","cyan","salmon","rosybrown4","plum4","orchid2","orangered3","olivedrab2","indianred","grey63","brown","aquamarine3","seagreen2")
## mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")


a.filter<-c( "chromosome_name", "start" , "end")
fil.vals<-list(as.character(the.chr), low.cut, high.cut)
fil.vals
### below for latest version
 exons<-getBM(attributes = c("ensembl_gene_id","external_gene_id","5_utr_start","5_utr_end","ensembl_exon_id","exon_chrom_start","exon_chrom_end","3_utr_start","3_utr_end","rank","strand","gene_biotype","chromosome_name","start_position","end_position"), filters = a.filter, values=fil.vals, mart = mart)

## fil = listFilters(mart)
## fil[grep("exon",fil[,1]),]

unique.genes<-unique(exons[,1])
unique.genes
col.array<-colors.ori[1:length(unique.genes)]
names(col.array)<-unique.genes
unique(exons[,"gene_biotype"])
keep<-exons[,"gene_biotype"] %in% c("protein_coding","miRNA")
exons<-exons[keep,]
## plus<-getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","5_utr_start","3_utr_start","ensembl_exon_id","exon_chrom_start","exon_chrom_end","rank","strand","gene_biotype"), filters = a.filter, values=fil.vals, mart = mart)


fil.vals<-list(as.character(the.chr),  low.cut, high.cut)
ann<-getBM(attributes = c( "ensembl_gene_id","external_gene_id","chromosome_name","start_position","end_position","strand","hgnc_symbol","gene_biotype"), filters = a.filter, values=fil.vals, mart = mart)
ann
keep<-ann[,"gene_biotype"] %in% c("protein_coding","miRNA")
ann<-ann[keep,]
ann

the.order<-order(exons[,"ensembl_gene_id"],exons[,"rank"])
exons<-exons[the.order,]

#################### read recombination file downloaded from   http://ftp.hapmap.org/recombination/latest/rates/
  setwd(rates.dir)
options(show.error.messages = TRUE)
file<-paste("genetic_map_chr",the.chr,"_b36.txt",sep="")
  file
recomb<-try(scan(file,what=numeric(3),skip=1,sep=" ",fill=TRUE))
num.lines<-length(recomb)/3
dim(recomb)<-c(3,num.lines)
recomb<-t(recomb)
colnames(recomb)<-c("position","rate","map")
recomb[1:5,]
the.range<-(as.numeric(recomb[,"position"])  >= low.cut) & (as.numeric(recomb[,"position"]) <= high.cut)
recomb<-recomb[the.range,]

recomb.low<-0

 setwd(plot.dir)
## recomb.high<-max(recomb[,"rate"])
## if(file2=="ANTXR2_LD_results.txt"){recomb.high<-10}
## if(file2=="ARTS1_LD_results.txt"){recomb.high<-16}

##################################

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#########################################################DO PLOT DO PLOT

if(with.conservation){
if( left.side==TRUE ){
nf<-layout(matrix(c(1,1,2,2,3,2,2,2,4,4,5,5),6,2,byrow=TRUE),heights=c(1.0,0.15,1.0,4.35,1.5,2),widths=c(1.4,2))
} else {
nf<-layout(matrix(c(1,1,2,2,2,3,2,2,4,4,5,5),6,2,byrow=TRUE),heights=c(1.0,0.15,1.0,4.35,1.5,2),widths=c(2,1.4))
}
}else{

 if( left.side==TRUE ){ 
nf<-layout(matrix(c(1,1,2,2,3,2,2,2,4,4),5,2,byrow=TRUE),heights=c(1.0,0.15,1.0,4.35,2),widths=c(1.4,2))
}else {
nf<-layout(matrix(c(1,1,2,2,2,3,2,2,4,4),5,2,byrow=TRUE),heights=c(1.0,0.15,1.0,4.35,2),widths=c(2,1.4))
}
}
  

## nf<-layout(matrix(c(1,1,5,5,2,2,6,6,2,3,6,7,2,2,6,6,4,4,8,8,   9,9,13,13,10,10,14,14,10,11,14,15,10,10,14,14,12,12,16,16)
## ,10,4,byrow=TRUE),heights=c(1.0,0.15,1.0,4.35,1.5,   1.0,0.15,1.0,4.35,1.5),widths=c(2,1.4,  2,1.4))
## nf<-layout(matrix(c(1,1,5,5,2,2,6,6,2,3,6,7,2,2,6,6,4,4,8,8,   9,9,13,13,10,10,14,14,10,11,14,15,10,10,14,14,12,12,16,16)
##,10,4,byrow=TRUE),heights=c(1.0,0.25,1.1,2.0,2.0,   1.0,0.25,1.1,2.0,2.0),widths=c(1.5,2,  1.5,2)) # 0.3 to 0.35
##  #test<-c(1,1,5,5,2,2,6,6,2,3,6,7,2,2,6,6,4,4,8,8)


layout.show(nf) # for use in mon-floating scale bar

##  nf<-layout(matrix(c(1,2,3,4),4,byrow=TRUE),heights=c(1.0,4.5,1.25,1.6))
#layout.show(nf)
par.opts.start<-par(no.readonly=TRUE)
  ############# DO BIOMART QUERYIES FIRST
require(lodplot)
par(mar=c(0,5.5,2.1,5.1),mgp=c(3,1,0)) #c(bottom, left, top, right) ## left

data(chrom.bands)
pos<-0
width<-0.25
lwd<-2
plot(c(0,max(chrom.bands[chrom.bands[,"chr"]==the.chr,"bases.bot"])),c(-1.75*width,0),pch="",axes=FALSE,xlab="",ylab="",main=paste("Chromosome",the.chr,sep=" "),cex.main=2.0)
paint.chromosomeBP(the.chr, pos = pos, width = width, bands = "major") ## function defines below:
segments(0,-1.5*width,min(ann[,"start_position"]),-1*width,lwd=lwd,col="red")
segments(min(ann[,"start_position"]),-1*width,min(ann[,"start_position"]),pos,lwd=lwd,col="red")
segments(max(ann[,"end_position"]),-1*width,max(ann[,"end_position"]),pos,lwd=lwd,col="red")
segments(max(ann[,"end_position"]),-1*width,max(chrom.bands[chrom.bands[,"chr"]==the.chr,"bases.bot"]),-1.5*width,lwd=lwd,col="red")
#layout()
#layout(matrix(1:3),heights=c(5,1,1))
#layout.show(1)    # Specify layout
#layout.show(2)    # Try specifying 
#layout.show(3)
par(mar=c(0,5.5,0,5.1),mgp=c(3,1,0))

             
the.plot<-plot(imput[,"position"],imput[,"Pval"],pch=pch.array,ylab=expression(bold(-log[10](P[val]))),xlab="",col=colorss[round(imput[,"R2"]*10,0)+1],main="",xlim=c(low.cut,high.cut),ylim=c(0,ymax.range),axes=FALSE,cex=2.0,cex.lab=2.5,font=2,font.lab=2,lwd=2) #  ,bg=colorss[round(imput[,"R2"]*10,0)+1] , xlim=c(low.cut,high.cut)
axis(2,lty=1,lwd=2,cex.axis=2.0,font=2) 



box()

################ recombination

par(mar=c(0,5.5,0,5.1),mgp=c(3,1,0),new=TRUE)  ## par(fig=c(2/3,1,2/3,1), new=T) ## also works

plot(recomb[,"position"],recomb[,"rate"],type="l",lty=2,xlab="",col="darkred",main="",axes=FALSE,cex=1.5,cex.lab=2,xlim=c(low.cut,high.cut),ylim=c(recomb.low,recomb.high),lwd=2,font.lab=2,ylab="") #  ,bg=colorss[round(imput[,"R2"]*10,0)+1]

axis(4,lty=1,lwd=2,cex.axis=2.0,col="darkred",col.axis="darkred",font=2)
mtext("Recombination Rate (cM/Mb)",side=4,padj=2.0,col="darkred",cex=1.75,las=0,font=2)          


############### ALTERNATIVE PLACEMNET METHOD
#fig=c(2/3,1,2/3,1) #NDC coordinates'c(x1, x2, y1, y2)
## shift<-0.55 # shift to right
## shift<-0    # shidt tp left
## par(fig=c(0.60-shift,1-shift,0.75,0.95),mar=c(4,0.5,3,5.5),mgp=c(3,1,0),new=TRUE) #c(bottom, left, top, right)

## par(fig=c(0.60,1,0.75,0.95),mar=c(4,0.5,3,5.5),mgp=c(3,1,0),new=TRUE)
## temp2<-matrix(0:10)  #temp2<-matrix(1:10,4,5) paste(expression(R^2),"for","rs999999",sep=" ")
## #plot(c(0,1),c(0,1),add=TRUE)
## #box()
## image(temp2, col = colorss,xlab=bquote(R^2~"with"~.(first)),main=expression("\u25CF"~~textstyle(Geneotyped~SNPs)~~~~bold("\u25C7")~~textstyle(Imputed~SNPs)),axes=FALSE,cex=1,cex.lab=2.5,cex.main=2,lwd=3,add=FALSE)
#### proble is that image not plooting in desired location

## ### NOTE USE APPLICATION -> CHARACTER MAP TO GET UNICODE ID
## axis(1,at=seq(0,1,0.1),labels=as.character(seq(0,1,0.1)),cex.axis=1.5,font=2 )
## box()


################### STOP AND CAHNGE SNP NAME
if( left.side==TRUE ){
par(mar=c(4,6.0,2,0.5),mgp=c(1.5,0,0)) #c(bottom, left, top, right) ## left
} else {
par(mar=c(4,0.5,2,5.5),mgp=c(1.5,0,0)) #c(bottom, left, top, right) ## rightpar(mar=c(4,0.5,3,5.5),mgp=c(3.1,1,0))
}


temp2<-matrix(0:10)  #temp2<-matrix(1:10,4,5) paste(expression(R^2),"for","rs999999",sep=" ")
#image(temp2, col = colorss,xlab=bquote(bold(R^2~"with"~.(first))),main="",axes=FALSE,cex=1,cex.lab=1.75,cex.main=1.5,lwd=3)
#### unicode only seems to work with 2.9.2
## image(temp2, col = colorss,xlab=bquote(bold(R^2~"with"~.(the.snp))),main=expression('{//ZapfDingbats \165}'~~bold(Genotyped~SNPs)~~~~bold("\u2666")~~bold(Imputed~SNPs)),axes=FALSE,cex=1,cex.lab=2.0,cex.main=2.0,lwd=3)


if(got.ld){
image(temp2, col = colorss[1:11],xlab=bquote(bold(R^2~"with"~.(the.snp))),main=expression("\u25CF"~~bold(Genotyped~SNPs)~~~~bold("\u25C7")~~bold(Imputed~SNPs)),axes=FALSE,cex=1,cex.lab=2.0,cex.main=2.0,lwd=3)
text(seq(0,1,0.1),0.1,labels=as.character(seq(0,1,0.1)),cex=2.0,font=2 )
box()
}else{
image(temp2, col = "white",xlab="",main=expression("\u25CF"~~bold(Genotyped~SNPs)~~~~bold("\u25C7")~~bold(Imputed~SNPs)),axes=FALSE,cex=1,cex.lab=2.0,cex.main=2.0,lwd=3)
}
##main=expression("\u25CF"~~bold(Genotyped~SNPs)~~~~bold("\u25C7")~~bold(Imputed~SNPs))
## image(matrix(0:10), col = colorss,xlab=bquote(bold(.(target))),main=expression(bold(Scale)),axes=FALSE,cex=1,cex.lab=1.75,cex.main=1.75,lwd=3)
##
## text(seq(0,1,0.1),0.1,labels=as.character(seq(0,1,0.1)),cex=2.0,font=2 )

### NOTE USE APPLICATION -> CHARACTER MAP TO GET UNICODE ID
#axis(1,at=seq(0,1,0.1),labels=as.character(seq(0,1,0.1)),cex.axis=1.0,font=2,padj=-1.0 )


############################# exons #####################
  if(with.conservation){
par(mar=c(0, 5.5,0,5.1),mgp=c(3,1,0)) #c(bottom, left, top, right)
}else{
par(mar=c(5.1,5.5,0,5.1),mgp=c(3,1,0))
}


 if(with.conservation){ 
bar.height<-0.05
bar.center<-0.08
}else{
bar.height<-0.05
bar.center<-0.08
}
  

exons.ori<-exons
ann.ori<-ann
exons<-exons.ori


#################### warning


 
utr.5.end.plus<-!is.na(exons[,"5_utr_end"]) & exons[,"strand"]==1
utr.5.end.minus<-!is.na(exons[,"5_utr_end"]) & exons[,"strand"]== -1
exons[utr.5.end.plus,"exon_chrom_start"]<-exons[utr.5.end.plus,"5_utr_end"]
exons[utr.5.end.minus,"exon_chrom_end"]<-exons[utr.5.end.minus,"5_utr_start"]

utr.3.start.plus<-!is.na(exons[,"3_utr_start"])  & exons[,"strand"]==1
utr.3.start.minus<-!is.na(exons[,"3_utr_start"])  & exons[,"strand"]== -1
exons[utr.3.start.plus,"exon_chrom_end"]<-exons[utr.3.start.plus,"3_utr_start"]
exons[utr.3.start.minus,"exon_chrom_start"]<-exons[utr.3.start.minus,"3_utr_end"]

no.nas<-!is.na(exons[,"exon_chrom_end"]) & !is.na(exons[,"exon_chrom_start"] )
exon.centers<-abs(exons[no.nas,"exon_chrom_end"]+exons[no.nas,"exon_chrom_start"])/2
y.vals<-bar.center*exons[no.nas,"strand"]

#plot(x=c(low.cut,high.cut) , y= c(-1*(bar.height+bar.center),(bar.height+bar.center)) ,xlim=c(low.cut,high.cut),type="n",axes=TRUE,ylab="",xlab="") ## protype
                                       # original is below
#plot(x=c(min(imput[,"position"]),max(imput[,"position"])) , y= c(-1*(bar.height+bar.center),(bar.height+bar.center)) ,xlim=c(low.cut,high.cut),type="n",axes=FALSE,ylab="",xlab="") ##run

if(with.conservation){
plot(x=c(min(imput[,"position"]),max(imput[,"position"])) , y= c(-1*(bar.height+bar.center),(bar.height+bar.center)) ,xlim=c(low.cut,high.cut),type="n",cex=2.0,cex.lab=2.4,lwd=3,cex.axis=2.0,axes=F,ylab="",xlab="",font=2,font.lab=2)
else{
plot(x=c(min(imput[,"position"]),max(imput[,"position"])) , y= c(-1*(bar.height+bar.center),(bar.height+bar.center)) ,xlim=c(low.cut,high.cut),type="n",cex=2.0,cex.lab=2.4,lwd=3,cex.axis=2.0,axes=FALSE,ylab="",xlab="Position (bp)",font=2,font.lab=2)
axis(1,lty=1,lwd=2,cex.axis=2.0,font=2) 

###################################################################################






 symbols(exon.centers, y.vals, rectangles=cbind(abs(exons[no.nas,"exon_chrom_end"]-exons[no.nas,"exon_chrom_start"]),bar.height ), inches=FALSE, fg =col.array[exons[no.nas,"ensembl_gene_id"]],bg =col.array[exons[no.nas,"ensembl_gene_id"]],add=TRUE,)  ### exons
############# 5 utr
no.nas<-!is.na(exons[,"5_utr_end"]) & !is.na(exons[,"5_utr_start"])
exon.centers<-abs(exons[no.nas,"5_utr_end"]+exons[no.nas,"5_utr_start"])/2
y.vals<-bar.center*exons[no.nas,"strand"]
symbols(exon.centers, y.vals, rectangles=cbind(abs(exons[no.nas,"5_utr_end"]-exons[no.nas,"5_utr_start"]),bar.height/4 ), inches=FALSE,  ,bg =col.array[exons[no.nas,"ensembl_gene_id"]],fg =col.array[exons[no.nas,"ensembl_gene_id"]],add=TRUE)  ### utr  fg =col.array[exons[no.nas,"ensembl_gene_id"]]
##########3 utr
no.nas<-!is.na(exons[,"3_utr_end"]) & !is.na(exons[,"3_utr_start"])
exon.centers<-abs(exons[no.nas,"3_utr_end"]+exons[no.nas,"3_utr_start"])/2
y.vals<-bar.center*exons[no.nas,"strand"]
symbols(exon.centers, y.vals, rectangles=cbind(abs(exons[no.nas,"3_utr_end"]-exons[no.nas,"3_utr_start"]),bar.height/4 ), inches=FALSE,  ,bg =col.array[exons[no.nas,"ensembl_gene_id"]],fg =col.array[exons[no.nas,"ensembl_gene_id"]],add=TRUE)  ### utr  fg =col.array[exons[no.nas,"ensembl_gene_id"]]

gene.centers<-abs(ann[,"start_position"]+ann[,"end_position"])/2
y.vals<-bar.center*ann[,"strand"]
the.labels<-ann[,"hgnc_symbol"]
############################### choose default label
## if(length(the.labels)==1 & is.na(the.labels[1])){the.labels<-ann[,"ensembl_gene_id"]}
## the.labels[the.labels==""]<-ann[the.labels=="","ensembl_gene_id"]

if(length(the.labels)==1 & is.na(the.labels[1])){the.labels<-ann[,"external_gene_id"]}
the.labels[the.labels==""]<-ann[the.labels=="","external_gene_id"]

keep<-rep(TRUE,times=length(the.labels))
#offsets<-rep(2.2,times=length(the.labels))
poss<-rep(3,times=length(the.labels))
to.shift.low<-FALSE
to.shift.amount.low<-0
to.shift.amount.high<-0
to.shift.high<-FALSE
set.below<-{}
to.ignore<-{}

  ############# change for a specific case below


###
## if(the.gene==2){
## to.shift.amount.high<-30000
## to.shift.high<-TRUE
## set.below<-{}
## to.ignore<-{}
## }

## if(the.gene==6){
## to.shift.low<-FALSE
## to.shift.amount.low<-30000
## to.shift.high<-FALSE
## set.below<-c(2)
## to.ignore<-c(3)
## }

## if(the.gene==1){
## to.shift.amount.high<-100000
## to.shift.high<-TRUE
## set.below<-c(8,6,3)
## to.ignore<-{}
## }

## if(the.gene==11){
## to.shift.amount.high<-0
## to.shift.high<-FALSE
## set.below<-c(6,8,3,4)
## to.ignore<-{}
## }

## if(the.gene=="6a"){
## }

if(the.gene==12){
to.shift.amount.high<-0
to.shift.high<-FALSE
set.below<-c(3)
to.ignore<-{}
}

if(the.gene==13){
to.shift.amount.low<-310000
to.shift.low<-TRUE
to.shift.high<-FALSE
set.below<-c(2)
to.ignore<-{}
}

if(the.gene==15){
to.shift.amount.high<-50000
to.shift.low<-FALSE
to.shift.high<-TRUE
set.below<-{}
to.ignore<-{}
}

the.labels
which.max(gene.centers)
gene.centers[which.max(gene.centers)]
gene.centers[which.max(gene.centers)]-100000

which.min(gene.centers)
gene.centers[which.min(gene.centers)]
gene.centers[which.min(gene.centers)]+310000

keep[to.ignore]<-FALSE  
poss[set.below]<- 1
gene.centers<-gene.centers[keep]
if(to.shift.low==TRUE){gene.centers[which.min(gene.centers)]<-gene.centers[which.min(gene.centers)]+to.shift.amount.low}
if(to.shift.high==TRUE){gene.centers[which.max(gene.centers)]<-gene.centers[which.max(gene.centers)]-to.shift.amount.high}

text(gene.centers,y.vals[keep],labels=the.labels[keep],pos=poss[keep],offset=2.0,col=col.array[ann[keep,"ensembl_gene_id"]],lwd=2,font=2,font.lab=2,cex=1.7)

num.arrows<-25
x <- seq(from=low.cut,to=high.cut,by =(high.cut-low.cut)/num.arrows)
y<- rep(bar.center,times=length(x))
s <- seq(length(x)-1)# one shorter than data
arrows(x[s], y[s], x[s+1], y[s+1], col="grey",angle=20,code=2,length=0.15,lwd=1.5)
x <- seq(from=high.cut,to=low.cut,by =(low.cut-high.cut)/num.arrows)
y<- rep(-1*bar.center,times=length(x))
arrows(x[s], y[s], x[s+1], y[s+1], col="grey",angle=20,code=2,length=0.15,lwd=1.5)
box()


################ Conservation ###########
#GENOME TABLE BROWER use region below
#paste(paste("chr",imput[1,"chrom"],sep=""),paste(min(imput[,"position"]),max(imput[,"position"]),sep="-"),sep=":")
# Comparative genomics ; 17-way Most Cons ; 
#Output as data points, set filter of points > 0 and 10,000,000 lines

if(with.conservation){                           
file<-paste("chr",the.gene,".conservation",".txt",sep="")
conserv<-read.delim(file,header=F,skip=10,sep="",fill=TRUE,stringsAsFactors=FALSE)
#conserv<-conserv[conserv[,2]!=0,]

#par(mar=c(0,5.5,0,2.1),mgp=c(3,1,0)) #c(bottom, left, top, right)
par(mar=c(5.1,5.5,0,5.1),mgp=c(3,1,0)) #c(bottom, left, top, right)


plot(x=c(low.cut,high.cut) , y= c(min(conserv[,6]),max(conserv[,6])) ,xlim=c(low.cut,high.cut),type="n",cex=2.0,cex.lab=2.4,lwd=3,cex.axis=2.0,axes=TRUE,ylab="Conservation",xlab="Position (bp)",font=2,font.lab=2)

symbols((conserv[,3]+conserv[,4])/2,conserv[,6]/2, rectangles=cbind(abs((conserv[,3]-conserv[,4])),conserv[,6] ), inches=FALSE,  ,bg ="darkblue",fg ="darkblue",add=TRUE)  ### utr
}


savePlot("chr10B.27pos.tiff",type="tiff")
savePlot("chr10B.27pos.jpeg",type="jpeg")


savePlot("chr6B.27pos.tiff",type="tiff")
savePlot("chr6B.27pos.jpeg",type="jpeg")



savePlot("chr21.tiff",type="tiff")
savePlot("chr21.jpeg",type="jpeg")

savePlot("CDH13.tiff",type="tiff")
savePlot("CDH13.jpeg",type="jpeg")



savePlot("chr16.tiff",type="tiff")
savePlot("chr16.jpeg",type="jpeg")

savePlot("PGR.tiff",type="tiff")
savePlot("PGR.jpeg",type="jpeg")


savePlot("ARCHGAP42.tiff",type="tiff")
savePlot("ARCHGAP42.jpeg",type="jpeg")


savePlot("rcl1.tiff",type="tiff")
savePlot("rcl1.jpeg",type="jpeg")


savePlot("chr9.tiff",type="tiff")
savePlot("chr9.jpeg",type="jpeg")


savePlot("chl1.tiff",type="tiff")
savePlot("chl1.jpeg",type="jpeg")

savePlot("slc44A5.tiff",type="tiff")
savePlot("slc44A5.jpeg",type="jpeg")



savePlot("chr4.tiff",type="tiff")
savePlot("chr4.jpeg",type="jpeg")


savePlot("RC3H1.tiff",type="tiff")
savePlot("RC3H1.jpeg",type="jpeg")





savePlot("PARK2.tiff",type="tiff")
savePlot("PARK2.jpeg",type="jpeg")



savePlot("IBSP.tiff",type="tiff")
savePlot("IBSP.jpeg",type="jpeg")


savePlot("CLCN7.tiff",type="tiff")
savePlot("CLCN7.jpeg",type="jpeg")

savePlot("SOX4.tiff",type="tiff")
savePlot("SOX4.jpeg",type="jpeg")

savePlot("LTBP3.tiff",type="tiff")
savePlot("LTBP3.jpeg",type="jpeg")

savePlot("TGFBR3.tiff",type="tiff")
savePlot("TGFBR3.jpeg",type="jpeg")

savePlot("rspo3.tiff",type="tiff")
savePlot("rspo3.jpeg",type="jpeg")

savePlot("galnt3.tiff",type="tiff")
savePlot("galnt3.jpeg",type="jpeg")
# savePlot("arts1.tiff",type="pdf")

save.image("erap_plot.RData")
dev.copy2pdf(file="erap.pdf")

savePlot("tnfr.tiff",type="tiff")
savePlot("tnfr.jpeg",type="jpeg")
# savePlot("arts1.tiff",type="pdf")

save.image("tnfr_plot.RData")
dev.copy2pdf(file="tnfr.pdf")


savePlot("stat3.tiff",type="tiff")
savePlot("stat3.jpeg",type="jpeg")
# savePlot("arts1.tiff",type="pdf")

save.image("stat3_plot.RData")
dev.copy2pdf(file="stat3.pdf")


##
pdf(file=)
dev.off









############################################END
################################## IDEOGRAN ##################################
 paint.chromosomeBP<-function (chrom, pos = 0, units = "cM", width = 0.4, bands = "major") 
{
    semicircle <- function(base.x, base.y, base.length, height = base.length, 
        side = 1, orientation = NULL, col = NULL) {
        radius <- base.length/2
        x <- radius * seq(-1, 1, length = 40)
        y <- height/radius * sqrt(radius^2 - x^2)
        if (is.null(orientation)) {
            co <- as.integer(cos(pi * (3 - side)/2))
            so <- as.integer(sin(pi * (3 - side)/2))
        }
        else {
            co <- cos(orientation)
            so <- sin(orientation)
        }
        tx <- co * x - so * y
        ty <- so * x + co * y
        if (is.null(orientation)) {
            if (side == 1 || side == 3) {
                base.x <- base.x + radius
            }
            else if (side == 2 || side == 4) {
                base.y <- base.y + radius
            }
        }
        x <- base.x + tx
        y <- base.y + ty
        polygon(x, y, col = col)
    }
    data(chrom.bands)
    chromdata <- subset(chrom.bands, chrom.bands$chr == chrom)
    lc <- nchar(chromdata$band)
    sel <- !(substr(chromdata$band, lc, lc) %in% letters)
    if (bands != "major") 
        sel <- !sel
    chromdata <- chromdata[sel, ]
    rm(lc, sel)
    bandcol <- gray(c(0.4, 0.6, 0.8, 0.8, 0.85))[match(chromdata$stain, 
        c("acen", "gneg", "gpos", "gvar", "stalk"))]
    n <- nrow(chromdata)
    centromere <- which(chromdata$arm[-n] != chromdata$arm[-1])
    idx <- c(2:(centromere - 1), (centromere + 2):(n - 1))
    rect(chromdata$bases.top[idx], pos, chromdata$bases.bot[idx], pos - 
        width, col = bandcol[idx])
    semicircle(chromdata$bases.bot[1], pos - width, width, chromdata$bases.bot[1] - 
        chromdata$bases.top[1], 2, col = bandcol[1])
    semicircle(chromdata$bases.top[n], pos - width, width, chromdata$bases.bot[n] - 
        chromdata$bases.top[n], 4, col = bandcol[n])
    semicircle(chromdata$bases.top[centromere], pos - width, width, 
        chromdata$bases.bot[centromere] - chromdata$bases.top[centromere], 
        4, col = bandcol[centromere])
    semicircle(chromdata$bases.bot[centromere + 1], pos - width, 
        width, chromdata$bases.bot[centromere + 1] - chromdata$bases.top[centromere + 
            1], 2, col = bandcol[centromere + 1])
    points(chromdata$bases.top[centromere], pos - 0.5 * width, col = "black", 
        cex = 3, pch = 16)
    points(chromdata$bases.top[centromere], pos - 0.5 * width, col = "white", 
        cex = 3, pch = 20)
}
########################################################################













############### BIOMART ###########
library(GenomeGraphs)

library(biomaRt)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
attributes = listAttributes(mart)
filters=listFilters(mart)


 attributes[grep("id",attributes[,1]),]
 hgnc_curated_gene_name,  hgnc_symbol hgnc_id  entrezgenegr

attributes[c(13,33,35,38,40,46,47,55,56,57,65,71,72,78), ]


           
##  "ensembl_gene_id"
## "ensembl_transcript_id"
## "ensembl_exon_id"
## "exon_chrom_start"
## "exon_chrom_end"
## "rank"
## "strand"
## "gene_biotype"
## "hgnc_symbol"
## "hgnc_id"-
### PLUS STRAND 5' -> 3'
                           name                   description
13              external_gene_id          Associated Gene Name
33 clone_based_ensembl_gene_name Clone based Ensembl gene name
35    clone_based_vega_gene_name    Clone based VEGA gene name
38                          embl             EMBL (Genbank) ID
40                          ottt  VEGA transcript ID(s) (OTTT)
46      hgnc_automatic_gene_name      HGNC automatic gene name/
47        hgnc_curated_gene_name        HGNC curated gene name/
55            mim_gene_accession            MIM Gene Accession
56          mim_gene_description          MIM Gene Description
57                       mirbase                       miRBase
65                          ucsc                       UCSC ID
71                   wikigene_id                   WikiGene ID
72                 wikigene_name                 WikiGene name
78                   dbass5_name              DBASS5 Gene Name


 fil.vals<-list(as.character(imput[1,"chrom"]), imput[1,"position"], imput[dim(imput)[1],"position"])
ann<-getBM(attributes = c( "ensembl_gene_id","external_gene_id","chromosome_name","start_position","end_position","strand","hgnc_symbol","gene_biotype","clone_based_ensembl_gene_name"),filters = a.filter, values=fil.vals, mart = mart)
ann

a.filter<-c( "chromosome_name", "start" , "end", "strand")
fil.vals<-list(as.character(imput[1,"chrom"]), imput[1,"position"], imput[dim(imput)[1],"position"],"1")
plus<-getBM(attributes = c("ensembl_gene_id","5_utr_start","5_utr_end","ensembl_exon_id","exon_chrom_start","exon_chrom_end","3_utr_start","3_utr_end","rank","strand","gene_biotype"), filters = a.filter, values=fil.vals, mart = mart)
 unique(plus[,1])
## plus<-getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","5_utr_start","3_utr_start","ensembl_exon_id","exon_chrom_start","exon_chrom_end","rank","strand","gene_biotype"), filters = a.filter, values=fil.vals, mart = mart)

fil.vals<-list(as.character(imput[1,"chrom"]), imput[1,"position"], imput[dim(imput)[1],"position"],"1")
ann.plus<-getBM(attributes = c( "ensembl_gene_id","chromosome_name","start_position","end_position","hgnc_symbol","gene_biotype"), filters = a.filter, values=fil.vals, mart = mart)

 require(stats); require(grDevices)
     x <- 1:10
     y <- sort(10*runif(10))
     z <- runif(10)
     z3 <- cbind(z, 2*runif(10), runif(10))
     symbols(x, y, thermometers=cbind(.5, 1, z), inches=.5, fg = 1:10)
     symbols(x, y, thermometers = z3, inches=FALSE)
     text(x,y, apply(format(round(z3, digits=2)), 1, paste, collapse = ","),
          adj = c(-.2,0), cex = .75, col = "purple", xpd=NA)

 arrows(x0, y0, x1, y1, length = 0.25, angle = 30, code = 2,
            col = par("fg"), lty = par("lty"), lwd = par("lwd"),
            ...)

 x <- stats::runif(12); y <- stats::rnorm(12)
     i <- order(x,y); x <- x[i]; y <- y[i]
     plot(x,y, main="arrows(.) and segments(.)")
     ## draw arrows from point to point :
     s <- seq(length(x)-1)# one shorter than data
     arrows(x[s], y[s], x[s+1], y[s+1], col= 1:3)
     s <- s[-length(s)]
     segments(x[s], y[s], x[s+2], y[s+2], col= 'pink')


##################### GENOME GRAPHS
plusStrand <- makeGeneRegion(chromosome = imput[1,"chrom"], start = imput[1,"position"], end =  imput[dim(imput)[1],"position"], strand = "+", biomart = mart,dp=DisplayPars(plotID=TRUE,idRotation=0,idColor="red",cex=0.5))

minStrand <- makeGeneRegion(chromosome = imput[1,"chrom"], start = imput[1,"position"], end =  imput[dim(imput)[1],"position"], strand = "-", biomart = mart)
ideogram <- makeIdeogram(chromosome = imput[1,"chrom"])
genomeAxis <- makeGenomeAxis(add53 = TRUE,add35=TRUE)
gdPlot(list(ideogram,plusStrand,genomeAxis, minStrand))


###2
plusStrand <- makeGeneRegion(chromosome = imput[1,"chrom"], start =  96130000, end =  96140000, strand = "+", biomart = mart, dp=DisplayPars(color="red",protein_coding="blue"))
  to <- makeTextOverlay("here is some text",xpos = 96133000, ypos = 0.5)
genomeAxis <- makeGenomeAxis(add53 = TRUE,add35=TRUE)
gdPlot(list(plusStrand,genomeAxis),overlay=c(to),add=TRUE)
text(1,1,"bob")
locator()
##############################

a <- new("GenomeAxis")
setPar(a, "size", 100) ### only gets dp parameters
gdPlot(a, minBase = 10, maxBase = 10000)


legend(x=96250000,10,legend= paste(expression(R^2),"=" ,as.character(seq(0,1,0.1)),sep=" "),col=colorss,pch=19)

> par()$mar
par(mar=c(5.1,4.1,4.1,2.1))
 par(mar=c(0,3,1,1))
     plot(x, y, xlim=xrange, ylim=yrange, xlab="", ylab="")
     par(mar=c(0,3,1,1))
     barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
     par(mar=c(3,0,1,1))


savePlot("Arts1.jpeg",type="jpeg")
print(the.plot,position = c(0, .3, 1, .9),sep=" "), more = TRUE)
     print(update(the.plot, aspect = "xy", main = "", xlab = "Year"),
           position = c(0, 0, 1, .3))
 
   print("TETS",position = c(0, 0, 1, .3))
]






########## get goodies for genome graphs
library(GenomeGraphs)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
minStrand <- makeGeneRegion(chromosome = 22,start = 30450000, end = 30550000,strand = "-", biomart = mart)
ideogram <- makeIdeogram(chromosome = 22)
genomeAxis <- makeGenomeAxis(add53 = TRUE,add35 = TRUE)
####################

########## subplot
ot(0:10, 0:10, type='n')
x <- rnorm(100)
plot(x)
par(fig=c(2/3,1,2/3,1), new=T)
hist(x, main="")
###################################

nf <- layout(matrix(c(1,2), byrow=TRUE))
layout.show(nf)
x<-1:100
y<-x^2
plot(x,y)
gdPlot(list(ideogram, plusStrand, genomeAxis, minStrand), minBase = 30450000, maxBase = 30550000)

#layout()
#layout(matrix(1:3),heights=c(5,1,1))
#layout.show(1)    # Specify layout
#layout.show(2)    # Try specifying 
#layout.show(3)
par(mar=c(1,5.5,4.1,2.1),mgp=c(3,1,0)) #c(bottom, left, top, right)
the.plot<-plot(imput[,"position"],imput[,"Pval"],pch=19,ylab=expression(-log[10](P[val])),xlab="",col=colorss[round(imput[,"R2"]*10,0)+1],main="ARTS1",axes=FALSE,cex=1.5,cex.lab=2)
axis(2,lty=1,lwd=2,cex.axis=1.75)
box()
par(mar=c(4,0,4,2.5),mgp=c(3,1,0))
temp2<-matrix(0:10)  #temp2<-matrix(1:10,4,5) paste(expression(R^2),"for","rs999999",sep=" ")
 image(temp2, col = colorss,xlab=expression(R^2~"for"~rs999999),main=expression("\u25CF"~~textstyle(Geneotyped~SNPs)~~~~"\u25C6"~~textstyle(Imputed~SNPs)),axes=FALSE,cex=1,cex.lab=2,cex.main=2,lwd=3)
### NOTE USE APPLICATION -> CHARACTER MAP TO GET UNICODE ID
axis(1,at=seq(0,1,0.1),labels=as.character(seq(0,1,0.1)),cex.axis=1.5 )
box()




axis(2, at=p, labels=format(p, decimal.mark="\u00B7"))




############### BIOMART ###########
library(GenomeGraphs)
mart <- useMart("ensembl")
mart <- useMart("ensembl_mart_50")
listMarts(archive=TRUE)
attributes =
  listDatasets(mart)
listDatasets(ensembl)
, dataset="hsapiens_gene_ensembl")
attributes = listAttributes(mart)
filters=listFilters(mart)

listMarts(archive=TRUE)
                       biomart                     version
1              ensembl_mart_51                  Ensembl 51
2                  snp_mart_51                      SNP 51
3                 vega_mart_51                     Vega 32
4              ensembl_mart_50                  Ensembl 50
5                  snp_mart_50                      SNP 50
6                 vega_mart_50                     Vega 32

Error in value[[3L]](cond) : 
  Request to BioMart web service failed. Verify if you are still connected to the internet.  Alternatively the BioMart web service is temporarily down.
In addition: Warning message:
In file(file, "r") : unable to resolve 'july2008.archive.ensembl.org'

mart<-  useMart("ensembl_mart_50",dataset="hsapiens_gene_ensembl",archive=TRUE)
mart<-  useMart("ensembl_mart_49",dataset="hsapiens_gene_ensembl",archive=TRUE)

 attributes[grep("utr",attributes[,1]),]
 hgnc_curated_gene_name,  hgnc_symbol hgnc_id  entrezgenegr

attributes[1:5, ]
##  "ensembl_gene_id"
## "ensembl_transcript_id"
## "ensembl_exon_id"
## "exon_chrom_start"
## "exon_chrom_end"
## "rank"
## "strand"
## "gene_biotype"
## "hgnc_symbol"
## "hgnc_id"-
### PLUS STRAND 5' -> 3'

a.filter<-c( "chromosome_name", "start" , "end", "strand")
fil.vals<-list(as.character(imput[1,"chrom"]), imput[1,"position"], imput[dim(imput)[1],"position"],"1")
plus<-getBM(attributes = c("ensembl_gene_id","5_utr_start","5_utr_end","ensembl_exon_id","exon_chrom_start","exon_chrom_end","3_utr_start","3_utr_end","rank","strand","gene_biotype"), filters = a.filter, values=fil.vals, mart = mart)
 unique(plus[,1])
## plus<-getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","5_utr_start","3_utr_start","ensembl_exon_id","exon_chrom_start","exon_chrom_end","rank","strand","gene_biotype"), filters = a.filter, values=fil.vals, mart = mart)

fil.vals<-list(as.character(imput[1,"chrom"]), imput[1,"position"], imput[dim(imput)[1],"position"],"1")
ann.plus<-getBM(attributes = c( "ensembl_gene_id","chromosome_name","start_position","end_position","hgnc_symbol","gene_biotype"), filters = a.filter, values=fil.vals, mart = mart)

 require(stats); require(grDevices)
     x <- 1:10
     y <- sort(10*runif(10))
     z <- runif(10)
     z3 <- cbind(z, 2*runif(10), runif(10))
     symbols(x, y, thermometers=cbind(.5, 1, z), inches=.5, fg = 1:10)
     symbols(x, y, thermometers = z3, inches=FALSE)
     text(x,y, apply(format(round(z3, digits=2)), 1, paste, collapse = ","),
          adj = c(-.2,0), cex = .75, col = "purple", xpd=NA)

##################### GENOME GRAPHS
plusStrand <- makeGeneRegion(chromosome = imput[1,"chrom"], start = imput[1,"position"], end =  imput[dim(imput)[1],"position"], strand = "+", biomart = mart,dp=DisplayPars(plotID=TRUE,idRotation=0,idColor="red",cex=0.5))

minStrand <- makeGeneRegion(chromosome = imput[1,"chrom"], start = imput[1,"position"], end =  imput[dim(imput)[1],"position"], strand = "-", biomart = mart)
ideogram <- makeIdeogram(chromosome = imput[1,"chrom"])
genomeAxis <- makeGenomeAxis(add53 = TRUE,add35=TRUE)
gdPlot(list(ideogram,plusStrand,genomeAxis, minStrand))


###2
plusStrand <- makeGeneRegion(chromosome = imput[1,"chrom"], start =  96130000, end =  96140000, strand = "+", biomart = mart, dp=DisplayPars(color="red",protein_coding="blue"))
  to <- makeTextOverlay("here is some text",xpos = 96133000, ypos = 0.5)
genomeAxis <- makeGenomeAxis(add53 = TRUE,add35=TRUE)
gdPlot(list(plusStrand,genomeAxis),overlay=c(to),add=TRUE)
text(1,1,"bob")
locator()
##############################

a <- new("GenomeAxis")
setPar(a, "size", 100) ### only gets dp parameters
gdPlot(a, minBase = 10, maxBase = 10000)


legend(x=96250000,10,legend= paste(expression(R^2),"=" ,as.character(seq(0,1,0.1)),sep=" "),col=colorss,pch=19)

> par()$mar
par(mar=c(5.1,4.1,4.1,2.1))
 par(mar=c(0,3,1,1))
     plot(x, y, xlim=xrange, ylim=yrange, xlab="", ylab="")
     par(mar=c(0,3,1,1))
     barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
     par(mar=c(3,0,1,1))


savePlot("Arts1.jpeg",type="jpeg")
print(the.plot,position = c(0, .3, 1, .9),sep=" "), more = TRUE)
     print(update(the.plot, aspect = "xy", main = "", xlab = "Year"),
           position = c(0, 0, 1, .3))
 
   print("TETS",position = c(0, 0, 1, .3))
]







library(GenomeGraphs)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
minStrand <- makeGeneRegion(chromosome = 17,start = 30450000, end = 30550000,strand = "-", biomart = mart)
ideogram <- makeIdeogram(chromosome = 17)
genomeAxis <- makeGenomeAxis(add53 = TRUE,add35 = TRUE)


gdPlot(list(ideogram, plusStrand, genomeAxis, minStrand), minBase = 30450000, maxBase = 30550000)


####################################### useful examples for plottting
X \264 ; Delta; lozenge; bullet \267; diamond \250
########### raw standard plot:
plot.new();
plot.window(c(0,10), c(0,10))
k=200
for(i in 1:10){
for(j in 1:10){
 text(i, j, expression(symbol(paste("\\",k,sep=""))))
}}

plot.new(); plot.window(c(0,4), c(15,1))
text(1, 1, "universal", adj=0); text(2.5, 1,  "\\042"); text(3, 1, expression(symbol("\250")))

the.plot<-plot(imput[,"position"],imput[,"Pval"],pch=19,col=colorss[round(imput[,"R2"]*10,0)+1],main="ARTS1",ylab=expression(-log[10](Pval)),xlab="Chromosome Position",)
legend(x=96250000,10,legend= paste(expression(R^2),"=" ,as.character(seq(0,1,0.1)),sep=" "),col=colorss,pch=19)



legend("topright",col=c("red","green","blue"), pch=c(15,16,17),
+ lty=c(1,2,3), legend=c("series 1", "series 2", "series 3"),
+ inset=0.05, bg='white')


 demo(plotmath)


dev.copy2pdf(file="blackbody.pdf")

#

matrix(1:2)
matrix(1:4)              # 4x1
matrix(1:4,2,2)          # 2x2
matrix(1:6,3,2)          # 3x2 ordered by columns
matrix(1:6,3,2,byrow=T)  # 3x2 ordered by rows
		

# To view the graphical layout, the following will show the borders of the sub-panels and the number identying each one:

layout(matrix(1:3))
layout.show(1)    # Specify layout for 4 panels, for the defined layout
layout.show(2)    # Try specifying just 2 instead
layout.show(3)		

# Now fill the layout with 4 plots:

x <- 1:10
plot(x,x)
plot(x,x^2)
plot(x,sqrt(x))
plot(x,log10(x))
curve(log10,add=T)  # Adds to last panel plotted

#####
?symbold plot symbols

The "heights" and "widths" arguments to "layout" are vectors of relative heights and widths of the matrix rows and columns, respectively.

# Specifying panels of different sizes:

 plot <- xyplot(sunspot.year ~ 1700:1988, xlab = "", type = "l",
                    scales = list(x = list(alternating = 2)),
                    main = "Yearly Sunspots")
     print(plot, position = c(0, .3, 1, .9), more = TRUE)
     print(update(plot, aspect = "xy", main = "", xlab = "Year"),
           position = c(0, 0, 1, .3))
 


library(sp)
data(meuse.grid)
coordinates(meuse.grid) <- c("x", "y")
gridded(meuse.grid) <- TRUE

With lattice graphics:

l1 <- list("SpatialPolygonsRescale", layout.scale.bar(),
  offset = c(180500,329800), scale = 500, fill=c("transparent","black"))
l2 = list("sp.text", c(180500,329900), "0")
l3 = list("sp.text", c(181000,329900), "500 m")
spplot(meuse.grid, "dist", col.regions=grey.colors(20),
  sp.layout=list(l1, l2, l3))

With base graphics:

image(meuse.grid, "dist", col=grey.colors(20))
SpatialPolygonsRescale(layout.scale.bar(), offset = c(180500, 329800),
  scale = 500, fill=c("transparent","black"), plot.grid=FALSE)
text(180500, 329900, "0")
text(181000, 329900, "500 m")


layout(matrix(1:4,2,2),heights=c(2,1)); layout.show(4)
replicate(4,plot(x,x))  # Repeat plot 4 times
with ( possum , plot ( density ( totlngth [ here ]) , type = " l " ))

    *  Spline interpolation of data

      x <- 1:20
      y <- jitter(x^2,factor=20)    # Add some noise to vector
      sf <- splinefun(x,y)          # Perform cubic spline interpolation
      # - note that "sf" is a function:
      sf(1.5)                       # Evaluate spline function at x=1.5
      plot(x,y); curve(sf(x),add=T) # Plot data points & spline curve

    * Scatter plot smoothing

      #--Using above data frame ("A"):
      plot(A$z,A$Tx)
      #--Return (X &) Y values of locally-weighted polynomial regression
      lowess(A$z,A$Tx)
      #--Plot smoothed data as line on graph:
      lines(lowess(A$z,A$Tx),col="red")

    * Numerical integration of a function (see "Functions" section below)

      # Create simple function:
      fun <- function(x,norm,index) norm*x^index
      # see "?integrate" for details. Note that:
      #  1) the names of arguments in "fun" must not match those of arguments
      #     in "integrate()" itself.
      #  2) fun must be able to return a vector, if supplied with a vector
      #     (i.e. not just a single value)
      # 
      i <- integrate(fun,lower=1,upper=10,norm=0.5,index=2)
      > i
      166.5 with absolute error < 1.8e-12
      > i[0]
      list()                         # Note that integrate() returns a list
      > names(i)                     # Show names components in list
      [1] "value"        "abs.error"    "subdivisions" "message"      "call"
      > i$value                      # If you want the integral value alone
      [1] 166.5



















