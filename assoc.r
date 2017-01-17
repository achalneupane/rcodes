









### het get het and homo counts
plink --bfile data  --allow-no-sex --hardy


sed -i 's/\tRS/\trs/g' test.bim  ## RS to rsls

cut -f 1,2 -d " " ANTXR2_sanger.fam >  samples.sanger.txt ## get out list of family members.

plink --bfile  AS_ALL_common --keep  samples.sanger.txt --make-bed --out tasc


 plink --bfile AS_ALL --allow-no-sex --extract common_snp_base.txt --make-bed --out AS_ALL_common ## keep / remove for samples
                                                                                                  ## extract  / exclude for SNPs

plink --bfile tasc --extract ANTXR2_snps.csv --make-bed --out ANTXR2.tasc
cut -f 1,2 -d " " ANTXR2.tasc.fam >  samples.tasc.txt 
cut -f 2  ANTXR2.tasc.bim > ANTXR2_snps_tasc 

plink --bfile ANTXR2_sanger --extract ANTXR2_snps_tasc --keep samples.tasc.txt --make-bed --out ANTXR2_sanger_on_tasc

plink --bfile ANTXR2.tasc --bmerge ANTXR2_sanger_on_tasc.bed ANTXR2_sanger_on_tasc.bim   ANTXR2_sanger_on_tasc.fam  --make-bed --out tasc_sanger_merge


plink --bfile ANTXR2_sanger_on_tasc --flip tasc_sanger_merge.missnp --make-bed --out  

plink --bfile ANTXR2_sanger_on_tasc_F --recode --out ANTXR2_sanger_on_tasc_F
plink --bfile ANTXR2.tasc --recode --out ANTXR2.tasc


plink --bfile ALL --update-ids recode.txt --make-bed --out  All_tascIDs  ## used match.sample.ids.r to make a 4 colume recode file of sample ids 

plink --bfile All_tascIDs --extract ANTXR2_snps.csv --keep  samples.sanger.txt --make-bed --out All_on_sanger



plink --bfile WTTC2 --update-ids recode_all.txt --make-bed --out  WTCC2_WTCCIDs 
plink --bfile WTCC2_WTCCIDs  --extract ANTXR2_snps.csv --keep  samples.sanger.txt --make-bed --out WTCC2_on_sanger


plink --bfile ALL --update-ids recode_all.txt --make-bed --out  All_WTCCIDs 
plink --bfile All_WTCCIDs  --extract ANTXR2_snps.csv --keep  samples.sanger.txt --make-bed --out All_on_sanger

 plink --bfile All_on_sanger --bmerge WTCC2_on_sanger.bed WTCC2_on_sanger.bim   WTCC2_on_sanger.fam  --make-bed --out All_WTCC2_merge ##check is a flip needed
plink --bfile All_on_sanger --recode --out All_on_sanger
plink --bfile WTCC2_on_sanger --recode --out WTCC2_on_sanger  ## convert to ped


 plink --bfile All_on_sanger --bmerge ANTXR2_sanger.bed ANTXR2_sanger.bim   ANTXR2_sanger.fam  --make-bed --out Allcat _sanger_merge

plink --bfile ANTXR2_sanger --flip All_sanger_merge.missnp --make-bed --out ANTXR2_sanger_F
 plink --bfile ANTXR2_sanger_F --recode --out ANTXR2_sanger_F



plink --bfile All_WTCCIDs  --extract ANTXR2_snps.csv --keep samples_WTCC2ids --make-bed --out All_on_WTCC
plink --bfile WTCC2_WTCCIDs  --extract ANTXR2_snps.csv --keep samples_WTCC2ids --make-bed --out WTCC2_on_WTCC
plink --bfile All_on_WTCC --recode --out All_on_WTCC
plink --bfile WTCC2_on_WTCC --recode --out WTCC2_on_WTCC

plink --bfile total_common_QC --update-ids recode_all.txt --make-bed --out  QCed_WTCCIDs 
sed -i 's/\tRS/\trs/g' QCed_WTCCIDs.bim 


plink --bfile QCed_WTCCIDs  --extract ANTXR2_snps.csv --keep samples_WTCC2ids --make-bed --out Qced_on_WTCC
plink --bfile All_on_sanger --bmerge Qced_on_WTCC.bed Qced_on_WTCC.bim   Qced_on_WTCC.fam  --make-bed --out All_QCed_merge
plink --bfile Qced_on_WTCC --flip All_QCed_merge.missnp --make-bed --out Qced_on_WTCC_F
plink --bfile Qced_on_WTCC_F --recode --out Qced_on_WTCC_F


plink --bfile ALL --extract ANTXR2_snps.csv --make-bed --out ALL_antrx

sed -i 's/\tRS/\trs/g' data.bim
plink --bfile data --extract ANTXR2_snps.csv --make-bed --out BC_antrx

plink --bfile ALL_antrx --bmerge BC_antrx.bed BC_antrx.bim   BC_antrx.fam  --make-bed --out All_BC
plink --bfile  BC_antrx --recode --out  BC_antrx
plink --bfile  ALL_antrx --recode --out  ALL_antrx

sed -i 's/\tRS/\trs/g' total_common_QC.bim
plink --bfile total_common_QC --extract test_snps.csv --make-bed --out Qced_test
plink --bfile ALL --extract test_snps.csv --make-bed --out All_test

plink --bfile All_test --bmerge Qced_test.bed Qced_test.bim   Qced_test.fam  --make-bed --out All_QCed_merge
plink --bfile Qced_test --flip All_QCed_merge.missnp --make-bed --out Qced_test_F
plink --bfile  Qced_test_F --recode --out  Qced_test_F
plink --bfile  All_test --recode --out  All_test

plink --bfile WTCC2_WTCCIDs  --extract ANTXR2_snps.csv --keep  samples.sanger.txt --make-bed --out WTCC2_on_sanger


plink --bfile AS_CAN  --update-map hg18_650Y_chr.txt --update-chr --make-bed --out AS_CAN_OR




plink --bfile AS_European_after_hapmap_nofp_diag_check2_wAlleles  --bmerge total_common_QC.bed  total_common_QC.bim  total_common_QC.fam --make-bed --out Ichip_QC_merge
plink --bfile AS_European_after_hapmap_nofp_diag_check2_wAlleles  --flip Ichip_QC_merge.missnp  --make-bed --out Ichip_F

plink --bfile Ichip_F  --bmerge total_common_QC.bed  total_common_QC.bim  total_common_QC.fam --make-bed --out Ichip_QC_merge
plink --bfile Ichip_QC_merge --geno 0.05 --hwe 0.00001 --make-bed --out Ichip_QC_merge_f --noweb
plink --bfile  Ichip_QC_merge_f --genome --genome-full --out Ichip_QC_merge_f.genome

perl -ne 'BEGIN { $i=0} ;@w=split;$i++; if ( $w[9] > 0.25 | $i == 1) {print $w[1]."\t".$w[3]."\t".$w[9]; print "\n"} ' Ichip_QC_merge_f.genome.genome >  Ichip_QC_merge_f.related
perl -ne 'BEGIN { $i=0} ;@w=split;$i++; if ( $w[2] > 0.98 | $i == 1) {print $w[0]."\t".$w[1]."\t".$w[2]; print "\n"} ' Ichip_QC_merge_f.related > Ichip_QC_merge_f.identical

5439670050_R04C02 5439670050_R04C02 0 0 0 1 ### control
UK_CA_B1183 UK_CA_B1183 0 0 0 2             ### control become a case            
plink --bfile Ichip_F --remove samples.1.txt  --make-bed --out  Ichip_FC  # remove control->case event


plink --bfile Ichip_FC --update-ids recode_Ichip_all.txt --make-bed --out  Ichip_F_QCids 


plink --bfile Ichip_F_QCids  --update-map  snp_recode_ichip.txt --update-name --make-bed --out  Ichip_FCS_QCids ### fix the immunochip ids.

plink --bfile total_common_QC  --extract test_snps.csv  --keep samples_common_ichip_QC.txt  --make-bed --out total_common_QCids_test
cut -f 2 total_common_QCids_test.bim > test_snps_inQC.csv

plink --bfile Ichip_FCS_QCids  --extract test_snps_inQC.csv  --keep samples_common_ichip_QC.txt  --make-bed --out Ichip_FCS_QCids_test

plink --bfile Ichip_FCS_QCids --bmerge total_common_QCids_test.bed total_common_QCids_test.bim total_common_QCids_test.fam --make-bed --out merge
plink --bfile Ichip_FCS_QCids_test --flip merge.missnp --make-bed --out Ichip_FFCS_QCids_test

plink --bfile  total_common_QCids_test --recode --out  total_common_QCids_test
plink --bfile  Ichip_FFCS_QCids_test --recode --out  Ichip_FFCS_QCids_test


plink --bfile Ichip_F_QCids  --extract Ichip_duplicate_snps.txt --make-bed --out  Ichip_FCS_QCids
plink --bfile Ichip_F_QCids  --update-map  snp_recode_ichip.txt --update-name --make-bed --out  Ichip_FCS_QCids
plink --bfile Ichip_FCS_QCids   --keep samples_common_ichip_QC.txt  --make-bed --out Ichip_FCS_QCids_UV





	--bfile ASAUvsControls
	--logistic
	--covar ASAU_covariates.txt
	--ci 0.95
	--allow-no-sex
	--out ASAU_vs_Controls

--standard-beta

cut -f 1,2 -d " " ASAU_NEGvsControls.fam > samples_ASAU_NEGvsControls.txt


plink --bfile total_common_QC --keep samples_ASAU_POSvsControls.txt --make-bed --out ASAU_POS
plink --bfile total_common_QC --keep samples_ASAU_NEGvsControls.txt --make-bed --out ASAU_NEG


plink	--bfile ASAU_POS  --logistic --standard-beta --covar ASAU_POScovariates.txt --ci 0.95 --allow-no-sex --out ASAU_POS_stB
plink	--bfile ASAU_POS  --logistic --covar ASAU_POScovariates.txt --ci 0.95 --allow-no-sex --out ASAU_POS
plink	--bfile ASAU_POS  --assoc  --ci 0.95 --allow-no-sex --out ASAU_POS


#plink	--bfile ASAU_POS  --logistic --standard-beta --covar ASAU_POScovariates.txt --ci 0.95 --allow-no-sex --out ASAU_POS_stB
plink	--bfile ASAU_NEG  --logistic --covar ASAU_NEGcovariates.txt --ci 0.95 --allow-no-sex --out ASAU_NEG
plink	--bfile ASAU_NEG  --assoc  --ci 0.95 --allow-no-sex --out ASAU_NEG



plink	--bfile ASAU_POS --remove strange.QC.samples.txt --make-bed --out ASAU_POS_trim
plink	--bfile ASAU_POS_trim  --logistic --covar ASAU_POScovariatesTrim.txt  --ci 0.95 --allow-no-sex --out ASAU_POS_trim
plink	--bfile ASAU_POS_trim  --assoc  --ci 0.95 --allow-no-sex --out ASAU_POS_trim

plink	--bfile ASAU_NEG --remove strange.QC.samples.txt --make-bed --out ASAU_NEG_trim
plink	--bfile ASAU_NEG_trim  --logistic --covar ASAU_NEGcovariatesTrim.txt  --ci 0.95 --allow-no-sex --out ASAU_NEG_trim
plink	--bfile ASAU_NEG_trim  --assoc  --ci 0.95 --allow-no-sex --out ASAU_NEG_trim



/media/scratch2/AS_CANADIAN_GWAS
prahman@mun.ca
Guangju.Zhai@med.mun.ca

http://di-genetics-wiki.di.uq.edu.au/mediawiki/index.php/Stratification
#### not that the AS CAN gwas was on hg19 but the 650Y are on hg18 this is how you update the map:

cut -f 2,4 ALL_eth_650Y.bim >  hg18_650Y_map.txt
cut -f 2 ALL_eth_650Y.bim >  the_rs.txt
cut -f 1 ALL_eth_650Y.bim >  the_chr.txt

paste -d "\t" the_rs.txt the_chr.txt > hg18_650Y_chr.txt

plink --bfile AS_case_control_542114_merged_final_clean --update-map hg18_650Y_map.txt --make-bed --out AS_CAN 

plink --bfile AS_CAN  --update-map hg18_650Y_chr.txt --update-chr --make-bed --out AS_CAN_OR

plink --bfile AS_CAN_OR --make-bed --out AS_CAN

# plink --bfile AS_CAN  --make-bed --out AS_CAN_OR
# sed -i 's/rs/SNP/g' $SCRATCH/haplo$CHR.assoc.head
plink 

########################

plink --bfile AS_CAN_OR  --bmerge ALL_eth_650Y.bed  ALL_eth_650Y.bim  ALL_eth_650Y.fam --make-bed --out AS_CAN_650

plink --bfile AS_CAN_650 --geno 0.03 --hwe 0.0000001 --make-bed --out AS_CAN_650_f --noweb

plink --bfile AS_CAN_650_f --exclude Price2008.txt --range --make-bed --out AS_CAN_650_f_ld --noweb

### NOW LD TRIM USE THIS:
http://di-genetics-wiki.di.uq.edu.au/mediawiki/index.php/Stratification
 plink --file data --indep 50 5 1.5 --out ldtrimset  ## identify markers to priuse
plink --file data --exclude ldtrimset.prune.out ## remove tose markers
plink --file data --thin .3  ## additionall if need thin out randomly (removes 30% more)

## NOW NEED TO DO A QC SET AGAIN 

plink --bfile  AS_CAN_650_f_ld --allow-no-sex --maf 0.05 --max-maf 0.95 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out  AS_CAN_650_f_ld_f

nohup /media/Bioinform-D/Research/EIGEN/bin/smartpca.perl -i AS_CAN_650_f_ld_f.bed -a AS_CAN_650_f_ld_f.bim -b AS_CAN_650_f_ld_f.fam -o AS_CAN_650_f_ld_f.pca -p AS_CAN_650_f_ld_f.plot -e AS_CAN_650_f_ld_f.eval -m 0 -l AS_CAN_650_f_ld_f.log  > AS_CAN_650_f_ld_f_m0.log

plink --bfile  AS_CAN_650_f_ld_f --genome --genome-full --out AS_CAN_650_f_ld_f.genome


#####UViteis
nohup /media/Bioinform-D/Research/EIGEN/bin/smartpca.perl -i data_no_lrld_clustering3f.bed -a data_no_lrld_clustering3f.bim -b data_no_lrld_clustering3f.fam -o data_no_lrld_clustering3f.pca -p data_no_lrld_clustering3f.plot -e data_no_lrld_clustering3f.eval -m 0 -l data_no_lrld_clustering3f.log  > output.lstxt



RELATEDNESS 

plink --bfile data --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --allow-no-sex --genome --genome-full


perl -ne '@w=split; if ($w[8] > 0.5) {print $_} ' plink.genome
perl -ne 'BEGIN { $i=0} ;@w=split;$i++; if ( $w[9] > 0.5 | $i == 1) {print $w[1]."\t".$w[3]."\t".$w[9]; print "\n"} ' plink.genome


perl -ne 'BEGIN { $i=0} ;@w=split;$i++; if ( $w[2] > 0.98 | $i == 1) {print $w[1]."\t".$w[2]."\t".$w[3]; print "\n"} ' test
#################### NEW WORK ##########################


################# flit starnds of icontrolsDB , to_flip filtes made in R /matt brown/TOPBOT
################# from illumina annottaions updated to build 36
################# bed dumped from BCSnp using 128 snp build
#Files in D:\Research\Matt Brown\WTCCC

plink --bfile Controls_300v1_07 --flip to_flip_300v1.txt --make-bed --allow-no-sex --out Controls_300v1_07_f
plink --bfile Controls_300v1lip to_flip_300v1.txt --make-bed --allow-no-sex --out Controls_300v1_08_f
plink --bfile Controls_1Mv1 --flip to_flip_1Mv1.txt --make-bed --allow-no-sex --out Controls_1Mv1_f
plink --bfile Controls_550v1 --flip to_flip_550v1.txt --make-bed --allow-no-sex --out Controls_550v1_f
plink --bfile Controls_550v3 --flip to_flip_550v3.txt --make-bed --allow-no-sex --out Controls_550v3_f
plink --bfile Controls_550v3duo --flip to_flip_550v3.txt --make-bed --allow-no-sex --out Controls_550v3duo_f
plink --bfile Controls_610quad --flip to_flip_610quad.txt --make-bed --allow-no-sex --out Controls_610quad_f
plink --bfile hapmap --flip to_flip_300v1.txt --make-bed --allow-no-sex --out hapmap_f   #add hapmap

########### merge controls
plink --bfile Controls_WTCCC --bmerge Controls_550v3_f.bed  Controls_550v3_f.bim  Controls_550v3_f.fam --merge-mode 7
plink --bfile Controls_ALL --bmerge hapmap_f.bed  hapmap_f.bim  hapmap_f.fam  --allow-no-sex --make-bed --out Controls_ALL_hapmap
### Note: Merge Mode:     7    Report mismatching non-missing calls (diff mode -- do not merge)

##########Below should merge... missing snps in one setA comapered to set B set as missing in Controls_ALL
plink --merge-list merge_controls_list.txt --allow-no-sex --make-bed --out Controls_ALL

plink --bfile Controls_ALL_hapmap --bmerge cases.bed  cases.bim  cases.fam  --allow-no-sex --make-bed --out AS_ALL


######################################################################################################
######################################################################################################
###################  ALL ethnicity on 650 Y chips
plink --bfile Controls_WTCCC --bmerge ALL_eth_650Y_f.bed  ALL_eth_650Y_f.bim  ALL_eth_650Y_f.fam  --allow-no-sex --make-bed --out ALL_eth_650Y_WTCCC
plink --bfile  ALL_eth_650Y --flip ALL_eth_650Y_WTCCC.missnp --make-bed --allow-no-sex --out ALL_eth_650Y_WTCCC_f

plink --bfile ALL_eth_650Y_WTCCC --extract common_snp_base.txt --make-bed --out ALL_eth_650Y_WTCCC_C
plink --bfile ALL_eth_650Y_WTCCC_C --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out ALL_eth_650Y_WTCCC_C_F


nohup ~/bin/smartpca.perl -i ALL_eth_650Y_WTCCC_C_F.bed -a ALL_eth_650Y_WTCCC_C_F.bim -b ALL_eth_650Y_WTCCC_C_F.fam -m 0 -o ALL_eth_650Y_WTCCC_C_F.pca -p ALL_eth_650Y_WTCCC_C_F.plot -e ALL_eth_650Y_WTCCC_C_F.eval -l ALL_eth_650Y_WTCCC_C_F.log  > output.txt
 plink --bfile  ALL_eth_650Y_WTCCC --extract ethnicity_snps_avaliable.txt --allow-no-sex --make-bed --out ALL_eth_650Y_WTCCC_eth

############ put in for fo reading from R:
plink --bfile  ALL_eth_650Y_WTCCC_eth --recodeAD --out ALL_eth_650Y_WTCCC_eth_R
 head -1 ALL_eth_650Y_WTCCC_eth_R.raw > first.line.txt
 wc -w first.line.txt # 104: 2x(49 Snps in list) +6
 perl -pe 'BEGIN { @index=(1); for($i=6;$i<=104;$i=$i+2){push(@index,$i)} }   @w=split; $_=join("\t",@w[@index]); print "\n" ' ALL_eth_650Y_WTCCC_eth_R.raw > ALL_eth_650Y_WTCCC_eth_R_calls.txt
######################################################################################################
######################################################################################################
plink --bfile Controls_WTCCC --bmerge Aff_get_58_Chiamo_ALL.bed  Aff_get_58_Chiamo_ALL.bim  Aff_get_58_Chiamo_ALL.fam  --allow-no-sex --make-bed --out Aff_get_58_Chiamo_ALL_WTCCC
plink --bfile  Aff_get_58_Chiamo_ALL --flip Aff_get_58_Chiamo_ALL_WTCCC.missnp --make-bed --allow-no-sex --out Aff_get_58_Chiamo_ALL_f
plink --bfile Aff_get_58_Chiamo_ALL_f --exclude double_snp.txt --make-bed  --out Aff_get_58_Chiamo_ALL_REV
plink --bfile  Aff_get_58_Chiamo_ALL  --exclude double_snp.txt --make-bed  --out Aff_get_58_Chiamo_ALL_FWD


plink --bfile Controls_WTCCC --bmerge Aff_get_58_Chiamo_ALL_REV.bed  Aff_get_58_Chiamo_ALL_REV.bim  Aff_get_58_Chiamo_ALL_REV.fam  --make-bed --out Aff_get_58_Chiamo_ALL_WTCCC


#  Merging 8 samples, final sample contains 5675 individuals and 1069443 markers
# Before frequency and genotyping pruning, there are 1069443 SNPs
# 5675 founders and 0 non-founders found
# 5399 SNPs with no founder genotypes observed
# Warning, MAF set to 0 for these SNPs (see --nonfounders)
# Writing list of these SNPs to [ Controls_ALL.nof ]
# Total genotyping rate in remaining individuals is nan
# 0 SNPs failed missingness test ( GENO > 1 )
# 0 SNPs failed frequency test ( MAF < 0 )
# After frequency and genotyping pruning, there are 1069443 SNPs
# After filtering, 0 cases, 5675 controls and 0 missing
# After filtering, 0 males, 0 females, and 5675 of unspecified sex
# Writing pedigree information to [ Controls_ALL.fam ]
# Writing map (extended format) information to [ Controls_ALL.bim ] 
# Writing genotype bitfile to [ Controls_ALL.bed ] 
# Using (default) SNP-major mode

Analysis finished: Fri Sep 12 16:02:25 2008


############# filter    common_snp_base has 301866 SNPs have 5675 controls and 2181 cases
  plink --bfile AS_ALL --allow-no-sex --extract common_snp_base.txt --make-bed --out AS_ALL_common
  
  ############## find genome
 nohup plink --bfile AS_ALL_common --allow-no-sex --genome --out AS_ALL_common_G > AS_ALL_common_G_lines &

 nohup plink --bfile AS_ALL_common --keep test_genome.txt --allow-no-sex  --genome --out test2_genome_G > test2_genome_lines &
#  /home/pleo/EIGENSTRAT/MERGE ALL 08-08/AS_ALL_common.fam
# /home/pleo/EIGENSTRAT/MERGE ALL 08-08/AS_ALL_common_G.genome
# /home/pleo/EIGENSTRAT/MERGE ALL 08-08/AS_ALL_common_G.nosex
# /home/pleo/EIGENSTRAT/MERGE ALL 08-08/AS_ALL_common.bed
# /home/pleo/EIGENSTRAT/MERGE ALL 08-08/AS_ALL_common.bim
# /home/pleo/EIGENSTRAT/MERGE ALL 08-08/AS_ALL_common_G.log
# /home/pleo/EIGENSTRAT/MERGE ALL 08-08/AS_ALL_common_G_line

perl -pe '@w=split; $_=$w[1]."\t".$w[3]."\t".$w[7]; print "\n" ' AS_ALL_common_G.genome  > AS_ALL_common_G_piHAT.txt

perl -ne '@w=split; if ($w[2] > 0.05) {print $w[0]."\t".$w[1]."\t".$w[2], "\n"} ' AS_ALL_common_G_piHAT.txt > AS_ALL_common_G_related.txt
perl -ne '@w=split; if ($w[7] > 0.04) {print $_} ' AS_ALL_common_G.genome  >  AS_ALL_common_G_related_all_cols.txt
perl -ne '@w=split; if ($w[7] > 2.5) {print $_} ' AS_ALL_common_G.genome  >  test_original_genome.txt

perl -ne '@w=split; if ($w[9] > 0.01) {print $_} ' AS_ALL_common_G.genome  >  AS_ALL_common_G_related_all_cols.txt
#perl -ne '@w=split; @w= sort { $a <=> $b } @w ; {print $w[0], "\n"} ' 


# bzip2 --keep --best    AS_ALL_common_G.genome
# bunzip2 --
# 
# 
#                        Z0           Z1           Z2           PI_hat
# MZ                        0              0              1              1
# Sibs                        0.25        0.5          0.25        0.5
# Parent-offspring              0              1              0              0.5
# Half-sibs, g-parent/g-child, aunts/uncles               0.5          0.5          0              0.25
# 1st cousins           0.75        0.25        0              0.125
#
#   Ratio > 2.5 will be used

## use remove related.r -> related_samples.txt        416 samples in here


#plink --bfile AS_ALL_common --allow-no-sex --exclude p6_remove.txt  --remove related_samples.txt --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --make-bed --out merge_all_final_cSNP_no6pF_mafF3



plink --bfile AS_ALL_common --allow-no-sex  --remove related_samples.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out AS_ALL_common.noR.F

## 301866 markers to
#convert f2c library to static   sudo ln -s /usr/lib/libf2c.a /usr/local/lib/libf2c.a
# then changed makile -lf2c to /usr/local/lib/libf2c.a seemed to work...


## run smartpca
nohup ~/bin/smartpca.perl -i AS_ALL_common.noR.F.bed -a AS_ALL_common.noR.F.bim -b AS_ALL_common.noR.F.fam -o AS_ALL_common.noR.F.pca -p AS_ALL_common.noR.F.plot -e AS_ALL_common.noR.F.eval -l AS_ALL_common.noR.F.log  > output.txt


# ** For gPLINK compatibility, do not use '.' in --out **
# Reading map (extended format) from [ AS_ALL_common.bim ] 
# 301866 markers to be included from [ AS_ALL_common.bim ]
# Reading pedigree information from [ AS_ALL_common.fam ]
# 8028 individuals read from [ AS_ALL_common.fam ] 
# 8028 individuals with nonmissing phenotypes
# Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)
# Missing phenotype value is also -9
# 2181 cases, 5847 controls and 0 missing
# 0 males, 0 females, and 8028 of unspecified sex
# Warning, found 8028 individuals with ambiguous sex codes
# Writing list of these individuals to [ AS_ALL_common.noR.F.nosex ]
# Reading genotype bitfile from [ AS_ALL_common.bed ] 
# Detected that binary PED file is v1.00 SNP-major mode
# Reading individuals to remove [ related_samples.txt ] ... 416 read
# 416 individuals removed with --remove option
# Before frequency and genotyping pruning, there are 301866 SNPs
# 7612 founders and 0 non-founders found
# 0 of 7612 individuals removed for low genotyping ( MIND > 0.1 )
# 11414 markers to be excluded based on HWE test ( p <= 1e-07 )
#         9019 markers failed HWE test in cases
#         11414 markers failed HWE test in controls
# Total genotyping rate in remaining individuals is 0.996359
# 2327 SNPs failed missingness test ( GENO > 0.05 )
# 134 SNPs failed frequency test ( MAF < 0.01 )
# After frequency and genotyping pruning, there are 288526 SNPs
# After filtering, 2119 cases, 5493 controls and 0 missing
# After filtering, 0 males, 0 females, and 7612 of unspecified sex
# Writing pedigree information to [ AS_ALL_common.noR.F.fam ]
# Writing map (extended format) information to [ AS_ALL_common.noR.F.bim ] 
# Writing genotype bitfile to [ AS_ALL_common.noR.F.bed ] 
# Using (default) SNP-major mode
# 
# Analysis finished: Fri Oct 17 18:54:46 2008                   plink --bfile AS_ALL_common.noR.F --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1  --out AS_ALL_common.noR.F.Assoc


 ## 301866 markers to
#convert f2c library to static   sudo ln -s /usr/lib/libf2c.a /usr/local/lib/libf2c.a
# then changed makile -lf2c to /usr/local/lib/libf2c.a seemed to work...
# libf2c call libg2c now
# eigenstat on red hat covered llpack and blas the same way

## run smartpca
nohup ~/bin/smartpca.perl -i AS_ALL_common.noR.F.bed -a AS_ALL_common.noR.F.bim -b AS_ALL_common.noR.F.fam -o AS_ALL_common.noR.F.pca -p AS_ALL_common.noR.F.plot -e AS_ALL_common.noR.F.eval -l AS_ALL_common.noR.F.log  > output.txt

###run plink to get snp list identift 6p snps
plink --bfile AS_ALL_common.noR.F --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1  --out AS_ALL_common.noR.F.Assoc
#chr 6 0-60000000 -   used excell    - p6_remove.txt
# added 12 known bad genothype (from  BC) - bad.snps
# combined is p6_remove_bas_snps.txt

## to get list of all samples:
plink --bfile AS_ALL_common --test-missing --out   AS_ALL_common_missingCC     #  8028 founders and 0 non-founders found
#cp  AS_ALL_common_missingCC.nosex AS_ALL_samples.txt   #### AS_ALL_samples.txt comtains list of all samples



 /usr/lib64/libblas.so.3.0.3

 libf2c2 - Shared libraries for use with FORTRAN applications
libf2c2-dev - Development libraries for use with f2c


plink --bfile AS_ALL_common --allow-no-sex  --exclude  p6_remove_bad_snps.txt --remove related_samples.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out AS_ALL_common.noR.no6p.F
plink --bfile AS_ALL_common --allow-no-sex  --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out AS_ALL_common.no6p.F


## run smartpca on GA
nohup ~/bin/smartpca.perl -i AS_ALL_common.no6p.F.bed -a AS_ALL_common.no6p.F.bim -b AS_ALL_common.no6p.F.fam -o AS_ALL_common.no6p.F.pca -p AS_ALL_common.no6p.F.plot -e AS_ALL_common.no6p.F.eval -l AS_ALL_common.no6p.F.log

### snapmax03
nohup ~/bin/smartpca.perl -i AS_ALL_common.noR.no6p.F.bed -a AS_ALL_common.noR.no6p.F.bim -b AS_ALL_common.noR.no6p.F.fam -o AS_ALL_common.noR.no6p.F.pca -p AS_ALL_common.noR.no6p.F.plot -e AS_ALL_common.noR.no6p.F.eval -l AS_ALL_common.noR.no6p.F.log  > output.noR.no6p.F.txt

### snapmax03- use no strat filters
nohup ~/bin/smartpca.perl -i AS_ALL_common.no6p.F.bed -a AS_ALL_common.no6p.F.bim -b AS_ALL_common.no6p.F.fam -m 0 -o AS_ALL_common.no6p.F.ALL2.pca -p AS_ALL_common.no6p.F.ALL2.plot -e AS_ALL_common.no6p.F.ALL2.eval -l AS_ALL_common.no6p.F.ALL2.log  > output.no6p.F.ALL2.txt

nohup ~/bin/smartpca.perl -i AS_ALL_common.no6p.F.bed -a AS_ALL_common.no6p.F.bim -b AS_ALL_common.no6p.F.fam  -o AS_ALL_common.no6p.F.pca -p AS_ALL_common.no6p.F.plot -e AS_ALL_common.no6p.F.eval -l AS_ALL_common.no6p.F.log  > output.no6p.F.txt


 
#~/bin/smartpca.perl -i AS_ALL_common.no6p.F.bed -a AS_ALL_common.no6p.F.bim -b AS_ALL_common.no6p.F.fam -m 0 -o AS_ALL_common.no6p.F.ALL.pca -p AS_ALL_common.no6p.F.ALL.plot -e AS_ALL_common.no6p.F.ALL.eval -l AS_ALL_common.no6p.F.ALL.log  > output.no6p.F.ALL.txt




plink --bfile AS_ALL_common --allow-no-sex  --remove related_miss_samples.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out AS_ALL_common.noR.noHM.F
plink --bfile AS_ALL_common --allow-no-sex  --remove related_miss_samples.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out AS_ALL_common.noRHM.no6p.F
plink --bfile AS_ALL_common --allow-no-sex  --remove related_miss_samples.txt  --maf 0 --max-maf 1 --geno 1 --hwe 0 --mind 1 --make-bed --out AS_ALL_common.noRHM

## 301866 markers to
#convert f2c library to static   sudo ln -s /usr/lib/libf2c.a /usr/local/lib/libf2c.a
# then changed makile -lf2c to /usr/local/lib/libf2c.a seemed to work...

###################################RUNS WITH OR sample filters below ####################

## run smartpca
nohup ~/bin/smartpca.perl -i AS_ALL_common.noR.noHM.F.bed -a AS_ALL_common.noR.noHM.F.bim -b AS_ALL_common.noR.noHM.F.fam -o AS_ALL_common.noR.noHM.F.pca -p AS_ALL_common.noR.noHM.F.plot -e AS_ALL_common.noR.noHM.F.eval -l AS_ALL_common.noR.noHM.F.log  > output.noR.noHM.F.txt


# tions in effect:
#         --bfile AS_ALL_common
#         --allow-no-sex
#         --remove related_miss_samples.txt
#         --maf 0.01
#         --max-maf 1
#         --geno 0.05
#         --hwe 0.0000001
#         --mind 0.1
#         --make-bed
#         --out AS_ALL_common.noR.noHM.F
# 
# ** For gPLINK compatibility, do not use '.' in --out **
# Reading map (extended format) from [ AS_ALL_common.bim ] 
# 301866 markers to be included from [ AS_ALL_common.bim ]
# Reading pedigree information from [ AS_ALL_common.fam ] 
# 8028 individuals read from [ AS_ALL_common.fam ] 
# 8028 individuals with nonmissing phenotypes
# Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)
# Missing phenotype value is also -9
# 2181 cases, 5847 controls and 0 missing
# 0 males, 0 females, and 8028 of unspecified sex
# Warning, found 8028 individuals with ambiguous sex codes
# Writing list of these individuals to [ AS_ALL_common.noR.noHM.F.nosex ]
# Reading genotype bitfile from [ AS_ALL_common.bed ] 
# Detected that binary PED file is v1.00 SNP-major mode
# Reading individuals to remove [ related_miss_samples.txt ] ... 597 read
# 597 individuals removed with --remove option
# Before frequency and genotyping pruning, there are 301866 SNPs
# 7431 founders and 0 non-founders found
# 0 of 7431 individuals removed for low genotyping ( MIND > 0.1 )
# 11337 markers to be excluded based on HWE test ( p <= 1e-07 )
#         9014 markers failed HWE test in cases
#         11337 markers failed HWE test in controls
# Total genotyping rate in remaining individuals is 0.996664
# 2226 SNPs failed missingness test ( GENO > 0.05 )
# 135 SNPs failed frequency test ( MAF < 0.01 )
# After frequency and genotyping pruning, there are 288662 SNPs
# After filtering, 2090 cases, 5341 controls and 0 missing
# After filtering, 0 males, 0 females, and 7431 of unspecified sex
# Writing pedigree information to [ AS_ALL_common.noR.noHM.F.fam ] 
# Writing map (extended format) information to [ AS_ALL_common.noR.noHM.F.bim ] 
# Writing genotype bitfile to [ AS_ALL_common.noR.noHM.F.bed ] 
# Using (default) SNP-major mode
#
# Analysis finished: Mon Nov  3 17:33:55 2008




 plink --bfile AS_ALL_common --allow-no-sex  --remove related_samples.txt   --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --het   --out AS_ALL_common.noR.F.Het
 pleo@snpmax03:~/MergeAll_08_08$ head AS_ALL_common.noR.F.Het.het
              FID               IID       O(HOM)       E(HOM)        N(NM)            F
       WTCCC66069        WTCCC66069       186469    1.863e+05       288133     0.001907
       WTCCC66071        WTCCC66071       187583    1.864e+05       288369       0.0113
       WTCCC66073        WTCCC66073       186317    1.865e+05       288457    -0.001614
awk 'NR > 1 {print $0, ($4-$3)/$4}'  AS_ALL_common.noR.F.Het > temp    # From : het=1-Hom
 
 plink --bfile AS_ALL_common --allow-no-sex  --remove related_samples.txt   --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --missing   --out AS_ALL_common.noR.F.Miss
 plink --bfile AS_ALL_common --allow-no-sex  --remove related_miss_samples.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out AS_ALL_common.noRHM.no6p.F
 plink --bfile AS_ALL_common --allow-no-sex  --remove related_miss_samples.txt --exclude  bad.snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out AS_ALL_common.noRHM.with6p.F

plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_1.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out AS_ALL_common.keep1.no6p.F

Options in effect:
	--bfile AS_ALL_common
	--allow-no-sex
	--keep keep_num_1.txt
	--exclude p6_remove_bad_snps.txt
	--maf 0.01
	--max-maf 1
	--geno 0.05
	--hwe 0.0000001
	--mind 0.1
	--make-bed
	--out AS_ALL_common.keep1.no6p.F

** For gPLINK compatibility, do not use '.' in --out **
Reading map (extended format) from [ AS_ALL_common.bim ] 
301866 markers to be included from [ AS_ALL_common.bim ]
Reading pedigree information from [ AS_ALL_common.fam ] 
8028 individuals read from [ AS_ALL_common.fam ] 
8028 individuals with nonmissing phenotypes
Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)
Missing phenotype value is also -9
2181 cases, 5847 controls and 0 missing
0 males, 0 females, and 8028 of unspecified sex
Warning, found 8028 individuals with ambiguous sex codes
Writing list of these individuals to [ AS_ALL_common.keep1.no6p.F.nosex ]
Reading genotype bitfile from [ AS_ALL_common.bed ] 
Detected that binary PED file is v1.00 SNP-major mode
Reading list of SNPs to exclude [ p6_remove_bad_snps.txt ] ... 7914 read
Reading individuals to keep [ keep_num_1.txt ] ... 7195 read
833 individuals removed with --keep option
Before frequency and genotyping pruning, there are 293953 SNPs
7195 founders and 0 non-founders found
0 of 7195 individuals removed for low genotyping ( MIND > 0.1 )
9367 markers to be excluded based on HWE test ( p <= 1e-07 )
	8887 markers failed HWE test in cases
	9367 markers failed HWE test in controls
Total genotyping rate in remaining individuals is 0.996623
2237 SNPs failed missingness test ( GENO > 0.05 )
193 SNPs failed frequency test ( MAF < 0.01 )
After frequency and genotyping pruning, there are 282562 SNPs
After filtering, 2053 cases, 5142 controls and 0 missing
After filtering, 0 males, 0 females, and 7195 of unspecified sex
Writing pedigree information to [ AS_ALL_common.keep1.no6p.F.fam ] 
Writing map (extended format) information to [ AS_ALL_common.keep1.no6p.F.bim ] 
Writing genotype bitfile to [ AS_ALL_common.keep1.no6p.F.bed ] 
Using (default) SNP-major mode

Analysis finished: Thu Mar 26 20:40:52 2009

#plink --bfile AS_ALL_common --allow-no-sex  --remove related_miss_samples.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out AS_ALL_common.noRHM.no6p.F

pleo@snpmax03:~/MergeAll_08_08$ head  AS_ALL_common.noR.F.Miss.imiss    ### F_miss is the number wanted
              FID               IID MISS_PHENO   N_MISS   N_GENO   F_MISS
       WTCCC66069        WTCCC66069          N      468   301866  0.00155
       WTCCC66071        WTCCC66071          N      302   301866    0.001


#### run with keep1 +2 no need to prune the iterations
nohup ~/bin/smartpca.perl -i AS_ALL_common.keep1.no6p.F.bed -a AS_ALL_common.keep1.no6p.F.bim -b AS_ALL_common.keep1.no6p.F.fam -m 0 -o AS_ALL_common.keep1.no6p.F.pca -p AS_ALL_common.keep1.no6p.F.plot -e AS_ALL_common.keep1.no6p.F.eval -l AS_ALL_common.keep1.no6p.F.log  > output.noRHMkeep1.no6p.F.ALL.txt &


nohup ~/bin/smartpca.perl -i AS_ALL_common.noRHM.no6p.F.bed -a AS_ALL_common.noRHM.no6p.F.bim -b AS_ALL_common.noRHM.no6p.F.fam -o AS_ALL_common.noRHM.no6p.F.pca -p AS_ALL_common.noRHM.no6p.F.plot -e AS_ALL_common.noRHM.no6p.F.eval -l AS_ALL_common.noRHM.no6p.F.log  > output.noRHM.no6p.F.txt

nohup ~/bin/smartpca.perl -i AS_ALL_common.noRHM.no6p.F.bed -a AS_ALL_common.noRHM.no6p.F.bim -b AS_ALL_common.noRHM.no6p.F.fam -m 0 -o AS_ALL_common.noRHM.no6p.F.ALL.pca -p AS_ALL_common.noRHM.no6p.F.ALL.plot -e AS_ALL_common.noRHM.no6p.F.ALL.eval -l AS_ALL_common.noRHM.no6p.F.ALL.log  > output.noRHM.no6p.F.ALL.txt

nohup ~/bin/smartpca.perl -i AS_ALL_common.noRHM.with6p.F.bed -a AS_ALL_common.noRHM.with6p.F.bim -b AS_ALL_common.noRHM.with6p.F.fam -o AS_ALL_common.noRHM.with6p.F2.pca -p AS_ALL_common.noRHM.with6p.F2.plot -e AS_ALL_common.noRHM.with6p.F2.eval -l AS_ALL_common.noRHM.with6p.F2.log  > output.noRHM.with6p.F2.txt  &

plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_1.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --assoc --out AS_ALL_common.1


plink --bfile AS_ALL_common.noRHM.no6p.F --allow-no-sex  --keep keep_num_8.txt --maf 0 --max-maf 1 --geno 1 --hwe 0 --mind 0.1 --assoc --out AS_ALL_common.noRHM.no6p.F.4Meta
Options in effect:
	--bfile AS_ALL_common.noRHM.no6p.F
	--allow-no-sex
	--keep keep_num_8.txt
	--maf 0
	--max-maf 1
	--geno 1
	--hwe 0
	--mind 0.1
	--assoc
	--out AS_ALL_common.noRHM.no6p.F.4Meta

** For gPLINK compatibility, do not use '.' in --out **
Reading map (extended format) from [ AS_ALL_common.noRHM.no6p.F.bim ] 
280897 markers to be included from [ AS_ALL_common.noRHM.no6p.F.bim ]
Reading pedigree information from [ AS_ALL_common.noRHM.no6p.F.fam ] 
7430 individuals read from [ AS_ALL_common.noRHM.no6p.F.fam ] 
7430 individuals with nonmissing phenotypes
Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)
Missing phenotype value is also -9
2090 cases, 5340 controls and 0 missing
0 males, 0 females, and 7430 of unspecified sex
Warning, found 7430 individuals with ambiguous sex codes
Writing list of these individuals to [ AS_ALL_common.noRHM.no6p.F.4Meta.nosex ]
Reading genotype bitfile from [ AS_ALL_common.noRHM.no6p.F.bed ] 
Detected that binary PED file is v1.00 SNP-major mode
Reading individuals to keep [ keep_num_8.txt ] ... 5790 read
1640 individuals removed with --keep option
Before frequency and genotyping pruning, there are 280897 SNPs
5790 founders and 0 non-founders found
0 of 5790 individuals removed for low genotyping ( MIND > 0.1 )
0 markers to be excluded based on HWE test ( p <= 0 )
	0 markers failed HWE test in cases
	0 markers failed HWE test in controls
Total genotyping rate in remaining individuals is 0.997411
0 SNPs failed missingness test ( GENO > 1 )
0 SNPs failed frequency test ( MAF < 0 )
After frequency and genotyping pruning, there are 280897 SNPs
After filtering, 1905 cases, 3885 controls and 0 missing
After filtering, 0 males, 0 females, and 5790 of unspecified sex
Writing main association results to [ AS_ALL_common.noRHM.no6p.F.4Meta.assoc ] 






plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_1.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --model --out AS_ALL_common.1
plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_2.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --model --out AS_ALL_common.2
plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_3.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --model --out AS_ALL_common.3
plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_4.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --model --out AS_ALL_common.4
plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_5.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --model --out AS_ALL_common.5
plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_6.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --model --out AS_ALL_common.6
plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_7.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --model --out AS_ALL_common.7
plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_8.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --model --out AS_ALL_common.8

plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_9.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --model --out AS_ALL_common.9

plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_10.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --model --out AS_ALL_common.10


 plink --bfile AS_ALL_common   --keep test_genome.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1  --recode



 plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_1.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out AS_ALL_common.keep1

 plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_8.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out AS_ALL_common.keep8

#/home/pleo/plink-1.02-x86_64/plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_8.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --model --out AS_ALL_common.8.2

plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_8.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --model --out AS_ALL_common.8.final

plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_8.txt --exclude  p6_remove.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --model --out AS_ALL_common.8.no6p

plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_8.txt --exclude  p6_remove_bad_snps_IL23R_ARTS.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --model --out AS_ALL_common.8.final.IL23.ARTS

plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_8.txt  --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --model --out AS_ALL_common.8.ALL


################################ Control vs Control to filter SNPs
  cp AS_ALL_common.fam AS_ALL_common.TRUE.fam
   cp AS_ALL_common_wtccc_cases.fam AS_ALL_common.fam

nohup plink --bfile AS_ALL_common --allow-no-sex  --keep keep_8_nocases.txt --exclude  p6_remove_bad_snps.txt  --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --model --out AS_ALL_common.8.controlUKvUS

 > mean(assoc.file[assoc.file[,"TEST"]=="ALLELIC","CHISQ"],na.rm=TRUE)
[1] 1.404701

 cp AS_ALL_common.TRUE.fam  AS_ALL_common.fam 

plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_8.txt --exclude  p6_remove_bad_control_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --assoc --out AS_ALL_common.8.snpF




####################################### eigenstrat run
nohup ~/bin/smartpca.perl -i AS_ALL_common.noRHM.no6p.F.bed -a AS_ALL_common.noRHM.no6p.F.bim -b AS_ALL_common.noRHM.no6p.F.ind -m 0 -o AS_ALL_common.noRHM.no6p.F.ALL.pca -p AS_ALL_common.noRHM.no6p.F.ALL.plot -e AS_ALL_common.noRHM.no6p.F.ALL.eval -l AS_ALL_common.noRHM.no6p.F.ALL.log  > output.noRHM.no6p.F.ALL.txt

nohup ~/bin/smartpca.perl -i AS_ALL_common.noRHM.no6p.F.bed -a AS_ALL_common.noRHM.no6p.F.bim -b AS_ALL_common.noRHM.no6p.F.fam -o AS_ALL_common.noRHM.no6p.F2.pca -p AS_ALL_common.noRHM.no6p.F.plot -e AS_ALL_common.noRHM.no6p.F.eval -l AS_ALL_common.noRHM.no6p.F.log  > output.noRHM.no6p.F.txt

nohup ~/bin/smartpca.perl -i AS_ALL_common.noRHM.no6p.F.bed -a AS_ALL_common.noRHM.no6p.F.bim -b AS_ALL_common.noRHM.no6p.F.fam -o AS_ALL_common.noRHM.no6p.F3.pca -p AS_ALL_common.noRHM.no6p.F3.plot -e AS_ALL_common.noRHM.no6p.F3.eval -s 3 -l AS_ALL_common.noRHM.no6p.F3.log  > output.noRHM.no6p.F3.txt

 plink --bfile AS_ALL_common.noRHM.no6p.F --assoc --out  AS_ALL_common.noRHM.no6p.F

../bin/convertf -p ../CONVERTF/par.BED.EIGENSTRAT3    # convert types this way does not work  ls *
#../bin/convertf -p ../CONVERTF/par.PACKEDPED.PACKEDPED   # convert types this way does not work
../CONVERTF/ind2pheno.perl AS_ALL_common.noRHM.no6p.F.ind AS_ALL_common.noRHM.no6p.F.pheno
evec2pca.perl 10  AS_ALL_common.noRHM.no6p.F.pca.evec AS_ALL_common.noRHM.no6p.F.fam AS_ALL_common.noRHM.no6p.F.pca
evec2pca_with_FAM.perl 10  AS_ALL_common.noRHM.no6p.F.pca.evec AS_ALL_common.noRHM.no6p.F.fam AS_ALL_common.noRHM.no6p.F.pca ##to fix name:name problem in the fam file

../bin/eigenstrat.big.perl -i  AS_ALL_common.noRHM.no6p.F.eigenstratgeno -j  AS_ALL_common.noRHM.no6p.F.pheno -p AS_ALL_common.noRHM.no6p.F.pca -o AS_ALL_common.noRHM.no6p.F2.chisq

eigenstrat -i  AS_ALL_common.noRHM.no6p.F.eigenstratgeno -j  AS_ALL_common.noRHM.no6p.F.pheno -p AS_ALL_common.noRHM.no6p.F.pca -o AS_ALL_common.noRHM.no6p.F.chisq

plink --bfile AS_ALL_common --allow-no-sex --keep keep_wtccc_cases_RHM.txt --exclude  p6_remove_bad_control_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --assoc --out AS_ALL_common_WTCCC_CASE


plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_1.txt --exclude  p6_remove_bad_control_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out AS_ALL_common.1.no6pRHM
plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_8.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out AS_ALL_common.8.no6pRHM
plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_1.txt --chr 1 --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out AS_ALL_common_1_chr1_no6pRHM
plink --bfile AS_ALL_common_1_chr1_no6pRHM --exclude  p6_remove_bad_control_snps.txt --make-bed --out AS_ALL_common_1_chr1_C_no6pRHM
######################################################################################################
######################################################################################################
############################################################################################################################################################################################################
######################################################################################################

################################### Shared controls

##in /media/Bioinform-D/Research/Matt Brown/EIGENSOFT/Merge ALL + HapMap Aug 2008/FINAL BED files
 /media/Bioinform-D/Research/Matt Brown/EIGENSOFT/Merge ALL + HapMap Aug 2008/FINAL BED files

######## files built using:
############ filtered
plink --bfile AS_ALL_common --allow-no-sex  --remove related_miss_samples.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out AS_ALL_common.noRHM.no6p.F
# 2053 cases, 5142 controls there are 290361 SNPs

############ startification corrected with and without LD regions
########################## LD trimmed
plink --bfile AS_ALL_common --allow-no-sex  --keep eigenstrat_keep_samples.txt --extract ld.trimmed.snp.list.txt --make-bed --out AS_ALL_common.EIGEN.noLD.F

#############
plink --bfile AS_ALL_common --allow-no-sex  --keep eigenstrat_keep_samples.txt --extract full.snp.list.txt --make-bed --out AS_ALL_common.FULL.noLD.F

############## controls given to PEX
plink --bfile AS_ALL_common.FULL.noLD.F --filter-controls --make-bed --out  ALL_common.FULL.noLD.F


######################################################################################################
######################################################################################################
############################################################################################################################################################################################################
######################################################################################################

################################### MONEY
################################### MONEY
################################### MONEY
################################### MONEY
################################### MONEY

#################  EIGENSTRAT RUN WITH 6p
#### first build genotype file with required samples

######################################################################################################
######################################################################################################
############################################################################################################################################################################################################
######################################################################################################

### eigenstart 3 max need to use:
/bin/smartpca.perl -i main_data_02.eigenstratgeno -a main_data_02.bim -b
main_data_02.ind -k 10 -o main_data_02.pca -p main_data_02.plot -e
main_data_02.eval -l main_data_02.log
 
################################### MONEY
################################### MONEY
################################### MONEY
################################### MONEY
################################### MONEY

#################  EIGENSTRAT RUN WITH 6p
#### first build genotype file with required samples 
###eigenstrat_keep_samples.txt same as keep1 +2 
plink --bfile AS_ALL_common --allow-no-sex  --keep eigenstrat_keep_samples.txt --exclude  bad.snps.more.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --make-bed --out AS_ALL_common.EIGEN.with6p.F
# 2053 cases, 5142 controls there are 290361 SNPs

########################## LD trimmed
plink --bfile AS_ALL_common --allow-no-sex  --keep eigenstrat_keep_samples.txt --extract ld.trimmed.snp.list.txt --make-bed --out AS_ALL_common.EIGEN.noLD.F

nohup ~/bin/smartpca.perl -i AS_ALL_common.EIGEN.noLD.F.bed -a AS_ALL_common.EIGEN.noLD.F.bim -b AS_ALL_common.EIGEN.noLD.F.fam -o AS_ALL_common.EIGEN.noLD.F.pca -p AS_ALL_common.EIGEN.noLD.F.plot -e AS_ALL_common.EIGEN.noLD.F.eval -l AS_ALL_common.EIGEN.noLD.F.log  > output.noLD.txt &

########### just using 3 sigma :

nohup ~/bin/smartpca.perl -i AS_ALL_common.EIGEN.noLD.F.bed -a AS_ALL_common.EIGEN.noLD.F.bim -b AS_ALL_common.EIGEN.noLD.F.fam  -o AS_ALL_common.EIGEN.noLD.F.sig3.pca -p AS_ALL_common.EIGEN.noLD.F.sig3.plot -e AS_ALL_common.EIGEN.noLD.F.sig3.eval -l AS_ALL_common.EIGEN.noLD.F.sig3.log -s 3  > output.noLD.sig3.txt &


../bin/convertf -p ../CONVERTF/par.BED.EIGENSTRAT10
../CONVERTF/ind2pheno.perl AS_ALL_common.EIGEN.noLD.F.ind AS_ALL_common.EIGEN.noLD.F.pheno
evec2pca.perl 10  AS_ALL_common.EIGEN.noLD.F.pca.evec AS_ALL_common.EIGEN.noLD.F.fam AS_ALL_common.EIGEN.noLD.F.pca


nohup ../bin/eigenstrat.big.perl -i  AS_ALL_common.EIGEN.noLD.F.eigenstratgeno -j   AS_ALL_common.EIGEN.noLD.F.pheno -p   AS_ALL_common.EIGEN.noLD.F.pca -o  AS_ALL_common.EIGEN.noLD.F.eig10.chisq &

nohup ../bin/eigenstrat.big.perl -i  AS_ALL_common.EIGEN.noLD.F.eigenstratgeno -j   AS_ALL_common.EIGEN.noLD.F.pheno -p   AS_ALL_common.EIGEN.noLD.F.pca -l 2  -o  AS_ALL_common.EIGEN.noLD.F.eig2.chisq &

nohup ../bin/eigenstrat.big.perl -i  AS_ALL_common.EIGEN.noLD.F.eigenstratgeno -j   AS_ALL_common.EIGEN.noLD.F.pheno -p   AS_ALL_common.EIGEN.noLD.F.pca -l 3  -o  AS_ALL_common.EIGEN.noLD.F.eig3.chisq &

nohup ../bin/eigenstrat.big.perl -i  AS_ALL_common.EIGEN.noLD.F.eigenstratgeno -j   AS_ALL_common.EIGEN.noLD.F.pheno -p   AS_ALL_common.EIGEN.noLD.F.pca -l 4  -o  AS_ALL_common.EIGEN.noLD.F.eig4.chisq &

#############  USE ALL SNPSs
plink --bfile AS_ALL_common --allow-no-sex  --keep eigenstrat_keep_samples.txt --extract full.snp.list.txt --make-bed --out AS_ALL_common.FULL.noLD.F

../bin/convertf -p ../CONVERTF/par.BED.EIGENSTRAT11
../CONVERTF/ind2pheno.perl AS_ALL_common.FULL.noLD.F.ind AS_ALL_common.FULL.noLD.F.pheno
evec2pca.perl 10  AS_ALL_common.EIGEN.noLD.F.pca.evec AS_ALL_common.FULL.noLD.F.fam AS_ALL_common.FULL.noLD.F.pca

nohup ../bin/eigenstrat.big.perl -i  AS_ALL_common.FULL.noLD.F.eigenstratgeno -j   AS_ALL_common.FULL.noLD.F.pheno -p   AS_ALL_common.FULL.noLD.F.pca -o  AS_ALL_common.FULL.noLD.F.eig10.chisq &

nohup ../bin/eigenstrat.big.perl -i  AS_ALL_common.FULL.noLD.F.eigenstratgeno -j   AS_ALL_common.FULL.noLD.F.pheno -p   AS_ALL_common.FULL.noLD.F.pca -l 2  -o  AS_ALL_common.FULL.noLD.F.eig2.chisq &

nohup ../bin/eigenstrat.big.perl -i  AS_ALL_common.FULL.noLD.F.eigenstratgeno -j   AS_ALL_common.FULL.noLD.F.pheno -p   AS_ALL_common.FULL.noLD.F.pca -l 3  -o  AS_ALL_common.FULL.noLD.F.eig3.chisq &

nohup ../bin/eigenstrat.big.perl -i  AS_ALL_common.FULL.noLD.F.eigenstratgeno -j   AS_ALL_common.FULL.noLD.F.pheno -p   AS_ALL_common.FULL.noLD.F.pca -l 4  -o  AS_ALL_common.FULL.noLD.F.eig4.chisq &



#################### run equivalent plink run ### note change name of file else .pheno file causes default to qualitative traithead 
plink --bfile AS_ALL_common --allow-no-sex  --keep keep_num_8.txt --extract full.snp.list.txt  --assoc --out AS_ALL_common.8.BEST

# plink --bfile AS_ALL_common.EIGEN.with6p.F --allow-no-sex  --keep keep_num_8.txt --exclude  p6_remove_bad_snps.txt --maf 0.01 --max-maf 1 --geno 0.05 --hwe 0.0000001 --mind 0.1 --assoc --out AS_ALL_common.8.BESTe


../bin/convertf -p ../CONVERTF/par.BED.EIGENSTRAT8    # convert types this way does not work  ls *
#../bin/convertf -p ../CONVERTF/par.PACKEDPED.PACKEDPED   # convert types this way does not work
../CONVERTF/ind2pheno.perl AS_ALL_common.EIGEN.with6p.F.ind AS_ALL_common.EIGEN.with6p.F.pheno
### get new pca file have to reduce the size of  AS_ALL_common.noRHM.no6p.F.pca it contains ALL samples
evec2pca.perl 10  AS_ALL_common.noRHM.no6p.F.pca.evec AS_ALL_common.EIGEN.with6p.F.fam AS_ALL_common.EIGEN.with6p.F.pca

## run with new PCA file ############# try differenet eigenvectors
nohup ../bin/eigenstrat.big.perl -i AS_ALL_common.EIGEN.with6p.F.eigenstratgeno -j  AS_ALL_common.EIGEN.with6p.F.pheno -p  AS_ALL_common.EIGEN.with6p.F.pca -o AS_ALL_common.EIGEN.with6p.F.eig10.chisq &

nohup ../bin/eigenstrat.big.perl -i AS_ALL_common.EIGEN.with6p.F.eigenstratgeno -j  AS_ALL_common.EIGEN.with6p.F.pheno -p  AS_ALL_common.EIGEN.with6p.F.pca -l 1  -o AS_ALL_common.EIGEN.with6p.F.eig1.chisq &

nohup ../bin/eigenstrat.big.perl -i AS_ALL_common.EIGEN.with6p.F.eigenstratgeno -j  AS_ALL_common.EIGEN.with6p.F.pheno -p  AS_ALL_common.EIGEN.with6p.F.pca -l 2  -o AS_ALL_common.EIGEN.with6p.F.eig2.chisq &




##################   noibdnoIBD  sampless
plink --bfile AS_ALL_common --allow-no-sex  --keep no.ibd.keep_samples.txt --extract full.snp.list.txt --make-bed --out AS_ALL_common.noIBD.with6p.F

../bin/convertf -p ../CONVERTF/par.BED.EIGENSTRAT9    # convert types this way does not work  ls *
../CONVERTF/ind2pheno.perl AS_ALL_common.noIBD.with6p.F.ind AS_ALL_common.noIBD.with6p.F.pheno
evec2pca.perl 10  AS_ALL_common.EIGEN.noLD.F.pca.evec AS_ALL_common.noIBD.with6p.F.fam AS_ALL_common.noIBD.with6p.F.pca

## run with original PCA file
../bin/eigenstrat.big.perl -i AS_ALL_common.noIBD.with6p.F.eigenstratgeno -j  AS_ALL_common.noIBD.with6p.F.pheno -p AS_ALL_common.noIBD.with6p.F.pca  -o AS_ALL_common.noIBD.with6p.F.chisq

../bin/eigenstrat.big.perl -i AS_ALL_common.noIBD.with6p.F.eigenstratgeno -j  AS_ALL_common.noIBD.with6p.F.pheno -p AS_ALL_common.noIBD.with6p.F.pca -l 4 -o AS_ALL_common.noIBD.with6p.F.eig4.chisq

../bin/eigenstrat.big.perl -i AS_ALL_common.noIBD.with6p.F.eigenstratgeno -j  AS_ALL_common.noIBD.with6p.F.pheno -p AS_ALL_common.noIBD.with6p.F.pca -l 4 -o AS_ALL_common.noIBD.with6p.F.eig2.chisq
######################## Get EIGENSTRAT working ### SOLN make sure that pca file does not contain zero

nohup ~/bin/smartpca.perl -i AS_ALL_common.1.no6pRHM.bed -a AS_ALL_common.1.no6pRHM.bim -b AS_ALL_common.1.no6pRHM.fam  -o AS_ALL_common.1.no6pRHM.pca -p AS_ALL_common.1.no6pRHM.plot -e AS_ALL_common.1.no6pRHM.eval -l AS_ALL_common.1.no6pRHM.log -m 0 > output.AS_ALL_common.1.no6pRHM.txt
nohup ~/bin/smartpca.perl -i AS_ALL_common.1.no6pRHM.bed -a AS_ALL_common.1.no6pRHM.bim -b AS_ALL_common.1.no6pRHM.fam  -o AS_ALL_common.1.no6pRHM2.pca -p AS_ALL_common.1.no6pRHM2.plot -e AS_ALL_common.1.no6pRHM2.eval -l AS_ALL_common.1.no6pRHM2.log > output.AS_ALL_common.1.no6pRHM2.txt
nohup ~/bin/smartpca.perl -i AS_ALL_common.8.no6pRHM.bed -a AS_ALL_common.8.no6pRHM.bim -b AS_ALL_common.8.no6pRHM.fam  -o AS_ALL_common.8.no6pRHM.pca -p AS_ALL_common.8.no6pRHM.plot -e AS_ALL_common.8.no6pRHM.eval -l AS_ALL_common.8.no6pRHM.log > output.AS_ALL_common.8.no6pRHM.txt &
nohup ~/bin/smartpca.perl -i AS_ALL_common.8.no6pRHM.bed -a AS_ALL_common.8.no6pRHM.pedsnp -b AS_ALL_common.8.no6pRHM.pedind  -o AS_ALL_common.8.no6pRHM2.pca -p AS_ALL_common.8.no6pRHM2.plot -e AS_ALL_common.8.no6pRHM2.eval -l AS_ALL_common.8.no6pRHM2.log > output.AS_ALL_common.8.no6pRHM2.txt &
nohup ~/bin/smartpca.perl -i AS_ALL_common.8.no6pRHM.bed -a AS_ALL_common_8_no6pRHM.snp -b AS_ALL_common_8_no6pRHM.ind  -o AS_ALL_common.8.no6pRHM3.pca -p AS_ALL_common.8.no6pRHM3.plot -e AS_ALL_common.8.no6pRHM3.eval -m 0 -l AS_ALL_common.8.no6pRHM3.log > output.AS_ALL_common.8.no6pRHM3.txt &

../bin/convertf -p ../CONVERTF/par.BED.EIGENSTRAT7

../bin/convertf -p ../CONVERTF/par.BED.EIGENSTRAT4    # convert types this way does not work  ls *

mv  AS_ALL_common.1.no6pRHM.eigenstratgeno  AS_ALL_common_1_no6pRHM.eigenstratgeno
mv AS_ALL_common.1.no6pRHM.snp AS_ALL_common_1_no6pRHM.snp
mv  AS_ALL_common.1.no6pRHM.ind  AS_ALL_common_1_no6pRHM.ind
mv AS_ALL_common.1.no6pRHM.pca AS_ALL_common_1_no6pRHM.pca

../CONVERTF/ind2pheno.perl AS_ALL_common_1_no6pRHM.ind AS_ALL_common_1_no6pRHM.pheno
eigenstrat -i  AS_ALL_common_1_no6pRHM.eigenstratgeno -j  AS_ALL_common_1_no6pRHM.pheno -p  AS_ALL_common_1_no6pRHM.pca -o AS_ALL_common_1_no6pRHM.chisq


mv
../bin/eigenstrat.big.perl -i  AS_ALL_common_1_no6pRHM.eigenstratgeno -j  AS_ALL_common.1.no6pRHM.pheno -p AS_ALL_common.1.no6pRHM.pca -o AS_ALL_common.1.no6pRHM.chisq



 plink --bfile AS_ALL_common.1.no6pRHM --allow-no-sex   --maf 1 --max-maf 1 --geno 1 --hwe 0 --mind 1 --assoc


../bin/convertf -p ../CONVERTF/par.BED.EIGENSTRAT5    # convert types this way does not work  ls *
../CONVERTF/ind2pheno.perl merge_ALL_hapmap_no6pF.ind merge_ALL_hapmap_no6pF.pheno
eigenstrat -i merge_ALL_hapmap_no6pF.eigenstratgeno -j merge_ALL_hapmap_no6pF.pheno -p merge_ALL_hapmap_no6pF.pca -o merge_ALL_hapmap_no6pF.chisq


../bin/convertf -p ../CONVERTF/par.BED.EIGENSTRAT6    # convert types this way does not work  ls *
../CONVERTF/ind2pheno.perl AS_ALL_common_1_chr1_C_no6pRHM.ind AS_ALL_common_1_chr1_C_no6pRHM.pheno
eigenstrat -i  AS_ALL_common_1_chr1_C_no6pRHM.eigenstratgeno -j  AS_ALL_common_1_chr1_C_no6pRHM.pheno -p  AS_ALL_common_1_no6pRHM.pca -o AS_ALL_common_1_chr1_C_no6pRHM.chisq


pleo@snpmax03:~/MergeAll_08_08$ wc -l  merge_ALL_hapmap_no6pF.eigenstratgeno
277909 merge_ALL_hapmap_no6pF.eigenstratgeno
pleo@snpmax03:~/MergeAll_08_08$ wc -L  merge_ALL_hapmap_no6pF.eigenstratgeno
6051 merge_ALL_hapmap_no6pF.eigenstratgeno
pleo@snpmax03:~/MergeAll_08_08$ wc -l  AS_ALL_common.1.no6pRHM.eigenstratgeno
266125 AS_ALL_common.1.no6pRHM.eigenstratgeno
pleo@snpmax03:~/MergeAll_08_08$ wc -L  AS_ALL_common.1.no6pRHM.eigenstratgeno
7193 AS_ALL_common.1.no6pRHM.eigenstratgeno


 

evec2pca.perl 10  merge_ALL_hapmap_no6pF.pca.evec merge_ALL_hapmap_no6pF.fam test.pca
evec2pca.perl 10 AS_ALL_common.8.no6pRHM.pca.evec AS_ALL_common.8.no6pRHM.fam test2.pca
############################################################
############################################################
############################################################

plink --bfile AS_ALL_common --allow-no-sex --keep wtccc_good.txt --extract UK_confirmation_set_SNPs_ALL.txt --maf 0 --max-maf 1 --geno 1 --hwe 0 --mind 1  --make-bed --out WTCCC_conf

plink --bfile UK_conf_cases --bmerge WTCCC_conf.bed  WTCCC_conf.bim  WTCCC_conf.fam  --allow-no-sex --make-bed --out UK_WTCCC_conf

plink --bfile UK_WTCCC_conf --allow-no-sex    --assoc 
plink --bfile UK_WTCCC_conf --allow-no-sex    --model

   1   RS11209026    A    G    TREND        99/2069       188/2680        9.166    1     0.002466
   1   RS11209026    A    G  ALLELIC        99/2069       188/2680        9.086    1     0.002576


plink --bfile UK_conf_cases --recode --out UK_conf_cases_ATGC
plink --bfile US_conf_cases --recode --out US_conf_cases_ATGC



plink --bfile UK_conf_cases --flip  UK_conf_flips.txt --make-bed --allow-no-sex --out UK_conf_cases_f
plink --bfile US_conf_cases --flip  US_conf_flips.txt --make-bed --allow-no-sex --out US_conf_cases_f


plink --bfile US_conf_cases_f --bmerge WTCCC_conf.bed  WTCCC_conf.bim  WTCCC_conf.fam  --allow-no-sex --make-bed --out US_WTCCC_conf
plink --bfile UK_conf_cases_f --bmerge WTCCC_conf.bed  WTCCC_conf.bim  WTCCC_conf.fam  --allow-no-sex --make-bed --out UK_WTCCC_conf

plink --bfile UK_conf_cases --recode --out UK_conf_cases_ATGC  # now on correct strand.
plink --bfile US_conf_cases --recode --out US_conf_cases_ATGC

plink --bfile US_WTCCC_conf --allow-no-sex    --assoc


plink --bfile US_conf_cases_f --allow-no-sex    --model
   1   RS11209026    A    G    TREND        94/1988       180/3204        1.685    1       0.1943
   1   RS11209026    A    G  ALLELIC        94/1988       180/3204        1.751    1       0.1858
   1   RS11209026    A    G      DOM         91/950       171/1521           NA   NA           NA
   1   RS11209026    A    G      REC         3/1038         9/1683           NA   NA           NA
plink --bfile US_conf_cases_f --remove US.conf.cut.txt --allow-no-sex    --model
   1   RS11209026    A    G    TREND        91/1803       171/3003        0.795    1       0.3726
   1   RS11209026    A    G  ALLELIC        91/1803       171/3003        0.822    1       0.3646
   1   RS11209026    A    G      DOM         88/859       163/1424           NA   NA           NA
   1   RS11209026    A    G      REC          3/944         8/1579           NA   NA           NA
plink --bfile US_conf_cases_f --remove US.conf.cut.txt --exclude ethnicity_snps_avaliable.txt --allow-no-sex  --model   


   

plink --bfile US_conf_cases_f --filter-cases   --make-bed --out US_conf_casesONLY_f
plink --bfile US_conf_cases_f --filter-controls   --make-bed --out US_conf_controlONLY_f

plink --bfile US_conf_casesONLY_f --bmerge WTCCC_conf.bed  WTCCC_conf.bim  WTCCC_conf.fam  --allow-no-sex --make-bed --out US_only_WTCCC_conf
plink --bfile US_only_WTCCC_conf --allow-no-sex    --model
 1   RS11209026    A    G    TREND        94/1988       188/2680         9.21    1     0.002407
   1   RS11209026    A    G  ALLELIC        94/1988       188/2680        9.346    1     0.002235
   1   RS11209026    A    G      DOM         91/950       181/1253           NA   NA           NA
   1   RS11209026    A    G      REC         3/1038         7/1427           NA   NA           NA
plink --bfile US_only_WTCCC_conf  --remove US.conf.cut.txt --allow-no-sex    --model
  1   RS11209026    A    G     GENO         0/1/59         1/7/67           NA   NA           NA
   1   RS11209026    A    G    TREND          1/119          9/141        4.268    1      0.03885
   1   RS11209026    A    G  ALLELIC          1/119          9/141         4.99    1       0.0255
   
plink --bfile US_conf_cases_f --bmerge UK_conf_cases_f.bed  UK_conf_cases_f.bim  UK_conf_cases_f.fam  --allow-no-sex --make-bed --out US_UK_conf

plink --bfile US_UK_conf --remove conf.cut.txt --exclude ethnicity_snps_avaliable.txt  --allow-no-sex    --assoc
 1   RS11209026    A    G    TREND       191/3921       171/3063        1.577    1       0.2092
 1   RS11209026    A    G  ALLELIC       191/3921       171/3063        1.596    1       0.2065

plink --bfile US_conf_casesONLY_f --bmerge UK_conf_cases_f.bed  UK_conf_cases_f.bim  UK_conf_cases_f.fam  --allow-no-sex --make-bed --out US_UK_conf_cases
plink --bfile US_UK_conf_cases --bmerge WTCCC_conf.bed  WTCCC_conf.bim  WTCCC_conf.fam  --allow-no-sex --make-bed --out US_UK_cases_WTCCC
plink --bfile US_UK_cases_WTCCC --remove conf.cut.txt --allow-no-sex   --exclude ethnicity_snps_avaliable.txt  --assoc
 1   RS11209026    A    G    TREND       191/3921       188/2680        12.02    1     0.000525
 1   RS11209026    A    G  ALLELIC       191/3921       188/2680        12.01    1    0.0005305




plink --bfile  AS_ALL_common.noRHM --extract ethnicity_snps_avaliable.txt --allow-no-sex --make-bed --out All_samples_ethnicity
plink --bfile  US_conf_cases_f --extract ethnicity_snps_avaliable.txt --allow-no-sex --make-bed --out US_conf_cases_f_eth
plink --bfile  UK_conf_cases_f --extract ethnicity_snps_avaliable.txt --allow-no-sex --make-bed --out UK_conf_cases_f_eth

 plink --bfile  AS_ALL_common.noRHM --extract ethnicity_snps_avaliable.txt  --out All_ethnicity_snps_avaliable
 plink --bfile AS_ALL_common.noRHM  --flip to_flip_300v1.txt --make-bed --allow-no-sex --out Controls_300v1_07_f

 plink --merge-list ethnicity_list.txt --allow-no-sex --make-bed --out ALL_eth_samples_UK_US
 plink --bfile  AS_ALL_common.noRHM --extract ethnicity_snps_avaliable.txt --allow-no-sex --make-bed --out All_samples_ethnicity

nohup ~/bin/smartpca.perl -i ALL_samples_UK_US.bed -a ALL_samples_UK_US.bim -b ALL_samples_UK_US.fam -m 0 -o ALL_samples_UK_US.pca -p ALL_samples_UK_US.plot -e ALL_samples_UK_US.eval -l ALL_samples_UK_US.log  > out_ALL_samples_UK_US.txt

nohup ~/bin/smartpca.perl -i ALL_eth_samples_UK_US.bed -a ALL_eth_samples_UK_US.bim -b ALL_eth_samples_UK_US.fam -m 0 -o ALL_eth_samples_UK_US.pca -p ALL_eth_samples_UK_US.plot -e ALL_eth_samples_UK_US.eval -l ALL_eth_samples_UK_US.log  > out_ALL_eth_samples_UK_US.txt   &


































Miptables -I RH-Firewall-1-INPUT -m state --state NEW -p tcp --destination-port 5901 -j ACCEPT
use 130.102.227.233



you have already excluded all the individuals you want (with the per-individual genotyping threshold option), then setting

     --mind 1 

will skip the step where per-individual genotyping rates are calculated, which can reduce the time
taken to load the file. Note, the command --all is equivalent to specifying --mind 1 --geno 1 --maf 0 (i.e. do not apply any filters).


 
 
 ##flip
  plink --bfile US_control_550a300 --flip mergeUS_550a300.missnp --make-bed --allow-no-sex --out US_control_550a300pf
 ##merge
   plink --bfile US_combined --bmerge US_control_550a300.bed US_control_550a300.bim US_control_550a300.fam --allow-no-sex --make-bed --out mergeUS_550a300
 ##run
 plink --bfile mergeUS_550a300pf --allow-no-sex --assoc --maf0.01 --max-maf 1 --geno 0.10 --hwe 0.00001 --mind 0.1


   plink --mergeUS_300pf --assoc -perm -p2

   plink --bfile mergeALL --exclude to_exclude_SNPs.txt --make-bed --allow-no-sex --out mergeALLcSNPs
   plink --bfile mergeALLcSNPs --allow-no-sex --genome --out mergeALLcSNPs_genome &
   plink --bfile mergeALLcSNPs --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00001 --mind 0.10 --out mergeALLcSNPsAssoc

    #### process and include HapMap data
    plink --bfile hapmap --exclude to_exclude_SNPs.txt  --make-bed --allow-no-sex --out hapmap_cSNPs
    plink --bfile mergeALLcSNPs --bmerge hapmap_cSNPs.bed hapmUSap_cSNPs.bim hapmap_cSNPs.fam --allow-no-sex --make-bed --out mergeALLcSNPsHapMap
#  Found 153495 SNPs that do not match in terms of strand
# (nb. flipped A/T and C/G SNPs will of course go undetected here)
# Writing identifiable strand problem SNPs to [ mergeALLcSNPsHapMap.missnp ]

    plink --bfile hapmap_cSNPs --flip mergeALLcSNPsHapMap.missnp --make-bed --allow-no-sex --out hapmap_cSNpsf

     plink --bfile mergeALLcSNPs --bmerge hapmap_cSNPsf.bed hapmap_cSNPsf.bim hapmap_cSNPsf.fam --allow-no-sex --make-bed --out mergeALLcSNPsHapMapf
      plink --bfile mergeALLcSNPsHapMapf --allow-no-sex --genome --out mergeALLcSNPsHapMap_G &

     ./../bin/smartpca.perl -i hapmap_cSNPsf.bed -a hapmap_cSNPsf.bim -b hapmap_cSNPsf.fam -o merge.pca -p merge.plot -e merge.eval -l merge.log

 ./../bin/convertf -p ../CONVERTF/par.BED.EIGENSTRAT   # run convert on a bed file

./../bin/smartpca.perl -i US_control_550v3_550v1_300_cSNPs.eigenstratgeno -a US_control_550v3_550v1_300_cSNPs.snp -b US_control_550v3_550v1_300_cSNPs.ind -o US_control_550v3_550v1_300_cSNPs.pca -p US_control_550v3_550v1_300_cSNPs.plot -e US_control_550v3_550v1_300_cSNPs.eval -l US_control_550v3_550v1_300_cSNPs.log   &

 ./plink --bfile mergeUS_550v3_550v1_300 --exclude to_exclude_SNPs.txt --make-bed --allow-no-sex --out mergeUS_control_cSNPs
 ./plink --bfile mergeUS --exclude to_exclude_SNPs.txt --make-bed --allow-no-sex --out mergeUS_WTCCC_cSNPs
./plink --bfile mergeUS_WTCCC_cSNPs --freq

### STrand analysis:
Recombine icontrol DB data without strand flips:
  ./plink --bfile US_control_550v3 --exclude to_exclude_SNPs.txt --make-bed --allow-no-sex --out US_control_550v3_cSNPs
  ./plink --bfile US_control_550v1 --exclude to_exclude_SNPs.txt --make-bed --allow-no-sex --out US_control_550v1_cSNPs
  ./plink --bfile US_control_300 --exclude to_exclude_SNPs.txt --make-bed --allow-no-sex --out US_control_300_cSNPs

   ./plink --bfile US_control_550v3_cSNPs --bmerge US_control_550v1_cSNPs.bed US_control_550v1_cSNPs.bim US_control_550v1_cSNPs.fam --allow-no-sex --make-bed --out US_control_550v3_550v1_cSNPs
    ./plink --bfile US_control_550v3_550v1_cSNPs --bmerge US_control_300_cSNPs.bed US_control_300_cSNPs.bim US_control_300_cSNPs.fam --allow-no-sex --make-bed --out US_control_550v3_550v1_300_cSNPs
      s
 ##
   ./plink --bfile US_control_550v3_550v1_cSNPs --bmerge hapmap_cSNPs.bed hapmap_cSNPs.bim hapmap_cSNPs.fam --allow-no-sex --make-bed --out US_control_550v3_550v1_300_hapmap_cSNPs
     
### RUN eigenstrat with US controls + hapmap   non-flipped
      ./../bin/convertf -p ../CONVERTF/par.BED.EIGENSTRAT2
     ./../bin/smartpca.perl -i US_control_550v3_550v1_300_hapmap_cSNPs.eigenstratgeno -a US_control_550v3_550v1_300_hapmap_cSNPs.snp -b US_control_550v3_550v1_300_hapmap_cSNPs.ind -o US_control_550v3_550v1_300_hapmap_cSNPs.pca -p US_control_550v3_550v1_300_hapmap_cSNPs.plot -e US_control_550v3_550v1_300_hapmap_cSNPs.eval -l US_control_550v3_550v1_300_hapmap_cSNPs.log   &

    ########## Get which strands need to be changed: GOTO topbot.r
   ### test true_flips.txt
./plink --bfile US_control_550v3_550v1_300_hapmap_cSNPs --flip true_flips.txt --make-bed --allow-no-sex --out US_control_550v3_550v1_300_hapmap_cSNPspf
./plink --bfile US_combined --exclude to_exclude_SNPs.txt --make-bed --allow-no-sex --out US_combined_cSNPs
./plink --bfile US_control_550v3_550v1_300_hapmap_cSNPspf --bmerge US_combined_cSNPs.bed US_combined_cSNPs.bim US_combined_cSNPs.fam --allow-no-sex --make-bed --out merge_US_550v3_550v1_300_hapmap_cSNPspf


### RUN eigenstrat with US data + US controls + hapmap   flipped data
./../bin/convertf -p ../CONVERTF/par.BED.EIGENSTRAT3
./../bin/smartpca.perl -i merge_US_550v3_550v1_300_hapmap_cSNPspf.eigenstratgeno -a merge_US_550v3_550v1_300_hapmap_cSNPspf.snp -b merge_US_550v3_550v1_300_hapmap_cSNpspf.ind -o merge_US_550v3_550v1_300_hapmap_cSNPspf.pca -p merge_US_550v3_550v1_300_hapmap_cSNPspf.plot -e merge_US_550v3_550v1_300_hapmap_cSNPspf.eval -l merge_US_550v3_550v1_300_hapmap_cSNPspf.log   &

## moved to SNPMAX03 are finish other merges: meregUk has UK and WTCCC data:
plink --bfile mergeUK --exclude to_exclude_SNPs.txt --make-bed --allow-no-sex --out mergeUK_cSNPs
plink --bfile merge_US_550v3_550v1_300_hapmap_cSNPspf --bmerge mergeUK_cSNPs.bed mergeUK_cSNPs.bim mergeUK_cSNPs.fam --allow-no-sex --make-bed --out merge_ALL_hapmap


### RUN eigenstrat with UK data + US data + WTCCC controls + US controls + hapmap   flipped data
./../bin/convertf -p ../CONVERTF/par.BED.EIGENSTRAT4
./../bin/smartpca.perl -i merge_ALL_hapmap.eigenstratgeno -a merge_ALL_hapmap.snp -b merge_ALL_hapmap.ind -o merge_ALL_hapmap.pca -p merge_ALL_hapmap.plot -e merge_ALL_hapmap.eval -l merge_ALL_hapmap.log   &




 ./plink --bfile merge_ALL_hapmap --remove Hap_6sig_delete_MergeAll.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1


### remove 6p arm      7920 snps                                  ls
./plink --bfile merge_ALL_hapmap --exclude p6_remove.txt --make-bed --allow-no-sex --out merge_ALL_hapmap_no6p

### below remove entries based on plink criteria:
./plink --bfile merge_ALL_hapmap_no6p  --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --make-bed --allow-no-sex --out merge_ALL_hapmap_no6pF
./../bin/convert*f -p ../CONVERTF/par.BED.EIGENSTRAT6
./../bin/./smartpca.perl -i merge_ALL_hapmap_no6pF.eigenstratgeno -a merge_ALL_hapmap_no6pF.snp -b merge_ALL_hapmap_no6pF.ind -o merge_ALL_hapmap_no6pF.pca -p merge_ALL_hapmap_no6pF.plot -e merge_ALL_hapmap_no6pF.eval -l merge_ALL_hapmap_no6pF.log   &


./plink --bfile merge_ALL_hapmap --remove Hap_6sig_delete_MergeAll.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1


############ run assoc  The perm run excludes addiotion points marked as bad in the pca for controls
 ./plink --bfile merge_ALL_hapmap --keep include_samples.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --perm --out merge_ALL_hapmap_pop1_perm
 ./plink --bfile merge_ALL_hapmap --keep include_samples.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --out merge_ALL_hapmap_pop1

./../bin/smartpca.perl -i merge_ALL_hapmap_no6p.eigenstratgeno -a merge_ALL_hapmap_no6p.snp -b merge_ALL_hapmap_no6p.ind -o merge_ALL_hapmap_no6p.pca -p merge_ALL_hapmap_no6p.plot -e merge_ALL_hapmap_no6p.eval -l merge_ALL_hapmap_no6p.log   &


### eigenstrat run
./../CONVERTF/ind2pheno.perl merge_ALL_hapmap_no6pF.ind  merge_ALL_hapmap_no6pF.pheno
./../bin/eigenstrat -i merge_ALL_hapmap_no6pF.eigenstratgeno -j merge_ALL_hapmap_no6pF.pheno -p merge_ALL_hapmap_no6pF.pca -o merge_ALL_hapmap_no6pF.chisq   &


#### Multidimensional scaling with plink
 nohup ./plink --bfile merge_ALL_hapmap_no6p --allow-no-sex --read-genome merge_ALL_hapmap_no6p_G.genome --cluster --mds-plot 8 --out merge_ALL_hapmap_no6p.mds &
 nohup ./plink --bfile merge_ALL_hapmap_no6p --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --read-genome merge_ALL_hapmap_no6p_G.genome --cluster --mds-plot 8 --out merge_ALL_hapmap_no6pF.mds &

############# different SNP lists


## 01/08 file analysis
############## final data merge ##### ## Does not include HAPHAP except for CEU from with icontrolDB data
./plink --bfile US_control_550v3_550v1_300_hapmap_cSNPspf --bmerge all_affected.bed all_affected.bim6pF all_affected.fam --allow-no-sex --make-bed --out merge_all_final_cSNPs_noWTCCC
./plink --bfile WTCCC --exclude to_exclude_SNPs.txt --make-bed --allow-no-sex --out WTCCC_cSNPs
./plink --bfile merge_all_final_cSNPs_noWTCCC --bmerge WTCCC_cSNPs.bed WTCCC_cSNPs.bim WTCCC_cSNPs.fam --allow-no-sex --make-bed --out merge_all_final_cSNPs_aSNPs &
./plink --bfile merge_all_final_cSNPs_aSNPs --exclude p6_remove.txt --make-bed --allow-no-sex --out merge_all_final_cSNP_no6p
                ### filter data for bad samples and snps - 6938 originally #####
./plink --bfile merge_all_final_cSNP_no6p  --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --make-bed --allow-no-sex --out merge_all_final_cSNP_no6pF
                ### 278016 snps left - see log file
./../bin/convertf -p ../CONVERTF/par.BED.EIGENSTRAT7
nohup ./../bin/./smartpca.perl -i merge_all_final_cSNP_no6pF.eigenstratgeno -a merge_all_final_cSNP_no6pF.snp -b merge_all_final_cSNP_no6pF.ind -o merge_all_final_cSNP_no6pF.pca -p merge_all_final_cSNP_no6pF.plot -e merge_all_final_cSNP_no6pF.eval -l merge_all_final_cSNP_no6pF.log   &
nohup ./../bin/./smartpca.perl -i merge_all_final_cSNP_no6pF.eigenstratgeno -a merge_all_final_cSNP_no6pF.snp -b merge_all_final_cSNP_no6pF.ind -m 0 -o merge_all_final_cSNP_no6pF_nR.pca -p merge_all_final_cSNP_no6pF_nR.plot -e merge_all_final_cSNP_no6pF_nR.eval -l merge_all_final_cSNP_no6pF_nR.log   &

nohup plink --bfile merge_all_final_cSNP_no6pF --allow-no-sex --genome --out merge_all_final_cSNP_no6pF_G &

#### Merged control file:   ## for later use
nohup ./plink --bfile US_control_550v3_550v1_300_hapmap_cSNPspf --bmerge WTCCC_cSNPs.bed WTCCC_cSNPs.bim WTCCC_cSNPs.fam --allow-no-sex --make-bed --out merge_all_controls
nohup ./plink --bfile US_control_550v3_550v1_300_hapmap_cSNPspf --bmerge WTCCC_cSNPs.bed WTCCC_cSNPs.bim WTCCC_cSNPs.fam --allow-no-sex --make-bed --out merge_all_controls

### cluster controls
plink --bfile merge_all_controls --exclude p6_remove.txt --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --make-bed --allow-no-sex --out merge_all_controls_cSNP_no6pF
../../eigensoft/bin/convertf -p ../../eigensoft/CONVERTF/par.BED.EIGENSTRAT2
nohup ../../eigensoft/bin/smartpca.perl -i merge_all_controls_cSNP_no6pF.eigenstratgeno -a merge_all_controls_cSNP_no6pF.snp -b merge_all_controls_cSNP_no6pF.ind -o merge_all_controls_cSNP_no6pF.pca -p merge_all_controls_cSNP_no6pF.plot -e merge_all_controls_cSNP_no6pF.eval -l merge_all_controls_cSNP_no6pF.log   &
nohup ../../eigensoft/bin/smartpca.perl -i merge_all_controls_cSNP_no6pF.eigenstratgeno -a merge_all_controls_cSNP_no6pF.snp -b merge_all_controls_cSNP_no6pF.ind -m 0 -o merge_all_controls_cSNP_no6pF_nR.pca -p merge_all_controls_cSNP_no6pF_nR.plot -e merge_all_controls_cSNP_no6pF_nR.eval -l merge_all_controls_cSNP_no6pF_nR.log   &



#######


#############ls Cluster just the affected samples:
plink --bfile all_affected  --exclude p6_remove.txt --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --make-bed --allow-no-sex --out all_affected_F
../../eigensoft/bin/convertf -p ../../eigensoft/CONVERTF/par.BED.EIGENSTRAT
../../eigensoft/bin/smartpca.perl -i all_affected_F.eigenstratgeno -a all_affected_F.snp -b all_affected_F.ind -o all_affected_F.pca -p all_affected_F.plot -e all_affected_F.eval -l all_affected_F.log   &
plink --bfile all_affected --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --out all_affected_F_all
plink --bfile merge_all_controls --keep include_samples.txt --allow-no-sex --make-bed --out merge_all_controls_pop1
plink --bfile merge_all_controls_pop1 --bmerge all_affected_F.bed all_affected_F.bim all_affected_F.fam --allow-no-sex --make-bed --out prelim_run
plink --bfile prelim_run --remove exclude_samples.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --out merge_aff_selcontrols_6sig
plink --bfile prelim_run --remove exclude_samples.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --out merge_aff_selcontrols_6sig_0p021
plink --bfile prelim_run --remove exclude_samples.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --out merge_aff_selcontrols_6sig_0p08

#############ls Cluster just the control samples:
./plink --bfile merge_all_final_cSNP_no6pF --remove HapMap_samples.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --out merge_all_final_cSNP_no6pF_noHapMap

./plink --bfile merge_all_final_cSNP_no6pF --keep include_samples.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --out merge_all_final_cSNP_no6pF_keep2
./plink --bfile merge_all_final_cSNP_no6pF --keep include_samples.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --out merge_all_final_cSNP_no6pF_cut2
./plink --bfile merge_all_final_cSNP_no6pF --keep include_samples.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --out merge_all_final_cSNP_no6pF_keep3
./plink --bfile merge_all_final_cSNP_no6pF --keep include_samples.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --out merge_all_final_cSNP_no6pF_cut3
./plink --bfile merge_all_final_cSNP_no6pF --keep include_samples.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --out merge_all_final_cSNP_no6pF_cut4
 ./plink --bfile merge_all_final_cSNP_no6pF --remove exclude_samples_6sig.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --out merge_all_final_cSNP_no6pF_6sig

 ./plink --bfile merge_all_final_cSNPs_aSNPs --keep include_samples.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --out merge_all_final_cSNPs_aSNPs_cut2

### eigenstrat run
./../CONVERTF/ind2pheno.perl merge_all_final_cSNP_no6pF.ind  merge_all_final_cSNP_no6pF.pheno
./../bin/eigenstrat -i merge_all_final_cSNP_no6pF.eigenstratgeno -j merge_all_final_cSNP_no6pF.pheno -p merge_all_final_cSNP_no6pF.pca -o merge_all_final_cSNP_no6pF.chisq   &

#cluster with plink
nohup plink --bfile merge_all_final_cSNP_no6pF --allow-no-sex --genome --out merge_all_final_cSNP_no6pF_G &

nohup ./plink --bfile merge_all_final_cSNP_no6pF --allow-no-sex --read-genome merge_all_final_cSNP_no6pF_G.genome --cluster --mds-plot 8    &
	
./plink --bfile merge_all_final_cSNP_no6pF --allow-no-sex --exclude strange_snps.txt  --remove related_samples.txt --make-bed --out merge_all_final_cSNP_no6pF_mafF2

nohup ./plink --bfile merge_all_final_cSNP_no6pF_mafF2 --allow-no-sex --read-genome merge_all_final_cSNP_no6pF_G.genome --cluster --mds-plot 8

./plink --bfile merge_all_final_cSNP_no6pF --allow-no-sex --exclude swap_snps.txt  --remove related_samples.txt --make-bed --out merge_all_final_cSNP_no6pF_mafF3

nohup ./plink --bfile merge_all_final_cSNP_no6pF_mafF3 --allow-no-sex --read-genome merge_all_final_cSNP_no6pF_G.genome --cluster --mds-plot 8

############# plink chromosome 2 clustering
nohup ./plink --bfile merge_all_chr2F  --genome --out merge_all_chr2F_G > nohup_chr2.out &
nohup ./plink --bfile merge_all_chr2F_noswaps  --genome --out merge_all_chr2F_noswaps_G > nohup_chr2_noswaps.out &


./plink --bfile merge_all_chr2F --exclude swap_snps_CHR2.txt --make-bed --out merge_all_chr2F_noswaps

 nohup ./plink --bfile merge_all_chr2F  --read-genome merge_all_chr2F_G.genome --cluster --mds-plot 4    &
 nohup ./plink --bfile merge_all_chr2F_noswaps  --read-genome merge_all_chr2F_noswaps_G.genome --cluster --mds-plot 4 >nuhup_noswaps.out   &

################plink MDS check

nohup  .././plink --bfile merge_ALL_HapMap --keep subset_mds_samples.txt --allow-no-sex --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --make-bed --out merge_subset  > nohup.gt0.out     &
.././plink --bfile merge_subset --exclude p6_remove.txt  --make-bed --out merge_subset_no6p
.././plink --bfile merge_subset --exclude p6_remove_all.txt  --make-bed --out merge_subset_no6pALL

###### re-run the big job -7 days run time
.././plink --bfile  merge_ALL_HapMap --exclude p6_remove_all.txt --allow-no-sex --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --make-bed --out  merge_ALL_HapMap_no6pF
nohup .././plink --bfile merge_ALL_HapMap_no6pF  --genome --out merge_ALL_HapMap_no6pF > nohup_big.out &
 nohup .././plink --bfile merge_ALL_HapMap_no6pF  --read-genome merge_ALL_HapMap_no6pF.genome --cluster --mds-plot 4  > big_cluster.out  &
 
 #### frequencies for all and population modling runs
.././plink --bfile  merge_ALL_HapMap  --allow-no-sex  --freq --out merge_ALL_HapMap_FREQS > junk.out &
.././plink --bfile  merge_ALL_HapMap  --allow-no-sex --assoc --maf 0.0 --geno 1.0 --out merge_ALL_HapMap_ASSOC 

.././plink --bfile merge_subset_NF --exclude p6_remove.txt  --make-bed --out merge_subset_no6pNF
 nohup .././plink --bfile merge_subset_no6pNF  --genome --out merge_subset_no6pNF_G > nohup_6p.out &


.././plink --bfile  merge_ALL_HapMap  --keep hapmap_set.txt --recodeAD --out hapmap_pop

http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#recode











 ## run everything on SNPMAX03 
 ## small job on SNPMAX01 to test version problem.
  ./plink --bfile merge_ALL_hapmap_no6p --allow-no-sex --genome --out merge_ALL_hapmap_no6p_G &
 nohup ./plink --bfile merge_ALL_hapmap_no6p --allow-no-sex --read-genome merge_ALL_hapmap_no6p_G.genome --cluster --mds-plot 8 --out merge_ALL_hapmap_no6p.mds &



#### try no filtering
 nohup  .././plink --bfile merge_ALL_HapMap --keep subset_mds_samples.txt --allow-no-sex --make-bed --out merge_subset_NF  > nohup.NF.out     &
 .././plink --bfile merge_subset_NF --exclude p6_remove.txt  --make-bed --out merge_subset_no6pNF
 nohup .././plink --bfile merge_subset_no6pNF --noweb  --genome --out merge_subset_no6pNF_G > nohup_6p.out &
nohup  .././plink --bfile merge_subset_no6pNF   --read-genome merge_subset_no6pNF_G.genome    --cluster --mds-plot 4  --out merge_subset_no6pNF_C > nohupNF.out  &

 ############## ON SNPMAX01   plink/plink-0.99.../MDS_check
nohup plink --bfile merge_subset_no6pNF --noweb  --genome --out merge_subset_no6pNFold_G > nohup_60.out &
nohup plink --bfile merge_subset_no6pNF --noweb  --read-genome merge_subset_no6pNFold_G.genome    --cluster --mds-plot 4  --out merge_subset_no6pNF_C > nohupNF.out  &

#### trp with 6p included and no filtering
nohup plink --bfile merge_subset_NF --noweb  --genome --out merge_subset_NF_G > nohup_8.out &
nohup plink --bfile merge_subset_NF --noweb  --read-genome  merge_subset_NF_G.genome    --cluster --mds-plot 4  --out merge_subset_NF_C > nohup9.out  &


##!!!!!!!!!!   clustering not reproducred with filter/nofilter of MAF etc; 6p/no6p/incomplete 6p removal; plink version; combinations of 6p removal and filters
####

nohup .././plink --bfile merge_subset_no6p  --genome --out merge_subset_no6pNF_G > nohup_6p.out &
nohup .././plink --bfile merge_subset_no6pALL  --genome --out merge_subset_no6pALL_G > nohup_6pALL.out &

nohup ../.././plink --bfile merge_subset_no6pALL  --read-genome merge_subset_no6pALL_G.genome --cluster --mds-plot 4  --out merge_subset_no6pALL_C > nohup2.out  &
nohup .././plink --bfile merge_subset_no6p     --read-genome merge_subset_no6p_G.genome    --cluster --mds-plot 4  --out merge_subset_no6p_C > nohup1.out  &





./plink --bfile merge_all_final_cSNP_no6pF --keep include_top_cluster_samples.txt --allow-no-sex --freq --out top_cluster


 perl -ne 'BEGIN{$i=0}; $i=$i+1 ;if($i % 10 == 0){print " $_"} ' test.txt

perl -ne 'BEGIN{$i=0}; $i=$i+1 ;if($i % 10 == 0){print $i "\n"} ' test.txt

perl -ne 'BEGIN{$i=0}; $i=$i+1 ;print $i "\n" ' test.txt


 ######### split the raw file wc -l
 perl -pe 'BEGIN { @index=(1); for($i=6;$i<=20;$i=$i+2){@index=(@index,$i)} }   @w=split; $_=join("\t",@w[@index]); print "\n" ' test.txt
 perl -pe 'BEGIN { @index=(1); for($i=6;$i<=642404;$i=$i+2){@index=(@index,$i)} }   @w=split; $_=join("\t",@w[@index]); print "\n" ' test.txt

 perl -pe 'BEGIN {@index<-for ($i=6;i<=642404;i=i+2}   @w=split; $_=join("\t",@w[1,2,3,4,5,6,7]); print "\n" ' test.txt

######### get info from genome file...
 cut -d' ' -f1,3,8 test.txt >> result.txt");
perl -pe '@w=split; $_=$w[1]; print "\n" ' test.txt   # extract column from a file
perl -pe '@w=split; $_=$w[1]."\t".$w[3]."\t".$w[7]; print "\n" ' test.txt # extract columns from a file

perl -pe '@w=split; $_=$w[1]."\t".$w[3]."\t".$w[7]; print "\n" ' merge_all_final_cSNP_no6pF_G.genome  > merge_all_final_cSNP_no6pF_piHAT.txt
perl -pe '@w=split; $_=$w[1]."\t".$w[3]."\t".$w[7]; print "\n" ' test.txt  > merge_all_final_cSNP_no6pF_piHAT.txt

 head merge_all_final_cSNP_no6pF_G.genome  >test.txt

perl -ne '@w=split; if ($w[2] > 0.05) {print $w[0]."\t".$w[1]."\t".$w[2], "\n"} ' merge_all_final_cSNP_no6pF_piHAT.txt > merge_all_final_related.txt

perl -ne '@w=split; {print $w[2],"\n"} ' test.txt
perl -ne '@w=split; if($w[2]> 0){print $w[1],"\n"} ' test.txt    # must be 'ne' NOT en  'pe' prints $_  'ane' splits automatically into @F

perl -ne '@w=split; if($w[2]> 0){print $w[1],"\n"} ' test.txt    # must be 'ne' NOT en  'pe' prints $_  'ane' splits automatically into @F

perl -ne '@w=split; if ($w[7] > 0.05) {print $w[0]."\t".$w[1]."\t".$w[3]."\t".$w, "\n"} ' merge_all_final_cSNP_no6pF_G.genome > related.txt

perl -ne '@w=split; if ($w[7] > 0.05) {print $_} ' merge_all_final_cSNP_no6pF_G.genome > related.txt

perl -ne '@w=split; if ($w[7] > 0.005) {print $_} ' test.txt  >out.txt
perl -ne '@w=split; if ($w[7] > 0.005) {print $w[0]."\t".$w[1]."\t".join("\t",$w[5 .. 17]),"\t"} ' test.txt  >out.txt
 ######## REMEMBER PERL STARTS FROM INDEX 0
perl -ne 'BEGIN { $i=0} if ($i<9) then {print $_} else { @w=split; if ($w[7] > 0.005) {print $w[0]."\t".$w[1]."\t".join("\t",$w[5 .. 17]),"\t"} ' test.txt  >out.txt

perl -pe 'BEGIN {$i=0} if ($i <  10){$i++; chomp($_); print "\n" }else{  @w=split; $_=split(//,$w[0]) ; print "\n"}  ' test.txt
perl -ne 'BEGIN {$i=0} if ($i <  10){$i++; chomp($_); print "\n" }else{  @w=split; $count=split(//,$w[0]) ; if($count>8 & $w[1].eq."PEX0020"){ print $w[0]."\t".$w[1]."\t".$w[2]."\t".$w[3], "\n"}}  ' test.txt

perl -ne 'BEGIN {$i=0} if ($i <  10){$i++; chomp($_); print $_,"\n" }else{  @w=split; $count=split(//,$w[0]) ; if(($count>9) && ($w[1] eq "PEX_0020")){ print $w[0]."\t".$w[1]."\t".$w[2]."\t".$w[3], "\n" }}  ' test.txt

perl -ne 'BEGIN {$i=0} if ($i <  10){$i++; chomp($_); print $_,"\n" }else{  @w=split; $count=split(//,$w[0]) ; if(($count>16) && ($w[1] eq "PEX_0020")){ print $w[0]."\t".$w[1]."\t".$w[2]."\t".$w[3], "\n" }}  ' final_report.txt  > counts.txt

perl -ne '@w=split; if ($w[7] > 0.05) {print $_} ' PAX_G.genome > PAX_related.txt


./plink --bfile merge_all_final_cSNP_no6pF --keep include_samples.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --out merge_all_final_cSNP_no6pF_cut2_no

#permutation in plik                                       hed
nohup  ./plink --bfile merge_all_final_cSNP_no6pF --keep include_samples.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --perm --out merge_all_final_cSNP_no6pF_cut2_perm > nohup_perm.log &


# eigenstrat with 6p arm in place,
./plink --bfile merge_all_final_cSNPs_aSNPs --remove exclude_samples_6sig.txt --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --make-bed --allow-no-sex --out merge_all_final_cSNPs_F
./plink --bfile merge_all_final_cSNPs_aSNPs --remove exclude_samples_6sig.txt --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --make-bed --allow-no-sex --out merge_all_final_cSNPs_F2

                ### 278016 snps left - see log file
./../bin/convertf -p ../CONVERTF/par.BED.EIGENSTRAT8
nohup ./../bin/./smartpca.perl -i merge_all_final_cSNPs_F.eigenstratgeno -a merge_all_final_cSNPs_F.snp -b merge_all_final_cSNPs_F.ind -o merge_all_final_cSNPs_F.pca -p merge_all_final_cSNPs_F.plot -e merge_all_final_cSNPs_F.eval -l merge_all_final_cSNPs_F.log   &
nohup ./../bin/./smartpca.perl -i merge_all_final_cSNPs_F.eigenstratgeno -a merge_all_final_cSNPs_F.snp -b merge_all_final_cSNPs_F.ind -m 0 -o merge_all_final_cSNPs_F_nR.pca -p merge_all_final_cSNPs_F_nR.plot -e merge_all_final_cSNPs_F_nR.eval -l merge_all_final_cSNPs_F_nR.log   &
### eigenstrat run
./../bin/convertf -p ../CONVERTF/par.BED.EIGENSTRAT9    # remove bad snps  6p and 8 invert
./../CONVERTF/ind2pheno.perl merge_all_final_cSNPs_F2.ind  merge_all_final_cSNPs_F2.pheno  #new pheno
./../bin/eigenstrat -i merge_all_final_cSNPs_F2.eigenstratgeno -j merge_all_final_cSNPs_F2.pheno -p merge_all_final_cSNPs_F.pca -o merge_all_final_cSNPs_F2.chisq   &

#### plink run with MHC:
./plink --bfile merge_all_final_cSNPs_F --keep include_samples.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --out merge_all_final_cSNPs_F_cut2


########################REDO WITH REMOVE SAMPLES
#### remove related samples # has 6p arm #  6sig removed # related removed
#./plink --bfile merge_all_final_cSNPs_F --remove related_samples.txt --allow-no-sex --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --out shit
./plink --bfile merge_all_final_cSNPs_F --remove exclude_samples_6sig_related.txt --allow-no-sex --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --make-bed --out merge_all_final_cSNPs_F_NR


### First WAY
./../bin/convertf -p ../CONVERTF/par.BED.EIGENSTRAT10
nohup ./../bin/./smartpca.perl -i merge_all_final_cSNPs_F_NR.eigenstratgeno -a merge_all_final_cSNPs_F_NR.snp -b merge_all_final_cSNPs_F_NR.ind -o merge_all_final_cSNPs_F_NR2.pca -p merge_all_final_cSNPs_F_NR2.plot -e merge_all_final_cSNPs_F_NR2.eval -l merge_all_final_cSNPs_F_NR2.log > nohup_NR2.out  &
#nohup ./../bin/./smartpca.perl -i merge_all_final_cSNPs_F_NR.eigenstratgeno -a merge_all_final_cSNPs_F_NR.snp -b merge_all_final_cSNPs_F_NR.ind -m 0 -o merge_all_final_cSNPs_F_nR.pca -p merge_all_final_cSNPs_F_nR.plot -e merge_all_final_cSNPs_F_nR.eval -l merge_all_final_cSNPs_F_nR.log   &
 ### eigenstrat run
./../bin/convertf -p ../CONVERTF/par.BED.EIGENSTRAT10    # remove bad snps  6p and 8 invert
./../CONVERTF/ind2pheno.perl merge_all_final_cSNPs_F_NR.ind  merge_all_final_cSNPs_F_NR.pheno  #new pheno
./../bin/eigenstrat -i merge_all_final_cSNPs_F_NR.eigenstratgeno -j merge_all_final_cSNPs_F_NR.pheno -p merge_all_final_cSNPs_F_NR2.pca -o merge_all_final_cSNPs_F_NR2.chisq   &
########################REDO WITH REMOVE SAMPLES



### Second WAY
### this combination failed ... thee is no reson why it shoudl have

##YES what happens is the the removal of SNPs has fucked it up
### smart pca with packed data no removeal of 6sig data as already removed
nohup ./../bin/./smartpca.perl -i merge_all_final_cSNPs_F_NR.bed -a merge_all_final_cSNPs_F_NR.bim -b merge_all_final_cSNPs_F_NR.fam -m 0 -o merge_all_final_cSNPs_F_NR.pca -p merge_all_final_cSNPs_F_NR.plot -e merge_all_final_cSNPs_F_NR.eval -l merge_all_final_cSNPs_F_NR.log   & cl_related.out
 ### eigenstrat run
./../bin/convertf -p ../CONVERTF/par.BED.EIGENSTRAT10    # remove bad snps  6p and 8 invert
./../CONVERTF/ind2pheno.perl merge_all_final_cSNPs_F_NR.ind  merge_all_final_cSNPs_F_NR.pheno  #new pheno
./../bin/eigenstrat -i merge_all_final_cSNPs_F_NR.eigenstratgeno -j merge_all_final_cSNPs_F_NR.pheno -p merge_all_final_cSNPs_F_NR.pca -o merge_all_final_cSNPs_F_NR.chisq   &
########################REDO WITH REMOVE SAMPLES








nohup  ./plink --bfile merge_all_final_cSNPs_F --keep include_samples.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --perm --out merge_all_final_cSNPs_F_cut2_perm > nohup_perm.log &
nohup  ./plink --bfile merge_all_final_cSNPs_F --keep include_samples.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --perm --out merge_all_final_cSNPs_F_cut2_noRelated_perm > nohup_perm.log &


#nohup ./../bin/./smartpca.perl -i merge_all_final_cSNP_F.eigenstratgeno -a merge_all_final_cSNP_F.snp -b merge_all_final_cSNP_F.ind -m 0 -o merge_all_final_cSNP_F_nR.pca -p merge_all_final_cSNP_F_nR.plot -e merge_all_final_cSNP_F_nR.eval -l merge_all_final_cSNP_F_nR.log   &
nohup plink --bfile merge_all_final_cSNP_F --allow-no-sex --genome --out merge_all_final_cSNP_no6pF_G &


#### chr 8 invesion clusters
./plink --bfile merge_all_final_cSNP_no6pF --keep include_top_cluster_samples.txt --allow-no-sex --freq --out top_cluster
./plink --bfile merge_all_final_cSNP_no6pF --keep include_mid_cluster_samples.txt --allow-no-sex --freq --out mid_cluster
./plink --bfile merge_all_final_cSNP_no6pF --keep include_bot_cluster_samples.txt --allow-no-sex --freq --out _cluster


#### mds cluster problem  ### first remove known pop stratification and related samples
## this will hoepfully ensure pop stratification is not encountered.
 ./plink --bfile merge_all_final_cSNP_no6pF --keep include_samples.txt --allow-no-sex --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --make-bed --out merge_all_final_cSNP_no6pF_cut2

nohup ./plink --bfile merge_all_final_cSNP_no6pF_cut2 --keep get_mds_set_lt0.txt --allow-no-sex --freq --out lt0_cluster_cut2  > nohup.lt0.out  &
nohup  ./plink --bfile merge_all_final_cSNP_no6pF_cut2 --keep get_mds_set_gt0.txt --allow-no-sex --freq --out gt0_cluster_cut2  > nohup.gt0.out     &

nohup ./plink --bfile merge_all_final_cSNP_no6pF_cut2 --keep get_mds_set_lt0.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --out lt0_cluster_cut2_assoc  > nohup.lt0a.out  &
nohup  ./plink --bfile merge_all_final_cSNP_no6pF_cut2 --keep get_mds_set_gt0.txt --allow-no-sex --assoc --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --out gt0_cluster_cut2_assoc  > nohup.gt0a.out     &

nohup ./plink --bfile merge_all_final_cSNP_no6pF --keep get_mds_set_lt0.txt --allow-no-sex --freq --out lt0_cluster  > nohup.lt0.out  &
nohup  ./plink --bfile merge_all_final_cSNP_no6pF --keep get_mds_set_gt0.txt --allow-no-sex --freq --out gt0_cluster  > nohup.gt0.out     &

nohup ./plink --bfile merge_all_final_cSNP_no6pF --keep gt0_split_1.txt --allow-no-sex --freq --out gt0_split_1  > nohup.gt0.1.out  &
nohup  ./plink --bfile merge_all_final_cSNP_no6pF --keep gt0_split_2.txt --allow-no-sex --freq --out gt0_split_2  > nohup.gt0.2.out     &

















##### remove pop-strat that might cause MAF to change between groups (as one group has WTCCC which has no pop-strat
 ./plink --bfile merge_all_final_cSNP_no6pF --keep include_samples_003.txt --allow-no-sex --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --make-bed --out merge_all_final_cSNP_no6pF_003

nohup ./plink --bfile merge_all_final_cSNP_no6pF_003 --keep get_mds_set_lt0.txt --allow-no-sex --freq --out lt0_cluster_003  > nohup.lt0_003.out  &
nohup  ./plink --bfile merge_all_final_cSNP_no6pF_003 --keep get_mds_set_gt0.txt --allow-no-sex --freq --out gt0_cluster_003  > nohup.gt0_003.out     &

nohup ./plink --bfile merge_all_final_cSNP_no6pF_003 --keep gt0_split_1.txt --allow-no-sex --freq --out gt0_split_1_003  > nohup.gt0.1_003.out  &
nohup  ./plink --bfile merge_all_final_cSNP_no6pF_003 --keep gt0_split_2.txt --allow-no-sex --freq --out gt0_split_2_003  > nohup.gt0.2_003.out     &

nohup ./plink --bfile merge_all_final_cSNP_no6pF  --allow-no-sex --freq --out  merge_all_final_cSNP_no6pF_nR


################## PEX RUNS
 ./plink --bfile WTCCC --bmerge PAX.bed PAX.bim PAX.fam --allow-no-sex --make-bed --out PEX_WTCCC
 ./plink --bfile PEX_WTCCC --bmerge hapmap.bed hapmap.bim hapmap.fam --allow-no-sex --make-bed --out PEX_WTCCC_hapmap
 ./plink --bfile hapmap --flip PEX_WTCCC_hapmap.missnp --make-bed --allow-no-sex --out hapmap_f
./plink --bfile PEX_WTCCC --bmerge hapmap_f.bed hapmap_f.bim hapmap_f.fam --allow-no-sex --make-bed --out PEX_WTCCC_hapmapf
./plink --bfile PEX_WTCCC_hapmapf  --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --make-bed --allow-no-sex --out PEX_WTCCC_hapmapf_F

  ./plink --bfile WTCCC --exclude to_exclude_SNPs.txt --make-bed --allow-no-sex --out mergeUS_control_cSNPs

 ./plink --bfile PAX --allow-no-sex --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --assoc --out PAX

 nohup ./../bin/./smartpca.perl -i PAX.bed -a PAX.bim -b PAX.fam -o PAX.pca -p PAX.plot -e PAX.eval -l PAX.log > nohip_PAX.out  &
 nohup ./plink --bfile PAX --allow-no-sex --genomepwd --out PAX_G > nohup_genome.out &

 nohup ./../bin/./smartpca.perl -i PEX_WTCCC_hapmapf_F.bed -a PEX_WTCCC_hapmapf_F.bim -b PEX_WTCCC_hapmapf_F.fam -o PEX_WTCCC_hapmapf_F.pca -p PEX_WTCCC_hapmapf_F.plot -e PEX_WTCCC_hapmapf_F.eval -l PEX_WTCCC_hapmapf_F.log > nohup_PEX_WTCCC_hapmapf_F.out  &

./plink --bfile merge_all_final_cSNPs_aSNPs --remove exclude_samples_6sig.txt --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00005 --mind 0.1 --make-bed --allow-no-sex --out merge_all_final_cSNPs_F


  .././plink --bfile PEX_WTCCC_cSNP --allow-no-sex --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00001 --mind 0.1 --assoc --out PEX_WTCCC_cSNPf
    .././plink --bfile PEX_WTCCC_cSNP --allow-no-sex --maf 0.01 --max-maf 1 --geno 0.10 --hwe 0.00001 --mind 0.1 --make-bed --out PEX_WTCCC_cSNPf
 nohup ./../../bin/./smartpca.perl -i PEX_WTCCC_cSNPf.bed -a PEX_WTCCC_cSNPf.bim -b PEX_WTCCC_cSNPf.fam -o PEX_WTCCC_cSNPf.pca -p PEX_WTCCC_cSNPf.plot -e PEX_WTCCC_cSNPf.eval -l PEX_WTCCC_cSNPf.log > nohip_PEX.out  &
 
./../../bin/convertf -p ../../CONVERTF/par.BED.EIGENSTRAT11    # remove bad snps  6p and 8 invert
./../../CONVERTF/ind2pheno.perl PEX_WTCCC_cSNPf.ind  PEX_WTCCC_cSNPf.pheno  #new pheno
./../../bin/eigenstrat -i PEX_WTCCC_cSNPf.eigenstratgeno -j PEX_WTCCC_cSNPf.pheno -p PEX_WTCCC_cSNPf.pca -o PEX_WTCCC_cSNPf.chisq   &

############ convert plink binary to bcos format
./plink_converter.sh PAX.bed PAX.fam PAX.bim 999 PAX_bcos.txt



c300v1<-read.delim("HumanHap300_v1_FinalReport.txt",header=T,skip=10,sep="",fill=TRUE)
c550v1<-read.delim("HumanHap550_v1_FinalReport.txt",header=T,skip=10,sep="",fill=TRUE)
c550v3<-read.delim("HumanHap550_v3_FinalReport.txt",header=T,skip=10,sep="",fill=TRUE)

c300v1<-c300v1[,-2]
c550v1<-c550v1[,-2]
c550v3<-c550v3[,-2]

in300.550v3<-intersect(c300v1,c550v3)
length(c300v1)
[1] 317502
length(test)
[1] 311397

in300.550v3.550v1<-intersect(in300.550v3,c550v1)
> length(in300.550v3.550v1)
[1] 307797


> length(c300v1)
[1] 317502
> length(c550v1)
[1] 555351
> length(c550v3)
[1] 561465


> all.snps<-union(c300v1,c550v1)
> length(all.snps)
[1] 559349
> length(c550v1)
[1] 555351
> all.snps<-union(c550v1,c300v1)
> length(all.snps)
[1] 559349
> all.snps<-union(c550v3,all.snps)
> length(all.snps)
[1] 572135
all.snps<-gsub("rs","RS",all.snps)

to_exclude_snps<-setdiff(all.snps,in300.550v3.550v1)
> length(to_exclude_snps)
[1] 264338
ls

  new<-data.frame( MARKER=in300.550v3.550v1 ,
                       MAF="",
                       CALLRATE="",
                       HWE=""
                        )
                        
  new2<-data.frame( MARKER=to_exclude_snps ,
                       MAF="",
                       CALLRATE="",
                       HWE=""
                        )



 in300.550v3.550v1<-gsub("rs","RS",in300.550v3.550v1)
write.table(in300.550v3.550v1,"common_SNPs.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
> lost_SNPs.txt<-setdiff(c300v1,in300.550v3.550v1)
> length(lost_SNPs.txt)
[1] 9705

 lost_SNPs.txt<-gsub("rs","RS",lost_SNPs.txt)
write.table(lost_SNPs.txt,"lost_SNPs.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

write.table(new2,"excluded_SNPs_SAMPLES.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

we loose about 10000 Snps this way


> mean(UK)
         CHR        CHISQ            P           OR         posn 
8.959151e+00 1.050088e+00 4.954470e-01 1.006591e+00 7.917611e+07 
> mean(US)
         CHR        CHISQ            P           OR         posn 
8.964583e+00 1.568465e+00 4.322679e-01 1.018406e+00 7.914919e+07 
> mean(UK_US)
         CHR        CHISQ            P           OR         posn
8.964548e+00 1.272669e+00 4.645669e-01 1.011915e+00 7.914953e+07



##########################  Q-Q plots
#### the null hypothysis for p-vales follows a uniform distribution
   runif(n, min=0, max=1)  # n is the number of observations

 rchisq(n, df, ncp=0)  # - chisg df=1 from plink
   pchisq(0.2471,1,lower.tail=FALSE) # convert to p-values
 qchisq(0.2106,df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) # convert p values to chisq



indata<-UK
indata.all<-UK.all

indata<-US
indata.all<-US.all

indata<-UK_US
indata.all<-UK_US.all

null<-runif(length(indata[,"P"]))
UKqq<-  qqplot(-log10(null),-log10(indata[,"P"]),plot.it=FALSE)

null<-runif(length(indata.all[,"P"]))
UKqq.all<-  qqplot(-log10(null),-log10(indata.all[,"P"]),plot.it=FALSE)


p1<-plot(c(-1,18),c(-1,119),ylim=c(0,20),xlim=c(0,7),ylab=expression(paste("Observed (",-log[10]," of P value)" ) ),
xlab=expression(paste("Expected (",-log[10]," of P value)" ) ),main="Q-Q plot: UK affected vs. WTCCC controls")
points(UKqq$x,UKqq$y,col="blue",pch=21)
points(UKqq.all$x,UKqq.all$y,col="red",pch=21)
#lines(UKqq.all$x,UKqq.all$x,col="green")
abline(0,1)
legend(4,20,c("6p-arm included","6p-arm excluded"),col=c("red","blue"),
            pch = c(21, 21))

##########



p1<-plot(c(-1,18),c(-1,119),ylim=c(0,20),xlim=c(0,7),ylab=expression(paste("Observed (",-log[10]," of P value)" ) ),
xlab=expression(paste("Expected  (",-log[10]," of P value)" ) ),main="Q-Q plot: US affected vs. US controls")
points(UKqq$x,UKqq$y,col="blue",pch=21)
points(UKqq.all$x,UKqq.all$y,col="red",pch=21)
abline(0,1)
legend(0,20,c("6p-arm included","6p-arm excluded"),col=c("red","blue"),
            pch = c(21, 21))

# qqplot(-log10(null),-log10(indata[,"P"]),plot.it=TRUE)
# qqline(-log10(indata[,"P"]))
# qqline(-log10(null))
 USp<- UKqq
  USp.all<- UKqq.all


#############
p1<-plot(c(-1,18),c(-1,119),ylim=c(0,20),xlim=c(0,7),ylab=expression(paste("Observed (",-log[10]," of P value)" ) ),
xlab=expression(paste("Expected  (",-log[10]," of P value)" ) ),main="Q-Q plot: UK+US affected vs. WTCCC controls")
points(UKqq$x,UKqq$y,col="blue",pch=21)
points(UKqq.all$x,UKqq.all$y,col="red",pch=21)
abline(0,1)
legend(4,20,c("6p-arm included","6p-arm excluded"),col=c("red","blue"),
            pch = c(21, 21))
            
            
 UK_USp<- UKqq
  UK_USp.all<- UKqq.all

#############   ALL


p1<-plot(c(-1,18),c(-1,119),ylim=c(0,20),xlim=c(0,7),ylab=expression(paste("Observed (",-log[10]," of P value)" ) ),
xlab=expression(paste("Expected (",-log[10]," of P value)" ) ),main="Q-Q plot: UK and US preliminary analysis")
points(UKqq$x,UKqq$y,col="blue",pch=20)
points(UKqq.all$x,UKqq.all$y,col="light blue",pch=20)
abline(0,1,col="black")
legend(0,20,c("UK+US 6p-arm excluded","UK+US 6p-arm included","US 6p-arm excluded","US 6p-arm included","UK 6p-arm excluded","UK 6p-arm included"),
col=c("blue","light blue","dark red","red","green","light green"),pch = c(20, 20,20,20,20,20))
points(USp$x,USp$y,col="dark red",pch=20)
points(USp.all$x,USp.all$y,col="red",pch=20)
points(UKp$x,UKp$y,col="green",pch=20)
points(UKp.all$x,UKp.all$y,col="light green",pch=20)


            
legend(6,20,c("6p-arm included","6p-arm excluded"),col=c("red","blue"),
            lty = c(1, 1, 1), merge = TRUE, bg='gray90')

qqplot(-log10(null)[1:100],-log10(UK[,"P"])[1:100],ylab=expression(Observed (-log^{10} of P value)),xlab="Expected (-log_{10} of P value)"  )
qqline(null)


p1<-qqplot(-log10(null)[1:200],-log10(UK[,"P"])[1:200],ylab=expression(paste("Observed (",-log[10]," of P value)" ) ),
xlab=expression(paste("Observed (",-log[10]," of P value)" ) ),main="UK sample vs. WTCCC controls" ,col="green" )


p1<-qqplot(-log10(null)[1:100],-log10(UK[,"P"])[1:100],ylab=expression(paste("Observed (",-log[10]," of P value)" ) ),
xlab=expression(paste("Observed (",-log[10]," of P value)" ) ))

################ Heat plot used :
plot(range(korn_genes.pca$x[,1]),range(korn_genes.pca$x[,2]),xlab="PCA1",ylab="PCA2",main="Spectral Clustering Genes (KORN)")
colours <-  rainbow(7)
#points(sam_genes.pca$x[,1],sam_genes.pca$x[,2],col=colours[sam_genes.cl$cluster],pch=sam_genes.cl$cluster,cex=1.0,bg=colours[sam_genes.cl$cluster]) #colours and symbols
text(korn_genes.pca$x[,1],korn_genes.pca$x[,2],label=rownames(korn_genes.pca$x),col=colours[korn_genes.cl$cluster],cex=0.75) #colours and symbols

 plot(p1$x,p1$y)



 #################################### ASSOCIATION ANALYSIS ######################################

# NCBI build 35=UCSC hg17
# dbSNP 125 is (hg17)
# dbSNP 126 is hg18 
#
# covert locatiosn use UCSC liftOver utility    to convert locatiosn to the previous build
#
snps<-read.delim("ILMN_HumanHap300_SNPlist.txt",header=T,sep="",nrows=5,fill=TRUE)
rownames(snps)<-snps[,1]



file<-"UK.assoc"
file<-"UK_US_aff.assoc"
file<-"US_WTCCC.assoc"
file<-"US_WTCCC_only.assoc"
file<-"mergeUS_550a300.assoc"
file<-"Merge_All_HapMap.assoc"
file<-"Merge_All_HapMap_pop1.assoc"

basfi<-read.delim(file,header=T,sep="",fill=TRUE)
#basfi[,1]<-gsub("RS","rs",basfi[,2])
rownames(basfi)<-basfi[,2]
basfi<-basfi[,-c(2:6)]

merge_all<-basfi
UK.all<-basfi
UK_US.all<-basfi
US.all<-basfi       # Us.all and US.only.all are the same obtained in separate runs
US.only.all<-basfi

file<-"exclude_UK.txt"
file<-"exclude_US.txt"
file<-"exclude_UK_US.txt"

temp<-read.delim(file,header=T,sep="",fill=TRUE)
temp<-temp[temp[,"Check"]==0,]


UK.bad<-as.character(temp[,"SNP"])
US.bad<-as.character(temp[,"SNP"])
UK_US.bad<-as.character(temp[,"SNP"])

bad<-UK.bad
tab<-UK.all
      alist<- apply(as.matrix(bad),1,function(x) grep(x,rownames(tab),,fixed=TRUE)
      locations<-unlist(lapply(alist,max))
      swiss_prot<-RG$genes[,"swiss_prot"]
      geneIDs<-david[locations,"ENTREZ_GENE_ID"]

 #assoc w alls saved  here


UK.all<-UK.all[-c(match(UK.bad,rownames(UK.all))),]
US.all<-US.all[-match(US.bad,rownames(US.all)),]
UK_US.all<-UK_US.all[-match(UK_US.bad,rownames(UK_US.all)),]

bads<-match(US.bad,rownames(US.all))
bads<-bads[!is.na(bads)]
US.all<-US.all[-bads,]

########### Use BioMart ###########
# got.snps<-rownames(UK)
# got.snps<-gsub("RS","rs",got.snps)
# ### Below do to with R
# library(biomaRt)
# listMarts()           # get the data set info
# test=useMart("snp")
# listDatasets(test)
# ######################
# human_snp_126 = useMart("snp",dataset="hsapiens_snp")
# 
# snp_ann<-getBM(attributes=c("refsnp_id","chr_name","chrom_start","chrom_strand"),
# filters=c("refsnp"),
# values=c(got.snps), mart=human_snp_126)



# grap from illumina ############# http://www.illumina.com/pages.ilmn?ID=153
snp.ann<-read.delim("ill_300v2_annotations.txt",header=T,sep="\t",fill=TRUE)
snp.ann[,"Name"]<-gsub("rs","RS",snp.ann[,"Name"])
rownames(snp.ann)<-snp.ann[,"Name"]
snp.ann<-snp.ann[,-c(1,2)]


#######
UK.all[,5]<-snp.ann[rownames(UK.all),"MapInfo"]
colnames(UK.all)[5]<-"posn"

UK_US.all[,5]<-snp.ann[rownames(UK_US.all),"MapInfo"]
colnames(UK_US.all)[5]<-"posn"

US.all[,5]<-snp.ann[rownames(US.all),"MapInfo"]
colnames(US.all)[5]<-"posn"
##########

lenth.names<-nchar(rownames(flip.strands))
alist<- apply(as.matrix(add.flips),1,function(x) flip.strands[(rownames(flip.strands)==x) &  (nchar(x)==lenth.names),] )

merge_all[,5]<-snp.ann[match(rownames(merge_all),rownames(snp.ann)),"MapInfo"]    #merge_all posns in snp.ann
merge_all[,6]<-snp.ann[match(rownames(merge_all),rownames(snp.ann)),"Chr"] 
colnames(merge_all)[5]<-"posn"
colnames(merge_all)[6]<-"CHR.test"


p6.arm<-UK.all[,"CHR"]==6 & UK.all[,"posn"]<60000000
UK<-UK.all[!p6.arm,]

p6.arm<-UK_US.all[,"CHR"]==6 & UK_US.all[,"posn"]<60000000
UK_US<-UK_US.all[!p6.arm,]

p6.arm<-US.all[,"CHR"]==6 & US.all[,"posn"]<60000000
US<-US.all[!p6.arm,]

p6.arm<-merge_all[,"CHR"]==6 & merge_all[,"posn"]<60000000
merge<-merge_all[!p6.arm,]

6p.arm.snps<-merge_all[p6.arm,]
p6.remove<-rownames(p6.arm.snps)
write.table(p6.remove,"p6_remove.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)









                       m



### order
orders<-sort(UK[,"P"],index.return=TRUE)$ix
UK<-UK[orders,]

orders<-sort(UK_US[,"P"],index.return=TRUE)$ix
UK_US<-UK_US[orders,]

orders<-sort(US[,"P"],index.return=TRUE)$ix
US<-US[orders,]


##### remove 6p arm

top<-UK[1:50,]
top[,2]<-US[rownames(top),"P"]
top[,4]<-UK_US[rownames(top),"P"]
colnames(top)<-c("CHR","US","UK","UK_US","posn")


top<-US[1:100,]
top[,2]<-UK[rownames(top),"P"]
top[,4]<-UK_US[rownames(top),"P"]
colnames(top)<-c("CHR","UK","US","UK_US","posn")



write.table(top[,c(1,5,3,2,4)],"top50UK.txt",sep="\t",col.names=TRUE,row.names=TRUE,quote=FALSE)
write.table(UK[1:100,],"UK.txt",sep="\t",col.names=TRUE,row.names=TRUE,quote=FALSE)
write.table(US[1:100,],"US.txt",sep="\t",col.names=TRUE,row.names=TRUE,quote=FALSE)
write.table(UK_US[1:100,],"UK_US.txt",sep="\t",col.names=TRUE,row.names=TRUE,quote=FALSE)

library(annaffy)
############ a table
filename<- "topUS50.html"
title <- "Top 100 in US-UScontrols : with corresponding UK+US-WTCCC and UK-WTCCC for comparison "

newURLs<-lapply(sub("RS","",rownames(top)),function(x) paste('<a href="http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=',x,'">', sep = "", collapse = NULL) )
newURLs<-paste(newURLs,rownames(top),sep="",collapse=NULL)
newURLs<-unlist(newURLs)
acc_table<-aafTable("SNP",items=list(newURLs),colnames="SNP")


chr_table <-aafTable("Chr"=as.character(top[,"CHR"]))
loc_table <-aafTable("Posn"=as.character(top[,"posn"]))
UK_table <-aafTable("UK"=signif(as.numeric(top[,"UK"]),3 ))
UK_US_table<- aafTable("UK_US"=signif(as.numeric(top[,"UK_US"]),3 ))
US_table<- aafTable("US"=signif(as.numeric(top[,"US"]),3 ))

data_table <- merge(acc_table,chr_table)
data_table <- merge(data_table,loc_table)
data_table <- merge(data_table,UK_table)
data_table <- merge(data_table,UK_US_table)
data_table <- merge(data_table,US_table)
saveHTML(data_table,filename,title)
saveText(data_table, sub("html","txt",filename) )


## second/third time
############ a table  </td>


results<- common_all
filename<- "CAR_7-CAR interestion KER_7-KER .html"
title <- "CAR_7-CAR interestion KER_7-KER "


newURLs<-lapply(results[,"swiss_prot"],function(x) paste('<a href="http://au.expasy.org/uniprot/',x,'">', sep = "", collapse = NULL) )
newURLs<-paste(newURLs,results[,"acc"],sep="",collapse=NULL)
newURLs<-unlist(newURLs)
acc_table<-aafTable("Accession",items=list(newURLs),colnames="Accession")

fold_table <-aafTable("fold_ker"=signif(as.numeric(results[,"logFC_ker"]),3),signed=TRUE  )
pval_table<- aafTable("p-val_ker"=signif(as.numeric(results[,"p-val_ker"]),3 ))
Bval_table<- aafTable("B-val_ker"=signif(exp(as.numeric(results[,"B_ker"])),3 ))



 ### first time
data_table <- merge(acc_table,fold_table)
data_table <- merge(data_table,pval_table)
data_table <- merge(data_table,Bval_table)

fold_table <-aafTable("fold_car"=signif(as.numeric(results[,"logFC_car"]),3),signed=TRUE  )
pval_table<- aafTable("p-val_car"=signif(as.numeric(results[,"p-val_car"]),3 ))
Bval_table<- aafTable("B-val_car"=signif(exp(as.numeric(results[,"B_car"])),3 ))

data_table <- merge(data_table,fold_table)
data_table <- merge(data_table,pval_table)
data_table <- merge(data_table,Bval_table)

saveHTML(data_table,filename,title)








    <table border="2">
<tr>
<th>Accession</th>
<th>UK</th>
<th>UK_US</th>
<th>US</th>
</tr>
<tr>
<td class="character"><a href="http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=785702">RS785702</td>
<td class="numeric">3.22e-253</td>
<td class="numeric">0.00463</td>
<td class="numeric">0.00617</td>
</tr>
<tr>
<td class="character"><a href="http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=258415">RS258415</td>
<td class="numeric">1.11e-199</td>
<td class="numeric">NA</td>
<td class="numeric">NA</td>
</tr>
<tr>
<td class="character"><a href="http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=7006137">RS7006137</td>
<td class="numeric">4.77e-197</td>
<td class="numeric">5.18e-205</td>
<td class="numeric">9.32e-119</td>
</tr>
<tr>


















6p.arm.posns<-snp.got[,"CHR"]==6 & snp.got[,"posn"]<60000000



#snp.ann<-scan("NCBI126_from_BC_sm.txt",what=c(character(1),numeric(2)),skip=1,sep="\t",fill=TRUE,na.strings="")

#snp.ann<-scan("test.txt",what=c(character(1),numeric(2)),skip=1,nlines=5,sep=",",fill=TRUE,na.strings="")
num.snps<-length(snp.ann)/(3)
dim(snp.ann)<-c(3,num.snps)
snp.ann<-t(snp.ann)
rownames(snp.ann)<-snp.ann[,1]
snp.ann<-snp.ann[,-1]
colnames(snp.ann)<-c("chromo","posn")
snp.got<-snp.ann[rownames(UK),]

dim(snp.got)    # check not missing values
dim(UK)

snp.got<-as.numeric(snp.got)
6p.arm.posns<-snp.got[,"chromo"]==6 & snp.got[,"posn"]>60000000   # posn of 6p arm in got snps

6p.arm<-rownames(snp.got)[6p.arm.posns] # 6p-arm RS snps

 ### remove 6p from all arrahes
all<-UK
orders<-sort(all[,"P"],index.return=TRUE)$ix

all<-all[orders,]
all




#############################################################
######################## START ########################
basfi<-read.delim("WTCCC BASFI.dat",header=T,sep="\t")
basfi[,1]<-gsub("RS","rs",basfi[,1])
rownames(basfi)<-basfi[,1]

basfi[,1]<-NA
basfi[,3]<-NA

### Below do to with R
library(biomaRt)
listMarts()           # get the data set info
test=useMart("snp")
listDatasets(test)
######################
human_snp_126 = useMart("snp",dataset="hsapiens_snp")

snp_ann<-getBM(attributes=c("refsnp_id","chr_name","chrom_start","chrom_strand"),
filters=c("refsnp"),
values=c(rs.snps), mart=human_snp_126)

## This contains some info of different contigs that I want to remove.
snp_ann[,5]<-as.numeric(snp_ann[,"chr_name"]) < 23 | snp_ann[,"chr_name"]=="X" |  snp_ann[,"chr_name"]=="Y"
snp_ann<-snp_ann[ as.numeric(snp_ann[,"chr_name"]) < 23 | snp_ann[,"chr_name"]=="X" |  snp_ann[,"chr_name"]=="Y",]
snp_ann<-snp_ann[!is.na(snp_ann[,1]),]
##### snp_ann cleaned up
#tapply(snp_ann[,1],snp_ann[,1],length)

rownames(snp_ann)<-snp_ann[,1]

missing<-setdiff(rownames(basfi),rownames(snp_ann))



# works but is slow and unreliable...... direct interrigation of web site.
missing_ann<-data.frame( name=missing,
                         chromo=NA,
                         chr_posn=NA)
rs.snps<-as.numeric(gsub("rs","",missing_ann[,"name"]))
for(i in 1:length(rs.snps)) {
url<-sprintf('http://gvs.gs.washington.edu/GVS/PopStatsServlet?searchBy=dbsnp+rsID&target=%d&upstreamSize=0&downstreamSize=0&x=&y=',rs.snps[i])
results<-htmlTreeParse(url)
test<-unlist(results$children$html[2])
line<-test[grep("chromoStart",test)]
if(length(line)==0){next}
line2<-strsplit(line,"&")
#missing_ann[i,1]<-missing[i]
missing_ann[i,2]<-as.numeric(strsplit(line2[[1]][1],"=")[[1]][2])  # chromosome
missing_ann[i,3]<-as.numeric(strsplit(line2[[1]][2],"=")[[1]][2])  # chromosome posn
}

lost.snps<-missing_ann[is.na(missing_ann[,"chromo"] ),]
missing_ann<-missing_ann[!is.na(missing_ann[,"chromo"] ),]

snp.length<-length(snp_ann[,1])         #undo: snp_ann<-snp_ann[1:snp.length,]
for(i in 1:length(missing_ann[,1])) {
snp_ann[i+snp.length,"refsnp_id"]<-as.character(missing_ann[i,"name"])
snp_ann[i+snp.length,"chr_name"]<-as.character(missing_ann[i,"chromo"])
snp_ann[i+snp.length,"chrom_start"]<-as.numeric(missing_ann[i,"chr_posn"])}

rownames(snp_ann)<-snp_ann[,1]

# add the data to snp_annotation:
snp_ann[rownames(snp_ann),4]<-basfi[rownames(snp_ann),"P"]
reorder<-order(snp_ann[,"chr_name"],snp_ann[,"chrom_start"])
snp_ann<-snp_ann[reorder,]
###################### Annotation file constructed #############
################################################################
################# Write the output file #######################
snp_ann[,5]<-snp_ann[,"chrom_start"]-c(1,snp_ann[1:(length(snp_ann[,"chrom_start"])-1),"chrom_start"]) -1
snp_ann[,6]<-c(snp_ann[2:length(snp_ann[,"chrom_start"]),"chrom_start"],snp_ann[length(snp_ann[,"chrom_start"]),"chrom_start"])-snp_ann[,"chrom_start"] -1

outfile<-file(description="wiggle_p.txt", open="at")
text<-"browser position chr1:1-247,249,719\nbrowser hide all\nbrowser dense refGene encodeRegions snp126 multiz28way"
writeLines(text, con = outfile, sep = "\n")
text<-'track type=wiggle_0 name="variableStep2" description="variableStep format" visibility=full autoScale=off viewLimits=-0.5:0.95 color=255,200,0 yLineMark=0.0 yLineOnOff=on priority=10'
writeLines(text, con = outfile, sep = "\n")   # flush(outfile) need to see text
#snp_ann[,"chr_name"]<-paste("chr",snp_ann[,"chr_name"],sep="")
current_chr<- snp_ann[1,"chr_name"]
#width<-1
for(i in 1:dim(snp_ann)[1]) {
width_high<-1
width_low<-1
if(current_chr==snp_ann[i,"chr_name"]){width_high=min(width,snp_ann[i,6])
                                         width_low=min(width,snp_ann[i,5])  }
text<-paste( paste("chr",snp_ann[i,"chr_name"],sep="") ,snp_ann[i,"chrom_start"]-width_low, snp_ann[i,"chrom_start"]+width_high,snp_ann[i,4],sep=" ")
writeLines(text, con = outfile, sep = "\n")
}
close(outfile)









track type=wiggle_0 name="variableStep" description="variableStep format" visibility=full autoScale=off viewLimits=0.0:1.0 color=100,50,0 useScore=1 altColor=0,100,200 yLineMark=0.05 yLineOnOff=on priority=10
browser position chr1:1-59311000
browser hide all
browser dense refGene encodeRegions snp126 multiz28way 





write(x, file = "wiggle_p.txt",
      ncolumns = if(is.character(x)) 1 else 5,
      append = TRUE, sep = " ")

  snp_ann[,5]<-snp_ann[,"chrom_start"]-c(1,snp_ann[1:(length(snp_ann[,"chrom_start"])-1),"chrom_start"])
  snp_ann[,6]<-c(snp_ann[2:length(snp_ann[,"chrom_start"]),"chrom_start"],snp_ann[length(snp_ann[,"chrom_start"]),"chrom_start"])-snp_ann[,"chrom_start"]

> test1<-test[-1]
> test2<-test[-length(test)]
> test<-test1-test2
  test<-c(1,test)
> snp_ann[,5]<-test
 > snp_ann[,5]<-test
## works tto get get more missing than when using the above code
# missing_ann<-data.frame( name=missing,
#                          chromo=NA,
#                          chr_posn=NA)
# #rs.snps<-as.numeric(gsub("rs","",missing_ann[,"name"]))
# for(i in 1:length(missing_ann[,1])) {                        # but loop dies if it can't find the SNP
# asnp<-try(SNPinfo(missing_ann[i,1]),silent=TRUE)
# if(asnp[1]!=missing_ann[i,1]){next}
# missing_ann[i,2]<-gsub("chr","",asnp[2])  # chromosome
# missing_ann[i,3]<-as.numeric(asnp[3])  # chromosome posn
# }

    test[rownames(test),5]<-basfi[rownames(snp_ann),"P"]


url<-'http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=1050150#map'
#<a href="Navigate?chrNumber=1&chromoStart=1110294&chromoStop=1110294&comingFrom=mapButton">    # return from line above


output<-basfi[,c("SNP","P")]
output<-basfi[,c("SNP")]
write.table(output,"the_snps.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

header<-'track type=wiggle_0 name="variableStep" description="variableStep format" visibility=full autoScale=off viewLimits=0.0:25.0 color=255,200,0 yLineMark=11.76 yLineOnOff=on priority=10'



























url<-'http://genome.ucsc.edu/cgi-bin/das/hg16/features?segment=1:1,100000;type=snp'    # works
url<-'http://genome.ucsc.edu/cgi-bin/das/hg16/features?segment=1:1,10000;type=snp'     #works
results <- xmlRoot(xmlTreeParse(url))           # but the passig looks terrible

http://genome.ucsc.edu/cgi-bin/das/hg18/entry_points
http://genome.ucsc.edu/cgi-bin/das/hg18/types


 library(ACME)
samps <- data.frame(SampleID = 1:3, Sample = paste("Sample",
+ 1:3))
 dat <- matrix(rnorm(90000), nc = 3)
 colnames(dat) <- 1:3
 pos <- rep(seq(1, 10000 * 100, 100), 3)
 chrom <- rep(paste("chr", 1:3, sep = ""), each = 10000)
 annot <- data.frame(Chromosome = chrom, Location = pos)
 a <- new("aGFF", data = dat, annotation = annot, samples = samps)

 calc <- do.aGFF.calc(a, window = 1000, thresh = 0.95)
























1[CHR] AND 67394757:67508250[CHRPOS] AND human[Organism]

Use SNP and then export usng dbSNP batch report (make sure resitrict to species)

Cut and paste the llist or use the chromosome report


snps<-read.delim("SNP_10kb.txt",header=F,sep="\t")
snps<-read.delim("SNP_gene.txt",header=F,sep="\t")


snp_gene<-as.vector(snps[,1]) snp_gene
all_gene<-SNP_ann[snp_gene,]
 founds<-!is.na(all_gene[,"MapInfo"])
found_gene<-all_gene[founds,]
true_names<-snp_gene[founds]
found_genes<-found_gene[rownames(found_gene)==true_names,]


snps<-read.delim("SNP_10kb.txt",header=F,sep="\t",fill=TRUE)


snp_10kb<-as.vector(snps_10kb[,1]) # cause its a data frame
all<-SNP_ann[snp_10kb,]

found_10kb<-all[!is.na(all[,"MapInfo"]),]



  snps_nei<-read.delim("neighbours.txt",header=F,sep="\t",fill=TRUE)
  snp_nei<-snps_nei[,"V4"]
SNP_ann[as.vector(snp_nei),]

rs790631    1 67388943 in intron 6
SNP_ann[c(snps),]







path<-"C:\\Research\\Matt Brown\\Linkage Chip investigation"
datadir <- system.file(path)
readLines(paste(path, "forpaul.csv", sep = "\\"))

sample_descriptions<-read.delim("forpaul.csv",header=T,skip=7,sep=",",fill=TRUE)
f_report<-read.delim("ForPaul-GTcomparison_FinalReport.csv",header=F,skip=10,sep="\t")
colnames(f_report)<-c("SNP_name","Sample_ID","GC_Score","Allele1","Allele2","GT_score","X_Raw","Y_Raw")
samples<-sample_descriptions[,1]


posns<-grep(samples[1],f_report[,2])
> length(posns)
[1] 317503

locations<-seq(0,length(posns)*length(samples),length(posns))
locations[1]<-1

#alist_of_locations<-apply(as.matrix(samples),1,function(x,column) grepalist(x,column),column=f_report[,2])

#f_report[(locations[site]-3):(locations[site]+3),]



 site<-1
calls<-f_report[(locations[site]+1):locations[site+1],]
rownames(calls)<-f_report[(locations[site]+1):locations[site+1],1]
colnames(calls)<-samples # just happens to be the same length
for (site in 1:length(samples)){ calls[,site]<-paste(f_report[(locations[site]+1):locations[site+1],4],f_report[(locations[site]+1):locations[site+1],5],sep="") }
 calls<-calls[,sort(colnames(calls))]   # order together

rbind(f_report[(locations[site]):(locations[site]+5),c(1,2,4,5)] ,  f_report[(locations[site+1]):(locations[site+1]+5),c(1,2,4,5)] #check


## get the SNP annotations
SNP_ann<-read.delim("ILMN_HumanHap300_SNPlist.txt",header=T,sep="\t")

rownames(SNP_ann)<-SNP_ann[,1]
SNP_ann[1:5,]
#                Name GenomeBuild Chr MapInfo
# rs3934834 rs3934834          35   1 1045729
# rs3737728 rs3737728          35   1 1061338
# rs6687776 rs6687776          35   1 1070488
# rs9651273 rs9651273          35   1 1071463
# rs4970405 rs4970405          35   1 1088878
SNP_ann<-SNP_ann[,c(-1,-2,-4)]
SNP_ann[1:5,]

### order SNP_ann the same as in calls
 SNP_ann<-SNP_ann[rownames(calls),] ## takes a LONG time TO RUN
write.table(SNP_ann,"SNP_ann.txt",sep="\t",col.names=TRUE,row.names=TRUE)

 ## Make comparisons
GB_comp<-calls[,1:4]

 colnames(GB_comp)<-c("HI_S_12454","LO_S_2516","UK330.3","UK366.3")

for (site in 1:4){GB_comp[,site]<-calls[,(2*site-1)]!=calls[,(2*site)] }


> apply(GB_comp,2,sum)
HI_S_12454  LO_S_2516    UK330.3    UK366.3
     20703      13318      35949      33367

 ## Different calls between GOOD and BAD arrays
>  100-100*apply(GB_comp,2,sum)/dim(GB_comp)[1]
HI_S_12454  LO_S_2516    UK330.3    UK366.3
  93.47943   95.80539   88.67759   89.49081 

b_UK330<-calls[GB_comp[,3],c(5,6)]



  control_sig<-tapply(temp[,colnames(temp)== "Signal.Median.Control" ],acc,mean)

SNPS_diff<-apply(GB_comp,2,function(x) rownames(x)[x])


### count the number not genotyped
missing<-apply(calls,2,function(x) x=="--")
100-100*apply(missing,2,sum)/dim(missing)[1]

> apply(missing,2,sum)
 HI_S_12454_bad HI_S_12454_good   LO_S_2516_bad  LO_S_2516_good     UK330.3_bad    UK330.3_good
          20064            3324           14158            3462           36566            3706
    UK366.3_bad    UK366.3_good
          32345            2947


# check sexes
SNP_male<-SNP_ann[,1]=="XY"
 > calls[SNP_male,]
           HI_S_12454_bad HI_S_12454_good LO_S_2516_bad LO_S_2516_good UK330.3_bad UK330.3_good UK366.3_bad
rs17842869             AA              AA            BB             BB          BB           BB          BB
rs4933045              AB              BB            AB             AB          --           BB          AB
           UK366.3_good
rs17842869           BB
rs4933045            AB
>
  par(mfrow=c(1,1))
toplot<-apply(missing,2,sum)
barplot(toplot,xlab="Sample",ylab="Number not genotyped",col=rainbow(20),cex.names=0.6,space=0,las=2)
title("Number not genotyped in each sample")
grid(ny=7,nx=0)

toplot<-100*apply(GB_comp,2,sum)/dim(GB_comp)[1]
barplot(toplot,xlab="Sample",ylab="Good/Bad genotype differences",col=rainbow(4),cex.names=0.75,space=0,las=2,ylim=c(0,12))
title("% of genotype differences")
grid(ny=12,nx=0)


labs<-c(1:22,"X","XY")



for(i in 1:8) {
toplot<-100*tapply(missing[,i],chrs,sum)/totals
barplot(toplot[labs],xlab="Chromosome",ylab="% missing per chr",col=rainbow(20),cex.names=0.75,space=0,las=2,ylim=c(0,10))
title(as.character(colnames(missing)[i]))
grid(ny=5,nx=0)      }

"UK366.3_good"
UK330.3_good
LO_S_2516_good
HI_S_12454_good
UK366.3_bad
UK330.3_bad
LO_S_2516_bad
HI_S_12454_bad

SNPdata <- read.SnpSetIllumina(paste(path, "forpaul.csv" ,sep = "\\"), reportpath=path,reportfile="ForPaul-GTcomparison_FinalReport.csv",verbose=TRUE)
> SNPdata

test<-calls[SNP_X,site]
test.f<-factor(test)
tapply(test,test.f,length)

"UK366.3_good"
  --   AA   AB   BB
 359 4131    5 4676 

 "UK366.3_bad"
  --   AA   AB   BB
1445 3537  147 4042

"HI_S_12454_good"
  --   AA   AB   BB
  93 2821 3205 3052

par(mfrow=c(4,1))
  for(i in 1:4) {
toplot<-100*tapply(GB_comp[,i],chrs,sum)/totals
barplot(toplot[labs],xlab="Chromosome",ylab="% different per chr",col=rainbow(20),cex.names=0.75,space=0,las=2,ylim=c(0,15))
title(as.character(colnames(GB_comp)[i]))
grid(ny=5,nx=0)      }

GTcalls<-f_report[(locations[site]+1):locations[site+1],]
 for (site in 1:length(samples)){ GTcalls[,site]<-f_report[(locations[site]+1):locations[site+1],6]}
 colnames(GTcalls)<-samples
 GTcalls<-GTcalls[,sort(colnames(GTcalls))]   # order together


 par(mfrow=c(1,2))
  test<-GTcalls[SNP_X,2]
 hist(test,xlab="X Chromosome: GT Score",ylab="Number",main="",col=rainbow(24))
title(as.character(colnames(GTcalls)[2]))





######## QC control in X chromosome in MEM : must be heterozygous

GTcallX<-GTcalls[SNP_X,]
> dim(GTcallX)
[1] 9171    8
callsX<-calls[SNP_X,]
> dim(callsX)
[1] 9171    8

callsX_AB<-apply(callsX,2,function(x) x=="AB")
apply(callsX_AB,2,sum)

#How many:
> apply(callsX_AB,2,sum)
 HI_S_12454_bad HI_S_12454_good   LO_S_2516_bad  LO_S_2516_good     UK330.3_bad    UK330.3_good 
           3044            3205            3133            3255             119               7 
    UK366.3_bad    UK366.3_good
            147               5

# just males       only UK samples    :
callsX<-callsX[,5:8]
callsX_AB<-callsX_AB[,5:8]
GTcallX<-GTcallX[,5:8]

# number of wrounf calls that is AB call make when should be homo
apply(callsX_AB,2,sum)
 UK330.3_bad UK330.3_good  UK366.3_bad UK366.3_good 
         119            7          147            5

> GTcallX[callsX_AB[,2],2]
 "UK330.3_good" [1] 0.8087 0.8123 0.6839 0.6719 0.6866 0.8173 0.7778 mean:  0.7512143
> GTcallX[callsX_AB[,4],4]
"UK366.3_good": [1] 0.8065 0.8068 0.7885 0.7046 0.7983  mean: 0.78094

 "UK330.3_bad" mean: 0.8291639
 "UK366.3_bad" mean: 0.8234279

rbind(GTcallX[callX_AB[,1],1],callX[callX_AB[,1],1])
rbind(GTcallX[callsX_AB[,1],1],callsX[callsX_AB[,1],1])



> apply(callsX_AB,2,sum)



 ######## GC score

 GCcalls<-f_report[(locations[site]+1):locations[site+1],]
 for (site in 1:length(samples)){ GCcalls[,site]<-f_report[(locations[site]+1):locations[site+1],3]}
 colnames(GCcalls)<-samples
 rownames(GCcalls)<-f_report[(locations[site]+1):locations[site+1],1]
 GCcalls<-GCcalls[,sort(colnames(GCcalls))]   # order together

 GCcallX<-GCcalls[SNP_X,]
 GCcallX<-GCcallX[,5:8]

 GCcallX[callsX_AB[,2],2]
 "UK330.3_good"  [1] 0.7390 0.8099 0.5712 0.2640 0.6136 0.5709 0.7799   mean: 0.6212143
 
  GCcallX[callsX_AB[,4],4]
 "UK366.3_good": 0.8098 0.6018 0.6322 0.6481 0.8129    mean:  0.70096
 

   > mean( GCcallX[callsX_AB[,1],1])
[1] 0.675226
> mean( GCcallX[callsX_AB[,3],3])
[1] 0.6708905
 
     > median( GCcallX[callsX_AB[,1],1])
[1] 0.7157
> median( GCcallX[callsX_AB[,3],3])
[1] 0.7366

### if GC score was < 0.81 then remove.      for any genotype except the empties


passed<-apply(GCcalls[,c(2,4,6,8)],2,function(x)  x>0.8  | x==0.0000)

miss_callsX<-apply(callsX,2,function(x) x=="--")
> apply(miss_callsX,2,sum)     # numners not assigned
 UK330.3_bad UK330.3_good  UK366.3_bad UK366.3_good 
        2756          449         1445          359

 total number f SNPs in X: 9171 ### 9171-449=total called 7/8722 * 100 _>error rate about 0.08%
 
 ## average all (9171-449)+ (9171-359)= toal good calls  =17534
 ## total errors 5+7 = 0.07% error rate.

 #the number lost:        SNP_
  passed<-apply(GCcalls[,c(2,4,6,8)],2,function(x)  x>0.8  | x==0.0000)
> apply(passed,2,sum)/dim(passed)[1]
  (1-apply(passed,2,sum)/dim(passed)[1])*100
  
   HI_S_12454_good  LO_S_2516_good    UK330.3_good    UK366.3_good
       13.05720        13.17342        13.23263        12.96145

   > passed<-apply(GCcalls[,c(2,4,6,8)],2,function(x)  x>0.71  | x==0.0000)
> (1-apply(passed,2,sum)/dim(passed)[1])*100
HI_S_12454_good  LO_S_2516_good    UK330.3_good    UK366.3_good 
       4.933812        5.012866        5.094755        4.860112


       > (1-apply(passed,2,sum)/dim(passed)[1])*100
HI_S_12454_good  LO_S_2516_good    UK330.3_good    UK366.3_good
       2.944224        3.002176        3.070837        2.866115 
x<- seq(0,1,0.01)
y<-x
for (i in  1:length(x) ){
passed<-apply(GCcalls[,c(2,4,6,8)],2,function(x,c1,c2)  x>c1[c2]  | x==0.0000,c1=x,c2=i)
y[i]<-mean((1-apply(passed,2,sum)/dim(passed)[1])*100 )
}
## basically in order to remove no more than 0.08% of calls th GC Score can be no higher than 0.48. 
## This would have only removed one bad genotyping error  from the AB X-chromosome test in males.







