#!/bin/bash

###To run script call ./QC_cleaning.bash path_to_input_plink_binary_files output_prefix 
# QC_cleaning.bash /40/AD/GWAS_data/Source_Plink/202102_DystoniaCoalitionExternal/01_pre_QC/emory/05_DYS_WT_919samples_SNPs_updated_alleles_no_unkown_pos_no_dups_cleaned DYS_WT_919samples

plink_data=$1
name=$2
Y_code=$3
HWE=$4

rm report_2_QC.txt
# remve the Y chr for missingness
echo "Removing Y chromosome ..." >> report_2_QC.txt
plink1.9 --bfile ${plink_data} --allow-no-sex --keep-allele-order --not-chr 24 --make-bed --out 06_${name}_preQC_cleaned_no_Ychr
echo "number of variants after removing Y chromosome: " >> report_2_QC.txt
wc -l 06_${name}_preQC_cleaned_no_Ychr.bim >> report_2_QC.txt



#also generate a file that only has chrY if we need to add them back at some point
echo "Generating a Y chromosome only dataset if sex chr present..." >> report_2_QC.txt
plink1.9 --bfile ${plink_data} --allow-no-sex --keep-allele-order --chr ${Y_code} --make-bed --out 06_${name}_preQC_cleaned_ONLY_Ychr
echo "Number of variants in Y chr: " >> report_2_QC.txt
wc -l 06_${name}_preQC_cleaned_ONLY_Ychr.bim >> report_2_QC.txt







											         ############################################################
							                                         ## 1) Investigate missingness per individual and per SNP  ##
												 ############################################################



echo "     Filter for missingness " >> report_2_QC.txt

######## Relaxed #######
#always SNPs first, then individuals

# Delete SNPs with missingness >0.2
plink1.9 --bfile 06_${name}_preQC_cleaned_no_Ychr --allow-no-sex --keep-allele-order --geno 0.2 --make-bed --out 07_${name}_preQC_cleaned_no_Ychr_geno0.2
echo "Number of variants after filetering for SNPs with missingness >0.2 : " >> report_2_QC.txt
wc -l 07_${name}_preQC_cleaned_no_Ychr_geno0.2.bim >> report_2_QC.txt

#Delete individuals with missingness >0.2
plink1.9 --bfile 07_${name}_preQC_cleaned_no_Ychr_geno0.2 --allow-no-sex --keep-allele-order --mind 0.2 --make-bed --out 08_${name}_preQC_cleaned_no_Ychr_geno0.2_mind0.2
echo "Number of individuals after filetering for individuals with missingness >0.2 : " >> report_2_QC.txt
wc -l 08_${name}_preQC_cleaned_no_Ychr_geno0.2_mind0.2.fam >> report_2_QC.txt



######### stringent ##########
# Delete SNPs with missingness >0.2
plink1.9 --bfile 08_${name}_preQC_cleaned_no_Ychr_geno0.2_mind0.2 --allow-no-sex --keep-allele-order --geno 0.02 --make-bed --out 09_${name}_preQC_cleaned_no_Ychr_geno0.02
echo "Number of variants after filetering for SNPs with missingness >0.02 : " >> report_2_QC.txt
wc -l 09_${name}_preQC_cleaned_no_Ychr_geno0.02.bim >> report_2_QC.txt

#Delete individuals with missingness >0.2
plink1.9 --bfile 09_${name}_preQC_cleaned_no_Ychr_geno0.02 --allow-no-sex --keep-allele-order --mind 0.02 --make-bed --out 10_${name}_preQC_cleaned_no_Ychr_geno0.02_mind0.02
echo "Number of individuals after filetering for individuals with missingness >0.02 : " >> report_2_QC.txt
wc -l 10_${name}_preQC_cleaned_no_Ychr_geno0.02_mind0.02.fam >> report_2_QC.txt



											         ############################################################
							                                         ## 2) Check for sex discrepancy and delete those with discrepancy  ##
												 ############################################################

echo "     Check for sex discrepancy   " >> report_2_QC.txt


if [ ${Y_code} != NA ]
then
# we have to remove the samples (indv) that were removed due to missingness from the Y chr subset before merging. otherwise it will add those samples back
#find the 6 samples that have to be removed.  I would use diff to compare the two files and then output only lines that are in the left file but not in right one. Such lines are flagged by diff with < so it suffices to grep that symbol at the beginning of the line
	diff 06_${name}_preQC_cleaned_ONLY_Ychr.fam 10_${name}_preQC_cleaned_no_Ychr_geno0.02_mind0.02.fam  | grep \^\< > 10_list_samples_to_remove_from_chY.txt
	cut -d " " -f 2,3 10_list_samples_to_remove_from_chY.txt > 10_samples_to_remove.txt

# remove 6 sample from chrY
	plink1.9 --bfile 06_${name}_preQC_cleaned_ONLY_Ychr --allow-no-sex --keep-allele-order --remove 10_samples_to_remove.txt --make-bed --out ONLY_Ychr_preQC_cleaned_min0.02_samples_removed

#merge the chrY with the missingness filtered file
	plink1.9 --bfile 10_${name}_preQC_cleaned_no_Ychr_geno0.02_mind0.02 --allow-no-sex --keep-allele-order --bmerge ONLY_Ychr_preQC_cleaned_min0.02_samples_removed --recode --out 11_${name}_preQC_cleaned_geno0.02_mind0.02_chrY_added
	echo "Number of variants after merging no Y chr data with Y chr data: " >> report_2_QC.txt
	wc -l 11_${name}_preQC_cleaned_geno0.02_mind0.02_chrY_added.bim >> report_2_QC.txt
	echo "Number of individuals after merging no Y chr data with Y chr data: " >> report_2_QC.txt
	wc -l 11_${name}_preQC_cleaned_geno0.02_mind0.02_chrY_added.fam >> report_2_QC.txt

# do the sex check
	plink1.9 --bfile 11_${name}_preQC_cleaned_geno0.02_mind0.02_chrY_added --allow-no-sex --keep-allele-order --check-sex 

	echo "Number of individuals with sex discrepancy: ">> report_2_QC.txt
	grep -c "PROBLEM" plink.sexcheck >> report_2_QC.txt

# Impute
# This imputes the sex based on the genotype information into your data set
	plink1.9 --bfile 11_${name}_preQC_cleaned_geno0.02_mind0.02_chrY_added --allow-no-sex --keep-allele-order --impute-sex --make-bed --out 12_${name}_preQC_cleaned_geno0.02_mind0.02_imputed_sex


### use R scripts to plot gender data:
echo "Plotting sexcheck data in R ... The plots are in files Gender_check.pdf, Men_check.pdf, and Women_check.pdf." >> report_2_QC.txt
path=$(pwd)
Rscript sexcheck_plot.R $path plink.sexcheck


else 
	echo "Can not check for sex discrepancy, NO sex chromosome in the file" >> report_2_QC.txt
fi


											         ############################################################
							                                         ##        3) Filter to keep only AUTOSOMAL SNPs      ##

												  ############################################################

echo "      Filter to keep only AUTOSOMAL SNPs      " >> report_2_QC.txt


if [ ${Y_code} != NA ]
	then
#### WE DO NOT FILTER FOR MAF IN THIS LAB
# Select autosomal SNPs only (i.e., from chromosomes 1 to 22).

#	echo " Number of autosomal variants:" >> report_2_QC.txt
#	awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' 12_${name}_preQC_cleaned_geno0.02_mind0.02_imputed_sex.bim > 12_variants_chr_1_22.txt
#	wc -l 12_variants_chr_1_22.txt >> report_2_QC.txt
	
	plink1.9 --bfile 12_${name}_preQC_cleaned_geno0.02_mind0.02_imputed_sex --allow-no-sex --keep-allele-order --autosome --make-bed --out 13_${name}_preQC_cleaned_geno0.02_mind0.02_Autosomal
	echo "Number of autosomal variants using plink --autosome filter:">> report_2_QC.txt
	wc -l 13_${name}_preQC_cleaned_geno0.02_mind0.02_Autosomal.bim >> report_2_QC.txt
	
	plink1.9 --bfile 12_${name}_preQC_cleaned_geno0.02_mind0.02_imputed_sex --chr 23,24,25,26 --allow-no-sex --keep-allele-order --make-bed --out 13_${name}_DC_preQC_cleaned_geno0.02_mind0.02_Not_autosomal
	echo "Number of non-autosomal variants after Plink filter: ">> report_2_QC.txt
	wc -l 13_${name}_DC_preQC_cleaned_geno0.02_mind0.02_Not_autosomal.bim >> report_2_QC.txt

else
	plink1.9 --bfile 10_${name}_preQC_cleaned_no_Ychr_geno0.02_mind0.02 --allow-no-sex --keep-allele-order --autosome --make-bed --out 13_${name}_preQC_cleaned_geno0.02_mind0.02_Autosomal
	echo "Number of autosomal variants using plink --autosome filter: ">> report_2_QC.txt
	wc -l 13_${name}_preQC_cleaned_geno0.02_mind0.02_Autosomal.bim >> report_2_QC.txt

	plink1.9 --bfile 10_${name}_preQC_cleaned_no_Ychr_geno0.02_mind0.02 --chr 23,24,25,26 --allow-no-sex --keep-allele-order --make-bed --out 13_${name}_DC_preQC_cleaned_geno0.02_mind0.02_Not_autosomal
	echo "Number of non-autosomal variants after Plink filter: ">> report_2_QC.txt
	wc -l 13_${name}_DC_preQC_cleaned_geno0.02_mind0.02_Not_autosomal.bim >> report_2_QC.txt

fi


											         ############################################################
							                                         ##         4) Delete SNPs which are not in HWE            ##
												 ############################################################

echo "      Check for HWE      " >> report_2_QC.txt
###use 1e-6 for both case and controls
#--hwe's 'midp' modifier applies the mid-p adjustment described in Graffelman J, Moreno V (2013) The mid p-value in exact tests for Hardy-Weinberg equilibrium. The mid-p adjustment tends to bring the null rejection rate in line with the nominal p-value, and also reduces the filter's tendency to favor retention of variants with missing data. We recommend its use.

plink1.9 --bfile 13_${name}_preQC_cleaned_geno0.02_mind0.02_Autosomal --hwe ${HWE} midp --keep-allele-order --allow-no-sex --make-bed --out 14_${name}_DC_preQC_cleaned_geno0.02_mind0.02_Autosomal_hwe
echo "Number of variants after filtering for HWE : " >> report_2_QC.txt
wc -l 14_${name}_DC_preQC_cleaned_geno0.02_mind0.02_Autosomal_hwe.bim >> report_2_QC.txt

											         ############################################################
							                                         ##         5)        Check for heterozygosity           ##
												 ############################################################

echo "      Check for heterozygosity      " >> report_2_QC.txt


# extract the list of SNPs that are not in LD
plink1.9 --bfile 14_${name}_DC_preQC_cleaned_geno0.02_mind0.02_Autosomal_hwe --keep-allele-order --allow-no-sex --indep-pairwise 50 5 0.2 --out 14_indepSNP
#indepSNP.prune.in contains the list of SNPs that are not in LD and in selected inversion regions
echo "Number of independent variants after pruning (--indep-pairwise 50 5 0.2) for hetro and IBD analyses:" >> report_2_QC.txt
wc -l 14_indepSNP.prune.in >> report_2_QC.txt


plink1.9  --bfile 14_${name}_DC_preQC_cleaned_geno0.02_mind0.02_Autosomal_hwe --keep-allele-order --extract 14_indepSNP.prune.in --het --out 14_R_check_hetero

path=$(pwd)
echo "plotting heterozygosity rate and identifying outliers usig check_heterozygosity_rate.R. The plot is in heterozygosity.pdf and the outliers are in fail-het-qc.txt. Samples whose het rate varies 3sd from the mean are considered outliers. " >> report_2_QC.txt

Rscript check_heterozygosity_rate.R $path 

# selecting only the first two columns.
awk '{print$1, $2}' fail-het-qc.txt > het_fail_ind.txt

# Remove heterozygosity rate outliers. Note that the outliers are removed from the file that is filtered for HWE and only has autosomal chr(.
plink1.9 --bfile 14_${name}_DC_preQC_cleaned_geno0.02_mind0.02_Autosomal_hwe --allow-no-sex --keep-allele-order --remove het_fail_ind.txt --make-bed --out 15_${name}_preQC_cleaned_geno0.02_mind0.02_Autosomal_hwe_noHet
echo "Number of individuals after filtering for heterozygosity : ">> report_2_QC.txt
wc -l 15_${name}_preQC_cleaned_geno0.02_mind0.02_Autosomal_hwe_noHet.fam >> report_2_QC.txt



											         ############################################################
							                                         ##              6)cryptic relatedness                     ##
												 ############################################################

# ONLY USE VARIANTS ON AUTOSOMAL CHR

# Check for relationships between individuals with a pihat > 0.2. (-- genome gives IBD based on IBS)
echo "     Check for cryptic relatedness      " >> report_2_QC.txt
echo "Extract individuals with pihat > 0.2 (--min 0.2) after filtering for --geno 0.0001 --hwe 0.05 --maf 0.15" >> report_2_QC.txt
plink1.9 --bfile 15_${name}_preQC_cleaned_geno0.02_mind0.02_Autosomal_hwe_noHet --allow-no-sex --keep-allele-order --extract 14_indepSNP.prune.in --genome --min 0.2 --geno 0.0001 --hwe 0.05 --maf 0.15 --out pihat_min0.2



###############plot the data
awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome > zoom_pihat.genome

echo "number of sample pairs with PI-HAT > 0.2 :" >> report_2_QC.txt
wc -l pihat_min0.2.genome >> report_2_QC.txt
#37

echo "number of sample pairs with PI-HAT > 0.9 :" >> report_2_QC.txt
wc -l zoom_pihat.genome >> report_2_QC.txt


path=$(pwd)
Rscript Relatedness.R $path


#stats- you can run these after you get the pihat freq table and plots. you have to see what are the break points. usually there should be one above 0.8 which are duplicates, and another cluster between 0.45-0.55 but the exact number could be retrived from the break points in the output table from R script in bash
#are there any duplicate IDs (expected duplicats)
#cut -f 1 pihat_min0.2.genome | sort | uniq -D >> report_2_QC.txt


#how many identical samples/duplicate (pihat =1)
#awk '($10 >0.98 ) {count++ } END { print count }' pihat_min0.2.genome
#awk '{ if ($10 >= 0.98) print $1,$2,$3,$4 }' pihat_min0.2.genome > identical_samples_pihat_0.98.txt
#wc -l identical_samples_pihat_0.98.txt


# how many indv are first degree relatives (pihat =0.45-0.55)
#get the list
#awk '{ if ($10 >= 0.45 && $10 <= 0.55) print $1,$2,$3,$4 }' pihat_min0.2.genome > first_degree_relatives_pihat_0.45_0.55.txt
#wc -l first_degree_relatives_pihat_0.45_0.55.txt




											         ############################################################
							                                         ##      7)Recombine all chromosomes if applicable         ##
												 ############################################################

echo "       Recombine all chromosomes if applicable        "
# we need to match the indv before merging the variants
# get the list of final indv for plink filtering
awk '{print$1, $2}' 15_${name}_preQC_cleaned_geno0.02_mind0.02_Autosomal_hwe_noHet.fam > 15_list_final_indv.txt


#keep only indv in list from non-autosomal data
plink1.9 --bfile 13_${name}_DC_preQC_cleaned_geno0.02_mind0.02_Not_autosomal --allow-no-sex --keep-allele-order --keep 15_list_final_indv.txt --make-bed --out 13_${name}_preQC_cleaned_geno0.02_mind0.02_Not_autosomal_samples_removed

plink1.9 --bfile 15_${name}_preQC_cleaned_geno0.02_mind0.02_Autosomal_hwe_noHet --allow-no-sex --keep-allele-order --bmerge 13_${name}_preQC_cleaned_geno0.02_mind0.02_Not_autosomal_samples_removed --recode --out 16_${name}_preQC_cleaned_geno0.02_mind0.02_hwe_noHet_chr_recombined
echo "Number of variants in recombined autosomal and non-autosomal files (if applicable): " >> report_2_QC.txt
wc -l 16_${name}_preQC_cleaned_geno0.02_mind0.02_hwe_noHet_chr_recombined.bim >> report_2_QC.txt

echo "Number of individuals in recombined autosomal and non-autosomal files (if applicable): " >> report_2_QC.txt
wc -l 16_${name}_preQC_cleaned_geno0.02_mind0.02_hwe_noHet_chr_recombined.fam >> report_2_QC.txt



