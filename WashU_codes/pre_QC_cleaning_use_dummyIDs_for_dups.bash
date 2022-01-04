#!/bin/bash

###To run script call ./pre_QC_cleaning.bash path_to_input_plink_binary_files path_to_manifest_file

plink_data=$1
name=$2
manifest=$3
dup_cutoff=$4

rm report_1_pre_QC.txt

# 1) check the phenotype col
echo "num of indv/ pheno:" > report_1_pre_QC.txt
awk '{print $6}' ${plink_data}.fam | sort | uniq -c  >> report_1_pre_QC.txt

#1.5) number of variants
echo "num of variants:" >> report_1_pre_QC.txt
wc -l ${plink_data}.bim  >> report_1_pre_QC.txt


#2) check if sex code is missing
echo "number/ sex code(0=missing, 1=male, 2=female):" >> report_1_pre_QC.txt
awk '{print $5}' ${plink_data}.fam | sort | uniq -c >> report_1_pre_QC.txt

#3) check if we have duplicate sample IDs
echo "number of duplicate sample IDs:" >> report_1_pre_QC.txt
cut -f 1 ${plink_data}.fam | sort | uniq -D >> report_1_pre_QC.txt

#4) check what chr are in the data
echo " number of variants on each chr:" >> report_1_pre_QC.txt
awk '{print $1}' ${plink_data}.bim | sort | uniq -c >> report_1_pre_QC.txt


							#################################################
							############ Update alleles      ##############
							#################################################


echo '                THIS IS THE START OF UPDATE ALLELE PROCESS         '

#Find SNPs that have both alleles as 0
awk '{if($5 =="0"  && $6 == "0") print $2}' ${plink_data}.bim > SNPs_double0_remove.txt
echo "num of variants with both alleles as missing:" >> report_1_pre_QC.txt
wc -l SNPs_double0_remove.txt >> report_1_pre_QC.txt



#exclude the SNPs with both alleles as 0
plink1.9 --bfile  ${plink_data} --keep-allele-order --allow-no-sex --exclude SNPs_double0_remove.txt --make-bed --out 01_${name}_SNPs_double0_removed
echo "number of variants after filtering those with both missing allele:" >> report_1_pre_QC.txt
wc -l 01_${name}_SNPs_double0_removed.bim >> report_1_pre_QC.txt



############# Need to take care of all variants with one zero allele before looking at duplicates

awk '{if($5 == "0" && $6 != "0") print $0}' 01_${name}_SNPs_double0_removed.bim  > SNPs_one_allele0.txt
echo "num of variants with one missing allele:" >> report_1_pre_QC.txt
wc -l SNPs_one_allele0.txt >> report_1_pre_QC.txt 


# check if we have allele zero in col 6 (allele2)
awk '{if($5 != "0" && $6 == "0") print $0}' 01_${name}_SNPs_double0_removed.bim > test.txt
echo "num of variants with missing allele in col6:" >> report_1_pre_QC.txt
wc -l test.txt >> report_1_pre_QC.txt


### use R scripts to update alleles:
echo "Updating the allele info with R script. Check the R_report_allele_update.txt file for report." >> report_1_pre_QC.txt

rm R_report_allele_update.txt
path=$(pwd)
Rscript allele_update.R $path SNPs_one_allele0.txt $manifest


echo "num of variants with one missing allele in array that do not exist in manifes:" >> report_1_pre_QC.txt
wc -l variants_not_in_manifest_to_remove.txt >> report_1_pre_QC.txt

##  update the alleles in plink
plink1.9 --bfile 01_${name}_SNPs_double0_removed --update-alleles allele_change.txt  --allow-no-sex --keep-allele-order --make-bed --out 02_${name}_SNPs_updated_alleles

## remove variants with one missing allele and not in manifest file
plink1.9 --bfile 02_${name}_SNPs_updated_alleles --allow-no-sex --keep-allele-order  --exclude variants_not_in_manifest_to_remove.txt --make-bed --out 02_${name}_SNPs_updated_alleles2

# check again the zero alleles
awk '{if($5 == "0" && $6 != "0") print $0}' 02_${name}_SNPs_updated_alleles2.bim > test.txt
echo "Number of variants with zero alleles after allele update and removal of those not in manifest (should be zero): " >> report_1_pre_QC.txt
wc -l test.txt >> report_1_pre_QC.txt





                                                        #################################################
                                                        ####  Update chr, pos if applicable    ##########
                                                        #################################################

##################### update chr,pos if available ( if there are not unknown pos alleles, this will still work because file to exclude or update those info will be empty)
#get th list of variants with unkown pos
echo "           UPDATE UNKOWN VARIANT MAP IF APLICABLE               "

awk '{if($1 == "0" ) print $2}' 02_${name}_SNPs_updated_alleles2.bim > 02_list_variants_unknown_pos.txt
echo "Number of variants with unkown position: " >> report_1_pre_QC.txt
wc -l 02_list_variants_unknown_pos.txt >> report_1_pre_QC.txt

## are these in manifest file? 
awk -F "," 'FNR==NR{a[$1];next} ($2 in a)' 02_list_variants_unknown_pos.txt $manifest > 02_unknown_pos_variants_manifest_info.txt

## add pos for those missing
awk -F "," '{ print $2,$10,$11}' 02_unknown_pos_variants_manifest_info.txt > 02_unknown_pos_variants_positions_from_manifest.txt

# select those without zero pos
awk '{if($2 != "0" ) print $0}' 02_unknown_pos_variants_positions_from_manifest.txt > 02_unknown_pos_variants_manifest_nozero.txt
echo "Number of unknown posiotion variants with known pos in manifest: " >> report_1_pre_QC.txt
wc -l 02_unknown_pos_variants_manifest_nozero.txt >> report_1_pre_QC.txt

#update chr and bp pos in plink
plink1.9 --bfile 02_${name}_SNPs_updated_alleles2 --update-chr 02_unknown_pos_variants_manifest_nozero.txt 2 1 --update-map 02_unknown_pos_variants_manifest_nozero.txt 3 1 --allow-no-sex --keep-allele-order --make-bed --out 03_${name}_SNPs_updated_alleles_chr_map

#extract the list of variants with no chr, bp info after updating the chr/map
awk '{if($1 == "0" ) print $2}' 03_${name}_SNPs_updated_alleles_chr_map.bim > 03_list_variants_unknown_pos_to_remove.txt

#exclude variants with no pos with plink
plink1.9 --bfile 03_${name}_SNPs_updated_alleles_chr_map --exclude 03_list_variants_unknown_pos_to_remove.txt --allow-no-sex --keep-allele-order --make-bed --out 04_${name}_SNPs_updated_alleles_no_unkown_pos
echo "Plink:Number of variants after removing those with unknown pos in manifest: " >> report_1_pre_QC.txt
wc -l 04_${name}_SNPs_updated_alleles_no_unkown_pos.bim >> report_1_pre_QC.txt



                                                        #################################################
                                                        ############ Update Indels      #############
                                                        #################################################
echo '                THIS IS THE START OF UPDATE INDEL PROCESS         '

#extract indels
awk '{if($5 =="I" || $5 == "D" || $6 == "I" || $6 == "D") print $0}' 04_${name}_SNPs_updated_alleles_no_unkown_pos.bim > 04_list_indels.txt
echo "Number of variants with InDels coded as D/I: " >> report_1_pre_QC.txt 
wc -l 04_list_indels.txt >> report_1_pre_QC.txt

### use R code
rm R_report_InDel_update.txt
echo " InDels are updated using indel_update.R.... The report can be found in R_report_InDel_update.txt">> report_1_pre_QC.txt 
path=$(pwd)
Rscript indel_update.R $path 04_list_indels.txt $manifest


#update indel
plink1.9 --bfile 04_${name}_SNPs_updated_alleles_no_unkown_pos --update-alleles allele_change_Indels.txt  --allow-no-sex --keep-allele-order --make-bed --out 04_${name}_SNPs_updated_alleles_no_unkown_pos_indel_update


#remove any indels not in manifest file (found by R script)
plink1.9 --bfile 04_${name}_SNPs_updated_alleles_no_unkown_pos_indel_update --allow-no-sex --keep-allele-order --exclude InDels_not_in_manifest_to_remove.txt --make-bed --out 04_${name}_SNPs_updated_alleles_no_unkown_pos_indel_update2


#double check if everything is updated
awk '{if($5 =="I" || $5 == "D" || $6 == "I" || $6 == "D") print $0}' 04_${name}_SNPs_updated_alleles_no_unkown_pos_indel_update.bim > test.txt
echo "Re-check the number of variants with D/I coding after update and removing those not in manifest file: " >> report_1_pre_QC.txt 
wc -l test.txt >> report_1_pre_QC.txt 

                                                       ##################################################  
                                                        ############ remove duplicate vars  #############
                                                        #################################################

echo '                THIS IS THE START OF DUPLICATE VARIANT REMOVAL         '


############### ficates. add a col whith new names for all markers based on chr and pos
awk -F "\t" 'OFS="\t" { print $1,$2,$3,$4,$5,$6,$1"_"$4}' FS='\t' 04_${name}_SNPs_updated_alleles_no_unkown_pos_indel_update2.bim > 04_dummy_names_chr_pos.bim


#Check if there are any duplicate SNPs 
cut -f 7 04_dummy_names_chr_pos.bim | sort | uniq -D   > 04_dummy_names_chr_pos_dups.txt
echo "Number of variants with duplicate dummy IDs (chr_pos): " >> report_1_pre_QC.txt >> report_1_pre_QC.txt
wc -l 04_dummy_names_chr_pos_dups.txt >> report_1_pre_QC.txt


#extract these data from the bim file
awk 'FNR==NR{a[$1];next} ($7 in a)' 04_dummy_names_chr_pos_dups.txt 04_dummy_names_chr_pos.bim > 04_dummy_names_chr_pos_dups_bim.txt
cut -f 2 04_dummy_names_chr_pos_dups_bim.txt > list_duplicates_to_extract.txt

# extract thes duplicates, to calculate missingness (call rate)
plink1.9 --bfile  04_${name}_SNPs_updated_alleles_no_unkown_pos_indel_update2 --allow-no-sex --keep-allele-order --extract list_duplicates_to_extract.txt --make-bed --out 04_duplicates_SNPs

#calculate the missing rate (call rate)
plink1.9 --bfile 04_duplicates_SNPs --keep-allele-order --allow-no-sex --missing --out 04_duplicates_SNPs

#convert to vcf
plink1.9 --bfile 04_duplicates_SNPs --no-pheno --recode vcf-iid --keep-allele-order --out 04_duplicates_SNPs

#remove the # before CHR to read as headr into R
awk '/^#CHROM/{sub(/^#/,"")} 1' 04_duplicates_SNPs.vcf > 04_duplicates_SNPs_edited.vcf


#add the dummy IDs to the dup only bim file to read in R
awk -F "\t" 'OFS="\t" { print $1,$2,$3,$4,$5,$6,$1"_"$4}' FS='\t' 04_duplicates_SNPs.bim > 04_duplicates_SNPs_dummyID_for_R.txt

rm R_report_dup_var_removal.txt
echo "Duplicate variants are being removed using dup_var_removal.R. This will take a while, please be patient! the report is in file R_report_dup_var_removal.txt" >> report_1_pre_QC.txt
path=$(pwd)
Rscript dup_var_removal.R $path 04_duplicates_SNPs_dummyID_for_R.txt 04_duplicates_SNPs.lmiss 04_duplicates_SNPs_edited.vcf $dup_cutoff

#exclude the duplicates
plink1.9 --bfile 04_${name}_SNPs_updated_alleles_no_unkown_pos_indel_update2 --exclude  dups_to_be_excluded.txt --allow-no-sex --keep-allele-order --make-bed --out 05_${name}_SNPs_updated_alleles_no_unkown_pos_indel_update_no_dups
echo "Num of variants after duplicates (same pos, same alleles (or RC alleles)) removed: " >> report_1_pre_QC.txt
wc -l 05_${name}_SNPs_updated_alleles_no_unkown_pos_indel_update_no_dups.bim >> report_1_pre_QC.txt


# double check again with plink
plink1.9 --bfile 05_${name}_SNPs_updated_alleles_no_unkown_pos_indel_update_no_dups --allow-no-sex --keep-allele-order --list-duplicate-vars --out dups_check_after_removal

echo "number of duplicate pairs found by Plink after duplicate removal: " >> report_1_pre_QC.txt
line_num=$( cat dups_check_after_removal.dupvar| wc -l)
expr $line_num - 1 >> report_1_pre_QC.txt 
 my manifest file did not have all variants with one allele zero could not be resolved, those will be removed for the rest of the analysis
# find variants with one allele =0 

awk '{if($5 == "0" && $6 != "0") print $2}' 05_${name}_SNPs_updated_alleles_no_unkown_pos_indel_update_no_dups.bim > 05_list_SNPs_one_allele0.txt
echo "If certain variants with one missing allele do not exist in manifest file, find them and remove them before next step. The number of these variants is: " >> report_1_pre_QC.txt
wc -l 05_list_SNPs_one_allele0.txt >> report_1_pre_QC.txt

## exclude the variants with one 0 allele that could not be fixed with manifest file
plink1.9 --bfile 05_${name}_SNPs_updated_alleles_no_unkown_pos_indel_update_no_dups --exclude  05_list_SNPs_one_allele0.txt --allow-no-sex --keep-allele-order --make-bed --out 05_${name}_SNPs_updated_alleles_no_unkown_pos_indel_update_no_dups_cleaned
echo "the final number of variants for the pre_QC step is: " >> report_1_pre_QC.txt
wc -l 05_${name}_SNPs_updated_alleles_no_unkown_pos_indel_update_no_dups_cleaned.bim >> report_1_pre_QC.txt

echo "the final number of individuals for the pre_QC step is: " >> report_1_pre_QC.txt
wc -l 05_${name}_SNPs_updated_alleles_no_unkown_pos_indel_update_no_dups_cleaned.fam >> report_1_pre_QC.txt


echo "the distribution of variants on chrs in the final pre_QC file is: " >> report_1_pre_QC.txt
awk '{print $1}' 05_${name}_SNPs_updated_alleles_no_unkown_pos_indel_update_no_dups_cleaned.bim | sort | uniq -c >> report_1_pre_QC.txt


rm test.txt
echo "Pre_QC is completed! "


