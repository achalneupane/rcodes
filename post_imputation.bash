#!/bin/bash

#create "imputation_sever_results" folder inside "05_post_imputation" folder and save all the downloaded files from TOPmed server here
#save the password in file "passwd.txt" in directory 05_post_imputation

name=$1
pre_imputation=$2
dup_cutoff=$3

cd imputation_sever_results

#check the md5sum check
md5sum *.zip | sort -k 2 | sed 's/  / /g'  > downloaded_zips_md5sum.txt
sort -k 2 results.md5 > results_sorted.md5
echo "The difference between md5sum reported by TOPmed and downloaded files" > ../report_post_imputation.txt
diff downloaded_zips_md5sum.txt results_sorted.md5 | grep \^\< > compare_md5sum.txt >> ../report_post_imputation.txt
cd ..

#Unzipping files in the main directory
7za e "./imputation_sever_results/*.zip" -p$(cat passwd.txt)


echo "my current folder is:" >> report_post_imputation.txt
pwd >> report_post_imputation.txt

#Converting vcfs to plink
echo "Converting vcfs to plink format, this may take a while... " >> report_post_imputation.txt
parallel -j 2 plink1.9 --vcf {} --double-id  --allow-no-sex --keep-allele-order --make-bed --out '{= s:\.[^.]+$::;s:\.[^.]+$::;s:\.[^.]+$::;=}' ::: chr*.dose.vcf.gz

#create a list of filenames to be merged
ls *.bed > file_list.txt
sed -i -e 's/.bed//g' file_list.txt

#Merge the imputed chromosomes together
echo "Merge the imputed chromosomes together , this may take a while... " >> report_post_imputation.txt
plink1.9 --merge-list file_list.txt --keep-allele-order --allow-no-sex --make-bed --out ${name}_hg38imputed

#Running Hardy-Weinberg threshold 1e-6
plink1.9 --bfile ${name}_hg38imputed --hwe 1e-6 midp --keep-allele-order --allow-no-sex --make-bed --out ${name}_hg38imputed_hwe
echo "Filtered IMPUTED file for HWE threshold 1e-6. The number of variants after HWE filter is : " >> report_post_imputation.txt
wc -l ${name}_hg38imputed_hwe.bim >> report_post_imputation.txt

#Converting SNP renamed final pre_imputation file with all chr from vcf to plink 
echo "Converting SNP renamed final pre_imputation file with all chr from vcf to plink... " >> report_post_imputation.txt
plink1.9 --vcf ${pre_imputation} --double-id --keep-allele-order --allow-no-sex --make-bed --out ${name}_Pre_imputation_final_No_chrM_ref_annotated_renamed

#Remove invalid alleles, maybe the ones that mismatch from the SNP renamed file
grep "Invalid" ./imputation_sever_results/snps-excluded.txt | awk '{print $1}' > SNPs_remove.txt
grep "mismatch" ./imputation_sever_results/snps-excluded.txt | awk '{print $1}' >> SNPs_remove.txt

echo "Number of invalid alleles and mismatches to be removed from pre_imputation file:" >> report_post_imputation.txt
wc -l SNPs_remove.txt  >> report_post_imputation.txt

plink1.9 --bfile ${name}_Pre_imputation_final_No_chrM_ref_annotated_renamed --exclude SNPs_remove.txt --allow-no-sex --keep-allele-order --make-bed --out ${name}_Pre_imputation_final_No_chrM_ref_annotated_renamed_clean
echo "Number of variants in the pre_imputation file after removing the invalid and mismatch alleles:" >> report_post_imputation.txt
wc -l ${name}_Pre_imputation_final_No_chrM_ref_annotated_renamed_clean.bim >> report_post_imputation.txt


#Merge imputed file to clean SNP renamed file (merge-mode 2 :Only overwrite calls which are missing in the original file.)
echo "Merging imputed file with cleaned SNP renamed file (merge-mode 2 :Only overwrite calls which are missing in the original file.). This may take a while..." >> report_post_imputation.txt
plink1.9 --bfile ${name}_Pre_imputation_final_No_chrM_ref_annotated_renamed_clean --bmerge ${name}_hg38imputed_hwe --allow-no-sex --keep-allele-order --merge-mode 2 --make-bed --out ${name}_hg38_imputed_merged_temp
echo "Number of variants after the merge is:" >> report_post_imputation.txt
wc -l ${name}_hg38_imputed_merged_temp.bim >> report_post_imputation.txt


################ Remove duplicate variants that are generated during the merge. Mostly are the same pos with inverted alleles.

######### approach 1: use dummyIDs genome wide:
### this gives a very large set but it is more comprehensive. I keep it here in case it will be needed.
#awk -F "\t" 'OFS="\t" { print $1,$2,$3,$4,$5,$6,$1"_"$4}' FS='\t' ${name}_hg38_imputed_merged_temp.bim > dummy_names_chr_pos.bim

#Check if there are any duplicate SNPs 
#cut -f 7 dummy_names_chr_pos.bim | sort | uniq -D   > dummy_names_chr_pos_dups.txt

#extract these data from the bim file
#awk 'FNR==NR{a[$1];next} ($7 in a)' dummy_names_chr_pos_dups.txt dummy_names_chr_pos.bim > dummy_names_chr_pos_dups_bim.txt
#cut -f 2 dummy_names_chr_pos_dups_bim.txt > list_duplicates_to_extract.txt


# extract these duplicates, to calculate missingness (call rate)
#plink1.9 --bfile ${name}_hg38_imputed_merged_temp --allow-no-sex --keep-allele-order --extract list_duplicates_to_extract.txt --make-bed --out duplicate_SNPs



############approach 2: use dups identified by plink (it does detect vars with inversted alleles as dups)
plink1.9 --bfile ${name}_hg38_imputed_merged_temp --allow-no-sex --keep-allele-order --list-duplicate-vars

cut -f 4 plink.dupvar > dups_from_plink.txt
awk 'gsub(/ /,"\n")' dups_from_plink.txt > dups_from_plink_1col.txt
echo "Number of duplicate variables identified by Plink --list-duplicate-vars: "  >> report_post_imputation.txt
wc -l dups_from_plink_1col.txt  >> report_post_imputation.txt

# extract these duplicates, to calculate missingness (call rate)
plink1.9 --bfile ${name}_hg38_imputed_merged_temp --allow-no-sex --keep-allele-order --extract dups_from_plink_1col.txt --make-bed --out duplicate_SNPs

############## This part is the same for both approaches
#calculate the missing rate (call rate)
plink1.9 --bfile duplicate_SNPs --keep-allele-order --allow-no-sex --missing --out duplicate_SNPs

#convert to vcf
plink1.9 --bfile duplicate_SNPs --no-pheno --recode vcf-iid --keep-allele-order --out duplicate_SNPs

#Remove vcf header
sed '/^#/d' duplicate_SNPs.vcf > duplicate_SNPs_noheader.vcf

#prep bim for R
awk -F "\t" 'OFS="\t" { print $1,$2,$3,$4,$5,$6,$1"_"$4}' FS='\t' duplicate_SNPs.bim > duplicate_SNPs_dummyID_for_R.txt

## identify duplicates in R and python and use the list to exclude in plink ( dups_to_be_excluded.txt)
rm R_report_dup_var_removal.txt
echo "Duplicate variants are being removed using dup_var_removal.R. This will take a while, please be patient! the report is in file R_report_dup_var_removal.txt" >> report_post_imputation.txt

path=$(pwd)
Rscript dup_var_removal_python_call_cor.R $path duplicate_SNPs_dummyID_for_R.txt duplicate_SNPs.lmiss $dup_cutoff

#Removing problematic SNPs and SNPs with call rates under 98%
plink1.9 --bfile ${name}_hg38_imputed_merged_temp --allow-no-sex --keep-allele-order --exclude dups_to_be_excluded.txt --geno 0.02 --make-bed --out ${name}_hg38_imputed_final

# double check again with plink
plink1.9 --bfile ${name}_hg38_imputed_final --allow-no-sex --keep-allele-order --list-duplicate-vars --out dups_check_after_removal
echo "number of duplicate pairs found by Plink after duplicate removal: " >> report_post_imputation.txt
line_num=$( cat dups_check_after_removal.dupvar| wc -l)
expr $line_num - 1 >> report_post_imputation.txt 

##report final numbers
echo "Final number of variants after removing duplicates:" >> report_post_imputation.txt
wc -l ${name}_hg38_imputed_final.bim >> report_post_imputation.txt

echo "Final number of indviduals:" >> report_post_imputation.txt
wc -l ${name}_hg38_imputed_final.fam >> report_post_imputation.txt







