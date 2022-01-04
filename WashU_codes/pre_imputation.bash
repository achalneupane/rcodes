#!/bin/bash

clean_plink=$1
name=$2

m report_pre_imputation.txt

#Converting chr 23,24 to X, Y, removed M, removed indels
plink1.9 --bfile ${clean_plink} --output-chr M --chr 1-22, X,Y,XY --merge-x 'no-fail' --snps-only no-DI --make-bed --out ${name}_no_chrM
echo "Converted chr 23,24 to X, Y, removed M, removed indels..." > report_pre_imputation.txt

#you can check the chr (1-22, + x and Y)
echo "Number of variants on each chr after removing M and updating chr coding:" >> report_pre_imputation.txt
awk '{print $1}' ${name}_no_chrM.bim | sort | uniq -c >> report_pre_imputation.txt


#Create a symbolic link to the reference and index file
ln -s /data/GATK_pipeline/20190522_bundle/Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta

../vcfCooker --in-bfile ${name}_no_chrM --ref Homo_sapiens_assembly38.fasta --out ${name}_no_chrM_ref_annotated.vcf.gz --write-vcf --bgzf
echo "Corrected ref and alt and strand based on reference using vcfCooker." >> report_pre_imputation.txt

#BCFtools snp name change to chr:pos:ref:alt
tabix -p vcf ${name}_no_chrM_ref_annotated.vcf.gz
bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%ALT' ${name}_no_chrM_ref_annotated.vcf.gz -Oz  > ${name}_no_chrM_ref_annotated_renamed.vcf.gz
tabix -p vcf ${name}_no_chrM_ref_annotated_renamed.vcf.gz
echo "Renamed and Indexed VCF file." >> report_pre_imputation.txt


#split vcf into indv chr and index
parallel -j 6 bcftools view -r {} ${name}_no_chrM_ref_annotated_renamed.vcf.gz -Oz -o {}_hg38_ref_annotatedSNP.vcf.gz :::  $(cut -f1 chr_list.txt)
parallel -j 6 bcftools sort {} -Oz -o sorted_{} ::: chr*.vcf.gz 
parallel -j 6 tabix -p vcf {} ::: sorted_chr*.vcf.gz 
echo "VCF file is split into chrs and indexed" >> report_pre_imputation.txt

rm chr*.vcf.gz
rm Homo_sapiens_assembly38.fasta


