#!/bin/bash

###To run script call ./Liftover.bash Clean_plink_file_path ouput_prefix liftover_builds

clean_plink=$1
name=$2 
liftover=$3
rm report_liftover.txt

#Converting chr 23,24,25,26 to X, Y, XY and M
# --merge-x changes chromosome codes of all XY variants back to X (and 'no-fail' has the same effect)
#--output-chr lets you specify a different coding scheme by providing the desired human mitochondrial code; supported options are '26' (default), 'M', 'MT', '0M', 'chr26', 'chrM', and 'chrMT'. (PLINK 1.9 correctly interprets all of these encodings in input files.)
echo "Converting chr 23,24,25,26 to X, Y, XY and M ..." > report_liftover.txt
plink1.9 --bfile ${clean_plink} --output-chr MT --merge-x 'no-fail' --allow-no-sex --keep-allele-order --make-bed --out ${name}_chrM


#Creating bed file for liftover - had to do the start position -1 because the bed is 0-based according with UCSC
awk '{OFS="\t"}; {print "chr"$1,$4-1,$4,$2}'  ${name}_chrM.bim > ${name}_tolift.bed

$liftover ${name}_tolift.bed
echo "Lifted over the positions." >> report_liftover.txt

#this will show all chr variation in this file: awk '{print $1}' QC_cleaned_DYS_WT_tolift.bed-hg38 | sort | uniq -c

echo "Flagging problematic variants for removal: " >> report_liftover.txt

rm to_remove.txt
grep -v  '#Deleted' ${name}_tolift.bed-unmapped | awk '{print $4}' > to_remove.txt
echo "Number of unmapped variants (Deleted) :" >> report_liftover.txt
wc -l to_remove.txt >> report_liftover.txt

grep '_alt'  ${name}_tolift.bed-hg38 | awk '{print $4}'  >> to_remove.txt
echo "Number of chr_alt variants :" >> report_liftover.txt
grep -c '_alt'  ${name}_tolift.bed-hg38 >> report_liftover.txt

grep '_random'  ${name}_tolift.bed-hg38 | awk '{print $4}' >> to_remove.txt
echo "Number of chrN_random variants :" >> report_liftover.txt
grep -c '_random'  ${name}_tolift.bed-hg38 >> report_liftover.txt

grep 'chrUn_'  ${name}_tolift.bed-hg38 | awk '{print $4}' >> to_remove.txt
echo "Number of chrUn variants :" >> report_liftover.txt
grep -c 'chrUn_'  ${name}_tolift.bed-hg38 >> report_liftover.txt

#Exclude variants in alternative or random chromosomes and that failed liftover

plink1.9 --bfile ${clean_plink} --exclude to_remove.txt --make-bed --out Temp
echo "Excluding problematic variants using Plink" >> report_liftover.txt

#Creating files to update plink file to hg38

awk '{print $4, $3}'  ${name}_tolift.bed-hg38 > update_pos.txt

awk '{print $4, $1}'  ${name}_tolift.bed-hg38 > update_chr.txt

#Updating chr and position to hg38

plink1.9 --bfile Temp --update-chr update_chr.txt --update-map update_pos.txt --allow-no-sex --keep-allele-order --make-bed --out ${name}_hg38
echo "Updated chr and map in plink." >> report_liftover.txt

echo " Number of variants in the liftover file:" >> report_liftover.txt
wc -l ${name}_hg38.bim >> report_liftover.txt

echo "done"

#mkdir Plink_log_files
#mv *.log ./Plink_log_files
rm *_chrM*
rm *Temp*



