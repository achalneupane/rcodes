#!/bin/bash
### generate_assoc.sh - Modified by Achal on 2021-07-29 from  2017-0-20 - vifehe
### script to put together all pvals from 
### USAGE SET="geneset-MAF-0.005" CHR="" ./generate_assoc.sh

SET="${SET}"  ## set being geneset-MAF1P, geneset-CADD20
chr="${CHR}"
TEST="${TEST}"  ## test being: SKAT // SKAT-COHORT_covar // SKAT-COHORT_APOE_covar

	## ADD header to file (so all will have same amount of lines)
	cut -f1 chr${chr}/${SET}_chr${chr}_${TEST}.pval  > ${SET}_chr${chr}_${TEST}_GENES.tmp;
	#echo -e "${SET}_GENE ${SET}_SNPs" > a; cat a ${SET}_chr${chr}_GENE-SNPs > b ;  mv b ${SET}_chr${chr}_GENE-SNPs;
	# Create a temporary SNP file from which to extract CHR POS REF ALT for assoc file
	cut -f 2 chr${chr}/${SET}_chr${chr}_${TEST}.pval > ${SET}_chr${chr}_${TEST}_SNPs.tmp;
	echo "CHR" > ${SET}_chr${chr}_${TEST}_CHR.tmp ; tail -n+2 ${SET}_chr${chr}_${TEST}_SNPs.tmp | cut -d ":"  -f 1  >> ${SET}_chr${chr}_${TEST}_CHR.tmp;
	# sed -i "s/^X/23/g" ${SET}_chr${chr}_${TEST}_CHR.tmp;
	# sed -i "s/^Y/24/g" ${SET}_chr${chr}_${TEST}_CHR.tmp;
	echo "POS" > ${SET}_chr${chr}_${TEST}_POS.tmp ; tail -n+2 ${SET}_chr${chr}_${TEST}_SNPs.tmp | cut -d ":"  -f 2  >> ${SET}_chr${chr}_${TEST}_POS.tmp;
	echo "REF" > ${SET}_chr${chr}_${TEST}_REF.tmp ; tail -n+2 ${SET}_chr${chr}_${TEST}_SNPs.tmp | cut -d ":"  -f 3  >> ${SET}_chr${chr}_${TEST}_REF.tmp;
	echo "ALT" > ${SET}_chr${chr}_${TEST}_ALT.tmp ; tail -n+2 ${SET}_chr${chr}_${TEST}_SNPs.tmp | cut -d ":"  -f 4  >> ${SET}_chr${chr}_${TEST}_ALT.tmp;
	echo "N0" > ${SET}_chr${chr}_${TEST}_N0.tmp ; tail -n+2 chr${chr}/${SET}_chr${chr}_${TEST}.pval | cut -f 3  >> ${SET}_chr${chr}_${TEST}_N0.tmp;
	echo "N1" > ${SET}_chr${chr}_${TEST}_N1.tmp ; tail -n+2 chr${chr}/${SET}_chr${chr}_${TEST}.pval | cut -f 4  >> ${SET}_chr${chr}_${TEST}_N1.tmp;
	# EXTRACT PVALS COLUMNS
	cut -f 5 chr${chr}/${SET}_chr${chr}_${TEST}.pval > ${SET}_chr${chr}_${TEST}_SKAT.tmp;
	cut -f 6 chr${chr}/${SET}_chr${chr}_${TEST}.pval > ${SET}_chr${chr}_${TEST}_SKATO.tmp;
	cut -f 7 chr${chr}/${SET}_chr${chr}_${TEST}.pval > ${SET}_chr${chr}_${TEST}_SKAT125.tmp;
	cut -f 8 chr${chr}/${SET}_chr${chr}_${TEST}.pval > ${SET}_chr${chr}_${TEST}_SKAT_C.tmp;
	
	## CREATE FILES WITH 0 for values F_A F_U OR CHISQ
	# Count how many genes per set:
	END=`cat chr${chr}/${SET}_chr${chr}_${TEST}.pval | wc -l`;
	#for i in $(seq 1 $END); do echo $i; done
	# CREATE COLUMns AS long as the number of genes ( IT AHS TO BE 2 units less thant the TOTAL number of genes present)
	echo "OR" >  ${SET}_chr${chr}_${TEST}_OR.tmp ; for x in $(seq 2 $END); do echo -e "0"; done >>  ${SET}_chr${chr}_${TEST}_OR.tmp;
	echo "AF" >  ${SET}_chr${chr}_${TEST}_AF.tmp ; for x in $(seq 2 $END); do echo -e "0"; done >>  ${SET}_chr${chr}_${TEST}_AF.tmp;
	echo "CHISQ" >  ${SET}_chr${chr}_${TEST}_CHISQ.tmp ; for x in $(seq 2 $END); do echo -e "0"; done >>  ${SET}_chr${chr}_${TEST}_CHISQ.tmp;
	echo "FU" >  ${SET}_chr${chr}_${TEST}_FU.tmp ; for x in $(seq 2 $END); do echo -e "0"; done >>  ${SET}_chr${chr}_${TEST}_FU.tmp;
	echo "FA" >  ${SET}_chr${chr}_${TEST}_FA.tmp ; for x in $(seq 2 $END); do echo -e "0"; done >>  ${SET}_chr${chr}_${TEST}_FA.tmp;
	## Put intermediate files together in an assoc file
	paste  ${SET}_chr${chr}_${TEST}_CHR.tmp  ${SET}_chr${chr}_${TEST}_GENES.tmp ${SET}_chr${chr}_${TEST}_POS.tmp  ${SET}_chr${chr}_${TEST}_REF.tmp  ${SET}_chr${chr}_${TEST}_N0.tmp  ${SET}_chr${chr}_${TEST}_N1.tmp  ${SET}_chr${chr}_${TEST}_ALT.tmp  ${SET}_chr${chr}_${TEST}_CHISQ.tmp  ${SET}_chr${chr}_${TEST}_SKAT.tmp  ${SET}_chr${chr}_${TEST}_OR.tmp >  ${SET}_chr${chr}_${TEST}_SKAT.assoc;
	paste  ${SET}_chr${chr}_${TEST}_CHR.tmp  ${SET}_chr${chr}_${TEST}_GENES.tmp ${SET}_chr${chr}_${TEST}_POS.tmp  ${SET}_chr${chr}_${TEST}_REF.tmp  ${SET}_chr${chr}_${TEST}_N0.tmp  ${SET}_chr${chr}_${TEST}_N1.tmp  ${SET}_chr${chr}_${TEST}_ALT.tmp  ${SET}_chr${chr}_${TEST}_CHISQ.tmp  ${SET}_chr${chr}_${TEST}_SKATO.tmp  ${SET}_chr${chr}_${TEST}_OR.tmp >  ${SET}_chr${chr}_${TEST}_SKATO.assoc;
	paste  ${SET}_chr${chr}_${TEST}_CHR.tmp  ${SET}_chr${chr}_${TEST}_GENES.tmp ${SET}_chr${chr}_${TEST}_POS.tmp  ${SET}_chr${chr}_${TEST}_REF.tmp  ${SET}_chr${chr}_${TEST}_N0.tmp  ${SET}_chr${chr}_${TEST}_N1.tmp  ${SET}_chr${chr}_${TEST}_ALT.tmp  ${SET}_chr${chr}_${TEST}_CHISQ.tmp  ${SET}_chr${chr}_${TEST}_SKAT125.tmp  ${SET}_chr${chr}_${TEST}_OR.tmp >  ${SET}_chr${chr}_${TEST}_SKAT125.assoc;
	paste  ${SET}_chr${chr}_${TEST}_CHR.tmp  ${SET}_chr${chr}_${TEST}_GENES.tmp ${SET}_chr${chr}_${TEST}_POS.tmp  ${SET}_chr${chr}_${TEST}_REF.tmp  ${SET}_chr${chr}_${TEST}_N0.tmp  ${SET}_chr${chr}_${TEST}_N1.tmp  ${SET}_chr${chr}_${TEST}_ALT.tmp  ${SET}_chr${chr}_${TEST}_CHISQ.tmp  ${SET}_chr${chr}_${TEST}_SKAT_C.tmp  ${SET}_chr${chr}_${TEST}_OR.tmp >  ${SET}_chr${chr}_${TEST}_SKAT_C.assoc;


	## CHECK THEY AL KEEP SAME AMOUNT OF LINES
	wc -l  ${SET}_chr${chr}_${TEST}_*;
	rm *.tmp


