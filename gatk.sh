srun -N 1 --cpus-per-task=10 --time=5:00:00 --partition=compute --pty bash


# download data
wget https://de.cyverse.org/dl/d/3CE425D7-ECDE-46B8-AB7F-FAF07048AD42/samples.tar.gz
  tar xvzf samples.tar.gz
  rm samples.tar.gz


mkdir GATK_tutorial
cd GATK_tutorial


mkdir dbSNP
mkdir SAM
mkdir BAM
mkdir sortedBAM
mkdir SNPs
mkdir ${SNPs}MERGED
mkdir ${SNPs}persample

GATK="${PWD}/"
dbSNP="${GATK}/dbSNP/"
SAM="${GATK}/SAM/"
BAM="${GATK}/BAM/"
sortedBAM="${GATK}/sortedBAM/"
SNPs="${GATK}/SNPs/"
MERGED="${SNPs}MERGED/"
outdir="${SNPs}persample/"


#load modules
module load samtools
module load bwa
module load R

Rscript -e "install.packages('ggplot2', contriburl=contrib.url('http://cran.r-project.org/'))"

Rscript -e "install.packages('gplots', contriburl=contrib.url('http://cran.r-project.org/'))"

Rscript -e "install.packages('reshape', contriburl=contrib.url('http://cran.r-project.org/'))"

Rscript -e "install.packages('gsalib', contriburl=contrib.url('http://cran.r-project.org/'))"

Rscript -e "install.packages('Biobase', contriburl=contrib.url('http://bioconductor.org/packages/release/bioc/'))"



# download reference genome
wget https://de.cyverse.org/dl/d/A9330898-FC54-42A5-B205-B1B2DC0E91AE/dog_chr5.fa.gz
gunzip dog_chr5.fa.gz
bwa index -a bwtsw dog_chr5.fa


#########RUN THIS: ADD read Group!!!!#############

for R1 in *_R1_001.pe.fq.gz;do
    SM=$(echo $R1 | cut -d"_" -f1)                                          ##sample ID
    LB=$(echo $R1 | cut -d"_" -f1,2)                                        ##library ID
    PL="Illumina"                                                           ##platform (e.g. illumina, solid)
    RGID=$(zcat $R1 | head -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)       ##read group identifier 
    PU=$RGID.$LB                                                            ##Platform Unit
    echo -e "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU"

    R2=$(echo $R1 | sed 's/_R1_/_R2_/')
    echo $R1 $R2
    bwa mem -t 4 -M -R "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU" dog_chr5.fa $R1 $R2 > ${SAM}${R1%_R1_001.pe.fq.gz}.sam
  done



# Generate BAM file
for samfile in ${SAM}*.sam;do
  sample=$(basename "${samfile%.*}")
  echo "Doing: " $sample
  samtools view -bS -o ${BAM}${sample}.bam $samfile
  echo "Created: " ${BAM}${sample}.bam
  samtools sort ${BAM}${sample}.bam -o ${sortedBAM}${sample}.sorted.bam
  echo "Created: " ${sortedBAM}${sample}.sorted.bam
done







wget https://github.com/broadinstitute/picard/releases/download/2.9.4/picard.jar
chmod u+x *




wget https://de.cyverse.org/dl/d/6177B1E0-718A-4F95-A83B-C3B88E23C093/GenomeAnalysisTK-3.7-0.tar.bz2
tar xjf GenomeAnalysisTK-3.7-0.tar.bz2



java -Xmx10g -jar ${GATK}picard.jar CreateSequenceDictionary R=dog_chr5.fa O=dog_chr5.dict
 samtools faidx dog_chr5.fa


# Merge BAM replicates

java  -jar ${GATK}picard.jar MergeSamFiles I="${sortedBAM}BD143_TGACCA_L005.sorted.bam" I="${sortedBAM}BD143_TGACCA_L006.sorted.bam" OUTPUT="${sortedBAM}BD143_TGACCA_merged.sorted.bam"




#Mark duplicates


for sample in ${sortedBAM}*.sorted.bam;do
  #name=${sample%.sorted.bam}
  #name=$(basename "${sample%.*}")
  name=$(basename "${sample%.sorted.bam}")
  echo "Doing: " $name
  java  -Xmx10g -jar ${GATK}picard.jar MarkDuplicates INPUT=$sample OUTPUT=${sortedBAM}${name}.dedup.bam METRICS_FILE=$name.metrics.txt;
done

# cd sortedBAM
samtools view -H ${sortedBAM}BD143_TGACCA_L005.sorted.bam
samtools view -H ${sortedBAM}BD143_TGACCA_L006.sorted.bam
samtools view -H ${sortedBAM}BD143_TGACCA_merged.sorted.bam
samtools view -H ${sortedBAM}BD143_TGACCA_merged.dedup.bam

rm ${sortedBAM}BD143_TGACCA_L00*.sorted.bam



http://m.ensembl.org/Canis_familiaris/Info/Annotation#assembly

# wget 'ftp://ftp.ensembl.org/pub/release-89/variation/vcf/canis_familiaris/Canis_familiaris.vcf.gz'
wget 'http://ftp.ensembl.org/pub/release-99/variation/vcf/canis_familiaris/canis_familiaris.vcf.gz' -O ${dbSNP}canis_familiaris.vcf.gz

# mv Canis_familiaris.vcf.gz canis_familiaris.vcf.gz

gunzip -c ${dbSNP}canis_familiaris.vcf.gz > ${dbSNP}canis_familiaris.vcf
grep "^#" ${dbSNP}canis_familiaris.vcf > ${dbSNP}canis_fam_chr5.vcf
grep "^5" ${dbSNP}canis_familiaris.vcf | sed 's/^5/chr5/' >> ${dbSNP}canis_fam_chr5.vcf
# Run Recalibration
# BQSR stands for Base Quality Score Recalibration.


for sample in ${sortedBAM}*.dedup.bam;do
  #name=${sample%.dedup.bam}
  name=$(basename "${sample%.dedup.bam}")
  echo "Doing: " $name
  samtools index $sample
  java -Xmx10g -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R ${GATK}dog_chr5.fa -I $sample -knownSites ${dbSNP}canis_fam_chr5.vcf -o ${SNPs}${name}.1st.table
  java -Xmx10g -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R ${GATK}dog_chr5.fa -I $sample -knownSites ${dbSNP}canis_fam_chr5.vcf -BQSR ${SNPs}${name}.1st.table -o ${SNPs}${name}.2nd.table
  java -Xmx10g -jar GenomeAnalysisTK.jar -T PrintReads -R ${GATK}dog_chr5.fa -I $sample -BQSR ${SNPs}${name}.2nd.table -o ${SNPs}${name}.recal.bam
  java -Xmx10g -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R ${GATK}dog_chr5.fa -before ${SNPs}${name}.1st.table -after ${SNPs}${name}.2nd.table -plots ${SNPs}${name}.BQSR.pdf
done




for sample in ${SNPs}*.recal.bam;do
  name=$(basename "${sample%.recal.bam}")
  java -Xmx10g -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R ${GATK}dog_chr5.fa --dbsnp ${dbSNP}canis_fam_chr5.vcf -I $sample --emitRefConfidence GVCF -nct 3 -o ${outdir}${name}.g.vcf
done


#Combine call
java -Xmx10g -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R ${GATK}dog_chr5.fa --dbsnp ${dbSNP}canis_fam_chr5.vcf \
 --variant ${outdir}BD143_TGACCA_merged.g.vcf \
 --variant ${outdir}BD174_CAGATC_L005.g.vcf \
 --variant ${outdir}BD225_TAGCTT_L007.g.vcf \
 -o ${MERGED}raw_variants.vcf


 # split variants into SNPs and indels
mkdir ${MERGED}/SNPs
mkdir ${MERGED}/INDELs
java -Xmx10g -jar ${GATK}GenomeAnalysisTK.jar -T SelectVariants -R ${GATK}dog_chr5.fa -V ${MERGED}raw_variants.vcf -selectType SNP -o "${MERGED}/SNPs/"raw_SNP.vcf 
java -Xmx10g -jar ${GATK}GenomeAnalysisTK.jar -T SelectVariants -R ${GATK}dog_chr5.fa -V ${MERGED}raw_variants.vcf -selectType INDEL -o "${MERGED}/INDELs/"raw_INDEL.vcf 


# Distribution of variants
cd $MERGED
mkdir both
cp INDELs/* ./both/
cp SNPs/* ./both/
cd "$MERGED/both"
wget https://raw.githubusercontent.com/drtamermansour/angus/2017/densityCurves.R
for var in "SNP" "INDEL";do
 for ann in "QD" "MQRankSum" "FS" "SOR" "ReadPosRankSum";do
  annFile=$var.$ann; echo $annFile;
  awk -v k="$ann=" '!/#/{n=split($8,a,";"); for(i=1;i<=n;i++) if(a[i]~"^"k) {sub(k,$3" ",a[i]); print a[i]}}' raw_$var.vcf > $annFile
  grep -v "^\." $annFile > known.$annFile
  grep "^\." $annFile > novel.$annFile
  Rscript densityCurves.R "$annFile"
  rm $annFile known.$annFile novel.$annFile
done; done

# Apply filters

java -Xmx10g -jar GenomeAnalysisTK.jar -T VariantFiltration -R dog_chr5.fa -V raw_SNP.vcf \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
--filterName "snp_filter" \
-o filtered_SNP.vcf

java -Xmx10g -jar GenomeAnalysisTK.jar -T VariantFiltration -R dog_chr5.fa -V raw_INDEL.vcf \
--filterExpression "QD < 2.0 || FS > 200.0" \
--filterName "indel_filter" \
-o filtered_INDEL.vcf


#NEXT>>>>>> R Programming
GWAS
Family-Based/Gene-based Study
Ethnicity Mapping
Cancer gene discovery/Clonal evolution
Personalized medicine





