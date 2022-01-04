## https://meyer-lab-cshl.github.io/plinkQC/articles/AncestryCheck.html
# https://meyer-lab-cshl.github.io/plinkQC/articles/plinkQC.html
# cd /40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/01-Aquilla-preQC/03-PLINK-QC-files3/plinkQC
## sudo docker run -w $PWD -ti --rm achalneupane/plinkqc
# sudo  docker run -ti --rm -v $PWD:/usr/local/src/  achalneupane/plinkqc bash




# install.packages("devtools")
# vim.tiny myscript.R
library(devtools)
# install_github("meyer-lab-cshl/plinkQC")
# install.packages("plinkQC")
library(plinkQC)



# Introduction
# Genotyping arrays enable the direct measurement of an individuals genotype at
# thousands of markers. Subsequent analyses such as genome-wide association
# studies rely on the high quality of these marker genotypes.

# Anderson and colleagues introduced a protocol for data quality control in
# genetic association studies heavily based on the summary statistics and
# relatedness estimation functions in the PLINK software [1]. PLINK is a
# comprehensive, open-source command-line tool for genome-wide association
# studies (GWAS) and population genetics research [2]. It’s main functionalities
# include data management, computing individual- and marker- level summary
# statistics, identity-by-state estimation and association analysis.

# Integration with R is achieved through its R plugin or PLINK/SEQ R Package
# [4]. While the plugin is limited to operations yielding simple genetic marker
# vectors as output, the PLINK/SEQ R Package is limited in the functionalities
# it can access.

# plinkQC facilitates genotype quality control for genetic association studies
# as described by [1]. It wraps around plink basic statistics (e.g. missing
# genotyping rates per individual, allele frequencies per genetic marker) and
# relationship functions and generates a per-individual and per-marker quality
# control report. Individuals and markers that fail the quality control can
# subsequently be removed with plinkQC to generate a new, clean dataset. Removal
# of individuals based on relationship status is optimised to retain as many
# individuals as possible in the study.

# plinkQC depends on the PLINK (version 1.9), which has to be manually installed
# prior to the usage of plinkQC. It assumes the genotype have already been
# determined from the original probe intensity data of the genotype array and is
# available in plink format.

# Note: PLINK 2.0 is still in alpha status with many re-implementations and
# updates, including output file changes. While these changes are ongoing,
# plinkQC will rely on users using PLINK 1.9. I will monitor PLINK 2.0 changes
# and check compatibility with plinkQC. If users require specific PLINK 2.0
# output, I recommend using the Step-by-Step approach described below and
# manually saving to the output files expected from PLINK 1.9.

# The protocol is implemented in three main functions, the per-individual
# quality control (perIndividualQC), the per-marker quality control
# (perMarkerQC) and the generation of the new, quality control dataset
# (cleanData):
  
  # Per-individual quality control
# The per-individual quality control with perIndividualQC wraps around these
# functions: (i) check_sex: for the identification of individuals with
# discordant sex information, (ii) check_heterozygosity_and_missingness: for the
# identification of individuals with outlying missing genotype and/or
# heterozygosity rates, (iii) check_relatedness: for the identification of
# related individuals, (iv) check_ancestry: identification of individuals of
# divergent ancestry.

# Per-marker quality control
# The per-marker quality control with perMarkerQC wraps around these functions:
# (i) check_snp_missingnes: for the identifying markers with excessive missing
# genotype rates, (ii) check_hwe: for the identifying markers showing a
# significant deviation from Hardy-Weinberg equilibrium (HWE), (iii) check_maf:
# for the removal of markers with low minor allele frequency (MAF).

# Clean data
# cleanData takes the results of perMarkerQC and perIndividualQC and creates a
# new dataset with all individuals and markers that passed the quality control
# checks.

# Workflow
# In the following, genotype quality control with plinkQC is applied on a small
# example dataset with 200 individuals and 10,000 markers (provided with this
# package). The quality control is demonstrated in three easy steps,
# per-individual and per-marker quality control followed by the generation of
# the new dataset. In addition, the functionality of each of the functions
# underlying perMarkerQC and perIndividualQC is demonstrated at the end of this
# vignette.

package.dir <- find.package('plinkQC')
indir <- file.path(package.dir, 'extdata')
# qcdir <- tempdir()
qcdir <- "/usr/local/src"
name <- 'data'
# path2plink <- "/Users/hannah/bin/plink"
path2plink <- "./plink"
# Per-individual quality control
# For perIndividualQC, one simply specifies the directory where the data is
# stored (qcdir) and the prefix of the plink files (i.e. prefix.bim, prefix.bed,
# prefix.fam). In addition, the names of the files containing information about
# the reference population and the merged dataset used in check_ancestry have to
# be provided: refSamplesFile, refColorsFile and prefixMergedDataset. Per
# default, all quality control checks will be conducted.

# In addition to running each check, perIndividualQC writes a list of all fail
# individual IDs to the qcdir. These IDs will be removed in the computation of
# the perMarkerQC. If the list is not present, perMarkerQC will send a message
# about conducting the quality control on the entire dataset.

# NB: To reduce the data size of the example data in plinkQC, data.genome has
# already been reduced to the individuals that are related. Thus the relatedness
# plots in C only show counts for related individuals only.

# NB: To demonstrate the results of the ancestry check, the required eigenvector
# file of the combined study and reference datasets have been precomputed and
# for the purpose of this example will be copied to the qcdir. In practice, the
# qcdir will often be the same as the indir and this step will not be required.

file.copy(file.path(package.dir, 'extdata/data.HapMapIII.eigenvec'), qcdir)
#> [1] TRUE
# perIndividualQC displays the results of the quality control steps in a multi-panel plot.

# fail_individuals <- perIndividualQC(indir=indir, qcdir=qcdir, name=name,
#                                     refSamplesFile=paste(indir, "/HapMap_ID2Pop.txt",
#                                                          sep=""), 
#                                     refColorsFile=paste(indir, "/HapMap_PopColors.txt",
#                                                         sep=""),
#                                     prefixMergedDataset="data.HapMapIII",
#                                     path2plink=path2plink, do.run_check_ancestry = FALSE,
#                                     interactive=TRUE, verbose=TRUE)


fail_individuals <- perIndividualQC(indir=indir, qcdir=qcdir, name=name,
                                    refSamplesFile=paste(indir, "/HapMap_ID2Pop.txt",
                                                         sep=""), 
                                    refColorsFile=paste(indir, "/HapMap_PopColors.txt",
                                                        sep=""),
                                    prefixMergedDataset="data.HapMapIII",
                                    path2plink=path2plink, do.run_check_ancestry = FALSE,
                                    interactive=F, verbose=TRUE)


# overviewperIndividualQC depicts overview plots of quality control failures and
# the intersection of quality control failures with ancestry exclusion.

overview_individuals <- overviewPerIndividualQC(fail_individuals,
                                                interactive=FALSE)


# Per-marker quality control
# perMarkerQC applies its checks to data in the specified directory (qcdir),
# starting with the specified prefix of the plink files (i.e. prefix.bim,
# prefix.bed, prefix.fam). Optionally, the user can specify different thresholds
# for the quality control checks and which check to conduct. Per default, all
# quality control checks will be conducted. perMarkerQC displays the results of
# the QC step in a multi-panel plot.

fail_markers <- perMarkerQC(indir=indir, qcdir=qcdir, name=name,
                            path2plink=path2plink,
                            verbose=TRUE, interactive=TRUE,
                            showPlinkOutput=FALSE)


# overviewPerMarkerQC depicts an overview of the marker quality control failures
# and their overlaps.

overview_marker <- overviewPerMarkerQC(fail_markers, interactive=TRUE)


# Create QC-ed dataset
# After checking results of the per-individual and per-marker quality control,
# individuals and markers that fail the chosen criteria can automatically be
# removed from the dataset with cleanData, resulting in the new dataset
# qcdir/data.clean.bed,qcdir/data.clean.bim, qcdir/data.clean.fam. For
# convenience, cleanData returns a list of all individuals in the study split
# into keep and remove individuals.

Ids  <- cleanData(indir=indir, qcdir=qcdir, name=name, path2plink=path2plink,
                  verbose=TRUE, showPlinkOutput=FALSE)
# Step-by-step
# Individuals with discordant sex information
# The identification of individuals with discordant sex information helps to
# detect sample mix-ups and samples with very poor genotyping rates. For each
# sample, the homozygosity rates across all X-chromosomal genetic markers are
# computed and compared with the expected rates (typically $<$0.2 for females
# and $>$0.8 for males). For samples where the assigned sex (PEDSEX in the .fam
# file) contradicts the sex inferred from the homozygosity rates (SNPSEX), it
# should be checked that the sex was correctly recorded (genotyping often occurs
# at different locations as phenotyping and misrecording might occur). Samples
# with discordant sex information that is not accounted for should be removed
# from the study. Identifying individuals with discordant sex information is
# implemented in check_sex. It finds individuals whose SNPSEX != PEDSEX.
# Optionally, an extra data.frame with sample IDs and sex can be provided to
# double check if external and PEDSEX data (often processed at different
# centers) match. If a mismatch between PEDSEX and SNPSEX was detected, by
# SNPSEX == Sex, PEDSEX of these individuals can optionally be updated.
# check_sex depicts the X-chromosomal heterozygosity (SNPSEX) of the samples
# split by their (PEDSEX).

fail_sex <- check_sex(indir=indir, qcdir=qcdir, name=name, interactive=TRUE,
                      verbose=TRUE, path2plink=path2plink)


# Individuals with outlying missing genotype and/or heterozygosity rates
# The identification of individuals with outlying missing genotype and/or
# heterozygosity rates helps to detect samples with poor DNA quality and/or
# concentration that should be excluded from the study. Typically, individuals
# with more than 3-7% of their genotype calls missing are removed. Outlying
# heterozygosity rates are judged relative to the overall heterozygosity rates
# in the study, and individuals whose rates are more than a few standard
# deviations (sd) from the mean heterozygosity rate are removed. A typical
# quality control for outlying heterozygosity rates would remove individuals who
# are three sd away from the mean rate. Identifying related individuals with
# outlying missing genotype and/or heterozygosity rates is implemented in
# check_het_and_miss. It finds individuals that have genotyping and
# heterozygosity rates that fail the set thresholds and depicts the results as a
# scatter plot with the samples’ missingness rates on x-axis and their
# heterozygosity rates on the y-axis.

fail_het_imiss <- check_het_and_miss(indir=indir, qcdir=qcdir, name=name,
                                     interactive=TRUE, path2plink=path2plink)


# Related individuals
# Depending on the future use of the genotypes, it might required to remove any
# related individuals from the study. Related individuals can be identified by
# their proportion of shared alleles at the genotyped markers (identity by
# descend, IBD). Standardly, individuals with second-degree relatedness or
# higher will be excluded. Identifying related individuals is implemented in
# check_relatedness. It finds pairs of samples whose proportion of IBD is larger
# than the specified highIBDTh. Subsequently, for pairs of individual that do
# not have additional relatives in the dataset, the individual with the greater
# genotype missingness rate is selected and returned as the individual failing
# the relatedness check. For more complex family structures, the unrelated
# individuals per family are selected (e.g. in a parents-offspring trio, the
# offspring will be marked as fail, while the parents will be kept in the
# analysis).

# NB: To reduce the data size of the example data in plinkQC, data.genome has
# already been reduced to the individuals that are related. Thus the relatedness
# plots in C only show counts for related individuals only.

exclude_relatedness <- check_relatedness(indir=indir, qcdir=qcdir, name=name,
                                         interactive=TRUE,
                                         path2plink=path2plink)


# Individuals of divergent ancestry
# The identification of individuals of divergent ancestry can be achieved by
# combining the genotypes of the study population with genotypes of a reference
# dataset consisting of individuals from known ethnicities (for instance
# individuals from the Hapmap or 1000 genomes study [5]). Principal component
# analysis on this combined genotype panel can be used to detect population
# structure down to the level of the reference dataset (for Hapmap and 1000
# Genomes, this is down to large-scale continental ancestry). Identifying
# individuals of divergent ancestry is implemented in check_ancestry. Currently,
# check ancestry only supports automatic selection of individuals of European
# descent. It uses information from principal components 1 and 2 to find the
# center of the European reference samples. All study samples whose euclidean
# distance from the centre falls outside a specified radius are considered
# non-European. check_ancestry creates a scatter plot of PC1 versus PC2
# color-coded for samples of the reference populations and the study population.

exclude_ancestry <- check_ancestry(indir=indir, qcdir=qcdir, name=name,
                                   refSamplesFile=paste(indir, "/HapMap_ID2Pop.txt",
                                                        sep=""), 
                                   refColorsFile=paste(indir, "/HapMap_PopColors.txt",
                                                       sep=""),
                                   prefixMergedDataset="data.HapMapIII",
                                   path2plink=path2plink, run.check_ancestry = FALSE,
                                   interactive=TRUE)


# Markers with excessive missingness rate
# Markers with excessive missingness rate are removed as they are considered
# unreliable. Typically, thresholds for marker exclusion based on missingness
# range from 1%-5%. Identifying markers with high missingness rates is
# implemented in snp_missingness. It calculates the rates of missing genotype
# calls and frequency for all variants in the individuals that passed the
# perIndividualQC.

fail_snpmissing <- check_snp_missingness(indir=indir, qcdir=qcdir, name=name,
                                         interactive=TRUE,
                                         path2plink=path2plink, 
                                         showPlinkOutput=FALSE)


# Markers with deviation from HWE
# Markers with strong deviation from HWE might be indicative of genotyping or
# genotype-calling errors. As serious genotyping errors often yield very low
# p-values (in the order of 10−50), it is recommended to choose a reasonably low
# threshold to avoid filtering too many variants (that might have slight,
# non-critical deviations). Identifying markers with deviation from HWE is
# implemented in check_hwe. It calculates the observed and expected heterozygote
# frequencies per SNP in the individuals that passed the perIndividualQC and
# computes the deviation of the frequencies from Hardy-Weinberg equilibrium
# (HWE) by HWE exact test.

fail_hwe <- check_hwe(indir=indir, qcdir=qcdir, name=name, interactive=TRUE,
                      path2plink=path2plink, showPlinkOutput=FALSE)


# Markers with low minor allele frequency
# Markers with low minor allele count are often removed as the actual genotype
# calling (via the calling algorithm) is very difficult due to the small sizes
# of the heterozygote and rare-homozygote clusters. Identifying markers with low
# minor allele count is implemented in check_maf. It calculates the minor allele
# frequencies for all variants in the individuals that passed the
# perIndividualQC.

fail_maf <- check_maf(indir=indir, qcdir=qcdir, name=name, interactive=TRUE,
                      path2plink=path2plink, showPlinkOutput=FALSE)



######################################################







#######################################################








#####################################################


# Ancestry estimation
# The identification of individuals of divergent ancestry can be achieved by
# combining the genotypes of the study population with genotypes of a reference
# dataset consisting of individuals from known ethnicities (for instance
# individuals from the Hapmap or 1000 genomes study [1]). Principal component
# analysis (PCA) on this combined genotype panel can then be used to detect
# population structure down to the level of the reference dataset (for Hapmap
# and 1000 Genomes, this is down to large-scale continental ancestry).

# In the following, the workflow for combining a study dataset with the
# reference samples, conducting PCA and estimating ancestry is demonstrated. The
# study dataset consists of 200 individuals and 10,000 genetic markers and is
# provided with plinkQC in file.path(find.package('plinkQC'),'extdata').

# Workflow
# Download reference data
# A suitable reference dataset should be downloaded and if necessary,
# re-formated into PLINK format. Vignettes ‘Processing HapMap III reference data
# for ancestry estimation’ and ‘Processing 1000Genomes reference data for
# ancestry estimation’, show the download and processing of the HapMap phase III
# and 1000Genomes phase III dataset, respectively. In this example, we will use
# the HapmapIII data as the reference dataset.

# Set-up
# We will first set up some bash variables and create directories needed;
# storing the names and directories of the reference and study will make it easy
# to use updated versions of the reference or new datasets in the future. Is is
# also useful to keep the PLINK log-files for future reference. In order to keep
# the data directory tidy, we’ll create a directory for the log files and move
# them to the log directory here after each analysis step.

qcdir='~/qcdir'
refdir='~/reference'
name='data'
refname='HapMapIII'

mkdir -r $qcdir/plink_log
# Match study genotypes and reference data
# In order to compute joint principal components of the reference and study
# population, we’ll need to combine the two datasets. The plink –merge function
# enables this merge, but requires the variants in the datasets to be matching
# by chromosome, position and alleles. The following sections show how to
# extract the relevant data from the reference and study dataset and how to
# filter matching variants.

# Filter reference and study data for non A-T or G-C SNPs
# We will use an awk script to find A→T and C→G SNPs. As these SNPs are more
# difficult to align and only a subset of SNPs is required for the analysis, we
# will remove them from both the reference and study data set.

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
|| $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
$qcdir/$name.bim  > \
$qcdir/$name.ac_gt_snps

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
|| $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
$refdir/$refname.bim  > \
$qcdir/$refname.ac_gt_snps

plink --bfile  $refdir/$refname \
--exclude $qcdir/$refname.ac_gt_snps \
--make-bed \
--out $qcdir/$refname.no_ac_gt_snps
mv  $qcdir/$refname.no_ac_gt_snps.log $qcdir/plink_log/$refname.no_ac_gt_snps.log

plink --bfile  $qcdir/$name \
--exclude $qcdir/$name.ac_gt_snps \
--make-bed \
--out $qcdir/$name.no_ac_gt_snps
mv  $qcdir/$name.no_ac_gt_snps.log $qcdir/plink_log/$name.no_ac_gt_snps.log
# Prune study data
# We will conduct principle component analysis on genetic variants that are
# pruned for variants in linkage disequilibrium (LD) with an r2>0.2 in a 50kb
# window. The LD-pruned dataset is generated below, using plink –indep-pairwise
# to compute the LD-variants; additionally exclude range is used to remove
# genomic ranges of known high-LD structure. This file was originally provided
# by [6] and is available in
# file.path(find.package('plinkQC'),'extdata','high-LD-regions.txt').

plink --bfile  $qcdir/$name.no_ac_gt_snps \
--exclude range  $refdir/$highld \
--indep-pairwise 50 5 0.2 \
--out $qcdir/$name.no_ac_gt_snps
mv  $qcdir/$name.prune.log $qcdir/plink_log/$name.prune.log

plink --bfile  $qcdir/$name.no_ac_gt_snps \
--extract $qcdir/$name.no_ac_gt_snps.prune.in \
--make-bed \
--out $qcdir/$name.pruned
mv  $qcdir/$name.pruned.log $qcdir/plink_log/$name.pruned.log
# Filter reference data for the same SNP set as in study
# We will use the list of pruned variants from the study sample to reduce the
# reference dataset to the size of the study samples:
  
  plink --bfile  $refdir/$refname \
--extract $qcdir/$name.prune.in \
--make-bed \
--out $qcdir/$refname.pruned
mv  $qcdir/$refname.pruned.log $qcdir/plink_log/$refname.pruned.log
# Check and correct chromosome mismatch
# The following section uses an awk-script to check that the variant IDs of the
# reference data have the same chromosome ID as the study data. For computing
# the genetic PC, the annotation is not important, however, merging the files
# via PLINK will only work for variants with perfectly matching attributes. For
# simplicity, we update the pruned reference dataset. Note, that sex chromosomes
# are often encoded differently and might make the matching more difficult.
# Again, for simplicity and since not crucial to the final task, we will ignore
# XY-encoded sex chromosomes (via sed -n '/^[XY]/!p').

awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
($2 in a && a[$2] != $1)  {print a[$2],$2}' \
$qcdir/$name.pruned.bim $qcdir/$refname.pruned.bim | \
sed -n '/^[XY]/!p' > $qcdir/$refname.toUpdateChr

plink --bfile $qcdir/$refname.pruned \
--update-chr $qcdir/$refname.toUpdateChr 1 2 \
--make-bed \
--out $qcdir/$refname.updateChr
mv $qcdir/$refname.updateChr.log $qcdir/plink_log/$refname.updateChr.log
# Position mismatch
# Similar to the chromosome matching, we use an awk-script to find variants with
# mis-matching chromosomal positions.

awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
($2 in a && a[$2] != $4)  {print a[$2],$2}' \
$qcdir/$name.pruned.bim $qcdir/$refname.pruned.bim > \
$qcdir/${refname}.toUpdatePos
# Possible allele flips
# Unlike chromosomal and base-pair annotation, mismatching allele-annotations
# will not only prevent the plink –merge, but also mean that it is likely that
# actually a different genotype was measured. Initially, we can use the
# following awk-script to check if non-matching allele codes are a simple case
# of allele flips.

awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
$qcdir/$name.pruned.bim $qcdir/$refname.pruned.bim > \
$qcdir/$refname.toFlip
# Upate positions and flip alleles
# We use plink to update the mismatching positions and possible allele-flips
# identified above.

plink --bfile $qcdir/$refname.updateChr \
--update-map $qcdir/$refname.toUpdatePos 1 2 \
--flip $qcdir/$refname.toFlip \
--make-bed \
--out $qcdir/$refname.flipped
mv $qcdir/$refname.flipped.log $qcdir/plink_log/$refname.flipped.log
# Remove mismatches
# Any alleles that do not match after allele flipping, are identified and
# removed from the reference dataset.

awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
$qcdir/$name.pruned.bim $qcdir/$refname.flipped.bim > \
$qcdir/$refname.mismatch

plink --bfile $qcdir/$refname.flipped \
--exclude $qcdir/$refname.mismatch \
--make-bed \
--out $qcdir/$refname.clean
mv $qcdir/$refname.clean.log $qcdir/plink_log/$refname.clean.log
# Merge study genotypes and reference data
# The matching study and reference dataset can now be merged into a combined
# dataset with plink –bmerge. If all steps outlined above were conducted
# successfully, no mismatch errors should occur.

plink --bfile $qcdir/$name.pruned  \
--bmerge $qcdir/$refname.clean.bed $qcdir/$refname.clean.bim \
$qcdir/$refname.clean.fam  \
--make-bed \
--out $qcdir/$name.merge.$refname
mv $qcdir/$name.merge.$refname.log $qcdir/plink_log
# PCA on the merged data
# We can now run principal component analysis on the combined dataset using
# plink –pca which returns a .eigenvec file with the family and individual ID in
# columns 1 and 2, followed by the first 20 principal components.

plink --bfile $qcdir/$name.merge.$refname \
--pca \
--out $qcdir/$name.$reference
mv $qcdir/$name.$reference.log $qcdir/plink_log
Check ancestry
We can use the .eigenvec file to estimate the ancestry of the study samples. Identifying individuals of divergent ancestry is implemented in check_ancestry. Currently, check ancestry only supports automatic selection of individuals of European descent. It uses principal components 1 and 2 to find the center of the known European reference samples. All study samples whose Euclidean distance from the centre falls outside the radius specified by the maximum Euclidean distance of the reference samples multiplied by the chosen europeanTh are considered non-European. check_ancestry shows the result of the ancestry analysis in a scatter plot of PC1 versus PC2 colour-coded for samples of the reference populations and the study population. From within R, run the following command to the ancestry check:
  
  library(plinkQC)
indir <- system.file("extdata", package="plinkQC")
name <- 'data'
refname <- 'HapMapIII'
prefixMergedDataset <- paste(name, ".", refname, sep="")

exclude_ancestry <-
  evaluate_check_ancestry(indir=indir, name=name,
                          prefixMergedDataset=prefixMergedDataset,
                          refSamplesFile=paste(indir, "/HapMap_ID2Pop.txt",
                                               sep=""), 
                          refColorsFile=paste(indir, "/HapMap_PopColors.txt",
                                              sep=""),
                          interactive=TRUE)