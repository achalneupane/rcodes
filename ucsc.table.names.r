#####################################variables for 
############# Get column title from UCSC track  column labels keep just first column (ID) and the last column (description)
##dgv "field	example	SQL type 	info  description"

#### id the database name contains a ++ then use name without tthis cat++dog _> catdog as variable names can't have a ++
Gerp.labels.wanted<-c("chrom","chromStart","chromEnd","score")                   
Gerp.labels<-c("bin","Indexing field to speed chromosome range queries.",
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"score","Score from 0-1000")

segdup.labels.wanted<-c("chrom","chromStart","chromEnd","otherChrom","otherEnd","otherEnd","fracMatch")
segdup.labels<-c("bin", "Indexing field to speed chromosome range queries",
"chrom", "Reference sequence chromosome or scaffold",
"chromStart",     "Start position in chromosome",
"chromEnd","End position in chromosome",
"name",  "Other chromosome involved",
"score",  "Score based on the raw BLAST alignment score. Set to 0 and not used in later versions",
"strand", "Value should be + or -",
"otherChrom",   "Other chromosome or scaffold",
"otherStart", "Start in other sequence",
"otherEnd",  "End in other sequence",
"otherSize", 	"Total size of other sequence (otherEnd - otherStart)",
"uid",	"Unique id shared by the query and subject",
"posBasesHit", " For future use",
"testResult", " For future use",
"verdict", "For future use",
"chits", "For future use",
"ccov", "For future use",
"alignfile", "alignment file path",
"alignL", 	"spaces/positions in alignment",
"indelN", "number of indels",
"indelS", 	"indel spaces",
"alignB", "bases Aligned",
"matchB", "aligned bases that match",
"mismatchB", "aligned bases that do not match",
"transitionsB", "number of transitions",
"transversionsB", "number of transversions",
"fracMatch",	"fraction of matching bases",
"fracMatchIndel", "fraction of matching bases with indels",
"jcK", 	"K-value calculated with Jukes-Cantor",
"k2K", "Kimura K")

##dgv "field	example	SQL type 	info  description"
dgv.labels.wanted<-c("chrom","chromStart","chromEnd","varType","sample")                   
dgv.labels<-c("bin","Indexing field to speed chromosome range queries.",
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"name","Name of item",
"score","Score from 0-1000",
"strand","+ or -",
"thickStart","Same as chromStart (placeholder for BED 9+ format)",
"thickEnd","Same as chromEnd (placeholder for BED 9+ format)",
"itemRgb","Item R,G,B color.",
"landmark","Genomic marker near the variation locus",
"varType","Type of variation",
"reference","Literature reference for the study that included this variant",
"pubMedId","For linking to pubMed abstract of reference",
"method","Brief description of method/platform",
"sample","Description of sample population for the study")

gerpelem.labels.wanted<-c("chrom","chromStart","chromEnd","score")                   
gerpelem.labels<-c(
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"score","Score from 0-1000",
"p-value","Significance")

jaxQtlAsIs.labels.wanted<-c("chrom","chromStart","chromEnd","name")                   
jaxQtlAsIs.labels<-c(
"bin","Indexing field to speed chromosome range queries.",                     
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"name","Jax name")

## ljb_gerp.labels.wanted<-c("chrom","chromStart","chromEnd","score")                   
## ljb_gerp.labels<-c(
## "chrom","Reference sequence chromosome or scaffold",
## "chromStart","Start position in chromosome",
## "chromEnd","End position in chromosome",
## "score","Score from 0-1000",
## "p-value","Significance")


tfbsConsSites.labels.wanted<-c("chrom","chromStart","chromEnd","name","zScore")
tfbsConsSites.labels<-c(
"bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"name","Name of item",
"score","Score from 0-1000",
"strand","+ or -",
"zScore","zScore")

tfbs.labels.wanted<-c("chrom","chromStart","chromEnd","name","zScore")
tfbs.labels<-c(
"bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"name","Name of item",
"score","Score from 0-1000",
"strand","+ or -",
"zScore","zScore")

mce28way.labels.wanted<-c("chrom","chromStart","chromEnd","score")
mce28way.labels<-c(
"bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"name","Name of item",
"score","Score from 0-1000")

mce44way.labels.wanted<-c("chrom","chromStart","chromEnd","score")
mce44way.labels<-c(
"bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"name","Name of item",
"score","Score from 0-1000")

mce30way.labels.wanted<-c("chrom","chromStart","chromEnd","score")
mce30way.labels<-c(
"bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"name","Name of item",
"score","Score from 0-1000")


mirnatarget.labels.wanted<-c("chrom","chromStart","chromEnd","name","score")
mirnatarget.labels<-c(
"bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"name","Name of item",
"score","Score from 0-1000",
"strand","+ or -" )

mirna.labels.wanted<-c("chrom","chromStart","chromEnd","name","score","type")
mirna.labels<-c(
"bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"name","Name of item",
"score","Score from 0-1000",
"strand","+ or -",
"thickStart","Start codon",
"thickEnd","stop codon",
"type","Types of RNA")

miRNA.labels.wanted<-c("chrom","chromStart","chromEnd","name")
miRNA.labels<-c(
"bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"name","Name of item",
"score","Score from 0-1000",
"strand","+ or -"
)

jaxPhenotype.labels.wanted<-c("chrom","chromStart","chromEnd","source")
jaxPhenotype.labels<-c(
"bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"name","Name of item",
"score","Score from 0-1000",
"strand","+ or -",
"thickStart","Start codon",
"thickEnd","stop codon",
"reserved","Used as itemRgb as of 2004-11-22",
"blockCount","Number of blocks",
"blockSizes","Comma separated list of block sizes",
"chromStarts","Start positions relative to chromStart",
"source","source of item" )


rnaGene.labels.wanted<-c("chrom","chromStart","chromEnd","name","score","type")
rnaGene.labels<-c(
"bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"name","Name of item",
"score","Score from 0-1000",
"strand","+ or -",
"source","detection method",
"type","Types of RNA",
"fullScore","Eddies Score",
"isPseudo","Is Pseudo gene" )

rgdQtl.labels.wanted<-c("chrom","chromStart","chromEnd","name")
rgdQtl.labels<-c(
"bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"name","Name of item")

## decription fro above
rgdQtlLink.labels.wanted<-c("name","description")
rgdQtlLink.labels<-c(
"id","RGD QTL ID",
"name","symbolic name",
"description","QTL description")

rgdQtl.labels.wanted<-c("chrom","chromStart","chromEnd","name")
rgdQtl.labels<-c(
"bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"name","Name of item")

#### Encode require -scorecolm 5
cpgIslandExt.labels.wanted<-c("chrom","chromStart","chromEnd","length","name")
cpgIslandExt.labels<-c(
"bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"name","Name of item",
"length","Island Length",
"cpgNum","Number of CpGs in island",
"gcNum","Number of C and G in island",
"perCpg","Percentage of island that is CpG",
"perGc","Percentage of island that is C or G",
"obsExp","Ratio of observed(cpgNum) to expected(numC*numG/length) CpG in island")

wgEncodeRegDnaseClustered.labels.wanted<-c("chrom","chromStart","chromEnd","name","score")
wgEncodeRegDnaseClustered.labels<-c(
"bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"name","Name of item",
"score","Score from 0-1000")


wgEncodeRegTfbsClustered.labels.wanted<-c("chrom","chromStart","chromEnd","name","score")
wgEncodeRegTfbsClustered.labels<-c(
"bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"name","Name of item",
"score","Score from 0-1000",
"strand","+ or -",
"thickStart","Start codon",
"thickEnd","stop codon",
"reserved","Used as itemRgb as of 2004-11-22",
"blockCount","Number of blocks",
"blockSizes","Comma separated list of block sizes",
"chromStarts","Start positions relative to chromStart",
"expCount","Number of experiment values",
"expIds","Comma separated list of experiment IDs",
"expScores","Comma separated list of experiment scores" )


gwasCatalog.labels.wanted<-c("name","trait","pValue","orOrBeta","cnv")
gwasCatalog.labels<-
c("bin","Indexing field to speed chromosome range queries.",
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"name","ID of SNP associated with trait",
"pubMedID","PubMed ID of publication of the study",
"author","First author of publication",
"pubDate","Date of publication",
"journal","Journal of publication",
"title","Title of publication",
"trait","Disease or trait assessed in study",
"initSample","Initial sample size",
"replSample","Replication sample size",
"region","Chromosome band / region of SNP",
"genes","Reported Gene(s)",
"riskAllele","Strongest SNP-Risk Allele",
"riskAlFreq","Risk Allele Frequency",
"pValue","p-Value",
"pValueDesc","p-Value Description",
"orOrBeta","Odds ratio or beta",
"ci95","95% Confidence Interval",
"platform","Platform and [SNPs passing QC]",
"cnv","Copy Number Variant")

## snp132.labels<-c("bin","Indexing field to speed chromosome range queries",
## "chrom","Reference sequence chromosome or scaffold",
## "chromStart","Start position in chrom",
## "chromEnd","End position in chrom",
## "name","dbSNP Reference SNP (rs) identifier",
## "score","Not used",
## "strand","Which DNA strand contains the observed alleles",
## "refNCBI","Reference genomic sequence from dbSNP",
## "refUCSC","Reference genomic sequence from UCSC lookup of chrom,chromStart,chromEnd",
## "observed","The sequences of the observed alleles from rs-fasta files (C/T)",
## "molType","Sample type from exemplar submitted SNPs ('unknown', 'genomic', 'cDNA')",
## "class","('unknown', 'single', 'in-del', 'het', 'microsatellite', 'named', 'mixed', 'mnp', 'insertion', 'deletion'",
## "valid","('unknown', 'by-cluster', 'by-frequency', 'by-submitter', 'by-2hit-2allele', 'by-hapmap', 'by-1000genomes') 	Validation status of the SNP",
## "avHet","Average heterozygosity from all observations. Note: may be computed on small number of samples",
## "avHetSE","Standard Error for the average heterozygosity",
## "func","Functional category of the SNP (coding-synon, coding-nonsynon, intron, etc.)",
## "locType","Type of mapping inferred from size on reference; may not agree with class",
## "weight","The quality of the alignment: 1 = unique mapping, 2 = non-unique, 3 = many matches",
## "exceptions","Unusual conditions noted by UCSC that may indicate a problem with the data",
## "submitterCount","Number of distinct submitter handles for submitted SNPs for this ref SNP",
## "submitters","List of submitter handles",
## "alleleFreqCount","Number of observed alleles with frequency data",
## "alleles","Observed alleles for which frequency data are available",
## "alleleNs","Count of chromosomes (2N) on which each allele was observed. Note: this is extrapolated by dbSNP from submitted frequencies and total sample 2N, and is not always an integer",
## "alleleFreqs","Allele frequencies",
## "bitfields","( clinically-assoc ,  maf-5-some-pop ,  maf-5-all-pops ,  has-omim-omia ,  microattr-tpa ,  submitted-by-lsdb ,  genotype-conflict ,  rs-cluster-nonoverlapping-alleles ,  observed-mismatch ) 	SNP attributes extracted from dbSNP s SNP_bitfield table")

snp130.labels.wanted<-c("name","observed","class","func","avHet","avHetSE")
snp130.labels<-c("bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStar","Start position in chrom",
"chromEnd","End position in chrom",
"name","dbSNP Reference SNP identifier",
"score","Not used",
"strand","Which DNA strand contains the observed alleles",
"refNCBI","Reference genomic sequence from dbSNP",
"refUCSC","Reference genomic sequence from UCSC lookup of chrom,chromStart,chromEnd",
"observed","The sequences of the observed alleles from rs-fasta files",
"molType","Sample type from exemplar submitted sequence (ss)",
"class","Class of variant (single, in-del, named, mixed, etc.)",
"valid","Validation status of the SNP",
"avHet","Average heterozygosity from all observations",
"avHetSE","Standard Error for the average heterozygosity",
"func","Functional category of the SNP (coding-synon, coding-nonsynon, intron, etc.)",
"locType","Type of mapping inferred from size on reference; may not agree with class",
"weight","The quality of the alignment: 1 = unique mapping, 2 = non-unique, 3 = many matches")


snp131.labels.wanted<-c("name","observed","class","func","avHet","avHetSE")
snp131.labels<-c("bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStar","Start position in chrom",
"chromEnd","End position in chrom",
"name","dbSNP Reference SNP identifier",
"score","Not used",
"strand","Which DNA strand contains the observed alleles",
"refNCBI","Reference genomic sequence from dbSNP",
"refUCSC","Reference genomic sequence from UCSC lookup of chrom,chromStart,chromEnd",
"observed","The sequences of the observed alleles from rs-fasta files",
"molType","Sample type from exemplar submitted sequence (ss)",
"class","Class of variant (single, in-del, named, mixed, etc.)",
"valid","Validation status of the SNP",
"avHet","Average heterozygosity from all observations",
"avHetSE","Standard Error for the average heterozygosity",
"func","Functional category of the SNP (coding-synon, coding-nonsynon, intron, etc.)",
"locType","Type of mapping inferred from size on reference; may not agree with class",
"weight","The quality of the alignment: 1 = unique mapping, 2 = non-unique, 3 = many matches")


#####this is for snp132 on hg18 downloaded from annovar
snp132.labels.wanted<-c("name","observed","A1_A2","class","func","avHet","avHetSE")
snp132.labels<-c("bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStar","Start position in chrom",
"chromEnd","End position in chrom",
"name","dbSNP Reference SNP identifier",
"score","Not used",
"strand","Which DNA strand contains the observed alleles",
"A1","allele 1 can be long",
"A2","allele 2 can be long",
"A1_A2","A1 and A2",
"molType","Sample type from exemplar submitted sequence (ss)",
"class","Class of variant (single, in-del, named, mixed, etc.)",
"Submitted","Submitter",
"avHet","Average heterozygosity from all observations",
"avHetSE","Standard Error for the average heterozygosity",
"func","Functional category of the SNP (coding-synon, coding-nonsynon, intron, etc.)",
"locType","Type of mapping inferred from size on reference; may not agree with class",
"weight","The quality of the alignment: 1 = unique mapping, 2 = non-unique, 3 = many matches")


#####this is for snp132 on hg18 downloaded from annovar
snp132.labels.wanted<-c("name","observed","A1_A2","class","func","avHet","avHetSE")
snp132.labels<-c("bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStar","Start position in chrom",
"chromEnd","End position in chrom",
"name","dbSNP Reference SNP identifier",
"score","Not used",
"strand","Which DNA strand contains the observed alleles",
"A1","allele 1 can be long",
"A2","allele 2 can be long",
"A1_A2","A1 and A2",
"molType","Sample type from exemplar submitted sequence (ss)",
"class","Class of variant (single, in-del, named, mixed, etc.)",
"Submitted","Submitter",
"avHet","Average heterozygosity from all observations",
"avHetSE","Standard Error for the average heterozygosity",
"func","Functional category of the SNP (coding-synon, coding-nonsynon, intron, etc.)",
"locType","Type of mapping inferred from size on reference; may not agree with class",
"weight","The quality of the alignment: 1 = unique mapping, 2 = non-unique, 3 = many matches")

#####this is for snp128 for mm9
snp137.labels.wanted<-c("name","observed","refNCBI","refUCSC","class","func","avHet","avHetSE")
snp137.labels<-c("bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStar","Start position in chrom",
"chromEnd","End position in chrom",
"name","dbSNP Reference SNP identifier",
"score","Not used",
"strand","Which DNA strand contains the observed alleles",
"refNCBI","NCBI allele",
"refUCSC","UCSC allele",
"observed","Alleles",
"molType","Sample type from exemplar submitted sequence (ss)",
"class","Class of variant (single, in-del, named, mixed, etc.)",
"valid","type",
"avHet","Average heterozygosity from all observations",
"avHetSE","Standard Error for the average heterozygosity",
"func","Functional category of the SNP (coding-synon, coding-nonsynon, intron, etc.)",
"locType","Type of mapping inferred from size on reference; may not agree with class",
"weight","The quality of the alignment: 1 = unique mapping, 2 = non-unique, 3 = many matches")

#####this is for snp128 for mm9
snp138.labels.wanted<-c("name","observed","refNCBI","refUCSC","class","func","avHet","avHetSE")
snp138.labels<-c("bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStar","Start position in chrom",
"chromEnd","End position in chrom",
"name","dbSNP Reference SNP identifier",
"score","Not used",
"strand","Which DNA strand contains the observed alleles",
"refNCBI","NCBI allele",
"refUCSC","UCSC allele",
"observed","Alleles",
"molType","Sample type from exemplar submitted sequence (ss)",
"class","Class of variant (single, in-del, named, mixed, etc.)",
"valid","type",
"avHet","Average heterozygosity from all observations",
"avHetSE","Standard Error for the average heterozygosity",
"func","Functional category of the SNP (coding-synon, coding-nonsynon, intron, etc.)",
"locType","Type of mapping inferred from size on reference; may not agree with class",
"weight","The quality of the alignment: 1 = unique mapping, 2 = non-unique, 3 = many matches")

#####this is for snp128 for mm9
snp138NonFlagged.labels.wanted<-c("name","observed","refNCBI","refUCSC","class","func","avHet","avHetSE")
snp138NonFlagged.labels<-c("bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStar","Start position in chrom",
"chromEnd","End position in chrom",
"name","dbSNP Reference SNP identifier",
"score","Not used",
"strand","Which DNA strand contains the observed alleles",
"refNCBI","NCBI allele",
"refUCSC","UCSC allele",
"observed","Alleles",
"molType","Sample type from exemplar submitted sequence (ss)",
"class","Class of variant (single, in-del, named, mixed, etc.)",
"valid","type",
"avHet","Average heterozygosity from all observations",
"avHetSE","Standard Error for the average heterozygosity",
"func","Functional category of the SNP (coding-synon, coding-nonsynon, intron, etc.)",
"locType","Type of mapping inferred from size on reference; may not agree with class",
"weight","The quality of the alignment: 1 = unique mapping, 2 = non-unique, 3 = many matches")

omim.desc.labels<-c("disease","genes","omim::name","band") ## complementary table hg19_omimMorbidMap.txt that describes diseases
omimGene.labels.wanted<-c("name")
omimGene.labels<-c("bin","Indexing field to speed chromosome range queries",
"chrom","Reference sequence chromosome or scaffold",
"chromStart","Start position in chromosome",
"chromEnd","End position in chromosome",
"name","Name of item")





