



data<-summary ### summary could be a subset of a.indel - . a.indel[wanted,core.ann.gr)
 #'elementMetadata' cannot use "seqnames", "ranges", "strand", "seqlengths", "start", "end", "width", or "element" as column names
data.gr<-GRanges(seqnames =data[,"chr"],ranges = IRanges(start=as.numeric(data[,"start"]),end=as.numeric(data[,"end"])),strand=data[,"strand"])

reserved.gr.names<-c("seqnames", "ranges", "strand", "seqlengths", "isCircular", "start", "end", "width", "element")
names.have<-colnames(data)[!(colnames(data) %in% core.ann.gr)]
names.have[names.have %in% reserved.gr.names]<-paste(names.have[names.have %in% reserved.gr.names],1,sep=".")
colnames(data)[!(colnames(data) %in% core.ann.gr)]<-names.have
values(data.gr)<-data[,(!(colnames(data) %in% core.ann.gr)) & !remove.cols.from.values ]

which<-  data.gr
which


setwd(bam.dir)

files<-dir(getwd())
files[grepl("bam",files)]
files[match(bam.file,files)]  ### check  bam file must exist in specified locations


params<-ScanBamParam(which=which,flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=NA,isNotPassingQualityControls=FALSE),simpleCigar = FALSE,reverseComplement = FALSE,what=c("qname","flag","rname","seq","strand","pos","qwidth","cigar","qual","mapq") )  ### NOTE isValidVendorRead=FALSE shoudl be TRUE

## params<-ScanBamParam(which=which,flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=NA,isValidVendorRead=FALSE),simpleCigar = FALSE,reverseComplement = FALSE,what=c("seq","strand","pos","cigar") )
params
which
## params<-ScanBamParam(which=which,flag=scanBamFlag(isPaired=TRUE,isUnmappedQuery=FALSE,hasUnmappedMate=FALSE,isFirstMateRead=TRUE),what=c("qname","flag","strand","pos","qwidth","cigar","mrnm","mpos","isize") )

param.pile<-PileupParam(max_depth=2500, min_base_quality=10, min_mapq=13,
min_nucleotide_depth=1, min_minor_allele_depth=0,
distinguish_strands=TRUE, distinguish_nucleotides=TRUE,
ignore_query_Ns=TRUE, include_deletions=TRUE, include_insertions=FALSE,
cycle_bins=numeric() )

test<-pileup(bam.file,scanBamParam=params,pileupParam=param.pile)
