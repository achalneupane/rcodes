# PUBLICATION: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5525021/

###R function to implement the clustering algorithm### pos(numeric) is a vector
#containing the variant positions of the respective region sorted in an
#increasing order names(character) is a vector with the corresponding variant
#names R-function to implement the clustering algorithm described in Fier et al.
#(2017), On the association analysis of genome-sequencing data: A spatial
#clustering approach for partitioning the entire genome into non- overlapping
#windows, Genetic Epidemiology




# From a biological point of view, there are several arguments, why variants in
# a distinct genetic region could be part of the same disease mechanism: First,
# variants within the same protein functional domain are likely to be located in
# close proximity in the DNA sequence and could have similar impact on disease
# risk (Krebs, Goldstein, & Kilpatrick, 2014). Second, gene regulatory elements
# tend to cluster in certain genomic locations, such as the promoter regions
# (Raab & Kamakaka, 2010). Finally, most recently it has been shown that local
# recombination rates affect the accumulation of disease-associated mutations
# along the chromosome (Hussin et al., 2015). Given these findings, the
# identification of groups of variants that are closely located to each other
# (i.e., cluster) and similarly affect the disease susceptibility is essential
# for WGS analysis, as suitable sets of variants are required for region-based
# association analyses (Fier et al., 2012, Lee et al., 2014, Lin & Tang, 2011,
# Wu et al., 2011)

# Due to the likely presence of linkage disequilibrium (LD) between single
# nucleotide polymorphisms (SNPs) in one region, it is important to use a
# suitable association test that accounts for LD between the SNPs, for example,
# for family data the region-based rare variant family-based association test
# (RV-FBAT) (De, Yip, Ionita-Laza, & Laird, 2013), and for case-control data
# Sequence Kernel Association Test (SKAT) (Ionita-Laza, Lee, Makarov, Buxbaum, &
# Lin, 2013, Wu et al., 2011).

#window-size(numeric) and overlap(numeric) are two parameters of the algorithm , specifying the number of variants in the sliding windows and the overlap. Default is 100 for window-size and 50 for overlap
#return_what (“names”,”positions”): Either returns the “positions” or the “names” of the variants in the clusters, default is “names”
#write_output(TRUE/FALSE): Specifies whether the output should be written in a file, default=TRUE
#path_output(name of output file): The name of the output file

#####################################################################################
# names <- MAP$names
# pos <- MAP$pos
# window_size=100
# overlap=50
# return_what="names"
  
cluster<-function(pos,names,window_size=100,overlap=50,return_what="names", write_output=TRUE, path_output=NA){
    
    if (!(return_what %in% c("names","positions"))){
      stop('Please enter "names" OR "positions" as argument for return_what')
    }
    
    if (write_output & is.na(path_output)){
      stop('Please specify a valid output name in path_output','\n')
    }
    
    cat('Using window size == ',window_size,'\n')
    
    require('zoo')
    
    b<-as.numeric(pos)
    comb<-cbind(b,names)
    comb<-comb[order(b),]
    b<-as.numeric(comb[,1])
    names<-comb[,2]
    bdiff<-diff(b)
    
    # calculate lambda under the assumption that distances follow exp distribution
    
    lambda<-log(2)/median(bdiff)
    
    br<-sort(b,decreasing=T)
    brdiff<-abs(diff(br))
    
    #observed mean values in the windows
    
    ra<-rollapply(bdiff,100,mean,by=50,align="left")
    
    #theoretical mean value given medians in the windows
    ws<-as.numeric(window_size)
    sh<-as.numeric(overlap)
    
    ra2<-rollapply(bdiff,ws,function(x) 1/(log(2)/median(x)),by=sh,align="left")
    
    #actual positions
    rap<-rollapply(b,ws+1,function(x) x,by=sh,align="left")
    
    maximas<-vector("list",length(ra))
    for (i in 1:length(ra)){
      if (ra[i]<=ra2[i]) {
        maximas[[i]]<-NA
      } else {
        ord<-order(diff(rap[i,]),decreasing=T)
        j<-1
        excl<-mean(diff(rap[i,])[-c(ord[1])])
        while (ra2[i]<excl){
          j<-j+1
          excl<-mean(diff(rap[i,])[-c(ord[1:j])])
        }
        j<-j-1
        maximas[[i]]<-rap[i,][c(ord[1:j])+1]
      }
    }
    maxima<-unlist(maximas)
    maxima<-maxima[!is.na(maxima)]
    maximaf<-sort(unique(maxima))
    
    
    if (!is.na(path_output)) file.create(path_output)
    if (return_what=="positions") {
      lo_final<-split(b,findInterval(b,maximaf))
      ll<-as.numeric(lapply(lo_final,length))
      if (return_what=="positions" & write_output==TRUE) lapply(lo_final, write, path_output, append=T, ncolumns=max(ll) )
    } else {
      lo_final<-split(names,findInterval(b,maximaf))
      ll<-as.numeric(lapply(lo_final,length))
      if (return_what=="names" & write_output==TRUE)  lapply(lo_final, write, path_output, append=T, ncolumns=max(ll) )
    }
    
    return(lo_final)
    
  }
  
#####################################################################################
    
setwd("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/FBAT")    
## MAF Filtered
# ALL_MAP <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/FBAT/1803_tanzi_ped_recode12_for_fbat.map", header = FALSE)    
## No MAF filtered
ALL_MAP <- read.table("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/FBAT/1803_tanzi_ped_spatial_clustering_all_variants_recode12_for_fbat.map", header = FALSE)    

ALL_MAP <- ALL_MAP[c("V1", "V4")]
colnames(ALL_MAP) <- c("names", "pos")
ALL_MAP$names <- as.character(ALL_MAP$names)

## RUN on each CHR
CHR_pattern <- paste0("^",c(1:22, "X", "Y"),":")

for (i in 1:length(CHR_pattern)){
CHR<- gsub("\\^|:", "", CHR_pattern[i])
MAP <- ALL_MAP[grepl(CHR_pattern[i], ALL_MAP$name),]
# dim(MAP)
# head(MAP)
print(paste0("Doing::CHR", CHR))
## MAF filtered
# OUTFILE <- paste0("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/FBAT/fbat_spatial_clusters_discovery_set/", "1803_Spatial_clustering_window_100_overlap_50_CHR_", CHR,".txt")
## No MAF (All variants)
OUTFILE <- paste0("/40/AD/AD_Seq_Data/05.-Analyses/06-Aquilla_202101/09-Tanzi-replication/01-familial/03-PLINK-QC-files/FBAT/fbat_spatial_clusters_discovery_set/", "1803_Spatial_clustering_window_100_overlap_50_CHR_all_variants_", CHR,".txt")
cluster(pos= MAP$pos,names = MAP$names, window_size=100,overlap=50,return_what="names", write_output=TRUE, path_output=OUTFILE)
}






