UQCCG.data<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015"
project<-"2015-07-14_MODY_NSGC_DISH_PCC"
project.name <- project
analysis.dir<-paste(UQCCG.data,project,"Analysis/2015-07-14_MODY_NSGC_DISH_PCC.BY-CHR",sep="/")


UQCCG.data<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015"
project<-"2015-06-29_MODY_NSGC"
project.name <- project
analysis.dir<-paste(UQCCG.data,project,"Analysis",sep="/")

#/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-29_MODY_NSGC/Analysis/2015-06-29_MODY_NSGC.chr5.DISH.analysis.txt

#project.extension<-"analysis-maf-filtered.txt"
project.extension<-"analysis.txt"
#".analysis-maf-filtered.txt"   ".analysis-maf-filtered.txt.geno.all.txt"  ## just the exterion not fam.extension! ".analysis-maf-filtered.txt.geno.all.txt" #
fam<-c(".NSGC-FAM-49.") 


  project.files<- "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-07-14_MODY_NSGC_DISH_PCC/Analysis/2015-07-14_MODY_NSGC_DISH_PCC.NSGC-FAM-49.All-maf-filtered.txt"

project.files<- "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-29_MODY_NSGC/Analysis/2015-06-29_MODY_NSGC.NSGC-FAM-49.All-maf-filtered.txt"

    ## "/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-07-14_MODY_NSGC_DISH_PCC/Analysis/2015-07-14_MODY_NSGC_DISH_PCC.BY-CHR/2015-07-14_MODY_NSGC_DISH_PCC.chr1.NSGC-FAM-49.analysis-maf-filtered.txt"
    ## /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-07-14_MODY_NSGC_DISH_PCC/Analysis/2015-07-14_MODY_NSGC_DISH_PCC.BY-CHR/2015-07-14_MODY_NSGC_DISH_PCC.chr1.NSGC-FAM-49.analysis.txt

analysis.dir<-dirname(project.files)


    /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-07-14_MODY_NSGC_DISH_PCC/Analysis/2015-07-14_MODY_NSGC_DISH_PCC.NSGC-FAM-49.wanted.All-maf-filtered.txt
    /media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-07-14_MODY_NSGC_DISH_PCC/Analysis/2015-07-14_MODY_NSGC_DISH_PCC.NSGC-FAM-49.All-maf-filtered.txt

   
/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-06-29_MODY_NSGC/Analysis/2015-06-29_MODY_NSGC.chr1.NSGC-FAM-49.analysis.txt





code.dir<-"/mnt/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/"
setwd(code.dir)
source("ucsc.table.names.r")   # load in the UCSC tables these use the db file names and not their lable-names 
source("annotate_SNPs_subroutines.r")
source("ucsc.table.names.processor.r")





#################################################################################
#### assume has format project.chr.fam.extension or chr.project.fam.extension
setwd(analysis.dir)
getwd()
files<-dir(analysis.dir)
the.extension<-paste(project.extension,"$",sep="")
files<-files[grepl(the.extension ,files)]
files
if(fam=="ALL" | fam=="All" | fam=="all" ){
  fam<-unique(unlist( mapply(function(x){x[length(x)]}, strsplit(gsub(the.extension,"",files),split=".",fixed=TRUE)   )))
}

fam

#
ifam<-1
  
  the.extension<-paste(fam[ifam],project.extension,"$",sep="")
  project.files<-files[grepl(the.extension ,files)]
  print(sort(paste("Doing: ",project.files,sep=""))) # project.files<-project.files[1:22]
  

 
  data.summary<-{} # used if doing genotype counts
  project.files
  # ichr<-1
 
  project.files<-project.files[!grepl("chrALL",project.files)] ## remove old combined runs 

  if(length(project.files)!=24){
    print("########################################### WARNING #################################")
    print("less that 24 chromosomes detected")
    print(fam[ifam])
    print("########################################### WARNING #################################") 
  }
  project.files  #  project.files <- project.files[-23]
   ichr<-20


#############################################
  indels<-{}
  the.col<-{}


ichr<-1
project.files
setwd(analysis.dir)
getwd()
  for(ichr in 1:length(project.files)){
    print(project.files[ichr])
    

    
    ################## fast read ###########
    column.labels<-read.delim(project.files[ichr],header=F,nrows=1,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
    num.vars<-dim(column.labels)[2]
    a.indel<-scan(project.files[ichr],what=character(num.vars),skip=1,sep="\t",fill=TRUE,na.strings="",quote="") # quote="\""
    num.lines<-length(a.indel)/(num.vars)
    dim(a.indel)<-c(num.vars,num.lines)
    a.indel<-t(a.indel)
    colnames(a.indel)<-column.labels

    print(table(a.indel[,"chr"]))
    ########################################

        if(is.null(dim(indels))){
      indels<-a.indel
      the.col<-colnames(a.indel)
    }else{
      
      if(sum(!(the.col %in% colnames(a.indel)))>0  | sum(!(colnames(a.indel) %in% the.col))>0){
        print("error colnames don't match")
        print(the.col[!(the.col %in% colnames(a.indel))])
        print(colnames(a.indel)[!(colnames(a.indel) %in% the.col)])
        next
      } # columns don't match
      
      if(sum(!(colnames(a.indel) %in% the.col))>0 ){
        a.indel<-a.indel[,the.col]     } # reduce a.indels }
      
      indels<-rbind(indels,a.indel[,the.col])
    } ## is null so  not first

}


if(length(project.files)>1){
    a.indel<-indels
    rm(indels)
}
    
print(sort(table(a.indel[,"chr"])))


summary<-table(a.indel[,"chr"])

#a.indel[1:5,1:50]

table((a.indel[,"FILTER"]))
sum(grepl("^snp",a.indel[,"TYPE"]))
sum(grepl("^indel",a.indel[,"TYPE"]))


table((a.indel[,"MAF.lt:0"]))
sum(as.logical(a.indel[,"MAF.lt:0"]))
sum(as.logical(a.indel[,"MAF.lt:0.001"]))
sum(as.logical(a.indel[,"MAF.lt:0.01"]))

the.samples<-colnames(a.indel)[grepl(".GT$",colnames(a.indel))]
the.samples
## MAF.lt:0.005 MAF.lt:0.01 MAF.lt:0.025

dim(a.indel)
sort(summary)
sort(names(summary))

> dim(a.indel)
[1] 5631  248
> sort(summary)

  Y  18  21  13  20  22   X   5  14  15   4   8  16   9   6   7  10  12   2  17 
  2  70  81  82 101 120 125 134 136 144 146 150 221 223 301 309 312 322 340 348 
 19  11   3   1 
354 443 555 612 
> sort(names(summary))
 [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "21" "22"
[16] "3"  "4"  "5"  "6"  "7"  "8"  "9"  "X"  "Y"

>
project.files


#############################################################

 a.indel.maf.06<-a.indel # a.indel<-a.indel.maf.07


                         PASS  VQSRTrancheINDEL99.00to99.90 
                        17000                          2803 
VQSRTrancheINDEL99.90to100.00   VQSRTrancheSNP99.90to100.00 
                          235                          1654 
> sum(grepl("^snp",a.indel[,"TYPE"]))
[1] 8401
> sum(grepl("^indel",a.indel[,"TYPE"]))
[1] 13291
> table((a.indel[,"MAF.lt:0"]))

FALSE  TRUE 
12291  9401 
> sum(as.logical(a.indel[,"MAF.lt:0"]))
[1] 9401
> sum(as.logical(a.indel[,"MAF.lt:0.001"]))
[1] 18362
> sum(as.logical(a.indel[,"MAF.lt:0.01"]))
[1] 19094
> the.samples<-colnames(a.indel)[grepl(".GT$",colnames(a.indel))]
> the.samples
[1] "NSGC-49.3.GT"
> dim(a.indel)
[1] 21692   242
> sort(summary)

   Y   21   18   22   13   20    X   15   14    8   16    4    5    9   10   19 
  22  315  373  403  462  518  537  656  716  738  786  830  867  889  987 1077 
   7   11   17    6   12    2    3    1 
1206 1240 1243 1259 1281 1513 1638 2136 




















#############################################################

a.indel.maf.07<-a.indel # a.indel<-a.indel.maf.07

colname(


                       PASS VQSRTrancheSNP99.90to100.00 
                       6750                        1653 
> sum(grepl("^snp",a.indel[,"TYPE"]))
[1] 8403
> sum(grepl("^indel",a.indel[,"TYPE"]))
[1] 0
> table((a.indel[,"MAF.lt:0"]))

FALSE  TRUE 
 7212  1191 
> sum(as.logical(a.indel[,"MAF.lt:0"]))
[1] 1191
> sum(as.logical(a.indel[,"MAF.lt:0.001"]))
[1] 5710
> sum(as.logical(a.indel[,"MAF.lt:0.01"]))
[1] 6301
> the.samples<-colnames(a.indel)[grepl(".GT$",colnames(a.indel))]
> the.samples
[1] "NSGC-49.3.GT"
> dim(a.indel)
[1] 8403  248
> sort(summary)

  Y  18  13  21  22   X  14  20   5   4  15   8  16   9  10   6  12  19   7   2 
 11 120 135 164 183 190 201 220 229 236 238 243 317 351 421 442 452 455 470 515 
 17  11   3   1 
540 559 764 947 
> sort(names(summary))
 [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "21" "22"
[16] "3"  "4"  "5"  "6"  "7"  "8"  "9"  "X"  "Y" 
> 

##############
a.indel.maf.chr.07<-a.indel
> table((a.indel[,"FILTER"]))

                         PASS  VQSRTrancheINDEL99.00to99.90 
                        89642                         10086 
VQSRTrancheINDEL99.90to100.00   VQSRTrancheSNP99.90to100.00 
                         1052                          3414 
> sum(grepl("^snp",a.indel[,"TYPE"]))
[1] 32464
> sum(grepl("^indel",a.indel[,"TYPE"]))
[1] 71730
> table((a.indel[,"MAF.lt:0"]))

FALSE  TRUE 
48265 55929 
> sum(as.logical(a.indel[,"MAF.lt:0"]))
[1] 55929
> sum(as.logical(a.indel[,"MAF.lt:0.001"]))
[1] 84310
> sum(as.logical(a.indel[,"MAF.lt:0.01"]))
[1] 89602
> the.samples<-colnames(a.indel)[grepl(".GT$",colnames(a.indel))]
> the.samples
[1] "NSGC-49.3.GT"
> dim(a.indel)
[1] 104194    232
> sort(summary)

    Y    21    18    22    13    20    14     X    15    16     8     4     9 
  104  1256  1869  2132  2267  2730  3403  3482  3557  3741  3759  4020  4282 
   10    19     5     7    11    12     6    17     3     2     1 
 4592  4872  4960  5359  5683  5870  5958  6002  6635  7432 10229 
> sort(names(summary))
 [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "21" "22"
[16] "3"  "4"  "5"  "6"  "7"  "8"  "9"  "X"  "Y" 
> 









################################################
a.indel.07<-a.indel # a.indel<-a.indel.07



> table((a.indel[,"FILTER"]))

                         PASS  VQSRTrancheINDEL99.00to99.90 
                       492922                         21412 
VQSRTrancheINDEL99.90to100.00   VQSRTrancheSNP99.90to100.00 
                         2097                         11071 
> sum(grepl("^snp",a.indel[,"TYPE"]))
[1] 291974
> sum(grepl("^indel",a.indel[,"TYPE"]))
[1] 235528
> table((a.indel[,"MAF.lt:0"]))

 FALSE   TRUE 
344610 182892 
> sum(as.logical(a.indel[,"MAF.lt:0"]))
[1] 182892
> sum(as.logical(a.indel[,"MAF.lt:0.001"]))
[1] 273899
> sum(as.logical(a.indel[,"MAF.lt:0.01"]))
[1] 290198
> the.samples<-colnames(a.indel)[grepl(".GT$",colnames(a.indel))]
> the.samples
[1] "NSGC-49.3.GT"
> dim(a.indel)
[1] 527502    232
> sort(summary)

    Y    21    18    13    22     X    20    14    15     8    16     4     9 
  200  7671  8925 11366 12832 13954 14350 17059 17975 18493 20981 21196 21480 
   10     5     7    11    12    17     3     6    19     2     1 
22740 23587 26475 27328 29116 29614 30101 30365 33208 36798 51688 
> sort(names(summary))
 [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "21" "22"
[16] "3"  "4"  "5"  "6"  "7"  "8"  "9"  "X"  "Y" 
> a.indel.07<-a.indel # a.indel<-a.indel.07 
################################################


a.indel.06<-a.indel

identical(a.indel,a.indel.06)

> table((a.indel[,"FILTER"]))

                         PASS  VQSRTrancheINDEL99.00to99.90 
                       372583                         13830 
VQSRTrancheINDEL99.90to100.00   VQSRTrancheSNP99.90to100.00 
                         1122                         11071 
> sum(grepl("^snp",a.indel[,"TYPE"]))
[1] 291974
> sum(grepl("^indel",a.indel[,"TYPE"]))
[1] 106632


 FALSE   TRUE 
285377 113229 
> sum(as.logical(a.indel[,"MAF.lt:0"]))
[1] 113229
> sum(as.logical(a.indel[,"MAF.lt:0.001"]))
[1] 177101
> sum(as.logical(a.indel[,"MAF.lt:0.01"]))
[1] 189277
> the.samples<-colnames(a.indel)[grepl(".GT$",colnames(a.indel))]
> the.samples
[1] "NSGC-49.3.GT"
> dim(a.indel)
[1] 398606    226
> sort(summary)

    Y    21    18    13    22     X    20    14    15     8     4    16     9 
  138  5993  6725  8535 10310 10715 11374 12624 13328 13877 15473 16295 16321 
   10     5     7    11    12     3     6    17    19     2     1 
17047 17282 19830 20780 21565 22228 22709 22900 26305 26922 39330 
> sort(names(summary))
 [1] "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "21" "22"
[16] "3"  "4"  "5"  "6"  "7"  "8"  "9"  "X"  "Y" 
>
save.image("monogenic_comparison.RData")
load("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-07-14_MODY_NSGC_DISH_PCC/Analysis/2015-07-14_MODY_NSGC_DISH_PCC.BY-CHR/monogenic_comparison.RData")


sum(grepl("^snp",a.indel.06[,"TYPE"]))

sum(grepl("^snp",a.indel.07[,"TYPE"]))

sum(grepl("^indel",a.indel.06[,"TYPE"]))

sum(grepl("^indel",a.indel.07[,"TYPE"]))


sum(grepl("^snp",a.indel.maf.06[,"TYPE"]))

sum(grepl("^snp",a.indel.maf.07[,"TYPE"]))

sum(grepl("^snp",a.indel.maf.chr.07[,"TYPE"]))


dim(a.indel.maf.07)
dim(a.indel.maf.chr.07)
dim(a.indel.maf.06)

colnames(a.indel.maf.07)[grepl("^MAF.lt",colnames(a.indel.maf.07))]


sum(as.logical(a.indel.maf.07[,"MAF.lt:0"]))


a.snp<-grepl("^snp",a.indel.maf.06[,"TYPE"])
rare<-as.logical(a.indel.maf.06[,"MAF.lt:0"])


a.snp<-grepl("^snp",a.indel.maf.07[,"TYPE"])
rare<-as.logical(a.indel.maf.07[,"MAF.lt:0"])


a.snp<-grepl("^snp",a.indel.maf.chr.07[,"TYPE"])
rare<-as.logical(a.indel.maf.chr.07[,"MAF.lt:0"])

dim(a.indel.maf.07)
dim(a.indel.maf.07)
dim(a.indel.maf.07)

length(a.snp)
sum(a.snp)

length(rare)
sum(rare)

sum(a.snp & rare)

sum(as.logical(a.indel.maf.06[a.snp,"MAF.lt:0"]))
