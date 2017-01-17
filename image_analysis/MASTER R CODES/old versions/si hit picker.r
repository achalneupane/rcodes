############################################### Ca screen 2 " same as above but the very latest
setwd( "/media/Bioinform-D/Research/Cellomics/Ca screen/Ca screen 2")
load("Ca2_ann.RData")  ### renames the SUMMARY files using the map file
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","G1","G2","S","Gr-S+G2/G1","Gr-S+G2+ >4N/G1","Gr-<2N","Gr->4N","Gr-G1","Gr-G2","Gr-S","NtGr-S+G2/G1","NtGr-S+G2+ >4N/G1","NtGr-<2N","NtGr->4N","NtGr-G1","NtGr-G2","NtGr-S" )


the.screen<-"Ca siRNA2"
files<- paste("plate_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Ca2_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Ca siRNA2.NOTGREEN3.Wells.txt"
well.type<-384
row.type<-16
col.type<-24

map.file<-"cellomics.to.expt.map.csv" #
exported.well.data<-"well_based_data.txt"
normalized.file<-paste(the.screen,"NORMALIZED","txt",sep=".")

targets<-c("percent_R", "ak_leakage","resazurin","Num-cells"  )
## low.cut<-c(23.3,84,0.54,1200)# version 2 not used
low.cut<-c(21,84,0.33,760) 
## high.cut<-c(41,360,1.10,3500)# version 2 not used
high.cut<-c(41,550,1.10,3500)
names(low.cut)<-targets
names(high.cut)<-targets
#######################################################
## [1] "Rep hits low cut : percent_R : 154"
## [1] "Rep hits high cut : percent_R : 0"
## [1] "Rep hits low cut : ak_leakage : 8"
## [1] "Rep hits high cut : ak_leakage : 93"
## [1] "Rep hits low cut : resazurin : 96"
## [1] "Rep hits high cut : resazurin : 24"
## [1] "Rep hits low cut : Num-cells : 93"
## [1] "Rep hits high cut : Num-cells : 17"

###############################################Ca screen 3  ### latest DNA and flexible annotation
setwd("/media/Bioinform-D/Research/Cellomics/Ca screen/Ca screen 3")
load("Ca3_ann.RData")
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann


core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","S")

the.screen<-"Ca siRNA3"
files<- paste("plate_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Ca3_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Ca siRNA3.NOTGREEN3.Wells.txt" # name now after adding well data
well.type<-384
row.type<-16
col.type<-24

map.file<-"cellomics.to.expt.map.csv" #
exported.well.data<-"well_based_data.txt"
normalized.file<-paste(the.screen,"NORMALIZED","txt",sep=".")

targets<-c("percent_R", "ak_leakage","resazurin","Num-cells"  )
## low.cut<-c(39.6,890,11500,600) ## version 2 not used
low.cut<-c(37,890,10300,600) # first sent to Greg
## high.cut<-c(52,1380,32000,1500) ## version 2 not used
high.cut<-c(52,1380,32000,1500) # first sent to Greg
names(low.cut)<-targets
names(high.cut)<-targets
#######################################################
## [1] "Rep hits low cut : percent_R : 219"
## [1] "Rep hits high cut : percent_R : 2"
## [1] "Rep hits low cut : ak_leakage : 2"
## [1] "Rep hits high cut : ak_leakage : 195"
## [1] "Rep hits low cut : resazurin : 15"
## [1] "Rep hits high cut : resazurin : 2"
## [1] "Rep hits low cut : Num-cells : 10"
## [1] "Rep hits high cut : Num-cells : 2"
############################################### millian  ### latest DNA and flexible annotation
setwd("/media/Bioinform-D/Research/Cellomics/millian")
load("siFB_ann.RData")  ### "plate" must match "plate" in file  cellomics to expt map.csv
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","G1","G2","S","Gr-S+G2/G1","Gr-S+G2+ >4N/G1","Gr-<2N","Gr->4N","Gr-G1","Gr-G2","Gr-S","NtGr-S+G2/G1","NtGr-S+G2+ >4N/G1","NtGr-<2N","NtGr->4N","NtGr-G1","NtGr-G2","NtGr-S"
             )

the.screen<-"siFB"
files<- paste("plate_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"siFB_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"siFB.NOTGREEN3.Wells.txt"
well.type<-384
row.type<-16
col.type<-24
exported.well.data<-"well_based_data.txt"

normalized.file<-paste(the.screen,"NORMALIZED","txt",sep=".")


targets<-c("S+G2_over_G1","percent_R", "ak_leakage","resazurin","Num-cells")

##"S+G2_over_G1" > 40 look crappy
low.cut<-c(1.05,42,3600,25500,4700)
high.cut<-c(2.5,65.2,8600,53000,7900)
names(low.cut)<-targets
names(high.cut)<-targets
#######################################################






###for single color siRNA
pos.control.list<-c("PLK1")
norm.control.list<-c("NT","LAMA","GFP","ATP2C1","No Lipid") # controls that could be used for  normalization

#### new annotation way
normalized<-read.delim(normalized.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
labels<-c(core.ann,core.vars,colnames(ann)[!place.core.ann])## 
colnames(normalized)[1:length(labels)]<-labels 
plates.all<-normalized
labels<-gsub("%","percent_",colnames(plates.all))
labels<-gsub("/","_over_",labels)
colnames(plates.all)<-labels
colnames(plates.all)<-gsub("/","_",colnames(plates.all))

the.screen
colnames(plates.all)
keep<-c("plate","row","col","Symbol","percent_R","ak_leakage","resazurin","Num-cells","S+G2_over_G1","S+G2+ >4N_over_G1","<2N",">4N","MEAN_Area","SD_Area","MEAN_P2A","SD_P2A","MEAN_LWR","SD_LWR")
norm<-plates.all[,keep]
targets<-c("percent_R", "ak_leakage","resazurin","Num-cells"  )

hits.low<-matrix(FALSE,nrow=dim(norm)[1],ncol=length(targets))
hits.high<-matrix(FALSE,nrow=dim(norm)[1],ncol=length(targets))
colnames(hits.low)<-targets
colnames(hits.high)<-targets
for(i in 1:length(targets)){
  the.target<-targets[i]
print(the.target)
  hits.low[,the.target]<-(as.numeric(norm[,the.target]) <= low.cut[the.target]) & !is.na(norm[,the.target])
  hits.high[,the.target]<-(as.numeric(norm[,the.target]) >= high.cut[the.target]) & !is.na(norm[,the.target])

print(paste("low cut",the.target,sum(hits.low[,the.target]),sep=" : "))
print(paste("high cut",the.target,sum(hits.high[,the.target]),sep=" : "))
}

a.screen<-{}
hits.low.rep<-matrix(0,nrow=dim(norm)[1],ncol=length(targets))
hits.high.rep<-matrix(0,nrow=dim(norm)[1],ncol=length(targets))
colnames(hits.low.rep)<-targets
colnames(hits.high.rep)<-targets
for(i in 1:length(targets)){
 the.target<-targets[i]
 
 low.counts<- sort(tapply(norm[hits.low[,the.target],"Symbol"],norm[hits.low[,the.target],"Symbol"],length))
 high.counts<- sort(tapply(norm[hits.high[,the.target],"Symbol"],norm[hits.high[,the.target],"Symbol"],length))


  posns<-apply(as.matrix(names(low.counts)),1,function(x) grep(paste("^",x,"$",sep=""),norm[,"Symbol"]))
   for(ii in 1:length(posns)){
     x<-posns[[ii]]
     hits.low.rep[x,the.target]<-low.counts[ii]
   }

   posns<-apply(as.matrix(names(high.counts)),1,function(x) grep(paste("^",x,"$",sep=""),norm[,"Symbol"]))
   for(ii in 1:length(posns)){
     x<-posns[[ii]]
     hits.high.rep[x,the.target]<-high.counts[ii]
   }
 
 
 a.test<-list(low.counts=low.counts,high.counts=high.counts)
 a.screen[i]<- list(a.test)

 
 print(paste("Rep hits low cut",the.target,sum(low.counts>1),sep=" : "))
 print(paste("Rep hits high cut",the.target,sum(high.counts>1),sep=" : "))
}
names(a.screen)<-targets
lapply(a.screen,function(x) lapply(x,length))
## hits.rep<-matrix(0,nrow=dim(norm)[1],ncol=length(targets))
## colnames(hits.rep)<-targets

###### fix cases wher genes in both extremes often a spacial defect majority rules
for(i in 1:length(targets)){
 the.target<-targets[i]
 high.variation<-(hits.high.rep[,the.target] !=0 & hits.low.rep[,the.target] !=0)

 print(sum(hits.high.rep[,the.target] !=0 & hits.low.rep[,the.target] !=0))
 
 even.call<-hits.high.rep[high.variation,the.target] == hits.low.rep[high.variation,the.target] # kill both
 hits.high.rep[high.variation,][even.call,the.target]<-0
 hits.high.rep[high.variation,][even.call,the.target]<-0

  high.call<-hits.high.rep[high.variation,the.target] > hits.low.rep[high.variation,the.target]
  hits.low.rep[high.variation,][high.call,the.target]<-0 # kill low

  low.call<-hits.high.rep[high.variation,the.target] < hits.low.rep[high.variation,the.target]
  hits.high.rep[high.variation,][low.call,the.target]<-0 # kill high
 
 
 print(sum(hits.high.rep[,the.target] !=0 & hits.low.rep[,the.target] !=0))
}
## ## test<-(hits.high.rep[,the.target] !=0 & hits.low.rep[,the.target] !=0)

####ca2 examples
## hits.low.rep[high.variation,][low.call,]
## hits.high.rep[high.variation,][low.call,]
## norm[high.variation,][low.call,]  ##IGSF2 plate 3.* 2 low one high  GNB5 two low one high

## hits.low.rep[high.variation,][high.call,]
## hits.high.rep[high.variation,][high.call,]
## norm[high.variation,][high.call,]  ##none

## norm[test,]
## cbind(hits.high.rep[test,],hits.low.rep[test,])
hits.rep<-hits.high.rep+hits.low.rep  ###case all fixed none in both extremes
colnames(hits.rep)<-paste("Reps.",colnames(hits.rep),sep="")



class<-matrix("x",nrow=dim(norm)[1],ncol=length(targets))
colnames(class)<-targets
for(i in 1:length(targets)){
 the.target<-targets[i]
 if(sum(hits.high[,the.target] & hits.low[,the.target])>0){print(paste("error in ", the.target,sep=""))}
 class[hits.low[,the.target],the.target]<-"-"
 class[hits.high[,the.target],the.target]<-"+"
}

class.string<-apply(class,1,toString)

summary<-cbind(norm,hits.rep,class.string)

colnames(summary)[dim(summary)[2]]<-"per.R,AK,Rez,Num.cells"
summary[1:5,]
dim(summary)
write.table(summary,paste(the.screen,"SUMMARY3","txt",sep="."),row.names=FALSE,sep="\t")

the.screen
a.screen


############ choose one 
a.screen.ca2<-a.screen
summary.ca2<-summary
hits.high.rep.ca2<-hits.high.rep
hits.low.rep.ca2<-hits.low.rep
#save.image("working3.ca2.RData")
setwd( "/media/Bioinform-D/Research/Cellomics/Ca screen/Ca screen 2")
load("working3.ca2.RData")

a.screen.ca3<-a.screen
summary.ca3<-summary
hits.high.rep.ca3<-hits.high.rep
hits.low.rep.ca3<-hits.low.rep

#save.image("working3.ca3.RData")
setwd( "/media/Bioinform-D/Research/Cellomics/Ca screen/Ca screen 3")

###use the rep hit status to see if is a good hits
### collect all good hits

all.hits<-{}
a.hit.ca2.high<-apply(hits.high.rep.ca2,1,function(x) sum(x>1)>0)
a.hit.ca2.low<-apply(hits.low.rep.ca2,1,function(x) sum(x>1)>0)
a.hit.ca2<-a.hit.ca2.high | a.hit.ca2.low
all.hits<-unique(summary.ca2[a.hit.ca2,"Symbol"])

##   [1] "No Cells" "OCRL"     "EDN2"     "ADRA2A"   "CCR7"     "CCR10"    "PLK1"     "GRK4"     "empty"    "EGFR"    
##  [11] "IRS1"     "SYK"      "PRND"     "GPER"     "NT"       "PTCRA"    "PMPCA"    "AK5"      "GNG3"     "MARCKSL1"
##  [21] "GABBR2"   "ATP2C1"   "TRHR"     "NAT1"     "RGS19"    "GFP"      "ADRA1A"   "F2R"      "C5"       "CCR5"    
##  [31] "CXCR3"    "GNAI1"    "MCHR2"    "SFTPA1B"  "PIK3CG"   "CXCR4"    "PPP1CA"   "SH3BP5"   "EDN1"     "PTPN11"  
##  [41] "AMH"      "BCL2L1"   "CAMP"     "CD72"     "P2RX3"    "CCL27"    "KITLG"    "CD24"     "WNT7A"    "GJA1"    
##  [51] "LIF"      "CCL3L1"   "GP1BA"    "CLSTN1"   "IL5"      "CSF1"     "AMPD1"    "TNFSF14"  "LHB"      "NCR3"    
##  [61] "CCL14"    "LRPAP1"   "FGF2"     "MCL1"     "FGG"      "HNF4A"    "CD19"     "CORT"     "GRIK4"    "PLAUR"   
##  [71] "IL10"     "ITGAL"    "NMU"      "NTF5"     "TLR2"     "TRPM3"    "CACNB4"   "RGS17"    "ADM"      "CCL21"   
##  [81] "CD44"     "HBEGF"    "GRIN2D"   "TRPV4"    "TTR"      "NOS3"     "IGSF2"    "GNB5"     "PRLH"     "TREML1"  
##  [91] "PKD1"     "PLA2G7"   "DEFA3"    "NPSR1"    "PAPOLA"   "ADCY6"    "CACNB1"   "KLRC2"    "TP53"     "ATP5F1"  
## [101] "EIF2S1"   "SYPL2"    "ATP5G2"   "BCL2"     "HMGCS2"   "PPAP2A"   "GNB2"     "NT5C2"    "OTOP1"    "TRDN"    
## [111] "PILRA"    "CD22"     "ATP6V0E2" "KLK2"     "TYRP1"    "ATP6V0C"  "ACSS2"    "ATP6V1C1" "CNP"      "SLN"     
## [121] "PLCZ1"    "HOXA3"    "KL"       "ATP5L"    "FREQ"     "DEF6"     "ATP6V1E2" "ATP6V0D2" "CHGA"     "NUCB2"   
## [131] "ADSL"     "NT5C3"    "ATP6V1B1" "ATP6V0B"  "GPRIN1"   "MMP3"     "LPA"


a.hit.ca3.high<-apply(hits.high.rep.ca3,1,function(x) sum(x>1)>0)
a.hit.ca3.low<-apply(hits.low.rep.ca3,1,function(x) sum(x>1)>0)
a.hit.ca3<-a.hit.ca3.high | a.hit.ca3.low
all.hits.ca3<-unique(summary.ca3[a.hit.ca3,"Symbol"])

##   [1] "No Lipid" "C4A"      "CD5"      "PLK1"     "APOE"     "CAST"     "DRD2"     "ABCC3"    "BCL2L1"   "CCL26"   
##  [11] "ITGB3"    "VDR"      "CD8A"     "PMCH"     "CD24"     "CACNA1A"  "BDNF"     "GP1BA"    "RYR2"     "LILRA6"  
##  [21] "PTHLH"    "ALMS1"    "CCL1"     "RASGRP3"  "FCGR2B"   "ATP2C1"   "TRPC1"    "CACNA1F"  "LAIR1"    "ESR1"    
##  [31] "CACNA1H"  "NCR3"     "RELA"     "ATP5C1"   "CACNA1G"  "CD79B"    "CORT"     "IGF2"     "IL10"     "TRPM3"   
##  [41] "CACNB4"   "RGS17"    "TRPV4"    "TTR"      "EDN2"     "GRK4"     "NPR1"     "empty"    "IRS1"     "PRND"    
##  [51] "RAPGEF3"  "NPY2R"    "TACR1"    "TRAT1"    "GPR35"    "RAP2B"    "CAMK2G"   "CCL17"    "PTGIR"    "RARRES2" 
##  [61] "FGR"      "JAK2"     "CMKLR1"   "WFS1"     "AVPR1B"   "IL8RA"    "MARCKSL1" "CX3CR1"   "GHSR"     "GABBR2"  
##  [71] "TBXA2R"   "AKAP1"    "OXTR"     "TRHR"     "PTK2B"    "EDG7"     "CKM"      "CKMT1B"   "CSF1R"    "GNA14"   
##  [81] "SAA1"     "INSR"     "PIK3CG"   "ADCY1"    "ATP4B"    "ATP6V1G2" "CYP3A4"   "SRI"      "PRLH"     "SLC26A9" 
##  [91] "TREML1"   "GCM2"     "MYC"      "ATP5J"    "PLA2G5"   "ATP5A1"   "KIR2DS4"  "PPY"      "PKD1"     "WIPF1"   
## [101] "ATP5B"    "ATP6AP1"  "ABCC4"    "ADCY4"    "TMOD1"    "PAPOLA"   "DIO2"     "OPRS1"    "S100A1"   "ATP5F1"  
## [111] "EIF2S1"   "KLRC3"    "SERPINE2" "CACNB3"   "LAT2"     "BCL2"     "ITGAV"    "S100P"    "NUTF2"    "SYNJ1"   
## [121] "ATP5I"    "ASPH"     "FCGR2C"   "CD52"     "SLC34A1"  "CHERP"    "SLC9A3R2" "FOS"      "A4GALT"   "PLCG1"   
## [131] "ATP6V1C1" "CNP"      "F7"       "ATP2A3"   "PLCG2"    "ATP2B4"   "MMP1"     "ATP6V1B1" "ITPR2"    "REG3A"   
## [141] "MMP3"     "ATP6V1H"  "LPA"
all.hits<-c(all.hits,all.hits.ca3)
all.hits<-unique(all.hits)

##   [1] "No Cells" "OCRL"     "EDN2"     "ADRA2A"   "CCR7"     "CCR10"    "PLK1"     "GRK4"     "empty"    "EGFR"    
##  [11] "IRS1"     "SYK"      "PRND"     "GPER"     "NT"       "PTCRA"    "PMPCA"    "AK5"      "GNG3"     "MARCKSL1"
##  [21] "GABBR2"   "ATP2C1"   "TRHR"     "NAT1"     "RGS19"    "GFP"      "ADRA1A"   "F2R"      "C5"       "CCR5"    
##  [31] "CXCR3"    "GNAI1"    "MCHR2"    "SFTPA1B"  "PIK3CG"   "CXCR4"    "PPP1CA"   "SH3BP5"   "EDN1"     "PTPN11"  
##  [41] "AMH"      "BCL2L1"   "CAMP"     "CD72"     "P2RX3"    "CCL27"    "KITLG"    "CD24"     "WNT7A"    "GJA1"    
##  [51] "LIF"      "CCL3L1"   "GP1BA"    "CLSTN1"   "IL5"      "CSF1"     "AMPD1"    "TNFSF14"  "LHB"      "NCR3"    
##  [61] "CCL14"    "LRPAP1"   "FGF2"     "MCL1"     "FGG"      "HNF4A"    "CD19"     "CORT"     "GRIK4"    "PLAUR"   
##  [71] "IL10"     "ITGAL"    "NMU"      "NTF5"     "TLR2"     "TRPM3"    "CACNB4"   "RGS17"    "ADM"      "CCL21"   
##  [81] "CD44"     "HBEGF"    "GRIN2D"   "TRPV4"    "TTR"      "NOS3"     "IGSF2"    "GNB5"     "PRLH"     "TREML1"  
##  [91] "PKD1"     "PLA2G7"   "DEFA3"    "NPSR1"    "PAPOLA"   "ADCY6"    "CACNB1"   "KLRC2"    "TP53"     "ATP5F1"  
## [101] "EIF2S1"   "SYPL2"    "ATP5G2"   "BCL2"     "HMGCS2"   "PPAP2A"   "GNB2"     "NT5C2"    "OTOP1"    "TRDN"    
## [111] "PILRA"    "CD22"     "ATP6V0E2" "KLK2"     "TYRP1"    "ATP6V0C"  "ACSS2"    "ATP6V1C1" "CNP"      "SLN"     
## [121] "PLCZ1"    "HOXA3"    "KL"       "ATP5L"    "FREQ"     "DEF6"     "ATP6V1E2" "ATP6V0D2" "CHGA"     "NUCB2"   
## [131] "ADSL"     "NT5C3"    "ATP6V1B1" "ATP6V0B"  "GPRIN1"   "MMP3"     "LPA"

############# remove unique to one screen
posns<-match(all.hits,summary.ca2[,"Symbol"])
missing<-is.na(posns)
all.hits<-all.hits[!missing]

posns<-match(all.hits,summary.ca3[,"Symbol"])
missing<-is.na(posns)
all.hits<-all.hits[!missing]
##########################################

a.hit.gene.ca2<-summary.ca2[,"Symbol"] %in% all.hits
a.hit.gene.ca3<-summary.ca3[,"Symbol"] %in% all.hits


##########################choose one
a.screen<-a.screen.ca2
summary<-summary.ca2
hits.high.rep<-hits.high.rep.ca2
hits.low.rep<-hits.low.rep.ca2
a.hit.gene<-a.hit.gene.ca2
the.cells<-"MCF10a"
## a.hit.ca2.high<-a.hit.ca2.high
## a.hit.ca2.low

a.screen<-a.screen.ca3
summary<-summary.ca3
hits.high.rep<-hits.high.rep.ca3
hits.low.rep<-hits.low.rep.ca3
a.hit.gene<-a.hit.gene.ca3
the.cells<-"MDA-231"

###################

summary[summary[,"Symbol"] %in% all.hits[i],]
i<-i+1

sum.class<-{}
sum.median<-{}
num.reps<-{}
for(i in 1:length(all.hits)){
the.test<-summary[,"Symbol"] %in% all.hits[i]
num.reps<-c(num.reps,sum(the.test))
summary[summary[,"Symbol"] %in% all.hits[i],]

data<-summary[summary[,"Symbol"] %in% all.hits[i],c("per.R,AK,Rez,Num.cells",targets)]
class<-as.character(data[,"per.R,AK,Rez,Num.cells"])
data<-data[,targets]
dim(class)<-c(sum(the.test),1)
class
class<-apply(as.matrix(class),1,function(x) strsplit(x,split=", "))
class<-unlist(class)
dim(class)<-c(length(targets),sum(the.test))
class<-t(class)
class<-apply(class,2,function(x){ x<-gsub("^x$",0,x)
                                  x<-gsub("^-$",-1,x)
                                  x<-gsub("^\\+$",1,x)
                                as.numeric(x)
                                })

##### report hit means and regular mean of non-hits
## the.number<-apply(class,2,function(x) abs(sum(x)) )
the.median<-apply(abs(data*class),2,function(x) median(x[x!=0],na.rm=TRUE))
## the.mean<-the.sum/the.number
the.hits<-!is.na(the.median)
all.median<-apply(data,2,function(x) median(x,na.rm=TRUE))
the.median[!the.hits]<-all.median[!the.hits]


if(i==1){sum.class<-apply(class,2,function(x) sum(as.numeric(x)))}else{
sum.class<-rbind(sum.class,apply(class,2,function(x) sum(as.numeric(x)))  )
}


if(i==1){sum.median<-the.median}else{sum.median<-rbind(sum.median,the.median)}  


}
rownames(sum.class)<-all.hits
sum.class<-cbind(the.cells,all.hits,sum.class,num.reps,sum.median)
colnames(sum.class)<-c("Cell Line","Gene",targets,"num reps",paste("<",targets,">",sep=""))

sum.class[1:15,]

sum.class.ca3<-sum.class
sum.class.ca2<-sum.class

sum.class.ca3[1:15,]
sum.class.ca2[1:15,]

sum.class.ca2[,"Gene"]
sum.class.ca3[,"Gene"]

############combine so are alternating
dim(sum.class.ca2)
dim(sum.class.ca3)

both.screens<-rbind(sum.class.ca2,sum.class.ca3)
the.order<-as.numeric(unlist(apply(as.matrix(unique(both.screens[,"Gene"])),1,function(x) grep(paste("^",x,"$",sep=""),both.screens[,"Gene"]))))

both.screens<-both.screens[the.order,]
both.screens[1:6,]

getwd()
write.table(both.screens,"Ca2 and Ca3 SUMMARY3.txt",row.names=FALSE,sep="\t")
#save.image("working_final3.RData")
setwd( "/media/Bioinform-D/Research/Cellomics/Ca screen/Ca screen 2")
load("working_final2.RData")


AA1 , SLC26A9,SLN found
ADCY6

my.hits<-c("EDN2","NPR1","TACR1","ATP2C1","BDKRB2","CD8A","CD24","CCL3L1","GP1BA","CSF1","TRPC1","FGA","TRPM3","TRPV4","GNB5","TREML1","DEFA1","KIR3DL1","DEFA5","EIF2S1","SERPINE2","S100P","NUTF2","F7","ADCY6","SAA1","SLC26A9","SLN") 
my.2nd.hits<-c("CCR7","RAPGEF3","ADRA1B","IL10","NTF5","THY1","PAPOLA","DIO2","ATP5I","FOS","HSPA5")
conf.hits<-c("PLK1","BCL2L1","IRS1","EGFR","MCL1")
########################Clustering of genes

 ######Heatmap
 ####advanced heatmap
 library("gplots")


par.opts.start<-par(no.readonly=TRUE)$mar
[1] 5.1 4.1 4.1 2.1
par(mar=c(5.1, 4.1,4.1,2.1),mgp=c(3,1,0),las=2)  #c(bottom, left, top, right)
 par(font.axis=2)
 par(font=2)

both.screens[1:6,]
labels<-paste(both.screens[,"Cell Line"],both.screens[,"Gene"],sep="-")
rownames(both.screens)<-labels
 "Cell Line"    "Gene"         "percent_R"    "ak_leakage"   "resazurin"    "Num-cells"    "num reps"     "<percent_R>"  "<ak_leakage>""<resazurin>"  "<Num-cells>"

data2<-sum.class.ca2[,targets]
data2[1:5,]
the.dim2<-dim(data2)
data2<-as.numeric(data2)
dim(data2)<-the.dim2
rownames(data2)<-rownames(sum.class.ca2)
colnames(data2)<-targets
data2[1:5,]
data2<-data2/as.numeric(sum.class.ca2[,"num reps"])


data3<-sum.class.ca3[,targets]
data3[1:5,]
the.dim3<-dim(data3)
data3<-as.numeric(data3)
dim(data3)<-the.dim3
rownames(data3)<-rownames(sum.class.ca3)
colnames(data3)<-targets
data3[1:5,]
data3<-data3/as.numeric(sum.class.ca3[,"num reps"])


data.all<-both.screens[,targets]
data.all[1:5,]
the.dim.all<-dim(data.all)
data.all<-as.numeric(data.all)
dim(data.all)<-the.dim.all
rownames(data.all)<-rownames(both.screens)
colnames(data.all)<-targets
data.all<-data.all/as.numeric(both.screens[,"num reps"])

########################## choose one
data<-data.all

data<-data3-data2

data<-abs(data3-data2)
data[1:5,]


AVPR1B
GJB1     -0.50000000  0.33333333  0.000000000 -0.08333333
AVPR1B   -0.83333333  0.00000000  0.000000000  0.12500000
CD8A     -0.75000000   
ColSideColors=rep("black",times=dim(data)[1])

ColSideColors[rownames(data) %in% conf.hits]<-"green"
## ColSideColors[rownames(data) %in% "PLK1"]<-"green"
ColSideColors[rownames(data) %in% my.2nd.hits]<-"orange"
ColSideColors[rownames(data) %in% my.hits]<-"red"
rownames(data)[rownames(data) %in% my.hits]
               
data[rownames(data) %in% conf.hits,]
data2[rownames(data2) %in% conf.hits,]
data3[rownames(data3) %in% conf.hits,]

## colorss.ori<-colorss
## #colorss<-bluered(max(data)*4)
## color = colorRampPalette(c("blue","red"))
## colorss<-color(max(data)*4) #  color = colorRampPalette(c("blue","red"))
## colorss<-colorss[c(-2,-3,-4,-5,-6)]

the.heat<-heatmap.2(data,margins=c(8,8))
the.heat<-heatmap.2(t(data),margins=c(5,9),trace="row",ColSideColors=ColSideColors,keysize=0.8,cex.lab=2.5)
heatmap.2(as.matrix(data),Colv=patients.reorder,col=colorss,dendrogram="none",key=TRUE,keysize=1.47,scale="none", symkey=FALSE, density.info="none", trace="none",Rowv=FALSE,main="Citrullinated Antigen",ylab="Cytokine",margins=c(7,5),xlab="Patient :: CCP Status :: Number of Alleles")

savePlot("neive clustering.jpeg",type="jpeg")
savePlot("distance clustering.jpeg",type="jpeg")
savePlot("distance clustering.tiff",type="tiff")

SAA1 , SLC26A9,SLN found
ADCY6


the.pca <- prcomp(data,scale = TRUE) # for hits/genes
attributes(the.pca )
dim(the.pca$x)

the.pca.var <- round(the.pca$sdev^2 / sum(the.pca$sdev^2)*100,2)
plot(c(1:length(the.pca.var)),the.pca.var,type="b",xlab="# components",ylab="% variance",main="Scree Plot for Hits",col="red",cex=1.5,cex.lab=1.5)
savePlot("scree plot.jpeg",type="jpeg")

SLCcenters<-25
the.cl<-kmeans(the.pca$x[,1:2],centers=centers,iter.max=500) #Do kmeans

colours <-  rainbow(centers)

plot(range(the.pca$x[,1]),range(the.pca$x[,2]),xlab="PCA1",ylab="PCA2",main="Spectral clustering of differential hits")
text(the.pca$x[,1],the.pca$x[,2],label=rownames(the.pca$x),col=colours[the.cl$cluster],cex=0.75)




dotchart(the.pca$x[,1:3],labels=as.character(rownames(the.pca$x)))
################### choose one
the.test<-my.hits
color<-"red"

the.test<-my.2nd.hits
color<-"orange"

the.test<-conf.hits
color<-"green"
###########################

wanted<-rownames(the.pca$x) %in% the.test


points(the.pca$x[wanted,1],the.pca$x[wanted,2],col=color,cex=5.0)


savePlot("spectral.jpeg",type="jpeg")
savePlot("spectral.tiff",type="tiff")


library(scatterplot3d)



## scatterplot3d(the.pca$x[,1],the.pca$x[,2],the.pca$x[,3],angle=100,pch=as.character(rownames(the.pca$x)))

## scatterplot3d(the.pca$x[,1],the.pca$x[,2],the.pca$x[,3],color=colours[the.cl$cluster],pch=as.character(rownames(the.pca$x)))
s3d<-scatterplot3d(range(the.pca$x[,1]),range(the.pca$x[,2]),range(the.pca$x[,3]),xlab="PCA1",ylab="PCA2",zlab="PCA3",main="Spectral clustering of differential hits",angle=120)
## text(s3d$xyz.convert(range(the.pca$x[wanted,1]),range(the.pca$x[wanted,2]),range(the.pca$x[wanted,3])),col="red",cex=5.0)
text(s3d$xyz.convert(the.pca$x[,1],the.pca$x[,2],the.pca$x[,3]),label=rownames(the.pca$x),col=colours[the.cl$cluster],cex=0.75)
## points(the.pca$x[wanted,1],the.pca$x[wanted,2],col=color,cex=5.0)
points(s3d$xyz.convert(the.pca$x[wanted,1],the.pca$x[wanted,2],the.pca$x[wanted,3]),col=color,cex=5.0)


savePlot("3d clustering.jpeg",type="jpeg")
savePlot("3d clustering.tiff",type="tiff")


######## cleck a cluster with a knoWn member
the.member<-"TREML1"
the.member<-"SAA1"
group<-the.cl$cluster[the.cl$cluster==the.cl$cluster[names(the.cl$cluster)==the.member]]
group
the.pca$x[names(group),]

###################################

the.test<-my.hits
color<-"red"

the.test<-my.2nd.hits
color<-"orange"

the.test<-conf.hits
color<-"green"

wanted.hits<-rownames(the.pca$x) %in% the.test
wanted.2nd.hits<-rownames(the.pca$x) %in% the.test
wanted.control<-rownames(the.pca$x) %in% the.test


the.spin<-seq(from=0,to=360,by=10)
 Sys.sleep(5)
for (i in 1:length(the.spin)){
  print(i)
s3d<-scatterplot3d(range(the.pca$x[,1]),range(the.pca$x[,2]),range(the.pca$x[,3]),xlab="PCA1",ylab="PCA2",zlab="PCA3",main="Spectral clustering of differential hits",angle=the.spin[i])
## text(s3d$xyz.convert(range(the.pca$x[wanted,1]),range(the.pca$x[wanted,2]),range(the.pca$x[wanted,3])),col="red",cex=5.0)
text(s3d$xyz.convert(the.pca$x[,1],the.pca$x[,2],the.pca$x[,3]),label=rownames(the.pca$x),col=colours[the.cl$cluster],cex=0.75)
## points(the.pca$x[wanted,1],the.pca$x[wanted,2],col=color,cex=5.0)
points(s3d$xyz.convert(the.pca$x[wanted.hits,1],the.pca$x[wanted.hits,2],the.pca$x[wanted.hits,3]),col="red",cex=5.0)
points(s3d$xyz.convert(the.pca$x[wanted.2nd.hits,1],the.pca$x[wanted.2nd.hits,2],the.pca$x[wanted.2nd.hits,3]),col="orange",cex=5.0)
points(s3d$xyz.convert(the.pca$x[wanted.control,1],the.pca$x[wanted.control,2],the.pca$x[wanted.control,3]),col="green",cex=5.0)
Sys.sleep(1)
  savePlot(paste("3d clustering angle=",the.spin[i],".jpeg",sep=""),type="jpeg")
}

savePlot("3d clustering.jpeg",type="jpeg")
savePlot("3d clustering.tiff",type="tiff")



hits<-rep("",times=dim(the.pca$x)[1])
names(hits)<-rownames(the.pca$x)

hits[wanted.hits]<-"HIT"
hits[wanted.2nd.hits]<-"Possible HIT"
hits[wanted.control]<-"Control"

both.screens[1:6,]

ann[1:5,]
unique.genes<-unique(ann[,"Symbol"])
posns<-match(unique.genes,ann[,"Symbol"])
sum(is.na(posns))
use.ann<-ann[posns, c("Symbol","Entrez_ID","Cell_loc","Description","Gene.ID","Accession.Number","GI.Number","Library.Type")]
rownames(use.ann)<-use.ann[,"Symbol"]
use.ann[1:5,]

both.summary<-cbind(both.screens,the.cl$cluster[both.screens[,"Gene"]],the.pca$x[both.screens[,"Gene"],],hits[both.screens[,"Gene"]],use.ann[both.screens[,"Gene"],])
colnames(both.summary)[c(12,17)]<-c("cluster ID","Hit Class")
both.summary[1:15,]

getwd()
write.table(both.summary,"Ca2 and Ca3 hit SUMMARY.txt",row.names=FALSE,sep="\t")
#save.image("working_final2.RData")
setwd( "/media/Bioinform-D/Research/Cellomics/Ca screen/Ca screen 2")
load("working_final2.RData")








library(scatterplot3d)


 x <- seq(range(the.pca$x[,1])[1], range(the.pca$x[,1])[2], length= 100)
y <- seq(range(the.pca$x[,2])[1], range(the.pca$x[,2])[2], length= 100)
z <- seq(range(the.pca$x[,3])[1], range(the.pca$x[,3])[2], length= 100)
z <- outer(x, y, z)
z<-the.pca$x[,1:3]
persp(z)












    dotchart(VADeaths, main = "Death Rates in Virginia - 1940")
     op <- par(xaxs="i")# 0 -- 100%
     dotchart(t(VADeaths), xlim = c(0,100),
              main = "Death Rates in Virginia - 1940")
     par(op)
     

    require(grDevices) # for trans3d
     ## More examples in  demo(persp) !!
     ##                   -----------
     
     # (1) The Obligatory Mathematical surface.
     #     Rotated sinc function.
     
     x <- seq(-10, 10, length= 30)
     y <- x
     f <- function(x,y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }
     z <- outer(x, y, f)
     z[is.na(z)] <- 1
     op <- par(bg = "white")
     persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue")
     persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
           ltheta = 120, shade = 0.75, ticktype = "detailed",
           xlab = "X", ylab = "Y", zlab = "Sinc( r )"
     ) -> res
     round(res, 3)
     
 xE <- c(-10,10); xy <- expand.grid(xE, xE)
     points(trans3d(xy[,1], xy[,2], 6, pmat = res), col = 2, pch =16)
     lines (trans3d(x, y=10, z= 6 + sin(x), pmat = res), col = 3)
     
     phi <- seq(0, 2*pi, len = 201)
     r1 <- 7.725 # radius of 2nd maximum
     xr <- r1 * cos(phi)
     yr <- r1 * sin(phi)
     lines(trans3d(xr,yr, f(xr,yr), res), col = "pink", lwd = 2)
     ## (no hidden lines)
     
     # (3) Visualizing a simple DEM model
     
     z <- 2 * volcano        # Exaggerate the relief
     x <- 10 * (1:nrow(z))   # 10 meter spacing (S to N)
     y <- 10 * (1:ncol(z))   # 10 meter spacing (E to W)
     ## Don't draw the grid lines :  border = NA
     par(bg = "slategray")
     persp(x, y, z, theta = 135, phi = 30, col = "green3", scale = FALSE,
           ltheta = -120, shade = 0.75, border = NA, box = FALSE)
     
     # (4) Surface colours corresponding to z-values
     
     par(bg = "white")
     x <- seq(-1.95, 1.95, length = 30)
     y <- seq(-1.95, 1.95, length = 35)
     z <- outer(x, y, function(a,b) a*b^2)
     nrz <- nrow(z)
     ncz <- ncol(z)
     # Create a function interpolating colors in the range of specified colors
     jet.colors <- colorRampPalette( c("blue", "green") ) 
     # Generate the desired number of colors from this palette
     nbcol <- 100
     color <- jet.colors(nbcol)
     # Compute the z-value at the facet centres
     zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
     # Recode facet z-values into color indices
     facetcol <- cut(zfacet, nbcol)
     persp(x, y, z, col=color[facetcol], phi=30, theta=-30)
     
     par(op)
     


 require(grDevices) # for colours
     x <- y <- seq(-4*pi, 4*pi, len=27)
     r <- sqrt(outer(x^2, y^2, "+"))
     image(z = z <- cos(r^2)*exp(-r/6), col=gray((0:32)/32))
     image(z, axes = FALSE, main = "Math can be beautiful ...",
           xlab = expression(cos(r^2) * e^{-r/6}))
     contour(z, add = TRUE, drawlabels = FALSE)
     
     # Volcano data visualized as matrix. Need to transpose and flip
     # matrix horizontally.
     image(t(volcano)[ncol(volcano):1,])
     
     # A prettier display of the volcano
     x <- 10*(1:nrow(volcano))
     y <- 10*(1:ncol(volcano))
     image(x, y, volcano, col = terrain.colors(100), axes = FALSE)
     contour(x, y, volcano, levels = seq(90, 200, by = 5),
             add = TRUE, col = "peru")
     axis(1, at = seq(100, 800, by = 100))
     axis(2, at = seq(100, 600, by = 100))
     box()
     title(main = "Maunga Whau Volcano", font.main = 4)
     





#points(sam_genes.pca$x[,1],sam_genes.pca$x[,2],col=colours[sam_genes.cl$cluster],pch=sam_genes.cl$cluster,cex=1.0,bg=colours[sam_genes.cl$cluster]) #colours and symbols
text(sam_genes.pca$x[,1],sam_genes.pca$x[,2],label=rownames(sam_genes.pca$x),col=colours[sam_genes.cl$cluster],cex=0.75) #colours and symbols

rownames(sam_genes.pca$x)<- diff_expressed_sam_SYM





 diff_expressed_sam_SYM<-mget(rownames(diff_expressed_sam),env=rgu34aSYMBOL)
 diff_expressed_sam_GENE<-mget(rownames(diff_expressed_sam),env=rgu34aGENENAME)

# if no symbol use probe name
diff_expressed_sam_GENE[diff_expressed_sam_GENE=="NA"]<-c(labels(diff_expressed_sam_GENE)[diff_expressed_sam_GENE=="NA"])

diff_expressed_sam_SYM[diff_expressed_sam_SYM=="NA"]<-c(labels(diff_expressed_sam_SYM)[diff_expressed_sam_SYM=="NA"])
diff_expressed_sam_SYM[diff_expressed_sam_SYM=="na"]<-c(labels(diff_expressed_sam_SYM)[diff_expressed_sam_SYM=="na"])

 #rownames(diff_expressed_sam)<-diff_expressed_sam_SYM
 
 ###To undo
# rownames(diff_expressed_sam)<-labels(diff_expressed_sam_SYM)

## truncate long names for heatmap
diff_expressed_sam_GENEt<- substr(as.character(diff_expressed_sam_GENE),1,45)
rownames(diff_expressed_sam)<-diff_expressed_sam_GENEt
 colnames(diff_expressed_sam)<-samples
 ##heatmap with gene names
 heatmap(as.matrix(diff_expressed_sam),col=cm.colors(16),RowSideColors=colours[sam_genes.cl$cluster],main="Heat map: Class and Genes (SAM)",xlab="Class")

 ##heatmap with gene symbols
 rownames(diff_expressed_sam)<-diff_expressed_sam_SYM
 heatmap(as.matrix(diff_expressed_sam),col=cm.colors(16),main="Heat map: Class and Genes (SAM)",xlab="Class")
  ### purple is higher expressed


  heatmap(as.matrix(diff_expressed_sam),col=cm.colors(16), col=colours[sam_genes.cl$cluster],main="Heat map: Class and Genes (SAM)",xlab="Class")

 ##############Spectral decomposition

 ##PCA first
 #PCS on genes    SCALE set since expression values differ by 3-4 orders of magnitude

#### QUESTION sam_class.pca <- prcomp(diff_expressed_sam,scale = TRUE)
sam_genes.pca <- prcomp(diff_expressed_sam,scale = TRUE)
> attributes(sam_genes.pca )
$names
[1] "sdev"     "rotation" "center"   "scale"    "x"

$class
[1] "prcomp"

dim(sam_genes.pca$x)
[1] 45 15

 dim(sam_genes.pca$x)  #with transpose
[1] 15 15


####scree plot

sam_genes.pca.var <- round(sam_genes.pca$sdev^2 / sum(sam_genes.pca$sdev^2)*100,2)
plot(c(1:length(sam_genes.pca.var)),sam_genes.pca.var,type="b",xlab="# components",ylab="% variance",main="Scree plot-Genes (SAM)",col="red")

 ###first 2 component have most or varaiblity so use K-means on frist two components
 ##from Heat map 7-8 clusters available

sam_genes.cl<-kmeans(sam_genes.pca$x[,1:2],centers=7,iter.max=500) #Do kmeans



################ Heat plot used :
plot(range(sam_genes.pca$x[,1]),range(sam_genes.pca$x[,2]),xlab="PCA1",ylab="PCA2",main="Spectral clustering Genes (SAM)")
colours <-  rainbow(7)
#points(sam_genes.pca$x[,1],sam_genes.pca$x[,2],col=colours[sam_genes.cl$cluster],pch=sam_genes.cl$cluster,cex=1.0,bg=colours[sam_genes.cl$cluster]) #colours and symbols
text(sam_genes.pca$x[,1],sam_genes.pca$x[,2],label=rownames(sam_genes.pca$x),col=colours[sam_genes.cl$cluster],cex=0.75) #colours and symbols

 rownames(sam_genes.pca$x)<- diff_expressed_sam_SYM
 plot(range(sam_genes.pca$x[,1]),range(sam_genes.pca$x[,2]),xlab="PCA1",ylab="PCA2",main="Spectral clustering Genes (SAM)")
text(sam_genes.pca$x[,1],sam_genes.pca$x[,2],label=rownames(sam_genes.pca$x),col=colours[sam_genes.cl$cluster],cex=0.75)



rownames(diff_expressed_sam)<-labels(diff_expressed_sam_SYM)
 diff_expressed_sam_ID<-mget(rownames(diff_expressed_sam),env=rgu34aLOCUSID)

rownames(diff_expressed_korn)<-labels(diff_expressed_korn_SYM)
 diff_expressed_korn_ID<-mget(rownames(diff_expressed_korn),env=rgu34aLOCUSID)



