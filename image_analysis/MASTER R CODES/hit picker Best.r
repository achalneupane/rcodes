


############################################### Jane VALIDATION  ### latest DNA and flexible annotation
## setwd("/media/Bioinform-D/Research/Cellomics/Jane/Validation/Green 128 run")
setwd("/media/Bioinform-D/Research/Cellomics/Jane/Validation") ## this is thegreen 32 run
load("Jane_validation_ann.RData")  ### "plate" must match "plate" in file  cellomics to expt map.csv
core.ann<-c("plate", "row", "col", "Symbol")
place.core.ann<-colnames(ann) %in% core.ann

core.vars<-c("%R","Trans_%","Prolif_%","Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","RG","RG_expected","RG_sd","RG_Low","RG_Low_expected","RG_Mid","RG_Mid_expected","RG_High","RG_High_expected","R","RnotG","G","nG","G_Low","G_Mid","G_High","Relative S+G2","Relative S+G2+>4N","Relative <2N","Relative >4N","Relative S","G2/G1 for fit","G1 posn", "Fit quality G","Fit quality nG","Num-cells","score","RG/G","RnG/nG","nR","gt4N","Total","<cells_field>","sd_cells_field","S+G2/G1","S+G2+ >4N/G1","<2N",">4N","G1","G2","S","Gr-S+G2/G1","Gr-S+G2+ >4N/G1","Gr-<2N","Gr->4N","Gr-G1","Gr-G2","Gr-S","NtGr-S+G2/G1","NtGr-S+G2+ >4N/G1","NtGr-<2N","NtGr->4N","NtGr-G1","NtGr-G2","NtGr-S"
             )

the.screen<-"Jane_validation"
files<- paste("plate_LIZ_HITS_",unique(ann[,"plate"]),"_SUMMARY",".RData",sep="")
field.output.file<-"Jane_validation_field_summary_NOTGREEN_DNA.txt"
well.output.file<-"Jane_validation_summary_NOTGREEN_DNA.txt"
well.type<-96
row.type<-8
col.type<-12

normalized.file<-paste(the.screen,"NORMALIZED","txt",sep=".")

targets<-c("Num-cells","Z-all Green","RG_over_G","RG_over_G_Z","Gr-G1","RG_over_G_over_EtoH")
validation<-TRUE
load("a.screen.cuts.RData")
load("a.screen.stats.RData")
load("a.screen.counts.RData")
num.reps.to.pass<-1  ### to initially be recognised
percent.reps.to.be.hit<-0.4  ### 0.4 GREATER THAN 40 of the clones must be a hit : reject: 1/3, 2/5  kepp  2/3 , 2/4 ,3/5

pos.control.list<-c("CCND1","C-MYC")
norm.control.list<-c("plv411") # controls that could be used for  normalization
confounders.list<-c("UTMOCK","MOCK")
#######################################################


############################################BEGIN FIRST need to process all the independent screen and collate the hits for acah screen
## ###for single color siRNA
## pos.control.list<-c("PLK1")
## norm.control.list<-c("NT","LAMA","GFP","ATP2C1","No Lipid") # controls that could be used for  normalization

#### new annotation way
the.screen.name<-the.screen
normalized<-read.delim(normalized.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
col<-colnames(normalized) %in% core.ann
labels<-c(colnames(normalized)[col],core.vars,colnames(ann)[!place.core.ann])## 
colnames(normalized)[1:length(labels)]<-labels

plates.all<-normalized
labels<-gsub("%","percent_",colnames(plates.all))
labels<-gsub("/","_over_",labels)
labels<-gsub(">","_gt_",labels)
colnames(plates.all)<-labels
colnames(plates.all)<-gsub("/","_",colnames(plates.all))

the.screen


toString(colnames(plates.all))

extras.for.indiv.filtering<-c("RG", "RG_expected", "RG_sd", "Total","G")
keep<-c(colnames(ann),targets,extras.for.indiv.filtering)

## keep<-c("plate","row","col","Symbol","Num-cells","Z-all Green",  "Relative S+G2+_gt_4N","Gr-S+G2+ _gt_4N_over_G1","RG_over_G_Z","Gr-G1","RG_over_G","RG_over_G_over_EtoH","Z-Low Green", "Z-Mid Green", "Z-High Green", "RG", "RG_expected", "RG_sd", "Total","G","core.plate", "Accessions", "GeneID", "MusIDT", "M.position", "TREAT" ) ##for Jane

## keep<-c("plate","row","col","Symbol","percent_R","ak_leakage","resazurin","Num-cells","S+G2_over_G1","S+G2+ >4N_over_G1","<2N",">4N","MEAN_Area","SD_Area","MEAN_P2A","SD_P2A","MEAN_LWR","SD_LWR") ### for sica2 and 3



## targets<-c("percent_R", "ak_leakage","resazurin","Num-cells"  ) # defined above

################################################### if norm not a file for one cell line then breakup now!
plates.all.ori<-plates.all




########################RESTART to do individual screens
plates.all.ori[1:5,]
##  norm<-norm.all

a.screen.cuts


the.screens<-names(a.screen.cuts)
the.screens

##### sometimes you might not want to run all the screens at once:
## to.exclude<-c("4a","4b","4c") # to.exclude<-c("a","b","c")
to.exclude<-c()

run.screens<-the.screens[!(the.screens %in% to.exclude)]
run.screens
targets

## recount<-FALSE ### use TRUE if have existing thresholds:  reget the statistics but do no mofify thresholds
## cut.sd.posn<-1.64 ## 1.64 for  p=0.5 for upper limit  ## cut.sd.posn<-2 for p=0.03 in upper tail
## reps.for.hit<-2 ## reps.for.hit<-


for(k in 1:length(run.screens)){

############### to do for each screen
the.screen<-run.screens[k]

wanted.screen<-grepl(the.screen,plates.all[,"plate"])
and.exclude<-!(plates.all[,"plate"] %in% to.exclude)

sum(wanted.screen & and.exclude)

#### one off
if(the.screen.name=="Jane_validation"){
if(the.screen=="a" | the.screen=="b" | the.screen=="c" ) {and.exclude<-!(plates.all[,"plate"] %in% c("4a","4b","4c"))}
if(the.screen=="4a" | the.screen=="4b" | the.screen=="4c" ) {and.exclude<-!(plates.all[,"plate"] %in% c("1a","1b","1c","2a","2b","2c","3a","3b","3c"))}
}
###

################## norm holds the screen data
norm<-plates.all.ori[wanted.screen & and.exclude,keep]

print(unique(norm[,"plate"]))

##################################
## hits.high and hits.low contain the position on individual high and low hits in boolean format
## hiths.high.all and hits.low a.ll contain the position of calif wells
## hits.high.rep/hights.low.rep  label the total counts high or low that are hits, bit that are provided at verey clone postion (even those that are not hits;

low.cut<-a.screen.cuts[[the.screen]]$low.cut
high.cut<-a.screen.cuts[[the.screen]]$high.cut

hits.low<-matrix(FALSE,nrow=dim(norm)[1],ncol=length(targets))
hits.high<-matrix(FALSE,nrow=dim(norm)[1],ncol=length(targets))
hits.low.all<-matrix(FALSE,nrow=dim(norm)[1],ncol=length(targets))
hits.high.all<-matrix(FALSE,nrow=dim(norm)[1],ncol=length(targets))
colnames(hits.low)<-targets
colnames(hits.high)<-targets
colnames(hits.low.all)<-targets
colnames(hits.high.all)<-targets

a.screen<-{}
hits.low.rep<-matrix(0,nrow=dim(norm)[1],ncol=length(targets))
hits.high.rep<-matrix(0,nrow=dim(norm)[1],ncol=length(targets))
hits.full.rep<-matrix(0,nrow=dim(norm)[1],ncol=length(targets))
colnames(hits.low.rep)<-targets
colnames(hits.high.rep)<-targets
colnames(hits.full.rep)<-targets

###############

  
  
for(i in 1:length(targets)){

the.target<-targets[i]

if(the.screen=="a" & the.target=="RG_over_G_over_EtoH"){next} ### cause it's all 1 by defination
if(the.screen=="4a" & the.target=="RG_over_G_over_EtoH"){next}

print(the.screen)



print(the.target)
norm[1:5,the.target]

norm[1:5,]
####sica2 (ann in one)
dim(norm)
low.cut
high.cut


filter.low<-!is.na(norm[,the.target])
filter.high<-filter.low  ### both shoudl not contain and NA for that attribute

########### Note low cutoff migh be different to high cutoff
### low rations ok with low numbers but high ratios with low numbers are probably artifacts
if(grepl("percent_R",the.target)){filter.high<- filter.low & norm[,"R"]>50} ### stop high percentages due so small counts (high num cells and low red ok)
if(grepl("RG_over_G",the.target)){filter.high<- filter.low & norm[,"G"]>=50} ### low ratio ok with small numbers of green cells
if(grepl("RG_over_G_Z",the.target)){filter.high<- filter.low & norm[,"G"]>=50}
if(grepl("Z-all Green",the.target)){filter.high<- filter.low & norm[,"RG"]>=35} ## have used RG-expected also
if(grepl("Num-cells",the.target)){filter.high<- filter.low & norm[,"Total"]>=100} ## true numer of cells must be large enough
if(grepl("Gr-G1",the.target)){filter.high<- filter.low & norm[,"G"]>=50}
if(grepl("RG_over_G_over_EtoH",the.target)){filter.high<- filter.low & norm[,"G"]>=50}
###################################################
filter.low<-filter.high # for jane validation
  
hits.low[,the.target]<-(as.numeric(norm[,the.target]) <= low.cut[the.target] ) & filter.low

hits.high[,the.target]<-(as.numeric(norm[,the.target]) >= high.cut[the.target] ) & filter.high

hits.high.all[,the.target]<- filter.high

hits.low.all[,the.target]<- filter.low

print(paste("low cut",the.target,sum(hits.low[,the.target]),sep=" : "))
print(paste("high cut",the.target,sum(hits.high[,the.target]),sep=" : "))


###------------ move to above

###------------move to above
full.counts<-sort(tapply(norm[(filter.high & filter.low),"Symbol"],norm[(filter.high & filter.low),"Symbol"],length)) ## cchanges from target to target cause norm changes
## lost.counts<-sort(tapply(is.na(norm[,the.target]),norm[,"Symbol"],sum))
## lost.counts<-lost.counts[names(full.counts)]
## lost.counts[is.na(lost.counts)]<-0
## full.counts<-full.counts-lost.counts


### these are the correct counts filtered for NAs AND other filtering casue hits.high and hits.ow are specif to each well
low.counts<- sort(tapply(norm[hits.low[,the.target],"Symbol"],norm[hits.low[,the.target],"Symbol"],length))
high.counts<- sort(tapply(norm[hits.high[,the.target],"Symbol"],norm[hits.high[,the.target],"Symbol"],length))
 
## norm[hits.low[,the.target],]
## norm[hits.high[,the.target],]
 
if(length(low.counts)< 1){posns<-{}}else{
  ## posns<-apply(as.matrix(names(low.counts)),1,function(x) grep(paste("^",x,"$",sep=""),norm[,"Symbol"])) # warning will return a matrix if square
  posns<-lapply(as.list(names(low.counts)),function(x) grep(paste("^",x[[1]],"$",sep=""),norm[,"Symbol"]))
}
   for(ii in 1:length(posns)){
     x<-posns[[ii]]
     hits.low.rep[x,the.target]<-low.counts[ii] ## hits.low.rep[,the.target]
   }

 if(length(high.counts)< 1){posns<-{}}else{
   ## posns<-apply(as.matrix(names(high.counts)),1,function(x) grep(paste("^",x,"$",sep=""),norm[,"Symbol"])) # warning will return a matrix if square
   posns<-lapply(as.list(names(high.counts)),function(x) grep(paste("^",x[[1]],"$",sep=""),norm[,"Symbol"]))
 }
   for(ii in 1:length(posns)){
     x<-posns[[ii]]
     hits.high.rep[x,the.target]<-high.counts[ii]
   }

if(length(full.counts)< 1){posns<-{}}else{
  ## posns<-apply(as.matrix(names(low.counts)),1,function(x) grep(paste("^",x,"$",sep=""),norm[,"Symbol"])) # warning will return a matrix if square
  posns<-lapply(as.list(names(full.counts)),function(x) grep(paste("^",x[[1]],"$",sep=""),norm[,"Symbol"]))
}
   for(ii in 1:length(posns)){
     x<-posns[[ii]]
     hits.full.rep[x,the.target]<-full.counts[ii] ## hits.low.rep[,the.target]
   }

 
 ## all.counts<-full.counts[unique(c(names(high.counts),names(low.counts)))]
 a.test<-list(low.counts=low.counts,high.counts=high.counts,full.counts=full.counts)

#####$$$$$$$$$$$$$$$
 a.screen[i]<- list(a.test)

 
 print(paste("Rep hits low cut",the.target,sum(low.counts>1),sep=" : "))
 print(paste("Rep hits high cut",the.target,sum(high.counts>1),sep=" : "))


names(a.screen)[i]<-the.target
a.screen[i]

## hits.rep<-matrix(0,nrow=dim(norm)[1],ncol=length(targets))
## colnames(hits.rep)<-targets

###### fix cases wher genes in both extremes often a spacial defect majority rules

### 0...........0.33. 0.4.........0.6..0.66........1

 high.variation<-(hits.high.rep[,the.target] !=0 & hits.low.rep[,the.target] !=0)

## norm[high.variation,]
hits.high.rep[high.variation,the.target]
hits.low.rep[high.variation,the.target]
norm[high.variation,"Symbol"]


 
 print(paste("Number to resolve->",sum(hits.high.rep[,the.target] !=0 & hits.low.rep[,the.target] !=0),sep=""))
 ## even.call<-hits.high.rep[high.variation,the.target] == hits.low.rep[high.variation,the.target] # kill both

 sum.hits<-hits.high.rep[high.variation,the.target] + hits.low.rep[high.variation,the.target]
 even.call<-(hits.high.rep[high.variation,the.target]/sum.hits >= 0.4) & (hits.high.rep[high.variation,the.target]/sum.hits <= 0.6)

 norm[high.variation,][even.call,"Symbol"] #  norm[high.variation,][even.call,]
 even.genes<-unique(norm[high.variation,][even.call,"Symbol"])

 hits.high.rep[high.variation,][even.call,the.target]<-0
 hits.low.rep[high.variation,][even.call,the.target]<-0

 hits.high[high.variation,][even.call,the.target]<-FALSE   ## sum(hits.high[high.variation,][even.call,the.target]) # Some may already be false ok just kill all
 hits.low[high.variation,][even.call,the.target]<-FALSE   ## sum(hits.low[high.variation,][even.call,the.target])
 ## remove hits from a.screen  names(a.screen[[the.target]]$"low.counts")
 a.screen[[the.target]][["low.counts"]]<-a.screen[[the.target]][["low.counts"]][!(names(a.screen[[the.target]][["low.counts"]]) %in% even.genes)]
 a.screen[[the.target]][["high.counts"]]<-a.screen[[the.target]][["high.counts"]][!(names(a.screen[[the.target]][["high.counts"]]) %in% even.genes)] 
 


  high.call<-hits.high.rep[high.variation,the.target] > hits.low.rep[high.variation,the.target]
  high.genes<-unique(norm[high.variation,][high.call,"Symbol"])
  hits.low.rep[high.variation,][high.call,the.target]<-0 # is high ;  kill low # hits.high.rep[high.variation,][high.call,the.target]
  hits.low[high.variation,][high.call,the.target]<-FALSE 
  a.screen[[the.target]][["low.counts"]]<-a.screen[[the.target]][["low.counts"]][!(names(a.screen[[the.target]][["low.counts"]]) %in% high.genes)]

  low.call<-hits.high.rep[high.variation,the.target] < hits.low.rep[high.variation,the.target]
  low.genes<-unique(norm[high.variation,][low.call,"Symbol"])
  hits.high.rep[high.variation,][low.call,the.target]<-0 # is low ; kill high
  hits.high[high.variation,][low.call,the.target]<-FALSE  
  a.screen[[the.target]][["high.counts"]]<-a.screen[[the.target]][["high.counts"]][!(names(a.screen[[the.target]][["high.counts"]]) %in% low.genes)]
 
 print(paste("Done; this should be zero->",sum(hits.high.rep[,the.target] !=0 & hits.low.rep[,the.target] !=0),sep=""))


} # end loop over targets


lapply(a.screen,function(x) lapply(x,length))

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

colnames(summary)[dim(summary)[2]]<-toString(targets)
summary[1:5,]
dim(summary)
write.table(summary,paste(the.screen,"SUMMARY","txt",sep="."),row.names=FALSE,sep="\t")
save(list=c("the.screen","a.screen","summary","hits.high.rep","hits.low.rep","hits.full.rep","hits.low","hits.high","hits.low.all","hits.high.all"),file=paste("data.",the.screen,".RData",sep=""))

} ## end loop over screens


##------loop over screen


##############################################################################
#############################################################################

 data.names<-paste("data.",the.screens,".RData",sep="")


## data.names<-  c("data.MUS.a.RData","data.MUS.b.RData","data.MUS.c.RData")

 ## data.names<-c("data.Ca siRNA2.RData","data.Ca siRNA3.RData")





## ### collect all good hits
data.names
all.hits<-{}
summary.genes<-list()
 for(i in 1:length(data.names)){
   load(data.names[i])
   print(paste("In Directory: ",getwd(),sep=" "))
   print(paste("doing: ",data.names[i],sep=" "))
 hits.full.rep  
a.hit.high<-apply(hits.high.rep/hits.full.rep,1,function(x) sum( (x> percent.reps.to.be.hit) & is.finite(x)  )>0) # test<-hits.high.rep/hits.full.rep
a.hit.low<-apply(hits.low.rep/hits.full.rep,1,function(x)  sum( (x> percent.reps.to.be.hit) & is.finite(x)  )>0)
a.hit<-a.hit.high | a.hit.low
all.hits<-c(all.hits,summary[a.hit,"Symbol"])
all.hits<-unique(all.hits)
## data<-eval(as.name(data.names[i]))
summary.genes[[i]]<-summary[,"Symbol"] ##this is just all the genes in the screen
names(summary.genes)[i]<-the.screen

## assign(data.names[i],value=data)
}


############## remove form all hits genes that occur in only one screen:
length(all.hits)
for(i in 1:length(summary.genes)){
posns<-match(all.hits,summary.genes[[i]])
missing<-is.na(posns)
paste("Removing::")
print(all.hits[missing])
all.hits<-all.hits[!missing]
}
length(all.hits)


###################c("the.screen","a.screen","summary","hits.high.rep","hits.low.rep","hits.full.rep") the.test<-summary[,"Symbol"] %in% "DDX21"
## all.screens<-{}
## for(ii in 1:length(data.names)){
## load(data.names[ii])
   
## sum.class<-{}
## sum.median<-{}
## num.reps<-{}
## for(i in 1:length(all.hits)){
## the.test<-summary[,"Symbol"] %in% all.hits[i]
## hits.full.rep[the.test,]
## num.reps<-c(num.reps,sum(the.test))
## summary[the.test,]
## hits.low.rep[the.test,]
## hits.high.rep[the.test,]

## data<-summary[summary[,"Symbol"] %in% all.hits[i],c(toString(targets),targets)]
## ## data<-summary[summary[,"Symbol"] %in% all.hits[i],c("per.R,AK,Rez,Num.cells",targets)]
## class<-as.character(data[,toString(targets)])
## data<-data[,targets]
## dim(class)<-c(sum(the.test),1)
## class
## class<-apply(as.matrix(class),1,function(x) strsplit(x,split=", "))
## class<-unlist(class)
## dim(class)<-c(length(targets),sum(the.test))
## class<-t(class)
## class<-apply(class,2,function(x){ x<-gsub("^x$",0,x)
##                                   x<-gsub("^-$",-1,x)
##                                   x<-gsub("^\\+$",1,x)
##                                 as.numeric(x)
##                                 })

## ##### report hit means and regular mean of non-hits
## ## the.number<-apply(class,2,function(x) abs(sum(x)) )
## the.median<-apply(abs(data*class),2,function(x) median(x[x!=0],na.rm=TRUE))
## ## the.mean<-the.sum/the.number
## the.hits<-!is.na(the.median)
## all.median<-apply(data,2,function(x) median(x,na.rm=TRUE))
## the.median[!the.hits]<-all.median[!the.hits]

## dim(class)<-c(num.reps[i],length(targets))
## if(i==1){
##   sum.class<-apply(class,2,function(x) sum(as.numeric(x)))

##        }else{
## sum.class<-rbind(sum.class,apply(class,2,function(x) sum(as.numeric(x)))  )
## }


## if(i==1){sum.median<-the.median}else{sum.median<-rbind(sum.median,the.median)}  


## } ### end loop over all hits 
## rownames(sum.class)<-all.hits
## sum.class<-cbind(the.screen,all.hits,sum.class,num.reps,sum.median)
## colnames(sum.class)<-c("Cell Line","Gene",targets,"num reps",paste("<",targets,">",sep=""))


## all.screens<-rbind(all.screens,sum.class)

## } # loop over screens

################################################collate all screens hits and the median scores
all.screens<-{}

for(ii in 1:length(data.names)){
load(data.names[ii])

   
sum.class<-{}
for(i in 1:length(targets)){
  the.target<-targets[i]
high<-tapply(summary[hits.high[,the.target],the.target],summary[hits.high[,the.target],"Symbol"],median) ## summary[hits.high[,the.target],]
low<-tapply(summary[hits.low[,the.target],the.target],summary[hits.low[,the.target],"Symbol"],median)
all<-tapply(summary[(hits.high.all[,the.target] | hits.low.all[,the.target]) ,the.target],summary[(hits.high.all[,the.target] | hits.low.all[,the.target]) ,"Symbol"],median) ##hits.low.all just the ones with NAs this ok I think

the.median<-all
the.median[names(high)]<-high
the.median[names(low)]<-low
  
class.high<-tapply(summary[hits.high[,the.target],the.target],summary[hits.high[,the.target],"Symbol"],length)
class.low<--1*tapply(summary[hits.low[,the.target],the.target],summary[hits.low[,the.target],"Symbol"],length) ## -1 to denote is a negative hit
class.sum<-c(class.high,class.low) # the two groups have no common members
num.reps<-tapply(summary[(hits.high.all[,the.target] | hits.low.all[,the.target]),the.target],summary[(hits.high.all[,the.target] | hits.low.all[,the.target]),"Symbol"],length)
  
    if(i==1){
sum.class<-cbind(the.screen,all.hits,class.sum[all.hits],num.reps[all.hits],the.median[all.hits])
rownames(sum.class)<-all.hits
colnames(sum.class)<-c("Cell Line","Gene",the.target,paste("num reps",the.target,sep="."),paste("<",the.target,">",sep=""))
  
sum.class[is.na(sum.class[,the.target]),the.target]<-0
sum.class[is.na(sum.class[,the.target]),paste("<",the.target,">",sep="")]<-0
sum.class[is.na(sum.class[,paste("num reps",the.target,sep=".")]),paste("num reps",the.target,sep=".")]<-0

       }else{
sum.class<-cbind(sum.class,class.sum[all.hits],num.reps[all.hits],the.median[all.hits])
colnames(sum.class)<-c(colnames(sum.class)[1:(dim(sum.class)[2]-3)],the.target,paste("num reps",the.target,sep="."),paste("<",the.target,">",sep=""))
  
sum.class[is.na(sum.class[,the.target]),the.target]<-0
sum.class[is.na(sum.class[,paste("<",the.target,">",sep="")]),paste("<",the.target,">",sep="")]<-0
sum.class[is.na(sum.class[,paste("num reps",the.target,sep=".")]),paste("num reps",the.target,sep=".")]<-0
}

  
}

  all.screens<-rbind(all.screens,sum.class)

} # loop over screens

all.screens[,the.target]
all.screens[,paste("<",the.target,">",sep="")]


wanted.order<-c("Cell Line","Gene",targets,paste("num reps",targets,sep="."),paste("<",targets,">",sep=""))
all.screens<-all.screens[,wanted.order]

### note means of a different sign to before...



the.order<-as.numeric(unlist(apply(as.matrix(unique(all.screens[,"Gene"])),1,function(x) grep(paste("^",x,"$",sep=""),all.screens[,"Gene"]))))
all.screens<-all.screens[the.order,]
all.screens[1:6,]

## if(the.screen.name=="Jane_validation"){


## }


getwd()
## write.table(both.screens,"Ca2 and Ca3 SUMMARY2.txt",row.names=FALSE,sep="\t")

write.table(all.screens,"Combined_Hits.txt",row.names=FALSE,sep="\t")
save(list=c("all.screens","targets","a.screen.cuts","the.screen.name"),file=paste(the.screen.name,"Combined Hits.RData",sep=" "))



######################################## set up for PCA analyis 

my.hits<-c("EDN2","NPR1","TACR1","ATP2C1","BDKRB2","CD8A","CD24","CCL3L1","GP1BA","CSF1","TRPC1","FGA","TRPM3","TRPV4","GNB5","TREML1","DEFA1","KIR3DL1","DEFA5","EIF2S1","SERPINE2","S100P","NUTF2","F7","ADCY6","SAA1","SLC26A9","SLN") 
my.2nd.hits<-c("CCR7","RAPGEF3","ADRA1B","IL10","NTF5","THY1","PAPOLA","DIO2","ATP5I","FOS","HSPA5")
conf.hits<-c("PLK1","BCL2L1","IRS1","EGFR","MCL1")


########################Clustering of genes
my.hits<-c("TGFBR1","STK16","IRAK3","PRKCN","KSR2","PRPSAP2","MAGI-3","GRK4","AURKC","KALRN","DGKE","TPRXL","AURKB","AURKA","MAP3K9","MAPK12","PCTK1","PCTK3","MERTK","PIK3C2A","MYLK2","BCR","DCK","ACVRL1","ADRB2","GSG2","HUNK","AVPR1B")
my.2nd.hits<-c("DAPK1","TEX14","DUSP10")
conf.hits<-c("PLK1","CSNK1D","IRS1","MAP3K7","SNRK","WEE1","STK39","PRKAA2","PRKAR2B","JAK2","COPB2","ERBB4","CNKSR1","CHEK1")

######JANE
my.hits<-c("SFRS7","DDX21","C-MYC") # groth in bothin both rel and RG scores
my.hits.rel("HEATR1","IFRD2","SNX5","SH3RF2","GSTTP1","CAPN3","S100G","CYB5A","RLN2","RP1-21O18.1")
my.hits.RG("MSI2","SLC25A32","PWP2","LARP4","DUSP16","NUP25","HNRNPAB","SLC25A24","BCL2","CHEK1","TSSK6","CDK9","NLN")
my.2nd.hits<-c("ABAT","GAR1","GLA","DTL","WISP2","FARSB","KIF12","PTP4A2","C4orf16","STC2","GNB1","UMPS","MINA","UCHL5","METTL1","DKC1","SLC16A6","HNRNPD","ZWINT") # growth in EX not in ICI in both rel and RG scores
conf.hits<-c("HSPCN","CCND1")

######JANE VALIDATION
my.hits<-c("SFRS7","DDX21","C-MYC") # groth in bothin both rel and RG scores
## my.hits.rel("HEATR1","IFRD2","SNX5","SH3RF2","GSTTP1","CAPN3","S100G","CYB5A","RLN2","RP1-21O18.1")
## my.hits.RG("MSI2","SLC25A32","PWP2","LARP4","DUSP16","NUP25","HNRNPAB","SLC25A24","BCL2","CHEK1","TSSK6","CDK9","NLN")
## my.2nd.hits<-c("ABAT","GAR1","GLA","DTL","WISP2","FARSB","KIF12","PTP4A2","C4orf16","STC2","GNB1","UMPS","MINA","UCHL5","METTL1","DKC1","SLC16A6","HNRNPD","ZWINT")
                                        # growth in EX not in ICI in both rel and RG scores
conf.hits<-c("HSPCN","CCND1")
 ######Heatmap

 ######Heatmap
 ####advanced heatmap  biocLite(c("gplots"))
 library("gplots")


par.opts.start<-par(no.readonly=TRUE)$mar
## [1] 5.1 4.1 4.1 2.1
par(mar=c(5.1, 4.1,4.1,2.1),mgp=c(3,1,0),las=2)  #c(bottom, left, top, right)
 par(font.axis=2)
 par(font=2)




## all.screens[1:5,]


## dim(data.all)
## data.all[1:5,]
## length((all.screens[,"Cell Line"]))

## the.cell.lines<-unique(all.screens[,"Cell Line"])
## the.cell.lines


## data.cell.lines<-{}
## for( i in 1:length(the.cell.lines)){
## temp<- data.all[all.screens[,"Cell Line"] %in% the.cell.lines[i],]
## colnames(temp)<-paste(the.cell.lines[i],colnames(temp),sep="_")
## data.cell.lines<-cbind(data.cell.lines,temp)
## }
## data.cell.lines[1:5,]  ## hit number NOT numerical

a.test<-"ALG8"
all.screens[all.screens[,"Gene"]==a.test,]
wanted<-all.screens[,"Gene"]==a.test

data.names
ii<-3
load(data.names[ii])
 the.test<-summary[,"Symbol"] %in% a.test

hits.full.rep[the.test,]
num.reps<-c(num.reps,sum(the.test))
summary[the.test,]
hits.low.rep[the.test,]
hits.high.rep[the.test,]


#############################################
#################  MAKE numerical data out of hit data:
old.names<-names(a.screen.cuts)
## names(a.screen.cuts)<-c("EtOH","ICI","EX")

cell.lines<-all.screens[,"Cell Line"]
genes<-all.screens[,"Gene"]

##### turn the medians to numeric put in data.num
all.screens[1:5,]
data.num<-all.screens[,paste("<",targets,">",sep="")]
data.num[1:5,]
the.dim.num<-dim(data.num)
data.num<-as.numeric(data.num)
dim(data.num)<-the.dim.num
rownames(data.num)<-cell.lines
colnames(data.num)<-targets

##### turn the counts to numeric put in data.all , scale by the num reps
data.all<-all.screens[,targets]
data.all[1:5,]
the.dim.all<-dim(data.all)
data.all<-as.numeric(data.all)
dim(data.all)<-the.dim.all
rownames(data.all)<-cell.lines
colnames(data.all)<-targets

data.nreps<-all.screens[,paste("num reps",targets,sep=".")]
data.nreps[1:5,]
the.dim.all<-dim(data.nreps)
data.nreps<-as.numeric(data.nreps)
dim(data.nreps)<-the.dim.all
rownames(data.nreps)<-cell.lines
colnames(data.nreps)<-targets

data.all<-data.all/data.nreps

#####data.nreps can be zero making NaN so set these to zero data.all MUST also be zero in those cases

for(i in 1:length(targets)){
  the.target<-targets[i]
  data.all[!is.finite(data.all[,the.target]),the.target]<-0
  data.num[!is.finite(data.num[,the.target]),the.target]<-0
}
sum(!is.finite(data.num[,6]))
sum(is.na(data.num[,6]))
## rownames(data.all)<-cell.lines
## rownames(data.all)<-genes


data.num[1:5,]
data.all[1:5,]
## data.num<-data.num*data.all #### not perfect as not the true average. really only works here as 

pos.hits<-data.all>0 ## positive hits
neg.hits<-data.all<0 ## negative hits (recall I  *1 to the number of negaitve hits (low.counts above)
NA.hits<-is.na(data.num)
pos.hits[1:5,]
neg.hits[1:5,]
NA.hits[1:5,]
data.new<-matrix(data=0,nrow=dim(data.num)[1],ncol=dim(data.num)[2])
colnames(data.new)<-colnames(data.num)
rownames(data.new)<-rownames(data.num)

data.high<-data.new
data.low<-data.new

the.row<-rownames(data.num)  ## the rownames of data.num are the screen names as per same as a.screen.cuts
the.col<-colnames(data.num)
for( i in 1:dim(data.num)[1]){
  for( j in 1:dim(data.num)[2]){
    if(pos.hits[i,j]){
      ## a.screen.cuts[[the.row[i]]]$high.cut[the.col[j]]
      ## data.num[i,j]
      data.high[i,j]<-abs(data.num[i,j]-a.screen.cuts[[the.row[i]]]$high.cut[the.col[j]])
    }
     if(neg.hits[i,j]){
      ## a.screen.cuts[[the.row[i]]]$low.cut[the.col[j]]
      ## data.num[i,j]
      data.low[i,j]<--1*abs(a.screen.cuts[[the.row[i]]]$low.cut[the.col[j]]-data.num[i,j])
    }
        if(NA.hits[i,j]){
      ## data.num[i,j]
      data.new[i,j]<-0
    }
    
  } #cols
} #rows
      
## rownames(data.high)<-genes
## rownames(data.low)<-genes
## rownames(data.new)<-genes
data.all[1:5,]
data.high[1:5,]
data.low[1:5,]
data.num[1:5,]

data.high.thresh<-data.high ## use for unscaled calculations  data.high<-data.high.thresh
data.low.thresh<-data.low ## use for unscaled calculations   data.low<-data.low.thresh


####################################get a feel for sizes and any data issues x<-data.high[rownames(data.high)==the.screen,6]

run.screens<-names(a.screen.cuts)
for(i in 1:length(run.screens)){
  the.screen<-run.screens[i]
  print( the.screen)
    print("High")
print( apply(data.high[rownames(data.high)==the.screen,],2,function(x){ test<-median(x[x!=0]); if(sum(x!=0)==0){test<-0}; test } )  ) ## test is all zeros then get trival NAs
  print("Low")
print( apply(data.low[rownames(data.high)==the.screen,],2, function(x){ test<-median(x[x!=0]); if(sum(x!=0)==0){test<-0}; test } )  )
  print("-------------------------------------------------")
}

sum(is.na(data.high[,5]))

for(i in 1:length(targets)){
  the.target<-targets[i]
  print(the.target)
      print("High")
print( apply(data.high,2,function(x) median(x[x!=0])) )
  print("Low")
print( apply(data.low,2,function(x) median(x[x!=0]))  )
    print("-------------------------------------------------")
  
}


apply(data.high,2,function(x) median(x[x!=0]))
apply(data.low,2,function(x) median(x[x!=0]))

apply(data.high,2,function(x) max(x[x!=0]))
apply(data.high,2,function(x) min(x[x!=0]))

apply(data.low,2,function(x) max(x[x!=0]))
apply(data.low,2,function(x) min(x[x!=0]))

############# PLOT data righ click to move to the next plot
for(i in 1:length(targets)){
  the.target<-targets[i]
  plot(1:dim(data.high)[1],sort(data.high[,the.target]),xlab="Rank",ylab=the.target,xlim=(c(0,(dim(data.high)[1]+ dim(data.high)[1]/7))))
  selected.data<-identify(1:dim(data.high)[1],sort(data.high[,the.target]),labels=paste(all.screens[,"Gene"],all.screens[,"Cell Line"],sep=":"),col="red",cex=1,atpen='TRUE')

     }

for(i in 1:length(targets)){
  the.target<-targets[i]
  plot(1:dim(data.low)[1],sort(data.low[,the.target]),xlab="Rank",ylab=the.target,xlim=(c(0,(dim(data.high)[1]+ dim(data.high)[1]/7))))
  selected.data<-identify(1:dim(data.low)[1],sort(data.low[,the.target]),labels=paste(all.screens[,"Gene"],all.screens[,"Cell Line"],sep=":"),col="red",cex=1,atpen='TRUE')

     }

################## LIMT RANGE if too extreme if needed
## data.high[data.high[,"Relative S+G2+_gt_4N"]>5,"Relative S+G2+_gt_4N"]<-5
## data.high[data.high[,"Gr-S+G2+ _gt_4N_over_G1"]>5,"Gr-S+G2+ _gt_4N_over_G1"]<-5
data.high[data.high[,"RG_over_G_over_EtoH"]>5,"RG_over_G_over_EtoH"]<-5
data.high[data.high[,"Gr-G1"]>0.1,"Gr-G1"]<-0.1
data.high[data.high[,"Z-all Green"]> 15,"Z-all Green"]<-15
data.high[data.high[,"RG_over_G"]> 0.5,"RG_over_G"]<-0.5
data.high[data.high[,"RG_over_G_Z"]> 30,"RG_over_G_Z"]<-30

data.low[data.low[,"Z-all Green"]< -3,"Z-all Green"]<--3
data.low[data.low[,"RG_over_G"]< -0.05,"RG_over_G"]<- -0.05
data.low[data.low[,"RG_over_G_Z"]< -2,"RG_over_G_Z"]<- -2
data.low[data.low[,"Gr-G1"]< -0.45,"Gr-G1"]<- -0.45
data.low[data.low[,"RG_over_G_over_EtoH"]< -0.2,"RG_over_G_over_EtoH"]<- -0.2
########################## RESCALE THE DATA 0-1 ans 0- (-1) 


## data.high
the.medians<-apply(data.high,2,function(x) median(x[x!=0]))
the.sd<-apply(data.high,2,function(x) sd(x[x!=0]))
the.max<-apply(data.high,2,function(x) max(x[x!=0]))
the.min<-apply(data.high,2,function(x) min(x[x!=0]))
the.range<-the.max-the.min


for(i in 1:length(targets)){
data.high[,targets[i]]<-(data.high[,targets[i]]-the.min[targets[i]]) / the.range[targets[i]]
}
data.high[!pos.hits]<-0

data.high[1:5,]



## data.low
the.medians<-apply(data.low,2,function(x) median(x[x!=0]))
the.sd<-apply(data.low,2,function(x) sd(x[x!=0]))
the.max<-apply(data.low,2,function(x) max(x[x!=0]))
the.min<-apply(data.low,2,function(x) min(x[x!=0]))
the.range<-the.max-the.min

for(i in 1:length(targets)){
## data.low<-(data.low[,names(the.medians)[i]]-the.medians[names(the.medians)[i]]) / the.sd[names(the.medians)[i]]
data.low[,targets[i]]<-(data.low[,targets[i]]-the.max[targets[i]]) / the.range[targets[i]]
}
data.low[!neg.hits]<-0

data.low[1:5,]
data.low.thresh[1:5,]


data.high[1:5,]
data.low[1:5,]
#################### RESTART :  FOR MULTI PARAMETER  chose which one to use
data.new<-data.high+data.low   ### numerical data w.r.t thrshold  AND scaled

## data.new<-data.high
## data.new[,"Gr-G1"]<--1*data.low[,"Gr-G1"] # cause low Gr-G1 is high proliferation
## rownames(data.new)<-cell.lines
rownames(data.new)<-genes
do.pca<-TRUE
################### data.new contains the numerical data


########################## RESTART : FOR ONE ATTRIBUTE ONLY
######### Get the single attribute and the reshape to 
data.new<-data.high+data.low  ### numerical data w.r.t thrshold WITH scaled none hits as zero #ued this for Liz
## data.new<-data.high.thresh+data.low.thresh  ### numerical data w.r.t thrshold NOT scaled non hits as zero

data.new.cell.lines<-{}
run.screens<-names(a.screen.cuts)
for(i in 1:length(run.screens)){
  the.screen<-run.screens[i]
temp<- data.new[all.screens[,"Cell Line"] %in% the.screen,]
colnames(temp)<-paste(the.screen,colnames(temp),sep="_")
data.new.cell.lines<-cbind(data.new.cell.lines,temp)
}
all.screens[1:7,]
data.new[1:7,]
data.new.cell.lines[1:2,]



targets
data.new.cell.lines[1:5,]
wanted.target<-targets[8]
wanted.target
data<-data.new.cell.lines[,grepl(paste(wanted.target,"$",sep=""),colnames(data.new.cell.lines))]
colnames(data)<-gsub(paste("_",wanted.target,sep=""),"",colnames(data))
data[1:5,]

do.pca<-TRUE ## use if only 3 variables
############################################### RESTART with new targets


############ OR if no scalling is needed that is no rescale to cutoffs which I think is mor appropriate here- note non hits are zero
## wanted<-targets[7]
## wanted
## data.cell.lines[1:5,]

## data<-data.cell.lines[,grepl(wanted,colnames(data.cell.lines))]
## colnames(data)<-gsub(paste("_",wanted,sep=""),"",colnames(data))
## data[1:5,]
## no.pca<-TRUE
########################################straight up distance calculations:
## genes<-unique(rownames(data.all))
## posns<-match(genes,rownames(data.all))
## the.dists<-{}
## for( i in 1:length(targets)){
##   wanted<-grepl(targets[i],colnames(data.all))
##   adist<-dist(data.all[,wanted])
##   adist<-as.matrix(adist)
##   adist<- diag(adist[posns+1,posns])
##   if(i==1){the.dists<-adist}else{the.dists<-cbind(the.dists,adist)}
## }
## colnames(the.dists)<-targets
## ### this is absolute distance getween genes


########################## choose one
the.cell.lines<-names(a.screen.cuts)
the.cell.lines

rownames(data.new)<-genes # rownames(data.new)<-cell.lines

data3<- data.all[all.screens[,"Cell Line"] %in% cell.lines[3],]
data2<- data.all[all.screens[,"Cell Line"] %in% cell.lines[2],]
data


####jane
cell.lines
[1] "EtOH" "ICI"  "EX"  
data3<- data.all[all.screens[,"Cell Line"] %in% cell.lines[3],]
data2<- data.all[all.screens[,"Cell Line"] %in% cell.lines[2],]
data1<- data.all[all.screens[,"Cell Line"] %in% cell.lines[1],]

data.new[1:15,]
the.cell.lines
[1] "EtOH-a" "ICI-b"  "EX-c"
> cell.lines
[1] "a"  "b"  "c"  "4a" "4b" "4c"
data3<- data.new[all.screens[,"Cell Line"] %in% the.cell.lines[3],] # c
data2<- data.new[all.screens[,"Cell Line"] %in% the.cell.lines[2],] # b 
data1<- data.new[all.screens[,"Cell Line"] %in% the.cell.lines[1],] # a

new.genes<-genes[all.screens[,"Cell Line"] %in% the.cell.lines[1]]

data<-data3+data2 # Jane
data[1:5,]
data<-data[,c(1,2,5,6)] # jane multi-parameter
data<-data[,c("Num-cells", "Z-all Green", "RG_over_G", "RG_over_G_Z")]


data<-data3-data2  ## Fawau caski-c33q()

data<-abs(data3-data2)


data3[1:5,]
data2[1:5,]
data1[1:5,]
max(data[,"Relative S+G2+_gt_4N"])
min(data[,"Relative S+G2+_gt_4N"])

######################################################################### BEGIN DIMENSIONAL ANALYSIS

ColSideColors=rep("black",times=dim(data)[1])
ColSideColors[rownames(data) %in% conf.hits]<-"green"
## ColSideColors[rownames(data) %in% "PLK1"]<-"green"
ColSideColors[rownames(data) %in% my.2nd.hits]<-"orange"
ColSideColors[rownames(data) %in% my.hits]<-"red"
rownames(data)[rownames(data) %in% my.hits]

colsep<-grep(FALSE, ColSideColors=="black")
sepcolor<-ColSideColors[colsep]
               
data[rownames(data) %in% conf.hits,]
data2[rownames(data2) %in% conf.hits,]
data3[rownames(data3) %in% conf.hits,]

## colorss.ori<-colorss
## #colorss<-bluered(max(data)*4)
par.opts.start<-par(no.readonly=TRUE)$mar
## [1] 5.1 4.1 4.1 2.1
par(mar=c(5.1, 4.1,4.1,2),mgp=c(3,1,0),las=2)  #c(bottom, left, top, right)
 par(font.axis=2)
 par(font=2)

########################################################################## DO HEATMAP
##Choose  number of clusters
centers<-16
colours<-rainbow(centers)

if(!do.pca){

gene.cl<-hclust(dist((data)))

gene.cl.members<-cutree(gene.cl,k=centers)
the.heat<-heatmap.2(t(data),margins=c(5,9),trace="row",ColSideColors=colours[gene.cl.members],keysize=1.0,cex.lab=1.5,main=paste("Heatmap of hits for:",wanted.target,sep=" "))
########################################################################
}else{
the.heat<-heatmap.2(t(data),margins=c(10,15),trace="row",ColSideColors=ColSideColors,keysize=1.0,cex.lab=1.5,main="Heatmap of differential hits")

#### to color side bar with colors from done PCA below
## the.heat<-heatmap.2(t(data),margins=c(10,15),trace="row",ColSideColors=colours[the.cl$cluster],keysize=1.0,cex.lab=1.5,main="Heatmap of differential hits")


}

plot.name<-"Six paramenter heatmap"
plot.name<-"Growth Six paramenter heatmap"
plot.name<-"four paramenter heatmap"
plot.name<-"high only no DNA clustering MCF10a new"
plot.name<-"Gr-G1 clustering by cell line"
plot.name<-"RG_over_G clustering by cell line"
plot.name<-"RG_over_G_over_EtoH clustering by cell line"
savePlot("neive clustering.jpeg",type="jpeg")

savePlot(paste(plot.name,"jpeg",sep="."),type="jpeg")
savePlot(paste(plot.name,"tiff",sep="."),type="tiff")

par(mar=c(5.1, 4.1,4.1,2.1),mgp=c(3,1,0),las=2) # reset plots
if(!do.pca){
ii<-2
jj<-3
colours <-  rainbow(centers)
plot(range(data[,ii]),range(data[,jj]),xlab=colnames(data)[ii],ylab=colnames(data)[jj],main=paste("Spectral clustering of hits for:",wanted.target,sep=" "))

plot(range(data[,ii]),range(data[,jj]),xlab=colnames(data)[ii],ylab=colnames(data)[jj],xlim=c(-0.2,0.2),ylim=c(-0.2,0.2),main=paste("Spectral clustering of hits for:",wanted.target,sep=" "))

text(data[,ii],data[,jj],label=rownames(data),col=colours[gene.cl.members],cex=0.75,font=2)

}else{
  
the.pca <- prcomp(data,scale = TRUE) # for hits/genes
attributes(the.pca )
dim(the.pca$x)

the.pca.var <- round(the.pca$sdev^2 / sum(the.pca$sdev^2)*100,2)
plot(c(1:length(the.pca.var)),the.pca.var,type="b",xlab="# components",ylab="% variance",main="Scree Plot for Hits",col="red",cex=1.5,cex.lab=1.5)
savePlot(paste("Scree Plot for",plot.name,"jpeg",sep="."),type="jpeg")

the.cl<-kmeans(the.pca$x[,1:2],centers=centers,iter.max=1000,nstart=50) #Do kmeans

#########full scale
plot(range(the.pca$x[,1]),range(the.pca$x[,2]),xlab="PCA1",ylab="PCA2",main="Spectral clustering of differential hits",col="white",font=2)

#########zoom in
## plot(range(the.pca$x[,1]),range(the.pca$x[,2]),xlim=c(-3,2),ylim=c(-1.5,1.5),xlab="PCA1",ylab="PCA2",main="Spectral clustering of differential hits",font=2)
##  plot(range(the.pca$x[,1]),range(the.pca$x[,2]),xlim=c(-1.5,2),ylim=c(-1.25,1.5),xlab="PCA1",ylab="PCA2",main="Spectral clustering of differential hits",font=2)

text(the.pca$x[,1],the.pca$x[,2],label=rownames(the.pca$x),col=colours[the.cl$cluster],cex=0.85,font=2)
}

########################################  HIGHTLIGH some gene if you wish
## dotchart(the.pca$x[,1:3],labels=as.character(rownames(the.pca$x)))
################### choose one
the.test<-my.hits
color<-"red"

the.test<-my.2nd.hits
color<-"orange"

the.test<-conf.hits
color<-"green"
###########################

if(!do.pca){
wanted<-rownames(data) %in% the.test
}else{
wanted<-rownames(the.pca$x) %in% the.test
}

if(!do.pca){
points(data[wanted,ii],data[wanted,jj],col=color,cex=5.0)
}else{
points(the.pca$x[wanted,1],the.pca$x[wanted,2],col=color,cex=5.0)
}

####################################
plot.name<-"Growth Six parameter spectral"
plot.name<-"Six parameter spectral"
plot.name<-"Four parameter spectral"
plot.name<-"Gr-G1 spectral"
plot.name<-"RG_over_G-over EtoH spectral"
plot.name<-"RG_over_G"

plot.name<-"Growth Six parameter spectral zoom in"
plot.name<-"Four parameter spectral zoom in"
plot.name<-"Gr-G1 spectral zoom in"
plot.name<-"RG_over_G zoom in"
plot.name<-"RG_over_G-over EtoH spectral zoom in"

savePlot(paste(plot.name,"jpeg",sep="."),type="jpeg")
savePlot(paste(plot.name,"tiff",sep="."),type="tiff")

## load((paste("Six parameter spectral","RData",sep="."))
save.image(paste(plot.name,"RData",sep=".")) ## keep so have the cluster data

library(scatterplot3d)  # biocLite("scatterplot3d")
#c(-0.2,0.4),c(-0.5,0.5),c(0,0.25) RG-over G
if(!do.pca){
s3d<-scatterplot3d(range(data[,1]),range(data[,2]),range(data[,3]),xlab=colnames(data)[1],ylab=colnames(data)[2],zlab=colnames(data)[3],main=paste("Spectral clustering of hits for:",wanted.target,sep=" "),angle=120, cex.lab=1.5,cex.axis=1.5)


s3d<-scatterplot3d(c(-0.2,0.2),c(-0.2,0.2),c(-0.2,0.2),xlab=colnames(data)[1],ylab=colnames(data)[2],zlab=colnames(data)[3],main=paste("Spectral clustering of hits for:",wanted.target,sep=" "),angle=120, cex.lab=1.5,cex.axis=1.5,grid=TRUE)

## text(s3d$xyz.convert(range(data[wanted,1]),range(data[wanted,2]),range(data[wanted,3])),col="red",cex=5.0)xlim=c(-0.2,0.5),ylim=c(-0.02,0.1)
text(s3d$xyz.convert(data[,1],data[,2],data[,3]),label=rownames(data),col=colours[gene.cl.members],cex=1.0)
## points(data[wanted,1],data[wanted,2],col=color,cex=5.0)
points(s3d$xyz.convert(data[wanted,1],data[wanted,2],data[wanted,3]),col=color,cex=5.0)

}else{
## scatterplot3d(the.pca$x[,1],the.pca$x[,2],the.pca$x[,3],angle=100,pch=as.character(rownames(the.pca$x)))

######## ALL DATA
s3d<-scatterplot3d(range(the.pca$x[,1]),range(the.pca$x[,2]),range(the.pca$x[,3]),xlab="PCA1",ylab="PCA2",zlab="PCA3",main="Spectral clustering of differential hits",angle=120,color="white")


##########ZOOM IN 
## s3d<-scatterplot3d(c(-4,4),c(-2,1),c(0,2),xlab="PCA1",ylab="PCA2",zlab="PCA3",main="Spectral clustering of differential hits",angle=120,color="white")


text(s3d$xyz.convert(the.pca$x[,1],the.pca$x[,2],the.pca$x[,3]),label=rownames(the.pca$x),col=colours[the.cl$cluster],cex=1.0)
## points(the.pca$x[wanted,1],the.pca$x[wanted,2],col=color,cex=5.0)
points(s3d$xyz.convert(the.pca$x[wanted,1],the.pca$x[wanted,2],the.pca$x[wanted,3]),col=color,cex=7.0)

}
plot.name<-"3D Growth Six parameter spectral"
plot.name<-"3D Four parameter spectral"
plot.name<-"3D Gr-G1 spectral"
plot.name<-"3D RG_over_G spectral"
plot.name<-"3D RG_over_G_over_EtOH spectral"

plot.name<-"3D Growth Six parameter spectral zoom in"
plot.name<-"3D Four parameter spectral zoom in"
plot.name<-"3D Gr-G1 spectral zoom in"
plot.name<-"3D RG_over_G spectral zoom in"
plot.name<-"3D RG_over_G_over_EtOH spectral"

savePlot(paste(plot.name,"jpeg",sep="."),type="jpeg")
savePlot(paste(plot.name,"tiff",sep="."),type="tiff")


###############  cleck a cluster with a ONE  member
the.member<-"TREML1"
the.member<-"CLK2"
if(!do.pca){
group<-gene.cl.members[gene.cl.members==gene.cl.members[names(gene.cl.members)==the.member]]
group
the.pca$x[names(group),]
}else{
group<-the.cl$cluster[the.cl$cluster==the.cl$cluster[names(the.cl$cluster)==the.member]]
group
the.pca$x[names(group),]
}
###################################

the.test<-my.hits
color.hit<-"red"
wanted.hits<-rownames(the.pca$x) %in% the.test

the.test<-my.2nd.hits
color.2nd.hit<-"orange"
wanted.2nd.hits<-rownames(the.pca$x) %in% the.test

the.test<-conf.hits
color.conf<-"green"
wanted.control<-rownames(the.pca$x) %in% the.test

load("Four parameter.RData")
load("Gr-G1 new .RData")
load("RG_over_G new .RData")
load("RG_over_G_over_EtOH new .RData")


the.spin<-seq(from=0,to=360,by=5)
 Sys.sleep(5)

if(!do.pca){
  for (i in 1:length(the.spin)){
  print(i)
s3d<-scatterplot3d(range(data[,1]),range(data[,2]),range(data[,3]),xlab=colnames(data)[1],ylab=colnames(data)[2],zlab=colnames(data)[3],main=paste("Spectral clustering of hits for:",wanted.target,sep=" "),angle=the.spin[i])
## text(s3d$xyz.convert(range(data[wanted,1]),range(data[wanted,2]),range(data[wanted,3])),col="red",cex=5.0)
text(s3d$xyz.convert(data[,1],data[,2],data[,3]),label=rownames(data),col=colours[gene.cl.members],cex=0.75)
## points(data[wanted,1],data[wanted,2],col=color,cex=5.0)
points(s3d$xyz.convert(data[wanted.hits,1],data[wanted.hits,2],data[wanted.hits,3]),col=color.hit,cex=5.0)
points(s3d$xyz.convert(data[wanted.2nd.hits,1],data[wanted.2nd.hits,2],data[wanted.2nd.hits,3]),col=color.2nd.hit,cex=5.0)
points(s3d$xyz.convert(data[wanted.control,1],data[wanted.control,2],data[wanted.control,3]),col=color.conf,cex=5.0)
Sys.sleep(1)
  ## savePlot(paste("3d clustering angle=",the.spin[i],".jpeg",sep=""),type="jpeg")
}
  
}else{

for (i in 1:length(the.spin)){
  print(i)
s3d<-scatterplot3d(range(the.pca$x[,1]),range(the.pca$x[,2]),range(the.pca$x[,3]),xlab="PCA1",ylab="PCA2",zlab="PCA3",main="Spectral clustering of differential hits",angle=the.spin[i])
## text(s3d$xyz.convert(range(the.pca$x[wanted,1]),range(the.pca$x[wanted,2]),range(the.pca$x[wanted,3])),col="red",cex=5.0)
text(s3d$xyz.convert(the.pca$x[,1],the.pca$x[,2],the.pca$x[,3]),label=rownames(the.pca$x),col=colours[the.cl$cluster],cex=0.75)
## points(the.pca$x[wanted,1],the.pca$x[wanted,2],col=color,cex=5.0)
points(s3d$xyz.convert(the.pca$x[wanted.hits,1],the.pca$x[wanted.hits,2],the.pca$x[wanted.hits,3]),col=color.hit,cex=5.0)
points(s3d$xyz.convert(the.pca$x[wanted.2nd.hits,1],the.pca$x[wanted.2nd.hits,2],the.pca$x[wanted.2nd.hits,3]),col=color.2nd.hit,cex=5.0)
points(s3d$xyz.convert(the.pca$x[wanted.control,1],the.pca$x[wanted.control,2],the.pca$x[wanted.control,3]),col=color.conf,cex=5.0)
Sys.sleep(1)
  ## savePlot(paste("3d clustering angle=",the.spin[i],".jpeg",sep=""),type="jpeg")
}

}


savePlot("3d clustering.jpeg",type="jpeg")
savePlot("3d clustering.tiff",type="tiff")



####################################define specific hits if you want
if(!do.pca){
hits<-rep("",times=dim(data)[1])
names(hits)<-rownames(data)
}else{
hits<-rep("",times=dim(the.pca$x)[1])
names(hits)<-rownames(the.pca$x)

best.names<-strsplit(names(hits),split=".",fixed=TRUE)
best.names<-unlist(lapply(best.names,function(x) x[1]))
names(hits)<-best.names

}


## JAne
hits[wanted.hits]<-"ICI & EX  HIT"
hits[wanted.2nd.hits]<-"EX HIT"
hits[wanted.control]<-"CONTROL HIT"

### JAne validation
screen.hits<-read.delim("selected hits.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ## JANe validation
for(i in 1:dim(screen.hits)){
  print(i)
  screen.hits[,i]<-gsub("*","",screen.hits[,i],fixed=TRUE)
  temp<-strsplit(screen.hits[,i],split=" ")
   screen.hits[,i]<-unlist(lapply(temp,function(x) x[1]))
}
screen.hits[1:5,] # screen.hits[,i]
hits[screen.hits[,1]]<-paste(hits[screen.hits[,1]],"4 Param",sep="|")  ### ICI & EX  HIT"
hits[screen.hits[,2]]<-paste(hits[screen.hits[,2]],"Ratio EtOh",sep="|")
hits[screen.hits[,3]]<-paste(hits[screen.hits[,3]],"%G1",sep="|")
hits[screen.hits[,4]]<-paste(hits[screen.hits[,4]],"RG/R",sep="|") 

hits

hits<-hits[best.names]
######################################################

all.screens[1:6,]
data.high[1:5,]

ann[1:5,]
unique.genes<-unique(ann[,"Symbol"])
posns<-match(unique.genes,ann[,"Symbol"])
sum(is.na(posns))
use.ann<-ann[posns, !(colnames(ann) %in% core.ann[core.ann !="Symbol"]) ]

use.ann<-as.matrix(use.ann) ## so can use dupcalted row names when making all.summary
rownames(use.ann)<-use.ann[,"Symbol"]
use.ann[1:5,]

####################
##can add more annotation here

best.names<-strsplit(use.ann[,"Symbol"],split=".",fixed=TRUE)
best.names<-unlist(lapply(best.names,function(x) x[1]))



use.ann[1:5,]
summary<-use.ann  ## so can use annote with drug.r
colnames(summary)[colnames(summary)=="Symbol"]<-"Gene Symbol"

summary[,"Gene Symbol"]<-best.names
summary[1:5,]
col.names<-colnames(summary)
summary<-cbind(summary,NA,NA,NA,NA,NA,NA,NA)
colnames(summary)<-c(col.names,"Description","Known Cancer Gene","Known in Breast","COSMIC-Primary Hist","COSMIC-Hist SubType","COSMIC-Primary Site","COSMIC-Combinations")

### GO TO OTHER CODE

summary[1:5,]

colnames(summary)[colnames(summary)=="Gene Symbol"]<-"Symbol"
use.ann<-cbind(use.ann,summary[,c("Description","Known Cancer Gene","Known in Breast","COSMIC-Primary Hist","COSMIC-Hist SubType","COSMIC-Primary Site","COSMIC-Combinations")])

use.ann[1:5,]
###################




rownames(use.ann)
########################################################################  BUILD SUMMARY
all.summary[1:5,]
all.screens[1:5,]

if(!do.pca){
all.summary<-cbind(all.screens,gene.cl.members[all.screens[,"Gene"]],hits[all.screens[,"Gene"]],data,data[all.screens[,"Gene"],],use.ann[all.screens[,"Gene"],])
colnames(all.summary)
colnames(all.summary)<-c(colnames(all.screens),"cluster ID","Hit Class",paste("SCALED",colnames(data.high),sep="-"),colnames(data),colnames(use.ann))
all.summary[1:5,]
}else{
all.summary<-cbind(all.screens,the.cl$cluster[all.screens[,"Gene"]],hits[all.screens[,"Gene"]],data[all.screens[,"Gene"],],the.pca$x[all.screens[,"Gene"],],use.ann[all.screens[,"Gene"],])
colnames(all.summary)
colnames(all.summary)<-c(colnames(all.screens),"cluster ID","Hit Class",paste("SCALED",colnames(data.high),sep="-"),colnames(the.pca$x),colnames(use.ann))
all.summary[1:2,]
}

     dim(all.summary)
getwd()
file.name<-"Combined_Hits_SUMMARY 6 Growth  parameters.txt"
file.name<-"Combined_Hits_SUMMARY 6 parameters.txt"
file.name<-"Combined_Hits_SUMMARY 4 parameters.txt"
file.name<-"Combined_Hits_SUMMARY Gr-G1.txt"
file.name<-"Combined_Hits_SUMMARY RG_over_G.txt"
file.name<-"Combined_Hits_SUMMARY RG_over_G_over EtOH.txt"
write.table(all.summary,file.name,row.names=FALSE,sep="\t")
save.image("working_final growth.RData")
## setwd( "/media/Bioinform-D/Research/Cellomics/Ca screen/Ca screen 2")


all.symmary.6.parameter<-all.summary
all.symmary.gr.g1<-all.summary
all.symmary.RG.over.G<-all.summary
all.symmary.RG.over.G.over.EtOH<-all.summary

all.summary[all.summary[,"cluster ID"]=="12",1:5]

all.summary.1st[1:5,]
all.summary<-read.delim("Combined_Hits_SUMMARY 4 parameters.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

save.image("all work.RData")

data.names<-c("all.symmary.4.parameter","all.symmary.gr.g1","all.symmary.RG.over.G","all.symmary.RG.over.G.over.EtOH")

clusters<-{}
for(i in 1:length(data.names)){
data<-eval(as.name(data.names[i]))
print(sum(data[,"Gene"]!=genes))
clusters<-cbind(clusters,data[,"cluster ID"])
## assign(data.names[i],value=data)
}
colnames(clusters)<-paste("Cluster ID:", gsub("all.symmary.","",data.names),sep="")

combined.summary<-cbind(all.summary[,1:19],clusters,all.summary[,21:41])
write.table(combined.summary,"combined.summary.txt",row.names=FALSE,sep="\t")


all.summary[1:10,"Cell Line"]
all.summary[,"Cell Line"]<-gsub("a","Etoh",all.summary[,"Cell Line"])
all.summary[,"Cell Line"]<-gsub("b","ICI",all.summary[,"Cell Line"])
all.summary[,"Cell Line"]<-gsub("c","EX",all.summary[,"Cell Line"])
all.summary[,"Cell Line"]<-gsub("4","T47D ",all.summary[,"Cell Line"])
                                
################################END 





































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


##################################heatmap and rescale
rowclust <- hclust(mydDist(mat, "myMetric"), "myLinkage")
colclust <- hclust(myDist(mat, "anotherMetric"), "anotherLinkage")
stdizedmat <- t(scale(t(mat)))
clipper <- 3
stdizedmat[stdizedmat < -clipper]<- -clipper
stdidezmat[stdizedmat > clipper] <- clipper

heatmap(stdizedmat,
       Rowv=as.dendrogram(rowclust),
       Colv=as.dendrogram(colclust),
       scale='none', zlim=c(-clipper, clipper))
###########################
