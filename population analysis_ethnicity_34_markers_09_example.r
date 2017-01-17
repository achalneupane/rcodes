
##################################################################################################
############################## START HERE  all.markers
############################## START HERE  hapmap.calls train.extended tain.ann.extended
############################## START HERE
############################## START HERE
############################## START HERE
############################## START HERE
/media/old-scratch/media/Bioinform-D/Research/Matt Brown/Popluation Statification SNPs/population analysis_ethnicity_34_markers_09_example.r
/media/old-scratch/media/Bioinform-D/Research/Matt Brown/Popluation Statification SNPs/Ethnicity.prediction.input.RData
library(MASS)
library(randomForest) 
library(ipred)   
library(e1071)

working.directory<-"/media/old-scratch/media/Bioinform-D/Research/Matt Brown/Popluation Statification SNPs"
setwd(working.directory) 
load("Ethnicity.prediction.input.RData")  ## data with known ethnicitys to train


#all.samples<-read.delim("FHCRC_TASCOG-Plate01.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
all.samples<-read.delim("FHCRC_Plate01.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
#sum(tapply(all.samples[,"sample"],all.samples[,"sample"] ,length)!=1 )   # should be zero
all.samples[1:5,]


############################ data in sequnome format
##         sample rs1923416 rs1113480 rs953035 rs1890191 rs1831024 rs718686
## 1 H-2401900139       F F       V V      F F       V F       V V      F F
## 2 H-2401990803       F F       V V      F F       V F       V V      F F
## 3 H-2401992129       F F       V V      F F       F F       V F      V F
## 4 H-2401992790       F F       V F      F F       V F       V V      F F
## 5 H-2401993010       F F       V V      F F       F F       V V      V F

##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
################################# this section to reformat data

rownames(all.samples)<-all.samples[,1]
all.samples<-all.samples[,-1]
colnames(all.samples)<-gsub("rs","RS",colnames(all.samples))

new.order<-names(hapmap.calls[colnames(all.samples)])
all.samples<-all.samples[,new.order]
markers34[colnames(all.samples),3:5]
markers34<-markers34[colnames(all.samples),]


dim(all.samples)
all.samples[1:5,]
i<-1
all.samples[i,]
all.samples[,i]
##### chnge from VF to allelles as defined by markers:for each snp SO IN SAME STRAND AS markers34
the.snps<-colnames(all.samples)
for(i in 1:length(the.snps)){
one.col<-gsub("V",markers34[the.snps[i],"a1"],all.samples[,the.snps[i]])
one.col<-gsub("F",markers34[the.snps[i],"a2"],one.col)
all.samples[,the.snps[i]]<-one.col
}
############# all.samples  Now in allele format on same strand as markers34
############ now check stand using markers as above must have an A


##### THIS ONLY WORKS CAUSE i HAVE ILLUMINA MARKERS OF SNP TYPE 1
to.flip<-unlist(apply(as.matrix(markers34[,"context"]),1,function(x) length(grep("A",unlist(strsplit(x,split="/"))) ) ==0   ))
sum(rownames(markers34)!=colnames(all.samples)) # MUST BE ZERO check in same order
names(to.flip)<-rownames(markers34)
the.snps<-colnames(all.samples)
for(i in 1:length(the.snps)){
  if(to.flip[the.snps[i]]){
one.col<-gsub("T","A",all.samples[,the.snps[i]])
one.col<-gsub("C","Cfirst",one.col)
one.col<-gsub("G","C",one.col)
one.col<-gsub("Cfirst","G",one.col)

all.samples[,the.snps[i]]<-one.col
}
}
#####################
################ IN TOP STRAND

############### NOW count alleles
all.samples[1:5,]

the.snps<-colnames(all.samples)
for(i in 1:length(the.snps)){
snp.counts<-unlist(apply(as.matrix(all.samples[,the.snps[i]]),1,function(x) length(grep(hapmap.calls[the.snps[i]],unlist(strsplit(x,split=" ")))) ) )
missing<-unlist(apply(as.matrix(all.samples[,the.snps[i]]),1,function(x) length(unlist(strsplit(x,split=" ")))  ==0  ) ) 
## snp.counts<-unlist(apply(as.matrix(all.samples[,the.snps[i]]),1,function(x) length(grep(hapmap.calls[the.snps[i]],unlist(strsplit(x,split="")))) ) )  ### note missing vales are zero here
## missing<-unlist(apply(as.matrix(all.samples[,the.snps[i]]),1,function(x) length(unlist(strsplit(x,split="")))  ==0  ) )
snp.counts[missing]<-NA
all.samples[,the.snps[i]]<-snp.counts
}
sum(apply(all.samples,1,function(x) sum(x,na.rm=TRUE))==0) ## should not be zero unless ALL SNPs all missing for 1 sample
sum(apply(all.samples,2,function(x) sum(x,na.rm=TRUE))==0) ## should not be zero unless 1 SNP all missing for ALL samples



#############################TESTS ######################
#############################TESTS ######################
#############################TESTS ######################
#############################TESTS ######################
################# locate problem sample or SNP
all.samples[apply(all.samples,1,function(x) sum(x,na.rm=TRUE))==0,] #plate 2 HB0663 no genotypes


############# chect for samples with missing data
rownames(all.samples[apply(all.samples,1,function(x) sum(x,na.rm=TRUE))==0,] )
#all.samples[apply(all.samples,1,function(x) sum(!is.na(x))==0),]  "H-2422930043" "H-2482940881" "H-2482970762"


#all.samples[,apply(all.samples,2,function(x) sum(x,na.rm=TRUE))==0] #plate 4 RS953035 is all one genotype!



##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
################################# END  section to reformat data



####################################### READY FOR TRAINING ######################


########################### Now since loaded the old data set there train.extened are alreday built with
########################## extened CAU and YRI,CHB ethnicities
 train<-train.extended
 train.ann<-train.ann.extended
#########################################

#######################################  Do training
 tapply(train.ann[,"Ethnicity"],train.ann[,"Ethnicity"],length)
##  1  2  3   # CUE Asian African
## 89 96 45
################ 3 ethnic groups to use

 sum(is.na(train.ann[,"Ethnicity"]))
#[1] 0

 
################ for mat 34 market set
#all.samples<-geno34.FHCRC.plate2.snps
keep.snps<-colnames(all.samples)
train<-train[,keep.snps]
##################################

 train<-cbind(train,train.ann[,"Ethnicity"])
 colnames(train)[dim(train)[2]]<-"Ethnicity"
  tapply(train[,"Ethnicity"],train[,"Ethnicity"],length)

#lda.hapmap<-lda(train,as.factor(train.ann[,"Ethnicity"]),CV=FALSE)  # if true cant use predict (does not return a lda object
lda.hapmap<-lda(Ethnicity~.,data=train,CV=FALSE)
forest.hapmap<-randomForest(Ethnicity~.,data=train,mtry=2,importance=TRUE)

z<-scale(train[,1:dim(train)[2]-1],scale=FALSE) %*% lda.hapmap$scaling
#plot(z[,1],pch=as.character(train.ann[,"Ethnicity"]),col=c("blue","green"))

y.col<-as.character(train.ann[,"Ethnicity"])
y.col[y.col==1]<-"blue"
y.col[y.col==2]<-"green"
y.col[y.col==3]<-"brown"
plot(z[,1],z[,2],pch=as.character(train.ann[,"Ethnicity"]),col=y.col,xlab="LDA1",ylab="LDA2",main="LDA HapMap trained")


################### need this for lda
 mypredict.lda <- function(object, newdata) predict(object, newdata = newdata)$class


 ###############################################Misclassification error:  0.02 


######################################################################
################################ DO A NEW MODEL #####################
######################################################################
###############################TEST MODEL ON REAL DATA     all.eth.650.snps<-all.samples   all.eth.650.calls<-all.samples.calls


test.snps<-all.samples

dim(all.samples) #eth 650Y


####################### just for discovery    #SNPs in extened model)

# sort(   apply(test.snps,2,function(x){sum(is.na(x))})  )    # snps
# sort(   apply(test.snps,1,function(x){sum(is.na(x))})  )    # samples

test.good.snps<-test.snps[apply(test.snps,1,function(x){sum(is.na(x))})==0,]
test.NAs.snps<-test.snps[apply(test.snps,1,function(x){sum(is.na(x))})!=0,]


 predictor.lda<-predict(lda.hapmap,test.good.snps,na.action=na.omit)
 predictor.forest<-predict(forest.hapmap,test.good.snps,na.action=na.omit)

 plot(predictor.lda$x[,1],predictor.lda$x[,2],pch=as.character(predictor.lda$class))

 #error.data<-list(name="ErrorRate",name="SNPs")
 error.data<-cbind(errorest(Ethnicity ~ .,data=train,model=randomForest,mtry=2)$error , NA)
 #error.data<-cbind(errorest(Ethnicity ~ .,data=train,model=randomForest,mtry=2)$error , NA)

 colnames(error.data)<-c("ErrorRate","SNPs")

 class.results<-predictor.lda$x

class.results[,2]<-predictor.forest
class.results[,1]<-predictor.lda$class
class.results<-cbind(class.results,0,as.numeric(error.data[dim(error.data)[1],"ErrorRate"]),1)
colnames(class.results)<-c("lda","forest","missing","error","posn.in.error.dat")


##################### decided not to do this just recalculate the modle #########
##################### model with fewer SNPs #####################
##################### order by throwinh away less important first   #####################
##################### Some infor about the model    #####################
#important.snps<-sort( forest.hapmap$importance[,4],dec=FALSE)  ## smaller numbers less important
#new.order<-sort(rank.with.missing*important.snps,dec=FALSE)
#rank.with.missing<-missing.snps[labels(important.snps)]
###########################################################

test.NAs.snps.ori<-test.NAs.snps

#### WARNING assume SNPS in same order as in train     rownames(test.NAs.snps)
test.NAs.snps<-test.NAs.snps[,colnames(train)[1:dim(train)[2]-1]]

count<-0
min.snps<-25
num.missing.snps<-0
while(dim(test.NAs.snps)[1]> 0 & (num.missing.snps < min.snps)  ){


 print(paste("Loop counter at ",count,sep=":"))
 count<-count+1
###Calculate / Re-Calculate positions of
missing.snps<-sort(apply(test.NAs.snps,1,function(x){sum(is.na(x))})  )

test.NAs.snps<-test.NAs.snps[labels(missing.snps),]
posns.of.NAs<-  apply(test.NAs.snps,1,function(x){ grep("TRUE",is.na(x))} )
if(class(posns.of.NAs)=="matrix"){
posns.of.NAs<-t(posns.of.NAs)
#dim(posns.of.NAs)= c(dim(test.NAs.snps)[1], (dim(posns.of.NAs)[1]*dim(posns.of.NAs)[2]/dim(test.NAs.snps)[1])  )
#rownames(posns.of.NAs) <- rownames(test.NAs.snps)
tmp<-apply(posns.of.NAs,1,function(x) list(x))
posns.of.NAs<-tmp
}

posns.of.NAs.lengths<-lapply(posns.of.NAs,function(x) length(unlist(x)) )

to.get<-lapply(posns.of.NAs,function(x){length(intersect(unlist(posns.of.NAs[1]),unlist(x)))})
to.get<-(length(unlist(posns.of.NAs[1]))==unlist(to.get)) & ( unlist(posns.of.NAs.lengths)==unlist(to.get))

to.test<- test.NAs.snps[to.get,]
#dim(test.NAs.snps)
#dim(to.test)
snps.to.excluded<-unlist(posns.of.NAs[1])
 
num.missing.snps <- length(snps.to.excluded)
print(paste("Number of missing SNPs ",num.missing.snps,sep=":"))


if( num.missing.snps > min.snps ){
    next
    ## print("More than 1/2 of SNps missing - no prediction made")
    ## class.results.part<-cbind("NA","NA",num.missing.snps ,1,count+1)
    ## predictor.forest.part<-"NA"
    ## predictor.lda.part.class<-"NA"
}else{
 
lda.hapmap.part<-lda(Ethnicity~.,data=train[,-snps.to.excluded],CV=FALSE)
forest.hapmap.part<-randomForest(Ethnicity~.,data=train[,-snps.to.excluded],mtry=2,importance=TRUE)    # train new model

predictor.lda.part<-predict(lda.hapmap.part,to.test,na.action=na.omit)     # run prediction LDA
predictor.forest.part<-predict(forest.hapmap.part,to.test,na.action=na.omit)  # run prediction RandomForest

error.data<-rbind(error.data,c(errorest(Ethnicity ~ .,data=train[,-snps.to.excluded],model=randomForest,mtry=2)$error, list(colnames(test.NAs.snps)[unlist(posns.of.NAs[1])] )   ) ) #estimate error

test.NAs.snps<-test.NAs.snps[!to.get,]
print(paste("Number of remaining samples ",dim(test.NAs.snps)[1],sep=":"))
#test.NAs.snps

class.results.part<-predictor.lda.part$x
class.results.part<-cbind(class.results.part,num.missing.snps,as.numeric(error.data[dim(error.data)[1],"ErrorRate"]),count+1)
predictor.lda.part.class<-predictor.lda.part$class
}

 

colnames(class.results.part)<-c("lda","forest","missing","error","posn.in.error.dat")
 class.results.part[,2]<-predictor.forest.part
class.results.part[,1]<-predictor.lda.part.class

#class.results.part
#dim(class.results.part)
#dim(test.NAs.snps)

class.results<-rbind(class.results,class.results.part)
#dim(class.results)

}

dim(class.results) #is not correct use below
dim(all.samples) ## difference due to low genotypeing fix later


dim(test.snps)

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
####################### FINISHED ethnicity predictions

######################################################################################################
######################################################################################################

######################################################################################################
##############################   clean up and write results:
class.results[1:5,]


class.results[class.results[,"forest"]!=1,]

rownames(class.results[class.results[,"forest"]!=1,])

not.ceu<-class.results[class.results[,"forest"]!=1,]
not.ceu.samples<-rownames(class.results[class.results[,"forest"]!=1,])

not.ceu<-not.ceu[,c("forest","error")]
colnames(not.ceu)<-c("ethnicity_Prediction","False_Discovery_rate_for_model")


not.ceu[not.ceu[,1]==2,1]<-"Asian"
not.ceu[not.ceu[,1]==3,1]<-"African"
not.ceu[,2]<-round(as.numeric(not.ceu[,2]),3)

tapply(not.ceu[,"ethnicity_Prediction"],not.ceu[,"ethnicity_Prediction"],length)
100*dim(not.ceu)[1]/dim(all.samples)[1] # not ceu rate


################## identify samples with less than min.snps where no prediction was made
if (length(setdiff(rownames(test.snps),rownames(class.results)))>0 ){
low.genotypes<-cbind(setdiff(rownames(test.snps),rownames(class.results)),"low genotyping",NA)
rownames(low.genotypes)<-low.genotypes[,1]
low.genotypes<-low.genotypes[,-1]
colnames(low.genotypes)<-c("ethnicity_Prediction","False_Discovery_rate_for_model")

not.ceu<-rbind(not.ceu,low.genotypes)}

###################################################################################

dim(class.results) #is not correct use below
dim(all.samples)
not.ceu

################# write data that are NOT CEU
write.table(not.ceu,"DATA NOT CEU table.txt",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t")

class.results[1:5,]
class.report<-class.results[,1:4]
class.report[class.report[,1]==1,1]<-"CEU"
class.report[class.report[,1]==2,1]<-"Asian"
class.report[class.report[,1]==3,1]<-"African"

class.report[class.report[,2]==1,2]<-"CEU"
class.report[class.report[,2]==2,2]<-"Asian"
class.report[class.report[,2]==3,2]<-"African"

colnames(class.report)<-c("ethnicity_Prediction_LDA","ethnicity_Prediction_Ran.Forest","SNPs.missing","False_Discovery_rate_for_model")
class.report[1:5,]

write.table(class.report,"Prediction_report.txt",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")



######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
