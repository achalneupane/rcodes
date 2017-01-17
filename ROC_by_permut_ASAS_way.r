



source("/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/functions.R")
library(akima)
library(ROCR)
library(nlme)

## load clinical data

## table(data[,"PHEN"])

## table(data[,"PHEN"])

##    0    1 
## 4294  137



## table(data$PHEN)  ## must be 0 and 1 not (1 and 2 else routines below give problems)
## data[1:5,]
## ##  0  1 
## 63 59
model<-"test"

data<-test[,c("pheno",snps)]
dim(data)

## Split data for CV nchooseK
Nsplits <- 300


phen.col="pheno"
data[1:5,]


## Nsplits=dim(data)[1]
## (1-1/Nsplits)*Nsplits
## data.list <- get_N_data_splits(data=data,trai.p=(1-1/Nsplits),N=Nsplits,phen.col='PHEN')


data.list <- get_N_data_splits(data=data,trai.p=.5,N=Nsplits,phen.col=phen.col) # using 0.7 as a slip made no difference
## highter numbers inclde more samples in training less in test

## Model to try
names(data.list[[1]])
#data.list[[1]]$train[,phen.col]
sum(data.list[[1]]$train[,phen.col])
sum(data.list[[1]]$test[,phen.col])

####
formula <- form


############################################ RUN the formulae
#PHEN ~ something + (1|ethnicity) + (1|site)  
predictions0 <- list()   
labels <- list()
problem.boot<-{}
top.performers<-{}

i<-1

for(i in 1:Nsplits) {
  train <- data.list[[i]]$train
  test <- data.list[[i]]$test
  ## table(train$PHEN)
  ## table(test$PHEN)
  ## table(train$imm_16_28236248_A)

  fit1 <- glm(formula,data=train, family=binomial(logit),control=glm.control(epsilon=1e-16,maxit=1000))
 #fit1 <- lm(formula,data=train) # gives the same results as above
 #  fit1 <- glmmPQL(formula , data=train, random= (~1 |ethnicity) , family=binomial(logit),control=glm.control(epsilon=1e-16,maxit=1000)) 
#  fit1 <- glmmPQL(formula , data=train, random= list((~1 | ethnicity),(~1 | HLAB27)) , family=binomial(logit),control=glm.control(epsilon=1e-16,maxit=1000))
 ##  print("---------------------------------------")
 ## print(i)
 ## print(summary(fit1))
 ##   print("---------------------------------------")
  
  pred <- predict(fit1,newdata=test,type='response')
  predictions0[[i]] <- pred

  
  top.performers<-c(top.performers,names(sort(summary(fit1)$coefficients[,"Pr(>|z|)"],decreasing=FALSE))[1:5])
 if(sum(abs(summary(fit1)$coefficients[,"z value"])< 1e-6)){if(length(problem.boot)==0){problem.boot<-i}else{problem.boot<-c(problem.boot,i)}}

  labels[[i]] <-  test[,phen.col]

  
}

length(labels)
length(predictions0)


sort(tapply(top.performers,top.performers,length),decreasing=TRUE)
names(sort(tapply(top.performers,top.performers,length),decreasing=TRUE))
gsub(", "," + ",toString(names(sort(tapply(top.performers,top.performers,length),decreasing=TRUE))))

## > sort(tapply(top.performers,top.performers,length),decreasing=TRUE)
##   (Intercept)    HLA_C_0303   HLA_DRB1_13 HLA_DQA1_0102 HLA_DRB1_1301 
##           300           206           193           113            74 
##   HLA_DQB1_03    HLA_B_1501 HLA_DQB1_0603   HLA_DRB1_15 HLA_DRB1_0401 
##            72            65            63            59            53 
## HLA_DQB1_0301   HLA_DRB1_04 HLA_DQB1_0602 HLA_DQA1_0103 HLA_DRB1_1501 
##            46            35            34            31            28 
##      HLA_C_03      HLA_A_03    HLA_A_0301      HLA_B_44    HLA_B_4402 
##            23            20            20            19            18 
##    HLA_C_0702      HLA_B_15 
##            16            12 
## class(error)
length((predictions0[[1]]))
length((predictions0[[1]]))
## Get plots
length(predictions0)
length(problem.boot)
dim(data.list[[1]]$train)

if(length(problem.boot)>0){ ### these indicated instances where one covaraiate not polymorphic/ variable for a givem predict/train combination they do not affect the requst
predictions0<-predictions0[ -1*problem.boot ]
labels<-labels[ -1*problem.boot ]
}
predict.gene.definite<-predictions0
length(predictions0)




formula
model


pred0 <- prediction(predictions0,labels)
perf0 <- performance(pred0, measure = "tpr", x.measure = "fpr")

length(pred0)
plot(perf0,avg="threshold",lwd=2,col=1,colorize=F)
plot(perf0,lwd=2,col=1,colorize=F,add=T)

## # it allows two different plots in the same frame
## par(mfrow = c(1,2))
## # plot a ROC curve for a single prediction run
## # and color the curve according to cutoff.
## library(ROCR)
## data(ROCR.simple)
## pred <- prediction(ROCR.simple$predictions, ROCR.simple$labels)
## perf <- performance(pred,"tpr", "fpr")
## plot(perf,colorize = TRUE)
## # plot a ROC curve for a single prediction run
## # with CI by bootstrapping and fitted curve
## library(verification)
## roc.plot(ROCR.simple$labels,ROCR.simple$predictions, xlab = "False positive rate",
## ylab = "True positive rate", main = NULL, CI = T, n.boot = 100, plot = "both", binormal = TRUE)
## (auc <- as.numeric(performance(pred, measure = "auc", x.measure = "cutoff")@y.values))


par(mar=c(4.5,4.5,2.5,2),mgp=c(3,1,0)) #c(bottom, left, top, right)
perf0.acc <- performance(pred0, "acc")
plot(perf0.acc, avg= "vertical", spread.estimate="boxplot", lwd=3,col='blue', show.spread.at= seq(0.1, 0.9, by=0.1),main= "",cex=1.5,cex.lab=1.75,font=2,xaxis.cex.axis=2,yaxis.cex.axis=2,xlab.cex.axis=2,colorkey.relwidth=0.5)
                                        # main= "Accuracy across the range of possible cutoffs"
savePlot(paste(model,"Accuracy_vs_Cutoff.png",sep=""),type="png")

plot(perf0, avg= "threshold", colorize=T, lwd= 5, main= "",colorkey.relwidth=0.5,colorkey.pos="right",cex.lab=1.75,font=2,xaxis.cex.axis=2,yaxis.cex.axis=2,xlab.cex.axis=2)
#print.cutoffs.at=c(0.4,0.5,0.6)
auc<-performance(pred0,'auc')
auc<- unlist(auc@y.values)
the.auc<-mean(auc,trim=0.1)
the.auc
sd.auc<-sd(auc)
sd.auc


x<-0.6 # x<-0.8

sen<-performance(pred0,'sens')
datax<-sen@x.values
datay<-sen@y.values

##
temp<-lapply(datax,length)
tapply(unlist(temp),unlist(temp),length)
test<-tapply(unlist(temp),unlist(temp),length)
to.keep<-max(test)
to.keep<-as.numeric(names(test)[test==max(test)])
to.keep
if(length(to.keep)>1){to.keep<-to.keep[1]} # in case multiple ones of the same length
datax<-datax[unlist(temp) ==to.keep]
datay<-datay[unlist(temp) ==to.keep]
length(datax)




av.x<-apply(do.call("rbind",datax),2,function(x) mean(x,trim=0.5))
av.y<-apply(do.call("rbind",datay),2,function(x) mean(x,trim=0.5)) # need trim to remove outlyers
sd.y<-apply(do.call("rbind",datay),2,function(x) sd(x))

trim<-1
the.length<-length(av.x)
av.sens<-aspline(as.numeric(av.x)[(1+trim):(the.length-trim)], as.numeric(av.y)[(1+trim):(the.length-trim)],x,method="improved",degree=3)$y
sd.sens<-aspline(as.numeric(av.x)[(1+trim):(the.length-trim)], as.numeric(sd.y)[(1+trim):(the.length-trim)],x,method="improved",degree=3)$y

## plot(av.x,av.y)
## plot(av.x,sd.y)
## length(as.numeric(av.x)[(1+trim):(the.length-trim)])
## length(as.numeric(av.y)[(1+trim):(the.length-trim)])


spec<-performance(pred0,'spec')
datax<-spec@x.values
datay<-spec@y.values

temp<-lapply(datax,length)
tapply(unlist(temp),unlist(temp),length)
test<-tapply(unlist(temp),unlist(temp),length)
to.keep<-max(test)
to.keep<-as.numeric(names(test)[test==max(test)])
if(length(to.keep)>1){to.keep<-to.keep[1]}
datax<-datax[unlist(temp) ==to.keep]
datay<-datay[unlist(temp) ==to.keep]
length(datax)


av.x<-apply(do.call("rbind",datax),2,function(x) mean(x,trim=0.5))
av.y<-apply(do.call("rbind",datay),2,function(x) mean(x,trim=0.5)) # need trim to remove outlyers
sd.y<-apply(do.call("rbind",datay),2,function(x) sd(x))
trim<-1
the.length<-length(av.x)
av.spec<-aspline(as.numeric(av.x)[(1+trim):(the.length-trim)], as.numeric(av.y)[(1+trim):(the.length-trim)],x,method="improved",degree=3)$y
sd.spec<-aspline(as.numeric(av.x)[(1+trim):(the.length-trim)], as.numeric(sd.y)[(1+trim):(the.length-trim)],x,method="improved",degree=3)$y




## plot(av.x,av.y,col="blue",cex=3)
## auc<- unlist(auc@y.values)
## the.auc<-mean(auc,trim=0.1)
## the.auc
## sd.auc<-sd(auc)
## sd.auc

leg.txt<-c(paste("Av. AUC        : ",signif(the.auc,digits=2), expression("\u00B1") ,signif(sd.auc,digits=1),sep=" "), # use character map to define
           paste("Av. Sensitivity: ",signif(av.sens,digits=2), expression("\u00B1") ,signif(abs(sd.sens),digits=1),sep=" "),
              paste("Av. Specificity: ",signif(av.spec,digits=2), expression("\u00B1") ,signif(abs(sd.spec),digits=1),sep=" ") )

leg.txt
# leg.txt<-gsub("0.7 ","0.70 ",leg.txt)
# leg.txt<-gsub("0.6 ","0.60 ",leg.txt)
#text(0.1,0.2,expression(paste("Av. AUC: ",eval(val),sep="") %+-% 1),cex=2)
#legend(-0.1,1.0,legend=leg.txt,cex=1.5,bty="n")

legend(0.3,1.0,legend=leg.txt,cex=1.5,bty="n")
#legend(0.2,0.8,legend=leg.txt,cex=1.5,bty="n")
#text(0,0.8,labels=leg.txt,cex=1.5,pos=4)
paste(model,"_ROC.png",sep="")
getwd()
savePlot(paste(model,"_ROC.png",sep=""),type="png")
savePlot(paste(model,"_ROC.tiff",sep=""),type="tiff")




###################################### SAVE a hi resolution epx figure###############
X11(width = 7, height = 7) 
plot(perf0, avg= "threshold", colorize=T, lwd= 5, main= "",colorkey.relwidth=0.5,colorkey.pos="right",cex.lab=1.75,font=2,xaxis.cex.axis=2,yaxis.cex.axis=2,xlab.cex.axis=2)

#legend(-0.1,1.0,legend=leg.txt,cex=1.5,bty="n")

legend(0.3,0.6,legend=leg.txt,cex=1.5,bty="n")

dev.copy2eps(file=paste(model,"_ROC.eps",sep=""))
dev.off()


##################### straight ROC no color
plot(perf0, avg= "threshold", colorize=F, lwd= 5, main= "",cex.lab=1.75,font=2,xaxis.cex.axis=2,yaxis.cex.axis=2,xlab.cex.axis=2)


#legend(-0.1,1.0,legend=leg.txt,cex=1.5,bty="n")
legend(0.3,0.6,legend=leg.txt,cex=1.5,bty="n")

savePlot(paste(model,"_ROC_noColour.png",sep=""),type="png")
savePlot(paste(model,"_ROC_noColour.tiff",sep=""),type="tiff")
dev.off()
########## EPS

X11(width = 7, height = 7) 
plot(perf0, avg= "threshold", colorize=F, lwd= 5, main= "",cex.lab=1.75,font=2,xaxis.cex.axis=2,yaxis.cex.axis=2,xlab.cex.axis=2)

#legend(-0.1,1.0,legend=leg.txt,cex=1.5,bty="n")
legend(0.3,0.6,legend=leg.txt,cex=1.5,bty="n")

dev.copy2eps(file=paste(model,"_ROC_noColour.eps",sep=""))
dev.off()

##########################################################################



