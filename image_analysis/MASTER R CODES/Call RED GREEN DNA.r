 setwd("/media/Bioinform-D/Research/Cellomics/Hugo screen/")


 update.packages(lib="/home/pleo/R_latest/library")

############################################## START REQUIRED FUNCTIONs ###################
############################################## START REQUIRED FUNCTIONs ###################
############################################## START REQUIRED FUNCTIONs ###################
a.model<-function(x,A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static){
              ak.curve<-aspline(static$x, static$y,x,method="improved",degree=3)
              ak.curve$y[!is.finite(ak.curve$y)]<-0.0
              A*dnorm(x,mean=g1.peak.posn,sd=g1.sd, log=FALSE) + B*dnorm(x,mean=g2.peak.posn,sd=g2.sd, log=FALSE) + abs(ak.curve$y*(( 1*pnorm( sqrt(2)*( ((x-g1.peak.posn)/g1.sd.inter)-k1), lower.tail=TRUE, log.p=FALSE)-1) -(1*(pnorm( sqrt(2)*( ((x-g2.peak.posn)/g2.sd.inter)+k2), lower.tail=TRUE, log.p=FALSE))-1 ))) + ak.curve$y*(pnorm( sqrt(2)*( ((x-g2.peak.posn)/g2.sd)-k3), lower.tail=TRUE, log.p=FALSE)) +  -ak.curve$y*(pnorm( sqrt(2)*( ((x-g1.peak.posn)/g1.sd)+k0), lower.tail=TRUE, log.p=FALSE)-1)                                                                               }

 a.model.S<-function(x,A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static){
              ak.curve<-aspline(static$x, static$y,x,method="improved",degree=3)
               ak.curve$y[!is.finite(ak.curve$y)]<-0.0
             abs(ak.curve$y*(( 1*pnorm( sqrt(2)*( ((x-g1.peak.posn)/g1.sd.inter)-k1), lower.tail=TRUE, log.p=FALSE)-1) -(1*(pnorm( sqrt(2)*( ((x-g2.peak.posn)/g2.sd.inter)+k2), lower.tail=TRUE, log.p=FALSE))-1 ))) }

  a.model.gt4<-function(x,A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static){
              ak.curve<-aspline(static$x, static$y,x,method="improved",degree=3)
              ak.curve$y[!is.finite(ak.curve$y)]<-0.0
              ak.curve$y*(pnorm( sqrt(2)*( ((x-g2.peak.posn)/g2.sd)-k3), lower.tail=TRUE, log.p=FALSE))
              }

 a.model.lt2<-function(x,A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static){
              ak.curve<-aspline(static$x, static$y,x,method="improved",degree=3)
               ak.curve$y[!is.finite(ak.curve$y)]<-0.0
               -ak.curve$y*(pnorm( sqrt(2)*( ((x-g1.peak.posn)/g1.sd)+k0), lower.tail=TRUE, log.p=FALSE)-1)   }

  a.model.G1andG2<-function(x,A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static){
                A*dnorm(x,mean=g1.peak.posn,sd=g1.sd, log=FALSE)+B*dnorm(x,mean=g2.peak.posn,sd=g2.sd, log=FALSE)       }

 a.model.true<-function(x,A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static){
                 ak.curve<-aspline(static$x, static$y,x,method="improved",degree=3)
                 ak.curve$y[!is.finite(ak.curve$y)]<-0.0
                 ak.curve$y  }

a.model.G1<-function(x,A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static){
               A*dnorm(x,mean=g1.peak.posn,sd=g1.sd, log=FALSE)              }

a.model.G2<-function(x,A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static){
                B*dnorm(x,mean=g2.peak.posn,sd=g2.sd, log=FALSE)       }

###############################
peaks <- function(series, span = 3, do.pad = TRUE) {
    if((span <- as.integer(span)) %% 2 != 1) stop("'span' must be odd")
    s1 <- 1:1 + (s <- span %/% 2)
    if(span == 1) return(rep.int(TRUE, length(series)))
    z <- embed(series, span)
    v <- apply(z[,s1] > z[, -s1, drop=FALSE], 1, all)
    if(do.pad) {
        pad <- rep.int(FALSE, s)
        c(pad, v, pad)
    } else v
}

peaksign <- function(series, span = 3, do.pad = TRUE)
{
    ## Purpose: return (-1 / 0 / 1) if series[i] is ( trough / "normal" / peak )
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 25 Nov 2005

    if((span <- as.integer(span)) %% 2 != 1 || span == 1)
        stop("'span' must be odd and >= 3")
    s1 <- 1:1 + (s <- span %/% 2)
    z <- embed(series, span)
    d <- z[,s1] - z[, -s1, drop=FALSE]
    ans <- rep.int(0:0, nrow(d))
    ans[apply(d > 0, 1, all)] <- as.integer(1)
    ans[apply(d < 0, 1, all)] <- as.integer(-1)
    if(do.pad) {
        pad <- rep.int(0:0, s)
        c(pad, ans, pad)
    } else ans
}


check.pks <- function(y, span = 3)
    stopifnot(identical(peaks( y, span), peaksign(y, span) ==  1),
              identical(peaks(-y, span), peaksign(y, span) == -1))

for(y in list(1:10, rep(1,10), c(11,2,2,3,4,4,6,6,6))) {
    for(sp in c(3,5,7))
        check.pks(y, span = sp)
    stopifnot(peaksign(y) == 0)
}


############################################## END REQUIRED FUNCTIONs ###################
############################################## END REQUIRED FUNCTIONs ###################
############################################## END REQUIRED FUNCTIONs ###################
#setwd("/media/Bioinform-D/Data/Cellomics/Inhibitor_data")

library(robust)
library(robustbase)
library(nls2)
library(akima)


## plate.numbers<-c(1:15,17:20,25:34)  #20
## file.list=paste("Kiril_Plate",plate.numbers,sep="")
## file.name.list<-file.list

## file.list<-c("plate_1","plate_2","plate_3","plate_4","plate_5","plate_6","plate_7","plate_8","plate_9","plate_10","plate_11","plate_12","plate_13","plate_14","plate_15","plate_16","plate_17","plate_18","plate_21")
## file.name.list<-file.list

file.list<-c("Ca_1","Ca_2","Ca_3")
file.name.list<-file.list

do.trace<-FALSE # tarce for non-linear fit
robust.shaver<-1  # used 0.1 before
max.cells.per.field<-3000
min.green.cells<-100
min.green.threshold<-50
max.g1.posn<-12000
min.g1.posn<-5000
expected.g1.posn<-8000
max.ObjectTotalIntenCh1<-30000
double.exposure<-TRUE
two.color<- TRUE
red.thresh<-50


for (ifile in 1:length(file.list)) {
#  for (ifile in 20:20) {
file<-file.list[ifile]
file.name<-file.name.list[ifile]
print(file.name)
### set file_name and jump here for single file


  
### get the samples in the column names
reads<-100000
if(double.exposure){keep<-c(2,3,5,6,10,17:27)}else{keep<-c(2,3,5,6,10,17:24)}  # columns to keep
#if(double.exposure){keep<-c(3,5,6,10,17:26)}else{keep<-c(3,5,6,10,17:24)} # for kiril

options(show.error.messages = TRUE)
chromo<-try(read.delim(paste(file,".TXT",sep=""),header=T,nrows=1,sep="\t",fill=TRUE))
num.vars<-dim(chromo)[2]
vars.names<-colnames(chromo)[1:dim(chromo)[2]]
vars.names<-sub("TargetActivationV3Cell.","",vars.names)
##########################   dim(chromo)<-c(num.lines,num.vars)
header.lines<-1
num.lines<-1
cells<-{}
################################### read one plate in one go
chromo<-try(scan(paste(file,".TXT",sep=""),what=character(num.vars),skip=header.lines,sep="\t",fill=TRUE))
num.lines<-length(chromo)/(num.vars)
dim(chromo)<-c(num.vars,num.lines)
chromo<-t(chromo)
cells<-chromo[,keep]

###################################to a read in a lrage file
# counter<- -1
# while (num.lines >0  ){
# counter<-counter+1
# counter
# chromo<-try(scan(#paste(file,"TXT",sep="."),what=character(num.vars),skip=(reads*counter)+header.lines,nlines=reads,sep="\t",fill=TRUE))
# num.lines<-length(chromo)/(num.vars) # -1 cause of ContrilCase0
# dim(chromo)<-c(num.vars,num.lines)
# chromo<-t(chromo)
# cells<-rbind(cells,chromo[,keep])
# }   # while rad in one data file

colnames(cells)<-vars.names[keep]
cellcells<-cells[-dim(cells)[1],]

######################### remap for Kiril 2 exposure settings
if(double.exposure){
colnames(cells)[colnames(cells)=="TotalIntenCh3"]<-"TotalIntenCh2b"
colnames(cells)[colnames(cells)=="AvgIntenCh3"]<-"AvgIntenCh2b"
colnames(cells)[colnames(cells)=="VarIntenCh3"]<-"VarIntenCh2b"
colnames(cells)[colnames(cells)=="TotalIntenCh4"]<-"TotalIntenCh3"
colnames(cells)[colnames(cells)=="AvgIntenCh4"]<-"AvgIntenCh3"
colnames(cells)[colnames(cells)=="VarIntenCh4"]<-"VarIntenCh3"
}

if(two.color){
  cells[,"TotalIntenCh3"]<-cells[,"TotalIntenCh2"]
  cells[,"AvgIntenCh3"]<-cells[,"AvgIntenCh2"]
  cells[,"VarIntenCh3"]<-cells[,"VarIntenCh2"]
}

#####################
#loop if have multiple plate in the file
order.in.file<-unique(cells[,"BarCode"])
sizes<-tapply(cells[,"BarCode"],cells[,"BarCode"],length)
sizes<-sizes[order.in.file]
barCodes<-unique(cells[,"BarCode"])
for (iBarCodes in 1:length(barCodes)){
  print(barCodes[iBarCodes])
the.cells<-cells[cells[,"BarCode"]==barCodes[iBarCodes],]
 the.cells<-the.cells[the.cells[,"Row"]!="",] #strip out crap
 file.plate.name<-paste(file.name,barCodes[iBarCodes],sep="_")


## chk.num.plates<-length(unique(barCodes))
## print(paste("FILE: ",file.name,"  BARCODE:",the.cells[1,"BarCode"],"  UNIQUE:",chk.num.plates,sep=""))
## chk.barcode1<-gsub("OCL1030000","Kiril_Plate",the.cells[1,"BarCode"])
## chk.barcode2<-gsub("OCL103000","Kiril_Plate",the.cells[1,"BarCode"])

## if(chk.barcode1 !=file.name){
##   if(chk.barcode2 !=file.name){print (paste("WARNING","ERROR","BARCODE MISMATCH",sep=" "))}}

number.cols<-c(2:dim(the.cells)[2])
cells.num<-as.numeric(the.cells[,number.cols])
dim(cells.num)<-dim(the.cells[,number.cols] )
colnames(cells.num)<-colnames(the.cells)[number.cols]
cells.num<-as.data.frame(cells.num)
rm(the.cells)
dim(cells.num)

rows<-tapply(cells.num[,"Row"],cells.num[,"Row"],length)
cols<-tapply(cells.num[,"Col"],cells.num[,"Col"],length)


# dim(green.c)<-c(length(rows),length(cols))
# dim(red.c)<-c(length(rows),length(cols))
# dim(notGreen.c)<-c(length(rows),length(cols))
# dim(notRed.c)<-c(length(rows),length(cols))
col.index<-c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
red.c<-c(rep.int(0,length(rows)*length(cols)))
dim(red.c)<-c(length(rows),length(cols))
rownames(red.c)<-col.index[c(as.integer(names(rows)))+1] # +1 cause index starts at zero
colnames(red.c)<-as.character(c(as.integer(names(cols)))+1   )
green.c<-red.c
notGreen.c<-red.c
notRed.c <-red.c
redAndGreen <-red.c
redAndGreenLow <-red.c
redAndGreenMid <-red.c
redAndGreenHigh <-red.c
redLow <-red.c
redMid <-red.c
redHigh <-red.c
greenLow <-red.c
greenMid <-red.c
greenHigh <-red.c
notRedAndGreen <-red.c
redAndNotGreen  <-red.c
all.c <-red.c
all.using <-red.c
big.c <-red.c

data.store<-c(rep(NA,length(rows)*length(cols)))
dim(data.store)<-c(length(rows),length(cols))
rownames(data.store)<-col.index[c(as.integer(names(rows)))+1] # +1 cause index starts at zero
colnames(data.store)<-as.character(c(as.integer(names(cols)))+1   )
lm.red.slope<-data.store
lm.green.slope<-data.store
lm.red.inter<-data.store
lm.green.inter<-data.store
lm.red.coverg  <-data.store
lm.green.coverg  <-data.store
cells.P.field<-data.store
R.P.field  <-data.store
G.P.field  <-data.store

 DNA.G1<-red.c
 DNA.G2<-red.c
 DNA.G1andG2<-list(adenG=red.c,adenNG=red.c,aden=red.c)
 DNA.aboveG1<-list(adenG=red.c,adenNG=red.c,aden=red.c)
 DNA.gt4 <-list(adenG=red.c,adenNG=red.c,aden=red.c)
 DNA.lt2<-list(adenG=red.c,adenNG=red.c,aden=red.c)
 DNA.S<-list(adenG=red.c,adenNG=red.c,aden=red.c)
 DNA.fit.success<-list(adenG=red.c,adenNG=red.c,aden=red.c)
 DNA.A<-list(adenG=red.c,adenNG=red.c,aden=red.c)
 DNA.B<-list(adenG=red.c,adenNG=red.c,aden=red.c)

############RUN THE LOOP

for(i in 1:length(rows)) {
for(j in 1:length(cols)) {
# print(i)
# print(j)
 row.label<-as.numeric(names(rows[i]))    #0 to n label
 col.label<-as.numeric(names(cols[j]))    #0 to n label

 row.label.letter<- col.index[c(as.integer(names(rows[i])))+1]
 col.label.number<- as.character(c(as.integer(names(cols[j])))+1 )
#####load data


test<-cells.num[ (cells.num[,"Row"]==row.label & cells.num[,"Col"]==col.label),]
 
if(dim(test)[1]>200){
 ############# remove high density fields OR FIELDS IN GENERAL
 num.per.field<-tapply(test[,"FieldIndex"],test[,"FieldIndex"],length)
 fields.to.remove<-names(num.per.field[num.per.field> max.cells.per.field])
 posns<-unlist(apply(as.matrix(fields.to.remove),1, function(x) grep(paste("^",x,"$",sep=""),test[,"FieldIndex"])))
 if(length(posns)>0){test<-test[-posns,]}
 ################################

 ### assumes it is symmetric
all.c[i,j]<-dim(test)[1]
#### check if have enough cells


#so have index with respect to wells in HCSView
  rownames(test)<-c(1:dim(test)[1])
big<- (test[,"ObjectTotalIntenCh1"] > max.ObjectTotalIntenCh1 )
big.c[i,j]<-sum(big)
all.using[i,j]<-sum(!big)

test<-test[!big,] #exclude big i.e. gt 4N cells from analysis


xaxis<-"ObjectTotalIntenCh1"
yaxis<-"TotalIntenCh3"
 
the.model <-ltsReg(TotalIntenCh3~ObjectTotalIntenCh1,data=test) # need mcd to get weights but not as good at fitting
#resid.quant<-quantile(the.model$residuals,c(0.90))
#red<- ( (test[,yaxis]>= (coef(the.model)[2]* test[,xaxis]+ coef(the.model)[1]+resid.quant))  )
red<-the.model$lts.wt==0 &  residuals(the.model)>0 & test[,"AvgIntenCh3" ]>= red.thresh
red.c[i,j]<-sum(red)
names(red)<-rownames(test)

##  abline( b=coef(the.model)[2] ,a=(coef(the.model)[1]) ,lty=10,col="green",lwd=2)
## abline( b=coef(the.model)[2] ,a=(coef(the.model)[1]+resid.quant) ,lty=10,col="red",lwd=2) # a intercept b is slope

the.model<-ltsReg(TotalIntenCh3~ObjectTotalIntenCh1,data=test,subset=c(1:dim(test)[1])[!red])
red.trimmed<- the.model$lts.wt==0 &  residuals(the.model)>0 & test[!red,"AvgIntenCh3" ]>= red.thresh
red[rownames(test[!red,])[red.trimmed]]<-TRUE
red.c[i,j]<-sum(red)

##  the.model<-lmrob(TotalIntenCh3~ ObjectTotalIntenCh1,data=test)
## #the.model<-lmrob(TotalIntenCh3~ I(0.5454545*ObjectTotalIntenCh1-125.4) ,data=test)
## #print(as.character(the.model$converged))
## lm.red.slope[i,j]<-coef(the.model)[2]
## lm.red.inter[i,j]<-coef(the.model)[1]
## lm.red.coverg[i,j]<-the.model$converged
## #red<- (abs(weights(the.model)) < 0.1/length(weights(the.model))) &  residuals(the.model)>0 & test[,"AvgIntenCh3"] >= red.thresh
## red<- (abs(weights(the.model)) <  robust.shaver/length(weights(the.model))) &  residuals(the.model)>0 & test[,"AvgIntenCh3"] >= red.thresh
## red.c[i,j]<-sum(red)
## names(red)<-rownames(test)

 the.model<-lmrob(TotalIntenCh3~ObjectTotalIntenCh1,data=test,subset=c(1:dim(test)[1])[!red])
#print(as.character(the.model$converged))
lm.red.slope[i,j]<-coef(the.model)[2]
lm.red.inter[i,j]<-coef(the.model)[1]
lm.red.coverg[i,j]<-the.model$converged
red.trimmed<- (abs(weights(the.model)) < robust.shaver/length(weights(the.model))) &  residuals(the.model)>0 & test[!red,"AvgIntenCh3"] >= red.thresh
red[rownames(test[!red,])[red.trimmed]]<-TRUE
red.c[i,j]<-sum(red)

#######PLOTS
xaxis<-"ObjectTotalIntenCh1"
yaxis<-"TotalIntenCh3"
 #c(bottom, left, top, right)
jpeg( paste(paste(file.plate.name,row.label.letter,col.label.number,sep="_"),"jpeg",sep="."), width = 1280, height = 1024, units = "px",quality = 100 )
par(mar=c(3,5.5,2.5,5.1),mgp=c(2,1,0))
layout(matrix(c(1,2,3),3,byrow=TRUE),heights=c(1,1,1))
#layout.show(nf)
plot((test[,xaxis]),(test[,yaxis]),pch=20,cex=0.1,xlab="Total Intensity of DAPI",ylab="EDU Intensity",main=paste(file.name," Well:",row.label.letter,col.label.number,sep=" "),font.lab=2,cex.lab=1.5,cex.main=2.0)
abline(coef(the.model),lty=10,col="red",lwd=2)
points(test[red,xaxis],test[red,yaxis],col="red",pch=19,cex=0.4)
#dev.off()
#######$$$$ PLOTS
green_thresh<-min.green.threshold
xaxis<-"AvgIntenCh2"
yaxis<-"VarIntenCh2"
#backgroud OR stauration
well.rows.sub<-c(1:dim(test)[1])[ test[,xaxis]>green_thresh & test[,xaxis]<3500 &  test[,yaxis]>0 ]
  
if(length(well.rows.sub) > min.green.cells) {
#print(length(well.rows.sub))
test[,xaxis]<-log2(test[,xaxis])
test[,yaxis]<-log2(test[,yaxis])


the.model <-ltsreg(VarIntenCh2~AvgIntenCh2+I(AvgIntenCh2^2),data=test,subset=well.rows.sub )
resid.quant<-quantile(the.model$residuals,c(0.85))
mostly.green<- ( (test[,yaxis]<= (coef(the.model)[3]* test[,xaxis]^2+coef(the.model)[2]* test[,xaxis]+ coef(the.model)[1]+resid.quant)) & is.finite(test[,xaxis])  &  is.finite(test[,yaxis]) & test[,xaxis]>=log2(green_thresh) )
the.model2 <-lmrob(VarIntenCh2~AvgIntenCh2+I(AvgIntenCh2^2),data=test,subset=c(1:dim(test)[1])[mostly.green] )


#print(as.character(the.model2$converged))
lm.green.slope[i,j]<-toString(signif(coef(the.model),3))
lm.green.coverg[i,j]<-the.model2$converged

green<-(test[,yaxis]<=  ((coef(the.model)[3]+1*sqrt(diag(the.model2$cov))[3])* test[,xaxis]^2  +  (coef(the.model)[2]+0.15*sqrt(diag(the.model2$cov))[2])* test[,xaxis]+ coef(the.model)[1]+ 1*sqrt(diag(the.model2$cov))[1]    ) & is.finite(test[,xaxis])  &  is.finite(test[,yaxis]) & test[,xaxis]>=log2(green_thresh)  )
if( sum(mostly.green) < sum(green) ){green<-mostly.green} ##flipped to < for kiril was > for joseph

############## select green cells in ranges of a third
ghist<-hist(test[green,xaxis],breaks=50,plot=FALSE)
xvals<-ghist$breaks[2:length(ghist$breaks)]
yvals<- cumsum(ghist$counts)
cuts<-quantile(yvals,c(0.333,0.666))
#green.cut.low<- max(xvals[yvals<=cuts[1]])
#green.cut.high<- min(xvals[yvals>=cuts[2]])
green.cut.low<- max(xvals[yvals<= yvals[length(yvals)]/4  ])
green.cut.mid<- max(xvals[yvals<= 2*yvals[length(yvals)]/4  ])
green.cut.high<- min(xvals[yvals>= 3*yvals[length(yvals)]/4 ])

green.c[i,j]<-sum(green)

#######PLOTS
#jpeg( paste(paste(file.name,"GREEN",row.label.letter,col.label.number,sep="_"),"jpeg",sep=".") )
par(mar=c(3,5.5,1.1,5.1),mgp=c(2,1,0)) #c(bottom, left, top, right)
plot((test[,xaxis]),(test[,yaxis]),pch=20,cex=0.1,xlab="log2( Average Green Signal )",ylab="log2( Varience Green Signal )",main="",font.lab=2,cex.lab=1.5) 

points(test[green,xaxis],test[green,yaxis],col="green",pch=20,cex=0.1)
points(test[red,xaxis],test[red,yaxis],col="red",pch=20,cex=0.1)
order.in<-order(test[well.rows.sub,xaxis])
lines(test[well.rows.sub[order.in],xaxis],the.model$fit[order.in],col="magenta")
#dev.off()
########PLOTS

}else{green<-rep(FALSE,dim(test)[1])
      green.c[i,j]<-0
      #jpeg( paste(paste(file.name,"GREEN",row.label.letter,col.label.number,sep="_"),"jpeg",sep=".") )
      par(mar=c(3,5.5,1.1,5.1),mgp=c(2.5,1,0))
      plot(log2(test[,xaxis]),log2(test[,yaxis]),pch=20,cex=0.1,xlab="log2( Average Green Signal )",ylab="log2( Varience Green Signal )",font.lab=2,cex.lab=1.5 )
     # dev.off()
      } #less than 20 green objects detected

######################## do DNA histgrams #############################
  xaxis<-"ObjectTotalIntenCh1"
  aden<-density(test[,xaxis])
  adenNG<-density(test[!green,xaxis])
  if(sum(green)<min.green.cells){adenG<-adenNG}else{adenG<-density(test[green,xaxis])}
  par(mar=c(3.5,5.5,1.1,5.1),mgp=c(2.5,1,0)) #c(bottom, left, top, right)
 the.max.range<-max(aden$y,adenG$y,adenNG$y)
  plot(adenNG,lwd=2,col="grey50",main="",font.lab=2,cex.lab=1.5,ylim=c(0,the.max.range))
  lines(adenG$x,adenG$y,col="green",lwd=2)
  lines(aden$x,aden$y,col="black",lwd=2)
 

  a.peak<-peaks(aden$y,span=55) 
  points(aden$x[a.peak],aden$y[a.peak],pch=23,col="red") ### found peak
 potential.peaks<-aden$y[a.peak]
  highest.place<-order(potential.peaks,decreasing=TRUE)[1]  ## incase find g2 peak first
  g1.peak.posn<-aden$x[a.peak][highest.place]
  g1.peak.height<-aden$y[a.peak][highest.place]
  g1.peak.place<-sum(aden$x <= g1.peak.posn)
  if(g1.peak.posn>max.g1.posn | g1.peak.posn< min.g1.posn){
     highest.place<-order(abs(potential.peaks-expected.g1.posn),decreasing=FALSE)[1]
      g1.peak.posn<-aden$x[a.peak][highest.place]
      g1.peak.height<-aden$y[a.peak][highest.place]
      g1.peak.place<-sum(aden$x <= g1.peak.posn)
      if(g1.peak.posn>max.g1.posn | g1.peak.posn< min.g1.posn){ # still not a good start
        g1.peak.place<-sum(aden$x <= expected.g1.posn)
        g1.peak.posn<-aden$x[g1.peak.place]
        g1.peak.height<-aden$y[g1.peak.place]}
                                 }

 
 # aden$x[g1.peak.place]
  g1.sd<-g1.peak.posn/3.5 # initial guess
  
  A<-g1.peak.height/dnorm(g1.peak.posn,mean=g1.peak.posn,sd=g1.sd, log=FALSE) # initial guess
  g1.region<-aden$x<=(g1.peak.posn+0.15*g1.peak.posn)
  to.fit<-data.frame(y=aden$y[g1.region], x=aden$x[g1.region])
  the.fit<-nls(y~A*dnorm(x,mean=g1.peak.posn,sd=g1.sd, log=FALSE),
               data=to.fit,
               start=list(A=A,g1.peak.posn=g1.peak.posn,g1.sd=g1.sd),
               ,lower=c(A/10, min.g1.posn, g1.sd/5)
               ,upper=c(10, max.g1.posn, g1.sd*5),
                ,trace=do.trace
               ,algorithm="port"
               ,control=list(maxiter=1000, minFactor=1/4048,tol=1e-4, warnOnly = TRUE))
  the.coefs<-coef(the.fit)
  A<-as.numeric(the.coefs["A"])
  g1.peak.posn<-as.numeric(the.coefs["g1.peak.posn"])
  g1.sd<-as.numeric(the.coefs["g1.sd"])
  g1.peak.height<-A*dnorm( g1.peak.posn,mean=g1.peak.posn,sd=g1.sd, log=FALSE)
   points(g1.peak.posn,g1.peak.height,pch=23,col="blue")

   second.highest.place<-order(potential.peaks,decreasing=TRUE)[2]
  g2.peak.posn<-aden$x[a.peak][second.highest.place]
  g2.peak.height<-aden$y[a.peak][second.highest.place]
  g2.peak.place<-sum(aden$x < g2.peak.posn)


  #### for fit range of allowed vales is 2.35 to 1.9 * gi.peak.posn
  if(g2.peak.posn>=2.35*g1.peak.posn |  g2.peak.posn <= 1.9*g1.peak.posn | is.na(aden$x[a.peak][second.highest.place]) ){  ## can't get a good peak
    g2.peak.place<-sum(aden$x < 2.1*g1.peak.posn)
    g2.peak.posn<-aden$x[g2.peak.place]
    g2.peak.height<-aden$y[g2.peak.place]
                                   }
    points(g2.peak.posn,g2.peak.height,pch=23,col="blue")
  g2.sd<-1*g1.sd # initial guess
  
  B<-g2.peak.height/dnorm(g2.peak.posn,mean=g2.peak.posn,sd=g2.sd, log=FALSE) # initial guess
  k0<-2
  k1<-0.7
  k2<-0.7
  k3<-0.7
  g1.sd.inter<-g1.sd
  g2.sd.inter<-g2.sd
  
 ##  mid.place<-g1.peak.place+floor(g2.peak.place-g1.peak.place)/2
##   mid.height<-aden$y[mid.place]
##   mid.posn<-aden$x[mid.place]
 static<-data.frame(x=aden$x,y=aden$y)

to.fit<-data.frame(y=aden$y, x=aden$x)
the.fit<-nls(y~(a.model(x,A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static) )
          ,data=to.fit
          ,start=list(A=A,B=B,g1.peak.posn=g1.peak.posn,g1.sd=g1.sd,g2.peak.posn=g2.peak.posn,g2.sd=g2.sd,k0=k0,k1=k1,k2=k2,k3=k3,g1.sd.inter=g1.sd.inter, g2.sd.inter=g2.sd.inter)
          ,lower=c(A-0.15*A, B-0.3*B, g1.peak.posn-0.5*g1.sd, g1.sd/2, g1.peak.posn*1.9,  g2.sd/2,0,  0 , 0, 0,g1.sd/5,g2.sd/5)
          ,upper=c(A+0.15*A, B+0.5*B, g1.peak.posn+0.5*g1.sd, g1.sd*2, g1.peak.posn*2.35,  g2.sd*2,10, 3,  3, 3,g1.sd*5,  g2.sd*5)
          ,trace=do.trace
          ,algorithm="port"
          ,control=list(maxiter=1000, minFactor=1/2048,tol=1e-4, warnOnly = TRUE) )


the.coefs<-coef(the.fit)
  
  
## the.coefs<-c(7.807640e-01,  3.113320e-01 , 5.947625e-06 , 1.553233e+04 , 4.024419e+03,  3.319354e+04 , 5.603872e+03 , 8.992209e+03, -7.230066e+03)
## the.coefs<-c( 0.720911 ,0.173168, 2.86649e-05 , 15621.2 , 4154.37 , 34761.3 , 4168.79, 0.000169886 ,3.96350e-05 )
## names(the.coefs)<-c("A","B","C","g1.peak.posn","g1.sd","g2.peak.posn","g2.sd","g1.sd.inter","g2.sd.inter","k1","k2")
  
  A<-the.coefs["A"]
  B<-the.coefs["B"]
  g1.peak.posn<-the.coefs["g1.peak.posn"]
  g1.sd<-the.coefs["g1.sd"]
  g2.peak.posn<-the.coefs["g2.peak.posn"]
  g2.sd<-the.coefs["g2.sd"]
  g1.sd.inter<-the.coefs["g1.sd.inter"]
  g2.sd.inter<-the.coefs["g2.sd.inter"]
  k0<-the.coefs["k0"]
  k1<-the.coefs["k1"]
  k2<-the.coefs["k2"]
  k3<-the.coefs["k3"]
DNA.G1[i,j]<- g1.peak.posn
DNA.G2[i,j]<- g2.peak.posn

 curve(A*dnorm(x,mean=g1.peak.posn,sd=g1.sd, log=FALSE),add=TRUE, col="red",lwd=2,lty="dashed")
 curve(B*dnorm(x,mean=g2.peak.posn,sd=g2.sd, log=FALSE),add=TRUE, col="orange3",lwd=2,lty="dashed")
 curve( a.model.S(x,A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static), add=TRUE, col="purple",lwd=2,lty="dashed")
 curve( a.model.gt4(x,A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static), add=TRUE, col="turquoise4",lwd=2,lty="dashed")
   curve( a.model.lt2(x,A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static), add=TRUE, col="turquoise",lwd=2,lty="dashed")
 curve( a.model(x,A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static), add=TRUE, col="blue",lwd=2,lty="dashed")

profiles<-c("adenG","adenNG","aden") # densities want to analyse
the.colors<-c("green","grey34")
for(iprofile in 1:length(profiles)){
  
   static<-data.frame(x=eval(as.name(profiles[iprofile]))$x,y=eval(as.name(profiles[iprofile]))$y)
   
    s.int<-integrate(a.model.S,lower=0,upper=max(static$x),subdivisions=length(static$x),A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static)

 gt4.int<-integrate(a.model.gt4,lower=0,upper=max(static$x),subdivisions=length(static$x),A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static)
    
 lt2.int<-integrate(a.model.lt2,lower=0,upper=max(static$x),subdivisions=length(static$x),A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static)
  
  g1ANDg2<-a.model.true(static$x,A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static)-a.model.S(static$x,A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static)-a.model.gt4(static$x,A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static)-a.model.lt2(static$x,A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static)
  
#lines(static$x, g1ANDg2, col="pink",lwd=2)
  
new.to.fit<-data.frame(y=g1ANDg2, x=static$x) ### fit just A and B to the new data 
new.the.fit<-nls(y~(a.model.G1andG2(x,A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static) )
          ,data=new.to.fit
          ,start=list(A=A,B=B)
          ,lower=c(A/10, B/10)
          ,upper=c(1.1, 1.1)
          ,trace=do.trace
          ,algorithm="port"
          ,control=list(maxiter=1000, minFactor=1/2048,tol=1e-4, warnOnly = TRUE) )

  new.the.coefs<-coef(new.the.fit)
  new.A<-as.numeric(new.the.coefs["A"])
  new.B<-as.numeric(new.the.coefs["B"])
   #### failure can offure if outside paramneter range very rare get A-1.01 sometimes for example
   if(is.na(new.A)){new.A<-A}
    if(is.na(new.B)){new.B<-B}
  
new.g1.int<-integrate(a.model.G1,lower=0,upper=max(static$x),subdivisions=length(aden$x),A=new.A,B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static)


new.g2.int<-integrate(a.model.G2,lower=0,upper=max(static$x),subdivisions=length(aden$x),A,B=new.B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static)


 the.fit.success<- new.g1.int$val +  new.g2.int$val + s.int$val + gt4.int$val + lt2.int$val
 curve( a.model(x,new.A,new.B,g1.peak.posn,g2.peak.posn,g1.sd,g2.sd,k0,k1,k2,k3,g1.sd.inter,g2.sd.inter,static), add=TRUE, col=the.colors[iprofile],lwd=2,lty="dotted")

DNA.G1andG2[[eval(profiles[iprofile])]][i,j]<- (new.g2.int$val + s.int$val)/new.g1.int$val
DNA.aboveG1[[eval(profiles[iprofile])]][i,j]<- (new.g2.int$val + s.int$val + gt4.int$val)/new.g1.int$val
DNA.gt4[[eval(profiles[iprofile])]][i,j]<- gt4.int$val
DNA.lt2[[eval(profiles[iprofile])]][i,j]<- lt2.int$val
DNA.S[[eval(profiles[iprofile])]][i,j]<- s.int$val
DNA.fit.success[[eval(profiles[iprofile])]][i,j]<-  the.fit.success
DNA.A[[eval(profiles[iprofile])]][i,j]<- new.A
DNA.B[[eval(profiles[iprofile])]][i,j]<- new.B
 }
  
leg.txt<-c(paste("Red=",sum(red),sep=""),
             paste("Green=",sum(green),sep=""),
             paste("NotGreen=",sum(!green),sep=""),
             paste("Red and Green=",sum(green & red),sep=""),
             "--------------------------",
             "RATIOs ARE FOR  Green/NotGreen",
               paste("(S+G2)/G1 =",round((DNA.G1andG2$adenG/DNA.G1andG2$adenNG)[i,j],2),sep=""),
                paste("(S+G2+ >4N)/G1 =",round(( DNA.aboveG1$adenG/DNA.aboveG1$adenNG)[i,j],2),sep=""),
               paste("S =",round((  DNA.S$adenG/DNA.S$adenN)[i,j],2),sep=""),
               paste("> 4N =",round((DNA.gt4$adenG/DNA.gt4$adenN)[i,j],2),sep=""),
               paste("< 2N =",round((  DNA.lt2$adenG/DNA.lt2$adenN)[i,j],2),sep=""),
               "--------------------------",
               paste("Fit Quality: Green[NotGreen] =",round((DNA.fit.success$adenG)[i,j],2),"[",round((DNA.fit.success$adenNG)[i,j],2) ,"]" ," : G2/G1=",round(DNA.G2[i,j]/DNA.G1[i,j],2),sep=""))
legend(x=expected.g1.posn*2.5,the.max.range,legend=leg.txt,text.col=c("black","black","black","black","black","black","red","orange4","black","black","black"),pch="",bty="n",cex=1.25)

dev.off()

cells.P.field.vec <-(tapply(test[,"FieldIndex"],test[,"FieldIndex"],length))
R.P.field.vec <- (tapply(test[red,"FieldIndex"],test[red,"FieldIndex"],length))
G.P.field.vec<- (tapply(test[green,"FieldIndex"],test[green,"FieldIndex"],length))

cells.P.field[i,j] <-toString(cells.P.field.vec) # don't need to peocess
G.P.field.vec<- G.P.field.vec[names(cells.P.field.vec)]
G.P.field.vec[is.na(G.P.field.vec)]<-0
G.P.field[i,j]<- toString(G.P.field.vec)

R.P.field.vec<- R.P.field.vec[names(cells.P.field.vec)]
R.P.field.vec[is.na(R.P.field.vec)]<-0
R.P.field[i,j]<- toString(R.P.field.vec)
      
notGreen.c[i,j]<-sum(!green)
notRed.c[i,j] <-sum(!red)
redAndGreen[i,j] <-sum(red & green)

redAndGreenLow[i,j] <-sum(red & green & (test[,"AvgIntenCh2"] >= green.cut.low) )
redAndGreenMid[i,j] <-sum(red & green & (test[,"AvgIntenCh2"] >= green.cut.mid ))
redAndGreenHigh[i,j] <-sum(red & green & (test[,"AvgIntenCh2"] >= green.cut.high) )

#redLow[i,j] <-sum(red & (test[,"AvgIntenCh2"] >= green.cut.low) )
#redMid[i,j] <-sum(red & (test[,"AvgIntenCh2"] >= green.cut.low) & (test[,"AvgIntenCh2"] <= green.cut.high) )
#redHigh[i,j]<-sum(red & (test[,"AvgIntenCh2"] > green.cut.high) )

greenLow[i,j] <-sum( green & (test[,"AvgIntenCh2"] >= green.cut.low) )
greenMid[i,j] <-sum(green & (test[,"AvgIntenCh2"] >= green.cut.mid)  )
greenHigh[i,j]<-sum(green & (test[,"AvgIntenCh2"] >= green.cut.high)  )

notRedAndGreen[i,j] <-sum(!red & green)
redAndNotGreen[i,j]  <-sum(red & !green)

  } 
#id test big enough

 }} #loop over roes and columns

score<-redAndGreen*notGreen.c/(redAndNotGreen*green.c)   ###( red AND Green/Green)  / (red and NOT Green/ NOT Green)
mean(as.numeric(redAndGreen*notGreen.c/(redAndNotGreen*green.c)),na.rm=TRUE)

save(list=c("score","greenLow","greenMid","greenHigh","redAndGreenLow","redAndGreenMid","redAndGreenHigh","cells.P.field","R.P.field","G.P.field","big.c","red.c","green.c","notRed.c","redAndGreen","notGreen.c","notRedAndGreen","redAndNotGreen","all.c","lm.red.slope","lm.green.slope","lm.red.inter","lm.green.inter","lm.red.coverg","lm.green.coverg","DNA.G1andG2","DNA.aboveG1","DNA.gt4","DNA.lt2","DNA.S","DNA.fit.success","DNA.A","DNA.B","DNA.G1","DNA.G2"),file=paste(file.plate.name,"_SUMMARY",".RData",sep=""))

} ## loop over plates in one file

} ## loop over files

[1] "Ca_1"
Read 34536510 items
[1] "ARVEC02_14:07:31"
[1] "Ca_2"
Read 78717420 items
[1] "ARVEC02_04:45:25"
[1] "siCa1_3.1"
[1] "siCa1_3.2"
[1] "Ca_3"
Read 166233708 items
[1] "ARVEC02_15:58:29"
[1] "ARVEC02_17:46:35"
[1] "siCa1_2.1"
[1] "ARVEC02_21:24:00"
[1] "ARVEC02_23:13:01"
order run:


 DNA.G1andG2,DNA.aboveG1,DNA.gt4,DNA.lt2,DNA.S,DNA.fit.success,DNA.A,DNA.B,


save(list=c("score","greenLow","greenMid","greenHigh","redAndGreenLow","redAndGreenMid","redAndGreenHigh","cells.P.field","R.P.field","G.P.field","big.c","red.c","green.c","notRed.c","redAndGreen","notGreen.c","notRedAndGreen","redAndNotGreen","all.c","lm.red.slope","lm.green.slope","lm.red.inter","lm.green.inter","lm.red.coverg","lm.green.coverg"),file="hugo_004_SUMMARY3.RData")
save(list=c("score","greenLow","greenMid","greenHigh","redAndGreenLow","redAndGreenMid","redAndGreenHigh","cells.P.field","R.P.field","G.P.field","big.c","red.c","green.c","notRed.c","redAndGreen","notGreen.c","notRedAndGreen","redAndNotGreen","all.c","lm.red.slope","lm.green.slope","lm.red.inter","lm.green.inter","lm.red.coverg","lm.green.coverg"),file="hugo_002_SUMMARY3.RData")

save(list=c("score","greenLow","greenMid","greenHigh","redAndGreenLow","redAndGreenMid","redAndGreenHigh","cells.P.field","R.P.field","G.P.field","big.c","red.c","green.c","notRed.c","redAndGreen","notGreen.c","notRedAndGreen","redAndNotGreen","all.c","lm.red.slope","lm.green.slope","lm.red.inter","lm.green.inter","lm.red.coverg","lm.green.coverg"),file="hugo_001_SUMMARY3.RData")





save(list=c("score","big.c","red.c","green.c","notRed.c","redAndGreen","notGreen.c","notRedAndGreen","redAndNotGreen","all.c","lm.red.slope","lm.green.slope","lm.red.inter","lm.green.inter","lm.red.coverg","lm.green.coverg"),file="Plate2_V3.6_RESCAN_FINAL.RData")

save(list=c("score","big.c","red.c","green.c","notRed.c","redAndGreen","notGreen.c","notRedAndGreen","redAndNotGreen","all.c","lm.red.slope","lm.green.slope","lm.red.inter","lm.green.inter","lm.red.coverg","lm.green.coverg"),file="Plate4_V3.6_RESCAN_FINAL.RData")

save(list=c("score","big.c","red.c","green.c","notRed.c","redAndGreen","notGreen.c","notRedAndGreen","redAndNotGreen","all.c","lm.red.slope","lm.green.slope","lm.red.inter","lm.green.inter","lm.red.coverg","lm.green.coverg"),file="Plate1_V3.6_RESCAN_FINAL.RData")


 signif(redAndGreen*notGreen.c/(redAndNotGreen*green.c),2)

 signif(score ,2)
 signif(red.c  ,2)
 signif(green.c ,2)
 signif(notRed.c  ,2)
 signif(redAndGreen ,2)
 signif(notGreen.c ,2)
 signif(notRedAndGreen ,2)
 signif(redAndNotGreen  ,2)
 signif( lm.red.slope ,2)
 signif(lm.green.slope ,2)
 signif(lm.red.inter ,2)
 signif(lm.green.inter ,2)
  signif(big.c/all.c  ,1)
lm.red.coverg
lm.green.coverg

 rm(cells.num)
 save(list=c("score","red.c","green.c","notRed.c","redAndGreen","notGreen.c","notRedAndGreen","redAndNotGreen","all.c","lm.red.slope","lm.green.slope","lm.red.inter","lm.green.inter","lm.red.coverg","lm.green.coverg"),file="Plate4_V3_ANA.RData")
 save(list=c("score","red.c","green.c","notRed.c","redAndGreen","notGreen.c","notRedAndGreen","redAndNotGreen"),file="Plate1_V3_ANA.RData")

 score
 red.c
 green.c
 notRed.c
 redAndGreen
 notGreen.c
 notRedAndGreen
 redAndNotGreen
lm.red.slope
lm.green.slope
lm.red.inter
lm.green.inter
lm.red.coverg
lm.green.coverg  

 DE V3
>  red.c
     8   9
D    0 435
E 2206   0
>  notRed.c
      8    9
D     0 6672
E 14331    0
>  redAndGreen
     8   9
D    0 285
E 1057   0
>  notGreen.c
      8    9
D     0 1671
E 10558    0
>  notRedAndGreen
     8    9
D    0 5151
E 4922    0
>  redAndNotGreen
     8   9
D    0 150
E 1149   0

#######DE V4 no cut on Ch1 totalInt
>  red.c
     8   9
D    0 446
E 2302   0
>  notRed.c
      8    9
D     0 6740
E 14790    0
>  redAndGreen
     8   9
D    0 280
E 1067   0
>  notGreen.c
      8    9
D     0 1870
E 11106    0
>  notRedAndGreen
     8    9
D    0 5036
E 4919    0
>  redAndNotGreen
     8   9
D    0 166
E 1235   0




     [,1] [,2]
[1,]    0 1671
[2,]    0    0
> notRed.c
     [,1] [,2]
[1,]    0 6672
[2,]    0    0
> redAndGreen
     [,1] [,2]
[1,]    0  285
[2,]    0    0
> notGreen.c
     [,1] [,2]                   s
[1,]    0 1671
[2,]    0    0
> notRedAndGreen
     [,1] [,2]
[1,]    0 5151
[2,]    0    0
> redAndNotGreen
     [,1] [,2]
[1,]    0  150
[2,]    0    0


> red.c
     8 9
D    0 0
E 2206 0
> notRed.c
      8 9
D     0 0
E 14331 0
> redAndGreen
     8 9
D    0 0
E 1057 0
> notGreen.c
      8 9
D     0 0
E 10558 0
> notRedAndGreen
     8 9
D    0 0
E 4922 0
>  redAndNotGreen
     8 9
D    0 0
E 1149 0



>  red.c
    1    2   3   4   5    6    7    8    9   10  11  12
A 876  603 715 507 173  323  478  329  417 1240 891 779
B 746  797 835 890 977 1041 1122 1253  866  281 226 156
C 929  702 727 607 674  770  257  171  359  434 668 132
D 423  438 924 149 207  228  755  465  435 1004 905 666
E 312  473 616 271 228  279 2571 2206 1171  421 741 522
F 645 1204 658 665 616  767  694  904  451  351 417 102
G 111  860 750 437 386  482  569  482  789  773 386 239
H 347  801 741 996 809 1413  672  764  786  486 358  71
>  notRed.c
      1     2     3     4     5     6     7     8     9    10    11    12
A 15568  6712  9328  7892  3607  6593  5269  6166  5607  5992  7334 10345
B  8429 14515  6790  6719  4778  7504 11382 10926  8415  3924  4000  3812
C  8759  6937  6375  7859  5806  9766 20300 24179 17317  5574  8554  3599
D  4687 16671 12422 21855 20063 19810  6119  5576  6672 14833 10395  8591
E  4636 20161  6225 19638 18924 19810 18945 14331 20189 17105 21266 17231
F  9367 17547  6223  9808  8996  6592  7350  8203  6045  4313  6190  2365
G  2649  6790  5725 15395 16474 21815  6853  5450  7246  6859  4679  4071
H 10632  9008  7647 13072 10054 13692  7474  9361  8810 15297 14394  6265
>  redAndGreen
    1   2   3   4   5   6    7    8   9  10  11  12
A  12 331 247 187  85 145  184  146 193 994 673 518
B 191  99 381 105  89 107  283  300 312  98  96  76
C 215 226 301 301 421 347    5    6   6 169 223  74
D  65  22 266   0   0   0  471  328 285 383 414 355
E 122   0 345  10   6   9 1062 1057 321   6  24  16
F 197  30 355 114 149 127  101  118  65 118 140  72
G  56 429 425   4   1   3   13    7  12 291 135  69
H   4 367 346   4   4   4  305  383 393   0   0   0
>  notGreen.c
      1     2    3     4     5     6     7     8     9    10    11    12
A 16036  2687 5206  3602   763  2890  2354  2567  2517  1843  2908  5043
B  5818 12595 2969  6282  4568  7090  8521  7950  5077  1525  1241  1097
C  6793  3901 2862  3546  1624  5286 20360 24073 17506  2933  4607  1209
D  3375 15746 9199     1     1     1  1705  1078  1671  9608  5527  3786
E  1934     1 2279 19508 18773 19713 17476 10558 19399 17357 21474 17067
F  5703 18179 2246  8685  7024  5293  6425  7434  4952  1613  2454   421
G   577  2496 1836 15652 16402 22047  7113  5557  7705  3665  2381  1771
H 10400  4657 3487 14018 10816 14998  3434  4555  4082     1     1     1
>  notRedAndGreen
     1    2    3    4    5    6    7    8    9   10   11   12
A  396 4297 4590 4610 2932 3881 3209 3782 3314 4395 4644 5563
B 3166 2618 4275 1222 1098 1348 3700 3929 3892 2582 2889 2795
C 2680 3512 3939 4619 4435 4903  192  271  164 2906 4392 2448
D 1670 1341 3881    0    0    0 4698 4635 5151 5846 5359 5116
E 2892    0 4217  391  373  367 2978 4922 1640  163  509  670
F 4112  542 4280 1674 2439 1939 1518 1555 1479 2933 4013 1974
G 2127 4725 4214  176  457  247  296  368  318 3676 2549 2470
H  575 4785 4555   46   43  103 4407 5187 5121    0    0    0
>  redAndNotGreen
    1    2   3   4   5    6    7    8   9  10  11  12
A 864  272 468 320  88  178  294  183 224 246 218 261
B 555  698 454 785 888  934  839  953 554 183 130  80
C 714  476 426 306 253  423  252  165 353 265 445  58
D 358  416 658 149 207  228  284  137 150 621 491 311
E 190  473 271 261 222  270 1509 1149 850 415 717 506
F 448 1174 303 551 467  640  593  786 386 233 277  30
G  55  431 325 433 385  479  556  475 777 482 251 170
H 343  434 395 992 805 1409  367  381 393 486 358  71
>

signif(redAndGreen/redAndNotGreen,2 )
rat<-redAndGreen/redAndNotGreen
m<-mean(as.numeric(rat))
psd<-sd(as.numeric(rat))

signif(((redAndGreen/redAndNotGreen)-m)/psd,2 )


