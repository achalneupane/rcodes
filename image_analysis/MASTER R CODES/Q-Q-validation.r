

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
#######################################################




#### new annotation way
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

#####old annotatio way
## plates.all<-read.delim(well.output.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##Unnormalized
## colnames(plates.all)<-sum.vars
## normalized<-read.delim(paste(the.screen,"NORMALIZED","txt",sep="."),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
## colnames(normalized)<-sum.vars


#################################################################################
#normalized.file<-paste(the.screen,"NORMALIZED","txt",sep=".")
## targets<-c("Z-all Green","Z-Low Green","Z-Mid Green","Z-High Green","%R","Num-cells","MEAN_Area","MEAN_LWR","MEAN_P2A","MEAN_Total_I","SD_LWR","SD_P2A")
## targets<-c("%R","S+G2_G1","S","ak_leakage","resazurin","Cell_count","MEAN_Area","MEAN_LWR","MEAN_P2A","MEAN_Avg_Inten","MEAN_Total_Inten","SD_LWR","SD_P2A")
## targets<-c("percent_R", "ak_leakage","resazurin","Num-cells"  ) # ca  2 and 3
## targets<-c("S+G2_over_G1","percent_R", "ak_leakage","resazurin","Num-cells") # fisal

sort(tapply(plates.all[,"Symbol"],plates.all[,"Symbol"],length))
##two.color




########Jane validation
pos.control.list<-c("CCND1","C-MYC")
pos.control.list2<-c("HSPCN")
norm.control.list<-c("plv411") # controls that could be used for  normalization
confounders.list<-c("UTMOCK","MOCK")


#############################################################################################
#############################################################################################
#############################################################################################
################################# restart from written files  Q_Q PLOTS
#plates.all<-plates.all.ori
plates.all.ori<-plates.all

## for(i in 1:
## a.screen.cuts.ori<-a.screen.cuts
a.screen.cuts<-list()

low.cut<-rep(0,times=length(targets))
high.cut<-rep(0,times=length(targets))
names(low.cut)<-targets
names(high.cut)<-targets


## the.screens<-c("si FB1","si FB2","si FB3")
## the.screens<-c("MUS.a","MUS.b","MUS.c")
the.screens<-c("a","b","c","4a","4b","4c")



for(i in 1:length(the.screens)){
cuts<-list(low.cut=low.cut,high.cut=high.cut)
 a.screen.cuts[the.screens[i]]<- list(cuts)
}
a.screen.cuts




########## add extra parameters to an existing screen:
## the.screens<-c("MUS.a","MUS.b","MUS.c")
## the.screens<-c("a","b","c")
## for(i in 1:length(the.screens)){
##   for(j in 7:length(targets)){
## a.screen.cuts[[the.screens[i]]]$low.cut[targets[j]]<-0
## a.screen.cuts[[the.screens[i]]]$high.cut[targets[j]]<-0
## }
## }

#######################
a.screen.stats<-a.screen.cuts
a.screen.counts<-a.screen.cuts

the.screens
to.exclude<-c("4a","4b","4c") # to.exclude<-c("a","b","c")

run.screens<-the.screens[!(the.screens %in% to.exclude)]
run.screens
targets

recount<-FALSE ### use TRUE if have existing thresholds:  reget the statistics but do no mofify thresholds
cut.sd.posn<-1.64 ## 1.64 for  p=0.5 for upper limit  ## cut.sd.posn<-2 for p=0.03 in upper tail
reps.for.hit<-2 ## reps.for.hit<-1

for(k in 1:length(run.screens)){
  
for(i in 1:length(targets)){

    
plates.all<-plates.all.ori

the.screen<-run.screens[k]
the.score<-targets[i]

if(the.screen=="a" & the.score=="RG_over_G_over_EtoH"){next} ### cause it's all 1 by defination
if(the.screen=="4a" & the.score=="RG_over_G_over_EtoH"){next}

print(the.screen)
wanted.screen<-grepl(the.screen,plates.all[,"plate"])
and.exclude<-!(plates.all[,"plate"] %in% to.exclude)

sum(wanted.screen & and.exclude)
plates.all<-plates.all[wanted.screen & and.exclude,]
print(unique(plates.all[,"plate"]))


print(the.score)
plates.all[1:5,the.score]


dim(plates.all)

confounders<-plates.all[,"Symbol"] %in% confounders.list

are.finite<-is.finite(plates.all[,the.score])

filter<-are.finite & !confounders

if(the.score=="percent_R"){filter<- filter & plates.all[,"R"]>=60} # need mre than 60 red cells to make a judgement
if(the.score=="RG_over_G"){filter<- filter & plates.all[,"G"]>=50}
if(the.score=="RG_over_G_Z"){filter<- filter & plates.all[,"G"]>=50}
if(the.score=="Gr-G1"){filter<- filter & plates.all[,"G"]>=50}
if(the.score=="RG_over_G_over_EtoH"){filter<- filter & plates.all[,"G"]>=50 & plates.all[,"RG_over_G_over_EtoH"]<=10 }
if(the.score=="Z-all Green"){filter<- filter & plates.all[,"RG"]>=35}
if(the.score=="Num-cells"){filter<- filter & plates.all[,"RG"]>=70}

dim(plates.all[filter,]) # plates.all[filter & (plates.all[,"Symbol"] %in% pos.control.list ),1:15]

## plates.all[(plates.all[,"Symbol"] %in% pos.control.list ),1:15]


labels<-as.character(paste(plates.all[filter,"plate"],"-",plates.all[filter,"row"],plates.all[filter,"col"],"-",plates.all[filter,"Symbol"],sep=""))
xlabels<-as.character(paste(plates.all[filter,"plate"],"-",plates.all[filter,"row"],plates.all[filter,"col"],sep=""))
symbols<-plates.all[filter,"Symbol"]
plates<-plates.all[filter,"plate"]
cols<-plates.all[filter,"col"]
data.in<- plates.all[filter,the.score]

pos.control<-symbols %in% pos.control.list
norm.control<-symbols %in% norm.control.list

if(validation){
the.mean<-mean(data.in[norm.control])
the.sd<-sd(data.in[norm.control])
norm.test<-shapiro.test(data.in[norm.control])
print(norm.test)
}else{
the.mean<-mean(data.in[!pos.control])
the.sd<-sd(data.in[!pos.control])
norm.test<-shapiro.test(data.in[pos.control])
print(norm.test)
}

print(the.mean)
print(the.sd)

range(data.in)

if(recount){
the.low.cut<-a.screen.cuts[[the.screen]]$low.cut[the.score]
the.high.cut<-a.screen.cuts[[the.screen]]$high.cut[the.score]
}else{
the.low.cut<-the.mean-cut.sd.posn*the.sd
sum(data.in<the.low.cut)/length(data.in)

the.high.cut<-the.mean+cut.sd.posn*the.sd
sum(data.in>the.high.cut)/length(data.in)

a.screen.cuts[[the.screen]]$low.cut[the.score]<-the.low.cut
a.screen.cuts[[the.screen]]$high.cut[the.score]<-the.high.cut
}

print(sort(tapply(symbols[data.in>=the.high.cut],symbols[data.in>=the.high.cut],length),decreasing=TRUE))

genes.high.cut<-sort(tapply(symbols[data.in>=the.high.cut],symbols[data.in>=the.high.cut],length),decreasing=TRUE)
high.passing<-sum(genes.high.cut>=reps.for.hit)
no.control.high<-genes.high.cut[names(genes.high.cut) %in% norm.control.list]
print(high.passing)


print(sort(tapply(symbols[data.in<=the.low.cut],symbols[data.in<=the.low.cut],length),decreasing=TRUE))

genes.low.cut<-sort(tapply(symbols[data.in<=the.low.cut],symbols[data.in<=the.low.cut],length),decreasing=TRUE)
low.passing<-sum(genes.low.cut>=reps.for.hit)
no.control.low<-genes.low.cut[names(genes.low.cut) %in% norm.control.list]
print(low.passing)

a.screen.stats[[the.screen]]$low.cut[the.score]<-toString(paste("Hits:",low.passing,"| Num Controls:",no.control.low,"| Shapiro:",signif(norm.test$p.value,digits=2),sep=" "))
a.screen.stats[[the.screen]]$high.cut[the.score]<-toString(paste("Hits:",high.passing,"| Num Controls:",no.control.high,"| Shapiro:",signif(norm.test$p.value,digits=2),sep=" "))

a.screen.counts[[the.screen]]$low.cut[the.score]<-toString(paste(names(genes.low.cut),genes.low.cut,sep=":"))
a.screen.counts[[the.screen]]$high.cut[the.score]<-toString(paste(names(genes.high.cut),genes.high.cut,sep=":"))

print("------------------------------------")
}
}

save(list=c("a.screen.cuts"),file="a.screen.cuts.RData")
save(list=c("a.screen.counts"),file="a.screen.counts.RData")
save(list=c("a.screen.stats"),file="a.screen.stats.RData")
##################################### end setting cutoffs
a.screen.cuts.ori<-a.screen.cuts
a.screen.counts.ori<-a.screen.counts
a.screen.stats.ori<-a.screen.stats

a.screen.stats # a.screen.stats<-a.screen.stats.ori
a.screen.counts

a.screen.stats.ori
a.screen.counts.ori

a.screen.stats.ori<-a.screen.stats

######################################one off- q-q plots are below

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################


targets
names(a.screen.cuts)

the.score<-"Gr-G1"
the.screen<-"b"

the.score<-"RG_over_G_Z"
the.score<-"Z-all Green"
the.screen<-"b"

run.screens<-the.screens

for(k in 1:length(run.screens)){
  
for(i in 1:length(targets)){

  
## plates.all<-plates.all.ori

the.screen<-run.screens[k]
the.score<-targets[i]
if(the.screen=="a" & the.score=="RG_over_G_over_EtoH"){next} ### cause it's all 1 by defination
if(the.screen=="4a" & the.score=="RG_over_G_over_EtoH"){next}

wanted.screen<-grepl(the.screen,plates.all.ori[,"plate"])
and.exclude<-!(plates.all.ori[,"plate"] %in% to.exclude)
plates.all<-plates.all.ori[wanted.screen & and.exclude,]

a.screen.cuts[[the.screen]]$high.cut[the.score] ## HIGH
a.screen.stats[[the.screen]]$high.cut[the.score] ## HIGH

a.screen.cuts[[the.screen]]$low.cut[the.score]  ## LOW
a.screen.stats[[the.screen]]$low.cut[the.score] ## LOW

########################### Z score to P-value ##########################
##  ## erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
##  erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE) ### 1-erf:  2 tailed sigma=1 68% p=0.32; sigma=2 95% p=0.05;  sigma=3 99.7% p=0.003;  = 2*erfc(x)
## z<-0
## erfc(z/sqrt(2))
##  ### this is 2 tailed (eg. higher and lower than 2 sd if x=2 )
## 1-pnorm(x, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE) ### this is 1 tailed (eg. higher  2 sd if x=2 )
## need to abs(z)
## x<-1:6
## signif(2*(pnorm(x, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)),digits=3) ## is the 2 tailed one to USE!
## x<-1.64  # p=0.05 for upper tail  ## x=2 p=0.03 in upper tail
## signif(
##        (pnorm(x, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))

##        ,digits=3) 
#########################################################################


confounders.list
confounders<-plates.all[,"Symbol"] %in% confounders.list
are.finite<-is.finite(plates.all[,the.score])

filter<-are.finite & !confounders
if(the.score=="percent_R"){filter<- filter & plates.all[,"R"]>=60} # need mre than 60 red cells to make a judgement
if(the.score=="RG_over_G"){filter<- filter & plates.all[,"G"]>=50}
if(the.score=="RG_over_G_Z"){filter<- filter & plates.all[,"G"]>=50}
if(the.score=="Gr-G1"){filter<- filter & plates.all[,"G"]>=50}
if(the.score=="RG_over_G_over_EtoH"){filter<- filter & plates.all[,"G"]>=50 & plates.all[,"RG_over_G_over_EtoH"]<=10 }
if(the.score=="Z-all Green"){filter<- filter & plates.all[,"RG"]>=35}
if(the.score=="Num-cells"){filter<- filter & plates.all[,"RG"]>=70}

dim(plates.all[filter,]) # plates.all[filter & (plates.all[,"Symbol"] %in% pos.control.list ),1:15]

## plates.all[(plates.all[,"Symbol"] %in% pos.control.list ),1:15]


labels<-as.character(paste(plates.all[filter,"plate"],"-",plates.all[filter,"row"],plates.all[filter,"col"],"-",plates.all[filter,"Symbol"],sep=""))
xlabels<-as.character(paste(plates.all[filter,"plate"],"-",plates.all[filter,"row"],plates.all[filter,"col"],sep=""))
symbols<-plates.all[filter,"Symbol"]
plates<-plates.all[filter,"plate"]
cols<-plates.all[filter,"col"]
data.in<- plates.all[filter,the.score]

sort(tapply(symbols[data.in>=a.screen.cuts[[the.screen]]$high.cut[the.score] ],symbols[data.in>=a.screen.cuts[[the.screen]]$high.cut[the.score] ],length),decreasing=TRUE)
sum(sort(tapply(symbols[data.in>=a.screen.cuts[[the.screen]]$high.cut[the.score] ],symbols[data.in>=a.screen.cuts[[the.screen]]$high.cut[the.score] ],length),decreasing=TRUE)>=2)

pos.control<-symbols %in% pos.control.list
norm.control<-symbols %in% norm.control.list

if(validation){
the.mean<-mean(data.in[norm.control])
the.sd<-sd(data.in[norm.control])
}else{
the.mean<-mean(data.in[!pos.control])
the.sd<-sd(data.in[!pos.control])
}
the.mean
the.sd
the.mean+cut.sd.posn*the.sd
the.mean-cut.sd.posn*the.sd
range(data.in)


qq<- qq.data(data.in,distribution="norm",the.mean=the.mean,the.sd=the.sd,plot.it=FALSE)

my.qq.plot(data.in,distribution="norm",col="blue",xlab="Expected Score",ylab="Observed score",xlim=range(qq$x), ylim=range(data.in),main=paste("Screen:",the.screen,"with 95% confidence intervals for",":",the.score,sep=" "),the.mean=the.mean,the.sd=the.sd,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex=1.5)

## my.qq.plot(data.in,distribution="norm",col="blue",xlab="Expected Score",ylab="Observed score", xlim=c(1.1,2.2), ylim=c(1,1.4),main=paste(the.screen,"with 95% confidence intervals",":",the.score,sep=" "),the.mean=the.mean,the.sd=the.sd,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex=1.5)


posns<-c(1:length(symbols))[symbols[qq$ord] %in% norm.control.list]
points(qq$x[posns],qq$y[posns],col="cyan",pch=19,cex=0.85)

posns<-c(1:length(symbols))[symbols[qq$ord] %in% pos.control.list ]
points(qq$x[posns],qq$y[posns],col="red",pch=19,cex=1.0)
 
abline(h=a.screen.cuts[[the.screen]]$high.cut[the.score],col="red",lwd=2)
abline(h=a.screen.cuts[[the.screen]]$low.cut[the.score],col="red",lwd=2)
text(range(qq$x)[1],a.screen.cuts[[the.screen]]$high.cut[the.score],labels=paste("hits>",signif(a.screen.cuts[[the.screen]]$high.cut[the.score],digits=3),sep=""),pos=3,col="red",cex=1.5)
text(range(qq$x)[1],a.screen.cuts[[the.screen]]$low.cut[the.score],labels=paste("hits<",signif(a.screen.cuts[[the.screen]]$low.cut[the.score],digits=3),sep=""),pos=1,col="red",cex=1.5)

legend(x=range(qq$x)[1],y=range(data.in)[2],legend=c(toString(norm.control.list),toString(pos.control.list),toString(pos.control.list2)),col=c("cyan","red","pink"),pch=c(20,20,20,21),cex=1.5)
selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="red",cex=1,atpen='TRUE') ##plate row col symbol
############ if cutoffs are known
savePlot(filename=paste("Q-Q",the.screen,the.score,"jpeg",sep="."),type="jpeg")


}} # end loop over all plots


#####annotate curve
selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="red",cex=1,atpen='TRUE') ##plate row col symbol
selected.data<-identify(qq$x,qq$y,labels=labels[qq$ord],col="red",cex=1,atpen='TRUE') ## sybmol
selected.data<-identify(qq$x,qq$y,labels=as.character(round(data.in[qq$ord],2)),col="forestgreen",cex=1.25,atpen='TRUE') # observed score
#####


savePlot(filename=paste("Q-Q",the.screen,the.score,"jpeg",sep="."),type="jpeg")

###########redefine


the.low.cut<-the.mean-2*the.sd
the.high.cut<-the.mean+4*the.sd


the.low.cut<-4.96

the.high.cut<-4.2

the.high.cut/the.sd

a.screen.cuts[[the.screen]]$low.cut[the.score]<-the.low.cut

a.screen.cuts[[the.screen]]$high.cut[the.score]<-the.high.cut




##########################
my.qq.plot<-function (x, distribution = "norm", ylab = deparse(substitute(x)),
    xlab = paste(distribution, "quantiles"), main = NULL, las = par("las"), 
    envelope = 0.95, labels = FALSE, col = palette()[2], lwd = 2, the.mean=0,the.sd=1,cex.lab=2,
    pch = 1, cex = 1, line = c("quartiles", "robust", "none"),xlim=c(0,100),ylim=c(0,20),font.lab=2,font.axis=2,font.main=2,cex.axis=1,cex.main=1,
    ...)
{
    result <- NULL
    line <- match.arg(line)
    good <- !is.na(x)
    ord <- order(x[good])
    ord.x <- x[good][ord]
    q.function <- eval(parse(text = paste("q", distribution, 
        sep = "")))
    d.function <- eval(parse(text = paste("d", distribution, 
        sep = "")))
    n <- length(ord.x)
    P <- ppoints(n)
    z <- q.function(P, mean=the.mean, sd=the.sd, ...)
    plot(z, ord.x, xlab = xlab, ylab = ylab, main = main, las = las, 
        col = col, pch = pch,cex = cex,xlim=xlim,ylim=ylim,cex.lab=cex.lab,font.lab=font.lab,font.axis=font.axis,font.main=font.main,cex.main=cex.main,cex.axis=cex.axis)
    if (line == "quartiles") {
        Q.x <- quantile(ord.x, c(0.25, 0.75))
        Q.z <- q.function(c(0.25, 0.75), mean=the.mean, sd=the.sd, ...)
        b <- (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
        a <- Q.x[1] - b * Q.z[1]
        abline(a, b, col = "red", lwd = lwd)
    }
    if (line == "robust") {
        if (!require("MASS")) 
            stop("MASS package not available")
        coef <- coefficients(rlm(ord.x ~ z))
        a <- coef[1]
        b <- coef[2]
        abline(a, b)
    }         ###################  Envelope function
    if (line != "none" & envelope != FALSE) {
        zz <- qnorm(1 - (1 - envelope)/2)
        SE <- (b/d.function(z, mean=the.mean, sd=the.sd, ...)) * sqrt(P * (1 - P)/n)
        fit.value <- a + b * z
        upper <- fit.value + zz * SE
        lower <- fit.value - zz * SE
        lines(z, upper, lty = 2, lwd = lwd/2, col = "red")
        lines(z, lower, lty = 2, lwd = lwd/2, col = "red")
    }       #####################
    if (labels[1] == TRUE & length(labels) == 1)
        labels <- seq(along = z)
    if (labels[1] != FALSE) {
        selected <- identify(z, ord.x, labels[good][ord])
        result <- seq(along = x)[good][ord][selected]
    }
    if (is.null(result)) 
        invisible(result)
    else sort(result)
}






      qq.data<- function (x, plot.it = TRUE, distribution = "norm", df=1, the.mean=0, the.sd=1,  xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)) , ...)
{
    good <- !is.na(x)
    ord <- order(x[good])
    ord.x <- x[good][ord]
    q.function <- eval(parse(text = paste("q", distribution, 
        sep = "")))
    n <- length(ord.x)
    P <- ppoints(n)
    z <- q.function(P, mean=the.mean, sd=the.sd, ...)

    if (plot.it)
        plot(z, ord.x, xlab = xlab, ylab = ylab, ...)
    invisible(list(x = z, y = ord.x, ord=ord))
}


#####################################################################################################################












