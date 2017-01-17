



cluster.definition.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/Final_FANC_clusters.csv"
clusters<-read.delim(cluster.definition.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
cluster.names<-colnames(clusters)
clinical.genes<-clusters[,"Clinical"]
clinical.genes<-clinical.genes[clinical.genes!=""]



working.dir<-c("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis")
project.extension<-"_analysis-maf-filtered.txt.txt"
project.name<-"2015-08-14_AML_mixedAligners."
fam<-c(".ALL.ALL_GENOTYPES")



#snap.file<- "NEWwBIASwQG.WITH.RECOVERY" # name of the run
snap.file<- "NEWwBIASwQG" # name of the run



do.MAFS<-c(0.001,0.01) # do.MAFS<-c(0.001,0.005,0.01)
#do.MAFS<-c(0.001)

ido.mafs<-1
for (ido.mafs in 1:length(do.MAFS)){
p<-do.MAFS[ido.mafs]
types<-c("Single Point","non-coding","bad.effect","sliding.window","coding")# ,"sliding.window")
#types<-c("coding") 

itypes<-1

       
for(itypes in 1:length(types)){
top.gene.list<-{}
the.test<-paste(types[itypes],p,sep=".")
#############
prefixes<-c(paste("Burden.ALL.",snap.file,the.test,".2015-08-14_AML_mixedAligners.",sep=""),paste("SKATO.ALL.",snap.file,the.test,".2015-08-14_AML_mixedAligners.",sep=""),paste("EVERYTHING.",p,".GENOTYPE.conponents..",snap.file,the.test,".2015-08-14_AML_mixedAligners.",sep=""))
the.extension<-".ALL.ALL_GENOTYPES_analysis-maf-filtered.txt.txt"

print(prefixes)
ip<-2

for(ip in 1:length(prefixes)){

    prefix<-prefixes[ip]

    if(grepl("^SKATO",prefix) & types[itypes]=="Single Point"){print("SKIP");next} # no skatO for single point
#    if(grepl("^EVERYTHING.",prefix) & !types[itypes]=="coding"){print("SKIP");next} # only one everything do on coding so fits the name


print(prefix)
print(the.extension)

setwd(working.dir)
files<-dir(getwd())

project.files<-files[grepl(the.extension,files)]

project.files<-project.files[grepl(paste0("^",prefix),project.files)]


print(sort(paste("Doing: ",project.files,sep=""))) # project.files<-project.files[1:22]

indels<-{}
the.col<-{}
project.files


  if(length(project.files)!=24){
    print("########################################### WARNING #################################")
    print("less that 24 chromosomes detected")
    print("########################################### WARNING #################################") 
  }

  project.files
length(project.files)

##########################

indels<-{}
ichr<-1

for(ichr in 1:length(project.files)){
    
    
   # print(project.files[ichr])
    
    ################## fast read ###########
    column.labels<-read.delim(project.files[ichr],header=F,nrows=1,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
    num.vars<-dim(column.labels)[2]
    a.indel<-scan(project.files[ichr],what=character(num.vars),skip=1,sep="\t",fill=TRUE,na.strings="",quote="\"")
    num.lines<-length(a.indel)/(num.vars)
    dim(a.indel)<-c(num.vars,num.lines)
    a.indel<-t(a.indel)
    colnames(a.indel)<-column.labels
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
    
    
    rm(a.indel)

     
  } # ichr

indels[1:20,]
dim(indels)

are.clusters<-indels[,"gene"] %in% cluster.names
indels<-indels[!are.clusters,]

order.by<-order(indels[,"p"])
indels<-indels[order.by,]
indels<-indels[!(indels[,"gene"] %in% cluster.names),]
enum<-1:dim(indels)[1]
indels<-cbind(enum,indels)

indels[1:20,]

if(length(top.gene.list)==0){top.gene.list<-indels[1:200,"gene"]}else{top.gene.list<-unique(c(top.gene.list,indels[1:200,"gene"]))}

dim(indels)
paste(prefix,"data_summary.txt",sep="")
write.table(indels,file=paste(prefix,"data_summary.txt",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


  if(grepl("^EVERYTHING.",prefix)){
 write.table(indels[indels[,"Gene.Names"] %in% top.gene.list,],file=paste(prefix,"TOP_200_summary.txt",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
}


  } # ip prefix

} # types

} # maf
  ################## write per family data
  
  




######################################################################################################


present<- indels[,"Gene.Names"] %in% meta.results.burden[1:100,"gene"]
sum(present)

  write.table(indels[present,],file=paste(prefix,"data_100.txt",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



cluster.definition.file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/Final_FANC_clusters.csv"
clusters<-read.delim(cluster.definition.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
cluster.names<-colnames(clusters)
clinical.genes<-clusters[,"Clinical"]
clinical.genes<-clinical.genes[clinical.genes!=""]

## file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/Burden.ALL.coding.0.001.2015-08-14_AML_mixedAligners.data_summary.txt"
## file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/Burden.ALL.No.missing.0.05coding.0.001.2015-08-14_AML_mixedAligners.data_summary.txt"
## file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/Burden.ALL.ORIGINALwBIAScoding.0.001.2015-08-14_AML_mixedAligners.data_summary.txt"
## file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/Burden.ALL.NEWwBIAScoding.0.001.2015-08-14_AML_mixedAligners.data_summary.txt"
## file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/Burden.ALL.NEWwBIASwQGcoding.0.001.2015-08-14_AML_mixedAligners.data_summary.txt"



file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/NO_RECOVERY_PUBLICATION/Burden.ALL.NEWwBIASwQGcoding.0.001.2015-08-14_AML_mixedAligners.data_summary.txt"
file<-"/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/NO_RECOVERY_PUBLICATION/Burden.ALL.NEWwBIASwQGnon-coding.0.001.2015-08-14_AML_mixedAligners.data_summary.txt"


setwd("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/NO_RECOVERY_PUBLICATION")
## prefixes<-c(paste("Burden.ALL.",snap.file,the.test,".2015-08-14_AML_mixedAligners.",sep=""),paste("SKATO.ALL.",snap.file,the.test,".2015-08-14_AML_mixedAligners.",sep=""),paste("EVERYTHING.",p,".GENOTYPE.conponents..",snap.file,the.test,".2015-08-14_AML_mixedAligners.",sep=""))
## the.extension<-".data_summary.txt"
 file.prefix<-gsub(".txt","",basename(file))

meta.results.burden<-read.delim(file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

meta.results.burden[1:15,]

top.500<-meta.results.burden$gene[1:500]
## posns<-match(meta.results.burden[,"gene"],rownames(poss.model))
## meta.results.burden<-cbind(poss.model[posns,],meta.results.burden)
## meta.results.burden[1:5,1:10]
## write.table(meta.results.burden,file="containamination loci.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)




clusters.wanted<-cluster.names

dups<-duplicated(meta.results.burden[,"gene"])
sum(dups)
meta.results.burden[dups,]

meta.results.burden<-meta.results.burden[!dups,]

meta.results.burden[1:5,]
       dim(meta.results.burden)
## bad.genes<-c(bad.genes,clusters.wanted)
subset<-meta.results.burden[ !(meta.results.burden[,"gene"] %in%  clusters.wanted) ,]
       subset[1:10,]

sum(!is.na(as.numeric(subset[ ,"p"])))
subset<-subset[!is.na(as.numeric(subset[ ,"p"])),]

## subset[1:10,]
## no.weight.sig.genes<-subset[as.numeric(subset[ ,"p"])<1.2e-6,"gene"]
## tight.weight.sig.genes<-subset[as.numeric(subset[ ,"p"])<1.2e-6,"gene"]

symbols<-subset[,"gene"]

## subset[1:10,]
## posns<-match(interesting,subset[,"gene"])
## names(posns)<-interesting
## posns
## write.table(subset[1:350,],file="string.txt",,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## dups<-duplicated(subset[,"gene"])
## sum(dups)

## write.table(subset[1:350,],file="string.txt",,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
z1<-qchisq(subset[ ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)
z<-qchisq(subset[ ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) ## if have no chisq valuse
#z0<-qchisq(meta.results.skatO[ !(meta.results.skatO[,"gene"] %in% genes.and.clusters ) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE) 
z[1:5]
subset[1:10,]

length(z)
dim(subset)
exclude.from.lambda<-subset[,"gene"] %in% clinical.genes
sum(exclude.from.lambda)

p.value<-as.numeric(subset[ ,"p"])
par(mfrow=c(1,1),font=2,font.lab=2,font.axis=2,mgp=c(3.5,1,0),mar=c(5,5,4,2)+0.1)


##  z<-rchisq(length(p.val), df=1, ncp = 0) ## null test


median(z[!exclude.from.lambda],na.rm=TRUE)/0.456  #1.071491
#median(z0,na.rm=TRUE)/0.456  #1.071491


sum(is.na(p.val))
################## p-values
z0=qnorm(p.value[!is.na(p.value)]/2)
lambda = round(median(z0[!exclude.from.lambda]^2)/0.454,3)
lambda


plot(x=c(1:10))
## source("http://bioconductor.org/biocLite.R") 
## biocLite("GWASTools")
## setRep
## qq.Plot(pvals)


help[1:5,]
## Reads data
## S <- read.table(input,header=F)
## if (stat_type == "Z")
##    z=S[,1]
## if (stat_type == "CHISQ")
##    z=sqrt(S[,1])
## if (stat_type == "PVAL")
##    z0=qnorm(meta.results.skatO[,"p"]/2)
## ## calculates lambda
lambda = round(median(z0^2)/.454,3)
## lambda
signif(range(z),digits=3)

z[!is.na(z)]
the.plot<-my.qq.plot(z,dist="chisq",df=1,ncp=0,col="blue",ylab="Observed chi-squared value",xlab="Expected chi-squared value",main="",cex=1,xlim=c(0,18),ylim=c(0,130),cex.lab=2.0,cex.axis=2.0,font.lab=2,font.axis=2,lwd=2,line="robust",plot.it=TRUE) # function defined below

the.plot<-my.qq.plot(z,dist="chisq",df=1,ncp=0,col="blue",ylab="Observed chi-squared value",xlab="Expected chi-squared value",main="",cex=1,xlim=c(7,18),ylim=c(0,130),cex.lab=2.0,cex.axis=2.0,font.lab=2,font.axis=2,lwd=2,line="robust",plot.it=TRUE)


## z.all<-qchisq(meta.results.burden[ !(meta.results.burden[,"gene"] %in% clusters.wanted) ,"p"],df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)
## range(z.all)
qq<-  qq.data(z,plot.it=FALSE)       ## qq plot used same method as in car library
#points(qq$x,qq$y,col="magenta",pch=21)

symbols<-subset[,"gene"]



#selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="forestgreen",cex=1.3,offset=1,atpen='TRUE')



interesting<-clinical.genes # c("ASXL1", "TET2", "FLT3", "IDH1", "RUNX1", "NRAS", "IDH2", "DNMT3A", "NPM1")
show<-symbols[qq$ord] %in% interesting
points(qq$x[show],qq$y[show],col="forestgreen",pch=20,cex=2,lwd=2)

selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="forestgreen",cex=1.0,offset=1,atpen='TRUE')
selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="blue",cex=1.0,offset=1,atpen='TRUE')

selected.data<-identify(qq$x,qq$y,labels=symbols[qq$ord],col="blue",cex=1.3,offset=3,atpen='TRUE')

toString(symbols[qq$ord][selected.data])

bonferroni<-0.05/dim(subset)[1]

abline(h=qchisq(bonferroni,df=1,ncp=0,lower.tail=FALSE,log.p=FALSE),col="blue",lw=2)


leg.txt<-c(expression(paste("Genes")),"95% confidence intervals")
legend(8,30,leg.txt,col=c("blue","red"),lty=c(-1,2),pch=c(1,-1),cex=2.0)


## leg.txt<-c(expression(paste("Genes")),"95% confidence intervals","Genes with multiple mutations","Genes with focal mutations")
## legend(2,100,leg.txt,col=c("blue","red","magenta","forestgreen"),lty=c(-1,2,-1,-1),pch=c(1,-1,19,19),cex=2.0)


#fig.prefix<-"Non-Coding_HIts"

fig.prefix<-paste(file.prefix,"Q-Q_plot",sep=".")

savePlot(filename=paste(fig.prefix,".png",sep=''),type="png")
savePlot(filename=paste(fig.prefix,".jpeg",sep=''),type="jpeg")
savePlot(filename=paste(fig.prefix,".tiff",sep=''),type="tiff")
savePlot(filename=paste(fig.prefix,".bmp",sep=''),type="bmp")
dev.print(svg,paste(fig.prefix,".svg",sep=''))









#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################

my.qq.plot<-function (x, distribution = "chisq",df=1,ncp=0, ylab = deparse(substitute(x)),
    xlab = paste(distribution, "quantiles"), main = NULL, las = par("las"), 
    envelope = 0.95, labels = FALSE, col = palette()[2], lwd = 2, 
    pch = 1, cex = 1, line = c("quartiles", "robust", "none"),xlim=c(0,100),ylim=c(0,20),font.lab=2,font.axis=2,font.main=2,cex.lab=2.5,cex.axis=1.0,plot.it=TRUE, ...){
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
    z <- q.function(P,df=df,ncp=ncp, ...)
    if(plot.it){
    plot(z, ord.x, xlab = xlab, ylab = ylab, main = main, las = las, 
        col = col, pch = pch,cex = cex,xlim=xlim,ylim=ylim,font.lab=font.lab,font.axis=font.axis,font.main=2,cex.lab=cex.lab,cex.axis=cex.axis)}
    if (line == "quartiles") {
        Q.x <- quantile(ord.x, c(0.25, 0.75))
        Q.z <- q.function(c(0.25, 0.75),df=df,ncp=ncp, ...)
        b <- (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
        a <- Q.x[1] - b * Q.z[1]
        if(plot.it){
        abline(a, b, col = "red", lwd = lwd)}
    }
    if (line == "robust") {
        if (!require("MASS")) 
            stop("MASS package not available")
        coef <- coefficients(rlm(ord.x ~ z))
        a <- coef[1]
        b <- coef[2]
          if(plot.it){
        abline(a, b,col="red")}
    }         ###################  Envelope function
    if (line != "none" & envelope != FALSE) {
        zz <- qnorm(1 - (1 - envelope)/2)
        SE <- (b/d.function(z,df=df,ncp=ncp, ...)) * sqrt(P * (1 - P)/n)
        fit.value <- a + b * z
        upper <- fit.value + zz * SE
        lower <- fit.value - zz * SE
          if(plot.it){
        lines(z, upper, lty = 2, lwd = lwd, col = "red")
        lines(z, lower, lty = 2, lwd = lwd, col = "red")}
    }       #####################
    if (labels[1] == TRUE & length(labels) == 1)
        labels <- seq(along = z)
    if (labels[1] != FALSE) {
        selected <- identify(z, ord.x, labels[good][ord])
        result <- seq(along = x)[good][ord][selected]
    }
    if (is.null(result)) 
        invisible(list(result=result,a=a,b=b,x=z,y = ord.x,ord=ord,upper=upper,lower=lower))
    else {sort(result)
           invisible(list(result=result,a=a,b=b,x=z,y = ord.x,ord=ord,upper=upper,lower=lower))}
}



      qq.data<- function (x, plot.it = TRUE, distribution = "chisq", df=1,ncp=0, xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)) , ...)
{
    good <- !is.na(x)
    ord <- order(x[good])
    ord.x <- x[good][ord]
    q.function <- eval(parse(text = paste("q", distribution, 
        sep = "")))
    n <- length(ord.x)
    P <- ppoints(n)
    z <- q.function(P,df=df,ncp=ncp, ...)

    if (plot.it)
        plot(z, ord.x, xlab = xlab, ylab = ylab, ...)
    invisible(list(x = z, y = ord.x, ord=ord))
} ##ord is the order if use identify


######################################### END SECTION

## qq<- qq.data(data.in,distribution="norm",the.mean=the.mean,the.sd=the.sd,plot.it=FALSE)

## my.qq.plot(data.in,distribution="norm",col="blue",xlab="Expected Score",ylab="Observed score",xlim=range(qq$x), ylim=range(data.in),main=paste("Screen:",the.screen,"with 95% confidence intervals for",":",the.score,sep=" "),the.mean=the.mean,the.sd=the.sd,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex=1.5)
########################## USE FUNTIONS BELOW IF YOU REFERNCE FUNCYION IS A NORMAL DISTRIBUTION NOT A CHISQ
my.qq.plot.mean<-function (x, distribution = "norm", ylab = deparse(substitute(x)),
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






      qq.data.mean<- function (x, plot.it = TRUE, distribution = "norm", df=1, the.mean=0, the.sd=1,  xlab = deparse(substitute(x)),
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






