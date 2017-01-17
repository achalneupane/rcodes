


plates.all<-read.delim("Hugo_summary4_50.txt",header=T,sep="\t",fill=TRUE)
fields<-read.delim("Hugo_field_summary4_50.txt",header=T,sep="\t",fill=TRUE)
the.expt<-"Hugo"

plates.all<-read.delim("Ivan_summary4_50.txt",header=T,sep="\t",fill=TRUE)
fields<-read.delim("Ivan_field_summary4_50.txt",header=T,sep="\t",fill=TRUE)
the.expt<-"Ivan"

plates.all<-read.delim("Joseph_summary.txt",header=T,sep="\t",fill=TRUE)
fields<-read.delim("Joseph_field_summary.txt",header=T,sep="\t",fill=TRUE)
the.expt<-"Joseph"

plates.all<-read.delim("Joseph_summary_NOTGREEN.txt",header=T,sep="\t",fill=TRUE)
fields<-read.delim("Joseph_field_summary_NOTGREEN.txt",header=T,sep="\t",fill=TRUE)
the.expt<-"Joseph"

the.score<-"Z.Low.Green"
the.score<-"Z.Mid.Green"
the.score<-"Z.High.Green" 
the.score<-"Z.all.Green"


################### PRELIMIANATY q-q
data.in<- (plates.all[is.finite(plates.all[,the.score]),the.score])
hist(data.in,col=col.array[i],breaks=100,main="Joseph screen RAW",xlab=the.score)

cut.high<- 2
cut.low<- -2
shapiro.test(data.in[data.in < cut.high & data.in > cut.low])
shapiro.test(data.in)
data.in<- scale(data.in )
hist(data.in,col=col.array[i],breaks=100,main="Joseph screen RAW",xlab=the.score)

labels<-as.character(paste(plates.all[is.finite(plates.all[,the.score]),"plate"],plates.all[is.finite(plates.all[,the.score]),"row"],plates.all[is.finite(plates.all[,the.score]),"col"],"-",plates.all[is.finite(plates.all[,the.score]),"Symbol"],sep=""))

my.qq.plot(data.in,distribution="norm",col="blue",ylab="Observed score",xlim=c(-15,15), ylim=c(-15,25),main=paste("Joseph UNNORMALIZED  expt: ",the.score,"  With 95% confidence intervals",sep=" "),the.mean=mean(data.in),the.sd=sd(data.in))
savePlot(filename=paste(the.expt,"Unnormalized Q-Q Plot for",the.score,"jpeg",sep="."),type="jpeg")

qq<- qq.data(data.in,distribution="norm",the.mean=mean(data.in),the.sd=sd(data.in),plot.it=FALSE)
selected.data<-identify(qq$x,qq$y,labels=labels[qq$ord],col="red",cex=0.75,atpen='TRUE')
selected.data<-identify(qq$x,qq$y,labels=as.character(round(data.in[qq$ord],3)),col="forestgreen",cex=1.0,atpen='TRUE')
points(qq$x,qq$y,col="magenta",pch=21)

#####################Normalization ###############
### do q-q plot to get estimates of cut-ss to apply
#plates.all[!is.finite(plates.all[,the.score]),the.score]<-NA
plate.list<-levels(as.factor(plates.all[,"plate"]))
cut.RG<-50
cut.low<- -14 # for  well estimates
cut.high<- 2 # for  well estimates


cut.low<- -20 # for Not Green well estimates
cut.high<- 2 # for  NOt Green  well estimates

samples<-plates.all[,"RG"] > cut.RG & plates.all[,the.score]>cut.low & plates.all[,the.score] < cut.high

plates.all[,"Z.High.Green"] > cut.low & plates.all[,"Z.High.Green"] < cut.high

 centers.all<-tapply(plates.all[samples,"Z.all.Green"],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE))
 centers.Low<-tapply(plates.all[samples,"Z.Low.Green"],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE))
 centers.Mid<-tapply(plates.all[samples,"Z.Mid.Green"],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE))
 centers.High<-tapply(plates.all[samples,"Z.High.Green"],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE))

tapply(plates.all[samples,the.score],plates.all[samples,"plate"],length)
### I think I should median center and then and then get scale  
#scale(plates.all[plates.all[,"plate"]==plate.list[i],the.score]  
hist(plates.all[samples,the.score],breaks=100,main="Joseph screen CHECK",xlab=the.score)
shapiro.test(plates.all[samples,the.score])


Kolmogorov-Smirnov test, taken from manual page ("?ks.test")

x <- rnorm(50)   # 50 Gaussian random numbers
y <- runif(30)   # 30 uniform random numbers
# Do x and y come from the same distribution?
ks.test(x, y)
######################### for scaled values
plates.all.ori<-plates.all  # hugo<-hugo.ori 
for(i in 1:length(plate.list)){
plates.all[plates.all[,"plate"]==plate.list[i],"Z.all.Green"] <-plates.all[plates.all[,"plate"]==plate.list[i],"Z.all.Green"]-centers.all[plate.list[i]]
plates.all[plates.all[,"plate"]==plate.list[i],"Z.Low.Green"] <-plates.all[plates.all[,"plate"]==plate.list[i],"Z.Low.Green"]-centers.Low[plate.list[i]]
plates.all[plates.all[,"plate"]==plate.list[i],"Z.Mid.Green"] <-plates.all[plates.all[,"plate"]==plate.list[i],"Z.Mid.Green"]-centers.Mid[plate.list[i]]
plates.all[plates.all[,"plate"]==plate.list[i],"Z.High.Green"] <-plates.all[plates.all[,"plate"]==plate.list[i],"Z.High.Green"]-centers.High[plate.list[i]]
}

      
scale.all<-tapply(plates.all[samples,"Z.all.Green"],plates.all[samples,"plate"],function(x) sqrt(sum(x[!is.na(x)]^2))/(length(x[!is.na(x)])-1)  )
scale.Low<-tapply(plates.all[samples,"Z.Low.Green"],plates.all[samples,"plate"],function(x) sqrt(sum(x[!is.na(x)]^2))/(length(x[!is.na(x)])-1)  )
scale.Mid<-tapply(plates.all[samples,"Z.Mid.Green"],plates.all[samples,"plate"],function(x) sqrt(sum(x[!is.na(x)]^2))/(length(x[!is.na(x)])-1)   )
scale.High<-tapply(plates.all[samples,"Z.High.Green"],plates.all[samples,"plate"],function(x) sqrt(sum(x[!is.na(x)]^2))/(length(x[!is.na(x)])-1)  )

for(i in 1:length(plate.list)){
plates.all[plates.all[,"plate"]==plate.list[i],"Z.all.Green"] <-plates.all[plates.all[,"plate"]==plate.list[i],"Z.all.Green"]*scale.all[plate.list[i]]
plates.all[plates.all[,"plate"]==plate.list[i],"Z.Low.Green"] <-plates.all[plates.all[,"plate"]==plate.list[i],"Z.Low.Green"]*scale.Low[plate.list[i]]
plates.all[plates.all[,"plate"]==plate.list[i],"Z.Mid.Green"] <-plates.all[plates.all[,"plate"]==plate.list[i],"Z.Mid.Green"]*scale.Mid[plate.list[i]]
plates.all[plates.all[,"plate"]==plate.list[i],"Z.High.Green"] <-plates.all[plates.all[,"plate"]==plate.list[i],"Z.High.Green"]*scale.High[plate.list[i]]
}      
hist(plates.all[samples,the.score],breaks=100,main="Joseph screen CHECK",xlab=the.score)      
shapiro.test(as.numeric(plates.all[samples  ,the.score]))
# p-value = 6.49e-06 unnormalized but filtered9
# p-value = 0.0001199 unnormalized but filtered 12 and 0
#  p-value = 0.001775  centered
# p-value = 5.178e-05 centered and scaled

## NOT GREEN
#p-value = 0.06337 centered
#p-value = 0.01608 centered and scaled
hist(plates.all[samples,the.score],col=col.array[i],breaks=100,main="Joseph screen z-Normalized",xlab=the.score)      
hist(plates.all[,the.score],col=col.array[i],breaks=100,main="Joseph screen z-Normalized",xlab=the.score)
savePlot("Z-all-Green Normalized Not Green.jpeg",type="jpeg")
i<-18
hist(plates.all[plates.all[,"plate"]==plate.list[i] & samples ,the.score],breaks=10,col=col.array[i])
write.table(plates.all,"Joseph_Norm_summary-NotGreen.txt",row.names=FALSE,sep="\t")

the.score<-"Z.Low.Green"
the.score<-"Z.Mid.Green"
the.score<-"Z.High.Green" 
the.score<-"Z.all.Green"
      
data.in<- (plates.all[is.finite(plates.all[,the.score]),the.score])

labels<-as.character(paste(plates.all[is.finite(plates.all[,the.score]),"plate"],plates.all[is.finite(plates.all[,the.score]),"row"],plates.all[is.finite(plates.all[,the.score]),"col"],"-",plates.all[is.finite(plates.all[,the.score]),"Symbol"],sep=""))
plate.labels<-as.character(plates.all[is.finite(plates.all[,the.score]),"plate"])

my.qq.plot(data.in,dist="norm",col="blue",ylab="Observed score",xlim=c(-10,10), ylim=c(-15,15),main=paste("Joseph NORMALIZED expt: ",the.score,"  With 95% confidence intervals",sep=" "),the.mean=mean(data.in),the.sd=sd(data.in))
savePlot(filename=paste(the.expt,"Q-Q Plot for NotGreen",the.score,"jpeg",sep="."),type="jpeg")

my.qq.plot(data.in,dist="norm",col="blue",ylab="Observed score",xlim=c(5,10), ylim=c(5,15),main=paste("UPPER ARM Plates.All expt: ",the.score,"  With 95% confidence intervals",sep=" "),the.mean=mean(data.in),the.sd=sd(data.in))
savePlot(filename=paste(the.expt,"Upper arm Q-Q Plot for",the.score,"png",sep="."),type="png")

my.qq.plot(data.in,dist="norm",col="blue",ylab="Observed score",xlim=c(-4,-1), ylim=c(-10,-5),main=paste("LOWER ARM Plates.All expt: ",the.score,"  With 95% confidence intervals",sep=" "))
savePlot(filename=paste(the.expt,"Lower arm Q-Q Plot for",the.score,"png",sep="."),type="png")


qq<- qq.data(data.in,distribution="norm",the.mean=mean(data.in),the.sd=sd(data.in),plot.it=FALSE)
selected.data<-identify(qq$x,qq$y,labels=labels[qq$ord],col="red",cex=0.75,atpen='TRUE')
selected.data<-identify(qq$x,qq$y,labels=as.character(round(data.in[qq$ord],2)),col="forestgreen",cex=1.0,atpen='TRUE')

col.array<-rainbow(length(plate.list))
points(qq$x,qq$y,col="magenta",pch=plate.labels[qq$ord],cex=1.5)
points(qq$x,qq$y,col=col.array[as.numeric(plate.labels[qq$ord])],pch=23,cex=1.5)

      
      
                                                                                       
########################## scale uses mean centered and  root-mean-square scaling
data.in<- (plates.all[is.finite(plates.all[,"Trans_."]),])
samples<-data.in[,"Symbol"] != "MOCK"


col.array<-rainbow(length(plate.list))
par(mar=c(5.5,5.5,2.5,5.1),mgp=c(3,1,0))
hist(data.in[samples,"Trans_."],breaks=75,col="green",main="Percent Transduction in Joseph",xlab="Percent Tranduction",cex=2,,font.lab=2.5,cex.lab=2.5,cex.main=2.0,font.axis=2,cex.axis=2.0)
savePlot("percent Tramsduction.jpeg",type="jpeg")
###################### NOrmalization EXPLORE ###############################
plate.list<-unique(plates.all[,"plate"])
tapply(plates.all[,the.score],plates.all[,"plate"],function(x) median(x,na.rm=TRUE))
tapply(plates.all[samples,the.score],plates.all[samples,"plate"],length)

plates.all[,"plate"],function(x) median(x,na.rm=TRUE))
       1        2        3        4        5        6        7        8 
 -7.0500  -7.1830  -7.1960  -5.6720  -4.5535  -6.7810  -8.4020  -9.3090 
       9       10       11       12       13       14       15       17 
 -8.7600 -11.1230  -9.5420  -8.6355  -8.7475  -6.8040  -7.8035  -7.8150 
      18       21 
 -7.7610  -4.1255
  ########### DO a q-q plot outlyers are >2 and < -14 ## core ones with enough RG and are not hits
samples<-plates.all[,"RG"] > 50 & plates.all[,the.score]>-14 & plates.all[,the.score]<2
samples<-plates.all[,"RG"] > 50 & plates.all[,the.score]>- 14 & plates.all[,the.score] < 2 & plates.all[,"Z.High.Green"] > -14 & plates.all[,"Z.High.Green"] < 2

samples<-plates.all[,"RG"] > 50 & plates.all[,the.score]>- 12 & plates.all[,the.score] < 0 & plates.all[,"Z.High.Green"] > -12 & plates.all[,"Z.High.Green"] < 0 ### about same as before 
  tapply(plates.all[samples,the.score],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE)) ## Z-high Green & Z-all- Green
        1        2        3        4        5        6        7        8 ##-14, 2
 -7.1970  -6.6010  -7.5200  -8.3630  -4.8495  -8.4840  -8.3220  -9.6680 
       9       10       11       12       13       14       15       17 
 -9.0170 -10.2370  -8.5510  -8.8410  -7.8870  -6.6410  -7.0100  -6.6010 
      18       21 
 -6.7685  -3.4800

       1        2        3        4        5        6        7        8 ##-12 0
 -7.3235  -6.3605  -7.1740  -8.2830  -4.8495  -7.4970  -8.1310  -9.1660 
       9       10       11       12       13       14       15       17 
 -8.9980 -10.1720  -8.6345  -8.7610  -7.8800  -6.7385  -7.2270  -6.5695 
      18       21 
 -6.4480  -3.4790 

hist(plates.all[samples,the.score],breaks=105,col=col.array[i])
i<-1
hist(plates.all[plates.all[,"plate"]==plate.list[i] & samples ,the.score],breaks=10,col=col.array[i])



tapply(plates.all[samples,"plate"],plates.all[samples,"plate"],length)
1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 17 18 21 ##-14, 2
37 54 65 38 36 46 65 58 62 63 52 57 73 52 68 60 45 38

 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 17 18 21 ##-12 0
32 52 57 32 36 34 55 45 56 56 48 47 65 50 52 56 42 31

 tapply(plates.all[samples,the.score],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE)) ## Z-high Green
      1       2       3       4       5       6       7       8       9      10 
-5.1730 -6.0510 -5.8850 -4.8860 -4.8620 -5.6515 -4.6020 -5.8645 -5.6915 -7.4310 
     11      12      13      14      15      17      18      21 
-6.7655 -5.5110 -5.4930 -5.1450 -5.2320 -6.3365 -6.2400 -4.1575
 tapply(plates.all[samples,the.score],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE)) ## Z-high Green
# Looks  more consistent than controls alone
shapiro.test(as.numeric(plates.all[samples ,the.score]))
# p-value = 6.49e-06 unnormalized but filtered9
# p-value = 0.0001199 unnormalized but filtered 12 and 0


samples<-plates.all

hist(plates.all[samples,the.score],breaks=105,col=col.array[i])
i<-21
hist(plates.all[plates.all[,"plate"]==plate.list[i] & samples ,the.score],breaks=10,col=col.array[i])


       1        2        3        4        5        6        7        8 
-10.4550 -10.5920 -10.8580  -9.8245  -6.1545 -10.9340 -12.5820 -13.4250 
       9       10       11       12       13       14       15       17 
-11.7510 -14.0260 -13.3950 -11.9935 -11.9240 -10.0430  -9.8705 -10.4390 
      18       21 
 -9.8340  -5.0850





samples<-plates.all[,"Symbol"]=="PLV101"& plates.all[,"RG"] > 50
tapply(plates.all[samples,the.score],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE))
> tapply(plates.all[samples,the.score],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE))
      1       2       3       4       5       6       7       8       9      10 
 -1.353  -7.063  -8.725  -9.748  -3.649 -10.559  -8.305 -10.814  -9.845 -10.903 
     11      12      13      14      15      17      18      21 
-11.496 -10.527  -8.706 -10.095  -9.523  -6.283  -5.715  -2.595
######### using only with well where > 50 RG cells
samples<-plates.all[,"RG"] > 50
tapply(plates.all[samples,the.score],plates.all[samples,"plate"],length)
 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 17 18 21 
48 57 73 59 36 62 76 72 64 64 58 64 79 65 74 63 45 42


tapply(plates.all[samples,the.score],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE))
       1        2        3        4        5        6        7        8 
 -7.3015  -8.0180  -7.7600  -6.7890  -5.3785  -7.7780  -8.5310  -9.6170 
       9       10       11       12       13       14       15       17 
 -9.0940 -11.1420 -10.0410  -9.1885  -8.6740  -6.8040  -7.8035  -7.8200 
      18       21 
 -7.9710  -4.4605
######### using only with well where > 50 RG cells
samples<-plates.all[,"RG"] > 100
tapply(plates.all[samples,the.score],plates.all[samples,"plate"],length)
tapply(plates.all[samples,the.score],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE))
 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 17 18 21 
31 40 48 42 13 51 58 47 43 36 34 50 64 29 58 45 25 30 
> tapply(plates.all[samples,the.score],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE))
       1        2        3        4        5        6        7        8 
 -5.5560  -7.6930  -7.7370  -6.6545  -6.8100  -6.1230  -8.5310  -9.6080 
       9       10       11       12       13       14       15       17 
 -9.4430 -11.3130 -10.3440  -9.1115  -8.6900  -6.8040  -7.4195  -7.8150 
      18       21 
 -7.0440  -4.3720



tapply(plates.all[samples,"X.R"],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE))
     1      2      3      4      5      6      7      8      9     10     11 
 4.538 11.360  9.836 16.384 10.999 14.837 15.217 16.445 15.423 16.149 15.493 
    12     13     14     15     17     18     21 
16.837 16.844 15.252 16.473 16.300 15.449 16.710

tapply(plates.all[samples,"Total"],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE))
   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   17 
 673 1916 3130 4458 6684 5338 4800 6121 6295 5524 5374 5534 7378 4130 4252 6074 
  18   21 
7541 9075


> tapply(plates.all[,"Total"],plates.all[,"plate"],function(x) sd(x,na.rm=TRUE))
       1        2        3        4        5        6        7        8 
2582.250 2746.694 2402.255 3860.152 2896.970 3683.041 2871.537 2312.270 
       9       10       11       12       13       14       15       17 
2446.584 2701.131 2546.627 2482.391 2223.013 2598.105 2808.478 2765.582 
      18       21 
2950.014 4013.994 
> tapply(plates.all[samples,"Total"],plates.all[samples,"plate"],function(x) sd(x,na.rm=TRUE))
         1          2          3          4          5          6          7 
  96.55724  509.49063  505.93893  974.14116 1813.26576  704.85483  373.14921 
         8          9         10         11         12         13         14 
 277.82249  528.39979  819.57692 1474.30696  944.02348  687.97844 1077.82389 
        15         17         18         21 
 928.13271 1568.11553 1134.63496  672.82712 



tapply(plates.all[samples,"X.R"],plates.all[samples,"plate"],function(x) sd(x,na.rm=TRUE))
        1         2         3         4         5         6         7         8 
1.6273582 1.7592144 3.3047083 2.1522375 0.8068856 0.5433461 0.1205156 0.4873729 
        9        10        11        12        13        14        15        17 
0.7561484 1.0882621 2.2443532 0.9438482 2.9162895 1.9935346 0.8064498 2.7115986 
       18        21 
2.1723720 0.5104913

check=sqrt(as.numeric(plates.all[,"Total"])*(plates.all[,"X.R"]/100)*(1-(plates.all[,"X.R"]/100)))

tapply(check[samples],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE))
tapply(check,plates.all[,"plate"],function(x) median(x,na.rm=TRUE))

> tapply(check[samples],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE))
        1         2         3         4         5         6         7         8 
 5.174183 14.737680 15.161541 24.712956 26.081658 25.814616 24.782695 29.251499 
        9        10        11        12        13        14        15        17 
28.606247 27.505293 23.908450 27.836652 32.188239 23.104850 24.905804 26.531440 
       18        21 
33.294995 35.649349 
> tapply(check,plates.all[,"plate"],function(x) median(x,na.rm=TRUE))
       1        2        3        4        5        6        7        8 
21.27909 24.97753 25.13392 30.91770 24.00315 27.48460 27.66900 23.73302 
       9       10       11       12       13       14       15       17 
28.11241 26.26289 23.71153 28.15789 26.08540 21.61618 25.85133 24.31030 
      18       21 
27.78374 37.71699


tapply(plates.all[samples,"Total"],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE))
> tapply(plates.all[samples,"Total"],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE))
      1       2       3       4       5       6       7       8       9      10 
 4261.0  5385.0  5648.0  7690.5  6682.0  6021.0  4645.0  4580.0  5955.0  4933.5 
     11      12      13      14      15      17      18      21 
 4449.0  5990.0  4864.0  3700.0  4976.0  4877.0  6741.0 11296.0
tapply(plates.all[samples,"Total"],plates.all[samples,"plate"],function(x) sd(x,na.rm=TRUE))
       1        2        3        4        5        6        7        8 
1805.579 2033.554 2189.243 2919.187 1643.149 2983.470 2175.194 2174.159 
       9       10       11       12       13       14       15       17 
2077.419 2810.453 1514.201 1958.105 1918.507 1707.228 2196.840 2303.147 
      18       21 
2236.619 2774.365

> tapply(plates.all[samples,"Trans_."],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE))
      1       2       3       4       5       6       7       8       9      10 
37.2230 36.3660 32.1720 22.8630 19.5995 29.5970 37.3170 36.0160 31.3340 29.9200 
     11      12      13      14      15      17      18      21 
35.6630 29.5650 37.1420 38.4940 34.9050 36.4010 26.1010 11.6940 
> tapply(plates.all[samples,"Trans_."],plates.all[samples,"plate"],function(x) mean(x,na.rm=TRUE))
       1        2        3        4        5        6        7        8 
37.10824 37.22642 33.74609 23.74648 21.82472 29.86945 37.37507 38.71374 
       9       10       11       12       13       14       15       17 
30.52114 31.36492 37.78525 30.41488 36.92820 38.31963 33.63130 35.83087 
      18       21 
25.79327 16.14338 
tapply(plates.all[samples,"Trans_."],plates.all[samples,"plate"],function(x) sd(x,na.rm=TRUE))
       1         2         3         4         5         6         7         8 
14.843948 11.754386 13.458241 10.672068 10.463499 11.870782 12.543951 12.357003 
        9        10        11        12        13        14        15        17 
12.194360 12.802979 16.266918 11.826850 11.195500 12.945871 12.149010 14.678545 
       18        21 
 8.967137  9.67974
 tapply(plates.all[samples,"X.R"],plates.all[samples,"plate"],function(x) median(x,na.rm=TRUE))
      1       2       3       4       5       6       7       8       9      10 
20.0000 18.8040 18.1920 19.5500 13.6310 20.1860 24.4340 21.8800 19.0320 20.9900 
     11      12      13      14      15      17      18      21 
21.7860 19.7600 20.8990 17.7420 19.6935 20.7380 16.6240 13.6700 
> tapply(plates.all[samples,"X.R"],plates.all[samples,"plate"],function(x) mean(x,na.rm=TRUE))
       1        2        3        4        5        6        7        8 
20.37407 18.72960 17.58305 18.90886 13.76708 20.04388 24.67332 22.75095 
       9       10       11       12       13       14       15       17 
18.73702 21.64362 22.14965 20.17174 21.33893 19.17007 20.21063 20.08997 
      18       21 
16.82947 13.93073 
> tapply(plates.all[samples,"X.R"],plates.all[samples,"plate"],function(x) sd(x,na.rm=TRUE))
       1        2        3        4        5        6        7        8 
5.195297 4.155517 5.519607 3.275432 3.509765 5.188146 4.430769 3.728062 
       9       10       11       12       13       14       15       17 
3.221360 4.104403 4.813066 4.748213 4.075727 4.983271 4.427129 5.615670 
      18       21 
4.493462 3.406489
tapply(plates.all[samples,"RG"],plates.all[samples,"plate"],function(x) sd(x,na.rm=TRUE))
> tapply(plates.all[samples,"RG"],plates.all[samples,"plate"],function(x) sd(x,na.rm=TRUE))
        1         2         3         4         5         6         7         8 
 65.12918 145.38636  56.81956 103.00155  56.48351  84.08123 125.79473  87.37889 
        9        10        11        12        13        14        15        17 
 75.33735  50.09736  68.93327  73.78803  98.10328  56.39693  87.87238  81.34496 
       18        21 
 70.97544 157.6750





hist(plates.all[,the.score],breaks=100,xlab=the.score,main="Unnormalized plate values",col="blue")
col.array<-rainbow(length(plate.list))
plot(c(-20,15),c(0,20),pch="")
for(i in 1:length(plate.list)){
#for(i in 1:5){
a.hist<-hist(plates.all[plates.all[,"plate"]==plate.list[i],the.score],breaks=25,plot=FALSE)
lines(a.hist$mids,a.hist$counts,col=col.array[i],lwd=2,type="b")
}

shapiro.test(as.numeric(plates.all[plates.all[,the.score]> -15 & plates.all[,the.score]< 2.0 ,the.score]))
i<-17
hist(plates.all[plates.all[,"plate"]==plate.list[i],the.score],breaks=15,col=col.array[i])
samples<- plates.all[,"RG"] > 50
hist(plates.all[plates.all[,"plate"]==plate.list[i] & plates.all[,"RG"] > 50 ,the.score],breaks=15,col=col.array[i])
barplot(sort(plates.all[plates.all[,"plate"]==plate.list[i] & plates.all[,"RG"] > 50 ,the.score]))

low<-min(plates.all[is.finite(plates.all[,the.score]),the.score])
high<-max(plates.all[is.finite(plates.all[,the.score]),the.score])
high<-15
low<--10



  


for(i in 1:length(plate.list)){
i<-plate.list[i]

plates<-plates.all[is.finite(plates.all[,the.score]),]
plates<-plates[plates[,"plate"]==i,]

col.array<-rep("blue",times=dim(plates)[1])
col.labels<-rep("black",times=dim(plates)[1])
col.labels[plates[,"Z.all.Green"] > 0]<-"red"



col.array[plates[,"Symbol"]=="PLV101"]<-"green"
order.score<-order(plates[,the.score])

the.plot<-barplot(plates[order.score,the.score],names.arg="",col=col.array[order.score],las=2,offset=0,ylab=the.score, ylim=c(low,high),main=paste("Expt: ",the.expt,",  Plate: ", i,",  Using Score:", the.score))
#the.plot<-barplot(plates[order.score,"Prob"],names.arg=plates[order.score,"Symbol"],col=col.array[order.score],las=2,offset=0,bg="grey")
text(the.plot[,1],0.5,labels=paste(plates[order.score,"plate"],plates[order.score,"row"],plates[order.score,"col"],"-",plates[order.score,"Symbol"],sep=" "),srt=90,adj=0,col=col.labels[order.score],cex=0.8) #at=the.plot[,1],
savePlot(filename=paste(the.expt,i,the.score,"png",sep="."),type="png")
}
#mtext(plates[order.score,"Symbol"],side=1,las=2,at=the.plot[,1],adj=-0.9) #at=the.plot[,1],

######################### all plates
###############################################

the.score<-"Z.Low.Green"
the.score<-"Z.Mid.Green"
the.score<-"Z.High.Green" 
the.score<-"Z.all.Green"

######################### for scaled values
plates.all.ori<-plates.all  # hugo<-hugo.ori 
for(i in 1:clength(plate.list)){
plates.all[plates.all[,"plate"]==plate.list[i],the.score] <-scale(plates.all[plates.all[,"plate"]==plate.list[i],the.score]
  )}

plate.list<-c("1","2","3")
for(i in 1:length(plate.list)){
plates.all[plates.all[,"plate"]==plate.list[i],the.score] <-scale(plates.all[plates.all[,"plate"]==plate.list[i],the.score]
  )}
##########################  



low<-min(plates.all[is.finite(plates.all[,the.score]),the.score])
high<-max(plates.all[is.finite(plates.all[,the.score]),the.score])
high<-15
low<--10

plates<-plates.all[is.finite(plates.all[,the.score]),]
col.array<-rep("blue",times=dim(plates)[1])
col.labels<-rep("black",times=dim(plates)[1])
col.labels[plates[,"Z.all.Green"] > 0]<-"red"
col.labels[plates[,the.score] > 0]<-"red"
col.array[plates[,"Symbol"]=="PLV101"]<-"green"
order.score<-order(plates[,the.score])

the.plot<-barplot(plates[order.score,the.score],names.arg="",col=col.array[order.score],las=2,offset=0,ylab=the.score, ylim=c(low,high),main=paste("Expt: ",the.expt,",  Plate: ", i,",  Using Score:", the.score),space=0.1,beside=TRUE)
#the.plot<-barplot(plates[order.score,"Prob"],names.arg=plates[order.score,"Symbol"],col=col.array[order.score],las=2,offset=0,bg="grey")
text(the.plot[,1],0.5,labels=paste(plates[order.score,"plate"],plates[order.score,"row"],plates[order.score,"col"],"-",plates[order.score,"Symbol"],sep=" "),srt=90,adj=0,col=col.labels[order.score],cex=0.4) #at=the.plot[,1],
savePlot(filename=paste(the.expt,"ALL Plates SCALED",the.score,"png",sep="."),type="png")

hugo.scaled<-plates.all
save.image("Hugo scaled.RData")

############################### NORMAL?:
#### test normality


the.score<-"Z.Low.Green"
the.score<-"Z.Mid.Green"
the.score<-"Z.High.Green" 
the.score<-"Z.all.Green"

shapiro.test(as.numeric(hugo.scaled[,the.score]))
data:  as.numeric(hugo.scaled[, the.score]) 
W = 0.8478, p-value < 2.2e-16

hist(hugo.scaled[,the.score],breaks=50)
shapiro.test(as.numeric(hugo.scaled[hugo.scaled[,the.score]> -1 & hugo.scaled[,the.score]< 1.0 ,the.score]))
# scaling sux a bit for plate 4 
################# 
shapiro.test(as.numeric(hugo[,the.score]))

data:  as.numeric(hugo[, the.score]) 
W = 0.8478, p-value < 2.2e-16
 hist(hugo.scaled[,the.score],breaks=50)### plot and check for outlyers
shapiro.test(as.numeric(hugo[hugo[,the.score]> -11 & hugo[,the.score]< 2.5 ,the.score]))
	Shapiro-Wilk normality test

data:  as.numeric(hugo[hugo[, the.score] > -11 & hugo[, the.score] >      2.5, the.score]) 
W = 0.8671, p-value = 0.2865
########## score is IS normall distributed  when trimmed:

 hist(hugo[,the.score],breaks=50)

shapiro.test(rnorm(1000, mean = 0, sd = 3))     ### if significant then did NOT come from a normal distribution
                                                  ## This is said in Royston (1995) to be adequate for p.value < 0.1
        Shapiro-Wilk normality test

data:  rnorm(100, mean = 5, sd = 3) 
W = 0.9929, p-value = 0.8837

shapiro.test(runif(100, min = 2, max = 4))

        Shapiro-Wilk normality test

data:  runif(100, min = 2, max = 4)
W = 0.9413, p-value = 0.0002307
hist(runif(100, min = 2, max = 4))
hist(rnorm(1000, mean = 0, sd=2))


 #dim(all.prob)<-dim(red.c)
data.in<- (hugo[is.finite(hugo[,the.score]),the.score])
null<-rnorm(length(hugo[,the.score]),mean=mean(data.in),sd=sd(data.in))
qq<-  qqplot((null),(hugo[,the.score]),plot.it=TRUE)
abline(0,1)
qq<-  qqplot((null),(summ[,"Prob"]),plot.it=TRUE)


data.in<- (plates.all[is.finite(plates.all[,the.score]),the.score])

labels<-as.character(paste(plates.all[is.finite(plates.all[,the.score]),"plate"],plates.all[is.finite(plates.all[,the.score]),"row"],plates.all[is.finite(plates.all[,the.score]),"col"],"-",plates.all[is.finite(plates.all[,the.score]),"Symbol"],sep=""))

my.qq.plot(data.in,dist="norm",col="blue",ylab="Observed score",xlim=c(-25,15), ylim=c(-20,15),main=paste("Joseph UNNORMALIZED  expt: ",the.score,"  With 95% confidence intervals",sep=" "),the.mean=mean(data.in),the.sd=sd(data.in))
savePlot(filename=paste(the.expt,"Q-Q Plot for",the.score,"png",sep="."),type="png")

my.qq.plot(data.in,dist="norm",col="blue",ylab="Observed score",xlim=c(1.5,3.5), ylim=c(-2,20),main=paste("UPPER ARM Plates.All expt: ",the.score,"  With 95% confidence intervals",sep=" "))
savePlot(filename=paste(the.expt,"Upper arm Q-Q Plot for",the.score,"png",sep="."),type="png")

my.qq.plot(data.in,dist="norm",col="blue",ylab="Observed score",xlim=c(-4,-1), ylim=c(-10,-5),main=paste("LOWER ARM Plates.All expt: ",the.score,"  With 95% confidence intervals",sep=" "))
savePlot(filename=paste(the.expt,"Lower arm Q-Q Plot for",the.score,"png",sep="."),type="png")


qq<- qq.data(data.in,distribution="norm",the.mean=mean(data.in),the.sd=sd(data.in),plot.it=FALSE)
selected.data<-identify(qq$x,qq$y,labels=labels[qq$ord],col="red",cex=0.75,atpen='TRUE')
selected.data<-identify(qq$x,qq$y,labels=as.character(data.in[qq$ord]),col="forestgreen",cex=1.0,atpen='TRUE')
points(qq$x,qq$y,col="magenta",pch=21)





my.qq.plot(-log10(summ[,"Prob"]),dist="unif",col="blue",ylab="Observed chi-squared value",xlim=c(0,1),ylim=c(0,30))
xlab="Expected chi-squared value",main="Q-Q plot: 1905 Cases 5790 Controls",cex=1,xlim=c(0,22),ylim=c(0,30),font.lab=2,font.axis=2)











############################ compare expts:
setwd("/media/Bioinform-D/Research/Cellomics/Hugo screen/")
plates.all<-read.delim("Hugo_summary4_50.txt",header=T,sep="\t",fill=TRUE)
fields<-read.delim("Hugo_field_summary4_50.txt",header=T,sep="\t",fill=TRUE)
hugo<-plates.all

setwd("/media/Bioinform-D/Research/Cellomics/Ivan screen/")
plates.all<-read.delim("Ivan_summary4_50.txt",header=T,sep="\t",fill=TRUE)
fields<-read.delim("Ivan_field_summary4_50.txt",header=T,sep="\t",fill=TRUE)


ivan<-plates.all


> colnames(hugo)
 [1] "plate"          "row"            "col"            "Symbol"        
 [5] "Z.all.Green"    "Z.Low.Green"    "Z.Mid.Green"    "Z.High.Green"  
 [9] "RG"             "RG_expected"    "X.R"            "R"             
[13] "RG_Low"         "RG_Mid"         "RG_High"        "G"             
[17] "Num.cells"      "score"          "RG.G"           "RnG.nG"        
[21] "RnG"            "nR"             "nG"             "gt4N"          
[25] "Total"          "Trans_."        "X.cells_field." "sd_cells_field"

[29] "Description"    "GenBank_Accn"   "dbEST_Id"       "MGC.Accession"
var<-"X.R"
well.labels<-apply(hugo,1,function(x) as.character(paste(x[1],x[2],x[3],sep="")))

ramp<-colorRamp(c("red", "green"))
well.col <- rgb( ramp(seq(0, 1, length = dim(hugo)[1])), max = 255)

 well.col<-rainbow(dim(hugo)[1])


the.var<-"Trans_."
the.var<-"RG"
the.var<-"Total"
the.var<-"X.cells_field."
the.var<-"sd_cells_field"
the.var<-"RG_expected"
the.var<-"row"
the.var<-"col"
the.var<-"RnG.nG"
the.var<-"Z.High.Green"
the.var<-"Z.Mid.Green"
the.var<-"Z.Low.Green"
the.var<-"Z.all.Green" 
the.var<-"R"
the.order<-order(as.numeric(hugo[,the.var]))

plot(hugo[the.order,"X.R"],ivan[the.order,"X.R"],pch="",cex=0.5,xlab="hugo %R",ylab="ivan %R",main=paste("Comparion of Hugo and Ivan aginst: ",the.var, " : red(low)->green(high)",sep= " "))
text(hugo[the.order,"X.R"],ivan[the.order,"X.R"],labels=as.character(well.labels)[the.order],cex=0.7,col=well.col)
controls<-hugo[,"Symbol"]=="MOCK"
plv101<-hugo[,"Symbol"]=="PLV101"
points(hugo[controls,"X.R"],ivan[controls,"X.R"],pch=21,cex=1.9)
points(hugo[plv101,"X.R"],ivan[plv101,"X.R"],pch=22,cex=2.5)
legend(0,14,c("MOCK","PLV101"),pch=c(21,22))
abline(c(0,1))
savePlot(filename=paste("hugo and Ivan % R and ",the.var,"png",sep="."),type="png")

the.model<-lmrob(TotalIntenCh3~ObjectTotalIntenCh1,data=test)

 cbind(hugo[,1:3],hugo[,"X.R"],ivan[,"X.R"],as.character(well.labels))
 test<-cbind(hugo[the.order,1:3],hugo[the.order,"X.R"],ivan[the.order,"X.R"],(well.labels)[the.order],well.col[the.order])

 test<-cbind(hugo[,1:3],hugo[,"X.R"],ivan[,"X.R"],(well.labels),well.col[the.order])
colnames(test)<-c(1:7)
barplot(test[,4],col=test[,7])


barplot(hugo[,"X.R"],col=(well.col[the.order]))

barplot(hugo[,"X.R"])
tester<-well.col[the.order]
barplot(hugo[,"X.R"],col=tester)

barplot(hugo[the.order,"X.R"],col=well.col)

barplot(rep(1,length(hugo[,"X.R"])),col=the.order[well.col])

barplot(rep(1,length(hugo[,"X.R"])),col=well.col)

barplot(x[the.order],col=well.col)
barplot(x,col=well.col[the.order])
barplot(x[the.order],col=well.col)
barplot(x,col=well.col[the.order])


hist(hugo[controls | plv101,"X.R"])
median(hugo[controls | plv101,"X.R"])
mean(hugo[controls | plv101,"X.R"])
sd(hugo[controls | plv101,"X.R"])

#######################################################################
#####################################################################
#######################################################################
#####################################################################
#################################################### plate avergaes
the.score <- "Z.all.Green"
the.score <- "Z.Low.Green"
the.score <- "Z.Mid.Green"
the.score <- "Z.Total.Green"

hugo[controls | plv101,the.score]
iplate<-1
iplate2<-1
median(hugo[hugo[,"plate"]==iplate,the.score],na.rm=TRUE)
mean(hugo[hugo[,"plate"]==iplate,the.score],na.rm=TRUE)
sd(hugo[hugo[,"plate"]==iplate,the.score],na.rm=TRUE)

plot(hugo[hugo[,"plate"]==iplate,the.score],hugo[hugo[,"plate"]==iplate2,the.score])
hist(hugo[hugo[,"plate"]==iplate,the.score])
test<-scale(hugo[hugo[,"plate"]==iplate,the.score])

#######################################################################
#####################################################################
#######################################################################
#####################################################################
#################################################### Z-factor score

)       1. \text{Zfactor} = 1 - {3 \times (\sigma_p + \sigma_n) \over | \mu_p - \mu_n |}

An alternative but equivalent definition of Z-factor is calculated from the Sum of Standard Deviations (SSD) divided by the range of the assay (R):

       1. SSD = σp + σn
       2. R = | μp − μn |
       3. \text{Zfactor} = 1 - 3 \times {\text{SSD} \over R}

The following interpretations for the Z-factor were taken from Zhang, et. al. 1999:
Z-factor 	Interpretation
1.0 	Ideal. This is approached when you have a huge dynamic range with tiny standard deviations. Z-factors can never actually equal 1.0 and can certainly never be greater than 1.0.

between 0.5 and 1.0 	An excellent assay. Note that if σp = σn, 0.5 is equivalent to a separation of 12 standard deviations between μp and μn.

between 0 and 0.5 	A marginal assay.
less than 0 	The signal from the positive and negative controls overlap, making the assay essentially useless for screening purposes.

The Z-factor is no

################## controls
> median(hugo[controls | plv101,"X.R"])
[1] 3.601
> mean(hugo[controls | plv101,"X.R"])
[1] 4.075042
> sd(hugo[controls | plv101,"X.R"])
[1] 2.416469
################## positives
CCNE1 CD48 PLK3 1H7-CyclD2
pos<-hugo[,"Symbol"]=="CCNE1"
pos<-hugo[,"Symbol"]=="1H7-CyclD2"

the.score <- "Z.all.Green"
the.score <- "Z.Low.Green"
the.score <- "Z.Mid.Green"
the.score <- "Z.Total.Green"

hugo[controls | plv101,the.score]

median(hugo[controls | plv101,the.score],na.rm=TRUE)
neg.median <- mean(hugo[controls | plv101,the.score],na.rm=TRUE)
neg.mean <- mean(hugo[controls | plv101,the.score],na.rm=TRUE)
neg.sd <- sd(hugo[controls | plv101,the.score],na.rm=TRUE)
neg.median
 neg.mean
neg.sd

pos.median <- mean(hugo[pos,the.score],na.rm=TRUE)
pos.mean <- mean(hugo[pos,the.score],na.rm=TRUE)
pos.sd <- sd(hugo[pos,the.score],na.rm=TRUE)

pos.median
pos.mean
pos.sd

1-3*(pos.sd+neg.sd)/(abs(pos.mean-neg.mean))

########################################################################
#####################################################################



ramp<-colorRamp(c("red", "green"))
well.col <- rgb( ramp(seq(0, 1, length = length(x))), max = 255)
x<-c(8,6,7,6,10,4,3,1)
the.order<-order(x)
the.order

well.col<-c("darkred","red","darkorange1","orange","grey","lightgreen","green","darkgreen")
x<-c(8,6,7,6,10,4,3,1)
the.order<-order(x)
barplot(x[the.order],col=well.col)
barplot(x,col=well.col[the.order])
cbind(x,well.col[the.order],the.order,well.col)






 ##########################
my.qq.plot<-function (x, distribution = "norm", ylab = deparse(substitute(x)),
    xlab = paste(distribution, "quantiles"), main = NULL, las = par("las"), 
    envelope = 0.95, labels = FALSE, col = palette()[2], lwd = 2, the.mean=0,the.sd=1,
    pch = 1, cex = 1, line = c("quartiles", "robust", "none"),xlim=c(0,100),ylim=c(0,20),font.lab=2,font.axis=2,font.main=2,
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
        col = col, pch = pch,cex = cex,xlim=xlim,ylim=ylim,font.lab=font.lab,font.axis=font.axis,font.main=2)
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


