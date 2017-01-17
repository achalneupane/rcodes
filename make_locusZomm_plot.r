#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#########################################################DO PLOT DO PLOT


#colorss<-bluered(max(data)*4)
color = colorRampPalette(c("blue","cyan","green","red"),space="Lab",interpolate="linear")
colorss<-color(11)
colorss<-c(colorss,"black")                                       #  color = colorRampPalette(c("blue","red"))
temp2<-matrix(0:10)

## colorss<-c("blue","slategrey","green","orange","red","red") #,"purple")
## num.col<-length(colorss)
## colorss<-c(colorss,"purple")
## col.step=1/(num.col-1)
## temp2<-matrix(1:(num.col-1))
## rownames(temp2)<-temp2[,1]


 # for use in mon-floating scale bar

##  nf<-layout(matrix(c(1,2,3,4),4,byrow=TRUE),heights=c(1.0,4.5,1.25,1.6))
#layout.show(nf)
par.opts.start<-par(no.readonly=TRUE)
  ############# DO BIOMART QUERYIES FIRST

par(mar=c(0,5.5,2.1,5.1),mgp=c(3,1,0)) #c(bottom, left, top, right) ## left

data(chrom.bands)
pos<-0
width<-0.25
lwd<-2
plot(c(0,max(chrom.bands[chrom.bands[,"chr"]==the.chr,"bases.bot"])),c(-1.75*width,0),pch="",axes=FALSE,xlab="",ylab="",main=paste("Chromosome",the.chr,sep=" "),cex.main=2.0)
paint.chromosomeBP(the.chr, pos = pos, width = width, bands = "major") ## function defines below:
segments(0,-1.5*width,min(ann[,"start_position"]),-1*width,lwd=lwd,col="red")
segments(min(ann[,"start_position"]),-1*width,min(ann[,"start_position"]),pos,lwd=lwd,col="red")
segments(max(ann[,"end_position"]),-1*width,max(ann[,"end_position"]),pos,lwd=lwd,col="red")
segments(max(ann[,"end_position"]),-1*width,max(chrom.bands[chrom.bands[,"chr"]==the.chr,"bases.bot"]),-1.5*width,lwd=lwd,col="red")
#layout()
#layout(matrix(1:3),heights=c(5,1,1))
#layout.show(1)    # Specify layout
#layout.show(2)    # Try specifying 
#layout.show(3)
par(mar=c(0,5.5,0,5.1),mgp=c(3,1,0))

             
the.plot<-plot(imput[,"position"],imput[,"Pval"],pch=pch.array,ylab=expression(bold(-log[10](P[val]))),xlab="",col=colorss[round(imput[,"R2"]*10,0)+1],main="",xlim=c(low.cut,high.cut),ylim=c(0,ymax.range),axes=FALSE,cex=2.0,cex.lab=2.5,font=2,font.lab=2,lwd=2) #  ,bg=colorss[round(imput[,"R2"]*10,0)+1] , xlim=c(low.cut,high.cut)

highlight<-imput[,"SNP"] %in% label.snps
if(sum(highlight)>1){
text(imput[highlight,"position"],imput[highlight,"Pval"],labels=as.character(imput[highlight,"SNP"]),cex=2.0,font=2,pos=4 )
}

label.snps
## the.plot<-plot(imput[,"position"],imput[,"Pval"],pch=pch.array,ylab=expression(bold(-log[10](P[val]))),xlab="",col=colorss[floor(imput[,"R2"]/col.step)+1],main="",xlim=c(low.cut,high.cut),ylim=c(0,ymax.range),axes=FALSE,cex=2.0,cex.lab=2.5,font=2,font.lab=2,lwd=2)

label.snps
axis(2,lty=1,lwd=2,cex.axis=2.0,font=2) 



box()

################ recombination

par(mar=c(0,5.5,0,5.1),mgp=c(3,1,0),new=TRUE)  ## par(fig=c(2/3,1,2/3,1), new=T) ## also works

plot(recomb[,"position"],recomb[,"rate"],type="l",lty=2,xlab="",col="darkred",main="",axes=FALSE,cex=1.5,cex.lab=2,xlim=c(low.cut,high.cut),ylim=c(recomb.low,recomb.high),lwd=2,font.lab=2,ylab="") #  ,bg=colorss[round(imput[,"R2"]*10,0)+1]

axis(4,lty=1,lwd=2,cex.axis=2.0,col="darkred",col.axis="darkred",font=2)
mtext("Recombination Rate (cM/Mb)",side=4,padj=2.0,col="darkred",cex=1.75,las=0,font=2)          


############### ALTERNATIVE PLACEMNET METHOD
#fig=c(2/3,1,2/3,1) #NDC coordinates'c(x1, x2, y1, y2)
## shift<-0.55 # shift to right
## shift<-0    # shidt tp left
## par(fig=c(0.60-shift,1-shift,0.75,0.95),mar=c(4,0.5,3,5.5),mgp=c(3,1,0),new=TRUE) #c(bottom, left, top, right)

## par(fig=c(0.60,1,0.75,0.95),mar=c(4,0.5,3,5.5),mgp=c(3,1,0),new=TRUE)
## temp2<-matrix(0:10)  #temp2<-matrix(1:10,4,5) paste(expression(R^2),"for","rs999999",sep=" ")
## #plot(c(0,1),c(0,1),add=TRUE)
## #box()
## image(temp2, col = colorss,xlab=bquote(R^2~"with"~.(first)),main=expression("\u25CF"~~textstyle(Geneotyped~SNPs)~~~~bold("\u25C7")~~textstyle(Imputed~SNPs)),axes=FALSE,cex=1,cex.lab=2.5,cex.main=2,lwd=3,add=FALSE)
#### proble is that image not plooting in desired location

## ### NOTE USE APPLICATION -> CHARACTER MAP TO GET UNICODE ID
## axis(1,at=seq(0,1,0.1),labels=as.character(seq(0,1,0.1)),cex.axis=1.5,font=2 )
## box()


################### STOP AND CAHNGE SNP NAME
if( left.side==TRUE ){
par(mar=c(4,6.0,2,0.5),mgp=c(1.5,0,0)) #c(bottom, left, top, right) ## left
} else {
par(mar=c(4,0.5,2,5.5),mgp=c(1.5,0,0)) #c(bottom, left, top, right) ## rightpar(mar=c(4,0.5,3,5.5),mgp=c(3.1,1,0))
}




image(temp2, col = colorss[1:11],xlab=bquote(bold(R^2~"with"~.(the.snp))),main=expression("\u25CF"~~bold(Genotyped~SNPs)~~~~bold("\u25C7")~~bold(Imputed~SNPs)),axes=FALSE,cex=1,cex.lab=2.0,cex.main=2.0,lwd=3)
text(seq(0,1,0.1),0.1,labels=as.character(seq(0,1,0.1)),cex=2.0,font=2 )



#### unicode only seems to work with 2.9.2
## image(temp2, col = colorss,xlab=bquote(bold(R^2~"with"~.(the.snp))),main=expression('{//ZapfDingbats \165}'~~bold(Genotyped~SNPs)~~~~bold("\u2666")~~bold(Imputed~SNPs)),axes=FALSE,cex=1,cex.lab=2.0,cex.main=2.0,lwd=3)


## for hand done colors
## image(temp2, col = colorss[1:(num.col-1)],xlab=bquote(bold(R^2~"with"~.(the.snp))),main=expression("\u25CF"~~bold(Genotyped~SNPs)~~~~bold("\u25C7")~~bold(Imputed~SNPs)),axes=FALSE,cex=1,cex.lab=2.0,cex.main=2.0,lwd=3)
## text(seq(0,0.9,0.2),0.1,labels=as.character(seq(0,0.9,0.2)),cex=2.0,font=2 )

## offset<-0.125
## width<-1/(num.col-2)
## at<-c(seq(-offset,1+offset,width))
## labels<-seq(0,1,col.step)
## at<-at[1:length(labels)]
## at
## ##axis(1,at=at,labels=labels,cex.lab=1.7,cex.axis=3.0,font.axis=2)
## text(at,0.1,labels=labels,cex=3.0,font=2,cex.lab=1.7,cex.axis=3.0)

box()
## }else{
## image(temp2, col = "white",xlab="",main=expression("\u25CF"~~bold(Genotyped~SNPs)~~~~bold("\u25C7")~~bold(Imputed~SNPs)),axes=FALSE,cex=1,cex.lab=2.0,cex.main=2.0,lwd=3)
## }


##main=expression("\u25CF"~~bold(Genotyped~SNPs)~~~~bold("\u25C7")~~bold(Imputed~SNPs))
## image(matrix(0:10), col = colorss,xlab=bquote(bold(.(target))),main=expression(bold(Scale)),axes=FALSE,cex=1,cex.lab=1.75,cex.main=1.75,lwd=3)
##
## text(seq(0,1,0.1),0.1,labels=as.character(seq(0,1,0.1)),cex=2.0,font=2 )

### NOTE USE APPLICATION -> CHARACTER MAP TO GET UNICODE ID
#axis(1,at=seq(0,1,0.1),labels=as.character(seq(0,1,0.1)),cex.axis=1.0,font=2,padj=-1.0 )


############################# exons #####################
  if(with.conservation){
par(mar=c(0, 5.5,0,5.1),mgp=c(3,1,0)) #c(bottom, left, top, right)
}else{
par(mar=c(5.1,5.5,0,5.1),mgp=c(3,1,0))
}


 if(with.conservation){ 
bar.height<-0.05
bar.center<-0.08
}else{
bar.height<-0.05
bar.center<-0.08
}
  

exons.ori<-exons
ann.ori<-ann
exons<-exons.ori


#################### warning


 
utr.5.end.plus<-!is.na(exons[,"5_utr_end"]) & exons[,"strand"]==1
utr.5.end.minus<-!is.na(exons[,"5_utr_end"]) & exons[,"strand"]== -1
exons[utr.5.end.plus,"exon_chrom_start"]<-exons[utr.5.end.plus,"5_utr_end"]
exons[utr.5.end.minus,"exon_chrom_end"]<-exons[utr.5.end.minus,"5_utr_start"]

utr.3.start.plus<-!is.na(exons[,"3_utr_start"])  & exons[,"strand"]==1
utr.3.start.minus<-!is.na(exons[,"3_utr_start"])  & exons[,"strand"]== -1
exons[utr.3.start.plus,"exon_chrom_end"]<-exons[utr.3.start.plus,"3_utr_start"]
exons[utr.3.start.minus,"exon_chrom_start"]<-exons[utr.3.start.minus,"3_utr_end"]

no.nas<-!is.na(exons[,"exon_chrom_end"]) & !is.na(exons[,"exon_chrom_start"] )
exon.centers<-abs(exons[no.nas,"exon_chrom_end"]+exons[no.nas,"exon_chrom_start"])/2
y.vals<-bar.center*exons[no.nas,"strand"]

#plot(x=c(low.cut,high.cut) , y= c(-1*(bar.height+bar.center),(bar.height+bar.center)) ,xlim=c(low.cut,high.cut),type="n",axes=TRUE,ylab="",xlab="") ## protype
                                       # original is below
#plot(x=c(min(imput[,"position"]),max(imput[,"position"])) , y= c(-1*(bar.height+bar.center),(bar.height+bar.center)) ,xlim=c(low.cut,high.cut),type="n",axes=FALSE,ylab="",xlab="") ##run

if(with.conservation){
plot(x=c(min(imput[,"position"]),max(imput[,"position"])) , y= c(-1*(bar.height+bar.center),(bar.height+bar.center)) ,xlim=c(low.cut,high.cut),type="n",cex=2.0,cex.lab=2.4,lwd=3,cex.axis=2.0,axes=F,ylab="",xlab="",font=2,font.lab=2)
}else{
plot(x=c(min(imput[,"position"]),max(imput[,"position"])) , y= c(-1*(bar.height+bar.center),(bar.height+bar.center)) ,xlim=c(low.cut,high.cut),type="n",cex=2.0,cex.lab=2.4,lwd=3,cex.axis=2.0,axes=FALSE,ylab="",xlab="Position (bp)",font=2,font.lab=2)
axis(1,lty=1,lwd=2,cex.axis=2.0,font=2)
}

###################################################################################
## rectangles[1:5,]
## y.vals[1:5]
## exon.centers[1:5]
## 50311864 
## symbols(50300000, -0.08, rectangles=cbind(10,0.05 ), inches=FALSE, fg ="red",bg ="red",add=TRUE)

 symbols(exon.centers, y.vals, rectangles=cbind(abs(exons[no.nas,"exon_chrom_end"]-exons[no.nas,"exon_chrom_start"]),bar.height ), inches=FALSE, fg =col.array[exons[no.nas,"ensembl_gene_id"]],bg =col.array[exons[no.nas,"ensembl_gene_id"]],add=TRUE,)  ### exons
############# 5 utr
no.nas<-!is.na(exons[,"5_utr_end"]) & !is.na(exons[,"5_utr_start"])
exon.centers<-abs(exons[no.nas,"5_utr_end"]+exons[no.nas,"5_utr_start"])/2
y.vals<-bar.center*exons[no.nas,"strand"]
symbols(exon.centers, y.vals, rectangles=cbind(abs(exons[no.nas,"5_utr_end"]-exons[no.nas,"5_utr_start"]),bar.height/4 ), inches=FALSE,  ,bg =col.array[exons[no.nas,"ensembl_gene_id"]],fg =col.array[exons[no.nas,"ensembl_gene_id"]],add=TRUE)  ### utr  fg =col.array[exons[no.nas,"ensembl_gene_id"]]
##########3 utr
no.nas<-!is.na(exons[,"3_utr_end"]) & !is.na(exons[,"3_utr_start"])
exon.centers<-abs(exons[no.nas,"3_utr_end"]+exons[no.nas,"3_utr_start"])/2
y.vals<-bar.center*exons[no.nas,"strand"]
symbols(exon.centers, y.vals, rectangles=cbind(abs(exons[no.nas,"3_utr_end"]-exons[no.nas,"3_utr_start"]),bar.height/4 ), inches=FALSE,  ,bg =col.array[exons[no.nas,"ensembl_gene_id"]],fg =col.array[exons[no.nas,"ensembl_gene_id"]],add=TRUE)  ### utr  fg =col.array[exons[no.nas,"ensembl_gene_id"]]

gene.centers<-abs(ann[,"start_position"]+ann[,"end_position"])/2
y.vals<-bar.center*ann[,"strand"]
the.labels<-ann[,"hgnc_symbol"]
############################### choose default label
## if(length(the.labels)==1 & is.na(the.labels[1])){the.labels<-ann[,"ensembl_gene_id"]}
## the.labels[the.labels==""]<-ann[the.labels=="","ensembl_gene_id"]

if(length(the.labels)==1 & is.na(the.labels[1])){the.labels<-ann[,"external_gene_name"]}
the.labels[the.labels==""]<-ann[the.labels=="","external_gene_name"]

keep<-rep(TRUE,times=length(the.labels))
#offsets<-rep(2.2,times=length(the.labels))
poss<-rep(3,times=length(the.labels))
to.shift.low<-FALSE
to.shift.amount.low<-0
to.shift.amount.high<-0
to.shift.high<-FALSE
set.below<-{}
to.ignore<-{}

  ############# change for a specific case below


## if(the.gene==12){
## to.shift.amount.high<-0
## to.shift.high<-FALSE
## set.below<-c(3)
## to.ignore<-{}
## }

## if(the.gene==13){
## to.shift.amount.low<-310000
## to.shift.low<-TRUE
## to.shift.high<-FALSE
## set.below<-c(2)
## to.ignore<-{}
## }

## if(the.gene==15){
## to.shift.amount.high<-50000
## to.shift.low<-FALSE
## to.shift.high<-TRUE
## set.below<-{}
## to.ignore<-{}
## }




the.labels
which.max(gene.centers)
gene.centers[which.max(gene.centers)]
gene.centers[which.max(gene.centers)]-100000

which.min(gene.centers)
gene.centers[which.min(gene.centers)]
gene.centers[which.min(gene.centers)]+310000

keep[to.ignore]<-FALSE  
poss[set.below]<- 1
gene.centers<-gene.centers[keep]
if(to.shift.low==TRUE){gene.centers[which.min(gene.centers)]<-gene.centers[which.min(gene.centers)]+to.shift.amount.low}
if(to.shift.high==TRUE){gene.centers[which.max(gene.centers)]<-gene.centers[which.max(gene.centers)]-to.shift.amount.high}


########## alternate gene name above and beo the line
lower<-ann[keep,"strand"]==1 ## boolean array every sesond one what to shift down by 0.08
sum(lower)
shifter<-rep(c(0,bar.center),times=length(lower)) # longer than needed
shifter<-shifter[1:sum(lower)]
y.vals[keep][lower]<-y.vals[keep][lower]-shifter


lower<-ann[keep,"strand"]==-1 ## boolean array every sesond one what to shift down by 0.08
sum(lower)
shifter<-rep(c(0,bar.center),times=length(lower)) # longer than needed
shifter<-shifter[1:sum(lower)]
y.vals[keep][lower]<-y.vals[keep][lower]-shifter
###########################################################

text(gene.centers,y.vals[keep],labels=the.labels[keep],pos=poss[keep],offset=2.0,col=col.array[ann[keep,"ensembl_gene_id"]],lwd=2,font=2,font.lab=2,cex=1.2)





#cbind(gene.centers,y.vals[keep],labels=the.labels[keep])

num.arrows<-25
x <- seq(from=low.cut,to=high.cut,by =(high.cut-low.cut)/num.arrows)
y<- rep(bar.center,times=length(x))
s <- seq(length(x)-1)# one shorter than data
arrows(x[s], y[s], x[s+1], y[s+1], col="grey",angle=20,code=2,length=0.15,lwd=1.5)
x <- seq(from=high.cut,to=low.cut,by =(low.cut-high.cut)/num.arrows)
y<- rep(-1*bar.center,times=length(x))
arrows(x[s], y[s], x[s+1], y[s+1], col="grey",angle=20,code=2,length=0.15,lwd=1.5)
box()


################ Conservation ###########
#GENOME TABLE BROWER use region below
#paste(paste("chr",imput[1,"chrom"],sep=""),paste(min(imput[,"position"]),max(imput[,"position"]),sep="-"),sep=":")
# Comparative genomics ; 17-way Most Cons ; 
#Output as data points, set filter of points > 0 and 10,000,000 lines

if(with.conservation){                           
file<-paste("chr",the.gene,".conservation",".txt",sep="")
conserv<-read.delim(file,header=F,skip=10,sep="",fill=TRUE,stringsAsFactors=FALSE)
#conserv<-conserv[conserv[,2]!=0,]

#par(mar=c(0,5.5,0,2.1),mgp=c(3,1,0)) #c(bottom, left, top, right)
par(mar=c(5.1,5.5,0,5.1),mgp=c(3,1,0)) #c(bottom, left, top, right)


plot(x=c(low.cut,high.cut) , y= c(min(conserv[,6]),max(conserv[,6])) ,xlim=c(low.cut,high.cut),type="n",cex=2.0,cex.lab=2.4,lwd=3,cex.axis=2.0,axes=TRUE,ylab="Conservation",xlab="Position (bp)",font=2,font.lab=2)

symbols((conserv[,3]+conserv[,4])/2,conserv[,6]/2, rectangles=cbind(abs((conserv[,3]-conserv[,4])),conserv[,6] ), inches=FALSE,  ,bg ="darkblue",fg ="darkblue",add=TRUE)  ### utr
}


######################################################### end plot







############################################END
################################## IDEOGRAN ##################################
##  paint.chromosomeBP<-function (chrom, pos = 0, units = "cM", width = 0.4, bands = "major") 
## {
##     semicircle <- function(base.x, base.y, base.length, height = base.length, 
##         side = 1, orientation = NULL, col = NULL) {
##         radius <- base.length/2
##         x <- radius * seq(-1, 1, length = 40)
##         y <- height/radius * sqrt(radius^2 - x^2)
##         if (is.null(orientation)) {
##             co <- as.integer(cos(pi * (3 - side)/2))
##             so <- as.integer(sin(pi * (3 - side)/2))
##         }
##         else {
##             co <- cos(orientation)
##             so <- sin(orientation)
##         }
##         tx <- co * x - so * y
##         ty <- so * x + co * y
##         if (is.null(orientation)) {
##             if (side == 1 || side == 3) {
##                 base.x <- base.x + radius
##             }
##             else if (side == 2 || side == 4) {
##                 base.y <- base.y + radius
##             }
##         }
##         x <- base.x + tx
##         y <- base.y + ty
##         polygon(x, y, col = col)
##     }
##     data(chrom.bands)
##     chromdata <- subset(chrom.bands, chrom.bands$chr == chrom)
##     lc <- nchar(chromdata$band)
##     sel <- !(substr(chromdata$band, lc, lc) %in% letters)
##     if (bands != "major") 
##         sel <- !sel
##     chromdata <- chromdata[sel, ]
##     rm(lc, sel)
##     bandcol <- gray(c(0.4, 0.6, 0.8, 0.8, 0.85))[match(chromdata$stain, 
##         c("acen", "gneg", "gpos", "gvar", "stalk"))]
##     n <- nrow(chromdata)
##     centromere <- which(chromdata$arm[-n] != chromdata$arm[-1])
##     idx <- c(2:(centromere - 1), (centromere + 2):(n - 1))
##     rect(chromdata$bases.top[idx], pos, chromdata$bases.bot[idx], pos - 
##         width, col = bandcol[idx])
##     semicircle(chromdata$bases.bot[1], pos - width, width, chromdata$bases.bot[1] - 
##         chromdata$bases.top[1], 2, col = bandcol[1])
##     semicircle(chromdata$bases.top[n], pos - width, width, chromdata$bases.bot[n] - 
##         chromdata$bases.top[n], 4, col = bandcol[n])
##     semicircle(chromdata$bases.top[centromere], pos - width, width, 
##         chromdata$bases.bot[centromere] - chromdata$bases.top[centromere], 
##         4, col = bandcol[centromere])
##     semicircle(chromdata$bases.bot[centromere + 1], pos - width, 
##         width, chromdata$bases.bot[centromere + 1] - chromdata$bases.top[centromere + 
##             1], 2, col = bandcol[centromere + 1])
##     points(chromdata$bases.top[centromere], pos - 0.5 * width, col = "black", 
##         cex = 3, pch = 16)
##     points(chromdata$bases.top[centromere], pos - 0.5 * width, col = "white", 
##         cex = 3, pch = 20)
## }
## ########################################################################













## ############### BIOMART ###########
## library(GenomeGraphs)

## library(biomaRt)
## mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
## attributes = listAttributes(mart)
## filters=listFilters(mart)


##  attributes[grep("id",attributes[,1]),]
##  hgnc_curated_gene_name,  hgnc_symbol hgnc_id  entrezgenegr

## attributes[c(13,33,35,38,40,46,47,55,56,57,65,71,72,78), ]


           
## ##  "ensembl_gene_id"
## ## "ensembl_transcript_id"
## ## "ensembl_exon_id"
## ## "exon_chrom_start"
## ## "exon_chrom_end"
## ## "rank"
## ## "strand"
## ## "gene_biotype"
## ## "hgnc_symbol"
## ## "hgnc_id"-
## ### PLUS STRAND 5' -> 3'
##                            name                   description
## 13              external_gene_id          Associated Gene Name
## 33 clone_based_ensembl_gene_name Clone based Ensembl gene name
## 35    clone_based_vega_gene_name    Clone based VEGA gene name
## 38                          embl             EMBL (Genbank) ID
## 40                          ottt  VEGA transcript ID(s) (OTTT)
## 46      hgnc_automatic_gene_name      HGNC automatic gene name/
## 47        hgnc_curated_gene_name        HGNC curated gene name/
## 55            mim_gene_accession            MIM Gene Accession
## 56          mim_gene_description          MIM Gene Description
## 57                       mirbase                       miRBase
## 65                          ucsc                       UCSC ID
## 71                   wikigene_id                   WikiGene ID
## 72                 wikigene_name                 WikiGene name
## 78                   dbass5_name              DBASS5 Gene Name


##  fil.vals<-list(as.character(imput[1,"chrom"]), imput[1,"position"], imput[dim(imput)[1],"position"])
## ann<-getBM(attributes = c( "ensembl_gene_id","external_gene_id","chromosome_name","start_position","end_position","strand","hgnc_symbol","gene_biotype","clone_based_ensembl_gene_name"),filters = a.filter, values=fil.vals, mart = mart)
## ann

## a.filter<-c( "chromosome_name", "start" , "end", "strand")
## fil.vals<-list(as.character(imput[1,"chrom"]), imput[1,"position"], imput[dim(imput)[1],"position"],"1")
## plus<-getBM(attributes = c("ensembl_gene_id","5_utr_start","5_utr_end","ensembl_exon_id","exon_chrom_start","exon_chrom_end","3_utr_start","3_utr_end","rank","strand","gene_biotype"), filters = a.filter, values=fil.vals, mart = mart)
##  unique(plus[,1])
## ## plus<-getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","5_utr_start","3_utr_start","ensembl_exon_id","exon_chrom_start","exon_chrom_end","rank","strand","gene_biotype"), filters = a.filter, values=fil.vals, mart = mart)

## fil.vals<-list(as.character(imput[1,"chrom"]), imput[1,"position"], imput[dim(imput)[1],"position"],"1")
## ann.plus<-getBM(attributes = c( "ensembl_gene_id","chromosome_name","start_position","end_position","hgnc_symbol","gene_biotype"), filters = a.filter, values=fil.vals, mart = mart)

##  require(stats); require(grDevices)
##      x <- 1:10
##      y <- sort(10*runif(10))
##      z <- runif(10)
##      z3 <- cbind(z, 2*runif(10), runif(10))
##      symbols(x, y, thermometers=cbind(.5, 1, z), inches=.5, fg = 1:10)
##      symbols(x, y, thermometers = z3, inches=FALSE)
##      text(x,y, apply(format(round(z3, digits=2)), 1, paste, collapse = ","),
##           adj = c(-.2,0), cex = .75, col = "purple", xpd=NA)

##  arrows(x0, y0, x1, y1, length = 0.25, angle = 30, code = 2,
##             col = par("fg"), lty = par("lty"), lwd = par("lwd"),
##             ...)

##  x <- stats::runif(12); y <- stats::rnorm(12)
##      i <- order(x,y); x <- x[i]; y <- y[i]
##      plot(x,y, main="arrows(.) and segments(.)")
##      ## draw arrows from point to point :
##      s <- seq(length(x)-1)# one shorter than data
##      arrows(x[s], y[s], x[s+1], y[s+1], col= 1:3)
##      s <- s[-length(s)]
##      segments(x[s], y[s], x[s+2], y[s+2], col= 'pink')


## ##################### GENOME GRAPHS
## plusStrand <- makeGeneRegion(chromosome = imput[1,"chrom"], start = imput[1,"position"], end =  imput[dim(imput)[1],"position"], strand = "+", biomart = mart,dp=DisplayPars(plotID=TRUE,idRotation=0,idColor="red",cex=0.5))

## minStrand <- makeGeneRegion(chromosome = imput[1,"chrom"], start = imput[1,"position"], end =  imput[dim(imput)[1],"position"], strand = "-", biomart = mart)
## ideogram <- makeIdeogram(chromosome = imput[1,"chrom"])
## genomeAxis <- makeGenomeAxis(add53 = TRUE,add35=TRUE)
## gdPlot(list(ideogram,plusStrand,genomeAxis, minStrand))


## ###2
## plusStrand <- makeGeneRegion(chromosome = imput[1,"chrom"], start =  96130000, end =  96140000, strand = "+", biomart = mart, dp=DisplayPars(color="red",protein_coding="blue"))
##   to <- makeTextOverlay("here is some text",xpos = 96133000, ypos = 0.5)
## genomeAxis <- makeGenomeAxis(add53 = TRUE,add35=TRUE)
## gdPlot(list(plusStrand,genomeAxis),overlay=c(to),add=TRUE)
## text(1,1,"bob")
## locator()
## ##############################

## a <- new("GenomeAxis")
## setPar(a, "size", 100) ### only gets dp parameters
## gdPlot(a, minBase = 10, maxBase = 10000)


## legend(x=96250000,10,legend= paste(expression(R^2),"=" ,as.character(seq(0,1,0.1)),sep=" "),col=colorss,pch=19)

## > par()$mar
## par(mar=c(5.1,4.1,4.1,2.1))
##  par(mar=c(0,3,1,1))
##      plot(x, y, xlim=xrange, ylim=yrange, xlab="", ylab="")
##      par(mar=c(0,3,1,1))
##      barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
##      par(mar=c(3,0,1,1))


## savePlot("Arts1.jpeg",type="jpeg")
## print(the.plot,position = c(0, .3, 1, .9),sep=" "), more = TRUE)
##      print(update(the.plot, aspect = "xy", main = "", xlab = "Year"),
##            position = c(0, 0, 1, .3))
 
##    print("TETS",position = c(0, 0, 1, .3))
## ]






## ########## get goodies for genome graphs
## library(GenomeGraphs)
## mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
## minStrand <- makeGeneRegion(chromosome = 22,start = 30450000, end = 30550000,strand = "-", biomart = mart)
## ideogram <- makeIdeogram(chromosome = 22)
## genomeAxis <- makeGenomeAxis(add53 = TRUE,add35 = TRUE)
## ####################

## ########## subplot
## ot(0:10, 0:10, type='n')
## x <- rnorm(100)
## plot(x)
## par(fig=c(2/3,1,2/3,1), new=T)
## hist(x, main="")
## ###################################

## nf <- layout(matrix(c(1,2), byrow=TRUE))
## layout.show(nf)
## x<-1:100
## y<-x^2
## plot(x,y)
## gdPlot(list(ideogram, plusStrand, genomeAxis, minStrand), minBase = 30450000, maxBase = 30550000)

## #layout()
## #layout(matrix(1:3),heights=c(5,1,1))
## #layout.show(1)    # Specify layout
## #layout.show(2)    # Try specifying 
## #layout.show(3)
## par(mar=c(1,5.5,4.1,2.1),mgp=c(3,1,0)) #c(bottom, left, top, right)
## the.plot<-plot(imput[,"position"],imput[,"Pval"],pch=19,ylab=expression(-log[10](P[val])),xlab="",col=colorss[round(imput[,"R2"]*10,0)+1],main="ARTS1",axes=FALSE,cex=1.5,cex.lab=2)
## axis(2,lty=1,lwd=2,cex.axis=1.75)
## box()
## par(mar=c(4,0,4,2.5),mgp=c(3,1,0))
## temp2<-matrix(0:10)  #temp2<-matrix(1:10,4,5) paste(expression(R^2),"for","rs999999",sep=" ")
##  image(temp2, col = colorss,xlab=expression(R^2~"for"~rs999999),main=expression("\u25CF"~~textstyle(Geneotyped~SNPs)~~~~"\u25C6"~~textstyle(Imputed~SNPs)),axes=FALSE,cex=1,cex.lab=2,cex.main=2,lwd=3)
## ### NOTE USE APPLICATION -> CHARACTER MAP TO GET UNICODE ID
## axis(1,at=seq(0,1,0.1),labels=as.character(seq(0,1,0.1)),cex.axis=1.5 )
## box()




## axis(2, at=p, labels=format(p, decimal.mark="\u00B7"))




## ############### BIOMART ###########
## library(GenomeGraphs)
## mart <- useMart("ensembl")
## mart <- useMart("ensembl_mart_50")
## listMarts(archive=TRUE)
## attributes =
##   listDatasets(mart)
## listDatasets(ensembl)
## , dataset="hsapiens_gene_ensembl")
## attributes = listAttributes(mart)
## filters=listFilters(mart)

## listMarts(archive=TRUE)
##                        biomart                     version
## 1              ensembl_mart_51                  Ensembl 51
## 2                  snp_mart_51                      SNP 51
## 3                 vega_mart_51                     Vega 32
## 4              ensembl_mart_50                  Ensembl 50
## 5                  snp_mart_50                      SNP 50
## 6                 vega_mart_50                     Vega 32

## Error in value[[3L]](cond) : 
##   Request to BioMart web service failed. Verify if you are still connected to the internet.  Alternatively the BioMart web service is temporarily down.
## In addition: Warning message:
## In file(file, "r") : unable to resolve 'july2008.archive.ensembl.org'

## mart<-  useMart("ensembl_mart_50",dataset="hsapiens_gene_ensembl",archive=TRUE)
## mart<-  useMart("ensembl_mart_49",dataset="hsapiens_gene_ensembl",archive=TRUE)

##  attributes[grep("utr",attributes[,1]),]
##  hgnc_curated_gene_name,  hgnc_symbol hgnc_id  entrezgenegr

## attributes[1:5, ]
## ##  "ensembl_gene_id"
## ## "ensembl_transcript_id"
## ## "ensembl_exon_id"
## ## "exon_chrom_start"
## ## "exon_chrom_end"
## ## "rank"
## ## "strand"
## ## "gene_biotype"
## ## "hgnc_symbol"
## ## "hgnc_id"-
## ### PLUS STRAND 5' -> 3'

## a.filter<-c( "chromosome_name", "start" , "end", "strand")
## fil.vals<-list(as.character(imput[1,"chrom"]), imput[1,"position"], imput[dim(imput)[1],"position"],"1")
## plus<-getBM(attributes = c("ensembl_gene_id","5_utr_start","5_utr_end","ensembl_exon_id","exon_chrom_start","exon_chrom_end","3_utr_start","3_utr_end","rank","strand","gene_biotype"), filters = a.filter, values=fil.vals, mart = mart)
##  unique(plus[,1])
## ## plus<-getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","5_utr_start","3_utr_start","ensembl_exon_id","exon_chrom_start","exon_chrom_end","rank","strand","gene_biotype"), filters = a.filter, values=fil.vals, mart = mart)

## fil.vals<-list(as.character(imput[1,"chrom"]), imput[1,"position"], imput[dim(imput)[1],"position"],"1")
## ann.plus<-getBM(attributes = c( "ensembl_gene_id","chromosome_name","start_position","end_position","hgnc_symbol","gene_biotype"), filters = a.filter, values=fil.vals, mart = mart)

##  require(stats); require(grDevices)
##      x <- 1:10
##      y <- sort(10*runif(10))
##      z <- runif(10)
##      z3 <- cbind(z, 2*runif(10), runif(10))
##      symbols(x, y, thermometers=cbind(.5, 1, z), inches=.5, fg = 1:10)
##      symbols(x, y, thermometers = z3, inches=FALSE)
##      text(x,y, apply(format(round(z3, digits=2)), 1, paste, collapse = ","),
##           adj = c(-.2,0), cex = .75, col = "purple", xpd=NA)

## ##################### GENOME GRAPHS
## plusStrand <- makeGeneRegion(chromosome = imput[1,"chrom"], start = imput[1,"position"], end =  imput[dim(imput)[1],"position"], strand = "+", biomart = mart,dp=DisplayPars(plotID=TRUE,idRotation=0,idColor="red",cex=0.5))

## minStrand <- makeGeneRegion(chromosome = imput[1,"chrom"], start = imput[1,"position"], end =  imput[dim(imput)[1],"position"], strand = "-", biomart = mart)
## ideogram <- makeIdeogram(chromosome = imput[1,"chrom"])
## genomeAxis <- makeGenomeAxis(add53 = TRUE,add35=TRUE)
## gdPlot(list(ideogram,plusStrand,genomeAxis, minStrand))


## ###2
## plusStrand <- makeGeneRegion(chromosome = imput[1,"chrom"], start =  96130000, end =  96140000, strand = "+", biomart = mart, dp=DisplayPars(color="red",protein_coding="blue"))
##   to <- makeTextOverlay("here is some text",xpos = 96133000, ypos = 0.5)
## genomeAxis <- makeGenomeAxis(add53 = TRUE,add35=TRUE)
## gdPlot(list(plusStrand,genomeAxis),overlay=c(to),add=TRUE)
## text(1,1,"bob")
## locator()
## ##############################

## a <- new("GenomeAxis")
## setPar(a, "size", 100) ### only gets dp parameters
## gdPlot(a, minBase = 10, maxBase = 10000)


## legend(x=96250000,10,legend= paste(expression(R^2),"=" ,as.character(seq(0,1,0.1)),sep=" "),col=colorss,pch=19)

## > par()$mar
## par(mar=c(5.1,4.1,4.1,2.1))
##  par(mar=c(0,3,1,1))
##      plot(x, y, xlim=xrange, ylim=yrange, xlab="", ylab="")
##      par(mar=c(0,3,1,1))
##      barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
##      par(mar=c(3,0,1,1))


## savePlot("Arts1.jpeg",type="jpeg")
## print(the.plot,position = c(0, .3, 1, .9),sep=" "), more = TRUE)
##      print(update(the.plot, aspect = "xy", main = "", xlab = "Year"),
##            position = c(0, 0, 1, .3))
 
##    print("TETS",position = c(0, 0, 1, .3))
## ]







## library(GenomeGraphs)
## mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
## minStrand <- makeGeneRegion(chromosome = 17,start = 30450000, end = 30550000,strand = "-", biomart = mart)
## ideogram <- makeIdeogram(chromosome = 17)
## genomeAxis <- makeGenomeAxis(add53 = TRUE,add35 = TRUE)


## gdPlot(list(ideogram, plusStrand, genomeAxis, minStrand), minBase = 30450000, maxBase = 30550000)


## ####################################### useful examples for plottting
## X \264 ; Delta; lozenge; bullet \267; diamond \250
## ########### raw standard plot:
## plot.new();
## plot.window(c(0,10), c(0,10))
## k=200
## for(i in 1:10){
## for(j in 1:10){
##  text(i, j, expression(symbol(paste("\\",k,sep=""))))
## }}

## plot.new(); plot.window(c(0,4), c(15,1))
## text(1, 1, "universal", adj=0); text(2.5, 1,  "\\042"); text(3, 1, expression(symbol("\250")))

## the.plot<-plot(imput[,"position"],imput[,"Pval"],pch=19,col=colorss[round(imput[,"R2"]*10,0)+1],main="ARTS1",ylab=expression(-log[10](Pval)),xlab="Chromosome Position",)
## legend(x=96250000,10,legend= paste(expression(R^2),"=" ,as.character(seq(0,1,0.1)),sep=" "),col=colorss,pch=19)



## legend("topright",col=c("red","green","blue"), pch=c(15,16,17),
## + lty=c(1,2,3), legend=c("series 1", "series 2", "series 3"),
## + inset=0.05, bg='white')


##  demo(plotmath)


## dev.copy2pdf(file="blackbody.pdf")

## #

## matrix(1:2)
## matrix(1:4)              # 4x1
## matrix(1:4,2,2)          # 2x2
## matrix(1:6,3,2)          # 3x2 ordered by columns
## matrix(1:6,3,2,byrow=T)  # 3x2 ordered by rows
		

## # To view the graphical layout, the following will show the borders of the sub-panels and the number identying each one:

## layout(matrix(1:3))
## layout.show(1)    # Specify layout for 4 panels, for the defined layout
## layout.show(2)    # Try specifying just 2 instead
## layout.show(3)		

## # Now fill the layout with 4 plots:

## x <- 1:10
## plot(x,x)
## plot(x,x^2)
## plot(x,sqrt(x))
## plot(x,log10(x))
## curve(log10,add=T)  # Adds to last panel plotted

## #####
## ?symbold plot symbols

## The "heights" and "widths" arguments to "layout" are vectors of relative heights and widths of the matrix rows and columns, respectively.

## # Specifying panels of different sizes:

##  plot <- xyplot(sunspot.year ~ 1700:1988, xlab = "", type = "l",
##                     scales = list(x = list(alternating = 2)),
##                     main = "Yearly Sunspots")
##      print(plot, position = c(0, .3, 1, .9), more = TRUE)
##      print(update(plot, aspect = "xy", main = "", xlab = "Year"),
##            position = c(0, 0, 1, .3))
 


## library(sp)
## data(meuse.grid)
## coordinates(meuse.grid) <- c("x", "y")
## gridded(meuse.grid) <- TRUE

## With lattice graphics:

## l1 <- list("SpatialPolygonsRescale", layout.scale.bar(),
##   offset = c(180500,329800), scale = 500, fill=c("transparent","black"))
## l2 = list("sp.text", c(180500,329900), "0")
## l3 = list("sp.text", c(181000,329900), "500 m")
## spplot(meuse.grid, "dist", col.regions=grey.colors(20),
##   sp.layout=list(l1, l2, l3))

## With base graphics:

## image(meuse.grid, "dist", col=grey.colors(20))
## SpatialPolygonsRescale(layout.scale.bar(), offset = c(180500, 329800),
##   scale = 500, fill=c("transparent","black"), plot.grid=FALSE)
## text(180500, 329900, "0")
## text(181000, 329900, "500 m")


## layout(matrix(1:4,2,2),heights=c(2,1)); layout.show(4)
## replicate(4,plot(x,x))  # Repeat plot 4 times
## with ( possum , plot ( density ( totlngth [ here ]) , type = " l " ))

##     *  Spline interpolation of data

##       x <- 1:20
##       y <- jitter(x^2,factor=20)    # Add some noise to vector
##       sf <- splinefun(x,y)          # Perform cubic spline interpolation
##       # - note that "sf" is a function:
##       sf(1.5)                       # Evaluate spline function at x=1.5
##       plot(x,y); curve(sf(x),add=T) # Plot data points & spline curve

##     * Scatter plot smoothing

##       #--Using above data frame ("A"):
##       plot(A$z,A$Tx)
##       #--Return (X &) Y values of locally-weighted polynomial regression
##       lowess(A$z,A$Tx)
##       #--Plot smoothed data as line on graph:
##       lines(lowess(A$z,A$Tx),col="red")

##     * Numerical integration of a function (see "Functions" section below)

##       # Create simple function:
##       fun <- function(x,norm,index) norm*x^index
##       # see "?integrate" for details. Note that:
##       #  1) the names of arguments in "fun" must not match those of arguments
##       #     in "integrate()" itself.
##       #  2) fun must be able to return a vector, if supplied with a vector
##       #     (i.e. not just a single value)
##       # 
##       i <- integrate(fun,lower=1,upper=10,norm=0.5,index=2)
##       > i
##       166.5 with absolute error < 1.8e-12
##       > i[0]
##       list()                         # Note that integrate() returns a list
##       > names(i)                     # Show names components in list
##       [1] "value"        "abs.error"    "subdivisions" "message"      "call"
##       > i$value                      # If you want the integral value alone
##       [1] 166.5



















