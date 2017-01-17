
 library(gplots)

load("clustering.RData")
exclude<-c(14,28)  #14_C and 14J excluded from analysis
gene.list2<-read.delim("BRB 465 gene list2.txt",header=T,sep="\t",fill=TRUE)
gene.list2<-as.character(gene.list2[,"Unique.id"])


data<-as.matrix(selDataMatrix[gene.list,-exclude])
data<-as.matrix(selDataMatrix[gene2.list,-exclude])

colnames(data)<-some.labels[-exclude]

# HCLUST method: the agglomeration method to be used. This should be (an
#           unambiguous abbreviation of) one of '"ward"', '"single"',
#           '"complete"', '"average"', '"mcquitty"', '"median"' or
#           '"centroid"'.
### centered correlation average linkage
##cor use: method = c("pearson", "kendall", "spearman")





###################heatmaps for paper
 plot(hclust(dist(t(data))), ColSideColors=coloursC[colour.posns]   )
 
patients<-hclust(dist(t(data)))
patients_dd <- as.dendrogram(patients)
heatm<-heatmap(as.matrix(selDataMatrix[gene.list,-exclude]),col=heat.colors(8),labCol=some.labels[-exclude],labRow="",ColSideColors=coloursC[colour.posns],xlab="Sample",margins=c(6,6))

##tests
plot(patients_dd)
plot(kinase_dd)
patients.reorder<-patients_dd

patients.reorder <- reorder(patients_dd,heatm$colInd)  # note sure about this from
                                                        # plot(patients_dd) & order.dendrogram(patients_dd)
plot(patients_dd)
 plot(patients.reorder)

jpeg("AS_heatmap_dendro.jpeg",quality=100)



##########################################
mydist<-as.dist((1-cor(data,method="kendall")))

plot(hclust(mydist,method="ward"))     # 3 normals in AS group

patients<-hclust(mydist,method="correlation")
plot(hclust(mydist,method="average"),main="Hierarchial Clustering Dendrogram with Filtered data")

heatmap(data, col=rbg, distfun = function(c){as.dist(1 - cor(c,method = "pearson"))},hclustfun= function (d) {hclust(d, method="ward")} )


heatmap(data)

arrays.t <- t(arrays)
arrays.dist <-dist(arrays.t)
plot(hclust(arrays.dist),labels=arrays,main="Hierarchial Clustering Dendrogram with Filtered data")

 patients<-hclust(dist(t(data),method="euclidean"))
 kinases<- hclust(dist(data))


exclude<-c(14,28)
gene.list<-read.delim("ClassComparison 1000 10 80.txt",header=T,sep="",fill=TRUE)
gene.list<-gene.list[,1]
gene.list<-as.character(gene.list)

  library(pamr)
data<-pamr.listgenes(trained, pam.data, threshold=2.789,genenames=TRUE)
data<-read.delim("gene list 70.txt",header=T,fill=TRUE)
 gene.list<-data[,"id"]
 gene.list<-as.character(gene.list)

> colnames(data)
[1] "id"        "name"      "AS.score"  "Nor.score"

##heatmap with gene names
 #heatmap(as.matrix(diff_expressed_korn),col=cm.colors(16),main="Class and Genes (KORN)",xlab="Class")


############## Col colours 
coloursC <-  rainbow(2)

colour.posns<-class.labels[-exclude]
colour.posns[ colour.posns=="AS"]<-1
colour.posns[ colour.posns=="Nor"]<-2
colour.posns<-as.numeric(colour.posns)    

############## Row colours use if have too many genes
genes.pca <- prcomp(selDataMatrix[gene.list,-exclude],scale = TRUE)
genes.pca.var <- round(genes.pca$sdev^2 / sum(genes.pca$sdev^2)*100,2)
plot(c(1:length(genes.pca.var)),genes.pca.var,type="b",xlab="# components",ylab="% variance",main="Scree plot-Genes ",col="blue")






 ###first 2 component have most or varaiblity so use K-means on frist two components
 ##from Heat map 7-8 clusters available

num.clusters<-16
coloursR <-  rainbow(num.clusters)
genes.cl<-kmeans(genes.pca$x[,1:2],centers=num.clusters,iter.max=500) #Do kmeans


################ Spectral clustering:
plot(range(genes.pca$x[,1]),range(genes.pca$x[,2]),xlab="PCA1",ylab="PCA2",main="Spectral Clustering of Genes")

#points(sam_genes.pca$x[,1],sam_genes.pca$x[,2],col=colours[sam_genes.cl$cluster],pch=sam_genes.cl$cluster,cex=1.0,bg=colours[sam_genes.cl$cluster]) #colours and symbols
text(genes.pca$x[,1],genes.pca$x[,2],label=rownames(genes.pca$x),col=coloursR[genes.cl$cluster],cex=0.75) #colours and symbols
text(genes.pca$x[,1],genes.pca$x[,2],label=as.character(genes.cl$cluster),col=coloursR[genes.cl$cluster],cex=0.75)
 median(abs(gene.means[labels(genes.cl$cluster[genes.cl$cluster==1])]))

#  >  median((gene.means[labels(genes.cl$cluster[genes.cl$cluster==8])]))
# [1] 8.689024
# >  median((gene.means[labels(genes.cl$cluster[genes.cl$cluster==7])]))
# [1] 10.63088
# >  median((gene.means[labels(genes.cl$cluster[genes.cl$cluster==6])]))
# [1] 11.34963
# >  median((gene.means[labels(genes.cl$cluster[genes.cl$cluster==5])]))
# [1] 12.97176
# >  median((gene.means[labels(genes.cl$cluster[genes.cl$cluster==4])]))
# [1] 9.199757
# >  median((gene.means[labels(genes.cl$cluster[genes.cl$cluster==3])]))
# [1] 9.938235
# >  median((gene.means[labels(genes.cl$cluster[genes.cl$cluster==2])]))
# [1] 7.387274
# >  median((gene.means[labels(genes.cl$cluster[genes.cl$cluster==1])]))
##############################################################################################


############ dendrogram
expressed.t <- t(selDataMatrix[gene.list,-exclude])    #samples
expressed.t <- (selDataMatrix[gene.list,-exclude])    #genes

expressed.dist <-dist(expressed.t)
plot(hclust(expressed.dist),hang=-1,main="Hierarchial Clustering Dendrogram with Diff. Exp. Genes ")   #labels=dimnames(selDataMatrix[gene.list,-exclude])


hc<-hclust(expressed.dist)
num.clusters<-45
memb <- cutree(hc, k = num.clusters)
coloursR <-  rainbow(num.clusters)


########### annotations for dendrogram:
memb.ann<-memb[hc$order] # order as in heatmap
locations<-match(labels(memb.ann),ann[,"Array_Address_Id"])
#LL<- as.character(ann[locations,"Entrez_Gene_ID"])    # for gene ids
LL<- as.character(ann[locations,"Symbol"])     # for gene symbols
# LL<- as.character(ann[locations,"Accession"])
##############################################################################################
14_C 14_J
 "AS::03_A" ,  "AS::03_C" ,  "AS::03_E"
,"AS::03_G",   "AS::03_I"  , "AS::03_K"
 ,"AS::10_A",   "AS::10_C"  , "AS::10_E"
 ,"AS::10_G",  "AS::10_I" , "AS::10_K"
 ,"AS::14_A",  "AS::14_C" , "AS::14_E"
 ,"AS::14_G",  "AS::14_I",  "AS::14_K"
 ,"Nor::10_J", "Nor::10_F", "Nor::10_B"
 ,"Nor::03_J", "Nor::14_D" ,"Nor::14_F"
 ,"Nor::10_D", "Nor::10_L", "Nor::03_F"
 ,"Nor::14_J", "Nor::14_H" ,"Nor::03_B"
,"Nor::03_D", "Nor::10_H", "Nor::14_L"
 ,"Nor::03_L", "Nor::14_B" ,"Nor::03_H"           )
 
some.labels<-c(rep("AS",18),rep("Control",18)) 


## bases on PCA so expression level
#heatmap(as.matrix(selDataMatrix[gene.list,-exclude]),col=cm.colors(16),RowSideColors=coloursR[genes.cl$cluster],ColSideColors=coloursC[colour.posns],main="Differentially expressed genes",xlab="Class")

# rainbow(n, s = 1, v = 1, start = 0, end = max(1,n - 1)/n,
#         gamma = 1, alpha = 1)
# heat.colors(n, alpha = 1)
# terrain.colors(n, alpha = 1)
# topo.colors(n, alpha = 1)
# cm.colors(n, alpha = 1)
#
##bases on tree cut
par(font=2,font.lab=2,font.axis=2)
heatmap(as.matrix(selDataMatrix[gene.list,-exclude]),labRow="",labCol=some.labels[-exclude],col=topo.colors(8),RowSideColors=coloursR[memb],ColSideColors=coloursC[colour.posns],main="Differentially expressed genes",xlab="Sample",margins=c(6,6))

gene.names<-as.character(data[,"name"])


gene.names[54]<-"AL713633"
heatmap(as.matrix(selDataMatrix[gene.list,-exclude]),col=heat.colors(8),labCol=some.labels[-exclude],labRow=gene.names,ColSideColors=coloursC[colour.posns],main="70 PAM genes",xlab="Sample",margins=c(6,6))

heatmap(as.matrix(selDataMatrix[gene.list,-exclude]),col=heat.colors(8),labCol=some.labels[-exclude],labRow="",ColSideColors=coloursC[colour.posns],main="Differentially expressed genes",xlab="Sample",margins=c(6,6))

heatmap(as.matrix(selDataMatrix[gene.list,-exclude]),col=heat.colors(8),labCol=some.labels[-exclude],labRow="",ColSideColors=coloursC[colour.posns],xlab="Sample",margins=c(6,6))

plot(hclust(as.matrix(selDataMatrix[gene.list,-exclude]) ))


jpeg(filename = "Rplot%03d.jpeg",
          width = 480, height = 480, units = "px",
          pointsize = 12, quality = 75, bg = "white", res = NA, ...,
          type = c("cairo", "Xlib", "quartz"), antialias)


  RowSideColors:
 ##heatmap with gene symbols
 rownames(diff_expressed_korn)<-diff_expressed_korn_SYM
 heatmap(as.matrix(diff_expressed_korn),col=cm.colors(16),RowSideColors=colours[korn_genes.cl$cluster],main="Heat map: Class and Genes (KORN)",xlab="Class")
  ### purple is higher expressed
  
  ########################)####################################
 patients<-hclust(dist(t(data),method="euclidean"))
kinase<-hclust(dist(dat),method="euclidean"))
 weights_pep<-apply(data,1,sum)
 data["P1a","AF::CCP+ 2"]<-"0.1"
 data[2,12]<-1.3
 data[6,19]<-1.2

 weights<-c(rep(200,12),rep(100,9),rep(50,2),rep(0,4))
  weights[grep("2",colnames(data))]<-200
   weights[grep("1",colnames(data))]<-100
   weights[grep("0",colnames(data))]<-50             appl
   weights[grep("- 0",colnames(data))]<-0

kinase_dd<- as.dendrogram(kinases)
patients_dd <- as.dendrogram(patients)

##tests
plot(patients_dd)
plot(kinase_dd)

patients.reorder <- reorder(patients_dd,c(1:6,27:7))  # note sure about this from
                                                        # plot(patients_dd) & order.dendrogram(patients_dd)
plot(patients.reorder)                                  # inspect to fix order

#patients.reorder <- rev(reorder(patients_dd, weights,agglo.FUN= mean))
kinase.reorder <- reorder(kinase_dd, weights_pep,agglo.FUN= mean)


# breaks <- seq(from = 0, to = max(data), length = max(data)+1)
# breaks <- seq(from = 0, to = 4, length = 4+1)

#heatmap(as.matrix(data),Colv=patients_dd,Rowv=kinase_dd,main="Heat map: Total Antigen",ylab="Kinase")

colorss.ori<-colorss  # colorss<-colorss.ori
#colorss<-bluered(max(data)*4)
color = colorRampPalette(c("blue","red"))
colorss<-color(max(data)*4) #  color = colorRampPalette(c("blue","red"))


################# for CIT vs other
colorss<-colorss[c(-2,-3,-4,-5,-6)]
data<-t(apply(data,1,as.numeric))
colnames(data)<- acc
order<-c("P11a",  "P11b", "P12a", "P12b", "P1a", "P1b", "P2a", "P2b",    "P3a", "P3b")
order<-c("P11a", "P12a", "P11b",  "P12b", "P1a", "P1b", "P2a", "P2b",    "P3a", "P3b")
data<-data[order,]

heatmap.2(as.matrix(data),Colv=patients_dd,col=colorss,dendrogram="column",key=TRUE,keysize=1.47,scale="none", symkey=FALSE, density.info="none", trace="none",Rowv=FALSE,main="Peptides",ylab="Peptide",margins=c(7,5),xlab="Patient :: CCP Status :: Number of Alleles")

heatmap.2(as.matrix(data),Colv=FALSE,col=colorss,dendrogram="none",key=TRUE,keysize=1.47,scale="none", symkey=FALSE, density.info="none", trace="none",Rowv=FALSE,main="Peptides",ylab="Peptide",margins=c(7,5),xlab="Patient :: CCP Status :: Number of Alleles")
###########################################

