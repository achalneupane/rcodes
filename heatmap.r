
library(gplots)
getwd()
#setwd( "/media/old-scratch/media/Bioinform-D/Research/MicroArray/Ran's AS gene expression study")
load("/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/clustering.RData")
exclude<-c(14,28)  #14_C and 14J excluded from analysis
data<-as.matrix(selDataMatrix[gene.list,-exclude])
colnames(data)<-some.labels[-exclude]

# HCLUST method: the agglomeration method to be used. This should be (an
#           unambiguous abbreviation of) one of '"ward"', '"single"',
#           '"complete"', '"average"', '"mcquitty"', '"median"' or
#           '"centroid"'.
### centered correlation average linkage
##cor use: method = c("pearson", "kendall", "spearman")

################### BASIC heatmaps for paper
plot(hclust(dist(t(data))), ColSideColors=coloursC[colour.posns]   )
 
patients<-hclust(dist(t(data)))
patients_dd <- as.dendrogram(patients)
heatm<-heatmap(as.matrix(selDataMatrix[gene.list,-exclude]),col=heat.colors(8),labCol=some.labels[-exclude],labRow="",ColSideColors=coloursC[colour.posns],xlab="Sample",margins=c(6,6))

##tests
plot(patients_dd)
patients.reorder <- reorder(patients_dd,heatm$colInd)  # note sure about this from
                                                        # plot(patients_dd) & order.dendrogram(patients_dd)
plot(patients_dd)
plot(patients.reorder)

jpeg("AS_heatmap_dendro.jpeg",quality=100) ## print


dev.off()
##########################################
mydist<-as.dist((1-cor(data,method="pearson")))
plot(hclust(mydist))     # ## CORRELATION

patients<-hclust(mydist,method="average") ## avergae
plot(hclust(mydist,method="average"),main="Hierarchial Clustering Dendrogram with Filtered data")






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
##############################################################################################         )
 
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

######################### different versions

heatmap(as.matrix(selDataMatrix[gene.list,-exclude]),col=heat.colors(8),labCol=some.labels[-exclude],labRow=gene.names,ColSideColors=coloursC[colour.posns],main="70 PAM genes",xlab="Sample",margins=c(6,6))

heatmap(as.matrix(selDataMatrix[gene.list,-exclude]),col=heat.colors(8),labCol=some.labels[-exclude],labRow="",ColSideColors=coloursC[colour.posns],main="Differentially expressed genes",xlab="Sample",margins=c(6,6))

heatmap(as.matrix(selDataMatrix[gene.list,-exclude]),col=heat.colors(8),labCol=some.labels[-exclude],labRow="",ColSideColors=coloursC[colour.posns],xlab="Sample",margins=c(6,6))

plot(hclust(as.matrix(selDataMatrix[gene.list,-exclude]) ))


jpeg(filename = "Rplot%03d.jpeg",
          width = 480, height = 480, units = "px",
          pointsize = 12, quality = 75, bg = "white", res = NA, ...,
          type = c("cairo", "Xlib", "quartz"), antialias)













####################### BELOW FOR PCA
############## Col colours 
coloursC <-  rainbow(2)
colour.posns<-class.labels[-exclude]
colour.posns[ colour.posns=="AS"]<-1
colour.posns[ colour.posns=="Nor"]<-2
colour.posns<-as.numeric(colour.posns)    

############## Row colours use if have too many genes



#### DO A PCA 
genes.pca <- prcomp(selDataMatrix[gene.list,-exclude],scale = TRUE)


### THIS IS A SCREE PLOT THAT TELLS YOU NOW MANT PCA's to use
dev.off()
genes.pca.var <- round(genes.pca$sdev^2 / sum(genes.pca$sdev^2)*100,2)
plot(c(1:length(genes.pca.var)),genes.pca.var,type="b",xlab="# components",ylab="% variance",main="Scree plot",col="blue",cex.lab=2)

fig.prefix<-"figure_1"
savePlot(filename="scree_plot",type="png")
savePlot(filename=paste(fig.prefix,".madhatt",.png",sep=''),type="png")
savePlot(filename=paste(fig.prefix,".madhatt".jpeg",sep=''),type="jpeg")
savePlot(filename=paste(fig.prefix,".madhatt.tiff",sep=''),type="tiff")
savePlot(filename=paste(fig.prefix,".madhatt.bmp",sep=''),type="bmp")
dev.print(svg,paste(fig.prefix,".madhatt.svg",sep=''))





 ###first 2 component have most or varaiblity so use K-means on frist two components
 ##from Heat map 7-8 clusters available
dev.off()
num.clusters<-16
coloursR <-  rainbow(num.clusters)
genes.cl<-kmeans(genes.pca$x[,1:2],centers=num.clusters,iter.max=500) #Do kmeans


################ Spectral clustering:
plot(range(genes.pca$x[,1]),range(genes.pca$x[,2]),xlab="PCA1",ylab="PCA2",main="Spectral Clustering of Genes")

#points(sam_genes.pca$x[,1],sam_genes.pca$x[,2],col=colours[sam_genes.cl$cluster],pch=sam_genes.cl$cluster,cex=1.0,bg=colours[sam_genes.cl$cluster]) #colours and symbols
text(genes.pca$x[,1],genes.pca$x[,2],label=rownames(genes.pca$x),col=coloursR[genes.cl$cluster],cex=0.75) #colours and symbols
text(genes.pca$x[,1],genes.pca$x[,2],label=as.character(genes.cl$cluster),col=coloursR[genes.cl$cluster],cex=0.75)
 median(abs(gene.means[labels(genes.cl$cluster[genes.cl$cluster==1])]))


##############################################################################################
