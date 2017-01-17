

pca.file<-"MND_shell2.evecs"
fam.file<-"ch_mnd_650_f_ld_f.fam"
working.directory<-"/media/Bioinform-D/Research/shellfish/shellfish"


setwd(working.directory)  ## change working directory



the.pca<-read.table(pca.file,header=F,skip=0,fill=TRUE,stringsAsFactors=FALSE)
dim(the.pca)
the.pca<-t(the.pca)
dim(the.pca)


the.fam<-read.table(fam.file,header=F,skip=0,fill=TRUE,stringsAsFactors=FALSE)
dim(the.fam)

the.fam[1:5,]
the.pca[1:5,]

the.pca<-cbind(the.fam[,1],the.pca)

the.pca[1:5,]
colnames(the.pca)<-c("Sample",paste("evec",1:(dim(the.pca)[2]-1),sep=":"))

write.table(the.pca,file="new.pca.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


y.axis<-"evec:1"
x.axis<-"evec:2"

plot(the.pca[,x.axis],the.pca[,y.axis],type="p")
