########################################################################################################################
########################################################################################################################
########################################################################################################################


geno.all<-{}

target<-"AOGCHBMFam"
the.samples<-sample.sheet.full[sample.sheet.full[,"SampleProject"]==target,"ParticipantCode"]
the.samples
genotypes<-a.indel[,paste(the.samples,"GT",sep=".")]
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")
summary.geno[1:5,]

geno.all<-summary.geno

target<-"AOGCHGM"
the.samples<-sample.sheet.full[sample.sheet.full[,"SampleProject"]==target,"ParticipantCode"]
the.samples
genotypes<-a.indel[,paste(the.samples,"GT",sep=".")]
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")
summary.geno[1:5,]

geno.all<-cbind(geno.all,summary.geno)

target<-"AOGCHBM"
the.samples<-sample.sheet.full[sample.sheet.full[,"SampleProject"]==target,"ParticipantCode"]
the.samples
genotypes<-a.indel[,paste(the.samples,"GT",sep=".")]
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")
summary.geno[1:5,]

geno.all<-cbind(geno.all,summary.geno)

target<-"AOGCLBM"
the.samples<-sample.sheet.full[sample.sheet.full[,"SampleProject"]==target,"ParticipantCode"]
the.samples
genotypes<-a.indel[,paste(the.samples,"GT",sep=".")]
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),target,sep=".")
summary.geno[1:5,]

geno.all<-cbind(geno.all,summary.geno)

geno.all[1:5,]

#insert about line after col 48:
setwd(analysis.dir)
write.table(cbind(a.indel[,1:6],geno.all),file=paste(project.files[ichr],"geno.all.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
