########################################################################################################################
########################################################################################################################
########################################################################################################################


geno.all<-{}



# all.genotypes1<-filtered.genotype.summary(indels,all.samples,prefix="",suffix="",1,0,1.1,1.1,0,0,0) ## no filtering checked ok
all.genotypes<-full.genotype.summary(indels,all.samples,prefix="",suffix="")
geno.all<-cbind(all.genotypes$"summary.geno.group",all.genotypes$"summary.depths.group")
print("got unfiltered")

## sum(all.genotypes$"summary.geno.group"[,3]!=all.genotypes1$"summary.geno.group"[,3])
## sum(all.genotypes$"summary.depths.group"[,2]!=all.genotypes1$"summary.depths.group"[,2])

all.genotypes<-filtered.genotype.summary(indels,all.samples,prefix="",suffix=".filtered",30,0.1,0.90,0.05,0.95,10,5) # most often used
print("got filtered UQCCG")
geno.all<-cbind(geno.all,all.genotypes$"summary.geno.group",all.genotypes$"summary.depths.group")

all.genotypes<-filtered.genotype.summary(indels,all.samples,prefix="",suffix=".25p.16x.filtered",1,0.25,0.75,0.05,0.95,16,16) # most often used
print("got filtered peter mac")
geno.all<-cbind(geno.all,all.genotypes$"summary.geno.group",all.genotypes$"summary.depths.group")


#insert about line after col 48:
## dim(geno.all)
## geno.all[1:5,]

remove.from.all.samples.local<-c("GT","AD","DP")
to.remove.all.local<-expand.labels.to.samples(remove.from.all.samples.local,all.possible.samples)
to.remove.all.local<-unique(to.remove.all.local)

indels<-indels[,colnames(indels)[!(colnames(indels) %in% to.remove.all.local)]]

print(dim(indels))

setwd(analysis.dir)
write.table(cbind(indels,geno.all),file=paste(gsub(".analysis.txt$","",project.files[ichr],perl=TRUE),"controls_mac.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(cbind(indels,geno.all[,1:11]),file=paste(gsub(".analysis.txt$","",project.files[ichr],perl=TRUE),"controls.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
