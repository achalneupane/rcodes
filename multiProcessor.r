

code.dir<-"/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")



######## below line can be run to demonstrante the process
#foreach(genotypes.bit=iter(genotypes,by='row',chunksize=as.integer(dim(genotypes)[1]/num.bits)), .combine='cbind', .multicombine=TRUE,.inorder=TRUE) %dopar%  print(dim(genotypes.bit))
library(doMC)
num.cores<-5
registerDoMC(cores=num.cores)

###### EXAMPLE 1

#################  that takes matrix genotypes and applies function logscan to each ROW
################ result is a matrix where processed rows are in the same order as the oroginal


if(dim(genotypes)[1]<200){
num.bits<-1; print(paste("Only",dim(genotypes)[1],"genotypes remaining",sep=" ")) }else{num.bits<-num.cores}  # this line just in case the matrix is small

res<-foreach(genotypes.bit=iter(genotypes,by='row',chunksize=as.integer(dim(genotypes)[1]/num.bits)), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% logscan(genotypes.bit,a.fam,"pheno",num.samples)

write.table(res,file=paste("rescale",traits[i],projects[ip],sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


####3 EXAMPLE 2
######################### another example
fil.genotypes<-foreach(a.indel.bit=iter(a.indel,by='row',chunksize=5000 ), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype(a.indel.bit,the.samples.fam,prefix="",suffix="",20,0.2,0.80,0.05,0.95,10,5)




### EXAMPLE 3
################# a more complicated exaple proceccing a matrix all.genotypes$"genotypes", here the iter function operates using the row numbers and the reults are added together..( don't this it need to be this complicated as function requires the same "chunk" from many matrices

Aij<-foreach(genotypes.bit=iter(as.matrix(c(1:(dim(all.genotypes$"genotypes")[1]))),by='row',chunksize=as.integer(dim(all.genotypes$"genotypes")[1]/num.bits)), .combine='+', .multicombine=FALSE,.inorder=FALSE) %dopar% genetic.sample.QC.accumilate(all.genotypes$"genotypes"[as.integer(genotypes.bit),],p[as.integer(genotypes.bit)],position.filter.QC[as.integer(genotypes.bit)],the.samples,on.X[as.integer(genotypes.bit)],parX[as.integer(genotypes.bit)])
