
#clusters.wanted.subset<-c("FANC_complex.all", "Clinical" ,"BRCA.Genes","RAD51.Paralogues","BLM.Complex","Checkpoint.Proteins","citric","Citric_final","BLM.Complex_AND_Checkpoint","FANCD2_minimal_mono_ubi","MLH_cluster","Richard", "C1Alpha","C1Delta","C1Beta","C2","C3","C4","C5","MRC","MRC_IDH","C1","NDUFA","NDUFB","NDUFS","NDUFCV")
## to.unwind<-c("FGFR3") #, "MCM7", "RNPC3")
## interesting.gene<-c("BCL2L12","CCDC61","STK19","KNSTRN","TRHDE","FREM2","EBNA1BP2","PHACTR3","DCLK1","LRRIQ1","PHACTR3","CSMD3")
## ## to.unwind<-c("FANC_complex.all") # to.unwind<-meta.results.burden[8,"gene"]
##  to.unwind<-c("BCL2L12","CCDC61","STK19","KNSTRN","TRHDE","FREM2","EBNA1BP2","PHACTR3","DCLK1","LRRIQ1","PHACTR3","CSMD3") #"test.set"

## to.unwind<-c("TP53","NOTCH1","NOTCH2","FAT4","STK19","ISX","TRHDE","ARHGAP35","PREX1","KL","PIK3CA","KNSTRN","BCL2L12")
## to.unwind<-c("TP53","NOTCH1","NOTCH2","ATM","ACD","ASIP","BAP1","CASP8","CCND1","CDK4","MC1R","MITF","MTAP","MX2","OCA2","PARP1","PLA2G6","POT1","SLC45A2","TERF2IP","TERT","TYR","TYRP1","VDR","BCL2L12","KNSTRN","ISX","CDKN2A","BCL2L11","STK19","FJX1","TRHDE")

#to.unwind<-c("CDKN2A") #,"TCEA1","POLR2A","CTDP1")
## to.unwind<-c("Clinical")

## to.unwind<-c("FANCD2_minimal_mono_ubi")
## to.unwind<-c("BLM.Complex")
## to.unwind<-c("BLM.Complex_AND_Checkpoint")


#to.unwind<-c(clusters.wanted.subset)
#to.unwind<-c("FANC_complex.all","Citric_final","citric","Clinical")


#sum(!(the.genes %in% wanted.genes))
## to.unwind<-c("Citric_final")
## to.unwind<-c("citric")

## to.unwind<-c("Ubin.proteo","lipid_raft","caveolae","Citric")
#to.unwind<-c(clusters.wanted[!(clusters.wanted %in% c("Ubin.proteo","lipid_raft","caveolae","Checkpoint_extendedx1","Checkpoint_extendedx2"))])
#grep(to.unwind,meta.results.burden[,"gene"])
#to.unwind
#to.unwind.name<-to.unwind[1]
#to.unwind.name<-"TOP_550_LENGTH_CONTROL"
#match(net,meta.results.burden[,"gene"])
# to.unwind.name<-"SYNON_test"
# to.unwind.name<-"Pathways"
# to.unwind.name<-"ALL_significant"
# to.unwind.name<-"ALL_significant"



to.unwind<-c(meta.results.burden[,"gene"])# ,meta.results.skatO[1:the.top,"gene"])

to.unwind.name<-"ALL"


snpinfo.ex<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,]
loci<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,"Name"] # this is IDH1 not IDH1 in cluster # are the snp.names
loci<-unique(loci)
the.genes<-unique(snpinfo.ex[,"cluster"])
the.genes<-unique(snpinfo.ex[,"gene"])
the.genes<-the.genes[!(the.genes %in% clusters.wanted)]

sort(the.genes) #245 ### if used a cluster name need to do back up to (**) the.genes<-c(the.genes,"STAG2")

################ use below to get data for genes

## the.genes.burden<-meta.results.burden[meta.results.burden[,"gene"] %in% the.genes,]

## the.genes.burden
## write.table(the.genes.burden,file=paste(to.unwind.name,"conponents:","Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

## the.genes.burden<-meta.results.skatO[meta.results.skatO[,"gene"] %in% the.genes,]

## write.table(the.genes.burden,file=paste(paste(to.unwind.name,collapse="."),"conponents:","SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



########### single point only
## length(loci)
## meta.results.burden[1:5,]
## loci<-meta.results.burden[1:550,"Name"]
###############


dim(genotypes)
genotypes[1:5,1:5]
genotypes.ex<-genotypes[,loci]


dim(genotypes.ex)
genotypes.ex[is.na(genotypes.ex)]<-0
dim(genotypes.ex)






snpinfo.ex<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,]
dim(snpinfo.ex)
dim(genotypes.ex)
dim(pheno)
snpinfo.ex[1:5,]


########### single point only
## snpinfo.ex<-snpinfo[snpinfo[,"Name"] %in% loci,]
## dim(snpinfo.ex)
## meta.results.burden.ex<-meta.results.burden[1:550,]
             ###############
# summary.geno.extra[loci,]
#high.missing[loci,]
#sum(are.in.repeats[loci])

    
# qual[loci,]
## snpinfo[1:5,]
## qual[1:5,c("FILTER_PASS", "FILTER_100" )]

cohort.seq.ex <- skatCohort(genotypes.ex, formula, SNPInfo = snpinfo.ex, data=pheno,aggregateBy = "Name",verbose=FALSE)
## meta.results.skat.ex<-skatMeta(cohort.seq,SNPInfo = snpinfo)
meta.results.burden.ex<-burdenMeta(cohort.seq.ex,wts=1,mafRange = c(0,1),SNPInfo = snpinfo.ex,aggregateBy = "Name")
meta.results.burden.ex[1:5,]
pheno[1:5,]
dim(meta.results.burden.ex)
dim(genotypes.ex)

cohort.seq.test <- skatCohort(genotypes.ex, formula, SNPInfo = snpinfo.ex, data=pheno,aggregateBy = "cluster",verbose=FALSE)

meta.results.burden.test<-burdenMeta(cohort.seq.test,wts=1,mafRange = c(0,1),SNPInfo = snpinfo.ex,aggregateBy = "cluster")
#meta.results.burden.test

## meta.results.skat.ex<-skatMeta(cohort.seq,SNPInfo = snpinfo)
#meta.results.skatO.test<-skatOMeta(cohort.seq.test,burden.wts =1,SNPInfo = snpinfo.ex,aggregateBy="cluster")
#meta.results.skatO.test
figure<- match(loci,key)

#genotypes.PD<-a.indel[figure, c("LPH-001-27_PD.GT",paste(pheno.ori[pheno.ori[,"PD"],"SAMPLE"],".GT",sep="")) ]
PDs<-pheno.ori[pheno.ori[,"AML-Child"] | pheno.ori[,"Asian-AML-Child"] | pheno.ori[,"Asian-AML"]  | pheno.ori[,"AML-NotDiagnosis-Child"] | pheno.ori[, "Asian-AML-NotDiagnosis-Child"],"SAMPLE"]
genotypes.PD<-a.indel[figure, c(paste(PDs,".GT",sep="")) ]
genotypes.PD<-t(genotypes.PD)

genotypes.PD[genotypes.PD=="NA"]<-NA
genotypes.PD[genotypes.PD=="0/0"]<-0
genotypes.PD[genotypes.PD=="0/1"]<-1
genotypes.PD[genotypes.PD=="1/1"]<-2

rownames(genotypes.PD)<-gsub(".GT","",rownames(genotypes.PD))

dim(genotypes.ex)
dim(genotypes.PD)

options(max.print=200)
muts.in.PD<-apply(genotypes.PD,2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
muts.in.cases<-apply(genotypes.ex[pheno[,"AML"],],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
muts.in.controls<-apply(genotypes.ex[pheno[,"Control"],],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})


 controls<- paste(pheno[pheno[,"Control"],"SAMPLE"],".GT",sep="")
## a.indel[figure,controls]
##   table(a.indel[figure,controls][2,])  

## muts.in.cases<-apply(genotypes.ex[pheno[,"SampleProject"]==1,],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
## muts.in.controls<-apply(genotypes.ex[pheno[,"SampleProject"]==0,],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})




########################################################
check<-16

quality.cases<-rep("",times=length(loci))
quality.controls<-rep("",times=length(loci))
quality.PD<-rep("",times=length(loci))

depth.cases<-rep("",times=length(loci))
depth.fad.cases<-rep("",times=length(loci))
dup.cases<-rep("",times=length(loci))

depth.controls<-rep("",times=length(loci))
depth.fad.controls<-rep("",times=length(loci))
dup.controls<-rep("",times=length(loci))

depth.PD<-rep("",times=length(loci))
depth.fad.PD<-rep("",times=length(loci))
dup.PD<-rep("",times=length(loci))

a.indel.sub<-a.indel[figure,]
a.indel.stats.sub<-a.indel.stats[figure,]


somatic.matrix.desc.full.sub<-somatic.matrix.desc.full[figure,]
somatic.matrix.p.full.sub<-somatic.matrix.p.full[figure,]

somatic.cases<-rep("",times=length(loci))
somatic.PD<-rep("",times=length(loci))

somatic.p.cases<-rep("",times=length(loci))
somatic.p.PD<-rep("",times=length(loci))



# a.indel.stats.sub[1:5,1:20]
# dim(a.indel.sub)

check<-1
for(check in 1:length(loci)){
# print(check)
#check<-"chr11:130066457:130066457:-:A:indel"
# posn<-grep(loci[check],key)
posn<-check


if(muts.in.PD[check]!=""){
#the.gt<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GT",sep=".")
the.gq<-paste(unlist(strsplit(muts.in.PD[check],split=",")),"GQ",sep=".")
quality.PD[check]<-paste(a.indel.sub[posn,the.gq],collapse=",")

the.ad<-paste(unlist(strsplit(muts.in.PD[check],split=",")),"AD",sep=".")
depth.PD[check]<-paste(a.indel.sub[posn,the.ad],collapse=";")

the.fad<-paste(unlist(strsplit(muts.in.PD[check],split=",")),"FAD",sep=".")
depth.fad.PD[check]<-paste(a.indel.stats.sub[posn,the.fad],collapse=";")

the.dup<-paste(unlist(strsplit(muts.in.PD[check],split=",")),"DUP",sep=".")
dup.PD[check]<-paste(a.indel.stats.sub[posn,the.dup],collapse=";")

the.ad.soma<-paste(unlist(strsplit(muts.in.PD[check],split=",")),"GT",sep=".")
somatic.PD[check]<-paste(somatic.matrix.desc.full.sub[posn,the.ad.soma],collapse=";")

the.ad.soma<-paste(unlist(strsplit(muts.in.PD[check],split=",")),"GT",sep=".")
somatic.p.PD[check]<-paste(signif(somatic.matrix.p.full.sub[posn,the.ad.soma],digits=4),collapse=";")


a.indel[posn,the.gq]
## a.indel[posn,the.gt]
## a.indel[posn,the.dp]
}




if(muts.in.cases[check]!=""){
#the.gt<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GT",sep=".")
the.gq<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GQ",sep=".")
quality.cases[check]<-paste(a.indel.sub[posn,the.gq],collapse=",")

the.ad<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"AD",sep=".")
depth.cases[check]<-paste(a.indel.sub[posn,the.ad],collapse=";")

the.fad<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"FAD",sep=".")
depth.fad.cases[check]<-paste(a.indel.stats.sub[posn,the.fad],collapse=";")

the.dup<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"DUP",sep=".")
dup.cases[check]<-paste(a.indel.stats.sub[posn,the.dup],collapse=";")


the.ad.soma<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GT",sep=".")
somatic.cases[check]<-paste(somatic.matrix.desc.full.sub[posn,the.ad.soma],collapse=";")

the.ad.soma<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GT",sep=".")
somatic.p.cases[check]<-paste(signif(somatic.matrix.p.full.sub[posn,the.ad.soma],digits=4),collapse=";")



a.indel[posn,the.gq]
## a.indel[posn,the.gt]
## a.indel[posn,the.dp]
}

if(muts.in.controls[check]!=""){
#the.gt<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"GT",sep=".")
the.gq<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"GQ",sep=".")
quality.controls[check]<-paste(a.indel.sub[posn,the.gq],collapse=",")

the.ad<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"AD",sep=".")
depth.controls[check]<-paste(a.indel.sub[posn,the.ad],collapse=";")

the.fad<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"FAD",sep=".")
depth.fad.controls[check]<-paste(a.indel.stats.sub[posn,the.fad],collapse=";")

the.dup<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"DUP",sep=".")
dup.controls[check]<-paste(a.indel.stats.sub[posn,the.dup],collapse=";")

a.indel[posn,the.gq]
## a.indel[posn,the.gt]
## a.indel[posn,the.dp]
}

} # end check
##########################################################################

                             
#figure
length(figure)
dim(meta.results.burden.ex)
length(muts.in.cases)
length(muts.in.controls)
#pass[figure]
#help[figure,]

toString(colnames(a.indel)[c(1:6,8,11,16,28,7,30,34,35,36,37:42,43,14,32,33)])
colnames(a.indel)[1:60]

 ann.cols<-c("chr","start","end","REF","ALT","TYPE","refGene::type","knownGene::type","Gene.Names","Genes.mentioned.at.ASH","refGene::location","OMIM (Gene::Status::OMIM::description::disease)","Consequence.Embl","Uploaded_variation.Embl","Gene.Embl","Feature.Embl", "Protein_position.Embl", "Amino_acids.Embl" , "ensGene::type","ID::maf","FILTER")# ,"rs.id")


annotations<-a.indel[,ann.cols]
dim(annotations)
dim(help)
dim(summary.geno.extra)
dim(a.indel)
dim(poss.model)
length(quality.cases)
length(figure)
dim(meta.results.burden.ex)
gerp.scores<-a.indel[,"gerp.scores"]
                             
#sum(meta.results.burden.ex[,"gene"]!=loci)
## colnames(a.indel)[1:50]

## key[grep("chr17",key)[1:100]]
## grep("chr17:41197708",key)
## key[grep("10088407",key)]
#out<-cbind(meta.results.burden.ex,a.indel[figure,c(1:6,16,28,7,30,34,37:42,43)],summary.geno.extra[figure,],high.missing[figure,],help[figure,])
## out<-cbind(meta.results.burden.ex,a.indel[figure,c(1:6,16,28,7,30,34,37:42,43,14,32,33)],summary.geno.extra[figure,c("GENO.AML","GENO.Control","GENO.AML.filt","GENO.Control.filt")],high.missing[figure,])
## summary.geno.extra[figure,]
## annotations[figure,]
## help[figure,]

dim(meta.results.burden.ex)
if(!exists("summary.geno.extra.ori")){summary.geno.extra.ori<-summary.geno.extra}
if(!exists("summary.geno.extra.ori")){summary.geno.extra.ori<-summary.geno.extra}
## if(!exists("pass.old")){pass.old<-pass}
## if(!exists("pass.new")){pass.new<-pass}
#out<-cbind(meta.results.burden.ex,a.indel[figure,c(1:6,16,43,28,7,30,34,37:42)],summary.geno.extra[figure,c("GENO.AML","GENO.Control","GENO.AML.filt","GENO.Control.filt")],help[figure,],muts.in.cases,muts.in.controls)
a.functions<-a.indel[,c("PolyPhen.scores","SIFT.scores","PolyPhen.desc","SIFT.desc")]


posns<-match(key,filt.key)
missing<-is.na(posns)
sum(missing)
filt.sub<-filt[posns,]



 filt.sub[figure,c("FILTER_SUMMARY","SUMMARY_CALLED","SUMMARY_NOT_CALLED")]
                             
out<-cbind(meta.results.burden.ex,a.functions[figure,],gerp.scores[figure],annotations[figure,],maf.lt.all[figure,],is.benign.missense[figure],annotations[figure,],summary.geno.extra[figure,colnames(summary.geno.extra)[grep("^GENO",colnames(summary.geno.extra))]], filt.sub[figure,c("FILTER_SUMMARY","SUMMARY_CALLED","SUMMARY_NOT_CALLED")],pass.single[figure],pass.coding[figure],pass.sliding[figure],help[figure,],high.missing.table[figure,],poss.model[figure,],poss.model.lib[figure,],muts.in.cases,somatic.cases,somatic.p.cases,quality.cases,depth.fad.cases,depth.cases,dup.cases,muts.in.PD,somatic.PD,somatic.p.PD,quality.PD,depth.fad.PD,depth.PD,muts.in.controls,quality.controls,depth.fad.controls,depth.controls,dup.controls,summary.geno.extra.ori[figure,colnames(summary.geno.extra.ori)[grep("^GENO",colnames(summary.geno.extra.ori))]])

## out<-cbind(meta.results.burden.ex,a.functions[figure,],gerp.scores[figure],annotations[figure,],maf.lt.all[figure,],is.benign.missense[figure],annotations[figure,],summary.geno.extra[figure,colnames(summary.geno.extra)[grep("^GENO",colnames(summary.geno.extra))]], filt.sub[figure,c("FILTER_SUMMARY","SUMMARY_CALLED","SUMMARY_NOT_CALLED")],pass.lose.filters[figure],pass.0.01.use[figure],pass.0.001.use[figure],pass.0.01.bad.loc[figure],pass.0.001.bad.loc[figure],pass.0.01.old[figure],pass.0.001.old[figure],help[figure,],high.missing.table[figure,],validated.posn[figure],poss.model[figure,],poss.model.lib[figure,],muts.in.cases,somatic.cases,somatic.p.cases,quality.cases,depth.fad.cases,depth.cases,dup.cases,muts.in.PD,somatic.PD,somatic.p.PD,quality.PD,depth.fad.PD,depth.PD,muts.in.controls,quality.controls,depth.fad.controls,depth.controls,dup.controls,summary.geno.extra.ori[figure,colnames(summary.geno.extra.ori)[grep("^GENO",colnames(summary.geno.extra.ori))]])



#all.data[figure,]
#out<-cbind(meta.results.burden.ex,annotations[figure,],muts.in.cases,muts.in.controls)
dim(out)
out[,1:13]



## table(out[,"refGene::location"])
## table(out[,"Consequence.Embl"]) # to.unwind.name<-"IDH"
getwd()
setwd(analysis.dir)
paste(paste(to.unwind,collapse="."))
paste(to.unwind.name,collapse=".")
  paste(paste(to.unwind.name,collapse="."),"GENOTYPE.conponents.","SkatO","clusters",snap.file,"txt",sep=".")

order.by<-order(out[,"p"],decreasing=FALSE)
#enum<-1:dim(meta.results.burden.ex)[1]
out[order.by,][1:10,1:10]
setwd(analysis.dir)
write.table(out[order.by,],file=paste(paste(to.unwind.name,collapse="."),"GENOTYPE.conponents.","Burden","clusters_FINAL",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

getwd()
