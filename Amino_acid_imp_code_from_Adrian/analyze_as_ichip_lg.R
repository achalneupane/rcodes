

rm(list=ls())
options(stringsAsFactors=F)

#####################
## Initialize Data ##
#####################

host = system("hostname",intern=T)

if( identical(host,"bioinform06") ) {
  wd <- "/home/acortes/Projects2/Immunochip_Utrecht/LG_ANALYSIS_PLINK"
} else {
  wd <- "/home/acortes/Projects/Ichip_Utrecht/Analysis/LG_ANALYSIS_PLINK"
}
setwd(wd)
source("scripts/aa_functions.R")

do.snp.tests = TRUE
do.aa.tests = TRUE
do.hla.tests = TRUE
do.cond = TRUE



## Samples to remove based on dosage ##
samples.out.dosage <- read.table("../imputed_data/AS_MHC_RESULTS/QC/dosage_samples_outliers.txt")

## Load the PCA data and those samples ##
## to be included in the analysis      ##
famW.file = paste("../pca_data/AS_European_after_hapmap_nofp_diag_check2b.fam",sep='')
famW <- read.table(famW.file)
famPCA = read.table(paste("../pca_data/European_all_nohighLD_trimset.fam",sep=''))
famPCA[grepl("/",famPCA$V1),'V1'] = gsub("/","_",famPCA[grepl("/",famPCA$V1),'V1'])
evecs = read.table(paste("../pca_data/European_PCA_no_assoc.evecs",sep=''))
evecs = t(evecs)
PCAdata = cbind(famPCA,evecs)
PCAdata = PCAdata[,c(1,7:16)]
colnames(PCAdata) <- c("ID",paste("PC",1:10,sep=''))

## Non-cofounding covariates
PCA.covariates = PCAdata[ PCAdata$ID %in% famW$V1 , ]

## remove samples with outlying dosage
PCA.covariates = PCA.covariates[ ! PCA.covariates$ID %in% samples.out.dosage$V1 , ]

fam.file = "../imputed_data/AS_MHC_RESULTS/all_dosage_files/all_groups.fam"
bim.file = "../imputed_data/AS_MHC_RESULTS/GROUP4/ankyl_spond.08102012.group4_IMPUTED.bim"
dosage.file = "../imputed_data/AS_MHC_RESULTS/all_dosage_files/all_groups.dosage"

allele.freq.thr = 0.01

fam = read.table(fam.file,header=F)
fam$PHEN = fam$V6 - 1
bim = read.table(bim.file,header=F)

which.rs = grepl("^rs",bim$V2)
which.1kg = grepl("^1kg",bim$V2)
which.AA = grepl("^AA",bim$V2)
which.SNP = grepl("^SNP",bim$V2)
which.HLA = grepl("^HLA",bim$V2)
which.SNPS = which.rs | which.1kg | which.SNP

otherids = c('imm_6_30184807', '6_30223299', 'INS_C_295x296_31345779_VLAVLA', 
	'INS_B_294x295_31430920', 'imm_6_31976744', 'gw_060011', 'gw_060012',
	'INS_DQB1_226x227_32736002_PQGPPPAG')

which.notclass = bim$V2 %in% otherids

yougot = which.rs | which.1kg | which.AA | which.SNP | which.HLA | which.notclass

all.pos = unique(bim$V4)
all.pos = all.pos[order(all.pos)]

####################
## Set Parameters ##
####################

cond1 = get.aa.data(
  pos = 31432180,
  bim = bim,
  dosage.file = dosage.file,
  fam = fam)

controls = subset(cond1,PHEN==0)
means = colMeans(controls[,c(3:8)])
col2rm = which.max(means) + 2
cond1 = cond1[,-col2rm]

## cond1 = cond1[,-2]
## colnames(cond1) = c("ID","COND1","COND2","COND3","COND4","COND5")

cond.data = cond1[,-2]
colnames(cond.data) = c("ID","COND1","COND2","COND3","COND4","COND5")

file.out = "results/cond1_97omnibus2.Rdata"

############################################


#####################################################
## Do biallelic analysis on all markers imputed.   ##
## Akin the suggested methods in the SNP2HLA site. ##
#####################################################

source("scripts/biallelic_tests.R")

#####################################################

###############################
## Do position omnibus tests ##
###############################


aa.results = data.frame()
snp.results = data.frame()
hla.results = data.frame()

print.count = 100

## all.pos = all.pos[ all.pos >= 31000000 & all.pos <= 31800000 ]

for( i in 1:length(all.pos) ) {

  pos = all.pos[i]
  pos.bim = bim[bim$V4 %in% pos,]

  ## exceptions
  if( all(pos.bim$V2 %in% c("rs2844750","rs17354619")) ) {
    pos.bim = pos.bim[1,]
  }
  if( all(pos.bim$V2 %in% c("rs3819284","SNP_B_31430745")) ) {
    pos.bim = pos.bim[1,]
  }
  if( all(pos.bim$V2 %in% c("rs1059615","SNP_DRB1_32657541")) ) {
    pos.bim = pos.bim[1,]
  }
  if( all(pos.bim$V2 %in% c("rs1048087","SNP_DQA1_32717264")) ) {
    pos.bim = pos.bim[1,]
  }

  ## check if you have amino acids
  have.aa = any( pos.bim$V2 %in% bim$V2[which.AA] )
  have.snp = any( pos.bim$V2 %in% bim$V2[which.SNPS] )
  have.hla = any( pos.bim$V2 %in% bim$V2[which.HLA] )

  if( do.hla.tests ) {
    if( have.hla ) {

      ids = pos.bim$V2[pos.bim$V2 %in% bim$V2[which.HLA]]
      row.idx = which(bim$V2 %in% ids)

      ## make sure ids are coded consequentially
      if( (max(row.idx) - min(row.idx)) != (length(ids)-1) ) {
        stop("may have problems getting the correct data! ",i,"\n")
      }
      tmp.file = tempfile()
      cmd = paste("sed -n ",min(row.idx),",",max(row.idx),"p ",
        dosage.file," > ",tmp.file,sep='')
      system(cmd,intern=FALSE)
      data = read.table(tmp.file,header=F,sep='\t',
        colClasses=c(rep("character",3),rep("numeric",nrow(fam))))
      system(paste("rm ",tmp.file,sep=''))

      hla.gens = data[,4:ncol(data)]
      hla.gens2 = t(hla.gens)

      colnames(hla.gens2) = data[,1]

      ids.codes = sapply(strsplit(colnames(hla.gens2),split='_'),function(x) return(x[3]))
      are.2.digit = nchar(ids.codes) == 2
      are.4.digit = nchar(ids.codes) == 4

      two.digit.data = cbind(fam[,c(1,7)],hla.gens2[,are.2.digit])
      four.digit.data = cbind(fam[,c(1,7)],hla.gens2[,are.4.digit])
      
      rowsums.2 = rowSums(two.digit.data[,3:ncol(two.digit.data)])
      allelefreqs.2 = colSums(two.digit.data[,3:ncol(two.digit.data)])/(2*nrow(fam))

      rowsums.4 = rowSums(four.digit.data[,3:ncol(four.digit.data)])
      allelefreqs.4 = colSums(four.digit.data[,3:ncol(four.digit.data)])/(2*nrow(fam))

      ## iterate through the 2 digit alleles.
      for( j in 1:length(which(are.2.digit)) ) {
        
        allele.an = colnames(two.digit.data)[2+j]
        allele.in = two.digit.data[,allele.an]
        allele.out = two.digit.data[,-c(1,2,
          which( colnames(two.digit.data) %in% allele.an ))]
        allele.out = rowSums(allele.out)
        assoc.data = data.frame(
          ID=two.digit.data[,1],
          PHEN=two.digit.data$PHEN,
          AL1=allele.in,
          AL2=allele.out)

        assoc.data = assoc.data[ assoc.data$ID %in% PCA.covariates$ID,]
        assoc.data = merge(assoc.data,PCA.covariates)

        form0 = paste("PHEN ~ ",
          paste(paste("PC",1:10,sep=''),collapse=" + "),sep='')

        if( do.cond ) {
          n.cond = ncol(cond.data) - 1
          form0 = paste(form0," + ",
            paste(paste("COND",1:n.cond,sep=''),collapse=" + "),sep='')
          assoc.data = assoc.data = merge(assoc.data,cond.data)
        }
        
        form1 = paste(form0," + AL1",sep='')
        
        form0 = as.formula(form0)
        form1 = as.formula(form1)

        fit0 = glm(form0,data=assoc.data,family=binomial("logit"))
        fit1 = glm(form1,data=assoc.data,family=binomial("logit"))

        require(lmtest)
        lrtest.o = lrtest(fit1,fit0)
        lrtest.pval = lrtest.o$Pr[2]
        lrtest.z = lrtest.o$Chisq[2]
        lrtest.df = lrtest.o$Df[2]
        
        hla.line = data.frame(
          POS=pos,
          TYPE="HLA2",
          PVAL=lrtest.pval,
          Z=lrtest.z,
          DF=lrtest.df,
          NOTE=allele.an)
        
        hla.results = rbind(hla.results,hla.line)
        
      }

      ## iterate through the 4 digit alleles.
      for( j in 1:length(which(are.4.digit)) ) {
        allele.an = colnames(four.digit.data)[2+j]
        allele.in = four.digit.data[,allele.an]
        allele.out = four.digit.data[,
          -c(1,2,which( colnames(four.digit.data) %in% allele.an ))]
        allele.out = rowSums(allele.out)
        assoc.data = data.frame(
          ID=four.digit.data[,1],
          PHEN=four.digit.data$PHEN,
          AL1=allele.in,
          AL2=allele.out)

        assoc.data = assoc.data[ assoc.data$ID %in% PCA.covariates$ID,]
        assoc.data = merge(assoc.data,PCA.covariates)
        
        form0 = paste("PHEN ~ ",
          paste(paste("PC",1:10,sep=''),collapse=" + "),sep='')

        if( do.cond ) {
          n.cond = ncol(cond.data) - 1
          form0 = paste(form0," + ",
            paste(paste("COND",1:n.cond,sep=''),collapse=" + "),sep='')
          assoc.data = merge(assoc.data,cond.data)
        }          
        
        form1 = paste(form0," + AL1",sep='')

        form0 = as.formula(form0)
        form1 = as.formula(form1)

        fit0 = glm(form0,data=assoc.data,family=binomial("logit"))
        fit1 = glm(form1,data=assoc.data,family=binomial("logit"))

        require(lmtest)
        lrtest.o = lrtest(fit1,fit0)
        lrtest.pval = lrtest.o$Pr[2]
        lrtest.z = lrtest.o$Chisq[2]
        lrtest.df = abs(lrtest.o$Df[2])
        
        hla.line = data.frame(
          POS=pos,
          TYPE="HLA4",
          PVAL=lrtest.pval,
          Z=lrtest.z,
          DF=lrtest.df,
          NOTE=allele.an)
        
        hla.results = rbind(hla.results,hla.line)

      }
    }
  }
  
  if( do.snp.tests ) {
    if( have.snp ) {
      ids = pos.bim$V2[pos.bim$V2 %in% bim$V2[which.SNPS]]
      
      if(length(ids) == 1) {

        row.idx = which(bim$V2 %in% ids)
        tmp.file = tempfile()
        cmd = paste("sed -n ",row.idx,"p ",dosage.file," > ",tmp.file,sep='')
        system(cmd,intern=FALSE)
        data = read.table(tmp.file,header=F,sep='\t',
          colClasses=c(rep("character",3),rep("numeric",nrow(fam))))
        system(paste("rm ",tmp.file,sep=''))
        a1 = bim[bim$V2 %in% ids,5]
        a2 = bim[bim$V2 %in% ids,6]
        
        maf = sum(data[,4:ncol(data)])/(2*nrow(fam))
        
        if( maf < allele.freq.thr ) {
          snp.line = data.frame(
            POS=pos,
            TYPE="SNP",
            PVAL=NA,
            Z=NA,
            DF=NA,
            NOTE=paste("MAF=",round(maf,4),sep=''))
        } else {
          fam$A1 = t(data[,c(4:ncol(data))])
          assoc.data = PCA.covariates
          assoc.data$PHEN = fam[ match(assoc.data$ID,fam$V1), 'PHEN']
          assoc.data$A1 = fam[ match(assoc.data$ID,fam$V1), 'A1']

          form0 = paste("PHEN ~ ",
            paste(paste("PC",1:10,sep=''),collapse=" + "),sep='')

          if( do.cond ) {
            n.cond = ncol(cond.data) - 1

            form0 = paste(form0," + ",
              paste(paste("COND",1:n.cond,sep=''),collapse=" + "),sep='')

            assoc.data = merge(assoc.data,cond.data)
            
          }
          
          form1 = paste(form0," + A1",sep='')
          
          form0 = as.formula(form0)
          form1 = as.formula(form1)
          
          fit0 = glm(form0,data=assoc.data,family=binomial("logit"))
          fit1 = glm(form1,data=assoc.data,family=binomial("logit"))
                
          require(lmtest)
          lrtest.o = lrtest(fit1,fit0)
          lrtest.pval = lrtest.o$Pr[2]
          lrtest.z = lrtest.o$Chisq[2]
          lrtest.df = abs(lrtest.o$Df[2])
                
          snp.line = data.frame(
            POS=pos,
            TYPE="SNP",
            PVAL=lrtest.pval,
            Z=lrtest.z,
            DF=lrtest.df,
            NOTE=paste("Alleles:",a1,a2,sep='')
            )
        }
      } else { ## A multiallelic SNP

        pos.snp.types = sapply(strsplit(ids,split='_'),function(x) return(x[4]))
        ids.to.keep = nchar(pos.snp.types) == 1
        ids = ids[ids.to.keep]
        row.idx = which(bim$V2 %in% ids)
        tmp.file = tempfile()
        cmd = paste("sed -n ",min(row.idx),",",
          max(row.idx),"p ",dosage.file," > ",tmp.file,sep='')
        system(cmd,intern=FALSE)
        data = read.table(tmp.file,header=F,sep='\t',
          colClasses=c(rep("character",3),rep("numeric",nrow(fam))))
        system(paste("rm ",tmp.file,sep=''))
        colsums = colSums(data[,4:ncol(data)])
        allelefreqs = rowSums(data[,4:ncol(data)])/(2*nrow(fam))
        
        if( mean(colsums) > 2.1 | mean(colsums) < 1.9 ) {
          stop("check colsums (snp analysis) ",i,"\n")
        }
        if( sum(allelefreqs) > 1.1 | sum(allelefreqs) < 0.9 ) {
          stop("check rowsums (snp analysis) ",i,"\n")
        }
              
        data2 = t(data[,4:ncol(data)])
        data2 = cbind(fam[,c(1,7)],data2)

        ## remove those alleles below the frequency threshold
        which.remove.freq = which( allelefreqs < allele.freq.thr )
        ids2 = ids
        if(length(which.remove.freq)>0) {
          data2 = data2[,-(2+which.remove.freq)]				
          ids2 = ids[-which.remove.freq]
          allelefreqs = allelefreqs[-which.remove.freq]
        }
        ids2 = sapply(strsplit(ids2,split='_'),function(x) return(x[4]))
        colnames(data2) = c("ID","PHEN",ids2)
                           
        if(ncol(data2) <= 3) {
          snp.line = data.frame(
            POS=pos,
            TYPE="SNP",
            PVAL=NA,
            Z=NA,
            DF=NA,
            NOTE="almost_monomorphic")
        } else {

          data2 = data2[,-c(2+which.max(allelefreqs))]
          data2 = data2[ data2$ID %in% PCA.covariates$ID ,]

          ids2 = ids2[-which.max(allelefreqs)]
          data2 = merge(data2,PCA.covariates,by=1)
          
          form0 = paste("PHEN ~ ",
            paste(paste("PC",1:10,sep=''),collapse=" + "),sep='')

          if( do.cond ) {

            n.cond = ncol(cond.data) - 1

            form0 = paste(form0," + ",
              paste(paste("COND",1:n.cond,sep=''),collapse=" + "),sep='')

            data2 = merge(data2,cond.data)
            
          }
          
          form1 = paste(form0," + ",
            paste(ids2,collapse=" + "),sep='')
          
          form0 = as.formula(form0)
          form1 = as.formula(form1)
          
          fit0 = glm(form0,data=data2,family=binomial("logit"))
          fit1 = glm(form1,data=data2,family=binomial("logit"))
          
          require(lmtest)
          lrtest.o = lrtest(fit1,fit0)
          lrtest.pval = lrtest.o$Pr[2]
          lrtest.z = lrtest.o$Chisq[2]
          lrtest.df = abs(lrtest.o$Df[2])
          
          snp.line = data.frame(
            POS=pos,
            TYPE="SNP",
            PVAL=lrtest.pval,
            Z=lrtest.z,
            DF=lrtest.df,
            NOTE=NA)
        }
      }
            
      snp.results = rbind(snp.results,snp.line)
    }
  }
  
  if( do.aa.tests ) {
    if( have.aa ) {

      ids = pos.bim$V2[pos.bim$V2 %in% bim$V2[which.AA]]

      if(length(ids)==1) { 

        aa.line = do_aa_assoc_biallelic(
          pos = pos,
          ids = ids,
          bim = bim,
          dosage.file = dosage.file,
          allele.freq.thr = allele.freq.thr,
          fam = fam,
          PCA.covariates = PCA.covariates,
          do.cond = do.cond,
          cond.data = cond.data
          )
        
      } else {

        aa.line = do_aa_assoc_multiallelic(
          pos = pos,
          ids = ids,
          bim = bim,
          dosage.file = dosage.file,
          fam = fam,
          allele.freq.thr = allele.freq.thr,
          PCA.covariates = PCA.covariates,
          do.cond = do.cond,
          cond.data = cond.data
          )

      }
      aa.results = rbind(aa.results,aa.line)
    }
  }

  if(exists(as.character(substitute(snp.line)))) rm(snp.line)
  if(exists(as.character(substitute(aa.line)))) rm(aa.line)
  if(exists(as.character(substitute(hla.line)))) rm(hla.line)


  if( ! i %% print.count ) {
    cat(".")
  }
  if( ! i %% 1000 ) {
    cat("\n")
  }
  
}

save.image(file=file.out)








