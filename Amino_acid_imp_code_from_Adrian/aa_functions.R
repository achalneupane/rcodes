
do_aa_assoc_biallelic <- function
(
 pos = pos,
 ids = ids,
 bim = bim,
 dosage.file = dosage.file,
 allele.freq.thr = allele.freq.thr,
 fam = fam,
 PCA.covariates = PCA.covariates,
 do.cond = do.cond,
 cond.data = cond.data
 ) {

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
    aa.line = data.frame(
      POS=pos,
      TYPE="AA",
      PVAL=NA,
      Z=NA,
      DF=NA,
      NOTE=paste("MAF=",round(maf,4),sep=''))
  } else {
    assoc.data = data.frame(
      ID=fam$V1,
      PHEN=fam$PHEN,
      A1=t(data[,c(4:ncol(data))]))
        
    assoc.data = assoc.data[assoc.data$ID %in% PCA.covariates$ID,]
    assoc.data = merge(assoc.data,PCA.covariates)
        
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
    
    aa.line = data.frame(
      POS=pos,
      TYPE="AA",
      PVAL=lrtest.pval,
      Z=lrtest.z,
      DF=lrtest.df,
      NOTE=paste("Alleles:",a1,a2,sep='')
      )
  }

  return(aa.line) 
                                  
}

do_aa_assoc_multiallelic <- function
(
 pos = pos,
 ids = ids,
 bim = bim,
 dosage.file = dosage.file,
 fam = fam,
 allele.freq.thr = allele.freq.thr,
 PCA.covariates = PCA.covariates,
 do.cond = do.cond,
 cond.data = cond.data
 ) {

  pos.aa.types = sapply(strsplit(ids,split='_'),function(x) return(x[5]))
  ids.to.keep = nchar(pos.aa.types) == 1

  ids = ids[ids.to.keep]
  row.idx = which(bim$V2 %in% ids)
      
  ## make sure ids are coded consequentially
  if( (max(row.idx) - min(row.idx)) != (length(ids)-1) ) {
    stop("may have problems getting the correct data! ",i,"\n")
  }

  tmp.file = tempfile()
  cmd = paste("sed -n ",min(row.idx),",",max(row.idx),
    "p ",dosage.file," > ",tmp.file,sep='')
  system(cmd,intern=FALSE)
  data = read.table(tmp.file,header=F,sep='\t',
    colClasses=c(rep("character",3),rep("numeric",nrow(fam))))
  system(paste("rm ",tmp.file,sep=''))
  colsums = colSums(data[,4:ncol(data)])
  allelefreqs = rowSums(data[,4:ncol(data)])/(2*nrow(fam))


  ### check dosages
  if( mean(colsums) > 2.1 | mean(colsums) < 1.9 ) {
    stop("check colsums (aa analysis) ",i,"\n")
  }
  if( sum(allelefreqs) > 1.1 | sum(allelefreqs) < 0.9 ) {
    stop("check rowsums (aa analysis) ",i,"\n")
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

  ids2 = sapply(strsplit(ids2,split='_'),function(x) return(x[5]))
  colnames(data2) = c("ID","PHEN",ids2)
      
  data2 = data2[,-c(2+which.max(allelefreqs))]
  ### defined reference allele use minor alleles in fitting
  ids2 = ids2[-which.max(allelefreqs)]

  ### removed too much causes error momomorphic
  if(ncol(data2) <= 3) {
    aa.line = data.frame(
      POS=pos,
      TYPE="AA",
      PVAL=NA,
      Z=NA,
      DF=NA,
      NOTE="almost_monomorphic")
  } else {
    data2 = data2[data2$ID %in% PCA.covariates$ID,]
    data2 = merge(data2,PCA.covariates,by=1)
        
    form0 = paste("PHEN ~ ",
      paste(paste("PC",1:10,sep=''),collapse=" + "),sep='')

    if( do.cond ) {          
      n.cond = ncol(cond.data) - 1
      form0 = paste(form0," + ",
        paste(paste("COND",1:n.cond,sep=''),collapse=" + "),sep='')
      data2 = merge(data2,cond.data)          
    }

    ### here add all minor alleles
    form1 = paste(form0," + ",
      paste(ids2,collapse=" + "),sep='')
    
    form0 = as.formula(form0)
    form1 = as.formula(form1)
    
    fit0 = glm(form0,data=data2,family=binomial("logit"))
    fit1 = glm(form1,data=data2,family=binomial("logit"))

    ##### need to get comnbined effect of A_1 + A_2+ A_3
    require(lmtest)
    lrtest.o = lrtest(fit1,fit0)
    lrtest.pval = lrtest.o$Pr[2]
    lrtest.z = lrtest.o$Chisq[2]
    lrtest.df = abs(lrtest.o$Df[2])
        
    aa.line = data.frame(
      POS=pos,
      TYPE="AA",
      PVAL=lrtest.pval,
      Z=lrtest.z,
      DF=lrtest.df,
      NOTE=NA
      )
  }

  return(aa.line)

}

get.hla.class.allele.data <- function
(
 pos = NA,
 allele = NA,
 reso = 2,
 bim = bim,
 dosage.file = dosage.file,
 fam = fam
 ) {
  bim.pos = which( bim$V4 %in% pos )
  bim.pos.ids = bim[bim.pos,2]
  
  ty = sapply(strsplit(bim.pos.ids,split="_"),function(x) return(x[3]))
  keep = nchar(ty) == reso
  ids = bim.pos.ids[keep]

  row.idx = which(bim$V2 %in% bim.pos.ids)
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

  reso.data = cbind(fam[,c(1,7)],hla.gens2[,keep])

  allele.in = reso.data[,allele]
  
  data = data.frame(
    ID=fam$V1,
    PHEN=fam$PHEN,
    A1 = allele.in)

  return(data)
}


get.snp.data <- function
(
 pos = NA,
 bim = bim,
 dosage.file = dosage.file,
 fam = fam
 ) {
  bim.pos = which( bim$V4 %in% pos )
  bim.pos.ids = bim[bim.pos,2]
  bim.pos.ids = bim.pos.ids[grepl("^SNP",bim.pos.ids) |
                            grepl("^rs",bim.pos.ids) |
                            grepl("^1kg",bim.pos.ids)]

  if( length(bim.pos.ids) == 1 ) {
    row.idx = which(bim$V2 %in% bim.pos.ids)
    tmp.file = tempfile()
    cmd = paste("sed -n ",row.idx,"p ",dosage.file," > ",tmp.file,sep='')
    system(cmd,intern=FALSE)
    data = read.table(tmp.file,header=F,sep='\t',
      colClasses=c(rep("character",3),rep("numeric",nrow(fam))))
    system(paste("rm ",tmp.file,sep=''))

    data2 = data.frame(
      ID=fam$V1,
      PHEN=fam$PHEN,
      A1=t(data[,c(4:ncol(data))])
      )

  } else {

    stop("Multiallelic SNP. Not yet implemented\n")

  }

  return(data2)

}


  fam.file = "../imputed_data/AS_MHC_RESULTS/all_dosage_files/all_groups.fam"
bim.file = "../imputed_data/AS_MHC_RESULTS/GROUP4/ankyl_spond.08102012.group4_IMPUTED.bim"
dosage.file = "../imputed_data/AS_MHC_RESULTS/all_dosage_files/all_groups.dosage"

get.aa.data <- function
(
 pos = NA,
 bim = bim,
 dosage.file = dosage.file,
 fam = fam
 ) {

  bim.pos = which( bim$V4 %in% pos )
  bim.pos.ids = bim[bim.pos,2]
  bim.pos.ids = bim.pos.ids[grepl("^AA",bim.pos.ids)]

  ## need to add and if statement if length of ids is 1...
  if( length(bim.pos.ids) == 1 ) {
    row.idx = which(bim$V2 %in% bim.pos.ids)
    tmp.file = tempfile()
    cmd = paste("sed -n ",row.idx,"p ",dosage.file," > ",tmp.file,sep='')
    system(cmd,intern=FALSE)
    data = read.table(tmp.file,header=F,sep='\t',
      colClasses=c(rep("character",3),rep("numeric",nrow(fam))))
    system(paste("rm ",tmp.file,sep=''))

    data2 = data.frame(
      ID=fam$V1,
      PHEN=fam$PHEN,
      A1=t(data[,c(4:ncol(data))])
      )
      
  } else {
  
    types = sapply(strsplit(bim.pos.ids,split='_'),function(x) return(x[5]))
    ids.to.keep = nchar(types) == 1
    ids = bim.pos.ids[ids.to.keep]
    types = types[ids.to.keep]

    row.idx = which(bim$V2 %in% ids)
    tmp.file = tempfile()
    cmd = paste("sed -n ",min(row.idx),",",max(row.idx),
      "p ",dosage.file," > ",tmp.file,sep='')
    system(cmd,intern=FALSE)
    data = read.table(tmp.file,header=F,sep='\t',
      colClasses=c(rep("character",3),rep("numeric",nrow(fam))))

    data2 = t(data[,4:ncol(data)])
    data2 = cbind(fam[,c(1,7)],data2)
    colnames(data2) = c("ID","PHEN",types)

  }

  return(data2)
}



get.aa.data.res <- function
(
 pos = NA,
 res = NA,
 bim = bim,
 dosage.file = dosage.file,
 fam = fam
 ) {

  bim.pos = which( bim$V4 %in% pos )
  bim.pos.ids = bim[bim.pos,2]
  bim.pos.ids = bim.pos.ids[grepl("^AA",bim.pos.ids)]

  ## need to add and if statement if length of ids is 1...
  if( length(bim.pos.ids) == 1 ) {
    row.idx = which(bim$V2 %in% bim.pos.ids)
    tmp.file = tempfile()
    cmd = paste("sed -n ",row.idx,"p ",dosage.file," > ",tmp.file,sep='')
    system(cmd,intern=FALSE)
    data = read.table(tmp.file,header=F,sep='\t',
      colClasses=c(rep("character",3),rep("numeric",nrow(fam))))
    system(paste("rm ",tmp.file,sep=''))

    data2 = data.frame(
      ID=fam$V1,
      PHEN=fam$PHEN,
      A1=t(data[,c(4:ncol(data))])
      )
      
  } else {
  
    types = sapply(strsplit(bim.pos.ids,split='_'),function(x) return(x[5]))
    ids.to.keep = nchar(types) == 1
    ids = bim.pos.ids[ids.to.keep]
    types = types[ids.to.keep]

    row.idx = which(bim$V2 %in% ids)
    tmp.file = tempfile()
    cmd = paste("sed -n ",min(row.idx),",",max(row.idx),
      "p ",dosage.file," > ",tmp.file,sep='')
    system(cmd,intern=FALSE)
    data = read.table(tmp.file,header=F,sep='\t',
      colClasses=c(rep("character",3),rep("numeric",nrow(fam))))

    data2 = t(data[,4:ncol(data)])
    data2 = cbind(fam[,c(1,7)],data2)
    colnames(data2) = c("ID","PHEN",types)

  }

  return(data2)
}



do_aa_assoc_multiallelic_return_fit1 <- function
(
 pos = pos,
 ids = ids,
 bim = bim,
 dosage.file = dosage.file,
 fam = fam,
 allele.freq.thr = allele.freq.thr,
 PCA.covariates = PCA.covariates,
 do.cond = do.cond,
 cond.data = cond.data
 ) {

  pos.aa.types = sapply(strsplit(ids,split='_'),function(x) return(x[5]))
  ids.to.keep = nchar(pos.aa.types) == 1

  ids = ids[ids.to.keep]
  row.idx = which(bim$V2 %in% ids)
      
  ## make sure ids are coded consequentially
  if( (max(row.idx) - min(row.idx)) != (length(ids)-1) ) {
    stop("may have problems getting the correct data! ",i,"\n")
  }

  tmp.file = tempfile()
  cmd = paste("sed -n ",min(row.idx),",",max(row.idx),
    "p ",dosage.file," > ",tmp.file,sep='')
  system(cmd,intern=FALSE)
  data = read.table(tmp.file,header=F,sep='\t',
    colClasses=c(rep("character",3),rep("numeric",nrow(fam))))
  system(paste("rm ",tmp.file,sep=''))
  colsums = colSums(data[,4:ncol(data)])
  allelefreqs = rowSums(data[,4:ncol(data)])/(2*nrow(fam))

  if( mean(colsums) > 2.1 | mean(colsums) < 1.9 ) {
    stop("check colsums (aa analysis) ",i,"\n")
  }
  if( sum(allelefreqs) > 1.1 | sum(allelefreqs) < 0.9 ) {
    stop("check rowsums (aa analysis) ",i,"\n")
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

  ids2 = sapply(strsplit(ids2,split='_'),function(x) return(x[5]))
  colnames(data2) = c("ID","PHEN",ids2)
      
  data2 = data2[,-c(2+which.max(allelefreqs))]
  ids2 = ids2[-which.max(allelefreqs)]
      
  if(ncol(data2) <= 3) {
    aa.line = data.frame(
      POS=pos,
      TYPE="AA",
      PVAL=NA,
      Z=NA,
      DF=NA,
      NOTE="almost_monomorphic")

    out = NA
  } else {
    data2 = data2[data2$ID %in% PCA.covariates$ID,]
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
        
    aa.line = data.frame(
      POS=pos,
      TYPE="AA",
      PVAL=lrtest.pval,
      Z=lrtest.z,
      DF=lrtest.df,
      NOTE=NA
      )

    out = fit1
  }

  return(out)

}



do_aa_assoc_multiallelic_return_fit1_all_aa <- function
(
 pos = pos,
 ids = ids,
 bim = bim,
 dosage.file = dosage.file,
 fam = fam,
 allele.freq.thr = allele.freq.thr,
 PCA.covariates = PCA.covariates,
 do.cond = do.cond,
 cond.data = cond.data
 ) {

  pos.aa.types = sapply(strsplit(ids,split='_'),function(x) return(x[5]))
  ids.to.keep = nchar(pos.aa.types) == 1

  ids = ids[ids.to.keep]
  row.idx = which(bim$V2 %in% ids)
      
  ## make sure ids are coded consequentially
  if( (max(row.idx) - min(row.idx)) != (length(ids)-1) ) {
    stop("may have problems getting the correct data! ",i,"\n")
  }

  tmp.file = tempfile()
  cmd = paste("sed -n ",min(row.idx),",",max(row.idx),
    "p ",dosage.file," > ",tmp.file,sep='')
  system(cmd,intern=FALSE)
  data = read.table(tmp.file,header=F,sep='\t',
    colClasses=c(rep("character",3),rep("numeric",nrow(fam))))
  system(paste("rm ",tmp.file,sep=''))
  colsums = colSums(data[,4:ncol(data)])
  allelefreqs = rowSums(data[,4:ncol(data)])/(2*nrow(fam))

  if( mean(colsums) > 2.1 | mean(colsums) < 1.9 ) {
    stop("check colsums (aa analysis) ",i,"\n")
  }
  if( sum(allelefreqs) > 1.1 | sum(allelefreqs) < 0.9 ) {
    stop("check rowsums (aa analysis) ",i,"\n")
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

  ids2 = sapply(strsplit(ids2,split='_'),function(x) return(x[5]))
  colnames(data2) = c("ID","PHEN",ids2)
      
  # data2 = data2[,-c(2+which.max(allelefreqs))]
  # ids2 = ids2[-which.max(allelefreqs)]
      
  if(ncol(data2) <= 3) {
    aa.line = data.frame(
      POS=pos,
      TYPE="AA",
      PVAL=NA,
      Z=NA,
      DF=NA,
      NOTE="almost_monomorphic")

    out = NA
  } else {
    data2 = data2[data2$ID %in% PCA.covariates$ID,]
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
        
    aa.line = data.frame(
      POS=pos,
      TYPE="AA",
      PVAL=lrtest.pval,
      Z=lrtest.z,
      DF=lrtest.df,
      NOTE=NA
      )

    out = fit1
  }

  return(out)

}



do_aa_assoc_multiallelic_return_fit1_all_aa_all_alleles_vs_one <- function
(
 pos = pos,
 ids = ids,
 bim = bim,
 dosage.file = dosage.file,
 fam = fam,
 allele.freq.thr = allele.freq.thr,
 PCA.covariates = PCA.covariates,
 do.cond = do.cond,
 cond.data = cond.data
 ) {

  pos.aa.types = sapply(strsplit(ids,split='_'),function(x) return(x[5]))
  ids.to.keep = nchar(pos.aa.types) == 1

  ids = ids[ids.to.keep]
  row.idx = which(bim$V2 %in% ids)
      
  ## make sure ids are coded consequentially
  if( (max(row.idx) - min(row.idx)) != (length(ids)-1) ) {
    stop("may have problems getting the correct data! ",i,"\n")
  }

  tmp.file = tempfile()
  cmd = paste("sed -n ",min(row.idx),",",max(row.idx),
    "p ",dosage.file," > ",tmp.file,sep='')
  system(cmd,intern=FALSE)
  data = read.table(tmp.file,header=F,sep='\t',
    colClasses=c(rep("character",3),rep("numeric",nrow(fam))))
  system(paste("rm ",tmp.file,sep=''))
  colsums = colSums(data[,4:ncol(data)])
  allelefreqs = rowSums(data[,4:ncol(data)])/(2*nrow(fam))

  if( mean(colsums) > 2.1 | mean(colsums) < 1.9 ) {
    stop("check colsums (aa analysis) ",i,"\n")
  }
  if( sum(allelefreqs) > 1.1 | sum(allelefreqs) < 0.9 ) {
    stop("check rowsums (aa analysis) ",i,"\n")
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

  ids2 = sapply(strsplit(ids2,split='_'),function(x) return(x[5]))
  colnames(data2) = c("ID","PHEN",ids2)
      
  # data2 = data2[,-c(2+which.max(allelefreqs))]
  # ids2 = ids2[-which.max(allelefreqs)]
      
  if(ncol(data2) <= 3) {
    aa.line = data.frame(
      POS=pos,
      TYPE="AA",
      PVAL=NA,
      Z=NA,
      DF=NA,
      NOTE="almost_monomorphic")

    out = NA
  } else {
    data2 = data2[data2$ID %in% PCA.covariates$ID,]
    data2 = merge(data2,PCA.covariates,by=1)

    aa.res = data.frame()
    
    for( aa in ids2 ) {

      have.aa  = data2[,aa]
      have.not = data2[,ids2[ ! ids2 %in% aa ]]
      have.not2 = rowSums(have.not)

      data3 = data2[,colnames(data2)[ ! colnames(data2) %in% ids2]]
      data3$AA = have.aa
      data3$AAnot = have.not2

      form0 = paste("PHEN ~ ",
        paste(paste("PC",1:10,sep=''),collapse=" + "),sep='')
      form1 = paste(form0," + ","AA",sep='')

      form0 = as.formula(form0)
      form1 = as.formula(form1)
    
      fit0 = glm(form0,data=data3,family=binomial("logit"))
      fit1 = glm(form1,data=data3,family=binomial("logit"))

      fit1.aa.effect = coefficients(fit1)['AA']
      fit1.aa.pval = summary(fit1)$coefficients['AA',4]
           
      require(lmtest)
      lrtest.o = lrtest(fit1,fit0)
      lrtest.pval = lrtest.o$Pr[2]
      lrtest.z = lrtest.o$Chisq[2]
      lrtest.df = abs(lrtest.o$Df[2])

      aa.line = data.frame(
        POS=pos,
        TYPE="AA",
        AA = aa,
        AA_BETA = round(fit1.aa.effect,3),
        AA_PVAL = fit1.aa.pval,
        PVAL=lrtest.pval,
        Z=lrtest.z,
        DF=lrtest.df,
        NOTE=NA
        )
      
      aa.res = rbind(aa.res,aa.line)
      
    }

    out = fit1
  }

  return(out)

}


get_aa_biallelic_freqs_info<- function
(
 pos = pos,
 ids = ids,
 bim = bim,
 dosage.file = dosage.file,
 allele.freq.thr = allele.freq.thr,
 fam = fam,
 PCA.covariates = PCA.covariates
 ) {

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
  assoc.data = data.frame(
      ID=fam$V1,
      PHEN=fam$PHEN,
      A1=t(data[,c(4:ncol(data))]))
        
  assoc.data = assoc.data[assoc.data$ID %in% PCA.covariates$ID,]
  cas = subset(assoc.data,PHEN==1)
  con = subset(assoc.data,PHEN==0)
  
  a1.freq.cas = sum(cas$A1)/(2*nrow(cas))
  a1.freq.con = sum(con$A1)/(2*nrow(con))

  a2.freq.cas = 1-a1.freq.cas
  a2.freq.con = 1-a1.freq.con

  out = data.frame(
    POS = pos,
    alleles = paste(a1,a2,sep="/"),
    Freq_Cas = paste(round(a1.freq.cas,2),round(a2.freq.cas,2),sep='/'),
    Freq_Con = paste(round(a1.freq.con,2),round(a2.freq.con,2),sep='/')
    )

  return(out)
}



get_aa_multiallelic_freqs_info <- function
(
 pos = pos,
 ids = ids,
 bim = bim,
 dosage.file = dosage.file,
 fam = fam,
 allele.freq.thr = allele.freq.thr,
 PCA.covariates = PCA.covariates
 ) {

  pos.aa.types = sapply(strsplit(ids,split='_'),function(x) return(x[5]))
  ids.to.keep = nchar(pos.aa.types) == 1

  ids = ids[ids.to.keep]
  row.idx = which(bim$V2 %in% ids)
  tmp.file = tempfile()
  cmd = paste("sed -n ",min(row.idx),",",max(row.idx),
    "p ",dosage.file," > ",tmp.file,sep='')
  system(cmd,intern=FALSE)
  data = read.table(tmp.file,header=F,sep='\t',
    colClasses=c(rep("character",3),rep("numeric",nrow(fam))))
  system(paste("rm ",tmp.file,sep=''))
  colsums = colSums(data[,4:ncol(data)])
  allelefreqs = rowSums(data[,4:ncol(data)])/(2*nrow(fam))
  data2 = t(data[,4:ncol(data)])
  data2 = cbind(fam[,c(1,7)],data2)   
  ids2 = ids
  ids2 = sapply(strsplit(ids2,split='_'),function(x) return(x[5]))
  colnames(data2) = c("ID","PHEN",ids2)
  data2 = data2[data2$ID %in% PCA.covariates$ID,]
  
  cas = subset(data2,PHEN==1)
  con = subset(data2,PHEN==0)

  cas.freqs = colSums(cas[,3:ncol(cas)])/(2*nrow(cas))
  con.freqs = colSums(con[,3:ncol(con)])/(2*nrow(con))

  out = data.frame(
    POS = pos,
    alleles = paste(ids2,collapse='/'),
    Freq_Cas = paste(round(cas.freqs,2),collapse='/'),
    Freq_Con = paste(round(con.freqs,2),collapse='/')
    )

  return(out)
}
