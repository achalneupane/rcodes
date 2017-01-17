
biallelic.tab.res = data.frame()
dosage.con = file(dosage.file,open='r')
j = 0
stay.in.loop = TRUE

while(stay.in.loop) {
  oneline <- readLines(dosage.con,n=1,warn=FALSE)
  j=j+1
  if (j %% 100 == 0) {
    cat ("biallelic number: ", j, "\n")    
  }
   if(length(oneline)==0) {
    stay.in.loop=FALSE
    break;
  }
  
  pos = bim[j,4]
  a1 = bim[j,5]
  a2 = bim[j,6]
  
  data <- strsplit(oneline,split='\t')[[1]]

  maf <- sum(as.numeric(data[4:length(data)]))/(2*nrow(fam))
  
  if( maf < allele.freq.thr ) {
    biallelic.line = data.frame(
      POS=pos,
      TYPE="biallelic",
      PVAL=NA,
      Z=NA,
      DF=NA,
      NOTE=paste("MAF=",round(maf,4),":ID=",bim[j,2],sep=''))
  } else {
    fam$A1 = as.numeric(data[c(4:length(data))])
    assoc.data = PCA.covariates
    assoc.data$PHEN = fam[ match(assoc.data$ID,fam$V1), 'PHEN']
    assoc.data$A1 = fam[ match(assoc.data$ID,fam$V1), 'A1']

    form0 = paste("PHEN ~ ",
      paste(paste("PC",1:4,sep=''),collapse=" + "),sep='')

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
                
    biallelic.line = data.frame(
      POS=pos,
      TYPE="SNP",
      PVAL=lrtest.pval,
      Z=lrtest.z,
      DF=lrtest.df,
      NOTE=paste("Alleles:",a1,a2,":ID=",bim[j,2],sep='')
      )
  }
  
  biallelic.tab.res = rbind(biallelic.tab.res,biallelic.line)
}



