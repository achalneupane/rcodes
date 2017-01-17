# This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
# Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of 
# Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000  

# NOTE: return code of -1.0 signals an error condition


## genotype_counts_r.txt:

##           HET   HOM1  HOM2
## MARKER_1  100    3    5
## MARKER_2   57   14   50
## MARKER_3   31   32   51
## MARKER_4   47    3    5
## MARKER_5  150   32   55
## MARKER_6  122    7   32
## MARKER_7   99    3   14
## MARKER_8  146   13   54
## MARKER_9  177  100   57
## MARKER_10 184   57  155


## Calling code from run_marker_tests.r:

## data <- read.table("genotype_counts_r.txt")                
## n_markers <- 10
## p_values <-rep(0, n_markers)
   
## test_markers <- function()
##    {
##    for (i in 1:n_markers)
##       {
##       hets  <- data[i, 1]
##       hom_1 <- data[i, 2]
##       hom_2 <- data[i, 3]   
##       p_values[i] <- SNPHWE(hets, hom_1, hom_2)
##       print(p_values[i])
##       }
##    }
   
## test_markers()




 HWExactMat.pval<- function(X, ...) 
{
        out <- HWExactMat(X, ...)
    return(out$pvalvec)
}


getHWE <- function(geno,num.cores=8) #  geno<-summary.geno.extra[,"GENO.Control"] geno<-geno[1:5]
{
library(doMC)
registerDoMC(cores=num.cores)
## print(num.cores)
#  http://cran.r-project.org/web/packages/HardyWeinberg/HardyWeinberg.pdf
  require("HardyWeinberg")
geno<-strsplit(geno,split=",")
AA<-as.numeric(unlist(lapply(geno,function(x) x[1])))
BB<-as.numeric(unlist(lapply(geno,function(x) x[3])))
AB<-as.numeric(unlist(lapply(geno,function(x) x[2])))
x<-cbind(AA,AB,BB)
## x<-rbind(x,c(0,1,914),c(0,399,550),c(0,110,577),c(0,266,569))
## print(num.bits)
## print(num.cores)
 if(dim(x)[1]<1000){num.bits<-1}else{num.bits<-num.cores} 
 values<-foreach(x.bit=iter(x, by='row',chunksize=as.integer(dim(x)[1]/num.bits)), .combine=c,.multicombine=TRUE,.inorder=TRUE) %dopar% HWExactMat.pval(x.bit)

 return(values)
}

#HWAlltests(x[20,])



get.ONE.HWE <- function(geno) #  geno<-summary.geno.extra[,"GENO.Control"] 
{
geno<-strsplit(geno,split=",")
obs_hom1<-as.numeric(unlist(lapply(geno,function(x) x[1])))
obs_hom2<-as.numeric(unlist(lapply(geno,function(x) x[3])))
obs_hets<-as.numeric(unlist(lapply(geno,function(x) x[2])))


   if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0)
      return(-1.0)

   # total number of genotypes
   N <- obs_hom1 + obs_hom2 + obs_hets
   
   # rare homozygotes, common homozygotes
   obs_homr <- min(obs_hom1, obs_hom2)
   obs_homc <- max(obs_hom1, obs_hom2)

   # number of rare allele copies
   rare  <- obs_homr * 2 + obs_hets

   # Initialize probability array
   probs <- rep(0, 1 + rare)

   # Find midpoint of the distribution
   mid <- floor(rare * ( 2 * N - rare) / (2 * N))
   if ( (mid %% 2) != (rare %% 2) ) mid <- mid + 1

   probs[mid + 1] <- 1.0
   mysum <- 1.0

   # Calculate probablities from midpoint down 
   curr_hets <- mid
   curr_homr <- (rare - mid) / 2
   curr_homc <- N - curr_hets - curr_homr

   while ( curr_hets >=  2)
      {
      probs[curr_hets - 1]  <- probs[curr_hets + 1] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0)  * (curr_homc + 1.0))
      mysum <- mysum + probs[curr_hets - 1]

      # 2 fewer heterozygotes -> add 1 rare homozygote, 1 common homozygote
      curr_hets <- curr_hets - 2
      curr_homr <- curr_homr + 1
      curr_homc <- curr_homc + 1
      }    

   # Calculate probabilities from midpoint up
   curr_hets <- mid
   curr_homr <- (rare - mid) / 2
   curr_homc <- N - curr_hets - curr_homr
   
   while ( curr_hets <= rare - 2)
      {
      probs[curr_hets + 3] <- probs[curr_hets + 1] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
      mysum <- mysum + probs[curr_hets + 3]
         
      # add 2 heterozygotes -> subtract 1 rare homozygtote, 1 common homozygote
      curr_hets <- curr_hets + 2
      curr_homr <- curr_homr - 1
      curr_homc <- curr_homc - 1
      }    
 
    # P-value calculation
    target <- probs[obs_hets + 1]

    #plo <- min(1.0, sum(probs[1:obs_hets + 1]) / mysum)

    #phi <- min(1.0, sum(probs[obs_hets + 1: rare + 1]) / mysum)

    # This assignment is the last statement in the fuction to ensure 
    # that it is used as the return value
    p <- min(1.0, sum(probs[probs <= target])/ mysum)
    }


HWE <- function(obs_hets, obs_hom1, obs_hom2)
   {
   if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0)
      return(-1.0)

   # total number of genotypes
   N <- obs_hom1 + obs_hom2 + obs_hets
   
   # rare homozygotes, common homozygotes
   obs_homr <- min(obs_hom1, obs_hom2)
   obs_homc <- max(obs_hom1, obs_hom2)

   # number of rare allele copies
   rare  <- obs_homr * 2 + obs_hets

   # Initialize probability array
   probs <- rep(0, 1 + rare)

   # Find midpoint of the distribution
   mid <- floor(rare * ( 2 * N - rare) / (2 * N))
   if ( (mid %% 2) != (rare %% 2) ) mid <- mid + 1

   probs[mid + 1] <- 1.0
   mysum <- 1.0

   # Calculate probablities from midpoint down 
   curr_hets <- mid
   curr_homr <- (rare - mid) / 2
   curr_homc <- N - curr_hets - curr_homr

   while ( curr_hets >=  2)
      {
      probs[curr_hets - 1]  <- probs[curr_hets + 1] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0)  * (curr_homc + 1.0))
      mysum <- mysum + probs[curr_hets - 1]

      # 2 fewer heterozygotes -> add 1 rare homozygote, 1 common homozygote
      curr_hets <- curr_hets - 2
      curr_homr <- curr_homr + 1
      curr_homc <- curr_homc + 1
      }    

   # Calculate probabilities from midpoint up
   curr_hets <- mid
   curr_homr <- (rare - mid) / 2
   curr_homc <- N - curr_hets - curr_homr
   
   while ( curr_hets <= rare - 2)
      {
      probs[curr_hets + 3] <- probs[curr_hets + 1] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
      mysum <- mysum + probs[curr_hets + 3]
         
      # add 2 heterozygotes -> subtract 1 rare homozygtote, 1 common homozygote
      curr_hets <- curr_hets + 2
      curr_homr <- curr_homr - 1
      curr_homc <- curr_homc - 1
      }    
 
    # P-value calculation
    target <- probs[obs_hets + 1]

    #plo <- min(1.0, sum(probs[1:obs_hets + 1]) / mysum)

    #phi <- min(1.0, sum(probs[obs_hets + 1: rare + 1]) / mysum)

    # This assignment is the last statement in the fuction to ensure 
    # that it is used as the return value
    p <- min(1.0, sum(probs[probs <= target])/ mysum)
    }
