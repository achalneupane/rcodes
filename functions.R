get_N_data_splits <- function(data=data,trai.p=trai.p,N=1,phen.col='PHENOTYPE') {
  
  data_list <- lapply(1:N,function(x) return(split.data.samerate(data=data,trai.p=trai.p,
                                                                 phen.col=phen.col)))
  return(data_list)

}
split.data.samerate <- function(data,trai.p=0.8,phen.col='PHENOTYPE') {

  rate <- length(which(data[,phen.col]==1))/nrow(data)
  
  data.pos <- data[data[,phen.col]==1,]
  data.neg <- data[data[,phen.col]==0,]
  
  data.pos <- split.data.uniformly(data.pos,trai.p=trai.p)
  data.neg <- split.data.uniformly(data.neg,trai.p=trai.p)
  
  training <- rbind(data.pos$training,data.neg$training)
  test <- rbind(data.pos$test,data.neg$test)
  
  rate.trai <- length(which(training[,phen.col]==1))/nrow(training)
  rate.test <- length(which(test[,phen.col]==1))/nrow(test)
  
  return( list( training=training, test=test))
  
}
split.data.uniformly <- function(data,trai.p=0.8) {

  order <- order(runif(nrow(data)))
  data2 <- data[order,]
  
  training <- data2[ 1:floor(nrow(data)*trai.p) , ]
  test <- data2[ (floor(nrow(data)*trai.p)+1):nrow(data), ]
  
  return( list( training=training, test=test))
}
