
## sample.labels  ### this are in individual samples in ALL
## indels[1:5,]

names(indiv.quality.type)<-indiv.quality.labs                       
names(indiv.quality)<-indiv.quality.names


###############################
quality.labs<-expand.labels.to.samples(indiv.quality.labs,sample.labels)
quality.names<-expand.labels.to.samples(indiv.quality.names,sample.labels)
quality<-expand.quantity.to.samples(indiv.quality,sample.labels)
quality.type<-expand.quantity.to.samples(indiv.quality.type,sample.labels)
quality.dirn<-expand.quantity.to.samples(indiv.quality.dirn,sample.labels)

########################################################################
######################   Define global attributes ##################
######## Bulit ALL sample filter
quality.names<-c(quality.names,global.quality.names)  # lables in the range list quality.names<-c("ID","DP","QUAL","QD")  #"QUAL < 30.0 || AB > 0.75 && DP > 25 || QD < 5.0 || HRun > 5 || SB > -0.10"
quality.dirn<-c(quality.dirn,global.quality.dirn)

quality.labs<-c(quality.labs,global.quality.labs)  # what to call them in the filter
quality.cut<-c(quality,global.quality) ### remove the "." to keep the SNPs intact
names(quality.cut)<-quality.labs
names(quality.dirn)<-quality.labs
quality.type<-c(quality.type,global.quality.type)
names(quality.type)<-quality.names


avail<-names(quality.cut) %in% colnames(indels)
names(quality.cut)[!avail]
#### these two are in the same ORDER
## quality.cut
## quality.type
## quality.dirn
###################################

#########use RUN_FILTERING.r !!!!!!!!!
quality.thresh<-position.quality.filter(indels,quality.cut[avail],quality.type[avail],quality.dirn[avail])


## dim(quality.thresh)
## quality.thresh[1:5,]
rownames(quality.thresh)<-key.indels

