

well.type<-96
row.type<-8
col.type<-12

row.index<-c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")

plate.ids<-unique(c("plate1","plate2")) # for example
plate.ids<-sort(plate.ids)
for(i in 1:length(plate.ids)){
 print(plate.ids[i])
  print(tapply( ann[ann[,"plate"]==plate.ids[i],"row"],ann[ann[,"plate"]==plate.ids[i],"row"],length)) }

############### test plate nummbering
## plate.ids<-unique(ann[,"plate"])
plate.ids.t<-rep(plate.ids,well.type)
dim(plate.ids.t)<-c(length(plate.ids),well.type)
plate.ids.t<-t(plate.ids.t)
plate.ids.t<-as.vector(plate.ids.t)

############### 
row.ids<-row.index[1:row.type]
row.ids.t<-rep(row.ids,col.type)
dim(row.ids.t)<-c(length(row.ids),col.type)
row.ids.t<-t(row.ids.t)
row.ids.t<-as.vector(row.ids.t)
row.ids.t<-rep(row.ids.t,length(plate.ids))

###################

col.ids<-1:col.type
col.ids.t<-rep(col.ids,row.type)
col.ids.t<-rep(col.ids.t,length(plate.ids))

complete.plates<-cbind(plate.ids.t,row.ids.t,col.ids.t)

complete.plates[1:10,] # what a complete 2 plate would look like in vector form


############################################################
# make up some fake data
ann<-complete.plates
ann<-cbind(ann,NA)
colnames(ann)<-c("plate","row","col","Symbol")
ann[,"Symbol"]<-1:dim(ann)[1]

ann[1:10,]
#### make an error
swap<-ann[100,]
ann[100,]<-ann[10,]
ann[10,]<-swap


## ann[10,]
## ann[10,]<-c("plate3","P","45","shit") ### insert an error
#####################################


order.wanted<-paste(plate.ids.t,row.ids.t,col.ids.t,sep="::")
order.have<-paste(ann[,"plate"],ann[,"row"],ann[,"col"],sep="::")



tapply(ann[,"plate"],ann[,"plate"],length)
tapply(ann[,"row"],ann[,"row"],length)
tapply(ann[,"col"],ann[,"col"],length)
tapply(ann[,"Symbol"],ann[,"Symbol"],length)


#####rows ot of order:
sum(row.ids.t!=ann[,"row"]) # must be zero
test<-row.ids.t!=ann[,"row"]
ann[test,1:4]


####columns out of order
sum(col.ids.t!=ann[,"col"]) # must be zero
test<-col.ids.t!=ann[,"col"]
ann[test,1:4]

### test to see in all orders ok
sum(plate.ids.t!=ann[,"plate"])  # must be zero


######## fix the incorrecting ordering:
order.have<-paste(ann[,"plate"],ann[,"row"],ann[,"col"],sep="::")
posns<-match(order.wanted,order.have)
sum(is.na(posns)) # must be zero

ann<-ann[posns,]
################## now go back to above and run checks again

############################################################
