


options(width=2000)
# IBD check
related<-read.table("/media/ga-apps/UQCCG-Analysis/AOGC_exome_chip/IBD/related_samples.txt",header=T,fill=TRUE,stringsAsFactors=FALSE)
colnames(related)<-c("ID1","ID2","P(IBD=0)","P(IBD=1)","P(IBD=2)","PI_HAT","Call rate ID1","Call rate ID2")

call.rates.form.imiss<-FALSE ## FALSE is special case of data from adrians pipline

#### here related conatins the Call rates use imiss otherwise
if(call.rates.form.imiss){
the.call.rates<-read.table("/media/scratch2/AOGC-NGS/ExomeChip/UQDI.imiss",header=T,fill=TRUE,stringsAsFactors=FALSE)
call.rates<-cbind(the.call.rates[,"IID"], 1-as.numeric(the.call.rates[,"F_MISS"]))
colnames(call.rates)<-c("IID","call.rate")
}else{

########################### get "call.rates IID vs "call.rate"
rel1<-related[,c(1,7)]
rel2<-related[,c(2,8)]

colnames(rel1)<-c("IID","call.rate")
colnames(rel2)<-c("IID","call.rate")
                  
call.rates<-rbind(rel1,rel2)
the.dups<-duplicated(call.rates[,"IID"])
call.rates<-call.rates[!the.dups,]
}
##############################################

dim(related)
length(unique(related[,1]))
length(unique(related[,2]))

related[1:5,]
call.rates[1:5,]
                     
all<-c(related[,"ID1"],related[,"ID2"])
#all<-related[,1]
counts<-sort(tapply(all,all,length),decreasing=TRUE)
multi<-counts[counts>1]
multi
##### trim multi based on all rate

to.remove.samples<-{}
related.ori<-related


################# recursively remove the sample that appear more than once with the lowerest call rate till all only singles

while(!is.null(dim(multi))){


multi  
posns<-match(names(multi),call.rates[,"IID"])
order.by<-order(multi,(1-as.numeric(call.rates[posns,"call.rate"])),decreasing=TRUE)
multi<-multi[order.by]
to.remove<-names(multi)[1]

to.remove.samples<-c(to.remove.samples,to.remove)


#### find out which column ha the most entries
col1.counts<-sum(related[,"ID1"]==to.remove)
col2.counts<-sum(related[,"ID2"]==to.remove)

if(col1.counts >=col2.counts){
related<-related[!related[,"ID1"]==to.remove,]   
#related[related[,"ID1"]==to.remove,"ID1"]<-paste("dummy",related[related[,"ID1"]==to.remove,"ID1"],"col1",1:sum((related[,"ID1"]==to.remove)),sep="_")
}else{
related<-related[!related[,"ID2"]==to.remove,] 
}

all<-c(related[,"ID1"],related[,"ID2"])
counts<-sort(tapply(all,all,length),decreasing=TRUE)
multi<-counts[counts>1]
print(multi)
print(to.remove.samples)
                       }


######### now remove half of the singles chosing ones with the worst call rate. 
######### related now only contains unique samples
related[1:8,]

call.rates[1:8,]

posns1<-match(related[,"ID1"],call.rates[,"IID"])
posns2<-match(related[,"ID2"],call.rates[,"IID"])


#to.swap<-as.numeric(related[,"Call rate ID1"]) > as.numeric(related[,"Call rate ID2"])
to.swap<-as.numeric(call.rates[posns1,"call.rate"]) > as.numeric(call.rates[posns2,"call.rate"]) 

#cbind(related,call.rates[posns1,"call.rate"],as.numeric(call.rates[posns2,"call.rate"]),to.swap)[1:10,]

related[to.swap,"ID1"]<-related[to.swap,"ID1"]

to.remove.samples<-c(to.remove.samples,related[,"ID1"])
length(to.remove.samples)

write.table(cbind(to.remove.samples,to.remove.samples),file="related.to.remove.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
length(to.remove.samples)
################################################# Check frequency for calling algorithms


getwd()




############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
