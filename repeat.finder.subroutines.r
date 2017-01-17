indentify.IN.repeat<-function(indel,looking="forward",bases.about=6,di.run.max=3,homo.run.max=5,genome="BSgenome.Hsapiens.UCSC.hg19"){
  # requires "chr","start","end","REF","ALT","TYPE" columns
  
require(genome,character.only=TRUE)
require("Biostrings")
require("IRanges")
#Hsapiens

## test<-indel[chk.in.repeat,1:10]
## test["chr5:172385267:172385267:-:AA:indel",]
## indel<-a.indel[chk.in.repeat,]
## print(bases.about)


if( (bases.about< (homo.run.max+1)) | (bases.about< (di.run.max*2)) ){bases.about<-max(c((homo.run.max+1),(di.run.max*2)))} # look to look past the length of the repeart run
## chk.in.repeat<-large.indel & !are.repeats
## sum(chk.in.repeat)


if(is.null(rownames(indel))){core.ann<-c("chr","start","end","REF","ALT","TYPE") ;key.small<-build.key(indel,core.ann,add.chr.label=FALSE)}else{key.small<-rownames(indel)}
the.chrom<-indel[,"chr"]

### look forward
if(looking=="forward"){
starts<-as.numeric(indel[,"end"]) 
ends<-as.numeric(indel[,"end"]) +   bases.about
}else{ ## looking back
starts<-as.numeric(a.indel[chk.in.repeat,"start"]) -  bases.about 
ends<-as.numeric(a.indel[chk.in.repeat,"start"])   
}
##### stay within the existing region 
## starts<-pmax.int(starts,peaks[,"start"])
## ends<-pmin.int(ends,peaks[,"end"])
############################ found if I did not do this I got a single extra region  15496:  16 89084207 89084461  Krtap16-7


large.indels.test<-getSeq(Hsapiens, the.chrom, starts, ends) # large.indels.test<-all.genomic
## length(all.genomic)
## system.time(x<- XStringViews(all.genomic, "DNAString"))
## x.labels<-paste(the.chrom,starts,ends,sep=":")
## names(x)<-x.labels



di<-dinucleotideFrequency(large.indels.test,step=1)
widths<-width(large.indels.test)
#homo.run<-pmax(di[,colnames(di) %in% c("AA","TT","GG","CC")])
homo.run<-pmax(di[,"AA"],di[,"TT"],di[,"GG"],di[,"CC"])
homo.run[homo.run>0]<-homo.run[homo.run>0]+1  
#di.run<-pmax(di[,colnames(di) %in% c("AC","AG","AT","CG","CT","GT")])
di.run<-pmax(di[,"AC"],di[,"AG"],di[,"AT"],di[,"CG"],di[,"CT"],di[,"GT"])
di.run<-di.run*2 # number of bases in repeats

## di.run[1:50]
## homo.run[1:50]
## widths[1:50]

repeats<-(di.run==widths & (di.run/2)>=di.run.max) | (homo.run>=(widths-1) &  homo.run>=homo.run.max)


## repeats[1:5]
## grep(TRUE,repeats)[1:5]
## large.indels.test[grep(TRUE,repeats)[1:5]]
## length(repeats)
## sum(repeats)
## key.small[repeats][1:5]
## homo.run[repeats][1:5]
## di.run[repeats][1:5]

return(repeats)

## remove.repeats<-key.small[repeats]
## length(remove.repeats)
## remove.repeats[1:20]

## are.in.repeats.forward<- key %in% remove.repeats

## are.in.repeats.forward[1:20]
## sum(are.in.repeats.forward)

} # indentify.IN.repeat

###############################################################

identify.repeats<-function(indel,di.run.max=3,homo.run.max=5){
  # requires "chr","start","end","REF","ALT","TYPE" columns 
require(Biostrings)
## poly.morphic<-grepl(":\\d*$",key,perl=T)
## key[poly.morphic][1:50]
if(is.null(rownames(indel))){core.ann<-c("chr","start","end","REF","ALT","TYPE") ;key<-build.key(indel,core.ann,add.chr.label=FALSE)}else{key<-rownames(indel)}
REF.length<-nchar(as.character(indel[,"REF"]))
ALT.length<-nchar(as.character(indel[,"ALT"]))

large.indel<-REF.length>1 | ALT.length>1

#sum(large.indel)

the.large.indels<-indel[,"REF"][large.indel]
use.alt<-ALT.length[large.indel]>REF.length[large.indel]
#the.large.indels[1:5]
the.large.indels[use.alt]<-indel[,"ALT"][large.indel][use.alt]



large.indels.test<-DNAStringSet(the.large.indels)
#large.indels.test[1:50]
di<-dinucleotideFrequency(large.indels.test,step=1)
widths<-width(large.indels.test)
#homo.run<-pmax(di[,colnames(di) %in% c("AA","TT","GG","CC")])
homo.run<-pmax(di[,"AA"],di[,"TT"],di[,"GG"],di[,"CC"])
homo.run[homo.run>0]<-homo.run[homo.run>0]+1  
#di.run<-pmax(di[,colnames(di) %in% c("AC","AG","AT","CG","CT","GT")])
di.run<-pmax(di[,"AC"],di[,"AG"],di[,"AT"],di[,"CG"],di[,"CT"],di[,"GT"])
di.run<-di.run*2 # number of bases in repeats

## di.run[1:50]
## homo.run[1:50]
## widths[1:50]

repeats<-(di.run==widths & (di.run/2)>=di.run.max) | (homo.run>=(widths-1) &  homo.run>=homo.run.max)

## length(repeats)
## sum(repeats)
## the.large.indels[repeats][1:5]
## homo.run[repeats][1:5]
## di.run[repeats][1:5]

remove.repeats<-key[large.indel][repeats]
## length(remove.repeats)
## remove.repeats[1:50]

are.repeats<- key %in% remove.repeats
## sum(are.repeats)

## length(large.indel) 
## length(are.repeats)
dim(indel)

return(are.repeats)
}

