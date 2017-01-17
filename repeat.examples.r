####################################
################Poly morphic SITE TESTS

genome<-"BSgenome.Hsapiens.UCSC.hg19"
require(genome,character.only=TRUE)
require("Biostrings")
require("IRanges")

##########
load("test.RData")
source("repeat.finder.subroutines.r")

test[1:5,]
################ 


a.indel<-test
key<-rownames(a.indel)

REF.length<-nchar(as.character(a.indel[,"REF"]))
ALT.length<-nchar(as.character(a.indel[,"ALT"]))

large.indel<-REF.length>1 | ALT.length>1

### check to see IF REF or ALT is a homopolymer run length > 5 or dinucleitide repeat - > 6pb 
are.repeats<-identify.repeats(a.indel,di.run.max=3,homo.run.max=5) 

length(large.indel) 
length(are.repeats)
sum(are.repeats)
rownames(a.indel)[are.repeats][1:20] # example
rownames(a.indel)[!are.repeats][1:5] #example

#################### CHECK to see if they ARE in a repeat region in a genomic context:

chk.in.repeat<-large.indel & !are.repeats
sum(chk.in.repeat)
## check 6 bases forward:
are.sub.repeat<-indentify.IN.repeat(a.indel[chk.in.repeat,],looking="forward",bases.about=6,di.run.max=3,homo.run.max=5,genome="BSgenome.Hsapiens.UCSC.hg19")
remove.repeats<-key[chk.in.repeat][are.sub.repeat]
are.in.repeats.forward<- key %in% remove.repeats

remove.repeats[1:20] ### test ones have repats in forward direction
sum(are.in.repeats.forward)
## [1] 6988

###################### in repeats looking back

sum(chk.in.repeat)
chk.in.repeat<-large.indel & !are.repeats & !are.in.repeats.forward

are.sub.repeat<-indentify.IN.repeat(a.indel[chk.in.repeat,],looking="back",bases.about=6,di.run.max=3,homo.run.max=5,genome="BSgenome.Hsapiens.UCSC.hg19")
remove.repeats<-key[chk.in.repeat][are.sub.repeat]
are.in.repeats.back<- key %in% remove.repeats

remove.repeats[1:20]
sum(are.in.repeats.back) ### test ones have repeats in backward direction
## [1] 3224

are.in.repeats<- are.in.repeats.back | are.in.repeats.forward

length(are.in.repeats)
sum(are.in.repeats)

