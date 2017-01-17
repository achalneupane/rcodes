##########
if(exists("regionanno.DB")){
for(i in 1:length(regionanno.DB)){
# print(regionanno.DB[i])
the.region.DB<-regionanno.DB[i]
  if(grepl("+",the.region.DB)){the.region.DB<-gsub("+","",the.region.DB,fixed=TRUE)}

the.labels<-eval(as.name(paste(the.region.DB,"labels",sep=".")))
the.labels.wanted<-eval(as.name(paste(the.region.DB,"labels.wanted",sep=".")))

dim(the.labels)<-c(2,length(the.labels)/2)
the.labels<-t(the.labels)
colnames(the.labels)<-c("column","description")
# print(the.labels)

the.label.posns<-match(the.labels.wanted,the.labels[,"column"])
the.label.posns<-the.label.posns[!is.na(the.label.posns)]
the.label.posns<-as.integer(the.label.posns)-1 ## first is always the bin
# print(the.label.posns)

assign(paste(regionanno.DB[i],"labels",sep="."),value=the.labels[,"column"])
assign(paste(regionanno.DB[i],"labels.posns",sep="."),value=the.label.posns)
}
} # check exists


if(exists("regionanno.DB.score")){
for(i in 1:length(regionanno.DB.score)){
# print(regionanno.DB.score[i])
the.labels<-eval(as.name(paste(regionanno.DB.score[i],"labels",sep=".")))
the.labels.wanted<-eval(as.name(paste(regionanno.DB.score[i],"labels.wanted",sep=".")))

dim(the.labels)<-c(2,length(the.labels)/2)
the.labels<-t(the.labels)
colnames(the.labels)<-c("column","description")
# print(the.labels)

the.label.posns<-match(the.labels.wanted,the.labels[,"column"])
the.label.posns<-the.label.posns[!is.na(the.label.posns)]
the.label.posns<-as.integer(the.label.posns)-1 ## first is always the bin
# print(the.label.posns)
#if(exists( paste(regionanno.DB.score[i],"labels",sep=".") )){next}
assign(paste(regionanno.DB.score[i],"labels",sep="."),value=the.labels[,"column"])
assign(paste(regionanno.DB.score[i],"labels.posns",sep="."),value=the.label.posns)
}
} # check exists

############################# Note for geneanno MUST use the names "location" and "type"  as these are used above
############################# Geneanno have 2 files
if(exists("geneanno.DB")){
for(i in 1:length(geneanno.DB)){
# print(geneanno.DB[i])

the.labels<-c("location","type",core.ann)
the.labels.wanted<-c("location","type")
# print(the.labels)
# print(the.labels.wanted)

the.region.DB<-geneanno.DB[i]
if(grepl("+",the.region.DB)){the.region.DB<-gsub("+","",the.region.DB,fixed=TRUE)}
assign(paste(the.region.DB,"labels.variant_function",sep="."),value=the.labels)
assign(paste(the.region.DB,"labels.wanted.variant_function",sep="."),value=the.labels.wanted)

the.labels<-c("line","location","type",core.ann)
the.labels.wanted<-c("location","type")
# print(the.labels)
# print(the.labels.wanted)
#if(exists( paste(geneanno.DB[i],"labels",sep=".") )){next}
assign(paste(the.region.DB,"labels.exonic_variant_function",sep="."),value=the.labels)
assign(paste(the.region.DB,"labels.wanted.exonic_variant_function",sep="."),value=the.labels.wanted)
}
}

if(exists("filter.DB")){
for(i in 1:length(filter.DB)   ){
# print(filter.DB[i])
the.labels<-c("DB","maf",core.ann)
the.labels.wanted<-c("maf")
# print(the.labels)
# print(the.labels.wanted)
#if(exists( paste(filter.DB[i],"labels",sep=".") )){next} 
## assign(paste(filter.DB[i],"labels",sep="."),value=the.labels)
## assign(paste(filter.DB[i],"labels.wanted",sep="."),value=the.labels.wanted)

the.region.DB<-filter.DB[i]
if(grepl("+",the.region.DB)){the.region.DB<-gsub("+","",the.region.DB,fixed=TRUE)}
assign(paste(the.region.DB,"labels",sep="."),value=the.labels)
assign(paste(the.region.DB,"labels.wanted",sep="."),value=the.labels.wanted)
}
}

if(exists("generic.filter.DB")){
for(i in 1:length(generic.filter.DB)   ){
# print(generic.filter.DB[i])
the.labels<-c("DB","maf",core.ann)
the.labels.wanted<-c("maf")
# print(the.labels)
# print(the.labels.wanted)
#if(exists( paste(generic.filter.DB[i],"labels",sep=".") )){next}
the.region.DB<-generic.filter.DB[i]
if(grepl("+",the.region.DB)){the.region.DB<-gsub("+","",the.region.DB,fixed=TRUE)}
assign(paste(the.region.DB,"labels",sep="."),value=the.labels)
assign(paste(the.region.DB,"labels.wanted",sep="."),value=the.labels.wanted)
}
}

if(exists("function.filter.DB")){
for(i in 1:length(function.filter.DB)   ){
# print(function.filter.DB[i])
the.labels<-c("DB","score",core.ann)
the.labels.wanted<-c("score")
# print(the.labels)
# print(the.labels.wanted)
#if(exists( paste(function.filter.DB[i],"labels",sep=".") )){next} 
## assign(paste(function.filter.DB[i],"labels",sep="."),value=the.labels)
## assign(paste(function.filter.DB[i],"labels.wanted",sep="."),value=the.labels.wanted)

the.region.DB<-function.filter.DB[i]
if(grepl("+",the.region.DB)){the.region.DB<-gsub("+","",the.region.DB,fixed=TRUE)}
assign(paste(the.region.DB,"labels",sep="."),value=the.labels)
assign(paste(the.region.DB,"labels.wanted",sep="."),value=the.labels.wanted)
}
}
