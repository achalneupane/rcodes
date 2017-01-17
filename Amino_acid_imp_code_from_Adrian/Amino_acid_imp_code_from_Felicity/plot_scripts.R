
all.results = rbind(snp.results, aa.results,hla.results)
# 
write.table(snp.results,paste(file.stem,"snp.txt",sep="."),sep="\t",row.names=F)
write.table(hla.results,paste(file.stem,"hla.txt",sep="."),sep="\t",row.names=F)
write.table(aa.results,paste(file.stem,"aa.txt",sep="."),sep="\t",row.names=F)
write.table(all.results,paste(file.stem,"all.txt",sep="."),sep="\t",row.names=F)
# 
all.results = all.results[!is.na(all.results$Z),]

all.results = all.results[order(all.results$POS,decreasing=F),]
all.results$logp = -log10(all.results$PVAL)
#plot(all.results$POS,all.results$logp)

res1 = all.results

d = res1

d$ID = bim[ match(d$POS,bim$V4), 2]

d.alleles = d[d$TYPE %in% c("HLA2","HLA4"),]

genes.pos = list(
  HLAA=30019970,
  HLAC=31346171,
  HLAB=31431272,
  HLADRB1=32660042,
  HLADQA1=32716284,
  HLADQB1=32739039,
  HLADPA1=33145064,
  HLADPB1=33157346)

 t = d
 t$R2 = 0
 


for( i in 1:nrow(t) ) {
   t.pos = t[i,'POS']
   t.type = t[i,'TYPE']
 
   if( ! t.type %in% "SNP" ) next;
 
#   ld.b = ld[ld$BP_B == t.pos,]
#   ld.b = ld.b[ ! grepl("^AA",ld.b$SNP_B),]
#   ld.b = ld.b[ ! grepl("^HLA",ld.b$SNP_B),]
# 
#   r2 = max(ld.b$R2)

   t[i,'R2'] = 0.8
 }
# cond = 0

ymax=17
require(RColorBrewer)
blue.colors = colorRampPalette(brewer.pal(9,"Blues"))(100)
blue.colors = blue.colors[30:100]
size.pchs = seq(1,2,by=.01)

fig.height=5

svg(file=paste(fig.prefix,"_snps_aa.svg",sep=''),
     height=fig.height,width=14)

par(mar=c(5,5,2,2))

x.range = range(d$POS)
y.range = c(0,ymax)

toplot = d$TYPE == "SNP"

x = t$POS[toplot]
y = t$logp[toplot]
cols = 'blue'

plot(
  x,
  y,
  col=cols,
  pch=20,
  cex=1.3,
  xlim=x.range,ylim=y.range,
  ylab=expression(-log[10](italic(p))),
  xlab="chromosome 6 position (Mb)",
  bty='l',xaxt='n',cex.lab=1.7,cex.axis=1.4)

axis(1,at=seq(29e6,34e6,by=0.5e6),labels=seq(29,34,by=0.5),
     cex.axis=1.4)

toplot = d$TYPE == "AA"
points(t$POS[toplot],
       t$logp[toplot],
       col='red',pch=23,cex=2,lwd=2)



dev.off()

##########################################
## Now Plot Results for Imputed Alleles With labels##
##########################################


d2 = d[ d$TYPE %in% c("HLA2","HLA4"),]

alleles = c("HLA_A","HLA_C","HLA_B","HLA_DRB1","HLA_DQA1","HLA_DQB1","HLA_DPA1","HLA_DPB1")
alleles.p.s = c(1,4,7,10,13,16,19,22)


svg(file=paste(fig.prefix,"_hla.label.svg",sep=''),
    height=fig.height,width=14)

par(mar=c(5,5,2,2))

plot(0,
     type='n',
     axes=F,
     xlab="",
     ylab=expression(-log[10](italic(p))),
     xlim=c(1,24),
     ylim=c(0,ymax),cex.lab=1.7)

axis(2,cex.lab=1.7,cex.axis=1.4)

p.thr.lwd = 5
p.lwd = 5

for( i in 1:length(alleles) ) {
  allele = alleles[i]
  t = d2[grepl(allele,d2$NOTE) & d2$TYPE %in% "HLA2",]
  for( j in 1:nrow(t) ) {
    lwd = 1
    if( t[j,'logp'] > p.thr.lwd ) {
      lwd = p.lwd   
      txt.cx = 0.85
      allele.id = sapply(strsplit(as.character(t[j,'NOTE']),split='_'),function(x) return(x[3]))
      text(alleles.p.s[i] - 0.2, t[j,'logp'],labels=allele.id[1],cex=txt.cx);
    }
    segments(x0 = alleles.p.s[i] + 0.05,
             y0 = t[j,'logp'],
             x1 = alleles.p.s[i] + 0.95,
             lwd = lwd
             )
  }
  t = d2[grepl(allele,d2$NOTE) & d2$TYPE %in% "HLA4",]
  for( j in 1:nrow(t) ) {
    lwd = 1
    if( t[j,'logp'] > p.thr.lwd ) {
      lwd = p.lwd 
      txt.cx = 0.85
      allele.id = sapply(strsplit(as.character(t[j,'NOTE']),split='_'),function(x) return(x[3]))
      text(alleles.p.s[i] + 2.35, t[j,'logp'],labels=allele.id[1],cex=txt.cx);
    }
    segments(x0 = alleles.p.s[i] + 1.05,
             y0 = t[j,'logp'],
             x1 = alleles.p.s[i] + 1.95,
             lwd = lwd
             )
  }
}


txt.cx = 1.5
mtext(text="HLA-A",side=1,line=0.5,at=2,cex=txt.cx)
mtext(text="HLA-C",side=1,line=0.5,at=5,cex=txt.cx)
mtext(text="HLA-B",side=1,line=0.5,at=8,cex=txt.cx)
mtext(text="HLA-DRB1",side=1,line=0.5,at=11,cex=txt.cx)
mtext(text="HLA-DQA1",side=1,line=0.5,at=14,cex=txt.cx)
mtext(text="HLA-DQB1",side=1,line=0.5,at=17,cex=txt.cx)
mtext(text="HLA-DPA1",side=1,line=0.5,at=20,cex=txt.cx)
mtext(text="HLA-DPB1",side=1,line=0.5,at=23,cex=txt.cx)

dev.off()

##########################################
## Now Plot Results for Imputed Alleles Without labels##
##########################################


d2 = d[ d$TYPE %in% c("HLA2","HLA4"),]

alleles = c("HLA_A","HLA_C","HLA_B","HLA_DRB1","HLA_DQA1","HLA_DQB1","HLA_DPA1","HLA_DPB1")
alleles.p.s = c(1,4,7,10,13,16,19,22)


svg(file=paste(fig.prefix,"_hla.svg",sep=''),
    height=fig.height,width=14)

par(mar=c(5,5,2,2))

plot(0,
     type='n',
     axes=F,
     xlab="",
     ylab=expression(-log[10](italic(p))),
     xlim=c(1,24),
     ylim=c(0,ymax),cex.lab=1.7)

axis(2,cex.lab=1.7,cex.axis=1.4)

p.thr.lwd = 5
p.lwd = 5

for( i in 1:length(alleles) ) {
  allele = alleles[i]
  t = d2[grepl(allele,d2$NOTE) & d2$TYPE %in% "HLA2",]
  for( j in 1:nrow(t) ) {
    lwd = 1
    if( t[j,'logp'] > p.thr.lwd ) {
      lwd = p.lwd      
    }
    segments(x0 = alleles.p.s[i] + 0.05,
             y0 = t[j,'logp'],
             x1 = alleles.p.s[i] + 0.95,
             lwd = lwd
    )
  }
  t = d2[grepl(allele,d2$NOTE) & d2$TYPE %in% "HLA4",]
  for( j in 1:nrow(t) ) {
    lwd = 1
    if( t[j,'logp'] > p.thr.lwd ) {
      lwd = p.lwd 
    }
    segments(x0 = alleles.p.s[i] + 1.05,
             y0 = t[j,'logp'],
             x1 = alleles.p.s[i] + 1.95,
             lwd = lwd
    )
  }
}


txt.cx = 1.5
mtext(text="HLA-A",side=1,line=0.5,at=2,cex=txt.cx)
mtext(text="HLA-C",side=1,line=0.5,at=5,cex=txt.cx)
mtext(text="HLA-B",side=1,line=0.5,at=8,cex=txt.cx)
mtext(text="HLA-DRB1",side=1,line=0.5,at=11,cex=txt.cx)
mtext(text="HLA-DQA1",side=1,line=0.5,at=14,cex=txt.cx)
mtext(text="HLA-DQB1",side=1,line=0.5,at=17,cex=txt.cx)
mtext(text="HLA-DPA1",side=1,line=0.5,at=20,cex=txt.cx)
mtext(text="HLA-DPB1",side=1,line=0.5,at=23,cex=txt.cx)

dev.off()

