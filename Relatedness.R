#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
setwd(args[1])

# We are plotting Z0 (IBD =0) vs Z1 (IBD=1) and color code based on RT col which is the relationship (UN=unrelated and PO=parent-off spring)
pdf("relatedness.pdf")
relatedness = read.table("pihat_min0.2.genome", header=T)
par(pch=16, cex=1)
with(relatedness,plot(Z0,Z1, xlim=c(0,1), ylim=c(0,1), type="n"))
with(subset(relatedness,RT=="PO") , points(Z0,Z1,col="forestgreen"))
with(subset(relatedness,RT=="UN") , points(Z0,Z1,col="blue"))
legend(1,1, xjust=1, yjust=1, legend=levels(relatedness$RT), pch=16, col=c("blue","forestgreen"))
dev.off()

pdf("relatedness_ggplot.pdf")
ggplot(relatedness, aes(Z0,Z1, color=RT)) +
  geom_point(alpha=0.5, size=2) +
  labs(y="Z1", x="Z0", subtitle="Scatter plot - RT color") + xlim(0,1)+ylim(0,1)



pdf("zoom_relatedness.pdf")
relatedness_zoom = read.table("zoom_pihat.genome", header=T)
par(pch=16, cex=1)
with(relatedness_zoom,plot(Z0,Z1, xlim=c(0,0.02), ylim=c(0.98,1), type="n"))
with(subset(relatedness_zoom,RT=="PO") , points(Z0,Z1,col="forestgreen"))
with(subset(relatedness_zoom,RT=="UN") , points(Z0,Z1,col="blue"))
legend(0.02,1, xjust=1, yjust=1, legend=levels(relatedness$RT), pch=16, col=c("blue","forestgreen"))
dev.off()

pdf("hist_relatedness.pdf")
relatedness = read.table("pihat_min0.2.genome", header=T)
hist(relatedness[,10],main="Histogram relatedness", xlab= "Pihat")  
dev.off()

### use my histo with ggplot2

pdf("hist_relatedness_ggplot.pdf")
ggplot(relatedness, aes(x=PI_HAT)) + geom_histogram(bins=50) +
geom_vline(aes(xintercept=mean(PI_HAT)),color="blue", linetype="dashed", size=1) +
  scale_x_continuous(breaks=seq(0,1,0.05)) +
  theme(text = element_text(size=10),axis.text.x = element_text(angle=90, hjust=1)) 
dev.off()


#create a table of frequency for pi-hat distributions to see the vreak points
aa <- as.data.frame(table(relatedness$PI_HAT))
colnames(aa) <- c("Pi_hat", "freq")
write.table(aa,"15_Table_of_pihat_freq.txt", col.names=T, row.names=F, quote=F, sep="\t")


