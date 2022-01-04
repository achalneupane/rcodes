#!/usr/bin/Rscript --vanilla --slave --no-save

### Achal - 2020-01-13
#### GENERATE QC graphs from GATK-table file


#### LOAD R
library(tools)
require(data.table)
library(data.table)
library(ggplot2)

cat("\nreading passed by arguments\n")
args <- commandArgs(TRUE)
INFILE<-args[1]

print(paste0("Reading table for set ",INFILE))
TABLE<-read.table(INFILE, header=T)
#### OUTPUT SUMMARY TABLE
TABLEname<-file_path_sans_ext(INFILE)
print(paste0("generating summary for ",TABLEname))
out<-capture.output(summary(TABLE))
cat(paste0(TABLEname,".summary"),out,file=paste0(TABLEname,".summary"), sep="\n", append=TRUE)
for (PARAM in c("QD","FS","SOR","InbreedingCoeff","MQ","ReadPosRankSum","MQRankSum","DP")){
        print(paste0("generating PLOTS for ",PARAM))
        ggplot(TABLE, aes_string(x=PARAM, fill="FILTER")) + geom_density(alpha=0.2) + labs(title=paste0(TABLEname,PARAM))
        title<-paste0(TABLEname,"-",PARAM,".jpg")
        ggsave(title, plot=last_plot(), device = NULL, scale = 1, width = 16, height = 9, dpi = 300, limitsize = TRUE )
}
