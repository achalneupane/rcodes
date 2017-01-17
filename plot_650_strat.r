
## shellfish.py --file <.geno and .map prefix> --numpcs 10 --evecs <.evecs
## file> --snpload --out <out file prefix>

## could try admixture as well
## http://www.genetics.ucla.edu/software/admixture/
####################CASE 1#############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-FALSE
plot.dir<-"/media/Bioinform-D/Research/Matt Brown/Popluation Statification SNPs"
plot.pca.file<-"ALL_eth_650Y_WTCCC_C_F.pca.evec"
fam.file<-"ALL_eth_650Y_WTCCC_C_F.fam" # not required if used.shellfish<-FALSE

plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/Matt Brown/Popluation Statification SNPs"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"MND_shell2.evecs" # "MND_shell2.evecs"
fam.file<-"ch_mnd_650_f_ld_f.fam" ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"AOGC_HBM_650_shell_more.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"hbm_aogc_merge_f_c_ld_650_f.fam" ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"AOGC_HBM_shell_more.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"hbm_aogc_merge_f_c_ld_f_f.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"AOGC_HBM_650_TRUNC_shell.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"hbm_aogc_merge_f_c_ld_650_f_PCA.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"AOGC_HBM_650_TRUNC_S02_shell.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"hbm_aogc_merge_f_c_ld_650_f_PCA_1.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################


#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## /media/Bioinform-D/Research/GWAS extreme regression/AOGC_HBM/AOGC_sequencing_shell.evecs_pca_eigenvectors.txt
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"AOGC_sequencing_shell.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"hbm_aogc_merge_f_c_ld_SEQ_PCA.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## /media/Bioinform-D/Research/GWAS extreme regression/AOGC_HBM/AOGC_sequencing_shell.evecs_pca_eigenvectors.txt
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"AOGC_seq_2.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"hbm_aogc_seq_pca_f_c_ld_f.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################



#######################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"AOGC_sequencing_650_shell.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"hbm_aogc_merge_f_c_ld_SEQ_650_PCA.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-TRUE
plot.dir<-"/media/scratch2/cervical cancer final"
plot.pca.file<-"all_samples_shell.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"all_samples_qc_shellfish.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"as_controls_650_thin5_shell.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"as_controls_650_thin5.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################
#######################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"as_controls_650_thin6_shell.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"as_controls_650_thin6.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"scl_controls_ld_650_thin3_shell.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"scl_controls_ld_650_thin3.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################


#######################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"GBS_wtccc_pca3.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"GBS_wtccc_commom_ld_f_thin_f_clean_650.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################


#######################CASE 2##############################
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"GBS_wtccc_pca2.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"GBS_wtccc_commom_ld_f_thin_f_clean.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################
#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## /media/Bioinform-D/Research/GWAS extreme regression/AOGC_HBM/AOGC_sequencing_shell.evecs_pca_eigenvectors.txt
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"AOGC_seq_2.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"hbm_aogc_seq_pca_f_c_ld_f.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################
#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## first one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.GWAS.650Y.f.1.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.650Y.f.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.GWAS.650Y.f.2.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.650Y.f.1.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################
#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.GWAS.HBM.650Y.f.2.ld.c.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.HBM.650Y.f.1.ld.c.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.fam"  ### ** fam ** file used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################



#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.2.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1.fam"  ### ** fam ** file used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################



#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.fam"  ### ** fam ** file used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.3.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.fam"  ### ** fam ** file used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################


#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.nPB.1.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.nPB.fam"  ### ** fam ** file used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################
 


#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.nPB.2.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1.nCB.2.nPB2.fam"  ### ** fam ** file used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"sanity.1.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"sanity.fam"  ### ** fam ** file used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################







#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"test.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.650Y.f.1.ld.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################


#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.GWAS.650Y.f.2.ld.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.650Y.f.1.ld.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"good.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.650Y.f.1.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"good.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.650Y.f.1.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################





#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.2.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.2s.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1s.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################


#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.2s2.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1s2.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf1.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.nl.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.nl.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.5"
eig2<-"e.6"
#################################################

####################### WORKING FOR AOGC ##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.nl23.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.HBM.650Y.ld.f.b.1s2sf.nl23.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################
#exomeChip.FINAL.1.fam has a chr5 loaded split 
#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.FINAL.2.p.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.FINAL.2.p.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"


color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################
exomeChip.GWAS.HBM.650Y.strat.FINAL

#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.GWAS.HBM.650Y.strat.FINAL.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.GWAS.HBM.650Y.strat.FINAL.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"



color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#######################CASE 2##############################  THIS ONE USED FOR SEQUENCING - second attempt after found PCA has a associtaion
#work.dir<-"/media/scratch2/DORITHs_WORK/PCA STRANGE"
## ## final one for exomechip
used.shellfish<-TRUE
plot.dir<-"/media/Bioinform-D/Research/shellfish/shellfish"
plot.pca.file<-"exomeChip.FINAL.BEST.evecs" # "AOGC_HBM_shell_2.evecs"  "AOGC_HBM_650_shell_more.evecs"
fam.file<-"exomeChip.FINAL.BEST.fam"  ### fam files used with PCA analysis
plot.ann.file<-{} # leave as {} if set to default othersise tab delimited file with column names "Sample","status"
extra.ann.cols.wanted<-{} # list list like c("Call Rate","Weight") etc

ann.dir<-"/media/Bioinform-D/Research/GWAS_Figures"
ann.650Y.file<-"Sample_Detail_All_ethnicity.txt"



color.with<-"Ethnicity"

eig1<-"e.1"
eig2<-"e.2"
#################################################

#53800 wrks as well as 172000 see /media/Bioinform-D/Research/GWAS_Figures AOGC_HBM_650_shell.evecs vs  AOGC_HBM_650_shell_more.evecs
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
##############
#core.allowed.ann<-c("Sample","origin","Study","Product","Gender","Ethnicity")

#allowed.ann<-c(core.allowed.ann,extra.ann.cols.wanted)


################ get PCA to lot loaded with 650Y samples #######
if(exists("vec")){rm(vec)}
if(!used.shellfish){
vec<-read.delim(paste(plot.dir,plot.pca.file,sep="/"),header=T,sep="",fill=TRUE,stringsAsFactors=FALSE)  ## pca.evec file from eigenstrat
rownames(vec)<-unlist(lapply(strsplit(rownames(vec),split=":"),function(x) x[1]))   # fix colnum names
colnames(vec)<-c(paste("e.",1:10,sep=""),"status")   #,"origin","color","points")
}else{
  
vec<-read.delim(paste(plot.dir,plot.pca.file,sep="/"),header=F,sep="",fill=TRUE,stringsAsFactors=FALSE)  ## pca.evec file from eigenstrat
vec<-t(vec)
fam<-read.delim(paste(plot.dir,fam.file,sep="/"),header=F,sep="",fill=TRUE,stringsAsFactors=FALSE)

fam[fam[,6]==2,6]<-"Case"
fam[fam[,6]==1,6]<-"Control"

if(dim(fam)[1]!=dim(vec)[1]){
  print("Error fam and evec file had different number of samples that MUST BE the same if shellfish was used")}

rownames(vec)<-fam[,1]   # fix colnum names
vec<-cbind(vec,fam[,6])
colnames(vec)<-c(paste("e.",1:10,sep=""),"status")
pca.for.logistic<-cbind(rownames(vec),rownames(vec),vec[,c(paste("e.",1:10,sep=""))])
colnames(pca.for.logistic)<-c("FID","IID",c(paste("PCA",1:10,sep="")))
paste(plot.pca.file,"pca_eigenvectors.txt",sep="_")
getwd()
setwd(plot.dir)
write.table(pca.for.logistic,file=paste(plot.pca.file,"pca_eigenvectors.txt",sep="_"),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
}


##  plot(as.numeric(vec[,eig1]),as.numeric(vec[,eig2]))
## identify(as.numeric(vec[,eig1]),as.numeric(vec[,eig2]),labels=as.character(rownames(vec)))
## savePlot(paste("AOGC_by_self","ZOOM_Ethnicity2.png",sep="_"),type="png")

if(!is.null(plot.ann.file)){
  plot.ann<-read.delim(paste(ann.dir,plot.ann.file,sep="/"),header=T,sep="",fill=TRUE) ### annotation file pca.evec samplea not in 650Y
}else{
  if(!("status" %in% colnames(vec))){ vec[,"status"]<-"unknown"}
  plot.ann<-cbind(rownames(vec),as.character(vec[,"status"]))
  colnames(plot.ann)<-c("Sample","status")
}
###############

ann.650Y<-read.delim(paste(ann.dir,ann.650Y.file,sep="/"),header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ### annotation file for 650 Y
ann.650Y[1:5,]
dim(ann.650Y)
############# make combine sample annotation file this has the 650Y and the pca.evec samples now.
dim(plot.ann)
dim(ann.650Y)


sample.ann<-merge(ann.650Y,plot.ann,by.x="Sample",by.y="Sample",all=TRUE,sort=FALSE)
dim(sample.ann)

#sample.ann<-ann.650Y
                  
rownames(sample.ann)<-sample.ann[,"Sample"]


##### bothe the same size put in order
dim(sample.ann)
dim(vec)

posns<-match(rownames(sample.ann),rownames(vec))
missing<-is.na(posns)
if(sum(missing)>0){print("ERROR missing annotations")}
vec<-vec[posns[!missing],]
sample.ann<-sample.ann[!missing,]

dim(vec)


sample.ann[1:5,]
dim(sample.ann)


## cases<- fam
## controls<-fam
## fam<-rbind(cases,controls)
## posns<-match(rownames(sample.ann),fam[,1])
## missing<-is.na(posns)
## sample.ann[!missing,color.with]<-fam[posns[!missing],6]

## posns<-match(fam[,1],controls[,1])
## missing<-is.na(posns)
## sum(missing)
## fam[!missing,6]<-1
## table(fam[,6])
## write.table(fam,file=fam.file,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

#################################################33
sample.ann[is.na(sample.ann[,color.with]),color.with]<-as.character(sample.ann[is.na(sample.ann[,color.with]),"status"])
sample.ann[is.na(sample.ann[,color.with]),color.with]<-"Data"

sort(tapply(sample.ann[,"Ethnicity"],sample.ann[,"Ethnicity"],length))
unique(sample.ann[,"Ethnicity"])

## posns<-match(rownames(sample.ann),rownames(vec))
## missing<-is.na(posns)
## sum(missing) # problem is sample id not found  as vec samples are a sunset of all known samples
## vec<-vec[posns,]
#vec<-vec[!missing,]

# posns<-match(rownames(sample.ann),fam[,1)
## missing<-is.na(posns)
## sum(missing) # problem is sample id not found  as vec samples are a sunset of all known samples
## vec<-vec[posns,]
#vec<-vec[!missing,]

tapply(rownames(sample.ann),sample.ann[,color.with],length) # coulsnts
unique(ann.650Y[,color.with])


############################## Get orginal colors for 650Y ###################
classes<-levels(as.factor(ann.650Y[,color.with]))  #order ok with the above
color.set<-rainbow(length(classes))          #get auto colors
color.set[10]<-"yellow4"
color.set[20]<-"forestgreen"
color.set[17]<-"lightblue"

pch.set<-c(0:25,33:45,47:58) # 650Y ethnicity
pch.set[21]<-pch.set[5]
pch.set[35]<-pch.set[2] #palastian
pch.set[36]<-pch.set[5]
pch.set[28]<-pch.set[4]
pch.set[31]<-pch.set[3]
pch.set[34]<-pch.set[6]
pch.set[32]<-pch.set[16]
pch.set[33]<-pch.set[1]
pch.set[30]<-pch.set[9]
pch.set[38:51]<-pch.set[1:14]

classes.full<-levels(as.factor(sample.ann[,color.with]))
classes.extra<-classes.full[!(classes.full %in% classes)]

classes<-c(classes.extra,classes)

gray.steps<-seq(from=0,to=0.75,by=(1-0.25)/length(classes.extra))
color.extra<-gray(gray.steps[1:length(classes.extra)])
pch.extra<-rep(20,times=length(classes.extra))

color.set<-c(color.extra,color.set)
pch.set<-c(pch.extra,pch.set) # 650Y ethnicity

names(color.set)<-classes
names(pch.set)<-classes

color.array<-color.set[sample.ann[,color.with]]
pch.array<-pch.set[sample.ann[,color.with]]
length(color.array)
length(pch.array)

range.eig1<-range(as.numeric(vec[,eig1]))
range.eig2<-range(as.numeric(vec[,eig2]))

include<-c(1:dim(vec)[1])
#include<-match(rownames(class.results),rownames(vec))
#include<-include[!is.na(include)]
 par(mfrow=c(1,1),font=2,font.lab=2,font.axis=2,mgp=c(3.1,1,0),mar=c(5,5,4,2)+0.1)
######## all data : plot(as.numeric(vec[,eig1]),as.numeric(vec[,eig2]))
bp<-plot(as.numeric(vec[include,eig1]),as.numeric(vec[include,eig2]),pch=pch.array, cex=1.35,lwd=1.0,col = color.array[include], xlim=range.eig1,ylim=range.eig2,main="",xlab=eig1, ylab=eig2,cex.lab=2.0,cex.axis=2.0)

legend(range.eig1[1],range.eig2[2],classes,col=color.set,pch=pch.set,cex=0.85,bty="n",xjust=0.3)
setwd(plot.dir)

savePlot(paste(plot.pca.file,"Ethnicity.FINAL.BEST.tiff",sep="_"),,type="tiff")
savePlot(paste(plot.pca.file,"Ethnicity.FINAL.BEST.png",sep="_"),type="png")
# savePlot(filename=paste("ALL ethnicity 650Y","png",sep="."),type="png")
test<-identify(as.numeric(vec[include,eig1]),as.numeric(vec[include,eig2]),labels=as.character(rownames(vec)[include]))
as.character(rownames(vec)[include])[test]
vec[test,]
identify(as.numeric(vec[include,eig1]),as.numeric(vec[include,eig2]),labels=as.character(rownames(vec)[include]),col="Black")

abline(v=0.055)
savePlot(paste(plot.pca.file,"Ethnicity_labels.FINAL.BEST.tiff",sep="_"),,type="tiff")
savePlot(paste(plot.pca.file,"Ethnicity_labels.FINAL.BEST.png",sep="_"),type="png")
abline(h=0.0)
################################################## zoom in options



####################################################
####################################################
####################################################
####################################################
####################################################

################## KEEP DATA SAMPLES samples  REMOVE UNWAHTED AND ANNOTAION SAMPLES
## dim(vec)
## vec[test,]
## sample.ann[1:5,]

## abline(h=-0.1)
## abline(v=0.05)

## dim(vec)
## dim(sample.ann)

## tapply(vec[,"status"],vec[,"status"],length)
## sort(tapply(sample.ann[,"Ethnicity"],sample.ann[,"Ethnicity"],length))


slice1<- as.numeric(as.character(vec[,"e.1"]))<= 0.0035 ## choose WANTED
slice2<- as.numeric(as.character(vec[,"e.2"]))> -0.2  ## choose WANTED
current.data<-sample.ann[,"Ethnicity"]=="Data" | sample.ann[,"Ethnicity"]=="Case" | sample.ann[,"Ethnicity"]=="Control" | sample.ann[,"Ethnicity"]=="-9"

dim(sample.ann)
sum(current.data)
sum(slice1)
sum(slice2)

wanted<-current.data & slice1 & slice2
sum(current.data)
sum(wanted)
sum(current.data)-sum(wanted) # 47 removed in phase 1

sum(!wanted)
table(sample.ann[wanted,"Ethnicity"]) #    558     430  for scl

getwd()
paste(plot.pca.file,"keep.txt",sep="_")
write.table(cbind(rownames(vec)[wanted],rownames(vec)[wanted]),file=paste(plot.pca.file,"keep_6SD.txt",sep="_"),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)


slice1<- as.numeric(as.character(vec[,"e.1"]))<= 0.01 ## choose WANTED
slice2<- as.numeric(as.character(vec[,"e.2"]))> -0.04  ## choose WANTED
current.data.w.Spikes<- slice1 & slice2
table(sample.ann[current.data.w.Spikes,"Ethnicity"])
paste(plot.pca.file,"wSpikes.01.keep.txt",sep="_")
write.table(cbind(rownames(vec)[current.data.w.Spikes],rownames(vec)[current.data.w.Spikes]),file=paste(plot.pca.file,"wSpikes.keep.01.txt",sep="_"),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)


           -9        Adygei        Basque       Bedouin Biaka_Pygmies 
         7592            16            24            11             1 
         Case       Control         Druze        French       Italian 
          781           444            46            29            13 
     Orcadian   Palestinian       Russian     Sardinian        Tuscan 
           16            13            25            28             8


           -9        Adygei       Balochi        Basque       Bedouin 
         7604            17            18            24            47 
Biaka_Pygmies        Brahui          Case       Control         Druze 
            1            20           781           445            47 
       French       Italian        Kalash       Makrani      Mozabite 
           29            13            21            14             4 
     Orcadian   Palestinian        Pathan       Russian     Sardinian 
           16            49             8            25            28 
       Sindhi        Tuscan 
            2             8 
> 
####################################################
####################################################
####################################################
####################################################
abline(v=0.001)

include<- as.numeric(as.character(vec[,"e.1"]))<= 0.03
#hist(as.numeric(vec[include,eig1]))
shapiro.test(sample(as.numeric(vec[include,eig1]),3500)) #  shapiro.test(sample(rnorm(8000, mean = 5, sd = 3),3500))
mean(as.numeric(vec[include,eig1]))
sd(as.numeric(vec[include,eig1]))
mean(as.numeric(vec[include,eig1]))+sd(as.numeric(vec[include,eig1]))*6
points(rnorm(mean=mean(as.numeric(vec[include,eig1])),sd=sd(as.numeric(vec[include,eig1]))

-0.002815501
0.0009711256*6
ship

# 426  and 4724
4724+426


vec[1:5,]
vec[include,][2962,]
test<-"hi_m_1788"
test<-"S12-F05-P01"
test<-"AOGC-02-0066"
test<-"S02-F01-P01"

"AOGC-08-0258" # not genotyped
outlyers<- 2  21 134 233 305 421 432 440 490 508 526 554 681 683 690 762 814 848 852 858 882 898 924
test<-"AOGC-14-3994"
test<-"AOGC-02-0400"
test<-"AOGC-14-2686"
test<-"AOGC-14-1118"
test<-"AOGC-02-0210"
test<-"AOGC-14-3359"
test<-"AOGC-02-0162"
points(vec[test,eig1],vec[test,eig2],col="cyan") # points(pca[test,2],pca[test,3],col="cyan")
abline(v=0.3)

strat<-rownames(vec) %in% strat.sequenced
strat<-rownames(vec) %in% ped[missing,"PATIENT"]

sum(strat) rownames(vec)[strat]
points(vec[strat,eig1],vec[strat,eig2],col="cyan",cex=4)
vec[strat,]

range.eig1
range.eig2

range.eig1<-c(-0.0062,0.005)
range.eig2<-c(-0.01,0.01)
  
bp<-plot(as.numeric(vec[include,eig1]),as.numeric(vec[include,eig2]),pch=pch.array, cex=1.35,lwd=1.0,col = color.array[include], xlim=range.eig1,ylim=range.eig2,main="",xlab="Eigenvector 1", ylab="Eigenvector 2",cex.lab=2.0,cex.axis=2.0)
legend(range.eig1[1],range.eig2[2],classes,col=color.set,pch=pch.set,cex=0.85,bty="n",xjust=0.3)


abline(v=0.0035)
abline(h=-0.002)

savePlot(paste(plot.pca.file,"ZOOM_Ethnicity.BEST.6SDtiff",sep="_"),,type="tiff")
savePlot(paste(plot.pca.file,"ZOOM_Ethnicity.BEST.6SD.png",sep="_"),type="png")


