#!/bin/bash
declare -A projects=(["Broad_WGS"]="174" ["DIAN"]="12" ["Genentech_WGS"]="47" ["MAPT_A152T"]="21" ["Macrogen_WGS"]="20" ["dbGaP_phs000572_ADSP"]="117" ["Genentech_WES"]="92" ["LOAD_WES"]="33" ["Otogenetics_WES"]="834" ["TGI_WES"]="298" ["MGI_201605"]="424" )
{ echo "Project Completed Total"; echo "------- -------- -----"; for i in \
    Broad_WGS \
    DIAN \
    Genentech_WGS \
    MAPT_A152T \
    Macrogen_WGS \
    dbGaP_phs000572_ADSP \
    Genentech_WES \
    LOAD_WES \
    Otogenetics_WES \
    TGI_WES \
    MGI_201605; do\
    echo -n "${i} "; ls /80/AD/AD_Seq_Data/02.-BWA_GATK/GRCh37/*${i} -d 2>/dev/null | wc -l | tr -d "\n"; echo " ${projects[${i}]}"; done; } | column -t

DATE=$(date "+%Y%m%d")
find /80/AD/AD_Seq_Data/02.-BWA_GATK/GRCh37/ -name "*.g.vcf.gz" | xargs ls -ltr --full-time > "gvcfgz_dates_${DATE}.csv"
spaces_to_one_tab.sh "gvcfgz_dates_${DATE}.csv" | sort -k6,6 | cut -f6 > "gvcfgz_justdates_${DATE}.csv"
wc -l "gvcfgz_justdates_${DATE}.csv"

exit 0

