#!/bin/bash

LOOKUP="/30/dbGaP/6109/sra/lookup.csv"

if [ $# -lt 1 ]; then
    echo "Error, must specify SRR number"
    exit 2
fi

SRR="${1}"
# if [ ! -r "${SRR}.sra" ]; then
#     echo "Error, file ${SRR}.sra does not exist or cannot be read"
#     exit 1
# fi

IFS=$'\n'
RGLINES=($(sam-dump ${SRR} | sed -n '/^[^@]/!p;//q' | grep ^@RG))
: '
args=(tee)
for RGLINE in ${RGLINES[@]}; do
  unset IFS
  RG=(${RGLINE})
#  args+=(\>\(grep -A3 --no-group-separator \"\\.${RG[1]#ID:}/\" \> ${SRR}.${RG[1]#ID:}.fastq-dump.split.defline.z.tee.fq\))
  args+=(\>\(grep -A3 --no-group-separator \"\\.${RG[1]#ID:}/\" \| gzip \> ${SRR}.${RG[1]#ID:}.fastq-dump.split.defline.z.tee.fq.gz\))
done
args+=(\>/dev/null)

echo "Splitting ${SRR}.sra into ${#RGLINES[@]} ReadGroups"
/usr/local/genome/bin/fastq-dump --split-3 --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' -Z "${SRR}.sra" | eval ${args[@]}

if [ $? -ne 0 ]; then
    echo "Error running fastq-dump, exiting."
    exit 1
fi

# Validate the .fq.gz that are produced
if [ $(ls ${SRR}*fq 2>/dev/null | wc -l) -eq 0 ]; then
    if [ $(ls ${SRR}*fq.gz 2>/dev/null | wc -l) -eq 0 ]; then
	echo "Error, cannot find any .fq or .fq.gz files for ${SRR}"
	exit 1
    else
	MODE=gz
    fi
else
    MODE=fq
fi

echo "Validating .fq.gz created from ${SRR}"
exitcode=0
'
# IFS=$'\n'
: '
SRAINFO=($(/usr/local/genome/bin/sra-stat --quick "${SRR}.sra"))
for line in ${SRAINFO[@]}; do
    IFS="|"
    split1=(${line})
    RG=${split1[1]}
    IFS=":"
    split2=(${split1[2]})
    READS=${split2[0]}
    ((READS*=8))
    unset IFS
    echo -n "Checking ${SRR} ReadGroup ${RG}, expect ${READS} lines..."
    if [ ${MODE} = "gz" ]; then
	LINES=$(zcat ${SRR}.${RG}.*.fq.gz | wc -l)
    elif [ ${MODE} = "fq" ]; then
	LINES=$(wc -l ${SRR}.${RG}.*.fq)
    fi
    echo "found ${LINES}"
    if [ ${READS} -ne ${LINES} ]; then
	((exitcode+=1))
    fi
done

if [ ${exitcode} -eq 0 ]; then
  echo "All spots from ${SRR} are represented in associated .fq or .fq.gz files"
  if [ -f "${SRR}.sra" ]; then
    rm -fv "${SRR}.sra"
    rm -fv "${SRR}.sra.vdbcache.cache" 2>/dev/null
    # This file gets created by sra-stat, even though it already exists in subdir
    rm -fv "/30/dbGaP/6109/sra/${SRR}.sra.vdbcache.cache" 2>/dev/null
  fi
'
  SM=$(grep -m1 "^${SRR}," "${LOOKUP}" | cut -d, -f2)
  PR=$(grep -m1 "^${SRR}," "${LOOKUP}" | cut -d, -f3)
  DNA="${SRR}"
  FULLSM="${SM}^${DNA}^${PR}"
  mkdir -p "ready/${FULLSM}"
  IFS=$'\n'
  for RGLINE in ${RGLINES[@]}; do
    OLD_RGID=$(echo ${RGLINE} | grep -o "ID:[^[:space:]]*" | sed 's/ID://g')
    NEW_RGID=$(echo ${OLD_RGID} | sed 's/\./^/g;s/_/^/g')
    mv -vi "${SRR}.${OLD_RGID}.fastq-dump.split.defline.z.tee.fq.gz" "ready/${FULLSM}/${FULLSM}.${NEW_RGID}.fq.gz"
    echo ${RGLINE} | sed "s@\bSM:\([^[:space:]]*\)\([[:space:]]\)@SM:${SM}\2@g;s/\t/\\\t/g" > "ready/${FULLSM}/${FULLSM}.${NEW_RGID}.rgfile"
  done
# fi
unset IFS

exit $?
