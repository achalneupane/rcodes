#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Error, must specify SRR number"
    exit 2
fi

SRR="${1}"
if [[ ! ${SRR} =~ ^SRR ]]; then
    echo "Error, SRR number must begin with 'SRR'"
    exit 1
fi

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

exitcode=0
IFS=$'\n'
SRAINFO=($(/usr/local/genome/bin/sra-stat --quick ${SRR}))
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

echo "All spots from ${SRR} are represented in associated .fq or .fq.gz files"
exit ${exitcode}
