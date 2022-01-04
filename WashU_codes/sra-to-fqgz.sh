#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Error, must specify SRR number"
    exit 2
fi

SRR="${1}"
if [ ! -r "${SRR}.sra" ]; then
    echo "Error, file ${SRR}.sra does not exist or cannot be read"
    exit 1
fi

IFS=$'\n'
RGLINES=($(sam-dump ${SRR}.sra | sed -n '/^[^@]/!p;//q' | grep ^@RG))
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

exit $?
