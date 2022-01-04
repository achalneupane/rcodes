#!/bin/bash
if [ $# -lt 3 ]; then
    echo "ERROR, must specify SRR#####, SampleID, and Project Name"
    echo "e.g. $0 SRR1391346 84467 dbGaP_phs000572_ADSP"
    exit 2
fi

SRR="${1}"
DNA="${SRR}"
SM="${2}"
PR="${3}"

# Rename to canonical
FULLSM="${SM}^${DNA}^${PR}"
mv -vi "${SRR}" "${FULLSM}" || { echo "ERROR, cannot change directory ${SRR} to ${FULLSM}, exiting"; exit 1; }
cd "${FULLSM}"
for BAMFILE in *.bam; do
  RGID=$(samtools view -H "${BAMFILE}" | grep "^@RG" | grep -o "ID:[^[:space:]]*" | sed 's/ID://g;s/\./^/g;s/_/^/g')
  samtools view -H "${BAMFILE}" | grep "^@RG" | sed "s@\bSM:\([^[:space:]]*\)\([[:space:]]\)@SM:${SM}\2DS:${FULLSM}\t@g;s/\t/\\\t/g" > ${FULLSM}.${RGID}.rgfile
  mv -vi "${BAMFILE}" "${FULLSM}.${RGID}.bam"
done
cd ..

exit 0
