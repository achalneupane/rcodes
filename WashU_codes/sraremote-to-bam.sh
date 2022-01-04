#!/bin/bash
if [ $# -lt 1 ]; then
  echo "Error, must provide SRR run number"
  exit 2
fi

if [[ ! ${1} =~ ^SRR ]] ; then
  echo "Error, SRR run reference must begin with SRR"
  exit 1
fi

JAVA="/usr/java/jre1.8.0_91/bin/java"
JAVAOPTS="-Xms4g -Xmx8g -XX:ParallelGCThreads=1 -Djava.io.tmpdir=/data/tmp"
PICARD="/usr/local/genome/picard-tools-2.4.1/picard.jar"
OUT="${1}.bam"

exitcode=0
/usr/local/genome/sratoolkit.2.8.1-centos_linux64/bin/sam-dump --unaligned "${1}" \
  | ${JAVA} ${JAVAOPTS} -jar ${PICARD} SamFormatConverter I=/dev/stdin O="${OUT}"
for i in ${PIPESTATUS[@]}; do ((exitcode+=i)); done

exit ${exitcode}
