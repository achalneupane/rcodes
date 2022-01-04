#!/bin/bash
if [ $# -lt 1 ]; then
  echo "Error, must provide .sra file"
  exit 2
fi

if [ ! -r "${1}" ]; then
  echo "Error, ${1} does not exist or cannot be read"
  exit 1
fi

if [[ ! ${1} =~ \.sra$ ]] ; then
  echo "Error, file must be an .sra file"
  exit 1
fi

JAVA="/usr/java/jre1.8.0_91/bin/java"
JAVAOPTS="-Xms4g -Xmx8g -XX:ParallelGCThreads=1 -Djava.io.tmpdir=/data/tmp"
PICARD="/usr/local/genome/picard-tools-2.4.1/picard.jar"
OUT="${1%.sra}.bam"

exitcode=0
/usr/local/genome/sratoolkit.2.8.0-centos_linux64/bin/sam-dump --unaligned --spot-group "${1}" \
  | ${JAVA} ${JAVAOPTS} -jar ${PICARD} SamFormatConverter I=/dev/stdin O="${OUT}"
for i in ${PIPESTATUS[@]}; do ((exitcode+=i)); done

if [ ${exitcode} -eq 0 ]; then
  rm -fv "${1}" && mv -vi "${1}.vdbcache" done/ && mv -vi "${OUT}" ready/
  for i in ${PIPESTATUS[@]}; do ((exitcode+=i)); done
elif [ ${exitcode} -ne 0 ]; then
  rm -fv "${OUT}"
  ((exitcode+=$?))
fi

exit ${exitcode}
