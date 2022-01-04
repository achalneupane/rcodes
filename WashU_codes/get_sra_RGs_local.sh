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

sam-dump "${SRR}.sra" | sed -n '/^[^@]/!p;//q' | grep ^@RG

exit $?
