#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Error, must specify SRR number"
    exit 2
fi

SRR="${1}"
sam-dump "${1}" | sed -n '/^[^@]/!p;//q' | grep ^@RG

exit $?
