#!/bin/bash

if [ $# -lt 1 ]
then
    echo ""
    echo "Usage: $0 <output prefix>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../mamba/bin:${PATH}

OUTP=${1}
python ${BASEDIR}/../scripts/qc.py -p ${OUTP} > ${OUTP}.qc.summary

