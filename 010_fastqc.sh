#! /bin/sh
#fastqc:0.11.9

if [ $# -ne 2 ]; then
  echo "Error: 2 arguments must be contained."
  exit 1;
fi
if [ ! -e ${1} ]; then
  echo "Error: Input file does not exist."
  exit 1;
fi
if [ ! -d ${2} ]; then
  mkdir ${2}
fi

fastqc --nogroup -o ${2} ${1}
