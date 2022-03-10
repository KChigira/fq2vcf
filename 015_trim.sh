#! /bin/bash
#trimmomatic-0.39.jar must be installed at home directory

if [ $# -ne 5 ]; then
  echo "Error: 5 arguments must be contained."
  exit 1;
fi
if [ ! -e ${1} ] || [ ! -e ${2} ]; then
  echo "Error: Input file does not exist."
  exit 1;
fi
if [ ! -d ${3} ]; then
  mkdir ${3}
fi
if expr "${5}" : "[0-9]*$" >&/dev/null; then
	echo "OK"
else
  echo "Error: Argument 5 must be numeric."
	exit 1;
fi

java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar \
					PE -threads 8 -phred33 -trimlog ${3}/${4}_trimlog.txt \
					${1} ${2} \
					${3}/${4}_trimed_1.fastq.gz \
					${3}/${4}_trimed_unpaired_1.fastq.gz \
					${3}/${4}_trimed_2.fastq.gz \
					${3}/${4}_trimed_unpaired_2.fastq.gz \
					CROP:${5} \
					LEADING:20 \
					TRAILING:20 \
					SLIDINGWINDOW:4:15 \
					MINLEN:50

rm ${3}/${4}_trimlog.txt
