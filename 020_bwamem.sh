#! /bin/sh
#bwa:0.7.17
#gatk4-4.2.4.0
#gatk4 does not work in jdk-11.
#use jdk-8, choosing by "sudo update-alternatives --config java"
#picard-2.26.10
#CollectAlignmentSummaryMetrics requires R.
#R must be installed in local, not in conda environment.
#If r-base is installed in conda environment, It must be uninstalled.

if [ $# -ne 5 ]; then
  echo "Error: 5 arguments must be contained."
  exit 1;
fi
if [ ! -e ${1} ] || [ ! -e ${2} ] || [ ! -e ${3} ]; then
  echo "Error: Input file does not exist."
  exit 1;
fi
if [ ! -d ${4} ]; then
  mkdir ${4}
fi

REF=$1
IN_1=$2
IN_2=$3
OUTDIR=$4
NAME=$5

bwa mem \
    -t 10 -K 100000000 \
    -Y \
    -R "@RG\tID:${NAME}\tLB:${NAME}\tPL:ILLUMINA\tSM:${NAME}" \
    ${REF} ${IN_1} ${IN_2} \
    > ${OUTDIR}/${NAME}_aligned_reads.sam

gatk MarkDuplicatesSpark \
     -I ${OUTDIR}/${NAME}_aligned_reads.sam \
     -M ${OUTDIR}/${NAME}_dedup_metrics.txt \
     -O ${OUTDIR}/${NAME}_sorted_reads.bam

rm ${OUTDIR}/${NAME}_aligned_reads.sam

java -jar ~/picard/picard.jar CollectAlignmentSummaryMetrics \
       -R ${REF} \
       -I ${OUTDIR}/${NAME}_sorted_reads.bam \
       -O ${OUTDIR}/${NAME}_alignment_metrics.txt

java -jar ~/picard/picard.jar CollectInsertSizeMetrics \
       -I ${OUTDIR}/${NAME}_sorted_reads.bam \
       -O ${OUTDIR}/${NAME}_insert_metrics.txt \
       -H ${OUTDIR}/${NAME}_insert_size_histgram.pdf
