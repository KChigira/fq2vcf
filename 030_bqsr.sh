#! /bin/sh

if [ $# -ne 4 ]; then
  echo "Error: 4 arguments must be contained."
  exit 1;
fi

REF=$1
INPUT=$2
OUTDIR=$3
NAME=$4

if [ ! -e ${REF} ] || [ ! -e ${INPUT} ]; then
  echo "Error: Input file does not exist."
  exit 1;
fi

if [ ! -d ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi

gatk HaplotypeCaller \
       -R ${REF} \
       -I ${INPUT} \
       -O ${OUTDIR}/${NAME}_raw_variants.vcf

gatk SelectVariants \
       -R ${REF} \
       -V ${OUTDIR}/${NAME}_raw_variants.vcf \
       -select-type SNP \
       -O ${OUTDIR}/${NAME}_raw_snps.vcf

gatk SelectVariants \
       -R ${REF} \
       -V ${OUTDIR}/${NAME}_raw_variants.vcf \
       -select-type INDEL \
       -O ${OUTDIR}/${NAME}_raw_indels.vcf

gatk VariantFiltration \
       -R ${REF} \
       -V ${OUTDIR}/${NAME}_raw_snps.vcf \
       -O ${OUTDIR}/${NAME}_filtered_snps.vcf \
       -filter-name "QD_filter" -filter "QD < 20.0" \
       -filter-name "FS_filter" -filter "FS > 60.0" \
       -filter-name "MQ_filter" -filter "MQ < 40.0" \
       -filter-name "SOR_filter" -filter "SOR > 4.0" \
       -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
       -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \

gatk VariantFiltration \
       -R ${REF} \
       -V ${OUTDIR}/${NAME}_raw_indels.vcf \
       -O ${OUTDIR}/${NAME}_filtered_indels.vcf \
       -filter-name "QD_filter" -filter "QD < 20.0" \
       -filter-name "FS_filter" -filter "FS > 200.0" \
       -filter-name "SOR_filter" -filter "SOR > 10.0"

gatk SelectVariants --exclude-filtered  \
       -V ${OUTDIR}/${NAME}_filtered_snps.vcf \
       -O ${OUTDIR}/${NAME}_bqsr_snps.vcf \

gatk SelectVariants --exclude-filtered  \
       -V ${OUTDIR}/${NAME}_filtered_indels.vcf \
       -O ${OUTDIR}/${NAME}_bqsr_indels.vcf \

gatk BaseRecalibrator \
       -R ${REF} \
       -I ${INPUT} \
       --known-sites ${OUTDIR}/${NAME}_bqsr_snps.vcf \
       --known-sites ${OUTDIR}/${NAME}_bqsr_indels.vcf \
       -O ${OUTDIR}/${NAME}_recal_data.table

gatk ApplyBQSR \
       -R ${REF} \
       -I ${INPUT} \
       -bqsr ${OUTDIR}/${NAME}_recal_data.table \
       -O ${OUTDIR}/${NAME}_recal_reads.bam

gatk BaseRecalibrator \
       -R ${REF} \
       -I ${OUTDIR}/${NAME}_recal_reads.bam \
       --known-sites ${OUTDIR}/${NAME}_bqsr_snps.vcf \
       --known-sites ${OUTDIR}/${NAME}_bqsr_indels.vcf \
       -O ${OUTDIR}/${NAME}_post_recal_data.table

gatk AnalyzeCovariates \
       -before ${OUTDIR}/${NAME}_recal_data.table \
       -after ${OUTDIR}/${NAME}_post_recal_data.table \
       -plots ${OUTDIR}/${NAME}_recalibration_plots.pdf \
       -csv ${OUTDIR}/${NAME}_recalibration_plots.csv
