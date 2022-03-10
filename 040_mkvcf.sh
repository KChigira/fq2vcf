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
       -O ${OUTDIR}/${NAME}_raw_variants_final.vcf

gatk SelectVariants \
       -R ${REF} \
       -V ${OUTDIR}/${NAME}_raw_variants_final.vcf \
       -select-type SNP \
       -O ${OUTDIR}/${NAME}_raw_snps_final.vcf

gatk SelectVariants \
       -R ${REF} \
       -V ${OUTDIR}/${NAME}_raw_variants_final.vcf \
       -select-type INDEL \
       -O ${OUTDIR}/${NAME}_raw_indels_final.vcf

gatk VariantFiltration \
       -R ${REF} \
       -V ${OUTDIR}/${NAME}_raw_snps_final.vcf \
       -O ${OUTDIR}/${NAME}_final_snps.vcf \
       -filter-name "QD_filter" -filter "QD < 20.0" \
       -filter-name "FS_filter" -filter "FS > 60.0" \
       -filter-name "MQ_filter" -filter "MQ < 40.0" \
       -filter-name "SOR_filter" -filter "SOR > 4.0" \
       -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
       -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \

gatk VariantFiltration \
       -R ${REF} \
       -V ${OUTDIR}/${NAME}_raw_indels_final.vcf \
       -O ${OUTDIR}/${NAME}_final_indels.vcf \
       -filter-name "QD_filter" -filter "QD < 20.0" \
       -filter-name "FS_filter" -filter "FS > 200.0" \
       -filter-name "SOR_filter" -filter "SOR > 10.0"

#
bgzip -c ${OUTDIR}/${NAME}_final_snps.vcf \
      > ${OUTDIR}/${NAME}_final_snps.vcf.gz
bcftools index ${OUTDIR}/${NAME}_final_snps.vcf.gz

bgzip -c ${OUTDIR}/${NAME}_final_indels.vcf \
      > ${OUTDIR}/${NAME}_final_indels.vcf.gz
bcftools index ${OUTDIR}/${NAME}_final_indels.vcf.gz

bcftools concat -o ${OUTDIR}/${NAME}_final_variants.vcf.gz \
                -a -O z \
                ${OUTDIR}/${NAME}_final_snps.vcf.gz \
                ${OUTDIR}/${NAME}_final_indels.vcf.gz
