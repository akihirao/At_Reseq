#!/bin/bash -i
#Pipe.10.VariantFiltration.sh
#by HIRAO Akira

set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0) && pwd)

CPU=8

target_ID=AT48

reference_folder=/zfs/Arabidopsis/Reference_v1.1
main_folder=/zfs/Arabidopsis/work/At_Reseq
work_folder=$main_folder/vcf_out


module load samtools/1.10
module load gatk/4.1.7.0
module load bedtools2/2.27.1
module load bcftools/1.9


cd $work_folder


#Selecting AT01-AT03 snp VCFs
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.snp.DPfilterNoCall.vcf.gz\
 --sample-name AT01\
 -O AT01.snp.DPfilterNoCall.vcf

gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.snp.DPfilterNoCall.vcf.gz\
 --sample-name AT02\
 -O AT02.snp.DPfilterNoCall.vcf

gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.snp.DPfilterNoCall.vcf.gz\
 --sample-name AT03\
 -O AT03.snp.DPfilterNoCall.vcf


#Selecting AT01-AT03 indel VCFs
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.indel.DPfilterNoCall.vcf.gz\
 --sample-name AT01\
 -O AT01.indel.DPfilterNoCall.vcf

gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.indel.DPfilterNoCall.vcf.gz\
 --sample-name AT02\
 -O AT02.indel.DPfilterNoCall.vcf

gatk SelectVariants\
 -V $target_ID.indel.DPfilterNoCall.vcf.gz\
 --sample-name AT03\
 -O AT03.indel.DPfilterNoCall.vcf


cd $SCRIPT_DIR

module unload samtools/1.10
module unload gatk/4.1.7.0
module unload bedtools2/2.27.1
module unload bcftools/1.9

