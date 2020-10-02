#!/bin/bash -i
#Pipe.09.VariantFiltration.sh
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


gatk VariantFiltration\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.nonvariant.vcf.gz\
 -G-filter "GQ < 20"\
 -G-filter-name "lowGQ"\
 -G-filter "DP < 10 || DP > 200"\
 -G-filter-name "DP_10-200"\
 -O $target_ID.nonvariant.DPfilterPASSED.vcf.gz


#Set filtered sites to no call:SNP
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.nonvariant.DPfilterPASSED.vcf.gz\
 --set-filtered-gt-to-nocall\
 -O $target_ID.nonvariant.DPfilterNoCall.vcf.gz


cd $SCRIPT_DIR

module unload samtools/1.10
module unload gatk/4.1.7.0
module unload bedtools2/2.27.1
module unload bcftools/1.9

