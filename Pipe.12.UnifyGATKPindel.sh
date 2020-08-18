#!/bin/bash -i
#Pipe.12.UnifyGATKPindel.sh
#by HIRAO Akira

set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0) && pwd)

CPU=6

module load gatk/4.1.7.0
module load vcftools/0.1.15


main_folder=/zfs/Arabidopsis/work/At_Reseq
reference_folder=/zfs/Arabidopsis/Reference_v1.1
BioAlcidaeJdk_path=/usr/local/jvarkit/dist
vcf_folder=$main_folder/vcf_out
work_folder=$main_folder/pindel_out

cd $work_folder

#compare indel variants identified by gatk with those by ppindel
bcftools isec $vcf_folder/M2.indel.unique.vcf.gz AT.M2.unique.pindel.indel.vcf.gz -p AT.M2.gatk.pindel.common -n=2

bcftools isec $vcf_folder/M2.indel.unique.vcf.gz AT.M2.unique.pindel.indel.vcf.gz -p AT.M2.gatk.unique -C
perl $SCRIPT_DIR/Vcf2BED_chr_start_end.pl < $work_folder/AT.M2.gatk.unique/0000.vcf > $work_folder/AT.M2.filterout.pindel.gatk.indel.bed

cd $SCRIPT_DIR

	
module unload gatk/4.1.7.0
module unload vcftools/0.1.15
