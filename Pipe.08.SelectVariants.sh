#!/bin/bash -i
#Pipe.08.SelectVariants.sh
#by HIRAO Akira

set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0) && pwd)

CPU=8

target_ID=AT48

reference_folder=/zfs/Arabidopsis/Reference_v1.1
main_folder=/zfs/Arabidopsis/work/At_Reseq
work_folder=$main_folder/vcf_out


module load gatk/4.1.7.0


cd $work_folder

#----------------------------------------------------------------------------------
#Spliting SNPv
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $work_folder/$target_ID.vcf.gz\
 -select-type SNP\
 --restrict-alleles-to BIALLELIC\
 -O $work_folder/$target_ID.snp.vcf.gz

#Spliting INDEL (only biallelic)
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $work_folder/$target_ID.vcf.gz\
 -select-type INDEL\
 --restrict-alleles-to BIALLELIC\
 -O $work_folder/$target_ID.indel.vcf.gz
#----------------------------------------------------------------------------------


cd $SCRIPT_DIR


module unload gatk/4.1.7.0
