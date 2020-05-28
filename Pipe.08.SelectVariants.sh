#!/bin/bash -i
#Pipe.08.SelectVariants.sh
#by HIRAO Akira

set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)

CPU=16

target_ID=AT48

reference_folder=/zfs/Arabidopsis/Reference_v1.1
maing_folder=/zfs/Arabidopsis/work/At_Reseq
work_folder=$main_folder/SNV_call


module load gatk/4.1.7.0

GQ_threshold=20

cd $output_folder

#----------------------------------------------------------------------------------
#Output mendelian violation site only: --pedigree & --mendelian-violation
gatk SelectVariants \
-R $reference_folder/TAIR10.fa \
-V $work_folder/$target_ID.vcf.gz \
--pedigree $SCRIPT_DIR/AT48.ped \
--mendelian-violation \
--mendelian-violation-qual-threshold $GQ_threshold \
-O $work_folder/$target_ID.mutation_candidates.vcf.gz

#Spliting SNP
gatk SelectVariants \
-R $reference_folder/TAIR10.fa \
-V $work_folder/$target_ID.mutation_candidates.vcf.gz \
-select-type SNP \
-O $work_folder/$target_ID.mutation_candidates.snp.vcf.gz

#Spliting INDEL
gatk SelectVariants \
-R $reference_folder/TAIR10.fa \
-V $work_folder/$target_ID.mutation_candidates.vcf.gz \
-select-type INDEL \
-O $work_folder/$target_ID.mutation_candidates.indel.vcf.gz
#----------------------------------------------------------------------------------

#----------------------------------------------------------------------------------
#Output all sites
#Spliting SNP
gatk SelectVariants \
-R $reference_folder/TAIR10.fa \
-V $work_folder/$target_ID.vcf.gz \
-select-type SNP \
-O $work_folder/$target_ID.snp.vcf.gz

#Spliting INDEL
gatk SelectVariants \
-R $reference_folder/TAIR10.fa \
-V $work_folder/$target_ID.vcf.gz \
-select-type INDEL \
-O $work_folder/$target_ID.indel.vcf.gz
#----------------------------------------------------------------------------------


cd $SCRIPT_DIR

module unload gatk/4.1.7.0


