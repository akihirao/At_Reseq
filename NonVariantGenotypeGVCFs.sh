#!/bin/bash -i
#NonVariantGenotypeGVCFs.sh
#For assessing all of callable regions including non variant sites
#by HIRAO Akira

set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0) && pwd)

CPU=8

module load gatk/4.1.7.0


main_folder=/zfs/Arabidopsis/work/At_Reseq
reference_folder=/zfs/Arabidopsis/Reference_v1.1


target_ID=AT48
output_folder=$main_folder/vcf_out
mkdir -p $output_folder

genomicsDB_name=genomicsDB.$target_ID
DB_path=$main_folder/gDB/$genomicsDB_name


cd $output_folder

gatk GenotypeGVCFs \
-R $reference_folder/TAIR10.fa \
-V gendb://$DB_path \
--include-non-variant-sites \
-O $output_folder/$target_ID.nonvariant.vcf.gz

cd $SCRIPT_DIR

	
module unload gatk/4.1.7.0

