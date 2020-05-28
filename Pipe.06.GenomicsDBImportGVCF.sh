#!/bin/bash -i
#Pipe.06.GenomicsDBImportGVCF.sh
#by HIRAO Akira

set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0) && pwd)

CPU=8

module load gatk/4.1.7.0


main_folder=/zfs/Arabidopsis/work/At_Reseq
reference_folder=/zfs/Arabidopsis/Reference_v1.1

#-----------------------------------------------------
# defining argument of samples for GenomicsDBImports
input_samples=""
option_lab="-V "
TAIR10_gvcf_lab=".TAIR10.g.vcf.gz"
one_space=" "
slash_lab="/"

while read sample; do

	echo $sample
	gvcf_folder=$main_folder/SNV_call/$sample$slash_lab
		
	input_samples=$input_samples$option_lab$gvcf_folder$sample$TAIR10_gvcf_lab$one_space

done < $SCRIPT_DIR/sample_list.txt #list of samples

echo $input_samples

#-----------------------------------------------------


main_folder=/zfs/Arabidopsis/work/At_Reseq
reference_folder=/zfs/Arabidopsis/Reference_v1.1


target_ID=AT48
output_folder=$main_folder/gDB
mkdir -p $output_folder

genomicsDB_name=genomicsDB.$target_ID
DB_path=$output_folder/$genomicsDB_name


cd $output_folder


gatk GenomicsDBImport \
$input_samples \
--genomicsdb-workspace-path  $DB_path \
--intervals $SCRIPT_DIR/AT_Chr1-5.list \
--reader-threads $CPU


cd $SCRIPT_DIR


	
module unload gatk/4.1.7.0

