#!/bin/bash -i
#Pipe.15.intansv.sh
#by HIRAO Akira

set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0) && pwd)

CPU=6

target_ID=AT48


module load R


main_folder=/zfs/Arabidopsis/work/At_Reseq
reference_folder=/zfs/Arabidopsis/Reference_v1.1
vcf_folder=$main_folder/vcf_out
work_folder=$main_folder/intansv_out


cd $work_folder

while read sample; do

	mkdir -p $sample
	cd $sample

	Rscript $SCRIPT_DIR/Do.intansv.R $sample

	cd ../
#done < $SCRIPT_DIR/sample_ID.list #list of samples
#done < $SCRIPT_DIR/sample_ID_test.list #list of samples
done < $SCRIPT_DIR/sample_ID_test_AT01.list #list of samples

cd $SCRIPT_DIR

	
module unload gatk/4.1.7.0
module unload vcftools/0.1.15
