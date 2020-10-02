#!/bin/bash -i
#Pipe.13.Pindel.indiv.sh
#by HIRAO Akira


set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)

CPU=12


reference_folder=/zfs/Arabidopsis/Reference_v1.1
main_folder=/zfs/Arabidopsis/work/At_Reseq
bwa_folder=$main_folder/bwa_out
work_folder=$main_folder/pindel_out

slash_lab="/"
tab_lab="	"
insert_length=450

mkdir -p $work_folder

module load miniconda2
module load vcftools/0.1.15
module load gatk/4.1.7.0

cd $work_folder

while read sample; do

	mkdir -p $sample
	cd $sample

	config_contents=$bwa_folder$slash_lab$sample$slash_lab$sample.TAIR10.bqsr.RG.bam$tab_lab$insert_length$tab_lab$sample
	echo $config_contents > $sample.config.txt

	#Identifying structural variants with using pindel
	#six types of structural variants（D = deletion、SI = short insertion、INV = inversion、TD = tandem duplication、LI = large insertion、BP = unassigned breakpoints）
	#if Option "-c ALL" does not work well in the environment, please work for each of chromosomes   We changed widonws size from 5 to 1 (e.g. "-w 1") due to high momory usage 
	/usr/local/miniconda2/bin/pindel -T $CPU -f $reference_folder/TAIR10.fa -i $work_folder/$sample/$sample.config.txt -c ALL -w 1 -o $work_folder/$sample/$sample.pindel.out

cd ../
done < $SCRIPT_DIR/sample_ID.list #list of samples




cd $SCRIPT_DIR

module unload miniconda2
module unload vcftools/0.1.15
module unload gatk/4.1.7.0


