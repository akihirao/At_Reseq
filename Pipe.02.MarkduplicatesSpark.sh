#!/bin/bash -i
#Pipe.02.MarkduplicatesSpark.sh
#by HIRAO Akira


set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0) && pwd)

CPU=16

reference_folder=/zfs/Arabidopsis/Reference_v1.1
main_folder=/zfs/Arabidopsis/work/At_Reseq

module load samtools/1.10
module load gatk/4.1.7.0


while read sample fastq_R1 fastq_R2; do

	work_folder=$main_folder/bwa_out/$sample
	cd $work_folder

	#marking PCR duplicates
	gatk MarkDuplicatesSpark -I $sample.TAIR10.bam -M $sample.TAIR10.metrics.txt -O $sample.TAIR10.MarkDup.bam
	#removing multiple-mapped reads
	samtools view -@ $CPU -b -F 256 -q 4 $sample.TAIR10.MarkDup.bam > $sample.TAIR10.filtered.bam
	samtools index -@ $CPU $sample.TAIR10.filtered.bam

done < $SCRIPT_DIR/fastq_list.txt  #list of MIDs and folder deposited their read fastq.gz


cd $SCRIPT_DIR


module unload samtools/1.10
module unload gatk/4.1.7.0
