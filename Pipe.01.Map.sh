#!/bin/bash -i
#Pipe.01.Map.sh
#by HIRAO Akira

set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)

CPU=16

reference_folder=/zfs/Arabidopsis/Reference_v1.1
main_folder=/zfs/Arabidopsis/work/At_Reseq


module load bwa #Version: 0.7.17-r1188
module load samtools/1.10


while read sample fastq_R1 fastq_R2; do

	work_folder=$main_folder/bwa_out/$sample
	mkdir -p $work_folder
	cd $work_folder
		
	#setting RG: @RG\tID:AT48\tSM:$sample\tPL:Illumina
	specific_ID="AT48"
	tag_read_group_part1="@RG\tID:"
	tag_read_group_part2="\tSM:"
	tag_read_group_part3="\tPL:Illumina"
	tag_read_group=$tag_read_group_part1$specific_ID$tag_read_group_part2$sample$tag_read_group_part3

	bwa mem -t $CPU -M -R $tag_read_group $reference_folder/TAIR10.fa $fastq_R1 $fastq_R2 | samtools view -@ $CPU -b | samtools sort -@ $CPU > $sample.TAIR10.bam
	samtools index -@ $CPU $sample.TAIR10.bam
	
done < $SCRIPT_DIR/fastq_list.txt #list of MIDs and folder deposited their read fastq.gz

cd $SCRIPT_DIR


module unload bwa
module unload samtools/1.10
