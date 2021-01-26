#!/bin/bash -i
#Pipe.01.fastp.sh
#by HIRAO Akiras

set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)


CPU=8

raw_fastq_folder=/zfs/Arabidopsis/data/raw_fastq
QC_folder=/zfs/Arabidopsis/data/AT48

while read sample; do

	if [ ! -e $QC_folder/$sample ]; then
		mkdir $QC_folder/$sample
	fi
	
	R1_tag="_R1"
	R2_tag="_R2"
	sample_R1=$sample$R1_tag
	sample_R2=$sample$R2_tag
	
	fastp -i $raw_fastq_folder/$sample_R1.fastq.gz -I $raw_fastq_folder/$sample_R2.fastq.gz -3 \
	-o $QC_folder/$sample/$sample_R1.trimQ30.fastq.gz -O $QC_folder/$sample/$sample_R2.trimQ30.fastq.gz \
	-h $QC_folder/$sample/$sample.fastp.trimQ30.report.html -j $QC_folder/$sample/$sample.fastq.trimQ30.report.json -q 30 \
	-w $CPU

done < $SCRIPT_DIR/sample_ID.list #list of ID
