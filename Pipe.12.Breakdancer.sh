#!/bin/bash -i
#Pipe.12.Breakdancer.sh
#by HIRAO Akira


set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)

CPU=12


reference_folder=/zfs/Arabidopsis/Reference_v1.1
main_folder=/zfs/Arabidopsis/work/At_Reseq
bwa_folder=$main_folder/bwa_out
work_folder=$main_folder/breakdancer_out


#-----------------------------------------------------
# defining the argument for 48 samples
samples_48=()

while read sample; do
	echo $sample
	samples_48+=($sample)
done < $SCRIPT_DIR/sample_ID.list #list of samples

echo $samples_48

mother_vec=(${samples_48[@]:0:12})
M2_vec=(${samples_48[@]:12:36})
echo ${mother_vec[@]}
echo ${M2_vec[@]}
#-----------------------------------------------------


mkdir -p $work_folder

module load vcftools/0.1.15
module load gatk/4.1.7.0
module load samtools/1.10
module load breakdancer

cd $bwa_folder

while read sample; do

	cd $sample
#	gatk AddOrReplaceReadGroups --INPUT $sample.TAIR10.bqsr.bam --OUTPUT $sample.TAIR10.bqsr.RG.bam --RGID AT48 --RGLB $sample --RGPL Illumina --RGSM $sample --RGPU cell
#	samtools index -@ $CPU $sample.TAIR10.bqsr.RG.bam
	bam2cfg.pl $sample.TAIR10.bqsr.RG.bam > $work_folder/$sample.breakdancer.cfg
	breakdancer-max  $work_folder/$sample.breakdancer.cfg >  $work_folder/$sample.breakdancer.out
cd ../
done < $SCRIPT_DIR/sample_ID.list #list of samples
#done < $SCRIPT_DIR/sample_ID_test.list #list of samples

cd $SCRIPT_DIR





#breakdancer-max  $work_folder/AT_breakdancer_config.txt >  $work_folder/AT_breakdancer_output.txt
#-----------------------------------------------------



cd $SCRIPT_DIR


module unload vcftools/0.1.15
module unload gatk/4.1.7.0
module unload samtools/1.10
module unload breakdancer


