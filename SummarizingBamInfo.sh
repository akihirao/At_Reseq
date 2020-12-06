#!/bin/bash -i
#SummarizingBamInfo.sh
#by HIRAO Akira


set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)


reference_folder=/zfs/Arabidopsis/Reference_v1.1
main_folder=/zfs/Arabidopsis/work/At_Reseq
vcf_folder=$main_folder/vcf_out
bwa_folder=$main_folder/bwa_out
work_folder=$main_folder/bwa_out


echo -n >| $SCRIPT_DIR/BAM.info.summary.txt

#bamtoools. 2.5.1
module load bamtools

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

cd $work_folder


for target_sample in ${M2_vec[@]}
do

	cd $target_sample

	#Sumamarizing BAM info
	bamtools stats -in $target_sample.TAIR10.bam > $target_sample.TAIR10.bam.stats.txt
	mapping_rate=$(grep "Mapped reads" $target_sample.TAIR10.bam.stats.txt)
	both_mapping_rate=$(grep "Both pairs mapped" $target_sample.TAIR10.bam.stats.txt)

	echo $target_sample $mapping_rate $both_mapping_rate >> $SCRIPT_DIR/BAM.info.summary.txt
	
	cd ../

done

cd $SCRIPT_DIR

module unload bamtools

