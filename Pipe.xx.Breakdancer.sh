#!/bin/bash -i
#Pipe.xx.Breakdancer.sh
#by HIRAO Akira


set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)

CPU=16


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

module load miniconda2
module load vcftools/0.1.15
module load gatk/4.1.7.0
module load breakdancer

cd $work_folder


perl /usr/local/breakdancer/lib/breakdancer-maxunstable/bam2cfg.pl $bwa_folder/
#-----------------------------------------------------



cd $SCRIPT_DIR

module unload miniconda2
module unload vcftools/0.1.15
module unload gatk/4.1.7.0


