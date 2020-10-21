#!/bin/bash -i
#Pipe.11.Annotation.sh
#by HIRAO Akira


set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)


reference_folder=/zfs/Arabidopsis/Reference_v1.1
main_folder=/zfs/Arabidopsis/work/At_Reseq
pindel_folder=$main_folder/pindel_out
work_folder=$main_folder/vcf_out


echo -n >| snpEff_all_samples_summary_file.txt

module load java
module load R
module load gatk/4.1.7.0

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

	gatk SelectVariants\
	 -R $reference_folder/TAIR10.fa\
	 -V $target_sample.all.vcf.gz\
	 -XL $pindel_folder/AT.M2.filterout.pindel.gatk.indel.bed\
	 -O $target_sample.pindel_gatk_common.mutants.vcf

	java -Xmx4g -jar /usr/local/snpEff/snpEff.jar Arabidopsis_thaliana $target_sample.pindel_gatk_common.mutants.vcf > $target_sample.final.mutants.snpeff.vcf

	Rscript $SCRIPT_DIR/snpEff_effect_summaryzing.R $target_sample

	cat $target_sample.snpEff_effect_summary.txt | tail -n +2 >> $SCRIPT_DIR/snpEff_all_samples_summary_file.txt

	cd ../

done

cd $SCRIPT_DIR

module unload java
module unload R
module unload gatk/4.1.7.0
