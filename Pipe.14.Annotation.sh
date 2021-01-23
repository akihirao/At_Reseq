#!/bin/bash -i
#Pipe.14.Annotation.sh
#by HIRAO Akira
#requirement Vcf_SnpEff_AT_Gene.pl (making annotation table corresponding with M2.mutations.full.list.csv)


set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)


reference_folder=/zfs/Arabidopsis/Reference_v1.1
main_folder=/zfs/Arabidopsis/work/At_Reseq
pindel_folder=$main_folder/pindel_out
vcf_folder=$main_folder/vcf_out
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
	 -V $vcf_folder/M2.snp.indel.gatk.pindl.common.homo.hetero.familyclustered.vcf.gz\
	 --set-filtered-gt-to-nocall\
	 --sample-name $target_sample\
	 -O  $vcf_folder/$target_sample/$target_sample.homo.hetero.familyclustered.pre.vcf.gz

	gatk SelectVariants\
	 -R $reference_folder/TAIR10.fa\
	 -V $vcf_folder/$target_sample/$target_sample.homo.hetero.familyclustered.pre.vcf.gz\
	 -select "vc.getGenotype('${target_sample}').isHomVar() || vc.getGenotype('${target_sample}').isHet() "\
	 -O $vcf_folder/$target_sample/$target_sample.homo.hetero.familyclustered.vcf

	java -Xmx4g -jar /usr/local/snpEff/snpEff.jar Arabidopsis_thaliana $vcf_folder/$target_sample/$target_sample.homo.hetero.familyclustered.vcf >$vcf_folder/$target_sample/$target_sample.final.mutations.snpeff.vcf

	Rscript $SCRIPT_DIR/snpEff_effect_summaryzing.R $target_sample

	cat $target_sample.snpEff_effect_summary.txt | tail -n +2 >> $SCRIPT_DIR/snpEff_all_samples_summary_file.txt

	cd ../

done

cd $SCRIPT_DIR

perl $SCRIPT_DIR/Vcf_SnpEff_AT_Gene.pl
(head -n +1 SnpEff.AT.gene.unsorted.txt && tail -n +2 SnpEff.AT.gene.unsorted.txt | sort -k 1,1 -k 2n,2) | uniq > SnpEff.AT.gene.txt
rm SnpEff.AT.gene.unsorted.txt

module unload java
module unload R
module unload gatk/4.1.7.0
