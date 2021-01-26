#!/bin/bash -i
#Pipe.13.UnifyGATKPindel.sh
#by HIRAO Akira

set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0) && pwd)

CPU=6

target_ID=AT48


module load gatk/4.1.7.0
module load vcftools/0.1.15


main_folder=/zfs/Arabidopsis/work/At_Reseq
reference_folder=/zfs/Arabidopsis/Reference_v1.1
BioAlcidaeJdk_path=/usr/local/jvarkit/dist
vcf_folder=$main_folder/vcf_out
pindel_folder=$main_folder/pindel_out
work_folder=$main_folder/pindel_out

BioAlcidaeJdk_path=/usr/local/jvarkit/dist

cd $work_folder

#compare indel variants identified by gatk with those by pindel: all
bcftools isec $vcf_folder/M2.indel.all.vcf.gz AT.M2.unique.pindel.D.SI.vcf.gz -p AT.M2.gatk.pindel.common.all -n=2
perl $SCRIPT_DIR/Vcf2BED_chr_start_end.pl < $work_folder/AT.M2.gatk.pindel.common.all/0000.vcf > $work_folder/AT.M2.indel.common.pindel.gatk.all.bed

cat $vcf_folder/M2.snp.all.bed $work_folder/AT.M2.indel.common.pindel.gatk.all.bed | sort -k 1,1 -k 2n,2 >  $work_folder/M2.snp.indel.gatk.pindl.common.all.bed


#select back variants with all.bed 
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $vcf_folder/$target_ID.mu.snp.indel.DPfilterNoCall.vcf.gz\
 -L $work_folder/M2.snp.indel.gatk.pindl.common.all.bed\
 --exclude-sample-name $SCRIPT_DIR/Mother_ID.list\
 -O $vcf_folder/M2.snp.indel.gatk.pindl.common.hetero.homo.familyclustered.vcf
bgzip -c $vcf_folder/M2.snp.indel.gatk.pindl.common.hetero.homo.familyclustered.vcf > $vcf_folder/M2.snp.indel.gatk.pindl.common.hetero.homo.familyclustered.vcf.gz
tabix -f -p vcf $vcf_folder/M2.snp.indel.gatk.pindl.common.hetero.homo.familyclustered.vcf.gz
	

perl $SCRIPT_DIR/Vcf2List_R_glm.pl < $vcf_folder/M2.snp.indel.gatk.pindl.common.hetero.homo.familyclustered.vcf > $SCRIPT_DIR/M2.mutations.full.list.csv

cd $SCRIPT_DIR

	
module unload gatk/4.1.7.0
module unload vcftools/0.1.15

