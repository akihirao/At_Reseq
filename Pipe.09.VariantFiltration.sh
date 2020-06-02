#!/bin/bash -i
#Pipe.09.VariantFiltration.sh
#by HIRAO Akira

#modified to restrict to biallelic @2020530
#modified for filtering out ExcessHet P < 0.05

set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0) && pwd)

CPU=8

target_ID=AT48

reference_folder=/zfs/Arabidopsis/Reference_v1.1
main_folder=/zfs/Arabidopsis/work/At_Reseq
work_folder=$main_folder/vcf_out


module load samtools/1.10
module load gatk/4.1.7.0
#!!! module load bcftools centOS　environmnt moduleの設定

cd $work_folder

#defining ExcessHet parameter
#"ExcessHet > 13.0" means excess of heterozygosity with p value 0.05 
ExcessHet_P=0.05
ExcessHet_Q=`echo "scale=5; -10 * l($ExcessHet_P) /l(10)" |bc -l | xargs printf %.1f`
ExcessHet_param="ExcessHet > ${ExcessHet_Q}"
echo $ExcessHet_param


#VariantFiltration for SNP
gatk VariantFiltration\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.snp.vcf.gz \
 --filter-expression "QD < 2.0" --filter-name "QDlt2"\
 --filter-expression "FS > 60.0" --filter-name "FSgt60"\
 --filter-expression "MQ < 40.0" --filter-name "MQlt40"\
 --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSumltnegative12.5"\
 --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSumltnegative8"\
 --filter-expression "SOR > 4.0" --filter-name "SORgt4"\
 --filter-expression "${ExcessHet_param}" --filter-name "ExHet"\
 -O $target_ID.snp.filter.vcf

grep -E '^#|PASS' $target_ID.snp.filter.vcf | bgzip > $target_ID.snp.filterPASSED.vcf.gz
tabix -f -p vcf $target_ID.snp.filterPASSED.vcf.gz


#VariantFiltration for INDEL
gatk VariantFiltration\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.indel.vcf.gz\
 --filter-expression "QD < 2.0" --filter-name "QDlt2"\
 --filter-expression "FS > 200.0" --filter-name "FSgt200"\
 --filter-expression "ReadPosRankSum < -20.0" --filter-name "RPRSltnagative20"\
 --filter-expression "SOR > 10.0" --filter-name "SORgt10"\
 --filter-expression "${ExcessHet_param}" --filter-name "ExHet"\
 -O $target_ID.indel.filter.vcf

grep -E '^#|PASS' $target_ID.indel.filter.vcf | bgzip > $target_ID.indel.filterPASSED.vcf.gz
tabix -f -p vcf $target_ID.indel.filterPASSED.vcf.gz


#---------------------------------------------------------------------------------------------------------------
#DepthFiltering for SNP: DP < 10 & DP > 200 & GQ <20 (Low Genotype Quality: less than 99%)
gatk VariantFiltration\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.snp.filterPASSED.vcf.gz\
 -G-filter "GQ < 20"\
 -G-filter-name "lowGQ"\
 -G-filter "DP < 10 || DP > 200"\
 -G-filter-name "DP_10-200"\
 -O $target_ID.snp.DPfilterPASSED.vcf.gz

#DepthFiltering for INDEL: DP < 10 & DP > 200 & GQ < 20 (Low Genotype Quality: less than 99%) 
gatk VariantFiltration\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.indel.filterPASSED.vcf.gz\
 -G-filter "GQ < 20"\
 -G-filter-name "lowGQ"\
 -G-filter "DP < 10 || DP > 200"\
 -G-filter-name "DP_10-200"\
 -O $target_ID.indel.DPfilterPASSED.vcf.gz

#Set filtered sites to no call:SNP
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.snp.DPfilterPASSED.vcf.gz\
 --set-filtered-gt-to-nocall\
 -O $target_ID.snp.DPfilterNoCall.vcf.gz
bcftools index -f $target_ID.snp.DPfilterNoCall.vcf.gz
bcftools view $target_ID.snp.DPfilterNoCall.vcf.gz -Ov -o $target_ID.snp.DPfilterNoCall.vcf

#Set filtered sites to no call:INDEL
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.indel.DPfilterPASSED.vcf.gz\
 --set-filtered-gt-to-nocall\
 -O $target_ID.indel.DPfilterNoCall.vcf.gz
bcftools index -f $target_ID.indel.DPfilterNoCall.vcf.gz
bcftools view $target_ID.indel.DPfilterNoCall.vcf.gz -Ov -o $target_ID.indel.DPfilterNoCall.vcf
#=====================================================================


#Merge SNPs and INDELs vcf files into a SNV vcf file
gatk MergeVcfs\
 -I $target_ID.snp.DPfilterNoCall.vcf.gz\
 -I $target_ID.indel.DPfilterNoCall.vcf.gz\
 -O $target_ID.snp.indel.DPfilterNoCall.vcf.gz
bcftools index -f $target_ID.snp.indel.DPfilterNoCall.vcf.gz
bcftools view  $target_ID.snp.indel.DPfilterNoCall.vcf.gz -Ov -o  $target_ID.snp.indel.DPfilterNoCall.vcf


cd $SCRIPT_DIR

module unload samtools/1.10
module unload gatk/4.1.7.0


