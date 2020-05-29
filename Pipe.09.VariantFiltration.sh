#!/bin/bash -i
#Pipe.09.VariantFiltration.sh
#by HIRAO Akira

set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)

CPU=8

target_ID=AT48

reference_folder=/zfs/Arabidopsis/Reference_v1.1
maing_folder=/zfs/Arabidopsis/work/At_Reseq
work_folder=$main_folder/vcf_out

module load samtools/1.10
module load gatk/4.1.7.0


cd $work_folder


#=====================================================================
#<< COMMENTOUT1
#Only medelian violation site only
#VariantFiltration for SNP
gatk VariantFiltration \
-R $reference_folder/TAIR10.fa \
-V $target_ID.family.mendelian.snp.vcf.gz \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0" \
--filter-name "basic_snp_filter" \
-O $target_ID.mendelian.snp.filter.vcf

#grep -E '^#|PASS' $target_ID.mendelian.snp.filter.vcf > $target_ID.mendelian.snp.filterPASSED.vcf
grep -E '^#|PASS' $target_ID.mendelian.snp.filter.vcf > $target_ID.mendelian.snp.filterPASSED.vcf

#VariantFiltration for INDEL
gatk VariantFiltration \
-R $reference_folder/TAIR10.fa \
-V $target_ID.mendelian.indel.vcf.gz \
--filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" \
--filter-name "basic_indel_filter" \
-O $target_ID.mendelian.indel.filter.vcf

#grep -E '^#|PASS' $target_ID.mendelian.indel.filter.vcf > $target_ID.mendelian.indel.filterPASSED.vcf
grep -E '^#|PASS' $target_ID.mendelian.indel.filter.vcf > $target_ID.mendelian.indel.filterPASSED.vcf

#COMMENTOUT1
#---------------------------------------------------------------------------------------------------------------
#all sites
#VariantFiltration for SNP
#<< COMMENTOUT2
gatk VariantFiltration \
-R $reference_folder/TAIR10.fa \v
-V $target_ID.snp.vcf.gz \
--filter-expression "QD < 2.0 || FvS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0" \
--filter-name "basic_snp_filter" \
-O $target_ID.snp.filter.vcf

grep -E '^#|PASS' $target_ID.snp.filter.vcf > $target_ID.snp.filterPASSED.vcf

#VariantFiltration for INDEL
gatk VariantFiltration \
-R $reference_folder/TAIR10.fa \
-V $target_ID.indel.vcf.gz \
--filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" \
--filter-name "basic_indel_filter" \
-O $target_ID.indel.filter.vcf

grep -E '^#|PASS' $target_ID.indel.filter.vcf > $target_ID.indel.filterPASSED.vcf

#COMMENTOUT2
#=====================================================================


#=====================================================================
#<< COMMENTOUT3

#Only medelian violation site only
#DepthFiltering for SNP: DP < 10 & DP > 200 & GQ <20 (Low Genotype Quality: less than 99%)
gatk VariantFiltration \
-R $reference_folder/TAIR10.fa \
-V $target_ID.mendelian.snp.filterPASSED.vcf \
-G-filter "GQ < 20" \
-G-filter-name "lowGQ" \
-G-filter "DP < 10 || DP > 200" \
-G-filter-name "DP_10-200" \
-O $target_ID.mendelian.snp.DPfilterPASSED.vcf

#DepthFiltering for INDEL: DP < 10 & DP > 200 & GQ < 20 (Low Genotype Quality: less than 99%) 
gatk VariantFiltration \
-R $reference_folder/TAIR10.fa \
-V $target_ID.mendelian.indel.filterPASSED.vcf \
-G-filter "GQ < 20" \
-G-filter-name "lowGQ" \
-G-filter "DP < 10 || DP > 200" \
-G-filter-name "DP_10-200" \
-O $target_ID.mendelian.indel.DPfilterPASSED.vcf

#Set filtered sites to no call:SNP
gatk SelectVariants \
-R $reference_folder/TAIR10.fa \
-V $target_ID.mendelian.snp.DPfilterPASSED.vcf \
--set-filtered-gt-to-nocall \
-O $target_ID.mendelian.snp.DPfilterNoCall.vcf
bgzip -c $target_ID.mendelian.snp.DPfilterNoCall.vcf > $target_ID.mendelian.snp.DPfilterNoCall.vcf.gz
tabix -f -p vcf $target_ID.mendelian.snp.DPfilterNoCall.vcf.gz
bcftools index -f $target_ID.mendelian.snp.DPfilterNoCall.vcf.gz

#Set filtered sites to no call:INDEL
gatk SelectVariants \
-R $reference_folder/TAIR10.fa \
-V $target_ID.mendelian.indel.DPfilterPASSED.vcf \
--set-filtered-gt-to-nocall \
-O $target_ID.mendelian.indel.DPfilterNoCall.vcf
bgzip -c $target_ID.mendelian.indel.DPfilterNoCall.vcf > $target_ID.mendelian.indel.DPfilterNoCall.vcf.gz
tabix -f -p vcf $target_ID.mendelian.indel.DPfilterNoCall.vcf.gz
bcftools index -f $target_ID.mendelian.indel.DPfilterNoCall.vcf.gz

#COMMENTOUT3

#---------------------------------------------------------------------------------------------------------------
#all sites
#DepthFiltering for SNP: DP < 10 & DP > 200 & GQ <20 (Low Genotype Quality: less than 99%)
gatk VariantFiltration \
-R $reference_folder/TAIR10.fa \
-V $target_ID.snp.filterPASSED.vcf \
-G-filter "GQ < 20" \
-G-filter-name "lowGQ" \
-G-filter "DP < 10 || DP > 200" \
-G-filter-name "DP_10-200" \
-O $target_ID.snp.DPfilterPASSED.vcf

#DepthFiltering for INDEL: DP < 10 & DP > 200 & GQ < 20 (Low Genotype Quality: less than 99%) 
gatk VariantFiltration \
-R $reference_folder/TAIR10.fa \
-V $target_ID.indel.filterPASSED.vcf \
-G-filter "GQ < 20" \
-G-filter-name "lowGQ" \
-G-filter "DP < 10 || DP > 200" \
-G-filter-name "DP_10-200" \
-O $target_ID.indel.DPfilterPASSED.vcf

#Set filtered sites to no call:SNP
gatk SelectVariants \
-R $reference_folder/TAIR10.fa \
-V $target_ID.family.snp.DPfilterPASSED.vcf \
--set-filtered-gt-to-nocall \
-O $target_ID.snp.DPfilterNoCall.vcf
bgzip -c $target_ID.snp.DPfilterNoCall.vcf > $target_ID.snp.DPfilterNoCall.vcf.gz
tabix -f -p vcf $target_ID.snp.DPfilterNoCall.vcf.gz
bcftools index -f $target_ID.snp.DPfilterNoCall.vcf.gz

#Set filtered sites to no call:INDEL
gatk SelectVariants \
-R $reference_folder/TAIR10.fa \
-V $target_ID.indel.DPfilterPASSED.vcf \
--set-filtered-gt-to-nocall \
-O $target_ID.indel.DPfilterNoCall.vcf
bgzip -c $target_ID.indel.DPfilterNoCall.vcf > $target_ID.indel.DPfilterNoCall.vcf.gz
tabix -f -p vcf $target_ID.indel.DPfilterNoCall.vcf.gz
bcftools index -f $target_ID.indel.DPfilterNoCall.vcf.gz

#=====================================================================


cd $SCRIPT_DIR

module unload samtools/1.10
module unload gatk/4.1.7.0


