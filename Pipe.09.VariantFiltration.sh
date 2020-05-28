#!/bin/bash -i
#Pipe.09.VariantFiltration.sh
#by HIRAO Akira


set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)

CPU=16

reference_folder=/zfs/Arabidopsis/Reference_v1.1

module load samtools/1.10
module load gatk/4.1.7.0



target_ID=AT48

reference_folder=/zfs/Arabidopsis/Reference_v1.1
input_folder=/zfs/Arabidopsis/work/FamilySNV/$target_ID.family
output_folder=/zfs/Arabidopsis/work/FamilySNV/$target_ID.family

cd $output_folder



#=====================================================================
#<< COMMENTOUT1
#Only medelian violation site only
#VariantFiltration for SNP
gatk VariantFiltration \
-R $reference_folder/TAIR10.fa \
-V $target_ID.family.mutation_candidates.snp.vcf.gz \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0" \
--filter-name "basic_snp_filter" \
-O $target_ID.family.mutation_candidates.snp.filter.vcf

#grep -E '^#|PASS' $target_ID.family.mutation_candidates.snp.filter.vcf > $target_ID.family.mutation_candidates.snp.filterPASSED.vcf
grep -E '^#|PASS' $target_ID.family.mutation_candidates.snp.filter.vcf > $target_ID.family.snp.filterPASSED.vcf

#VariantFiltration for INDEL
gatk VariantFiltration \
-R $reference_folder/TAIR10.fa \
-V $target_ID.family.mutation_candidates.indel.vcf.gz \
--filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" \
--filter-name "basic_indel_filter" \
-O $target_ID.family.mutation_candidates.indel.filter.vcf

#grep -E '^#|PASS' $target_ID.family.mutation_candidates.indel.filter.vcf > $target_ID.family.mutation_candidates.indel.filterPASSED.vcf
grep -E '^#|PASS' $target_ID.family.mutation_candidates.indel.filter.vcf > $target_ID.family.indel.filterPASSED.vcf

#COMMENTOUT1
#---------------------------------------------------------------------------------------------------------------
#all sites
#VariantFiltration for SNP
<< COMMENTOUT2
gatk VariantFiltration \
-R $reference_folder/TAIR10.fa \
-V $target_ID.family.snp.vcf.gz \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0" \
--filter-name "basic_snp_filter" \
-O $target_ID.family.snp.filter.vcf

grep -E '^#|PASS' $target_ID.family.snp.filter.vcf > $target_ID.family.snp.filterPASSED.vcf

#VariantFiltration for INDEL
gatk VariantFiltration \
-R $reference_folder/TAIR10.fa \
-V $target_ID.family.indel.vcf.gz \
--filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" \
--filter-name "basic_indel_filter" \
-O $target_ID.family.indel.filter.vcf

grep -E '^#|PASS' $target_ID.family.indel.filter.vcf > $target_ID.family.indel.filterPASSED.vcf

COMMENTOUT2
#=====================================================================


#=====================================================================
<< COMMENTOUT3

#Only medelian violation site only
#DepthFiltering for SNP: DP < 10 & DP > 200 & GQ <20 (Low Genotype Quality: less than 99%)
gatk VariantFiltration \
-R $reference_folder/TAIR10.fa \
-V $target_ID.family.mutation_candidates.snp.filterPASSED.vcf \
-G-filter "GQ < 20" \
-G-filter-name "lowGQ" \
-G-filter "DP < 10 || DP > 200" \
-G-filter-name "DP_10-200" \
-O $target_ID.family.mutation_candidates.snp.DPfilterPASSED.vcf

#DepthFiltering for INDEL: DP < 10 & DP > 200 & GQ < 20 (Low Genotype Quality: less than 99%) 
gatk VariantFiltration \
-R $reference_folder/TAIR10.fa \
-V $target_ID.family.mutation_candidates.indel.filterPASSED.vcf \
-G-filter "GQ < 20" \
-G-filter-name "lowGQ" \
-G-filter "DP < 10 || DP > 200" \
-G-filter-name "DP_10-200" \
-O $target_ID.family.mutation_candidates.indel.DPfilterPASSED.vcf

#Set filtered sites to no call:SNP
gatk SelectVariants \
-R $reference_folder/TAIR10.fa \
-V $target_ID.family.mutation_candidates.snp.DPfilterPASSED.vcf \
--set-filtered-gt-to-nocall \
-O $target_ID.family.mutation_candidates.snp.DPfilterNoCall.vcf
bgzip -c $target_ID.family.mutation_candidates.snp.DPfilterNoCall.vcf > $target_ID.family.mutation_candidates.snp.DPfilterNoCall.vcf.gz
tabix -f -p vcf $target_ID.family.mutation_candidates.snp.DPfilterNoCall.vcf.gz
bcftools index -f $target_ID.family.mutation_candidates.snp.DPfilterNoCall.vcf.gz

#Set filtered sites to no call:INDEL
gatk SelectVariants \
-R $reference_folder/TAIR10.fa \
-V $target_ID.family.mutation_candidates.indel.DPfilterPASSED.vcf \
--set-filtered-gt-to-nocall \
-O $target_ID.family.mutation_candidates.indel.DPfilterNoCall.vcf
bgzip -c $target_ID.family.mutation_candidates.indel.DPfilterNoCall.vcf > $target_ID.family.mutation_candidates.indel.DPfilterNoCall.vcf.gz
tabix -f -p vcf $target_ID.family.mutation_candidates.indel.DPfilterNoCall.vcf.gz
bcftools index -f $target_ID.family.mutation_candidates.indel.DPfilterNoCall.vcf.gz

COMMENTOUT3

#---------------------------------------------------------------------------------------------------------------
#all sites
#DepthFiltering for SNP: DP < 10 & DP > 200 & GQ <20 (Low Genotype Quality: less than 99%)
gatk VariantFiltration \
-R $reference_folder/TAIR10.fa \
-V $target_ID.family.snp.filterPASSED.vcf \
-G-filter "GQ < 20" \
-G-filter-name "lowGQ" \
-G-filter "DP < 10 || DP > 200" \
-G-filter-name "DP_10-200" \
-O $target_ID.family.snp.DPfilterPASSED.vcf

#DepthFiltering for INDEL: DP < 10 & DP > 200 & GQ < 20 (Low Genotype Quality: less than 99%) 
gatk VariantFiltration \
-R $reference_folder/TAIR10.fa \
-V $target_ID.family.indel.filterPASSED.vcf \
-G-filter "GQ < 20" \
-G-filter-name "lowGQ" \
-G-filter "DP < 10 || DP > 200" \
-G-filter-name "DP_10-200" \
-O $target_ID.family.indel.DPfilterPASSED.vcf

#Set filtered sites to no call:SNP
gatk SelectVariants \
-R $reference_folder/TAIR10.fa \
-V $target_ID.family.snp.DPfilterPASSED.vcf \
--set-filtered-gt-to-nocall \
-O $target_ID.family.snp.DPfilterNoCall.vcf
bgzip -c $target_ID.family.snp.DPfilterNoCall.vcf > $target_ID.family.snp.DPfilterNoCall.vcf.gz
tabix -f -p vcf $target_ID.family.snp.DPfilterNoCall.vcf.gz
bcftools index -f $target_ID.family.snp.DPfilterNoCall.vcf.gz

#Set filtered sites to no call:INDEL
gatk SelectVariants \
-R $reference_folder/TAIR10.fa \
-V $target_ID.family.indel.DPfilterPASSED.vcf \
--set-filtered-gt-to-nocall \
-O $target_ID.family.indel.DPfilterNoCall.vcf
bgzip -c $target_ID.family.indel.DPfilterNoCall.vcf > $target_ID.family.indel.DPfilterNoCall.vcf.gz
tabix -f -p vcf $target_ID.family.indel.DPfilterNoCall.vcf.gz
bcftools index -f $target_ID.family.indel.DPfilterNoCall.vcf.gz

#=====================================================================


cd $SCRIPT_DIR

module unload samtools/1.10
module unload gatk/4.1.7.0


