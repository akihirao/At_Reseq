#!/bin/bash -i
#Pipe.09.VariantFiltration.sh
#by HIRAO Akira

#modified to restrict to biallelic @2020530
#modified for filtering out ExcessHet P < 0.05

set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)

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

#=====================================================================
#<< COMMENTOUT1
#Only medelian violation site only
#VariantFiltration for SNP
#/Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java -jar /usr/local/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar VariantFiltration\
gatk VariantFiltration\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.mendelian.snp.vcf.gz\
 --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0"\
 --filter-name "basic_snp_filter"\
 --filter-expression "${ExcessHet_param}"\
 --filter-name "ExHet"\
 -O $target_ID.mendelian.snp.filter.vcf.gz

#grep -E '^#|PASS' $target_ID.mendelian.snp.filter.vcf > $target_ID.mendelian.snp.filterPASSED.vcf
gzcat $target_ID.mendelian.snp.filter.vcf.gz | grep -E '^#|PASS' |bgzip > $target_ID.mendelian.snp.filterPASSED.vcf.gz
tabix -f -p vcf $target_ID.mendelian.snp.filterPASSED.vcf.gz

#VariantFiltration for INDEL
#/Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java -jar /usr/local/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar VariantFiltration\
gatk VariantFiltration\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.mendelian.indel.vcf.gz \
 --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0"\
 --filter-name "basic_indel_filter"\
 --filter-expression "${ExcessHet_param}"\
 --filter-name "ExHet"\
 -O $target_ID.mendelian.indel.filter.vcf.gz

#grep -E '^#|PASS' $target_ID.mendelian.indel.filter.vcf > $target_ID.mendelian.indel.filterPASSED.vcf
gzcat $target_ID.mendelian.indel.filter.vcf.gz | grep -E '^#|PASS' | bgzip > $target_ID.mendelian.indel.filterPASSED.vcf.gz
tabix -f -p vcf $target_ID.mendelian.indel.filterPASSED.vcf.gz

#COMMENTOUT1
#---------------------------------------------------------------------------------------------------------------
#all sites
#VariantFiltration for SNP
#<< COMMENTOUT2
#/Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java -jar /usr/local/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar VariantFiltration\
gatk VariantFiltration\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.snp.vcf.gz \
 --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0"\
 --filter-name "basic_snp_filter"\
 --filter-expression "${ExcessHet_param}"\
 --filter-name "ExHet"\
 -O $target_ID.snp.filter.vcf.gz

#grep -E '^#|PASS' $target_ID.snp.filter.vcf > $target_ID.snp.filterPASSED.vcf
gzcat $target_ID.snp.filter.vcf.gz | grep -E '^#|PASS' | bgzip > $target_ID.snp.filterPASSED.vcf.gz
tabix -f -p vcf $target_ID.snp.filterPASSED.vcf.gz


#VariantFiltration for INDEL
#/Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java -jar /usr/local/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar VariantFiltration\
gatk VariantFiltration\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.indel.vcf.gz\
 --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0"\
 --filter-name "basic_indel_filter"\
 --filter-expression "${ExcessHet_param}"\
 --filter-name "ExHet"\
 -O $target_ID.indel.filter.vcf.gz

#grep -E '^#|PASS' $target_ID.indel.filter.vcf > $target_ID.indel.filterPASSED.vcf
gzcat $target_ID.indel.filter.vcf.gz | grep -E '^#|PASS' | bgzip > $target_ID.indel.filterPASSED.vcf.gz
tabix -f -p vcf $target_ID.indel.filterPASSED.vcf.gz

#COMMENTOUT2
#=====================================================================


#=====================================================================
#<< COMMENTOUT3

#Only medelian violation site only
#DepthFiltering for SNP: DP < 10 & DP > 200 & GQ <20 (Low Genotype Quality: less than 99%)
#/Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java -jar /usr/local/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar VariantFiltration\
gatk VariantFiltration\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.mendelian.snp.filterPASSED.vcf.gz\
 -G-filter "GQ < 20"\
 -G-filter-name "lowGQ"\
 -G-filter "DP < 10 || DP > 200"\
 -G-filter-name "DP_10-200"\
 -O $target_ID.mendelian.snp.DPfilterPASSED.vcf.gz

#DepthFiltering for INDEL: DP < 10 & DP > 200 & GQ < 20 (Low Genotype Quality: less than 99%) 
#/Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java -jar /usr/local/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar VariantFiltration\
gatk VariantFiltration\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.mendelian.indel.filterPASSED.vcf.gz\
 -G-filter "GQ < 20"\
 -G-filter-name "lowGQ"\
 -G-filter "DP < 10 || DP > 200"\
 -G-filter-name "DP_10-200"\
 -O $target_ID.mendelian.indel.DPfilterPASSED.vcf.gz

#Set filtered sites to no call:SNP
#--restrict-alleles-to BIALLELIC
#/Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java -jar /usr/local/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar SelectVariants\
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.mendelian.snp.DPfilterPASSED.vcf.gz\
 --restrict-alleles-to BIALLELIC\
 --set-filtered-gt-to-nocall\
 -O $target_ID.mendelian.snp.DPfilterNoCall.vcf.gz
#bgzip -c $target_ID.mendelian.snp.DPfilterNoCall.vcf > $target_ID.mendelian.snp.DPfilterNoCall.vcf.gz
#tabix -f -p vcf $target_ID.mendelian.snp.DPfilterNoCall.vcf.gz
bcftools index -f $target_ID.mendelian.snp.DPfilterNoCall.vcf.gz
bcftools view $target_ID.mendelian.snp.DPfilterNoCall.vcf.gz -Ov -o $target_ID.mendelian.snp.DPfilterNoCall.vcf

#Set filtered sites to no call:INDEL
#/Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java -jar /usr/local/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar SelectVariants\
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.mendelian.indel.DPfilterPASSED.vcf.gz\
 --restrict-alleles-to BIALLELIC\
 --set-filtered-gt-to-nocall\
 -O $target_ID.mendelian.indel.DPfilterNoCall.vcf.gz
#bgzip -c $target_ID.mendelian.indel.DPfilterNoCall.vcf > $target_ID.mendelian.indel.DPfilterNoCall.vcf.gz
#tabix -f -p vcf $target_ID.mendelian.indel.DPfilterNoCall.vcf.gz
bcftools index -f $target_ID.mendelian.indel.DPfilterNoCall.vcf.gz
bcftools view $target_ID.mendelian.indel.DPfilterNoCall.vcf.gz -Ov -o $target_ID.mendelian.indel.DPfilterNoCall.vcf

#COMMENTOUT3

#---------------------------------------------------------------------------------------------------------------
#all sites
#DepthFiltering for SNP: DP < 10 & DP > 200 & GQ <20 (Low Genotype Quality: less than 99%)
#/Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java -jar /usr/local/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar VariantFiltration\
gatk VariantFiltration\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.snp.filterPASSED.vcf.gz\
 -G-filter "GQ < 20"\
 -G-filter-name "lowGQ"\
 -G-filter "DP < 10 || DP > 200"\
 -G-filter-name "DP_10-200"\
 -O $target_ID.snp.DPfilterPASSED.vcf.gz

#DepthFiltering for INDEL: DP < 10 & DP > 200 & GQ < 20 (Low Genotype Quality: less than 99%) 
#/Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java -jar /usr/local/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar VariantFiltration\
gatk VariantFiltration\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.indel.filterPASSED.vcf.gz\
 -G-filter "GQ < 20"\
 -G-filter-name "lowGQ"\
 -G-filter "DP < 10 || DP > 200"\
 -G-filter-name "DP_10-200"\
 -O $target_ID.indel.DPfilterPASSED.vcf.gz

#Set filtered sites to no call:SNP
#/Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java -jar /usr/local/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar SelectVariants\
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.snp.DPfilterPASSED.vcf.gz\
 --restrict-alleles-to BIALLELIC\
 --set-filtered-gt-to-nocall\
 -O $target_ID.snp.DPfilterNoCall.vcf.gz
#bgzip -c $target_ID.snp.DPfilterNoCall.vcf > $target_ID.snp.DPfilterNoCall.vcf.gz
#tabix -f -p vcf $target_ID.snp.DPfilterNoCall.vcf.gz
bcftools index -f $target_ID.snp.DPfilterNoCall.vcf.gz
bcftools view $target_ID.snp.DPfilterNoCall.vcf.gz -Ov -o $target_ID.snp.DPfilterNoCall.vcf

#Set filtered sites to no call:INDEL
#/Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java -jar /usr/local/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar SelectVariants\
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.indel.DPfilterPASSED.vcf.gz\
 --restrict-alleles-to BIALLELIC\
 --set-filtered-gt-to-nocall\
 -O $target_ID.indel.DPfilterNoCall.vcf.gz
#bgzip -c $target_ID.indel.DPfilterNoCall.vcf > $target_ID.indel.DPfilterNoCall.vcf.gz
#tabix -f -p vcf $target_ID.indel.DPfilterNoCall.vcf.gz
bcftools index -f $target_ID.indel.DPfilterNoCall.vcf.gz
bcftools view $target_ID.indel.DPfilterNoCall.vcf.gz -Ov -o $target_ID.indel.DPfilterNoCall.vcf
#=====================================================================


cd $SCRIPT_DIR

module unload samtools/1.10
module unload gatk/4.1.7.0


