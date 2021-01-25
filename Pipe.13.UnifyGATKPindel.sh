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


#compare indel variants identified by gatk with those by pindel: hetero
#bcftools isec $vcf_folder/M2.indel.hetero.vcf.gz AT.M2.unique.pindel.D.SI.vcf.gz -p AT.M2.gatk.pindel.common.hetero -n=2
#perl $SCRIPT_DIR/Vcf2BED_chr_start_end.pl < $work_folder/AT.M2.gatk.pindel.common.hetero/0000.vcf > $work_folder/AT.M2.indel.common.pindel.gatk.hetero.bed

#compare indel variants identified by gatk with those by pindel: homo
#bcftools isec $vcf_folder/M2.indel.homo.vcf.gz AT.M2.unique.pindel.D.SI.vcf.gz -p AT.M2.gatk.pindel.common.homo -n=2
#perl $SCRIPT_DIR/Vcf2BED_chr_start_end.pl < $work_folder/AT.M2.gatk.pindel.common.homo/0000.vcf > $work_folder/AT.M2.indel.common.pindel.gatk.homo.bed

#compare indel variants identified by gatk with those by pindel: family-clustered
#bcftools isec $vcf_folder/M2.indel.familyclustered.vcf.gz AT.M2.unique.pindel.D.SI.vcf.gz -p AT.M2.gatk.pindel.common.familyclustered -n=2
#perl $SCRIPT_DIR/Vcf2BED_chr_start_end.pl < $work_folder/AT.M2.gatk.pindel.common.familyclustered/0000.vcf > $work_folder/AT.M2.indel.common.pindel.gatk.familyclustered.bed

#compare indel variants identified by gatk with those by pindel: all
bcftools isec $vcf_folder/M2.indel.all.vcf.gz AT.M2.unique.pindel.D.SI.vcf.gz -p AT.M2.gatk.pindel.common.all -n=2
perl $SCRIPT_DIR/Vcf2BED_chr_start_end.pl < $work_folder/AT.M2.gatk.pindel.common.all/0000.vcf > $work_folder/AT.M2.indel.common.pindel.gatk.all.bed


#at $vcf_folder/M2.snp.hetero.bed $work_folder/AT.M2.indel.common.pindel.gatk.hetero.bed | sort -k 1,1 -k 2n,2 >  $work_folder/M2.snp.indel.gatk.pindl.common.hetero.bed
#cat $vcf_folder/M2.snp.homo.bed $work_folder/AT.M2.indel.common.pindel.gatk.homo.bed | sort -k 1,1 -k 2n,2 >  $work_folder/M2.snp.indel.gatk.pindl.common.homo.bed
#cat $vcf_folder/M2.snp.hetero.bed $work_folder/AT.M2.indel.common.pindel.gatk.familyclustered.bed | sort -k 1,1 -k 2n,2 >  $work_folder/M2.snp.indel.gatk.pindl.common.familyclustered.bed

cat $vcf_folder/M2.snp.all.bed $work_folder/AT.M2.indel.common.pindel.gatk.all.bed | sort -k 1,1 -k 2n,2 >  $work_folder/M2.snp.indel.gatk.pindl.common.all.bed


#select back variants with hetero.bed 
#gatk SelectVariants\
# -R $reference_folder/TAIR10.fa\
# -V $vcf_folder/$target_ID.mu.snp.indel.DPfilterNoCall.vcf.gz\
# -L $work_folder/M2.snp.indel.gatk.pindl.common.hetero.bed\
# --exclude-sample-name $SCRIPT_DIR/Mother_ID.list\
# -O $vcf_folder/M2.snp.indel.gatk.pindl.common.hetero.vcf

#select back variants with homo.bed 
#gatk SelectVariants\
# -R $reference_folder/TAIR10.fa\
# -V $vcf_folder/$target_ID.mu.snp.indel.DPfilterNoCall.vcf.gz\
# -L $work_folder/M2.snp.indel.gatk.pindl.common.homo.bed\
# --exclude-sample-name $SCRIPT_DIR/Mother_ID.list\
# -O $vcf_folder/M2.snp.indel.gatk.pindl.common.homo.vcf

#select back variants with familyclustered.bed 
#gatk SelectVariants\
# -R $reference_folder/TAIR10.fa\
# -V $vcf_folder/$target_ID.mu.snp.indel.DPfilterNoCall.vcf.gz\
# -L $work_folder/M2.snp.indel.gatk.pindl.common.familyclustered.bed\
# --exclude-sample-name $SCRIPT_DIR/Mother_ID.list\
# -O $vcf_folder/M2.snp.indel.gatk.pindl.common.familyclustered.vcf

#select back variants with all.bed 
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $vcf_folder/$target_ID.mu.snp.indel.DPfilterNoCall.vcf.gz\
 -L $work_folder/M2.snp.indel.gatk.pindl.common.all.bed\
 --exclude-sample-name $SCRIPT_DIR/Mother_ID.list\
 -O $vcf_folder/M2.snp.indel.gatk.pindl.common.hetero.homo.familyclustered.vcf
bgzip -c $vcf_folder/M2.snp.indel.gatk.pindl.common.hetero.homo.familyclustered.vcf > $vcf_folder/M2.snp.indel.gatk.pindl.common.hetero.homo.familyclustered.vcf.gz
tabix -f -p vcf $vcf_folder/M2.snp.indel.gatk.pindl.common.hetero.homo.familyclustered.vcf.gz
	

#java -jar $BioAlcidaeJdk_path/bioalcidaejdk.jar -e 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0);println(V.getContig()+" "+V.getStart()+" "+V.getReference()+" "+g.getSampleName()+" "+g.getAlleles());});' $vcf_folder/M2.snp.indel.gatk.pindl.common.hetero.vcf > $vcf_folder/M2.snp.indel.gatk.pindl.common.hetero.list.txt
#java -jar $BioAlcidaeJdk_path/bioalcidaejdk.jar -e 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0);println(V.getContig()+" "+V.getStart()+" "+V.getReference()+" "+g.getSampleName()+" "+g.getAlleles());});' $vcf_folder/M2.snp.indel.gatk.pindl.common.all.vcf > $vcf_folder/M2.snp.indel.gatk.pindl.common.hetero.homo.list.txt


#perl $SCRIPT_DIR/BioalcidaejdkList2Rform.pl < $vcf_folder/M2.snp.indel.gatk.pindl.common.hetero.list.txt > $vcf_folder/M2.snp.indel.gatk.pindl.common.hetero.R.list.txt
#perl $SCRIPT_DIR/BioalcidaejdkList2Rform.pl < $vcf_folder/M2.snp.indel.gatk.pindl.common.hetero.homo.list.txt > $vcf_folder/M2.snp.indel.gatk.pindl.common.hetero.homo.R.list.txt


#Merge homo + hetero vcfs + familyclustered vcf
#gatk MergeVcfs\
# -I $vcf_folder/M2.snp.indel.gatk.pindl.common.hetero.vcf\
# -I $vcf_folder/M2.snp.indel.gatk.pindl.common.homo.vcf\
# -I $vcf_folder/AT.M2.family.clustered.mu.non_neighbor.vcf.gz\
# -O $vcf_folder/M2.snp.indel.gatk.pindl.common.homo.hetero.familyclustered.vcf.gz
#bcftools index -f $vcf_folder/M2.snp.indel.gatk.pindl.common.homo.hetero.familyclustered.vcf.gz
#bcftools view  $vcf_folder/M2.snp.indel.gatk.pindl.common.homo.hetero.familyclustered.vcf.gz -Ov -o $vcf_folder/M2.snp.indel.gatk.pindl.common.homo.hetero.familyclustered.vcf


perl $SCRIPT_DIR/Vcf2List_R_glm.pl < $vcf_folder/M2.snp.indel.gatk.pindl.common.hetero.homo.familyclustered.vcf > $SCRIPT_DIR/M2.mutations.full.list.csv

cd $SCRIPT_DIR

	
module unload gatk/4.1.7.0
module unload vcftools/0.1.15

