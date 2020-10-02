#!/bin/bash -i
#Pipe.12.UnifyGATKPindel.sh
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
work_folder=$main_folder/pindel_out

BioAlcidaeJdk_path=/usr/local/jvarkit/dist

cd $work_folder

#compare indel variants identified by gatk with those by pindel
bcftools isec $vcf_folder/M2.indel.unique.vcf.gz AT.M2.unique.pindel.indel.vcf.gz -p AT.M2.gatk.pindel.common -n=2
perl $SCRIPT_DIR/Vcf2BED_chr_start_end.pl < $work_folder/AT.M2.gatk.pindel.common/0000.vcf > $work_folder/AT.M2.indel.common.pindel.gatk.bed


bcftools isec $vcf_folder/M2.indel.unique.vcf.gz AT.M2.unique.pindel.indel.vcf.gz -p AT.M2.gatk.unique -C
perl $SCRIPT_DIR/Vcf2BED_chr_start_end.pl < $work_folder/AT.M2.gatk.unique/0000.vcf > $work_folder/AT.M2.filterout.pindel.gatk.indel.bed

cat $vcf_folder/M2.snp.unique.bed $work_folder/AT.M2.indel.common.pindel.gatk.bed | sort -k 1,1 -k 2n,2 >  $work_folder/M2.snp.indel.unique.gatk.pindl.common.bed

#select back variants with unique.bed 
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $vcf_folder/$target_ID.mu.snp.indel.DPfilterNoCall.vcf.gz\
 -L $work_folder/M2.snp.indel.unique.gatk.pindl.common.bed\
 --exclude-sample-name $SCRIPT_DIR/Mother_ID.list\
 -O $vcf_folder/M2.snp.indel.unique.gatk.pindl.common.vcf

java -jar $BioAlcidaeJdk_path/bioalcidaejdk.jar -e 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0);println(V.getContig()+" "+V.getStart()+" "+V.getReference()+" "+g.getSampleName()+" "+g.getAlleles());});' $vcf_folder/M2.snp.indel.unique.gatk.pindl.common.vcf > $vcf_folder/M2.snp.indel.unique.gatk.pindl.common.list.txt

perl $SCRIPT_DIR/BioalcidaejdkList2Rform.pl < $vcf_folder/M2.snp.indel.unique.gatk.pindl.common.list.txt > $vcf_folder/M2.snp.indel.unique.gatk.pindl.common.R.list.txt

cd $SCRIPT_DIR

	
module unload gatk/4.1.7.0
module unload vcftools/0.1.15
