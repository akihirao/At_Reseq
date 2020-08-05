#!/bin/bash -i
#Run.MutationIdentification.pindel_vcf.sh
#by HIRAO Akira

#requirement:
#*FilteringVcfNeighborSNVs.pl: filtering out combined variants (SNPs and INDELs) around neighborhood
#*VariantFilteredAF.pl: filtering out mutation sites whrere proportion of mutant reads < 25% && > 80% && GQ < 99
#*MakeMulationList.pl: output mutation list

set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)

module load gatk/4.1.7.0
module load vcftools/0.1.15
module load bedops/2.4.39


target_ID=AT48

reference_folder=/zfs/Arabidopsis/Reference_v1.1
main_folder=/zfs/Arabidopsis/work/At_Reseq
work_folder=$main_folder/pindel_out

BioAlcidaeJdk_path=/usr/local/jvarkit/dist

#mutation_summary_file="mutation_summary"
#mutation_summary_file_homo="mutation_summary.homo"
#mutation_list_file="M2.mutation_list"

#echo -n >| $mutation_summary_file.txt
#echo -n >| $mutation_summary_file_homo.txt
#echo -n >| $mutation_list_file.txt


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

#six types of structural variants（D = deletion、SI = short insertion、INV = inversion、TD = tandem duplication、LI = large insertion、BP = unassigned breakpoints）
bgzip -c $target_ID.pindel_out.D.vcf > $target_ID.pindel_out.D.vcf.gz
tabix -f -p vcf $target_ID.pindel_out.D.vcf.gz

bgzip -c $target_ID.pindel_out.SI.vcf > $target_ID.pindel_out.SI.vcf.gz
tabix -f -p vcf $target_ID.pindel_out.SI.vcf.gz


#filtering out neighborhood mutations 
perl $SCRIPT_DIR/FilteringVcfNeighborSNVs.pl < $target_ID.pindel_out.D.vcf > $target_ID.pindel_out.D.non_neighbor.vcf
perl $SCRIPT_DIR/FilteringVcfNeighborSNVs.pl < $target_ID.pindel_out.SI.vcf > $target_ID.pindel_out.SI.non_neighbor.vcf

#identifying unique SNVs with using bioalcidaejdk.jar
#See detail for https://www.biostars.org/p/329423/#329742
java -jar $BioAlcidaeJdk_path/bioalcidaejdk.jar -e 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0);println(V.getContig()+" "+V.getStart()+" "+V.getReference()+" "+g.getSampleName()+" "+g.getAlleles());});' $target_ID.pindel_out.D.non_neighbor.vcf > $target_ID.unique.pindel.D.list.txt
java -jar $BioAlcidaeJdk_path/bioalcidaejdk.jar -e 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0);println(V.getContig()+" "+V.getStart()+" "+V.getReference()+" "+g.getSampleName()+" "+g.getAlleles());});' $target_ID.pindel_out.SI.non_neighbor.vcf > $target_ID.unique.pindel.SI.list.txt


perl $SCRIPT_DIR/BioalcidaejdkOut2BED.pl < $target_ID.unique.pindel.D.list.txt > $target_ID.unique.pindel.D.bed
perl $SCRIPT_DIR/BioalcidaejdkOut2BED.pl < $target_ID.unique.pindel.SI.list.txt > $target_ID.unique.pindel.SI.bed
cat $target_ID.unique.pindel.D.bed $target_ID.unique.pindel.SI.bed | sort -k 1,1 -k 2n,2 > $target_ID.unique.pindel.indel.bed

#select back variants with unique.bed 
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.pindel_out.D.vcf.gz\
 -L $target_ID.unique.pindel.D.bed\
 --exclude-sample-name $SCRIPT_DIR/Mother_ID.list\
 --exclude-non-variants\
 -O AT.M2.unique.pindel.D.vcf
bgzip -c AT.M2.unique.pindel.D.vcf > AT.M2.unique.pindel.D.vcf.gz
tabix -f -p vcf AT.M2.unique.pindel.D.vcf.gz

gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.pindel_out.SI.vcf.gz\
 -L $target_ID.unique.pindel.SI.bed\
 --exclude-sample-name $SCRIPT_DIR/Mother_ID.list\
 --exclude-non-variants\
 -O AT.M2.unique.pindel.SI.vcf
bgzip -c AT.M2.unique.pindel.SI.vcf > AT.M2.unique.pindel.SI.vcf.gz
tabix -f -p vcf AT.M2.unique.pindel.SI.vcf.gz

vcf-concat AT.M2.unique.pindel.D.vcf.gz AT.M2.unique.pindel.SI.vcf.gz > AT.M2.unique.pindel.indel.unsorted.vcf
vcf-sort AT.M2.unique.pindel.indel.unsorted.vcf > AT.M2.unique.pindel.indel.vcf
bgzip -c AT.M2.unique.pindel.indel.vcf > AT.M2.unique.pindel.indel.vcf.gz
tabix -f -p vcf AT.M2.unique.pindel.indel.vcf.gz
bcftools index AT.M2.unique.pindel.indel.vcf.gz

#compare indel variants identified by gatk with those by ppindel
bcftools isec ../vcf_out/M2.indel.unique.vcf.gz AT.M2.unique.pindel.D.vcf.gz -p AT.M2.gatk.pindel.common -n=2


cd $SCRIPT_DIR

module unload gatk/4.1.7.0
module unload vcftools/0.1.15
module unload bedops/2.4.39
