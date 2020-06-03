#!/bin/bash -i
#Pipe.10.MutationIdentification.sh
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
work_folder=$main_folder/vcf_out

BioAlcidaeJdk_path=/usr/local/jvarkit/dist

mutation_summary_file="mutation_summary.txt"
mutation_summary_file_homo="mutation_summary.homo.txt"
mutation_list_file="M2.mutation_list"

echo -n >| $work_folder/$mutation_summary_file
echo -n >| $work_folder/$mutation_summary_file_homo
echo -n >| $work_folder/$mutation_list_file.txt


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

#----------------------------------------------------------------------------------
#Output mendelian violation site only: --pedigree & --mendelian-violation
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.snp.indel.DPfilterNoCall.vcf.gz\
 --pedigree $SCRIPT_DIR/AT48.ped\
 --mendelian-violation\
 --mendelian-violation-qual-threshold 30\
 -O $work_folder/$target_ID.mu.snp.indel.DPfilterNoCall.vcf.gz
bcftools index -f $work_folder/$target_ID.mu.snp.indel.DPfilterNoCall.vcf.gz
bcftools view $work_folder/$target_ID.mu.snp.indel.DPfilterNoCall.vcf.gz\
 -Ov -o $work_folder/$target_ID.mu.snp.indel.DPfilterNoCall.vcf


#identifying unique SNVs with using bioalcidaejdk.jar
#See detail for https://www.biostars.org/p/329423/#329742
java -jar $BioAlcidaeJdk_path/bioalcidaejdk.jar -e 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0);println(V.getContig()+" "+V.getStart()+" "+V.getReference()+" "+g.getSampleName()+" "+g.getAlleles());});' $work_folder/$target_ID.mu.snp.indel.DPfilterNoCall.vcf > $work_folder/$target_ID.unique.snp.indel.list.txt

perl $SCRIPT_DIR/BioalcidaejdkOut2BED.pl < $work_folder/$target_ID.unique.snp.indel.list.txt > $work_folder/$target_ID.unique.bed

gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $work_folder/$target_ID.mu.snp.indel.DPfilterNoCall.vcf.gz\
 -L $work_folder/$target_ID.unique.bed\
 --exclude-sample-name $SCRIPT_DIR/Mother_ID.list\
 -O $work_folder/AT.M2.unique.vcf.gz

#mark homozygous mutation sites
gatk VariantFiltration\
 -R $reference_folder/TAIR10.fa\
 -V $work_folder/AT.M2.unique.vcf.gz\
 -G-filter "isHomVar==1"\
 -G-filter-name "homozygous_mutation"\
 -G-filter "isHomRef==1"\
 -G-filter-name "homozygous_ref"\
 -O $work_folder/AT.M2.unique.hetero_marked.vcf.gz

#for target_sample in ${M2_vec[@]}
for target_sample in A011 A012 A014
do

	mkdir -p $target_sample

	gatk SelectVariants\
	 -R $reference_folder/TAIR10.fa\
	 -V $work_folder/AT.M2.unique.hetero_marked.vcf.gz\
	 --set-filtered-gt-to-nocall\
	 --sample-name $target_sample\
	 -O $work_folder/$target_sample/$target_sample.hetero.filter.vcf.gz

	gatk SelectVariants\
	 -R $reference_folder/TAIR10.fa\
	 -V $work_folder/$target_sample/$target_sample.hetero.filter.vcf.gz\
	 --max-nocall-number 0\
	 -O $work_folder/$target_sample/$target_sample.hetero.vcf

	bgzip -c $work_folder/$target_sample/$target_sample.hetero.vcf > $work_folder/$target_sample/$target_sample.hetero.vcf.gz
	tabix -f -p vcf $work_folder/$target_sample/$target_sample.hetero.vcf.gz

	#select homozygous mutation sites per samples
	gatk SelectVariants\
	 -R $reference_folder/TAIR10.fa\
	 -V $work_folder/AT.M2.unique.vcf.gz\
	 -select 'vc.getGenotype("${target_sample}").isHomRef() || vc.getGenotype("${target_sample}").isHomVar()'\
	 --sample-name $target_sample\
	 -O $work_folder/$target_sample/$target_sample.homo.vcf
	bgzip -c $work_folder/$target_sample/$target_sample.homo.vcf > $work_folder/$target_sample/$target_sample.homo.vcf.gz
	tabix -f -p vcf $work_folder/$target_sample/$target_sample.homo.vcf.gz


	#filtering out combined mutations 
	perl $SCRIPT_DIR/FilteringVcfNeighborSNVs.pl < $work_folder/$target_sample/$target_sample.hetero.vcf > $work_folder/$target_sample/$target_sample.non_neighbor.vcf
	bgzip -c $work_folder/$target_sample/$target_sample.non_neighbor.vcf > $work_folder/$target_sample/$target_sample.non_neighbor.vcf.gz
	tabix -f -p vcf $work_folder/$target_sample/$target_sample.non_neighbor.vcf.gz

	#picking up combined mutations
	vcf-isec -c $work_folder/$target_sample/$target_sample.hetero.vcf.gz $work_folder/$target_sample/$target_sample.non_neighbor.vcf.gz > $work_folder/$target_sample/$target_sample.combined.vcf

	#filtering out mutation sites whrere proportion of mutant reads < 25% or GQ < 99
	perl $SCRIPT_DIR/VariantFilteredAF.pl < $work_folder/$target_sample/$target_sample.non_neighbor.vcf > $work_folder/$target_sample/$target_sample.final.mutants.vcf
	perl $SCRIPT_DIR/VariantFilteredAF.pl < $work_folder/$target_sample/$target_sample.combined.vcf > $work_folder/$target_sample/$target_sample.final.combined.vcf

	#output heterozygous mutation.list 
	perl $SCRIPT_DIR/MakeMulationList.pl $work_folder/$target_sample/$target_sample.final.mutants.vcf $target_sample >> $work_folder/$mutation_list_file.txt
	#output homozygous mutation.list 
	perl $SCRIPT_DIR/MakeMulationList.pl $work_folder/$target_sample/$target_sample.homo.vcf $target_sample >> $work_folder/$mutation_list_file.homo.txt


	vcf2bed --snvs <  $work_folder/$target_sample/$target_sample.final.mutants.vcf > $work_folder/$target_sample/$target_sample.snp.bed
	vcf2bed --insertions < $work_folder/$target_sample/$target_sample.final.mutants.vcf > $work_folder/$target_sample/$target_sample.insertion.bed
	vcf2bed --deletions <  $work_folder/$target_sample/$target_sample.final.mutants.vcf > $work_folder/$target_sample/$target_sample.deletion.bed
	vcf2bed <  $work_folder/$target_sample/$target_sample.final.combined.vcf > $work_folder/$target_sample/$target_sample.final.combined.bed
  
	vcf2bed --snvs <  $work_folder/$target_sample/$target_sample.homo.vcf > $work_folder/$target_sample/$target_sample.snp.homo.bed
	vcf2bed --insertions < $work_folder/$target_sample/$target_sample.homo.vcf > $work_folder/$target_sample/$target_sample.insertion.homo.bed
	vcf2bed --deletions <  $work_folder/$target_sample/$target_sample.homo.vcf > $work_folder/$target_sample/$target_sample.deletion.homo.bed

	no_snp=`wc -l $work_folder/$target_sample/$target_sample.snp.bed |awk '{print $1}'`
	no_insertion=`wc -l $work_folder/$target_sample/$target_sample.insertion.bed |awk '{print $1}'`
	no_deletion=`wc -l $work_folder/$target_sample/$target_sample.deletion.bed |awk '{print $1}'`
	no_indel=$((no_insertion + no_deletion))

	no_snp_homo=`wc -l $work_folder/$target_sample/$target_sample.snp.homo.bed |awk '{print $1}'`
	no_insertion_homo=`wc -l $work_folder/$target_sample/$target_sample.insertion.homo.bed |awk '{print $1}'`
	no_deletion_homo=`wc -l $work_folder/$target_sample/$target_sample.deletion.homo.bed |awk '{print $1}'`
	no_indel_homo=$((no_insertion_homo + no_deletion_homo))
	 
	tab_lab="	"
	output_info=$target_sample$tab_lab$no_snp$tab_lab$no_indel$tab_lab$no_insertion$tab_lab$no_deletion
	echo $output_info >> $work_folder/$mutation_summary_file.txt

	output_info_homo=$target_sample$tab_lab$no_snp_homo$tab_lab$no_indel_homo$tab_lab$no_insertion_homo$tab_lab$no_deletion_homo
	echo $output_info_homo >> $work_folder/$mutation_summary_file_homo.txt

done


#Sort SNPs + INDELs heterozygous mutation list 
cat $work_folder/$mutation_list_file.txt | sort -k 1,1 -k 2,2 >  $work_folder/$mutation_list_file.position_sort.txt
#Extract SNPs heterozygous mutation list
awk '/SNP/' $work_folder/$mutation_list_file.position_sort.txt > $work_folder/M2.snp.list.txt
#Extract INDELs heterozygous mutation list
awk '/Insertion/ || /Deletion/' $work_folder/$mutation_list_file.position_sort.txt > $work_folder/M2.indel.list.txt

#Sort SNPs + INDELs homozygous mutation list 
cat $work_folder/$mutation_list_file.homo.txt | sort -k 1,1 -k 2,2 >  $work_folder/$mutation_list_file.position_sort.homo.txt
#Extract SNPs homozygous mutation list
awk '/SNP/' $work_folder/$mutation_list_file.position_sort.homo.txt > $work_folder/M2.snp.homo.list.txt
#Extract INDELs homozygous mutation list
awk '/Insertion/ || /Deletion/' $work_folder/$mutation_list_file.position_sort.homo.txt > $work_folder/M2.indel.homo.list.txt




cd $SCRIPT_DIR

module unload gatk/4.1.7.0
module unload vcftools/0.1.15
module unload bedops/2.4.39
