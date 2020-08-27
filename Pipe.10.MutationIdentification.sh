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

mutation_summary_file="mutation_summary.gatk"
mutation_summary_file_homo="mutation_summary.gatk.homo"
mutation_list_file="M2.mutation_list.gatk"

echo -n >| $mutation_summary_file.txt
echo -n >| $mutation_summary_file_homo.txt
echo -n >| $mutation_list_file.txt


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
 -O $target_ID.mu.snp.indel.DPfilterNoCall.vcf.gz
bcftools index -f $target_ID.mu.snp.indel.DPfilterNoCall.vcf.gz
bcftools view $target_ID.mu.snp.indel.DPfilterNoCall.vcf.gz\
 -Ov -o $target_ID.mu.snp.indel.DPfilterNoCall.vcf


#identifying unique SNVs with using bioalcidaejdk.jar
#See detail for https://www.biostars.org/p/329423/#329742
java -jar $BioAlcidaeJdk_path/bioalcidaejdk.jar -e 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0);println(V.getContig()+" "+V.getStart()+" "+V.getReference()+" "+g.getSampleName()+" "+g.getAlleles());});' $target_ID.mu.snp.indel.DPfilterNoCall.vcf > $target_ID.unique.snp.indel.list.txt

perl $SCRIPT_DIR/BioalcidaejdkOut2BED.pl < $target_ID.unique.snp.indel.list.txt > $target_ID.unique.bed

#select back variants with unique.bed 
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.mu.snp.indel.DPfilterNoCall.vcf.gz\
 -L $target_ID.unique.bed\
 --exclude-sample-name $SCRIPT_DIR/Mother_ID.list\
 --max-nocall-fraction 0.1\
 --exclude-filtered\
 -O AT.M2.unique.vcf.gz

#mark homozygous mutation sites
gatk VariantFiltration\
 -R $reference_folder/TAIR10.fa\
 -V AT.M2.unique.vcf.gz\
 -G-filter "isHomVar==1"\
 -G-filter-name "homozygous_mutation"\
 -G-filter "isHomRef==1"\
 -G-filter-name "homozygous_ref"\
 -O AT.M2.unique.hetero_marked.vcf.gz


for target_sample in ${M2_vec[@]}
do

	mkdir -p $target_sample

	gatk SelectVariants\
	 -R $reference_folder/TAIR10.fa\
	 -V AT.M2.unique.hetero_marked.vcf.gz\
	 --set-filtered-gt-to-nocall\
	 --sample-name $target_sample\
	 -O $target_sample/$target_sample.hetero.filter.vcf.gz

	gatk SelectVariants\
	 -R $reference_folder/TAIR10.fa\
	 -V $target_sample/$target_sample.hetero.filter.vcf.gz\
	 --max-nocall-number 0\
	 -O $target_sample/$target_sample.hetero.vcf

	bgzip -c $target_sample/$target_sample.hetero.vcf > $target_sample/$target_sample.hetero.vcf.gz
	tabix -f -p vcf $target_sample/$target_sample.hetero.vcf.gz

	#select homozygous mutation sites per samples with filtering out QG < 99
	gatk SelectVariants\
	 -R $reference_folder/TAIR10.fa\
	 -V AT.M2.unique.vcf.gz\
	 -select "vc.getGenotype('${target_sample}').isHomVar() && vc.getGenotype('${target_sample}').getGQ() == 99 "\
	 --sample-name $target_sample\
	 -O $target_sample/$target_sample.homo.vcf
	bgzip -c $target_sample/$target_sample.homo.vcf > $target_sample/$target_sample.homo.vcf.gz
	tabix -f -p vcf $target_sample/$target_sample.homo.vcf.gz


	#filtering out combined mutations 
	perl $SCRIPT_DIR/FilteringVcfNeighborSNVs.pl < $target_sample/$target_sample.hetero.vcf > $target_sample/$target_sample.non_neighbor.vcf
	bgzip -c $target_sample/$target_sample.non_neighbor.vcf > $target_sample/$target_sample.non_neighbor.vcf.gz
	tabix -f -p vcf $target_sample/$target_sample.non_neighbor.vcf.gz

	perl $SCRIPT_DIR/FilteringVcfNeighborSNVs.pl < $target_sample/$target_sample.homo.vcf > $target_sample/$target_sample.non_neighbor.homo.vcf
	cp $target_sample/$target_sample.non_neighbor.homo.vcf $target_sample/$target_sample.final.mutants.homo.vcf

	#picking up combined mutations
	vcf-isec -c $target_sample/$target_sample.hetero.vcf.gz $target_sample/$target_sample.non_neighbor.vcf.gz > $target_sample/$target_sample.combined.vcf

	#filtering out mutation sites whrere proportion of mutant reads < 25% or GQ < 99
	perl $SCRIPT_DIR/VariantFilteredAF.pl < $target_sample/$target_sample.non_neighbor.vcf > $target_sample/$target_sample.final.mutants.vcf
	perl $SCRIPT_DIR/VariantFilteredAF.pl < $target_sample/$target_sample.combined.vcf > $target_sample/$target_sample.final.combined.vcf

	#output heterozygous mutation.list 
	perl $SCRIPT_DIR/MakeMulationList.pl $target_sample/$target_sample.final.mutants.vcf $target_sample >> $mutation_list_file.unsort.txt
	#output homozygous mutation.list 
	perl $SCRIPT_DIR/MakeMulationList.pl $target_sample/$target_sample.final.mutants.homo.vcf $target_sample >> $mutation_list_file.homo.unsort.txt


	vcf2bed --snvs <  $target_sample/$target_sample.final.mutants.vcf > $target_sample/$target_sample.snp.bed
	vcf2bed --insertions < $target_sample/$target_sample.final.mutants.vcf > $target_sample/$target_sample.insertion.bed
	vcf2bed --deletions <  $target_sample/$target_sample.final.mutants.vcf > $target_sample/$target_sample.deletion.bed
	vcf2bed <  $target_sample/$target_sample.final.mutants.vcf > $target_sample/$target_sample.bed
	
	vcf2bed <  $target_sample/$target_sample.final.combined.vcf > $target_sample/$target_sample.combined.bed
  
	vcf2bed --snvs <  $target_sample/$target_sample.final.mutants.homo.vcf > $target_sample/$target_sample.snp.homo.bed
	vcf2bed --insertions < $target_sample/$target_sample.final.mutants.homo.vcf > $target_sample/$target_sample.insertion.homo.bed
	vcf2bed --deletions <  $target_sample/$target_sample.final.mutants.homo.vcf > $target_sample/$target_sample.deletion.homo.bed
	vcf2bed  <  $target_sample/$target_sample.final.mutants.homo.vcf > $target_sample/$target_sample.homo.bed
	

	no_snp=`wc -l $target_sample/$target_sample.snp.bed |awk '{print $1}'`
	no_insertion=`wc -l $target_sample/$target_sample.insertion.bed |awk '{print $1}'`
	no_deletion=`wc -l $target_sample/$target_sample.deletion.bed |awk '{print $1}'`
	no_indel=$((no_insertion + no_deletion))

	no_snp_homo=`wc -l $target_sample/$target_sample.snp.homo.bed |awk '{print $1}'`
	no_insertion_homo=`wc -l $target_sample/$target_sample.insertion.homo.bed |awk '{print $1}'`
	no_deletion_homo=`wc -l $target_sample/$target_sample.deletion.homo.bed |awk '{print $1}'`
	no_indel_homo=$((no_insertion_homo + no_deletion_homo))
	 
	tab_lab="	"
	output_info=$target_sample$tab_lab$no_snp$tab_lab$no_indel$tab_lab$no_insertion$tab_lab$no_deletion
	echo $output_info >> $mutation_summary_file.txt

	output_info_homo=$target_sample$tab_lab$no_snp_homo$tab_lab$no_indel_homo$tab_lab$no_insertion_homo$tab_lab$no_deletion_homo
	echo $output_info_homo >> $mutation_summary_file_homo.txt

done


#Sort SNPs + INDELs heterozygous mutation list 
cat $mutation_list_file.unsort.txt | sort -k 1,1 -k 2n,2 >  $mutation_list_file.txt
rm $mutation_list_file.unsort.txt
#Extract SNPs heterozygous mutation list
awk '/SNP/' $mutation_list_file.txt > M2.snp.list.txt
#transform to bed for final vcf in pipe.12
perl $SCRIPT_DIR/BioalcidaejdkOut2BED.pl < M2.snp.list.txt > M2.snp.unique.bed
#Extract INDELs heterozygous mutation list
awk '/Insertion/ || /Deletion/' $mutation_list_file.txt > M2.indel.list.txt

#Sort SNPs + INDELs homozygous mutation list 
cat $mutation_list_file.homo.unsort.txt | sort -k 1,1 -k 2n,2 >  $mutation_list_file.homo.txt
rm $mutation_list_file.homo.unsort.txt
#Extract SNPs homozygous mutation list
awk '/SNP/' $mutation_list_file.homo.txt > M2.snp.homo.list.txt
#Extract INDELs homozygous mutation list
awk '/Insertion/ || /Deletion/' $mutation_list_file.homo.txt > M2.indel.homo.list.txt


#Output for M2.indel.unique.vcf
perl $SCRIPT_DIR/BioalcidaejdkOut2BED.pl < M2.indel.list.txt > M2.indel.unique.bed
#select back variants with unique.bed 
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.mu.snp.indel.DPfilterNoCall.vcf.gz\
 -L M2.indel.unique.bed\
 --exclude-sample-name $SCRIPT_DIR/Mother_ID.list\
 -O M2.indel.unique.vcf.gz

 bcftools index M2.indel.unique.vcf.gz



cd $SCRIPT_DIR

module unload gatk/4.1.7.0
module unload vcftools/0.1.15
module unload bedops/2.4.39
