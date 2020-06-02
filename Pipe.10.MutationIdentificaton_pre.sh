#!/bin/bash -i
#Pipe.10.MutationIdentificaton.sh
#by HIRAO Akira

#requirement: FilteringUniqueSNVs.pl, summarizing_Family.unique_snp.R

set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)

module load R/3.6.0
module load gatk/4.1.7.0
module load vcftools/0.1.15


target_ID=AT48

reference_folder=/zfs/Arabidopsis/Reference_v1.1
main_folder=/zfs/Arabidopsis/work/At_Reseq
work_folder=$main_folder/vcf_out

BioAlcidaeJdk_path=/usr/local/jvarkit/dist
#-----------------------------------------------------
# defining the argument for 48 samples
samples_48=()

while read sample; do

	echo $sample
	samples_48+=($sample)

done < $SCRIPT_DIR/sample_list.txt #list of samples

echo $samples_48

#-----------------------------------------------------
mother_vec=(${samples_48[@]:0:12})
M2_vec=(${samples_48[@]:12:36})

echo ${mother_vec[@]}
echo ${M2_vec[@]}

cd $work_folder
#mkdir -p vcf_compare


#----------------------------------------------------------------------------------
#Output mendelian violation site only: --pedigree & --mendelian-violation
GQ_threshold=20
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.snv.DPfilterNoCall.vcf.gz\
 --pedigree $SCRIPT_DIR/AT48.ped\
 --mendelian-violation\
 --mendelian-violation-qual-threshold $GQ_threshold\
 -O $work_folder/$target_ID.mendelian.snv.DPfilterNoCall.vcf.gz
bcftools index -f $work_folder/$target_ID.mendelian.snv.DPfilterNoCall.vcf.gz
bcftools view $work_folder/$target_ID.mendelian.snv.DPfilterNoCall.vcf.gz\
 -Ov -o $work_folder/$target_ID.mendelian.snv.DPfilterNoCall.vcf


#identifying unique SNVs with using bioalcidaejdk.jar
#See detail for https://www.biostars.org/p/329423/#329742
java -jar $BioAlcidaeJdk_path/bioalcidaejdk.jar -e 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0);println(V.getContig()+" "+V.getStart()+" "+V.getReference()+" "+g.getSampleName()+" "+g.getAlleles());});' $work_folder/$target_ID.mendelian.snv.DPfilterNoCall.vcf > $work_folder/$target_ID.unique.snv.list.txt

perl BioalcidaejdkOut2BED.pl < $work_folder/$target_ID.unique.snv.list.txt > $work_folder/$target_ID.unique.bed


gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $work_folder/$target_ID.mendelian.snv.DPfilterNoCall.vcf\
 -L $work_folder/$target_ID.unique.bed\
 --exclude-sample-name $SCRIPT_DIR/Mother_ID.list
 -O $work_folder/AT.M2.unique.vcf




<< COMMENTOUT2

#addition@2020,3.30
#extraction of unique variants: SNPs
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.snp.DPfilterNoCall.vcf.gz\
 -select 'set =="Intersection";'\
 -invertSelect\
 -O $work_folder/$target_ID.snp.unique.vcf

#extraction of unique variants: INDELs
gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.indel.DPfilterNoCall.vcf.gz\
 -select 'set =="Intersection";'\
 -invertSelect\
 -O $work_folder/$target_ID.indel.unique.vcf



perl $SCRIPT_DIR/FilteringUniqueSNVs.pl < $work_folder/$target_ID.unique_snps_list.txt > $work_folder/$target_ID.unique_snps_list.filtered.txt
perl $SCRIPT_DIR/FilteringUniqueSNVs.pl < $work_folder/$target_ID.unique_indels_list.txt > $work_folder/$target_ID.unique_indels_list.filtered.txt
perl $SCRIPT_DIR/FilteringUniqueSNVs.pl < $work_folder/$target_ID.unique_snvs_list.txt > $work_folder/$target_ID.unique_snvs_list.filtered.txt

perl $SCRIPT_DIR/FilteredText2BED.pl < $work_folder/$target_ID.unique_snvs_list.filtered.txt > $work_folder/$target_ID.unique_snvs_list.filtered.bed



cp $SCRIPT_DIR/target_ID.list $work_folder/target_ID.list
#extraction of unique variants: SNPs
#--max-nocall-fraction 0.2 全サンプルの2割のGTが“./."の場合はフィルターアウト
gatk SelectVariants \
-R $reference_folder/TAIR10.fa \
-V $work_folder/$target_ID.snv.DPfilterNoCall.vcf \
-exclude-sample-expressions target_ID.list \
--max-nocall-fraction 0.2 \
-L $work_folder/$target_ID.unique_snvs_list.filtered.bed \
-select-type SNP \
-O $work_folder/$target_ID.M2.unique_snps.vcf

#extraction of unique variants: INDELs
gatk SelectVariants \
-R $reference_folder/TAIR10.fa \
-V $work_folder/$target_ID.snv.DPfilterNoCall.vcf \
-exclude-sample-expressions target_ID.list \
--max-nocall-fraction 0.2 \
-L $work_folder/$target_ID.unique_snvs_list.filtered.bed \
-select-type INDEL \
-O $work_folder/$target_ID.M2.unique_indels.vcf


mutation_summary_file="mutation_summary.txt"
echo -n >| $work_folder/$mutation_summary_file
echo "SampleID"	"SNP_mutations"	"INDEL_mutations" >> $work_folder/$mutation_summary_file

mkdir -p unique_vcfs
mkdir -p mutations_vcf
for target_sample in ${M2_vec[@]}
do

	#Filtering SNPs
	gatk SelectVariants \
	-R $reference_folder/TAIR10.fa \
	-V $work_folder/$target_ID.M2.unique_snps.vcf \
	--sample-name $target_sample \
	-O $work_folder/unique_vcfs/$target_sample.M2.unique_snps.vcf

	#Marking to Filter out homozygous SNP mutations
	gatk VariantFiltration \
	-R $reference_folder/TAIR10.fa \
	-V $work_folder/unique_vcfs/$target_sample.M2.unique_snps.vcf \
	-G-filter "isHomVar==1" \
	-G-filter-name "homozygous_mutation" \
	-G-filter "isHomRef==1" \
	-G-filter-name "homozygous_ref" \
	-O $work_folder/unique_vcfs/$target_sample.M2.unique_snps.hetero_marked.vcf

	#Filtering out homozygous SNP mutations
	gatk SelectVariants \
	-R $reference_folder/TAIR10.fa \
	-V $work_folder/unique_vcfs/$target_sample.M2.unique_snps.hetero_marked.vcf \
	--set-filtered-gt-to-nocall \
	-O $work_folder/unique_vcfs/$target_sample.M2.unique_snps.hetero_nocall.vcf

	perl $SCRIPT_DIR/FilteredNoCallVcf.pl < $work_folder/unique_vcfs/$target_sample.M2.unique_snps.hetero_nocall.vcf > $work_folder/mutations_vcf/$target_sample.M2.snps.mutation.vcf
	bgzip -c $work_folder/mutations_vcf/$target_sample.M2.snps.mutation.vcf > $work_folder/mutations_vcf/$target_sample.M2.snps.mutation.vcf.gz
	tabix -f -p vcf $work_folder/mutations_vcf/$target_sample.M2.snps.mutation.vcf.gz
	no_snp_mutation=$(grep -v "#" $work_folder/mutations_vcf/$target_sample.M2.snps.mutation.vcf | wc -l)


	#Filtering INDELs
	gatk SelectVariants \
	-R $reference_folder/TAIR10.fa \
	-V $work_folder/$target_ID.M2.unique_indels.vcf \
	--sample-name $target_sample \
	-O $work_folder/unique_vcfs/$target_sample.M2.unique_indels.vcf

	#Marking to Filter out homozygous INDEL mutations
	gatk VariantFiltration \
	-R $reference_folder/TAIR10.fa \
	-V $work_folder/unique_vcfs/$target_sample.M2.unique_indels.vcf \
	-G-filter "isHomVar==1" \
	-G-filter-name "homozygous_mutation" \
	-G-filter "isHomRef==1" \
	-G-filter-name "homozygous_ref" \
	-G-filter "GQ < 99" \
	-G-filter-name "GQ99" \
	-O $work_folder/unique_vcfs/$target_sample.M2.unique_indels.hetero_marked.vcf

	#Filtering out homozygous INDEL mutations
	gatk SelectVariants \
	-R $reference_folder/TAIR10.fa \
	-V $work_folder/unique_vcfs/$target_sample.M2.unique_indels.hetero_marked.vcf \
	--set-filtered-gt-to-nocall \
	-O $work_folder/unique_vcfs/$target_sample.M2.unique_indels.hetero_nocall.vcf

	perl $SCRIPT_DIR/FilteredNoCallVcf.pl < $work_folder/unique_vcfs/$target_sample.M2.unique_indels.hetero_nocall.vcf > $work_folder/mutations_vcf/$target_sample.M2.indels.mutation.vcf
	bgzip -c $work_folder/mutations_vcf/$target_sample.M2.indels.mutation.vcf > $work_folder/mutations_vcf/$target_sample.M2.indels.mutation.vcf.gz
	tabix -f -p vcf $work_folder/mutations_vcf/$target_sample.M2.indels.mutation.vcf.gz
	no_indel_mutation=$(grep -v "#" $work_folder/mutations_vcf/$target_sample.M2.indels.mutation.vcf | wc -l)

	echo "${target_sample}"	"${no_snp_mutation}"	"${no_indel_mutation}" >> $work_folder/$mutation_summary_file

done




vcf-merge \
$work_folder/mutations_vcf/$sample_No_13.M2.snps.mutation.vcf.gz \
...
$work_folder/mutations_vcf/$sample_No_48.M2.snps.mutation.vcf.gz | bgzip -c > $work_folder/$target_ID.M2.snps.mutation.vcf.gz


vcf-merge \
$work_folder/mutations_vcf/$sample_No_13.M2.indels.mutation.vcf.gz \
...

$work_folder/mutations_vcf/$sample_No_48.M2.indels.mutation.vcf.gz | bgzip -c > $work_folder/$target_ID.M2.indels.mutation.vcf.gz


gunzip -c $work_folder/$target_ID.M2.snps.mutation.vcf.gz> $work_folder/$target_ID.M2.snps.mutation.vcf
gunzip -c $work_folder/$target_ID.M2.indels.mutation.vcf.gz> $work_folder/$target_ID.M2.indels.mutation.vcf


#identifying unique SNVs with using bioalcidaejdk.jar
java -jar /usr/local/jvarkit/dist/bioalcidaejdk.jar -e 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0);println(V.getContig()+" "+V.getStart()+" "+V.getReference()+" "+g.getSampleName()+" "+g.getAlleles());});' $work_folder/$target_ID.M2.snps.mutation.vcf > $work_folder/$target_ID.M2.snps.mutation.list.txt
java -jar /usr/local/jvarkit/dist/bioalcidaejdk.jar -e 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0);println(V.getContig()+" "+V.getStart()+" "+V.getReference()+" "+g.getSampleName()+" "+g.getAlleles());});' $work_folder/$target_ID.M2.indels.mutation.vcf > $work_folder/$target_ID.M2.indels.mutation.list.txt


#cd $SCRIPT_DIR/$work_folder

cp $SCRIPT_DIR/summarizing_Family.unique_snp.R $work_folder/summarizing_Family.unique_snp.R
cp $SCRIPT_DIR/summarizing_Family.unique_indel.R $work_folder/summarizing_Family.unique_indel.R
cp $SCRIPT_DIR/summarizing_Family.unique_snv.R $work_folder/summarizing_Family.unique_snv.R

R --no-save $target_ID $sample_No_01 $sample_No_02 $sample_No_03 $sample_No_04 $sample_No_05 $sample_No_06 $sample_No_07 $sample_No_08 $sample_No_09 $sample_No_10 $sample_No_11 $sample_No_12 $sample_No_13 $sample_No_14 $sample_No_15 $sample_No_16 $sample_No_17 $sample_No_18 $sample_No_19 $sample_No_20 $sample_No_21 $sample_No_22 $sample_No_23 $sample_No_24 $sample_No_25 $sample_No_26 $sample_No_27 $sample_No_28 $sample_No_29 $sample_No_30 $sample_No_31 $sample_No_32 $sample_No_33 $sample_No_34 $sample_No_35 $sample_No_36 $sample_No_37 $sample_No_38 $sample_No_39 $sample_No_40 $sample_No_41 $sample_No_42 $sample_No_43 $sample_No_44 $sample_No_45 $sample_No_46 $sample_No_47 $sample_No_48 < summarizing_Family.unique_snp.R 
R --no-save $target_ID $sample_No_01 $sample_No_02 $sample_No_03 $sample_No_04 $sample_No_05 $sample_No_06 $sample_No_07 $sample_No_08 $sample_No_09 $sample_No_10 $sample_No_11 $sample_No_12 $sample_No_13 $sample_No_14 $sample_No_15 $sample_No_16 $sample_No_17 $sample_No_18 $sample_No_19 $sample_No_20 $sample_No_21 $sample_No_22 $sample_No_23 $sample_No_24 $sample_No_25 $sample_No_26 $sample_No_27 $sample_No_28 $sample_No_29 $sample_No_30 $sample_No_31 $sample_No_32 $sample_No_33 $sample_No_34 $sample_No_35 $sample_No_36 $sample_No_37 $sample_No_38 $sample_No_39 $sample_No_40 $sample_No_41 $sample_No_42 $sample_No_43 $sample_No_44 $sample_No_45 $sample_No_46 $sample_No_47 $sample_No_48 < summarizing_Family.unique_indel.R 
R --no-save $target_ID $sample_No_01 $sample_No_02 $sample_No_03 $sample_No_04 $sample_No_05 $sample_No_06 $sample_No_07 $sample_No_08 $sample_No_09 $sample_No_10 $sample_No_11 $sample_No_12 $sample_No_13 $sample_No_14 $sample_No_15 $sample_No_16 $sample_No_17 $sample_No_18 $sample_No_19 $sample_No_20 $sample_No_21 $sample_No_22 $sample_No_23 $sample_No_24 $sample_No_25 $sample_No_26 $sample_No_27 $sample_No_28 $sample_No_29 $sample_No_30 $sample_No_31 $sample_No_32 $sample_No_33 $sample_No_34 $sample_No_35 $sample_No_36 $sample_No_37 $sample_No_38 $sample_No_39 $sample_No_40 $sample_No_41 $sample_No_42 $sample_No_43 $sample_No_44 $sample_No_45 $sample_No_46 $sample_No_47 $sample_No_48 < summarizing_Family.unique_snv.R 


COMMENTOUT2

cd $SCRIPT_DIR


module unload R/3.6.0
module unload gatk/4.1.7.0
module unload vcftools/0.1.15