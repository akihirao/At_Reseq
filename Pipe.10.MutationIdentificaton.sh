#!/bin/bash -i
#Pipe.10.MutationIdentificaton.sh
#by HIRAO Akira

#requirement: FilteringUniqueSNVs.pl, summarizing_Family.unique_snp.R

module load R
module load gatk
module load vcftools

set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)


mother_ID=all
list_file_tag=fastq_list_
list_file=$list_file_tag$mother_ID
sample_list=($(grep '^A' $SCRIPT_DIR/$list_file.txt)) 
#sample_list=($(grep '^A' $SCRIPT_DIR/fastq_list_AT31.txt)) 
#sample_list=($(grep '^AT[0-9]\{2,3\}' $SCRIPT_DIR/fastq_list_AT01.txt)) 
sample_No_01=${sample_list[0]}
sample_No_02=${sample_list[3]}
sample_No_03=${sample_list[6]}
sample_No_04=${sample_list[9]}
sample_No_05=${sample_list[12]}
sample_No_06=${sample_list[15]}
sample_No_07=${sample_list[18]}
sample_No_08=${sample_list[21]}
sample_No_09=${sample_list[24]}
sample_No_10=${sample_list[27]}
sample_No_11=${sample_list[30]}
sample_No_12=${sample_list[33]}
sample_No_13=${sample_list[36]}
sample_No_14=${sample_list[39]}
sample_No_15=${sample_list[42]}
sample_No_16=${sample_list[45]}
sample_No_17=${sample_list[48]}
sample_No_18=${sample_list[51]}
sample_No_19=${sample_list[54]}
sample_No_20=${sample_list[57]}
sample_No_21=${sample_list[60]}
sample_No_22=${sample_list[63]}
sample_No_23=${sample_list[66]}
sample_No_24=${sample_list[69]}
sample_No_25=${sample_list[72]}
sample_No_26=${sample_list[75]}
sample_No_27=${sample_list[78]}
sample_No_28=${sample_list[81]}
sample_No_29=${sample_list[84]}
sample_No_30=${sample_list[87]}
sample_No_31=${sample_list[90]}
sample_No_32=${sample_list[93]}
sample_No_33=${sample_list[96]}
sample_No_34=${sample_list[99]}
sample_No_35=${sample_list[102]}
sample_No_36=${sample_list[105]}
sample_No_37=${sample_list[108]}
sample_No_38=${sample_list[111]}
sample_No_39=${sample_list[114]}
sample_No_40=${sample_list[117]}
sample_No_41=${sample_list[120]}
sample_No_42=${sample_list[123]}
sample_No_43=${sample_list[126]}
sample_No_44=${sample_list[129]}
sample_No_45=${sample_list[132]}
sample_No_46=${sample_list[135]}
sample_No_47=${sample_list[138]}
sample_No_48=${sample_list[141]}

M2_ID_vec=(${sample_No_13} ${sample_No_14} ${sample_No_15} ${sample_No_16}  ${sample_No_17}  ${sample_No_18}  ${sample_No_19}  ${sample_No_20}  ${sample_No_21}  ${sample_No_22}  ${sample_No_23}  ${sample_No_24}  ${sample_No_25}  ${sample_No_26}  ${sample_No_27}  ${sample_No_28}  ${sample_No_29}  ${sample_No_30}  ${sample_No_31}  ${sample_No_32}  ${sample_No_33}  ${sample_No_34}  ${sample_No_35}  ${sample_No_36}  ${sample_No_37}  ${sample_No_38}  ${sample_No_39} ${sample_No_40} ${sample_No_41} ${sample_No_42} ${sample_No_43} ${sample_No_44} ${sample_No_45} ${sample_No_46} ${sample_No_47} ${sample_No_48})
#echo "${M2_ID_vec[0]}","${M2_ID_vec[1]}","${M2_ID_vec[2]}","${M2_ID_vec[3]}" 

reference_folder=/zfs/Arabidopsis/Reference_v1.1
input_folder=/zfs/Arabidopsis/work/FamilySNV/$mother_ID.family
output_folder=/zfs/Arabidopsis/work/FamilySNV/$mother_ID.family

cd $output_folder
mkdir -p vcf_compare

cp $SCRIPT_DIR/summarizing_Family.unique_snp.R $output_folder/summarizing_Family.unique_snp.R
cp $SCRIPT_DIR/summarizing_Family.unique_indel.R $output_folder/summarizing_Family.unique_indel.R
cp $SCRIPT_DIR/summarizing_Family.unique_snv.R $output_folder/summarizing_Family.unique_snv.R


#addition@2020,3.30
#extraction of unique variants: SNPs
#gatk SelectVariants \
#-R $reference_folder/TAIR10.fa \
#-V $mother_ID.family.snp.DPfilterNoCall.vcf.gz \
#-select 'set =="Intersection";' -invertSelect \
#-O $output_folder/$mother_ID.family.snp.unique.vcf

#extraction of unique variants: INDELs
#gatk SelectVariants \
#-R $reference_folder/TAIR10.fa \
#-V $mother_ID.family.indel.DPfilterNoCall.vcf.gz \
#-select 'set =="Intersection";' -invertSelect \
#-O $output_folder/$mother_ID.family.indel.unique.vcf


#snpsindels
gunzip -c $input_folder/$mother_ID.family.snp.DPfilterNoCall.vcf.gz > $input_folder/$mother_ID.family.snp.DPfilterNoCall.vcf
gunzip -c $input_folder/$mother_ID.family.indel.DPfilterNoCall.vcf.gz > $input_folder/$mother_ID.family.indel.DPfilterNoCall.vcf
gunzip -c $input_folder/$mother_ID.family.snv.DPfilterNoCall.vcf.gz > $input_folder/$mother_ID.family.snv.DPfilterNoCall.vcf


#identifying unique SNVs with using bioalcidaejdk.jar
#See detail for https://www.biostars.org/p/329423/#329742
java -jar /usr/local/jvarkit/dist/bioalcidaejdk.jar -e 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0);println(V.getContig()+" "+V.getStart()+" "+V.getReference()+" "+g.getSampleName()+" "+g.getAlleles());});' $input_folder/$mother_ID.family.snp.DPfilterNoCall.vcf > $output_folder/$mother_ID.family.unique_snps_list.txt
java -jar /usr/local/jvarkit/dist/bioalcidaejdk.jar -e 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0);println(V.getContig()+" "+V.getStart()+" "+V.getReference()+" "+g.getSampleName()+" "+g.getAlleles());});' $input_folder/$mother_ID.family.indel.DPfilterNoCall.vcf > $output_folder/$mother_ID.family.unique_indels_list.txt
java -jar /usr/local/jvarkit/dist/bioalcidaejdk.jar -e 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0);println(V.getContig()+" "+V.getStart()+" "+V.getReference()+" "+g.getSampleName()+" "+g.getAlleles());});' $input_folder/$mother_ID.family.snv.DPfilterNoCall.vcf > $output_folder/$mother_ID.family.unique_snvs_list.txt


perl $SCRIPT_DIR/FilteringUniqueSNVs.pl < $output_folder/$mother_ID.family.unique_snps_list.txt > $output_folder/$mother_ID.family.unique_snps_list.filtered.txt
perl $SCRIPT_DIR/FilteringUniqueSNVs.pl < $output_folder/$mother_ID.family.unique_indels_list.txt > $output_folder/$mother_ID.family.unique_indels_list.filtered.txt
perl $SCRIPT_DIR/FilteringUniqueSNVs.pl < $output_folder/$mother_ID.family.unique_snvs_list.txt > $output_folder/$mother_ID.family.unique_snvs_list.filtered.txt

perl $SCRIPT_DIR/FilteredText2BED.pl < $output_folder/$mother_ID.family.unique_snvs_list.filtered.txt > $output_folder/$mother_ID.family.unique_snvs_list.filtered.bed



cp $SCRIPT_DIR/Mother_ID.list $output_folder/Mother_ID.list
#extraction of unique variants: SNPs
#--max-nocall-fraction 0.2 全サンプルの2割のGTが“./."の場合はフィルターアウト
gatk SelectVariants \
-R $reference_folder/TAIR10.fa \
-V $output_folder/$mother_ID.family.snv.DPfilterNoCall.vcf \
-exclude-sample-expressions Mother_ID.list \
--max-nocall-fraction 0.2 \
-L $output_folder/$mother_ID.family.unique_snvs_list.filtered.bed \
-select-type SNP \
-O $output_folder/$mother_ID.family.M2.unique_snps.vcf

#extraction of unique variants: INDELs
gatk SelectVariants \
-R $reference_folder/TAIR10.fa \
-V $output_folder/$mother_ID.family.snv.DPfilterNoCall.vcf \
-exclude-sample-expressions Mother_ID.list \
--max-nocall-fraction 0.2 \
-L $output_folder/$mother_ID.family.unique_snvs_list.filtered.bed \
-select-type INDEL \
-O $output_folder/$mother_ID.family.M2.unique_indels.vcf


mutation_summary_file="mutation_summary.txt"
echo -n >| $output_folder/$mutation_summary_file
echo "SampleID"	"SNP_mutations"	"INDEL_mutations" >> $output_folder/$mutation_summary_file

mkdir -p unique_vcfs
mkdir -p mutations_vcf
for target_sample in A011	A012	A014	A021	A023	A024	A031	A033	A034	A113	A114	A115	A121	A123	A125	A131	A132	A135	A212	A213	A215	A221	A222	A225	A231	A233	A234	A311	A313	A314	A321	A323	A324	A331	A332	A334
do

	#Filtering SNPs
	gatk SelectVariants \
	-R $reference_folder/TAIR10.fa \
	-V $output_folder/$mother_ID.family.M2.unique_snps.vcf \
	--sample-name $target_sample \
	-O $output_folder/unique_vcfs/$target_sample.M2.unique_snps.vcf

	#Marking to Filter out homozygous SNP mutations
	gatk VariantFiltration \
	-R $reference_folder/TAIR10.fa \
	-V $output_folder/unique_vcfs/$target_sample.M2.unique_snps.vcf \
	-G-filter "isHomVar==1" \
	-G-filter-name "homozygous_mutation" \
	-G-filter "isHomRef==1" \
	-G-filter-name "homozygous_ref" \
	-O $output_folder/unique_vcfs/$target_sample.M2.unique_snps.hetero_marked.vcf

	#Filtering out homozygous SNP mutations
	gatk SelectVariants \
	-R $reference_folder/TAIR10.fa \
	-V $output_folder/unique_vcfs/$target_sample.M2.unique_snps.hetero_marked.vcf \
	--set-filtered-gt-to-nocall \
	-O $output_folder/unique_vcfs/$target_sample.M2.unique_snps.hetero_nocall.vcf

	perl $SCRIPT_DIR/FilteredNoCallVcf.pl < $output_folder/unique_vcfs/$target_sample.M2.unique_snps.hetero_nocall.vcf > $output_folder/mutations_vcf/$target_sample.M2.snps.mutation.vcf
	bgzip -c $output_folder/mutations_vcf/$target_sample.M2.snps.mutation.vcf > $output_folder/mutations_vcf/$target_sample.M2.snps.mutation.vcf.gz
	tabix -f -p vcf $output_folder/mutations_vcf/$target_sample.M2.snps.mutation.vcf.gz
	no_snp_mutation=$(grep -v "#" $output_folder/mutations_vcf/$target_sample.M2.snps.mutation.vcf | wc -l)


	#Filtering INDELs
	gatk SelectVariants \
	-R $reference_folder/TAIR10.fa \
	-V $output_folder/$mother_ID.family.M2.unique_indels.vcf \
	--sample-name $target_sample \
	-O $output_folder/unique_vcfs/$target_sample.M2.unique_indels.vcf

	#Marking to Filter out homozygous INDEL mutations
	gatk VariantFiltration \
	-R $reference_folder/TAIR10.fa \
	-V $output_folder/unique_vcfs/$target_sample.M2.unique_indels.vcf \
	-G-filter "isHomVar==1" \
	-G-filter-name "homozygous_mutation" \
	-G-filter "isHomRef==1" \
	-G-filter-name "homozygous_ref" \
	-G-filter "GQ < 99" \
	-G-filter-name "GQ99" \
	-O $output_folder/unique_vcfs/$target_sample.M2.unique_indels.hetero_marked.vcf

	#Filtering out homozygous INDEL mutations
	gatk SelectVariants \
	-R $reference_folder/TAIR10.fa \
	-V $output_folder/unique_vcfs/$target_sample.M2.unique_indels.hetero_marked.vcf \
	--set-filtered-gt-to-nocall \
	-O $output_folder/unique_vcfs/$target_sample.M2.unique_indels.hetero_nocall.vcf

	perl $SCRIPT_DIR/FilteredNoCallVcf.pl < $output_folder/unique_vcfs/$target_sample.M2.unique_indels.hetero_nocall.vcf > $output_folder/mutations_vcf/$target_sample.M2.indels.mutation.vcf
	bgzip -c $output_folder/mutations_vcf/$target_sample.M2.indels.mutation.vcf > $output_folder/mutations_vcf/$target_sample.M2.indels.mutation.vcf.gz
	tabix -f -p vcf $output_folder/mutations_vcf/$target_sample.M2.indels.mutation.vcf.gz
	no_indel_mutation=$(grep -v "#" $output_folder/mutations_vcf/$target_sample.M2.indels.mutation.vcf | wc -l)

	echo "${target_sample}"	"${no_snp_mutation}"	"${no_indel_mutation}" >> $output_folder/$mutation_summary_file

done




vcf-merge \
$output_folder/mutations_vcf/$sample_No_13.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_14.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_15.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_16.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_17.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_18.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_19.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_20.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_21.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_22.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_23.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_24.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_25.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_26.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_27.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_28.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_29.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_30.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_31.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_32.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_33.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_34.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_35.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_36.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_37.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_38.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_39.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_40.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_41.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_42.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_43.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_44.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_45.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_46.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_47.M2.snps.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_48.M2.snps.mutation.vcf.gz | bgzip -c > $output_folder/$mother_ID.family.M2.snps.mutation.vcf.gz


vcf-merge \
$output_folder/mutations_vcf/$sample_No_13.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_14.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_15.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_16.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_17.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_18.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_19.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_20.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_21.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_22.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_23.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_24.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_25.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_26.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_27.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_28.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_29.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_30.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_31.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_32.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_33.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_34.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_35.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_36.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_37.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_38.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_39.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_40.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_41.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_42.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_43.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_44.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_45.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_46.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_47.M2.indels.mutation.vcf.gz \
$output_folder/mutations_vcf/$sample_No_48.M2.indels.mutation.vcf.gz | bgzip -c > $output_folder/$mother_ID.family.M2.indels.mutation.vcf.gz


gunzip -c $output_folder/$mother_ID.family.M2.snps.mutation.vcf.gz> $output_folder/$mother_ID.family.M2.snps.mutation.vcf
gunzip -c $output_folder/$mother_ID.family.M2.indels.mutation.vcf.gz> $output_folder/$mother_ID.family.M2.indels.mutation.vcf


#identifying unique SNVs with using bioalcidaejdk.jar
java -jar /usr/local/jvarkit/dist/bioalcidaejdk.jar -e 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0);println(V.getContig()+" "+V.getStart()+" "+V.getReference()+" "+g.getSampleName()+" "+g.getAlleles());});' $output_folder/$mother_ID.family.M2.snps.mutation.vcf > $output_folder/$mother_ID.family.M2.snps.mutation.list.txt
java -jar /usr/local/jvarkit/dist/bioalcidaejdk.jar -e 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0);println(V.getContig()+" "+V.getStart()+" "+V.getReference()+" "+g.getSampleName()+" "+g.getAlleles());});' $output_folder/$mother_ID.family.M2.indels.mutation.vcf > $output_folder/$mother_ID.family.M2.indels.mutation.list.txt


#cd $SCRIPT_DIR/$output_folder

R --no-save $mother_ID $sample_No_01 $sample_No_02 $sample_No_03 $sample_No_04 $sample_No_05 $sample_No_06 $sample_No_07 $sample_No_08 $sample_No_09 $sample_No_10 $sample_No_11 $sample_No_12 $sample_No_13 $sample_No_14 $sample_No_15 $sample_No_16 $sample_No_17 $sample_No_18 $sample_No_19 $sample_No_20 $sample_No_21 $sample_No_22 $sample_No_23 $sample_No_24 $sample_No_25 $sample_No_26 $sample_No_27 $sample_No_28 $sample_No_29 $sample_No_30 $sample_No_31 $sample_No_32 $sample_No_33 $sample_No_34 $sample_No_35 $sample_No_36 $sample_No_37 $sample_No_38 $sample_No_39 $sample_No_40 $sample_No_41 $sample_No_42 $sample_No_43 $sample_No_44 $sample_No_45 $sample_No_46 $sample_No_47 $sample_No_48 < summarizing_Family.unique_snp.R 
R --no-save $mother_ID $sample_No_01 $sample_No_02 $sample_No_03 $sample_No_04 $sample_No_05 $sample_No_06 $sample_No_07 $sample_No_08 $sample_No_09 $sample_No_10 $sample_No_11 $sample_No_12 $sample_No_13 $sample_No_14 $sample_No_15 $sample_No_16 $sample_No_17 $sample_No_18 $sample_No_19 $sample_No_20 $sample_No_21 $sample_No_22 $sample_No_23 $sample_No_24 $sample_No_25 $sample_No_26 $sample_No_27 $sample_No_28 $sample_No_29 $sample_No_30 $sample_No_31 $sample_No_32 $sample_No_33 $sample_No_34 $sample_No_35 $sample_No_36 $sample_No_37 $sample_No_38 $sample_No_39 $sample_No_40 $sample_No_41 $sample_No_42 $sample_No_43 $sample_No_44 $sample_No_45 $sample_No_46 $sample_No_47 $sample_No_48 < summarizing_Family.unique_indel.R 
R --no-save $mother_ID $sample_No_01 $sample_No_02 $sample_No_03 $sample_No_04 $sample_No_05 $sample_No_06 $sample_No_07 $sample_No_08 $sample_No_09 $sample_No_10 $sample_No_11 $sample_No_12 $sample_No_13 $sample_No_14 $sample_No_15 $sample_No_16 $sample_No_17 $sample_No_18 $sample_No_19 $sample_No_20 $sample_No_21 $sample_No_22 $sample_No_23 $sample_No_24 $sample_No_25 $sample_No_26 $sample_No_27 $sample_No_28 $sample_No_29 $sample_No_30 $sample_No_31 $sample_No_32 $sample_No_33 $sample_No_34 $sample_No_35 $sample_No_36 $sample_No_37 $sample_No_38 $sample_No_39 $sample_No_40 $sample_No_41 $sample_No_42 $sample_No_43 $sample_No_44 $sample_No_45 $sample_No_46 $sample_No_47 $sample_No_48 < summarizing_Family.unique_snv.R 


cd $SCRIPT_DIR

module unload R
module unload gatk
module unload vcftools