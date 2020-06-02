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






#identifying unique SNVs with using bioalcidaejdk.jar
#See detail for https://www.biostars.org/p/329423/#329742
#java -jar $BioAlcidaeJdk_path/bioalcidaejdk.jar -e 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0);println(V.getContig()+" "+V.getStart()+" "+V.getReference()+" "+g.getSampleName()+" "+g.getAlleles());});' $work_folder/$target_ID.mendelian.snp.DPfilterNoCall.vcf > $work_folder/$target_ID.unique_snps_list.txt
#java -jar $BioAlcidaeJdk_path/bioalcidaejdk.jar -e 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0);println(V.getContig()+" "+V.getStart()+" "+V.getReference()+" "+g.getSampleName()+" "+g.getAlleles());});' $work_folder/$target_ID.mendelian.indel.DPfilterNoCall.vcf > $work_folder/$target_ID.unique_indels_list.txt
java -jar $BioAlcidaeJdk_path/bioalcidaejdk.jar -e 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0);println(V.getContig()+" "+V.getStart()+" "+V.getReference()+" "+g.getSampleName()+" "+g.getAlleles());});' $work_folder/$target_ID.mendelian.snv.DPfilterNoCall.vcf > $work_folder/$target_ID.unique.snv.list.txt

gatk SelectVariants\
 -R $reference_folder/TAIR10.fa\
 -V $target_ID.snv.DPfilterNoCall.vcf.gz\
 -select 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0)'\
 -O $work_folder/$target_ID.snv.unique.text.vcf




cd $SCRIPT_DIR


module unload R/3.6.0
module unload gatk/4.1.7.0
module unload vcftools/0.1.15