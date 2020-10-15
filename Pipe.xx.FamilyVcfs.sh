#!/bin/bash -i
#Pipe.xx.FamilyVcfs.sh
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
vcf_family_folder=$main_folder/vcf_out_family

BioAlcidaeJdk_path=/usr/local/jvarkit/dist


cd $vcf_family_folder


#-----------------------------------------------------
# defining the argument for 48 samples

while read family sample1 sample2 sample3; do
	echo $sample1
	echo $sample2
	echo $sample3

	gatk SelectVariants\
	 -R $reference_folder/TAIR10.fa\
	 -V $main_folder/vcf_out/$target_ID.mu.snp.indel.DPfilterNoCall.vcf.gz\
	 --sample-name $sample1 --sample-name $sample2 --sample-name $sample3\
	 -O $family.mu.snp.indel.DPfilterNoCall.vcf.gz

	gunzip -c $family.mu.snp.indel.DPfilterNoCall.vcf.gz > $family.mu.snp.indel.DPfilterNoCall.vcf

	java -jar $BioAlcidaeJdk_path/bioalcidaejdk.jar -e 'stream().forEach(V->{final List<Genotype> L=V.getGenotypes().stream().filter(G->G.isHet() || G.isHomVar()).collect(Collectors.toList());if(L.size()!=1) return;final Genotype g=L.get(0);println(V.getContig()+" "+V.getStart()+" "+V.getReference()+" "+g.getSampleName()+" "+g.getAlleles());});' $family.mu.snp.indel.DPfilterNoCall.vcf > $family.mu.snp.indel.uniuqe.DPfilterNoCall.list

	perl $SCRIPT_DIR/BioalcidaejdkOut2BED.pl < $family.mu.snp.indel.uniuqe.DPfilterNoCall.list > $family.mu.snp.indel.unique.bed


	gatk SelectVariants\
	 -R $reference_folder/TAIR10.fa\
	 -V $family.mu.snp.indel.DPfilterNoCall.vcf.gz\
	 -XL $family.mu.snp.indel.unique.bed\
	 -O $family.mu.snp.indel.nonunique.DPfilterNoCall.vcf.gz

	gunzip -c $family.mu.snp.indel.nonunique.DPfilterNoCall.vcf.gz > $family.mu.snp.indel.nonunique.DPfilterNoCall.vcf

	vcf2bed < $family.mu.snp.indel.nonunique.DPfilterNoCall.vcf > $family.mu.snp.indel.nonunique.DPfilterNoCall.bed
	
	perl $SCRIPT_DIR/Vcf2Position.pl <  $family.mu.snp.indel.nonunique.DPfilterNoCall.vcf >  $family.mu.snp.indel.nonunique.DPfilterNoCall.position.txt

	perl $SCRIPT_DIR/BioalcidaejdkOut2BED.pl < $family.mu.snp.indel.uniuqe.DPfilterNoCall.list > $family.mu.snp.indel.unique.bed


done < $SCRIPT_DIR/family.list #list of samples


#bcftools isec A010.mu.snp.indel.nonunique.DPfilterNoCall.vcf.gz\
# A010.mu.snp.indel.nonunique.DPfilterNoCall.vcf.gz\
# A020.mu.snp.indel.nonunique.DPfilterNoCall.vcf.gz\
# A030.mu.snp.indel.nonunique.DPfilterNoCall.vcf.gz\
# A110.mu.snp.indel.nonunique.DPfilterNoCall.vcf.gz\
# A120.mu.snp.indel.nonunique.DPfilterNoCall.vcf.gz\
# A130.mu.snp.indel.nonunique.DPfilterNoCall.vcf.gz\
# A210.mu.snp.indel.nonunique.DPfilterNoCall.vcf.gz\
# A220.mu.snp.indel.nonunique.DPfilterNoCall.vcf.gz\
# A230.mu.snp.indel.nonunique.DPfilterNoCall.vcf.gz\
# A310.mu.snp.indel.nonunique.DPfilterNoCall.vcf.gz\
# A320.mu.snp.indel.nonunique.DPfilterNoCall.vcf.gz\
# A330.mu.snp.indel.nonunique.DPfilterNoCall.vcf.gz\
# -p famiy.unique -n-1




cd $SCRIPT_DIR





module unload gatk/4.1.7.0
module unload vcftools/0.1.15
module unload bedops/2.4.39
