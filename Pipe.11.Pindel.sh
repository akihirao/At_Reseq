#!/bin/bash -i
#Pipe.11.Pindel.sh
#by HIRAO Akira


set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)

CPU=12

reference_folder=/zfs/Arabidopsis/Reference_v1.1
main_folder=/zfs/Arabidopsis/work/At_Reseq
bwa_folder=$main_folder/bwa_out
work_folder=$main_folder/pindel_out
BioAlcidaeJdk_path=/usr/local/jvarkit/dist
bwa_path=/zfs/Arabidopsis/work/At_Reseq/bwa_out
slash_lab="/"
tab_lab="	"
insert_length=450

mkdir -p $work_folder
cd $work_folder

module load miniconda2
module load vcftools/0.1.15


while read sample fastq_R1 fastq_R2; do

	mkdir -p $sample
	cd $sample

	config_contents=$bwa_path$slash_lab$sample$slash_lab$sample.TAIR10.bqsr.bam$tab_lab$insert_length$tab_lab$sample
	echo $config_contents > $sample.config.txt

	#Identifying structural variants with using pindel
	#six types of structural variants（D = deletion、SI = short insertion、INV = inversion、TD = tandem duplication、LI = large insertion、BP = unassigned breakpoints）
	#if Option "-c ALL" does not work well in the environment, please work for each of chromosomes   
	/usr/local/miniconda2/bin/pindel -T $CPU -f $reference_folder/TAIR10.fa -i $sample.config.txt -c ALL -o $sample.pindel.out

	#conversion to vcf: six types of structural variants（D = deletion、SI = short insertion、INV = inversion、TD = tandem duplication、LI = large insertion、BP = unassigned breakpoints）
	/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/$sample/$sample.pindel.out_D -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/$sample/$sample.pindel_out.D.vcf -e 10
	/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/$sample/$sample.pindel.out_SI -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/$sample/$sample.pindel_out.SI.vcf -e 10
	/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/$sample/$sample.pindel.out_INV -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/$sample/$sample.pindel_out.INV.vcf -e 10
	/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/$sample/$sample.pindel.out_TD -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/$sample/$sample.pindel_out.TD.vcf -e 10
	/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/$sample/$sample.pindel.out_LI -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/$sample/$sample.pindel_out.LI.vcf -e 10
	/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/$sample/$sample.pindel.out_BP -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/$sample/$sample.pindel_out.BP.vcf -e 10
	/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/$sample/$sample.pindel.out_RP -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/$sample/$sample.pindel_out.RP.vcf -e 10

	#filtering out neighborhood mutations 
	perl $SCRIPT_DIR/FilteringVcfNeighborSNVs.pindel.pl < $sample.pindel_out.D.vcf > $sample.pindel_out.D.non_neighbor.vcf
	perl $SCRIPT_DIR/FilteringVcfNeighborSNVs.pindel.pl < $sample.pindel_out.SI.vcf > $sample.pindel_out.SI.non_neighbor.vcf


	#-----------------------------------------------------
	#Identification of mutations
	#six types of structural variants（D = deletion、SI = short insertion、INV = inversion、TD = tandem duplication、LI = large insertion、BP = unassigned breakpoints）
	bgzip -c $sample.pindel_out.D.non_neighbor.vcf > $sample.pindel_out.D.non_neighbor.vcf.gz
	tabix -f -p vcf $sample.pindel_out.D.non_neighbor.vcf.gz

	bgzip -c $sample.pindel_out.SI.non_neighbor.vcf > $sample.pindel_out.SI.non_neighbor.vcf.gz
	tabix -f -p vcf $sample.pindel_out.SI.non_neighbor.vcf.gz

	vcf-concat $sample.pindel_out.D.non_neighbor.vcf.gz $sample.pindel_out.SI.non_neighbor.vcf.gz > $sample.pindel.indel.unsorted.vcf
	vcf-sort $sample.pindel.indel.unsorted.vcf > $sample.pindel.indel.vcf
	bgzip -c $sample.pindel.indel.vcf > $sample.pindel.indel.vcf.gz
	tabix -f -p vcf $sample.pindel.indel.vcf.gz
	bcftools index $sample.pindel.indel.vcf.gz

	cd $work_folder

done < $SCRIPT_DIR/fastq_list_non1.txt #list of MIDs and folder deposited their read fastq.gz
#one < $SCRIPT_DIR/fastq_list.txt #list of MIDs and folder deposited their read fastq.gz


cd $SCRIPT_DIR

module unload miniconda2
module unload vcftools/0.1.15


