#!/bin/bash -i
#Pipe.11.pindel.sh
#by HIRAO Akira


set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)

CPU=16

reference_folder=/zfs/Arabidopsis/Reference_v1.1
main_folder=/zfs/Arabidopsis/work/At_Reseq
bwa_folder=$main_folder/bwa_out
work_folder=$main_folder/pindel_out
BioAlcidaeJdk_path=/usr/local/jvarkit/dist


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


mkdir -p $work_folder

module load miniconda2
module load vcftools/0.1.15

#Identifying structural variants with using pindel
#six types of structural variants（D = deletion、SI = short insertion、INV = inversion、TD = tandem duplication、LI = large insertion、BP = unassigned breakpoints）
#if Option "-c ALL" does not work well in the environment, please work for each of chromosomes   
/usr/local/miniconda2/bin/pindel -T $CPU -f $reference_folder/TAIR10.fa -i $SCRIPT_DIR/AT48.pindel.config.txt -c ALL -o $work_folder/AT48.pindel.out
#/usr/local/miniconda2/bin/pindel -T $CPU -f $reference_folder/TAIR10.fa -i $SCRIPT_DIR/AT48.pindel.config.txt -c Chr1 -o $work_folder/AT48.pindel.chr1.out
#/usr/local/miniconda2/bin/pindel -T $CPU -f $reference_folder/TAIR10.fa -i $SCRIPT_DIR/AT48.pindel.config.txt -c Chr2 -o $work_folder/AT48.pindel.chr2.out
#/usr/local/miniconda2/bin/pindel -T $CPU -f $reference_folder/TAIR10.fa -i $SCRIPT_DIR/AT48.pindel.config.txt -c Chr3 -o $work_folder/AT48.pindel.chr3.out
#/usr/local/miniconda2/bin/pindel -T $CPU -f $reference_folder/TAIR10.fa -i $SCRIPT_DIR/AT48.pindel.config.txt -c Chr4 -o $work_folder/AT48.pindel.chr4.out
#/usr/local/miniconda2/bin/pindel -T $CPU -f $reference_folder/TAIR10.fa -i $SCRIPT_DIR/AT48.pindel.config.txt -c Chr5 -o $work_folder/AT48.pindel.chr5.out


#conversion to vcf: six types of structural variants（D = deletion、SI = short insertion、INV = inversion、TD = tandem duplication、LI = large insertion、BP = unassigned breakpoints）
/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/AT48.pindel.out_D -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/AT48.pindel_out.D.vcf -e 10
/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/AT48.pindel.out_SI -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/AT48.pindel_out.SI.vcf -e 10
/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/AT48.pindel.out_INV -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/AT48.pindel_out.INV.vcf -e 10
/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/AT48.pindel.out_TD -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/AT48.pindel_out.TD.vcf -e 10
/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/AT48.pindel.out_LI -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/AT48.pindel_out.LI.vcf -e 10
/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/AT48.pindel.out_BP -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/AT48.pindel_out.BP.vcf -e 10
/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/AT48.pindel.out_RP -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/AT48.pindel_out.RP.vcf -e 10


#-----------------------------------------------------
#Identification of mutations
#six types of structural variants（D = deletion、SI = short insertion、INV = inversion、TD = tandem duplication、LI = large insertion、BP = unassigned breakpoints）
bgzip -c $target_ID.pindel_out.D.vcf > $target_ID.pindel_out.D.vcf.gz
tabix -f -p vcf $target_ID.pindel_out.D.vcf.gz

bgzip -c $target_ID.pindel_out.SI.vcf > $target_ID.pindel_out.SI.vcf.gz
tabix -f -p vcf $target_ID.pindel_out.SI.vcf.gz


#filtering out neighborhood mutations 
perl $SCRIPT_DIR/FilteringVcfNeighborSNVs.pindel.pl < $target_ID.pindel_out.D.vcf > $target_ID.pindel_out.D.non_neighbor.vcf
perl $SCRIPT_DIR/FilteringVcfNeighborSNVs.pindel.pl < $target_ID.pindel_out.SI.vcf > $target_ID.pindel_out.SI.non_neighbor.vcf

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
#-----------------------------------------------------



cd $SCRIPT_DIR

module unload miniconda2
module unload vcftools/0.1.15

