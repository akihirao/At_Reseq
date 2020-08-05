#!/bin/bash -i
#Run.pindel2vcf.sh
#by HIRAO Akira


set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)

CPU=1

reference_folder=/zfs/Arabidopsis/Reference_v1.1
main_folder=/zfs/Arabidopsis/work/At_Reseq
bwa_folder=$main_folder/bwa_out
work_folder=$main_folder/pindel_out

mkdir -p $work_folder

module load miniconda2

#six types of structural variants（D = deletion、SI = short insertion、INV = inversion、TD = tandem duplication、LI = large insertion、BP = unassigned breakpoints）
/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/AT48.pindel.out_D -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/AT48.pindel_out.D.vcf -e 10
/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/AT48.pindel.out_SI -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/AT48.pindel_out.SI.vcf -e 10
/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/AT48.pindel.out_INV -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/AT48.pindel_out.INV.vcf -e 10
/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/AT48.pindel.out_TD -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/AT48.pindel_out.TD.vcf -e 10
/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/AT48.pindel.out_LI -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/AT48.pindel_out.LI.vcf -e 10
/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/AT48.pindel.out_BP -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/AT48.pindel_out.BP.vcf -e 10
/usr/local/miniconda2/bin/pindel2vcf -p $work_folder/AT48.pindel.out_RP -r $reference_folder/TAIR10.fa -R TAIR10 -d 20101117 -v $work_folder/AT48.pindel_out.RP.vcf -e 10



cd $SCRIPT_DIR

module unload miniconda2

