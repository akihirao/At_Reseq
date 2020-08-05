#!/bin/bash -i
#Run.pindel.sh
#by HIRAO Akira


set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)

CPU=16

reference_folder=/zfs/Arabidopsis/Reference_v1.1
main_folder=/zfs/Arabidopsis/work/At_Reseq
bwa_folder=$main_folder/bwa_out
work_folder=$main_folder/pindel_out

mkdir -p $work_folder

module load miniconda2

#Identifying structural variants with using pindelコマンド
#six types of structural variants（D = deletion、SI = short insertion、INV = inversion、TD = tandem duplication、LI = large insertion、BP = unassigned breakpoints）
/usr/local/miniconda2/bin/pindel -T $CPU -f $reference_folder/TAIR10.fa -i $SCRIPT_DIR/AT48.pindel.config.txt -c ALL -o $work_folder/AT48.pindel.out

cd $SCRIPT_DIR

module load miniconda2

