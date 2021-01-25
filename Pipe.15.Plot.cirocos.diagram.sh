#!/bin/bash -i
#Plot.cirocos.diagram.sh
#by HIRAO Akiras

set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)

cd circos_working

perl Make.circos.input.pl < ../M2.mutations.full.list.csv


module load miniconda2

circos -conf circos.AT.mutations.conf -outputfile AT.mutations.png

cd ../

module unload miniconda2
