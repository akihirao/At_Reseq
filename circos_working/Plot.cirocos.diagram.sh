#!/bin/bash -i
#Plot.cirocos.diagram.sh
#by HIRAO Akiras

set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)


module load miniconda2

circos -conf circos.AT.mutations.conf -outputfile AT.mutations.png

module unload miniconda2
