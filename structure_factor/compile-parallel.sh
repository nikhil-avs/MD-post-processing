#! /bin/bash

## PREREQUISITES ##
## gfortran (version > 9.0)
## LIBGMXFORT (https://zenodo.org/record/1170728)

#source ~/.bash-gmxfort
module load libgmxfort/default


if [[ $# -ne 1 ]]; then
    echo "script should be used like -- $0 <EXEC_SUFFIX> "
    exit
fi

suf_exec="$1"
incl_path="${GMXFORT_PATH}/include"

gfortran -fopenmp -lgmxfort -I ${GMXFORT_PATH}/include/ -O3 sofq.F08 -o sofq-${suf_exec}.x
