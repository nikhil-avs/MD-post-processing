#! /bin/bash

## PREREQUISITES ##
## gfortran (version > 9.0)
## LIBGMXFORT (https://zenodo.org/record/1170728)

#source ~/.bash-gmxfort
module load libgmxfort/default

if [[ -d data-inter ]]; then
    echo "output directory already exists!"
    exit
fi

if [[ -d data-agnostic ]]; then
    echo "output directory already exists!"
    exit
fi

mkdir data-inter data-agnostic

echo " -- GOFR CODE V1.0 -- log file" 

tempv=$(date)
echo "RUN-INFO) Program started on : " ${tempv}  

tempv=$(whoami)
echo "RUN-INFO) Program run on     : " ${tempv}  

tempv=$(hostname)
echo "RUN-INFO) Program run by     : " ${tempv}  

tempv=$(pwd)
echo "RUN-INFO) Program run at     : " ${tempv}  

./gofr-parallel.x gofr.inp 2> junk.dat
