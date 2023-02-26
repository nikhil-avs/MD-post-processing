#! /bin/bash

## PREREQUISITES ##
## gfortran (version > 9.0)
## LIBGMXFORT (https://zenodo.org/record/1170728)

#source ~/.bash-gmxfort
module load libgmxfort/default

# performance compile on server
gfortran -fopenmp -lgmxfort -I ${GMXFORT_PATH}/include/ -O3 gofr.F08 -o gofr-parallel.x
