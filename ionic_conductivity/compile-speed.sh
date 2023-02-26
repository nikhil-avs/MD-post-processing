#! /bin/bash

## PREREQUISITES ##
## gfortran (version > 9.0)
## LIBGMXFORT (https://zenodo.org/record/1170728)

# load the prerequisites
module load libgmxfort/default

# performance compile on server
gfortran -fopenmp -lgmxfort -I ${GMXFORT_PATH}/include/ -O3 ic-msdmj-v1.F08 -o ic-msdmj-speed.x
