#! /bin/bash

## PREREQUISITES ##
## gfortran (version > 9.0)
## LIBGMXFORT (https://zenodo.org/record/1170728)

#source ~/.bash-gmxfort
module load libgmxfort/default

# debug compile
#gfortran -fcheck=all -lgmxfort -I /home/nikhil/projects/9_srimayee_project/1-hopping/codes/gmxfort/include/ sdf-bmim-cl-ntf2-v1.F08 -o sdf-bmim-cl-ntf2-v1-test.x

# performance compile on local machine
#gfortran -O3 -fopenmp -lgmxfort -I /home/nikhil/projects/9_srimayee_project/1-hopping/codes/gmxfort/include/ sdf-bmim-cl-ntf2-v1.F08 -o sdf-bmim-cl-ntf2-v1.x

# performance compile on server
gfortran -O3 -fopenmp -lgmxfort -I ${GMXFORT_PATH}/include/ sdf-bmim-cl-ntf2-v3.F08 -o sdf-bmim-cl-ntf2-v3.x
