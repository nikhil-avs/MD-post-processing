#! /bin/bash

# LAST MODIFIED ON 13 FEB 2020

## PREREQUISITES ##
## gfortran (version > 9.0)
## LIBGMXFORT (https://zenodo.org/record/1170728)

#source ~/.bash-gmxfort
module load libgmxfort/default

# check if you have two command line arguments 
if [[ $# -ne 3 ]]; then
    echo "script should be used like -- nohup $0 <EXEC_NAME> <INPFILE> <NUM_THREADS> &"
    exit
fi

# set some names using command line arguments
exec_name=$1
inpfile=$2
export OMP_NUM_THREADS=$3

# assuming the input file to be called sofq.inp extract suffix name
suffix=$(tail -n1 ${inpfile})
echo "USING SUFFIX: "${suffix}

# check if the output directories exist already
if [[ -d sabq-${suffix} ]]; then
    echo "output directory already exists! Remove it and rerun"
    exit
fi

if [[ -d sabs-${suffix} ]]; then
    echo "output directory already exists! Remove it and rerun"
    exit
fi

if [[ -d sqpart-${suffix} ]]; then
    echo "output directory already exists! Remove it and rerun"
    exit
fi

if [[ -d sspart-${suffix} ]]; then
    echo "output directory already exists! Remove it and rerun"
    exit
fi

if [[ -d ffs-${suffix} ]]; then
    echo "output directory already exists! Remove it and rerun"
    exit
fi

# create output directories
mkdir sabq-${suffix}  sabs-${suffix} ffs-${suffix} sqpart-${suffix}  sspart-${suffix}

# name the log files
logfile="sofq-"${suffix}".log"
junkfile="junk-"${suffix}".dat"

# write essential run info to log file
echo " -- STRUCTURE FACTOR CODE V1.0 -- log file" > ${logfile}

tempv=$(date)
echo "RUN-INFO) Program started on : " ${tempv} >> ${logfile}
 
tempv=$(whoami)
echo "RUN-INFO) Program run on     : " ${tempv} >> ${logfile}

tempv=$(hostname)
echo "RUN-INFO) Program run by     : " ${tempv} >> ${logfile}

tempv=$(pwd)
echo "RUN-INFO) Program run at     : " ${tempv} >> ${logfile}

# finally run the damn thing!
"./"${exec_name} ${inpfile} 1>> ${logfile} 2> ${junkfile}
