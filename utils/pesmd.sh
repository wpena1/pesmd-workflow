#!/bin/bash
input=$1
outputprefix=$2
calc=$3
dim=$4
#conda_dir=$5

if [ -z "$input" ];then
echo "Usage: $0 input_file outputprefix"
exit 1
fi

#this includes putting plumed in the path

#source ${conda_dir}/etc/profile.d/conda.sh
#conda activate plumedgmx
chk=$(which plumed)
if [ -z "$chk" ];then
echo "plumed not found in path: make sure you conda env was installed and activated"
exit 1
else
    plumed --no-mpi pesmd < $input
    if (( $(echo "${calc} == 0" | bc -l) )); #FES case
    then
        if (( $(echo "${dim} == 1" | bc -l) )); # slip model case
        then
        plumed sum_hills  --mintozero --bin 100 --min 0.0 --max 25.0 --hills ${outputprefix}_pesmd.hills --outfile ${outputprefix}.fes.dat
        fi
        if (( $(echo "${dim} == 2" | bc -l) )); # catch model case
        then
        plumed sum_hills  --mintozero --bin 300,300 --min -5.0,-5.0 --max 31.0,12.0 --hills ${outputprefix}_pesmd.hills --idw d1.x --kt 2.49 --outfile ${outputprefix}.fesx.dat
        plumed sum_hills  --mintozero --bin 300,300 --min -5.0,-5.0 --max 31.0,12.0 --hills ${outputprefix}_pesmd.hills --idw d1.y --kt 2.49 --outfile ${outputprefix}.fesy.dat
        plumed sum_hills  --mintozero --bin 300,300 --min -5.0,-5.0 --max 31.0,12.0 --hills ${outputprefix}_pesmd.hills --outfile ${outputprefix}.fes.dat
        fi
    fi
fi