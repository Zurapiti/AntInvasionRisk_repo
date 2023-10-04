#!/bin/bash
##
## author:     Yazmin Zurapiti
## title:      11_deigo_launchSlurm_*.sh
## aim:        schedule scripts to run in deigo
## version:    0.1: template for deigo scheduler
##             2022.01.11: 1.1: Try to run scripts iteratively, and separate run code from copying to bucket
## call script:  sh $dir_flash_scripts/11_deigo_launchSlurm_*.sh [r_script] [save_date] [target_region] [bioset] [file_species]
##  

## define common directories from script


local_analysis="FALSE"
source ./AntInvasionRisk/scripts/06_defineDirectories.sh

##### main

## define the arguments for the job
slurm_script=11_deigo_2023-04-21_launchRscript_SpeciesComputation.slurm

r_script=$1

save_date=$2
target_region=$3
bioset=$4
algor=$5
run_date=$6
sp_status=$7
file_species=$8

### iterate over the species in the species file of the target_region

while read -r line; do
#    echo -e "$line\n"
#    remove quotes ("...") 
this_species=$(sed -e 's/^"//' -e 's/"$//' <<<"$line")
echo $this_species

#echo $dir_flash_scripts/$slurm_script $r_script $save_date $target_region $bioset $this_species

## batch the job
sbatch $dir_flash_scripts/$slurm_script $r_script $save_date $target_region $bioset $algor $run_date $sp_status $this_species


done <$file_species


##### end



