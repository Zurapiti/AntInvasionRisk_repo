#!/bin/bash
##
## author:     Yazmin Zurapiti
## title:      11_deigo_launchSlurm_RegionRichnessSources.sh
## aim:        schedule scripts to run in deigo
## version:    0.1: template for deigo scheduler
##             2022.01.11: 1.1: Try to run scripts iteratively, and separate run code from copying to bucket
## call script:  sh $dir_flash_scripts/33_20221214_launchSlurm_invaderSpeciesMESSvalue.sh [file] [save_date] [bioset]
##  

## define common directories from script

local_analysis="FALSE"
source ./AntInvasionRisk/scripts/06_defineDirectories.sh

##### main

## define the arguments for the job
slurm_script=11_deigo_launchRscript_RegionRichnessSource.slurm

r_script=$1

save_date=$2
#target_region="Japan"
bioset=$3
#algor='maxent.jar' 
run_date=$4

### iterate over the target_region

#for target_region in Japan # Taiwan Madagascar New-Zealand Cuba
for sp_status in invader #potential 
do
    for model_clean_criteria in fc0 fc1 fc2
    do
        for mask in full mess OR10 both
        do
        
#echo $dir_flash_scripts/$slurm_script $r_script $save_date $bioset $run_date $sp_status $model_clean_criteria "prefecture" $mask

sbatch $dir_flash_scripts/$slurm_script $r_script $save_date $bioset $run_date $sp_status $model_clean_criteria "prefecture" $mask

        done
    done
done




##### end



