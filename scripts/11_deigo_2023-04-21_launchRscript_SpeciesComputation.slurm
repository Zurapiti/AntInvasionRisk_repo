#!/bin/bash
#SBATCH -p compute
#SBATCH -t 03:00:00
#SBATCH --mem=100G
#SBATCH -c 16
#SBATCH --job-name=Run_AntInv
#SBATCH --output=Run_AntInv-%j.out
#SBATCH --mail-user=yazmin.zurita@oist.jp
#SBATCH --mail-type=FAIL
##
## author:     Yazmin Zurapiti
## title:      11_deigo_launchRscript_Computation.slurm
## aim:        schedule scripts to run in deigo
## version:    0.1: template for deigo scheduler
##  

##### modules and directories

source ./AntInvasionRisk/scripts/08_deigo_modules.sh

##### main

# get arguments

r_script=$1

save_date=$2
target_region=$3
bioset=$4
algor=$5
run_date=$6
sp_status=$7
this_species=$8


## run script

echo $dir_flash_scripts/$r_script $save_date $target_region $bioset $algor $run_date $sp_status $this_species

ls -l /bucket/EconomoU/Zurapiti/AntInvasionRisk/data/invasionRiskData_5km_2023-02-09/model_inputs/$this_species

Rscript $dir_flash_scripts/$r_script $save_date $target_region $bioset $algor $run_date $sp_status $this_species


##### end



