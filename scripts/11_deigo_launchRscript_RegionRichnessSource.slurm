#!/bin/bash
#SBATCH -p compute
#SBATCH -t 03:00:00
#SBATCH --mem=200G
#SBATCH -c 16
#SBATCH --job-name=Run_RegRichSrc_trg
#SBATCH --output=Run_RegRichSrc_trg-%j.out
#SBATCH --mail-user=yazmin.zurita@oist.jp
#SBATCH --mail-type=FAIL
##
## author:     Yazmin Zurapiti
## title:      11_deigo_launchRscript_RegionRichnessSources.slurm
## aim:        schedule scripts to run in deigo
## version:    0.1: template for deigo schedule
##  

##### modules and directories

source ./AntInvasionRisk/scripts/08_deigo_modules.sh

##### main

# get arguments

r_script=$1

save_date=$2
target_region="Japan"
bioset=$3
algor='maxent.jar' 
run_date=$4
sp_status=$5
model_clean_criteria=$6
div_level=$7 
mask=$8

## run script
#echo $dir_flash_scripts/$r_script $save_date $target_region $bioset $algor $run_date $sp_status $model_clean_criteria $div_level $mask

Rscript $dir_flash_scripts/$r_script $save_date $target_region $bioset $algor $run_date $sp_status $model_clean_criteria $div_level $mask






##### end



