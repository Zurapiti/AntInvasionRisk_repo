#!/bin/bash
#SBATCH -p compute
#SBATCH -t 05:00:00
#SBATCH --mem=10G
#SBATCH -c 8
#SBATCH --job-name=ScheduledJobs
#SBATCH --output=ScheduledJob-%j.out
#SBATCH --error=ScheduledJob-%j.err
#SBATCH --mail-user=yazmin.zurita@oist.jp
#SBATCH --mail-type=END,FAIL
##
## author:     Yazmin Zurapiti
## title:      21_deigo_launchR_calcMESS_perSpecies.slurm
## aim:        schedule scripts to run in deigo
## version:    0.1: template for deigo scheduler
##  

#####  directories

## assuming you used 00_deigoLogin.sh to log into deigo:
## call the other script to define all directories
#local_analysis="FALSE"
#source $dir_flash/$proj/"scripts"/"06_defineDirectories.sh"


##### load modules

## to launch an interactive job
#srun -p compute -t 0-6 --mem=100G -c 16 --pty bash
#srun -p compute -t 0-3 --mem=50G -c 16 --pty bash



module use /apps/unit/EconomoU/.modulefiles/amd
## here are some of the most popular loads:
module load geos/3.8.0
module load gdal/3.0.1
module load PROJ/6.1.1
module load R/4.2.1

## by doing this, I kind of log in to a different computer, so
## need to reload directories
user="EconomoU/Zurapiti"
dir_bucket="/bucket/"$user
dir_flash="/flash/"$user
proj="AntInvasionRisk"

cd $dir_flash

local_analysis="FALSE"
source $dir_flash/$proj/"scripts"/"06_defineDirectories.sh"





#R



##### end



