## 
## author:  Yazmin Zurapiti
## title:   14_GetSpeciesOccurenceData_RemoveSpaces.sh
## aim:     record all copies I do deigo-chontali or between deigo e.g. Jamie-Jazz
## version: 1.0: 2022.04.04 onwards, copy data to deigo
##               2022.04.26 I am editting this to change some structure in the names of files and folders to make it more comprehensible
##               2022.05.16 Jamie updated data as there was a bug and we removed "indoor introduced"
##               2022.07.28 Rerunning whole analysis, thus running this again
## run:     run from chontali
##

deigo="raven-zurapiti@deigo.oist.jp"
ssh $deigo

user="EconomoU/Zurapiti"
dir_bucket="/bucket/"$user
dir_flash="/flash/"$user
proj="AntInvasionRisk"

cd $dir_flash

# to use the script that defines directories
local_analysis="FALSE"
## call the other script to define all directories
source ./AntInvasionRisk/scripts/06_defineDirectories.sh  


####


## enter a computing node
#srun -p compute -t 0-6 --mem=200G -c 16 --pty bash  ## heavy computation
srun -p compute -t 0-1 --mem=20G -c 4 --pty bash

## call all modules and directories needed:
source ./AntInvasionRisk/scripts/08_deigo_modules.sh

## crop the rasters for the target region
trg=Japan
Rscript $dir_flash_scripts/010_prepData_clim_013_20230530_trg_Environment.R $trg

## farther down the line Japan analysis needs a very specific structure of subregions:
## only Japan needs this script as we have custom subdivisions
Rscript $dir_flash_scripts/010_prepData_clim_013_20230530_Japan_Subregions.R 

## exit the computing node to synchornize to bucket:
exit

## syncronize to bucket
dirf=$dir_flash_data/"EnvironmentalRasters"
#ls $dirf/*/*
dirb=$dir_bucket_data/"EnvironmentalRasters"
#ls $dirb/*/*
rsync -crt $dirf/ $dirb
#ls $dirb/*/*

clear

#ls $dirf/*/*
#ls $dirb/*/*

rm -r $dirf/*        




##############################################
###### species occurrence data

save_date="2023-05-31"


cd $dir_flash

## split native, invader, potential of each target_region

srun -p compute -t 0-3 --mem=20G -c 8 --pty bash

source ./AntInvasionRisk/scripts/08_deigo_modules.sh


reg=Japan
Rscript $dir_flash_scripts/010_prepData_ants_026_20230131_TargetRegion_SppStatusSeparation.R $save_date $reg "country"


dirf=$dir_flash_data  
#ls $dirf/*/*
dirb=$dir_bucket_data
rsync -crt $dirf/ $dirb
#ls $dirb/*/*
rm -r $dirf/*        

## all done for the species data


##############################################
###### variable selection

############################
#save_date="2023-05-31"

## enter a computing node: 
srun -p compute -t 0-6 --mem=200G -c 16 --pty bash
## heare asking for a lot of memory to avoid reaching the limit

source ./AntInvasionRisk/scripts/08_deigo_modules.sh


##### PART 1
## prepare data for PCA-correlation analyses and VIF. 
## data is large, so I first create it and then run the analysis to avoid crashing a job 
## this scripts produces data and plots:

Rscript $dir_flash_scripts/010_prepData_vars_031_20230208_BioclimSubsets_PrepData.R $save_date


## exit the computing node:
exit

## syncronize data to bucket
dirf=$dir_flash_data  
#ls $dirf/*/*
dirb=$dir_bucket_data
rsync -crt $dirf/ $dirb
#ls $dirb/*/*
rm -r $dirf/*        


##### PART 2
## creating the bioclimatic variable sets

## another computing node: (this time we don't need that many resources)
srun -p compute -t 0-3 --mem=20G -c 8 --pty bash

source ./AntInvasionRisk/scripts/08_deigo_modules.sh

Rscript $dir_flash_scripts/010_prepData_vars_032_20230207_BioclimSubsets_PCACorr.R $save_date


exit

## syncronize to bucket
dirf=$dir_flash_data  
#ls $dirf/*/*
dirb=$dir_bucket_data
rsync -crt $dirf/ $dirb
#ls $dirb/*/*
rm -r $dirf/*        


## syncronize figures to bucket
dirf=$dir_flash_figures 
#ls $dirf/*/*
dirb=$dir_bucket_figures
rsync -crt $dirf/ $dirb
#ls $dirb/*/*
rm -r $dirf/*        


## copy the variable set file (from chontali) to bucket


###### 
exit

## code ends








