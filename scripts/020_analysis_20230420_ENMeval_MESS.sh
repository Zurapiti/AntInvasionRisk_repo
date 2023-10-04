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

##################
## srun -p compute -t 0-6 --mem=200G -c 16 --pty bash
## source ./AntInvasionRisk/scripts/08_deigo_modules.sh

#---# this code launches the species-level jobs in deigo

#sp_sh_script=11_deigo_2023-04-21_launchSlurm_SpeciesComputation.sh
#sp_slurm_script=11_deigo_2023-04-21_launchRscript_SpeciesComputation.slurm

#---# this code launches the region-level jobs in deigo


#---# SDM and MESS calculations


######################
## testing files: 
## invader test file:    inv_test.csv
## potential test file:   pot_test.csv

sp_sh_script=11_deigo_2023-04-21_launchSlurm_SpeciesComputation.sh
r_script=020_analysis_20230428_ENMeval-MESS.R

save_date="2023-05-31"
target_region="Japan"
bioset="bioset05"
algor="maxent.jar"

#### TEST
run_date="2023-05-31_test"   ## testing date

sp_status="invader"
file_species=$dir_bucket_data/"EnvironmentalRasters"/$target_region/inv_test.csv
ls $file_species

sh $dir_flash_scripts/$sp_sh_script $r_script $save_date $target_region $bioset $algor $run_date $sp_status $file_species

sp_status="potential"
file_species=$dir_bucket_data/"EnvironmentalRasters"/$target_region/pot_test.csv
ls $file_species

sh $dir_flash_scripts/$sp_sh_script $r_script $save_date $target_region $bioset $algor $run_date $sp_status $file_species


###### READY TO RUN
run_date="2023-06-23"   ## testing date

sp_status="invader"
file_species=$dir_bucket_data/"EnvironmentalRasters"/$target_region/$sp_status".csv"
ls $file_species

sh $dir_flash_scripts/$sp_sh_script $r_script $save_date $target_region $bioset $algor $run_date $sp_status $file_species

##
sp_status="potential"
file_species=$dir_bucket_data/"EnvironmentalRasters"/$target_region/$sp_status".csv"
#file_species=$dir_bucket_data/"EnvironmentalRasters"/$target_region/pot_test.csv
ls $file_species

sh $dir_flash_scripts/$sp_sh_script $r_script $save_date $target_region $bioset $algor $run_date $sp_status $file_species

## sinchronize to bucket here!


######### COUNTRY LEVEL 

#### This needs to be done "by hand" as we need to analyse the plots and recovered patters
srun -p compute -t 0-6 --mem=200G -c 16 --pty bash
source ./AntInvasionRisk/scripts/08_deigo_modules.sh

#### analyse relations between "good transfers" and transfer/source metrics with
save_date="2023-05-31"
trg="Japan"
bioset="bioset05"
algor='maxent.jar' #"maxnet"
run_date="2023-06-23"

### how good our transfers are
r_script=020_analysis_20230620_transferabilityEvaluation.R
Rscript $dir_flash_scripts/$r_script $save_date $trg $bioset $algor $run_date

### some sort of "predictions" of what to expect
#r_script=trg_046_20230509_H0.R   
#sp_status="invader" # "potential" #"prediction" # 
#Rscript $dir_flash_scripts/$r_script $save_date $trg $bioset $algor $run_date $sp_status 



################
## richness map, clustering analysis and source regions

save_date="2023-05-31"
target_region="Japan"
bioset="bioset05"
algor='maxent.jar' 
run_date="2023-06-23"

## launch all combinations of the filters (for the invaders) through a script to analyse and decide what to use for exotics
sp_sh_script=11_deigo_launchSlurm_RegionRichnessSource.sh
r_script=trg_049_20230725_RichnessMap_Source.R
#trg_048_20230529_RichnessMap_plottingAndOtherApproaches_SourceFromSDM.R

sh $dir_flash_scripts/$sp_sh_script $r_script $save_date $bioset $run_date


## ## after we decided what were our best combinations:
slurm_script=11_deigo_launchRscript_RegionRichnessSource.slurm
r_script=trg_049_20230725_RichnessMap_Source.R

sp_status="potential" #"invader" # 
model_clean_criteria="fc2" ##  (occ-cbi combination)
div_level="prefecture"  ## cluster  ## the administrative level 
mask="both"  ## what mask to use for the prefecture prediction
sbatch $dir_flash_scripts/$slurm_script $r_script $save_date $bioset $run_date $sp_status $model_clean_criteria $div_level $mask

## hierarchical analysis
r_script=trg_052_20230721_HirerarchichalClusterAnalysis.R
sbatch $dir_flash_scripts/$slurm_script $r_script $save_date $bioset $run_date $sp_status $model_clean_criteria $div_level $mask

### rerun RichnessMap_Source for clusters 
r_script=trg_049_20230725_RichnessMap_Source.R
div_level="cluster"  ## cluster  ## the administrative level 
sbatch $dir_flash_scripts/$slurm_script $r_script $save_date $bioset $run_date $sp_status $model_clean_criteria $div_level $mask



#####

######

save_date="2023-05-31"
run_date="2023-06-23"

## syncronize to bucket
dir_sync=run-$run_date\_data-$save_date\_varset-bioset05/Japan

## models
dirf=$dir_flash_models/$dir_sync
dirb=$dir_bucket_models/$dir_sync
mkdir -vp $dirb
rsync -crt $dirf/ $dirb
#ls $dirb/*/*
#rm -r $dirf/* 

## figures 
dirf=$dir_flash_figures/$dir_sync
dirb=$dir_bucket_figures/$dir_sync
mkdir -vp $dirb
rsync -crt $dirf/ $dirb
#ls $dirb/*/*
#rm -r $dirf/* 


###### 
exit

## code ends








