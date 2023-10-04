#!/bin/bash
## 
## author:  Yazmin Zurapiti
## title:   06_defineDirectories.sh
## aim:     define pointers to the directory structure in bash to source them to other scripts to be consistent. 
##          this is quite tricky because in deigo we read data from bucket, but can only write on flash, therefore the data structure gets quite confusing. 
##          this script should have an accompanying R version.
## version: 2022.04.26: this will have the most basic file structure for now
##

## local_analysis=$false  example of how to declare the variable


# so, the most general directories needed:
user="EconomoU/Zurapiti"
dir_bucket="/bucket"/$user
dir_flash="/flash"/$user
proj="AntInvasionRisk"

# then the actual folders for this project
dir_bucket_proj=$dir_bucket/$proj
dir_flash_proj=$dir_flash/$proj

# 2022.07.28: I declared some folders to mimic deigo structure in chontali and make it easier to transfer from deigo to chontali and vice versa.
# I declare a boolean variable called local_analysis (before calling this code) to override the deigo directores with the local ones when T
if [[ "$local_analysis" == "TRUE" ]] ; then
user="/home/natureza"
proj="Documents/13.Thesis/AntInvasionRisk"
dir_bucket_proj=$user/$proj/"bucket"
dir_flash_proj=$user/$proj/"flash"
dir_bucket=$dir_bucket_proj
dir_flash=$dir_flash_proj
fi

## and we are finally ready to 
## declare the rest of the structure which should be 
## the same across computers/users

## bucket directories are reading directories:
dir_bucket_scripts=$dir_bucket_proj/"scripts"  ##  working/running dir
dir_bucket_data=$dir_bucket_proj/"data"
dir_bucket_models=$dir_bucket_proj/"models"
dir_bucket_figures=$dir_bucket_proj/"figures"

## flash directories are write directories:
dir_flash_scripts=$dir_flash_proj/"scripts"  ##  working/running dir
dir_flash_data=$dir_flash_proj/"data"        ##  global data
dir_flash_models=$dir_flash_proj/"models"    ##  output dir 
dir_flash_figures=$dir_flash_proj/"figures"      ##  figures (default general figures)


## code ends




