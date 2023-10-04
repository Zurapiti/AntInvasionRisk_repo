## 
## author:  Yazmin Zurapiti
## title:   06_defineDirectories.R
## aim:     define the directories and source them to other scripts to be consistent
## version: 2022.04.26: this will have the most basic file structure for now
##

## local_analysis=F  example of how to declare the variable

# declare directories in deigo which are super tricky because we read from bucket, but can only write in flash

# so, the most general directories needed:
user="EconomoU/Zurapiti"
dir_bucket=file.path("/bucket", user)
dir_flash=file.path("/flash", user)
proj="AntInvasionRisk"

# then the actual folders for this project
dir_bucket_proj=file.path(dir_bucket, proj)
dir_flash_proj=file.path(dir_flash, proj)

# If I want to use the same declaration but in my own computer where the deigo structure does not exits. 
# I declare a boolean variable called local_analysis (before calling this code) to override the deigo directores with the local ones when T
if(local_analysis==T){
  user="/home/natureza"
  proj="Documents/13.Thesis/AntInvasionRisk"
  dir_bucket_proj=file.path(user,proj,"bucket")
  dir_flash_proj=file.path(user,proj,"flash")
  dir_bucket=dir_bucket_proj
  dir_flash=dir_flash_proj
}


## and we are finally ready to 
## declare the rest of the structure which should be 
## the same across computers/users

## bucket directories:
dir_bucket_scripts=file.path(dir_bucket_proj,"scripts")
dir_bucket_data=file.path(dir_bucket_proj,"data")
dir_bucket_models=file.path(dir_bucket_proj, "models")
dir_bucket_figures=file.path(dir_bucket_proj, "figures")

## write in flash:
dir_flash_scripts=file.path(dir_flash_proj,"scripts")     ##  working/running dir
dir_flash_data=file.path(dir_flash_proj,"data")           ##  global data
dir_flash_models=file.path(dir_flash_proj,"models")       ##  output dir (default same as local data)
dir_flash_figures=file.path(dir_flash_proj,"figures")       ##  figures (default general figures)


