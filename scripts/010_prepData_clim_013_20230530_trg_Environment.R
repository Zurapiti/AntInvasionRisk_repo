## 
## author:  Yazmin Zurapiti
## title:   010_prepData_clim_013_20230417_TargetRegion_Subregions-Japan.R
## aim:     crop and mask the so far chosen target regions using natural earth data.
## version: 2022.04.21. I had already successfully produced regions that have no spaces in their names.
##          2022.04.27. Try again to produce New Zealand and the UK.
##                      Updating directories
##
## inputs:            needs to get the worldclim data 
##   1.inputdir: 
##   1.inputfile:
## outputs: 
##   1.outputdir: 
##   1.outputfile: 
## 

# install.packages("rnaturalearth")
# devtools::install_github("ropenscilabs/rnaturalearthdata")
# devtools::install_github("ropensci/rnaturalearthhires") 

library(rnaturalearth)
library(ggplot2)
library(raster)
library(sf)
library(terra) 
library(tidyterra)
library(fasterize)
library(dplyr)


## clear memory
rm(list = ls())

## declare directories
local_analysis=F
if(local_analysis){script06path="./06_defineDirectories.R"} else {
script06path="./AntInvasionRisk/scripts/06_defineDirectories.R"}
source(script06path)

## more directories

### main

#target_area_name=
trg="Japan" #gsub(" ", "-", target_area_name) 
#trg="Taiwan"

## get parameters
args <- commandArgs(trailingOnly = TRUE)   ## allows to pass parameters to the script
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("Not enough parameters", call.=FALSE)
}

trg=args[1]


## reading directories:
dirbd_env=file.path(dir_bucket_data, "EnvironmentalRasters")
dirbd_env_trg=file.path(dirbd_env, trg)  ## subregions' shapes

# writing directories:
dirfd_env_trg=file.path(dir_flash_data, "EnvironmentalRasters", trg)  ## subregions' shapes
if(!dir.exists(dirfd_env_trg)){dir.create(dirfd_env_trg, recursive = T)}

## get the data I need from natural earth

target_region_div=ne_states(country=trg)
## then convert target_region_div object to sf so we can operate on it more easily
target_region_div_sf=st_as_sf(target_region_div)

target_buffered = vect(target_region_div_sf) %>% aggregate() %>% buffer(1000)
ggplot(target_buffered) + theme_bw() + geom_spatvector()

filename=paste0("bufferedSilhouette",".rds")
saveRDS(target_buffered, file = file.path(dirfd_env_trg, filename))

##########

## worldclim data 
wc=list.files(file.path(dirbd_env, "WorldClim_5km"), full.names=TRUE)
worldclim=raster::stack(wc)
laynames=c("bio01", "bio02", "bio03", "bio04", "bio05", "bio06", "bio07", "bio08", "bio09",
           "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")
names(worldclim)=laynames

trg_env = terra::crop(x = rast(worldclim), y = target_buffered, mask=T)   ## crop to the silhouette
writeRaster(trg_env, filename = file.path(dirfd_env_trg, paste0(laynames,".tif")), overwrite = TRUE)

print("owaru, mata ne~")
### ## code ends


