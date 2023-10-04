## 
## author:  Yazmin Zurapiti 
## title:   vars_031_20230206_BioclimSubsets_PrepData.R
## aim:     run a correlation analysis over the climatic variables in the space that ants occupy
##          to reduce the set to a more independent set. Hopefully improving our predictions.
## version:   0.1 20220923: 
##                            
## call:    Rscript vars_031_20230206_BioclimSubsets_PrepData.R [save_date]
## inputs: 
##   1.inputdir: 
##   1.inputfile:
## outputs: 
##   1.outputdir: 
##   1.outputfile: 
## 
###### -------------

## libraries
library(stringr)  ## for leading zeros
library(terra)
library(sf)
library(sp)
library(tidyverse)
library(rasterVis)
#library(maptools)
library(viridis) 

## clear memory
rm(list = ls())

# declare directories 
local_analysis=F
if(local_analysis){script06path="./06_defineDirectories.R"} else {
  script06path="./AntInvasionRisk/scripts/06_defineDirectories.R"}
source(script06path)

# save_date="2023-01-16"
## collect true parameters:
args <- commandArgs(trailingOnly = TRUE)   ##  parameters to the script
if (length(args)<1) {stop("Not enough parameters", call.=FALSE)}
save_date=args[1]

newdata=paste0("invasionRiskData_5km_",save_date)

dirb_ants=file.path(dir_bucket_data, newdata, "model_inputs")
dirb_env=file.path(dir_bucket_data, "EnvironmentalRasters")

dirf_pca=file.path(dir_flash_data, newdata, "BioclimaticSubsets")   ## to save the PCA and cormatrix
if(!dir.exists(dirf_pca)){dir.create(dirf_pca, recursive = T)}
dirf_plots=file.path(dir_flash_figures, newdata, "BioclimaticSubsets")
if(!dir.exists(dirf_plots)){dir.create(dirf_plots, recursive=T)}


### load climatic data: 
print("read data")
laynames=c("bio01", "bio02", "bio03", "bio04", "bio05", "bio06", "bio07", "bio08", "bio09",
"bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")
source_area_name="World"
dir_source_area=file.path(dirb_env,"WorldClim_5km")
sc=list.files(dir_source_area, full.names=TRUE)
source_area=rast(sc)
names(source_area) = laynames

## load ants data:
dir_ants_points=file.path(dirb_ants, "exotics_gabidb_list_thin.rds")
ant_points=readRDS(dir_ants_points)

########### main

print("rbind all points, make them one single huge data frame")
all_points=do.call(rbind, ant_points) %>% as.data.frame() %>% select("lon_opt", "lat_opt") %>% distinct()

##
if(local_analysis==TRUE){
  n1=sample.int(n=nrow(all_points), size=50, replace=F)
  all_points=all_points[n1,]
}
###

print("create big buffer")

print("sf solution to the previously done by geos all species buffer")
x= all_points %>% st_as_sf(coords=c("lon_opt", "lat_opt"))
sp_buff3=st_buffer(x, dist=2) ## independent tiny polygons
sp_buff4=st_combine(sp_buff3)  ## make them a single multipolygon
sp_buff5=st_buffer(sp_buff4, dist=1) # re-buffer the multipolygon to have a single shape

# plot(x, pch=16, cex=0.3)
# plot(sp_buff3, add=T, border="green")
# plot(sp_buff4, add=T, border="red")
# plot(sp_buff5, add=T, border="blue")

print("put things back to terra objects")

buf6=vect(sp_buff5) 

print("saving shape to file")
writeVector(buf6, file.path(dirf_pca, "allAnts_buffer.shp"), overwrite=T)

print("rasterize and mask world with it")
buf7=rasterize(buf6, source_area[[1]])
writeRaster(buf7, file.path(dirf_pca, "allAnts_raster.tif"), overwrite=T)

###############################
## plot such buffer/raster over the worldclim data 
## to see if we hit discontinuities 

print("plotting buffer over worldclim layers to check for discontinuities")


## basic terra plotting approach
for (lay in 1:19) {
  pdf(file=file.path(dirf_plots, paste0("bioclim", str_pad(lay, 2, pad = "0"), ".pdf")))
  plot(source_area[[lay]], ylim=c(-90, 90))
  plot(buf6, add=T, border="grey", alpha=125)
  dev.off()
}

## using levelplot to create the figures
for (lay in 1:19) {
  r=source_area[[lay]]
  u=minmax(r) 
  u = t(minmax(r)) %>% as.data.frame()# %>% class()
  
  ### this plots look nice:
  p1=levelplot(r, margin=FALSE,                       
               colorkey=list(space='bottom'),       
               col.regions=viridis, at=seq(u$min, u$max, len=201))#
  
  p2=levelplot(buf7, margin=FALSE, colorkey=F,
               col.regions='darkgray', alpha.regions=0.5)
  
  #p1
  #p2
  pdf(file=file.path(dirf_plots, paste0("bioclim", str_pad(lay, 2, pad = "0"), "_levelplot", ".pdf")))
  print(p1+p2)
  dev.off()
  png(filename=file.path(dirf_plots, paste0("bioclim",  str_pad(lay, 2, pad = "0"), "_levelplot", ".png")))
  print(p1+p2)
  dev.off()
  
}


###############################
## then extract the actual values within those cells to run PCA and VIF analysis

## then I mask the world data with buf4, so we keep only where both rasters intersect
buf8=mask(source_area, buf7)  ## this is a rasterBrick
# plot(buf8)

print("extract values of non-NA cells")
## get the RasterStack data as a matrix with
# raster cells as rows and raster layers as columns
antsValuesMatrix=as.data.frame(buf8, xy = TRUE, na.rm=T)
#names(antsValuesMatrix) = c("lon","lat", names(buf8))

antsValuesMatrix %>% dim()   ## ok, looks good
antsValuesMatrix %>% head()   ## ok, looks good

print("write such data to file")
## PCA and VIF and all those analyses should be ran on this data
write.csv(antsValuesMatrix, file=file.path(dirf_pca, paste0("antsVals",".csv")), row.names = F)


print("reached the end of this script and all looks good")
print("bye bye")
#####     #####


### end
