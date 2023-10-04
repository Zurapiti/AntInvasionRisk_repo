## 
## author:  Yazmin Zurapiti 
## title:   vars_032_20230207_BioclimSubsets_PCACorr.R
## aim:     run a correlation analysis over the climatic variables in the space that ants occupy
##          to reduce the set to a more independent set. Hopefully improving our predictions.
## version:   0.1 20220923: 
##                            
## call:    vars_032_20230207_BioclimSubsets_PCACorr.R [save_date]
## inputs: 
##   1.inputdir: 
##   1.inputfile:
## outputs: 
##   1.outputdir: 
##   1.outputfile: 
## 

## libraries

### 2022.09.28: 

###### -------------

library(ggplot2)
library(tidyverse)
library(ggbiplot)
library(corrplot)

## clear memory
rm(list = ls())

# declare directories 
local_analysis=F
if(local_analysis){script06path="./06_defineDirectories.R"} else {
  script06path="./AntInvasionRisk/scripts/06_defineDirectories.R"}
source(script06path)


## collect true parameters:
#trg="Japan"

args <- commandArgs(trailingOnly = TRUE)   ##  parameters to the script
if (length(args)<1) {stop("Not enough parameters", call.=FALSE)}
save_date=args[1]
#trg=args[2]

newdata=paste0("invasionRiskData_5km_",save_date)

dirb_pca=file.path(dir_bucket_data, newdata, "BioclimaticSubsets") 
#dirb_data_trg=file.path(dir_bucket_data,"TargetRegionData",trg)  ## subregions' shapes

dirf_pca=file.path(dir_flash_data, newdata, "BioclimaticSubsets")   ## to save the PCA and cormatrix
if(!dir.exists(dirf_pca)){dir.create(dirf_pca, recursive = T)}
dirf_plots=file.path(dir_flash_figures, newdata, "BioclimaticSubsets")   ## to save the PCA and cormatrix
if(!dir.exists(dirf_plots)){dir.create(dirf_plots, recursive = T)}


print("read data")
#save_date="2023-01-26"  #"2023-01-30"
AntsVals=read.csv(file.path(dirb_pca, paste0("antsVals",".csv")))
AntsVals %>% head()

## this function runs a correlation analysis and PCA to the specified variable set
## in y, and saves to file some objects with the name specified by varnames.
basicCorrPCA=function(y, varnames){
  antsCor=cor(y)   ## correlations will not change regardless subsets
  antsPCA=prcomp(y, center = TRUE, scale. = TRUE)   # run the pc
  
  #antsCor=antsCor %>% class()
  print("save to file")
  
  write.csv(antsCor, file.path(dirf_pca, paste0("antsCor_", varnames ,".csv")))
  saveRDS(antsPCA, file=file.path(dirf_pca, paste0("antsPCA_", varnames ,".rds")))
  write.csv(antsPCA$rotation, file.path(dirf_pca, paste0("antsPCARotation_", varnames ,".csv")))
  
  print("then some plots")
  
  antsPCAX=antsPCA$x %>% as.data.frame()
  
  biplot1 = ggbiplot(antsPCA, color = "lightgrey", alpha = 0.01, varname.size=5, choices=c(1,2)) + theme_bw()
  biplot2 = ggbiplot(antsPCA, color = "lightgrey", alpha = 0.01, varname.size=5, choices=c(2,3)) + theme_bw()
  
  pdf(file.path(dirf_plots, paste0("AntsPCA_", varnames , "_biplot1", ".pdf")))
  print(biplot1)
  dev.off()
  pdf(file.path(dirf_plots, paste0("AntsPCA_", varnames ,"_biplot2", ".pdf")))
  print(biplot2)
  dev.off()
  
  pdf(file.path(dirf_plots, paste0("AntsCor_", varnames, ".pdf")))
  corrplot(antsCor)
  dev.off()
  pdf(file.path(dirf_plots, paste0("AntsPCARotation_", varnames, ".pdf")))
  corrplot(antsPCA$rotation)
  dev.off()
}


#set0=names(y) 


## correlation analysis
print("run correlation and PCA")

## remove lon, lat columns for the analysis
y=AntsVals[,c(-1,-2)]

varnames="allvars"
basicCorrPCA(y=y, varnames=varnames)

## remove interactive species because in our global-wide study these may cause 
## problematic models due to their discontinuities
## this decision resulted from the plots of vars_031_*.R
y = y %>% dplyr::select(-bio08, -bio09, -bio18, -bio19)
varnames="nonInteractive"
basicCorrPCA(y=y, varnames=varnames)


#####     #####


### end
