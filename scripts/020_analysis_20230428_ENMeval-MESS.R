## 
## author:  Yazmin Zurapiti
## title:   020_analysis_20230428_ENMeval-MESS.R
## aim:     this script calculates the ENMEval procedure, selects the "best" model, 
##          calculates the MESS with the variables and reports metrics that will be used later. 
## version:   0.1 2022.02.18: merge suitability and mess analysis
##                2022.03.23: removed full=T to get only one layer of mess
##            1.1 2022.04.28, fix directories, remove duplicated coordinate points in data
##            2.0 2022.08.01, adjust directories after 07.28 archiving
##            2.1 2022.08.01, use thinned points directly, region_admin_rank argument is no longer needed
##            2.3, source=world, target!=world version
##            3.0 2023.02.14, a bunch of changes in the data have happened, so preparing all this again.
##                            
## call:  lauch from /flash/EconomoU/Zurapiti as there is where the R packages are installed
##        Rscript 020_analysis_20230428_ENMeval-MESS.R [save_date] [trg] [bioset] [this_species] [algor] 
## inputs: 
##   1.inputdir: 
##   1.inputfile:
## outputs: 
##   1.outputdir: 
##   1.outputfile: 
## 

# Load packages -- the order here is important because some pkg functions overwrite others.
options(java.parameters = "-Xmx8000m")
## Picked up _JAVA_OPTIONS: "-Xss2560k"  ## this is what deigo "picks"

# install.packages("remotes")
# remotes::install_github("sjevelazco/flexsdm")

library(sp)
library(sf)
library(terra)
library(raster)
library(rnaturalearth)
library(ggplot2)
library(tidyterra)
library(stringi)  ## to make a matrix out of a list of uneven number of rows 
library(stringr)  ## for leading zeros #str_pad(lay, 2, pad = "0")
library(purrr)
library(rlist)
library(dismo)
library(ENMeval)
library(dplyr)
library(rJava)

## clear memory
rm(list = ls())
set.seed(1)

# testing parameters
save_date="2023-02-09"
trg="Japan"
bioset="bioset05"
algor='maxent.jar' #"maxnet"
run_date="2023-04-19"
sp_status="invader" # "potential" # 
this_species="Paratrechina_longicornis" 
# "Myrmica_rubra" 
# "Brachyponera_obscurans" # "Leptogenys_pavesii" # "Tapinoma_sessile" #
# "Tetramorium_sericeiventre" # 
# "Pheidole_megacephala" #"Tapinoma_melanocephalum" # 

## get parameters
args <- commandArgs(trailingOnly = TRUE)   ## allows to pass parameters to the script
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("Not enough parameters", call.=FALSE)
}

save_date=args[1]
trg=args[2]
bioset=args[3]
algor=args[4]
run_date=args[5]
sp_status=args[6]
this_species=args[7]

sessionInfo()

print("parameters:")
print(paste("data date =", save_date, ", target_region =", trg, 
            ", bioset =", bioset, ", algorithm =", algor, ", running date =", run_date,
            ", species type =", sp_status, ", this_species =", this_species, sep=" " ))

# declare directories 
print("setting up directories")

local_analysis=F
if(local_analysis){script06path="./06_defineDirectories.R"} else {
  script06path="./AntInvasionRisk/scripts/06_defineDirectories.R"}
source(script06path)

newdata=paste0("invasionRiskData_5km_", save_date)  ## dir_data
dirbd_ants=file.path(dir_bucket_data, newdata, "model_inputs")
dirbd_varsets=file.path(dir_bucket_data, newdata, "BioclimaticSubsets")

dirbd_env=file.path(dir_bucket_data, "EnvironmentalRasters")
dirbd_env_src=file.path(dirbd_env, "WorldClim_5km")
dirbd_env_trg=file.path(dirbd_env, trg)

datavar=paste0("run-", run_date, "_data-", save_date, "_varset-", bioset) ## dir_models
dirbm_datvar_trg=file.path(dir_bucket_models, datavar, trg)
dirbm_datvar_tsp=file.path(dirbm_datvar_trg, sp_status, this_species)   ## species directory (only declared here)

## save directories:
dirfm_datvar_trg=file.path(dir_flash_models, datavar, trg)
dirfm_datvar_tsp=file.path(dirfm_datvar_trg, sp_status, this_species)   ## species directory (only declared here)
if(!dir.exists(dirfm_datvar_tsp)){dir.create(dirfm_datvar_tsp, recursive = T)}

dirff_datvar_trg=file.path(dir_flash_figures, datavar, trg)
dirff_datvar_tsp=file.path(dirff_datvar_trg, sp_status, this_species)   ## species directory (only declared here)
if(!dir.exists(dirff_datvar_tsp)){dir.create(dirff_datvar_tsp, recursive = T)}

## functions

## test that the points I am trying to use don't "land" on NA cells
pointsNA = function(area, points) {
  spp_env = terra::extract(x=area, y=points, ID=F)
  naIndx = spp_env %>% rowSums() %>% is.na() %>% which()
  return(naIndx)
}

## calculate the area of a raster
RastArea = function(rast){
  csize=cellSize(rast, mask=T, unit="km")
  buff_area=as.data.frame(csize) %>% select(area) %>% sum()
  return(buff_area)
}

## this function helps to save pdf and png plots (from an object) at the same time
plot_pdf_png = function(the_plot, plot_dir, plot_name, pointsize=12, pdf_width=7){
  fname=paste0(plot_name,".pdf")
  pdf(file=file.path(plot_dir, fname), pointsize=pointsize, width = pdf_width)
  print(the_plot)
  dev.off()
  
  fname=paste0(plot_name,".png")
  png(file=file.path(plot_dir, fname), pointsize=pointsize, width = round(pdf_width*68.6))
  print(the_plot)
  dev.off()
}

# 2022-04-20: function from Jamie's package to calculate a threshold.
calc.10p.trainThresh <- function(pred.train) {
  n <- length(pred.train)
  if(n < 10) {pct90.train <- floor(n * 0.9)} else {pct90.train <- ceiling(n * 0.9)}
  pct10.train.thr <- rev(sort(pred.train))[pct90.train]
  return(pct10.train.thr)
}

## make a data frame out of the lambda text output
lambdasDF=function(mlambdas){
  mlambdas=mlambdas  %>%  strsplit(",") %>% stri_list2matrix(byrow = T) %>% as.data.frame()
  mlambdas=mlambdas[1:(nrow(mlambdas)-4) , ]
  m2lambdas=lapply(mlambdas[,2:4], as.numeric) %>% simplify2array()
  mlambdas[,2:4]=m2lambdas
  return(mlambdas)
}

## test if the models have at least one lambda different from zero
anyPossitiveLambda = function(mlambdas){
  mlambdas=lambdasDF(mlambdas)
  sumLambdas=mlambdas$V2 %>% as.numeric() %>% abs() %>% sum()
  return(sumLambdas>0)
}

### don't understand why, but some mess measures are Inf,
## so we make then NaN
Inf_to_NaN=function(x){
  x[x[]==Inf] = NaN
  return(x)
}




# 
# ## to prepare a binary map of where the species already occurs in the trg
# emptyPred = function(rast){
#   pred=rast[[1]]
#   pred[!is.na(pred)] = 0    ## prediction raster; non NA cells get value zero
#   names(pred) = "obs"
#   return(pred)
# }
# 
# 
# ## threshold a raster of the terra class
# ## specifically to threshold MESS rasters which will be binary rasters in the end
# thresholdMESSRaster=function(x, thr){
#   x[x[]==Inf] = NaN
#   rasminmax=terra::minmax(x)
#   m <- c(Inf, Inf, NaN,
#          rasminmax[1,], thr, 0,
#          thr, rasminmax[2,], 1)
#   rclmat <- matrix(m, ncol=3, byrow=TRUE)
#   rc1 <- terra::classify(x, rclmat, include.lowest=TRUE)
#   return(rc1)
# }

##### main
print("data preparation starts")

wcl=ne_coastline() %>% vect()

laynames=c("bio01", "bio02", "bio03", "bio04", "bio05", "bio06", "bio07", "bio08", "bio09",
           "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")

## Select the variable set
vars_sets=read.csv(file.path(dirbd_varsets, "BioclimSetsTest_HandMade.csv") )
biovars=vars_sets %>% dplyr::select(all_of(bioset)) %>% t() #%>% paste()
biovars=biovars[biovars != ""]

## Environmental data for target_region:
print("target")
tg_files=list.files(dirbd_env_trg, pattern=".tif", full.names=TRUE)
target_area=raster::stack(tg_files)    ## target area needs to be a Raster* for dismo
names(target_area) = laynames
target_area=subset(target_area, subset=biovars)

## we also need the buffered silhouette to remove the points in the target area
trg_div=readRDS(file = file.path(dirbd_env_trg, "Subregions.rds")) ## sf class
trg_silh = trg_div %>% vect() %>% aggregate() ## terra class, for plotting
#shp_file="bufferedSilhouette.rds"
#target_silh=readRDS(file.path(dirbd_env_trg, shp_file))

## Environmental data for source_region:
print("source")
sp_rg=list.files(file.path(dirbd_ants, this_species), pattern=".tif", full.names=T)
print(sp_rg)
total_area = rast(sp_rg)
names(total_area) = laynames
total_area=subset(total_area, subset=biovars)  ## the species' shape
print("subsetted to biovars")
total_sqm=RastArea(rast=total_area[[1]])
print("area calculated")

rastIntersects=relate(total_area[[1]], rast(target_area[[1]]), relation="intersects")

## species polygonal shape
spp_silh=readRDS(file.path(dirbd_ants, "exotics_species_hull_list.rds"))  ## all species polygons
spp_silh = spp_silh[[this_species]]
polygon_shape=names(spp_silh)[1]
spp_silh = spp_silh[[1]] 

## species occurrence data
print("species occurrence data")
spp_data = readRDS(file.path(dirbd_ants, "exotics_gabidb_list_thin.rds"))   ## all species 
spp_data = spp_data[[this_species]] %>% as.data.frame() # this particular species
spp_coords_df = spp_data %>% select("lon_opt", "lat_opt") %>% rename(lon=lon_opt, lat=lat_opt) # data.frame
spp_tot_occ = nrow(spp_coords_df)
spp_occ  = spp_coords_df %>% vect(geom=c("lon", "lat"))  # terra object
crs(spp_occ) = "+proj=longlat +datum=WGS84 +no_defs"

occindexNA=pointsNA(total_area, spp_occ)
remove_occpoint = (length(occindexNA) > 0) 
if(remove_occpoint){
  removed_coords = spp_coords_df[occindexNA,] 
  spp_occ = spp_occ[-occindexNA,]
  report_NA_occ = data.frame(occindexNA, removed_coords)
  }
spp_noNA_occ = nrow(spp_occ)

## background:
spp_bgr_df = read.csv(file.path(dirbd_ants, this_species, "bg_pts.csv")) %>% rename(lon=x, lat=y) 
spp_tot_bgr = nrow(spp_bgr_df)
spp_bgr=spp_bgr_df %>% vect(geom=c("lon", "lat"))

bgrindexNA=pointsNA(total_area, spp_bgr)
remove_bgrpoint = (length(bgrindexNA) > 0) 
if(remove_bgrpoint){
  removed_bgr = spp_bgr_df[bgrindexNA,]
  spp_bgr = spp_bgr[-bgrindexNA,]
  report_NA_bgr = data.frame(bgrindexNA, removed_bgr)
} 
spp_noNA_bgr = nrow(spp_bgr)

################
#### points in target region?: 
print("separate target from source")

## define source_region by erasing the target_region
spp_occ_source  = terra::erase(x=spp_occ, y=trg_silh)
spp_bgr_source  = erase(x=spp_bgr, y=trg_silh)
spp_scr_occ=nrow(spp_occ_source)  
spp_scr_bgr=nrow(spp_bgr_source)

## define test_data with the points inside the target_region
spp_occ_target  = crop(x=spp_occ, y=trg_silh)
spp_bgr_target  = crop(x=spp_bgr, y=trg_silh)
spp_trg_occ=nrow(spp_occ_target)  
spp_trg_bgr=nrow(spp_bgr_target)

if(rastIntersects){
  ## source = total - target
  spp_area_source = mask(x=total_area, mask=trg_silh, inverse=T)
  source_sqm=RastArea(rast=spp_area_source[[1]])
  ## what is the occupied area in the target area?
  spp_area_target = crop(x=total_area, y=trg_silh, snap="near", mask=T)
  target_sqm=RastArea(rast=spp_area_target[[1]])
} else {
  spp_area_source = total_area
  source_sqm=total_sqm
  target_sqm=0
}

recordnames = c(  "species", "sp_status", "rastIntersects", "total_sqm", "polygon_shape", 
                  "spp_tot_occ", "remove_occpoint", "spp_noNA_occ",
                  "spp_tot_bgr", "remove_bgrpoint", "spp_noNA_bgr",
                  "source_sqm", "spp_scr_occ", "spp_scr_bgr", "target_sqm", "spp_trg_occ", "spp_trg_bgr")
recorddata = data.frame(this_species, sp_status, rastIntersects, total_sqm, polygon_shape, 
                        spp_tot_occ, remove_occpoint, spp_noNA_occ,
                        spp_tot_bgr, remove_bgrpoint, spp_noNA_bgr,
                        source_sqm, spp_scr_occ, spp_scr_bgr, target_sqm, spp_trg_occ, spp_trg_bgr)
#recorddata
### collect data from the species to be used later for other analysis: 
write.table(recorddata, file=file.path(dirfm_datvar_trg, "raw_data.csv"), col.names=FALSE, sep=",", append = TRUE)

if(exists(quote(report_NA_occ)) || exists(quote(report_NA_bgr))){
  lost_points=data.frame(this_species, sp_status)
  if(exists(quote(report_NA_occ))){lost_points=cbind(lost_points, report_NA_occ)} 
  if(exists(quote(report_NA_bgr))){lost_points=cbind(lost_points, report_NA_bgr)} 
  #lost_points
  write.table(lost_points, file=file.path(dirfm_datvar_trg, "raw_data_NApoints.csv"), col.names=FALSE, sep=",", append = TRUE)
}

##################
#spp_range=(total_area$bio15 >= 0) %>% na.omit()

species_range=
  ggplot() + theme_bw() +
  geom_spatvector(data = wcl, fill = "grey85", col="grey90") +
  geom_spatvector(data = vect(spp_silh), fill = NA, col="darkorchid4", lwd=0.5) +
  geom_spatvector(data = spp_occ, fill = NA, col="darkorchid1", cex=0.2, alpha=0.5) +
  geom_spatvector(data = trg_silh, fill = NA, col="cyan", lwd=0.2) 
pname=paste0("species_known_total_range")
plot_pdf_png(species_range, dirff_datvar_tsp, pname, pointsize = 4, pdf_width = 14.5)


print("data preparation finishes")

############################################################################

print("start calculations")

## dismo requires objects of the class raster, or plain matrices or data frames
occ_source  = as.data.frame(spp_occ_source, geom="XY") 
bgr_source  = as.data.frame(spp_bgr_source, geom="XY")
area_source = raster::stack(spp_area_source)

print("calculating enmeval")
enmeval_models <- ENMevaluate(occs = occ_source, bg = bgr_source, envs = area_source,
                              algorithm = algor, ## 'maxnet',   ## in deigo it should be 'maxent.jar'
                              partitions = 'block', 
                              tune.args = list(fc = c("LQH"), rm = c(2:5)))  ## 

## save the enmeval object to file
print("save enmeval to file")
filename=file.path(dirfm_datvar_tsp, paste0("enmeval_", algor, "_allmodels", ".rds")) 
saveRDS(enmeval_models, filename)

print("finished calculating ENMeval")

######################

print("model selection and evaluation")
## model selection following Kass, et al. (2022)
# 1) filtered out all models without non-zero coefficients
# posL = lapply(enmeval_models@models, function(m){anyPossitiveLambda(m@lambdas)}) %>% bind_rows() %>% t()
# candidate_results=enmeval_models@results[posL,]

candidate_results = enmeval_models@results %>% filter(ncoef > 0)

# 2) filtered out models that performed poorly (≤0 or NA) for the Continuous Boyce Index 
if(nrow(candidate_results)>1){
  c1=candidate_results %>% filter(cbi.val.avg>0)
  if(nrow(c1)>0){ candidate_results=c1 } else { candidate_results=candidate_results }
  # 3) selected those with the lowest 10 percentile omission rate
  min_or=candidate_results$or.10p.avg %>% min()
  candidate_results=candidate_results %>% filter(or.10p.avg==min_or)
  # 4) choosing the model with the highest validation area under the ROC curve (AUC)
  max_auc=candidate_results$auc.val.avg %>% max()
  candidate_results=candidate_results %>% filter(auc.val.avg==max_auc)
  # 5) choosing the model with the least number of parameters
  min_par=candidate_results$ncoef %>% min()
  candidate_results=candidate_results %>% filter(ncoef==min_par)
}
bm_metadata=candidate_results

## the best models after the filtering criteria
bm_name=bm_metadata$tune.args %>% as.character()

## get the best model information from the lists:
best_model=enmeval_models@models[[bm_name]]                 ## the selected model
bm_varImp=enmeval_models@variable.importance[[bm_name]]     ## variable importance
#bm_lambdas=best_model[["betas"]]               ### for algor=maxnet 
#bm_lambdas_used=best_model[["samplemeans"]] %>% names()  ### for algor=maxnet 
bm_lambdas=best_model@lambdas %>% lambdasDF()               ## lambdas file
bm_lambdas_used=regmatches(bm_lambdas$V1,
                           regexpr("bio[0-9][0-9]+", bm_lambdas$V1)) %>% unique()  ## these are to run the mess with them
bm_prediction=enmeval_models@predictions[[bm_name]]        ## prediction of the model

## save a bunch of stuff of the best_model to file
filename=file.path(dirfm_datvar_tsp, paste0("variable_importance.csv")) 
write.table(bm_varImp, file=filename, col.names = F, row.names = F)
filename=file.path(dirfm_datvar_tsp, paste0("lambdas_raw.csv")) 
write.table(best_model@lambdas, file=filename, col.names = F, row.names = F)
filename=file.path(dirfm_datvar_tsp, paste0("lambdas_dataframe.csv")) 
write.table(bm_lambdas, file=filename, col.names = F, row.names = F)
filename=file.path(dirfm_datvar_tsp, paste0("lambdas_used.csv")) 
write.table(bm_lambdas_used, file=filename, col.names = F, row.names = F)
filename=file.path(dirfm_datvar_tsp, paste0("suitability_prediction_source",".tif")) 
writeRaster(bm_prediction, filename=filename, overwrite=TRUE)

source_plot=
  ggplot() + theme_bw() +
  geom_spatvector(data = wcl, fill = "grey85", col="grey90") +
  geom_spatraster(data = rast(bm_prediction)) +
  scale_fill_whitebox_c(palette = "deep") +
  labs(fill = "suitability") + 
  theme(legend.position = "top", legend.title = element_text(size=9), legend.text = element_text(size=7),
        legend.key.height= unit(0.7, 'cm'), legend.key.width= unit(2, 'cm')) 

pname=paste0("suitability_prediction_source")
plot_pdf_png(source_plot, dirff_datvar_tsp, pname, pointsize = 4, pdf_width = 12.7)


fname=paste0("ModelResponse.pdf")
pdf(file=file.path(dirff_datvar_tsp, fname))
#plot(best_model)
dismo::response(best_model) 
dev.off()


print("measure how good the best model is for the source_area")

best_model_partitions=enmeval_models@results.partitions %>% filter(tune.args == bm_name)

### worst performing partitions:
auc.val.Pmin=min(best_model_partitions$auc.val)
cbi.val.Pmin=min(best_model_partitions$cbi.val)
or.10p.Pmax=max(best_model_partitions$or.10p)

print("external calculation of the ten-percentile omission rate, and few other measures")

occ_suit=raster::extract(bm_prediction, occ_source)
occ_suit=occ_suit[!is.na(occ_suit)]

## 2022-04-20: the 10 percentil omission rate as threshold 
## with Jamie's function
calc.10p.enmeval=calc.10p.trainThresh(occ_suit)

#### my previous way to calculate it
src_suit_min=min(occ_suit)
src_suit_mean=mean(occ_suit)
src_suit_med=median(occ_suit)
src_suit_max=max(occ_suit)

src_suit_q10=quantile(occ_suit, prob=0.1)
low_suit_n = occ_suit[occ_suit<src_suit_q10] %>% length()
src_or_10p=low_suit_n /spp_scr_occ

bm_metadata=data.frame(this_species, sp_status, candidate_results, 
                       auc.val.Pmin, cbi.val.Pmin, or.10p.Pmax, calc.10p.enmeval,
                       src_suit_min, src_suit_mean, src_suit_med, src_suit_max, 
                       src_suit_q10, src_or_10p)
bm_metadata_names = c("this_species", "sp_status", "fc", "rm", "tune.args", 
                      "auc.train", "cbi.train", "auc.diff.avg", "auc.diff.sd", "auc.val.avg", 
                      "auc.val.sd", "cbi.val.avg", "cbi.val.sd", "or.10p.avg", "or.10p.sd", 
                      "or.mtp.avg", "or.mtp.sd", "AICc", "delta.AICc", "w.AIC", "ncoef",
                      "auc.val.Pmin", "cbi.val.Pmin", "or.10p.Pmax", "calc.10p.enmeval",
                      "src_suit_min", "src_suit_mean", "src_suit_med", "src_suit_max", 
                      "src_suit_q10", "src_or_10p")

write.table(bm_metadata, file=file.path(dirfm_datvar_trg, "selected_model_metadata.csv"), col.names=FALSE, sep=",", append = TRUE)

#############

print("transfer model to target region")

## make predictions to target area
bm_transfer = dismo::predict(best_model, target_area) ## suitability values for target_region points
filename=file.path(dirfm_datvar_tsp, paste0("suitability_prediction_transfer",".tif")) 
writeRaster(bm_transfer, filename=filename, overwrite=TRUE)

bm_transfer_rast = rast(bm_transfer)

transfer_plot=
  ggplot() + theme_bw() +
  geom_spatvector(data = trg_div, fill = "grey85", col="grey90") +
  geom_spatraster(data = bm_transfer_rast) +
  geom_spatvector(data = spp_occ_target, fill = NA, col="darkorchid1", cex=0.2, alpha=0.5) +
    scale_fill_whitebox_c(palette = "deep") +
  labs(fill = "suitabitily") + 
  theme(legend.position = "top", legend.title = element_text(size=9), legend.text = element_text(size=7),
        legend.key.height= unit(0.7, 'cm'), legend.key.width= unit(2, 'cm')) 

pname=paste0("suitability_prediction_transfer")
plot_pdf_png(transfer_plot, dirff_datvar_tsp, pname, pointsize = 4)

#spp_occ_target
trg_suit=as.data.frame(bm_transfer_rast)[[1]]

trg_suit_min=min(trg_suit)
trg_suit_mean=mean(trg_suit)
trg_suit_med=median(trg_suit, na.rm = T)
trg_suit_max=max(trg_suit)

## decide threshold for species
#suit_thr=bm_metadata$or.10p.avg
true_thr=max(trg_suit_min, calc.10p.enmeval)
or_thr=(true_thr == calc.10p.enmeval)

#############
## calculate the mess to check that there are no weird extrapolations 
print("calculating MESS")

## subset only the lambdas that were used in the best_model
#bm_lambdas_used=biovars
spp_area_source_usedL=subset(spp_area_source, subset=bm_lambdas_used)  ## the species' shape
target_area_usedL=subset(target_area, subset=bm_lambdas_used)

spp_envvals = terra::extract(x=spp_area_source_usedL, y=spp_occ_source, ID=F)
mess_raster=dismo::mess(target_area_usedL, spp_envvals, full=T) # 2022.03.23 to get only one layer, full=TRUE)
names(mess_raster) = c(biovars, "mess")
mess_raster=rast(mess_raster)
mess_raster$mess=Inf_to_NaN(mess_raster$mess)

## save the raster to file
print("save mess to file")
filename=file.path(dirfm_datvar_tsp, paste0("MESS_full",".tif"))
writeRaster(mess_raster, filename=filename, overwrite=TRUE)
#mess_raster=rast(filename)
mess_minmax=t(minmax(mess_raster$mess))

print("finished MESS calc and save")

bm_transdata=data.frame(this_species, sp_status, trg_suit_min, trg_suit_mean, 
                        trg_suit_med, trg_suit_max, or_thr, true_thr, mess_minmax)

#############

if(sp_status=="invader"){
## evaluate model (can be done only for invaders)
print("evaluating how good the transfer is")

## evaluate the model on the expected presences and absences
## presence, exact points where species has been collected (in the target region)
## absence, randomly sampled points in the out_area of the species
occ_target=as.data.frame(spp_occ_target, geom="XY") %>% select("x", "y")  
bgr_target=dismo::randomPoints(target_area, n=1000, p=occ_target)
transfer_eval=dismo::evaluate(p=occ_target, a=bgr_target, model=best_model, x=target_area)
trg_auc=transfer_eval@auc
trg_cor=transfer_eval@cor

filename=file.path(dirfm_datvar_tsp, paste0("bm_transfer_evaluation", ".rds")) 
saveRDS(transfer_eval, filename)

print("calculate boyce index for target")

fname=paste0("predicted_expected_ratio_target_ecospatboyce.pdf")
pdf(file=file.path(dirff_datvar_tsp, fname))
transferBoyce=ecospat::ecospat.boyce(fit=bm_transfer, obs=occ_target, nclass=0,
                                    window.w="default", res=100, PEplot = T, method="spearman")
dev.off()
ext_trg_boyce=transferBoyce$cor

# #trg_suit
# trg_occ_suit=terra::extract(x=bm_transfer_rast, y=spp_occ_target, ID=F)
# fname=paste0("predicted_expected_ratio_target_ecospatboyce.pdf")
# pdf(file=file.path(dirff_datvar_tsp, fname))
# transferBoyce=ecospat::ecospat.boyce(fit=trg_suit, obs=trg_occ_suit[1,], nclass=0, 
#                                      window.w="default", res=100, PEplot = T, method="spearman")
# dev.off()
# ext_trg_boyce=transferBoyce$cor

bm_transdata=data.frame(bm_transdata, trg_auc, trg_cor, ext_trg_boyce)
}

print("finished evaluating how good the transfer is")

bm_transdata_names = c("this_species", "sp_status", "trg_suit_min", "trg_suit_mean", 
                       "trg_suit_med", "trg_suit_max", "or_thr", "true_thr", 
                       "mess_min", "mess_max", "trg_auc", "trg_cor", "ext_trg_boyce")

write.table(bm_transdata, file=file.path(dirfm_datvar_trg, "selected_model_transferdata.csv"), col.names=FALSE, sep=",", append = TRUE)


# ## external evaluation of the boyce index: 
# print("calculate boyce index at source region")
# 
# fname=paste0("predicted_expected_ratio_source_ecospatboyce.pdf")
# pdf(file=file.path(dirff_datvar_tsp, fname))
# #try(
# ecospatBoyce=ecospat::ecospat.boyce(fit=bm_prediction, obs=occ_source, nclass=0, 
#                                     window.w="default", res=100, PEplot = T, method="spearman")
# #)
# dev.off()
# 
# ext_src_boyce=ifelse(exists("ecospatBoyce"), ecospatBoyce$cor, "NA")
# 
# bm_src_boyce=data.frame(this_species, sp_status, ext_src_boyce)
# bm_src_boyce_names = c("this_species", "sp_status", "ext_src_boyce")
# 
# write.table(bm_src_boyce, file=file.path(dirfm_datvar_trg, "selected_model_extscr_boyce.csv"), col.names=FALSE, sep=",", append = TRUE)
# 

print("finished this script")
print("bye bye")

############


##### code end  




