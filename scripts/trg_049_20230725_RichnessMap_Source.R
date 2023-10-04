## 
## author:  Yazmin Zurapiti
## title:   trg_047_20230509_RichnessMap_invaders.R
## aim:     We need to calculate occupancy at country and subregion level, 
##          and source for the potential species
## version: 20230511
## call:  lauch from /flash/EconomoU/Zurapiti as there is where the R packages are installed
##        trg_043_20230306_invaders_agg_source.R [save_date] [trg] [bioset] [algor] 
## inputs: 
##   1.inputdir: 
##   1.inputfile:
## outputs: 
##   1.outputdir: 
##   1.outputfile: 
## 

# Load packages -- the order here is important because some pkg functions overwrite others.

#install.packages("rlist")
#install.packages("tidyterra")
library(sp)
library(sf)
library(terra)
library(raster)
library(rnaturalearth)
library(ggplot2)
library(tidyterra)
# library(stringi)  ## to make a matrix out of a list of uneven number of rows 
# library(stringr)  ## for leading zeros #str_pad(lay, 2, pad = "0")
library(purrr)
library(rlist)
# library(dismo)
# library(ENMeval)
library(dplyr)

## clear memory
rm(list = ls())

# declare directories 
local_analysis=F
if(local_analysis){script06path="./06_defineDirectories.R"} else 
  {script06path="./AntInvasionRisk/scripts/06_defineDirectories.R"}
source(script06path)

# testing parameters
# save_date="2023-02-09"  # chontali
save_date="2023-05-31"  # deigo
trg="Japan"
bioset="bioset05"
algor='maxent.jar' #"maxnet"
## run_date="2023-05-08"  # chontali
run_date="2023-06-23"  # deigo
sp_status="potential" #"invader" # "prediction" # 
model_clean_criteria="fc2" # fc0, fc1, postmodeling filtering criteria, (occ/cbi/auc combinations)
div_level="cluster"  ## cluster # prefecture # block # parallel  ## the administrative level for the risk of establishment and source predictions
mask="both"  ## what mask to use for the prefecture prediction


# fc1 = cbi.val.Pmin and auc 
# fc2 = occ and cbi.val.Pmin

## thresholding values
spp_scr_occ_thr = 100
cbi.val.Pmin_thr = 0  #-0.07
auc.diff.sd_thr =  0.13

## 
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
model_clean_criteria=args[7]
div_level=args[8]
mask = args[9]


sessionInfo()
print("parameters:")
print(paste("data date =", save_date, ", target_region =", trg, 
            ", bioset =", bioset, ", algorithm =", algor, ", running date =", run_date,
            ", species type =", sp_status, # ", this_species =", this_species, 
            #", apply filter =", apply_model_filter, 
            ", model cleaning criteria =", model_clean_criteria,
            ", division level =", div_level, " mask =", mask, 
            sep=" " ))

newdata=paste0("invasionRiskData_5km_", save_date)
dirbd_ants=file.path(dir_bucket_data, newdata, "model_inputs")
dirbd_varsets=file.path(dir_bucket_data, newdata, "BioclimaticSubsets")

dirbd_env=file.path(dir_bucket_data, "EnvironmentalRasters")
dirbd_env_src=file.path(dirbd_env, "WorldClim_5km")
dirbd_env_trg=file.path(dirbd_env, trg)

datavar=paste0("run-", run_date, "_data-", save_date, "_varset-", bioset) ## dir_models
dirbm_datvar_trg=file.path(dir_bucket_models, datavar, trg)
dirbm_datvar_st=file.path(dirbm_datvar_trg, sp_status)   ## species status directory (only declared here)

### A directory to save more info

## save directories:
dirfm_datvar_trg=file.path(dir_flash_models, datavar, trg)
dirfm_datvar_st=file.path(dirfm_datvar_trg, div_level, model_clean_criteria, sp_status)   
#dirfm_datvar_h0=file.path(dirfm_datvar_trg, "predictions") 
if(!dir.exists(dirfm_datvar_st)){dir.create(dirfm_datvar_st, recursive = T)}

dirff_datvar_trg=file.path(dir_flash_figures, datavar, trg)
dirff_datvar_st=file.path(dirff_datvar_trg, div_level, model_clean_criteria, sp_status) 
#dirff_datvar_h0=file.path(dirff_datvar_trg, "predictions") 
if(!dir.exists(dirff_datvar_st)){dir.create(dirff_datvar_st, recursive = T)}

## functions
mess_reclass_zero = function(rr, thr){
  min=minmax(rr)[1] - 1   ## extending a bit the limits
  max=minmax(rr)[2] + 1   ## to prevent leaving out values
  if(max<thr){mat = c(min, max, 0)} else {
    mat <- c(min, thr, 0, ## 
             thr, max, 1)
  }
  rclmat <- matrix(mat, ncol=3, byrow=TRUE)
  reclass <- classify(rr, rclmat, include.lowest=TRUE)
  return(reclass)
}

#####
pred_reclass_zero = function(rr, thr){
  mat <- c(0, thr, 0,  ## 
           thr, 1, 1)
  rclmat <- matrix(mat, ncol=3, byrow=TRUE)
  reclass <- classify(rr, rclmat, include.lowest=TRUE)
  return(reclass)
}

mess_reclass_NA = function(rr, thr){
  min=minmax(rr)[1] - 1   ## extending a bit the limits
  max=minmax(rr)[2] + 1   ## to prevent leaving out values
  if(max<thr){mat = c(min, max, NA)} else {
    mat <- c(min, thr, 0, 
             thr, max, 1)
  }
  rclmat <- matrix(mat, ncol=3, byrow=TRUE)
  reclass <- classify(rr, rclmat, include.lowest=TRUE)
  return(reclass)
}

#####
pred_reclass_NA = function(rr, thr){
  mat <- c(0, thr, NA,  
           thr, 1, 1)
  rclmat <- matrix(mat, ncol=3, byrow=TRUE)
  reclass <- classify(rr, rclmat, include.lowest=TRUE)
  return(reclass)
}


### this functions aggregates the prefectures into a larger (sub)region
### and calculates the mean coordinates that represent the subregions  
### to arrange them by latitude
prepareSubregion = function(trg_div, div_level) {
  use_div = trg_div %>% select(all_of(div_level)) #%>% summarise()
  use_div = use_div %>% group_by(across(all_of(div_level))) %>% summarise()
  names(use_div) = c("subregion", "geometry")  ## rename columns
  subreg_names = use_div$subregion %>% unique()  ## get the subregion names
  
  coords_reg=list()
  for(j in 1:length(subreg_names)){
    print(j)
    geo=use_div[use_div$subregion==subreg_names[j],]$geometry 
    geo_flat=list.flatten(geo)
    geo_flat=lapply(geo_flat, as.data.frame)
    coords=bind_rows(geo_flat)
    colnames(coords)=c("longitude", "latitude")
    coords_reg[[j]]=colMeans(coords)%>% as.data.frame() %>% t()
  }
  names(coords_reg) = subreg_names
  coords_reg = coords_reg %>% purrr::map_dfr( ~as_tibble(.), .id="subregion" ) 
  subreg_order = use_div %>% full_join(coords_reg) %>% arrange(latitude) 
  subreg_order = subreg_order %>% 
    mutate(subregion=factor(subreg_order$subregion, levels = subreg_order$subregion))
  
  return(subreg_order)
}

### this makes a dataframe pairing the subregion with the species that could 
### inhabit given the occ_r passed by parameter
makeDF_SubregSpecies =  function(occ_r, use_div, div_level){
  names(use_div) = c("subregion", "geometry")  ## rename columns
  subreg_names = use_div$subregion %>% unique()  ## get the subregion names
  
  subreg_spp_list=list()  # the list that will return
  #j=2
  for(j in 1:length(subreg_names)){
    subreg_name = subreg_names[[j]]
    print(paste0(j, ": ", subreg_name))
    
    subregion = use_div %>% filter(subregion==subreg_name) ## regional level
    subregion_pred = crop(x=occ_r, y=subregion)
    subregion_minmax=minmax(subregion_pred) ## without zeros
    
    ## take the average suitability of the region
    ## add zeros to the prediction, where NaN of the region  ## this is new 20230725
# subregion_pred = subst(subregion_pred, NaN, 0)   ## this doesn't work!
# print(plot(subregion_pred))
    subregion_pred_df=as.data.frame(subregion_pred, na.rm=NA)
    sum = colSums(subregion_pred_df, na.rm=T)
    mean = sum/nrow(subregion_pred_df)
 
    RegionalPotential=rbind(subregion_minmax, mean) %>% t() %>% as.data.frame()
    RegionalPotential$species = rownames(RegionalPotential)  
    subreg_spp_list[[j]]=RegionalPotential
  }
  names(subreg_spp_list)=subreg_names
  
  ## make the list a tibble
  reg_spp=subreg_spp_list %>% purrr::map_dfr( ~as_tibble(.), .id="subregion" ) %>% filter(max>0)
  return(reg_spp)
}

### arrange the dataframe and count the number of subregions per species and 
### and species per subregion to make a dataframe ready to plot
countsForPlot_SubregSpecies = function(reg_spp){
  
  subreg_names = unique(reg_spp$subregion)
  ## count in how many subregions the species occur
  spp_occ = reg_spp %>% group_by(species) %>% summarise(spp_nReg=n()) %>% arrange(-spp_nReg) 
  spp_occ$species_nReg = paste0(spp_occ$species, " (", spp_occ$spp_nReg, ")")
  
  ## count number of species per subregion
  reg_occ = reg_spp %>% group_by(subregion) %>% summarise(srg_nSpp=n()) %>% 
    mutate(subregion = factor(x=subregion, levels=subreg_names)) %>% arrange(subregion)
  reg_occ$subreg_nSpp = paste0(reg_occ$subregion, " (", reg_occ$srg_nSpp, ")")
  
  reg_spp_occ = reg_spp %>% full_join(spp_occ, by="species") %>% 
    full_join(reg_occ, by="subregion") %>%  ## join info
    mutate(subreg_nSpp = factor(x=subreg_nSpp, levels=unique(reg_occ$subreg_nSpp))) %>% 
    mutate(species_nReg = factor(x=species_nReg, levels=unique(spp_occ$species_nReg)))
  return(reg_spp_occ)
}

## species - mean occurrence tile plot 
sppOccTilePlot = function(spp_occ_mtx, level){
  gSS=ggplot(data=spp_occ_mtx, aes(y = subreg_nSpp, x = species_nReg)) +
    theme() + ylab(level) + xlab("species") +
    theme(legend.position = "top", legend.text = element_text(size=6),
          legend.key.height= unit(0.7, 'cm'), legend.key.width= unit(2, 'cm'),
          axis.text = element_text(size=6), axis.text.x = element_text(angle = 60, hjust=0.95)) +
    geom_tile(aes(fill=mean)) + scale_fill_whitebox_c("mean suitability", palette = "deep")
  return(gSS)
}

## to plot species source 
species_source_plot=function(raster_stack, labname){
  ggplot() + theme_bw() + scale_fill_whitebox_c(palette = "deep") + labs(fill = labname) + 
    theme(legend.position = "top", legend.title = element_text(size=9), legend.text = element_text(size=7),
          legend.key.height= unit(0.7, 'cm'), legend.key.width= unit(2, 'cm')) +
    geom_spatvector(data = wcl, fill = "grey85", col="grey90") +
    geom_spatraster(data = raster_stack)
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

##### main
print("data preparation starts")

## one raster from the environmental variables to get the resolution to rasterize the spp source shapes
worldTemplate=rast(file.path(dir_bucket_data, "EnvironmentalRasters","WorldClim_5km","bio12.tif"))
wcl=vect(ne_coastline()) ## plotting purposes

trg_div = readRDS(file = file.path(dirbd_env_trg, "Subregions.rds")) ## sf class
trg_silh = trg_div %>% vect()  ## terra class, for plotting

pHC_file=file.path(dirbm_datvar_trg, "PrefectureClusterColor.csv")
if(file.exists(pHC_file)){ 
  pref_cluster_df = read.csv(pHC_file, row.names = 1, header = T)
  trg_div = trg_div %>% full_join(pref_cluster_df, by = join_by(prefecture)) 
}

## species data

recordnames = c(  "species", "sp_status", "rastIntersects", "total_sqm", "polygon_shape", 
                  "spp_tot_occ", "remove_occpoint", "spp_noNA_occ",
                  "spp_tot_bgr", "remove_bgrpoint", "spp_noNA_bgr",
                  "source_sqm", "spp_scr_occ", "spp_scr_bgr", "target_sqm", "spp_trg_occ", "spp_trg_bgr")
recorddata=data.frame( read.table(file=file.path(dirbm_datvar_trg, "raw_data.csv"), sep=",") ) %>% select(-V1)
colnames(recorddata) = recordnames

bm_metadata_names = c("species", "sp_status", "fc", "rm", "tune.args", 
                      "auc.train", "cbi.train", "auc.diff.avg", "auc.diff.sd", "auc.val.avg", 
                      "auc.val.sd", "cbi.val.avg", "cbi.val.sd", "or.10p.avg", "or.10p.sd", 
                      "or.mtp.avg", "or.mtp.sd", "AICc", "delta.AICc", "w.AIC", "ncoef",
                      "auc.val.Pmin", "cbi.val.Pmin", "or.10p.Pmax", "calc.10p.enmeval",
                      "src_suit_min", "src_suit_mean", "src_suit_med", "src_suit_max", 
                      "src_suit_q10", "src_or_10p")#, "ext_src_boyce")
bm_metadata=data.frame( read.table(file=file.path(dirbm_datvar_trg, "selected_model_metadata.csv"), sep=",") ) %>% select(-V1)
colnames(bm_metadata) = bm_metadata_names

# bm_src_boyce_names = c("species", "sp_status", "ext_src_boyce")
# bm_src_boyce=data.frame( read.table(file=file.path(dirbm_datvar_trg, "selected_model_extscr_boyce.csv"), sep=",") ) %>% select(-V1)
# colnames(bm_src_boyce) = bm_src_boyce_names

bm_transdata_names = c("species", "sp_status", "trg_suit_min", "trg_suit_mean", 
                       "trg_suit_med", "trg_suit_max", "or_thr", "true_thr", 
                       "mess_min", "mess_max", "trg_auc", "trg_cor", "ext_trg_boyce")
bm_transdata=data.frame( read.table(file=file.path(dirbm_datvar_trg, "selected_model_transferdata.csv"), sep=",", fill=T) ) %>% select(-V1)
colnames(bm_transdata) = bm_transdata_names

#info_all = bm_metadata %>% full_join(bm_src_boyce, by=c("species", "sp_status"))
info_all = recorddata %>% full_join(bm_metadata, by=c("species", "sp_status")) %>% 
  full_join(bm_transdata, by=c("species", "sp_status"))

## apply (or not) the filter by spp_status
spp_list = info_all %>% filter(sp_status == UQ(sp_status)) %>% arrange(species) 
prefilter_n = nrow(spp_list)

#if(apply_model_filter){
  # apply or not the threshold model values found by the trg_045 code
print(paste0("apply postmodeling filtering (transferability) criteria: ", model_clean_criteria))
if(model_clean_criteria == "fc0"){remove = spp_list %>% filter(!(species %in% spp_list$species) ) %>% select(species)}
if(model_clean_criteria == "fc1"){remove = spp_list %>% filter(auc.diff.sd >= auc.diff.sd_thr | cbi.val.Pmin <= cbi.val.Pmin_thr) %>% select(species)} ## species to be removed
if(model_clean_criteria == "fc2"){remove = spp_list %>% filter(spp_scr_occ < spp_scr_occ_thr | cbi.val.Pmin <= cbi.val.Pmin_thr) %>% select(species)} ## species to be removed
spp_list = spp_list %>% filter( !(species %in% remove$species) ) ## species to be included
#}
postfilter_n = nrow(spp_list)
sp_names = spp_list$species 
#if(local_analysis){sp_names=sp_names[c(1,6,7,10,15)]}
## 32 invader species, 26 retained after "model_filter"; 6 filtered
## 
######

## this part reads and thresholds the species 
## suitability prediction in Japan and mess rasters
print("read source, transfer and mess rasters and prepare them for analysis")
source_glb = list()
transf_full = list()
transf_mess = list()
transf_OR10 = list()
transf_both = list()
#i=86
for(i in 1:length(sp_names)){
  this_species=sp_names[i]
  print(paste0(i, ": ", this_species))
  
  spp_data = spp_list %>% filter(species == this_species) #%>% select(true_thr)
  source_thr = spp_data$calc.10p.enmeval  ## source OR10 threshold
  transf_thr = spp_data$true_thr  ## transfer OR10 (or min transfer value, whichever smaller) threshold
  
  dirbm_datvar_tsp = file.path(dirbm_datvar_st, this_species)   ## species directory (only declared here)
   # source_file = file.path(dirbm_datvar_tsp, paste0("source_prediction.tif"))   # chontali file
   # pred_file=file.path(dirbm_datvar_tsp, paste0("prediction_", "full",".tif"))  # chontali file
  source_file = file.path(dirbm_datvar_tsp, paste0("suitability_prediction_source.tif"))
  transf_file = file.path(dirbm_datvar_tsp, paste0("suitability_prediction_transfer.tif"))
    mess_file = file.path(dirbm_datvar_tsp, paste0("MESS_full",".tif"))
  
  if(file.exists(source_file) & file.exists(transf_file) & file.exists(mess_file)){
    print(paste0("files exist"))
    ## file exist, so we read the files
    source_f = rast(source_file)
    transf_f = rast(transf_file)
      mess_f = rast(mess_file)$mess
    ## create the threshold predictions
    source_tNA = pred_reclass_NA(source_f, source_thr) # NA, 1 prediction
    transf_tNA = pred_reclass_NA(transf_f, transf_thr)
      mess_tNA = mess_reclass_NA(mess_f, 0) 
    source_tze = pred_reclass_zero(source_f, source_thr)
    transf_tze = pred_reclass_zero(transf_f, transf_thr)
      mess_tze = mess_reclass_zero(mess_f, 0) 
    
    ## correct predictions: threshold the predictions and extend
     source_glb[[i]] = extend((source_f * source_tNA), ext(worldTemplate))
    transf_full[[i]] = transf_f
    transf_mess[[i]] = (transf_f * mess_tze)
    transf_OR10[[i]] = (transf_f * transf_tze)
    transf_both[[i]] = (transf_f * transf_tze * mess_tze) 
  } else {print(paste0("some files don't exist"))}
}
names(source_glb)=sp_names
names(transf_full)=sp_names
names(transf_mess)=sp_names
names(transf_OR10)=sp_names
names(transf_both)=sp_names

 source_glb_C = compact(source_glb)
transf_full_C = compact(transf_full)
transf_mess_C = compact(transf_mess)
transf_OR10_C = compact(transf_OR10)
transf_both_C = compact(transf_both)

### make a list of all post-modelling cleaning scenarios
masks_n=c("full", "mess", "OR10", "both")
masks_d=list(transf_full_C, transf_mess_C, transf_OR10_C, transf_both_C)
names(masks_d) = masks_n

##########   the tile plot at prefecture level
subdiv = prepareSubregion(trg_div, div_level) # as before
use_div = subdiv %>% select(subregion, geometry)

#for(k in 1:length(filters_n)){
print(paste0("prediction mask type: ", mask))
## some species might not have prediction
transf_rast = rast(Filter(Negate(is.null), masks_d[[mask]])) ## remove null predictions

## country wide richness prediction
rich_rast = transf_rast %>% sum(na.rm=T)     ## richness prediction
pname=paste0("richness","_MF-", model_clean_criteria,"_mask-",mask)
writeRaster(rich_rast, file.path(dirfm_datvar_st, paste0(pname, ".tif")), overwrite=T)

min=minmax(rich_rast)[1] %>% floor()
max=minmax(rich_rast)[2] %>% ceiling()
rich_rast_NA = pred_reclass_NA(rich_rast, 0)
richness_plot=
  ggplot() + theme_bw() +
  geom_spatvector(data = trg_silh, fill = "grey85", col="grey90") +
  geom_spatraster(data = rich_rast_NA) +
  scale_fill_whitebox_c(palette = "deep") +
  labs(fill = "number of species") + 
  theme(legend.position = "top", legend.title = element_text(size=9), legend.text = element_text(size=7),
        legend.key.height= unit(0.7, 'cm'), legend.key.width= unit(2, 'cm')) 
richness_plot
plot_pdf_png(richness_plot, dirff_datvar_st, pname, pointsize = 4)

reg_spp = makeDF_SubregSpecies(transf_rast, use_div, div_level)
## write this table to file
fname=paste0("Species-",div_level,"_MF-",model_clean_criteria,"_mask-",mask,".csv")
write.csv(reg_spp, file=file.path(dirfm_datvar_st, fname))

reg_spp_occ = countsForPlot_SubregSpecies(reg_spp)
gSS=sppOccTilePlot(reg_spp_occ, div_level)
gSS
pname=paste0("Species-", div_level,"_MF-", model_clean_criteria, "_mask-",mask)
plot_pdf_png(gSS, dirff_datvar_st, pname, pointsize = 4)

###################
## source for the same admin level as tile plot above


### HERE: almost done the above procedure, 
## still making sure we are getting zeros and NAs in the correct places so I can 
## continue with the source and also be sure we are getting the correct prediction


#   ## same plot by block and parallel 
#   source_level = "block"
#   #source_level = "parallel"
#   srclvl_div = trg_div %>% mutate(subregion = get(eval(div_level))) %>% as_tibble()
#   
#   newgroup=c("species", eval(source_level))
#   srclvl_meanOcc = reg_spp_occ %>% full_join(srclvl_div) %>% group_by(across(all_of(newgroup))) %>% 
#     summarise(mean = mean(mean)) %>% mutate(subregion = get(eval(source_level)))
#   newreg_spp_occ = countsForPlot_SubregSpecies(srclvl_meanOcc) %>% na.omit()
#   fname=paste0("Species-",source_level,"_MF-",model_clean_criteria,"_mask-",filter_name,".csv")
#   write.csv(newreg_spp_occ, file=file.path(dirfm_datvar_st, fname))
#   #newreg_spp_occ = read.csv(file=file.path(dirbm_datvar_st, fname))
#   
#   gSS=sppOccTilePlot(newreg_spp_occ, source_level)
#   pname=paste0("Species-", source_level,"_MF-",model_clean_criteria,"_mask-",filter_name)
#   plot_pdf_png(gSS, dirff_datvar_st, pname, pointsize = 4)
# 
# }

## filter out the species that had zero suitability for the whole country
## this is the last reg_spp_occ calculated; this is with both masks for the current code
sp_names_suit = unique(reg_spp$species)
establish_n = length(sp_names_suit) 

print(paste0("potential total: ", prefilter_n))
print(paste0("potential post cleaning: ", postfilter_n))
print(paste0("potential establishment: ", establish_n))

source_glb_list = source_glb[sp_names_suit]
length(source_glb_list)

## country-level source prediction  (no weighting)
source_glb_sum=rast(source_glb_list) %>% sum(na.rm=T)  ## stack and sum the list, 
source_glb_sum=mask(source_glb_sum, mask=worldTemplate) ## remove water bodies

pname=paste0("Source_", trg,"_MF-", model_clean_criteria,"_mask-", mask)
writeRaster(source_glb_sum, file.path(dirfm_datvar_st, paste0(pname, ".tif")), overwrite=T)
source_plot=species_source_plot(source_glb_sum, "sum of species suitability")
source_plot
plot_pdf_png(source_plot, dirff_datvar_st, pname, pointsize = 4, pdf_width = 12.7)


### subregion-level source prediction


# ### first the species that are "everywhere"
# ### "ubiquitous" species
# subreg_name = "ubiquituos"
# 
# # widely spread species (prefecture)
# Nreg = trg_div$prefecture %>% unique() %>% length()
# Nreg90 = ceiling( Nreg*0.9 )
# ubireg = reg_spp_occ %>% filter(spp_nReg >= Nreg90)
# srg_species=ubireg$species %>% unique()
# srg_species %>% length()
# 
# ## widely spread species (block)
# Nblock = trg_div$block %>% unique() %>% length()
# Nblock90 = ceiling( Nblock*0.9 )
# ubiblock = newreg_spp_occ %>% filter(spp_nReg >= Nblock90) 
# sbl_species=ubiblock$species %>% unique()
# sbl_species %>% length()
# 
# ubi_spp=intersect(srg_species, sbl_species)
# ubi_spp %>% length()
# fname=paste0("UbiquituousSpecies-",source_level,"_MF-",model_clean_criteria,"_mask-",filter_name,".csv")
# write.csv(ubi_spp, file=file.path(dirfm_datvar_st, fname))
# 
# srg_sp_list=species_shape_list2[ubi_spp]
# species_shape_stack=rast(srg_sp_list) 
# species_shape_stack=sum(species_shape_stack, na.rm=T)  ## stack and sum the list, 
# species_shape_stack=mask(species_shape_stack, mask=worldTemplate) ## remove water bodies
# 
# pname=paste0("Source_", source_level, "-", subreg_name,"_MF-", model_clean_criteria,"_mask-",filter_name)
# writeRaster(species_shape_stack, file.path(dirfm_datvar_st, paste0(pname, ".tif")))
# subreg_source_plot=species_source_plot(species_shape_stack, "sum of species suitability")
# plot_pdf_png(subreg_source_plot, dirff_datvar_st, pname, pointsize = 4, pdf_width = 12.7)
# 
# ### then the "rest"
# subreg_name = "rest"

# rest = reg_spp_occ %>% filter(!(species %in% ubi_spp)) %>% distinct(species) 
# subreg_spp = srclvl_meanOcc %>% filter(species %in% rest$species)
# subreg_names = subreg_spp$subregion %>% unique() 



#srg_sp_list=species_shape_list2[rest$species]
# species_shape_stack=rast(srg_sp_list) 
# species_shape_stack=sum(species_shape_stack, na.rm=T)  ## stack and sum the list, 
# species_shape_stack=mask(species_shape_stack, mask=worldTemplate) ## remove water bodies

# pname=paste0("Source_", source_level, "-", subreg_name,"_MF-", model_clean_criteria,"_mask-",filter_name)
# writeRaster(species_shape_stack, file.path(dirfm_datvar_st, paste0(pname, ".tif")))
# subreg_source_plot=species_source_plot(species_shape_stack, "sum of species suitability")
# plot_pdf_png(subreg_source_plot, dirff_datvar_st, pname, pointsize = 4, pdf_width = 12.7)
# 
# j=2
# 
# reg_spp  ## this have the region, species and weights
# source_glb_list ## this is the source prediction for all exotics species

## get the unique subregions (admin level specified by the div_level parameter)
subreg_names=unique(reg_spp$subregion)  

#j=4
for(j in 1:length(subreg_names)){
  subreg_name = subreg_names[[j]]
  print(paste0(j, ": ", subreg_name))
  # select species and weights for the species in this subregion
  srg_species = reg_spp %>% filter(subregion == subreg_name) %>% na.omit()

if(nrow(srg_species) > 0){
 #   srg_sppNames = srg_species$species 
#    i=4
  weightedSource = list()
  for(i in 1:nrow(srg_species)){
  spName = srg_species$species[i]
  spWeight = srg_species$mean[i]
  print(paste0(spName, " weight: ", spWeight))
  weightedSource[[i]]=(source_glb_list[[spName]] * spWeight) 
  }
names(weightedSource) = srg_species$species
}

source_reg_stack=rast(weightedSource) 
source_reg_sum=sum(source_reg_stack, na.rm=T)  ## stack and sum the list, 
source_reg_sum=mask(source_reg_sum, mask=worldTemplate) ## remove water bodies

pname=paste0("Source_", div_level, "-", subreg_name,"_MF-", model_clean_criteria,"_mask-",mask)
writeRaster(source_reg_sum, file.path(dirfm_datvar_st, paste0(pname, ".tif")), overwrite=T)
subreg_source_plot=species_source_plot(source_reg_sum, "sum of species weighted suitability")
subreg_source_plot
plot_pdf_png(subreg_source_plot, dirff_datvar_st, pname, pointsize = 4, pdf_width = 12.7)
  
}




#### 
print("finished aggregation")

############
print("owaru")
print("bye bye")


##### code end  




