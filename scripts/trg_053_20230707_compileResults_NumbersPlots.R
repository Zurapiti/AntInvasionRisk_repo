## 
## author:  Yazmin Zurapiti
## title:   trg_043_20230306_invaders_agg_source.R
## aim:     We need to calculate occupancy at country and subregion level, 
##          and source for the potential species
## version: 3.0
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
library(tidyr)
library(sp)
library(sf)
library(terra)
library(raster)
library(rnaturalearth)
library(ggplot2)
library(ggpubr)
library(tidyterra)
library(ggdendro)

#library(ggtree)
library(aplot)
#library("dendsort")
library(dendextend)

# library(stringi)  ## to make a matrix out of a list of uneven number of rows 
# library(stringr)  ## for leading zeros #str_pad(lay, 2, pad = "0")
# library(dismo)
# library(ENMeval)
library(grid)
library(purrr)
library(rlist)
library(dplyr)

## clear memory
rm(list = ls())

# testing parameters
#save_date="2023-02-09"
save_date="2023-05-31"
trg="Japan"
bioset="bioset05"
algor='maxent.jar' #"maxnet"
run_date="2023-06-23"
#run_date="2023-05-08"
#sp_status="potential" # "invader" # 
#this_species="Myrmica_rubra" # "Tetramorium_sericeiventre" # 
#sp_status="invader" # "potential" # 
#this_species="Pheidole_megacephala" # "Paratrechina_longicornis" # "Tapinoma_melanocephalum" # 

# ## get parameters
# args <- commandArgs(trailingOnly = TRUE)   ## allows to pass parameters to the script
# # test if there is at least one argument: if not, return an error
# if (length(args)<1) {
#   stop("Not enough parameters", call.=FALSE)
# }
# save_date=args[1]
# trg=args[2]
# bioset=args[3]
# algor=args[4]
# run_date=args[5]
# #sp_status=args[6]
# #this_species=args[7]

sessionInfo()

print("parameters:")
print(paste("data date =", save_date, ", target_region =", trg, 
            ", bioset =", bioset, ", algorithm =", algor, ", running date =", run_date,
            #", species type =", sp_status, ", this_species =", this_species, 
            sep=" " ))


# declare directories 
print("setting up directories")

local_analysis=T
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
dirbm_datvar_trg_cl=file.path(dirbm_datvar_trg, "clustering")

#dirbm_datvar_tsp=file.path(dirbm_datvar_trg, sp_status, this_species)   ## species directory (only declared here)

## save directories:
dirfm_datvar_trg=file.path(dir_flash_models, datavar, trg)
#dirfm_datvar_tsp=file.path(dirfm_datvar_trg, sp_status, this_species)   ## species directory (only declared here)
if(!dir.exists(dirfm_datvar_trg)){dir.create(dirfm_datvar_trg, recursive = T)}

dirff_datvar_trg=file.path(dir_flash_figures, datavar, trg)
dirff_datvar_trg_fM = file.path(dirff_datvar_trg, "figsManuscript")
#dirff_datvar_tsp=file.path(dirff_datvar_trg, sp_status, this_species)   ## species directory (only declared here)
if(!dir.exists(dirff_datvar_trg_fM)){dir.create(dirff_datvar_trg_fM, recursive = T)}

## functions

lmfitmydata = function(minimal_data){
  colnames(minimal_data) = c("x", "y")
  minimal_data %>% na.omit()
  lm=lm(y ~ x, data=minimal_data)
  print(summary(lm))
  minimal_data = cbind(minimal_data, predict(lm, interval = 'confidence') )
  return(minimal_data)
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

#species_source_plot=function(raster_stack, labname){
#}

## this function helps to save pdf and png plots (from an object) at the same time
plot_pdf_png = function(the_plot, plot_dir, plot_name, pointsize=12, pdf_width=7, pdf_height = 7, bg="transparent"){
  fname=paste0(plot_name,".pdf")
  pdf(file=file.path(plot_dir, fname), pointsize=pointsize, width = pdf_width, height = pdf_height, bg="transparent")
  print(the_plot)
  dev.off()
  
  fname=paste0(plot_name,".png")
  png(file=file.path(plot_dir, fname), pointsize=pointsize, width = round(pdf_width*68.6), height = round(pdf_height*68.6), bg="transparent")
  print(the_plot)
  dev.off()
}

##### main
print("read data")
## one raster from the environmental variables to get the resolution to rasterize the spp source shapes
worldTemplate=rast(file.path(dir_bucket_data, "EnvironmentalRasters","WorldClim_5km","bio12.tif"))
wcl=vect(ne_coastline()) ## plotting purposes

trg_div = readRDS(file = file.path(dirbd_env_trg, "Subregions.rds")) ## sf class
trg_silh = trg_div %>% vect()  ## terra class, for plotting
trg_jp = ne_countries(scale = 50, country = "japan") %>% vect()

##### P1: report all numbers to be written in the manuscript 
##        and double check that are coherent across files
sp_list <- readRDS(file.path(dirbd_ants, "gabidb_list.rds"))
n_gabispp = length(sp_list)
sp_list <- readRDS(file.path(dirbd_ants, "exotics_gabidb_list_thin.rds"))
n_gabiexo = length(sp_list)

names_allexo = read.csv(file=file.path(dirbd_ants, "AllExoticNames.csv"), header = T)
n_allexo = nrow(names_allexo)

names_nat = read.csv(file=file.path(dirbd_env_trg, "native.csv"), header = F)
colnames(names_nat) = "species"
n_nat = nrow(names_nat)
names_inv = read.csv(file=file.path(dirbd_env_trg, "invader.csv"), header = F)
colnames(names_inv) = "species"
n_inv = nrow(names_inv)
names_pot = read.csv(file=file.path(dirbd_env_trg, "potential.csv"), header = F)
colnames(names_pot) = "species"
n_pot = nrow(names_pot)
n_exo = n_nat + n_inv + n_pot
n_ana = n_inv + n_pot

##### 

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

bm_transdata_names = c("species", "sp_status", "trg_suit_min", "trg_suit_mean", 
                       "trg_suit_med", "trg_suit_max", "or_thr", "true_thr", 
                       "mess_min", "mess_max", "trg_auc", "trg_cor", "ext_trg_boyce")
bm_transdata=data.frame( read.table(file=file.path(dirbm_datvar_trg, "selected_model_transferdata.csv"), sep=",", fill=T) ) %>% select(-V1)
colnames(bm_transdata) = bm_transdata_names

#####

info_all = recorddata %>% full_join(bm_metadata, by=c("species", "sp_status")) %>% 
  full_join(bm_transdata, by=c("species", "sp_status"))
n_aall = nrow(info_all)

info_all_supmat = info_all %>% select(species, sp_status, polygon_shape, spp_tot_occ, spp_scr_occ, 
                                      tune.args, auc.train, cbi.train, 
                                      cbi.val.Pmin, auc.diff.sd, trg_auc, trg_cor, ext_trg_boyce)

write.table(info_all_supmat, file.path(dirfm_datvar_trg,"SupplementalTable.csv"), row.names = F)

###




info_inv = info_all %>% filter(species %in% names_inv$species)
info_pot = info_all %>% filter(species %in% names_pot$species)
n_ainv = nrow(info_inv)
n_apot = nrow(info_pot)


######  check that all numbers make sense:
print(paste0("total number of ant species in gabi: ", n_gabispp))
print(paste0("total number of all identified exotic ant species in gabi: ", n_allexo))
print(paste0("total number of cleaned exotic ant (20+ records): ", n_gabiexo))
print(paste0("total number of alien ant species: ", n_exo))
print(paste0("number of analyzed native species: ", n_nat))
print(paste0("total number of analyzed species: ", n_ana))
print(paste0("number of analyzed invader species: ", n_inv))
print(paste0("number of analyzed exotic species: ", n_pot))

print(paste0("is the total number of global exotic species congruent? ", n_gabiexo == n_exo))
print(paste0("is the total number of analyzed species congruent? ", n_ana  == n_aall))
print(paste0("pre and post analysis invader species congruent? ", n_inv  == n_ainv))
print(paste0("pre and post analysis potential species congruent? ", n_pot  == n_apot))

print(paste0("total number of records for all species: ", sum(info_all$spp_tot_occ)))
print(paste0("total number of records for invader species: ", sum(info_inv$spp_tot_occ)))
print(paste0("number of source records: ", sum(info_inv$spp_scr_occ)))
print(paste0("number of target records: ", sum(info_inv$spp_trg_occ)))
print(paste0("total number of records for exotic species: ", sum(info_pot$spp_tot_occ)))

######
info_all %>% distinct(polygon_shape)
info_all %>% select(auc.train) %>% range()
info_all %>% select(auc.val.avg) %>% range()
info_all %>% select(cbi.train) %>% range()
info_all %>% select(cbi.val.avg) %>% range()


## now our evaluation on how well these transfers happened:
### Invaders data:

oppo = c("Tetramorium_smithi", "Cardiocondyla_itsukii", "Tetramorium_kraepelini", "Cardiocondyla_kagutsuchi")
over = c("Ooceraea_biroi", "Pheidole_parva", "Hypoponera_ragusai" ,"Nylanderia_amia", "Tetraponera_allaborans")
under = c("Tapinoma_melanocephalum", "Tetramorium_bicarinatum")
n_oppo = length(oppo)
n_over = length(over)
n_good = n_inv - n_oppo - n_over

print(paste0("results: invaders opposite models: ", n_oppo ))
print(paste0("results: invaders overpredicted models: ", n_over))
print(paste0("results: invaders expected models: ", n_good))

paste0("emph{", oppo[order(oppo)], "}")
paste0("emph{", over[order(over)], "}")

###
info_inv = info_inv %>%  
  mutate(transfer=ifelse(species %in% over, "extreme", 
                         ifelse(species %in% oppo, "opposite", "expected"))) 

## artificially change the NA values for the cbi for 0 to be able to visualize them
info_inv = info_inv %>%  
  mutate(cbi_edit = ifelse(is.na(ext_trg_boyce), 0, ext_trg_boyce)) %>% 
  mutate(transfer_edit = ifelse(species %in% under, "underpredicted",
                                ifelse(is.na(ext_trg_boyce), paste0(transfer,"_cbiNA"), transfer)))


##### Figure: transfer examples for invader species
examples = c("Anoplolepis_gracilipes", "Tetramorium_smithi", "Hypoponera_ragusai")
sp_status = "invader"

transfer_plot=list()
#ss = 3
for(ss in 1:3){
  this_species = examples[ss]
  print(this_species)
  dirbm_datvar_tsp=file.path(dirbm_datvar_trg, sp_status, this_species)   ## species directory (only declared here)

  spp_data = sp_list[[this_species]] %>% as.data.frame() # this particular species
  spp_coords_df = spp_data %>% select("lon_opt", "lat_opt") %>% rename(lon=lon_opt, lat=lat_opt) # data.frame
  spp_occ  = spp_coords_df %>% vect(geom=c("lon", "lat"))  # terra object
  crs(spp_occ) = "+proj=longlat +datum=WGS84 +no_defs"
  spp_occ_target  = crop(x=spp_occ, y=trg_silh)
  filename=file.path(dirbm_datvar_tsp, paste0("suitability_prediction_transfer",".tif")) 
  bm_transfer_rast = rast(filename)

  transfer2 = crop(bm_transfer_rast, c(123, 147, 24, 46))  
  
transfer_plot[[ss]]=
  ggplot() + theme_bw() +
  #  geom_spatvector(data = trg_silh, fill = "grey85", col="grey90") +
    geom_spatraster(data = transfer2) +
    geom_spatvector(data = spp_occ_target, fill = NA, col="darkorchid1", pch=8, cex=3, alpha=0.7) +
    scale_fill_whitebox_c(palette = "deep", expand=c(0,1)) +
    labs(fill = "suitability") + 
    theme(legend.position = "top", 
          legend.title = element_text(size=22, vjust = 0.75), 
          legend.text = element_text(size=22),
          axis.text = element_text(size=18),
          legend.key.height= unit(1, 'cm'), legend.key.width= unit(5, 'cm')) 
  
}

leg=get_legend(transfer_plot[[3]])
ptrex = ggarrange(transfer_plot[[1]], transfer_plot[[2]], transfer_plot[[3]], 
          nrow=1, ncol=3, legend.grob = leg, 
          labels=c("a: expected", "b: opposite", "c: extreme"),
          font.label = list(size = 20))
plot_pdf_png(ptrex, dirff_datvar_trg_fM, "F2_transferExamples_GoodBadUgly", pointsize = 4, pdf_width=15)

##### Figure: selection of training metrics
## transfer metrics: cor vs AUC
p3maintext=
  ggplot() + theme_bw() + labs(x="(transfer) cor", y="(transfer) AUC", col="transfer") + 
  theme(legend.position = "top", text = element_text(size=22), 
        legend.text = element_text(size=22),
        axis.text = element_text(size=18)) + 
  geom_point(data=info_inv, aes(trg_cor, trg_auc, col=transfer), cex=2.5)  +
  geom_hline(yintercept = 0.5, col="grey60", lty=2, lwd=1) +  # theoretical values
  geom_vline(xintercept = 0.0, col="grey60", lty=2, lwd=1)# +  # theoretical values
 # geom_hline(yintercept = 0.71, col="grey20", lty=3, lwd=1) + 
  #geom_vline(xintercept = 0.08, col="grey20", lty=3, lwd=1) 
p3maintext

## training metrics: transfer AUC vs training points
p_occu=
  ggplot() + theme_bw() + labs(y="(transfer) AUC", x="(training) occurrence points", col="transfer") + 
  theme(legend.position = "bottom", text = element_text(size=16), legend.text = element_text(size=14)) +
  geom_point(data=info_inv, aes(y=trg_auc, x=spp_scr_occ, col=transfer), cex=2.5) + 
  geom_smooth(data=info_inv, aes(y=trg_auc, x=spp_scr_occ), span=0.9, lwd=0.5, fill="blue", alpha=0.07)  +
  geom_hline(yintercept = 0.5, col="grey60", lty=2, lwd=1) +  # theoretical values
  geom_hline(yintercept = 0.71, col="grey60", lty=3, lwd=1) + 
  #geom_vline(xintercept = 20, col="grey20", lty=4, lwd=1) +
  geom_vline(xintercept = 100, col="grey20", lty=3, lwd=1) 
p_occu

## training metrics: transfer AUC vs training min CBI
tt_p1 = 
  ggplot() + theme_bw() + labs(x=paste0("(training) ", "cbi.val.Pmin"), y="(transfer) AUC", col="transfer") + 
  theme(legend.position = "bottom", text = element_text(size=16), legend.text = element_text(size=14)) +
  geom_point(data=info_inv, aes(x=cbi.val.Pmin, y=trg_auc, col=transfer), cex=2.5)  +
  geom_smooth(data=info_inv, aes(x=cbi.val.Pmin, y=trg_auc), lwd=0.5, fill="blue", alpha=0.07)  +
  geom_hline(yintercept = 0.5, col="grey60", lty=2, lwd=1) +  # theoretical values
  geom_hline(yintercept = 0.71, col="grey60", lty=3, lwd=1) + 
  geom_vline(xintercept = 0, col="grey20", lty=2, lwd=1) + # theoretical values
  geom_vline(xintercept = -0.07, col="grey20", lty=3, lwd=1)   ## remove opposite + some extreme
tt_p1

## add potential invaders
p_occu_cbi_pot=
  ggplot() + theme_bw() + labs(y=paste0("cbi.val.min"), x="source occurrence points", col="invader transfer") + 
  theme(legend.position = "top", text = element_text(size=22), 
        legend.text = element_text(size=22),
        axis.text = element_text(size=18)) +
  annotate('rect', xmin=100, xmax=2000, ymin=0, ymax=1, fill='lightcyan', alpha=0.5) +
  geom_point(data=info_pot, aes(y=cbi.val.Pmin, x=spp_scr_occ), col="grey", cex=2.5, alpha=0.4) + 
  geom_point(data=info_inv, aes(y=cbi.val.Pmin, x=spp_scr_occ, col=transfer), cex=2.5) + 
  geom_hline(yintercept = 0, col="grey60", lty=2, lwd=1)  + ## remove opposite + some extreme
  geom_vline(xintercept = 100, col="grey30", lty=3, lwd=1) 
p_occu_cbi_pot

info_inv %>% filter(spp_scr_occ < 100 & cbi.val.Pmin < 0.2 & transfer == "expected") %>% select(species)

leg=get_legend(p_occu_cbi_pot)
midpanel = ggarrange(p_occu, tt_p1,  nrow=2, ncol=1, legend = "none", labels=c("b", "c")) 
midpanel 

ptransftrain = ggarrange(p3maintext, midpanel, p_occu_cbi_pot, nrow=1, ncol=3, 
              legend.grob = leg, labels = c("a", "", "d"), widths=c(1.2,1,1.2))
ptransftrain
plot_pdf_png(ptransftrain, dirff_datvar_trg_fM, "F3_TransferTrainingRelations_AUC_Occ_minCBI_potential", pointsize = 4, pdf_width=21)

##
ptransftrain = ggarrange(p3maintext, p_occu_cbi_pot, nrow=1, ncol=2, 
                       legend.grob = leg, labels = c("a","b"), 
                       font.label = list(size = 20), widths=c(1,1))
ptransftrain
plot_pdf_png(ptransftrain, dirff_datvar_trg_fM, "F3_TransferTraining", pointsize = 4, pdf_width=14)

### post-modelling filters, numbers related to each criteria
source_sqm_thr = 9700000
spp_scr_occ_thr = 100
cbi.val.Pmin_thr = 0  #-0.07
auc.diff.sd_thr =  0.13

filters = c("sqm", "occ", "cbi", "auc", "sqm_occ", "cbi_auc", "occ_auc", "occ_cbi", "occ_auc_cbi")  
filter_val = c(source_sqm_thr, spp_scr_occ_thr, cbi.val.Pmin_thr, auc.diff.sd_thr, "", "", "", "", "")

filterout = info_all %>% 
  mutate(filter_sqm = ifelse(source_sqm <= source_sqm_thr, 1, NA)) %>% 
  mutate(filter_occ = ifelse(spp_scr_occ <= spp_scr_occ_thr, 1, NA)) %>% 
  mutate(filter_cbi = ifelse(cbi.val.Pmin <= cbi.val.Pmin_thr, 1, NA)) %>% 
  mutate(filter_auc = ifelse(auc.diff.sd >= auc.diff.sd_thr, 1, NA)) %>% 
  mutate(filter_sqm_occ = ifelse((source_sqm <= source_sqm_thr | spp_scr_occ <= spp_scr_occ_thr), 1, NA)) %>% 
  mutate(filter_cbi_auc = ifelse((auc.diff.sd >= auc.diff.sd_thr | cbi.val.Pmin <= cbi.val.Pmin_thr), 1, NA)) %>% 
  mutate(filter_occ_auc = ifelse((spp_scr_occ <= spp_scr_occ_thr | auc.diff.sd >= auc.diff.sd_thr ), 1, NA)) %>% 
  mutate(filter_occ_cbi = ifelse((spp_scr_occ <= spp_scr_occ_thr | cbi.val.Pmin <= cbi.val.Pmin_thr), 1, NA)) %>% 
  mutate(filter_occ_auc_cbi = ifelse((spp_scr_occ <= spp_scr_occ_thr | auc.diff.sd >= auc.diff.sd_thr | cbi.val.Pmin <= cbi.val.Pmin_thr), 1, NA)) 

fout_inv = filterout %>% filter(sp_status == "invader") %>% 
  select(species, filter_sqm, filter_occ, filter_cbi, filter_auc, filter_sqm_occ, filter_cbi_auc, filter_occ_auc, filter_occ_cbi, filter_occ_auc_cbi)  
fout_inv_n = fout_inv %>% select(-species) %>% colSums(na.rm=T)
fout_inv_p = fout_inv_n/nrow(info_inv)

fout_pot = filterout %>% filter(sp_status == "potential") %>% 
  select(species, filter_sqm, filter_occ, filter_cbi, filter_auc, filter_sqm_occ, filter_cbi_auc, filter_occ_auc, filter_occ_cbi, filter_occ_auc_cbi)  
fout_pot_n = fout_pot %>% select(-species) %>% colSums(na.rm=T)
fout_pot_p = fout_pot_n/nrow(info_pot)

all_filters = data.frame(filters, filter_val,
             fout_inv_n, fout_inv_p, fout_pot_n, fout_pot_p) %>% 
    mutate(fin_inv_n = n_inv - fout_inv_n) %>% mutate(fin_inv_p = 1 - fout_inv_p) %>% 
    mutate(fin_pot_n = n_pot - fout_pot_n) %>% mutate(fin_pot_p = 1 - fout_pot_p) 
all_filters
all_filters %>% filter(filters == "occ_cbi")


##### P2: whole country suitability and prefecture heatmap:
#require(patchwork)

model_clean_criteria="fc2"
sp_status="potential"
mask="both"

## just for this file we need the prefecture level:
div_level="prefecture"
dirbm_datvar_st=file.path(dirbm_datvar_trg, div_level, model_clean_criteria, sp_status)   
suit_pref = read.csv(file.path(dirbm_datvar_st, paste0("Species-",div_level,"_MF-",model_clean_criteria,"_mask-",mask,".csv")), row.names = 1) %>%
  rename(!!div_level:=subregion)

## then we change to cluster level:
div_level="cluster"
dirbm_datvar_st=file.path(dirbm_datvar_trg, div_level, model_clean_criteria, sp_status)   
richness = rast(file.path(dirbm_datvar_st, paste0("richness","_MF-",model_clean_criteria,"_mask-",mask,".tif")))
suit_clust = read.csv(file.path(dirbm_datvar_st, paste0("Species-",div_level,"_MF-",model_clean_criteria,"_mask-",mask,".csv")), row.names = 1) %>%
  rename(!!div_level:=subregion)

pref_dendro_data = readRDS(file.path(dirbm_datvar_trg_cl, "pref_dendro_data.rds"))
spp_dendro_data = readRDS(file.path(dirbm_datvar_trg_cl, "spp_dendro_data.rds"))
pref_clustCols=read.csv(file.path(dirbm_datvar_trg_cl,"PrefectureClusterColor.csv"), row.names = 1)
spp_groupCols = read.csv(file.path(dirbm_datvar_trg_cl,"SpeciesGroupColor.csv"), row.names = 1)
clust_col = pref_clustCols %>% distinct(cluster, cl_colour)
group_col = spp_groupCols %>% distinct(group, gr_colour)
trg_div2 = trg_div %>% full_join(pref_clustCols, by = join_by(prefecture)) 

richness_plot = 
  ggplot() + theme_bw() +
  geom_spatvector(data = trg_div2, fill = "grey85", col="grey90") +
  geom_spatraster(data = richness) +
  # geom_spatvector(data = spp_occ_target, fill = NA, col="darkorchid1", cex=0.2, alpha=0.5) +
  scale_fill_whitebox_c(palette = "deep") +
  labs(fill = "number of species") + 
  theme(legend.position = "top", legend.title = element_text(size=11), legend.text = element_text(size=9),
        legend.key.height= unit(1, 'cm'), legend.key.width= unit(1.7, 'cm'),
        panel.background = element_rect(fill='transparent')) 
richness_plot
plot_pdf_png(richness_plot, dirff_datvar_trg_fM, "F4_1_richness_countrywide", pointsize = 4, pdf_width=7, pdf_height = 7, bg="transparent")

#r2 = richness_plot + xlim(123, 147)
#plot_pdf_png(r2, dirff_datvar_trg_fM, "F4_1_richness_plot_2", pointsize = 4, pdf_width=7, pdf_height = 7, bg="transparent")

minmax(richness)

### 2023-09-27 : try to split cluster 1 to plot separately
## Jamie thinks it is not easy to see that prediction.

richness   ## the country-wide raster
## there is an island of Ogasawara which is so far that no matter the resolution
## never gets to plot nicely, so we will crop it out for visual purposes
richness2 = crop(richness, c(123, 147, 24, 46))  
plot(richness2)

## first remove the cluster I want to plot separately from the rest of the country
c1c = trg_div2 %>% filter(cluster != "cluster_1") 
rich_c1c = crop(richness2, c1c, mask=T)  
richness_plot = 
  ggplot() + theme_bw() +
  geom_spatvector(data = c1c, fill = "grey85", col="grey90") +
  geom_spatraster(data = rich_c1c) +
  # geom_spatvector(data = spp_occ_target, fill = NA, col="darkorchid1", cex=0.2, alpha=0.5) +
  scale_fill_whitebox_c(palette = "deep") +
  labs(fill = "number of species") + 
  theme(legend.position = "top", legend.title = element_text(size=11), legend.text = element_text(size=9),
        legend.key.height= unit(1, 'cm'), legend.key.width= unit(1.7, 'cm'),
        panel.background = element_rect(fill='transparent')) 
richness_plot
plot_pdf_png(richness_plot, dirff_datvar_trg_fM, "F4_1_richness_c1complement", pointsize = 4, pdf_width=7, pdf_height = 7, bg="transparent")

richness2 = crop(richness, c(123, 143, 24, 46))  
c1 = trg_div2 %>% filter(cluster == "cluster_1")  ## cluster I want to split
rich_c1 = crop(richness2, c1, mask=T)  
richness_plot = 
  ggplot() + theme_bw() +
 # geom_spatvector(data = c1, fill = "grey85", col="grey90") +
  geom_spatraster(data = rich_c1) +
  # geom_spatvector(data = spp_occ_target, fill = NA, col="darkorchid1", cex=0.2, alpha=0.5) +
  scale_fill_whitebox_c(palette = "deep") +
  labs(fill = "number of species") + 
  theme(legend.position = "top", legend.title = element_text(size=11), legend.text = element_text(size=9),
        legend.key.height= unit(1, 'cm'), legend.key.width= unit(1.7, 'cm'),
        panel.background = element_rect(fill='transparent')) 
richness_plot
plot_pdf_png(richness_plot, dirff_datvar_trg_fM, "F4_1_richness_c1", pointsize = 4, pdf_width=7, pdf_height = 7, bg="transparent")

#plot(rich_c1)
for(pref in c1$prefecture){
  print(pref)
  ppp = trg_div2 %>% filter(prefecture == pref)  ## cluster I want to split
  rich_ppp = crop(richness2, ppp, mask=T)  

richness_plot = 
  ggplot() + theme_bw() +
#  geom_spatvector(data = ppp, fill = "grey85", col="grey90") +
  geom_spatraster(data = rich_ppp) +
  # geom_spatvector(data = spp_occ_target, fill = NA, col="darkorchid1", cex=0.2, alpha=0.5) +
  scale_fill_whitebox_c(palette = "deep") +
  labs(fill = "number of species") + 
  theme(legend.position = "top", legend.title = element_text(size=11), legend.text = element_text(size=9),
        legend.key.height= unit(1, 'cm'), legend.key.width= unit(1.7, 'cm'),
        panel.background = element_rect(fill='transparent')) 
richness_plot
plot_pdf_png(richness_plot, dirff_datvar_trg_fM, paste0("F4_1_richness_c1_",pref), pointsize = 4, pdf_width=7, pdf_height = 7, bg="transparent")
}


## the clusters with different colour each
clust_jap = 
  ggplot() + theme_bw() + 
  theme(legend.position = c(0.15, 0.78), legend.title = element_blank(), 
     #   legend.text = element_text(size=20), legend.key.size = unit(5,"line"),
        axis.text = element_blank(), 
        panel.background = element_rect(fill='transparent')) +
  geom_spatvector(data = trg_div2, col="grey20") + 
  geom_spatvector(data = trg_div2, mapping = aes(fill = cluster)) + 
  scale_fill_manual(values = clust_col$cl_colour) +
  guides(fill=guide_legend(nrow=4, byrow=TRUE, reverse=T))
clust_jap 
plot_pdf_png(clust_jap, dirff_datvar_trg_fM, "F4_2_clust_jap", pointsize = 4, pdf_width = 7, bg="transparent")

regs_jap = 
  ggplot() + theme_bw() + 
  theme(legend.position = c(0.15, 0.78), legend.title = element_blank(), 
     #   legend.text = element_text(size=20), legend.key.size = unit(5,"line"),
        axis.text = element_blank(), 
        panel.background = element_rect(fill='transparent')) +
#  geom_spatvector(data = trg_div2, col=region_official) + 
  geom_spatvector(data = trg_div2, mapping = aes(fill = region))# + 
 # scale_fill_manual(values = clust_col$cl_colour) #+
 # guides(fill=guide_legend(nrow=4, byrow=TRUE, reverse=T))
regs_jap 

# richness_clusters = richness_plot + inset_element(clust_jap, left = 0.02, bottom = 0.54, right = 0.32, top = 1.01)
# richness_clusters
# plot_pdf_png(richness_clusters, dirff_datvar_trg_fM, "richness_clusters", pointsize = 4, pdf_width=7, pdf_height = 7)

### HERE: I already like the figures up to this point.

pref_dendro_ggdendrogram = 
  ggdendrogram(data = pref_dendro_data, rotate = TRUE) + 
  theme(axis.text.y = element_text(size = 10, colour = pref_dendro_data$labels$cl_colour))
pref_dendro_ggdendrogram
plot_pdf_png(pref_dendro_ggdendrogram, dirff_datvar_trg_fM, "F4_3_pref_dendro_ggdendrogram", pointsize = 4, pdf_width=7, pdf_height = 7)

spp_dendro_ggdendrogram = 
  ggdendrogram(data = spp_dendro_data) + 
  theme(axis.text.x = element_text(size = 10)) #, colour = spp_dendro_data$labels$gr_colour))
spp_dendro_ggdendrogram
plot_pdf_png(spp_dendro_ggdendrogram, dirff_datvar_trg_fM, "F4_4_spp_dendro_ggdendrogram", pointsize = 4, pdf_width=14, pdf_height = 7)

## species - prefectures
## ordered as prefectures and species as the dendrograms and add color info
suit_pref_plotOrd = suit_pref %>% 
  mutate(species = factor(species, levels = spp_dendro_data$labels$species)) %>%
  mutate(prefecture = factor(prefecture, levels=pref_dendro_data$labels$prefecture))

gSS_pref=
  ggplot(data=suit_pref_plotOrd, aes(y = prefecture, x = species)) +
  #theme() + 
  ylab(NULL) + xlab(NULL) +
  theme(legend.position = "bottom", legend.text = element_text(size=9),
        legend.key.height= unit(0.8, 'cm'), legend.key.width= unit(1.7, 'cm'),
        axis.text.y = element_text(size=9, colour = pref_dendro_data$labels$cl_colour), 
        axis.text.x = element_text(size=9, angle = 90, hjust=0, vjust=2)) + #, colour = spp_dendro_data$labels$gr_colour
  scale_y_discrete(position = "right") + scale_x_discrete(position = "top") +
  geom_tile(aes(fill=mean)) + scale_fill_whitebox_c("mean suitability", palette = "deep")
gSS_pref
plot_pdf_png(gSS_pref, dirff_datvar_trg_fM, "F4_5_gSS_pref_names", pointsize = 4, pdf_width=(7*1.5), pdf_height = 7)

gSS_pref=
  ggplot(data=suit_pref_plotOrd, aes(y = prefecture, x = species)) +
  #theme() + 
  ylab(NULL) + xlab(NULL) +
  theme(legend.position = "bottom", legend.text = element_text(size=9),
        legend.key.height= unit(0.8, 'cm'), legend.key.width= unit(1.7, 'cm'),
        axis.text.y = element_blank(),
        axis.text.x = element_blank()) + 
  scale_y_discrete(position = "right") + scale_x_discrete(position = "top") +
  geom_tile(aes(fill=mean)) + scale_fill_whitebox_c("mean suitability", palette = "deep")
gSS_pref
plot_pdf_png(gSS_pref, dirff_datvar_trg_fM, "F4_5_gSS_pref_NOnames", pointsize = 4, pdf_width=14, pdf_height = 7)


# library(heatmaply)
# test01 = suit_pref_plotOrd %>% select(-min, -max) %>% pivot_wider(names_from = prefecture, values_from = mean) %>% 
#   as.data.frame() 
# rownames(test01) = test01$species
# test01 = test01 %>% select(-species) %>% as.matrix() %>% t()
# test01[is.na(test01)] <- as.double("NA")
# 
# #heatmapr(suit_pref_plotOrd)
# ggheatmap(test01)
# 
# library(pheatmap)
# pheatmap(test01, na_col = "grey", cluster_rows=F, cluster_cols=F)
# 
# library(gplots)
# heatmap.2(test01)#, trace="none", na.color = "black", scale="none")
# #gSS_pref %>% insert_left(pref_dendro_ggdendrogram)
# 
# #ggarrange(richness_plot, gSS_pref, pref_dendro_ggdendrogram,  widths = c(1.6, 0.7, 0.7), ncol=3)
# 
# #plot_pdf_png(, dirff_datvar_trg_fM, "Richness_Heatmap", pointsize = 4, pdf_width=7)


##### P3: source plot: each cluster one at the time
suit_clust2 = suit_clust %>% mutate(subregion = cluster)
suit_clust2 = countsForPlot_SubregSpecies(suit_clust2)
suit_clust2 %>% group_by(cluster) %>% summarise(mean(mean)) 
suit_clust2 %>% group_by(species) %>% summarise(n()) 
suit_clust2 %>% select(species, cluster) %>% group_by(cluster) %>% summarise(n()) 

source_plots = list()
for(j in seq(1:nrow(clust_col))){
  print(j)
  
  clustj_source_rast=rast(file.path(dirbm_datvar_st, paste0("Source_", div_level,"-",clust_col$cluster[j],"_MF-",model_clean_criteria,"_mask-",mask,".tif")))
#  clustj_source = species_source_plot(clustj_source_rast, "weighted richness")
  raster_stack = clustj_source_rast
  labname = "weighted richness"
  
  clustj_source = ggplot() + theme_bw() + scale_fill_whitebox_c(palette = "deep") + labs(fill = labname) + 
    theme(legend.position = c(0.65, 0.075), legend.direction = "horizontal", 
          legend.title = element_text(size=16, vjust=0.8), 
          legend.text = element_text(size=15),
          legend.key.height= unit(0.7, 'cm'), legend.key.width= unit(2.5, 'cm')) +
    geom_spatvector(data = wcl, col="grey85") +
   # geom_spatvector(data = trg_jp, col="grey20", fill="grey85") + 
    geom_spatraster(data = raster_stack)
  clustj_source 
  
  clustj_jap = 
    ggplot() + theme_bw() + 
    theme(legend.position = c(0.5, 1.15), axis.text = element_blank(), plot.background = NULL,
          legend.text = element_text(size=16), legend.title = element_blank()) +
    geom_spatvector(data = trg_div2, col="grey20") + 
    geom_spatvector(data = subset(trg_div2, cluster == clust_col$cluster[j]), mapping = aes(fill = cluster), col="grey20") + 
    scale_fill_manual(values = clust_col$cl_colour[j]) #+
#  clustj_jap
  
  gg3 = clustj_source + inset_element(clustj_jap, left = 0.0, bottom = 0.01, right = 0.27, top = 0.45)
#  gg3
  source_plots[[j]] = gg3
}
#gg3
cluster_sources = ggarrange(plotlist = source_plots, ncol= 2, nrow=2)
cluster_sources 
plot_pdf_png(cluster_sources, dirff_datvar_trg_fM, "SourceHotspots_PClusters", pointsize = 4, pdf_width=(7*1.5*2), pdf_height = (7*1.8))





############
print("owaru")
print("bye bye")


##### code end  




