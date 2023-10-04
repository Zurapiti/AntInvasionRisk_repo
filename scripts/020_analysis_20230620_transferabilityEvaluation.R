## 
## author:  Yazmin Zurapiti
## title:   020_analysis_20230620_transferabilityEvaluation.R
## aim:     Here I analyse what variables are the most relevant to explain good or bad transfers
## version: 3.0
## call:  lauch from /flash/EconomoU/Zurapiti as there is where the R packages are installed
##        020_analysis_20230620_transferabilityEvaluation.R [save_date] [trg] [bioset] [algor] 
## inputs: 
##   1.inputdir: 
##   1.inputfile:
## outputs: 
##   1.outputdir: 
##   1.outputfile: 
## 

# Load packages -- the order here is important because some pkg functions overwrite others.
library(sp)
library(sf)
library(terra)
library(raster)
library(rnaturalearth)
library(ggplot2)
library(ggpubr)

# library(stringi)  ## to make a matrix out of a list of uneven number of rows 
# library(stringr)  ## for leading zeros #str_pad(lay, 2, pad = "0")
# library(dismo)
# library(ENMeval)
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

## get parameters
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
#dirbm_datvar_tsp=file.path(dirbm_datvar_trg, sp_status, this_species)   ## species directory (only declared here)

## save directories:
dirfm_datvar_trg=file.path(dir_flash_models, datavar, trg)
#dirfm_datvar_tsp=file.path(dirfm_datvar_trg, sp_status, this_species)   ## species directory (only declared here)
if(!dir.exists(dirfm_datvar_trg)){dir.create(dirfm_datvar_trg, recursive = T)}

dirff_datvar_trg=file.path(dir_flash_figures, datavar, trg)
dirff_datvar_trg_trans = file.path(dirff_datvar_trg, "transferabilityEvaluation")
#dirff_datvar_tsp=file.path(dirff_datvar_trg, sp_status, this_species)   ## species directory (only declared here)
if(!dir.exists(dirff_datvar_trg_trans)){dir.create(dirff_datvar_trg_trans, recursive = T)}

## functions

lmfitmydata = function(minimal_data){
  colnames(minimal_data) = c("x", "y")
  minimal_data %>% na.omit()
  lm=lm(y ~ x, data=minimal_data)
  print(summary(lm))
  minimal_data = cbind(minimal_data, predict(lm, interval = 'confidence') )
  return(minimal_data)
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
print("read data")

# trg_div=readRDS(file = file.path(dirbd_env_trg, "Subregions.rds")) ## sf class
# trg_silh = trg_div %>% vect()  ## terra class, for plotting
# block_plot=
#   ggplot() + theme_bw() +
#   geom_spatvector(data = trg_div, mapping=aes(fill=block)) + #, fill = "region", col="grey90") +
#   theme(legend.position = "top", legend.title = element_text(size=12), legend.text = element_text(size=12),
#         legend.key.height= unit(0.6, 'cm'), legend.key.width= unit(0.6, 'cm'))
# plot_pdf_png(block_plot, dirff_datvar_trg, "Japan_blocks", pointsize = 4)

names_nat = read.csv(file=file.path(dirbd_env_trg, "native.csv"), header = F)
colnames(names_nat) = "species"
names_inv = read.csv(file=file.path(dirbd_env_trg, "invader.csv"), header = F)
colnames(names_inv) = "species"
names_pot = read.csv(file=file.path(dirbd_env_trg, "potential.csv"), header = F)
colnames(names_pot) = "species"

n_wexo_spp = nrow(names_nat) + nrow(names_inv) + nrow(names_pot)
n_ana_spp = nrow(names_inv) + nrow(names_pot)
n_inv_spp = nrow(names_inv) 
n_pot_spp = nrow(names_pot)
print(paste0("total number of alien ant species: ", n_wexo_spp))
print(paste0("total number of analyzed species: ", n_ana_spp))
print(paste0("number of analyzed invader species: ", n_inv_spp))
print(paste0("number of analyzed exotic species: ", n_pot_spp))

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

# recorddata
# bm_metadata
# bm_src_boyce
# bm_transdata
##########################
## some basic data

info_all = recorddata %>% full_join(bm_metadata, by=c("species", "sp_status")) %>% 
  full_join(bm_transdata, by=c("species", "sp_status"))

info_inv = info_all %>% filter(species %in% names_inv$species)
info_pot = info_all %>% filter(species %in% names_pot$species)

print(paste0("number of species in models files: ", nrow(info_all)))
print(paste0("number of invaders in models files: ", nrow(info_inv)))
print(paste0("number of exotics in models files: ", nrow(info_pot)))

print(paste0("total number of records for all species: ", sum(info_all$spp_tot_occ)))
print(paste0("total number of records for invader species: ", sum(info_inv$spp_tot_occ)))
print(paste0("total number of records for exotic species: ", sum(info_pot$spp_tot_occ)))

#######################
## now our evaluation on how well these transfers happened:

oppo = c("Tetramorium_smithi", "Cardiocondyla_itsukii", "Tetramorium_kraepelini", "Cardiocondyla_kagutsuchi")
over = c("Ooceraea_biroi", "Pheidole_parva", "Hypoponera_ragusai" ,"Nylanderia_amia", "Tetraponera_allaborans")
under = c("Tapinoma_melanocephalum", "Tetramorium_bicarinatum")

info_inv = info_inv %>%  
  mutate(transfer=ifelse(species %in% over, "extreme", 
                                ifelse(species %in% oppo, "opposite", "expected"))) 

## artificially change the NA values for the cbi for 0 to be able to visualize them
info_inv = info_inv %>%  
  mutate(cbi_edit = ifelse(is.na(ext_trg_boyce), 0, ext_trg_boyce)) %>% 
  mutate(transfer_edit = ifelse(species %in% under, "underpredicted",
                                ifelse(is.na(ext_trg_boyce), paste0(transfer,"_cbiNA"), transfer)))

## run 2023-06-23 has an extra species which is overpredicted. ## previous values were: 
## transfer variable selection: ext_trg_boyce ~ cbi_edit = 0.42, trg_cor = 0.08, trg_auc = 0.72
## training variable selection: filter(auc.diff.sd >= 0.13 | cbi.val.Pmin <= -0.2 ) 

### here I evaluate what are the relationships between the transfer metrics:
### we can also try to find values that distinguish the "good" from the "bad" predictions.

#### the cbi edited version, no NA 
d = info_inv %>% select(species, transfer, transfer_edit, cbi_edit, trg_cor) #%>% rename(x=cbi_edit, y=trg_cor) #%>% na.omit()
#lm=lm(y ~ x, data=d)
#summary(lm)
#d = cbind(d, predict(lm, interval = 'confidence') )
t1=
  ggplot() + theme_bw() + labs(x="transfer CBI", y="transfer correlation", col="transfer") + 
  theme(legend.position = "bottom", text = element_text(size=16)) + 
  geom_point(data=info_inv, aes(x=cbi_edit, y=trg_cor, col=transfer_edit), cex=2.5) +  # col=transfer removes the cbi = NA 
  #  geom_line(data=d, aes(x, fit), col="blue2") +
  #  geom_ribbon(data=d, aes(x,ymin=lwr,ymax=upr), fill="blue2", alpha=0.3) +
  geom_hline(yintercept = 0.0, col="grey60", lty=2, lwd=1) +   ## theoretical value
  geom_vline(xintercept = 0.0, col="grey60", lty=2, lwd=1) +   ## theoretical value
  geom_hline(yintercept = 0.08 , col="grey20", lty=3, lwd=1) + 
  geom_vline(xintercept = -0.42, col="grey20", lty=3, lwd=1) 
t1
plot_pdf_png(t1, dirff_datvar_trg_trans, "transferEval_EcospatCBI_DismoCor_NoNA", pointsize = 4)

## no NA in cbi
d = info_inv %>% select(species, transfer, transfer_edit, cbi_edit, trg_auc) #%>% rename(x=cbi_edit, y=trg_auc) #%>% na.omit()
#lm=lm(y ~ x, data=d)
#summary(lm)
#d = cbind(d, predict(lm, interval = 'confidence') )
t2=
  ggplot() + theme_bw() + labs(x="transfer CBI", y="transfer AUC", col="transfer") + 
  theme(legend.position = "bottom", text = element_text(size=16)) + 
  geom_point(data=info_inv, aes(x=cbi_edit, y=trg_auc, col=transfer_edit), cex=2.5)  +
  #  geom_line(data=d, aes(x, fit), col="blue2") +
  #  geom_ribbon(data=d, aes(x,ymin=lwr,ymax=upr), fill="blue2", alpha=0.3)  +
  geom_hline(yintercept = 0.5, col="grey60", lty=2, lwd=1) +   ## theoretical value
  geom_vline(xintercept = 0.0, col="grey60", lty=2, lwd=1) +   ## theoretical value
  geom_hline(yintercept = 0.71, col="grey20", lty=3, lwd=1) +
  geom_vline(xintercept = -0.42, col="grey20", lty=3, lwd=1)
t2
plot_pdf_png(t2, dirff_datvar_trg_trans, "transferEval_EcospatCBI_DismoAUC_NoNA", pointsize = 4)

####
d = info_inv %>% select(species, transfer, transfer_edit, trg_cor, trg_auc) #%>% rename(x=trg_cor, y=trg_auc) #%>% na.omit()
#lm=lm(y ~ x, data=d)
#summary(lm)
#d = cbind(d, predict(lm, interval = 'confidence') )
t3=
  ggplot() + theme_bw() + labs(x="transfer correlation", y="transfer AUC", col="transfer") + 
  theme(legend.position = "bottom", text = element_text(size=16)) + 
  geom_point(data=info_inv, aes(trg_cor, trg_auc, col=transfer_edit), cex=2.5)  +
#  geom_line(data=d, aes(x, fit), col="blue2") +
#  geom_ribbon(data=d, aes(x,ymin=lwr,ymax=upr), fill="blue2", alpha=0.1) +
  geom_hline(yintercept = 0.5, col="grey60", lty=2, lwd=1) +  # theoretical values
  geom_vline(xintercept = 0.0, col="grey60", lty=2, lwd=1) +  # theoretical values
  geom_hline(yintercept = 0.71, col="grey20", lty=3, lwd=1) + 
  geom_vline(xintercept = 0.08, col="grey20", lty=3, lwd=1) #+
   # ylim(c(-0.05,1.2))
t3
plot_pdf_png(t3, dirff_datvar_trg_trans, "transferEval_DismoCor_DismoAUC_NoNA", pointsize = 4)

## plot for supplementary information
t=ggarrange(t1,t2,t3, nrow=1, ncol=3, common.legend = T, labels=c("a", "b", "c"))
t
plot_pdf_png(t, dirff_datvar_trg_trans, "transferEval_metricsCompare", pointsize = 4, pdf_width=21)

p3maintext=
  ggplot() + theme_bw() + labs(x="(transfer) cor", y="(transfer) AUC", col="transfer") + 
  theme(legend.position = "bottom", text = element_text(size=16)) + 
  geom_point(data=info_inv, aes(trg_cor, trg_auc, col=transfer), cex=2.5)  +
  #  geom_line(data=d, aes(x, fit), col="blue2") +
  #  geom_ribbon(data=d, aes(x,ymin=lwr,ymax=upr), fill="blue2", alpha=0.1) +
  geom_hline(yintercept = 0.5, col="grey60", lty=2, lwd=1) +  # theoretical values
  geom_vline(xintercept = 0.0, col="grey60", lty=2, lwd=1) +  # theoretical values
  geom_hline(yintercept = 0.71, col="grey20", lty=3, lwd=1) + 
  geom_vline(xintercept = 0.08, col="grey20", lty=3, lwd=1) #+
# ylim(c(-0.05,1.2))
p3maintext

## because the CBI = NA models were expected patterns, 
## we can use "transfer", instead of "transfer_edit" for further plots

#plot_pdf_png(p, dirff_datvar_trg, "transferMetrics-Relationship", pointsize = 4)
info_inv %>% filter(trg_cor <= 0.08) %>% select(species, trg_cor)
info_inv %>% filter(trg_auc <= 0.71) %>% select(species, trg_auc)

Ndistinguish=length(oppo)+length(over)
Ntot=info_inv %>% nrow()

length(oppo)/Ntot
Ndistinguish/Ntot
(Ntot - Ndistinguish)/Ntot

###############################################################################
### evaluating the transfer against the raw data

### here it seems like the transfer metrics are not very correlated to the 
## raw data, which should be a good sign, as that means that 
## we cannot predict a "bad" transfer by the number of source points or area occupied.

# let's plot the two slightly more correlated
d = info_inv %>% select(species, transfer, source_sqm, spp_scr_occ, trg_auc) 
#d = lmfitmydata(d)
p_area=
  ggplot() + theme_bw() + labs(y="transfer AUC", x="source area", col="transfer") + 
  theme(legend.position = "bottom", text = element_text(size=16)) +
  geom_point(data=info_inv, aes(y=trg_auc, x=source_sqm, col=transfer), cex=2.5)  +
  geom_smooth(data=info_inv, aes(y=trg_auc, x=source_sqm), lwd=0.5, fill="blue", alpha=0.07)  +
#  geom_line(data=d, aes(x, fit), col="blue2") +
#  geom_ribbon(data=d, aes(x, ymin=lwr, ymax=upr), fill="blue2", alpha=0.1) +
  geom_hline(yintercept = 0.5, col="grey60", lty=2, lwd=1) +  # theoretical values
  geom_hline(yintercept = 0.71, col="grey60", lty=3, lwd=1) + 
  geom_vline(xintercept = 9700000, col="grey20", lty=5, lwd=1) 
p_area
plot_pdf_png(p_area, dirff_datvar_trg_trans, "sourceProp_sourceArea", pointsize = 4)

## perhaps the number of source points has some effects

#d = info_inv %>% select(species, transfer, spp_scr_occ, trg_auc) 
#d = lmfitmydata(d)
p_occu=
  ggplot() + theme_bw() + labs(y="(transfer) AUC", x="(training) occurrence points", col="transfer") + 
  theme(legend.position = "bottom", text = element_text(size=16)) +
  geom_point(data=info_inv, aes(y=trg_auc, x=spp_scr_occ, col=transfer), cex=2.5) + 
  geom_smooth(data=info_inv, aes(y=trg_auc, x=spp_scr_occ), span=0.9, lwd=0.5, fill="blue", alpha=0.07)  +
#  geom_line(data=d, aes(x, fit), col="blue2") +
#  geom_ribbon(data=d, aes(x,ymin=lwr,ymax=upr), fill="blue2", alpha=0.1) +
  geom_hline(yintercept = 0.5, col="grey60", lty=2, lwd=1) +  # theoretical values
  geom_hline(yintercept = 0.71, col="grey60", lty=3, lwd=1) + 
  #geom_vline(xintercept = 20, col="grey20", lty=4, lwd=1) +
  geom_vline(xintercept = 100, col="grey20", lty=3, lwd=1) 
p_occu
plot_pdf_png(p_occu, dirff_datvar_trg_trans, "sourceProp_SourceOccurrencePoints_NoNA", pointsize = 4)
## zoom into the lower range of this plot
q = p_occu + xlim(20, 125) 
q
plot_pdf_png(q, dirff_datvar_trg_trans, "sourceProp_SourceOccurrencePoints_NoNA_LOW", pointsize = 4)

d = info_inv %>% select(species, transfer, spp_trg_occ, trg_auc) 
#d = lmfitmydata(d)
p_trgocc=
  ggplot() + theme_bw() + labs(y="transfer AUC", x="(tranfer) occurrence points", col="transfer") + 
  theme(legend.position = "bottom", text = element_text(size=16)) +
  geom_point(data=info_inv, aes(y=trg_auc, x=spp_trg_occ, col=transfer), cex=2.5) + 
  geom_smooth(data=info_inv, aes(y=trg_auc, x=spp_trg_occ), span=0.9, lwd=0.5, fill="blue", alpha=0.07)  +
  #  geom_line(data=d, aes(x, fit), col="blue2") +
  #  geom_ribbon(data=d, aes(x,ymin=lwr,ymax=upr), fill="blue2", alpha=0.1) +
  geom_hline(yintercept = 0.5, col="grey60", lty=2, lwd=1) +  # theoretical values
  geom_hline(yintercept = 0.71, col="grey60", lty=3, lwd=1)
p_trgocc
plot_pdf_png(p_trgocc, dirff_datvar_trg_trans, "sourceProp_TargetOccurrencePoints_NoNA", pointsize = 4)

### maaaybe there was something to see there in the end... 
## number of points seems to have a strong impact in the quality of the training transfer

#source_sqm <= 9700000 # remove
#spp_scr_occ <= 100 # remove

d = info_inv %>% select(spp_scr_occ, source_sqm) 
d = lmfitmydata(d)
p_occu_area=
  ggplot() + theme_bw() + labs(y="source area (sqm)", x="source occurrence points", col="transfer") + 
  theme(legend.position = "bottom", text = element_text(size=16)) +
  geom_point(data=info_inv, aes(y=source_sqm, x=spp_scr_occ, col=transfer), cex=2.5)  +
#  geom_smooth(data=info_inv, aes(y=source_sqm, x=spp_scr_occ), span=0.9, lwd=0.5, fill="blue", alpha=0.07) +
  geom_line(data=d, aes(x, fit), col="blue2") +
  geom_ribbon(data=d, aes(x,ymin=lwr,ymax=upr), fill="blue2", alpha=0.1) +
  geom_point(data=info_pot, aes(y=source_sqm, x=spp_scr_occ), col="grey", cex=2.5, alpha=0.3) +
  # geom_hline(yintercept = 0.5, col="grey60", lty=2, lwd=1)  # theoretical values
  geom_hline(yintercept = 9700000, col="grey60", lty=3, lwd=1) + 
  # geom_vline(xintercept = 20, col="grey20", lty=4, lwd=1) +
  geom_vline(xintercept = 100, col="grey60", lty=3, lwd=1) 

p_occu_area

###############################################################################
### evaluating the transfer against the model metrics

## cbi.val.Pmin 
d = info_inv %>% select(species, transfer, trg_auc, cbi.val.Pmin, auc.diff.sd) 
#d = lmfitmydata(d)
tt_p1 = 
ggplot() + theme_bw() + labs(x=paste0("(training) ", "cbi.val.Pmin"), y="(transfer) AUC", col="transfer") + 
  theme(legend.position = "bottom", text = element_text(size=16)) +
  geom_point(data=info_inv, aes(x=cbi.val.Pmin, y=trg_auc, col=transfer), cex=2.5)  +
  geom_smooth(data=info_inv, aes(x=cbi.val.Pmin, y=trg_auc), lwd=0.5, fill="blue", alpha=0.07)  +
  # geom_line(data=d, aes(x, fit), col="blue2") +
  # geom_ribbon(data=d, aes(x, ymin=lwr, ymax=upr), fill="blue2", alpha=0.1)+
  geom_hline(yintercept = 0.5, col="grey60", lty=2, lwd=1) +  # theoretical values
  geom_hline(yintercept = 0.71, col="grey60", lty=3, lwd=1) + 
  geom_vline(xintercept = 0, col="grey20", lty=2, lwd=1) + # theoretical values
  #geom_vline(xintercept = -0.2, col="grey0", lty=4, lwd=0.7) +  ## remove only opposite
  geom_vline(xintercept = -0.07, col="grey20", lty=3, lwd=1)   ## remove opposite + some extreme
tt_p1
plot_pdf_png(tt_p1, dirff_datvar_trg_trans, "traningRel_DismoAUC_PartMinCBI_noNA", pointsize = 4)

## auc.diff.sd 
tt_p2=
  ggplot() + theme_bw() + labs(x=paste0("(training) ", "auc.diff.sd"), y="(transfer) AUC", col="transfer") + 
  theme(legend.position = "bottom", text = element_text(size=16)) +
  geom_point(data=info_inv, aes(x=auc.diff.sd, y=trg_auc, col=transfer), cex=2.5)  +
  geom_smooth(data=info_inv, aes(x=auc.diff.sd, y=trg_auc), lwd=0.5, fill="blue", alpha=0.07)  +
  # geom_line(data=d, aes(x, fit), col="blue2") +
  # geom_ribbon(data=d, aes(x,ymin=lwr,ymax=upr), fill="blue2", alpha=0.1)+
  geom_hline(yintercept = 0.5, col="grey60", lty=2, lwd=1) +  # theoretical values
  geom_hline(yintercept = 0.71, col="grey60", lty=3, lwd=1) +
#  geom_vline(xintercept = 0.1, col="grey20", lty=2, lwd=1) + # "theoretical" values
  #geom_vline(xintercept = 0.165, col="grey0", lty=4, lwd=0.7) +  ## remove only opposite
  geom_vline(xintercept = 0.13, col="grey20", lty=3, lwd=1)   ## remove opposite + some extreme
tt_p2
## auc.diff.sd >= 0.13; 4 or 5 strange species
## auc.diff.sd >= 0.105; 5-6 strange species, but we remove one with good transfer
plot_pdf_png(tt_p2, dirff_datvar_trg_trans, "traningRel_DismoAUC_AUCDiffSD_NoNA", pointsize = 4)

ptt=ggarrange(p_occu, tt_p1, tt_p2, nrow=1, ncol=3, common.legend = T, labels=c("a", "b", "c"))
ptt
plot_pdf_png(ptt, dirff_datvar_trg_trans, "transferEval_RawModelMetrics", pointsize = 4, pdf_width=21)




##### combinations of the selected metrics

######
###  occurrence vs cbi
d = info_inv %>% select(species, transfer, auc.diff.sd, spp_scr_occ) 
#d = lmfitmydata(d)
p_occu_acu_inv=
  ggplot() + theme_bw() + labs(y=paste0("(training) ", "auc.diff.sd"), x="source occurrence points", col="transfer") + 
  theme(legend.position = "top", text = element_text(size=16)) +
  geom_point(data=info_inv, aes(y=auc.diff.sd, x=spp_scr_occ, col=transfer), cex=2.5) + 
  # geom_smooth(data=info_inv, aes(y=cbi.val.Pmin, x=spp_scr_occ), span=0.9, lwd=0.5, fill="blue", alpha=0.07)  +
  #  geom_point(data=info_pot, aes(y=cbi.val.Pmin, x=spp_scr_occ), col="grey", cex=2.5, alpha=0.4) + 
  #  geom_line(data=d, aes(x, fit), col="blue2") +
  #  geom_ribbon(data=d, aes(x,ymin=lwr,ymax=upr), fill="blue2", alpha=0.1) +
  geom_hline(yintercept = 0.13, col="grey60", lty=2, lwd=1)  + ## remove opposite + some extreme
  # geom_hline(yintercept = 0.5, col="grey60", lty=2, lwd=1) +  # theoretical values
  # geom_hline(yintercept = 0.71, col="grey60", lty=3, lwd=1) + 
  # geom_vline(xintercept = 20, col="grey20", lty=4, lwd=1) +
  geom_vline(xintercept = 100, col="grey30", lty=3, lwd=1) 
p_occu_acu_inv

p_occu_acu_pot=
  ggplot() + theme_bw() + labs(y=paste0("(training) ", "auc.diff.sd"), x="source occurrence points", col="transfer") + 
  theme(legend.position = "top", text = element_text(size=16)) +
  geom_point(data=info_inv, aes(y=auc.diff.sd, x=spp_scr_occ, col=transfer), cex=2.5) + 
  # geom_smooth(data=info_inv, aes(y=cbi.val.Pmin, x=spp_scr_occ), span=0.9, lwd=0.5, fill="blue", alpha=0.07)  +
  geom_point(data=info_pot, aes(y=auc.diff.sd, x=spp_scr_occ), col="grey", cex=2.5, alpha=0.4) + 
  #  geom_line(data=d, aes(x, fit), col="blue2") +
  #  geom_ribbon(data=d, aes(x,ymin=lwr,ymax=upr), fill="blue2", alpha=0.1) +
  geom_hline(yintercept = 0.13, col="grey60", lty=2, lwd=1)  + ## remove opposite + some extreme
  # geom_hline(yintercept = 0.5, col="grey60", lty=2, lwd=1) +  # theoretical values
  # geom_hline(yintercept = 0.71, col="grey60", lty=3, lwd=1) + 
  # geom_vline(xintercept = 20, col="grey20", lty=4, lwd=1) +
  geom_vline(xintercept = 100, col="grey30", lty=3, lwd=1) 
p_occu_acu_pot


###  occurrence vs cbi
d = info_inv %>% select(species, transfer, cbi.val.Pmin, spp_scr_occ) 
#d = lmfitmydata(d)
p_occu_cbi_inv=
  ggplot() + theme_bw() + labs(y=paste0("(training) ", "cbi.val.Pmin"), x="source occurrence points", col="transfer") + 
  theme(legend.position = "top", text = element_text(size=16)) +
  geom_point(data=info_inv, aes(y=cbi.val.Pmin, x=spp_scr_occ, col=transfer), cex=2.5) + 
  # geom_smooth(data=info_inv, aes(y=cbi.val.Pmin, x=spp_scr_occ), span=0.9, lwd=0.5, fill="blue", alpha=0.07)  +
  #  geom_point(data=info_pot, aes(y=cbi.val.Pmin, x=spp_scr_occ), col="grey", cex=2.5, alpha=0.4) + 
  #  geom_line(data=d, aes(x, fit), col="blue2") +
  #  geom_ribbon(data=d, aes(x,ymin=lwr,ymax=upr), fill="blue2", alpha=0.1) +
  geom_hline(yintercept = 0, col="grey60", lty=2, lwd=1)  + ## remove opposite + some extreme
  # geom_hline(yintercept = 0.5, col="grey60", lty=2, lwd=1) +  # theoretical values
  # geom_hline(yintercept = 0.71, col="grey60", lty=3, lwd=1) + 
  # geom_vline(xintercept = 20, col="grey20", lty=4, lwd=1) +
  geom_vline(xintercept = 100, col="grey30", lty=3, lwd=1) 
p_occu_cbi_inv

plot_pdf_png(p_occu_cbi_inv, dirff_datvar_trg_trans, "sourceTrainProp_SourceOccu_TrainCBI_inv_NoNA", pointsize = 4)
## zoom into the lower range of this plot
q = p_occu_cbi_inv + xlim(20, 150) 
q
plot_pdf_png(q, dirff_datvar_trg_trans, "sourceTrainProp_SourceOccu_TrainCBI_inv_NoNA_LOW", pointsize = 4)

## add potential invaders
p_occu_cbi_pot=
  ggplot() + theme_bw() + labs(y=paste0("(training) ", "cbi.val.Pmin"), x="source occurrence points", col="transfer") + 
  theme(legend.position = "top", text = element_text(size=16)) +
  #geom_smooth(data=info_inv, aes(y=cbi.val.Pmin, x=spp_scr_occ), span=0.9, lwd=0.5, fill="blue", alpha=0.07)  +
  geom_point(data=info_pot, aes(y=cbi.val.Pmin, x=spp_scr_occ), col="grey", cex=2.5, alpha=0.4) + 
  geom_point(data=info_inv, aes(y=cbi.val.Pmin, x=spp_scr_occ, col=transfer), cex=2.5) + 
  #  geom_line(data=d, aes(x, fit), col="blue2") +
  #  geom_ribbon(data=d, aes(x,ymin=lwr,ymax=upr), fill="blue2", alpha=0.1) +
  geom_hline(yintercept = 0, col="grey60", lty=2, lwd=1)  + ## remove opposite + some extreme
  geom_vline(xintercept = 100, col="grey30", lty=3, lwd=1) 

p_occu_cbi_pot
plot_pdf_png(p_occu_cbi_pot, dirff_datvar_trg_trans, "sourceTrainProp_SourceOccu_TrainCBI_pot_NoNA", pointsize = 4)
## zoom into the lower range of this plot
q = p_occu_cbi_pot + xlim(20, 150) 
q
plot_pdf_png(q, dirff_datvar_trg_trans, "sourceTrainProp_SourceOccu_TrainCBI_pot_NoNA_LOW", pointsize = 4)




######################### some numbers 

print(paste0("number of analyzed exotic species: ", n_pot_spp))

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

## here I also want to calculate how many of the expected vs others we get to filter out:
info_inv_vis = info_inv %>% select(species, transfer)
fout_inv_n = fout_inv %>% full_join(info_inv_vis) %>% group_by(transfer) %>% summarise(n = n())
fout_inv_tcat = fout_inv %>% full_join(info_inv_vis) %>% group_by(transfer) %>% 
summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) #%>% t()
fout_inv_tcat = full_join(fout_inv_n, fout_inv_tcat)

fout_pot = filterout %>% filter(sp_status == "potential") %>% 
  select(species, filter_sqm, filter_occ, filter_cbi, filter_auc, filter_sqm_occ, filter_cbi_auc, filter_occ_auc, filter_occ_cbi, filter_occ_auc_cbi)  
fout_pot_n = fout_pot %>% select(-species) %>% colSums(na.rm=T)
fout_pot_p = fout_pot_n/nrow(info_pot)


all_filters = data.frame(filters, filter_val,
                         fout_inv_n, fout_inv_p, fout_pot_n, fout_pot_p) %>% 
  mutate(fin_inv_n = n_inv_spp - fout_inv_n) %>% mutate(fin_inv_p = 1 - fout_inv_p) %>% 
  mutate(fin_pot_n = n_pot_spp - fout_pot_n) %>% mutate(fin_pot_p = 1 - fout_pot_p) 
all_filters

write.table(all_filters, file=file.path(dirfm_datvar_trg, "Transferability_Filters.csv"), col.names = T, row.names = F)
write.table(fout_inv_tcat, file=file.path(dirfm_datvar_trg, "Transferability_Filters_InvadersTransferCategories.csv"), col.names = T, row.names = F)

### conclusion, we will use the occ_cbi filter:
all_filters %>% filter(filters == "occ_cbi")



############
print("owaru")
print("bye bye")


##### code end  




