## 
## author:  Yazmin Zurapiti
## title:   trg_052_20230721_HierarchicalClusteringAnalysis.R
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
library(factoextra)
library(cluster)
#library(NbClust)

library(ggplot2)
library(ggpubr)
library(tidyterra)
library(ggdendro)

#library(ggtree)
#library(aplot)
#library(dendsort)
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
sp_status="potential" #"invader" # "prediction" # 
model_clean_criteria="fc2" # fc0, fc1, postmodeling filtering criteria, (occ/cbi/auc combinations)
div_level="cluster"  ## cluster # prefecture # block # parallel  ## the administrative level for the risk of establishment and source predictions
mask="both"  ## what mask to use for the prefecture prediction

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
#dirbm_datvar_tsp=file.path(dirbm_datvar_trg, sp_status, this_species)   ## species directory (only declared here)

## save directories:
dirfm_datvar_trg=file.path(dir_flash_models, datavar, trg, "clustering")
#dirfm_datvar_tsp=file.path(dirfm_datvar_trg, sp_status, this_species)   ## species directory (only declared here)
if(!dir.exists(dirfm_datvar_trg)){dir.create(dirfm_datvar_trg, recursive = T)}

dirff_datvar_trg=file.path(dir_flash_figures, datavar, trg, "clustering")
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

species_source_plot = function(raster_stack, labname){
  ggplot() + theme_bw() + scale_fill_whitebox_c(palette = "deep") + labs(fill = labname) + 
    theme(legend.position = "top", legend.title = element_text(size=18, vjust=0.7), legend.text = element_text(size=17),
          legend.key.height= unit(1, 'cm'), legend.key.width= unit(3, 'cm')) +
    geom_spatvector(data = wcl, fill = "grey85", col="grey90") +
    geom_spatraster(data = raster_stack)
}

## this function helps to save pdf and png plots (from an object) at the same time
plot_pdf_png = function(the_plot, plot_dir, plot_name, pointsize=12, pdf_width=7, pdf_height = 7){
  fname=paste0(plot_name,".pdf")
  pdf(file=file.path(plot_dir, fname), pointsize=pointsize, width = pdf_width, height = pdf_height)
  print(the_plot)
  dev.off()
  
  fname=paste0(plot_name,".png")
  png(file=file.path(plot_dir, fname), pointsize=pointsize, width = round(pdf_width*68.6), height = round(pdf_height*68.6))
  print(the_plot)
  dev.off()
}

##### main
print("read data")

trg_div = readRDS(file = file.path(dirbd_env_trg, "Subregions.rds")) ## sf class
trg_silh = trg_div %>% vect()  ## terra class, for plotting

##### P2: suitability in Japan, arrange species by hierarchical grouping

### how many of those species could establish populations in Japan?
fname = paste0("Species-prefecture_MF-",model_clean_criteria,"_mask-",mask,".csv")
suit_pref = read.csv(file.path(dirbm_datvar_trg, model_clean_criteria, sp_status, fname), row.names = 1) %>% 
  mutate(prefecture=subregion) %>% select(-subregion)
# number of species
n_pot_Jp = suit_pref$species %>% unique() %>% length()
print(paste0("potential: establishing in prefectures ", n_pot_Jp))

##### define orders that we have defined before:
# latitudinal order of prefectures
pref_latup = trg_div %>% distinct(prefecture, latitude) %>% arrange(latitude)
#pref_latdown = trg_div %>% distinct(prefecture, latitude) %>% arrange(-latitude)

# #suit_block = read.csv(file.path(dirbm_datvar_trg, "fc2", "potential", "Species-block_MF-fc2_mask-both.csv"), row.names = 1)
# bm_transfer_rast=rast(file.path(dirbm_datvar_trg, "fc2", "potential", "richness_MF-fc2_mask-both.tif"))
# minmax(bm_transfer_rast)  # min and max spp per cell
# 
## reference to decide on how to preform the hierarchical clustering
## https://www.statology.org/hierarchical-clustering-in-r/
## https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/

# create all species - prefecture combinations and fill in with zeros where NA
# because hierarchical clustering analysis does not accept NA
# we are assuming that those NA are places where species could not establish a 
# population, thus, 0 suitability prefectures:
exo_names = suit_pref %>% distinct(species)
sp_pf_combo = expand_grid(species = exo_names$species, prefecture = pref_latup$prefecture) 
sp_pf_mean = suit_pref %>% distinct(species, prefecture, mean) 
sp_pf_suit = sp_pf_combo %>% full_join(sp_pf_mean) %>% mutate(mean = ifelse(is.na(mean), 0, mean))
# 

## some references about dendrograms
## https://stackoverflow.com/questions/8045538/labelling-ggdendro-leaves-in-multiple-colors
## https://jcoliver.github.io/learn-r/008-ggplot-dendrograms-and-heatmaps.html

## matrices for hierarchical clustering analysis
suit_pref_matrix = sp_pf_suit %>% select(species, prefecture, mean) %>% 
  pivot_wider(names_from = prefecture, values_from = mean) %>% as.data.frame()
rownames(suit_pref_matrix) = suit_pref_matrix$species
suit_pref_matrix = suit_pref_matrix %>% select(-species) %>% as.matrix()

## what you want to cluster has to be in the rows, and you cluster by the columns
m1 = suit_pref_matrix    # species in rows
m2 = t(suit_pref_matrix) # prefectures in rows

################### determine number of clusters

# silhouette method
opt_clust_silh=fviz_nbclust(m2, kmeans, method = "silhouette") + labs(subtitle = "Silhouette method: kmeans")
opt_clust_silh
# # Gap statistic
# set.seed(49)
# opt_clust_gap=fviz_nbclust(m2, kmeans, nstart = 25, method = "gap_stat", nboot = 100) + labs(subtitle = "Gap statistic method: kmeans")

## get the number of clusters that makes sense: 
## the one that maximized the silhouette metrics
clust_silh_data = opt_clust_silh$data
n_clust_d = clust_silh_data %>% filter( y == max(clust_silh_data$y)) #%>% select(clusters)
n_clust=n_clust_d$clusters %>% as.numeric()
n_clust

km_res <- kmeans(m2, centers = n_clust, nstart = 20)
sil <- silhouette(km_res$cluster, dist(m2))
fviz_silhouette(sil)
d12=fviz_cluster(km_res, m2, ellipse.type = "norm", main="")
d23=fviz_cluster(km_res, m2, ellipse.type = "norm", main="" , axes=c(2,3))

d123=ggarrange(d12, d23,  nrow=1, ncol=2, common.legend = T)
d123

clust_d123=ggarrange(opt_clust_silh, d123, nrow=1, ncol=2, widths = c(1,2))
plot_pdf_png(clust_d123, file.path(dirff_datvar_trg), "HowManyClusters_pref", pointsize = 4, pdf_width=21, pdf_height = 7)
clust_d123

#################################################
## some colours for the clusters
clust_col = tibble(cluster = paste0("cluster_", 1:n_clust),    # prefecture clusters
                   cl_colour = scales::hue_pal(l=75)(n_clust))
species_col = tibble(group = paste0("group_", 1:n_clust),      # species groups
                     gr_colour = scales::hue_pal(l=35)(n_clust))

##### create the clusters for the analysis

## first we cluster the prefectures, according to their exotic species communities ; 
## after this we will cluster the species to be able to arrange them in the plots by groups and risk to the country
pref_dist = dist(m2, method="euclidean")    ## distance matrix
pref_clust = hclust(pref_dist, method = "ward.D2")   # I think this is the best grouping; I tried all, but this one is the one that looks more reasonable.
pref_dendro_str = as.dendrogram(pref_clust)
## this dendrogram has the leaves in very unituitive orders, 
## so we rotate the leaves so they are organised by latitude as much as possible
pref_dendro_str_lat = pref_dendro_str %>% rotate(pref_latup$prefecture) #%>% as.dendrogram()  ## reorder the dendogram by prefecture
#pref_dendro_str_lat %>% plot()

#####
pref_cluster = cutree(pref_dendro_str_lat, n_clust) # cutree to create prefecture clusters
## organise in tibble and join to colour code
pref_cluster_df = tibble(pref_cluster) %>% mutate(prefecture = names(pref_cluster)) %>% 
  mutate(cluster = paste0("cluster_", pref_cluster)) %>% full_join(clust_col, by = join_by(cluster)) 
pref_dendro_data = dendro_data(pref_dendro_str_lat)  # get dendrogram as a list 
## vector with the "final" order of prefectures: mix between clustering and latitude 
pref_order = label(pref_dendro_data)$label 
## join to the tibble above and make sure prefectures are a factor
pref_labs = label(pref_dendro_data) %>% mutate(prefecture = label) %>%  
  full_join(pref_cluster_df) %>% mutate(prefecture = factor(prefecture, levels=pref_order))
pref_dendro_data$labels = pref_labs  ## reassign info to dendrogram
saveRDS(pref_dendro_data, file.path(dirfm_datvar_trg, "pref_dendro_data.rds"))
# pref_dendro_data = readRDS(file.path(dirbm_datvar_trg, "pref_dendro_data.rds"))

## plot the dendrogram
pref_dendro_ggdendrogram =
  ggdendrogram(data = pref_dendro_data, rotate = TRUE) +
  theme(axis.text.y = element_text(size = 7, colour = pref_dendro_data$labels$cl_colour))
pref_dendro_ggdendrogram
plot_pdf_png(pref_dendro_ggdendrogram, file.path(dirff_datvar_trg), "Pref_dendro", pointsize = 4, pdf_width=7/2, pdf_height = 7)

# join prefecture clusters to the political division and visualize
trg_div2 = trg_div %>% full_join(pref_cluster_df, by = join_by(prefecture)) 
pref_block_plot = ggplot() + theme_bw() +  theme(legend.position = "top") +
  geom_spatvector(data = trg_div2, mapping = aes(fill = cluster)) + 
  scale_fill_manual(values = clust_col$cl_colour)
pref_block_plot 

#### 

div_level = "cluster"   ## define division level name
newgroup = c("species", eval(div_level))

## make the dataframe for the tile plot:
srclvl_meanOcc = suit_pref %>%
  full_join(trg_div2) %>% group_by(across(all_of(newgroup))) %>%
  summarise(mean = mean(mean)) %>% mutate(subregion = get(eval(div_level)))

suit_prefClust = countsForPlot_SubregSpecies(srclvl_meanOcc) # %>% na.omit()

##### here: trying to order species by "risk": number or clusters, number of prefectures and average country suitability
sp1 = sp_pf_mean %>% group_by(species) %>% summarise(spp_npref=n()) %>% arrange(-spp_npref)
sp2 = sp_pf_suit %>% group_by(species) %>% summarise(spp_Mean=mean(mean)) %>% arrange(-spp_Mean)
sp_measures = suit_prefClust %>% distinct(species, spp_nReg) %>% rename(spp_nclust = spp_nReg) %>% 
  full_join(sp1) %>% full_join(sp2) %>%  mutate(clust_pref_mean = (spp_nclust*spp_npref*spp_Mean) ) %>% arrange(-clust_pref_mean)

## now cluster species according to the prefectures they can inhabit
spp_dist = dist(m1,  method="euclidean")    ## distance matrix
spp_clust = hclust(spp_dist, method = "ward.D2")   # using ward.D2 inherited from prefectures
spp_dendro_str = as.dendrogram(spp_clust)
## dendrogram order might be a bit weird, 
## so we reorganise leaves by "risk" ~ clust_pref_mean
spp_dendro_str_lat = spp_dendro_str %>% rotate(sp_measures$species) #%>% as.dendrogram()  ## reorder the dendogram by prefecture
spp_group = cutree(spp_dendro_str_lat, n_clust)  ## create species groups
spp_group_df = tibble(spp_group) %>% mutate(species = names(spp_group)) %>% ## join to colours
  mutate(group = paste0("group_", spp_group)) %>% full_join(species_col, by = join_by(group))
spp_dendro_data = dendro_data(spp_dendro_str_lat)

spp_order = label(spp_dendro_data)$label  ## the order mix between grouping and risk
## as for prefectures, join tibble and assure species is a factor
spp_labs = label(spp_dendro_data) %>% mutate(species = label) %>% 
  full_join(spp_group_df) %>% mutate(species = factor(species, levels=spp_order))
spp_dendro_data$labels = spp_labs    ## reassign info to dendrogram
saveRDS(spp_dendro_data, file.path(dirfm_datvar_trg, "spp_dendro_data.rds"))

## plot dendrogram:
spp_dendro_ggdendrogram =
  ggdendrogram(data = spp_dendro_data) +
  # scale_y_continuous(expand = c(0, 0.5), labels = levels(spp_dendro_data$labels$species),
  #                    breaks = 1:length(spp_dendro_data$labels)) +
  theme(axis.text.x = element_text(size = 7, colour = spp_dendro_data$labels$gr_colour))
spp_dendro_ggdendrogram
plot_pdf_png(spp_dendro_ggdendrogram, file.path(dirff_datvar_trg), "Spp_dendro", pointsize = 4, pdf_width=7/2, pdf_height = 7)

### now we make sure that both prefecture and cluster heatmaps are in the same order:
suit_prefClust = suit_prefClust %>% full_join(spp_labs, by = join_by(species)) %>%
  mutate(species = factor(species, levels = spp_labs$species)) #%>% arrange(species)

suit_pref_plotOrd = suit_pref %>% 
  full_join(pref_labs, by = join_by(prefecture)) %>% 
  mutate(prefecture = factor(prefecture, levels=pref_labs$prefecture)) %>% 
  full_join(spp_labs, by = join_by(species)) %>%
  mutate(species = factor(species, levels = spp_labs$species))

gSS_pref =
  ggplot(data=suit_pref_plotOrd, aes(y = prefecture, x = species)) +
  theme() + ylab(NULL) + xlab(NULL) +
  theme(legend.position = "none", legend.text = element_text(size=9),
        legend.key.height= unit(0.8, 'cm'), legend.key.width= unit(1.7, 'cm'),
        axis.text.y = #element_blank(),
          element_text(size = 6.2, colour = pref_labs$cl_colour),
        axis.text.x = element_blank()
        #  element_text(size=7, angle = 90, hjust=1.05, vjust=0.5, colour = spp_labs$gr_colour)) + #
  ) +
  geom_tile(aes(fill=mean)) + scale_fill_whitebox_c("mean suitability", palette = "deep") +
  scale_y_discrete(position = "right") #+ scale_x_discrete(position = "top")
gSS_pref

gSS_cluster =
  ggplot(data=suit_prefClust, aes(y = cluster, x = species)) +
  theme() + ylab(NULL) + xlab("species") +
  theme(legend.position = "bottom", legend.text = element_text(size=9),
        legend.key.height= unit(0.8, 'cm'), legend.key.width= unit(1.7, 'cm'),
        axis.text.y = element_text(size = 11, colour = clust_col$cl_colour),
        axis.text.x = #element_blank()
          element_text(size=7, angle = 90, hjust=1.05, vjust=0.5, colour = spp_labs$gr_colour)  #
  ) +
  geom_tile(aes(fill=mean)) + scale_fill_whitebox_c("mean suitability", palette = "deep") +
  scale_y_discrete(position = "right") #+ scale_x_discrete(position = "top")
gSS_cluster
leg=get_legend(gSS_cluster)
heatmaps = ggarrange(gSS_pref, gSS_cluster, nrow=2, legend.grob = leg, heights = c(2,0.7))
heatmaps 
plot_pdf_png(heatmaps , file.path(dirff_datvar_trg), paste0("Heatmaps_PrefClusters"), pointsize = 4, pdf_width=7, pdf_height=(7*1.5))


pref_cluster_df = pref_cluster_df %>% mutate(prefecture = factor(prefecture, levels=pref_order)) %>% arrange(prefecture)
write.csv(pref_cluster_df, file.path(dirfm_datvar_trg, "PrefectureClusterColor.csv"))

spp_group_df = spp_group_df %>% mutate(species = factor(species, levels=spp_order)) %>% arrange(species)
write.csv(spp_group_df, file.path(dirfm_datvar_trg, "SpeciesGroupColor.csv"))

######################################
############
print("owaru")
print("bye bye")


##### code end  




