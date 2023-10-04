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
library(raster)
library(sf)
library(terra) 
library(fasterize)
library(dplyr)
library(ggplot2)
library(tidyterra)

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

# ## get parameters
# args <- commandArgs(trailingOnly = TRUE)   ## allows to pass parameters to the script
# # test if there is at least one argument: if not, return an error
# if (length(args)<1) {
#   stop("Not enough parameters", call.=FALSE)
# }
# 
# trg=args[1]


## reading directories:
dirbd_env=file.path(dir_bucket_data, "EnvironmentalRasters")
dirbd_env_trg=file.path(dirbd_env, trg)  ## subregions' shapes

# writing directories:
dirfd_env_trg=file.path(dir_flash_data, "EnvironmentalRasters", trg)  ## subregions' shapes
if(!dir.exists(dirfd_env_trg)){dir.create(dirfd_env_trg, recursive = T)}

## get the data I need from natural earth

target_region_div=ne_states(country=trg, returnclass = "sf")
# target_region_div is an SpatialPolygonsDataFrame describing Japan at prefecture level
#plot(target_region_div)

## then convert target_region_div object to sf so we can operate on it more easily
#target_region_div_sf=st_as_sf(target_region_div)

## select only what I need
## name_alt and name_local will be to add Izu and Ogasawara Islands later
## wikipedia is added only because it gives a blank silhouette
target_region_div_data = target_region_div %>% mutate(name=woe_name) %>% 
  select(name, name_alt, name_local, region, latitude, longitude, type, type_en, provnum_ne, labelrank, wikipedia, geometry)

## name is going to be used to reorganize the prefectures as we need
target_region_div_data = target_region_div_data %>% mutate(prefecture_official = name) ## prefectures of Japan

## now find the data and organise it as I need it:
target_region_div_data$region   ## larger regions of Japan (subregion)
## Saga and Nagasaki prefectures don't have region name, they belong to Kyushu, asign it
target_region_div_data$region[target_region_div_data$name %in% c("Saga", "Nagasaki")]="Kyushu"

## Okinawa officially belongs to Kyushu, but this data has it separate, 
## so I add two more fields to register that without losing current structure:
target_region_div_data$region_official=target_region_div_data$region
target_region_div_data$region_official[target_region_div_data$name %in% c("Okinawa")] = "Kyushu"

## name_alt is not an occupied field, so I will use it to split find small islands (later)
target_region_div_data$name_local = target_region_div_data$prefecture

## rearrange a bit the columns
target_region_div_data = target_region_div_data %>% 
  select(region, name, region_official, prefecture_official, longitude, latitude, everything())

## split Izu islands and Ogasawara from Tokyo prefecture
# first, select Tokyo from the dataset
TokyoPrefecture=target_region_div_data[target_region_div_data$prefecture_official=="Tokyo",]
## we can split (disaggregate) the prefecture into its individual polygons (kind of each island)
## but first we need to put it back as a SpatialPolygonsDataFrame:
tokyo=as(TokyoPrefecture, "Spatial")
tokyo_islands=disaggregate(tokyo)
## they we make it a dataframe again for easy manipulation
TokyoPrefectureIslands=st_as_sf(tokyo_islands)

### change info that will not make much sense:
TokyoPrefectureIslands$name_alt = NA     ## this will be the archipelago name: Izu, Ogasawara
TokyoPrefectureIslands$name_local = NA   ## this will be the island name: whatever they are called...
#TokyoPrefectureIslands$longitude = NA  
#TokyoPrefectureIslands$latitude = NA  

TokyoPrefectureIslands$name_alt[1] = "Tokyo-Mainland"
TokyoPrefectureIslands$name_local[1] = "Tokyo-toshobu"
TokyoPrefectureIslands$type[1] = "Toshobu"
TokyoPrefectureIslands$type_en[1] = "Metropolis"

TokyoPrefectureIslands$type[2:nrow(TokyoPrefectureIslands)] = "Shima"
TokyoPrefectureIslands$type_en[2:nrow(TokyoPrefectureIslands)] = "Island"
TokyoPrefectureIslands$labelrank[2:nrow(TokyoPrefectureIslands)] = 2  ## same as Okinawa

## index k for each subsubregion (island)
### I will go one at the time and figure out their names and archipelago or something like that 

# k=23
# TokyoPrefectureIslands[k,] %>% extent()           ## to know where it is, look in google maps
# TokyoPrefectureIslands[c(1,2,5,k),"wikipedia"] %>% plot()  ## to know how far from metropoli is.
# TokyoPrefectureIslands[k,"wikipedia"] %>% plot()  ## to know the shape of the island

# TokyoPrefectureIslands[c(2,3,4,9,10,12,13,14,15,23),"wikipedia"] %>% plot()  ## try to plot Owasawara together
# TokyoPrefectureIslands[c(1,5,6,7,8,16,17,18,19,20,21),"wikipedia"] %>% plot()  ## try to plot Izu together

TokyoPrefectureIslands[23,]$name_local="Nakano/Sasayo" ##? 
TokyoPrefectureIslands[22,]$name_local="Minami-Tori" ##again?
TokyoPrefectureIslands[22:23,]$name_alt="Ogasawara"

TokyoPrefectureIslands[21,]$name_local="Mikura"   ## still Izu
TokyoPrefectureIslands[20,]$name_local="Unnamed-rock"
TokyoPrefectureIslands[19,]$name_local="Aogashima"
TokyoPrefectureIslands[18,]$name_local="MiddleOfNowhere2"       ## Beyonesu? middle of nowhere literally, but still Aogashima archi
TokyoPrefectureIslands[17,]$name_local="MiddleOfNowhere1"       ## Sumisu? cannot find it on google maps, but still Aogashima archi
TokyoPrefectureIslands[16,]$name_local="Tori"   ## Izu
TokyoPrefectureIslands[16:21,]$name_alt="Izu"

TokyoPrefectureIslands[15,]$name_local="Kitano"  ## is the closest in google maps, but the shape is quite different
TokyoPrefectureIslands[14,]$name_local="Nakodo"
TokyoPrefectureIslands[13,]$name_local="Ototo"
TokyoPrefectureIslands[12,]$name_local="Kita-Iwo"
TokyoPrefectureIslands[11,]$name_local="Minami-Tori-2"   ## one of the farthest island in Owasawara, but shape is different
TokyoPrefectureIslands[10,]$name_local="Minami-Iwo"
TokyoPrefectureIslands[9,]$name_local="Nishino"    # Ogasawara
TokyoPrefectureIslands[9:15,]$name_alt="Ogasawara"

TokyoPrefectureIslands[8,]$name_local="Ooshima"
TokyoPrefectureIslands[7,]$name_local="Nii"   # Izu
TokyoPrefectureIslands[6,]$name_local="Miyake"   # Izu
TokyoPrefectureIslands[5,]$name_local="Hachijo"   # Izu
TokyoPrefectureIslands[5:8,]$name_alt="Izu"

TokyoPrefectureIslands[4,]$name_local="Chichi"
TokyoPrefectureIslands[3,]$name_local="Haha"
TokyoPrefectureIslands[2,]$name_local="Iwo"
TokyoPrefectureIslands[2:4,]$name_alt="Ogasawara"

#####
##### try to get a good latitude, longitude value for the islands

## Ogasawara
geo=TokyoPrefectureIslands[TokyoPrefectureIslands$name_alt=="Ogasawara",]$geometry
coords=NULL
for(k in 1:length(geo)){
  coords=rbind(coords, geo[[k]][[1]])
}
colnames(coords)=c("longitude", "latitude")
coords_mean=colMeans(coords) 

TokyoPrefectureIslands[TokyoPrefectureIslands$name_alt=="Ogasawara",]$longitude=coords_mean["longitude"]
TokyoPrefectureIslands[TokyoPrefectureIslands$name_alt=="Ogasawara",]$latitude=coords_mean["latitude"]
TokyoPrefectureIslands[TokyoPrefectureIslands$name_alt=="Ogasawara",]$region = "Ogasawara"

## Izu
geo=TokyoPrefectureIslands[TokyoPrefectureIslands$name_alt=="Izu",]$geometry
coords=NULL
for(k in 1:length(geo)){
  coords=rbind(coords, geo[[k]][[1]])
}
colnames(coords)=c("longitude", "latitude")
coords_mean=colMeans(coords) 

TokyoPrefectureIslands[TokyoPrefectureIslands$name_alt=="Izu",]$longitude=coords_mean["longitude"]
TokyoPrefectureIslands[TokyoPrefectureIslands$name_alt=="Izu",]$latitude=coords_mean["latitude"]
TokyoPrefectureIslands[TokyoPrefectureIslands$name_alt=="Izu",]$region = "Izu"

# ## these are the two archipelagos I need to aggregate:
# TokyoPrefectureIslands[c(2,3,4,9,10,12,13,14,15,23),"wikipedia"] %>% plot()  ## try to plot Owasawara together
# TokyoPrefectureIslands[c(1,5,6,7,8,16,17,18,19,20,21),"wikipedia"] %>% plot()  ## try to plot Izu together

## this groups/"dissolves" all islands with the same name_alt into a single multipolygon
Tokyo3sects=TokyoPrefectureIslands %>% 
  group_by(region_official, region, prefecture_official, name_alt, latitude, longitude, 
           type, type_en, provnum_ne, labelrank, wikipedia) %>% summarise()

Tokyo3sects$name=Tokyo3sects$name_alt  ## put the "name-alt" of the region back to name
Tokyo3sects$name_local=Tokyo3sects$name_alt

## Tokyo is ready now
## so now we remove the Tokyo record from the original data set
NoTokyo=target_region_div_data[target_region_div_data$prefecture_official!="Tokyo",]
#NoTokyo[NoTokyo$region == "Kanto",]

## now I can bind it as before:
Japan_StudySubregions=rbind(NoTokyo,Tokyo3sects)
Japan_StudySubregions=Japan_StudySubregions %>% mutate(prefecture = name) 

### I will create two more "admin" groups:
# 1) island, which will be literal the natural island groups 
# 2) block, something in between island and region to group few regions together 
##          mainly to avoid having many similarly looking plots

Japan_StudySubregions=Japan_StudySubregions %>% 
  mutate(island = region) %>% mutate(block = region) 

Honshu=c("Tohoku", "Kanto", "Chubu", "Chugoku", "Kinki")
#Honshu_N = c("Tohoku")
Honshu_C = c("Kanto", "Chubu") 
Honshu_S = c("Chugoku", "Kinki")
Shi_Kyu = c("Shikoku", "Kyushu")

Japan_StudySubregions=Japan_StudySubregions %>% 
  mutate(island=ifelse(region %in% Honshu, "Honshu", region)) 

#Japan_StudySubregions[Japan_StudySubregions$region %in% Honshu_N,]$block = "Honshu_N"
Japan_StudySubregions[Japan_StudySubregions$region %in% Honshu_C,]$block = "Honshu_C"
Japan_StudySubregions[Japan_StudySubregions$region %in% Honshu_S,]$block = "Honshu_S"
Japan_StudySubregions[Japan_StudySubregions$region %in% Shi_Kyu,]$block = "Shi_Kyu"

## parallel: A40; A35B40; A30B35; Arc_Izu; Arc_Ogasawara, Arc_Okinawa
Japan_StudySubregions=Japan_StudySubregions %>% 
  mutate(parallel = ifelse(block %in% c("Okinawa", "Izu", "Ogasawara"), "Archipelagoes", 
                           ifelse(block %in% c("Honshu_S", "Shi_Kyu"), "HonS_Shi_Kyu",
                                  ifelse(block %in% c("Honshu_C", "Tohoku"), "Honshu_CN", block))))

Japan_StudySubregions$region %>% unique() %>% length()  ## no. of regions
Japan_StudySubregions$block %>% unique() %>% length()  ## no. of blocks
Japan_StudySubregions$island %>% unique() %>% length()  ## no. of islands
Japan_StudySubregions$parallel %>% unique() %>% length()  ## no. of latitudinal groups

## and make factors by latitude
Japan_StudySubregions=Japan_StudySubregions %>% arrange(latitude, provnum_ne, name)
pref_names=Japan_StudySubregions$name %>% unique 
reg_names=Japan_StudySubregions$region %>% unique 
block_names=Japan_StudySubregions$block %>% unique 
paral_names=Japan_StudySubregions$parallel %>% unique 

Japan_StudySubregions=Japan_StudySubregions %>% 
  mutate(name = factor(Japan_StudySubregions$name, levels=pref_names)) %>% 
  ## add a new field called prefecture with "our" prefecture needs
  mutate(prefecture = factor(Japan_StudySubregions$name, levels=pref_names)) %>% 
  mutate(region = factor(Japan_StudySubregions$region, levels=reg_names)) %>% 
  mutate(block = factor(Japan_StudySubregions$block, levels=block_names)) %>% 
  mutate(parallel = factor(Japan_StudySubregions$parallel, levels=paral_names)) %>% 
  select(-name_alt, -wikipedia) %>% select(name, prefecture, region, block, parallel, island, everything())

## check product:
jp = vect(Japan_StudySubregions)

ggplot(data = jp) + theme_bw() +
  geom_spatvector(aes(fill=region)) + 
  theme(legend.position = "top", legend.title = element_text(size=9), legend.text = element_text(size=7)) 

ggplot(data = jp) + theme_bw() +
  geom_spatvector(aes(fill=block)) + 
  theme(legend.position = "top", legend.title = element_text(size=9), legend.text = element_text(size=7)) 

ggplot(data = jp) + theme_bw() +
  geom_spatvector(aes(fill=parallel)) + 
  theme(legend.position = "top", legend.title = element_text(size=9), legend.text = element_text(size=7)) 

## save product
filename=paste0("Subregions",".rds")
saveRDS(Japan_StudySubregions, file = file.path(dirfd_env_trg, filename))

print("owaru, mata ne~")

#############################################

### ## code ends


