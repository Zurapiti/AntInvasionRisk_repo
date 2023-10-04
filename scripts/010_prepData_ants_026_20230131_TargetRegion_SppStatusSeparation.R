## 
## author:  Yazmin Zurapiti
## title:   ants_026_20230131_TargetRegion_SppStatusSeparation.R
## aim:     make separate files to list for each region: native, invader, potential
## version: 0.4 2022.04.27: use snapped points and update directories.
##                          save this as data instead of models
##              2022.08.01: minimal adjusting to the new files.  
##              2023.01.31: now accepting save_date as parameter
##
## inputs: 
##   1.inputdir: 
##   1.inputfile:
## outputs: 
##   1.outputdir: 
##   1.outputfile: 
## parameters: receives 2: target_area_name (e.g. Japan), area_admin_rank (in country, bentity2)
## call: ants_026_20230131_TargetRegion_SppStatusSeparation.R [save_date] [target_area_name] [area_admin_rank] 


## load libraries
# library()
library(dplyr)
library(tidyr)
library(purrr)

## clear memory
rm(list = ls())

local_analysis=F
if(local_analysis){script06path="./06_defineDirectories.R"} else {
  script06path="./AntInvasionRisk/scripts/06_defineDirectories.R"}
source(script06path)

args <- commandArgs(trailingOnly = TRUE)   ## allows to pass parameters to the script
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("Not enough parameters", call.=FALSE)
}
## testing params
# save_date="2023-01-16"
# target_area_name = "Madagascar"
# area_admin_rank = "country"

## passed as args
save_date=args[1]
target_area_name = args[2]
area_admin_rank = args[3]


### ## main

## saving directories:
## if it does not exist yet, make a new directory for this analysis
dir_sepStatus=file.path(dir_flash_data, "EnvironmentalRasters")
if(!file_test("-d", dir_sepStatus)){ dir.create(dir_sepStatus, recursive=TRUE) }

## read directories
### post parameter directories
newdata=paste0("invasionRiskData_5km_", save_date)
dirb_newdata=file.path(dir_bucket_data, newdata, "model_inputs")   ## read

print("reading data")

## read cleaned species
filename=paste0("ExoticNames","_Clean.csv")
spp_names=read.table(file.path(dirb_newdata, filename))
names(spp_names)="species"

filename="exotics_gabidb_list_thin.rds"
spp_data=readRDS(file.path(dirb_newdata, filename))

print("processing data")

spp_data_tibb=spp_data %>% #spp_data[names(spp_data) %in% spp_names$species] %>% 
  map_dfr( ~as_tibble(.), .id="species")  %>% 
  rename(status = exotic) %>% mutate(status=ifelse(is.na(status), "native", "exotic"))
print("number of species in exotics_gabidb_list_thin.rds")
spp_data_tibb %>% distinct(species) %>% count()    ## ok 245; with the new list including Benoit's considerations


## we need the species, area_admin_rank and status info
# this !! rlang::parse_expr( par ) allows to evaluate the value of par and use that as the value passed to filter
spp_region_status = spp_data_tibb %>% distinct(species, !! rlang::parse_expr(area_admin_rank), status) 

print("number of species")
spp_region_status %>% distinct(species) %>% nrow()   ## number of species 245
print("number of admin regions")
spp_region_status %>% distinct(!! rlang::parse_expr(area_admin_rank)) %>% nrow()   # number of admin regions 244

#j=1
#for(j in 1:length(regions_goodNames)){
#  trg=regions_goodNames[j]

trg=target_area_name
target_area_name=gsub("-", " ", trg)
print(target_area_name)

print("saving data")

# make a directory for each area
dir_trg=file.path(dir_sepStatus, trg)
if(!file_test("-d", dir_trg)){ dir.create(dir_trg, recursive=TRUE) }

### then, we need to separate the (species/records?) into three categories:
# native to target_area - we actually don't care about these ones
# exotic to target_area - without records -[potential]- can be stacked/modeled as they are
#                         with records -[invaders]-  need some special modeling or extra attention

# native = native to target_area_name, should not be included in any further analysis for that region
native = spp_region_status %>% filter(!! rlang::parse_expr(area_admin_rank) == target_area_name & status == "native") %>% select(species)

## 20220725: remove native species from the species set
rest = spp_region_status %>% filter(!(species %in% native$species))
## to remove species like Lasius sakagamii which is native to North of Japan, but invasive to Okinawa

# invader = exotic and established in target_area_name, need extra processing and are the "test cases"
#           this will show us how good our predictions could be
invader = rest %>% filter(!! rlang::parse_expr(area_admin_rank) == target_area_name & status == "exotic") %>% select(species)
## 20220725: remove invader species from the species set
rest = rest %>% filter(!(species %in% invader$species))

# potential = would be exotic to the region but haven't reached yet
#             these are the ones we are interested in assessing the risk for
## these are the potential new invaders
potential = rest %>% distinct(species)

## both so we can take the ones that are not in the region
in_region = native %>% full_join(invader, by="species")
## so... interesting that there are some species that are listed as native and invasive in the same country.
# maybe it is something decided at a lowe scale resolution. 
nrow(spp_names) == nrow(in_region) + nrow(potential)


## now need to store in file native, invader and potential
filename=paste0(deparse(substitute(native)), ".csv")
write.table(x=native, file=file.path(dir_trg, filename), col.names = F,  row.names = F, sep=",")

filename=paste0(deparse(substitute(invader)), ".csv")
write.table(x=invader, file=file.path(dir_trg, filename), col.names = F,  row.names = F, sep=",")

filename=paste0(deparse(substitute(potential)), ".csv")
write.table(x=potential, file=file.path(dir_trg, filename), col.names = F,  row.names = F, sep=",")

print(paste0("number of ", deparse(substitute(native)), " ", nrow(native)))
print(paste0("number of ", deparse(substitute(invader)), " ", nrow(invader)))
print(paste0("number of ", deparse(substitute(potential)), " ", nrow(potential)))
#} # end of for

## code ends




