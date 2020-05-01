# Species map for Phillipine bats, primates, and rodents
# Author: Anna Willoughby
# Date: 24 April 2020

# setup library
library(rredlist) # for accessing IUCN data
## data/code cleaning 
library(dplyr) 
library(tidyr)
library(magrittr)
## mapping packages
library(raster)
library(rgeos)
library(fasterize)
library(sf)
library(sp)
library(maptools)
# color library
library(viridis)

## IUCN API token <ideal to make your own at https://apiv3.iucnredlist.org/api/v3/token 
token <- "303c8cdd6b8792e3bf4e740719940a81310f0a04297faf9d10fd1d5691ca9484"

# ISO code for the Phillipines is PH (2) and PHL (3)
# search IUCN for species list for PHL using rredlist. this outputs species present in PH
PH_species <- rl_sp_country("PH", key = token)

# loop species through search to retrieve upper taxonomy
taxon_info <- rl_search(PH_species$result$scientific_name[[1]], key = token)
df <- taxon_info$result
for(i in PH_species$result$scientific_name){
  taxon_info_loop <- rl_search(i, key = token)
  df_loop <- taxon_info_loop$result
  df <- rbind(df, df_loop)
}
# get rid of duplicate first row 
df <- unique(df) 

# filter to bats, primates, and rodents 
PH_target_taxa <- dplyr::filter(df, order %in% c("RODENTIA", "PRIMATES", "CHIROPTERA"))
# remove humans 
PH_target_taxa <- dplyr::filter(PH_target_taxa, scientific_name != "Homo sapiens")

# load in spatial files 
terr = shapefile("data/TERRESTRIAL_MAMMALS_20-04-2020/TERRESTRIAL_MAMMALS.shp", verbose = T)
# check that all species have a map; all species have a map 
PH_target_taxa$IUCN_map <- ifelse(PH_target_taxa$scientific_name %in% terr@data$binomial, "TRUE", "No")
write.csv(PH_target_taxa, file = "data/PH_target-species.csv")

# load in PH ADM, downloaded from gadm 
PH_gadm <- readRDS("data/PhillipinesADM/gadm36_PHL_2_sp.rds")
PH_gadm_proj <- spTransform(PH_gadm, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # make sure projection is to WGS84
PH_bbox <- c(117.17427453, 5.58100332277, 126.537423944, 18.5052273625) # for croping the output window

# filter shapefiles to species of interest
PH_shapes <- subset(terr, terr@data$binomial %in% PH_target_taxa$scientific_name)


# create Phillipines species richness map (not distinguishing orders)
PH_sf <- st_as_sf(PH_shapes) # change file from sp to sf for fasterize function
r <- raster(ncol = 1001, nrow = 1001)
extent(r)
PH_rs <- fasterize(PH_sf, r, field = "origin", fun = "sum")
if (require(rgdal)) {
  PH_raster <- writeRaster(PH_rs, filename="PH_species.tif", format="GTiff", overwrite=TRUE)
}

# clip the raster to only Phillipines
PH_sp_map <- mask(PH_rs, PH_gadm) # changes raster coverage
PH_sp_map_crop <- crop(PH_sp_map, PH_gadm) # changes raster extent 

# plot raster + PH adm boundaries 
plot(PH_sp_map_crop)
plot(PH_gadm, add = TRUE)

# extract species richness values per PH ADM 3 vector objects
species_per_ADM2 <- extract(x = PH_sp_map_crop, y = PH_gadm, fun = max, df = TRUE) # this values for each pixel
colnames(species_per_ADM2)[2] <-"target_sp_richness"

# filter shapefiles for each taxa of interest 
# Chiroptera species richness
PH_bats <- subset(PH_shapes, PH_shapes@data$order_  == "CHIROPTERA")
PH_bats_sf <- st_as_sf(PH_bats) # change file from sp to sf for fasterize function
PH_bats_rs <- fasterize(PH_bats_sf, r, field = "origin", fun = "sum")
PH_bats_map <- mask(PH_bats_rs, PH_gadm) # clip to Phillipines
PH_bats_map_crop <- crop(PH_bats_map, PH_gadm) # changes raster extent 
# calculate values by ADM2 
bats_per_ADM2 <- extract(x = PH_bats_map_crop, y = PH_gadm, fun = max, df = TRUE) # this values for each pixel
colnames(bats_per_ADM2)[2] <-"bats_sp_richness"

# Chiroptera order map 
PH_bats_order <- union(PH_bats)

# Primates species richness
PH_primates <- subset(PH_shapes, PH_shapes@data$order_  == "PRIMATES")
PH_primates_sf <- st_as_sf(PH_primates) # change file from sp to sf for fasterize function
PH_primates_rs <- fasterize(PH_primates_sf, r, field = "origin", fun = "sum")
PH_primates_map <- mask(PH_primates_rs, PH_gadm) # clip to Phillipines
PH_primates_map_crop <- crop(PH_primates_map, PH_gadm) # changes raster extent 

# calculate values by ADM2 
primates_per_ADM2 <- extract(x = PH_primates_map_crop, y = PH_gadm, fun = max, df = TRUE) # this values for each pixel
colnames(primates_per_ADM2)[2] <-"primate_sp_richness"

# Primate order map
PH_primate_order <- union(PH_primates)
  
# Rodents species richness
PH_rats <- subset(PH_shapes, PH_shapes@data$order_ == "RODENTIA")
PH_rats_sf <- st_as_sf(PH_rats) # change file from sp to sf for fasterize function
PH_rats_rs <- fasterize(PH_rats_sf, r, field = "origin", fun = "sum")
PH_rats_map <- mask(PH_rats_rs, PH_gadm) # clip to Phillipines
PH_rats_map_crop <- crop(PH_rats_map, PH_gadm) # changes raster extent 
# calculate values by ADM2 
rats_per_ADM2 <- extract(x = PH_rats_map_crop, y = PH_gadm, fun = max, df = TRUE) # this values for each pixel
colnames(rats_per_ADM2)[2] <-"rats_sp_richness"

# Rodent order map
PH_rats_order <- union(PH_rats)
  
# merge raster values by pixel 
species_per_ADM2 <- cbind(species_per_ADM2, bats_per_ADM2, primates_per_ADM2, rats_per_ADM2)
species_per_ADM2[,3] <- NULL # remove duplicate ID
species_per_ADM2[,4] <- NULL # remove duplicate ID
species_per_ADM2[,5] <- NULL # remove duplicate ID
species_per_ADM2[is.na(species_per_ADM2)] <- 0 # fill in zeros
species_per_ADM2$species_sum <- species_per_ADM2$primate_sp_richness + species_per_ADM2$bats_sp_richness + species_per_ADM2$rats_sp_richness# fill in zeros
# confirm pixels line up (there are issues here)
species_per_ADM2$pixel_aligned <- species_per_ADM2$target_sp_richness - species_per_ADM2$species_sum

pal <- colorNumeric(c("#440154FF", "#46337EFF" ,"#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF"), values(PH_sp_map),
                    na.color = "transparent")

country_map <- leaflet() %>% addTiles() 

country_map %>%
        addPolygons(data=PH_gadm_proj, weight = 2, fillColor = "grey") %>%
        addRasterImage(PH_sp_map_crop, colors = pal, opacity = 0.8) %>%
        addLegend("bottomright", pal = pal, values = ~PH_sp_map_crop,
            title = "Target Order Species Richness",
            opacity = 1)


png("PH_species.png", width = 800, height = 480)
{levelplot(PH_sp_map, 
      col = viridis(8), 
      breaks = c(0,10,20,30,40,50,60,70)) 
  plot(PH_gadm, 
       add = TRUE)}
dev.off()


