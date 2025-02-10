library(auk)
library(lubridate)
library(sf)
library(gridExtra)
library(tidyverse)
library(raster)
library(terra) ## confirm the GDAL version being used
library(tiff)
library(viridis)
library(exactextractr)
library(PresenceAbsence)
library(verification)
library(ebirdst)
library(fields)
library(tidyverse)
library(dggridR)
library(raster)
library(ncdf4)
library(foreign)
library(fs)
library(tools)
setwd("/Users/liyingnceas/GitHub/gulf_of_mexico_slr_birds/")
wd <-"/Users/liyingnceas/GitHub/gulf_of_mexico_slr_birds/"
#load the study area

# Path to the .gpkg file
shp_path <- "data/gis-data-gcoos.gpkg"

# List all layers in the GeoPackage
gcoos_extent <- shp_path %>%  
  st_read(layer = "gcoos") %>% 
# crop, buffer gcoos by 10 km to provide a little wiggly room
  st_buffer(dist = 10000) %>% 
  st_transform(crs = crs(elev))


# Open the NetCDF file as a SpatRaster (terra object)
elev <- rast("/Users/liyingnceas/GitHub/gulf_of_mexico_slr_birds/data/northern_gulf_1_navd88_2010.nc")
# Assign CRS to the raster
crs(elev) <- "EPSG:4326"

shp_dir <- "data/NOS80K_NGoMshoreline/nos80k_NGoM.shp"
ngom<-vect(shp_dir)
# Plot the first raster (e.g., elev)
plot(elev, main = "Elevation and NGOM")
# Add the second raster (ngom) using the 'add' argument to overlay
plot(ngom, add = TRUE, col = terrain.colors(100))  # Adjust color palette if needed

# List the land cover files
landcover_files <- list.files("data/landcover", pattern = "ccap_landcover_20200311", full.names = TRUE)

# Load the rasters as SpatRaster objects
# only pick the years wanted
landcover_rasters <- lapply(landcover_files[3:5], rast)

# Get the extents of each raster
extents <- lapply(landcover_rasters, ext)

# Calculate the minimum and maximum values for x and y (extent coordinates)
xmin <- max(sapply(extents, function(e) e[1]))  # maximum of xmin
xmax <- min(sapply(extents, function(e) e[2]))  # minimum of xmax
ymin <- max(sapply(extents, function(e) e[3]))  # maximum of ymin
ymax <- min(sapply(extents, function(e) e[4]))  # minimum of ymax

# Create a SpatExtent object from these values
common_extent <- ext(xmin, xmax, ymin, ymax)
# Define bounding box coordinates
bbox_coords <- c(xmin = -98, xmax = -79, ymin = 22, ymax = 32)
# Now, crop each raster to the common extent
cropped_rasters <- lapply(landcover_rasters, crop, common_extent)

# Stack the cropped rasters
landcover_stack <- rast(cropped_rasters)

# Check the stack
print(landcover_stack)



# label layers with year
landcover <- names(landcover_stack) %>% 
  str_extract("[0-9]{4}") %>% 
  paste0("y", .) %>% 
  setNames(landcover_stack, .)

# with a radius equal to 200 or 800 c-cap cells 
neighborhood_radius <- 200 * ceiling(max(res(landcover))) / 2
agg_factor <- round(2 * neighborhood_radius / res(landcover))

r <- raster(landcover) %>% 
  aggregate(agg_factor) 
r <- ngom %>% 
  st_transform(crs = projection(r)) %>% 
  rasterize(r, field = 1) %>% 
  # remove any empty cells at edges
  trim()
r <- writeRaster(r, filename = "data/prediction-surface-gcoos.tif", overwrite = TRUE)
r <- rast("data/prediction-surface-gcoos.tif")

species<-"blkski"
# load ebird data
ebird <- read_csv(paste0("data/", species, "/", species, "_19962006_bbox_zf.csv"))

# buffer to create neighborhood (6 km by 6 km or 24km by 24km) around each checklist point 


ebird_buff <- ebird %>% 
  distinct(year = format(observation_date, "%Y"),
           locality_id, latitude, longitude) %>% 
  # match the ebird year with corresponding year landcover data
  mutate(year_lc = paste0("y", year))  %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  # transform to cdl projection
  st_transform(crs = crs(landcover)) %>% 
  st_buffer(dist = neighborhood_radius) %>% 
  # nest by year to  match to land cover data from the corresponding year
  nest(data = c(year, locality_id, geometry))



# iterate over all years extracting landcover for all checklists in each
lc_extract <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year
  regions <- ebird_buff$data[[which(yr == ebird_buff$year_lc)]]
  # get landcover values within each buffered checklist area
  ee <- exact_extract(landcover[[yr]], regions, progress = FALSE)
  # count the number of each landcover class for each checklist buffer
  ee_count <- map(ee, ~ count(., landcover = value))
  # attach the year and locality id back to the checklists
  ee_summ <- tibble(st_drop_geometry(regions), data = ee_count) %>% 
    unnest(data)
  # bind to results
  lc_extract <- bind_rows(lc_extract, ee_summ)
}
#calculate PLAND: the proportion of each land cover class within the neighborhood.
pland <- lc_extract %>% 
  # calculate proporiton
  group_by(locality_id, year) %>% 
  mutate(pland = n / sum(n)) %>% 
  ungroup() %>% 
  dplyr::select(-n) %>% 
  # remove NAs after tallying so pland is relative to total number of cells
  filter(!is.na(landcover))

ccap_classes <- read.csv(file = "data/ccap_landcover_classes.csv")


pland_ccap <- inner_join(pland, ccap_classes,  by = "landcover")

# tranform to wide format, filling in implicit missing values with 0s%>% 


pland_wide <- pland_ccap %>% 
  select(-landcover) %>%
  group_by(year, locality_id, landclass) %>% 
  summarise(pland = sum (pland)) %>%     # Summarize pland to calculate total per landclass
  ungroup() %>%                             
  pivot_wider(names_from = landclass, 
              values_from = pland, 
              values_fill = list(pland = 0))

# save

write_csv(pland_wide, paste0("data/pland_location_6km_", species, "_meds.csv"))

# ------------------------------------------
# add elevation to the data
# ------------------------------------------
locs <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- ebird_buff$data[[which(yr == ebird_buff$year_lc)]]
  locs_cell <- st_set_geometry(regions, NULL) 
  elev_cell <- exact_extract(elev, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(elevation_median = mean(.$value, na.rm = TRUE),
                     elevation_sd = sd(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_cell, .)
  # bind to results
  locs <- bind_rows(locs, elev_cell)
}
# prediction surface covariates
pland_elev <- inner_join(pland_wide, locs, by = c("locality_id", "year"))
write_csv(pland_elev, "data/pland-elev_location_6km_brnpel_MEDS.csv")


# ------------------------------------------
# add northness to checklists and pred surface
# ------------------------------------------
#repeat the same process for northness covariates
nort <- raster("~/GitHub/gulf_of_mexico_slr_birds/data/northness_1KMmn_GMTEDmd.tif")
nort <-  nort %>% 
  crop(., gcoos_extent) %>% 
  projectRaster(crs = crs(landcover)) 
  



#repeat for prediction surface
nort_locs <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- ebird_buff$data[[which(yr == ebird_buff$year_lc)]]
  locs_nort_cell <- st_set_geometry(regions, NULL)
  nort_cell <- exact_extract(nort, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(northness_median = mean(.$value, na.rm = TRUE),
                     northness_sd = sd(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_nort_cell, .)
  # bind to results
  nort_locs <- bind_rows(nort_locs, nort_cell)
}


# prediction surface covariates

pland_elev_nort <- inner_join(pland_elev, nort_locs,  by = c("locality_id", "year"))
write_csv(pland_elev_nort, paste0("data/pland-elev_nort_location_", species, "_MEDS.csv"))


# ------------------------------------------
# add eastness to checklists 
# ------------------------------------------

#repeat with eastness
east <- raster("~/GitHub/gulf_of_mexico_slr_birds/data/eastness_1KMmn_GMTEDmd.tif")
east <-  east %>% 
  crop(., gcoos_extent) %>% 
  projectRaster(crs = crs(landcover))



#repeat for prediction surface
east_locs <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- ebird_buff$data[[which(yr == ebird_buff$year_lc)]]
  locs_east_cell <- st_set_geometry(regions, NULL) 
  east_cell <- exact_extract(east, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(eastness_median = mean(.$value, na.rm = TRUE),
                     eastness_sd = sd(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_east_cell, .)
  # bind to results
  east_locs <- bind_rows(east_locs, east_cell)
}


# prediction surface covariates

pland_elev_nort_east <- inner_join(pland_elev_nort, east_locs,  by = c("locality_id", "year"))
write_csv(pland_elev_nort_east, paste0("data/pland-elev_nort_east_location_", species, "_MEDS.csv"))

# ------------------------------------------
# add slope to checklists 
# ------------------------------------------

#repeat with slope
slope <- raster("~/GitHub/gulf_of_mexico_slr_birds/data/slope_1KMmn_GMTEDmd.tif")
slope <-  slope %>% 
  crop(., gcoos_extent) %>% 
  projectRaster(crs = crs(landcover))

slope_locs <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- ebird_buff$data[[which(yr == ebird_buff$year_lc)]]
  locs_slope_cell <- st_set_geometry(regions, NULL) 
  slope_cell <- exact_extract(slope, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(slope_median = mean(.$value, na.rm = TRUE),
                     slope_sd = sd(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_slope_cell, .)
  # bind to results
  slope_locs <- bind_rows(slope_locs, slope_cell)
}



# covariates


pland_elev_nort_east_slope <- inner_join(pland_elev_nort_east, slope_locs,  by = c("locality_id", "year"))
write_csv(pland_elev_nort_east_slope, paste0("data/pland-elev_nort_east_slope_location_", species, "_MEDS.csv"))


# ------------------------------------------
# start to add covariates to prediction surface 
# ------------------------------------------


r <- raster("data/prediction-surface-gcoos.tif")
# buffer to create neighborhood (6 km by 6 km) around each checklist point 
# with a radius equal to 100/2 CDL cells 
neighborhood_radius <- 200 * ceiling(max(res(landcover))) / 2
neighborhood_radius

cell_buff <- NULL
for (yr in years) {
  r_centers <- rasterToPoints(r, spatial = TRUE) %>% 
    st_as_sf() %>% 
    transmute(id = row_number()) %>% 
    mutate(year=as.character(yr), year_lc = paste0("y", yr))
  r_cells <- st_buffer(r_centers, dist = neighborhood_radius)
  cell_buff <- bind_rows(cell_buff, r_cells)
  # nest by year to  match to land cover data from the corresponding year
}
cell_buff <- nest(cell_buff, data = c(year, id, geometry))

lc_extract_pred <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year
  regions <- cell_buff$data[[which(yr == cell_buff$year_lc)]]
  # get landcover values within each buffered checklist area
  cc <- exact_extract(landcover[[yr]], regions, progress = FALSE)
  # count the number of each landcover class for each checklist buffer
  cc_count <- map(cc, ~ count(., landcover = value))
  # attach the year and locality id back to the checklists
  cc_summ <- tibble(st_drop_geometry(regions), data = cc_count) %>% 
    unnest(data)
  # bind to results
  lc_extract_pred <- bind_rows(lc_extract_pred, cc_summ)
}
# ------------------------------------------
# creat pland with c-cap first
# ------------------------------------------

# calculate the percent for each landcover class
pland_pred <- lc_extract_pred %>% 
  group_by(id, year) %>% 
  mutate(pland = n / sum(n)) %>% 
  ungroup() %>% 
  dplyr::select(-n) %>% 
  # remove NAs after tallying so pland is relative to total number of cells
  filter(!is.na(landcover))

#change landcover class name to ccap landuse names
pland_pred_ccap <- inner_join(pland_pred, ccap_classes,  by = "landcover")

# tranform to wide format, filling in implicit missing values with 0s%>% 


pland_pred_wide <- pland_pred_ccap %>% 
  select(-landcover) %>%
  group_by(year, id, landclass) %>% 
  summarise(pland = sum (pland)) %>%     # Summarize pland to calculate total per landclass
  ungroup() %>%                             
  pivot_wider(names_from = landclass, 
              values_from = pland, 
              values_fill = list(pland = 0))


# save

write_csv(pland_pred_wide, "data/covariates/pland_pred_6km_MEDS.csv")


# ------------------------------------------
# add elevation to the prediction surface
# ------------------------------------------

locs_pred <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff$data[[which(yr == cell_buff$year_lc)]]
  locs_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  elev_cell <- exact_extract(elev, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(elevation_median = mean(.$value, na.rm = TRUE),
                     elevation_sd = sd(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_cell, .)
  # bind to results
  locs_pred <- bind_rows(locs_pred, elev_cell)
}
# prediction surface covariates
pland_elev_pred <- inner_join(pland_pred_wide, locs_pred, by = c("id", "year"))
write_csv(pland_elev_pred, "data/pland-elev_prediction-surface_6km_MEDS.csv")

# ------------------------------------------
# add canopy to checklists and pred surface
# ------------------------------------------



# Load the archive package
library(archive)

# Open the ZIP file (stream it)
# Define the folder containing the canopy ZIP files
zip_folder <- "data/canopy"  # Replace with your folder path
# Create a folder to extract the files
unzip_folder <- file.path(zip_folder, "unzipped")
dir_create(unzip_folder)  # Create the unzip folder if it doesn't exist

# List all ZIP files in the folder
zip_files <- dir(zip_folder, pattern = "\\.zip$", full.names = TRUE)

for (zip_file in zip_files) {
  # Unzip using the system's unzip command
  system(paste("unzip", zip_file, "-d", unzip_folder))
}

# Initialize an empty list to store raster files
raster_list <- list()


# Find the .tif files inside the unzipped folder
tif_files <- dir(unzip_folder, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
# Load all the rasters
rasters <- lapply(tif_files, rast)

# Create a mosaic from the rasters
mosaic_raster <- do.call(mosaic, rasters)



# Optionally save the raster stack as a new file
writeRaster(mosaic_raster, "data/canopy_rasters.tif", overwrite = TRUE)


# ------------------------------------------
# add northness to pred surface
# ------------------------------------------

#repeat for prediction surface
nort_pred <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff$data[[which(yr == cell_buff$year_lc)]]
  locs_nort_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  nort_cell <- exact_extract(nort, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(northness_median = mean(.$value, na.rm = TRUE),
                     northness_sd = sd(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_nort_cell, .)
  # bind to results
  nort_pred <- bind_rows(nort_pred, nort_cell)
}




# prediction surface covariates

pland_elev_nort_pred <- inner_join(pland_elev_pred, nort_pred,  by = c("id", "year"))
write_csv(pland_elev_nort_pred, "data/pland-elev_nort_prediction-surface_MEDS.csv")

# ------------------------------------------
# add eastness to pred surface
# ------------------------------------------

#repeat for prediction surface
east_pred <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff$data[[which(yr == cell_buff$year_lc)]]
  locs_east_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  east_cell <- exact_extract(east, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(eastness_median = mean(.$value, na.rm = TRUE),
                     eastness_sd = sd(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_east_cell, .)
  # bind to results
  east_pred <- bind_rows(east_pred, east_cell)
}

# prediction surface covariates


pland_elev_nort_east_pred <- inner_join(pland_elev_nort_pred, east_pred,  by = c("id", "year"))
write_csv(pland_elev_nort_east_pred, "data/pland-elev_nort_east_prediction-surface_MEDS.csv")

# ------------------------------------------
# add slope to pred surface
# ------------------------------------------


slope_pred <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff$data[[which(yr == cell_buff$year_lc)]]
  pred_slope_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  slope_cell <- exact_extract(slope, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(slope_median = mean(.$value, na.rm = TRUE),
                     slope_sd = sd(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(pred_slope_cell, .)
  # bind to results
  slope_pred <- bind_rows(slope_pred, pred_slope_cell)
}



# prediction surface covariates


pland_elev_nort_east_slope_pred <- inner_join(pland_elev_nort_east_pred, slope_pred,  by = c("id", "year"))
write_csv(pland_elev_nort_east_slope_pred, "data/pland-elev_nort_east_slope_prediction-surface_MEDS.csv")




#repeat with canopy height
canopy <- raster(paste0(wd, "data/canopy_height_2019.tif"))
canopy <-  canopy %>% 
  crop(., cv_extent) %>% 
  projectRaster(crs = projection(landcover))


canopy_pred_cdl <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff_cdl$data[[which(yr == cell_buff_cdl$year_lc)]]
  locs_canopy_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  canopy_cell <- exact_extract(canopy, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(canopy = mean(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_canopy_cell, .)
  # bind to results
  canopy_pred_cdl <- bind_rows(canopy_pred_cdl, canopy_cell)
}
#replace nan with 0
canopy_pred_cdl <- as.data.frame(canopy_pred_cdl)
canopy_pred_cdl <- rapply(canopy_pred_cdl, f=function(x) ifelse(is.nan(x),0,x), how="replace" )

# Get the list of all .tif file paths in the folder
tif_files <- list.files(canopy_folder, pattern = "\\.tif$", full.names = TRUE)

# Create a list of raster objects using `rast()`
canopy_list <- lapply(tif_files, rast)

# Optionally, you can name the list with the file names (without the path and extension)
names(canopy_list) <- basename(tif_files)


canopy_pred <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff$data[[which(yr == cell_buff$year_lc)]]
  locs_canopy_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  for (canopy in canopy_list) {
  canopy_cell <- exact_extract(canopy, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(canopy = mean(.$value, na.rm = FALSE))) 
  }
  canopy_cell <- canopy_cell %>% 
    # join to lookup table to get id
    bind_cols(locs_canopy_cell, .)
  # bind to results
  canopy_pred <- bind_rows(canopy_pred, canopy_cell)
}
#replace nan with 0
canopy_pred <- as.data.frame(canopy_pred)
canopy_pred <- rapply(canopy_pred, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
# prediction surface covariates

pland_elev_nort_east_pred_crop<-read.csv("data/pland-elev_nort_east_prediction-surface_pmp_dream.csv")
pland_elev_nort_east_canopy_pred_cdl <- inner_join(pland_elev_nort_east_pred_crop, canopy_pred_cdl,  by = c("id", "year"))
write_csv(pland_elev_nort_east_canopy_pred_cdl, "data/pland-elev_nort_east_canopy_prediction-surface_pmp_dream.csv")

# ------------------------------------------
# ------------------------------------------
# perturb with pmp dust scenario
# ------------------------------------------
# ------------------------------------------




crop_change_pmp <- st_read("/Users/lily/Library/CloudStorage/Box-Box/calvin/gwregion_crop_change_BBAU.shp")
pland_geometry_10years <-NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year
  regions <- cell_buff_cdl$data[[which(yr == cell_buff_cdl$year_lc)]]
  pland_geometry_year <- inner_join(pland_pred_crop_category_cdl_wide, regions, by=c("id", "year"))
  pland_geometry_10years <-bind_rows(pland_geometry_10years, pland_geometry_year)
}
pland_geometry_10years<-st_as_sf(pland_geometry_10years)
crop_change_pmp <- st_transform(crop_change_pmp, crs=4326)
pland_geometry_10years <- st_transform(pland_geometry_10years, crs=4326)

st_is_valid(pland_geometry_10years) #test geometry values are valid
# Perform intersection to extract values from the shapefile based on spatial data frame's geometry
intersection_result <- st_intersection(pland_geometry_10years, crop_change_pmp)
names(intersection_result)
crop_change<-NULL
for (cp in names(crop_change_pmp)[2:8]) {
  crop_change <- intersection_result %>%
    mutate(!!sym(cp) := !!sym(cp) * .[[paste0(cp,".", 1)]])
  
}

crop_change_tb<-st_drop_geometry(crop_change)
write.csv(crop_change_tb, file="crop_change_pmp_10years_pred_bbau.csv")
crop_change_tb$year <- as.character(crop_change_tb$year)
crop_change_pmp <-NULL
for (yr in 2011:2020) {
  # get the buffered checklists for a given year
  regions <- cell_buff_cdl$data[[yr-2010]]
  crop_change_year <-inner_join(crop_change_tb, regions, by=c("id", "year"))
  crop_change_pmp <-bind_rows(crop_change_pmp, crop_change_year)
}
crop_change_pmp_slim <- crop_change_pmp[-18:-25]






# join in coordinates
r <- raster("data/prediction-surface-cv.tif")


crop_pmp_coords <- NULL
for (yr in 2011:2020) {
  pland_year <- crop_change_pmp_slim %>% filter(year == as.character(yr))
  pland_coords_year <- pland_year %>% st_as_sf()%>%st_transform(crs = 4326) %>% 
    st_coordinates() %>% 
    as.data.frame() %>% 
    rename(longitude = X, latitude = Y)  %>% 
    cbind(id=pland_year$id)%>%
    inner_join(pland_year, by = "id")
  crop_pmp_coords <- bind_rows(crop_pmp_coords, pland_coords_year)
}


elev <- raster(paste0(wd, "data/elevation_1KMmd_GMTEDmd.tif"))
locs_pred_pmp <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff_cdl$data[[which(yr == cell_buff_cdl$year_lc)]]
  locs_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  elev_cell <- exact_extract(elev, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(elevation_median = mean(.$value, na.rm = TRUE),
                     elevation_sd = sd(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_cell, .)
  # bind to results
  locs_pred_pmp <- bind_rows(locs_pred_pmp, elev_cell)
}
# prediction surface covariates
pland_elev_pred_pmp <- inner_join(crop_pmp_coords, locs_pred_pmp, by = c("id", "year"))
nrow(pland_elev_pred_pmp) 
pland_elev_pred_pmp_1 <- pland_elev_pred_pmp[1:7500000,]
pland_elev_pred_pmp_2 <-pland_elev_pred_pmp[7500001:14985850,]
write_csv(pland_elev_pred_pmp_1, "data/pland-elev_prediction-surface_pmp_10year_bbau_1.csv")
write_csv(pland_elev_pred_pmp_2, "data/pland-elev_prediction-surface_pmp_10year_bbau_2.csv")

# ------------------------------------------
# add northness to checklists and pred surface
# ------------------------------------------
#repeat the same process for northness covariates
nort <- raster(paste0(wd, "data/northness_1KMmd_GMTEDmd.tif"))
nort <-  nort %>% 
  crop(., cv_extent) %>% 
  projectRaster(crs = projection(landcover))



#repeat for prediction surface
nort_pred <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff_cdl$data[[which(yr == cell_buff_cdl$year_lc)]]
  locs_nort_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  nort_cell <- exact_extract(nort, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(northness_median = mean(.$value, na.rm = TRUE),
                     northness_sd = sd(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_nort_cell, .)
  # bind to results
  nort_pred <- bind_rows(nort_pred, nort_cell)
}




# prediction surface covariates

pland_elev_nort_pred_crop <- inner_join(pland_elev_pred_pmp, nort_pred,  by = c("id", "year"))


# ------------------------------------------
# add eastness to checklists and pred surface
# ------------------------------------------

#repeat with eastness
east <- raster(paste0(wd, "data/eastness_1KMmd_GMTEDmd.tif"))
east <-  east %>% 
  crop(., cv_extent) %>% 
  projectRaster(crs = projection(landcover))



#repeat for prediction surface
east_pred <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff_cdl$data[[which(yr == cell_buff_cdl$year_lc)]]
  locs_east_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  east_cell <- exact_extract(east, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(eastness_median = mean(.$value, na.rm = TRUE),
                     eastness_sd = sd(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_east_cell, .)
  # bind to results
  east_pred <- bind_rows(east_pred, east_cell)
}



# prediction surface covariates

pland_elev_nort_east_pred_crop <- inner_join(pland_elev_nort_pred_crop, east_pred,  by = c("id", "year"))

#repeat with canopy height
canopy <- raster(paste0(wd, "data/canopy_height_2019.tif"))
canopy <-  canopy %>% 
  crop(., cv_extent) %>% 
  projectRaster(crs = projection(landcover))


canopy_pred_cdl <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff_cdl$data[[which(yr == cell_buff_cdl$year_lc)]]
  locs_canopy_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  canopy_cell <- exact_extract(canopy, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(canopy = mean(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_canopy_cell, .)
  # bind to results
  canopy_pred_cdl <- bind_rows(canopy_pred_cdl, canopy_cell)
}
#replace nan with 0
canopy_pred_cdl <- as.data.frame(canopy_pred_cdl)
canopy_pred_cdl <- rapply(canopy_pred_cdl, f=function(x) ifelse(is.nan(x),0,x), how="replace" )



# prediction surface covariates
# pland_elev_nort_east_pred_crop <- read.csv("data/pland-elev_nort_east_prediction-surface_pmp_dust.csv")
# canopy_pred_cdl$id <- as.character(canopy_pred_cdl$id)
# pland_elev_nort_east_pred_crop$year<-as.character(pland_elev_nort_east_pred_crop$year)
pland_elev_nort_east_canopy_pred_cdl <- inner_join(pland_elev_nort_east_pred_crop, canopy_pred_cdl,  by = c("id", "year"))

pland_elev_nort_east_canopy_pred_cdl_1 <- pland_elev_nort_east_canopy_pred_cdl[1:7500000,]
pland_elev_nort_east_canopy_pred_cdl_2 <- pland_elev_nort_east_canopy_pred_cdl[7500001:14985850,]


write_csv(pland_elev_nort_east_canopy_pred_cdl_1, "data/pland-elev_nort_east_canopy_prediction-surface_pmp_bbau_1.csv")
write_csv(pland_elev_nort_east_canopy_pred_cdl_2, "data/pland-elev_nort_east_canopy_prediction-surface_pmp_bbau_2.csv")



pland_elev_nort_east_pred_crop <-read.csv("data/pland-elev_nort_east_prediction-surface_pmp_dust.csv")
pland_elev_nort_east_canopy_pred_cdl <- inner_join(pland_elev_nort_east_pred_crop, canopy_pred_cdl,  by = c("id", "year"))
write_csv(pland_elev_nort_east_canopy_pred_cdl, "data/pland-elev_nort_east_canopy_prediction-surface_pmp_dust.csv")


# -----------------------------------------------------
# create maps of both cdl and dust pland for comparison
# -----------------------------------------------------
crop_change_tb <- read.csv("crop_change_pmp_10years.csv")
crop_cdl <- read.csv("data/pland_elev_greyel.csv")


#Convert data frame format to spatial format of pland metrics
land_cover_pmp <- crop_pmp_coords %>% 
  filter(year == 2015) %>%
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(landcover)) 
land_cover_cdl <- crop_cdl_coords %>% 
  filter(year == 2015) %>%
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(landcover)) 

# rasterize points
#names(rice_cover)[names(rice_cover)=="3"] <- "rice"
land_cover_pmp_raster <-  rasterize(land_cover_pmp, r, field = c("rice", "corn", "pasture","grain", "alfalfa", "field_row", "perennial"))
land_cover_cdl_raster <-  rasterize(land_cover_cdl, r, field = c("rice", "corn", "pasture","grain", "alfalfa", "field_row", "perennial"))

# Extract raster layers from the stacks (replace with the layer indices you want to plot)

# # Create a blank plot to hold all the histograms
# plot.new()
# 
# # Set up the plotting layout
# par(mfrow = c(4, 2), mar = c(3, 3, 3, 3), oma = c(3, 3, 3, 3))
# 
# for (i in 1:7){
# pmp <- land_cover_pmp_raster[[i]]
# cdl <- land_cover_cdl_raster[[i]]
# 
# # Extract values from the raster layers
# values1 <- values(pmp)
# values2 <- values(cdl)
# 
# # Manually define breaks for histograms
# breaks1 <- seq(0, 1, by = 0.1)  # Adjust the breaks as needed
# breaks2 <- seq(0, 1, by = 0.1)    # Adjust the breaks as needed
# 
# # Create histograms with manual breaks
# hist1 <- hist(values1, breaks = breaks1, plot = FALSE)
# hist2 <- hist(values2, breaks = breaks2, plot = FALSE)
# # Determine the common x-axis limits for both histograms
# #x_limits <- range(c(hist1$breaks, hist2$breaks))
# 
# # Combine the histogram data
# hist_data <- rbind(hist1$counts, hist2$counts)
# 
# 
# 
# # Set the size of the plot device (adjust width and height as needed)
# #png("histogram_plot.png", width = 1000, height = 600, res = 96, units = "px")
# 
# # Set up the plot with custom size
# # Adjust margins (bottom, left, top, right)
# # Set up the plot
# barplot(hist_data, beside = TRUE, col = c("blue", "red"), xlim = c(0, ncol(hist_data) + 15),
#         names.arg = hist1$mids, legend.text = c("pmp", "cdl"),
#         main = names(land_cover_pmp_raster)[i], xlab = "Land cover class proportion in each cell", ylab = "Frequency")
# 
#   # Reset to the default layout
# 
# }
# # Add an overall title to the entire plot
# title("Histograms of Land Cover Proportions")
# 
# # Reset the plotting layout
# par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
# 
# # Save or display the plot
# #dev.off()  # If saving to a file
# 
# # Reset the plotting parameters
# par(mfrow = c(1, 1))





# Create a blank plot to hold all the histograms
plot.new()


for (i in 1:7){
  pmp <- land_cover_pmp_raster[[i]]
  cdl <- land_cover_cdl_raster[[i]]
  
  # Extract values from the raster layers
  values1 <- values(pmp)
  values2 <- values(cdl)
  
  # Manually define breaks for histograms
  breaks1 <- seq(0, 1, by = 0.1)  # Adjust the breaks as needed
  breaks2 <- seq(0, 1, by = 0.1)    # Adjust the breaks as needed
  
  # Create histograms with manual breaks
  hist1 <- hist(values1, breaks = breaks1, plot = FALSE)
  hist2 <- hist(values2, breaks = breaks2, plot = FALSE)
  # Determine the common x-axis limits for both histograms
  #x_limits <- range(c(hist1$breaks, hist2$breaks))
  
  # Combine the histogram data
  hist_data <- rbind(hist1$counts, hist2$counts)
  
  # Set up the plotting layout
  par(mfrow = c(1, 1), mar = c(3, 3, 3, 3), oma = c(3, 3, 3, 3))
  
  
  # Set the size of the plot device (adjust width and height as needed)
  #png("histogram_plot.png", width = 1000, height = 600, res = 96, units = "px")
  
  # Set up the plot with custom size
  # Adjust margins (bottom, left, top, right)
  # Set up the plot
  barplot(hist_data, beside = TRUE, col = c("blue", "red"), xlim = c(0, ncol(hist_data) + 15),
          names.arg = hist1$mids, legend.text = c("pmp", "cdl"),
          main = names(land_cover_pmp_raster)[i], xlab = "Land cover class proportion in each cell", ylab = "Frequency")
  
  # Reset to the default layout
  
}
# Add an overall title to the entire plot
title("Histograms of Land Cover Proportions")

# Reset the plotting layout
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))

# Save or display the plot
#dev.off()  # If saving to a file

# Reset the plotting parameters
par(mfrow = c(1, 1))



# project to albers equal-area for mapping
land_cover_pmp_raster <- land_cover_pmp_raster %>% projectRaster(crs = "ESRI:102003") %>% 
  trim()
land_cover_cdl_raster <- land_cover_cdl_raster %>% projectRaster(crs = "ESRI:102003") %>% 
  trim()

# load gis data for making maps
map_proj <- map_proj <- "ESRI:102003"
ne_land <- read_sf("data/gis-data.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
cv <- read_sf("data/gis-data.gpkg", "cv") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_country_lines <- read_sf("data/gis-data.gpkg", "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

pmp_points <- crop_pmp_coords %>% 
  filter(year == 2015) %>%
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  dplyr::select(rice, corn, pasture, alfalfa, perennial)
cdl_points <- crop_cdl_coords %>% 
  filter(year == 2015) %>%
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  dplyr::select(rice, corn, pasture, alfalfa, perennial)

both_pts<- list(pmp=pmp_points, cdl = cdl_points)

# # Source add.alpha function from Github to brew alpha color
# require(RCurl)
# source(textConnection(getURL("https://gist.github.com/mages/5339689/raw/576263b8f0550125b61f4ddba127f5aa00fa2014/add.alpha.R")))
# myColours = c(1, "steelblue", "#FFBB00", rgb(0.4, 0.2, 0.3))
# myColoursAlpha <- add.alpha(myColours, alpha=0.4)
### "#00000066" "#4682B466" "#FFBB0066" "#66334D66" 

# map
p <- par(mfrow = c(1, 2))
for (i in seq_along(both_pts)) {
  par(mar = c(0.25, 0.25, 0.25, 0.25))
  # set up plot area
  # plot(cv, col = NA, border = NA)
  # plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
  plot(st_geometry(both_pts[[i]]), col = NA)
  #borders
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
  plot(cv, col = "#cccccc", border = "#000000", lwd = 1, add = TRUE)
  plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
  
  
  # ebird observations
  # rice
  plot(filter(both_pts[[i]],rice>0.3),
       pch = 19, cex = 0.1, col = alpha("#4682B466", 0.25),
       add = TRUE)
  # corn
  plot(filter(both_pts[[i]],corn>0.3) %>% st_geometry(),
       pch = 19, cex = 0.3, col = alpha("#4daf4a", 0.5),
       add = TRUE)
  # pasture
  plot(filter(both_pts[[i]],pasture>0.3) %>% st_geometry(),
       pch = 19, cex = 0.3, col = alpha("#FFBB0066", 0.5),
       add = TRUE)
  # alfalfa
  plot(filter(both_pts[[i]],alfalfa>0.3) %>% st_geometry(),
       pch = 19, cex = 0.3, col = alpha("#66334D66", 0.5),
       add = TRUE)
  # perennial
  plot(filter(both_pts[[i]],perennial>0.3) %>% st_geometry(),
       pch = 19, cex = 0.3, col = alpha("#00000066", 0.5),
       add = TRUE)
  # legend
  
  legend("bottomright", bty = "n",
         col = c("#4682B466", "#4daf4a", "#FFBB0066", "#66334D66","#00000066"),
         legend = c("Rice", "Corn", "Pasture", "Alfalfa", "Perennial"),
         pch = 19)
  
  box()
  par(new = TRUE, mar = c(0, 0, 3, 0))
  if (names(both_pts)[i] == "PMP") {
    title("PMP")
  } else {
    title("CDL")
  }
}
par(p)




both_ras<- stack(land_cover_pmp_raster, land_cover_cdl_raster)
names(land_cover_dust_raster)
# map
p <- par(mfrow = c(5, 2))
par(mar = c(0.25, 0.5, 0.25, 0.25))
for(i in 1:5) {
  plot(land_cover_dust_raster[[i]], col= viridis(10), axes = FALSE, box = FALSE)
  plot(land_cover_cdl_raster[[i]], col= viridis(10),  axes = FALSE, box = FALSE)
  title(main=names(land_cover_dust_raster[[i]]), line=-8*i, outer=TRUE)
  
  
}
par(p)


levelplot(both_ras, par.settings=rasterTheme(viridis_pal(option = "D")(255)),
          colorkey=list(labels=list(at=seq(0.1,1,0.1), 
                                    labels=seq(0.1,1,0.1)),
                        space="bottom", height=0.6),
          margin=FALSE, contour = FALSE, 
          scales=list(draw=FALSE))




# # Define the raster extent and resolution
# raster_extent <- extent(pland_geometry)
# raster_resolution <- max(res(landcover))  # Adjust the resolution as needed
# 
# # Create an empty raster with the specified extent and resolution
# empty_raster <- raster(raster_extent, res = raster_resolution)
# 
# # Choose the field to use for rasterizing
# fields_to_rasterize <- names(crop_change_pmp)[2:8] # Replace with the actual field name
# 
# # Create a list to store the rasterized fields
# raster_list <- list()
# 
# # Rasterize each field and add it to the list
# for (field in fields_to_rasterize) {
#   rasterized <- rasterize(crop_change_pmp, empty_raster, field)
#   raster_list[[field]] <- rasterized
# }
# 
# # Create a raster stack from the list of rasterized fields
# crop_change_raster <- stack(raster_list)
# proj4string(crop_change_raster) <- crs(landcover)
# 
# # Save the raster stack to the specified file
# writeRaster(crop_change_raster, filename = "crop_change_raster.tif", format = "GTiff", overwrite = TRUE)
# 
# min(pland_geometry$year)
# ebird_buff_pmp <- NULL
# for (cp in names(crop_change_raster)) {
#   ebird_buff_crop <- ebird %>% 
#     distinct(crop = cp,
#              locality_id, latitude, longitude)
#   ebird_buff_pmp <- rbind(ebird_buff_pmp, ebird_buff_crop)
# }
# 
# 
# 
# 
# #define diameter
# neighborhood_radius <- 100 * ceiling(max(res(landcover))) / 2
# 
# 
# # convert to spatial features
# ebird_buff_pmp <- ebird_buff_pmp %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
#   # transform to cdl projection
#   st_transform(crs = projection(crop_change_raster)) %>% 
#   st_buffer(dist = neighborhood_radius) %>% 
#   # nest by year to  match to land cover data from the corresponding year
#   nest(data = c(locality_id, geometry))
# 
# # iterate over all species extracting abundance for all checklists in each buff
# dland <- NULL
# for (cp in names(crop_change_raster)) {
#   # get the buffered checklists for a given year
#   regions <- ebird_buff_pmp$data[[which(cp == ebird_buff_pmp$crop)]]
#   #locs <- st_set_geometry(regions, NULL) 
#   # get landcover values within each buffered checklist area
#   dl <- exact_extract(crop_change_raster[[cp]], regions, progress = FALSE)
#   # count the number of each landcover class for each checklist buffer
#   dl_checklists <- dl %>% 
#     map_dfr(~ tibble(dland = mean(.$value, na.rm = TRUE), crop=cp)) %>% 
#     # join to lookup table to get locality_id
#     bind_cols(regions, .)
#   # attach the year and locality id back to the checklists
#   #sp_summ <- tibble(st_drop_geometry(regions), data = sp_checklists) %>% 
#   #  unnest(data)
#   # bind to results
#   dland <- bind_rows(dland, dl_checklists)
#   
# }
# 
# write.csv(dland, file="extract_dland_table.csv", row.names = FALSE) 
# dland <- st_drop_geometry(dland)
# dland<-read_csv("extract_dland_table.csv")
# dland_rmna <- na.omit(dland)
# dland$dland <- ifelse(is.na(dland$dland), 1, dland$dland)
# # Sample data frames
# pland_crop_change <- NULL
# 
# for(cp in names(crop_change_raster)) {
#   df1 <- dland%>%filter(crop==cp)%>%
#     dplyr::select(-crop)
#   
#   df2 <- pland_crop_category_wide
#   
#   # Merge data frames based on the common 'ID' column
#   dcrop <- pland_crop_category_wide %>% merge(df1, by = "locality_id", all = TRUE) %>%
#     map_dfr(~ tibble(cp=dland*pland_crop_category_wide$cp)) %>%
#  
# summary(dland$dland)
# }