library(auk)
library(lubridate)
library(sf)
library(gridExtra)
library(tidyverse)
library(raster)
library(terra) ## confirm the GDAL version being used
library(viridis)
library(PresenceAbsence)
library(verification)
library(fields)
library(dggridR)
library(raster)
library(ncdf4)
library(foreign)
library(fs)
library(tools)
library(dplyr)
library(rhdf5)
setwd("/Users/liyingnceas/GitHub/gulf_of_mexico_slr_birds/")
select <- dplyr::select

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rhdf5")
#-------------------------------------
#load all the environmental covariates
#-------------------------------------

# Open the NetCDF file as a SpatRaster (terra object)
elev <- rast("data/northern_gulf_1_navd88_2010.nc")
# Assign CRS to the raster
crs(elev) <- "EPSG:4326"


# Get the bounding box as an extent object
bbox_extent <- ext(elev)

# Optionally, convert the extent to a vector
bbox_vector <- as.vector(bbox_extent)
names(bbox_vector) <- c("xmin", "xmax", "ymin", "ymax")
bbox_vector["ymin"] <- 30

ext(elev)<-ext(bbox_vector)

# Extract the tar.gz file
untar("~/R/0170338.1.1.tar.gz", exdir = "/home/lli/R/data")




# Define a new color palette for the elevation raster (e.g., "viridis" or custom)
elev_colors <- viridis(100)  # Yellow-Green-Blue palette

# Plot the elevation raster with the new color palette
plot(elev, main = "Elevation and NGOM", col = elev_colors)

# Add the NGOM shapefile on top with transparent fill and boundary only
plot(gcoos, add=TRUE, col = "transparent", lwd = 2)  # Transparent fill with blue boundary

# Add legend for elevation raster
legend("topright", legend = "Elevation (meters)", fill = elev_colors, title = "Elevation")

# Add legend for NGOM shapefile (transparent fill, boundary shown)
legend("bottomright", legend = "NGOM shoreline", col = elev_colors, lty = 1, lwd = 2, title = "Shoreline Boundary")

# List the land cover files
#landcover_files <- list.files("data/landcover", pattern = "ccap_landcover_20200311", full.names = TRUE)
landcover_files <- list.files("data/landcover", pattern = "2016_ccap_landcover_20200311", full.names = TRUE)
landcover_files 
# Load the rasters as SpatRaster objects
landcover_rasters <- lapply(landcover_files, rast)
# #------------use this if involving multiple years landcover---------------------
# # Get the extents of each raster
# extents <- lapply(landcover_rasters, ext)
# 
# # Calculate the minimum and maximum values for x and y (extent coordinates)
# xmin <- max(sapply(extents, function(e) e[1]))  # maximum of xmin
# xmax <- min(sapply(extents, function(e) e[2]))  # minimum of xmax
# ymin <- max(sapply(extents, function(e) e[3]))  # maximum of ymin
# ymax <- min(sapply(extents, function(e) e[4]))  # minimum of ymax
# 
# # Create a SpatExtent object from these values
# common_extent <- ext(xmin, xmax, ymin, ymax)
# # Define bounding box coordinates
# bbox_coords <- c(xmin = -98, xmax = -79, ymin = 22, ymax = 32)
# # Now, crop each raster to the common extent
# cropped_rasters <- lapply(landcover_rasters, crop, common_extent)
# # Stack the cropped rasters
# landcover_stack <- rast(cropped_rasters)
# 
# # label layers with years of interest
# landcover <- setNames(landcover_stack, c("y2019", "y2020", "y2021"))


# get 2016 land cover used for 2020, the closest ccap landcover available
l <- landcover_rasters[[1]]
landcover<-setNames(l, "y2020")




#load gcoos boundary and crop elev 
gcoos <-read_sf("data/gis-data-gcoos.gpkg", "gcoos")%>%
  # project to the cdl projection
  st_transform(crs = 4326)

# with a radius equal to 100  c-cap cells 
neighborhood_radius <- 100 * ceiling(max(res(landcover))) / 2

agg_factor <- round(2 * neighborhood_radius / res(landcover))




elev <-  elev %>% 
  crop(., gcoos_extent) 



extent(gcoos)

# Get the bounding box
bbox <- st_bbox(gcoos)

# Modify ymax
bbox["ymax"] <- 31.16

# Update the sf object with the modified bbox
gcoos <- st_crop(gcoos, bbox)

# crop, buffer cv_cov by 10 km to provide a little wiggly room
gcoos_extent <- gcoos %>% 
  st_buffer(dist = 10000) %>% 
  st_transform(crs = crs(elev))
gcoos_crop<-st_transform(gcoos, crs= crs(landcover))

# Plot the landcover raster
plot(landcover_cropped)

# Overlay the gcoos_crop polygon boundary for reference
plot(gcoos, add=TRUE, border = 'red')

# Crop landcover to the bounding box of gcoos_extent
landcover_cropped <- crop(landcover, gcoos_crop)
writeRaster(landcover_cropped, "data/landcover_cropped_no_blank.tif", overwrite = TRUE)



#repeat the same process for northness covariates
nort <- raster("data/northness_1KMmn_GMTEDmd.tif")
nort <-  nort %>% 
  crop(., gcoos_extent) %>% 
  projectRaster(crs = crs(elev)) 

temperature <- raster("data/t_raster.tif")
temperature <-  temperature %>% 
  crop(., gcoos_extent) %>% 
  projectRaster(crs = crs(elev)) 

trange <- raster("data/t_range_raster.tif")
trange <-  trange %>% 
  crop(., gcoos_extent) %>% 
  projectRaster(crs = crs(elev)) 

bio1 <- raster("data/bio1_585_raster.tif")
bio1 <-  bio1 %>% 
  crop(., gcoos_extent) %>% 
  projectRaster(crs = crs(elev)) 

bio7 <- raster("data/bio7_585_raster.tif")
bio7 <-  bio7 %>% 
  crop(., gcoos_extent) %>% 
  projectRaster(crs = crs(elev)) 


#repeat with eastness
east <- raster("data/eastness_1KMmn_GMTEDmd.tif")
east <-  east %>% 
  crop(., gcoos_extent) %>% 
  projectRaster(crs = crs(elev))

slope <- raster("data/slope_1KMmn_GMTEDmd.tif")
slope <-  slope %>% 
  crop(., gcoos_extent) %>% 
  projectRaster(crs = crs(elev))

poc <- raster("data/normalized_poc_raster.tif")
poc <-  poc %>% 
  crop(., gcoos_extent) %>% 
  projectRaster(crs = crs(elev)) 

chl <- raster("data/normalized_chl_raster.tif")
chl <-  chl %>% 
  crop(., gcoos_extent) %>% 
  projectRaster(crs = crs(elev)) 

#load npp rasters and average across years
# List all HDF files in the folder
npp_files <- list.files("data/poc_viirs", pattern = "\\.nc$", full.names = TRUE)
npp<-rast("GEDI_L4B_Gridded_Biomass_V2_1_2299/data/GEDI04_B_MW019MW223_02_002_02_R01000M_MU.tif")


# Initialize an empty list to store rasters
raster_list <- list()

# Loop through each HDF file
for (file in npp_files) {
  # Read the raster from the HDF file
  raster <- rast(file)
  
  # If there are multiple layers, select the layer you want or just take the first one
  raster_list <- append(raster_list, list(raster))
}



# Stack all rasters together (if they have the same extent and resolution)
raster_stack <- rast(raster_list)
npp_ocean <- mean(raster_stack, na.rm = TRUE)


npp<-project(npp, crs(npp_ocean)) 
npp_resampled <- resample(npp, npp_ocean)
# Crop the rasters
npp_resampled_cropped <- crop(npp_resampled, gcoos_extent)
npp_ocean_cropped <- crop(npp_ocean, gcoos_extent)


min_npp_value <- as.numeric(global(npp_resampled_cropped, "min", na.rm = TRUE))
max_npp_value <- as.numeric(global(npp_resampled_cropped, "max", na.rm = TRUE))

npp_normalized <- app(npp_resampled_cropped, fun = function(x) {
  (x - min_npp_value) / (max_npp_value - min_npp_value)
})

min_npp_ocean_value <- as.numeric(global(npp_ocean_cropped, "min", na.rm = TRUE))
max_npp_ocean_value <- as.numeric(global(npp_ocean_cropped, "max", na.rm = TRUE))

npp_ocean_normalized <- app(npp_ocean_cropped, fun = function(x) {
  (x - min_npp_ocean_value) / (max_npp_ocean_value - min_npp_ocean_value)
})

npp_poc_raster <- app(
  c(npp_normalized, npp_ocean_normalized),
  fun = function(x) max(x, na.rm = TRUE),
  filename = "normalized_npp_raster.tif",
  overwrite = TRUE
)

# Replace -Inf (or other negative values) with NA
combined_raster_masked <- app(npp_poc_raster, fun = function(x) {
  x[is.infinite(x)] <- NA  # Replace -Inf with NA
  return(x)
})


# Assuming combined_raster is a SpatRaster object from terra
combined_raster_filled <- focal(combined_raster_masked, w = matrix(1, 3, 3), fun = function(x) {
  x[which.min(is.na(x))]  # Return the value of the nearest neighbor
})

plot(combined_raster_filled)

# Save the filled raster to a file
writeRaster(combined_raster_filled, "data/normalized_poc_raster.tif", overwrite = TRUE)



#Chlorophyll
#unRAR
# Path to your RAR file
rar_file <- "data/2020_LCC.rar"

# Path where you want to extract the contents
output_dir <- "data/"

system(paste("/opt/homebrew/bin/unar", shQuote(rar_file), "-o", shQuote(output_dir)))

#filter the files to cover GOM only
# Load required library
library(fs)  # For working with file paths

# Define the root directory containing the 12 folders
root_dir <- "data/2020_LCC"

# Define the tile identifiers to search for
tiles_to_find <- c("h09v05", "h09v06", "h10v05", "h10v06")

# Get all folder names in the root directory
all_folders <- dir(root_dir, full.names = TRUE)
# Extract month from folder names
month_indices <- substr(basename(all_folders), 5, 6)

# Initialize an empty list to store the results
monthly_tile_files <- vector("list", 6)

# Loop through each month (01 to 5) for the months all 4 tiles available
for (month in sprintf("%02d", 1:5)) {
  # Identify folders for the current month
  month_folders <- all_folders[month_indices == month]
  
  # Initialize a list for files in this month
  monthly_files <- list()
  
  # Loop through each folder for the current month
  for (folder in month_folders) {
    # Get all .tif files in the current folder
    all_files <- dir(folder, full.names = TRUE, pattern = "\\.tif$")
    
    # Find files matching the tiles
    matching_files <- lapply(tiles_to_find, function(tile) {
      grep(tile, all_files, value = TRUE)
    })
    
    # Append the matching files to the list
    monthly_files <- c(monthly_files, unlist(matching_files))
  }
  
  # Store the files for the current month
  monthly_tile_files[[as.numeric(month)]] <- monthly_files
}

# Add names to the list for clarity
names(monthly_tile_files) <- sprintf("Month%02d", 1:6)

# View the result
print(monthly_tile_files)

# Load required libraries
library(terra)
library(leaflet)


# Define output directory for mosaics
output_dir <- "data/LCC_mosaic"
dir.create(output_dir, showWarnings = FALSE)

# List of months
months <- sprintf("%02d", 1:5)

# Initialize a list to store Leaflet maps
leaflet_maps <- list()

# Loop through each month
for (month in months) {
  # Get the list of tiles for this month
  tile_files <- monthly_tile_files[[as.numeric(month)]]
  
  if (length(tile_files) > 0) {
    # Read all tiles as SpatRaster objects
    tiles <- lapply(tile_files, rast)
    
    # Mosaic the tiles
    mosaic_raster <- do.call(mosaic, c(tiles, list(fun = mean)))  # Use mean to handle overlaps
    
    # Define the output file name
    output_file <- file.path(output_dir, paste0("mosaic_", month, ".tif"))
    
    # Save the mosaic raster to disk
    writeRaster(mosaic_raster, output_file, overwrite = TRUE)
    
    
  }
}


# List all GeoTIFF files for the months
tif_files <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE)

# Load all the rasters as a stack
monthly_rasters <- rast(tif_files)

# Compute the average raster
LCC_average_raster <- mean(monthly_rasters, na.rm = TRUE)

# Save the average raster to a new GeoTIFF file
average_raster_file <- file.path(output_dir, "LCC_average_raster.tif")
crs(land_bio1_aligned) == crs(ocean_t_aligned)


#superinpose CHL with LCC

#load npp rasters and average across years
# List all HDF files in the folder
chl_files <- list.files("data/snpp_viirs", pattern = "\\.nc$", full.names = TRUE)

# Initialize an empty list to store rasters
raster_list <- list()

# Loop through each HDF file
for (file in chl_files) {
  # Read the raster from the HDF file
  raster <- rast(file)
  
  # If there are multiple layers, select the layer you want or just take the first one
  raster_list <- append(raster_list, list(raster))
}



# Stack all rasters together (if they have the same extent and resolution)
raster_stack <- rast(raster_list)
chl_ocean <- mean(raster_stack, na.rm = TRUE)


LCC<-project(LCC_average_raster, crs(chl_ocean)) 
LCC_resampled <- resample(LCC, chl_ocean)
# Crop the rasters
LCC_resampled_cropped <- crop(LCC_resampled, gcoos_extent)
chl_ocean_cropped <- crop(chl_ocean, gcoos_extent)

library(terra)

# Step 1: Compute Mean and Standard Deviation for LCC_resampled_cropped
mean_lcc <- global(LCC_resampled_cropped, fun = "mean", na.rm = TRUE)[1]
sd_lcc <- global(LCC_resampled_cropped, fun = "sd", na.rm = TRUE)[1]
mean_lcc_value <- as.numeric(mean_lcc)
sd_lcc_value <- as.numeric(sd_lcc)

# Ensure raster math for standardization
LCC_standardized <- app(LCC_resampled_cropped, fun = function(x) {
  (x - mean_lcc_value) / sd_lcc_value
})

# Step 1: Compute Mean and Standard Deviation for LCC_resampled_cropped
mean_chl <- global(chl_ocean_cropped, fun = "mean", na.rm = TRUE)[1]
sd_chl <- global(chl_ocean_cropped, fun = "sd", na.rm = TRUE)[1]
mean_chl_value <- as.numeric(mean_chl)
sd_chl_value <- as.numeric(sd_chl)

# Ensure raster math for standardization
chl_standardized <- app(chl_ocean_cropped, fun = function(x) {
  (x - mean_chl_value) / sd_chl_value
})



chl_combined_raster <- app(
  c(LCC_standardized, chl_standardized),
  fun = function(x) max(x, na.rm = TRUE),
  filename = "standardized_chl_raster.tif",
  overwrite = TRUE
)

plot(chl_combined_raster)

# Save the filled raster to a file
writeRaster(chl_combined_raster, "standardized_chl_raster.tif", overwrite = TRUE)



min_lcc_value <- as.numeric(global(LCC_resampled_cropped, "min", na.rm = TRUE))
max_lcc_value <- as.numeric(global(LCC_resampled_cropped, "max", na.rm = TRUE))

LCC_normalized <- app(LCC_resampled_cropped, fun = function(x) {
  (x - min_lcc_value) / (max_lcc_value - min_lcc_value)
})

min_chl_value <- as.numeric(global(chl_ocean_cropped, "min", na.rm = TRUE))
max_chl_value <- as.numeric(global(chl_ocean_cropped, "max", na.rm = TRUE))

chl_normalized <- app(LCC_resampled_cropped, fun = function(x) {
  (x - min_lcc_value) / (max_lcc_value - min_lcc_value)
})

chl_normalized_raster <- app(
  c(LCC_normalized, chl_normalized),
  fun = function(x) max(x, na.rm = TRUE),
  filename = "standardized_chl_raster.tif",
  overwrite = TRUE
)

plot(chl_normalized_raster)

# Save the filled raster to a file
writeRaster(chl_normalized_raster, "normalized_chl_raster.tif", overwrite = TRUE)




world_future<-rast("data/bioclim/wc2.1_30s_bioc_BCC-CSM2-MR_ssp245_2081-2100.tif")

#load rasters and average across months
# List all tif files in the folder
bioclim_files <- list.files("data/bioclim/wc2.1_30s_tavg/", pattern = "\\.tif$", full.names = TRUE)

# Initialize an empty list to store rasters
raster_list <- list()

# Loop through each HDF file
for (file in bioclim_files) {
  # Read the raster from the HDF file
  raster <- rast(file)
  
  # If there are multiple layers, select the layer you want or just take the first one
  raster_list <- append(raster_list, list(raster))
}



# Stack all rasters together (if they have the same extent and resolution)
raster_stack <- rast(raster_list)
land_bio1 <- mean(raster_stack, na.rm = TRUE)

ocean_t <-rast("data/bioclim/tas_baseline_2000_2020_depthsurf_0103_d100_33a3_U1734391611301.nc")



land_bio1<-project(land_bio1, crs(elev)) 
ocean_t<-project(ocean_t, crs(elev))
land_bio1_resampled <- resample(land_bio1, ocean_t)

# Crop the rasters
land_bio1_resampled<- crop(land_bio1_resampled, gcoos_extent)
plot(ocean_t) <- crop(ocean_t, gcoos_extent)


# Step 1: Compute Mean and Standard Deviation for LCC_resampled_cropped
min_land_bio1 <- global(land_bio1_resampled, fun = "min", na.rm = TRUE)[1]
max_land_bio1 <- global(land_bio1_resampled, fun = "max", na.rm = TRUE)[1]
min_land_bio1_value <- as.numeric(min_land_bio1)
max_land_bio1_value <- as.numeric(max_land_bio1)

land_bio1_normalized <- app(land_bio1_resampled, fun = function(x) {
  (x - min_land_bio1_value) / (max_land_bio1_value - min_land_bio1_value)
})

# Step 1: Compute Mean and Standard Deviation for LCC_resampled_cropped
min_ocean_t <- global(ocean_t, fun = "min", na.rm = TRUE)[1]
max_ocean_t <- global(ocean_t, fun = "max", na.rm = TRUE)[1]
min_ocean_t_value <- as.numeric(min_ocean_t)
max_ocean_t_value <- as.numeric(max_ocean_t)

ocean_t_normalized <- app(ocean_t, fun = function(x) {
  (x - min_ocean_t_value) / (max_ocean_t_value - min_ocean_t_value)
})



t_raster <- app(
  c(land_bio1_normalized, ocean_t_normalized),
  fun = function(x) max(x, na.rm = TRUE),
  filename = "t_raster.tif",
  overwrite = TRUE
)

#fill the no-value cells
t_raster <- app(t_raster, fun = function(x) {
  x[is.infinite(x)] <- NA  # Replace -Inf with NA
  return(x)
})


# fill NA
t_raster<- focal(t_raster, w = matrix(1, 3, 3), fun = function(x) {
  x[which.min(is.na(x))]  # Return the value of the nearest neighbor
})


writeRaster(t_raster, "data/t_ave_raster.tif", overwrite=TRUE)


#------------------future t average---------------------------------

#load rasters and average across months
# List all tif files in the folder
bioclim_files <- list.files("data/bioclim/wc2.1_30s_tavg/", pattern = "\\.tif$", full.names = TRUE)

# Initialize an empty list to store rasters
raster_list <- list()

# Loop through each HDF file
for (file in bioclim_files) {
  # Read the raster from the HDF file
  raster <- rast(file)
  
  # If there are multiple layers, select the layer you want or just take the first one
  raster_list <- append(raster_list, list(raster))
}



# Stack all rasters together (if they have the same extent and resolution)
raster_stack <- rast(raster_list)
land_bio1 <- mean(raster_stack, na.rm = TRUE)

land_t<-rast("data/bioclim/wc2.1_30s_bioc_BCC-CSM2-MR_ssp585_2081-2100.tif")
# Compute the average across layers
BCC585_t_bio1 <- land_t$`wc2.1_30s_bioc_BCC-CSM2-MR_ssp585_2081-2100_1`
BCC585_t_bio7 <- land_t$`wc2.1_30s_bioc_BCC-CSM2-MR_ssp585_2081-2100_7`
land_t1 <-rast("data/bioclim/wc2.1_30s_bioc_MIROC6_ssp585_2081-2100.tif")
# Define the function to convert Fahrenheit to Celsius
fahrenheit_to_celsius <- function(temp_f) {
  (temp_f - 32) * 5/9
}

# Apply the function to the SpatRaster object (land_t1)
MICRO585_bio1<- land_t1$`wc2.1_30s_bioc_MIROC6_ssp585_2081-2100_1`
MICRO585_bio7 <- land_t1$`wc2.1_30s_bioc_MIROC6_ssp585_2081-2100_7`


land_t2 <-rast("data/bioclim/wc2.1_30s_bioc_MRI-ESM2-0_ssp585_2081-2100.tif")
# Compute the average across layers
MRI585_bio1 <- land_t2$`wc2.1_30s_bioc_MRI-ESM2-0_ssp585_2081-2100_1`
MRI585_bio7<-land_t2$`wc2.1_30s_bioc_MRI-ESM2-0_ssp585_2081-2100_7`

land_t3 <-rast("data/bioclim/wc2.1_30s_bioc_BCC-CSM2-MR_ssp245_2081-2100.tif")
# Compute the average across layers
BCC245_bio1 <- land_t3$`wc2.1_30s_bioc_BCC-CSM2-MR_ssp245_2081-2100_1`
BCC245_bio7<-land_t3$`wc2.1_30s_bioc_BCC-CSM2-MR_ssp245_2081-2100_7`

land_t4 <-rast("data/bioclim/wc2.1_30s_bioc_MRI-ESM2-0_ssp245_2081-2100.tif")
# Compute the average across layers
MRI245_bio1 <- land_t4$`wc2.1_30s_bioc_MRI-ESM2-0_ssp245_2081-2100_1`
MRI245_bio7<-land_t4$`wc2.1_30s_bioc_MRI-ESM2-0_ssp245_2081-2100_7`

ocean_t <-rast("data/bioclim/tas_ssp585_2020_2100_depthsurf_0bec_f749_1d2f_U1734391577148.nc")
ocean_t1 <- rast("data/bioclim/tas_ssp585_2020_2100_depthsurf_c6fd_e300_ff0c_U1734391581734.nc")
ocean_t2 <-rast("data/bioclim/tas_ssp245_2020_2100_depthsurf_0bec_f749_1d2f_U1735770179064.nc")
ocean_t3 <- rast("data/bioclim/tas_ssp245_2020_2100_depthsurf_c6fd_e300_ff0c_U1735770169897.nc")



land_t585 <- c(MICRO585_bio1, MRI585_bio1)

bio1_585<-app(land_t585, fun=mean, rm.na=TRUE)
land_t245 <- c( BCC245_bio1, MRI245_bio1)

bio1_245<-app(land_t245, fun=mean, rm.na=TRUE)
bio1_245_resampled <- resample(bio1_245, ocean_t2)

bio1_585<-project(bio1_585, crs(elev)) 

bio1_585_resampled <- resample(bio1_585, ocean_t)

land585_bio7 <- c(MICRO585_bio7, MRI585_bio7)

bio7_585<-app(land585_bio7, fun=mean, rm.na=TRUE)

land245_bio7 <- c( BCC245_bio7, MRI245_bio7)

bio7_245<-app(land_t245, fun=mean, rm.na=TRUE)

bio7_585<-project(bio7_585, crs(elev)) 

bio7_585_resampled <- resample(bio7_585, ocean_t1)

bio7_245_resampled <- resample(bio7_245, ocean_t3)
# Crop the rasters
bio1_585_resampled<- crop(bio1_585_resampled, gcoos_extent)
bio7_585_resampled<- crop(bio7_585_resampled, gcoos_extent)
ocean_t <- crop(ocean_t, gcoos_extent)
ocean_t<-project(ocean_t, crs(elev))

ocean_t1 <- crop(ocean_t1, gcoos_extent)
ocean_t1<-project(ocean_t1, crs(elev))
ocean_t1 <-   273.15+ocean_t1

ocean_t2 <- crop(ocean_t2, gcoos_extent)
ocean_t2<-project(ocean_t2, crs(elev))

ocean_t3 <- crop(ocean_t3, gcoos_extent)
ocean_t3<-project(ocean_t3, crs(elev))
ocean_t3 <-   273.15+ocean_t3
bio1_245_resampled<- crop(bio1_245_resampled, gcoos_extent)
bio7_245_resampled<- crop(bio7_245_resampled, gcoos_extent)

# Step 1: Compute Mean and Standard Deviation for LCC_resampled_cropped
min_land_bio1 <- global(bio1_585_resampled, fun = "min", na.rm = TRUE)[1]
max_land_bio1 <- global(bio1_585_resampled, fun = "max", na.rm = TRUE)[1]
min_land_bio1_value <- as.numeric(min_land_bio1)
max_land_bio1_value <- as.numeric(max_land_bio1)

land_bio1_normalized <- app(bio1_585_resampled, fun = function(x) {
  (x - min_land_bio1_value) / (max_land_bio1_value - min_land_bio1_value)
})

# Step 1: Compute Mean and Standard Deviation for LCC_resampled_cropped
min_ocean_t <- global(ocean_t, fun = "min", na.rm = TRUE)[1]
max_ocean_t <- global(ocean_t, fun = "max", na.rm = TRUE)[1]
min_ocean_t_value <- as.numeric(min_ocean_t)
max_ocean_t_value <- as.numeric(max_ocean_t)

ocean_t_normalized <- app(ocean_t, fun = function(x) {
  (x - min_ocean_t_value) / (max_ocean_t_value - min_ocean_t_value)
})



bio1_raster <- app(
  c(land_bio1_normalized, ocean_t_normalized),
  fun = function(x) max(x, na.rm = TRUE),
  filename = "t_raster.tif",
  overwrite = TRUE
)

#fill the no-value cells
bio1_raster <- app(bio1_raster, fun = function(x) {
  x[is.infinite(x)] <- NA  # Replace -Inf with NA
  return(x)
})


# fill NA
bio1_raster<- focal(bio1_raster, w = matrix(1, 3, 3), fun = function(x) {
  x[which.min(is.na(x))]  # Return the value of the nearest neighbor
})


writeRaster(bio1_raster, "data/bio1_585_raster.tif", overwrite=TRUE)

# Step 1: Compute Mean and Standard Deviation for LCC_resampled_cropped
min_land_bio1 <- global(bio1_245_resampled, fun = "min", na.rm = TRUE)[1]
max_land_bio1 <- global(bio1_245_resampled, fun = "max", na.rm = TRUE)[1]
min_land_bio1_value <- as.numeric(min_land_bio1)
max_land_bio1_value <- as.numeric(max_land_bio1)

land_bio1_normalized <- app(bio1_245_resampled, fun = function(x) {
  (x - min_land_bio1_value) / (max_land_bio1_value - min_land_bio1_value)
})

# Step 1: Compute Mean and Standard Deviation for LCC_resampled_cropped
min_ocean_t <- global(ocean_t2, fun = "min", na.rm = TRUE)[1]
max_ocean_t <- global(ocean_t2, fun = "max", na.rm = TRUE)[1]
min_ocean_t_value <- as.numeric(min_ocean_t)
max_ocean_t_value <- as.numeric(max_ocean_t)

ocean_t_normalized <- app(ocean_t2, fun = function(x) {
  (x - min_ocean_t_value) / (max_ocean_t_value - min_ocean_t_value)
})



bio1_raster <- app(
  c(land_bio1_normalized, ocean_t_normalized),
  fun = function(x) max(x, na.rm = TRUE),
  filename = "t_raster.tif",
  overwrite = TRUE
)

#fill the no-value cells
bio1_raster <- app(bio1_raster, fun = function(x) {
  x[is.infinite(x)] <- NA  # Replace -Inf with NA
  return(x)
})


# fill NA
bio1_raster<- focal(bio1_raster, w = matrix(1, 3, 3), fun = function(x) {
  x[which.min(is.na(x))]  # Return the value of the nearest neighbor
})


writeRaster(bio1_raster, "data/bio1_245_raster.tif", overwrite=TRUE)


# Step 1: Compute Mean and Standard Deviation for LCC_resampled_cropped
min_land_bio7 <- global(bio7_585_resampled, fun = "min", na.rm = TRUE)[1]
max_land_bio7 <- global(bio7_585_resampled, fun = "max", na.rm = TRUE)[1]
min_land_bio7_value <- as.numeric(min_land_bio7)
max_land_bio7_value <- as.numeric(max_land_bio7)

land_bio7_normalized <- app(bio7_585_resampled, fun = function(x) {
  (x - min_land_bio7_value) / (max_land_bio7_value - min_land_bio7_value)
})

# Step 1: Compute Mean and Standard Deviation for LCC_resampled_cropped
min_ocean_t1 <- global(ocean_t1, fun = "min", na.rm = TRUE)[1]
max_ocean_t1 <- global(ocean_t1, fun = "max", na.rm = TRUE)[1]
min_ocean_t1_value <- as.numeric(min_ocean_t1)
max_ocean_t1_value <- as.numeric(max_ocean_t1)

ocean_t1_normalized <- app(ocean_t1, fun = function(x) {
  (x - min_ocean_t1_value) / (max_ocean_t1_value - min_ocean_t1_value)
})



bio7_raster <- app(
  c(land_bio7_normalized, ocean_t1_normalized),
  fun = function(x) max(x, na.rm = TRUE),
  filename = "t_raster.tif",
  overwrite = TRUE
)

#fill the no-value cells
bio7_raster <- app(bio7_raster, fun = function(x) {
  x[is.infinite(x)] <- NA  # Replace -Inf with NA
  return(x)
})


# fill NA
plot(bio7_raster)<- focal(bio7_raster, w = matrix(1, 3, 3), fun = function(x) {
  x[which.min(is.na(x))]  # Return the value of the nearest neighbor
})


writeRaster(bio7_raster, "data/bio7_585_raster.tif", overwrite=TRUE)


# Step 1: Compute Mean and Standard Deviation for LCC_resampled_cropped
min_land_bio7 <- global(bio7_245_resampled, fun = "min", na.rm = TRUE)[1]
max_land_bio7 <- global(bio7_245_resampled, fun = "max", na.rm = TRUE)[1]
min_land_bio7_value <- as.numeric(min_land_bio7)
max_land_bio7_value <- as.numeric(max_land_bio7)

land_bio7_normalized <- app(bio7_245_resampled, fun = function(x) {
  (x - min_land_bio7_value) / (max_land_bio7_value - min_land_bio7_value)
})

# Step 1: Compute Mean and Standard Deviation for LCC_resampled_cropped
min_ocean_t1 <- global(ocean_t3, fun = "min", na.rm = TRUE)[1]
max_ocean_t1 <- global(ocean_t3, fun = "max", na.rm = TRUE)[1]
min_ocean_t1_value <- as.numeric(min_ocean_t1)
max_ocean_t1_value <- as.numeric(max_ocean_t1)

ocean_t1_normalized <- app(ocean_t3, fun = function(x) {
  (x - min_ocean_t1_value) / (max_ocean_t1_value - min_ocean_t1_value)
})



bio7_raster <- app(
  c(land_bio7_normalized, ocean_t1_normalized),
  fun = function(x) max(x, na.rm = TRUE),
  filename = "t_raster.tif",
  overwrite = TRUE
)

#fill the no-value cells
bio7_raster <- app(bio7_raster, fun = function(x) {
  x[is.infinite(x)] <- NA  # Replace -Inf with NA
  return(x)
})


# fill NA
bio7_raster<- focal(bio7_raster, w = matrix(1, 3, 3), fun = function(x) {
  x[which.min(is.na(x))]  # Return the value of the nearest neighbor
})


writeRaster(bio7_raster, "data/bio7_245_raster.tif", overwrite=TRUE)


#load rasters and average across months
# List all tif files in the folder
mint_files <- list.files("data/bioclim/wc2.1_30s_tmin/", pattern = "\\.tif$", full.names = TRUE)

# Initialize an empty list to store rasters
raster_list <- list()

# Loop through each HDF file
for (file in mint_files) {
  # Read the raster from the HDF file
  raster <- rast(file)
  
  # If there are multiple layers, select the layer you want or just take the first one
  raster_list <- append(raster_list, list(raster))
}



# Stack all rasters together (if they have the same extent and resolution)
raster_stack <- rast(raster_list)
mint <- mean(raster_stack, na.rm = TRUE)

#load rasters and average across months
# List all tif files in the folder
maxt_files <- list.files("data/bioclim/wc2.1_30s_tmax/", pattern = "\\.tif$", full.names = TRUE)

# Initialize an empty list to store rasters
raster_list <- list()

# Loop through each HDF file
for (file in maxt_files) {
  # Read the raster from the HDF file
  raster <- rast(file)
  
  # If there are multiple layers, select the layer you want or just take the first one
  raster_list <- append(raster_list, list(raster))
}



# Stack all rasters together (if they have the same extent and resolution)
raster_stack <- rast(raster_list)
maxt <- mean(raster_stack, na.rm = TRUE)

trange<-maxt-mint

ocean_t_range <-rast("data/bioclim/tas_baseline_2000_2020_depthsurf_7db6_0749_5630_U1734473815048.nc")
ocean_t_range <-   273.15+ocean_t_range


trange<-project(trange, crs(elev)) 
ocean_t_range<-project(ocean_t_range, crs(elev))
trange_resampled <- resample(trange, ocean_t_range)

# Crop the rasters
trange_resampled<- crop(land_bio1_resampled, gcoos_extent)
ocean_t_range <- crop(ocean_t_range, gcoos_extent)


# Step 1: Compute Mean and Standard Deviation for LCC_resampled_cropped
min_trange <- global(trange_resampled, fun = "min", na.rm = TRUE)[1]
max_trange <- global(trange_resampled, fun = "max", na.rm = TRUE)[1]
min_trange_value <- as.numeric(min_trange)
max_trange_value <- as.numeric(max_trange)

trange_normalized <- app(trange_resampled, fun = function(x) {
  (x - min_trange_value) / (max_trange_value - min_trange_value)
})

# Step 1: Compute Mean and Standard Deviation for LCC_resampled_cropped
min_ocean_t_range <- global(ocean_t_range, fun = "min", na.rm = TRUE)[1]
max_ocean_t_range <- global(ocean_t_range, fun = "max", na.rm = TRUE)[1]
min_ocean_trange_value <- as.numeric(min_ocean_t_range)
max_ocean_trange_value <- as.numeric(max_ocean_t_range)

ocean_trange_normalized <- app(ocean_t_range, fun = function(x) {
  (x - min_ocean_trange_value) / (max_ocean_trange_value - min_ocean_trange_value)
})



trange_raster <- app(
  c(trange_normalized, ocean_trange_normalized),
  fun = function(x) max(x, na.rm = TRUE),
  filename = "t_raster.tif",
  overwrite = TRUE
)

#fill the no-value cells
trange_raster <- app(trange_raster, fun = function(x) {
  x[is.infinite(x)] <- NA  # Replace -Inf with NA
  return(x)
})


# fill NA
trange_raster<- focal(trange_raster, w = matrix(1, 3, 3), fun = function(x) {
  x[which.min(is.na(x))]  # Return the value of the nearest neighbor
})

plot(trange_raster)

writeRaster(trange_raster, "data/t_range_raster.tif", overwrite=TRUE)


# # Set CRS, extent, and resolution for each raster in the stack
# raster_stack <- lapply(raster_stack, function(r) {
#   crs(r) <- "EPSG:4326"  # Set CRS
#   ext(r) <- c(-180, 180, -90, 90)  # Global extent
#   res(r) <- c(0.0833, 0.0833)  # Set rsolution for latitude and longitude
#   return(r)
# })
# 
# 
# # 
# # Convert back to SpatRaster stack
# raster_stack <- rast(raster_stack)
# 
# # Crop and project each raster and then stack them together
# processed_stack <- lapply(raster_stack, function(r) {
#   r <- crop(r, gcoos_extent)  # Crop based on gcoos_extent
#   r <- project(r, crs(elev))  # Project raster to match the CRS of 'elev'
#   return(r)
# })
# # 
# # Convert list of rasters back into a single stack
# processed_stack <- rast(processed_stack)
# 
# # Compute the mean raster across the stack
# average_raster <- mean(processed_stack, na.rm = TRUE)
# plot(npp_mean)<-average_raster


npp_ocean<-rast()
library(terra)

# Reproject the raster to EPSG:4326
plot(npp_4326) <- project(npp, "EPSG:4326")







# load ebird data
ebird <- read_csv(paste0("data/ebird/all_species/",  "ebd_gom", "_bbox_zf.csv"))
length(levels(as.factor(ebird$scientific_name)))
# buffer to create neighborhood (6 km by 6 km or 24km by 24km) around each checklist point 
library(dplyr)

# Filter rows where latitude >= 31.16
filtered_ebird <- ebird %>%
  filter(latitude <= 31.16)


scientific_to_common <-read.csv("data/ebird-taxonomy.csv")
scientific_to_common <- scientific_to_common %>% select(scientific_name, common_name)
library(dplyr)
library(stringr)


# Create a function to sanitize the scientific name for use as a file name
sanitize_filename <- function(name) {
  # Replace spaces with underscores and remove non-alphanumeric characters
  name <- gsub("[^[:alnum:][:space:]]", "", name)  # Remove special characters
  name <- gsub(" ", "_", name)  # Replace spaces with underscores
  return(name)
}

# Loop over each unique scientific name
unique_scientific_name <- unique(filtered_ebird$scientific_name)




# Define batch size
batch_size <- 10
data_with_common<-NULL
# Get the total number of unique scientific names
total_scientific_names <- length(unique_scientific_names)

# Process in batches
for (i in seq(1, total_scientific_names, by = batch_size)) {
  # Get the current batch of scientific names
  current_batch <- unique_scientific_names[i:min(i + batch_size - 1, total_scientific_names)]
  
  # Filter and join for the entire batch
  data_batch <- filtered_ebird %>%
    filter(scientific_name %in% current_batch) %>%
    left_join(scientific_to_common, by = "scientific_name")
  
  # Bind the processed batch to the result
  data_with_common$scientific_name <- bind_rows(data_with_common, data_batch)
}

#data_with_common["common_name"][is.na(data_with_common["common_name"])]<- "Cooper's Hawk"
    # Save the subset to a CSV file


t(data_with_common)
# Pivot data
ebird_with_commonname <- data_with_common %>%
  select(locality_id, latitude, longitude, common_name, observation_count, species_observed, protocol_type, day_of_year, time_observations_started,
         duration_minutes, effort_distance_km, number_observers) 
write.csv(ebird_with_commonname, file = "data/ebird/all_species/ebird_with_commonname.csv", row.names = FALSE)

ebird_with_commonname <- read_csv("data/ebird/all_species/cleaned_ebird_data.csv")

ebird_count <- ebird_with_commonname %>% 
  mutate(protocol_type = factor(protocol_type, 
                                levels = c("Stationary" , "Traveling"))) %>%
  # remove observations with no count
  filter(!is.na(observation_count))
# habitat covariates with crops categorized 
habitat <- read_csv("data/esrd/pland_elev_north_east_slope_ocean_chl_poc_t_trange_3km.csv") %>% 
  mutate(year = as.integer(year))

ebird_habitat <-inner_join(ebird_count, habitat, by = c("locality_id")) 

ebird_habitat_new <- bind_cols(ebird_habitat[,1:27], ebird_habitat[,30:40])

write.csv(ebird_habitat_new, file = "data/ebird/all_species/ebird_habitat.csv", row.names = FALSE)

ebird_habitat<-read_csv("data/ebird/all_species/ebird_habitat_float.csv")
library(tidyr)

# Assuming your data frame is named df
pivoted_df <- ebird_habitat %>%
  pivot_wider(
    names_from = common_name, 
    values_from = observation_count,
    values_fill = list(observation_count = 0)
  )

# View the pivoted data frame
print(pivoted_df)




ebird_buff <- ebird_with_commonname %>% 
  distinct(locality_id, latitude, longitude) %>% 
  # match the ebird year with corresponding year landcover data
  mutate(year = 2020, year_lc = paste0("y", "2020"))  %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  # transform to cdl projection
  st_transform(crs = crs(elev)) %>% 
  st_buffer(dist = neighborhood_radius) %>% 
  # nest by year to  match to land cover data from the corresponding year
  nest(data = c(year, geometry))

ccap_classes <- read.csv(file = "data/ccap_landcover_classes.csv")



library(dplyr)
library(tidyr)
library(sf)

# # Initialize parameters
# chunk_size <- 1000  # Define the number of rows per chunk
# lc_extract_crs <- NULL  # Initialize an empty object to store results
# 
# # Iterate over unique years in the landcover data
# for (yr in unique(ebird_buff$year_lc)) {
#   # Filter rows for the current year
#   year_data <- ebird_buff %>%
#     filter(year_lc == yr)
#   # Get the buffered checklists for the given year
#   regions <- bind_rows(lapply(year_data$data, st_as_sf))
#   
#   # Process data in chunks
#   total_rows <- nrow(regions)
#   
#   for (start in seq(1, total_rows, by = chunk_size)) {
#     # Extract a chunk of rows
#     chunk <- regions %>% slice(start:(start + chunk_size - 1))
#     
#     # Extract landcover values within the buffered checklist areas
#     ee <- extract(landcover[[yr]], st_transform(chunk, crs(landcover[[yr]])), progress = FALSE)
#     
#     # Count the number of each landcover class for each checklist buffer
#     ee_count <- ee %>%
#       as_tibble() %>%
#       group_by(ID, landcover = !!sym(yr), drop = FALSE) %>%
#       summarise(count = n(), .groups = 'drop')
#     # Convert extracted values into a tibble and replace NA landcover values with 26
#     
#     # Calculate PLAND: the proportion of each land cover class within the neighborhood
#     pland <- ee_count %>%
#       group_by(ID) %>%
#       mutate(pland = count / sum(count)) %>%
#       ungroup() %>%
#       select(-count)
#     
#     # Add landcover class names
#     pland_ccap <- inner_join(pland, ccap_classes, by = "landcover")
#     
#     # Transform to wide format, filling implicit missing values with 0
#     pland_wide <- pland_ccap %>%
#       select(-landcover) %>%
#       group_by(ID, landclass) %>%
#       summarise(pland = sum(pland), .groups = 'drop') %>%
#       pivot_wider(
#         names_from = landclass,
#         values_from = pland,
#         values_fill = 0
#       )
#     
#     # Attach the year and locality ID back to the checklists
#     pland_summ <- tibble(st_drop_geometry(chunk), data = pland_wide)
#     
#     # Append the processed chunk to the results
#     lc_extract_crs <- bind_rows(lc_extract_crs, pland_summ)
#   }
# }


# iterate over all years extracting landcover for all checklists in each
lc_extract <- NULL
for (yr in "y2020") {
  # get the buffered checklists for a given year
  # Filter rows for the current year
  year_data <- ebird_buff %>%
    filter(year_lc == yr)
  # Get the buffered checklists for the given year
  regions <- bind_rows(lapply(year_data$data, st_as_sf))
  coords_year <- regions %>% st_as_sf()%>%st_transform(crs = 4326) %>% st_centroid() %>%
    st_coordinates() %>% 
    as.data.frame() %>% 
    rename(longitude = X, latitude = Y)  
  # get landcover values within each buffered checklist area
  ee <- extract(landcover[[yr]], regions, progress = FALSE)
  # count the number of each landcover class for each checklist buffer
  ee_count <- ee %>%
    as_tibble() %>%
    group_by(ID, landcover = !!sym(yr), .drop = FALSE) %>%
    summarise(count = n(), .groups = 'drop')
  
  #calculate PLAND: the proportion of each land cover class within the neighborhood.
  pland <- ee_count %>% 
    # calculate proporiton
    group_by(ID) %>% 
    mutate(pland = count / sum(count)) %>% 
    ungroup() %>% 
    dplyr::select(-count) 
  pland_ccap <- inner_join(pland, ccap_classes,  by = "landcover")
  # tranform to wide format, filling in implicit missing values with 0s%>% 
  pland_wide <- pland_ccap %>% 
    select(-landcover) %>% 
    group_by(ID, landclass)%>%
    summarize(pland=sum(pland)) %>%
    pivot_wider(names_from = landclass, 
                values_from = pland, 
                values_fill =  0)
  # attach the year and locality id back to the checklists
  pland_sum <- tibble(st_drop_geometry(regions), coords_year) 
  
  pland_sumn <- bind_cols(ebird_buff[, 'locality_id'], pland_sum, pland_wide)
  pland_sumn <-
  lc_extract <- bind_row(lc_extract, pland_sumn)
}


na_count_B <- sum(is.na(lc_extract$ocean))
print(paste("Column B has", na_count_B, "NA values."))


# coords_ebird <- ebird_buff$data[[1]] %>% st_as_sf()%>%st_transform(crs = 4326) %>% st_centroid() %>%
#   st_coordinates() %>% 
#   as.data.frame() %>% 
#   rename(longitude = X, latitude = Y)  
# 
# # Final result is stored in `lc_extract`
# lc_extract_crs<-lc_extract_crs %>% 
#   unnest(, cols = data) %>%
#   select(-ID) %>% 
#   bind_cols(coords_ebird)



# Replace non-1 values with 0
lc_extract$ocean_fact <- ifelse(lc_extract$ocean == 1, 1, 0)

# Convert the column to a factor (categorical)
lc_extract$ocean_fact <- as.factor(lc_extract$ocean_fact)


# Initialize an empty tibble to store results
all_locs <- tibble()

# Define the chunk size
chunk_size <- 100  # Adjust based on your memory and data size

# Iterate over unique years in the landcover data
for (yr in unique(ebird_buff$year_lc)) {
  # Filter rows for the current year
  year_data <- ebird_buff %>%
    filter(year_lc == yr)
  # Get the buffered checklists for the given year
  regions <- bind_rows(lapply(year_data$data, st_as_sf))
  
  if (is.null(regions) || nrow(regions) == 0) next  # Skip if no data for this year
  
  # Process data in chunks
  total_rows <- nrow(regions)
  
  for (start in seq(1, total_rows, by = chunk_size)) {
    # Define the slicing range
    end <- min(start + chunk_size - 1, total_rows)
    
    # Extract a chunk of rows
    chunk <- regions %>% slice(start:end)
    
    # Extract locality information without geometry
    locs_cell <- st_set_geometry(chunk, NULL)
    
    # Extract elevation values and calculate median and sd
    elev_cell <- extract(elev, chunk, progress = FALSE) %>%
      as_tibble() %>%
      group_by(ID) %>%
      summarise(
        elevation_median = mean(Band1, na.rm = TRUE),
        elevation_sd = sd(Band1, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      # Combine with locality info
      bind_cols(locs_cell) %>%
      select(-ID)
    
    # Append chunk results to the combined dataset
    all_locs <- bind_rows(all_locs, elev_cell)
  }
  
  # Force garbage collection to free up memory
  rm(regions)
  gc()
}

# Perform the final inner join after all rows are bound
pland_elev <- bind_cols(lc_extract, all_locs[,1:2])

# Extract the values
elevation_values <- pland_elev$elevation_median

# Perform min-max normalization
elevation_normalized <- (elevation_values - min(elevation_values, na.rm = TRUE)) / 
  (max(elevation_values, na.rm = TRUE) - min(elevation_values, na.rm = TRUE))

# Add the normalized values back to the data frame
pland_elev$elevation_median <- elevation_normalized

# Extract the values
elevation_values <- pland_elev$elevation_sd

# Perform min-max normalization
elevation_normalized <- (elevation_values - min(elevation_values, na.rm = TRUE)) / 
  (max(elevation_values, na.rm = TRUE) - min(elevation_values, na.rm = TRUE))

# Add the normalized values back to the data frame
pland_elev$elevation_sd <- elevation_normalized


# Save or inspect results
write.csv(pland_elev, "data/esrd/pland_elev_all_species_ocean.csv", row.names = FALSE)




nort_locs <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  # Filter rows for the current year
  year_data <- ebird_buff %>%
    filter(year_lc == yr)
  # Get the buffered checklists for the given year
  regions <- bind_rows(lapply(year_data$data, st_as_sf))
  locs_cell <- st_set_geometry(regions, NULL) 
  nort_cell <- extract(nort, regions, progress = FALSE) %>%
    map_dfr(~ tibble(northness_median = mean(., na.rm = TRUE),
                     northness_sd = sd(., na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_cell, .)  
    
  # bind to results
  nort_locs <- bind_rows(nort_locs, nort_cell)
}


pland_elev_north <- bind_cols(pland_elev, nort_locs[,2:3])

east_locs <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  # Filter rows for the current year
  year_data <- ebird_buff %>%
    filter(year_lc == yr)
  # Get the buffered checklists for the given year
  regions <- bind_rows(lapply(year_data$data, st_as_sf))
  locs_cell <- st_set_geometry(regions, NULL) 
  east_cell <- extract(east, regions, progress = FALSE) %>%
    map_dfr(~ tibble(eastness_median = mean(., na.rm = TRUE),
                     eastness_sd = sd(., na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_cell, .) 
  # bind to results
  east_locs <- bind_rows(east_locs, east_cell)
}
# prediction surface covariates
pland_elev_north_east <- bind_cols(pland_elev_north, east_locs[,2:3])


slope_locs <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  # Filter rows for the current year
  year_data <- ebird_buff %>%
    filter(year_lc == yr)
  # Get the buffered checklists for the given year
  regions <- bind_rows(lapply(year_data$data, st_as_sf))
  locs_cell <- st_set_geometry(regions, NULL) 
  slope_cell <- extract(slope, regions, progress = FALSE) %>%
    map_dfr(~ tibble(slopeness_median = mean(., na.rm = TRUE),
                     slopeness_sd = sd(., na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_cell, .) 
  # bind to results
  slope_locs <- bind_rows(slope_locs, slope_cell)
}
# prediction surface covariates
pland_elev_north_east_slope <- bind_cols(pland_elev_north_east, slope_locs[,2:3])

# save

write_csv(pland_elev_north_east_slope, paste0("data/esrd/", "pland_elev_north_east_slope_ocean_3km.csv"))

#repeat the same process for northness covariates
poc <- raster("data/normalized_poc_raster.tif")
poc <-  poc %>% 
  crop(., gcoos_extent) %>% 
  projectRaster(crs = crs(elev)) 

chl <- raster("normalized_chl_raster.tif")
chl <-  chl %>% 
  crop(., gcoos_extent) %>% 
  projectRaster(crs = crs(elev)) 

chl_locs <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  # Filter rows for the current year
  year_data <- ebird_buff %>%
    filter(year_lc == yr)
  # Get the buffered checklists for the given year
  regions <- bind_rows(lapply(year_data$data, st_as_sf))
  locs_cell <- st_set_geometry(regions, NULL) 
  chl_cell <- extract(chl, regions, progress = FALSE) %>%
    map_dfr(~ tibble(chl_mean = mean(., na.rm = TRUE),
                     chl_sd = sd(., na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_cell, .) 
  # bind to results
  chl_locs <- bind_rows(chl_locs, chl_cell)
}



# prediction surface covariates
pland_elev_north_east_slope_chl <- bind_cols(pland_elev_north_east_slope, chl_locs[,2])

poc_locs <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  # Filter rows for the current year
  year_data <- ebird_buff %>%
    filter(year_lc == yr)
  # Get the buffered checklists for the given year
  regions <- bind_rows(lapply(year_data$data, st_as_sf))
  locs_cell <- st_set_geometry(regions, NULL) 
  poc_cell <- extract(poc, regions, progress = FALSE) %>%
    map_dfr(~ tibble(poc_mean = mean(., na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_cell, .) 
  # bind to results
  poc_locs <- bind_rows(poc_locs, poc_cell)
}


# prediction surface covariates
pland_elev_north_east_slope_chl_poc <- bind_cols(pland_elev_north_east_slope_chl, poc_locs[,2])

# save

write_csv(pland_elev_north_east_slope_chl_poc, paste0("data/esrd/", "pland_elev_north_east_slope_ocean_chl_poc_3km.csv"))


t_locs <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  # Filter rows for the current year
  year_data <- ebird_buff %>%
    filter(year_lc == yr)
  # Get the buffered checklists for the given year
  regions <- bind_rows(lapply(year_data$data, st_as_sf))
  locs_cell <- st_set_geometry(regions, NULL) 
  t_cell <- extract(temperature, regions, progress = FALSE) %>%
    map_dfr(~ tibble(t_mean = mean(., na.rm = TRUE),
                     t_sd = sd(., na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_cell, .) 
  # bind to results
  t_locs <- bind_rows(t_locs, t_cell)
}

# prediction surface covariates
pland_elev_north_east_slope_chl_poc_t <- bind_cols(pland_elev_north_east_slope_chl_poc, t_locs[,2:3])

# save

write_csv(pland_elev_north_east_slope_chl_poc_t, paste0("data/esrd/", "pland_elev_north_east_slope_ocean_chl_poc_t_3km.csv"))

trange_locs <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  # Filter rows for the current year
  year_data <- ebird_buff %>%
    filter(year_lc == yr)
  # Get the buffered checklists for the given year
  regions <- bind_rows(lapply(year_data$data, st_as_sf))
  locs_cell <- st_set_geometry(regions, NULL) 
  trange_cell <- extract(trange_raster, regions, progress = FALSE) %>%
    map_dfr(~ tibble(trange_mean = mean(., na.rm = TRUE),
                     trange_sd = sd(., na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_cell, .) 
  # bind to results
  trange_locs <- bind_rows(trange_locs, trange_cell)
}

# prediction surface covariates
pland_elev_north_east_slope_t_poc_chl <- bind_cols(pland_elev_north_east_slope, t_locs[2], trange_locs[2], poc_locs[2], chl_locs[2])
# save

write_csv(pland_elev_north_east_slope_t_poc_chl, paste0("data/esrd/", "pland_elev_north_east_slope_ocean_chl_poc_t_trange_3km.csv"))


#create prediction grid with 3km*3km resolution
pts <- st_sample(gcoos_extent, 100000) %>% 
  st_sf(as.data.frame(st_coordinates(.)), geometry = .) %>% #converts the sampled points into an sf object with a geometry column representing the spatial points. 
  rename(lat = Y, lon = X)
dggs <- dgconstruct(spacing=3)
pts$cell <- dgGEO_to_SEQNUM(dggs, pts$lon, pts$lat)$seqnum #his assigns a cell identifier to each point in the pts dataset. 
#It uses the dgGEO_to_SEQNUM() function to convert the latitude and longitude coordinates of the points to cell identifiers (seqnum) based on the DGGS grid.
hexagon <- dgcellstogrid(dggs, unique(pts$cell)) %>%  #unique(pts$cell): This extracts the unique cell identifiers from the pts dataset.
  st_as_sf() #This converts the hexagonal grid polygons to an sf object.
#The resulting pts and hexagon objects will contain the sampled points and the hexagonal grid polygons
nrow(hexagon)
ggplot() +
  geom_sf(data = hexagon, colour = "red", fill = NA ) +
  theme_bw()

library(terra)
library(sf)

# Plot the SpatRaster (elevation data)
plot(elev, main = "Hexagon Boundaries with Elevation Raster", col = terrain.colors(100))
ext(elev)
extent(hexagon)
# Overlay hexagon boundaries
plot(st_geometry(hexagon), add = TRUE, col = NA, border = "red", lwd = 0.5)

hexagon_multipolygon <- st_cast(hexagon, "MULTIPOLYGON") %>% st_transform(crs=4326)
# Save to a GeoJSON file
st_write(hexagon, "data/hexagon.geojson", append=FALSE)




##cannot use this becuase the crs not desired
# # with a radius equal to 200 or 800 c-cap cells 
# neighborhood_radius <- 100 * ceiling(max(res(landcover))) / 2
# agg_factor <- round(2 * neighborhood_radius / res(landcover))


# r <- raster(landcover_cropped) %>% 
#   aggregate(agg_factor) 
# 
# r <- gcoos %>% 
#   st_transform(crs = projection(r)) %>% 
#   rasterize(r, field = 1) %>% 
#   # remove any empty cells at edges
#   trim()

#writeRaster(r, filename = "data/prediction-surface-drnet.tif", overwrite = TRUE)


# buffer to create neighborhood (6 km by 6 km) around each checklist point 
# with a radius equal to 100/2 ccap cells 1500m
neighborhood_radius <- 100 * ceiling(max(res(landcover))) / 2

hexagon<-st_read("data/hexagon_multipolygon.geojson")
cell_buff <- NULL
for (yr in "2020") {
  r_centers <- st_centroid(hexagon) %>% 
    st_as_sf() %>% 
    transmute(id = row_number()) %>% 
    mutate(year=as.character(yr), year_lc = paste0("y", yr))
  r_cells <- st_buffer(r_centers, dist = neighborhood_radius)
  
  cell_buff <- bind_rows(cell_buff, r_cells)
  # nest by year to  match to land cover data from the corresponding year
}
#cell_buff <- nest(cell_buff, data = c(year, id, geometry))

# ------------------------------------------
# add elevation to the surface
# ------------------------------------------


locs_pred <- NULL
na.omit(locs_pred)

for (yr in "y2020") {
  regions <- cell_buff
  
  # Process in smaller chunks
  for (chunk in split(regions, ceiling(seq_along(regions$id) / 1000))) {
    locs_cell <- st_set_geometry(chunk, NULL)
    elev_cell <- extract(elev, chunk, progress = FALSE) %>%
      as_tibble() %>%
      group_by(ID) %>%
      summarise(
        elevation_median = median(Band1, na.rm = TRUE),
        elevation_sd = sd(Band1, na.rm = TRUE)
      ) %>% 
      select(-ID)
    
    locs_pred <- bind_rows(locs_pred, bind_cols(locs_cell, elev_cell))
  }
}

# prediction surface covariates


locs_cleaned<-na.omit(locs_pred)

# Filter cell_buff to keep only rows with IDs present in valid_ids
cell_buff_short <- cell_buff %>% 
  filter(id %in% locs_cleaned$id)

cell_buff_short <- cell_buff_short %>% 
  transmute(id = row_number())


st_write(cell_buff_short, "data/cell_buff_short.shp", delete_layer = TRUE)
cell_buff_short <-st_read("data/cell_buff_short.shp")

# save


ccap_classes <- read.csv(file = "data/ccap_landcover_classes.csv")

# iterate over all years extracting landcover for all checklists in each
lc_extract_pred <- NULL
for (yr in "y2020") {
  # get the buffered checklists for a given year
  regions <- cell_buff_short
  coords_year <- regions %>% st_as_sf()%>%st_transform(crs = 4326) %>% st_centroid() %>%
    st_coordinates() %>% 
    as.data.frame() %>% 
    rename(longitude = X, latitude = Y)  
  # get landcover values within each buffered checklist area
  ee <- extract(landcover[[yr]], regions, progress = FALSE)
  # count the number of each landcover class for each checklist buffer
  ee_count <- ee %>%
    as_tibble() %>%
    group_by(ID, landcover = !!sym(yr), .drop = FALSE) %>%
    summarise(count = n(), .groups = 'drop')
  
  #calculate PLAND: the proportion of each land cover class within the neighborhood.
  pland_pred <- ee_count %>% 
    # calculate proporiton
    group_by(ID) %>% 
    mutate(pland = count / sum(count)) %>% 
    ungroup() %>% 
    dplyr::select(-count) 
  pland_pred_ccap <- inner_join(pland_pred, ccap_classes,  by = "landcover")
  # tranform to wide format, filling in implicit missing values with 0s%>% 
  pland_pred_wide <- pland_pred_ccap %>% 
    select(-landcover) %>% 
    group_by(ID, landclass)%>%
    summarize(pland=sum(pland)) %>%
    pivot_wider(names_from = landclass, 
                values_from = pland, 
                values_fill =  0)
  # attach the year and locality id back to the checklists
  pland_pred_sum <- tibble(st_drop_geometry(regions), coords_year) 
  
  pland_pred_sumn <- pland_pred_wide %>% rename(id = ID) %>% inner_join(pland_pred_sum, by="id")
  # bind to results
  lc_extract_pred <- bind_rows(lc_extract_pred, pland_pred_sumn)
}



# Replace non-1 values with 0
lc_extract_pred$ocean_fact <- ifelse(lc_extract_pred$ocean == 1, 1, 0)

# Convert the column to a factor (categorical)
lc_extract_pred$ocean_fact <- as.factor(lc_extract_pred$ocean_fact)

library(terra)

# Check extent of hexagon
hex_extent <- ext(gcoos)
print(hex_extent)  # Verify the xmin, xmax, ymin, ymax

# Set desired resolution in degrees (3km for EPSG:4326)
res <- 0.027

# Compute the number of rows and columns
ncol <- ceiling((hex_extent$xmax - hex_extent$xmin) / res)
nrow <- ceiling((hex_extent$ymax - hex_extent$ymin) / res)

# Create the raster template
r_template <- rast(nrows = nrow, ncols = ncol, 
                   xmin = hex_extent$xmin, xmax = hex_extent$xmax, 
                   ymin = hex_extent$ymin, ymax = hex_extent$ymax, 
                   crs = 4326)


# fill the raster with 1s inside the study region
r <- rasterize(gcoos, r_template) |> 
  setNames("study_region")

# save for later use
r <- writeRaster(r, "data/prediction-grid_gom.tif",
                 overwrite = TRUE,
                 gdal = "COMPRESS=DEFLATE")
# buffer to create neighborhood (6 km by 6 km) around each checklist point 
# with a radius equal to 100/2 CDL cells 
neighborhood_radius <- 1500


cell_buff <- NULL
for (yr in "2020") {
  r_centers <- terra::xyFromCell(r, 1:terra::ncell(r)) %>% 
    as.data.frame() %>%
    rename(longitude = x, latitude = y) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(r)$wkt) %>%
    transmute(id = row_number()) %>% 
    mutate(year = as.character(yr), year_lc = paste0("y", yr))
  r_cells <- st_buffer(r_centers, dist = neighborhood_radius)
  cell_buff <- bind_rows(cell_buff, r_cells)
  # nest by year to  match to land cover data from the corresponding year
}
cell_buff <- nest(cell_buff, data = c(year, id, geometry))


# Convert to spatial features
ocean_cover_sf <- lc_extract_pred |> 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Rasterize points
ocean_cover <- rasterize(ocean_cover_sf, r, field = "ocean", fun = mean, background = NA)



# Make a map
par(mar = c(0.25, 0.25, 2, 0.25))
plot(ocean_cover, 
     axes = FALSE, box = FALSE, col = viridis(10), 
     main = "Ocean (% cover)")

pland_elev_pred <- bind_cols(lc_extract_pred, locs_cleaned[4:5])

write_csv(pland_elev_pred, "data/pred_surface/pred_pland_elev_20192021_3km.csv")

# ------------------------------------------
# add northness to  pred surface
# ------------------------------------------

#repeat for prediction surface
nort_pred <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff_short
  locs_nort_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  nort_cell <- extract(nort, regions, progress = FALSE) %>% 
     map_dfr(~ tibble(northness_median = mean(., na.rm = TRUE),
                     northness_sd = sd(., na.rm = TRUE))) %>% 
    # join to lookup table to get id
  bind_cols(locs_nort_cell, .)
  # bind to results
  nort_pred <- bind_rows(nort_pred, nort_cell)
}




# prediction surface covariates

pland_elev_nort_pred <- bind_cols(pland_elev_pred, nort_pred[2:3])

# ------------------------------------------
# add eastness to checklists and pred surface
# ------------------------------------------

#repeat with eastness

east_pred <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff_short
  locs_east_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  east_cell <- extract(east, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(eastness_median = mean(., na.rm = TRUE),
                     eastness_sd = sd(., na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_east_cell, .)
  # bind to results
  east_pred <- bind_rows(east_pred, east_cell)
}



# prediction surface covariates


pland_elev_nort_east_pred <- bind_cols(pland_elev_nort_pred, east_pred[2:3])
write_csv(pland_elev_nort_east_pred, "data/pred_surface/pred_pland_elev_nort_east_20192021_3km.csv")


slope_pred <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff_short
  locs_slope_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  slope_cell <- extract(slope, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(slopeness_median = mean(., na.rm = TRUE),
                     slopeness_sd = sd(., na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_slope_cell, .) 
  # bind to results
  slope_pred <- bind_rows(slope_pred, slope_cell)
}

# prediction surface covariates
pland_elev_nort_east_slope_pred <- bind_cols(pland_elev_nort_east_pred, slope_pred[2:3])
write_csv(pland_elev_nort_east_slope_pred, "data/pred_surface/pred_pland_elev_nort_east_slope_20192021_3km.csv")



temperature_pred <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff_short
  locs_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  bio1_cell <- extract(temperature, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(t_mean = mean(., na.rm = TRUE)
    )) %>% 
    # join to lookup table to get id
    bind_cols(locs_cell, .) 
  # bind to results
  temperature_pred <- bind_rows(temperature_pred, bio1_cell)
}
pland_elev_nort_east_slope_t_pred <- bind_cols(pland_elev_nort_east_slope_pred, temperature_pred[2])
write_csv(pland_elev_nort_east_slope_t_pred, "data/pred_surface/pred_pland_elev_nort_east_slope_t_20192021_3km.csv")

trange_pred <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff_short
  locs_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  bio7_cell <- extract(trange, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(trange_mean = mean(., na.rm = TRUE)
    )) %>% 
    # join to lookup table to get id
    bind_cols(locs_cell, .) 
  # bind to results
  trange_pred <- bind_rows(trange_pred, bio7_cell)
}

poc_pred <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff_short
  locs_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  poc_cell <- extract(poc, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(poc_mean = mean(., na.rm = TRUE)
    )) %>% 
    # join to lookup table to get id
    bind_cols(locs_cell, .) 
  # bind to results
  poc_pred <- bind_rows(poc_pred, poc_cell)
}

chl_pred <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff_short
  locs_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  chl_cell <- extract(chl, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(chl_mean = mean(., na.rm = TRUE)
    )) %>% 
    # join to lookup table to get id
    bind_cols(locs_cell, .) 
  # bind to results
  chl_pred <- bind_rows(chl_pred, chl_cell)
}

pland_elev_nort_east_slope_t_poc_chl_pred <- bind_cols(pland_elev_nort_east_slope_t_pred, trange_pred[2], poc_pred[2], chl_pred[2])
write_csv(pland_elev_nort_east_slope_t_poc_chl_pred, "data/pred_surface/pred_pland_elev_nort_east_slope_t_poc_chl_20192021_3km.csv")


bio1_pred <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff_short
  locs_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  bio1_cell <- extract(bio1, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(bio1 = mean(., na.rm = TRUE)
                     )) %>% 
    # join to lookup table to get id
    bind_cols(locs_cell, .) 
  # bind to results
  bio1_pred <- bind_rows(bio1_pred, bio1_cell)
}
bio7_pred <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff_short
  locs_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  bio7_cell <- extract(bio7, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(bio7 = mean(., na.rm = TRUE)
    )) %>% 
    # join to lookup table to get id
    bind_cols(locs_cell, .) 
  # bind to results
  bio7_pred <- bind_rows(bio7_pred, bio7_cell)
}
# prediction surface covariates
pland_elev_nort_east_slope_t_pred <- bind_cols(pland_elev_nort_east_slope_pred, bio1_pred[2], bio7_pred[2])
write_csv(pland_elev_nort_east_slope_t_pred, "data/pred_surface/pred_pland_elev_nort_east_slope_t_20192021_3km.csv")

# ------------------------------------------
# modify elevation
# ------------------------------------------

elev_init <- rast("data/storm_surge/EESLR_Storm-Surge_1.2_MOM_Initial_Conditions-004/StormSurge_MOM_2100_Initial_WSE.tif")
#reproject to match the crs of all other variables
elev_init <- project(elev_init, elev)
elev_high <- rast("data/storm_surge/EESLR_Storm-Surge_1.1_MOM_High_2m-001/StormSurge_MOM_2100_High_WSE.tif")
#reproject to match the crs of all other variables

idag_init <- rast("data/storm_surge/EESLR_Storm-Surge_1.2_MOM_Initial_Conditions-004/StormSurge_MOM_2100_Initial_IDAG.tif")
#reproject to match the crs of all other variables
idag_init <- project(idag_init, elev)
plot(elev_high)
plot(idag_init)



idag_pred <- NULL
for (yr in "y2020") {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff$data[[which(yr == cell_buff$year_lc)]]
  locs_cell <- st_set_geometry(regions, NULL) 
  idag_cell <- extract(idag_init, regions, progress = FALSE)%>%
    as_tibble() %>%
    group_by(ID) %>%
    summarize(idag_max = mean(StormSurge_MOM_2100_Initial_IDAG, na.rm = TRUE),
              idag_sd = sd(StormSurge_MOM_2100_Initial_IDAG, na.rm = TRUE)) %>% 
    ungroup() %>%
    # join to lookup table to get id
    bind_cols(locs_cell) %>%
    select(-ID) 
  # bind to results
  idag_pred <- bind_rows(idag_pred, idag_cell)
}


locs_init <- locs_pred %>%
  left_join(
    idag_pred %>%
      replace_na(list(idag_max = 0)) %>%  # Replace NA in idag_max with 0
      select(id, idag_max), 
    by = "id"
  ) %>%
  mutate(elevation_median = elevation_median + idag_max) %>%
  select(-idag_max)



# prediction surface covariates
pland_elev_pred <- inner_join(lc_extract_pred, locs_init, by = c("id"))
write_csv(pland_elev_pred, "data/pred_surface/pred_pland_elev_init_20192021_gcoos_6km.csv")




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




#define diameter
neighborhood_radius <- 100 * ceiling(max(res(landcover))) / 2

# # min(pland_geometry$year)
# ebird_buff_pmp <- NULL
# for (cp in names(crop_change_raster)) {
#   ebird_buff_crop <- ebird %>%
#     distinct(crop = cp,
#              locality_id, latitude, longitude)
#   ebird_buff_pmp <- rbind(ebird_buff_pmp, ebird_buff_crop)
# }
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