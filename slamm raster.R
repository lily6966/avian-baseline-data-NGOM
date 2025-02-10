library(sf)
library(rnaturalearth)
library(dplyr)
library(terra)
library(raster)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

select <- dplyr::select

# CV boundary sf

wd <- ("/Users/liyingnceas/GitHub/gulf_of_mexico_slr_birds/")
setwd(wd)
# file to save spatial data
gpkg_dir <- "data"
if (!dir.exists(gpkg_dir)) {
  dir.create(gpkg_dir)
}
f_ne <- file.path(gpkg_dir, "gis-data.gpkg")

# List all raster files in the folder with extension .asc
raster_slamms <- list.files(path = "data/slamms/", 
                           pattern = "\\.ASC$",  # Adjust this pattern for other extensions like ".tif"
                           full.names = TRUE)  # Get the full file paths

# Load all rasters as a list
raster_list <- lapply(raster_slamms, rast)  # Use rast() for 'terra', raster() for 'raster'



# Initialize with the extent of the first raster
combined_extent <- ext(raster_list[[1]])

# Loop through the remaining rasters and compute the union of their extents
for (i in 2:3) {
  combined_extent <- union(combined_extent, ext(raster_list[[i]]))
}


# Combine the list of polygons into a single SpatialPolygons object
combined_polygons <- as.polygons(combined_extent)

# Convert to an sf object for easier manipulation if needed
combined_sf <- st_as_sf(combined_polygons)

# Step 1: Check the CRS of both polygons
crs_ngom <- crs(ngom)
crs_combined <- crs(combined_polygons)
if (is.na(crs_combined) || crs_combined == "") {
  # Assign the correct CRS (replace 'EPSG:YOUR_CODE_HERE' with the actual EPSG code)
  crs(combined_polygons) <- "EPSG:4326"
}
# Step 2: If the CRS is different, reproject 'combined_polygons' to match 'ngom'

combined_polygons <- project(combined_polygons, crs(ngom))


# Step 3: Perform the intersection
intersection_result <- intersect(ngom, combined_polygons)

# Step 4: Plot the result to visualize the intersection
plot(intersection_result, col = "lightblue", border = "darkblue")

# Step 2: Extend all rasters to match the combined extent
extended_rasters <- lapply(raster_list[1:3], function(r) {
  extend(r, combined_extent)
})

# Step 3: Resample all rasters to match the resolution of raster_list[[1]], keeping the new extent
resampled_rasters <- lapply(extended_rasters, function(r) {
  resample(r, extended_rasters[[1]], method = "bilinear")  # or "ngb" for categorical
})

# Step 2: Resample each raster to match resolution but use the full combined extent
aligned_rasters <- lapply(raster_list, function(r) {
  # Extend each raster to match the combined extent
  r_extended <- extend(r, combined_extent)
  
  # Resample to match resolution (but not extent)
  resample(r_extended, raster_list[[1]], method = "bilinear")  # or "ngb" for categorical data
})
# Step 3: Now mosaic the aligned rasters to combine their extents
combined_raster <- do.call(mosaic, c(resampled_rasters, fun = mean))

# Plot the combined raster to see the result
plot(combined_raster)


# Convert the raster to polygons
polygons_whole <- as.polygons(combined_raster, dissolve = TRUE)

# Plot the resulting polygons
plot(polygons_whole, main = "Combined Raster Boundaries", col = "lightblue", border = "darkblue")

shp_dir <- "data/NOS80K_NGoMshoreline/nos80k_NGoM.shp"
ngom<-vect(shp_dir)
plot(ngom)
shp_dir <- "data/The_GCOOS_Region_3195904373465877595/GCOOS_Region_US-EEZ.shp"
GCOOS<-vect(shp_dir)
plot(GCOOS)


# Save as a shapefile
writeVector(ngom_polygon, "ngom_polygon.shp", filetype = "ESRI Shapefile")

list.files("ngom_polygon.shp", "ngom", ignore.case = TRUE, full.names = TRUE) %>% 
  unlink()

# Load world map data for countries
world <- ne_countries(scale = "medium", returnclass = "sf")

# Load state boundaries
states <- ne_states(returnclass = "sf")

# Load coastline data
coastline <- ne_download(scale = "medium", type = "coastline", category = "physical", returnclass = "sf")

# Convert SpatVector to sf object
ngom_sf <- sf::st_as_sf(ngom)

# Repair geometries
ngom_sf <- st_make_valid(ngom_sf)
world <- st_make_valid(world)
states <- st_make_valid(states)
coastline <- st_make_valid(coastline)
# Crop world map by ngom
world_cropped <- st_intersection(world, ngom_sf)

# Crop state boundaries by ngom
states_cropped <- st_intersection(states, ngom_sf)

# Crop coastline data by ngom
coastline_cropped <- st_intersection(coastline, ngom_sf)

# output
unlink(f_ne)
write_sf(states_cropped, f_ne, "ne_land")
write_sf(coastline_cropped, f_ne, "ne_coast")
write_sf(world_cropped, f_ne, "ne_country_lines")
write_sf(ngom_sf, f_ne, "ngom")

# Read the 'ngom' data
ngom_sf <- read_sf("data/gis-data.gpkg", layer = "ngom")

# Plot the 'ngom' data
plot(ngom_sf, main = "NGOM Data")



# Plot the bounding box
plot(bbox_polygon, col = "lightgrey", main = "Bounding Box Polygon")



shp_dir <- "data/NOS80K_NGoMshoreline/nos80k_NGoM.shp"
ngom<-vect(shp_dir)
plot(ngom)



# file to save spatial data
gpkg_dir <- "data"
if (!dir.exists(gpkg_dir)) {
  dir.create(gpkg_dir)
}
f_ne_bcr <- file.path(gpkg_dir, "gis-data-bcr.gpkg")

#download bcrs
tmp_dir <- normalizePath(tempdir())
tmp_bcr <- file.path(tmp_dir, "bcr.zip")
paste0("https://www.birdscanada.org/download/gislab/", 
       "bcr_terrestrial_shape.zip") %>% 
  download.file(destfile = tmp_bcr)
unzip(tmp_bcr, exdir = tmp_dir)
list.files(tmp_dir)

bcr <- file.path(tmp_dir, "BCR_Terrestrial/BCR_Terrestrial_master_International.shp") %>% 
  read_sf() %>% 
  select(bcr_code = BCR, bcr_name = Label) %>% 
  filter(bcr_code == 27)
# clean up
# clean up
list.files(tmp_dir, "bcr", ignore.case = TRUE, full.names = TRUE) %>% 
  unlink()

# political boundaries
# land border with lakes removed
ne_land <- ne_download(scale = 50, category = "cultural",
                       type = "admin_0_countries_lakes",
                       returnclass = "sf") %>%
  filter(CONTINENT == "North America") %>%
  st_set_precision(1e6) %>%
  st_union()
# country lines
# downloaded globally then filtered to north america with st_intersect()
ne_country_lines <- ne_download(scale = 50, category = "cultural",
                                type = "admin_0_boundary_lines_land",
                                returnclass = "sf") %>% 
  st_geometry()
ne_country_lines <- st_intersects(ne_country_lines, ne_land, sparse = FALSE) %>%
  as.logical() %>%
  {ne_country_lines[.]}
# states, north america
ne_state_lines <- ne_download(scale = 50, category = "cultural",
                              type = "admin_1_states_provinces_lines",
                              returnclass = "sf") %>%
  filter(ADM0_A3 %in% c("USA", "CAN")) %>%
  mutate(iso_a2 = recode(ADM0_A3, USA = "US", CAN = "CAN")) %>% 
  select(country = ADM0_NAME)

# Crop world map by ngom
ne_land_cropped <- st_intersection(ne_land, ngom)

# Crop state boundaries by ngom
ne_state_lines_cropped <- st_intersection(ne_state_lines, ngom)

# Crop coastline data by ngom
ne_country_lines_cropped <- st_intersection(ne_country_lines, ngom)


# output
unlink(f_ne_bcr)
write_sf(ne_land, f_ne_bcr, "ne_land")
write_sf(ne_country_lines, f_ne_bcr, "ne_country_lines")
write_sf(ne_state_lines, f_ne_bcr, "ne_state_lines")
write_sf(bcr, f_ne_bcr, "bcr")
