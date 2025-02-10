library(sf)
library(rnaturalearth)
library(dplyr)
library(terra)
library(raster)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
install.packages("devtools")
select <- dplyr::select

# CV boundary sf

wd <- ("/Users/liyingnceas/GitHub/gulf_of_mexico_slr_birds/")
setwd(wd)
# file to save spatial data
gpkg_dir <- "data"
if (!dir.exists(gpkg_dir)) {
  dir.create(gpkg_dir)
}
f_ne <- file.path(gpkg_dir, "gis-data-ngom.gpkg")


shp_dir <- "data/NOS80K_NGoMshoreline/nos80k_NGoM.shp"
ngom<-vect(shp_dir)
plot(ngom)
shp_dir <- "data/The_GCOOS_Region_3195904373465877595/GCOOS_Region_US-EEZ.shp"
GCOOS<-vect(shp_dir)
plot(GCOOS)
# Open the NetCDF file as a SpatRaster (terra object)
elev <- rast("/Users/liyingnceas/GitHub/gulf_of_mexico_slr_birds/data/northern_gulf_1_navd88_2010.nc")
# Assign CRS to the raster
crs(elev) <- "EPSG:4326"
# Define bounding box coordinates

# Convert the raster extent to a polygon get elev extent
elev_extent_sf <- st_as_sfc(st_bbox(elev))


# Clip the 'ngom' polygons using the raster extent
ngom_clipped <- st_intersection(ngom_sf, elev_extent_sf)


# Check validity of the 'ngom' geometries
ngom_sf <- st_as_sf(ngom)
# make geometry valid
ngom_sf <- st_make_valid(ngom_sf)  # Fix the geometry if it's not valid
# Check the validity of each geometry in 'ngom_sf'
validity_ngom <- st_is_valid(ngom_sf)
# Print out the indices of invalid geometries
invalid_indices_ngom <- which(!validity_ngom)
if (length(invalid_indices_ngom) > 0) {
  print(paste("Invalid geometries in ngom_sf at indices:", paste(invalid_indices_ngom, collapse = ", ")))
}

#ignore the invalid geometry
ngom_clipped <-st_intersection(ngom_sf[1:9990,], elev_extent_sf)
# Save as a shapefile
st_write(ngom_clipped, "ngom_clipped_polygon.shp", filetype = "ESRI Shapefile")

list.files("ngom_clipped_polygon.shp", "ngom", ignore.case = TRUE, full.names = TRUE) %>% 
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
ngom_clipped <- st_make_valid(ngom_clipped)
world <- st_make_valid(world)
states <- st_make_valid(states)
coastline <- st_make_valid(coastline)
# Crop world map by ngom
world_cropped <- st_intersection(world, ngom_clipped)

# Crop state boundaries by ngom
states_cropped <- st_intersection(states, ngom_clipped)

# Crop coastline data by ngom
coastline_cropped <- st_intersection(coastline, ngom_clipped)

# output
unlink(f_ne)
write_sf(states_cropped, f_ne, "ne_land")
write_sf(coastline_cropped, f_ne, "ne_coast")
write_sf(world_cropped, f_ne, "ne_country_lines")
write_sf(ngom_clipped, f_ne, "ngom")


# Plot the bounding box
plot(bbox_polygon, col = "lightgrey", main = "Bounding Box Polygon")



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
