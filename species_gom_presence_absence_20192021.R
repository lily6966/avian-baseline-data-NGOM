#setup
install.packages("remotes")
remotes::install_github("mstrimas/ebppackages")
setwd("/Users/liyingnceas/GitHub/gulf_of_mexico_slr_birds/")

library(auk)
library(lubridate)
library(sf)
library(gridExtra)
library(tidyverse)
library(raster)
library(rgdal) ## confirm the GDAL version being used
library(tiff)
library(viridis)
library(exactextractr)

#set ebird directory
auk::auk_set_ebd_path("/home/lli/R", overwrite = TRUE)
setwd("/home/lli/R")
getwd()
# ------------------------------------------
# extract ebird data
# ------------------------------------------
select <- dplyr::select

# Open the NetCDF file as a SpatRaster (terra object)
elev <- rast("data/northern_gulf_1_navd88_2010.nc")
# Assign CRS to the raster
crs(elev) <- "EPSG:4326"
plot(elev)
  


species<-c("blkski", "brnpel", "greegr", "laugul", "motduc", "royter", "seaspa", "clarai")
species_long<-c("Black Skimmer", "Brown Pelican", "Great Egret", "Laughing Gull", "Mottled Duck","Royal Tern", "Seaside Sparrow", "Clapper Rail")


for (i in 2:2){
# resolve namespace conflicts


# # setup data directory, produce presence-only data
# dir.create("data/ebird/ameavo", showWarnings = FALSE)

ebd <- auk_ebd(paste0("data/ebird/", species[i],"/ebd_US_",species[i],"_201901_202112_smp_relSep-2024.txt"))
ebd <- auk_ebd("data/ebd_US_smp_relSep-2024/ebd_US_relSep-2024.txt")

# #define the filters that you want to apply to the EBD
# ebd_filters <- ebd %>% 
#   # restrict to the standard traveling and stationary count protocols and
#   #keep complete checklists with auk_complete()
#   auk_protocol(protocol = c("Stationary", "Traveling")) %>% 
#   auk_complete()

# #output files produce presence-only data
# data_dir <- "data/ebird/ameavo"
# if (!dir.exists(data_dir)) {
#   dir.create(data_dir)
# }
# f_ebd1 <- file.path(data_dir, "data/ebird/ameavo/ebd_US-CA_ameavo_201012_202011_relMar-2024_protocal.txt")
# 
# # only run if the files don't already exist
# if (!file.exists(f_ebd)) {
#   auk_filter(ebd_filters, file = f_ebd)
# }


#produce zero-filled, detection/non-detection data (also called presence/absence data)
#dir.create("data/ebird/whimbr", showWarnings = FALSE)

# define an EBD checklists reference and a set of filters
ebd_filters <- auk_ebd("data/ebd_US_smp_relSep-2024/ebd_US_relSep-2024.txt",
                        #use the same sampling event data for the -2023 ebird data
                            file_sampling = "data/ebd_US_smp_relSep-2024/ebd_US_relSep-2024_sampling.txt") %>% 
  # species: common and scientific names can be mixed
  #auk_species(species = c(species_long[i])) %>%
  # bbox: long and lat in decimal degrees
  # formatted as `c(lng_min, lat_min, lng_max, lat_max)`
  auk_bbox(bbox = c(-90.7504,28.49986,-84.99986,31.25014)) %>%   #ext(elev)
  # date: use standard ISO date format `"YYYY-MM-DD"`
  auk_date(date = c("2019-01-01", "2021-12-30")) %>%
  # restrict to the standard traveling and stationary count protocols and
  #keep complete checklists with auk_complete()
  auk_protocol(protocol = c("Stationary", "Traveling")) %>% 
  auk_complete()

data_dir <- paste0("data/ebird/","gom")
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}
f_ebd <- file.path(data_dir, paste0("ebd_gom_","pa","_201901_202112.txt"))
f_sampling <- file.path(data_dir, paste0("ebd_gom_","checklists", "_201901_202112.txt"))
# Reading a text file
data_text <- read.table("data/ebird/gom/ebd_gom_checklists_201901_202112.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
lines <- readLines("data/ebird/gom/ebd_gom_checklists_201901_202112.txt", n = 20)
cat(lines[15:20], sep = "\n")  # View lines near the issue



head(data_text)  # Preview the data


# Note that  `ebd_filters`  have been set. At this point, we've only defined the filters, 
#not applied them to the EBD. The last step is to use `auk_filter()` to compile the filters 
#into an AWK script and run it to produce two output files: one for the EBD and one for the SED. 
#**This step typically takes several hours to run since the files are so large.** 
#As a result, it's wise to wrap this in an `if` statement, so it's only run once. 
#As noted in the [Introduction](intro-setup-software), Windows users will need to [install Cygwin](https://www.cygwin.com/) for this next step to work.


if (!file.exists(f_ebd)) {
  auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling, overwrite=TRUE)
}

#Importing (reading) data and zero-filling

ebd_zf <- auk_zerofill(f_ebd, f_sampling, collapse = TRUE,overwrite=TRUE)

# function to convert time observation to hours since midnight
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}


ebd_zf <- ebd_zf %>% 
# clean up variables
 
  mutate(
    # convert X to NA
    observation_count = if_else(observation_count == "X", 
                                NA_character_, observation_count),
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type != "Traveling", 
                                 0, effort_distance_km),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )
#impose some more consistent structure on the data by filtering observations on the effort variables.
#restricting checklists to less than 5 hours long and 5 km in length, and with 10 or fewer observers. 
ebd_zf_filtered <- ebd_zf %>% 
  filter(
    # effort filters
    duration_minutes <= 5 * 60,
    effort_distance_km <= 5,
    # 10 or fewer observers
    number_observers <= 10)

#remove redundant variables and create csv file

ebird_pa <- ebd_zf_filtered %>%  #write ebird presence/absence data
  select(checklist_id, observer_id, sampling_event_identifier,
         scientific_name,
         observation_count, species_observed, 
         state_code, locality_id, latitude, longitude,
         protocol_type, all_species_reported,
         observation_date, year, day_of_year,
         time_observations_started, 
         duration_minutes, effort_distance_km,
         number_observers)

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
unique_scientific_names <- unique(ebird_pa$scientific_name)




# Define batch size
batch_size <- 10

# Get the total number of unique scientific names
total_scientific_names <- length(unique_scientific_names)

# Process in batches
for (i in seq(1, total_scientific_names, by = batch_size)) {
  # Get the current batch of scientific names
  current_batch <- unique_scientific_names[i:min(i + batch_size - 1, total_scientific_names)]
  
  # Loop over the current batch
  for (sci_name in current_batch) {
    # Filter the data for that specific scientific name
    data_subset <- ebird_pa %>% filter(scientific_name == sci_name)
    
    # Join with scientific_to_common to get the common name
    data_subset_with_common <- data_subset %>% left_join(scientific_to_common, by = "scientific_name")
    
    # Create the file name by sanitizing the scientific name
    file_name <- paste0("data/ebird/species_by_name/", sanitize_filename(sci_name), ".csv")
    
    # Save the subset to a CSV file
    write.csv(data_subset, file = file_name, row.names = FALSE)
    
    # Create the file name by sanitizing the common name
    common_name <- data_subset_with_common$common_name[1]
    if (is.na(common_name) || common_name == "") {
      common_name <- sanitize_filename(sci_name) # Fallback to scientific name
    }
    file_name <- paste0("data/ebird/species_by_common_name/", sanitize_filename(common_name), ".csv")
    
    # Save the subset to a CSV file
    write.csv(data_subset_with_common, file = file_name, row.names = FALSE)
  }
  
  
}


library(dplyr)



# Define a chunk size for processing rows
chunk_size <- 1000000

# Function to process a single scientific name
process_scientific_name <- function(sci_name, ebird_pa, scientific_to_common, chunk_size) {
  tryCatch({
    # Filter for rows corresponding to the scientific name
    sci_data <- ebird_pa %>% filter(scientific_name == sci_name)
    
    if (nrow(sci_data) > chunk_size) {
      # Process in chunks if the dataset is too large
      num_chunks <- ceiling(nrow(sci_data) / chunk_size)
      
      for (chunk_idx in seq_len(num_chunks)) {
        sci_chunk <- sci_data %>% slice((chunk_idx - 1) * chunk_size + 1:min(chunk_idx * chunk_size, n()))
        
        # Join with the scientific_to_common
        sci_chunk_with_common <- sci_chunk %>%
          left_join(scientific_to_common, by = "scientific_name")
        
        # Common name extraction
        common_name <- sci_chunk_with_common$common_name[1]
        if (is.na(common_name) || common_name == "") {
          common_name <- sanitize_filename(sci_name) # Fallback to scientific name
        }
        
        # Create file names
        file_name_sci <- paste0("data/ebird/species_by_name/", sanitize_filename(sci_name), "_chunk_", chunk_idx, ".csv")
        file_name_common <- paste0("data/ebird/species_by_common_name/", sanitize_filename(common_name), "_chunk_", chunk_idx, ".csv")
        
        # Save the chunk data
        write.csv(sci_chunk, file = file_name_sci, row.names = FALSE)
        write.csv(sci_chunk_with_common, file = file_name_common, row.names = FALSE)
      }
    } else {
      # Process normally if data size is manageable
      sci_data_with_common <- sci_data %>%
        left_join(scientific_to_common, by = "scientific_name")
      
      # Common name extraction
      common_name <- sci_data_with_common$common_name[1]
      if (is.na(common_name) || common_name == "") {
        common_name <- sanitize_filename(sci_name) # Fallback to scientific name
      }
      
      # Create file names
      file_name_sci <- paste0("data/ebird/species_by_name/", sanitize_filename(sci_name), ".csv")
      file_name_common <- paste0("data/ebird/species_by_common_name/", sanitize_filename(common_name), ".csv")
      
      # Save the data
      write.csv(sci_data, file = file_name_sci, row.names = FALSE)
      write.csv(sci_data_with_common, file = file_name_common, row.names = FALSE)
    }
    
    cat("Processed:", sci_name, "\n")
  }, error = function(e) {
    # Handle errors for a specific scientific name
    cat("Error processing", sci_name, ":", e$message, "\n")
  })
}

# Loop through unique scientific names and process each
unique_scientific_names <- unique(ebird_pa$scientific_name)
for (sci_name in unique_scientific_names) {
  process_scientific_name(sci_name, ebird_pa, scientific_to_common, chunk_size)
}

unique_scientific_names <- unique(ebird_pa$scientific_name) %>%
  sort() %>%
  .[grepl("^[M-Z]", ., ignore.case = TRUE)]


# Check the output files in your working directory




library(data.table)
# Example: Load large data and split into manageable chunks
chunk_size <- 1e5  # Number of rows per chunk (adjust as needed)
file_path <- paste0("data/ebird/all_species/",  "ebd_gom", "_bbox_zf.csv")

# Open connection to write CSV
fwrite(head(ebird_pa, 0), file = file_path)  # Write header first

# Save chunks of data
for (start_row in seq(1, nrow(ebird_pa), by = chunk_size)) {
  end_row <- min(start_row + chunk_size - 1, nrow(ebird_pa))
  fwrite(ebird_pa[start_row:end_row, ], file = file_path, append = TRUE)
}








# how many species?
distinct_count <- n_distinct(ebird_pa$scientific_name)

cat("Number of distinct scientific names:", distinct_count, "\n")






#this goes back to line to clean up presence-only data) 

# ebd_filters_r <- read_ebd(f_ebd1)
# 
# 
# 
# # clean up variables
# ebd_filters_r <- ebd_filters_r %>% 
#   mutate(
#     # convert X to NA
#     observation_count = if_else(observation_count == "X", 
#                                 NA_character_, observation_count),
#     observation_count = as.integer(observation_count),
#     # effort_distance_km to 0 for non-travelling counts
#     effort_distance_km = if_else(protocol_type != "Traveling", 
#                                  0, effort_distance_km),
#     # convert time to decimal hours since midnight
#     time_observations_started = time_to_decimal(time_observations_started),
#     # split date into year and day of year
#     year = year(observation_date),
#     day_of_year = yday(observation_date)
#   )

# #additional filtering accounting for variation in detectability by filtering observations on the effort variables
# ebd_filtered <- ebd_filters_r %>% 
#   filter(
#     # effort filters
#     duration_minutes <= 5 * 60,
#     effort_distance_km <= 5,
#     # 10 or fewer observers
#     number_observers <= 10)
# 
# 
# #remove redundant variables, keeping only the variables we want for modelling. 
# 
# ebird <- ebd_filtered %>% 
#   select(checklist_id, observer_id, sampling_event_identifier,
#          scientific_name,
#          observation_count,species_observed, 
#          state_code, locality_id, latitude, longitude,
#          protocol_type, all_species_reported,
#          observation_date, year, day_of_year,
#          time_observations_started, 
#          duration_minutes, effort_distance_km,
#          number_observers)
# write_csv(ebird, "data/ebird/least sandpiper/ebd_US-CA_leasan_201012_202011_relSep-2022_protocal.csv", na = "")
# 






#Exploratory analysis and visualization of the presence/absence data
# load and project gis data
map_proj <- "ESRI:102003"


ne_land <- read_sf("data/gis-data.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
cv <- read_sf("data/gis-data.gpkg", "cv") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_country_lines <- read_sf("data/gis-data.gpkg", "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

# ------------------------------------------
# Map ebird observation Data 
# ------------------------------------------

# prepare ebird data for mapping
ebird_pa_sf <- ebird_pa %>% 
  # convert to spatial points
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = map_proj)

# map
par(mar = c(0.25, 0.25, 0.25, 0.25))
# set up plot area
plot(st_geometry(ebird_sf), col = NA)
# contextual gis data
plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
plot(cv, col = "#cccccc", border = NA, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
# ebird observations
plot(filter(ebird_pa_sf, species_observed) %>% st_geometry(),
     pch = 19, cex = 0.3, col = alpha("#4daf4a", 1),
     add = TRUE)

# not observed
plot(st_geometry(ebird_pa_sf),
     pch = 19, cex = 0.1, col = alpha("#555555", 0.25),
     add = TRUE)

# legend
legend("bottomright", bty = "n",
       col = c("#555555", "#4daf4a"),
       legend = c("eBird checklists", "least sandpiper"),
       pch = 19)
box()
par(new = TRUE, mar = c(0, 0, 3, 0))
title("least sandpiper eBird Observations \nDec 2010-Nov 2020, Central Valley_CA")

# ------------------------------------------
# Compiling environmental covariates 
# ------------------------------------------


wd <- ("/Users/lily/Box/Thesis/SDM/code/data/")
elev <- raster(paste0(wd, "CV_20101117_gmted_med075.tif"))



# load the landcover data
landcover <- list.files("data/landcover", "^cdl_cv_albers_", 
                        full.names = TRUE) %>% 
  stack()
landcover
# label layers with year
landcover <- names(landcover) %>% 
  str_extract("[0-9]{4}") %>% 
  paste0("y", .) %>% 
  setNames(landcover, .)
#load cv boundary and crop elev 
cv <-read_sf("data/gis-data.gpkg", "cv") %>% 
  # project to the cdl projection
  st_transform(crs = projection(landcover))

# crop, buffer cv_cov by 10 km to provide a little wiggly room
cv_extent <- cv %>% 
  st_buffer(dist = 10000) %>% 
  st_transform(crs = projection(elev))

# load the dust landcover data
landcover_dust <- list.files("data/landcover_dust", "^scn418.sc.it1.ts", 
                        full.names = TRUE) %>% 
  stack()
# label layers with year
landcover_dust <- names(landcover_dust) %>% 
  str_extract("[0-9]{4}") %>% 
  paste0("y", .) %>% 
  setNames(landcover_dust, .)


# Set the extent of each landcover_dust raster stack to same extent and projection as landcover
landcover_dust <- landcover_dust %>% setExtent(extent(landcover)) %>% 
  projectRaster(crs = projection(landcover))

# Check the result
extent(landcover_dust)
extent(landcover)  
  
  
elev <-  elev %>% 
  crop(., cv_extent) %>% 
  projectRaster(crs = projection(landcover))




# load ebird data

ebird <- read_csv("data/ebd_lobcur_bbox_zf.csv")

# buffer to create neighborhood (3 km by 3 km) around each checklist point 
# with a radius equal to 100/2 CDL cells 
neighborhood_radius <- 100 * ceiling(max(res(landcover))) / 2
neighborhood_radius
ebird_buff <- ebird %>% 
  distinct(year = format(observation_date, "%Y"),
           locality_id, latitude, longitude) %>% 
  # match the ebird year with corresponding year landcover data
  mutate(year_lc = paste0("y", year))  %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  # transform to cdl projection
  st_transform(crs = projection(landcover)) %>% 
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

crop_category <- read.csv(file = "crop_category_lucas.csv")
pland_crop_category <- inner_join(pland, crop_category,  by = "landcover")

# tranform to wide format, filling in implicit missing values with 0s%>% 


pland_crop_category_wide <- pland_crop_category %>% 
  # dplyr::select(-landcover, -crop_name) %>%
  group_by(locality_id, year, landuse) %>% 
  mutate(pland_crop = sum (pland)) %>% 
  dplyr::select(-landcover, -crop_name, -pland) %>%
  group_by(locality_id, year) %>%
  distinct(landuse, .keep_all = TRUE) %>% 
  # dplyr::select(-pland, -landcover, -crop_name) %>% 
  pivot_wider(names_from = landuse, 
              values_from = pland_crop, 
              values_fill = list(pland_crop = 0))




# make a map
par(mar = c(0.25, 0.25, 2, 0.25))
t <- str_glue("Proportion of Rice\n",
              "CDL Landcover in 2012")
plot(rice_cover_raster, axes = FALSE, box = FALSE, col = viridis(10), main = t)


# add additional environmental covariates 
# extract elevation values


elev_checklists <- NULL
for (yr in names(landcover)) {
    # get the buffered checklists for a given year and extract elevation values and calculate median and sd
      regions <- ebird_buff$data[[which(yr == ebird_buff$year_lc)]]
      locs <- st_set_geometry(regions, NULL) 
      elev_checklists_year <- exact_extract(elev, regions, progress = FALSE) %>% 
          map_dfr(~ tibble(elevation_median = mean(.$value, na.rm = TRUE),
                             +                      elevation_sd = sd(.$value, na.rm = TRUE))) %>% 
         # join to lookup table to get locality_id
         bind_cols(locs, .)
     # bind to results
        elev_checklists <- bind_rows(elev_checklists, elev_checklists_year)
       }

# checklist covariates
pland_elev_checklist <- inner_join(pland_crop_category_wide, elev_checklists, by = c("locality_id", "year"))
write_csv(pland_elev_checklist, "data/pland_elev_lobcur.csv")
#repeat this process to calculate the elevation covariates for the prediction surface
# create prediction surface using landcover buffered neighborhood
# to confine the prediction to the study area

library(sf)
# agg_factor <- round(2 * neighborhood_radius / res(landcover))
# agg_factor
# r <- raster(landcover) %>% 
#   aggregate(agg_factor)  #a template raster with cells equal in size to the neighborhoods
# 
# #load cv shapefile
# cv_dir <- "data/central_valley_aq.shp" 
# cv_cov <- readOGR(cv_dir) 
# 
# r <- cv_cov %>% 
#   st_transform(crs = projection(r)) %>% 
#   rasterize(r, field = 1) %>% 
#   # remove any empty cells at edges
#   trim()
# r <- writeRaster(r, filename = "data/prediction-surface-cv.tif", overwrite = FALSE)
#load prediction surface
# ------------------------------------------
# create cell_buff for both cdl and dust
# ------------------------------------------
r <- raster("data/prediction-surface-cv.tif")

cell_buff <- NULL
for (yr in 2051:2060) {
  r_centers <- rasterToPoints(r, spatial = TRUE) %>% 
    st_as_sf() %>% 
    transmute(id = row_number()) %>% 
    mutate(year=as.character(yr), year_lc = paste0("y", yr))
  r_cells <- st_buffer(r_centers, dist = neighborhood_radius)
  cell_buff <- bind_rows(cell_buff, r_cells)
  # nest by year to  match to land cover data from the corresponding year
}
cell_buff <- nest(cell_buff, data = c(year, id, geometry))

cell_buff_cdl <- NULL
for (yr in 2011:2020) {
  r_centers <- rasterToPoints(r, spatial = TRUE) %>% 
    st_as_sf() %>% 
    transmute(id = row_number()) %>% 
    mutate(year=as.character(yr), year_lc = paste0("y", yr))
  r_cells <- st_buffer(r_centers, dist = neighborhood_radius)
  cell_buff_cdl <- bind_rows(cell_buff_cdl, r_cells)
  # nest by year to  match to land cover data from the corresponding year
}
cell_buff_cdl <- nest(cell_buff_cdl, data = c(year, id, geometry))

lc_extract_pred_cdl <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year
  regions <- cell_buff_cdl$data[[which(yr == cell_buff_cdl$year_lc)]]
  # get landcover values within each buffered checklist area
  ee <- exact_extract(landcover[[yr]], regions, progress = FALSE)
  # count the number of each landcover class for each checklist buffer
  ee_count <- map(ee, ~ count(., landcover = value))
  # attach the year and locality id back to the checklists
  ee_summ <- tibble(st_drop_geometry(regions), data = ee_count) %>% 
    unnest(data)
  # bind to results
  lc_extract_pred_cdl <- bind_rows(lc_extract_pred_cdl, ee_summ)
}
# ------------------------------------------
# creat pland with cdl first
# ------------------------------------------


# calculate the percent for each landcover class
pland_pred_cdl <- lc_extract_pred_cdl %>% 
  group_by(id, year) %>% 
  mutate(pland = n / sum(n)) %>% 
  ungroup() %>% 
  dplyr::select(-n) %>% 
  # remove NAs after tallying so pland is relative to total number of cells
  filter(!is.na(landcover))

#change landcover class name to crop and other landuse names

pland_pred_crop_category_cdl <- inner_join(pland_pred_cdl, crop_category,  by = "landcover")

# tranform to wide format, filling in implicit missing values with 0s%>% 


pland_pred_crop_category_cdl_wide <- pland_pred_crop_category_cdl %>% 
  # dplyr::select(-landcover, -crop_name) %>%
  group_by(id, landuse) %>% 
  mutate(pland_crop = sum (pland)) %>% 
  dplyr::select(-landcover_dust, -pland) %>%
  group_by(id) %>%
  distinct(landuse, .keep_all = TRUE) %>% 
  # dplyr::select(-pland, -landcover, -crop_name) %>% 
  pivot_wider(names_from = landuse, 
              values_from = pland_crop, 
              values_fill = list(pland_crop = 0))


# map landcover classes using pland
# join in coordinates

# map landcover classes using pland
# join in coordinates
pland_coords <- NULL
for (yr in 2051:2060) {
  pland_pred_year <- pland_pred_crop_category_cdl_wide %>% filter(year == as.character(yr))
  pland_coords_year <- st_transform(r_centers, crs = 4326) %>% 
    st_coordinates() %>% 
    as.data.frame() %>% 
    cbind(id = r_centers$id, .) %>% 
    rename(longitude = X, latitude = Y) %>% 
    inner_join(pland_pred_year, by = "id")
  pland_coords <- bind_rows(pland_coords, pland_coords_year)
}

landcover_dust[['y2055']]
# ------------------------------------------
# create pland with dust
# ------------------------------------------
# iterate over all years extracting landcover_dust for all landcover_dust pixels in each cell
lc_pred_extract <- NULL
for (yr in names(landcover_dust)) {
  # get the buffered checklists for a given year
  regions <- cell_buff$data[[which(yr == cell_buff$year_lc)]]
  # get landcover values within each buffered checklist area
  ee <- exact_extract(landcover_dust[[yr]], regions, progress = FALSE)
  # count the number of each landcover class for each checklist buffer
  ee_count <- map(ee, ~ count(., landcover = value))
  # attach the year and locality id back to the checklists
  ee_summ <- tibble(st_drop_geometry(regions), data = ee_count) %>% 
    unnest(data)

  # bind to results
  lc_pred_extract <- bind_rows(lc_pred_extract, ee_summ)
}

summary(lc_pred_extract)
#calculate PLAND: the proportion of each land cover class within the neighborhood.
pland_pred <- lc_pred_extract %>% 
  # calculate proporiton
  group_by(id, year) %>% 
  mutate(pland = n / sum(n)) %>% 
  ungroup() %>% 
  dplyr::select(-n) %>% 
  # remove NAs after tallying so pland is relative to total number of cells
  filter(!is.na(landcover))


#change landcover class name to crop and other landuse names
crop_category_dust <- read.csv(file = "crop_category_dust.csv")
pland_pred_crop_category <- inner_join(pland_pred, crop_category_dust,  by = "landcover")

# tranform to wide format, filling in implicit missing values with 0s%>% 


pland_pred_crop_category_wide <- pland_pred_crop_category %>% 
  # dplyr::select(-landcover, -crop_name) %>%
  group_by(id, landuse) %>% 
  mutate(pland_crop = sum (pland)) %>% 
  dplyr::select(-landcover_dust, -pland) %>%
  group_by(id) %>%
  distinct(landuse, .keep_all = TRUE) %>% 
  # dplyr::select(-pland, -landcover, -crop_name) %>% 
  pivot_wider(names_from = landuse, 
              values_from = pland_crop, 
              values_fill = list(pland_crop = 0))


# map landcover classes using pland
# join in coordinates

# map landcover classes using pland
# join in coordinates
pland_coords_dust <- NULL
for (yr in 2051:2060) {
  pland_pred_year <- pland_pred_crop_category_wide %>% filter(year == as.character(yr))
  pland_coords_year <- st_transform(r_centers, crs = 4326) %>% 
    st_coordinates() %>% 
    as.data.frame() %>% 
    cbind(id = r_centers$id, .) %>% 
    rename(longitude = X, latitude = Y) %>% 
    inner_join(pland_pred_year, by = "id")
  pland_coords_dust <- bind_rows(pland_coords_dust, pland_coords_year)
}

# -----------------------------------------------------
# create maps of both cdl and dust pland for comparison
# -----------------------------------------------------
#Convert data frame format to spatial format of pland metrics
rice_cover <- pland_coords_dust %>% 
  filter(year == 2052) %>%
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) 
# rasterize points
#names(rice_cover)[names(rice_cover)=="3"] <- "rice"
rice_cover_raster <-  rasterize(rice_cover, r, field = "rice") 
# project to albers equal-area for mapping
rice_cover_raster <- rice_cover_raster %>% projectRaster(crs = "ESRI:102003") %>% 
  trim()

rice_cover_raster


# ------------------------------------------
# adding elevation to pred surface
# ------------------------------------------



locs_pred <- NULL
for (yr in names(landcover_dust)) {
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

pland_elev_pred_crop <- inner_join(pland_coords_crop, locs_pred, by = c("id", "year"))

write_csv(pland_elev_pred_crop, "data/pland-elev_prediction-surface_lobcur_10year.csv")
# ------------------------------------------
# add northness to checklists and pred surface
# ------------------------------------------
#repeat the same process for northness covariates
nort <- raster(paste0(wd, "northness_1KMmd_GMTEDmd.tif"))
nort <-  nort %>% 
  crop(., cv_extent) %>% 
  projectRaster(crs = projection(landcover))
nort 
nort_checklists <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- ebird_buff$data[[which(yr == ebird_buff$year_lc)]]
  locs_nort <- st_set_geometry(regions, NULL) 
  nort_checklists_year <- exact_extract(nort, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(northness_median = mean(.$value, na.rm = TRUE),
                     northness_sd = sd(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get locality_id
    bind_cols(locs_nort, .)
  # bind to results
  nort_checklists <- bind_rows(nort_checklists, nort_checklists_year)
}

#repeat for prediction surface
nort_pred <- NULL
for (yr in names(landcover_dust)) {
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


# checklist covariates
pland_elev_nort_checklist <- inner_join(pland_elev_checklist, nort_checklists, c("locality_id", "year"))
write_csv(pland_elev_nort_checklist, "data/pland-elev_nort_lobcur.csv")

# prediction surface covariates

pland_elev_nort_pred_crop <- inner_join(pland_elev_pred_crop, nort_pred,  by = c("id", "year"))
write_csv(pland_elev_nort_pred_crop, "data/pland-elev_nort_prediction-surface_lobcur.csv")

# ------------------------------------------
# add eastness to checklists and pred surface
# ------------------------------------------

#repeat with eastness
east <- raster(paste0(wd, "eastness_1KMmd_GMTEDmd.tif"))
east <-  east %>% 
  crop(., cv_extent) %>% 
  projectRaster(crs = projection(landcover))

east_checklists <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- ebird_buff$data[[which(yr == ebird_buff$year_lc)]]
  locs_east <- st_set_geometry(regions, NULL) 
  east_checklists_year <- exact_extract(east, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(eastness_median = mean(.$value, na.rm = TRUE),
                     eastness_sd = sd(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get locality_id
    bind_cols(locs_east, .)
  # bind to results
  east_checklists <- bind_rows(east_checklists, east_checklists_year)
}

#repeat for prediction surface
east_pred <- NULL
for (yr in names(landcover_dust)) {
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



# checklist covariates
pland_elev_nort_east_checklist <- inner_join(pland_elev_nort_checklist, east_checklists, c("locality_id", "year"))
write_csv(pland_elev_nort_east_checklist, "data/pland-elev_nort_east_lobcur1.csv")

# prediction surface covariates

pland_elev_nort_east_pred_crop <- inner_join(pland_elev_nort_pred_crop, east_pred,  by = c("id", "year"))
write_csv(pland_elev_nort_east_pred_crop, "data/pland-elev_nort_east_prediction-surface_lobcur_10year.csv")

typeof(east_pred$id)
#repeat with canopy height
canopy <- raster(paste0(wd, "canopy_height_2019.tif"))
canopy <-  canopy %>% 
  crop(., cv_extent) %>% 
  projectRaster(crs = projection(landcover))

canopy_checklists <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- ebird_buff$data[[which(yr == ebird_buff$year_lc)]]
  locs_canopy <- st_set_geometry(regions, NULL) 
  canopy_checklists_year <- exact_extract(canopy, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(canopy = mean(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get locality_id
    bind_cols(locs_canopy, .)
  # bind to results
  canopy_checklists <- bind_rows(canopy_checklists, canopy_checklists_year)
}
#replace nan with 0
canopy_checklists <- as.data.frame(canopy_checklists)
canopy_checklists <- rapply(canopy_checklists, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
#repeat for prediction surface
canopy_pred <- NULL
for (yr in names(landcover)) {
  # get the buffered checklists for a given year and extract elevation values and calculate median and sd
  regions <- cell_buff$data[[which(yr == cell_buff$year_lc)]]
  locs_canopy_cell <- st_set_geometry(regions, NULL) %>% 
    mutate(id = row_number())
  canopy_cell <- exact_extract(canopy, regions, progress = FALSE) %>% 
    map_dfr(~ tibble(canopy = mean(.$value, na.rm = TRUE))) %>% 
    # join to lookup table to get id
    bind_cols(locs_canopy_cell, .)
  # bind to results
  canopy_pred <- bind_rows(canopy_pred, canopy_cell)
}
#replace nan with 0
canopy_pred <- as.data.frame(canopy_pred)
canopy_pred <- rapply(canopy_pred, f=function(x) ifelse(is.nan(x),0,x), how="replace" )

# checklist covariates

pland_elev_nort_east_canopy_checklist_crop <- inner_join(pland_elev_nort_east_checklist_crop, canopy_checklists, c("locality_id", "year"))
 write_csv(pland_elev_nort_east_canopy_checklist_crop, "data/pland-elev_nort_east_canopy_location-10year-crop.csv")
 # prediction surface covariates

pland_elev_nort_east_canopy_pred_crop <- inner_join(pland_elev_nort_east_pred_crop, canopy_pred,  rm.na= TRUE, by = c("id", "year"))
write_csv(pland_elev_nort_east_canopy_pred_crop, "data/pland-elev_nort_east_canopy_prediction-surface-crop.csv")
pland_elev_nort_east_canopy_pred_crop