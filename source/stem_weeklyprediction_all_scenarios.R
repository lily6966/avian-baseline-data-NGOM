setwd("/home/lli/R/")
parent.dir <- "/home/lli/R/"

# -----------------------
source.dir <- paste(parent.dir,"source/",sep="")
source(paste(source.dir,"vMR.library.R",sep=""))
source(paste(source.dir,"vMR_baseModels.R",sep=""))
#------------
library(fields)
require(rpart)
library(mgcv)
library(sf)
library(raster)
library(lubridate)
library(scam)
library(PresenceAbsence)
library(verification)
library(fields)
library(gridExtra)
library(tidyverse)
library(sf)
library(raster)
library(lubridate)
library(ranger)
library(scam)
library(PresenceAbsence)
library(verification)
library(fields)
library(gridExtra)
library(tidyverse)
library(randomForest)
library(dggridR)
library(mlr)
library(dplyr)
library(ISOweek)
library(gbm)
# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection
#------------

# set random number seed to insure fully repeatable results
set.seed(1)



# ----------------------------------------------
# 1
# ----------------------------------------------
species_collection1 <- c("blkski", "brnpel", "greegr", "laugul", "motduc", "royter", "seaspa", "clarai")

scenarios<-c("baseline","initial","low", "high")
species<-species_collection1[2]
for (species in species_collection1){
  for (i in 1:1){
    scenario <- scenarios[i]

# ebird data
ebird <- read_csv(paste0("data/ebd_gcoos", species, "_bbox_zf.csv")) %>% 
  # year required to join to habitat data
  mutate(year = year(observation_date))
ebird_count <- ebird %>% 
  mutate(protocol_type = factor(protocol_type, 
                                levels = c("Stationary" , "Traveling"))) %>%
  # remove observations with no count
  filter(!is.na(observation_count))
# habitat covariates with crops categorized 
habitat <- read_csv(paste0("data/covariates/", species, "_pland_elev_north_east_slope_location-20192021_GCOOS_6km.csv")) %>% 
  mutate(year = as.integer(year))

# combine ebird_count and habitat data
ebird_habitat_count <- inner_join(ebird_count, habitat, by = c("locality_id", "year"))

# prediction surface

r <- raster("data/prediction-surface-gcoos.tif")

# load gis data for making maps
map_proj <- map_proj <- "ESRI:102003"
ne_land <- read_sf("data/gis-data-gcoos.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
gcoos <- read_sf("data/gis-data-gcoos.gpkg", "gcoos") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_country_lines <- read_sf("data/gis-data-gcoos.gpkg", "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

# generate hexagonal grid with ~ 10 km betweeen cells
dggs <- dgconstruct(spacing = 10)
# get hexagonal cell id and week number for each checklist
checklist_cell_count <- ebird_habitat_count %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
         year = year(observation_date),
         week = week(observation_date))

# sample one checklist per grid cell per week
# sample detection/non-detection independently 
ebird_ss_count <- checklist_cell_count %>% 
  group_by(species_observed, year, week, cell) %>% 
  sample_n(size = 1) %>% 
  ungroup() 

ebird_split_count <- ebird_ss_count %>% 
  # select only the columns to be used in the model
  dplyr::select(-checklist_id, -observer_id, -sampling_event_identifier,
                -scientific_name, -species_observed, -state_code,
                -locality_id, -all_species_reported, -observation_date,   
                -cell) %>% 
  drop_na() 


# split 80/20
ebird_split_count <- ebird_split_count %>% 
  split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))
map_int(ebird_split_count, nrow)
# ------------------------------------------
# STEM Data structures 
# ------------------------------------------
train.data.time <- NULL      #add a time column make a time column to count from the first day to the last day of the 10 years
for (i in 2019:2021) {
  if (i == 2019) {
    train.data.time.year <-ebird_split_count$train %>% 
      filter(year == i)  %>% 
      mutate(time = day_of_year + (i-2019)*365)
  } else if (i>=2020 ) {
    train.data.time.year <-ebird_split_count$train %>% 
      filter(year == i)  %>% 
      mutate(time = day_of_year + (i-2019)*365 +1)
  } 
  
  
  train.data.time <- bind_rows(train.data.time, train.data.time.year)
}



range01 <- function(x){(x-min(x))/(max(x)-min(x))} #function to normalize longitude, latitude, time

train.data <- train.data.time%>% 
  mutate(x=range01(longitude), y=range01(latitude), t=range01(time)) %>% 
  select(-longitude, -latitude, -time)

train.locations <- train.data.time %>% dplyr::select(longitude, latitude, time) %>% 
  mutate(x=range01(longitude), y=range01(latitude), t=range01(time)) %>% 
  select(-longitude, -latitude, -time)

# ------------------------------------------
# repeat STEM Data structures for test data
# ------------------------------------------
test.locs <- NULL
for (i in 2019:2021) {
  if (i == 2019) {
    test.locs.year <-ebird_split_count$test %>% 
      filter(year == i) %>% 
      mutate(time = day_of_year + (i-2019)*365)
  } else if (i>=2020) {
    test.locs.year <- ebird_split_count$test %>% 
      filter(year == i) %>% 
      mutate(time = day_of_year + (i-2019)*365 +1)
  } 
  test.locs <- bind_rows(test.locs, test.locs.year)
}

test.data <- test.locs %>% 
  mutate(x=range01(longitude), y=range01(latitude), t=range01(time)) %>% 
  select(-longitude, -latitude, -time)

test.locations <- test.locs %>% dplyr::select(longitude, latitude, time) %>% 
  mutate(x=range01(longitude), y=range01(latitude), t=range01(time)) %>% 
  select(-longitude, -latitude, -time)


train.data.y <- train.data %>% mutate(y=observation_count) %>% select(-observation_count, -week)


##uncomment if want to use GBM_train instead of GBM_NAremove_train
d6.stem <- stem(
  data = train.data,
  locations = train.locations,
  width = c(0.3, 0.3, 0.008),
  n = 10,
  baseModel = GBM_train,
  # No BaseModel parameters
  # BaseModel parameters
  bm.formula = "observation_count ~ slopeness_median+slopeness_sd+year+day_of_year+time_observations_started+duration_minutes+
    effort_distance_km+number_observers+time_observations_started+estuarine_wetland+palustrine_wetland+NA+pasture+agriculture+elevation_median+elevation_sd+
    northness_median+northness_sd+eastness_median+eastness_sd+water+urban+barren+grassland+forest",
  min.data.size = 100,
  distribution = "gaussian",
  n.trees=1000,
  bag.fraction=.95,
  interaction.depth=2,
  verbose=T)

d6.stem <- stem( 	
  data = train.data.y,
  locations = train.locations,
  width = c(0.3, 0.3, 0.008),
  n = 10,
  baseModel = GBM_NAremove_train, 
  # No BaseModel parameters
  min.data.size = 100,
  distribution = "gaussian",
  n.trees=1000,
  bag.fraction=.95,
  interaction.depth=2,
  verbose=T)
pred6.stem <- predict.stem(
  stem.obj = d6.stem, 
  new.data = test.data,
  new.locations = test.locations,
  baseModel = GBM_predict,
  n.trees = 1000)		


#regression plot of observed and predicted
obs_count <- select(test.data, obs = observation_count)

m_gaussian_pred <-
  tibble(species = species, pred = pred6.stem$mean) %>%
  bind_cols(obs_count) %>%
  filter(pred >= 0)

correlation <- cor(m_gaussian_pred$obs, m_gaussian_pred$pred)
pseudo_r_squared <- correlation^2
# plot predicted vs. observed
ticks <- c(0, 1, 10, 100, 1000)
mx <- round(max(m_gaussian_pred$obs))
plot<-ggplot(m_gaussian_pred) +
  aes(x = log10(obs + 1),
      y = log10(pred + 1)) +
  geom_jitter(alpha = 0.2, height = 0) +
  # y = x line
  geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
  # area where counts off by a factor of 10
  geom_area(data = tibble(x = log10(seq(0, mx-1)+1),
                          y = log10(seq(0, mx-1) / 10+1)),
            mapping = aes(x = x, y = y),
            fill = "red", alpha = 0.2) +
  # loess fit
  geom_smooth(method = "loess",
              method.args = list(span = 2 / 3, degree = 1)) +
  scale_x_continuous(breaks = log10(ticks + 1), labels = ticks) +
  scale_y_continuous(breaks = log10(ticks + 1), labels = ticks) +
  labs(x = "Observed count",
       y = "Predicted count") +
  facet_wrap(~ species, nrow = 1) +
  geom_text(aes(label = paste("Pseudo R-squared =", round(pseudo_r_squared, 3))), 
            x = 0, y = 0, vjust = -32, hjust = 0)#vjust = -5, hjust = -2 to postion bottom right
ggsave(paste0(species, scenario, "_obs_pred.pdf"), plot = plot, width = 5, height = 5, units = "in")
# ----------------------------------------------------------------------
# Predictor importance of GBM 
# ---------------------------------------------------------------------- 

# # Create a dataframe
# pi_stem <- data.frame(predictor = character(), importance = numeric(), stringsAsFactors = FALSE)
# predictors <- c("canopy", "year", "day_of_year", "time_observations_started", "duration_minutes",
#                 "effort_distance_km", "number_observers","field_row", "rice", "perennial",
#                 "pasture", "grain", "elevation_median", "elevation_sd",
#                 "northness_median", "northness_sd", "eastness_median", "eastness_sd", "alfalfa", "seasonal_wetlands",
#                 "semi_permanent_wetlands", "developed", "barren", "corn", "water", "shrub", "forest")
# #for (i in names(train.data)){
# pi_stem_i<-stem.gbm.predictor.importance(
#   d6.stem,
#   partition.number = 1,
#   n.trees=1000, 
#   predictor.names = predictors) 
# # Add each predictor importance
# 
# #  pi_stem <- rbind(pi_stem, c(i, pi_stem_i))
# #}
# pi_stem_i_narm <- na.omit(pi_stem_i)
# pi_stem <- pi_stem_i_narm %>% select(-year)%>%colMeans()
# pi_stem_d <- data.frame(predictor=names(pi_stem[4:29]), importance=pi_stem[4:29], row.names = NULL)
# 
# write_csv(pi_stem_d, paste(species, "_pi",".csv",sep=""))



# prediction surface
pred_surface_count <- read_csv(paste0("data/pred_surface/pland-elev_nort_east_slope_prediction-surface_", scenario, "_20192021_6km.csv"))


# add effort covariates to prediction 
r <- raster("data/prediction-surface-cv.tif")
pred_surface_eff_count <- pred_surface_count %>% 
  mutate(
    time_observations_started = 6.114894,
    duration_minutes = 60,
    effort_distance_km = 1,
    number_observers = 1,
    protocol_type = "Traveling") 

range01 <- function(x){(x-min(x))/(max(x)-min(x))} #function to normalize longitude, latitude, time

# ------------------------------------------------------------------------------------------------------
#start to prediction week by week and average across the 10 year for each week of the 52 weeks of a year
# --------------------------------------------------------------------------------------------------

# create directory to save the raster

tif_dir <- paste0(species, "_new_abundance_30days_", scenario)
if (!dir.exists(tif_dir)) {
  dir.create(tif_dir)
}


###jjj
for (wk in 1:52){
  #divide a year into a loop of 52 weeks
  #create a place holder for each week's 10 year assembly
  pred_surface_gbm_count_10year <- NULL
  for (yr in 2051:2060){
    
    # Get the start date of the week
    s <- paste(yr, "-","W", sprintf("%02d", wk),"-", 1, sep = "")
    start_date <- ISOweek2date(s)
    # start_date <- sprintf("%04d-W%02s-1", yr, wk) #another way
    # Get the end date of the week
    e <- paste(yr, "-","W", sprintf("%02d", wk),"-", 7, sep = "")
    end_date <- ISOweek2date(e)
    
    # Create a sequence of dates for the entire week
    week_dates <- seq(start_date, end_date, by = "day")
    
    
    #select each of the 10 year prediction surface data
    pred_surface_eff_count_year <- pred_surface_eff_count %>%
      filter(year==yr) %>% mutate(year=year-40) #align the year with train data year range
    
    # Create a sequence of dates for the entire week
    week_dates <- seq(start_date, end_date, by = "day")
    pred_surface_gbm_count_week <- NULL
    for (k in 1:7) {
      #create a place holder for the 7 days in each of the 52 weeks
      
      #add a day_of_year co variate
      pred_surface_eff_count_day  <- pred_surface_eff_count_year%>%
        mutate(obeservation_date = week_dates[k]-years(40),  #align the year with train data year range
               day_of_year = yday(week_dates[k])
        )
      
      #add a time column based on day_of_year and scale the day count to a ten-year span
      if (yr>=2051 & yr<=2052) {
        pred_surface_day.time <- pred_surface_eff_count_day  %>%
          mutate(time = day_of_year + (yr-2051)*365)
      } else if (yr>2052 & yr<=2056) {
        pred_surface_day.time <- pred_surface_eff_count_day  %>%
          mutate(time = day_of_year + (yr-2051)*365 +1)
      } else {
        pred_surface_day.time <- pred_surface_eff_count_day  %>%
          mutate(time = day_of_year + (yr-2051)*365 +2)
      }
      # #create dummy years match with train data
      # 
      #       if (yr==2052 | yr==2056 ) {
      #         pred_surface_day.locations <- pred_surface_day.time %>% dplyr::select(longitude, latitude, day_of_year) %>%
      #           mutate(x=range01(longitude), y=range01(latitude), t=(day_of_year-1)/(366-1)) %>%
      #           select(-longitude, -latitude, -day_of_year)
      # 
      #         pred_surface_day.data <- pred_surface_day.time %>%
      #           mutate(x=range01(longitude), y=range01(latitude), t=(day_of_year-1)/(366-1)) %>%
      #           select(-longitude, -latitude)
      #       #normalize the geometry and time for stem
      # 
      #         } else {
      #          pred_surface_day.locations <- pred_surface_day.time %>% dplyr::select(longitude, latitude, day_of_year) %>%
      #           mutate(x=range01(longitude), y=range01(latitude), t=(day_of_year-1)/(365-1)) %>%
      #             select(-longitude, -latitude, -day_of_year)
      # 
      #           pred_surface_day.data <- pred_surface_day.time %>%
      #             mutate(x=range01(longitude), y=range01(latitude), t=(day_of_year-1)/(365-1)) %>%
      #             select(-longitude, -latitude)
      #       }
      
      
      pred_surface_day.locations <- pred_surface_day.time %>% dplyr::select(longitude, latitude, time) %>%
        mutate(x=range01(longitude), y=range01(latitude), t=(time-min(train.data.time$time))/(3648-min(train.data.time$time))) %>%
        select(-longitude, -latitude, -time)
      
      pred_surface_day.data <- pred_surface_day.time %>%
        mutate(x=range01(longitude), y=range01(latitude), t=(time-min(train.data.time$time))/(3648-min(train.data.time$time))) %>%
        select(-longitude, -latitude, -time)
      
      #make prediction on the single day
      pred2_surface_day.stem <- predict.stem(
        stem.obj = d6.stem, 
        new.data = pred_surface_day.data, 
        new.locations = pred_surface_day.locations,
        baseModel = GBM_predict,
        n.trees = 1000,
        na.action = na.omit)
      
      #grab the mean of local neiborhood base model prediction
      pred_surface_gbm_count_day <- pred2_surface_day.stem$mean
      #combine daily prediction to a dataframe
      pred_surface_gbm_count_week <- bind_cols(pred_surface_gbm_count_week,  pred_surface_gbm_count_day)
      
      
      
    }
    # Calculate the mean of each column
    weekly_mean <- rowMeans(pred_surface_gbm_count_week, na.rm =TRUE)
    # Add coordinates to abundance
    pred_surface_gbm_count_10year <- bind_cols(pred_surface_gbm_count_10year,  weekly_mean)
    # 
  } 
  
  if (!any(complete.cases(pred_surface_gbm_count_10year))) {
    print(paste("no date week",wk))
    # If all values are NA, exit and proceed to the next loop iteration
    next
  }
  mean_10year <- rowMeans(pred_surface_gbm_count_10year, na.rm = TRUE)
  #add coordinates to weekly abundance averaged across 10 years
  pred_surface_gbm_er_count_week <- bind_cols(pred_surface_day.locations, abundance = mean_10year) %>% 
    # convert to spatial featurest
    select(y, x, abundance) %>% 
    #denomalize
    mutate(latitude=y*(max(pred_surface_day.time$latitude)-min(pred_surface_day.time$latitude))+min(pred_surface_day.time$latitude),  #reverse coordinate normalization
           longitude=x*(max(pred_surface_day.time$longitude)-min(pred_surface_day.time$longitude))+min(pred_surface_day.time$longitude),
           abundance = pmin(pmax(abundance, 0), 1)) %>% select(-y,-x)
  #then rasterize the points 
  r_pred_surface_count_week <- pred_surface_gbm_er_count_week %>% 
    na.omit() %>% #remove abundance NA value to rasterize
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
    st_transform(crs = projection(r)) %>%
    # rasterize
    rasterize(r)
  r_pred_surface_count <- r_pred_surface_count_week[[-1]]
  # ----------------------------------------------
  # 3
  # ----------------------------------------------  
  
  # save the raster
  
  raster_write <- paste(species, "_new_abundance_30days_",scenario,"/week",wk,".tif",sep="")
  writeRaster(r_pred_surface_count, raster_write,  
              overwrite = TRUE)
}
###jjj

  }
  # prediction surface
  pred_surface_count_cdl <- read_csv("data/pland-elev_nort_east_canopy_prediction-surface_greyel_cdl.csv")
  
  
  
  
  # add effort covariates to prediction 
  r <- raster("data/prediction-surface-cv.tif")
  pred_surface_eff_count_cdl <- pred_surface_count_cdl %>% 
    mutate(
      time_observations_started = 6.114894,
      duration_minutes = 60,
      effort_distance_km = 1,
      number_observers = 1,
      protocol_type = "Traveling") 
  
  range01 <- function(x){(x-min(x))/(max(x)-min(x))} #function to normalize longitude, latitude, time
  
  # ----------------------------------------------
  # 4
  # ----------------------------------------------
  
  # create directory to save the raster
  tif_dir <- paste0(species, "_abundance_30days_cdl")
  if (!dir.exists(tif_dir)) {
    dir.create(tif_dir)
  }
  
  
  ###jjj
  for (wk in 1:52){
    #divide a year into a loop of 52 weeks
    #create a place holder for each week's 10 year assembly
    pred_surface_gbm_count_10year <- NULL
    for (yr in 2011:2020){
      
      # Get the start date of the week
      s <- paste(yr, "-","W", sprintf("%02d", wk),"-", 1, sep = "")
      start_date <- ISOweek2date(s)
      # start_date <- sprintf("%04d-W%02s-1", yr, wk) #another way
      # Get the end date of the week
      e <- paste(yr, "-","W", sprintf("%02d", wk),"-", 7, sep = "")
      end_date <- ISOweek2date(e)
      
      # Create a sequence of dates for the entire week
      week_dates <- seq(start_date, end_date, by = "day")
      
      
      #select each of the 10 year land cover data
      pred_surface_eff_count_year <- pred_surface_eff_count_cdl %>%
        filter(year==yr)
      
      # Create a sequence of dates for the entire week
      week_dates <- seq(start_date, end_date, by = "day")
      pred_surface_gbm_count_week <- NULL
      for (k in 1:7) {
        #create a place holder for the 7 days in each of the 52 weeks
        
        #add a day_of_year co variate
        pred_surface_eff_count_day  <- pred_surface_eff_count_year%>%
          mutate(obeservation_date = week_dates[k]-years(40),
                 day_of_year = yday(week_dates[k])
          )
        
        #add a time column based on day_of_year and scale the day count to a ten-year span
        if (yr>=2011 & yr<=2012) {
          pred_surface_day.time <- pred_surface_eff_count_day  %>%
            mutate(time = day_of_year + (yr-2011)*365)
        } else if (yr>2012 & yr<=2016) {
          pred_surface_day.time <- pred_surface_eff_count_day  %>%
            mutate(time = day_of_year + (yr-2011)*365 +1)
        } else {
          pred_surface_day.time <- pred_surface_eff_count_day  %>%
            mutate(time = day_of_year + (yr-2011)*365 +2)
        }
        # #create dummy years match with train data
        
        # if (i==2052 | i==2056 ) {
        #   pred_surface_day.locations <- pred_surface_day.time %>% dplyr::select(longitude, latitude, day_of_year) %>%
        #     mutate(x=range01(longitude), y=range01(latitude), t=(day_of_year-1)/(366-1)) %>%
        #     select(-longitude, -latitude, -day_of_year)
        #   
        #   pred_surface_day.data <- pred_surface_day.time %>%
        #     mutate(x=range01(longitude), y=range01(latitude), t=(day_of_year-1)/(366-1)) %>%
        #     select(-longitude, -latitude)
        # #normalize the geometry and time for stem 
        # 
        #   } else {
        #    pred_surface_day.locations <- pred_surface_day.time %>% dplyr::select(longitude, latitude, day_of_year) %>%
        #     mutate(x=range01(longitude), y=range01(latitude), t=(day_of_year-1)/(365-1)) %>%
        #       select(-longitude, -latitude, -day_of_year)
        #     
        #     pred_surface_day.data <- pred_surface_day.time %>%
        #       mutate(x=range01(longitude), y=range01(latitude), t=(day_of_year-1)/(365-1)) %>%
        #       select(-longitude, -latitude)
        # }
        # 
        # 
        pred_surface_day.locations <- pred_surface_day.time %>% dplyr::select(longitude, latitude, time) %>%
          mutate(x=range01(longitude), y=range01(latitude), t=(time-min(train.data.time$time))/(max(train.data.time$time)-min(train.data.time$time))) %>%
          select(-longitude, -latitude, -time)
        
        pred_surface_day.data <- pred_surface_day.time %>%
          mutate(x=range01(longitude), y=range01(latitude), t=(time-min(train.data.time$time))/(max(train.data.time$time)-min(train.data.time$time))) %>%
          select(-longitude, -latitude, -time)
        
        #make prediction on the single day
        pred2_surface_day.stem <- predict.stem(
          stem.obj = d6.stem, 
          new.data = pred_surface_day.data, 
          new.locations = pred_surface_day.locations,
          baseModel = GBM_predict,
          n.trees = 1000,
          na.action = na.omit)
        
        #grab the mean of local neiborhood base model prediction
        pred_surface_gbm_count_day <- pred2_surface_day.stem$mean
        #combine daily prediction to a dataframe
        pred_surface_gbm_count_week <- bind_cols(pred_surface_gbm_count_week,  pred_surface_gbm_count_day)
        
        
        
      }
      # Calculate the mean of each column
      weekly_mean <- rowMeans(pred_surface_gbm_count_week, na.rm =TRUE)
      # Add coordinates to abundance
      pred_surface_gbm_count_10year <- bind_cols(pred_surface_gbm_count_10year,  weekly_mean)
      # 
    } 
    
    if (!any(complete.cases(pred_surface_gbm_count_10year))) {
      print(paste("no date week",wk))
      # If all values are NA, exit and proceed to the next loop iteration
      next
    }
    mean_10year <- rowMeans(pred_surface_gbm_count_10year, na.rm = TRUE)
    #add coordinates to weekly abundance averaged across 10 years
    pred_surface_gbm_er_count_week <- bind_cols(pred_surface_day.locations, abundance = mean_10year) %>% 
      # convert to spatial featurest
      select(y, x, abundance) %>% 
      #denomalize
      mutate(latitude=y*(max(pred_surface_day.time$latitude)-min(pred_surface_day.time$latitude))+min(pred_surface_day.time$latitude),  #reverse coordinate normalization
             longitude=x*(max(pred_surface_day.time$longitude)-min(pred_surface_day.time$longitude))+min(pred_surface_day.time$longitude),
             abundance = pmin(pmax(abundance, 0), 1)) %>% select(-y,-x)
    #then rasterize the points 
    r_pred_surface_count_week <- pred_surface_gbm_er_count_week %>% 
      na.omit() %>% #remove abundance NA value to rasterize
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
      st_transform(crs = projection(r)) %>%
      # rasterize
      rasterize(r)
    r_pred_surface_count <- r_pred_surface_count_week[[-1]]
    # ----------------------------------------------
    # 3
    # ----------------------------------------------  
    
    # save the raster
    raster_write <- paste(species, "_abundance_30days_cdl/week",wk,".tif",sep="")
    writeRaster(r_pred_surface_count, raster_write, 
                overwrite = TRUE)
  }
}