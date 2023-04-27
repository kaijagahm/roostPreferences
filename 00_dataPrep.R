# Prepare data for movement analysis
# Setup -------------------------------------------------------------------
## Load packages -----------------------------------------------------------
library(vultureUtils)
library(sf)
library(tidyverse)
library(move)
library(feather)
library(readxl)

## Download data from movebank
# download data from movebank (just a subset of the times for now)
base::load("movebankCredentials/pw.Rda")
MB.LoginObject <- move::movebankLogin(username = "kaijagahm", password = pw)
rm(pw)

dat <- vultureUtils::downloadVultures(loginObject = MB.LoginObject, removeDup = T, dfConvert = T, quiet = T, dateTimeStartUTC = "2020-09-01 00:00", dateTimeEndUTC = "2021-11-30 11:59") # 2020 through the 2022 nb season. We don't have the full 2023 breeding season yet.
write_feather(dat, "data/dat.feather")
dat <- read_feather("data/dat.feather")
## fix trackId
dat <- dat %>%
  mutate(trackId = as.character(trackId),
         trackId = case_when(trackId == "E03" ~ "E03w",
                             TRUE ~ trackId))

# Add the Nili_id
ww <- read_excel("data/whoswho_vultures_20230315_new.xlsx", sheet = "all gps tags")[,1:35] %>%
  dplyr::select(Nili_id, Movebank_id) %>%
  distinct()

all(dat$trackId %in% ww$Movebank_id) # true

dat2 <- left_join(dat, ww, by = c("trackId" = "Movebank_id"))

## annotate the data with periods to remove
toRemove <- read_excel("data/whoswho_vultures_20230315_new.xlsx", sheet = "periods_to_remove")

toRemove <- toRemove %>%
  dplyr::select(Nili_id,
                "trackId" = Movebank_id,
                remove_start,
                remove_end,
                reason) %>%
  mutate(across(c(remove_start, remove_end), .fns = function(x){
    lubridate::ymd(x)
  })) %>%
  dplyr::filter(!is.na(remove_end))

toRemove_long <- toRemove %>%
  group_by(Nili_id, reason) %>%
  # sequence of daily dates for each corresponding start, end elements
  dplyr::mutate(dateOnly = map2(remove_start, remove_end, seq, by = "1 day")) %>%
  # unnest the list column
  unnest(cols = c(dateOnly)) %>%
  # remove any duplicate rows
  distinct() %>%
  dplyr::select(-c(remove_start, remove_end)) %>%
  rename("status" = reason)

# Join to the original data
datAnnot <- dat2 %>%
  left_join(toRemove_long, by = c("Nili_id", "dateOnly")) %>%
  mutate(status = replace_na(status, "valid"))
nrow(datAnnot) == nrow(dat2) #T

# Clean the data
## Region masking, downsampling, removal of speed outliers, setting altitude outliers to NA, etc.
mask <- sf::st_read("data/CutOffRegion.kml")
datAnnotCleaned <- vultureUtils::cleanData(dataset = datAnnot, mask = mask, inMaskThreshold = 0.33, removeVars = F, idCol = "Nili_id", downsample = F)
save(datAnnotCleaned, file = "data/datAnnotCleaned.Rda")

## Load data ---------------------------------------------------------------
load("data/datAnnotCleaned.Rda")

# Fix time zone so dates make sense ---------------------------------------
## Overwrite the dateOnly column from the new times
datAnnotCleaned <- datAnnotCleaned %>%
  mutate(timestampIsrael = lubridate::with_tz(timestamp, tzone = "Israel"),
         dateOnly = lubridate::date(timestampIsrael))

# Split into seasons ------------------------------------------------------
datAnnotCleaned <- datAnnotCleaned %>%
  mutate(month = lubridate::month(timestampIsrael),
         year = lubridate::year(timestampIsrael),
         season = case_when(month %in% 7:11 ~ "nb",
                            month %in% c(12, 1:6) ~ "b"),
         seasonUnique = case_when(season == "nb" ~ paste(year, season, sep = "_"),
                                  season == "b" & month == 12 ~ 
                                    paste(year+1, season, sep = "_"),
                                  TRUE ~ paste(year, season, sep = "_")))

# Separate the seasons -----------------------------
seasons <- datAnnotCleaned %>%
  group_by(seasonUnique) %>%
  group_split(.keep = T)
seasonNames <- map_chr(seasons, ~.x$seasonUnique[1])

# Restrict to southern individuals ----------------------------------------
# Based on previous investigations for the 2022 breeding and non-breeding seasons, have found that a good cutoff for southern vs. non-southern is 3550000 (in ITM)
## Transform to SF object, so we can get centroids
seasonsSF <- map(seasons, ~.x %>%
                   sf::st_as_sf(coords = c("location_long", "location_lat"), remove = F) %>%
                   sf::st_set_crs("WGS84") %>%
                   sf::st_transform(32636))

## Get centroids, so we can see who's "southern" for that season.
centroids <- map(seasonsSF, ~.x %>%
                   group_by(Nili_id) %>%
                   summarize(geometry = st_union(geometry)) %>%
                   st_centroid())

## Examine a histogram of centroid latitudes 
walk(centroids, ~hist(st_coordinates(.x)[,2])) # looks like 3550000 is generally a good cutoff point here.

## Get southern individuals for each season, so we can filter the data
southernIndivs <- map(centroids, ~.x %>%
                        filter(st_coordinates(.)[,2] < 3550000) %>%
                        pull(Nili_id))

## Remove individuals not in the south
seasons <- map2(.x = seasons, .y = southernIndivs, ~.x %>% filter(Nili_id %in% .y))

# XX note: not removing individuals not tracked for long enough, or individuals without a lot of points per day (so, this script differs from the one in MvmtSoc in that way). I don't think we need to subset individuals this way for the purpose of roost preference analysis.

# Get roosts for each season ----------------------------------------------
roosts_seasons <- purrr::map(seasons, ~vultureUtils::get_roosts_df(df = .x, id = "Nili_id")) 

roosts_seasons <- roosts_seasons %>%
  map(., ~st_as_sf(.x, crs = "WGS84", coords = c("location_long", "location_lat")))

# Export the data
save(seasons, file = "data/seasons.Rda")
save(roosts_seasons, file = "data/roosts_seasons.Rda")
