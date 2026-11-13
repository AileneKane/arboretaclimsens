# This extracts data from the General GBIF occurence records

rm(list=ls())
library(readxl)
library(rgbif)
library(dplyr)
library(openxlsx)

species_name <- "Abies koreana";   # Define a species name
species_lookup <- name_backbone(name = species_name)
my_target_key <- species_lookup$usageKey

# Filter parameters for info relevant to provenance
filters <- c("scientificName", "eventDate", "recordedBy", "stateProvince",
             "municipality", "county", "higherGeography","verbatimLocality",
             "locality", "verbatimElevation", "elevation", "decimalLatitude", 
             "decimalLongitude")

# Sample from occurences (parameters such as country, stateProvince, recordedBy
# used to narrow down results depending on available info from arboretum)
collector_data <- occ_search(
  taxonKey = my_target_key, 
  hasCoordinate = TRUE,
  #country = "",
  #stateProvince = "",
  #recordedBy = "",
  limit = 120)

# To subset relevant info (related to provenance): 
collector_data <- collector_data$data

cols <- filters[which(filters %in% names(collector_data))]
entry <- collector_data %>% select(all_of(cols))

# Write to excel sheet
write.xlsx(entry, "gbif_search.xlsx")
