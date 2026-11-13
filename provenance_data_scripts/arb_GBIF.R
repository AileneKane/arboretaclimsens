# This extracts data from the Harvard Herbaria Records on GBIF: 
# https://www.gbif.org/dataset/861e6afe-f762-11e1-a439-00145eb45e9a
library(rgbif)
library(dplyr)
library(openxlsx)

# Define species name
species_name <- "Abies koreana";
species_lookup <- name_backbone(name = species_name)
my_target_key <- species_lookup$usageKey

# Filter Parameters
filters <- c("scientificName", "eventDate", "recordedBy", "stateProvince",
             "municipality", "county", "higherGeography","verbatimLocality",
             "locality", "verbatimElevation", "elevation", "decimalLatitude", 
             "decimalLongitude")

# Sample from occurences dataset
harvard_sample <- occ_search(
  datasetKey = "861e6afe-f762-11e1-a439-00145eb45e9a",
  taxonKey = my_target_key,
  #country = "",
  #stateProvince = "",
  #recordedBy = "",
  limit = 120
)

# To subset relevant info (related to provenance):
harvard_sample <- harvard_sample$data

cols <- filters[which(filters %in% names(harvard_sample))]
entry <- harvard_sample %>% select(all_of(cols))

# Write to excel sheet
write.xlsx(entry, "gbif_search.xlsx")
