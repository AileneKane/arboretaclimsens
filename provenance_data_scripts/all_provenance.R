# This script is used to combine provenance records for which the direct locality is known
# with sampled coordinates from GBIF for the remaining trees

# Output: Write a list of all trees and their provenance coordinates to all_provenances

library(readxl)
library(openxlsx)
library(terra)
library(dplyr)

filepath <- "C:/Temporal Ecology Lab/arboretaclimsens/provenance_trees.xlsx"
provenance_trees <- read_excel(filepath, sheet = "exact_locality")
coords_samples <- read_excel(filepath, sheet = "GBIF_sample")

# Save all trees with precise coordinates in entries
entries <- na.omit(provenance_trees)

# Save trees without coordinate info, will draw from GBIF samples for these
sample_cores <- provenance_trees %>% filter(is.na(Latitude))
set.seed(23)

for (i in 1:nrow(sample_cores)){
  x <- coords_samples %>% filter(Country == sample_cores[i,]$Country & 
                                   Species == sample_cores[i,]$SciName)
  index <- sample(1:nrow(x), 1) # Select a random index from list of coords
  sample_cores[i,]$Latitude <- x[index,]$Latitude
  sample_cores[i,]$Longitude <- x[index,]$Longitude
}

# Combine both sets of data into one data frame, and write to all_provenances
provenance_trees <- bind_rows(entries, sample_cores)
write.xlsx(provenance_trees, "all_provenances.xlsx")

