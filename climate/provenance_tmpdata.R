# ** Extracting tmp Data from CRU dataset for provenances and arboretum **
library(terra)
library(dplyr)
library(readxl)

# Importing Excel Data (Modify filepath depending on where file is stored)
filepath <- "C:/Temporal Ecology Lab/arboretaclimsens/all_provenances.xlsx"
provenance_trees <- read_excel(filepath)

# TEMPORARY: Removing rows with missing provenance info
entries <- na.omit(provenance_trees)

# Extract list of coordinates from the data frame
target_coords <- entries %>% select("Longitude", "Latitude")
arb_coords <- matrix(c(-71.12589, 42.29793), ncol=2)

# Raster for CRU tmp. data
cru_filepath <- "C:/Temporal Ecology Lab/arboretaclimsens/cru_ts4.09.1901.2024.tmp.dat.nc/cru_ts4.09.1901.2024.tmp.dat.nc"
tmp_raster <- rast(cru_filepath)

# Starting Month: 1950-01
# Ending Month: 2009-12
lyr_start <- (1950-1901)*12+1
lyr_end <- (2010-1901)*12

dates <- time(tmp_raster)[lyr_start:lyr_end]
tmp_data <- tmp_raster[[lyr_start:lyr_end]]

# terra:extract to obtain values, remove 'ID' column and assign dates
cell_tmp <- terra::extract(tmp_data, target_coords)
cell_tmp <- cell_tmp[, -1]
names(cell_tmp) <- dates

# Sum all of the mean temps for each month (columns) 
# over the 30 year period, for each row (site)
avg_1950_1979 <- rowSums(cell_tmp[,1:360])
avg_1980_2009 <- rowSums(cell_tmp[,361:720])

# Divide by total number of months
avg_1950_1979 <- avg_1950_1979 / 360
avg_1980_2009 <- avg_1980_2009 / 360

# Add result back to data frame
entries$avg_1950_1979 <- avg_1950_1979
entries$avg_1980_2009 <- avg_1980_2009

# Repeat for Arboretum Site
arb_tmp <- terra::extract(tmp_data, arb_coords)
names(arb_tmp) <- dates

arb_1950_1979 <- (rowSums(arb_tmp[, 1:360])) / 360
arb_1980_2009 <- (rowSums(arb_tmp[, 361:720])) / 360

# *** Write to excel file for easier access ***
#library(openxlsx)
#write.xlsx(entries, "climaticData.xlsx")

