# ** Extracting Precip. Data from CRU dataset for provenances and arboretum **
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

# Raster for CRU precip. data
cru_filepath <- "C:/Temporal Ecology Lab/arboretaclimsens/cru_ts4.09.1901.2024.pre.dat.nc/cru_ts4.09.1901.2024.pre.dat.nc"
pre_raster <- rast(cru_filepath)

# Starting Month: 1950-01
# Ending Month: 2009-12
lyr_start <- (1950-1901)*12+1
lyr_end <- (2010-1901)*12

dates <- time(pre_raster)[lyr_start:lyr_end]
pre_data <- pre_raster[[lyr_start:lyr_end]]

# terra:extract to obtain values, remove 'ID' column and assign dates
cell_pre <- terra::extract(pre_data, target_coords)
cell_pre <- cell_pre[, -1]
names(cell_pre) <- dates

# Subset data into the two 30-year periods
pre_1950_1979 <- cell_pre[,1:360]
pre_1980_2009 <- cell_pre[,361:720]

# Sum the monthly precip. for each year to get annual precip.
yr <- format(as.Date(colnames(pre_1950_1979)), "%Y")
annual_1950_1979 <- sapply(unique(yr), function(y) {
  rowSums(pre_1950_1979[, yr == y], na.rm = TRUE)
})

yr <- format(as.Date(colnames(pre_1980_2009)), "%Y")
annual_1980_2009 <- sapply(unique(yr), function(y) {
  rowSums(pre_1980_2009[, yr == y], na.rm = TRUE)
})

# Take the 30-year mean for annual precip.
means_1950_1979 <- rowMeans(annual_1950_1979, na.rm = FALSE)
means_1980_2009 <- rowMeans(annual_1980_2009, na.rm = FALSE)

entries$pre_1950_1979 <- means_1950_1979
entries$pre_1980_2009 <- means_1980_2009

# Repeat for Arboretum Site
arb_pre <- terra::extract(pre_data, arb_coords)
names(arb_pre) <- dates

arb_1950_1979 <- arb_pre[, 1:360]
arb_1980_2009 <- arb_pre[, 361:720]

yr <- format(as.Date(colnames(pre_1950_1979)), "%Y")
arb_1950_1979 <- sapply(unique(yr), function(y) {
  rowSums(arb_1950_1979[, yr == y], na.rm = TRUE)})

yr <- format(as.Date(colnames(pre_1980_2009)), "%Y")
arb_1980_2009 <- sapply(unique(yr), function(y) {
  rowSums(arb_1980_2009[, yr == y], na.rm = TRUE)})

arb_1950_1979 <- mean(arb_1950_1979, na.rm = FALSE)
arb_1980_2009 <- mean(arb_1980_2009, na.rm = FALSE)


