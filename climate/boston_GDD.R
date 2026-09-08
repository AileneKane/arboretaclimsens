# This is for calculating GDD from PRISM Data for the arboretum 
# to use for modeling annual climatic variation
library(terra)

# Arboretum Coordinates
lat <- 42.29793
lon <- -71.12589
arboretum_coors <- vect(cbind(lon, lat),crs = "EPSG:4269")

# Folder containing TIFFs in server
prism_dir <- "/data/climate/prism/800m/us/tmean"

# Find all tif files
files <- list.files(prism_dir,
                    pattern = "\\.tif$",
                    recursive = TRUE,
                    full.names = TRUE
)

r <- rast(files)
tmean_data <- terra::extract(r, arboretum_coors)

next_months <- seq(as.Date(paste0(1895, "-02-01")),
                   as.Date(paste0(2025, "-01-01")),
                   by = "month")

# Compute the last day of each month (to get number of days in each month)
last_days <- next_months - 1

tmean_vals <- data.frame(
  yr = as.numeric(format(last_days, "%Y")),
  dy = as.numeric(format(last_days, "%d")),
  tmean = as.numeric(tmean_data[1, -1])
)

# For each month: (Tmean - 5)*(N days in month)
tmean_vals$weighted <- (tmean_vals$tmean - 5) * tmean_vals$dy

# Filter out negative values (e.g. Tmean was below 5 deg C)
tmean_vals <- subset(tmean_vals, weighted >= 0)

# Sum total for each month
annual_GDDs <- tapply(tmean_vals$weighted,
                      tmean_vals$yr,
                      sum,
                      na.rm = TRUE)

# Post processing: GDD = GDD * 0.982 + 136.05
annual_GDDs <- annual_GDDs * 0.982 + 136.05

# Converting to data frame just for copy paste!
#output <- stack(annual_GDDs)
#print(output)


