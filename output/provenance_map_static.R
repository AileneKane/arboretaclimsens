# Static Map for visualizing provenance data
# Zoom to see full figure with legend
library(dplyr)
library(readxl)
library(sf)
library(ggplot2)
library(rnaturalearth)

# Replace with file path to where provenance_trees file is downloaded
filepath <- "C:/Temporal Ecology Lab/arboretaclimsens/provenance_trees.xlsx"
provenance_trees <- read_excel(filepath)

# TEMPORARY: Removing rows with missing provenance info
entries <- na.omit(provenance_trees)

# Filter Variables (Currently displays all)
# Ex: replace with genus <- c("Abies", "Betula")
a_g <- c("Gymnosperm","Angiosperm")
genus <- sort(unique(entries$Genus))
species <- sort(unique(entries$Species))

# Adjust range for map span (e.g North America: -130, -60, 25, 60)
# Default for world map (-180, 180, -90, 90)
min_long <- -180
max_long <- 180
min_lat <- -90
max_lat <- 90

# Filter dataset according to filter variables
entries_filtered <- entries %>% filter(G_or_A %in% a_g, Genus %in% genus, Species %in% species)

# Filter dataset according to location constraints
entries_local <- entries_filtered %>%
  filter(
    Longitude >= min_long,
    Longitude <= max_long,
    Latitude >= min_lat,
    Latitude <= max_lat
  )

# Draw Base Map
basemap_sf <- ne_countries(
  scale = 110,
  returnclass = "sf"
)

# Convert data frame to spatial object for ggplot2
trees_sf <- st_as_sf(
  entries_local,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)

# Plot data points on map
ggplot() +
  geom_sf(data = basemap_sf, fill = "grey95") +
  geom_sf(data = trees_sf, aes(color = Genus), size = 2, show.legend = TRUE) +
  coord_sf(
    xlim = c(min_long, max_long),
    ylim = c(min_lat, max_lat)
  ) +
  theme_minimal() +
  labs(
    title = "Provenance Locations of Cored Trees"
  ) + theme(
    legend.position = "left"
  )
