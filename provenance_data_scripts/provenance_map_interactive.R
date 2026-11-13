# Libraries for Map
library(leaflet)
library(dplyr)
library(shiny)
library(RColorBrewer)

# Importing Excel Data (Replace excel file path with the one from your download)
library(readxl)
filepath <- "C:/Temporal Ecology Lab/arboretaclimsens/provenance_trees.xlsx"
provenance_trees <- read_excel(filepath)

# TEMPORARY: Removing rows with missing provenance info
entries <- na.omit(provenance_trees)

# Setting up color palette for map
pal <- colorFactor(
  palette = colorRampPalette(
    brewer.pal(11, "Spectral")
  )(length(unique(entries$Genus))),
  domain = sort(unique(entries$Genus))
)

# Setup of main map and navigation bar layout 
ui <- fluidPage(
  titlePanel("Provenance Locations of Cored Trees"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "genus",
        "Select Genus",
        choices = c("All", sort(unique(entries$Genus))),
        selected = "All",
        multiple = FALSE
      ),
      selectInput(
        "G_or_A",
        "Select Gymnosperm or Angiosperm",
        choices = c("Both", "Gymnosperm", "Angiosperm"),
        selected = "Both",
        multiple = FALSE
      )
    ),
    
    mainPanel(
      
      leafletOutput("map", height = 700)
      
    )
  )
)

# Server responds to user selection
server <- function(input, output, session) {
  
  # Reactive filtered dataset
  filtered_data <- reactive({
    
    if (input$genus == "All") {
      if (input$G_or_A == "Both") {
        entries
      } else {
        entries %>%
          filter(G_or_A == input$G_or_A)
      }
    } else {
      entries %>%
        filter(Genus == input$genus)
        
    }
    
  })
  
  # Initial map render
  output$map <- renderLeaflet({
    
    leaflet() %>%
      addTiles(options = tileOptions(noWrap = TRUE)) %>%
      setView(
        lng = mean(entries$Longitude, na.rm = TRUE),
        lat = mean(entries$Latitude, na.rm = TRUE),
        zoom = 10
      )
    
  })
  
  # Update markers without redrawing map
  observe({
    
    leafletProxy("map", data = filtered_data()) %>%
      
      clearGroup("trees") %>%
      
      addCircleMarkers(
        lng = ~Longitude,
        lat = ~Latitude,
        
        fillColor = ~pal(Genus),
        fillOpacity = 0.7,
        group = "trees",
        
        color = "black",
        weight = 1,
        radius = 6,
        
        clusterOptions = markerClusterOptions(
          spiderfyOnMaxZoom = TRUE,
          spiderfyDistanceMultiplier = 2,
          showCoverageOnHover = FALSE,
          zoomToBoundsOnClick = TRUE
        ),
        
        label = ~lapply(
          paste0(
            "<b>Tree ID:</b> ", TreeID,
            "<br><b>Genus:</b> ", Genus,
            "<br><b>Species:</b> ", Species,
            "<br>", G_or_A
          ),HTML)
      )
    
  })
  
}

shinyApp(ui, server)


