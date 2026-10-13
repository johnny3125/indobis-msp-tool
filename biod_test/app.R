library(shiny)
library(leaflet)
library(robis)
library(dplyr)
library(sf)
library(rnaturalearth)
library(zip)

# Function to compute diversity indices safely following OBIS indicator standards
calc_diversity <- function(sp_vec) {
  counts <- table(sp_vec)
  N <- sum(counts)
  p <- counts / N
  
  # Shannon Index (H')
  shannon <- -sum(p * log(p))
  
  # Simpson Index (D)
  simpson <- sum(p^2)
  
  # Species Richness (S)
  S <- length(counts)
  
  # Pielou's Evenness (J')
  evenness <- if (S > 1) shannon / log(S) else 0
  
  # Hurlbert's ES50 (Estimated species richness for sample of 50 individuals)
  # Standard OBIS methodology requires N >= 50
  es50 <- if (N >= 50) {
    sum(1 - exp(lchoose(N - counts, 50) - lchoose(N, 50)))
  } else {
    NA_real_
  }
  
  list(
    shannon = round(shannon, 3),
    simpson = round(simpson, 3),
    evenness = round(evenness, 3),
    es50 = if (is.na(es50)) NA_real_ else round(es50, 3)
  )
}

# --- DEFINE UI ---
ui <- fluidPage(
  titlePanel("IndOBIS Marine Spatial Planning (MSP) Tool - Hexagonal Grid & Indicators"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Set coordinates for your Northern Indian Ocean region of interest:"),
      
      numericInput("xmin", "Min Longitude (Xmin):", value = 60.0, min = -180, max = 180),
      numericInput("xmax", "Max Longitude (Xmax):", value = 100.0, min = -180, max = 180),
      numericInput("ymin", "Min Latitude (Ymin):", value = 5.0, min = -90, max = 90),
      numericInput("ymax", "Max Latitude (Ymax):", value = 25.0, min = -90, max = 90),
      
      actionButton("fetchBtn", "1. Fetch & Mask Land Records", class = "btn-primary", style = "width:100%; margin-top: 10px;"),
      
      hr(),
      helpText("Step 2: Choose metric and calculate/update map locally."),
      
      selectInput("selected_metric", "Select Map Layer / Metric:",
                  choices = c(
                    "Species Richness (S)" = "species_richness",
                    "Shannon-Wiener Index (H')" = "shannon",
                    "Simpson Index (D)" = "simpson",
                    "Pielou Evenness (J')" = "evenness",
                    "Hurlbert ES50" = "es50",
                    "IUCN Threatened Species" = "threatened_count",
                    "Total Observation Records" = "total_records"
                  ),
                  selected = "species_richness"),
      
      actionButton("calcBtn", "2. Calculate & Update Map", class = "btn-success", style = "width:100%; margin-top: 5px;"),
      
      hr(),
      downloadButton("downloadCSV", "Export Data (CSV)", style = "width:100%; margin-bottom: 5px;"),
      downloadButton("downloadSHP", "Export Spatial Grid (Zipped Shapefile)", style = "width:100%;")
      ),
    
    mainPanel(
      leafletOutput("map", height = "650px"),
      verbatimTextOutput("summaryText")
    )
  )
)

# --- DEFINE SERVER LOGIC ---
server <- function(input, output, session) {
  
  reactiveData <- reactiveValues(
    joined = NULL,
    grid_sf = NULL,
    grid = NULL,
    summary = "Ready. Set bounding box coordinates, then click '1. Fetch & Mask Land Records'."
  )
  
  # Pre-load world land polygon once for spatial clipping
  world_land <- ne_countries(scale = "medium", returnclass = "sf")
  land_union <- st_union(st_geometry(world_land))
  
  # Initial Leaflet Map
  output$map <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%
      setView(lng = 77.0, lat = 15.0, zoom = 5)
  })
  
  # STEP 1: Fetch raw OBIS data & apply land mask
  observeEvent(input$fetchBtn, {
    req(input$xmin, input$xmax, input$ymin, input$ymax)
    
    showNotification("Querying OBIS API and applying spatial land mask...", type = "message", duration = NULL, id = "fetching")
    
    tryCatch({
      wkt_geometry <- sprintf(
        "POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))",
        input$xmin, input$ymin,
        input$xmax, input$ymin,
        input$xmax, input$ymax,
        input$xmin, input$ymax,
        input$xmin, input$ymin
      )
      
      occ <- occurrence(geometry = wkt_geometry)
      
      if (is.null(occ) || nrow(occ) == 0) {
        removeNotification("fetching")
        showNotification("No OBIS records found for this region. Try a different bounding box.", type = "error")
        return()
      }
      
      occ_clean <- occ %>% filter(!is.na(decimalLongitude), !is.na(decimalLatitude), !is.na(species))
      
      if (nrow(occ_clean) == 0) {
        removeNotification("fetching")
        showNotification("Retrieved records lacked valid species/coordinates.", type = "error")
        return()
      }
      
      occ_sf <- st_as_sf(occ_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
      
      # Spatial Land Masking
      land_intersects <- st_intersects(occ_sf, land_union, sparse = FALSE)
      occ_sf_marine <- occ_sf[!land_intersects[, 1], ]
      
      if (nrow(occ_sf_marine) == 0) {
        removeNotification("fetching")
        showNotification("All fetched points fell on land. Please choose an offshore marine area.", type = "error")
        return()
      }
      
      # OBIS Standard MSP Grid: Generate 0.1-degree HEXAGONAL cells (square = FALSE)
      grid_polygon <- st_make_grid(
        st_as_sfc(st_bbox(occ_sf_marine)), 
        cellsize = c(0.25, 0.25), 
        square = FALSE, 
        what = "polygons"
      )
      grid_sf <- st_sf(grid_id = 1:length(grid_polygon), geometry = grid_polygon)
      
      joined <- st_join(occ_sf_marine, grid_sf, join = st_intersects)
      
      reactiveData$joined <- joined
      reactiveData$grid_sf <- grid_sf
      reactiveData$summary <- paste("Successfully fetched", nrow(occ_sf_marine), "marine records into hexagonal cells. Click '2. Calculate & Update Map'.")
      
      removeNotification("fetching")
      showNotification("Marine records fetched and gridded into hexagons!", type = "default")
      
    }, error = function(e) {
      removeNotification("fetching")
      showNotification(paste("Error fetching/masking data:", e$message), type = "error")
    })
  })
  
  # STEP 2: Calculate OBIS indicators locally and update map
  observeEvent(input$calcBtn, {
    req(reactiveData$joined, reactiveData$grid_sf)
    
    showNotification("Calculating OBIS biodiversity indicators locally...", type = "message", duration = NULL, id = "calculating")
    
    tryCatch({
      joined <- reactiveData$joined
      grid_sf <- reactiveData$grid_sf
      
      grid_indicators <- joined %>%
        st_drop_geometry() %>%
        group_by(grid_id) %>%
        summarise(
          total_records = n(),
          species_richness = n_distinct(species),
          shannon = calc_diversity(species)$shannon,
          simpson = calc_diversity(species)$simpson,
          evenness = calc_diversity(species)$evenness,
          es50 = calc_diversity(species)$es50,
          threatened_count = {
            cnt <- 0
            if ("category" %in% names(.)) {
              cnt <- cnt + sum(.data$category %in% c("VU", "EN", "CR", "Vulnerable", "Endangered", "Critically Endangered"), na.rm = TRUE)
            }
            if ("iucnRedListCategory" %in% names(.)) {
              cnt <- cnt + sum(.data$iucnRedListCategory %in% c("VU", "EN", "CR", "Vulnerable", "Endangered", "Critically Endangered"), na.rm = TRUE)
            }
            cnt
          }
        )
      
      final_grid <- left_join(grid_sf, grid_indicators, by = "grid_id") %>% filter(!is.na(total_records))
      reactiveData$grid <- final_grid
      
      reactiveData$summary <- paste("Calculation complete across", nrow(final_grid), "active hexagonal marine cells. Map updated.")
      
      removeNotification("calculating")
      showNotification("Indicators calculated successfully!", type = "default")
      
    }, error = function(e) {
      removeNotification("calculating")
      showNotification(paste("Error during calculation:", e$message), type = "error")
    })
  })
  
  # Render map hexagons and dynamic tooltips
  observe({
    req(reactiveData$grid, input$selected_metric)
    grid_data <- reactiveData$grid
    metric <- input$selected_metric
    
    val_vector <- grid_data[[metric]]
    valid_vals <- val_vector[!is.na(val_vector)]
    if (length(valid_vals) == 0) valid_vals <- c(0, 1)
    
    pal <- colorNumeric(palette = "viridis", domain = valid_vals, na.color = "#808080")
    
    leafletProxy("map", data = grid_data) %>%
      clearShapes() %>%
      clearControls() %>%
      addPolygons(
        fillColor = ~pal(val_vector),
        weight = 0.5,
        opacity = 0.9,
        color = "white",
        fillOpacity = 0.7,
        popup = ~paste(
          "<b>Hex ID:</b>", grid_id, "<br/>",
          "<b>Selected Metric (", metric, "):</b>", val_vector, "<br/>",
          "<hr style='margin:3px 0;'>",
          "<b>Species Richness (S):</b>", species_richness, "<br/>",
          "<b>Shannon Index (H'):</b>", shannon, "<br/>",
          "<b>Simpson Index (D):</b>", simpson, "<br/>",
          "<b>Pielou Evenness (J'):</b>", evenness, "<br/>",
          "<b>Hurlbert ES50:</b>", ifelse(is.na(es50), "N/A (<50 records)", es50), "<br/>",
          "<b>Threatened Species:</b>", threatened_count, "<br/>",
          "<b>Total Observations:</b>", total_records
        )
      ) %>%
      addLegend("bottomright", pal = pal, values = val_vector,
                title = input$selected_metric,
                opacity = 0.8)
  })
  
  output$summaryText <- renderText({
    reactiveData$summary
  })
  
  # 1. Export CSV with all metrics and centroids
  output$downloadCSV <- downloadHandler(
    filename = function() {
      paste("indobis_msp_indices_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(reactiveData$grid)
      centroids <- st_centroid(reactiveData$grid) %>% st_coordinates()
      export_df <- reactiveData$grid %>%
        st_drop_geometry() %>%
        mutate(
          longitude = centroids[,1],
          latitude = centroids[,2]
        ) %>%
        select(
          grid_id, longitude, latitude, species_richness, 
          shannon, simpson, evenness, es50, threatened_count, 
          total_records
        )
      write.csv(export_df, file, row.names = FALSE)
    }
  )
  
  # 2. Export Spatial Hexagonal Grid as a Zipped Shapefile for QGIS / ArcGIS
  output$downloadSHP <- downloadHandler(
    filename = function() {
      paste("indobis_msp_hex_grid_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      req(reactiveData$grid)
      
      temp_dir <- tempdir()
      shp_name <- "indobis_msp_hex_grid"
      shp_path <- file.path(temp_dir, paste0(shp_name, ".shp"))
      
      suppressWarnings({
        file.remove(list.files(temp_dir, pattern = shp_name, full.names = TRUE))
      })
      
      st_write(reactiveData$grid, shp_path, driver = "ESRI Shapefile", quiet = TRUE)
      
      shp_files <- list.files(temp_dir, pattern = shp_name, full.names = TRUE)
      zip::zip(zipfile = file, files = shp_files, mode = "cherry-pick")
    },
    contentType = "application/zip"
  )
}

shinyApp(ui, server)