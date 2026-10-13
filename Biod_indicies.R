#install.packages(c("robis", "dplyr", "sf", "ggplot2")))
library(robis)
library(dplyr)
library(sf)
library(ggplot2)

# --- STEP 1: Define a Bounding Box (e.g., a region around the Galápagos or custom sea area) ---
# Format: c(xmin, ymin, xmax, ymax)
# Example: Western Mediterranean or a coastal box
my_bbox <- c(xmin = 5, ymin = 80, xmax = 23, ymax = 94)

#cat("Fetching OBIS data for the bounding box...\n")
# Fetch occurrence records from OBIS API within the geometry
# (robis handles pagination automatically)
#occ <- occurrence(geometry = sf::st_as_sfc(sf::st_bbox(my_bbox, crs = 4326)))

# 2. Convert bounding box to an sf polygon and then to a WKT string
bbox_polygon <- st_as_sfc(st_bbox(my_bbox, crs = 4326))
wkt_geometry <- st_as_text(bbox_polygon)

# 3. Query OBIS using the WKT string format
occ <- occurrence(geometry = wkt_geometry)

# Check results
print(nrow(occ))

cat(paste("Retrieved", nrow(occ), "records from OBIS.\n"))

# Quick data cleaning: filter out records with missing coordinates or invalid species names
occ_clean <- occ %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude), !is.na(species))

# --- STEP 2: Create a Spatial Grid & Aggregate ---
# Convert occurrences to an sf spatial object
occ_sf <- st_as_sf(occ_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# Create a regular square grid over the bounding box (e.g., 0.5 degree cells)
grid_polygon <- st_make_grid(
  st_as_sfc(st_bbox(occ_sf)),
  cellsize = c(0.5, 0.5),
  what = "polygons"
)

# Turn grid into an sf dataframe with an ID
grid_sf <- st_sf(grid_id = 1:length(grid_polygon), geometry = grid_polygon)

# --- STEP 3: Calculate Biodiversity Indices per Grid Cell ---
# Spatial join: map occurrences to grid cells
joined <- st_join(occ_sf, grid_sf, join = st_intersects)

# Calculate indicators: Species Richness (S) and Total Records (Abundance)
grid_indicators <- joined %>%
  st_drop_geometry() %>%
  group_by(grid_id) %>%
  summarise(
    total_records = n(),
    species_richness = n_distinct(species),
    threatened_count = sum(category %in% c("VU", "EN"), na.rm = TRUE)
  )

# Join indicators back to the spatial grid polygons
final_grid_map <- left_join(grid_sf, grid_indicators, by = "grid_id") %>%
  filter(!is.na(total_records)) # Remove empty cells with no data

# --- STEP 4: Quick Plot ---
ggplot(final_grid_map) +
  geom_sf(aes(fill = species_richness), color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "viridis", name = "Species Richness") +
  theme_minimal() +
  labs(
    title = "OBIS Species Richness Grid for Marine Spatial Planning",
    subtitle = "Aggregated from live OBIS database queries"
  )
