#Here the idea was to calculate the raw probability that fall within
#the collection regions using all points and then creating a buffer
#area. But given the whole North America extent, we then calculate
#the ratio considering the four regions sampled as a unique with a buffer. 

setwd("C:/Users/abraa/OneDrive/Documents/Empoasca/Isotope/2assign/buffer_check")

#load packages
library(terra)
library(sf)
library(readxl)
library(ggplot2)
library(tidyterra)
library(rnaturalearth)
library(geosphere)

#Here we put all the the geotiff files from the assignment.
#for the ratio, we need the continous ones.
#the binary will be used for presentation only.

#I organize the script using the Claude AI so it will be more easier to follow

#define the directory
probability_dir <- "C:/Users/abraa/OneDrive/Documents/Empoasca/Isotope/2assign/buffer_check/probability"  
binary_dir      <- "C:/Users/abraa/OneDrive/Documents/Empoasca/Isotope/2assign/buffer_check/binary"       
points_xlsx     <- "C:/Users/abraa/OneDrive/Documents/Empoasca/Isotope/2assign/buffer_check/point.xlsx"

#create the buffer for each point for later unique zone 
buffer_km <- 10

#I define two randomly, but later we will compare what is the ratio
#between those found by the minimum distance
enrichment_threshold <- 1      

#define the output
output_dir <- "C:/Users/abraa/OneDrive/Documents/Empoasca/Isotope/2assign/buffer_check/output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#define the extent same as in the assignment
MAP_EXTENT <- list(xmin = -105, xmax = -65, ymin = 24, ymax = 49)  # match main pipeline

#load the borders from the represented countries and their states
na_borders_sf <- ne_countries(scale = "large", returnclass = "sf",
                              country = c("Canada", "United States of America", "Mexico"))
us_states_sf  <- ne_states(country = "United States of America", returnclass = "sf")
ca_provs_sf   <- ne_states(country = "Canada", returnclass = "sf")

#ggplot2 set
map_theme <- function() {
  theme_void() +
    theme(
      plot.title      = element_text(size = 11, face = "bold", hjust = 0.5, margin = margin(b = 6)),
      legend.position = "bottom",
      plot.margin     = margin(10, 10, 10, 10)
    )
}

map_borders <- function() {
  list(
    geom_sf(data = na_borders_sf, fill = NA, color = "black",  linewidth = 0.4, inherit.aes = FALSE),
    geom_sf(data = us_states_sf,  fill = NA, color = "grey50", linewidth = 0.2, inherit.aes = FALSE),
    geom_sf(data = ca_provs_sf,   fill = NA, color = "grey50", linewidth = 0.2, inherit.aes = FALSE)
  )
}

map_coord <- function() {
  coord_sf(xlim = c(MAP_EXTENT$xmin, MAP_EXTENT$xmax),
           ylim = c(MAP_EXTENT$ymin, MAP_EXTENT$ymax), expand = FALSE)
}

#So since we sampled the same field multiple times, this will filter 
#unique points (i.e. if i point is repeated numerically it will be removed)

#build the reference zone based on the unique sample points
points_raw  <- as.data.frame(read_excel(points_xlsx))
unique_pts  <- unique(points_raw[, c("Latitude", "Longitude")])
cat(sprintf("Loaded %d unique collection points.\n", nrow(unique_pts)))
pts_sf <- st_as_sf(unique_pts, coords = c("Longitude", "Latitude"), crs = 4326)

# Project to a metric CRS — EPSG:3347 (Statistics Canada Lambert) avoids
# UTM zone-edge distortion since the sites span two UTM zones.
pts_proj <- st_transform(pts_sf, crs = 3347)

#Create a single single convex hull spanning ALL collection points together
hull_all        <- st_convex_hull(st_union(pts_proj))
local_zone_proj <- st_buffer(hull_all, dist = buffer_km * 1000)
local_zone_sf   <- st_transform(st_sf(geometry = local_zone_proj), crs = 4326)
local_zone_sv   <- vect(local_zone_sf)

st_write(local_zone_sf, file.path(output_dir, "local_reference_zone.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)
cat(sprintf("Local reference zone (one convex hull + %.0f km margin) saved: %s/local_reference_zone.gpkg\n",
            buffer_km, output_dir))

#zone zoom in
zone_bbox <- st_bbox(local_zone_sf)
zoom_pad_x <- as.numeric(zone_bbox["xmax"] - zone_bbox["xmin"]) * 0.25
zoom_pad_y <- as.numeric(zone_bbox["ymax"] - zone_bbox["ymin"]) * 0.25
zoom_xlim  <- c(zone_bbox["xmin"] - zoom_pad_x, zone_bbox["xmax"] + zoom_pad_x)
zoom_ylim  <- c(zone_bbox["ymin"] - zoom_pad_y, zone_bbox["ymax"] + zoom_pad_y)

p_zonein <- ggplot() +
  map_borders() +
  geom_sf(data = local_zone_sf, fill = "steelblue", alpha = 0.35, color = "steelblue4", linewidth = 0.5) +
  geom_point(data = unique_pts, aes(x = Longitude, y = Latitude),
             color = "black", size = 1.6, shape = 19) +
  coord_sf(xlim = zoom_xlim, ylim = zoom_ylim, expand = FALSE) +
  map_theme() +
  labs(title = sprintf("Local reference zone, n = %d sites",
                       buffer_km, nrow(unique_pts)))
p_zonein

ggsave(file.path(output_dir, "local.png"),
       p_zonein, width = 8, height = 6, dpi = 300, bg = "white")
cat(sprintf("Saved: %s/local_reference_zone_map.png\n", output_dir))

p_zone <- ggplot() +
  map_borders() +
  geom_sf(data = local_zone_sf, fill = "steelblue", alpha = 0.35, color = "steelblue4", linewidth = 0.5) +
  geom_point(data = unique_pts, aes(x = Longitude, y = Latitude),
             color = "black", size = 1.6, shape = 19) +
  map_coord() +
  map_theme() +
  labs(title = sprintf("Local reference zone (one convex hull + %.0f km margin, n = %d sites)",
                       buffer_km, nrow(unique_pts)))

ggsave(file.path(output_dir, "local_reference_zone_map.png"),
       p_zone, width = 8, height = 6, dpi = 300, bg = "white")
cat(sprintf("Saved: %s/local_reference_zone_map.png\n", output_dir))

#classify each sample
prob_files <- list.files(probability_dir, pattern = "_probability\\.tif$", full.names = TRUE)
cat(sprintf("\nFound %d probability rasters in %s\n", length(prob_files), probability_dir))
if (length(prob_files) == 0) stop("No *_probability.tif files found — check probability_dir path (see save_probability_raster_patch.R).")

# Zone area (km^2), computed once in a metric CRS
zone_area_km2 <- as.numeric(sum(st_area(st_transform(local_zone_sf, 3347)))) / 1e6

# Domain area (km^2), computed once from the first raster (same grid for all samples)
first_r <- rast(prob_files[1])
domain_area_km2 <- as.numeric(global(mask(cellSize(first_r, unit = "km"), first_r), "sum", na.rm = TRUE)[1, 1])
zone_area_fraction <- zone_area_km2 / domain_area_km2

cat(sprintf("Zone area: %.1f km^2 | Domain area: %.0f km^2 | Zone = %.4f%% of domain area\n",
            zone_area_km2, domain_area_km2, 100 * zone_area_fraction))

classify_sample <- function(prob_tif_path, zone_sv, zone_area_fraction, enrichment_threshold) {
  
  r <- rast(prob_tif_path)
  zone_sv_r <- zone_sv
  if (!identical(crs(zone_sv_r), crs(r))) crs(zone_sv_r) <- crs(r)
  total_mass <- as.numeric(global(r, "sum", na.rm = TRUE)[1, 1])
  zone_mass  <- as.numeric(global(mask(r, zone_sv_r), "sum", na.rm = TRUE)[1, 1])
  mass_fraction_in_zone <- if (is.finite(total_mass) && total_mass > 0) zone_mass / total_mass else NA_real_
  enrichment <- if (is.finite(mass_fraction_in_zone) && zone_area_fraction > 0) {
    mass_fraction_in_zone / zone_area_fraction
  } else NA_real_
  is_local_A <- isTRUE(is.finite(enrichment) && enrichment >= enrichment_threshold)
  
  list(mass_fraction_in_zone = mass_fraction_in_zone, enrichment = enrichment, is_local_A = is_local_A)
}

results <- vector("list", length(prob_files))

for (i in seq_along(prob_files)) {
  f <- prob_files[i]
  sample_id <- trimws(sub("_probability\\.tif$", "", basename(f)))
  
  cls <- classify_sample(f, local_zone_sv, zone_area_fraction, enrichment_threshold)
  
  results[[i]] <- data.frame(
    ID                     = sample_id,
    Mass_fraction_in_zone  = round(cls$mass_fraction_in_zone, 5),
    Enrichment_ratio       = round(cls$enrichment, 2),
    Classification_A       = ifelse(cls$is_local_A, "Local", "Migratory")
  )
  
  cat(sprintf("[%d/%d] %-20s | A: enrich=%.2fx -> %s\n",
              i, length(prob_files), sample_id,
              cls$enrichment, ifelse(cls$is_local_A, "LOCAL", "MIGRATORY")))
}

results_df <- do.call(rbind, results)
write.csv(results_df, file.path(output_dir, "local_migratory_classification.csv"), row.names = FALSE)


