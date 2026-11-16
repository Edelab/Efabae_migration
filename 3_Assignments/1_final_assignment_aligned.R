#This code was based on the ones from Reich et al. 2025 including the 
#pallete color adopted (Vanessa_cardui)

#M. S. Reich, et al. PNAS Nexus 4, pgae586 (2025).

#We used three isoscapes (_aligned files), and run the analyses considering
#the isotope combination for each sample. For each case, it will create 
#for each sample four files, two plots (continuous and binary maps), and
#two giff files continuous and binary posterior probability surfaces

#Assignment can take time specially in dual and triple assignments.

#The loops presented in this code was co-created using Claude and 
#i also used to organize each section of the script.

setwd("C:/Users/abraa/OneDrive/Documents/Empoasca/Isotope/2assign")
library(readxl)
library(terra)
library(assignR)
library(readxl)
library(rnaturalearth)
library(ggplot2)
library(tidyterra)
library(scico)
library(sf)
library(geosphere)   

# -----------------------------------------------------------------------------
# 1. USER SETTINGS
# -----------------------------------------------------------------------------

#In Reich paper they considered 100 km. While there is not such estimation for 
#minimum distance for local dispersal in E. fabae. 
#Here i used 25 km for this initial estimation to determine locals. 
#We will compare the minimum with enrichment ratio for the samples later
#lets see if they agree or not.

local_threshold_km <- 25
output_dir         <- "assignment_results"
dir.create(output_dir, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. LOAD DATA
# -----------------------------------------------------------------------------
isodata  <- read_excel("data.xlsx", sheet = 1)
location <- read_excel("data.xlsx", sheet = 2)

isodata$d2H <- as.numeric(isodata$d2H)
isodata$Sr  <- as.numeric(isodata$Sr)
isodata$S   <- as.numeric(isodata$S)

# Collection points as SpatVector (for wDist) and sf (for ggplot)
location_sv <- vect(location, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
location_sf <- st_as_sf(location, coords = c("Longitude", "Latitude"), crs = 4326)

# -----------------------------------------------------------------------------
# 3. LOAD ALIGNED ISOSCAPES
# -----------------------------------------------------------------------------
#This extent i took considering the main occurrence range of E. fabae
#E. fabae has been related in California and in other provinces such as Alberta
#But the main impact, following DeLong 1938 and Baker et al. 2015 
#still on this range

MAP_EXTENT <- list(
  xmin = -105,
  xmax =  -55,
  ymin =   24,
  ymax =   50 
)

domain_ext <- ext(MAP_EXTENT$xmin, MAP_EXTENT$xmax, MAP_EXTENT$ymin, MAP_EXTENT$ymax)

iso_H  <- crop(rast("H_aligned.tif"),  domain_ext)
iso_Sr <- crop(rast("Sr_aligned.tif"), domain_ext)
iso_S  <- crop(rast("S_aligned.tif"),  domain_ext)
res(iso_H)

#Check the isoscapes

par(mfrow = c(1, 3))
plot(iso_H[[1]],  main = "H (mean) — cropped extent")
plot(iso_Sr[[1]], main = "Sr (mean) — cropped extent")
plot(iso_S[[1]],  main = "S (mean) — cropped extent")
par(mfrow = c(1, 1))

# -----------------------------------------------------------------------------
# 4. MAP BORDER LAYERS
# -----------------------------------------------------------------------------

na_borders_sf <- ne_countries(scale   = "large", returnclass = "sf",
                              country = c("Canada", "United States of America", "Mexico"))
us_states_sf  <- ne_states(country = "United States of America", returnclass = "sf")
ca_provs_sf   <- ne_states(country = "Canada", returnclass = "sf")

# -----------------------------------------------------------------------------
# 5. HELPER FUNCTIONS
# -----------------------------------------------------------------------------

# -- Minimum great circle distance (km) from capture point to HP area --
min_gc_distance <- function(binary_layer, cap_lon, cap_lat) {
  hp <- as.data.frame(binary_layer, xy = TRUE, na.rm = TRUE)
  hp <- hp[hp[, 3] == 1, c("x", "y")]
  if (nrow(hp) == 0) return(NA_real_)
  dists <- distHaversine(cbind(cap_lon, cap_lat), as.matrix(hp))
  min(dists, na.rm = TRUE) / 1000
}

# -- Check if capture location falls inside high-probability area --
is_local_by_overlap <- function(binary_layer, cap_lon, cap_lat) {
  pt  <- vect(matrix(c(cap_lon, cap_lat), nrow = 1), crs = "EPSG:4326")
  val <- terra::extract(binary_layer, pt)[1, 2]
  isTRUE(val == 1)
}

# -- ggplot theme applied to all maps --
map_theme <- function() {
  theme_void() +
    theme(
      plot.title        = element_text(size = 11, face = "bold", hjust = 0.5,
                                       margin = margin(b = 6)),
      legend.position   = "bottom",
      legend.key.width  = unit(2, "cm"),
      legend.key.height = unit(0.4, "cm"),
      legend.title      = element_text(size = 9),
      legend.text       = element_text(size = 8),
      plot.margin       = margin(10, 10, 10, 10)
    )
}

# -- Map border layers list (added to every ggplot) --
map_borders <- function() {
  list(
    geom_sf(data = na_borders_sf, fill = NA, color = "black",  linewidth = 0.4,
            inherit.aes = FALSE),
    geom_sf(data = us_states_sf,  fill = NA, color = "grey50", linewidth = 0.2,
            inherit.aes = FALSE),
    geom_sf(data = ca_provs_sf,   fill = NA, color = "grey50", linewidth = 0.2,
            inherit.aes = FALSE)
  )
}

# -- coord_sf zoomed to MAP_EXTENT (applied to every map) --
map_coord <- function() {
  coord_sf(
    xlim   = c(MAP_EXTENT$xmin, MAP_EXTENT$xmax),
    ylim   = c(MAP_EXTENT$ymin, MAP_EXTENT$ymax),
    expand = FALSE
  )
}

# -----------------------------------------------------------------------------
# 6. CORE ASSIGNMENT FUNCTION
#
# Runs pdRaster(), qtlRaster(), wDist(), saves rasters, plots maps,
# and returns a distance summary dataframe — all for one combination.
#
# Arguments:
#   iso_stack  : isoStack or single SpatRaster for this combination
#   samples_df : subset of isodata for samples in this combination
#   combo_label: string e.g. "H", "H_Sr", "H_Sr_S"
# -----------------------------------------------------------------------------

run_assignment <- function(iso_stack, samples_df, combo_label) {
  
  cat(sprintf("\n========== %s assignment (%d samples) ==========\n",
              combo_label, nrow(samples_df)))
  
  combo_dir <- file.path(output_dir, combo_label)
  dir.create(combo_dir, showWarnings = FALSE)
  
  # -- 1. Run pdRaster() for all samples in this combination at once --
  prob_r   <- pdRaster(iso_stack, samples_df, genplot = FALSE)
  binary_r <- qtlRaster(prob_r, threshold = 0.333,
                        thresholdType = "prob", genplot = FALSE)
  
  # Record layer names as assigned by pdRaster (used for file naming)
  sample_ids <- names(prob_r)
  
  # -- 2. Individual binary .tif files are saved inside the per-sample loop --
  cat(sprintf("  Starting per-sample export (%d layers)\n", nlyr(binary_r)))
  
  # -- 3. wDist() — weighted distance and bearing for all samples at once --
  # sites must be a SpatVector matching the order of samples in prob_r
  sites_sv <- location_sv[match(samples_df$ID, location_sv$ID), ]
  
  wd <- tryCatch(
    wDist(prob_r, sites_sv),
    error = function(e) { warning(paste("wDist failed:", e$message)); NULL }
  )
  
  if (!is.null(wd)) {
    wdist_df        <- c(wd)
    wdist_df$Combo  <- combo_label
    write.csv(wdist_df,
              file.path(combo_dir, paste0(combo_label, "_wdist_summary.csv")),
              row.names = FALSE)
    cat(sprintf("  wDist summary saved (%d rows)\n", nrow(wdist_df)))
  }
  
  dist_results <- vector("list", length(sample_ids))
  
  # Terra extent object for cropping rasters to the map window
  map_ext_terra <- ext(MAP_EXTENT$xmin, MAP_EXTENT$xmax,
                       MAP_EXTENT$ymin, MAP_EXTENT$ymax)
  
  for (i in seq_along(sample_ids)) {
    
    sample_id  <- sample_ids[i]
    prob_layer <- prob_r[[i]]      # positional — safe regardless of layer name
    bin_layer  <- binary_r[[i]]    # positional — safe regardless of layer name
    
    # Get capture coordinates for this sample
    loc       <- location[location$ID == sample_id, ]
    cap_lon   <- loc$Longitude[1]
    cap_lat   <- loc$Latitude[1]
    sample_pt <- data.frame(Longitude = cap_lon, Latitude = cap_lat)
    
    # -- Distance metrics (computed on full-extent layers, before cropping) --
    min_dist      <- min_gc_distance(bin_layer, cap_lon, cap_lat)
    local_overlap <- is_local_by_overlap(bin_layer, cap_lon, cap_lat)
    local_dist    <- !is.na(min_dist) && min_dist < local_threshold_km
    putative_local <- local_overlap || local_dist
    
    dist_results[[i]] <- data.frame(
      ID               = sample_id,
      Combo            = combo_label,
      Capture_Lon      = round(cap_lon, 5),
      Capture_Lat      = round(cap_lat, 5),
      Min_Distance_km  = round(min_dist, 2),
      Local_by_overlap = local_overlap,
      Local_by_dist    = local_dist,
      Putative_local   = putative_local
    )
    
    prob_crop <- crop(prob_layer, map_ext_terra)
    bin_crop  <- crop(bin_layer,  map_ext_terra)
    
    # Normalized probability raster (0–1), from cropped layer
    prob_max   <- max(terra::values(prob_crop), na.rm = TRUE)
    norm_layer <- prob_crop / prob_max
    
    # Force clean 0/1 numeric raster; any other value → NA
    bin_num <- ifel(bin_crop == 1, 1, ifel(bin_crop == 0, 0, NA))
    
    # Diagnostic — confirm raw values
    uvals <- sort(unique(na.omit(as.vector(terra::values(bin_num)))))
    cat(sprintf("  [%d/%d] %s — binary unique values: {%s}\n",
                i, length(sample_ids), sample_id,
                paste(uvals, collapse = ", ")))
    
    # Save individual binary .tif for this sample (FIX 1: per-sample, not stack)
    writeRaster(bin_layer,
                file.path(combo_dir, paste0(sample_id, "_binary.tif")),
                overwrite = TRUE)
    
    #save the continuous probability raster
    writeRaster(prob_layer,
                file.path(combo_dir, paste0(sample_id, "_probability.tif")),
                overwrite = TRUE)
    
    cardui_cols <- colorRampPalette(c("#F2F2F2", "#a7d8d9", "#faab36", "#fd5900"))(10)
    
    p_norm <- ggplot() +
      geom_spatraster(data = norm_layer) +
      scale_fill_gradientn(
        colours  = cardui_cols,
        limits   = c(0, 1),
        breaks   = seq(0, 1, by = 0.2),
        na.value = "transparent",
        name     = "Proportion",
        guide    = guide_colorbar(
          title.position = "top", title.hjust = 0.5,
          barwidth       = unit(6, "cm"),
          barheight      = unit(0.4, "cm")
        )
      ) +
      map_borders() +
      geom_point(data        = sample_pt,
                 aes(x = Longitude, y = Latitude),
                 color       = "black", size = 2.5, shape = 19,
                 inherit.aes = FALSE) +
      map_coord() +
      map_theme() +
      labs(title = paste0(sample_id, "  |  ", combo_label))
    
    ggsave(file.path(combo_dir, paste0(sample_id, "_normalized.png")),
           p_norm, width = 8, height = 6, dpi = 300, bg = "white")
    
    p_bin <- ggplot() +
      geom_spatraster(data = bin_num) +
      scale_fill_gradientn(
        colours  = c("grey92", "#b03018"),
        values   = c(0, 1),
        limits   = c(0, 1),
        breaks   = c(0, 1),
        labels   = c("Low probability", "High probability\n(top 33%)"),
        na.value = "transparent",
        name     = NULL,
        guide    = guide_colorbar(
          title          = NULL,
          direction      = "horizontal",
          label.position = "bottom",
          barwidth       = unit(6, "cm"),
          barheight      = unit(0.4, "cm"),
          ticks          = FALSE
        )
      ) +
      map_borders() +
      geom_point(data        = sample_pt,
                 aes(x = Longitude, y = Latitude),
                 color       = "black", size = 2.5, shape = 19,
                 inherit.aes = FALSE) +
      map_coord() +
      map_theme() +
      labs(title = paste0(sample_id, "  |  ", combo_label, " binary (2:1 odds)"))
    
    ggsave(file.path(combo_dir, paste0(sample_id, "_binary.png")),
           p_bin, width = 8, height = 6, dpi = 300, bg = "white")
    
    cat(sprintf("             min dist: %.1f km | putative local: %s\n",
                ifelse(is.na(min_dist), NA, min_dist),
                putative_local))
    
    rm(prob_layer, bin_layer, prob_crop, bin_crop, bin_num, norm_layer,
       p_norm, p_bin)
    gc()
  }
  
  # Clean combination-level objects before returning
  rm(prob_r, binary_r)
  gc()
  
  do.call(rbind, dist_results)
}

# -----------------------------------------------------------------------------
# 7. DETECT COMBINATIONS AND SUBSET SAMPLES
# -----------------------------------------------------------------------------

isodata$combo <- paste0(
  ifelse(!is.na(isodata$d2H), "H",   ""),
  ifelse(!is.na(isodata$Sr),  "_Sr", ""),
  ifelse(!is.na(isodata$S),   "_S",  "")
)
isodata$combo <- sub("^_", "", isodata$combo)

cat("\nIsotope combinations detected:\n")
print(table(isodata$combo))

# -----------------------------------------------------------------------------
# 8. RUN ASSIGNMENTS PER COMBINATION
# -----------------------------------------------------------------------------

all_distances <- list()


# -- H only --
if (any(isodata$combo == "H")) {
  subs <- as.data.frame(isodata[isodata$combo == "H", c("ID", "d2H")])
  all_distances[["H"]] <- run_assignment(iso_H, subs, "H")
}

# -- Sr only --
if (any(isodata$combo == "Sr")) {
  subs <- as.data.frame(isodata[isodata$combo == "Sr", c("ID", "Sr")])
  all_distances[["Sr"]] <- run_assignment(iso_Sr, subs, "Sr")
}

# -- S only --
if (any(isodata$combo == "S")) {
  subs <- as.data.frame(isodata[isodata$combo == "S", c("ID", "S")])
  all_distances[["S"]] <- run_assignment(iso_S, subs, "S")
}

# -- H + Sr --
if (any(isodata$combo == "H_Sr")) {
  subs     <- as.data.frame(isodata[isodata$combo == "H_Sr", c("ID", "d2H", "Sr")])
  iso_H_Sr <- isoStack(iso_H, iso_Sr, clean = TRUE)
  all_distances[["H_Sr"]] <- run_assignment(iso_H_Sr, subs, "H_Sr")
  rm(iso_H_Sr); gc()
}

# -- H + S --
if (any(isodata$combo == "H_S")) {
  subs    <- as.data.frame(isodata[isodata$combo == "H_S", c("ID", "d2H", "S")])
  iso_H_S <- isoStack(iso_H, iso_S, clean = TRUE)
  all_distances[["H_S"]] <- run_assignment(iso_H_S, subs, "H_S")
  rm(iso_H_S); gc()
}

# -- Sr + S --
if (any(isodata$combo == "Sr_S")) {
  subs     <- as.data.frame(isodata[isodata$combo == "Sr_S", c("ID", "Sr", "S")])
  iso_Sr_S <- isoStack(iso_Sr, iso_S, clean = TRUE)
  all_distances[["Sr_S"]] <- run_assignment(iso_Sr_S, subs, "Sr_S")
  rm(iso_Sr_S); gc()
}

# -- H + Sr + S --
if (any(isodata$combo == "H_Sr_S")) {
  subs       <- as.data.frame(isodata[isodata$combo == "H_Sr_S", c("ID", "d2H", "Sr", "S")])
  iso_H_Sr_S <- isoStack(iso_H, iso_Sr, iso_S, clean = TRUE)
  all_distances[["H_Sr_S"]] <- run_assignment(iso_H_Sr_S, subs, "H_Sr_S")
  rm(iso_H_Sr_S); gc()
}

# -----------------------------------------------------------------------------
# 9. EXPORT COMBINED DISTANCE RESULTS
# -----------------------------------------------------------------------------

final_distances <- do.call(rbind, all_distances)
write.csv(final_distances,
          file.path(output_dir, "migration_distances.csv"),
          row.names = FALSE)
