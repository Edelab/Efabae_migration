#This code was based on the ones from Reich et al. 2025 including the 
#pallete color adopted (Vanessa_cardui)

#M. S. Reich, et al. PNAS Nexus 4, pgae586 (2025).

#The loops presented in this code was co-created using Claude and 
#i also used to organize each section of the script.
#Once we defined the putative locals, the binary geotiff files was transfer to a 
#unique file

library(terra)
library(tidyterra)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(patchwork)
library(readxl)
library(scales)




# -----------------------------------------------------------------------------
# 1. SETTINGS
# -----------------------------------------------------------------------------
base_dir   <- "C:/Users/abraa/OneDrive/Documents/Empoasca/Isotope/2assign/finalplots"
locals_dir <- file.path(base_dir, "5_Locals")            # all putative Local *_binary.tif files
points_xlsx  <- file.path(base_dir, "points.xlsx")       # collection sites (black dots)
points_sheet <- "Raw_data"
output_dir <- file.path(locals_dir, "locais")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#filter from the excel file the putative locals information
classification_xlsx <- file.path(base_dir, "local_migratory_classification.xlsx")
class_id_col    <- "ID"
class_final_col <- "Final_putative"
class_local_lab <- "Local"

#same extent since the assignments
MAP_EXTENT <- list(xmin = -105, xmax = -65, ymin = 24, ymax = 49)
months <- c("JUN", "JUL", "AUG", "SEP")                  # panel order of the per-month figure
n_cols <- 2                                              # panels per row in the per-month figure

# page sizes (mm)
total_w_mm <- 85;  total_h_mm <- 100                     # one quadrant
grid_w_mm  <- 170; grid_h_mm  <- 200                     # full page, per-month figure

# -----------------------------------------------------------------------------
# 2. MAP LAYERS
# -----------------------------------------------------------------------------
na_borders_sf <- ne_countries(scale = "large", returnclass = "sf",
                              country = c("Canada", "United States of America", "Mexico"))
us_states_sf  <- ne_states(country = "United States of America", returnclass = "sf")
ca_provs_sf   <- ne_states(country = "Canada", returnclass = "sf")
cardui_cols   <- colorRampPalette(c("#F2F2F2", "#a7d8d9", "#faab36", "#fd5900"))(10)
map_ext_terra <- ext(MAP_EXTENT$xmin, MAP_EXTENT$xmax, MAP_EXTENT$ymin, MAP_EXTENT$ymax)

points_raw <- if (file.exists(points_xlsx)) as.data.frame(read_excel(points_xlsx, sheet = points_sheet)) else NULL

# -----------------------------------------------------------------------------
# 3. FILES, IDs, MONTHS
# -----------------------------------------------------------------------------
binary_files <- list.files(locals_dir, pattern = "_binary\\.tif$", full.names = TRUE)
binary_files <- binary_files[!grepl("proportion_summed", binary_files)]
sample_ids   <- trimws(sub("_binary\\.tif$", "", basename(binary_files)))

if (!is.null(classification_xlsx) && file.exists(classification_xlsx)) {
  cl <- as.data.frame(read_excel(classification_xlsx))
  loc_ids <- trimws(cl[[class_id_col]][cl[[class_final_col]] == class_local_lab])
  drop <- !(sample_ids %in% loc_ids)
  if (any(drop)) cat("Ignored (not Local in classification):", paste(sample_ids[drop], collapse = ", "), "\n")
  absent <- setdiff(loc_ids, sample_ids)
  if (length(absent)) cat("Local in classification but no tif found:", paste(absent, collapse = ", "), "\n")
  binary_files <- binary_files[!drop]; sample_ids <- sample_ids[!drop]
}
if (length(binary_files) == 0) stop("No Local *_binary.tif files found in ", locals_dir)

id_month <- sub("^.*_[0-9]+([A-Za-z]{3})[0-9]{2}.*$", "\\1", sample_ids)
id_month <- toupper(id_month)
cat(sprintf("%d Locals found. By month: %s\n", length(sample_ids),
            paste(names(table(id_month)), table(id_month), sep = "=", collapse = ", ")))

# -----------------------------------------------------------------------------
# 4. HELPERS
# -----------------------------------------------------------------------------
# stack rasters that may sit on different grids (e.g., H-only maps): resample to the finest one
load_stack <- function(files) {
  rl  <- lapply(files, rast)
  tmpl <- rl[[which.max(sapply(rl, ncell))]]
  rl <- lapply(rl, function(r) if (isTRUE(compareGeom(r, tmpl, stopOnError = FALSE))) r else resample(r, tmpl, method = "near"))
  rast(rl)
}

prop_map <- function(files) {
  stk  <- load_stack(files)
  prop <- sum(stk, na.rm = TRUE) / nlyr(stk)
  valid <- ifel(sum(!is.na(stk)) > 0, 1, NA)
  crop(mask(prop, valid), map_ext_terra)
}

build_panel <- function(files, ids, label = NULL) {
  pr  <- prop_map(files)
  pts <- if (!is.null(points_raw))
    unique(points_raw[points_raw$Assignment_ID %in% ids, c("Latitude", "Longitude")]) else NULL
  p <- ggplot() +
    geom_spatraster(data = pr) +
    scale_fill_gradientn(colours = cardui_cols, limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
                         na.value = "transparent", name = "Proportion",
                         guide = guide_colorbar(barwidth = unit(3.5, "cm"), barheight = unit(0.2, "cm"),
                                                title.position = "left", title.vjust = 1)) +
    geom_sf(data = na_borders_sf, fill = NA, color = "black",  linewidth = 0.4, inherit.aes = FALSE) +
    geom_sf(data = us_states_sf,  fill = NA, color = "grey50", linewidth = 0.2, inherit.aes = FALSE) +
    geom_sf(data = ca_provs_sf,   fill = NA, color = "grey50", linewidth = 0.2, inherit.aes = FALSE) +
    { if (!is.null(pts) && nrow(pts) > 0)
      geom_point(data = pts, aes(x = Longitude, y = Latitude), color = "black", size = 1.5,
                 alpha = 0.6, shape = 19, inherit.aes = FALSE) } +
    coord_sf(xlim = c(MAP_EXTENT$xmin, MAP_EXTENT$xmax), ylim = c(MAP_EXTENT$ymin, MAP_EXTENT$ymax), expand = FALSE) +
    theme_void() +
    theme(legend.text = element_text(size = 6), legend.title = element_text(size = 7),
          plot.margin = margin(0, 0, 0, 0))
  if (!is.null(label))
    p <- p + annotate("label", x = MAP_EXTENT$xmin + 0.8, y = MAP_EXTENT$ymax - 0.6, label = label,
                      hjust = 0, vjust = 1, size = 2.6, fontface = "bold",
                      fill = alpha("white", 0.85), label.size = 0, label.padding = unit(0.15, "lines"))
  p
}

# -----------------------------------------------------------------------------
# 5. TOTAL MAP (all Locals, one panel, no title)
# -----------------------------------------------------------------------------
# vertical legend parked over the empty Atlantic (lower-right), clear of Florida and the coast
p_total <- build_panel(binary_files, sample_ids, label = sprintf("n = %d", length(sample_ids))) +
  scale_fill_gradientn(colours = cardui_cols, limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
                       na.value = "transparent", name = "Proportion",
                       guide = guide_colorbar(barwidth = unit(0.25, "cm"), barheight = unit(2.6, "cm"))) +
  theme(legend.position = c(0.985, 0.03), legend.justification = c(1, 0), legend.direction = "vertical",
        legend.background = element_blank(), legend.margin = margin(2, 4, 2, 4))
ggsave(file.path(output_dir, "Locals_total.png"), p_total, width = total_w_mm, height = total_h_mm,
       units = "mm", dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Locals_total.svg"), p_total, width = total_w_mm, height = total_h_mm,
       units = "mm", bg = "white")

# -----------------------------------------------------------------------------
# 6. PER-MONTH FIGURE (supplement): one panel per month, shared legend, A-D
# -----------------------------------------------------------------------------
panels <- list(); k <- 0
for (m in months) {
  sel <- which(id_month == m)
  if (length(sel) == 0) { cat(sprintf("  %s: no Locals -> panel skipped\n", m)); next }
  k <- k + 1
  panels[[k]] <- build_panel(binary_files[sel], sample_ids[sel],
                             label = sprintf("%s  %s  (n = %d)", LETTERS[k], m, length(sel)))
}
unknown <- setdiff(unique(id_month), months)
if (length(unknown)) cat("Months outside the panel list (not plotted):", paste(unknown, collapse = ", "), "\n")

fig <- wrap_plots(panels, ncol = n_cols) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave(file.path(output_dir, "Locals_by_month.png"), fig, width = grid_w_mm, height = grid_h_mm,
       units = "mm", dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Locals_by_month.svg"), fig, width = grid_w_mm, height = grid_h_mm,
       units = "mm", bg = "white")

write.csv(data.frame(ID = sample_ids, Month = id_month),
          file.path(output_dir, "locals_by_month_log.csv"), row.names = FALSE)
cat("\nDone. Files in:", output_dir, "\n")