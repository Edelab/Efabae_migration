#This script was based on Reich et al. 2025, includuing the pallete in the paper for Vanessa cardui.

#M. S. Reich, et al. PNAS Nexus 4, pgae586 (2025).

#I used Claude AI (Sonnet 5.5) to organize the sections to facilitate
#reproduction


# =============================================================================

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
base_dir <- "C:/Users/abraa/OneDrive/Documents/Empoasca/Isotope/2assign/finalplots"
cluster_csv <- file.path(base_dir, "cluster_assignment_all_transforms.csv")
id_col      <- "ID"
month_col   <- "Month"
cluster_col <- "Cluster_sqrt"   #cluster by square root mean provide the best silhoute

month_map   <- c("JUN" = "1_Jun", "JUL" = "2_Jul", "AUG" = "3_Aug", "SEP" = "4_Sep")
month_order <- c("JUN", "JUL", "AUG", "SEP")   # keep in this order to avoid aug being on the top

pool_months       <- c()               

swap_months       <- c()               
align_columns     <- TRUE              
n_pattern_cols    <- NULL              
mean_surface_dir  <- file.path(base_dir, sub("^Cluster_", "", cluster_col), "cluster_mean_surfaces")

bearing_xlsx   <- file.path(base_dir, "migration_distance33.xlsx")
bearing_id_col <- "Sample_ID"
bearing_col    <- "wMeanBear"

points_xlsx  <- file.path(base_dir, "points.xlsx")
points_sheet <- "Raw_data"

MAP_EXTENT <- list(xmin = -105, xmax = -65, ymin = 24, ymax = 49)

# rose inset position/size inside each map panel 
rose_box <- list(left = 0.62, bottom = 0.01, right = 0.98, top = 0.45)

rose_bin_width <- 3
rose_start     <- 0
rose_direction <- 1

#max size from PNAS guidelines
target_width_mm  <- 170
target_height_mm <- 200   

#output 
output_dir <- file.path(base_dir, "cluster_final")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 2. SHARED LAYERS + DATA
# -----------------------------------------------------------------------------
na_borders_sf <- ne_countries(scale = "large", returnclass = "sf",
                              country = c("Canada", "United States of America", "Mexico"))
us_states_sf  <- ne_states(country = "United States of America", returnclass = "sf")
ca_provs_sf   <- ne_states(country = "Canada", returnclass = "sf")
cardui_cols   <- colorRampPalette(c("#F2F2F2", "#a7d8d9", "#faab36", "#fd5900"))(10)
map_ext_terra <- ext(MAP_EXTENT$xmin, MAP_EXTENT$xmax, MAP_EXTENT$ymin, MAP_EXTENT$ymax)

points_raw <- if (file.exists(points_xlsx)) as.data.frame(read_excel(points_xlsx, sheet = points_sheet)) else NULL
cluster_df <- read.csv(cluster_csv, stringsAsFactors = FALSE)
bearing_df <- as.data.frame(read_excel(bearing_xlsx))

# ---- units to draw: one per (month, cluster); a pooled month is a single unit ----
units <- list()
for (m in month_order) {
  rows <- cluster_df[cluster_df[[month_col]] == m, ]
  if (nrow(rows) == 0) next
  if (m %in% pool_months) {
    units[[length(units) + 1]] <- list(month = m, orig = sort(unique(rows[[cluster_col]])), ids = rows[[id_col]])
  } else {
    for (k in sort(unique(rows[[cluster_col]])))
      units[[length(units) + 1]] <- list(month = m, orig = k, ids = rows[[id_col]][rows[[cluster_col]] == k])
  }
}
unit_month <- sapply(units, function(u) u$month)
unit_n     <- sapply(units, function(u) length(u$ids))
K <- max(table(unit_month))
cat("Clusters drawn per month:", paste(names(table(unit_month)), as.integer(table(unit_month)), sep = "=", collapse = ", "), "\n")

# ---- mean surface of a unit (weighted by cluster size for a pooled month) ----
read_unit_surface <- function(u) {
  rows <- cluster_df[cluster_df[[month_col]] == u$month, ]
  ws <- sapply(u$orig, function(k) sum(rows[[cluster_col]] == k))
  fs <- file.path(mean_surface_dir, sprintf("%s_cluster%d_mean.tif", u$month, u$orig))
  if (any(!file.exists(fs))) return(NULL)
  rl <- lapply(fs, rast)
  Reduce(`+`, Map(function(x, w) x * w, rl, ws)) / sum(ws)
}

# ---- column of each unit ----
col_of <- integer(length(units))
apply_swaps <- function(col_of, ncols) {                       # mirror the columns of the months in swap_months
  for (m in swap_months) { ui <- which(unit_month == m); col_of[ui] <- ncols + 1 - col_of[ui] }
  col_of
}
G <- if (is.null(n_pattern_cols)) K else n_pattern_cols
G <- min(max(G, K), length(units))      # never fewer patterns than clusters in the busiest month
header_prefix <- "Cluster"

if (align_columns) {
  surfs <- lapply(units, read_unit_surface)
  if (any(sapply(surfs, is.null))) {
    cat("Note: cluster-mean surfaces not found in", mean_surface_dir, "- columns by cluster size instead.\n")
    align_columns <- FALSE
  }
}

if (align_columns) {
  tmpl  <- surfs[[1]]
  surfs <- lapply(surfs, function(x) if (compareGeom(x, tmpl, stopOnError = FALSE)) x else resample(x, tmpl, method = "bilinear"))
  M  <- do.call(cbind, lapply(surfs, function(x) as.numeric(values(x))))
  ok <- complete.cases(M)
  M  <- M[ok, , drop = FALSE]
  cors    <- cor(M, method = "spearman")
  lat     <- yFromCell(tmpl, which(ok))
  cen_lat <- apply(M, 2, function(p) sum(p * lat) / sum(p))
  
  grp <- cutree(hclust(as.dist(1 - cors), method = "average"), k = G)
  grp_lat  <- tapply(cen_lat, grp, mean)
  grp_rank <- match(grp, as.integer(names(sort(grp_lat))))      # 1 = southernmost pattern
  
  assign_month <- function(S) {                                 # best one-to-one match of clusters to patterns
    nm <- nrow(S); Gc <- ncol(S)
    grid <- expand.grid(rep(list(seq_len(Gc)), nm))
    grid <- grid[apply(grid, 1, function(r) !anyDuplicated(r)), , drop = FALSE]
    best <- which.max(apply(grid, 1, function(r) sum(S[cbind(seq_len(nm), r)])))
    as.integer(unlist(grid[best, ]))
  }
  for (m in unique(unit_month)) {
    ui <- which(unit_month == m)
    S <- sapply(seq_len(G), function(g) sapply(ui, function(i) {
      mem <- setdiff(which(grp_rank == g), i)
      if (length(mem) == 0) 0 else mean(cors[i, mem])
    }))
    S <- matrix(S, nrow = length(ui), ncol = G)
    col_of[ui] <- assign_month(S)
  }
  col_of <- apply_swaps(col_of, G)
  header_prefix <- "Pattern"
  assign_tab <- data.frame(Month = unit_month, n = unit_n,
                           Original_cluster = sapply(units, function(u) paste(u$orig, collapse = "+")),
                           Column = col_of, Centroid_lat = round(cen_lat, 1))
  cat("\nColumn assignment (patterns ordered south -> north):\n"); print(assign_tab, row.names = FALSE)
  cat("\nSpearman correlation between cluster-mean surfaces:\n")
  cm <- round(cors, 2); dimnames(cm) <- list(paste(unit_month, sapply(units, function(u) paste(u$orig, collapse = "+")), sep = "_"),
                                             paste(unit_month, sapply(units, function(u) paste(u$orig, collapse = "+")), sep = "_"))
  print(cm)
  write.csv(assign_tab, file.path(output_dir, sprintf("column_assignment_%s.csv", cluster_col)), row.names = FALSE)
  write.csv(cm, file.path(output_dir, sprintf("cluster_surface_correlation_%s.csv", cluster_col)))
} else {
  for (m in unique(unit_month)) {                               # columns by size, largest first
    ui <- which(unit_month == m)
    col_of[ui] <- rank(-unit_n[ui], ties.method = "first")
  }
  col_of <- apply_swaps(col_of, K)
}
Kc <- max(col_of)

# -----------------------------------------------------------------------------
# 3. HELPERS
# -----------------------------------------------------------------------------
make_rose_plot <- function(bearings, n_label) {
  bearings <- bearings[!is.na(bearings)] %% 360
  df <- data.frame(bearing = bearings, ind = seq_along(bearings), height = 1)
  df$slot <- floor(df$bearing / rose_bin_width) * rose_bin_width + rose_bin_width / 2
  r_max   <- max(table(df$slot), 1)
  label_r <- r_max * 1.25
  cardinal_df <- data.frame(x = c(0, 90, 180, 270), y = label_r, lab = c("N", "E", "S", "W"))
  ggplot(df, aes(x = slot, y = height, group = ind)) +
    geom_col(position = "stack", fill = "#a9c9ef", color = "grey35", linewidth = 0.15, width = rose_bin_width * 0.9) +
    geom_text(data = cardinal_df, aes(x = x, y = y, label = lab), inherit.aes = FALSE,
              size = 2.2, fontface = "bold", color = "grey25") +
    coord_polar(theta = "x", start = rose_start, direction = rose_direction, clip = "off") +
    scale_x_continuous(limits = c(0, 360), breaks = c(0, 90, 180, 270), expand = c(0, 0)) +   # only N/E/S/W spokes, no diagonals
    scale_y_continuous(limits = c(0, r_max * 1.35), breaks = scales::pretty_breaks(n = 2), expand = c(0, 0)) +   # fewer concentric rings
    theme_void(base_size = 6) +
    theme(
      panel.grid.major = element_line(color = "grey85", linewidth = 0.12),
      plot.background  = element_rect(fill = NA, color = NA),
      plot.caption     = element_text(size = 7.5, hjust = 0.5, margin = margin(t = 1)),
      plot.margin      = margin(2, 2, 0, 2)
    ) +
    labs(caption = sprintf("n = %d", n_label))   # n= sits BELOW the rose, like the reference figure
}

build_panel <- function(month_label, ids, tag_letter) {
  folder_name <- month_map[[month_label]]
  month_dir   <- file.path(base_dir, folder_name)
  
  bin_all <- list.files(month_dir, pattern = "_binary\\.tif$", recursive = TRUE, full.names = TRUE)
  bin_ids <- sub("_binary\\.tif$", "", basename(bin_all))
  binary_files <- bin_all[match(ids, bin_ids)]
  found <- !is.na(binary_files)
  if (any(!found)) cat(sprintf("  %s: no binary tif for %s\n", month_label, paste(ids[!found], collapse = ", ")))
  ids_found    <- ids[found]
  binary_files <- binary_files[found]
  n_ind <- length(binary_files)
  
  group_pts <- if (!is.null(points_raw)) {
    unique(points_raw[points_raw$Assignment_ID %in% ids_found, c("Latitude", "Longitude")])
  } else NULL
  
  # binary tifs may sit on slightly different grids (e.g. H-only) -> bring to one template
  rl   <- lapply(binary_files, rast)
  tmpl <- rl[[which.max(sapply(rl, function(x) prod(dim(x)[1:2])))]]
  bad  <- !sapply(rl, function(x) compareGeom(x, tmpl, stopOnError = FALSE))
  if (any(bad)) rl[bad] <- lapply(rl[bad], function(x) resample(x, tmpl, method = "near"))
  bin_stack <- rast(rl)
  sum_r  <- sum(bin_stack, na.rm = TRUE)
  prop_r <- ifel(sum(!is.na(bin_stack)) > 0, sum_r / n_ind, NA)
  prop_crop <- crop(prop_r, map_ext_terra)
  
  p_map <- ggplot() +
    # grey land background FIRST (bottom layer) so it only shows through
    # where the raster is transparent/NA -- it must NOT be drawn after the
    # raster, or its solid fill paints right over the colored data
    geom_sf(data = na_borders_sf, fill = "grey96", color = NA, inherit.aes = FALSE) +
    geom_spatraster(data = prop_crop) +
    scale_fill_gradientn(colours = cardui_cols, limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
                         na.value = "transparent",
                         name = "Proportion of assignment probability",
                         guide = guide_colorbar(title.position = "bottom", title.hjust = 0.5,
                                                barwidth = unit(6, "cm"), barheight = unit(0.3, "cm"))) +
    # borders redrawn on top, unfilled, so they stay crisp over the raster
    geom_sf(data = na_borders_sf, fill = NA, color = "black",  linewidth = 0.3, inherit.aes = FALSE) +
    geom_sf(data = us_states_sf,  fill = NA, color = "grey50", linewidth = 0.15, inherit.aes = FALSE) +
    geom_sf(data = ca_provs_sf,   fill = NA, color = "grey50", linewidth = 0.15, inherit.aes = FALSE) +
    { if (!is.null(group_pts)) geom_point(data = group_pts, aes(x = Longitude, y = Latitude), color = "black", alpha = 0.6, size = 1.3, shape = 19, inherit.aes = FALSE) } +
    coord_sf(xlim = c(MAP_EXTENT$xmin, MAP_EXTENT$xmax), ylim = c(MAP_EXTENT$ymin, MAP_EXTENT$ymax), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(t = 11, r = 1, b = 1, l = 1),
          plot.tag = element_text(size = 10, face = "bold"),
          plot.tag.position = "topleft",
          plot.tag.location = "margin",   # renders the tag in the margin, OUTSIDE the panel border, not overlapping the map (requires ggplot2 >= 3.5.0)
          panel.border = element_blank()) +   # no box around the panel
    labs(tag = tag_letter)
  
  match_idx <- match(ids_found, bearing_df[[bearing_id_col]])
  bearings  <- bearing_df[[bearing_col]][match_idx]
  n_matched <- sum(!is.na(bearings))
  
  if (n_matched > 0) {
    p_rose <- make_rose_plot(bearings, n_matched)
    # bottom-right inset (moved from top-right)
    p_map <- p_map + inset_element(p_rose, left = rose_box$left, bottom = rose_box$bottom, right = rose_box$right, top = rose_box$top, align_to = "panel")
  }
  
  p_map
}

blank_plot <- function() ggplot() + theme_void()
text_plot  <- function(label, size = 5, angle = 0, face = "plain") {
  ggplot() +
    annotate("text", x = 0, y = 0, label = label, size = size, angle = angle, fontface = face) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

# -----------------------------------------------------------------------------
# 4. BUILD THE PANELS
# -----------------------------------------------------------------------------
panels <- list(); tag_i <- 0
for (m in month_order) {
  ui <- which(unit_month == m)
  for (i in ui[order(col_of[ui])]) {                            # letters follow reading order
    tag_i <- tag_i + 1
    panels[[paste0(m, "_", col_of[i])]] <- build_panel(m, units[[i]]$ids, LETTERS[tag_i])
  }
}

# -----------------------------------------------------------------------------
# 5. ASSEMBLE THE GRID: header row + 4 (label | K maps) rows
# -----------------------------------------------------------------------------
top_row <- Reduce("|", c(list(blank_plot()), lapply(seq_len(Kc), function(j) text_plot(paste(header_prefix, j), size = 5)))) +
  plot_layout(widths = c(0.05, rep(1, Kc)))

row_list <- lapply(month_order, function(m) {
  cells <- lapply(seq_len(Kc), function(j) {
    key <- paste0(m, "_", j)
    if (key %in% names(panels)) panels[[key]] else blank_plot()
  })
  Reduce("|", c(list(text_plot(m, size = 4.5, angle = 90)), cells)) + plot_layout(widths = c(0.05, rep(1, Kc)))
})

full_fig <- Reduce("/", c(list(top_row), row_list)) +
  plot_layout(heights = c(0.12, rep(1, length(month_order))), guides = "collect") &
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.box.just = "center",
        legend.title = element_text(size = 8, family = "Arial", colour = "black"),
        legend.text  = element_text(size = 7, family = "Arial", colour = "black"),
        legend.key.width = unit(1.1, "cm"),
        legend.key.height = unit(0.3, "cm"))

full_fig <- full_fig +
  plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1)))   # trims the outer blank margin around the whole assembled figure

# -----------------------------------------------------------------------------
# 6. SAVE
# -----------------------------------------------------------------------------
out_png <- file.path(output_dir, sprintf("isotope_figure_grid_%s.png", cluster_col))
out_svg <- file.path(output_dir, sprintf("isotope_figure_grid_%s.svg", cluster_col))
ggsave(out_png, full_fig, width = target_width_mm, height = target_height_mm, units = "mm", dpi = 300, bg = "white")
ggsave(out_svg, full_fig, width = target_width_mm, height = target_height_mm, units = "mm", bg = "white")
