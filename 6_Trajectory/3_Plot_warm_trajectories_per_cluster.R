#Plot of warm trajectories
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)

# -----------------------------------------------------------------------------
# 1. SETTINGS
# -----------------------------------------------------------------------------
retained_file <- "C:/Users/abraa/OneDrive/Documents/Empoasca/Wind/output/2_trajectory_warm_plots/retained_trajectories.rds"
cluster_file  <- "C:/Users/abraa/OneDrive/Documents/Empoasca/Wind/cluster_assignment_by_month_continuous.csv"
output_dir    <- "C:/Users/abraa/OneDrive/Documents/Empoasca/Wind/output/3_cluster_month_maps_lines"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


layout        <- "columns"   # "columns" = month rows x cluster columns; "overlay" = clusters together
ink_per_group <- 60          # total ink per cluster-month; larger = darker. Try 20-80
alpha_max     <- 0.35        # cap, so groups with few runs are not drawn too heavily
line_width    <- 0.25

cluster_cols <- c("1" = "#0072B2", "2" = "#D55E00")    # colour-blind-safe pair
MAP_EXTENT   <- list(xmin = -105, xmax = -55, ymin = 24, ymax = 55)
month_levels <- c("JUN", "JUL", "AUG", "SEP")
month_labels <- c("Jun", "Jul", "Aug", "Sep")

set.seed(42)   # shuffled draw order, so neither cluster always sits on top

# -----------------------------------------------------------------------------
# 2. LOAD AND PREPARE
# -----------------------------------------------------------------------------
traj <- readRDS(retained_file)

cl <- read.csv(cluster_file, stringsAsFactors = FALSE) %>%
  transmute(sample_id = ID, Month = toupper(trimws(Month)), Cluster = as.integer(Cluster))

all_ids <- unique(traj$sample_id)
cat(sprintf("Retained trajectories: %d samples | with a cluster assignment: %d\n",
            length(all_ids), sum(all_ids %in% cl$sample_id)))
cat("No assignment (not drawn):", paste(setdiff(all_ids, cl$sample_id), collapse = ", "), "\n\n")

traj <- traj %>%
  inner_join(cl, by = "sample_id") %>%
  mutate(run_key = paste(sample_id, run_id, sep = "|"),
         cl_f    = factor(Cluster),
         Month_f = factor(Month, levels = month_levels, labels = month_labels))

sample_runs  <- traj %>% distinct(sample_id, run_key) %>% count(sample_id, name = "n_runs")
group_sample <- traj %>% distinct(sample_id, Month_f, cl_f) %>% count(Month_f, cl_f, name = "n_samples_group")

traj <- traj %>%
  left_join(sample_runs,  by = "sample_id") %>%
  left_join(group_sample, by = c("Month_f", "cl_f")) %>%
  mutate(alpha_run = pmin(alpha_max, ink_per_group / (n_samples_group * n_runs)),
         run_key   = factor(run_key, levels = sample(unique(run_key)))) %>%   # shuffled draw order
  arrange(run_key, hour_along)                                               # each path in time order

group_n <- traj %>%
  distinct(sample_id, run_key, Month_f, cl_f) %>%
  group_by(Month_f, cl_f) %>%
  summarise(n_samples = n_distinct(sample_id), n_runs = n(), .groups = "drop")
print(as.data.frame(group_n))

# collection sites, snapped to 0.25 degrees so overlapping sites draw once.
# Month_f and cl_f are kept so each facet only shows the sites of ITS OWN samples
sites <- traj %>%
  filter(hour_along == 0) %>%
  mutate(lon = round(lon / 0.25) * 0.25, lat = round(lat / 0.25) * 0.25) %>%
  distinct(lon, lat, Month_f, cl_f)

# -----------------------------------------------------------------------------
# 3. MAP LAYERS
# -----------------------------------------------------------------------------
na_borders_sf <- ne_countries(scale = "large", returnclass = "sf",
                              country = c("Canada", "United States of America", "Mexico"))
us_states_sf  <- ne_states(country = "United States of America", returnclass = "sf")
ca_provs_sf   <- ne_states(country = "Canada", returnclass = "sf")

map_layers <- list(
  geom_sf(data = na_borders_sf, fill = NA, color = "grey35", linewidth = 0.35, inherit.aes = FALSE),
  geom_sf(data = us_states_sf,  fill = NA, color = "grey70", linewidth = 0.2,  inherit.aes = FALSE),
  geom_sf(data = ca_provs_sf,   fill = NA, color = "grey70", linewidth = 0.2,  inherit.aes = FALSE),
  coord_sf(xlim = c(MAP_EXTENT$xmin, MAP_EXTENT$xmax),
           ylim = c(MAP_EXTENT$ymin, MAP_EXTENT$ymax), expand = FALSE)
)

caption_txt <- paste0("Each line is one retained 80-h backward trajectory (complete, air temperature >= 9.5 C throughout). ",
                      "Each cluster has the same total ink within a month: darker = a larger share of that cluster's trajectories.")

# -----------------------------------------------------------------------------
# 4. FIGURE
# -----------------------------------------------------------------------------
p <- ggplot() + map_layers +
  geom_path(data = traj, aes(lon, lat, group = run_key, colour = cl_f, alpha = alpha_run),
            linewidth = line_width, lineend = "round") +
  scale_alpha_identity() +
  geom_point(data = sites, aes(lon, lat), shape = 21, size = 1.6, stroke = 0.7,
             fill = "white", colour = "black", inherit.aes = FALSE)

if (layout == "overlay") {
  ann <- group_n %>%
    mutate(label = sprintf("Cluster %s: %d samples, %d runs", as.character(cl_f), n_samples, n_runs),
           y     = MAP_EXTENT$ymax - 1 - (as.integer(as.character(cl_f)) - 1) * 2.2)
  p <- p +
    scale_colour_manual(values = cluster_cols, labels = function(x) paste("Cluster", x), name = NULL) +
    geom_text(data = ann, aes(x = MAP_EXTENT$xmin + 1, y = y, label = label, colour = cl_f),
              hjust = 0, vjust = 1, size = 3.2, fontface = "bold", show.legend = FALSE,
              inherit.aes = FALSE) +
    facet_wrap(~ Month_f, ncol = 2) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1.8)))
  fig_w <- 12; fig_h <- 9.5
} else {
  ann <- group_n %>% mutate(label = sprintf("%d samples, %d runs", n_samples, n_runs))
  p <- p +
    scale_colour_manual(values = cluster_cols, guide = "none") +
    facet_grid(Month_f ~ cl_f, switch = "y",
               labeller = labeller(cl_f = function(x) paste("Cluster", x)))
  fig_w <- 9; fig_h <- 12
}

p <- p +
  theme_void() +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 13),
        strip.text.y.left = element_text(face = "bold", size = 13, angle = 90),
        strip.placement = "outside")

out_file <- file.path(output_dir, sprintf("trajectory_lines_%s.png", layout))
ggsave(out_file, p, width = fig_w, height = fig_h, dpi = 300, bg = "white")
cat("\nDone:", out_file, "\n")
