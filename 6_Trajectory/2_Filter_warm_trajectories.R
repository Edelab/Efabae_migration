# =============================================================================
# For each of the candidate-day x hour trajectories (per height), this script walks
# BACKWARD from hour_along = 0 (the collection-site end) and finds the first
# point where air temperature drops below temp_threshold_C. Everything beyond
# that point is dropped -- what remains is the CONTINUOUS stretch, starting
# from the collection site, during which the parcel never got too cold for
# Empoasca to have been alive and flying' the whole way.
#
# temp_threshold_C defaults to 10 C, per your initial pass using the
# Carlson et al 1992 study's documented "ideal" transport temperature (10 C) as a
#
# Each complete warm trajectory is drawn in its own color, with position 
#markers every 3h (matching GDAS1 resolution)

##The loops in my R code were written with assistance 
#from Claude (Claude Sonnet 5.5, Anthropic, claude.ai). 
#I also used it to organize the steps, including the notes along with.

###################
library(readxl)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)

# -----------------------------------------------------------------------------
# 1. SETTINGS
# -----------------------------------------------------------------------------
traj_dir   <- "C:/Users/abraa/OneDrive/Documents/Empoasca/Wind/output/1_trajectories_per_sample"
output_dir <- "C:/Users/abraa/OneDrive/Documents/Empoasca/Wind/output/2_trajectory_warm_plots"

height_filter    <- c(150, 300, 600, 1000, 1500, 2000)
max_duration_h   <- 80      # a complete run has points at hour_along = 0, -1, ..., -80
temp_threshold_C <- 10      # threshold from the published study
temp_tolerance_C <- 0.5     # accept runs down to 9.5 C; set 0 for a strict 10 C
min_allowed_C    <- temp_threshold_C - temp_tolerance_C

sensitivity_thresholds_C <- c(10, 9.5, 9, 8, 5, 0)   # reported only, not used to filter

sample_ids_to_run <- NULL   # NULL = every sample in traj_dir, or c("A01_02JUL24", ...)
make_plots        <- TRUE

MAP_EXTENT <- list(xmin = -105, xmax = -55, ymin = 24, ymax = 55)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
plot_dir <- file.path(output_dir, "plots")
if (make_plots) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 2. MAP LAYERS (only if plotting)
# -----------------------------------------------------------------------------
if (make_plots) {
  na_borders_sf <- ne_countries(scale = "large", returnclass = "sf",
                                country = c("Canada", "United States of America", "Mexico"))
  us_states_sf  <- ne_states(country = "United States of America", returnclass = "sf")
  ca_provs_sf   <- ne_states(country = "Canada", returnclass = "sf")
}

# -----------------------------------------------------------------------------
# 3. FIND SAMPLES
# -----------------------------------------------------------------------------
sample_ids_all <- sub("_traj\\.xlsx$", "",
                      list.files(traj_dir, pattern = "_traj\\.xlsx$"))
if (length(sample_ids_all) == 0) stop("No *_traj.xlsx files found in ", traj_dir)

if (is.null(sample_ids_to_run)) {
  sample_ids <- sample_ids_all
} else {
  unmatched <- setdiff(sample_ids_to_run, sample_ids_all)
  if (length(unmatched) > 0) {
    stop("No traj file for: ", paste(unmatched, collapse = ", "))
  }
  sample_ids <- sample_ids_to_run
}
cat(sprintf("Found %d sample(s). Rule: complete %d-h runs with min temp >= %.2f C\n\n",
            length(sample_ids), max_duration_h, min_allowed_C))

# -----------------------------------------------------------------------------
# 4. PER-SAMPLE PROCESSING
# -----------------------------------------------------------------------------
process_sample <- function(sample_id) {
  
  traj <- read_excel(file.path(traj_dir, sprintf("%s_traj.xlsx", sample_id))) %>%
    filter(height_m %in% height_filter) %>%
    mutate(
      run_id     = cumsum(hour_along == 0),            # new run starts at each hour_along == 0
      air_temp_C = round(air_temp_K - 273.15, 2)
    )
  
  # sanity check: every run must begin at hour_along == 0 in row order
  bad <- traj %>% group_by(run_id) %>%
    summarise(first_ha = first(hour_along), .groups = "drop") %>%
    filter(first_ha != 0)
  if (nrow(bad) > 0) stop("rows are not ordered with hour_along = 0 first in each run")
  
  runs <- traj %>%
    group_by(run_id) %>%
    summarise(
      candidate_date = first(candidate_date),
      height_m       = first(height_m),
      start_hour     = first(hour),
      n_points       = n(),
      min_hour_along = min(hour_along),
      min_temp_C     = min(air_temp_C),
      hours_below    = sum(air_temp_C < min_allowed_C),
      .groups = "drop"
    ) %>%
    mutate(
      complete = n_points == max_duration_h + 1 & min_hour_along == -max_duration_h,
      retained = complete & min_temp_C >= min_allowed_C,
      sample_id = sample_id,
      .before = 1
    )
  
  points <- traj %>%
    semi_join(filter(runs, retained), by = "run_id") %>%
    mutate(sample_id = sample_id) %>%
    select(sample_id, run_id, candidate_date,
           any_of(c("days_before_collection")), height_m,
           any_of(c("year", "month", "day")), hour, hour_along,
           lat, lon, height, air_temp_C,
           any_of(c("rainfall", "mixdepth", "rh")))
  
  site <- traj %>% filter(hour_along == 0) %>% distinct(lat, lon) %>% slice(1)
  
  cat(sprintf("%s: %d runs | %d complete | %d retained (%.1f%%)\n",
              sample_id, nrow(runs), sum(runs$complete), sum(runs$retained),
              100 * sum(runs$retained) / nrow(runs)))
  
  list(runs = runs, points = points, site = site)
}

plot_sample <- function(sample_id, res) {
  pts <- res$points %>%
    mutate(height_f = factor(height_m, levels = height_filter,
                             labels = paste0(height_filter, " m AGL")),
           run_key  = paste(run_id),
           day      = factor(candidate_date))
  n_run <- nrow(res$runs); n_ret <- sum(res$runs$retained)
  
  p <- ggplot(pts, aes(lon, lat, group = run_key, color = day)) +
    geom_path(linewidth = 0.4, alpha = 0.85) +
    geom_point(data = filter(pts, hour_along %% 6 == 0), size = 1, alpha = 0.85) +
    geom_point(data = res$site, aes(lon, lat), inherit.aes = FALSE,
               shape = 8, size = 3, color = "black", stroke = 1.2) +
    scale_color_viridis_d(name = "Candidate day") +
    geom_sf(data = na_borders_sf, fill = NA, color = "black",  linewidth = 0.4, inherit.aes = FALSE) +
    geom_sf(data = us_states_sf,  fill = NA, color = "grey50", linewidth = 0.2, inherit.aes = FALSE) +
    geom_sf(data = ca_provs_sf,   fill = NA, color = "grey50", linewidth = 0.2, inherit.aes = FALSE) +
    coord_sf(xlim = c(MAP_EXTENT$xmin, MAP_EXTENT$xmax),
             ylim = c(MAP_EXTENT$ymin, MAP_EXTENT$ymax), expand = FALSE) +
    facet_wrap(~ height_f, ncol = 3, drop = FALSE) +
    theme_void() +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold")) +
    labs(
      title    = sprintf("%s: %d of %d runs retained", sample_id, n_ret, n_run),
      subtitle = sprintf("Complete %d-h runs with air temperature >= %.1f C throughout. Star = collection site; dots every 6 h.",
                         max_duration_h, min_allowed_C)
    )
  ggsave(file.path(plot_dir, sprintf("%s_retained_warm.png", sample_id)),
         p, width = 14, height = 11, dpi = 300, bg = "white")
}

# -----------------------------------------------------------------------------
# 5. RUN ALL SAMPLES
# -----------------------------------------------------------------------------
results <- list()
failed  <- character()

for (i in seq_along(sample_ids)) {
  sid <- sample_ids[i]
  cat(sprintf("[%d/%d] ", i, length(sample_ids)))
  tryCatch({
    res <- process_sample(sid)
    results[[sid]] <- res
    if (make_plots && nrow(res$points) > 0) plot_sample(sid, res)
  }, error = function(e) {
    failed <<- c(failed, sid)
    cat(sprintf("%s FAILED: %s\n", sid, conditionMessage(e)))
  })
}

# -----------------------------------------------------------------------------
# 6. COMBINED OUTPUTS
# -----------------------------------------------------------------------------
all_runs     <- bind_rows(lapply(results, function(r) r$runs))
retained_pts <- bind_rows(lapply(results, function(r) r$points))

sample_summary <- all_runs %>%
  group_by(sample_id) %>%
  summarise(
    n_runs            = n(),
    n_complete        = sum(complete),
    n_retained        = sum(retained),
    pct_retained      = round(100 * n_retained / n_runs, 1),
    n_days_retained   = n_distinct(candidate_date[retained]),
    n_heights_retained = n_distinct(height_m[retained]),
    .groups = "drop"
  ) %>%
  arrange(n_retained)

retained_counts_long <- all_runs %>%
  filter(retained) %>%
  count(sample_id, candidate_date, height_m, name = "n_retained")

sensitivity <- bind_rows(lapply(sensitivity_thresholds_C, function(th) {
  all_runs %>%
    group_by(sample_id) %>%
    summarise(n_retained = sum(complete & min_temp_C >= th), .groups = "drop") %>%
    mutate(threshold_C = th)
}))

write.csv(all_runs,              file.path(output_dir, "all_runs_summary.csv"),         row.names = FALSE)
write.csv(sample_summary,        file.path(output_dir, "sample_summary.csv"),           row.names = FALSE)
write.csv(retained_counts_long,  file.path(output_dir, "retained_counts_long.csv"),     row.names = FALSE)
write.csv(sensitivity,           file.path(output_dir, "sensitivity_by_threshold.csv"), row.names = FALSE)
write.csv(retained_pts,          file.path(output_dir, "retained_trajectories.csv"),    row.names = FALSE)
saveRDS(retained_pts,            file.path(output_dir, "retained_trajectories.rds"))

