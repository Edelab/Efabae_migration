#This runs the trajectory for each sampled point
#The GDAS files can be downloaded automatic in splitr package however, 
#it failed multiple times because the connection with the server was lost.
# a safe option is to have the files from may to september saved in the computer within a folder for each year
#these files can be downloaded from 
#https://www.ready.noaa.gov/data/archives/gdas1/
#each weeek has ~ 500 MB

#We also used the latest version of HYSPLIT model. Once downloaded, the file must be open
#and indicated in the respective directory. 
#https://www.ready.noaa.gov/HYSPLIT.php

#The loops in my R code were written with assistance 
#from Claude (Claude Sonnet 5.5, Anthropic, claude.ai). 
#I also used it to organize the steps for more clear reproducibility

setwd("C:/Users/abraa/OneDrive/Documents/Empoasca/Wind")

library(splitr)
library(dplyr)
library(readxl)
library(writexl)
library(lubridate)

# -----------------------------------------------------------------------------
# 1. SETTINGS
# -----------------------------------------------------------------------------
points_xlsx   <- "Points.xlsx"
points_sheet  <- "Raw_data"

#the met folder needs to be change when running each each year
met_folder    <- "C:/Users/abraa/OneDrive/Documents/wdata/gdas1/2024"
heights       <- c(150, 300, 600, 1000, 1500, 2000)
hours         <- c(0, 3, 6, 9, 12, 15, 18, 21)
window_days   <- 6
duration_h    <- 80

binary_dir    <- "C:/HYSPLIT/exec/"   #v5.4.2

#I split per year so the computer do not 'burn'

#2024
sample_ids_to_run <- c("A01_02JUL24", "A02_15JUL24", "A03_22JUL24", 
                       "A04_06AUG24", "A05_09SEP24", "A06_16SEP24", 
                       "A07_23SEP24", "A08_17JUN24", "A09_15JUL24",
                       "A10_22JUL24", "A11_26AUG24", "A12_3SEP24", 
                       "A13_23SEP24", "A14_19JUN24", "A15_10JUL24",
                       "A16_5AUG24", "A17_5SEP24", "A18_18JUN24", 
                       "A19_09JUL24", "A20_09JUL24", "A21_18JUL24",
                       "A22_7AUG24", "A23_27AUG24", "A24_04SEP24", 
                       "A25_18SEP24")


#2025
sample_ids_to_run <- c("A26_16JUN25", "A27_16JUN25","A28_23JUN25", 
                       "A29_23JUN25", "A30_30JUN25", "A31_21JUL25",
                       "A32_21JUL25", "A33_21JUL25", "A34_04AUG25",
                       "A34_04AUG25", "A35_04AUG25", "A36_25AUG25",
                       "A37_22SEP25", "A38_24SEP25", "A39_15JUL25",
                       "A40_22JUL25", "A41_02JUL25", "A42_16JUL25",
                       "A43_05JUN25", "A44_05JUN25", "A45_12JUN25", 
                       "A46_12JUN25", "A47_19JUN25", "A48_26JUN25",
                       "A49_26JUN25", "A50_03JUL25", "A51_03JUL25",
                       "A52_07AUG25", "A53_07AUG25", "A54_07AUG25", 
                       "A55_11AUG25")


#2026
sample_ids_to_run <- c("A56_16JUN26", "A57_29JUN26","A58_06JUL26", 
                       "A59_24JUL26", "A60_27JUL26", "A61_10AUG26",
                       "A62_31AUG26", "A63_14SEP26", "A64_29JUN26",
                       "A65_10AUG26", "A66_31AUG26", "A67_14SEP26",
                       "A68_28JUL26", "A69_01SEP26", "A70_01SEP26",
                       "A71_17JUN26", "A72_22JUL26", "A73_11JUN26",
                       "A74_24JUN26", "A75_02JUL26", "A76_24JUN26", 
                       "A77_02JUL26", "A78_02SEP26", "A79_09SEP26",
                       "A80_02SEP26", "A81_09SEP26", "A82_02JUL26",
                       "A83_24JUN26", "A84_02SEP26", "A85_02SEP26", 
                       "A86_09SEP26")

#save each trajectory (excel files)
output_dir <- "output/1_trajectories_per_sample"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#This function was adopted because split R did not return meteorological data
#of the trajectories
# Reads a raw HYSPLIT extended-met tdump file directly (splitr's own reader
# assumes an old line-wrapping format the v5.4.2 binary doesn't use, and
# throws "Column `year` doesn't exist" on this newer output -- this function
# reads the same file correctly instead).

read_traj_extended <- function(path) {
  lines <- readLines(path, warn = FALSE)
  header_idx <- grep("PRESSURE", lines)[1]
  if (is.na(header_idx)) return(NULL)
  data_lines <- trimws(lines[(header_idx + 1):length(lines)])
  data_lines <- data_lines[data_lines != ""]
  fields <- strsplit(data_lines, "\\s+")
  fields <- fields[sapply(fields, length) == 22]
  if (length(fields) == 0) return(NULL)
  mat <- matrix(as.numeric(unlist(fields)), ncol = 22, byrow = TRUE)
  yr <- mat[,3]; yr <- ifelse(yr < 50, yr + 2000, yr + 1900)
  data.frame(
    year = yr, month = mat[,4], day = mat[,5], hour = mat[,6],
    hour_along = mat[,9], lat = mat[,10], lon = mat[,11], height = mat[,12],
    pressure = mat[,13], theta = mat[,14],
    air_temp_K  = mat[,15],   # Kelvin -- 10C = 283.15K
    rainfall = mat[,16], mixdepth = mat[,17], rh = mat[,18],
    sp_humidity = mat[,19], h2o_mixrate = mat[,20],
    terr_msl = mat[,21], sun_flux = mat[,22]
  )
}

# as before but reads the raw output itself instead of splitr's reader.
run_traj_extended_met <- function(lat, lon, height, duration, days, daily_hours,
                                  direction, met_dir, binary_path, exec_dir) {
  before_top   <- list.files(exec_dir, full.names = TRUE)
  before_files <- list.files(exec_dir, recursive = TRUE, full.names = TRUE)
  
  suppressWarnings(try(
    hysplit_trajectory(
      lat = lat, lon = lon, height = height, duration = duration,
      days = days, daily_hours = daily_hours, direction = direction,
      met_type = "gdas1", met_dir = met_dir, extended_met = TRUE,
      binary_path = binary_path, exec_dir = exec_dir, clean_up = FALSE
    ),
    silent = TRUE
  ))
  
  after_top   <- list.files(exec_dir, full.names = TRUE)
  after_files <- list.files(exec_dir, recursive = TRUE, full.names = TRUE)
  
  new_files <- setdiff(after_files, before_files)   # actual files (dirs excluded automatically)
  new_top   <- setdiff(after_top, before_top)       # new top-level entries
  new_dirs  <- new_top[dir.exists(new_top)]         # just the new subfolder(s), for cleanup
  
  result <- dplyr::bind_rows(lapply(new_files, read_traj_extended))
  
  file.remove(new_files[file.exists(new_files)])
  unlink(new_dirs, recursive = TRUE)   # clean up the now-empty subfolder
  result
}


# -----------------------------------------------------------------------------
# 2. LOAD SAMPLES 
# -----------------------------------------------------------------------------
samples_all <- read_excel(points_xlsx, sheet = points_sheet) %>%
  filter(!is.na(Assignment_ID)) %>%
  transmute(
    id  = Assignment_ID,
    lat = Latitude,
    lon = Longitude,
    collection_date = as_date(`Collection date`)
  )

cat(sprintf("Loaded %d samples total in %s.\n", nrow(samples_all), points_xlsx))

if (is.null(sample_ids_to_run)) {
  samples <- samples_all
} else {
  unmatched <- setdiff(sample_ids_to_run, samples_all$id)
  if (length(unmatched) > 0) {
    stop(sprintf("sample_ids_to_run contains ID(s) not found in %s: %s",
                 points_xlsx, paste(unmatched, collapse = ", ")))
  }
  samples <- samples_all %>% filter(id %in% sample_ids_to_run)
}

#check if the samples match
cat(sprintf("This run will process %d sample(s): %s\n", nrow(samples), paste(samples$id, collapse = ", ")))

# -----------------------------------------------------------------------------
# 3. RUN THE TRAJECTORY FOR THE SELECTED SAMPLE(S) - this can take a bit of time
# -----------------------------------------------------------------------------
for (i in seq_len(nrow(samples))) {
  
  s <- samples[i, ]
  out_file <- file.path(output_dir, sprintf("%s_traj.xlsx", s$id))
  
  if (file.exists(out_file)) {
    cat(sprintf("[%d/%d] %-15s already done -- skipping.\n", i, nrow(samples), s$id))
    next
  }
  
  candidate_dates <- seq(s$collection_date - window_days, s$collection_date, by = "day")
  
  cat(sprintf("[%d/%d] %-15s | %d candidate day(s) x %d height(s) x %d hour(s)\n",
              i, nrow(samples), s$id, length(candidate_dates), length(heights), length(hours)))
  
  run_list <- vector("list", length(candidate_dates) * length(heights))
  k <- 1
  
  for (d in seq_along(candidate_dates)) {
    for (h in seq_along(heights)) {
      
      traj <- tryCatch(
        run_traj_extended_met(
          lat = s$lat,
          lon = s$lon,
          height = heights[h],
          duration = duration_h,
          days = as.character(candidate_dates[d]),
          daily_hours = hours,
          direction = "backward",
          met_dir = met_folder,
          binary_path = binary_dir,
          exec_dir = getwd()          # same folder as before (your setwd() above)
        ),
        error = function(e) {
          cat(sprintf("    [error] %s | date=%s height=%d: %s\n",
                      s$id, candidate_dates[d], heights[h], conditionMessage(e)))
          NULL
        }
      )

      if (!is.null(traj)) {
        traj$Assignment_ID <- s$id
        traj$candidate_date <- as.character(candidate_dates[d])
        traj$days_before_collection <- as.numeric(s$collection_date - candidate_dates[d])
        traj$height_m <- heights[h]
        run_list[[k]] <- traj
      }
      k <- k + 1
    }
  }
  
  sample_traj <- bind_rows(run_list)
  
  is_datetime_col <- vapply(sample_traj, function(x) inherits(x, "POSIXct") || inherits(x, "Date"), logical(1))
  sample_traj[is_datetime_col] <- lapply(sample_traj[is_datetime_col], as.character)
  
  write_xlsx(sample_traj, out_file)
  cat(sprintf("    Saved: %s (%d trajectory points)\n", out_file, nrow(sample_traj)))
}

#the trajectories will be saved as a excel files
