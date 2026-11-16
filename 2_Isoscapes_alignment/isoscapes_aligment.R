#Final assignment test build with the aid from Claude#

setwd("C:/Users/abraa/OneDrive/Documents/Empoasca/Isotope/2assign")

#load packages
library(terra)
library(assignR)
library(rnaturalearth)

#The S file can be obtained from Dargent et al. 2026
path_Sr <- "C:/Users/abraa/OneDrive/Documents/Empoasca/Isotope/2assign/Sr.tif"  
path_S  <- "C:/Users/abraa/OneDrive/Documents/Empoasca/Isotope/2assign/S.tif"  

output_dir <- "aligned_isoscapes"
dir.create(output_dir, showWarnings = FALSE)

# Bounding box for Empoasca fabae range (WGS84)
lon_min <- -101; lon_max <- -60
lat_min <-   20; lat_max <-  55

bbox_wgs84 <- vect(data.frame(
  lon = c(lon_min, lon_max, lon_max, lon_min, lon_min),
  lat = c(lat_min, lat_min, lat_max, lat_max, lat_min)
), crs = "EPSG:4326")

# North America land mask
na_sf   <- ne_countries(scale = "large", returnclass = "sf",
                        country = c("Canada", "United States of America", "Mexico"))
na_vect <- vect(na_sf)  # WGS84

#Hydrogen data
d2h_global <- getIsoscapes(isoType = "GlobalPrecipGS", timeout = 1200)

# Extract H layers 
d2h_mean_raw <- d2h_global$d2h
d2h_sd_raw   <- d2h_global$d2h.se

# Mask to North America land, then crop to species bounding box
na_wgs84      <- project(na_vect, crs(d2h_mean_raw))
d2h_mean_mask <- mask(crop(d2h_mean_raw, na_wgs84), na_wgs84)
d2h_sd_mask   <- mask(crop(d2h_sd_raw,   na_wgs84), na_wgs84)

bbox_h        <- project(bbox_wgs84, crs(d2h_mean_raw))
d2h_mean_crop <- crop(d2h_mean_mask, bbox_h)
d2h_sd_crop   <- crop(d2h_sd_mask,   bbox_h)

# Stack into 2-layer SpatRaster (mean + SD) for calRaster()
isoscape_H_raw <- c(d2h_mean_crop, d2h_sd_crop)
names(isoscape_H_raw) <- c("d2h", "d2h.se")

#H CALIBRATION — monarch tissue regression (VSMOW_H scale)
class(knownOrig$sites)

monarch_known <- subOrigData(marker    = "d2H",
                             taxon     = "Danaus plexippus",
                             ref_scale = "VSMOW_H",
                             genplot   = FALSE)

cal_H <- calRaster(known        = monarch_known,
                   isoscape     = isoscape_H_raw,
                   interpMethod = 2,
                   genplot      = FALSE,
                   verboseLM    = TRUE)


# cal_H is a 'rescale' object — isoStack accepts this directly
cat(sprintf("  Calibration: slope=%.4f, intercept=%.4f\n",
            cal_H$lm.model$coefficients[2],
            cal_H$lm.model$coefficients[1]))

#SR AND S ISOSCAPES — load only (isoStack handles reprojection)

r_Sr <- rast(path_Sr)
r_S  <- rast(path_S)
cat(sprintf("  Sr: CRS=%s, res=%.0f m\n", crs(r_Sr, describe=TRUE)$name, res(r_Sr)[1]))
cat(sprintf("  S:  CRS=%s, res=%.0f m\n", crs(r_S,  describe=TRUE)$name, res(r_S)[1]))

#STACK AND ALIGN with isoStack()
iso_stack <- isoStack(cal_H, r_Sr, r_S, clean = TRUE)
print(iso_stack)

#check
plot(iso_stack)

#Export and align to Sulfur
iso_names <- c("H", "Sr", "S")
for (i in seq_along(iso_stack)) {
  out_path <- file.path(output_dir, paste0(iso_names[i], "_aligned.tif"))
  writeRaster(iso_stack[[i]], out_path, overwrite = TRUE)
  cat(sprintf("  Saved: %s\n", out_path))
}
