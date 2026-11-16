#Sr isoscapes for assignment using 87Sr/86Sr plant isoscape 

#Strontium data was compiled from Reich et al. 2021 and Le Corre et al. 2025
#the other two is not available yet as manuscripts still in preparation

#Script was adapted from Reich et al. 2021 and Bataille et al. 2020. 

#References
#M. S. Reich et al. Methods Ecol Evol 12, 2445–2457 (2021).
#M. Le Corre, et al.FACETS 10, 1–17 (2025).
#C. P. Bataille, et al. Palaeogeography, Palaeoclimatology, Palaeoecology 555, 109849 (2020)


#list of packages----
library(parallel)
library(doParallel)
library(terra)
library(tidyverse)
library(sf)
library(viridis)
library(randomForest)
library(readxl)
library(rnaturalearth)
library(rnaturalearthdata)
library(caret)
library(VSURF)                                              
library(ranger)
library(corrplot)
library(raster)
library(mapmisc)
library(patchwork)
library(ggplot2)
library(pdp)
library(caret)
library(ggpubr)
library(ggspatial)

#set my directory and create folders
mydir<-setwd("C:/Users/abraa/OneDrive/Documents/Empoasca/Isotope/Stront")

os_path <- function(folder_path,filename){
  if (os=="win") return(paste(folder_path,"/",filename,sep="")) else
    return(paste(folder_path,"/",filename,sep="")) } 

Fig_Path <- "Figure"  
Out_Path <- "Output"
Ras_Path <- "Projected_rasters/Projected_rasters" #should already exist
OutRas_Path <- "Projected_rasters/Projected_rasters_NA"

#Step 1----

#import the Sr Plants data
sr_orig <- readxl::read_excel("Sr_Plants_Isotope_Database.xlsx",col_names = TRUE, na = "NA", sheet = "Sheet2")

sr_orig <- as.data.frame(sr_orig)

#Clean the training dataset removing missing values in latitude, longitude and 87Sr/86Sr
sr_orig_0 <- sr_orig[!is.na(sr_orig$collection_decimal_longitude),] #remove rows with missing XY coordinates
sr_orig_1 <- sr_orig_0[!is.na(sr_orig_0$collection_decimal_latitude),]
sr_orig_2 <- sr_orig_1[!is.na(sr_orig_1$'87Sr_86Sr'),]

target_crs_string <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

# Make spatial files of the isotope data
sr_proj <- st_as_sf(sr_orig_2, coords = c("collection_decimal_longitude", "collection_decimal_latitude"), crs = st_crs(4326)) %>% 
  st_transform(crs = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
sr_proj$collection_decimal_longitude <- sr_orig_2$collection_decimal_longitude
sr_proj$collection_decimal_latitude <- sr_orig_2$collection_decimal_latitude

#Plot the Sr points in the dataset


canada_sf <- ne_states(country = "Canada", returnclass = "sf")
usa_sf <- ne_states(country = "United States of America", returnclass = "sf")
mex_sf <- ne_states(country = "Mexico", returnclass = "sf")


sr_proj <- sr_proj %>%
  rename(Sr_ratio = `87Sr_86Sr`)

mapSr <- ggplot() +
  geom_sf(data = canada_sf, 
          fill = "white", color = "gray50", size = 0.2) +
  geom_sf(data = usa_sf, 
          fill = "white", color = "gray50", size = 0.2) +
  geom_sf(data = mex_sf, 
          fill = "white", color = "gray50", size = 0.2) +
  geom_point(data = sr_proj,
             aes(x = collection_decimal_longitude,
                 y = collection_decimal_latitude,
                 fill = Sr_ratio),
             size = 2, shape = 21, alpha = 0.8) +
  coord_sf(xlim = c(-160, -50), ylim = c(20, 70)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(16) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "right",
        legend.text = element_text(size = 16),
        panel.background = element_rect(fill = "#deebf7", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  annotation_scale(location = "br")

mapSr

#from this dataset I noticed one point in Alaska and other in the Island

#Since my focus is on the Contiguous USA Part I will remove those points
sr_proj <- sr_proj %>%
  filter(collection_decimal_longitude >= -140 & collection_decimal_longitude <= -50) %>%
  filter(collection_decimal_latitude >= 20 & collection_decimal_latitude <= 80)


mapSr2 <- ggplot() +
  geom_sf(data = canada_sf, 
          fill = "white", color = "gray50", size = 0.2) +
  geom_sf(data = usa_sf, 
          fill = "white", color = "gray50", size = 0.2) +
  geom_sf(data = mex_sf, 
          fill = "white", color = "gray50", size = 0.2) +
    geom_point(data = sr_proj,
             aes(x = collection_decimal_longitude,
                 y = collection_decimal_latitude,
                 fill = Sr_ratio),
             size = 2, shape = 21, alpha = 0.8) +
  coord_sf(xlim = c(-130, -40), ylim = c(25, 70)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(16) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "right",
        legend.text = element_text(size = 16),
        panel.background = element_rect(fill = "#deebf7", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  annotation_scale(location = "br")

mapSr2

breaks <- seq(floor(min(sr_proj$Sr_ratio)*1000)/1000, 
              max(sr_proj$Sr_ratio) + 0.003, 
              by = 0.003)

table(cut(sr_proj$Sr_ratio, breaks = breaks, include.lowest = TRUE), useNA = "ifany")

max(breaks)
max(sr_proj$Sr_ratio)

mapSr2 <- ggplot() +
  geom_sf(data = canada_sf, 
          fill = "white", color = "gray50", size = 0.2) +
  geom_sf(data = usa_sf, 
          fill = "white", color = "gray50", size = 0.2) +
  geom_sf(data = mex_sf, 
          fill = "white", color = "gray50", size = 0.2) +
  geom_point(data = sr_proj,
             aes(x = collection_decimal_longitude,
                 y = collection_decimal_latitude,
                 fill = cut(Sr_ratio, breaks = breaks, include.lowest = TRUE)),  # Cut into bins
             size = 2, shape = 21, alpha = 0.8) +
  coord_sf(xlim = c(-130, -40), ylim = c(25, 70)) +
  scale_fill_viridis_d(
    name = "Sr Ratio",
    drop=TRUE,
    option = "H",
    na.value = "grey50",
    labels = function(x) {
      gsub("\\(|\\]", "", x) %>%
        gsub(",", "-", .)
    }
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(16) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "right",
        legend.text = element_text(size = 14),  # Slightly smaller for longer labels
        panel.background = element_rect(fill = "#deebf7", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  annotation_scale(location = "br")

mapSr2

ggsave('mapSr2.png', width = 240, height = 200, units = "mm", dpi = 300)


breaks <- seq(floor(min(sr_proj$Sr_ratio)*1000)/1000, 
              ceiling(max(sr_proj$Sr_ratio)*1000)/1000, 
              by = 0.003)

# see if any values fall outside the break range
range(sr_proj$Sr_ratio)
range(breaks)

# see exactly which rows turn NA after cut()
test_cut <- cut(sr_proj$Sr_ratio, breaks = breaks)
sr_proj[is.na(test_cut), ]


#Step 2---- 

#The layers used in this study can are free for download. Please refer to the 
#supporting information for that. They couldnt be deposited given the large size
#of the total files (~11 GB)

#Crop to North America and keep layers in the same resolution (1 km)

# Define file paths
files <- list.files(path = file.path(mydir, Ras_Path), pattern = "\\.tif$", full.names = TRUE)
outfiles <- file.path(mydir, OutRas_Path, paste0("p_", basename(files)))

# Define the extent
efalim <- ext(-14000000, -2000000, 2000000, 8374837)
template_ext <- ext(efalim)

template_raster <- rast(
  crs = target_crs_string,
  resolution = 1000,
  extent = template_ext 
)

#fixed to force the same extent
for (i in seq_along(files)) {
  r <- rast(files[i])
  r1 <- crop(r, efalim, snap = "near")
  r1_projected <- project(r1, template_raster)
  r2 <- resample(r1_projected, template_raster, method = "bilinear")
  writeRaster(r2, outfiles[i], filetype = "GTiff", overwrite = TRUE)
}


#now load the raster files from the NA folder
#lets check the first and the last layer
#if the resolution is correct and the coordinate system


###Global median, Q1 and Q3 bedrock model 87Sr/86Sr (Bataille et al. 2020)
r.m1<-terra::rast("Projected_rasters/Projected_rasters_NA/p_rm1_reproj.tif")
plot(r.m1)
crs(r.m1)
res(r.m1)

r.srsrq1<-terra::rast("Projected_rasters/Projected_rasters_NA/p_srsrq1.tif")
plot(r.srsrq1)

r.srsrq3<-terra::rast("Projected_rasters/Projected_rasters_NA/p_srsrq3.tif")
plot(r.srsrq3)


###Global mean, min and max geological age in Myrs from the GLiM Hartman and 
#Moosdorf 2012 and Bataille et al. 2020 for text to numerica age field converstion
r.maxage_geol<-terra::rast("Projected_rasters/Projected_rasters_NA/p_agemax.tif")
plot(r.maxage_geol)

r.meanage_geol<-terra::rast("Projected_rasters/Projected_rasters_NA/p_agemean.tif")
plot(r.meanage_geol)

r.minage_geol<-terra::rast("Projected_rasters/Projected_rasters_NA/p_agemin.tif")
plot(r.minage_geol)


###Global mean of terrane age in Myrs (Mooney et al. 1988)
r.age<-terra::rast("Projected_rasters/Projected_rasters_NA/p_basement_age_reproj.tif")
plot(r.age)

###Global unconsolidated sediments type (see publications for each class) 
#(Börker et al., 2018) 
r.GUM<-terra::rast("Projected_rasters/Projected_rasters_NA/p_gum_mask3.tif")
plot(r.GUM)

###Global bouger gravity anomaly map (Balmino et al., 2012)
r.bouger <- terra::rast("Projected_rasters/Projected_rasters_NA/p_bouger_reproj.tif")
plot(r.bouger)

###Global elevation in m ((Jarvis et al., 2008) )
r.elevation<-terra::rast("Projected_rasters/Projected_rasters_NA/p_elevation_reproj.tif")
plot(r.elevation)

###Mean annual temperature (°C) (Harris et al., 2020) 
r.mat <- terra::rast("Projected_rasters/Projected_rasters_NA/p_mat_reproj.tif")
plot(r.mat)

### Mean annual precipitation (mm/yr) (Harris et al., 2020) 
r.map<-terra::rast("Projected_rasters/Projected_rasters_NA/p_map_reproj.tif")
plot(r.map)

###Global Aridity Index Aridity Index <-  Mean annual precipitation÷ Mean annual evapotranspiration (Trabuco & Zomer, 2019) 
r.ai<-terra::rast("Projected_rasters/Projected_rasters_NA/p_ai_reproj.tif")
plot(r.ai)

### Global Potential Evapo-Transpiration (mm day-1) (Trabuco & Zomer, 2019) 
r.pet<-terra::rast("Projected_rasters/Projected_rasters_NA/p_pet_reproj.tif")
plot(r.pet)

###World clim data
r.wind  <- terra::rast("Projected_rasters/Projected_rasters_NA/p_r.wind.tif")
plot(r.wind)

r.vapor <- terra::rast("Projected_rasters/Projected_rasters_NA/p_r.vapor.tif")
plot(r.vapor)

r.tmin  <- terra::rast("Projected_rasters/Projected_rasters_NA/p_r.tmin.tif")
plot(r.tmin)

r.tmax  <- terra::rast("Projected_rasters/Projected_rasters_NA/p_r.tmax.tif")
plot(r.tmax)

r.solar <- terra::rast("Projected_rasters/Projected_rasters_NA/p_r.solar.tif")
plot(r.solar)

###Multi-models average of dust deposition (g/m²/year) (Chien et al., 2016; Brahney et al., 2015) 
r.dust<-terra::rast("Projected_rasters/Projected_rasters_NA/p_dust_reproj.tif")
plot(r.dust)

###Simulation of sea salt aerosol deposition (g/m2/year) (Chien et al., 2016) 
r.salt<-terra::rast("Projected_rasters/Projected_rasters_NA/p_salt_reproj.tif")
plot(r.salt)

###Simulation of distance
r.distance <- terra::rast("Projected_rasters/Projected_rasters_NA/p_distance.tif")
plot(r.distance)

###Volcanos (S) (kg/m2/s) (Brahney et al., 2015) 
r.volc <- terra::rast("Projected_rasters/Projected_rasters_NA/p_volc_dep_reproj.tif")
plot(r.volc)

###Black carbon deposition (kg/m2/s) (Chien et al., 2016) 
r.fire <- terra::rast("Projected_rasters/Projected_rasters_NA/p_fire_dep_reproj.tif")
plot(r.fire)

### Fossil Fuel (Black Carbon) (kg/m2/s) (Chien et al., 2016) 
r.foss <- terra::rast("Projected_rasters/Projected_rasters_NA/p_fossilfuels_dep_reproj.tif")
plot(r.foss)

###Clay fraction (%), ph, Cation Exchange Capacity (cmol+ /kg) , bulk density, organic carbon content (Hengl et al., 2021) 
r.clay<-terra::rast("Projected_rasters/Projected_rasters_NA/p_rclay_reproj.tif")
plot(r.clay)

r.ph<-terra::rast("Projected_rasters/Projected_rasters_NA/p_rph_reproj.tif")
plot(r.ph)

r.cec<-terra::rast("Projected_rasters/Projected_rasters_NA/p_rcec_reproj.tif")
plot(r.cec)

r.bulk<-terra::rast("Projected_rasters/Projected_rasters_NA/p_rbulk_reproj.tif")
plot(r.bulk)

r.ocs<-terra::rast("Projected_rasters/Projected_rasters_NA/p_rocr_reproj.tif")
plot(r.ocs)

#Global Nitrogen fertilization 
r.nfert<-terra::rast("Projected_rasters/Projected_rasters_NA/p_nfert_reproj.tif")
plot(r.nfert)

r.w0010 <- terra::rast("Projected_rasters/Projected_rasters_NA/p_r.w0010.tif")
plot(r.w0010)

r.w0033 <- terra::rast("Projected_rasters/Projected_rasters_NA/p_r.w0033.tif")
plot(r.w0033)

r.w1500 <- terra::rast("Projected_rasters/Projected_rasters_NA/p_r.w1500.tif")
plot(r.w1500)

r.w0010_30 <- terra::rast("Projected_rasters/Projected_rasters_NA/p_r.w0010_30.tif")
plot(r.w0010_30)

###Soil organic content, organic carbon density and nitrogen content between 5-15cm in soil (Hengl et al., 2023) 
r.soc <- terra::rast("Projected_rasters/Projected_rasters_NA/p_r.soc.tif")
plot(r.soc)

r.ocd <- terra::rast("Projected_rasters/Projected_rasters_NA/p_r.ocd.tif")
plot(r.ocd)

r.nitrogen <- terra::rast("Projected_rasters/Projected_rasters_NA/p_r.nitrogen.tif")
plot(r.nitrogen)

###Permafrost cover and mean annual ground temperature (Obu et al., 2021) 
r.perma <- terra::rast("Projected_rasters/Projected_rasters_NA/p_r.perma.tif")
plot(r.perma)

r.magt  <- terra::rast("Projected_rasters/Projected_rasters_NA/p_r.magt.tif")
plot(r.magt)


###Global soil moisture map(Guevara et al., 2020) 
r.sw <- terra::rast("Projected_rasters/Projected_rasters_NA/p_r.sw.tif")
plot(r.sw)

###Global wetland(Fluet-Chouinard et al., 2023) 
r.wetland <- terra::rast("Projected_rasters/Projected_rasters_NA/p_r.wetland.tif") 
plot(r.wetland)

#GLiM 2nd lithological class attribute
r.litho<-terra::rast("Projected_rasters/Projected_rasters_NA/p_litho.tif")
plot(r.litho)

#Sand (weight %)
r.sand<-terra::rast("Projected_rasters/Projected_rasters_NA/p_sand_reproj.tif")
plot(r.sand)

#Sand (weight %)
r.silt<-terra::rast("Projected_rasters/Projected_rasters_NA/p_silt_reproj.tif")
plot(r.silt)

#Soil thickness (m)
r.soilthick<-terra::rast("Projected_rasters/Projected_rasters_NA/p_soilthick_reproj.tif")
plot(r.soilthick)

#Multi-model average sea salt wet+dry deposition (kg/ha/yr)
r.ssa<-terra::rast("Projected_rasters/Projected_rasters_NA/p_ssa.tif")
plot(r.ssa)

#Multi-model average sea salt wet deposition (kg/ha/yr)
r.ssaw<-terra::rast("Projected_rasters/Projected_rasters_NA/p_ssaw.tif")
plot(r.ssaw)

#GLiM 1st lithological class attribute
r.xx<-terra::rast("Projected_rasters/Projected_rasters_NA/p_xx.tif")
plot(r.xx)
crs(r.xx)
res(r.xx)


#Total variables = 49 variables

#Step 3----


#Extract the variables values for the random forest model
s_proj <- sr_proj

#load the function extract function
source("Extract_function.R")

#following the function
agexy <- extract_nearest_non_na(s_proj, r.age)
aixy<- extract_nearest_non_na(s_proj, r.ai)
bougerxy<- extract_nearest_non_na(s_proj, r.bouger)
bulkxy<- extract_nearest_non_na(s_proj, r.bulk)
cecxy<- extract_nearest_non_na(s_proj, r.cec)
clayxy<- extract_nearest_non_na(s_proj, r.clay)
distancexy<- extract_nearest_non_na(s_proj, r.distance)
dustxy<- extract_nearest_non_na(s_proj, r.dust)
elevationxy<- extract_nearest_non_na(s_proj, r.elevation)
firexy<- extract_nearest_non_na(s_proj, r.fire)
fossxy<- extract_nearest_non_na(s_proj, r.foss)
GUMxy<- extract_nearest_non_na(s_proj, r.GUM)
lithoxy<- extract_nearest_non_na(s_proj, r.litho)
m1xy<- extract_nearest_non_na(s_proj, r.m1)
magtxy<- extract_nearest_non_na(s_proj, r.magt)
mapxy<- extract_nearest_non_na(s_proj, r.map)
matxy<- extract_nearest_non_na(s_proj, r.mat)
maxage_geolxy<- extract_nearest_non_na(s_proj, r.maxage_geol)
meanage_geolxy<- extract_nearest_non_na(s_proj, r.meanage_geol)
minage_geolxy<- extract_nearest_non_na(s_proj, r.minage_geol)
nfertxy<- extract_nearest_non_na(s_proj, r.nfert)
nitrogenxy<- extract_nearest_non_na(s_proj, r.nitrogen)
ocdxy<- extract_nearest_non_na(s_proj, r.ocd)
ocsxy<- extract_nearest_non_na(s_proj, r.ocs)
permaxy<- extract_nearest_non_na(s_proj, r.perma)
petxy<- extract_nearest_non_na(s_proj, r.pet)
phxy<- extract_nearest_non_na(s_proj, r.ph)
saltxy<- extract_nearest_non_na(s_proj, r.salt)
sandxy<- extract_nearest_non_na(s_proj, r.sand)
siltxy<- extract_nearest_non_na(s_proj, r.silt)
socxy<- extract_nearest_non_na(s_proj, r.soc)
soilthickxy<- extract_nearest_non_na(s_proj, r.soilthick)
solarxy<- extract_nearest_non_na(s_proj, r.solar)
srsrq1xy<- extract_nearest_non_na(s_proj, r.srsrq1)
srsrq3xy<- extract_nearest_non_na(s_proj, r.srsrq3)
ssaxy<- extract_nearest_non_na(s_proj, r.ssa)
ssawxy<- extract_nearest_non_na(s_proj, r.ssaw)
swxy<- extract_nearest_non_na(s_proj, r.sw)
tmaxxy<- extract_nearest_non_na(s_proj, r.tmax)
tminxy<- extract_nearest_non_na(s_proj, r.tmin)
vaporxy<- extract_nearest_non_na(s_proj, r.vapor)
volcxy<- extract_nearest_non_na(s_proj, r.volc)
w0010xy<- extract_nearest_non_na(s_proj, r.w0010)
w0010_30xy<- extract_nearest_non_na(s_proj, r.w0010_30)
w0033xy<- extract_nearest_non_na(s_proj, r.w0033)
w1500xy<-extract_nearest_non_na(s_proj, r.w1500)
wetlandxy<- extract_nearest_non_na(s_proj, r.wetland)
windxy<- extract_nearest_non_na(s_proj, r.wind)
xxxy<- extract_nearest_non_na(s_proj, r.xx)

s_xy<- data.frame(st_coordinates(s_proj))

s_proj_xy <- data.frame(s_proj$collection_decimal_latitude, s_proj$collection_decimal_longitude, 
                        s_proj$'87Sr_86Sr',s_xy, agexy , aixy, bougerxy, bulkxy, cecxy, clayxy, distancexy, dustxy, elevationxy, firexy, fossxy, GUMxy, lithoxy,
                        m1xy, magtxy, mapxy, matxy, maxage_geolxy, meanage_geolxy, minage_geolxy, nfertxy, nitrogenxy, ocdxy,ocsxy,
                        permaxy, petxy,  phxy, saltxy, sandxy, siltxy, socxy, soilthickxy, solarxy, srsrq1xy, srsrq3xy, ssaxy, ssawxy, 
                        swxy, tmaxxy, tminxy, vaporxy, volcxy, w0010xy, w0010_30xy, w0033xy,w1500xy,wetlandxy,  windxy,  xxxy)

colnames(s_proj_xy)<-c('Latitude','Longitude','X87Sr_86Sr',
                       "X","Y",'r.age', 'r.ai', 'r.bouger', 'r.bulk', 'r.cec', 'r.clay', 'r.distance', 'r.dust', 
                       'r.elevation', 'r.fire', 'r.foss', 'r.GUM', 'r.litho',
                       'r.m1', 'r.magt', 'r.map', 'r.mat', 'r.maxage_geol', 'r.meanage_geol', 'r.minage_geol', 
                       'r.nfert', 'r.nitrogen', 'r.ocd','r.ocs',
                       'r.perma', 'r.pet',  'r.ph', 'r.salt', 'r.sand', 'r.silt', 'r.soc', 'r.soilthick', 'r.solar', 
                       'r.srsrq1', 'r.srsrq3', 'r.ssa', 'r.ssaw', 
                       'r.sw', 'r.tmax', 'r.tmin', 'r.vapor', 'r.volc', 'r.w0010', 'r.w0010_30', 
                       'r.w0033','r.w1500','r.wetland',  'r.wind',  'r.xx')




#Save
write.table(s_proj_xy, "s_proj_xy.txt", sep = "\t",col.names = T, row.names = F, quote = F)

#load
sr_proj_xy <-read.table("s_proj_xy.txt",sep="\t",h=T)

# remove NA values cases if any from the dataset
sr_proj_xy <- sr_proj_xy[!is.na(sr_proj_xy$r.age),]

# Optional. Filter 'sr_proj_xy' to exclude rows with super high or zero values
sr_proj_xy1 <- sr_proj_xy[sr_proj_xy$X87Sr_86Sr <= 0.8, ]
sr_proj_xy2 <- sr_proj_xy1[!(sr_proj_xy1$X87Sr_86Sr == 0), ]

#no changes in my dataset
#so far 2657 obs.

#Aggregate all samples by latitude and longitude to avoid having multiple samples per site
sr_agg1 <- aggregate(sr_proj_xy2,by = list(sr_proj_xy2$Latitude,sr_proj_xy2$Longitude), FUN = median, na.rm = TRUE)
sr_agg1 <- sr_agg1[,5:ncol(sr_agg1)]

#now 1104 obs
write.table(sr_agg1, "sr_agg1.txt", sep = "\t",col.names = T, row.names = F, quote = F)

#check correlation between covariates and remove redundant variables
cor <- round(cor(sr_agg1[,-(1:3)]),2)
data2 <- sr_agg1[,-(1:3)] 

# Compute correlation matrix
cor_matrix <- cor(data2, use = "complete.obs")

high_corr <- findCorrelation(cor_matrix, cutoff = 0.9, names = TRUE, exact = TRUE)
high_corr


#13 variables presented > 0.90
#"r.pet"          "r.tmax"         "r.solar"        "r.mat"         
#"r.magt"         "r.tmin"         "r.volc"         "r.maxage_geol" 
# "r.meanage_geol" "r.ssaw"         "r.ssa"          "r.sand"        
# "r.m1"


#remove it from the dataset
sr_agg <- sr_agg1[, !names(sr_agg1) %in% high_corr]

ncol(sr_agg1)
ncol(sr_agg)  
setdiff(names(sr_agg1), names(sr_agg))  

# Remove incomplete cases with NA in covariates
sr_agg <- sr_agg[complete.cases(sr_agg),]

# Save table on directory
write.table(sr_agg,"sr_agg_all.txt",sep = "\t",col.names = T, row.names = F, quote = F)

#load
sr_agg <- read.table("sr_agg_all.txt", sep = "\t", h = T)

eck4_crs <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

# convert sr_agg's Eckert IV X/Y back to WGS84 lon/lat
sr_agg_sf <- st_as_sf(sr_agg, coords = c("X", "Y"), crs = eck4_crs, remove = FALSE)
sr_agg_wgs84 <- st_transform(sr_agg_sf, crs = 4326)

# extract the lon/lat as plain columns
coords <- st_coordinates(sr_agg_wgs84)
sr_agg$lon <- coords[, 1]
sr_agg$lat <- coords[, 2]

canada_sf <- ne_states(country = "Canada", returnclass = "sf")
usa_sf    <- ne_states(country = "United States of America", returnclass = "sf")
mex_sf    <- ne_states(country = "Mexico", returnclass = "sf")

breaks <- seq(floor(min(sr_agg$X87Sr_86Sr)*1000)/1000, 
              max(sr_agg$X87Sr_86Sr) + 0.003, 
              by = 0.003)

sr_agg$Sr_bin <- cut(sr_agg$X87Sr_86Sr, breaks = breaks, include.lowest = TRUE)
sum(is.na(sr_agg$Sr_bin))  # should be 0

mapSr_training <- ggplot() +
  geom_sf(data = canada_sf, fill = "white", color = "gray50", size = 0.2) +
  geom_sf(data = usa_sf,    fill = "white", color = "gray50", size = 0.2) +
  geom_sf(data = mex_sf,    fill = "white", color = "gray50", size = 0.2) +
  geom_point(data = sr_agg,
             aes(x = lon, y = lat, fill = Sr_bin),
             size = 2, shape = 21, alpha = 0.8) +
  coord_sf(xlim = c(-130, -40), ylim = c(25, 70)) +
  scale_fill_viridis_d(
    name = "Sr Ratio",
    drop = TRUE,
    option = "H",
    na.value = "grey50",
    labels = function(x) {
      gsub("\\(|\\]", "", x) %>%
        gsub(",", "-", .)
    }
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(16) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "right",
        legend.text = element_text(size = 14),
        panel.background = element_rect(fill = "#deebf7", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  annotation_scale(location = "br")

mapSr_training

ggsave('mapSr_training.png', width = 240, height = 200, units = "mm", dpi = 300)


#Step 4----

#Variable selection and RF model

set.seed(1)

sr_agg0 <- sr_agg

sr_agg1.vsurf <- VSURF(
  sr_agg0[, 4:ncol(sr_agg0)],
  sr_agg0$X87Sr_86Sr,
  mtry = max(floor(ncol(sr_agg0)/3), 1),
  ntree.thres = 3000,
  nfor.thres = 20,
  nmin = 1,
  ntree.interp = 3000,
  nfor.interp = 10,
  nsd = 1,
  ntree.pred = 3000,
  nfor.pred = 10,
  nmj = 1,
  RFimplem = "ranger",
  parallel = TRUE,
  ncores = detectCores() - 1,
  clusterType = "PSOCK"
)

sr_agg1.vsurf$varselect.pred

#Plot VSURF results
plot(sr_agg1.vsurf)

#subset only the selected
sr_agg1.sub <- sr_agg0[, 4:ncol(sr_agg0)]
sr_agg1_VSURF <- sr_agg1.sub[c(sr_agg1.vsurf$varselect.pred)]

write.table(sr_agg1_VSURF,"sr_agg1_VSURF.txt",sep = "\t",col.names = T, row.names = F, quote = F)

sr_agg1_VSURF <- read.table("sr_agg1_VSURF.txt", sep = "\t", h = T)

sr_agg1_VSURF <- cbind(sr_agg[, 1:3], sr_agg1_VSURF)

training_rf <- sr_agg1_VSURF[, -c(2, 3)]

cluster <- parallel::makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)

fitControl <- trainControl(
  method = "repeatedcv",  # 10-fold Crossvalidation
  number = 10,
  repeats = 5,  # Repeated five times
  verboseIter = FALSE,
  returnResamp = "final",
  savePredictions = "all",
  allowParallel = TRUE
)

bestmtry <- tuneRF(training_rf, training_rf$X87Sr_86Sr, stepFactor = 1, 
                   improve = 1e-7, ntree = 3000)
mtry <- bestmtry[1]

tunegrid <- expand.grid(.mtry = mtry)

metric <- "Accuracy"

names(training_rf)

set.seed(1)


#first model with the 7 predictors

RF_all <- train(X87Sr_86Sr ~ r.w0010+r.srsrq3+r.dust+r.ocs+r.sw+r.w0010_30+r.cec,
                data = training_rf, ntree = 3000, method = "rf", importance=TRUE, 
                tuneGrid=tunegrid,trControl= fitControl)

save(RF_all, file = "RF_all")

RF_all

load("RF_all")
summary(RF_all)
print(RF_all)


#General model
p1 <- ggplot(RF_all$pred, aes(x = pred, y = obs)) +
  geom_point(shape = 21, size = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.001))+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001))+
  labs(y = expression(""^{87}*"Sr/"^{86}*"Sr observed"))+
  labs(x = expression(""^{87}*"Sr/"^{86}*"Sr predicted"))+
  theme_bw(12)+
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))
        

p1

ggsave("p1.tiff", width = 120, height = 80, units = "mm",dpi=600, compression = "lzw")
ggsave("p1.pdf", width = 120, height = 80, units = "mm",dpi=600)


#Ten fold cross-validation

p2 <- ggplot(RF_all$pred, aes(x = obs, y = pred)) +
  geom_bin_2d(binwidth = 0.0005) +
  geom_abline(slope = 1, intercept = 0, linetype = 3) +
  geom_smooth(method = "lm", se=FALSE, colour = "black") +
  scale_fill_viridis_c(option = "viridis", direction = -1, trans = "sqrt") +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.001))+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001))+
  labs(x = expression(""^{87}*"Sr/"^{86}*"Sr observed"))+
  labs(y = expression(""^{87}*"Sr/"^{86}*"Sr predicted"))+
  labs(fill = "Count") +
  theme_bw(12)+
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

p2

ggsave("p2.tiff", width = 120, height = 80, units = "mm",dpi=600, compression = "lzw")
ggsave("p2.pdf", width = 120, height = 80, units = "mm",dpi=600)


#Residuals

p3 <- ggplot(RF_all$pred, aes(x = obs, y = pred - obs)) +
  geom_bin_2d(binwidth = 0.0005, show.legend = FALSE) +
  scale_fill_viridis_c(option = "viridis", direction = -1, trans = "sqrt") +
  labs(x = expression(""^{87}*"Sr/"^{86}*"Sr observed"))+
  scale_x_continuous(labels = scales::label_number(accuracy = 0.001))+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001))+
  ylab("Residuals")+
  theme_bw(12)+
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

p3

ggsave("p3.tiff", width = 120, height = 80, units = "mm",dpi=600, compression = "lzw")
ggsave("p3.pdf", width = 120, height = 80, units = "mm",dpi=600)


#variable importance
var_imp <- varImp(RF_all, scale = TRUE)

imp_df <- as.data.frame(var_imp$importance)

rf_model <- RF_all$finalModel

varImpPlot(rf_model)

nodeplot_data <- data.frame(randomForest::importance(rf_model))
nodeplot_data$Variable <- row.names(nodeplot_data)
purity_importance <- nodeplot_data[, c("Variable", "IncNodePurity")]
purity_importance <- purity_importance[order(purity_importance$IncNodePurity, decreasing = FALSE), ]
purity_importance$Variable <- factor(purity_importance$Variable, 
                                     levels = purity_importance$Variable)


p4 <- ggplot(purity_importance, aes(x = IncNodePurity, y = Variable)) +
  geom_point(stat = "identity",color="blue", fill="blue", size = 2.5, shape=21) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.0001))+
  labs(title = "Variable importance",
    x = "Increase in node purity",
    y = "Predictor") +
    theme_bw(12)+
    theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))


p4

ggsave("p4.tiff", width = 130, height = 80, units = "mm",dpi=600, compression = "lzw")
ggsave("p4.pdf", width = 130, height = 80, units = "mm",dpi=600)


# List of predictors to plot
predictors <- c('r.w0010', 'r.srsrq3', 'r.dust', 
                'r.ocs', 'r.sw', 'r.w0010_30', 'r.cec')

pdp_1r.w0010 <- partial(RF_all, pred.var = "r.w0010")
pdp_2r.dust<- partial(RF_all, pred.var = "r.dust")
pdp_3r.srsrq3<- partial(RF_all, pred.var = "r.srsrq3")
pdp_4r.cec<- partial(RF_all, pred.var = "r.cec")
pdp_5r.sw<- partial(RF_all, pred.var = "r.sw")
pdp_6r.w0010_30<- partial(RF_all, pred.var = "r.w0010_30")
pdp_7r.ocs<- partial(RF_all, pred.var = "r.ocs")


# Plot the partial dependence for each predictor
pdp1<-autoplot(pdp_1r.w0010) + 
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1))+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.0001))+
  labs(y = expression(""^{87}*"Sr/"^{86}*"Sr Predicted"))+
  labs(x = 'Soil water content at -10 kPa (5-15 cm)')+
  theme_bw(12)+
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))
  
pdp1

ggsave(path=Out_Path, filename= "pdp1.tiff", width = 120, height = 80, units = "mm",dpi=600, compression = "lzw")


pdp2<-autoplot(pdp_2r.dust) + 
  scale_x_continuous(labels = scales::label_number(accuracy = 0.00000001))+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.0001))+
  labs(y = expression(""^{87}*"Sr/"^{86}*"Sr Predicted"))+
  labs(x = expression("Dust deposition (g.m"^{2}*".year)"))+
  theme_bw(12)+
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

pdp2

ggsave(path=Out_Path, filename= "pdp2.tiff", width = 120, height = 80, units = "mm",dpi=600, compression = "lzw")

pdp3<-autoplot(pdp_3r.srsrq3) + 
  scale_x_continuous(labels = scales::label_number(accuracy = 0.0001))+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.0001))+
  labs(y = expression(""^{87}*"Sr/"^{86}*"Sr Predicted"))+
  labs(x = expression("Q3 bedrock model"^{87}*"Sr/"^{86}*"Sr"))+
  theme_bw(12)+
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

pdp3

ggsave(path=Out_Path, filename= "pdp3.tiff", width = 120, height = 80, units = "mm",dpi=600, compression = "lzw")

pdp4<-autoplot(pdp_4r.cec) + 
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1))+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.0001))+
  labs(y = expression(""^{87}*"Sr/"^{86}*"Sr Predicted"))+
  labs(x = "Cation exchange capacity (cmol/kg)")+
  theme_bw(12)+
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

pdp4

ggsave(path=Out_Path, filename= "pdp4.tiff", width = 120, height = 80, units = "mm",dpi=600, compression = "lzw")


pdp5<-autoplot(pdp_5r.sw) + 
  #scale_x_continuous(labels = scales::label_number(accuracy = 0.001))+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.0001))+
  labs(y = expression(""^{87}*"Sr/"^{86}*"Sr Predicted"))+
  labs(x = expression("Soil moisture m"^3*"/m"^3))+
  theme_bw(12)+
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))
pdp5

ggsave(path=Out_Path, filename= "pdp5.tiff", width = 120, height = 80, units = "mm",dpi=600, compression = "lzw")


pdp6<-autoplot(pdp_6r.w0010_30) + 
  #scale_x_continuous(labels = scales::label_number(accuracy = 0.001))+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.0001))+
  labs(y = expression(""^{87}*"Sr/"^{86}*"Sr Predicted"))+
  labs(x = "Soil water capacity at −10 kPa (15-30 cm)")+
  theme_bw(12)+
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

pdp6

ggsave(path=Out_Path, filename= "pdp6.tiff", width = 120, height = 80, units = "mm",dpi=600, compression = "lzw")

pdp7<-autoplot(pdp_7r.ocs) + 
  #scale_x_continuous(labels = scales::label_number(accuracy = 0.001))+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.0001))+
  labs(y = expression(""^{87}*"Sr/"^{86}*"Sr Predicted"))+
  labs(x = "Organic carbon stock (0-30 cm) t/ha")+
  theme_bw(12)+
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

pdp7

ggsave(path=Out_Path, filename= "pdp7.tiff", width = 120, height = 80, units = "mm",dpi=600, compression = "lzw")

part10<-ggarrange(pdp1, pdp2, pdp3, pdp4,
         pdp5, pdp6, pdp7,
         ncol=2, nrow=4)

part10

ggsave(path=Out_Path, filename= "part10.tiff", width = 210, height = 220, units = "mm",dpi=600, compression = "lzw")

#load the predictors in case they are not loaded
r.srsrq3<-terra::rast("Projected_rasters/Projected_rasters_NA/p_srsrq3.tif")
r.w0010 <- terra::rast("Projected_rasters/Projected_rasters_NA/p_r.w0010.tif")
r.w0010_30 <- terra::rast("Projected_rasters/Projected_rasters_NA/p_r.w0010_30.tif")
r.ocs<-terra::rast("Projected_rasters/Projected_rasters_NA/p_rocr_reproj.tif")
r.sw <- terra::rast("Projected_rasters/Projected_rasters_NA/p_r.sw.tif")
r.cec<-terra::rast("Projected_rasters/Projected_rasters_NA/p_rcec_reproj.tif")
r.dust<-terra::rast("Projected_rasters/Projected_rasters_NA/p_dust_reproj.tif")


predictors <- c('r.w0010', 'r.srsrq3', 'r.dust', 
                'r.ocs', 'r.sw', 'r.w0010_30', 'r.cec')


#stack the predictors

stack_all <- c(r.w0010, r.srsrq3, r.dust, 
               r.ocs, r.sw, r.w0010_30, r.cec)


crs(stack_all) <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

names(stack_all) <- c('r.w0010', 'r.srsrq3', 'r.dust', 
                      'r.ocs', 'r.sw', 'r.w0010_30', 'r.cec')


###to define my extent
lat_long_bbox <- st_sfc(
  st_polygon(list(cbind(
    c(-130, -40, -40, -130, -130),
    c(25, 25, 55, 55, 25)
  ))),
  crs = "EPSG:4326"
)

eck4_crs <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
reprojected_bbox <- st_transform(lat_long_bbox, crs = eck4_crs)
new_eck4_extent <- st_bbox(reprojected_bbox)
print(new_eck4_extent)

ec <- extent(c(-11757406, -4604003 , 3232157, 6553182))
stack_all_cropped <- crop(stack_all, ec, snap='near')
NA.grid1 <- stack_all_cropped[[1]] / stack_all_cropped[[1]]

rf1 <- terra::predict(stack_all_cropped, RF_all, ext=NA.grid1, na.rm=TRUE, overwrite=TRUE)
writeRaster(rf1, filename="Output/Sr_mean.tif", overwrite=TRUE)
plot(rf1)

#create a mask to avoid extrapolation beyond our data range
training_rf2 <- subset(training_rf, select = -c(X87Sr_86Sr))
# Get the range of each predictor in the training data
min_values <- apply(training_rf2, 2, min)
max_values <- apply(training_rf2, 2, max)
# Create a raster with the same extent and resolution as your original stack
template_raster <- rast(stack_all_cropped[[1]])  # Creating a raster with the same extent and resolution
template_raster[] <- 1  # Set all values to 1

# Create a mask by checking if values in stack_all fall within the training range
mask <- app(stack_all_cropped, function(x) {
  all(x >= min_values & x <= max_values)
})

# Convert the mask to a binary raster
mask_raster <- ifel(mask, 1, NA)
mask_raster<-mask
writeRaster(mask_raster,"Output/mask.tif", overwrite=TRUE)

#Step5----
#quantile 

#this step took 6-7 hours on my computer.

### spatial uncertainty using Quantile regression forest using quantile 
#random forest 0.87 and 0.13 to calculate the spatial standard deviation. 

training_qrf<-read.table("sr_agg_all.txt",sep="\t",h=T)

#my model for reference
RF_all <- train(X87Sr_86Sr ~ r.w0010+r.srsrq3+r.dust+r.ocs+r.sw+r.w0010_30+r.cec,
                data = training_rf, ntree = 3000, method = "rf", importance=TRUE, 
                tuneGrid=tunegrid,trControl= fitControl)

qrf_all <- ranger(X87Sr_86Sr ~ r.w0010+r.srsrq3+r.dust+r.ocs+r.sw+r.w0010_30+r.cec,
                  data = training_qrf, quantreg=TRUE, num.trees=3000,mtry=2)   

# dataset need to be divided in several part to avoid memory problems
raster::rasterOptions(maxmemory = 12e+09)
stack_all2<-raster::stack(stack_all_cropped)
NA_pxl<-as(stack_all2,"SpatialPixelsDataFrame")
newdat0<-as.data.frame(NA_pxl)
cc<-which(complete.cases(NA_pxl@data))
lcc<-length(cc)
cc_div<-c(seq(0,lcc,100000),lcc) 

tab_pred<-as.data.frame(matrix(data=NA,ncol=3,nrow=0))

for (i in 1:(length(cc_div)-1))
{
  print(paste(i,(length(cc_div)-1),Sys.time(),sep=" / "))
  cc_sub<-cc[(cc_div[i]+1):cc_div[i+1]]
  newdat1<-newdat0[cc_sub,]
  s.rfd_low <- predict(qrf_all, data=newdat1, type="quantiles", quantiles=0.159, fun = function(model, ...) predict(model, ...)$predictions)
  gc()
  s.rfd_high <- predict(qrf_all, data=newdat1, type="quantiles", quantiles=0.841, fun = function(model, ...) predict(model, ...)$predictions)
  gc()
  pred<-cbind(cc_sub,s.rfd_low$predictions,s.rfd_high$predictions)
  tab_pred<-rbind(tab_pred,pred)
  rm(s.rfd_low)
  rm(s.rfd_high)
  rm(pred)
}

colnames(tab_pred)<-c("pxlnum","pred_low","pred_high")
tab_pred$error<-(tab_pred[,3]-tab_pred[,2])/2
NA_pxl@data$se.s<-NA
NA_pxl@data$se.s[cc]<-tab_pred$error
rast.se<-raster(NA_pxl[,"se.s"])
writeRaster(rast.se, filename="Output/Sr_sd.tif", format="GTiff", overwrite=TRUE)

#Step6----

#stack both layers for the assignR 

#load the layer (original data)
Sr.mean<-terra::rast("Sr_mean_filled.tif")
Sr.sd<-terra::rast("Sr_sd_filled.tif")  

plot(Sr.mean)
plot(Sr.sd)

stack_sr <- c(Sr.mean, Sr.sd)
crs(stack_sr) <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
names(stack_sr) <- c('Sr.mean', 'Sr.sd')
writeRaster(stack_sr, filename="Sr_2.tif", overwrite=TRUE)


#plot
setwd("C:/Users/abraa/OneDrive/Documents/Empoasca/Isotope/Stront/Output")

#plot Sr mean
Sr<-terra::rast("Sr_mean.tif")
plot(Sr)

sr_m1<-terra::project(Sr,"EPSG:5937",method="near")

plot(sr_m1)

plot(is.na(sr_m1), main = "NA cells (gaps + lakes)")
sr_filled <- terra::focal(sr_m1, w = 5, fun = "mean", na.policy = "only", na.rm = TRUE)
plot(sr_filled)
terra::writeRaster(sr_filled, "Sr_mean_filled.tif", overwrite = TRUE)

#plot
Sr<-terra::rast("Sr_mean_filled.tif")

# get country boundaries (medium scale is usually a good balance of detail vs. file size)
canada_sf <- ne_countries(country = "Canada", scale = "medium", returnclass = "sf")
usa_sf    <- ne_countries(country = "United States of America", scale = "medium", returnclass = "sf")
mex_sf    <- ne_countries(country = "Mexico", scale = "medium", returnclass = "sf")

# reproject to match sr_final's CRS
canada_sf2 <- st_transform(canada_sf, crs = "EPSG:5937")
usa_sf2    <- st_transform(usa_sf,    crs = "EPSG:5937")
mex_sf2    <- st_transform(mex_sf,    crs = "EPSG:5937")

# turn the SpatRaster into a data frame with x, y, value columns
sr_df <- as.data.frame(Sr, xy = TRUE, na.rm = TRUE)
names(sr_df)[3] <- "Sr_ratio"   # rename the value column if needed; check actual name first with names(sr_df)

breaks <- seq(floor(min(sr_df$Sr_ratio)*1000)/1000, 
              max(sr_df$Sr_ratio) + 0.003, 
              by = 0.003)

sr_df$Sr_bin <- cut(sr_df$Sr_ratio, breaks = breaks, include.lowest = TRUE)

# sanity check - no NAs
sum(is.na(sr_df$Sr_bin))

mapSr_final <- ggplot() +
  geom_raster(data = sr_df, aes(x = x, y = y, fill = Sr_bin)) +
  geom_sf(data = canada_sf2, fill = NA, color = "gray30", size = 0.3) +
  geom_sf(data = usa_sf2,    fill = NA, color = "gray30", size = 0.3) +
  geom_sf(data = mex_sf2,    fill = NA, color = "gray30", size = 0.3) +
  coord_sf(xlim = range(sr_df$x), ylim = range(sr_df$y), expand = FALSE) +
  scale_fill_viridis_d(
    name = "Sr Ratio",
    option = "H",
    drop = TRUE,
    labels = function(x) {
      gsub("\\(|\\]", "", x) %>%
        gsub(",", "-", .)
    }
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(16) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "right",
        legend.text = element_text(size = 12),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

mapSr_final

ggsave("mapSr_final.tiff", width = 200, height = 200, units = "mm",dpi=600, compression = "lzw")

#plot Sr quantiles
Sr<-terra::rast("Sr_sd.tif")

sr_sd<-terra::project(Sr,"EPSG:5937",method="near")

plot(sr_sd)

plot(is.na(sr_sd), main = "NA cells (gaps + lakes)")
sr_filled <- terra::focal(sr_sd, w = 5, fun = "mean", na.policy = "only", na.rm = TRUE)
plot(sr_sd)
terra::writeRaster(sr_filled, "Sr_sd_filled.tif", overwrite = TRUE)


#plot
Sr<-terra::rast("Sr_sd_filled.tif")
plot(Sr)

# get country boundaries (medium scale is usually a good balance of detail vs. file size)
canada_sf <- ne_countries(country = "Canada", scale = "medium", returnclass = "sf")
usa_sf    <- ne_countries(country = "United States of America", scale = "medium", returnclass = "sf")
mex_sf    <- ne_countries(country = "Mexico", scale = "medium", returnclass = "sf")

# reproject to match sr_final's CRS
canada_sf2 <- st_transform(canada_sf, crs = "EPSG:5937")
usa_sf2    <- st_transform(usa_sf,    crs = "EPSG:5937")
mex_sf2    <- st_transform(mex_sf,    crs = "EPSG:5937")

# turn the SpatRaster into a data frame with x, y, value columns
sr_df <- as.data.frame(Sr, xy = TRUE, na.rm = TRUE)
names(sr_df)[3] <- "Sr_sd"   # rename the value column if needed; check actual name first with names(sr_df)


mapSr_sd <- ggplot() +
  geom_raster(data = sr_df, aes(x = x, y = y, fill = Sr_sd)) +
  geom_sf(data = canada_sf2, fill = NA, color = "gray30", size = 0.3) +
  geom_sf(data = usa_sf2,    fill = NA, color = "gray30", size = 0.3) +
  geom_sf(data = mex_sf2,    fill = NA, color = "gray30", size = 0.3) +
  coord_sf(xlim = range(sr_df$x), ylim = range(sr_df$y), expand = FALSE) +
  scale_fill_viridis_c(
    name = "Sr (sd)",
    option = "H",
    breaks = c(min(sr_df$Sr_sd), mean(range(sr_df$Sr_sd)), max(sr_df$Sr_sd)),
    labels = scales::label_number(accuracy = 0.000001)
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(16) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "right",
        legend.text = element_text(size = 12),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

mapSr_sd

ggsave("mapSr_sd.tiff", width = 200, height = 200, units = "mm",dpi=600, compression = "lzw")



