# Function to extract covariables at the sampling sites or at the closest of the sampling site
# if the site fall within an empty cell
extract_nearest_non_na <- function(coordinates, raster) {
  # Convert the coordinates data.frame to an sf points object
  points <- st_as_sf(coordinates, coords = c("X", "Y"), crs = st_crs(raster))
  
  # Extract raster values under the points
  values <- terra::extract(raster, points,method="simple")
  values <- values[,2]
  
  # Identify NA values
  na_indices <- is.na(values)
  
  # Find the indices of NA values
  na_indices <- which(na_indices)
  
  # Loop over each NA value
  for (i in na_indices) {
    # Get the current point
    point <- points[i, ]
    st_crs(point) <- st_crs(raster)
    # Create a buffer around the point with a 50km radius
    buffer <- st_buffer(point, dist = 50000)
    
    # Convert the buffer to the same CRS as the raster
    st_crs(buffer) <- st_crs(raster)
    
    # Crop the raster to the buffered extent
    cropped_raster <- crop(raster, buffer)
    
    # Extract the non-NA values within the buffer
    non_na_values <- cropped_raster[]
    non_na_values <- non_na_values[!is.na(non_na_values)]
    
    # Find the nearest non-NA value
    nearest_value <- non_na_values[which.min(st_distance(point, buffer))]
    
    # Assign the nearest non-NA value to replace the NA value
    values[i] <- nearest_value
  }
  
  # Return the extracted values
  return(data.frame(values))
}
