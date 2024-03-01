# Clean the environment
rm(list = ls())
gc()


# Importing libraries
library(terra)
library(dplyr)
library(tidyterra)
library(lme4)
library(sf)
library(ggplot2)
library(MuMIn)
library(stringr)
library(purrr)
install.packages("MuMIn")

# Set working directory
setwd('I:/biocon/Emmanuel_Oceguera/projects/ChamoisHabitatConnectivityEurope_output')
getwd()

# Load in vector of European nations and ReWilding grid of Europe
studyarea = st_read("I:\\biocon\\Jeremy_Dertien\\NaturaConnect_Analysis\\Polygons\\NC_Countries_reproj.shp")
spgrid = vect("I:\\biocon\\Emmanuel_Oceguera\\projects\\2024_01_Mammals_species_distribution_DarwingCore\\data\\occurrances_data\\sp_occu_20190829.shp")

# Assuming spgrid and other required libraries are already loaded

# Extract genus part from column names
genus_names <- names(spgrid)[4:20] %>% 
  str_extract(pattern = "(?<=_)[a-z]{3}") %>%
  na.omit() %>%
  unique()


# Isolate species by genus and sum into one column per genus
spgrid_by_genus <- list()

for (genus in genus_names) {
  # Identify columns belonging to this genus
  cols <- grep(paste0("_", genus), names(spgrid), value = TRUE)
  
  # Select and sum columns for each genus, replace 9 with 0
  spgrid_by_genus[[genus]] <- spgrid %>%
    select(cellcode, eoforigin, noforigin, all_of(cols)) %>%
    rowwise() %>%
    mutate(genus_sum = sum(c_across(all_of(cols)), na.rm = TRUE),
           genus_sum = if_else(genus_sum == 9, 0, genus_sum)) %>%
    ungroup() %>%
    select(cellcode, eoforigin, noforigin, genus_sum)
  
  # Rename the summed column to reflect the genus
  names(spgrid_by_genus[[genus]])[4] <- paste0("m_", genus, "spp")
}


# Convert all the spatVect to sf object
sf_mammals_eu_lt <- lapply(spgrid_by_genus, sf::st_as_sf)
names(sf_mammals_eu_lt) <- paste0("m_", genus_names, "spp")


# Select the only presence squares with only 1 values
sf_mammals_eu2_lt <- lapply(sf_mammals_eu_lt, function(sf_obj){
  #Fiterl rows where any numeric column has a values of 1
  sf_obj %>% 
    filter(if_any(where(is.numeric), ~ .x == 1))
})

# check
sf_mammals_eu2_lt
# Import the rupssp to check it
m_rupspp_eu2 <- sf_mammals_eu2_lt[["m_rupspp"]]
st_write(m_rupspp_eu2, 'eu/m_rupspp_eu2.shp')



# Function to generate points only at the center of each square
st_centroid_within_poly <- function (poly) {
  
  #Check if centroid is in polygon
  ctrd <- st_centroid(poly, of_largest_polygon = TRUE)
  in_poly <- diag(st_within(ctrd, poly, sparse = F))
  
  #Replace geometries that are not within polygon with st_point_on_surface()
  st_geometry(ctrd[!in_poly,]) <- st_geometry(st_point_on_surface(poly[!in_poly,]))
  
  return(ctrd)
  
}


##Create points for presence polygons
mammals_eu_pts2_lt <- lapply(sf_mammals_eu2_lt, st_centroid_within_poly)

# Import the rupssp to check it
m_rupspp_pts2 <- mammals_eu_pts2_lt[["m_rupspp"]]
st_write(m_rupspp_pts2, 'eu/m_rupspp_pts2.shp')


##Buffer around the presence squares to produce habitat availability points
# Define the custom function to apply the operations
process_genus_sf <- function(sf_object) {
  
  # Buffer around the presence points to produce habitat availability points, buffer 40 km
  buff_40km <- st_buffer(sf_object, dist = 40000)
  
  # Dissolve the buffer
  buff_40km_union <- st_as_sf(data.frame(id = 1), geometry = st_union(buff_40km))
  
  # Optional: Fix invalid geometries if necessary
  buff_40km_union_valid <- st_make_valid(buff_40km_union)
  
  result <- buff_40km_union_valid
  
  return(result)
}

# Assuming genus_layers_no_na is your list of sf objects
mammals_eu_pts2_buff_40km_lt <- lapply(mammals_eu_pts2_lt, process_genus_sf)

# Import the rupssp to check it
m_rupspp_buff_sf <- mammals_eu_pts2_buff_40km_lt[["m_rupspp"]]
st_write(m_rupspp_buff_sf, 'eu/m_rupspp_pts2_buff_40km.shp')


# Define the st_erase function
st_erase <- function(x, y) {
  st_difference(x, st_union(st_combine(y)))
}

# Erase the presence squares to obtain a fully dissolved single polygon on the outside
# Initialize an empty list to store the results
mammals_eu_buff_40km_dif_lt <- list()

# Loop through the names of the first list
for (name in names(mammals_eu_pts2_buff_40km_lt)) {
  # Check if the name exists in the second list
  if (name %in% names(sf_mammals_eu2_lt)) {
    # Perform the erase operation using matched names
    mammals_eu_buff_40km_dif_lt[[name]] <- st_erase(mammals_eu_pts2_buff_40km_lt[[name]], sf_mammals_eu2_lt[[name]])
  } else {
    # If no matching name, just a placeholder to indicate, or you can choose to skip
    mammals_eu_buff_40km_dif_lt[[name]] <- NULL
    print(paste("No matching sf object for:", name))
  }
}

# Import the rupssp to check it
m_rupspp_buff_40km_dif <- mammals_eu_buff_40km_dif_lt[["m_rupspp"]]
st_write(m_rupspp_buff_40km_dif, 'eu/m_rupspp_buff_40km_dif.shp')


###Clip buffer by nations and full grid by the available buffer area
# Initialize an empty list to store the intersection results
mammals_eu_buff_country_lt  <- list()

# Loop through each sf object in the 'm_rupspp_buff_40km_dif' list
for (name in names(mammals_eu_buff_40km_dif_lt)) {
  # Perform intersection with 'studyarea' for each sf object
  mammals_eu_buff_country_lt[[name]] <- st_intersection(mammals_eu_buff_40km_dif_lt[[name]], studyarea)
}

# Import the rupssp to check it
m_rupspp_buff_country_clip <- mammals_eu_buff_country_lt[["m_rupspp"]]
st_write(m_rupspp_buff_country_clip, 'eu/m_rupspp_buff_country_clip.shp')


# Initialize an empty list to store the clipped results
mammals_eu_buff_country_clip_lt <- list() 

# Correctly use 'name' as the loop variable
for(name in names(mammals_eu_buff_country_lt)){
  
  # Check if 'name' exists in 'sf_mammals_eu_lt' to avoid errors
  if(name %in% names(sf_mammals_eu_lt)) {
    mammals_eu_buff_country_clip_lt[[name]] <- st_intersection(sf_mammals_eu_lt[[name]], 
                                                               mammals_eu_buff_country_lt[[name]])
  } else {
    warning(paste(name, 'not found in sf_mammals_eu_lt'))
  }
}

names(mammals_eu_buff_country_clip_lt)

# # Import the rupssp to check it
# m_rupspp_buff_country_clip_clip <- mammals_eu_buff_country_clip_lt[["m_rupspp"]]
# st_write(m_rupspp_buff_country_clip_clip, 'eu/m_rupspp_buff_country_clip_clip.shp')


# Create point form the 'availables' squares
##Create points for presence polygons
mammals_eu_avail_pts <- lapply(mammals_eu_buff_country_clip_lt, function(sf_objt){
  #Selected the first 4 columns along with the geometry columns
  selected_sf_obj <- sf_objt %>% 
    select(1:4, geometry)
  
  st_centroid_within_poly(selected_sf_obj)
  
  })


# # Import the rupssp to check it
m_rupspp_avail_pts <- mammals_eu_avail_pts[["m_rupspp"]]
st_write(m_rupspp_avail_pts, 'eu/m_rupspp_avail_pts.shp')


# Filter only available points
filtered_mammals_eu_avail_pts <- lapply(mammals_eu_avail_pts, function(sf_obj){
  #Fiterl rows where any numeric column has a values of 1
  sf_obj %>% 
    filter(if_any(where(is.numeric), ~ .x == 0))
})


filtered_mammals_eu_avail_pts[[1]]
mammals_eu_avail_pts[[1]]

## Final combination of presence and available point 
eu_mammals <- list() 
  
for (spp_names in names(mammals_eu_pts2_lt)) {
  # Retrieve the matching sf object from both list
  sf_objt1 <- mammals_eu_pts2_lt[[spp_names]]
  sf_objt2 <- filtered_mammals_eu_avail_pts[[spp_names]]
  
  # Combine the sf object usign rbind, ensuring thee have the same structure
  combined_sf <- rbind(sf_objt1, sf_objt2)
  
  # Store the combined sf object in the new list
  eu_mammals[[spp_names]] <- combined_sf 
} 


# # Import the rupssp to check it
m_rupspp_pres_avail_comb <- eu_mammals[["m_rupspp"]]
st_write(m_rupspp_pres_avail_comb, 'eu/m_rupspp_pres_avail_comb.shp')


##Loading in covariates
ungcovs = rast("I:\\biocon\\Jeremy_Dertien\\NaturaConnect_Analysis\\Covariates\\LULC_ZMeans\\Kilo1\\UngCovsRastStack2.grd")        
biogeo = rast("I:\\biocon\\Jeremy_Dertien\\NaturaConnect_Analysis\\Covariates\\LULC_ZMeans\\Kilo1\\Biogeographic2016.tif")
# plot(ungcovs)
# plot(biogeo)

# Resmaple bioregion so it has the same ext
biogeo = resample(biogeo, ungcovs, method = "bilinear")

# Add bio region to the covariants and rename it
mammalscovs = c(ungcovs, biogeo) %>% 
  rename(biogeo = short_name)

summary(mammalscovs)

# Extracting valus for all the covariants at each point
eu_mammals_cov_values <- lapply(eu_mammals, function(sf_obj){
  # Convert sf object to SpatVector for compatibility with terra::extract()
  spat_vector <- vect(sf_obj)
  #Extract covariate values for points in the SpatVector
  df <- terra::extract(mammalscovs, spat_vector, bind = TRUE) |>
    as.data.frame()
   
  # Optionally convert to a data frame if extract() doesn't automatically and if there are some Na we replace them
  df <- df %>% 
    mutate(across(where(is.numeric), ~replace_na(., 0)))
  
  return(df)
})

# # Import the rupssp to check the values
m_rupspp_cov_values <- eu_mammals_cov_values[["m_rupspp"]]
write.table(m_rupspp_cov_values, 'eu/m_rupspp_cov_values.txt', sep = ',')

########TEST###########
# Placeholder list to store models
models_list <- list()

# Loop through each species data frame in your list
models_list <- lapply(names(eu_mammals_cov_values), function(species_name) {
  # Access the data frame for the current species
  df <- eu_mammals_cov_values[[species_name]]
  
  # Construct the model formula dynamically if needed
  # mofify the morfula?
  formula <- as.formula(paste(species_name, "~ RdDen10k_Zmean + Forest_ZMean + TRI10km_reSamp + (1|biogeo)", sep = " "))
  
  # Fit the model using glmer
  model <- glmer(formula, family = "binomial", data = df)
  
  # Return the fitted model
  return(model)
})

# Assign names to the models in the list to match the species names
names(models_list) <- names(eu_mammals_cov_values)
summary(models_list[["m_rupspp"]])








# # Ensure the list is not empty before proceeding
# if (length(spgrid_by_genus) > 0) {
#   # Start with the first SpatVector in the list
#   combined_spatvector <- spgrid_by_genus[[1]]
#   
#   # If there are more SpatVectors in the list, append them one by one
#   if (length(spgrid_by_genus) > 1) {
#     for (i in 2:length(spgrid_by_genus)) {
#       combined_spatvector <- rbind(combined_spatvector, spgrid_by_genus[[i]])
#     }
#   }
#   
#   # Print or check the structure of the combined SpatVector
#   print(combined_spatvector)
# } else {
#   message("The list spgrid_by_genus is empty.")
# }
# 
# 
# 
# # Check the combined data frame
# head(combined_spatvector)

# Convert to sf
# sf_mammals_eu <- as_sf(combined_spatvector)

# # Convert all the Nan to 0
# sf_mammals_eu <- sf_mammals_eu %>% 
#   mutate(across(names(sf_mammals_eu)[4:14], ~replace(., is.nan(.), 0)))

# sf_mammals_eu2_0 <- sf_mammals_eu %>%
#   mutate_if(is.numeric, ~replace_na(., 0))

# # Assuming sf_mammals_eu is your sf object
# sf_mammals_eu2 <- sf_mammals_eu %>%
#   # Filter out rows that have only NA values in columns 4:20
#   filter(rowSums(is.na(select(., 4:14))) < (14-4+1))
# 

# # gettin names
# names_genus <- grep("^m_", names(sf_mammals_eu2), value = TRUE) # Identify genus columns
# # List to store each genus layer
# sf_mammals_eu2_lt <- list()
# 
# for (genus in names_genus) {
#   # Select the required columns for each genus
#   genus_layer <- sf_mammals_eu2 %>%
#     select(cellcode, eoforigin, noforigin, geometry, all_of(genus)) %>% 
#     
#     filter(!is.na(!!sym(genus)))
#   # Add the new layer to the list
#   sf_mammals_eu2_lt[[genus]] <- genus_layer
# }
# 

# 
# # List to store each genus layer with 0
# sf_mammals_eu2_0_lt <- list()
# 
# for (genus in names_genus) {
#   # Select the required columns for each genus
#   genus_layer <- sf_mammals_eu2_0 %>%
#     select(cellcode, eoforigin, noforigin, geometry, all_of(genus))
#     
#   # Add the new layer to the list
#   sf_mammals_eu2_0_lt[[genus]] <- genus_layer
# }
# 
# # Select only the presence squares cham with only 1 values
# sf_mammals_eu2 <- sf_mammals_eu %>% 
#   filter(names(sf_mammals_eu)[4:14] == 1)


# # Function to generate points only at the center of each square
# st_centroid_within_poly <- function (poly) {
#   
#   #Check if centroid is in polygon
#   ctrd <- st_centroid(poly, of_largest_polygon = TRUE)
#   in_poly <- diag(st_within(ctrd, poly, sparse = F))
#   
#   #Replace geometries that are not within polygon with st_point_on_surface()
#   st_geometry(ctrd[!in_poly,]) <- st_geometry(st_point_on_surface(poly[!in_poly,]))
#   
#   return(ctrd)
#   
# }

# ##Create points for presence polygons
# mammals_eu_pts2 <-  st_centroid_within_poly(sf_mammals_eu2)
# 
# # print(mammals_eu_pts2)
# 
# mammals_eu_pts2_lt = lapply(sf_mammals_eu2_lt, st_centroid_within_poly)
# # Import the rupssp to check it
# m_rupspp_pts2 <- mammals_eu_pts2_lt[["m_rupspp"]]
# st_write(m_rupspp_pts2, 'eu/m_rupspp_pts2.shp')
# 
# 


# Now, genus_layers_centroids is a list where each element is an sf object
# with points only at the center of each polygon.

# ##Buffer around the presence squares to produce habitat availability points
# # Define the custom function to apply the operations
# process_genus_sf <- function(sf_object) {
#   
#   # Buffer around the presence points to produce habitat availability points, buffer 40 km
#   buff_40km <- st_buffer(sf_object, dist = 40000)
#   
#   # Dissolve the buffer
#   buff_40km_union <- st_as_sf(data.frame(id = 1), geometry = st_union(buff_40km))
#   
#   # Optional: Fix invalid geometries if necessary
#   buff_40km_union_valid <- st_make_valid(buff_40km_union)
#   
#   result <- buff_40km_union_valid
#   
#   return(result)
# }
# 
# # Assuming genus_layers_no_na is your list of sf objects
# mammals_eu_buff_40km_lt <- lapply(mammals_eu_pts2_lt, process_genus_sf)
# 
# # Import the rupssp to check it
# m_rupspp_feature_sf <- mammals_eu_buff_40km_lt[["m_rupspp"]]
# st_write(m_rupspp_feature_sf, 'eu/m_rupspp_buff_40km.shp')
# 
# 
# # Define the st_erase function
# st_erase <- function(x, y) {
#   st_difference(x, st_union(st_combine(y)))
# }
# # 
# # # Proceed with the combination if all elements are sf objects
# # combined_sf <- dplyr::bind_rows(sf_mammals_eu2_lt)
# 
# # Initialize an empty list to store the results
# mammals_eu_buff_40km_dif_lt <- list()
# 
# # Loop through the names of the first list
# for (name in names(mammals_eu_buff_40km_lt)) {
#   # Check if the name exists in the second list
#   if (name %in% names(sf_mammals_eu2_lt)) {
#     # Perform the erase operation using matched names
#     mammals_eu_buff_40km_dif_lt[[name]] <- erase_sf(mammals_eu_buff_40km_lt[[name]], sf_mammals_eu2_lt[[name]])
#   } else {
#     # If no matching name, just a placeholder to indicate, or you can choose to skip
#     mammals_eu_buff_40km_dif_lt[[name]] <- NULL
#     print(paste("No matching sf object for:", name))
#   }
# }

# # Corrected erase_sf function to apply st_erase
# erase_sf <- function(sf_object, sf_ref_object) {
#   # Apply st_erase using sf_object and a reference sf object (sf_ref_object)
#   buff_dif <- st_difference(sf_object, st_union(sf_ref_object))
#   # Attempt to make the result valid
#   buff_dif_valid <- st_make_valid(buff_dif)
#   # Return the valid difference
#   return(buff_dif_valid)
# }
# 
# # Use lapply to apply erase_sf to each sf object in mammals_eu_buff_40km_lt
# # Ensure sf_mammals_eu2_lt is correctly defined as a single sf object for this operation
# mammals_eu_buff_40km_dif <- lapply(mammals_eu_buff_40km_lt, function(x) erase_sf(x, combined_sf))

# # Import the rupssp to check it
# m_rupspp_buff_40km_dif <- mammals_eu_buff_40km_dif_lt[["m_rupspp"]]
# st_write(m_rupspp_buff_40km_dif, 'eu/m_rupspp_buff_40km_dif.shp')

preprocessing


###Clip buffer by nations and full grid by the available buffer area
# Initialize an empty list to store the intersection results
# buf_country  <- list()
# 
# # Loop through each sf object in the 'results' list
# for (name in names(results)) {
#   # Perform intersection with 'studyarea' for each sf object in 'results'
#   buf_country[[name]] <- st_intersection(results[[name]], studyarea)
# }

# ##Buffer around the presence squares to produce habitat availability points
# # buffer 40 km 
# mammals_eu_buff_40km <- st_buffer(mammals_eu_pts2, dist = 40000)
# # dissolve the buffer
# mammals_eu_buff_40km_union <- st_as_sf(data.frame(id=1), geometry = st_union(mammals_eu_buff_40km))
# st_write(mammals_eu_buff_40km_union, 'eu/mammals_eu_buff_40km_union.shp')
# # # Fix invalid geometries
# # mammals_eu_buff_40km_union <- st_make_valid(mammals_eu_buff_40km_union)
# # sf_mammals_eu2 <- st_make_valid(sf_mammals_eu2)
# # Erase the presence squares to obtain a fully dissolved single polygon on the outside
# st_erase = function(x, y) st_difference(x, st_union(st_combine(y)))
# mammals_eu_buff_40km_diss <- st_erase(mammals_eu_buff_40km_union, st_make_valid(sf_mammals_eu2))



