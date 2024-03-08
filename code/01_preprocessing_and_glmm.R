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
setwd('I:\\biocon\\Emmanuel_Oceguera\\projects\\2024_03_EruoMammalsConnectivity\\outcomes')
getwd()

# Load in vector of European nations and ReWilding grid of Europe
studyarea = st_read("I:\\biocon\\Jeremy_Dertien\\NaturaConnect_Analysis\\Polygons\\NC_Countries_reproj.shp")
spgrid = vect("I:\\biocon\\Emmanuel_Oceguera\\projects\\2024_01_Mammals_species_distribution_DarwingCore\\data\\occurrances_data\\sp_occu_20190829.shp")


# Extract genus part from column names
genus_names <- names(spgrid)[4:20] %>% 
  str_extract(pattern = "(?<=_)[a-z]{3}") %>%
  na.omit() %>%
  unique()


# Isolate species by genus and sum into one column per genus

eu_mammals_by_genus <- list()

for (genus in genus_names) {
  
  # Identify columns belonging to this genus
  cols <- grep(paste0("_", genus), names(spgrid), value = TRUE)
  
  # Select and sum columns for each genus, replace 9 with 0
  eu_mammals_by_genus[[genus]] <- spgrid %>%
    select(cellcode, eoforigin, noforigin, all_of(cols)) %>%
    rowwise() %>%
    mutate(genus_sum = sum(c_across(all_of(cols)), na.rm = TRUE),
           genus_sum = if_else(genus_sum == 9, 0, genus_sum)) %>%
    ungroup() %>%
    select(cellcode, eoforigin, noforigin, genus_sum)
  
  # Rename the summed column to reflect the genus
  names(eu_mammals_by_genus[[genus]])[4] <- paste0("m_", genus, "spp")
}


# Convert all the spatVect to sf object
sf_mammals_eu_lt <- lapply(eu_mammals_by_genus, sf::st_as_sf)
names(sf_mammals_eu_lt) <- paste0("m_", genus_names, "spp")


# Select the only presence squares with only 1 values
sf_mammals_eu2_lt <- lapply(sf_mammals_eu_lt, function(sf_obj){
  #Fiterl rows where any numeric column has a values of 1
  sf_obj %>% 
    filter(if_any(where(is.numeric), ~ .x == 1))
})

# # check
# sf_mammals_eu2_lt
# # Import the rupssp to check it
# m_rupspp_eu2 <- sf_mammals_eu2_lt[["m_rupspp"]]
# st_write(m_rupspp_eu2, 'm_rupspp_eu2.shp')



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
eu_mammals_pts2_lt <- lapply(sf_mammals_eu2_lt, st_centroid_within_poly)


# # Import the rupssp to check it
# m_rupspp_pts2 <- eu_mammals_pts2_lt[["m_rupspp"]]
# st_write(m_rupspp_pts2, 'm_rupspp_pts2.shp')


##Buffer around the presence squares to produce habitat availability points
# function to accept the presence column name and buffer distances as arguments
st_buffer_presence_squeres <- function(sf_object, presence_column, buffer_distances) {
  # Assuming the presence column effectively indicates the species
  # And that each sf object is for a single species, so taking the first value is sufficient
  species_name <- names(buffer_distances)[which(names(buffer_distances) == presence_column)]
  
  # Extract the buffer distance for this species/presence column
  dist <- buffer_distances[presence_column]
  
  # Proceed with buffering, dissolving, and fixing geometries
  buff <- st_buffer(sf_object, dist = dist)
  buff_union <- st_as_sf(data.frame(id = 1), geometry = st_union(buff))
  buff_union_valid <- st_make_valid(buff_union)
  
  return(buff_union_valid)
}

# Define buffer distances for each species based on the presence column name
buffer_distances <- c("m_canspp" = 100000, "m_ursspp" = 50000, "m_lynspp" = 60000, 
                      "m_gulspp" = 40000, "m_alcspp" = 40000, "m_capspp" = 30000,
                      "m_cerspp" = 40000, "m_rupspp" = 40000, "m_susspp" = 70000,
                      "m_bisspp" = 50000, "m_ranspp" = 60000)


# Dynamically determine the presence column for each sf object
# Applying the function with a manually specified column 
eu_mammals_pts2_buff_varied_lt <- lapply(names(eu_mammals_pts2_lt), function(x) {
  sf_object <- eu_mammals_pts2_lt[[x]]
  presence_column <- x # Assuming the key in your list somehow indicates the presence column or species
  st_buffer_presence_squeres(sf_object, presence_column, buffer_distances)
})


# Assigning object names
names(eu_mammals_pts2_buff_varied_lt) <- paste0("m_", genus_names, "spp")

# Import the rupssp and canspp to check the buffers
# m_rupspp_buff_40km <- eu_mammals_pts2_buff_varied_lt[["m_rupspp"]]
# m_canspp_buff_40km <- eu_mammals_pts2_buff_varied_lt[["m_canspp"]]
# 
# st_write(m_rupspp_buff_40km, 'm_rupspp_pts2_buff_40km.shp')
# st_write(m_canspp_buff_40km, 'm_canspp_pts2_buff_40km.shp')
# 


# Buffer function with unique value -------------------------------

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
# mammals_eu_pts2_buff_40km_lt <- lapply(mammals_eu_pts2_lt, process_genus_sf)
# 
# # Import the rupssp to check it
# m_rupspp_buff_sf <- mammals_eu_pts2_buff_40km_lt[["m_rupspp"]]
# st_write(m_rupspp_buff_sf, 'eu/m_rupspp_pts2_buff_40km.shp')


# Continues ---------------------------------------------------------------


# Define the st_erase function
st_erase <- function(x, y) {
  st_difference(x, st_union(st_combine(y)))
}

# Erase the presence squares to obtain a fully dissolved single polygon on the outside
eu_mammals_buff_diff_lt <- list()

# Loop through the names of the first list
for (name in names(eu_mammals_pts2_buff_varied_lt)) {
  # Check if the name exists in the second list
  if (name %in% names(sf_mammals_eu2_lt)) {
    # Perform the erase operation using matched names
    eu_mammals_buff_diff_lt[[name]] <- st_erase(eu_mammals_pts2_buff_varied_lt[[name]], sf_mammals_eu2_lt[[name]])
  } else {
    # If no matching name, just a placeholder to indicate, or you can choose to skip
    eu_mammals_buff_diff_lt[[name]] <- NULL
    print(paste("No matching sf object for:", name))
  }
}

# Import the rupssp and canspp to check the buffers
# m_rupspp_buff_dif <- eu_mammals_buff_diff_lt[["m_rupspp"]]
# m_canspp_buff_dif <- eu_mammals_buff_diff_lt[["m_canspp"]]
# 
# st_write(m_rupspp_buff_dif, 'm_rupspp_buff_dif.shp')
# st_write(m_canspp_buff_dif, 'm_canspp_buff_dif.shp')


###Clip buffer by nations and full grid by the available buffer area
eu_mammals_buff_diff_ctry_lt <- list()

# Loop through each sf object in the 'm_rupspp_buff_40km_dif' list
for (name in names(eu_mammals_buff_diff_lt)) {
  # Perform intersection with 'studyarea' for each sf object
  eu_mammals_buff_diff_ctry_lt[[name]] <- st_intersection(eu_mammals_buff_diff_lt[[name]], studyarea)
}

# Import the rupssp and canspp to check the buffers clipped with countries
# m_rupspp_buff_dif_ctry <- eu_mammals_buff_diff_ctry_lt[["m_rupspp"]]
# m_canspp_buff_dif_ctry <- eu_mammals_buff_diff_ctry_lt[["m_canspp"]]
# 
# st_write(m_rupspp_buff_dif_ctry, 'm_rupspp_buff_dif_ctry.shp')
# st_write(m_canspp_buff_dif_ctry, 'm_canspp_buff_dif_ctry.shp')


# Clip the buffer area with the inicial grid with 0 and 1
eu_mammals_buff_diff_ctry_clip_lt <- list() 

# Correctly use 'name' as the loop variable
for(name in names(eu_mammals_buff_diff_ctry_lt)){
  
  # Check if 'name' exists in 'sf_mammals_eu_lt' to avoid errors
  if(name %in% names(sf_mammals_eu_lt)) {
    eu_mammals_buff_diff_ctry_clip_lt[[name]] <- st_intersection(sf_mammals_eu_lt[[name]], 
                                                               eu_mammals_buff_diff_ctry_lt[[name]])
  } else {
    warning(paste(name, 'not found in sf_mammals_eu_lt'))
  }
}

###Create point form the 'availables' squares
#Create points for presence polygons
eu_mammals_avail_pts <- lapply(eu_mammals_buff_diff_ctry_clip_lt, function(sf_objt){
  #Selected the first 4 columns along with the geometry columns
  selected_sf_obj <- sf_objt %>% 
    select(1:4, geometry)
  
  st_centroid_within_poly(selected_sf_obj)
  
  })

# Filter only available points, so all the 0
eu_mammals_avail_pts_cero <- lapply(eu_mammals_avail_pts, function(sf_obj){
  #Filter rows where any numeric column has a values of 1
  sf_obj %>% 
    filter(if_any(where(is.numeric), ~ .x == 0))
})

## Final combination of presence and available point 
eu_mammals_pres_avail_pts <- list() 
  
for (name in names(eu_mammals_pts2_lt)) {
  # Retrieve the matching sf object from both list
  sf_objt1 <- eu_mammals_pts2_lt[[name]]
  sf_objt2 <- eu_mammals_avail_pts_cero[[name]]
  
  # Combine the sf object usign rbind, ensuring thee have the same structure
  combined_sf <- rbind(sf_objt1, sf_objt2)
  
  # Store the combined sf object in the new list
  eu_mammals_pres_avail_pts[[name]] <- combined_sf 
} 

#summary(eu_mammals_pres_avail_pts[["m_rupspp"]])

# Import the rupssp and canspp to check the buffers clipped with countries
# m_rupspp_pres_avail_pts <- eu_mammals_pres_avail_pts[["m_rupspp"]]
# m_canspp_pres_avail_pts <- eu_mammals_pres_avail_pts[["m_canspp"]]
# 
# st_write(m_rupspp_pres_avail_pts, 'm_rupspp_pres_avail_pts.shp')
# st_write(m_canspp_pres_avail_pts, 'm_canspp_pres_avail_pts.shp')
# 


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


# Extracting valus for all the covariants at each point
eu_mammals_covs_values <- lapply(eu_mammals_pres_avail_pts, function(sf_obj){
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


# Import the rupssp and canspp to check the buffers clipped with countries
# m_rupspp_covs_values <- eu_mammals_covs_values[["m_rupspp"]]
# m_canspp_covs_values <- eu_mammals_covs_values[["m_canspp"]]
# 
# write.table(m_rupspp_covs_values, 'm_rupspp_covs_values.txt', sep = ',')
# write.table(m_canspp_covs_values, 'm_canspp_covs_values.txt', sep = ',')
# 

# Reference formula GLMM 
# full_model <- lme4::glmer(m_rupspp ~ Dvlpd_ZMean + RdDen10k_Zmean + Forest_ZMean + TRI10km_reSamp + I(TRI10km_reSamp^2) + Maxtmp10k_resamp + Wland_ZMean + Shrub10k_resamp + (1|biogeo),
#                           data = chamoisdf1, family = binomial)

# Function to fit models and perform model selection for each species
fit_best_model_for_species <- function(df, species_name) {
  
  options(na.action = "na.fail")
  # Dynamically create the full model formula
  # Assuming species_name is used as the response variable; adjust if your structure is different
  formula <- as.formula(paste0(species_name, " ~ Dvlpd_ZMean + RdDen10k_Zmean + 
                               Forest_ZMean + TRI10km_reSamp + I(TRI10km_reSamp^2) + 
                               Maxtmp10k_resamp + Wland_ZMean + Shrub10k_resamp + (1|biogeo)"))
  
  # Fit the full model
  full_model <- glmer(formula, data = df, family = binomial, control = glmerControl(optimizer = "bobyqa"))
  
  # Perform model selection using dredge from the MuMIn package
  model_set <- MuMIn::dredge(full_model)
  
  # Select the best model based on the lowest AIC
  best_model <- MuMIn::get.models(model_set, 1)[[1]]
  
  # Optionally, perform model averaging
  evg_model <- model.avg(model_set, subset = delta < 4) # Adjust the delta according to your criteria
  
  # Return the best model and the averaged model
  return(list(best_model = best_model, evg_model = evg_model))
}


# Loop through each species data frame in your list
models_list <- lapply(names(eu_mammals_covs_values), function(species_name) {
  df <- eu_mammals_covs_values[[species_name]]
  
  # Fit the model and perform model selection for this species
  models <- fit_best_model_for_species(df, species_name)
  
  # Return the fitted models
  return(models)
})


# Parallelize Raster Prediction Generation
eu_mammals_pred_rast <- future_map(names(models_list), function(spp_name) {
  best_model <- models_list[[spp_name]]$best_model
  pred_rast <- generate_raster_prediction(best_model, mammalscovs, spp_name)
  return(pred_rast)
}, .options = furrr_options(globals = c("generate_raster_prediction", "models_list", "mammalscovs")))






