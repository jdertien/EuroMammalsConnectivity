# Clean the environment
rm(list = ls())
gc()

# Importing necessary libraries
library(terra)    # For spatial data manipulation
library(tidyterra)# For data manipulation spatVector and spatRaster
library(dplyr)    # For data manipulation
library(lme4)     # For fitting generalized linear mixed models (GLMM)
library(sf)       # For handling spatial data in simple features format
library(MuMIn)    # For model selection and averaging
library(stringr)  # For string manipulation

# Utility Functions
# -----------------
# Function to generate points only at the center of each square
st_centroid_within_poly <- function(poly) {
  # Check if centroid is in polygon
  ctrd <- st_centroid(poly, of_largest_polygon = TRUE)
  in_poly <- diag(st_within(ctrd, poly, sparse = FALSE))
  # Replace geometries that are not within polygon
  st_geometry(ctrd[!in_poly,]) <- st_geometry(st_point_on_surface(poly[!in_poly,]))
  
  return(ctrd)
}
# Function to apply buffers and dissolve with variable distances
st_buffer_presence_squeres <- function(sf_object, presence_column, buffer_distances) { # nolint: line_length_linter.
  # Assuming the presence column effectively indicates the species and that each sf object is for a single species, # nolint: line_length_linter.
  # so taking the first values is sufficient
  species_name <- names(buffer_distances)[which(names(buffer_distances) == presence_column)] # nolint: line_length_linter.
  # Extract the buffer distance got this species/presence column
  dist <- buffer_distances[presence_column]
  # Proceed with buffering, dissolving and fixing geometries
  buff <- st_buffer(sf_object, dist = dist)
  buff_union <- st_as_sf(data.frame(id = 1), geometry = st_union(buff))
  buff_union_valid <- st_make_valid(buff_union)
  
  return(buff_union_valid)
}
# Function to erase
st_erase <- function(x, y) {
  st_difference(x, st_union(st_combine(y)))
}
# Function for model fitting and selection
fit_best_model_for_species <- function(df, species_name) {
  options(na.action = "na.fail")
  # Dinamically create the full model formula
  formula <- as.formula(paste0(species_name, " ~ Dvlpd_ZMean + RdDen10k_Zmean + 
                               Forest_ZMean + TRI10km_reSamp + I(TRI10km_reSamp^2) +  # nolint: line_length_linter.
                               Maxtmp10k_resamp + Wland_ZMean + Shrub10k_resamp + (1|biogeo)")) # nolint: line_length_linter.
  # Fit the full model
  full_model <- glmer(formula, data = df, family = binomial, control = glmerControl(optimizer = "bobyqa")) # nolint: line_length_linter.
  # Model selection table of models with combinations
  model_set <- dredge(full_model)
  # Get the best model with the lowest AIC
  best_model <- get.models(model_set, 1)[[1]]
  
  return(list(best_model = best_model))
}
# Prediction Function
predict_for_species <- function(fitted_model, species_name, prediction_grid) {
  # Fitted model is a list containing the best model for each specie
  best_model <- fitted_model$best_model
  # Generate predictions
  species_predictions <- terra::predict(prediction_grid, best_model, type = "response", re.form = NA) # nolint: line_length_linter.
  # Set CRS
  crs(species_predictions) <- "epsg:3035"
  # Save the results as TIF
  file_path <- paste0(species_name, "_predictions.tif")
  writeRaster(species_predictions, filename = file_path, overwrite = T)
  return(file_path)
}

# Main Script
# -----------
# Set working directory (user should adjust this path)
# setwd('path/to/your/project/directory')
setwd('I:\\biocon\\Emmanuel_Oceguera\\projects\\2024_03_EruoMammalsConnectivity\\outcomes\\prediction_test_03072024')

# Load and preprocess data
# User should replace the path placeholders with actual paths to their data
# studyarea <- st_read("path/to/NC_Countries_reproj.shp")
# spgrid <- vect("path/to/sp_occu_20190829.shp")
studyarea = st_read("I:\\biocon\\Jeremy_Dertien\\NaturaConnect_Analysis\\Polygons\\NC_Countries_reproj.shp")
spgrid = vect("I:\\biocon\\Emmanuel_Oceguera\\projects\\2024_01_Mammals_species_distribution_DarwingCore\\data\\occurrances_data\\sp_occu_20190829.shp")

# Extract genus part from column names
genus_names <- names(spgrid)[4:20] %>% 
  str_extract(pattern = "(?<=_)[a-z]{3}") %>%
  na.omit() %>%
  unique()


# Isolate species by genus and sum into one column per genus into a list
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

# Convert all the spatVect to sf object and set the name for each object in the list  # nolint
sf_eu_mammals_lst <- lapply(eu_mammals_by_genus, sf::st_as_sf)
names(sf_eu_mammals_lst) <- paste0("m_", genus_names, "spp")

# Select the only presence squares
sf_eu_mammals_pres_lst <- lapply(sf_eu_mammals_lst, function(sf_obj){
  #Fiterl rows where any numeric column has a values of 1
  sf_obj %>% 
    filter(if_any(where(is.numeric), ~ .x == 1))
})

# Create points for the presence squares
eu_mammals_pres_pts_lst <- lapply(sf_eu_mammals_pres_lst, st_centroid_within_poly) # nolint: line_length_linter.

# Create buffer with different buffer distances 
# Define the buffer distance
buffer_distances <- c("m_canspp" = 100000, "m_ursspp" = 50000, "m_lynspp" = 60000, # nolint: line_length_linter.
                      "m_gulspp" = 40000, "m_alcspp" = 40000, "m_capspp" = 30000, # nolint: line_length_linter.
                      "m_cerspp" = 40000, "m_rupspp" = 40000, "m_susspp" = 70000, # nolint: line_length_linter.
                      "m_bisspp" = 50000, "m_ranspp" = 60000)

# Create the buffer taking into acoount the buffer distances
eu_mammals_pres_pts_buff_lst <- lapply(names(eu_mammals_pres_pts_lst), function(x){ # nolint: line_length_linter, paren_body_linter.
  sf_object <- eu_mammals_pres_pts_lst[[x]]
  presence_column <- x # Assing the ley in the list somehow indicates the presence column # nolint: line_length_linter.
  st_buffer_presence_squeres(sf_object,presence_column, buffer_distances)
})

# Assigning objects names
names(eu_mammals_pres_pts_buff_lst) <- paste0("m_", genus_names, "spp")

# Erase the presence squares to obtain a fully dissolved single polygon on the outside # nolint: line_length_linter.
# Create a list to store it
eu_mammals_buff_diff_lst <- list()
# Loop through the names of the first list
for (name in names(eu_mammals_pres_pts_buff_lst)){
  # Chek the names 
  if (name %in% names(sf_eu_mammals_pres_lst)) {
    # If true Perfome the erase function
    eu_mammals_buff_diff_lst[[name]] <- st_erase(eu_mammals_pres_pts_buff_lst[[name]], sf_eu_mammals_pres_lst[[name]])  # nolint: line_length_linter.
  } else {
    # If not true we set NULL
    eu_mammals_buff_diff_lst[[name]] <- NULL
    print(paste("No matchinf sf object for:", name))
  }
}

# Clip the available buffer area by the nations
eu_mammals_buff_diff_ctry_lst <- list()
# Loop throught each sf object
for(name in names(eu_mammals_buff_diff_lst)){
  # Intersection with the nations "study area"
  eu_mammals_buff_diff_ctry_lst[[name]] <- st_intersection(eu_mammals_buff_diff_lst[[name]], studyarea) # nolint: line_length_linter.
}

# Clip the buffer area with the inical grid that contains o and 1
eu_mammals_buff_diff_ctry_clip_lst <- list() # nolint: object_length_linter.
# Loop to each object in the list
for(name in names(eu_mammals_buff_diff_ctry_lst)){
  # Check if 'name' exists in inical grid sf_eu_mammals_lst
  if(name %in% names(sf_eu_mammals_lst)) {
    eu_mammals_buff_diff_ctry_clip_lst[[name]] <- st_intersection(sf_eu_mammals_lst[[name]], eu_mammals_buff_diff_ctry_lst[[name]]) # nolint
  } else {
    warning(paste(name, 'not found in sf_eu_mammals_lst')) # nolint
  }
}

# Create point from available squares
eu_mammals_avail_pts <- lapply(eu_mammals_buff_diff_ctry_clip_lst, function(sf_object){ # nolint
  # Select the first 4 columns along with the geometry column
  selected_sf_object <- sf_object  %>%
    select(1:4, geometry)
  st_centroid_within_poly(selected_sf_object)
})

eu_mammals_avail_selected_pts <- lapply(eu_mammals_avail_pts, function(sf_object){ # nolint
  # Filter rows, where any numeric columns has a 0 values
  sf_object %>%
    filter(if_any(where(is.numeric), ~ .x == 0))
})

# Final result combine the presence and available points
eu_mammals_pres_avail_pts <- list()
for (name in names(eu_mammals_pres_pts_lst)){
  # Retrive the matching sf objects from both list
  sf_objec1 <- eu_mammals_pres_pts_lst[[name]]
  sf_objec2 <- eu_mammals_avail_selected_pts[[name]]
  # Combine the objects
  combined_sf <- rbind(sf_objec1, sf_objec2)
  eu_mammals_pres_avail_pts[[name]] <- combined_sf 
}

## glmmm and predictions
# covariants <- rast("path/to/covariants.grd")
# biogeo <- rast("path/to/biogeographic.tif")
covariants = rast("I:\\biocon\\Jeremy_Dertien\\NaturaConnect_Analysis\\Covariates\\LULC_ZMeans\\Kilo1\\UngCovsRastStack2.grd")
biogeo = rast("I:\\biocon\\Jeremy_Dertien\\NaturaConnect_Analysis\\Covariates\\LULC_ZMeans\\Kilo1\\Biogeographic2016.tif")

# Resample Biogeo
biogeo <- resample(biogeo, covariants, method = "bilinear")

# Add Biogeo to the covariants and rename it
mammals_covariants <- c(covariants, biogeo) %>%
  rename(biogeo = short_name)
# Extracting covariants values for all at each point
eu_mammals_covs_values <- lapply(eu_mammals_pres_avail_pts, function(sf_object){
  # Convert sf object to SpatVector for compatibility with terra::extract()
  spat_vector <- vect(sf_object)
  #Extract covariate values for points in the SpatVector
  df <- terra::extract(mammals_covariants, spat_vector, bind = TRUE) |>
    as.data.frame()
  # if there are some Na we replace them
  df <- df %>% 
    mutate(across(where(is.numeric), ~replace_na(., 0)))
  return(df)
})

### glmm 
## Model fitting and selection (using parallel processing) for all the species 
library(parallel)
# Determinate the number of cores
no_cores <- detectCores() - 1  # Leave one core free
# Initiate cluster
cl <- makeCluster(no_cores)
# Export the fit_species_model function and any required data/objects to cluster
clusterExport(cl, varlist = c("fit_best_model_for_species", "eu_mammals_covs_values"))
clusterEvalQ(cl, {
  library(lme4)
  library(MuMIn)
})
# Prepare a list of species names to iterate over
species_names <- names(eu_mammals_covs_values)
# Run models in parallel
species_models <- parLapply(cl, species_names, function(species_name) {
  species_data <- eu_mammals_covs_values[[species_name]]
  fit_best_model_for_species(species_data, species_name)
})
stopCluster(cl)
# Set names to the list of models
names(species_models) <- species_names
# Check all the models in a loop (Optional)
# for(i in seq_along(species_models)) {
#   cat("Summary for model", names(species_models)[i], ":\n")
#   print(summary(species_models[[i]]$best_model))
#   cat("\n\n")
# }

## Predictions with covariants 
# Before ruinning, we have to alight the rasters covariants to the covariants we used in the model
# Create a raste stack
mammals_covariants_pred_grid <- mammals_covariants 

# Set the names of covariantes
model_covariants <- c("Dvlpd_ZMean", "Forest_ZMean", "Maxtmp10k_resamp", "RdDen10k_Zmean", "Shrub10k_resamp",
                      "TRI10km_reSamp", "Wland_ZMean", "biogeo")
# subset only the covariantes we've used in the model
prediction_grid_sub <- terra::subset(mammals_covariants_pred_grid, model_covariants)

### I think is not necesary
# # Convert to a data frame
# prediction_grid_sub_df <- terra::as.data.frame(prediction_grid_sub, xy = TRUE, na.rm = TRUE)
# # convert back to a rast, although I think we don't have to do this.
# eu_mammals_covs <- terra::rast(prediction_grid_sub_df, type="xyz")

# spp predictions in a loop
spp_predition_rastPaths <- list()
for (species_name in names(species_models)) {
  fitted_model <- species_models[[species_name]]
  spp_predition_rastPaths[[species_name]] <- predict_for_species(fitted_model, species_name, prediction_grid_sub)
}



