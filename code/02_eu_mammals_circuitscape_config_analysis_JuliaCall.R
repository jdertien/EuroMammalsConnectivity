
# Clear all objects from the current R environment 
rm(list = ls())
gc()

# Set the working directory to the main project folder.
# This should be the root directory where your data and scripts are organized.
setwd("path/to_your/project_folder")
setwd()

# Load the JuliaCall package for interfacing with Julia from R.
# Tools is used for handling file paths and names.
library(JuliaCall)
library(tools)
# install_julia() # only necessary if Julia is not installed. 

# Initialize Julia. Replace the path with your actual Julia installation directory.
julia_home <- "path/to_your/julia/bin/"
JuliaCall::julia_setup(julia_home)

# Install and load the Circuitscape package in Julia.
julia_command("using Pkg; Pkg.add(\"Circuitscape\")")
julia_library("Circuitscape")

# Define paths to the input and output directories.
# resistance_maps_dir should contain ASCII grid files (.asc) representing resistance maps.
resistance_maps_dir <- "path/to_your/input/resistance_maps"
output_dir <- "path/to_your/output"
config_dir <- "path/to_your/config_files"

# List all resistance map files in the specified directory.
resistance_files <- list.files(resistance_maps_dir, pattern = "\\ .asc$", full.names = TRUE)
# Extract species names from the file names by removing file extensions and unwanted suffixes.
species_names <- sapply(resistance_files, function(x) tools::file_path_sans_ext(basename(x)))

# Store paths to configuration files to be used later.
config_files <- list()

# Generate configuration files for each species based on their specific resistance maps.
for(i in seq_along(resistance_files)) {
  species_name <- gsub("_resistance", "", species_names[[i]])
  sp_resistance_name <- resistance_files[i]
  output_file_path <- sprintf("%s/%s_circuitscape.out", output_dir, species_name)
  config_file_path <- sprintf("%s/%s.ini", config_dir, species_name)
  config_files[[species_name]] <- config_file_path
  
  # Configure Circuitscape options using species-specific details.
  config_options <- sprintf("
[Calculation options]
low_memory_mode = False
parallelize = True
max_parallel = 8
solver = cg+amg
print_timings = True
preemptive_memory_release = True
print_rusages = False

[Output options]
set_null_currents_to_nodata = True
set_focal_node_currents_to_zero = False
set_null_voltages_to_nodata = True
compress_grids = False
write_volt_maps = False
write_cur_maps = False
output_file = %s
write_cum_cur_map_only = True
log_transform_maps = False
write_max_cur_maps = True

[Version]
version = 5.0.0

[Logging Options]
log_level = critical
log_file = None
profiler_log_file = None
screenprint_log = False

[Options for pairwise and one-to-all and all-to-one modes]
included_pairs_file = path/to_your/node_pairs_file.txt
use_included_pairs = True
point_file = path/to_your/points_file.asc

[Connection scheme for raster habitat data]
connect_using_avg_resistances = True
connect_four_neighbors_only = False

[Habitat raster or graph]
habitat_map_is_resistances = True
habitat_file = %s

[Circuitscape mode]
data_type = raster
scenario = pairwise
", output_file_path, sp_resistance_name)

  # Save the configuration options to a file
  if (!file.exists(config_file_path)) {
    writeLines(config_options, config_file_path)
  }
  cat("Config file created at:", config_file_path, "
")
}

# Run Circuitscape for each species using their respective configuration files.
for (config_file_path in config_files) {
  julia_command <- sprintf("Circuitscape.compute(\"%s\")", config_file_path)
  cat("Executing command in Julia: ", julia_command, "
")
  tryCatch({   # Execute the command in Julia with tryCatch to handle potential errors
    julia_eval(julia_command)
  }, error = function(e) {
    cat("Error in executing Circuitscape for config:", config_file_path, " Error Message:", e$message, "
")
  })
}
