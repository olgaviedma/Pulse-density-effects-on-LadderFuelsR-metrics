
#**Impact of LiDAR Pulse Density on Forest Fuels Metrics derived using LadderFuelsR**
  
#**Authors: Olga Viedma and JM Moreno**

#This study evaluates the impact of LiDAR pulse density thinning on forest structure characterization using the LadderFuelsR framework. High-resolution LiDAR data were collected from Mediterranean forest sites with diverse structures and fire histories and systematically thinned to simulate varying pulse densities. Key forest metrics, including leaf area density (LAD), leaf area index (LAI), canopy base height, depth, height of fuel layers, and inter-layer distances, were derived to assess changes at different thinning levels.
#This first section deals with pre-processing normalized height LiDAR data: thinning, canopy height models (CHM), tree segmentation, crown metrics, cropping .las files and, calculating LAD profiles. 

# Getting Started

## Installation

remotes::install_github("Jean-Romain/lidRplugins")

#The CRAN version:
install.packages("LadderFuelsR")

# The development version:
install.packages("remotes")
library(remotes)
install_github("https://github.com/olgaviedma/LadderFuelsR", dependencies = TRUE)

install.packages("rgdal", dependencies = TRUE)
library(remotes)
library(rgdal)
install_github("r-spatial/sf")


## Required libraries

library(devtools)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(stringi)
library(purrr)
library(rlang)
library(tidyverse)
library(magrittr)
library(rlist)

library(sp)
library(sf)
library(raster)
library(data.table)
library(rgdal)

library(lidR)
library(leafR)
library(lidRplugins)
library(LadderFuelsR)


#SECTION 1. THINNIG, CHM, SEGMENTATION, LAD PROFILES

############################
#0. READ ZENODO FILES
############################

## María Olga Viedma Sillero. (2025). DATA.zip. Zenodo. Available at: https://zenodo.org/record/14786024

# Load necessary library
library(httr)

# Step 1: Define URLs and Local Paths
zenodo_url <- "https://zenodo.org/record/14786024/files/DATA.zip?download=1"  # Zenodo download link
zip_file <- "DATA.zip"       # Name of the downloaded ZIP file
output_dir <- "DATA"         # Final directory for extracted files
temp_dir <- "temp_extracted" # Temporary folder for intermediate extraction

# Step 2: Download the ZIP File (if not already downloaded)
if (!file.exists(zip_file)) {
  message("Downloading DATA.zip from Zenodo...")
  GET(zenodo_url, write_disk(zip_file, overwrite = TRUE))
} else {
  message("ℹ️  DATA.zip already exists. Skipping download.")
}

# Step 3: Validate the Downloaded ZIP File
file_size <- file.info(zip_file)$size
if (is.na(file_size) || file_size < 1e6) {  # Ensure the file is >1 MB
  stop("Download failed: ZIP file is too small or empty!")
} else {
  message(paste("ZIP file size:", round(file_size / (1024 * 1024), 2), "MB"))
}

# Step 4: Extract the ZIP File to a Temporary Folder
message("Extracting the ZIP file to a temporary folder...")
unzip(zip_file, exdir = temp_dir, overwrite = TRUE)

# Step 5: Move Contents to Final Output Directory
if (file.exists(file.path(temp_dir, "DATA"))) {
  # If ZIP contains an extra "DATA" folder, move its contents
  message("Moving extracted contents to the DATA directory...")
  file.rename(from = file.path(temp_dir, "DATA"), to = output_dir)
} else {
  # Otherwise, directly move the extracted files
  message("Moving extracted files to the DATA directory...")
  file.rename(from = temp_dir, to = output_dir)
}

# Step 6: Clean Up Temporary Files
unlink(temp_dir, recursive = TRUE)  # Delete the temporary folder

# Step 7: Verify the Final Extracted Files
message("Veifying the final structure of extracted files...")
extracted_files <- list.files(output_dir, recursive = TRUE)
if (length(extracted_files) == 0) {
  stop("Unzipping failed: No files were extracted.")
} else {
  message("Extraction successful! Files extracted:")
  print(extracted_files)
}

# Step 8 (Optional): Print a Preview of Extracted Files
message("First few extracted files:")
print(head(extracted_files))

############################
#0. FUNCTION CROWN METRICS
############################

custom_crown_metrics <- function(z, i) { # user-defined function
  metrics <- list(
    dz = 1, 
    th = 1,
    z_max = max(z),   # max height
    z_min = min(z),   # min height
    z_mean = mean(z),   # mean height
    z_sd = sd(z), # vertical variability of points
    z_q1=quantile(z, probs = 0.01),
    z_q5=quantile(z, probs = 0.05),
    z_q25=quantile(z, probs = 0.25),
    z_q50=quantile(z, probs = 0.50),
    z_q75=quantile(z, probs = 0.75),
    z_q95=quantile(z, probs = 0.95),
    crr=(mean(z)-min(z))/(max(z)-min(z))
  )
  return(metrics) # output
}
ccm = ~custom_crown_metrics(z = Z, i = Intensity)

############################
#1.THINNING LIDARPOD FILES (CAUTION THE THINNED DATA CAN BE CHANGE IF YOU USE LIDR OR LASTOOLS AND IF USE -HIGHEST or -LOWEST parameter in LasTools)
############################

### USING LASTOOLS (in a loop for checking thinning process consistency)

# Add LAStools bin directory to the PATH
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "C:/LAStools/bin", sep = ";"))

# Directories
las_dir <- "DATA"  
las_folders <- list.dirs(las_dir, full.names = TRUE, recursive = FALSE)
las_subfolders <- list.dirs(las_folders, full.names = TRUE, recursive = FALSE)

output_base_dir <-las_dir

# Thinning densities
target_densities <- c(100, 50, 25, 10, 5, 4, 3, 2, 1)

# Number of iterations
num_iterations <- 5

# Similarity threshold settings
cor_threshold <- 0.99  # High correlation implies consistency
rmse_threshold <- 0.1  # Adjust based on acceptable variation

# Function to compute similarity between two point clouds
compute_similarity <- function(file1, file2) {
  library(lidR)  # Ensure lidR is installed for point cloud processing
  pc1 <- readLAS(file1)
  pc2 <- readLAS(file2)
  
  if (is.empty(pc1) || is.empty(pc2)) {
    return(NA)  # Avoid errors on empty files
  }
  
  # Compute density histogram
  dens1 <- density(pc1@data$Z, n = 100)  # Elevation density
  dens2 <- density(pc2@data$Z, n = 100)
  
  # Compute correlation
  cor_value <- cor(dens1$y, dens2$y, use = "complete.obs")
  
  # Compute RMSE
  rmse_value <- sqrt(mean((dens1$y - dens2$y)^2))
  
  return(list(correlation = cor_value, RMSE = rmse_value))
}

# Loop through each subfolder
for (i in seq_along(las_subfolders)) {
  folder_path <- las_subfolders[i]
  site_name <- basename(folder_path)
  
  # List LAS files in the current folder
  las_list <- list.files(folder_path, pattern = glob2rx("*.las$"), full.names = TRUE)
  
  for (las_file in las_list) {
    for (density in target_densities) {
      # Calculate step size for thinning
      step <- sqrt(1 / density)
      
      for (iteration in 1:num_iterations) {
        # Define output directory and file
        thinning_folder <- file.path(output_base_dir, paste0("THINNED_", density, "P"), site_name)
        if (!dir.exists(thinning_folder)) dir.create(thinning_folder, recursive = TRUE)
        
        output_file <- file.path(thinning_folder, paste0(site_name, "_thinned_", density, "point_iter", iteration, ".las"))
        
        # Construct LAStools lasthin command
        cmd <- sprintf(
          'lasthin -i "%s" -o "%s" -step %.3f -random -seed 232', #lowest point
          las_file, output_file, step
        )
        
        # Debugging: Print command
        cat("Executing command:", cmd, "\n")
        
        # Execute the command
        tryCatch({
          result <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
          cat("Command output:\n", result, "\n")
        }, error = function(e) {
          cat("Error executing command for file:", las_file, "\n", e$message, "\n")
        })
        
        # Compute similarity if at least two iterations exist
        if (iteration > 1) {
          prev_file <- file.path(thinning_folder, paste0(site_name, "_thinned_", density, "point_iter", iteration - 1, ".las"))
          if (file.exists(prev_file)) {
            similarity <- compute_similarity(prev_file, output_file)
            cat(sprintf("Iteration %d vs %d for %s (Density %d): Correlation = %.3f, RMSE = %.3f\n",
                        iteration - 1, iteration, site_name, density, similarity$correlation, similarity$RMSE))
          }
        }
      }
      
      ### **Check Consistency Before Deleting Extra Iterations**
      # List all iteration files using correct pattern
      iteration_files <- list.files(thinning_folder, pattern = paste0("_thinned_", density, "point_iter[1-5]\\.las$"), full.names = TRUE)
      
      # Debug: Print detected files
      cat("Detected iteration files in", thinning_folder, ":", paste(iteration_files, collapse = ", "), "\n")
      
      if (length(iteration_files) < 2) {
        cat("Skipping", thinning_folder, "- Not enough iterations found.\n")
        next
      }
      
      # Check similarity for all iterations
      consistent <- TRUE  # Assume consistency unless proven otherwise
      for (iter in 2:5) {
        # Find the exact file names
        file1 <- list.files(thinning_folder, pattern = paste0("_thinned_", density, "point_iter", iter - 1, "\\.las$"), full.names = TRUE)
        file2 <- list.files(thinning_folder, pattern = paste0("_thinned_", density, "point_iter", iter, "\\.las$"), full.names = TRUE)
        
        if (length(file1) > 0 && length(file2) > 0) {
          similarity <- compute_similarity(file1[1], file2[1])  # Use first match
          
          if (!is.na(similarity$correlation) && !is.na(similarity$RMSE)) {
            cat(sprintf("Checking Iteration %d vs %d for %s (Density %d): Correlation = %.3f, RMSE = %.3f\n",
                        iter - 1, iter, site_name, density, similarity$correlation, similarity$RMSE))
            
            # If any iteration falls below threshold, mark as inconsistent
            if (similarity$correlation < cor_threshold || similarity$RMSE > rmse_threshold) {
              cat("WARNING: Inconsistent thinning detected in", thinning_folder, "\n")
              consistent <- FALSE
              break  # Stop checking if inconsistency is found
            }
          }
        } else {
          cat(sprintf("Skipping consistency check for %s, iteration %d: File not found.\n", thinning_folder, iter))
        }
      }
      
      # If all iterations are consistent, delete iterations 2-5
      if (consistent) {
        files_to_remove <- list.files(thinning_folder, pattern = paste0("_thinned_", density, "point_iter[2-5]\\.las$"), full.names = TRUE)
        
        if (length(files_to_remove) > 0) {
          file.remove(files_to_remove)
          cat("Deleted extra iterations in:", thinning_folder, "\n")
        }
      } else {
        cat("Skipping deletion in", thinning_folder, "- Iterations were not consistent.\n")
      }
    }
  }
}


############################
#2. CHM FOR LAS FILES
############################

library(raster)
library(lidR)

# Create corresponding structure in chm_dir

# Define LAS and crowns directories
las_dir <- "DATA"
chm_dir <- "CHM"

# List subfolders in LAS directory
las_folders <- list.dirs(las_dir, full.names = TRUE, recursive = FALSE)
las_subfolders <- list.dirs(las_folders, full.names = TRUE, recursive = FALSE)

# Remove trailing slash from crowns_dir (if exists)
chm_dir <- gsub("/$", "", chm_dir)

# Create corresponding structure in crowns_dir
for (las_folder in las_subfolders) {
  # Generate relative path
  relative_path <- gsub(paste0("^", las_dir, "/?"), "", las_folder)
  
  # Create corresponding folder in chm_dir
  chm_folder <- file.path(chm_dir, relative_path)
  dir.create(chm_folder, recursive = TRUE, showWarnings = FALSE)
}

# Check created folder structure
cat("Folders successfully created in:", chm_dir, "\n")


for (i in seq_along(las_subfolders)) {
  
  las_list <- list.files(las_subfolders[i], pattern = "*.las", full.names = TRUE, recursive = FALSE)
  
  for (j in seq_along(las_list)) {
    
    las_file <- lidR::readLAS(las_list[j])
    
    chm_pitfree <- grid_canopy(las_file, res = 0.5, pitfree(c(0, 2, 5, 10, 15, 20, 25, 30), c(0, 1.5), subcircle = 0.15))
    chm_pitfree[chm_pitfree > 40] <- NA
    chm_pitfree[chm_pitfree < 0] <- 0
    
    output_file <- paste0(chm_dir, gsub(las_dir, "", tools::file_path_sans_ext(las_list[j])), "_CHM05m.tif")
    
    # Write the current raster to its corresponding output file
    writeRaster(chm_pitfree, filename = output_file, format = "GTiff", overwrite = TRUE)
    
  }
}

############################
#3. CROWNS WATERSHED 0.25 FOR UNTHINNED LAS FILES
############################

library(lidRplugins)
library(lidR)
library(raster)

#########  CREATE NEW FOLDER FOR METRIC

las_dir <- "DATA/NO_THINNED"
las_folders <- list.dirs(las_dir, full.names = TRUE, recursive = FALSE)
las_list <- list.files(las_folders, pattern = "*.las", full.names = TRUE)

chm_dir <- "CHM/NO_THINNED"
chm_folders <- list.dirs(chm_dir, full.names = TRUE, recursive = FALSE)
chm_list <- list.files(chm_folders, pattern = glob2rx("*.tif"), full.names = TRUE)

crowns_dir <- "WATER_LAS/NO_THINNED"
# Remove trailing slash if exists
crowns_dir <- gsub("/$", "", crowns_dir)


# Create the same folder structure in crowns_dir1
for (las_folder in las_folders) {
  relative_path <- gsub(las_dir, "", las_folder)  # Get the relative path
  # Remove leading slash if exists
  relative_path <- gsub("^/", "", relative_path)
  crowns_folder <- file.path(crowns_dir, relative_path)  # Create corresponding folder in crowns_dir1
  dir.create(crowns_folder, recursive = TRUE, showWarnings = FALSE)  # Create folder if it doesn't exist
}

# CROWNS WATERSHED 0.25

for (i in seq_along(las_folders)) {
  lidar_list <- list.files(las_subfolders[i], pattern = glob2rx("*.las"), full.names = TRUE, recursive = FALSE)
  chm_list <- list.files(chm_folders[i], pattern = glob2rx("*.tif"), full.names = TRUE, recursive = FALSE)
  
  if (length(lidar_list) != length(chm_list)) {
    stop("Number of LAS files does not match number of CHM files.")
  }
  
  for (j in seq_along(lidar_list)) {
    # Load CHM and LAS files
    chm <- raster(chm_list[j])
    las_file <- lidR::readLAS(lidar_list[j])
    
    # Define the output file path
    output_file <- file.path(
      crowns_dir,
      paste0(gsub(las_dir, "", tools::file_path_sans_ext(lidar_list[j])), "_wat025.las")
    )
    
    # Apply watershed segmentation algorithm
    algo_watershed <- watershed(chm, th_tree = 4, tol = 0.25, ext = 1)
    crowns_water_las <- segment_trees(las_file, algo_watershed, attribute = "treeID", uniqueness = "incremental")
    
    # Filter points with valid treeID
    crowns_water_las_filtered <- filter_poi(crowns_water_las, !is.na(treeID))
    
    # Write the segmented LAS file with treeID
    lidR::writeLAS(crowns_water_las_filtered, output_file)
    
    cat("Processed and saved segmented LAS:", output_file, "\n")
  }
}


############################
#4 CROWNS METRICS FOR NO THINNED LAS FILES 
############################

library(lidR)
library(sf)
library(dplyr)

########  CREATE NEW FOLDER FOR METRIC

las_dir <- "WATER_LAS/NO_THINNED"
las_folders <- list.dirs(las_dir, full.names = TRUE, recursive = FALSE)
las_list <- list.files(las_folders, pattern = "*.las", full.names = TRUE)

metrics_dir <- "WATER_LAS/METRICS"
# Remove trailing slash if exists
metrics_dir <- gsub("/$", "", metrics_dir)


# Create the same folder structure in metrics_dir
for (las_folder in las_folders) {
  
  relative_path <- gsub(las_dir, "", las_folder)  # Get the relative path
  # Remove leading slash if exists
  relative_path <- gsub("^/", "", relative_path)
  metrics_folder <- file.path(metrics_dir, basename(las_dir), relative_path)  # Create corresponding folder in crowns_dir1
  dir.create(metrics_folder, recursive = TRUE, showWarnings = FALSE)  # Create folder if it doesn't exist
}

########  CALCULATE CROWN METRICS

# Define directories
las_dir <- "WATER_LAS/NO_THINNED"
las_folders <- list.dirs(las_dir, full.names = TRUE, recursive = FALSE)
las_list <- list.files(las_folders, pattern = "*_wat025.las", full.names = TRUE)

# Define the metrics directory
metrics_dir <- "WATER_LAS/METRICS"

# Ensure the metrics directory exists
if (!dir.exists(metrics_dir)) dir.create(metrics_dir, recursive = TRUE)

# Loop through LAS files
for (j in seq_along(las_list)) {
  # Get the relative path of the LAS file (relative to its base directory)
  relative_path <- sub("^.*?WATER_LAS/", "", las_list[j])
  
  # Remove the file name from the relative path to get the subdirectory
  sub_dir <- dirname(relative_path)
  
  # Construct the full output directory in metrics_dir
  output_dir <- file.path(metrics_dir, sub_dir)
  
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Get the base file name (without extension) from the LAS file
  file_name <- basename(tools::file_path_sans_ext(las_list[j]))
  
  # Construct the output file path
  output_file <- file.path(output_dir, paste0(file_name, "_metrics.shp"))
  
  # Debugging: print the constructed path
  print(output_file)
  
  # Your existing code to process LAS files and write outputs
  las_file1 <- lidR::readLAS(las_list[j])
  las_file2 <- filter_poi(las_file1, !is.na(treeID))
  las_file3 <- filter_poi(las_file2, Z >= 1)
  
  metrics1 <- crown_metrics(las_file3, func = .stdtreemetrics, geom = "convex")
  crown_diam <- data.frame(sqrt(metrics1$convhull_area / pi) * 2)
  names(crown_diam) <- "crown_diam"
  metrics2 <- crown_metrics(las_file3, func = ccm, geom = "convex")
  metrics_all <- dplyr::bind_cols(list(metrics1, crown_diam, metrics2))
  metrics_all1 <- metrics_all[, c(1:4, 6, 10:21)]
  names(metrics_all1) <- c("treeID", "Z", "npoints", "convhull_area", "crown_diam", 
                           "z_max", "z_min", "z_mean", "z_sd", "z_q1", "z_q5", 
                           "z_q25", "z_q50", "z_q75", "z_q95", "crr", "geometry")
  metrics_all2 <- st_as_sf(metrics_all1)
  
  # Save the output as a shapefile
  st_write(metrics_all2, output_file, driver = 'ESRI Shapefile', append = FALSE)
  
  # Print a success message
  cat("Metrics written to:", output_file, "\n")
}


############################
#5. CROP ALL THINNED LAS FILES WITH TREE POLYGONS derived from BENCHMARK LAS FILE (UNTHINNED)
############################

library(lidR)
library(sf)

# Define directories for segmented crowns

tree_dir <- "WATER_LAS/METRICS/NO_THINNED"
tree_folders <- list.dirs(tree_dir, full.names = TRUE, recursive = FALSE)
tree_list <- list.files(tree_folders, pattern = "*_metrics.shp", full.names = TRUE)

las_dir <- "DATA"
las_folders <- list.dirs(las_dir, full.names = TRUE, recursive = FALSE)
las_subfolders <- list.dirs(las_folders, full.names = TRUE, recursive = FALSE)
lidar_list <- list.files(las_subfolders, pattern = "*.las", full.names = TRUE)

crowns_dir1 <- "WATER_LAS/CROWNS_LIDARPOD"
# Remove trailing slash if it exists
crowns_dir1 <- gsub("/$", "", crowns_dir1)

# Create the same folder structure in crowns_dir1
for (las_subfolder in las_subfolders) {
  # Get the relative path by removing the base directory (las_dir)
  relative_path <- gsub(paste0("^", las_dir, "/?"), "", las_subfolder)
  
  # Construct the corresponding folder path in crowns_dir1
  crowns_folder <- file.path(crowns_dir1, relative_path)
  
  # Create the directory if it doesn't exist
  if (!dir.exists(crowns_folder)) dir.create(crowns_folder, recursive = TRUE)
  
  # Debugging: Print created folders for verification
  cat("Created folder:", crowns_folder, "\n")
}


### CROPPING 

# Define directories
las_dir <- "DATA"
crowns_dir2 <- "WATER_LAS/CROWNS_LIDARPOD"
tree_dir <- "WATER_LAS/METRICS/NO_THINNED"

# List directories
las_folders <- list.dirs(las_dir, full.names = TRUE, recursive = FALSE)
las_subfolders <- list.dirs(las_folders, full.names = TRUE, recursive = FALSE)
crowns_folders <- list.dirs(crowns_dir2, full.names = TRUE, recursive = FALSE)
tree_subfolders <- list.dirs(tree_dir, full.names = TRUE, recursive = FALSE)

# List LAS files
las_files <- list.files(
  path = las_subfolders,
  pattern = glob2rx("*.las"),
  full.names = TRUE,
  recursive = TRUE
)

cat("Total LAS files:", length(las_files), "\n")

# Iterate through LAS files
for (las_file in las_files) {
  # Extract folder and site information
  thinning_folder <- basename(dirname(dirname(las_file)))  # Parent folder name
  las_folder <- basename(dirname(las_file))  # Subfolder name
  las_base <- tools::file_path_sans_ext(basename(las_file))
  site_prefix <- sub("_.*", "", las_base)  # Extract SITE prefix like "SITE*"
  
  # Match tree subfolder
  tree_subfolder <- tree_subfolders[basename(tree_subfolders) == site_prefix]
  
  if (length(tree_subfolder) == 0) {
    cat("No matching tree subfolder for SITE:", site_prefix, "\n")
    next
  }
  
  # List tree shapefiles in the subfolder
  tree_files <- list.files(
    path = tree_subfolder,
    pattern = glob2rx("*_metrics.shp"),
    full.names = TRUE,
    recursive = FALSE
  )
  
  if (length(tree_files) == 0) {
    cat("No tree shapefiles found in tree subfolder for SITE:", site_prefix, "\n")
    next
  }
  
  tree_file <- tree_files[1]
  cat("Matched LAS file:", las_file, "\n")
  cat("Matched Tree file:", tree_file, "\n")
  
  # Read LAS and shapefile
  las <- lidR::readLAS(las_file)
  if (is.null(las)) {
    cat("Error: Unable to read LAS file:", las_file, "\n")
    next
  }
  
  tree_data <- sf::st_read(tree_file)
  if (nrow(tree_data) == 0) {
    cat("Warning: Tree shapefile is empty:", tree_file, "\n")
    next
  }
  
  # Ensure the crown output folder exists
  crown_folder <- crowns_folders[basename(crowns_folders) == thinning_folder]
  if (length(crown_folder) == 0) {
    cat("No matching crown folder for LAS folder:", thinning_folder, "\n")
    next
  }
  
  site_crown_folder <- file.path(crown_folder[1], las_folder)
  if (!dir.exists(site_crown_folder)) {
    dir.create(site_crown_folder, recursive = TRUE)
    cat("Created SITE folder in crown directory:", site_crown_folder, "\n")
  }
  
  # Process each polygon in the shapefile
  trees_ID <- tree_data %>% dplyr::select(treeID)
  
  for (k in seq_len(nrow(trees_ID))) {
    tree_polygon <- trees_ID[k, ]
    clipped_las <- clip_roi(las, tree_polygon)
    
    if (!is.null(clipped_las) && nrow(clipped_las@data) > 0) {
      # Extract the actual treeID for naming
      tree_id <- as.character(tree_polygon$treeID)
      
      # Create output filename using treeID
      output_filename <- file.path(
        site_crown_folder,
        paste0(tree_id, "_", thinning_folder, "_", site_prefix, "_CROWNS_LIDARPOD.las")
      )
      
      # Save the clipped LAS file
      lidR::writeLAS(clipped_las, file = output_filename)
      cat("Clipped LAS file saved:", output_filename, "\n")
    } else {
      cat("Warning: No points found within the polygon for LAS file:", las_file, "\n")
    }
  }
  
}


#### EXTRACT SUMMARIES ABOUT POINT DENSITY FOR EACH TREE


# Define the paths

las_dir <- "WATER_LAS/CROWNS_LIDARPOD"
las_folders <- list.dirs(las_dir, full.names = TRUE, recursive = FALSE)
las_subfolders <- list.dirs(las_folders, full.names = TRUE, recursive = FALSE)

density_list <- list()

# Loop through each subfolder
for (j in seq_along(las_subfolders)) {
  las_files <- list.files(las_subfolders[j], pattern = "\\.las$", full.names = TRUE, recursive = TRUE)
  las_folder_name <- basename(las_subfolders[j])
  
  # Create a data frame to store results for this subfolder
  density_results <- data.frame(File = character(), Density = numeric(), stringsAsFactors = FALSE)
  
  # Loop through LAS files and calculate density
  for (file in las_files) {
    las <- readLAS(file)
    
    if (is.empty(las)) {
      warning(paste("File", file, "is empty."))
      next
    }
    
    # Calculate point density
    density <- density(las)
    
    density_results <- rbind(density_results, data.frame(File = basename(file), Density = density))
  }
  
  # Save results or print them for the current subfolder
  if (nrow(density_results) > 0) {
    print(density_results)
    density_list[[length(density_list) + 1]] <- density_results
  }
}


# Combine all results into a single data frame
all_density_results <- do.call(rbind, density_list)

print("All density results combined")

# Filter out Inf and NaN values
cleaned_density_results <- all_density_results %>%
  filter(is.finite(Density))  # Keeps only finite numbers

# Save the combined results
write.csv(cleaned_density_results, file.path(las_dir, "all_density_results.csv"), row.names = FALSE)

summary_stats <- cleaned_density_results %>%
  summarise(
    Mean_Density = mean(Density, na.rm = TRUE),
    Median_Density = median(Density, na.rm = TRUE),
    SD_Density = sd(Density, na.rm = TRUE),
    Min_Density = min(Density, na.rm = TRUE),
    Max_Density = max(Density, na.rm = TRUE)
  )

print(summary_stats)

write.csv(summary_stats, file.path(las_dir, "summary_all_density_results.csv"), row.names = FALSE)


###############################################
#6. LAI-LAD METRICS BY TREE and point density statistics
###############################################

library(rlist)
library(lidR)
library(leafR)
library(dplyr)

# Define directories for LAD PROFILES

las_dir <- "DATA"
las_folders <- list.dirs(las_dir, full.names = TRUE, recursive = FALSE)
las_subfolders <- list.dirs(las_folders, full.names = TRUE, recursive = FALSE)
las_subfolders1 <- las_subfolders[!grepl("CROWNS_LIDARPOD", las_subfolders)]
lidar_list <- list.files(las_subfolders1, pattern = "*.las", full.names = TRUE)

crowns_dir1 <- "WATER_LAI"
# Remove trailing slash if it exists
crowns_dir1 <- gsub("/$", "", crowns_dir1)

# Create the same folder structure in crowns_dir1
for (las_subfolder in las_subfolders1) {
  # Get the relative path by removing the base directory (las_dir)
  relative_path <- gsub(paste0("^", las_dir, "/?"), "", las_subfolder)
  
  # Construct the corresponding folder path in crowns_dir1
  crowns_folder <- file.path(crowns_dir1, relative_path)
  
  # Create the directory if it doesn't exist
  if (!dir.exists(crowns_folder)) dir.create(crowns_folder, recursive = TRUE)
  
  # Debugging: Print created folders for verification
  cat("Created folder:", crowns_folder, "\n")
}


###  LAD PROFILES FOR EACH THINNED TREE 

crowns_dir1 <- "WATER_LAI"
crowns_folders <- list.dirs(crowns_dir1, full.names = TRUE, recursive = FALSE)
crowns_subfolders <- list.dirs(crowns_folders, full.names = TRUE, recursive = FALSE)  

las_dir <- "WATER_LAS/CROWNS_LIDARPOD"
las_folders <- list.dirs(las_dir, full.names = TRUE, recursive = FALSE)
las_subfolders <- list.dirs(las_folders, full.names = TRUE, recursive = FALSE)

for (i in seq_along(las_subfolders)) {
  
  las_list <- list.files(las_subfolders[i], pattern = glob2rx("*.las"), full.names = TRUE, recursive = TRUE)
  last_subfolder_name <- basename(las_subfolders[i])
  
  crowns_subfolders1 <-crowns_subfolders [i]
  
  files_with_more_than_10_points <- c()
  
  for (j in las_list) {
    las_data <- lidR::readLAS(j)
    las_data1 <- filter_poi(las_data, Z >= 1)
    
    if (!is.null(las_data1) && las_data1@header$`Number of point records` > 10) {
      files_with_more_than_10_points <- c(files_with_more_than_10_points, j)
    }
  }
  
  profile_list <- list()
  lidar_lai_list <- list()
  understory_lai_list <- list()
  LAHV_metric_list <- list()
  
  for (j in seq_along(files_with_more_than_10_points)) {
    tryCatch({
      las_file <- files_with_more_than_10_points[j]
      short_name_base <- stri_sub(las_file, 1, -5)
      short_name <- gsub(".*/", "", short_name_base)
      # Extract characters before the first underscore
      short_name1 <- gsub("^(.*?)_.*", "\\1", short_name)
      
      VOXELS_LAD <- leafR::lad.voxels(las_file, grain.size = 2)
      lad_profile <- leafR::lad.profile(VOXELS_LAD, relative = FALSE)
      lai_tot <- leafR::lai(lad_profile)
      understory_lai <- leafR::lai(lad_profile, min = 0.3, max = 2.5)
      LAHV_metric <- leafR::LAHV(lad_profile, LAI.weighting = FALSE, height.weighting = FALSE)
      
      lad_profile_df <- data.frame(lad_profile, treeID = as.factor(short_name))
      lai_tot_df <- data.frame(lai_tot, treeID = as.factor(short_name))
      understory_lai_df <- data.frame(understory_lai, ttreeID = as.factor(short_name))
      LAHV_metric_df <- data.frame(LAHV_metric, treeID = as.factor(short_name))
      
      profile_list[[j]] <- lad_profile_df
      lidar_lai_list[[j]] <- lai_tot_df
      understory_lai_list[[j]] <- understory_lai_df
      LAHV_metric_list[[j]] <- LAHV_metric_df
    }, error = function(e) {
      cat("Error processing file:", files_with_more_than_10_points[j], "\n")
    })
  }
  
  # Combine the lists into data frames
  profile_list1 <- do.call(rbind, profile_list)
  lidar_lai_list1 <- do.call(rbind, lidar_lai_list)
  understory_lai_list1 <- do.call(rbind, understory_lai_list)
  LAHV_metric_list1 <- do.call(rbind, LAHV_metric_list)
  
  # Debugging: Print the output_folder and file paths
  cat("Output Folder:", crowns_subfolders1, "\n")
  cat("File Paths:",
      file.path(crowns_subfolders1, paste0("alltrees_LAD_profile_voxels2m_", last_subfolder_name, ".txt")), "\n",
      file.path(crowns_subfolders1, paste0("alltrees_LAI_profile_voxels2m_", last_subfolder_name, ".txt")), "\n",
      file.path(crowns_subfolders1, paste0("alltrees_understory_profile_voxels2m_",last_subfolder_name, ".txt")), "\n",
      file.path(crowns_subfolders1, paste0("alltrees_LAHV_profile_voxels2m_", last_subfolder_name, ".txt")), "\n")
  
  # Adjusted the output file names to include the folder structure
  write.table(profile_list1, file = file.path(crowns_subfolders1, paste0("alltrees_LAD_profile_voxels2m_", last_subfolder_name , ".txt")), sep = "\t", row.names = FALSE)
  write.table(lidar_lai_list1, file =  file.path(crowns_subfolders1, paste0("alltrees_LAI_profile_voxels2m_", last_subfolder_name , ".txt")), sep = "\t", row.names = FALSE)
  write.table(understory_lai_list1, file = file.path(crowns_subfolders1, paste0("alltrees_understory_profile_voxels2m_",last_subfolder_name , ".txt")), sep = "\t", row.names = FALSE)
  write.table(LAHV_metric_list1, file = file.path(crowns_subfolders1, paste0("alltrees_LAHV_profile_voxels2m_", last_subfolder_name , ".txt")), sep = "\t", row.names = FALSE)
}

###############################################
#7.DEPURATING TREE LAD PROFILES (\>= 5 HEIGHT VALUES)
###############################################

LAD_dir <- "WATER_LAI"
LAD_folders <- list.dirs(LAD_dir, full.names = TRUE, recursive = FALSE)
LAD_subfolders <- list.dirs(LAD_folders, full.names = TRUE, recursive = FALSE)

for (k in seq_along(LAD_folders)) {
  
  LAD_folders1 <- list.dirs(LAD_folders[k], full.names = TRUE, recursive = FALSE)
  
  profiles_list3 <- list()
  
  for (m in seq_along(LAD_folders1)) {
    
    LAD_subfolders <- list.dirs(LAD_folders1[m], full.names = TRUE, recursive = TRUE)
    
    LAD_list <- list.files(LAD_subfolders, pattern = glob2rx("alltrees_LAD_profile_voxels2m_*.txt"), full.names = TRUE, recursive = TRUE)
    LAD_files <- lapply(LAD_list, function (X) read.table(X, sep = "\t", header = TRUE))
    
    profiles_list1 <- list()  
    
    for (i in seq_along(LAD_files)) {
      
      LAD_profiles <- LAD_files[[i]]
      LAD_profiles$treeID <- factor(LAD_profiles$treeID)
      
      trees_name1 <- as.character(LAD_profiles$treeID)
      trees_name2 <- factor(unique(trees_name1))
      
      profiles_list <- list()
      
      for (j in levels(trees_name2)) {
        tree2 <- LAD_profiles %>%
          filter(treeID == j) %>%
          mutate(lad = ifelse(is.na(lad), 0.0001, lad)) %>%
          filter(any(lad > 0.0001)) %>%  # Remove rows with lad values less than 0.0001
          filter(sum(lad != 0) > 0.0001)  # Keep trees with at least one non-zero LAD value
        
        profiles_list[[j]] <- tree2
      }
      
      profiles_list1 <- bind_rows(profiles_list)
      
      cases <- data.frame(table(profiles_list1$treeID))
      cases1 <- cases[cases$Freq >= 5, ]
      names(cases1) <- c("treeID", "Freq")
      
      profiles_list2 <- profiles_list1[profiles_list1$treeID %in% cases1$treeID, ]
      
    }
    
    profiles_list3[[m]] <- profiles_list2
    
    # Construct the output folder path based on the input folder path
    output_folder <- LAD_folders1[m]
    
    # Get short name from the folder path
    short_name2 <- basename(LAD_folders1[m])
    
    # Create the output folder if it doesn't exist
    if (!file.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    
    # Adjusted the output file names to include the folder structure
    write.table(profiles_list3[[m]], file = paste(output_folder, paste0("/", "depurated_LAD_profiles_", short_name2, ".txt"), sep = ""), sep = "\t", row.names = FALSE)
  }
}






