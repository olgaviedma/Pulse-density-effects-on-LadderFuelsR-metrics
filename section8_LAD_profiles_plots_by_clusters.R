
#**Impact of LiDAR Pulse Density on Forest Fuels Metrics derived using LadderFuelsR**

#**Authors: Olga Viedma and JM Moreno**

#This study evaluates the impact of LiDAR pulse density thinning on forest structure characterization using the LadderFuelsR framework. 
#High-resolution LiDAR data were collected from Mediterranean forest sites with diverse structures and fire histories and systematically thinned to simulate varying pulse densities. 
#Key forest metrics, including leaf area density (LAD), leaf area index (LAI), canopy base height, depth, height of fuel layers, and inter-layer distances, were derived to assess changes at different thinning levels.
#This section deals with plotting LAD profiles with LaddereFuelsR metrics.

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

#####################################################
#INPUTS FOR FIG. 6: PLOTS EFFECTIVE METRICS by CLUSTERS
#Plots of fuel layers with LAD percentage greater than 25 and the canopy base height (CBH) based on the maximum LAD percentage (get_plots_effective FUNCTION).
# Users can change the plotting function using other ones from LadderFuelsR
#####################################################

library(LadderFuelsR)
library(lattice)
library(dplyr)
library(ggplot2)
library(stringr)  # For string manipulation functions


cluster_dir <- "RESULTS/CLUSTER"
cluster_data <- read.table(file.path (cluster_dir, "ALL_METRICS_all_sites_CLUST_PCA4_CLUST_NT.txt"), sep = "\t", header = TRUE)

# Ensure factors are set correctly
cluster_data$treeID <- factor(cluster_data$treeID)
cluster_data$clust_nt <- factor(cluster_data$clust_nt)
cluster_data$cluster <- factor(cluster_data$cluster)
cluster_data$thin_level <- factor(cluster_data$thin_level)
cluster_data$sites <- factor(cluster_data$sites)

cluster_data$clust_change <- factor (paste0(cluster_data$clust_nt, "_",cluster_data$cluster))
cluster_data$clust_change <- factor(cluster_data$clust_change)

# Define the order of thin_level levels
new_order <- c("NO_THINNED",
               "THINNED_100P",
               "THINNED_50P",
               "THINNED_25P",
               "THINNED_10P",
               "THINNED_5P",
               "THINNED_4P",
               "THINNED_3P",
               "THINNED_2P",
               "THINNED_1P")
# Reorder factor levels
cluster_data$thin_level <- factor(cluster_data$thin_level, levels = new_order)


# PLOTTING 


library(ggplot2)
library(dplyr)

# Define directories
LAD_dir <- "WATER_LAI"
LAD_folders <- list.dirs(LAD_dir, full.names = TRUE, recursive = FALSE)
LAD_folders <- LAD_folders[!grepl("DATA_NO_THINNED", LAD_folders)]
LAD_subfolders <- list.dirs(LAD_folders, full.names = TRUE, recursive = FALSE)

figures_dir<-"FIGURES"
# Create the directory (including parent directories if they do not exist)
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
  cat("Directory created successfully:", figures_dir, "\n")
} else {
  cat("Directory already exists:", figures_dir, "\n")
}


# Main loop for saving plots
for (index in 1:nrow(cluster_data)) {
  
  # Extract necessary information
  current_treeID <- cluster_data$treeID[index]
  current_site <- cluster_data$sites[index]  
  current_cluster <- cluster_data$clust_change[index]
  current_thinning <- cluster_data$thin_level[index]
  
  # Print log for treeID consistency status
  cat("Processing treeID:", current_treeID, "with clust_change:", current_cluster, "in thinning level:", current_thinning, "\n")
  
  # Loop through each thinning level folder
  for (folder in LAD_folders) {
    
    # Check if the folder name matches the current thinning level
    if (grepl(current_thinning, folder, ignore.case = TRUE)) {
      
      # List subdirectories to find site matches
      LAD_subfolders <- list.dirs(folder, full.names = TRUE, recursive = FALSE)
      matched_subfolders <- LAD_subfolders[grepl(current_site, LAD_subfolders, ignore.case = TRUE)]
      
      # If there's a matching site folder, proceed
      if (length(matched_subfolders) > 0) {
        
        for (matched_folder in matched_subfolders) {
          print(paste("Matched site folder:", matched_folder, "for treeID:", current_treeID))
          
          # List the files for depurated LAD profiles and fuels LAD
          LAD_list <- list.files(matched_folder, pattern = glob2rx("depurated_LAD_profiles_*.txt"), full.names = TRUE)
          fuels_LAD_list <- list.files(matched_folder, pattern = glob2rx("7_fuels_lad_gt10perc_*.txt"), full.names = TRUE, ignore.case = TRUE)
          
          # Define the output directory based on consistent or individual cluster
          output_directory <- file.path(figures_dir, paste0("FIGURES_CLUST_", current_cluster, "_", current_thinning))
          
          # Create the output directory if it doesn't exist
          if (!dir.exists(output_directory)) {
            dir.create(output_directory, recursive = TRUE)
          }
          
          # Load and filter the data if files are found
          if (length(LAD_list) > 0 && length(fuels_LAD_list) > 0) {
            depurated_LAD_profiles <- lapply(LAD_list, function(X) read.table(X, sep = "\t", header = TRUE))
            fuels_LAD_files <- lapply(fuels_LAD_list, function(X) read.table(X, sep = "\t", header = TRUE))
            
            # Filter the data for the current treeID
            tree_data <- depurated_LAD_profiles[[1]] %>% dplyr::filter(treeID == current_treeID)
            fuels_data <- fuels_LAD_files[[1]] %>% dplyr::filter(treeID == current_treeID)
            
            if (nrow(tree_data) == 0 || nrow(fuels_data) == 0) {
              warning(paste("No data found for treeID:", current_treeID, "in thinning level:", current_thinning))
              next
            }
            
            # Generate plots for gaps and fbhs
            plots_trees_LAD <- get_plots_effective(tree_data, fuels_data, min_height = 1.5)
            
            # Save each plot immediately
            if (all(sapply(plots_trees_LAD, function(x) inherits(x, "ggplot")))) {
              for (name in names(plots_trees_LAD)) {
                plot <- plots_trees_LAD[[name]]
                if (inherits(plot, "ggplot")) {
                  # Define the output file path
                  output_file <- file.path(output_directory, paste0(name, ".png"))
                  
                  # Save the plot
                  tryCatch({
                    ggsave(output_file, plot = plot, width = 7, height = 7, units = "in")
                    cat("Saved plot for treeID:", name, "to", output_file, "\n")
                  }, error = function(e) {
                    warning(paste("Failed to save plot for", name, ":", e$message))
                  })
                } else {
                  warning(paste("Expected ggplot object but got", class(plot), "for treeID:", name))
                }
              }
            } else {
              warning(paste("Non-ggplot object found for treeID:", current_treeID))
            }
          } else {
            warning(paste("No files found for treeID:", current_treeID, "in folder:", matched_folder))
          }
        }
      }
    }
  }
}

```

