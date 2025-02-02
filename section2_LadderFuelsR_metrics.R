
#**Impact of LiDAR Pulse Density on Forest Fuels Metrics derived using LadderFuelsR**

#**Authors: Olga Viedma and JM Moreno**

#This study evaluates the impact of LiDAR pulse density thinning on forest structure characterization using the LadderFuelsR framework. High-resolution LiDAR data were collected from Mediterranean forest sites with diverse structures and fire histories and systematically thinned to simulate varying pulse densities. Key forest metrics, including leaf area density (LAD), leaf area index (LAI), canopy base height, depth, height of fuel layers, and inter-layer distances, were derived to assess changes at different thinning levels.
#This second section deals with the calculation of  LadderFuelsR metrics on the LAD profiles. 

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

#########################################
#7.LAD PERCENTAGE OF EACH HEIGHT VALUE
#########################################

library(dplyr)

LAD_dir <- "WATER_LAI"
LAD_folders <- list.dirs(LAD_dir, full.names = TRUE, recursive = FALSE)

for (k in seq_along(LAD_folders)) {
  
  LAD_folders1 <- list.dirs(LAD_folders[k], full.names = TRUE, recursive = FALSE)
  
  gaps_perc2_list2 <- list()  # Initialize outside the loop
  short_name2 <- NULL
  
  for (m in seq_along(LAD_folders1)) {
    
    LAD_subfolders <- list.dirs(LAD_folders1[m], full.names = TRUE, recursive = TRUE)
    
    LAD_list <- list.files(LAD_subfolders, pattern = glob2rx("depurated_LAD_profiles_*.txt"), full.names = TRUE, recursive = TRUE)
    LAD_files <- lapply(LAD_list, function (X) read.table(X, sep = "\t", header = TRUE))
    
    gaps_perc2_list1 <- list()  # Initialize inside the loop
    
    for (j in seq_along(LAD_files)) {
      LAD_profiles <- LAD_files[[j]]
      LAD_profiles$treeID <- factor(LAD_profiles$treeID)
      
      trees_name1 <- as.character(LAD_profiles$treeID)
      trees_name2 <- factor(unique(trees_name1))
      
      gaps_perc2_list <- list()
      
      for (i in levels(trees_name2)) {
        tree1 <- LAD_profiles |> filter(treeID == i)
        percentiles <- calculate_gaps_perc(tree1)
        gaps_perc2_list[[i]] <- percentiles
      }
      
      gaps_perc2 <- bind_rows(gaps_perc2_list)
      gaps_perc2_list1[[j]] <- gaps_perc2
    }
    
    # Construct the output folder path based on the input folder path
    output_folder <- LAD_folders1[m]
    
    # Get short name from the folder path
    short_name2 <- basename(LAD_folders1[m])
    
    # Create the output folder if it doesn't exist
    if (!file.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    gaps_perc2_list2[[m]] <- do.call(rbind, gaps_perc2_list1)
    
    # Adjusted the output file names to include the folder structure
    write.table(gaps_perc2_list2[[m]], file = paste(output_folder, paste0("/", "1a_percentile_metrics_", short_name2, ".txt"), sep = ""), sep = "\t", row.names = FALSE)
  }
}


head(gaps_perc2_list2[[1]])


#########################################
#8.GAPS AND FUEL LAYERS BASE HEIGHT (FBH)
#########################################

LAD_dir <- "WATER_LAI"
LAD_folders <- list.dirs(LAD_dir, full.names = TRUE, recursive = FALSE)

for (k in seq_along(LAD_folders)) {
  
  LAD_folders1 <- list.dirs(LAD_folders[k], full.names = TRUE, recursive = FALSE)
  
  metrics_all_percentile_list1 <- list()
  short_name2 <- NULL
  
  for (m in seq_along(LAD_folders1)) {
    
    LAD_subfolders <- list.dirs(LAD_folders1[m], full.names = TRUE, recursive = TRUE)
    
    LAD_list <- list.files(LAD_subfolders, pattern = glob2rx("depurated_LAD_profiles_*.txt"), full.names = TRUE, recursive = TRUE)
    LAD_files <- lapply(LAD_list, function (X) read.table(X, sep = "\t", header = TRUE))
    
    
    metrics_all_percentile_list <- list()
    tree_metrics_filtered <- NULL
    
    for (i in seq_along(LAD_files)) {
      
      depurated_LAD_profiles <- LAD_files[[i]]
      depurated_LAD_profiles$treeID <- factor(depurated_LAD_profiles$treeID)
      
      trees_name1 <- as.character(depurated_LAD_profiles$treeID)
      trees_name2 <- factor(unique(trees_name1))
      
      metrics_percentile_list <- list()
      
      for (j in levels(trees_name2)) {
        tree2 <- depurated_LAD_profiles %>% dplyr::filter(treeID == j)
        metrics_percentile <- get_gaps_fbhs(tree2, step=1, min_height=1.5, perc_gap= 25, perc_base= 25, verbose=TRUE)
        metrics_percentile_list[[j]] <- metrics_percentile
      }
      
      metrics_all_percentile <- dplyr::bind_rows(metrics_percentile_list)
      
      # Remove rows with all NA or zero values
      tree_metrics_filtered <- metrics_all_percentile[rowSums(is.na(metrics_all_percentile) | metrics_all_percentile == 0, na.rm = TRUE) < ncol(metrics_all_percentile), ]
      
    }
    
    # Construct the output folder path based on the input folder path
    output_folder <- LAD_folders1[m]
    
    # Get short name from the folder path
    short_name2 <- basename(LAD_folders1[m])
    
    # Create the output folder if it doesn't exist
    if (!file.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    metrics_all_percentile_list1[[m]] <- tree_metrics_filtered
    
    # Adjusted the output file names to include the folder structure
    write.table(metrics_all_percentile_list1[[m]], file = paste(output_folder, paste0("/", "1_gaps_fbhs_metrics_", short_name2, ".txt"), sep = ""), sep = "\t", row.names = FALSE)
  }
}


head(metrics_all_percentile_list1[[1]])

#########################################
#9.DISTANCE BETWEEN FUEL LAYERS
#########################################

LAD_dir <- "WATER_LAI"
LAD_folders <- list.dirs(LAD_dir, full.names = TRUE, recursive = FALSE)


for (k in seq_along(LAD_folders)) {
  
  LAD_folders1 <- list.dirs(LAD_folders[k], full.names = TRUE, recursive = FALSE)
  
  metrics_distance_list2 <- list()
  short_name2<-NULL  
  
  for (m in seq_along(LAD_folders1)) {
    
    LAD_subfolders <- list.dirs(LAD_folders1[m], full.names = TRUE, recursive = TRUE)
    
    gaps_fbhs_list <- list.files(LAD_subfolders, pattern = glob2rx("1_gaps_fbhs_metrics_*.txt"), full.names = TRUE, recursive = FALSE)
    gaps_fbhs_files <- lapply(gaps_fbhs_list, function (X) read.table(X, sep = "\t", header = TRUE))
    
    gaps_perc_list <- list.files(LAD_subfolders, pattern = glob2rx("1a_percentile_metrics_*.txt"), full.names = TRUE, ignore.case = TRUE)
    gaps_perc_files <- lapply(gaps_perc_list, function (X) read.table(X, sep = "\t", header = TRUE))
    
    
    metrics_distance_list1 <- list()
    
    for (i in seq_along(gaps_fbhs_files)) {
      
      gaps_fbh1 <- gaps_fbhs_files[[i]]
      gaps_fbh1$treeID <- factor(gaps_fbh1$treeID)
      
      gaps_perc1 <- gaps_perc_files[[i]]
      gaps_perc1$treeID <- factor(gaps_perc1$treeID)
      
      trees_name1 <- as.character(gaps_fbh1$treeID)
      trees_name2 <- factor(unique(trees_name1))
      
      metrics_distance_list <- list()
      
      for (j in levels(trees_name2)) {
        # Filter data for each tree
        tree2 <- gaps_fbh1 %>% dplyr::filter(treeID == j)
        tree3 <- gaps_perc1 %>% dplyr::filter(treeID == j)
        
        # Check if there is data for the tree
        if (nrow(tree2) == 0 || nrow(tree3) == 0) {
          cat("No data for tree:", j, "\n")
          metrics_distance_list[[j]] <- NULL
          next  # Skip to the next iteration if there is no data
        }
        
        # Get distance metrics for each tree
        metrics_distance <- get_distance(tree2,tree3, min_height=1.5, step=1, verbose = TRUE)
        metrics_distance_list[[j]] <- metrics_distance
      }
      
      # Combine the individual data frames
      metrics_distance_all <- dplyr::bind_rows(metrics_distance_list)
      #metrics_distance_all <- metrics_distance_all[, order(names(metrics_distance_all))]
      metrics_distance_list1[[i]] <- metrics_distance_all
    }
    
    # Construct the output folder path based on the input folder path
    output_folder <- LAD_folders1[m]
    
    # Get short name from the folder path
    short_name2 <- basename(LAD_folders1[m])
    
    # Create the output folder if it doesn't exist
    if (!file.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    metrics_distance_list2[[m]] <- do.call(rbind, metrics_distance_list1)
    
    # Adjusted the output file names to include the folder structure
    write.table(metrics_distance_list2[[m]], file = paste(output_folder, paste0("/", "2_distance_metrics_", short_name2, ".txt"), sep = ""), sep = "\t", row.names = FALSE)
  }
}

head(metrics_distance_list2[[1]])


#########################################
#10.FUEL LAYERS DEPTH
#########################################

library(dplyr)
library(magrittr)

LAD_dir <- "WATER_LAI"
LAD_folders <- list.dirs(LAD_dir, full.names = TRUE, recursive = FALSE)


for (k in seq_along(LAD_folders)) {
  
  LAD_folders1 <- list.dirs(LAD_folders[k], full.names = TRUE, recursive = FALSE)
  
  metrics_depth_list2 <- list()
  short_name2<-NULL
  
  for (m in seq_along(LAD_folders1)) {
    
    LAD_subfolders <- list.dirs(LAD_folders1[m], full.names = TRUE, recursive = TRUE)
    
    LAD_list <- list.files(LAD_subfolders, pattern = glob2rx("depurated_LAD_profiles_*.txt"), full.names = TRUE, recursive = FALSE)
    LAD_files <- lapply(LAD_list, function (X) read.table(X, sep = "\t", header = TRUE))
    
    distance_list <- list.files(LAD_subfolders, pattern = glob2rx("2_distance_metrics_*.txt"), full.names = TRUE, ignore.case = TRUE)
    distance_files <- lapply(distance_list, function (X) read.table(X, sep = "\t", header = TRUE))
    
    
    metrics_depth_list1 <- list()
    
    for (i in seq_along(distance_files)) {
      
      distances <- distance_files[[i]]
      distances$treeID <- factor(distances$treeID)
      
      depurated_LAD_profiles <-LAD_files[[i]]
      depurated_LAD_profiles$treeID <- factor(depurated_LAD_profiles$treeID)
      
      trees_name1 <- as.character(distances$treeID)
      trees_name2 <- factor(unique(trees_name1))
      
      metrics_depth_list <- list()
      
      for (j in levels(trees_name2)){
        
        tree1 <- depurated_LAD_profiles |> dplyr::filter(treeID == j)
        tree2 <- distances |> dplyr::filter(treeID == j)
        
        # Get depths for each tree
        metrics_depth <- get_depths(tree1, tree2, step= 1,min_height= 1.5, verbose=TRUE)
        metrics_depth_list[[j]] <- metrics_depth
      }
      
      # Combine the individual data frames
      metrics_all_depth <- dplyr::bind_rows(metrics_depth_list)
      metrics_all_depth <- metrics_all_depth[, order(names(metrics_all_depth))]
      
      metrics_depth_list1[[i]] <- metrics_all_depth
      
    }
    
    # Construct the output folder path based on the input folder path
    output_folder <- LAD_folders1[m]
    
    # Get short name from the folder path
    short_name2 <- basename(LAD_folders1[m])
    
    # Create the output folder if it doesn't exist
    if (!file.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    metrics_depth_list2[[m]] <- do.call(rbind, metrics_depth_list1)
    
    # Adjusted the output file names to include the folder structure
    write.table(metrics_depth_list2[[m]], file = paste(output_folder, paste0("/", "3_depth_metrics_", short_name2, ".txt"), sep = ""), sep = "\t", row.names = FALSE)
  }
}


head(metrics_depth_list2[[1]])


#########################################
#11.FUEL LAYERS BASE HEIGHT (FBH) AFTER REMOVING DISTANCES = 1
#########################################

library(SSBtools)
library(dplyr)
library(magrittr)

LAD_dir <- "WATER_LAI"
LAD_folders <- list.dirs(LAD_dir, full.names = TRUE, recursive = FALSE)

for (k in seq_along(LAD_folders)) {
  
  LAD_folders1 <- list.dirs(LAD_folders[k], full.names = TRUE, recursive = FALSE)
  
  fbh_corr_list2 <- list()
  short_name2 <- NULL
  
  for (m in seq_along(LAD_folders1)) {
    
    LAD_subfolders <- list.dirs(LAD_folders1[m], full.names = TRUE, recursive = TRUE)
    
    depth_list <- list.files(LAD_subfolders, pattern = glob2rx("3_depth_metrics_*.txt"), full.names = TRUE, ignore.case = TRUE)
    depth_files <- lapply(depth_list, function(X) read.table(X, sep = "\t", header = TRUE))
    
    
    fbh_corr_list1 <- list()
    
    for (i in seq_along(depth_files)) {
      
      depths <- depth_files[[i]]
      depths$treeID <- factor(depths$treeID)
      
      trees_name1 <- as.character(depths$treeID)
      trees_name2 <- factor(unique(trees_name1))
      
      fbh_corr_list <- list()
      
      for (j in levels(trees_name2)){
        
        # Filter data for each tree
        tree3 <- depths |> dplyr::filter(treeID == j)
        fbh_corr <- get_real_fbh(tree3, step= 1, number_steps = 1, min_height=1.5,verbose = TRUE)
        fbh_corr_list[[j]] <- fbh_corr
      }
      
      # Combine fbh values for all trees
      fbh_corr_all <- dplyr::bind_rows(fbh_corr_list)
      fbh_corr_all$treeID <- factor(fbh_corr_all$treeID)
      
      # Reorder columns
      original_column_names <- colnames(fbh_corr_all)
      # Specify prefixes
      prefixes <- c("treeID", "Hdist", "Hcbh", "Hdepth", "dist", "depth", "max_height")
      # Initialize vector to store new order
      new_order <- c()
      
      # Loop over prefixes
      for (prefix in prefixes) {
        # Find column names matching the current prefix
        matching_columns <- grep(paste0("^", prefix), original_column_names, value = TRUE)
        # Append to the new order
        new_order <- c(new_order, matching_columns)
      }
      
      # Reorder values
      fbh_corr_all <- fbh_corr_all[, new_order]
      fbh_corr_list1[[i]] <- fbh_corr_all
      
    }
    
    # Construct the output folder path based on the input folder path
    output_folder <- LAD_folders1[m]
    
    # Get short name from the folder path
    short_name2 <- basename(LAD_folders1[m])
    
    # Create the output folder if it doesn't exist
    if (!file.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    fbh_corr_list2[[m]] <- do.call(rbind, fbh_corr_list1)
    
    # Adjusted the output file names to include the folder structure
    write.table(fbh_corr_list2[[m]], file = paste(output_folder, paste0("/", "4_fbh_metrics_corr_", short_name2, ".txt"), sep = ""), sep = "\t", row.names = FALSE)
  }
}

head(fbh_corr_list2[[1]])


#########################################
#12.FUEL LAYERS DEPTH AFTER REMOVING DISTANCES = 1
#########################################

library(dplyr)
library(magrittr)
library(tidyr)


LAD_dir <- "WATER_LAI"
LAD_folders <- list.dirs(LAD_dir, full.names = TRUE, recursive = FALSE)

for (k in seq_along(LAD_folders)) {
  
  LAD_folders1 <- list.dirs(LAD_folders[k], full.names = TRUE, recursive = FALSE)
  
  depth_corr_list2 <- list()
  short_name2 <- NULL
  
  for (m in seq_along(LAD_folders1)) {
    
    LAD_subfolders <- list.dirs(LAD_folders1[m], full.names = TRUE, recursive = TRUE)
    
    fbhcor_list <- list.files(LAD_subfolders, pattern = glob2rx("4_fbh_metrics_corr_*.txt"), full.names = TRUE, ignore.case = TRUE)
    fbhcor_files <- lapply(fbhcor_list, function(X) read.table(X, sep = "\t", header = TRUE))
    
    depth_corr_list1 <- list()
    
    for (i in seq_along(fbhcor_files)) {
      
      fbhcor <- fbhcor_files[[i]]
      fbhcor$treeID <- factor(fbhcor$treeID)
      
      trees_name1 <- as.character(fbhcor$treeID)
      trees_name2 <- factor(unique(trees_name1))
      
      depth_corr_list <- list()
      
      for (j in levels(trees_name2)){
        
        # Filter data for each tree
        tree3 <- fbhcor|> dplyr::filter(treeID == j)
        depth_corr<- get_real_depths(tree3, step=1, min_height=1.5)
        depth_corr_list[[j]] <- depth_corr
      }
      
      # Combine fbh values for all trees
      depth_corr_all <- dplyr::bind_rows(depth_corr_list)
      depth_corr_all$treeID <- factor(depth_corr_all$treeID)
      
      # Reorder columns
      original_column_names <- colnames(depth_corr_all)
      # Specify prefixes
      desired_order <- c("treeID", "Hcbh", "dptf", "dist", "Hdist", "Hdptf", "max_height")
      
      # Identify unique prefixes
      prefixes <- unique(sub("^([a-zA-Z]+).*", "\\1", original_column_names))
      # Initialize vector to store new order
      new_order <- c()
      
      # Loop over desired order of prefixes
      for (prefix in desired_order) {
        # Find column names matching the current prefix
        matching_columns <- grep(paste0("^", prefix), original_column_names, value = TRUE)
        # Append to the new order
        new_order <- c(new_order, matching_columns)
      }
      
      # Reorder values
      depth_corr_all <- depth_corr_all[, new_order]
      depth_corr_list1[[i]] <- depth_corr_all
      
    }
    
    # Construct the output folder path based on the input folder path
    output_folder <- LAD_folders1[m]
    
    # Get short name from the folder path
    short_name2 <- basename(LAD_folders1[m])
    
    # Create the output folder if it doesn't exist
    if (!file.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    depth_corr_list2[[m]] <- do.call(rbind, depth_corr_list1)
    
    # Adjusted the output file names to include the folder structure
    write.table(depth_corr_list2[[m]], file = paste(output_folder, paste0("/", "5_depth_metrics_corr_", short_name2, ".txt"), sep = ""), sep = "\t", row.names = FALSE)
  }
}

head(depth_corr_list2[[1]])


#########################################
#13.FUEL LAYERS DISTANCES (\> 1 M)
#########################################

library(dplyr)
library(magrittr)
library(stringr)

LAD_dir <- "WATER_LAI"
LAD_folders <- list.dirs(LAD_dir, full.names = TRUE, recursive = FALSE)

for (k in seq_along(LAD_folders)) {
  
  LAD_folders1 <- list.dirs(LAD_folders[k], full.names = TRUE, recursive = FALSE)
  
  dist_corr_list2 <- list()
  short_name2 <- NULL
  
  for (m in seq_along(LAD_folders1)) {
    
    LAD_subfolders <- list.dirs(LAD_folders1[m], full.names = TRUE, recursive = TRUE)
    
    depth_corr_list <-list.files(LAD_subfolders, pattern = glob2rx("5_depth_metrics_corr_*.txt"), full.names = TRUE, ignore.case = TRUE)
    depth_corr_files <- lapply(depth_corr_list, function(X) read.table(X, sep = "\t", header = TRUE))
    
    dist_corr_list1 <- list()
    
    for (i in seq_along(depth_corr_files)) {
      
      depthcor <- depth_corr_files[[i]]
      depthcor$treeID <- factor(depthcor$treeID)
      
      trees_name1 <- as.character(depthcor$treeID)
      trees_name2 <- factor(unique(trees_name1))
      
      dist_corr_list <- list()
      
      for (j in levels(trees_name2)){
        
        # Filter data for each tree
        tree3 <- depthcor |> dplyr::filter(treeID == j)
        dist_corr <- get_effective_gap(tree3,number_steps = 1,min_height=1.5)
        dist_corr_list[[j]] <- dist_corr
      }
      
      # Combine depth values for all trees
      dist_corr_all <- dplyr::bind_rows(dist_corr_list)
      # =======================================================================#
      # REORDER COLUMNS:
      # =======================================================================#
      # Get original column names
      original_column_names <- colnames(dist_corr_all)
      
      # Specify prefixes
      desired_order <- c("treeID", "Hcbh", "dptf","effdist","dist", "Hdist", "Hdptf", "max_")
      
      # Identify unique prefixes
      prefixes <- unique(sub("^([a-zA-Z]+).*", "\\1", original_column_names))
      # Initialize vector to store new order
      new_order <- c()
      
      # Loop over desired order of prefixes
      for (prefix in desired_order) {
        # Find column names matching the current prefix
        matching_columns <- grep(paste0("^", prefix), original_column_names, value = TRUE)
        # Append to the new order
        new_order <- c(new_order, matching_columns)
      }
      
      # Reorder values
      dist_corr_all1 <- dist_corr_all[, new_order]
      
      dist_corr_list1[[i]] <- dist_corr_all1
      
    }
    
    # Construct the output folder path based on the input folder path
    output_folder <- LAD_folders1[m]
    
    # Get short name from the folder path
    short_name2 <- basename(LAD_folders1[m])
    
    # Create the output folder if it doesn't exist
    if (!file.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    
    dist_corr_list2[[m]] <- dist_corr_list1
    
    # Adjusted the output file names to include the folder structure
    write.table(dist_corr_list2[[m]], file = paste(output_folder, paste0("/", "6_dist_metrics_corr_", short_name2, ".txt"), sep = ""), sep = "\t", row.names = FALSE)
  }
}

head(dist_corr_list2[[1]])


#########################################
#14.FUELS LAD PERCENTAGE (\> 10 %)
#########################################

library(dplyr)
library(magrittr)

LAD_dir <- "WATER_LAI"
LAD_folders <- list.dirs(LAD_dir, full.names = TRUE, recursive = FALSE)

for (k in seq_along(LAD_folders)) {
  
  LAD_folders1 <- list.dirs(LAD_folders[k], full.names = TRUE, recursive = FALSE)
  
  LAD_metrics1_list2 <- list()  
  LAD_metrics2_list2 <- list()  
  short_name2<-NULL 
  
  for (m in seq_along(LAD_folders1)) {
    
    LAD_subfolders <- list.dirs(LAD_folders1[m], full.names = TRUE, recursive = TRUE)
    LAD_list <- list.files(LAD_subfolders, pattern = "^depurated_.*\\.txt$", full.names = TRUE, ignore.case = TRUE) 
    LAD_files <- lapply(LAD_list, function(X) read.table(X, sep = "\t", header = TRUE))
    
    dist_corr_list <- list.files(LAD_subfolders, pattern = "^6_dist_metrics_corr_.*\\.txt$", full.names = TRUE, ignore.case = TRUE)
    dist_corr_files <- lapply(dist_corr_list, function(X) read.table(X, sep = "\t", header = TRUE))
    
    LAD_metrics1_list1 <- list()  
    LAD_metrics2_list1 <- list()  
    
    for (i in seq_along(dist_corr_files)) {
      
      distcor <- dist_corr_files[[i]]
      distcor$treeID <- factor(distcor$treeID)
      
      LAD_profile1 <- LAD_files[[i]]
      LAD_profile1$treeID <- factor(LAD_profile1$treeID)
      
      trees_name1 <- as.character(distcor$treeID)
      trees_name2 <- factor(unique(trees_name1))
      
      LAD_metrics1_list <- list()  
      LAD_metrics2_list <- list()  
      
      for (j in levels(trees_name2)) {
        # Filter data for each tree
        tree1<-LAD_profile1 |> dplyr::filter(treeID ==j)
        tree2 <- distcor |> dplyr::filter(treeID ==j)
        
        # Get LAD metrics for each tree
        LAD_metrics <- get_layers_lad (tree1, tree2, threshold=10, step = 1,min_height=1.5,verbose=TRUE)
        LAD_metrics1_list[[j]] <- LAD_metrics$df1
        LAD_metrics2_list[[j]] <- LAD_metrics$df2
      }
      
      LAD_metrics_all1 <- dplyr::bind_rows(LAD_metrics1_list)
      LAD_metrics_all2 <- dplyr::bind_rows(LAD_metrics2_list)
      
      # List of data frames
      LAD_metrics_list <- list(LAD_metrics_all1, LAD_metrics_all2)
      
      # Initialize an empty list to store reordered data frames
      reordered_list <- list()
      
      # Loop over each data frame
      for (h in seq_along(LAD_metrics_list)) {
        
        LAD_metrics_all <- LAD_metrics_list[[h]]
        
        # Get original column names
        original_column_names <- colnames(LAD_metrics_all)
        
        # Specify prefixes
        desired_order <- c("treeID", "Hcbh", "dptf","effdist", "Hdist", "Hdptf", "max_","nlayers")
        
        # Identify unique prefixes
        prefixes <- unique(sub("^([a-zA-Z]+).*", "\\1", original_column_names))
        # Initialize vector to store new order
        new_order <- c()
        
        # Loop over desired order of prefixes
        for (prefix in desired_order) {
          # Find column names matching the current prefix
          matching_columns <- grep(paste0("^", prefix), original_column_names, value = TRUE)
          # Append to the new order
          new_order <- c(new_order, matching_columns)
        }
        
        # Reorder columns
        LAD_metrics_all <- LAD_metrics_all[, new_order]
        # Store the reordered data frame in the list
        reordered_list[[h]] <- LAD_metrics_all
      }
      
      LAD_metrics1_list1[[i]] <- reordered_list[[1]]
      LAD_metrics2_list1[[i]] <- reordered_list[[2]]
      
      # Construct the output folder path based on the input folder path
      output_folder <- LAD_folders1[m]
      
      # Get short name from the folder path
      short_name2 <- basename(LAD_folders1[m])
      
      # Create the output folder if it doesn't exist
      if (!file.exists(output_folder)) {
        dir.create(output_folder, recursive = TRUE)
      }
      
      
      
      write.table(LAD_metrics1_list1[[i]] , file = paste(output_folder, paste0("/", "7_fuels_lad_all_",short_name2,".txt"), sep = ""), sep = "\t", row.names = FALSE)
      write.table(LAD_metrics2_list1[[i]] , file = paste(output_folder, paste0("/", "7_fuels_lad_gt10perc_",short_name2, ".txt"), sep = ""), sep = "\t", row.names = FALSE)
      
    }  }
}

head(LAD_metrics2_list1[[1]])


#########################################
#15. CBH BASED ON DIFFERENT CRITERIA: MAX LAD, MAXIMUM AND LAST DISTANCE (EFFECTIVE FUEL LAYERS)
#########################################

library(dplyr)
library(magrittr)

LAD_dir <- "WATER_LAI"
LAD_folders <- list.dirs(LAD_dir, full.names = TRUE, recursive = FALSE)

for (k in seq_along(LAD_folders)) {
  
  LAD_folders1 <- list.dirs(LAD_folders[k], full.names = TRUE, recursive = FALSE)
  
  cbh_dist_list2 <- list()
  short_name2 <- NULL
  
  
  for (m in seq_along(LAD_folders1)) {
    
    LAD_subfolders <- list.dirs(LAD_folders1[m], full.names = TRUE, recursive = TRUE)
    
    fuels_LAD_list <- list.files(LAD_subfolders, pattern = glob2rx("7_fuels_lad_gt10perc_*.txt"), full.names = TRUE, ignore.case = TRUE)
    fuels_LAD_files <- lapply(fuels_LAD_list, function(X) read.table(X, sep = "\t", header = TRUE))
    
    #effective_LAD<-tree3
    
    cbh_dist_list1 <- list()
    
    for (i in seq_along(fuels_LAD_files)) {
      
      fuels_lad1 <- fuels_LAD_files[[i]]
      fuels_lad1$treeID <- factor(fuels_lad1$treeID)
      
      trees_name1 <- as.character(fuels_lad1$treeID)
      trees_name2 <- factor(unique(trees_name1))
      
      cbh_dist_list <- list()
      
      for (j in levels(trees_name2)){
        
        # Filter data for each tree
        tree3 <- fuels_lad1 |> dplyr::filter(treeID == j)
        cbh_dist <- get_cbh_metrics (tree3, min_height= 1.5,hdepth1_height = 2.5, verbose=TRUE)
        
        # Check if the column 'maxlad1_Hcbh' exists in the cbh_dist data frame
        if ("maxlad1_Hcbh" %in% colnames(cbh_dist)) {
          print(paste("maxlad1_Hcbh column exists in tree:", j, "of file:", i))
          
          cbh_dist$maxlad_Hcbh <- ifelse(!is.na(cbh_dist$maxlad1_Hcbh), cbh_dist$maxlad1_Hcbh, cbh_dist$maxlad_Hcbh)
        }
        
        cbh_dist_list[[j]] <- cbh_dist
      }
      
      # Combine depth values for all trees
      cbh_dist_all <- dplyr::bind_rows(cbh_dist_list)
      
      # Get original column names
      original_column_names <- colnames(cbh_dist_all)
      
      # Specify prefixes
      desired_order <- c("treeID", "Hcbh", "dptf","effdist","dist", "Hdist", "Hdptf","maxlad_","maxlad1_","max_","last_","nlayers")
      
      # Identify unique prefixes
      prefixes <- unique(sub("^([a-zA-Z]+).*", "\\1", original_column_names))
      # Initialize vector to store new order
      new_order <- c()
      
      # Loop over desired order of prefixes
      for (prefix in desired_order) {
        # Find column names matching the current prefix
        matching_columns <- grep(paste0("^", prefix), original_column_names, value = TRUE)
        # Append to the new order
        new_order <- c(new_order, matching_columns)
      }
      
      # Reorder columns
      cbh_dist_all1 <- cbh_dist_all[, new_order]
      
      cbh_dist_list1[[i]] <- cbh_dist_all1
      
    }
    
    # Construct the output folder path based on the input folder path
    output_folder <- LAD_folders1[m]
    
    # Get short name from the folder path
    short_name2 <- basename(LAD_folders1[m])
    
    # Create the output folder if it doesn't exist
    if (!file.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    
    cbh_dist_list2[[m]] <- cbh_dist_list1
    
    # Adjusted the output file names to include the folder structure
    write.table(cbh_dist_list2[[m]], file = paste(output_folder, paste0("/", "8_cbhs_metrics_ladgt10_", short_name2, ".txt"), sep = ""), sep = "\t", row.names = FALSE)
  }
}

head(cbh_dist_list2[[1]])


#########################################
#16. CBH BASED ON THE BREAKING POINT METHOD AND LAD PERCENTAGE
#########################################

library(dplyr)
library(magrittr)

LAD_dir <- "WATER_LAI"
LAD_folders <- list.dirs(LAD_dir, full.names = TRUE, recursive = FALSE)

for (k in seq_along(LAD_folders)) {
  
  LAD_folders1 <- list.dirs(LAD_folders[k], full.names = TRUE, recursive = FALSE)
  
  cum_LAD_metrics_list2 <- list()
  short_name2 <- NULL
  
  for (m in seq_along(LAD_folders1)) {
    
    LAD_subfolders <- list.dirs(LAD_folders1[m], full.names = TRUE, recursive = TRUE)
    
    LAD_list <- list.files(LAD_subfolders, pattern = glob2rx("depurated_LAD_profiles_*.txt"), full.names = TRUE, recursive = FALSE)
    LAD_files <- lapply(LAD_list, function (X) read.table(X, sep = "\t", header = TRUE))
    
    cbhs_list <- list.files(LAD_subfolders, pattern = glob2rx("8_cbhs_metrics_ladgt10_*.txt"), full.names = TRUE, ignore.case = TRUE)
    cbhs_files <- lapply(cbhs_list, function(X) read.table(X, sep = "\t", header = TRUE))
    
    
    cum_LAD_metrics_list1 <- list()  
    
    for (i in seq_along(cbhs_files)) {
      
      cbh1 <- cbhs_files[[i]]
      cbh1$treeID <- factor(cbh1$treeID)
      
      LAD_profile1 <- LAD_files[[i]]
      LAD_profile1$treeID <- factor(LAD_profile1$treeID)
      
      trees_name1 <- as.character(cbh1$treeID)
      trees_name2 <- factor(unique(trees_name1))
      
      cum_LAD_metrics_list <- list()
      
      for (j in levels(trees_name2)) {
        # Filter data for each tree
        tree1 <- LAD_profile1 |> dplyr::filter(treeID == j)
        tree2 <- cbh1 |> dplyr::filter(treeID == j)
        
        # Get cumulative LAD metrics for each tree
        cum_LAD_metrics <- get_cum_break(tree1, tree2, threshold =75, min_height= 1.5,verbose=T)
        cum_LAD_metrics_list[[j]] <- cum_LAD_metrics
      }
      
      # Combine the individual data frames
      cum_LAD_metrics_all <- dplyr::bind_rows(cum_LAD_metrics_list) 
      
      # =======================================================================#
      # REORDER COLUMNS
      # =======================================================================#
      
      # Get original column names
      original_column_names <- colnames(cum_LAD_metrics_all)
      
      # Specify prefixes (adjust accordingly)
      prefixes <- c("treeID", "Hcbh", "below", "above", "bp", "max", "cumlad")
      
      # Initialize vector to store new order
      new_order <- c()
      
      # Loop over prefixes
      for (prefix in prefixes) {
        # Find column names matching the current prefix
        matching_columns <- grep(paste0("^", prefix), original_column_names, value = TRUE)
        
        # Extract numeric suffixes and order the columns based on these suffixes
        numeric_suffixes <- as.numeric(gsub(paste0("^", prefix), "", matching_columns))
        matching_columns <- matching_columns[order(numeric_suffixes)]
        
        # Append to new order
        new_order <- c(new_order, matching_columns)
      }
      
      # Reorder columns
      cum_LAD_metrics_order <- cum_LAD_metrics_all[, new_order]
      
      cum_LAD_metrics_list1[[i]] <- cum_LAD_metrics_order
    }
    
    # Construct the output folder path based on the input folder path
    output_folder <- LAD_folders1[m]
    
    # Get short name from the folder path
    short_name2 <- basename(LAD_folders1[m])
    
    # Create the output folder if it doesn't exist
    if (!file.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    
    cum_LAD_metrics_list2[[m]] <- cum_LAD_metrics_list1
    
    # Adjusted the output file names to include the folder structure
    write.table(cum_LAD_metrics_list2[[m]], file = paste(output_folder, paste0("/", "9_cbh_breaking_point_", short_name2, ".txt"), sep = ""), sep = "\t", row.names = FALSE)
  }
}

head(cum_LAD_metrics_list2[[1]]) 


#########################################
#17.LAI METRICS
#########################################

metrics_crowns_dir <- "WATER_LAI"
metrics_crowns_folders <- list.dirs(metrics_crowns_dir, full.names = TRUE, recursive = FALSE)

for (k in seq_along(metrics_crowns_folders)) {
  
  metrics_crowns_subfolders <- list.dirs(metrics_crowns_folders[k], full.names = TRUE, recursive = F)
  
  metrics_tot1 <- NULL
  short_name2<-NULL
  
  for (m in seq_along(metrics_crowns_subfolders)) {
    
    metrics_crowns_list <-list.files(metrics_crowns_subfolders[m], pattern = glob2rx("alltrees_*.txt"), full.names = TRUE, ignore.case = TRUE)
    exclude_files <- list.files(metrics_crowns_subfolders[m], pattern = "alltrees_LAD_.*\\.txt$", full.names = TRUE, ignore.case = TRUE)
    metrics_crowns_list1 <- setdiff(metrics_crowns_list, exclude_files)
    metrics_crowns_files <- lapply(metrics_crowns_list1, function(X) read.table(X,sep = "\t", header = TRUE))
    
    
    short_name2 <- gsub(".*/", "", metrics_crowns_subfolders[m])
    
    metrics_file1 <- metrics_crowns_files[[1]]
    metrics_file1$treeID <- factor(metrics_file1$treeID)
    
    metrics_file2 <- metrics_crowns_files[[2]]
    metrics_file2$treeID <- factor(metrics_file2$treeID)
    
    metrics_file3 <- metrics_crowns_files[[3]]
    metrics_file3$treeID <- factor(metrics_file3$ttreeID)
    
    # Use left_join to merge based on treeID1
    metrics_tot<- metrics_file1 |> left_join(metrics_file2, by = "treeID")
    metrics_tot1<- metrics_tot |> left_join(metrics_file3, by = "treeID")
    
    output_folder <- metrics_crowns_subfolders[m]
    # Debugging: Print the output_folder and file paths
    cat("Output Folder:", output_folder, "\n")
    cat("File Paths:",
        file.path(output_folder, paste0("LAI_metrics_", short_name2, ".txt")), "\n")
    
    write.table(metrics_tot1, file.path(output_folder, paste0("LAI_metrics_", short_name2, ".txt")), sep="\t", row.names = FALSE)
    
  }
}



#########################################
#18. JOINING LADDER FUEL PROPERTIES WITH CROWN POLYGONS
#########################################

LAD_dir <- "WATER_LAI"
LAD_folders <- list.dirs(LAD_dir, full.names = TRUE, recursive = FALSE)

crowns_dir <- "WATER_LAS/METRICS/NO_THINNED"
crown_folders <- list.dirs(crowns_dir, full.names = TRUE, recursive = FALSE)

for (k in seq_along(LAD_folders)) {
  
  LAD_subfolders <- list.dirs(LAD_folders[k], full.names = TRUE, recursive = FALSE)
  crown_folders1 <- list.dirs(crowns_dir[1], full.names = TRUE, recursive = FALSE)
  
  short_name3 <- basename(LAD_folders[k])
  
  merged_maps_list2 <- list()
  short_name2 <- NULL
  
  for (m in seq_along(LAD_subfolders)) {
    
    # List files for each type of data
    bp_list <- list.files(LAD_subfolders[m], pattern = glob2rx("9_cbh_breaking_point_*.txt"), full.names = TRUE, ignore.case = TRUE)
    cbh_list <- list.files(LAD_subfolders[m], pattern = glob2rx("8_cbhs_metrics_ladgt10_*.txt"), full.names = TRUE, ignore.case = TRUE)
    lai_list <- list.files(LAD_subfolders[m], pattern = glob2rx("LAI_metrics_*.txt"), full.names = TRUE, ignore.case = TRUE)
    crowns_list <- list.files(crown_folders1[m], pattern = glob2rx("*_metrics.shp"), full.names = TRUE, ignore.case = TRUE)
    
    # Skip subfolders with missing files
    if (length(bp_list) == 0 || length(cbh_list) == 0 || length(lai_list) == 0 || length(crowns_list) == 0) {
      cat("Skipping subfolder due to missing files:", LAD_subfolders[m], "\n")
      next
    }
    
    # Read data files
    bp_files <- lapply(bp_list, function(X) read.table(X, sep = "\t", header = TRUE))
    cbhs_files <- lapply(cbh_list, function(X) read.table(X, sep = "\t", header = TRUE))
    lai_files <- lapply(lai_list, function(X) read.table(X, sep = "\t", header = TRUE))
    crowns_files <- lapply(crowns_list, function(X) st_read(X))
    
    merged_maps_list1 <- list()
    
    for (i in seq_along(bp_list)) {
      
      # Process individual files
      cbhs_file1 <- cbhs_files[[i]]
      cbhs_file1$treeID1 <- factor(cbhs_file1$treeID1)
      
      bp_file1 <- bp_files[[i]]
      bp_file1$treeID1 <- factor(bp_file1$treeID1)
      
      lai_file1 <- lai_files[[i]]
      lai_file1$treeID1 <- gsub("^(.*?)_.*", "\\1", as.character(lai_file1$treeID))
      lai_file1$treeID1 <- factor(lai_file1$treeID1)
      
      crowns_files1 <- crowns_files[[i]]
      crowns_files1$treeID1 <- factor(crowns_files1$treeID)
      
      
      # Output folder and file writing
      output_folder <- LAD_subfolders[m]
      short_name2 <- basename(LAD_subfolders[m])
      
      if (!file.exists(output_folder)) {
        dir.create(output_folder, recursive = TRUE)
      }
      
      # Merge data files
      files_to_merge <- list(cbhs_file1, bp_file1, lai_file1)
      merge_two <- function(x, y) merge(x, y, by = "treeID1", all = TRUE)
      merged_data <- Reduce(merge_two, files_to_merge)
      
      # Remove unwanted columns
      cols_to_remove <- grep("^treeID.y|max_height.y|treeID.x|ttreeID", names(merged_data))
      merged_data1 <- merged_data[, -cols_to_remove]
      
      # Merge with crowns
      crowns_files2 <- crowns_files1 |> left_join(merged_data1, by = "treeID1")
      tree_cols <- grep("^tree", names(crowns_files2), value = TRUE)
      crowns_files3 <- crowns_files2[, -which(names(crowns_files2) %in% tree_cols)]
      
      trees <- st_drop_geometry(crowns_files2[, c(tree_cols)])
      names(trees) <- c("treeIDf", "treeID1", "treeID")
      crowns_files4 <- cbind(trees, crowns_files3)
      crowns_files4$thin_level <- rep(short_name3, nrow(crowns_files4))
      crowns_files4$sites <- rep(short_name2, nrow(crowns_files4))
      # Append to results list
      merged_maps_list1[[i]] <- crowns_files4
    }
    
    # Combine results for all files
    if (length(merged_maps_list1) > 0) {
      merged_maps_list2[[m]] <- do.call(rbind, merged_maps_list1)
    } else {
      cat("No valid data to merge for subfolder:", LAD_subfolders[m], "\n")
      next
    }
    
    
    st_write(merged_maps_list2[[m]], file.path(output_folder, paste0("ladderfuels_metrics_", short_name3, "_", short_name2, ".shp")), driver = "ESRI Shapefile", append = F)
  }
}


#########################################
#19.ALL METRICS and DIFFERENCES II: (NO THINNED - REST) 
#########################################

library(dplyr)
library(sf)

#### create the FOLDER: "DATA_NO_THINNED" AND COPY THE FOLDER "NO_THINNED" FROM "WATER_LAI"

# Define the source and destination paths
source_folder <- "WATER_LAI/NO_THINNED"
destination_folder <- "WATER_LAI/DATA_NO_THINNED"

# Create the destination folder if it does not exist
if (!dir.exists(destination_folder)) {
  dir.create(destination_folder, recursive = TRUE)
}

# List all files and directories within the source folder
files_and_dirs <- list.files(source_folder, full.names = TRUE, recursive = TRUE)

# Copy each file while maintaining the original structure
for (item in files_and_dirs) {
  # Compute the relative path from the source folder
  relative_path <- sub(paste0("^", source_folder, "/"), "", item)
  
  # Define the destination path while preserving structure
  destination_path <- file.path(destination_folder, "NO_THINNED", relative_path)
  
  # Check if it is a directory or a file
  if (dir.exists(item)) {
    # Create directory if it does not exist
    if (!dir.exists(destination_path)) {
      dir.create(destination_path, recursive = TRUE)
    }
  } else {
    # Create the parent directory if it does not exist
    parent_dir <- dirname(destination_path)
    if (!dir.exists(parent_dir)) {
      dir.create(parent_dir, recursive = TRUE)
    }
    
    # Copy file
    file.copy(item, destination_path, overwrite = TRUE)
  }
}

# Print success message
cat("Folder 'NO_THINNED' successfully copied to 'DATA_NO_THINNED' with structure preserved.\n")

###################3
extract_site_number <- function(folder_path) {
  site_match <- regexpr("SITE\\d+", folder_path)
  if (site_match == -1) {  # Check if site_match is -1
    return(NA)
  } else {
    site_string <- substr(folder_path, site_match, site_match + attr(site_match, "match.length") - 1)
    site_number <- as.numeric(gsub("\\D", "", site_string))
    return(site_number)
  }
}

crown_dir<- "WATER_LAI"

metrics_crowns_dir <- "WATER_LAI"
metrics_crowns_folders <- list.dirs(metrics_crowns_dir, full.names = TRUE, recursive = FALSE)
metrics_crowns_subfolders <- list.dirs(metrics_crowns_folders, full.names = TRUE, recursive = FALSE)
metrics_crowns_subfolders1 <- metrics_crowns_subfolders[!grepl("DATA_NO_THINNED", metrics_crowns_subfolders)]

metrics_NO_THINNED_dir <- "WATER_LAI/DATA_NO_THINNED/NO_THINNED"
metrics_NO_THINNED_folders <- list.dirs(metrics_NO_THINNED_dir, full.names = TRUE, recursive = FALSE)  # Corrected

metrics_crowns_subfolders <- metrics_crowns_subfolders[order(sapply(metrics_crowns_subfolders, extract_site_number))]
metrics_NO_THINNED_folders <- metrics_NO_THINNED_folders[order(sapply(metrics_NO_THINNED_folders, extract_site_number))]


## DIFFERENCES BETWEEN NO-THINNED LEVEL AND THE REST OF THINNING LEVELS FOR ALL SITES

# Initialize lists to store results
all_differences_list <- list()
all_differences_list1 <- list()
alldata1_list<- list()
alldata2_list<- list()

# Helper function to extract SITE name from file paths
get_site_name <- function(path) {
  basename(dirname(path))  # Adjust this if the SITE name is in a different part of the path
}

# Create the list of all METRICS_CROWNS files and group them by SITE
metrics_crowns_list <- list.files(metrics_crowns_subfolders1, 
                                  pattern = glob2rx("ladderfuels_metrics_*.shp"), 
                                  full.names = TRUE, 
                                  recursive = TRUE)

# Group METRICS_CROWNS files by SITE
metrics_crowns_grouped <- split(metrics_crowns_list, sapply(metrics_crowns_list, get_site_name))

# Create the list of all NO_THINNED files
metrics_NO_THINNED_list <- list.files(metrics_NO_THINNED_dir, 
                                      pattern = glob2rx("ladderfuels_metrics_*.shp"), 
                                      full.names = TRUE, 
                                      recursive = TRUE)

# Loop through each NO_THINNED file
for (m in seq_along(metrics_NO_THINNED_list)) {
  # Get the current NO_THINNED file and its SITE
  NO_THINNED_file <- metrics_NO_THINNED_list[m]
  NO_THINNED_site <- get_site_name(NO_THINNED_file)
  
  # Fetch all METRICS_CROWNS files for this SITE
  crowns_files_for_site <- metrics_crowns_grouped[[NO_THINNED_site]]
  
  if (is.null(crowns_files_for_site)) {
    message("No matching METRICS_CROWNS files found for SITE: ", NO_THINNED_site)
    next
  }
  
  # Load the NO_THINNED file
  NO_THINNED_metrics <- st_read(NO_THINNED_file, quiet = TRUE)
  NO_THINNED_metrics$treeID1 <- factor(NO_THINNED_metrics$treeID1)
  
  # Compare with each METRICS_CROWNS file in the site
  for (crowns_file in crowns_files_for_site) {
    # Load the crowns file
    crowns_metrics <- st_read(crowns_file, quiet = TRUE)
    crowns_metrics$treeID1 <- factor(crowns_metrics$treeID1)
    
    outdir <- dirname(crowns_file)
    thin_level <- basename(dirname(dirname(crowns_file)))
    
    # Find common variables
    common_vars <- intersect(names(crowns_metrics), names(NO_THINNED_metrics))
    
    # Subset to common variables
    crowns_metrics_common <- crowns_metrics[, common_vars, drop = FALSE]
    NO_THINNED_common <- NO_THINNED_metrics[, common_vars, drop = FALSE]
    
    # Filter out rows with NA in Hcbh1
    crowns_metrics_common <- crowns_metrics_common[!is.na(crowns_metrics_common$Hcbh1), ]
    NO_THINNED_common <- NO_THINNED_common[!is.na(NO_THINNED_common$Hcbh1), ]
    
    # Calculate differences
    data1 <- st_drop_geometry(crowns_metrics_common)
    data2 <- st_drop_geometry(NO_THINNED_common)
    
    data1$Hcbh_br[is.na(data1$Hcbh_br)] <- data1$bp_Hcbh[is.na(data1$Hcbh_br)]
    data2$Hcbh_br[is.na(data2$Hcbh_br)] <- data2$bp_Hcbh[is.na(data2$Hcbh_br)]
    
    data1$treeID1 <-as.character(data1$treeID1)
    data2$treeID1<- as.character(data2$treeID1)
    
    common_tree_ids <- intersect(data1$treeID1, data2$treeID1)
    data1 <- data1[data1$treeID1 %in% common_tree_ids, ]
    data2 <- data2[data2$treeID1 %in% common_tree_ids, ]
    
    
    # Define tree-related columns to exclude
    tree_columns <- c("treeID", "treeIDf", "treeID1", "thn_lvl", "sites")
    numeric_columns <- setdiff(names(data1), tree_columns)
    
    data1_numeric <- data1[, numeric_columns, drop = FALSE]
    data2_numeric <- data2[, numeric_columns, drop = FALSE]
    
    
    # Ensure data frames have matching columns
    if (!all(names(data1_numeric) == names(data2_numeric))) {
      stop("Column names in data1_numeric and data2_numeric do not match!")
    }
    
    # Calculate differences directly
    differences <- data1_numeric - data2_numeric
    
 
    # Combine with tree-related columns
    tree_columns_data <- data1[, tree_columns, drop = FALSE]
    differences_tot <- cbind(tree_columns_data, differences)
    
    
    # Define the breaks
    breaks_height <- c(-Inf, -3, -2, -1,0,1,2,3, Inf)
    
    # Use cut() to create a factor variable based on the breaks
    differences_tot$height_factor <- cut(differences_tot$mx_hgh_, breaks = breaks_height, labels = c("<= -3", "-3 to -2", "-2 to -1","-1 to 0",
                                                                                                     "0 to 1","1 to 2","2 to 3","> 3"))
    
    # Subset the dataframe to keep only rows with absolute difference less than 3
    differences_tot1 <- differences_tot[abs(differences_tot$mx_hgh_) < 3, ]
    differences_tot1_height_factor <- differences_tot1[,c("treeID","height_factor")]
    
    crowns_metrics$treeID <- as.character(crowns_metrics$treeID)
    differences_tot1_height_factor$treeID <- as.character(differences_tot1_height_factor$treeID)
    
   
    differences_shapes<-merge(crowns_metrics_common[,c("treeID")], differences_tot, by="treeID", all = TRUE)   
    differences_shapes1<-merge(crowns_metrics_common[,c("treeID")], differences_tot1, by="treeID", all = TRUE)
    
    st_write(differences_shapes, file.path(outdir, paste0("ladderfuels_diff_NO_THINNED_", thin_level,"_",NO_THINNED_site, "_tot.shp")), driver = "ESRI Shapefile", append=F)
    st_write(differences_shapes1, file.path(outdir, paste0("ladderfuels_diff_NO_THINNED_", thin_level,"_",NO_THINNED_site, "_depurat.shp")), driver = "ESRI Shapefile", append=F)
    
    write.table(differences_tot, file.path(outdir, paste0("ALL_DIFFERENCES_NO_THINNED_", thin_level,"_", NO_THINNED_site, ".txt")), sep = "\t", row.names = FALSE) 
    write.table(differences_tot1, file.path(outdir, paste0("ALL_DIFFERENCES_NO_THINNED_", thin_level,"_", NO_THINNED_site, "_depurat.txt")), sep = "\t", row.names = FALSE) 
    
    # Store results
    all_differences_list[[length(all_differences_list) + 1]] <- differences_tot
    all_differences_list1[[length(all_differences_list1) + 1]] <- differences_tot1
    
    
    data1_numeric$thin_level <- rep(thin_level, nrow(data1_numeric))
    data1_numeric$sites <- rep(NO_THINNED_site, nrow(data1_numeric))
    alldata <- cbind(tree_columns_data, data1_numeric)
    
    trees1 <- unique(alldata$treeID1)
    trees2 <- unique(differences_tot$treeID1)
    common_tree_ids <- intersect(trees1, trees2)
    all_data_wide1 <- alldata[alldata$treeID1 %in% common_tree_ids, ]
    
    trees3 <- unique(differences_tot1$treeID1)
    common_tree_ids1 <- intersect(trees1, trees3)
    all_data_wide2 <- alldata[alldata$treeID1 %in% common_tree_ids1, ]
    
    all_data_wide1 <-st_drop_geometry(all_data_wide1)
    all_data_wide2<-st_drop_geometry(all_data_wide2)
    
    write.table(all_data_wide1, file.path(outdir,paste0("ALL_METRICS_together_long_NO_THINNED_", thin_level,"_",NO_THINNED_site, "_tot.txt")), sep = "\t", row.names = FALSE)
    write.table(all_data_wide2, file.path(outdir,paste0("ALL_METRICS_together_long_NO_THINNED_", thin_level,"_",NO_THINNED_site, "_depurat.txt")), sep = "\t", row.names = FALSE)
    
    
    alldata1_list[[length(alldata1_list) + 1]] <- all_data_wide1
    alldata2_list[[length(alldata2_list) + 1]] <- all_data_wide2
    
  }
}

# Combine all results
all_differences_final <- RbindAll (all_differences_list)
all_differences_final1 <- RbindAll (all_differences_list1)
all_data_tot <- RbindAll (alldata1_list)
all_data_depurat <- RbindAll (alldata2_list)

crown_dir<- "WATER_LAI"

all_differences_final1$tree_site <-paste0(all_differences_final1$treeID1,"_", all_differences_final1$sites)

# Count the frequency of each treeID1 per site
treeID_counts <- all_differences_final1 %>%
  dplyr::group_by(tree_site) %>%
  tally(name = "frequency")

# Display treeID counts for each site
print(treeID_counts)

# Filter treeID1 replicated exactly 10 times across all sites
selected_treeIDs <- treeID_counts %>%
  dplyr::group_by(tree_site) %>%
  summarise(total_count = sum(frequency)) %>%
  dplyr::filter(total_count == 10) %>%
  pull(tree_site)

# Subset the data for the selected treeIDs
all_differences_final2 <- all_differences_final1 %>%
  dplyr::filter(tree_site %in% selected_treeIDs) %>%
  droplevels()


# Save the final result
write.table(all_differences_final, file.path(crown_dir, "ALL_DIFFERENCES_NO_THINNED_together_long_tot.txt"), sep = "\t", row.names = FALSE)
write.table(all_differences_final2, file.path(crown_dir, "ALL_DIFFERENCES_NO_THINNED_together_long_depurated.txt"), sep = "\t", row.names = FALSE)

all_data_depurat$tree_site <-paste0(all_data_depurat$treeID1,"_", all_data_depurat$sites)

# Count the frequency of each treeID1 per site
treeID_counts <- all_data_depurat %>%
  dplyr::group_by(tree_site) %>%
  tally(name = "frequency")

# Display treeID counts for each site
print(treeID_counts)

# Filter treeID1 replicated exactly 10 times across all sites
selected_treeIDs <- treeID_counts %>%
  dplyr::group_by(tree_site) %>%
  summarise(total_count = sum(frequency)) %>%
  dplyr::filter(total_count == 10) %>%
  pull(tree_site)

# Subset the data for the selected treeIDs
all_data_depurat1 <- all_data_depurat %>%
  dplyr::filter(tree_site %in% selected_treeIDs) %>%
  droplevels()

table(all_data_depurat1$sites, all_data_depurat1$thin_level)


# Write combined dataframe to a file
write.table(all_data_tot, file.path(crown_dir, "ALL_METRICS_together_long_NO_THINNED_tot.txt"), sep = "\t", row.names = FALSE)
write.table(all_data_depurat1, file.path(crown_dir, "ALL_METRICS_together_long_NO_THINNED_depurated.txt"), sep = "\t", row.names = FALSE)



## JOINING ALL DATA + DIFFERENCES BY SITES

crowns_dir <- "WATER_LAI"
metrics_file1 <- read.table(file.path(crown_dir,"ALL_METRICS_together_long_NO_THINNED_depurated.txt"), sep = "\t", header = TRUE)

metrics_file1$treeID <- as.character(metrics_file1$treeID)
metrics_file1$thin_level <- as.character(metrics_file1$thn_lvl)
metrics_file1$sites <- as.character(metrics_file1$sites)

table(metrics_file1$sites, metrics_file1$thin_level)

diff_file1 <- read.table(file.path(crown_dir,"ALL_DIFFERENCES_NO_THINNED_together_long_depurated.txt"),  sep = "\t", header = TRUE)
colnames(diff_file1) <- paste0(colnames(diff_file1), "_difNTH")

diff_file1$treeID <- as.character(diff_file1$treeID_difNTH)
diff_file1$thin_level <- as.character(diff_file1$thn_lvl_difNTH)
diff_file1$sites <- as.character(diff_file1$sites_difNTH)   

table(diff_file1$sites, diff_file1$thin_level)

# Define the function to reorder levels
reorder_levels <- function(group, thin_level_order) {
  group$thin_level <- factor(group$thin_level, levels = thin_level_order)
  return(group)
}

thin_level_order <- c("NO_THINNED",
                      "THINNED_100P",
                      "THINNED_50P",
                      "THINNED_25P",
                      "THINNED_10P",
                      "THINNED_5P",
                      "THINNED_4P",
                      "THINNED_3P",
                      "THINNED_2P",
                      "THINNED_1P")

check_thin_level_order <- function(group, thin_level_order) {
  all_levels <- levels(group$thin_level)
  all_levels_in_order <- all(all_levels == thin_level_order)
  return(all_levels_in_order)
}


metrics_file1$thin_level <- factor(metrics_file1$thin_level)
diff_file1$thin_level <- factor(diff_file1$thin_level)

# Reorder levels of 'thin_level' variable
metrics_file2 <- reorder_levels(metrics_file1, thin_level_order)
diff_file2 <- reorder_levels(diff_file1, thin_level_order)

metrics_file2$tree_site <- as.character(metrics_file2$tree_site)
diff_file2$tree_site <- as.character(diff_file2$tree_site)

metrics_file2$thin_level <- as.character(metrics_file2$thin_level)
diff_file2$thin_level <- as.character(diff_file2$thin_level)

all_data_diff <- metrics_file2 |> left_join(diff_file2, by = c("tree_site", "thin_level"))

cols_to_drop <- c("treeID1.y", "treeID_y", "sites.1", "sites.y")
all_data_diff <- all_data_diff[, !names(all_data_diff) %in% cols_to_drop]
names(all_data_diff) <- sub("\\.x$", "", names(all_data_diff))

table(all_data_diff$sites, all_data_diff$thin_level)

# Load necessary library
library(dplyr)

# Count the frequency of each treeID1 per site
treeID_counts <- all_data_diff %>%
  group_by(treeID1, sites) %>%
  tally(name = "frequency")

# Display treeID counts for each site
print(treeID_counts)

# Filter treeID1 replicated exactly 10 times across all sites
selected_treeIDs <- treeID_counts %>%
  group_by(treeID1, sites) %>%
  summarise(total_count = sum(frequency)) %>%
  filter(total_count == 10) %>%
  pull(treeID1)

# Subset the data for the selected treeIDs
metrics_file_selec <- all_data_diff %>%
  filter(treeID1 %in% selected_treeIDs) %>%
  droplevels()

# Ensure treeID is a factor and levels with 0 cases are removed
metrics_file_selec$treeID1 <- factor(metrics_file_selec$treeID1)

table(metrics_file_selec$sites, metrics_file_selec$thin_level)

write.table(metrics_file_selec, file.path(crowns_dir,paste0("ALL_METRICS_AND_DIFFERENCES_NO_THINNED_ALL_SITES_depurated.txt")), sep = "\t", row.names = FALSE) 



## ADDING DIFFERENCES INFORMATION TO LADDERFUELS SHAPEFILES

library(sf)

# Directory paths
crown_dir <- "WATER_LAI"
metrics_diff <- read.table(
  file.path(crown_dir, "ALL_DIFFERENCES_NO_THINNED_together_long_depurated.txt"),
  sep = "\t", header = TRUE
)

colnames(metrics_diff) <- paste0(colnames(metrics_diff), "_df")
table(metrics_diff$sites, metrics_diff$thn_lvl)
table(metrics_diff$tree_site)

# Get the unique site names from the metrics data
sites_list <- unique(metrics_diff$sites)

# Get all crown subfolders, excluding "DATA_NO_THINNED"
metrics_crowns_subfolders <- list.dirs(crown_dir, full.names = TRUE, recursive = FALSE)
metrics_crowns_subfolders <- metrics_crowns_subfolders[!grepl("DATA_NO_THINNED", metrics_crowns_subfolders)]
metrics_crowns_subfolders1 <- list.dirs(metrics_crowns_subfolders, full.names = TRUE, recursive = FALSE)

# Loop over each site
for (site in sites_list) {
  cat("Processing site:", site, "\n")
  
  # Convert `sites` to a factor to ensure consistency
  metrics_diff$sites <- factor(metrics_diff$sites)
  
  # Filter metrics data for the current site
  metrics_file_site <- subset(metrics_diff, sites == site)
  
  if (nrow(metrics_file_site) == 0) {
    cat("No data found in metrics_file_site for site:", site, "\n")
    next
  }
  
  # Ensure consistent column types
  metrics_file_site$treeID <- as.character(metrics_file_site$treeID_df)
  metrics_file_site$treeID1 <- as.character(metrics_file_site$treeID1_df)
  metrics_file_site$treeIDf <- as.character(metrics_file_site$treeIDf_df)
  metrics_file_site$thin_level <- as.character(metrics_file_site$thn_lvl_df)
  metrics_file_site$tree_site <- as.character(metrics_file_site$tree_site_df)
  
  # Loop over each subfolder to find shapefiles related to the current site
  for (subfolder in metrics_crowns_subfolders1) {
    # List all shapefiles in the current subfolder
    crown_files <- list.files(
      subfolder,
      pattern = glob2rx("ladderfuels_metrics_*.shp"),
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    )
    
    if (length(crown_files) == 0) next  # Skip if no shapefiles are found
    
    for (crown_file in crown_files) {
      cat("Processing file:", crown_file, "\n")
      
      # Read the crown metric shapefile
      crown_metrics <- st_read(crown_file, quiet = TRUE)
      
      if (nrow(crown_metrics) == 0) {
        cat("No data in crown_metrics for file:", crown_file, "\n")
        next
      }
      
      # Add sites field if missing or incorrect
      if (!"sites" %in% colnames(crown_metrics)) {
        crown_metrics$sites <- gsub(".*_(SITE[0-9]+).*", "\\1", crown_file)
        cat("Inferred site from filename:", crown_metrics$sites[1], "\n")
      }
      
      # Ensure sites is a factor and filter for the current site
      crown_metrics$sites <- factor(crown_metrics$sites)
      crown_metrics <- subset(crown_metrics, sites == site)
    
      
      # Validate geometry
      if (!all(st_is_valid(crown_metrics))) {
        crown_metrics <- st_make_valid(crown_metrics)
      }
      
      crown_metrics$tree_site <- paste0(crown_metrics$treeID1, "_", crown_metrics$sites)
      crown_metrics$thin_level <- as.character(crown_metrics$thn_lvl)
      
      # Perform the join
      joined_data <- merge(
        crown_metrics, metrics_file_site,
        by = c("tree_site", "thin_level"),
        all.x = TRUE
      )
      
      if (nrow(joined_data) == 0) {
        cat("No joined data for site:", site, "in file:", crown_file, "\n")
        next
      }
      
      # Drop redundant columns
      cols_to_drop <- c("treeIDf.y", "treeID1.y", "sites.y", "thin_level.y", "height_factor.y")
      joined_data <- joined_data[, !names(joined_data) %in% cols_to_drop]
      names(joined_data) <- sub("\\.x$", "", names(joined_data))
      
      # Remove rows with NA in the tree_site2 column
      if ("tree_site_df" %in% colnames(joined_data)) {
        joined_data1 <- subset(joined_data, !is.na(tree_site_df))
      } else {
        cat("Warning: 'tree_site_df' column not found in joined_data. Skipping NA removal.\n")
      }
      
      # Convert to sf object
      joined_data_sf <- st_as_sf(joined_data1, crs = st_crs(crown_metrics))
      
      # Validate geometry
      if (!all(st_is_valid(joined_data_sf))) {
        joined_data_sf <- st_make_valid(joined_data_sf)
      }
      
      # Save the output to a new shapefile
      output_file <- gsub(".shp$", paste0("_diff.shp"), crown_file)
      st_write(joined_data_sf, output_file, quiet = TRUE, append = FALSE)
      cat("Saved joined file to:", output_file, "\n")
    }
  }
}







