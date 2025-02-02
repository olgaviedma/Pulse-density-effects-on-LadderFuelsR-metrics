
#**Impact of LiDAR Pulse Density on Forest Fuels Metrics derived using LadderFuelsR**

#**Authors: Olga Viedma and JM Moreno**

#This study evaluates the impact of LiDAR pulse density thinning on forest structure characterization using the LadderFuelsR framework. High-resolution LiDAR data were collected from Mediterranean forest sites with diverse structures and fire histories and systematically thinned to simulate varying pulse densities. Key forest metrics, including leaf area density (LAD), leaf area index (LAI), canopy base height, depth, height of fuel layers, and inter-layer distances, were derived to assess changes at different thinning levels.
#This forth section deals with the calculation of the median of LadderFuelsR metrics at different levels.

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
#23.CALCULATE MEDIAN I VALUES OF EACH VARIABLE BY TREES
#####################################################

# Define the directory path
crowns_dir <- "RESULTS/MODELS"

# Create the directory (including parent directories if they do not exist)
if (!dir.exists(crowns_dir)) {
  dir.create(crowns_dir, recursive = TRUE)
  cat("Directory created successfully:", crowns_dir, "\n")
} else {
  cat("Directory already exists:", crowns_dir, "\n")
}


all_data_lidarpod <- read.table("WATER_LAI/ALL_METRICS_together_long_NO_THINNED_depurated.txt", sep = "\t", header = TRUE)


all_data_lidarpod$treeID <- factor(all_data_lidarpod$treeID,
                                   levels=unique(all_data_lidarpod$treeID))
all_data_lidarpod$thin_level <- factor(all_data_lidarpod$thin_level,
                                       levels=unique(all_data_lidarpod$thin_level))
all_data_lidarpod$sites <- factor(all_data_lidarpod$sites,
                                  levels=unique(all_data_lidarpod$sites))

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
all_data_lidarpod$thin_level <- factor(all_data_lidarpod$thin_level, levels = new_order)

# Check the levels after reordering
levels(all_data_lidarpod$thin_level)
table(all_data_lidarpod$thin_level)

levels(all_data_lidarpod$sites)
table(all_data_lidarpod$sites)


#################################################3

all_data_lidarpod1 <- all_data_lidarpod[, c("treeID","treeID1","thin_level" ,"sites", "Hcbh1","dptf1","Hcb1_H1","Hdptf1",  "mxld_Hc","bp_Hcbh","Hcbh_br","mxld_dp","bp_dptf","mxld_ff","bp_ffds", "mxld_Hdp", "bp_Hdpt","mx_hgh_", "lai_tot", "undrst_")]

names(all_data_lidarpod1) <- c("treeID","treeID1","thin_level" ,"sites","F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh", "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist", "MXLAD_Hdepth","BP_Hdepth", "MAX_height", "total_LAI", "understory_LAI")


responses <-c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
              "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
              "MXLAD_Hdepth","BP_Hdepth", "MAX_height", "total_LAI", "understory_LAI")

write.table (all_data_lidarpod1, file.path(crowns_dir, "ALL_METRICS_together_long_NO_THINNED_depurated_names.txt"), sep = "\t", row.names =  F)


reference_values <- all_data_lidarpod1 %>%
  dplyr::select(all_of(c(responses, "treeID", "treeID1","sites"))) %>%
  group_by(treeID1,sites) %>%
  dplyr::summarize(across(where(is.numeric),
                          list(median = ~median(.x, na.rm = TRUE),
                               sd = ~sd(.x, na.rm = TRUE),
                               lower_ci = ~quantile(.x, 0.5, na.rm = TRUE),
                               upper_ci = ~quantile(.x, 0.95, na.rm = TRUE)),
                          .names = "{.col}_{.fn}"),
                   .groups = 'drop') 



####  RESHAPE TO LONG FORMAT ########3

# Select only the statistic columns for pivoting
reference_values0 <- reference_values %>%
  dplyr::select(treeID1,sites, ends_with(c("_median", "_sd", "_lower_ci", "_upper_ci")))


reference_values0_long <- reference_values0 %>%
  pivot_longer(
    cols = -c(treeID1,sites),                     # All columns except 'cluster_nt'
    names_to = c("Variable", "ValueType"),   # Split the variable name and the type (median/sd)
    names_pattern = "(.*)_(median|sd|lower_ci|upper_ci)",      # Capture variable name and either median or sd
    values_to = "Value"                      # Name the new column for values
  )

write.table (reference_values0_long, file.path(crowns_dir, "MEDIAN_VARIABLES_BY_TREES_quantiles0595_LONG_all_sites.txt"), sep = "\t", row.names =  F)


reference_values0_wide <- reference_values0_long %>%
  pivot_wider(
    names_from = ValueType,  # Move 'median' and 'sd' to new columns
    values_from = Value      # Values from the 'Value' column
  )

# Check the result
reference_values0_wide

write.table (reference_values0_wide, file.path(crowns_dir, "MEDIAN_VARIABLES_BY_TREES_quantiles0595_WIDE_all_sites.txt"), sep = "\t", row.names =  F)



#####################################################
#24. CALCULATE MEDIAN II VALUES OF EACH VARIABLE DIRECTLY
#####################################################

crowns_dir <- "RESULTS/MODELS"

all_data_lidarpod1 <- read.table(file.path(crowns_dir, "ALL_METRICS_together_long_NO_THINNED_depurated_names.txt"), sep = "\t", header = TRUE)


all_data_lidarpod1$treeID <- factor(all_data_lidarpod1$treeID,
                                    levels=unique(all_data_lidarpod1$treeID))
all_data_lidarpod1$thin_level <- factor(all_data_lidarpod1$thin_level,
                                        levels=unique(all_data_lidarpod1$thin_level))
all_data_lidarpod1$sites <- factor(all_data_lidarpod1$sites,
                                   levels=unique(all_data_lidarpod1$sites))

table(all_data_lidarpod1$sites, all_data_lidarpod1$thin_level)

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
all_data_lidarpod1$thin_level <- factor(all_data_lidarpod1$thin_level, levels = new_order)

responses <-c("F1_cbh", "F1_depth", "F1_LAD","F1_Hdepth","MXLAD_cbh", "BP_cbh", "BR_cbh",
              "MXLAD_depth", "BP_depth", "MXLAD_dist", "BP_dist",
              "MXLAD_Hdepth","BP_Hdepth", "MAX_height", "total_LAI", "understory_LAI")




## CALCULATE MEDIAN VALUES OF EACH VARIABLE BY THINNED LEVELS (OVER ALL TREES )


filtered_data <- all_data_lidarpod1 %>%
  dplyr::select(all_of(c(responses, "thin_level"))) %>%
  group_by(thin_level) %>%
  dplyr::summarize(across(where(is.numeric),
                          list(median = ~median(.x, na.rm = TRUE),
                               sd = ~sd(.x, na.rm = TRUE),
                               lower_ci = ~quantile(.x, 0.05, na.rm = TRUE),
                               upper_ci = ~quantile(.x, 0.95, na.rm = TRUE)),
                          .names = "{.col}_{.fn}"),
                   .groups = 'drop') 

# Select only the median columns for pivoting
filtered_data1 <- filtered_data %>%
  dplyr::select(thin_level, ends_with("_median"))

write.table (filtered_data1, file.path(crowns_dir, "MEDIAN_VARIABLES_BY_THIN_LEVELS_quantiles0595_WIDE_all_sites.txt"), sep = "\t", row.names =  F)




## CALCULATE MEDIAN VALUES OF EACH VARIABLE  (OVER ALL TREES AND OVER ALL THINNED LEVELS)


filtered_data2 <- all_data_lidarpod1 %>%
  dplyr::select(all_of(c(responses))) %>%
  dplyr::summarize(across(where(is.numeric),
                          list(median = ~median(.x, na.rm = TRUE),
                               sd = ~sd(.x, na.rm = TRUE),
                               lower_ci = ~quantile(.x, 0.05, na.rm = TRUE),
                               upper_ci = ~quantile(.x, 0.95, na.rm = TRUE)),
                          .names = "{.col}_{.fn}"),
                   .groups = 'drop') 



# Selecting only the median and related columns
filtered_data_long1 <- filtered_data2 %>%
  dplyr::select(ends_with(c("_median", "_sd", "_lower_ci", "_upper_ci")))

# Check the column names to ensure the selection was correct
print(colnames(filtered_data_long1))

# Pivoting longer
filtered_data_long2 <- filtered_data_long1 %>%
  pivot_longer(
    cols = everything(),  # Use everything() to select all columns
    names_to = c("Variable", "ValueType"),  # Split the variable name and the type (median/sd)
    names_pattern = "(.*)_(median|sd|lower_ci|upper_ci)",  # Corrected the pattern
    values_to = "Value"  # Name the new column for values
  )

# Display the result
print(filtered_data_long2)

filtered_data_sensor_wide <- filtered_data_long2 %>%
  pivot_wider(
    names_from = ValueType,  # Move 'median' and 'sd' to new columns
    values_from = Value      # Values from the 'Value' column
  )

# Check the result
print(filtered_data_sensor_wide)

write.table (filtered_data_sensor_wide, file.path(crowns_dir, "MEDIAN_VARIABLES_quantiles0595_WIDE_all_sites.txt"), sep = "\t", row.names =  F)










